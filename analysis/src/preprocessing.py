'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import numpy as np
import yaml
from abc import ABC, abstractmethod

from ROOT import TFile, TLorentzVector

import sys
sys.path.append('..')
sys.path.append('../..')

from torchic import AxisSpec, Dataset
from torchic.physics import ITS
from torchic.physics.calibration import cluster_size_parametrisation
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

np_cluster_size_parametrisation = np.vectorize(cluster_size_parametrisation)

from utils.particles import ParticleMasses, ParticlePID, ParticleLabels
PIDlabels = {int(ParticlePID[key]): ParticleLabels[key] for key in ParticlePID.keys()}

class Preprocessor(ABC):

    def __init__(self, dataset: Dataset) -> None:
        '''
            - dataset (pd.DataFrame):  data to be preprocessed
            - recoDataset (pd.DataFrame): reconstructed tracks
            - dataset['full'] (pd.DataFrame): generated and reconstructed particles

        '''
        
        print()
        self.dataset = dataset

    # PUBLIC METHODS
    
    @abstractmethod 
    @timeit
    def define_variables(self) -> None:
        
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Defining variables')

        # cut in pseudorapidity
        self.dataset.query('-0.9 < fEtaHe3 < 0.9', inplace=True)
        self.dataset.query('-0.9 < fEtaPr < 0.9', inplace=True)

        ## definition of reconstructed variables
        self.dataset.eval('fSignedPtHe3 = fPtHe3', inplace=True)
        self.dataset.eval('fSignedPtPr = fPtPr', inplace=True)
        self.dataset.eval('fPtHe3 = abs(fPtHe3)', inplace=True)
        self.dataset.eval('fPtPr = abs(fPtPr)', inplace=True)
        self.dataset.eval('fSignHe3 = fSignedPtHe3/fPtHe3', inplace=True)
        self.dataset.eval('fSignPr = fSignedPtPr/fPtPr', inplace=True)

        self.dataset.eval('fPxHe3 = fPtHe3 * cos(fPhiHe3)', inplace=True)
        self.dataset.eval('fPxPr = fPtPr * cos(fPhiPr)', inplace=True)
        self.dataset.eval('fPyHe3 = fPtHe3 * sin(fPhiHe3)', inplace=True)
        self.dataset.eval('fPyPr = fPtPr * sin(fPhiPr)', inplace=True)
        self.dataset.eval('fPzHe3 = fPtHe3 * sinh(fEtaHe3)', inplace=True)
        self.dataset.eval('fPzPr = fPtPr * sinh(fEtaPr)', inplace=True)
        self.dataset.eval('fPHe3 = fPtHe3 * cosh(fEtaHe3)', inplace=True)
        self.dataset.eval('fPPr = fPtPr * cosh(fEtaPr)', inplace=True)
        self.dataset.eval(f'fEHe3 = sqrt(fPHe3**2 + {ParticleMasses["He"]}**2)', inplace=True)
        self.dataset.eval(f'fEPr = sqrt(fPPr**2 + {ParticleMasses["Pr"]}**2)', inplace=True)
        self.dataset.eval('fDeltaEta = fEtaHe3 - fEtaPr', inplace=True)
        self.dataset.eval('fDeltaPhi = fPhiHe3 - fPhiPr', inplace=True)

        self.dataset.eval('fInnerParamTPCHe3 = fInnerParamTPCHe3 * 2', inplace=True)
        self.dataset.eval('fSignedInnerParamTPCHe3 = fInnerParamTPCHe3 * fSignHe3', inplace=True)
        self.dataset.eval('fSignedInnerParamTPCPr = fInnerParamTPCPr * fSignPr', inplace=True)
        
        self.dataset.eval('fCosLambdaHe3 = 1/cosh(fEtaHe3)', inplace=True)
        self.dataset.eval('fCosLambdaPr = 1/cosh(fEtaPr)', inplace=True)
        
        self.dataset['fClSizeITSMeanHe3'], self.dataset['fNHitsITSHe3'] = ITS.average_cluster_size(self.dataset['fItsClusterSizeHe3'])
        self.dataset['fClSizeITSMeanPr'], self.dataset['fNHitsITSPr'] = ITS.average_cluster_size(self.dataset['fItsClusterSizePr'])
        self.dataset.eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 * fCosLambdaHe3', inplace=True)
        self.dataset.eval('fClSizeITSCosLamPr = fClSizeITSMeanPr * fCosLambdaPr', inplace=True)

        # beta*gamma (for Bethe-Bloch formula)
        self.dataset.eval(f'fBetaGammaHe3 = fInnerParamTPCHe3 * 2 / {ParticleMasses["He"]}', inplace=True)
        self.dataset.eval(f'fBetaGammaPr = fInnerParamTPCPr / {ParticleMasses["Pr"]}', inplace=True)

        # invariant mass 
        self.dataset.eval('fPtLi = sqrt(fPtHe3**2 + fPtPr**2 + 2*fPtHe3*fPtPr*cos(fDeltaPhi))', inplace=True)
        self.dataset.eval('fSignedPtLi = fPtLi * fSignHe3', inplace=True)    
        self.dataset.eval('fMassInvLi = sqrt( (fEHe3 + fEPr)**2 - ((fPxHe3 + fPxPr)**2 + (fPyHe3 + fPyPr)**2 + (fPzHe3 + fPzPr)**2) )', inplace=True)

        # TOF nsigma
        expTOFmassPr = 0.9487
        self.dataset['fExpTOFMassPr'] = expTOFmassPr
       
        # exponential
        resolution_params = {
            'kp0': 1.22204e-02,
            'kp1': 7.48467e-01,
        }
        self.dataset.eval(f'fSigmaTOFMassPr = ({resolution_params["kp0"]} * exp({resolution_params["kp1"]} * abs(fPtPr))) * {expTOFmassPr}', inplace=True)
        self.dataset.eval(f'fNSigmaTOFPr = (fMassTOFPr - fExpTOFMassPr) / fSigmaTOFMassPr', inplace=True)

        # separate matter and antimatter
        self.dataset.add_subset('antimatter', self.dataset._data['fSignHe3'] < 0)
        self.dataset.add_subset('matter', self.dataset._data['fSignHe3'] > 0)
        self.dataset.add_subset('unlike-sign', self.dataset._data['fIsBkgLS'] == True) # this is actually isBkgUS
        self.dataset.add_subset('like-sign', self.dataset._data['fIsBkgLS'] == False) # this is actually isBkgUS
    
    @staticmethod
    def compute_kstar(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
        '''
            Kstar: variable used to study the correlation between two particles. It is computed as half of the difference between the 
            momentum of the two particles in the center of mass reference frame.

            NOTE: this is not vectorized

            Parameters:
            - pt1, eta1, phi1, e1: 4-momentum of the first particle
            - pt2, eta2, phi2, e2: 4-momentum of the second particle
        '''

        p1mu = TLorentzVector(0, 0, 0, 0)
        p1mu.SetPtEtaPhiM(pt1, eta1, phi1, m1)
        p2mu = TLorentzVector(0, 0, 0, 0)
        p2mu.SetPtEtaPhiM(pt2, eta2, phi2, m2)
        Pboost = (p1mu + p2mu).BoostVector()

        p1mu_star = TLorentzVector(0, 0, 0, 0)
        p1mu_star.SetPtEtaPhiM(pt1, eta1, phi1, m1)
        p2mu_star = TLorentzVector(0, 0, 0, 0)
        p2mu_star.SetPtEtaPhiM(pt2, eta2, phi2, m2)

        p1mu_star.Boost(-Pboost)
        p2mu_star.Boost(-Pboost)

        kmu_star = p1mu_star - p2mu_star
        kstar = 0.5 * kmu_star.P()

        del p1mu, p2mu, Pboost, p1mu_star, p2mu_star, kmu_star
        return kstar
    
    np_compute_kstar = np.vectorize(compute_kstar)

    @timeit
    def define_kstar(self):
        '''
            Kstar: variale used for the study of correlation of p-He3
        '''

        self.dataset.query('fPtPr < 10 and fPtHe3 < 10', inplace=True)
        self.dataset['fKstar'] = self.np_compute_kstar(self.dataset['fPtHe3'], self.dataset['fEtaHe3'], self.dataset['fPhiHe3'], ParticleMasses['He'], 
                                                  self.dataset['fPtPr'], self.dataset['fEtaPr'], self.dataset['fPhiPr'], ParticleMasses['Pr'])
    
    @abstractmethod
    def correct_pt_H3_hp(self) -> None: 
        '''
            Corrected pT for He3 identified as H3 in tracking.
        '''

        print(tc.RED+'[ERROR]: '+tc.RESET+'Method not implemented - .loc method needs to be implemented in Dataset!')

        curve_params = {'kp0': -0.233625,
                       'kp1': 0.12484,
                       'kp2': -0.015673
                       }
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Correcting pT for He3 identified as H3 in tracking')
        print('Using pol2 correction')
        print('Parameters:', curve_params)

        # change values only to rows where fPIDtrkHe3 == 6
        # pol1 correction
        # self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] * (1 + curve_params['kp0'] + curve_params['kp1'] * self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'])
        
        # pol2 correction
        self._dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] + curve_params['kp0'] + curve_params['kp1'] * self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] + curve_params['kp2'] * self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3']**2
      
    def selections_He3(self) -> None:
        '''
            Selections for He3
        '''

        ItsClSizeParams = {
            'kp1': 2.781,
            'kp2': 1.159,
            'kp3': 5.116,
        }
        
        self.dataset['fExpClSizeITSHe3'] = np_cluster_size_parametrisation(self.dataset['fBetaGammaHe3'], *ItsClSizeParams.values())
        its_resolution = 0.11 # 11% resolution
        self.dataset.eval(f'fSigmaClSizeCosLHe3 = fExpClSizeITSHe3 * {its_resolution}', inplace=True)
        self.dataset.eval('fNSigmaITSHe3 = (fClSizeITSCosLamHe3 - fExpClSizeITSHe3) / fSigmaClSizeCosLHe3', inplace=True) 

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on He3')
        #self.dataset.query('fClSizeITSCosLamHe3 > 5.5', inplace=True)
        self.dataset.query('fNSigmaITSHe3 > -1.5', inplace=True)
        self.dataset.query('-2 < fNSigmaTPCHe3 < 2', inplace=True)
        self.dataset.query('0.5 < fChi2TPCHe3 < 4', inplace=True)
        self.dataset.query('abs(fInnerParamTPCHe3) > 1.6', inplace=True)
        self.dataset.query('abs(fDCAxyHe3) < 0.1', inplace=True)

    def selections_Pr(self) -> None:
        '''
            Selections for Pr
        '''

        ItsClSizeParams = {
            'kp1': 0.903,
            'kp2': 2.014,
            'kp3': 2.440,
        }
        
        self.dataset['fExpClSizeITSPr'] = np_cluster_size_parametrisation(self.dataset['fBetaGammaPr'], *ItsClSizeParams.values())
        its_resolution = 0.2 # 20% resolution
        self.dataset.eval(f'fSigmaClSizeCosLPr = fExpClSizeITSPr * {its_resolution}', inplace=True)
        self.dataset.eval('fNSigmaITSPr = (fClSizeITSCosLamPr - fExpClSizeITSPr) / fSigmaClSizeCosLPr', inplace=True)         

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on Pr')
        self.dataset.query('0.5 < fChi2TPCPr < 4', inplace=False)
        self.dataset.query('fSignedPtPr > -3.', inplace=False)

        self.dataset.query('(fPtPr < 0.8) or (-2 < fNSigmaTOFPr < 2)', inplace=True)
        self.dataset.query('fNSigmaITSPr > -1.5', inplace=True)
        self.dataset.query('abs(fNSigmaTPCPr) < 2', inplace=True)
   
class DataPreprocessor(Preprocessor):

    def __init__(self, dataset: Dataset) -> None:
        '''
            Method inheritated from Processor class. Data is stored in the 'full' key.
            'reco' key is not used.
        '''
        super().__init__(dataset)
        self.dataset.add_subset('full', self.dataset._data['fPtPr'] > -999.) # dummy check

    def define_variables(self) -> None:
        '''
            Definition of variables for data.
        '''
        super().define_variables()

    def correct_pt_H3_hp(self) -> None:
        '''
            Corrected pT for He3 identified as H3 in tracking using parameters from a dedicated curve
            (see studies.py H3inTrkStudy class)
        '''

        super().correct_pt_H3_hp()

    def visualize(self, outputFilePath, config) -> None:
        ''' 
            Visualization of data.
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Visualizing')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Creating output file '+tc.UNDERLINE+tc.CYAN+outputFilePath+tc.RESET)
        out_file = TFile(outputFilePath, 'recreate')
        out_dirs = {}
        for dir in config['outDirs']:   out_dirs[dir] = out_file.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outDirs':                continue
            
            for part in cfg['particle']:

                opt = cfg.get('opt', 'full')
                    
                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset.columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable']+part,'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.dataset.build_hist(cfg['xVariable']+part, axis_spec_x, subset=opt)
                    
                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset.columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable']+part,'not present in dataset!')
                        continue
                    elif cfg['yVariable']+part not in self.dataset.columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['yVariable']+part,'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axis_spec_y = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.dataset.build_hist(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y, subset=opt)
                    #if cfg['xVariable'] == 'fPIDtrk':   hist = self.hist_handler[opt].setLabels(hist, PIDlabels, 'x')
                    #if cfg['yVariable'] == 'fPIDtrk':   hist = self.hist_handler[opt].setLabels(hist, PIDlabels, 'y')

                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()
        
        out_file.Close()


class MCPreprocessor(Preprocessor):

    def __init__(self, dataset: Dataset) -> None:
        '''
            Method inheritated from Processor class. Full sample is stored in the 'full' key.
            Reconstructed tracks are stored in the 'reco' key.
        '''

        super().__init__(dataset)        
        self.dataset.add_subset('reco', self.dataset._data['fPIDtrkHe3'] != 0xFFFFF)
    
    # PUBLIC METHODS

    def define_variables(self) -> None:
        '''
            Definition of variables for MC data.
        '''
        super().define_variables()

        ## MC variables
        self.dataset.eval('fPtMCLi = fSignedPtMC', inplace=True)
        self.dataset.eval('fMassMCLi = fMassMC', inplace=True)
        self.dataset.eval('fSignLi = (-1)**(fPtMCLi < 0)', inplace=True)

        self.dataset.eval('fSignedPtMCHe3 = fPtMCHe3', inplace=True)
        self.dataset.eval('fSignedPtMCPr = fPtMCPr', inplace=True)
        self.dataset.eval('fSignedPtMCLi = fPtMCLi', inplace=True)

        self.dataset.eval('fPtMCHe3 = abs(fPtMCHe3)', inplace=True)
        self.dataset.eval('fPtMCPr = abs(fPtMCPr)', inplace=True)
        self.dataset.eval('fPtMCLi = abs(fPtMCLi)', inplace=True)

        self.dataset.eval('fPMCHe3 = fPtMCHe3 * cosh(fEtaMCHe3)', inplace=True)
        self.dataset.eval('fPMCPr = fPtMCPr * cosh(fEtaMCPr)', inplace=True)
        #self.dataset['fPMCLi'] = self.dataset['fPtMCLi'] * np.cosh(self.dataset['fEtaLi'])

        self.dataset.eval('fDeltaEtaMC = fEtaMCHe3 - fEtaMCPr', inplace=True)
        self.dataset.eval('fDeltaPhiMC = fPhiMCHe3 - fPhiMCPr', inplace=True)

        # momentum resolution
        self.define_resolution()

        ## reconstructed variables
        self.dataset.eval(f'fEMCHe3 = sqrt(fPMCHe3**2 + {ParticleMasses["He"]}**2)', inplace=True)
        self.dataset.eval(f'fEMCPr = sqrt(fPMCPr**2 + {ParticleMasses["Pr"]}**2)', inplace=True)
        self.dataset.eval('fMassInvMCLi = sqrt((fEMCHe3 + fEMCPr)**2 - (fPtMCHe3 * cos(fPhiMCHe3) + fPtMCPr * cos(fPhiMCPr))**2 - (fPtMCHe3 * sin(fPhiMCHe3) + fPtMCPr * sin(fPhiMCPr))**2 - (fPtMCHe3 * sinh(fEtaMCHe3) + fPtMCPr * sinh(fEtaMCPr))**2 )', inplace=True)

    def define_resolution(self) -> None:
        '''
            Definition of resolution in the dataframe
        '''

        ## pt resolution
        self.dataset.eval('fPtResHe3 = (fPtMCHe3 - fPtHe3) / fPtMCHe3', inplace=True)
        self.dataset.eval('fPtResPr = (fPtMCPr - fPtPr) / fPtMCPr', inplace=True)
        self.dataset.eval('fPtResLi = (fPtMCLi - fPtLi) / fPtMCLi', inplace=True)
        self.dataset.eval('fPtResNotNormHe3 = (fPtMCHe3 - fPtHe3)', inplace=True)
        self.dataset.eval('fPtResNotNormPr = (fPtMCPr - fPtPr)', inplace=True)
        self.dataset.eval('fPtResNotNormLi = (fPtMCLi - fPtLi)', inplace=True)

        ## p resolution
        self.dataset.eval('fPResHe3 = (fPMCHe3 - fInnerPTPCHe3) / fPMCHe3', inplace=True)
        self.dataset.eval('fPResPr = (fPMCPr - fInnerPTPCPr) / fPMCPr', inplace=True)

    def define_kstar(self):
        super().define_kstar()

    def correct_pt_H3_hp(self) -> None:
        '''
            Corrected pT for He3 identified as H3 in tracking using parameters from a dedicated curve
            (see studies.py H3inTrkStudy class)
        '''

        super().correct_pt_H3_hp()
        self.define_resolution()

    def visualize(self, config) -> None:
        ''' 
            Visualization of MC.
        '''

        super().visualize(config)

        print(tc.GREEN+'[INFO]:'+tc.RESET+' visualizing')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print(tc.GREEN+'[INFO]:'+tc.RESET+' creating output file '+tc.UNDERLINE+tc.CYAN+f'{config["outputFilePath"]}'+tc.RESET)
        out_file = TFile(config['outputFilePath'], 'recreate')
        out_dirs = {}
        for dir in config['outDirs']:   out_dirs[dir] = out_file.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outputFilePath':         continue
            if key == 'studiesOutputFilePath':  continue
            if key == 'out_dirs':                continue
            
            for part in cfg['particle']:

                opt = cfg.get('opt', 'full')

                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    hist = self.dataset.build_hist(cfg['xVariable']+part, axis_spec_x, subset=opt)
                    
                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:
                    
                    if cfg['xVariable']+part not in self.dataset.columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue
                    elif cfg['yVariable']+part not in self.dataset.columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['yVariable'],'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axis_spec_y = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')
                    hist = self.dataset.build_hist(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y, subset=opt)

                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()

                if 'TEfficiency' in cfg['type']:
                    print(tc.RED+'[ERROR]:'+tc.RESET+' TEfficiency not implemented yet!')

                #    axisSpecNum = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+'Reco'+part, cfg['title']+f' {part}')
                #    axisSpecDen = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                #    hNumerator = self.hist_handler['reco'].buildTH1(cfg['numVariable']+part, axisSpecNum)
                #    hDenominator = self.hist_handler['full'].buildTH1(cfg['denVariable']+part, axisSpecDen)
#
                #    eff = self.hist_handler['reco'].buildEfficiency(hNumerator, hDenominator)
                #    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                #    else:                       out_file.cd()  
                #    eff.Write()
        #
        out_file.Close()