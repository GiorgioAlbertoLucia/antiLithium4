'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import numpy as np
import yaml
from copy import deepcopy
from abc import ABC, abstractmethod

from ROOT import TFile, TLorentzVector
#from ROOT.Math import Vector4D
from ROOT.Math import LorentzVector, PtEtaPhiMVector, Boost

import sys
sys.path.append('..')
sys.path.append('../..')

from framework.src.axis_spec import AxisSpec
from framework.src.hist_handler import *
from framework.src.data_handler import *
from framework.utils.timeit import timeit

from utils.particles import ParticleMasses, ParticlePID, ParticleLabels
PIDlabels = {int(ParticlePID[key]): ParticleLabels[key] for key in ParticlePID.keys()}

class Preprocessor(ABC):

    def __init__(self, data_handler: DataHandler) -> None:
        '''
            - dataset (pd.DataFrame):  data to be preprocessed
            - recoDataset (pd.DataFrame): reconstructed tracks
            - dataset['full'] (pd.DataFrame): generated and reconstructed particles

        '''
        
        print()
        self.data_handler = data_handler
        tmp_dataset = data_handler.inData.copy()

        self._dataset = {'full': tmp_dataset,
                        'reco': tmp_dataset,
                        'matter': tmp_dataset,
                        'antimatter': tmp_dataset}

    def _update_antimatter(self) -> None:
        '''
            Filter out particles with negative charge.
        '''
        if 'fSignHe3' not in self._dataset['full'].columns or 'fSignPr' not in self._dataset['full'].columns:
            print(tc.RED+'[ERROR]: '+tc.RESET+'Columns fSignHe3 or fSignPr not present in dataset')
        self._dataset['antimatter'] = self._dataset['full'].query('fSignHe3 < 0 and fSignPr < 0', inplace=False)
        self._dataset['matter'] = self._dataset['full'].query('fSignHe3 > 0 and fSignPr > 0', inplace=False)
    
    def update_antimatter(func):
        def wrapper(self, *args, **kwargs):
            result = func(self, *args, **kwargs)
            self._update_antimatter()
            return result
        return wrapper

    @property
    def dataset(self):
        return self._dataset

    # PUBLIC METHODS
    
    @abstractmethod 
    @update_antimatter
    @timeit
    def define_variables(self) -> None:
        
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Defining variables')

        # cut in pseudorapidity
        self._dataset['full'] = self._dataset['full'].query('-0.9 < fEtaHe3 < 0.9', inplace=False)
        self._dataset['full'] = self._dataset['full'].query('-0.9 < fEtaPr < 0.9', inplace=False)

        # convert rigidity to momentum (for He3, if not done in the task)
        # self._dataset['full']['fPtHe3'] = self._dataset['full']['fPtHe3'] * 2.

        ## definition of reconstructed variables
        self._dataset['full']['fSignedPtHe3'] = self._dataset['full']['fPtHe3']
        self._dataset['full']['fSignedPtPr'] = self._dataset['full']['fPtPr']
        self._dataset['full']['fPtHe3'] = abs(self._dataset['full']['fPtHe3'])
        self._dataset['full']['fPtPr'] = abs(self._dataset['full']['fPtPr'])
        self._dataset['full'].eval('fSignHe3 = fSignedPtHe3/fPtHe3', inplace=True)
        self._dataset['full'].eval('fSignPr = fSignedPtPr/fPtPr', inplace=True)

        self._dataset['full'].eval('fPxHe3 = fPtHe3 * cos(fPhiHe3)', inplace=True)
        self._dataset['full'].eval('fPxPr = fPtPr * cos(fPhiPr)', inplace=True)
        self._dataset['full'].eval('fPyHe3 = fPtHe3 * sin(fPhiHe3)', inplace=True)
        self._dataset['full'].eval('fPyPr = fPtPr * sin(fPhiPr)', inplace=True)
        self._dataset['full'].eval('fPzHe3 = fPtHe3 * sinh(fEtaHe3)', inplace=True)
        self._dataset['full'].eval('fPzPr = fPtPr * sinh(fEtaPr)', inplace=True)
        self._dataset['full'].eval('fPHe3 = fPtHe3 * cosh(fEtaHe3)', inplace=True)
        self._dataset['full'].eval('fPPr = fPtPr * cosh(fEtaPr)', inplace=True)
        self._dataset['full'].eval(f'fEHe3 = sqrt(fPHe3**2 + {ParticleMasses["He"]}**2)', inplace=True)
        self._dataset['full'].eval(f'fEPr = sqrt(fPPr**2 + {ParticleMasses["Pr"]}**2)', inplace=True)
        self._dataset['full'].eval('fDeltaEta = fEtaHe3 - fEtaPr', inplace=True)
        self._dataset['full'].eval('fDeltaPhi = fPhiHe3 - fPhiPr', inplace=True)

        self._dataset['full']['fInnerPTPCHe3'] = self._dataset['full']['fInnerParamTPCHe3'] * 2
        self._dataset['full']['fInnerPTPCPr'] = self._dataset['full']['fInnerParamTPCPr']

        self._dataset['full'].eval('fSignedPTPCHe3 = fInnerPTPCHe3 * fSignHe3', inplace=True)
        self._dataset['full'].eval('fSignedPTPCPr = fInnerPTPCPr * fSignPr', inplace=True)
        
        # ITS cluster size 
        for layer in range(7):
            read4Bits = lambda x, layer: (x >> layer*4) & 0b1111
            self._dataset['full'][f'fClSizeITS{layer}He3'] = self._dataset['full']['fItsClusterSizeHe3'].apply(read4Bits, args=(layer,))
            self._dataset['full'][f'fClSizeITS{layer}Pr'] = self._dataset['full']['fItsClusterSizePr'].apply(read4Bits, args=(layer,))
        
        self._dataset['full'].eval('fCosLambdaHe3 = 1/cosh(fEtaHe3)', inplace=True)
        self._dataset['full'].eval('fCosLambdaPr = 1/cosh(fEtaPr)', inplace=True)
        
        self._dataset['full'].eval('fNHitsITSHe3 = (fClSizeITS0He3 > 0) * 1. + (fClSizeITS1He3 > 0) * 1. + (fClSizeITS2He3 > 0) * 1. + (fClSizeITS3He3 > 0) * 1. + (fClSizeITS4He3 > 0) * 1. + (fClSizeITS5He3 > 0) * 1. + (fClSizeITS6He3 > 0) * 1.', inplace=True)
        self._dataset['full'].eval('fNHitsITSPr = (fClSizeITS0Pr > 0) * 1. + (fClSizeITS1Pr > 0) * 1. + (fClSizeITS2Pr > 0) * 1. + (fClSizeITS3Pr > 0) * 1. + (fClSizeITS4Pr > 0) * 1. + (fClSizeITS5Pr > 0) * 1. + (fClSizeITS6Pr > 0) * 1.', inplace=True)
        self._dataset['full'].eval('fClSizeITSMeanHe3 = (fClSizeITS0He3 + fClSizeITS1He3 + fClSizeITS2He3 + fClSizeITS3He3 + fClSizeITS4He3 + fClSizeITS5He3 + fClSizeITS6He3) / fNHitsITSHe3', inplace=True)
        self._dataset['full'].eval('fClSizeITSMeanPr = (fClSizeITS0Pr + fClSizeITS1Pr + fClSizeITS2Pr + fClSizeITS3Pr + fClSizeITS4Pr + fClSizeITS5Pr + fClSizeITS6Pr) / fNHitsITSPr', inplace=True)
        self._dataset['full'].eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 * fCosLambdaHe3', inplace=True)
        self._dataset['full'].eval('fClSizeITSCosLamPr = fClSizeITSMeanPr * fCosLambdaPr', inplace=True)

        # beta*gamma (for Bethe-Bloch formula)
        self._dataset['full'].eval(f'fBetaGammaHe3 = fInnerParamTPCHe3 * 2 / {ParticleMasses["He"]}', inplace=True)
        self._dataset['full'].eval(f'fBetaGammaPr = fInnerParamTPCPr / {ParticleMasses["Pr"]}', inplace=True)

        # invariant mass 
        self._dataset['full'].eval('fPtLi = sqrt(fPtHe3**2 + fPtPr**2 + 2*fPtHe3*fPtPr*cos(fDeltaPhi))', inplace=True)
        self._dataset['full'].eval('fSignedPtLi = fPtLi * fSignHe3', inplace=True)    
        self._dataset['full'].eval('fMassInvLi = sqrt( (fEHe3 + fEPr)**2 - ((fPxHe3 + fPxPr)**2 + (fPyHe3 + fPyPr)**2 + (fPzHe3 + fPzPr)**2) )', inplace=True)
    
    def compute_kstar(self, p1x, p1y, p1z, e1, p2x, p2y, p2z, e2):
        '''
            Kstar: variable used to study the correlation between two particles. It is computed as half of the difference between the 
            momentum of the two particles in the center of mass reference frame.

            NOTE: this is not vectorized

            Parameters:
            - p1x, p1y, p1z, e1: 4-momentum of the first particle
            - p2x, p2y, p2z, e2: 4-momentum of the second particle
        '''

        p1mu = TLorentzVector(p1x, p1y, p1z, e1)
        p2mu = TLorentzVector(p2x, p2y, p2z, e2)
        Pmu = p1mu + p2mu

        beta = Pmu.Beta()
        betax = beta * np.cos(Pmu.Phi()) * np.sin(Pmu.Theta())
        betay = beta * np.sin(Pmu.Phi()) * np.sin(Pmu.Theta())
        betaz = beta * np.cos(Pmu.Theta())

        p1mu_star, p2muStar = TLorentzVector(p1mu), TLorentzVector(p2mu)
        p1mu_star.SetXYZM(p1x, p1y, p1z, p1mu.M())
        p2mu_star.SetXYZM(p2x, p2y, p2z, p2mu.M())
        
        boost = Boost(-betax, -betay, -betaz)
        p1mu_star = boost(p1mu_star)
        p2mu_star = boost(p2mu_star)

        kmu_star = 0.5 * (p1mu_star - p2mu_star)
        kstar = kmu_star.P()
        del p1mu, p2mu, Pmu, beta, betax, betay, betaz, p1mu_star, p2mu_star, boost, kmu_star
        return kstar
    
    def compute_kstar2(self, pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
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

    @update_antimatter
    @timeit
    def define_kstar(self):
        '''
            Kstar: variale used for the study of correlation of p-He3
        '''

        self._dataset['full'].query('fPtPr < 10 and fPtHe3 < 10', inplace=True)
        self._dataset['full']['fKstar'] = self._dataset['full'].apply(lambda x: self.compute_kstar2(x['fPtHe3'], x['fEtaHe3'], x['fPhiHe3'], ParticleMasses['He'], 
                                                                                                 x['fPtPr'], x['fEtaPr'], x['fPhiPr'], ParticleMasses['Pr']), axis=1)
    
    @abstractmethod
    @update_antimatter
    def correct_pt_H3_hp(self) -> None: 
        '''
            Corrected pT for He3 identified as H3 in tracking.
        '''

        curve_params = {'kp0': -0.233625,
                       'kp1': 0.12484,
                       'kp2': -0.015673
                       }
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Correcting pT for He3 identified as H3 in tracking')
        print('Using pol2 correction')
        print('Parameters:', curve_params)

        # change values only to rows where fPIDtrkHe3 == 6
        # pol1 correction
        # self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] = self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] * (1 + curve_params['kp0'] + curve_params['kp1'] * self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'])
        # pol2 correction
        self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] = self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] + curve_params['kp0'] + curve_params['kp1'] * self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] + curve_params['kp2'] * self._dataset['full'].loc[self._dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3']**2
      
    def visualize(self, config) -> None:
        
        self.hist_handler = None # reset histHandler
        self.hist_handler = {key: HistHandler.createInstance(self._dataset[key]) for key in self._dataset.keys()}

    def remove_non_reco(self) -> None:
        '''
            Remove particles that were not reconstructed.
        '''
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Removing particles that were not reconstructed')
        self._dataset['full'].query('fPIDtrkHe3 != 0xFFFFF and fPIDtrkPr != 0xFFFFF', inplace=True)


    def ItsClSizeFunc(betagamma, ItsClSizeParams):
        return ItsClSizeParams['kp1'] / betagamma ** ItsClSizeParams['kp2'] + ItsClSizeParams['kp3']
    ItsClSizeFunc = np.vectorize(ItsClSizeFunc)

    @update_antimatter
    def selections_He3(self) -> None:
        '''
            Selections for He3
        '''

        ItsClSizeParams = {
            'kp1': 2.781,
            'kp2': 1.159,
            'kp3': 5.116,
        }
        
        self._dataset['full'].eval('fExpClSizeITSHe3 = @self.ItsClSizeFunc(fBetaGammaHe3, @ItsClSizeParams)', inplace=True)
        its_resolution = 0.11 # 11% resolution
        self._dataset['full'].eval('fSigmaClSizeCosLHe3 = fExpClSizeITSHe3 * @its_resolution', inplace=True)
        self._dataset['full'].eval('fNSigmaITSHe3 = (fClSizeITSCosLamHe3 - fExpClSizeITSHe3) / fSigmaClSizeCosLHe3', inplace=True) 

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on He3')
        #self._dataset['full'].query('fClSizeITSCosLamHe3 > 5.5', inplace=True)
        self._dataset['full'].query('fNSigmaITSHe3 > 0', inplace=True)
        self._dataset['full'].query('-2 < fNSigmaTPCHe3 < 2', inplace=True)
        self._dataset['full'].query('0.5 < fChi2TPCHe3 < 4', inplace=True)
        self._dataset['full'].query('fInnerPTPCHe3 > 1.6', inplace=True)

    @update_antimatter
    def selections_Pr(self) -> None:
        '''
            Selections for Pr
        '''

        ItsClSizeParams = {
            'kp1': 0.903,
            'kp2': 2.014,
            'kp3': 2.440,
        }
        
        self._dataset['full'].eval('fExpClSizeITSPr = @self.ItsClSizeFunc(fBetaGammaPr, @ItsClSizeParams)', inplace=True)
        its_resolution = 0.2 # 20% resolution
        self._dataset['full'].eval('fSigmaClSizeCosLPr = fExpClSizeITSPr * @its_resolution', inplace=True)
        self._dataset['full'].eval('fNSigmaITSPr = (fClSizeITSCosLamPr - fExpClSizeITSPr) / fSigmaClSizeCosLPr', inplace=True)         

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on Pr')
        self._dataset['full'].query('fNSigmaITSPr > -1', inplace=True)
        self._dataset['full'].query('0.5 < fChi2TPCPr < 4', inplace=True)
        self._dataset['full'].query('fSignedPtPr > -3.', inplace=True)

        expTOFmassPr = 0.9487
        self._dataset['full']['fExpTOFMassPr'] = expTOFmassPr
        # pol2
        #resolution_params = {
        #    'kp0': 3.80467e-02,
        #    'kp1': 3.54744e-02,
        #    'kp2': 2.21394e-02
        #}
        #self._dataset['full']['fSigmaTOFMassPr'] = (resolution_params['kp0'] + resolution_params['kp1'] * self._dataset['full']['fSignedPtPr'] + resolution_params['kp2'] * self._dataset['full']['fSignedPtPr']**2) * expTOFmassPr

        # exponential
        resolution_params = {
            'kp0': 1.22204e-02,
            'kp1': 7.48467e-01,
        }
        self._dataset['full']['fSigmaTOFMassPr'] = (resolution_params['kp0'] * np.exp(resolution_params['kp1'] * np.abs(self._dataset['full']['fPtPr']))) * expTOFmassPr
        self._dataset['full']['fNSigmaTOFPr'] = (self._dataset['full']['fMassTOFPr'] - self._dataset['full']['fExpTOFMassPr']) / self._dataset['full']['fSigmaTOFMassPr']

        self._dataset['full'].query('(fPtPr < 0.8) or (-1 < fNSigmaTOFPr < 1)', inplace=True)
        self._dataset['full'].query('(fPtPr > 0.8) or (-5 < fNSigmaTOFPr < 5)', inplace=True)

    @update_antimatter
    def single_track_id(self) -> None:
        '''
            Select every track only once.
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Selecting single track events')
        self._dataset['full'] = self._dataset['full'].drop_duplicates(subset=['fTrackIDPr'], keep='first').reset_index(drop=True)
        #self._dataset['full'] = self._dataset['full'].drop_duplicates(subset=['fTrackIDHe3'], keep='first').reset_index(drop=True)
   
class DataPreprocessor(Preprocessor):

    def __init__(self, data_handler: DataHandler) -> None:
        '''
            Method inheritated from Processor class. Data is stored in the 'full' key.
            'reco' key is not used.
        '''
        super().__init__(data_handler)

    @property
    def dataset(self):
        return super().dataset

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

        super().visualize(config)
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

                    if cfg['xVariable']+part not in self.dataset[opt].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable']+part,'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.hist_handler[opt].buildTH1(cfg['xVariable']+part, axis_spec_x)
                    if cfg['xVariable'] == 'fPIDtrk':   self.hist_handler[opt].setLabels(hist, PIDlabels, 'x')
                    
                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset[opt].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable']+part,'not present in dataset!')
                        continue
                    elif cfg['yVariable']+part not in self.dataset[opt].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['yVariable']+part,'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axis_spec_y = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.hist_handler[opt].buildTH2(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y)
                    #if cfg['xVariable'] == 'fPIDtrk':   hist = self.hist_handler[opt].setLabels(hist, PIDlabels, 'x')
                    #if cfg['yVariable'] == 'fPIDtrk':   hist = self.hist_handler[opt].setLabels(hist, PIDlabels, 'y')

                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()
        
        out_file.Close()


class MCPreprocessor(Preprocessor):

    def __init__(self, data_handler: DataHandler) -> None:
        '''
            Method inheritated from Processor class. Full sample is stored in the 'full' key.
            Reconstructed tracks are stored in the 'reco' key.
        '''

        super().__init__(data_handler)

    def __update_reco_tracks(self) -> None:
            
        self._dataset['reco'] = self._dataset['full'].copy().query('fPIDtrkHe3 != 0xFFFFF', inplace=False)
    
    # PUBLIC METHODS

    def define_variables(self) -> None:
        '''
            Definition of variables for MC data.
        '''
        super().define_variables()

        ## MC variables
        self._dataset['full']['fPtMCLi'] = self._dataset['full']['fSignedPtMC'] 
        self._dataset['full']['fMassMCLi'] = self._dataset['full']['fMassMC']

        # pt MC stores sign
        self._dataset['full'].eval('fSignLi = (-1)**(fPtMCLi < 0)', inplace=True)

        self._dataset['full']['fSignedPtMCHe3'] = self._dataset['full']['fPtMCHe3']
        self._dataset['full']['fSignedPtMCPr'] = self._dataset['full']['fPtMCPr']
        self._dataset['full']['fSignedPtMCLi'] = self._dataset['full']['fPtMCLi']

        self._dataset['full']['fPtMCHe3'] = abs(self._dataset['full']['fPtMCHe3'])
        self._dataset['full']['fPtMCPr'] = abs(self._dataset['full']['fPtMCPr'])
        self._dataset['full']['fPtMCLi'] = abs(self._dataset['full']['fPtMCLi'])

        self._dataset['full'].eval('fPMCHe3 = fPtMCHe3 * cosh(fEtaMCHe3)', inplace=True)
        self._dataset['full'].eval('fPMCPr = fPtMCPr * cosh(fEtaMCPr)', inplace=True)
        #self._dataset['full']['fPMCLi'] = self._dataset['full']['fPtMCLi'] * np.cosh(self._dataset['full']['fEtaLi'])

        self._dataset['full'].eval('fDeltaEtaMC = fEtaMCHe3 - fEtaMCPr', inplace=True)
        self._dataset['full'].eval('fDeltaPhiMC = fPhiMCHe3 - fPhiMCPr', inplace=True)

        # momentum resolution
        self.define_resolution()

        ## reconstructed variables
        self._dataset['full'].eval(f'fEMCHe3 = sqrt(fPMCHe3**2 + {ParticleMasses["He"]}**2)', inplace=True)
        self._dataset['full'].eval(f'fEMCPr = sqrt(fPMCPr**2 + {ParticleMasses["Pr"]}**2)', inplace=True)
        self._dataset['full'].eval('fMassInvMCLi = sqrt((fEMCHe3 + fEMCPr)**2 - (fPtMCHe3 * cos(fPhiMCHe3) + fPtMCPr * cos(fPhiMCPr))**2 - (fPtMCHe3 * sin(fPhiMCHe3) + fPtMCPr * sin(fPhiMCPr))**2 - (fPtMCHe3 * sinh(fEtaMCHe3) + fPtMCPr * sinh(fEtaMCPr))**2 )', inplace=True)

        self.__update_reco_tracks()

    def define_resolution(self) -> None:
        '''
            Definition of resolution in the dataframe
        '''

        ## pt resolution
        self._dataset['full']['fPtResHe3'] = (self._dataset['full']['fPtMCHe3'] - self._dataset['full']['fPtHe3']) / self._dataset['full']['fPtMCHe3']
        self._dataset['full']['fPtResPr'] = (self._dataset['full']['fPtMCPr'] - self._dataset['full']['fPtPr']) / self._dataset['full']['fPtMCPr']
        self._dataset['full']['fPtResLi'] = (self._dataset['full']['fPtMCLi'] - self._dataset['full']['fPtLi']) / self._dataset['full']['fPtMCLi']
        self._dataset['full']['fPtResNotNormHe3'] = (self._dataset['full']['fPtMCHe3'] - self._dataset['full']['fPtHe3'])
        self._dataset['full']['fPtResNotNormPr'] = (self._dataset['full']['fPtMCPr'] - self._dataset['full']['fPtPr'])
        self._dataset['full']['fPtResNotNormLi'] = (self._dataset['full']['fPtMCLi'] - self._dataset['full']['fPtLi'])

        ## p resolution
        self._dataset['full']['fPResHe3'] = (self._dataset['full']['fPMCHe3'] - self._dataset['full']['fInnerPTPCHe3']) / self._dataset['full']['fPMCHe3']
        self._dataset['full']['fPResPr'] = (self._dataset['full']['fPMCPr'] - self._dataset['full']['fInnerPTPCPr']) / self._dataset['full']['fPMCPr']
        #self._dataset['full']['fPResLi'] = (self._dataset['full']['fPMCLi'] - self._dataset['full']['fInnterPTPCLi']) / self._dataset['full']['fPMCLi']

    def define_kstar(self):
        super().define_kstar()
        self.__update_reco_tracks()

    def correct_pt_H3_hp(self) -> None:
        '''
            Corrected pT for He3 identified as H3 in tracking using parameters from a dedicated curve
            (see studies.py H3inTrkStudy class)
        '''

        super().correct_pt_H3_hp()
        self.define_resolution()
        self.__update_reco_tracks()

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

                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    
                    hist = self.hist_handler[cfg['opt']].buildTH1(cfg['xVariable']+part, axis_spec_x)
                    if cfg['xVariable'] == 'fPIDtrk':   self.hist_handler[cfg['opt']].setLabels(hist, PIDlabels, 'x')
                    
                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:
                    
                    if cfg['xVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue
                    elif cfg['yVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['yVariable'],'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axis_spec_y = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.hist_handler[cfg['opt']].buildTH2(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y)
                    if cfg['xVariable'] == 'fPIDtrk':   self.hist_handler[cfg['opt']].setLabels(hist, PIDlabels, 'x')
                    if cfg['yVariable'] == 'fPIDtrk':   self.hist_handler[cfg['opt']].setLabels(hist, PIDlabels, 'y')

                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()
                    hist.Write()

                if 'TEfficiency' in cfg['type']:

                    axisSpecNum = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+'Reco'+part, cfg['title']+f' {part}')
                    axisSpecDen = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    hNumerator = self.hist_handler['reco'].buildTH1(cfg['numVariable']+part, axisSpecNum)
                    hDenominator = self.hist_handler['full'].buildTH1(cfg['denVariable']+part, axisSpecDen)

                    eff = self.hist_handler['reco'].buildEfficiency(hNumerator, hDenominator)
                    if cfg['dir'] != 'None':    out_dirs[cfg['dir']].cd()
                    else:                       out_file.cd()  
                    eff.Write()
        
        out_file.Close()