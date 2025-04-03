'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import ctypes
import numpy as np
import yaml
from abc import ABC, abstractmethod

from ROOT import TFile
from ROOT.Math import PtEtaPhiMVector, Boost

import sys
sys.path.append('..')
sys.path.append('../..')

from torchic import AxisSpec, Dataset
from torchic.physics.ITS import average_cluster_size, sigma_its, expected_cluster_size
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

from utils.particles import ParticleMasses, ParticlePID, ParticleLabels
from utils.parametrisations import ITS as ITS_params
from utils.parametrisations import TOF as TOF_params
from utils.parametrisations import PIDforTracking as PID_params
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
        self.available_subsets = ['full']

    def add_subset(self, name: str, condition: np.ndarray) -> None:
        '''
            Add a new subset to the dataset.
        '''
        self.dataset.add_subset(name, condition)
        self.available_subsets.append(name)

    # DEFINE VARIABLES      
    
    @abstractmethod 
    @timeit
    def define_variables(self) -> None:
        '''
        # conversion from the older version of the task
        self.dataset.eval('fPtHad = fPtPr', inplace=True)
        self.dataset.eval('fEtaHad = fEtaPr', inplace=True)
        self.dataset.eval('fPhiHad = fPhiPr', inplace=True)
        self.dataset.eval('fInnerParamTPCHad = fInnerParamTPCPr', inplace=True)
        self.dataset.eval('fDCAxyHad = fDCAxyPr', inplace=True)
        self.dataset.eval('fDCAzHad = fDCAzPr', inplace=True)
        self.dataset.eval('fChi2TPCHad = fChi2TPCPr', inplace=True)
        self.dataset.eval('fItsClusterSizeHad = fItsClusterSizePr', inplace=True)
        self.dataset.eval('fPIDtrkHad = fPIDtrkPr', inplace=True)
        self.dataset.eval('fNSigmaTPCHad = fNSigmaTPCPr', inplace=True)
        self.dataset.eval('fMassTOFHad = fMassTOFPr', inplace=True)
        self.dataset.eval('fIsBkgUS = fIsBkgLS', inplace=True)
        '''
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Defining variables')

        ## definition of reconstructed variables
        self.dataset.eval('fSignedPtHe3 = fPtHe3', inplace=True)
        self.dataset.eval('fSignedPtHad = fPtHad', inplace=True)
        self.dataset.eval('fPtHe3 = abs(fPtHe3)', inplace=True)
        self.dataset.eval('fPtHad = abs(fPtHad)', inplace=True)
        self.dataset.eval('fSignHe3 = fSignedPtHe3/fPtHe3', inplace=True)
        self.dataset.eval('fSignHad = fSignedPtHad/fPtHad', inplace=True)

        if 'fIsBkgUS' not in self.dataset.columns:
            self.dataset['fIsBkgUS'] = True
            self.dataset.loc[self.dataset['fSignHe3'] == self.dataset['fSignHad'], 'fIsBkgUS'] = False

        self.correct_pt_H3_hp()
        self.dataset.eval('fSignedPtHe3 = fPtHe3 * fSignHe3', inplace=True)

        self.dataset.eval(f'fEHe3 = sqrt((fPtHe3 * cosh(fEtaHe3))**2 + {ParticleMasses["He"]}**2)', inplace=True)
        self.dataset.eval(f'fEHad = sqrt((fPtHad * cosh(fEtaHad))**2 + {ParticleMasses["Pr"]}**2)', inplace=True)
        self.dataset.eval('fDeltaEta = fEtaHe3 - fEtaHad', inplace=True)
        self.dataset.eval('fDeltaPhi = fPhiHe3 - fPhiHad', inplace=True)

        self.dataset.eval('fInnerParamTPCHe3 = fInnerParamTPCHe3 * 2', inplace=True)
        self.dataset.eval('InnerParamTPCHe3 = fInnerParamTPCHe3 * fSignHe3', inplace=True)
        self.dataset.eval('InnerParamTPCHad = fInnerParamTPCHad * fSignHad', inplace=True)
        
        self.dataset['fClSizeITSMeanHe3'], self.dataset['fNHitsITSHe3'] = average_cluster_size(self.dataset['fItsClusterSizeHe3'])
        self.dataset['fClSizeITSMeanHad'], self.dataset['fNHitsITSHad'] = average_cluster_size(self.dataset['fItsClusterSizeHad'])
        self.dataset.eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 / cosh(fEtaHe3)', inplace=True)
        self.dataset.eval('fClSizeITSCosLamHad = fClSizeITSMeanHad / cosh(fEtaHad)', inplace=True)

        # beta*gamma (for Bethe-Bloch formula)
        self.dataset.eval(f'fBetaGammaHe3 = abs(fInnerParamTPCHe3) * 2 / {ParticleMasses["He"]}', inplace=True)
        self.dataset.eval(f'fBetaGammaHad = abs(fInnerParamTPCHad) / {ParticleMasses["Pr"]}', inplace=True)

        # invariant mass 
        self.dataset.eval('fPxLi = fPtHe3 * cos(fPhiHe3) + fPtHad * cos(fPhiHad)', inplace=True)
        self.dataset.eval('fPyLi = fPtHe3 * sin(fPhiHe3) + fPtHad * sin(fPhiHad)', inplace=True)
        self.dataset.eval('fPzLi = fPtHe3 * sinh(fEtaHe3) + fPtHad * sinh(fEtaHad)', inplace=True)
        self.dataset.eval('fELi = fEHe3 + fEHad', inplace=True)
        self.dataset.eval('fPLi = sqrt(fPxLi**2 + fPyLi**2 + fPzLi**2)', inplace=True)

        self.dataset.eval('fPtLi = sqrt(fPxLi**2 + fPyLi**2)', inplace=True)
        self.dataset.eval('fKt = fPtLi / 2', inplace=True)
        self.dataset.eval('fEtaLi = arccosh(fPLi / fELi)', inplace=True)
        self.dataset.eval('fPhiLi = arctan2(fPyLi, fPxLi)', inplace=True)
        self.dataset.eval('fSignedPtLi = fPtLi * fSignHe3', inplace=True)
        self.dataset.eval('fMassInvLi = sqrt(fELi**2 - fPLi**2)', inplace=True)
        self.dataset.eval('fMassTLi = sqrt(fELi**2 - fPtLi**2)', inplace=True)

        #self.dataset.query('fMassInvLi > 3.755', inplace=True)

        # remove unnecessary columns
        self.dataset._data.drop(columns=[
                                   'fEHe3', 'fInnerParamTPCHe3', 'fClSizeITSMeanHe3', 'fItsClusterSizeHe3',
                                   'fEHad', 'fInnerParamTPCHad', 'fClSizeITSMeanHad', 'fItsClusterSizeHad',
                                   'fPxLi', 'fPyLi', 'fPzLi', 'fELi', 'fPLi', 
                                   #'fPtLi', 
                                   'fEtaLi', 'fPhiLi',
                                   ], inplace=True)

        # separate matter and antimatter
        self.add_subset('antimatter', self.dataset['fSignHe3'] < 0)
        self.add_subset('matter', self.dataset['fSignHe3'] > 0)
        
        #if 'fCentralityFT0C' in self.dataset.columns:
        #    self.add_subset('antimatter-cent0-10', (self.dataset['fSignHe3'] < 0) & (0 < self.dataset['fCentralityFT0C']) & (self.dataset['fCentralityFT0C'] < 10))
        #    self.add_subset('matter-cent0-10', (self.dataset['fSignHe3'] > 0) & (0 < self.dataset['fCentralityFT0C']) & (self.dataset['fCentralityFT0C'] < 10))
        
        if 'fIsBkgUS' in self.dataset.columns:
            self.dataset['fIsBkgUSfloat'] = self.dataset['fIsBkgUS'].astype(float)
            self.add_subset('unlike-sign', self.dataset['fIsBkgUSfloat'] == 1)
            self.add_subset('like-sign', self.dataset['fIsBkgUSfloat'] == 0)
            if 'fCentralityFT0C' in self.dataset.columns:
                self.add_subset('unlike-sign-cent0-10', (self.dataset['fIsBkgUSfloat'] == 1) & (0 < self.dataset['fCentralityFT0C']) & (self.dataset['fCentralityFT0C'] < 10))
        
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

        p1mu = PtEtaPhiMVector(pt1, eta1, phi1, m1)
        p2mu = PtEtaPhiMVector(pt2, eta2, phi2, m2)
        P_beta_vector = (p1mu + p2mu).BoostToCM()
        P_bx, P_by, P_bz = ctypes.c_double(0), ctypes.c_double(0), ctypes.c_double(0)
        P_beta_vector.GetCoordinates(P_bx, P_by, P_bz)
        P_boost = Boost(P_bx, P_by, P_bz)

        p1mu_star = PtEtaPhiMVector(pt1, eta1, phi1, m1)
        p2mu_star = PtEtaPhiMVector(pt2, eta2, phi2, m2)

        p1mu_star = P_boost(p1mu_star)
        p2mu_star = P_boost(p2mu_star)

        kmu_star = p1mu_star - p2mu_star
        kstar = 0.5 * kmu_star.P()

        del p1mu, p2mu, P_boost, p1mu_star, p2mu_star, kmu_star
        return kstar
    
    np_compute_kstar = np.vectorize(compute_kstar)

    @timeit
    def define_kstar(self):
        '''
            Kstar: variale used for the study of correlation of p-He3
        '''
        
        self.dataset.query('fPtHad < 10 and fPtHe3 < 10', inplace=True)
        self.dataset['fKstar'] = self.np_compute_kstar(self.dataset['fPtHe3'], self.dataset['fEtaHe3'], self.dataset['fPhiHe3'], ParticleMasses['He'], 
                                                  self.dataset['fPtHad'], self.dataset['fEtaHad'], self.dataset['fPhiHad'], ParticleMasses['Pr'])
        self.add_subset('anti-kstar0-200', (self.dataset['fKstar'] < 0.2) & (self.dataset['fSignHe3'] < 0))
        self.add_subset('matter-kstar0-200', (self.dataset['fKstar'] < 0.2) & (self.dataset['fSignHe3'] > 0))
    
    def define_nsigmaTOF_Pr(self) -> None:
        self.dataset['fExpTOFMassHad'] = ParticleMasses['Pr']
        self.dataset.eval(f'fSigmaTOFMassHad = ({TOF_params.pr_res_params["res1"]} * exp({TOF_params.pr_res_params["res2"]} * abs(fPtHad))) * fExpTOFMassHad', inplace=True)
        self.dataset.eval(f'fNSigmaTOFHad = (fMassTOFHad - fExpTOFMassHad) / fSigmaTOFMassHad', inplace=True)
        self.dataset.drop(columns=['fExpTOFMassHad', 'fSigmaTOFMassHad'], inplace=True)
    
    def define_nsigmaITS_Pr(self) -> None:
        its_params = list(ITS_params.pr_exp_params.values()) + list(ITS_params.pr_res_params.values())
        self.dataset['fExpClSizeITSHad'] = expected_cluster_size(self.dataset['fBetaGammaHad'], its_params)
        self.dataset['fSigmaClSizeCosLHad'] = sigma_its(self.dataset['fBetaGammaHad'], its_params)
        self.dataset.eval('fNSigmaITSHad = (fClSizeITSCosLamHad - fExpClSizeITSHad) / fSigmaClSizeCosLHad', inplace=True)
        self.dataset.drop(columns=['fExpClSizeITSHad', 'fSigmaClSizeCosLHad'], inplace=True)
    
    def define_nsigmaITS_He3(self) -> None:
        its_params = list(ITS_params.he_exp_params.values()) + list(ITS_params.he_res_params.values())
        self.dataset['fExpClSizeITSHe3'] = expected_cluster_size(self.dataset['fBetaGammaHe3'], its_params)
        self.dataset['fSigmaClSizeCosLHe3'] = sigma_its(self.dataset['fBetaGammaHe3'], its_params)
        self.dataset.eval('fNSigmaITSHe3 = (fClSizeITSCosLamHe3 - fExpClSizeITSHe3) / fSigmaClSizeCosLHe3', inplace=True)
        self.dataset.drop(columns=['fExpClSizeITSHe3', 'fSigmaClSizeCosLHe3'], inplace=True)

    # CUTS - DEPRECATED
    
    def general_selections(self) -> None:
        '''
            Selections for all particles.
            Needed columns are those of the original tree.
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying general selections')
        
        # cut in pseudorapidity
        self.dataset.query('-0.9 < fEtaHe3 < 0.9', inplace=True)
        self.dataset.query('-0.9 < fEtaHad < 0.9', inplace=True)

    def selections_He3(self) -> None:
        '''
            Selections for He3
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on He3')
        #self.dataset.query('fPIDtrkHe3 == 7', inplace=True)
        self.dataset.query('fPtHe3 < 2.5 or fPIDtrkHe3 == 7', inplace=True)
        
        #self.dataset.query('abs(fPtHe3) > 1.6', inplace=True)
        #self.dataset.query('0.5 < fChi2TPCHe3 < 4', inplace=True)
        self.dataset.query('fChi2TPCHe3 < 4', inplace=True)
        self.dataset.query('-2 < fNSigmaTPCHe3 < 2', inplace=True)
        self.dataset.query('abs(fDCAxyHe3) < 0.1', inplace=True)
        self.dataset.query('abs(fDCAzHe3) < 1.0', inplace=True)
        
        self.define_nsigmaITS_He3()
        #self.dataset.query('fNSigmaITSHe3 > -1.5', inplace=True)

    def selections_Pr(self) -> None:
        '''
            Selections for Pr
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on Had')
        #self.dataset.query('0.5 < fChi2TPCHad < 4', inplace=True)
        self.dataset.query('fChi2TPCHad < 4', inplace=True)
        self.dataset.query('abs(fNSigmaTPCHad) < 2', inplace=True)

        self.define_nsigmaTOF_Pr()
        self.dataset.query('(fPtHad < 0.8) or (-2 < fNSigmaTOFHad < 2)', inplace=True)
        
        self.define_nsigmaITS_Pr() 
        #self.dataset.query('fNSigmaITSHad > -1.5', inplace=True)

    def run_selections(self) -> None:
        '''
            Run selections on both He3 and Pr
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on He3')
        #self.dataset.query('fPIDtrkHe3 == 7', inplace=True)
        self.dataset.query('fPtHe3 < 2.5 or fPIDtrkHe3 == 7', inplace=True)
        
        #self.dataset.query('abs(fPtHe3) > 1.6', inplace=True)
        self.dataset.query('0.5 < fChi2TPCHe3 < 4', inplace=True)
        #self.dataset.query('fChi2TPCHe3 < 4', inplace=True)
        self.dataset.query('-2 < fNSigmaTPCHe3 < 2', inplace=True)
        self.dataset.query('abs(fDCAxyHe3) < 0.1', inplace=True)
        self.dataset.query('abs(fDCAzHe3) < 1.0', inplace=True)
        
        self.define_nsigmaITS_He3()
        self.dataset.query('fNSigmaITSHe3 > -1.5', inplace=True)


        ############################
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on Had')
        self.dataset.query('0.5 < fChi2TPCHad < 4', inplace=True)
        #self.dataset.query('fChi2TPCHad < 4', inplace=True)
        self.dataset.query('abs(fNSigmaTPCHad) < 2', inplace=True)

        self.define_nsigmaTOF_Pr()
        self.dataset.query('(fPtHad < 0.8) or (-2 < fNSigmaTOFHad < 2)', inplace=True)
        
        self.define_nsigmaITS_Pr()        
        self.dataset.query('fNSigmaITSHad > -1.5', inplace=True)

    # CUTS

    def apply_cut(self, col_exp:str) -> None:
        '''
            Apply a cut on the dataset
            Parameters:
            - col_exp (str): column: expression to be evaluated (e.g. 'fPtHe3:fPtHe3 > 1.6')
        '''
        column, expression = col_exp.split(':')
        if column not in self.dataset.columns:
            if column == 'fNSigmaITSHad':   self.define_nsigmaITS_Pr()
            elif column == 'fNSigmaITSHe3': self.define_nsigmaITS_He3()
            elif column == 'fNSigmaTOFHad': self.define_nsigmaTOF_Pr()
            else:
                print(tc.MAGENTA+'[WARNING]:'+tc.RESET+' Column',column,'not present in dataset!')
                return
        self.dataset.query(expression, inplace=True)

    # OTHERS

    @abstractmethod
    def correct_pt_H3_hp(self) -> None: 
        '''
            Corrected pT for He3 identified as H3 in tracking.
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Correcting pT for He3 identified as H3 in tracking')
        print(PID_params.he_params)
        self.dataset.loc[(self.dataset['fPIDtrkHe3'] == 6) & (self.dataset['fPtHe3'] < PID_params.he_params['ptmax']), 
                         'fPtHe3'] = self.dataset['fPtHe3'] - self.dataset['fPtHe3']*(PID_params.he_params['kp1'] + PID_params.he_params['kp2'] * self.dataset['fPtHe3'])
        #old code
        #self.dataset.loc[self.dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset['fPtHe3'] * (1 + PID_params.he_params['kp0'] + PID_params.he_params['kp1'] * self.dataset['fPtHe3'])

    def close_pair_selection(self) -> None:
        '''
            Close pair rejection
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying close pair rejection')
        sigma_delta_eta = 0.008
        sigma_delta_phi = 0.017
        self.dataset.query(f'(fDeltaPhi/(3*{sigma_delta_phi}))**2+ (fDeltaEta/(3*{sigma_delta_eta}))**2 > 1 ', inplace=True)

class DataPreprocessor(Preprocessor):

    def __init__(self, dataset: Dataset) -> None:
        '''
            Method inheritated from Processor class. Data is stored in the 'full' key.
            'reco' key is not used.
        '''
        super().__init__(dataset)
        self.dataset.add_subset('full', self.dataset['fPtHe3'] > -999.) # dummy check

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
                if opt not in self.available_subsets:
                    print(tc.MAGENTA+'[WARNING]:'+tc.RESET+' Subset',opt,'not available!')
                    continue
                #if 'cent' in opt:
                #    print(tc.RED+'[ERROR]:'+tc.RESET+' Centrality selection not implemented yet!')
                #    continue
                    
                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset.columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable']+part,'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.dataset.build_th1(cfg['xVariable']+part, axis_spec_x, subset=opt)
                    
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

                    hist = self.dataset.build_th2(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y, subset=opt)
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
        self.dataset.eval('fSignedPtMCHad = fPtMCHad', inplace=True)
        self.dataset.eval('fSignedPtMCLi = fPtMCLi', inplace=True)

        self.dataset.eval('fPtMCHe3 = abs(fPtMCHe3)', inplace=True)
        self.dataset.eval('fPtMCHad = abs(fPtMCHad)', inplace=True)
        self.dataset.eval('fPtMCLi = abs(fPtMCLi)', inplace=True)

        self.dataset.eval('fPMCHe3 = fPtMCHe3 * cosh(fEtaMCHe3)', inplace=True)
        self.dataset.eval('fPMCHad = fPtMCHad * cosh(fEtaMCHad)', inplace=True)
        #self.dataset['fPMCLi'] = self.dataset['fPtMCLi'] * np.cosh(self.dataset['fEtaLi'])

        self.dataset.eval('fDeltaEtaMC = fEtaMCHe3 - fEtaMCHad', inplace=True)
        self.dataset.eval('fDeltaPhiMC = fPhiMCHe3 - fPhiMCHad', inplace=True)

        # momentum resolution
        self.define_resolution()

        ## reconstructed variables
        self.dataset.eval(f'fEMCHe3 = sqrt(fPMCHe3**2 + {ParticleMasses["He"]}**2)', inplace=True)
        self.dataset.eval(f'fEMCHad = sqrt(fPMCHad**2 + {ParticleMasses["Had"]}**2)', inplace=True)
        self.dataset.eval('fMassInvMCLi = sqrt((fEMCHe3 + fEMCHad)**2 - (fPtMCHe3 * cos(fPhiMCHe3) + fPtMCHad * cos(fPhiMCHad))**2 - (fPtMCHe3 * sin(fPhiMCHe3) + fPtMCHad * sin(fPhiMCHad))**2 - (fPtMCHe3 * sinh(fEtaMCHe3) + fPtMCHad * sinh(fEtaMCHad))**2 )', inplace=True)

    def define_resolution(self) -> None:
        '''
            Definition of resolution in the dataframe
        '''

        ## pt resolution
        self.dataset.eval('fPtResHe3 = (fPtMCHe3 - fPtHe3) / fPtMCHe3', inplace=True)
        self.dataset.eval('fPtResHad = (fPtMCHad - fPtHad) / fPtMCHad', inplace=True)
        self.dataset.eval('fPtResLi = (fPtMCLi - fPtLi) / fPtMCLi', inplace=True)
        self.dataset.eval('fPtResNotNormHe3 = (fPtMCHe3 - fPtHe3)', inplace=True)
        self.dataset.eval('fPtResNotNormHad = (fPtMCHad - fPtHad)', inplace=True)
        self.dataset.eval('fPtResNotNormLi = (fPtMCLi - fPtLi)', inplace=True)

        ## p resolution
        self.dataset.eval('fPResHe3 = (fPMCHe3 - fInnerPTPCHe3) / fPMCHe3', inplace=True)
        self.dataset.eval('fPResHad = (fPMCHad - fInnerPTPCHad) / fPMCHad', inplace=True)

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
                if opt not in self.available_subsets:
                    print(tc.RED+'[WARNING]:'+tc.RESET+' Subset',opt,'not available!')
                    continue

                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue

                    axis_spec_x = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    hist = self.dataset.build_th1(cfg['xVariable']+part, axis_spec_x, subset=opt)
                    
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
                    hist = self.dataset.build_th2(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y, subset=opt)

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