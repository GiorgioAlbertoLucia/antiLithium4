'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import ctypes
import numpy as np
import polars as pl
import yaml
from abc import ABC, abstractmethod

import uproot
from ROOT import TFile
from ROOT.Math import PtEtaPhiMVector, Boost

import sys
sys.path.append('..')
sys.path.append('../..')

from torchic import AxisSpec, PlDataset
from torchic.physics.ITS import average_cluster_size, sigma_its, expected_cluster_size
from torchic.physics import py_BetheBloch
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

from utils.particles import ParticleMasses, ParticlePID, ParticleLabels
from utils.parametrisations import ITS as ITS_params
from utils.parametrisations import TPC as TPC_params
from utils.parametrisations import TOF as TOF_params
from utils.parametrisations import PIDforTracking as PID_params
PIDlabels = {int(ParticlePID[key]): ParticleLabels[key] for key in ParticlePID.keys()}

class Preprocessor(ABC):

    def __init__(self, dataset: PlDataset) -> None:
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
        
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Defining variables')
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

        ## definition of reconstructed variables
        self.dataset.with_columns( pl.col('fPtHe3').alias('fSignedPtHe3'))
        self.dataset.with_columns( pl.col('fPtHad').alias('fSignedPtHad'))
        self.dataset.with_columns( abs(pl.col('fPtHe3')).alias('fPtHe3'))
        self.dataset.with_columns( abs(pl.col('fPtHad')).alias('fPtHad'))
        self.dataset.with_columns( (pl.col('fSignedPtHe3')/pl.col('fPtHe3')).alias('fSignHe3'))
        self.dataset.with_columns( (pl.col('fSignedPtHad')/pl.col('fPtHad')).alias('fSignHad'))

        #if 'fIsBkgUS' not in self.dataset.columns:
        #    self.dataset['fIsBkgUS'] = True
        #    self.dataset.loc[self.dataset['fSignHe3'] == self.dataset['fSignHad'], 'fIsBkgUS'] = False

        self.correct_pt_H3_hp()
        self.dataset.with_columns((pl.col('fPtHe3') * pl.col('fSignHe3')).alias('fSignedPtHe3'))

        self.dataset.with_columns(np.sqrt((pl.col('fPtHe3') * np.cosh(pl.col('fEtaHe3')))**2 + ParticleMasses["He"]**2).alias(f'fEHe3'))
        self.dataset.with_columns(np.sqrt((pl.col('fPtHad') * np.cosh(pl.col('fEtaHad')))**2 + ParticleMasses["Pr"]**2).alias(f'fEHad'))
        self.dataset.with_columns((pl.col('fEtaHe3') - pl.col('fEtaHad')).alias('fDeltaEta'))
        self.dataset.with_columns((pl.col('fPhiHe3') - pl.col('fPhiHad')).alias('fDeltaPhi'))

        self.dataset.with_columns((pl.col('fInnerParamTPCHe3') * 2).alias('fInnerParamTPCHe3'))
        self.dataset.with_columns((pl.col('fInnerParamTPCHe3') * pl.col('fSignHe3')).alias('InnerParamTPCHe3'))
        self.dataset.with_columns((pl.col('fInnerParamTPCHad') * pl.col('fSignHad')).alias('InnerParamTPCHad'))
        
        fClSizeITSMeanHe3, fNHitsITSHe3 = average_cluster_size(self.dataset['fItsClusterSizeHe3'])
        self.dataset.with_columns(pl.Series(values=fClSizeITSMeanHe3, name='fClSizeITSMeanHe3'))
        self.dataset.with_columns(pl.Series(values=fNHitsITSHe3, name='fNHitsITSHe3'))
        fClSizeITSMeanHad, fNHitsITSHad = average_cluster_size(self.dataset['fItsClusterSizeHad'])
        self.dataset.with_columns(pl.Series(values=fClSizeITSMeanHad, name='fClSizeITSMeanHad'))
        self.dataset.with_columns(pl.Series(values=fNHitsITSHad, name='fNHitsITSHad'))
        self.dataset.with_columns((pl.col('fClSizeITSMeanHe3') / np.cosh(pl.col('fEtaHe3'))).alias('fClSizeITSCosLamHe3'))
        self.dataset.with_columns((pl.col('fClSizeITSMeanHad') / np.cosh(pl.col('fEtaHad'))).alias('fClSizeITSCosLamHad'))

        # beta*gamma (for Bethe-Bloch formula)
        self.dataset.with_columns((abs(pl.col('fInnerParamTPCHe3')) / ParticleMasses["He"]).alias(f'fBetaGammaHe3'))
        self.dataset.with_columns((abs(pl.col('fInnerParamTPCHad')) / ParticleMasses["Pr"]).alias(f'fBetaGammaHad'))

        # invariant mass 
        self.dataset.with_columns((pl.col('fPtHe3') * np.cos(pl.col('fPhiHe3')) + pl.col('fPtHad') * np.cos(pl.col('fPhiHad'))).alias('fPxLi'))
        self.dataset.with_columns((pl.col('fPtHe3') * np.sin(pl.col('fPhiHe3')) + pl.col('fPtHad') * np.sin(pl.col('fPhiHad'))).alias('fPyLi'))
        self.dataset.with_columns((pl.col('fPtHe3') * np.sinh(pl.col('fEtaHe3')) + pl.col('fPtHad') * np.sinh(pl.col('fEtaHad'))).alias('fPzLi'))
        self.dataset.with_columns((pl.col('fEHe3') + pl.col('fEHad')).alias('fELi'))
        self.dataset.with_columns(np.sqrt(pl.col('fPxLi')**2 + pl.col('fPyLi')**2 + pl.col('fPzLi')**2).alias('fPLi'))

        self.dataset.with_columns(np.sqrt(pl.col('fPxLi')**2 + pl.col('fPyLi')**2).alias('fPtLi'))
        self.dataset.with_columns((pl.col('fPtLi') / 2).alias('fKt'))
        self.dataset.with_columns(np.arccosh(pl.col('fPLi') / pl.col('fELi')).alias('fEtaLi'))
        self.dataset.with_columns(np.arctan2(pl.col('fPyLi'), pl.col('fPxLi')).alias('fPhiLi'))
        self.dataset.with_columns((pl.col('fPtLi') * pl.col('fSignHe3')).alias('fSignedPtLi'))
        self.dataset.with_columns(np.sqrt(pl.col('fELi')**2 - pl.col('fPLi')**2).alias('fMassInvLi'))
        self.dataset.with_columns(np.sqrt(pl.col('fELi')**2 - pl.col('fPtLi')**2).alias('fMassTLi'))

        #self.dataset.query('fMassInvLi > 3.755', inplace=True)

        # remove unnecessary columns
        self.dataset.drop(columns=[
                                   'fEHe3', 'fInnerParamTPCHe3', 'fClSizeITSMeanHe3', 'fItsClusterSizeHe3',
                                   'fEHad', 'fInnerParamTPCHad', 'fClSizeITSMeanHad', 'fItsClusterSizeHad',
                                   'fPxLi', 'fPyLi', 'fPzLi', 'fELi', 'fPLi', 
                                   #'fPtLi', 
                                   'fEtaLi', 'fPhiLi',
                                   ])

        # separate matter and antimatter
        self.add_subset('antimatter', pl.col('fSignHe3') < 0)
        self.add_subset('matter', pl.col('fSignHe3') > 0)
        
        #if 'fCentralityFT0C' in self.dataset.columns:
        #    self.add_subset('antimatter-cent0-10', (self.dataset['fSignHe3'] < 0) & (0 < self.dataset['fCentralityFT0C']) & (self.dataset['fCentralityFT0C'] < 10))
        #    self.add_subset('matter-cent0-10', (self.dataset['fSignHe3'] > 0) & (0 < self.dataset['fCentralityFT0C']) & (self.dataset['fCentralityFT0C'] < 10))
        
        if 'fIsBkgUS' in self.dataset.columns:
            #self.dataset['fIsBkgUSfloat'] = self.dataset['fIsBkgUS'].astype(float)
            self.add_subset('unlike-sign', pl.col('fIsBkgUS') == True)
            self.add_subset('like-sign', pl.col('fIsBkgUS') == False)
            if 'fCentralityFT0C' in self.dataset.columns:
                self.add_subset('unlike-sign-cent0-10', (pl.col('fIsBkgUSfloat') == 1) & (0 < pl.col('fCentralityFT0C')) & (pl.col('fCentralityFT0C') < 10))
        
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
    
    @timeit
    def define_kstar(self):
        '''
            Kstar: variale used for the study of correlation of p-He3
        '''
        
        np_compute_kstar = np.vectorize(self.compute_kstar)
        self.dataset.query('fPtHad < 10 and fPtHe3 < 10', inplace=True)
        fKstar = np_compute_kstar(self.dataset['fPtHe3'], self.dataset['fEtaHe3'], self.dataset['fPhiHe3'], ParticleMasses['He'], 
                                                  self.dataset['fPtHad'], self.dataset['fEtaHad'], self.dataset['fPhiHad'], ParticleMasses['Pr'])
        self.dataset.with_columns(pl.Series(values=fKstar, name='fKstar'))
        self.add_subset('anti-kstar0-200', (pl.col('fKstar') < 0.2) & (pl.col('fSignHe3') < 0))
        self.add_subset('matter-kstar0-200', (pl.col('fKstar') < 0.2) & (pl.col('fSignHe3') > 0))
    
    def define_nsigmaTOF_Pr(self) -> None:
        self.dataset.with_columns(pl.Series(values=np.ones(self.dataset.shape[0])*ParticleMasses['Pr'], name='fExpTOFMassHad'))
        self.dataset.with_columns(((TOF_params.pr_res_params['res1'] * np.exp(TOF_params.pr_res_params['res2'] * np.abs(pl.col('fPtHad')))) * pl.col('fExpTOFMassHad')).alias('fSigmaTOFMassHad'))
        self.dataset.with_columns(((pl.col('fMassTOFHad') - pl.col('fExpTOFMassHad')) / pl.col('fSigmaTOFMassHad')).alias('fNSigmaTOFHad'))
        self.dataset.drop(columns=['fExpTOFMassHad', 'fSigmaTOFMassHad'])
    
    def define_nsigmaITS_Pr(self) -> None:
        its_params = list(ITS_params.pr_exp_params.values()) + list(ITS_params.pr_res_params.values())

        fExpClSizeITSHad = expected_cluster_size(self.dataset['fBetaGammaHad'], its_params)
        self.dataset.with_columns(pl.Series(values=fExpClSizeITSHad, name='fExpClSizeITSHad'))

        fSigmaClSizeCosLHad = sigma_its(self.dataset['fBetaGammaHad'], its_params)
        self.dataset.with_columns(pl.Series(values=fSigmaClSizeCosLHad, name='fSigmaClSizeCosLHad'))

        self.dataset.with_columns((( pl.col('fClSizeITSCosLamHad') - pl.col('fExpClSizeITSHad')) / pl.col('fSigmaClSizeCosLHad')).alias('fNSigmaITSHad'))
        self.dataset.drop(columns=['fExpClSizeITSHad', 'fSigmaClSizeCosLHad'])
    
    def define_nsigmaITS_He3(self) -> None:
        its_params = list(ITS_params.he_exp_params.values()) + list(ITS_params.he_res_params.values())

        fExpClSizeITSHe3 = expected_cluster_size(self.dataset['fBetaGammaHe3'], its_params)
        self.dataset.with_columns(pl.Series(values=fExpClSizeITSHe3, name='fExpClSizeITSHe3'))

        fSigmaClSizeCosLHe3 = sigma_its(self.dataset['fBetaGammaHe3'], its_params)
        self.dataset.with_columns(pl.Series(values=fSigmaClSizeCosLHe3, name='fSigmaClSizeCosLHe3'))

        self.dataset.with_columns(((pl.col('fClSizeITSCosLamHe3') - pl.col('fExpClSizeITSHe3')) / pl.col('fSigmaClSizeCosLHe3')).alias('fNSigmaITSHe3'))
        self.dataset.drop(columns=['fExpClSizeITSHe3', 'fSigmaClSizeCosLHe3'])

    def define_nsigmaTPC_He3(self) -> None:
        np_bethe_bloch = np.vectorize(py_BetheBloch)
        fExpTPCSignalHe3 = np_bethe_bloch(self.dataset['fBetaGammaHe3'], *TPC_params.he_exp_params.values())
        self.dataset.with_columns(pl.Series(values=fExpTPCSignalHe3, name='fExpTPCSignalHe3'))
        fSigmaTPCHe3 = TPC_params.he_res_params['res'] * self.dataset['fExpTPCSignalHe3']
        self.dataset.with_columns(pl.Series(values=fSigmaTPCHe3, name='fSigmaTPCHe3'))
        self.dataset.with_columns(((pl.col('fSignalTPCHe3') - pl.col('fExpTPCSignalHe3')) / pl.col('fSigmaTPCHe3')).alias('fNSigmaTPCHe3'))
        self.dataset.drop(columns=['fExpTPCSignalHe3', 'fSigmaTPCHe3'])

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
        self.dataset.with_columns(
            pl.when((pl.col("fPIDtrkHe3") == 6) & (pl.col("fPtHe3") < PID_params.he_params["ptmax"]))
              .then(
                  pl.col("fPtHe3") - pl.col("fPtHe3") * (
                      PID_params.he_params["kp1"] + PID_params.he_params["kp2"] * pl.col("fPtHe3")
                  )
              )
              .otherwise(pl.col("fPtHe3"))
              .alias("fPtHe3")
)

    def close_pair_selection(self) -> None:
        '''
            Close pair rejection
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying close pair rejection')
        sigma_delta_eta = 0.008
        sigma_delta_phi = 0.017
        self.dataset.query(f'(fDeltaPhi/(3*{sigma_delta_phi}))**2+ (fDeltaEta/(3*{sigma_delta_eta}))**2 > 1 ', inplace=True)

class DataPreprocessor(Preprocessor):

    def __init__(self, dataset: PlDataset) -> None:
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

    def visualize_boost(self, outputFilePath, config) -> None:
        ''' 
            Visualization of data.
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Visualizing')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Creating output file '+tc.UNDERLINE+tc.CYAN+outputFilePath+tc.RESET)
        out_file = uproot.recreate(outputFilePath)

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
                    hist = self.dataset.build_boost1(cfg['xVariable']+part, axis_spec_x, subset=opt)
                    
                    hist_name = None
                    if cfg['dir'] != 'None':    hist_name = f'{cfg["dir"]}/'
                    else:                       hist_name = f''
                    out_file[hist_name+axis_spec_x.name] = hist

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
                    hist = self.dataset.build_boost2d(cfg['xVariable']+part, cfg['yVariable']+part, axis_spec_x, axis_spec_y, subset=opt)

                    hist_name = None
                    if cfg['dir'] != 'None':    hist_name = f'{cfg["dir"]}/'
                    else:                       hist_name = f''
                    out_file[hist_name+axis_spec_x.name] = hist
        
        out_file.Close()

    def save_df(self, outputFilePath: str, columns: list) -> None:
        '''
            Save the dataframe to a ROOT file.
            Parameters:
            - outputFilePath (str): path to the output file
            - columns (list): list of columns to be saved
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Saving dataframe to',outputFilePath)
        if 'fIsMatter' not in self.dataset.columns:
            self.dataset['fIsMatter'] = self.dataset['fSignHe3'] > 0
        with uproot.recreate(outputFilePath) as outfile:
            outfile['outTree'] = self.dataset[columns]

class MCPreprocessor(Preprocessor):

    def __init__(self, dataset: PlDataset) -> None:
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