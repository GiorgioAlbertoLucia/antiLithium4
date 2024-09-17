'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import numpy as np
import yaml
from abc import ABC, abstractmethod

from ROOT import TFile, TLorentzVector
from ROOT.Math import Boost

import sys
sys.path.append('..')
sys.path.append('../..')

from framework.src.axis_spec import AxisSpec
from framework.src.hist_handler import *
from framework.src.data_handler import *

from utils.particles import ParticleMasses, ParticlePID, ParticleLabels
PIDlabels = {int(ParticlePID[key]): ParticleLabels[key] for key in ParticlePID.keys()}

class Preprocessor(ABC):

    def __init__(self, dataHandler: DataHandler) -> None:
        '''
            - dataset (pd.DataFrame):  data to be preprocessed
            - recoDataset (pd.DataFrame): reconstructed tracks
            - dataset['full'] (pd.DataFrame): generated and reconstructed particles

        '''
        
        print()
        self.dataHandler = dataHandler
        tmpDataset = dataHandler.inData.copy()

        self.dataset = {'full': tmpDataset,
                        'reco': tmpDataset}

    # PUBLIC METHODS
    
    @abstractmethod
    def defineVariables(self) -> None:
        
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Defining variables')

        ## definition of reconstructed variables
        self.dataset['full']['fSignedPtHe3'] = self.dataset['full']['fPtHe3']
        self.dataset['full']['fSignedPtPr'] = self.dataset['full']['fPtPr']
        self.dataset['full']['fPtHe3'] = abs(self.dataset['full']['fPtHe3'])
        self.dataset['full']['fPtPr'] = abs(self.dataset['full']['fPtPr'])
        self.dataset['full']['fSignHe3'] = self.dataset['full']['fSignedPtHe3'] / self.dataset['full']['fPtHe3']
        self.dataset['full']['fSignPr'] = self.dataset['full']['fSignedPtPr'] / self.dataset['full']['fPtPr']

        self.dataset['full']['fPHe3'] = self.dataset['full']['fPtHe3'] * np.cosh(self.dataset['full']['fEtaHe3'])
        self.dataset['full']['fPPr'] = self.dataset['full']['fPtPr'] * np.cosh(self.dataset['full']['fEtaPr'])
        self.dataset['full']['fEHe3'] = np.sqrt(self.dataset['full']['fPHe3']**2 + ParticleMasses['He']**2)
        self.dataset['full']['fEPr'] = np.sqrt(self.dataset['full']['fPPr']**2 + ParticleMasses['Pr']**2)
        self.dataset['full']['fAlpha'] = self.dataset['full']['fPhiHe3'] - self.dataset['full']['fPhiPr']  # separation angle
        self.dataset['full']['fDeltaEta'] = self.dataset['full']['fEtaHe3'] - self.dataset['full']['fEtaPr']
        self.dataset['full']['fDeltaPhi'] = self.dataset['full']['fPhiHe3'] - self.dataset['full']['fPhiPr']

        self.dataset['full']['fInnerPTPCHe3'] = self.dataset['full']['fInnerParamTPCHe3'] * 2
        self.dataset['full']['fInnerPTPCPr'] = self.dataset['full']['fInnerParamTPCPr']

        self.dataset['full']['fSignedPTPCHe3'] = self.dataset['full']['fInnerPTPCHe3'] * self.dataset['full']['fSignHe3']
        self.dataset['full']['fSignedPTPCPr'] = self.dataset['full']['fInnerPTPCPr'] * self.dataset['full']['fSignPr']
        
        # ITS cluster size 
        for layer in range(7):
            read4Bits = lambda x, layer: (x >> layer*4) & 0b1111
            self.dataset['full'][f'fClSizeITS{layer}He3'] = self.dataset['full']['fItsClusterSizeHe3'].apply(read4Bits, args=(layer,))
            self.dataset['full'][f'fClSizeITS{layer}Pr'] = self.dataset['full']['fItsClusterSizePr'].apply(read4Bits, args=(layer,))
        
        self.dataset['full']['fCosLambdaHe3'] = 1/np.cosh(self.dataset['full']['fEtaHe3'])
        self.dataset['full']['fCosLambdaPr'] = 1/np.cosh(self.dataset['full']['fEtaPr'])
        
        self.dataset['full'].eval('fNHitsITSHe3 = (fClSizeITS0He3 > 0) * 1. + (fClSizeITS1He3 > 0) * 1. + (fClSizeITS2He3 > 0) * 1. + (fClSizeITS3He3 > 0) * 1. + (fClSizeITS4He3 > 0) * 1. + (fClSizeITS5He3 > 0) * 1. + (fClSizeITS6He3 > 0) * 1.', inplace=True)
        self.dataset['full'].eval('fNHitsITSPr = (fClSizeITS0Pr > 0) * 1. + (fClSizeITS1Pr > 0) * 1. + (fClSizeITS2Pr > 0) * 1. + (fClSizeITS3Pr > 0) * 1. + (fClSizeITS4Pr > 0) * 1. + (fClSizeITS5Pr > 0) * 1. + (fClSizeITS6Pr > 0) * 1.', inplace=True)
        self.dataset['full'].eval('fClSizeITSMeanHe3 = (fClSizeITS0He3 + fClSizeITS1He3 + fClSizeITS2He3 + fClSizeITS3He3 + fClSizeITS4He3 + fClSizeITS5He3 + fClSizeITS6He3) / fNHitsITSHe3', inplace=True)
        self.dataset['full'].eval('fClSizeITSMeanPr = (fClSizeITS0Pr + fClSizeITS1Pr + fClSizeITS2Pr + fClSizeITS3Pr + fClSizeITS4Pr + fClSizeITS5Pr + fClSizeITS6Pr) / fNHitsITSPr', inplace=True)
        self.dataset['full'].eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 * fCosLambdaHe3', inplace=True)
        self.dataset['full'].eval('fClSizeITSCosLamPr = fClSizeITSMeanPr * fCosLambdaPr', inplace=True)

        # beta*gamma (for Bethe-Bloch formula)
        self.dataset['full']['fBetaGammaHe3'] = self.dataset['full']['fInnerParamTPCHe3'] / ParticleMasses['He'] * 2.
        self.dataset['full']['fBetaGammaPr'] = self.dataset['full']['fInnerParamTPCPr'] / ParticleMasses['Pr']

        # invariant mass 
        self.dataset['full']['fPtLi'] = np.sqrt( self.dataset['full']['fPtHe3']**2 + self.dataset['full']['fPtPr']**2 + 
                                        2*self.dataset['full']['fPtHe3']*self.dataset['full']['fPtPr']*np.cos(self.dataset['full']['fAlpha']) )
        self.dataset['full']['fSignedPtLi'] = self.dataset['full']['fPtLi'] * np.sign(self.dataset['full']['fPtHe3'])
        self.dataset['full']['fMassInvLi'] = np.sqrt((self.dataset['full']['fEHe3'] + self.dataset['full']['fEPr'])**2 - 
                                  (self.dataset['full']['fPtHe3'] * np.cos(self.dataset['full']['fPhiHe3']) + self.dataset['full']['fPtPr'] * np.cos(self.dataset['full']['fPhiPr']))**2 -
                                  (self.dataset['full']['fPtHe3'] * np.sin(self.dataset['full']['fPhiHe3']) + self.dataset['full']['fPtPr'] * np.sin(self.dataset['full']['fPhiPr']))**2 -
                                  (self.dataset['full']['fPtHe3'] * np.sinh(self.dataset['full']['fEtaHe3']) + self.dataset['full']['fPtPr'] * np.sinh(self.dataset['full']['fEtaPr']))**2 )
    
    def computeKstar(self, p1, eta1, phi1, e1, p2, eta2, phi2, e2):
        '''
            Kstar: variable used to study the correlation between two particles. It is computed as half of the difference between the 
            momentum of the two particles in the center of mass reference frame.

            NOTE: this is not vectorized

            Parameters:
            - p1: modulus of the 3-momentum of particle 1
            - eta1: pseudorapidity of particle 1
            - phi1: angle on the transverse plane for particle 1
            - e1: energy for particle 1
            - p2: modulus of the 3-momentum of particle 2
            - eta2: pseudorapidity of particle 2
            - phi2: angle on the transverse plane for particle 2
            - e2: energy for particle 2
        '''

        p1x = p1/np.sinh(eta1)*np.cos(phi1)
        p1y = p1/np.sinh(eta1)*np.sin(phi1)
        p1z = p1/np.cosh(eta1)
        p2x = p2/np.sinh(eta2)*np.cos(phi2)
        p2y = p2/np.sinh(eta2)*np.sin(phi2)
        p2z = p2/np.cosh(eta2)

        p1mu = TLorentzVector(p1x, p1y, p1z, e1)
        p2mu = TLorentzVector(p2x, p2y, p2z, e2)
        Pmu = p1mu + p2mu

        beta = Pmu.Beta()
        betax = beta * np.cos(Pmu.Phi()) * np.sin(Pmu.Theta())
        betay = beta * np.sin(Pmu.Phi()) * np.sin(Pmu.Theta())
        betaz = beta * np.cos(Pmu.Theta())

        p1muStar, p2muStar = TLorentzVector(p1mu), TLorentzVector(p2mu)
        p1muStar.SetXYZM(p1x, p1y, p1z, p1mu.M())
        p2muStar.SetXYZM(p2x, p2y, p2z, p2mu.M())
        
        boost = Boost(-betax, -betay, -betaz)
        p1muStar = boost(p1muStar)
        p2muStar = boost(p2muStar)

        Kstarmu = 0.5 * (p1muStar - p2muStar)
        return Kstarmu.P()

    def defineKstar(self):
        '''
            Kstar: variale used for the study of correlation of p-He3
        '''
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Defining Kstar')
        vComputeKstar = np.vectorize(self.computeKstar)
        self.dataset['full']['fKstar'] = vComputeKstar(self.dataset['full']['fPHe3'], self.dataset['full']['fEtaHe3'], self.dataset['full']['fPhiHe3'], self.dataset['full']['fEHe3'],
                                                       self.dataset['full']['fPPr'], self.dataset['full']['fEtaPr'], self.dataset['full']['fPhiPr'], self.dataset['full']['fEPr'])
        #self.dataset['full']['fPtHe3'].apply(lambda x: self.computeKstar(x, self.dataset['full']['fEtaHe3'], self.dataset['full']['fPhiHe3'], self.dataset['full']['fEHe3'],
                                         #                                                               self.dataset['full']['fPtPr'], self.dataset['full']['fEtaPr'], self.dataset['full']['fPhiPr'], self.dataset['full']['fEPr']))
    
    @abstractmethod
    def correctPtH3hp(self) -> None: 
        '''
            Corrected pT for He3 identified as H3 in tracking.
        '''

        curveParams = {'kp0': -0.233625,
                       'kp1': 0.12484,
                       'kp2': -0.015673
                       }
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Correcting pT for He3 identified as H3 in tracking')
        print('Using pol2 correction')
        print('Parameters:', curveParams)

        # change values only to rows where fPIDtrkHe3 == 6
        # pol1 correction
        # self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] * (1 + curveParams['kp0'] + curveParams['kp1'] * self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'])
        # pol2 correction
        self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] + curveParams['kp0'] + curveParams['kp1'] * self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] + curveParams['kp2'] * self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3']**2
      
    def visualize(self, config) -> None:
        
        self.histHandler = None # reset histHandler
        self.histHandler = {'full': HistHandler.createInstance(self.dataset['full']),
                            'reco': HistHandler.createInstance(self.dataset['reco'])}

    def filterAntimatter(self) -> None:
        '''
            Filter out particles with negative charge.
        '''
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Selecting only particles with negative charge')
        self.dataset['full'].query('fSignHe3 == -1 and fSignPr == -1', inplace=True)

    def selectionsHe3(self) -> None:
        '''
            Selections for He3
        '''

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Applying selections on He3')
        self.dataset['full'] = self.dataset['full'].query('fClSizeITSCosLamHe3 > 4.7', inplace=False)
        self.dataset['full'] = self.dataset['full'].query('0.5 < fChi2TPCHe3 < 4', inplace=False)

    def selectionsPr(self) -> None:
        '''
            Selections for Pr
        '''

        BB_params = {
            'kp1': -0.031712,
            'kp2': -45.0275,
            'kp3': -0.997645,
            'kp4': 1.68228,
            'kp5': 0.0108484
        }

        def BBfunc(x):
            x = np.abs(x)
            beta = x / np.sqrt(1 + x**2)
            aa = beta**BB_params['kp4']
            bb = (1/x)**BB_params['kp5']
            if (BB_params['kp3'] + bb) < 0:
                return -1.
            bb = np.log(BB_params['kp3'] + bb)
            return (BB_params['kp2'] - aa - bb) * BB_params['kp1'] / aa
        
        self.dataset['full']['fExpClSizeITSCosLamPr'] = self.dataset['full']['fBetaGammaPr']
        self.dataset['full']['fExpClSizeITSCosLamPr'].apply(BBfunc)
        sigma_params = {
            'kp0': 0.418451,
            'kp1': -0.040885
        }
        self.dataset['full']['fSigmaClSizeCosLPr'] = self.dataset['full']['fExpClSizeITSCosLamPr'] * (sigma_params['kp0'] + sigma_params['kp1'] * self.dataset['full']['fExpClSizeITSCosLamPr'])
        self.dataset['full']['fNSigmaPr'] = (self.dataset['full']['fClSizeITSCosLamPr'] - self.dataset['full']['fExpClSizeITSCosLamPr']) / self.dataset['full']['fSigmaClSizeCosLPr']
        self.dataset['full'] = self.dataset['full'].query('fNSigmaPr > 0', inplace=False)
        
class DataPreprocessor(Preprocessor):

    def __init__(self, dataHandler: DataHandler) -> None:
        '''
            Method inheritated from Processor class. Data is stored in the 'full' key.
            'reco' key is not used.
        '''
        super().__init__(dataHandler)

    def defineVariables(self) -> None:
        '''
            Definition of variables for data.
        '''
        super().defineVariables()

    def correctPtH3hp(self) -> None:
        '''
            Corrected pT for He3 identified as H3 in tracking using parameters from a dedicated curve
            (see studies.py H3inTrkStudy class)
        '''

        super().correctPtH3hp()

    def visualize(self, config, output_suffix='') -> None:
        ''' 
            Visualization of data.
        '''

        super().visualize(config)
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Visualizing')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        outputFilePath = config['outputFilePath'].replace('.root', output_suffix+'.root')
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Creating output file '+tc.UNDERLINE+tc.CYAN+outputFilePath+tc.RESET)
        outFile = TFile(outputFilePath, 'recreate')
        outDirs = {}
        for dir in config['outDirs']:   outDirs[dir] = outFile.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outputFilePath':         continue
            if key == 'studiesOutputFilePath':  continue
            if key == 'outDirs':                continue
            
            for part in cfg['particle']:

                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset['full'].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    
                    hist = self.histHandler['full'].buildTH1(cfg['xVariable']+part, axisSpecX)
                    #if cfg['xVariable'] == 'fPIDtrk':   self.histHandler['full'].setLabels(hist, PIDlabels, 'x')
                    
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset['full'].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue
                    elif cfg['yVariable']+part not in self.dataset['full'].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['yVariable'],'not present in dataset!')
                        continue

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.histHandler['full'].buildTH2(cfg['xVariable']+part, cfg['yVariable']+part, axisSpecX, axisSpecY)
                    #if cfg['xVariable'] == 'fPIDtrk':   hist = self.histHandler['full'].setLabels(hist, PIDlabels, 'x')
                    #if cfg['yVariable'] == 'fPIDtrk':   hist = self.histHandler['full'].setLabels(hist, PIDlabels, 'y')

                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()

        
        outFile.Close()



class MCPreprocessor(Preprocessor):

    def __init__(self, dataHandler: DataHandler) -> None:
        '''
            Method inheritated from Processor class. Full sample is stored in the 'full' key.
            Reconstructed tracks are stored in the 'reco' key.
        '''

        super().__init__(dataHandler)

    def __updateRecoTracks(self) -> None:
            
        self.dataset['reco'] = self.dataset['full'].copy().query('fPIDtrkHe3 != 0xFFFFF', inplace=False)
    
    # PUBLIC METHODS

    def defineVariables(self) -> None:
        '''
            Definition of variables for MC data.
        '''
        super().defineVariables()

        ## MC variables
        self.dataset['full']['fPtMCLi'] = self.dataset['full']['fSignedPtMC'] 
        self.dataset['full']['fMassMCLi'] = self.dataset['full']['fMassMC']

        # pt MC stores sign
        self.dataset['full'].eval('fSignLi = (-1)**(fPtMCLi < 0)', inplace=True)

        self.dataset['full']['fSignedPtMCHe3'] = self.dataset['full']['fPtMCHe3']
        self.dataset['full']['fSignedPtMCPr'] = self.dataset['full']['fPtMCPr']
        self.dataset['full']['fSignedPtMCLi'] = self.dataset['full']['fPtMCLi']

        self.dataset['full']['fPtMCHe3'] = abs(self.dataset['full']['fPtMCHe3'])
        self.dataset['full']['fPtMCPr'] = abs(self.dataset['full']['fPtMCPr'])
        self.dataset['full']['fPtMCLi'] = abs(self.dataset['full']['fPtMCLi'])

        self.dataset['full']['fPMCHe3'] = self.dataset['full']['fPtMCHe3'] * np.cosh(self.dataset['full']['fEtaMCHe3'])
        self.dataset['full']['fPMCPr'] = self.dataset['full']['fPtMCPr'] * np.cosh(self.dataset['full']['fEtaMCPr'])
        #self.dataset['full']['fPMCLi'] = self.dataset['full']['fPtMCLi'] * np.cosh(self.dataset['full']['fEtaLi'])

        self.dataset['full']['fDeltaEtaMC'] = self.dataset['full']['fEtaMCHe3'] - self.dataset['full']['fEtaMCPr']
        self.dataset['full']['fDeltaPhiMC'] = self.dataset['full']['fPhiMCHe3'] - self.dataset['full']['fPhiMCPr']

        # momentum resolution
        self.defineResolution()

        ## reconstructed variables
        self.dataset['full']['fEMCHe3'] = np.sqrt(self.dataset['full']['fPMCHe3']**2 + ParticleMasses['He']**2)
        self.dataset['full']['fEMCPr'] = np.sqrt(self.dataset['full']['fPMCPr']**2 + ParticleMasses['Pr']**2)
        self.dataset['full']['fMassInvMCLi'] = np.sqrt((self.dataset['full']['fEMCHe3'] + self.dataset['full']['fEMCPr'])**2 - 
                                  (self.dataset['full']['fPtMCHe3'] * np.cos(self.dataset['full']['fPhiMCHe3']) + self.dataset['full']['fPtMCPr'] * np.cos(self.dataset['full']['fPhiMCPr']))**2 -
                                  (self.dataset['full']['fPtMCHe3'] * np.sin(self.dataset['full']['fPhiMCHe3']) + self.dataset['full']['fPtMCPr'] * np.sin(self.dataset['full']['fPhiMCPr']))**2 -
                                  (self.dataset['full']['fPtMCHe3'] * np.sinh(self.dataset['full']['fEtaMCHe3']) + self.dataset['full']['fPtMCPr'] * np.sinh(self.dataset['full']['fEtaMCPr']))**2 )

        self.__updateRecoTracks()

    def defineResolution(self) -> None:
        '''
            Definition of resolution in the dataframe
        '''

        ## pt resolution
        self.dataset['full']['fPtResHe3'] = (self.dataset['full']['fPtMCHe3'] - self.dataset['full']['fPtHe3']) / self.dataset['full']['fPtMCHe3']
        self.dataset['full']['fPtResPr'] = (self.dataset['full']['fPtMCPr'] - self.dataset['full']['fPtPr']) / self.dataset['full']['fPtMCPr']
        self.dataset['full']['fPtResLi'] = (self.dataset['full']['fPtMCLi'] - self.dataset['full']['fPtLi']) / self.dataset['full']['fPtMCLi']
        self.dataset['full']['fPtResNotNormHe3'] = (self.dataset['full']['fPtMCHe3'] - self.dataset['full']['fPtHe3'])
        self.dataset['full']['fPtResNotNormPr'] = (self.dataset['full']['fPtMCPr'] - self.dataset['full']['fPtPr'])
        self.dataset['full']['fPtResNotNormLi'] = (self.dataset['full']['fPtMCLi'] - self.dataset['full']['fPtLi'])

        ## p resolution
        self.dataset['full']['fPResHe3'] = (self.dataset['full']['fPMCHe3'] - self.dataset['full']['fInnerPTPCHe3']) / self.dataset['full']['fPMCHe3']
        self.dataset['full']['fPResPr'] = (self.dataset['full']['fPMCPr'] - self.dataset['full']['fInnerPTPCPr']) / self.dataset['full']['fPMCPr']
        #self.dataset['full']['fPResLi'] = (self.dataset['full']['fPMCLi'] - self.dataset['full']['fInnterPTPCLi']) / self.dataset['full']['fPMCLi']

    def defineKstar(self):
        super().defineKstar()
        self.__updateRecoTracks()

    def filterAntimatter(self) -> None:

        '''
            Filter out particles with negative charge.
        '''

        super().filterAntimatter()
        self.__updateRecoTracks()

    def correctPtH3hp(self) -> None:
        '''
            Corrected pT for He3 identified as H3 in tracking using parameters from a dedicated curve
            (see studies.py H3inTrkStudy class)
        '''

        super().correctPtH3hp()
        self.defineResolution()
        self.__updateRecoTracks()

    def visualize(self, config) -> None:
        ''' 
            Visualization of MC.
        '''

        super().visualize(config)

        print(tc.GREEN+'[INFO]:'+tc.RESET+' visualizing')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print(tc.GREEN+'[INFO]:'+tc.RESET+' creating output file '+tc.UNDERLINE+tc.CYAN+f'{config["outputFilePath"]}'+tc.RESET)
        outFile = TFile(config['outputFilePath'], 'recreate')
        outDirs = {}
        for dir in config['outDirs']:   outDirs[dir] = outFile.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outputFilePath':         continue
            if key == 'studiesOutputFilePath':  continue
            if key == 'outDirs':                continue
            
            for part in cfg['particle']:

                if 'TH1' in cfg['type']:

                    if cfg['xVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    
                    hist = self.histHandler[cfg['opt']].buildTH1(cfg['xVariable']+part, axisSpecX)
                    if cfg['xVariable'] == 'fPIDtrk':   self.histHandler[cfg['opt']].setLabels(hist, PIDlabels, 'x')
                    
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:
                    
                    if cfg['xVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['xVariable'],'not present in dataset!')
                        continue
                    elif cfg['yVariable']+part not in self.dataset[cfg['opt']].columns:
                        print(tc.MAGENTA+'[WARNING]:'+tc.RESET,cfg['yVariable'],'not present in dataset!')
                        continue

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.histHandler[cfg['opt']].buildTH2(cfg['xVariable']+part, cfg['yVariable']+part, axisSpecX, axisSpecY)
                    if cfg['xVariable'] == 'fPIDtrk':   self.histHandler[cfg['opt']].setLabels(hist, PIDlabels, 'x')
                    if cfg['yVariable'] == 'fPIDtrk':   self.histHandler[cfg['opt']].setLabels(hist, PIDlabels, 'y')

                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()

                if 'TEfficiency' in cfg['type']:

                    axisSpecNum = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+'Reco'+part, cfg['title']+f' {part}')
                    axisSpecDen = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    hNumerator = self.histHandler['reco'].buildTH1(cfg['numVariable']+part, axisSpecNum)
                    hDenominator = self.histHandler['full'].buildTH1(cfg['denVariable']+part, axisSpecDen)

                    eff = self.histHandler['reco'].buildEfficiency(hNumerator, hDenominator)
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()  
                    eff.Write()
        
        outFile.Close()