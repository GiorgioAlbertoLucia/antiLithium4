'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import pandas as pd
import numpy as np
import yaml
from abc import ABC, abstractmethod

from hipe4ml.tree_handler import TreeHandler

from ROOT import TFile, TH1F, TH2F, TEfficiency

import sys
sys.path.append('..')
sys.path.append('../..')

from framework.src.axisSpec import AxisSpec
from framework.src.histHandler import *
from framework.src.dataHandler import *

from utils.particle import particleMass, PIDlabels

class Preprocessor(ABC):

    def __init__(self, dataHandler: DataHandler) -> None:
        '''
            - dataset (pd.DataFrame):  data to be preprocessed
            - recoDataset (pd.DataFrame): reconstructed tracks
            - dataset['full'] (pd.DataFrame): generated particles

        '''
        
        self.dataHandler = dataHandler
        tmpDataset = dataHandler.inData.copy()

        self.dataset = {'full': tmpDataset,
                        'reco': tmpDataset}

        self.histHandler = {'full': HistHandler.createInstance(self.dataset['full']),
                            'reco': HistHandler.createInstance(self.dataset['reco'])}

    # PUBLIC METHODS
    
    @abstractmethod
    def defineVariables(self) -> None:
        
        print('\nDefining variables...')

        ## definition of reconstructed variables
        self.dataset['full'].eval('fSignHe3 = (-1)**(fPtHe3 < 0)', inplace=True)
        self.dataset['full'].eval('fSignPr = (-1)**(fPtPr < 0)', inplace=True)
        self.dataset['full']['fSignedPtHe3'] = self.dataset['full']['fPtHe3']
        self.dataset['full']['fSignedPtPr'] = self.dataset['full']['fPtPr']
        self.dataset['full']['fPtHe3'] = abs(self.dataset['full']['fPtHe3'])
        self.dataset['full']['fPtPr'] = abs(self.dataset['full']['fPtPr'])

        self.dataset['full']['fPHe3'] = self.dataset['full']['fPtHe3'] * np.cosh(self.dataset['full']['fEtaHe3'])
        self.dataset['full']['fPPr'] = self.dataset['full']['fPtPr'] * np.cosh(self.dataset['full']['fEtaPr'])
        self.dataset['full']['fEHe3'] = np.sqrt(self.dataset['full']['fPHe3']**2 + particleMass['He3']**2)
        self.dataset['full']['fEPr'] = np.sqrt(self.dataset['full']['fPPr']**2 + particleMass['Pr']**2)
        self.dataset['full']['fAlpha'] = self.dataset['full']['fPhiHe3'] - self.dataset['full']['fPhiPr']  # separation angle

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
        
        self.dataset['full'].eval('fClSizeITSMeanHe3 = (fClSizeITS0He3 + fClSizeITS1He3 + fClSizeITS2He3 + fClSizeITS3He3+ fClSizeITS4He3 + fClSizeITS5He3 + fClSizeITS6He3) / 7', inplace=True)
        self.dataset['full'].eval('fClSizeITSMeanPr = (fClSizeITS0Pr + fClSizeITS1Pr + fClSizeITS2Pr + fClSizeITS3Pr+ fClSizeITS4Pr + fClSizeITS5Pr + fClSizeITS6Pr) / 7', inplace=True)
        self.dataset['full'].eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 * fCosLambdaHe3', inplace=True)
        self.dataset['full'].eval('fClSizeITSCosLamPr = fClSizeITSMeanPr * fCosLambdaPr', inplace=True)

        # beta*gamma (for Bethe-Bloch formula)
        self.dataset['full']['fBetaGammaHe3'] = self.dataset['full']['fInnerParamTPCHe3'] / particleMass['He3'] * 2.
        self.dataset['full']['fBetaGammaPr'] = self.dataset['full']['fInnerParamTPCPr'] / particleMass['Pr']

        # invariant mass 
        self.dataset['full']['fPtLi'] = np.sqrt( self.dataset['full']['fPtHe3']**2 + self.dataset['full']['fPtPr']**2 + 
                                        2*self.dataset['full']['fPtHe3']*self.dataset['full']['fPtPr']*np.cos(self.dataset['full']['fAlpha']) )
        self.dataset['full']['fSignedPtLi'] = self.dataset['full']['fPtLi'] * np.sign(self.dataset['full']['fPtHe3'])
        self.dataset['full']['fMassInvLi'] = np.sqrt((self.dataset['full']['fEHe3'] + self.dataset['full']['fEPr'])**2 - 
                                  (self.dataset['full']['fPtHe3'] * np.cos(self.dataset['full']['fPhiHe3']) + self.dataset['full']['fPtPr'] * np.cos(self.dataset['full']['fPhiPr']))**2 -
                                  (self.dataset['full']['fPtHe3'] * np.sin(self.dataset['full']['fPhiHe3']) + self.dataset['full']['fPtPr'] * np.sin(self.dataset['full']['fPhiPr']))**2 -
                                  (self.dataset['full']['fPtHe3'] * np.sinh(self.dataset['full']['fEtaHe3']) + self.dataset['full']['fPtPr'] * np.sinh(self.dataset['full']['fEtaPr']))**2 )
    
    @abstractmethod
    def correctPtH3hp(self) -> None: 
        '''
            Corrected pT for He3 identified as H3 in tracking.
        '''

        curveParams = {'kp0': -0.87643,
                       'kp1': 0.79767,
                       'kp2': -0.19314
                       }

        # change values only to rows where fPIDtrkHe3 == 6
        # pol1 correction
        # self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] * (1 + curveParams['kp0'] + curveParams['kp1'] * self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'])
        # pol2 correction
        self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] = self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] * (1 + curveParams['kp0'] + curveParams['kp1'] * self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3'] + curveParams['kp2'] * self.dataset['full'].loc[self.dataset['full']['fPIDtrkHe3'] == 6, 'fPtHe3']**2)

        
    @abstractmethod
    def visualize(self, config) -> None:
        return NotImplemented


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

    def visualize(self, config) -> None:
        ''' 
            Visualization of data.
        '''

        print('\nVisualizing...')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print('Creating output file '+tc.UNDERLINE+tc.CYAN+f'{config["outputFilePath"]}'+tc.RESET+'...')
        outFile = TFile(config['outputFilePath'], 'recreate')
        outDirs = {}
        for dir in config['outDirs']:   outDirs[dir] = outFile.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outputFilePath':         continue
            if key == 'studiesOutputFilePath':  continue
            if key == 'outDirs':                continue
            
            for part in cfg['particle']:

                if 'TH1' in cfg['type']:

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    
                    hist = self.histHandler['full'].buildTH1(cfg['xVariable']+part, axisSpecX)
                    if cfg['xVariable'] == 'fPIDtrk':   self.histHandler['full'].setLabels(hist, PIDlabels, 'x')
                    
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+part, cfg['title']+f' {part}')

                    hist = self.histHandler[cfg['opt']].buildTH2(cfg['xVariable']+part, cfg['yVariable']+part, axisSpecX, axisSpecY)
                    if cfg['xVariable'] == 'fPIDtrk':   self.histHandler['full'].setLabels(hist, PIDlabels, 'x')
                    if cfg['yVariable'] == 'fPIDtrk':   self.histHandler['full'].setLabels(hist, PIDlabels, 'y')

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

        self.dataset['full']['fPtMCHe3'] = abs(self.dataset['full']['fPtMCHe3'])
        self.dataset['full']['fPtMCPr'] = abs(self.dataset['full']['fPtMCPr'])
        self.dataset['full']['fPtMCLi'] = abs(self.dataset['full']['fPtMCLi'])

        self.dataset['full']['fPMCHe3'] = self.dataset['full']['fPtMCHe3'] * np.cosh(self.dataset['full']['fEtaHe3'])
        self.dataset['full']['fPMCPr'] = self.dataset['full']['fPtMCPr'] * np.cosh(self.dataset['full']['fEtaPr'])
        #self.dataset['full']['fPMCLi'] = self.dataset['full']['fPtMCLi'] * np.cosh(self.dataset['full']['fEtaLi'])

        # momentum resolution
        self.defineResolution()

        ## reconstructed variables
        self.dataset['full']['fEMCHe3'] = np.sqrt(self.dataset['full']['fPMCHe3']**2 + particleMass['He3']**2)
        self.dataset['full']['fEMCPr'] = np.sqrt(self.dataset['full']['fPMCPr']**2 + particleMass['Pr']**2)

        self.dataset['full']['fMassInvMCLi'] = np.sqrt((self.dataset['full']['fEMCHe3'] + self.dataset['full']['fEMCPr'])**2 - 
                                  (self.dataset['full']['fPtMCHe3'] * np.cos(self.dataset['full']['fPhiHe3']) + self.dataset['full']['fPtMCPr'] * np.cos(self.dataset['full']['fPhiPr']))**2 -
                                  (self.dataset['full']['fPtMCHe3'] * np.sin(self.dataset['full']['fPhiHe3']) + self.dataset['full']['fPtMCPr'] * np.sin(self.dataset['full']['fPhiPr']))**2 -
                                  (self.dataset['full']['fPtMCHe3'] * np.sinh(self.dataset['full']['fEtaHe3']) + self.dataset['full']['fPtMCPr'] * np.sinh(self.dataset['full']['fEtaPr']))**2 )

        self.__updateRecoTracks()

    def defineResolution(self) -> None:
        '''
            Definition of resolution in the dataframe
        '''

        ## pt resolution
        self.dataset['full']['fPtResHe3'] = (self.dataset['full']['fPtMCHe3'] - self.dataset['full']['fPtHe3']) / self.dataset['full']['fPtMCHe3']
        self.dataset['full']['fPtResPr'] = (self.dataset['full']['fPtMCPr'] - self.dataset['full']['fPtPr']) / self.dataset['full']['fPtMCPr']
        self.dataset['full']['fPtResLi'] = (self.dataset['full']['fPtMCLi'] - self.dataset['full']['fPtLi']) / self.dataset['full']['fPtMCLi']

        ## pt resolution
        self.dataset['full']['fPResHe3'] = (self.dataset['full']['fPMCHe3'] - self.dataset['full']['fInnerPTPCHe3']) / self.dataset['full']['fPMCHe3']
        self.dataset['full']['fPResPr'] = (self.dataset['full']['fPMCPr'] - self.dataset['full']['fInnerPTPCPr']) / self.dataset['full']['fPMCPr']
        #self.dataset['full']['fPResLi'] = (self.dataset['full']['fPMCLi'] - self.dataset['full']['fInnterPTPCLi']) / self.dataset['full']['fPMCLi']

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

        print('\nVisualizing...')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print('Creating output file '+tc.UNDERLINE+tc.CYAN+f'{config["outputFilePath"]}'+tc.RESET+'...')
        outFile = TFile(config['outputFilePath'], 'recreate')
        outDirs = {}
        for dir in config['outDirs']:   outDirs[dir] = outFile.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outputFilePath':         continue
            if key == 'studiesOutputFilePath':  continue
            if key == 'outDirs':                continue
            
            for part in cfg['particle']:

                if 'TH1' in cfg['type']:

                    axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+part, cfg['title']+f' {part}')
                    
                    hist = self.histHandler[cfg['opt']].buildTH1(cfg['xVariable']+part, axisSpecX)
                    if cfg['xVariable'] == 'fPIDtrk':   self.histHandler[cfg['opt']].setLabels(hist, PIDlabels, 'x')
                    
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()

                if 'TH2' in cfg['type']:

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