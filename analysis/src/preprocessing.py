'''
    Class to read and preprocess data. A first visualization is also implemented.
'''

import pandas as pd
import numpy as np
import yaml

from hipe4ml.tree_handler import TreeHandler

from ROOT import TFile, TH1F, TH2F, TEfficiency

import sys
sys.path.append('..')

from analysisFramework.src.axisSpec import AxisSpec
from analysisFramework.src.histHandler import *
from analysisFramework.src.dataHandler import *

from utils.particle import particleMass, PIDlabels

class Preprocessor:

    def __init__(self, dataHandler: DataHandler) -> None:
        '''
            - dataset (pd.DataFrame):  data to be preprocessed
            - recoDataset (pd.DataFrame): reconstructed tracks
            - dataset['gen'] (pd.DataFrame): generated particles

        '''
        
        self.dataHandler = dataHandler
        __dataset = dataHandler.inData.copy()

        self.dataset = {'gen': __dataset,
                        'reco': __dataset}

        self.histHandler = {'gen': HistHandler.createInstance(self.dataset['gen']),
                            'reco': HistHandler.createInstance(self.dataset['reco'])}


    def __selectRecoTracks(self) -> None:
            
        self.dataset['reco'] = self.dataset['gen'].copy().query('fPIDtrkHe3 != 0xFFFFF', inplace=False)

    # public methods
    def defineVariables(self) -> None:

        self.dataset['gen']['fPtMCLi'] = self.dataset['gen']['fSignedPtMC'] 
        self.dataset['gen']['fMassMCLi'] = self.dataset['gen']['fMassMC']

        # pt MC stores sign
        self.dataset['gen'].eval('fSignHe3 = (-1)**(fPtMCHe3 < 0)', inplace=True)
        self.dataset['gen'].eval('fSignPr = (-1)**(fPtMCPr < 0)', inplace=True)
        self.dataset['gen'].eval('fSignLi = (-1)**(fPtMCLi < 0)', inplace=True)

        self.dataset['gen']['fSignedPtHe3'] = self.dataset['gen']['fPtHe3'] * self.dataset['gen']['fSignHe3']
        self.dataset['gen']['fSignedPtPr'] = self.dataset['gen']['fPtPr'] * self.dataset['gen']['fSignPr']

        self.dataset['gen']['fPtMCHe3'] = abs(self.dataset['gen']['fPtMCHe3'])
        self.dataset['gen']['fPtMCPr'] = abs(self.dataset['gen']['fPtMCPr'])
        self.dataset['gen']['fPtMCLi'] = abs(self.dataset['gen']['fPtMCLi'])

        self.dataset['gen']['fPMCHe3'] = self.dataset['gen']['fPtMCHe3'] * np.cosh(self.dataset['gen']['fEtaHe3'])
        self.dataset['gen']['fPMCPr'] = self.dataset['gen']['fPtMCPr'] * np.cosh(self.dataset['gen']['fEtaPr'])
        #self.dataset['gen']['fPMCLi'] = self.dataset['gen']['fPtMCLi'] * np.cosh(self.dataset['gen']['fEtaLi'])

        self.dataset['gen']['fPHe3'] = self.dataset['gen']['fPtHe3'] * np.cosh(self.dataset['gen']['fEtaHe3'])
        self.dataset['gen']['fPPr'] = self.dataset['gen']['fPtPr'] * np.cosh(self.dataset['gen']['fEtaPr'])
        self.dataset['gen']['fEHe3'] = np.sqrt(self.dataset['gen']['fPHe3']**2 + particleMass['He3']**2)
        self.dataset['gen']['fEPr'] = np.sqrt(self.dataset['gen']['fPPr']**2 + particleMass['Pr']**2)
        self.dataset['gen']['fAlpha'] = self.dataset['gen']['fPhiHe3'] - self.dataset['gen']['fPhiPr']            # separation angle

        self.dataset['gen']['fInnerPTPCHe3'] = self.dataset['gen']['fInnerParamTPCHe3'] * 2
        self.dataset['gen']['fInnerPTPCPr'] = self.dataset['gen']['fInnerParamTPCPr']
        
        self.dataset['gen']['fSignedPTPCHe3'] = self.dataset['gen']['fInnerPTPCHe3'] * self.dataset['gen']['fSignHe3']
        self.dataset['gen']['fSignedPTPCPr'] = self.dataset['gen']['fInnerPTPCPr'] * self.dataset['gen']['fSignPr']

        # beta*gamma (for Bethe-Bloch formula)
        self.dataset['gen']['fBetaGammaHe3'] = self.dataset['gen']['fInnerParamTPCHe3'] / particleMass['He3'] * 2.
        self.dataset['gen']['fBetaGammaPr'] = self.dataset['gen']['fInnerParamTPCPr'] / particleMass['Pr']

        # invariant mass 
        self.dataset['gen']['fPtLi'] = np.sqrt( self.dataset['gen']['fPtHe3']**2 + self.dataset['gen']['fPtPr']**2 + 2*self.dataset['gen']['fPtHe3']*self.dataset['gen']['fPtPr']*np.cos(self.dataset['gen']['fAlpha']) )
        self.dataset['gen']['fMassLi'] = np.sqrt(particleMass['He3']**2 + particleMass['Pr']**2 + 
                                  2 * (self.dataset['gen']['fEHe3'] * self.dataset['gen']['fEPr'] - self.dataset['gen']['fPHe3'] * self.dataset['gen']['fPPr'] * np.cos(self.dataset['gen']['fAlpha'])) )
        
        self.dataset['gen']['fPtResHe3'] = (self.dataset['gen']['fPtMCHe3'] - self.dataset['gen']['fPtHe3']) / self.dataset['gen']['fPtMCHe3']
        self.dataset['gen']['fPtResPr'] = (self.dataset['gen']['fPtMCPr'] - self.dataset['gen']['fPtPr']) / self.dataset['gen']['fPtMCPr']
        self.dataset['gen']['fPtResLi'] = (self.dataset['gen']['fPtMCLi'] - self.dataset['gen']['fPtLi']) / self.dataset['gen']['fPtMCLi']

        self.dataset['gen']['fPResHe3'] = (self.dataset['gen']['fPMCHe3'] - self.dataset['gen']['fInnerPTPCHe3']) / self.dataset['gen']['fPMCHe3']
        self.dataset['gen']['fPResPr'] = (self.dataset['gen']['fPMCPr'] - self.dataset['gen']['fInnerPTPCPr']) / self.dataset['gen']['fPMCPr']
        #self.dataset['gen']['fPResLi'] = (self.dataset['gen']['fPMCLi'] - self.dataset['gen']['fInnterPTPCLi']) / self.dataset['gen']['fPMCLi']

        # cluster size on different layers
        for layer in range(7):
            read4Bits = lambda x, layer: (x >> layer*4) & 0b1111
            self.dataset['gen'][f'fClSizeITS{layer}He3'] = self.dataset['gen']['fItsClusterSizeHe3'].apply(read4Bits, args=(layer,))
            self.dataset['gen'][f'fClSizeITS{layer}Pr'] = self.dataset['gen']['fItsClusterSizePr'].apply(read4Bits, args=(layer,))
        
        self.dataset['gen']['fCosLambdaHe3'] = 1/(np.cosh(self.dataset['gen']['fEtaHe3'])**2)
        self.dataset['gen']['fCosLambdaPr'] = 1/(np.cosh(self.dataset['gen']['fEtaPr'])**2)
        
        self.dataset['gen'].eval('fClSizeITSMeanHe3 = (fClSizeITS0He3 + fClSizeITS1He3 + fClSizeITS2He3 + fClSizeITS3He3+ fClSizeITS4He3 + fClSizeITS5He3 + fClSizeITS6He3) / 7', inplace=True)
        self.dataset['gen'].eval('fClSizeITSMeanPr = (fClSizeITS0Pr + fClSizeITS1Pr + fClSizeITS2Pr + fClSizeITS3Pr+ fClSizeITS4Pr + fClSizeITS5Pr + fClSizeITS6Pr) / 7', inplace=True)
        self.dataset['gen'].eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 * fCosLambdaHe3', inplace=True)
        self.dataset['gen'].eval('fClSizeITSCosLamPr = fClSizeITSMeanPr * fCosLambdaPr', inplace=True)

        self.__selectRecoTracks()

    def visualize(self, config) -> None:

        print('Visualizing...')
        with open(config, 'r') as file:     config = yaml.safe_load(file)

        print(f'Creating output file {config["outputFilePath"]}...')
        outFile = TFile(config['outputFilePath'], 'recreate')
        outDirs = {}
        for dir in config['outDirs']:   outDirs[dir] = outFile.mkdir(dir)

        for key, cfg in config.items():
            
            if key == 'outputFilePath': continue
            if key == 'outDirs':        continue
            
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
                    hDenominator = self.histHandler['gen'].buildTH1(cfg['denVariable']+part, axisSpecDen)

                    eff = self.histHandler['reco'].buildEfficiency(hNumerator, hDenominator)
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()  
                    eff.Write()
        
        outFile.Close()
