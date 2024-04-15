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
from utils.particle import particleMass
from utils.plotUtils import setPIDlabels

class Preprocessor:

    def __init__(self) -> None:
        '''
            - dataset (pd.DataFrame):  data to be preprocessed
            - recoDataset (pd.DataFrame): reconstructed tracks
            - genDataset (pd.DataFrame): generated particles

        '''
        
        self.genDataset = pd.DataFrame()
        self.recoDataset = pd.DataFrame()




    def open(self, inputFile, **kwargs):
        
        # check extension
        if inputFile.endswith('.root'):
            
            treeName = None
            dirPrefix = None

            for key, value in kwargs.items():
                if key == 'treeName':   treeName = value
                if key == 'dirPrefix':  dirPrefix = value
            
            print(f'Opening {inputFile}...')
            print(f'Using tree {treeName} and directory prefix {dirPrefix}')
            th = TreeHandler(inputFile, treeName, folder_name=dirPrefix)
            self.genDataset = th.get_data_frame()

        else:   print('Error: file extension not supported')

    def __selectRecoTracks(self) -> None:
            
        self.recoDataset = self.genDataset.query('fPIDtrkHe3 != 0xFFFFF', inplace=False)

    # private methods
    def __buildHist(self, variable, name, title, nBins, xMin, xMax, opt:str='reco') -> TH1F:
        
        hist = TH1F(name, title, nBins, xMin, xMax)
        if opt == 'gen':
            for value in self.genDataset[variable]:    hist.Fill(value)
        elif opt == 'reco':
            for value in self.recoDataset[variable]:    hist.Fill(value)
        else:
            raise ValueError(f'opt can only be gen or reco. Empty histogram {name} is returned.')
        return hist
    
    def __buildHist2(self, xVariable, yVariable, name, title, nXBins, xMin, xMax, nYBins, yMin, yMax, opt:str='reco') -> TH2F:

        hist2 = TH2F(name, title, nXBins, xMin, xMax, nYBins, yMin, yMax)
        if opt == 'gen':
            for x, y in zip(self.genDataset[xVariable], self.genDataset[yVariable]):    hist2.Fill(x, y)
        elif opt == 'reco':
            for x, y in zip(self.recoDataset[xVariable], self.recoDataset[yVariable]):    hist2.Fill(x, y)
        else:
            raise ValueError(f'opt can only be gen or reco. Empty histogram {name} is returned.')
        return hist2
    
    def __buildEfficiency(self, num, den, title, nBins, xMin, xMax) -> TEfficiency:

        hist1 = TH1F(num, title, nBins, xMin, xMax)
        hist2 = TH1F(den, title, nBins, xMin, xMax)

        for value in self.recoDataset[num]:    hist1.Fill(value)
        for value in self.genDataset[den]:    hist2.Fill(value)

        eff = TH1F(num+'_eff', title, nBins, xMin, xMax)
        for i in range(1, nBins+1):
            if hist2.GetBinContent(i) != 0:     eff.SetBinContent(i, hist1.GetBinContent(i) / hist2.GetBinContent(i))
            else:                               eff.SetBinContent(i, 0)

        return eff


    # public methods
    def defineVariables(self) -> None:

        self.genDataset['fPtMCLi'] = self.genDataset['fSignedPtMC'] 
        self.genDataset['fMassMCLi'] = self.genDataset['fMassMC']

        # pt MC stores sign
        self.genDataset.eval('fSignHe3 = (-1)**(fPtMCHe3 < 0)', inplace=True)
        self.genDataset.eval('fSignPr = (-1)**(fPtMCPr < 0)', inplace=True)
        self.genDataset.eval('fSignLi = (-1)**(fPtMCLi < 0)', inplace=True)

        self.genDataset['fSignedPtHe3'] = self.genDataset['fPtHe3'] * self.genDataset['fSignHe3']
        self.genDataset['fSignedPtPr'] = self.genDataset['fPtPr'] * self.genDataset['fSignPr']

        self.genDataset['fPtMCHe3'] = abs(self.genDataset['fPtMCHe3'])
        self.genDataset['fPtMCPr'] = abs(self.genDataset['fPtMCPr'])
        self.genDataset['fPtMCLi'] = abs(self.genDataset['fPtMCLi'])

        self.genDataset['fPMCHe3'] = self.genDataset['fPtMCHe3'] * np.cosh(self.genDataset['fEtaHe3'])
        self.genDataset['fPMCPr'] = self.genDataset['fPtMCPr'] * np.cosh(self.genDataset['fEtaPr'])
        #self.genDataset['fPMCLi'] = self.genDataset['fPtMCLi'] * np.cosh(self.genDataset['fEtaLi'])

        self.genDataset['fPHe3'] = self.genDataset['fPtHe3'] * np.cosh(self.genDataset['fEtaHe3'])
        self.genDataset['fPPr'] = self.genDataset['fPtPr'] * np.cosh(self.genDataset['fEtaPr'])
        self.genDataset['fEHe3'] = np.sqrt(self.genDataset['fPHe3']**2 + particleMass['He3']**2)
        self.genDataset['fEPr'] = np.sqrt(self.genDataset['fPPr']**2 + particleMass['Pr']**2)
        self.genDataset['fAlpha'] = self.genDataset['fPhiHe3'] - self.genDataset['fPhiPr']            # separation angle

        self.genDataset['fInnerPTPCHe3'] = self.genDataset['fInnerParamTPCHe3'] * 2
        self.genDataset['fInnerPTPCPr'] = self.genDataset['fInnerParamTPCPr']
        
        self.genDataset['fSignedPTPCHe3'] = self.genDataset['fInnerPTPCHe3'] * self.genDataset['fSignHe3']
        self.genDataset['fSignedPTPCPr'] = self.genDataset['fInnerPTPCPr'] * self.genDataset['fSignPr']

        # beta*gamma (for Bethe-Bloch formula)
        self.genDataset['fBetaGammaHe3'] = self.genDataset['fInnerParamTPCHe3'] / particleMass['He3'] * 2.
        self.genDataset['fBetaGammaPr'] = self.genDataset['fInnerParamTPCPr'] / particleMass['Pr']

        # invariant mass 
        self.genDataset['fPtLi'] = np.sqrt( self.genDataset['fPtHe3']**2 + self.genDataset['fPtPr']**2 + 2*self.genDataset['fPtHe3']*self.genDataset['fPtPr']*np.cos(self.genDataset['fAlpha']) )
        self.genDataset['fMassLi'] = np.sqrt(particleMass['He3']**2 + particleMass['Pr']**2 + 
                                  2 * (self.genDataset['fEHe3'] * self.genDataset['fEPr'] - self.genDataset['fPHe3'] * self.genDataset['fPPr'] * np.cos(self.genDataset['fAlpha'])) )
        
        self.genDataset['fPtResHe3'] = (self.genDataset['fPtMCHe3'] - self.genDataset['fPtHe3']) / self.genDataset['fPtMCHe3']
        self.genDataset['fPtResPr'] = (self.genDataset['fPtMCPr'] - self.genDataset['fPtPr']) / self.genDataset['fPtMCPr']
        self.genDataset['fPtResLi'] = (self.genDataset['fPtMCLi'] - self.genDataset['fPtLi']) / self.genDataset['fPtMCLi']

        self.genDataset['fPResHe3'] = (self.genDataset['fPMCHe3'] - self.genDataset['fInnerPTPCHe3']) / self.genDataset['fPMCHe3']
        self.genDataset['fPResPr'] = (self.genDataset['fPMCPr'] - self.genDataset['fInnerPTPCPr']) / self.genDataset['fPMCPr']
        #self.genDataset['fPResLi'] = (self.genDataset['fPMCLi'] - self.genDataset['fInnterPTPCLi']) / self.genDataset['fPMCLi']

        # cluster size on different layers
        for layer in range(7):
            read4Bits = lambda x, layer: (x >> layer*4) & 0b1111
            self.genDataset[f'fClSizeITS{layer}He3'] = self.genDataset['fItsClusterSizeHe3'].apply(read4Bits, args=(layer,))
            self.genDataset[f'fClSizeITS{layer}Pr'] = self.genDataset['fItsClusterSizePr'].apply(read4Bits, args=(layer,))
        
        self.genDataset['fCosLambdaHe3'] = np.sqrt( 1 - 1/(np.cosh(self.genDataset['fEtaHe3'])**2) )
        self.genDataset['fCosLambdaPr'] = np.sqrt( 1 - 1/(np.cosh(self.genDataset['fEtaPr'])**2) )
        
        self.genDataset.eval('fClSizeITSMeanHe3 = (fClSizeITS0He3 + fClSizeITS1He3 + fClSizeITS2He3 + fClSizeITS3He3+ fClSizeITS4He3 + fClSizeITS5He3 + fClSizeITS6He3) / 7', inplace=True)
        self.genDataset.eval('fClSizeITSMeanPr = (fClSizeITS0Pr + fClSizeITS1Pr + fClSizeITS2Pr + fClSizeITS3Pr+ fClSizeITS4Pr + fClSizeITS5Pr + fClSizeITS6Pr) / 7', inplace=True)
        self.genDataset.eval('fClSizeITSCosLamHe3 = fClSizeITSMeanHe3 * fCosLambdaHe3', inplace=True)
        self.genDataset.eval('fClSizeITSCosLamPr = fClSizeITSMeanPr * fCosLambdaPr', inplace=True)

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
            
            if cfg['particle'] is None:
                if 'TH1' in cfg['type']:
                    hist = self.__buildHist(cfg['xVariable'], cfg['name'], cfg['title'], cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['opt'])
                    if cfg['xVariable'] == 'fPIDtrk':   setPIDlabels(hist, 'x')
                    if cfg['dir'] != 'None':  outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist.Write()
                if 'TH2' in cfg['type']:
                    hist2 = self.__buildHist2(cfg['xVariable'], cfg['yVariable'], cfg['name'], cfg['title'], cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['opt'])
                    if cfg['xVariable'] == 'fPIDtrk':   setPIDlabels(hist2, 'x')
                    if cfg['yVariable'] == 'fPIDtrk':   setPIDlabels(hist2, 'y')
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    hist2.Write()
                if 'TEfficiency' in cfg['type']:
                    eff = self.__buildEfficiency(cfg['numVariable'], cfg['denVariable'], cfg['title'], cfg['nXBins'], cfg['xMin'], cfg['xMax'])
                    if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                    else:                       outFile.cd()
                    eff.Write()
                continue
            else:
                for part in cfg['particle']:
                    if 'TH1' in cfg['type']:
                        hist = self.__buildHist(cfg['xVariable']+part, cfg['name']+part, cfg['title']+f' {part}', cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['opt'])
                        if cfg['xVariable'] == 'fPIDtrk':   setPIDlabels(hist, 'x')
                        if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                        else:                       outFile.cd()
                        hist.Write()
                    if 'TH2' in cfg['type']:
                        hist2 = self.__buildHist2(cfg['xVariable']+part, cfg['yVariable']+part, cfg['name']+part, cfg['title']+f' {part}', cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['opt'])
                        if cfg['xVariable'] == 'fPIDtrk':   setPIDlabels(hist2, 'x')
                        if cfg['yVariable'] == 'fPIDtrk':   setPIDlabels(hist2, 'y')
                        if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                        else:                       outFile.cd()
                        hist2.Write()
                    if 'TEfficiency' in cfg['type']:
                        eff = self.__buildEfficiency(cfg['numVariable']+part, cfg['denVariable']+part, cfg['title']+f' {part}', cfg['nXBins'], cfg['xMin'], cfg['xMax'])
                        if cfg['dir'] != 'None':    outDirs[cfg['dir']].cd()
                        else:                       outFile.cd()  
                        eff.Write()
        
        outFile.Close()
