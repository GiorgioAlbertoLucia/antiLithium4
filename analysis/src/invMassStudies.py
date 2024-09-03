'''
    Classes for invariant mass studies
'''
import os
import numpy as np
import pandas as pd
import yaml

from ROOT import TFile, TH1F, TH2F, TF1, TCanvas, gInterpreter, TObjArray

from .studies import Study

import sys
sys.path.append('..')
from utils.particle import PID
from .preprocessing import Preprocessor

sys.path.append('../..')
from framework.src.axisSpec import AxisSpec
from framework.src.histHandler import HistHandler
from framework.utils.terminalColors import TerminalColors as tc

class InvariantMassStudy(Study):

    def __init__(self, preprocessor: Preprocessor, config):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(preprocessor, config)
        self.dir = InvariantMassStudy.outFile_shared.mkdir('invariantMass')

        cfg = self.config['InvMass']
        self.axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name'], cfg['title'])
        self.hBkg = TH1F('InvMassLiBkg', 'Invariant mass; m_{inv} (GeV/c^{2}); Counts', self.axisSpecX.nbins, self.axisSpecX.xmin, self.axisSpecX.xmax)
        self.hSignal = TH1F('InvMassLiSignal', 'Invariant mass; m_{inv} (GeV/c^{2}); Counts', self.axisSpecX.nbins, self.axisSpecX.xmin, self.axisSpecX.xmax)

    def generalSelections(self) -> None:

        selections = [
            '-2 < fNSigmaTPCHe3 < 2',
            '-2 < fNSigmaTPCPr < 2',
        ]

        for selection in selections: self.dataset['full'].query(selection, inplace=True)

    def cutList(self) -> None:
        '''
            Create a dictionary of cuts to apply on the dataset.

            Returns:
                cuts: dictionary with the cuts to be applied on the dataset
                  (keys: variable name; values: [min, max, step])

        '''
        cuts = {
            'fClSizeITSCosLamHe3': [3.5, 4.5, 0.1, 'g'],    # minimum cl size for He3
            'fClSizeITSCosLamPr': [2.8, 3.2, 0.1, 'l'],      # maximum cl size for Pr
            #'fDCAxyHe3': [0.05, 0.015, 0.001, 'labs'],      # maximum abs(DCAxy) for He3
            #'fDCAxyPr': [0.05, 0.015, 0.001, 'labs'],       # maximum abs(DCAxy) for Pr
            'fDCAzHe3': [0.75, 1.25, 0.05, 'labs'],         # maximum abs(DCAz) for He3
            'fDCAzPr': [0.75, 1.25, 0.05, 'labs'],          # maximum abs(DCAz) for Pr
        }

        return cuts
    
    def applyCut(self, cutVariable:str, cutValue:list, cutMethod:str):
        '''
            Apply a cut on the dataset.

            Args:
                cutVariable: variable name to apply the cut
                cutValue: list with the cut values [min, max, step]
                cutMethod: method to apply the cut ('labs' for lesser than absolute,
                                                    'gabs' for greater than absolute,
                                                    'l' for lesser than,
                                                    'g' for greater than)

            Returns:
                cutData: dataset with the cut applied
        '''
        data = self.dataset['full']
        
        if cutMethod == 'labs': cutData = data.query(f'abs({cutVariable}) < {cutValue}', inplace=False)
        elif cutMethod == 'gabs': cutData = data.query(f'abs({cutVariable}) > {cutValue}', inplace=False)
        elif cutMethod == 'l': cutData = data.query(f'{cutVariable} < {cutValue}', inplace=False)
        elif cutMethod == 'g': cutData = data.query(f'{cutVariable} > {cutValue}', inplace=False)
        else: raise ValueError('Invalid cut method. Use "labs", "gabs", "l" or "g".')

        return cutData
        
    def invariantMass(self) -> None:
        '''
            Create a 2D histogram with the invariant mass distribution for different cuts.
        '''
        
        cuts = self.cutList()
        operation = {'labs': '< abs(', 'gabs': '> abs(', 'l': '<', 'g': '>'}

        for cutVar, cutRange in cuts.items():

            cutHists = {}

            for cutValue in np.arange(cutRange[0], cutRange[1], cutRange[2]):
                
                cutData = self.applyCut(cutVar, cutValue, cutRange[3])
                
                title = f'{cutVar} {operation} {cutValue:.2f};m_{{inv}} (GeV/c^{{2}});Counts'
                axisSpecX = AxisSpec(100, 3.7, 3.9, f'invariantMass_{cutVar}_{cutValue}', title)
                
                histHandler = HistHandler.createInstance(cutData)
                hist = histHandler.buildTH1('fMassInvLi', axisSpecX)
                cutHists[cutValue] = hist
            
            hist = TH2F(f'minv_{cutVar}', f'Invariant mass ({cutVar});m_{{inv}} (GeV/c^{{2}});{cutVar}', 100, 3.7, 3.9, len(cuts), cutRange[0]-0.5*cutRange[2], cutRange[1]+0.5*cutRange[2])
            for cutValue, cutHist in cutHists.items(): 
                for bin in range(1, hist.GetNbinsX()+1): hist.SetBinContent(bin, hist.GetYaxis().FindBin(cutValue), cutHist.GetBinContent(bin))

            self.dir.cd()
            hist.Write()    

    def normalizeEventMixingBkg(self, signalPath:str, signalName:str, lowInvMass:float=3.78, upperInvMass:float=3.85)  -> None:
        '''
            Normalize the event mixing background to the signal.

            Args:
                signalPath: path to the signal file
                signalName: name of the signal histogram
                lowInvMass: lower limit of the invariant mass
                upperInvMass: upper limit of the invariant mass
        '''
        signalFile = TFile(signalPath, 'read')
        hSignal = signalFile.Get(signalName)
        signalIntegral = hSignal.Integral(hSignal.FindBin(lowInvMass), hSignal.FindBin(upperInvMass))

        for x in self.dataset['full']['fMassInvLi']: self.hBkg.Fill(x)
        bkgIntegral = self.hBkg.Integral(self.hBkg.FindBin(lowInvMass), self.hBkg.FindBin(upperInvMass))

        self.hBkg.Scale(signalIntegral/bkgIntegral)

        self.dir.cd()
        self.hBkg.Write('InvMassLiNormalized')

    def bkgSubtraction(self, signalPath:str, signalName:str) -> None:
        '''
            Subtract the background from the signal.

            Args:
                signalPath: path to the signal file
                signalName: name of the signal histogram
        '''
        signalFile = TFile(signalPath, 'read')
        hSignal = signalFile.Get(signalName)
        bkgHist = self.dir.Get('InvMassLiNormalized')

        hSubtracted = hSignal.Clone()
        hSubtracted.Add(bkgHist, -1)

        self.dir.cd()
        hSubtracted.Write('InvMassLiSubtracted')
        