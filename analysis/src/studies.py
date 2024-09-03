'''
    Class to visualize different variables with particular selections on the dataset
'''
import os
import numpy as np
import pandas as pd
import yaml

from ROOT import TFile, TH1F, TH2F, TF1, TCanvas, gInterpreter, TObjArray

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BETHEBLOCH_DIR = os.path.join(CURRENT_DIR, '..', 'include', 'BetheBloch.hh')
gInterpreter.ProcessLine(f'#include "{BETHEBLOCH_DIR}"')
from ROOT import BetheBloch


import sys
sys.path.append('..')
from utils.particle import PID
from .preprocessing import Preprocessor

sys.path.append('../..')
from framework.src.axisSpec import AxisSpec
from framework.src.histHandler import HistHandler
from framework.utils.terminalColors import TerminalColors as tc

class Study:

    # static variable to check if a shared file is open between classes
    isFileOpen_shared = False
    outFile_shared = None

    def __init__(self, preprocessor: Preprocessor, config):

        self.dataset = preprocessor.dataset
        with open(config, 'r') as file:     self.config = yaml.safe_load(file)

        self.outFilePath = self.config['studiesOutputFilePath']
        if Study.isFileOpen_shared == False:  self.openFile(self.outFilePath)
    
    @classmethod
    def openFile(cls, outFilePath) -> None:
        
        Study.outFile_shared = TFile(outFilePath, 'recreate')
        print('Creating output file '+tc.UNDERLINE+tc.CYAN+f'{outFilePath}'+tc.RESET+'...')
        Study.isFileOpen_shared = True

    @classmethod
    def close(cls) -> None:

        if Study.isFileOpen_shared:   Study.outFile_shared.Close()


class PIDforTrkStudy(Study):

    def __init__(self, preprocessor: Preprocessor, config):

        super().__init__(preprocessor, config)
        self.dir = PIDforTrkStudy.outFile_shared.mkdir('PIDforTracking')

    def pidVsPtRes(self):

        cfg = self.config['ptRes']
        for part in ['He3', 'Pr']:
            for pidHP in np.delete(self.dataset['reco'][f'fPIDtrk{part}'].unique(), np.where(self.dataset['reco'][f'fPIDtrk{part}'].unique() == 0xFFFFF)):
                
                title = cfg['title'].split(';', 1)[0] + f' {part} - PID trk hp {PID[pidHP]["label"]};' + cfg['title'].split(';', 1)[1]
                data = self.dataset['reco'].query(f'fPIDtrk{part} == {pidHP}', inplace=False)
                
                axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+f'{part}_pid_{PID[pidHP]["label"]}', title)
                axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+f'{part}_pid_{PID[pidHP]["label"]}', title)

                histHandler = HistHandler.createInstance(data)
                hist = histHandler.buildTH2(cfg['xVariable']+f'{part}', cfg['yVariable']+f'{part}', axisSpecX, axisSpecY)
                self.dir.cd()
                hist.Write()

class PtResolutionStudy(Study):

    def __init__(self, preprocessor: Preprocessor, config):

        super().__init__(preprocessor, config)
        self.dir = PtResolutionStudy.outFile_shared.mkdir('ptResolution')

    def ptResolution(self) -> None:

        cfg = self.config['ptRes']

        for part in cfg['particle']:

            title = cfg['title'].split(';', 1)[0] + f' {part} (Reconstructed only);' + cfg['title'].split(';', 1)[1]
            axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+f'{part}', title)
            axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+f'{part}', title)
            
            histHandler = HistHandler.createInstance(self.dataset['reco'])
            hist = histHandler.buildTH2(cfg['xVariable']+f'{part}', cfg['yVariable']+f'{part}', axisSpecX, axisSpecY)
            self.dir.cd()
            hist.Write()

class H3inTrkStudy(Study):

    def __init__(self, preprocessor: Preprocessor, config):

        super().__init__(preprocessor, config)
        self.dir = H3inTrkStudy.outFile_shared.mkdir('H3identified')

        H3pidHp = 6
        self.H3data = self.dataset['reco'].query(f'fPIDtrkHe3 == {H3pidHp}', inplace=False)

        '''
        self.H3fitPtResParams = {'kp0': -1, 
                                 'kp1': 2, 
                                 #'kp2': 1.3
                                 }
        '''
        self.H3fitPtResParams = {'kp0': -1, 
                                 'kp1': 2, 
                                 'kp2': 0.01
                                 }

    def fitH3pt(self) -> None:

        cfg = self.config['ptResNotNorm']

        axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name'], cfg['title'])
        axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name'], cfg['title'])
        
        histHandler = HistHandler.createInstance(self.H3data)
        self.hPtRes = histHandler.buildTH2(cfg['xVariable']+'He3', cfg['yVariable']+'He3', axisSpecX, axisSpecY)

        xMin = self.H3data[cfg['xVariable']+'He3'].min()
        xMax = self.H3data[cfg['xVariable']+'He3'].max()
        # self.H3fitPtRes = TF1('fitH3PtRes', '[0] + exp(-[1]*(x-[2]))', xMin, xMax) # exponential fit
        # self.H3fitPtRes = TF1('fitH3PtRes', '[0] + [1]*x', xMin, xMax) # pol1 fit
        self.H3fitPtRes = TF1('fitH3PtRes', '[0] + [1]*x + [2]*x^2', xMin, xMax) # pol2 fit

        gaus = TF1('gaus', 'gaus')
        results = TObjArray()
        self.hPtRes.FitSlicesY(gaus, 0, -1, 0, 'Q', results)

        self.H3hist = results[1]
        H3res = results[2]

        # error on points = 10% resolution
        for i in range(1, self.H3hist.GetNbinsX()+1):    self.H3hist.SetBinError(i, 0.2*H3res.GetBinContent(i))
        for ipar, param in self.H3fitPtResParams.items():    self.H3fitPtRes.SetParameter(ipar, param)
        self.H3hist.Fit(self.H3fitPtRes, 'RM+')

        print(tc.BOLD+'H3 resolution fit parameters:'+tc.RESET)
        for i, (parName, param) in enumerate(self.H3fitPtResParams.items()):    
            self.H3fitPtResParams[parName] = self.H3fitPtRes.GetParameter(i)
        
            print(tc.RED+f'     {parName}:\t'+tc.RESET+f'{self.H3fitPtRes.GetParameter(i)}')
        print(tc.BOLD+'     chi2/ndf:\t'+tc.RESET+str(self.H3fitPtRes.GetChisquare())+'/'+str(self.H3fitPtRes.GetNDF()))
        print()

    def drawH3pt(self) -> None:

        cfg = self.config['ptResNotNorm']

        canvas = TCanvas(f'H3pidfortrk', f'He3 pt resolution')
        hframe = canvas.DrawFrame(cfg['xMin'], cfg['yMin'], cfg['xMax'], cfg['yMax'], 'He3 pt resolution - H3 in PID for tracking; p_{T} (GeV/#it{c}); p_{T}^{true} - p_{T}^{reco} (a.u.)')

        self.dir.cd()    
        self.H3hist.SetName('H3pidfortrk_hist')
        self.H3hist.Write()
        self.H3fitPtRes.Write()
        self.hPtRes.Write()
            
        canvas.cd()
        self.hPtRes.Draw()
        self.H3fitPtRes.Draw('same')
        self.dir.cd()
        canvas.Write()

class BetheBlochStudy(Study):
    
    def __init__(self, preprocessor: Preprocessor, config):

        super().__init__(preprocessor, config)
        self.dir = BetheBlochStudy.outFile_shared.mkdir('BetheBloch')

        # Bethe-Bloch parameters
        '''
        # Parameters from data
        self.BetheBlochParams = {'kp1':  -170.929,
                            'kp2':  -0.02554,
                            'kp3':  1.58698,
                            'kp4':  0.97078,
                            'kp5':  2.91625,
                            #'resolution':  0.09
                }
        '''
        # Parameters from MC
        self.BetheBlochParams = {'kp1':  -187.9255,
                            'kp2':  -0.26878,
                            'kp3':  1.16252,
                            'kp4':  1.15149,
                            'kp5':  2.47912,
                            #'resolution':  0.09
                }

    def fitBetheBloch(self, **kwargs) -> None:

        cfg = self.config['BetheBloch']
        
        axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name']+f'He3', cfg['title'])
        axisSpecY = AxisSpec(cfg['nYBins'], cfg['yMin'], cfg['yMax'], cfg['name']+f'He3', cfg['title'])
            
        histHandler = HistHandler.createInstance(self.dataset['reco'])
        self.dEdx = histHandler.buildTH2(cfg['xVariable']+'He3', cfg['yVariable']+'He3', axisSpecX, axisSpecY)

        self.BBcurve = TF1(f'BetheBlochHe3', BetheBloch, cfg['xMin'], cfg['xMax'], len(self.BetheBlochParams.values()))

        for i, (parName, param) in enumerate(self.BetheBlochParams.items()):    
            self.BBcurve.SetParameter(i, param)
            self.BBcurve.SetParName(i, parName)
        
        if 'rebin' in kwargs:    self.dEdx.RebinX(kwargs['rebin'])
        
        gaus = TF1('gaus', 'gaus')
        results = TObjArray()
        self.dEdx.FitSlicesY(gaus, 0, -1, 0, 'Q', results)
            
        self.BBhist = results[1]
        BBres = results[2]
        # for i in range(1, self.BBhist.GetNbinsX()+1):    self.BBhist.SetBinError(i, BBres.GetBinContent(i))
        for i in range(1, self.BBhist.GetNbinsX()+1):    self.BBhist.SetBinError(i, 0.09*self.BBhist.GetBinContent(i)) # assume 9% resolution
        self.BBcurve.SetRange(kwargs.get('xMinFit', 0.25), kwargs.get('xMaxFit', 0.78))
        self.BBhist.Fit(self.BBcurve, 'RM+')

        self.dir.cd()
    
        print(tc.BOLD+'Bethe Bloch parameters:'+tc.RESET)
        for i, (parName, param) in enumerate(self.BetheBlochParams.items()):    
            self.BetheBlochParams[parName] = self.BBcurve.GetParameter(i)
            print(tc.RED+f'     {parName}:\t'+tc.RESET+f'{self.BBcurve.GetParameter(i)}')
        print(tc.BOLD+'     chi2/ndf:\t'+tc.RESET+str(self.BBcurve.GetChisquare())+'/'+str(self.BBcurve.GetNDF()))
        print()
        
    
    def drawBetheBloch(self) -> None:

        cfg = self.config['BetheBloch']

        canvas = TCanvas(f'BBHe3', f'Bethe Bloch curve - He3')
        hframe = canvas.DrawFrame(cfg['xMin'], cfg['yMin'], cfg['xMax'], cfg['yMax'], f'Bethe Bloch He3; #beta #gamma; #frac{{dE}}{{dX}} (a.u.)')

        self.dir.cd()    
        self.BBhist.SetName('BBhistHe3')
        self.BBhist.Write()
        self.BBcurve.Write()
        self.dEdx.Write()
            
        canvas.cd()
        self.dEdx.Draw()
        self.BBcurve.Draw('same')
        self.dir.cd()
        canvas.Write()
