'''
    Class to visualize different variables with particular selections on the dataset
'''
import numpy as np

from ROOT import TF1, TCanvas, TObjArray

import sys
sys.path.append('..')
from utils.particle import PID
from ..src.preprocessing import Preprocessor
from .studies import Study

sys.path.append('../..')
from framework.src.axis_spec import AxisSpec
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

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
