'''
    Class to visualize different variables with particular selections on the dataset
'''
import os
from ROOT import TF1, TCanvas, gInterpreter, TObjArray

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BETHEBLOCH_DIR = os.path.join(CURRENT_DIR, '..', 'include', 'BetheBloch.hh')
gInterpreter.ProcessLine(f'#include "{BETHEBLOCH_DIR}"')
from ROOT import BetheBloch

import sys
sys.path.append('..')
from .studies import StandaloneStudy

sys.path.append('../..')
from framework.src.axis_spec import AxisSpec
from framework.src.hist_info import HistLoadInfo
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

class BetheBlochStudy(StandaloneStudy):
    
    def __init__(self, config, h2_dedx_info:HistLoadInfo):

        super().__init__(config)
        self.dir = BetheBlochStudy.outFile_shared.mkdir('BetheBloch')

        # Bethe-Bloch parameters
    
        # Parameters from data
        self.BetheBlochParams = {'kp1':  -321.34,
                                 'kp2':  0.6539,
                                 'kp3':  1.591,
                                 'kp4':  0.8225,
                                 'kp5':  2.363,
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
        '''
        self.dEdx = HistHandler.loadHist(h2_dedx_info)

    def rebinx(self, rebin_factor:int=2) -> None:
        self.dEdx.RebinX(rebin_factor)

    def fit_BetheBloch(self, **kwargs) -> None:

        xMin = kwargs.get('xMin', 0.0)
        xMax = kwargs.get('xMax', 5.0)
        
        self.BBcurve = TF1(f'BetheBlochHe3', BetheBloch, xMin, xMax, len(self.BetheBlochParams.values()))

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
        self.BBcurve.SetRange(kwargs.get('xMinFit', 0.25), kwargs.get('xMaxFit', 5.))
        self.BBhist.Fit(self.BBcurve, 'RM+')

        self.dir.cd()
    
        print(tc.BOLD+'Bethe Bloch parameters:'+tc.RESET)
        for i, (parName, param) in enumerate(self.BetheBlochParams.items()):    
            self.BetheBlochParams[parName] = self.BBcurve.GetParameter(i)
            print(tc.RED+f'     {parName}:\t'+tc.RESET+f'{self.BBcurve.GetParameter(i)}')
        print(tc.BOLD+'     chi2/ndf:\t'+tc.RESET+str(self.BBcurve.GetChisquare())+'/'+str(self.BBcurve.GetNDF()))
        print()
        
    
    def draw(self) -> None:

        canvas = TCanvas(f'BBHe3', f'Bethe Bloch curve - He3')
        #hframe = canvas.DrawFrame(cfg['xMin'], cfg['yMin'], cfg['xMax'], cfg['yMax'], f'Bethe Bloch He3; #beta #gamma; #frac{{dE}}{{dX}} (a.u.)')

        self.dir.cd()    
        self.BBhist.SetName('BBhistHe3')
        self.BBhist.Write()
        self.BBcurve.Write()
        self.dEdx.Write()
            
        canvas.cd()
        self.dEdx.Draw('colz same')
        self.BBcurve.Draw('same')
        self.dir.cd()
        canvas.Write()
