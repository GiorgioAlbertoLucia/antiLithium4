'''
    Class to visualize different variables with particular selections on the dataset
'''
from ROOT import RooRealVar

import sys
sys.path.append('..')
from .studies import StandaloneStudy

from torchic import Roofitter, HistLoadInfo
from torchic.core.histogram import load_hist
from torchic.physics.calibration import bethe_bloch_calibration

class BetheBlochStudy(StandaloneStudy):
    
    def __init__(self, config, outputFile, h2_dedx_info:HistLoadInfo):

        super().__init__(config, outputFile)
        self.dir = self.outFile.mkdir('BetheBloch')

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
        self.dEdx = load_hist(h2_dedx_info)

    def rebinx(self, rebin_factor:int=2) -> None:
        self.dEdx.RebinX(rebin_factor)

    def fit_BetheBloch(self, **kwargs) -> None:

        xMin = kwargs.get('xMin', 0.55)
        xMax = kwargs.get('xMax', 4.0)
        yMin = kwargs.get('yMin', 0.0)
        yMax = kwargs.get('yMax', 3000.0)

        bb_options = {'first_bin_fit_by_slices': self.dEdx.GetXaxis().FindBin(xMin),
                      'last_bin_fit_by_slices': self.dEdx.GetXaxis().FindBin(xMax),
                      'output_dir': self.dir,
                      'signal_func_name': 'gaus_0',
                      'signal_range': (yMin, yMax),
                     }
        dEdx = RooRealVar('x', 'x', yMin, yMax)
        fitter = Roofitter(dEdx, ['gaus'])
        fitter.init_param('gaus_0_mean', 500., 200., 800.)
        fitter.init_param('gaus_0_sigma', 50, 0., 500.)
    
        self.BetheBlochParams = bethe_bloch_calibration(self.dEdx, self.dir, fitter, **bb_options)