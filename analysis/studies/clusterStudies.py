'''
    Class to visualize different variables with particular selections on the dataset
'''

from ROOT import RooRealVar, TCanvas, TF1

import sys
sys.path.append('..')
from .studies import StandaloneStudy

from torchic import Roofitter, HistLoadInfo
from torchic.core.histogram import load_hist
from torchic.physics.calibration import cluster_size_calibration

class ClusterSizeParamStudy(StandaloneStudy):

    def __init__(self, config, h2_cluster_info:HistLoadInfo):

        super().__init__(config)
        self.dir = self.outFile.mkdir('ClusterSizeParametrisation')

        # Bethe-Bloch parameters
        
        # Parameters from data
        self.clusterSizeParams = {'kp1': 0.903,
                                 'kp2': 2.014,
                                 'kp3': 2.440,
                                 'charge': 1.,
                                 'kp4': 1.,
                                },

        self.h2_clsize = load_hist(h2_cluster_info)

    def rebinx(self, rebin_factor:int=2) -> None:
        self.bb_param.h2.RebinX(rebin_factor)

    def fit(self, **kwargs) -> None:
        
        xMin = kwargs.get('xMin', 0.5)
        xMax = kwargs.get('xMax', 2.7)
        yMin = kwargs.get('yMin', 1.0)
        yMax = kwargs.get('yMax', 9.5)

        cl_options = {'first_bin_fit_by_slices': self.h2_clsize.GetXaxis().FindBin(xMin),
                      'last_bin_fit_by_slices': self.h2_clsize.GetXaxis().FindBin(xMax),
                      'output_dir': self.dir,
                      'fit_range': [yMin, yMax],
                      'simil_bethe_bloch_pars': {'kp1': 0.903,
                                                 'kp2': 2.014,
                                                 'kp3': 2.440,
                                                 'charge': 1.,
                                                 'kp4': 1.,
                                                },
                      'signal_func_name': 'exp_mod_gaus_0',
        }

        x = RooRealVar('x', 'x', 1., 9.5)
        fitter = kwargs.get('fitter', None)
        if fitter is None:
            fitter = Roofitter(x, ['exp_mod_gaus'])
            fitter.init_param('exp_mod_gaus_0_mean', 5., 2., 6.)
            fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0., 10.)
            fitter.init_param('exp_mod_gaus_0_tau', 1., -10., 10.)    
        self.clusterSizeParams, its_resolution_pars = cluster_size_calibration(self.h2_clsize, self.dir, fitter, **cl_options)
    
    def draw(self, canvas_name:str='c_ClSizeAndCurve') -> None:

        simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
        for ipar, par in enumerate(self.clusterSizeParams.values()):
            simil_bethe_bloch_func.SetParameter(ipar, par)
        c_clsize = TCanvas(canvas_name, 'canvas', 800, 600)
        self.h2_clsize.Draw('colz')
        simil_bethe_bloch_func.Draw('same')

        self.dir.cd()
        self.h2_clsize.Write()
        simil_bethe_bloch_func.Write()
        c_clsize.Write()