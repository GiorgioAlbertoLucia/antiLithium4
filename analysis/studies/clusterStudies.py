'''
    Class to visualize different variables with particular selections on the dataset
'''

from ROOT import TH2F

import sys
sys.path.append('..')
from src.preprocessing import Preprocessor
from src.bethe_bloch_parametrisation import BetheBlochParametrisation
from .studies import StandaloneStudy

sys.path.append('../..')
from framework.src.axis_spec import AxisSpec
from framework.src.hist_info import HistLoadInfo
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

class ClusterSizeParamStudy(StandaloneStudy):

    def __init__(self, config, bb_config, h2_cluster_info:HistLoadInfo):

        super().__init__(config)
        self.dir = ClusterSizeParamStudy.outFile_shared.mkdir('ClusterSizeParametrisation')

        # Bethe-Bloch parameters
        
        # Parameters from data
        self.BetheBlochParams = {'kp1': -0.031712,
                                 'kp2': -45.0275,
                                 'kp3': -0.997645,
                                 'kp4': 1.68228,
                                 'kp5': 0.0108484
                                 #'resolution':  0.09
                                }

        self.bb_param = BetheBlochParametrisation()
        self.bb_param.load_config(bb_config)
        self.bb_param._set_output_dir(self.dir)
        self.bb_param.select_fit_particle('Pr')
        self.bb_param.reset_fit_results()
        self.bb_param.init_config('betagamma')

        h2 = HistHandler.loadHist(h2_cluster_info)
        self.bb_param.upload_h2(h2)
        self.bb_param.upload_bethe_bloch_params(self.BetheBlochParams)

    def rebinx(self, rebin_factor:int=2) -> None:
        self.bb_param.h2.RebinX(rebin_factor)

    def fitBetheBloch(self, **kwargs) -> None:

        self.bb_param.generate_bethe_bloch_points()
        out_fit_results = kwargs.get('out_fit_results', None)
        if out_fit_results:
            self.bb_param.save_fit_results(out_fit_results)
        self.bb_param.fit_bethe_bloch()
    
    def drawBetheBloch(self, canvas_name) -> None:

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Drawing Bethe Bloch with current parameters')
        self.bb_param.draw_bethe_bloch_fit('Pr', canvas_name=canvas_name)
    
    def print_results(self) -> None:

        params = self.bb_param.BetheBloch_params
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'Fit results:')
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'kp1: {params["kp1"]}')
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'kp2: {params["kp2"]}')
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'kp3: {params["kp3"]}')
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'kp4: {params["kp4"]}')
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'kp5: {params["kp5"]}')
