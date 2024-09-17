'''
    Class to visualize different variables with particular selections on the dataset
'''

from ROOT import TH2F

import sys
sys.path.append('..')
from src.preprocessing import Preprocessor
from src.bethe_bloch_parametrisation import BetheBlochParametrisation
from .studies import Study

sys.path.append('../..')
from framework.src.axis_spec import AxisSpec
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

class ClusterSizeParamStudy(Study):

    def __init__(self, preprocessor: Preprocessor, config, bb_config):

        super().__init__(preprocessor, config)
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

        h2 = TH2F('test', 'Bethe Bloch Fit; #beta#gamma; #LT Cluster size #GT #times #LT cos#lambda #GT; Counts', 200, 0.3, 5.0, 90, 0, 15)
        for x, y in zip(self.dataset['reco']['fBetaGammaPr'], self.dataset['reco']['fClSizeITSCosLamPr']):
            h2.Fill(x, y)
        self.bb_param.upload_h2(h2)
        self.bb_param.upload_bethe_bloch_params(self.BetheBlochParams)

    def fitBetheBloch(self, **kwargs) -> None:

        self.bb_param.generate_bethe_bloch_points()
        out_fit_results = kwargs.get('out_fit_results', None)
        if out_fit_results:
            self.bb_param.save_fit_results(out_fit_results)
        self.bb_param.fit_bethe_bloch()
    
    def drawBetheBloch(self, canvas_name) -> None:

        print(tc.GREEN+'[INFO]: '+tc.RESET+'Drawing Bethe Bloch with current parameters')
        self.bb_param.draw_bethe_bloch_fit('Pr', canvas_name=canvas_name)
