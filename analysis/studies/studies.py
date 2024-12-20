'''
    Class to visualize different variables with particular selections on the dataset
'''
import yaml
from ROOT import TFile

import sys
sys.path.append('..')
from src.preprocessing import Preprocessor

sys.path.append('../..')
from framework.utils.terminal_colors import TerminalColors as tc

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
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Creating output file '+tc.UNDERLINE+tc.CYAN+f'{outFilePath}'+tc.RESET)
        Study.isFileOpen_shared = True

    @classmethod
    def close(cls) -> None:

        if Study.isFileOpen_shared:   Study.outFile_shared.Close()

class StandaloneStudy:

    def __init__(self, config, outFile: TFile):

        with open(config, 'r') as file:     self.config = yaml.safe_load(file)
        self.outFile = outFile