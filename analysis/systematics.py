'''
    Analysis routine to study the effects of systematics in the correlation function using the cut variation technique.
'''

import yaml
import argparse

from torchic import Dataset
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

import sys
sys.path.append('/home/galucia/antiLithium4/analysis')
from src.preprocessing import DataPreprocessor

def systematics(args):

    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.BOLD+tc.WHITE+'\t\tRunning analysis'+tc.RESET)
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET)

    with open(args.cfgInputFile, 'r') as file:     
        cfgInput = yaml.safe_load(file)
    inFilePath = cfgInput['inFilePath']
    outFilePath = cfgInput['outFilePath']
    outQaFilePath = cfgInput['outQaFilePath']
    cfgVisualFile = cfgInput['visualFilePath']
    cfgQaFile = cfgInput['qaFilePath']

    dataset_SE = 
    preprocessor 
