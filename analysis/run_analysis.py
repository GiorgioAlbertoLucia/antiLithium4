import yaml
import numpy as np
import argparse

from torchic import Dataset
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

import sys
sys.path.append('..')
from analysis.src.preprocessing import DataPreprocessor

@timeit
def preprocessing(args) -> DataPreprocessor:

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

    tree_names = cfgInput['treeNames']
    folder_name = cfgInput.get('folderName', None)

    dataset = None
    for treeName in tree_names:
        if dataset is None:
            dataset = Dataset.from_root(inFilePath, tree_name=treeName, folder_name=folder_name)
        else:
            dataset = dataset.concat(Dataset.from_root(inFilePath, tree_name=treeName, folder_name=folder_name), axis=1)
    
    print(tc.GREEN+'[INFO]: '+tc.RESET+f'{dataset.columns=}')

    if args.USonly: 
        dataset.query('fIsBkgUS == 1', inplace=True)
    if args.LSonly: 
        dataset.query('fIsBkgUS == 0', inplace=True)
    preprocessor = DataPreprocessor(dataset)
    preprocessor.define_variables()
    #preprocessor.visualize(outFilePath, cfgVisualFile)

    preprocessor.selections_He3()
    #preprocessor.define_kstar()
    #preprocessor.visualize(outFilePath.replace('.root', '_selectionsHe3.root'), cfgVisualFile)
    preprocessor.selections_Pr()
    preprocessor.define_kstar()
    preprocessor.visualize(outFilePath.replace('.root', '_selectionsPr.root'), cfgVisualFile)
    if args.qa: preprocessor.visualize(outQaFilePath.replace('.root', '_selectionsPr.root'), cfgQaFile)

    return preprocessor

@timeit
def close_pair_rejection(args) -> DataPreprocessor:

    with open(args.cfgInputFile, 'r') as file:     cfgInput = yaml.safe_load(file)
    inFilePath = cfgInput['inFilePath']
    outFilePath = cfgInput['outFilePath']
    cfgVisualFile = cfgInput['visualFilePath']

    dataset = Dataset(inFilePath, tree_name='O2lithium4table', folder_name='DF*')
    preprocessor = DataPreprocessor(dataset)
    preprocessor.define_variables(args.isBkgUS)
    preprocessor.close_pair_selection()
    preprocessor.visualize(outFilePath.replace('.root', '_closePair.root'), cfgVisualFile)

    return preprocessor


if __name__ == '__main__':

    print()
    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--config-file', dest='cfgInputFile',
                        help='path to the YAML file with configuration.', default='')
    parser.add_argument('--QA', dest='qa', help='Draw QA plots.', action='store_true')
    parser.add_argument('--isBkgUS', dest='isBkgUS', help='Flag for US background.', action='store_true')
    parser.add_argument('--USonly', dest='USonly', help='Flag to analyse only US pairs.', action='store_true')
    parser.add_argument('--LSonly', dest='LSonly', help='Flag to analyse only LS pairs.', action='store_true')
    parser.add_argument('--closePair', dest='closePair', help='Flag for close pair rejection.', action='store_true')
    args = parser.parse_args()

    if args.cfgInputFile == '':
        print(tc.RED+'[ERROR]: '+tc.RESET+'No config file provided, exiting.')
        exit(1)

    if args.closePair:
        preprocessor = close_pair_rejection(args)
    else:
        preprocessor = preprocessing(args)
