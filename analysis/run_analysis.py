import yaml
import argparse

from torchic import Dataset
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

import sys
sys.path.append('..')
from analysis.src.preprocessing import DataPreprocessor

@timeit
def preprocessing(args) -> DataPreprocessor:

    with open(args.cfgInputFile, 'r') as file:     cfgInput = yaml.safe_load(file)
    inFilePath = cfgInput['inFilePath']
    outFilePath = cfgInput['outFilePath']
    cfgVisualFile = cfgInput['visualFilePath']

    dataset = Dataset(inFilePath, tree_name='O2lithium4table', folder_name='DF*')
    preprocessor = DataPreprocessor(dataset)
    preprocessor.define_variables(args.isBkgUS)
    #preprocessor.visualize(outFilePath, cfgVisualFile)

    preprocessor.selections_He3()
    preprocessor.define_kstar()
    preprocessor.visualize(outFilePath.replace('.root', '_selectionsHe3.root'), cfgVisualFile)
    preprocessor.selections_Pr()
    preprocessor.visualize(outFilePath.replace('.root', '_selectionsPr.root'), cfgVisualFile)

    return preprocessor

if __name__ == '__main__':

    print()
    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--config-file', dest='cfgInputFile',
                        help='path to the YAML file with configuration.', default='')
    parser.add_argument('--isBkgUS', dest='isBkgUS', help='Flag for US background.', action='store_true')
    args = parser.parse_args()

    if args.cfgInputFile == '':
        print(tc.RED+'[ERROR]: '+tc.RESET+'No config file provided, exiting.')
        exit(1)

    preprocessor = preprocessing(args)
