import yaml
import argparse
import sys
sys.path.append('..')
from framework.src.data_handler import TableHandler, TaskHandler
from framework.utils.timeit import timeit
from framework.utils.terminal_colors import TerminalColors as tc

from src.preprocessing import DataPreprocessor


@timeit
def preprocessing(cfgInputFile) -> DataPreprocessor:

    with open(cfgInputFile, 'r') as file:     cfgInput = yaml.safe_load(file)
    inFilePath = cfgInput['inFilePath']
    outFilePath = cfgInput['outFilePath']
    cfgVisualFile = cfgInput['visualFilePath']
    antimatterOnly = cfgInput.get('antimatterOnly', False)

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4table', dirPrefix='DF*')
    preprocessor = DataPreprocessor(dataHandler)
    preprocessor.define_variables()
    preprocessor.visualize(outFilePath, cfgVisualFile)

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
    args = parser.parse_args()

    if args.cfgInputFile == '':
        print(tc.RED+'[ERROR]: '+tc.RESET+'No config file provided, exiting.')
        exit(1)

    preprocessor = preprocessing(args.cfgInputFile)
