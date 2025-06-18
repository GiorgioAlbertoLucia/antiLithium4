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

    with open(args.cfg_input_file, 'r') as file:     
        cfgInput = yaml.safe_load(file)

    tree_names = cfgInput['treeNames']
    folder_name = cfgInput.get('folderName', None)

    dataset = None
    for treeName in tree_names:
        if dataset is None:
            dataset = Dataset.from_root(cfgInput['inFilePath'], tree_name=treeName, folder_name=folder_name)
        else:
            dataset = dataset.concat(Dataset.from_root(cfgInput['inFilePath'], tree_name=treeName, folder_name=folder_name), axis=1)
    
    print(tc.GREEN+'[INFO]: '+tc.RESET+f'{dataset.columns=}')
    dataset['fItsClusterSizeHe3'] = np.array(dataset['fItsClusterSizeHe3'], dtype=np.int32)
    dataset['fItsClusterSizeHad'] = np.array(dataset['fItsClusterSizeHad'], dtype=np.int32)

    if args.us_only: 
        dataset.query('(fPtHe3 > 0 and fPtHad < 0) or (fPtHe3 < 0 and fPtHad > 0)', inplace=True)
    if args.ls_only: 
        dataset.query('(fPtHe3 > 0 and fPtHad > 0) or (fPtHe3 < 0 and fPtHad < 0)', inplace=True)
        #dataset.query('fIsBkgUS == 0', inplace=True)
    preprocessor = DataPreprocessor(dataset)
    preprocessor.define_variables()
    preprocessor.define_nsigmaTPC_He3()

    for selection in cfgInput['selections']:
        preprocessor.apply_cut(selection)

    #if args.qa: 
    #    preprocessor.define_nsigmaTOF_Pr()
    #    preprocessor.visualize(outQaFilePath.replace('.root', '_purity.root'), cfgInput['qaFilePath'])
    
    preprocessor.define_kstar()
    if args.save_df:
        preprocessor.save_df(cfgInput['dfFilePath'], ['fKstar', 'fCentralityFT0C', 'fIsMatter'])
        return preprocessor
    #preprocessor.visualize(cfgInput['outFilePath'], cfgInput['visualFilePath'])
    preprocessor.visualize_boost(cfgInput['outFilePath'], cfgInput['visualFilePath'])
    if args.qa: preprocessor.visualize(cfgInput['outQaFilePath'], cfgInput['qaFilePath'])
    
    return preprocessor


if __name__ == '__main__':

    print()
    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--config-file', dest='cfg_input_file',
                        help='path to the YAML file with configuration.', default='')
    parser.add_argument('--QA', dest='qa', help='Draw QA plots.', action='store_true')
    parser.add_argument('--USonly', dest='us_only', help='Flag to analyse only US pairs.', action='store_true')
    parser.add_argument('--LSonly', dest='ls_only', help='Flag to analyse only LS pairs.', action='store_true')
    parser.add_argument('--saveDf', dest='save_df', help='Flag for close pair rejection.', action='store_true')
    args = parser.parse_args()

    if args.cfg_input_file == '':
        print(tc.RED+'[ERROR]: '+tc.RESET+'No config file provided, exiting.')
        exit(1)

    preprocessor = preprocessing(args)
