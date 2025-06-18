'''
    Analysis routine to study the effects of systematics in the correlation function using the cut variation technique.
'''

import numpy as np
import yaml
import argparse
from alive_progress import alive_bar

from ROOT import TFile

from torchic import Dataset, AxisSpec
from torchic.utils import timeit
from torchic.utils import TerminalColors as tc

import sys
sys.path.append('/home/galucia/antiLithium4/analysis')
from src.preprocessing import DataPreprocessor
from studies.correlationStudies import CorrelationStudy

def add_variation(condition: bool, var_name: str, var_value, dataset: Dataset) -> bool:
    condition = True
    if 'fEta' in var_name:  condition = (abs(dataset[var_name]) < var_value)
    if 'fNSigmaTPC' in var_name:  condition = (abs(dataset[var_name]) < var_value)
    if 'fNSigmaTOFHad' in var_name:  condition = (abs(dataset[var_name]) < var_value) & (abs(dataset['fPtHad']) > 0.8)

    return condition

def load_cut_dict(cfgInput):
    '''
        Load the cut variation dictionary from the configuration file.
        The dictionary is expected to have the following structure:
        {
            'cut_name': [min_value, max_value, n_steps],'
            ...
        }
    '''

    cut_variation_dict_cfg = cfgInput['cutVariationDict']
    cut_variation_dict = {}
    for var_name, var_range in cut_variation_dict_cfg.items():
        var_min, var_max, n_steps = var_range
        cut_variation_dict[var_name] = np.linspace(var_min, var_max, n_steps)
    return cut_variation_dict

def _load_from_file(file_path, tree_names, folder_name=None):
    '''
        Load the dataset from a ROOT file.
        If the folder_name is provided, it will load the dataset from the specified folder.
        If the folder_name is None, it will load the dataset from the root of the file.
    '''
    dataset = None
    for tree_name in tree_names:
        if dataset is None:
            dataset = Dataset.from_root(file_path, folder_name=folder_name, tree_name=tree_name)
        else:
            dataset = dataset.concat(Dataset.from_root(file_path, folder_name=folder_name, tree_name=tree_name), axis=1)
    return dataset

def _drop_columns(dataset, columns_to_keep):
    '''
        Drop the columns that are not in the columns_to_keep list.
    '''
    columns_to_drop = [col for col in dataset.columns if col not in columns_to_keep]
    dataset.drop(columns=columns_to_drop, inplace=True)
    return dataset

def load_data(cfgInput):
    '''
        Load the data from the ROOT files and preprocess it.
    '''

    dataset_SE = _load_from_file(cfgInput['SEFilePath'], cfgInput['SETreeName'], cfgInput.get('SEFolderName', None))
    
    process_SE = DataPreprocessor(dataset_SE)
    process_SE = DataPreprocessor(dataset_SE)
    process_SE.define_variables()
    process_SE.define_nsigmaTPC_He3()
    for selection in cfgInput['selections']:
        process_SE.apply_cut(selection)
    
    dataset_SE = process_SE.dataset
    dataset_SE = _drop_columns(dataset_SE, cfgInput['columnsToKeep'])
    print(tc.GREEN+'[INFO]: '+tc.RESET+f'{dataset_SE.columns=}')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')

    dataset_ME = _load_from_file(cfgInput['MEFilePath'], cfgInput['METreeName'], cfgInput.get('MEFolderName', None))
    
    process_ME = DataPreprocessor(dataset_ME)
    process_ME.define_variables()
    process_ME.define_nsigmaTPC_He3()
    for selection in cfgInput['selections']:
        process_ME.apply_cut(selection)
    
    dataset_ME = process_ME.dataset
    dataset_ME = _drop_columns(dataset_ME, cfgInput['columnsToKeep'])
    print(tc.GREEN+'[INFO]: '+tc.RESET+f'{dataset_ME.columns=}')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')

    return dataset_SE, dataset_ME

def run_systematic_iteration(dataset_SE, dataset_ME, cut_variation_dict, iter, corr_study):
    '''
        Run a single iteration of the systematic analysis.
        This function is called for each variation of the cuts.
    '''
    condition_SE = True
    condition_ME = True
    for var_name, cut_values in cut_variation_dict.items():
        cut_value = np.random.choice(cut_values)
        condition_SE = condition_SE & add_variation(condition_SE, var_name, cut_value, dataset_SE)
        condition_ME = condition_ME & add_variation(condition_ME, var_name, cut_value, dataset_ME)

    dataset_SE.add_subset(f'cut_{iter}', condition_SE)
    dataset_ME.add_subset(f'cut_{iter}', condition_ME)

    hist_SE = dataset_SE.build_th1('fKstar', 
                                    AxisSpec(40, 0.0, 0.8, f'fKstar_SE_{iter}', '; #it{k}* (GeV/#it{c}); counts'),
                                    subset=f'cut_{iter}')
    hist_ME = dataset_ME.build_th1('fKstar', 
                                    AxisSpec(40, 0.0, 0.8, f'fKstar_ME_{iter}', '; #it{k}* (GeV/#it{c}); counts'),
                                    subset=f'cut_{iter}')

    corr_study.set_same_event(hist_SE)
    corr_study.set_mixed_event(hist_ME)
    corr_study.normalize(low=0.05, high=0.75)
    hist_corr = corr_study.correlation_function()

    return hist_corr

def systematics(args):

    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.BOLD+tc.WHITE+'\t\tRunning systematics'+tc.RESET)
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET)

    with open(args.cfgInputFile, 'r') as file:     
        cfgInput = yaml.safe_load(file)
    

    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'\t\tCut variation:')
    
    cut_variation_dict = load_cut_dict(cfgInput)
    dataset_SE, dataset_ME = load_data(cfgInput)

    hist_corr_funcs = []
    corr_study = CorrelationStudy()

    N_VARIATIONS = cfgInput.get('nVariations', 100)
    with alive_bar(N_VARIATIONS, title=tc.GREEN+'[INFO]: '+tc.RESET) as bar:
        for iter in range(N_VARIATIONS):
            
            hist_corr = run_systematic_iteration(dataset_SE, dataset_ME, cut_variation_dict, iter, corr_study)
            hist_corr_funcs.append(hist_corr)
            bar()

    point_variations = {}   # dict of lists with the values for each point of the CF for each cut
                            # the systematic on the point is evaluated as the RMS on the corresponding histogram
    
    for ibin in range(1, hist_corr_funcs[0].GetNbinsX()+1):
        for hist in hist_corr_funcs:
            point_variations[ibin] = point_variations.get(ibin, []) + [hist.GetBinContent(ibin)]

    hist_systematics = hist_corr_funcs[0].Clone('hist_systematics')
    hist_systematics.Reset()
    for ibin in range(1, hist_systematics.GetNbinsX()+1):
        point_systematics = np.std(point_variations[ibin])
        hist_systematics.SetBinContent(ibin, point_systematics)

    outfile_path = cfgInput['outFilePath']
    outfile = TFile.Open(outfile_path, 'RECREATE')

    hist_systematics.Write()
    corr_dir = outfile.mkdir('correlation_functions')
    corr_dir.cd()
    for hist in hist_corr_funcs:
        hist.Write()
    outfile.Close()

    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Systematics analysis')
    parser.add_argument('--config-file', dest='cfgInputFile', type=str, help='Path to the configuration file')
    args = parser.parse_args()

    systematics(args)