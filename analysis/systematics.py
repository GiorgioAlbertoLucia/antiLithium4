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
    add_condition = True
    if 'fEta' in var_name:  condition = (abs(dataset[var_name]) < var_value)
    if 'fNSigmaTPC' in var_name:  condition = (abs(dataset[var_name]) < var_value)
    if 'fNSigmaTOFHad' in var_name:  condition = (abs(dataset[var_name]) < var_value) & (abs(dataset['fPtHad']) > 0.8)

    return condition


def systematics(args):

    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.BOLD+tc.WHITE+'\t\tRunning analysis'+tc.RESET)
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET)

    with open(args.cfgInputFile, 'r') as file:     
        cfgInput = yaml.safe_load(file)
    SE_file_path = cfgInput['SEFilePath']
    SE_tree_name = cfgInput['SETreeName']
    SE_folder_name = cfgInput.get('SEFolderName', None)
    ME_file_path = cfgInput['MEFilePath']
    ME_tree_name = cfgInput['METreeName']
    ME_folder_name = cfgInput.get('MEFolderName', None)

    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'\t\tCut variation:')
    # cut variation dictionary: {var_name: [var_min, var_max, n_steps]}
    cut_variation_dict_cfg = cfgInput['cutVariationDict']
    cut_variation_dict = {}
    for var_name, var_range in cut_variation_dict_cfg.items():
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'{var_name}: {var_range}')
        var_min, var_max, n_steps = var_range
        cut_variation_dict[var_name] = np.linspace(var_min, var_max, n_steps)
    print(tc.GREEN+'[INFO]: '+tc.RESET)

    dataset_SE = None
    for tree_name in SE_tree_name:
        if dataset_SE is None:
            dataset_SE = Dataset.from_root(SE_file_path, folder_name=SE_folder_name, tree_name=tree_name)
        else:
            dataset_SE = dataset_SE.concat(Dataset.from_root(SE_file_path, folder_name=SE_folder_name, tree_name=tree_name), axis=1)
    dataset_SE.query('fPIDtrkHe3 != 6', inplace=True)
    process_SE = DataPreprocessor(dataset_SE)
    process_SE.define_variables()
    process_SE.define_nsigmaTOF_Pr()
    process_SE.define_kstar()
    dataset_SE = process_SE.dataset
    columns_to_keep = [ 'fPtHad',
                        'fEtaHe3',
                        'fEtaHad',
                        'fZVertex',
                        'fDCAxyHe3',
                        'fDCAzHe3',
                        'fDCAxyHad',
                        'fDCAzHad',
                        'fNSigmaTPCHe3',
                        'fNSigmaTPCHad',
                        'fNSigmaTOFHad',
                        'fCentralityFT0C',
                        'fKstar',
                        'fMassInvLi']
    columns_to_drop_SE = [col for col in dataset_SE.columns if col not in columns_to_keep]
    dataset_SE.drop(columns=columns_to_drop_SE, axis=1, inplace=True)
    print(tc.GREEN+'[INFO]: '+tc.RESET+f'{dataset_SE.columns=}')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')

    dataset_ME = None
    for tree_name in ME_tree_name:
        if dataset_ME is None:
            dataset_ME = Dataset.from_root(ME_file_path, tree_name=tree_name)
        else:
            dataset_ME = dataset_ME.concat(Dataset.from_root(ME_file_path, tree_name=tree_name), axis=1)
    dataset_ME.query('fPIDtrkHe3 != 6', inplace=True)
    process_ME = DataPreprocessor(dataset_ME)
    process_ME.define_variables()
    process_ME.define_nsigmaTOF_Pr()
    process_ME.define_kstar()
    dataset_ME = process_ME.dataset
    columns_to_drop_ME = [col for col in dataset_ME.columns if col not in columns_to_keep]
    dataset_ME.drop(columns=columns_to_drop_ME, inplace=True)
    print(tc.GREEN+'[INFO]: '+tc.RESET+f'{dataset_ME.columns=}')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'------------------------------------------------------')

    hist_corr_funcs = []
    corr_study = CorrelationStudy()

    N_VARIATIONS = cfgInput.get('nVariations', 100)
    with alive_bar(N_VARIATIONS, title=tc.GREEN+'[INFO]: '+tc.RESET) as bar:
        for iter in range(N_VARIATIONS):
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