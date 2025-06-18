'''
    Macro to fit the p-He3 correlation function
'''

import argparse
from torchic.core.histogram import HistLoadInfo, load_hist
from torchic.utils.terminal_colors import TerminalColors as tc

from ROOT import TFile, TH1F

import sys
sys.path.append('..')
from femto.fit_workflow_cf import FitWorkflowCF

SUFFIX_DICT = {
    'dot':
    {
        'integrated': '',
        '0-10': '_cent0.0_10.0',
        '10-30': '_cent10.0_30.0',
        '30-50': '_cent30.0_50.0',
    },
    'undot':
    {
        'integrated': '',
        '0-10': '_cent0_10',
        '10-30': '_cent10_30',
        '30-50': '_cent30_50',
    }
}

CORAL_DICT = {
    'integrated': '3.6',
    '0-10': '4.2',
    '10-30': '3.4',
    '30-50': '2.6',

}

def multiply_hist(hist_to_multiply: TH1F, hist: TH1F) -> TH1F:
    '''
        Multiply two histograms
    '''
    _hist_to_multiply = hist_to_multiply.Clone()
    
    for ibin in range(1, hist_to_multiply.GetNbinsX() + 1):
        original_bin_content = _hist_to_multiply.GetBinContent(ibin)
        original_bin_error = _hist_to_multiply.GetBinError(ibin)
        bin_to_multiply = hist.FindBin(_hist_to_multiply.GetBinCenter(ibin))
        value_to_multiply = hist.GetBinContent(bin_to_multiply)
        error_to_multiply = hist.GetBinError(bin_to_multiply)

        bin_content = original_bin_content * value_to_multiply
        bin_error = 0
        if original_bin_content > 0 and value_to_multiply > 0:
            bin_error = bin_content * (original_bin_error / original_bin_content + error_to_multiply / value_to_multiply)

        hist_to_multiply.SetBinContent(ibin, bin_content)
        hist_to_multiply.SetBinError(ibin, bin_error)
    del _hist_to_multiply
    return hist_to_multiply

def load_bkg_template(centrality_opt: str, model: str) -> TH1F:
    '''
        Load the Coulomb template
    '''
    h_coulomb = None
    if model == 'CATS':
        cent = SUFFIX_DICT['undot'][centrality_opt]
        h_coulomb = load_hist(HistLoadInfo(f'/home/galucia/antiLithium4/analysis/output/CATS/CATS{cent}_new.root', 
                                                           'hHe3_p_Coul_CF_LS'))
    elif model == 'CorAL':
        h_coulomb = load_hist(HistLoadInfo(f'/home/galucia/antiLithium4/analysis/output/CorAL/sqwell_correlation.root', 
                                                           f'radius_{CORAL_DICT[centrality_opt]}fm/CF_{CORAL_DICT[centrality_opt]}fm'))
    elif model == 'CorALCoul':
        h_coulomb = load_hist(HistLoadInfo(f'/home/galucia/antiLithium4/analysis/output/CorAL/coulomb_correlation.root', 
                                                           f'radius_{CORAL_DICT[centrality_opt]}fm/CF_{CORAL_DICT[centrality_opt]}fm'))
    return h_coulomb

def init_parser():
    '''
        Initialize the parser for the command line arguments
    '''
    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--sign', dest='sign',
                        help='Perform the fit on either matter (Matter) or antimatter (Anti).', default='Anti')
    parser.add_argument('--model', dest='model',
                        help='Model to use for the Coulomb template.', default='CATS')
    parser.add_argument('--kstar', dest='kstar',
                        help='Fit the kstar distribution instead of the correlation function?', default=False, action='store_true')
    return parser

def main(workflow: FitWorkflowCF, centrality: str):

    cent = SUFFIX_DICT['undot'][centrality]
    cent_dot = SUFFIX_DICT['dot'][centrality]
    workflow.make_directory(f'dir{cent}')

    h_signal =      load_hist(HistLoadInfo(
                    '/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root',
                    f'Correlations/fKstar{args.sign}'
                    ))
    h_mixed_event = load_hist(HistLoadInfo(
                    '/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
                    f'Correlations/fKstar{args.sign}'
                    ))
    if not args.kstar: 
        h_signal.Divide(h_mixed_event)
    workflow.prepare_signal_fit(h_signal)

    h_bkg = load_bkg_template(centrality, args.model)
    if args.kstar:
        h_bkg = multiply_hist(h_bkg, h_mixed_event)
    workflow.prepare_bkg_template(h_bkg)

    # p-He3 correlation function
    data_name = 'Same' if args.kstar else 'Correlation'
    print(f'Correlation{args.sign}/h{data_name}_kstar{cent_dot}')
    h_data =    load_hist(HistLoadInfo(
                '/home/galucia/antiLithium4/analysis/output/PbPb/studies.root', 
                f'Correlation{args.sign}/h{data_name}_kstar{cent_dot}'
                ))
    workflow.load_hist_data(h_data)
    
    workflow.pull_from_bkg(h_bkg)
    logy = args.kstar
    workflow.fit_CF(logy)
    workflow.integral_CF()

    workflow.evaluate_pvalue_and_significance(cent)
    if args.kstar:
        workflow.evaluate_pvalue_and_significance_counts(cent)
    workflow.scan_pvalue_and_significance()
    workflow.clear_data()

if __name__ == '__main__':

    parser = init_parser()
    args = parser.parse_args()
    mode = 'kstar' if args.kstar else 'CF'

    outfile = TFile(f'output/{args.sign}Lithium4Fit{mode}{args.model}.root', 'RECREATE')
    pdf_path = f'output/{args.sign}Lithium4Fit{mode}{args.model}.pdf'
    workflow = FitWorkflowCF(outfile, pdf_path)

    for centrality_opt in ['integrated', '0-10', '10-30', '30-50']:
        print(f'\n\n{tc.BOLD+tc.GREEN}{centrality_opt=}{tc.RESET}')
        main(workflow, centrality_opt)

    workflow.close_canvas()
    outfile.Close()
