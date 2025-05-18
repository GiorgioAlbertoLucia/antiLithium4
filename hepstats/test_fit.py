from enum import Enum
import uproot
from hist import Hist, storage, axis
import boost_histogram as bh

from fit_workflow import FitWorkflowCF
from core.utils import match_hist_binning

from torchic.core.histogram import HistLoadInfo, load_hist
from torchic.core.histogram import HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc

import sys
sys.path.append('/home/galucia/antiLithium4/analysis')
from studies.correlationStudies import CorrelationStudy

from ROOT import TFile, TH1F, TH2F
import numpy as np

class Centrality(Enum):
    CENT_0_10 = 0
    CENT_10_30 = 1
    CENT_30_50 = 2
    INTEGRATED = 3

SUFFIX_DICT = {
    'dot': 
    [
        '_cent0.0_10.0',
        '_cent10.0_30.0',
        '_cent30.0_50.0',
        '',
    ],
    'undot':
    [
        '_cent0_10',
        '_cent10_30',
        '_cent30_50',
        '',
    ]
}

CORAL_DICT = [
    '4.2',
    '3.4',
    '2.6',
    '3.6',
]

def load_bkg_template(centrality_opt: int, model: str) -> Hist:
    '''
        Load the Coulomb template
    '''
    h_bkg = None
    if model == 'CATS':
        cent = SUFFIX_DICT['undot'][centrality_opt]
        h_bkg = uproot.open(f'/home/galucia/antiLithium4/analysis/output/CATS/CATS{cent}_new.root')\
                                                           ['hHe3_p_Coul_CF_LS']
    elif model == 'CorAL':
        h_bkg = uproot.open(f'/home/galucia/antiLithium4/analysis/output/CorAL/sqwell_correlation.root')\
                                                           [f'radius_{CORAL_DICT[centrality_opt]}fm/CF_{CORAL_DICT[centrality_opt]}fm']
    elif model == 'CorALCoul':
        h_bkg = uproot.open(f'/home/galucia/antiLithium4/analysis/output/CorAL/coulomb_correlation.root')\
                                                           [f'radius_{CORAL_DICT[centrality_opt]}fm/CF_{CORAL_DICT[centrality_opt]}fm']
        
    return h_bkg.to_hist()

def compute_correlation_functions(bin_edges: np.ndarray, sign: str) -> Hist:
    '''
        Compute the correlation function. 
        The correlation function is defined as:
            C(k*) = h_same_event / h_mixed_event

        Args:
            h_same_event_centrality (TH2F): Same event distribution for centrality vs k* 
            h_mixed_event_centrality (TH2F): Mixed event distribution for centrality vs k*
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Computing the correlation function.')
    h_same_event_centrality = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/data_visual_selectionsPr.root',
                                            f'Correlations/fKstarCentrality{sign}')
    h_mixed_event_centrality = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
                                            f'Correlations/fKstarCentrality{sign}')

    study = CorrelationStudy(config=None, outputFile=None, 
                             sameEvent=h_same_event_centrality,
                             mixedEvent=h_mixed_event_centrality, 
                             opt='anti')

    study.custom_binning(bin_edges)
    study.normalize(low=0.4, high=0.8)

    CENTRALITY_BIN_EDGES = np.array([0., 10., 30., 50.], dtype=np.float32)
    study.CentralityBinEdges = CENTRALITY_BIN_EDGES
    study.correlation_function_centrality(low_value_norm=0.05, high_value_norm=0.75, bin_edges=bin_edges)
    h_correlation_functions = []

    for centrality in Centrality:
        if centrality == Centrality.INTEGRATED:
            th1_correlation_function = study.hCorrelation
        else:
            th1_correlation_function = study.hCorrelationCent[centrality.value]

        nbins   = th1_correlation_function.GetNbinsX()  
        edges   = [th1_correlation_function.GetBinLowEdge(1) + i * th1_correlation_function.GetBinWidth(1) for i in range(nbins+1)]
        values    = [th1_correlation_function.GetBinContent(i) for i in range(1, nbins+1)]
        variances = [th1_correlation_function.GetBinError(i)**2    for i in range(1, nbins+1)]

        h_correlation_function = Hist(
            axis.Regular(nbins, edges[0], edges[-1], name=th1_correlation_function.GetName()),
            storage=storage.Weight()
        )
        h_correlation_function[...] = np.stack([values, variances], axis=-1)
        h_correlation_functions.append(h_correlation_function)

    return h_correlation_functions


if __name__ == "__main__":

    sign = 'Anti'
    h_correlation_functions = []

    with uproot.recreate('output/test_fit.root') as file:
        for centrality in Centrality:

            print(f'{centrality.value=}')
            print(f'{type(centrality.value)=}')

            cent = SUFFIX_DICT['undot'][centrality.value]
            cent_dot = SUFFIX_DICT['dot'][centrality.value]
            fit_workflow = FitWorkflowCF(outfile=file, pdf_outpath=f'output/plots{cent}.pdf')
            fit_workflow.add_directory(f'dir{cent}')

            GLOBAL_REBIN = 2
            LOW_VALUE_KSTAR = 0.0j     # lower and upper limits of the k* range
            HIGH_VALUE_KSTAR = 0.80j

            h_signal_mc = uproot.open('/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root')\
                                    [f'Correlations/fKstar{sign}'].to_hist()\
                                    [LOW_VALUE_KSTAR:HIGH_VALUE_KSTAR:]
            h_mixed_event = uproot.open('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root')\
                                    [f'Correlations/fKstar{sign}'].to_hist()\
                                    [LOW_VALUE_KSTAR:HIGH_VALUE_KSTAR:]
            h_bkg_template = load_bkg_template(centrality.value, 'CATS')
            h_bkg_template = h_bkg_template[LOW_VALUE_KSTAR:HIGH_VALUE_KSTAR:bh.rebin(5*GLOBAL_REBIN)] # rebin

            bin_edges = h_bkg_template.axes[0].edges
            if len(h_correlation_functions) < centrality.value + 1:
                h_correlation_functions = compute_correlation_functions(bin_edges=bin_edges, sign=sign)
            fit_workflow.load_hist_CF(h_correlation_functions[centrality.value])

            fit_workflow.prepare_bkg_template(h_bkg_template)
            fit_workflow.prefit_signal(h_signal_mc, h_mixed_event)
            fit_workflow.fit_CF(get_upper_limit=False,
                                get_pnull_significance=True,
                                get_pnull_significance_per_bin=False,)
            
            fit_workflow.clear_data()
