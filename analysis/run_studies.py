'''
    Run the studies
'''
import os
import numpy as np
import argparse
import yaml

from ROOT import TFile

from studies.correlationStudies import CorrelationStudy
from studies.invMassStudies import InvariantMassStudy

from torchic import HistLoadInfo
from torchic.core.histogram import load_hist
from torchic.utils.terminal_colors import TerminalColors as tc


# temporary
def add_hist(hist1, hist2):
    assert hist1.GetNbinsX() == hist2.GetNbinsX()
    hist_sum = hist1.Clone(hist1.GetName()+'_sum')
    for ibin in range(1, hist1.GetNbinsX()+1):
        hist_sum.SetBinContent(ibin, hist1.GetBinContent(ibin) + hist2.GetBinContent(ibin))
        hist_sum.SetBinError(ibin, np.sqrt(hist1.GetBinError(ibin)**2 + hist2.GetBinError(ibin)**2))
    return hist_sum

def run_correlation_study(config, outputFile:TFile) -> float:
    '''
        Run the correlation study
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET)
    print(tc.GREEN+'[INFO]: '+tc.RESET+'----------------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'\t\tRunning the correlation study.')
    print(tc.GREEN+'[INFO]: '+tc.RESET)

    cfg = yaml.safe_load(open(config, 'r'))
    
    #for opt in ['Matter', 'Anti', 'Both', 'US']:
    for opt in ['Matter', 'Anti']:
        print(tc.GREEN+f'[INFO]: {opt}'+tc.RESET)
        if opt == 'Both':
            hSameEventMatter = load_hist(HistLoadInfo(cfg['SEFileCorrelation'], 'Correlations/fKstarMatter'))
            hSameEventAnti = load_hist(HistLoadInfo(cfg['SEFileCorrelation'], 'Correlations/fKstarAnti'))
            sameEventLoad = add_hist(hSameEventMatter, hSameEventAnti)
            print('type same event:', type(sameEventLoad))

            hMixedEventMatter = load_hist(HistLoadInfo(cfg['MEFileCorrelation'], 'Correlations/fKstarMatter'))
            hMixedEventAnti = load_hist(HistLoadInfo(cfg['MEFileCorrelation'], 'Correlations/fKstarAnti'))
            mixedEventLoad = add_hist(hMixedEventMatter, hMixedEventAnti)
            print('type mixed event:', type(mixedEventLoad))
        else:
            sameEventLoad = HistLoadInfo(cfg['SEFileCorrelation'],
                                        f'Correlations/fKstarCentrality{opt}')
            mixedEventLoad = HistLoadInfo(cfg['MEFileCorrelation'],
                                        f'Correlations/fKstarCentrality{opt}')
            sameEventMassTLoad = HistLoadInfo(cfg['SEFileCorrelation'],
                                        f'Correlations/fKstarMassT{opt}')
            mixedEventMassTLoad = HistLoadInfo(cfg['MEFileCorrelation'],
                                        f'Correlations/fKstarMassT{opt}')
            #sameEventLoad = HistLoadInfo(cfg['SEFileCorrelation'],
            #                f'Correlations/fKstar{opt}')
            #mixedEventLoad = HistLoadInfo(cfg['MEFileCorrelation'],
            #                            f'Correlations/fKstar{opt}')
        study = CorrelationStudy(config, outputFile, sameEvent=sameEventLoad, mixedEvent=mixedEventLoad, opt=opt)
        study.save('_prerebin')
        #bin_edges = np.array([0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0], dtype=np.float32)
        #bin_edges = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0], dtype=np.float32)
        bin_edges = np.linspace(0.0, 0.8, 61)
        study.custom_binning(bin_edges)
        #study.rebin(5)
        #study.self_normalize()
        norm_factor = study.normalize(low=0.05, high=0.75)
        #norm_factor = 1.0
        # study.correlation_function()

        #study.correlation_function_massT(sameEvent=sameEventMassTLoad, mixedEvent=mixedEventMassTLoad, low_value_norm=0.05, high_value_norm=0.75, bin_edges=bin_edges)

        if 'Centralities' in cfg:
            centrality_bin_edges = np.array(cfg['Centralities']['CentralityBinEdges'], dtype=np.float32)
            study.CentralityBinEdges = centrality_bin_edges
            genuine_loads = []
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Running centrality study. Centrality classes are: ', end='')
            for icent in range(len(centrality_bin_edges)-1):
                print(f'[{centrality_bin_edges[icent]} - {centrality_bin_edges[icent+1]}]%, ', end=' ')
                if 'GenuineFileCorrelation' in cfg and 'GenuineNameCorrelation' in cfg:
                    genuine_file = os.path.splitext(cfg['GenuineFileCorrelation'])[0] + f'_cent{int(centrality_bin_edges[icent])}_{int(centrality_bin_edges[icent+1])}' + '.root'
                    genuine_loads.append(HistLoadInfo(genuine_file, cfg['GenuineNameCorrelation']))
            print()

            study.correlation_function_centrality(low_value_norm=0.05, high_value_norm=0.75, bin_edges=bin_edges)
            if 'GenuineFileCorrelation' in cfg and 'GenuineNameCorrelation' in cfg:
                study.pull_distribution_centrality(genuine_loads)
            
        
        if 'GenuineFileCorrelation' in cfg and 'GenuineNameCorrelation' in cfg:
            genuine_load = HistLoadInfo(cfg['GenuineFileCorrelation'], cfg['GenuineNameCorrelation'])
        else:
            genuine_load = None
        study.pull_distribution(genuine_load)


        study.save()
        study.produce_plot(cfg['figuresFolderPath']+f'correlationFunction{opt}.pdf')
    
    print(tc.GREEN+'[INFO]: '+tc.RESET+'----------------------------------------------------------------')

    return norm_factor

def run_invariant_mass_study(config, outputFile:TFile):
    '''
        Run the invariant mass study
    '''
    print(tc.GREEN+'[INFO]: '+tc.RESET)
    print(tc.GREEN+'[INFO]: '+tc.RESET+'----------------------------------------------------------------')
    print(tc.GREEN+'[INFO]: '+tc.RESET+'\t\tRunning the invariant mass study.')
    print(tc.GREEN+'[INFO]: '+tc.RESET)

    cfg = yaml.safe_load(open(config, 'r'))
    for opt in ['Matter', 'Anti']:
        sameEventLoad = HistLoadInfo(cfg['SEFileInvMass'],
                                     f'InvMass/InvMassCentrality{opt}')
        mixedEventLoad = HistLoadInfo(cfg['MEFileInvMass'],
                                      f'InvMass/InvMassCentrality{opt}')
        #correctionLoad = HistLoadInfo(cfg['CorrectionFileInvMass'],
        #                                f'hHe3_p_Coul_InvMass')
        #hCorrection = load_hist(correctionLoad)

        study = InvariantMassStudy(config, outputFile, sameEvent=sameEventLoad, mixedEvent=mixedEventLoad, opt=opt)
        study.save('_prerebin')
        #bin_edges = np.array([3.747, 3.749, 3.751, 3.753, 3.755, 3.757, 3.759, 3.761, 3.763, 3.765, 3.767, 3.769, 3.771, 3.773, 3.775, 3.777, 3.779, 3.784, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85], dtype=np.float32)#, 3.86, 3.88, 3.9], dtype=np.float32)
        #bin_edges = np.array([3.747, 3.749, 3.751, 3.753, 3.755, 3.757, 3.759, 3.761, 3.763, 3.765, 3.767, 3.769, 3.771, 3.773, 3.775, 3.777, 3.779, 3.784, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.88, 3.9], dtype=np.float32)
        #bin_edges = np.array([3.747, 3.751, 3.755, 3.759, 3.763, 3.767, 3.771, 3.775, 3.779, 3.784, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85], dtype=np.float32)
        #study.custom_binning(bin_edges)
        study.rebin(4)
        #study.self_normalize()

        if 'Centralities' in cfg:
            centrality_bin_edges = np.array(cfg['Centralities']['CentralityBinEdges'], dtype=np.float32)
            study.CentralityBinEdges = centrality_bin_edges
            correction_hists = []
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Running centrality study. Centrality classes are: ', end='')
            for icent in range(len(centrality_bin_edges)-1):
                print(f'[{centrality_bin_edges[icent]} - {centrality_bin_edges[icent+1]}]%, ', end=' ')
                if 'CorrectionFileInvMass' in cfg and 'CorrectionNameInvMass' in cfg:
                    # correct, needs the invariant mass added to the centrality file
                    correction_file = os.path.splitext(cfg['CorrectionFileInvMass'])[0] + f'_cent{int(centrality_bin_edges[icent])}_{int(centrality_bin_edges[icent+1])}' + '_new.root'
                    #correction_file = cfg['CorrectionFileInvMass']
                    correction_hists.append(load_hist(HistLoadInfo(correction_file, cfg['CorrectionNameInvMass'])))
            print()

            if 'CorrectionFileInvMass' in cfg and 'CorrectionNameInvMass' in cfg:
                study.set_corrections_centrality(correction_hists)
            study.invariant_mass_centrality(low_value_norm=3.78, high_value_norm=3.799, bin_edges=[], rebin_factor=2)

        #study.correct_mixed_event(hCorrection)
        study.normalize(low=3.78, high=3.84)
        study.bkg_subtraction()
        study.pull_distribution(low_edge_fit=None, high_edge_fit=3.84)
        study.ratio_distribution()
        study.save()
        study.produce_plot(cfg['figuresFolderPath']+f'invariantMass{opt}.pdf')
    
    print(tc.GREEN+'[INFO]: '+tc.RESET+'----------------------------------------------------------------')

def main():

    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--all', dest='run_all',
                        help='Run all the studies.', default=False, action='store_true')
    parser.add_argument('--correlation', dest='run_correlation',
                        help='Run the correlation study.', default=False, action='store_true')
    parser.add_argument('--invMass', dest='run_invMass',
                        help='Run the invariant mass study.', default=False, action='store_true')
    args = parser.parse_args()

    #config = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/config/cfg_studies.yml'
    config = '/home/galucia/antiLithium4/analysis/config/cfg_studies.yml'
    norm_factor = 1.0
    cfg = yaml.safe_load(open(config, 'r'))
    outputFile = TFile.Open(cfg['studiesOutputFilePath'], 'recreate')

    if args.run_all or args.run_correlation:        norm_factor = run_correlation_study(config, outputFile)
    if args.run_all or args.run_invMass:            run_invariant_mass_study(config, outputFile)
    
    outputFile.Close()

if __name__ == '__main__':

    main()