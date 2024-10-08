'''
    Run the studies
'''
import numpy as np
import argparse
import yaml

from studies.studies import StandaloneStudy
from studies.correlationStudies import CorrelationStudy
from studies.invMassStudies import InvariantMassStudy
from studies.clusterStudies import ClusterSizeParamStudy
from studies.betheBlochStudies import BetheBlochStudy
from studies.comparisonStudies import ComparisonStudy
from studies.TOFselectionStudies import TOFselectionStudy

import sys
sys.path.append('..')
from framework.src.hist_info import HistLoadInfo
from framework.utils.terminal_colors import TerminalColors as tc

def run_cluster_size_param_study(config):
    '''
        Run the cluster size param study
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Running the cluster size param study.')

    bb_config = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/config/cfg_BetheBloch.yml'
    h2_cluster_info = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/data_visual_selectionsHe3.root',
                                   'ClusterSize/ClusterSizevsPPr')
    study = ClusterSizeParamStudy(config, bb_config, h2_cluster_info)
    study.rebinx(2)
    study.drawBetheBloch('OLD_PARAMETERS_h2_cBB_Pr')    # fit the Bethe Bloch curve with the old parameters
    study.fitBetheBloch()
    study.print_results()

def run_tpc_calibration_study(config):
    '''
        Run the TPC calibration study
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Running the TPC calibration study.')

    cfg = yaml.safe_load(open(config, 'r'))
    sameEventLoad = HistLoadInfo(cfg['SEFileTPCcalibration'],
                                 'TPC/dEdXvsBetaGammaHe3')

    study = BetheBlochStudy(config, sameEventLoad)
    study.rebinx(2)
    study.fit_BetheBloch(xMinFit=0.75, xMaxFit=2.)
    study.draw()

def run_correlation_study(config) -> float:
    '''
        Run the correlation study
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Running the correlation study.')

    cfg = yaml.safe_load(open(config, 'r'))
    sameEventLoad = HistLoadInfo(cfg['SEFileCorrelation'],
                                 'Correlations/fKstar')
    mixedEventLoad = HistLoadInfo(cfg['MEFileCorrelation'],
                                  'Correlations/fKstar')

    study = CorrelationStudy(config, sameEvent=sameEventLoad, mixedEvent=mixedEventLoad)
    study.save('_prerebin')
    #bin_edges = np.array([0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0], dtype=np.float32)
    bin_edges = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0], dtype=np.float32)
    study.custom_binning(bin_edges)
    #study.rebin(5)
    #study.self_normalize()
    norm_factor = study.normalize(low=0.5, high=0.8)
    #norm_factor = 1.0
    study.correlation_function()
    study.save()
    study.produce_plot(cfg['figuresFolderPath']+'correlationFunction.pdf')

    return norm_factor

def run_comparison_study(config, scaling_factor=1.0):
    '''
        Run the comparison study
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Running the comparison study.')

    cfg = yaml.safe_load(open(config, 'r'))
    study = ComparisonStudy(config)
    print(scaling_factor)

    for var in cfg['comparisonVariables']:
        sameEventLoad = HistLoadInfo(cfg['SEFileComparison'], var)
        mixedEventLoad = HistLoadInfo(cfg['SEFileComparison'], var)
        study.load_same_event(sameEventLoad)
        study.load_mixed_event(mixedEventLoad)  
        study.self_normalize()  
        #study.scale_mixing(scaling_factor)
        study.save()

def run_invariant_mass_study(config):
    '''
        Run the invariant mass study
    '''

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Running the invariant mass study.')

    cfg = yaml.safe_load(open(config, 'r'))
    for opt in ['Matter', 'Anti']:
        sameEventLoad = HistLoadInfo(cfg['SEFileInvMass'],
                                     f'InvMass/InvMass{opt}Li')
        mixedEventLoad = HistLoadInfo(cfg['MEFileInvMass'],
                                      f'InvMass/InvMass{opt}Li')

        study = InvariantMassStudy(config, sameEvent=sameEventLoad, mixedEvent=mixedEventLoad, opt=opt)
        study.save('_prerebin')
        bin_edges = np.array([3.7, 3.74, 3.742, 3.744, 3.746, 3.748, 3.75, 3.752, 3.754, 3.756, 3.758, 3.76, 3.762, 3.764, 3.766, 3.768, 3.77, 3.772, 3.774, 3.776, 3.778, 3.78, 3.79, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.88, 3.9], dtype=np.float32)
        study.custom_binning(bin_edges)
        #study.rebin(2)
        study.self_normalize()
        study.normalize(low=3.8, high=3.89)
        study.bkg_subtraction()
        study.save()
        study.produce_plot(cfg['figuresFolderPath']+f'invariantMass{opt}.pdf')

def run_TOF_selection_study(config):

    print(tc.GREEN+'[INFO]: '+tc.RESET+'Running the TOF selection study.')

    prLoad = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/event_mixing_visual_selectionsHe3.root',
                         'TOF/TOFmassvsPtPr')
    study = TOFselectionStudy(config, prLoad)
    study.rebinx(2)
    study.fit_slices()
    study.fit_resolution()
    study.fit_mass()
    study.save()

def main():

    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--all', dest='run_all',
                        help='Run all the studies.', default=False, action='store_true')
    parser.add_argument('--correlation', dest='run_correlation',
                        help='Run the correlation study.', default=False, action='store_true')
    parser.add_argument('--invMass', dest='run_invMass',
                        help='Run the invariant mass study.', default=False, action='store_true')
    parser.add_argument('--clusterSize', dest='run_clusterSize',
                        help='Run the cluster size param study.', default=False, action='store_true')
    parser.add_argument('--comparison', dest='run_comparison',
                        help='Run the comparison study.', default=False, action='store_true')
    parser.add_argument('--TOF', dest='run_TOF',
                        help='Run the TOF selection study.', default=False, action='store_true')
    parser.add_argument('--TPCcalibration', dest='run_TPCcalibration',
                        help='Run the TPC calibration study.', default=False, action='store_true')
    args = parser.parse_args()

    #config = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/config/cfg_studies.yml'
    config = '/home/galucia/antiLithium4/analysis/config/cfg_studies.yml'
    norm_factor = 1.0
    if args.run_all or args.run_correlation:    norm_factor = run_correlation_study(config)
    if args.run_all or args.run_invMass:        run_invariant_mass_study(config)
    if args.run_all or args.run_clusterSize:    run_cluster_size_param_study(config)
    if args.run_all or args.run_TOF:            run_TOF_selection_study(config)
    if args.run_all or args.run_TPCcalibration: run_tpc_calibration_study(config)
    if args.run_all or args.run_comparison:     run_comparison_study(config, norm_factor)
    StandaloneStudy.close()

if __name__ == '__main__':

    main()