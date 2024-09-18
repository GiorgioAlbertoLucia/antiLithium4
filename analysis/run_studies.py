'''
    Run the studies
'''
import numpy as np
from argparse import ArgumentParser

from studies.studies import StandaloneStudy
from studies.correlationStudies import CorrelationStudy
from studies.invMassStudies import InvariantMassStudy
from studies.clusterStudies import ClusterSizeParamStudy

import sys
sys.path.append('..')
from framework.src.hist_info import HistLoadInfo

def run_cluster_size_param_study(args, config):
    '''
        Run the cluster size param study
    '''

    bb_config = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/src/config/cfgBetheBloch.yml'
    h2_cluster_info = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/data_visual_selectionsHe3.root',
                                   'ClusterSize/ClusterSizevsPPr')
    study = ClusterSizeParamStudy(config, bb_config, h2_cluster_info)
    study.rebinx(2)
    study.drawBetheBloch('OLD_PARAMETERS_h2_cBB_Pr')    # fit the Bethe Bloch curve with the old parameters
    study.fitBetheBloch()
    study.print_results()

def run_correlation_study(args, config):
    '''
        Run the correlation study
    '''

    sameEventLoad = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/data_visual_selectionsHe3.root',
                                 'Correlations/fKstar')
    mixedEventLoad = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/event_mixing_visual_selectionsHe3.root',
                                  'Correlations/fKstar')

    study = CorrelationStudy(config, sameEvent=sameEventLoad, mixedEvent=mixedEventLoad)
    study.save('_prerebin')
    bin_edges = np.array([0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0], dtype=np.float32)
    study.custom_binning(bin_edges)
    #study.rebin(5)
    study.normalize()
    study.correlation_function()
    study.save()
    study.produce_plot('/Users/glucia/Projects/ALICE/antiLithium4/analysis/figures/correlationFunction.pdf')

def run_invariant_mass_study(args, config):
    '''
        Run the invariant mass study
    '''

    sameEventLoad = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/data_visual_selectionsHe3.root',
                                 'InvMass/InvMassLi')
    mixedEventLoad = HistLoadInfo('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24/event_mixing_visual_selectionsHe3.root',
                                  'InvMass/InvMassLi')

    study = InvariantMassStudy(config, sameEvent=sameEventLoad, mixedEvent=mixedEventLoad)
    study.normalize()
    study.bkg_subtraction()
    study.save()
    #study.produce_plot('/Users/glucia/Projects/ALICE/antiLithium4/analysis/figures/invariantMass.pdf')

def main():

    config = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/src/config/cfgStudies.yml'
    run_correlation_study(None, config)
    run_invariant_mass_study(None, config)
    run_cluster_size_param_study(None, config)
    StandaloneStudy.close()

if __name__ == '__main__':

    main()