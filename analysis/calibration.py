'''
    Code to run calibration of ITS and TPC parametrisations
'''

import numpy as np
from ROOT import TFile
from ROOT import RooRealVar, RooCrystalBall, RooExponential
from torchic import Dataset, AxisSpec

from torchic.physics import cluster_size_calibration, bethe_bloch_calibration
from torchic.physics.ITS import average_cluster_size

from utils.particles import ParticleMasses

if __name__ == '__main__':

    particle = 'Pr'
    #particle = 'He'


    infile_path = '/Users/glucia/Projects/ALICE/data/lithium/MC/LHC25a4.root'
    folder_name = 'DF*'
    outfile_path = f'/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/MC/tpc_calibration_{particle}.root'

    if particle == 'Pr':
        dataset = Dataset.from_root(infile_path, 'O2he3hadtable', folder_name, columns=['fPtHad', 'fEtaHad', 'fSignalTPCHad', 'fItsClusterSizeHad'])
        dataset['AvgClusterSize'], __ = average_cluster_size(dataset['fItsClusterSizeHad'])
        dataset['AvgClSizeCosLam'] = dataset['AvgClusterSize'] / np.cosh(dataset['fEtaHad'])
        dataset['BetaGamma'] = dataset['fPtHad'] * np.cosh(dataset['fEtaHad']) / ParticleMasses['Pr']
    else:
        dataset = Dataset.from_root(infile_path, 'O2he3hadtable', folder_name, columns=['fPtHe3', 'fEtaHe3', 'fSignalTPCHe3', 'fItsClusterSizeHe3'])
        dataset['AvgClusterSize'], __ = average_cluster_size(dataset['fItsClusterSizeHe3'])
        dataset['AvgClSizeCosLam'] = dataset['AvgClusterSize'] / np.cosh(dataset['fEtaHe3'])
        dataset['BetaGamma'] = dataset['fPtHe3'] * np.cosh(dataset['fEtaHe3']) / ParticleMasses['He']

    outfile = TFile.Open(outfile_path, 'recreate')
    axis_spec_bg = AxisSpec(100, 0, 5, 'bg', ';;')
    axis_spec_clsize = AxisSpec(75, 0, 15, 'cluster_size_cal', ';#beta#gamma;#LT ITS Cluster Size #GT #times cos #LT #lambda #GT')
    axis_spec_tpcsignal = AxisSpec(200, 0, 400, 'cluster_size_cal', ';#beta#gamma;TPC signal') # pr
    #axis_spec_tpcsignal = AxisSpec(200, 0, 2000, 'tpc_size_cal', ';#beta#gamma;TPC signal') # he

    clsize_dir = outfile.mkdir('clsize_dir')
    h2_clsize = dataset.build_th2('BetaGamma', 'AvgClSizeCosLam', axis_spec_bg, axis_spec_clsize)
    clsize_dir.cd()
    h2_clsize.Write()

    # RooFit initialization
    clsize = RooRealVar('fITSClSizeCosLam', 'cluster size', 0., 15.)
    clsize_cb_mean = RooRealVar('clsize_cb_mean', 'cb_mean', 0., 15, 'GeV/c^{2}')
    clsize_cb_sigma = RooRealVar('clsize_cb_sigma', 'cb_sigma', 0., 2., 'GeV/c^{2}')
    clsize_cb_alpha = RooRealVar('clsize_cb_a', 'cb_a', -10., 0.)
    clsize_cb_n = RooRealVar('clsize_cb_n', 'cb_n', 0., 30.)
    clsize_signal = RooCrystalBall('clsize_cb', 'cb', clsize, clsize_cb_mean, clsize_cb_sigma, clsize_cb_alpha, clsize_cb_n, doubleSided=False)
    #b0 = RooRealVar('b0', 'b0', -1., 1.)
    #bkg = RooExponential('bkg', 'bkg', clsize, b0)
    #sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
    clsize_fit_params = {'clsize_cb_mean': clsize_cb_mean,
                         'clsize_cb_sigma': clsize_cb_sigma,
                         'clsize_cb_alpha': clsize_cb_alpha,
                         'clsize_cb_n': clsize_cb_n,
                         #'b0': b0
                        }

    clsize_first_bin_fit = h2_clsize.GetXaxis().FindBin(0.7)
    clsize_last_bin_fit = h2_clsize.GetXaxis().FindBin(3.5)    # pr
    #clsize_last_bin_fit = h2_clsize.GetXaxis().FindBin(2.5)     # he

    n_low_bg_bins, n_high_bg_bins = (0, 0)
    for ibin in range(clsize_first_bin_fit, clsize_last_bin_fit+1):
        if h2_clsize.GetXaxis().GetBinCenter(ibin) < 1.:    n_low_bg_bins += 1
        else:                                               n_high_bg_bins += 1
    
    fit_range_low = [[0., 15.] for ibin in range(n_low_bg_bins)]
    fit_range_high = [[0., 15.]for ibin in range(n_high_bg_bins)]
    fit_ranges = fit_range_low + fit_range_high

    clsize_kwargs = {'first_bin_fit_by_slices': clsize_first_bin_fit,
                     'last_bin_fit_by_slices': clsize_last_bin_fit,
                     'fit_ranges': fit_ranges,
                     'output_dir': clsize_dir,
                     'simil_bethe_bloch_pars': {'kp1': 1.18941,
                                                'kp2': 1.53792,
                                                'kp3': 1.69961
                                               },
                     'signal_func_name': 'clsize_cb',
                     #'signal_func_name': 'exp_mod_gaus_0',
                    }    

    clsize_pars, its_resolution = cluster_size_calibration(clsize, h2_clsize, clsize_dir, clsize_signal, clsize_fit_params, 
                                                           #fit_mc=True, 
                                                           **clsize_kwargs)

    input()
    ## TPC

    tpc_dir = outfile.mkdir('tpc_dir')
    if particle == 'Pr':
        h2_tpcsignal = dataset.build_th2('BetaGamma', 'fSignalTPCHad', axis_spec_bg, axis_spec_tpcsignal)
    else:
        h2_tpcsignal = dataset.build_th2('BetaGamma', 'fSignalTPCHe3', axis_spec_bg, axis_spec_tpcsignal)
    tpc_dir.cd()
    h2_tpcsignal.Write()

    # RooFit initialization
    tpcsignal = RooRealVar('fSignalTPC', 'tpc signal', 0., 500.)                   # pr     
    tpc_cb_mean = RooRealVar('tpc_cb_mean', 'cb_mean', 0., 500, 'GeV/c^{2}')       # pr
    tpc_cb_sigma = RooRealVar('tpc_cb_sigma', 'cb_sigma', 0., 50., 'GeV/c^{2}')    # pr
    #tpcsignal = RooRealVar('fSignalTPC', 'tpc signal', 0., 2000.)                   # he
    #tpc_cb_mean = RooRealVar('tpc_cb_mean', 'cb_mean', 0., 2000, 'GeV/c^{2}')       # he
    #tpc_cb_sigma = RooRealVar('tpc_cb_sigma', 'cb_sigma', 0., 200., 'GeV/c^{2}')    # he
    tpc_cb_alpha = RooRealVar('tpc_cb_a', 'cb_a', -10., 0.)
    tpc_cb_n = RooRealVar('tpc_cb_n', 'cb_n', 0., 30.)
    tpc_signal = RooCrystalBall('tpc_cb', 'cb', tpcsignal, tpc_cb_mean, tpc_cb_sigma, tpc_cb_alpha, tpc_cb_n, doubleSided=False)
    #b0 = RooRealVar('b0', 'b0', -1., 1.)
    #bkg = RooExponential('bkg', 'bkg', tpcsignal, b0)
    #sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
    tpc_fit_params = {'tpc_cb_mean': tpc_cb_mean,
                      'tpc_cb_sigma': tpc_cb_sigma,
                      'tpc_cb_alpha': tpc_cb_alpha,
                      'tpc_cb_n': tpc_cb_n,
                      #'b0': b0
                     }
    
    tpc_first_bin_fit = h2_tpcsignal.GetXaxis().FindBin(0.5)   # pr
    #tpc_first_bin_fit = h2_tpcsignal.GetXaxis().FindBin(0.8)    # he
    tpc_last_bin_fit = h2_tpcsignal.GetXaxis().FindBin(3.5)

    n_low_bg_bins, n_high_bg_bins = (0, 0)
    for ibin in range(tpc_first_bin_fit, tpc_last_bin_fit+1):
        if h2_tpcsignal.GetXaxis().GetBinCenter(ibin) < 1.5:    n_low_bg_bins += 1
        else:                                                   n_high_bg_bins += 1
    
    fit_range_low = [[0., 400.] for ibin in range(n_low_bg_bins)]  # pr
    fit_range_high = [[0., 100.]for ibin in range(n_high_bg_bins)] # pr
    #fit_range_low = [[0., 2000.] for ibin in range(n_low_bg_bins)]  # he
    #fit_range_high = [[0., 800.]for ibin in range(n_high_bg_bins)]  # he
    fit_ranges = fit_range_low + fit_range_high

    tpc_kwargs = {'first_bin_fit_by_slices': tpc_first_bin_fit,
                    'last_bin_fit_by_slices': tpc_last_bin_fit,
                    'fit_ranges': fit_ranges,
                    'output_dir': tpc_dir,
                    'simil_bethe_bloch_pars': {'kp1': 0.903,
                                                 'kp2': 2.014,
                                                 'kp3': 2.440
                                                },
                    'signal_func_name': 'tpc_cb',
                    #'signal_func_name': 'exp_mod_gaus_0',
                     }
    
    tpc_pars = bethe_bloch_calibration(tpcsignal, h2_tpcsignal, tpc_dir, tpc_signal, tpc_fit_params, **tpc_kwargs)

    outfile.Close()
    