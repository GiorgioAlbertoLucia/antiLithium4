'''
    Code to run calibration of ITS and TPC parametrisations
'''

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from ROOT import TFile, TCanvas, TF1
from ROOT import RooRealVar, RooCrystalBall, RooGExpModel, RooAddPdf, RooGaussian
from torchic import Dataset, AxisSpec

from torchic.physics.ITS import expected_cluster_size, sigma_its, average_cluster_size
from torchic.physics import BetheBloch, py_BetheBloch

from utils import calibration_fit_slice, create_graph, initialize_means_and_sigmas
import sys
sys.path.append('..')
from analysis.utils.parametrisations import ITS as ITS_params
from analysis.utils.particles import ParticleMasses

DATASET_COLUMN_NAMES = {
    'P': 'fInnerParamTPC',
    'Pt': 'fPt',
    'Eta': 'fEta',
    'TpcSignal': 'fSignalTPC',
    'Chi2TPC': 'fChi2TPC',
    'ItsClusterSize': 'fItsClusterSize',
}

def TPC_calibration_He(infile_path:str, folder_name:str, tree_name:str, outfile:TFile, column_names:dict=DATASET_COLUMN_NAMES):

    dataset = Dataset.from_root(infile_path, tree_name, folder_name, columns=[col for col in column_names.values()])

    dataset.query(f'0.5 < {column_names["Chi2TPC"]} < 4', inplace=True)
    dataset.eval(f'fSign = {column_names["Pt"]} / abs({column_names["Pt"]})', inplace=True)
    dataset.eval(f'fBetaGamma = fSign * {column_names["P"]} * 2 / {ParticleMasses["He"]}', inplace=True)

    axis_spec_betagamma = AxisSpec(160, -8, 8, 'beta_gamma', ';#beta#gamma;dE/dx (a.u.)')
    axis_spec_tpcsignal = AxisSpec(100, 0, 1200, 'tpc_signal', ';#beta#gamma;dE/dx (a.u.)')
    h2_tpc = dataset.build_th2('fBetaGamma', column_names['TpcSignal'], axis_spec_betagamma, axis_spec_tpcsignal)

    bg_min = 0.6
    bg_max = 4.

    tpc_dir = outfile.mkdir('tpc_dir')
    
    tpc_signal = RooRealVar('fSignalTPC', 'dE/dx (a.u.)', 0., 1200.)
    signal_pars = {
        'mean': RooRealVar('mean', 'mean', 800., 0., 1200., 'GeV/c^{2}'),
        'sigma': RooRealVar('sigma', 'sigma', 0.01, 1000, 'GeV/c^{2}'),
    }
    signal = RooGaussian('signal', 'signal', tpc_signal, signal_pars['mean'], signal_pars['sigma'])
    bkg_pars = {
        'mean_1': RooRealVar('mean_1', 'mean_1', 0., 0., 1200., 'GeV/c^{2}'),
        'sigma_1': RooRealVar('sigma_1', 'sigma_1', 60, 20., 100, 'GeV/c^{2}'),
        'rlife_1': RooRealVar('rlife_1', 'rlife_1', 0., 10.),
    }
    bkg_1 = RooGaussian('bkg_1', 'bkg_1', tpc_signal, bkg_pars['mean_1'], bkg_pars['sigma_1']) #, bkg_pars['rlife_1'])

    fit_params = None
    g_mean = None

    for sign in ['matter', 'antimatter']:

        if sign == 'matter':
            slice_range = [bg_min, bg_max]
        else:
            slice_range = [-bg_max, -bg_min]

        model = None
        fit_results_df = None

        bg_bin_min = h2_tpc.GetXaxis().FindBin(slice_range[0])
        bg_bin_max = h2_tpc.GetXaxis().FindBin(slice_range[1])
        for bg_bin in range(bg_bin_min, bg_bin_max+1):
            
            tpc_signal.setRange(0., 1200.)
            bg = h2_tpc.GetXaxis().GetBinCenter(bg_bin)
            bg_low_edge = h2_tpc.GetXaxis().GetBinLowEdge(bg_bin)
            bg_high_edge = h2_tpc.GetXaxis().GetBinLowEdge(bg_bin+1)
            
            h_tpc = h2_tpc.ProjectionY(f'tpc_signal_{bg:.2f}', bg_bin, bg_bin, 'e')
            if np.abs(bg) < 0.8:
                sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
                model = RooAddPdf('model', 'model', [signal, bkg_1], [sig_frac])
                
                means, sigmas = initialize_means_and_sigmas(h_tpc, 2)
                signal_pars['mean'].setVal(means[1])
                signal_pars['sigma'].setVal(sigmas[1])
                bkg_pars['mean_1'].setVal(means[0])
                bkg_pars['sigma_1'].setVal(sigmas[0])
            else:
                model = signal

            frame, fit_results = calibration_fit_slice(model, h_tpc, tpc_signal, signal_pars, bg_low_edge, bg_high_edge)
            fit_results['bg'] = np.abs(bg)
            if fit_results_df is None:
                fit_results_df = pd.DataFrame.from_dict([fit_results])
            else:
                fit_results_df = pd.concat([fit_results_df, pd.DataFrame.from_dict([fit_results])], ignore_index=True)

            canvas = TCanvas(f'cNSigmaTPC_{bg:.2f}', f'cNSigmaTPC_{bg:.2f}', 800, 600)
            frame.Draw()
            tpc_dir.cd()
            canvas.Write()

        g_mean = create_graph(fit_results_df, 'bg', 'mean', 0, 'mean_err', 
                                    f'g_mean_{sign}', ';#beta#gamma;#LT m_{TOF} #GT (GeV/#it{c}^{2})')
        f_mean = TF1('f_mean', BetheBloch, bg_min, bg_max, 5)
        f_mean.SetParameters(-241.4902, 0.374245, 1.397847, 1.0782504, 2.048336)
        g_mean.Fit(f_mean, 'RMS+')
        fit_params = [f_mean.GetParameter(i) for i in range(5)]
        g_resolution = create_graph(fit_results_df, 'bg', 'resolution', 0, 'resolution_err', 
                                    f'g_resolution_{sign}', ';#beta#gamma;#sigma_{m_{TOF}} / #LT m_{TOF} #GT')
        #f_resolution = TF1('f_resolution', '[0] + [1]/x^2', bg_min, bg_max)
        f_resolution = TF1('f_resolution', '[0]', bg_min, bg_max)
        f_resolution.SetParameters(0.07)
        g_resolution.Fit(f_resolution, 'RMS+')
        
        tpc_dir.cd()
        g_mean.Write()
        g_resolution.Write()

    np_bethe_bloch = np.vectorize(py_BetheBloch)
    print(f'{fit_params=}')
    dataset['fExpTpcSignal'] = np_bethe_bloch(np.abs(dataset['fBetaGamma']), *fit_params)
    dataset['fNSigmaTPC'] = (dataset[column_names['TpcSignal']] - dataset['fExpTpcSignal']) / (dataset['fExpTpcSignal'] * 0.09)

    axis_spec_nsigmatpc = AxisSpec(100, -5, 5, 'nsigma_tpc', ';#beta#gamma;#LT m_{TOF} #GT (GeV/#it{c}^{2})')
    h2_nsigmatpc = dataset.build_th2('fBetaGamma', 'fNSigmaTPC', axis_spec_betagamma, axis_spec_nsigmatpc)
    h2_exptpc = dataset.build_th2('fBetaGamma', 'fExpTpcSignal', axis_spec_betagamma, axis_spec_tpcsignal)

    tpc_dir.cd()
    h2_tpc.Write()
    h2_nsigmatpc.Write()
    h2_exptpc.Write('exp_tpc_signal')
    canvas = TCanvas('cNSigmaTPC', 'cNSigmaTPC', 800, 600)
    h2_tpc.Draw('colz')
    g_mean.Draw('same')
    canvas.Write()


if __name__ == '__main__':

    
    infile_path = '/data/galucia/lithium_local/same/LHC23_PbPb_pass4_long_same_lsus.root'
    #infile_path = '/data/galucia/lithium_local/same/LHC24as_pass1_same.root'
    folder_name = 'DF*'
    tree_name = 'O2he3hadtable'
    outfile = TFile('/home/galucia/antiLithium4/calibration/output/TPC_calibration_23.root', 'recreate')

    columns_names = {key: value+'He3' for key, value in DATASET_COLUMN_NAMES.items()}
    TPC_calibration_He(infile_path, folder_name, tree_name, outfile, columns_names)

    outfile.Close()
    