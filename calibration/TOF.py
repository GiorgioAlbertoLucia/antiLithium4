'''
    Code to run calibration of ITS and TPC parametrisations
'''

import numpy as np
import pandas as pd
from ROOT import TFile, TCanvas, TF1
from ROOT import RooRealVar, RooCrystalBall, RooGExpModel, RooAddPdf, RooGaussian
from torchic import Dataset, AxisSpec

from torchic.physics.ITS import expected_cluster_size, sigma_its, average_cluster_size

from utils import calibration_fit_slice, create_graph
import sys
sys.path.append('..')
from analysis.utils.parametrisations import ITS as ITS_params
from analysis.utils.particles import ParticleMasses

DATASET_COLUMN_NAMES = {
    'Pt': 'fPt',
    'Eta': 'fEta',
    'TofMass': 'fMassTOF',
    'Chi2TPC': 'fChi2TPC',
    'NSigmaTPC': 'fNSigmaTPC',
    'ItsClusterSize': 'fItsClusterSize',
}

def TOF_calibration_Pr(infile_path:str, folder_name:str, tree_name:str, outfile:TFile, column_names:dict=DATASET_COLUMN_NAMES):

    dataset = Dataset.from_root(infile_path, tree_name, folder_name, columns=[col for col in column_names.values()])
    # proton selections
    dataset.query(f'0.5 < {column_names["Chi2TPC"]} < 4', inplace=True)
    #dataset.query(f'abs({column_names["NSigmaTPC"]}) < 2', inplace=True)
    #dataset.eval(f'fBetaGamma = abs({column_names["Pt"]})*cosh({column_names["Eta"]})/{ParticleMasses["Pr"]}', inplace=True)
    #its_params = list(ITS_params.pr_exp_params.values()) + list(ITS_params.pr_res_params.values())
    #dataset['fClSizeITSCosLam'], _ = average_cluster_size(dataset[column_names['ItsClusterSize']])
    #dataset['fExpClSizeITS'] = expected_cluster_size(dataset['fBetaGamma'], its_params)
    #dataset['fSigmaClSizeCosL'] = sigma_its(dataset['fBetaGamma'], its_params)
    #dataset.eval('fNSigmaITS = (fClSizeITSCosLam - fExpClSizeITS) / fSigmaClSizeCosL', inplace=True)
    #dataset.query('fNSigmaITS > -1.5', inplace=True)
    #dataset.drop(columns=['fExpClSizeITS', 'fSigmaClSizeCosL', 'fBetaGamma', 'fClSizeITSCosLam', 'fNSigmaITS', column_names['Eta'],
    #                      column_names['ItsClusterSize'], column_names['Chi2TPC'], column_names['NSigmaTPC']], inplace=True)

    axis_spec_pt = AxisSpec(100, -5, 5, 'pt', ';;')
    axis_spec_tofmass = AxisSpec(50, 0, 2, 'tof_mass', ';#it{p}_{T} (GeV/#it{c});m_{TOF} (GeV/#it{c}^{2})')
    h2_tof = dataset.build_th2(column_names['Pt'], column_names['TofMass'], axis_spec_pt, axis_spec_tofmass)

    pt_min = 0.7
    pt_max = 2

    tof_dir = outfile.mkdir('tof_dir')
    
    # RooFit initialization
    tof_mass = RooRealVar('fMassTOF', 'm_{TOF} (GeV/#it{c}^{2})', 0., 2.)
    signal_pars = {
        'mean': RooRealVar('mean', 'mean', 0.7, 1.2, 'GeV/c^{2}'),
        'sigma': RooRealVar('sigma', 'sigma', 0.01, 0.1, 'GeV/c^{2}'),
        'aL': RooRealVar('aL', 'aL', 0.7, 30.),
        'nL': RooRealVar('nL', 'nL', 0.3, 30.),
        'aR': RooRealVar('aR', 'aR', 0.7, 30.),
        'nR': RooRealVar('nR', 'nR', 0.3, 30.),
    }
    signal = RooCrystalBall('signal', 'signal', tof_mass, signal_pars['mean'], signal_pars['sigma'],
                            signal_pars['aL'], signal_pars['nL'], signal_pars['aR'], signal_pars['nR'])
    bkg_pars = {
        'mean_1': RooRealVar('mean_1', 'mean_1', 0., 0.3, 'GeV/c^{2}'),
        'sigma_1': RooRealVar('sigma_1', 'sigma_1', 0.01, 0.1, 'GeV/c^{2}'),
        'rlife_1': RooRealVar('rlife_1', 'rlife_1', 0., 10.),
        'mean_2': RooRealVar('mean_2', 'mean_2', 0.4, 0.6, 'GeV/c^{2}'),
        'sigma_2': RooRealVar('sigma_2', 'sigma_2', 0.01, 0.1, 'GeV/c^{2}'),
        'rlife_2': RooRealVar('rlife_2', 'rlife_2', 0., 10.),
        #'mean_3': RooRealVar('mean_3', 'mean_3', 0., 0.1, 'GeV/c^{2}'),
        #'sigma_3': RooRealVar('sigma_3', 'sigma_3', 0.01, 0.04, 'GeV/c^{2}'),
        #'rlife_3': RooRealVar('rlife_3', 'rlife_3', 0., 10.),
    }
    bkg_1 = RooGaussian('bkg_1', 'bkg_1', tof_mass, bkg_pars['mean_1'], bkg_pars['sigma_1']) #, bkg_pars['rlife_1'])
    bkg_2 = RooGaussian('bkg_2', 'bkg_2', tof_mass, bkg_pars['mean_2'], bkg_pars['sigma_2']) #, bkg_pars['rlife_2'])
    #bkg_3 = RooGaussian('bkg_3', 'bkg_3', tof_mass, bkg_pars['mean_3'], bkg_pars['sigma_3']) #, bkg_pars['rlife_3'])


    for sign in ['matter', 'antimatter']:

        if sign == 'matter':
            slice_range = [pt_min, pt_max]
        else:
            slice_range = [-pt_max, -pt_min]

        model = None
        fit_results_df = None

        pt_bin_min = h2_tof.GetXaxis().FindBin(slice_range[0])
        pt_bin_max = h2_tof.GetXaxis().FindBin(slice_range[1])
        for pt_bin in range(pt_bin_min, pt_bin_max+1):
            
            tof_mass.setRange(0.1, 2.)
            pt = h2_tof.GetXaxis().GetBinCenter(pt_bin)
            pt_low_edge = h2_tof.GetXaxis().GetBinLowEdge(pt_bin)
            pt_high_edge = h2_tof.GetXaxis().GetBinLowEdge(pt_bin+1)
            if np.abs(pt) < 0.9:
                sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
                bkg1_frac = RooRealVar('bkg1_frac', 'bkg1_frac', 0.5, 0., 1.)
                #bkg2_frac = RooRealVar('bkg2_frac', 'bkg2_frac', 0.5, 0., 1.)
                model = RooAddPdf('model', 'model', [signal, bkg_1, bkg_2],#, bkg_3], 
                                  [sig_frac, bkg1_frac]) #, bkg2_frac])
            else:
                model = signal

            h_tof = h2_tof.ProjectionY(f'tof_mass_{pt:.2f}', pt_bin, pt_bin, 'e')
            frame, fit_results = calibration_fit_slice(model, h_tof, tof_mass, signal_pars, pt_low_edge, pt_high_edge)
            fit_results['pt'] = np.abs(pt)
            if fit_results_df is None:
                fit_results_df = pd.DataFrame.from_dict([fit_results])
            else:
                fit_results_df = pd.concat([fit_results_df, pd.DataFrame.from_dict([fit_results])], ignore_index=True)

            canvas = TCanvas(f'cNSigmaTPC_{pt:.2f}', f'cNSigmaTPC_{pt:.2f}', 800, 600)
            frame.Draw()
            tof_dir.cd()
            canvas.Write()

        g_mean = create_graph(fit_results_df, 'pt', 'mean', 0, 'mean_err', 
                                    f'g_mean_{sign}', ';#it{p}_{T} (GeV/#it{c};#LT m_{TOF} #GT (GeV/#it{c}^{2})')
        f_mean = TF1('f_mean', '[0] + [1]*x', pt_min, pt_max)
        f_mean.SetParameters(1, -0.01)
        g_mean.Fit(f_mean, 'RMS+')
        g_resolution = create_graph(fit_results_df, 'pt', 'resolution', 0, 'resolution_err', 
                                    f'g_resolution_{sign}', ';#it{p}_{T} (GeV/#it{c};#sigma_{m_{TOF}} / #LT m_{TOF} #GT')
        f_resolution = TF1('f_resolution', '[0] + [1]*x', pt_min, pt_max)
        f_resolution.SetParameters(0., 0.04)
        g_resolution.Fit(f_resolution, 'RMS+')
        
        tof_dir.cd()
        g_mean.Write()
        g_resolution.Write()

    tof_dir.cd()
    h2_tof.Write()


if __name__ == '__main__':

    
    infile_path = '/data/galucia/lithium_local/same/LHC23_PbPb_pass4_long_same_lsus.root'
    folder_name = 'DF*'
    tree_name = 'O2he3hadtable'
    outfile = TFile('/home/galucia/antiLithium4/calibration/output/TOF_calibration.root', 'recreate')

    columns_names = {key: value+'Had' for key, value in DATASET_COLUMN_NAMES.items()}
    TOF_calibration_Pr(infile_path, folder_name, tree_name, outfile, columns_names)

    outfile.Close()
    