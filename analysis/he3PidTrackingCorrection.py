from torchic import Dataset, AxisSpec, Roofitter
from ROOT import TFile, RooRealVar, TDirectory, TGraphErrors, TF1
import pandas as pd
import numpy as np

def parametrise_resolution(h2_res, xmin, xmax, ymin, ymax, output_file: TDirectory):

    fits_dir = output_file.mkdir('fits')
    fit_results = pd.DataFrame()

    x = RooRealVar('x', 'x', ymin, ymax)
    #roofitter = Roofitter(x, ['exp_mod_gaus', 'crystal_ball'])
    roofitter = Roofitter(x, ['exp_mod_gaus'])
    roofitter.init_param('exp_mod_gaus_0_mean', 0.0, -0.5, 0.5)
    roofitter.init_param('exp_mod_gaus_0_sigma', 0.01, 0.0, 1.0)
    roofitter.init_param('exp_mod_gaus_0_tau', 0.1, -10.0, 10.0)
    
    low_bin = h2_res.GetXaxis().FindBin(xmin)
    high_bin = h2_res.GetXaxis().FindBin(xmax)

    for ibin in range(low_bin, high_bin+1):
        h1 = h2_res.ProjectionY(f'h1_{ibin}', ibin, ibin)
        fractions = roofitter.fit(h1, xmin, xmax)
        roofitter.plot(fits_dir, canvas_name=f'c_{h1.GetName()}', funcs_to_plot=list(roofitter.pdfs.keys()))
        
        bin_fit_results = pd.DataFrame.from_dict({key: [val] for key, val in roofitter.fit_results.items()})
        bin_fit_results['bin_center'] = h2_res.GetXaxis().GetBinCenter(ibin)
        bin_fit_results['unnorm_integral'] = h1.Integral() * bin_fit_results['integral'] * fractions[0].getVal()

        fit_results = pd.concat([fit_results, bin_fit_results], ignore_index=True)
    
    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'].values[0])/2.
    fit_results['bin_error'] = bin_error

    fit_results['mean_error'] = fit_results['exp_mod_gaus_0_sigma'] / fit_results['unnorm_integral']
    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center'], dtype=np.float32), np.array(fit_results['exp_mod_gaus_0_mean'], dtype=np.float32), np.array(fit_results['bin_error'], dtype=np.float32), np.array(fit_results['mean_error'], dtype=np.float32))
    graph_mean.SetName('gResolutionMean')
    graph_mean.SetTitle(';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}')

    xmin_fit = 1.6
    xmax_fit = 2.4
    # fit_res = TF1('fit_res', '[0] / (1 + exp((x - [1])/[2]) )', xmin_fit, xmax_fit)
    # fit_res.SetParameter(0, 0.1)
    # fit_res.SetParameter(1, 1.9)
    # fit_res.SetParameter(2, 0.1)
    # fit_res = TF1('fit_res', '[0] + [1] * x + [2] * x^2', xmin_fit, xmax_fit)
    fit_res = TF1('fit_res', '[0] + [1] * x', xmin_fit, xmax_fit)
    fit_res.SetParameter(0, 0.1)
    fit_res.SetParameter(1, -1.)
    # fit_res.SetParameter(2, 0.1)
    graph_mean.Fit(fit_res, 'RMS+')

    # fit_res_results = [fit_res.GetParameter(i) for i in range(3)]
    fit_res_results = [fit_res.GetParameter(i) for i in range(2)]

    output_file.cd()
    graph_mean.Write()
    fit_res.Write()

    return fit_res_results

if __name__ == '__main__':

    dataset = Dataset.from_root('/data/galucia/lithium_local/MC/LHC24b2c_nucleispectra_injectedhe3.root', 
                                tree_name='O2nucleitable',
                                folder_name='DF*')

    print(f'{dataset.columns=}')

    # select he3
    dataset.query('abs(fPDGcode) == 1000020030', inplace=True)
    readPidTracking = lambda x: (x >> 12) & 0x1F
    readPidTracking = np.vectorize(readPidTracking)
    dataset['fPidTracking'] = readPidTracking(dataset['fFlags'])
    dataset.add_subset('H3', dataset['fPidTracking'] == 6)
    dataset.add_subset('He4', dataset['fPidTracking'] == 8)
    dataset.loc[dataset['fPDGcode'] < 0, 'fPt'] = - dataset['fPt']

    dataset['fPt'] = dataset['fPt'] * 2.0
    dataset['fPtDiff'] = dataset['fPt'] - dataset['fgPt']
    dataset['fPtRes'] = dataset['fPtDiff'] / dataset['fPt']

    output_file_path = '/home/galucia/antiLithium4/analysis/output/MC/he3PidTrackingCorrection.root'
    axis_spec_pt = AxisSpec(100, 0., 5., 'fPt', ';#it{p}_{T}^{rec} (GeV/#it{c})')
    axis_spec_gpt = AxisSpec(100, 0., 5., 'fgPt', ';#it{p}_{T}^{gen} (GeV/#it{c})')

    h_pt = dataset.build_hist('fPt', axis_spec_pt)
    h_gpt = dataset.build_hist('fgPt', axis_spec_gpt)
    h2_res = dataset.build_hist('fPt', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtRes', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'))
    h2_resh3 = dataset.build_hist('fPt', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResH3', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'), subset='H3')
    h2_reshe4 = dataset.build_hist('fPt', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResHe4', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'), subset='He4')

    output_file = TFile.Open(output_file_path, 'recreate')
    h_pt.Write()
    h_gpt.Write()
    h2_res.Write()
    h2_resh3.Write()
    h2_reshe4.Write()

    fit_res_results = parametrise_resolution(h2_resh3, 1.4, 2.35, -0.5, 0.5, output_file)

    # pol1 correction
    dataset.loc[dataset['fPidTracking'] == 6, 'fPt'] = dataset['fPt'] - dataset['fPt']*(fit_res_results[0] + fit_res_results[1] * dataset['fPt'])
    # pol2 correction
    # dataset.loc[dataset['fPidTracking'] == 6, 'fPt'] = dataset['fPt'] - dataset['fPt']*(fit_res_results[0] + fit_res_results[1] * dataset['fPt'] + fit_res_results[2] * dataset['fPt']**2)
    dataset['fPtDiff'] = dataset['fPt'] - dataset['fgPt']
    dataset['fPtRes'] = dataset['fPtDiff'] / dataset['fPt']
    h2_res_corrected = dataset.build_hist('fPt', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResCorrected', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'))
    h2_resh3_corrected = dataset.build_hist('fPt', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResH3Corrected', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'), subset='H3')

    output_file.cd()
    h2_res_corrected.Write()
    h2_resh3_corrected.Write()
