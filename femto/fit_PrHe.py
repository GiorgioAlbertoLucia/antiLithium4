'''
    Macro to fit the p-He3 correlation function
'''

import numpy as np
from torchic.core.histogram import HistLoadInfo, load_hist
from torchic.utils.terminal_colors import TerminalColors as tc
from ROOT import RooRealVar, RooDataHist, RooHistPdf, RooCrystalBall, RooFit, RooAddPdf, RooGenericPdf, RooProdPdf, RooStats, RooWorkspace
from ROOT import TFile, TH1F, TPaveText, TCanvas, TDirectory, TGraph

def prepare_coulomb_CF(hist: TH1F, outfile: TFile) -> RooHistPdf:

    kstar = RooRealVar('kstar', 'kstar', 0, 0.8)
    data_hist = RooDataHist('coulomb', 'coulomb', [kstar], Import=hist)
    pdf = RooHistPdf('pdf', 'pdf', {kstar}, data_hist, 0)

    frame = kstar.frame()
    data_hist.plotOn(frame)
    pdf.plotOn(frame)

    outfile.cd()
    frame.Write('coulomb_CF')

    return pdf

def prepare_peak_CF(hist: TH1F, me_hist: TH1F, outfile: TFile, use_mixed:bool=False) -> dict:

    if use_mixed:
        hist.Add(me_hist)
        tmp_hist = hist.Clone()
        for ibin in range(1, hist.GetNbinsX() + 1):
            value = hist.GetBinContent(ibin) / me_hist.GetBinContent(ibin) - 1
            error = value * np.sqrt(hist.GetBinError(ibin)**2 / hist.GetBinContent(ibin)**2 + me_hist.GetBinError(ibin)**2 / me_hist.GetBinContent(ibin)**2)
            tmp_hist.SetBinContent(ibin, value)
            tmp_hist.SetBinError(ibin, error)
        del hist
        hist = tmp_hist

    kstar = RooRealVar('kstar', 'kstar', 0., 0.25)
    fit_params = {
        'cb_mean': RooRealVar('cb_mean', 'cb_mean', 0., 0.25, 'GeV/c^{2}'),
        'cb_sigma': RooRealVar('cb_sigma', 'cb_sigma', 0., 0.1, 'GeV/c^{2}'),
        'cb_aL': RooRealVar('cb_aL', 'cb_aL', 0.1, 10.),
        'cb_nL': RooRealVar('cb_nL', 'cb_nL', 0.1, 30.),
        'cb_aR': RooRealVar('cb_aR', 'cb_aR', 0.1, 10.),
        'cb_nR': RooRealVar('cb_nR', 'cb_nR', 0.1, 30.),
    }
    signal = RooCrystalBall('cb', 'cb', kstar, fit_params['cb_mean'], fit_params['cb_sigma'], fit_params['cb_aL'], 
                            fit_params['cb_nL'], fit_params['cb_aR'], fit_params['cb_nR'])

    datahist = RooDataHist('datahist', 'datahist', [kstar], hist)
    sumw2error_bool = not use_mixed
    signal.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.2), SumW2Error=sumw2error_bool)
    frame = kstar.frame()
    datahist.plotOn(frame)
    signal.plotOn(frame)

    text = TPaveText(0.7, 0.7, 0.9, 0.9, 'NDC')
    for name, param in fit_params.items():
        text.AddText(f'{name} = {param.getVal():.4f} +/- {param.getError():.4f}')
    frame.addObject(text)
    frame.Draw()
    outfile.cd()
    frame.Write(f'li4_peak')

    return fit_params

def fit_CF(hist: TH1F, outfile, li4_peak_params: dict, sign: str, cent: str = ''):

    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    print(tc.RED+f'Fit the CF'+tc.RESET)

    outfile.cd()
    hist.Write('data_hist')

    kstar = RooRealVar('kstar', '#it{k}* (GeV/#it{c})', 0., 0.8)
    
    infile_path_coulomb = f'/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent{cent}_new.root'
    if cent == '':
        infile_path_coulomb = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_new.root'
    coulomb_template_hist = load_hist(HistLoadInfo(infile_path_coulomb, 'hHe3_p_Coul_CF_LS'))
    data_hist_coulomb = RooDataHist('coulomb', 'coulomb', [kstar], Import=coulomb_template_hist)
    coulomb_pdf = RooHistPdf('coulomb_pdf', 'coulomb_pdf', {kstar}, data_hist_coulomb, 0)

    data_hist = RooDataHist('CF', 'CF', [kstar], Import=hist)
    
    par_vals = {name: param.getVal() for name, param in li4_peak_params.items()}
    mean_val = li4_peak_params['cb_mean'].getVal()
    mean_err = li4_peak_params['cb_mean'].getError()
    fit_params = {
                    'li4_mean': RooRealVar('li4_mean', '#mu', par_vals['cb_mean'], mean_val - 3 * mean_err, mean_val + 3 * mean_err, 'GeV/c^{2}'),
                    'li4_sigma': RooRealVar('li4_sigma', '#sigma', par_vals['cb_sigma'], par_vals['cb_sigma'], 4 * par_vals['cb_sigma'], 'GeV/c^{2}'),
                    'li4_cb_aL': RooRealVar('li4_cb_aL', 'a_{L}', par_vals['cb_aL'], 0.1, 10.),
                    'li4_cb_nL': RooRealVar('li4_cb_nL', 'n_{L}', par_vals['cb_nL'], 0.1, 30.),
                    'li4_cb_aR': RooRealVar('li4_cb_aR', 'a_{R}', par_vals['cb_aR'], 0.1, 10.),
                    'li4_cb_nR': RooRealVar('li4_cb_nR', 'n_{R}', par_vals['cb_nR'], 0.1, 30.),
                }

    #fit_params['li4_mean'].setConstant(True)
    fit_params['li4_cb_aL'].setConstant(True)
    fit_params['li4_cb_nL'].setConstant(True)
    fit_params['li4_cb_aR'].setConstant(True)
    fit_params['li4_cb_nR'].setConstant(True)

    li4_signal = RooCrystalBall('li4_signal', 'li4_signal', kstar, fit_params['li4_mean'], fit_params['li4_sigma'], fit_params['li4_cb_aL'], 
                                fit_params['li4_cb_nL'], fit_params['li4_cb_aR'], fit_params['li4_cb_nR'])

    li4_counts = RooRealVar('li4_counts', 'li4_counts', 0.5, 0., 1e3)
    coulomb_counts = RooRealVar('coulomb_counts', 'coulomb_counts', 0.9, 0., 1e3)
    model = RooAddPdf('model', 'model', [li4_signal, coulomb_pdf], [li4_counts, coulomb_counts])
    
    kstar.setRange('prefit_range', 0.3, 0.8)
    li4_counts.setVal(0)
    li4_counts.setConstant(True)
    model.fitTo(data_hist, Range='prefit_range', PrintLevel=-1, Extended=True, SumW2Error=False)
    coulomb_counts.setConstant(True)
    li4_counts.setConstant(False)
    #model = RooAddPdf('model', 'model', [li4_signal], [li4_counts])
    model.fitTo(data_hist, PrintLevel=-1, Extended=True, SumW2Error=False)

    frame = kstar.frame()
    frame.SetYTitle('C(#it{k}*)')
    data_hist.plotOn(frame)
    model.plotOn(frame, LineStyle=1, LineColor=2)
    model.plotOn(frame, Components={li4_signal}, LineStyle=2, LineColor=4)
    model.plotOn(frame, Components={coulomb_pdf}, LineStyle=2, LineColor=3)

    # fit significance
    ws = RooWorkspace('workspace')
    model_config = RooStats.ModelConfig('model_config', ws)
    model_config.SetPdf(model)
    model_config.SetParametersOfInterest({li4_counts})
    model_config.SetObservables({kstar})
    getattr(ws, 'import')(model_config)
    getattr(ws, 'import')(data_hist)
    outfile.cd()
    ws.Write('workspace')

    kstar_max = 0.2
    kstar.setRange('int_range', 0., kstar_max)
    signal_integral = li4_signal.createIntegral(kstar, kstar, 'int_range').getVal() * li4_counts.getVal()
    model_integral = model.createIntegral(kstar, kstar, 'int_range').getVal() * (li4_counts.getVal() + coulomb_counts.getVal())
    significance = signal_integral / np.sqrt(model_integral)
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    print(tc.RED+f'Fit significance'+tc.RESET)
    print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
    print(f'S = S_R / sqrt(S_R + B_R)')
    print(f'S_R = {signal_integral:.3f}')
    print(f'B_R = {model_integral - signal_integral:.3f}')
    print(f'Significance = {significance:.3f}')
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)

    # compute li4 yield
    kstar_li4 = fit_params['li4_mean'].getVal()
    kstar_max = kstar_li4 + 5 * fit_params['li4_sigma'].getVal()
    kstar.setRange('int_range', 0., kstar_max)
    kstar2_func = RooGenericPdf('kstar2_func', 'kstar2_func', 'kstar * kstar', [kstar])
    S_R_func = RooProdPdf('S_R_func', 'S_R_func', [li4_signal, kstar2_func])
    #S_R_func.plotOn(frame, LineStyle=2, LineColor=6)
    S_R = S_R_func.createIntegral(kstar, kstar, 'int_range').getVal()       # normalized to 1
    S_R *= li4_counts.getVal()                                            # normalised to fraction of the signal
    S_R *= hist.GetEntries()                                                # normalised to the number of entries in the histogram
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    print(tc.RED+f'Li4 yield'+tc.RESET)
    print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
    print(f'dN(Li4)/(d3k* d3P) = 3/4 S_R(k*max) dN(p)/d3p dN(He3)/d3p')
    print(f'S_R(k*max) is the integral of the Li4 signal in the CF in the range [0, k*max]')
    print(f'k*max = k*(Li4) + 3 * sigma')
    print(f'\nk*(Li4) = {kstar_li4:.4f}')
    print(f'sigma = {fit_params["li4_sigma"].getVal():.4f}')
    print(f'k*max = {kstar_max:.4f}')
    print(f'\nS_R(k*max) = {S_R:.2f}')
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)

    text = TPaveText(0.55, 0.25, 0.75, 0.65, 'NDC')
    for name, param in fit_params.items():
        tmp_name = name.replace('li4_', '')
        tmp_name = tmp_name.replace('cb_', '')
        if (tmp_name != 'mean') and (tmp_name != 'sigma'):
            text.AddText(f'{param.GetTitle()} = {param.getVal():.4f} (fixed)')
        else:
            text.AddText(f'{param.GetTitle()} = {param.getVal():.4f} #pm {param.getError():.4f}')
    #text.AddText(f'li4_counts = {li4_counts.getVal():.4f} #pm {li4_counts.getError():.4f}')
    #text.AddText(f'coulomb_counts = {coulomb_counts.getVal():.4f} #pm {coulomb_counts.getError():.4f}')
    text.AddText(f'#chi^{{2}} / NDF = {frame.chiSquare():.4f}')
    text.SetBorderSize(0)
    text.SetFillStyle(0)
    #text.AddText(f'S_R = {S_R:.2f}')
    frame.addObject(text)

    hpull = frame.pullHist()
    frame_pull = kstar.frame()
    frame_pull.addPlotable(hpull, 'P')
    frame_pull.SetYTitle('Pull')
    
    outfile.cd()
    frame.Write('fit_CF')
    frame_pull.Write('pull_CF')

    canvas = TCanvas('canvas', 'canvas', 800, 800)
    frame.Draw()
    canvas.SaveAs(f'/home/galucia/antiLithium4/femto/output/fit_PrHe_{sign}_{cent}.pdf')

def compute_significance(outdir: TDirectory, sign: str, cent: str = ''):
    
    ws = outdir.Get(f'workspace')

    CF_datahist = ws.data('CF')
    model = ws.obj('model_config')
    poi = model.GetParametersOfInterest().first()
    model.SetSnapshot({poi})
    # create the b-only model
    bkg_model = model.Clone()
    bkg_model.SetName('coulomb_pdf_config')
    
    old_poi_val = poi.getVal()
    poi.setVal(0)
    bkg_model.SetSnapshot(poi)
    poi.setVal(old_poi_val)

    bkg_model.Print()
    ws.var('li4_mean').setConstant(True)
    ws.var('li4_sigma').setConstant(True)

    asymp_calc = RooStats.AsymptoticCalculator(CF_datahist, model, bkg_model)
    asymp_calc.SetPrintLevel(0)
    asymp_calc_result = asymp_calc.GetHypoTest()
    p_value = asymp_calc_result.NullPValue()
    p_value_err = asymp_calc_result.NullPValueError()
    significance = asymp_calc_result.Significance()
    significance_error = asymp_calc_result.SignificanceError()
    
    ### perform a scan in kstar and compute the significance
    '''
    kstars = []
    p0_values = []
    p0_values_expected = []
    kstar_array = np.linspace(ws.var('kstar').getMin(), ws.var('kstar').getMax(), 40)
    for kstar in kstar_array:

        ws.var('kstar').setVal(kstar)
        ws.var('kstar').setConstant(True)
        asymp_calc_scan = RooStats.AsymptoticCalculator(CF_datahist, model, bkg_model)
        asymp_calc_scan.SetOneSidedDiscovery(True)
        asym_calc_result_scan = asymp_calc_scan.GetHypoTest()
        null_p_value_scan = asym_calc_result_scan.NullPValue()
        kstars.append(kstar)
        p0_values.append(null_p_value_scan)

        print(f'k*: {kstar} GeV/c, p0: {null_p_value_scan:.10f}')

    ## create a graph with the p0 values
    local_pvalue_graph = TGraph(len(kstars), np.array(kstars), np.array(p0_values))
    local_pvalue_graph.SetName('p0_values')
    local_pvalue_graph.SetTitle('; #it{k}* (GeV/#it{c}); Local p-value')
    # log Y axis
    local_pvalue_graph.SetMarkerStyle(20)
    local_pvalue_graph.SetMarkerColor(418)
    local_pvalue_graph.SetMarkerSize(0)
    local_pvalue_graph.SetLineColor(418)
    local_pvalue_graph.SetLineWidth(2)
    
    outdir.cd()
    local_pvalue_graph.Write()
    '''

    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    print(tc.RED+f'Fit significance'+tc.RESET)
    print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
    print(f'p-value = {p_value:.10f} +/- {p_value_err:.10f}')
    print(f'significance = {significance:.10f} +/- {significance_error:.10f}')
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    


if __name__ == '__main__':
    
    #sign = 'Matter'
    sign = 'Anti'

    outfile = TFile.Open(f'/home/galucia/antiLithium4/femto/output/fit_PrHe{sign}.root', 'recreate')

    infile_path_li4 = '/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root'
    li4_peak_hist = load_hist(HistLoadInfo(infile_path_li4, f'Correlations/fKstar{sign}'))
    me_hist = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
                                f'Correlations/fKstar{sign}'))
    li4_peak_params = prepare_peak_CF(li4_peak_hist, me_hist, outfile, use_mixed=True)

    outdir = outfile.mkdir(f'cent_integrated')
    infile_path_coulomb = f'/home/galucia/antiLithium4/analysis/output/CATS/CATS_new.root'
    coulomb_template_hist = load_hist(HistLoadInfo(infile_path_coulomb, 'hHe3_p_Coul_CF_LS'))
    coulomb_pdf = prepare_coulomb_CF(coulomb_template_hist, outdir)

    infile_path = '/home/galucia/antiLithium4/analysis/output/PbPb/studies.root'
    hist = load_hist(HistLoadInfo(infile_path, f'Correlation{sign}/hCorrelation_kstar'))
    fit_CF(hist, outdir, li4_peak_params, sign)
    compute_significance(outdir, sign)

    centralities_dotted = ['0.0_10.0', '10.0_30.0', '30.0_50.0']
    centralities = ['0_10', '10_30', '30_50']

    for cent, cent_dot in zip(centralities, centralities_dotted):

        outdir = outfile.mkdir(f'cent{cent}')
        infile_path_coulomb = f'/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent{cent}_new.root'
        coulomb_template_hist = load_hist(HistLoadInfo(infile_path_coulomb, 'hHe3_p_Coul_CF_LS'))
        coulomb_pdf = prepare_coulomb_CF(coulomb_template_hist, outdir)

        infile_path = '/home/galucia/antiLithium4/analysis/output/PbPb/studies.root'
        hist = load_hist(HistLoadInfo(infile_path, f'Correlation{sign}/hCorrelation_kstar_cent{cent_dot}'))
        fit_CF(hist, outdir, li4_peak_params, sign, cent)
        compute_significance(outdir, sign, cent)

    outfile.Close()
