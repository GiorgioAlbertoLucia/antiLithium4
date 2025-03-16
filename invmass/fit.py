'''
    Macro to fit the p-He3 correlation function
'''

import numpy as np
from torchic.core.histogram import HistLoadInfo, load_hist
from torchic.utils.terminal_colors import TerminalColors as tc
from ROOT import RooRealVar, RooDataHist, RooHistPdf, RooCrystalBall, RooFit, RooAddPdf, RooGenericPdf, RooProdPdf, RooStats, RooWorkspace
from ROOT import TFile, TH1F, TPaveText, TCanvas, TDirectory

def prepare_peak_invmass(hist: TH1F, outfile: TFile) -> dict:

    mass = RooRealVar('mass', 'm(p + ^{3}He) (GeV/#it{c}^{2})', 3.747, 3.787)
    fit_params = {
        'cb_mean': RooRealVar('cb_mean', 'cb_mean', 3.751, 3.747, 3.787, 'GeV/c^{2}'),
        'cb_sigma': RooRealVar('cb_sigma', 'cb_sigma', 0.004, 0.001, 0.1, 'GeV/c^{2}'),
        'cb_aL': RooRealVar('cb_aL', 'cb_aL', 1, 0., 15),
        'cb_nL': RooRealVar('cb_nL', 'cb_nL', 3, 0., 15),
        'cb_aR': RooRealVar('cb_aR', 'cb_aR', 1, 0., 15),
        'cb_nR': RooRealVar('cb_nR', 'cb_nR', 3, 0., 15),
    }
    signal = RooCrystalBall('cb', 'cb', mass, fit_params['cb_mean'], fit_params['cb_sigma'], fit_params['cb_aL'], 
                            fit_params['cb_nL'], fit_params['cb_aR'], fit_params['cb_nR'])

    datahist = RooDataHist('datahist', 'datahist', [mass], hist)
    signal.fitTo(datahist, RooFit.Save(), RooFit.Range(3.747, 3.787), Extended=True)
    frame = mass.frame()
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

def fit_invmass(hist: TH1F, hist_me: TH1F, outfile, li4_peak_params: dict, sign: str, cent: str = ''):

    outfile.cd()
    hist.Write('data_hist')

    mass = RooRealVar('mass', 'm(p + ^{3}He) (GeV/#it{c}^{2})', 3.747, 3.847)
    data_hist_me = RooDataHist('mixed_event', 'mixed_event', [mass], Import=hist_me)
    bkg_pdf = RooHistPdf('bkg_pdf', 'bkg_pdf', {mass}, data_hist_me, 0)

    data_hist = RooDataHist('CF', 'CF', [mass], Import=hist)
    
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

    li4_signal = RooCrystalBall('li4_signal', 'li4_signal', mass, fit_params['li4_mean'], fit_params['li4_sigma'], fit_params['li4_cb_aL'], 
                                fit_params['li4_cb_nL'], fit_params['li4_cb_aR'], fit_params['li4_cb_nR'])

    li4_fraction = RooRealVar('li4_fraction', 'li4_fraction', 0.5, 0., 1e3)
    bkg_fraction = RooRealVar('coulomb_fraction', 'coulomb_fraction', 0.9, 0., 1e3)
    #model = RooAddPdf('model', 'model', [li4_signal, bkg_pdf], [li4_fraction, bkg_fraction])
    model = RooAddPdf('model', 'model', [bkg_pdf], [bkg_fraction])
    #model = RooAddPdf('model', 'model', [li4_signal], [li4_fraction])
    model.fitTo(data_hist, PrintLevel=-1, Extended=True)

    ##
    model = RooAddPdf('model', 'model', [li4_signal, bkg_pdf], [li4_fraction, bkg_fraction])

    frame = mass.frame()
    frame.SetYTitle('Counts')
    data_hist.plotOn(frame)
    model.plotOn(frame, LineStyle=1, LineColor=2)
    model.plotOn(frame, Components={li4_signal}, LineStyle=2, LineColor=4)
    model.plotOn(frame, Components={bkg_pdf}, LineStyle=2, LineColor=3)

    # fit significance
    ws = RooWorkspace('workspace')
    model_config = RooStats.ModelConfig('model_config', ws)
    model_config.SetPdf(model)
    model_config.SetParametersOfInterest({li4_fraction})
    model_config.SetObservables({mass})
    getattr(ws, 'import')(model_config)
    getattr(ws, 'import')(data_hist)
    outfile_path = outfile.GetFile().GetName()
    path, ext = outfile_path.split('.') 
    ws.writeToFile(f'{path}_ws_{sign}_{cent}.{ext}', True)

    text = TPaveText(0.55, 0.25, 0.75, 0.65, 'NDC')
    for name, param in fit_params.items():
        tmp_name = name.replace('li4_', '')
        tmp_name = tmp_name.replace('cb_', '')
        if (tmp_name != 'mean') and (tmp_name != 'sigma'):
            text.AddText(f'{param.GetTitle()} = {param.getVal():.4f} (fixed)')
        else:
            text.AddText(f'{param.GetTitle()} = {param.getVal():.4f} #pm {param.getError():.4f}')
    #text.AddText(f'li4_fraction = {li4_fraction.getVal():.4f} #pm {li4_fraction.getError():.4f}')
    #text.AddText(f'coulomb_fraction = {coulomb_fraction.getVal():.4f} #pm {coulomb_fraction.getError():.4f}')
    text.AddText(f'#chi^{{2}} / NDF = {frame.chiSquare():.4f}')
    text.SetBorderSize(0)
    text.SetFillStyle(0)
    #text.AddText(f'S_R = {S_R:.2f}')
    frame.addObject(text)

    outfile.cd()
    frame.Write('fit_invmass')

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    frame.Draw()
    canvas.SaveAs(f'/home/galucia/antiLithium4/invmass/output/fit_PrHe_{sign}_{cent}.pdf')

def compute_significance(infile: TDirectory, sign: str, cent: str = ''):
    
    outfile_path = infile.GetFile().GetName()
    path, ext = outfile_path.split('.') 
    ws_file = TFile.Open(f'{path}_ws_{sign}_{cent}.{ext}')
    ws = ws_file.Get('workspace')

    CF_datahist = ws.data('CF')
    model = ws.obj('model_config')
    poi = model.GetParametersOfInterest().first()
    model.SetSnapshot({poi})
    # create the b-only model
    bkg_model = model.Clone()
    bkg_model.SetName('bkg_pdf_config')
    poi.setVal(0)
    bkg_model.SetSnapshot(poi)
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
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    print(tc.RED+f'Fit significance'+tc.RESET)
    print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
    print(f'p-value = {p_value:.10f} +/- {p_value_err:.10f}')
    print(f'significance = {significance:.10f} +/- {significance_error:.10f}')
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    


if __name__ == '__main__':
    
    #sign = 'Matter'
    sign = 'Anti'

    outfile = TFile.Open(f'/home/galucia/antiLithium4/invmass/output/fit_PrHe{sign}.root', 'recreate')

    infile_path_li4 = '/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root'
    li4_peak_hist = load_hist(HistLoadInfo(infile_path_li4, f'InvMass/InvMass{sign}Li'))
    li4_peak_params = prepare_peak_invmass(li4_peak_hist, outfile)

    outdir = outfile.mkdir(f'cent_integrated')
    infile_path_me = f'/home/galucia/antiLithium4/analysis/output/PbPb/studies.root'
    hist_me = load_hist(HistLoadInfo(infile_path_me, f'InvariantMass{sign}/hMixed_invMass'))

    infile_path = '/home/galucia/antiLithium4/analysis/output/PbPb/studies.root'
    hist = load_hist(HistLoadInfo(infile_path, f'InvariantMass{sign}/hSame_invMass'))
    fit_invmass(hist, hist_me, outdir, li4_peak_params, sign)
    compute_significance(outdir, sign)

    centralities_dotted = ['0.0_10.0', '10.0_30.0', '30.0_50.0']
    centralities = ['0_10', '10_30', '30_50']

    for cent, cent_dot in zip(centralities, centralities_dotted):

        outdir = outfile.mkdir(f'cent{cent}')
        #infile_path_me = f'/home/galucia/antiLithium4/analysis/output/PbPb/studies.root'
        #hist_me = load_hist(HistLoadInfo(infile_path_me, f'InvariantMass{sign}/hMixed_invMass_{cent_dot.split("_")[0]}'))
        infile_path_me = f'/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root'
        hist_me = load_hist(HistLoadInfo(infile_path_me, f'InvMass/InvMass{sign}Li'))
        infile_path_coulomb = f'/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent{cent}_new.root'
        coulomb = load_hist(HistLoadInfo(infile_path_coulomb, 'hHe3_p_Coul_CF_LS'))
        for ibin in range(1, hist_me.GetNbinsX() + 1):
            imass = hist_me.GetBinCenter(ibin)
            ibin_coulomb = coulomb.FindBin(imass)
            hist_me.SetBinContent(ibin, hist_me.GetBinContent(ibin) * coulomb.GetBinContent(ibin_coulomb))
            hist_me.SetBinError(ibin, hist_me.GetBinError(ibin) * coulomb.GetBinContent(ibin_coulomb))

        infile_path = '/home/galucia/antiLithium4/analysis/output/PbPb/studies.root'
        hist = load_hist(HistLoadInfo(infile_path, f'InvariantMass{sign}/hSame_invMass_{cent_dot.split("_")[0]}'))
        fit_invmass(hist, hist_me, outdir, li4_peak_params, sign, cent)
        compute_significance(outdir, sign, cent)

    outfile.Close()
