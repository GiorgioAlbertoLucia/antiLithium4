'''
    Macro to fit the p-He3 correlation function
'''

from torchic.core.histogram import HistLoadInfo, load_hist
from ROOT import RooRealVar, RooDataHist, RooHistPdf, RooCrystalBall, RooFit, RooAddPdf, RooArgList
from ROOT import TFile, TH1F, TPaveText

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

def prepare_peak_CF(hist: TH1F, outfile: TFile) -> dict:

    kstar = RooRealVar('kstar', 'kstar', 0., 0.25)
    fit_params = {
        'cb_mean': RooRealVar('cb_mean', 'cb_mean', 0., 0.25, 'GeV/c^{2}'),
        'cb_sigma': RooRealVar('cb_sigma', 'cb_sigma', 0., 0.1, 'GeV/c^{2}'),
        'cb_aL': RooRealVar('cb_aL', 'cb_aL', -10., 0.),
        'cb_nL': RooRealVar('cb_nL', 'cb_nL', 0., 30.),
        'cb_aR': RooRealVar('cb_aR', 'cb_aR', -10., 0.),
        'cb_nR': RooRealVar('cb_nR', 'cb_nR', 0., 30.),
    }
    signal = RooCrystalBall('cb', 'cb', kstar, fit_params['cb_mean'], fit_params['cb_sigma'], fit_params['cb_aL'], 
                            fit_params['cb_nL'], fit_params['cb_aR'], fit_params['cb_nR'])

    datahist = RooDataHist('datahist', 'datahist', [kstar], hist)
    signal.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.2))
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

def fit_CF(hist: TH1F, outfile, coulomb_pdf: RooHistPdf, li4_peak_params: dict):

    outfile.cd()
    hist.Write('data_hist')

    kstar = RooRealVar('kstar', 'kstar', 0., 0.8)
    
    infile_path_coulomb = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent0_10_new.root'
    coulomb_template_hist = load_hist(HistLoadInfo(infile_path_coulomb, 'hHe3_p_Coul_CF_LS'))
    data_hist_coulomb = RooDataHist('coulomb', 'coulomb', [kstar], Import=coulomb_template_hist)
    coulomb_pdf = RooHistPdf('pdf', 'pdf', {kstar}, data_hist_coulomb, 0)

    data_hist = RooDataHist('CF', 'CF', [kstar], Import=hist)
    
    par_vals = {name: param.getVal() for name, param in li4_peak_params.items()}
    fit_params = {
                    'li4_mean': RooRealVar('li4_mean', 'li4_mean', par_vals['cb_mean'], 0., 0.25, 'GeV/c^{2}'),
                    'li4_sigma': RooRealVar('li4_sigma', 'li4_sigma', par_vals['cb_sigma'], par_vals['cb_sigma'], 4 * par_vals['cb_sigma'], 'GeV/c^{2}'),
                    'li4_cb_aL': RooRealVar('li4_cb_aL', 'li4_cb_aL', par_vals['cb_aL'], -10., 0.),
                    'li4_cb_nL': RooRealVar('li4_cb_nL', 'li4_cb_nL', par_vals['cb_nL'], 0., 30.),
                    'li4_cb_aR': RooRealVar('li4_cb_aR', 'li4_cb_aR', par_vals['cb_aR'], -10., 0.),
                    'li4_cb_nR': RooRealVar('li4_cb_nR', 'li4_cb_nR', par_vals['cb_nR'], 0., 30.),
                }

    fit_params['li4_mean'].setConstant(True)
    fit_params['li4_cb_aL'].setConstant(True)
    fit_params['li4_cb_nL'].setConstant(True)
    fit_params['li4_cb_aR'].setConstant(True)
    fit_params['li4_cb_nR'].setConstant(True)

    li4_signal = RooCrystalBall('cb', 'cb', kstar, fit_params['li4_mean'], fit_params['li4_sigma'], fit_params['li4_cb_aL'], 
                                fit_params['li4_cb_nL'], fit_params['li4_cb_aR'], fit_params['li4_cb_nR'])

    li4_fraction = RooRealVar('li4_fraction', 'li4_fraction', 0.1, 0., 1.)
    model = RooAddPdf('model', 'model', [li4_signal, coulomb_pdf], [li4_fraction])
    #model = RooAddPdf('model', 'model', [li4_signal], [li4_fraction])
    model.fitTo(data_hist, PrintLevel=-1)

    frame = kstar.frame()
    data_hist.plotOn(frame)
    model.plotOn(frame, Components={coulomb_pdf, li4_signal}, LineStyle=2, LineColor=2)

    text = TPaveText(0.7, 0.2, 0.9, 0.5, 'NDC')
    for name, param in fit_params.items():
        text.AddText(f'{name} = {param.getVal():.4f} +/- {param.getError():.4f}')
    frame.addObject(text)

    outfile.cd()
    frame.Write('fit_CF')
    


if __name__ == '__main__':
    
    outfile = TFile.Open('/home/galucia/antiLithium4/femto/output/fit_PrHe.root', 'recreate')

    infile_path_li4 = '/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root'
    li4_peak_hist = load_hist(HistLoadInfo(infile_path_li4, 'Correlations/fKstarAnti'))
    li4_peak_params = prepare_peak_CF(li4_peak_hist, outfile)

    centralities_dotted = ['0.0_10.0', '10.0_30.0', '30.0_50.0']
    centralities = ['0_10', '10_30', '30_50']

    for cent, cent_dot in zip(centralities, centralities_dotted):

        outdir = outfile.mkdir(f'cent{cent}')
        infile_path_coulomb = f'/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent{cent}_new.root'
        coulomb_template_hist = load_hist(HistLoadInfo(infile_path_coulomb, 'hHe3_p_Coul_CF_LS'))
        coulomb_pdf = prepare_coulomb_CF(coulomb_template_hist, outdir)

        infile_path = '/home/galucia/antiLithium4/analysis/output/PbPb/studies_noH3.root'
        hist = load_hist(HistLoadInfo(infile_path, f'CorrelationAnti/hCorrelation_kstar_cent{cent_dot}'))
        fit_CF(hist, outdir, coulomb_pdf, li4_peak_params)

    outfile.Close()
