'''
    Macro to fit the correlation function
'''

from ROOT import TPaveText, TFile
from ROOT import RooRealVar, RooCrystalBall, RooDataHist, RooArgList, RooFit
from torchic.core.histogram import load_hist, HistLoadInfo

def template_from_mc(hist_load_info: HistLoadInfo, output_file: TFile):

    hist = load_hist(hist_load_info)

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

    datahist = RooDataHist('datahist', 'datahist', RooArgList(kstar), hist)
    #signal.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.2))
    signal.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.2))
    frame = kstar.frame()
    datahist.plotOn(frame)
    signal.plotOn(frame)

    text = TPaveText(0.7, 0.7, 0.9, 0.9, 'NDC')
    for name, param in fit_params.items():
        text.AddText(f'{name} = {param.getVal():.4f} +/- {param.getError():.4f}')
    frame.addObject(text)
    frame.Draw()
    output_file.cd()
    frame.Write(f'frame')

    return fit_params


if __name__ == '__main__':

    mc_hist_info = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root', 
                                'Correlations/fKstarAnti')
    output_file = TFile.Open('/home/galucia/antiLithium4/femto/output/li4_peak.root', 'RECREATE')
    fit_params = template_from_mc(mc_hist_info, output_file)
    output_file.Close()

    