'''
    Macro to fit the correlation function
'''

import numpy as np
from ROOT import TPaveText, TFile, TCanvas
from ROOT import RooRealVar, RooCrystalBall, RooDataHist, RooArgList, RooFit, RooGenericPdf, RooAddPdf
from torchic.core.histogram import load_hist, HistLoadInfo

def template_from_mc(hist_mc_load_info: HistLoadInfo, hist_me_load_info: HistLoadInfo, output_file: TFile):

    hist_mc = load_hist(hist_mc_load_info)
    hist_me = load_hist(hist_me_load_info)
    hist_mc.Add(hist_me)
    hist = hist_mc.Clone()
    for ibin in range(1, hist_mc.GetNbinsX() + 1):
        value = hist_mc.GetBinContent(ibin) / hist_me.GetBinContent(ibin) #- 1
        error = value * np.sqrt(hist_mc.GetBinError(ibin)**2 / hist_mc.GetBinContent(ibin)**2 + hist_me.GetBinError(ibin)**2 / hist_me.GetBinContent(ibin)**2)
        hist.SetBinContent(ibin, value)
        hist.SetBinError(ibin, error)

    kstar = RooRealVar('kstar', '#it{k}* (GeV/#it{c})', 0., 0.25)
    fit_params = {
        'cb_mean': RooRealVar('cb_mean', '#mu', 0., 0.25, 'GeV/#it{c}'),
        'cb_sigma': RooRealVar('cb_sigma', '#sigma', 0., 0.1, 'GeV/#it{c}'),
        'cb_aL': RooRealVar('cb_aL', 'a_{L}', 0.1, 10., ''),
        'cb_nL': RooRealVar('cb_nL', 'n_{L}', 0.1, 30., ''),
        'cb_aR': RooRealVar('cb_aR', 'a_{R}', 0.1, 10., ''),
        'cb_nR': RooRealVar('cb_nR', 'n_{R}', 0.1, 30., ''),
    }
    signal = RooCrystalBall('cb', 'cb', kstar, fit_params['cb_mean'], fit_params['cb_sigma'], fit_params['cb_aL'], 
                            fit_params['cb_nL'], fit_params['cb_aR'], fit_params['cb_nR'])

    datahist = RooDataHist('datahist', 'datahist', RooArgList(kstar), hist)
    #signal.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.2))
    #signal.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.2))
    frame = kstar.frame()
    datahist.plotOn(frame)
    #signal.plotOn(frame)

    #text = TPaveText(0.5, 0.5, 0.8, 0.8, 'NDC')
    #for name, param in fit_params.items():
    #    text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
    ##text.AddText(f'#chi^{{2}} / NDF = {frame.chiSquare():.2f}')
    #text.SetFillColor(0)
    #text.SetBorderSize(0)
    #frame.addObject(text)
    frame.Draw()

    output_file.cd()
    frame.Write(f'frame')
    hist_mc.Write(f'hist_mc')
    hist_me.Write(f'hist_me')
    hist.Write(f'hist')

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    canvas.DrawFrame(0., 0.1, 0.25, 4e2, ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
    frame.Draw('same')
    canvas.SetLogy()
    canvas.SaveAs('/home/galucia/antiLithium4/femto/output/li4_peak_before_subtraction.pdf')

    return fit_params


if __name__ == '__main__':

    mc_hist_info = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root', 
                                'Correlations/fKstarAnti')
    me_hist_info = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
                                'Correlations/fKstarAnti')
    output_file = TFile.Open('/home/galucia/antiLithium4/femto/output/li4_peak.root', 'RECREATE')
    fit_params = template_from_mc(mc_hist_info, me_hist_info, output_file)
    output_file.Close()

    