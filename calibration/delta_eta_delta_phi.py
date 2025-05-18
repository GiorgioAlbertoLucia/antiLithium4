'''
    Calibration for the Delta eta vs Delta phi cut
'''

from torchic import Dataset, AxisSpec
from ROOT import TFile, TH1F, TF1, TCanvas, TPaveText

def fit_depletion(hist: TH1F, outfile: TFile):

    fit_func = TF1('fit_func', '[0] - gaus(1)', -0.05, 0.05)
    fit_func.SetParameter(0, 100)
    fit_func.SetParameter(1, 50)
    fit_func.SetParameter(2, 0.0)
    fit_func.SetParameter(3, 0.003)
    hist.Fit(fit_func, 'R')

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    hist.Draw()
    fit_func.Draw('same')

    text = TPaveText(0.7, 0.3, 0.9, 0.5, 'NDC')
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.AddText(f'Fit function: {fit_func.GetName()}')
    text.AddText(f'baseline = {fit_func.GetParameter(0):.3f} #pm {fit_func.GetParError(0):.3f}')
    text.AddText(f'N = {fit_func.GetParameter(1):.3f} #pm {fit_func.GetParError(1):.3f}')
    text.AddText(f'#mu = {fit_func.GetParameter(2):.3f} #pm {fit_func.GetParError(2):.3f}')
    text.AddText(f'#sigma = {fit_func.GetParameter(3):.3f} #pm {fit_func.GetParError(3):.3f}')
    text.Draw()

    outfile.cd()
    canvas.Write(hist.GetName()+ '_fit')

    return fit_func.GetParameter(3)


def delta_eta_delta_phi_calibration(dataset:Dataset, outfile:TFile):

    dataset['fDeltaEta'] = dataset['fEtaHe3'] - dataset['fEtaHad']
    dataset['fDeltaPhi'] = dataset['fPhiHe3'] - dataset['fPhiHad']

    axis_spec_deltaeta = AxisSpec(100, -0.05, 0.05, 'deltaeta', ';#Delta#eta;#Delta#phi')
    axis_spec_deltaphi = AxisSpec(100, -0.05, 0.05, 'deltaphi', ';#Delta#eta;#Delta#phi')

    h2 = dataset.build_th2('fDeltaEta', 'fDeltaPhi', axis_spec_deltaeta, axis_spec_deltaphi)
    h_delta_phi = h2.ProjectionY('delta_phi', h2.GetXaxis().FindBin(-0.005), h2.GetXaxis().FindBin(0.005))
    h_delta_phi.SetTitle(';#Delta#phi;')
    sigma_dphi = fit_depletion(h_delta_phi, outfile)
    h_delta_eta = h2.ProjectionX('delta_eta', h2.GetYaxis().FindBin(-0.005), h2.GetYaxis().FindBin(0.005))
    h_delta_eta.SetTitle(';#Delta#eta;')
    sigma_deta = fit_depletion(h_delta_eta, outfile)

    outfile.cd()
    h2.Write('delta_eta_delta_phi')
    h_delta_phi.Write('delta_phi')
    h_delta_eta.Write('delta_eta')

    return sigma_deta, sigma_dphi

def apply_cut(dataset:Dataset, sigma_deta:float, sigma_dphi:float, outfile:TFile):

    cut = f'(fDeltaEta/(2*{sigma_deta}))**2 + (fDeltaPhi/(2*{sigma_dphi}))**2 > 1'
    dataset.query(cut, inplace=True)

    axis_spec_deltaeta = AxisSpec(100, -0.05, 0.05, 'deltaeta', ';#Delta#eta;#Delta#phi')
    axis_spec_deltaphi = AxisSpec(100, -0.05, 0.05, 'deltaphi', ';#Delta#eta;#Delta#phi')
    h2 = dataset.build_th2('fDeltaEta', 'fDeltaPhi', axis_spec_deltaeta, axis_spec_deltaphi)
    h2.SetTitle(';#Delta#eta;#Delta#phi')
    outfile.cd()
    h2.Write('delta_eta_delta_phi_cut')

if __name__ == '__main__':

    infile_path = '/data/galucia/lithium_local/same/LHC23_PbPb_pass4_long_same_lsus.root'
    folder_name = 'DF*'
    tree_name = 'O2he3hadtable'
    dataset = Dataset.from_root(infile_path, folder_name=folder_name, tree_name=tree_name, columns=['fEtaHe3', 'fEtaHad', 'fPhiHe3', 'fPhiHad'])
    
    outfile = TFile('/home/galucia/antiLithium4/calibration/output/delta_eta_delta_phi_calibration.root', 'recreate')
    sigma_deta, sigma_dphi = delta_eta_delta_phi_calibration(dataset, outfile)
    apply_cut(dataset, sigma_deta, sigma_dphi, outfile)
    outfile.Close()