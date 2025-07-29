'''
    Macro to fit the correlation function
'''

import numpy as np
from ROOT import TPaveText, TFile, TCanvas, gStyle, TF1, TH1F, gInterpreter
from ROOT import RooRealVar, RooCrystalBall, RooDataHist, RooArgList, RooFit, RooGenericPdf, RooAddPdf, RooGaussian, RooHistPdf, RooFFTConvPdf, RooNumConvPdf, RooVoigtian
from torchic.core.histogram import load_hist, HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc

import os
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
INCLUDE_DIR = os.path.join(CURRENT_DIR, 'include')
gInterpreter.ProcessLine(f'#include "{INCLUDE_DIR}/RooCustomVoigtian.h"')
gInterpreter.ProcessLine(f'#include "{INCLUDE_DIR}/RooSillPdf.h"')
from ROOT import RooCustomVoigtian, RooSillPdf

def combine_histograms(h1: TH1F, h2: TH1F):
    '''
        Bin by bin addition
    '''

    h_combined = h1.Clone()
    for ibin in range(1, h1.GetNbinsX() + 1):
        value = h1.GetBinContent(ibin) + h2.GetBinContent(ibin)
        error = np.sqrt(h1.GetBinError(ibin)**2 + h2.GetBinError(ibin)**2)
        h_combined.SetBinContent(ibin, value)
        h_combined.SetBinError(ibin, error)
    return h_combined

def correlation_function(h_signal_mc: TH1F, h_mixed_event: TH1F, substract_one: bool = False):

    h_same_event = h_signal_mc.Clone()
    h_same_event = combine_histograms(h_signal_mc, h_mixed_event)
    h_correlation_function = h_same_event.Clone()
    for ibin in range(1, h_same_event.GetNbinsX() + 1):
        value = h_same_event.GetBinContent(ibin) / h_mixed_event.GetBinContent(ibin)
        if substract_one:
            value -= 1
        check = h_same_event.GetBinContent(ibin) > 0 and h_mixed_event.GetBinContent(ibin) > 0
        error = value * np.sqrt(h_same_event.GetBinError(ibin)**2 / h_same_event.GetBinContent(ibin)**2 + h_mixed_event.GetBinError(ibin)**2 / h_mixed_event.GetBinContent(ibin)**2) if check else 0
        h_correlation_function.SetBinContent(ibin, value)
        h_correlation_function.SetBinError(ibin, error)

    return h_correlation_function, h_same_event

def draw_same_mixed_event(h_signal_mc: TH1F, h_mixed_event: TH1F, h_same_event: TH1F, output_file: TFile):

    canvas_hist = TCanvas('canvas_hist', 'canvas_hist', 800, 600)
    hframe = canvas_hist.DrawFrame(0., 1e-4, 0.25, 0.12, ';#it{k}* (GeV/#it{c});Counts (a.u.)')
    canvas_hist.SetLogy()

    h_same_event.SetLineColor(797)
    h_same_event.SetMarkerColor(797)
    h_same_event.SetMarkerStyle(22)
    h_same_event.Draw('hist e0 same')
    h_signal_mc.SetLineColor(601)
    h_signal_mc.SetMarkerColor(601)
    h_signal_mc.SetMarkerStyle(21)
    h_signal_mc.Draw('hist e0 same')
    h_mixed_event.SetLineColor(418)
    h_mixed_event.SetMarkerColor(418)
    h_mixed_event.SetMarkerStyle(23)
    h_mixed_event.Draw('hist e0 same')
    canvas_hist.SaveAs('/home/galucia/antiLithium4/femto/output/li4_peak_hist.pdf')

    output_file.cd()
    h_same_event.Write(f'h_same_event')
    h_signal_mc.Write(f'h_signal_mc')
    h_mixed_event.Write(f'h_mixed_event')

def template_from_mc(h_signal_mc: TH1F, h_mixed_event: TH1F, output_file: TFile):

    gStyle.SetOptStat(0)
    h_signal_mc.Rebin(5) # matching rebin
    h_correlation_function, h_same_event = correlation_function(h_signal_mc, h_mixed_event)

    draw_same_mixed_event(h_signal_mc, h_mixed_event, h_same_event, output_file)
    output_file.cd()
    h_correlation_function.Write(f'h_correlation_function')

    #h_correlation_function_substracted, _ = correlation_function(h_signal_mc, h_mixed_event, substract_one=True)
    

def fit_mc_signal(h_signal_mc_load_info: HistLoadInfo, output_file: TFile = None) -> dict:

    print(f'\n{tc.BOLD+tc.GREEN}Fitting MC signal with Crystal Ball function{tc.RESET}')
    mass = RooRealVar('mass', 'm (GeV/#it{c}^{2})', 3.727, 3.807)
    h_mc_inv_mass = load_hist(h_signal_mc_load_info)
    mc_invmass_datahist = RooDataHist('bkg_dh', 'bkg_dh', [mass], Import=h_mc_inv_mass)

    mc_signal_pars = {
        'mean': RooRealVar('mean', 'mean', 3.751, 3.747, 3.767),
        'sigma': RooRealVar('sigma', 'sigma', 0.003, 0.0001, 0.1),
        'aL': RooRealVar('aL', 'aL', 0.01, 10.),
        'nL': RooRealVar('nL', 'nL', 0.1, 10.),
        'aR': RooRealVar('aR', 'aR', 0.01, 10.),
        'nR': RooRealVar('nR', 'nR', 0.1, 10.)
    }
    mc_signal_fit = RooCrystalBall('mc_signal', 'mc_signal', mass, mc_signal_pars['mean'], mc_signal_pars['sigma'],
                                    mc_signal_pars['aL'], mc_signal_pars['nL'],
                                    #doubleSided=True)
                                    mc_signal_pars['aR'], mc_signal_pars['nR'])
    mc_signal_fit.fitTo(mc_invmass_datahist, RooFit.Save())

    frame = mass.frame()
    mc_invmass_datahist.plotOn(frame, MarkerStyle=20, MarkerColor=418, LineColor=418, DrawOption='e0')
    mc_signal_fit.plotOn(frame, LineColor=601, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    mc_signal_fit.paramOn(frame, Label='Crystal Ball fit', Format='NEU', ShowConstants=True, Layout=(0.15, 0.45, 0.7))
    frame.SetTitle(';m (GeV/#it{c}^{2});Counts (a.u.)')
    canvas = TCanvas('canvas', 'canvas', 800, 600)
    gStyle.SetOptStat(0)
    canvas.cd()
    frame.Draw('same')
    canvas.SetLogy()
    output_file.cd('sampling')
    canvas.Write('mc_signal_fit_canvas')

    del mc_invmass_datahist, h_mc_inv_mass, mc_signal_fit
    return mc_signal_pars

def fit_mc_signal_with_custom_voigtian(h_signal_mc_load_info: HistLoadInfo, output_file: TFile = None) -> dict:

    print(f'\n{tc.BOLD+tc.GREEN}Fitting MC signal with Custom Voigtian function{tc.RESET}')
    mass = RooRealVar('mass', 'm (GeV/#it{c}^{2})', 3.727, 3.807)
    h_mc_inv_mass = load_hist(h_signal_mc_load_info)
    mc_invmass_datahist = RooDataHist('bkg_dh', 'bkg_dh', [mass], Import=h_mc_inv_mass)

    voigtian_pars = {
        'mean': RooRealVar('mean', 'mean', 3.751, 3.747, 3.767),
        'sigma': RooRealVar('sigma', 'sigma', 0.003, 0.0001, 0.1),
        'aL': RooRealVar('aL', 'aL', 2, 0.01, 10.),
        'nL': RooRealVar('nL', 'nL', 0.095, 0.01, 10.),
    }
    Gamma = RooRealVar('Gamma', 'Gamma', 0.)  # Intrinsic width
    Gamma.setConstant(True)  # Set intrinsic width to zero for the Voigtian fit

    voigtian_fit = RooCustomVoigtian('voigtian_fit', 'voigtian_fit', mass, voigtian_pars['mean'], Gamma, voigtian_pars['sigma'],
                                     voigtian_pars['aL'], voigtian_pars['nL'])
    voigtian_fit.fitTo(mc_invmass_datahist, RooFit.Save())

    frame = mass.frame()
    mc_invmass_datahist.plotOn(frame, MarkerStyle=20, MarkerColor=418, LineColor=418, DrawOption='e0')
    voigtian_fit.plotOn(frame, LineColor=601, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    voigtian_fit.paramOn(frame, Label='Voigtian fit', Format='NEU', ShowConstants=True, Layout=(0.15, 0.45, 0.7))
    frame.SetTitle(';m (GeV/#it{c}^{2});Counts (a.u.)')
    canvas = TCanvas('canvas', 'canvas', 800, 600)
    gStyle.SetOptStat(0)
    canvas.cd()
    frame.Draw('same')
    canvas.SetLogy()
    output_file.cd('sampling')
    canvas.Write('voigtian_fit')

    del mc_invmass_datahist, h_mc_inv_mass, voigtian_fit
    return voigtian_pars

def fit_mc_signal_with_sill(h_signal_mc_load_info: HistLoadInfo, output_file: TFile = None) -> dict:

    print(f'\n{tc.BOLD+tc.GREEN}Fitting MC signal with Custom Voigtian function{tc.RESET}')
    mass = RooRealVar('mass', 'm (GeV/#it{c}^{2})', 3.727, 3.807)
    h_mc_inv_mass = load_hist(h_signal_mc_load_info)
    mc_invmass_datahist = RooDataHist('bkg_dh', 'bkg_dh', [mass], Import=h_mc_inv_mass)

    M_PR = 0.938272 # GeV/c^2
    M_HE = 2.80839 # GeV/c^2
    sill_pars = {
        'mean': RooRealVar('mean', 'mean', 3.751, 3.747, 3.767),
        'gamma': RooRealVar('gamma', 'gamma', 0.003, 0.0001, 0.1),
        'eth': RooRealVar('eth', 'eth', M_HE+M_PR),  # Eth is the threshold energy
    }
    sill_pars['eth'].setConstant(True)  # Set threshold energy to constant

    sill_fit = RooSillPdf('sill_fit', 'sill_fit', mass, sill_pars['mean'], sill_pars['gamma'],
                                     sill_pars['eth'])
    sill_fit.fitTo(mc_invmass_datahist, RooFit.Save())

    frame = mass.frame()
    mc_invmass_datahist.plotOn(frame, MarkerStyle=20, MarkerColor=418, LineColor=418, DrawOption='e0')
    sill_fit.plotOn(frame, LineColor=601, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    sill_fit.paramOn(frame, Label='Voigtian fit', Format='NEU', ShowConstants=True, Layout=(0.15, 0.45, 0.7))
    frame.SetTitle(';m (GeV/#it{c}^{2});Counts (a.u.)')
    canvas = TCanvas('canvas', 'canvas', 800, 600)
    gStyle.SetOptStat(0)
    canvas.cd()
    frame.Draw('same')
    canvas.SetLogy()
    output_file.cd('sampling')
    canvas.Write('sill_fit')

    del mc_invmass_datahist, h_mc_inv_mass, sill_fit
    return sill_pars

def mc_with_intrinsic_width(h_signal_mc_load_info: HistLoadInfo, output_file: TFile = None):

    mc_signal_pars = fit_mc_signal(h_signal_mc_load_info, output_file)
    voigtian_pars = fit_mc_signal_with_custom_voigtian(h_signal_mc_load_info, output_file)
    sill_pars = fit_mc_signal_with_sill(h_signal_mc_load_info, output_file)
    mass = RooRealVar('mass', 'm (GeV/#it{c}^{2})', 3.6, 3.8)
    mass.setBins(1000, 'cache')

    mc_signal_shape = RooCrystalBall('mc_signal_shape', 'mc_signal_shape', mass, mc_signal_pars['mean'], mc_signal_pars['sigma'],
                                            mc_signal_pars['aL'], mc_signal_pars['nL'], 
                                            #doubleSided=True)
                                            mc_signal_pars['aR'], mc_signal_pars['nR'])
    
    for param in mc_signal_pars.values():
        param.setConstant(True)
    
    LI_INTRINSIC_WIDTH = 0.003 # GeV/c^2 
    LI_MASS = 3.751 # GeV/c^2
    mean = RooRealVar('mean', 'mean', LI_MASS, LI_MASS - 0.2, LI_MASS + 0.2)
    sigma = RooRealVar('sigma', 'sigma', LI_INTRINSIC_WIDTH, 0.001, 0.1)
    mean.setConstant(True)
    sigma.setConstant(True)
    gaussian_smearing = RooGaussian('gauss', 'gauss', mass, mean, sigma)
    
    convoluted_pdf = RooVoigtian('conv_pdf', 'conv_pdf', mass, mc_signal_pars['mean'], sigma, mc_signal_pars['sigma'])
    convoluted_pdf_custom = RooCustomVoigtian('conv_pdf', 'conv_pdf', mass, mc_signal_pars['mean'], sigma, mc_signal_pars['sigma'], mc_signal_pars['aL'], 
                                               mc_signal_pars['nL'])
    
    sill_pdf = RooSillPdf('sill_pdf', 'sill_pdf', mass, sill_pars['mean'], sill_pars['gamma'], sill_pars['eth'])

    _sigma = mc_signal_pars['sigma'].getVal()
    mc_signal_pars['sigma'].setVal(np.sqrt(_sigma*_sigma + LI_INTRINSIC_WIDTH*LI_INTRINSIC_WIDTH))
    sill_pars['gamma'].setVal(np.sqrt(sill_pars['gamma'].getVal()**2 + LI_INTRINSIC_WIDTH**2))
    
    sample_size = 100_000
    sample_crystal_ball = mc_signal_shape.generate([mass], sample_size)
    sample_convolution = convoluted_pdf_custom.generate([mass], sample_size)
    sample_sill = sill_pdf.generate([mass], sample_size)

    h_sample_convolution = sample_convolution.createHistogram('h_sample_convolution', mass, RooFit.Binning(1000))
    h_sample_crystal_ball = sample_crystal_ball.createHistogram('h_sample_crystal_ball', mass, RooFit.Binning(1000))
    h_sample_sill = sample_sill.createHistogram('h_sample_sill', mass, RooFit.Binning(1000))

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    gStyle.SetOptStat(0)
    frame = mass.frame()
    gaussian_smearing.plotOn(frame, LineColor=418, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    mc_signal_shape.plotOn(frame, LineColor=601, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    convoluted_pdf.plotOn(frame, LineColor=797, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    convoluted_pdf_custom.plotOn(frame, LineColor=800, LineWidth=2, MarkerStyle=20, DrawOption='l e0')
    frame.SetTitle(';m (GeV/#it{c}^{2});Counts (a.u.)')

    canvas.cd()
    frame.Draw('same')
    canvas.SetLogy()
    
    canvas.SaveAs('/home/galucia/antiLithium4/femto/output/li4_intrinsic_width.pdf')

    output_file.cd('sampling')
    convoluted_pdf.Write('conv_pdf')
    canvas.Write('canvas')
    frame.Write('frame_fit')
    h_sample_convolution.Write('h_sample_convolution')
    h_sample_crystal_ball.Write('h_sample_crystal_ball')
    h_sample_sill.Write('h_sample_sill')

    #return h_sample_convolution
    #return h_sample_crystal_ball
    return h_sample_sill

def mass_to_kstar(mass: float) -> float:
    """
    Convert k* to mass using the relation m^2 = k*^2 + m^2_pi.
    """
    M_PR = 0.938272 # GeV/c^2
    M_HE = 2.80839 # GeV/c^2
    lambda_func = lambda x, y, z: x**2 + y**2 + z**2 -2*x*y - 2*x*z - 2*y*z
    numerator = np.sqrt(lambda_func(mass**2, M_PR**2, M_HE**2)) if lambda_func(mass**2, M_PR**2, M_HE**2) >= 0 else 0
    denominator = 2 * mass
    return numerator / denominator

def kstar_resonance_sampling(h_sampled_mass: TH1F, output_file: TFile = None, sample_size: int = 100_000):
    """
    Sample k* resonance from a given histogram.
    """

    h_sample_kstar = TH1F('h_sample_kstar', ';k^{*} (GeV/#it{c}); Counts (a.u.)', 400, 0, 0.4)
    for iter in range(sample_size):
        mass = h_sampled_mass.GetRandom()
        kstar = mass_to_kstar(mass)
        h_sample_kstar.Fill(kstar)
    
    output_file.cd('sampling')
    h_sample_kstar.Write('h_sample_kstar')
    return h_sample_kstar




if __name__ == '__main__':

    mc_hist_info = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root', 
                                'Correlations/fKstarAnti')
    mc_invmass_hist_info = HistLoadInfo('/home/galucia/antiLithium4/root_dataframe/output/mc.root',
                                        'InvariantMassAntimatter/hInvariantMassAntimatter')
    me_hist_info = HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
                                'Correlations/fKstarAnti')
    
    #output_file = TFile.Open('/home/galucia/antiLithium4/femto/output/li4_peak.root', 'RECREATE')
    #fit_params = template_from_mc(mc_hist_info, me_hist_info, output_file)
    #output_file.Close()

    output_file_intrinsic = TFile.Open('/home/galucia/antiLithium4/femto/output/li4_intrinsic_width.root', 'RECREATE')
    output_file_intrinsic.mkdir('sampling')
    h_sample = mc_with_intrinsic_width(mc_invmass_hist_info, output_file_intrinsic)
    h_sample_kstar = kstar_resonance_sampling(h_sample, output_file_intrinsic, 10_000_000)

    h_mixed_event = load_hist(me_hist_info)
    h_mixed_event.Scale(1. / h_mixed_event.Integral())
    h_sample_kstar.Scale(1. / (1e2 * h_sample_kstar.Integral()))
    template_from_mc(h_sample_kstar, h_mixed_event, output_file_intrinsic)

    output_file_intrinsic.Close()