'''
    Macro to fit the p-He3 correlation function
'''

import numpy as np
import argparse
from contextlib import contextmanager
from torchic.core.histogram import HistLoadInfo, load_hist

from torchic.utils.terminal_colors import TerminalColors as tc
from ROOT import RooRealVar, RooDataHist, RooHistPdf, RooCrystalBall, RooFit, RooAddPdf, RooGenericPdf, RooProdPdf, RooStats, RooWorkspace, RooExtendPdf
from ROOT import TFile, TH1F, TPaveText, TCanvas, TDirectory, TGraph

from fit_utils import write_params_to_text

SUFFIX_DICT = {
    'dot':
    {
        'integrated': '',
        '0-10': '_cent0.0_10.0',
        '10-30': '_cent10.0_30.0',
        '30-50': '_cent30.0_50.0',
    },
    'undot':
    {
        'integrated': '',
        '0-10': '_cent0_10',
        '10-30': '_cent10_30',
        '30-50': '_cent30_50',
    }
}

CORAL_DICT = {
    'integrated': '3.6',
    '0-10': '4.2',
    '10-30': '3.4',
    '30-50': '2.6',

}

class PdfCanvas:
    def __init__(self, filename):
        self.filename = filename
        self.canvas = TCanvas('canvas', 'canvas', 800, 600)

    def __enter__(self):
        self.canvas.Print(f'{self.filename}(')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.canvas.Clear()
        self.canvas.Print(f'{self.filename})')

    def draw(self, obj):
        self.canvas.Clear()
        self.canvas.cd()
        obj.Draw()
        self.canvas.Print(self.filename)

def prepare_signal_fit(roows: RooWorkspace, h_signal: TH1F, h_mixed_event: TH1F, outfile: TFile):
    '''
        Fit the signal peak with a Crystal Ball function
    '''

    h_signal.Add(h_mixed_event)
    h_signal.Divide(h_mixed_event)
    h_unif = h_signal.Clone()
    for ibin in range(1, h_signal.GetNbinsX() + 1): h_unif.SetBinContent(ibin, 1.)
    h_signal.Add(h_unif, -1.)

    kstar = roows.var('kstar')
    fit_params = {
        'cb_mean': RooRealVar('cb_mean', '#mu', 0., 0.25, 'GeV/c'),
        'cb_sigma': RooRealVar('cb_sigma', '#sigma', 0.001, 0.1, 'GeV/c'),
        'cb_aL': RooRealVar('cb_aL', 'a_{L}', 1.3, 0.1, 10.),
        'cb_nL': RooRealVar('cb_nL', 'n_{L}', 2.7, 0.1, 10.),
        'cb_aR': RooRealVar('cb_aR', 'a_{R}', 1.1, 0.1, 10.),
        'cb_nR': RooRealVar('cb_nR', 'n_{R}', 5.7, 0.1, 10.),
    }
    signal_pdf = RooCrystalBall('signal_pdf', 'signal_pdf', kstar, fit_params['cb_mean'], fit_params['cb_sigma'], fit_params['cb_aL'], 
                            fit_params['cb_nL'], fit_params['cb_aR'], fit_params['cb_nR'])

    datahist = RooDataHist('signal_dh', 'signal_dh', [kstar], Import=h_signal)
    signal_pdf.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.25))
    frame = kstar.frame()
    datahist.plotOn(frame)
    signal_pdf.plotOn(frame)

    text = write_params_to_text(fit_params.values(), coordinates=(0.7, 0.3, 0.9, 0.5))
    frame.addObject(text)

    outfile.cd()
    frame.Write(f'signal_mc')
    getattr(roows, 'import')(signal_pdf)
    getattr(roows, 'import')(datahist)

def load_coulomb_template(centrality_opt: str, model: str) -> TH1F:
    '''
        Load the Coulomb template
    '''
    h_coulomb = None
    if model == 'CATS':
        cent = SUFFIX_DICT['undot'][centrality_opt]
        h_coulomb = load_hist(HistLoadInfo(f'/home/galucia/antiLithium4/analysis/output/CATS/CATS{cent}_new.root', 
                                                           'hHe3_p_Coul_CF_LS'))
    elif model == 'CorAL':
        h_coulomb = load_hist(HistLoadInfo(f'/home/galucia/antiLithium4/analysis/output/CorAL/sqwell_correlation.root', 
                                                           f'radius_{CORAL_DICT[centrality_opt]}fm/CF_{CORAL_DICT[centrality_opt]}fm'))
    return h_coulomb

def prepare_coulomb_template(roows: RooWorkspace, h_coulomb: TH1F, outfile: TDirectory):
    '''
        Get the Coulomb template for the p-He3 correlation function
    '''

    datahist = RooDataHist('coulomb_dh', 'coulomb_dh', [kstar], Import=h_coulomb)
    coulomb_pdf = RooHistPdf('coulomb_pdf', 'coulomb_pdf', [kstar], datahist)

    frame = kstar.frame()
    datahist.plotOn(frame)
    coulomb_pdf.plotOn(frame)

    outfile.cd()
    frame.Write(f'coulomb_pdf')
    getattr(roows, 'import')(coulomb_pdf)

def pull_from_coulomb(h_CF: TH1F, h_coulomb: TH1F, outdir: TDirectory):
    '''
        Pull the Coulomb template from the histogram
    '''

    pull = h_CF.Clone('pull')
    for ibin in range(1, h_CF.GetNbinsX() + 1):
        ikstar = h_CF.GetBinCenter(ibin)
        imeasured = h_CF.GetBinContent(ibin)
        ierror = h_CF.GetBinError(ibin)
        iexpected = h_coulomb.GetBinContent(h_coulomb.FindBin(ikstar))
        ipull = 0. if ierror < 1.e-9 else (imeasured - iexpected) / ierror
        pull.SetBinContent(ibin, ipull)
        pull.SetBinError(ibin, 1.)

    pull.SetTitle(';#it{k}* (GeV/#it{c});Pull')
    pull.SetMarkerStyle(20)
    pull.SetMarkerColor(418)

    outdir.cd()
    pull.Write(f'pull')

def prefit_CF(roows: RooWorkspace, model, h_CF: TH1F, outdir: TDirectory) -> RooRealVar:
    '''
        Fit the p-He3 correlation function
    '''

    kstar = roows.var('kstar')

    datahist = RooDataHist('datahist', 'datahist', [kstar], Import=h_CF)
    model.fitTo(datahist, RooFit.Save(), RooFit.Range(0.25, 0.8))

    frame = kstar.frame()
    datahist.plotOn(frame)
    model.plotOn(frame, Components={'coulomb_pdf'}, LineStyle=2, LineColor=3)

    text = write_params_to_text(model.getParameters(datahist), coordinates=(0.7, 0.3, 0.9, 0.5))
    frame.addObject(text)

    outdir.cd()
    frame.Write(f'prefit_CF')

def prepare_model_for_significance(roows: RooWorkspace, model: RooAddPdf, sig_counts: RooRealVar, kstar: RooRealVar):
    '''
        Prepare the model for significance estimation
    '''

    model_config = RooStats.ModelConfig('model_config', roows)
    model_config.SetPdf(model)
    model_config.SetParametersOfInterest({sig_counts})
    model_config.SetObservables({kstar})
    getattr(roows, 'import')(model_config)

def fit_CF(roows: RooWorkspace, h_CF: TH1F, outdir: TDirectory, pdf_canvas: PdfCanvas):
    '''
        Fit the p-He3 correlation function
    '''

    print(tc.GREEN + f'Fitting the correlation function' + tc.RESET)

    kstar = roows.var('kstar')
    
    bkg_pdf = roows.pdf('coulomb_pdf')
    signal_pdf = roows.pdf('signal_pdf')
    signal_dh = roows.data('signal_dh')

    sig_params = signal_pdf.getParameters(signal_dh)
    for param in sig_params:
        if 'mean' in param.GetName():
            #param.setRange(param.getVal() - 8 * param.getError(), param.getVal() + 5 * param.getError())
            param.setRange(0.06, param.getVal() + 5 * param.getError())
        elif 'sigma' in param.GetName():
            param.setRange(param.getVal(), 5 * param.getVal())
        else:
            param.setConstant(True)
        print(f'{param.GetName()} = {param.getVal()}')

    datahist = RooDataHist('datahist', 'datahist', [kstar], Import=h_CF)
    getattr(roows, 'import')(datahist)
    
    sig_counts = RooRealVar('sig_counts', 'sig_counts', 0.05, 0., 3)#1e3)
    bkg_counts = RooRealVar('bkg_counts', 'bkg_counts', h_CF.Integral(), 0., 1e6)
    model = RooAddPdf('model', 'model', [signal_pdf, bkg_pdf], [sig_counts, bkg_counts])

    sig_counts.setVal(0.)
    sig_counts.setConstant(True)
    prefit_CF(roows, model, h_CF, outdir)
    bkg_counts.setConstant(True)
    sig_counts.setConstant(False)

    model.fitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.8))

    frame = kstar.frame()
    datahist.plotOn(frame)
    model.plotOn(frame)
    model.plotOn(frame, Components={'signal_pdf'}, LineStyle=2, LineColor=2)
    model.plotOn(frame, Components={'coulomb_pdf'}, LineStyle=2, LineColor=3)
    
    text = write_params_to_text(sig_params, coordinates=(0.5, 0.2, 0.8, 0.5))
    text.AddText(f'sig counts = ({sig_counts.getVal():.4f} #pm {sig_counts.getError():.4f})')
    text.AddText(f'bkg counts = ({bkg_counts.getVal():.4f} #pm {bkg_counts.getError():.4f})')
    text.AddText(f'#chi^{{2}} / NDF = {frame.chiSquare():.4f}')
    frame.addObject(text)

    prepare_model_for_significance(roows, model, sig_counts, kstar)

    pdf_canvas.draw(frame)

    outdir.cd()
    frame.Write(f'fit_CF')

def evaluate_pvalue_and_significance(roows: RooWorkspace, outdir: TDirectory):

    datahist = roows.data('datahist')
    model = roows.obj('model_config')
    poi = model.GetParametersOfInterest().first()
    model.SetSnapshot({poi})
    # create the b-only model
    bkg_model = model.Clone()
    bkg_model.SetName('bkg_model')
    
    old_poi_val = poi.getVal()
    poi.setVal(0)
    bkg_model.SetSnapshot(poi)
    poi.setVal(old_poi_val)

    bkg_model.Print()
    roows.var('cb_mean').setConstant(True)
    roows.var('cb_sigma').setConstant(True)

    asymp_calc = RooStats.AsymptoticCalculator(datahist, model, bkg_model)
    asymp_calc.SetPrintLevel(0)
    asymp_calc_result = asymp_calc.GetHypoTest()
    p_value = asymp_calc_result.NullPValue()
    p_value_err = asymp_calc_result.NullPValueError()
    significance = asymp_calc_result.Significance()
    significance_error = asymp_calc_result.SignificanceError()
    
    ### perform a scan in kstar and compute the significance

    kstars = []
    p0_values = []
    p0_values_expected = []
    kstar_array = np.linspace(roows.var('kstar').getMin(), roows.var('kstar').getMax(), 40)
    for kstar in kstar_array:

        roows.var('kstar').setVal(kstar)
        roows.var('kstar').setConstant(True)
        asymp_calc_scan = RooStats.AsymptoticCalculator(datahist, model, bkg_model)
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


    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    print(tc.RED+f'Fit significance'+tc.RESET)
    print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
    print(f'p-value = {p_value:.10f} +/- {p_value_err:.10f}')
    print(f'significance = {significance:.10f} +/- {significance_error:.10f}')
    print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
    input()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
    parser.add_argument('--sign', dest='sign',
                        help='Perform the fit on either matter (Matter) or antimatter (Anti).', default='Anti')
    parser.add_argument('--model', dest='model',
                        help='Model to use for the Coulomb template.', default='CATS')
    args = parser.parse_args()

    outfile = TFile(f'output/{args.sign}Lithium4Fit{args.model}.root', 'RECREATE')

    with PdfCanvas(f'output/{args.sign}Lithium4Fit{args.model}.pdf') as pdf_canvas:
        for centrality_opt in ['integrated', '0-10', '10-30', '30-50']:

            cent = SUFFIX_DICT['undot'][centrality_opt]
            cent_dot = SUFFIX_DICT['dot'][centrality_opt]
            outdir = outfile.mkdir(f'dir{cent}')

            roows = RooWorkspace(f'ws{cent}')
            kstar = RooRealVar('kstar', '#it{k}* (GeV/#it{c})', 0.02, 0.8)
            getattr(roows, 'import')(kstar)

            # signal template - shared in all centralities
            h_signal = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/MC/data_visual_selectionsPr.root',
                                                f'Correlations/fKstar{args.sign}'))
            h_mixed_event = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
                                                    f'Correlations/fKstar{args.sign}'))
            prepare_signal_fit(roows, h_signal, h_mixed_event, outdir)

            # coulomb template
            h_coulomb = load_coulomb_template(centrality_opt, args.model)
            prepare_coulomb_template(roows, h_coulomb, outdir)

            # p-He3 correlation function
            h_CF = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/PbPb/studies.root', 
                                          f'Correlation{args.sign}/hCorrelation_kstar{cent_dot}'))
            pull_from_coulomb(h_CF, h_coulomb, outdir)
            fit_CF(roows, h_CF, outdir, pdf_canvas)
            evaluate_pvalue_and_significance(roows, outdir)

            del roows

    outfile.Close()
