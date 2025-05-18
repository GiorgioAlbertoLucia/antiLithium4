import numpy as np
import pandas as pd
from ROOT import TFile, TGraphErrors, TCanvas, TPaveText, TF1, TH2F
from ROOT import RooRealVar, RooCrystalBall, RooGenericPdf, RooAddPdf, RooDataHist, RooArgList, RooChebychev, RooGaussian
from torchic import HistLoadInfo
from torchic import histogram

from utils.root_setter import obj_setter
from utils.utils import write_params_to_text, create_graph, PdfCanvas

def build_signal_model(x):
    # signal function
    signal_pars = {
        'mean': RooRealVar('mean', '#mu', 0., -1, 1),
        'sigma': RooRealVar('sigma', '#sigma', 1, 0.3, 2),
        'aL': RooRealVar('aL', '#alpha', 1., 0.7, 10.),
        'nL': RooRealVar('nL', 'n', 1., 0.1, 30.)
    }
    signal = RooCrystalBall('signal', 'signal', x, 
                            signal_pars['mean'], signal_pars['sigma'], 
                            signal_pars['aL'], signal_pars['nL'], doubleSided=True)
    return signal, signal_pars

def build_background_model(x):
    # background function
    bkg_exp_pars = {
        'bkg_exp_alpha': RooRealVar('alpha', 'alpha', 4., 1., 10.),
        'bkg_exp_offset': RooRealVar('offset', 'offset', 0.5, 0., 10.)
    }
    bkg_exp = RooGenericPdf('bkg_exp', 'bkg_exp', 'exp(-alpha * (x - offset))', 
                            [x,
                             bkg_exp_pars['bkg_exp_alpha'],
                             bkg_exp_pars['bkg_exp_offset']])
    bkg_gaus_pars = {
        'bkg_gaus_mean': RooRealVar('mean_gaus', 'mean_gaus', 0., -6, -2),
        'bkg_gaus_sigma': RooRealVar('sigma_gaus', 'sigma_gaus', 0.5, 1e-3, 1e3)
    }
    bkg_gaus = RooGaussian('bkg_gaus', 'bkg_gaus', x, 
                            bkg_gaus_pars['bkg_gaus_mean'], bkg_gaus_pars['bkg_gaus_sigma'])
    
    bkg_const = RooGenericPdf('bkg_const', 'bkg_const', '1', [x])
    bkg_pol1_pars = {
        'bkg_pol1_c0': RooRealVar('c0', 'c_{0}', 1., -10., 10.),
        'bkg_pol1_c1': RooRealVar('c1', 'c_{1}', 1., -10., 10.),
    }
    bkg_pol1 = RooGenericPdf('bkg_pol1', 'bkg_pol1', f'c0 + c1 * {x.GetName()}', [x, 
                                                                     bkg_pol1_pars['bkg_pol1_c0'],
                                                                     bkg_pol1_pars['bkg_pol1_c1']])
    return bkg_const, bkg_pol1, bkg_pol1_pars

def purity_fit_slice(model, signal_pdf, bkg_pdf, sig_frac, hist, x, signal_pars, pt_low_edge, pt_high_edge, nsigma_low_int, nsigma_high_int):

    pt = 0.5 * (pt_high_edge + pt_low_edge)
    dh = RooDataHist(f'dh_{pt:.2f}', 'dh', [x], Import=hist)
    model.fitTo(dh, PrintLevel=-1)

    frame = x.frame(Title=f'{pt_low_edge:.2f} < #it{{p}}_{{T}} < {pt_high_edge:.2f} GeV/#it{{c}}')
    dh.plotOn(frame)
    model.plotOn(frame, LineColor=2)
    model.plotOn(frame, Components={signal_pdf}, LineColor=3, LineStyle='--')
    model.plotOn(frame, Components={bkg_pdf}, LineColor=4, LineStyle='--')
    #for ipdf, bkg_pdf in enumerate(bkg_pdfs, start=4):
    #    model.plotOn(frame, Components={bkg_pdf}, LineColor=ipdf, LineStyle='--')

    x.setRange('integral_range', nsigma_low_int, nsigma_high_int)
    integral_signal = signal_pdf.createIntegral(x, x, 'integral_range').getVal() * sig_frac.getVal()
    integral_signal_unnorm = integral_signal * hist.Integral(hist.FindBin(nsigma_low_int), hist.FindBin(nsigma_high_int))
    
    #integral_bkg = 0
    #for bkg_pdf in bkg_pdfs:
    #    integral_bkg += bkg_pdf.createIntegral(x, x, 'integral_range').getVal() * sig_frac.getVal() ## WIP
    integral_bkg = bkg_pdf.createIntegral(x, x, 'integral_range').getVal() * (1 - sig_frac.getVal())

    ifit_results = {'pt': np.abs(pt), 'int_sig': integral_signal, 'sig_mean': signal_pars['mean'].getVal(), 
                   'sig_sigma': signal_pars['sigma'].getVal(), 
                   'sig_mean_err': signal_pars['sigma'].getVal()/np.sqrt(integral_signal_unnorm),
                   'sig_sigma_err': signal_pars['sigma'].getError(),
                   'purity': integral_signal / (integral_signal + integral_bkg)}

    return frame, ifit_results

def draw_results(h2_nsigmaTPC: TH2F, fit_results: pd.DataFrame, pt_min, pt_max, pdf_canvas: PdfCanvas, out_dir: TFile, sign: str):
    '''
        Draw the results of the fit
    '''
    fit_results['pt_err'] = h2_nsigmaTPC.GetXaxis().GetBinWidth(1)/2.
    graph_signal = create_graph(fit_results, 'pt', 'int_sig', 'pt_err', 0., name=f'g_signal_{sign}', title=';#it{p}_{T} (GeV/c); int (n#sigma_{TPC})')

    graph_signal_mean = create_graph(fit_results, 'pt', 'sig_mean', 'pt_err', 0., name=f'g_signal_mean_{sign}', title=';#it{p}_{T} (GeV/c); #mu(n#sigma_{TPC})')
    fit_mean = TF1('fit_mean', '[0]', pt_min, pt_max)
    graph_signal_mean.Fit(fit_mean, 'RSM+')
    graph_signal_mean.SetMarkerStyle(20)
    pdf_canvas.draw_object(graph_signal_mean, logy=False, draw_option='ap')

    graph_signal_sigma = create_graph(fit_results, 'pt', 'sig_sigma', 'pt_err', 0., name=f'g_signal_sigma_{sign}', title=';#it{p}_{T} (GeV/c); #sigma(n#sigma_{TPC})')
    fit_sigma = TF1('fit_sigma', 'pol0', pt_min, pt_max)
    graph_signal_sigma.Fit(fit_sigma, 'RSM+')
    graph_signal_sigma.SetMarkerStyle(20)
    pdf_canvas.draw_object(graph_signal_sigma, logy=False, draw_option='ap')

    graph_purity = create_graph(fit_results, 'pt', 'purity', 'pt_err', 0., name=f'g_purity_{sign}', title=';#it{p}_{T} (GeV/c); purity')

    out_dir.cd()
    graph_signal.Write(f'g_signal_{sign}')
    graph_signal_mean.Write(f'g_signal_mean_{sign}')
    graph_signal_sigma.Write(f'g_signal_sigma_{sign}')
    graph_purity.Write(f'g_purity_{sign}')

def purity_he3_TPC(h2_nsigmaTPC: TH2F, output_file: TFile):

    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)
    
    pt_min, pt_max = 1.6, 3.5
    nsigma_min, nsigma_max = -3, 3
    nsigma_low_int, nsigma_high_int = -2, 2
    
    nsigma = RooRealVar('nsigma', 'n#sigma_{TPC}', nsigma_min, nsigma_max)

    signal, signal_pars = build_signal_model(nsigma)
    bkg_const, bkg_pol1, bkg_pol1_pars = build_background_model(nsigma)

    sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
    #model_anti = RooAddPdf('model_anti', 'signal + bkg_const', [signal, bkg_const], [sig_frac])
    model_anti = RooAddPdf('model_anti', 'signal + bkg_pol1', [signal, bkg_pol1], [sig_frac])

    out_dir = output_file.mkdir('He_TPC')

    for sign in ['matter', 'antimatter']:

        output_pdf_path = output_file.GetFile().GetName().split('.')[0] + '_He.pdf'
        output_pdf_path = output_pdf_path.replace('.pdf', f'_{sign}.pdf')
        with PdfCanvas(output_pdf_path) as pdf_canvas:
        
            fit_results = None

            if sign == 'antimatter':
                pt_min, pt_max = -3.5, -1.6
                model = model_anti
            else:
                pt_min, pt_max = 1.6, 3.5
                model = model_anti

            for pt_bin in range(h2_nsigmaTPC.GetXaxis().FindBin(pt_min), h2_nsigmaTPC.GetXaxis().FindBin(pt_max)):

                pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)
                h_nsigmaTPC = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt:.2f}', pt_bin, pt_bin, 'e')
                h_nsigmaTPC.Rebin(4)

                pt_low_edge = h2_nsigmaTPC.GetXaxis().GetBinLowEdge(pt_bin)
                pt_up_edge = h2_nsigmaTPC.GetXaxis().GetBinLowEdge(pt_bin+1)
                nsigma_frame, ifit_results = purity_fit_slice(model, signal, bkg_pol1, sig_frac, h_nsigmaTPC, 
                                                              nsigma, signal_pars, pt_low_edge, pt_up_edge, 
                                                              nsigma_low_int, nsigma_high_int)
                
                if fit_results is None:
                    fit_results = pd.DataFrame.from_dict([ifit_results])
                else:
                    fit_results = pd.concat([fit_results, pd.DataFrame.from_dict([ifit_results])], ignore_index=True)

                text = write_params_to_text({**signal_pars, **bkg_pol1_pars}.values(), coordinates=[0.4, 0.45, 0.63, 0.75])
                #text.AddText(f'#chi^{{2}} / NDF = {nsigma_frame.chiSquare():.2f}')
                nsigma_frame.addObject(text)
                pdf_canvas.draw_object(nsigma_frame, logy=True)

                out_dir.cd()
                nsigma_frame.Write(f'nsigma_frame_{sign}_{pt:.2f}')

            draw_results(h2_nsigmaTPC, fit_results, pt_min, pt_max, pdf_canvas, out_dir, sign)
                
    h2_nsigmaTPC.Write()

if __name__ == '__main__':
    
    output_file = TFile('output/LHC24PbPb/purity_new.root', 'recreate')
    #output_file = TFile('output/LHC23PbPb/purity.root', 'recreate')

    #infile_path_TPC = '/Users/glucia/Projects/ALICE/data/lithium/same/AnalysisResults_LHC24ar_pass1_same.root'
    #infile_path_TPC = '/data/galucia/lithium_local/same/AnalysisResults_LHC24ar_pass1_same.root'
    #infile_path_TOF = 'output/LHC24PbPb/qa_purity.root'
    infile_path_TPC = 'output/LHC24PbPb/qa_selectionsPr.root'
    #infile_path_TOF = 'output/LHC23PbPb/qa_purity.root'

    h2_nsigmaTPC_info = HistLoadInfo(infile_path_TPC, 'TPC/NSigmaTPCvsPtHe3')
    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)

    purity_he3_TPC(h2_nsigmaTPC, output_file)
    output_file.Close()