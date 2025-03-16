import numpy as np
from ROOT import TFile, TGraphErrors, TCanvas, TPaveText
from ROOT import RooRealVar, RooCrystalBall, RooGenericPdf, RooAddPdf, RooDataHist, RooArgList, RooChebychev, RooGaussian
from torchic import HistLoadInfo
from torchic import histogram
from torchic.core.roofitter import Roofitter

def purity_proton_TPC(infile_path: str, output_file: TFile):

    h2_nsigmaTPC_info = HistLoadInfo(infile_path, 'he3hadronfemto/QA/h2NsigmaHadronTPC_preselection')
    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)
    
    pt_min = 0.2
    pt_max = 0.8
    nsigma_min = -8.
    nsigma_max = 8.

    nsigma = RooRealVar('nsigma', 'n#sigma_{TPC}', nsigma_min, nsigma_max)

    # signal function
    signal_pars = {
        'mean': RooRealVar('mean', 'mean', 0., -2.5, 2.5),
        'sigma': RooRealVar('sigma', 'sigma', 0.5, 1e-3, 1e3),
        'aL': RooRealVar('aL', 'aL', 1., 0.7, 30.),
        'nL': RooRealVar('nL', 'nL', 1., 0.1, 30.),
        'aR': RooRealVar('aR', 'aR', 1., 0.7, 30.),
        'nR': RooRealVar('nR', 'nR', 1., 0.1, 30.)
    }
    signal = RooCrystalBall('signal', 'signal', nsigma, 
                            signal_pars['mean'], signal_pars['sigma'], 
                            signal_pars['aL'], signal_pars['nL'], signal_pars['aR'], signal_pars['nR'])
    bkg_exp_pars = {
        'alpha': RooRealVar('alpha', 'alpha', 0.5, 0., 10.),
        'offset': RooRealVar('offset', 'offset', 0.5, 0., 10.)
    }
    bkg_exp = RooGenericPdf('bkg_exp', 'bkg_exp', 'exp(-alpha * (nsigma - offset))', 
                            RooArgList(nsigma,
                             bkg_exp_pars['alpha'],
                             bkg_exp_pars['offset']))

    sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
    model = RooAddPdf('model', 'signal + bkg_frac', [signal, bkg_exp], [sig_frac])

    out_dir = output_file.mkdir('Pr_TPC')

    for sign in ['matter', 'antimatter']:

        pts = []
        sig = []
        bkg = []
        tot = []
        purity = []

        if sign == 'antimatter':
            pt_min = -0.8
            pt_max = -0.2
        else:
            pt_min = 0.2
            pt_max = 0.8

        for pt_bin in range(h2_nsigmaTPC.GetXaxis().FindBin(pt_min), h2_nsigmaTPC.GetXaxis().FindBin(pt_max)):

            pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)
            h_nsigmaTPC = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt:.2f}', pt_bin, pt_bin, 'e')
            h_nsigmaTPC.Rebin(4)

            dh = RooDataHist('dh', f'dh_{pt:.2f}', [nsigma], Import=h_nsigmaTPC)
            model.fitTo(dh, PrintLevel=-1)

            pt_low_edge = h2_nsigmaTPC.GetXaxis().GetBinLowEdge(pt_bin)
            pt_up_edge = h2_nsigmaTPC.GetXaxis().GetBinLowEdge(pt_bin+1)
            nsigma_frame = nsigma.frame(Title=f'{pt_low_edge:.2f} < p_{{T}} < {pt_up_edge:.2f} GeV/#it{{c}}')
            dh.plotOn(nsigma_frame)
            model.plotOn(nsigma_frame, LineColor=2)
            model.paramOn(nsigma_frame)
            model.plotOn(nsigma_frame, Components={signal}, LineColor=3, LineStyle='--')
            model.plotOn(nsigma_frame, Components={bkg_exp}, LineColor=4, LineStyle='--')

            nsigma.setRange('integral_range', -2, 2)
            integral_signal = signal.createIntegral(nsigma, nsigma, 'integral_range').getVal() * sig_frac.getVal()
            integral_bkg = bkg_exp.createIntegral(nsigma, nsigma, 'integral_range').getVal() * (1 - sig_frac.getVal())

            pts.append(np.abs(pt))
            sig.append(integral_signal)
            bkg.append(integral_signal)
            tot.append(integral_signal + integral_bkg)
            purity.append(integral_signal / (integral_signal + integral_bkg))

            canvas = TCanvas(f'cNSigmaTPC_{pt:.2f}', f'cNSigmaTPC_{pt_bin}', 800, 600)
            nsigma_frame.Draw()
            out_dir.cd()
            canvas.Write()

        pt_bin_width = h2_nsigmaTPC.GetXaxis().GetBinWidth(1)/2.
        graph_signal = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_bkg = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_tot = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(tot, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_purity = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(purity, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))

        out_dir.cd()
        graph_signal.Write(f'g_signal_{sign}')
        graph_bkg.Write(f'g_bkg_{sign}')
        graph_tot.Write(f'g_tot_{sign}')
        graph_purity.Write(f'g_purity_{sign}')
    h2_nsigmaTPC.Write()

def purity_proton_TOF(infile_path: str, output_file: TFile):

    h2_nsigmaTOF_info = HistLoadInfo(infile_path, 'TOF/NSigmaTOFvsPtHad')
    h2_nsigmaTOF = histogram.load_hist(h2_nsigmaTOF_info)
    output_pdf_path = output_file.GetFile().GetName()
    output_pdf_path = output_pdf_path.split('.')[0] + '_Pr_TOF.pdf'

    pt_min = 0.7
    pt_max = 2
    nsigma_min = -8
    nsigma_max = 8

    nsigma = RooRealVar('nsigma', 'n#sigma_{TOF}', nsigma_min, nsigma_max)

    # signal function
    signal_pars = {
        'mean': RooRealVar('mean', 'mean', 0., -2.5, 2.5),
        'sigma': RooRealVar('sigma', 'sigma', 0.5, 1e-3, 1e3),
        'aL': RooRealVar('aL', 'aL', 1., 0.7, 10.),
        'nL': RooRealVar('nL', 'nL', 1., 2., 40.),
        'aR': RooRealVar('aR', 'aR', 1., 0.7, 10.),
        'nR': RooRealVar('nR', 'nR', 1., 2., 30.)
    }
    signal = RooCrystalBall('signal', 'signal', nsigma, 
                            signal_pars['mean'], signal_pars['sigma'], 
                            signal_pars['aL'], signal_pars['nL'], signal_pars['aR'], signal_pars['nR'])
    #bkg_exp_pars = {
    #    'alpha': RooRealVar('alpha', 'alpha', 0.5, 0., 10.),
    #    'offset': RooRealVar('offset', 'offset', 0.5, 0., 10.)
    #}
    #bkg_exp = RooGenericPdf('bkg_exp', 'bkg_exp', 'exp(-alpha * (nsigma - offset))', 
    #                        [nsigma,
    #                         bkg_exp_pars['alpha'],
    #                         bkg_exp_pars['offset']])
    pol_pars = {
        'c0': RooRealVar('c0', 'c0', 1., 0., 10.),
        'c1': RooRealVar('c1', 'c1', 1., -10., 10.),
    }
    bkg_pol = RooGenericPdf('bkg_pol', 'bkg_pol', 'c0 + c1 * nsigma', 
                            [nsigma,
                             pol_pars['c0'],
                             pol_pars['c1']])
    sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
    #model = RooAddPdf('model', 'signal + bkg_frac', [signal, bkg_exp], [sig_frac])
    model = RooAddPdf('model', 'signal + bkg_frac', [signal, bkg_pol], [sig_frac])

    out_dir = output_file.mkdir('Pr_TOF')

    for sign in ['antimatter']: #['matter', 'antimatter']:

        pts = []
        sig = []
        sig_mean = []
        sig_mean_err = []
        sig_sigma = []
        sig_sigma_err = []
        bkg = []
        tot = []
        purity = []

        if sign == 'antimatter':
            pt_min = -2
            pt_max = -0.7
        else:
            pt_min = 0.7
            pt_max = 2

        pt_bin_min = h2_nsigmaTOF.GetXaxis().FindBin(pt_min)
        pt_bin_max = h2_nsigmaTOF.GetXaxis().FindBin(pt_max)
        for ibin, pt_bin in enumerate(range(pt_bin_min, pt_bin_max)):

            pt = h2_nsigmaTOF.GetXaxis().GetBinCenter(pt_bin)
            h_nsigmaTOF = h2_nsigmaTOF.ProjectionY(f'h_nsigmaTOF_{pt:.2f}', pt_bin, pt_bin, 'e')
            h_nsigmaTOF.Rebin(4)

            dh = RooDataHist(f'dh_{pt:.2f}', 'dh', [nsigma], Import=h_nsigmaTOF)
            model.fitTo(dh, PrintLevel=-1)

            pt_low_edge = h2_nsigmaTOF.GetXaxis().GetBinLowEdge(pt_bin)
            pt_up_edge = h2_nsigmaTOF.GetXaxis().GetBinLowEdge(pt_bin+1)
            nsigma_frame = nsigma.frame(Title=f'{pt_low_edge:.2f} < p_{{T}} < {pt_up_edge:.2f} GeV/#it{{c}}')
            dh.plotOn(nsigma_frame)
            model.plotOn(nsigma_frame, LineColor=2)
            #model.paramOn(nsigma_frame)
            model.plotOn(nsigma_frame, Components={signal}, LineColor=3, LineStyle='--')
            #model.plotOn(nsigma_frame, Components={bkg_exp}, LineColor=4, LineStyle='--')
            model.plotOn(nsigma_frame, Components={bkg_pol}, LineColor=4, LineStyle='--')

            nsigma.setRange('integral_range', -2, 2)
            integral_signal = signal.createIntegral(nsigma, nsigma, 'integral_range').getVal() * sig_frac.getVal()
            #integral_bkg = bkg_exp.createIntegral(nsigma, nsigma, 'integral_range').getVal() * (1 - sig_frac.getVal())
            integral_bkg = bkg_pol.createIntegral(nsigma, nsigma, 'integral_range').getVal() * (1 - sig_frac.getVal())

            pts.append(np.abs(pt))
            sig.append(integral_signal)
            sig_mean.append(signal_pars['mean'].getVal())
            sig_mean_err.append(signal_pars['sigma'].getVal()/np.sqrt(integral_signal / (integral_signal + integral_bkg) * h_nsigmaTOF.Integral(
                h_nsigmaTOF.FindBin(-2), h_nsigmaTOF.FindBin(2)
            ) ))
            sig_sigma.append(signal_pars['sigma'].getVal())
            sig_sigma_err.append(signal_pars['sigma'].getError())
            bkg.append(integral_signal)
            tot.append(integral_signal + integral_bkg)
            purity.append(integral_signal / (integral_signal + integral_bkg))

            canvas = TCanvas(f'cNSigmaTPC_{pt:.2f}', f'cNSigmaTPC_{pt_bin}', 800, 600)
            nsigma_frame.Draw()
            text = TPaveText(0.6, 0.55, 0.85, 0.85, 'NDC')
            for param in signal_pars.values():
                text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
            for param in pol_pars.values():
                text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
            #text.AddText(f'#chi^{{2}} / NDF = {nsigma_frame.chiSquare():.2f}')
            text.SetFillColor(0)
            text.SetBorderSize(0)
            text.Draw()
            out_dir.cd()
            canvas.SetLogy()
            canvas.Write()

            path, extension = output_pdf_path.split('.')
            if ibin == 0:
                canvas.Print(f'{path}_{sign}.{extension}(')
            elif ibin == pt_bin_max - pt_bin_min - 1:
                canvas.Print(f'{path}_{sign}.{extension})')
            else:
                canvas.Print(f'{path}_{sign}.{extension}')


        pt_bin_width = h2_nsigmaTOF.GetXaxis().GetBinWidth(1)/2.
        graph_signal = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_signal_mean = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig_mean, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array(sig_mean_err, dtype=np.float32))
        graph_signal_sigma = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig_sigma, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array(sig_sigma_err, dtype=np.float32))
        graph_bkg = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_tot = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(tot, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_purity = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(purity, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))

        out_dir.cd()
        graph_signal.Write(f'g_signal_{sign}')
        graph_signal_mean.Write(f'g_signal_mean_{sign}')
        graph_signal_sigma.Write(f'g_signal_sigma_{sign}')
        graph_bkg.Write(f'g_bkg_{sign}')
        graph_tot.Write(f'g_tot_{sign}')
        graph_purity.Write(f'g_purity_{sign}')
    h2_nsigmaTOF.Write()

def purity_he3_TPC(infile_path: str, output_file: TFile):

    h2_nsigmaTPC_info = HistLoadInfo(infile_path, 'he3hadronfemto/QA/h2NsigmaHe3TPC_preselection')
    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)
    output_pdf_path = output_file.GetFile().GetName()
    output_pdf_path = output_pdf_path.split('.')[0] + '_He.pdf'
    
    pt_min = 1.6
    pt_max = 3.5
    nsigma_min = -6
    nsigma_max = 6
    nsigma_low_int = -2
    nsigma_high_int = 2

    # 2024
    # nsigma_mean = 0.35
    # nsigma_sigma = 0.7
    #nsigma_low_int = -1.05
    #nsigma_high_int = 1.75

    # 2023
    # nsigma_mean = 0.5
    # nsigma_sigma = 0.7
    #nsigma_low_int = -0.90
    #nsigma_high_int = 1.9 
    
    nsigma = RooRealVar('nsigma', 'n#sigma_{TPC}', nsigma_min, nsigma_max)

    # signal function
    signal_pars = {
        'mean': RooRealVar('mean', 'mean', 0., -0.5, 2.5),
        'sigma': RooRealVar('sigma', 'sigma', 1, 0.3, 2),
        'aL': RooRealVar('aL', 'aL', 1., 0.7, 10.),
        'nL': RooRealVar('nL', 'nL', 1., 0.1, 30.),
        #'aR': RooRealVar('aR', 'aR', 1., 0.1, 10.),
        #'nR': RooRealVar('nR', 'nR', 1., 0.1, 30.)
    }
    #signal = RooCrystalBall('signal', 'signal', nsigma, 
    #                        signal_pars['mean'], signal_pars['sigma'], 
    #                        signal_pars['aL'], signal_pars['nL'], signal_pars['aR'], signal_pars['nR'])
    signal = RooCrystalBall('signal', 'signal', nsigma, 
                            signal_pars['mean'], signal_pars['sigma'], 
                            signal_pars['aL'], signal_pars['nL'], doubleSided=True)
    bkg_exp_pars = {
        'alpha': RooRealVar('alpha', 'alpha', 4., 0., 10.),
        'offset': RooRealVar('offset', 'offset', 0.5, 0., 10.)
    }
    bkg_exp = RooGenericPdf('bkg_exp', 'bkg_exp', 'exp(-alpha * (nsigma - offset))', 
                            [nsigma,
                             bkg_exp_pars['alpha'],
                             bkg_exp_pars['offset']])
    bkg_gaus_pars = {
        'mean': RooRealVar('mean_gaus', 'mean_gaus', 0., -6, -2),
        'sigma': RooRealVar('sigma_gaus', 'sigma_gaus', 0.5, 1e-3, 1e3)
    }
    bkg_gaus = RooGaussian('bkg_gaus', 'bkg_gaus', nsigma, 
                            bkg_gaus_pars['mean'], bkg_gaus_pars['sigma'])
    sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
    gaus_frac = RooRealVar('gaus_frac', 'gaus_frac', 0.5, 0., 1.)
    model_anti = RooAddPdf('model_anti', 'signal + bkg_exp', [signal, bkg_exp], [sig_frac])
    model_matter = RooAddPdf('model_matter', 'signal + bkg_exp + bkg_gaus', [signal, bkg_gaus, bkg_exp], [sig_frac, gaus_frac])

    out_dir = output_file.mkdir('He_TPC')

    for sign in ['matter', 'antimatter']:

        pts = []
        sig = []
        sig_mean = []
        sig_mean_err = []
        sig_sigma = []
        sig_sigma_err = []
        bkg = []
        tot = []
        purity = []

        if sign == 'antimatter':
            pt_min = -3.5
            pt_max = -1.6
            model = model_anti
        else:
            pt_min = 1.6
            pt_max = 3.5
            model = model_matter

        pt_bin_min = h2_nsigmaTPC.GetXaxis().FindBin(pt_min)
        pt_bin_max = h2_nsigmaTPC.GetXaxis().FindBin(pt_max)

        for ibin, pt_bin in enumerate(range(h2_nsigmaTPC.GetXaxis().FindBin(pt_min), h2_nsigmaTPC.GetXaxis().FindBin(pt_max))):

            pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)
            if pt > 2.6:
                model = model_anti # no more H3 bkg
            h_nsigmaTPC = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt:.2f}', pt_bin, pt_bin, 'e')
            h_nsigmaTPC.Rebin(4)

            dh = RooDataHist(f'dh_{pt:.2f}', 'dh', [nsigma], Import=h_nsigmaTPC)
            model.fitTo(dh, PrintLevel=-1)

            pt_low_edge = h2_nsigmaTPC.GetXaxis().GetBinLowEdge(pt_bin)
            pt_up_edge = h2_nsigmaTPC.GetXaxis().GetBinLowEdge(pt_bin+1)
            nsigma_frame = nsigma.frame(Title=f'{pt_low_edge:.2f} < p_{{T}} < {pt_up_edge:.2f} GeV/#it{{c}}')
            dh.plotOn(nsigma_frame)
            model.plotOn(nsigma_frame, LineColor=2)
            #model.paramOn(nsigma_frame)
            model.plotOn(nsigma_frame, Components={signal}, LineColor=3, LineStyle='--')
            model.plotOn(nsigma_frame, Components={bkg_exp}, LineColor=4, LineStyle='--')
            model.plotOn(nsigma_frame, Components={bkg_gaus}, LineColor=5, LineStyle='--')

            nsigma.setRange('integral_range', nsigma_low_int, nsigma_high_int)
            integral_signal = signal.createIntegral(nsigma, nsigma, 'integral_range').getVal() * sig_frac.getVal()
            integral_bkg_exp = bkg_exp.createIntegral(nsigma, nsigma, 'integral_range').getVal() * (1 - sig_frac.getVal() - gaus_frac.getVal())
            integral_bkg_gaus = bkg_gaus.createIntegral(nsigma, nsigma, 'integral_range').getVal() * gaus_frac.getVal()
            integral_bkg = integral_bkg_exp + integral_bkg_gaus

            pts.append(np.abs(pt))
            sig.append(integral_signal)
            sig_mean.append(signal_pars['mean'].getVal())
            sig_mean_err.append(signal_pars['sigma'].getVal()/np.sqrt(integral_signal / (integral_signal + integral_bkg) * h_nsigmaTPC.Integral(
                h_nsigmaTPC.FindBin(nsigma_low_int), h_nsigmaTPC.FindBin(nsigma_high_int)
            ) ))
            sig_sigma.append(signal_pars['sigma'].getVal())
            sig_sigma_err.append(signal_pars['sigma'].getError())
            bkg.append(integral_signal)
            tot.append(integral_signal + integral_bkg)
            purity.append(integral_signal / (integral_signal + integral_bkg))

            canvas = TCanvas(f'cNSigmaTPC_{pt:.2f}', f'cNSigmaTPC_{pt:.2f}', 800, 600)
            nsigma_frame.Draw()
            text = TPaveText(0.6, 0.55, 0.85, 0.85, 'NDC')
            for param in signal_pars.values():
                text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
            for param in bkg_gaus_pars.values():
                text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
            for param in bkg_exp_pars.values():
                text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
            #text.AddText(f'#chi^{{2}} / NDF = {nsigma_frame.chiSquare():.2f}')
            text.SetFillColor(0)
            text.SetBorderSize(0)
            text.Draw()
            out_dir.cd()
            canvas.SetLogy()
            canvas.Write()

            path, extension = output_pdf_path.split('.')
            if ibin == 0:
                canvas.Print(f'{path}_{sign}.{extension}(')
            elif ibin == pt_bin_max - pt_bin_min - 1:
                canvas.Print(f'{path}_{sign}.{extension})')
            else:
                canvas.Print(f'{path}_{sign}.{extension}')

        pt_bin_width = h2_nsigmaTPC.GetXaxis().GetBinWidth(1)/2.
        graph_signal = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_signal_mean = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig_mean, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array(sig_mean_err, dtype=np.float32))
        graph_signal_sigma = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(sig_sigma, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array(sig_sigma_err, dtype=np.float32))
        graph_bkg = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_tot = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(tot, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
        graph_purity = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(purity, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))

        out_dir.cd()
        graph_signal.Write(f'g_signal_{sign}')
        graph_signal_mean.Write(f'g_signal_mean_{sign}')
        graph_signal_sigma.Write(f'g_signal_sigma_{sign}')
        graph_bkg.Write(f'g_bkg_{sign}')
        graph_tot.Write(f'g_tot_{sign}')
        graph_purity.Write(f'g_purity_{sign}')
    h2_nsigmaTPC.Write()

if __name__ == '__main__':
    
    output_file = TFile('output/LHC24PbPb/purity.root', 'recreate')
    #output_file = TFile('output/LHC23PbPb/purity.root', 'recreate')

    #infile_path_TPC = '/Users/glucia/Projects/ALICE/data/lithium/same/AnalysisResults_LHC24ar_pass1_same.root'
    infile_path_TPC = '/data/galucia/lithium_local/same/AnalysisResults_LHC24ar_pass1_same.root'
    infile_path_TOF = 'output/LHC24PbPb/qa_purity.root'
    #infile_path_TPC = '/data/galucia/lithium_local/same/AnalysisResults_LHC23_PbPb_pass4_same.root'
    #infile_path_TOF = 'output/LHC23PbPb/qa_purity.root'

    purity_proton_TPC(infile_path_TPC, output_file)
    purity_proton_TOF(infile_path_TOF, output_file)
    purity_he3_TPC(infile_path_TPC, output_file)
    output_file.Close()