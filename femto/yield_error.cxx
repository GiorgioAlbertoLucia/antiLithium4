/**
 * This code computes the chi2 disitrbution for the null hypothesis 
*/

#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>

#include <RooRealVar.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>

#include "crystal_ball.hh"

typedef struct CrystalBallPars {
    RooRealVar mean = RooRealVar("mean", "mean", 0.081, 0.06, 0.1);
    RooRealVar sigma = RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1);
    RooRealVar aL = RooRealVar("aL", "aL", 1.3, 0.1, 10.0);
    RooRealVar nL = RooRealVar("nL", "nL", 2.7, 0.1, 10.0);
    RooRealVar aR = RooRealVar("aR", "aR", 1.1, 0.1, 10.0);
    RooRealVar nR = RooRealVar("nR", "nR", 5.7, 0.1, 10.0);

    void print() {
        std::cout << "Mean: " << mean.getVal() << " +/- " << mean.getError() << std::endl;
        std::cout << "Sigma: " << sigma.getVal() << " +/- " << sigma.getError() << std::endl;
        std::cout << "aL: " << aL.getVal() << " +/- " << aL.getError() << std::endl;
        std::cout << "nL: " << nL.getVal() << " +/- " << nL.getError() << std::endl;
        std::cout << "aR: " << aR.getVal() << " +/- " << aR.getError() << std::endl;
        std::cout << "nR: " << nR.getVal() << " +/- " << nR.getError() << std::endl;
    }

    void setConstantTail() {
        aL.setConstant(true);
        nL.setConstant(true);
        aR.setConstant(true);
        nR.setConstant(true);
    }

} CrystalBallPars;


void load_same_mixed(TH1F *& h_same, TH1F *& h_mixed) {

    TFile *file = TFile::Open("/home/galucia/antiLithium4/analysis/output/PbPb/studies.root");

    //const char * sign = "Matter";
    const char * sign = "Anti";
    
    h_same = (TH1F*)file->Get(Form("Correlation%s/hSame_kstar", sign));
    h_same->SetDirectory(0);

    h_mixed = (TH1F*)file->Get(Form("Correlation%s/hMixed_kstar", sign));
    h_mixed->SetDirectory(0);

    file->Close();
}

void compute_correlation_function(TH1F *h_same, TH1F *h_mixed, TH1F* h_correlation) {

    for  (int ibin = 1; ibin <= h_same->GetNbinsX()+1; ++ibin) {
        double same_value = h_same->GetBinContent(ibin);
        double mixed_value = h_mixed->GetBinContent(ibin);
        double same_error = std::sqrt(same_value);
        double mixed_error = std::sqrt(mixed_value);

        if (mixed_value > 0) {
            double correlation_value = same_value / mixed_value;
            h_correlation->SetBinContent(ibin, correlation_value);
            h_correlation->SetBinError(ibin, correlation_value*std::sqrt((same_error/same_value)*(same_error/same_value) + (mixed_error/mixed_value)*(mixed_error/mixed_value)) );
        } else {
            h_correlation->SetBinContent(ibin, 1e-12);
            h_correlation->SetBinError(ibin, 1e-12);
        }
    }
}

void prepare_signal_fit(CrystalBallPars &pars, TFile *outfile) {

    using namespace RooFit;

    TFile *file_mc = TFile::Open("/home/galucia/antiLithium4/femto/output/li4_intrinsic_width.root");
    auto h_signal = (TH1F*)file_mc->Get("sampling/h_sample_kstar");
    h_signal->SetDirectory(0);
    file_mc->Close();

    TFile *file_mixed = TFile::Open("/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root");
    auto h_mixed = (TH1F*)file_mixed->Get("Correlations/fKstarAnti");
    h_mixed->SetDirectory(0);
    file_mixed->Close();

    TH1F* h_correlation_signal = (TH1F*)h_signal->Clone("h_correlation_signal");
    compute_correlation_function(h_signal, h_mixed, h_correlation_signal);

    RooRealVar kstar("kstar", "kstar", 0.02, 0.25);
    RooCrystalBall signal("signal", "signal", kstar, pars.mean, pars.sigma, pars.aL, pars.nL, pars.aR, pars.nR);
    RooDataHist data("data", "data", kstar, Import(*h_correlation_signal));

    RooFitResult* result = signal.fitTo(data, PrintLevel(-1), Save(), SumW2Error(true));
    std::cout << "Crystal Ball parameters:" << std::endl;
    pars.print();

    const double mean0 = pars.mean.getVal();
    const double mean0_error = pars.mean.getError();
    const double sigma0 = pars.sigma.getVal();
    const double sigma0_error = pars.sigma.getError();

    pars.mean.setRange(mean0 - 3 * mean0_error, mean0 + 3 * mean0_error);
    //pars.sigma.setRange(0.5 * sigma0, 3 * sigma0);
    pars.sigma.setRange(sigma0 - 3 * sigma0_error, sigma0 + 3 * sigma0_error);

    TCanvas canvas("canvas_fit", "Signal Fit", 800, 600);
    RooPlot* frame = kstar.frame();
    data.plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    signal.plotOn(frame, LineColor(kRed), LineWidth(2));
    signal.paramOn(frame, Layout(0.6, 0.95, 0.95), Format("NEU", AutoPrecision(2)));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
    canvas.cd();
    frame->Draw();

    outfile->mkdir("signal_fit");
    outfile->cd("signal_fit");
    h_signal->Write("h_signal");
    h_mixed->Write("h_mixed");
    canvas.Write();
    
    delete h_correlation_signal;
    delete h_signal;
    delete frame;
    delete result;
}

void prepare_background_model(RooRealVar &kstar, RooDataHist* &background_data, RooHistPdf* &background_pdf, TFile *outfile) {

    using namespace RooFit;

    TFile *file = TFile::Open("/home/galucia/antiLithium4/analysis/output/CATS/CATS_converted.root");
    auto h_correlation_background = (TH1F*)file->Get("hHe3_p_Coul_CF");
    h_correlation_background->SetDirectory(0);
    file->Close();

    background_data = new RooDataHist("data", "data", kstar, Import(*h_correlation_background));
    background_pdf = new RooHistPdf("background_pdf", "background_pdf", kstar, *background_data);

    RooPlot* frame = kstar.frame();
    background_data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    background_pdf->plotOn(frame, LineColor(kBlue), LineWidth(2));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
    TCanvas canvas("canvas_background_fit", "Background Fit", 800, 600);
    canvas.cd();
    frame->Draw();

    outfile->mkdir("background_fit");
    outfile->cd("background_fit");
    h_correlation_background->Write("h_correlation_background");
    canvas.Write();

    delete h_correlation_background;
    delete frame;
}

void prefit_background(RooRealVar &kstar, RooDataHist* &data, RooAddPdf* &model, TFile *outfile, 
                       const char* background_name = "background_pdf", bool do_drawing = false) {

    using namespace RooFit;

    model->chi2FitTo(*data, PrintLevel(-2), SumW2Error(true), Range(0.2, 0.4), 
                     Extended(true), Save());

    if (!do_drawing) {
        return;
    }

    RooPlot* frame = kstar.frame();
    data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    model->plotOn(frame, LineColor(kBlue), LineWidth(2));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");

    TCanvas canvas("canvas_prefit_background", "Prefit Background", 800, 600);
    canvas.cd();
    frame->Draw();

    outfile->cd("prefit_background");
    canvas.Write(background_name);
    delete frame;
}

double fit(RooRealVar &kstar, RooDataHist* &data, RooAddPdf* &model, TFile *outfile, 
          const char* fit_name = "fit_result", bool do_drawing = false) {

    using namespace RooFit;

    RooFitResult* result = model->chi2FitTo(*data, PrintLevel(-2), SumW2Error(true), Range(0.02, 0.4), 
                                        Extended(true), Save());

    RooPlot* frame = kstar.frame();
    data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    model->plotOn(frame, LineColor(kBlue), LineWidth(2));
    model->paramOn(frame, Layout(0.45, 0.85, 0.35), Format("NEU", AutoPrecision(2)));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");

    double chi2 = frame->chiSquare();

    if (do_drawing) {
        TCanvas canvas("canvas_fit", "Fit Result", 800, 600);
        canvas.cd();
        frame->Draw();

        outfile->cd("fit_results");
        canvas.Write(fit_name);
    }
    
    delete frame;
    delete result;

    return chi2;
}

double compute_yield(CrystalBallPars & signal_pars, const double nsig, TH1F*& h_mixed_event, TFile *outfile, const char * yield_name = "yield_result") {

    TF1* f_signal = new TF1("f_signal", CrystalBall, 0.02, 0.4, 6);
    f_signal->SetParameters(signal_pars.mean.getVal(), signal_pars.sigma.getVal(), 
                            signal_pars.aL.getVal(), signal_pars.nL.getVal(), 
                            signal_pars.aR.getVal(), signal_pars.nR.getVal());

    TH1F* h_same_event_signal = (TH1F*)h_mixed_event->Clone("h_same_event_signal");
    for (int ibin = 1; ibin <= h_same_event_signal->GetNbinsX(); ++ibin) {
        const double x = h_same_event_signal->GetBinCenter(ibin);
        const double signal_value = f_signal->Eval(x) * nsig * h_mixed_event->GetBinContent(ibin);
        h_same_event_signal->SetBinContent(ibin, signal_value);
    }

    double yield = h_same_event_signal->Integral(1, h_same_event_signal->GetNbinsX());

    outfile->cd("yield_results");
    h_same_event_signal->Write(yield_name);

    return yield;

}

void compute_yield_error() {

    //const int N_ITERATIONS = 1000000; // 1 million
    const int N_ITERATIONS = 10000; // 10k

    TH1F* h_same, * h_mixed;
    load_same_mixed(h_same, h_mixed);
    auto h_correlation = (TH1F*)h_same->Clone("h_correlation");
    std::cout << "Cloned histogram for correlation function." << std::endl;
    compute_correlation_function(h_same, h_mixed, h_correlation);

    std::cout << "Loaded histograms: " << h_same->GetName() << " and " << h_mixed->GetName() << std::endl;

    auto h_same_iter = (TH1F*)h_same->Clone("h_same_iter");
    auto h_mixed_iter = (TH1F*)h_mixed->Clone("h_mixed_iter");
    auto h_correlation_iter = (TH1F*)h_same_iter->Clone("h_correlation_iter");
    std::cout << "Cloned histograms for iterations." << std::endl;

    auto outfile = TFile::Open("/home/galucia/antiLithium4/femto/output/yield_error_Anti_mc.root", "RECREATE");
    //auto outfile = TFile::Open("/home/galucia/antiLithium4/femto/output/yield_error_Matter_mc.root", "RECREATE");

    CrystalBallPars pars;
    prepare_signal_fit(pars, outfile);
    pars.setConstantTail();
    const double mean_mc = pars.mean.getVal();
    const double sigma_mc = pars.sigma.getVal();
    
    auto h_raw_yield = new TH1F("h_raw_yield", "Raw yield distribution;#it{N}_{raw}(^{4}Li);Counts", 500, 0, 10000);
    auto h_chi2_fit = new TH1F("h_chi2_fit", "Chi2 distribution;#chi^{2};Counts", 500, 0, 100);

    RooRealVar kstar("kstar", "kstar", 0.02, 0.4);
    RooCrystalBall signal_pdf("signal_pdf", "signal", kstar, pars.mean, pars.sigma, pars.aL, pars.nL, pars.aR, pars.nR);
    RooDataHist* background_data = nullptr;
    RooHistPdf* background_pdf = nullptr;

    prepare_background_model(kstar, background_data, background_pdf, outfile);

    RooRealVar nsig("nsig", "Signal Yield", 0, 0, 1e3);
    RooRealVar nbkg("nbkg", "Background Yield", 1, 0, 1e6);
    RooAddPdf* model = new RooAddPdf("model", "Signal + Background", RooArgList(signal_pdf, *background_pdf), RooArgList(nsig, nbkg));
    
    // Example fit to the correlation function
    RooDataHist* data = new RooDataHist("data", "data", kstar, RooFit::Import(*h_correlation));
    nsig.setVal(0);
    nsig.setConstant(true);
    prefit_background(kstar, data, model, outfile);
    nsig.setConstant(false);
    nbkg.setConstant(true);

    outfile->mkdir("prefit_background");
    outfile->mkdir("fit_results");
    outfile->mkdir("yield_results");
    outfile->mkdir("high_yield_results");
    bool do_drawing = true; // Set to false to skip drawing

    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        if (iter % 100 == 0) {
            std::cout << "Processing iteration: " << iter << "/" << N_ITERATIONS << std::endl;
            do_drawing = true;
        } else {
            do_drawing = false;
        }

        for (int ibin = 1; ibin <= h_same->GetNbinsX(); ++ibin) {
            h_same_iter->SetBinContent(ibin, gRandom->Poisson(h_same->GetBinContent(ibin)));
            h_mixed_iter->SetBinContent(ibin, gRandom->Poisson(h_mixed->GetBinContent(ibin)));
        }
        compute_correlation_function(h_same_iter, h_mixed_iter, h_correlation_iter);

        RooDataHist* data_iter = new RooDataHist("data_iter", "data_iter", kstar, RooFit::Import(*h_correlation_iter));
        
        nsig.setVal(0);
        nsig.setConstant(true);
        prefit_background(kstar, data_iter, model, outfile, Form("background_pdf_%d", iter), do_drawing);
        nsig.setConstant(false);
        nbkg.setConstant(true);

        pars.mean.setVal(mean_mc);
        pars.sigma.setVal(sigma_mc);

        double chi2 = fit(kstar, data_iter, model, outfile, Form("fit_result_%d", iter), do_drawing);
        h_chi2_fit->Fill(chi2);

        double yield  = compute_yield(pars, nsig.getVal(), h_mixed_iter, outfile, Form("yield_result_%d", iter));
        h_raw_yield->Fill(yield);

        if (yield > 1000) {
            RooPlot* frame_high_yield = kstar.frame();
            data_iter->plotOn(frame_high_yield, RooFit::MarkerStyle(20), RooFit::MarkerSize(0.5));
            model->plotOn(frame_high_yield, RooFit::LineColor(kBlue), RooFit::LineWidth(2));
            model->paramOn(frame_high_yield, RooFit::Layout(0.45, 0.85, 0.35), RooFit::Format("NEU", RooFit::AutoPrecision(2)));
            frame_high_yield->SetTitle(Form("High Yield Fit Result - Iteration %d;#it{k}* (GeV/#it{c});C(#it{k}*)", iter));
            TCanvas canvas_high_yield(Form("canvas_high_yield_%d", iter), "High Yield Fit Result", 800, 600);
            canvas_high_yield.cd();
            frame_high_yield->Draw();
            outfile->cd("high_yield_results");
            canvas_high_yield.Write(Form("canvas_high_yield_%d", iter));
            delete frame_high_yield;
        }
        
        delete data_iter;

    }
    
    outfile->cd();
    h_raw_yield->Write();
    h_chi2_fit->Write();
    outfile->Close();

}