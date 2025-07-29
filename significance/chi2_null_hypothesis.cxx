/**
 * This code computes the chi2 disitrbution for the null hypothesis 
*/

#include <iostream>
#include <vector>
#include <deque>
#include <numeric>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TMath.h>


void load_same_mixed(TH1F *& h_same, TH1F *& h_mixed) {

    TFile *file_mixed = TFile::Open("/home/galucia/antiLithium4/analysis/output/PbPb/studies.root");
    h_mixed = (TH1F*)file_mixed->Get("CorrelationAnti/hMixed_kstar");
    h_mixed->SetDirectory(0);
    file_mixed->Close();

    TFile *file_correction = TFile::Open("/home/galucia/antiLithium4/analysis/output/CATS/CATS_new.root");
    auto h_correction = (TH1F*)file_correction->Get("hHe3_p_Coul_CF_LS");
    h_correction->SetDirectory(0);
    file_correction->Close();

    h_same = (TH1F*)h_mixed->Clone("h_same");
    h_same->SetDirectory(0);

    for (int ibin = 1; ibin <= h_same->GetNbinsX()+1; ++ibin) {
        double mixed_value = h_mixed->GetBinContent(ibin);
        double bin_center = h_mixed->GetBinCenter(ibin);
        double correction_value = h_correction->GetBinContent(h_correction->FindBin(bin_center));
        h_same->SetBinContent(ibin, mixed_value * correction_value);
    }

    delete h_correction;
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
            h_correlation->SetBinContent(ibin, 0.0);
        }
    }
}

void run_mc_chi2(TH1F *& h_chi2, TH1F *& h_chi2_far_from_signal,
                 std::vector<TH1F *>& running_chi2_histograms,
                 std::vector<TH1F *>& window_chi2_histograms,
                 float kstar_bin_centers[],
                 TFile * outfile,
                 const int N_BINS = 40 /* 40 bins */,
                 const int N_ITERATIONS = 1000000 /* 1 mln */) {

    TH1F* h_same, * h_mixed;
    load_same_mixed(h_same, h_mixed);
    auto h_correlation = (TH1F*)h_same->Clone("h_correlation");
    std::cout << "Cloned histogram for correlation function." << std::endl;
    compute_correlation_function(h_same, h_mixed, h_correlation);

    std::cout << "Loaded histograms: " << h_same->GetName() << " and " << h_mixed->GetName() << std::endl;

    auto h_same_iter = (TH1F*)h_same->Clone("h_same_iter");
    auto h_mixed_iter = (TH1F*)h_mixed->Clone("h_mixed_iter");
    auto h_correlation_iter = (TH1F*)h_same_iter->Clone("h_correlation_iter");

    const int N_WINDOW_BINS = running_chi2_histograms.size() - window_chi2_histograms.size();
    std::deque<float> chi2_deque;
    chi2_deque.resize(N_WINDOW_BINS);

    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        if (iter % 10000 == 0) {
            std::cout << "Processing iteration: " << iter << "/" << N_ITERATIONS << std::endl;
        }
        h_same_iter->Reset();
        h_mixed_iter->Reset();

        for (int ibin = 1; ibin <= h_same->GetNbinsX(); ++ibin) {

            h_same_iter->SetBinContent(ibin, gRandom->Poisson(h_same->GetBinContent(ibin)));
            h_mixed_iter->SetBinContent(ibin, gRandom->Poisson(h_mixed->GetBinContent(ibin)));
        }
        compute_correlation_function(h_same_iter, h_mixed_iter, h_correlation_iter);

        double chi2 = 0.0, kstar = 0.0, expected = 0.0, observed = 0.0, error = 0.0;
        double chi2_cumulated = 0.0;
        double chi2_limited = 0;
        double chi2_far_from_signal = 0;
        //const int FIRST_BIN = h_correlation->FindBin(0.01);
        const int FIRST_BIN = h_correlation->FindBin(0.02);

        for (int ibin = FIRST_BIN; ibin <= h_same_iter->GetNbinsX(); ++ibin) {
            kstar = h_same_iter->GetBinCenter(ibin);
            expected = h_correlation_iter->GetBinContent(ibin);  
            observed = h_correlation->GetBinContent(ibin);
            error = h_correlation_iter->GetBinError(ibin);

            if (error > 0) {
                chi2 = (observed - expected) * (observed - expected) / (error * error);
            }
            chi2_cumulated += chi2;

            if (kstar < 0.15) {
                chi2_limited = chi2_cumulated;
                //chi2_limited = chi2_cumulated / (ibin+1); // reduced chi2
            } else {
                chi2_far_from_signal += chi2;
            }
            running_chi2_histograms[ibin-1]->Fill(chi2_cumulated);
            
            chi2_deque.pop_front();
            chi2_deque.push_back(chi2);
            if (ibin > N_WINDOW_BINS) {
                double chi2_window = std::accumulate(chi2_deque.begin(), chi2_deque.end(), 0., std::plus<double>());
                window_chi2_histograms[ibin - N_WINDOW_BINS - 1]->Fill(chi2_window);
            }
            //running_chi2_histograms[ibin-1]->Fill(chi2 / (ibin+1)); // reduced chi2

            if (iter == 0) {
                kstar_bin_centers[ibin-1] = h_correlation->GetBinCenter(ibin);
            }
        }

        h_chi2->Fill(chi2_limited);
        h_chi2_far_from_signal->Fill(chi2_far_from_signal);
    }

    outfile->cd();
    h_same->Write();
    h_mixed->Write();
    h_correlation->Write();
    h_same_iter->Write();
    h_mixed_iter->Write();
    h_correlation_iter->Write();

    delete h_same;
    delete h_mixed;
    delete h_correlation;
    delete h_same_iter;
    delete h_mixed_iter;
    delete h_correlation_iter;

}

void display_running_result(TH1F *& h_chi2,
                            std::vector<TH1F *> & running_chi2_histograms,
                            TFile * outfile,
                            float kstar_bin_centers[],
                            const int N_BINS = 40 /* 40 bins */,
                            const int N_ITERATIONS = 1000000 /* 1 mln */) {


    // antimatter
    //float running_chi2[] = {1.653, 1.944, 2.056, 2.061, 3.552,
    //                        5.436, 9.618, 16.31, 16.50, 16.68, 
    //                        18.60, 19.07, 19.48, 22.14, 24.63, 
    //                        24.79, 30.41, 30.89, 31.42, 38.91, 
    //                        39.77, 43.52, 43.68, 44.08, 44.65, 
    //                        54.65, 55.15, 55.21, 55.24, 55.46, 
    //                        58.89, 59.68, 63.17, 72.19, 72.26, 
    //                        72.62, 73.18, 73.21, 74.14, 74.15};
    // matter
    float running_chi2[] = {0.607, 5.360, 8.214, 8.295, 8.32, 
                            8.595, 11.86, 12.43, 13.86, 15.39, 
                            15.50, 15.60, 16.77, 16.87, 17.09, 
                            19.78, 20.08, 20.09, 20.23, 20.40, 
                            20.63, 21.31, 21.36, 21.79, 23.29, 
                            24.35, 24.92, 25.26, 25.82, 26.88, 
                            28.10, 29.81, 32.56, 34.21, 34.79, 
                            34.99, 35.48, 35.49, 36.29, 36.37};
    // matter 0.02
    //float running_chi2[] = {2.854, 2.935, 2.962, 3.235, 6.504, 7.067, 8.496, 10.025, 10.144, 10.240, 11.411, 11.504, 11.726, 14.415, 14.723, 14.734, 14.872, 15.035, 15.274, 15.946, 16.001, 16.434, 17.926, 18.991, 19.561, 19.896, 20.464, 21.515, 22.739, 24.446, 27.200, 28.845, 29.431, 29.626, 30.118, 30.132, 30.925, 31.011};


    outfile->cd();
    h_chi2->Write();

    outfile->mkdir("running_chi2");
    outfile->cd("running_chi2");

    TGraph *g_running_pvalue = new TGraph(N_BINS);
    g_running_pvalue->SetTitle("Running P-value;#it{k}* (GeV/#it{c});P-value");
    TGraph *g_running_significance = new TGraph(N_BINS);
    g_running_significance->SetTitle("Running Significance;#it{k}* (GeV/#it{c});Significance");

    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        running_chi2_histograms[ibin]->Write();
        float kstar = kstar_bin_centers[ibin];
        float chi2_value = running_chi2[ibin];
        //float chi2_value = running_chi2[ibin] / (ibin+1); // reduced chi2
        float pvalue = running_chi2_histograms[ibin]->Integral(running_chi2_histograms[ibin]->FindBin(chi2_value), running_chi2_histograms[ibin]->GetNbinsX()+1) / N_ITERATIONS;
        float significance = TMath::NormQuantile(1. - pvalue/2.);

        g_running_pvalue->SetPoint(ibin, kstar, pvalue);
        g_running_significance->SetPoint(ibin, kstar, significance);
    }
    g_running_pvalue->Write("g_running_pvalue");
    g_running_significance->Write("g_running_significance");

}

void display_window_result(TH1F *& h_chi2,
                           std::vector<TH1F *> & window_chi2_histograms,
                           TFile * outfile,
                           float kstar_bin_centers[],
                           const int N_BINS = 40 /* 40 bins */,
                           const int N_BINS_WINDOW = 6 /* 6 bins for the window */,
                           const int N_ITERATIONS = 1000000 /* 1 mln */) {

    auto infile_window = TFile::Open("/home/galucia/antiLithium4/femto/output/NEW_AntiLithium4FitCFCATS.root", "READ");
    auto h_window_chi2_data = (TH1F*)infile_window->Get("dir/running_chi2");
    h_window_chi2_data->SetDirectory(0);
    infile_window->Close();

    outfile->mkdir("window_chi2");
    outfile->cd("window_chi2");
    TGraph *g_window_pvalue = new TGraph(N_BINS - N_BINS_WINDOW);
    g_window_pvalue->SetTitle("Window P-value;#it{k}* (GeV/#it{c});P-value");
    TGraph *g_window_significance = new TGraph(N_BINS - N_BINS_WINDOW);
    g_window_significance->SetTitle("Window Significance;#it{k}* (GeV/#it{c});Significance");
    
    for (int ibin = 0; ibin < N_BINS - N_BINS_WINDOW; ++ibin) {

        float chi2_value = h_window_chi2_data->GetBinContent(ibin + N_BINS_WINDOW / 2 + 1);
        float pvalue = window_chi2_histograms[ibin]->Integral(window_chi2_histograms[ibin]->FindBin(chi2_value), window_chi2_histograms[ibin]->GetNbinsX()+1) / N_ITERATIONS;
        float significance = TMath::NormQuantile(1. - pvalue/2.);

        float kstar = kstar_bin_centers[ibin + N_BINS_WINDOW / 2];
        g_window_pvalue->SetPoint(ibin, kstar, pvalue);
        g_window_significance->SetPoint(ibin, kstar, significance);

        window_chi2_histograms[ibin]->Write(Form("h_window_chi2_%d", ibin));
    }
    g_window_pvalue->Write("g_window_pvalue");
    g_window_significance->Write("g_window_significance");
}

void compute_chi2() {

    const int N_ITERATIONS = 1000000; // 1 million
    const int N_BINS = 40; // 40 bins
    const int N_BINS_WINDOW = 4; // 6 bins for the window

    auto outfile = TFile::Open("/home/galucia/antiLithium4/significance/chi2_output_matter_0p02.root", "RECREATE");
    auto h_chi2 = new TH1F("h_chi2", "Chi2 Distribution;#chi^{2};Counts", 1000, 0, 100);
    auto h_chi2_far_from_signal = new TH1F("h_chi2_far_from_signal", "Chi2 Distribution (far from signal);#chi^{2};Counts", 1000, 0, 100);

    std::vector<TH1F*> running_chi2_histograms, window_chi2_histograms;
    running_chi2_histograms.reserve(N_BINS);
    window_chi2_histograms.reserve(N_BINS - N_BINS_WINDOW);
    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        std::string name_running = "h_running_chi2_" + std::to_string(ibin);
        auto h_running_chi2 = new TH1F(name_running.c_str(), Form("Running Chi2 %d ;#chi^{2};Counts", ibin), 1000, 0, 100);
        running_chi2_histograms.emplace_back(h_running_chi2);

        if (ibin > N_BINS_WINDOW / 2 && ibin <= N_BINS - (N_BINS_WINDOW / 2)) {
            std::string name_window = "h_window_chi2_" + std::to_string(ibin - N_BINS_WINDOW);
            auto h_window_chi2 = new TH1F(name_window.c_str(), Form("Window Chi2 %d ;#chi^{2};Counts", ibin - N_BINS_WINDOW), 1000, 0, 100);
            window_chi2_histograms.emplace_back(h_window_chi2);
        }
    }

    float kstar_bin_centers[N_BINS];
    run_mc_chi2(h_chi2, h_chi2_far_from_signal, running_chi2_histograms, window_chi2_histograms,
                kstar_bin_centers, outfile, N_BINS, N_ITERATIONS);

    std::cout << std::endl;
    // chi2 computed for kstar < 0.15 GeV/c
    const float chi2_data_antimatter_MARIO = 23.765; // k* (0.01, 0.15) 
    //const float chi2_data_antimatter_MARIO = 24.627; 
    //const float chi2_data_antimatter_MARIO = 24.627 / 13; // reduced chi2
    const float pvalue_null_antimatter_MARIO = h_chi2->Integral(h_chi2->FindBin(chi2_data_antimatter_MARIO), h_chi2->GetNbinsX()+1) / N_ITERATIONS;
    std::cout << "P-value for null hypothesis  (MARIO): " << pvalue_null_antimatter_MARIO << std::endl;
    std::cout << "Significance  (MARIO): " << TMath::NormQuantile(1. - pvalue_null_antimatter_MARIO/2.) << std::endl;

    const float chi2_data_matter_MARIO = 12.361;
    //const float chi2_data_matter_MARIO = 12.361 / 13; // reduced chi2
    const float pvalue_null_matter_MARIO = h_chi2->Integral(h_chi2->FindBin(chi2_data_matter_MARIO), h_chi2->GetNbinsX()+1) / N_ITERATIONS;
    std::cout << "P-value for null hypothesis: " << pvalue_null_matter_MARIO << std::endl;
    std::cout << "Significance: " << TMath::NormQuantile(1. - pvalue_null_matter_MARIO/2.) << std::endl;

    const float chi2_data_antimatter_far_from_signal = 49.5236;
    //const float chi2_data_antimatter_far_from_signal = 49.5236 / 25; // reduced chi2
    const float pvalue_null_antimatter_far_from_signal = h_chi2_far_from_signal->Integral(h_chi2_far_from_signal->FindBin(chi2_data_antimatter_far_from_signal), h_chi2_far_from_signal->GetNbinsX()+1) / N_ITERATIONS;
    std::cout << "P-value for null hypothesis (far from signal): " << pvalue_null_antimatter_far_from_signal << std::endl;
    std::cout << "Significance (far from signal): " << TMath::NormQuantile(1. - pvalue_null_antimatter_far_from_signal/2.) << std::endl;

    display_running_result(h_chi2, running_chi2_histograms, outfile, kstar_bin_centers, N_BINS, N_ITERATIONS);

    display_window_result(h_chi2, window_chi2_histograms, outfile, kstar_bin_centers, N_BINS, N_BINS_WINDOW, N_ITERATIONS);

    outfile->Close();

}