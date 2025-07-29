#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TDecompChol.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TString.h>
#include <iostream>
#include <vector>

#include "crystal_ball.hh"

double F(double k, const TVectorD& p) {
    double p3 = 1.3011;
    double p4 = 2.7734;
    double p5 = 1.1105;
    double p6 = 5.7432;

    double xs[1] = {k};
    double pars[6] = {p[1], p[2], p3, p4, p5, p6};
    return p[0] * CrystalBall(xs, pars);
}

double integrate(const std::vector<double>& k, const std::vector<double>& S) {
    double result = 0.0;
    for (size_t i = 0; i < k.size() - 1; ++i) {
        double dk = k[i+1] - k[i];
        result += 0.5 * dk * (S[i] + S[i+1]);
    }
    return result;
}

void MonteCarloYield(const char* filename = "yourfile.root", const char* histname = "M_hist", const char* outfile = "output.root") {

    TFile* file = TFile::Open(filename);
    TH1F* h_Mixed = (TH1F*)file->Get(histname);
    h_Mixed->SetDirectory(0);  // Detach histogram from file
    file->Close();

    const int N_samples = 1000;
    const int N_params = 3;
    const int N_bins = h_Mixed->GetNbinsX();

    TVectorD param_mean(N_params);
    TMatrixD cov(N_params, N_params);
    
    // matter
    param_mean[0] = 2.8588e-1; // counts
    param_mean[1] = 8.1759e-2; // mean
    param_mean[2] = 1.3137e-2;  // sigma

    cov(0,0) = 5.565;       cov(0,1) = -1.435e-8;   cov(0,2) = -1.224e-8; 
    cov(1,0) = -1.435e-8;   cov(1,1) =  9.142e-10;  cov(1,2) = -6.303e-17;
    cov(2,0) = -1.224e-8;   cov(2,1) = -6.303e-17;  cov(2,2) = 3.982e-10;

    // antimatter
    //param_mean[0] = 5.0644e-1; // counts
    //param_mean[1] = 7.9246e-2; // mean
    //param_mean[2] = 1.401e-2;  // sigma

    //cov(0,0) = 6.631;       cov(0,1) = 6.629e-08;   cov(0,2) = -1.89e-06; 
    //cov(1,0) = 6.629e-08;   cov(1,1) = 7.143e-07;   cov(1,2) = 2.31e-13;
    //cov(2,0) = -1.89e-06;   cov(2,1) = 2.31e-13;    cov(2,2) = 1.474e-08;

    TDecompChol chol(cov);
    TMatrixD L = chol.GetU();  // upper triangle, take transpose for lower
    TMatrixD LT(N_params, N_params);
    LT.Transpose(L);

    TRandom3 rng(42);
    std::vector<double> Y_samples;
    TH1F* h_Signal = (TH1F*)h_Mixed->Clone("S_samples");

    TFile* outFile = TFile::Open(outfile, "RECREATE");
    h_Mixed->Write("M_mixed");
    auto outDir = outFile->mkdir("Iterations");
    outDir->cd();

    for (int i = 0; i < N_samples; ++i) {
        std::cout << "Sample " << i+1 << "/" << N_samples << std::endl;

        // Sample parameter vector
        TVectorD z(N_params);
        for (int j = 0; j < N_params; ++j)
            z[j] = rng.Gaus(0, 1);

        TVectorD params = param_mean;
        for (int j = 0; j < N_params; ++j)
            for (int k = 0; k < N_params; ++k)
                params[j] += LT(j,k) * z[k];

        // Sample M(k) and compute S(k)
        for (int ibin = 1; ibin < N_bins+1; ++ibin) {
            //double M_sample = h_Mixed->GetBinContent(ibin) + rng.Gaus(0, h_Mixed->GetBinError(ibin));
            double M_sample = h_Mixed->GetBinContent(ibin);
            h_Signal->SetBinContent(ibin, F(h_Mixed->GetBinCenter(ibin), params) * M_sample); // ibin+1 because ROOT histograms are 1-indexed
        }
        h_Signal->Write(Form("S_sample_%d", i+1));
        double yield = h_Signal->Integral(1, N_bins);

        Y_samples.push_back(yield);
    }

    // Compute mean and standard deviation
    double mean = 0.0;
    for (auto y : Y_samples) mean += y;
    mean /= N_samples;

    double var = 0.0;
    for (auto y : Y_samples) var += (y - mean)*(y - mean);
    var /= (N_samples - 1);

    std::cout << "Y = " << mean << " ± " << std::sqrt(var) << std::endl;

    TVectorD result(2);
    result[0] = mean; // mean yield
    result[1] = std::sqrt(var); // standard deviation
    outFile->cd();
    result.Write("Yield");

    outFile->Close();
}

void runMonteCarlo() {
    //MonteCarloYield("/home/galucia/antiLithium4/femto/output/AntiLithium4FitCFCATS.root", 
    //                "dir/h_mixed_event_norm1.0",
    //                "output/AntiLithium4Yield.root");

    MonteCarloYield("/home/galucia/antiLithium4/femto/output/MatterLithium4FitCFCATS.root", 
                    "dir/h_mixed_event_norm1.0",
                    "output/MatterLithium4Yield.root");
}
