//
//  Script to evaluate the efficiency and the transverse momentum resolution of the selection
//

#include <string>

#include <Riostream.h>

#include <TEfficiency.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>

/**
 * @brief Evaluate the efficiency of the selection. The efficiency is defined as Li4-reconstructed / Li4-total.
 * From the lithium4analtsis.cpp, the reconstructed candidates have specific values for some variables. 
 * Non-reconstructed candidates will have default values for these variables.
 * (default: fHe3NSigmaTPC = -10, fPrNSigmaTPC = -10)
 *  
*/
void EvaluateEfficiency(TTree & inTree, TDirectory & outDir, const char * dirName)
{

    TH1F lithium4("lithium4", "lithium4; Mass [GeV/#it{c}^{2}]; Counts", 20, 3.7, 3.9);                            // selected with flag fIsBkg_ == 0
    TH1F lithium4Total("lithium4Total", "lithium4Total; Mass [GeV/#it{c}^{2}]; Counts", 20, 3.7, 3.9);

    inTree.Draw("fMassMC >> lithium4", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");
    inTree.Draw("fMassMC >> lithium4Total");
    
    TEfficiency efficiency(lithium4, lithium4Total);
    efficiency.SetTitle("Efficiency; Mass [GeV/#it{c}^{2}]; Efficiency");
    outDir.cd();
    efficiency.Write("efficiency");
    lithium4.Write("lithium4");
    lithium4Total.Write("lithium4Total");
}

/**
 * Momentum resolution is defined as the difference between the reconstructed and the true momentum, 
 * divided by the true momentum.
*/
void EvaluatePtResolution(TTree & inTree, TDirectory & outDir, const char * dirName)
{
    const int nBinsX = 100;
    const int nBinsY = 100;
    const int nEntries = inTree.GetEntries();

    TH1F ptHe3("ptHe3", "ptHe3; p_{T, reco} [GeV/#it{c}]; Counts", nBinsX, 0., 10.);
    TH1F ptHe3True("ptHe3True", "ptHe3True; p_{T, true} [GeV/#it{c}]; Counts", nBinsX, 0., 10.);
    TH1F ptPr("ptPr", "ptPr; p_{T, reco} [GeV/#it{c}]; Counts", nBinsX, 0., 10.);
    TH1F ptPrTrue("ptPrTrue", "ptPrTrue; p_{T, true} [GeV/#it{c}]; Counts", nBinsX, 0., 10.);

    // resolution
    TH2F ptResHe3("ptResHe3", "ptResHe3; p_{T, true} [GeV/#it{c}]; (p_{T, true} - p_{T, reco}) / p_{T, true}", nBinsX, 0., 10., nBinsY, -.3, .3);
    TH2F ptResPr("ptResPr", "ptResPr; p_{T, true} [GeV/#it{c}]; (p_{T, true} - p_{T, reco}) / p_{T, true}", nBinsX, 0., 10., nBinsY, -.3, .3);

    inTree.Draw("fPtHe3 >> ptHe3", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");
    inTree.Draw("fPtTrueHe3 >> ptHe3True", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");
    inTree.Draw("fPtPr >> ptPr", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");
    inTree.Draw("fPtTruePr >> ptPrTrue", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");

    inTree.Draw("(fPtTrueHe3 - fPtHe3)/fPtTrueHe3:fPtTrueHe3 >> ptResHe3", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");
    inTree.Draw("(fPtTruePr - fPtPr)/fPtTruePr:fPtTruePr >> ptResPr", "(fHe3NSigmaTPC > -10) && (fPrNSigmaTPC > -10)");

    outDir.cd();
    
    ptHe3.Write("ptHe3");
    ptHe3True.Write("ptHe3True");
    ptResHe3.Write("ptResHe3");
    
    ptPr.Write("ptPr");
    ptPrTrue.Write("ptPrTrue");
    ptResPr.Write("ptResPr");
}

/*
void efficiency_ptresolution(const int folderIdx)
{   
    const char * inFilePath = "/Users/glucia/Projects/ALICE/antiLithium4/MCWorkflowLauncher/AO2D_lit_mc.root";
    const char * outFilePath = "/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/output/efficiency_ptresolution.root";

    TFile * inFile = TFile::Open(inFilePath);
    TKey * key = (TKey *) inFile->GetListOfKeys()->At(folderIdx);
    TDirectory * dir = (TDirectory *) key->ReadObj();
    std::string dirName = dir->GetName();

    TTree * tree = (TTree *) dir->Get("O2lithium4tablemc");
    //tree->Print();

    TFile * outFile = TFile::Open(outFilePath, "UPDATE");
    TDirectory * outDir = outFile->mkdir(dirName.c_str());

    EvaluateEfficiency(*tree, *outDir, dirName.c_str());
    EvaluatePtResolution(*tree, *outDir, dirName.c_str());

    outFile->Close();
    inFile->Close();

}
*/

void efficiency_ptresolution()
{
    const char * inFilePath = "/Users/glucia/Projects/ALICE/antiLithium4/MCWorkflowLauncher/AO2D_lit_mc.root";
    const char * outFilePath = "/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/output/efficiency_ptresolution.root";

    TFile * inFile = TFile::Open(inFilePath);
    TFile * outFile = TFile::Open(outFilePath, "RECREATE");
    const int nKeys = inFile->GetListOfKeys()->GetEntries();
    auto listKeys = inFile->GetListOfKeys();

    for (int i = 0; i < nKeys-1; i++)
    {
        TKey * key = (TKey *) listKeys->At(i);
        if(!key) continue;
        TDirectory * dir = (TDirectory *) key->ReadObj();
        std::string dirName = dir->GetName();

        TTree * tree = (TTree *) dir->Get("O2lithium4tablemc");
        //tree->Print();

        TDirectory * outDir = outFile->mkdir(dirName.c_str());

        EvaluateEfficiency(*tree, *outDir, dirName.c_str());
        EvaluatePtResolution(*tree, *outDir, dirName.c_str());

        delete tree;
    }

    outFile->Close();
    inFile->Close();
}
