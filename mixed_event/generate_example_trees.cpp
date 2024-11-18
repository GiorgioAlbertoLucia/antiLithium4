#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>

#include <vector>

/**
 * This function creates a file with two folders each containing two trees with a single branch.
*/
void generate_example_trees(const char* outputFileName) {
    // Create a new ROOT file
    TFile *file = TFile::Open(outputFileName, "RECREATE");

    std::vector<TDirectory*> directories;
    std::vector<TTree*> trees;
    for (int idir = 0; idir < 2; ++idir) {
        // Create a new directory
        TDirectory *dir = file->mkdir(Form("dir%d", idir));
        dir->cd();
        directories.push_back(dir);

        // Create new trees
        for (int itree = 0; itree < 2; ++itree) {
            TTree *tree = new TTree(Form("Tree%d", itree), Form("Tree%d", itree));
            int x;
            tree->Branch(Form("x%d", itree), &x);
            for (int i = 0; i < 10; ++i) {
                x = i + itree;
                tree->Fill();
            }
            tree->Write();
            trees.push_back(tree);
        }
    }

    // Close the file
    file->Close();

}