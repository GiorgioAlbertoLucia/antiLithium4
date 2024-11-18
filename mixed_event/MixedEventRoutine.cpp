#include <iostream>
#include <vector>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>

#include "TreeManager.cpp"

void MixedEventRoutine(const char* inputFileName, const char* outputFileName, std::vector<std::string>& treeNames) {
    
    ROOT::RDataFrame df = HStackTreeInDataFrame(inputFileName, treeNames);
    
    // Define the mixed event routine here
    // For example, the following code will print the number of entries in the dataframe
    std::cout << "Number of entries: " << df.Count().GetValue() << std::endl;
    
    // Save the dataframe to a file
    TFile * outputFile = TFile::Open(outputFileName, "RECREATE");
    df.Snapshot("mixed_event", outputFile->GetName());
    outputFile->Close();
    
    return;
}
