#include <string>
#include <TString.h>
#include <TSystem.h>

void load(TString myopt = "fast") {
    gSystem->AddIncludePath((std::string("-I ")+"/home/galucia/antiLithium4/mixed_event/build").c_str());
    TString opt;
    if(myopt.Contains("force"))   opt = "kfg";
    else                          opt = "kg";
  
    gSystem->CompileMacro(".L /home/galucia/antiLithium4/mixed_event/TreeManager.cpp+", opt.Data(), "", "build");
    gSystem->CompileMacro(".L /home/galucia/antiLithium4/mixed_event/MixedEventRoutine.cpp+", opt.Data(), "", "build");
}