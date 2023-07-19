#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int massCutter(std::string point) {
    // MC
    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/kskl2bgen.cpp");
    auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl" + point + "_Merged.root";
    auto file = TFile::Open((fname).c_str());
    gROOT->ProcessLine("kskl2bGen a(tr_ph_merged)");
    gROOT->ProcessLine(("a.Loop(\"KsKl_Smeared/MC" + point + ".root\")").c_str());
    
    // EXP
    // gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/ksklExp.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/KnPrelim/Prelim" + point + ".root";
    // auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("ksklExp a(tr_ph)");
    // gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/expKsKl/exp" + point + ".root\")").c_str());

    return 0;
}