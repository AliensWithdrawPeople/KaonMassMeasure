#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int TrackEffCutter(std::string point) {
    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/KsKl/ksklTrEff.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl" + point + "_Merged.root";
    std::string fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MultiHadron/tr_ph_run000030.root";
    auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("ksklTrEff a(tr_ph_merged)");
    gROOT->ProcessLine("ksklTrEff a(tr_ph)");

    gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/TrackEfficiency/MC/TrEff_MC" + point + "_mult.root\")").c_str());
    return 0;
}