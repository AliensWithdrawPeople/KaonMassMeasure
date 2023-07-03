#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int cutter_MC(std::string point) {
    // gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/prelim.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl" + point + "_Merged.root";
    // auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("prelim a(tr_ph_merged)");

    // gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KnPrelim/MC_Kn" + point + ".root\")").c_str());

    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/PhiToKn.cpp");
    auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KnPrelim/MC_Kn" + point + ".root";
    auto file = TFile::Open((fname).c_str());
    gROOT->ProcessLine("PhiToKn a(tr_ph_merged)");

    gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/MC/MC_Kn" + point + ".root\", " + point + ")").c_str());
    return 0;
}