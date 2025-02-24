#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int cutter_MC(std::string point) {
    // gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/prelim.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/xsec_cutted/MCGPJ_kskl" + point + "_Merged_XsecConv.root";
    // auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("prelim a(tr_ph_merged)");
    // gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/KnPrelim/MC/New formfactor/MC_Prelim" + point + ".root\")").c_str());

    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/PhiToKn/PhiToKn_MC.cpp");
    auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/KnPrelim/MC/New formfactor/MC_Prelim" + point + ".root";
    auto file = TFile::Open((fname).c_str());
    gROOT->ProcessLine("PhiToKn_MC a(tr_ph_merged)");

    gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/MC/MC_Kn" + point + ".root\", " + point + ")").c_str());

    return 0;
}

// C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MultiHadron/tr_ph_525_62741_multihadr_no_KsKl_2E1.0498_5.root
// a.Loop("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MultiHadron/tr_ph_525_62741_proc.root", 525.)