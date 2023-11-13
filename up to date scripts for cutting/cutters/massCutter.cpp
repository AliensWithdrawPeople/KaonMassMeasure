#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int massCutter(std::string point) {
    // MC
    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/KsKl/kskl2bgen.cpp");
    auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/xsec_cutted/MCGPJ_kskl" + point + "_Merged_XsecConv.root";
    auto file = TFile::Open((fname).c_str());
    gROOT->ProcessLine("kskl2bGen a(tr_ph_merged)");
    gROOT->ProcessLine(("a.Loop(\"KsKl_Smeared/New formfactor/XsecConv/MC" + point + "_XsecConv.root\")").c_str());

    // MC without xsection-energy spectrum convolution
    // gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/KsKl/kskl2bgen.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/MCGPJ_kskl" + point + "_Merged.root";
    // auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("kskl2bGen a(tr_ph_merged)");
    // gROOT->ProcessLine(("a.Loop(\"KsKl_Smeared/New formfactor/MC" + point + ".root\")").c_str());

    // MC with heterogeneous field
    // gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/KsKl/kskl2bgen.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl" + point + "_Merged_Field_New.root";
    // auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("kskl2bGen a(tr_ph_merged)");
    // gROOT->ProcessLine(("a.Loop(\"KsKl_Smeared/Field/MC" + point + "_Field.root\")").c_str());
    
    // EXP
    // gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/KsKl/ksklExp.cpp");
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/KnPrelim/Prelim" + point + ".root";
    // auto file = TFile::Open((fname).c_str());
    // gROOT->ProcessLine("ksklExp a(tr_ph)");
    // gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/expKsKl/exp" + point + ".root\")").c_str());

    return 0;
}