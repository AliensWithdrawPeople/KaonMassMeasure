#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int Xsec_EnergySpectrum_Cutter(std::string point)
{
    // Cut events so energy spectum is a convolution of xsection and gaussian distribution (energy smearing).
    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/KsKl/Xsec_EnergySpectrumCut/xsec_cut.cpp");
    auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/MCGPJ_kskl" + point + "_Merged.root";
    // auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 phi_width 4.5 MeV/MCGPJ_kskl" + point + "_Merged_phi_width_4.5_MeV.root";
    auto file = TFile::Open((fname).c_str());
    gROOT->ProcessLine("xsec_cut a(tr_ph_merged)");
    gROOT->ProcessLine(("a.Loop(\"MCGPJ_kskl" + point + "_Merged_XsecConv.root\")").c_str());
    return 0;
}