#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>

int cutter(std::string point) {
    gROOT->LoadMacro("C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/PhiToKn/PhiToKn.cpp");
    auto fname = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/KnPrelim/Prelim" + point + ".root";
    auto file = TFile::Open((fname).c_str());
    gROOT->ProcessLine("PhiToKn a(tr_ph)");

    gROOT->ProcessLine(("a.Loop(\"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn" + point + ".root\", " + point + ")").c_str());
    return 0;
}