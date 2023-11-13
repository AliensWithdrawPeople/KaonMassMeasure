#define xsec_cut_cxx
#include "xsec_cut.hpp"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>
#include <TLine.h>

#include <iostream>
#include <set>
#include <map>

#include "C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/KnXSec.cpp"

void xsec_cut::Loop(std::string histFileName)
{
    //   In a ROOT session, you can do:
    //      root> .L ksklCut_mcgpj.C
    //      root> ksklCut_mcgpj t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    // by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    auto xsec = new KnXSec();
    TFile *new_file = new TFile(("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/xsec_cutted/" + histFileName).c_str(), "recreate");
    auto tnew = fChain->CloneTree(0);

    auto energy_spectrum = new TH1D("energy_spectrum", "energy_spectrum", 3500, 500, 535);
    fChain->Draw("emeas>>energy_spectrum", "", "goff");
    std::set<double> energies = {};

    Long64_t nentries = fChain->GetEntriesFast();
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        fChain->GetEntry(jentry);
        energies.insert(emeas);
    }
    std::map<double, double> xsec_vals; 
    for(const auto &en : energies)
    { xsec_vals[en] = xsec->eval(en); }
    auto max_xsec = std::max_element(xsec_vals.begin(), xsec_vals.end(), 
        [](std::pair<double, double> a, std::pair<double, double> b)
        {
            return a.second < b.second;
        })->second;

    std::map<double, int> events_limits;
    std::map<double, int> events;
    for(auto& [en, xsection] : xsec_vals)
    {
        events_limits[en] = int(xsection / max_xsec * fChain->GetEntries(("fabs(emeas - " + std::to_string(en) + ") < 0.05").c_str()));
        events[en] = 0;
    }


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        if(events[emeas] < events_limits[emeas])
        {
            events[emeas]++;
            tnew->Fill();
        }

    }
    std::cout << tnew->GetEntries() << std::endl;
    new_file->Write();
    new_file->Save();
}