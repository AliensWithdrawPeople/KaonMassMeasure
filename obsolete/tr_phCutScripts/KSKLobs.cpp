#define KSKL_cxx
#include "KSKL.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>

void KSKL::Loop(std::string histFileName)
{
  //   In a ROOT session, you can do:
  //      root> .L KSKL.C
  //      root> KSKL t
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
  //by  b_branchname->GetEntry(ientry); //read only this branch
    TFile *top = new TFile(histFileName.c_str(), "recreate");
    TH1D *hNgood = new TH1D("hNgood", "Number of good tracks", 10, 0, 10);
    TH1D *hKsdPsi = new TH1D("hKsdPsi", "Ks dpsi", 960, 0, 3.2);
    TH1D *hEmeas = new TH1D("emeas", "Measured Beam energy", 1000, 509, 510);
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("nks", &nks, "nks/I");
    tNew->Branch("ksminv", ksminv, "ksminv[nks]/F");
    tNew->Branch("ksalign", ksalign, "ksalign[nks]/F");
    tNew->Branch("ksdpsi", ksdpsi, "ksdpsi[nks]/F");
    tNew->Branch("kslen", kslen, "kslen[nks]/F");

    double cutChi2r = 20.;
    int cutNhit = 10;
    double cutRmax = 6.0;
    double cutZtrack = 10.;
    double cutPtot = 1000;

    int NgoodTr = 0;

    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
        break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        for (int i = 0; i < nt_total; i++)
        {
            if(tptot[i] < cutPtot && trho[i] < cutRmax && fabs(tz[i]) < cutZtrack &&
             tchi2r[i] < cutChi2r && tnhit[i] > cutNhit)
            { NgoodTr++; }
        }

        if(NgoodTr == 2 && ksalign[0] > 0.8 && kstlen[0] <= 1.7 && ksminv[0] < 570 && ksminv[0] > 420 && ksptot[0] < 180 && ksptot[0] > 40)
        { 
            hKsdPsi->Fill(ksdpsi[0]); 
            hEmeas->Fill(emeas0);
            tNew->Fill();
        }

        hNgood->Fill(NgoodTr);
        NgoodTr = 0;
    }

    top->Write();
    top->Save();
}