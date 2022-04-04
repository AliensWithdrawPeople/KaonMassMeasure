#define KSKL_cxx
#include "KSKL.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>

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
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    Float_t dpsi; Float_t Y;
    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", &dpsi, "dpsi/F");
    tNew->Branch("Y", &Y, "Y/F");

    Double_t halfPi = TMath::Pi() / 2;

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
        
        if(nt == 2 && nks == 1 && is_coll!=1 && tnhit[0]>10 && tnhit[1]>10 && (tdedx[0]+tdedx[1])/2 < 5000 && kstlen[0]>0.1 && kstlen[0]<1.3 && ksalign[0]>0.9 &&
        abs(tth[0] - TMath::Pi()/2) <= 0.5 && abs(tth[1] - TMath::Pi()/2) <= 0.5 && tcharge[0]*tcharge[1] < 0 &&
        fabs(tz[0]) < 10 && tchi2r[0] < 8 && fabs(tz[1]) < 10 && tchi2r[1] < 8 && tchi2z[1] < 4 && tchi2z[1] < 4 &&
        sqrt(emeas*emeas - ksptot[0]*ksptot[0]) < 520 && sqrt(emeas*emeas - ksptot[0]*ksptot[0]) > 460)
        {
            if(tcharge[0] > 0)
            { Y = tptot[0] / tptot[1]; }
            else
            { Y = tptot[1] / tptot[0]; }
            dpsi = ksdpsi[0];
            tNew->Fill();
        }
    }

    top->Write();
    top->Save();
}