#define KChargedCut_cxx
#include "KChargedCut.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>

void KChargedCut::Loop(std::string histFileName)
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

    TFile *fNew = new TFile(histFileName.c_str(), "recreate");
    
    
    auto tNew = new TTree("kChargedTree", "Cutted for K+K- tr_ph");

    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("tptot", tptot, "tptot[2]/F");
    tNew->Branch("tdedx", tdedx, "tdedx[2]/F");
    
/*
   auto hKchRho = new TH1D("KchRho", "Kch Rho", 400, -3.0, 3.0);
   auto hKchTh = new TH1D("KchTheta", "Kch theta", 200, 0, 4.0);
   auto hKchPavg = new TH1D("KchPavg", "KchPavg", 1000, 0, 500);

   auto hRho = new TH1D("KsRho", "KsRho", 200, -1.0, 10.0);
   auto hdedx = new TH1D("dedx", "Ks dedx", 2000, 0.0, 25000.0);
   auto hKsTh = new TH1D("KsTheta", "Track theta (for Ks)", 200, 0, 3.5);
   auto hKsLen = new TH1D("KsLen", "KsLen", 1000, 0.0, 10.0);
   auto hAlign = new TH1D("ksalign", "Ks Align", 110, 0.0, 1.1);
   auto hKsPavg = new TH1D("KsPavg", "KsPavg", 1000, 0, 500);

   auto hCut1 = new TH2D("cut1", "dedx : ptot (for Kcharged)", 1000, 0.0, 1000, 1000, 0, 25000);
   auto hCut2 = new TH2D("cut2", "p1 : p2", 1000, 0, 1000, 1000, 0, 1000);
   */


   Double_t halfPi = TMath::Pi() / 2;

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

        if(nt == 2 && is_coll==1 && tnhit[0]>10 && tnhit[1]>10 && tcharge[0]*tcharge[1] < 0 &&
        TMath::Abs(trho[0]) < 0.3 && TMath::Abs(trho[1]) < 0.3 &&
        (tth[0] - tth[1] + TMath::Pi())/2 < TMath::Pi() && (tth[0] - tth[1] + TMath::Pi())/2 > 1 &&
        abs(tth[0] - TMath::Pi()/2) <= 0.4 && abs(tth[1] - TMath::Pi()/2) <= 0.4  && (tdedx[0]+tdedx[1])/2 > 7000 &&
        fabs(tptot[0]-tptot[1])/(tptot[0]+tptot[1]) < 0.3 && tptot[0] > 70 && tptot[1] > 70) 
        { tNew->Fill(); }
        /*
        if(nt == 2 && tcharge[0] * tcharge[1] < 0) 
        { 
            hKchRho->Fill(trho[0]); hKchRho->Fill(trho[1]);
            hKchPavg->Fill((tptot[0] + tptot[1]) / 2);
            hKchTh->Fill(tth[0]); hKchTh->Fill(tth[1]);
            hCut1->Fill((tptot[0] + tptot[1]) / 2, (tdedx[0] + tdedx[1]) / 2);
        }

        if(nks == 1 && nt == 2 && tcharge[0] * tcharge[1] < 0 && is_coll != 1 &&
        kstlen[0] > 0.05 && ) 
        { 
            hRho->Fill(kstlen[0]); hdedx->Fill(tdedx[0]);
            hKsTh->Fill(tth[1]); hKsTh->Fill(tth[0]);
            hAlign->Fill(ksalign[0]); hKsLen->Fill(kslen[0]);
            hKsPavg->Fill((tptot[1] + tptot[0]) / 2);
            hCut2->Fill(tptot[ksvind[0][0]], tptot[ksvind[0][1]]);
        }
        */
    }

    fNew->Write();
    fNew->Save();
}