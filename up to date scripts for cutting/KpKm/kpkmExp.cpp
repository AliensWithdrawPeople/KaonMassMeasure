#define kpkmExp_cxx
#include "kpkmExp.hpp"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TVector3.h>

void kpkmExp::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      root> .L kpkmExp.C
//      root> kpkmExp t
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
   if (fChain == 0)
        return;

   TFile *top = new TFile(("../tr_ph/MC/Kch/" + histFileName).c_str(), "recreate");
   auto tNew = new TTree("kChargedTree", "Cutted tr_ph (Kch Energy stability important data)");

   auto hKp = new TH2D("hKp", "K+", 5000, 0, 1000, 5000, 0, 40000);
   auto hKm = new TH2D("hKm", "K-", 5000, 0, 1000, 5000, 0, 40000);
   auto hMomentums =  new TH2D("hMoms", "K+ mom vs K- mom", 1000, 0, 1000, 1000, 0, 1000);

   Float_t dpsi;
   Float_t Y;
   Float_t eff;
   tNew->Branch("ebeam", &ebeam, "ebeam/F");
   tNew->Branch("emeas", &emeas0, "emeas/F");
   tNew->Branch("demeas", &demeas0, "demeas/F");
   tNew->Branch("runnum", &runnum, "runnum/I");
   tNew->Branch("tptot", tptotv, "tptot[2]/F");
   tNew->Branch("tdedx", tdedx, "tdedx[2]/F");

   Double_t halfPi = TMath::Pi() / 2;
   int NgoodTr = 0;
   int NgoodTrS = 0;
   const double cutChi2r = 15.;
   const double cutChi2z = 10.;
   const int cutNhitMin = 10;
   const double cutZtrack = 10.;
   const double cutPtot = 40;
   
   int counter = 0;

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

      for (int i = 0; i < nt; i++)
      {
         if (tptotv[i] > cutPtot && fabs(tz[i]) < cutZtrack &&
               tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin)
         { NgoodTr++; }
      }

      if (NgoodTr == 2)
      { NgoodTrS++; }

      if (NgoodTr == 2 && is_coll==1 && tnhit[0]>10 && tnhit[1]>10 && tcharge[0]*tcharge[1] < 0 &&
         tth[0] > 1 && tth[0] < TMath::Pi() - 1 &&
         tth[1] > 1 && tth[1] < TMath::Pi() - 1 &&
         tdedx[0] > 40 * exp(-(tptotv[0] - 60) / 40) + 7000 &&
         tdedx[1] > 40 * exp(-(tptotv[1] - 60) / 40) + 7000 &&
         fabs(tptotv[0]-tptotv[1])/(tptotv[0]+tptotv[1]) < 0.3 &&
         tptotv[0] > 60 && tptotv[0] < 150 &&
         tptotv[1] > 60 && tptotv[1] < 150)
      {
         tNew->Fill();
         counter++;
         if(tcharge[0] > 0)
         {
            hMomentums->Fill(tptotv[1], tptotv[0]);
            hKp->Fill(tptotv[0],tdedx[0]);
            hKm->Fill(tptotv[1],tdedx[1]);
         }
         else
         {
            hMomentums->Fill(tptotv[1], tptotv[0]);
            hKm->Fill(tptotv[0],tdedx[0]);
            hKp->Fill(tptotv[1],tdedx[1]);
         }
      }

      NgoodTr = 0;
   }
   std::cout << "Number of events in the tree = " << nentries << std::endl;
   std::cout << "Number of good tracks = " << NgoodTrS << std::endl;
   std::cout << "efficiency = " << double(counter) / nentries << std::endl;
   top->Write();
   top->Save();
}