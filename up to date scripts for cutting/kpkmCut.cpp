#define kpkmCut_cxx
#include "kpkmCut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void kpkmCut::Loop()
{
//   In a ROOT session, you can do:
//      root> .L kpkmCut.C
//      root> kpkmCut t
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

   TFile *top = new TFile("kchCut21May.root", "recreate");
   auto tNew = new TTree("kChargedTree", "Cutted tr_ph (Kch Energy stability important data)");

   auto hist1 = new TH2D("histPlus", "a", 1000, 0, 500, 1000, 0, 20000);
   auto hist2 = new TH2D("histMinus", "a", 1000, 0, 500, 1000, 0, 20000);

   Float_t dpsi;
   Float_t Y;
   Float_t eff;
   tNew->Branch("emeas", &emeas0, "emeas/F");
   tNew->Branch("demeas", &demeas0, "demeas/F");
   tNew->Branch("runnum", &runnum, "runnum/I");
   tNew->Branch("tptot", tptot, "tptot[2]/F");
   tNew->Branch("tdedx", tdedx, "tdedx[2]/F");

   Double_t halfPi = TMath::Pi() / 2;
   int NgoodTr = 0;
   int NgoodTrS = 0;
   double cutChi2r = 15.;
   double cutChi2z = 10.;
   int cutNhitMin = 10;
   int cutNhitMax = 30;
   double cutRmin = 0.1;
   double cutRmax = 6.0;
   double cutZtrack = 12.;
   double cutPtot = 40;

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
         if (tptot[i] > cutPtot && fabs(trho[i]) < cutRmax && fabs(tz[i]) < cutZtrack &&
               tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin && tnhit[i] < cutNhitMax)
         { NgoodTr++; }
      }

      if (NgoodTr == 2)
      { NgoodTrS++; }

      if (NgoodTr == 2 && is_coll==1 && tnhit[0]>10 && tnhit[1]>10 && tcharge[0]*tcharge[1] < 0 &&
         (tth[0]*tcharge[0] + tth[1]*tcharge[1] + TMath::Pi()) / 2 < TMath::Pi() - 1 && 
         (tth[0]*tcharge[0] + tth[1]*tcharge[1] + TMath::Pi()) / 2 > 1 && 
         abs(tth[1] - TMath::Pi()/2) <= 0.9  &&
         (tdedx[0]+tdedx[1])/2 > 7000 &&
         fabs(tptot[0]-tptot[1])/(tptot[0]+tptot[1]) < 0.3 &&
         tcharge[0] * tcharge[1] < 0)
      {
         tNew->Fill();

         if(tcharge[0] > 0)
         {
            hist1->Fill(tptot[0],tdedx[0]);
            hist2->Fill(tptot[1],tdedx[1]);
         }
         else
         {
            hist2->Fill(tptot[0],tdedx[0]);
            hist1->Fill(tptot[1],tdedx[1]);
         }
      }

      NgoodTr = 0;
   }

   top->Write();
   top->Save();
}
