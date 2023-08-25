#define pion_cxx
#include "pion.hpp"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include "TVector3.h"
#include "TLorentzVector.h"

void pion::Loop(std::string filename)
{
//   In a ROOT session, you can do:
//      root> .L pion.C
//      root> pion t
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
   if (fChain == 0) return;

   TFile *top = new TFile(filename.c_str(), "recreate");
   auto tNew = new TTree("pion", "Cutted tr_ph (pi+pi- events)");

   auto hMTracks = new TH1D("hMTracks", "hMTracks", 1000, -200, 200);
   auto hCalEnergy_bhabha = new TH1D("hCalEnergy_bhabha", "hCalEnergy_bhabha", 1000, 0, 1000);

   double phiPos;
   double thetaPos;
   double momPos;
   double phiNeg;
   double thetaNeg;
   double momNeg;
   double CalEnergyPos;
   double CalEnergyNeg;

   double phiPos_v;
   double thetaPos_v;
   double momPos_v;
   double phiNeg_v;
   double thetaNeg_v;
   double momNeg_v;

   tNew->Branch("emeas", &emeas, "emeas/F");

   tNew->Branch("phiPos", &phiPos, "phiPos/D");
   tNew->Branch("thetaPos", &thetaPos, "thetaPos/D");
   tNew->Branch("momPos", &momPos, "momPos/D");
   tNew->Branch("phiNeg", &phiNeg, "phiNeg/D");
   tNew->Branch("thetaNeg", &thetaNeg, "thetaNeg/D");
   tNew->Branch("momNeg", &momNeg, "momNeg/D");
   tNew->Branch("CalEnergyPos", &CalEnergyPos, "CalEnergyPos/D");
   tNew->Branch("CalEnergyNeg", &CalEnergyNeg, "CalEnergyNeg/D");

   tNew->Branch("phiPos_v", &phiPos_v, "phiPos_v/D");
   tNew->Branch("thetaPos_v", &thetaPos_v, "thetaPos_v/D");
   tNew->Branch("momPos_v", &momPos_v, "momPos_v/D");
   tNew->Branch("phiNeg_v", &phiNeg_v, "phiNeg_v/D");
   tNew->Branch("thetaNeg_v", &thetaNeg_v, "thetaNeg_v/D");
   tNew->Branch("momNeg_v", &momNeg_v, "momNeg_v/D");


   const double cutChi2r = 15.;
   const double cutChi2z = 10.;
   const int cutNhitMin = 10;
   const double cutRmin = 0.05;
   const double cutRmax = 0.3;
   const double cutZtrack = 12.;
   const double cutPtot = 40;
   const double cutTrackTheta = 0.7;

   int counter = 0;
   TVector3 pos3(1, 1, 1);
   TVector3 neg3(1, 1, 1);
   TLorentzVector pos(1, 1, 1, 1);
   TLorentzVector neg(1, 1, 1, 1);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(nt != 2)
      { continue; }

      for(int i = 0; i < nt; i++)
      {
         if (tptot[i] > cutPtot && fabs(trho[i]) < cutRmax && fabs(tz[i]) < cutZtrack &&
            tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin &&
            tdedx[i] < 5e3)
         { counter++; }


         if(counter > 2)
         { break; }
      }

      if(counter == 2 && is_coll == 1)
      {
         auto posIndex = tcharge[0] > 0? 0 : 1;
         auto negIndex = posIndex == 0? 1 : 0;
         phiPos = tphi[posIndex];
         thetaPos = tth[posIndex];
         momPos = tptot[posIndex];
         phiNeg = tphi[negIndex];
         thetaNeg = tth[negIndex];
         momNeg = tptot[negIndex]; 

         CalEnergyPos = ten[posIndex];
         CalEnergyNeg = ten[negIndex];

         phiPos_v = tphiv[posIndex];
         thetaPos_v = tthv[posIndex];
         momPos_v = tptotv[posIndex];
         phiNeg_v = tphiv[negIndex];
         thetaNeg_v = tthv[negIndex];
         momNeg_v = tptotv[negIndex]; 

         if(is_bhabha == 1)
         { hCalEnergy_bhabha->Fill((ten[posIndex] + ten[negIndex]) / 2); }

         pos3.SetMagThetaPhi(momPos_v, thetaPos_v, phiPos_v);
         neg3.SetMagThetaPhi(momNeg_v, thetaNeg_v, phiNeg_v);

         pos.SetVectM(pos3, 139.570);
         neg.SetVectM(neg3, 139.570);

         hMTracks->Fill((pos + neg).M() - 2 * emeas);

         if(fabs((pos + neg).M() - 2 * emeas) < 0.05 * 2 * emeas && is_bhabha != 1)
         { tNew->Fill(); }
      }

      counter = 0;
   }

   top->Write();
   top->Save();
}
