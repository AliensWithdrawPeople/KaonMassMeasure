#define pipiCut_cxx
#include "pipiCut.h"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVector3.h>
#include <Math/SpecFuncMathCore.h>

void pipiCut::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      root> .L pipiCut.C
//      root> pipiCut t
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

      TFile *top = new TFile(histFileName.c_str(), "recreate");
      auto tNew = new TTree("pionsTree", "Cutted tr_ph (pi+pi- events)");

    Float_t p1; Float_t p2;
    tNew->Branch("emeas", &emeas, "emeas/F");
    tNew->Branch("demeas", &demeas, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("tphi", tphi, "tphi[10]/F");
    tNew->Branch("tth", tth, "tth[10]/F");
    tNew->Branch("tptot", tptot, "tptot[10]/F");
    tNew->Branch("tcharge", tcharge, "tcharge[10]/I");

    int NgoodTr = 0;
    int NgoodTrS = 0;
    double cutChi2r = 15.;
    double cutChi2z = 10.;
    int cutNhitMin = 10;
    int cutNhitMax = 30;
    double cutRmin = 0.05;
    double cutRmax = 6;
    double cutZtrack = 12.;
    double cutPtot = 40;
    double cutTrackTheta = 0.7;

    double pionMomentum0 = 0;
    std::vector<int> nTracks;

    auto hDeltaMom = new TH1D("hDeltaMom", "Lorentz delta mom", 1000, -0.2, 0.2);
    auto hDeltaMomVsDeltaPhi = new TH2D("hDeltaMomVsDeltaPhi", "delta mom vs delta phi", 2000, -20, 20, 2000, -3.15, 3.15);
    auto hDeltaMomVsDeltaPhiCowboy = new TH2D("hDeltaMomVsDeltaPhiCowboy", "delta mom vs delta phi Cowboy", 2000, -20, 20, 2000, -3.15, 3.15);
    auto hDeltaMomVsDeltaPhiSailor = new TH2D("hDeltaMomVsDeltaPhiSailor", "delta mom vs delta phi Sailor", 2000, -20, 20, 2000, -3.15, 3.15);
    auto hDeltaPhiVsPhiCowboy = new TH2D("hDeltaPhiVsPhiCowboy", "delta phi vs phi (rec - gen) Cowboy", 1000, 0, 2 * 3.15, 1000, -1, 1);
    auto hDeltaPhiVsPhiSailor = new TH2D("hDeltaPhiVsPhiSailor", "delta phi vs phi (rec - gen) Sailor", 1000, 0, 2 * 3.15, 1000, -1, 1);
    auto hDeltaPhiRecVsGenCowboy = new TH2D("hDeltaPhiRecVsGenCowboy", "delta phi Rec vs Gen data Cowboy", 1000, 0, 6.3, 1000, 0, 6.3);
    auto hDeltaPhiRecVsGenSailor = new TH2D("hDeltaPhiRecVsGenSailor", "delta phi Rec vs Gen data Sailor", 2000, 0, 6.3, 2000, 0, 6.3);

    auto hPhiColGen = new TH1D("hPhiColGen", "", 5000, -0.15, 0.15);
    auto hPhiColRec = new TH1D("hPhiColRec", "", 1000, -0.02, 0.02);
    auto hPtotGen = new TH1D("hPtotGen", "PtotGen", 1000, 200, 300);
    auto hPtotRec = new TH1D("hPtotRec", "PtotRec", 1000, 200, 300);

    auto hDeltaPhiPos = new TH1D("hDeltaPhiPos", "tphiv - tphi pi+", 1000, -0.2, 0.2); 
    auto hDeltaPhiNeg = new TH1D("hDeltaPhiNeg", "tphiv - tphi pi-", 1000, -0.2, 0.2); 
    auto hDeltaPtotVsDeltaPhiPos = new TH2D("hDeltaPtotVsDeltaPhiPos", "tptotv - tptot vs tphiv - tphi pi+", 5000, -0.2, 0.2, 5000, -10, 10);
    auto hDeltaPtotVsDeltaPhiNeg = new TH2D("hDeltaPtotVsDeltaPhiNeg", "tptotv - tptot vs tphiv - tphi pi-", 5000, -0.2, 0.2, 5000, -10, 10);
    auto hDeltaPhiVsPhiPos = new TH2D("hDeltaPhiVsPhiPos", "tphiv - tphi vs tphiv pi+", 1000, 0, 2 * TMath::Pi(), 1000, -0.2, 0.2);
    auto hDeltaPhiVsPhiNeg = new TH2D("hDeltaPhiVsPhiNeg", "tphiv - tphi vs tphiv pi-", 1000, 0, 2 * TMath::Pi(), 1000, -0.2, 0.2);

    auto hDeltaPhiVsThetaPos = new TH2D("hDeltaPhiVsThetaPos", "tphiv - tphi vs tthv pi+", 1000, 0, TMath::Pi(), 1000, -0.2, 0.2);
    auto hDeltaPhiVsThetaNeg = new TH2D("hDeltaPhiVsThetaNeg", "tphiv - tphi vs tthv pi-", 1000, 0, TMath::Pi(), 1000, -0.2, 0.2);
    
    TVector3 piPos;
    TVector3 piNeg;
    TVector3 piPosRec;
    TVector3 piNegRec;
    TVector3 field(0., 0., 1.);
    int posTrackNumber = 0;
    int negTrackNumber = 0;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        for (int i = 0; i < nt; i++)
        {
            if (tptot[i] > cutPtot && fabs(trho[i]) < cutRmax && fabs(tz[i]) < cutZtrack &&
            tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin && tnhit[i] < cutNhitMax)
            { 
                NgoodTr++;
                nTracks.push_back(i);
            }
        }
        pionMomentum0 = sqrt(emeas * emeas - 139.57 * 139.57);
        if(NgoodTr == 2 && is_coll == 1 && tptotv[nTracks[0]] > pionMomentum0 - 40 && tptotv[nTracks[1]] > pionMomentum0 - 40 &&
                                            tptotv[nTracks[0]] < pionMomentum0 + 20 && tptotv[nTracks[1]] < pionMomentum0 + 20)
        {
            posTrackNumber = tcharge[nTracks[0]] > 0 ? 0 : 1;
            negTrackNumber = posTrackNumber == 1 ? 0 : 1;

            piPosRec.SetMagThetaPhi(tptotv[nTracks[posTrackNumber]], tthv[nTracks[posTrackNumber]], tphiv[nTracks[posTrackNumber]]);
            piNegRec.SetMagThetaPhi(tptotv[nTracks[negTrackNumber]], tthv[nTracks[negTrackNumber]], tphiv[nTracks[negTrackNumber]]);

            hDeltaPtotVsDeltaPhiPos->Fill(tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]], 
                                            tptotv[nTracks[posTrackNumber]] - tptot[nTracks[posTrackNumber]]);
            hDeltaPtotVsDeltaPhiNeg->Fill(tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]], 
                                            tptotv[nTracks[negTrackNumber]] - tptot[nTracks[negTrackNumber]]);

            hDeltaPhiPos->Fill(tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            hDeltaPhiNeg->Fill(tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);

            hDeltaPhiVsPhiPos->Fill(tphiv[nTracks[posTrackNumber]], tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            hDeltaPhiVsPhiNeg->Fill(tphiv[nTracks[negTrackNumber]], tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);

            hDeltaPhiVsThetaPos->Fill(tthv[nTracks[posTrackNumber]], tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            hDeltaPhiVsThetaNeg->Fill(tthv[nTracks[negTrackNumber]], tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);

            hPhiColRec->Fill(fabs(piPosRec.Phi() - piNegRec.Phi()) - TMath::Pi());
            hPtotRec->Fill(piPosRec.Mag());

            // for(int i = 0; i < nsim; i++)
            // {
            //     if(simtype[i] == 211)
            //     { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }

            //     if(simtype[i] == -211)
            //     { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
            // }
        
            // hPtotGen->Fill(piPos.Mag());
            // hPhiColGen->Fill(fabs(piPos.Phi() - piNeg.Phi()) - TMath::Pi());
            // int mult = piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()) < 0;
            // hDeltaMomVsDeltaPhi->Fill(fabs(piPos.Mag() / piNeg.Mag()) - fabs(piPosRec.Mag() / piNegRec.Mag()),
            //                         fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())) - fabs(piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()))
            //                         - 0 * TMath::Power(-1, mult) * TMath::Pi() );
            // //  Cowboy Type
            // if(piPosRec.Cross(field).XYvector().DeltaPhi(piNegRec.XYvector()) < TMath::Pi() / 2)
            // { 
            //     hDeltaPhiVsPhiCowboy->Fill(piPos.XYvector().Phi(), piPos.XYvector().DeltaPhi(piPosRec.XYvector()));
            //     hDeltaPhiRecVsGenCowboy->Fill(fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())), 
            //                                 fabs(piPosRec.XYvector().DeltaPhi(piNegRec.XYvector())) );

            //     // hPhi->Fill(piPos.Phi() - piPosRec.Phi()); 
            //     // hPhi->Fill(piNeg.Phi() - piNegRec.Phi()); 
            //     hDeltaMomVsDeltaPhiCowboy->Fill(fabs(piPos.Mag() - piNeg.Mag()) - fabs(piPosRec.Mag() - piNegRec.Mag()),
            //                         piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector())
            //                         - TMath::Power(-1, mult) * TMath::Pi() );
            // }
            // //  Sailor Type
            // if(piPosRec.Cross(field).XYvector().DeltaPhi(piNegRec.XYvector()) > TMath::Pi() / 2)
            // {  
            //     hDeltaPhiVsPhiSailor->Fill(piPos.XYvector().Phi(), piPos.XYvector().DeltaPhi(piPosRec.XYvector()));

            //     hDeltaPhiRecVsGenSailor->Fill(fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())), 
            //                                 fabs(piPosRec.XYvector().DeltaPhi(piNegRec.XYvector())) );
            //     hDeltaMomVsDeltaPhiSailor->Fill(fabs(piPos.Mag() - piNeg.Mag()) - fabs(piPosRec.Mag() - piNegRec.Mag()),
            //                         piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector())
            //                         - TMath::Power(-1, mult) * TMath::Pi() );
            // }

            tNew->Fill();
        }
        NgoodTr = 0;
        nTracks.clear();
        nTracks.shrink_to_fit();
    }

    auto hyperGauss = new TF1("HyperGauss", "[0] * exp([1] * (1 - sqrt(1 + (x - [2]) * (x - [2]) / [3])))", -1000, 1000);
    // 2 * K(1, 2) / K(2, 2) = 1.10234; [3]_0 = sigma^2 * 2 * K(1, 2) / K(2, 2)
    // hyperGauss->SetParameters(hPhiColGen->GetMaximum(), 4, 0, hPhiColGen->GetRMS() * hPhiColGen->GetRMS() * 1.10234);
    // hPhiColGen->Fit(hyperGauss, "SM", "", - 2 * hPhiColGen->GetRMS(), 2 * hPhiColGen->GetRMS());
    
    // hyperGauss->SetParameters(hPhiColRec->GetMaximum(), 4., 0, hPhiColRec->GetRMS() * hPhiColRec->GetRMS() * 1.10234);
    // hPhiColRec->Fit(hyperGauss, "SLQM", "", -2 * hPhiColRec->GetRMS(), 2 * hPhiColRec->GetRMS());

    // hyperGauss->SetParameters(hPtotRec->GetMaximum(), 2., hPtotRec->GetBinCenter(hPtotRec->GetMaximumBin()), 
    // hPtotRec->GetRMS() * hPtotRec->GetRMS() * 1.10234);
    // hPtotRec->Fit(hyperGauss, "SM", "", hPtotRec->GetMean() - hPtotRec->GetRMS(), hPtotRec->GetMean() + hPtotRec->GetRMS());

    top->Write();
    top->Save();
}