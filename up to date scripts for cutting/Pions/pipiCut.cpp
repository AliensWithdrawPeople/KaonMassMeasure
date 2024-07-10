#define pipiCut_cxx
#include "pipiCut.hpp"
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
#include "TProfile.h"
#include <Math/SpecFuncMathCore.h>

void printVector(std::vector<double> &vec, std::string name)
{
    std::cout << name + " = {";
    for(auto elem : vec)
    { std::cout << elem << ", "; }
    std::cout << "}" << std::endl;
}

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
    tNew->Branch("tdedx", tdedx, "tdedx[10]/F");
    tNew->Branch("trho", trho, "trho[10]/F");
    tNew->Branch("tcharge", tcharge, "tcharge[10]/I");

    int NgoodTr = 0;
    int NgoodTrS = 0;
    const double cutChi2r = 15.;
    const double cutChi2z = 10.;
    const int cutNhitMin = 10;
    const int cutNhitMax = 30;
    const double cutRmin = 0.05;
    const double cutRmax = 6;
    const double cutZtrack = 12.;
    const double cutPtot = 40;
    const double cutTrackTheta = 0.7;

    double pionMomentum0 = 0;
    std::vector<int> nTracks;
    tNew->Branch("nTracks", &nTracks);


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
    
    auto hDeltaThetaPos = new TH1D("hDeltaThetaPos", "tthv - tth pi+", 1000, -0.02, 0.02);
    auto hDeltaThetaNeg = new TH1D("hDeltaThetaNeg", "tthv - tth pi-", 1000, -0.02, 0.02);

    auto hDeltaPhiVsThetaPos = new TH2D("hDeltaPhiVsThetaPos", "tphiv - tphi vs tthv pi+", 1000, 0, TMath::Pi(), 1000, -0.2, 0.2);
    auto hDeltaPhiVsThetaNeg = new TH2D("hDeltaPhiVsThetaNeg", "tphiv - tphi vs tthv pi-", 1000, 0, TMath::Pi(), 1000, -0.2, 0.2);
    
    // Vector of histograms containing (deltaPhi / deltaPtot vs Phi) info. 
    // [0], [1] - pi+ and pi+ respectively; 
    std::vector<TH2D *> vDeltaPhiVsPhi{new TH2D("hDeltaPhiVsPhiPos", "tphiv - tphi vs tphiv pi+", 36, 0, 2 * TMath::Pi(), 1000, -0.2, 0.2), 
                                        new TH2D("hDeltaPhiVsPhiNeg", "tphiv - tphi vs tphiv pi-", 36, 0, 2 * TMath::Pi(), 1000, -0.2, 0.2)};
    std::vector<TH2D *> vDeltaPtotVsPhi{new TH2D("hDeltaPtotVsPhiPos", "tptotv - tptotv vs tphiv pi+", 36, 0, 2 * TMath::Pi(), 4000, -20, 20), 
                                        new TH2D("hDeltaPtotVsPhiNeg", "tptotv - tptotv vs tphiv pi-", 36, 0, 2 * TMath::Pi(), 4000, -20, 20)};

    std::vector<TProfile *> vDeltaPhiVsPhiProfPos;
    std::vector<TProfile *> vDeltaPhiVsPhiProfNeg;
    std::vector<TProfile *> vDeltaPtotVsPhiProfPos;
    std::vector<TProfile *> vDeltaPtotVsPhiProfNeg;
    auto hDeltaPhiVsPhiPosFull_pfx = new TProfile("hDeltaPhiVsPhiPos_pfx_Full2pi", "hDeltaPhiVsPhiPos_pfx_Full2pi", 18*150, 0, 2 * TMath::Pi(), -0.2, 0.2);
    auto hDeltaPhiVsPhiNegFull_pfx = new TProfile("hDeltaPhiVsPhiNeg_pfx_Full2pi", "hDeltaPhiVsPhiNeg_pfx_Full2pi", 18*150, 0, 2 * TMath::Pi(), -0.2, 0.2);
    auto hDeltaPtotVsPhiPosFull_pfx = new TProfile("hDeltaPtotVsPhiPos_pfx_Full2pi", "hDeltaPtotVsPhiPos_pfx_Full2pi", 18*150, 0, 2 * TMath::Pi(), -0.2, 0.2);
    auto hDeltaPtotVsPhiNegFull_pfx = new TProfile("hDeltaPtotVsPhiNeg_pfx_Full2pi", "hDeltaPtotVsPhiNeg_pfx_Full2pi", 18*150, 0, 2 * TMath::Pi(), -0.2, 0.2);
    
    int Nchambers = 18;
    for(int i = 0; i < Nchambers; i++)
    { 
        vDeltaPhiVsPhiProfPos.push_back(new TProfile(("hDeltaPhiVsPhiPos_pfx" + std::to_string(i + 1)).c_str(),
                                                    ("hDeltaPhiVsPhiPos_pfx" + std::to_string(i + 1)).c_str(), 
                                                    150, 2 * i * TMath::Pi() / Nchambers, 2 * (i + 1) * TMath::Pi() / Nchambers, -0.1, 0.1)); 
        vDeltaPhiVsPhiProfNeg.push_back(new TProfile(("hDeltaPhiVsPhiNeg_pfx" + std::to_string(i + 1)).c_str(),
                                                    ("hDeltaPhiVsPhiNeg_pfx" + std::to_string(i + 1)).c_str(), 
                                                    150, 2 * i * TMath::Pi() / Nchambers, 2 * (i + 1) * TMath::Pi() / Nchambers, -0.1, 0.1));

        vDeltaPtotVsPhiProfPos.push_back(new TProfile(("hDeltaPtotVsPhiPos_pfx" + std::to_string(i + 1)).c_str(),
                                                    ("hDeltaPtotVsPhiPos_pfx" + std::to_string(i + 1)).c_str(), 
                                                    150, 2 * i * TMath::Pi() / Nchambers, 2 * (i + 1) * TMath::Pi() / Nchambers, -10, 10)); 
        vDeltaPtotVsPhiProfNeg.push_back(new TProfile(("hDeltaPtotVsPhiNeg_pfx" + std::to_string(i + 1)).c_str(),
                                                    ("hDeltaPtotVsPhiNeg_pfx" + std::to_string(i + 1)).c_str(), 
                                                    150, 2 * i * TMath::Pi() / Nchambers, 2 * (i + 1) * TMath::Pi() / Nchambers, -10, 10));                                                                                      
    }

    TVector3 piPos;
    TVector3 piNeg;
    TVector3 piPosRecV;
    TVector3 piNegRecV;
    TVector3 piPosRec;
    TVector3 piNegRec;
    TVector3 field(0., 0., 1.);
    int posTrackNumber = 0;
    int negTrackNumber = 0;
    int indexPos = 0;
    int indexNeg = 0;
    double deltaPhiShiftedPos = 0;
    double deltaPhiShiftedNeg = 0;
    double phiShiftedPos = 0;
    double phiShiftedNeg = 0;

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
            
            piPosRecV.SetMagThetaPhi(tptotv[nTracks[posTrackNumber]], tthv[nTracks[posTrackNumber]], tphiv[nTracks[posTrackNumber]]);
            piNegRecV.SetMagThetaPhi(tptotv[nTracks[negTrackNumber]], tthv[nTracks[negTrackNumber]], tphiv[nTracks[negTrackNumber]]);

            piPosRec.SetMagThetaPhi(tptot[nTracks[posTrackNumber]], tth[nTracks[posTrackNumber]], tphi[nTracks[posTrackNumber]]);
            piNegRec.SetMagThetaPhi(tptot[nTracks[negTrackNumber]], tth[nTracks[negTrackNumber]], tphi[nTracks[negTrackNumber]]);

            hDeltaPtotVsDeltaPhiPos->Fill(tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]], 
                                            tptotv[nTracks[posTrackNumber]] - tptot[nTracks[posTrackNumber]]);
            hDeltaPtotVsDeltaPhiNeg->Fill(tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]], 
                                            tptotv[nTracks[negTrackNumber]] - tptot[nTracks[negTrackNumber]]);

            hDeltaPhiPos->Fill(tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            hDeltaPhiNeg->Fill(tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);
            
            hDeltaPhiVsThetaPos->Fill(tthv[nTracks[posTrackNumber]], tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            hDeltaPhiVsThetaNeg->Fill(tthv[nTracks[negTrackNumber]], tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);

            hDeltaThetaPos->Fill(tthv[nTracks[posTrackNumber]] - tth[nTracks[posTrackNumber]]);
            hDeltaThetaNeg->Fill(tthv[nTracks[negTrackNumber]] - tth[nTracks[negTrackNumber]]);

            hPhiColRec->Fill(fabs(piPosRecV.Phi() - piNegRecV.Phi()) - TMath::Pi());
            hPtotRec->Fill(piPosRecV.Mag());
            
            deltaPhiShiftedPos = 0;
            deltaPhiShiftedNeg = 0;
            phiShiftedPos = TVector2::Phi_0_2pi(tphiv[nTracks[posTrackNumber]]);
            phiShiftedNeg = TVector2::Phi_0_2pi(tphiv[nTracks[negTrackNumber]]);
            indexPos = int(phiShiftedPos / (2 * TMath::Pi() / Nchambers));
            indexNeg = int(phiShiftedNeg / (2 * TMath::Pi() / Nchambers));

            vDeltaPhiVsPhi[0]->Fill(phiShiftedPos, tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            vDeltaPhiVsPhi[1]->Fill(phiShiftedNeg, tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);
            
            vDeltaPtotVsPhi[0]->Fill(phiShiftedPos, tptotv[nTracks[posTrackNumber]] - tptot[nTracks[posTrackNumber]]);
            vDeltaPtotVsPhi[1]->Fill(phiShiftedNeg, tptotv[nTracks[negTrackNumber]] - tptot[nTracks[negTrackNumber]]);

            vDeltaPhiVsPhiProfPos[indexPos]->Fill(phiShiftedPos, tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            vDeltaPhiVsPhiProfNeg[indexNeg]->Fill(phiShiftedNeg, tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);
            
            vDeltaPtotVsPhiProfPos[indexPos]->Fill(phiShiftedPos, tptotv[nTracks[posTrackNumber]] - tptot[nTracks[posTrackNumber]]);
            vDeltaPtotVsPhiProfNeg[indexNeg]->Fill(phiShiftedNeg, tptotv[nTracks[negTrackNumber]] - tptot[nTracks[negTrackNumber]]);

            hDeltaPhiVsPhiPosFull_pfx->Fill(phiShiftedPos, tphiv[nTracks[posTrackNumber]] - tphi[nTracks[posTrackNumber]]);
            hDeltaPhiVsPhiNegFull_pfx->Fill(phiShiftedNeg, tphiv[nTracks[negTrackNumber]] - tphi[nTracks[negTrackNumber]]);
            
            hDeltaPtotVsPhiPosFull_pfx->Fill(phiShiftedPos, tptotv[nTracks[posTrackNumber]] - tptot[nTracks[posTrackNumber]]);
            hDeltaPtotVsPhiNegFull_pfx->Fill(phiShiftedNeg, tptotv[nTracks[negTrackNumber]] - tptot[nTracks[negTrackNumber]]);

            // for(int i = 0; i < nsim; i++)
            // {
            //     if(simtype[i] == 211)
            //     { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }

            //     if(simtype[i] == -211)
            //     { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
            // }
        
            // hPtotGen->Fill(piPos.Mag());
            // hPhiColGen->Fill(fabs(piPos.Phi() - piNeg.Phi()) - TMath::Pi());
            // int mult = piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRecV.XYvector().DeltaPhi(piNegRecV.XYvector()) < 0;
            // hDeltaMomVsDeltaPhi->Fill(fabs(piPos.Mag() / piNeg.Mag()) - fabs(piPosRecV.Mag() / piNegRecV.Mag()),
            //                         fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())) - fabs(piPosRecV.XYvector().DeltaPhi(piNegRecV.XYvector()))
            //                         - 0 * TMath::Power(-1, mult) * TMath::Pi() );
            // //  Cowboy Type
            // if(piPosRecV.Cross(field).XYvector().DeltaPhi(piNegRecV.XYvector()) < TMath::Pi() / 2)
            // { 
            //     hDeltaPhiVsPhiCowboy->Fill(piPos.XYvector().Phi(), piPos.XYvector().DeltaPhi(piPosRecV.XYvector()));
            //     hDeltaPhiRecVsGenCowboy->Fill(fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())), 
            //                                 fabs(piPosRecV.XYvector().DeltaPhi(piNegRecV.XYvector())) );

            //     // hPhi->Fill(piPos.Phi() - piPosRecV.Phi()); 
            //     // hPhi->Fill(piNeg.Phi() - piNegRecV.Phi()); 
            //     hDeltaMomVsDeltaPhiCowboy->Fill(fabs(piPos.Mag() - piNeg.Mag()) - fabs(piPosRecV.Mag() - piNegRecV.Mag()),
            //                         piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRecV.XYvector().DeltaPhi(piNegRecV.XYvector())
            //                         - TMath::Power(-1, mult) * TMath::Pi() );
            // }
            // //  Sailor Type
            // if(piPosRecV.Cross(field).XYvector().DeltaPhi(piNegRecV.XYvector()) > TMath::Pi() / 2)
            // {  
            //     hDeltaPhiVsPhiSailor->Fill(piPos.XYvector().Phi(), piPos.XYvector().DeltaPhi(piPosRecV.XYvector()));

            //     hDeltaPhiRecVsGenSailor->Fill(fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())), 
            //                                 fabs(piPosRecV.XYvector().DeltaPhi(piNegRecV.XYvector())) );
            //     hDeltaMomVsDeltaPhiSailor->Fill(fabs(piPos.Mag() - piNeg.Mag()) - fabs(piPosRecV.Mag() - piNegRecV.Mag()),
            //                         piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRecV.XYvector().DeltaPhi(piNegRecV.XYvector())
            //                         - TMath::Power(-1, mult) * TMath::Pi() );
            // }

            tNew->Fill();
        }
        NgoodTr = 0;
        nTracks.clear();
        nTracks.shrink_to_fit();
    }

    std::vector<TCanvas *> canvs;
    for(int i = 0; i < 3; i++)
    { 
        canvs.push_back(new TCanvas(("canv" + std::to_string(i + 1)).c_str(), ("deltaPhi vs Phi pi+ Part" + std::to_string(i + 1)).c_str() ) ); 
        canvs[canvs.size() - 1]->Divide(2, 3);
    }

    std::vector<double> deltaPhiPos;
    std::vector<double> deltaPhiNeg;
    std::vector<double> deltaPtotPos;
    std::vector<double> deltaPtotNeg;

    std::vector<double> deltaPhiPosErr;
    std::vector<double> deltaPhiNegErr;
    std::vector<double> deltaPtotPosErr;
    std::vector<double> deltaPtotNegErr;
    for(int i = 0; i < Nchambers; i++)
    {
        canvs[int(i / 6)]->cd(i - 6 * int(i / 6) + 1);
        vDeltaPhiVsPhiProfPos[i]->Draw();

        deltaPhiPos.push_back(vDeltaPhiVsPhiProfPos[i]->GetMean(2));
        deltaPhiNeg.push_back(vDeltaPhiVsPhiProfNeg[i]->GetMean(2));

        deltaPtotPos.push_back(vDeltaPtotVsPhiProfPos[i]->GetMean(2));
        deltaPtotNeg.push_back(vDeltaPtotVsPhiProfNeg[i]->GetMean(2));

        deltaPhiPosErr.push_back(vDeltaPhiVsPhiProfPos[i]->GetMeanError(2));
        deltaPhiNegErr.push_back(vDeltaPhiVsPhiProfNeg[i]->GetMeanError(2));

        deltaPtotPosErr.push_back(vDeltaPtotVsPhiProfPos[i]->GetMeanError(2));
        deltaPtotNegErr.push_back(vDeltaPtotVsPhiProfNeg[i]->GetMeanError(2));
    }
    for(auto canv : canvs)
    { canv->Write(); }

    
    printVector(deltaPhiPos, "deltaPhiPos");
    printVector(deltaPhiNeg, "deltaPhiNeg");
    printVector(deltaPtotPos, "deltaPtotPos");
    printVector(deltaPtotNeg, "deltaPtotNeg");

    printVector(deltaPhiPosErr, "deltaPhiPosErr");
    printVector(deltaPhiNegErr, "deltaPhiNegErr");
    printVector(deltaPtotPosErr, "deltaPtotPosErr");
    printVector(deltaPtotNegErr, "deltaPtotNegErr");

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
