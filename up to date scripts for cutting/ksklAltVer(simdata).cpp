#define ksklAltVer_cxx
#include "ksklAltVer.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>
#include "Math/Vector4D.h"

void ksklAltVer::Loop(std::string histFileName)
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

    TFile *top = new TFile(histFileName.c_str(), "recreate");
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    Float_t dpsi;
    Float_t Y;
    Float_t eff;
    Float_t p1; Float_t p2;
    tNew->Branch("emeas", &emeas, "emeas/F");
    tNew->Branch("demeas", &demeas, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", &dpsi, "dpsi/F");
    tNew->Branch("Y", &Y, "Y/F");

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

    double cosKlKs = 0;
    TVector3 ks;
    TVector3 kl;
    // Auxilary var (for dPhi).
    int foo = 0;
    // dPhi - angle between Ks and Kl
    double dPhi = 0;
    int tmpCounter = 0;
    std::vector<Int_t> ksCand = {};

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);
    auto histKlCands = new TH1D("KlCands", "number of Kl candidates for one Ks candidate", 5, 0, 5);
    auto histKsCands = new TH1D("KsCands", "number of Ks candidates", 5, 0, 5);
    // Theta = kl.Theta()
    auto hdPhiTheta = new TH2D("hdPhiTheta", "dPhi between Ks and Kl vs Theta  Kl", 600, 0, TMath::Pi(), 600, 0, 2*TMath::Pi());
    auto hdThetadPhi = new TH2D("hdThetadPhi", "dTheta vs dPhi between Ks and Kl", 600, 0, 2*TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hClEdPhi = new TH2D("hClEdPhi", "Cluster Energy vs dPhi", 600, 0, 2 * TMath::Pi(), 600, 0, 600);
    auto hPsi = new TH1D("hPsi", "", 628, 0, 6.28);

    hdThetadPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hdThetadPhi->GetYaxis()->SetTitle("#Delta#theta, rad");

    hdPhiTheta->GetXaxis()->SetTitle("#theta of Kl, rad");
    hdPhiTheta->GetYaxis()->SetTitle("#Delta#phi, rad");

    hPsi->GetXaxis()->SetTitle("Angle between motion vectors Ks and Kl, rad");

    hClEdPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hClEdPhi->GetYaxis()->SetTitle("Cluster Energy,  MeV");

    TVector3 piPos;
    TVector3 piNeg;
    double mass = 0;
    double massDif = 0;
    int foobarbaz = 0;
    int foobarbazDif = 0;
    int counter1 = 0;
    int counter2 = 0;
    ROOT::Math::PxPyPzEVector bar4Vec(0, 0, 0, emeas);
    auto hMass = new TH1D("hMass", "Mass", 1000, 495, 505);
    auto hMomentum = new TH1D("hMomentum", "Momentum", 1000, 50, 200);

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
        
        // if (NgoodTr == 2 && is_coll != 1 )
        // {
        //     for(int k = 0; k < nks; k++)
        //     {
        //         if(ksalign[k] > 0.85 && (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000 &&
        //             abs(kspith[k][0] - TMath::Pi() / 2) <= 0.7 && 
        //             abs(kspith[k][1] - TMath::Pi() / 2) <= 0.7 &&
        //             //20 - half of the linear size of Drift Chamber
        //             //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
        //             //kspipt[k][0] > 120 && kspipt[k][1] > 120 && 
        //             //kspipt[k][0] < 350 && kspipt[k][1] < 350 &&
        //             tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0) // Added kstype[k] == 0.
        //         {
        //             ks.SetMagThetaPhi(1, ksth[k], ksphi[k]);

        //             for(int j = 0; j < nph; j++)
        //             {
        //                 // phth0 and phphi0 are theta and phi respectively angles of KL candidate.
        //                 kl.SetMagThetaPhi(phrho[j], phth0[j], phphi0[j]);
        //                 // Shift to the center of phi decay. 
        //                 kl.SetX(kl.X() - xbeam);
        //                 kl.SetY(kl.Y() - ybeam);
        //                 kl.SetZ(kl.Z() - ksz0[k]);

        //                 hPsi->Fill(ks.Angle(kl));

        //                 foo = ks.Phi() - kl.Phi() < 0;
        //                 if(fabs(ks.Phi() - kl.Phi()) <= TMath::Pi())
        //                 { 
        //                     dPhi = ks.Phi() - kl.Phi() + foo * 2 * TMath::Pi();
        //                     hdThetadPhi->Fill(dPhi, ks.Theta() + kl.Theta() - TMath::Pi());
        //                     // phen0 - cluster energy. So this is the energy deposition of Kl candidate.
        //                     hClEdPhi->Fill(dPhi, phen0[j]); 
        //                     hdPhiTheta->Fill(kl.Theta(), dPhi);
        //                 }
        //                 else
        //                 { 
        //                     dPhi = std::pow(-1, foo) * 2 * TMath::Pi() - (ks.Phi() - kl.Phi()) + foo * 2 * TMath::Pi();
        //                     hdThetadPhi->Fill(dPhi, ks.Theta() + kl.Theta() - TMath::Pi()); 
        //                     hClEdPhi->Fill(dPhi, phen0[j]);
        //                     hdPhiTheta->Fill(kl.Theta(), dPhi);
        //                 }

        //                 if(fabs(dPhi - TMath::Pi()) < 0.5 && phen0[j] > 40)
        //                 { tmpCounter++; } 
        //             }        
        
        //             histKlCands->Fill(tmpCounter);
        //             if(tmpCounter > 0)
        //             { ksCand.push_back(k); }
        //             tmpCounter = 0;
        //         }
        //     }

        //     histKsCands->Fill(ksCand.size());
        //     if(ksCand.size() > 0)
        //     {
        //         if (tcharge[ksvind[0][0]] > 0)
        //         {
        //             Y = kspipt[0][0] / kspipt[0][1];
        //             hist->Fill(kspipt[0][0], kspipt[0][1]);
        //             p1 = kspipt[0][0]; p2 = kspipt[0][1];
        //         }
        //         else
        //         {
        //             Y = kspipt[0][1] / kspipt[0][0];
        //             hist->Fill(kspipt[0][1], kspipt[0][0]);
        //             p1 = kspipt[0][1]; p2 = kspipt[0][0];
        //         }

                
        //         for(int i = 0; i < nsim; i++)
        //         {
        //             if(simtype[i] == 211 && simorig[i] == 310)
        //             { 
        //                 piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]);

        //             }

        //             if(simtype[i] == -211 && simorig[i] == 310)
        //             { 
        //                 piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); 
                        
        //             }

        //             if(simtype[i] == 310)
        //             {
        //                 bar4Vec = ROOT::Math::PxPyPzEVector(simmom[i] * sin(simtheta[i]) * cos(simphi[i]), simmom[i] * sin(simtheta[i]) * sin(simphi[i]), 
        //                                                     simmom[i] * cos(simtheta[i]), emeas);
                        
        //                 hMass->Fill(bar4Vec.M());
        //                 hMomentum->Fill(simmom[i]);
        //             }
        //         }
        //         dpsi = piPos.Angle(piNeg);
        //         Y = piPos.Mag() / piNeg.Mag();
        //         tNew->Fill();
        //     }
        //     ksCand.clear();
        //     ksCand.shrink_to_fit();
        // }
        
        for(int i = 0; i < nsim; i++)
        {
            if(simtype[i] == 211 && simorig[i] == 310)
            { 
                piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); 
                counter1++;
            }

            if(simtype[i] == -211 && simorig[i] == 310)
            { 
                piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); 
                counter2++;
            }

            if(simtype[i] == 310)
            {
                bar4Vec = ROOT::Math::PxPyPzEVector(simmom[i] * sin(simtheta[i]) * cos(simphi[i]), simmom[i] * sin(simtheta[i]) * sin(simphi[i]), 
                                                    simmom[i] * cos(simtheta[i]), emeas);
                hMass->Fill(bar4Vec.M());
                mass += bar4Vec.M();
                foobarbaz++;
                hMomentum->Fill(simmom[i]);
            }
        }
        if(counter1 == 1 && counter2 == 1)
        {
            bar4Vec.SetE(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));
            massDif += bar4Vec.M();
            foobarbazDif++;
            dpsi = piPos.Angle(piNeg);
            Y = piPos.Mag() / piNeg.Mag();
            emeas = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()); 
            tNew->Fill();
        }
        counter1 = 0;
        counter2 = 0;

        NgoodTr = 0;
    }

    eff = (double) tNew->GetEntriesFast() / nentries;
    std::cout << "nentries " << nentries << std::endl;
    std::cout << "e_mc = " <<  eff << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;
    std::cout << "M_Ks = " << mass / foobarbaz << std::endl;
    std::cout << "M_Ks Alt = " << massDif / foobarbazDif << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    hist->Draw("COL");
    hMass->Draw("COL");
    hMomentum->Draw();

    top->Write();
    top->Save();
}