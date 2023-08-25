#define ksklTrEff_cxx
#include "ksklTrEff.hpp"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLine.h>


void ksklTrEff::Loop(std::string outFileName)
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    TFile *top = new TFile(outFileName.c_str(), "recreate");
    auto tNew = new TTree("tr_eff", "Cutted tr_ph (track efficiency KsKl related data)");

    TLorentzVector Kl(1, 1, 1, 1);
    TLorentzVector Tr1(1, 1, 1, 1);
    TLorentzVector Tr2(1, 1, 1, 1);
    TLorentzVector Lost(1, 1, 1, 1);
    TLorentzVector Tot(TVector3(0, 0, 0), 505);
    std::vector<int> trackNumber = {};
    std::vector<double> recoilMass_Kl = {};
    std::vector<double> recoilMass_2tracks = {};

    int trackCounter = 0;
    int clusterCounter = 0;
    int isKsKl = 0;
    tNew->Branch("tracks", &trackCounter, "tracks/I");
    tNew->Branch("clusters", &clusterCounter, "clusters/I");
    tNew->Branch("isKsKl", &isKsKl, "isKsKl/I");


    auto hClEnergy = new TH1D("hClEnergy", "hClEnergy", 1000, 0, 1000);
    hClEnergy->GetXaxis()->SetTitle("E_{no track cluster}, MeV");

    auto hTracks = new TH1D("hTracks", "hTracks", 5, 0, 5);
    auto hSimType = new TH1D("hSimType", "hSimType", 2000, -1000, 1000);
    auto hKsMass_Kl = new TH1D("hKsMass_Kl", "hKsMass_Kl", 4000, 0, 1000);
    auto hKchMass = new TH1D("hKchMass", "hKchMass", 20000, 0, 1000);
    auto hdEdx = new TH2D("hdEdx", "hdEdx", 1000, 0, 1000, 3e4, 0, 3e4);

    auto hMasses = new TH2D("hMasses", "hMasses", 4000, -1e6, -1e6, 4000, -1e6, -1e6);

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if(nt < 1 || nt > 2)
        { continue; }

        if(nt == 2 && is_coll == 1)
        { continue; }

        for(int i = 0; i < nt; i++)
        {
            if (fabs(trho[i]) < 6 && fabs(tz[i]) < 10 &&
                tchi2r[i] < 15 && tchi2z[i] < 10 && tnhit[i] > 6 && 
                0.2 < fabs(trho[i]) && 
                // tdedx[i] < 3e3 && 
                170 < tptot[i] && tptot[i] < 320 &&
                1.1 < tth[i] && tth[i] < TMath::Pi() - 1.1 )
            { 
                trackNumber.push_back(i);
                trackCounter++; 
            }
        }

        if(trackCounter > 0)
        {
            Tot.SetE(2 * emeas);
            TVector3 pi1(1, 1, 1);
            TVector3 pi2(1, 1, 1);
            pi1.SetMagThetaPhi(tptot[trackNumber[0]], tth[trackNumber[0]], tphi[trackNumber[0]]);
            Tr1.SetVectM(pi1, 139.570);
            if(trackCounter == 2)
            { 
                pi2.SetMagThetaPhi(tptot[trackNumber[1]], tth[trackNumber[1]], tphi[trackNumber[1]]);
                Tr2.SetVectM(pi2, 139.570); 
            }
            
            for(int i = 0; i < nph; i++)
            {
                if(isKsKl == 1)
                { break; }

                hClEnergy->Fill(phen0[i]);
                if(phen0[i] > 40 && 1.1 < phth0[i] && phth0[i] < TMath::Pi() - 1.1)
                { 
                    clusterCounter++; 
                    TVector3 kl(1, 1, 1);
                    kl.SetMagThetaPhi(sqrt(emeas * emeas - 497.614 * 497.614), phth0[i], phphi0[i]);
                    Kl.SetVect(kl);
                    Kl.SetE(emeas);

                    recoilMass_Kl.push_back((Tot - Kl).M2());
                    hKsMass_Kl->Fill((Tot - Kl - Tr1).M());
                    if(trackCounter == 2)
                    { hMasses->Fill((Tot - Kl - Tr1).M2(), (Tot - Kl - Tr2).M2()); }

                    auto recoilMass = (Tot - Kl - Tr1).M2();
                    // Tr1.SetVect(Tr1.Vect());
                    // Tr1.SetE(emeas);
                    // hKchMass->Fill(Tr1.M());
                    // if(isKsKl == 0 && fabs(Tr1.M() - 493.677) > 20 && recoilMass > 0 && 120 < sqrt(recoilMass) && sqrt(recoilMass) < 150)
                    // { 
                    //     isKsKl = 1;
                    //     hTracks->Fill(trackCounter);
                    //     for(int j = 0; j < nsim; j++)
                    //     { 
                    //         if(simorig[j] == 0)
                    //         { hSimType->Fill(simtype[j]); }
                    //     }
                    // }
                }
            }

            Tr1.SetVect(Tr1.Vect());
            Tr1.SetE(emeas);
            hKchMass->Fill(Tr1.M());

            if(isKsKl == 0 && fabs(Tr1.M() - 493.677) > 20)
            { 
                isKsKl = 1;
                hTracks->Fill(trackCounter);
                hdEdx->Fill(tptot[trackNumber[0]], tdedx[trackNumber[0]]);
                for(int j = 0; j < nsim; j++)
                { 
                    if(simorig[j] == 0)
                    { hSimType->Fill(simtype[j]); }
                }
            }
        }

        tNew->Fill();
        trackCounter = 0;
        clusterCounter = 0;
        isKsKl = 0;
        trackNumber.clear();
    }

    top->Write();
    top->Save();
}