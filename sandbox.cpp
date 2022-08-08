#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

int sandbox()
{   
    // TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/New folder/tr_ph_run2bGen_511.root");    
    TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/tr_ph_run2bGen510NoField.root");    
    auto tr1 = (TTree *)file->Get("tr_ph");

    Int_t nks;
    Int_t tcharge[50];
    Float_t piTheta[50];
    Float_t piPhi[50];    
    Int_t ksTracks[50][50];

    Float_t simmom[20];
    Float_t simphi[20];
    Float_t simtheta[20];
    Int_t nsim;
    Int_t simtype[20];
    Int_t simorig[20];
    

    tr1->SetBranchAddress("nks", &nks); 
    tr1->SetBranchAddress("tcharge", tcharge); 
    tr1->SetBranchAddress("kspiphi", piPhi); 
    tr1->SetBranchAddress("kspith", piTheta); 
    tr1->SetBranchAddress("ksvind", ksTracks); 
    tr1->SetBranchAddress("simphi", simphi); 
    tr1->SetBranchAddress("simtheta", simtheta); 
    tr1->SetBranchAddress("simtype", simtype); 
    tr1->SetBranchAddress("simorig", simorig); 
    tr1->SetBranchAddress("simmom", simmom); 
    tr1->SetBranchAddress("nsim", &nsim); 

    auto hDeltaPhi = new TH1D("hDeltaPhi", "#Delta#phi", 500, 0, 6.28);
    auto hDeltaPhiGen = new TH1D("hDeltaPhiGen", "#Delta#phi Gen", 250, 0, 6.28);
    auto hMomentumTotal = new TH1D("hMomentumTotal", "Total momentum", 250, 0, 20);

    double phi1 = 0; 
    double phi2 = 0; 
    int foo = 0;
    double dPhi = 0;

    TVector3 v1;
    TVector3 v2;

    TVector3 kl;
    TVector3 piPos;
    TVector3 piNeg;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;

    for(int i = 0; i < tr1->GetEntriesFast(); i++)
    {
        tr1->GetEntry(i);
        // std::cout << nks << std::endl;
        if(nks > 0)
        {
            if (tcharge[ksTracks[0][0]] > 0)
            { 
                v1.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[0]);
                v2.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[1]);
            }
            else
            { 
                v1.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[1]);
                v2.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[0]);
            }
            hDeltaPhi->Fill(v1.Angle(v2));

            // for(int j = 0; j < nsim; j++)
            // {
            //     if(simtype[j] == 211 && simorig[j] == 310)
            //     { 
            //         phi1 = simphi[j]; 
            //         v1.SetMagThetaPhi(1, TMath::Pi() / 2, simphi[j]);
            //     }

            //     if(simtype[j] == -211 && simorig[j] == 310)
            //     { 
            //         phi2 = simphi[j]; 
            //         v2.SetMagThetaPhi(1, TMath::Pi() / 2, phi2);
            //     }
            // }
            // hDeltaPhiGen->Fill(v1.Angle(v2));
        }

        for(int j = 0; j < nsim; j++)
        {
            if(simtype[j] == 211 && simorig[j] == 310)
            { 
                piPos.SetMagThetaPhi(simmom[j], simtheta[j], simphi[j]);
                phi1 = simphi[j]; 
                v1.SetMagThetaPhi(1, TMath::Pi() / 2, simphi[j]);
                count1 = 1;
            }

            if(simtype[j] == -211 && simorig[j] == 310)
            { 
                piNeg.SetMagThetaPhi(simmom[j], simtheta[j], simphi[j]);
                phi2 = simphi[j]; 
                v2.SetMagThetaPhi(1, TMath::Pi() / 2, phi2);
                count2 = 1;
            }

            if(simtype[j] == 130)
            { 
                kl.SetMagThetaPhi(simmom[j], simtheta[j], simphi[j]); 
                count3 = 1;
            }
        }
        // hDeltaPhiGen->Fill(v1.Angle(v2));
        hDeltaPhiGen->Fill(fabs(phi1 - phi2));
        if(count1 == 1 && count2 == 1 && count3 == 1)
        {
            hMomentumTotal->Fill((kl + piPos + piNeg).Mag());
            if((kl + piPos + piNeg).Mag() > 1e-4)
            {   std::cout << (kl + piPos + piNeg).Mag() << std::endl; }
        }
        count1 = 0;
        count2 = 0;
        count3 = 0;
    }
    hDeltaPhiGen->SetLineColor(kRed);
    // hDeltaPhi->Draw();
    hDeltaPhiGen->Draw("Same");
    hMomentumTotal->Draw();

    return 0;
}