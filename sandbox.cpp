#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TRandom.h"

int sandbox()
{   
    TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/New folder/mcgpj/tr_ph_run000013.root");    
    // TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/tr_ph_run2bGen510NoField.root");    
    auto tr1 = (TTree *)file->Get("tr_ph");

    Int_t nks;
    Int_t tcharge[50];
    Float_t ksdpsi[50];
    Float_t piMom[50][2];
    Float_t piTheta[50][2];
    Float_t piPhi[50][2];    
    Int_t ksTracks[50][50];

    Float_t simmom[20];
    Float_t simphi[20];
    Float_t simtheta[20];
    Int_t nsim;
    Int_t simtype[20];
    Int_t simorig[20];
    

    tr1->SetBranchAddress("nks", &nks); 
    tr1->SetBranchAddress("tcharge", tcharge); 
    tr1->SetBranchAddress("ksdpsi", ksdpsi); 
    tr1->SetBranchAddress("kspipt", piMom); 
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
    auto hDeltaThetaGen = new TH1D("hDeltaThetaGen", "#Delta#theta Gen", 100, -3.14, 3.14);
    auto hThetaKsGen = new TH1D("hThetaKsGen", "#theta Ks Gen", 100, 0, 3.14);
    auto hMomentumTotal = new TH1D("hMomentumTotal", "Total momentum", 250, 0, 20);
    auto hDeltaThetaVsDeltaPhiGen = new TH2D("hDeltaThetaVsDeltaPhiGen", "", 100, -3.14, 3.14, 100, 0, 6.28);
    auto hLorentzPhi = new TH1D("hLorentzPhi", "Lorentz #phi", 500, 6.28, 6.28);
    auto hPhi = new TH1D("hPhi", "Lorentz #phi", 4000, -0.2, 0.2);
    auto hDeltaMom = new TH1D("hDeltaMom", "Lorentz delta mom", 4000, -50, 50);

    double phi1 = 0; 
    double phi2 = 0; 
    double theta1 = 0; 
    double theta2 = 0; 
    int foo = 0;
    double dPhi = 0;

    TVector3 v1;
    TVector3 v2;

    TVector3 kl;
    TVector3 piPos;
    TVector3 piNeg;
    TVector3 piPosRec(0, 0, 0);
    TVector3 field(0., 0., 1.);
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
                v1.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[0][0]);
                v2.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[0][1]);
            }
            else
            { 
                v1.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[0][1]);
                v2.SetMagThetaPhi(1, TMath::Pi() / 2, piPhi[0][0]);
            }

            if(gRandom->Uniform(1) > 0.5)
            { 
                hDeltaPhi->Fill(piPhi[0][0]  - piPhi[0][1]  < 0? piPhi[0][0]  - piPhi[0][1] + 2*TMath::Pi() : piPhi[0][0]  - piPhi[0][1]);
            }
            else
            { 
                hDeltaPhi->Fill(piPhi[0][1]  - piPhi[0][0]  < 0? piPhi[0][1]  - piPhi[0][0] + 2*TMath::Pi() : piPhi[0][1]  - piPhi[0][0]);
            }
        


            for(int j = 0; j < nsim; j++)
            {
                if(simtype[j] == 211 && simorig[j] == 310)
                { 
                    piPos.SetMagThetaPhi(simmom[j], simtheta[j], simphi[j]);
                    phi1 = simphi[j]; 
                    theta1 = simtheta[j];
                    v1.SetMagThetaPhi(1, TMath::Pi() / 2, simphi[j]);
                    count1 = 1;
                }

                if(simtype[j] == -211 && simorig[j] == 310)
                { 
                    piNeg.SetMagThetaPhi(simmom[j], simtheta[j], simphi[j]);
                    phi2 = simphi[j]; 
                    theta2 = simtheta[j];
                    v2.SetMagThetaPhi(1, TMath::Pi() / 2, phi2);
                    count2 = 1;
                }

                if(simtype[j] == 130)
                { 
                    kl.SetMagThetaPhi(simmom[j], simtheta[j], simphi[j]); 
                    count3 = 1;
                }

                if(simtype[j] == 310)
                { 
                    hThetaKsGen->Fill(simtheta[j]);
                }
            }

            if(gRandom->Uniform(1) > 0.5)
            { 
                hDeltaThetaGen->Fill(theta1 + theta2 - TMath::Pi());
                // hDeltaPhiGen->Fill(std::fmod(2*TMath::Pi() - phi1 + phi2, 2*TMath::Pi()));
                hDeltaPhiGen->Fill(phi2  - phi1 < 0? phi2  - phi1 + 2*TMath::Pi() : phi2  - phi1);
            }
            else
            { 
                hDeltaThetaGen->Fill(theta2 + theta1 - TMath::Pi());
                // hDeltaPhiGen->Fill(std::fmod(2*TMath::Pi() - phi2 + phi1, 2*TMath::Pi()) );
                hDeltaPhiGen->Fill(phi1  - phi2 < 0? phi1  - phi2 + 2*TMath::Pi() : phi1  - phi2);
            }
        }
        // hDeltaPhiGen->Fill(v1.Angle(v2));

        

        hDeltaThetaVsDeltaPhiGen->Fill(theta2 + theta1 - TMath::Pi(), fabs(phi1 - phi2));


        if(count1 == 1 && count2 == 1 && count3 == 1)
        {
            hMomentumTotal->Fill((kl + piPos + piNeg).Mag());
        }

        
        if(piNeg.Phi() < -TMath::Pi() / 2 && piPos.Phi() < TMath::Pi() / 2)
        { hLorentzPhi->Fill(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector())); }


        if (tcharge[ksTracks[0][0]] > 0)
        { piPosRec.SetMagThetaPhi(piMom[0][0], TMath::Pi() / 2, piPhi[0][0]); }
        else
        { piPosRec.SetMagThetaPhi(piMom[0][1], TMath::Pi() / 2, piPhi[0][1]); }

        if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) > TMath::Pi() / 2)
        { 
            hPhi->Fill(piPos.Phi() - piPosRec.Phi()); 
            hDeltaMom->Fill(piPos.Mag() - piPosRec.Mag());
        }

        count1 = 0;
        count2 = 0;
        count3 = 0;
    }
    hDeltaPhiGen->SetLineColor(kRed);

    hDeltaMom->Draw();
    // hPhi->Draw();
    // hLorentzPhi->Draw();

    // hDeltaPhiGen->Draw();
    // hDeltaPhi->Draw("Same");

    // hThetaKsGen->Draw();
    // hDeltaThetaGen->Draw();
    // hMomentumTotal->Draw();
    // hDeltaThetaVsDeltaPhiGen->Draw("COL");

    return 0;
}