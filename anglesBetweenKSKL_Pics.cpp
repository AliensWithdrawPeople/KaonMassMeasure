#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

int anglesBetweenKSKL_Pics()
{
    //TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/scan2018_omphi_tr_ph_fc_e509.5_v8.root");    
    //TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/tr_ph_kskl_2bodygen600k.root");    
    TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/tr_ph_KsKlsim_mcgpj800k.root");    
    auto tr1 = (TTree *)file->Get("tr_ph");

    Float_t ksTheta[10]; Float_t klTheta[30]; Float_t ksZ[30];
    Float_t ksPhi[10]; Float_t klPhi[30]; Float_t klRho[30];
    Float_t klEn[30];
    Int_t nks; Int_t nph;
    Float_t xbeam; Float_t ybeam;

    tr1->SetBranchAddress("ksz0", ksZ);
    tr1->SetBranchAddress("ksth", ksTheta);
    tr1->SetBranchAddress("ksphi", ksPhi);
    tr1->SetBranchAddress("phth0", klTheta);
    tr1->SetBranchAddress("phphi0", klPhi);
    tr1->SetBranchAddress("phrho", klRho);
    tr1->SetBranchAddress("phen0", klEn);

    tr1->SetBranchAddress("nks", &nks); 
    tr1->SetBranchAddress("nph", &nph);
    tr1->SetBranchAddress("xbeam", &xbeam);
    tr1->SetBranchAddress("ybeam", &ybeam);

    
    // Theta = kl.Theta()
    auto hdPhiTheta = new TH2D("hdPhiTheta", "", 600, 0, TMath::Pi(), 600, 0, 2*TMath::Pi());
    auto hdThetadPhi = new TH2D("hdThetadPhi", "", 600, 0, 2*TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hClEdPhi = new TH2D("hClEdPhi", "", 600, 0, 2*TMath::Pi(), 600, 0, 600);
    auto hClEnergy = new TH1D("hClEnergy", "", 350, 0, 350);
    auto hPsi = new TH1D("hPsi", "", 628, 0, 6.28);
    auto hCosPsi = new TH1D("hCosPsi", "", 2000, -1.1, 1.1);
    int n = 0;

    auto c1 = new TCanvas("c1", "c1", 200, 10, 800, 600);    
    //n = tr1->Draw("phen0[0] >> hClEnergy", "nph == 1 && phen[0] > 100 && nks == 1");

    // Spatial angle between Ks and Kl vectors of movements.
    double cosPsi = 0;
    TVector3 ks;
    TVector3 kl;

    TVector3 ksPhiVec;
    TVector3 klPhiVec;
    double rotAngle = 0;
    int foo = 0;
    double dPhi = 0;

    for(int i = 0; i < tr1->GetEntriesFast(); i++)
    {
        tr1->GetEntry(i);
        
        for(int k = 0; k < nks; k++)
        {
            ks.SetMagThetaPhi(1, ksTheta[k], ksPhi[k]);
            for(int j = 0; j < nph; j++)
            {

                kl.SetMagThetaPhi(klRho[j], klTheta[j], klPhi[j]);
                kl.SetX(kl.X() - xbeam);
                kl.SetY(kl.Y() - ybeam);
                kl.SetZ(kl.Z() - ksZ[k]);

                ksPhiVec.SetMagThetaPhi(1, TMath::Pi() / 2, ks.Phi());
                klPhiVec.SetMagThetaPhi(1, TMath::Pi() / 2, kl.Phi());

                cosPsi = ks.Unit() * kl.Unit();

                if(fabs(cosPsi) > 0.8)
                {
                    for(int l = 0; l < nph; l++)
                    { hClEnergy->Fill(klEn[l]); }
                }
                
                hCosPsi->Fill(cosPsi);
                
                foo = ks.Phi() - kl.Phi() < 0;
                if(fabs(ks.Phi() - kl.Phi()) <= TMath::Pi())
                { 
                    dPhi = ks.Phi() - kl.Phi() + foo * 2 * TMath::Pi(); // wrong
                    dPhi = ksPhiVec.Angle(klPhiVec); // correct
                    hdThetadPhi->Fill(dPhi, ks.Theta() + kl.Theta() - TMath::Pi());
                    hClEdPhi->Fill(dPhi, klEn[j]); 
                    hdPhiTheta->Fill(kl.Theta(), dPhi);
                }
                else
                { 
                    dPhi = std::pow(-1, foo) * 2 * TMath::Pi() - (ks.Phi() - kl.Phi()) + foo * 2 * TMath::Pi(); // wrong
                    dPhi = ksPhiVec.Angle(klPhiVec); // correct
                    hdThetadPhi->Fill(dPhi, ks.Theta() + kl.Theta() - TMath::Pi()); 
                    hClEdPhi->Fill(dPhi, klEn[j]);
                    hdPhiTheta->Fill(kl.Theta(), dPhi);
                }

                hPsi->Fill(ks.Angle(kl));
            }
        }
    }
    n = hPsi->GetEntries();
    
    hClEnergy->SetTitle("Distribution of the enrgy of Kl candidates");
    hdThetadPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hdThetadPhi->GetYaxis()->SetTitle("#Delta#theta, rad");

    hdPhiTheta->GetXaxis()->SetTitle("#theta of Kl, rad");
    hdPhiTheta->GetYaxis()->SetTitle("#Delta#phi, rad");

    hPsi->GetXaxis()->SetTitle("Angle between motion vectors Ks and Kl, rad");

    hClEdPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hClEdPhi->GetYaxis()->SetTitle("Cluster Energy,  MeV");

    hClEnergy->Draw("COL");

    return n;
}