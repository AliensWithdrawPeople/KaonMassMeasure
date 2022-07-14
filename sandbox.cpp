#include "TH2D.h"
#include "TH1D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TLine.h"
#include "TArrow.h"
#include "TText.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVector3.h"

int sandbox()
{
    //TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/scan2018_omphi_tr_ph_fc_e509.5_v8.root");    
    TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/tr_ph_kskl_2bodygen600k.root");    
    //TFile *file = TFile::Open("E:/Science/BINP/Kaon Mass Measure/tr_ph/tr_ph_KsKlsim_mcgpj800k.root");    
    auto tr1 = (TTree *)file->Get("tr_ph");

    Float_t ksTheta[10]; Float_t klTheta[30]; Float_t ksZ[30];
    Float_t ksPhi[10]; Float_t klPhi[30]; Float_t klRho[30];
    Float_t klEn[30];
    Int_t nks; Int_t nph;
    Float_t xbeam; Float_t ybeam;

    
    tr1->SetBranchAddress("ksth", ksTheta);
    tr1->SetBranchAddress("phth0", klTheta);
    tr1->SetBranchAddress("ksphi", ksPhi);
    tr1->SetBranchAddress("phphi0", klPhi);
    tr1->SetBranchAddress("phrho", klRho);
    tr1->SetBranchAddress("ksz0", ksZ);

    tr1->SetBranchAddress("phen0", klEn);

    tr1->SetBranchAddress("nks", &nks); 
    tr1->SetBranchAddress("nph", &nph);
    tr1->SetBranchAddress("xbeam", &xbeam);
    tr1->SetBranchAddress("ybeam", &ybeam);
    

    auto hist = new TH2D("hist", "", 1000, 0, 500, 1000, 0, 20000);
    auto hist1 = new TH2D("hist1", "", 600, 0, 2*TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hist2 = new TH2D("hist2", "", 600, -TMath::Pi(), TMath::Pi(), 600, 0, 600);
    auto hist3 = new TH1D("hist3", "", 350, 0, 350);
    auto hist4 = new TH1D("hist4", "", 628, 0, 6.28);
    auto hist5 = new TH1D("hist5", "", 2000, -1.1, 1.1);
    int n = 0;

    auto c1 = new TCanvas("c1", "c1", 200, 10, 800, 600);
    //tr1->Draw("p>>hist4", "", "goff");
    
    //n = tr1->Draw("phen0[0] >> hist3", "nph == 1 && phen[0] > 100 && nks == 1");

    // Spatial angle between Ks and Kl vectors of movements.
    double cosPsi = 0;
    TVector3 ks;
    TVector3 kl;
    double tmp = 0;
    int foo = 0;

    for(int i = 0; i < tr1->GetEntriesFast(); i++)
    {
        tr1->GetEntry(i);
        for(int k = 0; k < nks; k++)
        {
            for(int j = 0; j < nph; j++)
            {
                ks.SetMagThetaPhi(1, ksTheta[j], ksPhi[j]);
                kl.SetMagThetaPhi(klRho[j], klTheta[j], klPhi[j]);
                kl.SetX(kl.X() - xbeam);
                kl.SetY(kl.Y() - ybeam);
                kl.SetZ(kl.Z() - ksZ[k]);

                cosPsi = ks.Unit() * kl.Unit();
                
                /*
                psi = TMath::ACos(TMath::Sin(ksTheta[k]) * TMath::Cos(ksPhi[k]) * TMath::Sin(klTheta[j]) * TMath::Cos(klPhi[j]) + 
                                  TMath::Sin(ksTheta[k]) * TMath::Sin(ksPhi[k]) * TMath::Sin(klTheta[j]) * TMath::Sin(klPhi[j]) +
                                  TMath::Cos(ksTheta[k]) * TMath::Cos(klTheta[j]) );
                */

                if(fabs(cosPsi) > 0.8)
                {
                    for(int l = 0; l < nph; l++)
                    { hist3->Fill(klEn[l]); }
                }
                
                hist5->Fill(cosPsi);

                // if(z-component of the cross product of Ks and Kl vectors < 0)
                /*
                tmp = TMath::ATan2(TMath::Sin(ksTheta[k]) * TMath::Sin(ksPhi[k]) + TMath::Cos(klTheta[j]), TMath::Cos(ksTheta[k]) + TMath::Sin(klTheta[j]) * TMath::Sin(klPhi[j]));
                if(cos(tmp) * (TMath::Sin(ksTheta[k])*TMath::Cos(ksPhi[k]) * TMath::Sin(klTheta[j])*TMath::Sin(klPhi[j]) - 
                   TMath::Sin(ksTheta[k]) * TMath::Sin(ksPhi[k]) * TMath::Sin(klTheta[j])*TMath::Cos(klPhi[j])) + 
                   sin(tmp) * (TMath::Cos(ksTheta[k]) * TMath::Sin(klTheta[j]) * TMath::Cos(klPhi[j]) - TMath::Cos(klTheta[j]) * TMath::Sin(ksTheta[k]) * TMath::Cos(ksPhi[k])) < 0)
                { psi = 2 * TMath::Pi() - psi; }
                */
                
                // d? = |?1 - ?2| - ?  è d? =(?1 + ?2) - ?;
                foo = ks.Phi() - kl.Phi() < 0;
                if(fabs(ks.Phi() - kl.Phi()) <= TMath::Pi())
                { hist1->Fill(ks.Phi() - kl.Phi() + foo * 2 * TMath::Pi(), ks.Theta() + kl.Theta() - TMath::Pi()); }
                else
                { hist1->Fill(std::pow(-1, foo) * 2 * TMath::Pi() - (ks.Phi() - kl.Phi()) + foo * 2 * TMath::Pi(), ks.Theta() + kl.Theta() - TMath::Pi()); }


                //hist1->Fill(ks.Phi() - kl.Phi(), ks.Theta() + kl.Theta() - TMath::Pi());
                hist4->Fill(ks.Angle(kl));

                hist2->Fill(fabs(ks.Phi() - kl.Phi()) - TMath::Pi(), klEn[j]);
            }
        }
    }
    n = hist4->GetEntries();
    
    //hist3->SetTitle("Distribution of the enrgy of Kl candidates");
    hist1->GetXaxis()->SetTitle("#Delta#phi, rad");
    hist1->GetYaxis()->SetTitle("#Delta#theta, rad");

    hist2->GetXaxis()->SetTitle("#Delta#phi, rad");
    hist2->GetYaxis()->SetTitle("Cluster Energy,  MeV");

    hist1->Draw("COL");


/*    
n = tr->Draw("tptot[0]:tptot[1]>>hist", 
    "nt == 2 && is_coll!=1 && kstype[0] == 0 && tnhit[ksvind[0][0]]>10 && tnhit[ksvind[0][1]]>10 &&\
    (tdedx[ksvind[0][0]]+tdedx[ksvind[0][1]])/2 < 5000 && kstlen[0]<2 && ksalign[0]>0.85 &&\
    abs(tth[ksvind[0][0]] - TMath::Pi()/2) <= 0.4 && abs(tth[ksvind[0][1]] - TMath::Pi()/2) <= 0.4", "COL");

    "nt == 2 && is_coll!=1 && kstype[0] == 0 && tnhit[ksvind[0][0]]>10 && tnhit[ksvind[0][1]]>10 &&\
    (tdedx[ksvind[0][0]]+tdedx[ksvind[0][1]])/2 < 5000 && kstlen[0]>0.1 && kstlen[0]<1.3 && ksalign[0]>0.85 &&\
    abs(tth[ksvind[0][0]] - TMath::Pi()/2) <= 0.5 && abs(tth[ksvind[0][1]] - TMath::Pi()/2) <= 0.5 &&\
    tcharge[ksvind[0][0]]*tcharge[ksvind[0][1]] < 0 && fabs(tz[ksvind[0][0]]) < 10 && fabs(tz[ksvind[0][1]]) < 10"

    tr->Draw("(tdedx[0] + tdedx[1])/2 : (tptot[0] + tptot[1])/2 >>hist2", 
    "nt == 2 && is_coll==1 && tnhit[0]>10 && tnhit[1]>10 && tcharge[0]*tcharge[1] < 0 &&\
    (tth[0] - tth[1] + TMath::Pi())/2 < TMath::Pi() && (tth[0] - tth[1] + TMath::Pi())/2 > 1 &&\
    abs(tth[0] - TMath::Pi()/2) <= 0.4 && abs(tth[1] - TMath::Pi()/2) <= 0.4  && (tdedx[0]+tdedx[1])/2 > 7000 && \
    fabs(tptot[0]-tptot[1])/(tptot[0]+tptot[1]) < 0.15 && tptot[0] > 70 && tptot[1] > 70", "COL");
*/

/*
    n = tr->Draw("fabs(tptot[0]-tptot[1])/(tptot[0]+tptot[1])>>hist3", 
    "nt == 2 && is_coll==1 && tnhit[0]>10 && tnhit[1]>10 && tnhit[0]<21 && tnhit[1]<21 && tcharge[0]*tcharge[1] < 0 &&\
    TMath::Abs(trho[0]) < 0.3 && TMath::Abs(trho[1]) < 0.3 &&\
    (tth[0] - tth[1] + TMath::Pi())/2 < TMath::Pi() && (tth[0] - tth[1] + TMath::Pi())/2 > 1 &&\
    abs(tth[0] - TMath::Pi()/2) <= 0.4 && abs(tth[1] - TMath::Pi()/2) <= 0.4  && (tdedx[0]+tdedx[1])/2 > 7000 &&\
    fabs(tptot[0]-tptot[1])/(tptot[0]+tptot[1]) < 0.15", "COL");
*/

    //tr->Draw("(tdedx[0]+tdedx[1])/2 : (tptot[0]+tptot[1])/2>> hist2", "nt == 2 && is_coll!=1 && tnhit[0]>10 && tnhit[1]>10 && nks==1 && kstlen[0]>0.1 && kstlen[0]<1.3", "COL");
    
    return n;
}