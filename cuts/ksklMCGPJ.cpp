#define ksklMCGPJ_cxx
#include "ksklMCGPJ.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>

void ksklMCGPJ::Loop(std::string histFileName)
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
    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", &dpsi, "dpsi/F");
    tNew->Branch("pi+ mom", &dpsi, "dpsi/F");
    tNew->Branch("pi- mom", &dpsi, "dpsi/F");
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
    int tmpCounter = 0;
    std::vector<Int_t> ksCand = {};

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);
    auto histKlCands = new TH1D("KlCands", "number of Kl candidates for one Ks candidate", 5, 0, 5);
    auto histKsCands = new TH1D("KsCands", "number of Ks candidates", 5, 0, 5);

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
        
        /*
        if (NgoodTr == 2 && nks == 1 && is_coll != 1 && ksalign[0] > 0.85 && (tdedx[ksvind[0][0]] + tdedx[ksvind[0][1]]) / 2 < 5000 &&
            abs(kspith[0][0] - TMath::Pi() / 2) <= 0.9 && abs(kspith[0][1] - TMath::Pi() / 2) <= 0.9 &&
            //20 - half of the linear size of Drift Chamber
            //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
            kspipt[0][0] > 120 && kspipt[0][1] > 120 && 
            kspipt[0][0] < 350 && kspipt[0][1] < 350 &&
            
            tcharge[ksvind[0][0]] * tcharge[ksvind[0][1]] < 0)
        {
            if (tcharge[ksvind[0][0]] > 0)
            {
                Y = kspipt[0][0] / kspipt[0][1];
                hist->Fill(kspipt[0][0], kspipt[0][1]);
                p1 = kspipt[0][0]; p2 = kspipt[0][1];
            }
            else
            {
                Y = kspipt[0][1] / kspipt[0][0];
                hist->Fill(kspipt[0][1], kspipt[0][0]);
                p1 = kspipt[0][1]; p2 = kspipt[0][0];
            }
            dpsi = ksdpsi[0];
            tNew->Fill();
        }
        */

        if (NgoodTr == 2 && is_coll != 1 )
        {
            for(int k = 0; k < nks; k++)
            {
               if(ksalign[k] > 0.85 && (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000 &&
                  abs(kspith[k][0] - TMath::Pi() / 2) <= 0.9 && abs(kspith[k][1] - TMath::Pi() / 2) <= 0.9 &&
                  //20 - half of the linear size of Drift Chamber
                  //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
                  kspipt[k][0] > 120 && kspipt[k][1] > 120 && 
                  kspipt[k][0] < 350 && kspipt[k][1] < 350 &&
                  tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0)
               {
                  for(int j = 0; j < nph; j++)
                  {
                     // phth0 and phphi0 are theta and phi respectively angles of KL candidate ,
                     cosKlKs =   TMath::Sin(ksth[k])*TMath::Cos(ksphi[k]) * TMath::Sin(phth0[j])*TMath::Cos(phphi0[j]) + 
                                 TMath::Sin(ksth[k]) * TMath::Sin(ksphi[k]) * TMath::Sin(phth0[j])*TMath::Sin(phphi0[j]) +
                                 TMath::Cos(ksth[k]) * TMath::Cos(phth0[j]);

                    if(fabs(cosKlKs) > 0.95)
                     { tmpCounter++; } 
                  }

                  
      
                  histKlCands->Fill(tmpCounter);
                  if(tmpCounter > 0)
                  { ksCand.push_back(k); }
                  tmpCounter = 0;
               }
            }

            histKsCands->Fill(ksCand.size());
            if(ksCand.size() > 0)
            {
               if (tcharge[ksvind[0][0]] > 0)
               {
                  Y = kspipt[0][0] / kspipt[0][1];
                  hist->Fill(kspipt[0][0], kspipt[0][1]);
                  p1 = kspipt[0][0]; p2 = kspipt[0][1];
               }
               else
               {
                  Y = kspipt[0][1] / kspipt[0][0];
                  hist->Fill(kspipt[0][1], kspipt[0][0]);
                  p1 = kspipt[0][1]; p2 = kspipt[0][0];
               }
               dpsi = ksdpsi[0];
               tNew->Fill();
            }
            ksCand.clear();
            ksCand.shrink_to_fit();
        }

        NgoodTr = 0;
    }

    eff = (double) tNew->GetEntriesFast() / nentries;
    std::cout << "nentries " << nentries << std::endl;
    std::cout << "e_mc = " <<  eff << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    hist->Draw("COL");

    top->Write();
    top->Save();
}
