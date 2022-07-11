#define KSKL_cxx
#include "KSKL.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>

void KSKL::Loop(std::string histFileName)
{
    //   In a ROOT session, you can do:
    //      root> .L KSKL.C
    //      root> KSKL t
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
    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
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

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);


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

        if(NgoodTr == 2)
        { NgoodTrS++; }

        if (NgoodTr == 2 && nks == 1 && is_coll != 1 && ksalign[0] > 0.85 && (tdedx[ksvind[0][0]] + tdedx[ksvind[0][1]]) / 2 < 5000 &&
            abs(kspith[0][0] - TMath::Pi() / 2) <= 0.9 && abs(kspith[0][1] - TMath::Pi() / 2) <= 0.9 &&
            //20 - half of the linear size of Drift Chamber
            //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
            kspipt[0][0] > 130 && kspipt[0][1] > 130 && 
            kspipt[0][0] < 320 && kspipt[0][1] < 320 &&
            ksminv[0] > 480 && ksminv[0] < 510 &&
            tcharge[ksvind[0][0]] * tcharge[ksvind[0][1]] < 0)
        {
            if (tcharge[ksvind[0][0]] > 0)
            {
                Y = kspipt[0][0] / kspipt[0][1];
                hist->Fill(kspipt[0][0], kspipt[0][1]);
            }
            else
            {
                Y = kspipt[0][1] / kspipt[0][0];
                hist->Fill(kspipt[0][1], kspipt[0][0]);
            }
            dpsi = ksdpsi[0];
            tNew->Fill();
        }

        NgoodTr = 0;
    }

    std::cout << "NgoodTrS = "<< NgoodTrS << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    hist->Draw("COL");

    top->Write();
    top->Save();
}
