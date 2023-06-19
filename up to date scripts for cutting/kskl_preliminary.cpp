#define kskl_preliminary_cxx
#include "kskl_preliminary.h"
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

void kskl_preliminary::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      root> .L kskl_preliminary.C
//      root> kskl_preliminary t
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
    auto tNew = new TTree("tr_ph_cutted", "Preliminary cut tree");
    int NgoodTr = 0;
    int NgoodTrS = 0;

    tNew->Branch("NgoodTr", &NgoodTr, "NgoodTr/F");
    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("xbeam", &xbeam, "xbeam/F");
    tNew->Branch("ybeam", &ybeam, "ybeam/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", ksdpsi, "ksdpsi[50]/F");
    tNew->Branch("ksvind", ksvind, "ksvind[50][50]/F");
    tNew->Branch("tptot", tptotv, "tptot[50]/F");
    tNew->Branch("tdedx", tdedx, "tdedx[50]/F");
    tNew->Branch("ksalign", ksalign, "ksalign[50]/F");
    tNew->Branch("kspith", kspith, "kspith[50][50]/F");
    tNew->Branch("kstlen", kstlen, "kstlen[50]/F");
    tNew->Branch("tcharge", tcharge, "tcharge[50]/F");
    tNew->Branch("kstype", kstype, "kstype[50]/F");
    tNew->Branch("ksth", ksth, "ksth[50]/F");
    tNew->Branch("ksphi", ksphi, "ksphi[50]/F");
    tNew->Branch("phrho", phrho, "phrho[50]/F");
    tNew->Branch("phth0", phth0, "phth0[50]/F");
    tNew->Branch("phphi0", phphi0, "phphi0[50]/F");
    tNew->Branch("phen0", phen0, "phen0[50]/F");
    tNew->Branch("ksz0", ksz0, "ksz0[50]/F");
    tNew->Branch("kspiphi", kspiphi, "kspiphi[50][50]/F");
    tNew->Branch("kspipt", kspipt, "kspipt[50][50]/F");
    tNew->Branch("kspith", kspith, "kspith[50][50]/F");
    tNew->Branch("ksminv", ksminv, "ksminv[50]/F");
    tNew->Branch("kstlen", kstlen, "kstlen[50]/F");
    tNew->Branch("tphi", tphi, "tphi[50]/F");
    tNew->Branch("tth", tth, "tth[50]/F");

    Double_t halfPi = TMath::Pi() / 2;
    double cutChi2r = 15.;
    double cutChi2z = 10.;
    int cutNhitMin = 10;
    int cutNhitMax = 30;
    // double cutRmin = 0.05;
    double cutRmax = 6;
    double cutZtrack = 12.;
    double cutPtot = 40;
    double cutTrackTheta = 0.3;

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

        if(NgoodTr > 1 && is_coll != 1)
        { tNew->Fill(); }
    }

    std::cout << "nentries " << nentries << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;

    top->Write();
    top->Save();
}
