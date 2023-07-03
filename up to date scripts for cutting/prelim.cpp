#define prelim_cxx
#include "prelim.h"
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


void prelim::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      root> .L prelim.C
//      root> prelim t
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
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    TFile *top = new TFile(histFileName.c_str(), "recreate");
    auto tNew = new TTree("ksPrelim", "Preliminary cut tree");

    tNew->Branch("emeas0", &emeas0, "emeas0/F");
    tNew->Branch("demeas0", &demeas0, "demeas0/F");
    tNew->Branch("emeas", &emeas, "emeas/F");
    tNew->Branch("demeas", &demeas, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/F");
    tNew->Branch("nks", &nks, "nks/F");

    tNew->Branch("tcharge", tcharge, "tcharge[15]/I");
    tNew->Branch("ksvind", ksvind, "ksvind[15][2]/I");
    tNew->Branch("kspiphi", kspiphi, "kspiphi[15][2]/F");
    tNew->Branch("kspith", kspith, "kspith[15][2]/F");
    tNew->Branch("kspipt", kspipt, "kspipt[15][2]/F");
    
    tNew->Branch("ksminv", ksminv, "ksminv[15]/F");
    tNew->Branch("kstlen", kstlen, "kstlen[15]/F");
    tNew->Branch("ksalign", ksalign, "ksalign[15]/F");

    Double_t halfPi = TMath::Pi() / 2;
    double cutChi2r = 15.;
    double cutChi2z = 10.;
    int cutNhitMin = 6;
    // double cutRmin = 0.05;
    double cutZtrack = 12.;
    double cutPtot = 40;
    double cutTrackTheta = 0.3;

    auto isGoodTrack = [&](int TrackNum) {
        return (tptot[TrackNum] > cutPtot  && 
                fabs(tz[TrackNum]) < cutZtrack && tnhit[TrackNum] > cutNhitMin && 
                tchi2z[TrackNum] < cutChi2z && tchi2r[TrackNum] < cutChi2r);
    };

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        for(int k = 0; k < nks; k++)
        {
            if( isGoodTrack(ksvind[k][0]) && isGoodTrack(ksvind[k][1]) &&
                (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000 &&
                1 < kspith[k][0] && kspith[k][0] < TMath::Pi() - 1 && 
                1 < kspith[k][1] && kspith[k][1] < TMath::Pi() - 1 && 
                kspipt[k][0] > 130 && kspipt[k][1] > 130 && 
                kspipt[k][0] < 320 && kspipt[k][1] < 320 &&
                tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0)
            { tNew->Fill(); }
        }
    }

    std::cout << "nentries = " << tNew->GetEntries() << std::endl;
    std::cout << "e_MC =  " << double(tNew->GetEntries()) / nentries << std::endl;

    top->Write();
    top->Save();
}
