#define PhiToKn_cxx
#include "PhiToKn.h"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TVector3.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

void PhiToKn::Loop(std::string output_fname, double energy)
{
//   In a ROOT session, you can do:
//      root> .L PhiToKn.C
//      root> PhiToKn t
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

    TVector3 piPos;
    TVector3 piNeg;
    TVector3 missingMom;

    double piPosEn = 0;
    double piNegEn = 0;
    double missingMass = 0;
    int posTrackNumber = 0;
    int negTrackNumber = 0;

    Double_t halfPi = TMath::Pi() / 2;
    double cutChi2r = 15.;
    double cutChi2z = 10.;
    int cutNhitMin = 6;
    int cutNhitMax = 30;
    // double cutRmin = 0.05;
    double cutRmax = 6;
    double cutZtrack = 12.;
    double cutPtot = 40;
    double cutTrackTheta = 0.3;

    double p1 = 0;
    double p2 = 0;

    std::vector<double> KsCandMasses;

    TFile *top = new TFile(output_fname.c_str(), "recreate");
    auto tNew = new TTree("Kn", "Cutted tr_ph (phi xsec related data)");

    auto hKsMass = new TH1D("hKsMass", "M^{(#pi^{+}#pi^{-})}_{inv}", 1000, 450, 550);
    auto hKsMom = new TH1D("hKsMom", "P_{K_{S}}", 2000, 0, 1000);
    auto hKstlen = new TH1D("hKstlen", "P_{K_{S}}", 2000, 0, 10);

    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksminv", &ksminv, "ksminv/F");
    tNew->Branch("kstlen", &kstlen, "kstlen/F");
    bool flag = false; 

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        for(int k = 0; k < nks; k++)
        { hKsMom->Fill(ksptot[k]); }
    }

    auto res = hKsMom->Fit("gaus", "SQME", "goff", sqrt(energy * energy - 497.614 * 497.614) - 15, sqrt(energy * energy - 497.614 * 497.614) + 15);
    res = hKsMom->Fit("gaus", "SQME", "goff", res->Parameter(1) - 2 * res->Parameter(2), res->Parameter(1) + 2 * res->Parameter(2));
    const double ksMomUpperBound = res->Parameter(1) + 5 * res->Parameter(2);
    const double ksMomLowerBound = res->Parameter(1) - 5 * res->Parameter(2);

    nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        std::vector<double> KsCandMasses = {};
        for(int k = 0; k < nks; k++)
        {
            if( (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000. &&
                1 < kspith[k][0] && kspith[k][0] < TMath::Pi() - 1 && 
                1 < kspith[k][1] && kspith[k][1] < TMath::Pi() - 1 &&
                ksalign[k] > 0.85 && kstlen[k] < 6. &&
                ksptot[k] > ksMomLowerBound && ksptot[k] < ksMomUpperBound
                ) 
            {
                posTrackNumber = tcharge[ksvind[k][0]] > 0 ? 0 : 1;
                negTrackNumber = posTrackNumber == 1 ? 0 : 1;

                p1 = kspipt[k][posTrackNumber];
                p2 = kspipt[k][negTrackNumber];

                piPos.SetMagThetaPhi(p1, kspith[k][posTrackNumber], kspiphi[k][posTrackNumber]);
                                        
                piNeg.SetMagThetaPhi(p2, kspith[k][negTrackNumber], kspiphi[k][negTrackNumber]);

                missingMom = -(piPos + piNeg);
                piPosEn = sqrt(139.57 * 139.57 + piPos.Mag2());
                piNegEn = sqrt(139.57 * 139.57 + piNeg.Mag2());
                missingMass = sqrt(4 * emeas0 * emeas0 + 2 * 139.57 * 139.57 
                                        - 2 * 2 * emeas0 * piPosEn - 2 * 2 * emeas0 * piNegEn
                                        + 2 * (piPosEn * piNegEn - piPos.Dot(piNeg)));

                KsCandMasses.push_back(sqrt(2 * 139.57 * 139.57 + 2 * (piPosEn * piNegEn - piPos.Dot(piNeg)) ) );
                flag = true;
                hKstlen->Fill(kstlen[k]);
            }
        }
        if(flag)
        { 
            double tmp = fabs(KsCandMasses[0] - 497.614);
            int candNum = 0;
            for(int i = 1; i < KsCandMasses.size(); i++) 
            { 
                if(fabs(KsCandMasses[i] - 497.614) < tmp)
                { 
                    tmp = fabs(KsCandMasses[i] - 497.614); 
                    candNum = i;
                }
            }

            hKsMass->Fill(KsCandMasses[candNum]);
            tNew->Fill(); 
        }

        KsCandMasses.clear();
        KsCandMasses.shrink_to_fit();
        flag = false; 
    }
    std::cout << "e_MC = " << double(tNew->GetEntries()) / nentries << std::endl;
    top->Write();
    // top->Save();
}