#define PhiToKn_cxx
#include "PhiToKn.hpp"
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <set>
#include <fstream>
#include <algorithm>

void PhiToKn::Loop(std::string output_fname, double energy0)
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

    double mass = 0;
    double z = 0;

    const Double_t halfPi = TMath::Pi() / 2;
    const double cutChi2r = 15.;
    const double cutChi2z = 10.;
    const int cutNhitMin = 6;
    const int cutNhitMax = 30;
    // const double cutRmin = 0.05;
    const double cutRmax = 6;
    const double cutZtrack = 12.;
    const double cutPtot = 40;
    const double cutTrackTheta = 0.3;
    const double max_ks_vertex_z = 8;
    const double min_kstlen = 0.7;
    const double max_kstlen = 6;

    double p1 = 0;
    double p2 = 0;

    double piThetaPos = 0;
    double piThetaNeg = 0;
    double ksCandtlen = 0;

    TFile *top = new TFile(output_fname.c_str(), "recreate");
    auto tNew = new TTree("Kn", "Cutted tr_ph (phi xsec related data)");

    auto hKsMass = new TH1D("hKsMass", "M^{(#pi^{+}#pi^{-})}_{inv}", 640, 420, 580);
    auto hPsi = new TH1D("hPsi", "Psi", 3142, 0, 3.142);
    auto hKsMom = new TH1D("hKsMom", "P_{K_{S}}", 2000, 0, 1000);
    auto hKstlen = new TH1D("hKstlen", "kstlen_{K_{S}}", 100, 0, 10);

    auto hKspipt = new TH1D("hKspipt", "hKspipt", 1000, 0, 1000);

    auto hMissingMass = new TH1D("hMissingMass", "M_{K_{L}}", 1000, 0, 1000);
    auto hKstlenVsMinv = new TH2D("hKstlenVsMinv", "", 100, 0, 10, 640, 420, 580);
    auto hKsMomVsMinv = new TH2D("hKsMomVsMinv", "", 1000, 0, 1000, 640, 420, 580);

    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("trigbits", &trigbits, "trigbits/I");

    tNew->Branch("mass", &mass, "mass/D");
    tNew->Branch("z", &z, "z/D");
    tNew->Branch("lumi", &lumoff, "lumiRun/F");
    tNew->Branch("lumiErr", &lumofferr, "lumiErr/F");
    // tNew->Branch("ksminv", ksminv, "ksminv[6]/F");
    // tNew->Branch("kstlen", kstlen, "kstlen[6]/F");
    tNew->Branch("kstlen", &ksCandtlen, "kstlen/D");
    tNew->Branch("kspipt", kspipt, "kspipt[6][2]/F");

    tNew->Branch("piThetaPos", &piThetaPos, "piThetaPos/D");
    tNew->Branch("piThetaNeg", &piThetaNeg, "piThetaNeg/D");


    bool flag = false; 
    int prevRunnum = -1;
    std::set<int> runNums = {};

    double lumi = 0;
    double lumi_err = 0;

    auto isGoodTrack = [&](int TrackNum) {
        return (tptot[TrackNum] > cutPtot  && 
                fabs(tz[TrackNum]) < cutZtrack && tnhit[TrackNum] > cutNhitMin && 
                tchi2z[TrackNum] < cutChi2z && tchi2r[TrackNum] < cutChi2r);
    };

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        for(int k = 0; k < nks; k++)
        { hKsMom->Fill(ksptot[k]); }
    }

    auto res = hKsMom->Fit("gaus", "SQME", "goff", sqrt(energy0 * energy0 - 497.614 * 497.614) - 15, sqrt(energy0 * energy0 - 497.614 * 497.614) + 15);
    res = hKsMom->Fit("gaus", "SQME", "goff", res->Parameter(1) - 2 * res->Parameter(2), res->Parameter(1) + 2 * res->Parameter(2));
    const double ksMomUpperBound = res->Parameter(1) + 5 * res->Parameter(2);
    const double ksMomLowerBound = energy0 < 513 ? res->Parameter(1) - 5 * res->Parameter(2) : 85.;
    std::cout << "Ks Momentum bounds = " <<  ksMomLowerBound << "--" << ksMomUpperBound << std::endl;
    const double min_psi = 0.8 * 2 * sqrt(TMath::ACos((energy0 * energy0 - 497.614 * 497.614) / (energy0 * energy0 - 4 * 139.57 * 139.57)));
    std::cout << "Min psi_{pions} = " <<  min_psi << std::endl;

    nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        double energy = emeas0;
        std::vector<int> KsCand = {};
        std::vector<double> KlCandMasses = {};
        std::vector<double> KsCandMasses = {};

        for(int k = 0; k < nks; k++)
        {
            if( isGoodTrack(ksvind[k][0]) && isGoodTrack(ksvind[k][1]) &&
                kspipt[k][0] > 130 && kspipt[k][1] > 130 && 
                kspipt[k][0] < 320 && kspipt[k][1] < 320 &&
                tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0 && 
                energy > 100 && 
                tdedx[ksvind[k][0]] < 5000. && 
                tdedx[ksvind[k][1]] < 5000. &&
                
                1.1 < kspith[k][0] && kspith[k][0] < TMath::Pi() - 1.1 && 
                1.1 < kspith[k][1] && kspith[k][1] < TMath::Pi() - 1.1 &&

                kstlen[k] < max_kstlen &&
                fabs(ksz0[k]) < max_ks_vertex_z &&
                ksptot[k] > ksMomLowerBound && ksptot[k] < ksMomUpperBound
            ) 
            {
                flag = true;
                hKspipt->Fill(kspipt[k][0]);
                hKspipt->Fill(kspipt[k][1]);
                posTrackNumber = tcharge[ksvind[k][0]] > 0 ? 0 : 1;
                negTrackNumber = posTrackNumber == 1 ? 0 : 1;

                p1 = kspipt[k][posTrackNumber];
                p2 = kspipt[k][negTrackNumber];

                piPos.SetMagThetaPhi(p1, kspith[k][posTrackNumber], kspiphi[k][posTrackNumber]);
                                        
                piNeg.SetMagThetaPhi(p2, kspith[k][negTrackNumber], kspiphi[k][negTrackNumber]);

                missingMom = -(piPos + piNeg);
                piPosEn = sqrt(139.57 * 139.57 + piPos.Mag2());
                piNegEn = sqrt(139.57 * 139.57 + piNeg.Mag2());

                KsCand.push_back(k);
                hPsi->Fill(piPos.Angle(piNeg));

                if(piPos.Angle(piNeg) < min_psi)
                { flag = false; }

                // if(KsCand.size() > 1)
                // {
                //     flag = false;
                //     break;
                // }
                hKstlen->Fill(kstlen[k]);

                p1 = tptot[ ksvind[k][posTrackNumber] ];
                p2 = tptot[ ksvind[k][negTrackNumber] ];
                piPos.SetMagThetaPhi(p1, tth[ ksvind[k][posTrackNumber] ], tphi[ ksvind[k][posTrackNumber] ]);
                piNeg.SetMagThetaPhi(p2, tth[ ksvind[k][negTrackNumber] ], tphi[ ksvind[k][negTrackNumber] ]);
                
                missingMom = -(piPos + piNeg);
                piPosEn = sqrt(139.57 * 139.57 + piPos.Mag2());
                piNegEn = sqrt(139.57 * 139.57 + piNeg.Mag2());

                TLorentzVector piPlus(piPos, piPosEn);
                TLorentzVector piMinus(piNeg, piNegEn);
                TLorentzVector pseudo(TVector3(0, 0, 0), 2 * energy);
                missingMass = (pseudo - piPlus - piMinus).M();
                KlCandMasses.push_back(missingMass);
                KsCandMasses.push_back((piPlus + piMinus).M());
                hMissingMass->Fill(missingMass);
            }
        }

        if(flag)
        { 
            int candNum = 0;
            {
                auto it = std::min_element(KsCandMasses.begin(), KsCandMasses.end(), [](const double &c1, const double &c2){
                    return fabs(c1 - 497.614) < fabs(c2 - 497.614);
                });
                candNum = std::distance(KsCandMasses.begin(), it);
            }
            
            if(missingMass > 350)
            {
                posTrackNumber = tcharge[ksvind[KsCand[candNum]][0]] > 0 ? 0 : 1;
                negTrackNumber = posTrackNumber == 1 ? 0 : 1;

                piThetaPos = kspith[KsCand[candNum]][posTrackNumber];
                piThetaNeg = kspith[KsCand[candNum]][negTrackNumber];
                ksCandtlen = kstlen[KsCand[candNum]];
                z = ksz0[KsCand[candNum]];

                if(runnum != prevRunnum)
                {
                    lumi += lumoff;
                    lumi_err += lumofferr * lumofferr;
                    prevRunnum = runnum;
                    runNums.insert(runnum);
                }
                
                mass = KsCandMasses[candNum];
                hKstlenVsMinv->Fill(kstlen[KsCand[candNum]], mass);
                hKsMomVsMinv->Fill(ksptot[KsCand[candNum]], mass);
                hKsMass->Fill(ksminv[KsCand[candNum]]);

                tNew->Fill(); 
            }
        }

        KsCandMasses.clear();
        KsCandMasses.shrink_to_fit();

        KlCandMasses.clear();
        KlCandMasses.shrink_to_fit();

        KsCand.clear();
        KsCand.shrink_to_fit();
        flag = false; 
    }
    int n_events = tNew->GetEntries();
    // if(energy0 < 506 || energy0 > 513)
    
    auto hMass = new TH1D("hMass", "Mass without background", 80, 420, 580);
    tNew->Draw("mass >> hMass", "", "goff");
    auto res1 = hMass->Fit("pol0", "SQMEL", "goff", 540, 576);
    // auto res2 = hMass->Fit("pol0", "SQMEL", "goff", 420, 465);
    // double bckgLevel = (res1->Parameter(0) + res2->Parameter(0)) / 2.;
    // double bckgLevelErr = sqrt(res1->ParError(0) * res1->ParError(0) + res2->ParError(0) * res2->ParError(0));
    double bckgLevel = res1->Parameter(0);
    double bckgLevelErr = res1->ParError(0);
    n_events = hMass->Integral() - bckgLevel * hMass->GetNbinsX();
    
    std::cout << "bckgLevel_Left = " << res1->Parameter(0) << "; bckgLevel_Right = " << res1->Parameter(0) << "; bckgLevel_Avg = " << bckgLevel << std::endl;
    std::cout << "N_bckg = " << bckgLevel * hMass->GetNbinsX() << std::endl;
    std::cout << "N_bckg_err = " <<  bckgLevelErr * hMass->GetNbinsX() << std::endl;
    std::cout << "res1 chi2 /ndf = " << res1->Chi2() / res1->Ndf() << std::endl;
    std::cout << "res2 chi2 /ndf = " << res1->Chi2() / res1->Ndf() << std::endl;
    

    std::cout << "luminosity = " << lumi << std::endl;
    std::cout << "luminosity_err = " << sqrt(lumi_err) << std::endl;
    std::cout << "entries = " << n_events << std::endl;
    std::cout << "runCounter = " << runNums.size() << std::endl;

    std::ofstream output_Runs("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Runs/runs_" + std::to_string(energy0) + ".txt");
    for(auto& run : runNums)
    { output_Runs << run << std::endl; }

    top->Write();
}