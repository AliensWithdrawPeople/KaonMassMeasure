#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVector3.h"


int CalcThetaSigmas()
{
    auto filename = [](std::string energyPoint) { return "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/pipi/pipiExp" + energyPoint + ".root"; };
    std::vector<std::string> points = {"210", "230", "250", "270", "274", "280", "290", "295"};

    auto tree = TFile::Open(filename("274").c_str())->Get<TTree>("pion");
    double phiPos;
    double thetaPos;
    double momPos;
    double phiNeg;
    double thetaNeg;
    double momNeg;

    double phiPos_v;
    double thetaPos_v;
    double momPos_v;
    double phiNeg_v;
    double thetaNeg_v;
    double momNeg_v;

    tree->SetBranchAddress("phiPos", &phiPos);
    tree->SetBranchAddress("thetaPos", &thetaPos);
    tree->SetBranchAddress("momPos", &momPos);
    tree->SetBranchAddress("phiNeg", &phiNeg);
    tree->SetBranchAddress("thetaNeg", &thetaNeg);
    tree->SetBranchAddress("momNeg", &momNeg);

    tree->SetBranchAddress("phiPos_v", &phiPos_v);
    tree->SetBranchAddress("thetaPos_v", &thetaPos_v);
    tree->SetBranchAddress("momPos_v", &momPos_v);
    tree->SetBranchAddress("phiNeg_v", &phiNeg_v);
    tree->SetBranchAddress("thetaNeg_v", &thetaNeg_v);
    tree->SetBranchAddress("momNeg_v", &momNeg_v);

    auto hDeltaThetaVsThetaPos = new TH2D("hDeltaThetaVsThetaPos", "hDeltaThetaVsThetaPos", 600, 0, 3.15, 20000, -1, 1);
    auto hAngleVsThetaPos = new TH2D("hAngleVsThetaPos", "hAngleVsThetaPos", 600, 0, 3.15, 20000, -1, 1);
    TVector3 pos(1, 0, 0);
    TVector3 neg(1, 0, 0);

    for(int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        pos.SetMagThetaPhi(1, thetaPos, phiPos);
        neg.SetMagThetaPhi(1, thetaNeg, phiNeg);

        hDeltaThetaVsThetaPos->Fill(thetaPos, thetaPos + thetaNeg - TMath::Pi());
        hAngleVsThetaPos->Fill(thetaPos, pos.Angle(neg) - TMath::Pi());
    }

    hDeltaThetaVsThetaPos->DrawClone("col");
    hAngleVsThetaPos->DrawClone("col");

    return 0;
}