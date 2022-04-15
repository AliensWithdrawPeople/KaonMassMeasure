#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include <vector>
#include <algorithm>


class Corrections
{
private:
    TTree *ksTr;

    // [0] = E
    TF1 *massFuncCrAngle;
    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57 * 139.57 / E^2]))]
    TF1 *massFuncFullRec;

    // sigma_psi - uncertainty of critical angle between pi+ and pi-(particles after decay of Ks).
    double sigmaPsiCrAngle;
    double avgPsiCrAngle;
    double sigmaPsiFullRec;
    double avgPsiFullRec;
    
    double energyPoint;
    std::vector<int> badRuns;

    Float_t emeas; Float_t demeas;
    Int_t runnum; Float_t ksdpsi;
    Float_t Y;

    void GetResolutions();
    double AvgCor(double massKs, double energy, double psiAvg, double sigmaPsi, TF1 *massFunc, bool isCriticalAngleMethod = false);

public:
    Corrections(std::string fKsKl, std::vector<int> badRuns, double energy);
    ~Corrections();
};

Corrections::Corrections(std::string fKsKl, std::vector<int> bRuns, double energy)
{
    TFile *file = TFile::Open(fKsKl.c_str());
    ksTr = (TTree *)file->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);

    energyPoint = energy;
    for(auto item : bRuns)
    { badRuns.push_back(item); }
    massFuncCrAngle = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    massFuncFullRec = new TF1("MassLnY", 
    "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
}

Corrections::~Corrections()
{
    delete ksTr;
    delete massFuncCrAngle;
    delete massFuncFullRec;
}

void Corrections::GetResolutions()
{
    auto hPsilnY = new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 200, 2, 4);
    auto hPsi = new TH1D("hPsi", "Psi", 200, 2.4, 3.2);

    double Yavg = 0;
    double eAvg = 0;
    int counter = 0;
    for(int i = 0; ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && fabs(Y - 1) > 1e-6 && fabs(log(Y)) < 0.4)
        {
            Yavg += Y;
            eAvg += emeas;
            counter++;
            hPsi->Fill(ksdpsi);
            hPsilnY->Fill(log(Y), ksdpsi);
        }
    }
    auto fermiStep = new TF1("fermiStep", "[2]/(exp(-(x-[0])/[1])+1)/((x-[3])^2 + [5])");
    TFitResultPtr res;
    for(int i = 0; i < 10; i++)
    {
        fermiStep->SetParameters(/* TO-DO!!! */);
        res = hPsi->Fit("fermiStep", "LSQE0", "", 2.54, 2.66);
        if (res->IsValid())
        { break; }
    }
    avgPsiCrAngle = res->Parameter(0);
    sigmaPsiCrAngle = res->ParError(0);
    Yavg /= counter;
    eAvg /= counter;

    massFuncCrAngle->SetParameter(0, eAvg);
    massFuncFullRec->SetParameter(eAvg, (1 - Yavg * Yavg) / (1 + Yavg * Yavg));
}


double Corrections::AvgCor(double massKs, double energy, double psiAvg, double sigmaPsi, TF1 *massFunc, bool isCriticalAngleMethod = false)
{
    double gamma2 = energy * energy / (massKs * massKs);
    // R = 4*M^2(pi+-)/M^2(Ks).
    double R = 4 * 139.57 * 139.57 / (massKs * massKs);
    double Mcorr = -1;
    if(isCriticalAngleMethod)
    { Mcorr = massKs * (1 - 1./8 * sigmaPsi * sigmaPsi * (R * gamma2 - 1)); }
    else
    {
        TF1 avgCorrected("avg non-linearity", 
        [&](double *x, double *p) 
        { return massFunc->Eval(x[0]) * 1./sqrt(2 * TMath::Pi() * p[0]) * exp(-(x[0] - p[1]) * (x[0] - p[1]) / (2 * p[0]) ); }, 0, 10000, 1);
        avgCorrected.SetParameters(sigmaPsi * sigmaPsi, psiAvg);
        Mcorr = avgCorrected.Integral(0, TMath::Pi());
    }
    return Mcorr;
}