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
#include <tuple>
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
    double sigmaPsi;
    double avgPsiFullRec;
    
    std::vector<int> badRuns;

    Float_t emeas; Float_t demeas;
    Int_t runnum; Float_t ksdpsi;
    Float_t Y;

    void CalcResolutions();

public:
    Corrections(std::string fKsKl, std::vector<int> badRuns);
    double CalcCor(double massKs, double energy, double psiAvg, double sigmaPsi, bool isCriticalAngleMethod = false);
    double GetCorrectedMass(double massKs, double energy, bool isCriticalAngleMethod = true);
    std::tuple<double, double> GetResolution(bool isCriticalAngleMethod = true);
    ~Corrections();
};

Corrections::Corrections(std::string fKsKl, std::vector<int> bRuns)
{
    TFile *file = TFile::Open(fKsKl.c_str());
    ksTr = (TTree *)file->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);
    
    for(auto item : bRuns)
    { badRuns.push_back(item); }
    massFuncCrAngle = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    massFuncFullRec = new TF1("MassLnY", 
    "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
    CalcResolutions();
}

Corrections::~Corrections()
{
    delete ksTr;
    delete massFuncCrAngle;
    delete massFuncFullRec;
}

void Corrections::CalcResolutions()
{
    auto hPsilnY = new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 200, 2.45, 2.8);
    auto hPsi = new TH1D("hPsi", "Psi", 200, 2.4, 3.2);

    double Yavg = 0;
    double eAvg = 0;
    int counter = 0;
    for(int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
        {
            eAvg += emeas;
            Yavg += Y;
            counter++;
            hPsi->Fill(ksdpsi);
            if(fabs(Y - 1) > 1e-6 && fabs(log(Y)) <= 0.2)
            { hPsilnY->Fill(log(Y), ksdpsi); }
        }
    }
    Yavg /=counter;
    eAvg /= counter;
    auto fermiStep = new TF1("fermiStep", "[2]/(exp(-(x-[0])/[1])+1)/((x-[3])^4 + (x-[4])^2)");
    TFitResultPtr res;
    for(int i = 0; i < 10; i++)
    {
        fermiStep->SetParameters(2.61082, 0.0102736, 4.36229, 2.29701);
        res = hPsi->Fit("fermiStep", "SQE0", "", 2.50, 2.65);
        if (res->IsValid())
        { break; }
    }
    avgPsiCrAngle = res->Parameter(0);
    sigmaPsiCrAngle = res->ParError(0);

    auto psiProjY = hPsilnY->ProjectionY();
    for(int i = 0; i < 10; i++)
    {
        res = psiProjY->Fit("gaus", "LSQE", "", 2.58, 2.65);
        if (res->IsValid())
        { break; }
    }
    avgPsiFullRec = res->Parameter(1); //2.61545; // DO IT CORRECTLY!!!
    sigmaPsi = res->Parameter(2); //0.0226; // DO IT CORRECTLY!!!

    massFuncCrAngle->SetParameter(0, eAvg);
    massFuncFullRec->SetParameters(eAvg, (1 - Yavg * Yavg) / (1 + Yavg * Yavg));
    
    delete fermiStep;
    delete hPsi;
    delete hPsilnY;
}

double Corrections::CalcCor(double massKs, double energy, double psiAvg, double sigmaPsi, bool isCriticalAngleMethod = false)
{
    double Mcorr = massKs - sigmaPsi*sigmaPsi / 2 * (isCriticalAngleMethod? massFuncCrAngle->Derivative2(psiAvg) : massFuncFullRec->Derivative2(psiAvg));
    return Mcorr;
}

std::tuple<double, double> Corrections::GetResolution(bool isCriticalAngleMethod = true)
{
    if(isCriticalAngleMethod)
    { return std::make_tuple(avgPsiCrAngle, sigmaPsiCrAngle); }
    else
    { return std::make_tuple(avgPsiFullRec, sigmaPsi); }
}

double Corrections::GetCorrectedMass(double massKs, double energy, bool isCriticalAngleMethod = true)
{
    double massCorrected = -1;
    if(isCriticalAngleMethod)
    { massCorrected = CalcCor(massKs, energy, avgPsiCrAngle, sigmaPsi, true); }
    else
    { massCorrected = CalcCor(massKs, energy, avgPsiFullRec, sigmaPsi, false); }
    return massCorrected;
}

void corrections()
{
    auto cor = new Corrections("hists and root files/cuts/kskl_2bgen600k(min_nthit == 11 min_rho = 0.1).root", {});
    std::cout << cor->GetCorrectedMass(497.579, 510, true) << std::endl;
    // FullRec 497.604
    // CrAngle 497.579
    delete cor;
}