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
    double avgPsi;
    double sigmaPsi;

    std::vector<int> badRuns;

    Float_t emeas; Float_t demeas;
    Int_t runnum; Float_t ksdpsi;
    Float_t Y;

    void CalcResolutions();

public:
    Corrections(std::string fKsKl, std::vector<int> badRuns);
    double CalcCor(double massKs, double energy, double psiAvg, double sigmaPsi, bool isCriticalAngleMethod = false);
    double GetCorrectedMass(double massKs, double energy, bool isCriticalAngleMethod = true);
    std::tuple<double, double> GetResolution();
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
            if(fabs(Y - 1) > 1e-9)
            { hPsilnY->Fill(log(Y), ksdpsi); }
        }
    }
    Yavg /=counter;
    eAvg /= counter;
    TFitResultPtr res;
    for(int i = 0; i < 10; i++)
    {
        res = hPsilnY->ProfileX()->Fit("pol6", "SQE0", "", -0.4, 0.4);
        if (res->IsValid())
        { break; }
    }
    avgPsi = res->Parameter(0);
    sigmaPsi = 1.64407e-02; //res->ParError(0);
    std::cout << "sigma_psi = " << sigmaPsi << std::endl;

    massFuncCrAngle->SetParameter(0, eAvg);
    massFuncFullRec->SetParameters(eAvg, (1 - Yavg * Yavg) / (1 + Yavg * Yavg));
    
    delete hPsi;
    delete hPsilnY;
}

double Corrections::CalcCor(double massKs, double energy, double psiAvg, double sigmaPsi, bool isCriticalAngleMethod = false)
{
    double Mcorr = sigmaPsi * sigmaPsi / 2 * (isCriticalAngleMethod? massFuncCrAngle->Derivative2(psiAvg) : massFuncFullRec->Derivative2(psiAvg));
    double sigmaOfSigmaPsi = 1.04220e-04;
    std::cout << "sigma_deltaM = " << sigmaOfSigmaPsi * 2 * abs(Mcorr) * 1000 << " keV" << std::endl;
    return Mcorr;
}

std::tuple<double, double> Corrections::GetResolution()
{ return std::make_tuple(avgPsi, sigmaPsi); }

double Corrections::GetCorrectedMass(double massKs, double energy, bool isCriticalAngleMethod = true)
{
    double massCorrected = -1;
    massCorrected = massKs - CalcCor(massKs, energy, avgPsi, sigmaPsi, isCriticalAngleMethod);
    return massCorrected;
}

void corrections()
{
    auto cor = new Corrections("hists and root files/cuts/kskl_2bgen600k(min_nthit == 11 min_rho = 0.1).root", {});
    std::cout << cor->GetCorrectedMass(497.601, 510, false) << " MeV" << std::endl;
    // FullRec 2body 497.602 +- 0.003 MeV
    // CrAngle 2body 497.623 +- 0.007 MeV
    // FullRec MCGPJ 497.724 +- 0.003 MeV
    // CrAngle MCGPJ 497.726 +- 0.006 MeV
    delete cor;
}