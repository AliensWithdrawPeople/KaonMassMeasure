#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include <vector>
#include <algorithm>

#include <chrono>
#include <ctime> 

class EnergyHandler
{
private:
    TTree *kTr;
    TTree *ksTr;
    TProfile *pfY;

    Float_t emeas; Float_t demeas;

    Int_t runnum;
    Float_t ksdpsi;

    // Momentum ratio = P1/P2, where P1 is the momentum of pi+, P2 is the momentum of pi-.
    Float_t Y;
    // Energy of beam calculated as sqrt(p^2 + m^2), where p is the momentum of Kch, m is the mass of Kch.
    std::vector<Float_t> e;
    // Error of energy of beam calculated as sqrt(p^2 + m^2), where p is the momentum of Kch, m is the mass of Kch.
    std::vector<Float_t> eErr;

    std::vector<Int_t> entryNum;
    std::vector<Int_t> runEntriesNum;

    std::vector<Float_t> rNum;
    std::vector<Float_t> runnumV;

    std::vector<Float_t> massEMeas;
    std::vector<Float_t> massEMeasErr;

    std::vector<Float_t> eMerge;
    std::vector<Float_t> eMergeErr;
    std::vector<Float_t> groupsStartRunnum;

    int nRunMax;
    std::vector<int> badRuns;
    // Number of groups of runs
    Int_t groupsAmount; 

    // Fills badRuns vector and returns number of good runs.
    int BadRunSearch();
    // Fills vectors e and eErr.
    void EnergyCalculation();
    // Weeds out bad runs.
    void EntryFilter();
    // Merges runs into groups, fills vector massHistMeas and fits each hist.
    size_t RunsDivider();
    double RadCor(TF1 *massFunc);
    void Draw();

public:
    EnergyHandler(std::string fChargedK, std::string fKsKl);
    void MassLnY();
    void MassCriticalAngle();
    ~EnergyHandler();
};

EnergyHandler::EnergyHandler(std::string fChargedK, std::string fKsKl)
{
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree *)file->Get("kChargedTree");
    auto h2EvsRun = new TH2D("dEdXvsPtot", "", 1000, 100, 120, 659, 60741, 61400);
    Int_t runnum_;
    Float_t tdedx[2];
    Float_t tptot[2];
    kTr->SetBranchAddress("runnum", &runnum_);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);

    kTr->Draw("runnum :(tptot[0] + tptot[1]) / 2 >> dEdXvsPtot",
              "(tptot[0] + tptot[1]) / 2 > 100 && (tptot[0] + tptot[1]) / 2 < 200 && (tdedx[0] + tdedx[1]) / 2 > 7e3 && (tdedx[0] + tdedx[1]) / 2 < 25000", "goff");
    pfY = h2EvsRun->ProfileY();

    TFile *file1 = TFile::Open(fKsKl.c_str());
    ksTr = (TTree *)file1->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);

    BadRunSearch();

    delete h2EvsRun;
}

int EnergyHandler::BadRunSearch()
{
    Float_t pAvg = 0;
    int goodRunsCounter = 0;
    for(int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if(emeas == 509.5 && demeas == 0 && std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
        { badRuns.push_back(runnum); }
    }
    for(int i = 0; i < pfY->GetNbinsX(); i++)
    {
        pAvg = pfY->GetBinContent(i);
        if(std::find(badRuns.begin(), badRuns.end(), int(pfY->GetBinCenter(i))) == badRuns.end())
        {
            if(pAvg > 50)
            { goodRunsCounter++; }
            else
            { badRuns.push_back(int(pfY->GetBinCenter(i))); }
        } 
    }
    return goodRunsCounter;
}

void EnergyHandler::EnergyCalculation()
{
    int counter = 0;
    Float_t pAvg = 0;
    Float_t pErr = 0;
    // Invariant mass of charged Kaon
    Float_t minv = 493.677;

    for (int i = 0; i < pfY->GetNbinsX(); i++)
    {
        pAvg = pfY->GetBinContent(i);
        if (std::find(badRuns.begin(), badRuns.end(), int(pfY->GetBinCenter(i))) == badRuns.end())
        {
            pErr = pfY->GetBinError(i);
            rNum.push_back(pfY->GetBinCenter(i) - 0.5);
            e.push_back(TMath::Sqrt(minv * minv + pAvg * pAvg));
            eErr.push_back(pErr * pAvg / *e.rbegin());
            counter++;
        }
    }
}

void EnergyHandler::EntryFilter()
{
    int tmp = 0;
    nRunMax = 0;
    int entriesCounter = 0;
    for (int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if (std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
        {
            if (tmp != runnum)
            {
                entryNum.push_back(i);
                runEntriesNum.push_back(entriesCounter);
                entriesCounter = 0;
                nRunMax++;
            }
            entriesCounter++;
        }
        tmp = runnum;
    }
}

size_t EnergyHandler::RunsDivider()
{
    auto massFunc = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    auto massErrFunc = new TF1("massErr1", "(1 - cos(x / 2) * cos(x / 2)) / TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    auto revMassFunc = new TF1("revMassFunc", "2*TMath::ACos(TMath::Sqrt((1 - x * x / [0] / [0])/(1-4*139.57 * 139.57 / [0] / [0])))");

    auto timeNow = std::time(nullptr);
    auto date = std::localtime(&timeNow);
    std::string rootFileName = "hists and root files/massMeas macro output/massHist" + std::to_string(date->tm_hour) + "h" + 
                                std::to_string(date->tm_min) + "m_"+std::to_string(date->tm_yday)+"day.root";
    TFile file2(rootFileName.c_str(), "recreate");
    std::vector<TH1D *> massHistsMeas;

    auto fermiStep = new TF1("fermiStep", "[2]/(exp(-(x-[0])/[1])+1)/((x-[3])^2+(x-[4])^4 + [5])");
    fermiStep->SetParameters(497.418, 0.335709, 66506, 399.187, 499.504, -9292.33);

    TFitResultPtr res;
    int j = 1; int m = 0; int l = 0; int tmp = 0;
    Float_t eSum = 0; Float_t eErrBackSum = 0;
    Float_t eSumMeas = 0; Float_t eErrBackSumMeas = 0;
    int mhCounter = 0; // counter for massHists
    int mhCounterPrev = -1;
    int sumEntries = 0;
    while (m < nRunMax)
    {
        if (mhCounter != mhCounterPrev)
        { massHistsMeas.push_back(new TH1D(("Meas" + std::to_string(mhCounter)).c_str(), "Ks mass", 150, 494, 510)); }
        j = 0; eSum = 0; eErrBackSum = 0; sumEntries = 0;
        while (m + j < nRunMax && sumEntries < 7000)
        {
            ksTr->GetEntry(entryNum[m + j]);
            sumEntries += runEntriesNum[m + j];
            eSum += e[m + j]; eErrBackSum += 1 / (eErr[m + j] * eErr[m + j]);
            l = entryNum[m + j];
            tmp = runnum;
            while (tmp == runnum && l < entryNum[nRunMax - 1] + runEntriesNum[nRunMax - 1])
            {
                massFunc->SetParameter(0, emeas);
                massHistsMeas[mhCounter]->Fill(massFunc->Eval(ksdpsi));
                massFunc->SetParameter(0, e[m + j]);
                l++;
                ksTr->GetEntry(l);
            }
            j++;
        }
        if (1 / eErrBackSum > 0.00001)
        {
            eMerge.push_back(eSum / j); eMergeErr.push_back(sqrt(1 / eErrBackSum));
            groupsStartRunnum.push_back(rNum[m]);
            for (int i = 0; i < 10; i++)
            {
                fermiStep->SetParameters(497.418, 0.335709, 66506, 399.187, 499.504, -9292.33);
                res = massHistsMeas[mhCounter]->Fit("fermiStep", "LSQE", "", 495.0, 501);
                if (res->IsValid())
                { break; }
            }
            massEMeas.push_back(res->Parameter(0));
            ksTr->GetEntry(entryNum[m + j]);
            massErrFunc->SetParameter(0, emeas);
            revMassFunc->SetParameter(0, emeas);
            massEMeasErr.push_back(res->ParError(0));

            massErrFunc->SetParameter(0, eSum / j);
            revMassFunc->SetParameter(0, eSum / j);

            mhCounter++;
        }
        m += j;
    }
    file2.Write();
    file2.Close();

    for(auto hist : massHistsMeas)
    { delete hist; }
    massHistsMeas.clear();
    
    delete fermiStep; delete massFunc;
    delete massErrFunc; delete revMassFunc; 

    return eMerge.size();
}

void EnergyHandler::Draw()
{
    std::vector<Float_t> emeasGrouped;
    std::vector<Float_t> demeasGrouped;

    for(int i = 0; i < groupsAmount; i++)
    {
        for(int j = 0; j < nRunMax; j++)
        {
            if(abs(rNum[j] - groupsStartRunnum[i]) < 0.5)
            { 
                ksTr->GetEntry(entryNum[j]);
                emeasGrouped.push_back(emeas); 
                demeasGrouped.push_back(demeas);
                break; 
            }
        }
    }

    std::cout<<"Number of groups: " << groupsAmount <<std::endl;

    std::vector<Float_t> zeroes(groupsAmount, 0.0);

    TGraphErrors *gEvsRunMeas = new TGraphErrors(groupsAmount, groupsStartRunnum.data(), emeasGrouped.data(), zeroes.data(), demeasGrouped.data());
    TGraphErrors *gEvsRunCalc = new TGraphErrors(groupsAmount, groupsStartRunnum.data(), eMerge.data(), zeroes.data(), eMergeErr.data());
    TGraphErrors *gMassMeasVsRun = new TGraphErrors(groupsAmount, groupsStartRunnum.data(), massEMeas.data(), zeroes.data(), massEMeasErr.data());
  
    auto c1 = new TCanvas("EnergyStab", "EvsRunMeas, EvsRunCalc, MassVsRunMeas, MassVsRunCalc", 200, 10, 600, 400);
    c1->Divide(2, 2);
    gEvsRunMeas->SetLineColor(2);

    c1->cd(1);
    gEvsRunMeas->SetTitle("E meas vs Run");
    gEvsRunMeas->SetMarkerStyle(kFullDotLarge);
    gEvsRunMeas->DrawClone("AP");

    c1->cd(3);
    gEvsRunCalc->SetTitle("E calc vs Run");
    gEvsRunCalc->SetMarkerStyle(kFullDotLarge);
    gEvsRunCalc->DrawClone("AP");

    c1->cd(2);
    gMassMeasVsRun->SetTitle("Mass (E meas) vs Run");
    gMassMeasVsRun->SetMarkerStyle(kFullDotLarge);
    gMassMeasVsRun->DrawClone("AP");

    delete gEvsRunMeas;
    delete gEvsRunCalc;
    delete gMassMeasVsRun;
}

void EnergyHandler::MassCriticalAngle()
{
    EnergyCalculation();
    EntryFilter();
    
    groupsAmount = RunsDivider();
    Draw();

    double bar = 0;
    double barErr = 0;
    for (int i = 0; i < massEMeas.size() - 1; i++)
    {
        bar += massEMeas[i];
        if (massEMeasErr[i] != 0)
        { barErr += 1 / (massEMeasErr[i] * massEMeasErr[i]); }
        else
        { std::cout<<"massEMeasErr["<< i <<"] = 0" << std::endl; }
    }
    std::cout << "mass(emeas) = " << bar / (massEMeas.size() - 1) << " +- " << 1 / sqrt(barErr) << std::endl;
    
}

void EnergyHandler::MassLnY()
{
    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57 * 139.57 / E^2]))]
    auto massF = new TF1("MassLnY", 
    "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
    
    auto hMlnY = new TH2D("hMlnY", "M(lnY)", 200, -1, 1, 200, 480, 520);
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-6 && Y < 5 && Y > 0.4)
        {
            massF->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            hMlnY->Fill(log(Y), massF->Eval(ksdpsi));
        }
    }
    auto canv = new TCanvas("MlnY","Mass(lnY)", 200, 10, 600, 400);
    hMlnY->DrawClone();

    delete massF;
    delete hMlnY;
}

EnergyHandler::RadCor(TF1 *massFunc, Float_t E, Float_t pRatio, Float_t psi)
{
    double kaonMass = 497.611; double electronMass = 0.511;
    double s = 4 * E * E;
    double beta_ = sqrt(1 - kaonMass * kaonMass / E / E);
    TF1 Li2("Li2", "-TMath::Log(1-t)/t");
    double L = 2 * TMath::Log(2*E / electronMass);
    double a = TMath::Pi()*TMath::Pi() / 6 - 1/4;
    double b = -1 + (1/beta_-1)*L/2 + 1/beta_*TMath::Log((1+beta_)/2) + 
                (1/beta_ + beta_)/2 *(-Li2.Integral(0,(beta_-1)/(1+beta_)) + Li2.Integral(0, (1-beta_)/(1+beta_)) - TMath::Pi()*TMath::Pi()/12 + 
                L*TMath::Log((1+beta_)/2)-2*L*TMath::Log(beta_) + 1.5*(TMath::Log((1+beta_)/2))^2 -0.5*(TMath::Log(beta_))^2 - 
                3*TMath::Log(beta_)*TMath::Log((1+beta_)/2) + L + 2*TMath::Log((1+beta_)/2) );

    /* 
    D(x,s) = D_gamma + D_e+e-. 
    [0] - b=2*alpha/pi * (L-1); [1] - L; [2] - 2*m_e/E.
    For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    */
    TF1 D("D func", "[0]/2*(1-x)^([0]/2-1) * ( 1+3*[0]/8+[0]*[0]/16*(9/8-TMath::Pi()*TMath::Pi()/3)) - [0]/4*(1+x)+ \
        [0]*[0]/32*(4*(1+x)*TMath::Log(1/(1-x)) + (1+3*x^2)/(1-x)*TMath::Log(1/x)-5-x) + \
        [0]/2*(1-x)^([0]/2-1)*(-[0]*[0]/288*(2*[1] - 15)) + (1/137/TMath::Pi())^2 * ( 1/12/(1-x) * \
        (1-x-[2])^([0]/2)*(TMath::Log((1-x)^2 /[2]/[2])-5/3)^2 *(1+x^2+[0]/6*(TMath::Log((1-x)^2 /[2]/[2])-5/3)) +\
        [1]*[1]/4*(2/3*(1-x^3)/x+1/2*(1-x) + (1+x)*TMath::Log(x)) ) * ((1-x-[2] > 0) ? 1 : 0) ");

    /*
    Cross-section of e+e- -> KsKl.
    For definition go to https://arxiv.org/abs/1604.02981v3 or https://doi.org/10.1142/S0217751X92001423.
    TO-DO!!!
    */
    TF1 sigma0("sigma0", "");

    // Krc(s, x1, x2); 
    // For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    TF2 Krc("K RadCor", [&](double* x, double*p) { 
        massFunc->SetParameters(E * (1-x[0]) * (1-x[1]), pRatio);
        double bForDFormula = 2/137/TMath::Pi()*(L - 1);
        D.SetParameters(bForDFormula, L, 2 * electronMass / E);
        sigma0.SetParameters(); // TO-DO!!!
        return (fabs(p[0] - 1) < 0.1 ? massFunc->Eval(psi) : 1) * (1 + 2/137/TMath::Pi() * (1+a+b)) * D.Eval(x[0]) * D.Eval(x[1]) * sigma0(s*(1-x[0]) * (1-x[1])); 
        });
    Krc.SetParameters(0);
    Double_t N = Krc.Integral(0., 1., 0., 1.);
    
    Krc.SetParameters(1);
    return Krc.Integral(0., 1., 0., 1.) / N;
}

EnergyHandler::~EnergyHandler()
{
    delete pfY;
    delete ksTr;
    delete kTr;
}

int massMeasRefactored()
{
    auto start = std::chrono::system_clock::now();
    
    //auto eHandler = new EnergyHandler("hists and root files/cuts/cutKch16Mar22_17h41m.root", "hists and root files/cuts/cutKs28Feb_23h37m.root");
    auto eHandler = new EnergyHandler("hists and root files/cuts/cutKch16Mar22_17h41m.root", "hists and root files/cuts/kskl_2bgen600k(min_nthit == 11 min_rho = 0.1).root");
    
    //eHandler->MassCriticalAngle();
    eHandler->MassLnY();
    delete eHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}