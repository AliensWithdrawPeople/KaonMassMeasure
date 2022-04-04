#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <vector>
#include <algorithm>

class EnergyHandler
{
private:
TCanvas *c1;
TTree *kTr; TTree *ksTr;
TH2D *eMeasVsCalc;

Float_t emeas; Float_t demeas; Int_t runnum; Float_t ksdpsi; Float_t Y;
Int_t *entryNum; Int_t *runEntriesNum;
Float_t *e; Float_t *eErr; Float_t *rNum;
Float_t *runnumV; Float_t *emeasV; Float_t *demeasV;
Float_t *massEMeas; Float_t *massECalc; Float_t *massEMeasErr; Float_t *massECalcErr; 
Float_t *eMerge; Float_t *eMergeErr; Float_t *groupsStartRunnum;

int nRunMax;
std::vector<int> badRuns;
Int_t groupsAmount; // amount of groups of runs 

int BadRunSearch(TProfile *pfY);
void EnergyCalculation(int goodRunsCounter, TProfile *pfY);
void EntryFilter();
size_t RunsDivider();
void DrawShort(TGraphErrors *tg, std::string title, int position = 0);

public:
    EnergyHandler(std::string fChargedK, std::string fKsKl);
    void MassViaY();
    void Draw();
    ~EnergyHandler();
};

EnergyHandler::EnergyHandler(std::string fChargedK, std::string fKsKl)
{
    c1 = new TCanvas("EnergyStab","EvsRunMeas, EvsRunCalc, MassVsRunMeas, MassVsRunCalc", 200, 10, 600, 400);
    
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree*)file->Get("kChargedTree");
    auto h2EvsRun = new TH2D("dEdXvsPtot", "", 1000, 100, 120, 659, 60741, 61400);
    Int_t runnum_; Float_t tdedx[2]; Float_t tptot[2];
    kTr->SetBranchAddress("runnum", &runnum_);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);

    kTr->Draw("runnum :(tptot[0] + tptot[1]) / 2 >> dEdXvsPtot", 
    "(tptot[0] + tptot[1]) / 2 > 100 && (tptot[0] + tptot[1]) / 2 < 200 && (tdedx[0] + tdedx[1]) / 2 > 7e3 && (tdedx[0] + tdedx[1]) / 2 < 25000", "goff");
    auto pfY = h2EvsRun->ProfileY(); 

    TFile *file1 = TFile::Open(fKsKl.c_str());
    ksTr = (TTree*)file1->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas); ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum); ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);
    
    int goodRunsCounter = BadRunSearch(pfY);
    EnergyHandler::EnergyCalculation(goodRunsCounter, pfY);
    EntryFilter();
    groupsAmount = RunsDivider();
}

int EnergyHandler::BadRunSearch(TProfile *pfY)
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

void EnergyHandler::EnergyCalculation(int goodRunsCounter, TProfile *pfY)
{
    e = new Float_t[goodRunsCounter]; eErr = new Float_t[goodRunsCounter];
    rNum = new Float_t[goodRunsCounter];
    
    int counter = 0;
    Float_t pAvg = 0; Float_t pErr = 0; Float_t minv = 493.677;

    for(int i = 0; i < pfY->GetNbinsX(); i++)
    {
        pAvg = pfY->GetBinContent(i);
        if(std::find(badRuns.begin(), badRuns.end(), int(pfY->GetBinCenter(i))) == badRuns.end())
        {
            pErr = pfY->GetBinError(i); rNum[counter] = pfY->GetBinCenter(i) - 0.5;
            e[counter] = TMath::Sqrt(minv * minv + pAvg * pAvg);
            eErr[counter] = pErr * pAvg / e[counter];
            counter++;
        } 
    }
}

void EnergyHandler::EntryFilter()
{
    Long64_t n = ksTr->Draw("emeas : runnum","","goff");
    entryNum = new Int_t[n]; runEntriesNum = new Int_t[n];
    int tmp = 0; nRunMax = 0; int entriesCounter = 0;
    for(int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
        {
            if(tmp != runnum)
            {
                entryNum[nRunMax] = i;
                runEntriesNum[nRunMax] = entriesCounter;
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
    eMeasVsCalc = new TH2D("eMeasVsCalc", "eMeasVsCalc", 10, 509.45, 509.6, 20, 505.7, 505.85);
    auto massFunc = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    auto massErrFunc = new TF1("massErr1", "(1 - cos(x / 2) * cos(x / 2)) / TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    auto revMassFunc = new TF1("revMassFunc", "2*TMath::ACos(TMath::Sqrt((1 - x * x / [0] / [0])/(1-4*139.57 * 139.57 / [0] / [0])))");
    
    std::vector<Float_t> foos; std::vector<Float_t> fooErrs;
    std::vector<Float_t> tres;
    std::vector<Float_t> mass1; std::vector<Float_t> mass2;
    std::vector<Float_t> mass1Err; std::vector<Float_t> mass2Err;

    int j = 1; int m = 0; int l = 0; int tmp = 0;
    Float_t eSum = 0; Float_t eErrBackSum = 0;

    TFile file2("hists and root files/massHist.root", "recreate");
    std::vector<TH1D*> massHistsCalc; std::vector<TH1D*> massHistsMeas;
    auto fermiStep = new TF1("fermiStep", "[2] / (exp(-(x - [0])/[1]) + 1) / (x-[3]) / (x - [3])");
    fermiStep->SetParameters(497.253, 0.35052, 8030.5, 490.363);
    TFitResultPtr res;
    int mhCounter = 0; // counter for massHists
    int mhCounterPrev = -1;
    int sumEntries = 0;
    while(m < nRunMax)
    {
        if(mhCounter != mhCounterPrev)
        {
            massHistsMeas.push_back(new TH1D(("Meas" + std::to_string(mhCounter)).c_str(), "Ks mass", 150, 494, 510)); 
            massHistsCalc.push_back(new TH1D(("Calc" + std::to_string(mhCounter)).c_str(), "Ks mass", 150, 490, 505));
        }
        j = 0; eSum = 0;
        eErrBackSum = 0; sumEntries = 0;
        while(m + j < nRunMax && sumEntries < 7000)
        {  
            sumEntries += runEntriesNum[m + j];
            eSum += e[m + j]; 
            eErrBackSum += 1/ (eErr[m + j] * eErr[m + j]);
            ksTr->GetEntry(entryNum[m + j]);
            l = entryNum[m + j];
            tmp = runnum;
            while(tmp == runnum && l < entryNum[nRunMax - 1] + runEntriesNum[nRunMax - 1])
            {
                massFunc->SetParameter(0, emeas);
                massHistsMeas[mhCounter]->Fill(massFunc->Eval(ksdpsi));
                massFunc->SetParameter(0, e[m + j]);
                massHistsCalc[mhCounter]->Fill(massFunc->Eval(ksdpsi));
                l++; ksTr->GetEntry(l);
            }
            j++;
        }
        if(1 / eErrBackSum > 0.00001)
        { 
            foos.push_back(eSum / j); fooErrs.push_back(sqrt(1 / eErrBackSum)); tres.push_back(rNum[m]);
            for(int i = 0; i < 10; i++)
            {
                fermiStep->SetParameters(497.253, 0.35052, 8030.5, 490.363);
                res = massHistsMeas[mhCounter]->Fit("fermiStep", "LSQ", "", 495.0, 502);
                if(res->IsValid())
                { break; }
            }
            mass1.push_back(res->Parameter(0));
            massErrFunc->SetParameter(0, emeas); revMassFunc->SetParameter(0, emeas);

            //std::cout << mhCounter << " : mass1 = " << res->Parameter(0) << " : mass1Err = " <<demeas * massErrFunc->Eval(revMassFunc->Eval(res->Parameter(0))) << std::endl;
            mass1Err.push_back(demeas * massErrFunc->Eval(revMassFunc->Eval(res->Parameter(0))));

            for(int i = 0; i < 10; i++)
            {
                fermiStep->SetParameters(493.731, 0.35017, 6833.16, 489.16);
                res = massHistsCalc[mhCounter]->Fit("fermiStep", "LSQ", "", 491.0, 498.0);
                if(res->IsValid())
                { break; }
            }
            mass2.push_back(res->Parameter(0));
            massErrFunc->SetParameter(0, eSum / j); revMassFunc->SetParameter(0, eSum / j);
            mass2Err.push_back(sqrt(1 / eErrBackSum) * massErrFunc->Eval(revMassFunc->Eval(res->Parameter(0))));

            mhCounter++;
        }
        m += j;
    }
    file2.Write();
    file2.Close();

    massEMeas = mass1.data(); massECalc = mass2.data();
    massEMeasErr = mass1Err.data(); massECalcErr = mass2Err.data();
    eMerge = foos.data(); eMergeErr = fooErrs.data();
    groupsStartRunnum = tres.data();

    emeasV = new Float_t[nRunMax]; demeasV = new Float_t[nRunMax];

    for(int i = 0; i < nRunMax; i++)
    {
        ksTr->GetEntry(entryNum[i]);
        emeasV[i] = emeas; demeasV[i] = demeas;
        eMeasVsCalc->Fill(emeas, e[i]);
    }

    double bar = 0; double barErr = 0;
    for(int i = 0; i < mass1.size() - 1; i++)
    {
        bar += mass1[i];
        if(mass1Err[i] != 0)
        { barErr += 1 / (mass1Err[i] * mass1Err[i]); }
    }
    std::cout << "mass(emeas) = " << bar / (mass1.size() - 1) << " +- " << 1 / sqrt(barErr) << std::endl;
    return foos.size();
}

void EnergyHandler::Draw()
{
    Float_t *zeroesVect1 = new Float_t[nRunMax]; Float_t *zeroesVect2 = new Float_t[groupsAmount];
    for(int i = 0; i < nRunMax; i++)
    { zeroesVect1[i] = 0; }
    for(int i = 0; i < groupsAmount; i++) 
    { zeroesVect2[i] = 0; }

    Float_t *eDiff = new Float_t[groupsAmount]; Float_t *eDiffErr = new Float_t[groupsAmount];
    Float_t *emeasVgrouped = new Float_t[groupsAmount]; Float_t *demeasVgrouped = new Float_t[groupsAmount];
    for(int i = 0; i < groupsAmount; i++)
    {
        for(int j = 0; j < nRunMax; j++)
        {
            if(abs(rNum[j] - groupsStartRunnum[i]) < 0.5)
            { 
                emeasVgrouped[i] = emeasV[j]; demeasVgrouped[i] = demeasV[j];
                eDiff[i] = emeasV[j] - eMerge[i];
                eDiffErr[i] = sqrt(demeasV[j] * demeasV[j] + eMergeErr[i] * eMergeErr[i]); 
                break; 
            }
        }
    }

    std::cout<<"grpsAmount: " << groupsAmount <<std::endl;

    //TGraphErrors *gEvsRunMeas = new TGraphErrors(nRunMax, rNum, emeasV, zeroesVect1, demeasV);
    TGraphErrors *gEvsRunMeas = new TGraphErrors(groupsAmount, groupsStartRunnum, emeasVgrouped, zeroesVect2, demeasVgrouped);
    TGraphErrors *gEvsRunCalc = new TGraphErrors(groupsAmount, groupsStartRunnum, eMerge, zeroesVect2, eMergeErr);
    TGraphErrors *gEDiffvsRun = new TGraphErrors(groupsAmount, groupsStartRunnum, eDiff, zeroesVect2, eDiffErr);
    TGraphErrors *gMassMeasVsRun = new TGraphErrors(groupsAmount, groupsStartRunnum, massEMeas, zeroesVect2, massEMeasErr);
    TGraphErrors *gMassCalcVsRun = new TGraphErrors(groupsAmount, groupsStartRunnum, massECalc, zeroesVect2, massECalcErr);
  
    c1->Divide(2, 3);
    gEvsRunMeas->SetLineColor(2);

    c1->cd(1);
    gEvsRunMeas->SetTitle("E meas vs Run");
    gEvsRunMeas->SetMarkerStyle(kFullDotLarge);
    gEvsRunMeas->DrawClone("AP");

    c1->cd(2);
    gMassCalcVsRun->SetTitle("Mass (E calc) vs Run");
    gMassCalcVsRun->SetMarkerStyle(kFullDotLarge);
    gMassCalcVsRun->DrawClone("AP");

    c1->cd(3);
    gEvsRunCalc->SetTitle("E calc vs Run");
    gEvsRunCalc->SetMarkerStyle(kFullDotLarge);
    gEvsRunCalc->DrawClone("AP");

    c1->cd(5);
    gEDiffvsRun->SetTitle("E diff vs Run");
    gEDiffvsRun->SetMarkerStyle(kFullDotLarge);
    gEDiffvsRun->DrawClone("AP");

    c1->cd(4);
    gMassMeasVsRun->SetTitle("Mass (E meas) vs Run");
    gMassMeasVsRun->SetMarkerStyle(kFullDotLarge);
    gMassMeasVsRun->DrawClone("AP");

    delete[] zeroesVect1; delete[] zeroesVect2;
    delete[] eDiff; delete[] eDiffErr; delete[] emeasVgrouped; delete[] demeasVgrouped;
}

EnergyHandler::~EnergyHandler()
{
    delete[] entryNum; delete[] runEntriesNum; delete[] e; delete[] eErr; 
    delete[] rNum; delete[] emeasV; delete[] demeasV; 
    delete[] massEMeas; delete[] massECalc; delete[] massEMeasErr; delete[] massECalcErr;
    delete[] eMerge; delete[] eMergeErr; delete[] groupsStartRunnum;
}

void EnergyHandler::MassViaY()
{
    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57 * 139.57 / E^2]))]
    auto massF = new TF1("massViaY", 
    "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
    std::vector<Float_t> mY; std::vector<Float_t> tmp;

    auto hMlnY = new TH2D("hMlnY", "M(lnY)", 200, -1, 1, 200, 480, 520);
    int countIt = 0;
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);

        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-6 && Y < 5 && Y > 0.4)
        {
            countIt++;
            massF->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            mY.push_back(massF->Eval(ksdpsi));
            tmp.push_back(log(Y));
            hMlnY->Fill(log(Y), massF->Eval(ksdpsi));
        }
    }
    Float_t *massY = mY.data(); Float_t *lnY = tmp.data();
    TGraph *gMassLnY = new TGraphErrors(groupsAmount, lnY, massY);
    std::cout << countIt << std::endl;
    std::cout << mY.size() << std::endl;
    auto canv = new TCanvas("MlnY","Mass(lnY)", 200, 10, 600, 400);
    gMassLnY->SetTitle("Mass(lnY)"); gMassLnY->SetMarkerStyle(kFullDotLarge);
    //gMassLnY->Draw("AP");

    hMlnY->DrawClone();
    //delete [] massY; delete [] lnY;
}


int massMeas()
{
    auto eHandler = new EnergyHandler("hists and root files/cuts/cutKch16Mar22_17h41m.root", "hists and root files/cuts/cutKs28Feb_23h37m.root");
    eHandler->Draw();
    eHandler->MassViaY();
    //delete eHandler;
    return 0;
}