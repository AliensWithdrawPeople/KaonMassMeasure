#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
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

    TF1 *massCrAngle;
    TF1 *massFullRec;

    Float_t emeas; Float_t demeas;
    Int_t runnum; Float_t ksdpsi;

    // Momentum ratio = P1/P2, where P1 is the momentum of pi+, P2 is the momentum of pi-.
    Float_t Y;
    // Energy of beam calculated as sqrt(p^2 + m^2), where p is the momentum of Kch, m is the mass of Kch.
    std::vector<Float_t> e;
    // Error of energy of beam calculated as sqrt(p^2 + m^2), where p is the momentum of Kch, m is the mass of Kch.
    std::vector<Float_t> eErr;

    // First entry of run; Last entry start with (ksTr->GetEntriesFast() - 1), because it's just a mark of the end.
    std::vector<Int_t> entryNum;
    // Number of entries in a run; Last run has 0 entries.
    std::vector<Int_t> runEntriesNum;

    std::vector<Float_t> rNum;
    std::vector<Float_t> runnumV;

    std::vector<Float_t> massEMeas;
    std::vector<Float_t> massEMeasErr;

    std::vector<Float_t> eMerge;
    std::vector<Float_t> eMergeErr;
    std::vector<Float_t> groupsStartRunnum;


    std::vector<Float_t> vSigma;

    int nRunMax;
    double sigmaPsiCrAngle;
    double avgPsiCrAngle;

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
    void Draw();

public:
    EnergyHandler(std::string fChargedK, std::string fKsKl, std::vector<Float_t> sigmas);
    void MassLnY(int drawOpt = 0);
    double GetMassCriticalAngle();

    std::vector<int> GetBadRunsList();

    ~EnergyHandler();
};

EnergyHandler::EnergyHandler(std::string fChargedK, std::string fKsKl, std::vector<Float_t> sigmas)
{
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree *)file->Get("kChargedTree");
    TH2D h2EvsRun("dEdXvsPtot", "", 1000, 100, 120, 659, 60741, 61400);
    Int_t runnum_;
    Float_t tdedx[2];
    Float_t tptot[2];

    vSigma = sigmas;

    // [0] = energy;
    massCrAngle = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57018 * 139.57018 / [0] / [0]) * cos(x / 2) * cos(x / 2))");

    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57 * 139.57 / E^2]))]
    massFullRec = new TF1("MassLnY", 
    "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57018 * 139.57018 / [0] / [0]) ) )/ [1] / [1] ) )");

    kTr->SetBranchAddress("runnum", &runnum_);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);

    kTr->Draw("runnum :(tptot[0] + tptot[1]) / 2 >> dEdXvsPtot",
              "(tptot[0] + tptot[1]) / 2 > 100 && (tptot[0] + tptot[1]) / 2 < 200 && (tdedx[0] + tdedx[1]) / 2 > 7e3 && (tdedx[0] + tdedx[1]) / 2 < 25000", "goff");
    pfY = h2EvsRun.ProfileY();

    TFile *file1 = TFile::Open(fKsKl.c_str());
    ksTr = (TTree *)file1->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);

    BadRunSearch();
}

std::vector<int> EnergyHandler::GetBadRunsList()
{ return badRuns; }

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
    int tmp = -1;
    nRunMax = 0;
    int entriesCounter = 0;
    for (int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if (std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
        {
            if (tmp != runnum || i == ksTr->GetEntriesFast() - 1)
            {
                if(entryNum.size() > 0)
                { runEntriesNum.push_back(entriesCounter); }

                entryNum.push_back(i);
                entriesCounter = 0;

                nRunMax++;
            }
            entriesCounter++;
        }
        tmp = runnum;
    }
    runEntriesNum.push_back(0);
}

size_t EnergyHandler::RunsDivider()
{
    // [0] - energy
    auto massErrFunc = new TF1("massErr1", "(1 - cos(x / 2) * cos(x / 2)) / TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    auto revMassFunc = new TF1("revMassFunc", "2*TMath::ACos(TMath::Sqrt((1 - x * x / [0] / [0])/(1-4*139.57 * 139.57 / [0] / [0])))");

    auto timeNow = std::time(nullptr);
    auto date = std::localtime(&timeNow);
    std::string rootFileName = "hists and root files/massMeas macro output/massHist" + std::to_string(date->tm_hour) + "h" + 
                                std::to_string(date->tm_min) + "m_"+std::to_string(date->tm_yday)+"day.root";
    TFile file2(rootFileName.c_str(), "recreate");
    std::vector<TH1D> massHistsMeas;

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
        { massHistsMeas.push_back(TH1D(("Meas" + std::to_string(mhCounter)).c_str(), "Ks mass", 150, 494, 510)); }
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
                massCrAngle->SetParameter(0, emeas);
                massHistsMeas[mhCounter].Fill(massCrAngle->Eval(ksdpsi));
                massCrAngle->SetParameter(0, e[m + j]);
                l++;
                ksTr->GetEntry(l);
            }
            j++;
        }
        if (1 / eErrBackSum > 0.00001 && runEntriesNum[m] != 0)
        {
            eMerge.push_back(eSum / j); eMergeErr.push_back(sqrt(1 / eErrBackSum));
            groupsStartRunnum.push_back(rNum[m]);
            for (int i = 0; i < 10; i++)
            {
                fermiStep->SetParameters(497.6, 0.335709, 66506, 399.187, 499.504, -9292.33);
                res = massHistsMeas[mhCounter].Fit("fermiStep", "SE", "", 495.0, 501);
                if (res->IsValid())
                { break; }
            }
            massEMeas.push_back(res->Parameter(0));
            ksTr->GetEntry(entryNum[m + j]);
            massErrFunc->SetParameter(0, eSum / j);
            revMassFunc->SetParameter(0, eSum / j);
            massEMeasErr.push_back(res->ParError(0));
            mhCounter++;
        }
        m += j;
    }

    massHistsMeas.pop_back();
    massHistsMeas.shrink_to_fit();
    nRunMax = nRunMax - 1;
    
    file2.Write();
    file2.Close();
    
    massHistsMeas.clear();
    delete fermiStep; delete massErrFunc; 
    delete revMassFunc; 
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

    std::vector<Float_t> eRatio;
    for(int i = 0; i < eMerge.size(); i++)
    {
        eMerge[i] += 3.8;
    }

    std::cout<<"Number of groups: " << groupsAmount <<std::endl;

    std::vector<Float_t> zeroes(groupsAmount, 0.0);

    TGraphErrors gEvsRunMeas(groupsAmount, groupsStartRunnum.data(), emeasGrouped.data(), zeroes.data(), demeasGrouped.data());
    TGraphErrors gEvsRunCalc(groupsAmount, groupsStartRunnum.data(), eMerge.data(), zeroes.data(), eMergeErr.data());
    TGraphErrors gMassMeasVsRun(groupsAmount, groupsStartRunnum.data(), massEMeas.data(), zeroes.data(), massEMeasErr.data());
    TGraph gRatio(groupsAmount-1, groupsStartRunnum.data(), eRatio.data());
    
    TMultiGraph *mg = new TMultiGraph();

    auto c1 = new TCanvas("EnergyStab", "EvsRunMeas, EvsRunCalc, MassVsRunMeas, MassVsRunCalc", 200, 10, 600, 400);
    //c1->Divide(2, 2);
    gEvsRunMeas.SetLineColor(kBlack);

    //c1->cd(1);
    gEvsRunMeas.SetTitle("E meas vs Run");
    gEvsRunMeas.GetXaxis()->SetTitle("Number of run");
    gEvsRunMeas.GetYaxis()->SetTitle("E_{beam}, MeV");
    //gEvsRunMeas.DrawClone("");

    //c1->cd(3);
    //gEvsRunCalc.SetTitle("E calc vs Run");
    gEvsRunCalc.SetMarkerColor(kBlue);
    //gEvsRunCalc.DrawClone("Same");
    
    mg->Add(&gEvsRunMeas);
    mg->Add(&gEvsRunCalc);
    mg->GetXaxis()->SetTitle("Number of run");
    mg->GetYaxis()->SetTitle("E_{beam}, MeV");
    mg->DrawClone("AP");

    delete mg;

/*
    c1->cd(2);
    gMassMeasVsRun.SetTitle("Mass (E meas) vs Run");
    gMassMeasVsRun.SetMarkerStyle(kFullDotLarge);
    gMassMeasVsRun.DrawClone("AP");
*/
}

double EnergyHandler::GetMassCriticalAngle()
{
    EnergyCalculation();
    EntryFilter();
    groupsAmount = RunsDivider();
    Draw();

    double bar = 0;
    double barErr = 0;
    for (int i = 0; i < massEMeas.size(); i++)
    {
        bar += massEMeas[i];
        if (massEMeasErr[i] != 0)
        { barErr += 1 / (massEMeasErr[i] * massEMeasErr[i]); }
        else
        { std::cout<<"massEMeasErr["<< i <<"] = 0" << std::endl; }
    }

    std::cout << "mass(emeas) = " << bar / massEMeas.size() << " +- " << 1 / sqrt(barErr) << std::endl;
    return bar / massEMeas.size();
}

void EnergyHandler::MassLnY(int drawOpt = 0)
{    
    auto hMlnY = new TH2D("hMlnY", "M(lnY)", 200, -0.4, 0.4, 200, 480, 520);
    auto hMPsi = new TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520);
    auto hM_CrAnglelnY = new TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 200, -0.3, 0.3, 200, 490, 515);
    auto hPsilnY = new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 200, 2.4, 3.3);

    auto hPsilnY1 = new TH2D("hPsilnY1", "Psi1(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY2 = new TH2D("hPsilnY2", "Psi2(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY3 = new TH2D("hPsilnY3", "Psi3(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY4 = new TH2D("hPsilnY4", "Psi4(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY5 = new TH2D("hPsilnY5", "Psi5(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY6 = new TH2D("hPsilnY6", "Psi6(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY7 = new TH2D("hPsilnY7", "Psi7(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);
    auto hPsilnY8 = new TH2D("hPsilnY8", "Psi8(lnY)", 400, -0.4, 0.4, 400, 2.4, 3.);

    auto hPsi = new TH1D("hPsi", "Psi", 200, 2.4, 3.1);

    double sigmaPsi = 0;
    double energy = 1;

    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9)
        {
            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massCrAngle->SetParameter(0, emeas);

            {
                if(abs(log(Y)) <= 0.05)
                { sigmaPsi = vSigma[0]; }
                else if (abs(log(Y)) <= 0.1)
                { sigmaPsi = vSigma[1]; }
                else if (abs(log(Y)) <= 0.15)
                { sigmaPsi = vSigma[2]; }
                else if (abs(log(Y)) <= 0.2)
                { sigmaPsi = vSigma[3]; }
                else if (abs(log(Y)) <= 0.25)
                { sigmaPsi = vSigma[4]; }
                else if (abs(log(Y)) <= 0.3)
                { sigmaPsi = vSigma[5]; }
                else if (abs(log(Y)) <= 0.35)
                { sigmaPsi = vSigma[6]; }
                else if (abs(log(Y)) <= 0.4)
                { sigmaPsi = vSigma[7]; }
            }

            // if(massFullRec->Eval(ksdpsi) > 490 && massFullRec->Eval(ksdpsi) < 505)
            // { hMlnY->Fill(log(Y), massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi)); }
            hMlnY->Fill(log(Y), massFullRec->Eval(ksdpsi));

            hM_CrAnglelnY->Fill(log(Y), massCrAngle->Eval(ksdpsi));
            if(massFullRec->Eval(ksdpsi) > 490 && massFullRec->Eval(ksdpsi) < 505)
            { 
                hPsilnY->Fill(log(Y), ksdpsi); 
                {
                    if(abs(log(Y)) <= 0.05)
                    { hPsilnY1->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.1)
                    { hPsilnY2->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.15)
                    { hPsilnY3->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.2)
                    { hPsilnY4->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.25)
                    { hPsilnY5->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.3)
                    { hPsilnY6->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.35)
                    { hPsilnY7->Fill(log(Y), ksdpsi); }
                    else if (abs(log(Y)) <= 0.4)
                    { hPsilnY8->Fill(log(Y), ksdpsi); }
                }  
            }

            hMPsi->Fill(ksdpsi, massFullRec->Eval(ksdpsi));
            hPsi->Fill(ksdpsi);
        }
    }
    std::cout << "!!!!!!! = " << massCrAngle->Eval(2.59085) << std::endl;
    TFitResultPtr r;
    std::vector<TH1D *> p;
    p.push_back(hPsilnY1->ProjectionY("py1"));
    p.push_back(hPsilnY2->ProjectionY("py2") );
    p.push_back(hPsilnY3->ProjectionY("py3") );
    p.push_back(hPsilnY4->ProjectionY("py4") );
    p.push_back(hPsilnY5->ProjectionY("py5") );
    p.push_back(hPsilnY6->ProjectionY("py6" ));
    p.push_back(hPsilnY7->ProjectionY("py7") );
    p.push_back(hPsilnY8->ProjectionY("py8") );
    std::cout << "Sigmas: " << std::endl;
    for(int i = 0; i < p.size(); i++)
    {
        r = p[i]->Fit("gaus", "SQE", "goff", p[i]->GetMean() - 0.04, p[i]->GetMean() + 0.04);
        std::cout << r->Parameter(2) << ", ";
    }
    std::cout<<std::endl;
    

    hMlnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnY->GetXaxis()->SetTitle("ln(Y)");

    r = hM_CrAnglelnY->ProfileX()->Fit("pol4", "SQE");
    std::cout << "Mass_CrAngle = " << r->Parameter(0) << " +/- " << r->ParError(0) << std::endl;
    r = hPsilnY->Fit("pol2", "SQE");
    std::cout << "Psi = " << r->Parameter(0) << " +/- " << r->ParError(0) << std::endl;
    std::cout << "emeas = " << emeas << std::endl;
    auto canv = new TCanvas("MlnY","Mass(lnY)", 200, 10, 600, 400);
    switch (drawOpt)
    {
    case 0:
        // hMlnY->GetXaxis()->SetRangeUser(-0.4, 0.4);
        // hMlnY->GetYaxis()->SetRangeUser(490., 505);
        // hMlnY->ProfileX("pfx")->DrawClone();
        hMlnY->DrawClone();
        break;
    case 1:
        hM_CrAnglelnY->ProfileX()->DrawClone();
        break;
    case 2:
        hPsilnY->DrawClone();
        break;
    case 3:
        hPsi->DrawClone();
        break;
    default:
        break;
    }

    delete hMlnY;
    delete hM_CrAnglelnY;
    delete hPsilnY;
    delete hMPsi;
    delete hPsi;

    delete hPsilnY1;
    delete hPsilnY2;
    delete hPsilnY3;
    delete hPsilnY4;
    delete hPsilnY5;
    delete hPsilnY6;
    delete hPsilnY7;
    delete hPsilnY8;
}

EnergyHandler::~EnergyHandler()
{
    delete pfY;
    delete ksTr;
    delete kTr;
    delete massCrAngle;
    delete massFullRec;
}

int massMeasRefactored()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    std::vector<Float_t> vSigma0 = {0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<Float_t> vSigma2bGen510 = {0.0161352, 0.0168335, 0.0174627, 0.0189506, 0.0203554, 0.0227744, 0.0256168, 0.0300882};
    std::vector<Float_t> vSigma2bGen511 = {0.0142948, 0.014196, 0.0139079, 0.0161541, 0.0177081, 0.0189663, 0.0221388, 0.0253829};
    std::vector<Float_t> vSigma2bGen514 = {1.37127e-02, 1.44228e-02, 1.45462e-02, 1.46740e-02, 1.72791e-02, 1.84825e-02,  2.11118e-02, 2.26650e-02};
    std::vector<Float_t> vSigmaMCGPJ505 = {1.39789e-02, 1.48976e-02, 1.63576e-02,  1.83967e-02, 2.14074e-02, 2.52082e-02, 2.92340e-02, 3.85429e-02};
    std::vector<Float_t> vSigmaMCGPJ510 = {0.0176344, 0.0187807, 0.0198675, 0.020309, 0.0222395, 0.024432, 0.027991, 0.0309448};
    std::vector<Float_t> vSigmaMCGPJ514 = {0.0423904, 0.0438116, 0.0445002, 0.0460169, 0.0480358, 0.052797, 0.0556817, 0.061061};
    
    //auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", "hists and root files/cuts/ksklCut_11May22.root");

    // auto eHandler = new EnergyHandler("hists and root files/cuts/kpkm_2bGen.root", "tr_ph/kskl2bGen514_withKLcut.root", vSigma0);
    // auto eHandler = new EnergyHandler("hists and root files/cuts/kpkm_2bGen.root", "tr_ph/tr_ph_kskl2bGen600k(diff theta cut).root", vSigma2bGen510);
    auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", "tr_ph/New folder/bonkXXL.root", vSigma0); // with Kl cut
    //auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", "tr_ph/tr_ph_ksklMCGPJ15Jul14h00m.root", vSigmaMCGPJ510); // with Kl cut
    
    //eHandler->GetMassCriticalAngle();
    eHandler->MassLnY(1);
    delete eHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}