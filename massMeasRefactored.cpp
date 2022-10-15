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
#include "TROOT.h"
#include "TLine.h"

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
    // Kaon energy calculated as invariant mass of pi+pi-.
    Float_t etrue;
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
    double energyCorrected;

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
    EnergyHandler(std::string fChargedK, std::string fKsKl, std::vector<Float_t> sigmas, double energyCorr = -1);
    void MassLnY(int drawOpt = 0);
    double GetMassCriticalAngle();

    std::vector<int> GetBadRunsList();

    ~EnergyHandler();
};

EnergyHandler::EnergyHandler(std::string fChargedK, std::string fKsKl, std::vector<Float_t> sigmas, double energyCorr = -1)
{
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree *)file->Get("kChargedTree");
    TH2D h2EvsRun("dEdXvsPtot", "", 1000, 100, 120, 659, 60741, 61400);
    Int_t runnum_;
    Float_t tdedx[2];
    Float_t tptot[2];

    vSigma = sigmas;
    energyCorrected = energyCorr;

    // [0] = energy;
    massCrAngle = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57018 * 139.57018 / [0] / [0]) * cos(x / 2) * cos(x / 2))");

    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57 * 139.57 / E^2]))]
    massFullRec = new TF1("MassLnY", 
    "sqrt( [0]*[0] * (1 - 1 / [1] / [1] * (1 + cos(x) * sqrt(1-[1]*[1])) * (1 - sqrt(1-[1]*[1] * (1 - 4 * 139.57*139.57 / [0] / [0]) ) ) ) )");

    kTr->SetBranchAddress("runnum", &runnum_);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);

    kTr->Draw("runnum :(tptot[0] + tptot[1]) / 2 >> dEdXvsPtot",
              "(tptot[0] + tptot[1]) / 2 > 100 && (tptot[0] + tptot[1]) / 2 < 200 && (tdedx[0] + tdedx[1]) / 2 > 7e3 && (tdedx[0] + tdedx[1]) / 2 < 25000", "goff");
    pfY = h2EvsRun.ProfileY();

    TFile *file1 = TFile::Open(fKsKl.c_str());
    ksTr = (TTree *)file1->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("etrue", &etrue);
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
        // In later versions of tr_ph (after v8) if energy was not measured during specific run, then emeas == -1.
        // But right now (version 8) in this case emeas == stake energy.
        // if(emeas == *stake energy* && demeas == 0 && std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
        if(emeas == -1 && std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end())
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
    
    auto hMlnY = new TH2D("hMlnY", "M(lnY)", 250, -1, 1, 40000, 480, 520);
    auto hMlnYpfx  = new TProfile("hMlnYpfx","Profile of M versus lnY", 40, -1, 1, 490, 505);
    auto hMPsi = new TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520);
    auto hM_CrAnglelnY = new TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 200, -0.4, 0.4, 40000, 490, 515);
    auto hPsilnY = new TH2D("hPsilnY", "Psi(lnY)", 250, -0.5, 0.5, 10000, 2.4, 3.3);
    auto hEnergySpectrum = new TH1D("hEnergySpectrum", "", 1000, emeas - 20, emeas + 5);
    // After profile cut
    auto hEnergySpectrumCut = new TH1D("hEnergySpectrumCut", "", 1000, emeas - 20, emeas + 5);

    std::vector<TH2D *> psilnYs;
    for(int i = 0; i < 8; i++)
    { psilnYs.push_back(new TH2D(("hPsilnY" + std::to_string(i + 1)).c_str(), ("Psi" + std::to_string(i + 1) + "(lnY)").c_str(), 400, -0.4, 0.4, 10000, 2.4, 3.3)); }

    auto hPsi = new TH1D("hPsi", "Psi", 200, 2.4, 3.3);

    double sigmaPsi = 0;
    double energy = 1;

    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9)
        {
            emeas = (energyCorrected == -1) ? emeas : energyCorrected;
            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massCrAngle->SetParameter(0, emeas);
            if(int((abs(log(Y)) + 1e-7) / 0.05) < vSigma.size())
            { sigmaPsi = vSigma[int((abs(log(Y)) + 1e-7) / 0.05)]; }

            hMlnY->Fill(log(Y), massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
            // hMlnYpfx->Fill(log(Y), massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
        
            hM_CrAnglelnY->Fill(log(Y), massCrAngle->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massCrAngle->Derivative2(ksdpsi));
            hPsilnY->Fill(log(Y), ksdpsi); 
            if(int((abs(log(Y)) + 1e-12) / 0.05) < psilnYs.size())
            { psilnYs[int((abs(log(Y)) + 1e-12) / 0.05)]->Fill(log(Y), ksdpsi); }

            if(massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi) > 490 && 
                massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi) < 505)
            { 
                hMlnYpfx->Fill(log(Y), massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
                // hPsilnY->Fill(log(Y), ksdpsi); 
                // if(int((abs(log(Y)) + 1e-12) / 0.05) < psilnYs.size())
                // { psilnYs[int((abs(log(Y)) + 1e-12) / 0.05)]->Fill(log(Y), ksdpsi); }  
                if(fabs(log(Y)) < 0.3)
                { hEnergySpectrumCut->Fill(etrue); }
            }

            hMPsi->Fill(ksdpsi, massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
            hPsi->Fill(ksdpsi);
            if(fabs(log(Y)) < 0.2)
            { hEnergySpectrum->Fill(etrue); }
        }
    }

    TFitResultPtr r;
    std::vector<TH1D *> psilnYprojYs;
    std::cout << "Sigmas: " << std::endl;
    for(int i = 0; i < psilnYs.size(); i++)
    {
        psilnYprojYs.push_back(psilnYs[i]->ProjectionY(("py" + std::to_string(i + 1)).c_str()));
        psilnYprojYs.back()->Rebin(40);
        r = psilnYprojYs.back()->Fit("gaus", "SEQ", "", psilnYprojYs.back()->GetXaxis()->GetBinCenter(psilnYprojYs.back()->GetMaximumBin()) - 3 * psilnYprojYs.back()->GetStdDev(), 
                                                psilnYprojYs.back()->GetXaxis()->GetBinCenter(psilnYprojYs.back()->GetMaximumBin()) + 3 * psilnYprojYs.back()->GetStdDev());
        r = psilnYprojYs.back()->Fit("gaus", "SEQ", "", r->Parameter(1) - 2 * r->Parameter(2), r->Parameter(1) + 2 * r->Parameter(2));
        std::cout << r->Parameter(2) << ", ";
    }
    std::cout<<std::endl;
    std::cout << "RMS: " << std::endl;
    for(auto hist : psilnYprojYs)
    { std::cout << hist->GetStdDev() << ", "; }
    std::cout<<std::endl;

    hMlnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnY->GetXaxis()->SetTitle("ln(Y)");
    hMlnYpfx->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnYpfx->GetXaxis()->SetTitle("ln(Y)");

    r = hM_CrAnglelnY->ProfileX()->Fit("pol4", "SMQE", "", -0.2, 0.2);
    std::cout << "Mass_CrAngle = " << r->Parameter(0) << " +/- " << r->ParError(0) 
                << "; chi2 / ndf = " << r->Chi2() << "/" << r->Ndf() << "; Prob = " << r->Prob() << std::endl;
    r = hPsilnY->ProfileX("pfxAng")->Fit("pol2", "SQME", "", -0.2, 0.2);
    std::cout << "Psi = " << r->Parameter(0) << " +/- " << r->ParError(0) << std::endl;
    auto canv = new TCanvas("MlnY","Mass(lnY)", 200, 10, 600, 400);
    std::cout << "Average E cut = " << hEnergySpectrumCut->GetMean() << " + / - " << hEnergySpectrumCut->GetMeanError() << std::endl;
    std::cout << "Average E  = " << hEnergySpectrum->GetMean() << " + / - " << hEnergySpectrum->GetMeanError() << std::endl;

    auto massCutHorizontal1 = new TLine(-0.8, 490, 0.8, 490);
    massCutHorizontal1->SetLineColor(kBlue);
    massCutHorizontal1->SetLineWidth(2);

    auto massCutHorizontal2 = new TLine(-0.8, 505, 0.8, 505);
    massCutHorizontal2->SetLineColor(kBlue);
    massCutHorizontal2->SetLineWidth(2);

    auto massCutVertical1 = new TLine(-0.4, 490, -0.4, 505);
    massCutVertical1->SetLineColor(kBlue);
    massCutVertical1->SetLineWidth(2);

    auto massCutVertical2 = new TLine(0.4, 490, 0.4, 505);
    massCutVertical2->SetLineColor(kBlue);
    massCutVertical2->SetLineWidth(2);

    auto massLine = new TLine(-0.3, 497.614, 0.3, 497.614);
    massLine->SetLineColor(kBlue);
    massLine->SetLineWidth(2);

    hM_CrAnglelnY->GetXaxis()->SetTitle("lnY");
    hM_CrAnglelnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");

    auto MCrAngleProf = hM_CrAnglelnY->ProfileX();
    MCrAngleProf->GetXaxis()->SetTitle("lnY");
    MCrAngleProf->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");

    switch (drawOpt)
    {
    case 0:
        hMlnYpfx->GetXaxis()->SetRangeUser(-0.5, 0.5);
        hMlnYpfx->DrawClone();
        massLine->Draw("same");
        break;
    case 1:
        MCrAngleProf->Fit("pol4", "SMQE", "", -0.2, 0.2);
        MCrAngleProf->DrawClone();
        break;
    case 2:
        hPsilnY->DrawClone();
        break;
    case 3:
        hPsi->DrawClone();
        break;
    case 4:
        hEnergySpectrumCut->SetLineColor(kBlue);
        hEnergySpectrum->Draw();
        hEnergySpectrumCut->Draw("same");
        break;
    default:
        break;
    }

    // delete hMlnY;
    // delete hM_CrAnglelnY;
    // delete hPsilnY;
    // delete hMPsi;
    // delete hPsi;

    delete massCutHorizontal1;
    delete massCutHorizontal2;
    delete massCutVertical1;
    delete massCutVertical2;
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

    std::vector<Float_t> vSigma0(8, 0);
    std::vector<Float_t> vSigmaMCGPJ = {0.0148821, 0.0151739, 0.0163453, 0.0178507, 0.0199841, 0.0231457, 0.0260142, 0.0300946};
    std::vector<Float_t> vSigmaExp509_5 = {0.0131431, 0.0148787, 0.014976, 0.0162632, 0.0167996, 0.0211273, 0.0235995, 0.0395321};
    std::vector<Float_t> vSigmaFit514MC = {0.0141507, 0.0145647, 0.0164458, 0.0163694, 0.0183365, 0.0201849, 0.0228716, 0.0266007};
    std::vector<Float_t> vSigmaFit = {0.0142348, 0.0147611, 0.0157042, 0.0187489, 0.0183891, 0.0225073, 0.0255563, 0.0311781};
    std::vector<Float_t> vSigmaRMS = {0.045531, 0.0475883, 0.0484802, 0.0514004, 0.051982, 0.0571173, 0.060207, 0.0640742};
    
    //auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", "hists and root files/cuts/ksklCut_11May22.root");

    auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", "tr_ph/MC508tmp.root", vSigmaFit, 507.908);
    // auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", "tr_ph/exp509_5_newtry.root", vSigmaExp509_5);
    eHandler->MassLnY(1);
    delete eHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}