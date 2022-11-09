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
    etrue = -1;
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
    auto hDeltaM = new TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1);
    auto hMlnYpfx  = new TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505);
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
    double energy = 0;
    double massFullRecWithEmeas = 0;
    double lnY = 0;
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9)
        {
            lnY = log(Y);
            if(int((abs(lnY) + 1e-7) / 0.05) < vSigma.size())
            { sigmaPsi = vSigma[int((abs(lnY) + 1e-7) / 0.05)]; }

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi);

            emeas = (energyCorrected == -1) ? emeas : energyCorrected;
            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massCrAngle->SetParameter(0, emeas);
        
            hMlnY->Fill(lnY, massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
            hM_CrAnglelnY->Fill(lnY, massCrAngle->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massCrAngle->Derivative2(ksdpsi));
            // hPsilnY->Fill(lnY, ksdpsi); 

            // if(int((abs(lnY) + 1e-12) / 0.05) < psilnYs.size())
            // { psilnYs[int((abs(lnY) + 1e-12) / 0.05)]->Fill(lnY, ksdpsi); }

            // if((massFullRecWithEmeas > 490 && massFullRecWithEmeas < 496) || (massFullRecWithEmeas > 500 && massFullRecWithEmeas < 505))
            // if(massFullRecWithEmeas > 496 && massFullRecWithEmeas < 500)
            if(massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
            { 
                if(int((abs(lnY) + 1e-12) / 0.05) < psilnYs.size())
                { psilnYs[int((abs(lnY) + 1e-12) / 0.05)]->Fill(lnY, ksdpsi); }
                hDeltaM->Fill(lnY, - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
                hMlnYpfx->Fill(lnY, massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
                hPsilnY->Fill(lnY, ksdpsi); 
                if(fabs(lnY) < 0.34)
                { hEnergySpectrumCut->Fill(etrue); }
            }

            hMPsi->Fill(ksdpsi, massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
            hPsi->Fill(ksdpsi);
            hEnergySpectrum->Fill(etrue);
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
    {
        // hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()) - 3 * hist->GetStdDev(),
        //                                 hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()) + 3 * hist->GetStdDev());
        std::cout << hist->GetStdDev() << ", "; 
    }
    std::cout<<std::endl;

    hMlnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnY->GetXaxis()->SetTitle("ln(Y)");
    hMlnYpfx->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnYpfx->GetXaxis()->SetTitle("ln(Y)");

    r = hM_CrAnglelnY->ProfileX()->Fit("pol4", "SMQE", "", -0.2, 0.2);
    std::cout << "Mass_CrAngle = " << r->Parameter(0) << " +/- " << r->ParError(0) 
                << "; chi2 / ndf = " << r->Chi2() << "/" << r->Ndf() << "; Prob = " << r->Prob() << std::endl;
    r = hMlnYpfx->Fit("pol0", "SMQE", "goff", -0.34, 0.34);
    std::cout << "Mass_FullRec = " << r->Parameter(0) << " +/- " << r->ParError(0) 
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
        hEnergySpectrumCut->SetMarkerColor(kBlue);
        hEnergySpectrum->Draw();
        hEnergySpectrumCut->Draw("same E0");
        break;
    case 5:
        hMlnY->DrawClone();
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
    std::vector<Float_t> vSigmaFit514MC = {0.0141507, 0.0145647, 0.0164458, 0.0163694, 0.0183365, 0.0201849, 0.0228716, 0.0266007};
    std::vector<Float_t> vSigmaFit = {0.0150457, 0.0157218, 0.0167191, 0.0180674, 0.0199512, 0.0224721, 0.0259625, 0.0309317};
    std::vector<Float_t> vSigmaRMS = {0.0210113, 0.0216181, 0.0220604, 0.0237359, 0.0255735, 0.0279509, 0.0309818, 0.0349982,};
    
    // std::vector<Float_t> vSigmaFit514MC2 = {0.0139102, 0.0159495, 0.0161846, 0.0163694, 0.0183365, 0.0201849, 0.0228716, 0.0266007};

    std::string fileName = "tr_ph/MC514tmp.root";
    // std::string fileName = "tr_ph/exp509_5Cowboy.root";
    // std::string fileName = "tr_ph/exp509_MissinMassCut1.root";
    double energy = 512.393;
    // auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", fileName, vSigma0);
    auto eHandler = new EnergyHandler("hists and root files/cuts/kchCut21May.root", fileName, vSigmaFit514MC, energy);
    eHandler->MassLnY(5);
    delete eHandler;

    // std::vector<double> y = {5, 7, 9, 10, 16, 53, 99, 142, 138, 138, 131, 153, 146, 146, 142, 169, 116, 119, 139, 143, 150, 133, 143, 138, 111, 85, 37, 16, 6, 6, 6};
    // std::vector<double> x = {-30, -28, -26, -24, -22, -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    // std::vector<double> yErr = {2.23606797749979, 2.6457513110645907, 3.0, 3.1622776601683795, 4.0, 7.280109889280518, 9.9498743710662, 11.916375287812984, 11.74734012447073, 11.74734012447073, 11.445523142259598, 12.36931687685298, 12.083045973594572, 12.083045973594572, 11.916375287812984, 13.0, 10.770329614269007, 10.908712114635714, 11.789826122551595, 11.958260743101398, 12.24744871391589, 11.532562594670797, 11.958260743101398, 11.74734012447073, 10.535653752852738, 9.219544457292887, 6.082762530298219, 4.0, 2.449489742783178, 2.449489742783178, 2.449489742783178};
    // std::vector<double> xErr(yErr.size(), 0.0);
    // auto gr = new TGraphErrors(y.size(), x.data(), y.data(), xErr.data(), yErr.data());
    // gr->DrawClone();

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}