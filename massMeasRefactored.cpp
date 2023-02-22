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
#include "Energy.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>

#include <chrono>
#include <ctime> 

class MassHandler
{
private:
    TTree *ksTr;

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


    std::vector<Float_t> vSigmaFit;
    std::vector<Float_t> vSigmaRMS;

    int nRunMax;
    double sigmaPsiCrAngle;
    double avgPsiCrAngle;
    double energyCorrected;

    std::vector<int> badRuns;
    // Number of groups of runs
    Int_t groupsAmount; 

    // Fills badRuns vector and returns number of good runs.
    int ReadBadRuns(std::string filename);
    int BadRunSearch();
    // Calculate sigma(lnY) via Fit and RMS
    void CalcSigmas();

public:
    /*
    * withFitSigma: true == use fit sigmas, false == use RMS;
    * drawOpt: 0 - M_FullRec vs lnY profile, 1 - M_CrAngle profile, 2 - Psi vs lnY, 
    * 3 - Psi distribution, 4 - EnergySpectrum (for MC only), 5 - M_FullRec vs lnY;
    */
    void MassLnY(std::map<int, Float_t> &energyDiff, bool isEnergyDiffCor, bool withFitSigma = true, int drawOpt = 0);
    std::vector<int> GetBadRunsList();

    MassHandler(std::string fKsKl, std::string badRunsFileName = "", double energyCorr = -1);
    ~MassHandler();
};

MassHandler::MassHandler(std::string fKsKl, std::string badRunsFileName, double energyCorr)
{
    TH2D h2EvsRun("dEdXvsPtot", "", 1000, 100, 120, 659, 60741, 61400);
    Int_t runnum_;
    Float_t tdedx[2];
    Float_t tptot[2];

    energyCorrected = energyCorr;

    // [0] = energy;
    massCrAngle = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57018 * 139.57018 / [0] / [0]) * cos(x / 2) * cos(x / 2))");

    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57 * 139.57 / E^2]))]
    massFullRec = new TF1("MassLnY", 
    "sqrt( [0]*[0] * (1 - 1 / [1] / [1] * (1 + cos(x) * sqrt(1-[1]*[1])) * (1 - sqrt(1-[1]*[1] * (1 - 4 * 139.57*139.57 / [0] / [0]) ) ) ) )");


    TFile *file1 = TFile::Open(fKsKl.c_str());
    ksTr = (TTree *)file1->Get("ksTree");
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("etrue", &etrue);
    etrue = -1;
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);
    
    ReadBadRuns(badRunsFileName);
    BadRunSearch();
}

std::vector<int> MassHandler::GetBadRunsList()
{ return badRuns; }

int MassHandler::ReadBadRuns(std::string filename)
{
    std::ifstream input(filename);
    if(input.is_open() && input.good())
    {
        std::istream_iterator<double> start(input), end;
        badRuns.insert(badRuns.end(), start, end);
    }
    return badRuns.size();
}

int MassHandler::BadRunSearch()
{
    int goodRunsCounter = 0;
    for(int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        // In later versions of tr_ph (after v8) if energy was not measured during specific run, then emeas0 == -1.
        // But right now (version 8) in this case emeas == stake energy.
        // if(emeas == *stake energy* && demeas == 0)
        if(false && (emeas == -1 || demeas == 0))
        { badRuns.push_back(runnum); }
    }

    return goodRunsCounter;
}

void MassHandler::CalcSigmas()
{
    std::vector<TH2D *> psilnYs;
    for(int i = 0; i < 16; i++)
    { psilnYs.push_back(new TH2D(("hPsilnY" + std::to_string(i + 1)).c_str(), ("Psi" + std::to_string(i + 1) + "(lnY)").c_str(), 400, -0.4, 0.4, 10000, 2.4, 3.3)); }

    double sigmaPsi = 0;
    double massFullRecWithEmeas = 0;
    double lnY = 0;
    double psiBinNum = 0;
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9)
        {
            lnY = log(Y);
            psiBinNum = lnY > -0.4 ? floor((lnY + 0.4 + 1e-12) / 0.05) : 16;
            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(ksdpsi);
    
            if(massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505 && psiBinNum < psilnYs.size())
            { psilnYs[psiBinNum]->Fill(lnY, ksdpsi); }
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
        vSigmaFit.push_back(r->Parameter(2));
    }
    std::cout<<std::endl;
    std::cout << "RMS: " << std::endl;
    for(auto hist : psilnYprojYs)
    {
        // hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()) - 3 * hist->GetStdDev(),
        //                                 hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()) + 3 * hist->GetStdDev());
        std::cout << hist->GetStdDev() << ", "; 
        vSigmaRMS.push_back(hist->GetStdDev());
    }
    std::cout<<std::endl;
}

void MassHandler::MassLnY(std::map<int, Float_t> &energyDiff, bool isEnergyDiffCor, bool withFitSigma = true, int drawOpt = 0)
{   
    CalcSigmas();

    auto hMlnY = new TH2D("hMlnY", "M(lnY)", 250, -1, 1, 40000, 480, 520);
    auto hDeltaM = new TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1);
    auto hMlnYpfx  = new TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505);
    auto hMPsi = new TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520);
    auto hM_CrAnglelnY = new TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 30, -0.4, 0.4, 40000, 490, 515);
    auto hPsilnY = new TH2D("hPsilnY", "Psi(lnY)", 250, -0.5, 0.5, 10000, 2.4, 3.3);
    auto hEnergySpectrum = new TH1D("hEnergySpectrum", "", 1000, emeas - 20, emeas + 5);
    // After profile cut
    auto hEnergySpectrumCut = new TH1D("hEnergySpectrumCut", "", 1000, emeas - 20, emeas + 5);
    auto hPsi = new TH1D("hPsi", "Psi", 200, 2.4, 3.3);

    double sigmaPsi = 0;
    double energy = 0;
    double massFullRecWithEmeas = 0;
    double lnY = 0;
    double psiBinNum = 0;
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9)
        {
            lnY = log(Y);
            psiBinNum = lnY > -0.4 ? floor((lnY + 0.4 + 1e-12) / 0.05) : 16;
            if(psiBinNum < vSigmaFit.size())
            { sigmaPsi = withFitSigma? vSigmaFit[psiBinNum] : vSigmaRMS[psiBinNum]; }

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi);

            emeas = (energyCorrected == -1) ? emeas :  energyCorrected + (isEnergyDiffCor? 1 : 0) * energyDiff[runnum];
            if(runnum >= 60790 && runnum <= 60921)
            { emeas = 15.7437 / (runnum-60763.1) / (runnum-60763.1) + 509.518; }

            if(runnum >= 60922 && runnum <= 61174)
            { emeas = 2441.92 / (runnum-60737.2) / (runnum-60737.2) + 509.497; }

            if(runnum >= 61175 && runnum <= 61379)
            { emeas = 631.041 / (runnum-61094.5) / (runnum-61094.5) + 509.527; }

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massCrAngle->SetParameter(0, emeas);
        
            hMlnY->Fill(lnY, massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
            hM_CrAnglelnY->Fill(lnY, massCrAngle->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massCrAngle->Derivative2(ksdpsi));
            // hPsilnY->Fill(lnY, ksdpsi); 

            if(massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
            { 
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

    hMlnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnY->GetXaxis()->SetTitle("ln(Y)");
    hMlnYpfx->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnYpfx->GetXaxis()->SetTitle("ln(Y)");

    TFitResultPtr r;
    r = hM_CrAnglelnY->ProfileX()->Fit("pol4", "SMQE", "", -0.2, 0.2);
    std::cout << "Mass_CrAngle = " << r->Parameter(0) << " +/- " << r->ParError(0) 
                << "; chi2 / ndf = " << r->Chi2() << "/" << r->Ndf() << "; Prob = " << r->Prob() << std::endl;
    r = hMlnYpfx->Fit("pol0", "SMQE", "goff", -0.3, 0.3);
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
        massLine->DrawClone("same");
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
        hEnergySpectrum->DrawClone();
        hEnergySpectrumCut->DrawClone("same E0");
        break;
    case 5:
        hMlnY->DrawClone();
    default:
        break;
    }

    delete hMlnY;
    delete hM_CrAnglelnY;
    delete hPsilnY;
    delete hMPsi;
    delete hPsi;

    delete massCutHorizontal1;
    delete massCutHorizontal2;
    delete massCutVertical1;
    delete massCutVertical2;
}

MassHandler::~MassHandler()
{
    delete ksTr;
    delete massCrAngle;
    delete massFullRec;
}

int massMeasRefactored()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    std::string fileName = "tr_ph/expKsKl/exp509.5_v9.root";
    // std::string fileName = "tr_ph/MC/MC514_v9.root";
    auto kchEnergyHandler = new Energy("tr_ph/expKpKm/kchExp510.root", 509.957, 0.005, 15, 3.709);
    // auto kchEnergyHandler = new Energy("tr_ph/expKpKm/kchExp509.5.root", 509.528, 0.01, 15, 3.751);
    auto energyDiff = kchEnergyHandler->GetEnergyDiff();
    delete kchEnergyHandler;
    auto massHandler = new MassHandler(fileName, "txt/BadRuns.txt");
    massHandler->MassLnY(energyDiff, false, true, 5);
    delete massHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}