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
#include <optional>

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
    std::optional<double> energyCorrected;
    std::optional<std::map<int, Float_t>> energyDiff;

    std::vector<int> badRuns;
    // Number of groups of runs
    Int_t groupsAmount; 

    // Fills badRuns vector and returns number of good runs.
    int ReadBadRuns(std::string filename);
    int BadRunSearch();
    // Calculate sigma(lnY) via Fit and RMS
    void CalcSigmas();
    /* Return either corrected energy for a given run according to a polynomial func with parameters from fit if there is one
    * or average energy / emeas for this run.
    */
    Float_t GetRightEnergy(int rnum);

public:
    /*
    * withFitSigma: true == use fit sigmas, false == use RMS;
    * drawOpt: 0 - M_FullRec vs lnY profile, 1 - M_CrAngle profile, 2 - Psi vs lnY, 
    * 3 - Psi distribution, 4 - EnergySpectrum (for MC only), 5 - M_FullRec vs lnY;
    * Returns mass from the fit. 
    */
    double MassLnY(double fitRange = 0.33, bool useEtrue = false, bool withFitSigma = true, int drawOpt = 0);
    std::vector<int> GetBadRunsList();

    MassHandler(std::string fKsKl, std::string badRunsFileName = "", 
                std::optional<double> energyCorr = std::nullopt, 
                std::optional<std::map<int, Float_t>> energyDiff = std::nullopt);
    ~MassHandler();
};

MassHandler::MassHandler(std::string fKsKl, std::string badRunsFileName, std::optional<double> energyCorr, std::optional<std::map<int, Float_t>> energyDiff)
{
    TH2D h2EvsRun("dEdXvsPtot", "", 1000, 100, 120, 659, 60741, 61400);
    Int_t runnum_;
    Float_t tdedx[2];
    Float_t tptot[2];

    energyCorrected = energyCorr;
    this->energyDiff = energyDiff;

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

Float_t MassHandler::GetRightEnergy(int rnum)
{
    Float_t energy = energyCorrected.value_or(emeas) + (energyDiff.has_value() ? energyDiff->operator[](rnum) : 0);
    
    // E = 508.5 MeV
    if(rnum >= 60196 && rnum <= 60259)
    { energy = 35.83 / (rnum - 60182.1) / (rnum - 60182.1) + 508.393; }

    if(rnum >= 60260 && rnum <= 60416)
    { energy = 769.991 / (rnum-60153.3) / (rnum-60153.3) + 508.361; }

    if(rnum >= 60417 && rnum <= 60498)
    { energy = 174994 / (rnum-59094.7) / (rnum-59094.7) + 508.288; }

    // E = 509 MeV
    if(rnum >= 60520 && rnum <= 60702)
    { energy = 0.0076 * sin(0.0950 * (rnum - 60584.6)) + 508.945; }

    // E = 509.5 MeV
    if(rnum >= 60790 && rnum <= 60921)
    { energy = 15.7437 / (rnum-60763.1) / (rnum-60763.1) + 509.518; }

    if(rnum >= 60922 && rnum <= 61174)
    { energy = 2441.92 / (rnum-60737.2) / (rnum-60737.2) + 509.497; }

    if(rnum >= 61175 && rnum <= 61378)
    { energy = 631.041 / (rnum-61094.5) / (rnum-61094.5) + 509.527; }
    
    // E = 510 MeV
    if(rnum >= 61380 && rnum <= 61461)
    { energy = 364.312 / (rnum - 61332.6) / (rnum - 61332.6) + 509.925; }

    if(rnum >= 61560 && rnum <= 61689)
    { energy = 598.053 / (rnum - 61479.6) / (rnum - 61479.6) + 509.928; }

    if(rnum >= 61689 && rnum <= 61856)
    { energy = 0.0101 * sin(0.0397 * (rnum-61704.6)) + 509.960; }

    // E = 510.5 MeV
    if(rnum >= 61859 && rnum <= 61958)
    { energy = 1587 / (rnum-61683) / (rnum-61683) + 510.434; }

    if(rnum >= 61958 && rnum <= 62075)
    { energy = -0.0109 * sin(0.068 * (rnum-62018)) + 510.448; }

    // E = 511 MeV
    if(rnum >= 62108 && rnum <= 62167)
    { energy = 50.2021 / (rnum-62077.4) / (rnum-62077.4) + 511.014; }

    if(rnum >= 62194 && rnum <= 62311)
    { energy = 14.23 / (rnum-62179) / (rnum-62179) + 511.031; }

    return energy;
}

double MassHandler::MassLnY(double fitRange, bool useEtrue, bool withFitSigma, int drawOpt )
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

    auto sigmas = std::vector<double>({0.0303571, 0.0244383, 0.0230856, 0.0203081, 0.0182954, 0.0188669, 0.0167881, 0.0151023, 0.0154767, 0.0172125, 0.0169003, 0.0184711, 0.0220203, 0.0231788, 0.0256208, 0.0287016});
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
            // { sigmaPsi = sigmas[psiBinNum]; }

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi);

            emeas = useEtrue? etrue : GetRightEnergy(runnum);

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
                if(fabs(lnY) < fitRange + 1e-6)
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
    r = hMlnYpfx->Fit("pol0", "SMQE", "goff", -fitRange, fitRange);
    double mass = r->Parameter(0);
    std::cout << "Mass_FullRec = " << r->Parameter(0) << " +/- " << r->ParError(0) 
                << "; chi2 / ndf = " << r->Chi2() << "/" << r->Ndf() << "; Prob = " << r->Prob() << std::endl;
    r = hPsilnY->ProfileX("pfxAng")->Fit("pol2", "SQME", "", -0.2, 0.2);
    std::cout << "Psi = " << r->Parameter(0) << " +/- " << r->ParError(0) << std::endl;
    auto canv = new TCanvas("MlnY","Mass(lnY)", 200, 10, 600, 400);
    std::cout << "Average E cut = " << hEnergySpectrumCut->GetMean() << " +/- " << hEnergySpectrumCut->GetMeanError() << std::endl;
    std::cout << "Average E  = " << hEnergySpectrum->GetMean() << " +/- " << hEnergySpectrum->GetMeanError() << std::endl;

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
        hMlnY->DrawClone("col");
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

    return mass;
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

    const std::vector<double> meanEnergies_vec = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};
    const std::vector<double> deltaM_RC_Smeared = {0.132, 0.086, 0.076, 0.062, 0.071, 0.106, 0.178, 0.325, 1.453};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    std::map<std::string, std::pair<double, double>> meanEnergies;
    std::map<std::string, double> radiativeCorrections;
    for(int i = 0; i < energyPoints.size(); i++)
    { 
        meanEnergies[energyPoints[i]] = std::make_pair(meanEnergies_vec[i], meanEnergiesErr[i]); 
        radiativeCorrections[energyPoints[i]] = deltaM_RC_Smeared[i]; 
    }

    std::string energyPoint = "510.5";
    // std::string fileName = "tr_ph/expKsKl/exp" + energyPoint + "_v9.root";
    // std::string fileName = "tr_ph/MC/MC510.5_v9.root";
    std::string fileName = "tr_ph/MC/MC510.5_Smeared.root";
    // std::string fadRunsFile = "tr_ph/MC/MC509.5_Smeared.root";

    auto kchEnergyHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoint + ".root", "txt/BadRuns.txt", 
                                                        meanEnergies[energyPoint].first, meanEnergies[energyPoint].second, 30, 0); 
    auto energyDiff = kchEnergyHandler->GetEnergyDiff();
    delete kchEnergyHandler;

    auto massHandler = new MassHandler(fileName, "txt/BadRuns.txt", meanEnergies[energyPoint].first);
    // auto massHandler = new MassHandler(fileName, "", 510.694);
    auto mass = massHandler->MassLnY(0.33, false) - radiativeCorrections[energyPoint];
    std::cout << "M_NCRC_smeared = " << mass << std::endl;
    delete massHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}