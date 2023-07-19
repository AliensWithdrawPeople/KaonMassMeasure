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
#include <memory>

#include <chrono>
#include <ctime> 

class MassHandler
{
private:
    std::unique_ptr<TTree> ksTr;

    std::unique_ptr<TF1> massCrAngle;
    std::unique_ptr<TF1> massFullRec;

    Float_t emeas; Float_t demeas;
    // Kaon energy calculated as invariant mass of pi+pi-.
    Float_t etrue;
    Float_t ksTheta;
    Float_t ksPhi;
    Float_t piThetaPos;
    Float_t piThetaNeg;
    Int_t runnum; 
    Float_t ksdpsi;
    Int_t nhitPos;
    Int_t nhitNeg;
    // Momentum ratio = P1/P2, where P1 is the momentum of pi+, P2 is the momentum of pi-.
    Float_t Y;

    std::vector<Float_t> vSigmaFit;
    std::vector<Float_t> vSigmaRMS;
    std::vector<Float_t> vY_RMS;
    std::vector<int> badRuns;

    std::unique_ptr<TH2D> hMlnY;
    std::unique_ptr<TH2D> hMPsi;
    std::unique_ptr<TH2D> hM_CrAnglelnY;
    std::unique_ptr<TH2D> hPsilnY;
    std::unique_ptr<TH2D> hMassVsKsTheta;
    std::unique_ptr<TH2D> hKsThetaVsLnY;

    std::unique_ptr<TProfile> hDeltaM;
    std::unique_ptr<TProfile> hMlnYpfx;

    std::unique_ptr<TH1D> hEnergySpectrum;
    // After profile cut
    std::unique_ptr<TH1D> hEnergySpectrumCut;
    std::unique_ptr<TH1D> hPsi;

    std::unique_ptr<TGraphErrors> grMassLnYFit;
    std::unique_ptr<TGraphErrors> grPsiLnYFit;


    std::optional<double> energyCorrected;
    std::optional<std::map<int, Float_t>> energyDiff;


    // Fills badRuns vector and returns number of good runs.
    int ReadBadRuns(std::string filename);
    int BadRunSearch();
    // Calculate sigma(lnY) (obtained from fit and as RMS).
    void CalcSigmas(bool verbose = true);
    void FillHists(double fitRange, double deltaE, bool withFitSigma, bool useEtrue);

    /* Return either corrected energy for a given run according to a polynomial func with parameters from fit if there is one
    * or average energy / emeas for this run.
    * Also it adds a deltaE correction (Ebergy = Energy0 + deltaE).
    */
    Float_t GetCorrectedEnergy(int rnum, double deltaE = 0.0);

    /*
    * This function takes 2D histogram, slices it along X-axis to nSlices, rebin obtained slices by rebinFactor  
    * and fit each slice with gaussian function in range (Mean - fitRange.first * Sigma, Mean + fitRange.second * Sigma) 
    * in two steps: at the firs step Mean = Slice.Mean, Sigma = Slice.RMS, at the second step it uses results of the previous fit. 
    * If isVerbose == true, results of fits will be shown. 
    * If saveSlices == true, the slices will be saved in the ROOT context (TDirectory), hence won't be deleted.
    */
    static std::unique_ptr<TGraphErrors> FitSlices(const std::unique_ptr<TH2D> &hist, int nSlices, std::pair<double, double> fitRange = std::make_pair(1.5, 0.8), 
                                                    std::optional<int> rebinFactor = std::nullopt, bool isVerbose = false);
    /* 
    * Draws graphs and hists related to Ks mass measurement.
    * drawOpt: 0 - M_FullRec vs lnY profile, 1 - M_CrAngle profile, 2 - Psi vs lnY, 
    * 3 - Psi distribution, 4 - EnergySpectrum (for MC only), 5 - M_FullRec vs lnY;
    */
    void DrawGraphs(int drawOpt);

public:

    /*
    * withFitSigma: true == use fit sigmas, false == use RMS.
    * drawOpt:      
    * 0 - M_FullRec vs lnY profile, 
    * 1 - M_CrAngle profile, 
    * 2 - Psi vs lnY, 
    * 3 - Psi distribution, 
    * 4 - EnergySpectrum (for MC only), 
    * 5 - M_FullRec vs lnY.
    * 
    * methodOpt:    
    * 0 - Mass_FullRec
    * 1 - MassVsLnY_Fit_FullRec
    * 2 - Mass_CrAngle
    * 3 - massCrAngle_Psi
    * 4 - massCrAngle_Psi_Fit.
    * 
    * Returns pair(mass, mass_error). 
    */
    std::pair<double, double> GetMass(double fitRange = 0.33, bool useEtrue = false, double deltaE = 0.0, bool withFitSigma = true, int drawOpt = 0, int methodOpt = 0, bool verbose = true);

    std::vector<int> GetBadRunsList();


    MassHandler(std::string fKsKl, std::string badRunsFileName = "", 
                std::optional<double> energyCorr = std::nullopt, 
                std::optional<std::map<int, Float_t>> energyDiff = std::nullopt,
                bool isVerbose = true);
    ~MassHandler();
};

MassHandler::MassHandler(std::string fKsKl, std::string badRunsFileName, std::optional<double> energyCorr, std::optional<std::map<int, Float_t>> energyDiff, bool isVerbose)
{
    energyCorrected = energyCorr;
    this->energyDiff = energyDiff;

    // [0] = energy;
    massCrAngle = std::make_unique<TF1>(TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57039 * 139.57039 / [0] / [0]) * cos(x / 2) * cos(x / 2))"));

    // [0] = E; [1] = eta = (1 - Y^2) / (1 + Y^2)
    // M = sqrt(E^2[1- 1/eta^2 (1+sqrt(1-eta^2)cosx) * (1 - sqrt(1-eta^2 * [1 - 4 * 139.57039 * 139.57039 / E^2]))]
    // massFullRec = new TF1("MassLnY", 
    // "sqrt( [0]*[0] * (1 - 1 / [1] / [1] * (1 + cos(x) * sqrt(1-[1]*[1])) * (1 - sqrt(1-[1]*[1] * (1 - 4 * 139.57039*139.57039 / [0] / [0]) ) ) ) )");

    massFullRec = std::make_unique<TF1>( TF1("MassLnY", 
    "sqrt( [0]*[0] * (1 - 1 / [1] / [1] * (1 + cos(x) * sqrt(1-[1]*[1])) * (1 - sqrt(1-[1]*[1] * (1 - 4 * 139.57039*139.57039 / [0] / [0]) ) ) ) )"));

    TFile *file = TFile::Open(fKsKl.c_str());
    ksTr = std::unique_ptr<TTree>(file->Get<TTree>("ksTree"));

    std::cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "KsKl filename: " << fKsKl << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;

    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("etrue", &etrue);
    etrue = -1;
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", &ksdpsi);
    ksTr->SetBranchAddress("Y", &Y);
    ksTr->SetBranchAddress("nhitPos", &nhitPos);
    ksTr->SetBranchAddress("nhitNeg", &nhitNeg);
    ksTr->SetBranchAddress("kstheta", &ksTheta);
    ksTr->SetBranchAddress("ksphi", &ksPhi);

    ksTr->SetBranchAddress("piThetaPos", &piThetaPos);
    ksTr->SetBranchAddress("piThetaNeg", &piThetaNeg);

    
    hMlnY = std::make_unique<TH2D>(TH2D("hMlnY", "M(lnY)", 300, -0.3, 0.3, 600, 480, 520));
    hDeltaM = std::make_unique<TProfile>(TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1));
    hMlnYpfx  = std::make_unique<TProfile>(TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505));
    hMPsi = std::make_unique<TH2D>(TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520));
    hM_CrAnglelnY = std::make_unique<TH2D>(TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 30, -0.4, 0.4, 40000, 490, 515));
    hPsilnY = std::make_unique<TH2D>(TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 10000, 2.4, 3.3));
    hPsi = std::make_unique<TH1D>(TH1D("hPsi", "Psi distr", 1000, 0, 3.15));
    hEnergySpectrum = std::make_unique<TH1D>(TH1D("hEnergySpectrum", "hEnergySpectrum", 6000, 480, 540));
    // After profile cut
    hEnergySpectrumCut = std::make_unique<TH1D>(TH1D("hEnergySpectrumCut", "hEnergySpectrumCut", 6000, 480, 540));

    hMassVsKsTheta = std::make_unique<TH2D>(TH2D("hMassVsKsTheta", "M vs KsTheta", 315, 0, 3.15, 600, 480, 520));
    hKsThetaVsLnY = std::make_unique<TH2D>(TH2D("hKsThetaVsLnY", "KsTheta vs lnY", 300, -0.6, 0.6, 315, 0, 3.15));

    hMlnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnY->GetXaxis()->SetTitle("ln(Y)");
    hMlnYpfx->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnYpfx->GetXaxis()->SetTitle("ln(Y)");

    hM_CrAnglelnY->GetXaxis()->SetTitle("lnY");
    hM_CrAnglelnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");

    ReadBadRuns(badRunsFileName);
    BadRunSearch();
    CalcSigmas(isVerbose);
}

MassHandler::~MassHandler()
{
    // delete ksTr;
    // delete massCrAngle;
    // delete massFullRec;
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
    input.close();
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

Float_t MassHandler::GetCorrectedEnergy(int rnum, double deltaE)
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

    return energy + deltaE;
}

void MassHandler::CalcSigmas(bool verbose)
{
    std::vector<std::unique_ptr<TH1D>> psilnYprojYs;
    std::vector<std::unique_ptr<TH1D>> Ys;
    for(int i = 0; i < 16; i++)
    { 
        psilnYprojYs.push_back(std::make_unique<TH1D>(TH1D(("py" + std::to_string(i + 1)).c_str(), ("Y" + std::to_string(i + 1)).c_str(), 10000, 2.4, 3.3)));
        Ys.push_back(std::make_unique<TH1D>(TH1D(("hY" + std::to_string(i + 1)).c_str(), ("Y" + std::to_string(i + 1)).c_str(), 400, 0, 4))); 
    }

    double massFullRecWithEmeas = 0;
    double lnY = 0;
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);

        if(auto bin = int((ksTheta + 1e-12) / 0.2); bin < Ys.size())
        { Ys[bin]->Fill(Y); }
        
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9 && nhitPos > 10 && nhitNeg > 10 &&
            1.1 < piThetaPos && piThetaPos < TMath::Pi() - 1.1 && 
            1.1 < piThetaNeg && piThetaNeg < TMath::Pi() - 1.1 )
        {
            lnY = log(Y);
            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(ksdpsi);
    
            if(auto psiBinNum = lnY > -0.4 ? int(floor((lnY + 0.4 + 1e-12) / 0.05)) : -1; 
                psiBinNum >= 0 && psiBinNum < psilnYprojYs.size() && 
                massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
            { psilnYprojYs[psiBinNum]->Fill(ksdpsi); }
        }
    }

    TFitResultPtr r;
    for(int i = 0; i < psilnYprojYs.size(); i++)
    {
        psilnYprojYs[i]->Rebin(40);
        r = psilnYprojYs[i]->Fit("gaus", "SEQ", "", psilnYprojYs[i]->GetXaxis()->GetBinCenter(psilnYprojYs[i]->GetMaximumBin()) - 3 * psilnYprojYs[i]->GetStdDev(), 
                                                psilnYprojYs[i]->GetXaxis()->GetBinCenter(psilnYprojYs[i]->GetMaximumBin()) + 3 * psilnYprojYs[i]->GetStdDev());
        r = psilnYprojYs[i]->Fit("gaus", "SEQ", "", r->Parameter(1) - 2 * r->Parameter(2), r->Parameter(1) + 2 * r->Parameter(2));
        vSigmaFit.push_back(r->Parameter(2));
        vSigmaRMS.push_back(psilnYprojYs[i]->GetStdDev());

        vY_RMS.push_back(Ys[i]->GetStdDev());
    }
    Ys[1]->DrawClone();
    auto printSigmas = [](std::string name, const std::vector<Float_t> &sigmas) {
        std::cout << "\n" << name << ": " << std::endl;
        for(const auto &sigma : sigmas)
        { std::cout << sigma << ", "; }
        std::cout << "\n" << std::endl;
    };

    if(verbose)
    { 
        printSigmas("Sigma", vSigmaFit);
        printSigmas("RMS", vSigmaRMS);
        printSigmas("Y RMS", vY_RMS);
        std::cout << "\n" <<std::endl; 
    }
}

void MassHandler::FillHists(double fitRange, double deltaE, bool withFitSigma, bool useEtrue)
{
    auto sigmas = std::vector<double>({0.0303571, 0.0244383, 0.0230856, 0.0203081, 0.0182954, 0.0188669, 0.0167881, 0.0151023, 0.0154767, 0.0172125, 0.0169003, 0.0184711, 0.0220203, 0.0231788, 0.0256208, 0.0287016});
    double sigmaPsi = 0;
    double sigmaY = 0;
    double energy = 0;
    double massFullRecWithEmeas = 0;
    double lnY = 0;
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9 && nhitPos > 10 && nhitNeg > 10 &&
            1.1 < piThetaPos && piThetaPos < TMath::Pi() - 1.1 && 
            1.1 < piThetaNeg && piThetaNeg < TMath::Pi() - 1.1 )
        {
            lnY = log(Y);
            if(auto psiBinNum = lnY > -0.4 ? int(floor((lnY + 0.4 + 1e-12) / 0.05)) : -1; psiBinNum >= 0 && psiBinNum < vSigmaFit.size())
            { sigmaPsi = withFitSigma? vSigmaFit[psiBinNum] : vSigmaRMS[psiBinNum]; }
            // { sigmaPsi = sigmas[psiBinNum]; }

            if(auto bin = int((ksTheta + 1e-12) / 0.2); bin < vY_RMS.size())
            { sigmaY = vY_RMS[bin]; }

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi);

            emeas = useEtrue? etrue : GetCorrectedEnergy(runnum, deltaE);

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massCrAngle->SetParameter(0, emeas);
        
            auto massFullRecVal = massFullRec->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi);
            hMlnY->Fill(lnY, massFullRecVal);
            hM_CrAnglelnY->Fill(lnY, massCrAngle->Eval(ksdpsi) - sigmaPsi * sigmaPsi / 2 * massCrAngle->Derivative2(ksdpsi));
            if(massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
            { 
                hDeltaM->Fill(lnY, - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(ksdpsi));
                hMlnYpfx->Fill(lnY, massFullRecVal); 
                hPsilnY->Fill(lnY, ksdpsi); 

                hMassVsKsTheta->Fill(ksTheta, massFullRecVal);
                hKsThetaVsLnY->Fill(lnY, ksTheta);
                if(fabs(lnY) < fitRange + 1e-6)
                { hEnergySpectrumCut->Fill(etrue); }
            }

            hMPsi->Fill(ksdpsi, massFullRecVal);
            hPsi->Fill(ksdpsi);
            hEnergySpectrum->Fill(etrue);
        }
    }

    grMassLnYFit = FitSlices(hMlnY, 10);
    grPsiLnYFit = FitSlices(hPsilnY, 16, std::make_pair(2, 1), 8);

    std::cout << std::endl << "Average E cut = " << hEnergySpectrumCut->GetMean() << " +/- " << hEnergySpectrumCut->GetMeanError() << std::endl;
    std::cout << "Average E  = " << hEnergySpectrum->GetMean() << " +/- " << hEnergySpectrum->GetMeanError() << "\n\n\n" << std::endl;
}

void MassHandler::DrawGraphs(int drawOpt)
{
    auto massLine = std::make_unique<TLine>(TLine(-0.3, 497.614, 0.3, 497.614));
    massLine->SetLineColor(kBlue);
    massLine->SetLineWidth(2);

    auto MCrAngleProf = std::unique_ptr<TProfile>(hM_CrAnglelnY->ProfileX());
    MCrAngleProf->GetXaxis()->SetTitle("lnY");
    MCrAngleProf->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    
    switch (drawOpt)
    {
    case 0:
        hMlnYpfx->GetXaxis()->SetRangeUser(-0.5, 0.5);
        hMlnYpfx->DrawClone();
        massLine->DrawClone("same");
        // hMassVsKsTheta->DrawClone("col");
        // hKsThetaVsLnY->DrawClone("col");
        break;
    case 1:
        hM_CrAnglelnY->DrawClone();
        break;
    case 2:
        // hPsilnY->ProfileX()->DrawClone();
        grPsiLnYFit->DrawClone();
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
        // grMlnYFit->DrawClone();
    default:
        break;
    }
}

std::unique_ptr<TGraphErrors> MassHandler::FitSlices(const std::unique_ptr<TH2D> &hist, int nSlices, std::pair<double, double> fitRange, std::optional<int> rebinFactor, bool isVerbose)
{
    auto grMassLnY_Fit = std::make_unique<TGraphErrors>(TGraphErrors());
    TFitResultPtr res;

    int nOldBins = hist->GetNbinsX();
    int firstBin = 1;
    int lastBin = int(nOldBins / nSlices);

    for (size_t i = 0; i < nSlices; i++)
    {
        std::string title = hist->GetName() + std::string("_py") + std::to_string(i);
        auto sliceProj = std::unique_ptr<TH1D>(hist->ProjectionY(title.c_str(), firstBin, lastBin));
        
        if(rebinFactor)
        { sliceProj->Rebin(rebinFactor.value()); }

        res = sliceProj->Fit("gaus", "SQME", "", sliceProj->GetBinCenter(sliceProj->GetMaximumBin()) - fitRange.first * sliceProj->GetRMS(), 
                                                    sliceProj->GetBinCenter(sliceProj->GetMaximumBin()) + fitRange.second * sliceProj->GetRMS());
        res = sliceProj->Fit("gaus", "SQME", "", res->Parameter(1) - fitRange.first * res->Parameter(2),  res->Parameter(1) + fitRange.second * res->Parameter(2));
        
        if(isVerbose)
        { 
            std::cout << title << " fit; (firstBin, lastBin) = (" << hist->GetXaxis()->GetBinCenter(firstBin) << ", " << hist->GetXaxis()->GetBinCenter(lastBin) << "); "<<
                        "chi2/ndf = " << res->Chi2() << "/" << res->Ndf() << "; Prob = " << res->Prob() <<  "; Mean = " << res->Parameter(1) << "; MeanErr = " << 
                        res->ParError(1) << "; Sigma = " << res->Parameter(2) << std::endl; 
        }

        grMassLnY_Fit->SetPoint(i, hist->GetXaxis()->GetBinCenter(int((lastBin + firstBin) / 2)), res->Parameter(1));
        grMassLnY_Fit->SetPointError(i, 0., res->ParError(1));
        
        firstBin = lastBin;
        lastBin = lastBin + int(nOldBins / nSlices);
    }

    if(isVerbose)
    { std::cout << std::endl; }

    return grMassLnY_Fit;
}

std::pair<double, double> MassHandler::GetMass(double fitRange, bool useEtrue, double deltaE, bool withFitSigma, int drawOpt, int methodOpt, bool verbose)
{   
    FillHists(fitRange, deltaE, withFitSigma, useEtrue);

    TFitResultPtr r, res0, res1, res2, res3, res4;
    double mass = 0;
    double massErr = 0;

    massCrAngle->SetParameter(0, useEtrue? hEnergySpectrumCut->GetMean() : GetCorrectedEnergy(runnum, deltaE));

    auto printFitRes = [&verbose, &methodOpt](const TFitResultPtr &res, std::string resName, int methodNum) {
        if(verbose || methodNum == methodOpt) {
            std::cout << "\n" << resName << " = " << res->Parameter(0) << " +/- " << res->ParError(0) 
            << "; chi2 / ndf = " << res->Chi2() << "/" << res->Ndf() << "; Prob = " << res->Prob() << std::endl;              
        }
    };

    res0 = hMlnYpfx->Fit("pol0", "SMQE", "goff", -fitRange, fitRange);
    printFitRes(res0, "Mass_FullRec", 0);

    // res1 = grMassLnYFit->Fit("pol0", "SMQE", "goff", -fitRange, fitRange);
    // printFitRes(res1, "MassVsLnY_Fit_FullRec", 1);     

    res2 = hM_CrAnglelnY->ProfileX()->Fit("pol4", "SMQE", "", -0.15, 0.15);
    printFitRes(res2, "Mass_CrAngle", 2);  

    res3 = hPsilnY->ProfileX("pfxAng")->Fit("pol2", "SQME", "goff", -0.15, 0.15);
    printFitRes(res3, "Psi", 3);  
    std::cout << "M_crAngle(Psi_cr pfx) = " << massCrAngle->Eval(res3->Parameter(0)) << " +/- " << massCrAngle->Derivative(res3->Parameter(0)) * res3->ParError(0) << "\n" << std::endl;

    res4 = grPsiLnYFit->Fit("pol2", "SQME", "goff", -0.15, 0.15);
    printFitRes(res4, "PsiFit", 4);
    std::cout << "M_crAngle(Psi_cr fit) = " << massCrAngle->Eval(res4->Parameter(0)) << " +/- " << massCrAngle->Derivative(res4->Parameter(0)) * res4->ParError(0) << std::endl;
    std::cout << std::endl;

    switch (methodOpt)
    {
    case 1:
        mass = res1->Parameter(0);
        massErr = res1->ParError(0);
        break;
        
    case 2:
        mass = res2->Parameter(0);
        massErr = res2->ParError(0);
        break;

    case 3:
        mass = massCrAngle->Eval(res3->Parameter(0));
        massErr = massCrAngle->Derivative(res3->Parameter(0)) * res3->ParError(0);
        break;

    case 4:
        mass = massCrAngle->Eval(res4->Parameter(0));
        massErr = massCrAngle->Derivative(res4->Parameter(0)) * res4->ParError(0);
        break;

    default:
        mass = res0->Parameter(0);
        massErr = res0->ParError(0);
        break;
    }

    DrawGraphs(drawOpt);
    return std::pair(mass, massErr);
}

int Scan(const std::map<std::string, std::pair<double, double>> &meanEnergies, double deltaE = 0, bool isExp = true)
{
    double lnYrange = 0.33;
    std::vector<double> vals = {};
    std::vector<double> errs = {};

    std::string fadRunsFile;
    std::string energyFilename;
    std::string filename;

    for(const auto&[energyPoint, meanEnergy] : meanEnergies)
    {
        if(energyPoint == "505" || energyPoint == "508") 
        { lnYrange = 0.2; }
        else 
        { lnYrange = 0.27; }

        if(isExp)
        {
            fadRunsFile = "txt/BadRuns.txt";
            filename = "tr_ph/expKsKl/exp" + energyPoint + ".root";
            energyFilename = "C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoint + ".root";
        }
        else
        {
            fadRunsFile = "";
            energyFilename = "C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoint + ".root";
            filename = "tr_ph/MC/KsKl_Smeared/MC" + energyPoint + ".root";
        }

        auto massHandler = new MassHandler(filename, fadRunsFile, meanEnergy.first);

        std::cout << "\n\n" << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << 
                    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << 
                    "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
        std::cout << "Energy Point = " << energyPoint << std::endl;
        auto [mass, massErr] = massHandler->GetMass(0.15, false, deltaE, true, 0, 4, false);

        vals.push_back(mass);
        errs.push_back(massErr);

        delete massHandler;
    }

    std::cout << "\n\n\nMasses = {";
    for(auto val : vals)
    { std::cout << val << ", "; }
    std::cout << "}" << std::endl;

    std::cout << "Errs = {";
    for(auto err : errs)
    { std::cout << err << ", "; }
    std::cout << "}" << std::endl;
    return 0;
}

int KsMassMeas()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    const std::vector<double> meanEnergies_vec = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    // Means of Ks energy spectrums
    const std::vector<double> meanEnergiesSpectrum_vec = {504.683, 507.762, 508.323, 508.885, 509.445, 509.841, 510.263, 510.694, 512.32};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};

    // RC without energy shift.
    // const std::vector<double> deltaM_RC_Smeared = {0.112, 0.081, 0.077, 0.068, 0.076, 0.108, 0.185, 0.324, 1.462};

    // RC with energy shifted by 134 keV
    const std::vector<double> deltaM_RC_Smeared = {0.1021, 0.0966, 0.087, 0.0759, 0.0808, 0.1336, 0.2191, 0.3726, 1.5223};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    std::map<std::string, std::pair<double, double>> meanEnergies;
    std::map<std::string, std::pair<double, double>> meanEnergiesSpectrum;
    std::map<std::string, double> radiativeCorrections;
    for(int i = 0; i < energyPoints.size(); i++)
    { 
        meanEnergies[energyPoints[i]] = std::make_pair(meanEnergies_vec[i], meanEnergiesErr[i]); 
        meanEnergiesSpectrum[energyPoints[i]] = std::make_pair(meanEnergiesSpectrum_vec[i], meanEnergiesErr[i]); 
        radiativeCorrections[energyPoints[i]] = deltaM_RC_Smeared[i]; 
    }

    std::string energyPoint = "508.5";
    double deltaE = 0.0;
    std::string fileName = "tr_ph/expKsKl/exp" + energyPoint + ".root";
    fileName = "tr_ph/MC/KsKl_Smeared/MC" + energyPoint + ".root";
    // fileName = "tr_ph/MC/KsKl_Smeared/Field/MC" + energyPoint + "_Field.root";

    std::string fadRunsFile = "txt/BadRuns.txt";

    // Scan(meanEnergies, deltaE, true);
    // Scan(meanEnergiesSpectrum, deltaE, false);
    
    auto massHandler = new MassHandler(fileName, fadRunsFile, meanEnergies[energyPoint].first);
    auto mass = massHandler->GetMass(0.27, false, 0., true, 0, 0).first;
    std::cout << "M_NCRC_smeared = " << mass << "\n\n" << std::endl ;
    delete massHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}