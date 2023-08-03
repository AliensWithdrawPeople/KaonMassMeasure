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
#include "TVector3.h"
#include "TSpline.h"
#include "TRandom.h"

#include "Energy.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <optional>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <chrono>
#include <ctime> 

double deltaPsi(double phi1, double theta1, double phi2, double theta2) 
{ return acos(sin(theta1) * sin(theta2) * cos(phi1 - phi2) + cos(theta1) * cos(theta2)); }

double dPsiCorrection(double ksTheta, std::string energy)
{
    auto theta = (ksTheta - TMath::Pi() / 2);
    if(energy == "501") { return -0.000672261 + 0.0160195 * std::pow(theta, 2); }

    if(energy == "503") { return -0.000831758 + 0.0126012 * std::pow(theta, 2); }

    if(energy == "505") { return -0.000967487 + 0.0101455 * std::pow(theta, 2); }

    if(energy == "508") { return -0.0011324 + 0.010266 * std::pow(theta, 2); }

    if(energy == "508.5") { return -0.00131057 + 0.0127361 * std::pow(theta, 2); }


    // if(energy == "509") { return -0.00134738 + 0.0140213 * std::pow(theta, 2); }
    if(energy == "509") { return 8.91930e-04 -1.15960e-02 * std::pow(theta, 2) + 6.31244e-03 * std::pow(theta, 4); }


    if(energy == "509.5") { return -0.00122879 + 0.0128762 * std::pow(theta, 2); }

    if(energy == "510") { return -0.00122787 + 0.0114095 * std::pow(theta, 2); }

    if(energy == "510.5") { return -0.0012652 + 0.0128049 * std::pow(theta, 2); }

    if(energy == "511") { return -0.0010293 + 0.00869987 * std::pow(theta, 2); }

    if(energy == "511.5") { return -0.00102064 + 0.010116 * std::pow(theta, 2); }

    if(energy == "514") { return -0.0012627 + 0.0105834 * std::pow(theta, 2); }

    if(energy == "517") { return -0.00101428 + 0.00586478 * std::pow(theta, 2); }

    if(energy == "520") { return -0.0011453 + 0.00921438 * std::pow(theta, 2); }

    if(energy == "525") { return -0.00149835 + 0.0104855 * std::pow(theta, 2); }

    if(energy == "530") { return -0.000983491 + 0.00858641 * std::pow(theta, 2); }

    return 0;
}

double Ycorrection(std::string energy)
{
    if(energy == "501") { return 0.000346437; }

    if(energy == "503") { return 0.00108152; }

    if(energy == "505") { return 0.00132092; }

    if(energy == "508") { return 0.00122104; }

    if(energy == "508.5") { return 0.00124719; }

    if(energy == "509") { return 0.0014927; }

    if(energy == "509.5") { return 0.00187755; }

    if(energy == "510") { return 0.00150114; }

    if(energy == "510.5") { return 0.00163287; }

    if(energy == "511") { return 0.00135194; }

    if(energy == "511.5") { return 0.00143905; }

    if(energy == "514") { return 0.00195023; }

    if(energy == "517") { return 0.00154664; }

    if(energy == "520") { return 0.00192852; }

    if(energy == "525") { return 0.00175594; }

    if(energy == "530") { return 0.00221627; }

    return 0.;
}

class MassHandler
{
private:
    std::unique_ptr<TTree> ksTr;

    std::unique_ptr<TF1> massCrAngle;
    std::unique_ptr<TF1> massFullRec;

    std::string energyPoint;

    Float_t emeas; Float_t demeas;
    // Kaon energy calculated as invariant mass of pi+pi-.
    Float_t etrue;
    Float_t ksTheta;
    Float_t ksPhi;

    Int_t nhitPos;
    Int_t nhitNeg;

    Int_t runnum; 

    // Momentum ratio = P1/P2, where P1 is the momentum of pi+, P2 is the momentum of pi-.
    Float_t Y;
    Float_t ksdpsi;
    Float_t piThetaPos;
    Float_t piThetaNeg;
    Float_t piPhiPos;
    Float_t piPhiNeg;

    Float_t Y_gen;
    Float_t ksdpsi_gen;
    Float_t piThetaPos_gen;
    Float_t piThetaNeg_gen;
    Float_t piPhiPos_gen;
    Float_t piPhiNeg_gen;

    // Sigma_psi(lnY, kstheta) from gaus fit.
    std::vector<std::vector<Float_t>> vSigmaMatrixFit;

    std::vector<Float_t> vSigmaFit;
    std::vector<Float_t> vSigmaRMS;
    std::vector<int> badRuns;

    std::unique_ptr<TH2D> hMlnY;
    std::unique_ptr<TH2D> hMPsi;
    std::unique_ptr<TH2D> hM_CrAnglelnY;
    std::unique_ptr<TH2D> hPsilnY;
    std::unique_ptr<TH2D> hMassVsKsTheta;
    std::unique_ptr<TH2D> hKsThetaVsLnY;

    std::unique_ptr<TProfile> hDeltaM;
    std::unique_ptr<TProfile> hMlnYpfx;

    TProfile* hMlnYpfx_cowboy;
    TProfile* hMlnYpfx_sailor;

    std::unique_ptr<TH1D> hEnergySpectrum;
    // After profile cut
    std::unique_ptr<TH1D> hEnergySpectrumCut;
    std::unique_ptr<TH1D> hPsi;

    std::unique_ptr<TGraphErrors> grMassLnYFit;
    std::unique_ptr<TGraphErrors> grPsiLnYFit;

    std::unique_ptr<TH2D> hThetaDiffRecGen;
    std::unique_ptr<TH2D> hDiffRecGen_cowboy;
    std::unique_ptr<TH2D> hDiffRecGen_sailor;


    std::optional<double> energyCorrected;
    std::optional<std::map<int, Float_t>> energyDiff;

    TSpline3* ksdpsi_RecGenDiff_Spline;


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

    /*  mode == 0 -- create spline, 
        mode == 1 -- create spline and save it to the file(spline_filename), 
        mode == 2 -- read spline from the file(spline_filename).
    */
    void CreateSplines(std::string spline_filename = "", int mode = 0);

    /* TO-DO!!!
        Returns angle between pions momentums psi with applied correction from  ksdpsi_RecGenDiff_Spline.
    */
    //double GetCorrectedPsi(double phiPos, double thetaPos, double phiNeg, double thetaNeg);
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

    /*
    *splineMode == 0 -- create spline, 
    *           == 1 -- create spline and save it to the file(spline_filename), 
    *           == 2 -- read spline from the file(spline_filename).
    */
    MassHandler(std::string fKsKl, std::string energyPoint, std::string badRunsFileName = "", 
                std::string splineFilename = "", int splineMode = 0,
                std::optional<double> energyCorr = std::nullopt, 
                std::optional<std::map<int, Float_t>> energyDiff = std::nullopt,
                bool isVerbose = true);
    ~MassHandler();
};

MassHandler::MassHandler(std::string fKsKl, std::string energyPoint, 
                        std::string badRunsFileName, std::string splineFilename, int splineMode,
                        std::optional<double> energyCorr, std::optional<std::map<int, Float_t>> energyDiff, bool isVerbose)
{
    this->energyPoint = energyPoint;
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
    ksTr->SetBranchAddress("piPhiPos", &piPhiPos);
    ksTr->SetBranchAddress("piPhiNeg", &piPhiNeg);

    ksTr->SetBranchAddress("ksdpsi_gen", &ksdpsi_gen);
    ksTr->SetBranchAddress("Y_gen", &Y_gen);
    ksTr->SetBranchAddress("piPhiPos_gen", &piPhiPos_gen);
    ksTr->SetBranchAddress("piPhiNeg_gen", &piPhiNeg_gen);
    ksTr->SetBranchAddress("piThetaPos_gen", &piThetaPos_gen);
    ksTr->SetBranchAddress("piThetaNeg_gen", &piThetaNeg_gen);
    
    hMlnY = std::make_unique<TH2D>(TH2D("hMlnY", "M(lnY)", 300, -0.3, 0.3, 600, 480, 520));
    hDeltaM = std::make_unique<TProfile>(TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1));
    hMlnYpfx  = std::make_unique<TProfile>(TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505));
    hMPsi = std::make_unique<TH2D>(TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520));
    hM_CrAnglelnY = std::make_unique<TH2D>(TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 30, -0.4, 0.4, 40000, 490, 515));
    hPsilnY = std::make_unique<TH2D>(TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 10000, 2.4, 3.3));
    hPsi = std::make_unique<TH1D>(TH1D("hPsi1", "Psi distr", 20000, -1, 6.3));
    hEnergySpectrum = std::make_unique<TH1D>(TH1D("hEnergySpectrum", "hEnergySpectrum", 6000, 480, 540));
    // After profile cut
    hEnergySpectrumCut = std::make_unique<TH1D>(TH1D("hEnergySpectrumCut", "hEnergySpectrumCut", 6000, 480, 540));

    hMassVsKsTheta = std::make_unique<TH2D>(TH2D("hMassVsKsTheta", "M vs KsTheta", 600, -1.57, 1.57, 40000, 480, 520));
    hKsThetaVsLnY = std::make_unique<TH2D>(TH2D("hKsThetaVsLnY", "KsTheta vs lnY", 300, -0.6, 0.6, 315, 0, 3.15));

    hThetaDiffRecGen = std::make_unique<TH2D>(TH2D("hThetaDiffRecGen", "hThetaDiffRecGen", 1200, -6., 6., 20000, -1, 1));
    hDiffRecGen_cowboy = std::make_unique<TH2D>(TH2D("hDiffRecGen_cowboy", "hDiffRecGen_cowboy", 1200, -6., 6., 20000, -1, 1));
    hDiffRecGen_sailor = std::make_unique<TH2D>(TH2D("hDiffRecGen_sailor", "hDiffRecGen_sailor", 1200, -6., 6., 20000, -1, 1));

    hMlnYpfx_cowboy = new TProfile("hMlnYpfx_cowboy","Profile of M versus lnY, cowboy", 30, -1, 1, 490, 505);
    hMlnYpfx_sailor = new TProfile("hMlnYpfx_sailor","Profile of M versus lnY, sailor", 30, -1, 1, 490, 505);

    hMlnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnY->GetXaxis()->SetTitle("ln(Y)");
    hMlnYpfx->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    hMlnYpfx->GetXaxis()->SetTitle("ln(Y)");

    hM_CrAnglelnY->GetXaxis()->SetTitle("lnY");
    hM_CrAnglelnY->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");

    ReadBadRuns(badRunsFileName);
    BadRunSearch();
    
    CreateSplines(splineFilename, splineMode);
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

void MassHandler::CreateSplines(std::string spline_filename, int mode)
{
    if(mode != 0 && mode != 1 && mode != 2)
    { std::cout << "Wrong mode in CreateSplines!!!" << std::endl; }

    if(mode == 0 || mode == 1)
    {
        std::vector<Double_t> psi_delta = {};
        std::vector<Double_t> ksTheta_bins = {};
        TProfile psi("psi", "deltaPsi", 50, 0, 3.14, -0.1, 0.1);

        for(int i = 0; i < ksTr->GetEntries(); i++)
        {
            ksTr->GetEntry(i);
            if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9 && 
                nhitPos > 10 && nhitNeg > 10 && fabs(ksTheta - TMath::Pi() / 2) < 0.5)
            { psi.Fill(ksTheta, ksdpsi - ksdpsi_gen); }
        }

        for(int i = 0; i < psi.GetNbinsX(); i++)
        {
            ksTheta_bins.push_back(psi.GetBinCenter(i));
            psi_delta.push_back(psi.GetBinContent(i));
        }

        ksdpsi_RecGenDiff_Spline = new TSpline3("ksdpsi_RecGenDiff_Spline", ksTheta_bins.data(), psi_delta.data(), psi_delta.size());
        ksdpsi_RecGenDiff_Spline->SetName("ksdpsi_RecGenDiff_Spline");
        if(mode == 1)
        { 
            auto splineFile = TFile::Open(spline_filename.c_str(), "recreate");
            ksdpsi_RecGenDiff_Spline->Write();
            splineFile->Close();
            delete splineFile;
        }
    }
    
    if(mode == 2)
    { ksdpsi_RecGenDiff_Spline = TFile::Open(spline_filename.c_str())->Get<TSpline3>("ksdpsi_RecGenDiff_Spline"); }
}

void MassHandler::CalcSigmas(bool verbose)
{
    std::vector<std::unique_ptr<TH1D>> psilnYprojYs;
    for(int i = 0; i < 16; i++)
    { 
        psilnYprojYs.push_back(std::make_unique<TH1D>(TH1D(("py" + std::to_string(i + 1)).c_str(), ("Y" + std::to_string(i + 1)).c_str(), 10000, 2.4, 3.3)));
    }

    // Distributions of ksdpsi (angle between pions' momentums) in different bins (lnY, kstheta).
    std::vector<std::vector<std::unique_ptr<TH1D>>> psiDistrs;
    for(int i = 0; i < 8; i++)
    {
        vSigmaMatrixFit.push_back({});
        psiDistrs.push_back(std::vector<std::unique_ptr<TH1D>>());
        for(int j = 0; j < 12; j++)
        {
            std::string histName = "py_bin_lnY" + std::to_string(i) + "_kstheta" + std::to_string(j); 
            psiDistrs[i].push_back(std::make_unique<TH1D>(TH1D(histName.c_str(), histName.c_str(), 2000, -3.1415, 3.1415))); 
        }
    }

    double massFullRecWithEmeas = 0;
    double lnY = 0;
    TVector3 piPos(1, 1, 1);
    TVector3 piNeg(1, 1, 1);

    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);

        piPos.SetMagThetaPhi(1, piThetaPos, piPhiPos);
        piNeg.SetMagThetaPhi(1, piThetaNeg, piPhiNeg);
        ksdpsi = piPos.Angle(piNeg);
        ksdpsi -= ksdpsi_RecGenDiff_Spline->Eval(ksTheta);

        lnY = log(Y);
        massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
        massFullRecWithEmeas = massFullRec->Eval(ksdpsi);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9 && nhitPos > 10 && nhitNeg > 10 &&
            1.1 < piThetaPos && piThetaPos < TMath::Pi() - 1.1 &&
            1.1 < piThetaNeg && piThetaNeg < TMath::Pi() - 1.1 &&
            massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
        {
            if(auto psiBinNum = lnY > -0.4 ? int(floor((lnY + 0.4 + 1e-12) / 0.05)) : -1; 
                psiBinNum >= 0 && psiBinNum < psilnYprojYs.size())
            { psilnYprojYs[psiBinNum]->Fill(ksdpsi); }
            
            auto binY = lnY > -0.3? int(floor((lnY + 0.3) / 0.075)) : -1;
            auto binKsTheta = ksTheta > 0.57? int(floor((ksTheta - 0.57) / 0.25)) : -1;
            if(binY != -1 && binKsTheta != -1 &&
               binY < psiDistrs.size() && binKsTheta < psiDistrs.back().size())
            { psiDistrs[binY][binKsTheta]->Fill(ksdpsi); } 
        }
    }

    if(energyPoint == "514")
    {
        vSigmaFit = std::vector<Float_t>({0.0303571, 0.0244383, 0.0230856, 0.0203081, 0.0182954, 0.0188669, 0.0167881, 0.0151023, 0.0154767, 0.0172125, 0.0169003, 0.0184711, 0.0220203, 0.0231788, 0.0256208, 0.0287016});
        return;
    }

    TFitResultPtr r;
    for(int i = 0; i < psilnYprojYs.size(); i++)
    {
        psilnYprojYs[i]->Rebin(40);

        auto leftBorder = psilnYprojYs[i]->GetXaxis()->GetBinCenter(psilnYprojYs[i]->GetMaximumBin()) - 3 * psilnYprojYs[i]->GetStdDev();
        auto rightBorder = psilnYprojYs[i]->GetXaxis()->GetBinCenter(psilnYprojYs[i]->GetMaximumBin()) + 3 * psilnYprojYs[i]->GetStdDev();
        r = psilnYprojYs[i]->Fit("gaus", "SEQ", "", leftBorder, rightBorder);
        r = psilnYprojYs[i]->Fit("gaus", "SEQ", "", r->Parameter(1) - 2 * r->Parameter(2), r->Parameter(1) + 2 * r->Parameter(2));
        vSigmaFit.push_back(r->Parameter(2));
        vSigmaRMS.push_back(psilnYprojYs[i]->GetStdDev());
    }

    for(int i = 0; i < 8; i++)
    {
        for(int j = 0; j < 8; j++)
        {
            if(psiDistrs[i][j]->GetEntries() < 150)
            { 
                // vSigmaMatrixFit[i].push_back(psiDistrs[i][j]->GetEntries() > 50 ? psiDistrs[i][j]->GetStdDev() : 0.0); 
                vSigmaMatrixFit[i].push_back(psiDistrs[i][j]->GetStdDev()); 
                continue;
            }

            if(psiDistrs[i][j]->GetEntries() < 600)
            { psiDistrs[i][j]->Rebin(2 * int(3 - psiDistrs[i][j]->GetEntries() / 300)); }

            auto leftBorder = psiDistrs[i][j]->GetBinCenter(psiDistrs[i][j]->GetMaximumBin()) - 2 * psiDistrs[i][j]->GetStdDev();
            auto rightBorder = psilnYprojYs[i]->GetBinCenter(psiDistrs[i][j]->GetMaximumBin()) + 2 * psiDistrs[i][j]->GetStdDev();
            r = psiDistrs[i][j]->Fit("gaus", "SEQ", "", leftBorder, rightBorder);
            r = psiDistrs[i][j]->Fit("gaus", "SEQ", "", r->Parameter(1) - 2 * r->Parameter(2), r->Parameter(1) + 2 * r->Parameter(2));
            vSigmaMatrixFit[i].push_back(r->Parameter(2));
        }
    }

    auto printSigmas = [](std::string name, const std::vector<Float_t> &sigmas) {
        std::cout << "\n" << name << ": " << std::endl;
        for(const auto &sigma : sigmas)
        { std::cout << sigma << ", "; }
        std::cout << "\n" << std::endl;
    };

    if(verbose)
    { 
        auto default_precision = std::cout.precision();
        std::cout << "vSigmaMatrixFit:" << std::endl;
        for(int i = 0; i < 8; i++)
        {
            for(int j = 0; j < 8; j++)
            { std::cout << std::left << std::setprecision(4) << std::setw(8) << vSigmaMatrixFit[i][j] << " "; }
            std::cout << std::endl;
        }
        std::cout << std::setprecision(default_precision) << std::endl;

        printSigmas("Sigma", vSigmaFit);
        printSigmas("RMS", vSigmaRMS);
        std::cout << "\n" <<std::endl; 
    }
}

void MassHandler::FillHists(double fitRange, double deltaE, bool withFitSigma, bool useEtrue)
{
    auto sigmas = std::vector<double>({0.0303571, 0.0244383, 0.0230856, 0.0203081, 0.0182954, 0.0188669, 0.0167881, 0.0151023, 0.0154767, 0.0172125, 0.0169003, 0.0184711, 0.0220203, 0.0231788, 0.0256208, 0.0287016});
    double sigmaPsi = 0;
    double energy = 0;
    double massFullRecWithEmeas = 0;
    double lnY = 0;

    TVector3 piPos(1, 0, 0);
    TVector3 piNeg(1, 0, 0);
    TVector3 field(0., 0., 1.);

    hPsi->Reset();
    for(int i = 0; i < ksTr->GetEntries(); i++)
    {
        ksTr->GetEntry(i);
        if(std::find(badRuns.begin(), badRuns.end(), runnum) == badRuns.end() && abs(Y - 1) > 1e-9 && nhitPos > 10 && nhitNeg > 10 &&
            1.1 < piThetaPos && piThetaPos < TMath::Pi() - 1.1 &&
            1.1 < piThetaNeg && piThetaNeg < TMath::Pi() - 1.1 
            // && fabs(ksTheta - TMath::Pi() / 2) > 0.3
            && fabs(ksTheta - TMath::Pi() / 2) < 0.7
            )
        {
            lnY = log(Y);

            auto binY = lnY > -0.3? int(floor((lnY + 0.3) / 0.075)) : -1;
            auto binKsTheta = ksTheta > 0.57? int(floor((ksTheta - 0.57) / 0.25)) : -1;
            if(binY != -1 && binKsTheta != -1 &&
               binY < vSigmaMatrixFit.size() && binKsTheta < vSigmaMatrixFit.back().size())
            { sigmaPsi = vSigmaMatrixFit[binY][binKsTheta]; } 

            piPos.SetMagThetaPhi(1, piThetaPos, piPhiPos);
            piNeg.SetMagThetaPhi(1, piThetaNeg, piPhiNeg);
            ksdpsi = piPos.Angle(piNeg);
            auto theta = ksTheta - TMath::Pi() / 2;

            auto dpsi = ksdpsi - (ksdpsi_RecGenDiff_Spline->Eval(ksTheta) + 3e-4);

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massFullRecWithEmeas = massFullRec->Eval(dpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(dpsi);

            emeas = useEtrue? etrue : GetCorrectedEnergy(runnum, deltaE);

            massFullRec->SetParameters(emeas, (1 - Y*Y) / (1 + Y*Y));
            massCrAngle->SetParameter(0, emeas);

            auto massFullRecVal = massFullRec->Eval(dpsi) - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(dpsi);
            hM_CrAnglelnY->Fill(lnY, massCrAngle->Eval(dpsi) - sigmaPsi * sigmaPsi / 2 * massCrAngle->Derivative2(dpsi));

            hMlnY->Fill(lnY, massFullRecVal);
            if(massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
            { 
                hDeltaM->Fill(lnY, - sigmaPsi * sigmaPsi / 2 * massFullRec->Derivative2(dpsi));
                hMlnYpfx->Fill(lnY, massFullRecVal);

                if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) < TMath::Pi() / 2)
                { 
                    hDiffRecGen_cowboy->Fill(ksdpsi, ksdpsi - ksdpsi_gen);
                    hMlnYpfx_cowboy->Fill(lnY, massFullRecVal); 
                }

                if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) > TMath::Pi() / 2)
                { 
                    hDiffRecGen_sailor->Fill(ksdpsi, ksdpsi - ksdpsi_gen);
                    hMlnYpfx_sailor->Fill(lnY, massFullRecVal); 
                }

                hPsilnY->Fill(lnY, dpsi); 

                if(fabs(lnY) < 0.27)
                { hMassVsKsTheta->Fill(ksTheta - TMath::Pi() / 2, massFullRecVal); }
                hKsThetaVsLnY->Fill(lnY, ksTheta);
                if(fabs(lnY) < fitRange + 1e-6)
                { hEnergySpectrumCut->Fill(etrue); }
                hPsi->Fill(dpsi);
            }

            hMPsi->Fill(ksdpsi, massFullRecVal);
            hEnergySpectrum->Fill(etrue);
        }
    }
    // grMassLnYFit = FitSlices(hMlnY, 10);
    // grPsiLnYFit = FitSlices(hPsilnY, 16, std::make_pair(2, 1), 8);

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
        hMassVsKsTheta->DrawClone();

        hMlnYpfx_cowboy->GetXaxis()->SetRangeUser(-0.5, 0.5);
        hMlnYpfx_sailor->GetXaxis()->SetRangeUser(-0.5, 0.5);
        hMlnYpfx_cowboy->DrawClone();
        hMlnYpfx_sailor->DrawClone();
        hDiffRecGen_cowboy->DrawClone("same");
        hDiffRecGen_sailor->DrawClone("same");
        break;
    case 1:
        hM_CrAnglelnY->DrawClone();
        break;
    case 2:
        hPsilnY->ProfileX()->DrawClone();
        // grPsiLnYFit->DrawClone();
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

    // res4 = grPsiLnYFit->Fit("pol2", "SQME", "goff", -0.15, 0.15);
    // printFitRes(res4, "PsiFit", 4);
    // std::cout << "M_crAngle(Psi_cr fit) = " << massCrAngle->Eval(res4->Parameter(0)) << " +/- " << massCrAngle->Derivative(res4->Parameter(0)) * res4->ParError(0) << std::endl;
    // std::cout << std::endl;

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

int Scan(const std::map<std::string, std::pair<double, double>> &meanEnergies, double deltaE = 0, bool isExp = true, bool trueEnergy = false)
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

        std::string splineFilename = "C:/work/Science/BINP/Kaon Mass Measure/ksdpsi_splines/spline_" + energyPoint + ".root";
        int splineMode;


        if(isExp)
        {
            fadRunsFile = "txt/BadRuns.txt";
            filename = "tr_ph/expKsKl/exp" + energyPoint + ".root";
            energyFilename = "C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoint + ".root";
            splineMode = 2;
        }
        else
        {
            splineMode = 0;
            fadRunsFile = "";
            energyFilename = "C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoint + ".root";
            filename = "tr_ph/MC/KsKl_Smeared/Field/MC" + energyPoint + "_Field.root";
        }

        auto massHandler = new MassHandler(filename, energyPoint, fadRunsFile, splineFilename, splineMode, meanEnergy.first);

        std::cout << "\n\n" << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << 
                    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << 
                    "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
        std::cout << "Energy Point = " << energyPoint << std::endl;
        auto [mass, massErr] = massHandler->GetMass(lnYrange, trueEnergy, deltaE, (meanEnergy.first < 512)? true : false, 0, 0, false);

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
    // const std::vector<double> meanEnergiesSpectrum_vec = {504.683, 507.762, 508.323, 508.885, 509.445, 509.841, 510.263, 510.694, 512.32};
    const std::vector<double> meanEnergiesSpectrum_vec = {504.683, 507.762, 508.323, 508.885, 509.445, 509.841, 510.263, 510.694, 512.297};
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

    std::string energyPoint = "509";
    double deltaE = 0.0;
    std::string fileName = "tr_ph/expKsKl/exp" + energyPoint + ".root";
    // fileName = "tr_ph/MC/KsKl_Smeared/MC" + energyPoint + ".root";
    // fileName = "tr_ph/MC/KsKl_Smeared/Field/MC" + energyPoint + "_Field.root";

    std::string fadRunsFile = "txt/BadRuns.txt";
    std::string splineFilename = "C:/work/Science/BINP/Kaon Mass Measure/ksdpsi_splines/spline_" + energyPoint + ".root";
    // TO-DO!!! Create enum witm splineMode.
    int splineMode = 2;


    // Scan(meanEnergies, deltaE, false, true);
    // Scan(meanEnergiesSpectrum, deltaE, false, false);
    
    auto massHandler = new MassHandler(fileName, energyPoint, fadRunsFile, splineFilename, splineMode, meanEnergies[energyPoint].first);
    auto mass = massHandler->GetMass(0.27, false, 0., true, 0, 0, true).first;
    std::cout << "M_NCRC_smeared = " << mass << "\n\n" << std::endl ;
    delete massHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}