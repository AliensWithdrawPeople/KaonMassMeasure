#ifndef HandlerMC_hpp
#define HandlerMC_hpp

#include "Evaluator.hpp"
#include "Tree.hpp"
#include "HistContainer.hpp"
#include "Spline.hpp"
#include "CalcSigmas.hpp"
#include "Funcs.hpp"
#include "Painter.hpp"

#include <memory>
#include <utility>
#include <iostream>

#include "TMath.h"
#include "TVector3.h"


class HandlerMC: Evaluator
{
private:
    std::string energyPoint;
    std::optional<double> meanEnergy;
    bool verbose;
    std::unique_ptr<Tree> tree;
    const Tree::Data* data; 
    HistContainer container = HistContainer();
    std::unique_ptr<Painter> painter;

    double fitRange;
    bool useTrueEnergy = false;
    bool useEnergySmearing = true;

    /// @brief Numbers of entries that satisfy all selection cuts.
    std::vector<int> goodEntries = {};
    /// @brief Sigma_psi(lnY, kstheta) from gaus fit.
    std::vector<std::vector<double>> vSigmaMatrixFit;
    
    /// @brief Correction to psi_reco (dpsi = psi_reco - psi_gen).
    Spline deltaPsi_RecGenDiff;
    /// @brief Resolution of theta as a function of track's theta.
    Spline sigmaTheta;

    /// @brief Creates deltaPsi_RecGenDiff -- correction to psi_reco (dpsi = psi_reco - psi_gen).
    /// @param tree reference to Tree object storing data about events,
    Spline CreateDeltaPsiSpline();

    void FillHists(bool useTrueEnergy);

public:
    HandlerMC(std::string fKsKl, std::string energyPoint, double fitRange, std::optional<double> meanEnergy, 
                bool useTrueEnergy, bool saveSplines = false, bool useEnergySmearing = true, bool isVerbose = true);

    std::pair<double, double> GetMass(double fitRange = 0.27);
    std::pair<double, double> GetEnergySpectrumMean();
    
    std::pair<double, double> Eval() override { return GetMass(); }
    void Draw(std::string name) override;
    void Draw(std::string name, Range range);
    void SaveHists(std::string output_filename);
    void SaveSplines(std::string spline_filename);
};

HandlerMC::HandlerMC(std::string fKsKl, std::string energyPoint, double fitRange, std::optional<double> meanEnergy, 
                        bool useTrueEnergy, bool saveSplines, bool useEnergySmearing, bool isVerbose): 
    energyPoint{energyPoint}, fitRange{fitRange}, meanEnergy{meanEnergy}, useTrueEnergy{useTrueEnergy}, useEnergySmearing{useEnergySmearing}, verbose{isVerbose} 
{
    tree = std::unique_ptr<Tree>(new Tree(fKsKl));
    painter = std::unique_ptr<Painter>(new Painter(&container));
    data = tree->data; 

    container.Add("hMlnY", new TH2D("hMlnY", "M(lnY)", 800, -0.8, 0.8, 600, 485, 515));
    container.Add("hDeltaM", new TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1));
    container.Add("hPsilnY", new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 10000, 2.4, 3.3));
    container.Add("hEnergySpectrum", new TH1D("hEnergySpectrum", "hEnergySpectrum", 6000, 480, 540));
    container.Add("hMassVsKsTheta", new TH2D("hMassVsKsTheta", "M vs KsTheta", 600, -1.57, 1.57, 40000, 480, 520));

    container.Add("hDiffRecGen", new TH2D("hDiffRecGen", "hDiffRecGen", 1200, -6., 6., 20000, -1, 1));
    container.Add("hDiffRecGen_cowboy", new TH2D("hDiffRecGen_cowboy", "hDiffRecGen_cowboy", 1200, -6., 6., 20000, -1, 1));
    container.Add("hDiffRecGen_sailor", new TH2D("hDiffRecGen_sailor", "hDiffRecGen_sailor", 1200, -6., 6., 20000, -1, 1));

    container.Add("hMlnYpfx", new TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505));
    container.Add("hMlnYpfx_cowboy", new TProfile("hMlnYpfx_cowboy", "Profile of M versus lnY, cowboy", 30, -1, 1, 490, 505));
    container.Add("hMlnYpfx_sailor", new TProfile("hMlnYpfx_sailor", "Profile of M versus lnY, sailor", 30, -1, 1, 490, 505));

    for(const auto &entry : *tree.get())
    {
        if(abs(data->Y - 1) > 3e-7 && fabs(data->ks.theta - TMath::Pi() / 2) < 1 &&
            data->piPos.nhit > 10 && data->piNeg.nhit > 10 && 
            1.1 < data->piPos.theta && data->piPos.theta < TMath::Pi() - 1.1 &&
            1.1 < data->piNeg.theta && data->piNeg.theta < TMath::Pi() - 1.1
        )
        {
            auto massFullRecWithEmeas = FullRecMassFunc::Eval(data->ksdpsi, tree->emeas, data->Y);    
            if(massFullRecWithEmeas > 490 && massFullRecWithEmeas < 505)
            { goodEntries.push_back(entry); }
        }
    }

    vSigmaMatrixFit = Sigmas::GetSigmaMatrix(tree, goodEntries, verbose);
    auto [ksThetaBinCenters, deltaTheta, deltaThetaError] = Sigmas::GetThetaSigmas(tree, goodEntries);
    sigmaTheta = Spline(ksThetaBinCenters, deltaTheta, "sigmaTheta_spline");

    std::string spline_filename = misc::GetSplineFilename(energyPoint);
    deltaPsi_RecGenDiff = CreateDeltaPsiSpline();

    if(saveSplines)
    { SaveSplines(spline_filename); }

    FillHists(useTrueEnergy);
}

Spline HandlerMC::CreateDeltaPsiSpline()
{
    std::vector<Double_t> psi_delta = {};
    std::vector<Double_t> ksTheta_bins = {};
    TProfile psi("psi", "deltaPsi", 50, 0, 3.14, -0.1, 0.1);

    for(const auto &entry : goodEntries)
    {
        tree->GetEntry(entry);
        if(abs(tree->reco.Y - 1) > 1e-9 && tree->reco.piPos.nhit > 10 && tree->reco.piNeg.nhit > 10 && fabs(tree->reco.ks.theta - TMath::Pi() / 2) < 0.5)
        { psi.Fill(tree->reco.ks.theta, tree->reco.ksdpsi - tree->gen.ksdpsi); }
    }

    for(int i = 0; i < psi.GetNbinsX(); i++)
    {
        ksTheta_bins.push_back(psi.GetBinCenter(i));
        psi_delta.push_back(psi.GetBinContent(i));
    }

    return Spline(ksTheta_bins, psi_delta, "ksdpsiRecGenDiff_spline");
}

void HandlerMC::FillHists(bool useTrueEnergy)
{
    auto sigmaPhi = Sigmas::GetPhiSigma();
    auto sigmaY = Sigmas::GetYSigma();

    for(const auto &entry : goodEntries)
    {
        tree->GetEntry(entry);
        auto lnY = log(data->Y);
        container["hMlnY"]->Fill(lnY, FullRecMassFunc::Eval(tree->reco.ksdpsi, tree->emeas, data->Y));

        // Ks in a good region + optional cut to eliminate energy smearing.
        if(fabs(data->ks.theta - TMath::Pi() / 2) > 0.3 || (useEnergySmearing? false : fabs(tree->emeas - meanEnergy.value_or(0)) > 20e-3))
        { continue; }
        

        auto psiCor =   sigmaTheta(data->piPos.theta - TMath::Pi()/2) * sigmaTheta(data->piPos.theta - TMath::Pi()/2) / 2 * PsiFunc::Derivative(PsiFunc::Var::thetaPos, data->piPos, data->piNeg, 2) + 
                        sigmaTheta(data->piNeg.theta - TMath::Pi()/2) * sigmaTheta(data->piNeg.theta - TMath::Pi()/2) / 2 * PsiFunc::Derivative(PsiFunc::Var::thetaNeg, data->piPos, data->piNeg, 2) + 
                        sigmaPhi * sigmaPhi / 2 * PsiFunc::Derivative(PsiFunc::Var::phiPos, data->piPos, data->piNeg, 2) +
                        sigmaPhi * sigmaPhi / 2 * PsiFunc::Derivative(PsiFunc::Var::phiNeg, data->piPos, data->piNeg, 2);
        auto dpsi = tree->reco.ksdpsi + psiCor;

        auto energy = useTrueEnergy? tree->etrue : meanEnergy.value_or(tree->emeas);

        auto [binY, binKsTheta] = Sigmas::GetSigmaMatrix_bin(lnY, data->ks.theta);
        auto sigmaPsi = (binY != -1 && binKsTheta != -1)? vSigmaMatrixFit[binY][binKsTheta] : 0.;

        auto massCorr = -sigmaPsi * sigmaPsi / 2 * 
                        FullRecMassFunc::Derivative(FullRecMassFunc::Var::psi, dpsi, energy, data->Y, 2);
        auto massCorrY = sigmaY * sigmaY / 2 * 
                        FullRecMassFunc::Derivative(FullRecMassFunc::Var::Y, dpsi, energy, data->Y, 2);

        auto mass = FullRecMassFunc::Eval(dpsi, energy, data->Y) + massCorr + massCorrY;

        container["hDeltaM"]->Fill(lnY, massCorr);
        container["hMlnYpfx"]->Fill(lnY, mass); 

        container["hDiffRecGen"]->Fill(data->ks.theta - TMath::Pi() / 2, dpsi - tree->gen.ksdpsi);
        container["hMassVsKsTheta"]->Fill(data->ks.theta - TMath::Pi() / 2, mass); 
        auto eventType = misc::GetEventType(data->piPos, data->piNeg);
        if(eventType == misc::EventType::cowboy)
        { 
            container["hDiffRecGen_cowboy"]->Fill(data->ks.theta - TMath::Pi() / 2, dpsi - tree->gen.ksdpsi);
            container["hMlnYpfx_cowboy"]->Fill(lnY, mass); 
        }

        if(eventType == misc::EventType::sailor)
        { 
            container["hDiffRecGen_sailor"]->Fill(data->ks.theta - TMath::Pi() / 2, dpsi - tree->gen.ksdpsi);
            container["hMlnYpfx_sailor"]->Fill(lnY, mass); 
        }
        container["hPsilnY"]->Fill(lnY, dpsi); 

        if(fabs(lnY) < 0.3)
        { container["hEnergySpectrum"]->Fill(tree->etrue); }
    }
}

std::pair<double, double> HandlerMC::GetMass(double fitRange)
{   
    auto printFitRes = [](const TFitResultPtr &res, std::string resName) {
        std::cout << "\n" << resName << " = " << res->Parameter(0) << " +/- " << res->ParError(0) 
            << "; chi2 / ndf = " << res->Chi2() << "/" << res->Ndf() << "; Prob = " << res->Prob() << std::endl;              
    };

    auto res = container["hMlnYpfx"]->Fit("pol0", "SMQE", "goff", -fitRange, fitRange);
    if(verbose)
    { printFitRes(res, "Mass_FullRec"); }

    double mass = res->Parameter(0);
    double massErr = res->ParError(0);

    return std::pair(mass, massErr);
}

std::pair<double, double> HandlerMC::GetEnergySpectrumMean()
{
    auto mean = container["hEnergySpectrum"]->GetMean();
    auto meanErr = container["hEnergySpectrum"]->GetMeanError();
    return std::make_pair(mean, meanErr);
}

void HandlerMC::Draw(std::string name)
{
    painter->Draw(name);
}

void HandlerMC::Draw(std::string name, Range range)
{
    painter->Draw(name, "", range);
}

void HandlerMC::SaveHists(std::string output_filename)
{ container.Save(output_filename); }

void HandlerMC::SaveSplines(std::string spline_filename)
{
    deltaPsi_RecGenDiff.Save(spline_filename, true);
    sigmaTheta.Save(spline_filename, false);
}

#endif