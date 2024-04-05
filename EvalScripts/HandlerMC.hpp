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
    Tree::Data* data; 
    // const Tree::Data* data; 
    HistContainer container = HistContainer();
    std::unique_ptr<Painter> painter;

    double fitRange;
    bool useTrueEnergy = false;
    bool useEnergySmearing = true;

    std::pair<double, double> piPos_corr = {0, 0};
    std::pair<double, double> piNeg_corr = {0, 0};
    /// @brief Numbers of entries that satisfy all selection cuts.
    std::vector<int> goodEntries = {};
    /// @brief Sigma_psi(lnY, kstheta) from gaus fit.
    std::vector<std::vector<double>> vSigmaMatrixFit;
    double pion_theta_covariance;

    /// @brief Correction to psi_reco (dpsi = psi_reco - psi_gen).
    Spline deltaPsi_RecGenDiff;
    /// @brief Resolution of theta as a function of track's theta.
    Spline sigmaTheta;

    /// @brief Creates deltaPsi_RecGenDiff -- correction to psi_reco (dpsi = psi_reco - psi_gen).
    /// @param tree reference to Tree object storing data about events,
    Spline CreateDeltaPsiSpline();
    double GetPionThetaCovariance();
    double GetPionThetaCorrection_(double KsTheta, bool isPos, bool isWriteSession = false);
    void FillHists(bool useTrueEnergy);

    static bool KsThetaCut(double ksTheta)
    { return fabs(ksTheta - TMath::Pi() / 2) < 0.3; }

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

    container.Add("hMlnY", new TH2D("hMlnY", "M(lnY)", 800, -0.8, 0.8, 1000, 485, 515));
    container.Add("hDeltaM", new TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1));
    container.Add("hPsilnY", new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 10000, 2.4, 3.3));
    container.Add("hEnergySpectrum", new TH1D("hEnergySpectrum", "hEnergySpectrum", 6000, 480, 540));
    container.Add("hMassVsKsTheta", new TH2D("hMassVsKsTheta", "M vs KsTheta", 600, -1.57, 1.57, 40000, 480, 520));

    container.Add("hDiffRecGen", new TH2D("hDiffRecGen", "hDiffRecGen", 1200, -6., 6., 20000, -1, 1));    
    container.Add("hDiffRecGen_cowboy", new TH2D("hDiffRecGen_cowboy", "hDiffRecGen_cowboy", 1200, -6., 6., 20000, -1, 1));
    container.Add("hDiffRecGen_sailor", new TH2D("hDiffRecGen_sailor", "hDiffRecGen_sailor", 1200, -6., 6., 20000, -1, 1));

    container.Add("hMassVsKstlen", new TProfile("hMassVsKstlen", "M vs tlen", 100, 0, 2, 490, 505));
    container.Add("hMlnYpfx", new TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505));
    container.Add("hMlnYpfx_cowboy", new TProfile("hMlnYpfx_cowboy", "Profile of M versus lnY, cowboy", 30, -1, 1, 490, 505));
    container.Add("hMlnYpfx_sailor", new TProfile("hMlnYpfx_sailor", "Profile of M versus lnY, sailor", 30, -1, 1, 490, 505));

    container.Add("hThetaPionDiffRecGen", new TH2D("hThetaPionDiffRecGen", "hThetaPionDiffRecGen", 2000, -1., 1., 2000, -1, 1));
    container.Add("hPhiPionDiffRecGen", new TH2D("hPhiPionDiffRecGen", "hPhiPionDiffRecGen", 2000, -1., 1., 2000, -1, 1));
    container.Add("hPhiThetaPionDiffRecGen", new TH2D("hPhiThetaPionDiffRecGen", "hPhiThetaPionDiffRecGen", 2000, -1., 1., 2000, -1, 1));

    container.Add("hThetaDiffRecGenVsTheta_PionPos", new TH2D("hThetaDiffRecGenVsTheta_PionPos", "hThetaDiffRecGenVsTheta_PionPos", 2000, 0., 3.14, 2000, -1, 1));
    container.Add("hThetaDiffRecGenVsTheta_PionNeg", new TH2D("hThetaDiffRecGenVsTheta_PionNeg", "hThetaDiffRecGenVsTheta_PionNeg", 2000, 0., 3.14, 2000, -1, 1));

    container.Add("hPhiDiffRecGenVsTheta_PionNeg", new TH2D("hPhiDiffRecGenVsTheta_PionNeg", "hPhiDiffRecGenVsTheta_PionNeg", 2000, 0., 3.14, 2000, -1, 1));
    container.Add("hPhiDiffRecGenVsTheta_PionPos", new TH2D("hPhiDiffRecGenVsTheta_PionPos", "hPhiDiffRecGenVsTheta_PionPos", 2000, 0., 3.14, 2000, -1, 1));

    container["hPhiDiffRecGenVsTheta_PionPos"]->GetYaxis()->SetTitle("#phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}, rad");
    container["hPhiDiffRecGenVsTheta_PionPos"]->GetXaxis()->SetTitle("#theta_{K_{S}}, rad");
    container["hPhiDiffRecGenVsTheta_PionNeg"]->GetYaxis()->SetTitle("#phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}, rad");
    container["hPhiDiffRecGenVsTheta_PionNeg"]->GetXaxis()->SetTitle("#theta_{K_{S}}, rad");

    container["hThetaDiffRecGenVsTheta_PionPos"]->GetYaxis()->SetTitle("#theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}, rad");
    container["hThetaDiffRecGenVsTheta_PionPos"]->GetXaxis()->SetTitle("#theta_{K_{S}}, rad");
    container["hThetaDiffRecGenVsTheta_PionNeg"]->GetYaxis()->SetTitle("#theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}, rad");
    container["hThetaDiffRecGenVsTheta_PionNeg"]->GetXaxis()->SetTitle("#theta_{K_{S}}, rad");

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
    pion_theta_covariance = GetPionThetaCovariance();
    std::string spline_filename = misc::GetSplineFilename(energyPoint);
    deltaPsi_RecGenDiff = CreateDeltaPsiSpline();
    auto tmp = GetPionThetaCorrection_(0., true, true);
    if(saveSplines)
    { SaveSplines(spline_filename); }

    FillHists(useTrueEnergy);
}

double HandlerMC::GetPionThetaCovariance()
{
    TH2D tmp_hist("tmp_hist", "tmp_hist", 2000, -1., 1., 2000, -1, 1);
    for(const auto &entry : goodEntries)
    {
        tree->GetEntry(entry);
        if(!KsThetaCut(data->ks.theta) || (useEnergySmearing? false : fabs(tree->emeas - meanEnergy.value_or(0)) > 20e-3))
        { continue; }
        tmp_hist.Fill(data->piPos.theta - tree->gen.piPos.theta, data->piNeg.theta - tree->gen.piNeg.theta); 
    }
    tmp_hist.SetAxisRange(-0.06, 0.06, "X");
    tmp_hist.SetAxisRange(-0.06, 0.06, "Y");
    return tmp_hist.GetCovariance();
}

double HandlerMC::GetPionThetaCorrection_(double ksTheta, bool isPos, bool isWriteSession)
{
    if(isWriteSession)
    {
        auto tmp_hThetaDiffRecGenVsTheta_PionPos = new TProfile("tmp_hThetaDiffRecGenVsTheta_PionPos", "hThetaDiffRecGenVsTheta_PionPos", 200, 0., 3.14, -1, 1);
        auto tmp_hThetaDiffRecGenVsTheta_PionNeg = new TProfile("tmp_hThetaDiffRecGenVsTheta_PionNeg", "hThetaDiffRecGenVsTheta_PionNeg", 200, 0., 3.14, -1, 1);

        TProfile tmp_hPhiDiffRecGenVsTheta_PionNeg("tmp_hPhiDiffRecGenVsTheta_PionNeg", "hPhiDiffRecGenVsTheta_PionNeg", 200, 0., 3.14, -1, 1);
        TProfile tmp_hPhiDiffRecGenVsTheta_PionPos("tmp_hPhiDiffRecGenVsTheta_PionPos", "hPhiDiffRecGenVsTheta_PionPos", 200, 0., 3.14, -1, 1);
        for(const auto &entry : goodEntries)
        {
            tree->GetEntry(entry);
            if((useEnergySmearing? false : fabs(tree->emeas - meanEnergy.value_or(0)) > 20e-3))
            { continue; }

            tmp_hPhiDiffRecGenVsTheta_PionPos.Fill(data->ks.theta, data->piPos.phi - tree->gen.piPos.phi); 
            tmp_hPhiDiffRecGenVsTheta_PionNeg.Fill(data->ks.theta, data->piNeg.phi - tree->gen.piNeg.phi); 
            tmp_hThetaDiffRecGenVsTheta_PionNeg->Fill(data->ks.theta, data->piPos.theta - tree->gen.piPos.theta); 
            tmp_hThetaDiffRecGenVsTheta_PionPos->Fill(data->ks.theta, data->piNeg.theta - tree->gen.piNeg.theta); 
        }
        auto res_pos = tmp_hThetaDiffRecGenVsTheta_PionPos->Fit("pol1", "SQME", "goff", TMath::Pi()/2 - 0.8, TMath::Pi()/2 + 0.8);
        auto res_neg = tmp_hThetaDiffRecGenVsTheta_PionNeg->Fit("pol1", "SQME", "goff", TMath::Pi()/2 - 0.8, TMath::Pi()/2 + 0.8);
        piPos_corr = std::make_pair(res_pos->Parameter(0), res_pos->Parameter(1));
        piNeg_corr = std::make_pair(res_neg->Parameter(0), res_neg->Parameter(1));
        std::cout << piPos_corr.first << " : " << piPos_corr.second << std::endl;
        std::cout << piNeg_corr.first << " : " << piNeg_corr.second << std::endl;
        return 0;
    }
    if(isPos)
    {
        return -(piPos_corr.first + piPos_corr.second * ksTheta);
    }
    return -(piNeg_corr.first + piNeg_corr.second * ksTheta);
}

Spline HandlerMC::CreateDeltaPsiSpline()
{
    std::vector<Double_t> psi_delta = {};
    std::vector<Double_t> ksTheta_bins = {};
    TProfile psi("psi", "deltaPsi", 50, 0, 3.14, -0.1, 0.1);

    for(const auto &entry : goodEntries)
    {
        tree->GetEntry(entry);
        if(abs(tree->reco.Y - 1) > 1e-9 && tree->reco.piPos.nhit > 10 && tree->reco.piNeg.nhit > 10 && KsThetaCut(tree->reco.ks.theta))
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
        // Ks in a good region + optional cut to eliminate energy smearing.
        if(!KsThetaCut(data->ks.theta) || (useEnergySmearing? false : fabs(tree->emeas - meanEnergy.value_or(0)) > 20e-3))
        { continue; }
        data->piPos.theta += GetPionThetaCorrection_(data->ks.theta, true);
        data->piNeg.theta += GetPionThetaCorrection_(data->ks.theta, false);
        auto psiCor =   sigmaTheta(data->piPos.theta - TMath::Pi()/2) * sigmaTheta(data->piPos.theta - TMath::Pi()/2) / 2 * PsiFunc::Derivative(PsiFunc::Var::thetaPos, data->piPos, data->piNeg, 2) + 
                        sigmaTheta(data->piNeg.theta - TMath::Pi()/2) * sigmaTheta(data->piNeg.theta - TMath::Pi()/2) / 2 * PsiFunc::Derivative(PsiFunc::Var::thetaNeg, data->piPos, data->piNeg, 2) +
                        pion_theta_covariance *  PsiFunc::MixedThetaDerivative(data->piPos, data->piNeg) +
                        sigmaPhi * sigmaPhi / 2 * PsiFunc::Derivative(PsiFunc::Var::phiPos, data->piPos, data->piNeg, 2) +
                        sigmaPhi * sigmaPhi / 2 * PsiFunc::Derivative(PsiFunc::Var::phiNeg, data->piPos, data->piNeg, 2);
        // data->piPos.theta = data->piPos.theta + GetPionThetaCorrection_(data->ks.theta, true, false);
        // data->piNeg.theta = data->piNeg.theta + GetPionThetaCorrection_(data->ks.theta, false, false);
        auto dpsi = PsiFunc::Eval(data->piPos.theta, data->piPos.phi, data->piNeg.theta,  data->piNeg.phi) + psiCor;

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
        container["hMlnY"]->Fill(lnY, mass);

        container["hDiffRecGen"]->Fill(data->ks.theta - TMath::Pi() / 2, dpsi - tree->gen.ksdpsi);
        if(fabs(lnY) < fitRange)
        { 
            container["hMassVsKsTheta"]->Fill(data->ks.theta - TMath::Pi() / 2, mass); 
            container["hMassVsKstlen"]->Fill(fabs(data->ks_len), mass); 
        }

        container["hThetaPionDiffRecGen"]->Fill(data->piPos.theta - tree->gen.piPos.theta, data->piNeg.theta - tree->gen.piNeg.theta);
        container["hPhiPionDiffRecGen"]->Fill(data->piPos.phi - tree->gen.piPos.phi, data->piNeg.phi - tree->gen.piNeg.phi); 
        container["hPhiThetaPionDiffRecGen"]->Fill(data->piPos.phi - tree->gen.piPos.phi, data->piPos.theta - tree->gen.piPos.theta); 
        
        container["hPhiDiffRecGenVsTheta_PionNeg"]->Fill(data->ks.theta, data->piPos.phi - tree->gen.piPos.phi); 
        container["hPhiDiffRecGenVsTheta_PionPos"]->Fill(data->ks.theta, data->piNeg.phi - tree->gen.piNeg.phi); 
        container["hThetaDiffRecGenVsTheta_PionPos"]->Fill(data->ks.theta, data->piPos.theta - tree->gen.piPos.theta); 
        container["hThetaDiffRecGenVsTheta_PionNeg"]->Fill(data->ks.theta, data->piNeg.theta - tree->gen.piNeg.theta); 

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