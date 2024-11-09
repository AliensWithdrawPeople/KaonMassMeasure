#ifndef HandlerExp_h
#define HandlerExp_h

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
#include <stdexcept>

#include "TMath.h"
#include "TVector3.h"


class HandlerExp: Evaluator
{
private:
    std::string energyPoint;
    std::optional<double> meanEnergy;
    std::optional<double> energyCorrection;
    std::optional<double> pion_theta_covariance;
    std::optional<std::pair<double, double>> piPos_correction; 
    std::optional<std::pair<double, double>> piNeg_correction; 

    bool verbose;
    std::unique_ptr<Tree> tree;
    Tree::Data* data; 
    HistContainer container = HistContainer();
    std::unique_ptr<Painter> painter;
    
    double fitRange;
    /// @brief Width in sigmas in which sigma^{(M_K)}_{psi} [lnY][ksTheta] is evaluated.
    double sigma_matrix_window_width;
    bool useCorrectedEnergy = true;

    /// @brief Numbers of entries that satisfy all selection cuts.
    std::vector<int> goodEntries = {};
    /// @brief Sigma_psi(lnY, kstheta) from gaus fit.
    std::vector<std::vector<double>> vSigmaMatrixFit;


    /// @brief Correction to psi_reco (dpsi = psi_reco - psi_gen). 
    Spline deltaPsi_RecGenDiff;
    /// @brief Resolution of theta as a function of track's theta.
    Spline sigmaTheta;

    /// @brief Load ksdpsiRecGenDiff_spline and sigmaTheta_spline from the file.
    /// @param spline_filename name of .root file where splines are saved. 
    /// @return std::pair<ksdpsiRecGenDiff_spline, sigmaTheta_spline>.
    std::pair<Spline, Spline> LoadSplines(std::string spline_filename);
    double GetPionThetaCorrection_(double KsTheta, bool isPos);

    void FillHists(bool useCorrectedEnergy);

    static bool KsThetaCut(double ksTheta)
    { return fabs(ksTheta - TMath::Pi() / 2) < 0.3; }

public:
    HandlerExp(std::string fKsKl, 
            std::string energyPoint, 
            double fitRange, 
            std::optional<double> meanEnergy, 
            std::optional<double> energyCorrection, 
            std::optional<double> pion_theta_covariance, 
            std::optional<std::pair<double, double>> piPos_correction, 
            std::optional<std::pair<double, double>> piNeg_correction,
            double sigma_matrix_window_width, 
            bool useCorrectedEnergy, 
            bool isVerbose = true);


    std::pair<double, double> GetMass(double fitRange = 0.27);
    
    std::pair<double, double> Eval() override { return GetMass(); }
    void Draw(std::string name) override;
    void Draw(std::string name, Range range);
    void SaveHists(std::string output_filename);
};

HandlerExp::HandlerExp(std::string fKsKl, 
                    std::string energyPoint, 
                    double fitRange, 
                    std::optional<double> meanEnergy, 
                    std::optional<double> energyCorrection, 
                    std::optional<double> pion_theta_covariance, 
                    std::optional<std::pair<double, double>> piPos_correction, 
                    std::optional<std::pair<double, double>> piNeg_correction, 
                    double sigma_matrix_window_width, 
                    bool useCorrectedEnergy, 
                    bool isVerbose): 
    energyPoint{energyPoint}, fitRange{fitRange}, meanEnergy{meanEnergy}, energyCorrection{energyCorrection}, 
    useCorrectedEnergy{useCorrectedEnergy}, verbose{isVerbose}, pion_theta_covariance{pion_theta_covariance},
    piPos_correction{piPos_correction}, piNeg_correction{piNeg_correction}, sigma_matrix_window_width{sigma_matrix_window_width}
{
    tree = std::unique_ptr<Tree>(new Tree(fKsKl, true));
    painter = std::unique_ptr<Painter>(new Painter(&container));

    tree->SetReco(true);
    data = tree->data; 

    container.Add("hMlnY", new TH2D("hMlnY", "M(lnY)", 300, -0.8, 0.8, 600, 480, 520));
    container.Add("hDeltaM", new TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1));
    container.Add("MPsi", new TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520));
    container.Add("hM_CrAnglelnY", new TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 30, -0.4, 0.4, 40000, 490, 515));
    container.Add("hPsilnY", new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 10000, 2.4, 3.3));
    container.Add("hMassVsKsTheta", new TH2D("hMassVsKsTheta", "M vs KsTheta", 600, -1.57, 1.57, 40000, 480, 520));

    container.Add("hMassVsKstlen", new TProfile("hMassVsKstlen", "M vs tlen", 100, 0, 2, 490, 505));
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
    vSigmaMatrixFit = Sigmas::GetSigmaMatrix(tree, goodEntries, sigma_matrix_window_width, sigma_matrix_window_width, verbose); 
    std::string spline_filename = misc::GetSplineFilename(energyPoint);
    std::tie(deltaPsi_RecGenDiff, sigmaTheta) = LoadSplines(spline_filename);
    FillHists(useCorrectedEnergy);
}

double HandlerExp::GetPionThetaCorrection_(double ksTheta, bool isPos)
{
    if(!(piPos_correction.has_value() && piNeg_correction.has_value()))
    { return 0.; }
    if(isPos)
    {
        return -(piPos_correction.value().first + piPos_correction.value().second * ksTheta);
    }
    return -(piNeg_correction.value().first + piNeg_correction.value().second * ksTheta);
}

std::pair<Spline, Spline> HandlerExp::LoadSplines(std::string spline_filename)
{ return {Spline(spline_filename, "ksdpsiRecGenDiff_spline"), Spline(spline_filename, "sigmaTheta_spline")}; }

void HandlerExp::FillHists(bool useCorrectedEnergy)
{
    auto sigmaPhi = Sigmas::GetPhiSigma();
    auto sigmaY = Sigmas::GetYSigma();

    for(const auto &entry : goodEntries)
    {
        tree->GetEntry(entry);
        if(!KsThetaCut(data->ks.theta))
        { continue; }

        auto lnY = log(data->Y);
        data->piPos.theta += GetPionThetaCorrection_(data->ks.theta, true);
        data->piNeg.theta += GetPionThetaCorrection_(data->ks.theta, false);
        auto psiCor =   sigmaTheta(data->piPos.theta - TMath::Pi()/2) * sigmaTheta(data->piPos.theta - TMath::Pi()/2) / 2 * PsiFunc::Derivative(PsiFunc::Var::thetaPos, data->piPos, data->piNeg, 2) + 
                        sigmaTheta(data->piNeg.theta - TMath::Pi()/2) * sigmaTheta(data->piNeg.theta - TMath::Pi()/2) / 2 * PsiFunc::Derivative(PsiFunc::Var::thetaNeg, data->piPos, data->piNeg, 2) + 
                        pion_theta_covariance.value_or(0.) *  PsiFunc::MixedThetaDerivative(data->piPos, data->piNeg) +
                        sigmaPhi * sigmaPhi / 2 * PsiFunc::Derivative(PsiFunc::Var::phiPos, data->piPos, data->piNeg, 2) +
                        sigmaPhi * sigmaPhi / 2 * PsiFunc::Derivative(PsiFunc::Var::phiNeg, data->piPos, data->piNeg, 2);
        auto dpsi = PsiFunc::Eval(data->piPos.theta, data->piPos.phi, data->piNeg.theta,  data->piNeg.phi) + psiCor;

        auto energy = meanEnergy.value_or(tree->emeas);
        energy = useCorrectedEnergy? misc::GetCorrectedEnergy(tree->runnum, energy) : energy;
        energy = energy - energyCorrection.value_or(0.0);

        auto [binY, binKsTheta] = Sigmas::GetSigmaMatrix_bin(lnY, data->ks.theta);
        auto sigmaPsi = (binY != -1 && binKsTheta != -1)? vSigmaMatrixFit[binY][binKsTheta] : 0.;

        auto massCorr = -sigmaPsi * sigmaPsi / 2 * 
                        FullRecMassFunc::Derivative(FullRecMassFunc::Var::psi, dpsi, energy, data->Y, 2);
        auto massCorrY = sigmaY * sigmaY / 2 * 
                        FullRecMassFunc::Derivative(FullRecMassFunc::Var::Y, dpsi, energy, data->Y, 2);
    
        // This correction was obtained using MC truth.
        auto eventType = misc::GetEventType(data->piPos, data->piNeg);
        auto event_type_mass_correction = eventType == misc::EventType::sailor? -0.010 : 0.006; // MeV
        // It can be applied like that:
        // auto mass = FullRecMassFunc::Eval(dpsi, energy, data->Y) + massCorr + massCorrY + event_type_mass_correction;
        // I did it to estimate the systematic error of these shifts. Do not use it to obtain the results!!!

        auto mass = FullRecMassFunc::Eval(dpsi, energy, data->Y) + massCorr + massCorrY;

        container["hDeltaM"]->Fill(lnY, massCorr);
        container["hMlnYpfx"]->Fill(lnY, mass); 
        container["hMlnY"]->Fill(lnY, mass); 
        container["hMassVsKstlen"]->Fill(fabs(data->ks_len), mass); 

        if(fabs(lnY) < fitRange)
        { 
            container["hMassVsKsTheta"]->Fill(data->ks.theta - TMath::Pi() / 2, mass); 
            container["hMassVsKstlen"]->Fill(fabs(data->ks_len), mass); 
        }

        if(eventType == misc::EventType::cowboy)
        { 
            container["hMlnYpfx_cowboy"]->Fill(lnY, mass); 
        }

        if(eventType == misc::EventType::sailor)
        { container["hMlnYpfx_sailor"]->Fill(lnY, mass); }
        container["hPsilnY"]->Fill(lnY, dpsi); 
    }
}

std::pair<double, double> HandlerExp::GetMass(double fitRange)
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

void HandlerExp::Draw(std::string name)
{
    painter->Draw(name);
}

void HandlerExp::Draw(std::string name, Range range)
{ painter->Draw(name, "", range); }

void HandlerExp::SaveHists(std::string output_filename)
{ container.Save(output_filename); }

#endif