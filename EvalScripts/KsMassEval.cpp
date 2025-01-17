#include <chrono>
#include <ctime> 
#include <iostream> 

#include "HandlerMC.hpp"
#include "HandlerExp.hpp"


int KsMassEval()
{
    auto start = std::chrono::system_clock::now();

    const std::vector<double> meanEnergies_vec = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};
    // Means of Ks energy spectrums
    const std::vector<double> meanEnergiesSpectrum_vec = {504.683, 507.762, 508.323, 508.885, 509.445, 509.841, 510.263, 510.694,  510.694, 512.297};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009, 0.009};

    const std::vector<double> pion_theta_covariances = {0.000128597, 0.00015998, 0.000175285, 0.000159971, 0.000152738, 0.000153349, 0.000162581, 0.000144861, 0.000159156, 0.000152202};
    const std::vector<double> piPos_correction_constant = {-0.00699377, -0.00645158, -0.00664103, -0.00600068, -0.00584225, -0.00508524, -0.00560359, -0.0046329, -0.0050488, -0.00538691};
    const std::vector<double> piPos_correction_slope = {0.00443298, 0.00408992, 0.00410851, 0.00383344, 0.00373523, 0.00322678, 0.00349512, 0.00288437, 0.00325238, 0.00341808};
    const std::vector<double> piNeg_correction_constant = {-0.00683841, -0.00749158, -0.00641094, -0.00569952, -0.00557661, -0.00490508, -0.00509329, -0.00451928, -0.00465828, -0.00478692};
    const std::vector<double> piNeg_correction_slope = {0.00437792, 0.00472957, 0.00398921, 0.00359215, 0.0036446, 0.00310455, 0.00322322, 0.00279679, 0.00298795, 0.00309261};
    // Energy shift of Compton laser system
    const double energy_shift = 0.006;
    const double sigma_matrix_window_width = 2;

    const std::vector<double> deltaE_RC_Smeared = {0.105504, 0.0760349, 0.069423, 0.0598971, 0.0769491, 0.118673, 0.197142, 0.336072, 0.467132, 1.54063};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"};

    std::map<std::string, std::pair<double, double>> meanEnergies;
    std::map<std::string, std::pair<double, double>> meanEnergiesSpectrum;
    std::map<std::string, double> radiativeCorrections;
    std::map<std::string, double> pion_theta_covariance;
    std::map<std::string, std::pair<double, double>> piPos_correction;
    std::map<std::string, std::pair<double, double>> piNeg_correction;
    for(int i = 0; i < energyPoints.size(); i++)
    {   
        meanEnergies[energyPoints[i]] = std::make_pair(meanEnergies_vec[i] + energy_shift, meanEnergiesErr[i]); 
        meanEnergiesSpectrum[energyPoints[i]] = std::make_pair(meanEnergiesSpectrum_vec[i], meanEnergiesErr[i]); 
        radiativeCorrections[energyPoints[i]] = deltaE_RC_Smeared[i]; 
        pion_theta_covariance[energyPoints[i]] = pion_theta_covariances[i]; 
        piPos_correction[energyPoints[i]] = std::make_pair(piPos_correction_constant[i], piPos_correction_slope[i]);
        piNeg_correction[energyPoints[i]] = std::make_pair(piNeg_correction_constant[i], piNeg_correction_slope[i]);
    }

    std::string energyPoint = "511.5";
    
    std::string fileNameMC = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KsKl_Smeared/New formfactor/XsecConv/MC" + energyPoint + "_XsecConv.root";
    // // std::string fileNameMC = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KsKl_Smeared/phi_width 4.5 MeV/MC" + energyPoint + "_XsecConv.root";
    auto handlerMC = new HandlerMC(fileNameMC, energyPoint, 0.27, meanEnergies[energyPoint].first, true, false, true, sigma_matrix_window_width);
    auto [mass, massErr] = handlerMC->Eval();
    std::cout << mass << " + " << massErr << std::endl;

    auto [meanEnergy_Spectrum, meanEnergyErr_Spectrum] = handlerMC->GetEnergySpectrumMean();
    std::cout << meanEnergy_Spectrum << " + " << meanEnergyErr_Spectrum << std::endl;
    // handlerMC->Draw("hDiffRecGen", {-0.8, 0.8});
    handlerMC->Draw("hMlnYpfx", {-0.5, 0.5});
    handlerMC->SaveHists("C:/work/Science/BINP/Kaon Mass Measure/hists/MC/Hists_MC" + energyPoint + ".root");
    handlerMC->SaveSplines("C:/work/Science/BINP/Kaon Mass Measure/splines/spline_" + energyPoint + ".root");
    // delete handlerMC;
    // std::string fileNameExp = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/expKsKl/exp" + energyPoint + ".root";
    auto handlerExp = new HandlerExp(fileNameExp, energyPoint, 0.27, meanEnergies[energyPoint].first, 
                                    radiativeCorrections[energyPoint], pion_theta_covariance[energyPoint], 
                                    piPos_correction[energyPoint], piNeg_correction[energyPoint], 
                                    sigma_matrix_window_width, true);
    // auto [mass, massErr] = handlerExp->Eval();
    // std::cout << mass << " +/- " << massErr << std::endl;
    // handlerExp->Draw("hMassVsKsTheta", {-0.8, 0.8});
    // handlerExp->Draw("hMlnYpfx", {-0.5, 0.5});
    // handlerExp->SaveHists("C:/work/Science/BINP/Kaon Mass Measure/hists/Exp/Hists_Exp" + energyPoint + ".root");
    // delete handlerExp;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}