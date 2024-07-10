#include "Energy.hpp"

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();
    const std::vector<double> meanEnergies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    std::string badRuns = "C://work/Science/BINP/Kaon Mass Measure/txt/BadRuns.txt";
    std::map<std::string, Energy*> enDict;
    // for(int i = 0; i < energyPoints.size(); i++)
    // { 
    //     enDict[energyPoints[i]] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoints[i] + ".root", 
    //                                             badRuns, meanEnergies[i], meanEnergiesErr[i], 30, 3.4);        
    // }

    // MC
    // enDict["510"] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/Kch_MC510_smeared.root", 
    //                                             "", 509.956, 0.007, 30, 0, false);

    // No beam point
    std::string enPointName = "514";  
    // exp
    enDict[enPointName] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/KpKm_NoVertex/kchExp" + enPointName +  "_noVertex.root", 
                                                badRuns, meanEnergies[8], meanEnergiesErr[8], 30, 0);

    // MC                              
    // enDict[enPointName] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/MC_Kch" + enPointName +  "_NoBeamPoint.root", 
    //                                             "", meanEnergies[8], meanEnergiesErr[8], 30, 0, false);

    enDict[enPointName]->DrawGraph();
    std::vector<double> kchMeanEnergy;
    std::vector<double> kchMeanEnergyErr;
    // for(auto &handler : enDict)
    // {
    //     auto p = handler.second->DrawGraph();
    //     kchMeanEnergy.push_back(p.first);
    //     kchMeanEnergyErr.push_back(p.second);
    // }

    // for(auto en : kchMeanEnergy)
    // { std::cout<< en << ", "; }
    // std::cout<< std::endl;
    // for(auto err : kchMeanEnergyErr)
    // { std::cout<< err << ", "; }
    // std::cout<< std::endl;

    for(auto &handler : enDict)
    { delete handler.second; }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}