#include "C:\work\Science\BINP\Kaon Mass Measure\Energy.h"

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();
    const std::vector<double> meanEnergies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    std::map<std::string, Energy*> enDict;
    for(int i = 0; i < energyPoints.size(); i++)
    { 
        enDict[energyPoints[i]] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoints[i] + ".root", 
                                                "C://work/Science/BINP/Kaon Mass Measure/txt/BadRuns.txt", meanEnergies[i], meanEnergiesErr[i], 30, 0);        
    }

    // MC
    // enDict["510"] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/Kch_MC510_smeared.root", 
    //                                             "", 509.956, 0.007, 30, 0, false);
    
    enDict["510"]->DrawGraph();
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