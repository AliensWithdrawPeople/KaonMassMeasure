#include "C:\work\Science\BINP\Kaon Mass Measure\Energy.h"

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();
    const std::vector<double> meanEnergies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp508.5.root", 508.397, 0.008, 15, 4.06);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp509.5.root", 509.529, 0.004, 15, 3.753);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp510.root", 509.957, 0.005, 15, 3.6);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/exкpKpKm/kchExp510.5.root", 510.454, 0.007, 15, 3.677);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp511.root", 511.035, 0.009, 30, 0);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/MC505_FSR.root", 504.8, 0.007, 15, 0);

    std::vector<Energy*> handlers;
    for(int i = 0; i < energyPoints.size(); i++)
    { handlers.push_back(new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/MC" + energyPoints[i] + "_FSR.root", meanEnergies[i], 0.007, 30, 0)); }

    std::cout << "{";
    for(auto &handler : handlers)
    {
        auto shift = handler->GetEnergyDiff();
        for (auto k : shift)
        { std::cout << -k.second << ", "; }
    }
    std::cout << "}" << std::endl;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}