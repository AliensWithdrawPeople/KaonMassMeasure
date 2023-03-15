#include "C:\work\Science\BINP\Kaon Mass Measure\Energy.h"

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp508.5.root", 508.397, 0.008, 15, 4.06);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp509.5.root", 509.529, 0.004, 15, 3.753);
    auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp510.root", 509.957, 0.005, 15, 3.6);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/exкpKpKm/kchExp510.5.root", 510.454, 0.007, 15, 3.677);
    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp511.root", 511.035, 0.009, 30, 0);

    // auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/MC_kch510.5.root", 510.454, 0.007, 15, 0);
    
    enHandler->DrawGraph();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}