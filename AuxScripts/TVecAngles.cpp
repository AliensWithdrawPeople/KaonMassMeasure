#include "TVector3.h"
#include "TSpline.h"
#include "TRandom.h"
#include "TH1D.h"

#include <iostream>

int TVecAngles()
{
    auto hist1 = new TH1D("hist", "Wo shift", 10000, 0, 3.15);
    auto hist2 = new TH1D("hist", "With shift", 10000, 0, 3.15);
    TVector3 v1(1, 0, 0);
    TVector3 v2(1, 0, 0);
    auto rand = new TRandom(1);

    for(int i = 0; i < 1e5; i++)
    {
        v1.SetMagThetaPhi(1, TMath::Pi() / 2, rand->Gaus(2.7, 0.1));
        v2.SetMagThetaPhi(1, TMath::Pi() / 2, rand->Gaus(0, 0.1));
        hist1->Fill(v1.Angle(v2));

        v1.SetPhi(v1.Phi() + 5e-3);
        hist2->Fill(v1.Angle(v2));
        // hist->Fill(acos(cos(v1.Phi() - v2.Phi())));
    }

    std::cout << "Wo shift Mean = " << hist1->GetMean() << " +/- " << hist1->GetMeanError() << std::endl;
    std::cout << "With shift Mean = " << hist2->GetMean() << " +/- " << hist2->GetMeanError() << std::endl;
    std::cout << "delta = " << hist2->GetMean() - hist1->GetMean() << std::endl; 
    // hist1->DrawClone();

    return 0;
}
