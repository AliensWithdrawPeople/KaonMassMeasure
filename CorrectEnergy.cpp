#include "TTree.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TH1D.h"

double CorrectEnergy()
{
    int entry = 24937;
    TFile *file = TFile::Open("tr_ph/New folder/tr_ph_run000000mcgpj.root");
    auto tr = (TTree *)file->Get("tr_ph");
    int nsim; Float_t emeas;
    Float_t mom[20]; Float_t phi[20]; Float_t theta[20];
    Int_t orig[20]; Int_t type[20];
    tr->SetBranchAddress("emeas", &emeas);
    tr->SetBranchAddress("nsim", &nsim);
    tr->SetBranchAddress("simmom", mom);
    tr->SetBranchAddress("simphi", phi);
    tr->SetBranchAddress("simtheta", theta);
    tr->SetBranchAddress("simorig", orig);
    tr->SetBranchAddress("simtype", type);

    auto hNph = new TH1D("hNph", "N photons", 20, 0, 20);
    auto hE = new TH1D("hE", "Corrected energy", 100, 480, 520);
    std::vector<ROOT::Math::PxPyPzEVector> ph;
    ROOT::Math::PxPyPzEVector phTot;
    ROOT::Math::PxPyPzEVector ks;
    ROOT::Math::PxPyPzEVector def; 
    ROOT::Math::PxPyPzEVector newCMsys; 

    int nph = 0;
    for(int j = 0; j < 10000; j++)
    {
        tr->GetEntry(j);
        for(int i = 0; i < nsim; i++)
        {
            // orig == 0 <=> from beam && type == 22 <=> photon
            if(orig[i] == 0 && type[i] == 22)
            {  
                nph++;
                ph.push_back(ROOT::Math::PxPyPzEVector(mom[i] * sin(theta[i]) * cos(phi[i]), mom[i] * sin(theta[i]) * sin(phi[i]), mom[i] * cos(theta[i]), mom[i]));
            }
            if(type[i] == 310)
            { ks = ROOT::Math::PxPyPzEVector(mom[i] * sin(theta[i]) * cos(phi[i]), mom[i] * sin(theta[i]) * sin(phi[i]), mom[i] * cos(theta[i]), 0); }
        }
        if(!ph.empty())
        {
            for(auto v : ph)
            { phTot += v; }
            def = ROOT::Math::PxPyPzEVector(0, 0, 0, 2 * emeas);
            ROOT::Math::Boost boostToCM((def - phTot).BoostToCM());
            newCMsys = def - phTot;

            hE->Fill(newCMsys.E() / 2);
            hNph->Fill(nph);
            nph = 0;
            phTot = ROOT::Math::PxPyPzEVector(0, 0, 0, 0);
        }
        ph.clear();
        ph.shrink_to_fit();
    }
    hE->Draw();    
    return 0;
}