#include "TTree.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TH1D.h"

double CorrectEnergy()
{
    TFile *file = TFile::Open("tr_ph/New folder/mcgpj/tr_ph_run000019.root");
    // TFile *file = TFile::Open("tr_ph/New folder/tr_ph_run2bGen505.root");
    auto tr = (TTree *)file->Get("tr_ph");
    int nsim; Float_t emeas;
    Float_t mom[50]; Float_t phi[50]; Float_t theta[50];
    Int_t orig[50]; Int_t type[50];
    tr->SetBranchAddress("emeas", &emeas);
    tr->SetBranchAddress("nsim", &nsim);
    tr->SetBranchAddress("simmom", mom);
    tr->SetBranchAddress("simphi", phi);
    tr->SetBranchAddress("simtheta", theta);
    tr->SetBranchAddress("simorig", orig);
    tr->SetBranchAddress("simtype", type);

    auto hNph = new TH1D("hNph", "N photons", 20, 0, 20);
    auto hE = new TH1D("hE", "Corrected energy", 50000, 480, 515);
    auto hMom = new TH1D("hMomPi", "Pi+ mom", 700, 84, 88);
    auto hMomRel = new TH2D("hMomRel", "Pi+ mom vs Pi- mom", 600, 0, 350, 600, 0, 350);
    std::vector<ROOT::Math::PxPyPzEVector> ph;
    ROOT::Math::PxPyPzEVector phTot;
    ROOT::Math::PxPyPzEVector ks;
    ROOT::Math::PxPyPzEVector def; 
    ROOT::Math::PxPyPzEVector newCMsys; 
    double piPosMom = 0;
    double piNegMom = 0;

    int nph = 0;
    for(int j = 0; j < tr->GetEntriesFast(); j++)
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
            { 
                ks = ROOT::Math::PxPyPzEVector(mom[i] * sin(theta[i]) * cos(phi[i]), mom[i] * sin(theta[i]) * sin(phi[i]), mom[i] * cos(theta[i]), 0);
                hMom->Fill(mom[i]);
            }

            if(type[i] == 211 && orig[i] == 310)
            { piPosMom = mom[i]; }

            if(type[i] == -211 && orig[i] == 310)
            { piNegMom = mom[i]; }
        }
        hMomRel->Fill(piPosMom, piNegMom);
        if(piNegMom != 0 && piPosMom != 0)
        { hE->Fill(sqrt(139.57 * 139.57 + piNegMom * piNegMom) + sqrt(139.57 * 139.57 + piPosMom * piPosMom)); }
        piPosMom = 0;
        piNegMom = 0;
        if(!ph.empty())
        {
            for(auto v : ph)
            { phTot += v; }
            def = ROOT::Math::PxPyPzEVector(0, 0, 0, 2 * emeas);
            ROOT::Math::Boost boostToCM((def - phTot).BoostToCM());
            newCMsys = def - phTot;

            // hE->Fill(newCMsys.E() / 2);
            hNph->Fill(nph);
            nph = 0;
            phTot = ROOT::Math::PxPyPzEVector(0, 0, 0, 0);
        }
        ph.clear();
        ph.shrink_to_fit();
    }
    hE->Draw();
    // hMomRel->Draw("Col");    
    // hMom->Draw();    
    std::cout << "Energy point: " << emeas << "; Mean Energy = " << hE->GetMean() << " +/- " << hE->GetMeanError() << std::endl;
    return hNph->GetMean();
}