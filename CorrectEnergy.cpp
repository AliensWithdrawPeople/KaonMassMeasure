#include "TTree.h"
#include "Math/Vector4D.h"
#include "TH1D.h"

double CorrectEnergy()
{
    int entry = 24937;
    TFile *file = TFile::Open("tr_ph/New folder/tr_ph_run000000.root");
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

    int nph = 0;
    for(int j = 0; tr->GetEntriesFast(); j++)
    {
        tr->GetEntry(j);
        for(int i = 0; i < nsim; i++)
        {
            if(orig[i] == 0 && type[i] == 22) // orig == 0 <=> from beam; type == 22 <=> photon
            {  
                nph++;
            }
        }
        hNph->Fill(nph);
    }
    hNph->Draw();
    /*
    double px = mom * sin(theta) * cos(phi);
    double py = mom * sin(theta) * sin(phi);
    double pz = mom * cos(theta);
    ROOT::Math::PtEtaPhiEVector v1(mom * sin(theta) * cos(phi), mom * sin(theta) * sin(phi), mom * cos(theta), mom);
    */
    return 0;
}