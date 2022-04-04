#include "TH1D.h"
#include "TH2D.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TGraphErrors.h"

int KsMassMeas()
{
    double eNom = 509.5; // Nominal energy of beam, MeV
    double MinvPi = 139.57; // Pi+ Mass, Mev

    TFile *file = TFile::Open("hists and root files/hist_kskl.root");
    auto ksTr = (TTree*)file->Get("ksTree");
    auto hKsMass = new TH1D("hKsMass", "Ks mass", 960, 450, 550);
    TF1 *gaus = new TF1("fGaus", "gaus", 495.52, 498.54);
    auto h2EvsRun = new TH2D("h2EvsRun", "Measured energy vs run number", 1000, 60700, 61400, 1000, 509.3, 509.7);
    
    Float_t emeas; Float_t demeas; Int_t runnum; Float_t ksdpsi[7];
    ksTr->SetBranchAddress("emeas", &emeas);
    ksTr->SetBranchAddress("demeas", &demeas);
    ksTr->SetBranchAddress("runnum", &runnum);
    ksTr->SetBranchAddress("ksdpsi", ksdpsi);

    auto massFunc = new TF1("mass", "[0] * TMath::Sqrt(1 - (1 - 4 * [1] * [1] / [0] / [0]) * cos(x / 2) * cos(x / 2))");

    Long64_t n = ksTr->Draw("emeas : runnum","","goff");
    Int_t *entryNum = new Int_t[n];
    int tmp = 0; int nRunMax = 0;
    for(int i = 0; i < ksTr->GetEntriesFast(); i++)
    {
        ksTr->GetEntry(i);
        if(tmp != runnum)
        {
            entryNum[nRunMax] = i;
            nRunMax++;
        }
        tmp = runnum;
    }
    Double_t *v0 = new Double_t[nRunMax];
    Double_t *v1 = new Double_t[nRunMax];
    Double_t *eV0 = new Double_t[nRunMax];
    Double_t *eV1 = new Double_t[nRunMax];

    for(int i = 0; i < nRunMax; i++)
    {
        ksTr->GetEntry(entryNum[i]);
        v0[i] = runnum;
        v1[i] = emeas;
        eV0[i] = 0;
        eV1[i] = demeas;
    }
    TGraph *gEvsRun = new TGraphErrors(nRunMax, v0, v1, eV0, eV1);
    gEvsRun->Draw();
    delete [] v0; delete [] v1; delete [] eV0; delete [] eV1;
    Long64_t nSuccess = 0;
    Long64_t nentries = ksTr->GetEntries();

    for(int i = 0; i < nentries; i++)
    {
        ksTr->GetEntry(i);
        if(fabs(emeas - eNom) < 5)
        {
            massFunc->SetParameter(0, emeas);
            massFunc->SetParameter(1, MinvPi);
            hKsMass->Fill(massFunc->Eval(ksdpsi[0]));
            nSuccess++;
        }
    }

    for(int i = 0; i < 3; i++)
    { gaus->SetParameter(i, 1); }
    hKsMass->Fit("fGaus", "S", "", 495.52, 498.54);
    hKsMass->Draw();
    gaus->Draw("same");


    double mass = gaus->GetParameter(1);
    double massErr = gaus->GetParError(1);

    std::cout << "Mass: " << mass << " +- " << massErr << std::endl;
    std::cout << "Number of points: " << nentries << std::endl;
    std::cout << "Number of success points: " << nSuccess << std::endl;
    return 0;
} 