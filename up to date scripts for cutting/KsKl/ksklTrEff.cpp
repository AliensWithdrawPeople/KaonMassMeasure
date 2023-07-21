#define ksklTrEff_cxx
#include "ksklTrEff.h"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>
#include <TLine.h>


void ksklTrEff::Loop(std::string outFileName)
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        
        if(nt != 1 && nt != 2 && is_coll == 1)
            continue;

        for(int i = 0; i < nt; i++)
        {
            if (tptot[i] > 40 && fabs(trho[i]) < 6 && fabs(tz[i]) < 10 &&
                tchi2r[i] < 15 && tchi2z[i] < 10 && tnhit[i] > 10)
            {
                
            }
        }

        
    }
}
