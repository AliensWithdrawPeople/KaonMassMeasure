#define PionTracksCutter_cxx
#include "PionTracksCutter.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>
#include <TLine.h>

#include <vector>
#include <iostream>

void PionTracksCutter::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    auto file = TFile::Open("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/pipi/Pion_TrackAnglesDiff_Raw03.root", "recreate");
    auto hThetaDiff_PosTr = new TH2D("hThetaDiff_PosTr", "#pi^{+}; #theta^{(rec)}_{#pi^{+}}; #theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}", 3200, 0, 3.2, 4000, -2, 2); 
    auto hThetaDiff_NegTr = new TH2D("hThetaDiff_NegTr", "#pi^{-}; #theta^{(rec)}_{#pi^{-}}; #theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}", 3200, 0, 3.2, 4000, -2, 2); 

    auto hPhiDiff_PosTr = new TH2D("hPhiDiff_PosTr", "#pi^{+}; #phi^{(rec)}_{#pi^{+}}; #phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}", 3200, -3.2, 3.2, 4000, -2, 2); 
    auto hPhiDiff_NegTr = new TH2D("hPhiDiff_NegTr", "#pi^{-}; #phi^{(rec)}_{#pi^{-}}; #phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}", 3200, -3.2, 3.2, 4000, -2, 2); 


    auto hThetaDiff_PosTrV = new TH2D("hThetaDiff_PosTrV", "#pi^{+} with a common vertex; #theta^{(rec)}_{#pi^{+}}; #theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}", 3200, 0, 3.2, 4000, -2, 2); 
    auto hThetaDiff_NegTrV = new TH2D("hThetaDiff_NegTrV", "#pi^{-} with a common vertex; #theta^{(rec)}_{#pi^{-}}; #theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}", 3200, 0, 3.2, 4000, -2, 2); 

    auto hPhiDiff_PosTrV = new TH2D("hPhiDiff_PosTrV", "#pi^{+} with a common vertex; #phi^{(rec)}_{#pi^{+}}; #phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}", 3200, -3.2, 3.2, 4000, -2, 2); 
    auto hPhiDiff_NegTrV = new TH2D("hPhiDiff_NegTrV", "#pi^{-} with a common vertex; #phi^{(rec)}_{#pi^{-}}; #phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}", 3200, -3.2, 3.2, 4000, -2, 2); 

    std::vector<int> goodTr = {};

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        for (int i = 0; i < nt; i++)
        {
            if (tptot[i] > 40 && fabs(tz[i]) < 10 &&
                tchi2r[i] < 15 && tchi2z[i] < 10 && tnhit[i] > 10)
            { goodTr.push_back(i); }

            if(goodTr.size() > 2)
            { break; }
        }
        
    
        if (goodTr.size() == 2 && is_coll == 1)
        {
            for (int i = 0; i < nt; i++)
            {
                if( tcharge[goodTr[0]] * tcharge[goodTr[1]] < 0) 
                {

                    int posTrackNumber = tcharge[goodTr[0]] > 0 ? 0 : 1;
                    int negTrackNumber = posTrackNumber == 1 ? 0 : 1;
                    bool hasPi0 = false;
                    for(int j = 0; j < nsim; j++)
                    {
                        if(simtype[j] == 111 && simorig[j] == 0)
                        { hasPi0 = true; }
                    }

                    for(int j = 0; j < nsim; j++)
                    {
                        if(simtype[j] == 211 && simorig[j] == 0)
                        {
                            hThetaDiff_PosTr->Fill(tth[goodTr[posTrackNumber]], tth[goodTr[posTrackNumber]] - simtheta[j]);
                            hPhiDiff_PosTr->Fill(tphi[goodTr[posTrackNumber]], tphi[goodTr[posTrackNumber]] - simphi[j]);

                            hThetaDiff_PosTrV->Fill(tthv[goodTr[posTrackNumber]], tthv[goodTr[posTrackNumber]] - simtheta[j]);
                            hPhiDiff_PosTrV->Fill(tphiv[goodTr[posTrackNumber]], tphiv[goodTr[posTrackNumber]] - simphi[j]);
                        }

                        if(simtype[j] == -211 && simorig[j] == 0)
                        {
                            hThetaDiff_NegTr->Fill(tth[goodTr[negTrackNumber]], tth[goodTr[negTrackNumber]] - simtheta[j]);
                            hPhiDiff_NegTr->Fill(tphi[goodTr[negTrackNumber]], tphi[goodTr[negTrackNumber]] - simphi[j]);

                            hThetaDiff_NegTrV->Fill(tthv[goodTr[negTrackNumber]], tthv[goodTr[negTrackNumber]] - simtheta[j]);
                            hPhiDiff_NegTrV->Fill(tphiv[goodTr[negTrackNumber]], tphiv[goodTr[negTrackNumber]] - simphi[j]);
                        }
                    }
                }
            }

        }

        goodTr.clear();
        goodTr.shrink_to_fit();
    }
    std::cout << "N_events: " << hThetaDiff_NegTr->GetEntries() << std::endl;
    file->Write();
    file->Save();
}
