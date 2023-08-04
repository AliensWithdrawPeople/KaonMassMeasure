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
#include <TProfile.h>
#include <TAxis.h>

#include <vector>
#include <iostream>

void PionTracksCutter::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    auto file = TFile::Open("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/pipi/Pion_TrackAnglesDiff1.root", "recreate");

    auto hThetaAvgDiff = new TH2D("hThetaAvgDiff", "; #theta^{(rec)}_{avg}; #theta^{(rec)}_{avg} - #theta^{(gen)}_{avg}", 3200, 3.2, 3.2, 4000, -2, 2); 
    auto hPhiAvgDiff = new TH2D("hPhiAvgDiff", "; #phi^{(rec)}_{avg}; #phi^{(rec)}_{avg} - #phi^{(gen)}_{avg}", 3200, 0, 6.4, 4000, -2, 2); 

    auto hThetaAvgDiffV = new TH2D("hThetaAvgDiffV", "common vertex fit; #theta^{(rec)}_{avg}; #theta^{(rec)}_{avg} - #theta^{(gen)}_{avg}", 3200, 3.2, 3.2, 4000, -2, 2); 
    auto hPhiAvgDiffV = new TH2D("hPhiAvgDiffV", "common vertex fit; #phi^{(rec)}_{avg}; #phi^{(rec)}_{avg} - #phi^{(gen)}_{avg}", 3200, 0, 6.4, 4000, -2, 2); 

    auto hThetaDiff_PosTr = new TH2D("hThetaDiff_PosTr", "#pi^{+}; #theta^{(rec)}_{#pi^{+}}; #theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}", 3200, 0, 3.2, 4000, -2, 2); 
    auto hThetaDiff_NegTr = new TH2D("hThetaDiff_NegTr", "#pi^{-}; #theta^{(rec)}_{#pi^{-}}; #theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}", 3200, 0, 3.2, 4000, -2, 2); 

    auto hPhiDiff_PosTr = new TH2D("hPhiDiff_PosTr", "#pi^{+}; #phi^{(rec)}_{#pi^{+}}; #phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}", 3200, 0, 6.4, 4000, -2, 2); 
    auto hPhiDiff_NegTr = new TH2D("hPhiDiff_NegTr", "#pi^{-}; #phi^{(rec)}_{#pi^{-}}; #phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}", 3200, 0, 6.4, 4000, -2, 2); 


    auto hThetaDiff_PosTrV = new TH2D("hThetaDiff_PosTrV", "#pi^{+} with a common vertex; #theta^{(rec)}_{#pi^{+}}; #theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}", 3200, 0, 3.2, 4000, -2, 2); 
    auto hThetaDiff_NegTrV = new TH2D("hThetaDiff_NegTrV", "#pi^{-} with a common vertex; #theta^{(rec)}_{#pi^{-}}; #theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}", 3200, 0, 3.2, 4000, -2, 2); 

    auto hPhiDiff_PosTrV = new TH2D("hPhiDiff_PosTrV", "#pi^{+} with a common vertex; #phi^{(rec)}_{#pi^{+}}; #phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}", 3200, 0, 6.4, 4000, -2, 2); 
    auto hPhiDiff_NegTrV = new TH2D("hPhiDiff_NegTrV", "#pi^{-} with a common vertex; #phi^{(rec)}_{#pi^{-}}; #phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}", 3200, 0, 6.4, 4000, -2, 2); 

    std::vector<int> goodTr = {};

    double thetaAvg, thetaAvgV;
    double phiAvg, phiAvgV;
    double thetaPos_MC, thetaNeg_MC;
    double phiPos_MC, phiNeg_MC;

    TVector3 piPos_Rec(1, 1, 1);
    TVector3 piPos_Gen(1, 1, 1);

    TVector3 piNeg_Rec(1, 1, 1);
    TVector3 piNeg_Gen(1, 1, 1);

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

                            piPos_Rec.SetMagThetaPhi(1, tth[goodTr[posTrackNumber]], tphi[goodTr[posTrackNumber]]);
                            piPos_Gen.SetMagThetaPhi(1, simtheta[j], simphi[j]);
                            thetaPos_MC = simtheta[j];
                            phiPos_MC = simphi[j];
                            hThetaDiff_PosTr->Fill(tth[goodTr[posTrackNumber]], tth[goodTr[posTrackNumber]] - simtheta[j]);
                            hPhiDiff_PosTr->Fill(tphi[goodTr[posTrackNumber]], piPos_Rec.XYvector().DeltaPhi(piPos_Gen.XYvector()));

                            hThetaDiff_PosTrV->Fill(tthv[goodTr[posTrackNumber]], tthv[goodTr[posTrackNumber]] - simtheta[j]);
                            hPhiDiff_PosTrV->Fill(tphiv[goodTr[posTrackNumber]], tphiv[goodTr[posTrackNumber]] - simphi[j]);
                        }

                        if(simtype[j] == -211 && simorig[j] == 0)
                        {
                            piNeg_Rec.SetMagThetaPhi(1, tth[goodTr[negTrackNumber]], tphi[goodTr[negTrackNumber]]);
                            piNeg_Gen.SetMagThetaPhi(1, simtheta[j], simphi[j]);

                            thetaNeg_MC = simtheta[j];
                            phiNeg_MC = simphi[j];

                            hThetaDiff_NegTr->Fill(tth[goodTr[negTrackNumber]], tth[goodTr[negTrackNumber]] - simtheta[j]);
                            hPhiDiff_NegTr->Fill(tphi[goodTr[negTrackNumber]], piNeg_Rec.XYvector().DeltaPhi(piNeg_Gen.XYvector()));

                            hThetaDiff_NegTrV->Fill(tthv[goodTr[negTrackNumber]], tthv[goodTr[negTrackNumber]] - simtheta[j]);
                            hPhiDiff_NegTrV->Fill(tphiv[goodTr[negTrackNumber]], tphiv[goodTr[negTrackNumber]] - simphi[j]);
                        }
                    }

                    thetaAvg = (tth[goodTr[posTrackNumber]] - tth[goodTr[negTrackNumber]]) / 2;
                    thetaAvgV = (tthv[goodTr[posTrackNumber]] - tthv[goodTr[negTrackNumber]]) / 2;
                    phiAvg = (tphi[goodTr[posTrackNumber]] + tphi[goodTr[negTrackNumber]]) / 2;
                    phiAvgV = (tphiv[goodTr[posTrackNumber]] + tphiv[goodTr[negTrackNumber]]) / 2;

                    hThetaAvgDiff->Fill(thetaAvg, thetaAvg - (thetaPos_MC - thetaNeg_MC) / 2);
                    hPhiAvgDiff->Fill(phiAvg, phiAvg - (phiPos_MC + phiNeg_MC) / 2);
                    hThetaAvgDiffV->Fill(thetaAvgV, thetaAvgV - (thetaPos_MC - thetaNeg_MC) / 2);
                    hPhiAvgDiffV->Fill(phiAvgV, phiAvgV - (phiPos_MC + phiNeg_MC) / 2);
                }
            }

        }

        goodTr.clear();
        goodTr.shrink_to_fit();
    }

    auto pfx_ThetaDiff_PosTr = hThetaDiff_PosTr->ProfileX();
    auto pfx_PhiDiff_PosTr =  hPhiDiff_PosTr->ProfileX();
    auto pfx_ThetaDiff_PosTrV = hThetaDiff_PosTrV->ProfileX();
    auto pfx_PhiDiff_PosTrV = hPhiDiff_PosTrV->ProfileX();
    
    pfx_ThetaDiff_PosTr->Rebin(8);
    pfx_PhiDiff_PosTr->Rebin(8);
    pfx_ThetaDiff_PosTrV->Rebin(8);
    pfx_PhiDiff_PosTrV->Rebin(8);

    pfx_ThetaDiff_PosTr->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV; #theta^{(rec)}_{#pi^{+}}; #theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}");
    pfx_PhiDiff_PosTr->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV; #phi^{(rec)}_{#pi^{+}}; #phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}");
    pfx_ThetaDiff_PosTrV->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV, common vertex reconstruction; #theta^{(rec)}_{#pi^{+}}; #theta^{(rec)}_{#pi^{+}} - #theta^{(gen)}_{#pi^{+}}");
    pfx_PhiDiff_PosTrV->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV, common vertex reconstruction; #phi^{(rec)}_{#pi^{+}}; #phi^{(rec)}_{#pi^{+}} - #phi^{(gen)}_{#pi^{+}}");

    pfx_ThetaDiff_PosTr->GetYaxis()->SetRangeUser(-0.02, 0.02);
    pfx_PhiDiff_PosTr->GetYaxis()->SetRangeUser(-0.02, 0.02);
    pfx_ThetaDiff_PosTrV->GetYaxis()->SetRangeUser(-0.02, 0.02);
    pfx_PhiDiff_PosTrV->GetYaxis()->SetRangeUser(-0.02, 0.02);



    auto pfx_ThetaDiff_NegTr = hThetaDiff_NegTr->ProfileX();
    auto pfx_PhiDiff_NegTr =  hPhiDiff_NegTr->ProfileX();
    auto pfx_ThetaDiff_NegTrV = hThetaDiff_NegTrV->ProfileX();
    auto pfx_PhiDiff_NegTrV = hPhiDiff_NegTrV->ProfileX();
    
    pfx_ThetaDiff_NegTr->Rebin(8);
    pfx_PhiDiff_NegTr->Rebin(8);
    pfx_ThetaDiff_NegTrV->Rebin(8);
    pfx_PhiDiff_NegTrV->Rebin(8);

    pfx_ThetaDiff_NegTr->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV; #theta^{(rec)}_{#pi^{-}}; #theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}");
    pfx_PhiDiff_NegTr->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV; #phi^{(rec)}_{#pi^{-}}; #phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}");
    pfx_ThetaDiff_NegTrV->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV, common vertex reconstruction; #theta^{(rec)}_{#pi^{-}}; #theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}");
    pfx_PhiDiff_NegTrV->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV, common vertex reconstruction; #phi^{(rec)}_{#pi^{-}}; #phi^{(rec)}_{#pi^{-}} - #phi^{(gen)}_{#pi^{-}}");

    pfx_ThetaDiff_NegTr->GetYaxis()->SetRangeUser(-0.02, 0.02);
    pfx_PhiDiff_NegTr->GetYaxis()->SetRangeUser(-0.02, 0.02);
    pfx_ThetaDiff_NegTrV->GetYaxis()->SetRangeUser(-0.02, 0.02);
    pfx_PhiDiff_NegTrV->GetYaxis()->SetRangeUser(-0.02, 0.02);



    auto pfx_ThetaAvgDiff = hThetaAvgDiff->ProfileX();
    auto pfx_hPhiAvgDiff =  hPhiAvgDiff->ProfileX();
    auto pfx_hThetaAvgDiffV = hThetaAvgDiffV->ProfileX();
    auto pfx_hPhiAvgDiffV = hPhiAvgDiffV->ProfileX();

    pfx_ThetaAvgDiff->Rebin(8);
    pfx_hPhiAvgDiff->Rebin(8);
    pfx_hThetaAvgDiffV->Rebin(8);
    pfx_hPhiAvgDiffV->Rebin(8);

    pfx_ThetaAvgDiff->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV; #theta^{(rec)}_{avg}; #theta^{(rec)}_{avg} - #theta^{(gen)}_{avg}");
    pfx_hPhiAvgDiff->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV; #phi^{(rec)}_{avg}; #phi^{(rec)}_{avg} - #phi^{(gen)}_{avg}");
    pfx_hThetaAvgDiffV->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV common vertex reconstruction; #theta^{(rec)}_{avg}; #theta^{(rec)}_{avg} - #theta^{(gen)}_{avg}");
    pfx_hPhiAvgDiffV->SetTitle("#pi^{+}#pi^{-}, MCGPJ, E_{beam} = 509 MeV common vertex reconstruction; #phi^{(rec)}_{avg}; #phi^{(rec)}_{avg} - #phi^{(gen)}_{avg}");


    pfx_ThetaDiff_PosTr->GetYaxis()->SetTitleOffset(1.);
    pfx_PhiDiff_PosTr->GetYaxis()->SetTitleOffset(1.);
    pfx_ThetaDiff_PosTrV->GetYaxis()->SetTitleOffset(1.);
    pfx_PhiDiff_PosTrV->GetYaxis()->SetTitleOffset(1.);
    pfx_ThetaDiff_NegTr->GetYaxis()->SetTitleOffset(1.);
    pfx_PhiDiff_NegTr->GetYaxis()->SetTitleOffset(1.);
    pfx_ThetaDiff_NegTrV->GetYaxis()->SetTitleOffset(1.);
    pfx_PhiDiff_NegTrV->GetYaxis()->SetTitleOffset(1.);
    pfx_ThetaAvgDiff->GetYaxis()->SetTitleOffset(1.);
    pfx_hPhiAvgDiff->GetYaxis()->SetTitleOffset(1.);
    pfx_hThetaAvgDiffV->GetYaxis()->SetTitleOffset(1.);
    pfx_hPhiAvgDiffV->GetYaxis()->SetTitleOffset(1.);

    std::cout << "N_events: " << hThetaDiff_NegTr->GetEntries() << std::endl;
    file->Write();
    file->Save();
}
