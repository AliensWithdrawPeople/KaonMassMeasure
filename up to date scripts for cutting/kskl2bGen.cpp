#define kskl2bGen_cxx
#include "kskl2bGen.h"
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

void kskl2bGen::Loop(std::string histFileName)
{
    //   In a ROOT session, you can do:
    //      root> .L ksklCut_mcgpj.C
    //      root> ksklCut_mcgpj t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    // by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    TFile *top = new TFile(("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/" + histFileName).c_str(), "recreate");
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    Float_t dpsi;
    Float_t ksTheta;
    Float_t ksPhi;
    Float_t Y;
    Float_t eff;
    Float_t posMom, negMom;
    Float_t etrue;
    tNew->Branch("emeas", &emeas, "emeas/F");
    tNew->Branch("etrue", &etrue, "etrue/F");
    tNew->Branch("demeas", &demeas, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", &dpsi, "dpsi/F");
    tNew->Branch("kstheta", &ksTheta, "kstheta/F");
    tNew->Branch("ksphi", &ksPhi, "ksphi/F");
    tNew->Branch("Y", &Y, "Y/F");

    Double_t halfPi = TMath::Pi() / 2;
    int NgoodTr = 0;
    int NgoodTrS = 0;
    const double cutChi2r = 15.;
    const double cutChi2z = 10.;
    const int cutNhitMin = 6;
    const int cutNhitMax = 30;
    // const double cutRmin = 0.1;
    const double cutRmax = 6.0;
    const double cutZtrack = 12.;
    const double cutPtot = 40;
    // const double cutTrackTheta = 0.7;

    TVector3 ks;
    TVector3 kl;
    TVector3 ksGen;
    TVector3 klGen;
    // dPhi and dTheta - phi and theta angle between Ks and Kl respectively
    double dPhi = 0;
    double dTheta = 0;
    int KlCounter = 0;
    std::vector<Int_t> ksCand = {};

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);
    auto hSumMomTr = new TH1D("hSumMomTr", "Sum of moms", 4000, 0, 1200);
    auto histKlCands = new TH1D("histKlCands", "number of Kl candidates for one Ks candidate", 5, 0, 5);
    auto histKsCands = new TH1D("histKsCands", "number of Ks candidates", 5, 0, 5);
    // Theta = kl.Theta()
    auto hdPhiTheta = new TH2D("hdPhiTheta", "", 600, 0, TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hdThetadPhi = new TH2D("hdThetadPhi", "", 1200, 0, 2*TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    auto hdThetadPhiGen = new TH2D("hdThetadPhiGen", "", 1200, 0, 2*TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    auto hClEdPhi = new TH2D("hClEdPhi", "", 600, 0, 2 * TMath::Pi(), 600, 0, 600);
    auto hPsiUncutted = new TH1D("hPsiUncutted", "", 628, 0, 3.15);
    auto hPsiCutted = new TH1D("hPsiCutted", "", 628, 0, 3.15);
    auto hRho = new TH1D("hRho", "rho of vertex", 1000, 0, 10);
    // Phi angle between pions
    auto hPhi = new TH1D("hPhi", "", 1000, -0.2, 0.2);
    // Angles between tracks.
    auto hTrackColl = new TH2D("hTrackColl", "Angles between pions", 1200, -TMath::Pi(), TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    auto hTrackCollCutted = new TH2D("hTrackCollCutted", "Angles between pions cutted", 1200, -TMath::Pi(), TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    auto hDeltaMom = new TH1D("hDeltaMom", "Lorentz delta mom", 1000, -0.2, 0.2);
    auto hDeltaMomVsDeltaPhi1 = new TH2D("hDeltaMomVsDeltaPhi1", "delta mom vs delta phi", 1000, -0.1, 0.1, 1000, -10, 10);
    auto hDeltaMomVsDeltaPhi2 = new TH2D("hDeltaMomVsDeltaPhi2", "delta mom vs delta phi Cowboy", 1000, -1, 1, 1000, -10, 10);
    auto hDeltaMomVsDeltaPhi3 = new TH2D("hDeltaMomVsDeltaPhi3", "delta mom vs delta phi Sailor", 1000, -1, 1, 1000, -10, 10);
    auto hDeltaPhiVsPhi = new TH2D("hDeltaPhiVsPhi", "delta phi vs phi (rec - gen)", 1000, 0, 2 * 3.15, 1000, -1, 1);
    auto hDeltaPhiRecVsGenCowboy = new TH2D("hDeltaPhiRecVsGenCowboy", "delta phi Rec vs Gen data Cowboy", 1000, 0, 6.3, 1000, 0, 6.3);
    auto hDeltaPhiRecVsGenSailor = new TH2D("hDeltaPhiRecVsGenSailor", "delta phi Rec vs Gen data Sailor", 2000, 0, 6.3, 2000, 0, 6.3);

    auto hKsKlPhi = new TH1D("hKsKlPhi", "", 1000, -1.5, 1.5);
    auto hKsKlPhiGen = new TH1D("hKsKlPhiGen", "", 6000, -1.5, 1.5);

    auto hKsKlTheta = new TH1D("hKsKlTheta", "", 1000, -1.5, 1.5);
    auto hKsKlThetaGen = new TH1D("hKsKlThetaGen", "", 6000, -1.5, 1.5);

    auto hEVsdTheta = new TH2D("hEVsdTheta", "", 6000, 480, 540, 600, -1, 1);
    auto hEVsdThetaUncutted = new TH2D("hEVsdThetaUncutted", "", 6000, 480, 420, 600, -1, 1);
    auto hEVsdThetaGen = new TH2D("hEVsdThetaGen", "", 6000, 480, 540, 600, -1, 1);
    auto hEVsdThetaGenCut = new TH2D("hEVsdThetaGenCut", "", 6000, 480, 540, 600, -1, 1);

    auto hEUncutted = new TH1D("hEUncutted", "Uncutted Corrected energy", 5000, 480, 540);
    auto hE = new TH1D("hE", "Corrected energy", 6000, 480, 540);

    auto hHit = new TH1D("hHit", "nhits", 40, 0, 40);
    auto hTrackTheta = new TH1D("hTrackTheta", "hTrackTheta", 250, -TMath::Pi() / 2, TMath::Pi() / 2);

    auto hMissingMass = new TH1D("hMissingMass", "hMissingMass", 10000, 0, 10000);
    auto hKsInvMass = new TH1D("hKsInvMass", "hKsInvMass", 5000, 0, 1000);
    auto hFinalStateId = new TH1D("hFinalStateId", "hFinalStateId", 1000, 0, 1000);

    hdThetadPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hdThetadPhi->GetYaxis()->SetTitle("#Delta#theta, rad");

    hTrackColl->GetXaxis()->SetTitle("|#phi_{1} - #phi_{2}| - #pi, rad");
    hTrackColl->GetYaxis()->SetTitle("#theta_{1} + #theta_{2} - #pi, rad");

    hTrackCollCutted->GetXaxis()->SetTitle("|#phi_{1} - #phi_{2}| - #pi, rad");
    hTrackCollCutted->GetYaxis()->SetTitle("#theta_{1} + #theta_{2} - #pi, rad");

    hdPhiTheta->GetXaxis()->SetTitle("#theta of Kl, rad");
    hdPhiTheta->GetYaxis()->SetTitle("#Delta#phi, rad");

    hPsiUncutted->GetXaxis()->SetTitle("Angle between motion vectors Ks and Kl, rad");

    hClEdPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hClEdPhi->GetYaxis()->SetTitle("Cluster Energy,  MeV");

    TVector3 piPos;
    TVector3 piNeg;
    TVector3 piPosRec;
    TVector3 piNegRec;
    TVector3 field(0., 0., 1.);

    auto dPhiVsEcutFunc = new TF1("dPhiVsEcutFunc", "[1] / ([0] - x) + [2]", 490, 520);
    double missingMass = 0;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        for (int i = 0; i < nt; i++)
        {
            hRho->Fill(fabs(trho[i]));
            if (tptot[i] > cutPtot && fabs(trho[i]) < cutRmax && fabs(tz[i]) < cutZtrack &&
                tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin && tnhit[i] < cutNhitMax)
            { NgoodTr++; }
        }

        if(nt == 2)
        { hTrackColl->Fill(fabs(tphi[0] - tphi[1]) - TMath::Pi(), tth[0] + tth[1] - TMath::Pi()); }

        if (NgoodTr == 2)
        { NgoodTrS++; }
        
    
        if (NgoodTr == 2 && is_coll != 1)
        {
            for(int k = 0; k < nks; k++)
            {
                if(ksalign[k] > 0.85 && (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000 &&
                // abs(kspith[k][0] - TMath::Pi() / 2) <= cutTrackTheta && 
                // abs(kspith[k][1] - TMath::Pi() / 2) <= cutTrackTheta &&
                1.1 < kspith[k][0] && kspith[k][0] < TMath::Pi() - 1.1 && 
                1.1 < kspith[k][1] && kspith[k][1] < TMath::Pi() - 1.1 &&
                //20 - half of the linear size of Drift Chamber
                //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
                // kspipt[k][0] > 120 && kspipt[k][1] > 120 && 
                // kspipt[k][0] < 350 && kspipt[k][1] < 350 &&
                tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0) // Added kstype[k] == 0.
                {
                    ks.SetMagThetaPhi(1, ksth[k], ksphi[k]);

                    for(int j = 0; j < nph; j++)
                    {
                        // if(phen0[j] < 40)
                        // { continue; }
                        // phth0 and phphi0 are theta and phi respectively angles of KL candidate.
                        kl.SetMagThetaPhi(phrho[j], phth0[j], phphi0[j]);
                        // Shift to the center of phi decay. 
                        kl.SetX(kl.X() - xbeam);
                        kl.SetY(kl.Y() - ybeam);
                        kl.SetZ(kl.Z() - ksz0[k]);

                        hPsiUncutted->Fill(ks.Angle(kl));

                        dPhi = ks.DeltaPhi(kl);
                        dTheta = ks.Theta() + kl.Theta() - TMath::Pi();
                        // hdThetadPhi->Fill(dPhi, dTheta);
                        // phen0 - cluster energy. So this is the energy deposition of Kl candidate.
                        // hClEdPhi->Fill(dPhi, phen0[j]); 
                        hdPhiTheta->Fill(kl.Theta(), dPhi);

                        for(int i = 0; i < nsim; i++)
                        {
                            if(simtype[i] == 211 && simorig[i] == 310)
                            { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                            if(simtype[i] == -211 && simorig[i] == 310)
                            { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }

                            if(simtype[i] == 310)
                            { ksGen.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                            if(simtype[i] == 130)
                            { klGen.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                        }
                        // hEVsdTheta->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()), dTheta);
                        
                        // [1] / ([0] - x) + [2]
                        dPhiVsEcutFunc->SetParameters(emeas - 0.5, 0.1, 0.1);
                        if(dPhi < 0)
                        { 
                            hdThetadPhi->Fill(dPhi + 2*TMath::Pi(), dTheta); 
                            hClEdPhi->Fill(dPhi + 2*TMath::Pi(), phen0[j]);
                        }
                        else
                        { 
                            hdThetadPhi->Fill(dPhi, dTheta); 
                            hClEdPhi->Fill(dPhi, phen0[j]);
                        }

                        etrue = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2());
                        hEVsdThetaUncutted->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()), dTheta);
                        // if((dPhi < -TMath::Pi() + 1 || dPhi > TMath::Pi() - 1) && (fabs(dTheta) - dPhiVsEcutFunc->Eval(etrue) < 0 || etrue > emeas - 0.5) && phen0[j] > 40)
                        if((dPhi < -TMath::Pi() + 1 || dPhi > TMath::Pi() - 1) && fabs(dTheta) < 1 && phen0[j] > 40)
                        { 
                            KlCounter++; 
                            hPsiCutted->Fill(ks.Angle(kl));
                            // hdThetadPhi->Fill(dPhi, dTheta);
                            if(dPhi < 0)
                            { hKsKlPhi->Fill(dPhi+ TMath::Pi()); }
                            else
                            { hKsKlPhi->Fill(dPhi - TMath::Pi()); }
                            hKsKlTheta->Fill(dTheta);
                            hEVsdThetaGenCut->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()), ksGen.Theta() + klGen.Theta() - TMath::Pi());
                            hEVsdTheta->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()), dTheta);
                        } 
                    }        
        
                    histKlCands->Fill(KlCounter);
                    if(KlCounter > 0) // '|| 1' if there is no Kl cut
                    { ksCand.push_back(k); }
                    KlCounter = 0;
                }
            }

            histKsCands->Fill(ksCand.size());
            if(ksCand.size() > 0)
            {
                for(int i = 0; i < nsim; i++)
                {
                    if(simtype[i] == 211 && simorig[i] == 310)
                    { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }

                    if(simtype[i] == -211 && simorig[i] == 310)
                    { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                }

                int posTrackNumber = tcharge[ksvind[ksCand[0]][0]] > 0 ? 0 : 1;
                int negTrackNumber = posTrackNumber == 1 ? 0 : 1;
            
                Y = kspipt[ksCand[0]][posTrackNumber] / kspipt[ksCand[0]][negTrackNumber];
                hist->Fill(kspipt[ksCand[0]][posTrackNumber], kspipt[ksCand[0]][negTrackNumber]);
                posMom = kspipt[ksCand[0]][posTrackNumber]; 
                negMom = kspipt[ksCand[0]][negTrackNumber];
                piPosRec.SetMagThetaPhi(kspipt[ksCand[0]][posTrackNumber], kspith[ksCand[0]][posTrackNumber], kspiphi[ksCand[0]][posTrackNumber]);
                piNegRec.SetMagThetaPhi(kspipt[ksCand[0]][negTrackNumber], kspith[ksCand[0]][negTrackNumber], kspiphi[ksCand[0]][negTrackNumber]);

                hSumMomTr->Fill(kspipt[ksCand[0]][negTrackNumber] + kspipt[ksCand[0]][posTrackNumber]);
                dpsi = ksdpsi[ksCand[0]];
                ksTheta = ksth[ksCand[0]];
                ksPhi = ksphi[ksCand[0]];
                etrue = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2());            
                
                missingMass = sqrt(4 * emeas * emeas + 2 * 139.57 * 139.57 
                - 2 * 2 * emeas * sqrt(139.57 * 139.57 + piPosRec.Mag2()) - 2 * 2 * emeas * sqrt(139.57 * 139.57 + piNegRec.Mag2()) 
                + 2 * (sqrt(139.57 * 139.57 + piPosRec.Mag2()) * sqrt(139.57 * 139.57 + piNegRec.Mag2()) - piPosRec.Dot(piNegRec)));

                double pipiInvMass = sqrt(2 * (139.57 * 139.57 + sqrt(139.57 * 139.57 + piPosRec.Mag2()) * sqrt(139.57 * 139.57 + piNegRec.Mag2()) - piPosRec.Dot(piNegRec)));
                hMissingMass->Fill(missingMass);
                hKsInvMass->Fill(pipiInvMass);
                hFinalStateId->Fill(finalstate_id);


                // Data from the generator
                // Y = piPos.Mag() / piNeg.Mag();
                // dpsi = piPos.Angle(piNeg);
                // emeas = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2());

                if(missingMass > 350)
                { 
                    tNew->Fill(); 
                    hTrackCollCutted->Fill(fabs(tphi[0] - tphi[1]) - TMath::Pi(), tth[0] + tth[1] - TMath::Pi());
                    hE->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));

                    //  Cowboy Type
                    if(piPosRec.Cross(field).XYvector().DeltaPhi(piNegRec.XYvector()) < TMath::Pi() / 2)
                    { 
                        // tNew->Fill(); 
                        // hE->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));
                        hDeltaPhiVsPhi->Fill(piPos.XYvector().Phi(), piPos.XYvector().DeltaPhi(piPosRec.XYvector()));
                        hDeltaPhiRecVsGenCowboy->Fill(fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())), 
                                                    fabs(piPosRec.XYvector().DeltaPhi(piNegRec.XYvector())) );

                        // hPhi->Fill(piPos.Phi() - piPosRec.Phi()); 
                        // hPhi->Fill(piNeg.Phi() - piNegRec.Phi()); 

                        hDeltaMomVsDeltaPhi2->Fill( (piPos.Mag() - piNeg.Mag()) - (piPosRec.Mag() - piNegRec.Mag()),
                                                    piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()));
                    }
                    //  Sailor Type
                    if(piPosRec.Cross(field).XYvector().DeltaPhi(piNegRec.XYvector()) > TMath::Pi() / 2)
                    {  
                        // tNew->Fill();
                        // hE->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));
                        // hDeltaPhiVsPhi->Fill(piPos.XYvector().Phi(), piPos.XYvector().DeltaPhi(piPosRec.XYvector()));

                        hDeltaPhiRecVsGenSailor->Fill(fabs(piPos.XYvector().DeltaPhi(piNeg.XYvector())), 
                                                    fabs(piPosRec.XYvector().DeltaPhi(piNegRec.XYvector())) );
                        hDeltaMomVsDeltaPhi3->Fill(piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()), 
                                                (piPos.Mag() - piNeg.Mag()) - (piPosRec.Mag() - piNegRec.Mag())); 
                    }
                }

                // hDeltaMomVsDeltaPhi2->Fill(piPos.Phi() - piPosRec.Phi(), piPos.Mag() - piPosRec.Mag());
                // hDeltaMomVsDeltaPhi3->Fill(piNeg.Phi() - piNegRec.Phi(), piNeg.Mag() - piNegRec.Mag());
                hDeltaMomVsDeltaPhi1->Fill(piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()), 
                                            (piPos.Mag() - piNeg.Mag()) - (piPosRec.Mag() - piNegRec.Mag()));

                hDeltaMom->Fill( (piPos.Mag() - piPosRec.Mag()) / piPos.Mag());
            }
            ksCand.clear();
            ksCand.shrink_to_fit();

            for(int i = 0; i < nsim; i++)
            {
                if(simtype[i] == 211 && simorig[i] == 310)
                { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                if(simtype[i] == -211 && simorig[i] == 310)
                { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
            }
            hEUncutted->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));
        }
        NgoodTr = 0;

        for(int i = 0; i < nsim; i++)
        {
            if(simtype[i] == 211 && simorig[i] == 310)
            { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
            if(simtype[i] == -211 && simorig[i] == 310)
            { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }

            if(simtype[i] == 310)
            { ksGen.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
            if(simtype[i] == 130)
            { klGen.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
        }
        // Y = piPos.Mag() / piNeg.Mag();
        // dpsi = piPos.Angle(piNeg);
        // emeas = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2());
        // etrue = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2());
        // tNew->Fill();

        if(klGen.DeltaPhi(ksGen) < 0)
        { hKsKlPhiGen->Fill(klGen.DeltaPhi(ksGen) + TMath::Pi()); }
        else
        { hKsKlPhiGen->Fill(klGen.DeltaPhi(ksGen) - TMath::Pi()); }
        hKsKlThetaGen->Fill(ksGen.Theta() + klGen.Theta() - TMath::Pi());
        hEVsdThetaGen->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()), ksGen.Theta() + klGen.Theta() - TMath::Pi());
        hdThetadPhiGen->Fill(klGen.DeltaPhi(ksGen) + TMath::Pi(), ksGen.Theta() + klGen.Theta() - TMath::Pi());
    }

    eff = (double) tNew->GetEntriesFast() / nentries;
    std::cout << "nentries " << nentries << std::endl;
    std::cout << "e_mc = " <<  eff << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;

    // Drawing hits with cuts
    fChain->Draw("tnhit>>hHit", "", "goff");

    auto hitCutLine1 = new TLine(10, 0, 10, 10000);
    hitCutLine1->SetLineColor(kBlue);
    hitCutLine1->SetLineWidth(4);
    auto hitCutLine2 = new TLine(30, 0, 30, 10000);
    hitCutLine2->SetLineColor(kBlue);
    hitCutLine2->SetLineWidth(4);

    hHit->GetXaxis()->SetTitle("Number of hits");
    // hHit->Draw();
    // hitCutLine1->Draw("same");
    // hitCutLine2->Draw("same");
    // Drawing hits with cuts 

    // ****************************************************

    // Drawing track theta with cuts' lines 
    fChain->Draw("tth - TMath::Pi() / 2>>hTrackTheta", "", "goff");

    auto thetaCutLine1 = new TLine(-0.7, 0, -0.7, 1600);
    thetaCutLine1->SetLineColor(kBlue);
    thetaCutLine1->SetLineWidth(4);
    auto thetaCutLine2 = new TLine(0.7, 0, 0.7, 1600);
    thetaCutLine2->SetLineColor(kBlue);
    thetaCutLine2->SetLineWidth(4);

    hTrackTheta->GetXaxis()->SetTitle("#theta - #frac{#pi}{2}, rad");
    // hTrackTheta->Draw();
    // thetaCutLine1->Draw("same");
    // thetaCutLine2->Draw("same");
    // Drawing track theta with cuts' lines 

    // ****************************************************

    // Drawing KsKl dThetadPhi with cuts' lines
    auto dThetadPhiCutLine1 = new TLine(TMath::Pi() - 1, -1, TMath::Pi() + 1, -1);
    auto dThetadPhiCutLine2 = new TLine(TMath::Pi() - 1, -1, TMath::Pi() - 1, 1);
    auto dThetadPhiCutLine3 = new TLine(TMath::Pi() + 1, -1, TMath::Pi() + 1, 1);
    auto dThetadPhiCutLine4 = new TLine(TMath::Pi() - 1, 1, TMath::Pi() + 1, 1);
    dThetadPhiCutLine1->SetLineWidth(4);
    dThetadPhiCutLine2->SetLineWidth(4);
    dThetadPhiCutLine3->SetLineWidth(4);
    dThetadPhiCutLine4->SetLineWidth(4);
    hdThetadPhi->Draw("col");
    dThetadPhiCutLine1->Draw("same");
    dThetadPhiCutLine2->Draw("same");
    dThetadPhiCutLine3->Draw("same");
    dThetadPhiCutLine4->Draw("same");
    // Drawing KsKl dThetadPhi with cuts' lines

    // ****************************************************

    // Drawing KsKl Cluster Energy vs dPhi with cuts' lines
    auto hClEdPhiCutLine1 = new TLine(TMath::Pi() - 1, 40, TMath::Pi() + 1, 40);
    auto hClEdPhiCutLine2 = new TLine(TMath::Pi() - 1, 40, TMath::Pi() - 1, 600);
    auto hClEdPhiCutLine3 = new TLine(TMath::Pi() + 1, 40, TMath::Pi() + 1, 600);
    hClEdPhiCutLine1->SetLineWidth(4);
    hClEdPhiCutLine2->SetLineWidth(4);
    hClEdPhiCutLine3->SetLineWidth(4);
    // hClEdPhi->Draw("col");
    // hClEdPhiCutLine1->Draw("same");
    // hClEdPhiCutLine2->Draw("same");
    // hClEdPhiCutLine3->Draw("same");
    // Drawing KsKl Cluster Energy vs dPhi with cuts' lines

    std::cout << "Energy point: " << emeas << "; Mean Energy = " << hE->GetMean() << " +/- " << hE->GetMeanError() << "; e_mc = " <<  eff << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    // hist->Draw("COL");
    // hdThetadPhi->Draw("COL");
    // hPsiUncutted->Draw();
    // hPsiCutted->SetLineColor(kGreen);
    // hPsiCutted->Draw("same");

    // auto tmpForm = new TF1("tmpForm", "[1] / ([0] - x) + [2]", 490, 520);
    // tmpForm->SetParameters(emeas - 0.5, 0.1, 0.1);
    // hEVsdTheta->Draw("COL");
    // tmpForm->Draw("same");

    // hKsKlThetaGen->SetLineColor(kGreen);
    // hKsKlPhiGen->DrawNormalized("", hKsKlPhi->Integral());
    // hKsKlPhiGen->Draw();
    // hKsKlTheta->Draw();


    // hDeltaMomVsDeltaPhi2->GetYaxis()->SetTitle("#DeltaP_{#pi^{+}} [MeV/c]");
    // hDeltaMomVsDeltaPhi1->Draw("COL");
    // auto c = new TCanvas("deltaP vs deltaPhi", "#pi^{-}", 200, 10, 600, 400);
    // c->Divide(2, 1);
    // c->cd(1);
    // hDeltaMomVsDeltaPhi2->Draw("COL");
    // c->cd(2);
    // hDeltaMomVsDeltaPhi3->Draw("COL");


    top->Write();
    top->Save();
}