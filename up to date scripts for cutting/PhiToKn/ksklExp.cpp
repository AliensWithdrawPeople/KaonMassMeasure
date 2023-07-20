#define ksklExp_cxx
#include "ksklExp.h"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVector3.h>

double PtotBinNumber(double mom)
{
    int nBin = 0;
    if(mom > 257)
    { nBin = 7; }
    else if (mom > 248.5)
    { nBin = 6; }
    else if (mom > 239)
    { nBin = 5; }
    else if (mom > 233)
    { nBin = 4; }
    else if (mom > 219)
    { nBin = 3; }
    else if (mom > 194.5)
    { nBin = 2; }
    else if (mom > 169.5)
    { nBin = 1; }
    return nBin;
}

void ksklExp::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      root> .L ksklExp.C
//      root> ksklExp t
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
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    TFile *top = new TFile(histFileName.c_str(), "recreate");
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    Float_t dpsi;
    Float_t ksTheta;
    Float_t ksPhi;
    Float_t Y;
    Float_t posMom, negMom;

    Float_t piThetaPos;
    Float_t piThetaNeg;

    Float_t piPhiPos;
    Float_t piPhiNeg;

    Float_t piMomPos;
    Float_t piMomNeg;

    tNew->Branch("emeas", &emeas0, "emeas/F");
    tNew->Branch("demeas", &demeas0, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", &dpsi, "dpsi/F");
    tNew->Branch("kstheta", &ksTheta, "kstheta/F");
    tNew->Branch("ksphi", &ksPhi, "ksphi/F");
    tNew->Branch("Y", &Y, "Y/F");

    tNew->Branch("piThetaPos", &piThetaPos, "piThetaPos/F");
    tNew->Branch("piThetaNeg", &piThetaNeg, "piThetaNeg/F");

    tNew->Branch("piPhiPos", &piPhiPos, "piPhiPos/F");
    tNew->Branch("piPhiNeg", &piPhiNeg, "piPhiNeg/F");

    tNew->Branch("piMomPos", &piMomPos, "piMomPos/F");
    tNew->Branch("piMomNeg", &piMomNeg, "piMomNeg/F");

    // pi+ pi- track reconstruction difference correction
    std::vector<double> Pavg = {157, 182, 207, 231, 235, 243, 254, 260};
    
    std::vector<std::vector<double>> dPtotPos{
        {0.119915, 0.0367675, 0.115981, 0.0150124, 0.0787089, -0.114334, 0.352348, 0.159914, 0.116781, -0.120692, 0.0468401, 0.113305, -0.112598, 0.0118255, -0.0153768, 0.107273, 0.107945, -0.0314952},
        {-0.0223037, 0.00114706, 0.196464, 0.0940865, 0.197963, 0.163502, 0.147584, 0.0980236, 0.290308, 0.281272, 0.213902, 0.0278824, 0.0536076, 0.0277654, 0.171705, 0.0652256, 0.188947, 0.070509},
        {0.173089, 0.0596282, 0.224753, 0.158979, 0.0770608, 0.22157, 0.163729, 0.0839866, 0.130302, 0.222514, 0.156945, 0.212422, 0.0910114, 0.161208, 0.337515, 0.0672538, 0.149937, 0.106422},
        {0.0549019, 0.105076, 0.321787, 0.108172, -0.00915358, 0.232584, 0.108496, 0.120291, 0.314675, 0.206091, 0.0668328, 0.218722, 0.143755, 0.161072, 0.277399, 0.104557, 0.0804615, 0.117139},
        {0.161264, 0.0623348, 0.422493, 0.274454, 0.0465108, 0.28003, 0.170985, 0.0986363, 0.310545, 0.20847, 0.156922, 0.188814, 0.232051, 0.177157, 0.295083, 0.155113, 0.0452571, 0.148004},
        {0.0635158, 0.2037, 0.446141, 0.0920958, 0.125514, 0.254835, 0.232825, 0.0101496, 0.317645, 0.281854, 0.360682, 0.161407, 0.249745, 0.311142, 0.303925, 0.0861456, 0.048336, 0.254039},
        {0.19121, 0.146413, 0.423835, 0.281012, 0.197206, 0.430928, 0.0705234, 0.212912, 0.451742, 0.137366, 0.25676, 0.206582, 0.278049, 0.242995, 0.477527, 0.20123, 0.254296, 0.226805},
        {0.313781, 0.119914, 0.522704, 0.502443, 0.141722, 0.449812, 0.368661, 0.159862, 0.644029, -0.0526354, 0.473181, 0.314891, 0.318761, 0.307549, 0.590559, 0.293004, 0.319808, 0.463133}
    };

    std::vector<std::vector<double>> dPtotNeg{
        {0.084882, 0.136052, -0.119107, 0.178941, -0.00785945, -0.021713, -0.0966209, 0.0251788, 0.0485461, 0.139593, 0.158504, -0.0383832, -0.0198412, -0.0900882, -0.191491, -0.0995606, 0.00327012, 0.0207491},
        {0.0432425, 0.0146009, -0.179285, 0.157402, 0.0607335, -0.0635736, 0.049901, 0.0903557, 0.0396047, 0.0975017, 0.120204, -0.00104125, 0.0641284, -0.0321337, 0.0570621, 0.0861944, 0.139065, 0.12177},
        {0.0728522, -0.0678266, 0.00923529, -0.0110787, 0.0770165, -0.00568652, -0.0746065, -0.028713, -0.0206241, 0.166844, 0.0727375, -0.0562728, -0.0263362, -0.12587, 0.0170237, -0.015987, 0.0324239, 0.00893105},
        {0.2979, -0.0906616, -0.0995954, 0.0931528, 0.0907167, -0.0342835, -0.0125813, 0.0551172, -0.0988879, 0.241792, 0.12736, -0.103945, 0.167163, 0.0427573, 0.0573838, -0.0314714, -0.0755831, -0.158469},
        {-0.00629191, -0.0698449, -0.134875, 0.00464847, 0.0793518, -0.0318978, -0.10969, -0.139721, -0.228363, 0.027629, 0.0464644, -0.159761, 0.0335226, 0.0666787, -0.0577746, -0.038914, -0.067291, -0.105659},
        {0.0434556, -0.0631131, -0.178364, -0.05619, 0.0212681, -0.0319327, -0.19772, -0.0821845, -0.111698, 0.00237257, 0.130398, -0.117076, 0.0331765, 0.169129, -0.084263, -0.135945, 0.00350779, -0.072382},
        {0.0242586, -0.0535996, -0.249788, 0.0284716, 0.0668667, -0.0919622, -0.108293, -0.11017, -0.371949, 0.0316157, 0.113404, -0.281787, -0.0309801, 0.148044, -0.0671834, -0.0814128, 0.0353451, -0.18092},
        {-0.0661348, 0.193507, -0.37129, -0.3405, -0.0807373, -0.213648, -0.204026, 0.00466489, -0.590979, 0.0867649, -0.179036, -0.237288, -0.0179135, -0.209035, -0.308377, -0.250435, -0.0544558, -0.263817}
    };

    std::vector<std::vector<double>> dPhiPos{
        {-0.00171841, -0.000713086, -0.00194934, -0.00075552, -0.000622023, 0.000402848, 0.000557336, -0.000497259, 1.04461e-06, -0.000114628, -0.0015742, 0.000219615, 6.55397e-05, -0.000503745, -0.000320512, -0.00148366, -0.00102118, 0.000244124},
        {-0.000264661, -5.97575e-05, -0.00189582, -0.000735752, -0.00112439, -0.00122268, -0.00177736, -0.000237151, -0.00219885, -0.00186689, -0.00105541, 0.000278961, -1.59819e-05, -0.00015473, -0.00112962, 7.87714e-05, -0.0013869, -0.000373716},
        {-0.000828786, -0.000190122, -0.00104632, -0.000443133, 0.000126639, -0.000447401, -0.000321892, -0.000546166, -0.000809326, -0.00092867, -0.00148553, -0.000634037, 0.000481242, 0.000386851, -0.00111419, -0.000285888, -0.000905787, -0.000874813},
        {-0.000312404, -0.000589547, -0.000951601, -0.000464232, -0.000176916, -0.00073326, -0.000296456, -8.5357e-05, -0.0011892, -0.000752748, -0.000346934, -0.000425162, 0.000173945, 0.000198806, -0.000905391, -0.000399599, -0.000453043, -0.000500869},
        {-0.000858511, -0.000369408, -0.00176218, -0.00113121, 3.53633e-05, -0.000774616, -0.000335391, -8.11299e-05, -0.00100314, -0.000641059, -0.000252993, -9.57728e-05, -0.000304755, -0.000101027, -0.00102666, -0.000607029, -0.000215155, -0.000633286},
        {-0.000642031, -0.000901121, -0.00185643, -0.000582681, -0.000312786, -0.000180728, 0.000175219, 0.000133407, -0.000949754, -0.00139783, -0.000528676, 0.000132331, -0.000835965, 9.59472e-05, -0.000522444, -0.000569269, -0.000342442, -0.000933475},
        {-0.000685818, -0.000575581, -0.00155276, -0.00100481, 0.000128186, -0.00124043, 0.000412472, -0.000264698, -0.00110254, -0.000277593, -0.000436719, 8.76531e-05, 9.07182e-05, -6.51789e-05, -0.00113332, -0.000514981, -0.000717286, -0.000874952},
        {-0.00108688, -0.000577309, -0.00115158, -0.00146797, -1.34959e-05, -0.000166657, -0.000260389, 1.39855e-05, -0.00117191, 6.83239e-05, -0.000572017, 1.63837e-05, 0.000480998, 0.000446467, -0.00129792, -0.00066396, -0.000904274, -0.00117598},
    };

    std::vector<std::vector<double>> dPhiNeg{
        {-0.000503283, 0.000509125, -0.000589746, 0.00190402, -0.000345595, -0.000738336, -0.000628465, 0.000903111, 0.000773863, 0.0019279, 0.000966168, -0.000239811, -0.000184621, -0.000226235, -0.000407359, -0.000336102, -0.00242461, -0.000723888},
        {0.000126477, 0.000277416, -0.00244242, 0.000725527, -0.00101011, -0.000506723, 0.000641548, 0.000671107, 0.000443427, 0.000868782, 0.000887614, -0.00034488, 0.000608085, -0.000878325, -0.000506825, 0.00021444, 0.00156808, 0.00188528},
        {-0.000184174, -0.000796152, -0.000778516, -0.000471322, 1.43354e-05, 0.00014908, -0.000659223, -0.000454644, -1.0506e-05, 0.0013642, 0.000462455, 0.000172649, -0.000635907, -0.00137824, -0.000289819, -7.06172e-05, -0.000567085, 0.00156351},
        {0.00115229, -0.00068634, -0.000889196, -3.89684e-05, -0.000217931, -0.000180165, 0.000248681, 0.00039728, -0.000247258, 0.000950679, 0.000707743, -0.000310605, 0.000713548, -2.31645e-05, -0.000127349, -3.94741e-05, -0.000187029, -0.000562813},
        {-0.000117577, -0.000782354, -0.000853202, -0.000159879, 2.12009e-05, -0.000156944, -0.000194829, -0.000419818, -0.000682926, 0.000180843, 0.000236089, -0.000598577, 3.05832e-05, -5.37666e-05, -0.000630708, -0.000276431, -0.000200059, -0.00017595},
        {-3.22916e-05, -0.000895992, -0.00172729, 0.000531201, -0.000671116, -0.000357941, -0.000380609, -0.000304787, -2.9849e-05, -7.2982e-05, 0.000608587, -0.000520108, 2.93506e-05, 0.000114721, -0.000429317, -0.00106561, -0.000268374, 0.000328358},
        {-0.000131757, -0.000374072, -0.00116368, -0.000321933, -0.000317392, 1.00203e-05, -0.000113744, -0.000295933, -0.000714535, 0.000150162, 0.000239532, -0.000633271, -2.38972e-05, -0.000366941, -6.70228e-05, -0.000186264, -0.000351156, -0.000547458},
        {-0.00072863, 2.69428e-05, -0.00173713, -0.0016593, -0.000602371, -0.000417179, -0.000221917, 0.000120655, -0.00113072, 7.74914e-05, -0.000172286, -0.000580373, -0.000205413, -0.000927578, -0.00136473, -0.00067639, -0.000363969, -0.000727758}
    };

    double dThetaPos = 0;
    double dThetaNeg = 0;

    Double_t halfPi = TMath::Pi() / 2;
    int NgoodTr = 0;
    int NgoodTrS = 0;
    double cutChi2r = 15.;
    double cutChi2z = 10.;
    int cutNhitMin = 6;
    int cutNhitMax = 30;
    double cutRmin = 0.05;
    double cutRmax = 6;
    double cutZtrack = 12.;
    double cutPtot = 40;
    double cutTrackTheta = 0.3;

    TVector3 ks;
    TVector3 kl;
    // dPhi and dTheta - phi and theta angle between Ks and Kl respectively
    double dPhi = 0;
    double dTheta = 0;
    int tmpCounter = 0;
    std::vector<Int_t> ksCand = {};

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);
    auto histKlCands = new TH1D("histKlCands", "number of Kl candidates for one Ks candidate", 5, 0, 5);
    auto histKsCands = new TH1D("histKsCands", "number of Ks candidates", 5, 0, 5);
    // Angles between tracks.
    auto hTrackColl = new TH2D("hTrackColl", "Angles between pions", 1200, -TMath::Pi(), TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    // Theta = kl.Theta()
    auto hdPhiTheta = new TH2D("hdPhiTheta", "", 600, 0, TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hdThetadPhi = new TH2D("hdThetadPhi", "", 1200, 0, 2*TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    auto hClEdPhi = new TH2D("hClEdPhi", "", 600, 0, 2 * TMath::Pi(), 600, 0, 600);
    auto hPsiUncutted = new TH1D("hPsiUncutted", "", 628, 0, 3.15);
    auto hPsiCutted = new TH1D("hPsiCutted", "", 628, 0, 3.15);
    // Phi angle between pions
    auto hKsKlPhi = new TH1D("hKsKlPhi", "", 1000, -1.5, 1.5);

    auto hKsKlTheta = new TH1D("hKsKlTheta", "", 1000, -1.5, 1.5);

    auto hHit = new TH1D("hHit", "nhits", 40, 0, 40);
    auto hTrackTheta = new TH1D("hTrackTheta", "hTrackTheta", 250, -TMath::Pi() / 2, TMath::Pi() / 2);
    auto hKstlenVsMinv = new TH2D("hKstlenVsMinv", "hKstlwnVsMinv", 1000, 420, 580, 1000, 0, 100);
    auto hTwoPhotonsTotEn = new TH1D("hTwoPhotonsTotEn", "hTwoPhotonsTotEn", 1000, 0, 1000);
    auto hGammasPi0Angles = new TH2D("hGammasPi0Angles", "hGammasPi0Angles", 2000, -TMath::Pi(), TMath::Pi(), 2000, -TMath::Pi(), TMath::Pi());
    auto hMissingMass = new TH1D("hMissingMass", "hMissingMass", 1000, 0, 1000);
    auto hMissingMom = new TH1D("hMissingMom", "hMissingMom", 700, 0, 700);

    hdThetadPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hdThetadPhi->GetYaxis()->SetTitle("#Delta#theta, rad");

    hTrackColl->GetXaxis()->SetTitle("|#phi_{1} - #phi_{2}| - #pi, rad");
    hTrackColl->GetYaxis()->SetTitle("#theta_{1} + #theta{2} - #pi, rad");

    hdPhiTheta->GetXaxis()->SetTitle("#theta of Kl, rad");
    hdPhiTheta->GetYaxis()->SetTitle("#Delta#phi, rad");

    hPsiUncutted->GetXaxis()->SetTitle("Angle between motion vectors Ks and Kl, rad");

    hClEdPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hClEdPhi->GetYaxis()->SetTitle("Cluster Energy, MeV");

    TVector3 piPos;
    TVector3 piNeg;
    TVector3 missingMom;
    TVector3 ph1Vec;
    TVector3 ph2Vec;
    TVector3 field(0., 0., 1.);

    double piPosEn = 0;
    double piNegEn = 0;
    double missingMass = 0;
    int posTrackNumber = 0;
    int negTrackNumber = 0;
    int phiIndexPos = 0;
    int phiIndexNeg = 0;
    const int Nchambers = 18;

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
	        if (tptot[i] > cutPtot && fabs(trho[i]) < cutRmax && fabs(tz[i]) < cutZtrack &&
                tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin && tnhit[i] < cutNhitMax)
            { NgoodTr++; }
        }
        //if(nt == 2)
        //{ hTrackColl->Fill(fabs(tphi[0] - tphi[1]) - TMath::Pi(), tth[0] + tth[1] - TMath::Pi()); }

        if (NgoodTr == 2)
        { NgoodTrS++; }
        
    
        if (NgoodTr == 2 && is_coll != 1 )
        {
            for(int k = 0; k < nks; k++)
            {
                if(ksalign[k] > 0.8 && (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000 &&
                1. < kspith[k][0] && kspith[k][0] < TMath::Pi() - 1. && 
                1. < kspith[k][1] && kspith[k][1] < TMath::Pi() - 1. &&
                //20 - half of the linear size of Drift Chamber
                //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
                // kspipt[k][0] > 120 && kspipt[k][1] > 120 && 
                // kspipt[k][0] < 350 && kspipt[k][1] < 350 &&
		        // (kspipt[k][0]+kspipt[k][1]) > 500 &&
		        kstlen[k] < 1.7 &&
                tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0  &&
                emeas0 > 100) // Added kstype[k] == 0.
                {
                    ks.SetMagThetaPhi(1, ksth[k], ksphi[k]);

                    for(int j = 0; j < nph; j++)
                    {
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
                        // phen0 - cluster energy. So it's a Kl candidate energy deposition.
                        // hClEdPhi->Fill(dPhi, phen0[j]); 
                        hdPhiTheta->Fill(kl.Theta(), dPhi);
                        
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

                        if((dPhi < -TMath::Pi() + 1 || dPhi > TMath::Pi() - 1) && fabs(dTheta) < 1 && phen0[j] > 40)
                        { 
                            tmpCounter++; 
                            hPsiCutted->Fill(ks.Angle(kl));
                            // hdThetadPhi->Fill(dPhi, dTheta);
                            if(dPhi < 0)
                            { hKsKlPhi->Fill(dPhi + TMath::Pi()); }
                            else
                            { hKsKlPhi->Fill(dPhi - TMath::Pi()); }
                            hKsKlTheta->Fill(dTheta);
                        } 
                    }        
        
                    histKlCands->Fill(tmpCounter);
                    if(tmpCounter > 0) // '|| 1' if there is no Kl cut
                    { ksCand.push_back(k); }
                    tmpCounter = 0;
                }
            }

            histKsCands->Fill(ksCand.size());
            if(ksCand.size() > 0)
            {
                posTrackNumber = tcharge[ksvind[ksCand[0]][0]] > 0 ? 0 : 1;
                negTrackNumber = posTrackNumber == 1 ? 0 : 1;

                phiIndexPos = int(kspiphi[ksCand[0]][posTrackNumber] / (2 * TMath::Pi() / Nchambers));
                phiIndexNeg = int(kspiphi[ksCand[0]][negTrackNumber] / (2 * TMath::Pi() / Nchambers));
                posMom = kspipt[ksCand[0]][posTrackNumber] + dPtotPos[PtotBinNumber(kspipt[ksCand[0]][posTrackNumber])][phiIndexPos];
                negMom = kspipt[ksCand[0]][negTrackNumber] + dPtotNeg[PtotBinNumber(kspipt[ksCand[0]][negTrackNumber])][phiIndexNeg];

                Y = posMom / negMom;
    
                piPos.SetMagThetaPhi(posMom, kspith[ksCand[0]][posTrackNumber], 
                                    kspiphi[ksCand[0]][posTrackNumber] + dPhiPos[PtotBinNumber(kspipt[ksCand[0]][posTrackNumber])][phiIndexPos]);
                piNeg.SetMagThetaPhi(negMom, kspith[ksCand[0]][negTrackNumber], 
                                        kspiphi[ksCand[0]][negTrackNumber] + dPhiNeg[PtotBinNumber(kspipt[ksCand[0]][negTrackNumber])][phiIndexNeg]);

                
                hKstlenVsMinv->Fill(ksminv[ksCand[0]], kstlen[ksCand[0]]); // New hist!!!!!!!!!!!!!!!!!!!!!!!!!
                
                dpsi = piPos.Angle(piNeg);
                ksTheta = ksth[ksCand[0]];
                ksPhi = ksphi[ksCand[0]];

                missingMom = -(piPos + piNeg);
                piPosEn = sqrt(139.57 * 139.57 + piPos.Mag2());
                piNegEn = sqrt(139.57 * 139.57 + piNeg.Mag2());
                missingMass = sqrt(4 * emeas0 * emeas0 + 2 * 139.57 * 139.57 
                                    - 2 * 2 * emeas0 * piPosEn - 2 * 2 * emeas0 * piNegEn
                                    + 2 * (piPosEn * piNegEn - piPos.Dot(piNeg)));
                hMissingMass->Fill(missingMass);

		        hTrackColl->Fill(fabs(tphi[0] - tphi[1]) - TMath::Pi(), tth[0] + tth[1] - TMath::Pi());
                if(missingMass > 350)
                { 
                    piThetaPos = piPos.Theta();
                    piThetaNeg = piNeg.Theta();

                    piPhiPos = piPos.Phi();
                    piPhiNeg = piNeg.Phi();

                    piMomPos = piPos.Mag();
                    piMomNeg = piNeg.Mag();

                    tNew->Fill(); 
                    hist->Fill(piPos.Mag(), piNeg.Mag()); 

                    // Cowboy
                    // if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) < TMath::Pi() / 2)
                    // { 
                    //     tNew->Fill(); 
                    //     hist->Fill(piPos.Mag(), piNeg.Mag());    
                    // }
                    
                    // // Sailor
                    // if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) > TMath::Pi() / 2)
                    // {  
                    //     tNew->Fill(); 
                    //     hist->Fill(piPos.Mag(), piNeg.Mag());    
                    // }
                }
            }
            ksCand.clear();
            ksCand.shrink_to_fit();
        }
        NgoodTr = 0;
    }

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
    auto dThetadPhiCutLine1 = new TLine(TMath::Pi() - 1, -0.3, TMath::Pi() + 1, -0.3);
    auto dThetadPhiCutLine2 = new TLine(TMath::Pi() - 1, -0.3, TMath::Pi() - 1, 0.3);
    auto dThetadPhiCutLine3 = new TLine(TMath::Pi() + 1, -0.3, TMath::Pi() + 1, 0.3);
    auto dThetadPhiCutLine4 = new TLine(TMath::Pi() - 1, 0.3, TMath::Pi() + 1, 0.3);
    dThetadPhiCutLine1->SetLineWidth(4);
    dThetadPhiCutLine2->SetLineWidth(4);
    dThetadPhiCutLine3->SetLineWidth(4);
    dThetadPhiCutLine4->SetLineWidth(4);
    // hdThetadPhi->Draw("col");
    // dThetadPhiCutLine1->Draw("same");
    // dThetadPhiCutLine2->Draw("same");
    // dThetadPhiCutLine3->Draw("same");
    // dThetadPhiCutLine4->Draw("same");
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

    std::cout << "nentries " << nentries << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    hist->Draw("COL");
    // hdThetadPhi->Draw("COL");
    // hPsiUncutted->Draw();
    // hPsiCutted->SetLineColor(kGreen);
    // hPsiCutted->Draw("same");


    // hKsKlTheta->Draw();

    top->Write();
    top->Save();
}