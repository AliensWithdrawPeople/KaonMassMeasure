// #include "PhiToK0_XSecBornFunc.hpp"
#include "TMath.h"
#include "TComplex.h"
#include "TF1.h"
#include "TSpline.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TFile.h"

Double_t RhoWidth(Double_t s);

Double_t OmgWidth(Double_t s);

Double_t PhiWidth(Double_t s, Double_t MPhi = 1019.461, Double_t WPhi = 4.249);

Double_t PKPKM(Double_t s);

Double_t PKLKS(Double_t s);

Double_t QPGamma(Int_t Mode, Double_t s);

Double_t PhitoK0(double x, const std::vector<double>& par);

Double_t eff1(double s);

Double_t PKPKM(Double_t s)
{
    Double_t MKC = 493.677;

    if ( s < 4.0*MKC*MKC ) return 0.0;

    return TMath::Power((s/4.0-MKC*MKC),1.5);
}

Double_t PKLKS(Double_t s)
{
    Double_t MK0 = 497.614;

    if ( s < 4.0*MK0*MK0 ) return 0.0;

    return TMath::Power((s/4.0-MK0*MK0),1.5);
}

Double_t QPGamma(Int_t Mode, Double_t s)
{
    Double_t PMass;

    if ( Mode == 0 ) PMass = 134.9768;
    if ( Mode == 1 ) PMass = 547.862;

    if ( s < PMass*PMass ) return 0.0;

    Double_t Fval =  TMath::Power(0.5*(1.0-PMass*PMass/s),3)*TMath::Power(s,1.5);

    return Fval;
}

Double_t FAS_ASPO(Double_t TwoE)
{
    Double_t Polinom(Int_t, Double_t, Double_t*);

    const Int_t K_max = 4;
    Int_t I_mark = 100;
    Int_t Mode = 1;

    bool Calc = true;

    Double_t A[K_max] = {-6.12394622, 25.0341405, -34.1311022, 15.5413717};
    Double_t B[K_max] = {5.29354148, -7.90990714, -2.26007613, 5.21453902};
    Double_t C[K_max] = {-0.56436115415016, 2.69953793, -4.32966739, 2.33116866};
    Double_t D[K_max] = {-0.0548334238,  0.31600391, -0.609523718, 0.393667808};

    Double_t A1[K_max] = {-4.91624401 , 19.8655606, -26.9136128 , 12.2412286};
    Double_t D1[K_max] = {-0.00794774189,0.0522269164,-0.114526409 , 0.0838126536};

    Double_t LM[6] = {1.1,0.875,0.75,0.62,0.52,0.46};

    Double_t POL, TEMP1, TEMP2; 

    if ( Calc ) {
        TEMP2 = D[1]*LM[4]+D[2]*LM[4]*LM[4]+D[3]*LM[4]*LM[4]*LM[4];
        TEMP1 = D1[0]+D1[1]*LM[4]+D1[2]*LM[4]*LM[4]+D1[3]*LM[4]*LM[4]*LM[4];
        D[0]  = TEMP1 - TEMP2;

        TEMP2 = C[1]*LM[3]+C[2]*LM[3]*LM[3]+C[3]*LM[3]*LM[3]*LM[3];
        TEMP1 = D[0]+D[1]*LM[3]+D[2]*LM[3]*LM[3]+D[3]*LM[3]*LM[3]*LM[3];
        C[0]  = TEMP1 - TEMP2;

        TEMP2 = A1[1]*LM[2]+A1[2]*LM[2]*LM[2]+A1[3]*LM[2]*LM[2]*LM[2];
        TEMP1 = C[0]+C[1]*LM[2]+C[2]*LM[2]*LM[2]+C[3]*LM[2]*LM[2]*LM[2];
        A1[0] = TEMP1 - TEMP2;

        TEMP1 = A1[0]+A1[1]*LM[1]+A1[2]*LM[1]*LM[1]+A1[3]*LM[1]*LM[1]*LM[1];
        TEMP2 = A[1]*LM[1]+A[2]*LM[1]*LM[1]+A[3]*LM[1]*LM[1]*LM[1];
        A[0]  = TEMP1 - TEMP2;
        
        TEMP1 = A[0] + A[1]*LM[0]+A[2]*LM[0]*LM[0]+A[3]*LM[0]*LM[0]*LM[0];
        TEMP2 = B[1]*LM[0]+B[2]*LM[0]*LM[0]+B[3]*LM[0]*LM[0]*LM[0];
        B[0]  = TEMP1 - TEMP2;

        Calc = false;
    }
    POL = 0.0;
    if ( TwoE >= LM[0] ) {
        for( Int_t i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*B[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*B[i];
        }
    }

    if ( TwoE >= LM[1] && TwoE < LM[0] ) {
        for( Int_t i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*A[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*A[i];
        }
    }

    if ( TwoE >= LM[2] && TwoE < LM[1] ) {
        for( Int_t i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*A1[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*A1[i];
        }
    }

    if ( TwoE >= LM[3] && TwoE < LM[2] ) {
        for( Int_t i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*C[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*C[i];
        }
    }

    if ( TwoE >= LM[4] && TwoE < LM[3] ) {
        for( Int_t i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*D[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*D[i];
        }
    }

    if ( TwoE >= LM[5] && TwoE < LM[4] ) {
        for( Int_t i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*D1[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*D1[i];
        }
    }
    if ( TwoE < LM[5] ) {
        if ( Mode == 1 ) POL = D1[0]+D1[1]*LM[5]+D1[2]*LM[5]*LM[5]+D1[3]*LM[5]*LM[5]*LM[5];
        if ( Mode == 2 ) POL = 3.0*LM[5]*LM[5]*D1[3]+2.0*LM[5]*D1[2]+D1[1]; 
    }


    Double_t Fval = (POL/0.393728*(0.00749/0.0361478));

    return Fval;
}

Double_t RhoWidth(Double_t s)
{
    Double_t MRho = 775.26;
    Double_t WRho = 149.10;

    Double_t BrKC = 0.491;
    Double_t BrK0 = 0.339;

    Double_t BrPi0Gamma = 0.00047;
    Double_t BrEtaGamma = 0.0003;

    Double_t Mpi = 139.57039;
    Double_t MPhi = 1019.461;

    // Double_t PhiWidth(Double_t);
    // Double_t PKPKM(Double_t);
    // Double_t PKLKS(Double_t);
    // Double_t QPGamma(Int_t, Double_t);

    Double_t Fval = WRho*MRho*MRho/s
                    * 0.9987386*TMath::Power((s/4.0-Mpi*Mpi)/(MRho*MRho/4.0-Mpi*Mpi),1.5)
        + 0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrKC*PKPKM(s)/PKPKM(MPhi*MPhi)/s
                    + 0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrK0*PKLKS(s)/PKLKS(MPhi*MPhi)/s
                    + WRho*BrPi0Gamma*QPGamma(0,s)/QPGamma(0,MRho*MRho)
                    + WRho*BrEtaGamma*QPGamma(1,s)/QPGamma(1,MRho*MRho);
    return Fval;
} 

Double_t OmgWidth(Double_t s)
{
    Double_t MOmg = 782.66;
    Double_t WOmg =   8.68;
    Double_t MPhi = 1019.461;
    Double_t WPhi = 4.249;

    Double_t MPi  = 139.57039;
    Double_t MPi0 = 134.9768;
    Double_t BrKC = 0.491;
    Double_t BrK0 = 0.339;
    Double_t BrEtaGamma = 0.00045;
    Double_t BrOmg3Pi   = 0.892;
    Double_t BrOmg2Pi   = 0.0153;
    Double_t BrOmgPiG   = 0.087;

    // Double_t PhiWidth(Double_t);
    // Double_t PKPKM(Double_t);
    // Double_t PKLKS(Double_t);
    // Double_t QPGamma(Int_t, Double_t);
    // Double_t FAS_ASPO(Double_t);

    Double_t Temp1, Temp2, Temp3, Temp4, Temp5, Temp6;
    Double_t X;

    Double_t TwoE;

    TwoE = TMath::Sqrt(s)/1000.0;

    if ( TMath::Sqrt(s) < 750.0 ) {
        Temp1 = 1.0E-6;
    } else {
        X = TMath::Sqrt(s)/1000.0 - 0.75;
        Temp1 = 1.0469 + 10.457*X + 53.009*X*X+591.66*X*X*X + 1487.7*X*X*X*X-7623.2*X*X*X*X*X+6949.5*X*X*X*X*X*X;
    }
    X = MOmg/1000.0 - 0.75;
    Temp2 = 1.0469 + 10.457*X + 53.009*X*X+591.66*X*X*X + 1487.7*X*X*X*X-7623.2*X*X*X*X*X+6949.5*X*X*X*X*X*X;

    Temp3 = TMath::Power(0.5*(1.0-MPi0*MPi0/s)*TMath::Sqrt(s),3.0);
    Temp4 = TMath::Power(0.5*(1.0-MPi0*MPi0/MOmg/MOmg)*TMath::Sqrt(MOmg*MOmg),3.0);

    Temp5 = TMath::Power((s/4.0-MPi*MPi),1.5);
    Temp6 = TMath::Power((MOmg*MOmg/4.0-MPi*MPi),1.5);

    Double_t Fval = WOmg*(BrOmg3Pi*FAS_ASPO(TwoE)/FAS_ASPO(MOmg/1000.0) +
        BrOmgPiG*QPGamma(0,s)/QPGamma(0,MOmg*MOmg)    +
        BrOmg2Pi*Temp5*MOmg*MOmg/(Temp6*s))           +
                            0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrKC*PKPKM(s)/PKPKM(MPhi*MPhi)/s +
                            0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrK0*PKLKS(s)/PKLKS(MPhi*MPhi)/s +
                            WOmg*BrEtaGamma*QPGamma(1,s)/QPGamma(1,MOmg*MOmg);

    return Fval;
}

Double_t PhiWidth(Double_t s, Double_t MPhi, Double_t WPhi)
{
    // Double_t MPhi = 1019.461;
    // Double_t WPhi = 4.249;
    Double_t MPi  = 139.57039;
    Double_t MPi2;

    Double_t Br[4] = {0.491, 0.339, 0.1532, 0.01309};
    Double_t Br2Pi = 0.000073;
    Double_t BrPi0Gamma = 0.00132;

    Double_t MPhi2 = 1019.461 * 1019.461;

    // Double_t PKPKM(Double_t);
    // Double_t PKLKS(Double_t);
    // Double_t QPGamma(Int_t, Double_t);

    // Double_t FAS_ASPO(Double_t);

    Double_t W1, W2;
    Double_t TwoE;

    TwoE = TMath::Sqrt(s)/1000.0;
    W1   = FAS_ASPO(TwoE);
    W2   = FAS_ASPO(MPhi/1000.0);

    MPi2 = MPi*MPi;


    Double_t Fval = MPhi2*WPhi*(Br[0]*PKPKM(s)/PKPKM(MPhi2) +
                                    Br[1]*PKLKS(s)/PKLKS(MPhi2) +
                    Br2Pi*TMath::Power((s/4.0-MPi2)/(MPhi2/4.0-MPi2),1.5))/s +
                                WPhi*(Br[2]*W1/W2+BrPi0Gamma*QPGamma(0,s)/QPGamma(0,MPhi2) +
                                    Br[3]*QPGamma(1,s)/QPGamma(1,MPhi2));

    return Fval;
}

Double_t PhitoK0(double *x, double *par)
{
    TComplex DPhi = TComplex(0.0,0.0);

    TComplex APhi = TComplex(0.0,0.0);
    TComplex ARho = TComplex(0.0,0.0);
    TComplex AOmg = TComplex(0.0,0.0);

    TComplex ATot = TComplex(0.0,0.0);

    Double_t RePart;
    Double_t ImPart;

    Double_t Smax, MPhi, WPhi;
    Double_t MPhi2;
    Double_t S;

    Double_t MRho = 775.26;
    Double_t MRho2;
    Double_t MOmg = 782.66;
    Double_t MOmg2;

    Double_t WRho = 149.10;
    Double_t WOmg =   8.68;

    Double_t BrRhoEE = 4.72E-5;
    Double_t BrOmgEE = 7.38E-5;

    Double_t BrPhiKK = 0.339;

    Double_t C = 0.389379292E12;

    S    = x[0] * x[0];

    Smax = par[0];
    MPhi = 1000.0 + par[1];
    WPhi = par[2];

    MPhi2 = MPhi*MPhi;
    MRho2 = MRho*MRho;
    MOmg2 = MOmg*MOmg;

    ARho = TMath::Sqrt(PhiWidth(MPhi2, MPhi, WPhi)*RhoWidth(MRho2)
        *	     MRho2*MRho*MPhi2*6.0*TMath::Pi()
        *             BrRhoEE*BrPhiKK / 2)
        / TComplex(MRho2-S,-TMath::Sqrt(S)*RhoWidth(S));

    AOmg = TMath::Sqrt(PhiWidth(MPhi2, MPhi, WPhi)*OmgWidth(MOmg2)
        *	     MOmg2*MOmg*MPhi2*6.0*TMath::Pi()
        *             BrOmgEE*BrPhiKK / 2)
        / TComplex(MOmg2-S,-TMath::Sqrt(S)*OmgWidth(S));

    // ARho = par[3] / TComplex(MRho2-S,-TMath::Sqrt(S)*RhoWidth(S));

    // AOmg = par[3] / TComplex(MOmg2-S,-TMath::Sqrt(S)*OmgWidth(S));

    // APhi = TMath::Sqrt(PhiWidth(MPhi2, MPhi, WPhi)*PhiWidth(MPhi2, MPhi, WPhi)
    //     *	     MPhi2*MPhi*MPhi2*6.0*TMath::Pi()
    //     *             par[0])
    //     / TComplex(MPhi2-S,-TMath::Sqrt(S)*PhiWidth(S, MPhi, WPhi)); 

    APhi = PhiWidth(MPhi2, MPhi, WPhi)*MPhi2*MPhi*TMath::Sqrt(Smax*MPhi/C)
        / TComplex(MPhi2-S,-TMath::Sqrt(S)*PhiWidth(S, MPhi, WPhi)); 

    ATot = - APhi - AOmg + ARho;

    RePart = ATot.Re();
    ImPart = ATot.Im();

    Double_t Fval = C*PKLKS(S)/PKLKS(MPhi2)/TMath::Power(S,2.5)
                    * (RePart*RePart+ImPart*ImPart);
    // if S < threshold then Sigma = 0
    double M_K0 = 497.614;
    // return (S < 4 * M_K0 * M_K0)? 0. : Fval;
    return Fval;
}

Double_t eff1(double s)
{
    std::vector<double> energy = {1001.859, 1005.96, 1009.6, 1015.724, 1016.808, 1017.914, 1019.056, 1019.912, 1020.916, 1022.07, 1022.896, 1027.728};
    std::vector<double> eff_data = {0.100742, 0.14724, 0.168875, 0.189286, 0.191459, 0.194718, 0.19831, 0.198671, 0.199915, 0.202736, 0.202236, 0.205726};
    auto efficiency = new TSpline3("efficiency", energy.data(), eff_data.data(), eff_data.size());
    return efficiency->Eval(TMath::Sqrt(s));
}

int PhiToK0_XSecBornFunc()
{
    std::vector<double> energies = {1004.066, 1010.466, 1012.955, 1015.068, 1016.105, 1017.155, 1017.156, 1018.046, 1019.118,  1019.214, 1019.421, 1019.902, 
                1021.222, 1021.309, 1022.078,  1022.744, 1023.264,  1025.320, 1027.956, 1029.090, 1033.907, 1040.028,  1049.864,  1050.862,  1059.947};
    std::vector<double> cross_sections = {6.87, 42.16, 96.74, 219.53, 366.33, 628.15, 624.76, 996.62, 1413.65,  1433.05, 1434.84, 1341.91,  833.20, 
                        807.54, 582.93,  443.71, 377.77,  199.26, 115.93, 96.96, 50.12, 31.27,  16.93,  17.47,  12.09};
    std::vector<double> xsec_err = {0.42, 0.47, 1., 5.02, 3.33, 2.95, 9.89, 4.28,  6.02, 15.03, 18.40, 4.74, 4.89, 10.36, 4.03, 4.38, 
                                    5.31, 4.97, 1.7, 3., 1.26, 1.01, 0.5, 0.94, 0.71};
    std::vector<double> zeroes(100, 0.0);

    // std::vector<double> energies = {1001.859, 1005.96, 1009.6, 1015.724, 1016.808, 1017.914, 1019.056, 1019.912, 1020.916, 1022.07, 1022.896, 1027.728};
    // std::vector<double> cross_sections = {5.6382, 11.0821, 17.8353, 200.678, 360.1538, 625.5156, 908.8247, 900.9489, 734.6828, 444.511, 347.8403, 135.666};
    // std::vector<double> xsec_err = {0.328, 0.248, 0.483, 1.441, 1.237, 2.017, 1.54, 1.735, 2.123, 1.821, 1.87, 1.156};
    // std::vector<double> zeroes(100, 0.0);
    
    TCanvas c("canv", "canv");
    TGraphErrors a(energies.size(), energies.data(), cross_sections.data(), zeroes.data(), xsec_err.data());
    auto bb = new TF1("born", PhitoK0, 1000., 1030., 3);

    bb->SetParameters(1431.08, 19.4572, 4.25065);
    // bb->SetParameters(10.04, 19., 4.26);

    // a.Fit(bb, "ME");
    // a.DrawClone("AP");

    TFile *file = TFile::Open("C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/bcs_BigChi2=935.root");
    auto gr = (TGraphErrors *)file->Get("bcs");
    gr->SetMarkerColor(kRed);
    gr->DrawClone("AP");
    // a.DrawClone("P same");
    bb->DrawClone("same");


    std::vector<double> energy = {1001.859, 1005.96, 1009.6, 1015.724, 1016.808, 1017.914, 1019.056, 1019.912, 1020.916, 1022.07, 1022.896, 1027.728};
    std::vector<double> eff_data = {0.100742, 0.14724, 0.168875, 0.189286, 0.191459, 0.194718, 0.19831, 0.198671, 0.199915, 0.202736, 0.202236, 0.205726};
    auto efficiency = new TSpline3("efficiency", energy.data(), eff_data.data(), eff_data.size());
    std::cout << efficiency->Eval(1000); 
    // efficiency->DrawClone();

    c.DrawClone();
    return 0;
}