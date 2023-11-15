/**
 * @author Parametrization and all values of parameters are taken from E.A.Kozyrev's article (https://arxiv.org/abs/1604.02981).
 * Code was written by E.A.Kozyrev (e.a.kozyrev@inp.nsk.su) and was updated in terms of constants (pdg live, July 2023) 
 * and readability by Daniil Ivanov (daniilivanov1606@gmail.com).
*/

#include <memory>
#include <cmath>
#include <complex>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>

#define mK 493.677
#define mKn 497.614
#define me 0.511
#define C 0.389379292E12 // 1/Mev^2 = C [barn]; Used in conversion from natural units to barns.
#define alphafss 0.0072973525628
#define pi 3.1415926535


/**
 * // Cross section: KnCrossSection(double *Ecms, double *par);
 * // Form factor: KnFormFactor(double *Ecms, double *par);
 * double M_phi = 1019.457; \\ par[0], MeV -- phi meson mass.
 * double W_phi = 4.240; \\ par[1], MeV -- phi meson width.
 * double W_phi_ee_B_phi_KsKl = 0.428e-03 // par[2], MeV -- \Gamma_{\phi->e^+e^-} * B_{\phi->KsKl} = 0.428 keV = 0.428e-03 MeV. 
 * double g_rho_omega = 0.8; \\ par[3] -- g_{\rho, \omega} multiplicative factor for coupling constants  g_{\rho KsKl}, g_{\omega KsKl} (read more in the article).
 * std::vector<double> params = {1019.457, 4.240, 0.428e-03, 0.8};
 * // pass params as params.data().
*/

Double_t PointLikeKnCrossSection(Double_t s);
Double_t RhoWidth(Double_t s);
Double_t OmgWidth(Double_t s);
Double_t PhiWidth(Double_t s);
Double_t PKPKM(Double_t s);
Double_t PKLKS(Double_t s);
Double_t QPGamma(Int_t Mode, Double_t s);
Double_t FAS_ASPO(Double_t TwoE);

/** @param charge For Neutral kaons == -1. 
*/ 
std::complex<double> KFF_total(double Ecms, const std::vector<double> par, double charge)
{
    Double_t s = Ecms*Ecms;

    Double_t MPhi = par[0];
    Double_t WPhi = par[1];
    Double_t gb = par[2]; // Gamma * Branching
        
    Double_t Momega = 782.66;
    Double_t Mrho = 775.26;
    double Wrho = 147.4;
    double Womega = 8.68;

    Double_t Momega1420 = 1410.;
    double Womega1420 = 290.;
    double GBomega1420 = 2.20257e-07;
    
    double Maphi1680 = par[3]; 
    double Wphi1680 = par[4];
    double GBphi1680 = par[5];
    
    double Mrho1450 = par[6];
    double Wrho1450 = par[7];
    double GBrho1450 = par[8];
    
    double Bphi_ee = 2.979e-04;
    double Bphi_KK = gb/Bphi_ee/WPhi;
    double Brho_ee = 0.0000472;
    double Bomega_ee = 0.0000738;
    
    //************************phi width*********************************************
    Double_t Br[4] = {0.491, Bphi_KK, 0.1532, 0.01309};
    if(charge == 1){ Br[0] = Bphi_KK; Br[1] = 0.339;}
    
    Double_t MPhi2 = MPhi*MPhi;
    
    Double_t PKPKM(Double_t);
    Double_t PKLKS(Double_t);
    Double_t QPGamma(Int_t, Double_t);
    
    Double_t FAS_ASPO(Double_t);
    
    Double_t W1, W2;
    Double_t TwoE;
    
    TwoE = TMath::Sqrt(s)/1000.0;
    W1   = FAS_ASPO(TwoE);
    W2   = FAS_ASPO(MPhi/1000.0);
    
    std::complex<double> ATot = std::complex<double>(0.0,0.0);
    Double_t BrPi0Gamma = 0.00132;
    Double_t WPhi_s =   WPhi*MPhi*MPhi/s*(
                            Br[0]*PKPKM(s)/PKPKM(MPhi2) +
                            Br[1]*PKLKS(s)/PKLKS(MPhi2)) + 
                        WPhi*(Br[2]*W1/W2+BrPi0Gamma*QPGamma(0,s)/QPGamma(0,MPhi2) +
                        Br[3]*QPGamma(1,s)/QPGamma(1,MPhi2));
                        //WPhi*(Br[2]*W1/W2 + 1. - Br[0] - Br[1] - Br[2]);
    //******************************************************************************

    std::complex<double> aphi = WPhi*pow(MPhi,3.)*sqrt(12.*TMath::Pi()*gb/MPhi/WPhi)/std::complex<double>(MPhi*MPhi-s,-sqrt(s)*WPhi_s);  
    
    std::complex<double> arho = sqrt(WPhi*Wrho*pow(Mrho,3.)*pow(MPhi,2.)*6.*TMath::Pi()*Bphi_KK*Brho_ee)/std::complex<double>(Mrho*Mrho-s,-sqrt(s)*RhoWidth(s));
    std::complex<double> arho1450 = sqrt(12.*TMath::Pi()*pow(Mrho1450,5.)*Wrho1450*GBrho1450)/std::complex<double>(Mrho1450*Mrho1450-s, -sqrt(s)*Wrho1450);
    
    std::complex<double> aomega = sqrt(WPhi*Womega*pow(Momega,3.)*pow(MPhi,2.)*6.*TMath::Pi()*Bphi_KK*Bomega_ee)/std::complex<double>(Momega*Momega-s,-sqrt(s)*OmgWidth(s));
    std::complex<double> aomega1420 = sqrt(12.*TMath::Pi()*pow(Momega1420,5.)*Womega1420*GBomega1420)/std::complex<double>(Momega1420*Momega1420-s, -sqrt(s)*Womega1420);

    std::complex<double> aphi1680 = sqrt(12.*TMath::Pi()*pow(Maphi1680,5.)*Wphi1680*GBphi1680)/std::complex<double>(Maphi1680*Maphi1680-s, -sqrt(s)*Wphi1680);
    
    return aphi + par[9]*(aomega + charge*arho)*std::complex<double>(cos(par[10]),sin(par[10])) - par[11]*(charge*arho1450 + aphi1680) - par[12]*std::complex<double>(cos(par[13]),sin(par[13]));   
}

/** @brief Set mode of KFF_total: -1 -- K0
 * @param Ecms = 2 * E_beam MeV
 * @param par Vector of parameters that are used in phi -> KsKl cross section parametrization.
*/
std::complex<double> KnFF_total(double Ecms, const std::vector<double> par)
{ return KFF_total(Ecms,par,-1.); }

/**
 * @note pass params as params.data().
 * @param Ecms = 2 * E_beam MeV
 * @param par Vector of parameters that are used in phi -> KsKl cross section parametrization.
 * @returns Cross section in nb.
*/
double KnCrossSection(double *Ecms, double *par)
{
    double s = Ecms[0]*Ecms[0];
    std::vector<double> par1;
    par1.push_back(par[0]);
    par1.push_back(par[1]);
    par1.push_back(par[2]);
    par1.push_back(1680.69);
    par1.push_back(150.586);
    par1.push_back(1.33238e-07);
    par1.push_back(1503.66);
    par1.push_back(387.094);
    par1.push_back(2.68137e-09);
    par1.push_back(par[3]);
    par1.push_back(0.0);
    par1.push_back(1.);
    par1.push_back(0.0);
    par1.push_back(0.0);
    
    std::complex<double> ATot = KnFF_total(Ecms[0], par1);
    double MPhi = par[0];
    double Fval = (std::norm(ATot))/pow(s,2.5)*C*pow(s/4.-mKn*mKn,1.5)/pow(MPhi*MPhi/4.-mKn*mKn,1.5); 

    // if S < threshold then Sigma = 0
    double M_K0 = 497.614;
    return (s < 4 * M_K0 * M_K0)? 0. : Fval;
}

/** @brief Evaluate form-factor of neutral kaon F(s), which is defined as 
 * sigma_{e+e- -> KsKl}(s) = sigma^{point-like}_{e+e- -> KsKl}(s) * |F(s)|^2.
 * @note pass params as params.data().  
 * @param Ecms = 2 * E_beam MeV
 * @param par Vector of parameters that are used in phi -> KsKl cross section parametrization.
 * @returns Form factor of neutral kaon.
*/
std::complex<double> KnFormFactor(double *Ecms, double *par)
{
    double s = Ecms[0]*Ecms[0];
    std::vector<double> par1;
    par1.push_back(par[0]);
    par1.push_back(par[1]);
    par1.push_back(par[2]);
    par1.push_back(1680.69);
    par1.push_back(150.586);
    par1.push_back(1.33238e-07);
    par1.push_back(1503.66);
    par1.push_back(387.094);
    par1.push_back(2.68137e-09);
    par1.push_back(par[3]);
    par1.push_back(0.0);
    par1.push_back(1.);
    par1.push_back(0.0);
    par1.push_back(0.0);
    
    std::complex<double> ATot = KnFF_total(Ecms[0] ,par1);
    double MPhi = par[0];
    auto ff =  ATot * sqrt(3 / (8. * alphafss * alphafss * pi * pow(MPhi*MPhi/4.-mKn*mKn, 1.5)));

    // if S < threshold then Sigma = 0
    double M_K0 = 497.614;
    return (s < 4 * M_K0 * M_K0)? 0. : ff;
}

// ------------ Auxiliary functions ------------

// Cross section of e+e- -> KsKl for point-like kaons.
Double_t PointLikeKnCrossSection(Double_t s)
{ 
    double beta = sqrt(1 - 4*mKn*mKn/s);
    return pi * alphafss * alphafss * pow(beta,3.)/(3*s)*C; 
}

Double_t RhoWidth(Double_t s)
{
    Double_t MRho = 775.49;
    Double_t WRho = 149.10;

    Double_t BrKC = 0.491;
    Double_t BrK0 = 0.339;

    Double_t BrPi0Gamma = 0.0006;
    Double_t BrEtaGamma = 0.0003;

    Double_t Mpi = 139.57018;
    Double_t MPhi = 1019.461;

    Double_t PhiWidth(Double_t);
    Double_t PKPKM(Double_t);
    Double_t PKLKS(Double_t);
    Double_t QPGamma(Int_t, Double_t);
    
    Double_t Fval = WRho/s*MRho*MRho*TMath::Power((s/4.0-Mpi*Mpi)/(MRho*MRho/4.0-Mpi*Mpi),1.5)
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
    
    Double_t MPi  = 139.56995;
    Double_t MPi0 = 134.9766;
    Double_t BrKC = 0.491;
    Double_t BrK0 = 0.339;
    Double_t BrEtaGamma = 0.00065;
    Double_t BrOmg3Pi   = 0.891;
    Double_t BrOmg2Pi   = 0.017;
    Double_t BrOmgPiG   = 0.087;

    Double_t PhiWidth(Double_t);
    Double_t PKPKM(Double_t);
    Double_t PKLKS(Double_t);
    Double_t QPGamma(Int_t, Double_t);
    Double_t FAS_ASPO(Double_t);

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

Double_t PhiWidth(Double_t s)
{
    Double_t MPhi = 1019.461;
    Double_t WPhi = 4.249;
    Double_t MPi  = 139.569995;
    Double_t MPi2;

    Double_t Br[4] = {0.491, 0.339, 0.1532, 0.01309};
    Double_t Br2Pi = 0.000073;
    Double_t BrPi0Gamma = 0.00124;

    Double_t MPhi2 = 1019.461*1019.461;

    Double_t PKPKM(Double_t);
    Double_t PKLKS(Double_t);
    Double_t QPGamma(Int_t, Double_t);

    Double_t FAS_ASPO(Double_t);

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

    if ( Mode == 0 ) PMass = 134.9766;
    if ( Mode == 1 ) PMass = 547.862;

    if ( s < PMass*PMass ) return 0.0;
    
    Double_t Fval =  TMath::Power(0.5*(1.0-PMass*PMass/s),3)*TMath::Power(s,1.5);

    return Fval;
}

Double_t FAS_ASPO(Double_t TwoE)
{
    // Double_t Polinom(Int_t, Double_t, Double_t*);

    const Int_t K_max = 4;
    Int_t I_mark = 100;
    Int_t Mode = 1;

    bool Calc = true;

    Double_t A[4] = {-6.12394622, 25.0341405, -34.1311022, 15.5413717};
    Double_t B[4] = {5.29354148, -7.90990714, -2.26007613, 5.21453902};
    double CC[4] = {-0.5643611, 2.69953793, -4.32966739, 2.33116866};
    Double_t D[4] = {-0.0548334238,  0.31600391, -0.609523718, 0.393667808};

    Double_t A1[4] = {-4.91624401 , 19.8655606, -26.9136128 , 12.2412286};
    Double_t D1[4] = {-0.00794774189,0.0522269164,-0.114526409 , 0.0838126536};
    
    Double_t LM[6] = {1.1,0.875,0.75,0.62,0.52,0.46};

    Double_t POL, TEMP1, TEMP2; 

    if ( Calc ) {
        TEMP2 = D[1]*LM[4]+D[2]*LM[4]*LM[4]+D[3]*LM[4]*LM[4]*LM[4];
        TEMP1 = D1[0]+D1[1]*LM[4]+D1[2]*LM[4]*LM[4]+D1[3]*LM[4]*LM[4]*LM[4];
        D[0]  = TEMP1 - TEMP2;
    
        TEMP2 = CC[1]*LM[3]+CC[2]*LM[3]*LM[3]+CC[3]*LM[3]*LM[3]*LM[3];
        TEMP1 = D[0]+D[1]*LM[3]+D[2]*LM[3]*LM[3]+D[3]*LM[3]*LM[3]*LM[3];
        CC[0]  = TEMP1 - TEMP2;
    
        TEMP2 = A1[1]*LM[2]+A1[2]*LM[2]*LM[2]+A1[3]*LM[2]*LM[2]*LM[2];
        TEMP1 = CC[0]+CC[1]*LM[2]+CC[2]*LM[2]*LM[2]+CC[3]*LM[2]*LM[2]*LM[2];
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
        if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*CC[i];
        if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*CC[i];
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

int KnXSec_func()
{
    // data from E.A.Kozyrev's article (https://arxiv.org/abs/1604.02981).
    std::vector<double> energies = {1004.066, 1010.466, 1012.955, 1015.068, 1016.105, 1017.155, 1017.156, 1018.046, 1019.118,  1019.214, 1019.421, 1019.902, 
                1021.222, 1021.309, 1022.078,  1022.744, 1023.264,  1025.320, 1027.956, 1029.090, 1033.907, 1040.028,  1049.864,  1050.862,  1059.947};
    std::vector<double> cross_sections = {6.87, 42.16, 96.74, 219.53, 366.33, 628.15, 624.76, 996.62, 1413.65,  1433.05, 1434.84, 1341.91,  833.20, 
                        807.54, 582.93,  443.71, 377.77,  199.26, 115.93, 96.96, 50.12, 31.27,  16.93,  17.47,  12.09};
    std::vector<double> xsec_err = {0.42, 0.47, 1., 5.02, 3.33, 2.95, 9.89, 4.28,  6.02, 15.03, 18.40, 4.74, 4.89, 10.36, 4.03, 4.38, 
                                    5.31, 4.97, 1.7, 3., 1.26, 1.01, 0.5, 0.94, 0.71};
    std::vector<double> zeroes(100, 0.0);
    
    TCanvas canv("canv", "canv");
    TGraphErrors xsec(energies.size(), energies.data(), cross_sections.data(), zeroes.data(), xsec_err.data());

    auto born_xsec = new TF1("born", KnCrossSection, 1000., 1060., 4);
    
    born_xsec->SetParameters(1019.457, 4.240, 0.428e-03, 0.8);
    born_xsec->SetNpx(1e4);
    // born_xsec->DrawClone();

    std::vector<double> delta_cross_sections = {};
    for(const auto& en : energies)
    {
        // delta_cross_sections.push_back(born_xsec->EvalPar([en], [1019.457, 4.270, 0.428e-03, 0.8) - ))
    }

    // xsec.Fit(born_xsec, "ME");
    // xsec.DrawClone("AP");
    // born_xsec->DrawClone("same");

    // born_xsec->DrawDerivative();
    born_xsec->SetParameters(1019.457, 4.270, 0.428e-03, 0.8);
    born_xsec->SetLineColor(kRed);
    // born_xsec->DrawDerivative("same");

    double_t pars1[4] = {1019.457, 4.240, 0.428e-03, 0.8};
    double_t pars2[4] = {1019.457, 4.270, 0.428e-03, 0.8};
    double_t pars_mine[4] = {1019.446, 4.225, 0.428e-03, 0.8};

    auto born_xsec_deriv_diff = new TF1("born", [&](double *x, double *pars) {
        return born_xsec->EvalPar(x, pars1) - born_xsec->EvalPar(x, pars2); 
        // return born_xsec->Derivative(x[0], pars1) - born_xsec->Derivative(x[0], pars2);
    }, 1000., 1060., 1);
    born_xsec_deriv_diff->SetNpx(1e4);
    born_xsec_deriv_diff->Draw();
    canv.DrawClone();

    return 0;
}

class KnXSec
{
private:
    std::unique_ptr<TF1> born_xsec;

public:
    KnXSec();
    ~KnXSec() = default;

    double eval(double energy)
    {
        return born_xsec->Eval(2 * energy);
    }
};

KnXSec::KnXSec()
{
    born_xsec = std::unique_ptr<TF1>(new TF1("born", KnCrossSection, 990., 1060., 4));
    born_xsec->SetParameters(1019.457, 4.240, 0.428e-03, 0.8);
    born_xsec->SetNpx(1e4);
}

