#ifndef _Phi_TO_K0_XSecBorn_Func_HPP_
#define _Phi_TO_K0_XSecBorn_Func_HPP_
// #include "TMath.h"
// #include "TComplex.h"
// #include "TF1.h"
// #include "TSpline.h"

namespace Ivanov_PhiToK0{
    Double_t RhoWidth(Double_t s);

    Double_t OmgWidth(Double_t s);

    Double_t PhiWidth(Double_t s, Double_t MPhi = 1019.461, Double_t WPhi = 4.249);

    Double_t PKPKM(Double_t s);

    Double_t PKLKS(Double_t s);

    Double_t QPGamma(Int_t Mode, Double_t s);

    Double_t PhitoK0(double x, const std::vector<double>& par);

    Double_t eff1(double s);
}

#endif