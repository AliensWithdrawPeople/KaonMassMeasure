#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"

int MassVsE()
{
    double massK = 497.614;
    std::vector<Float_t> zeroes(100, 0.0);

    std::vector<Float_t> vMPDG = {497.742, 497.661, 497.625, 497.634, 497.583, 497.607};
    std::vector<Float_t> vMPDGerr = {0.085, 0.033, 0.031, 0.024, 0.021, 0.017};
    std::vector<Float_t> vExpNum = {505, 506, 507, 508, 509, 510};
    TGraphErrors grMPDG(vMPDG.size(), vExpNum.data(), vMPDG.data(), zeroes.data(), vMPDGerr.data());
    // grMPDG.DrawClone("same AP");
    // massKline->DrawClone("Same");

    std::vector<Float_t> vPavg = {157, 182, 207, 231, 235, 243, 254, 260};

    std::vector<Float_t> vM_Exp = {497.541, 497.564, 497.536, 497.563, 497.549, 497.578, 497.589, 497.571, 497.601};
    std::vector<Float_t> vM_Exp2 = {497.569, 497.579, 497.541, 497.574, 497.566, 497.601, 497.591, 497.571, 497.529};
    std::vector<Float_t> vMerrExp = {0.010, 0.012, 0.013, 0.013, 0.011, 0.012, 0.013, 0.014, 0.024};
    std::vector<Float_t> vE = {505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 514};
    TGraphErrors grRCNC_exp(vM_Exp.size(), vE.data(), vM_Exp.data(), zeroes.data(), vMerrExp.data());
    TGraphErrors grRCNC_exp2(vM_Exp.size(), vE.data(), vM_Exp2.data(), zeroes.data(), vMerrExp.data());

    grRCNC_exp.SetTitle("Black -- with Energy controll, Blue -- without");
    grRCNC_exp.GetXaxis()->SetTitle("E_{beam}, MeV");
    grRCNC_exp.GetYaxis()->SetTitle("M^{(FullRec)}_{RC NC RecCor}, #frac{MeV}{c^{2}}");
    grRCNC_exp.SetName("EnControl");

    grRCNC_exp.DrawClone("AP");
    grRCNC_exp2.SetMarkerColor(kBlue);
    // grRCNC_exp2.DrawClone("P Same");

    return 0;
}