#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"

int drawGraphs()
{
    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};
    std::vector<Float_t> zeroes(18, 0.0);

    std::vector<Float_t> vM_FullRec = {497.703, 497.678, 497.661, 497.701, 497.937, 499.116};
    std::vector<Float_t> vMerr_FullRec = {.007, 0.008, 0.008, 0.008, 0.010, 0.019};
    TGraphErrors grRCNC_FullRec(vM_FullRec.size(), vE.data(), vM_FullRec.data(), zeroes.data(), vMerr_FullRec.data());
    grRCNC_FullRec.DrawClone("AP");

    std::vector<Float_t> vM_CrAngle = {497.685, 497.695, 497.67, 497.728, 497.914, 499.079};
    std::vector<Float_t> vMerr_CrAngle = {.012, 0.016, 0.016, 0.018, 0.021, 0.048};
    TGraphErrors grRCNC_CrAngle(vM_CrAngle.size(), vE.data(), vM_CrAngle.data(), zeroes.data(), vMerr_CrAngle.data());
    // grRCNC_CrAngle.DrawClone("P same");

    std::vector<Float_t> deltaPhiPos = {-0.000312404, -0.000589547, -0.000951601, -0.000464232, -0.000176916, -0.00073326, -0.000296456, -8.5357e-05, -0.0011892, -0.000752748, -0.000346934, -0.000425162, 0.000173945, 0.000198806, -0.000905391, -0.000399599, -0.000453043, -0.000500869};
    std::vector<Float_t> deltaPhiNeg = {0.00115229, -0.00068634, -0.000889196, -3.89684e-05, -0.000217931, -0.000180165, 0.000248681, 0.00039728, -0.000247258, 0.000950679, 0.000707743, -0.000310605, 0.000713548, -2.31645e-05, -0.000127349, -3.94741e-05, -0.000187029, -0.000562813};
    std::vector<Float_t> deltaPtotPos = {0.0549019, 0.105076, 0.321787, 0.108172, -0.00915358, 0.232584, 0.108496, 0.120291, 0.314675, 0.206091, 0.0668328, 0.218722, 0.143755, 0.161072, 0.277399, 0.104557, 0.0804615, 0.117139};
    std::vector<Float_t> deltaPtotNeg = {0.2979, -0.0906616, -0.0995954, 0.0931528, 0.0907167, -0.0342835, -0.0125813, 0.0551172, -0.0988879, 0.241792, 0.12736, -0.103945, 0.167163, 0.0427573, 0.0573838, -0.0314714, -0.0755831, -0.158469};
    std::vector<Float_t> deltaPhiPosErr =  {0.000149231, 0.000152889, 0.000223089, 0.000179193, 0.000206281, 0.000256816, 0.000262754, 0.000252833, 0.000210241, 0.000217408, 0.000231568, 0.000281043, 0.000281636, 0.000251983, 0.000192513, 0.000205022, 0.000159534, 0.000169987};
    std::vector<Float_t> deltaPhiNegErr = {0.000235513, 0.000271333, 0.000287681, 0.000301175, 0.000265405, 0.000178755, 0.000165038, 0.000162146, 0.000145813, 0.000160739, 0.000176715, 0.000160374, 0.000173009, 0.000233536, 0.000215724, 0.000259462, 0.000240372, 0.000213341};
    std::vector<Float_t> deltaPtotPosErr = {0.0440824, 0.0483945, 0.0751344, 0.0535507, 0.0488162, 0.0467227, 0.0437407, 0.0474076, 0.0436462, 0.0409416, 0.0464952, 0.0433387, 0.0449354, 0.0447475, 0.0470045, 0.047424, 0.0434876, 0.0451984};
    std::vector<Float_t> deltaPtotNegErr = {0.0560044, 0.0757615, 0.0492184, 0.0474808, 0.049744, 0.04239, 0.0431219, 0.0467077, 0.0413148, 0.0408802, 0.0459421, 0.041855, 0.043048, 0.0500148, 0.0442303, 0.0477249, 0.0517014, 0.0448925};


    std::vector<Float_t> vCellNums = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    auto grDeltaPtotPos1 = new TGraphErrors(deltaPtotPos.size(), vCellNums.data(), deltaPtotPos.data(), zeroes.data(), deltaPtotPosErr.data());
    auto grDeltaPtotNeg1 = new TGraphErrors(deltaPtotNeg.size(), vCellNums.data(), deltaPtotNeg.data(), zeroes.data(), deltaPtotNegErr.data());

    auto grDeltaPhiPos1 = new TGraphErrors(deltaPhiPos.size(), vCellNums.data(), deltaPhiPos.data(), zeroes.data(), deltaPhiPosErr.data());
    auto grDeltaPhiNeg1 = new TGraphErrors(deltaPhiNeg.size(), vCellNums.data(), deltaPhiNeg.data(), zeroes.data(), deltaPhiNegErr.data());
    
    grDeltaPtotPos1->SetMarkerColor(kRed);
    grDeltaPtotNeg1->SetMarkerColor(kBlue);

    grDeltaPhiPos1->SetMarkerColor(kRed);
    grDeltaPhiNeg1->SetMarkerColor(kBlue);

    grDeltaPtotPos1->SetTitle("#pi^{+} - red, #pi^{-} - blue;Cell num;#DeltaP_{tr}, #frac{MeV}{c}");
    // grDeltaPtotPos1->SetTitle("#DeltaP_{#pi^{+}} - #DeltaP_{#pi^{-}};Cell num;#DeltaP_{tr}, #frac{MeV}{c}");
    grDeltaPhiPos1->SetTitle("#pi^{+} - red, #pi^{-} - blue;Cell num;#Delta#phi, rad");
    // grDeltaPtotPos1->SetTitle("Red - #pi^{+}, Blue - #pi^{-}");
    grDeltaPtotPos1->DrawClone("AP");
    grDeltaPtotNeg1->DrawClone("P same");
    
    return 0;
}