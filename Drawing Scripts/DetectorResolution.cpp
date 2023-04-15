#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TF1.h" 

using vec = std::vector<double>;


int DetectorResolution()
{
/*
***************************************************
* Detector resolution vs lnY (E_point = 510.5 MeV):
***************************************************
*/
    vec zeroes(100, 0.0);
    vec lnY = {-0.375, -0.325, -0.275, -0.225, -0.175, -0.125, -0.075, -0.025, 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375};
    vec xRange(100, 0.025);

    vec exp510_5_res = {0.0293367, 0.0257106, 0.0230893, 0.0202052, 0.0192319, 0.0173616, 0.0166222, 0.0163964, 0.0162607, 0.0164867, 0.0170824, 0.0185048, 0.0201503, 0.0232465, 0.0253798, 0.0296679}; 
    vec exp510_5_resErr = {0.000403798, 0.000368311, 0.000351037, 0.000294031, 0.000277827, 0.000276993, 0.000276293, 0.00024957, 0.000257862, 0.000267686, 0.00026952, 0.000287248, 0.000301728, 0.000341351, 0.000354024, 0.0004338}; 
    vec exp510_5_resRMS = {0.0353734, 0.0324411, 0.0285595, 0.0267774, 0.0249735, 0.0233319, 0.0230339, 0.0217395, 0.0230499, 0.0224724, 0.0233741, 0.0253803, 0.0263835, 0.0288288, 0.0323426, 0.0357933};
    vec exp510_5_resRMSerr = {0.000317791, 0.000292179, 0.000262911, 0.000252549, 0.000235936, 0.000222198, 0.00022381, 0.000211352, 0.000224624, 0.000216925, 0.000225524, 0.000240574, 0.000248217, 0.000264029, 0.000293081, 0.000321641};
    
    vec MC510_5_res = {0.0287751, 0.0240572, 0.0222305, 0.018954, 0.0168579, 0.0178757, 0.014925, 0.0160412, 0.0148722, 0.0159862, 0.0161041, 0.0179021, 0.0185923, 0.020374, 0.0242676, 0.0285658};
    vec MC510_5_resErr = {0.00112921, 0.000952216, 0.00095722, 0.000758307, 0.000766868, 0.000768446, 0.000682797, 0.000783418, 0.000688001, 0.000773766, 0.000915849, 0.000803759, 0.000834932, 0.000820302, 0.00106843, 0.00118497};
    vec MC510_5_resRMS = {0.0338171, 0.0338478, 0.0269535, 0.0245841, 0.0226219, 0.0214614, 0.0203239, 0.0226417, 0.0214698, 0.0222647, 0.0236017, 0.0234514, 0.0261891, 0.0289492, 0.0316015, 0.0338511};
    vec MC510_5_resRMSerr = {0.000882619, 0.000861404, 0.000696401, 0.000656101, 0.000638314, 0.000598931, 0.000559822, 0.000651983, 0.000602456, 0.000630751, 0.000671328, 0.000646458, 0.000696953, 0.000768231, 0.00080528, 0.000828352};
    
    vec MC510_5Smeared_res = {0.029446, 0.0253757, 0.0226888, 0.0201915, 0.0184904, 0.0170543, 0.0162259, 0.0156777, 0.0160234, 0.016628, 0.0169716, 0.0183647, 0.0209968, 0.0225542, 0.0253334, 0.0300803};
    vec MC510_5Smeared_resErr = {0.000377701, 0.000324142, 0.000306126, 0.000270792, 0.00025596, 0.000248088, 0.000234514, 0.000219965, 0.000238442, 0.000237539, 0.000244739, 0.000257661, 0.000290025, 0.000308156, 0.000321102, 0.000375205};
    vec MC510_5Smeared_resRMS = {0.035772, 0.0324929, 0.0298521, 0.0265526, 0.0243349, 0.0227267, 0.0211992, 0.0212944, 0.0219155, 0.0216133, 0.0230264, 0.0242471, 0.0264392, 0.0297325, 0.0326509, 0.0351503};
    vec MC510_5Smeared_resRMSerr = {0.000289067, 0.000267018, 0.000251079, 0.000227452, 0.000213414, 0.000203518, 0.000190852, 0.000191523, 0.000198332, 0.000193269, 0.000204568, 0.000213104, 0.000228298, 0.000250517, 0.000268715, 0.000282116};
    
    TGraphErrors grExp(lnY.size(), lnY.data(), exp510_5_res.data(), xRange.data(), exp510_5_resErr.data());
    TGraphErrors grExp_RMS(lnY.size(), lnY.data(), exp510_5_resRMS.data(), xRange.data(), exp510_5_resRMSerr.data());

    TGraphErrors grMC(lnY.size(), lnY.data(), MC510_5_res.data(), xRange.data(), MC510_5_resErr.data());
    TGraphErrors grMC_RMS(lnY.size(), lnY.data(), MC510_5_resRMS.data(), xRange.data(), MC510_5_resRMSerr.data());

    TGraphErrors grMCsmeared(lnY.size(), lnY.data(), MC510_5Smeared_res.data(), xRange.data(), MC510_5Smeared_resErr.data());
    TGraphErrors grMCsmeared_RMS(lnY.size(), lnY.data(), MC510_5Smeared_resRMS.data(), zeroes.data(), MC510_5Smeared_resRMSerr.data());

    grExp.SetTitle("E_{point} = 510.5 MeV, Black - data, Red - MC, Blue - MC with energy smearing");
    grExp.GetXaxis()->SetTitle("lnY");
    grExp.GetYaxis()->SetTitle("#sigma_{#psi}, rad");

    grExp.SetMarkerColor(kBlack);
    grExp.SetMarkerSize(1.5);
    grExp.SetMarkerStyle(22);
    
    grExp_RMS.SetMarkerColor(kBlack);
    grExp_RMS.SetMarkerSize(1.5);
    grExp_RMS.SetMarkerStyle(22);

    grMC.SetMarkerColor(kRed);
    grMC.SetMarkerSize(1.5);
    grMC.SetMarkerStyle(22);

    grMC_RMS.SetMarkerColor(kBlack);
    grMC_RMS.SetMarkerSize(1.5);
    grMC_RMS.SetMarkerStyle(23);

    grMCsmeared.SetMarkerColor(kBlue);
    grMCsmeared.SetMarkerSize(1.5);

    grMCsmeared_RMS.SetMarkerColor(kBlack);
    grMCsmeared_RMS.SetMarkerSize(1.5);
    grMCsmeared_RMS.SetMarkerStyle(23);

    grExp.DrawClone("AP");
    grMC.DrawClone("P same");
    grMCsmeared.DrawClone("P same");

/*
*****************************************
* deltaM_NC vs lnY (E_point = 510.5 MeV):
*****************************************
*/
    vec deltaM_exp = {0.036, 0.027, 0.022, 0.017, 0.015, 0.0125, 0.0115, 0.011, 0.011, 0.0115, 0.0125, 0.015, 0.017, 0.022, 0.027, 0.036};
    vec deltaM_MC = {0.034, 0.024, 0.020, 0.0145, 0.0135, 0.0115, 0.0105, 0.0095, 0.0095, 0.0105, 0.0135, 0.0145, 0.0155, 0.021, 0.025, 0.034};
    vec deltaM_MCsmeared = {0.036, 0.026, 0.022, 0.017, 0.015, 0.0125, 0.0115, 0.011, 0.011, 0.0115, 0.0125, 0.015, 0.017, 0.022, 0.026, 0.036};

    TGraphErrors grDeltaM_exp(lnY.size(), lnY.data(), deltaM_exp.data(), xRange.data(), zeroes.data());
    TGraphErrors grDeltaM_MC(lnY.size(), lnY.data(), deltaM_MC.data(), xRange.data(), zeroes.data());
    TGraphErrors grDeltaM_MCsmeared(lnY.size(), lnY.data(), deltaM_MCsmeared.data(), xRange.data(), zeroes.data());
    
    grDeltaM_exp.SetName("grDeltaM_exp");
    grDeltaM_MC.SetName("grDeltaM_MC");
    grDeltaM_MCsmeared.SetName("grDeltaM_MCsmeared");

    grDeltaM_exp.SetTitle("E_{point} = 510.5 MeV, Black - data, Red - MC, Blue - MC with energy smearing");
    grDeltaM_exp.GetXaxis()->SetTitle("lnY");
    grDeltaM_exp.GetYaxis()->SetTitle("#DeltaM_{NC}, MeV");

    grDeltaM_exp.SetMarkerColor(kBlack);
    grDeltaM_exp.SetMarkerSize(1.5);
    grDeltaM_exp.SetMarkerStyle(22);

    grDeltaM_MC.SetMarkerColor(kRed);
    grDeltaM_MC.SetMarkerSize(1.5);
    grDeltaM_MC.SetMarkerStyle(23);

    grDeltaM_MCsmeared.SetMarkerColor(kBlue);
    grDeltaM_MCsmeared.SetMarkerSize(1.5);

    grDeltaM_exp.DrawClone("AP");
    grDeltaM_MC.DrawClone("P same");
    grDeltaM_MCsmeared.DrawClone("P same");

    return 0;
}