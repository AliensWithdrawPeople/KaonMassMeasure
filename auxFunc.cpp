#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"

int auxFunc()
{
    double pRatio = 0.99999;
    double psi = 2.56547;

    double energy = 514;
    double mass =  499.109;
    double massK = 497.614;

    std::vector<Float_t> zeroes(100, 0.0);


    auto massFunc = new TF1("MassLnY", "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
    massFunc->SetParameter(1, (1 - pRatio*pRatio) / (1 + pRatio*pRatio));

    //auto tmpFunc = new TF1("MassLnY", 
    //"sqrt(x * x * (1 - (1 + sqrt(1 - [1] *[1]) * cos([0]))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / x / x)))/ [0] / [0] ))", 480, 520);
    auto tmpFunc = new TF1("MassLnY", "sqrt(x * x * (1 - (1 + sqrt(1 - [1] *[1]) * cos([0]))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / x / x)))/ [1] / [1] ))");
    tmpFunc->SetParameters(psi, (1 - pRatio*pRatio) / (1 + pRatio*pRatio));
    massFunc->SetParameter(0, energy);
    double Emin = tmpFunc->GetX(490, 450, 550, 1e-5, 10000);
    double Emax = tmpFunc->GetX(505, 450, 550, 1e-5, 10000);
    std::cout<< Emin << " : " << Emax << std::endl;


    auto revMassFunc = new TF1("revMassFunc", "2*TMath::ACos(TMath::Sqrt((1 - x * x / [0] / [0])/(1-4*139.57 * 139.57 / [0] / [0])))");
    revMassFunc->SetParameter(0, energy);
    std::cout << "At E = " << energy << " and M = " << mass << " Psi = " << revMassFunc->Eval(mass) << std::endl;


    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};

    // Full reconstruction
    // Avg energy radcor
    std::vector<Float_t> vMRad2 = { 497.597, 497.593, 497.597, 497.583, 497.617, 497.642};
    std::vector<Float_t> vMRad2NC = {497.609, 497.617, 497.608, 497.611, 497.617, 497.612};
    std::vector<Float_t> vMerrRad2 = {0.007, 0.009, 0.008, 0.010, 0.010, 0.019};
    TGraphErrors grRC(vMRad2.size(), vE.data(), vMRad2.data(), zeroes.data(), vMerrRad2.data());
    TGraphErrors grRCNC(vMRad2NC.size(), vE.data(), vMRad2NC.data(), zeroes.data(), vMerrRad2.data());
    // grRC2.SetMarkerColor(kRed);
    
    // RMS
    // std::vector<Float_t> vM2 = { 497.597, 497.600, 497.597, 497.587, 497.606, 497.581};
    std::vector<Float_t> vM2 = {497.620, 497.618, 497.615, 497.610, 497.633, 497.704};
    std::vector<Float_t> vMerr2 = {0.007, 0.009, 0.008, 0.010, 0.010, 0.019};
    TGraphErrors grRCNCnew(vM2.size(), vE.data(), vM2.data(), zeroes.data(), vMerr2.data());
    grRCNCnew.SetMarkerColor(kBlack);

    std::vector<Float_t> vDeltaMRad = {94, 61, 53, 90, 320, 1504};
    std::vector<Float_t> vDeltaMRadErr = { 0.008, 0.007, 0.004, 0.010, 0.013, 0.017};
    TGraphErrors grRCdeltaM(vDeltaMRad.size(), vE.data(), vDeltaMRad.data(), zeroes.data(), vDeltaMRadErr.data());
    // grRCdeltaM.SetMarkerColor(kGreen);

    // deltaM_NC
    std::vector<Float_t> vDeltaM_NC_RMS = {0.023, 0.018, 0.018, 0.023, 0.027, 0.123}; //RMS sigma
    std::vector<Float_t> vDeltaM_NC_Fit = { 0.012, 0.017, 0.011, 0.024, 0.011, 0.031}; //Fit sigma
    std::vector<Float_t> vDeltaM_NC_FitErr = {0.0006, 0.0004, 0.0006, 0.0010, 0.0007, 0.0009};
    std::vector<Float_t> vDeltaM_NC_RMSerr = {0.0006, 0.0004, 0.0006, 0.0010, 0.0021, 0.0050};
    TGraphErrors grDeltaM_NC_RMS(vDeltaM_NC_RMS.size(), vE.data(), vDeltaM_NC_RMS.data(), zeroes.data(), vDeltaM_NC_RMSerr.data());
    TGraphErrors grDeltaM_NC_Fit(vDeltaM_NC_Fit.size(), vE.data(), vDeltaM_NC_Fit.data(), zeroes.data(), vDeltaM_NC_FitErr.data());

    // Critical angle
    std::vector<Float_t> vMCrAngle = {497.685, 497.695, 497.67, 497.728, 497.914, 499.079};
    std::vector<Float_t> vMCrAngleRCNC_Fit = {497.596, 497.615, 497.612, 497.621, 497.598, 497.582}; //Fit sigma
    std::vector<Float_t> vMCrAngleRCNC_RMS = {497.608, 497.643, 497.608, 497.634, 497.646, 497.8}; //RMS sigma

    std::vector<Float_t> vMCrAngleErr = {0.010, 0.013, 0.013, 0.014, 0.015, 0.034};
    std::vector<Float_t> vMCrAngleRCNC_FitErr = {0.014, 0.016, 0.016, 0.018, 0.021, 0.057};
    std::vector<Float_t> vMCrAngleRCNC_RMSErr = {0.010, 0.013, 0.013, 0.013, 0.015, 0.034};

    TGraphErrors grMCrAngle(vMCrAngle.size(), vE.data(), vMCrAngle.data(), zeroes.data(), vMCrAngleErr.data());
    TGraphErrors grMCrAngleRCNC_Fit(vMCrAngleRCNC_Fit.size(), vE.data(), vMCrAngleRCNC_Fit.data(), zeroes.data(), vMCrAngleRCNC_FitErr.data());
    TGraphErrors grMCrAngleRCNC_RMS(vMCrAngleRCNC_RMS.size(), vE.data(), vMCrAngleRCNC_RMS.data(), zeroes.data(), vMCrAngleRCNC_RMSErr.data());

    auto massKline = new TLine(505, massK, 514, massK);
    massKline->SetLineColor(kBlue);
    massKline->SetLineWidth(2);

    grMCrAngle.SetMarkerColor(kRed);
    grMCrAngleRCNC_RMS.SetMarkerColor(kRed);
    grMCrAngleRCNC_Fit.SetMarkerColor(kBlack);
    grRCNCnew.SetMarkerColor(kRed);

    // grRCNC.DrawClone("AP");
    // grRCNCnew.DrawClone("AP same");
    grDeltaM_NC_Fit.DrawClone("AP");
    grDeltaM_NC_RMS.DrawClone("P same");
    // massKline->DrawClone("Same");

    std::vector<Float_t> vSigmaFit = {0.0139273, 0.0142717, 0.0153109, 0.0176764, 0.0188618, 0.0220763, 0.025119, 0.0298211};
    std::vector<Float_t> vSigmaRMS = {0.0186562, 0.019872, 0.0214141, 0.0219208, 0.0253945, 0.0270087, 0.0295075, 0.0336299};
    std::vector<Float_t> vSigmaErr = {2.41164e-04, 2.59689e-04, 2.51916e-04, 3.13132e-04, 3.32211e-04, 4.86097e-04, 4.64792e-04, 6.67834e-04};
    std::vector<Float_t> vWidth = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375};
    TGraphErrors grSigmaFit(vSigmaFit.size(), vWidth.data(), vSigmaFit.data(), zeroes.data(), vSigmaErr.data());
    TGraphErrors grSigmaRMS(vSigmaRMS.size(), vWidth.data(), vSigmaRMS.data(), zeroes.data(), vSigmaErr.data());

    // grSigmaRMS.SetMarkerColor(kRed);
    // grSigmaRMS.DrawClone("AP"); 
    // grSigmaFit.DrawClone("P same"); 

    std::vector<Float_t> vAngle = { 2.73092, 2.65629, 2.63373, 2.62245, 2.61362, 2.59637, 2.55372};
    std::vector<Float_t> vEtmp = {505, 508, 509, 509.527, 510, 511, 514};
    TGraph grAngleVsE(vEtmp.size(), vEtmp.data(), vAngle.data());
    std::cout << "Psi angle vs Energy correlation factor = " << grAngleVsE.GetCorrelationFactor() << std::endl;

    std::vector<Float_t> vMPDG = {497.742, 497.661, 497.625, 497.634, 497.583, 497.607};
    std::vector<Float_t> vMPDGerr = {0.085, 0.033, 0.031, 0.024, 0.021, 0.017};
    std::vector<Float_t> vExpNum = {505, 506, 507, 508, 509, 510};
    TGraphErrors grMPDG(vMPDG.size(), vExpNum.data(), vMPDG.data(), zeroes.data(), vMPDGerr.data());

    // grMPDG.DrawClone("same AP");
    // massKline->DrawClone("Same");

    // gr2.DrawClone("AP");

    // Positive tracks
    std::vector<Float_t> vDeltaPtot = {0.076, 0.096, 0.126, 0.165, 0.178, 0.209, 0.248, 0.349};
    std::vector<Float_t> vDeltaPtotErr = {0.012, 0.013, 0.017, 0.015, 0.016, 0.025, 0.017, 0.018};

    std::vector<Float_t> vDeltaPhi = {-0.00076, -0.00089, -0.00072, -0.00055, -0.00064, -0.0006, -0.00065, -0.00083};
    std::vector<Float_t> vDeltaPhiErr = {0.00019, 0.00014, 0.00009, 0.00004, 0.00002, 0.000065, 0.00005, 0.00005};

    std::vector<Float_t> vPavg = {157, 182, 207, 231, 235, 243, 254, 260};

    // for(int i = 0; i < vDeltaPtot.size(); i++)
    // {
    //     vDeltaPtot[i] = vDeltaPtot[i] / vPavg[i];
    //     vDeltaPtotErr[i] = vDeltaPtotErr[i] / vPavg[i];
    // }
    
    TGraphErrors grDeltaPtotPos(vDeltaPtot.size(), vPavg.data(), vDeltaPtot.data(), zeroes.data(), vDeltaPtotErr.data());
    TGraphErrors grDeltaPhiPos(vDeltaPhi.size(), vPavg.data(), vDeltaPhi.data(), zeroes.data(), vDeltaPhiErr.data());

    // Negative tracks
    vDeltaPtot = {0.027, 0.059,  0.059, 0.038, 0.047, 0.063, 0.051, 0.15};
    vDeltaPtotErr = {0.014, 0.016, 0.015, 0.015, 0.002, 0.03, 0.018, 0.02};

    vDeltaPhi = {0.00019, 0.00034, 0.00022, 0.00008, -0.00019, -0.00023, -0.00019, -0.00041};
    vDeltaPhiErr = {0.00019, 0.00013, 0.00008, 0.00005, 0.00002, 0.00002, 0.00005, 0.00005};

    TGraphErrors grDeltaPtotNeg(vDeltaPtot.size(), vPavg.data(), vDeltaPtot.data(), zeroes.data(), vDeltaPtotErr.data());
    TGraphErrors grDeltaPhiNeg(vDeltaPhi.size(), vPavg.data(), vDeltaPhi.data(), zeroes.data(), vDeltaPhiErr.data());

    // grDeltaPtotPos.DrawClone();
    // grDeltaPtotNeg.DrawClone("same");
    grDeltaPhiPos.SetTitle("Black - #pi^{+}, Red - #pi^{-}");
    grDeltaPhiNeg.DrawClone();
    // grDeltaPhiNeg.DrawClone("same");


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

    std::vector<Float_t> vDeltaPtotPos;
    std::vector<Float_t> vDeltaPtotNeg;
    std::vector<Float_t> vDeltaPhiPos;
    std::vector<Float_t> vDeltaPhiNeg;
    
    int a = 3;
    for(int i = 0; i < dPtotPos.size(); i++)
    {
        vDeltaPtotPos.push_back(dPtotPos[i][a] - dPtotNeg[i][a]);
        vDeltaPtotNeg.push_back(dPtotNeg[i][a]);
        vDeltaPhiPos.push_back(dPhiPos[i][a] - dPhiNeg[i][a]);
        vDeltaPhiNeg.push_back(dPhiNeg[i][a]);
    }
    std::vector<Float_t> vCellNums = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    std::vector<Float_t> vErrs1 = {0.0475, 0.05, 0.061, 0.0864, 0.02, 0.04, 0.05, 0.06 };
    std::vector<Float_t> vErrs2 = {0.0002, 0.0003, 0.0005, 0.0009, 0.0001, 0.0002, 0.0003, 0.0004 };
    auto grDeltaPtotPos1 = new TGraphErrors(vDeltaPtotPos.size(), vPavg.data(), vDeltaPtotPos.data(), zeroes.data(), vErrs1.data());
    auto grDeltaPtotNeg1 = new TGraphErrors(vDeltaPtotNeg.size(), vPavg.data(), vDeltaPtotNeg.data(), zeroes.data(), zeroes.data());

    auto grDeltaPhiPos1 = new TGraphErrors(vDeltaPhiPos.size(), vPavg.data(), vDeltaPhiPos.data(), zeroes.data(), vErrs2.data());
    auto grDeltaPhiNeg1 = new TGraphErrors(vDeltaPhiNeg.size(), vPavg.data(), vDeltaPhiNeg.data(), zeroes.data(), zeroes.data());
    
    grDeltaPtotPos1->SetMarkerColor(kRed);
    grDeltaPtotNeg1->SetMarkerColor(kBlue);

    grDeltaPhiPos1->SetMarkerColor(kRed);
    grDeltaPhiNeg1->SetMarkerColor(kBlue);
    auto fstScaleDeltaPhi = grDeltaPhiPos1->GetYaxis();
    auto sndScaleDeltaPhi = new TGaxis(270, -0.003, 270, 0.001, -0.030, 0.01, 50010, "+");
    sndScaleDeltaPhi->SetTitleOffset(1.2);
    sndScaleDeltaPhi->SetLabelOffset(0.06);
    sndScaleDeltaPhi->SetTitle("#DeltaM^{(CrAngle)}, #frac{MeV}{c^{2}}");
    // grDeltaPtotPos1->SetTitle("#pi^{+} - red, #pi^{-} - blue;Cell num;#DeltaP_{tr}, #frac{MeV}{c}");
    // grDeltaPhiPos1->SetTitle("#pi^{+} - red, #pi^{-} - blue;Cell num;#Delta#phi, rad");
    grDeltaPtotPos1->SetTitle("#pi^{+} - #pi^{-};P_{avg}, #frac{MeV}{c};#DeltaP_{#pi^{+}} - #DeltaP_{#pi^{-}}, #frac{MeV}{c}");
    grDeltaPhiPos1->SetTitle("#pi^{+} - #pi^{-};P_{avg}, #frac{MeV}{c};#Delta#phi_{#pi^{+}} - #Delta#phi_{#pi^{-}}, rad");
    // grDeltaPtotPos1->SetTitle("Red - #pi^{+}, Blue - #pi^{-}");
    grDeltaPhiPos1->DrawClone("AP");
    sndScaleDeltaPhi->DrawClone("same");
    // grDeltaPtotNeg1->DrawClone("P same");

    std::vector<Float_t> vM_Exp = {497.59, 497.578, 497.569, 497.59, 497.573, 497.556};
    std::vector<Float_t> vMerrExp = {0.0258545, 0.00895756, 0.00402005, 0.00247398, 0.00670656, 0.0283649};
    vE = {505, 508, 509, 510, 511, 514};
    TGraphErrors grRCNC_exp(vM_Exp.size(), vE.data(), vM_Exp.data(), zeroes.data(), vMerrExp.data());
    // grRCNC_exp.DrawClone("AP");


    // grDeltaPtotPos1.DrawClone("AP");
    return 0;
}