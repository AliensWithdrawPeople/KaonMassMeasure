#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"

void delta_E_isr()
{
    TCanvas canvas("canv", "", 800, 600);
    std::vector<Float_t> zeroes(100, 0.0);
    std::vector<Float_t> RC= {0.115504, 0.0860349, 0.069423, 0.0598971, 0.0769491, 0.118673, 0.197142, 0.346072, 0.477132, 1.56063};
    std::vector<Float_t> energies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};

    TGraph grRC(RC.size(), energies.data(), RC.data());


    grRC.GetXaxis()->SetTitle("E, MeV");
    grRC.GetYaxis()->SetTitle("#DeltaE^{(RC)}, MeV");

    grRC.GetYaxis()->SetRangeUser(0, 1.7);
    
    grRC.DrawClone("AP");

    TLatex lat; 
    double x = 0.2;
    double y = 0.8;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    lat.DrawLatex(x, y, "CMD-3");

    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/delta_E_isr.pdf");
}