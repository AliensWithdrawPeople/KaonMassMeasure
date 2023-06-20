#include <iostream>
#include <vector>

#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TF1.h"


using std::vector;
using std::cout;
using std::endl;


/* It's a ROOT alternative to Shift class from EnergyShift.py.
*/
int Shift()
{
    vector<double> zeroes(100, 0.0);
    vector<double> vis = {497.677, 497.667, 497.622, 497.634, 497.615, 497.686, 497.78, 497.902, 499.057};
    vector<double> visErr = {0.036, 0.014, 0.011, 0.008, 0.007, 0.01, 0.013, 0.014, 0.029};
    vector<double> RC_data = {0.1, 0.099, 0.089, 0.08, 0.071, 0.116, 0.191, 0.334, 1.453};
    vector<double> E = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};

    auto Mvis = new TSpline3("Mvis", E.data(), vis.data(), vis.size());
    auto RC = new TSpline3("RC", E.data(), RC_data.data(), RC_data.size());

    auto fitFunc = new TF1("fitFunc", [&](double *x, double *par){ return par[0] + RC->Eval(x[0] + par[1]); }, 504, 515, 3);
    fitFunc->SetParameters(497.6, 0.0);

    TGraphErrors grMvis(vis.size(), E.data(), vis.data(), zeroes.data(), visErr.data());

    TFitResultPtr r = grMvis.Fit("fitFunc", "SME");
    cout << "M_RC = " << r->Parameter(0) << " +/- " << r->ParError(0) << endl;
    cout << "deltaE = " << r->Parameter(1) << " +/- " << r->ParError(1) << endl;
    return 0;
}