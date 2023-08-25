#ifndef misc_hpp
#define misc_hpp

#include "Tree.hpp"

#include <functional>
#include <vector>
#include <numeric>

#include <TVector3.h>


namespace misc {
    enum EventType {cowboy, sailor};
    auto GetEventType = [](const track &piPos, const track &piNeg) {
        TVector3 momPos(1, 0, 0);
        TVector3 momNeg(1, 0, 0);
        TVector3 field(0., 0., 1.);
        momPos.SetMagThetaPhi(1, piPos.theta, piPos.phi);
        momNeg.SetMagThetaPhi(1, piNeg.theta, piNeg.phi);
        return momPos.Cross(field).XYvector().DeltaPhi(momNeg.XYvector()) < TMath::Pi() / 2?  cowboy : sailor;
    };

    std::string GetSplineFilename(std::string energyPoint)
    { return "C:/work/Science/BINP/Kaon Mass Measure/splines/spline_" + energyPoint + ".root"; }

    std::vector<int> pascalsTriangle(int n) 
    {
        std::vector<int> triangle(n + 1, 1);
        for (int i = 1; i < n; i++) 
        {
            for (int j = i; j >= 1; j--) 
            { triangle[j] += triangle[j - 1]; }
        }
        return triangle;
    }

    // Derivative with Richardson extrapolation.
    double derivative(std::function<double(double)> func, double x, double dx, int n) 
    {
        std::vector<int> binomial = pascalsTriangle(n);
        auto deriv = [&binomial, &func, &x, &n](double dx) {
            double result = 0.0;
            for (int i = 0; i <= n; i++) 
            { result += pow(-1, i) * binomial[i] * func(x + (n / 2. - i) * dx); }
            return result / pow(dx, n);
        };
        
        return (4 * deriv(dx / 2) - deriv(dx)) / 3;
    }

    double derivative2(std::function<double(double)> func, double x, double dx, int n) 
    {
        return (func(x + dx) - 2 * func(x) + func(x - dx)) / dx / dx;
    }

    double GetCorrectedEnergy(int rnum, double energy)
    {        
        // E = 508.5 MeV
        if(rnum >= 60196 && rnum <= 60259)
        { energy = 35.83 / (rnum - 60182.1) / (rnum - 60182.1) + 508.393; }

        if(rnum >= 60260 && rnum <= 60416)
        { energy = 769.991 / (rnum-60153.3) / (rnum-60153.3) + 508.361; }

        if(rnum >= 60417 && rnum <= 60498)
        { energy = 174994 / (rnum-59094.7) / (rnum-59094.7) + 508.288; }

        // E = 509 MeV
        if(rnum >= 60520 && rnum <= 60702)
        { energy = 0.0076 * sin(0.0950 * (rnum - 60584.6)) + 508.945; }

        // E = 509.5 MeV
        if(rnum >= 60790 && rnum <= 60921)
        { energy = 15.7437 / (rnum-60763.1) / (rnum-60763.1) + 509.518; }

        if(rnum >= 60922 && rnum <= 61174)
        { energy = 2441.92 / (rnum-60737.2) / (rnum-60737.2) + 509.497; }

        if(rnum >= 61175 && rnum <= 61378)
        { energy = 631.041 / (rnum-61094.5) / (rnum-61094.5) + 509.527; }
        
        // E = 510 MeV
        if(rnum >= 61380 && rnum <= 61461)
        { energy = 364.312 / (rnum - 61332.6) / (rnum - 61332.6) + 509.925; }

        if(rnum >= 61560 && rnum <= 61689)
        { energy = 598.053 / (rnum - 61479.6) / (rnum - 61479.6) + 509.928; }

        if(rnum >= 61689 && rnum <= 61856)
        { energy = 0.0101 * sin(0.0397 * (rnum-61704.6)) + 509.960; }

        // E = 510.5 MeV
        if(rnum >= 61859 && rnum <= 61958)
        { energy = 1587 / (rnum-61683) / (rnum-61683) + 510.434; }

        if(rnum >= 61958 && rnum <= 62075)
        { energy = -0.0109 * sin(0.068 * (rnum-62018)) + 510.448; }

        // E = 511 MeV
        if(rnum >= 62108 && rnum <= 62167)
        { energy = 50.2021 / (rnum-62077.4) / (rnum-62077.4) + 511.014; }

        if(rnum >= 62194 && rnum <= 62311)
        { energy = 14.23 / (rnum-62179) / (rnum-62179) + 511.031; }

        return energy;
    }
}

#endif