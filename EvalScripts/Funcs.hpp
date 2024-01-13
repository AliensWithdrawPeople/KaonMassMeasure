#ifndef Funcs_hpp
#define Funcs_hpp

#include "misc.hpp"
#include "Tree.hpp"

#include <math.h>

namespace FullRecMassFunc
{
    // Mass of pi meson in MeV/c^2.
    constexpr double mass_pi = 139.57039; 

    /// @brief Variables of Ks mass function (full reconstruction method).
    enum Var {psi, energy, Y};

    /// @brief Evaluate Ks mass function
    /// @param dpsi Angle between two pions tracks in radians;
    /// @param energy Beam energy in MeV
    /// @param Y Ration of pions' momentums (Y = p_{pi^+} / p_{pi^-});
    /// @return Value of Ks mass function (full reconstruction method).
    double Eval(double dpsi, double energy, double Y)
    {
        auto E_2 = energy * energy; 
        auto eta_2 = std::pow((1 - Y * Y) / (1 + Y * Y), 2); 
        
        auto part1 = 1 / eta_2 * (1 + cos(dpsi) * sqrt(1 - eta_2));
        auto part2 = 1 - sqrt(1 - eta_2 * (1 - 4 * mass_pi * mass_pi / E_2));
        auto tot = E_2 * (1 - part1 * part2);
        return sqrt(tot);
    }

    /// @brief Evaluate numerical partial derivative of Ks Mass function with respect to var variable
    /// in point (dpsi, energy, Y) of corresponding order.
    /// @param var Evaluate derivative with respect to this variable;
    /// @param dpsi Angle between two pions tracks in radians;
    /// @param energy Beam energy in MeV
    /// @param Y Ration of pions' momentums (Y = p_{pi^+} / p_{pi^-});
    /// @param order Order of partial derivative.
    /// @return Numerical partial derivative of Ks Mass function (full reconstruction method).
    double Derivative(Var var, double dpsi, double energy, double Y, int order)
    {
        double res = 0;
        double dx = 1e-6;
        
        switch (var)
        {
        case Var::psi:
            res = misc::derivative( [&energy, &Y](double psi){ 
                return Eval(psi, energy, Y); 
                }, dpsi, dx, order);
            break;

        case Var::Y:
            res = misc::derivative( [&energy, &dpsi](double Y_val){ 
                return Eval(dpsi, energy, Y_val); 
                }, Y, dx, order);
            break;
        
        case Var::energy:
            res = misc::derivative( [&dpsi, &Y](double E){ 
                return Eval(dpsi, E, Y); 
                }, energy, dx, order);
            break;
        }   

        return res; 
    }
};

namespace PsiFunc
{
    /// @brief Variables of dpsi(thetaPos, phiPos, thetaNeg, phiNeg) where dpsi is angle between pions.
    enum Var {thetaPos, phiPos, thetaNeg, phiNeg};

    /// @param thetaPos theta angle of pi^+ momentum,
    /// @param phiPos phi angle of pi^+ momentum,
    /// @param thetaNeg theta angle of pi^- momentum,
    /// @param phiNeg phi angle of pi^- momentum.
    /// @return Angle between pi^+ and pi^-
    double Eval(double thetaPos, double phiPos, double thetaNeg, double phiNeg) 
    { return acos(sin(thetaPos) * sin(thetaNeg) * cos(phiPos - phiNeg) + cos(thetaPos) * cos(thetaNeg)); }
    
    double Derivative(Var var, double thetaPos, double phiPos, double thetaNeg, double phiNeg, int order)
    {
        double res = 0;
        double dx = 1e-6;

        switch (var)
        {
        case Var::thetaPos:
            res = misc::derivative( [&phiPos, &thetaNeg, &phiNeg](double thetaPos){ 
                    return Eval(thetaPos, phiPos, thetaNeg, phiNeg); 
                    }, thetaPos, dx, order);
            break;
        
        case Var::thetaNeg:
            res = misc::derivative( [&thetaPos, &phiPos, &phiNeg](double thetaNeg){ 
                    return Eval(thetaPos, phiPos, thetaNeg, phiNeg); 
                    }, thetaNeg, dx, order);
            break;

        case Var::phiPos:
            res = misc::derivative( [&thetaPos, &thetaNeg, &phiNeg](double phiPos){ 
                    return Eval(thetaPos, phiPos, thetaNeg, phiNeg); 
                    }, phiPos, dx, order);
            break;

        case Var::phiNeg:
            res = misc::derivative( [&thetaPos, &phiPos, &thetaNeg](double phiNeg){ 
                    return Eval(thetaPos, phiPos, thetaNeg, phiNeg); 
                    }, phiNeg, dx, order);
            break;
        }

        return res;
    }

    /// @brief Returns d^2(psi)/(d(thetaPos)d(thetaNeg)).
    double MixedThetaDerivative(double thetaPos, double phiPos, double thetaNeg, double phiNeg)
    {
        return (-cos(thetaPos)*cos(thetaNeg)*cos(phiPos - phiNeg) - sin(thetaPos) * sin(thetaNeg)) / 
        sqrt(1 - std::pow(sin(thetaPos) * sin(thetaNeg) * cos(phiPos - phiNeg) + cos(thetaPos) * cos(thetaNeg), 2)) + 
        (
            (sin(thetaPos) * cos(thetaNeg) * cos(phiPos - phiNeg) - cos(thetaPos) * sin(thetaNeg)) * 
            (sin(thetaPos) * cos(thetaNeg) - cos(thetaPos) * sin(thetaNeg) * cos(phiPos - phiNeg)) * 
            (sin(thetaPos) * sin(thetaNeg) * cos(phiPos - phiNeg) + cos(thetaPos) * cos(thetaNeg))
        ) / std::pow(1 - std::pow(sin(thetaPos) * sin(thetaNeg) * cos(phiPos - phiNeg) + cos(thetaPos) * cos(thetaNeg), 2), 3./2.);
    }

    /// @brief Returns d^2(psi)/(d(thetaPos)d(thetaNeg)).
    double MixedThetaDerivative(track piPos, track piNeg)
    { return MixedThetaDerivative(piPos.theta, piPos.phi, piNeg.theta, piNeg.phi); }

    /// @param piPos track data of pi^+
    /// @param piNeg track data of pi^-
    /// @return Angle between pi^+ and pi^-
    double Eval(track piPos, track piNeg)
    { return Eval(piPos.theta, piPos.phi, piNeg.theta, piNeg.phi); }
    
    double Derivative(Var var, track piPos, track piNeg, int order)
    { return Derivative(var, piPos.theta, piPos.phi, piNeg.theta, piNeg.phi, order); }
}

#endif