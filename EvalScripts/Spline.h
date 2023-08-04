#ifndef Spline_h
#define Spline_h

#include <optional>
#include <functional>
#include <stdexcept>
#include <iostream>

#include "TSpline.h"
#include "TFile.h"

/**
* @brief Spline class used for interpolating data.
* It is a TSpline3 class wrapper.
*/
class Spline
{
private:
    std::optional<TSpline3> deltaPsi_RecGenDiff;
    std::vector<double> dataX; 
    std::vector<double> dataY;

public:
    Spline(std::vector<double> xs, std::vector<double> ys, bool createSpline = true);
    Spline(): Spline({}, {}, false) {}

    /**
     * @brief Reads spline from .root file.
    */
    Spline(std::string inputFilename, std::string splineName = "deltaPsi_RecGenDiff");

    /**
     * @brief Creates spline from earlier passed data.
    */
    void CreateSpline();
    double operator()(double x) const;

    /**
     * @brief Save spline as TSpline3 to a \b filename .root file. File recreates.
     * @returns \b true if spline was saved and \b false otherwise.
    */
    bool Save(std::string filename);
};

Spline::Spline(std::vector<double> xs, std::vector<double> ys, bool createSpline)
{
    if(xs.size() != ys.size())
    { throw std::logic_error("Error during constructing Spline: xs doesn't have the same size as ys."); }

    deltaPsi_RecGenDiff = std::nullopt;
    dataX = xs;
    dataY = ys;

    if(createSpline)
    { CreateSpline(); }
}

Spline::Spline(std::string inputFilename, std::string splineName)
{
    dataX = {};
    dataY = {};
    
    auto file = TFile::Open(inputFilename.c_str());
    deltaPsi_RecGenDiff.emplace(TSpline3(*file->Get<TSpline3>("splineName")));
    file->Close();
    delete file;
}

void Spline::CreateSpline()
{
    if(deltaPsi_RecGenDiff.has_value())
    { std::cout << "Warning in Spline::CreateSpline(): Spline has been already created. \
                    It was recreated with x and y that you might not provided." << std::endl; }

    deltaPsi_RecGenDiff.emplace(TSpline3("ksdpsi_RecGenDiff_Spline", dataX.data(), dataY.data(), dataX.size()));
}

double Spline::operator()(double x) const
{ 
    if(!deltaPsi_RecGenDiff.has_value())
    { throw std::runtime_error("Error in Spline class in operator() call: spline was not created"); } 
    else
    { return deltaPsi_RecGenDiff.value().Eval(x); }
}

bool Spline::Save(std::string filename)
{
    if(deltaPsi_RecGenDiff.has_value())
    {
        auto file = TFile::Open(filename.c_str(), "recreate");
        deltaPsi_RecGenDiff.value().Write("deltaPsi_RecGenDiff");
        file->Close();
        delete file;
        return true;
    }
    
    return false;
}

#endif