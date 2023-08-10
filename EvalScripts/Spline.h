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
    std::optional<TSpline3> spline;
    std::string name;
    std::vector<double> dataX; 
    std::vector<double> dataY;

public:
    Spline(std::vector<double> xs, std::vector<double> ys, std::string name, bool createSpline = true);
    Spline(std::string name): Spline({}, {}, name, false) {}
    Spline(): Spline({}, {}, "spline", false) {}

    /// @brief Reads spline from .root file.
    Spline(std::string inputFilename, std::string splineName);

    /// @brief Creates spline from earlier passed data.
    void CreateSpline();
    double operator()(double x) const;

    /// @brief Save spline as TSpline3 to a filename .root file. File recreates
    /// @param filename Name of the file where spline will be saved,
    /// @param recreate Recreate file or not,
    /// @return true if spline was saved and false otherwise.
    bool Save(std::string filename, bool recreate);
};

Spline::Spline(std::vector<double> xs, std::vector<double> ys, std::string name, bool createSpline) : name(name)
{
    if(xs.size() != ys.size())
    { throw std::logic_error("Error during constructing Spline " + name + ": xs is not of the same size as ys."); }

    spline = std::nullopt;
    dataX = xs;
    dataY = ys;

    if(createSpline)
    { CreateSpline(); }
}

Spline::Spline(std::string inputFilename, std::string splineName): name(splineName)
{
    dataX = {};
    dataY = {};
    auto file = TFile::Open(inputFilename.c_str(), "read");
    spline.emplace(TSpline3(*file->Get<TSpline3>(name.c_str())));
    file->Close();
    delete file;
}

void Spline::CreateSpline()
{
    if(spline.has_value())
    { std::cout << "Warning in Spline::CreateSpline(): Spline has been already created. \
                    It was recreated with x and y that you might not provided." << std::endl; }

    spline.emplace(TSpline3(name.c_str(), dataX.data(), dataY.data(), dataX.size()));
}

double Spline::operator()(double x) const
{ 
    if(!spline.has_value())
    { throw std::runtime_error("Error in Spline class in operator() call: spline " + name + " was not created"); } 
    else
    { return spline.value().Eval(x); }
}

bool Spline::Save(std::string filename, bool recreate)
{
    if(spline.has_value())
    {
        auto file = TFile::Open(filename.c_str(), (recreate? "recreate" : "update"));
        spline.value().Write(name.c_str());
        file->Close();
        delete file;
        return true;
    }
    
    return false;
}

#endif