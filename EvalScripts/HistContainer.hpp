#ifndef HistContainer_hpp
#define HistContainer_hpp

#include <memory>
#include <map>
#include <iostream>
#include <optional>
#include <stdexcept>

#include "TFile.h"
#include "TObject.h"
#include "TH2D.h"
#include "TH1.h"
#include "TH1D.h"
#include "TProfile.h"


class HistContainer
{
private:
    std::map<std::string, TH1*> container{};

public:
    ~HistContainer();


    bool Add(std::string name, TH1* obj);
    bool Erase(std::string name);
    TH1* operator[](std::string name);
    bool Save(std::string filename);
};

HistContainer::~HistContainer()
{
    for(const auto &[key, val] : container)
    { delete val; }
}

bool HistContainer::Add(std::string name, TH1* obj)
{
    auto [it, isInserted] = container.try_emplace(name, obj);
    if(!isInserted)
    { std::cout << "Warning: Object with such name (" << name << ") already exists. The obj was not added." << std::endl; }
    
    return isInserted;
}

bool HistContainer::Erase(std::string name)
{
    auto removedElements = container.erase(name);
    return removedElements == 1? true : false;
}

TH1* HistContainer::operator[](std::string name)
{
    try
    { return container.at(name); }
    catch(const std::out_of_range& e)
    { throw std::out_of_range("Error in HistContainer operator[]: There is no object with name " + name + "."); }
}

bool HistContainer::Save(std::string filename)
{
    auto file = TFile::Open(filename.c_str(), "recreate");
    for(const auto &[name, val] : container)
    { val->Write(name.c_str()); }
    file->Close();
    delete file;

    return true;
}

#endif