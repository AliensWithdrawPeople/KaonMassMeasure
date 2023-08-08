#ifndef HistContainer_h
#define HistContainer_h

#include <memory>
#include <map>
#include <iostream>
#include <optional>

#include "TFile.h"
#include "TObject.h"
#include "TH2D.h"
#include "TH1.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphErrors.h" 

class HistContainer
{
    std::map<std::string, TH1*> container{};
public:
    ~HistContainer();

    bool Add(std::string name, TH1* obj);
    bool Erase(std::string name);
    std::optional<TH1*> Get(std::string name);
    bool Save(std::string filename);
};

HistContainer::~HistContainer()
{
    for(const auto &[key, val] : container)
    { delete val; }
}

bool HistContainer::Add(std::string name, TH1* obj)
{
    if(container.find(name) != container.end())
    { std::cout << "Warning: Object with such name (" << name << ") already exists. The obj was not added." << std::endl; }
    auto [it, isInserted] = container.try_emplace(name, obj);
    return isInserted;
}

bool HistContainer::Erase(std::string name)
{
    auto removedElements = container.erase(name);
    return removedElements == 1? true : false;
}

std::optional<TH1*> HistContainer::Get(std::string name)
{
    if(container.find(name) == container.end())
    { 
        std::cout << "Warning: There is no object with name " << name << ". " << std::endl; 
        return std::nullopt;
    }
    return container[name];
}

bool HistContainer::Save(std::string filename)
{
    auto file = TFile::Open(filename.c_str(), "recreate");
    for(const auto &[name, val] : container)
    { val->Write(name.c_str()); }
    file->Close();
    delete file;
}

#endif