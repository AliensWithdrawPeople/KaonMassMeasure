#ifndef Painter_h
#define Painter_h

#include "HistContainer.h"

struct Range{
    double min;
    double max;
};

class Painter
{
private:
    HistContainer* container;
public:
    Painter(HistContainer* container);

    void Draw(std::string name, std::string option = "");
    void Draw(std::string name, std::string option, Range range);
    void Draw(std::string name, std::string option, Range range, Style_t markerStyle, Color_t color);
};

Painter::Painter(HistContainer* container): container(container)
{ }

void Painter::Draw(std::string name, std::string option)
{
    auto hist = (*container)[name];
    hist->DrawClone(option.c_str());
}

void Painter::Draw(std::string name, std::string option, Range range)
{
    auto hist = (*container)[name];
    hist->GetXaxis()->SetRangeUser(range.min, range.max);
    hist->DrawClone(option.c_str());
}

void Painter::Draw(std::string name, std::string option, Range range, Style_t markerStyle, Color_t color)
{
    auto hist = (*container)[name];
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->GetXaxis()->SetRangeUser(range.min, range.max);
    hist->DrawClone(option.c_str());
}


#endif