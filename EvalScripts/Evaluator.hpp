#ifndef Evaluator_h
#define Evaluator_h

#include "Tree.hpp"

#include <string>
#include <utility>


class Evaluator
{
public:
    virtual std::pair<double, double> Eval() = 0;
    virtual void Draw(std::string name) = 0;
};

#endif