#ifndef Evaluator_h
#define Evaluator_h

#include "Tree.h"

#include <string>


class Evaluator
{
public:
    virtual double Eval() = 0;
    virtual void Draw(std::string key) = 0;
};

#endif