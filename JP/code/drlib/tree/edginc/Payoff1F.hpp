//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Model1F.hpp
//
//   Description : One factor model parent class.
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : February 8, 2002
//
//----------------------------------------------------------------------------

#ifndef PAYOFF1F_HPP
#define PAYOFF1F_HPP
#include "edginc/config.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL Payoff1F{
public:
    virtual void preCalc(int step, int idx) {};
    virtual void PayoffAtMat(const double * s, int step, int bot, int top,
                             int pStart, int pEnd, double * const * price) = 0;
    virtual void PayoffBeforeMat(const double * s, int step, int bot, int top,
                                 int pStart, int pEnd, double * const * price) = 0;

    // returns true for option payoff, false otherwise
    virtual bool Positive() = 0;

    virtual void postCalc(const double * s, int step, int bot, int top,
                          int pStart, int pEnd, double * const * price){};
};

DRLIB_END_NAMESPACE
#endif
