//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CTree1fPayoff.hpp
//
//   Description : base payoff class for tree1f.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef EDG_TREE1F_PAYOFF_H
#define EDG_TREE1F_PAYOFF_H

#include "edginc/Class.hpp"
#include "edginc/Payoff1F.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL TreePayoff1F: public Payoff1F{
public:
    virtual void preCalcTree(int step, int idx) {  
        preCalc(step, idx);   
    };
    
    virtual void PayoffAtMatTree(const double * s, int step, int bot, int top,
        int pStart, int pEnd, double * const * price) {
        PayoffAtMat(s, step, bot, top, pStart, pEnd, price);
    };
    
    virtual void PayoffBeforeMatTree(const double * s, int step, int bot, int top,
                                   int pStart, int pEnd, double * const * price) {
        PayoffBeforeMat(s, step, bot, top, pStart, pEnd, price);
    };

    virtual void postCalcTree(const double * s, int step, int bot, int top,
                              int pStart, int pEnd, double * const * price){
        postCalc(s, step, bot, top, pStart, pEnd, price);
    };
};


DRLIB_END_NAMESPACE
#endif
