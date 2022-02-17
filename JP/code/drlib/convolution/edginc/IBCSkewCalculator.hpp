//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 27-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IBCSKEWCALCULATOR_HPP
#define QLIB_IBCSKEWCALCULATOR_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Simple interface for objects capable to compute a Base Correlations adjusted
 * beta for a given set of (name historical beta, strike, time).
 * (See example of implementation in BCtrancheConvolutionLossGen.cpp).
 * */
class CONVOLUTION_DLL IBCSkewCalculator: public virtual VirtualDestructorBase
{
public:

    /**
     * Returns a Base Correlations adjusted beta for the given name historical beta, 
     * (adjusted) strike and time inputs.
     * */
    virtual double computeBCadjustedBeta(double nameHistBeta, double strike, const DateTime& time) const = 0;
};

DRLIB_END_NAMESPACE

#endif /*QLIB_IBCSKEWCALCULATOR_HPP*/
