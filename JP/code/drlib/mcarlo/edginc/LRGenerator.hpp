//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LRGenerator.hpp
//
//   Description : Generates 'Likelihood Ratios'. Used for calculating 
//                 greeks within MC simulation
//
//   Author      : Mark A Robson
//
//   Date        : 24 June 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_LRGENERATOR_HPP
#define EDR_LRGENERATOR_HPP
#include "edginc/Sensitivity.hpp"
#include "edginc/OutputName.hpp"

DRLIB_BEGIN_NAMESPACE
/** Used for calculating greeks within MC simulation using 'Likelihood
    Ratios' methodology. (The split into 2 classes makes the implementation
    easier at the path generator level.) */
class MCARLO_DLL LRGenerator{
public:
    virtual ~LRGenerator();

    class MCARLO_DLL GreekCalculator{
    public:
        virtual ~GreekCalculator();
        /** Invoked after each path is priced in the MC simulation. pathIdx
            indicates the path index whilst price is the value returned by
            the product's payoff function */
        virtual void processPath(int pathIdx, double price) = 0;

        /** Invoked after the MC simulation is complete. It should
            store the result for the greek corresponding to this
            GreekCalculator. The pvFactor supplied is that from the
            method pvFromPaymentDate on IMCProduct. The same factor is
            used on the fair value. */
        virtual void storeResult(Results*   results,
                                 double     pvFactor) = 0;
    };
    typedef refCountPtr<GreekCalculator> GreekCalculatorSP;
    typedef vector<GreekCalculatorSP> GreekCalculatorArray;

    /** Returns a GreekCalculator object for the requested
        sensitivity. The object will calculate the given greek using
        the Likelihood Ratio methodology. The greek will be calculated
        with respect to either the names supplied in the Sensitivity
        or a list of names that need to be tweaked based upon those
        within the MultiFactor object. Null will be returned if the
        particular greek is not supported. The calcStdErr indicates whether
        the standard error should be calculated as well */
    virtual GreekCalculator* greekCalculator(const Sensitivity* sens,
                                             bool               calcStdErr) = 0;
};

DRLIB_END_NAMESPACE
#endif
