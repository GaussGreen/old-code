//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedBSIR.hpp
//
//   Description : What a VolatilityBS can do for IR Vols
//
//   Author      : Mark A Robson
//
//   Date        : 24 Jun 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_PROCESSED_BS_IR_HPP
#define EDR_VOL_PROCESSED_BS_IR_HPP
#include "edginc/VolProcessedBS.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/IRGridPointAbs.hpp"

DRLIB_BEGIN_NAMESPACE
class IRVolBase;
class YieldCurve;

/** Defines what a Black-Scholes processed volatility can do for IR
    Vols. Can be thought of as a Vol Curve */
class MARKET_DLL VolProcessedBSIR: public CVolProcessedBS{
public:
    static CClassConstSP const TYPE; // defined in VolProcessedBS.cpp

    /** Get benchmark details */
    virtual void getBMDetails(DateTimeArray& swaptionExpiries,
                              DateTimeArray& swapStarts,
                              DateTimeArray& swapMats,
                              DoubleArray&   swaptionVols) const = 0;

    /** Get benchmark details */
    virtual void getSwaptionGrid(
        DateTimeArray& expiries,     // in years, offset from today
        IntArray& tenors,            // in months = 0; 
        vector< vector<double> >& marketVols) const = 0; 

    /** Returns the day count convention of the underlying swaps */
    virtual DayCountConventionSP getSwapDCC() const = 0;

    /** Returns the frequency of the underlying swaps (eg 3M, 1A etc) */
    virtual MaturityPeriodSP getSwapFrequency() const = 0;

    /** ir vega sensitive points when models just do simple interpolation.
        In particular dates should be only those actually used */
    virtual IRGridPointAbsArraySP sensitiveIRVolPoints(
        const DateTimeArray& dates) const = 0;

    /** Returns relevant points on matrix as identified by getBMDetails. */
    virtual IRGridPointAbsArraySP sensitiveIRVolPoints() const = 0;
                      
    /** Wrapper around sensitiveIRVolPoints(const DateTimeArray&) above */
    static IRGridPointAbsArraySP sensitiveIRVolPoints(
        const IRVolBase*     irVol,
        const YieldCurve*    yc,
        CVolRequest*         volRequest,
        const DateTimeArray& dates); // in VolProcessedBS.cpp

protected:
    VolProcessedBSIR(const CClassConstSP& clazz);
private:
    VolProcessedBSIR(const VolProcessedBSIR &rhs);
    VolProcessedBSIR& operator=(const VolProcessedBSIR& rhs);
    
};

typedef smartConstPtr<VolProcessedBSIR> VolProcessedBSIRConstSP;
typedef smartPtr<VolProcessedBSIR> VolProcessedBSIRSP;

DRLIB_END_NAMESPACE
#endif
