//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RateConversion.hpp
//
//   Description : rate conversion utilities
//
//   Author      : Andrew J Swain
//
//   Date        : 30 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef RATECONVERSION_HPP
#define RATECONVERSION_HPP

#include "edginc/DateTime.hpp"
#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL RateConversion {

public:

    // convert rate quoted as 'inBasis' to one quoted as 'outBasis'
    static double rateConvert(
        const DateTime&           startDate,     // (I) 
        const DateTime&           endDate,       // (I) 
        double                    inRate,        // (I) 
        const DayCountConvention *inDcc,         // (I) 
        int                       inBasis,       // (I) 
        const DayCountConvention *outDcc,        // (I) 
        int                       outBasis);     // (I) 

    static double discountToRate(
        double                    discount,      // (I) Discount factor 
        const DateTime&           startDate,     // (I) 
        const DateTime&           endDate,       // (I) 
        const DayCountConvention *dcc,           // (I) 
        int                       basis);        // (I)

    static double rateToDiscount(
        double                    rate,          // (I) 
        const DateTime&           startDate,     // (I) 
        const DateTime&           endDate,       // (I) 
        const DayCountConvention *dcc,           // (I) 
        int                       basis);        // (I) 

    static double discountToRateYearFrac(
        double  discount,           // (I) Discount factor 
        double  yearFraction,       // (I) Year fraction  
        int     basis);             // (I) Basis for the rate 

    static double rateToDiscountYearFrac(
        double  rate,           // (I) 
        double  yearFraction,   // (I) Year fraction 
        int     basis);         // (I) Basis for the rate 

    /** Convert a future to a forward rate (simple basis) */
    static double futureToForward(
        const DateTime&           today, 
        double                    futuresRate,   // (I) Annual futures rate from T1 to T2
        double                    volatility,    // (I) Vol of short term rate
        double                    shortTermRate, // (I) Annualized short term rate
        const DateTime&           startDate,     // (I) T1 of future
        const DateTime&           endDate,       // (I) T2 of future
        const DayCountConvention *dcc);          // (I) 

    /** Convert a forward rate to a future */
    static double forwardToFuture(
        const DateTime&           today, 
        double                    fwdRate,       // (I) 
        double                    volatility,    // (I) Vol of short term rate
        const DateTime&           startDate,     // (I) T1 of fwd
        const DateTime&           endDate,       // (I) T2 of fwd
        const DayCountConvention *dcc);          // (I)         

    /** Convert forward rates to spot rates */
    static DoubleArraySP forwardsToSpots(
        const DateTime&           today,
        const DateTimeArray&      forwardDates,
        const DoubleArray&        forwardRates,
        const DayCountConvention* dcc);

    /** Convert spot rates to forward rates */
    static DoubleArraySP spotsToForwards(
        const DateTime&           today,
        const DateTimeArray&      spotDates,
        const DoubleArray&        spotRates,
        const DayCountConvention* dcc);

private:
    RateConversion();
};

DRLIB_END_NAMESPACE

#endif


