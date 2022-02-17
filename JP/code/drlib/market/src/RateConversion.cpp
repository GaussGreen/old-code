//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RateConversion.cpp
//
//   Description : rate conversion utilities
//
//   Author      : Andrew J Swain
//
//   Date        : 30 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"

DRLIB_BEGIN_NAMESPACE

// convert rate quoted as 'inBasis' to one quoted as 'outBasis'
double  RateConversion::rateConvert(
    const DateTime&           startDate,     // (I) 
    const DateTime&           endDate,       // (I) 
    double                    inRate,        // (I) 
    const DayCountConvention *inDcc,         // (I) 
    int                       inBasis,       // (I) 
    const DayCountConvention *outDcc,        // (I) 
    int                       outBasis)      // (I) 
{
    static const string method = "RateConversion::rateConvert";

    try {
        // don't convert if it's already what you want
        if (inBasis == outBasis && (inDcc->getClass() == outDcc->getClass())) {
            return inRate;
        }

        double disc = RateConversion::rateToDiscount(inRate,
                                                     startDate,
                                                     endDate,
                                                     inDcc,
                                                     inBasis);

        double rate = RateConversion::discountToRate(disc,
                                                     startDate,
                                                     endDate,
                                                     outDcc,
                                                     outBasis);
        return rate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


double RateConversion::discountToRate(
    double                    discount,      // (I) 
    const DateTime&           startDate,     // (I) 
    const DateTime&           endDate,       // (I) 
    const DayCountConvention *dcc,           // (I)  
    int                       basis)         // (I) 
{
    static const string method = "RateConversion::discountToRate";

    try {
        double rate;
        double yearfrac;

        if (discount < 0.0) {
            throw ModelException(method, "discount factor < 0.0");
        }
    
        if (startDate.getDate() == endDate.getDate()) {
            rate = 0.0;
        }
        else {
            if (discount == 0.0) {
                throw ModelException(method, "discount factor = 0.0");
            }

            if (!dcc) {
                throw ModelException(method, "Day count convention is NULL");
            }

            yearfrac = dcc->years(startDate, endDate);
            
            rate = RateConversion::discountToRateYearFrac(discount,
                                                          yearfrac,
                                                          basis);
        }

        return rate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


double RateConversion::rateToDiscount(
    double                    rate,          // (I) 
    const DateTime&           startDate,     // (I) 
    const DateTime&           endDate,       // (I) 
    const DayCountConvention *dcc,           // (I)
    int                       basis)         // (I)
{
    static const string method = "RateConversion::rateToDiscount";

    try {
        double discount;
        double yearfrac;

        if (!dcc) {
            throw ModelException(method, "Day count convention is NULL");
        }

        yearfrac = dcc->years(startDate, endDate);

        discount = rateToDiscountYearFrac(rate,
                                          yearfrac,
                                          basis);

        return discount;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


double RateConversion::discountToRateYearFrac(
    double  discount,           // (I) Discount factor 
    double  yearFraction,       // (I) Year fraction  
    int     basis)              // (I) Basis for the rate 
{
    static const string method = "RateConversion::discountToRateYearFrac";

    try {
        double rate;

        // very small discount factors are valid - so don't check for isZero
        if (discount <= 0.0)  
        {
            throw ModelException(method, "discount " + 
                                 Format::toString(discount) + " <= 0");
        }

        if (Maths::isZero(yearFraction))
        {
            throw ModelException(method, "year fraction " + 
                                 Format::toString(yearFraction) + " <= 0");
        }

        if (basis == CompoundBasis::SIMPLE) {
            rate = (1.0/discount - 1.0)/yearFraction;
        }
        else if (basis == CompoundBasis::CONTINUOUS) {
            rate = (-log(discount)/yearFraction);
        }
        else {
            rate = basis * (pow(discount, -1.0/(basis*yearFraction)) - 1.0);
        }

        return rate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double RateConversion::rateToDiscountYearFrac(
    double  rate,           // (I) 
    double  yearFraction,   // (I) Year fraction 
    int     basis)          // (I) Basis for the rate 
{
    static const string method = "RateConversion::rateToDiscountYearFrac";

    try {
        double discount;

        if (basis == CompoundBasis::SIMPLE) {
            double denom = 1.0 + rate * yearFraction;

            if (denom <= 0.0 || Maths::isZero(denom)) {
                throw ModelException(method, "invalid simple interest rate " + 
                                     Format::toString(rate) + " <= 0");
            }
            discount = 1.0 / denom;
        }
        else if (basis == CompoundBasis::CONTINUOUS) {
            discount = exp(-rate*yearFraction);
        }
        else {
            double tmp = 1.0 + rate / (double)basis;
            /* Since pow(x,y) is not defined when x < 0 and y is not an integer,
             * check before calling it.
             */
            if (tmp <= 0.0 || Maths::isZero(tmp)) {
                throw ModelException(method, "invalid rate " + 
                                     Format::toString(rate) + " <= 0");
            }
            discount = pow( tmp, -basis*yearFraction);           
        }
    
        return discount;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Convert a future to a forward rate (simple basis) */
double RateConversion::futureToForward(
    const DateTime&           today, 
    double                    futuresRate,   // (I) Annual futures rate from T1 to T2
    double                    volatility,    // (I) Vol of short term rate
    double                    shortTermRate, // (I) Annualized short term rate
    const DateTime&           startDate,     // (I) T1 of future
    const DateTime&           endDate,       // (I) T2 of future
    const DayCountConvention *dcc)           // (I) 
{
    static const string method = "RateConversion::futureToForward";
    try {
        if (!dcc) {
            throw ModelException(method, "Day count convention is NULL");
        }

        if (Maths::isNegative(volatility)) {
            throw ModelException(method, "ir vol (" + 
                                 Format::toString(volatility) + ") is < 0");
        }
        

        double t1 = dcc->years(today, startDate);
        double t2 = dcc->years(today, endDate);

        if (Maths::isZero(volatility)) {
            return futuresRate;
        }

        double volAdj = volatility * log (1.0 + shortTermRate);

        double fwd = exp(log(1.0+futuresRate) - 0.5*volAdj*volAdj * t1*t2)-1.0;
        return fwd;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Convert a forward rate to a future */
double RateConversion::forwardToFuture(
    const DateTime&           today, 
    double                    fwdRate,       // (I) 
    double                    volatility,    // (I) Vol of short term rate
    const DateTime&           startDate,     // (I) T1 of fwd
    const DateTime&           endDate,       // (I) T2 of fwd
    const DayCountConvention *dcc)           // (I) 
{
    static const string method = "RateConversion::forwardToFuture";
    try {
        double futRate;

        if (!dcc) {
            throw ModelException(method, "Day count convention is NULL");
        }

        if (Maths::isNegative(volatility)) {
            throw ModelException(method, "ir vol (" + 
                                 Format::toString(volatility) + ") is < 0");
        }
        
        if (Maths::isPositive(volatility)) {
            double   loYears   = dcc->years(today, startDate);
            double   hiYears   = dcc->years(today, endDate);

            if (fwdRate < -0.5) {
                throw ModelException(method, "fwd rate (" + 
                                     Format::toString(fwdRate) + 
                                     ") is < -0.5");
            }

            if (Maths::isPositive(loYears) && Maths::isPositive(hiYears)) {
                double a = -0.5*loYears*hiYears*volatility*volatility;
                double b = 1.0;
                double c = -log(1.0 + fwdRate);

                double sq = b*b - 4*a*c;
                if (Maths::isNegative(sq)) {
                    throw ModelException(method,
                                         "non-real solution to quadratic equation");
                }

                double x = (-b + sqrt(sq)) / (2.0 * a);

                futRate = exp(x) - 1.0;
            }
            else {
                // be generous here
                futRate = fwdRate;
            }
        }
        else {
            futRate = fwdRate;
        }

        return futRate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Convert forward rates to spot rates */
DoubleArraySP RateConversion::forwardsToSpots(
    const DateTime&           today,
    const DateTimeArray&      forwardDates,
    const DoubleArray&        forwardRates,
    const DayCountConvention* dcc)
{
    static const string method = "RateConversion::forwardsToSpots";
    try
    {
        int numFwds = forwardDates.size();
        DoubleArraySP spots = DoubleArraySP(new DoubleArray(numFwds));

        DateTime prevDate = today;
        double yearFrac;
        double rateTimesYearFrac;
        double totalYearFrac = 0.;
        double totalRateTimesYearFrac = 0.;

        //convert forwards into spots
        //spot(n) = (sum i=1 to n) fwd(i) * yf(i)
        //          -----------------------------
        //          (sum i=1 to n) yf(i)
        for (int i = 0; i < numFwds; i++)
        {
            DateTime thisDate = forwardDates[i];

            yearFrac = dcc->years(prevDate, thisDate);
            rateTimesYearFrac = yearFrac * forwardRates[i];

            //update state
            totalYearFrac += yearFrac;
            totalRateTimesYearFrac += rateTimesYearFrac;

            (*spots)[i] = totalRateTimesYearFrac / totalYearFrac;

            //update loop variables for next pass
            prevDate = thisDate;
        }
        return spots;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** Convert spot rates to forward rates */
DoubleArraySP RateConversion::spotsToForwards(
    const DateTime&           today,
    const DateTimeArray&      spotDates,
    const DoubleArray&        spotRates,
    const DayCountConvention* dcc)
{
    static const string method = "RateConversion::spotsToForwards";
    try
    {
        //spot(n) = (sum i=1 to n) fwd(i) * yf(i)
        //          -----------------------------
        //          (sum i=1 to n) yf(i)
        //gives
        //fwd(n) = [spot(n) * (sum i=1 to n) yf(i)] - [(sum i= 1 to n-1) fwd(i) * yf(i)]
        //         ---------------------------------------------------------------------
        //                                         yf(n)

        int numSpots = spotDates.size();
        DoubleArraySP forwards = DoubleArraySP(new DoubleArray(numSpots));

        DateTime prevDate = today;
        double yearFrac;
        double totalYearFrac = 0.;
        double totalRateTimesYearFrac = 0.;

        for (int i=0; i<numSpots; i++)
        {
            DateTime thisDate = spotDates[i];
            yearFrac  = dcc->years(prevDate,thisDate);
            totalYearFrac += yearFrac;

            //convert from spot rate to forward
            (*forwards)[i] = ((spotRates[i] * totalYearFrac) - totalRateTimesYearFrac) / yearFrac;

            //update the loop variables
            totalRateTimesYearFrac += (*forwards)[i] * yearFrac;
            prevDate = thisDate;
        }

        return forwards;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE

