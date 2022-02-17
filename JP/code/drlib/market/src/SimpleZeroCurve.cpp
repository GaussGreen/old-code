//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : SimpleZeroCurve.cpp
//
//   Description : Zero curve with dates/rate explicitly defined.
//
//   Author      : Richard Appleton
//
//   Date        : 27 June 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SimpleZeroCurve.hpp"
#include "edginc/Maths.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/Atomic.hpp"
#include <math.h>

DRLIB_BEGIN_NAMESPACE


SimpleZeroCurve::SimpleZeroCurve(
    const DateTime&           pBaseDate, 
    int                       pBasis,
    const DayCountConvention* pDcc,
    const string&             pInterpolation,
    const ExpiryArray&        pDates,
    const DoubleArray&        pValues)
: ZeroCurve(TYPE), baseDate(pBaseDate), basis(pBasis), 
  interpolation(pInterpolation), loBound(0), hiBound(0)
{
    static const string method = "SimpleZeroCurve::SimpleZeroCurve";

    if (pDates.size() != pValues.size())
    {
        string msg = Format::toString(
            "dates (%d) and rates (%d) array sizes differ",
            pDates.size(),
            pValues.size());
        throw ModelException(method, msg);
    }

    if (!CString::equalsIgnoreCase(interpolation, "L")
     && !CString::equalsIgnoreCase(interpolation, "Linear", 6))
    {
        string msg = Format::toString(
            "Invalid interpolation type [%s]", 
            interpolation.c_str());
        throw ModelException(method, msg);
    }

    dayCountConv = DayCountConventionConstSP(pDcc ? pDcc : new Actual365F());

    // insert values
    dates.reserve(pDates.size());
    rates.reserve(pValues.size());
    
    for (int i = 0 ; i < pDates.size() ; i++)
    {
        // TBD!! permit DF's as well as rates
        // if (df)
        //  addDiscountFactor(date,pValues[i])
        addRate(pDates[i]->toDate(baseDate), pValues[i]);
    }
}


SimpleZeroCurve::~SimpleZeroCurve() 
{
}


/** strip out the rates and dates */
CashFlowArraySP SimpleZeroCurve::getRatesAndDates() const
{
    return CashFlow::createCashFlowArray(dates, rates);
}


/** strip out the dates */
DateTimeArray SimpleZeroCurve::getDates()const
{
    return DateTimeArray(dates);
}


/** Convert this risk free zero curve to a risky zero curve */
void SimpleZeroCurve::CDSriskyZeroCurve(double recovery, CashFlowArray& defaultRates)
{
    static const string method = "SimpleZeroCurve::CDSriskyZeroCurve";
    try
    {
        double riskFree;
        DateTime date;
        double interpDefRate;
        for (int i=0; i < dates.size(); i++)
        {
            riskFree = rates[i];
            date = dates[i];
            interpDefRate = CashFlow::interpolate(defaultRates, date, true);
            // convert the risk free rate to a risky rate
            rates[i] = pow(1.+interpDefRate, 1.-recovery)*(1.+riskFree) - 1.;
        }   
    }
    catch (exception &e) 
    {
        throw ModelException(e, method);
    }
}


// Calculates discount factor for a date
double SimpleZeroCurve::discountFactor(const DateTime& date) const
{
    static const string method = "SimpleZeroCurve::discountFactor";

    try
    {
        double discount;
        if (date.equals(baseDate))
        {
            discount = 1.0;
        }
        else
        {
            double rate = zeroCouponRate(date);
            discount = RateConversion::rateToDiscount
                (rate, baseDate, date, dayCountConv.get(), basis);
        }
        return discount;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


// Calculates an interpolated rate from a ZCurve at some date
// with option on how to handle dates past end of curve - either
// go flat or linearly extrapolate
double SimpleZeroCurve::zeroCouponRate(const DateTime& date) const
{
    static const string method = "simpleZeroCurve::zeroCouponRate";
    bool extendFlat = true;

    try
    {
        int size = dates.size();

        if (dates.getLength() < 1)
        {
            throw ModelException(method, "no points in zero curve");
        }
        int    dt = date.getDate();
        double rate;

        if (size == 1)
        {
            // can size be different from getLength() ? Yes
            rate = rates[0];
        } 
        else if (dates[0].getDate() >= dt)
        {
            /* Do *flat* extrapolation only when going backwards. This
             * is done so that swaps which have payments before the beginning
             * of the stub zero curve will still value to par. This can happen
             * very easily if there are swaps with front stubs.
             * We still permit forward non-flat extrapolation.
             */
            rate = rates[0];
        } 
        else  if (extendFlat && endDate().getDate() <= dt) 
        {
            // extrapolate flat off end of zero curve
            rate = rates[size-1];
        } 
        else 
        {
            int  lo;
            int  hi;
            const DateTime& loDate = dates[loBound];
            const DateTime& hiDate = dates[hiBound];

            if (dt >= loDate.getDate() && hiDate.getDate() >= dt) 
            {
                // we've already got the bounds
                lo = loBound;
                hi = hiBound;
            }
            else 
            {
                // Do a binary search to find the lower and upper bounds
                int mid;
                lo = 0;
                hi = dates.size() -1;
                
                while ((hi - lo) > 1) 
                {
                    mid = (hi + lo) >> 1;  // compute a mid point
                    if (dt >= dates[mid].getDate()) 
                    {
                        lo = mid;
                    }
                    else 
                    {
                        hi = mid;
                    }
                }
                
                // for next time
                loBound = lo;
                hiBound = hi;
            }
            
            if (dates[lo].getDate() == dt) 
            {
                rate = rates[lo];
            }
            else if (dates[hi].getDate() == dt) 
            {
                rate = rates[hi];
            }
            else 
            {
                rate = interpolate(date, dates[lo], rates[lo], dates[hi], rates[hi]);
            }
        }
        return rate;
    }
    catch (exception &e) 
    {
        throw ModelException(e, method);
    }
}


// Calculates an interpolated rate from a ZCurve at some date
double SimpleZeroCurve::interpolate(
    const DateTime& date,
    const DateTime& loDate, 
    double          loRate, 
    const DateTime& hiDate, 
    double          hiRate) const
{
    double rate;

    if (CString::equalsIgnoreCase(interpolation, "L")
     || CString::equalsIgnoreCase(interpolation, "Linear", 6))
    {
        double hi_lo = dayCountConv->years(loDate, hiDate);
        double dt_lo = dayCountConv->years(loDate, date);
        rate = loRate + ((hiRate - loRate)/hi_lo) * dt_lo;
    }
    else if (CString::equalsIgnoreCase(interpolation, "Flat Forward")
        || CString::equalsIgnoreCase(interpolation, "FlatForward")
        || CString::equalsIgnoreCase(interpolation, "Flat-Forward")
        || CString::equalsIgnoreCase(interpolation, "Flat_Forward"))
    {
        double dflo  = RateConversion::rateToDiscount(loRate, baseDate, loDate, dayCountConv.get(), basis);
        double dfhi  = RateConversion::rateToDiscount(hiRate, baseDate, hiDate, dayCountConv.get(), basis);
        double alpha = dayCountConv->years(loDate, date) / dayCountConv->years(loDate, hiDate);
        double df    = dflo * std::pow(dfhi/dflo, alpha);
        rate = RateConversion::discountToRate(df, baseDate, date, dayCountConv.get(), basis);
    }
    else
    {
        string msg = Format::toString(
            "Invalid interpolation type [%s]",
            interpolation.c_str());
        throw ModelException("SimpleZeroCurve::interpolate", msg);
    }

    return rate;
}


// what is the base date ?
const DateTime& SimpleZeroCurve::getBaseDate() const
{
    return baseDate;
}


// how long is the curve ?
int SimpleZeroCurve::length() const
{
    return dates.size();
}


// when is first date for which we have genuine information?
const DateTime& SimpleZeroCurve::firstDate() const
{
    return dates.empty() ? baseDate : dates[0];
}


// when does it end ?
const DateTime& SimpleZeroCurve::endDate() const
{
    return dates.back();
}


/* 'fast' linear interpolator. Assumes curve dates are in order and that
   supplied idx is optimum place to search from. 
   Returns an Annual Act/365F rate. No error checking of inputs */
double SimpleZeroCurve::fastInterpRate(
    int        date,       // (I) interp date 
    int*       idx) const  // (M) where to search from 
{   
    static const string routine("ZeroCurve::fastInterpRate");
    double       rate;
    int          size = dates.size();
 
    /* see if we're of the ends of the curve - if so, extrapolate flat */
    if (date <= dates[0].getDate())
    {
        rate  = rates[0];
    } 
    else if (date >= endDate().getDate()) 
    {
        rate = rates[size - 1];
    } 
    else  
    {
        int jdx = *idx;
        if (jdx >= size -1)
        {
            // idx is the index of the lower bound
            throw ModelException(routine, "Internal error");
        }
        /* search for bounding points (if necessary) */
        if (date > dates[jdx+1].getDate()) 
        {
            /* search forwards */
            for (jdx++; dates[jdx+1].getDate() < date; jdx++); /* empty */
            *idx = jdx;
        } 
        else if (date < dates[jdx].getDate()) 
        {
            /* search backwards */
            for (jdx--; dates[jdx].getDate() >= date; jdx--); /* empty */
            *idx = jdx;
        } 

        rate = interpolate(DateTime(date,0), dates[jdx], rates[jdx], dates[jdx+1], rates[jdx+1]);
        /** 
         *  Note:   interpolate() returns rates with whatever day count and compounding as in 
         *          the zero curve, here we need Act/365F anually compounding rate
         */
	}

    // convert to Act/365F anually compounding rate
    double df  = RateConversion::rateToDiscount(rate, baseDate, DateTime(date,0), dayCountConv.get(), basis);
    Actual365F act365f;
    rate = RateConversion::discountToRate(df, baseDate, DateTime(date,0), &act365f, 1);

    return rate;
}


// Add rate at a specified date
void SimpleZeroCurve::addRate(const DateTime& date, double rate)
{
    static const string method = "SimpleZeroCurve::addRate";

    try
    {
        if (date < baseDate)
        {
            string msg = Format::toString(
                "Cannot add rate at %s before base date %s",
                date.toString().c_str(),
                baseDate.toString().c_str()
                );
            throw ModelException(method, msg);
        }

        dates.push_back(date);
        rates.push_back(rate);
    }
    catch(exception& e )
    {
        throw ModelException(e, method);
    }
}


void SimpleZeroCurve::insertValue(const DateTime& date, double value)
{
    static const string method = "SimpleZeroCurve::insertValue";

    if (dates.size() > 0 && endDate().isGreaterOrEqual(date))
    {                   
        // make sure date not already in list */
        bool found = false;
        int i = 0;
        
        while (i < dates.size() && !found)
        {
            // match on date only
            found = date.getDate() == dates[i].getDate();
            i++;
        }

        // don't want multiple entries on same date
        if (found && date.getTime() != dates[i-1].getTime())
        {
            string msg = Format::toString(
                "rate at date %s already exists - trying to add at %s",
                dates[i-1].toString().c_str(),
                date.toString().c_str());
            throw ModelException(method, msg);
        }

        if (found && !Maths::equals(value, rates[i-1]))
        {
            string msg = Format::toString(
                "rate at date %s already exists",
                date.toString().c_str());
            throw ModelException(method, msg);
        }

        if (found)
        {
            // we have 2 entries on the same date, but they agree on
            // time and value so avoid creating duplicate point in 
            // zero curve. 
            return;
        }
    }

    dates.push_back(date);
    rates.push_back(value);
}



class SimpleZeroCurve::LogOfDiscFactorKey: public YieldCurve::IKey
{
public:
    LogOfDiscFactorKey(const SimpleZeroCurve& zc): 
        // note: loDate and hiDate initialized to 0 - forces linear search
        // from base date. To do: add firstTime flag and do binary search
        // first time in
        zc(zc), loDate(0), loIdx(0), loRateTT(0),
        hiDate(0), hiIdx(0), hiRateTT(0)
        {
        }

    double years(const DateTime& startDate, const DateTime& endDate) const
    {
        return zc.dayCountConv->years(startDate, endDate);
    }

    /** Calculates the continuous forward rate between the
        two dates - we are trying to solve for r (want -r*(t2-t1))
        e^r(t2-t1) = (1+r2)^t2/(1+r1)^t1 
        or alternatively r(t2-t1) = t2*ln(1+r2) - t1*ln(1+r1) */
    double calc(const DateTime& loDateTime, const DateTime&  hiDateTime)
    {
        int newLoDate = loDateTime.getDate();
        int newHiDate = hiDateTime.getDate();

        /* see if we can quit early */
        if (loDateTime == hiDateTime) 
        {
            return 0.0;
        } 
        // to do: if first time in do binary search instead
        if (newHiDate == loDate) 
        {
            hiDate = newHiDate;
            hiIdx  = loIdx;
            hiRateTT = loRateTT;
            loDate = newLoDate;
            double rate = zc.fastInterpRate(newLoDate, &loIdx);
            loRateTT = -years(loDateTime, zc.baseDate) * log(1.0 + rate);
            //// why we care about rate, which depends on many other things like day count, compound?
            //// let's use df which is unambiguous
            //loRateTT = zc.fastInterpLnDF(newLoDate, &loIdx);
        } 
        else if (newLoDate == hiDate) 
        {
            loDate = newLoDate;
            loIdx  = hiIdx;
            loRateTT = hiRateTT;
            hiDate = newHiDate;
            double rate = zc.fastInterpRate(hiDate, &hiIdx);
            hiRateTT = -years(hiDateTime, zc.baseDate) * log(1.0 + rate);
        } 
        else 
        {
            if (newHiDate != hiDate) 
            {
                if (newHiDate < loDate) 
                {
                    hiIdx = loIdx;
                }
                hiDate = newHiDate;
                double rate = zc.fastInterpRate(newHiDate, &hiIdx);
                hiRateTT = -years(hiDateTime, zc.baseDate) * log(1.0 + rate);
            }
            if (newLoDate != loDate) 
            {
                if (newLoDate > hiDate) 
                {
                    loIdx = hiIdx;
                }
                loDate = newLoDate;
                double rate = zc.fastInterpRate(newLoDate, &loIdx);
                loRateTT = -years(loDateTime, zc.baseDate) * log(1.0 + rate);
            }
        }
        return (loRateTT-hiRateTT); // = -1* (hiRateTT - loRateTT)
    }

private:
    const SimpleZeroCurve& zc;
    int           loDate;    /* last loDate */
    int           loIdx;     /* index on zero curve for lo date */
    double        loRateTT;  /* rate from base date to lo date * year frac */
    int           hiDate;    /* last hiDate */
    int           hiIdx;     /* index on zero curve for hi date */
    double        hiRateTT;  /* rate from base date to hi date * year frac */
};


/** Returns a key used to optimize repeated calculations of
    forward rates */
YieldCurve::IKey* SimpleZeroCurve::logOfDiscFactorKey() const
{
    return new LogOfDiscFactorKey(*this);
}


/*
 * Reflection support.
 */

SimpleZeroCurve::SimpleZeroCurve()
  : ZeroCurve(TYPE), loBound(0), hiBound(0)
{
}


class SimpleZeroCurveHelper
{
public:
    /** Invoked when Class is 'loaded' - we only really need the reflection
        info to serialize the class - currently only used by UntweakableYC */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic();
        REGISTER(SimpleZeroCurve, clazz);
        SUPERCLASS(ZeroCurve);
        EMPTY_SHELL_METHOD(defaultZeroCurve);
        FIELD(baseDate, "base date");
        FIELD(basis, "basis");
//        FIELD(type, "(R)ates, (D)iscount factors");
        FIELD(interpolation, "interpolation type");
        FIELD       (dayCountConv, "day count convention");
        FIELD(dates, "zero curve dates");
        FIELD(rates, "zero curve rates");
    }

    static IObject* defaultZeroCurve()
    {
        return new SimpleZeroCurve();
    }
};

CClassConstSP const SimpleZeroCurve::TYPE = CClass::registerClassLoadMethod(
    "SimpleZeroCurve", typeid(SimpleZeroCurve), SimpleZeroCurveHelper::load);


DRLIB_END_NAMESPACE

