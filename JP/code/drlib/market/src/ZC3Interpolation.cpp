//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3Interpolation.cpp
//
//   Description : Implementation for linear interpolation methods
//
//   Author      : Richard Appleton
//
//   Date        : 19th April 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/ZC3Interpolation.hpp"
#include "edginc/ZC3Iteration.hpp"
#include "edginc/ZC3RecurseInfo.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include <algorithm>


DRLIB_BEGIN_NAMESPACE

/*
 * Zero interpolation classes.
 */

const string ZC3ZeroInterpolation::LINEAR           = "Linear";
const string ZC3ZeroInterpolation::FLAT_FORWARDS    = "Flat";
const string ZC3ZeroInterpolation::SMOOTH_FORWARDS  = "Smooth";
const string ZC3ZeroInterpolation::CUBIC_SPLINE     = "Cubic Spline";
const string ZC3ZeroInterpolation::EXP_CUBIC_SPLINE = "Expoential Cubic Spline";


/* static */
bool ZC3ZeroInterpolation::lessThan(
    const ZC3ZeroInterpolationSP& first, 
    const ZC3ZeroInterpolationSP& second)
{
    // empty date = forever (which must be last segment!)
    if (first->date.empty())
    {
        return false;
    }
    else if (second->date.empty())
    {
        return true;
    }
    else
    {
        return first->date < second->date;
    }
}


/* static */
ZC3ZeroInterpolationConstSP ZC3ZeroInterpolation::make(
    const string&   type, 
    const DateTime& date,
    const bool      shaped)
{
    if (CString::equalsIgnoreCase(type, LINEAR, 1)
     || CString::equalsIgnoreCase(type, LINEAR, 6))
    {
        return ZC3ZeroInterpolationConstSP(new ZC3LinearZeroInterpolation(date));
    }

    if (CString::equalsIgnoreCase(type, FLAT_FORWARDS, 1)
     || CString::equalsIgnoreCase(type, FLAT_FORWARDS, 4))
    {
        return ZC3ZeroInterpolationConstSP(new ZC3FlatForwardsZeroInterpolation(date, true, shaped));
    }

    if (CString::equalsIgnoreCase(type, SMOOTH_FORWARDS, 1)
     || CString::equalsIgnoreCase(type, SMOOTH_FORWARDS, 6))
    {
        return ZC3ZeroInterpolationConstSP(new ZC3SmoothForwardsZeroInterpolation(date));
    }

    if (CString::equalsIgnoreCase(type, CUBIC_SPLINE, 1)
     || CString::equalsIgnoreCase(type, CUBIC_SPLINE, 5))
    {
        return ZC3ZeroInterpolationConstSP(new ZC3CubicSplineZeroInterpolation(date));
    }

    if (CString::equalsIgnoreCase(type, EXP_CUBIC_SPLINE, 1)
     || CString::equalsIgnoreCase(type, EXP_CUBIC_SPLINE, 3))
    {
        return ZC3ZeroInterpolationConstSP(new ZC3ExpCubicSplineZeroInterpolation(date));
    }

    string msg = Format::toString("Unknown interpolation type %s", type.c_str());
    throw ModelException("ZC3ZeroInterpolation::make", msg);
}


/* static */
bool ZC3ZeroInterpolation::validate(const string& interpolation)
{
    bool result = false;

    string::size_type pos = interpolation.find_last_of(';');
    if (pos != string::npos)
    {
        result = pos > 1
              && validate(interpolation.substr(pos + 1))
              && validate(interpolation.substr(0, pos));
    }
    else
    {
        try
        {
            pos = interpolation.find(',');
            if (pos != string::npos)
            {
                make(interpolation.substr(0, pos));
            }
            else
            {
                make(interpolation);
            }

            result = true;
        }
        catch(const ModelException&)
        {
            result = false;
        }
    }

    return result;
}


// ALIB: zcbuild3.c#242 (GtoStructuredInterpTypeNew)
/* static */
ZC3ZeroInterpolationArraySP ZC3ZeroInterpolation::createArray(
    const string&           spec,
    const DateTime&         baseDate,
    const Holiday&          holidays,
    const BadDayConvention& defaultBadDayConv)
{
    static const string method = "ZC3Interpolation::createArray";
    
    ZC3ZeroInterpolationArraySP results(new ZC3ZeroInterpolationArray());

    string::size_type pos1 = 0;
    string::size_type pos2 = 0;
    while (pos2 != string::npos)
    {
        pos2 = spec.find_first_of(';', pos1);
        string type = (pos2 == string::npos) ? spec.substr(pos1) 
                                             : spec.substr(pos1, pos2 - pos1);
        pos1 = pos2 + 1;

        string::size_type pos3 = type.find_first_of(',');

        // process (interval, bad day convention) [not for last segment]
        if (pos2 != string::npos && pos3 != string::npos)
        {
            BadDayConventionConstSP bdc(&defaultBadDayConv);
            string tenor = type.substr(pos3 + 1);
            type = type.substr(0, pos3);

            // check for day count convention override
            string::size_type pos4 = tenor.find_first_of(',');
            if (pos4 != string::npos)
            {
                string tmp = tenor.substr(pos4 + 1);
                bdc = BadDayConventionConstSP(BadDayConventionFactory::make(tmp));
                tenor = tenor.substr(0, pos4);
            }

            DateTime date = DateFwdThenAdjust(baseDate, tenor, 1, *bdc, holidays);

            // work around lack of arrays of smart const pointers
            ZC3ZeroInterpolationSP tmp(dynamic_cast<ZC3ZeroInterpolation*>(make(type,date)->clone()));
            results->push_back(tmp);
        }
        else
        {
            // ignore duration for the final section
            if (pos3 != string::npos)
            {
                type = type.substr(0, pos3);
            }

            // work around lack of arrays of smart const pointers
            ZC3ZeroInterpolationSP tmp(dynamic_cast<ZC3ZeroInterpolation*>(make(type)->clone()));
            results->push_back(tmp);
        }
    }

    // ensure array returned is sorted
    sort(results->begin(), results->end(), lessThan);

    /*
     * We have special processing for interpolation types which are used for
     * recursive builds. For these we check that there is only one recursive
     * section and that it is either the penultimate section with the last
     * section being flat or is the last section. 
     *
     * If the recursive section is the penultimate one we then collapse the
     * last two sections into a single section, where the benchmark flags are
     * changed to indicate which swaps are added in the final flat section.
     */
    // ALIB: zcrecurs.c#485 (GtoZCRecurseProcessStructuredInterp)
    int numSegments = results->size();

    // first we check that there is only one recursive section
    int numRecursive = 0;
    for (int j = 0 ; j < numSegments ; j++)
    {
        if ((*results)[j]->isRecursive())
        {
            numRecursive++;
        }
    }

    switch(numRecursive)
    {
    case 0: // nothing to do
        break;

    case 1:
        /*
         * We have only one recursive segment. It can either be the penultimate
         * or last section. If it is the last section there is nothing to do.
         */
        if ((*results)[numSegments-1]->isRecursive())
        {
            break;
        }

        if ((*results)[numSegments-2]->isRecursive())
        {
            if (!ZC3FlatForwardsZeroInterpolation::TYPE->isInstance(*(*results)[numSegments-1]))
            {
                string msg = "The interpolation following a recursive "
                             "section must be flat forwards";
                throw ModelException(method, msg);
            }
        }
        else
        {
            string msg = "Recursive interpolation types must be either "
                         "the last or the penultimate section";
            throw ModelException(method, msg);
        }

        break;

    default:
        throw ModelException(method, "Too many recursive sections in interpolation type");
    }

    if (results->size() > 1)
    {
        // one or two ALIB tests fail, so disable this feature
        throw ModelException(method, "Multiple interpolation types are not currently supported");
    }

    return results;
}


ZC3ZeroInterpolation::ZC3ZeroInterpolation(
    const DateTime&      pDate, 
	const bool           pFastInterp,
    const CClassConstSP& clazz)
  : CObject(clazz), date(pDate), fastInterp(pFastInterp)
{
}


ZC3ZeroInterpolation::~ZC3ZeroInterpolation()
{
}


// ALIB: szcswaps.c#5290 (SZCValidateInterpType)
void ZC3ZeroInterpolation::validateInterpType(
    const ZC3CouponInterpolation* couponInterp) const
{
}


const DateTime& ZC3ZeroInterpolation::until() const
{
    return date;
}


bool ZC3ZeroInterpolation::hasFastInterp() const
{
	return fastInterp;
}


bool ZC3ZeroInterpolation::isSmoothed() const
{
    return false;
}


bool ZC3ZeroInterpolation::isRecursive() const
{
    return false;
}


bool ZC3ZeroInterpolation::isShaped() const
{
    return false;
}


double ZC3ZeroInterpolation::discountToRate(
    const ZC3ZeroCurve& curve, 
    double              discount, 
    const DateTime&     start, 
    const DateTime&     end) const
{
    return curve.discountToRate(discount, start, end, true);
}


double ZC3ZeroInterpolation::smoothedDiscountFactor(
    const ZC3ZeroCurve&    zc3,
    const ZC3CurveSegment& segment,
    const DateTime&        startDate,
    double                 pvStart,
    const DateTime&        date
    ) const
{
    throw ModelException(getClass()->getName() + "::smoothedDiscountFactor", 
        "Program bug. Interp type cannot perform this operation");
}


/*
 * Linear interpolation.
 */

ZC3LinearZeroInterpolation::ZC3LinearZeroInterpolation(const DateTime& date)
  : ZC3ZeroInterpolation(date, true, TYPE)
{
}


ZC3LinearZeroInterpolation::~ZC3LinearZeroInterpolation()
{
}


// ALIB: zcrecurs.c#941
DoubleArraySP ZC3LinearZeroInterpolation::spline(
    const DateTimeArray& benchmarkDates, 
    const DateTimeArray& fittingDates, 
    const DoubleArray&   benchmarkInputs) const
{
    throw ModelException("ZC3LinearZeroInterpolation::spline",
        "Recursive fitting not available for this interp type.");
}


double ZC3LinearZeroInterpolation::interpolate(
    const ZC3ZeroCurve&    zc3, 
    const DateTime&        x, 
    const ZC3CurveSegment& p1, 
    const ZC3CurveSegment& p2) const
{
    double hi_lo = zc3.years(p1.date, p2.date);
    double dt_lo = zc3.years(p1.date, x);
    return p1.rate + ((p2.rate - p1.rate)/hi_lo) * dt_lo;
}


// ALIB: zcsmooth.c#3355 (InterpDiscountLinear)
double ZC3LinearZeroInterpolation::unsmoothedDiscountFactor(
        const ZC3ZeroCurve&    zc3,
        const ZC3CurveSegment& prevSegment,
        const ZC3CurveSegment& segment,
        const DateTime&        date,
        bool                   useSmoothing
        ) const
{
    const DateTime& startDate = prevSegment.date;
    const DateTime& endDate = segment.date;
    double startRate = prevSegment.rate;
    double endRate = segment.rate;

    /* 
     * No-one with any sense will ever enter the following block.  It is used
     * if there is a linear interpolated segment following a smooth segment.
     */
    if (useSmoothing && prevSegment.isSmoothed() && !segment.isSmoothed())
    {
        startRate = zc3.discountToRate(prevSegment.smoothDiscount, startDate);
    }

    /*
     * Linear interpolation.
     */
    double rate;
    double t1 = zc3.years(startDate, date);
    double t2 = zc3.years(startDate, endDate);

    if (t1 > t2)
    {
        /*
         * We are after the last date of the segment. Therefore we may be
         * in an extrapolation situation.
         */
        if (zc3.extrapDate <= endDate)
        {
            // No extrapolation allowed at all
            rate = endRate;
        }
        else
        {
            // Extrapolation occurs - either full or partial.
            if (date > zc3.extrapDate)
            {
                /*
                 * We extrapolate the endRate out to the extrapDate and then
                 * assume a flat zeroRate from that date. So recalculate t1
                 * and use that in the linear formula.
                 */
                t1 = zc3.years(startDate, zc3.extrapDate);
            }

            rate = startRate + (t1 / t2) * (endRate - startRate);
        }
    }
    else if (t1 < 0)
    {
        /*
         * Before the beginning of the segment - should be before the
         * start of the zero curve!
         */
        rate = startRate;
    }
    else
    {
        /*
         * The normal case - in the middle of a segment with simple
         * linear interpolation.
         */
        rate = startRate + (t1 / t2) * (endRate - startRate);
    }

    return zc3.rateToDiscount(rate,date);
}


// ALIB: szcbuild.c#1318
void ZC3LinearZeroInterpolation::solve(
    const CashStream& cs,
    ZC3ZeroCurve&     curve,
    const ZeroCurve*  discountCurve,
    const IRVolBase*  volModelIR,
    ZC3CurveSegment&  prevSegment,
    ZC3CurveSegment&  segment,
    double            presentValue,
    double            pvRemainder,
    double            pvKnown,
    bool              foundAdjAbsolute,
    bool              useSmoothing,
    double            prevDiscount
    ) const
{
    const static string method = "ZC3LinearZeroInterpolation::solve";

    if (prevSegment.date < curve.getBaseDate())
    {
        string msg = Format::toString(
            "Previous segment has date (%s) less than base date (%s)",
            prevSegment.date.toString().c_str(),
            curve.getBaseDate().toString().c_str());
        throw ModelException(method, msg);
    }
    else if (prevSegment.date == curve.getBaseDate())
    {
        /*
         * We are at the beginning of the curve, so we do not know the start
         * rate. Thus all rates are the same in this case. We use the second 
         * form of the linear objective function. Also we do not have a guess 
         * on the rate from the curve.
         */
        LinearObjFunc2 func(curve, cs, prevSegment, segment, presentValue - pvKnown, 
            discountCurve, volModelIR, foundAdjAbsolute, useSmoothing);
        segment.rate = func.solve(0.06);    // arbitrary 6% initial guess
    }
    else
    {
        /*
         * We are adding to an existing curve, so we know the start rate. Thus 
         * we need to provide the start rate to allow for linear interpolation 
         * in the objective function (we use the first form in this case).
         *
         * We will use the previous segment's rate as the initial guess.
         */

        LinearObjFunc1 func(curve, cs, prevSegment, segment, presentValue - pvKnown, 
            discountCurve, volModelIR, foundAdjAbsolute, useSmoothing);
        segment.rate = func.solve(prevSegment.rate);
    }

    /*
     * Having calculated the rate, we can now calculate the discount factor at 
     * the end of the segment.
     */
    segment.discount = curve.rateToDiscount(segment.rate, segment.date);

    /*
     * Finally we can calculate the average rate for this segment, since we 
     * have the total discount for this segment and the previous segment.
     */
    double discount = segment.discount / prevDiscount;
    segment.avgRate = curve.discountToRate(discount, prevSegment.date, segment.date, true);
}


/*
 * Flat forwards interpolation.
 */

ZC3FlatForwardsZeroInterpolation::ZC3FlatForwardsZeroInterpolation(
    const DateTime&      date, 
	const bool           pFastInterp,
    const bool           pShaped,
    const CClassConstSP& clazz)
  : ZC3ZeroInterpolation(date, pFastInterp, clazz), shaped(pShaped)
{
}


ZC3FlatForwardsZeroInterpolation::~ZC3FlatForwardsZeroInterpolation()
{
}


bool ZC3FlatForwardsZeroInterpolation::isShaped() const
{
    return shaped;
}


// ALIB: zcrecurs.c#941
DoubleArraySP ZC3FlatForwardsZeroInterpolation::spline(
    const DateTimeArray& benchmarkDates, 
    const DateTimeArray& fittingDates, 
    const DoubleArray&   benchmarkInputs) const
{
    throw ModelException("ZC3FlatForwardsZeroInterpolation::spline",
        "Recursive fitting not available for this interp type.");
}


/*
 * Flat forwards interpolation is linear interpolation in r*t when the rates
 * are continuously compounded.
 */
double ZC3FlatForwardsZeroInterpolation::interpolate(
    const ZC3ZeroCurve&    zc3, 
    const DateTime&        x, 
    const ZC3CurveSegment& p1, 
    const ZC3CurveSegment& p2) const
{
    const DateTime& dt = zc3.getBaseDate();
    const double r1 = p1.getContinuousRate();
    const double r2 = p2.getContinuousRate();
    const double t1 = zc3.years(p1.date, dt);
    const double t2 = zc3.years(p2.date, dt);
    double t = zc3.years(x, dt);

    if (t == 0)
    {
        // should never happen as have a segment at curve base date
        throw ModelException("ZC3FlatForwardsZeroInterpolation::interpolate","t == 0");
    }

    // Do the actual interpolation in r*t
    double rate = r1 * t1 + ((r2 * t2 - r1 * t1)/(t2-t1)) * (t-t1);
    rate /= t;
    
    // Convert the rate back from continuous to annual compounding
    return exp(rate) - 1;
}


// ALIB: zcsmooth.c#3246 (InterpDiscountFlatFwds)
double ZC3FlatForwardsZeroInterpolation::unsmoothedDiscountFactor(
        const ZC3ZeroCurve&    zc3,
        const ZC3CurveSegment& prevSegment,
        const ZC3CurveSegment& segment,
        const DateTime&        date,
        bool                   useSmoothing
        ) const
{
    double prevDiscount = (useSmoothing && prevSegment.isSmoothed() && !segment.isSmoothed()) 
        ? prevSegment.smoothDiscount : prevSegment.discount;

    /*
     * We can extrapolate out to zc3.extrapDate. For dates greater than
     * this date we compute the zeroRate at zc3.extrapDate and use that
     * instead.
     */
    if (date > segment.date && date > zc3.extrapDate)
    {
        double zeroRate;

        if (zc3.extrapDate <= segment.date)
        {
            /*
             * No extrapolation at all for this curve.
             */
            zeroRate = segment.rate;
        }
        else
        {
            /*
             * We are extrapolating the flat forward rate to extrapDate.  This
             * gives a discount factor for that date. We convert this back to
             * a zeroRate and use this zeroRate for the required date.
             */
            double avgRate = segment.avgRate;
            double t = zc3.years(prevSegment.date, zc3.extrapDate);
            double extrapDiscount = prevDiscount * exp(-avgRate * t);
            zeroRate = zc3.discountToRate(extrapDiscount,zc3.extrapDate);
        }

        /*
         * We have now calculated zeroRate by one of two methods.
         * Use it to calculate the discount factor to the desired date.
         */
        return zc3.rateToDiscount(zeroRate,date);
    }
    else
    {
        double avgRate = segment.avgRate;
        double t = zc3.years(prevSegment.date, date);
        return prevDiscount * exp(-avgRate * t);
    }
}


// ALIB: szcbuild.c#1449
void ZC3FlatForwardsZeroInterpolation::solve(
    const CashStream& cs,
    ZC3ZeroCurve&     curve,
    const ZeroCurve*  discountCurve,
    const IRVolBase*  volModelIR,
    ZC3CurveSegment&  prevSegment,
    ZC3CurveSegment&  segment,
    double            presentValue,
    double            pvRemainder,
    double            pvKnown,
    bool              foundAdjAbsolute,
    bool              useSmoothing,
    double            prevDiscount
    ) const
{
    const static string method = "ZC3FlatForwardsZeroInterpolation::solve";

    FlatFwdObjFunc func(curve, cs, prevSegment, segment, 
        pvRemainder, discountCurve, volModelIR, foundAdjAbsolute, useSmoothing);

    // use previous segment result or arbitrary 6% initial guess
    double guess = (curve.data.size() > 2) ? prevSegment.avgRate : 0.06;
    segment.avgRate = func.solve(guess);

    /*
     * Now that we have calculated the avgRate we can calculate the discount 
     * factor for this segment, and thus the cumulative discount factor since 
     * the base date. Then we can calculate the rate.
     */
    double discount = curve.rateToDiscount(segment.avgRate, prevSegment.date, segment.date, true);
    segment.discount = discount * prevDiscount;
    segment.rate = curve.discountToRate(segment.discount, segment.date);
}


/*
 * Smooth forwards interpolation.
 */

ZC3SmoothForwardsZeroInterpolation::ZC3SmoothForwardsZeroInterpolation(const DateTime& date)
  : ZC3FlatForwardsZeroInterpolation(date, false, false, TYPE)
{
}


ZC3SmoothForwardsZeroInterpolation::~ZC3SmoothForwardsZeroInterpolation()
{
}


// ALIB: szcswaps.c#5290 (SZCValidateInterpType)
void ZC3SmoothForwardsZeroInterpolation::validateInterpType(
    const ZC3CouponInterpolation* couponInterp) const
{
    if (couponInterp)
    {
        throw ModelException("Smooth forwards interpolation requires coupon interpolation of NONE");
    }
}


bool ZC3SmoothForwardsZeroInterpolation::isSmoothed() const
{
    return true;
}


// ALIB: zcsmooth.c#3486 [SmoothParabolicFwdDiscount]
double ZC3SmoothForwardsZeroInterpolation::smoothedDiscountFactor(
    const ZC3ZeroCurve&    szc,
    const ZC3CurveSegment& segment,
    const DateTime&        startDate,
    double                 startPv,
    const DateTime&        date
    ) const
{
    double discount = 0.0;

    if (date > segment.date)
    {
        /*
         * Extrapolation rules for the smooth zero curve are as follows:
         *
         * We calculate the discount factor and the continuously compounded
         * forward rate at segment->date.
         *
         * For the period out to szc->extrapDate we can use the forward rate
         * for discounting.
         *
         * For dates greater than szc->extrapDate, we compute the zeroRate at 
         * the extrapDate and use that for all subsequent dates we have a 
         * fixed zeroRate.
         */
        double t = szc.years(startDate, segment.date);
        double endDisc = segment.parabDiscount(startPv, t);                 // df at segment date
        double fwdRate = segment.a0 + segment.a1 * t + segment.a2 * t * t;  // fwd rate at segment date

        if (date <= szc.extrapDate)
        {
            t = szc.years(segment.date, date);
            discount = endDisc * exp(-fwdRate * t);
        }
        else
        {
            DateTime extrapDate = std::max(szc.extrapDate, segment.date);
            t = szc.years(segment.date, extrapDate);
            double extrapDisc = endDisc * exp(-fwdRate * t);
            double zeroRate = szc.discountToRate(extrapDisc, szc.getBaseDate(), extrapDate);
            discount = szc.rateToDiscount(zeroRate, szc.getBaseDate(), date);
        }
    }
    else
    {
        /*
         * Normal interpolation.
         */
        double t = szc.years(startDate, date);
        discount = segment.parabDiscount(startPv, t);
    }

    return discount;
}


/*
 * Common spline interpolation methods.
 */
ZC3SplineZeroInterpolation::ZC3SplineZeroInterpolation(
    const DateTime&      date, 
    const CClassConstSP& clazz)
  : ZC3FlatForwardsZeroInterpolation(date, false, false, clazz)
{
}


ZC3SplineZeroInterpolation::~ZC3SplineZeroInterpolation()
{
}


bool ZC3SplineZeroInterpolation::isRecursive() const
{
    return true;
}


// ALIB: sintrp1.inc#37
DoubleArraySP ZC3SplineZeroInterpolation::interpolate(
	const DateTimeArray& x,
	const DoubleArray&   f,
	const DateTimeArray& xDesired) const
{
	static const string method = "ZC3SplineZeroInterpolation::interpolate";

	if (x.size() != f.size())
	{
        string msg = Format::toString("Array sizes differ [x = %d, f = %d]", 
			x.size(), f.size());
		throw ModelException(method, msg);
	}

	if (x.size() < 1)
	{
        string msg = Format::toString("No array data [size = %d]", x.size());
		throw ModelException(method, msg);
	}

	int numInterps = xDesired.size();
	DoubleArraySP fInterp(new DoubleArray(numInterps));
	DoubleArraySP secondDeriv = secondDerivative(x, f);

	for (int i = 0 ; i < numInterps ; i++)
	{
		(*fInterp)[i] = interpolatePoint(*secondDeriv, x, f, xDesired[i]);
	}

	return fInterp;
}


// ALIB: sintrp1.inc#297
DoubleArraySP ZC3SplineZeroInterpolation::secondDerivative(
	const DateTimeArray& x,
	const DoubleArray&   f) const
{
	static const string method = "ZC3SplineZeroInterpolation::secondDerivative";

	double p, sig;
	double qn = 0;	// assume 2nd deriv = 0 at end
	double un = 0;	// assume 2nd deriv = 0 at end
	int N = x.size();

	// NB. validation done by caller

	DoubleArraySP f2(new DoubleArray(x.size()));
	DoubleArray temp(x.size());

	(*f2)[0] = 0.0;	// assume 2nd deriv = 0 at beginning
	temp[0] = 0.0;

	if (x[0] >= x[1])
	{
        string msg = Format::toString("x[0](%s) >= x[1](%s)", 
            x[0].toString().c_str(), x[1].toString().c_str());
		throw ModelException(method, msg);
	}

	for (int i = 1 ; i < N - 1 ; i++)
	{
		if (x[i] >= x[i+1])
		{
            string msg = Format::toString("x[%d](%s) >= x[%d](%s)", 
                i, x[i].toString().c_str(), i+1, x[i+1].toString().c_str());
			throw ModelException(method, msg);
		}

		sig = static_cast<double>(x[i].daysDiff(x[i-1])) / x[i+1].daysDiff(x[i-1]);
		p = sig * (*f2)[i-1] + 2.0;
		(*f2)[i] = (sig - 1.0) / p;

		// forward slope minus backwards slope
		temp[i] = (f[i+1] - f[i]) / x[i+1].daysDiff(x[i]) - (f[i] - f[i-1]) / x[i].daysDiff(x[i-1]);
		temp[i] = (6.0 * temp[i] / x[i+1].daysDiff(x[i-1]) - sig * temp[i-1]) / p;
	}

	(*f2)[N-1] = (un - qn * temp[N-2]) / (qn * (*f2)[N-2] + 1.0);

	for (int j = N-2 ; j >= 0 ; j--)
	{
		(*f2)[j] = (*f2)[j] * (*f2)[j+1] + temp[j];
	}

	return f2;
}


// ALIB: sintrp1.inc#214
double ZC3SplineZeroInterpolation::interpolatePoint(
	const DoubleArray&   secondDeriv,
	const DateTimeArray& x,
	const DoubleArray&   f,
	const DateTime&      xDesired) const
{
	static const string method = "ZC3SplineZeroInterpolation::interpolatePoint";

	if (xDesired < x[0] || xDesired > x[x.size() - 1])
	{
        string msg = Format::toString(
			"Cannot extrapolate: xDesired %s not in range [%s, %s]", 
			xDesired.toString().c_str(), x[0].toString().c_str(), 
            x[x.size() - 1].toString().c_str());
		throw ModelException(method, msg);
	}

	int lo = 0;
	int hi = 0;
	binarySearch(xDesired, x, lo, hi);

	double delta_X = x[hi].daysDiff(x[lo]);
    if (Maths::isZero(delta_X))	// happens if only 1 point
	{
		return f[lo];
	}

	double a = x[hi].daysDiff(xDesired) / delta_X;
	double b = 1.0 - a;

	// now calculate interpolated value
	return a * f[lo] + b * f[hi] + 
		((a*a*a-a) * secondDeriv[lo] + (b*b*b-b) * secondDeriv[hi]) 
		* delta_X * delta_X / 6;
}


// binary search function [bsearch.inc]
// TBD!! make generic templated binarySearch function, and use here
void ZC3SplineZeroInterpolation::binarySearch(
	const DateTime&      date, 
	const DateTimeArray& data, 
	int&                 lo, 
	int&                 hi) const
{
    static const string method = "ZC3SplineZeroInterpolation::binarySearch";

    size_t N = data.size();

    if (N < 2)
    {
        if (N < 1)
        {
            throw ModelException(method, "no data points");
        }
        else
        {
            // only 1 point - best guess is that point
            lo = 0;
            hi = 0;
            return;
        }
    }

    // extrapolate if desired date is before smallest in curve
    if (date <= data[0])
    {
        lo = 0;
        hi = 1;
        return;
    }

    // extrapolate if desired date is greater than biggest in curve
    if (date >= data[N-1])
    {
        lo = N - 2;
        hi = N - 1;
        return;
    }

    // do binary search to find pair of values which surround the desired segment
    int mid;
    size_t i;
    lo = 0;
    hi = N - 2;
    for (i = N ; i > 0 ; i--)
    {
        mid = (hi + lo) / 2;
        if (date < data[mid])
        {
            hi = mid - 1;
        }
        else if (date > data[mid+1])
        {
            lo = mid + 1;
        }
        else
        {
            break;
        }
    }

    if (i == 0)
    {
        throw ModelException(method, "Dates are not in increasing order");
    }

    // protect against a run of identical values
    lo = mid;
    hi = mid + 1;
    while (data[lo] == data[hi])
    {
        hi++;
    }
}


/*
 * Cubic spline interpolation.
 */

ZC3CubicSplineZeroInterpolation::ZC3CubicSplineZeroInterpolation(const DateTime& date)
  : ZC3SplineZeroInterpolation(date, TYPE)
{
}


ZC3CubicSplineZeroInterpolation::~ZC3CubicSplineZeroInterpolation()
{
}


// ALIB: zcrecurs.c#1204  [InterpolateSplineDF]
DoubleArraySP ZC3CubicSplineZeroInterpolation::spline(
    const DateTimeArray& benchmarkDates, 
    const DateTimeArray& fittingDates, 
    const DoubleArray&   benchmarkInputs) const
{
    return interpolate(benchmarkDates, benchmarkInputs, fittingDates);
}


/*
 * Exponential cubic spline interpolation.
 */

ZC3ExpCubicSplineZeroInterpolation::ZC3ExpCubicSplineZeroInterpolation(const DateTime& date)
  : ZC3SplineZeroInterpolation(date, TYPE)
{
}


ZC3ExpCubicSplineZeroInterpolation::~ZC3ExpCubicSplineZeroInterpolation()
{
}


// ALIB: zcrecurs.c#1134 [InterpolateSplineLogDF]
DoubleArraySP ZC3ExpCubicSplineZeroInterpolation::spline(
    const DateTimeArray& benchmarkDates, 
    const DateTimeArray& fittingDates, 
    const DoubleArray&   benchmarkInputs) const
{
	DoubleArray logInputs(benchmarkInputs.size());
	for (int i = 0 ; i < benchmarkInputs.size() ; i++)
	{
		logInputs[i] = log(benchmarkInputs[i]);
	}

	DoubleArraySP result = interpolate(benchmarkDates, logInputs, fittingDates);
	for(int j = 0 ; j < result->size() ; j++)
	{
		(*result)[j] = exp((*result)[j]);
	}

	return result;
}



/*
 * Coupon interpolation classes.
 */

const string ZC3CouponInterpolation::LINEAR     = "Linear";
const string ZC3CouponInterpolation::ANNUALIZED = "Annualized";
const string ZC3CouponInterpolation::NONE       = "None";


/* static */
bool ZC3CouponInterpolation::validate(const string& interpolation)
{
    try
    {
        make(interpolation);
        return true;
    }
    catch(const ModelException&)
    {
        return false;
    }
}


/* static */
ZC3CouponInterpolationConstSP ZC3CouponInterpolation::make(const string& type)
{
    if (CString::equalsIgnoreCase(type, LINEAR, 1)
     || CString::equalsIgnoreCase(type, LINEAR, 6))
    {
        return ZC3CouponInterpolationConstSP(new ZC3LinearCouponInterpolation(false));
    }

    if (CString::equalsIgnoreCase(type, ANNUALIZED, 1)
     || CString::equalsIgnoreCase(type, ANNUALIZED, 10))
    {
        return ZC3CouponInterpolationConstSP(new ZC3LinearCouponInterpolation(true));
    }

    if (CString::equalsIgnoreCase(type, NONE, 1)
     || CString::equalsIgnoreCase(type, NONE, 4))
    {
        return ZC3CouponInterpolationConstSP(   );
    }

    string msg = Format::toString("Unknown interpolation type %s", type.c_str());
    throw ModelException("ZC3CouponInterpolation::make", msg);
}


ZC3CouponInterpolation::ZC3CouponInterpolation(const CClassConstSP& clazz)
  : CObject(clazz)
{
}


ZC3CouponInterpolation::~ZC3CouponInterpolation()
{
}


ZC3LinearCouponInterpolation::ZC3LinearCouponInterpolation(bool pAnnualize)
: ZC3CouponInterpolation(TYPE), annualize(pAnnualize)
{
}


ZC3LinearCouponInterpolation::~ZC3LinearCouponInterpolation()
{
}


// ALIB: szcswaps.c#3143
ZC3SwapDataArraySP ZC3LinearCouponInterpolation::getInterpRates(
        const ZC3SwapDataArray& swapsData,
        const ZeroCurve&        discountCurve,
        const ZeroCurve*        estimatingCurve,
        const BadDayConvention& badDayConv,
        const StubPlacement&    stubPos,
        const Holiday&          holidays,
        const Holiday&          basisHolidays,
        const bool              convDelayAdj,
        const IRVolBase*        volModelIR,
        const bool              withBasis) const
{
    static const string method = "ZC3LinearCouponInterpolation::getInterpRates";

    return ZC3SwapData::getInterpRates(
        swapsData,
        discountCurve,
        estimatingCurve,
        badDayConv, 
        stubPos, 
        holidays, 
        basisHolidays, 
        convDelayAdj, 
        volModelIR,
        annualize,
        withBasis);
}


/*
 * Reflection support.
 */

static IObject* defaultCouponInterpolation()
{
    return new ZC3LinearCouponInterpolation(false);
}


/** Invoked when Class is 'loaded' */
void ZC3CouponInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3CouponInterpolation, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCouponInterpolation);
}


CClassConstSP const ZC3CouponInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3CouponInterpolation", typeid(ZC3CouponInterpolation), load);


static IObject* defaultLinearCouponInterpolation()
{
    return new ZC3LinearCouponInterpolation(false);
}


/** Invoked when Class is 'loaded' */
void ZC3LinearCouponInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3LinearCouponInterpolation, clazz);
    SUPERCLASS(ZC3CouponInterpolation);
    EMPTY_SHELL_METHOD(defaultLinearCouponInterpolation);
    FIELD_NO_DESC(annualize);
}


CClassConstSP const ZC3LinearCouponInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3LinearCouponInterpolation", typeid(ZC3LinearCouponInterpolation), load);



static IObject* defaultZeroInterpolation()
{
    return new ZC3FlatForwardsZeroInterpolation(DateTime());
}


/** Invoked when Class is 'loaded' */
void ZC3ZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3ZeroInterpolation, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultZeroInterpolation);

    //fields
    FIELD_NO_DESC(date);
	FIELD_NO_DESC(fastInterp);
	FIELD_MAKE_TRANSIENT(fastInterp);
}


CClassConstSP const ZC3ZeroInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3ZeroInterpolation", typeid(ZC3ZeroInterpolation), load);


DEFINE_TEMPLATE_TYPE(ZC3ZeroInterpolationArray);


static IObject* defaultLinearZeroInterpolation()
{
    return new ZC3LinearZeroInterpolation(DateTime());
}


/** Invoked when Class is 'loaded' */
void ZC3LinearZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3LinearZeroInterpolation, clazz);
    SUPERCLASS(ZC3ZeroInterpolation);
    EMPTY_SHELL_METHOD(defaultLinearZeroInterpolation);
}


CClassConstSP const ZC3LinearZeroInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3LinearZeroInterpolation", typeid(ZC3LinearZeroInterpolation), load);


static IObject* defaultFlatForwardsZeroInterpolation()
{
    return new ZC3FlatForwardsZeroInterpolation(DateTime());
}


/** Invoked when Class is 'loaded' */
void ZC3FlatForwardsZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3FlatForwardsZeroInterpolation, clazz);
    SUPERCLASS(ZC3ZeroInterpolation);
    EMPTY_SHELL_METHOD(defaultFlatForwardsZeroInterpolation);
}


CClassConstSP const ZC3FlatForwardsZeroInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3FlatForwardsZeroInterpolation", typeid(ZC3FlatForwardsZeroInterpolation), load);


static IObject* defaultSmoothForwardsZeroInterpolation()
{
    return new ZC3SmoothForwardsZeroInterpolation(DateTime());
}


/** Invoked when Class is 'loaded' */
void ZC3SmoothForwardsZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3SmoothForwardsZeroInterpolation, clazz);
    SUPERCLASS(ZC3FlatForwardsZeroInterpolation);
    EMPTY_SHELL_METHOD(defaultSmoothForwardsZeroInterpolation);
}


CClassConstSP const ZC3SmoothForwardsZeroInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3SmoothForwardsZeroInterpolation", typeid(ZC3SmoothForwardsZeroInterpolation), load);


/** Invoked when Class is 'loaded' */
void  ZC3SplineZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3SplineZeroInterpolation, clazz);
    SUPERCLASS(ZC3FlatForwardsZeroInterpolation);
    // no fields
}


CClassConstSP const  ZC3SplineZeroInterpolation::TYPE = CClass::registerClassLoadMethod(
    " ZC3SplineZeroInterpolation", typeid( ZC3SplineZeroInterpolation),  ZC3SplineZeroInterpolation::load);


static IObject* defaultCubicSplineZeroInterpolation()
{
    return new ZC3CubicSplineZeroInterpolation(DateTime());
}


/** Invoked when Class is 'loaded' */
void ZC3CubicSplineZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3CubicSplineZeroInterpolation, clazz);
    SUPERCLASS(ZC3SplineZeroInterpolation);
    EMPTY_SHELL_METHOD(defaultCubicSplineZeroInterpolation);
}


CClassConstSP const ZC3CubicSplineZeroInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3CubicSplineZeroInterpolation", typeid(ZC3CubicSplineZeroInterpolation), load);


static IObject* defaultExpCubicSplineZeroInterpolation()
{
    return new ZC3ExpCubicSplineZeroInterpolation(DateTime());
}


/** Invoked when Class is 'loaded' */
void ZC3ExpCubicSplineZeroInterpolation::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3ExpCubicSplineZeroInterpolation, clazz);
    SUPERCLASS(ZC3SplineZeroInterpolation);
    EMPTY_SHELL_METHOD(defaultExpCubicSplineZeroInterpolation);
}


CClassConstSP const ZC3ExpCubicSplineZeroInterpolation::TYPE = 
    CClass::registerClassLoadMethod("ZC3ExpCubicSplineZeroInterpolation", typeid(ZC3ExpCubicSplineZeroInterpolation), load);


DRLIB_END_NAMESPACE
