//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3CurveSegment.cpp
//
//   Description : Helper class for ALIB zero curve 3 bootstrapping method.
//
//   Author      : Richard Appleton
//
//   Date        : 13th May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3CurveSegment.hpp"
#include "edginc/ZC3Iteration.hpp"
#include "edginc/ZC3Interval.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/B30360.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/StubSimple.hpp"


DRLIB_BEGIN_NAMESPACE

static StubSimple stubType;


/* static */

bool ZC3CurveSegment::hasFastInterp(const ZC3CurveSegmentArray& segments)
{
	for (int i = 0 ; i < segments.size() ; i++)
	{
		if (!segments[i].interpolation->hasFastInterp())
			return false;
	}

	return true;
}

/* static */
double ZC3CurveSegment::unadjustedZeroCouponRate(
    const ZC3CurveSegmentArray& segments, 
    const ZC3ZeroCurve&         szc,
    const DateTime&             date,
    int&                        loBound,
    int&                        hiBound)
{
    int size = segments.size();

    if (size < 1)
    {
        throw ModelException("ZC3CurveSegment::unadjustedZeroCouponRate", 
            "no points in zero curve");
    }

    int    dt = date.getDate();
    double rate;

    if (size == 1)
    { // can size be different from getLength() ? Yes
        rate = segments.front().rate;
    } 
    else if (segments.front().date.getDate() >= dt) 
    {
        /* Do *flat* extrapolation only when going backwards. This
         * is done so that swaps which have payments before the beginning
         * of the stub zero curve will still value to par. This can happen
         * very easily if there are swaps with front stubs.
         * We still permit forward non-flat extrapolation.
         */
        rate = segments.front().rate;
    } 
    else  if (segments.back().date.getDate() <= dt) 
    {
        // extrapolate off end of zero curve
        rate = segments.back().rate;
    } 
    else 
    {
        int  lo;
        int  hi;
        const DateTime& loDate = segments[loBound].date;
        const DateTime& hiDate = segments[hiBound].date;

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
            hi = size -1;
            
            while ((hi - lo) > 1) 
            {
                mid = (hi + lo) >> 1;  // compute a mid point
                if (dt >= segments[mid].date.getDate()) 
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
        
        if (segments[lo].date.getDate() == dt) 
        {
            rate = segments[lo].rate;
        }
        else if (segments[hi].date.getDate() == dt) 
        {
            rate = segments[hi].rate;
        }
        else
        {
            return segments[lo].interpolation->interpolate
                (szc, date, segments[lo], segments[hi]);
        }
    }

    return rate;
}


/* static */
void ZC3CurveSegment::calculateAverageRates(
    ZC3CurveSegmentArray& segments, 
    const ZC3ZeroCurve&   szc,
    int                   numStubs)
{
    for (int j = 0 ; j < segments.size() ; j++)
    {
        segments[j].index = j;

        // We need to calculate the average rate for what was previously the
        // base date segment, which explains the condition.
        if (j <= numStubs)
        {
            segments[j].averageRate(szc, true);
        }
    }
}


/* static */
// ALIB: zcsmooth.c#1973 (GtoSZCDiscount)
double ZC3CurveSegment::discountFactor(
    const ZC3CurveSegmentArray& segments, 
    const ZC3ZeroCurve&         szc,
    const DateTime&             date, 
    bool                        useSmoothing)
{
    static const string method = "ZC3CurveSegment::discountFactor";

    /*
     * Special case is when we have only 1 segment, 
     * i.e. essentially the zero curve is empty.
     */
    if (segments.size() == 1)
    {
        if (szc.getBaseDate() != szc.firstDate())
        {
            throw ModelException(method, "Program bug! firstDate != baseDate");
        }

        return 1.0;
    }

    int lo;
    int hi;
    binarySearch(segments, lo, hi, date.getDate());

    double discount = 0.0;
    const ZC3CurveSegment& segment = segments[hi];
    const ZC3CurveSegment& prevSegment = segments[lo];

    /*
     * A little bit of logic to decide what calculations to perform.
     */
    if (!useSmoothing || !segment.isSmoothed())
    {
        /*
         * Calculate the unsmoothed discount - either because that was
         * requested, or it is the only one available for this segment.
         *
         * We can only support interpolation types which only rely on the begin
         * and end points of the segment.
         *
         * These are flat forwards and linear interpolation.
         */
        discount = segment.unsmoothedDiscountFactor
            (szc, prevSegment, date, useSmoothing);
    }
    else
    {
        /*
         * Calculate the smooth discount
         */
        discount = segment.smoothedDiscountFactor(szc, prevSegment, date);
    }

    return discount;
}


/* static */
// binary search function [bsearch.inc]
// TBD!! make generic templated binarySearch algorithm
void ZC3CurveSegment::binarySearch(
    const ZC3CurveSegmentArray& segments, 
    int&                        lo, 
    int&                        hi, 
    const int                   date)
{
    static const string method = "ZC3CurveSegment::binarySearch";

    size_t N = segments.size();

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
    if (date <= segments.front().date.getDate())
    {
        lo = 0;
        hi = 1;
        return;
    }

    // extrapolate if desired date is greater than biggest in curve
    if (date >= segments.back().date.getDate())
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
        if (date < segments[mid].date.getDate())
        {
            hi = mid - 1;
        }
        else if (date > segments[mid+1].date.getDate())
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
        throw ModelException(method, "Segments not in increasing order");
    }

    // protect against a run of identical values
    lo = mid;
    hi = mid + 1;
    while (segments[lo].date.getDate() == segments[hi].date.getDate())
    {
        hi++;
    }
}


ZC3CurveSegment::ZC3CurveSegment()
  :CObject(TYPE), interpolation(NULL), 
    discount(0.0), rate(0.0), smoothDiscount(0.0), avgRate(0.0),
    a0(0.0), a1(0.0), a2(0.0), duration(0.0), index(-1), contRate(-1)
{
}


ZC3CurveSegment::ZC3CurveSegment(
    const DateTime&             pDate,
    const ZC3ZeroInterpolation& pInterpolation,
    double                      pRate,
    double                      pDiscount,
    double                      pSmoothedDiscount,
    int                         pIndex)
  : CObject(TYPE), date(pDate), interpolation(&pInterpolation), 
    discount(pDiscount), rate(pRate), smoothDiscount(pSmoothedDiscount), avgRate(0.0),
    a0(0.0), a1(0.0), a2(0.0), duration(0.0), index(pIndex), contRate(-1)
{
}


ZC3CurveSegment::ZC3CurveSegment(
    const DateTime&             pDate,
    const ZC3ZeroInterpolation& pInterpolation,
    const ZC3CurveSegment&      prevSegment)
  : CObject(TYPE), date(pDate), interpolation(&pInterpolation), 
    discount(0.0), rate(0.0), smoothDiscount(0.0), avgRate(0.0), 
    a0(0.0), a1(0.0), a2(0.0), duration(0.0), index(prevSegment.index + 1), contRate(-1)
{
}


ZC3CurveSegment::ZC3CurveSegment(
    const DateTime&        discDate,
    const ZC3CurveSegment& prevSegment,
	double                 pDiscount,
	const ZC3ZeroCurve&    curve)
	: CObject(TYPE), date(discDate), interpolation(ZC3ZeroInterpolation::make("Flat")),
    discount(pDiscount), smoothDiscount(0.0), 
    a0(0.0), a1(0.0), a2(0.0), duration(0.0), index(prevSegment.index + 1), contRate(-1)
{
	// convert the discount to an average rate
	double time = curve.years(prevSegment.date, discDate);
	avgRate = -log(discount / prevSegment.discount) / time;
	rate = curve.discountToRate(discount, date);
}


bool ZC3CurveSegment::isSmoothed() const
{
    if (index == 0)
    {
        return false;
    }

    return interpolation.get() ? interpolation->isSmoothed() : false;
}


double ZC3CurveSegment::getDiscount() const
{
    return getDiscount(isSmoothed());
}


double ZC3CurveSegment::getDiscount(bool smoothed) const
{
    return smoothed ? smoothDiscount : discount;
}


double ZC3CurveSegment::getContinuousRate() const
{
    if (contRate == -1)
    {
        contRate = log(1 + rate);
    }

    return contRate;
}


double ZC3CurveSegment::unsmoothedDiscountFactor(
    const ZC3ZeroCurve&    szc, 
    const ZC3CurveSegment& prevSegment,
    const DateTime&        pDate,
    bool                   useSmoothing) const
{
    double discount = interpolation->unsmoothedDiscountFactor
        (szc, prevSegment, *this, pDate, useSmoothing);

    return discount;
}


double ZC3CurveSegment::smoothedDiscountFactor(
    const ZC3ZeroCurve&    szc, 
    const ZC3CurveSegment& prevSegment,
    const DateTime&        pDate) const
{
    DateTime startDate = prevSegment.date;
    double pvStart = prevSegment.getDiscount();

    double discount = interpolation->smoothedDiscountFactor
        (szc, *this, startDate, pvStart, pDate);

    return discount;
}


// ALIB: zcsmooth.c#3597 (GtoSZCSegmentAverageRate)
void ZC3CurveSegment::averageRate(
    const ZC3ZeroCurve& szc, 
    bool                forceRecalc)
{
    static const string method = "ZC3CurveSegment::averageRate";

    if (index == 0)
    {
        avgRate = 0.0; // meaningless quantity in this case
        return;
    }
    else if (!forceRecalc)
    {
        return;
    }
    else
    {
        const ZC3CurveSegment& previous = szc.data[index - 1];

        double t = szc.years(previous.date, date);
        if (!Maths::isPositive(t))
        {
            throw ModelException(method,"Non-increasing time points");
        }

        avgRate = log(previous.discount / discount) / t;
    }
}


// ALIB: szcbuild.c#2884 (SegmentDuration)
void ZC3CurveSegment::segmentDuration(const ZC3ZeroCurve& szc)
{
    static const string method = "ZC3CurveSegment::segmentDuration";

    if (duration > 0.0)
    {
        // do not recalculate it
        return;
    }

    B30360 dcc30360;
    double days30360 = dcc30360.days(szc.getBaseDate(), date);

    if (days30360 == 0)
    {
        duration = 0.0;
    }
    else
    {
        MaturityPeriod sixMonths("6M");
        BadDayNone none;
        HolidaySP holidays(Holiday::noHolidays());
        szc.setUseSmoothing(false);
        double parRate = SwapTool::parSwapRate(
            szc, 
            szc.getBaseDate(), 
            date, 
            sixMonths, 
            dcc30360, 
            stubType,
            false, 
            none,
            none,
            *holidays);
        szc.setUseSmoothing(true);

        Actual365F dccAct365F;
        double years = dccAct365F.years(szc.getBaseDate(), date);
        double freq = 2.0;

        duration = Maths::isZero(parRate) ? years 
            : (1.0 - pow(1.0 + parRate/freq, -years * freq)) / parRate;
    }
}


// ALIB: szcbuild.c#1861
void ZC3CurveSegment::smooth(
    ZC3ZeroCurve&    szc, 
    double           startFwdRate,
    const ZeroCurve* discountCurve,
    const IRVolBase* volModelIR)
{
    static const string method = "ZC3CurveSegment::smooth";

    /*
     * There are two possibilities.
     * 1. We are at the start of the curve. In this case, the segment was
     *    previously flat, and we need to re-calculate its start rate. In
     *    fact we go backwards from the end of this segment and force a
     *    segment which is flat at t=0.
     *
     * 2. We are not at the start of the curve. In this case, a0 is the
     *    start rate, and startFwdRate is the end rate for this segment.
     *
     * We must do this before acting on the final segment.
     */
    if (isSmoothed())
    {
        if (szc.data.size() <= 2)
        {
            /*
             * This is impossible, since if the array size is precisely 2 then
             * the previous segment must be the base date segment, which is
             * never smoothed.
             */
            string msg = Format::toString(
                              "Program bug. Array size must be > 2 [%s line %d]", 
                __FILE__, __LINE__);
            throw ModelException(method, msg);
        }

        ZC3CurveSegment& prevSegment = szc.data[index -1];
        DateTime startDate = prevSegment.date;
        bool atStart = (startDate == szc.getBaseDate());

        /*
         * Check to see if there any absolute adjustments in the interval
         *    [prevSegment2.date, prevSegment.date]
         * If they exist, then we will perform more complex calculations in
         * the objective functions.
         *
         * Also if they exist, it prevents closed form solutions.
         *
         * Finally if they exist, then they must be included wholly within
         * the segment.
         */
        bool foundAdjAbsPrev = szc.findAbsoluteAdjustments(prevSegment, *this);

        double t1 = szc.years(startDate, date);

        if (szc.cs->size() == 1 && (*szc.cs)[0].isFixed() && !foundAdjAbsPrev)
        {
            // no iteration required
            szc.setUseSmoothing(false);
            double fv = (*szc.cs)[0].getAmount(szc, volModelIR);
            szc.setUseSmoothing(true);
            double adjustment = szc.adjustmentFactor(prevSegment.date, (*szc.cs)[0].getPayDate());
            double discount = szc.PV / (fv * adjustment);

            if (atStart)
            {
                if (t1 > 0.0)
                {
                    a0 = -0.5 * (startFwdRate + 3.0 * log(discount) / t1);
                    a1 = 0.0;
                    a2 = (startFwdRate - a0) / (t1 * t1);
                }
                else
                {
                    // zero length segment so it doesn't matter what we do
                    a0 = startFwdRate;
                    a1 = 0.0;
                    a2 = 0.0;
                }
            }
            else
            {
                /**
                 * Calculate parabola coefficients from discount factor and 
                 * forward rates at the extremes.
                 */
                // ALIB: szcbuild.c#3922 (ParabCoeffs)
                // startRate = a0, endRate = startFwdRate
                a1 = -(6.0 * log(discount) + (a0 + a0 + startFwdRate) * 2.0 * t1) / (t1*t1);
                a2 = (startFwdRate - a0 - a1*t1) / (t1*t1);
            }
        }
        else
        {
            // iteration required
            if (atStart)
            {
                if (t1 > 0.0)
                {
                    // iterate to solve for a0
                    SmoothFirstObjFunc func(szc, *szc.cs, prevSegment, *this,
                        szc.PV, discountCurve, volModelIR, foundAdjAbsPrev, startFwdRate);
                    a0 = func.solve(prevSegment.a0);
                    a1 = 0.0;
                    a2 = (startFwdRate - a0) / (t1 * t1);
                }
                else
                {
                    // zero length segment so it doesn't matter what we do
                    a0 = startFwdRate;
                    a1 = 0.0;
                    a2 = 0.0;
                }
            }
            else
            {
                // iterate to solve for a1
                SmoothParabObjFunc func(szc, *szc.cs, prevSegment, *this,
                    szc.PV, discountCurve, volModelIR, foundAdjAbsPrev, startFwdRate);
                double guess = (startFwdRate - prevSegment.a0) / t1;
                a1 = func.solve(guess);
                a0 = func.r0;
                a2 = (func.r1 - func.r0 - a1 * t1) / (t1 * t1);
            }
        }

        smoothDiscount = parabDiscount(prevSegment.getDiscount(), t1);
    }
}


// ALIB: szcbuild.c#2164
void ZC3CurveSegment::recalculate(
    ZC3ZeroCurve&     szc, 
    ZC3CurveSegment&  prevSegment,
    const CashStream& cs,
    double            startFwdRate,
    double            presentValue,
    double            pvKnown,
    double            pvRemainder,
    bool              foundAdjAbs,
    const ZeroCurve*  discountCurve,
    const IRVolBase*  volModelIR)
{
    if (isSmoothed())
    {
        // perform smoothing
        double t1 = szc.years(prevSegment.date, date);
        a0 = startFwdRate;

        if (cs.size() == 1 && cs[0].isFixed() && !foundAdjAbs)
        {
            // closed form solution when we have single known cash flow
            szc.setUseSmoothing(false);
            double fv = cs[0].getAmount(szc, volModelIR);
            szc.setUseSmoothing(true);
            double adjustment = szc.adjustmentFactor(prevSegment.date, cs[0].getPayDate());
            a1 = 3.0 * (log(fv * adjustment/pvRemainder) - a0 * t1) / (t1 * t1);
        }
        else
        {
            // iteration required
            SmoothLastObjFunc func(szc, cs, prevSegment, *this, pvRemainder, discountCurve, volModelIR, foundAdjAbs, startFwdRate);
            a1 = func.solve(prevSegment.rate);
        }

        a2 = -0.5 * a1 / t1;
        smoothDiscount = parabDiscount(prevSegment.getDiscount(), t1);
    }
    else
    {
        /*
         * The rate, discount and avgRate for the segment must be recalculated
         * to take into account that the previous segment has now been smoothed.
         */
        if (cs.size() == 1 && cs[0].isFixed() && !foundAdjAbs)
        {
            // no iteration required
            szc.setUseSmoothing(false);
            double fv = cs[0].getAmount(szc, volModelIR);
            szc.setUseSmoothing(true);
            double pv = pvRemainder;//szc.pv;
            double adjustment = szc.adjustmentFactor(prevSegment.date, cs[0].getPayDate());
            double df = pv / (fv * adjustment);
            discount = prevSegment.getDiscount() * df;
            rate = szc.discountToRate(discount, date);
            avgRate = szc.discountToRate(df, prevSegment.date, date, true);
        }
        else
        {
            interpolation->solve(cs, szc, discountCurve, volModelIR, prevSegment, 
                *this, presentValue, pvRemainder, pvKnown, 
                foundAdjAbs, true, prevSegment.getDiscount());
        }
    }
}


// ALIB: szcbuild.c#3958 & zcsmooth.c#3664
double ZC3CurveSegment::parabDiscount(double prevDiscount, double time) const
{
    double rt = (a0 + (0.5 * a1 + a2 * time / 3.0) * time) * time;
    double discount = prevDiscount * exp(-rt);
    return discount;
}


/*
 * Reflection support
 */

static IObject* defaultSegment()
{
    return new ZC3CurveSegment();
}


/** Invoked when Class is 'loaded' */
void ZC3CurveSegment::load(CClassSP& clazz)
{
    clazz->setPrivate(); // make invisible to EAS/spreadsheet
    REGISTER(ZC3CurveSegment, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSegment);

    //fields
    FIELD(date,           "date at end of this segment");
    FIELD       (interpolation,  "simple interpolation type for this segment");
    FIELD(discount,       "unsmoothed discount factor at date");
    FIELD(rate,           "unsmoothed rate in basis and dcc");
    FIELD(smoothDiscount, "smoothed discount factor at date");
    FIELD(avgRate,        "unsmoothed continuously compounded rate for this segment");
    FIELD(a0,             "smoothing parameter");
    FIELD(a1,             "smoothing parameter");
    FIELD(a2,             "smoothing parameter");
    FIELD(duration,       "smoothing parameter");
    FIELD(index,          "position within data array");

    FIELD_NO_DESC(contRate);
    FIELD_MAKE_TRANSIENT(contRate);
}


CClassConstSP const ZC3CurveSegment::TYPE = 
    CClass::registerClassLoadMethod("ZC3CurveSegment", typeid(ZC3CurveSegment), load);


/**
 * Array support.
 */

DEFINE_TEMPLATE_TYPE(ZC3CurveSegmentArray);


/* static */
IObjectConstSP arrayObjectCast<ZC3CurveSegment>::toIObject(const ZC3CurveSegment& value)
{
    IObjectConstSP ptr(IObjectConstSP::attachToRef(&value));
    return ptr;
}


/* static */
IObjectSP arrayObjectCast<ZC3CurveSegment>::toIObject(ZC3CurveSegment& value)
{
    return IObjectSP::attachToRef(&value);
}


/* static */
ZC3CurveSegment arrayObjectCast<ZC3CurveSegment>::fromIObject(IObjectSP& value)
{
    ZC3CurveSegment* ptr = DYNAMIC_CAST(ZC3CurveSegment, value.get());
    if (!ptr)
    {
        throw ModelException("arrayObjectCast::fromIObject", "Object is not a ZC3CurveSegment");
    }

    return *ptr;
}


DRLIB_END_NAMESPACE
