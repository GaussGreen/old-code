//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3Iteration.cpp
//
//   Description : Implementation for iteration function objects.
//
//   Author      : Richard Appleton
//
//   Date        : 10th June 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SecantBrentRootFinder.hpp"
#include "edginc/ZC3Iteration.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/Atomic.hpp"


DRLIB_BEGIN_NAMESPACE

double IterData::UPPER_BOUND =  10000.0;
double IterData::LOWER_BOUND = -UPPER_BOUND;

const int SmoothFunc::NO_SMOOTH_DISCOUNT = -1;


IterData::IterData(
         ZC3ZeroCurve&     curve, 
         const CashStream& pCs,
         ZC3CurveSegment&  pPrevSegment,
         ZC3CurveSegment&  pSegment,
         double            pPv,
         const ZeroCurve*  pDiscountCurve,
         const IRVolBase*  pVolModelIR,
         bool              pUseAdjAbsolute,
         bool              pUseSmoothing)
: szc(curve), cs(pCs), prevSegment(pPrevSegment), segment(pSegment), pv(pPv), 
  discountCurve(pDiscountCurve), volModelIR(pVolModelIR), useAdjAbsolute(pUseAdjAbsolute), useSmoothing(pUseSmoothing)
{
}


double IterData::solve(const double guess) const
{
    SecantBrentRootFinder rootfinder;
    return rootfinder.solve(*this, guess, LOWER_BOUND, UPPER_BOUND);
}


// ALIB: szcbuild.c#3505, 3385, 3239
double IterData::operator()(double  x) const
{
    static const string method = "IterData::operator()";

    initIteration(x);
    double result = pv;

    if (useAdjAbsolute)
    {
        adjAbsolute();
    }

    if (discountCurve)
    {
        for (int i = 0 ; i < cs.size() ; i++)
        {
            // calculate cash stream amount (accounting for smoothing)
            szc.setUseSmoothing(useSmoothing);
            double amount = cs[i].getAmount(szc, volModelIR);
            szc.setUseSmoothing(true);

            // get discount factor from discount curve
            double discount = discountCurve->discountFactor(cs[i].getPayDate());

            result -= discount * amount;
        }
    }
    else
    {
         result -= noDiscountCurvePv(x);
    }

    return result;
}


/*
 * Linear interpolation (GtoSZCLinearObjFunc).
 */
LinearObjFunc::LinearObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing)
  : IterData(curve, cs, prevSegment, segment, pv, discountCurve, volModelIR, useAdjAbsolute, useSmoothing)
{
    // Perform some pre-calculations to save on work for each iteration.
    t0 = curve.years(curve.getBaseDate(), prevSegment.date);
    t1 = curve.years(curve.getBaseDate(), segment.date);
}


// ALIB: szcbuild.c#3005 SZCLinearAdjAbsolute
void LinearObjFunc::adjAbsolute() const
{
    static const string method = "LinearObjFunc::adjAbsolute";

    /*
     * We need to match one or more absolute adjustments within the current
     * segment.
     *
     * To do this we will compute the relative adjustment implied by the
     * current choice of rate.
     *
     * This calculation is independent of whether we have a discount curve or
     * not.
     */
    double prevRate = 0.0;

    if (useSmoothing && prevSegment.isSmoothed() && !segment.isSmoothed())
    {
        prevRate = szc.discountToRate(prevSegment.smoothDiscount, t0);
    }
    else
    {
        prevRate = prevSegment.rate;
    }

    for (int i = 0 ; i < szc.adjAbsolute.size() ; i++)
    {
        ZC3Adjustment& adjustment = szc.adjAbsolute[i];
        const DateTime& startDate = adjustment.getStartDate();
        const DateTime& endDate = adjustment.getEndDate();

        if (adjustment.startsWithin(prevSegment.date, segment.date))
        {
            // lack of overlap was validated externally
            // use of this step validated externally

            double r0 = prevRate;                                  // rate at start of segment
            double r1 = segment.rate;                              // rate at end of segment
            double m = (r1 - r0) / (t1 - t0);                      // rate gradient
            double ts = szc.years(szc.getBaseDate(), startDate);
            double te = szc.years(szc.getBaseDate(), endDate);
            double rs = r0 + m * (ts - t0);                        // rate at adjustment start date
            double re = r0 + m * (te - t0);                        // rate at adjustment end date
            double ds = szc.rateToDiscount(rs, ts);                // discount at adjustment start date
            double de = szc.rateToDiscount(re, te);                // discount at adjustment end date
            double rt = log(de / (ds * adjustment.getDiscount())); // integral of rt over adjustment interval
            szc.modifyAdjustment(adjustment, rt);
        }
    }
}


LinearObjFunc1::LinearObjFunc1(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing)
  : LinearObjFunc(curve, cs, prevSegment, segment, pv, discountCurve, volModelIR, useAdjAbsolute, useSmoothing)
{
}


// ALIB: szcbuild.c#3256
void LinearObjFunc1::initIteration(double rate) const
{
    // setting segment changes curve prior to calling GtoCashStreamPointAmountUsingSZC
    segment.rate = rate;
}


// ALIB: szcbuild.c#3300
double LinearObjFunc1::noDiscountCurvePv(double rate) const
{
    double result = 0.0;
    double r0 = prevSegment.rate;   // rate at start of segment
    double r1 = segment.rate;       // rate at end of segment

    if (useSmoothing && prevSegment.isSmoothed() && !segment.isSmoothed())
    {
        r0 = szc.discountToRate(prevSegment.smoothDiscount, t0);
    }

    // we have precalculated t0 and t1 to save time here
    double m = (r1 - r0) / (t1 - t0); // rate gradient

    for (int i = 0 ; i < cs.size() ; i++)
    {
        // calculate cash stream amount without smoothing
        szc.setUseSmoothing(useSmoothing);
        double amount = cs[i].getAmount(szc, volModelIR);

        // calculate adjustment for this date from the base date
        double adjustment = szc.adjustmentFactor(szc.getBaseDate(), cs[i].getPayDate());
        szc.setUseSmoothing(true);

        /*
         * For linear interpolation we will calculate times from baseDate in the 
         * cash stream. This is because this is how the discount factors will be 
         * calculated, rather than relative to the previous point in the zero curve
         * We therefore also need to compute the pv at baseDate rather than
         * prevSegment->date.
         */
        double time = szc.years(szc.getBaseDate(), cs[i].getPayDate());

        // calculate discount factor from linearly interpolated rate
        double r = r0 + m * (time - t0);
        double discount = szc.rateToDiscount(r, time);

        result += discount * amount * adjustment;
    }

    return result;
}


LinearObjFunc2::LinearObjFunc2(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing)
  : LinearObjFunc(curve, cs, prevSegment, segment, pv, discountCurve, volModelIR, useAdjAbsolute, useSmoothing)
{
}


// ALIB: szcbuild.c#3404
void LinearObjFunc2::initIteration(double rate) const
{
    // setting segment and prevSegment changes curve prior to calling GtoCashStreamPointAmountUsingSZC
    segment.rate = rate;
    prevSegment.rate = rate;
}


// ALIB: szcbuild.c#3449
double LinearObjFunc2::noDiscountCurvePv(double rate) const
{
    double result = 0.0;
    
    for (int i = 0 ; i < cs.size() ; i++)
    {
        // calculate cash stream amount without smoothing
        szc.setUseSmoothing(useSmoothing);
        double amount = cs[i].getAmount(szc, volModelIR);

        // calculate adjustment for this date from the base date
        double adjustment = szc.adjustmentFactor(szc.getBaseDate(), cs[i].getPayDate());
        szc.setUseSmoothing(true);

        /*
         * For linear interpolation we will calculate times from baseDate in the 
         * cash stream. This is because this is how the discount factors will be 
         * calculated, rather than relative to the previous point in the zero curve
         * We therefore also need to compute the pv at baseDate rather than
         * prevSegment->date.
         */
        double time = szc.years(szc.getBaseDate(), cs[i].getPayDate());

        // calculate discount factor from linearly interpolated rate
        double discount = szc.rateToDiscount(rate, time);

        result += discount * amount * adjustment;
    }

    return result;
}


/*
 * Flat forwards interpolation (GtoSZCFlatFwdObjFunc).
 */
FlatFwdObjFunc::FlatFwdObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing)
  : IterData(curve, cs, prevSegment, segment, pv, discountCurve, volModelIR, useAdjAbsolute, useSmoothing)
{
}


// ALIB: szcbuild.c#3523
void FlatFwdObjFunc::initIteration(double rate) const
{
    // setting segment changes curve prior to calling GtoCashStreamPointAmountUsingSZC
    segment.avgRate = rate;
}


// ALIB: szcbuild.c#3564
double FlatFwdObjFunc::noDiscountCurvePv(double rate) const
{
    double result = 0.0;

    for (int i = 0 ; i < cs.size() ; i++)
    {
        // calculate cash stream amount without smoothing
        szc.setUseSmoothing(useSmoothing);
        double amount = cs[i].getAmount(szc, volModelIR);

        // calculate adjustment for this date within this segment
        szc.setUseSmoothing(false);
        double adjustment = szc.adjustmentFactor(prevSegment.date, cs[i].getPayDate());
        szc.setUseSmoothing(true);

        /*
         * For flat forwards interpolation we will calculate times from
         * prevSegment->date in the cash stream. This is because we are simply 
         * trying to calculate the average rate since the last time point in the 
         * zero curve. Thus within iterData we can use the pv of the remaining 
         * cash flows at the previous time point in the zero curve.
         */
        double time = szc.years(prevSegment.date, cs[i].getPayDate());

        // calculate discount factor from flat forward rate
        double discount = exp(-rate * time);

        result += discount * amount * adjustment;
    }

    return result;
}


void FlatFwdObjFunc::adjAbsolute() const
{
	// as functionality is also needed in ZC3RecurseInfo delegate to curve
	szc.flatFwdAdjAbsolute(prevSegment, segment);
}


/*
 * Smoothing support (GtoSZCSmoothFunc)
 */
SmoothFunc::SmoothFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            pR0,
             double            pR1)
  : IterData(curve, cs, prevSegment, segment, pv, discountCurve, volModelIR, useAdjAbsolute, true),
    r0(pR0), r1(pR1)
{
    t1 = curve.years(prevSegment.date, segment.date);
}

SmoothFirstObjFunc::SmoothFirstObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            startFwdRate)
  : SmoothFunc(curve, cs, prevSegment, segment, pv, discountCurve,
                     volModelIR, useAdjAbsolute, 0.0, startFwdRate)
{
}

SmoothLastObjFunc::SmoothLastObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            startFwdRate)
  : SmoothFunc(curve, cs, prevSegment, segment, pv, discountCurve, 
                     volModelIR, useAdjAbsolute, startFwdRate, 0.0)
{
}


SmoothParabObjFunc::SmoothParabObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            startFwdRate)
  : SmoothFunc(curve, cs, prevSegment, segment, pv, discountCurve, 
                     volModelIR, useAdjAbsolute, segment.a0, startFwdRate)
{
}


// ALIB:: szcbuild.c#3627
void SmoothFirstObjFunc::initIteration(double x) const
{
    // setting segment changes curve prior to calling GtoCashStreamPointAmountUsingSZC
    segment.smoothDiscount = NO_SMOOTH_DISCOUNT;
    segment.a0 = x;
    segment.a1 = 0.0;
    segment.a2 = (r1 - x) / (t1 * t1);
}


// ALIB:: szcbuild.c#3732
void SmoothLastObjFunc::initIteration(double x) const
{
    // setting segment changes curve prior to calling GtoCashStreamPointAmountUsingSZC
    segment.smoothDiscount = NO_SMOOTH_DISCOUNT;
    segment.a0 = r0;
    segment.a1 = x;
    segment.a2 = -0.5 * x / t1;
}


// ALIB:: szcbuild.c#3838
void SmoothParabObjFunc::initIteration(double x) const
{
    // setting segment changes curve prior to calling GtoCashStreamPointAmountUsingSZC
    segment.smoothDiscount = NO_SMOOTH_DISCOUNT;
    segment.a0 = r0;
    segment.a1 = x;
    segment.a2 = (r1 - r0 - x * t1) / (t1 * t1);
}


// ALIB: szcbuild.c#3163 (SZCSmoothFwdAdjAbsolute)
void SmoothFunc::adjAbsolute() const
{
    for (int i = 0 ; i < szc.adjAbsolute.size() ; i++)
    {
        ZC3Adjustment& adjustment = szc.adjAbsolute[i];

        if (adjustment.startsWithin(prevSegment.date, segment.date))
        {
            // lack of overlap validated externally
            // use of this step validated externally
            double ts = szc.years(prevSegment.date, adjustment.getStartDate());
            double te = szc.years(prevSegment.date, adjustment.getEndDate());
            double ds = segment.parabDiscount(1.0, ts);   // unadjusted discount at ts
            double de = segment.parabDiscount(1.0, te);   // unadjusted discount at te
            double rt = log(de / (ds * adjustment.getDiscount()));

            // these lines change the input curve
            // (in contrast to the unsmoothed case we do not modify adjustment.rt)
            szc.modifyAdjustment(adjustment, rt);
        }
    }
}


// ALIB: szcbuild.c#3670 / 3775 / 3881
double SmoothFunc::noDiscountCurvePv(double x) const
{
    double result = 0.0;

    for (int i = 0 ; i < cs.size() ; i++)
    {
        // calculate cash stream amount with smoothing
        szc.setUseSmoothing(useSmoothing);
        double amount = cs[i].getAmount(szc, volModelIR);
        szc.setUseSmoothing(true);

        // calculate adjustment for this date
        double adjustment = szc.adjustmentFactor(prevSegment.date, cs[i].getPayDate());

        // We will calculate times from prevSegment date in the cash stream.
        double time = szc.years(prevSegment.date, cs[i].getPayDate());

        // calculate discount factor from parabola
        double discount = segment.parabDiscount(1.0, time);
        
        result += discount * amount * adjustment;
    }

    return result;
}


ShapedFlatFwdObjFunc::ShapedFlatFwdObjFunc(
     ZC3ZeroCurve&      pCurve, 
     const CashStream&  pCs,
     const ZeroCurve&   pDiscountCurve,
     const double       pPv,
     const IRVolBase*   pVolModelIR,
     const BoolArray&   pUseAdjAbsolutes,
     const int          pNumSegs,
     const DoubleArray& pShapeRates,
     ZC3CurveSegment&   pPrevSegment)
  : curve(pCurve), cs(pCs), discountCurve(pDiscountCurve), pv(pPv), 
  volModelIR(pVolModelIR), useAdjAbsolutes(pUseAdjAbsolutes),
  numSegs(pNumSegs), startIndex(pCurve.data.size() - pNumSegs),
  shapeRates(pShapeRates), prevSegment(pPrevSegment)
{
}


double ShapedFlatFwdObjFunc::solve(const double guess) const
{
    SecantBrentRootFinder rootfinder;
    return rootfinder.solve(*this, guess, IterData::LOWER_BOUND, IterData::UPPER_BOUND);
}


double ShapedFlatFwdObjFunc::operator()(double  x) const
{
    double pvDiff = pv;

    // populate segments with shape and spread
    curve.populateShapedCurve(startIndex, numSegs, shapeRates, x, useAdjAbsolutes, false);

    for (int i = 0 ; i < cs.size() ; i++)
    {
        double amount = cs[i].getAmount(curve, volModelIR);
        double discount = discountCurve.discountFactor(cs[i].getPayDate());
        pvDiff -= discount * amount;
    }

    return pvDiff;
}


void ShapedFlatFwdObjFunc::initIteration(double rate) const
{
}


DRLIB_END_NAMESPACE
