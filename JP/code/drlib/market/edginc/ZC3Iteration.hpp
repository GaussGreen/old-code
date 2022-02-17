//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3Iteration.hpp
//
//   Description : Defines iteration function objects.
//
//   Author      : Richard Appleton
//
//   Date        : 10th June 2005
//
//----------------------------------------------------------------------------

#ifndef ALIB_ITERATION_HPP
#define ALIB_ITERATION_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/RootFinder.hpp"


DRLIB_BEGIN_NAMESPACE

class ZeroCurve;
class IRVolBase;
class CashStream;
class ZC3ZeroCurve;
class ZC3CurveSegment;


class MARKET_DLL IterData : public Func1D::NoDeriv
{
public:
    static double LOWER_BOUND;
    static double UPPER_BOUND;

    double solve(const double guess) const;
    virtual double operator()(double  x) const;

protected:
    IterData(ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing);
    
    // template method functions for operator() implementation
    virtual void initIteration(double rate) const = 0;
    virtual double noDiscountCurvePv(double rate) const = 0;

    // SZCLinearAdjAbsolute, SZCFlatFwdAdjAbsolute, SZCSmoothFwdAdjAbsolute
    virtual void adjAbsolute() const = 0;

    ZC3ZeroCurve&     szc;
    const CashStream& cs;
    ZC3CurveSegment&  prevSegment;
    ZC3CurveSegment&  segment;
    double            pv;
    const ZeroCurve*  discountCurve;
    const IRVolBase*  volModelIR;

    // do we have any absolute adjustments in this segment?
    bool              useAdjAbsolute;

    // do we use smoothing during this operation?
    bool              useSmoothing;
};


class MARKET_DLL LinearObjFunc : public IterData
{
public:
    LinearObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing);
    ~LinearObjFunc() {}

protected:
    // pre-calculations
    double                     t0;
    double                     t1;    

private:
    virtual void adjAbsolute() const;
};


/**
 * Linear interpolation function for the case that the startRate at the 
 * beginning of the interval is known.
 */
class MARKET_DLL LinearObjFunc1 : public LinearObjFunc
{
public:
    LinearObjFunc1(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing);
    ~LinearObjFunc1() {}

private:
    void   initIteration(double rate) const;
    double noDiscountCurvePv(double rate) const;
};


/**
 * Linear interpolation function for the case that the startRate at the 
 * beginning of the interval is unknown, and hence we have a single zero rate
 * for the entire interval.
 */
class MARKET_DLL LinearObjFunc2 : public LinearObjFunc
{
public:
    LinearObjFunc2(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing);
    ~LinearObjFunc2() {}

private:
    void   initIteration(double rate) const;
    double noDiscountCurvePv(double rate) const;
};


class MARKET_DLL FlatFwdObjFunc : public IterData
{
public:
    FlatFwdObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             bool              useSmoothing);
    ~FlatFwdObjFunc() {}

private:
    void   adjAbsolute() const;
    void   initIteration(double rate) const;
    double noDiscountCurvePv(double rate) const;
};


// Smoothing iteration classes (GtoSZCSmoothFunc)

class MARKET_DLL SmoothFunc : public IterData
{
public:
    static const int NO_SMOOTH_DISCOUNT;

    SmoothFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            r0,
             double            r1);
    ~SmoothFunc() {}

protected:
    double t1;
    double r0;
    double r1;

private:
    friend class ZC3CurveSegment;

    /*
     * We need to match one or more absolute adjustments within the current
     * segment.
     *
     * To do this we will compute the relative adjustment implied by
     * the current choice of rate.
     *
     * This calculation is independent of whether we have a discount
     * curve or not.
     */
    void   adjAbsolute() const;

    double noDiscountCurvePv(double x) const;
};


/**
 * We constrain a0, a1, a2 so that we are flat at the beginning, i.e. a1=0.
 * We solve for a0 and this implies a2.
 */
class MARKET_DLL SmoothFirstObjFunc : public SmoothFunc
{
public:
    SmoothFirstObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            startFwdRate);
    ~SmoothFirstObjFunc() {}

private:
    void initIteration(double a0) const;
};


class MARKET_DLL SmoothParabObjFunc : public SmoothFunc
{
public:
    SmoothParabObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            startFwdRate);
    ~SmoothParabObjFunc() {}

private:
    void initIteration(double a1) const;
};


class MARKET_DLL SmoothLastObjFunc : public SmoothFunc
{
public:
    SmoothLastObjFunc(
             ZC3ZeroCurve&     curve, 
             const CashStream& cs,
             ZC3CurveSegment&  prevSegment,
             ZC3CurveSegment&  segment,
             double            pv,
             const ZeroCurve*  discountCurve,
             const IRVolBase*  volModelIR,
             bool              useAdjAbsolute,
             double            startFwdRate);
    ~SmoothLastObjFunc() {}

private:
    void initIteration(double a1) const;
};


/*
 * Interpolation function for adding shaped segments
 */
class MARKET_DLL ShapedFlatFwdObjFunc : public Func1D::NoDeriv
{
public:
    ShapedFlatFwdObjFunc(
             ZC3ZeroCurve&      curve, 
             const CashStream&  cs,
             const ZeroCurve&   discountCurve,
             const double       pv,
             const IRVolBase*   volModelIR,
             const BoolArray&   useAdjAbsolutes,
             const int          numSegs,
             const DoubleArray& shapeRates,
             ZC3CurveSegment&   prevSegment);
    ~ShapedFlatFwdObjFunc() {}

    double solve(const double guess) const;
    double operator()(double  x) const;

private:
    void initIteration(double rate) const;

    ZC3ZeroCurve&          curve;
    const CashStream&      cs;
    const ZeroCurve&       discountCurve;
    const double           pv;
    const IRVolBase*       volModelIR;
    const BoolArray&       useAdjAbsolutes;
    const int              numSegs;
    const int              startIndex;
    const DoubleArray&     shapeRates;
    const ZC3CurveSegment& prevSegment;
};


DRLIB_END_NAMESPACE

#endif // ALIB_ITERATION_HPP
