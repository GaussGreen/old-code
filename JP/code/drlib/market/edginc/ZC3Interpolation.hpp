//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3Interpolation.hpp
//
//   Description : Defines interface for curve interpolation methods
//
//   Author      : Richard Appleton
//
//   Date        : 19th April 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_INTERPOLATION_HPP
#define ZC3_INTERPOLATION_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ZC3Stub.hpp"
#include "edginc/ZC3CurveInstrument.hpp"


DRLIB_BEGIN_NAMESPACE

class BadDayConvention;
class Holiday;
class ZeroCurve;
class IRVolBase;
class CashStream;
class StubPlacement;
class ZC3ZeroCurve;
class ZC3CurveSegment;


class ZC3CouponInterpolation;
typedef smartConstPtr<ZC3CouponInterpolation>              ZC3CouponInterpolationConstSP;

class ZC3ZeroInterpolation;
typedef smartPtr<ZC3ZeroInterpolation>                     ZC3ZeroInterpolationSP;
typedef smartConstPtr<ZC3ZeroInterpolation>                ZC3ZeroInterpolationConstSP;
typedef array<ZC3ZeroInterpolationSP,ZC3ZeroInterpolation> ZC3ZeroInterpolationArray;    // really want array of const pointers
typedef smartPtr<ZC3ZeroInterpolationArray>                ZC3ZeroInterpolationArraySP;


/*
 * Coupon interpolation.
 */

class MARKET_DLL ZC3CouponInterpolation : public CObject
{
public:
    static CClassConstSP const TYPE;

    static const string LINEAR;
    static const string ANNUALIZED;
    static const string NONE;

    
    /**
     * Validate an interpolation format string.
     */
    static bool validate(const string& interpolation);

    /**
     * Return pointer to global singleton for desired interpolation type.
     */
    static ZC3CouponInterpolationConstSP make(const string& type);


    virtual ~ZC3CouponInterpolation();


    /**
     * Gets the interpolated coupon rates structure.
     *
     * This is done as follows:
     * 1. Generated the date list for the swap we are adding.
     * 2. For any coupon date which is greater than the last date in the curve
     *    but less than the maturity date of the swap we will interpolate the
     *    coupon rate.
     * 3. Interpolation will be between the coupon rate which applies for the
     *    last coupon date which is covered by the zero curve.
     *
     * Failures can occur as follows:
     * 1. If value floating is OFF then the price of this swap must be 1.
     *    Otherwise we do not know how to interpolate the price.
     *
     * Note that if no interpolated rates are required, then we return the
     * interpolated rates structure with zero rates.
     */
    virtual ZC3SwapDataArraySP getInterpRates(
        const ZC3SwapDataArray& swapsData,
        const ZeroCurve&        discountCurve,
        const ZeroCurve*        estimatingCurve,
        const BadDayConvention& badDayConv,
        const StubPlacement&    stubPos,
        const Holiday&          holidays,
        const Holiday&          basisHolidays,
        const bool              convDelayAdj,
        const IRVolBase*        volModelIR,
        const bool              withBasis) const = 0;

protected:
    ZC3CouponInterpolation(const CClassConstSP& clazz = TYPE);

private:
    static void load(CClassSP& clazz);
};


class MARKET_DLL ZC3LinearCouponInterpolation : public ZC3CouponInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3LinearCouponInterpolation(bool annualize);
    ~ZC3LinearCouponInterpolation();

    virtual ZC3SwapDataArraySP getInterpRates(
        const ZC3SwapDataArray& swapsData,
        const ZeroCurve&        discountCurve,
        const ZeroCurve*        estimatingCurve,
        const BadDayConvention& badDayConv,
        const StubPlacement&    stubPos,
        const Holiday&          holidays,
        const Holiday&          basisHolidays,
        const bool              convDelayAdj,
        const IRVolBase*        volModelIR,
        const bool              withBasis) const;

private:
    static void load(CClassSP& clazz);

    bool annualize;
};


/*
 * Zero interpolation.
 */

class MARKET_DLL ZC3ZeroInterpolation : public CObject
{
public:
    static bool lessThan(const ZC3ZeroInterpolationSP& first, const ZC3ZeroInterpolationSP& second);
    static CClassConstSP const TYPE;
    
    static const string LINEAR;
    static const string FLAT_FORWARDS;
    static const string SMOOTH_FORWARDS;
    static const string CUBIC_SPLINE;
    static const string EXP_CUBIC_SPLINE;

    /**
     * Return pointer to global singleton for desired interpolation type.
     */
    static ZC3ZeroInterpolationConstSP make(
        const string& type, 
        const DateTime& date = DateTime(),
        const bool      shaped = false);

    /**
     * Validate an interpolation format string.  This performs basic checks,
     * with full validation being done when the curve is bootstrapped (eg.
     * checks re. number of recursive interpolation types, or tenors).
     */
    static bool validate(const string& interpolation);

    /**
     * The format of the interpolation type specification string is 
     * 'interpType,duration,badDayConv;' for each section of the curve with a 
     * different interpolation type (obviously the duration and bad day 
     * convention would not be required for the final section).
     *
     * The bad day convention field is optional and defaults to the swapBadDay
     * input to ZERO_CURVE3.
     *
     * We have special processing for interpolation types which are used for
     * recursive builds. For these we check that there is only one
     * recursive section and that it is either the penultimate section (with
     * the last section being flat), or is the last section. 
     */
    static ZC3ZeroInterpolationArraySP createArray(
        const string& spec,
        const DateTime&         baseDate,
        const Holiday&          holidays,
        const BadDayConvention& defaultBadDayConv);


    virtual ~ZC3ZeroInterpolation();

    /**
     * Validates whether a combination of interpolation types is valid.
     */
    virtual void validateInterpType
        (const ZC3CouponInterpolation* couponInterp) const;

    /**
     * Get date until which this interpolation type is used.
     */
    const DateTime& until() const;

	/**
	 * Test if zero interpolation type can be used in fast interp function.
	 */
	bool hasFastInterp() const;

    /**
     * Test if zero interpolation type requires recursive bootstrap process.
     */
    // ALIB: zcbuild3.h#56
    virtual bool isRecursive() const;

    virtual bool isSmoothed() const;

    /**
     * Test if zero interpolation type uses the shape of a discount curve.
     * (cf. ALIB szcbuild.c#1038)
     */
    virtual bool isShaped() const;

    virtual double discountToRate(
        const ZC3ZeroCurve& curve, 
        double              discount, 
        const DateTime&     start, 
        const DateTime&     end) const;

    virtual double interpolate(
        const ZC3ZeroCurve&    zc3, 
        const DateTime&        x, 
        const ZC3CurveSegment& p1, 
        const ZC3CurveSegment& p2) const = 0;

    virtual double unsmoothedDiscountFactor(
        const ZC3ZeroCurve&    zc3,
        const ZC3CurveSegment& prevSegment,
        const ZC3CurveSegment& segment,
        const DateTime&        date,
        bool                   useSmoothing
        ) const = 0;

    virtual double smoothedDiscountFactor(
        const ZC3ZeroCurve&    zc3,
        const ZC3CurveSegment& segment,
        const DateTime&        startDate,
        double                 pvStart,
        const DateTime&        date
        ) const;

    virtual void solve(
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
        double            prevDiscount) const = 0;

    virtual DoubleArraySP spline(
        const DateTimeArray& benchmarkDates, 
        const DateTimeArray& fittingDates, 
        const DoubleArray&   benchmarkInputs) const = 0;

protected:
    ZC3ZeroInterpolation(const DateTime& date, bool fastInterp, const CClassConstSP& clazz = TYPE);

private:
    static void load(CClassSP& clazz);

    DateTime date;
	bool     fastInterp;
};


class MARKET_DLL ZC3LinearZeroInterpolation : public ZC3ZeroInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3LinearZeroInterpolation(const DateTime& date);
    ~ZC3LinearZeroInterpolation();

    virtual double interpolate(
        const ZC3ZeroCurve&    zc3,
        const DateTime&        x, 
        const ZC3CurveSegment& p1, 
        const ZC3CurveSegment& p2) const;

    virtual double unsmoothedDiscountFactor(
        const ZC3ZeroCurve&    zc3,
        const ZC3CurveSegment& prevSegment,
        const ZC3CurveSegment& segment,
        const DateTime&        date,
        bool                   useSmoothing) const;

    virtual void solve(
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
        double            prevDiscount) const;

    virtual DoubleArraySP spline(
        const DateTimeArray& benchmarkDates, 
        const DateTimeArray& fittingDates, 
        const DoubleArray&   benchmarkInputs) const;

private:
    static void load(CClassSP& clazz);
};


class MARKET_DLL ZC3FlatForwardsZeroInterpolation : public ZC3ZeroInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3FlatForwardsZeroInterpolation(
        const DateTime&      date,
		const bool           fastInterp = false,
        const bool           shaped = false, 
        const CClassConstSP& clazz = TYPE);
    ~ZC3FlatForwardsZeroInterpolation();

    virtual bool isShaped() const;
    void setShaped();

    virtual double interpolate(
        const ZC3ZeroCurve&    zc3,
        const DateTime&        x, 
        const ZC3CurveSegment& p1, 
        const ZC3CurveSegment& p2) const;

    virtual double unsmoothedDiscountFactor(
        const ZC3ZeroCurve&    zc3,
        const ZC3CurveSegment& prevSegment,
        const ZC3CurveSegment& segment,
        const DateTime&        date,
        bool                   useSmoothing) const;

    virtual void solve(
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
        double            prevDiscount) const;

    virtual DoubleArraySP spline(
        const DateTimeArray& benchmarkDates, 
        const DateTimeArray& fittingDates, 
        const DoubleArray&   benchmarkInputs) const;

private:
    static void load(CClassSP& clazz);
    bool shaped; // $unregistered
};


class MARKET_DLL ZC3SmoothForwardsZeroInterpolation : public ZC3FlatForwardsZeroInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3SmoothForwardsZeroInterpolation(const DateTime& date);
    ~ZC3SmoothForwardsZeroInterpolation();

    /**
     * Validates whether a combination of interpolation types is valid.
     */
    virtual void validateInterpType
        (const ZC3CouponInterpolation* couponInterp) const;

    bool isSmoothed() const;

    double smoothedDiscountFactor(
        const ZC3ZeroCurve&    szc,
        const ZC3CurveSegment& segment,
        const DateTime&        startDate,
        double                 startPv,
        const DateTime&        date
        ) const;

private:
    static void load(CClassSP& clazz);
};


class MARKET_DLL ZC3SplineZeroInterpolation : public ZC3FlatForwardsZeroInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3SplineZeroInterpolation(const DateTime& date, const CClassConstSP& clazz);
    ~ZC3SplineZeroInterpolation();

    virtual bool isRecursive() const;

protected:
	/**
	 * Performs one-dimensional cubic spline interpolation and extrapolation
	 * on an array of numbers.
	 *
	 * Returns the interpolated (or extrapolated) values.
	 */
	DoubleArraySP interpolate(
		const DateTimeArray& x,
		const DoubleArray&   f,
		const DateTimeArray& xDesired) const;

private:
    static void load(CClassSP& clazz);

	/**
	 * Performs one-dimensional cubic spline interpolation and extrapolation 
	 * for an array of points.  This routine is taken from "Numerical Recipes 
	 * in C", 1988, p.97. In the book, the routine name is "splint".
	 */
	double interpolatePoint(
		const DoubleArray&   secondDeriv,
		const DateTimeArray& x,
		const DoubleArray&   f,
		const DateTime&      xDesired) const;

	/**
	 * Computes second derivative of f at defined points. This must be called
	 * before the interpolation can take place. This routine is taken from
	 * "Numerical Recipes in C", 1988, pp. 96-97.
	 */
	DoubleArraySP secondDerivative(
		const DateTimeArray& x,
		const DoubleArray&   f) const;

	void binarySearch(
		const DateTime&      date, 
		const DateTimeArray& data, 
		int&                 lo, 
		int&                 hi) const;
};


class MARKET_DLL ZC3CubicSplineZeroInterpolation : public ZC3SplineZeroInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3CubicSplineZeroInterpolation(const DateTime& date);
    ~ZC3CubicSplineZeroInterpolation();

	/**
	 * This function does an interpolation of the discount factors.
	 */
    virtual DoubleArraySP spline(
        const DateTimeArray& benchmarkDates, 
        const DateTimeArray& fittingDates, 
        const DoubleArray&   benchmarkInputs) const;

private:
    static void load(CClassSP& clazz);
};


class MARKET_DLL ZC3ExpCubicSplineZeroInterpolation : public ZC3SplineZeroInterpolation
{
public:
    static CClassConstSP const TYPE;

    ZC3ExpCubicSplineZeroInterpolation(const DateTime& date);
    ~ZC3ExpCubicSplineZeroInterpolation();

	/**
	 * This function does an interpolation of the log of the discount factors.
	 */
    virtual DoubleArraySP spline(
        const DateTimeArray& benchmarkDates, 
        const DateTimeArray& fittingDates, 
        const DoubleArray&   benchmarkInputs) const;

private:
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE

#endif // ZC3_INTERPOLATION_HPP
