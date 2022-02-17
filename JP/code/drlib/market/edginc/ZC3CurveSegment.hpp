//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3CurveSegment.hpp
//
//   Description : Helper class for ALIB zero curve 3 bootstrapping method.
//
//   Author      : Richard Appleton
//
//   Date        : 13th May 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_CURVE_SEGMENT_HPP
#define ZC3_CURVE_SEGMENT_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ZC3Interpolation.hpp"



DRLIB_BEGIN_NAMESPACE
class DayCountConvention;
class ZC3ZeroCurve;
class IRVolBase;
class CashStream;
class ZeroCurve;


/*
 * Infrastructure support. 
 */
class ZC3CurveSegment;
typedef smartPtr<ZC3CurveSegment> ZC3CurveSegmentSP;
typedef array<ZC3CurveSegment>    ZC3CurveSegmentArray;


class MARKET_DLL ZC3CurveSegment : public CObject
{
public:
    static double discountFactor(
        const ZC3CurveSegmentArray& segments, 
        const ZC3ZeroCurve&         szc,
        const DateTime&             date, 
        bool                        useSmoothing);
    
    static double unadjustedZeroCouponRate(
        const ZC3CurveSegmentArray& segments, 
        const ZC3ZeroCurve&         szc,
        const DateTime&             date,
        int&                        loBound,
        int&                        hiBound);

    static void   calculateAverageRates(
        ZC3CurveSegmentArray& segments, 
        const ZC3ZeroCurve&   szc,
        int                   numStubs);

	static bool hasFastInterp(const ZC3CurveSegmentArray& segments);

    static CClassConstSP const TYPE;

    ZC3CurveSegment();

    ZC3CurveSegment(
        const DateTime&             date,
        const ZC3ZeroInterpolation& interpolation,
        double                      rate,
        double                      discount,
        double                      smoothedDiscount,
        int                         index);
 
    ZC3CurveSegment(
        const DateTime&             date,
        const ZC3ZeroInterpolation& interpolation,
        const ZC3CurveSegment&      prevSegment);

	ZC3CurveSegment(
		const DateTime&             discDate,
		const ZC3CurveSegment&      prevSegment,
		double                      discount,
		const ZC3ZeroCurve&         curve);

    /**
     * true if segment has been smoothed.
     */
    bool isSmoothed() const;

    /**
     * Get discount factor at date.
     */
    double getDiscount() const;

    /**
     * Get discount factor at date.
     */
    double getDiscount(bool smoothed) const;

    /**
     * Get continuous compounded rate.  Used for flat forwards fast interpolation.
     */
    double getContinuousRate() const;

    /**
     * Calculate the unsmoothed discount
     */
    double unsmoothedDiscountFactor(
        const ZC3ZeroCurve&    szc,
        const ZC3CurveSegment& prevSegment,
        const DateTime&        date,
        bool                   useSmoothing) const;

    /**
     * Calculate the smooth discount
     */
    double smoothedDiscountFactor(
        const ZC3ZeroCurve&    szc,
        const ZC3CurveSegment& prevSegment,
        const DateTime&        date) const;

    /**
     * Calculates the unsmoothed segment average rate (continuously compounded).
     */
    void averageRate(const ZC3ZeroCurve& szc, bool forceRecalc);

    /**
     * Calculates the duration of a segment of a smooth zero curve.
     *
     * The method chosen is somewhat artificial and agrees with the existing
     * implementation.
     */
    void segmentDuration(const ZC3ZeroCurve& szc);

    /**
     * Complete the smoothing for the segment.
     * This is where we calculate the parabola.
     */
    void smooth(
        ZC3ZeroCurve&    szc, 
        double           startFwdRate,
        const ZeroCurve* discountCurve,
        const IRVolBase* volModelIR);

    /*
     * Smooths segment.  Always recalculates rate, discount and avgRate as the 
     * previous segment may have been smoothed.
     */
    void recalculate(
        ZC3ZeroCurve&     szc, 
        ZC3CurveSegment&  prevSegment,
        const CashStream& cs,
        double            startFwdRate,
        double            presentValue,
        double            pvKnown,
        double            pvRemainder,
        bool              foundAdjAbs,
        const ZeroCurve*  discountCurve,
        const IRVolBase*  volModelIR);

    /**
     * Calculates discount from parabola coefficients.
     */
    double parabDiscount(double prevDiscount, double time) const;

    
private:
    friend class ZC3ZeroCurve;
    friend class ZC3LinearZeroInterpolation;
    friend class ZC3FlatForwardsZeroInterpolation;
    friend class ZC3SmoothForwardsZeroInterpolation;
    friend class LinearObjFunc;
    friend class LinearObjFunc1;
    friend class LinearObjFunc2;
    friend class FlatFwdObjFunc;
    friend class SmoothFunc;
    friend class SmoothFirstObjFunc;
    friend class SmoothLastObjFunc;
    friend class SmoothParabObjFunc;
	friend class ZC3RecurseInfo;
    friend class LMParabCurve;
    
    /** search for segments surrounding desired date */
    static void   binarySearch(
        const ZC3CurveSegmentArray& segments, 
        int&                        lo, 
        int&                        hi, 
        const int                   date);

    static void load(CClassSP& clazz);

    DateTime                    date;           // date at end of this segment
    ZC3ZeroInterpolationConstSP interpolation;  // simple interpolation type for this segment
    double                      discount;       // unsmoothed discount factor at date
    double                      rate;           // unsmoothed rate in basis and dcc
    /* Note that rate and discount are equivalent */
    double                      smoothDiscount; // smoothed discount factor at date
    double                      avgRate;        // unsmoothed continuously compounded rate for this segment
    double                      a0;             // smoothing parameter
    double                      a1;             // smoothing parameter
    double                      a2;             // smoothing parameter
    double                      duration;       // smoothing parameter
    int                         index;          // position within data array

    // transient field for flat forwards fast interpolation
    mutable double              contRate;       // continuous rate
};



/**
 * Specialization of arrayObjectCast (needed as arrays are not array of pointers)
 */
template <> class MARKET_DLL arrayObjectCast<ZC3CurveSegment>
{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ZC3CurveSegment& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ZC3CurveSegment& value);

    /** Turns the IObjectSP into a ZC3CurveSegment */
    static ZC3CurveSegment fromIObject(IObjectSP& value);
};


DRLIB_END_NAMESPACE
#endif // ZC3_CURVE_SEGMENT_HPP
