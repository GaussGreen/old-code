//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : TimePoint2D.hpp
//
//   Description : Identifies a (DateTime, double, double) point
//
//   Author      : Antoine Gregoire
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef TIME_POINT_2D_HPP
#define TIME_POINT_2D_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Identifies a (DateTime, double, double) point */
class UTIL_DLL TimePoint2D : public CObject {
public:
    /** TYPE for TimePoint2D */
    static CClassConstSP const TYPE;

    /** Virtual destructor */
    virtual ~TimePoint2D();

    /** Constructor */
    TimePoint2D(const DateTime& date, double coord1, double coord2);
    
    /** Comparison by 1) date 2) coord1 3) coord2 */
    bool lessThanDateC1C2(const TimePoint2D& point) const;

    /** Comparison by 1) coord1 2) coord2 3) date */
    bool lessThanC1C2Date(const TimePoint2D& point) const;
    
    /** Access to date */
    DateTime getDate() const;
    
    /** Access to first coordinate */
    double getCoord1() const;

    /** Access to second coordinate */
    double getCoord2() const;
    
    /**
     * STL strict weak ordering functor for TimePoint2D,
     * defining the following order on TimePoint2D:
     * 1) increasing maturities
     * 2) increasing coord1
     * 3) increasing coord2 
     * */
    struct CompareDateC1C2 {
        /** _TimePoint2DCompareDateC1C2 binary function */
        UTIL_DLL bool operator()(const TimePoint2D& p1,
                                 const TimePoint2D& p2) const;
    };

private:
    /** Private constructor (for reflection) */
    TimePoint2D();
    
    /** Default constructor */
    static IObject* defaultTimePoint2D();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    // Fields
    DateTime date;
    double coord1;
    double coord2;
};

// Support for smart pointers
typedef smartPtr<TimePoint2D> TimePoint2DSP;

// Support for arrays
typedef array<TimePoint2DSP, TimePoint2D> TimePoint2DArray;
typedef smartPtr<TimePoint2DArray> TimePoint2DArraySP;

DRLIB_END_NAMESPACE

#endif /*TIME_POINT_2D_HPP*/

