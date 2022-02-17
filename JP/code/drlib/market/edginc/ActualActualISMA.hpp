//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ActualActualISMA.hpp
//
//   Description : Act/Act ISMA Day count convention
//
//   Author      : Andrew J Swain
//
//   Date        : 23 April 2004
//
//
//----------------------------------------------------------------------------

#ifndef ACTUALACTUALISMA_HPP
#define ACTUALACTUALISMA_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/Act ISMA Day count convention */

class MARKET_DLL ActualActualISMA : public DayCountConvention {
public:
    static CClassConstSP const TYPE;

    /** how many days between two dates */
    virtual int days(const DateTime &lodate, const DateTime &hidate) const;
    
    // returns year fraction for accrued interest calculations.
    // regular coupon period
    double accruedFactor(const DateTime &aiDate,
                             const DateTime &prevDate,
                             const DateTime &nextDate,
                             int couponFreq,
                             bool goneEx) const;

    // returns year fraction for accrued interest calculations.
    // regular and irregular coupon periods
    double accruedFactor(const DateTime &aiDate, 
                                 const DateTime &prevDate, 
                                 const DateTime &nextDate, 
                                 int couponFreq,
                                 const string &couponType,
                                 bool eomAdjSec,
                                 bool goneEx) const;

    // derived classes must implement this to call the regular version 
    // above with the coupon frequency
    virtual double accruedFactor(const DateTime &aiDate,
                                 const DateTime &prevDate,
                                 const DateTime &nextDate,
                                 bool  eomAdjSec,
                                 bool  goneEx) const = 0;

    /** year fraction between two dates */
    /** only works for "regular" coupons */
    double years(const DateTime &lodate, 
                 const DateTime &hidate) const;
   
    /** returns a string description e.g. Act/Act (ISMA) */
    virtual string toString() const;

    /** create subclasses */
    static ActualActualISMA* make(int couponsPerYear);

protected:
    ActualActualISMA(CClassConstSP clazz);
    ActualActualISMA();
    friend class ActualActualISMAHelper;
};

DRLIB_END_NAMESPACE


#endif
