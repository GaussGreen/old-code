//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : B30360.hpp
//
//   Description : 30/360 Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 20 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef B30360_HPP
#define B30360_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** 30/360 Day count convention */

class MARKET_DLL B30360 : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class B30360Helper;

    B30360();
    virtual ~B30360();
    
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate,
                       bool  eomAdjSec,
                       bool  eomIgnoreLeapYear) const;
    
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;

    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate,
                         bool eomAdjSec) const;

    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
