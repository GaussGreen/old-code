//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Actual360.hpp
//
//   Description : Act/360 Day count convention
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef ACTUAL360_HPP
#define ACTUAL360_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/360 Day count convention */

class MARKET_DLL Actual360 : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class Actual360Helper;

    Actual360();
    virtual ~Actual360();
    
    /** how many days between two dates */
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    /** year fraction between two dates */
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;
    
    /** returns a string description e.g. Act/360 */
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
