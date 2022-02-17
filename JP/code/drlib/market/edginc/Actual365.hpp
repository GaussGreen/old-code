//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Actual365.hpp
//
//   Description : Act/365 Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 20 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef ACTUAL365_HPP
#define ACTUAL365_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/365 Day count convention */

class MARKET_DLL Actual365 : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class Actual365Helper;

    Actual365();
    virtual ~Actual365();
    
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
