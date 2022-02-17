//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Actual365FJ.hpp
//
//   Description : Act/365FJ Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef ACTUAL365FJ_HPP
#define ACTUAL365FJ_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/365FJ Day count convention */

class MARKET_DLL Actual365FJ : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class Actual365FJHelper;

    Actual365FJ();
    virtual ~Actual365FJ();
    
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
