//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Actual365F.hpp
//
//   Description : Act/365F Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef ACTUAL365F_HPP
#define ACTUAL365F_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/365F Day count convention */

class MARKET_DLL Actual365F : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class Actual365FHelper;

    Actual365F();
    virtual ~Actual365F();
    
    /** how many days between two dates */
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    /** year fraction between two dates */
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;
    
    /** returns a string description e.g. Act/360 */
    virtual string toString() const;

    /** year fraction between two dates */
    static double yearFraction(const DateTime &lodate, 
                               const DateTime &hidate);
};

DRLIB_END_NAMESPACE

#endif
