//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ActualActualAFB.hpp
//
//   Description : Act/Act AFB [Act/Act (Euro)] Day count convention
//
//   Author      : Andrew J Swain
//
//   Date        : 23 April 2004
//
//
//----------------------------------------------------------------------------

#ifndef ACTUALACTUALAFB_HPP
#define ACTUALACTUALAFB_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/Act AFB [Act/Act (Euro)] Day count convention */

class MARKET_DLL ActualActualAFB : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class ActualActualAFBHelper;

    ActualActualAFB();
    
    /** how many days between two dates */
    virtual int days(const DateTime &lodate, const DateTime &hidate) const;
    
    /** year fraction between two dates */
    virtual double years(const DateTime &lodate, const DateTime &hidate) const;
    
    /** returns a string description e.g. Act/Act (Euro) */
    virtual string toString() const;
};

DRLIB_END_NAMESPACE

#endif
