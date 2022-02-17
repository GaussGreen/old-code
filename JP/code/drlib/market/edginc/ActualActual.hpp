//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ActualActual.hpp
//
//   Description : Act/Act Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 20 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef ACTUALACTUAL_HPP
#define ACTUALACTUAL_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Act/Act Day count convention - same as Act/365 */

class MARKET_DLL ActualActual : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class ActualActualHelper;

    ActualActual();
    virtual ~ActualActual();
    
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
