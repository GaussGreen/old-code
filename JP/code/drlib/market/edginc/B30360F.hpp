//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : B30360F.hpp
//
//   Description : 30/360F Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 14 September 2001
//
//
//----------------------------------------------------------------------------

#ifndef B30360F_HPP
#define B30360F_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** 30/360F Day count convention */

class MARKET_DLL B30360F : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class B30360FHelper;

    B30360F();
    virtual ~B30360F();
    
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;
    
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
