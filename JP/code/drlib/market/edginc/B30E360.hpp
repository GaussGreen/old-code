//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : B30E360.hpp
//
//   Description : 30E/360 Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 20 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef B30E360_HPP
#define B30E360_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** 30E/360 Day count convention */

class MARKET_DLL B30E360 : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class B30E360Helper;

    B30E360();
    virtual ~B30E360();
    
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;
    
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
