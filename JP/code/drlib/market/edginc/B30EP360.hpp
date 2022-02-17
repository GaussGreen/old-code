//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : B30EP360.hpp
//
//   Description : 30E+/360 Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef B30EP360_HPP
#define B30EP360_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** 30E+/360 Day count convention */

class MARKET_DLL B30EP360 : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class B30EP360Helper;

    B30EP360();
    virtual ~B30EP360();
    
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;
    
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
