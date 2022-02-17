//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : B30E360I.hpp
//
//   Description : 30E/360I Day count convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef B30E360I_HPP
#define B30E360I_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** 30E/360I Day count convention */

class MARKET_DLL B30E360I : public DayCountConvention {
public:
    static CClassConstSP const TYPE;
    friend class B30E360IHelper;

    B30E360I();
    virtual ~B30E360I();
    
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const;
    
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;
    
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
