//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DayCountConventionFactory.hpp
//
//   Description : Factory class for building DayCountConventions
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef DAYCOUNTCONVENTIONFACTORY_HPP
#define DAYCOUNTCONVENTIONFACTORY_HPP

#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Factory class for building DayCountConventions */
class MARKET_DLL DayCountConventionFactory {
public:
    static DayCountConvention* make(const string& dcc);
    static DayCountConvention* clone(const DayCountConvention* dcc);
    
private:
    DayCountConventionFactory();  
};

DRLIB_END_NAMESPACE

#endif
