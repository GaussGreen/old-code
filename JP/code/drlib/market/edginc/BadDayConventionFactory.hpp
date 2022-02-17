//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadDayConventionFactory.hpp
//
//   Description : Factory class for building BadDayConventions
//
//   Author      : Andrew J Swain
//
//   Date        : 15 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef BADDAYCONVENTIONFACTORY_HPP
#define BADDAYCONVENTIONFACTORY_HPP

#include "edginc/BadDayConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Factory class for building BadDayConventions */
class MARKET_DLL BadDayConventionFactory {
public:
    static BadDayConvention* make(const string& convention);
    static BadDayConvention* clone(const BadDayConvention* convention);
    
private:
    BadDayConventionFactory();  
};

DRLIB_END_NAMESPACE

#endif
