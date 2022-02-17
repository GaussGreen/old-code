//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolatilityBS.hpp
//
//   Description : Interface for those vol's that return a CVolProcessedBS
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLATILITY_BS_HPP
#define VOLATILITY_BS_HPP
#include "edginc/Object.hpp"
DRLIB_BEGIN_NAMESPACE

/** Contains methods implemented by vol which can provide a Black-Scholes
    view of the world */
class MARKET_DLL IVolatilityBS{
public:
    static CClassConstSP const TYPE;
    // no explicit methods currently
};

DRLIB_END_NAMESPACE
#endif
