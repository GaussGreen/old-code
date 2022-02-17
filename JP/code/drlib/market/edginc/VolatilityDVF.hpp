//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolatilityDVF.hpp
//
//   Description : Interface for those vol's that return a CVolProcessedDVF
//
//   Author      : JNJ
//
//   Date        : 06 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLATILITY_DVF_HPP
#define VOLATILITY_DVF_HPP
#include "edginc/Object.hpp"
DRLIB_BEGIN_NAMESPACE

/** Contains methods implemented by vol which can provide a Deterministic Vol Function
    view of the world */
class MARKET_DLL IVolatilityDVF{
public:
    static CClassConstSP const TYPE;
    // no explicit methods currently
};

DRLIB_END_NAMESPACE
#endif
