//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Inflation.hpp
//
//   Description : Inflation - acts like a ContangoCommodity asset
//
//   Author      : Andrew J Swain
//
//   Date        : 9 March 2005
//
//
//----------------------------------------------------------------------------

#ifndef _INFLATION_HPP
#define _INFLATION_HPP

#include "edginc/ContangoCommodity.hpp"

DRLIB_BEGIN_NAMESPACE

// Inflation - acts like a ContangoCommodity asset
class MARKET_DLL Inflation: public ContangoCommodity {
public:
    static CClassConstSP const TYPE;

private:
    friend class InflationHelper;

    Inflation();
    Inflation(const Inflation& rhs);
    Inflation& operator=(const Inflation& rhs);
};

DRLIB_END_NAMESPACE
#endif
