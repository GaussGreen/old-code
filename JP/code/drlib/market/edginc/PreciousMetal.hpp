//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PreciousMetal.hpp
//
//   Description : Commodity asset for gold, silver etc.
//
//   Author      : Andrew J Swain
//
//   Date        : 30 September 2004
//
//
//----------------------------------------------------------------------------

#ifndef _PRECIOUSMETAL_HPP
#define _PRECIOUSMETAL_HPP

#include "edginc/ContangoCommodity.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL PreciousMetal: public ContangoCommodity {
public:
    static CClassConstSP const TYPE;

private:
    friend class PreciousMetalHelper;

    PreciousMetal();
    PreciousMetal(const PreciousMetal& rhs);
    PreciousMetal& operator=(const PreciousMetal& rhs);
};

DRLIB_END_NAMESPACE
#endif
