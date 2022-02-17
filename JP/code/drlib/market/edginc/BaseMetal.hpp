//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BaseMetal.hpp
//
//   Description : Base metal asset
//
//   Author      : Andrew McCleery
//
//   Date        : 14 July 2005
//
//
//----------------------------------------------------------------------------

#ifndef _BaseMetal_HPP
#define _BaseMetal_HPP

#include "edginc/ContangoCommodity.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL BaseMetal: public ContangoCommodity {
public:
    static CClassConstSP const TYPE;

private:
    friend class BaseMetalHelper;

    BaseMetal();
    BaseMetal(const BaseMetal& rhs);
    BaseMetal& operator=(const BaseMetal& rhs);
};

DRLIB_END_NAMESPACE
#endif
