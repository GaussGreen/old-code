//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CommodityIndex.hpp
//
//   Description : Commodity index asset
//
//   Author      : Andrew J Swain
//
//   Date        : 4 October 2004
//
//
//----------------------------------------------------------------------------

#ifndef _CommodityIndex_HPP
#define _CommodityIndex_HPP

#include "edginc/ContangoCommodity.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CommodityIndex: public ContangoCommodity {
public:
    static CClassConstSP const TYPE;

private:
    friend class CommodityIndexHelper;

    CommodityIndex();
    CommodityIndex(const CommodityIndex& rhs);
    CommodityIndex& operator=(const CommodityIndex& rhs);
};

DRLIB_END_NAMESPACE
#endif
