//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Commodity.hpp
//
//   Description : Commodity asset base class 
//
//   Author      : Andrew J Swain
//
//   Date        : 2 July 2004
//
//
//----------------------------------------------------------------------------

#ifndef _COMMODITY_HPP
#define _COMMODITY_HPP
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

// more of a marker interface to tell commodities apart from everything else
class MARKET_DLL Commodity: public CAsset {
public:
    static CClassConstSP const TYPE;

protected:
    Commodity(CClassConstSP clazz); 
private:
    friend class CommodityHelper;
    Commodity();
};

typedef smartPtr<Commodity> CommoditySP;
typedef smartConstPtr<Commodity> CommodityConstSP;

DRLIB_END_NAMESPACE
#endif
