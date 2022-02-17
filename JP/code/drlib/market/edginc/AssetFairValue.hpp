//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetFairValue.hpp
//
//   Description : What an asset must implement to calculate
//                 its fair value inside the AssetValue instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 7 September 2001
//
//
//----------------------------------------------------------------------------

#ifndef ASSETFAIRVALUE_HPP
#define ASSETFAIRVALUE_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for priceable assets */
class MARKET_DLL IAssetFairValue: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    /** return the fair value  */
    virtual double fairValue() const = 0;

protected:
    IAssetFairValue();
};


DRLIB_END_NAMESPACE
#endif
