//----------------------------------------------------------------------------
//
//   Filename    : FixingType.hpp
//
//   Description : Extra context for retrieving samples for ObservableHistory
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//                 We include the (trivial) Asset specific type AssetFixType
//
//   Author      : Ian Stares
//
//   Date        : February 2 2006
//
//
//----------------------------------------------------------------------------

#ifndef FIXINGTYPE_HPP
#define FIXINGTYPE_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL FixingType :  public CObject {
public:
    static CClassConstSP const TYPE;

    virtual bool equals(const FixingType& other) const;

protected:
    FixingType(CClassConstSP clazz);

private:
    FixingType();

    static void load(CClassSP& clazz);
    static IObject* defaultFixingType();
};
typedef smartPtr<FixingType> FixingTypeSP;

// extra (!) context for asset sample look ups
class MARKET_DLL AssetFixType :  public FixingType {
public:
    static CClassConstSP const TYPE;

    AssetFixType(string assetName);

    virtual bool equals(const FixingType& other) const;

private:
    AssetFixType();

    static void load(CClassSP& clazz);
    static IObject* defaultAssetFixType();

    string assetName;
};
typedef smartPtr<AssetFixType> AssetFixTypeSP;

DRLIB_END_NAMESPACE

#endif




