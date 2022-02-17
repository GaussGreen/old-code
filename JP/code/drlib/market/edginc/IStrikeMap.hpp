//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : IStrikeMap.hpp
//
//   Author      : Andrew McCleery
//
//   Description : Interface to map a strike into a different representation
//                 (e.g. forward/spot moneyness, delta, dollar div moneyness etc.)
//
//   Date        : October 25, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/Asset.hpp"

#ifndef ISTRIKEMAP_HPP
#define ISTRIKEMAP_HPP

DRLIB_BEGIN_NAMESPACE

/////////////////////////////////////////////////////////////////////////////////////
// STRIKE MAP INTERFACE /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

class MARKET_DLL IStrikeMap: public virtual IObject {
public:
    static  CClassConstSP const TYPE;

    // Map from a translated strike to an absolute strike
    virtual double toStrike(double pseudoStrike, const DateTime& maturity) = 0;

    // Map from an absolute strike to a translated strike
    virtual double fromStrike(double strike, const DateTime& maturity) = 0;
    
    virtual ~IStrikeMap();

protected:
    static void load(CClassSP& clazz);
    IStrikeMap();
};

typedef smartConstPtr<IStrikeMap> IStrikeMapConstSP;
typedef smartPtr<IStrikeMap> IStrikeMapSP;

/////////////////////////////////////////////////////////////////////////////////////
// FORWARD_MONEYNESS STRIKE MAP /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

class MARKET_DLL StrikeMapFwdMoneyness: public virtual IStrikeMap,
                                        public CObject {
public:
    static CClassConstSP const TYPE;

    // Map from forward moneynss to an absolute strike
    double toStrike(double fwdMoneyness, const DateTime& maturity);

    // Map from an absolute strike to forward moneyness
    double fromStrike(double strike, const DateTime& maturity);
        
    StrikeMapFwdMoneyness(const CAsset* curAsset);

    virtual ~StrikeMapFwdMoneyness() {};

protected:
    CAssetSP                    asset;          // Asset to use in the mapping

private:
    StrikeMapFwdMoneyness();
    static void load(CClassSP& clazz);
    static IObject* defaultStrikeMapFwdMoneyness();
};

/////////////////////////////////////////////////////////////////////////////////////
// DO NOTHING STRIKE MAP /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

class MARKET_DLL StrikeMapNone: public virtual IStrikeMap,
                                        public CObject {
public:
    static CClassConstSP const TYPE;

    // Map from absolute strike to absolute strike
    double toStrike(double strike, const DateTime& maturity);

    // Map from absolute strike to absolute strike
    double fromStrike(double strike, const DateTime& maturity);

    StrikeMapNone();
    //StrikeMapNone(const CAsset* curAsset);

    virtual ~StrikeMapNone() {};

private:

    static void load(CClassSP& clazz);
    static IObject* defaultStrikeMapNone();
};

DRLIB_END_NAMESPACE
#endif

