//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : IStrikeMap.cpp
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
#include "edginc/IStrikeMap.hpp"

#ifndef ISTRIKEMAP_CPP
#define ISTRIKEMAP_CPP

DRLIB_BEGIN_NAMESPACE

IStrikeMap::IStrikeMap()  {}

IStrikeMap::~IStrikeMap() 
{}
 
void IStrikeMap::load(CClassSP& clazz){
    REGISTER_INTERFACE(IStrikeMap, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IStrikeMap::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IStrikeMap", typeid(IStrikeMap), load);

/////////////////////////////////////////////////////////////////////////////////////
// FORWARD_MONEYNESS STRIKE MAP /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// Map from forward moneynss to an absolute strike
double StrikeMapFwdMoneyness::toStrike(double fwdMoneyness, const DateTime& maturity) {
    static const string method = "StrikeMapFwdMoneyness::toStrike";
    try {
        if (!asset) {
            throw ModelException(method, "Internal Error - Asset has not been initialised");
        }
        double strike = fwdMoneyness * asset->fwdValue(maturity);
        return strike;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Map from an absolute strike to forward moneyness
double StrikeMapFwdMoneyness::fromStrike(double strike, const DateTime& maturity) {
    static const string method = "StrikeMapFwdMoneyness::fromStrike";
    try {
        if (!asset) {
            throw ModelException(method, "Internal Error - Asset has not been initialised");
        }
        double fwdMoneyness = strike / asset->fwdValue(maturity);
        return fwdMoneyness;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

StrikeMapFwdMoneyness::StrikeMapFwdMoneyness(const CAsset* curAsset) : CObject(TYPE) {
    asset = AssetSP(dynamic_cast<Asset *>(curAsset->clone()));
}

StrikeMapFwdMoneyness::StrikeMapFwdMoneyness() : CObject(TYPE) {}

void StrikeMapFwdMoneyness::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(StrikeMapFwdMoneyness, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IStrikeMap);
    EMPTY_SHELL_METHOD(defaultStrikeMapFwdMoneyness);
    FIELD_NO_DESC(asset);
    FIELD_MAKE_TRANSIENT(asset);
}

IObject* StrikeMapFwdMoneyness::defaultStrikeMapFwdMoneyness(){
        return new StrikeMapFwdMoneyness();
}

CClassConstSP const StrikeMapFwdMoneyness::TYPE = CClass::registerClassLoadMethod(
    "StrikeMapFwdMoneyness", typeid(StrikeMapFwdMoneyness), StrikeMapFwdMoneyness::load);

/////////////////////////////////////////////////////////////////////////////////////
// DO NOTHING STRIKE MAP /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// Map from absolute strike to absolute strike
double StrikeMapNone::toStrike(double strike, const DateTime& maturity) {
    static const string method = "StrikeMapNone::toStrike";
    try {
        return strike;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Map from absolute strike to absolute strike
double StrikeMapNone::fromStrike(double strike, const DateTime& maturity) {
    static const string method = "StrikeMapNone::fromStrike";
    try {
        return strike;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//StrikeMapNone::StrikeMapNone() : CObject(TYPE) {
    // No asset needed
//}

StrikeMapNone::StrikeMapNone() : CObject(TYPE) {};

void StrikeMapNone::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(StrikeMapNone, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IStrikeMap);
    EMPTY_SHELL_METHOD(defaultStrikeMapNone);
}

IObject* StrikeMapNone::defaultStrikeMapNone(){
        return new StrikeMapNone();
}

CClassConstSP const StrikeMapNone::TYPE = CClass::registerClassLoadMethod(
    "StrikeMapNone", typeid(StrikeMapNone), StrikeMapNone::load);

DRLIB_END_NAMESPACE
#endif

