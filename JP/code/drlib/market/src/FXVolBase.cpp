//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVolBase.cpp
//
//   Description : Abstract base for FX vols
//
//   Author      : Mark A Robson
//
//   Date        : 5 Dec 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_FXVOLBASE_CPP
#include "edginc/FXVolBase.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/FlatFXVol.hpp"

DRLIB_BEGIN_NAMESPACE


FXVolBase::~FXVolBase() {}

FXVolBase::FXVolBase(const CClassConstSP& clazz): 
        CVolBase(clazz) {}


// centralised method to handle both legacy and current way of using
// FX vols for ccy protected/struck assets. Given model, market, the name
// and type of the FX vol you "think" you want, it gives you back the type 
// you should be using.
// the rule is:
// 1. look for a FlatFXVol of the given name
// 2. if it exists, return the class of FlatFXVol
// 3. otherwise return the input clazz
CClassConstSP FXVolBase::volClassProtStruck(const IModel*     model,
                                            const MarketData* market,
                                            const string&     name,
                                            CClassConstSP     clazz) {
    static const string method("FXVolBase::volClassProtStruck");
    try {
        if (market->hasData(name, FlatFXVol::TYPE)) {
            return FlatFXVol::TYPE;
        }
        return clazz;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


void FXVolBase::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FXVolBase, clazz);
    SUPERCLASS(CVolBase);
}

CClassConstSP const FXVolBase::TYPE = CClass::registerClassLoadMethod(
    "FXVolBase", typeid(FXVolBase), load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(FXVolBaseWrapper);

DRLIB_END_NAMESPACE
