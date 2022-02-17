//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVolBase.hpp
//
//   Description : Abstract base for FX vols
//
//   Author      : Mark A Robson
//
//   Date        : 5 Dec 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FX_VOLBASE_HPP
#define EDR_FX_VOLBASE_HPP

#include "edginc/VolBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** FXVols are like other vols, just that the greeks come out
    differently.  (parallel vega comes out as fx vega etc). This class
    just allows classes to insist upon an FXVol if they should so
    desire. The main derived class from this class is FXVol which
    wraps any CVolBase and takes care of the greeks. */
class MARKET_DLL FXVolBase: public CVolBase {
public:
    static CClassConstSP const TYPE;

    virtual ~FXVolBase();

    // centralised method to handle both legacy and current way of using
    // FX vols for ccy protected/struck assets. Given model, market, the name
    // and type of the FX vol you "think" you want, it gives you back the type 
    // you should be using.
    static CClassConstSP volClassProtStruck(const IModel*     model,
                                            const MarketData* market,
                                            const string&     name,
                                            CClassConstSP     clazz);
                                            
protected:
    FXVolBase(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);

    FXVolBase(const CVolBase &rhs);
    FXVolBase& operator=(const CVolBase& rhs);
};

typedef smartConstPtr<FXVolBase> FXVolBaseConstSP;
typedef smartPtr<FXVolBase> FXVolBaseSP;
#ifndef QLIB_FXVOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<FXVolBase>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<FXVolBase>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<FXVolBase>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<FXVolBase>);
#endif

// support for wrapper class
typedef MarketWrapper<FXVolBase> FXVolBaseWrapper;
#ifndef QLIB_FXVOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<FXVolBase>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<FXVolBase>);
#endif

DRLIB_END_NAMESPACE

#endif

