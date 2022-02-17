//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcherLN.cpp
//
//   Description : Helper class for LN models to get data out of market cache
//
//   Author      : Andrew J Swain
//
//   Date        : 1 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MARKETDATAFETCHERLN_CPP
#include "edginc/IRCalib.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/FlatFXVol.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/CurrencyBasis.hpp"


DRLIB_BEGIN_NAMESPACE

MarketDataFetcherLN::MarketDataFetcherLN(const string& volType, 
                                         const bool useCcyBasis): 
    volType(volType), volClass(0) {
    setRetrievalMode(CurrencyBasis::TYPE, useCcyBasis, NULL);
    initialise();
}

MarketDataFetcherLN::MarketDataFetcherLN(): 
    volType(IVolatilityBS::TYPE->getName()), volClass(0) {
    initialise();
}

//// common part of constructors
void MarketDataFetcherLN::initialise(){
    // turn off 'advanced' IRVol data - only needed for backward compatibility
    // (this class is/was used for swaptions priced without smile)
    setRetrievalMode(IRCalib::SmileBase::TYPE, false, NULL);
    setRetrievalMode(IRCalib::Model::TYPE, false, NULL);
    setRetrievalMode(IRCalib::TYPE, false, NULL);
    // look up class
    volClass = CClass::forName(volType);
    // check to see it's a BS vol
    if (!IVolatilityBS::TYPE->isAssignableFrom(volClass)){
        throw ModelException("MarketDataFetcherLN::intialise",
                             "Specified volatility type (" +
                             volType + ") is not a Log Normal vol");
    }
}

/** LN fetch method */
MarketObjectSP MarketDataFetcherLN::fetch(
    const MarketData*    market,
    const string&        name,
    const CClassConstSP& type,
    const IModel*        model) const {

    static const string  method("MarketDataFetcherLN::fetch");
    static CClassConstSP volPreferred = CClass::forName("VolPreferred");
    CClassConstSP        typeToUse = type; // default
    try {
        if (CVolBase::TYPE->isAssignableFrom(type)){ 
            if (FXVolBase::TYPE->isAssignableFrom(type)){
                // if it's not FlatFXVol look for an ordinary vol and
                // then convert it
                if (!FlatFXVol::TYPE->isAssignableFrom(type)){
                    typeToUse = CVolBase::TYPE; // for error message
                    // get an ordinary vol call ourself to sort out
                    // vol preferred have a problem if the volType is
                    // set to a specific parameterisation (rather than
                    // (say) surface) as it's unlikely we're going to
                    // have an FX vol of that type too.  Sidestep by
                    // creating a new fetcher, looking for
                    // VolPreferred
                    MarketDataFetcherLN fxFetcher("VolPreferred");
                    IObjectSP obj(fxFetcher.fetch(market, name, CVolBase::TYPE, 
                        model));
                    // then turn it into an fx vol
                    CObject::checkType(obj, type);
                    return MarketObjectSP::dynamicCast(obj);
                }
            } else if (volClass->isAssignableFrom(type)){
                /* if type requested is already of the type we want
                   (or is derived from the type we want) then let it through */
            } else {
                // use our specific type
                typeToUse = volClass;
            }

            // see if we looking for VolPreferred - as there's a real chance
            // that won't be in the cache, we collapse to VolSurface if it's
            // not there rather than fail
            if (volPreferred->isAssignableFrom(volClass) &&
                !FXVolBase::TYPE->isAssignableFrom(type) &&
                !IRVolBase::TYPE->isAssignableFrom(type)) {
                if (!market->hasData(name, volClass)) {
                    typeToUse = VolSurface::TYPE;
                } 
            }
        }
        // then just invoke parent with overridden type
        return MarketDataFetcher::fetch(market, name, typeToUse, model);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed to get market data " + name +
                             " of type " + typeToUse->getName());
    }
}


DRLIB_END_NAMESPACE
