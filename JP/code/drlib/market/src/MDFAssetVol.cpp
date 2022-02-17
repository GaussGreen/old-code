//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Helper class for models to get vols for 'CAssets'
//                 (ie 'spot' type assets) out of market cache
//
//   Author      : Mark A Robson
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/FlatFXVol.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/MarketDataConvert.hpp"

DRLIB_BEGIN_NAMESPACE


MDFAssetVol::MDFAssetVol(const string& volType): 
    volType(volType), fxVolType(CVolBase::TYPE->getName()), 
    allowVolTypeConversions(false), 
    volClass(0), fxVolClass(0), collapseToVolSurface(false) {
    initialise();
}

MDFAssetVol::MDFAssetVol(const string& volType,
                         const string& fxVolType): 
    volType(volType), fxVolType(fxVolType), 
    allowVolTypeConversions(false), 
    volClass(0), fxVolClass(0), collapseToVolSurface(false) {
    initialise();
}

void MDFAssetVol::setAllowVolTypeConversions(bool allowVolTypeConversions) {
    this->allowVolTypeConversions = allowVolTypeConversions;
}

void MDFAssetVol::initialise(){
    // turn off 'advanced' IRVol data - only needed for backward compatibility
    // (this class is/was used for swaptions priced without smile)
    setRetrievalMode(IRCalib::SmileBase::TYPE, false, NULL);
    setRetrievalMode(IRCalib::Model::TYPE, false, NULL);
    setRetrievalMode(IRCalib::TYPE, false, NULL);
}

/** key method - uses specifed type to refine search when looking for
    'equity' and fx vols. Does not default to VolSurface when VolPreferred 
    specified (cf MarketDataFetcherLN) unless collapseToVolSurface is set */
MarketObjectSP MDFAssetVol::fetch(
    const MarketData*    market,
    const string&        name,
    const CClassConstSP& type,
    const IModel*        model) const {

    static const string  method("MDFAssetVol::fetch");
    CClassConstSP        typeToUse = type; // default

    try {
        if (volClass == 0){
            volClass = CClass::forName(volType);
        }
        if (fxVolClass == 0){
            fxVolClass = CClass::forName(fxVolType);
        }
        // see if we're asking for a vol of some sort
        if (CVolBase::TYPE->isAssignableFrom(type)){ 
            if (FXVolBase::TYPE->isAssignableFrom(type)){
                // if it's not FlatFXVol look for an ordinary vol and 
                // then convert it
                if (!FlatFXVol::TYPE->isAssignableFrom(type)){
                    typeToUse = fxVolClass; // for error message
                    // get an ordinary vol
                    // call ourself to sort out vol preferred
                    // have a problem if the volType is set to a specific 
                    // parameterisation (rather than (say) surface) as it's
                    // unlikely we're going to have an FX vol of that type too.
                    // Sidestep by creating a new fetcher, looking for VolPreferred
                    MarketDataFetcherLN fxFetcher(fxVolType);
                    IObjectSP obj(fxFetcher.fetch(market, name, CVolBase::TYPE, model));
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

            // If collapseToVolSurface is set, see if we are looking for 
            // VolPreferred - as there's a real chance that it will not
            // be in the cache, we collapse to VolSurface if it's
            // not there rather than fail
            if (collapseToVolSurface &&
                CClass::forName("VolPreferred")->isAssignableFrom(volClass) &&
                !FXVolBase::TYPE->isAssignableFrom(type) &&
                !IRVolBase::TYPE->isAssignableFrom(type)) 
            {
                if (!market->hasData(name, volClass)) {
                    typeToUse = VolSurface::TYPE;
                } 
            }
        }

        MarketObjectSP obj;
        if(!allowVolTypeConversions) {
            // Previous behaviour
            obj = MarketDataFetcher::fetch(market, name, typeToUse, model);
        } else {
            // New behaviour with market data conversions
            obj = convertAndfetch(market, name, typeToUse, model);
        }
        
        return obj;
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed to get market data " + name +
                             " of type " + typeToUse->getName());
    }
}

/** Set the collapseToVolSurface attribute and return the previous value */
bool MDFAssetVol::setCollapseToVolSurface (bool collapse){
    bool oldValue = collapseToVolSurface;
    collapseToVolSurface = collapse;
    return oldValue;
}

/** Set the volType attribute and return the previous value */
string MDFAssetVol::setVolType(const string& newVolType){
    string oldVolType = volType;
    volType = newVolType;
    volClass = 0;
    return oldVolType;
}

MarketObjectSP MDFAssetVol::convertAndfetch(const MarketData*    market,
                                            const string&        name,
                                            const CClassConstSP& type,
                                            const IModel*        model) const {
    static const string method = "MDFAssetVol::convertAndfetch";
    
    try{
        // Support market data conversions
        if(market->hasData(name, type)) {
            // Do the usual stuff
            return MarketDataFetcher::fetch(market, name, type, model);
        } else {
            // Find conversion interface
            const MarketDataConvert::Reg& reg = MarketDataConvert::getConversionInfo(type);

            CClassConstSP iConvert = reg.iConvert;
            if(market->hasData(name, iConvert)) {
                // Find all QLib classes that implement the convert interface
                const CClassVec& types = iConvert->listAllConstructorClasses();;
                
                // Loop over all classes 
                MarketObjectSP objDesiredType;
                for(size_t iType = 0; iType < types.size(); iType++) {
                    if(market->hasData(name, types[iType])) {
                        // Fetch object that can be converted to desired type
                        MarketObjectSP objCandidate = MarketDataFetcher::fetch(market, name, types[iType], model);
                
                        // Invoke conversion
                        MarketObjectSP desiredCandidate = reg.method(objCandidate);
                        if(!!desiredCandidate) {
                            // Successful conversion
                            objDesiredType = desiredCandidate;
                            break;
                        }
                    }
                }

                if(!objDesiredType) {
                    throw ModelException(
                        "Failed for name" + name + 
                        ". No market object of type " + type->getName() + 
                        " or objects that can be converted to desired type.");
                }
                
                return objDesiredType;
            } else {
                throw ModelException("No market object of type " + type->getName() + 
                    " or objects that can be converted to desired type.");
            }
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


DRLIB_END_NAMESPACE
