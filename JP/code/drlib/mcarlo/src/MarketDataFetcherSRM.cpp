//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcherSRM.cpp
//
//   Description : Helper class for models to get data out of market cache
//
//   Author      : Mark A Robson
//
//   Date        : June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MarketDataFetcherSRM.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/Format.hpp"
#include "edginc/CRCalib.hpp"
#include "edginc/MCPathConfigSRM.hpp"

DRLIB_BEGIN_NAMESPACE
MarketDataFetcherSRM::~MarketDataFetcherSRM(){}

MarketDataFetcherSRM::MarketDataFetcherSRM(
    const string&        irCalibSmileType, 
    const string&        irCalibModelType,
    bool                 getSwaptionVols,
	bool				 useIRVolPair,
    const string&        volType,
    const StringArray&   fxVolType, 
    const string&        crCalibSmileType,
    const string&        cdsVolType):
    fxVolType(fxVolType)
{
    // get ir vol inside yield curves
    setRetrievalMode(IRVolBase::TYPE, IYieldCurve::TYPE, true, NULL);
    // get cds vol in ICDSParSpreads
    setRetrievalMode(ICDSVol::TYPE, ICDSParSpreads::TYPE, true, NULL);
    setSwaptionVolFlag(getSwaptionVols);
    setIrCalibSmileType(irCalibSmileType);
    setIrCalibModelType(irCalibModelType);
    // do crCalibSmileClass, if specified - currently no set method
    if (!crCalibSmileType.empty()){
        CClassConstSP crCalibSmileClass;
        try {
            crCalibSmileClass = CClass::forName(crCalibSmileType);
        } catch (exception& e){
            throw ModelException(e, "Name for CR Calib smile type not "
                                 "recognised");
        }
        setRetrievalMode(CRCalib::Smile::TYPE, true, crCalibSmileClass);
    }
    // repeat for equity vol type, if specified
    if (!volType.empty()){
        CClassConstSP eqVolClass;
        try {
            eqVolClass = CClass::forName(volType);
        } catch (exception& e){
            throw ModelException(e, "Name for Equity Vol type not recognised");
        }
        setRetrievalMode(CVolBase::TYPE, true, eqVolClass);
    }
    // repeat for cds vol type, if specified
    if (!cdsVolType.empty()){
        CClassConstSP cdsVolClass;
        try {
            cdsVolClass = CClass::forName(cdsVolType);
        } catch (exception& e){
            throw ModelException(e, "Name for CDS Vol type not recognised");
        }
        setRetrievalMode(ICDSVol::TYPE, ICDSParSpreads::TYPE, true, cdsVolClass);
    }
}

/** Pulls out irCalibSmileType and irCalibModelType objects from cache */
MarketObjectSP MarketDataFetcherSRM::fetch(
    const MarketData*    market,
    const string&        name,
    const CClassConstSP& type,
    const IModel*        model) const {
    static const string method("MarketDataFetcherSRM::fetch");
    try {
        CClassConstSP clazzToUse = type; // default

#if 0
//FIXME: don't know what we do here:
                        //clazzToUse = getSwaptionVols? IRVol::TYPE : IRCalib::TYPE;
            if ( getSwaptionVols ) {
                clazzToUse = useIRVolPair ? IRVolPair::TYPE : IRVol::TYPE;
            } else {
	            clazzToUse = IRCalib::TYPE;
            }
#endif
        if (FXVolBase::TYPE->isAssignableFrom(type)){
            // XXX to be revisited fxVolType should really be of type map<string, string>
            // where 'first' is the market name (eg EURUSD) and 'second' is the 
            // fx vol type for that particular market name.
            // for now simply validate that fxVolType is of size 1
            if (fxVolType.size() != 1){
                throw ModelException(method,
                                     "fxVolType should be of size 1; got "
                                     + Format::toString(fxVolType.size()));
            }
            clazzToUse = CClass::forName(fxVolType.front());
        }
        return MarketDataFetcher::fetch(market, name, clazzToUse, model);
    } catch (exception& e) {
        throw ModelException(e, 
                             method, 
                             "Failed to get market data " + name + " of type " +
                             type->getName());
    }
}

/** change the irCalibSmileType */
void MarketDataFetcherSRM::setIrCalibSmileType(const string& irCalibSmileType){
    if (!irCalibSmileType.empty()){
        CClassConstSP irCalibSmileClass;
        try {
            irCalibSmileClass = CClass::forName(irCalibSmileType);
        } catch (exception& e){
            throw ModelException(e, "Name for IR Calib smile "
                                 "type not recognised");
        }
        setRetrievalMode(IRCalib::SmileBase::TYPE, true, irCalibSmileClass);
    }
}

/** change the irCalibModelType */
void MarketDataFetcherSRM::setIrCalibModelType(const string& irCalibModelType){
    if (!irCalibModelType.empty()){
        CClassConstSP irCalibModelClass;
        try {
            irCalibModelClass = CClass::forName(irCalibModelType);
        } catch (exception& e){
            throw ModelException(e, "Name for IR Calib model "
                                 "type not recognised");
        }
        setRetrievalMode(IRCalib::Model::TYPE, true, irCalibModelClass);
    }
}

/** Setting this to true, means that a request for a IRVolBase is turned
    into a request for an IRVol otherwise it is turned into one for 
    an IRCalib */
void MarketDataFetcherSRM::setSwaptionVolFlag(bool getSwaptionVols){
    // get either IRVol or IRCalib depending on flag
    setRetrievalMode(IRVolBase::TYPE, IYieldCurve::TYPE, true, 
                     getSwaptionVols? IRVolCommon::TYPE : IRCalib::TYPE);
    // some products still have explicit IRVolBases in them ie outside of YCs
    setRetrievalMode(IRVolBase::TYPE, true, 
                     getSwaptionVols? IRVolCommon::TYPE : IRCalib::TYPE);
}

DRLIB_END_NAMESPACE
