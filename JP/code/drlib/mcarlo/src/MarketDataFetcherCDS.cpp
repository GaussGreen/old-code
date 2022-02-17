//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcherCDS.cpp
//
//   Description : Helper class for models that use CDS Par Spreads.
//                 It's in the mcarlo directory since we might need this for
//                 a credit mc
//
//   Author      : Mark A Robson
//
//   Date        : 29 Nov 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MarketDataFetcherCDS.hpp"
#include "edginc/QuantoCDSParSpreads.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/ICDSSpotVol.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/SRMFXVol.hpp"

DRLIB_BEGIN_NAMESPACE
MarketDataFetcherCDS::~MarketDataFetcherCDS(){}

//// default constructor
MarketDataFetcherCDS::MarketDataFetcherCDS(){
    initialise(false); // do default settings
}

/** Will fetch stochastic yield curves plus type of smile and model
    specified. Will probably need some sort of vol choice once we know
    what we're doing */
MarketDataFetcherCDS::MarketDataFetcherCDS(bool getSwaptionVols){
    initialise(getSwaptionVols); // do default settings
    // get either IRVol or IRCalib depending on flag
    setRetrievalMode(IRVolBase::TYPE, true, 
        getSwaptionVols? IRVolCommon::TYPE : IRCalib::TYPE);
}

/** Will fetch stochastic yield curves plus type of smile and model
    specified. Will probably need some sort of vol choice once we know
    what we're doing */
MarketDataFetcherCDS::MarketDataFetcherCDS(bool getSwaptionVols, 
                                           bool getCreditVols,
					   CClassConstSP creditVolType) {
    initialise(getSwaptionVols); // do default settings
    // get the vol inside ICDSParSpreads curves if flag is set
    if (getCreditVols){
      setRetrievalMode(ICDSVol::TYPE,
		       ICDSParSpreads::TYPE,
		       true,
		       creditVolType);
    }
}

void MarketDataFetcherCDS::initialise(bool getSwaptionVols){
    // turn on currency basis
    setRetrievalMode(CurrencyBasis::TYPE, true, NULL);
    // hard code IRCalib::Model to 1 factor
    setRetrievalMode(IRCalib::Model::TYPE, true, IRCalib::getModelType(1));
    // don't get IRCalib::SmileBase inside quanto cds par spread curves
    setRetrievalMode(IRCalib::SmileBase::TYPE, QuantoCDSParSpreads::TYPE,
                     false, NULL);
    // get ICDSSpotVols inside quanto cds par spread curves when asked for
    // ICDSVol
    setRetrievalMode(ICDSVol::TYPE, QuantoCDSParSpreads::TYPE,
                     true, ICDSSpotVol::TYPE);
    // get SRMFXVol::VOL_TYPE inside quanto cds par spread curves when asked for
    // FXVolBase
    setRetrievalMode(FXVolBase::TYPE, QuantoCDSParSpreads::TYPE,
                     true, SRMFXVol::VOL_TYPE);
    // get either IRVol or IRCalib depending on flag (inside quanto curve)
    // when asked for an IRVol.
    // This essentially overrides the "don't get IR vols in yield curves" since
    // the outer object (here the QuantoCDSParSpreads) as precedence
    setRetrievalMode(IRVolBase::TYPE, QuantoCDSParSpreads::TYPE,
                     true, 
                     getSwaptionVols? IRVolCommon::TYPE: IRCalib::TYPE);
    if (getSwaptionVols){
        // then one slight twist, to make sure we get IRCalibs within IRVols:
        // Because the above setRetrievalMode will not apply to IRCalibs [when
        // getSwaptionVols is true] and as IRCalibs are IRVols so we
        // will pick up the default behaviour of not getting IRVols within
        // yield curves. We really should stop IRCalib being an IRVol ...
        // Anyway, start by saying don't get IRCalibs in general
        setRetrievalMode(IRCalib::TYPE, false, NULL);
        // and then do get them inside QuantoCDSParSpreads
        setRetrievalMode(IRCalib::TYPE, QuantoCDSParSpreads::TYPE, true, NULL);
    }
}


DRLIB_END_NAMESPACE

