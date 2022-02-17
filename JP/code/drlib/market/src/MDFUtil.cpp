//
//   Author      : Mark A Robson
//
//   Date        : 7th August 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/LocalCorrSqueeze.hpp"
#include "edginc/CorrelationTerm.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/VarSwapBasis.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sets whether correlation skew shall be used or not. */
void MDFUtil::setCorrSkewMode(MarketDataFetcher& mdf,
                              bool               setUseCorrSkew){
    mdf.setRetrievalMode(CorrelationCategory::TYPE, setUseCorrSkew, NULL);    
    mdf.setRetrievalMode(LocalCorrSqueeze::TYPE, setUseCorrSkew, NULL);
}
    
/** Retrieves whether correlation skew is being used */
bool MDFUtil::useCorrSkew(const MarketDataFetcher& mdf){
    if (mdf.getRetrievalMode(LocalCorrSqueeze::TYPE)) {
        return true;
    } else {
        return false; 
    }
}

/** Sets whether correlation term structure shall be used or not. */
void MDFUtil::setCorrTermMode(MarketDataFetcher& mdf,
                              bool               setUseCorrTermStructure){
    mdf.setRetrievalMode(CorrelationTerm::TYPE,
                         setUseCorrTermStructure, NULL);
    mdf.setRetrievalMode(CorrelationCategory::TYPE,
                         setUseCorrTermStructure, NULL);
}
    
/** Retrieves whether correlation term structure is being used */
bool MDFUtil::useCorrTerm(const MarketDataFetcher& mdf){
    return mdf.getRetrievalMode(CorrelationTerm::TYPE);
}

/** Specific models can request the support of currency basis, if
    it is supplied */
void MDFUtil::setUseCurrencyBasis(MarketDataFetcher& mdf,
                                  bool               flag){
    mdf.setRetrievalMode(CurrencyBasis::TYPE, flag, NULL);
}

/** Sets whether VarSwapBasis objects are retrieved or not */
void MDFUtil::setVarSwapBasis(MarketDataFetcher& mdf,
                              bool               flag){
    mdf.setRetrievalMode(VarSwapBasis::TYPE, flag, NULL);
}

/** Configure MDF so that only IR Swaption vols are retrieved ie 
    redirects IRVolBase to IRVol and does not retrieve any IRCalib
    type parameters */
void MDFUtil::setUseSimpleIRVol(MarketDataFetcher& mdf){
    // just get swaption vols when asked for ir vols.
    mdf.setRetrievalMode(IRVolBase::TYPE, IYieldCurve::TYPE, true, IRVolCommon::TYPE);
    // and don't get any IRCalib stuff
    mdf.setRetrievalMode(IRCalib::SmileBase::TYPE, false, NULL);
    mdf.setRetrievalMode(IRCalib::Model::TYPE, false, NULL);
    mdf.setRetrievalMode(IRCalib::TYPE, false, NULL);
}

DRLIB_END_NAMESPACE

