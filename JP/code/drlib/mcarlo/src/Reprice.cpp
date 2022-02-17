//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Reprice.cpp
//
//   Description : Captures whether paths in monte carlo simulation can be 
//                 skipped when calculating sensitivities
//
//   Author      : Mark A Robson
//
//   Date        : 29 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetSpotGreek.hpp"
#include "edginc/MCPathGenerator.hpp"

DRLIB_BEGIN_NAMESPACE
IReprice::~IReprice(){}

/** Given the a ScalarShift (which contains an OutputName). This
    returns (|shift size|)/(1-|shift size|) *(spot of asset
    i)*(scaleFactors[i]), where asset i is the one which is being
    shifted because it is sensitive to OutputName. (The division by
    1-|shiftSize| is to allow for the supplied market data being
    tweaked already - it is conservative).If there is more than one
    asset sensitive to the OutputName, the sum of the above is used. If
    there are no assets sensitive, 0 is returned.  The SensControl
    must implement the IAssetSpotGreek interface */
double IReprice::calcMaxAssetShiftSize(const IMultiFactors* assets,
                                       const ScalarShift*   delta,
                                       const DoubleArray&   scaleFactors){
    static const string routine("IReprice::calcMaxAssetShiftSize");
    const IAssetSpotGreek* assetSpotGreek = 
        dynamic_cast<const IAssetSpotGreek*>(delta);
    if (!assetSpotGreek){
        throw ModelException(routine, "Method is only valid for shifts "
                             "implementing IAssetSpotGreek");
    }
    double shiftSize = fabs(delta->getShiftSize());
    double result;
    if (Maths::isZero(shiftSize)){
        result = 0.0;
    } else if (Maths::equals(shiftSize, 1.0)){
        throw ModelException(routine, "100% shift size not supported");
    } else {
        // find which assets are sensitive to it
        IntArray assetIdxs(assets->getSensitiveAssets(delta, false));
        // chain rule - so sum them
        double change = 0.0;
        for (int i = 0; i < assetIdxs.size(); i++){
            int iAsset = assetIdxs[i];
            change += assets->assetGetSpot(iAsset) *
                fabs(scaleFactors[iAsset]);
        }
        result = change * shiftSize/(1.0 - shiftSize);
    }
    return result;
}

/** Calculates threshold for quick greeks when doing a two sided shift
    (eg delta/gamma/cross gamma). If 
    B = sumProduct(maxDeltaScalingFactor, assets) then this returns the
    largest absolute amount that B will move by under the specified
    shifts */
double IReprice::baskTwoFactorGreekTolerance(
    const DoubleArray&                  maxDeltaScalingFactor,
    const IMultiFactors*                assets,
    const MCPathGenerator*              futurePathGen,
    const ScalarShiftArray&             sens){
    static const string method("IReprice::baskTwoFactorGreekTolerance");
    try{
        // We assume here that the sensitivity
        // is spot related - note that the calcMaxAssetShiftSize does
        // indeed validate for this so we are ok
        double tolerance = 0.0;
        DoubleArray deltaScaleFactors(assets->NbAssets());
        // we're going to calculate what the biggest change in value of
        // a single path is
        for (int iAsset = 0; iAsset < deltaScaleFactors.size(); iAsset++){
            // scale by 'delta' of forward and then by whatever scaling
            // factor the performance does
            deltaScaleFactors[iAsset] = 
                futurePathGen->maxDriftProduct(iAsset) * 
                maxDeltaScalingFactor[iAsset];
        }
        for (int j = 0; j < sens.size(); j++){
            const ScalarShift* shift = sens[j].get();
            // change in price can't be bigger than change in asset value *
            // participation/ref level
            double maxFactor = calcMaxAssetShiftSize(assets, shift,
                                                     deltaScaleFactors);
            // tolerance only relevant for2ndOrderGreek 
            tolerance += fabs(maxFactor);
        }
        return tolerance;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

IReprice::Vanilla::~Vanilla(){}

/** Constructor. 'store' should be set to true if the object needs
    to record values */
IReprice::Vanilla::Vanilla(int mode, int numPaths, double notional): 
    isStoringMode(mode != 0), notional(notional){
    if (isStoringMode){
        untweakedPrices = DoubleArraySP(new DoubleArray(0));
        untweakedPrices->reserve(numPaths);
    }
}
  
/** Returns a deep copy of this object. Storage data is shallow copied */
ISubReprice* IReprice::Vanilla::clone() const{
    Vanilla* copy = new Vanilla(0, 0, notional);
    copy->isStoringMode = isStoringMode;
    copy->tolerance = tolerance;
    copy->repriceForCmpt = ISubRepriceSP(!repriceForCmpt? 
                                         0: repriceForCmpt->clone());
    copy->untweakedPrices = untweakedPrices; // shallow copy
    return copy;
}

/** Same as Constructor above but treats reprice as a reprices object
    for the f in MAX(f-k,0) (otherwise f is assumed to have zero
    cross gamma eg sum of spots). Note that only a reference is
    taken to repriceForCmpt. The caller must ensure that the right
    methods (eg setForQuickGreeks() etc) are called at the right
    time for this object */
IReprice::Vanilla::Vanilla(const ISubRepriceSP& repriceForCmpt, 
                           int                  mode,
                           int                  numPaths,
                           double               notional): 
    isStoringMode(mode != 0), notional(notional),
    repriceForCmpt(repriceForCmpt){
    if (isStoringMode){
        untweakedPrices = DoubleArraySP(new DoubleArray(0));
        untweakedPrices->reserve(numPaths);
    }
}

/** Returns the repriceForCmpt as supplied in the constructor (may
    be null) */
ISubRepriceSP IReprice::Vanilla::getRepriceForCmpt(){
    return repriceForCmpt;
}

/** Sets object in mode ready to pricing for a [1st order] greek
    which is done by a one sided taylor expansion. Use NULL for params */
void IReprice::Vanilla::setForOneSidedGreek(
    IOneSidedGreekParamsSet* params){
    isStoringMode = false;
}

//// Primitive approach
IReprice::Vanilla::TwoSidedGreekParamsSet::TwoSidedGreekParamsSet(
    double tolerance): tolerance(tolerance){}

//// Does the work for you
IReprice::Vanilla::TwoSidedGreekParamsSet::TwoSidedGreekParamsSet(
            const DoubleArray&               maxCmptScalingFactor,
            const MCPathGenerator* futurePathGen,
            const IMultiFactors*                   assets, 
            const ScalarShiftArray&          sens){
    tolerance = baskTwoFactorGreekTolerance(maxCmptScalingFactor, assets,
                                             futurePathGen, sens);
}

IReprice::Vanilla::TwoSidedGreekParamsSet::~TwoSidedGreekParamsSet(){}

/** Sets object in mode ready to pricing for a [1st or 2nd order]
    greek which is done by a two sided taylor expansion. Use
    TwoSidedGreekParamsSet for an instance of
    ITwoSidedGreekParamsSet */
void IReprice::Vanilla::setForTwoSidedGreek(
    ITwoSidedGreekParamsSet* params){
    TwoSidedGreekParamsSet& myParams =
        dynamic_cast<TwoSidedGreekParamsSet&>(*params);
    isStoringMode = false;
    this->tolerance = myParams.tolerance;
}    

/** Sets object in mode ready to pricing for cross gamma.  Use
    TwoSidedGreekParamsSet for an instance of
    ITwoSidedGreekParamsSet */
void IReprice::Vanilla::setForXGamma(
    IXGammaParamsSet* params){
    TwoSidedGreekParamsSet& myParams =
        dynamic_cast<TwoSidedGreekParamsSet&>(*params);
    setForTwoSidedGreek(&myParams);
}

/** Stores the supplied value for the current path on initial
    pricing run and this allows methods below to be
    implemented. The method must be invoked for each simulated
    path and called in order. The price should be before the
    MAX with zero is applied and without any scaling by notional */
void IReprice::Vanilla::store(double preMaxPrice){
    if (isStoringMode){
        untweakedPrices->push_back(preMaxPrice);
    }
}

/** Returns true if original value < 0.0 */
bool IReprice::Vanilla::firstDerivZero(int path) const{
    double origVal = (*untweakedPrices)[path];
    if (origVal <= 0.0){
        return true;
    }
    // if have repriceForCmpt, then by chain rule ...
    if (repriceForCmpt.get()){
        return repriceForCmpt->firstDerivZero(path);
    }
    return false;
}

/** Returns true if original value <= -tolerance. */
bool IReprice::Vanilla::firstNumericalDerivZero(int path) const{
    double origVal = (*untweakedPrices)[path];
    if (origVal <= -tolerance){
        return true;
    }
    // if have repriceForCmpt, then by chain rule ...
    if (repriceForCmpt.get()){
        return repriceForCmpt->firstNumericalDerivZero(path);
    }
    return false;
}
   
/** Returns true if original value <= -tolerance or 
    original value >= tolerance */
bool IReprice::Vanilla::crossNumericalDerivZero(int path) const{
    double origVal = (*untweakedPrices)[path];
    if (origVal <= -tolerance){
        return true; // out of the money
    }
    if (origVal < tolerance){
        return false; // at the money
    }
    // in the money - cross deriv zero unless repriceForCmpt isn't
    if (repriceForCmpt.get()){
        return repriceForCmpt->crossNumericalDerivZero(path);
    }
    return true;
}

/** Returns the original value (ie on pricing run) for specified path */
double IReprice::Vanilla::originalPrice(int path) const{
    return (notional * Maths::max((*untweakedPrices)[path], 0.0));
}

IReprice::Spread::~Spread(){}

/** Constructor. mode is as per IQuickGreeksCommon::createOrigPrices().
    payoff=(isCall?1: -1)*notional*(MAX(P1,0)-MAX(P2, 0)) with 
    Pi = (isCall?1: -1)*(f-ki). k1 = loStrike < k2=hiStrike. If not doing
    crossGamma then the important part is that the optionality
    lies around k1 and k2 (the overall price being
    irrelevant). The repriceForCmpt is optional and can be
    null. If present it is treated as a reprices object for the f
    in MAX(f-k,0) (otherwise f is assumed to have zero cross gamma
    eg sum of spots). Note that only a reference is taken to
    repriceForCmpt. The caller must ensure that the right methods
    (eg setForQuickGreeks() etc) are called at the right time for
    this object */
IReprice::Spread::Spread(const ISubRepriceSP& repriceForCmpt, 
                         int                  mode, 
                         int                  numPaths,
                         bool                 isCall, 
                         double               notional,
                         double               loStrike,
                         double               hiStrike):
    isStoringMode(mode != 0), isCall(isCall), notional(notional),
    loStrike(loStrike), hiStrike(hiStrike), repriceForCmpt(repriceForCmpt){
    if (isStoringMode){
        basketPrices = DoubleArraySP(new DoubleArray(0));
        basketPrices->reserve(numPaths);
    }
}

/** Returns a deep copy of this object. Storage data is shallow copied */
ISubReprice* IReprice::Spread::clone() const{
    ISubRepriceSP subCopy(!repriceForCmpt? 0: repriceForCmpt->clone());
    Spread* copy = new Spread(subCopy, 0, 0, isCall, notional, 
                              loStrike, hiStrike);
    copy->isStoringMode = isStoringMode;
    copy->tolerance = tolerance;
    copy->basketPrices = basketPrices; // shallow copy
    return copy;
}

/** Sets object in mode ready to pricing for a [1st order] greek
    which is done by a one sided taylor expansion. Use NULL for params */
void IReprice::Spread::setForOneSidedGreek(
    IOneSidedGreekParamsSet* params){
    isStoringMode = false;
}

/** Sets object in mode ready to pricing for a [1st or 2nd order]
    greek which is done by a two sided taylor expansion. Use
    TwoSidedGreekParamsSet for an instance of
    ITwoSidedGreekParamsSet */
void IReprice::Spread::setForTwoSidedGreek(
    ITwoSidedGreekParamsSet* params){
    Vanilla::TwoSidedGreekParamsSet& myParams =
        dynamic_cast<Vanilla::TwoSidedGreekParamsSet&>(*params);
    isStoringMode = false;
    this->tolerance = myParams.tolerance;
}

/** Sets object in mode ready to pricing for cross gamma.  Use
    TwoSidedGreekParamsSet for an instance of
    ITwoSidedGreekParamsSet */
void IReprice::Spread::setForXGamma(IXGammaParamsSet* params){
    Vanilla::TwoSidedGreekParamsSet& myParams =
        dynamic_cast<Vanilla::TwoSidedGreekParamsSet&>(*params);
    setForTwoSidedGreek(&myParams);
}

/** Stores the supplied value for the current path on initial
    pricing run and this allows methods below to be
    implemented. The method must be invoked for each simulated
    path and called in order. The basketValue supplied should be as
    described in the constructor (ie the variable f) */
void IReprice::Spread::store(double basketPrice){
    if (isStoringMode){
        basketPrices->push_back(basketPrice);
    }
}    

/** Returns the repriceForCmpt as supplied in the constructor (may
    be null) */
ISubRepriceSP IReprice::Spread::getRepriceForCmpt(){
    return repriceForCmpt;
}

/** Returns true if path lies out of the money or is above the
    hiStrike */
bool IReprice::Spread::firstDerivZero(int path) const{
    double origVal = (*basketPrices)[path];
    if (origVal <= loStrike || origVal >= hiStrike){
        return true;
    }
    // if have repriceForCmpt, then by chain rule ...
    if (repriceForCmpt.get()){
        return repriceForCmpt->firstDerivZero(path);
    }
    return false;
}

/** Returns true if path is below loStrike-tolerance or is above
    the hiStrike+tolerance */
bool IReprice::Spread::firstNumericalDerivZero(int path) const{
    double origVal = (*basketPrices)[path];
    if (origVal <= loStrike-tolerance || origVal >= hiStrike+tolerance){
        return true;
    }
    // in the money - cross deriv zero unless repriceForCmpt isn't
    if (repriceForCmpt.get()){
        return repriceForCmpt->firstNumericalDerivZero(path);
    }
    return false;
}    
   
/** Returns true if path lies away from loStrike and hiStrike (using
    tolerance) */
bool IReprice::Spread::crossNumericalDerivZero(int path) const{
    double origVal = (*basketPrices)[path];
    if (origVal <= loStrike-tolerance || origVal >= hiStrike+tolerance){
        return true; // out of the money
    }
    if (origVal < loStrike+tolerance || origVal > hiStrike-tolerance){
        return false; // at point of optionality
    }
    // 'in the money' - cross deriv zero unless repriceForCmpt isn't
    if (repriceForCmpt.get()){
        return repriceForCmpt->crossNumericalDerivZero(path);
    }
    return true;
}

/** Returns the original value (ie on pricing run) for specified path */
double IReprice::Spread::originalPrice(int path) const{
    double origVal = (*basketPrices)[path];
    double opt1 = origVal-loStrike;
    double opt2 = origVal-hiStrike;
    if (!isCall){
        opt1 = -opt1;
        opt2 = -opt2;
    }
    double price = notional*(Maths::max(opt1, 0.0)-Maths::max(opt2, 0.0));
    if (!isCall){
        price = -price;
    }
    return price;
}

DRLIB_END_NAMESPACE
