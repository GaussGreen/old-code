//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StruckEquity.cpp
//
//   Description : Implement asset for ccy struck assets 
//
//   Author      : Mark A Robson
//
//   Date        : 16 Oct 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/StruckAsset.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/PDFDefaultLNStrike.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/ATMVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE
StruckAsset::~StruckAsset(){}

/** Constructor needed for case when instrument specifies
    underlying asset and currency treatment */
StruckAsset::StruckAsset(const string& name,         // name for struck asset
                         const string& assetName,    // asset to make struck
                         const string& fxAssetName): // name for fx asset
    CAsset(TYPE), name(name), asset(assetName), fx(fxAssetName),
    getCorrFromCache(true){
    validatePop2Object();
}

/** Validation */
void StruckAsset::validatePop2Object(){
    static const string method("StruckAsset::validatePop2Object");
    // Since need fx but may not have yet called StruckAsset::getMarket
    // need to shield 
    if (fx.get()) {
        if (name.empty()){
            // use the base ccy's name not its iso code
            name = CAsset::getSyntheticName(asset.getName(),
                                            CAsset::CCY_TREATMENT_STRUCK,
                                            fx->getYCName());
        } else if (name == asset.getName()){
            throw ModelException(method, 
                                 "Struck asset's name must be different to"
                                 " underlying asset's name");
        }
        if (asset.get()){
            if (fx->getBaseCcyIsoCode() == AssetUtil::assetCcy(asset.get())){
                throw ModelException(method,
                                     "Struck currency is identical to "
                                     "underlying currency (" +
                                     fx->getBaseCcyIsoCode() + ")");
            }
        }
    }
    getCorrFromCache = !corrEqFx;
    if (!getCorrFromCache){
        // hack for some tree products which incorrectly do not always
        // call getMarket. We'll guess the types for now - if getMarket
        // is called it will be done correctly
        corrEqFx->configureForSensitivities(SimpleEquity::TYPE, FXAsset::TYPE);
    }
    if (asset.get() && !Asset::IStruckable::TYPE->isInstance(asset.get())){
        throw ModelException(method, "Asset to be struck must implement the "
                             "Asset::IStruckable interface");
    }
}

/** Pull out the asset date and correlation from the market data */
void StruckAsset::getMarket(const IModel* model, const MarketData* market){
    // get the data we need from the cache
    market->GetReferenceDate(baseDate); // get/check ref date

    // essentially do asset.getData(model, market) but tell model that the
    // ccy is changing
    CAsset::getAssetInNewCurrency(model, market, asset);
    // see if want flat fx vol here - NB ought to pass asset->getYCName() down
    FXAsset::getMarketForStruck(fx, model, market);
    if (getCorrFromCache){
        string corrName = market->getCorrelationName(asset->getName(),
                                                     fx->getName());
        corrEqFx = CorrelationSP::dynamicCast(
            model->getCorrelation(corrName, asset->getClass(), fx->getClass(),
                                  Correlation::TYPE, market));
    } else {
        // configure it correctly ourselves
        corrEqFx->configureForSensitivities(asset->getClass(), fx->getClass());
    }
    
    // finally, validate the object again
    validatePop2Object();
}

/** returns the spot price */
double StruckAsset::getSpot() const{
    return asset->getSpot() * fx->getSpot();
}

/** returns the asset name */
string StruckAsset::getName() const{
    return name;
}

/** returns the equity's name */
string StruckAsset::getTrueName() const{
    return asset->getTrueName();
}

/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
CVolProcessed * StruckAsset::getProcessedVol(
    const CVolRequest* volRequest) const {
    const Asset::IStruckable& struckable = 
        dynamic_cast<const Asset::IStruckable&>(*asset.get());
    return struckable.getProcessedVol(volRequest, this, 
                                      fx.get(), corrEqFx.get());
}

/** Calculates the expected spot price of the asset at the given date */
double StruckAsset::fwdValue(const DateTime& date) const{
    return asset->fwdValue(date) * fx->fwdValue(date);
}

void StruckAsset::fwdValue(const DateTimeArray&     dates,
                           const FwdValueAlgorithm& algo,
                           CDoubleArray&            result) const{

    if (dates.size() != result.size()) {
        throw ModelException("StruckAsset::fwdValue", "Date list and result"
                             " array have different sizes");
    }

    asset->fwdValue(dates, algo, result); // do underlying
    CDoubleArray    fxFwds(dates.size());
    fx->fwdValue(dates, fxFwds); // do fx

    for (int i = 0; i < dates.size(); i++){
        result[i] *= fxFwds[i];
    }
}

void StruckAsset::fwdValue(const DateTimeArray& dates,
                           CDoubleArray&        result) const{
    FwdValueAlgorithm algo(false);
    fwdValue(dates, algo, result);
}

/** Returns the name (not the ISO code) of the asset ccy */
string StruckAsset::getYCName() const {
    return fx->getYCName();
}

/** Calculate the settlement date associated with a given trade date */
DateTime StruckAsset::settleDate(const DateTime& tradeDate) const{
    return asset->settleDate(tradeDate);
}

/** (IStruck interface) Returns the fx spot */
double StruckAsset::getFXSpot() const{
    return fx->getSpot();
}

/** (IStruck interface) Returns the fx forward value */
double StruckAsset::fxFwdValue(const DateTime& date) const{
    return fx->fwdValue(date);
}

/** Shifts the object using given shift. */
bool StruckAsset::sensShift(Theta* shift){
    baseDate = shift->rollDate(baseDate);
    return true; 
}

/** returns sensitive strikes for a given vol request */
void StruckAsset::getSensitiveStrikes(
    const CVolRequest*               volRequest,
    OutputNameConstSP                outputName,
    const SensitiveStrikeDescriptor& sensStrikeDesc,
    DoubleArraySP                    sensitiveStrikes) const
{
    // to do: alter getSensitiveStrikes signature
    const CVolRequestLN& request = 
        dynamic_cast<const CVolRequestLN&>(*(IObject*)(volRequest));
    if (!sensStrikeDesc.forwardOnly) {
        // must turn strikes into struck asset's ccy
        double fxSpot = fx->getSpot();
        if (Maths::isPositive(fxSpot)) {
            if ( outputName->equals(fx->getVolName()) ) {
                // Create an atm vol request for the FX side of things
                CVolRequestLNSP atmRequest(new ATMVolRequest());
                atmRequest->getSensitiveStrike(fx->getSpot(),
                                               sensitiveStrikes);
            } else {
                // Handle the asset part of the struck underlying
                // NB. There are no explicit name checks because the underlying
                // may be a composite, so assets must check for themselves.
                CVolRequestLNSP undCcyVolRequest(copy(&request));
                undCcyVolRequest->scale(1.0/fxSpot);
                // then just pass down to asset
                asset->getSensitiveStrikes(undCcyVolRequest.get(), outputName, 
                                           sensStrikeDesc, sensitiveStrikes);
            }
        } else {
            throw ModelException("StruckAsset::getSensitiveStrikes", 
                                 "FX rates must be strictly positive");
        }
    }
}

/** overrides CObject version to allow for easy default */
bool StruckAsset::accept(ICollector* collector) const{
    if (!CClass::invokeAcceptMethod(this, collector)){
        // if no method registered default by passing onto asset
        return asset->accept(collector);
    }
    return false;
}

void StruckAsset::acceptValueDateCollector(const StruckAsset*   asset, 
                                           CValueDateCollector* collector){
    collector->valueDateValidate(asset->baseDate, asset->getName());
}

void StruckAsset::acceptImntCcy(const StruckAsset* asset,
                                AssetCcyCollector* collector){
    collector->currencyValidate(asset->fx->getBaseCcyIsoCode(),
                                asset->getName());
    // then make sure our asset that we're making struck is ok
    AssetCcyCollector collect;
    asset->asset->accept(&collect);
}

void StruckAsset::acceptNameCollector(
    const StruckAsset* asset, AssetNameCollector* collector){
    collector->assetNameValidate(asset->getName());
    asset->asset->accept(collector);
}

/** given a current spot level, get the next strike on the vol surface where 
    the slope is non-differentiable */
double StruckAsset::getNextStrike(const double& strike,
                                  bool          isUp,
                                  bool&         offSurface) const {
    double adjStrike = strike;
    double fxSpot = fx->getSpot();
    if (!Maths::isPositive(fxSpot) )
    {
        throw ModelException("StruckAsset::getNextStrike", 
                             "FX rates must be strictly positive");
    }
    adjStrike /= fxSpot;

    const IObject* obj = asset.get();
    const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(obj);
    double volStrike;
    if (nextStrike) {
        volStrike = nextStrike->getNextStrike(adjStrike, isUp, offSurface);
    } else {
        volStrike  = 0.0;
        offSurface = true;
    }
    // volStrike is in underlyer's ccy, but we want it in struck ccy
    return (volStrike * fxSpot);
}

PDFCalculator* StruckAsset::pdfCalculator(const PDFRequest* request) const {
    const PDFRequestLNStrike* lnRequest =
        dynamic_cast<const PDFRequestLNStrike*>(request);
    if (lnRequest){
        return new PDFDefaultLNStrike(baseDate, this, lnRequest);
    }
    throw ModelException("StruckAsset::pdfCalculator",
                         "Request of type "+request->getClass()->getName()+
                         " not supported");
}

/** returns dividend list - to be retired.
    Needed because of crap method in AssetUtil that casts around desperately
    looking for a type that it knows about */
    /** returns dividend list - to be retired - it's a bit meaningless */
DividendListSP  StruckAsset::getAllDivsBetweenDates(const DateTime& start,
                                                    const DateTime& end) const{
    return AssetUtil::getAllDivsBetweenDates(asset.get(), start, end);
}

/** Returns fair value of asset */
double StruckAsset::fairValue() const {
    const IAssetFairValue* fv = dynamic_cast<const IAssetFairValue*>(asset.get());
    if (!fv) {
        throw ModelException("StruckAsset::fairValue", 
                             "asset (" + asset->getTrueName() + 
                             ") doesn't support fair value calculation");
    }
    return fv->fairValue() * fx->getSpot();;
}

// Tactical method to get (copy of) underlying (non-struck) asset 
CAssetConstSP  StruckAsset::getPlainAsset() const {
    return CAssetConstSP(dynamic_cast<const CAsset *>(asset.get()));
}

//// Returns the ccy treatment for the asset. Default returns 
string StruckAsset::getCcyTreatment() const{
    return CCY_TREATMENT_STRUCK;
}

void StruckAsset::checkStruckSourceAndObsType(const ObservationSource* source,
                                              const ObservationType* type) const {
    // throw an exception if the source/type is not suitable for a struck asset
    // i.e. it should have 2 sources for asset and fx
    if(!(source->getClass() == CClass::forName("StruckObservationSource"))) {
        throw ModelException("StruckAsset::checkStruckSourceAndObsType",
                    "Inconsistent source for a struck asset. Struck assets "
                    "require 2 sources - one for the asset and one for the FX");
    }
    if(!(type->getClass() == CClass::forName("StruckObservationType"))) {
        throw ModelException("StruckAsset::checkStruckSourceAndObsType",
                    "Inconsistent observation type for a struck asset. Struck assets "
                    "require 2 types - one for the asset and one for the FX");
    }
}

double StruckAsset::pastValue(const DateTime&             sampleDate,
                              const ObservationType*      obsType,
                              const ObservationSource*    source,
                              const FixingType*           fixType,
                              const IObservationOverride* overrides,
                              const SamplingConvention*   sampleRule) const{
    // first try and get the sample from the overrides
    // overrides will be for the struck asset levels
    if (overrides) {
        double pValue = 0.0;
        if(overrides->value(sampleDate, obsType, &pValue)) {
            return pValue;
        }
    }                                   
    // get equity sample and FX sample and multiply
    try {
        checkStruckSourceAndObsType(source, obsType);
        // note after checking we know these dynamic casts will work
        const StruckObservationSource* compoundSource = 
                dynamic_cast<const StruckObservationSource*>(source);
        const StruckObservationType* compoundType = 
                dynamic_cast<const StruckObservationType*>(obsType);

        //this assumes the sample will be available for the FX when the asset is available
        double assetValue = asset->pastValue(sampleDate, compoundType->getPrimaryObsType(), 
                                            source, fixType, 0, sampleRule);

        double fxValue = fx->pastValue(sampleDate, compoundType->getFXObsType(), 
                                       compoundSource->getFXSource(),
                                       fixType, 0, sampleRule);
        return assetValue * fxValue;
    } catch (ModelException& e) {
        throw ModelException(e, "StruckAsset::pastValue", 
                    "Failed to retrieve past sample for struck asset " + 
                    getTrueName());
    }
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool StruckAsset::observationDate(const DateTime&           sampleDate,
                                  const ObservationSource*  source,
                                  const SamplingConvention* sampleRule,
                                  DateTime*                 obsDate) const {
    // Currently assume the FX and the asset move together based on the asset
    return asset->observationDate(sampleDate, source, sampleRule, obsDate);
}

// the IMarketObservable interface for retrieving past samples events
double StruckAsset::addPastSampleEvent(const DateTime&      sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const {
    static const string method("StruckAsset::addPastSampleEvent");
    try {
        // first try and get the sample from the overrides
        if (overrides) {
            double level = 0.0;
            if (overrides->addPastSampleEvent(sampleDate, obsType, source,
                                            fixType, collector, level)) {
                return level;
            }
        }                                   
        checkStruckSourceAndObsType(source, obsType);
        // note after checking we know these dynamic casts will work
        const StruckObservationSource* compoundSource = 
                dynamic_cast<const StruckObservationSource*>(source);
        const StruckObservationType* compoundType = 
                dynamic_cast<const StruckObservationType*>(obsType);

        DateTime obsDate;
        if (!observationDate(sampleDate, source, sampleRule, &obsDate)) {
            throw ModelException(method, "Internal error "
                                "Should not attempt to get sampling "
                                "event for omitted date");
        }

        AssetFixTypeSP assetFixType = 
                AssetFixTypeSP(new AssetFixType(asset->getName()));
        double level = asset->addPastSampleEvent(sampleDate, 
                                                    compoundType->getPrimaryObsType(), 
                                                    source,
                                                    assetFixType.get(), 0, 
                                                    sampleRule, collector);
        AssetFixTypeSP fxFixType = 
                AssetFixTypeSP(new AssetFixType(fx->getName()));
        level *= fx->addPastSampleEvent(sampleDate, compoundType->getFXObsType(), 
                                    compoundSource->getFXSource(),
                                    fxFixType.get(), 0, sampleRule, collector);
        // then add the struck asset sample event
        collector->addSample(sampleDate, obsDate, source, obsType,
                            level, fixType, false);
        return level;
    } catch (exception& e) {
        throw ModelException(e, method, "Failed to retrieve past sample event "
                                        "for asset " + getName() + 
                                        " for sample date (" +
                                        sampleDate.toString() + ")");
    }
}

// the IMarketObservable interface for 
// is the given date a holiday for the relevant source
bool StruckAsset::isHoliday(const DateTime&             sampleDate,
                             const ObservationSource*   source) const {
    // Currently assume the FX and the asset move together based on the asset
    return asset->isHoliday(sampleDate, source);
}

/* for reflection */
StruckAsset::StruckAsset(): CAsset(TYPE), getCorrFromCache(true) {}
    
class StruckAssetHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(StruckAsset, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(Asset::IStruck);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(INextStrike);
        IMPLEMENTS(IAssetFairValue);
        EMPTY_SHELL_METHOD(defaultStruckAsset);
        FIELD(name, "Asset's name");
        FIELD_MAKE_OPTIONAL(name);
        FIELD(asset, "Asset to make struck");
        FIELD(fx, "What the stock is struck into");
        FIELD(corrEqFx, "correlation between equity and fx");
        FIELD_MAKE_OPTIONAL(corrEqFx);
        FIELD(baseDate,    "Value Date");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(getCorrFromCache, "cached flag");
        FIELD_MAKE_TRANSIENT(getCorrFromCache); // hide from dd interface
        ClassSetAcceptMethod(StruckAsset::acceptValueDateCollector);
        ClassSetAcceptMethod(StruckAsset::acceptImntCcy);
        ClassSetAcceptMethod(StruckAsset::acceptNameCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const StruckAsset* asset,
                                        DividendCollector* collector){
        collector->addDivs(asset->asset.get(), asset->fx.get());
    }

    static IObject* defaultStruckAsset(){
        return new StruckAsset();
    }
};

CClassConstSP const StruckAsset::TYPE = CClass::registerClassLoadMethod(
    "StruckAsset", typeid(StruckAsset), StruckAssetHelper::load);

DRLIB_END_NAMESPACE
