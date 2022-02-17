//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StruckEquity.cpp
//
//   Description : Implement asset for simple ccy struck equity 
//
//   Author      : Mark A Robson
//
//   Date        : 8 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/Maths.hpp"
#include "edginc/PDFDefaultLNStrike.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/VolDeltaShiftSize.hpp"

DRLIB_BEGIN_NAMESPACE

StruckEquity::~StruckEquity(){}

/** Constructor needed for case when instrument specifies underlying asset and 
    currency treatment */
StruckEquity::StruckEquity(const string& name, // name for struck asset
                           const Equity* equity,
                           const string& volName, // name for vol
                           const string& fxAssetName) // name for fx asset
    :CAsset(TYPE), name(name), equity(copy(equity)), vol(volName),
     fx(fxAssetName), getCorrFromCache(true){
    validatePop2Object();
}

StruckEquity::StruckEquity(const string& name, // name for struck asset
                           const Equity* equity,
                           const CVolBaseWrapper    vol, 
                           const FXAssetWrapper     fx,
                           const Correlation*       corrEqFx) 
    :CAsset(TYPE), name(name), equity(copy(equity)), vol(vol),
     fx(fx), corrEqFx(copy(corrEqFx)){
}
                           
/** Validation */
void StruckEquity::validatePop2Object(){
    static const string method("StruckEquity::validatePop2Object");
    // Since need fx but may not have yet called StruckEquity::getMarket
    // need to shield 
    if (fx.get()) {
        if (name.empty()){
            // use the base ccy's name not its iso code
            name = CAsset::getSyntheticName(equity->getName(),
                                            CAsset::CCY_TREATMENT_STRUCK,
                                            fx->getYCName());
        } else if (name == equity->getName()){
            throw ModelException(method, 
                                 "Struck asset's name must be different to"
                                 " equity's name");
        }
        if ( equity->getYCIsoCode() == 
             fx->getBaseCcyIsoCode()) {
            throw ModelException(method,
                                 "Struck currency is identical to "
                                 "underlying currency (" +
                                 equity->getYCIsoCode() + ")");
        }
    }
    getCorrFromCache = !corrEqFx;
    if (!getCorrFromCache){
        // hack for some tree products which incorrectly do not always
        // call getMarket. We'll guess the types for now - if getMarket
        // is called it will be done correctly
        corrEqFx->configureForSensitivities(SimpleEquity::TYPE, FXAsset::TYPE);
    }
}

/** Pull out the vol, fx asset and correlation from the
    market data */
void StruckEquity::getMarket(const IModel* model, const MarketData* market){
    // need to identify new 'domestic' currency
    const string& ycName = equity->getYCName();

    // get the data we need from the cache
    vol.getData(model, market, ycName);
    // see if  want flat fx vol here - NB ought to pass ycName down
    FXAsset::getMarketForStruck(fx, model, market);
    if (getCorrFromCache){
        string corrName = market->getCorrelationName(equity->getName(),
                                                     fx->getName());
        corrEqFx = CorrelationSP::dynamicCast(
            model->getCorrelation(corrName, SimpleEquity::TYPE, fx->getClass(),
                                  Correlation::TYPE, market));
    } else {
        // configure it correctly ourselves
        corrEqFx->configureForSensitivities(SimpleEquity::TYPE, fx->getClass());
    }
    
    // then get data that we already have to fill in missing bits
    //  - NB ought to pass ycName down
    equity->getMarket(model, market);

    // finally, validate the object again
    validatePop2Object();
}

/** get a (const) plain asset */
CAssetConstSP  StruckEquity::getPlainAsset() const{
    CAssetConstSP asset = CAssetConstSP(new SimpleEquity(equity.get(), vol.get()));
    return asset;
}

/** returns the spot price */
double StruckEquity::getSpot() const{
    return equity->spot() * fx->getSpot();
}

/** returns the asset name */
string StruckEquity::getName() const{
    return name;
}

/** returns the equity's name */
string StruckEquity::getTrueName() const{
    return equity->getName();
}

/** returns the name of the vol base object */
string StruckEquity::getVolName() const {
    return vol.getName();
}

/** Returns fair value of asset */
double StruckEquity::fairValue() const {
    return equity->fairValue() * fx->getSpot();;
}
   
/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
CVolProcessed * StruckEquity::getProcessedVol(
    const CVolRequest* volRequest) const
{
    if (!corrEqFx){
        throw ModelException("StruckEquity::getProcessedVol", "No correlation"
                             " supplied");
    }
    return vol->getProcessedVol(volRequest, this, fx.get(), corrEqFx.get());
}

/** Calculates the expected spot price of the asset at the given date */
double StruckEquity::fwdValue(const DateTime& date) const{
    return equity->fwdValue(date) * fx->fwdValue(date);
}

void StruckEquity::fwdValue(const DateTimeArray&     dates,
                            const FwdValueAlgorithm& algo,
                            CDoubleArray&            result) const{

    if ( dates.size() != result.size() )
    {
        throw ModelException("StruckEquity::fwdValue", 
            "date list and result array have different sizes");
    }
    CDoubleArray    equityFwds(dates.size());
    CDoubleArray    fxFwds(dates.size());

    equity->fwdValue(dates, algo, equityFwds);
    fx->fwdValue(dates, fxFwds);

    for (int i = 0; i < dates.size(); i++){
        result[i] = equityFwds[i] * fxFwds[i];
    }
}

void StruckEquity::fwdValue(const DateTimeArray& dates,
                            CDoubleArray&        result) const{
    FwdValueAlgorithm algo(false);
    fwdValue(dates, algo, result);
}

/** Returns the name (not the ISO code) of the asset ccy */
string StruckEquity::getYCName() const {
    return fx->getYCName();
}

// /** Calculates the expected spot price of the asset at the given date if
//         the spot price had the given value spot on spotDate */
// double StruckEquity::fwdFwd(const DateTime& spotDate,
//                             double          spot, 
//                             const DateTime& fwdDate) const{
//     throw ModelException("StruckEquity::fwdFwd", "fwdFwd calculation is not allowd for struck equities");
//     return 0;
// }

/** Calculate the settlement date associated with a given trade date */
DateTime StruckEquity::settleDate(const DateTime& tradeDate) const{
    return equity->settles(tradeDate);
}

/** (IStruck interface) Returns the fx spot */
double StruckEquity::getFXSpot() const{
    return fx->getSpot();
}

/** (IStruck interface) Returns the fx forward value */
double StruckEquity::fxFwdValue(const DateTime& date) const{
    return fx->fwdValue(date);
}

/** return the VolBaseWrapper */
const CVolBaseWrapper& StruckEquity::getVol()const
{
    return vol;
}

/** Is the Eq/FX correlation in the cache ? */
bool StruckEquity::useCorrFromCache()const
{
    return getCorrFromCache;
}

/** return the Eq/FX correlation */
const Correlation* StruckEquity::getCorrelation()const
{
    return corrEqFx.get();
}

/** return the FXVolBaseWrapper */
const FXAssetWrapper& StruckEquity::getFXAsset()const
{
    return fx;
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string StruckEquity::sensName(VegaSkewParallel* shift) const{
    return vol->getName();
}
    
bool StruckEquity::sensShift(VegaSkewParallel* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string StruckEquity::sensName(VegaSkewPointwise* shift) const{
    return vol->getName();
}
    
ExpiryArrayConstSP StruckEquity::sensExpiries(VegaSkewPointwise* shift) const{
    return ExpiryArrayConstSP(); // return null SP - the vol will do this
}

bool StruckEquity::sensShift(VegaSkewPointwise* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string StruckEquity::sensName(DeltaSurface* shift) const{
    return equity->getName();
}
    
bool StruckEquity::sensShift(DeltaSurface* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot(), equity->getName(), vol->getName());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string StruckEquity::sensName(VolRelativeShift* shift) const{
    return VolRelativeShift::IShift::TYPE->isInstance(vol.get())
        ? vol->getName(): "";
}
    
bool StruckEquity::sensShift(VolRelativeShift* shift){
    // store the underlying spot price - want the real equity
    // not its struck cousin
    SimpleEquity simple(equity.get(), vol.get());
    shift->setSpot(equity->getValueDate(), &simple);
    return true; // continue on and do the vol
}

/** returns sensitive strikes for a given vol request */
void StruckEquity::getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const
{
    // to do: alter getSensitiveStrikes signature
    const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
        *(IObject*)(volRequest));
    if (! sensStrikeDesc.forwardOnly) {
        double fxSpot = fx->getSpot();

        if (Maths::isPositive(fxSpot)) {
            if (outputName->equals(fx->getVolName())) {
                // Create an atm vol request for the FX side of things
                CVolRequestLNSP atmRequest(new ATMVolRequest());
                atmRequest->getSensitiveStrike(fx->getSpot(),
                                               sensitiveStrikes);
            } else {
                // Handle the equity part of the struck equity
                // NB. There are no explicit name checks because the equity
                // may have proxy vol, so assets must check for themselves.
                CVolRequestLNSP undCcyVolRequest(copy(&request));
                undCcyVolRequest->scale(1.0/fxSpot);
                getPlainAsset()->getSensitiveStrikes(undCcyVolRequest.get(), outputName, 
                                           sensStrikeDesc, sensitiveStrikes);
            }
        } else {
            throw ModelException("StruckEquity::getSensitiveStrikes", 
                                 "FX rates must be strictly positive");
        }
    }
}

void StruckEquity::acceptDeltaShift(const StruckEquity* asset, 
                                    ShiftSizeCollector* collector)
{
    const IObject* obj = asset->vol.get();
    const IVolDeltaShiftSize* volDelShift = 
                        dynamic_cast<const IVolDeltaShiftSize*>(obj);
    if (volDelShift) {
        volDelShift->adjustDeltaShiftSize(collector,
                                    asset->getName(),
                                    asset->equity->spot());
    }
}

void StruckEquity::acceptValueDateCollector(const StruckEquity*  asset, 
                                            CValueDateCollector* collector)
{
    asset->equity->accept(collector);
}

void StruckEquity::acceptImntCcy(const StruckEquity* asset,
                                 AssetCcyCollector*  collector)
{
  collector->currencyValidate(asset->fx->getBaseCcyIsoCode(),
                              asset->getName());
}

void StruckEquity::acceptHoliday(const StruckEquity* asset,
                                 HolidayCollector*   collector)
{
    collector->setHoliday(asset->equity->getMarketHolidays());
}

void StruckEquity::acceptCriticalDateCollector(const StruckEquity*    asset,
                                               CriticalDateCollector* collector)
{
    asset->equity->accept((ICollector*)collector);
}

/** given a current spot level, get the next strike on the vol surface where 
    the slope is non-differentiable */
double StruckEquity::getNextStrike(const double& strike,
                                   bool          isUp,
                                   bool&         offSurface) const
{
    double volStrike;
    double adjStrike = strike;
    double fxSpot = fx->getSpot();
    if ( !Maths::isPositive(fxSpot) )
    {
        throw ModelException("StruckEquity::getNextStrike", 
                             "FX rates must be strictly positive");
    }
    adjStrike /= fxSpot;

    const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(vol.get());
    if ( nextStrike )
    {
        volStrike = nextStrike->getNextStrike(adjStrike,
                                              isUp,
                                              offSurface);
    }
    else
    {
        volStrike  = 0.0;
        offSurface = true;
    }
    // volStrike is in underlyer's ccy, but we want it in struck ccy
    return (volStrike * fxSpot);
}

/** Returns the name of the stock/asset - used to determine
    whether to shift the object */
string StruckEquity::sensName(SpotLevelProbability* shift) const {
    return equity->getName();
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
bool StruckEquity::sensShift(SpotLevelProbability* shift) {
    // we want to set the equity level based on its "raw" vol
    // NOT on the struck vol
    SimpleEquity simple(equity.get(), vol.get());
    double spot = shift->spotLevel(equity->getValueDate(), &simple);
    SpotLevel level(spot);
    equity->sensShift(&level);  // set equity spot to our new level
    return false;  // all done;
}



FXAssetConstSP StruckEquity::getFX() const {
    return fx.getSP();
}

/** returns dividend list */
DividendListConstSP  StruckEquity::getDivList() const {
    return equity->getDivList();
}

PDFCalculator* StruckEquity::pdfCalculator(const PDFRequest* request) const {
    const PDFRequestLNStrike* lnRequest =
        dynamic_cast<const PDFRequestLNStrike*>(request);
    if (lnRequest){
        return new PDFDefaultLNStrike(equity->getValueDate(), this, lnRequest);
    }
    throw ModelException("StruckEquity::pdfCalculator",
                         "Request of type "+request->getClass()->getName()+
                         " not supported");
}

/** adds the credit spread to the asset's growth curve */
void StruckEquity::makeRisky(ICreditCurveSP creditSpreads,
                             const  DateTime *maturityDate)
{
    equity->makeRiskyEquity(creditSpreads,maturityDate);
}

/* for IEqVolNamePair */
bool StruckEquity::getNamePairs(string& eqName, string& volName) const {
    eqName = getTrueName();
    volName = getVolName();
    return false;   // no more assets inside
}

//// Returns the ccy treatment for the asset. Default returns 
string StruckEquity::getCcyTreatment() const{
    return CCY_TREATMENT_STRUCK;
}

void StruckEquity::checkStruckSourceAndObsType(const ObservationSource* source,
                                               const ObservationType* type) const {
    // throw an exception if the source/type is not suitable for a struck asset
    // i.e. it should have 2 sources for asset and fx
    if(!(source->getClass() == CClass::forName("StruckObservationSource"))) {
        throw ModelException("StruckEquity::checkStruckSourceAndObsType",
                    "Inconsistent source for a struck asset. Struck assets "
                    "require 2 sources - one for the asset and one for the FX");
    }
    if(!(type->getClass() == CClass::forName("StruckObservationType"))) {
        throw ModelException("StruckEquity::checkStruckSourceAndObsType",
                    "Inconsistent observation type for a struck asset. Struck assets "
                    "require 2 types - one for the asset and one for the FX");
    }
}


// the IMarketObservable interface for retrieving a single sample
double StruckEquity::pastValue(const DateTime&             sampleDate,
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

        double equityValue = equity->pastValue(sampleDate, compoundType->getPrimaryObsType(), 
                                            source, fixType, 0, sampleRule);

        double fxValue = fx->pastValue(sampleDate, compoundType->getFXObsType(), 
                                       compoundSource->getFXSource(),
                                       fixType, 0, sampleRule);
        return equityValue * fxValue;
    } catch (ModelException& e) {
        throw ModelException(e, "StruckEquity::pastValue", 
                    "Failed to retrieve past sample for struck equity " + 
                    getTrueName());
    }
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool StruckEquity::observationDate(const DateTime&           sampleDate,
                                   const ObservationSource*  source,
                                   const SamplingConvention* sampleRule,
                                   DateTime*                 obsDate) const {
    // Currently assume the FX and the equity move together based on the equity
    return equity->observationDate(sampleDate, source, sampleRule, obsDate);
}

// the IMarketObservable interface for retrieving past samples events
double StruckEquity::addPastSampleEvent(const DateTime&     sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const {
    static const string method("StruckEquity::addPastSampleEvent");
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

        AssetFixTypeSP equityFixType = 
                AssetFixTypeSP(new AssetFixType(equity->getName()));
        double level = equity->addPastSampleEvent(sampleDate, 
                                                    compoundType->getPrimaryObsType(), 
                                                    source,
                                                    equityFixType.get(), 0, 
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
bool StruckEquity::isHoliday(const DateTime& sampleDate,
                             const ObservationSource*   source) const{
    // Currently assume the FX and the equity move together based on the equity
    return equity->isHoliday(sampleDate, source);
}

/* for reflection */
StruckEquity::StruckEquity(): CAsset(TYPE){}
    
class StruckEquityHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(StruckEquity, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(CAsset::IStruck);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(SpotLevelProbability::Shift);
        IMPLEMENTS(ICanBeRisky);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(IEqVolNamePair);
        EMPTY_SHELL_METHOD(defaultStruckEquity);
        FIELD(name, "Asset's name");
        FIELD_MAKE_OPTIONAL(name);
        FIELD(equity, "Stock");
        FIELD(vol, "Stock's Volatility");
        FIELD(fx, "What the stock is struck into");
        FIELD(corrEqFx, "correlation between equity and fx");
        FIELD_MAKE_OPTIONAL(corrEqFx);
        FIELD(getCorrFromCache, "cached flag");
        FIELD_MAKE_TRANSIENT(getCorrFromCache); // hide from dd interface
        ClassSetAcceptMethod(StruckEquity::acceptValueDateCollector);
        ClassSetAcceptMethod(StruckEquity::acceptDeltaShift);
        ClassSetAcceptMethod(StruckEquity::acceptImntCcy);
        ClassSetAcceptMethod(StruckEquity::acceptHoliday);
        ClassSetAcceptMethod(StruckEquity::acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const StruckEquity* asset,
                                        DividendCollector* collector){
        collector->addDivs(asset->equity->getDivList(), asset->fx.get());
    }

    static IObject* defaultStruckEquity(){
        return new StruckEquity();
    }
};

CClassConstSP const StruckEquity::TYPE = CClass::registerClassLoadMethod(
    "StruckEquity", typeid(StruckEquity), StruckEquityHelper::load);

DRLIB_END_NAMESPACE
