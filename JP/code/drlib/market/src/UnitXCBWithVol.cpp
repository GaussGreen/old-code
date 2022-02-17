//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : UnitXCBWithVol.cpp
//
//   Description : Implement asset for unit weighted cross currency basket
//                 with vol override
//
//   Author      : Mark A Robson
//
//   Date        : 9 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UnitXCBWithVol.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/** Validation */
void UnitXCBWithVol::validatePop2Object(){
    static const string method("UnitXCBWithVol::validatePop2Object");
    try {
        AssetUtil::validateUnitXCB(assets, weights);

        // check all the sampling fields are okay, if they exist
        int numAssets = assets.size();
        if (sources.size() > 0) {
            if (sources.size() != numAssets) {
                throw ModelException(method,
                            "Must have as many sources (" + 
                            Format::toString(sources.size()) +
                            ") as component assets (" +
                            Format::toString(numAssets) +
                            "in this type of basket");
            }
            if (obsTypes.size() != numAssets) {
                throw ModelException(method,
                            "Must have as many observation types (" + 
                            Format::toString(obsTypes.size()) +
                            ") as component assets (" +
                            Format::toString(numAssets) +
                            "in this type of basket");
            }
            if (fxSources.size() != numAssets) {
                throw ModelException(method,
                            "Must have as many FX sources (" + 
                            Format::toString(fxSources.size()) +
                            ") as component assets (" +
                            Format::toString(numAssets) +
                            "in this type of basket");
            }
            if (fxObsTypes.size() != numAssets) {
                throw ModelException(method,
                            "Must have as many FX observation types (" + 
                            Format::toString(fxObsTypes.size()) +
                            ") as component assets (" +
                            Format::toString(numAssets) +
                            "in this type of basket");
            }
            // then build the final objects we're going to need for each component
            obsSourceObjs = ObservationSourceArray(assets.size());
            obsTypeObjs = ObservationTypeArray(assets.size());
            for (int i = 0; i < numAssets; i++) {
                if(CString::equalsIgnoreCase(ccyTreatments[i], CAsset::CCY_TREATMENT_STRUCK)) {
                    obsSourceObjs[i] = 
                        ObservationSourceSP(new StruckObservationSource(sources[i],
                                                                        fxSources[i]));
                    obsTypeObjs[i] =
                        ObservationTypeSP(new StruckObservationType(obsTypes[i].get(),
                                                                    fxObsTypes[i].get()));
                } else {
                    obsSourceObjs[i] = 
                        ObservationSourceSP(new ObservationSource(sources[i]));
                    obsTypeObjs[i] = obsTypes[i];
                }
            }
        }
    } 
    catch(exception& e){
        throw ModelException(e, method, "Failed for xcb "+getName());
    }
}

/** Pull out the component assets from the market data */
void UnitXCBWithVol::getMarket(const IModel* model, const MarketData* market){
    try {
        // need the ccy code of payoff ccy for validation
        if (basketCcyCode.empty()){ // in case we're called twice
            if (!assets[0]) {
                // if using market cache
                YieldCurveWrapper yc(basketYCName);
                yc.getData(model, market);
                basketCcyCode = yc->getCcy();
            }
            else {
                basketCcyCode = basketYCName;
            }
        }
        // get the data we need from the cache
        for (int i = 0; i < assets.size(); i++){
            CAsset::getAssetMarketData(model, market, ccyTreatments[i],
                                       basketYCName, assets[i]);
        }
        basketVol.getData(model, market);
        marketHols.getData(model, market);
        market->GetReferenceDate(baseDate);
    } catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::getMarket", 
                             "Failed for xcb "+getName());
    }
}

/** returns the spot price */
double UnitXCBWithVol::getSpot() const{
    double spot = 0.0;
    for (int i = 0; i < assets.size(); i++){
        spot += assets[i]->getSpot() * weights[i];
    }
    return spot;
}

/** returns the asset name */
string UnitXCBWithVol::getName() const{
    return name;
}

/** Helper for centralised sampling/isda adjustment for XCBs
gets a list of dates for which the given sample will finally be known for 
each component. Results are put into obsDates array and the returned bool
is false if the date is omitted (note can't omit some components and not others)*/
bool UnitXCBWithVol::getObsDates(const DateTime&           sampleDate,
                      const SamplingConvention* sampleRule,
                      DateTimeArray&            obsDates,
                      DateTime&                 finalObsDate) const {
    static const string routine("UnitXCBWithVol::getObsDates");
    bool sampled;
    try {
        // ISS - note this logic is fine if all assets roll in same
        // direction. Might need revisiting if we have asset level
        // overrides i.e. some roll fwd, some back
        // we shouldn't have components rolling individually with omit rule
        sampled = assets[0]->observationDate(sampleDate, 
                                             obsSourceObjs[0].get(),
                                             sampleRule,
                                             &(obsDates[0]));
        finalObsDate = obsDates[0];
        bool lastsampled;
        bool rollAll = sampleRule->rollAssetsTogether();
        bool omit = (rollAll && !sampled);
        for (int i = 1; i < assets.size() && !omit; i++) {
            lastsampled = sampled;
            // note for an XCB the source coming in is ignored. Sampling
            // is governed by the component sources defined in the XCB
            // function calculates the date at which the sample will finally be known
            sampled = assets[i]->observationDate(rollAll ? obsDates[i-1] : sampleDate, 
                                                 obsSourceObjs[i].get(),
                                                 sampleRule,
                                                 &(obsDates[i]));
            omit = (rollAll && !sampled);
            if (!rollAll && (lastsampled != sampled)) {
                throw ModelException("Cannot use the omit adjustment convention on an "
                                     "XCB with the components rolling separately");
            }
            // note we're after date at which all bits will finally be known
            if (obsDates[i] > finalObsDate) {
                finalObsDate = obsDates[i];
            }
        }
        // if assets roll together we need to go back and adjust things
        if (rollAll && !omit) {
            // all assets move to the last date
            for (int i = 0; i < assets.size()-1; i++) {
                obsDates[i] = obsDates[assets.size()-1];
            }
        }
    } catch (exception& e) {
        throw ModelException(e, routine, "Failed to retrieve observation date "
                            "for sample date (" + sampleDate.toString() + ") in "
                            "XCB " + name);
    }
    return sampled;
}

// the IMarketObservable interface for retrieving a single sample
double UnitXCBWithVol::pastValue(const DateTime&             sampleDate,
                      const ObservationType*      obsType,
                      const ObservationSource*    source,
                      const FixingType*           fixType,
                      const IObservationOverride* overrides,
                      const SamplingConvention*   sampleRule) const{
    static const string routine("UnitXCBWithVol::pastValue");
    try {
        // first try and get the sample from the overrides
        // overrides will be for the XCB levels
        if (overrides) {
            double pValue = 0.0;
            if(overrides->value(sampleDate, obsType, &pValue)) {
                return pValue;
            } 
        }                                   
        // note for an XCB the source and obsType coming in are ignored. Sampling
        // is governed by the component sources/obsTypes defined in the XCB
        DateTimeArray obsDates(assets.size());
        DateTime finalObsDate;
        // shouldn't be trying to sample on an omitted date
        if (getObsDates(sampleDate, sampleRule, obsDates, finalObsDate)) {
            double value = 0.0;
            for (int i = 0; i < assets.size(); i++){
                try {
                    value += assets[i]->pastValue(obsDates[i], obsTypeObjs[i].get(),
                                                obsSourceObjs[i].get(), fixType,
                                                0, sampleRule) 
                            * weights[i];
                } catch (exception& e) {
                    throw ModelException(e, 
                        "Failed to retrieve sample for component " +
                        assets[i]->getTrueName() + " for date (" +
                        obsDates[i].toString() + ")in XCB " + name);
                }
            }
            return value;
        } else {
            throw ModelException(routine, "Cannot sample on an omitted date for an XCB");
        }
    } catch (exception& e) {
        throw ModelException(e, routine,  "Failed to retrieve sample for XCB " +
                                          name + " for sample date (" +
                                          sampleDate.toString() + ")");
    }
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool UnitXCBWithVol::observationDate(const DateTime&           sampleDate,
                          const ObservationSource*  source,
                          const SamplingConvention* sampleRule,
                          DateTime*                 obsDate) const {
    // note for an XCB the source coming in is ignored. Sampling
    // is governed by the component sources defined in the XCB
    // function calculates the date at which the sample will finally be known
    static const string routine("UnitXCBWithVol::observationDate");
    try {
        DateTimeArray obsDates(assets.size());
        bool sampled = getObsDates(sampleDate, sampleRule, obsDates, *obsDate);
        return sampled;
    } catch (exception& e) {
        throw ModelException(e, routine, "Failed to retrieve observation date "
                            "for sample date (" + sampleDate.toString() + ") in "
                            "XCB " + name);
    }
}

// the IMarketObservable interface for retrieving past samples events
double UnitXCBWithVol::addPastSampleEvent(const DateTime&          sampleDate,
                            const ObservationType*      obsType,
                            const ObservationSource*    source,
                            const FixingType*           fixType,
                            const IObservationOverride* overrides,
                            const SamplingConvention*   sampleRule,
                            PastSamplesCollector*        collector) const {
    static const string method("UnitXCBWithVol:addPastSampleEvent");
    try {
        double level = 0.0;
        // first try and get the sample from the overrides
        if (overrides) {
            if (overrides->addPastSampleEvent(sampleDate, obsType, source,
                                            fixType, collector, level)) {
                return level;
            }
        }                                   
        DateTimeArray obsDates(assets.size());
        DateTime finalObsDate;
        bool sampled = getObsDates(sampleDate, sampleRule, 
                                   obsDates, finalObsDate);
        if (!sampled) {
            throw ModelException(method, "Internal error "
                                "Should not attempt to get sampling "
                                "event for omitted date");
        }
        for (int i = 0; i < assets.size(); i++){
            AssetFixTypeSP assetFixType = 
                    AssetFixTypeSP(new AssetFixType(assets[i]->getName()));
            level += assets[i]->addPastSampleEvent(obsDates[i], obsTypeObjs[i].get(),
                                        obsSourceObjs[i].get(), assetFixType.get(),
                                        0, sampleRule, collector) 
                    * weights[i];
        }
        // add in the top level sample event
        collector->addSample(sampleDate, finalObsDate, source, obsType,
                            level, fixType, false);
        return level;

    } catch (exception& e) {
        throw ModelException(e, method, "Failed to retrieve past sample event "
                                        "for XCB " + name + " for sample date (" +
                                         sampleDate.toString() + ")");
    }
}

// the IMarketObservable interface for 
// is the given date a holiday for the relevant source
bool UnitXCBWithVol::isHoliday(const DateTime&            sampleDate,
                             const ObservationSource*   source) const {
    static const string method("UnitXCBWithVol:isHoliday");
    try {
        for (int i = 0; i < assets.size(); i++){
            if (assets[i]->isHoliday(sampleDate, obsSourceObjs[i].get())) {
                return true;
            }
        }
        return false;
    } catch (exception& e) {
        throw ModelException(e, method, "Holiday lookup failed for XCB "
                                        + name + " and date (" +
                                        sampleDate.toString() + ")");
    }
}

/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
CVolProcessed * UnitXCBWithVol::getProcessedVol(
    const CVolRequest* volRequest) const{
    return basketVol->getProcessedVol(volRequest, this);
}

/** Calculates the expected spot price of the asset at the given date */
double UnitXCBWithVol::fwdValue(const DateTime& date) const{
    double fwd = 0.0;
    try{
        for (int i = 0; i < assets.size(); i++){
            fwd += assets[i]->fwdValue(date) * weights[i];
        }
    } catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::fwdValue", "Failed for XCB "+
                             getName());
    }
    return fwd;
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void UnitXCBWithVol::fwdValue(const DateTimeArray& dates,
                              CDoubleArray&        result) const{
    try{
        AssetUtil::fwdValue(assets, weights, dates, result);
    } catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::fwdValue", "Failed for XCB "+
                             getName());
    }
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void UnitXCBWithVol::fwdValue(const DateTimeArray&     dates,
                              const FwdValueAlgorithm& algo,
                              CDoubleArray&            result) const{
    try{
        AssetUtil::fwdValue(assets, weights, dates, algo, result);
    } catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::fwdValue", "Failed for XCB "+
                             getName());
    }
}

/** Returns the name (not the ISO code) of the asset ccy */
string UnitXCBWithVol::getYCName() const {
    return basketYCName;
}

// /** Calculates the expected spot price of the asset at the given date if
//     the spot price had the given value spot on spotDate */
// double UnitXCBWithVol::fwdFwd(const DateTime& spotDate,
//                               double          spot, 
//                               const DateTime& fwdDate) const{
//     throw ModelException("UnitXCBWithVol::fwdFwd", "Not yet done");
//     return 0;
// }

/** Calculate the settlement date associated with a given trade date */
DateTime UnitXCBWithVol::settleDate(const DateTime& tradeDate) const{
    return AssetUtil::settleDate(assets, tradeDate);
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string UnitXCBWithVol::sensName(VegaSkewParallel* shift) const{
    return basketVol->getName();
}
    
bool UnitXCBWithVol::sensShift(VegaSkewParallel* shift){
    // store the underlying spot price
    shift->setSpot(getSpot());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string UnitXCBWithVol::sensName(VegaSkewPointwise* shift) const{
    return basketVol->getName();
}
    
ExpiryArrayConstSP UnitXCBWithVol::sensExpiries(VegaSkewPointwise* shift)const{
    return ExpiryArrayConstSP(); // return null SP - the vol will do this
}

bool UnitXCBWithVol::sensShift(VegaSkewPointwise* shift){
    // store the underlying spot price
    shift->setSpot(getSpot());
    return true; // continue on and do the vol
}

/** Returns the name of the XCB for Delta */
string UnitXCBWithVol::sensName(const BasketSpot*) const{
    return name;
}

/** Shifts all components in the XCB for Delta */
TweakOutcome UnitXCBWithVol::sensShift(const PropertyTweak<BasketSpot>& shift){
    try{
        double oldSpot = getSpot();
        RiskProperty<Spot>().tweakAllSubjects(IObjectSP::attachToRef(this),
                                               VoidSP(),
                                               shift.coefficient);

        return TweakOutcome(oldSpot,
                            oldSpot * (1 + shift.coefficient), // FIXME need to return init val for compat AND distnace!!
                            false); // don't shift any further
    }
    catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::sensShift");
    }
}

/** Returns name identifying this object for SpotShift */
string UnitXCBWithVol::sensName(SpotShift* shift) const{
    return name;
}
/** Shifts the object using given shift (see SpotShift::Shift)*/
bool UnitXCBWithVol::sensShift(SpotShift* shift){
    try {
        // new shift with no market name, so shifts everything
        SpotShift newShift(shift->getShiftSize());
        for (int i = 0; i < assets.size(); i++){
            newShift.findAndShift(IObjectSP::attachToRef(&assets[i]),
                                  OutputNameConstSP());
        }
    } 
    catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::sensShift");
    }
    return false; // don't shift any further
}


/** Returns the name of the XCB - used to determine whether to tweak
    the object */
string UnitXCBWithVol::sensName(DeltaSurface* shift) const{
    return name;
}
    
bool UnitXCBWithVol::sensShift(DeltaSurface* shift){
    // have to do this ourselves a la BasketDelta
    try {
        // must set the initial value
        shift->setInitialValue(getSpot());

        // store the underlying spot price
        shift->setSpot(getSpot(), name, basketVol->getName());

        // now set ourselves up for basket style shift
        // may as well just move component spot (not vol)
        // as it make sno contribution here
        DeltaSP newShift(new Delta(shift->getShiftSize()));
        for (int i = 0; i < assets.size(); i++){
            newShift->findAndShift(IObjectSP::attachToRef(&assets[i]),
                                   OutputNameConstSP());
        }
    } 
    catch (exception& e){
        throw ModelException(e, "UnitXCBWithVol::sensShift");
    }

    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string UnitXCBWithVol::sensName(VolRelativeShift* shift) const{
    return VolRelativeShift::IShift::TYPE->isInstance(basketVol.get())
        ? basketVol->getName(): "";
}
    
bool UnitXCBWithVol::sensShift(VolRelativeShift* shift){
    // store the underlying spot price
    shift->setSpot(baseDate, this);
    return true; // continue on and do the vol
}

/** returns sensitive strikes for a givyen vol request */
void UnitXCBWithVol::getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const
{
    // to do: alter getSensitiveStrikes signature
    const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
        *(IObject*)(volRequest));
    if ( sensStrikeDesc.forwardOnly == false && 
         outputName->equals(basketVol->getName()) ) {
        // pull out the basket strike
        request.getSensitiveStrike(getSpot(), sensitiveStrikes);
    }

    /** for the components we only want the strikes which are relevant 
        for the forward price calculation */
    SensitiveStrikeDescriptor cmpStrikeDesc;
    cmpStrikeDesc.forwardOnly = true;

    // from the assets we only need 
    for (int i=0 ; i<assets.size() ; ++i) {
        assets[i]->getSensitiveStrikes(volRequest, 
                                          outputName, 
                                          cmpStrikeDesc, 
                                          sensitiveStrikes);
    }
}

/** record forwards at maturity*/
void UnitXCBWithVol::recordFwdAtMat(
                          OutputRequest*  request,
                          CResults*       results,
                          const DateTime& maturityDate) const
{
    double fwd = fwdValue(maturityDate);

    // wrap forward price into object
    CDoubleSP fwdObj(CDouble::create(fwd));

    results->storeRequestResult(request,
                                fwdObj,
                                OutputNameSP(new OutputName(getName())));

    // record forward at maturity for each asset
    for (int i = 0; i < assets.size() ; ++i)
    {
        assets[i]->recordFwdAtMat(request, results, maturityDate);
    }
}

void UnitXCBWithVol::acceptCollector(const UnitXCBWithVol* asset, 
                                     StartDateCollector* collector)
{
    // check start dates for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void UnitXCBWithVol::acceptNameCollector(const UnitXCBWithVol* asset, 
                                         AssetNameCollector* collector)
{
    collector->assetNameValidate(asset->getName());
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void UnitXCBWithVol::acceptFutureCollector(const UnitXCBWithVol* asset, 
                                           FutureExpiryCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void UnitXCBWithVol::acceptValueDateCollector(const UnitXCBWithVol* asset, 
                                              CValueDateCollector* collector)
{
    // check value date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void UnitXCBWithVol::acceptDeltaShift(const UnitXCBWithVol* asset, 
                                      ShiftSizeCollector* collector)
{
    for (int idx=0 ; idx<asset->assets.size() ; idx++)
    {
        asset->assets[idx]->accept(collector);
    }
}

void UnitXCBWithVol::acceptImntCcy(const UnitXCBWithVol* asset,
                                   AssetCcyCollector* collector)
{
    collector->currencyValidate(asset->basketCcyCode, asset->getName());
    
    try {
        for (int i = 0; i < asset->assets.size(); ++i) {
            IObjectConstSP component(
                IObjectConstSP::attachToRef(asset->assets[i].get()));
            
            AssetCcyCollector::validateAllCurrencies(component,
                                                     asset->basketCcyCode);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "UnitXCBWithVol::acceptImntCcy",
                             "failed for XCB " + asset->getName());
    }
}

void UnitXCBWithVol::acceptHoliday(const UnitXCBWithVol* asset,
                                   HolidayCollector* collector)
{
    collector->setHoliday(asset->marketHols.getSP());
}

void UnitXCBWithVol::acceptWrapperNameCollector(const UnitXCBWithVol* asset, 
                                                WrapperNameCollector* collector)
{
    collector->addName(asset->getName());
}

void UnitXCBWithVol::acceptCriticalDateCollector(const UnitXCBWithVol* asset, 
                                                 CriticalDateCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

  
/** Returns fair value of fund price */
double UnitXCBWithVol::fairValue() const {
    return getSpot(); // do we want settlement ? need a yield curve if we do
}

PDFCalculator* UnitXCBWithVol::pdfCalculator(const PDFRequest* request) const {
    if (IPDFCalculator::TYPE->isInstance(basketVol.get())) {
        const IPDFCalculator* pdf = dynamic_cast< const IPDFCalculator*>(basketVol.get());
        return pdf->getPDFCalculator(request, this);
    }
    throw ModelException("UnitXCBWithVol::pdfCalculator",
                         "vol of type (" + basketVol->getClass()->getName() +
                         ") has no pdf calculator");
}

/* for IEqVolNamePair */
bool UnitXCBWithVol::getNamePairs(string& eqName, string& volName) const {
    eqName  = name;
    volName = basketVol->getName();
    return false;  // no point in reporting our components
}
 
/* for reflection */
UnitXCBWithVol::UnitXCBWithVol(): CAsset(TYPE){}
    
class UnitXCBWithVolHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(UnitXCBWithVol, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(ITweakableWithRespectTo<BasketSpot>);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(SpotShift::Shift);
        //IMPLEMENTS(IEqVolNamePair);
        EMPTY_SHELL_METHOD(defaultUnitXCBWithVol);
        FIELD(name, "XCB's name");
        FIELD(basketVol, "XCB's Volatility");
        FIELD(assets, "Assets in XCB");
        FIELD(ccyTreatments, "Currency treatments");
        FIELD(basketYCName, "XCB's currency");
        FIELD(weights, "Units of each asset");
        FIELD(marketHols, "Market hols for XCB");
        FIELD(basketCcyCode, "XCB's currency code");
        FIELD_MAKE_TRANSIENT(basketCcyCode);
        FIELD_NO_DESC(baseDate);
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(sources, "The source for sampling of each component");
        FIELD_MAKE_OPTIONAL(sources);
        FIELD(obsTypes, "The observation types for sampling of each component");
        FIELD_MAKE_OPTIONAL(obsTypes);
        FIELD(fxSources, "The FX source for sampling of each component");
        FIELD_MAKE_OPTIONAL(fxSources);
        FIELD(fxObsTypes, "The FX observation types for sampling of each component");
        FIELD_MAKE_OPTIONAL(fxObsTypes);
        FIELD(obsSourceObjs, "The (possibly composite) sources for sampling of each component");
        FIELD_MAKE_TRANSIENT(obsSourceObjs);
        FIELD(obsTypeObjs, "The (possibly composite) observation types for sampling of each component");
        FIELD_MAKE_TRANSIENT(obsTypeObjs);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptCollector);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptNameCollector);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptFutureCollector);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptValueDateCollector);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptDeltaShift);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptImntCcy);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptHoliday);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptWrapperNameCollector);
        ClassSetAcceptMethod(UnitXCBWithVol::acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const UnitXCBWithVol* asset,
                                        DividendCollector*    collector){
        for (int idx=0 ; idx<asset->assets.size() ; idx++) {
            collector->processComponentAsset(asset->assets[idx].get(),
                                             asset->weights[idx]);
        }
    }
    static IObject* defaultUnitXCBWithVol(){
        return new UnitXCBWithVol();
    }
};

CClassConstSP const UnitXCBWithVol::TYPE = CClass::registerClassLoadMethod(
    "UnitXCBWithVol", typeid(UnitXCBWithVol), UnitXCBWithVolHelper::load);

DRLIB_END_NAMESPACE
