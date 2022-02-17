//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PercXCBWithVol.cpp
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
#include "edginc/PercXCBWithVol.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/** Validation */
void PercXCBWithVol::validatePop2Object(){
    static const string method("PercXCBWithVol::validatePop2Object");
    try {
        int numAssets = assets.size();
        if (numAssets < 1){
            throw ModelException(method, 
                                 "Must have at least one asset in this "
                                 "type of basket");
        }
        if (startDate.isGreater(baseDate)){
            // fwd starting - allocate memory
            spotsAtStart.resize(numAssets);
        }

        // roll value date forward by 0 days - this is to populate any samples
        // which should be set now, but are not populated yet, for example when
        // running overnight grids for instruments which have a SOD sample
        ThetaSP thetaShift = ThetaSP(new Theta(0, HolidaySP(Holiday::noHolidays())));
        sensShift(thetaShift.get());

        // validate everthing else
        AssetUtil::validatePercXCB(assets, pubWeights);

        // check all the sampling fields are okay, if they exist
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
    catch (exception& e){
        throw ModelException(e, method, "Failed for XCB "+ getName());
    }
}

/** Pull out the component assets from the market data */
void PercXCBWithVol::getMarket(const IModel* model, const MarketData* market){
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
        market->GetReferenceDate(baseDate);
        for (int i = 0; i < assets.size(); i++){
            CAsset::getAssetMarketData(model, market, ccyTreatments[i],
                                       basketYCName, assets[i]);
        }
        basketVol.getData(model, market);
        marketHols.getData(model, market);

        if (!startDate.isGreater(baseDate)){
            // only need spot at start if it's not forward starting
            bool hasSamplingInfo = sources.size() > 0;
            AssetUtil::validateXCBSpotAtStart(assets, spotsAtStart, startDate,
                                                hasSamplingInfo,
                                                obsSourceObjs, obsTypeObjs);
        }
            
        validatePop2Object();
    } 
    catch (exception& e){
        throw ModelException(e, "PercXCBWithVol::getMarket",
                             "Failed for XCB "+ getName());
    }
}

/** Returns a reference to the internal weights to use for combining assets */
const CDoubleArray& PercXCBWithVol::weightsRef() const {
    bool emptyWeights     = intWeights.empty();
    CDoubleArray& weights = intWeights;
    try{
        if (emptyWeights){
            weights = CDoubleArray(assets.size());
        }
        // always recalcuate the weights when fwd starting
        if (startDate.isGreater(baseDate)){
            /* forward starting percentage basket - always calculate
               weights (avoid problem of having to clear). Let the
               equity do the caching */
            for (int i = 0; i < assets.size(); i++){
                double fwdPrice = assets[i]->fwdValue(startDate);
                weights[i] = 100.0 * pubWeights[i]/fwdPrice;
            }
        } else if (emptyWeights){
            // weights not calculated
            for (int i = 0; i < assets.size(); i++){
                weights[i] = 100.0 * pubWeights[i]/spotsAtStart[i];
            }
        }
    } catch (exception& e){
        throw ModelException(e, "PercXCBWithVol::weightsRef", "Failed to "
                             "calculate weights for XCB "+getName());
    }
    return weights;
}

/** returns the spot price */
double PercXCBWithVol::getSpot() const{
    double spot = 0.0;
    const CDoubleArray& weights = weightsRef();
    for (int i = 0; i < assets.size(); i++){
        spot += assets[i]->getSpot() * weights[i];
    }
    return spot;
}

/** returns the asset name */
string PercXCBWithVol::getName() const{
    return name;
}

/** Helper for centralised sampling/isda adjustment for XCBs
gets a list of dates for which the given sample will finally be known for 
each component. Results are put into obsDates array and the returned bool
is false if the date is omitted (note can't omit some components and not others)*/
bool PercXCBWithVol::getObsDates(const DateTime&           sampleDate,
                      const SamplingConvention* sampleRule,
                      DateTimeArray&            obsDates,
                      DateTime&                 finalObsDate) const {
    static const string routine("PercXCBWithVol::getObsDates");
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
double PercXCBWithVol::pastValue(const DateTime&             sampleDate,
                      const ObservationType*      obsType,
                      const ObservationSource*    source,
                      const FixingType*           fixType,
                      const IObservationOverride* overrides,
                      const SamplingConvention*   sampleRule) const{
    static const string routine("PercXCBWithVol::pastValue");
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
            const CDoubleArray& weights = weightsRef();
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
bool PercXCBWithVol::observationDate(const DateTime&           sampleDate,
                          const ObservationSource*  source,
                          const SamplingConvention* sampleRule,
                          DateTime*                 obsDate) const {
    // note for an XCB the source coming in is ignored. Sampling
    // is governed by the component sources defined in the XCB
    // function calculates the date at which the sample will finally be known
    static const string routine("PercXCBWithVol::observationDate");
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
double PercXCBWithVol::addPastSampleEvent(const DateTime&          sampleDate,
                            const ObservationType*      obsType,
                            const ObservationSource*    source,
                            const FixingType*           fixType,
                            const IObservationOverride* overrides,
                            const SamplingConvention*   sampleRule,
                            PastSamplesCollector*        collector) const {
    static const string method("PercXCBWithVol:addPastSampleEvent");
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
        const CDoubleArray& weights = weightsRef();
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
bool PercXCBWithVol::isHoliday(const DateTime&            sampleDate,
                             const ObservationSource*   source) const {
    static const string method("PercXCBWithVol:isHoliday");
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
CVolProcessed * PercXCBWithVol::getProcessedVol(
    const CVolRequest* volRequest) const{
    return basketVol->getProcessedVol(volRequest, this);
}

/** Calculates the expected spot price of the asset at the given date */
double PercXCBWithVol::fwdValue(const DateTime& date) const{
    double fwd = 0.0;
    try{
        const CDoubleArray& weights = weightsRef();
        for (int i = 0; i < assets.size(); i++){
            fwd += assets[i]->fwdValue(date) * weights[i];
        }
    } catch (exception& e){
        throw ModelException(e, "PercXCBWithVol::fwdValue", "Failed for XCB "+
                             getName());
    }
    return fwd;
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void PercXCBWithVol::fwdValue(const DateTimeArray& dates,
                              CDoubleArray&        result) const{
    try{
        const CDoubleArray& weights = weightsRef();
        AssetUtil::fwdValue(assets, weights, dates, result);
    } catch (exception& e){
        throw ModelException(e, "PercXCBWithVol::fwdValue", "Failed for XCB "+
                             getName());
    }
}

/** Calculates the expected spot prices of the asset at the given dates
    respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
void PercXCBWithVol::fwdValue(const DateTimeArray&     dateList,
                              const FwdValueAlgorithm& algo,
                              CDoubleArray&            result) const{
    try{
        const CDoubleArray& weights = weightsRef();
        AssetUtil::fwdValue(assets, weights, dateList, algo, result);
    } catch (exception& e){
        throw ModelException(e, "PercXCBWithVol::fwdValue", "Failed for XCB "+
                             getName());
    }
}
    
/** Returns the name (not the ISO code) of the asset ccy */
string PercXCBWithVol::getYCName() const {
    return basketYCName;
}

// /** Calculates the expected spot price of the asset at the given date if
//     the spot price had the given value spot on spotDate */
// double PercXCBWithVol::fwdFwd(const DateTime& spotDate,
//                               double          spot, 
//                               const DateTime& fwdDate) const{
//     throw ModelException("PercXCBWithVol::fwdFwd", "Not yet done");
// }

/** Calculate the settlement date associated with a given trade date */
DateTime PercXCBWithVol::settleDate(const DateTime& tradeDate) const{
    return AssetUtil::settleDate(assets, tradeDate);
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string PercXCBWithVol::sensName(VegaSkewParallel* shift) const{
    return basketVol->getName();
}
    
bool PercXCBWithVol::sensShift(VegaSkewParallel* shift){
    // store the underlying spot price
    shift->setSpot(getSpot());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string PercXCBWithVol::sensName(VegaSkewPointwise* shift) const{
    return basketVol->getName();
}
    
ExpiryArrayConstSP PercXCBWithVol::sensExpiries(VegaSkewPointwise* shift)const{
    return ExpiryArrayConstSP(); // return null SP - the vol will do this
}

bool PercXCBWithVol::sensShift(VegaSkewPointwise* shift){
    // store the underlying spot price
    shift->setSpot(getSpot());
    return true; // continue on and do the vol
}

/** Returns the name of the XCB for (Basket-)Delta */
string PercXCBWithVol::sensName(const BasketSpot*) const{
    return name;
}

/** Shifts all components in the XCB for Delta */
TweakOutcome PercXCBWithVol::sensShift(const PropertyTweak<BasketSpot>& shift){
    try{
        double oldSpot = getSpot();
        RiskProperty<Spot>().tweakAllSubjects(IObjectSP::attachToRef(this),
                                              VoidSP(),
                                              shift.coefficient);

        return TweakOutcome(oldSpot,
                            oldSpot * (1 + shift.coefficient), // FIXME need to return init val for compat AND distnace!!
                            false); // don't shift any further

    } catch (exception& e){
        throw ModelException(e, "XCB::sensShift");
    }
}

/** Returns name identifying this object for SpotShift */
string PercXCBWithVol::sensName(SpotShift* shift) const{
    return name;
}
/** Shifts the object using given shift (see SpotShift::Shift)*/
bool PercXCBWithVol::sensShift(SpotShift* shift){
    try {
        // new shift with no market name, so shifts everything
        SpotShift newShift(shift->getShiftSize());
        for (int i = 0; i < assets.size(); i++){
            newShift.findAndShift(IObjectSP::attachToRef(&assets[i]), 
                                  OutputNameConstSP());
        }
    } 
    catch (exception& e){
        throw ModelException(e, "PercXCBWithVol::sensShift");
    }
    return false; // don't shift any further
}

/** Shifts the object using given shift (see Theta::Shift)*/
bool PercXCBWithVol::sensShift(Theta* shift) 
{
    try
    {
        DateTime rolledDate    = shift->rollDate(baseDate);
        bool     populateZeros = startDate.equals(baseDate);
       
        // was forward starting but not anymore
        if ((startDate.isGreater(baseDate) && rolledDate.isGreaterOrEqual(startDate)) ||
            populateZeros)
        {
            // THETA_FORWARD_SPOT
            if (shift->useAssetFwds())
            {
                // reset the spots to the fwd value at rolled date
                for (int i = 0; i < assets.size(); i++)
                {
                    if ( !populateZeros || Maths::isZero(spotsAtStart[i])) {
                       spotsAtStart[i] = assets[i]->fwdValue(rolledDate);
                    }
                }
            } else {
                // reset the spots to the fwd value at the start date
                for (int i = 0; i < assets.size(); i++)
                {
                    if ( !populateZeros || Maths::isZero(spotsAtStart[i])) {
                       spotsAtStart[i] = assets[i]->getSpot();
                    }
                }
            }
        }


        /* otherwise no need to do anything special since intWeights are 
           recalculated everytime they are needed, so just roll the baseDate*/
        baseDate = rolledDate; 
    }
    catch (exception& e)
    {
        throw ModelException(e, "PercXCBWithVol::sensShift (theta)"); 
    }
    return true; // tweak stuff inside
}

/** Returns the name of the XCB - used to determine whether to tweak
    the object */
string PercXCBWithVol::sensName(DeltaSurface* shift) const{
    return name;
}
    
bool PercXCBWithVol::sensShift(DeltaSurface* shift){
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
        throw ModelException(e, "PercXCBWithVol::sensShift");
    }

    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string PercXCBWithVol::sensName(VolRelativeShift* shift) const{
    return VolRelativeShift::IShift::TYPE->isInstance(basketVol.get())
        ? basketVol->getName(): "";
}
    
bool PercXCBWithVol::sensShift(VolRelativeShift* shift){
    // store the underlying spot price
    shift->setSpot(baseDate, this);
    return true; // continue on and do the vol
}


/** returns sensitive strikes for a given vol request */
void PercXCBWithVol::getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const
{
    // to do: alter getSensitiveStrikes signature
    const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
        *(IObject*)(volRequest));
    // get the sensitive strike for the basket if required
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
void PercXCBWithVol::recordFwdAtMat(
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

void PercXCBWithVol::acceptCollector(const PercXCBWithVol* asset, 
                                     StartDateCollector* collector)
{
    collector->startDateValidate(asset->startDate, 
                                 asset->getName(),
                                 false);
    // check start dates for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void PercXCBWithVol::acceptNameCollector(const PercXCBWithVol* asset, 
                                         AssetNameCollector* collector)
{
    collector->assetNameValidate(asset->getName());
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void PercXCBWithVol::acceptFutureCollector(const PercXCBWithVol* asset, 
                                           FutureExpiryCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void PercXCBWithVol::acceptValueDateCollector(const PercXCBWithVol* asset, 
                                              CValueDateCollector* collector)
{
    collector->valueDateValidate(asset->baseDate, asset->getName());
}

void PercXCBWithVol::acceptDeltaShift(const PercXCBWithVol* asset, 
                                      ShiftSizeCollector* collector)
{
    for (int idx=0 ; idx<asset->assets.size() ; idx++)
    {
        asset->assets[idx]->accept(collector);
    }
}

void PercXCBWithVol::acceptImntCcy(const PercXCBWithVol* asset,
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
        throw ModelException(e, "PercXCBWithVol::acceptImntCcy",
                             "failed for XCB " + asset->getName());
    }
}

void PercXCBWithVol::acceptHoliday(const PercXCBWithVol* asset,
                                   HolidayCollector* collector)
{
    collector->setHoliday(asset->marketHols.getSP());
}

void PercXCBWithVol::acceptWrapperNameCollector(const PercXCBWithVol* asset, 
                                                WrapperNameCollector* collector)
{
    collector->addName(asset->getName());
}

void PercXCBWithVol::acceptCriticalDateCollector(const PercXCBWithVol* asset, 
                                                 CriticalDateCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

 
/** Returns fair value of fund price */
double PercXCBWithVol::fairValue() const {
    return getSpot(); // do we want settlement ? need a yield curve if we do
}
  
PDFCalculator* PercXCBWithVol::pdfCalculator(const PDFRequest* request) const {
    if (IPDFCalculator::TYPE->isInstance(basketVol.get())) {
        const IPDFCalculator* pdf = dynamic_cast< const IPDFCalculator*>(basketVol.get());
        return pdf->getPDFCalculator(request, this);
    }
    throw ModelException("PercXCBWithVol::pdfCalculator",
                         "vol of type (" + basketVol->getClass()->getName() +
                         ") has no pdf calculator");
}

/* for IEqVolNamePair */
bool PercXCBWithVol::getNamePairs(string& eqName, string& volName) const {
    eqName  = name;
    volName = basketVol->getName();
    return false;  // no point in reporting our components
}

/** Can this asset physically settle? */
bool PercXCBWithVol::canPhysicallySettle() const {
    return false;
}

/* for reflection */
PercXCBWithVol::PercXCBWithVol(): CAsset(TYPE){}
    
class PercXCBWithVolHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PercXCBWithVol, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(ITweakableWithRespectTo<BasketSpot>);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(SpotShift::Shift);
        //IMPLEMENTS(IEqVolNamePair);
        EMPTY_SHELL_METHOD(defaultPercXCBWithVol);
        FIELD(name, "XCB's name");
        FIELD(basketVol, "XCB's Volatility");
        FIELD(assets, "Assets in XCB");
        FIELD(ccyTreatments, "Currency treatments");
        FIELD(basketYCName, "XCB's currency");
        FIELD(pubWeights, "Percentage of each asset");
        FIELD(marketHols, "Market hols for XCB");
        FIELD(baseDate, "today");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(startDate, "When XCB starts");
        FIELD(spotsAtStart, "Spots of asset at startDate");
        FIELD_MAKE_OPTIONAL(spotsAtStart);
        FIELD(basketCcyCode, "XCB's currency code");
        FIELD_MAKE_TRANSIENT(basketCcyCode);
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
        ClassSetAcceptMethod(PercXCBWithVol::acceptCollector);
        ClassSetAcceptMethod(PercXCBWithVol::acceptNameCollector);
        ClassSetAcceptMethod(PercXCBWithVol::acceptFutureCollector);
        ClassSetAcceptMethod(PercXCBWithVol::acceptValueDateCollector);
        ClassSetAcceptMethod(PercXCBWithVol::acceptDeltaShift);
        ClassSetAcceptMethod(PercXCBWithVol::acceptImntCcy);
        ClassSetAcceptMethod(PercXCBWithVol::acceptHoliday);
        ClassSetAcceptMethod(PercXCBWithVol::acceptWrapperNameCollector);
        ClassSetAcceptMethod(PercXCBWithVol::acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const PercXCBWithVol*    asset,
                                        DividendCollector* collector){
        const CDoubleArray& weights = asset->weightsRef();
        for (int idx=0 ; idx<asset->assets.size() ; idx++) {
            collector->processComponentAsset(asset->assets[idx].get(),
                                             weights[idx]);
        }
    }

    static IObject* defaultPercXCBWithVol(){
        return new PercXCBWithVol();
    }
};

CClassConstSP const PercXCBWithVol::TYPE = CClass::registerClassLoadMethod(
    "PercXCBWithVol", typeid(PercXCBWithVol), PercXCBWithVolHelper::load);

DRLIB_END_NAMESPACE
