//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XCB.cpp
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
#include "edginc/XCB.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/CompositeVol.hpp"
#include "edginc/Addin.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/PDFParamLNStrike.hpp"
#include "edginc/VolProcessedDispatch.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/CorrelationCategory.hpp"

DRLIB_BEGIN_NAMESPACE

const string XCB::NO_SMILE = "C";
const string XCB::FIXED_STRIKE_SMILE = "D";
const string XCB::FLOAT_STRIKE_SMILE = "E";

const int XCB::noSmile = 0;
const int XCB::fixedStrike = 1;
const int XCB::floatStrike = 2;

#define EDG_MINIMUM_DELTA_SHIFT 0.0001

XCB::~XCB(){}

/** Pull out the component assets & correlations from the market data */
void XCB::getMarket(const IModel* model, const MarketData* market){
    static const string method("XCB::getMarket");
    try{
        int numAssets = assets.size();
        /** correlations: if we are explicitly supplied with just the numbers,
            then corr obj were built in validatePop2Object */ 

        if (!correlations.empty()){
            int pos = 0;
            for (int i = 0; i < assets.size(); i++) {
                for (int j = i + 1; j < assets.size(); j++, pos++) {
#if 0
                    corrObjects[pos]->configureForSensitivities(
                        assets[i]->getClass(), assets[j]->getClass());
#else
                    // force it to be EQ-EQ
                    corrObjects[pos]->configureForSensitivities(
                        CAsset::TYPE, CAsset::TYPE);
#endif
                    corrObjects[pos]->getMarket(model, market);
                }
            }
        } else {
            int pos = 0;
            corrObjects.resize((numAssets * numAssets - numAssets)/2);
            for (int i = 0; i < assets.size(); i++) {
                for (int j = i + 1; j < assets.size(); j++, pos++) {
                    string corrName = 
                        market->getCorrelationName(assets[i].getName(),
                                                       assets[j].getName());
                    CorrelationBaseSP corr(model->getCorrelation(corrName,
#if 0
                        assets[i]->getClass(),assets[j]->getClass(),
#else
                        CAsset::TYPE, CAsset::TYPE, // force it to be EQ-EQ
#endif
                        Correlation::TYPE,market));
                     corrObjects[pos] = CorrelationSP::dynamicCast(corr);
                }
            }
        }

        if(useCorrelationTerm) {
        /** now: all corr category and corr term objects should be available.
            if not, then we wont fail, but assign a dummy CorrelationTerm */
            
            // allocate space
            corrTermArray.resize((numAssets * numAssets - numAssets)/2);

            int pos=0, i,j;
            for (i = 0; i < assets.size(); i++) {
                if (market->hasData(assets[i].getName(), CorrelationCategory::TYPE)) {
                    MarketObjectSP mo1 = market->GetData(assets[i].getName(), CorrelationCategory::TYPE);
                    CorrelationCategorySP category1 = CorrelationCategorySP::dynamicCast(mo1);
                    for (j = i+1; j < assets.size(); j++, pos++) {
                        if (market->hasData(assets[j].getName(), CorrelationCategory::TYPE)) {
                            MarketObjectSP mo2 = market->GetData(assets[j].getName(), CorrelationCategory::TYPE);
                            CorrelationCategorySP category2 = CorrelationCategorySP::dynamicCast(mo2);

                            const string& corrTermName = 
                                market->getCorrelationTermName(category1->getCategoryName(), category2->getCategoryName());
                            MarketObjectSP moTerm = market->GetData(corrTermName, CorrelationTerm::TYPE);
                            corrTermArray[pos] = CorrelationTermSP::dynamicCast(moTerm);
                            corrTermArray[pos]->getMarket(model, market);
                        } else {
                            // fill with dummy data
                            corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                        }
                    }
                } else {
                    for (j = i+1; j < assets.size(); j++, pos++) {
                        // fill with dummy data
                        corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                    }

                }
            }
        }   

        // set correlation skew parameters
        getCorrelationSkewParameters(model,market);

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
        // get the rest of the data we need from the cache
        market->GetReferenceDate(baseDate);
        for (int i = 0; i < assets.size(); i++){
            CAsset::getAssetMarketData(model, market, ccyTreatments[i],
                                       basketYCName, assets[i]);
        }

        marketHols.getData(model, market);
        timeMetric->getMarket(model, market); // time metric needs hols

        bool needInitialSpots = this->needInitialSpots();

        // roll value date forward by 0 days - this is to populate any samples
        // which should be set now, but are not populated yet, for example when
        // running overnight grids for instruments which have a SOD sample
        ThetaSP thetaShift = ThetaSP(new Theta(0, HolidaySP(Holiday::noHolidays())));
        sensShift(thetaShift.get());

        // Theta shift deals with future spot, this deals with spots retrieved
        // from the asset history
        if (needInitialSpots && startDate<baseDate) {
            // Attempt to populate missing spots from asset history
            bool hasSamplingInfo = sources.size() > 0;
            AssetUtil::validateXCBSpotAtStart(assets, spotsAtStart, startDate,
                                                hasSamplingInfo,
                                                obsSourceObjs, obsTypeObjs);
        }

        // validate everthing else
        if (unitWeights){
            AssetUtil::validateUnitXCB(assets, pubWeights);
        } else {
            AssetUtil::validatePercXCB(assets, pubWeights);
        }
            
        // calculate derived basketAtStart value
        if (needInitialSpots){
            if (!unitWeights){
                basketAtStart = 100.0; // by definition
            } else {
                basketAtStart = 0.0;
                for (int i = 0; i < assets.size(); i++){
                    basketAtStart += pubWeights[i] * spotsAtStart[i];
                }
            }
        }

    }
    catch (exception& e){
        throw ModelException(e, "XCB::getMarket", "Failed for XCB "+getName());
    }
}

/** Validation */
void XCB::validatePop2Object(){
    static const string method("XCB::validatePop2Object");
    try {
        int numAssets = assets.size();
        
        if (numAssets < 2){
            throw ModelException(method, "Must have at least two assets in this type of XCB");
        }
        // if the optional spotsAtStart wasn't specified, it needs to be
        // created at the correct size
        if (startDate.isGreater(baseDate) || spotsAtStart.size() == 0){
            spotsAtStart.resize(numAssets); // allocate memory
        }
        /** correlations: if we are explicitly supplied with just the numbers,
            then generate correlation objects */
        if (!correlations.empty()){
            // validation of correlations
            if (correlations.numRows() != numAssets ||
                correlations.numCols() != numAssets){
                throw ModelException(method,
                    "Correlations must be an n x n matrix "
                    "with n = number of assets in XCB");
            }
            correlations.checkSymmetric();
            // First allocate space
            corrObjects.resize((numAssets * numAssets - numAssets)/2);
            // then generate correlation objects
            int pos = 0, i,j;
            for (i = 0; i < assets.size(); i++) {
                for (j = i + 1; j < assets.size(); j++, pos++) {
                    // the Correlation constructor reorders the names into
                    // alphabetical order if no name supplied. For backwards
                    // compatibility we keep the original order
                    string corrName(assets[i].getName() + "_" + 
                                    assets[j].getName());
                    corrObjects[pos] = 
                        CorrelationSP(new Correlation(corrName,
                                                      assets[i].getName(),
                                                      assets[j].getName(),
                                                      correlations[i][j]));
                }
            }
        }

        // check smile
        if (smileType == NO_SMILE){
            smile = noSmile;
        } else if (smileType == FIXED_STRIKE_SMILE){
            smile = fixedStrike;
        } else if (smileType == FLOAT_STRIKE_SMILE){
            smile = floatStrike;
        } else {
            throw ModelException(method, "Unknown smile type "+smileType);
        }

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
    } catch (exception& e){
        throw ModelException(e, method, "Failed for XCB " + getName());
    }
}

/** Returns a reference to the internal weights to use for combining assets */
const CDoubleArray& XCB::weightsRef() const {
    if (unitWeights){
        return pubWeights;
    } else {
        bool emptyWeights  = intWeights.empty();
        CDoubleArray&  weights = intWeights;
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
            return weights;
        } catch (exception& e){
            throw ModelException(e, "XCB::weightsRef", "Failed to "
                                 "calculate weights for XCB "+getName());
        }
    }
}

/** returns the spot price */
double XCB::getSpot() const{
    double spot = 0.0;
    const CDoubleArray& weights = weightsRef();
    for (int i = 0; i < assets.size(); i++){
        spot += assets[i]->getSpot() * weights[i];
    }
    return spot;
}

/** returns the asset name */
string XCB::getName() const{
    return name;
}

/** Helper for centralised sampling/isda adjustment for XCBs
gets a list of dates for which the given sample will finally be known for 
each component. Results are put into obsDates array and the returned bool
is false if the date is omitted (note can't omit some components and not others)*/
bool XCB::getObsDates(const DateTime&           sampleDate,
                      const SamplingConvention* sampleRule,
                      DateTimeArray&            obsDates,
                      DateTime&                 finalObsDate) const {
    static const string routine("XCB::getObsDates");
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
double XCB::pastValue(const DateTime&             sampleDate,
                      const ObservationType*      obsType,
                      const ObservationSource*    source,
                      const FixingType*           fixType,
                      const IObservationOverride* overrides,
                      const SamplingConvention*   sampleRule) const{
    static const string routine("XCB::pastValue");
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
                    throw ModelException(e, "XCB::pastValue", 
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
bool XCB::observationDate(const DateTime&           sampleDate,
                          const ObservationSource*  source,
                          const SamplingConvention* sampleRule,
                          DateTime*                 obsDate) const {
    // note for an XCB the source coming in is ignored. Sampling
    // is governed by the component sources defined in the XCB
    // function calculates the date at which the sample will finally be known
    static const string routine("XCB::observationDate");
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
double XCB::addPastSampleEvent(const DateTime&          sampleDate,
                            const ObservationType*      obsType,
                            const ObservationSource*    source,
                            const FixingType*           fixType,
                            const IObservationOverride* overrides,
                            const SamplingConvention*   sampleRule,
                            PastSamplesCollector*        collector) const {
    static const string method("XCB:addPastSampleEvent");
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
bool XCB::isHoliday(const DateTime&            sampleDate,
                    const ObservationSource*   source) const {
    static const string method("XCB:isHoliday");
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

/** calculate the vol interp scale factor for smile - this is the scaling that
    should be used for a non percentage strike (general idea is to let the
    vol interp decide whether or not to scale its strike) */
double XCB::volInterpSmileScale(
    double        basketSpot,   /* (I) the current basket spot price */
    int           i) const      /* (I) index of asset in compAsset */
{
    static const string routine("XCB:volInterpSmileScale");
    double scale;
    if (useSmileVolInterp()){
        bool   bskFwdStart = startDate.isGreater(baseDate);
        if (smile == floatStrike || bskFwdStart){
            /* either E type of fwd starting basket */
            double spot = assets[i]->getSpot();
            scale = spot/basketSpot;
        } else if (smile == fixedStrike) {
            // D type, xcb already started
            scale = spotsAtStart[i]/basketAtStart;
        } else {
            throw ModelException(routine, "Unknown smile type");
        }
    } else {
        scale = 1.0; // this value should be irrelevant
    }
    return scale;
}

CVolRequestLNArraySP XCB::getComponentVolRequests(
    const CVolRequestLN* volRequest) const
{
    static const string  routine("XCB::getComponentVolRequests");
    CVolRequestLNConstSP origInterp(
        CVolRequestLNConstSP::attachToRef(volRequest));
    // create an array of vol requests - one for each asset
    CVolRequestLNArraySP interps(new CVolRequestLNArray(assets.size()));
    bool useSmile = useSmileVolInterp();
    if (!useSmile){
        // override the interp
        origInterp = CVolRequestLNSP(new ATMVolRequest());
    }
    double basketSpot = getSpot(); // cache the basket spot price
    /* note that we scale each of the vol requests regardless of whether the
       basket is 'C' type or if fwd starting - instead we rely on the vol
       request to be able to correctly scale itself. This is because the
       XCB might be ccy protected and, if the xcb is fwd starting, we need to
       scale any absolute strike levels coming in */
    for (int i = 0; i < assets.size(); i++) {
        /* copy each of the interps so that we can scale as needed */
        (*interps)[i] = CVolRequestLNSP(origInterp.clone());
        // ensure no negative forward variance.
        // Overriding the setting in the request means an exception WILL be
        // thrown if neg fwd var is encountered. This is less severe than
        // using checkNegativeFwdVarDisallowed() since it will let through
        // cases where the data is ok.
        (*interps)[i]->allowNegativeFwdVar(false);

        double scale = volInterpSmileScale(basketSpot, i);
        (*interps)[i]->scale(scale);
    }
    return interps;
}

/** Generate a vol surface for the XCB using 'default' strikes */
VolSurfaceSP XCB::defaultVolSurface() const{
    static const string method("XCB::defaultVolSurface");
    try{
        if (pdfStrikes.empty()){
            // to avoid instability in the greeks reuse the same strikes
            // from the original pricing run. We assume that this method
            // will be called from the pricing run if it is called during
            // the tweaks
            pdfStrikes = CVolProcessedBS::defaultStrikes(baseDate, this, this);
        }
        const DoubleArray& strikes = pdfStrikes;

        // gather data needed for composite vol:

        // calculate correlation skew parameters on the fly...
        double  correlationSkew;
        double  correlationSkewPower;
        computeCorrelationSkewAndPower(correlationSkew, correlationSkewPower);        

        CorrTermDataSP data = CorrelationTerm::getCorrelationTermSqueezesAndExpiries(
            baseDate, assets.size(),
            corrObjects,
            corrTermArray);

        int numAssets = assets.size();
        DoubleMatrix corrMatrix(numAssets, numAssets);
        int pos = 0, i;
        for (i = 0; i < assets.size(); i++) {
            corrMatrix[i][i] = 1.0;
            for (int j = i + 1; j < assets.size(); j++, pos++) {
                if (Correlation::TYPE->isInstance(corrObjects[pos].get()))
                {
                    CorrelationSP p = CorrelationSP::dynamicCast(corrObjects[pos]);
                    corrMatrix[i][j] = p->getCorrelation();
                    corrMatrix[j][i] = corrMatrix[i][j];
                }
                else
                {
                    // TO DO: handle the case the correlation is post calibrated
                    throw ModelException(method, "correlation is calibrated");
                }
            }
        }

        data->correlationMatrix = CDoubleMatrixSP(new DoubleMatrix(corrMatrix));

        const CDoubleArray& weights  = weightsRef();
        // use composite vol object to cache fwd prices and also allow use
        // of more flexible methods
        CompositeVol  compositeVol(assets);
        DoubleMatrix  vols; // will hold vols across strikes and benchmarks
        const ExpiryArray*   bmExpiries = 0;
        for (i = 0; i < strikes.size(); i++){
            // build the volRequest for this strike
            LinearStrikeVolRequest volRequest(strikes[i], baseDate, baseDate,
                                              false /* not fwd starting */);
            // get the vol requests for the components
            CVolRequestLNArraySP interps =
                getComponentVolRequests(&volRequest);
            DoubleArray volAtDates =
                compositeVol.getVolsAtBMDates(&volRequest,
                                              this,
                                              startDate.isGreater(baseDate),
                                              startDate,
                                              name,
                                              timeMetric.get(),
                                              (*interps),
                                              weights,
                                              *data->correlationMatrix,
                                              *data->benchmarkTermExpiries,
                                              correlationSkew,
                                              correlationSkewPower,
                                              *data->shortTermSqueeze,
                                              *data->shortTermExpiries,
                                              *data->longTermSqueeze,
                                              *data->longTermExpiries);
            if (i == 0){
                // record bm dates
                bmExpiries = &compositeVol.getBMExpiries();
                vols = DoubleMatrix(strikes.size(), bmExpiries->size());
            } else {
                // sanity check
                if (bmExpiries->size() != volAtDates.size()){
                    throw ModelException(method,
                                         "Number of bm dates has changed");
                }
            }
            // copy vols over to double matrix
            for (int j = 0; j < bmExpiries->size(); j++){
                vols[i][j] = volAtDates[j];
            }
        }

        // note this is only an internal representation so we aren't ever
        // going to roll it through time. Thus we don't care if some
        // absolute dates occur after relative ones. We'll convert
        // everything to absolute dates to avoid VolSurface construction failing
        ExpiryArraySP exps(compositeVol.fixedExpiries());

        return VolSurfaceSP(new VolSurface(name,  timeMetric.get(),
                                           strikes, vols, exps.get(),
                                           baseDate));
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Sums correlationSkew and correlationSkewPower */
void XCB::computeCorrelationSkewAndPower(double&  correlationSkew,
                                         double&  correlationSkewPower) const{
    correlationSkew         = 0.0;
    correlationSkewPower    = 0.0;
    for (int i = 0; i < corrSkewNames.size(); i++){
        correlationSkew += corrSkewWeights[i] *
            corrSkewNames[i]->getCorrelationSkew();
        correlationSkewPower += corrSkewWeights[i] *
            corrSkewNames[i]->getCorrelationSkewPower();
    }
}

/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest. Implementation for LogNormal vols */
CVolProcessed * XCB::getProcessedVolLN(
    const CVolRequestLN* volRequest) const{
    static const string  routine("XCB::getProcessedVol");
    try{
        // gather data needed for composite vol:
        const CDoubleArray& weights  = weightsRef();

        // calculate correlation skew parameters on the fly...
        double  correlationSkew;
        double  correlationSkewPower;
        computeCorrelationSkewAndPower(correlationSkew, correlationSkewPower);

        // get the vol requests for the components
        CVolRequestLNArraySP interps =  getComponentVolRequests(volRequest);        

        CorrTermDataSP data = CorrelationTerm::getCorrelationTermSqueezesAndExpiries(
            baseDate, assets.size(),
            corrObjects, 
            corrTermArray);

        int numAssets = assets.size();
        DoubleMatrix corrMatrix(numAssets, numAssets);
        int pos = 0;
        for (int i = 0; i < numAssets; i++) {
            corrMatrix[i][i] = 1.0;
            for (int j = i + 1; j < numAssets; j++, pos++) {
                if (Correlation::TYPE->isInstance(corrObjects[pos].get()))
                {
                    CorrelationSP p = CorrelationSP::dynamicCast(corrObjects[pos]);
                    corrMatrix[i][j] = p->getCorrelation();
                    corrMatrix[j][i] = corrMatrix[i][j];
                }
                else
                {
                    // TO DO: handle the case the correlation is post calibrated
                    throw ModelException("XCB::getProcessedVolLN", "correlation is calibrated");
                }
            }
        }

        data->correlationMatrix = CDoubleMatrixSP(new DoubleMatrix(corrMatrix));

        // create and interpolate composite vol
        return CompositeVol::getProcessedVolCorrelationSkewAndTerm(
                                volRequest,
                                this,
                                startDate.isGreater(baseDate),
                                startDate,
                                name,
                                timeMetric.get(),
                                assets,
                                (*interps),
                                weights,
                                *data->correlationMatrix,
                                *data->benchmarkTermExpiries,
                                correlationSkew,
                                correlationSkewPower,
                                *data->shortTermSqueeze,
                                *data->shortTermExpiries,
                                *data->longTermSqueeze,
                                *data->longTermExpiries);
    } catch (exception& e){
        throw ModelException(e, routine, "Failed for XCB "+ getName());
    }
}

/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest. Implementation for 'local vols' */
CVolProcessed * XCB::getProcessedVolDVF(
    const CVolRequestDVF* volRequest) const{
    // to support local vol we choose a set of strikes and build a
    // vol surface for those strikes. We then spline that surface
    // and use that spline to support the local vol.
    // Create a vol surface with these strike
    VolSurfaceSP myVolSurface(defaultVolSurface());
    // Create spline
    smartPtr<VolSpline> volSpline(new VolSpline(*myVolSurface, getSpot()));
    // Pass volRequest to spline
    return volSpline->getProcessedVol(volRequest, this);
}
/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
CVolProcessed * XCB::getProcessedVol(
    const CVolRequest* volRequest) const
{
    // we delegate the work of doing a 'double dispatch' to VolProcessedDispatch
    // Here 'double dispatch' means choosing the method based upon the type
    // of the vol and the request
    return VolProcessedDispatch::dispatch(this, volRequest);
}

/** Calculates the expected spot price of the asset at the given date */
double XCB::fwdValue(const DateTime& date) const{
    double fwd = 0.0;
    try{
        const CDoubleArray& weights = weightsRef();
        for (int i = 0; i < assets.size(); i++){
            fwd += assets[i]->fwdValue(date) * weights[i];
        }
    } catch (exception& e){
        throw ModelException(e, "XCB::fwdValue", "Failed for XCB "+
                             getName());
    }
    return fwd;
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void XCB::fwdValue(const DateTimeArray& dates,
                              CDoubleArray&        result) const{
    try{
        const CDoubleArray& weights = weightsRef();
        AssetUtil::fwdValue(assets, weights, dates, result);
    } catch (exception& e){
        throw ModelException(e, "XCB::fwdValue", "Failed for XCB "+
                             getName());
    }
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void XCB::fwdValue(const DateTimeArray&     dates,
                   const FwdValueAlgorithm& algo,
                   CDoubleArray&            result) const{
    try{
        const CDoubleArray& weights = weightsRef();
        AssetUtil::fwdValue(assets, weights, dates, algo, result);
    } catch (exception& e){
        throw ModelException(e, "XCB::fwdValue", "Failed for XCB "+
                             getName());
    }
}

/** Returns the name (not the ISO code) of the asset ccy */
string XCB::getYCName() const {
    return basketYCName;
}

// /** Calculates the expected spot price of the asset at the given date if
//     the spot price had the given value spot on spotDate */
// double XCB::fwdFwd(const DateTime& spotDate,
//                               double          spot,
//                               const DateTime& fwdDate) const{
//     throw ModelException("XCB::fwdFwd", "Not yet done");
//     return 0;
// }

/** Calculate the settlement date associated with a given trade date */
DateTime XCB::settleDate(const DateTime& tradeDate) const{
    return AssetUtil::settleDate(assets, tradeDate);
}

/** Returns true if basket needs initial spots */
bool XCB::needInitialSpots() const {
    bool needed;
    bool fwdStart = startDate.isGreater(baseDate);
    if (fwdStart) {
        /* it's fwd starting so can't possibly need what we don't know */
        needed = false;
    } else {
        /* if it's %age weighted, need 'em to set weights */
        if (!unitWeights) {
            needed = true;
        } else {
            /* by now we must be a started unit basket, so only need
             * initial spots if we're vol type D */
            needed = smile == fixedStrike ? true : false;
        }
    }
    return needed;
}

/** Do we use smile when interpolating the vol */
bool XCB::useSmileVolInterp() const{
    return (smile != noSmile);
}


/** Returns the name of the XCB for Delta */
string XCB::sensName(const BasketSpot*) const{
    return name;
}

/** Shifts all components in the XCB for Delta */
TweakOutcome XCB::sensShift(const PropertyTweak<BasketSpot>& shift){
    try{
        double oldSpot = getSpot();
        RiskProperty<Spot>().tweakAllSubjects(IObjectSP::attachToRef(this),
                                               VoidSP(),
                                               shift.coefficient);

        return TweakOutcome(oldSpot,
                            oldSpot * (1 + shift.coefficient),
                            false); // don't shift any further
    }
    catch (exception& e){
        throw ModelException(e, "XCB::sensShift");
    }
}

/** Returns name identifying this object for SpotShift */
string XCB::sensName(SpotShift* shift) const{
    return name;
}
/** Shifts the object using given shift (see SpotShift::Shift)*/
bool XCB::sensShift(SpotShift* shift){
    try {
        // new shift with no market name, so shifts everything
        SpotShift newShift(shift->getShiftSize());
        for (int i = 0; i < assets.size(); i++){
            newShift.findAndShift(IObjectSP::attachToRef(&assets[i]),
                                  OutputNameConstSP());
        }
    }
    catch (exception& e){
        throw ModelException(e, "XCB::sensShift");
    }
    return false; // don't shift any further
}

/** Shifts the object using given shift (see Theta::Shift)*/
bool XCB::sensShift(Theta* shift)
{
    try
    {
        pdfStrikes.clear(); // use fresh strikes when doing theta
        DateTime rolledDate = shift->rollDate(baseDate);
        bool  recalcBasket  = false;
        bool  populateZeros = startDate.equals(baseDate);

        // was forward starting but not anymore
        if ( (startDate.isGreater(baseDate) && rolledDate.isGreaterOrEqual(startDate)) ||
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
            recalcBasket = true;
        }

        // otherwise no need to do anything special since intWeights are
        // recalculated everytime they are needed, so just roll the baseDate
        baseDate = rolledDate;

        // calculate derived basketAtStart value, if we have run over the
        // start date. We couldn't do this inside the if, as the baseDate
        // has not been set at that point and needInitialSpots refers to
        // baseDate
        if ( recalcBasket && needInitialSpots() ) {
            if (!unitWeights){
                basketAtStart = 100.0; // by definition
            } else {
                basketAtStart = 0.0;
                for (int i = 0; i < assets.size(); i++){
                    basketAtStart += pubWeights[i] * spotsAtStart[i];
                }
            }
        }

    }
    catch (exception& e)
    {
        throw ModelException(e, "XCB::sensShift (theta)");
    }
    return true; // tweak stuff inside
}

// DeltaSurface - don't reallt want it, but need a zero for it to work
string XCB::sensName(DeltaSurface* shift) const {
    return name;
}

bool XCB::sensShift(DeltaSurface* shift) {
    // must set the intitial value
    shift->setInitialValue(getSpot());
    // route through BasketSpot so our spot moves
    return sensShift(PropertyTweak<BasketSpot>(shift->getShiftSize())).
               tweakMembers();
}

/** returns sensitive strikes for a given vol request */
void XCB::getSensitiveStrikes(const CVolRequest* volRequest,
                              OutputNameConstSP outputName,
                              const SensitiveStrikeDescriptor& sensStrikeDesc,
                              DoubleArraySP sensitiveStrikes) const
{
    // get the vol requests for the components - shouldn't this method take
    // a CVolRequestLN to begin with?
    const CVolRequestLN* request = DYNAMIC_CAST(CVolRequestLN, volRequest);
    CVolRequestLNArraySP interps = getComponentVolRequests(request);

    for (int i = 0; i < assets.size(); i++) {
        // let the component take care of it's sensitive strikes
        assets[i]->getSensitiveStrikes((*interps)[i].get(),
                                       outputName,
                                       sensStrikeDesc,
                                       sensitiveStrikes);
    }
}

/** record forwards at maturity*/
void XCB::recordFwdAtMat(OutputRequest*  request,
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

void XCB::acceptStartDateCollector(const XCB* asset,
                                   StartDateCollector* collector)
{
    collector->startDateValidate(asset->startDate,
                                 asset->getName(),
                                 true);
    // check start dates for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void XCB::acceptNameCollector(const XCB* asset, AssetNameCollector* collector)
{
    collector->assetNameValidate(asset->getName());
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void XCB::acceptFutureCollector(const XCB* asset,
                                FutureExpiryCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void XCB::acceptValueDateCollector(const XCB* asset,
                                   CValueDateCollector* collector)
{
    collector->valueDateValidate(asset->baseDate, asset->getName());
}

void XCB::acceptCriticalDateCollector(const XCB* asset, CriticalDateCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < asset->assets.size() ; ++i)
    {
        asset->assets[i]->accept(collector);
    }
}

void XCB::acceptDeltaShift(const XCB*                asset,
                           ShiftSizeCollector* collector)
{
    ShiftSizeCollector::TAdjustmentType adjType =
        collector->getAdjustmentType();
    if ( adjType == ShiftSizeCollector::FWD_START_ADJUSTMENT )
    {
        // only worried about spot price - component assets deal with this
        for (int idx=0 ; idx<asset->assets.size() ; idx++)
        {
            asset->assets[idx]->accept(collector);
        }
    } else if ( adjType == ShiftSizeCollector::SPOT_START_ADJUSTMENT ) {
        /* depending on our smile type, we might be moving across the vol
           surface */
        const Delta*  deltaCtrl  = collector->getSensControl().get();
        asset->alterDeltaShiftSize(collector, deltaCtrl);
    }
}

void XCB::acceptImntCcy(const XCB* asset,
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
        throw ModelException(e, "XCB::acceptImntCcy",
                             "failed for XCB " + asset->getName());
    }
}

void XCB::acceptHoliday(const XCB* asset,
                        HolidayCollector* collector)
{
    collector->setHoliday(asset->marketHols.getSP());
}

void XCB::acceptWrapperNameCollector(const XCB* asset,
                                     WrapperNameCollector* collector)
{
    collector->addName(asset->getName());
}

double XCB::calculateMaxShiftSize(int tweakedAssetIndex,
                                  int i,
                                  const double& strike,
                                  const double& volSurfaceStrike,
                                  const double& weight,
                                  const double& assetSpot,
                                  const double& basketSpot) const
{
    // calculate the maximum delta shift size
    double maxShift = 0.0;

    if ( tweakedAssetIndex == i ) {
        double divisor = (strike - volSurfaceStrike*weight);
        if ( !Maths::isZero(divisor) )
        {
            maxShift = fabs((strike*assetSpot - volSurfaceStrike * basketSpot)/
                            divisor);
        }
    } else {
        if ( !Maths::isZero(weight) )
        {
            maxShift = fabs((strike * assetSpot /volSurfaceStrike-basketSpot)/
                            weight);
        }
    }
    return maxShift;
}

void XCB::alterDeltaShiftSize(ShiftSizeCollector* collector,
                              const Delta*        sensControl) const{
    // if the XCB is a component of another XCB we don't change
    // the delta shift
    if (smile != floatStrike )
    {
        return;
    }
    double strike = collector->getStrike();


    double maxDeltaUp = 0.0;   // stop compiler warning
    double maxDeltaDown = 0.0; // stop compiler warning

    // get the basket spot price
    double basketSpot = getSpot();

    // get the basket weightings
    const CDoubleArray& weights = weightsRef();

    double originalShift = sensControl->getShiftSize();
    double newShift      = originalShift;

    // get the spot price of the tweaked equity
    double tweakedAssetSpot  = 0.0;
    int    tweakedAssetIndex = 0;

    // retrieve the asset's spot price for the equity we're currently tweaking
    // create new Delta control in case the names have been set in the
    // original delta control
    Delta myDeltaCtrl(originalShift);
    int i;
    for (i=0 ; i<assets.size() ; ++i)
    {
        OutputNameArrayConstSP names(myDeltaCtrl.names(assets[i].get()));

        if ( names->size() ==1 &&
             (*names)[0].get()->equals(sensControl->getMarketDataName().get()))
        {
            tweakedAssetSpot  = assets[i]->getSpot();
            tweakedAssetIndex = i;
            break;  // found what we were looking for
        }
    }

    if ( Maths::isPositive(tweakedAssetSpot) )
    {
        for (i=0 ; i<assets.size() ; ++i)
        {
            const IObject*     obj = assets[i].get();
            const INextStrike* nextStrike = 
                dynamic_cast<const INextStrike*>(obj);

            // do nothing if asset is not derived from NextStrike interface
            if ( !nextStrike )
                continue;

            // scale strike for component
            double assetStrike = strike * volInterpSmileScale(basketSpot,i);

            double assetSpot   = assets[i]->getSpot();

            bool   upperOffSurface = false;
            double upperVolStrike = nextStrike->getNextStrike(assetStrike,
                                                              true,
                                                              upperOffSurface);


            if (!upperOffSurface)
            {
                maxDeltaUp = calculateMaxShiftSize(tweakedAssetIndex,
                                                   i,
                                                   strike,
                                                   upperVolStrike,
                                                   weights[tweakedAssetIndex],
                                                   assetSpot,
                                                   basketSpot);
            }

            bool   lowerOffSurface = false;
            double lowerVolStrike = nextStrike->getNextStrike(assetStrike,
                                                              false,
                                                              lowerOffSurface);
            if (!lowerOffSurface)
            {
                maxDeltaDown =
                    calculateMaxShiftSize(tweakedAssetIndex,
                                          i,
                                          strike,
                                          lowerVolStrike,
                                          weights[tweakedAssetIndex],
                                          assetSpot,
                                          basketSpot);
            }

            if ( !upperOffSurface && !lowerOffSurface )
            {
                newShift = Maths::min(Maths::min(
                    (fabs(maxDeltaUp)- 2.0 * DBL_EPSILON)/(2*tweakedAssetSpot),
                    (fabs(maxDeltaDown) - 2.0 * DBL_EPSILON )/
                    (2*tweakedAssetSpot)), newShift);
            } else if ( !upperOffSurface ) {
                newShift = Maths::min(
                    (fabs(maxDeltaUp)- 2.0 * DBL_EPSILON)/(2*tweakedAssetSpot),
                    newShift);
            } else if ( !lowerOffSurface ) {
                newShift = Maths::min(
                    (fabs(maxDeltaDown) - 2.0 * DBL_EPSILON )/
                    (2*tweakedAssetSpot),
                    newShift);
            }
        }
    }

    /* We only want to use the new shift if it's at least the minimum
       shift - otherwise we use the original shift size */
    if ( newShift > Delta::MINIMUM_SHIFT  &&
         newShift < originalShift )
    {
        collector->setShiftSize(newShift);
    }
}

/** Returns fair value of fund price */
double XCB::fairValue() const {
    return getSpot(); // do we want settlement ? need a yield curve if we do
}

/** return a pdf calculator */
PDFCalculator* XCB::pdfCalculator(const PDFRequest* request) const {
    if (!PDFRequestLNStrike::TYPE->isInstance(request)){
        throw ModelException("XCB::pdfCalculator", "Only LN pdf supported");
    }
    const PDFRequestLNStrike& lnRequest =
        dynamic_cast<const PDFRequestLNStrike&>(*request);
    PDFRequestLNStrikeConstSP lnRequestSP(copyIfRef(&lnRequest));
    double spot = getSpot();
    // for performance we choose some strikes, calculate a vol surface for
    // these. Then spline the vol surface and use the PDFParamLNStrike class
    // to do the work. This avoids umpteen calculations of fwds and vols.

    // Step 1 - create a vol surface with default strike
    VolSurfaceSP myVolSurface(defaultVolSurface());
    // Step 2 - create spline
    smartPtr<VolSpline> volSpline(new VolSpline(*myVolSurface, spot));
    // Step 3 - use PDFParamLNStrike class
    // need a CVolParam though
    CVolParamConstSP volParam(volSpline->getVolParam());
    smartConstPtr<VolSpline> volSplineConst(volSpline);
    return (new PDFParamLNStrike(baseDate, CAssetConstSP::attachToRef(this),
                                 timeMetric, volSplineConst,
                                 volParam, lnRequestSP));
}

/** sets correlation skew parameters */
void XCB::getCorrelationSkewParameters(const IModel* model, const MarketData* market){
    static const string method("XCB::getCorrelationSkewParameters");
    try {
        if ( !(corrSkewNames.empty() && corrSkewWeights.empty()) ){
            if( !(corrSkewNames.empty() || corrSkewWeights.empty()) ){
                if( corrSkewNames.size() == corrSkewWeights.size() ){

                    // validate
                    AssetUtil::checkWeights(corrSkewWeights,corrSkewWeights.size());

                    // calculate
                    for (int i = 0; i < corrSkewNames.size(); i++){
                        corrSkewNames[i].getData(model, market);
                    }
                }
                else {
                    throw ModelException(method,
                        "Provided " +
                        Format::toString(corrSkewNames.size()) +
                        " names and " +
                        Format::toString(corrSkewWeights.size()) +
                        " weights, should be equal.");
                }
            }
            else {
                throw ModelException(method,
                        "Correlation skew names or weights missing.");
            }
        }
    }
    catch (exception& e){
        throw ModelException(e, method, "Failed for XCB " + getName());
    }
}

/** Can this asset physically settle? */
bool XCB::canPhysicallySettle() const {
    return unitWeights;  // can't deliver % weighted basket
}


class XCBGetMarketAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    CAssetSP xcb;
    CMarketDataSP market;
    IModelSP model;

    // return the processed vol as a handle
    static IObjectSP getXCBMarket(XCBGetMarketAddin* params){
        CAssetSP cloneXCB(copy(params->xcb.get()));

        cloneXCB->getMarket(params->model.get(), params->market.get());

        return cloneXCB;
    }


    XCBGetMarketAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XCBGetMarketAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXCBGetMarketAddin);
        FIELD(xcb, "XCB");
        FIELD(market, "market");
        FIELD(model, "model");

        Addin::registerClassObjectMethod("XCB_GET_MARKET",
                                         Addin::UTILITIES,
                                         "Populates a XCB with it's market data",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getXCBMarket);
    }

    static IObject* defaultXCBGetMarketAddin(){
        return new XCBGetMarketAddin();
    }

};

CClassConstSP const XCBGetMarketAddin::TYPE =
CClass::registerClassLoadMethod("XCBGetMarketAddin",
                                typeid(XCBGetMarketAddin), load);



/* for reflection */
XCB::XCB(): CAsset(TYPE), useCorrelationTerm(false),
            smile(0), basketAtStart(0.0), pdfBoundaryProb(CVolProcessedBS::LEGACY_PDF_BOUNDARY_PROB) {}

class XCBHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(XCB, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(ITweakableWithRespectTo<BasketSpot>);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(SpotShift::Shift);
        IMPLEMENTS(IPDFBoundaryProb);
        EMPTY_SHELL_METHOD(defaultXCB);
        FIELD(name, "XCB's name");
        FIELD(assets, "Assets in XCB");
        FIELD(ccyTreatments, "Currency treatments");
        FIELD(basketYCName, "XCB's currency");
        FIELD(unitWeights, "Is XCB unit weighted");
        FIELD(pubWeights, "Units/Percentage of each asset");
        FIELD(marketHols, "Market hols for XCB");
        FIELD(baseDate, "today");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(startDate, "When XCB starts");
        FIELD(spotsAtStart, "Spots of asset at startDate");
        FIELD_MAKE_OPTIONAL(spotsAtStart);
        FIELD(smileType, "Type of smile");
        FIELD(correlations, "Equity-Equity Correlations");
        FIELD_MAKE_OPTIONAL(correlations);
        FIELD(timeMetric, "Trading time for composite vol");
        FIELD(corrSkewNames, "Name of correlation skew object");
        FIELD_MAKE_OPTIONAL(corrSkewNames);
        FIELD(corrSkewWeights, "Weight of correlation skew object");
        FIELD_MAKE_OPTIONAL(corrSkewWeights);
        FIELD(smile, "cached smile type");
        FIELD_MAKE_TRANSIENT(smile);
        FIELD(basketAtStart, "cached spot at start");
        FIELD_MAKE_TRANSIENT(basketAtStart);
        FIELD(intWeights, "cached weights");
        FIELD_MAKE_TRANSIENT(intWeights);
        FIELD(basketCcyCode, "XCB's currency code");
        FIELD_MAKE_TRANSIENT(basketCcyCode);
        FIELD(pdfStrikes, "pdfStrikes");
        FIELD_MAKE_TRANSIENT(pdfStrikes);
        FIELD(corrTermArray, "corrTermArray");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrTermArray);
        FIELD(corrObjects, "corrObjects");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrObjects);
        FIELD(useCorrelationTerm, "useCorrelationTerm");
        FIELD_MAKE_OPTIONAL(useCorrelationTerm);
        FIELD_USING_ALIAS(corrTermArrayNotUsed, corrTermArray, "DO NOT USE");
        FIELD_MAKE_OPTIONAL(corrTermArray);
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

        FIELD(pdfBoundaryProb, "limits implied pdf");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(pdfBoundaryProb);

        // register the flavours of vol request that we support
        VolProcessedDispatch::registerAssetMethod(&XCB::getProcessedVolLN);
        VolProcessedDispatch::registerAssetMethod(&XCB::getProcessedVolDVF);

        ClassSetAcceptMethod(XCB::acceptStartDateCollector);
        ClassSetAcceptMethod(XCB::acceptNameCollector);
        ClassSetAcceptMethod(XCB::acceptFutureCollector);
        ClassSetAcceptMethod(XCB::acceptValueDateCollector);
        ClassSetAcceptMethod(XCB::acceptDeltaShift);
        ClassSetAcceptMethod(XCB::acceptImntCcy);
        ClassSetAcceptMethod(XCB::acceptHoliday);
        ClassSetAcceptMethod(XCB::acceptWrapperNameCollector);
        ClassSetAcceptMethod(XCB::acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const XCB*               asset,
                                        DividendCollector* collector){
        const CDoubleArray& weights = asset->weightsRef();
        for (int idx=0 ; idx<asset->assets.size() ; idx++) {
            collector->processComponentAsset(asset->assets[idx].get(),
                                             weights[idx]);
        }
    }

    static IObject* defaultXCB(){
        return new XCB();
    }
};

CClassConstSP const XCB::TYPE = CClass::registerClassLoadMethod(
    "XCB", typeid(XCB), XCBHelper::load);


double XCB::getPDFBoundaryProb() const
{
    return pdfBoundaryProb;
}


DRLIB_END_NAMESPACE

