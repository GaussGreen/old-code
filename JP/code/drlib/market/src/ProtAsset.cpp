//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ProtAsset.cpp
//
//   Description : Generalized protected asset
//
//   Author      : Andrew J Swain
//
//   Date        : 2 October 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ProtAsset.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FlatFXVol.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/FXVegaPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

#define EDG_CCY_PROT_ATM_FACTOR 1.1

/** overrides CObject version to allow for easy default */
bool ProtAsset::accept(ICollector* collector) const{
    if (!CClass::invokeAcceptMethod(this, collector)){
        // if no method registered default by passing onto asset
        return asset->accept(collector);
    }
    return false;
}


/** Pull out the vol, fx asset and correlation from the
    market data */
void ProtAsset::getMarket(const IModel* model, const MarketData* market) {
    static const string method = "ProtAsset::getMarket";
    try {
        // get the data we need from the cache
        market->GetReferenceDate(baseDate);
        /* note that we don't populate the protYC wrapper to avoid
           generating rho sensitivities */
        if (protYC.get()){
            protYCCode = protYC->getCcy();
            protYC.setObject(MarketObjectSP()); // clear out to avoid greeks
            /* If there are issues with this approach we can use the
               ObjectIteration::IOverride interface to switch off tweaking for
               greeks */
        } else {
            YieldCurveWrapper protYCtmp(protYC.getName());
            protYCtmp.getData(model, market);
            // then cache transient field
            protYCCode = protYCtmp->getCcy();
        }
        // essentially do asset.getData(model, market) but tell model that the
        // ccy is changing
        CAsset::getAssetInNewCurrency(model, market, asset);
        // can look this up (as we've got the asset now)
        const string& ycName = asset->getYCName();

        // validate ISO codes are different
        if (protYCCode == AssetUtil::assetCcy(asset.get())) {
            throw ModelException(method,
                                 "Protected currency is identical to "
                                 "underlying currency (" +
                                 protYCCode+ ")");
        }

        // decide what fx vol we want here
        CClassConstSP fxVolType = FXVolBase::volClassProtStruck(model,
                                                                market,
                                                                fxVol.getName(),
                                                                fxVol.getMOType());
                                                           
        fxVol.getData(model, market, fxVolType, ycName);

        // ideally we get the FX asset so we know what ATM means for interpolating
        // on the FX vol surface. In reality we have 1000s of test files with no
        // FX asset and only a scalar FX vol, so need to not make these break
        // use getCorrFromCache as a proxy for old style files
        if (getCorrFromCache && fx.isEmpty()){
            // now get the FX itself - need this to handle FX vol surfaces
            string fxName = market->getFXName(ycName, protYC.getName());
            fx.setName(fxName);

            // this FXAsset method needs renaming now
            FXAsset::getMarketForStruck(fx, model, market);
        }

        if (getCorrFromCache){
            string growYCName = asset->getYCName();
            string fxName = market->getFXName(growYCName, protYC.getName());
            string corrName = market->getCorrelationName(asset->getName(),
                                                         fxName);
            corrAssetFX = CorrelationSP::dynamicCast(
                model->getCorrelation(corrName, asset->getClass(), 
                                      FXAsset::TYPE,
                                      Correlation::TYPE, market));
            validatePop2Object();
            getCorrFromCache = true; // restore value
        } else {
            // configure it correctly ourselves
            corrAssetFX->configureForSensitivities(asset->getClass(), 
                                                   FXAsset::TYPE);
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** returns the spot price */
double ProtAsset::getSpot() const{
    return scale * asset->getSpot();
}

/** returns the asset name */
string ProtAsset::getName() const{
    return name;
}

/** returns the 'true' name of the asset that is being protected */
string ProtAsset::getTrueName() const{
    return asset->getTrueName();
}

// the IMarketObservable interface for retrieving a single sample
double ProtAsset::pastValue(const DateTime&                sampleDate,
                               const ObservationType*      obsType,
                               const ObservationSource*    source,
                               const FixingType*           fixType,
                               const IObservationOverride* overrides,
                               const SamplingConvention*   sampleRule) const{
    return asset->pastValue(sampleDate, obsType, source, fixType,
                            overrides, sampleRule);
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool ProtAsset::observationDate(const DateTime&           sampleDate,
                                const ObservationSource*  source,
                                const SamplingConvention* sampleRule,
                                DateTime*                 obsDate) const {
    return asset->observationDate(sampleDate, source, sampleRule, obsDate);
}

// the IMarketObservable interface for retrieving past samples events
double ProtAsset::addPastSampleEvent(const DateTime&        sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const {
    return asset->addPastSampleEvent(sampleDate, obsType, source, fixType,
                                      overrides, sampleRule, collector);
}

// the IMarketObservable interface for 
// is the given date a holiday for the relevant source
bool ProtAsset::isHoliday(const DateTime&            sampleDate,
                             const ObservationSource*   source) const {
    return asset->isHoliday(sampleDate, source);
}

/** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
CVolProcessed * ProtAsset::getProcessedVol(
    const CVolRequest* volRequest) const{
    return asset->getProcessedVol(volRequest);
}

/** Calculate the settlement date associated with a given trade date */
DateTime ProtAsset::settleDate(const DateTime& tradeDate) const{
    return asset->settleDate(tradeDate);
}

// calculate protected adjustments
CDoubleArraySP ProtAsset::protAdjustment(
    const DateTimeArray& dates,
    const double         spotPrice,
    const CDoubleArray&  unprotectedFwds) const
{
    static const string routine("ProtAsset::protAdjustment");
    double              approxKInterp ;
    double              approxVolAtKInterp;

    CDoubleArraySP ccyAdjust(new CDoubleArray(dates.size()));

    /*  asset/fx correlation */
    double correlation = corrAssetFX->getCorrelation();

    /* interpolate the fx vols atm */
    ATMVolRequestSP fxVolRequest(new ATMVolRequest());
    // interpolate the vol
    CVolProcessedSP  fxVolPS(getProcessedFXVol(fxVolRequest.get()));
    // cast to the type of vol we're expecting
    CVolProcessedBSSP fxVolBS = CVolProcessedBSSP::dynamicCast(fxVolPS);

    /* interpolate the asset vols atm */
    CVolRequestLNSP atmVolRequest(new ATMVolRequest());
                                                             
    // interpolate the vol
    CVolProcessedSP  atmVol(getProcessedVol(atmVolRequest.get()));
    // cast to the type of vol we're expecting
    CVolProcessedBSSP atmVolBS = CVolProcessedBSSP::dynamicCast(atmVol);

    /* 
     * create a vol interp object to interpolate at 
     * atm * EDG_CCY_PROT_ATM_FACTOR. Note we do not want forward starting 
     * vol interpolation. 
     */
    LinearStrikeVolRequestSP otmVolRequest(new LinearStrikeVolRequest(
                             spotPrice * EDG_CCY_PROT_ATM_FACTOR, 
                             baseDate, 
                             baseDate,
                             false));
    // interpolate the vol
    CVolProcessedSP  otmVol(getProcessedVol(otmVolRequest.get()));
    // cast to the type of vol we're expecting
    CVolProcessedBSSP otmVolBS = CVolProcessedBSSP::dynamicCast(otmVol);

    double logOfOTMFactor = log(EDG_CCY_PROT_ATM_FACTOR);

    double varianceATM       = 0.0;
    double varianceFX        = 0.0;
    double volATM            = 0.0;
    double volOTM            = 0.0;
    double varAtApproxStrike = 0.0;

    for (int i=0 ;i < dates.size() ; ++i) {
        if (dates[i].equals(baseDate)) {
            (*ccyAdjust)[i] = 1.0;
        }
        else {

            // calculate the variances and vols 
            varianceATM = atmVolBS->CalcVar(baseDate, dates[i]);
            varianceFX  = fxVolBS->CalcVar(baseDate, dates[i]);
            volATM      = atmVolBS->CalcVol(baseDate, dates[i]);
            volOTM      = otmVolBS->CalcVol(baseDate, dates[i]);

            // calculate approxKInterp - an approximation of the strike
            // we would ideally interpolate at 
            approxKInterp = unprotectedFwds[i] *
                exp(( -correlation * sqrt( varianceATM * varianceFX )) +
                    ( 0.5 * varianceATM ));

            // approximate the volatility at the ideal interpolation level 
            // note that log(x) computes the natural logarithm of x 
            approxVolAtKInterp = volATM +
                ((volOTM - volATM) *
                 (log(approxKInterp / spotPrice) / logOfOTMFactor));

            // floor vol to 0.0 as vol can go negative if skew is negative and
            // the ratio of approximate strike and spot is large
            approxVolAtKInterp = Maths::max(approxVolAtKInterp, 0.0);

            varAtApproxStrike = 
                approxVolAtKInterp * approxVolAtKInterp *
                atmVol->calcTradingTime(baseDate, dates[i]);

            // calculate the currency protection factor 
            (*ccyAdjust)[i] = exp( -correlation * sqrt(varianceFX * varAtApproxStrike));
        }
    }
    return ccyAdjust;
}

/** Calculates the expected spot price of the asset at the given date */
double ProtAsset::fwdValue(const DateTime& date) const{

    double fwdPrice = asset->fwdValue(date);
    double spot     = getSpot();

    // set up array to avoid having two different methods for protected 
    // adjustments which do exactly the same
    DateTimeArray dates(1);
    CDoubleArray  fwdPrices(1);
    dates[0]     = date;
    fwdPrices[0] = fwdPrice;

    CDoubleArraySP ccyAdjust = protAdjustment(dates, spot, fwdPrices);
    return (fwdPrice * (*ccyAdjust)[0] * scale);
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void ProtAsset::fwdValue(const DateTimeArray& dates,
                         CDoubleArray&        result) const {
    int i = 0;

    if (dates.size() != result.size()) {
        throw ModelException("ProtAsset::fwdValue", 
            "date list and result array have different sizes");
    }

    // first calculate unprotected forward prices
    asset->fwdValue(dates, result);

    // get spot which is required for vol interpolation
    double spot = getSpot();

    // calculate protected adjustments
    CDoubleArraySP protAdjustment = this->protAdjustment(dates, spot, result);

    for (i = 0 ;i < dates.size(); ++i) {
        result[i] *= (*protAdjustment)[i] * scale;
    }
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void ProtAsset::fwdValue(const DateTimeArray&     dates,
                         const FwdValueAlgorithm& algo,
                         CDoubleArray&            result) const{
    try {
        int i = 0;

        if (dates.size() != result.size()) {
            throw ModelException("ProtAsset::fwdValue", 
                                 "date list and result array have different sizes");
        }

        // first calculate unprotected forward prices
        asset->fwdValue(dates, algo, result);

        // get spot which is required for vol interpolation
        double spot = getSpot();

        // calculate protected adjustments
        CDoubleArraySP protAdjustment = this->protAdjustment(dates, spot, result);

        for (i = 0 ;i < dates.size(); ++i) {
            result[i] *= (*protAdjustment)[i] * scale;
        }

    } catch (exception& e){
        throw ModelException(e, "ProtAsset::fwdValue", "Failed for "+
                             getName());
    }
}

/** Returns the name (not the ISO code) of the asset ccy */
string ProtAsset::getYCName() const {
    return protYC.getName();
}


/** Calculates the expected spot price of the underlying asset at
    the given date ie no currency protection adjustment is made */
double ProtAsset::unadjustedFwdValue(const DateTime& date) const{
    return asset->fwdValue(date);
}

/** Calculates the expected spot price of the underlying asset at
    the given dates ie no currency protection adjustment is made */
void ProtAsset::unadjustedFwdValue(const DateTimeArray& dates,
                                   CDoubleArray&        result) const{
    asset->fwdValue(dates, result);
} 

/** return the Eq/FX correlation */
const Correlation* ProtAsset::getCorrelation() const {
    return corrAssetFX.get();
}

/** Returns a processed vol - which combines the vol market data with
 *  the instrument data in the volRequest */
IVolProcessed* ProtAsset::getProcessedFXVol( const CVolRequest* volRequest ) const {
    if (!fx.get()) {
        throw ModelException("ProtAsset::getProcessedFXVol",
                "FX asset is null for "+getTrueName());
    }
    return fxVol->getProcessedVol(volRequest, fx.get());
}

/** Constructor needed for case when instrument specifies underlying asset and 
    currency treatment */
ProtAsset::ProtAsset(const string& name,
                     const Asset*  asset,
                     const string& fxVolName,
                     const string& protCcyCode,
                     bool homoGreeks): 
    CAsset(TYPE), name(name), asset(copy(asset)),
    fxVol(fxVolName), protYC(protCcyCode), coherentGreeks(homoGreeks) {
    validatePop2Object();
}

/** validation code to be called after object construction */
void ProtAsset::validatePop2Object()
{
    const string method("ProtAsset::validatePop2Object");
    // initialised cached asset name
    if (name.empty()) {
        name = CAsset::getSyntheticName(asset.getName(),
                                        CAsset::CCY_TREATMENT_PROTECTED,
                                        protYC.getName());
    } else if (asset.get()){
        // if market data is present, check names are different
        if (name == asset->getName()) {
            throw ModelException(method,
                                 "Protected asset's name must be different to"
                                 " asset's name");
        }
    }
    getCorrFromCache = !corrAssetFX;
    if (!getCorrFromCache){
        // hack for some tree products which incorrectly do not always
        // call getMarket. We'll make a guess as to what type the equity is
        // Should getMarket be called it will be done properly.
        corrAssetFX->configureForSensitivities(SimpleEquity::TYPE,
                                               FXAsset::TYPE);
    }
    scale = 1.0; // needs to be common to all constructors
}


/** Shifts the object using given shift. */
bool ProtAsset::sensShift(Theta* shift)
{
    try  {
        DateTime newDate = shift->rollDate(baseDate);
        
        if (shift->useAssetFwds()) {
            /* THETA_FORWARD_SPOT - get the ccy adjustment coefficient */
            scale = fwdValue(newDate) / asset->fwdValue(newDate);
        }

        baseDate = newDate;
    }
    catch (exception& e)
    {
        throw ModelException(e, "ProtAsset::sensShift (theta)"); 
    }
    return true; 
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string ProtAsset::sensName(VegaSkewParallel* shift) const{
     throw ModelException("ProtAsset::sensName(VegaSkewParallel *)",
                          "Vol name mechanism has been overridden.");
}

/** Matches the name of the vol - used to determine whether to tweak
    the object */
bool ProtAsset::sensNameMatches(VegaSkewParallel* shift, const OutputName& name) const {
    // We have to set up the spot when the FX vol is tweaked
    return name.equals(fxVol->getName());
}

/** Add the names of the vols that this underlying is sensitive to */
void ProtAsset::sensAppendName(VegaSkewParallel* shift, OutputNameArray& namesList) const {
    // Only register the FX if it's a surface that supports this sensitivity
    if(VegaSkewParallel::IShift::TYPE->isInstance(fxVol)) {
        OutputNameSP fxName(new OutputName(fxVol->getName()));
        namesList.push_back(fxName);
    }
}

bool ProtAsset::sensShift(VegaSkewParallel* shift){
    // store the underlying spot price, according to which vol surface is being tweaked
    // The asset, if it supports VegaSkewParallel, will look after itself.
    if( shift->getMarketDataName()->equals(fxVol->getName()) ) {
        shift->setSpot(fx->getSpot());
    }

    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string ProtAsset::sensName(VegaSkewPointwise* shift) const{
     throw ModelException("ProtAsset::sensName(VegaSkewPointwise *)",
                          "Vol name mechanism has been overridden.");
}

/** Matches the name of the vol - used to determine whether to tweak
    the object */
bool ProtAsset::sensNameMatches(VegaSkewPointwise* shift, const OutputName& name) const {
    // We have to set up the spot when the FX vol is tweaked
    return name.equals(fxVol->getName());
}

/** Add the names of the vols that this underlying is sensitive to */
void ProtAsset::sensAppendName(VegaSkewPointwise* shift, OutputNameArray& namesList) const {
    // Only register the FX if it's a surface that supports this sensitivity
    if(VegaSkewPointwise::IShift::TYPE->isInstance(fxVol)) {
        OutputNameSP fxName(new OutputName(fxVol->getName()));
        namesList.push_back(fxName);
    }
}

ExpiryArrayConstSP ProtAsset::sensExpiries(VegaSkewPointwise* shift) const{
    return ExpiryArrayConstSP(); // return null SP - the vol will do this
}

bool ProtAsset::sensShift(VegaSkewPointwise* shift){
    // store the underlying spot price, according to which vol surface is being tweaked
    // The asset, if it supports VegaSkewPointwise, will look after itself.
    if( shift->getMarketDataName()->equals(fxVol->getName()) ) {
        shift->setSpot(fx->getSpot());
    }

    return true; // continue on and do the vol
}

/** record forwards at maturity*/
void ProtAsset::recordFwdAtMat(OutputRequest*  request,
                               CResults*       results,
                               const DateTime& maturityDate) const
{
    double fwd = fwdValue(maturityDate);

    // wrap forward price into object
    CDoubleSP fwdObj(CDouble::create(fwd));

    results->storeRequestResult(request,
                                fwdObj,
                                OutputNameSP(new OutputName(getName())));

    // record forward at maturity for raw asset
    asset->recordFwdAtMat(request, results, maturityDate);    
}
                         
void ProtAsset::acceptValueDateCollector(const ProtAsset*     asset, 
                                         CValueDateCollector* collector)
{
    collector->valueDateValidate(asset->baseDate, asset->getName());
}

void ProtAsset::acceptImntCcy(const ProtAsset*   asset,
                              AssetCcyCollector* collector)
{
    // check our ccy is ok
    collector->currencyValidate(asset->protYCCode, asset->getName());
    // then make sure our asset that we're protecting is ok
    AssetCcyCollector collect;
    asset->asset->accept(&collect);
}

void ProtAsset::acceptNameCollector(
    const ProtAsset* asset, AssetNameCollector* collector)
{
    collector->assetNameValidate(asset->getName());
    asset->asset->accept(collector);
}

/** returns sensitive strikes for a given vol request */
void ProtAsset::getSensitiveStrikes(
    const CVolRequest*               volRequest,
    OutputNameConstSP                outputName,
    const SensitiveStrikeDescriptor& sensStrikeDesc,
    DoubleArraySP                    sensitiveStrikes) const
{
    static const string method = "ProtAsset::getSensitiveStrikes";
    try {
        if( outputName->equals(fxVol->getName())) {
            if (sensStrikeDesc.forwardOnly == false){
                // Create an atm vol request for the FX side of things
                CVolRequestLNSP atmRequest(new ATMVolRequest());
                atmRequest->getSensitiveStrike(fx->getSpot(),
                                               sensitiveStrikes);
            }
        } else {
            if (sensStrikeDesc.forwardOnly == false){
                asset->getSensitiveStrikes(volRequest,
                                           outputName,
                                           sensStrikeDesc,
                                           sensitiveStrikes);
            }

            /** for the prot fwd adjustment we don't want the asset to stop things
                if the forward only piece is turned on */
            SensitiveStrikeDescriptor cmpStrikeDesc;
            cmpStrikeDesc.forwardOnly = false;

            // the strike which corresponds to the equity's ATM level
            ATMVolRequestSP atmVolRequest(new ATMVolRequest());
            // need to route this through our asset (eg it might be an xcb)
            asset->getSensitiveStrikes(atmVolRequest.get(),
                                       outputName,
                                       cmpStrikeDesc,
                                       sensitiveStrikes);
        

            // the strike which corresponds to the equity's 
            // ATM*EDG_CCY_PROT_ATM_FACTOR level
            LinearStrikeVolRequestSP otmVolRequest(new LinearStrikeVolRequest(
                asset->getSpot() * EDG_CCY_PROT_ATM_FACTOR, 
                baseDate, 
                baseDate,
                false));

            // need to route this through our asset (eg it might be an xcb)
            asset->getSensitiveStrikes(otmVolRequest.get(),
                                       outputName,
                                       cmpStrikeDesc,
                                       sensitiveStrikes);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}                       

/** given a current spot level, get the next strike on the vol surface where 
    the slope is non-differentiable */
double ProtAsset::getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const
{
    double volStrike;
    const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(asset.get());
    if (nextStrike) {
        volStrike = nextStrike->getNextStrike(strike,
                                              isUp,
                                              offSurface);
    }
    else {
        volStrike  = 0.0;
        offSurface = true;
    }
    return volStrike;
}

/** return a pdf calculator */
PDFCalculator* ProtAsset::pdfCalculator(const PDFRequest* request) const {
    return asset->pdfCalculator(request);
}


/** returns dividend list - to be retired.
    Needed because of crap method in AssetUtil that casts around desperately
    looking for a type that it knows about */
    /** returns dividend list - to be retired - it's a bit meaningless */
DividendListSP  ProtAsset::getAllDivsBetweenDates(const DateTime& start,
                                                  const DateTime& end) const{
    return AssetUtil::getAllDivsBetweenDates(asset.get(), start, end);
}

//// Returns the ccy treatment for the asset. Default returns 
string ProtAsset::getCcyTreatment() const{
    return CCY_TREATMENT_PROTECTED;
}


/** Can this asset physically settle? */
bool ProtAsset::canPhysicallySettle() const {
    return false;
}

// Get (const) underlying (non-protected) asset 
CAssetConstSP  ProtAsset::getPlainAsset() const {
    return CAssetConstSP(dynamic_cast<const CAsset *>(asset.get()));
}

static CFieldConstSP fxVolField;

bool ProtAsset::recurse(const CFieldConstSP& field,
                        const CClassConstSP& targetClass) const{
    // this gets called as part of the tweaking and allows us to specify 
    // whether the fields within this class should be tweaked or not. The
    // target class indicates what is being shifted

    // don't want to tweak any of our transient fields
    if (field->isTransientForIteration()) {
        return false;
    }

    if (coherentGreeks) {
        // Turn off FX-specific greeks
        if (field == fxVolField) {
            if (targetClass == FXVega::IShift::TYPE                          ||
                targetClass == FXVega::IRestorableShift::TYPE                ||
                targetClass == FXVegaPointwise::IShift::TYPE                 ||
                targetClass == FXVegaPointwise::IRestorableShift::TYPE) {
                return false;
            }
  
        }
    } else {
        if (field == fxVolField) {
        // Turn off homogeneous greeks
            if (targetClass == ITweakableWithRespectTo<VolParallel>::TYPE                    ||
                targetClass == IRestorableWithRespectTo<VolParallel>::TYPE          ||
                targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE    ||
                targetClass == IRestorableWithRespectTo<VolPointwise>::TYPE   ||
                targetClass == VegaMatrix::IShift::TYPE                     ||
                targetClass == VegaMatrix::IRestorableShift::TYPE           ||
                targetClass == RootTimeVega::IShift::TYPE                   ||
                targetClass == RootTimeVega::IRestorableShift::TYPE         ||
                targetClass == VegaSkewParallel::IShift::TYPE               ||
                targetClass == VegaSkewParallel::IRestorableShift::TYPE     ||
                targetClass == VegaSkewPointwise::IShift::TYPE              ||
                targetClass == VegaSkewPointwise::IRestorableShift::TYPE    ||
                targetClass == VolLevel::IShift::TYPE                       || 
                targetClass == VolParallelShift::IShift::TYPE               ||
                targetClass == VolBenchmarkShift::IShift::TYPE              ||
                targetClass == PowerVega::IShift::TYPE                      ||
                targetClass == DeltaSurface::IShift::TYPE                   ||
                targetClass == DeltaSurface::IRestorableShift::TYPE         ||
                targetClass == VolRelativeShift::IShift::TYPE               ||
                targetClass == VolAbsoluteShift::IShift::TYPE) {
                return false;
            }
        }
    }
    return true;
}

/* for reflection */
ProtAsset::ProtAsset(): CAsset(TYPE), coherentGreeks(false) {}

class ProtAssetHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ProtAsset, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(Asset::IQuanto);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(INextStrike);
        IMPLEMENTS(ObjectIteration::IOverride);
        FIELD(name, "Asset's name");
        FIELD_MAKE_OPTIONAL(name);
        EMPTY_SHELL_METHOD(defaultProtAsset);
        FIELD(asset,       "Asset");
        FIELD(fxVol,       "FX Volatility");
        FIELD(baseDate,    "Value Date");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(corrAssetFX,        "Asset FX Correlation");
        FIELD_MAKE_OPTIONAL(corrAssetFX);
        FIELD(protYC, "Protected Ccy");
        FIELD(protYCCode, "Protected Ccy iso code");
        FIELD_MAKE_TRANSIENT(protYCCode); // hide from dd interface
        FIELD(getCorrFromCache, "Cached flag");
        FIELD_MAKE_TRANSIENT(getCorrFromCache); // hide from dd interface
        FIELD(scale, "Cached ccy adjustment");
        FIELD_MAKE_TRANSIENT(scale);            // hide from dd interface
        FIELD(coherentGreeks, "Use homogeneous greeks for FX");
        FIELD_MAKE_TRANSIENT(coherentGreeks);
        FIELD(fx, "");
        FIELD_MAKE_TRANSIENT(fx);

        ClassSetAcceptMethod(ProtAsset::acceptValueDateCollector);
        ClassSetAcceptMethod(ProtAsset::acceptNameCollector);
        ClassSetAcceptMethod(ProtAsset::acceptImntCcy);
        ClassSetAcceptMethod(acceptDividendCollector);

        // look up field for use on recurse
        fxVolField = clazz->getDeclaredField("fxVol");
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const ProtAsset*   asset,
                                        DividendCollector* collector){
        collector->processComponentAsset(asset->asset.get(), 1.0);
    }

    static IObject* defaultProtAsset(){
        return new ProtAsset();
    }
};

CClassConstSP const ProtAsset::TYPE = CClass::registerClassLoadMethod(
    "ProtAsset", typeid(ProtAsset), ProtAssetHelper::load);

DRLIB_END_NAMESPACE
