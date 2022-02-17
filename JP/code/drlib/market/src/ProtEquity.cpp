//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketObject.cpp
//
//   Description : Base class for objects stored in market data cache
//
//   Author      : Mark A Robson
//
//   Date        : 19 Mar 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AllStrikes.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/FlatFXVol.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/FXVegaPointwise.hpp"
#include "edginc/VolDeltaShiftSize.hpp"

DRLIB_BEGIN_NAMESPACE

#define EDG_CCY_PROT_ATM_FACTOR 1.1

/** Pull out the vol, fx asset and correlation from the
    market data */
void ProtEquity::getMarket(const IModel* model, const MarketData* market){
    static const string method = "ProtEquity::getMarket";
    try{
        // get the data we need from the cache
        market->GetReferenceDate(baseDate);

        // need to identify new 'domestic' currency
        const string& ycName = equity->getYCName();

        vol.getData(model, market, ycName);

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
        if (getCorrFromCache){
            // now get the FX itself - need this to handle FX vol surfaces
            string fxName = market->getFXName(ycName, protCcyCode);
            fx.setName(fxName);

            // this FXAsset method needs renaming now
            FXAsset::getMarketForStruck(fx, model, market);
        }

        // need equity to get its yield curve before we can start querying it
        // Note that we really ought to switch the 'domestic yc name' in the 
        // model here (cf treatment of vols etc below) or pass the yc name
        // to the equity (although currently it doesn't need it)
        equity->getMarket(model, market);
        if (getCorrFromCache){
            // block ought to be outside of cache but we need a method to
            // see of yc exists in cache or not - to do
            YieldCurveWrapper yc(protCcyCode);
            yc.getData(model, market, ycName);
            protCcyIsoCode = yc->getCcy();
            // end of block

            string fxName = market->getFXName(ycName, protCcyCode);
            string corrName = market->getCorrelationName(equity->getName(),
                                                         fxName);
            corrEqFX = CorrelationSP::dynamicCast(
                model->getCorrelation(corrName, SimpleEquity::TYPE,
                                      FXAsset::TYPE,
                                      Correlation::TYPE, market));
            if ( protCcyIsoCode == equity->getYCIsoCode() ) {
                throw ModelException("method",
                                     "Protected currency is identical to "
                                     "underlying currency (" +
                                     equity->getYCIsoCode() + ")");
            }

        } else {
            // configure it correctly ourselves
            corrEqFX->configureForSensitivities(SimpleEquity::TYPE,
                                                FXAsset::TYPE);
        }
    } catch (exception& e){
        throw ModelException(e, method, "Failed for equity " + name);
    }
}

/** get a (copy of) plain asset */
CAssetConstSP  ProtEquity::getPlainAsset() const{
    CAssetConstSP asset = CAssetConstSP(new SimpleEquity(equity.get(), vol.get()));
    return asset;
}

// access the FX asset
const Asset* ProtEquity::getFXAsset() const {
    if (!fx.get()) {
        throw ModelException("ProtEquity::getFXAsset", "FX asset is null for "+
                             getTrueName());
    }
    return fx.get();
}

/** returns the spot price */
double ProtEquity::getSpot() const{
    return scale * equity->spot();
}

/** returns the asset name */
string ProtEquity::getName() const{
    return name;
}

/** returns the equity's name */
string ProtEquity::getTrueName() const{
    return equity->getName();
}

/** returns the name of the vol base object */
string ProtEquity::getVolName() const {
    return vol.getName();
}

/** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
CVolProcessed * ProtEquity::getProcessedVol(
    const CVolRequest* volRequest) const{
    return vol->getProcessedVol(volRequest, this);
}

/** Calculate the settlement date associated with a given trade date */
DateTime ProtEquity::settleDate(const DateTime& tradeDate) const{
    return equity->settles(tradeDate);
}

// calculate protected adjustments
CDoubleArraySP ProtEquity::protAdjustment(
    const DateTimeArray& dates,
    const double         spotPrice,
    const CDoubleArray&  unprotectedFwds) const
{
    static const string routine("ProtEquity::protAdjustment");
    double        approxKInterp ;
    double        approxVolAtKInterp;

    CDoubleArraySP ccyAdjust(new CDoubleArray(dates.size()));

    /*  equity/fx correlation */
    double correlation = corrEqFX->getCorrelation();

    /* interpolate the fx vols atm */
    ATMVolRequestSP fxVolRequest(new ATMVolRequest());
    // interpolate the vol
    // we use the FXAsset to define what ATM means
    if (fx.isEmpty() && getCorrFromCache){
        throw ModelException(routine, "FX asset is null for " + getTrueName());
    }
    // FX might be null, but that should only occur if FX vol is a scalar
    CVolProcessedSP  fxVolPS(fxVol->getProcessedVol(fxVolRequest.get(), fx.get()));
    // cast to the type of vol we're expecting
    CVolProcessedBSSP fxVolBS = CVolProcessedBSSP::dynamicCast(fxVolPS);

    /* interpolate the asset vols atm */
    ATMVolRequestSP atmVolRequest(new ATMVolRequest());
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

    for ( int i=0 ; i<dates.size() ; ++i )
    {
        double ttimeYearFrac = atmVol->calcTradingTime(baseDate, dates[i]);
        if ( dates[i].equals(baseDate) || Maths::isZero(ttimeYearFrac))
        {
            (*ccyAdjust)[i] = 1.0;
        }
        else
        {

            /* calculate the variances and vols */
            varianceATM = atmVolBS->CalcVar(baseDate, dates[i]);
            varianceFX  =  fxVolBS->CalcVar(baseDate, dates[i]);
            volATM      = sqrt(varianceATM/ttimeYearFrac);
            volOTM      = otmVolBS->CalcVol(baseDate, dates[i]);

            /* calculate approxKInterp - an approximation of the strike
               we would ideally interpolate at */
            approxKInterp = unprotectedFwds[i] *
                exp(( -correlation * sqrt( varianceATM * varianceFX )) +
                    ( 0.5 * varianceATM ));

            /* approximate the volatility at the ideal interpolation level 
               note that log(x) computes the natural logarithm of x */
            approxVolAtKInterp = volATM +
                ((volOTM - volATM) *
                 (log(approxKInterp / spotPrice) / logOfOTMFactor));

            // floor vol to 0.0 as vol can go negative if skew is negative and
            // the ratio of approximate strike and spot is large
            approxVolAtKInterp = Maths::max(approxVolAtKInterp, 0.0);

            varAtApproxStrike = 
                approxVolAtKInterp * approxVolAtKInterp * ttimeYearFrac;

            /* calculate the currency protection factor */
            (*ccyAdjust)[i] = exp( -correlation * 
                                   sqrt(varianceFX * varAtApproxStrike));
        }
    }
    return ccyAdjust;
}

/** Returns the name (not the ISO code) of the asset ccy */
string ProtEquity::getYCName() const {
    return protCcyCode;
}

// /** Calculates the expected spot price of the asset at the given date if
//         the spot price had the given value spot on spotDate */
// double ProtEquity::fwdFwd(const DateTime& spotDate,
//                           double          spot, 
//                           const DateTime& fwdDate) const{
//     throw ModelException("ProtEquity::FwdFwd", "Not yet done");
//     return 0;
// }


YieldCurveWrapper  ProtEquity::getYC() const{
    return equity->getYC();
}

/** Calculates the expected spot price of the asset at the given date */
double ProtEquity::fwdValue(const DateTime& date) const{

    double fwdPrice  = equity->fwdValue(date);
    double spot      = getSpot();

    // set up array to avoid having two different methods for protected 
    // adjustments which do exactly the same
    DateTimeArray    dates(1);
    CDoubleArray     fwdPrices(1);
    dates[0]     = date;
    fwdPrices[0] = fwdPrice;

    CDoubleArraySP ccyAdjust = protAdjustment(dates, spot, fwdPrices);
    return (fwdPrice * (*ccyAdjust)[0] * scale);
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void ProtEquity::fwdValue(const DateTimeArray& dates,
                          CDoubleArray&        result) const {
    FwdValueAlgorithm algo(false);
    fwdValue(dates, algo, result);
}

/** return the VolBaseWrapper */
const CVolBaseWrapper& ProtEquity::getVol()const
{
    return vol;
}
    
/** return the FXVolBaseWrapper */
const FXVolBaseWrapper& ProtEquity::getFXVol()const
{
    return fxVol;
}

/** Is the Eq/FX correlation in the cache ? */
bool ProtEquity::useCorrFromCache()const
{
    return getCorrFromCache;
}

/** return the Eq/FX correlation */
const Correlation* ProtEquity::getCorrelation() const
{
    return corrEqFX.get();
}

/** Returns a processed vol - which combines the vol market data with
 *  the instrument data in the volRequest */
IVolProcessed* ProtEquity::getProcessedFXVol( const CVolRequest* volRequest ) const
{
    if (!fx.get()) {
        throw ModelException("ProtEquity::getProcessedFXVol",
                "FX asset is null for "+getTrueName());
    }
    return fxVol->getProcessedVol(volRequest, fx.get());
}

/** Calculates the expected spot prices of the asset at the given dates
    respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
void ProtEquity::fwdValue(const DateTimeArray&     dates,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const{
    int     i = 0;

    if ( dates.size() != result.size() )
    {
        throw ModelException("Equity::fwdValue", 
                             "date list and result array have different sizes");
        return;
    }

    // first calculate unprotected forward prices
    equity->fwdValue(dates, algo, result);

    // get spot which is required for vol interpolation
    double spot = getSpot();

    // calculate protected adjustments
    CDoubleArraySP protAdjustment = this->protAdjustment(dates, spot, result);

    for ( i=0 ; i<dates.size() ; ++i )
    {
        result[i] *= (*protAdjustment)[i] * scale;
    }
}    

/** Calculates the expected spot price of the underlying asset at
    the given date ie no currency protection adjustment is made */
double ProtEquity::unadjustedFwdValue(const DateTime& date) const{
    return equity->fwdValue(date);
}

/** Calculates the expected spot price of the underlying asset at
    the given dates ie no currency protection adjustment is made */
void ProtEquity::unadjustedFwdValue(const DateTimeArray& dates,
                                    CDoubleArray&        result) const{
    equity->fwdValue(dates, result);
} 


/** Constructor needed for case when instrument specifies underlying asset and 
    currency treatment */
ProtEquity::ProtEquity(const string&        name,
                       const Equity*        equity,
                       const string&        eqVolName,
                       const string&        fxVolName,
                       const string&        protCcyCode,
                       bool                 homoGreeks): 
    CAsset(TYPE), name(name), equity(copy(equity)),
    vol(eqVolName), fxVol(fxVolName), protCcyCode(protCcyCode), coherentGreeks(homoGreeks) {
    validatePop2Object(); //sets scale
}


ProtEquity::ProtEquity(const string&        name,
                       const Equity*        equity,
                       const CVolBaseWrapper  vol,
                       const FXVolBaseWrapper fxVol,
                       const Correlation*         corrEqFX,
                       const DateTime&             baseDate,
                       const string&           protCcyIsoCode,
                       const string&        protCcyCode,
                       const double         scale,
                       bool                 homoGreeks): 
    CAsset(TYPE), name(name), equity(copy(equity)),
    vol(vol), fxVol(fxVol), corrEqFX(copy(corrEqFX)),
    baseDate(baseDate),protCcyCode(protCcyCode),scale(scale),
    protCcyIsoCode(protCcyIsoCode), coherentGreeks(homoGreeks) {}

/** validation code to be called after object construction */
void ProtEquity::validatePop2Object()
{
    // initialised cached asset name
    if (name.empty()){
        name = CAsset::getSyntheticName(equity->getName(),
                                        CAsset::CCY_TREATMENT_PROTECTED,
                                        protCcyCode);
    } else if (name == equity->getName()){
        throw ModelException("ProtEquity::validatePop2Object",
                             "Protected asset's name must be different from"
                             " equity's name");
    }
    if ( protCcyIsoCode.empty() )
    {
        protCcyIsoCode = protCcyCode;
    }

    getCorrFromCache = !corrEqFX;
    if (!getCorrFromCache){
        // hack for some tree products which incorrectly do not always
        // call getMarket. We'll guess the types for now - if getMarket
        // is called it will be done correctly
        corrEqFX->configureForSensitivities(SimpleEquity::TYPE, FXAsset::TYPE);
    }
    scale = 1.0;
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string ProtEquity::sensName(VegaSkewParallel* shift) const{
     throw ModelException("ProtEquity::sensName(VegaSkewParallel *)",
                          "Vol name mechanism has been overridden.");
}

/** Matches the name of the vol - used to determine whether to tweak
    the object */
bool ProtEquity::sensNameMatches(VegaSkewParallel* shift, const OutputName& name) const {
    return name.equals(vol->getName()) || name.equals(fxVol->getName());
}

/** Add the names of the vols that this underlying is sensitive to */
void ProtEquity::sensAppendName(VegaSkewParallel* shift, OutputNameArray& namesList) const {
    OutputNameSP eqName(new OutputName(vol->getName()));
    namesList.push_back(eqName);
    if(VegaSkewParallel::IShift::TYPE->isInstance(fxVol)) {
        OutputNameSP fxName(new OutputName(fxVol->getName()));
        namesList.push_back(fxName);
    }
}

bool ProtEquity::sensShift(VegaSkewParallel* shift){
    // store the underlying spot price, according to which vol surface is being tweaked
    if( shift->getMarketDataName()->equals(vol->getName()) ) {
        shift->setSpot(equity->spot());
    } else if( shift->getMarketDataName()->equals(fxVol->getName()) ) {
        shift->setSpot(fx->getSpot());
    }

    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string ProtEquity::sensName(VegaSkewPointwise* shift) const{
     throw ModelException("ProtEquity::sensName(VegaSkewPointwise *)",
                          "Vol name mechanism has been overridden.");
}

/** Matches the name of the vol - used to determine whether to tweak
    the object */
bool ProtEquity::sensNameMatches(VegaSkewPointwise* shift, const OutputName& name) const {
    return name.equals(vol->getName()) || name.equals(fxVol->getName());
}

/** Add the names of the vols that this underlying is sensitive to */
void ProtEquity::sensAppendName(VegaSkewPointwise* shift, OutputNameArray& namesList) const {
    OutputNameSP eqName(new OutputName(vol->getName()));
    namesList.push_back(eqName);
    if(VegaSkewPointwise::IShift::TYPE->isInstance(fxVol)) {
        OutputNameSP fxName(new OutputName(fxVol->getName()));
        namesList.push_back(fxName);
    }
}

ExpiryArrayConstSP ProtEquity::sensExpiries(VegaSkewPointwise* shift) const{
    return ExpiryArrayConstSP(); // return null SP - the vol will do this
}

bool ProtEquity::sensShift(VegaSkewPointwise* shift){
    // store the underlying spot price, according to which vol surface is being tweaked
    if( shift->getMarketDataName()->equals(vol->getName()) ) {
        shift->setSpot(equity->spot());
    } else if( shift->getMarketDataName()->equals(fxVol->getName()) ) {
        shift->setSpot(fx->getSpot());
    }

    return true; // continue on and do the vol
}

/** Shifts the object using given shift. */
bool ProtEquity::sensShift(Theta* shift)
{
    try
    {
        DateTime newDate = shift->rollDate(baseDate);

        if (shift->useAssetFwds())
        {
            /* THETA_FORWARD_SPOT - get the ccy adjustment coefficient */
            scale = fwdValue(newDate) / equity->fwdValue(newDate);
        }

        baseDate = newDate;
    }
    catch (exception& e)
    {
        throw ModelException(e, "ProtEquity::sensShift (theta)"); 
    }
    return true; 
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string ProtEquity::sensName(DeltaSurface* shift) const{
    return equity->getName();
}
    
bool ProtEquity::sensShift(DeltaSurface* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot(), equity->getName(), vol->getName());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string ProtEquity::sensName(VolRelativeShift* shift) const{
    return VolRelativeShift::IShift::TYPE->isInstance(vol.get())
        ? vol->getName(): "";
}
    
bool ProtEquity::sensShift(VolRelativeShift* shift){
    // store the underlying spot price
    shift->setSpot(baseDate, this);
    return true; // continue on and do the vol
}


/** returns sensitive strikes for a given vol request */
void ProtEquity::getSensitiveStrikes(
    const CVolRequest* volRequest,
    OutputNameConstSP outputName,
    const SensitiveStrikeDescriptor& sensStrikeDesc,
    DoubleArraySP sensitiveStrikes) const
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
            CAssetConstSP plainAsset = getPlainAsset();
            if (sensStrikeDesc.forwardOnly == false){
                plainAsset->getSensitiveStrikes(volRequest,
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
            plainAsset->getSensitiveStrikes(atmVolRequest.get(),
                                       outputName,
                                       cmpStrikeDesc,
                                       sensitiveStrikes);
        

            // the strike which corresponds to the equity's 
            // ATM*EDG_CCY_PROT_ATM_FACTOR level
            LinearStrikeVolRequestSP otmVolRequest(new LinearStrikeVolRequest(
                equity->spot() * EDG_CCY_PROT_ATM_FACTOR, 
                baseDate, 
                baseDate,
                false));

            // need to route this through our asset (eg it might be an xcb)
            plainAsset->getSensitiveStrikes(otmVolRequest.get(),
                                       outputName,
                                       cmpStrikeDesc,
                                       sensitiveStrikes);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}                                   

void ProtEquity::acceptValueDateCollector(const ProtEquity*    asset, 
                                          CValueDateCollector* collector)
{
    collector->valueDateValidate(asset->baseDate, asset->getName());
}

void ProtEquity::acceptDeltaShift(const ProtEquity*   asset, 
                                  ShiftSizeCollector* collector)
{
    const IObject* obj = asset->vol.get();
    const IVolDeltaShiftSize* volDelShift = 
                        dynamic_cast<const IVolDeltaShiftSize*>(obj);
    if (volDelShift) {
        volDelShift->adjustDeltaShiftSize(collector,
                                    asset->getName(),
                                    asset->getSpot());
    }
}

void ProtEquity::acceptImntCcy(const ProtEquity*  asset,
                               AssetCcyCollector* collector)
{
    collector->currencyValidate(asset->protCcyIsoCode,
                                asset->getName());
}

void ProtEquity::acceptHoliday(const ProtEquity* asset,
                               HolidayCollector* collector)
{
    collector->setHoliday(asset->equity->getMarketHolidays());
}

void ProtEquity::acceptCriticalDateCollector(const ProtEquity*      asset,
                                             CriticalDateCollector* collector)
{
/** Defining CORRECT_YIELD_DIVS will mean that yield divs are adjusted in the
    same way as the fwd price. Not switched on yet as role of critical date
    collector under way */
#define xCORRECT_YIELD_DIVS
#ifdef CORRECT_YIELD_DIVS
    DividendListConstSP eqDivs(asset->equity->getDivList());
    DividendListSP divs(eqDivs->getAllDivsBetweenDates(
        collector->getStartDate(), collector->getEndDate()));
    if (collector->requireDollarDivs()){
        //CorrelationSP corrBak = asset->corrEqFX;
        //asset->corrEqFX = CorrelationSP(new Correlation("", "", "", 0.0));
        divs->convertYieldDivs(asset);
        //asset->corrEqFX = corrBak;
    }
    const DividendArray& divArray = divs->getArray();
    for (int i = 0; i < divArray.size(); i++){
        collector->addDividend(divArray[i]);
    }
#else
    asset->equity->accept((ICollector*)collector);
#endif
}


/** given a current spot level, get the next strike on the vol surface where 
    the slope is non-differentiable */
double ProtEquity::getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const
{
    double volStrike;
    const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(vol.get());
    if ( nextStrike )
    {
        volStrike = nextStrike->getNextStrike(strike,
                                              isUp,
                                              offSurface);
    }
    else
    {
        volStrike  = 0.0;
        offSurface = true;
    }
    return volStrike;
}

/** Returns the name of the stock/asset - used to determine
    whether to shift the object */
string ProtEquity::sensName(SpotLevelProbability* shift) const {
    return equity->getName();
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
bool ProtEquity::sensShift(SpotLevelProbability* shift) {
    double spot = shift->spotLevel(baseDate, this);
    SpotLevel level(spot);
    equity->sensShift(&level);  // set equity spot to our new level
    return false;  // all done;
}

/** returns dividend list */
DividendListConstSP  ProtEquity::getDivList() const {
    return equity->getDivList();
}

PDFCalculator* ProtEquity::pdfCalculator(const PDFRequest* request) const {
    if (IPDFCalculator::TYPE->isInstance(vol.get())) {
        const IPDFCalculator* pdf = dynamic_cast< const IPDFCalculator*>(vol.get());
        return pdf->getPDFCalculator(request, this);
    }
    throw ModelException("ProtEquity::pdfCalculator",
                         "vol of type (" + vol->getClass()->getName() +
                         ") has no pdf calculator");
}

/* for IEqVolNamePair */
bool ProtEquity::getNamePairs(string& eqName, string& volName) const {
    eqName = getTrueName();
    volName = getVolName();
    return false;   //no more assets inside
}

//// Returns the ccy treatment for the asset. Default returns 
string ProtEquity::getCcyTreatment() const{
    return CCY_TREATMENT_PROTECTED;
}

/** Can this asset physically settle? */
bool ProtEquity::canPhysicallySettle() const {
    return false;
}

// the IMarketObservable interface for retrieving a single sample
double ProtEquity::pastValue(const DateTime&               sampleDate,
                               const ObservationType*      obsType,
                               const ObservationSource*    source,
                               const FixingType*           fixType,
                               const IObservationOverride* overrides,
                               const SamplingConvention*   sampleRule) const{
    return equity->pastValue(sampleDate, obsType, source, fixType,
                             overrides, sampleRule);
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool ProtEquity::observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const {
    return equity->observationDate(sampleDate, source, sampleRule, obsDate);
}

// the IMarketObservable interface for retrieving past samples events
double ProtEquity::addPastSampleEvent(const DateTime&       sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const {
    return equity->addPastSampleEvent(sampleDate, obsType, source, fixType,
                                      overrides, sampleRule, collector);
}

// the IMarketObservable interface for 
// is the given date a holiday for the relevant source
bool ProtEquity::isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const{
    return equity->isHoliday(sampleDate, source);
}

static CFieldConstSP fxVolField;

bool ProtEquity::recurse(const CFieldConstSP& field,
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
            if (targetClass == ITweakableWithRespectTo<VolParallel>::TYPE   ||
                targetClass == IRestorableWithRespectTo<VolParallel>::TYPE  ||
                targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE  ||
                targetClass == IRestorableWithRespectTo<VolPointwise>::TYPE ||
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
ProtEquity::ProtEquity(): CAsset(TYPE), coherentGreeks(false) {}

class ProtEquityHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ProtEquity, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(Asset::IQuanto);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(INextStrike);
        IMPLEMENTS(SpotLevelProbability::Shift);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(IEqVolNamePair);
        IMPLEMENTS(ObjectIteration::IOverride);
        FIELD(name, "Asset's name");
        FIELD_MAKE_OPTIONAL(name);
        EMPTY_SHELL_METHOD(defaultProtEquity);
        FIELD(vol,         "Volatility");
        FIELD(equity,             "Stock");
        FIELD(fxVol,       "FX Volatility");
        FIELD(baseDate,    "Value Date");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(corrEqFX,           "Equity FX Correlation");
        FIELD_MAKE_OPTIONAL(corrEqFX);
        FIELD(protCcyCode, "Protected Ccy code");
        FIELD(getCorrFromCache, "Cached flag");
        FIELD_MAKE_TRANSIENT(getCorrFromCache); // hide from dd interface
        FIELD(scale, "Cached ccy adjustment");
        FIELD_MAKE_TRANSIENT(scale);            // hide from dd interface
        FIELD(protCcyIsoCode, "Protected ISO Ccy code");
        FIELD_MAKE_TRANSIENT(protCcyIsoCode);  // hide from dd interface
        FIELD(coherentGreeks, "Use homogeneous greeks for FX");
        FIELD_MAKE_TRANSIENT(coherentGreeks);
        FIELD(fx, "");
        FIELD_MAKE_TRANSIENT(fx);

        ClassSetAcceptMethod(ProtEquity::acceptValueDateCollector);
        ClassSetAcceptMethod(ProtEquity::acceptDeltaShift);
        ClassSetAcceptMethod(ProtEquity::acceptImntCcy);
        ClassSetAcceptMethod(ProtEquity::acceptHoliday);
        ClassSetAcceptMethod(ProtEquity::acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);

        // look up field for use on recurse
        fxVolField = clazz->getDeclaredField("fxVol");
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const ProtEquity*  asset,
                                        DividendCollector* collector){
        collector->addDivs(asset->equity->getDivList(), 0 /* not struck */);
    }

    static IObject* defaultProtEquity(){
        return new ProtEquity();
    }
};

CClassConstSP const ProtEquity::TYPE = CClass::registerClassLoadMethod(
    "ProtEquity", typeid(ProtEquity), ProtEquityHelper::load);

DRLIB_END_NAMESPACE
