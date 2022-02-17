//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CVolProcessedBS.cpp
//
//   Description : What a VolatilityBS can do
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLPROCESSEDBS_CPP
#include <algorithm>
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/Addin.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/IPDFBoundaryProb.hpp"

DRLIB_BEGIN_NAMESPACE

const double CVolProcessedBS::LEGACY_PDF_BOUNDARY_PROB = 0.000000000001;

/** Calculates variance between a series of dates. If the dateList
        has n dates in it, n-1 vols will be calculated. */
void CVolProcessedBS::CalcVar(const DateTimeArray& dateList,
                              TCalcType            calcType, 
                              CDoubleArray&        vars) const
{
    int numDates = dateList.size();
    static char routine[] = "CVolProcessedBS::CalcVar";
    if (numDates < 2 ){
        throw ModelException(routine, "Must supply at least two dates");
    }
    if (vars.size() < numDates-1){
        throw ModelException(routine, "Double array too short");
    }
    for (int i = 0; i < numDates-1; i++){
        vars[i] = CalcVar(dateList[i], dateList[i+1]);
    }
    switch (calcType){
        int i; // MSVC is broken - doesn't support separate variables in loop
    case forward:
        // do nothing
        break;
    case fromFirst:
        for (i = 1; i < numDates-1; i++){
            vars[i] += vars[i-1];
        }
        break;
    case toLast:
        for (i = numDates-2; i >= 0; i--){
            vars[i] += vars[i+1];
        }
        break;
    default:
        throw ModelException(routine, "Unknown calculate type");
    }
    return;
}

/** Calculates variance beginning at dateFrom. If the dateList
    has n dates in it, n variances will be calculated. */
void CVolProcessedBS::CalcVar(const DateTime &dateFrom,
                              const DateTimeArray& datesTo,
                              TCalcType            calcType, 
                              CDoubleArray&        vols) const
{
    int numDates = datesTo.size();
    DateTimeArray dateSeries(datesTo.size()+1);

    dateSeries[0] = dateFrom;
    for(int i=0;i<numDates;++i)
    {
        dateSeries[i+1] = datesTo[i];
    }
    CalcVar(dateSeries, calcType, vols);
}


/** Calculates vols between a series of dates */
void CVolProcessedBS::CalcVol(const DateTimeArray& dateList, 
                              TCalcType            calcType, 
                              CDoubleArray&        vols) const
{
    int numDates = dateList.size();
    static char routine[] = "CVolProcessedBS::CalcVol";
    if (numDates < 2 && vols.size() > 0 ) {
        throw ModelException(routine, "Must supply at least two dates");
    }
    if (vols.size() < numDates-1){
        throw ModelException(routine, "Double array too short");
    }
    switch (calcType){
        int i; // MSVC is broken - doesn't support separate variables in loop
    case forward:
        for (i = 0; i < numDates-1; i++){
            vols[i] = CalcVol(dateList[i], dateList[i+1]);
        }
        break;
    case fromFirst:
        // probably not the most efficient
        for (i = 1; i < numDates; i++){
            vols[i-1] = CalcVol(dateList[0], dateList[i]);
        }
        break;
    case toLast:
        // probably not the most efficient
        for (i = 0; i < numDates-1; i++){
            vols[i] = CalcVol(dateList[i], dateList[numDates-1]);
        }
        break;
    default:
        throw ModelException(routine, "Unknown calculate type");
    }
    return;
}

/** Calculates variance beginning at dateFrom. If the dateList
    has n dates in it, n variances will be calculated. */
void CVolProcessedBS::CalcVol(const DateTime &dateFrom,
                              const DateTimeArray& datesTo,
                              TCalcType            calcType, 
                              CDoubleArray&        vols) const
{
    int numDates = datesTo.size();
    DateTimeArray dateSeries(datesTo.size()+1);

    dateSeries[0] = dateFrom;
    for(int i=0;i<numDates;++i)
    {
        dateSeries[i+1] = datesTo[i];
    }
    CalcVol(dateSeries, calcType, vols);
}

/** utility method - interpolates vol then calculates vols */
void CVolProcessedBS::calcVol(const CVolBase*       vol,
                              const CVolRequest*    request,
                              const CAsset*         asset,
                              const DateTime&       dateFrom,
                              const DateTimeArray&  datesTo,
                              CDoubleArray&         vols){
    CVolProcessedSP interpVol(vol->getProcessedVol(request, asset));
    CVolProcessedBS& interped = dynamic_cast<CVolProcessedBS&>(*interpVol);
    interped.CalcVol(dateFrom, datesTo, fromFirst, vols);
}

/** Generates an array of 'default' strikes given an asset. The
    strikes are chosen such that they are evenly distributed in
    probability space based upon the at-the-money 1 year vol. Strikes
    are also generated for probabilities "0+epsilon" and "1-epsilon". The
    maxNumBands is an indicative indicator for the number of strikes to
    generate. pdfBoundaryProb determines the range of strikes to return
    by limiting the probability distribution from above. */
DoubleArray CVolProcessedBS::defaultStrikes(int             maxNumBands,
                                            const DateTime& baseDate,
                                            const CAsset*   asset,
                                            const IPDFBoundaryProb* marketObj){


    double pdfBoundaryProb = marketObj->getPDFBoundaryProb();
    double spot = asset->getSpot(); // for ease
    // Then somewhat arbitrarily use 1yr atm for our distribution
    DateTime oneYear(baseDate.rollDate(365));
    double oneYearFwd = asset->fwdValue(oneYear); // one year forward
    ATMVolRequest atmVolRequest;
    CVolProcessedBSSP vol(asset->getProcessedVol(&atmVolRequest));
    double oneYearVol = vol->CalcVol(baseDate, oneYear);
    // for familiarity map to usual greek names
    double sigma = oneYearVol; // since time in years = 1
    double mu = log(oneYearFwd/spot) - 0.5 * oneYearVol * oneYearVol;
    // number of bands: eg 20 => every 5%
    const int numBands = Maths::max(maxNumBands, 2); 
    //const double boundaryProb = 0.000000000001; // fairly arbitrary, now replaced by customisable pdfBoundaryProb
    DoubleArray strikes(Maths::isZero(oneYearVol)? 
                        1: (numBands+1)); // includes 'end' points
    // loop through probabilities
    for (int i = 0; i < strikes.size(); i++){
        double prob;
        if (i == 0){
            prob = pdfBoundaryProb;
        } else if (i == numBands){
            prob = 1.0 - pdfBoundaryProb; // avoid infinite strike
        } else {
            prob = (double)i/numBands;
        }
        double logStrikeRatio = sigma*N1Inverse(prob) + mu;
        strikes[i] = spot * exp(logStrikeRatio);
    }
    return strikes;
}

/** As above but uses maxNumBands = 20 ie every 5% */
DoubleArray CVolProcessedBS::defaultStrikes(const DateTime& baseDate,
                                            const CAsset*   asset,
                                            const IPDFBoundaryProb* marketObj){
    const int DEFAULT_NUM_BANDS = 20; /* fairly arbitrary (but how strikes
                                         would a trader use?) */
    return defaultStrikes(DEFAULT_NUM_BANDS, baseDate, asset, marketObj);
}

static void volProcessedBSLoad(CClassSP& clazz){
    REGISTER(CVolProcessedBS, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
}

CVolProcessedBS::CVolProcessedBS(const CClassConstSP& clazz): CObject(clazz){}

CClassConstSP const CVolProcessedBS::TYPE = CClass::registerClassLoadMethod(
    "VolProcessedBS", typeid(CVolProcessedBS), volProcessedBSLoad);

/** Wrapper around 
    VolProcessedBSIR::sensitiveIRVolPoints(const DateTimeArray&) */
IRGridPointAbsArraySP VolProcessedBSIR::sensitiveIRVolPoints(
    const IRVolBase*     irVol,
    const YieldCurve*    yc,
    CVolRequest*         volRequest,
    const DateTimeArray& dates){
    CVolProcessedSP processedVol(irVol->getProcessedVol(volRequest, yc));
    if (!VolProcessedBSIR::TYPE->isInstance(processedVol)){
        throw ModelException("VolProcessedBSIR::sensitiveIRVolPoints",
                             "Processed vol of type "+
                             processedVol->getClass()->getName()+
                             " not supported");
    }
    VolProcessedBSIR* bsir = DYNAMIC_CAST(VolProcessedBSIR, processedVol.get());
    return bsir->sensitiveIRVolPoints(dates);
}

/** for VolProcessedBSIR class */
static void volProcessedBSIRLoad(CClassSP& clazz){
    REGISTER(VolProcessedBSIR, clazz);
    SUPERCLASS(CVolProcessedBS);
}

VolProcessedBSIR::VolProcessedBSIR(const CClassConstSP& clazz): 
    CVolProcessedBS(clazz){}

CClassConstSP const VolProcessedBSIR::TYPE = CClass::registerClassLoadMethod(
    "VolProcessedBSIR", typeid(VolProcessedBSIR), volProcessedBSIRLoad);

class LNVolAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    DateTime           baseDate; 
    DateTimeArraySP    dates;     
    CVolProcessedBSSP  procVol;

    // return the vols from the given date to each future date
    static IObjectSP calcVolsLN(LNVolAddin* params){
        static const string method = "VolProcessedBS::calcVolsLN";

        // check that the dates are in ascending order
        int numDates = params->dates->size();

        for (int i = 1; i < numDates; i++)
        {
            if ((*params->dates)[i-1].isGreaterOrEqual((*params->dates)[i]))
            {
                throw ModelException(method,
                                     "dates are not in ascending order. (" +
                                     (*params->dates)[i-1].toString() +
                                     ") comes before (" +
                                     (*params->dates)[i].toString() +
                                     ") in the dates array.");
            } 
        }

        // if the base date is before the first array date then 
        // insert the base date at the front of the dates array
        if (params->baseDate.isGreaterOrEqual((*params->dates)[0]))
        {
            throw ModelException(method, 
                                 "base date (" +
                                 params->baseDate.toString() +
                                 ") is not before first array date (" +
                                 (*params->dates)[0].toString() +
                                 ")");
        }

        vector<DateTime>::iterator iter = params->dates->begin();
        params->dates->insert(iter, params->baseDate);

        // create the output double array and get the vols
        CDoubleArraySP vols(new CDoubleArray(params->dates->size() -1));
        CVolProcessedBS::TCalcType calcType = CVolProcessedBS::fromFirst;
        params->procVol->CalcVol((*params->dates),
                                 calcType,
                                 *vols);
                                 
        return vols;
    }
    

    LNVolAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LNVolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLNVolAddin);
        FIELD(baseDate, "date from");
        FIELD(dates, "dates to");
        FIELD(procVol, "processed vol");
        Addin::registerClassObjectMethod("LN_VOL",
                                         Addin::MARKET,
                                         "Calculates vols from a given date "
                                         "to a set of future dates.",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)calcVolsLN);
    }
    
    static IObject* defaultLNVolAddin(){
        return new LNVolAddin();
    }
 
};

CClassConstSP const LNVolAddin::TYPE = CClass::registerClassLoadMethod(
    "LNVolAddin", typeid(LNVolAddin), load);

// make DummyAsset available to LN vol interp addins, but keep
// the namespace clear for any other DummyAsset that may arise
class LNVolInterpHelper {
public:
    /** This is the safest way to create an asset that should not be used 
        (better than passing null). Might want to move this somewhere more
        accessible to other similar addin functions */
    class DummyAsset: public CAsset{
    private:
        ModelException fail() const{
            return ModelException("DummyAsset", "Asset information not"
                                  " available");
        }
        static void load(CClassSP& clazz){
            // NB not public
            REGISTER(DummyAsset, clazz);
            SUPERCLASS(CAsset);
            EMPTY_SHELL_METHOD(defaultDummyAsset);
        }
        static IObject* defaultDummyAsset(){
            return new DummyAsset();
        }
        
    public:
        static CClassConstSP const TYPE;
        DummyAsset(): CAsset(TYPE){}
        virtual double getSpot() const{
            return 0.0; /* unfortunate but eg for parameterised vols there is a
                           call to getSensitiveStrike which needs the spot -
                           even if it's not used */
        }
        virtual string getName() const{ throw fail(); }
        virtual CVolProcessed* getProcessedVol(
            const CVolRequest* volRequest) const{ throw fail(); }
        virtual double fwdValue(
            const DateTime& date) const { throw fail(); }
        virtual void fwdValue(
            const DateTimeArray&             dateList,
            const CAsset::FwdValueAlgorithm& algo,
            CDoubleArray&                    result) const{ throw fail(); }
        virtual string getYCName() const{ throw fail(); }
        virtual DateTime settleDate(
            const DateTime& tradeDate) const{ throw fail(); }
        virtual PDFCalculator* pdfCalculator(
            const PDFRequest* request) const{ throw fail(); }
        // the IMarketObservable interface for retrieving a sampling event
        double addPastSampleEvent(const DateTime&     sampleDate,
                                  const ObservationType*      obsType,
                                  const ObservationSource*    source,
                                  const FixingType*           fixType,
                                  const IObservationOverride* overrides,
                                  const SamplingConvention*   sampleRule,
                                  PastSamplesCollector*        collector) const {
            throw fail();
        }
        // the IMarketObservable interface for retrieving a single sample
        virtual double pastValue(const DateTime&            sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule) const{
            throw fail();
        }

        // IMarketObservable - retrieve a single observation date
        // Returns false if obs is to be omitted
        virtual bool observationDate(const DateTime&           sampleDate,
                                     const ObservationSource*  source,
                                     const SamplingConvention* sampleRule,
                                     DateTime*                 obsDate) const {
            throw fail();
        }

        // the IMarketObservable interface for 
        // is the given date a holiday for the relevant source
        virtual bool isHoliday(const DateTime& sampleDate,
                               const ObservationSource*   source) const{
            throw fail();
        }
    };
};    

CClassConstSP const LNVolInterpHelper::DummyAsset::TYPE = 
CClass::registerClassLoadMethod(
    "LNVolInterpHelper::DummyAsset", typeid(LNVolInterpHelper::DummyAsset), 
    load);

// not another vol interp addin you cry ! Well this one is fully market based
// and driven off the vol name & type
class LNMktVolInterpAddin: public CObject, public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    MarketDataSP       market;
    string             name;
    double             strike;
    DateTimeArray      toDates;
    string             volType;  // optional
    DateTime           fromDate; // optional

    // EdrAction version of addin
    IObjectSP run() {
        return volInterp();
    }

    // return the vols from the given date to each future date
    DoubleArraySP volInterp(){
        static const string method("LNMktVolInterpAddin::volInterp");
        try {
            // finds out where we want vols from
            DateTime fromDate;
            if (this->fromDate.empty()) {
                market->GetReferenceDate(fromDate);
            }
            else {
                fromDate = this->fromDate;
            }
              
            // need a model to get data out of the market
            MarketDataFetcherLNSP mdf(
                new MarketDataFetcherLN(volType.empty() ?
                                        "VolPreferred" : volType));
            IModelSP model(new NonPricingModel(mdf));;

            MarketObjectSP volObj(market->GetData(model.get(),
                                                  name,
                                                  CVolBase::TYPE));
            
            CVolBase* vol = dynamic_cast<CVolBase*>(volObj.get());

            // create a simple vol request
            LinearStrikeVolRequest volRequest(strike, 
                                              DateTime(), DateTime(),
                                              false);
            // create a dummy asset
            LNVolInterpHelper::DummyAsset dummyAsset;
            // process the vol
            CVolProcessedSP procVol(vol->getProcessedVol(&volRequest,
                                                         &dummyAsset));
            CVolProcessedBS& volLN = dynamic_cast<CVolProcessedBS&>(*procVol);
            // then interpolate
            DoubleArraySP theVols(new DoubleArray(toDates.size()));
            volLN.CalcVol(fromDate, 
                          toDates,
                          CVolProcessedBS::fromFirst,
                          *theVols);
            return theVols;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }            
    }
    

    LNMktVolInterpAddin(): CObject(TYPE), strike(0.0), volType(""){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(LNMktVolInterpAddin, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultLNMktVolInterpAddin);
        FIELD(market, "market data cache");
        FIELD(name, "vol name");
        FIELD(strike, "strike");
        FIELD(toDates, "dates for vol");
        FIELD(volType, "vol type");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(fromDate, "date from");
        FIELD_MAKE_OPTIONAL(fromDate);
        Addin::registerObjectMethod("LN_VOL_INTERP",
                                    Addin::MARKET,
                                    "Interp vol for a given strike",
                                    false,
                                    Addin::expandSimple,
                                    &LNMktVolInterpAddin::volInterp);
    }
    
    static IObject* defaultLNMktVolInterpAddin(){
        return new LNMktVolInterpAddin();
    }
 
};

CClassConstSP const LNMktVolInterpAddin::TYPE=CClass::registerClassLoadMethod(
    "LNMktVolInterpAddin", typeid(LNMktVolInterpAddin), load);


class LNInterpVolAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    CVolBaseWrapper    vol;
    double             absStrike;
    DateTime           fromDate; 
    DateTime           toDate;     
    MarketDataSP       market; // optional
    string             volType; // optional

    // return the vols from the given date to each future date
    double calcVolsLN(){
        static const string method = "LNInterpVolAddin::calcVolsLN";
        // if market supplied use it
        if (market.get()){
            IModelSP npm(new NonPricingModel());
            // if vol type supplied use it
            if (!volType.empty()){
                MarketDataFetcherLNSP mdf(
                    new MarketDataFetcherLN(volType));
                npm = IModelSP(new NonPricingModel(mdf));
            }
            vol.getData(npm.get(), market.get());
        }
        // create a simple vol request
        LinearStrikeVolRequest volRequest(absStrike, 
                                          DateTime(), DateTime(),
                                          false);
        // create a dummy asset
        LNVolInterpHelper::DummyAsset dummyAsset;
        // process the vol
        CVolProcessedSP procVol(vol->getProcessedVol(&volRequest,
                                                             &dummyAsset));
        CVolProcessedBS& procVolLN = dynamic_cast<CVolProcessedBS&>(*procVol);
        // then interpolate
        return (procVolLN.CalcVol(fromDate, toDate));
    }
    

    LNInterpVolAddin():
    CObject(TYPE),
    volType(""){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LNInterpVolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLNInterpVolAddin);
        FIELD(vol, "vol supporting LN interpolation");
        FIELD(absStrike, "Absolute strike");
        FIELD(fromDate, "date from");
        FIELD(toDate, "date to");
        FIELD(market, "optional market data cache");
        FIELD_MAKE_OPTIONAL(market);
        FIELD(volType, "vol type");
        FIELD_MAKE_OPTIONAL(volType);
        Addin::registerDoubleMethod("LN_INTERP_VOL",
                                    Addin::MARKET,
                                    "Calculates vol between two dates",
                                    &LNInterpVolAddin::calcVolsLN);
    }
    
    static IObject* defaultLNInterpVolAddin(){
        return new LNInterpVolAddin();
    }
 
};

CClassConstSP const LNInterpVolAddin::TYPE = CClass::registerClassLoadMethod(
    "LNInterpVolAddin", typeid(LNInterpVolAddin), load);

class LNVarAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    DateTime           baseDate; 
    DateTimeArraySP    dates;     
    CVolProcessedBSSP  procVol;


    // return the variance from the given date to each future date
    static IObjectSP calcVarsLN(LNVarAddin* params){
        static const string method = "VolProcessedBS::calcVarsLN";

        // check that the dates are in ascending order
        int numDates = params->dates->size();

        for (int i = 1; i < numDates; i++)
        {
            if ((*params->dates)[i-1].isGreaterOrEqual((*params->dates)[i]))
            {
                throw ModelException(method,
                                     "dates are not in ascending order. (" +
                                     (*params->dates)[i-1].toString() +
                                     ") comes before (" +
                                     (*params->dates)[i].toString() +
                                     ") in the dates array.");
            } 
        }

        // if the base date is before the first array date then 
        // insert the base date at the front of the dates array
        if (params->baseDate.isGreaterOrEqual((*params->dates)[0]))
        {
            throw ModelException(method, 
                                 "base date (" +
                                 params->baseDate.toString() +
                                 ") is not before first array date (" +
                                 (*params->dates)[0].toString() +
                                 ")");
        }

        vector<DateTime>::iterator iter = params->dates->begin();
        params->dates->insert(iter, params->baseDate);

        // create the output double array and get the vols
        CDoubleArraySP vars(new CDoubleArray(params->dates->size() -1));
        CVolProcessedBS::TCalcType calcType = CVolProcessedBS::fromFirst;
        params->procVol->CalcVar((*params->dates),
                                 calcType,
                                 *vars);
                                 
        return vars;
    }
    

    LNVarAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LNVarAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLNVarAddin);
        FIELD(baseDate, "date from");
        FIELD(dates, "dates to");
        FIELD(procVol, "processed vol");
        Addin::registerClassObjectMethod("LN_VAR",
                                         Addin::MARKET,
                                         "Calculates variance from a given date "
                                         "to a set of future dates.",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)calcVarsLN);
    }
    
    static IObject* defaultLNVarAddin(){
        return new LNVarAddin();
    }
 
};

CClassConstSP const LNVarAddin::TYPE = CClass::registerClassLoadMethod(
    "LNVarAddin", typeid(LNVarAddin), load);

/** Vol Surface generator addin */
class VSGAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    CVolBaseWrapper    vol;
    DoubleArray        strikes;
    ExpiryArray        benchmarks;
    DateTime           today;
    string             volType;
    CMarketDataSP      market;
    bool               negativeFwdVarAllowed;
             
    // Return a double matrix at requested strikes and maturities
    static IObjectSP calcVolSurf(VSGAddin* params){
        return params->calcSurf();
    }

    IObjectSP calcSurf(){
        static const string method = "VSGAddin::calcSurf";
        // use valueDate in market if it's there
        market->GetReferenceDate(today);
        if (!vol){
            // Ideally like to instantiate a Model that will select what
            // type of vol to get - however not in right directory.
            // Need to manually pull the vol out of the cache since we would
            // need the model to specify the type of vol to get and there
            // may be many BS vols in the cache
            CClassConstSP volClass = CClass::forName(volType);
            vol.setObject(market->GetData(vol.getName(), volClass));
        }
        // create a non pricing model in order to get the market data
        NonPricingModel dummyModel;
        vol->getMarket(&dummyModel, market.get());
        // convert expiries into dates
        DateTimeArray  dates(benchmarks.size());
        for (int i = 0; i < dates.size(); i++){
            dates[i] = benchmarks[i]->toDate(today);
        }
        // check dates are increasing and not empty
        DateTime::ensureIncreasing(dates, "Benchmark dates", true);
        if (strikes.empty()){
            throw ModelException(method, "No strikes specified");
        }
        // create space for returned object
        CDoubleMatrixSP output(new DoubleMatrix(strikes.size(), dates.size()));
        DoubleMatrix&  matrix = *output;
        // scratchpad for doing calculations
        DoubleArray    vols(dates.size());
        // then loop across strikes
        for (int j = 0; j < strikes.size(); j++){
            // create a vol request - v simple one
            LinearStrikeVolRequest request(strikes[j], today, today, false);
            request.allowNegativeFwdVar(negativeFwdVarAllowed);
            // use this lovely little utitility function
            CVolProcessedBS::calcVol(vol.get(),
                                     &request,
                                     0, // asset
                                     today,
                                     dates,
                                     vols);
            // copy over vols to matrix
            for (int k = 0; k < dates.size(); k++){
                matrix[j][k] = vols[k];
            }
        }
        return output;
    }
    

    VSGAddin(): CObject(TYPE), volType(IVolatilityBS::TYPE->getName()),
                market(new MarketData()), negativeFwdVarAllowed(false) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VSGAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVSGAddin);
        FIELD(vol, "Volatility to generate surface from");
        FIELD(strikes, "Strikes for surface");
        FIELD(benchmarks, "Benchmark dates for surface");
        FIELD(today, "Base date for surface");
        FIELD_MAKE_OPTIONAL(today);
        FIELD(volType, "Type of the vol in the market data cache");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(market, "Market data cache");
        FIELD_MAKE_OPTIONAL(market);
        FIELD(negativeFwdVarAllowed, "Allow negative fwd variance");
        FIELD_MAKE_OPTIONAL(negativeFwdVarAllowed);
        Addin::registerClassObjectMethod("VOL_SURFACE_GENERATOR",
                                         Addin::MARKET,
                                         "Calculates vol surface from "
                                         "supplied vol",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)calcVolSurf);
    }
    
    static IObject* defaultVSGAddin(){
        return new VSGAddin();
    }
 
};

CClassConstSP const VSGAddin::TYPE = CClass::registerClassLoadMethod(
    "VSGAddin", typeid(VSGAddin), load);

/** Vol Surface generator addin using asset rather than vol */
class VSG2Addin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    CAssetWrapper       asset;
    DoubleArray         strikes;
    ExpiryArray         benchmarks;
    DateTime            today;
    DateTime            startDate;
    string              volType;
    CMarketDataSP       market;
    bool                negativeFwdVarAllowed;
             
    // Return a double matrix at requested strikes and maturities
    static IObjectSP calcVolSurf(VSG2Addin* params){

        // create a non pricing model in order to get the market data
        MarketDataFetcherSP mdf(new MarketDataFetcherLN(params->volType));
        NonPricingModel dummyModel(mdf);
        params->asset.getData(&dummyModel, params->market.get());

        return calcSurf(params);
    }

    static IObjectSP calcSurf(VSG2Addin* params){
        static const string method = "VSG2Addin::calcSurf";
        if (!(params->asset)) {
            // Asset should have been previously populated with market data
            throw ModelException(method, "Internal error: Asset is not populated with market data");
        } else {
            return params->calcSurf();
        }
    }
    
    IObjectSP calcSurf() {
        static const string method = "VSG2Addin::calcSurf";
        
        // use valueDate in market if it's there
        market->GetReferenceDate(today);
        bool fwdStart;
        if (startDate.empty()){
            fwdStart = false;
        } else {
           fwdStart = startDate.isGreater(today);
        }
        if (!fwdStart){
            startDate = today;
        }
     
        // convert expiries into dates
        DateTimeArray  dates(benchmarks.size());
        for (int i = 0; i < dates.size(); i++){
            dates[i] = benchmarks[i]->toDate(startDate);
        }
        // check dates are increasing and not empty
        DateTime::ensureIncreasing(dates, "Benchmark dates", true);
        if (strikes.empty()){
            throw ModelException(method, "No strikes specified");
        }
        // create space for returned object
        CDoubleMatrixSP output(new DoubleMatrix(strikes.size(), dates.size()));
        DoubleMatrix&  matrix = *output;
        // scratchpad for doing calculations
        DoubleArray    vols(dates.size());
        // then loop across strikes
        for (int j = 0; j < strikes.size(); j++){
            // create a vol request - v simple one
            LinearStrikeVolRequest request(strikes[j],
                                           startDate, 
                                           dates.back(), 
                                           fwdStart);
            request.allowNegativeFwdVar(negativeFwdVarAllowed);
            CVolProcessedBSSP vol(asset->getProcessedVol(&request));
            vol->CalcVol(startDate, dates, CVolProcessedBS::fromFirst, vols);
            // copy over vols to matrix
            for (int k = 0; k < dates.size(); k++){
                matrix[j][k] = vols[k];
            }
        }
        return output;
    }
    

    VSG2Addin(): CObject(TYPE), volType(IVolatilityBS::TYPE->getName()),
                market(new MarketData()), negativeFwdVarAllowed(false) {}

    VSG2Addin(CAssetWrapper       asset, 
              DoubleArray         strikes,
              ExpiryArray         benchmarks,
              DateTime            today,
              DateTime            startDate,
              string              volType,
              CMarketDataSP       market,
              bool                negativeFwdVarAllowed): CObject(TYPE),
        asset(asset), strikes(strikes), benchmarks(benchmarks), today(today),
        startDate(startDate), volType(volType), market(market), 
        negativeFwdVarAllowed(negativeFwdVarAllowed)
    {};
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VSG2Addin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVSG2Addin);
        FIELD(asset, "Asset to use to build vol surface");
        FIELD(strikes, "Strikes for surface");
        FIELD(benchmarks, "Benchmark dates for surface"
                     " relative to start date");
        FIELD(startDate, "When to calculate vols from");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(today, "Today");
        FIELD_MAKE_OPTIONAL(today);
        FIELD(volType, "Type of the vol in the market data cache");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(market, "Market data cache");
        FIELD_MAKE_OPTIONAL(market);
        FIELD(negativeFwdVarAllowed, "Allow negative fwd variance");
        FIELD_MAKE_OPTIONAL(negativeFwdVarAllowed);
        Addin::registerClassObjectMethod("VOL_SURFACE_GENERATOR2",
                                         Addin::MARKET,
                                         "Calculates vol surface from "
                                         "supplied asset",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)calcVolSurf);
    }
    
    static IObject* defaultVSG2Addin(){
        return new VSG2Addin();
    }
 
};

CClassConstSP const VSG2Addin::TYPE = CClass::registerClassLoadMethod(
    "VSG2Addin", typeid(VSG2Addin), load);

typedef smartPtr<VSG2Addin> VSG2AddinSP;

class ImpliedVolSurface: public CObject, virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // Wrapper around VSG2Addin to imply a vol surface
    // for general purposes - e.g. building a fund vol surface or
    // implying a strike-based surface from an asset with a
    // delta-based surface
    ExpiryArraySP                   benchmarks;
    DoubleArray                     strikes;
    CAssetWrapper                   asset;
    CMarketDataSP                   market;
    string                          volType;
    bool                            strikesInPct;
    string                          ccyTreatment;
    YieldCurveWrapper               discount; // needed if struck
    bool                            includeSurfaceStrikes;
    DoubleArraySP                   excludeStrikes;

    // EdrAction version of addin
    IObjectSP run() {
        return ImpliedVolSurfaceAddin(this);
    }

    static IObjectSP ImpliedVolSurfaceAddin(ImpliedVolSurface* params){
        static const string routine = "ImpliedVolSurface::ImpliedVolSurfaceAddin";
        try {
            return params->implyVols();
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** the 'addin function' */
    // Return a DoubleMatrix if includeSurfaceStrikes=false and there are no excludeStrikes
    // Otherwise return a VolSurface (so that the strikes are explicit).
    IObjectSP implyVols() {
        static const string method = "ImpliedVolSurface::implyVols";
        try {
            MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
            NonPricingModel dummyModel(mdf);

            if (ccyTreatment==CAsset::CCY_TREATMENT_NONE || 
                ccyTreatment==CAsset::CCY_TREATMENT_VANILLA) {
               // No discount curve needed
                asset.getData(&dummyModel, market.get());
            } else {
                if (discount.isEmpty()) {
                    throw ModelException(method, "Discount curve must be supplied if asset is struck or protected");
                } else {
                    discount.getData(&dummyModel, market);
                    CAsset::getAssetMarketData(&dummyModel, market.get(), ccyTreatment, 
                                            discount, asset);
                }
            }
            
            VolSurfaceSP oldSurface;
            // Whether to return a surface or matrix
            bool returnSurface =
                (includeSurfaceStrikes || (!!excludeStrikes && excludeStrikes->size()>0));
            if  (returnSurface) {
                // Extract the volSurface from the market data
                double arbitraryShiftSize = 0.0001;
                SensControlPerNameSP shift(new VegaMatrix(arbitraryShiftSize));
                OutputNameArrayConstSP sensNames(shift->names(asset.getMO().get()));
                if (sensNames->size() > 1) {
                    throw ModelException(method, "Found more than one vol surface for asset " + asset.getName());
                } else if (sensNames->size() == 0) {
                    throw ModelException(method, "No vol surface found for asset " + asset.getName());
                }
                MarketObjectSP getVolSurface;
                try {
                    getVolSurface = market->GetData((*sensNames)[0]->idGet(0), VolSurface::TYPE);
                } catch (exception& e) {
                   throw ModelException(e, method, "Failed to retrieve strikes from vol surface for "+
                        (*sensNames)[0]->idGet(0) + ".");
                }
                oldSurface = VolSurfaceSP::dynamicCast(getVolSurface);
                oldSurface->getMarket(&dummyModel, market.get());
            }

            // Get strikes for which vols are to be retrieved
            DoubleArraySP volStrikes(gatherStrikes(oldSurface));

            DateTime today;
            market->GetReferenceDate(today);

            VSG2AddinSP myParams(new VSG2Addin(asset,     
                                               *volStrikes,
                                               *benchmarks,
                                               today,     // today
                                               today,     // startDate
                                               volType,
                                               market,
                                               false));        // -ve FwdVarAllowed
            
            IObjectSP matrix = VSG2Addin::calcSurf(myParams.get());

            if (returnSurface) {
                VolSurfaceSP surface(new VolSurface(oldSurface->getName(),
                                                    oldSurface->getTimeMetric().get(),
                                                    *volStrikes,
                                                    *(dynamic_cast<DoubleMatrix *>(matrix.get())),
                                                    benchmarks.get(),
                                                    oldSurface->getBaseDate()));
                return surface;
            }
            else {
                return matrix;
            }
        } 
        catch (exception& e){
            throw ModelException(e, method, 
                                 "Failed for asset " + asset->getName());
        }
    }

private:

    DoubleArraySP gatherStrikes(VolSurfaceSP surface) {
        static const string method = "ImpliedVolSurface::gatherStrikes";
        try {
            int i;
            DoubleArraySP volStrikes(new DoubleArray(0));

            // Extract the strikes from the surface
            if (includeSurfaceStrikes && !!surface) {
                DoubleArray surfaceStrikes = surface->getStrikes();
                for (i=0; i<surfaceStrikes.size(); i++) {
                    volStrikes->push_back(surfaceStrikes[i]);
                }
            }

            // Add the strikes explicitly supplied 
            for (i=0; i<strikes.size(); i++) {
                volStrikes->push_back(strikes[i]);
            }

            // Exclude those requested
            if (!!excludeStrikes) {
                for (i=0; i<excludeStrikes->size(); i++) { 
                    // remove() returns the iterator at the end of the shortened array
                    // erase() truncates the array at this point
                    volStrikes->erase(
                        remove(volStrikes->begin(),volStrikes->end(),(*excludeStrikes)[i]), 
                        volStrikes->end());
                }
            }

            // Sort from lowest to highest
            Algorithm::shellSort(*volStrikes);
            for (i=0; i<volStrikes->size()/2; i++) {
                double tmp = (*volStrikes)[i];
                (*volStrikes)[i]=(*volStrikes)[volStrikes->size()-i-1];
                (*volStrikes)[volStrikes->size()-i-1]=tmp;
            }

            // Remove duplicates
            volStrikes->erase(
                unique(volStrikes->begin(), volStrikes->end()),
                volStrikes->end());

            //for (i=1; i<volStrikes->size(); i++) {

              //  if ((*volStrikes)[i] == (*volStrikes)[i-1]) {
               //     volStrikes->erase(volStrikes->begin() + i);

            // Translate relative to absolute
            if (strikesInPct) {
                if (includeSurfaceStrikes) {
                    throw ModelException(method, "Cannot specify percentage strikes if existing "
                                                "vol surface strikes are to be used");
                }
                double spot = asset->getSpot();
                for (int i=0; i<volStrikes->size(); i++) {
                    (*volStrikes)[i] *= spot;
                }
            }            
        
            return volStrikes;
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    ImpliedVolSurface():  CObject(TYPE), volType(IVolatilityBS::TYPE->getName()), 
                                         ccyTreatment(CAsset::CCY_TREATMENT_NONE),
                                         includeSurfaceStrikes(false) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedVolSurface, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        FIELD(benchmarks, "Benchmarks for surface");
        FIELD(strikes, "Additional strikes for surface");
        FIELD(asset, "Asset");
        FIELD(market, "Market data");
        FIELD(volType, "Type of the vol in the market data cache");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(strikesInPct, "Are the strikes a percentage of spot");
        FIELD(ccyTreatment, "Currency treatment");
        FIELD_MAKE_OPTIONAL(ccyTreatment);
        FIELD(discount, "Payoff currency if struck");
        FIELD_MAKE_OPTIONAL(discount);
        FIELD(includeSurfaceStrikes, "Include existing vol surface strikes");
        FIELD_MAKE_OPTIONAL(includeSurfaceStrikes);
        FIELD(excludeStrikes, "Strikes to exclude");
        FIELD_MAKE_OPTIONAL(excludeStrikes);
        EMPTY_SHELL_METHOD(defaultImpliedVolSurface);

        Addin::registerClassObjectMethod("IMPLIED_VOL_SURFACE",
                                         Addin::RISK,
                                         "Returns an implied vol surface",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)ImpliedVolSurfaceAddin);
    }

    static IObject* defaultImpliedVolSurface(){
        return new ImpliedVolSurface();
    }
};

CClassConstSP const ImpliedVolSurface::TYPE = CClass::registerClassLoadMethod(
    "ImpliedVolSurface", typeid(ImpliedVolSurface), load);

DRLIB_END_NAMESPACE
