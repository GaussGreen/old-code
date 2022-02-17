//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRM.hpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef MCPathConfigSRM_HPP
#define MCPathConfigSRM_HPP

#include "edginc/MCPathConfigParametric.hpp"
#include "edginc/MCPathGenerator.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketDataFetcherSRM.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/BetaCorrelation.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/BasisIndexCurve.hpp"
#include "edginc/CIDParameters.hpp"
#include "edginc/QMCStrata.hpp"
#include ext_hash_set
#include <set>

DRLIB_BEGIN_NAMESPACE

/** Methods for discretising an SDE. Simplest is Euler-Maruyama method */
enum EquityDiffusionStyle {
    eqDEFAULT = 0,
    eqEULER,
    eqSECOND_ORDER, 
    eqSTRICT_SECOND_ORDER,
    eqMAPPING_METHOD
};

/** way of collecting market factors that we will model */
class MCARLO_DLL MFId{
public:
    string        name;
    CClassConstSP clazz;

    // method
    MFId(): clazz(0) {} // to allow hash_set to work
    MFId(const string& name, CClassConstSP clazz): 
    name(name), clazz(clazz){
        // backwards compatibility for StochasticYieldCurves - want
        // to treat all yield curves the same (or rather ignore 
        // StochasticYieldCurves)
        if (IYieldCurve::TYPE->isAssignableFrom(clazz)){
            this->clazz = IYieldCurve::TYPE;
        }
    }
    size_t operator()(const MFId& id) const {
        return (hash_string(id.name.c_str()) ^ (size_t) clazz);
    }
    bool operator()(const MFId& id1, const MFId& id2) const{
        return (id1.clazz == id2.clazz && id1.name == id2.name);
    }
};

typedef hash_set<MFId, MFId, MFId> MFHashSet;

class MCARLO_DLL MCPathConfigSRM: public MCPathConfigParametric,
                       public virtual IModelFamily { // for model family name
    friend class SV;
    friend class MCPathConfigSRMGen;
    friend struct SrmCorrData;
    
public:

    static CClassConstSP const TYPE;

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] 
    model for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const;

    //// helper method: If isoCode empty returns 0 otherwise returns
    //// the index in isoCode. Throws an exception if not found
    int findIsoCodeIndex(const string& isoCodeToFind) const;
	int findEnergyIndex(const string & energyName) const;

    virtual void validatePop2Object();

    const string getVolType() const { return volType; }
    const StringArray getFxVolType() const { return fxVolType; }
    const string getCDSVolType() const { return cdsVolType; }

    /** gets cdsVolType from the CreditModelType */
    string cdsVolTypeFromCreditModelType(
        const string& _creditModelType) const;

    // For the IModelFamily 
    const string* getFamilyName() const {
        return &MCPathConfigSRM::SRM_MODEL_FAMILY_NAME;
    }

    /** whether we're doing 'isCarefulRandoms' ie consistent random numbers
    when doing theta */
    virtual bool carefulRandoms() const{ return isCarefulRandoms; }

    /** Returns the amount of storage space the path generator needs to
    save paths (ie if PATH_CACHE bit is set in the call to 
    futurePathGenerator()). It does not need to be exact. */
    virtual int storagePerPath(IMCProduct* /*product*/) const {
        // to be done
        return 0;
    }

    /** Creates a PathGenerator used for historic dates */
    virtual MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod );

    /** Creates a PathGenerator used for future dates. The cachingMode
    parameter indicates what should be cached between tweaks. For example,
    if a path generator is
    created using with "cachingMode & PATH_RANDOM_ACCESS = true", then 
    subsequently another path generator can be created (with
    "cachingMode & PATH_RANDOM_ACCESS = false") which will
    support the pathIdx parameter in the generatePath() method, ie the
    paths can be generated in any order. */
    virtual MCPathGeneratorSP futurePathGenerator(
        int                                cachingMode,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates ); // see below

    /** Having retrieved a market object from the cache (or if presupplied) it
    is necessary to ensure that it has all its market data. */
    void getComponentMarketData(const IModel*         model,
                                const MarketData*     market,
                                MarketObjectSP        mo,
                                MarketDataFetcherSRM* mdf) const;

    //// returns the calibration style for specific currency
    const string& getCalibrationStyle(const string& isoCode) const;

    //// returns the calibration maturity for specific currency
    const string& getCalibrationMaturity(const string& isoCode) const;
	//// returns the calibration maturity for specific currency
	const string& getCalibrationMaturityCMS(const string& isoCode) const;

    int getNumIRFactors(const string& isoCode) const;

    bool calibrateAgainstSwaptionVols(const string& isoCode) const;

    //// returns the ir model params for specific currency
    const string& getIRModelParams(const string& isoCode) const;

    //// returns the ir smile params for specific currency
    const string& getIRSmileParams(const string& isoCode) const;

    //// returns the cr smile params for specific credit name
    const string& getCRSmileParams(const string& nameToFind) const;

    const string& getFXVolBootStrapMode(const string& isoCode) const;

    double getFXCutOffLevel(const string& isoCode) const;

    const string& getFXVolType(const string& isoCode) const;

    const string& getEqVolBootStrapMode(int index) const;

    double getEqCutOffLevel(int index) const;

    /** Returns GROWTH or DISCOUNT etc for supplied iso code */
    const string& getZCToDiffuse(const string& isoCode) const;

    /** Returns correlation swap start for supplied iso code */
    const string& getCorrelationSwapStart(const string& isoCode) const;

    /** Returns correlation swap maturity for supplied iso code */
    const string& getCorrelationSwapMat(const string& isoCode) const;

    /** Returns correlation swap day count convention for supplied iso code */
    const string& getCorrelationSwapDCC(const string& isoCode) const;

    /** Returns correlation swap frequency for supplied iso code */
    const string& getCorrelationSwapFreq(const string& isoCode) const;

    /** Returns CID parameters if they were provided */
    const ICIDParameters* getCIDParameters() const;


    /** Returns a VolProcessedBSIRSP for the supplied stochastic yield curve
    appropriate for the specified calibration style. Note returns null
    for calibration style NIL */
    CVolProcessedSP getProcessedIRVol(
        const string& isoCode,
        IYieldCurveConstSP   stochYC) const;

    bool skipNegIRVols(const string& isoCode) const;

    /** Invoked after instrument has got its market data. Allows model to
    get any extra data required. */
    virtual void getMarket(
        const IModel*      model,
        const MarketData*  market,
        IInstrumentCollectionSP instrument);

    /** For some vol models this may not make sense */
    virtual bool vegaMatrixSupported() const { return false; }

    /** Given the discount YC, what YC do we diffuse */
    IYieldCurveConstSP getDiffusedYC(IYieldCurveConstSP discount);


    /** few more inline accessors */
    
    const EquityDiffusionStyle& getEquityDiffusionStyle() const { return equityDiffusionStyle; }
    double getFlatVolIr() const { return flatVolIr; }
    const string& getChoiceCutoff() const { return choiceCutoff; }
    double getCutoffValue() const { return cutoffValue; }
    const CorrelationCommonArray& getCorrObjArray() const { return corrObjArray; }
    bool isStrictCorr() const { return strictCorr; }
    const string& getEqProtAdjMethod() { return eqProtAdjMethod; }
    double getNumSDMaxEffRateIR() const { return numSDMaxEffRateIR; }
    double getNumSDMinEffRateIR() const { return numSDMinEffRateIR; }
    double getNumSDMaxEffRateCR() const { return numSDMaxEffRateCR; }
    double getNumSDMinEffRateCR() const { return numSDMinEffRateCR; }
    double getNumSDMaxEffRateSP() const { return numSDMaxEffRateSP; }
    double getNumSDMinEffRateSP() const { return numSDMinEffRateSP; }
    double getBetaCorrMax() const { return betaCorrMax; }
    int getTimePointsPerYear() const { return timePointsPerYear; }

    string getCreditModelType() const { return creditModelType; }
    string getRatesModelType() const { return ratesModelType; }

    bool isFullCreditMC() const { return creditModelType == CREDIT_MODEL_TYPE_CID || 
        creditModelType == CREDIT_MODEL_TYPE_CID_FULL; }

    bool isCreditCID() const { return creditModelType == CREDIT_MODEL_TYPE_CID || 
        creditModelType == CREDIT_MODEL_TYPE_CID_FULL || 
        creditModelType == CREDIT_MODEL_TYPE_CID_FAST;}

    string getAllSupportedCreditModels() const {
        return CREDIT_MODEL_TYPE_CIR + " " + 
               CREDIT_MODEL_TYPE_HJM + " " +
               CREDIT_MODEL_TYPE_LIBOR + " " +
               CREDIT_MODEL_TYPE_CID + " " +
               CREDIT_MODEL_TYPE_CID_FULL + " " + 
               CREDIT_MODEL_TYPE_CID_FAST;
    }

    const QMCStratificationRecordArray& getStratification() const { return stratification; }


	//// returns the number of factor of the energy model
	int getEnergyNumFactor(const string & energyName) const;

	//// returns the start date of the energy correlation instrument
	const string & getEnergyCorrInstrStart(const string & energyName) const;

	//// returns the maturity of the energy correlation instrument
	const string & getEnergyCorrInstrMaturity(const string & energyName) const;

    virtual void createUtil( 
        IQMCGenDiffusibleAsset* gen, 
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    virtual DateTimeArray getCriticalDates(
        IQMCGenDiffusibleAsset* gen,
        const DateTime& start,       // likely to be "today"
        const DateTime& finish); // the latest requested date

    virtual void setDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        IQMCGenDiffusibleAsset* gen,
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);

    static IObject* defaultConstructor();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static const string CREDIT_MODEL_TYPE_HJM;
    static const string CREDIT_MODEL_TYPE_CIR;
    static const string CREDIT_MODEL_TYPE_LIBOR;
    static const string CREDIT_MODEL_TYPE_CID;
    static const string CREDIT_MODEL_TYPE_CID_FULL;
    static const string CREDIT_MODEL_TYPE_CID_FAST;

    static const string RATES_MODEL_TYPE_HJM;
    static const string RATES_MODEL_TYPE_LIBOR;
    static const string RATES_MODEL_TYPE_DETERM;

    static const string COMMON_MARKET_FACTOR_FOR_CID; // a token used in correlation of regular market data with CID's CM

private: 
    class MDF;
    /** Essentially a pass through for futurePathGenerator except that the
    relevant caches are created/updated etc */

    MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                  prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates );

   MCPathConfigSRM();

    ///// static fields  ////
    static const string CALIB_NIL; // "NIL" - no calibration against swaptns vol
    static const string ZC_GROWTH; // "GROWTH" - use 'growth' zero curve for diffusion
    static const string ZC_DISCOUNT; // "DISCOUNT" - use 'discount' instead
    static const string SRM_MODEL_FAMILY_NAME;

    /////// fields //////
    string          volType;    // optional
    StringArray     fxVolType;  // optional (possibly to be revisited)
    string          cdsVolType;    // optional
    int             timePointsPerYear;
    int             numICERuns;     /* how many times to run ICE. Default 2 */
    bool            isCarefulRandoms;
    EquityDiffusionStyle equityDiffusionStyle;  //Euler / 2nd order/ mapping method etc

    StringArray     calibrationStyle;  // eg CMS
    bool            useIRVolPair;	// HACK: to use IRVolPair object as well as IRVol (for regression tests).
    StringArray     calibrationMaturity; // eg 10Y
    StringArray     calibrationMaturityCMS; // eg 10Y
	
	StringArray     isoCode; /* either blank to use same for all currencies 
                             or matches up with calibrationStyle/
                             calibrationMaturity and other Arrays */
    IntArray        numIRFactors; // 1, 2, or 3 per iso code (or same for all)
    StringArray     irModelParams; // choice of inputs for IR model
    StringArray     irSmileParams; // choice of inputs for IR smile
    StringArray     correlationSwapStart; // eg 10Y (from valueDate)
    StringArray     correlationSwapMat;  // eg 10Y
    StringArray     correlationSwapDCC; // eg 30/360
    StringArray     correlationSwapFreq; // eg A or S etc
    BoolArray       skipNegativeIRVols;  // optional, default false
    double          flatVolIr; // optional, default 1.0

	// target swap rate correlations
	double          corrLower,
		            corrUpper;

    string          choiceCutoff;
    double          cutoffValue; 
    double          numSDMaxEffRateIR; // optional, default 2.0
    double          numSDMinEffRateIR; // optional, default 2.0
    StringArray     zcToDiffuse; // zc to diffuse (1 or n)
    StringArray     fxVolBootStrapMode; /* how to bootstrap fx vol (per foreign ccy) */
    DoubleArray     fxCutOffLevel;      /* to do (per foreign currency) */

    string          eqProtAdjMethod;    /* "DEFAULT" means local fx vol from simulation */
    StringArray     eqVolBootStrapMode; /* how to bootstrap eq vol (per equity) */
    DoubleArray     eqCutOffLevel;      /* per equity */

    StringArray     crSmileParams;      /* choice of inputs for CR smile (per credit) (optional)*/
    double          numSDMaxEffRateCR; // optional, default 2.0
    double          numSDMinEffRateCR; // optional, default 2.0
    double          numSDMaxEffRateSP; // optional, default 2.0
    double          numSDMinEffRateSP; // optional, default 0.0

    bool            deICE;  /* legacy Flag. optional, default false. Turn off ICE if true */
    bool            matchLegacySRM3; // true: match FI SRM3 numbers
    bool            offCycleICE;     // legacy flag. Do ICE offcycle
    string          dependenceType;
    string          creditModelType; // optional, default CIR
    string          ratesModelType; // optional, default HJM

	StringArray		energyNames;	// name should match energy curve names
	IntArray		energyNbFactors; // 1 or 2 for energy model
	StringArray		energyCorrInstrStarts; // energy corr instrument start date
	StringArray		energyCorrInstrMaturities; // enrg corr instrument maturity

    // an optional field that defines the stratified sampling
    QMCStratificationRecordArray  stratification;

    StringArray     rateDriverTargets;  // optional, target rates for driver rates
    IntArray        rateDrivers;        // optional, driver rate for rates

    StringArray     creditDriverTargets;
    IntArray        creditDrivers;      

    StringArray     energyDriverTargets;
    IntArray        energyDrivers;      

    // a field with CID parameters -- this is likely to be restructured
    CIDParametersWrapper    wCIDParametersWrapper;


    // For Sampras, many assets are diffused but there are very few
    // correlations.  Storing 0 correlations results in few Gigs of
    // memory being needed.  The following parameter was added to
    // address this.
    // When strictCorr = false, missing correlations are ignored
    // so that the associated corrObjects[pos] will be null.
    // When strictCorr = true, missing correlations result in exceptions.
    // True is the default value as this is the previous behavior in Qlib.
    bool            strictCorr;

    // transient fields
    /*   calibFlag in IR_INPUT - seems to be inferred from vol data */
    mutable MFHashSet           marketFactors; // for getting correlations
    mutable set<string>         yieldCurveNames; // for each ccy we use
    CorrelationCommonArray      corrObjArray;
    CorrelationCommonArray      eqeqCorrObjArray;
    BetaCorrelationArray        betaCorrArray;
    double                      betaCorrMax; // to prevent too high beta correlations, defaulted to 0.8
    CorrelationTermArray        eqeqCorrTermObjArray;
    YieldCurveArray             diffYCs;    // YC's to diffuse if not discount
    CAssetArray                 diffAssets; // assets to diffuse (equity & FX)
    ICDSParSpreadsArray         diffCDSs;   // CDS's to diffuse
	EnergyFuturesCurveArray		diffEnrgs;	// Energy future curves to diffuse
    IBasisIndexCurveArray       diffBasisArray; // Basis ccy curves to diffuse

};

typedef smartPtr<MCPathConfigSRM> MCPathConfigSRMSP;
typedef smartConstPtr<MCPathConfigSRM> MCPathConfigSRMConstSP;


DRLIB_END_NAMESPACE
#endif

