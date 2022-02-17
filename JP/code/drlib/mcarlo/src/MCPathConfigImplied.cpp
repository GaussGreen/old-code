//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathGenerator.cpp
//
//   Description : Implied Monte Carlo Path Generator
//
//
//   Date        : March 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/ImpliedSample.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductClient.hpp"
#include "edginc/MCProductEngineClient.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/DependenceLocalCorr.hpp"

#include "edginc/SVGenDiscFactor.hpp"
#define STATE_VARIABLES
#define STATEVAR_CACHING

#ifdef STATE_VARIABLES
#include "edginc/SVGenBarrier.hpp"
#include "edginc/HitSample.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#endif

#define PATH_DEBUG 0
#if PATH_DEBUG
#endif


DRLIB_BEGIN_NAMESPACE

// Available methodologies for caching strikes
static const string ALWAYS_CACHE_STRIKES = string("Always");            // Always cache strikes at tweaks
static const string CACHE_IF_POSSIBLE_STRIKES = string("IfPossible");   // Regenerate at extreme tweaks
static const string NEVER_CACHE_STRIKES = string("Never");              // Always throw away strikes
static const string DEFAULT_CACHE_STRIKES = string("Default");

/** Throws error if methodology is not one of the above */
static void checkCacheStrikesMethod(const string& cacheStrikes) {
    static const string method = "checkCacheStrikesMethod";

    bool invalid = !CString::equalsIgnoreCase(cacheStrikes, ALWAYS_CACHE_STRIKES) &&
                   !CString::equalsIgnoreCase(cacheStrikes, CACHE_IF_POSSIBLE_STRIKES) &&
                   !CString::equalsIgnoreCase(cacheStrikes, NEVER_CACHE_STRIKES);

    if(invalid) {
        throw ModelException("Unrecognied cacheStrikes method: " +
                             cacheStrikes +
                             ". Needs to be one of: " +
                             DEFAULT_CACHE_STRIKES + ", " +
                             ALWAYS_CACHE_STRIKES + ", " +
                             CACHE_IF_POSSIBLE_STRIKES + ", " +
                             NEVER_CACHE_STRIKES);
    }
}


///////////////////////////////////////////////
// Additional parameter values for MCImplied
///////////////////////////////////////////////
class ImpliedParams: public CObject {
public:
    static CClassConstSP const TYPE;

    /** For reflection */
    ImpliedParams(): CObject(TYPE), spotVol(default_spotVol), numNormals(default_numNormals),
                     callSpreadWidth(PDFRequestLNStrike::default_callSpreadWidth),
                     accuracy(PDFRequestLNStrike::default_accuracy) {}

    /** Full constructor */
    ImpliedParams(bool spotVol,
                  int numNormals,
                  double callSpreadWidth,
                  double accuracy):
        CObject(TYPE), spotVol(spotVol), numNormals(numNormals), callSpreadWidth(callSpreadWidth), accuracy(accuracy) {
        static string routine = "ImpliedParams::ImpliedParams";
        try {
            validatePop2Object();
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Basic checking of parameters */
    void validatePop2Object() {
        static string routine = "ImpliedParams::validatePop2Object";
        try {
            // Positive CallSpread width
            if(!Maths::isPositive(callSpreadWidth)) {
                throw ModelException("CallSpreadWidth " +
                                     Format::toString(callSpreadWidth) +
                                     " must be positive");
            }

            // Positive accuracy for callspreads
            if(!Maths::isPositive(accuracy)) {
                throw ModelException("Accuracy " +
                                     Format::toString(accuracy) +
                                     " must be positive");
            }

            // Need at least two points to discritize the Cumulative Normal function
            if(numNormals < 2 ) {
                throw ModelException("NumNormals " +
                                     Format::toString(numNormals) +
                                     " must be at least 2");
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    bool           spotVol;
    int            numNormals;
    double         callSpreadWidth;
    double         accuracy;

    static bool    default_spotVol;
    static int     default_numNormals;

private:
    /** Empty shell method */
    static IObject* defaultImpliedParams() {
        return new ImpliedParams();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(ImpliedParams, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultImpliedParams);
        FIELD(spotVol,  "Spot Vol");
        FIELD_MAKE_OPTIONAL(spotVol);
        FIELD(numNormals,  "Size of grid for cumulative normal");
        FIELD_MAKE_OPTIONAL(numNormals);
        FIELD(callSpreadWidth,  "Callspread width");
        FIELD_MAKE_OPTIONAL(callSpreadWidth);
        FIELD(accuracy,  "Accuracy");
        FIELD_MAKE_OPTIONAL(accuracy);
    }
};

bool   ImpliedParams::default_spotVol         = false;
int    ImpliedParams::default_numNormals      = 1000;

typedef smartPtr<ImpliedParams> ImpliedParamsSP;
typedef smartConstPtr<ImpliedParams> ImpliedParamsConstSP;

CClassConstSP const ImpliedParams::TYPE = CClass::registerClassLoadMethod(
    "ImpliedParams", typeid(ImpliedParams), ImpliedParams::load);


//////////////////////////////////
// DRIVERS VOLS
//////////////////////////////////
class DriversSqrtVar: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    DriversSqrtVar(const CAssetConstSP&        asset,
                   const VolRequestLNStrikeSP& volRequest,
                   const DateTime&             dateFrom,
                   const DateTimeArray&        datesTo,
                   const DoubleArray&          fwds,
                   bool                        spotVol,
                   bool                        isVolFwdStarting):
        CObject(TYPE) {
        static string routine = "DriversSqrtVar::DriversSqrtVar";

        try {
            int numDates = datesTo.size();

            if(numDates == 0) {
                throw ModelException("DatesTo is empty.");
            }

            // Count if there are any dates prior to the volStartDate
            DateTimeArray futDatesTo = dateFrom.getFutureDates(datesTo);
            int nbNullDates = numDates - futDatesTo.size();

            // Temporary arrays of size future dates
            int nbNonNullDates = numDates - nbNullDates;
            DoubleArray futSqrtTotalVar(nbNonNullDates);
            DoubleArray futSqrtFwdVar(nbNonNullDates);
            IntArray    futZeroTTime(nbNonNullDates);

            if(nbNonNullDates) {
                if(spotVol) {
                    generateVols(asset,
                                 volRequest,
                                 dateFrom,
                                 futDatesTo,
                                 isVolFwdStarting,
                                 futSqrtTotalVar,
                                 futSqrtFwdVar,
                                 futZeroTTime);

                } else {
                    // Use only the forwards for the future dates
                    DoubleArray futFwds(futDatesTo.size());
                    for(int iStep = 0; iStep < futDatesTo.size(); iStep++) {
                        futFwds[iStep] = fwds[iStep + nbNullDates];
                    }

                    generateVolsAtFwd(asset,
                                      volRequest,
                                      dateFrom,
                                      futDatesTo,
                                      isVolFwdStarting,
                                      futFwds,
                                      futSqrtTotalVar,
                                      futSqrtFwdVar,
                                      futZeroTTime);
                }
            }

            // Now copy over info and add some 0.0 for null dates
            sqrtTotalVar = DoubleArraySP(new DoubleArray(numDates));
            sqrtFwdVar   = DoubleArraySP(new DoubleArray(numDates));
            zeroTTime    = IntArraySP(new IntArray(numDates));

            int iStep;
            for(iStep = 0; iStep < nbNullDates; iStep++) {
                (*sqrtTotalVar)[iStep] = 0.0;
                (*sqrtFwdVar)[iStep]   = 0.0;
                (*zeroTTime)[iStep]    = 1;
            }
            for(; iStep < numDates; iStep++) {
                int index = iStep - nbNullDates;
                (*sqrtTotalVar)[iStep] = futSqrtTotalVar[index];
                (*sqrtFwdVar)[iStep]   = futSqrtFwdVar[index];
                (*zeroTTime)[iStep]    = futZeroTTime[index];
            }

            // Some validation for internal errors
            for(iStep = 0; iStep < numDates; iStep++) {
                if((*zeroTTime)[iStep] && !Maths::isZero((*sqrtTotalVar)[iStep])) {
                    throw ModelException("Zero trading time with non-zero variance encountered "
                                         "for maturity " + datesTo[iStep].toString() +
                                         " for asset " + asset->getName());
                }

                if(!(*zeroTTime)[iStep] && Maths::isZero((*sqrtTotalVar)[iStep])) {
                    throw ModelException("Non-zero trading time with zero variance encountered "
                                         "for maturity " + datesTo[iStep].toString() +
                                         " for asset " + asset->getName());
                }
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    DoubleArraySP              sqrtTotalVar;    //! Driver's total variance
    DoubleArraySP              sqrtFwdVar;      //! Driver's fwd variance
    IntArraySP                 zeroTTime;       //!< Whether trading time is zero

private:
    /** ATM Vol Interpolation for the driver */
    void generateVols(const CAssetConstSP&        asset,
                      const VolRequestLNStrikeSP& volRequest,
                      const DateTime&             dateFrom,
                      const DateTimeArray&        dates,
                      bool                        isVolFwdStarting,
                      DoubleArray&                sqrtTotalVar,
                      DoubleArray&                sqrtFwdVar,
                      IntArray&                   zeroTTime) {

        static const string routine("DriversSqrtVar::generateVols");

        try{
            // override with an ATM vol
            VolRequestLNStrikeSP requestATM(volRequest.clone());
            double strike = isVolFwdStarting? 1.0 : asset->getSpot();
            // BEWARE: the volRequests allow for negative FwdVariance
            // Here we really want positive fwd variance as we are constructing
            // a normal driver
            requestATM->allowNegativeFwdVar(false);
            requestATM->setStrike(strike);
            CVolProcessedSP vol(asset->getProcessedVol(requestATM.get()));
            CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));

            // get Total and Fwd variance
            DoubleArray totalVar(dates.size());
            DoubleArray fwdVar(dates.size());

            volBS->CalcVar(dateFrom, dates, volBS->fromFirst, totalVar);
            volBS->CalcVar(dateFrom, dates, volBS->forward,   fwdVar);

            for(int iStep = 0; iStep < dates.size(); iStep++) {
                sqrtTotalVar[iStep] = sqrt(totalVar[iStep]);
                sqrtFwdVar[iStep]   = sqrt(fwdVar[iStep]);

                if(Maths::isZero(sqrtTotalVar[iStep])) {
                    zeroTTime[iStep] = 1;
                } else {
                    zeroTTime[iStep] = 0;
                }
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** ATMFwd Vol Interpolation for the driver */
    void generateVolsAtFwd(const CAssetConstSP&        asset,
                           const VolRequestLNStrikeSP& volRequest,
                           const DateTime&             dateFrom,
                           const DateTimeArray&        dates,
                           bool                        isVolFwdStarting,
                           const DoubleArray&          fwds,
                           DoubleArray&                sqrtTotalVar,
                           DoubleArray&                sqrtFwdVar,
                           IntArray&                   zeroTTime) {

        static const string routine("DriversSqrtVar::generateVolsAtFwd");
        // The first date in the interpolation schedule is either
        // i) the first day in the list for started and FwdStarting products or
        // ii) the first AvgIn date for AvgIn products

        // Note: we are constructing the vol interp step-by-step as the interpolation
        // level changes for each maturity (remember we interpolate at Fwd)
        try{
            // Create temporary arrays for computing variance (cleaner code)
            DoubleArray totalVar(sqrtTotalVar.size());
            DoubleArray fwdVar(sqrtFwdVar.size());

            DoubleArray vars(1);
            DateTimeArray datesTo(1);

            // Set up vol request
            VolRequestLNStrikeSP requestATM(volRequest.clone());
            requestATM->allowNegativeFwdVar(false);

            // Compute fwd at vol start
            double fwdAtVolStart = asset->fwdValue(dateFrom);

            // 1) Compute the total variance at Fwd from startDate to each maturity
            int iStep;

            for(iStep = 0; iStep < dates.size(); iStep++) {
                // decide on relative vs absolute level for interpolation
                double strike = isVolFwdStarting ?
                    fwds[iStep] / fwdAtVolStart:
                    fwds[iStep];

                requestATM->setStrike(strike);
                CVolProcessedSP vol(asset->getProcessedVol(requestATM.get()));
                CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));

                datesTo[0] = dates[iStep];
                volBS->CalcVar(dateFrom, datesTo, volBS->fromFirst, vars);

                totalVar[iStep] = vars[0];
            }

            // 2) Compute fwd variance between successive maturities by subtracting
            //    totalVariances and flooring it to 0.0
            fwdVar[0] = totalVar[0];
            for(iStep = 1; iStep < dates.size(); iStep++) {
                // definition of FwdVariance as increment in Total Variance
                double thisFwdVar = totalVar[iStep] - totalVar[iStep-1];
                if(!Maths::isPositive(thisFwdVar)) {
                    // Floor fwd variance to zero
                    thisFwdVar = Maths::max(thisFwdVar, 0.0);
                }
                fwdVar[iStep] = thisFwdVar;
            }

            // 3) Recompute total variance by adding up fwdVariances as these may
            //    have been affected by flooring. Also deduce zero trading time now
            //    that total variance has been finalized
            for(iStep = 0; iStep < dates.size(); iStep++) {
                double thisTotalVar = fwdVar[iStep];
                if(iStep > 0) {
                    // Plus previous total variance
                    thisTotalVar += totalVar[iStep-1];
                }

                // Now apply sqrt and update the output arrays
                sqrtTotalVar[iStep] = sqrt(totalVar[iStep]);
                sqrtFwdVar[iStep] = sqrt(fwdVar[iStep]);

                // Zero trading time if total variance is zero
                if(Maths::isZero(thisTotalVar)) {
                    zeroTTime[iStep] = 1;
                } else {
                    zeroTTime[iStep] = 0;
                }
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** For Reflection */
    DriversSqrtVar(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(DriversSqrtVar, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDriversSqrtVar);
        FIELD(sqrtTotalVar, "SqRoot of total variance");
        FIELD(sqrtFwdVar,   "SqRoot of fwd variance");
        FIELD(zeroTTime,    "Zero trading time");
    }

    /** Empty shell method */
    static IObject* defaultDriversSqrtVar(){
        return new DriversSqrtVar();
    }
};

typedef smartPtr<DriversSqrtVar> DriversSqrtVarSP;
typedef smartConstPtr<DriversSqrtVar> DriversSqrtVarConstSP;

CClassConstSP const DriversSqrtVar::TYPE = CClass::registerClassLoadMethod(
    "DriversSqrtVar", typeid(DriversSqrtVar), load);


//////////////////////////////////////////////////
// Cache per Asset:
// Contains ImpliedVolatility for driver and
// Array of implied samples
//////////////////////////////////////////////////
class AssetCache: public CObject {
public:
    static CClassConstSP const TYPE;

    /** For relflection */
    AssetCache(): CObject(TYPE), numPDFSteps(-1) { };

    /** Overrides clone method to copy the Lattices */
    IObject* clone() const {
        static const string routine = "AssetCache::clone";
        try {
            IObject* Copy = CObject::clone(); // call parent

            AssetCache* copy = dynamic_cast<AssetCache*>(Copy);
            if(!copy) {
                throw ModelException("Clone method failed");
            }
            if(strikes.get()) {
                // now fill in the strikes lattice
                copy->strikes = CLatticeDoubleSP(new CLatticeDouble(*strikes));
            }
            if(probs.get()) {
                // now fill in the probs lattice
                copy->probs   = CLatticeDoubleSP(new CLatticeDouble(*probs));
            }

            return copy;
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    int                numPDFSteps;     // number of steps for PDF

    DriversSqrtVarSP   sqrtVars;        // driver's vol

    CLatticeDoubleSP   strikes;         // lattice of strikes $unregistered
    StrikesPartitionSP partition;       // partition of strikes
    CLatticeDoubleSP   probs;           // lattice of probs $unregistered

    LinearImpliedSampleArraySP pdfSamples;

private:
    /** Empty shell method */
    static IObject* defaultAssetCache(){
        return new AssetCache();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(AssetCache, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAssetCache);
        FIELD(numPDFSteps, "Number of PDFSteps");
        FIELD_MAKE_TRANSIENT(numPDFSteps);
        FIELD(sqrtVars, "Driver's SqRoot Var");
        FIELD_MAKE_TRANSIENT(sqrtVars);
        FIELD(partition, "Strikes partition");
        FIELD_MAKE_TRANSIENT(partition);
        FIELD(pdfSamples, "Array of Linear Implied Samples");
        FIELD_MAKE_TRANSIENT(pdfSamples);
    }
};

typedef smartPtr<AssetCache> AssetCacheSP;
typedef smartConstPtr<AssetCache> AssetCacheConstSP;

CClassConstSP const AssetCache::TYPE = CClass::registerClassLoadMethod(
    "AssetCache", typeid(AssetCache), load);

typedef array<AssetCacheSP, AssetCache> AssetCacheArray;
typedef smartPtr<AssetCacheArray> AssetCacheArraySP;
typedef smartConstPtr<AssetCacheArray> AssetCacheArrayConstSP;

DEFINE_TEMPLATE_TYPE(AssetCacheArray);


//////////////////////////////////////////////////
// Cache per Cliquet:
// Contains array of caches per asset
//////////////////////////////////////////////////
class CliquetCache: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    CliquetCache(const DateTime& startDate,
                 int             numSteps,
                 int             numAssets): // = numFutSteps at first Pricing
        CObject(TYPE), startDate(startDate), numSteps(numSteps), numAssets(numAssets) {

        static const string routine = "CliquetCache::CliquetCache";

        try {
            assetCaches = AssetCacheArraySP(new AssetCacheArray(numAssets));
            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                (*assetCaches)[iAsset] = AssetCacheSP(new AssetCache());
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }


    // All data are public so that they can accessed by MCPathBaseImplied
    // and MCPathConfigImplied
    AssetCacheArraySP assetCaches; // asset specific information
    DateTime          startDate;   // startDate at pricing time
    int               numSteps;    // the number of steps in the PRICING timeline
    int               numAssets;   // the number of assets

private:
    /** For reflection */
    CliquetCache(): CObject(TYPE) {}

    static IObject* defaultCliquetCache(){
        return new CliquetCache();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(CliquetCache, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCliquetCache);
        FIELD(assetCaches, "Array of asset caches");
        FIELD_MAKE_TRANSIENT(assetCaches);
        FIELD(startDate, "Start date");
        FIELD(numSteps, "Number of future steps");
        FIELD(numAssets, "Number of assets");
    }
};


typedef smartPtr<CliquetCache> CliquetCacheSP;
typedef smartConstPtr<CliquetCache> CliquetCacheConstSP;

CClassConstSP const CliquetCache::TYPE = CClass::registerClassLoadMethod(
    "CliquetCache", typeid(CliquetCache), load);

typedef array<CliquetCacheSP, CliquetCache> CliquetCacheArray;
typedef smartPtr<CliquetCacheArray> CliquetCacheArraySP;
typedef smartConstPtr<CliquetCacheArray> CliquetCacheArrayConstSP;

DEFINE_TEMPLATE_TYPE(CliquetCacheArray);


///////////////////////////////////////////////////
// Contains general information stored across
// tweaks AND array of CliquetCaches
///////////////////////////////////////////////////
class MCImpliedCache: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    MCImpliedCache(int numCliquets):
    CObject(TYPE), numCliquets(numCliquets), cliquetCaches(new CliquetCacheArray(numCliquets)),
    isResultsExported(true) {}

    /** Routine that kills the driver's vols for all assets */
    void killVols() {
        if(cliquetCaches.get()) {
            for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                CliquetCacheSP& cliquetCache = (*cliquetCaches)[iCliquet];
                if(cliquetCache.get() && cliquetCache->assetCaches.get()) {
                    for(int iAsset = 0; iAsset < cliquetCache->numAssets; iAsset++) {
                        (*cliquetCache->assetCaches)[iAsset]->sqrtVars = DriversSqrtVarSP(   );
                    }
                }
            }
        }
    }

    /** Routine that kills the implied probs for all assets */
    void killProbs() {
        if(cliquetCaches.get()) {
            for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                CliquetCacheSP& cliquetCache = (*cliquetCaches)[iCliquet];
                if(cliquetCache.get() && cliquetCache->assetCaches.get()) {
                    for(int iAsset = 0; iAsset < cliquetCache->numAssets;
                        iAsset++) {
                        (*cliquetCache->assetCaches)[iAsset]->probs =
                            CLatticeDoubleSP();
                    }
                }
            }
        }
    }

    /** Routine that kills the implied probs for assets listed in array */
    void killProbs(const IntArray& sensitiveAssets) {
        if(cliquetCaches.get()) {
            for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                CliquetCacheSP& cliquetCache = (*cliquetCaches)[iCliquet];
                if(cliquetCache.get() && cliquetCache->assetCaches.get()) {
                    for(int i = 0; i < sensitiveAssets.size(); i++) {
                        int iAsset = sensitiveAssets[i];
                        (*cliquetCache->assetCaches)[iAsset]->probs =
                            CLatticeDoubleSP();
                    }
                }
            }
        }
    }

    /** Routine that kills the cliquet caches */
    void killCliquetCaches(int numCliquetsInput) {
        numCliquets = numCliquetsInput;
        cliquetCaches->resize(numCliquetsInput);
        if(cliquetCaches.get()) {
            for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                (*cliquetCaches)[iCliquet] = CliquetCacheSP(   );
            }
        }
    }

    // All data are public so that they can accessed by MCPathBaseImplied
    // and MCPathConfigImplied
    int                             numCliquets;        // number of cliquets
    CliquetCacheArraySP             cliquetCaches;      // caches per cliquet
    bool                            isResultsExported;

private:
    /** For reflection */
    MCImpliedCache(): CObject(TYPE) {}

    static IObject* defaultMCImpliedCache(){
        return new MCImpliedCache();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCImpliedCache, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMCImpliedCache);
        FIELD(numCliquets, "Number of cliquets");
        FIELD(cliquetCaches, "Array of cliquet caches");
        FIELD_MAKE_TRANSIENT(cliquetCaches);
        FIELD(isResultsExported, "Has info been exported to DEBUG packet");
    }
};


typedef smartPtr<MCImpliedCache> MCImpliedCacheSP;
typedef smartConstPtr<MCImpliedCache> MCImpliedCacheConstSP;

CClassConstSP const MCImpliedCache::TYPE = CClass::registerClassLoadMethod(
    "MCImpliedCache", typeid(MCImpliedCache), load);


/////////////////////////////////////////
// DebugInfo:
// Provides a blueprint for constructing
// a LinearImpliedSample
/////////////////////////////////////////
class DebugInfo: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    DebugInfo(const StrikesPartitionSP& partition,
              const DoubleArraySP&      adjustment) :
        CObject(TYPE), partition(partition), adjustment(adjustment) {}

    StrikesPartitionSP      partition;
    DoubleArraySP           adjustment;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(DebugInfo, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDebugInfo);
        FIELD(partition,  "Strikes Partition");
        FIELD(adjustment, "Fwd Adjustments");
    }

    DebugInfo(): CObject(TYPE) {}

    static IObject* defaultDebugInfo() {
        return new DebugInfo();
    }
};

typedef smartPtr<DebugInfo> DebugInfoSP;
typedef smartConstPtr<DebugInfo> DebugInfoConstSP;

CClassConstSP const DebugInfo::TYPE = CClass::registerClassLoadMethod(
    "DebugInfo", typeid(DebugInfo), load);

typedef array<DebugInfoSP, DebugInfo> DebugInfoArray;
typedef smartPtr<DebugInfoArray> DebugInfoArraySP;
typedef smartConstPtr<DebugInfoArray> DebugInfoArrayConstSP;

DEFINE_TEMPLATE_TYPE(DebugInfoArray);


//////////////////////////////////////////
// Max Drifts required for QuickGreeks
//////////////////////////////////////////
class ImpliedMaxDrifts {
public:
    /** Constructor */
    ImpliedMaxDrifts(
        const vector<ImpliedSampleArray>& samples,          // implied samples [iAsset][iStep]
        const CDoubleArray&               fwdAtSimStart):   // fwd at sim start date
    maxDrifts(samples.size()) {

        static string routine = "ImpliedMaxDrifts::ImpliedMaxDrifts";
        try {
            /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
                the 'drift' between simulation date i and simulation date i+1.
                (With d_0 = simulation start date). Range of dates is restricted to
                those in the future so the past path generator will
                return 1.0. The 'drift' is the biggest ratio of the 2 spots possible
                within some sort of probability interval. See  MCPathConfigLN */
            double criticalPercentile = 1.0 - N1(normalPercentile);
            for(unsigned int iAsset = 0; iAsset < samples.size(); iAsset++) {
                const ImpliedSampleSP& lastSample = samples[iAsset].back();
                double maxDrift = lastSample->sample(criticalPercentile);
                if(lastSample->isRelativeReturn() == false) {
                    maxDrift /= fwdAtSimStart[iAsset];
                }
                maxDrifts[iAsset] = maxDrift;
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Returns the maxDrift for each asset */
    double getMaxDrift(int iAsset) {
        return maxDrifts[iAsset];
    }

private:
    CDoubleArray  maxDrifts;
    static double normalPercentile;
};

DECLARE_REF_COUNT(ImpliedMaxDrifts);

double ImpliedMaxDrifts::normalPercentile = 3.0; // same as in MCPathConfigLN


///////////////////////////////////////////////////////
// Contains all the data required for each cliquet.
///////////////////////////////////////////////////////
class Cliquet {
    // some code that was moved from Cliquet() as solaris opt was crashing
    // when an exception got thrown from getPartition
    void helper(AssetCacheSP         assetCache,
                VolRequestLNStrikeSP volRequest,
                ImpliedParamsSP      impliedParams,
                const DateTime&      today,
                const DateTimeArray& PDFTimeline,
                PDFParamsSP          pdfParams,
                DoubleArraySP        pdfFwds,
                DoubleArraySP        pdfSqrtVar,
                CAssetConstSP        asset,
                const string&        cacheStrikes){
        try {
            PDFRequestLNStrikeSP pdfRequest(
                new PDFRequestLNStrike(volRequest.get(),
                                       startDate,
                                       impliedParams->callSpreadWidth,
                                       impliedParams->accuracy));

            // Check if we need to regenerate the strikes
            bool useExistingStrikes = true;
            if(CString::equalsIgnoreCase(cacheStrikes, NEVER_CACHE_STRIKES) &&
               assetCache->strikes.get()) {
                // If we never want to cache them, request regeneration anyway
                useExistingStrikes = false;
            } else if(CString::equalsIgnoreCase(cacheStrikes, CACHE_IF_POSSIBLE_STRIKES) &&
                      assetCache->strikes.get()) {
                // See whether we will use the existing strikes or regenerate them
                // Regeneration criterion is based on fwd

                for(int iStep = 0; iStep < assetCache->strikes->size(); iStep++) {
                    CSliceDouble& strk = (*assetCache->strikes)[iStep];
                    int nbStrikes = strk.size();

                    // If fwd is near lowStrike or highStrike within 10%
                    // we need to regenerate the strikes as the resulting
                    // distribution will be very weird using the original
                    // strikes. This happens when the Delta shift sizes are
                    // so big that the spot moves near the boundaries of the
                    // support of the original distribution
                    double pct = ( (*pdfFwds)[iStep]  - strk[0]) /
                                   (strk[nbStrikes-1] - strk[0]);

                    if( pct < 0.1 || pct > 0.9 ) {
                        // request regeneration of whole lattice of strikes
                        useExistingStrikes = false;
                        break;
                    }
                }
            }

            CLatticeDoubleSP   inStrikes;
            CLatticeDoubleSP   inProbs;
            StrikesPartitionSP inPartition;
            if(useExistingStrikes) {
                // Get them from cache
                inStrikes   = assetCache->strikes;
                inProbs     = assetCache->probs;
                inPartition = assetCache->partition;
            }

            LinearImpliedSampler assetSampler(asset,
                                              volRequest,
                                              today,
                                              startDate,
                                              PDFTimeline,
                                              *pdfParams,
                                              pdfRequest,
                                              pdfFwds,
                                              pdfSqrtVar,
                                              inStrikes,
                                              inProbs,
                                              inPartition);

            // extract information from the sampler
            CLatticeDoubleSP   outStrikes    = assetSampler.getStrikes();
            CLatticeDoubleSP   outProbs      = assetSampler.getProbs();
            StrikesPartitionSP outPartition;
            try {
                outPartition  = assetSampler.getPartition();
            } catch (exception&){
                throw;
            }
            LinearImpliedSampleArraySP outPdfSamples =
                assetSampler.getImpliedSamples();

            if(useExistingStrikes) {
                // record strikes in the cache if we did not regenerate them
                // otherwise leave original in there
                assetCache->strikes    = outStrikes;
                assetCache->probs      = outProbs;
                assetCache->partition  = outPartition;
            }
            // Always write the samples as these are used inside the simulation
            assetCache->pdfSamples = outPdfSamples;
        } catch(exception& e) {
            string message = "Failed to create Linear Implied Samples "
                " for the second pass for asset " +
                asset->getName();
            throw ModelException::addTextToException(e, message);
        }
    }

    void populate(int                                 iCliquet,
                  const vector<VolRequestLNStrikeSP>& volRequests,
                  const IMultiFactors*                mAsset,
                  const DateTime&                     today,
                  PDFParamsSP                         pdfParams,
                  ImpliedParamsSP                     impliedParams,
                  CliquetCacheSP&                     cliquetCache,
                  const string&                       cacheStrikes,
                  EquidistantLinearInterpolantNonVirtualSP tabulatedNorm) {

        static const string routine = "Cliquet::populate";

        try {
            checkCacheStrikesMethod(cacheStrikes);

            // get volStartDate and check that it is identical for all assets
            const DateTime& volStartDate = volRequests[0]->getStartDate();
            int iAsset;
            for(iAsset = 1; iAsset < numAssets; iAsset++) {
                if(volStartDate != volRequests[iAsset]->getStartDate()) {
                    throw ModelException("All vol requests must have same start date. Check asset " +
                                         Format::toString(iAsset + 1));
                }
            }

            // check that volStartDate is in the simulation dates list
            if(volStartDate > today) {
                int occurances = count(datesTo.begin(),
                                       datesTo.end(),
                                       volStartDate);

                if(occurances == 0 && volStartDate != startDate) {
                    throw ModelException(routine,
                                         "Vol start date " +
                                         volStartDate.toString() +
                                         " is not a member of the Simulation "
                                         "Dates for cliquet " +
                                         Format::toString(iCliquet) +
                                         ".");
                }
            }

            // Definition of FwdStartingVol
            bool isVolFwdStarting = volStartDate > today;
            // Definition of FwdStarting Product (According to the definition, AvgIn is not FwdStarting)
            bool isFwdStarting = startDate > today ? true : false;

            // initialize the cache, unless it exists from some previous tweak
            if(!cliquetCache) {
                cliquetCache = CliquetCacheSP(new CliquetCache(startDate,
                                                               numFutSteps,
                                                               numAssets));
            }

            // the date from which we compute the driver
            DateTime dateFrom = today >= volStartDate ? today : volStartDate;

            // Check if the timeline has changed. This can happen
            // when rolling over simulation dates at THETA type tweaks
            // Note we have already captured rolling over a cliquet
            // Here we check for rolling over a date within a cliquet
            bool changedTimeline = cliquetCache->numSteps != numFutSteps ||
                (startDate != cliquetCache->startDate &&
                 iCliquet > 0);
            int iStep;
            // Asset specific from now on
            DoubleArray fwdAtSimStart(numAssets);
            DoubleArray fwdAtVolStart(numAssets);
            for(iAsset = 0; iAsset < numAssets; iAsset++) {
                CAssetConstSP asset      = CAssetConstSP(&mAsset->getAsset(iAsset));
                AssetCacheSP  assetCache = (*cliquetCache->assetCaches)[iAsset];
                string        assetName  = asset->getName();

                // compute the forwards
                fwds[iAsset] = DoubleArray(datesTo.size());
                mAsset->factorFwdValues(iAsset, datesTo, fwds[iAsset]);
                fwdAtSimStart[iAsset] = mAsset->factorFwdValue(iAsset, startDate);
                fwdAtVolStart[iAsset] = mAsset->factorFwdValue(iAsset, dateFrom);

                if(isFwdStarting) {
                    relativeLevel[iAsset] = fwdAtSimStart[iAsset];
                } else if(isVolFwdStarting) {
                    relativeLevel[iAsset] = fwdAtVolStart[iAsset];
                } else {
                    relativeLevel[iAsset] = 1.0;    // it's absolute
                }

                // get the vol for the driver at common dates and find if zeroTTime
                if(!assetCache->sqrtVars) {
                    try {
                        assetCache->sqrtVars = DriversSqrtVarSP(new DriversSqrtVar(
                                                                    asset,
                                                                    volRequests[iAsset],
                                                                    dateFrom,
                                                                    datesTo,
                                                                    fwds[iAsset],
                                                                    impliedParams->spotVol,
                                                                    isVolFwdStarting));
                    } catch(exception& e) {
                        string message = "Failed to compute driver's vols for asset " + assetName;
                        throw ModelException::addTextToException(e, message);
                    }
                }

                // compute asset's timeline by dropping more dates from the commonTimeline
                // depending on zero Trading Time for example
                DateTimeArray PDFTimeline;
                IntArray      PDFIndex;
                IntArray      droppedIndex;
                DoubleArraySP pdfFwds(new DoubleArray(0));
                DoubleArraySP pdfSqrtVar(new DoubleArray(0));

                DoubleArray&  driverSqrtTotalVar = *assetCache->sqrtVars->sqrtTotalVar;
                for(iStep = 0; iStep < datesTo.size(); iStep++) {
                    // drop dates with zeroTTime
                    bool zeroTTime = (*assetCache->sqrtVars->zeroTTime)[iStep] ? true : false;
                    if(zeroTTime) {
                        droppedIndex.push_back(iStep);
                    } else {
                        PDFTimeline.push_back(datesTo[iStep]);
                        PDFIndex.push_back(iStep);

                        pdfFwds->push_back(fwds[iAsset][iStep]);
                        pdfSqrtVar->push_back(driverSqrtTotalVar[iStep]);
                    }
                }

                int numPDFSteps = PDFTimeline.size();
                if(assetCache->numPDFSteps < 0) {
                    assetCache->numPDFSteps = numPDFSteps;
                }

                // check if the pdf steps have changed
                // this could happen without rolling over a date but
                // when we approach very near a sample date so we have
                // zero trading time
                if(assetCache->strikes.get()) {
                    bool killStrikes = assetCache->numPDFSteps != numPDFSteps ||
                        changedTimeline;
                    if(killStrikes) {
                        assetCache->strikes   = CLatticeDoubleSP();
                        assetCache->partition = StrikesPartitionSP();
                        assetCache->probs     = CLatticeDoubleSP();
                    }
                }

                // If we have PDFDates, get strikes, partition, probs and samples
                if(numPDFSteps != 0) {
                    try{
                        helper(assetCache, volRequests[iAsset], impliedParams,
                               today, PDFTimeline, pdfParams, pdfFwds,
                               pdfSqrtVar, asset, cacheStrikes);
                    } catch (exception&){
                        throw;
                    }
                }

                // create ImpliedSample array for the whole timeLine.
                // It consists of:
                // i) LinearImpliedSamples
                // ii) DegenerateImplied samples
                ImpliedSampleArray&  assetSamples = samples[iAsset];
                ImpliedMappingArray& assetMappings = mappings[iAsset];
                assetSamples  = ImpliedSampleArray(numFutSteps);
                assetMappings = ImpliedMappingArray(numFutSteps);
                int iPdf = 0;
                int iDet = 0;
                for(iStep = 0; iStep < numFutSteps; iStep++) {
                    if(iPdf < numPDFSteps && PDFIndex[iPdf] == iStep) {
                        // samples from the Implied Distribution
                        LinearImpliedSampleSP linSample = (*assetCache->pdfSamples)[iPdf];
                        assetSamples[iStep]  = linSample;
                        assetMappings[iStep] = ImpliedMappingSP(new
                            GaussianInterpImpliedMapping(
                            tabulatedNorm,
                            linSample,
                            (*pdfSqrtVar)[iPdf]));
                        iPdf++;
                    } else if (droppedIndex.size() != 0) {
                        // deterministic Samples
                        int index = droppedIndex[iDet];

                        double fwd;
                        bool isRelative;
                        if(isVolFwdStarting) {
                            // relative returns except prior to volStartDate e.g. Today, 1st AvgIn day
                            if(datesTo[index] <= volStartDate) {
                                fwd = fwds[iAsset][index];
                                isRelative = false;
                            } else {
                                fwd = fwds[iAsset][index] / fwdAtVolStart[iAsset];
                                isRelative = true;
                            }
                        } else {
                            // absolute
                            fwd = fwds[iAsset][index];
                            isRelative = false;
                        }

                        DegenerateImpliedSampleSP degSample(new DegenerateImpliedSample(fwd, isRelative));
                        assetSamples[iStep] = DegenerateImpliedSampleSP(degSample);
                        assetMappings[iStep] = ImpliedMappingSP(new
                            DegenerateImpliedMapping(degSample));
                        iDet++;
                    } else {
                        throw ModelException("Internal error. Samples should either be "
                                             "linear or degenerate.");
                    }
                }
            }

            // Store whether it is relative return or absolute price
            for(iStep = 0; iStep < numFutSteps; iStep++) {
                isRelativeReturn[iStep] = samples[0][iStep]->isRelativeReturn();
            }

            // Record the information for quick greeks
            maxDrifts = ImpliedMaxDriftsSP(new ImpliedMaxDrifts(samples, fwdAtSimStart));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }
public:
    /** Constructor that initializes reference to fwds. To be used for
        cliquets that want a deterministic reference e.g. 1st cliquet
        with AvgIn uses Fwd at vol start date. 2nd, 3rd etc. want a
        stochastic reference i.e. spot at end of previous cliquet so
        should use the other constructor */
    Cliquet(int                                 iCliquet,        // which cliquet
            const vector<VolRequestLNStrikeSP>& volRequests,     // volReqests per asset
            const IMultiFactors*                mAsset,          // multi asset
            const DateTime&                     today,           // valuation day
            const DateTime&                     startDate,       // product start day
            const DateTimeArray&                datesTo,         // maturities for this cliquet
            const PDFParamsSP&                  pdfParams,       // pdf params
            const ImpliedParamsSP&              impliedParams,   // implied distribution params
            CliquetCacheSP&                     cliquetCache,    // cache for writing info
            const string&                       cacheStrikes,    // Whether to cache strikes between tweaks
            EquidistantLinearInterpolantNonVirtualSP tabulatedNorm):  // Tabulated normal distribution
        numAssets(mAsset->NbAssets()), numFutSteps(datesTo.size()), startDate(startDate),
        datesTo(datesTo), fwds(numAssets), samples(numAssets), mappings(numAssets),
        isRelativeReturn(numFutSteps), relativeLevel(numAssets) {

        static const string routine("Cliquet::Cliquet");
        try{
            // Populate the fields
            populate(iCliquet, volRequests, mAsset, today, pdfParams, impliedParams,
                cliquetCache, cacheStrikes, tabulatedNorm);

            // Set the return Reference to the relativeLevel i.e. deterministic
            // level equal to some type of fwd (simStart or volStart)
            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                returnReference.push_back(&relativeLevel[iAsset]);
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Constructor that initializes reference to point to some part
        of the path. To be used for cliquets that want a stochastic
        reference e.g. 2nd, 3rd etc. whose reference is spot at end
        of previous cliquet */
    Cliquet(int                                 iCliquet,        // which cliquet
            const vector<VolRequestLNStrikeSP>& volRequests,     // volReqests per asset
            const IMultiFactors*                mAsset,          // multi asset
            const DateTime&                     today,           // valuation day
            const DateTime&                     startDate,       // product start day
            const DateTimeArray&                datesTo,         // maturities for this cliquet
            const PDFParamsSP&                  pdfParams,       // pdf params
            const ImpliedParamsSP&              impliedParams,   // implied distribution params
            CliquetCacheSP&                     cliquetCache,    // cache for writing info
            const string&                       cacheStrikes,    // Whether to cache strikes between tweaks
            EquidistantLinearInterpolantNonVirtualSP tabulatedNorm, // Tabulated normal distribution
            const vector<const double*>&        stochasticRelativeLevel): // relative level
        numAssets(mAsset->NbAssets()), numFutSteps(datesTo.size()), startDate(startDate),
        datesTo(datesTo), fwds(numAssets), returnReference(stochasticRelativeLevel),
        samples(numAssets), mappings(numAssets), isRelativeReturn(numFutSteps),
        relativeLevel(numAssets) {

        static const string routine("Cliquet::Cliquet");
        try{
            // Populate the fields
            populate(iCliquet, volRequests, mAsset, today, pdfParams, impliedParams,
                cliquetCache, cacheStrikes, tabulatedNorm);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    int                         numAssets;           //!< Number of assets
    int                         numFutSteps;         //!< Number of maturities

    DateTime                    startDate;           //!< Product's start date
    DateTimeArray               datesTo;             //!< Maturities for this cliquet

    vector<DoubleArray>         fwds;                //!< Fwds at datesTo [iAsset][iStep]
    vector<const double *>      returnReference;     //!< Either input or relativeLevel (See below)

    vector<ImpliedSampleArray>  samples;             //!< Implied sample  [iAsset][iStep]
    vector<ImpliedMappingArray> mappings;            //!< Implied mapping [iAsset][iStep]
    IntArray                    isRelativeReturn;    //!< Whether date is relative return or absolute [iStep]
    ImpliedMaxDriftsSP          maxDrifts;           //!< Required for quickGreeks

private:
    DoubleArray                 relativeLevel;       //!< Relative level w.r.t. which returns are computed
};

DECLARE_REF_COUNT(Cliquet);

/////////////////////////////////////////
// MCPathConfigImplied
/////////////////////////////////////////
/** Class splits its work into three smaller classes. This is possibly more
    code but it makes it much simpler and easier to develop */
class MCPathConfigImplied: public MCPathConfig {
    friend class MCPathBaseImplied;
    MCCacheManager cacheMgr;        //!< Cache manager $unregistered
public:
    static CClassConstSP const TYPE;

    class Gen;              //!< Referee simulation object
    class PathGenSpot;      //!< Spot path simulator
    class PathGenHVBB;      //!< Brownian Bridge PathGen for hitting values

    friend class Gen;
    friend class PathGenSpot;
    friend class PathGenHVBB;

    typedef refCountPtr<Gen> GenSP;
    typedef refCountPtr<PathGenSpot> PathGenSpotSP;
    typedef refCountPtr<PathGenHVBB> PathGenHVBBSP;

protected:
    MCPathConfigImplied(CClassConstSP clazz,
                        const string& volType,
                        const string& dependenceType):
        MCPathConfig(clazz, DependenceMakerGaussTermSP(new DependenceMakerGaussTerm())),
        volType(volType),
        dependenceType("not used"),
        numMidStrikesPerLogStrike(PDFParams::default_numMidStrikesPerLogStrike),
        minMidStrikes(PDFParams::default_numMiddleStrikes),
        numTailStrikesPerLogStrike(PDFParams::default_numTailStrikesPerLogStrike),
        minTailStrikes(PDFParams::default_numTailStrikes),
        numMidStdDevs(PDFParams::default_numMiddleStdDevs),
        numTailStdDevs(PDFParams::default_numTailStdDevs),
        cacheStrikes(DEFAULT_CACHE_STRIKES),
        DEBUG_propMidStrikes(PDFParams::default_propMiddleStrikes),
        DEBUG_callSpreadWidth(PDFRequestLNStrike::default_callSpreadWidth),
        DEBUG_accuracy(PDFRequestLNStrike::default_accuracy),
        DEBUG_failProbs(PDFParams::default_failProbs),
        DEBUG_spotVol(ImpliedParams::default_spotVol),
        DEBUG_numNormals(ImpliedParams::default_numNormals),
        isCarefulRandoms(false) {
    }

    MCPathConfigImplied():
        MCPathConfig(TYPE, DependenceMakerGaussTermSP(new DependenceMakerGaussTerm())),
        volType("VolSurface"),
        dependenceType("not used"),
        numMidStrikesPerLogStrike(PDFParams::default_numMidStrikesPerLogStrike),
        minMidStrikes(PDFParams::default_numMiddleStrikes),
        numTailStrikesPerLogStrike(PDFParams::default_numTailStrikesPerLogStrike),
        minTailStrikes(PDFParams::default_numTailStrikes),
        numMidStdDevs(PDFParams::default_numMiddleStdDevs),
        numTailStdDevs(PDFParams::default_numTailStdDevs),
        cacheStrikes(DEFAULT_CACHE_STRIKES),
        DEBUG_propMidStrikes(PDFParams::default_propMiddleStrikes),
        DEBUG_callSpreadWidth(PDFRequestLNStrike::default_callSpreadWidth),
        DEBUG_accuracy(PDFRequestLNStrike::default_accuracy),
        DEBUG_failProbs(PDFParams::default_failProbs),
        DEBUG_spotVol(ImpliedParams::default_spotVol),
        DEBUG_numNormals(ImpliedParams::default_numNormals){}

private:
    static string getDefaultCacheStrikes() {
        return ALWAYS_CACHE_STRIKES;
    }

    // fields ////
    PDFParamsSP      pdfParams;
    ImpliedParamsSP  impliedParams;

    string           volType;
    string           dependenceType;

    int              numMidStrikesPerLogStrike;
    int              minMidStrikes;
    int              numTailStrikesPerLogStrike;
    int              minTailStrikes;
    double           numMidStdDevs;
    double           numTailStdDevs;

    string           cacheStrikes;

    double           DEBUG_propMidStrikes;
    double           DEBUG_callSpreadWidth;
    double           DEBUG_accuracy;
    bool             DEBUG_failProbs;
    bool             DEBUG_spotVol;
    int              DEBUG_numNormals;

    MCImpliedCacheSP                         impliedCache;      // holds information accross tweaks
    bool                                     isCarefulRandoms;
    EquidistantLinearInterpolantNonVirtualSP tabulatedNorm;     // tabulated normal distribution

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
     *  for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const {
        MarketDataFetcherSP mdf(new MarketDataFetcherLNSpline(volType));
        dependenceMaker->modifyMarketDataFetcher(mdf);
        if(skewMaker.get()) {
            skewMaker->modifyMarketDataFetcher(mdf);            
        }
        return mdf;
    }

public:
    virtual MCPathGeneratorSP makePathGenerator(
        bool                     cachingRequested,
        int                      numPaths,
        const MCPathGeneratorSP& pastPathGenerator,
        const IMCProduct*         prod,
        Control*                 control,
        Results*                 results,
        DateTimeArray&           simDates ); // defined below

    /** Creates a past path generator */
    MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod);

    int storagePerPath(IMCProduct* product) const;

        /** Creates a future path generator */
    MCPathGeneratorSP futurePathGenerator(
        int                      cachingMode,
        int                      numPaths,
        const MCPathGeneratorSP& pastPathGenerator,
        const IMCProduct*         prod,
        Control*                 control,
        Results*                 results,
        DateTimeArray&           simDates );

    virtual bool vegaMatrixSupported() const { return false; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * Irrelevant since it's B-S not parametric.  See
     * IModel::wantsRiskMapping().
     */

    IModel::WantsRiskMapping wantsRiskMapping() const {
        return IModel::riskMappingIrrelevant;
    }

    virtual bool carefulRandoms() const { return isCarefulRandoms; }

    virtual void validatePop2Object() {
        static const string routine("MCPathConfigImplied::validatePop2Object");

        try {
            // Get default methodology for caching strikes
            if(CString::equalsIgnoreCase(cacheStrikes, DEFAULT_CACHE_STRIKES)) {
                cacheStrikes = getDefaultCacheStrikes();
            }

            checkCacheStrikesMethod(cacheStrikes);

            // Construct a LinearInterpolant for the standard normal distribution
            LinearInterpolatorSP interpolator = LinearInterpolatorSP(new LinearInterpolator());

            // double normMax = numTailStdDevs + 1.0 / (double)(DefaultImp::default_numFineStdDevs);
            tabulatedNorm = EquidistantLinearInterpolantNonVirtualSP::constCast(
                LinearInterpolator::computeEquidInterpNV(
                *interpolator, -numTailStdDevs, numTailStdDevs, DEBUG_numNormals, N1));

            // Ideally these classes should be an input and the set of doubles, ints etc.
            // should be retired. Cannot do it though because it would require
            // changes in IMS
            pdfParams = PDFParamsSP(new PDFParams(numMidStrikesPerLogStrike,
                                                  numTailStrikesPerLogStrike,
                                                  minMidStrikes,
                                                  minTailStrikes,
                                                  DEBUG_propMidStrikes,
                                                  numMidStdDevs,
                                                  numTailStdDevs,
                                                  PDFParams::default_maxStdDevs,
                                                  PDFParams::default_maxFineStdDevs,
                                                  DEBUG_failProbs));

            impliedParams = ImpliedParamsSP(new ImpliedParams(DEBUG_spotVol,
                                                              DEBUG_numNormals,
                                                              DEBUG_callSpreadWidth,
                                                              DEBUG_accuracy));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }


    /** Does special things for Theta-type tweaks */
    virtual bool sensShift(Theta* shift) {
        MCPathConfig::sensShift(shift); // call parent's method
        if (impliedCache.get()) {
            // We delete the vols and probs for theta tweaks
            impliedCache->killVols();
            impliedCache->killProbs();
        }
        cacheMgr = MCCacheManager();    // Blank cache manager now

        return true;
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigImplied, clazz);
        SUPERCLASS(MCPathConfig);
        EMPTY_SHELL_METHOD(defaultMCPathConfigImplied);
        FIELD(pdfParams, "Parameters for PDF calculation");
        FIELD_MAKE_TRANSIENT(pdfParams);
        FIELD(impliedParams, "Additional parameters for MCImplied");
        FIELD_MAKE_TRANSIENT(impliedParams);
        FIELD(volType, "Type of vol to use");
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);
        FIELD(minMidStrikes,  "Number of percentile strikes in middle");
        FIELD_MAKE_OPTIONAL(minMidStrikes);
        FIELD(numMidStrikesPerLogStrike,  "Number of middle strikes per log strike");
        FIELD_MAKE_OPTIONAL(numMidStrikesPerLogStrike);
        FIELD(minTailStrikes,  "Number of equidistant strikes in tail");
        FIELD_MAKE_OPTIONAL(minTailStrikes);
        FIELD(numTailStrikesPerLogStrike,  "Number of tail strikes per log strike");
        FIELD_MAKE_OPTIONAL(numTailStrikesPerLogStrike);
        FIELD(DEBUG_propMidStrikes,  "Proportion of middle strikes for first pass");
        FIELD_MAKE_OPTIONAL(DEBUG_propMidStrikes);
        FIELD(numTailStdDevs,  "Number of Std Devs");
        FIELD_MAKE_OPTIONAL(numTailStdDevs);
        FIELD(numMidStdDevs,  "Number of Std Devs");
        FIELD_MAKE_OPTIONAL(numMidStdDevs);
        FIELD(DEBUG_callSpreadWidth, "PDF Call spread width");
        FIELD_MAKE_OPTIONAL(DEBUG_callSpreadWidth);
        FIELD(DEBUG_accuracy, "PDF Accuracy");
        FIELD_MAKE_OPTIONAL(DEBUG_accuracy);
        FIELD(DEBUG_failProbs, "Fail at negative probabilities");
        FIELD_MAKE_OPTIONAL(DEBUG_failProbs);
        FIELD(DEBUG_spotVol, "Spot Vol");
        FIELD_MAKE_OPTIONAL(DEBUG_spotVol);
        FIELD(DEBUG_numNormals,  "Size of normal table");
        FIELD_MAKE_OPTIONAL(DEBUG_numNormals);
        FIELD(cacheStrikes,  "Cache strikes between tweaks");
        FIELD_MAKE_OPTIONAL(cacheStrikes);
        FIELD(impliedCache, "MCImplied Cache");
        FIELD_MAKE_TRANSIENT(impliedCache);
        FIELD(isCarefulRandoms, "isCarefulRandoms");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms);
        FIELD(tabulatedNorm, "Tabulated normal distribution");
        FIELD_MAKE_TRANSIENT(tabulatedNorm);

        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigImplied(){
        return new MCPathConfigImplied(TYPE, IVolatilityBS::TYPE->getName(),"not used");
    }
};

CClassConstSP const MCPathConfigImplied::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigImplied", typeid(MCPathConfigImplied), load);

class MCPathBaseImplied: virtual public MCPathBase,
                         virtual public IMCRandom::Callbacks,
                         public DependenceMakerGauss::Support,
                         public DependenceMakerGaussTerm::Support,
                         public DependenceMakerLocalCorr::Support {
public:
    MCPathBaseImplied(MCPathConfigImplied*     pathConfig,
                      const MCPathGeneratorSP& pastPathGenerator,
                      const IMCProduct*         prod,
                      const MCImpliedCacheSP&  impliedCache,
                      Results*                 results);

    /** Required for quickGreeks */
    double maxDriftProduct(int iAsset) const;

    /** Obtains timeline object from base */
    MCProductTimelineConstSP getTimeline() const {
        return timeline;
    }

    MCImpliedCacheSP getCache() {
        return impliedCache;
    }

    /** Returns number of assets */
    int NbSimAssets() const;

    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath);

    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const;

    /** Returns the paths */
    double* Path(int iAsset, int iPath);

    /** Returns the number of paths per asset */
    int nbPaths(int iAsset) const {
        return 1;
    }

    /** Returns if it is a single path generator or not */
    bool isSinglePath() const {
        return true;
    }

    /** Generates an N-Factor MCImplied path */
    virtual void generatePath(int pathIdx, int iAsset, int iPath);

    virtual void drawRandomNumbers(int pathIdx);

    /** Configures pathGen for antithetics */
    virtual void configureAntithetics();

    /** Configures pathGen for nonAntithetics. We need a better name */
    virtual void configureNonAntithetics();

    /** Gets VolRequests, volStartDate etc. */
    void getVolRequests(const IMCProduct*                        prod,
                        const MCPathGeneratorSP&      pastPathGenerator,
                        DateTimeArray&                          cliquetTimeline,
                        vector< vector<VolRequestLNStrikeSP> >& volRequests,
                        int&                                    numCliquets);

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const;

    // implement methods since derivation from DependenceMakerGaussTerm::Support
    virtual DateTimeArray getSimDates() const;
    virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const;
    virtual const IMultiFactors* getMultiFactors() const;

    // implement methods since derivation from DependenceMakerLocalCorr::Support
    virtual DependenceMakerSP getDependenceMaker() const;
    virtual int getNbPastDates() const;
    virtual const IRandomSP getRandomGenerator() const; 
    virtual const MCPathConfig::RandomCacheSP getIdiosynFactorRandomCache() const;
    virtual const MCPathConfig::RandomCacheSP getMarketFactorRandomCache() const;
    virtual bool carefulRandoms() const;    

private:
    MCImpliedCacheSP      impliedCache;        // the MCImplied cache
    Results*              results;             // for writing information to the DEBUG packet

    DoubleMatrix          driverPath;          // [iAsset][iStep]  brownian samples
    DoubleMatrix          cumulatives;         // [iAsset][iStep]  cumulative of brownian samples

    vector<CliquetSP>     cliquetArray;        // information per cliquet
    int                   numCliquets;         // number of cliquets
    IntArray              cliquetDiffOffsets;  // marks the beginning of each cliquet

    IMCRandomSP           randomGen;           // Random number generator
    MCProductTimelineSP   timeline;            // Product timeline
    DateTimeArray         cliquetTimeline;     // = today + liveCliquetDates
    RefLevelDataSP        refData;             // RefLevel data

    const IMultiFactors*  mAsset;              // Multi asset interface
    int                   numAssets;           // number of assets
    vector<DoubleArray>   productPaths;        // [iAsset][iStep]
    EquidistantLinearInterpolantNonVirtualSP tabulatedNorm; //!< Tabulated normal distribution

    // Additional fields used for LocalCorr dependence maker
    DependenceMakerSP           dependenceMaker;
    IRandomSP                   rand;
    MCPathConfig::RandomCacheSP idiosynFactorRandomCache;
    MCPathConfig::RandomCacheSP marketFactorRandomCache;
    bool                        isCarefulRandoms;
};

#ifdef STATE_VARIABLES


/////////////////////////////////////////////////////////////////////////////

/** Referee class that distributes simulation to components
    e.g. spot, hitTime etc. */
class MCPathConfigImplied::Gen: virtual public MCPathGenerator, // For backward compatibility
                                virtual public MCStatelessPathGen,
                                virtual public IStateVariableGen::IStateGen {
public:
    /** Constructor */
    Gen(MCPathConfigImplied*     pathConfig,
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient,
        MCCacheManager&          cacheMgr,
        bool                     cachingRequested,
        const MCImpliedCacheSP&  impliedCache,
        Results*                 results,
        DateTimeArray&           simDates);

    /** MCpathGenerator methods */
    // Deprecated methods
    virtual int NbSimAssets() const;

    virtual const double* Path(int iAsset, int iPath) const;

    virtual double refLevel(int iAsset, int iPath) const;

    virtual double maxDriftProduct(int iAsset) const;

    virtual int begin(int iAsset) const;

    virtual int end(int iAsset) const;

    // Live methods
    virtual bool hasPast() const;

    virtual bool doingPast() const;

    virtual void generatePath(int pathIdx);

    virtual void advance();

    virtual void reset();

    virtual int getPathIndex() const;

    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen);

    MCImpliedCacheSP getCache() const;

private:
    StateVarDBase                      svDBase;        //!< Collection of Generators + statevars
    MCPathConfigImplied::PathGenSpotSP pathGenSpot;    //!< Spot generator
    MCPathConfigImplied::PathGenHVBBSP pathGenHVBB;    //!< Hit value Brownian Bridge generator
    int                                nowPathIdx;     //!< Current path idx
};


/////////////////////////////////////////////////////////////////////////////


/** Spot path generator for MCImplied process using state variable approach */
class MCPathConfigImplied::PathGenSpot: virtual public MCPathGen,
                                        virtual public IMCRandom::Callbacks,
                                        public DependenceMakerGauss::Support,
                                        public DependenceMakerGaussTerm::Support {
public:
    PathGenSpot(MCPathConfigImplied*     pathConfig,
                const PastPathGenSpotSP& pastPathGenSpot,
                const MCProductClient*   prodClient,
                const SVGenSpotArray&       spotGenArray,
                const DateTimeArray&     driverDates,
                MCCacheManager&          cacheMgr,
                bool                     cachingRequested,
                StateVarDBase&           svDBase,
                const MCImpliedCacheSP&  suppliedCache,
                Results*                 results,
                bool                     requestDriver,
                DateTimeArray&           simDates);

    /** Generates an N-Factor MCImplied path */
    virtual void generatePath(int pathIdx);

    virtual void advance();

    virtual void reset();

    /** Returns false always */
    virtual bool doingPast() const;

    /** Configures pathGen for antithetics */
    virtual void configureAntithetics();

    /** Configures pathGen for nonAntithetics. We need a better name */
    virtual void configureNonAntithetics();

    /** Required for quickGreeks */
    double maxDriftProduct(int iAsset) const;

    bool hasPast() const {
        return simHasPast;
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const;

    // implement methods since derivation from DependenceMakerGaussTerm::Support
    virtual DateTimeArray getSimDates() const;
    virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const;
    virtual const IMultiFactors* getMultiFactors() const;

    MCImpliedCacheSP getCache() const {
        return impliedCache;
    }
    /** Information passed to driver users i.e. driver
        values and vols for all simulation dates. */
    class Driver {
    public:
        Driver(const DateTimeArray&     cliquetTimeline,
               const DateTimeArray&     simTimeline,
               const DoubleMatrix&      driverPath,
               const CliquetCacheArray& cliquetCaches,
               const vector<CliquetSP>& cliquets):
        cliquetTimeline(cliquetTimeline), simTimeline(simTimeline), driverPath(driverPath) {
            int numDates = simTimeline.size();
            int numAssets = driverPath.numCols();
            refReturns = vector<vector<const double*> >(numAssets);
            isRelative = vector<vector<int> >(numAssets);

            // Setup dimensions
            driverVols = DoubleMatrix(numAssets, numDates);
            driverVars = DoubleMatrix(numAssets, numDates);
            mappings = vector<ImpliedMappingArray>(numAssets);
            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                mappings[iAsset] = ImpliedMappingArray(numDates);
                int index = 0;
                for(int iCliquet = 0; iCliquet < cliquetCaches.size(); iCliquet++) {
                    const CliquetCache& clCache = *cliquetCaches[iCliquet];
                    int numSteps = (*clCache.assetCaches)[iAsset]->sqrtVars->sqrtFwdVar->size();
                    for(int iStep = 0 ; iStep < numSteps; iStep++) {
                        // Copy over info
                        driverVols[iAsset][index] = (*((*clCache.assetCaches)[iAsset]->sqrtVars->sqrtFwdVar))[iStep];
                        driverVars[iAsset][index] = driverVols[iAsset][index]*driverVols[iAsset][index];
                        ImpliedMappingSP mapping = cliquets[iCliquet]->mappings[iAsset][iStep];
                        mappings[iAsset][index] = mapping;
                        refReturns[iAsset].push_back(cliquets[iCliquet]->returnReference[iAsset]);
                        isRelative[iAsset].push_back(mapping->isRelativeReturn());
                        index++;
                    }
                }
            }
        }

        ~Driver() {}

        /** Mapping function from spot to driver */
        double spotFromDriver(double driver, int iAsset, int iStep) const {
            return mappings[iAsset][iStep]->spotFromDriver(driver);
        }

        /** Mapping function from driver to spot */
        double driverFromSpot(double spot, int iAsset, int iStep) const {
            return mappings[iAsset][iStep]->driverFromSpot(spot);
        }

        const DateTimeArray&           cliquetTimeline;   //!< Today + liveCliquetDates
        const DateTimeArray&           simTimeline;       //!< Timeline of driver
        const DoubleMatrix&            driverPath;        //!< Driver paths
        DoubleMatrix                   driverVols;        //!< Driver vols: concatenates vols from cliquets
        DoubleMatrix                   driverVars;
        vector<ImpliedMappingArray>    mappings;          //!< Implied mapping
        vector<vector<const double*> > refReturns;        //!< Reference level for implied mappings
        vector<vector<int> >           isRelative;        //!< Whether it is a relative return or not
    };
    DECLARE_REF_COUNT(Driver);

    /** Returns a driver object */
    DriverConstSP getDriver() const {
        static const string routine = "MCPathConfigImplied::PathGenSpot::getDriver";

        // Make sure that we are refreshing the driver
        if(!requestDriver) {
            throw ModelException(routine, "Internal error: driver has not been requested");
        }

        return DriverConstSP(new Driver(cliquetTimeline, simTimeline, driverPaths,
            *impliedCache->cliquetCaches, cliquetArray));
    }

private:
    /** Generates driver and cumulatives */
    void generateDriverAndCumulatives();

    /** Gets VolRequests, volStartDate etc. */
    void getVolRequests(const MCProductClient*                  prodClient,
                        const DateTime&                         maturity,
                        DateTimeArray&                          cliquetTimeline,
                        vector< vector<VolRequestLNStrikeSP> >& volRequests,
                        int&                                    numCliquets);

    /** Generate path for specified path in simulation (pathIdx), for
        specified asset (iAsset) */
    void generatePath(int pathIdx, int iAsset);

    bool                  requestDriver;        //!< Whether to refresh driver all the time or try to save time
    MCImpliedCacheSP      impliedCache;         //!< The MCImplied cache
    const IMultiFactors*  mAsset;               //!< Multi asset interface
    int                   numAssets;            //!< Number of assets

    DateTimeArray         simTimeline;          //!< Today + strictly future merged dates
    DateTimeArray         cliquetTimeline;      //!< Today + liveCliquetDates
    bool                  simHasPast;           //!< Whether product has past

    DoubleArray           spotAtStart;          //!< Spot at simulation start

    vector<DoubleArray>   spotProdPaths;        //!< Spot price at asset specific dates [iAsset][iStep]
    DoubleMatrix          driverPaths;          //!< Driver        at simTimeline [iAsset][iStep]
    DoubleMatrix          cumulatives;          //!< Uniforms      at simTimeline [iAsset][iStep]
    vector<IntArray>      doSpot;               //!< Simulate spot at simTimeline [iAsset][iStep]
    IntArray              spotProdOffset;       //!< Starting point for future path per asset

    MCRandomGenCacheSP    randomGen;            //!< Random number generator for normals
    DependenceSP          dependence;           //!< Dependence object e.g. gaussian
    SVGenSpotPathCacheSP     spotPathCache;        //!< Cache of generated paths

    vector<CliquetSP>     cliquetArray;         //!< Information per cliquet
    int                   numCliquets;          //!< Number of cliquets
    IntArray              cliquetDiffOffsets;   //!< Marks the beginning of each cliquet

    EquidistantLinearInterpolantNonVirtualSP tabulatedNorm; //!< Tabulated normal distribution

    static const MCCache::KeyUtils::Key keyPathCache;  //!< Key for spot path cache
    static const MCCache::KeyUtils::Key keyRandCache;  //!< Key for spot path random numbers cache

    vector<IAdvanceableStateVariable*> advanceSvSet; // for stateless advance
};

// Initialize the keys to the caches
const MCCache::KeyUtils::Key MCPathConfigImplied::PathGenSpot::keyPathCache = MCCache::KeyUtils::getNewKey();
const MCCache::KeyUtils::Key MCPathConfigImplied::PathGenSpot::keyRandCache = MCCache::KeyUtils::getNewKey();


/////////////////////////////////////////////////////////////////////////////


MCPathConfigImplied::PathGenSpot::PathGenSpot(
    MCPathConfigImplied*     pathConfig,
    const PastPathGenSpotSP& pastPathGenSpot,
    const MCProductClient*   prodClient,
    const SVGenSpotArray&       spotGenArray,
    const DateTimeArray&     driverDates,
    MCCacheManager&          cacheMgr,
    bool                     cachingRequested,
    StateVarDBase&           svDBase,
    const MCImpliedCacheSP&  suppliedCache,
    Results*                 results,
    bool                     requestDriver,
    DateTimeArray&           simDates):
requestDriver(requestDriver),
impliedCache(suppliedCache), // might be null
mAsset(prodClient->getMultiFactors()),
numAssets(prodClient->getNumAssets()),
simHasPast(pastPathGenSpot->hasPast()),
spotAtStart(numAssets),
spotProdPaths(prodClient->getNumAssets()),
doSpot(prodClient->getNumAssets()),
spotProdOffset(prodClient->getNumAssets()),
tabulatedNorm(pathConfig->tabulatedNorm) {

    static const string routine("PathGenSpot::PathGenSpot");
    try{
        const DateTime& today = prodClient->getToday();
        const IPastValues* pastValues = prodClient->getMCPastValues();

        // Get driver and spot dates
        for(int iDate = 0; iDate < driverDates.size(); iDate++) {
            if(today > driverDates[iDate]) {
                throw ModelException("Driver dates must be strictly in the future.");
            }
        }
        DateTimeArraySP spotGenDates = MCPath::getAllDates(spotGenArray);

        // Enrich driver dates by spot dates and drop past dates
        DateTimeArray allDriverDates = DateTime::merge(driverDates, *spotGenDates);
        DateTimeArray futDriverDates = today.getFutureDates(allDriverDates);

        // Create simulation timeline
        simTimeline = futDriverDates;
        simTimeline.insert(simTimeline.begin(), today);
        int nbRandoms = simTimeline.size() - 1;

        if (dynamic_cast<const IMCStatelessProductClient*>(prodClient))
            simDates = simTimeline;

        // Allocate memory for matrices now that timeline is known
        driverPaths = DoubleMatrix(numAssets, simTimeline.size());
        cumulatives = DoubleMatrix(numAssets, simTimeline.size());

        // At this stage, future spot dates are a subset of driver dates.
        // So we will always simulate the driver and sometimes the spot.

        // Create spot paths and driver paths per asset
        DateTimeArrayArray spotDatesPerAsset(numAssets);
        vector<const double*> spotPtrs(numAssets);
        vector<int> spotBeginInd(numAssets);
        vector<int> spotEndInd(numAssets);

        int iAsset;
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            // 1) SPOT PATHS
            // Populate past values
            DateTimeArraySP spotAssetDates =
                MCPath::getAllDates(spotGenArray, iAsset);
            spotProdPaths[iAsset]   = DoubleArray(spotAssetDates->size());
            DateTimeArray assetPastDates = today.getPastDates(*spotAssetDates);
            DoubleArray assetPastValues =
                pastValues->getPastValues(assetPastDates, iAsset, today);
            spotProdOffset[iAsset] = assetPastValues.size();
            int iStep;
            for(iStep = 0; iStep < assetPastValues.size(); iStep++) {
                spotProdPaths[iAsset][iStep] = assetPastValues[iStep];
            }

            // Create spot mappings
            spotDatesPerAsset[iAsset] = *spotAssetDates;
            spotPtrs[iAsset] = &spotProdPaths[iAsset][0];
            spotBeginInd[iAsset] = assetPastDates.size();
            spotEndInd[iAsset]   = spotAssetDates->size();

            // 2) SPOT AND DRIVER FLAGS
            IntArray& doSpotPerAsset = doSpot[iAsset];
            doSpotPerAsset = IntArray(simTimeline.size(), 0);
            // Be careful not to use iStep = 0 because today might be
            // an asset date but it has been dealt with in the past
            for(iStep = 1; iStep < simTimeline.size(); iStep++) {
                int occurances = count(spotAssetDates->begin(),
                                       spotAssetDates->end(),
                                       simTimeline[iStep]);
                doSpotPerAsset[iStep] = occurances > 0 ? 1: 0;
            }
        }

#ifdef STATEVAR_CACHING
        // Caching of spot paths
        IntArray numDatesPerAsset(numAssets, 0);
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            numDatesPerAsset[iAsset] = spotEndInd[iAsset] - spotBeginInd[iAsset];
        }
        spotPathCache = SVGenSpotPathCache::createCache(
            cacheMgr, keyPathCache, numDatesPerAsset, cachingRequested);

        // Caching of ramdom numbers
        MCRandomCacheSP randCache = MCRandomCache::createCache(cacheMgr,
            keyRandCache, numAssets, nbRandoms, cachingRequested);
#endif

        vector<double> maxDrifts(numAssets);
        // Create SVGenSpot::IStateVars and put them in database
        MCPath::IStateVarArray spotSVArray(MCPath::createPaths(
            false,
            spotGenArray,
            spotDatesPerAsset,
            spotBeginInd,
            spotEndInd,
            spotPtrs,
            maxDrifts));

        unsigned int iVar;
        for(iVar = 0; iVar < spotGenArray.size(); iVar++) {
            svDBase.append(spotGenArray[iVar], spotSVArray[iVar]);
        }

        // populate advanceables
        svDBase.filterAdvanceables(advanceSvSet);

        // Deduce the cliquet timeline and obtain volRequests for each cliquet
        vector< vector<VolRequestLNStrikeSP> > volRequests; //[iCliquet][iAsset]

        getVolRequests(prodClient,
                       simTimeline.back(),
                       cliquetTimeline,
                       volRequests,
                       numCliquets);

        // Categorize the simulation dates on a per cliquet basis
        int iCliquet, iStep;
        vector<DateTimeArray> datesTo(numCliquets);

        for(iStep = 0; iStep < simTimeline.size(); iStep++) {
            const DateTime& date = simTimeline[iStep];
            // figure out on which cliquet the date belongs
            for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                DateTime& startDate = cliquetTimeline[iCliquet];
                DateTime& endDate   = cliquetTimeline[iCliquet+1];
                if(startDate <= date && date <= endDate) {
                    // record the date
                    datesTo[iCliquet].push_back(date);
                    break;
                }
            }
        }

        // Record indices of the timeline at the beginning of each cliquet
        // Required for simulating paths by looping over cliquets
        cliquetDiffOffsets = IntArray(numCliquets);
        for(iCliquet = 1; iCliquet < numCliquets; iCliquet++) {
            cliquetDiffOffsets[iCliquet] = datesTo[iCliquet-1].size();
        }

        // initialize the cache, unless it exists from some previous tweak
        if(!impliedCache) {
            impliedCache = MCImpliedCacheSP(
                new MCImpliedCache(numCliquets));
        }

        // Kill the cache if we have rolled over a cliquet
        if(numCliquets != impliedCache->numCliquets) {
            impliedCache->killCliquetCaches(numCliquets);
        }

        // construct array of cliquets
        cliquetArray = vector<CliquetSP>(numCliquets);
        for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            try {
                const DateTime& cliqStartDate = cliquetTimeline[iCliquet];
                CliquetCacheSP& cliquetCache = (*impliedCache->cliquetCaches)[iCliquet];
                CliquetSP cliquet;
                if(iCliquet == 0) {
                    // For 1st cliquet invoke constructor that will use
                    // deterministic reference
                    cliquet = CliquetSP(
                        new Cliquet(iCliquet,
                                    volRequests[iCliquet],
                                    mAsset,
                                    today,
                                    cliqStartDate,
                                    datesTo[iCliquet],
                                    pathConfig->pdfParams,
                                    pathConfig->impliedParams,
                                    cliquetCache,
                                    pathConfig->cacheStrikes,
                                    tabulatedNorm));
                } else {
                    // Create stochastic reference levels for subsequent cliquets.
                    // This is a pointer to the spot path at cliqStartDate
                    vector<const double*> stochasticReference;
                    for(iAsset = 0; iAsset < numAssets; iAsset++) {
                        const DateTimeArray& spotAssetDates = spotDatesPerAsset[iAsset];
                        bool found = false;
                        for(iStep = 0; iStep < spotAssetDates.size(); iStep++) {
                            if(cliqStartDate == spotAssetDates[iStep]) {
                                const double* ref = &spotProdPaths[iAsset][iStep];
                                stochasticReference.push_back(ref);
                                found = true;
                                break;
                            }
                        }
                        if(!found) {
                            throw ModelException(
                                "Spot price at cliquet start date " +
                                cliqStartDate.toString() +
                                " has not been requested for asset " +
                                Format::toString(iAsset+1) +
                                ". All assets must request spot path at all cliquet start dates.");
                        }
                    }
                    cliquet = CliquetSP(
                        new Cliquet(iCliquet,
                                    volRequests[iCliquet],
                                    mAsset,
                                    today,
                                    cliqStartDate,
                                    datesTo[iCliquet],
                                    pathConfig->pdfParams,
                                    pathConfig->impliedParams,
                                    cliquetCache,
                                    pathConfig->cacheStrikes,
                                    tabulatedNorm,
                                    stochasticReference));
                }
                cliquetArray[iCliquet] = cliquet;

            } catch(exception& e) {
                string message = "Failed to construct cliquet number " +
                    Format::toString(iCliquet + 1);
                throw ModelException::addTextToException(e, message);
            }
        }

        // set up dependence
        dependence = pathConfig->dependenceMaker->createDependence(this);

        // Initialize random number generator
        randomGen = MCRandomGenCacheSP(new MCRandomGenCache(
            this,
            dependence,
            pathConfig->getRandomGenerator(),
            randCache,
            pathConfig->carefulRandoms(),
            nbRandoms,
            numAssets,
            allDriverDates.size() - futDriverDates.size()));

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/** Generates an N-Factor MCImplied path */
void MCPathConfigImplied::PathGenSpot::generatePath(int pathIdx) {
    // Draw random numbers
    randomGen->generate(pathIdx);

    for (int iAsset = 0; iAsset < numAssets; iAsset++) {
#ifdef STATEVAR_CACHING
        int offset = spotProdOffset[iAsset];
        // int offset = commonDates->size() - numFutSteps;
        double* futPath = &spotProdPaths[iAsset][offset];
        if (spotPathCache->isValid(iAsset)){
            // Read paths from cache
            spotPathCache->read(iAsset, pathIdx, futPath);
        } else {
            // Diffuse paths and write them to cache
            generatePath(pathIdx, iAsset);
            if (spotPathCache->updateAllowed(iAsset)) {
                spotPathCache->write(iAsset, pathIdx, futPath);
            }
        }
#else
        generatePath(pathIdx, iAsset);
#endif
    }

    reset();

}

// for stateless payoff
void MCPathConfigImplied::PathGenSpot::reset() {
    for (size_t i = 0; i < advanceSvSet.size(); ++i)
        advanceSvSet[i]->reset();
}

void MCPathConfigImplied::PathGenSpot::advance() {
    for (size_t i = 0; i < advanceSvSet.size(); ++i)
        advanceSvSet[i]->advance();
}

/** Returns false always */
bool MCPathConfigImplied::PathGenSpot::doingPast() const {
    return false;
}


void MCPathConfigImplied::PathGenSpot::configureAntithetics() {
    // Notice that there is some overhead as the MCRandom class has
    // already negated the random numbers

    // In the longer run, we might want to disable fast antithetics
    // because we are losing the driver values. We will need them for BB.

    // Strictly speaking, this pathGen should not be using an MCRandom
    // random number generator but a MCRandomUniform (to be coded).
    // The latter should do antithetics by u(anti) = 1.0 - u
    if(requestDriver) {
        // Generate drivers and cumulatives
        generateDriverAndCumulatives();
    } else {
        // Just do u(antithetic) = 1.0 - u
        cumulatives.negate();
        cumulatives.scalarAdd(1.0);
    }
}

// implement getGaussData since derivation from DependenceMakerGauss::Support
CDoubleMatrixConstSP MCPathConfigImplied::PathGenSpot::getGaussData() const {
    static const string method("MCPathConfigImplied::PathGenSpot::getGaussData");
    try {
        return mAsset->factorsCorrelationMatrix();
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

// implement getGaussTermData since derivation from DependenceMakerGaussTerm::Support
DateTimeArray MCPathConfigImplied::PathGenSpot::getSimDates() const {
    return simTimeline;
}
DoubleMatrix MCPathConfigImplied::PathGenSpot::getFwdVarAtDates(bool interpolateAtmFwd) const {
    // SV approach is different from nonSV approach: today needs to be passed to Driver
    const CliquetCacheArray* tmpCCA = impliedCache->cliquetCaches.get();
    MCPathConfigImplied::PathGenSpot::Driver driver(cliquetTimeline,
                                                    simTimeline,
                                                    driverPaths,
                                                    *tmpCCA,
                                                    cliquetArray);
    int nbDates = simTimeline.size()-1;
    int nbAssets = driver.driverVars.numCols();
    DoubleMatrix fwdVarAtDates(nbAssets,nbDates);
    for (int iDate=0; iDate<nbDates; iDate++) {
        for (int iAsset=0; iAsset<nbAssets; iAsset++) {
            fwdVarAtDates[iAsset][iDate] = driver.driverVars[iAsset][iDate+1];
        }
    }
    return fwdVarAtDates;
}
const IMultiFactors* MCPathConfigImplied::PathGenSpot::getMultiFactors() const {
    return mAsset;
}

void MCPathConfigImplied::PathGenSpot::configureNonAntithetics() {
    // Generate drivers and cumulatives
    generateDriverAndCumulatives();
}


/** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
    the 'drift' between simulation date i and simulation date i+1.
    (With d_0 = simulation start date). Range of dates is restricted to
    those in the future so the past path generator will
    return 1.0. The 'drift' is the biggest ratio of the 2 spots possible
    within some sort of probability interval. See  MCPathConfigLN */
double MCPathConfigImplied::PathGenSpot::maxDriftProduct(int iAsset) const {
    static const string routine("PathGenSpot::maxDriftProduct");

    try {
        // Return the maximum maxDrift over all possible cliquets
        double maxDrift = 0.0;
        for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            if(cliquetArray[iCliquet]->maxDrifts.get()) {
                maxDrift = Maths::max(maxDrift,
                                      cliquetArray[iCliquet]->maxDrifts->getMaxDrift(iAsset));
            } else {
                throw ModelException("Unable to get maxDrifts for cliquet number "+
                                     Format::toString(iCliquet+1) +
                                     " for asset " +
                                     mAsset->assetGetName(iAsset));
            }
        }
        return maxDrift;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void MCPathConfigImplied::PathGenSpot::generateDriverAndCumulatives() {
    // Given the random numbers, compute the driver path
    // and cumulative distributions
    const DoubleMatrix& randoms = randomGen->getRandomNumbers();
    for (int iAsset = 0; iAsset < numAssets; iAsset++) {
        // for speed
        double* assetDriver        = driverPaths[iAsset];
        double* assetCumul         = cumulatives[iAsset];
        const double* assetRandoms = randoms[iAsset];

        for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            // Skip today in 1st cliquet
            int skip = (iCliquet == 0) ? 1 : 0;
            int numFutSteps = cliquetArray[iCliquet]->numFutSteps;

            // For easier access to caches
            const DriversSqrtVarSP& sqrtVars =
                (*(*impliedCache->cliquetCaches)[iCliquet]->assetCaches)[iAsset]->sqrtVars;

            // Some pointer arithmetic to ensure we fill in the path properly
            int cliquetOffset = cliquetDiffOffsets[iCliquet];
            assetDriver  += cliquetOffset + skip;   // Don't update driver at Today
            assetCumul   += cliquetOffset + skip;   // Don't update cumulative at Today
            assetRandoms += cliquetOffset;          // Note that we never had a random for Today so don't offset

            // Extract driver information for the cliquet
            const double* sqrtTotalVar = &(*sqrtVars->sqrtTotalVar)[0] + skip;
            const double* sqrtFwdVar   = &(*sqrtVars->sqrtFwdVar)[0] + skip;
            const int* zeroTT          = &(*sqrtVars->zeroTTime)[0] + skip;

            for(int iStep = 0; iStep < numFutSteps - skip; iStep++) {
                if(zeroTT[iStep]) {
                    // Zero trading time up to here so the Brownian motion is still at 0.0
                    assetCumul[iStep] = 0.5;
                } else {
                    // Diffuse driver and compute cumulatives
                    if(iStep == 0) {
                        assetDriver[iStep] = 0.0;
                    } else {
                        assetDriver[iStep] = assetDriver[iStep - 1];
                    }
                    assetDriver[iStep] += sqrtFwdVar[iStep] * assetRandoms[iStep];
                    assetCumul[iStep] = tabulatedNorm->value(assetDriver[iStep] /
                                                             sqrtTotalVar[iStep]);
                }
            }
            if(skip) {
                assetDriver  -= skip;
                assetCumul   -= skip;
                assetRandoms -= skip;
            }
        }
    }
}


void MCPathConfigImplied::PathGenSpot::getVolRequests(
    const MCProductClient*                  prodClient,
    const DateTime&                         maturity,
    DateTimeArray&                          cliquetTimeline,
    vector< vector<VolRequestLNStrikeSP> >& volRequests,
    int&                                    numCliquets) {
    static const string routine("PathGenSpot::GetVolRequests");
    try{
        // get vol interp from product at this stage. Later on we would like
        // to get the VolRequest from somewhere else. The user should, ideally,
        // pass a VolRequest to MCImplied.
        const IMCProductLN* prodLN = dynamic_cast<const IMCProductLN*>(prodClient);
        if (!prodLN){
            throw ModelException("Product does not implement "
                                 "IMCProductLN interface");
        }

        // Get from the IMCProductLN the first VolRequest for the first asset
        int iAsset = 0;
        CVolRequestLNArray prodVolRequests(prodLN->getVolInterp(0, iAsset)); // pass null past path gen
        if (prodVolRequests.empty()){
            throw ModelException(routine, "No vol requests specified");
        }
        CVolRequestLNSP volReq(prodVolRequests[0]);

        // Deduce if it a cliquetVolRequest or VolRequestLNStrike
        VolRequestLNStrike* volReqLNStr = dynamic_cast<VolRequestLNStrike*>(volReq.get());
        CliquetVolRequest*  cliquetVol  = dynamic_cast<CliquetVolRequest*>(volReq.get());

        bool isCliquet;
        if(volReqLNStr) {
            isCliquet = false;
        } else if(cliquetVol) {
            isCliquet = true;
        } else {
            throw ModelException(routine, "Vol interp is neither LNStrike "
                                 "type nor Cliquet type for asset " +
                                 mAsset->assetGetName(iAsset));
        }

        // The cliquet TimeLine is the intersection of:
        // i) StartDay (if not cliquet) or cliquet dates (if cliquet)
        // ii) last simulation date
        vector<DateTimeArray> dates(2);
        DateTime startDate = prodClient->getToday();
        int iCliquet;
        if(isCliquet) {
            dates[0] = cliquetVol->getCliqStartDates();     // additional cliquet dates
            // the first live cliquet date might be in the past (started product)
            // or in the future (1st cliquet AvgIn). To make consistent with
            // non-cliquet we set it to startDate
            dates[0][0] = startDate;
        } else {
            dates[0] = DateTimeArray(1, startDate);         // product's start day
        }
        DateTimeArray mat(1);
        mat[0] = maturity;
        dates[1] = DateTimeArray(mat);                      // expiry

        cliquetTimeline = DateTime::merge(dates);
        numCliquets = cliquetTimeline.size() - 1;
        if(!numCliquets) {
            // Something has gone wrong. It cannot be that we have zero number of cliquets
            // Maybe we are being overcautious here and the engine can cope with 1 AvgIn date
            // being equal to 1 AvgOut date. For the moment be defensive
            throw ModelException(
                "Simulation start date " + startDate.toString() +
                " and simulation end date " + maturity.toString() +
                " coincide. They have to be distinct.");
        }

        for(iCliquet = 1; iCliquet < numCliquets; iCliquet++) {
            if(dates[0][iCliquet] <= startDate) {
                throw ModelException("Failed to construct the cliquet timeline. "
                                     "More than 1 live cliquet start dates are before the start date.");
            }
        }

        // construct VolRequests per cliquet, asset
        volRequests = vector< vector<VolRequestLNStrikeSP> >(numCliquets);
        for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            volRequests[iCliquet] = vector<VolRequestLNStrikeSP>(numAssets);
        }

        // Get VolRequests or create them
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            CVolRequestLNArray prodVolRequests(prodLN->getVolInterp(0, iAsset));    // pass null past path gen
            if (prodVolRequests.empty()){
                throw ModelException(routine, "No vol requests specified");
            }
            CVolRequestLNSP volReq(prodVolRequests[0]);
            CliquetVolRequest*  cliquetVol  = dynamic_cast<CliquetVolRequest*>(volReq.get());
            VolRequestLNStrike* volReqLNStr = dynamic_cast<VolRequestLNStrike*>(volReq.get());

            if(cliquetVol) {
                CVolRequestLNArray productVolRequests = cliquetVol->getRequestsArray();
                for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                    VolRequestLNStrike* volLNStrike = dynamic_cast<VolRequestLNStrike*>(productVolRequests[iCliquet].get());
                    if(volLNStrike) {
                        volRequests[iCliquet][iAsset] = VolRequestLNStrikeSP(copy(volLNStrike));
                    } else {
                        throw ModelException(routine, "Vol interp is not LNStrike "
                                             "type for asset " +
                                             mAsset->assetGetName(iAsset));
                    }
                }
            } else if(volReqLNStr) {
                volRequests[0][iAsset] = VolRequestLNStrikeSP(copy(volReqLNStr));
            } else {
                throw ModelException(routine, "Vol interp is neither LNStrike "
                                     "type nor Cliquet type for asset " +
                                     mAsset->assetGetName(iAsset));
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/** Generate path for specified path in simulation (pathIdx), for
    specified asset (iAsset) and specified vol interp (iPath) */
void MCPathConfigImplied::PathGenSpot::generatePath(int pathIdx, int iAsset) {
    // populate paths field
    const double* assetCumul  = cumulatives[iAsset];
    const int* doSpotPerAsset = &doSpot[iAsset][0];

    // The path in particular must be offset for the past for
    // the first cliquet. The remaining offsets are relative
    // so we should only do it once for the first cliquet
    int spotOffset    = spotProdOffset[iAsset];
    double* assetPath = &spotProdPaths[iAsset][spotOffset];     // for ease/speed

    int modStep = 0;
    for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
        // extract information for the cliquet
        const Cliquet& cliquet              = *cliquetArray[iCliquet];
        int numFutSteps                     = cliquet.numFutSteps;
#ifdef USE_MAPPINGS
        const ImpliedMappingArray& mappings = cliquet.mappings[iAsset];
#else
        const ImpliedSampleArray& samples   = cliquet.samples[iAsset];
#endif

        // Some pointer arithmetic to ensure we fill in the path properly
        int cliquetOffset = cliquetDiffOffsets[iCliquet];
        assetCumul     += cliquetOffset;
        doSpotPerAsset += cliquetOffset;

        // Get the reference level in advance
        double assetRefLevel = *(cliquet.returnReference[iAsset]);

        // Loop over sample dates and get a sample
        const int* isRelative = &cliquet.isRelativeReturn[0];

        for(int iStep = 0; iStep < numFutSteps; iStep++) {
            if(doSpotPerAsset[iStep]) {
                // Sample and multiply by reference if needed
#ifdef USE_MAPPINGS
                // Use the mapping
                double assetSpot = mappings[iStep]->spotFromDriver(driverPaths[iAsset][iStep]);
#else
                // Use the sample
                double assetSpot = samples[iStep]->sample(assetCumul[iStep]);
#endif
                if(isRelative[iStep]) {
                    assetSpot *= assetRefLevel;
                }
                assetPath[modStep] = assetSpot;
                modStep++;
            }
        }
    }
}


/////////////////////////////////////////////////////////////////////////////


/** Hitting values path generator for LogNormal process
    using state variable approach */
class MCPathConfigImplied::PathGenHVBB: virtual public MCPathGen,
                                        public DependenceMakerGauss::Support,
                                        public DependenceMakerGaussTerm::Support {
public:
    /** Constructor */
    PathGenHVBB(MCPathConfigImplied*         pathConfig,
                const MCProductClient*       prodClient,
                IStateVariableGen::IStateGen*    pastPathGen,
                IStateVariableGen::IStateGen*    futPathGen,
                const PathGenSpotSP&         pathGenSpot,
                const SVGenBarrierHVBBArray& barrierHVBBGens,
                MCCacheManager&              cacheMgr,
                bool                         cachingRequested,
                StateVarDBase&               svDBase,
                EquidistantLinearInterpolantNonVirtualSP tabulatedNorm):
    numBarriers(barrierHVBBGens.size()),
    numAssets(prodClient->getNumAssets()),
    mAsset(prodClient->getMultiFactors()),
    hitNoHitPath(numAssets),
    doHitNoHit(numAssets),
    pathGenSpot(pathGenSpot),
    tabulatedNorm(tabulatedNorm) {

        static const string routine = "MCPathConfigImplied::PathGenHVBB::PathGenHVBB";

        try {
            // Ensure that past and future spot path generators exist
            if(!pathGenSpot) {
                throw ModelException("Spot path generator not initialized. Unable "
                                     "to simulate hitting values.");
            }
            if(!barrierHVBBGens.size()) {
                throw ModelException("No barriers to simulate. Internal error.");
            }

            // Validate the object and assume for the moment that barriers are disjoint
            validate(barrierHVBBGens);

            // Obtain the driver from the path gen
            driver = pathGenSpot->getDriver();

            // Join all barriers together in a common timeline for BBs.
            joinBarriers(pathGenSpot, pastPathGen, futPathGen, barrierHVBBGens, driver, svDBase);

#ifdef STATEVAR_CACHING
            // Caching of HitNoHit values
            IntArray numDatesPerAsset(numAssets, 0);
            for (int iAsset=0; iAsset < numAssets; iAsset++) {
                numDatesPerAsset[iAsset] = hitNoHitPath[iAsset].size();
            }
            hitValuePathCache = MCHitValuePathCache::createCache(
                cacheMgr, keyHitValueCache, numDatesPerAsset, cachingRequested);

            // Caching of ramdom numbers
            MCRandomCacheSP randCache = MCRandomCache::createCache(cacheMgr,
                keyRandCache, numAssets, numBBs, cachingRequested);
#endif

            // set up dependence at driver dates
            dependence = pathConfig->dependenceMaker->createDependence(this);

            // Get random number generator
            IRandomSP rand = pathConfig->getRandomGenerator();
            randomGen = MCRandomGenCacheSP(new MCRandomGenCache(
                0,
                dependence,
                rand,
                randCache,
                pathConfig->carefulRandoms(),
                numBBs,
                numAssets,
                0));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }


    virtual void generatePath(int pathIdx) {
        // Draw random numbers for all assets
        randomGen->generate(pathIdx);

        // Generate path per asset. Need to read / write from cache here
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
#ifdef STATEVAR_CACHING
            int* futPath = &hitNoHitPath[iAsset][0];
            if (hitValuePathCache->isValid(iAsset)){
                // Read paths from cache
                hitValuePathCache->read(iAsset, pathIdx, futPath);
            } else {
                // Diffuse paths and write them to cache
                generatePath(pathIdx, iAsset);
                if (hitValuePathCache->updateAllowed(iAsset)) {
                    hitValuePathCache->write(iAsset, pathIdx, futPath);
                }
            }
#else
            generatePath(pathIdx, iAsset);
#endif
        }
    }

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool doingPast() const {
        return false;
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigImplied::PathGenHVBB::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // implement methods since derivation from DependenceMakerGaussTerm::Support
    virtual DateTimeArray getSimDates() const {
        return driver->simTimeline;
    }
    virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const {
        // SV approach is different from nonSV approach: today is incl in driver data
        // driver is always available (created in constructor of PathGenHVBB)
        int nbDates = driver->simTimeline.size()-1;
        int nbAssets = driver->driverVars.numCols();
        DoubleMatrix fwdVarAtDates(nbAssets,nbDates);
        for (int iDate=0; iDate<nbDates; iDate++) {
            for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                fwdVarAtDates[iAsset][iDate] = driver->driverVars[iAsset][iDate+1];
            }
        }
        return fwdVarAtDates;
    }
    virtual const IMultiFactors* getMultiFactors() const {
        return mAsset;
    }

private:
    /** Simulates brownian bridges for specified asset. */
    void generatePath(int pathIdx, int iAsset) {
        // Obtain random numbers
        const DoubleMatrix& randoms = randomGen->getRandomNumbers();

        const double* assetRandoms              = randoms[iAsset];
        const double* assetDriver               = driver->driverPath[iAsset];
        int* assetHitNoHit                      = &hitNoHitPath[iAsset][0];
        const BoolArray& doHitNoHitAsset        = doHitNoHit[iAsset];

        // Loop on the path and do BBs
        for(int iStep = 0, modStep = 0; iStep < maxNumBBs; iStep++) {
            if(doHitNoHitAsset[iStep]) {
                const BridgeData::AssetBridgeData& assetBBData = *bridgeData[iStep]->assetData[iAsset];

                double barrierStart, barrierEnd;
                assetBBData.getBarrierLevels(barrierStart, barrierEnd);

                double driverStart = assetDriver[iStep];
                double driverEnd   = assetDriver[iStep + 1];

                // Convert problem to discontinuous process against
                // continuous barrier
                barrierEnd -= assetBBData.sumDivs;
                driverEnd  -= assetBBData.sumDivs;

                // Update BB, obtain uniform and get a sample
                HitNoHitBB& bb = *assetBBData.bbSample;
                double uniform = tabulatedNorm->value(assetRandoms[modStep]);
                assetHitNoHit[iStep] = bb.updateAndSample(barrierStart,
                                barrierEnd,
                                driverStart,
                                driverEnd,
                                uniform);
                modStep++;
            }
        }
    }

    /** Basic validation of the object */
    void validate(const SVGenBarrierHVBBArray& barrierHVBBGens) {
        static const string method = "MCPathConfigImplied::PathGenHVBB::validate";

        try {
            // Get the barrier data
            vector<BarrierDataConstSP> barriers;
            unsigned int iBarrier;
            for(iBarrier = 0; iBarrier < barrierHVBBGens.size(); iBarrier++) {
                barriers.push_back(barrierHVBBGens[iBarrier]->getBarrierData());
            }

            // First implementation: only allow disjoint regions for each asset
            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                DateTimeArray startDates, endDates;
                for(iBarrier = 0; iBarrier < barriers.size(); iBarrier++) {
                    const DateTimeArray& monDates = *barriers[iBarrier]->getAssetBarrier(iAsset)->monitoringDates;
                    // Start and end dates for this barrier for the corresponding asset
                    startDates.push_back(monDates.front());
                    endDates.push_back(monDates.back());
                }

                // Validate that startDates and EndDates of other barriers are mutually disjoint
                // for this asset
                for(iBarrier = 0; iBarrier < barriers.size(); iBarrier++) {
                    for(unsigned int jBarrier = iBarrier + 1; jBarrier < barriers.size(); jBarrier++) {
                        bool disjoint =
                            (endDates[iBarrier] < startDates[jBarrier]) ||
                            (endDates[jBarrier] < startDates[iBarrier]);
                        if(! disjoint) {
                            throw ModelException("Barrier objects must be disjoint.");
                        }
                    }
                }
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Create common timeline from N barriers */
    void joinBarriers(const PathGenSpotSP&         pathGenSpot,
                      IStateVariableGen::IStateGen*    pastPathGen,
                      IStateVariableGen::IStateGen*    futPathGen,
                      const SVGenBarrierHVBBArray& barrierHVBBGens,
                      PathGenSpot::DriverConstSP   driver,
                      StateVarDBase&               svDBase) {

        static const string method = "MCPathConfigImplied::PathGenHVBB::joinBarriers";

        try {
            int iAsset, iStep, iBarrier;

            // Offsets for reading hitNoHitPath in each SV
            vector<vector<IntArray> > offsets(numBarriers);
            for(iBarrier = 0; iBarrier < numBarriers; iBarrier++) {
                offsets[iBarrier] = vector<IntArray>(numAssets);
                for(iAsset = 0; iAsset < numAssets; iAsset++) {
                    offsets[iBarrier][iAsset] = IntArray(0);
                }
            }

            // Flag dates at which we need samples
            maxNumBBs = driver->simTimeline.size() - 1;
            if(maxNumBBs < 1) {
                throw ModelException("Need at least 2 dates.");
            }

            // Decide on number of divs to be used in the simulation
            int maxNumDivsPerYear = barrierHVBBGens[0]->getBarrierData()->maxNumDivsPerYear;
            const DateTime& valueDate = barrierHVBBGens[0]->getBarrierData()->valueDate;
            for(iBarrier = 1; iBarrier < numBarriers; iBarrier++) {
                if(maxNumDivsPerYear != barrierHVBBGens[iBarrier]->getBarrierData()->maxNumDivsPerYear) {
                    throw ModelException(
                        "Maximum number of divs per years has to be common between barriers."
                        "Check barrier objects number " +
                        Format::toString(iBarrier + 1) + " and " +
                        Format::toString(iBarrier));
                }
            }

            const DateTimeArray& simTimeline = driver->simTimeline;
            const DateTime& simStart = simTimeline.front();
            const DateTime& simEnd   = simTimeline.back();
            int totalNbDivs = maxNumDivsPerYear == -1 ? maxNumDivsPerYear :
                (int)(maxNumDivsPerYear * simStart.yearFrac(simEnd));
            if(maxNumDivsPerYear != -1 && maxNumDivsPerYear != 0) {
                // Use a minimum number of divs
                totalNbDivs = Maths::max(totalNbDivs, BarrierData::MIN_DIVS_PER_PERIOD);
            }

            // Brownian Bridge data
            bridgeData = vector<BridgeDataSP>(maxNumBBs);
            vector<TimeMetricConstSP> assetTimeMetrics(numAssets);
            vector<DividendListSP> driverDivs(numAssets);
            for(iAsset = 0; iAsset < numAssets; iAsset++) {
                doHitNoHit[iAsset]   = BoolArray(maxNumBBs);
                hitNoHitPath[iAsset] = IntArray(maxNumBBs);

                // Get asset time metric
                ATMVolRequestSP  volRequest(new ATMVolRequest());
                CVolProcessedSP vol(mAsset->factorGetProcessedVol(iAsset, volRequest.get()));
                assetTimeMetrics[iAsset] = vol->GetTimeMetric();

                // Get asset dividends and convert them to driver yields by multiplying by (-1)
                // swince driver is decreasing function of spot
                DividendListSP divs = AssetUtil::getDiscreteDivs(
                    &mAsset->getAsset(iAsset), valueDate, simStart, simEnd,
                    totalNbDivs, DividendCollector::DOLLAR_TO_YIELD);
                DividendArray& divArray = const_cast<DividendArray&>(divs->getArray());
                for(int iDiv = 0; iDiv < divArray.size(); iDiv++) {
                    divArray[iDiv].scale(-1.0);
                }
                driverDivs[iAsset] = divs;
            }

            // Flag if date is useful or not
            for(iStep = 0; iStep < maxNumBBs; iStep++) {
                int iNextStep = iStep + 1;
                const DateTime& lastDate = simTimeline[iStep];
                const DateTime& thisDate = simTimeline[iStep + 1];
                BridgeDataSP bbData(new BridgeData(lastDate, thisDate, numAssets));
                for(iAsset = 0; iAsset < numAssets; iAsset++) {
                    // Get asset time metric
                    ATMVolRequestSP  volRequest(new ATMVolRequest());
                    CVolProcessedSP vol(mAsset->factorGetProcessedVol(iAsset, volRequest.get()));
                    TimeMetricConstSP timeMetric(vol->GetTimeMetric());

                    for(iBarrier = 0; iBarrier < numBarriers; iBarrier++) {
                        BarrierDataConstSP data = barrierHVBBGens[iBarrier]->getBarrierData();
                        BarrierPerAssetDataConstSP assetBarrier = data->getAssetBarrier(iAsset);
                        DateTimeArraySP monitoringDates = assetBarrier->monitoringDates;

                        // Define the start and end date for this asset's barrier
                        DateTime startDate = monitoringDates->front();
                        const DateTime& endDate   = monitoringDates->back();
                        DateTimeArray firstFutureDailyMonDate = assetBarrier->getFirstFutureDailyMonDate(
                            mAsset->getAsset(iAsset) , valueDate, data->monitorType);
                        if(firstFutureDailyMonDate.size()) {
                            // This will make sure we will not do BBs before firstFutureDailyMonDate
                            startDate = firstFutureDailyMonDate[0];
                        }

                        if(startDate < thisDate && thisDate <= endDate) {
                            doHitNoHit[iAsset][iStep] = true; // sample the asset
                            // Get the reference level from the appropriate barrier object
                            IRefLevel::IStateVarGenSP refLevelGen = barrierHVBBGens[iBarrier]->
                                getRefLevelGen();
                            IRefLevel::IStateVarSP oldStateVar =
                                refLevelGen->getRefLevelSV(IRefLevel::IStateVarSP(   ), pastPathGen);
                            IRefLevel::IStateVarSP refLevelSV =
                                refLevelGen->getRefLevelSV(oldStateVar, futPathGen);
                            bool isRefLevelInPast = refLevelGen->getDates(iAsset).back() <= valueDate;

                            // Brownian Bridge object
                            // Note that driver is decreasing function of spot since the implied
                            // distribution is computed using call spreads as opposed to put spreads.
                            // So an Up barrier for the spot is equivalent to a down barrier for
                            // the driver
                            double driverVol = driver->driverVols[iAsset][iNextStep];
                            double driverVar = Maths::square(driverVol);
                            bool driverIsUp = !assetBarrier->isUp;
                            double bridgeLength = assetTimeMetrics[iAsset]->yearFrac(lastDate, thisDate);

                            // See if we need to do a barrier shift for daily monitoring
                            bool adjust = CString::equalsIgnoreCase(
                                data->monitorType, BarrierData::DAILY_MONITORING);
                            double barrierShift = 0.0;
                            if(adjust) {
                                // Get business days between dates using holidays
                                HolidayConstSP hols = AssetUtil::getHoliday(&mAsset->getAsset(iAsset));
                                int numBusDays = hols->businessDaysDiff(lastDate, thisDate);
                                if(numBusDays) {
                                    int isUp = driverIsUp ? 1 : -1;
                                    barrierShift = isUp * ADJUST_CONSTANT * driverVol / sqrt((double)(numBusDays));
                                }
                            }

                            // The BB data
                            ScheduleSP schedule = assetBarrier->levels;
                            double barrierStart = schedule->interpolate(lastDate);
                            double barrierEnd;
                            if(CString::equalsIgnoreCase(schedule->getInterp(), Schedule::INTERP_STAIRS)) {
                                // Step interpolation: use start barrier
                                barrierEnd = barrierStart;
                            } else {
                                // Other i.e. None of Linear: use end barrier
                                barrierEnd =  schedule->interpolate(thisDate);
                            }

                            // Obtain dividends for this brownian bridge
                            DividendListSP assetBridgeDivs(driverDivs[iAsset]->getAllDivsBetweenDates(lastDate, thisDate));
                            DoubleArraySP divTimes(new DoubleArray(assetBridgeDivs->getArray().size()));
                            DateTimeArrayConstSP exDivDates = assetBridgeDivs->getExDivDates();
                            assetTimeMetrics[iAsset]->yearFrac(lastDate, *exDivDates, *divTimes);

                            HitNoHitBBSP bb(new HitNoHitBB(driverVar, bridgeLength, driverIsUp,
                                            assetBridgeDivs->getArray(), *divTimes));

                            bbData->assetData[iAsset] = BridgeData::AssetBridgeDataSP(new
                                BridgeData::AssetBridgeData(iAsset,
                                                            iStep,
                                                            driver,
                                                            barrierStart,
                                                            barrierEnd,
                                                            refLevelSV,
                                                            isRefLevelInPast,
                                                            bb,
                                                            bridgeLength,
                                                            adjust,
                                                            barrierShift,
                                                            assetBridgeDivs,
                                                            divTimes));

                            // Append to offsets
                            offsets[iBarrier][iAsset].push_back(iStep);
                        }
                    }
                }
                bridgeData[iStep] = bbData;
            }

            // Count on how many dates we need to do something
            numBBs = 0;
            for(iStep = 0; iStep < maxNumBBs; iStep++) {
                for(iAsset = 0; iAsset < numAssets; iAsset++) {
                    if(doHitNoHit[iAsset][iStep]) {
                        numBBs++;
                        break;
                    }
                }
            }

            // Populate the map of stateVarGens, StateVars
            for(iBarrier = 0; iBarrier < numBarriers; iBarrier++) {
                const SVGenBarrierHVBB* barrierBBGen = barrierHVBBGens[iBarrier];
                svDBase.append(barrierBBGen, SVGenBarrierHVBB::StateVarSP(new
                    SVGenBarrierHVBB::StateVar(pastPathGen,
                                               futPathGen,
                                               barrierBBGen->getBarrierData(),
                                               barrierBBGen->getRefLevelGen(),
                                               barrierBBGen->getSpotSmoothGen(),
                                               barrierBBGen->getSpotAtMonStartGen(),
                                               barrierBBGen->getSpotAtMonEndGen(),
                                               offsets[iBarrier],
                                               hitNoHitPath)));
            }
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Contains N-factor data for BB between 2 dates */
    class BridgeData {
    public:
        /** Contains 1-factor data for BB between 2 dates */
        class AssetBridgeData {
        public:
            /** Full constructor */
            AssetBridgeData(int iAsset,
                            int iStep,
                            PathGenSpot::DriverConstSP driver,
                            double barrierStartPct,
                            double barrierEndPct,
                            IRefLevel::IStateVarSP refLevelSV,
                            bool isRefLevelInPast,
                            HitNoHitBBSP bbSample,
                            double bridgeLength,
                            bool adjust,
                            double barrierShift,
                            DividendListSP divList,
                            DoubleArraySP divTimes):
            iAsset(iAsset), iStep(iStep), driver(driver), barrierStartPct(barrierStartPct),
            barrierEndPct(barrierEndPct), refLevelSV(refLevelSV), isRefLevelInPast(isRefLevelInPast),
            bbSample(bbSample), bridgeLength(bridgeLength), adjust(adjust), barrierShift(barrierShift),
            divList(divList), divTimes(divTimes), barrierStart(0.0), barrierEnd(0.0) {
                // Compute weight for forward problem
                const DividendArray& divArray = divList->getArray();
                double meanTime = 0.0;
                sumDivs  = 0.0;
                for(int iDiv = 0; iDiv < divArray.size(); iDiv++) {
                    double tDiv = (*divTimes)[iDiv];
                    sumDivs  += divArray[iDiv].getDivAmount();
                    meanTime += divArray[iDiv].getDivAmount() * tDiv;
                }
                if(!Maths::isZero(sumDivs) && bridgeLength > 0.0) {
                    meanTime /= sumDivs;
                    fwdWeight = 1.0 - meanTime / bridgeLength;
                } else {
                    // Just use the forward problem
                    fwdWeight = 1.0;
                }

                if(isRefLevelInPast) {
                    computeBarrierLevels(barrierStart, barrierEnd);
                }
            }

            /** Computes or returns precomputed barrier */
            void getBarrierLevels(double& start, double& end) const {
                if(isRefLevelInPast) {
                    start = barrierStart;
                    end   = barrierEnd;
                } else {
                    computeBarrierLevels(start, end);
                }
            }

            /** Destructor */
            virtual ~AssetBridgeData() {}

            // Mandatory fields
            int                         iAsset;               //!< Asset index
            int                         iStep;                //!< Indicates iStep on the driver's timeline
            PathGenSpot::DriverConstSP  driver;               //!< Driver
            double                      barrierStartPct;      //!< Barrier % at start of BB
            double                      barrierEndPct;        //!< Barrier % at end of BB
            IRefLevel::IStateVarSP      refLevelSV;           //!< Reference level state variable
            bool                        isRefLevelInPast;     //!< Whether refLevel is in the past i.e. has been determined
            HitNoHitBBSP                bbSample;             //!< Brownian Bridge
            double                      bridgeLength;         //!< Length of bridge
            bool                        adjust;               //!< Whether requires adjustment for "D" monitoring
            double                      barrierShift;         //!< Barrier shift for "D" monitoring
            DividendListSP              divList;              //!< Asset dividends for BB
            DoubleArraySP               divTimes;             //!< Dividend time relative to start of bridge

            // Transient fields
            double                      sumDivs;              //!< Sum of dividends for this BB
            double                      fwdWeight;            //!< Weight for forward problem
            double                      barrierStart;         //!< Precomputes driver's barrier if possible
            double                      barrierEnd;           //!< Precomputes driver's barrier if possible

        private:
            /** Computes driver's barrier */
            void computeBarrierLevels(double& start, double& end) const {
                const vector<int>& relativeAsset        = driver->isRelative[iAsset];
                const vector<const double*>& refReturns = driver->refReturns[iAsset];

                double ref      = refLevelSV->refLevel(iAsset);
                double refStart = (relativeAsset[iStep] == 0) ? 1.0 : (*refReturns[iStep]);
                double refEnd   = (relativeAsset[iStep + 1] == 0) ? 1.0 : (*refReturns[iStep + 1]);

                // The spots that enter the mapping are
                double spotBarrierStart = ref * barrierStartPct / refStart;
                double spotBarrierEnd   = ref * barrierEndPct / refEnd;

                // Map barriers to driver barriers (they vary with the path due to refLevel)
                double barrierStart = driver->driverFromSpot(spotBarrierStart, iAsset, iStep);
                double barrierEnd   = driver->driverFromSpot(spotBarrierEnd, iAsset, iStep + 1);

                // Further adjust the data for daily monitoring
                if(adjust) {
                    barrierStart += barrierShift;
                    barrierEnd   += barrierShift;
                }

                // Return barriers
                start = barrierStart;
                end   = barrierEnd;
            }

            /** Disabled default constructor */
            AssetBridgeData();

        };

        typedef refCountPtr<AssetBridgeData> AssetBridgeDataSP;

        /** Full constructor */
        BridgeData(const DateTime& startDate,
                   const DateTime& endDate,
                   int numAssets):
        startDate(startDate), endDate(endDate), assetData(numAssets) {}

        /** Destructor */
        virtual ~BridgeData() {}

        DateTime                  startDate;          //!< Start date for the BB
        DateTime                  endDate;            //!< End date for the BB
        vector<AssetBridgeDataSP> assetData;          //!< Brownian Bridge Per asset data

    private:
        /** Disabled default constructor */
        BridgeData();
    };

    typedef refCountPtr<BridgeData> BridgeDataSP;


    // Fields
    int                     numBarriers;        //!< Number of barriers
    int                     numAssets;          //!< Number of assets
    const IMultiFactors*    mAsset;             //!< Multi asset interface
    int                     maxNumBBs;          //!< Maximum nb of BBs = driver dates - 1
    int                     numBBs;             //!< Number of steps for BB
    DependenceSP            dependence;         //!< Dependence object e.g. gaussian
    MCRandomGenCacheSP      randomGen;          //!< Random number generator for normals

    MCHitValuePathCacheSP   hitValuePathCache;  //!< Hit value path cache
    MCRandomCacheSP         randCache;          //!< Random number cache

    vector<IntArray>        hitNoHitPath;       //!< HitNoHit flag for BB [iAsset][iStep]
    vector<BoolArray>       doHitNoHit;         //!< Whether HitNoHit is required [iAsset][iStep]

    PathGenSpotSP           pathGenSpot;        //!< Spot path generator
    PathGenSpot::DriverConstSP driver;          //!< Driver of pathGenSpot

    vector<BridgeDataSP>    bridgeData;         //!< Data for Brownian Bridges [iStep]

    // tabulated normal distribution. Use fast non virtual inteprolant
    EquidistantLinearInterpolantNonVirtualSP tabulatedNorm; //!< Tabulated normal distribution

    static const MCCache::KeyUtils::Key keyHitValueCache;   //!< Key for hit no hit path cache
    static const MCCache::KeyUtils::Key keyRandCache;       //!< Key for hit no hit path random numbers cache
    static const double ADJUST_CONSTANT;                    //!< Barrier adjustment constant
};

// Initialize the keys to the caches
const MCCache::KeyUtils::Key MCPathConfigImplied::PathGenHVBB::keyHitValueCache = MCCache::KeyUtils::getNewKey();
const MCCache::KeyUtils::Key MCPathConfigImplied::PathGenHVBB::keyRandCache = MCCache::KeyUtils::getNewKey();
const double MCPathConfigImplied::PathGenHVBB::ADJUST_CONSTANT = 0.5826;



/////////////////////////////////////////////////////////////////////////////


MCPathConfigImplied::Gen::Gen(MCPathConfigImplied*     pathConfig,
                              const MCPathGeneratorSP& pastPathGenerator,
                              const MCProductClient*   prodClient,
                              MCCacheManager&          cacheMgr,
                              bool                     cachingRequested,
                              const MCImpliedCacheSP&  impliedCache,
                              Results*                 results,
                              DateTimeArray&           simDates):
nowPathIdx(0) {
    static const string routine = "MCPathConfigImplied::Gen::Gen";

    try {
        // Collect state variables from product and categorize them
        StateVariableCollectorSP svCollector(new StateVariableCollector());
        prodClient->collectStateVars(svCollector);
        IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

        // 1) Spot requests
        SVGenSpotArray spotGenArray = filterStateVars<SVGenSpot>(stateVarGenArray);
        if(!spotGenArray.size()) {
            throw ModelException("No spot paths specified.");
        }

        // 2a) Brownian Bridge Barrier hit value requests
        SVGenBarrierHVBBArray barrierHVBBGens =
            filterStateVars<SVGenBarrierHVBB>(stateVarGenArray);


        // 3) Create discount factor state variables
        SVGenDiscFactorArray discFactors =
            filterStateVars<SVGenDiscFactor>(stateVarGenArray);
        unsigned int iVar;
        for(iVar = 0; iVar < discFactors.size(); iVar++) {
            IStateVariableSP sv(discFactors[iVar]->
                             determinsticSV(false /* not doing past */));
            svDBase.append(discFactors[iVar], sv);
        }

        // 4) Create expected discount factor (ZCBs) state variables
        vector<const SVGenExpectedDiscFactor*>  expDiscFactors(
            filterStateVars<SVGenExpectedDiscFactor>(stateVarGenArray));
        const DateTime& today = prodClient->getToday();
        for(iVar = 0; iVar < expDiscFactors.size(); iVar++) {
            IStateVariableSP sv(expDiscFactors[iVar]->determinsticSV(
                                 today, false /* not doing past */));
            svDBase.append(expDiscFactors[iVar], sv);
        }

        // Nothing else is supported
        if(!stateVarGenArray.empty()) {
            throw ModelException("Unable to recognize all state variable types.");
        }

        // Create driver generators for BarrierHVs
        DateTimeArraySP allDriverDates =
            SVGenBarrierHVBB::BarrierDatesHelper::getBarrierTimeline(
            prodClient, barrierHVBBGens, simDates);

        // Create a spot path generator for the SVGenSpot and MCDriver
        PastPathGen* pastPathGen = dynamic_cast<PastPathGen*>(pastPathGenerator.get());
        if(!pastPathGen) {
            throw ModelException("Past path generator is not of PastPathGen type.");
        }

        // Do any modules require using the driver of the pathGenSpot ?
        bool requestDriver = barrierHVBBGens.size() ? true: false;

        pathGenSpot = MCPathConfigImplied::PathGenSpotSP(new
            MCPathConfigImplied::PathGenSpot(pathConfig, pastPathGen->getPathGenSpot(),
            prodClient, spotGenArray, *allDriverDates, cacheMgr, cachingRequested,
            svDBase, impliedCache, results, requestDriver, simDates));

        // Create a hit value path generator. Pass the driverGens and the
        // pathGenSpot so that the HV simulator can obtain the driver SVs
        if(barrierHVBBGens.size()) {
            // throw ModelException("Cannot do BBs yet");
            pathGenHVBB = MCPathConfigImplied::PathGenHVBBSP(new
            MCPathConfigImplied::PathGenHVBB(pathConfig, prodClient,
            pastPathGen, this, pathGenSpot, barrierHVBBGens,
            cacheMgr, cachingRequested, svDBase, pathConfig->tabulatedNorm));
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


int MCPathConfigImplied::Gen::NbSimAssets() const {
    static const string routine = "MCPathConfigImplied::Gen::NbSimAssets";
    throw ModelException(routine, "Method is retired for StateVars");
}

const double* MCPathConfigImplied::Gen::Path(int iAsset, int iPath) const {
    static const string routine = "MCPathConfigImplied::Gen::Path";
    throw ModelException(routine, "Method is retired for StateVars");
};

double MCPathConfigImplied::Gen::refLevel(int iAsset, int iPath) const {
    static const string routine = "MCPathConfigImplied::Gen::refLevel";
    throw ModelException(routine, "Method is retired for StateVars");
}

double MCPathConfigImplied::Gen::maxDriftProduct(int iAsset) const {
    static const string routine = "MCPathConfigImplied::Gen::maxDriftProduct";
    try {
        return pathGenSpot->maxDriftProduct(iAsset);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

int MCPathConfigImplied::Gen::begin(int iAsset) const {
    static const string routine = "MCPathConfigImplied::Gen::begin";
    throw ModelException(routine, "Method is retired for StateVars");
}

int MCPathConfigImplied::Gen::end(int iAsset) const{
    static const string routine = "MCPathConfigImplied::Gen::end";
    throw ModelException(routine, "Method is retired for StateVars");
}

// Live methods
bool MCPathConfigImplied::Gen::hasPast() const {
    return pathGenSpot->hasPast();
}

bool MCPathConfigImplied::Gen::doingPast() const {
    return false;
}

void MCPathConfigImplied::Gen::generatePath(int pathIdx) {
    nowPathIdx = pathIdx;
    pathGenSpot->generatePath(pathIdx);
    if(pathGenHVBB.get()) {
        pathGenHVBB->generatePath(pathIdx);
    }
}


void MCPathConfigImplied::Gen::advance() {
    pathGenSpot->advance();
}

void MCPathConfigImplied::Gen::reset() {
    pathGenSpot->reset();
}

int MCPathConfigImplied::Gen::getPathIndex() const {
    return nowPathIdx;
}

/** Returns the state variable corresponding to generator.
    Part of the IStateVariableGen::IStateGen IFace */
IStateVariableSP MCPathConfigImplied::Gen::create(const IStateVariableGen* svGen) {
    static const string routine = "MCPathConfigImplied::Gen::create";

    try {
        return svDBase.find(svGen);
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

MCImpliedCacheSP MCPathConfigImplied::Gen::getCache() const {
    return pathGenSpot->getCache();
}


/////////////////////////////////////////////////////////////////////////////


/** Creates a past path generator */
MCPathGeneratorSP MCPathConfigImplied::pastPathGenerator(const IMCProduct* prod) {
    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
    if(!prodClient){
        // Old approach
        return MCPathConfig::pastPathGenerator(prod);
    } else {
        // State variables approach
        return MCPathGeneratorSP(new PastPathGen(prodClient));
    }
}


/** Creates a PathGenerator used for future dates which supports
    'random access' to the paths. That is if a path generator is
    created using savePaths = true, then subsequently another path
    generator can be created (with savePaths = false) which will
    support the pathIdx parameter in the generatePath() method, ie the
    paths can be generated in any order. */
MCPathGeneratorSP MCPathConfigImplied::futurePathGenerator(
    int                      cachingMode,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control,
    Results*                 results,
    DateTimeArray&           simDates ) {

    static const string routine = "MCPathConfigImplied::futurePathGenerator";

    try {
        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
        if(!prodClient) {
            // Old approach
            return MCPathConfig::futurePathGenerator(cachingMode, numPaths,
                pastPathGenerator, prod, control, results, simDates);
        } else {
            pathCache = PathCacheSP(new PathCache());
            randomCache = RandomCacheSP(new RandomCache());

            // State variables approach
            bool isPricing = control->isPricing();
            bool cachePaths = (cachingMode & PATH_CACHE)? true: false;
            // cachePaths => cacheRandoms, so no need to save states
            bool mySaveStates = !cachePaths && (cachingMode & PATH_RANDOM_ACCESS);
            /* no point saving paths if no greeks */
            if (mySaveStates){
                throw ModelException("MCPathConfigImplied::futurePathGenerator",
                                     "QuickGreeks are not supported without path caching");
            }

            if (cachePaths || isPricing) {
                // Kill them
                cacheMgr = MCCacheManager();
            } else {

                const IMultiFactors* multiFactor = prod->getMultiFactors();
                int numAssets = prod->getNumAssets();
                SensitivitySP sens(control->getCurrentSensitivity());
                // getSensitiveAssets() does the real work
                IntArray sensitivePhiAssets(multiFactor->getSensitiveAssets(
                                            sens.get(), true)); // include phi etc
                if (sensitivePhiAssets.empty()){
                    IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                                             sens.get(), false)); // exclude phi etc
                    cacheMgr.disallowAssets(sensitiveAssets);
                } else {
                    // mark all paths as invalid
                    IntArray sensitiveAssets(numAssets);
                    for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                        sensitiveAssets[iAsset] = iAsset;
                    }
                    cacheMgr.disallowAssets(sensitiveAssets);
                }
            }

            // Bit comparison to deduce whether user has asked for caching
            bool cachingRequested = (cachingMode & IMCPathConfig::PATH_CACHE) ? true : false;

            MCPathGeneratorSP gen = makePathGenerator(
                cachingRequested, numPaths, pastPathGenerator, prod, control, results, simDates);

            if (cachePaths || isPricing) {
                // Caches have been built. Configure them now
                cacheMgr.configureCache(numPaths);
            }

            return gen;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


MCPathGeneratorSP MCPathConfigImplied::makePathGenerator(
    bool                     cachingRequested,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control,
    Results*                 results,
    DateTimeArray&           simDates ) {

    static const string routine = "MCPathConfigImplied::makePathGenerator";

    try {
        // get hold of the random number generator in parent
        IRandomSP randGen(getRandomGenerator());
        /* switch between general case and specific case where each asset
           has the same simulation dates */
        MCImpliedCacheSP copyOfCache;
        if(impliedCache.get()) {
            if (control->isPricing()){
                // important to reuse cache here for case when overall calc is
                // split into 'sub-blocks'
                copyOfCache = impliedCache;
            } else {
                // Work out which assets which are being tweaked
                // Start by gettting hold of MultiFactors
                const IMultiFactors* multiFactor = prod->getMultiFactors();
                // then which greek we're doing
                SensitivitySP sens(control->getCurrentSensitivity());
                // getSensitiveAssets() does the real work
                IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                                             sens.get(), false)); // exclude phi etc
                /* copy our cache - this is poor performance (need a
                   selective 'shallow' copy) */
                copyOfCache = MCImpliedCacheSP(impliedCache.clone());
                copyOfCache->killProbs(sensitiveAssets);
            }
        } else {
            copyOfCache = impliedCache; // NULL cache
        }

        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
        if(!prodClient){
            // Old approach

            // Create the MCPathBaseImplied class
            refCountPtr<MCPathBaseImplied> futurePathGenBaseImplied(
                new MCPathBaseImplied(this,
                                      pastPathGenerator,
                                      prod,
                                      copyOfCache,
                                      results));
            if (!impliedCache){
                // for new cache, save for future use
                impliedCache = futurePathGenBaseImplied->getCache();
            }

            // Construct the MCPathGenerator
            return MCPathBase::createPathGenerator(pastPathGenerator,
                                                   this,
                                                   prod,
                                                   futurePathGenBaseImplied);

        } else {
            // State variables approach
            GenSP pathGen(new MCPathConfigImplied::Gen(
                this,
                pastPathGenerator,
                prodClient,
                cacheMgr,
                cachingRequested,
                copyOfCache,
                results,
                simDates));

            if (!impliedCache){
                // for new cache, save for future use
                impliedCache = pathGen->getCache();
            }

            return pathGen;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


int MCPathConfigImplied::storagePerPath(IMCProduct* product) const{
    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(product);
    if(!prodClient){
        // Call parent
        return MCPathConfig::storagePerPath(product);
    } else {
        // State variables approach

        if (!product->hasFuture()) {
            return 0; // nothing to simulate
        }

        // Let's just create one on the fly...
        smartPtr<MCPathConfigImplied> pathConfig(copy(this)); // copy to avoid const problems
        MCPathGeneratorSP pastPathGenerator(pathConfig->pastPathGenerator(product));
        product->pathGenUpdated(pastPathGenerator.get());
        // Calls to getVolInterp expect payoff() to have been called for the past
        if (pastPathGenerator->hasPast()) {
            // create our IMCPrices object
            IMCPricesSP prices(product->createOrigPrices(1, 1, 0));
            pastPathGenerator->generatePath(0); // may well do nothing
            product->payoff(pastPathGenerator.get(), *prices);
        }
        SensitivityArrayConstSP    sens;
        OutputRequestArrayConstSP  request;
        Control control(sens, request, false, "");
        Results results;
        DateTimeArray simDates;
        MCPathGeneratorSP pathGen(pathConfig->futurePathGenerator(2, // cache
                                                                  2, // num paths
                                                                  pastPathGenerator,
                                                                  product,
                                                                  &control,
                                                                  &results,
                                                                  simDates ));

        int memRequirement = pathConfig->cacheMgr.storagePerRun();
        return memRequirement;
    }
}

#endif


/** "bread & butter" path config that can be captured in Pyramid using
    current IMS */
class MCPathConfigImpliedDefault: public MCPathConfigImplied {
public:
    static CClassConstSP const TYPE;

private:
    MCPathConfigImpliedDefault():MCPathConfigImplied(
        TYPE,
        IVolatilityBS::TYPE->getName(),
        "not used"){}

    void validatePop2Object() {
        static const string method = "MCPathConfigImpliedDefault::validatePop2Object";

        try {
            // Call parent method
            MCPathConfigImplied::validatePop2Object();
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigImpliedDefault, clazz);
        SUPERCLASS(MCPathConfigImplied);
        EMPTY_SHELL_METHOD(defaultMCPathConfigImpliedDefault);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigImpliedDefault(){
        return new MCPathConfigImpliedDefault();
    }
};

CClassConstSP const MCPathConfigImpliedDefault::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigImpliedDefault", typeid(MCPathConfigImpliedDefault), load);


// external symbol to allow class to be forced to be linked in
bool MCPathConfigImpliedLoad(){
    return (MCPathConfigImplied::TYPE != 0 && MCPathConfigImpliedDefault::TYPE != 0);
}



////////////////////////////////////////
// MCPathBaseImplied
////////////////////////////////////////
MCPathBaseImplied::MCPathBaseImplied(
    MCPathConfigImplied*     pathConfig,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    const MCImpliedCacheSP&  suppliedCache,
    Results*                 results):
    impliedCache(suppliedCache), // might be null
    results(results),
    mAsset(prod->getMultiFactors()),
    numAssets(prod->getNumAssets()),
    productPaths(prod->getNumAssets()),
    tabulatedNorm(pathConfig->tabulatedNorm),
    // need to pass on the following three fields to LocalCorr
    dependenceMaker(pathConfig->getDependenceMaker()),
    rand(pathConfig->getRandomGenerator()), 
    idiosynFactorRandomCache(pathConfig->getIdiosynFactorRandomCache()),
    marketFactorRandomCache(pathConfig->getMarketFactorRandomCache()), 
    isCarefulRandoms(pathConfig->carefulRandoms()) {

    static const string routine("MCPathBaseImplied::MCPathBaseImplied");
    try{
        // Obtain product timeline
        timeline = getProductTimeline(prod, pastPathGenerator);

        // Obtain market data
        IntArray nbPaths(numAssets, 1);
        refData = getRefData(timeline, nbPaths, mAsset, prod, pastPathGenerator);

        // Paths
        int iAsset;
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            productPaths[iAsset] = DoubleArray(timeline->totalNumSteps);
        }

        // Start initializing member variables
        driverPath  = DoubleMatrix(numAssets, timeline->numFutSteps);
        cumulatives = DoubleMatrix(numAssets, timeline->numFutSteps);

        // Deduce the cliquet timeline and obtain volRequests for each cliquet
        vector< vector<VolRequestLNStrikeSP> > volRequests; //[iCliquet][iAsset]
        getVolRequests(prod,
                       pastPathGenerator,
                       cliquetTimeline,
                       volRequests,
                       numCliquets);

        const IMCProductImplied* prodImplied =
            dynamic_cast<const IMCProductImplied*>(prod);
        if (prodImplied)
        {
            prodImplied->initialiseImplied(pastPathGenerator.get());
        }

        // Categorize the simulation dates on a per cliquet basis
        int iCliquet, iStep;
        vector<DateTimeArray> datesTo(numCliquets);
        for(iStep = 0; iStep < timeline->numFutSteps; iStep++) {
            const DateTime& date = timeline->futureDates[iStep+1];
            // figure out on which cliquet the date belongs
            for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                DateTime& startDate = cliquetTimeline[iCliquet];
                DateTime& endDate   = cliquetTimeline[iCliquet+1];
                if(startDate <= date && date <= endDate) {
                    // record the date
                    datesTo[iCliquet].push_back(date);
                    break;
                }
            }
        }

        // Record indices of the timeline at the beginning of each cliquet
        // Required for simulating paths by looping over cliquets
        cliquetDiffOffsets = IntArray(numCliquets);
        for(iCliquet = 1; iCliquet < numCliquets; iCliquet++) {
            cliquetDiffOffsets[iCliquet] = datesTo[iCliquet-1].size();
        }

        // initialize the cache, unless it exists from some previous tweak
        if(!impliedCache) {
            impliedCache = MCImpliedCacheSP(
                new MCImpliedCache(numCliquets));
        }

        // Kill the cache if we have rolled over a cliquet
        if(numCliquets != impliedCache->numCliquets) {
            impliedCache->killCliquetCaches(numCliquets);
        }

        // construct array of cliquets
        cliquetArray = vector<CliquetSP>(numCliquets);
        int offset = 0;
        const DateTime& today = prod->getToday();
        for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            try {
                DateTime cliqStartDate = cliquetTimeline[iCliquet];
                CliquetCacheSP& cliquetCache = (*impliedCache->cliquetCaches)[iCliquet];
                offset += cliquetDiffOffsets[iCliquet];
                CliquetSP cliquet;
                if(iCliquet == 0) {
                    // For 1st cliquet invoke constructor that will use
                    // deterministic reference
                    cliquet = CliquetSP(
                        new Cliquet(iCliquet,
                                    volRequests[iCliquet],
                                    mAsset,
                                    today,
                                    cliqStartDate,
                                    datesTo[iCliquet],
                                    pathConfig->pdfParams,
                                    pathConfig->impliedParams,
                                    cliquetCache,
                                    pathConfig->cacheStrikes,
                                    tabulatedNorm));
                } else {
                    // Create stochastic reference levels for subsequent cliquets
                    vector<const double*> stochasticReference;
                    for(iAsset = 0; iAsset < numAssets; iAsset++) {
                        const double* ref =
                            &productPaths[iAsset][timeline->numPastDates + offset - 1];
                        stochasticReference.push_back(ref);
                    }
                    cliquet = CliquetSP(
                        new Cliquet(iCliquet,
                                    volRequests[iCliquet],
                                    mAsset,
                                    today,
                                    cliqStartDate,
                                    datesTo[iCliquet],
                                    pathConfig->pdfParams,
                                    pathConfig->impliedParams,
                                    cliquetCache,
                                    pathConfig->cacheStrikes,
                                    tabulatedNorm,
                                    stochasticReference));
                }
                cliquetArray[iCliquet] = cliquet;
            } catch(exception& e) {
                string message = "Failed to construct cliquet number " +
                    Format::toString(iCliquet + 1);
                throw ModelException::addTextToException(e, message);
            }
        }

        // Write information to the DEBUG packet
        if(!impliedCache->isResultsExported) {
            // Record an array of infos for each asset (one per cliquet)
            CliquetCacheArraySP& cliquetCaches = impliedCache->cliquetCaches;
            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                string assetName;
                try {
                    assetName = mAsset->assetGetName(iAsset);
                    DebugInfoArraySP xmlOutput(new DebugInfoArray(numCliquets));

                    for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                        CliquetCacheSP& cliquetCache = (*cliquetCaches)[iCliquet];
                        AssetCacheSP&   assetCache   = (*cliquetCache->assetCaches)[iAsset];

                        // The Fwd adjustment
                        int numPDFSteps = assetCache->numPDFSteps;
                        DoubleArraySP adjustment(new DoubleArray(numPDFSteps));
                        for(iStep = 0; iStep < numPDFSteps; iStep++) {
                            (*adjustment)[iStep] = (*assetCache->pdfSamples)[iStep]->getAdjustment();
                        }

                        (*xmlOutput)[iCliquet] = DebugInfoSP(new DebugInfo(assetCache->partition, adjustment));
                    }

                    results->storeGreek(xmlOutput,
                                        Results::DEBUG_PACKET,
                                        OutputNameSP(new OutputName(assetName)));
                } catch(exception& e) {
                    string message = "Failed to write partitions and adjustments to Results for asset " + assetName;
                    throw ModelException::addTextToException(e, message);
                }
            }

            // switch off writing to the DEBUG
            impliedCache->isResultsExported = true;
        }
        
        if (pathConfig->skewMaker.get() && (numAssets>1)) {
            // set up dependence
            DependenceMakerLocalCorr* dpm = dynamic_cast<DependenceMakerLocalCorr*>(pathConfig->skewMaker.get());
            if (!dpm) {
                throw ModelException(routine, "Internal Error");
            }
            DependenceSP dependence(dpm->createDependence(this));
            // a special one IRandomSP
            randomGen = IMCRandomSP(new MCRandom(
                this,                               // path generator
                dependence,                         // dependence
                IRandomSP(),                        // random number generator
                pathConfig->getRandomCache(),       // randomCache
                true,                               // not used (isCarefulRandoms)
                timeline->numFutSteps,              // numDates
                numAssets,                          // numFactors
                timeline->totalNumPastDates));      // numPastDates
        } else {
            // set up dependence
            DependenceSP dependence(pathConfig->dependenceMaker->createDependence(this));
            // the classic one         
            randomGen = IMCRandomSP(new MCRandom(
                this,                               // path generator
                dependence,                         // dependence
                pathConfig->getRandomGenerator(),   // random number generator
                pathConfig->getRandomCache(),       // randomCache
                pathConfig->carefulRandoms(),       // isCarefulRandoms
                timeline->numFutSteps,              // numDates
                numAssets,                          // numFactors
                timeline->totalNumPastDates));      // numPastDates
        }

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void MCPathBaseImplied::getVolRequests(
    const IMCProduct*                        prod,
    const MCPathGeneratorSP&                pastPathGenerator,
    DateTimeArray&                          cliquetTimeline,
    vector< vector<VolRequestLNStrikeSP> >& volRequests,
    int&                                    numCliquets) {
    static const string routine("MCPathBaseImplied::GetVolRequests");
    try{
        // get vol interp from product at this stage. Later on we would like
        // to get the VolRequest from somewhere else. The user should, ideally,
        // pass a VolRequest to MCImplied.
        const IMCProductLN* prodLN = dynamic_cast<const IMCProductLN*>(prod);
        if (!prodLN){
            throw ModelException("Product does not implement "
                                 "IMCProductLN interface");
        }

        // Get from the IMCProductLN the first VolRequest for the first asset
        int iAsset = 0;
        CVolRequestLNArray prodVolRequests(prodLN->getVolInterp(pastPathGenerator.get(), iAsset));
        if (prodVolRequests.empty()){
            throw ModelException(routine, "No vol requests specified");
        }
        CVolRequestLNSP volReq(prodVolRequests[0]);

        // Deduce if it a cliquetVolRequest or VolRequestLNStrike
        VolRequestLNStrike* volReqLNStr = dynamic_cast<VolRequestLNStrike*>(volReq.get());
        CliquetVolRequest*  cliquetVol  = dynamic_cast<CliquetVolRequest*>(volReq.get());

        bool isCliquet;
        if(volReqLNStr) {
            isCliquet = false;
        } else if(cliquetVol) {
            isCliquet = true;
        } else {
            throw ModelException(routine, "Vol interp is neither LNStrike "
                                 "type nor Cliquet type for asset " +
                                 mAsset->assetGetName(iAsset));
        }

        // The cliquet TimeLine is the intersection of:
        // i) StartDay (if not cliquet) or cliquet dates (if cliquet)
        // ii) last simulation date
        vector<DateTimeArray> dates(2);
        DateTime startDate = prod->getEffectiveSimStartDate();
        int iCliquet;
        if(isCliquet) {
            dates[0] = cliquetVol->getCliqStartDates();                      // additional cliquet dates
            // the first live cliquet date might be in the past (started product)
            // or in the future (1st cliquet AvgIn). To make consistent with
            // non-cliquet we set it to startDate
            dates[0][0] = startDate;
        } else {
            dates[0] = DateTimeArray(1, startDate);                          // product's start day
        }
        dates[1] = DateTimeArray(1, prod->getSimSeries()->getLastDate()); // expiry

        cliquetTimeline = DateTime::merge(dates);
        numCliquets = cliquetTimeline.size() - 1;
        if(!numCliquets) {
            // Something has gone wrong. It cannot be that we have zero number of cliquets
            // Maybe we are being overcautious here and the engine can cope with 1 AvgIn date
            // being equal to 1 AvgOut date. For the moment be defensive
            throw ModelException(
                "Simulation start date " + startDate.toString() +
                " and simulation end date " + dates[1][0].toString() +
                " coincide. They have to be distinct.");
        }

        for(iCliquet = 1; iCliquet < numCliquets; iCliquet++) {
            if(dates[0][iCliquet] <= startDate) {
                throw ModelException("Failed to construct the cliquet timeline. "
                                     "More than 1 live cliquet start dates are before the start date.");
            }
        }

        // construct VolRequests per cliquet, asset
        volRequests = vector< vector<VolRequestLNStrikeSP> >(numCliquets);
        for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            volRequests[iCliquet] = vector<VolRequestLNStrikeSP>(numAssets);
        }

        // Get VolRequests or create them
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            CVolRequestLNArray prodVolRequests(prodLN->getVolInterp(pastPathGenerator.get(), iAsset));
            if (prodVolRequests.empty()){
                throw ModelException(routine, "No vol requests specified");
            }
            CVolRequestLNSP volReq(prodVolRequests[0]);
            CliquetVolRequest*  cliquetVol  = dynamic_cast<CliquetVolRequest*>(volReq.get());
            VolRequestLNStrike* volReqLNStr = dynamic_cast<VolRequestLNStrike*>(volReq.get());

            if(cliquetVol) {
                CVolRequestLNArray productVolRequests = cliquetVol->getRequestsArray();
                for(iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
                    VolRequestLNStrike* volLNStrike = dynamic_cast<VolRequestLNStrike*>(productVolRequests[iCliquet].get());
                    if(volLNStrike) {
                        volRequests[iCliquet][iAsset] = VolRequestLNStrikeSP(copy(volLNStrike));
                    } else {
                        throw ModelException(routine, "Vol interp is not LNStrike "
                                             "type for asset " +
                                             mAsset->assetGetName(iAsset));
                    }
                }
            } else if(volReqLNStr) {
                volRequests[0][iAsset] = VolRequestLNStrikeSP(copy(volReqLNStr));
            } else {
                throw ModelException(routine, "Vol interp is neither LNStrike "
                                     "type nor Cliquet type for asset " +
                                     mAsset->assetGetName(iAsset));
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

// implement getGaussData since derivation from DependenceMakerGauss::Support
CDoubleMatrixConstSP MCPathBaseImplied::getGaussData() const {
    static const string method("MCPathBaseImplied::getGaussData");
    try {
        return mAsset->factorsCorrelationMatrix();
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

// implement getSimDates since derivation from DependenceMakerGaussTerm::Support
DateTimeArray MCPathBaseImplied::getSimDates() const {
    return timeline->futureDates;
}

// implement getFwdVarAtDates since derivation from DependenceMakerGaussTerm::Support
DoubleMatrix MCPathBaseImplied::getFwdVarAtDates(bool interpolateAtmFwd) const {
    // non SV approach: do not pass today to Driver
    DateTime valueDate = timeline->futureDates.front();
    DateTimeArray toDates = valueDate.getFutureDates(timeline->futureDates);
    if (interpolateAtmFwd) {

        // Driver object retrieves sqrtFwdVars for driver from cache
        const CliquetCacheArray* tmpCCA = impliedCache->cliquetCaches.get();
        MCPathConfigImplied::PathGenSpot::Driver driver(valueDate.getFutureDates(cliquetTimeline),
                                                        toDates,
                                                        driverPath,
                                                        *tmpCCA,
                                                        cliquetArray);
        return DoubleMatrix(driver.driverVars);  // sqrtFwdVars at simDates
    } else {
        // interpolate atm spot as opposed to atm fwd
        DoubleArrayArray var(numAssets);
        for(int iAsset=0; iAsset<numAssets; iAsset++ ) {
            var[iAsset] = DoubleArray(toDates.size());
            const CAsset& thisAsset = mAsset->getAsset(iAsset);
            ATMVolRequestSP volRequest(new ATMVolRequest());
            volRequest->allowNegativeFwdVar(false);
            CVolProcessedBSSP volBS(thisAsset.getProcessedVol(volRequest.get()));
            volBS->CalcVar(valueDate,
                           toDates,
                           volBS->forward,
                           var[iAsset]);
        }
        return DoubleMatrix(var);
    }
}

// implement getMultiFactors since derivation from DependenceMakerGaussTerm::Support
const IMultiFactors* MCPathBaseImplied::getMultiFactors() const {
    return mAsset;
}

// implement getDependenceMaker since derivation from DependenceMakerLocalCorr::Support
DependenceMakerSP MCPathBaseImplied::getDependenceMaker() const {
    return dependenceMaker;
};

// implement getNbPastDates since derivation from DependenceMakerLocalCorr::Support
int MCPathBaseImplied::getNbPastDates() const {
    return timeline->totalNumPastDates;
}

// implement getRandomGenerator since derivation from DependenceMakerLocalCorr::Support
const IRandomSP MCPathBaseImplied::getRandomGenerator() const {
    return rand;
}

// implement getFactorRandomCache since derivation from DependenceMakerLocalCorr::Support
const MCPathConfig::RandomCacheSP MCPathBaseImplied::getIdiosynFactorRandomCache() const {
    return idiosynFactorRandomCache;
}

// implement getFactorRandomCache since derivation from DependenceMakerLocalCorr::Support
const MCPathConfig::RandomCacheSP MCPathBaseImplied::getMarketFactorRandomCache() const {
    return marketFactorRandomCache;
}

// implement
bool MCPathBaseImplied::carefulRandoms() const {
    return isCarefulRandoms;
}

/** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
    the 'drift' between simulation date i and simulation date i+1.
    (With d_0 = simulation start date). Range of dates is restricted to
    those in the future so the past path generator will
    return 1.0. The 'drift' is the biggest ratio of the 2 spots possible
    within some sort of probability interval. See  MCPathConfigLN */
double MCPathBaseImplied::maxDriftProduct(int iAsset) const {
    static const string routine("MCPathBaseImplied::maxDriftProduct");

    try {
        // Return the maximum maxDrift over all possible cliquets
        double maxDrift = 0.0;
        for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            if(cliquetArray[iCliquet]->maxDrifts.get()) {
                maxDrift = Maths::max(maxDrift,
                                      cliquetArray[iCliquet]->maxDrifts->getMaxDrift(iAsset));
            } else {
                throw ModelException("Unable to get maxDrifts for cliquet number "+
                                     Format::toString(iCliquet+1) +
                                     " for asset " +
                                     mAsset->assetGetName(iAsset));
            }
        }
        return maxDrift;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/** Returns number of assets */
int MCPathBaseImplied::NbSimAssets() const {
    return numAssets;
}

/** Returns the reference level for iAsset, iPath */
double& MCPathBaseImplied::refLevel(int iAsset, int iPath) {
    return refData->refLevels[iAsset][0];
}

/** Returns the reflevel path */
IRefLevel::IMCPathSP& MCPathBaseImplied::refLevelPath() const {
    return refData->refLevelPath;
}

/** Returns the paths */
double* MCPathBaseImplied::Path(int iAsset, int iPath) {
    return &productPaths[iAsset][0];
}



/** Generate path for specified path in simulation (pathIdx), for
    specified asset (iAsset) and specified vol interp (iPath) */
void MCPathBaseImplied::generatePath(
    int pathIdx, int iAsset, int iPath){
    // populate paths field
    double* assetCumul = cumulatives[iAsset];
    double* assetPath  = &productPaths[iAsset][0]; // for ease/speed

    for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
        // extract information for the cliquet
        const Cliquet& cliquet            = *cliquetArray[iCliquet];
        int numFutSteps                   = cliquet.numFutSteps;
        const ImpliedSampleArray& samples = cliquet.samples[iAsset];

        // Some pointer arithmetic to ensure we fill in the path properly
        int cliquetOffset = cliquetDiffOffsets[iCliquet];
        assetCumul += cliquetOffset;
        assetPath  += cliquetOffset;
        // The path in particular must be offset for the past for
        // the first cliquet. The remaining offsets are relative
        // so we should only do it once for the first cliquet
        if (iCliquet == 0) {
            assetPath += timeline->numPastDates;
        }

        // Get the reference level in advance
        double assetRefLevel = *(cliquet.returnReference[iAsset]);

        // Loop over sample dates and get a sample
        const int* isRelative = &cliquet.isRelativeReturn[0];
        for(int iStep = 0; iStep < numFutSteps; iStep++) {
            // Sample and multiply by reference if needed
            double assetSpot = samples[iStep]->sample(assetCumul[iStep]);
            if(isRelative[iStep]) {
                assetSpot *= assetRefLevel;
            }
            assetPath[iStep] = assetSpot;
        }
    }

#if PATH_DEBUG
    // prints the path in a file
    if (iAsset == numAssets){
        ofstream  debugpathsfile("c:/paths.txt", ios_base::app);
        debugpathsfile << "Iteration : " << pathIdx << "\n";
        for (int j = 0; j < numAssets; j++) {
            debugpathsfile << "Asset : " << j << "\n";
            for (int i = 0; i < this->numFutSteps; i++) {
                debugpathsfile << paths[j][0][i] << "\n";
            }
        }
    }
#endif
}

void MCPathBaseImplied::drawRandomNumbers(int pathIdx) {
    // Generate random numbers
    randomGen->generate(pathIdx);
}



void MCPathBaseImplied::configureAntithetics() {
    // Notice that there is some overhead as the MCRandom class has
    // already negated the random numbers

    // In the longer run, we might want to disable fast antithetics
    // because we are losing the driver values. We will need them for BB.

    // Strictly speaking, this pathGen should not be using an MCRandom
    // random number generator but a MCRandomUniform (to be coded).
    // The latter should do antithetics by u(anti) = 1.0 - u
    cumulatives.negate();
    cumulatives.scalarAdd(1.0);
}

void MCPathBaseImplied::configureNonAntithetics() {
    // Given the random numbers, compute the driver path
    // and cumulative distributions
    const DoubleMatrix& randoms = randomGen->getRandomNumbers();
    for (int iAsset = 0; iAsset < numAssets; iAsset++) {
        // for speed
        double* assetDriver        = driverPath[iAsset];
        double* assetCumul         = cumulatives[iAsset];
        const double* assetRandoms = randoms[iAsset];

        for(int iCliquet = 0; iCliquet < numCliquets; iCliquet++) {
            int numFutSteps = cliquetArray[iCliquet]->numFutSteps;

            // For easier access to caches
            const DriversSqrtVarSP& sqrtVars =
                (*(*impliedCache->cliquetCaches)[iCliquet]->assetCaches)[iAsset]->sqrtVars;

            // Some pointer arithmetic to ensure we fill in the path properly
            int cliquetOffset = cliquetDiffOffsets[iCliquet];
            assetDriver  += cliquetOffset;
            assetCumul   += cliquetOffset;
            assetRandoms += cliquetOffset;

            // Extract driver information for the cliquet
            const double* sqrtTotalVar = &(*sqrtVars->sqrtTotalVar)[0];
            const double* sqrtFwdVar   = &(*sqrtVars->sqrtFwdVar)[0];
            const int* zeroTT          = &(*sqrtVars->zeroTTime)[0];

            for(int iStep = 0; iStep < numFutSteps; iStep++) {
                if(zeroTT[iStep]) {
                    // Zero trading time up to here so the Brownian motion is still at 0.0
                    assetCumul[iStep] = 0.5;
                } else {
                    // Diffuse driver and compute cumulatives
                    if(iStep == 0) {
                        assetDriver[iStep] = 0.0;
                    } else {
                        assetDriver[iStep] = assetDriver[iStep - 1];
                    }
                    assetDriver[iStep] += sqrtFwdVar[iStep] * assetRandoms[iStep];
                    assetCumul[iStep] = tabulatedNorm->value(assetDriver[iStep] /
                                                             sqrtTotalVar[iStep]);
                }
            }
        }
    }
}



DRLIB_END_NAMESPACE
