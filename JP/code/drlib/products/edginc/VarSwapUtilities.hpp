//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VarianceSwapUtilities.hpp
//
//   Description : Variance Swap Utilities
//
//   Author      : Manos Venardos
//
//   Date        : 1 February 2006
//
//
//----------------------------------------------------------------------------

#ifndef VARSWAP_UTILITIES_HPP
#define VARSWAP_UTILITIES_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Function.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/Integrator.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/MatrixResult.hpp"
#include "edginc/ModelLN.hpp"

DRLIB_BEGIN_NAMESPACE

/** Main options recorder */
FORWARD_DECLARE(VanillaInfo);
FORWARD_DECLARE(VanillaContractsRecorder);
FORWARD_DECLARE(VegaParallel);
FORWARD_DECLARE(VegaMatrix);
FORWARD_DECLARE(VegaMatrixLite);


class PRODUCTS_DLL VanillaInfo: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Empty constructor */
    VanillaInfo();

    /** Full constructor */
    VanillaInfo(const DateTime& maturity,
                double nbContracts,
                bool isCall,
                double fwd,
                double strike,
                double vol,
                double yearFrac);

    /** Scales results */
    void scaleNbContracts(double scalingFactor);

    /** Compute price of vanilla portfolio */
    double price() const;
    
    /** Compute vega of vanilla option defined as
        Vega = (Price(vol + shift) - Price(vol)) / shift */
    double vega(double shift) const;

    /** Returns indices such that 
        1) x < xArray[0]            ---> lowerBound = upperBound = 0 
        2) x > xArray[N-1]          ---> lowerBound = upperBound = N-1
        3) x[0] < x < xArray[N-1]   ---> lowerBound, upperBound such that x[lowerBound] < x < x[upperBound] */
    static void bracketArray(DoubleArray& xArray, double x, int& lowerBound, int& upperBound);
    
    /** VolSurface-type Implied Vol interpolation at a benchmark strike, between benchmark maturities */
    static double volInterpMaturity(double t1, double t2, 
                                    double vol_t1, double vol_t2,
                                    double t);

    /** VolSurface-type Implied Vol interpolation at a benchmark maturity, between benchmark strikes */
    static double volInterpStrike(double k1, double k2, 
                                  double vol_k1, double vol_k2,
                                  double k);
    
    /** VolSurface-type Implied Vol interpolation between benchmark maturities, strikes*/
    static double volInterp(double t1, double t2, double k1, double k2,
                            double vol_t1_k1, double vol_t1_k2, double vol_t2_k1, double vol_t2_k2,
                            double t, double k);
    
    /** Compute price of vanilla portfolio */
    static double price(VanillaContractsRecorderSP recorder);

    /** Compute vega of vanilla portfolio */
    static double vegaParallel(VegaParallelSP sens, VanillaContractsRecorderSP recorder);
    
    /** Compute vega matrix of vanilla portfolio */
    static void storeVegaMatrix(VegaMatrixLiteSP sens,
                                CInstrument* inst,
                                CModelLN* model,
                                VanillaContractsRecorderSP recorder,
                                CAssetSP asset,
                                const DateTime& valueDate,
                                const DateTime& instExpiry,
                                CResults* results);

private:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Empty shell method */
    static IObject* defaultVanillaInfo();
    
    DateTime     maturity;                //!< Maturity
    double       nbContracts;             //!< Number of contracts
    
    bool         isCall;                  //!< True: call, False: put
    double       fwd;                     //!< Forward price
    double       strike;                  //!< Absolute strike
    double       vol;                     //!< Implied Vol
    double       yearFrac;                //!< Year fraction to maturity
};


/////////////////////////////////////////////////////////////////////////////////////////////////////


/** Main options recorder */
class PRODUCTS_DLL VanillaContractsRecorder: public CObject {
public:
    static CClassConstSP const TYPE;
    
    /** Creates a recorder depending on the control */
    static VanillaContractsRecorderSP createVanillaOptionRecorder(const Control* control);

    /** Full constructor */
    VanillaContractsRecorder(double contractsScalingFactor, bool record);
    
    /** Records entry */
    void recordContract(VanillaInfoSP optionInfo);

    /** Records (T, K, nbOptions) entry */
    void recordContract(const DateTime& maturity,
                        double nbContracts,
                        bool isCall,
                        double fwd,
                        double strike,
                        double vol,
                        double yearFrac);
        
    /** Gives accress to scaling factor */
    double getScalingFactor() const;
    
    /** Changes scaling factor of nbContracts */
    void setScalingFactor(double newScalingFactor);

    /** Scales results */
    void scaleNbContracts(double x);

    /** Get results */
    VanillaInfoArraySP getResults() const;

private:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Default constructor */
    VanillaContractsRecorder();

    /** Empty shell method */
    static IObject* defaultVanillaContractsArray();
    
    double                      scalingFactor;               //!< Scaling factor to multiply all entries
    bool                        record;                      //!< True: record entries, False: ignore entries
    
    VanillaInfoArraySP          options;                     //!< Option array
};


DECLARE(VanillaContractsRecorder);


//////////////////////////////////////////////////////////////////////////////////////////////////////


/** Wrapper around BlackSholes */
double BlackPrice(bool isCall, double fwd, double strike, double pv, double vol, double yearFrac,
    double nbContracts, const DateTime& maturity, VanillaContractsRecorderSP recorder);


//////////////////////////////////////////////////////////////////////////////////////////////////////


/** Portfolio of continuum of vanillas at a given maturity: puts below fwd, calls above fwd */
class PRODUCTS_DLL StaticReplicationIntegrand: public Function1DDouble {
public:
    /** Full constructor */
    StaticReplicationIntegrand(Function1DDoubleConstSP      relativeStrikeWeight,
                               VanillaContractsRecorderSP   recorder,
                               CAssetConstSP                asset,
                               const DateTime&              valueDate,
                               const DateTime&              maturity,
                               bool                         negativeVar);

    /** Computes w(k) * PutCall(k) */
    virtual double operator()(double relativeStrike) const;

private:
    Function1DDoubleConstSP     relativeStrikeWeight;   //!< Strike weights
    VanillaContractsRecorderSP  recorder;               //!< Option recorder
    double                      fwd;                    //!< Forward
    DateTime                    valueDate;              //!< Value date
    DateTime                    maturity;               //!< Maturity
    CAssetConstSP               asset;                  //!< Asset
    LinearStrikeVolRequestSP    volRequest;             //!< VolRequest
    double                      yearFrac;               //!< Year frac from valueDate to maturity
};

typedef refCountPtr<StaticReplicationIntegrand> StaticReplicationIntegrandSP;


//////////////////////////////////////////////////////////////////////////////////////////////////////


/** Trapezium rule with n-points */
class PRODUCTS_DLL Trapez1DSimpleRecorder: public Trapez1DSimple {
public:
    static CClassConstSP const TYPE;

    Trapez1DSimpleRecorder(int nbSteps, VanillaContractsRecorderSP recorder);

    virtual double integrate(const Function1DDouble& func) const;
    
    double integrateAndRecord(const StaticReplicationIntegrand& func) const;

protected:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultTrapez1DSimpleRecorder();

    Trapez1DSimpleRecorder();

    VanillaContractsRecorderSP      recorder;       //!< Option recorder
};


//////////////////////////////////////////////////////////////////////////////////////////////////////


class PRODUCTS_DLL VarSwapUtilities {
public:
    /** Returns an array of log-returns i.e. first difference of input sample series.
    Properties of output array:
        1) Same size as input array
        2) Zero first element
        3) Zero for strictly future returns
        4) log-return up to current spot price if return is half-between past and future and
           computeCurrentReturn is set to true. If false then a zero will be computed
    
    Indexing is as follows
    S(0), S(1), S(2), S(3), ..., S(N)
      0,  y(1), y(2),   0, ...,    0
    where y(i) = log(S(i)/S(i-1)) if S(i) < today and 0 otherwise */
    static DoubleArraySP computeHistoricLogReturns(const CAsset*        asset,
                                                   const DateTimeArray& obsDates,
                                                   const DoubleArray&   obsSamples,
                                                   const DateTime&      valueDate,
                                                   bool                 computeCurrentReturn,
                                                   bool                 dividendAdjusted,
                                                   bool                 divAdjOnExDate);

    static double realisedVar(const CAsset*        asset,
                              const DateTimeArray& obsDates,
                              const DoubleArray&   obsSamples,
                              const DateTime&      valueDate,
                              bool                 computeCurrentReturn,
                              bool                 dividendAdjusted,
                              bool                 divAdjOnExDate);

    /** Computes cumulative historic dividends between 2 dates */
    static double getHistoricDivsBetweendates(DividendList::IDivAdjuster* divAdjuster,
                                              const CAsset* asset,
                                              const DateTime& valueDate,
                                              const DateTime& startDate,
                                              const DateTime& endDate);

    /** Computes variance due to future dividends: assuming discrete yields */
    static double futureDivsVariance(const CAsset*   asset,
                                     const DateTime& valueDate,
                                     const DateTime& firstDate,
                                     const DateTime& lastDate);
    
    /** Splits a given timeline into past and future samples according to VarSwap convention 
        i.e. futureExpectedN is the expectedN for a spot starting VarSwap until maturity
        (based on holiday schedule) and futureExpectedN + pastExpectedN = expectedN */
    static void pastAndfutureExpN(const DateTimeArray&  obsDates,
                                  const HolidayWrapper& assetHols,
                                  const DateTime&       valueDate,
                                  int                   expectedN,
                                  int&                  pastExpectedN,
                                  int&                  futureExpectedN);

    /** Fills in fixings for Theta shifts */
    static void thetaShiftCashFlows(const Theta*         shift,
                                    const CAsset*        asset,
                                    const DateTime&      valueDate,
                                    const DateTimeArray& obsDates,
                                    DoubleArray&         obsSamples);

    /** Fills in historic fixings from the asset history */
    static void populateSamples(const CAsset*               asset,
                                const ObservationSource*    assetHistorySource,                        
                                const DateTime&             valueDate,
                                const DateTimeArray&        obsDates,
                                const ObservationTypeArray& obsTypes,
                                DoubleArray&                obsSamples,
                                SamplingConventionSP        sampleRule = 
                                        SamplingConventionSP());
    
    static double futureDiscreteDivsAdjustment(const CAsset* asset,
                                               const DateTime& valueDate,
                                               const DateTime& startDate,
                                               const DateTime& endDate,
                                               int observationsPerYear, 
                                               int numSamplesInSwap);
};

DRLIB_END_NAMESPACE

#endif
