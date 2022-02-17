//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CreditMetricsBaseCorrelation.hpp
//
//   Description : Credit Metrics model with Base Correlation
//
//   Author      : Antoine Gregoire
//
//   Date        : May 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDIT_METRICS_BASE_CORRELATION_HPP
#define QLIB_CREDIT_METRICS_BASE_CORRELATION_HPP

#include "edginc/CreditMetricsModel.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/CompressionRatioTweak.hpp"
#include "edginc/StrikeMappingOverrideTweak.hpp"
#include "edginc/BCStrikeMappingOverrideTweak.hpp"
#include "edginc/QuasiContractualBaseCorrelation.hpp"
#include "edginc/SkewSurface.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IndexWeights);

/** Credit Metrics model with Base Correlation */
class PRODUCTS_DLL CreditMetricsBaseCorrelation :
    public CreditMetricsModel,
    virtual public TweakableWith<CompressionRatioTwk>,
    virtual public TweakableWith<StrikeMappingOverrideTwk>,
    virtual public TweakableWith<BCStrikeMappingOverrideTwk>,
    virtual public QuasiContractualBaseCorrelation::IShift
{
    struct IndexData; // definition in source file
    typedef refCountPtr<IndexData> IndexDataSP;
    
    /** A (portfolio name, index weights) map*/
    typedef map<string, IndexWeightsConstSP> IndexWeightsMap;

    /** A (index name, IndexData) map */
    typedef map<string, IndexDataSP> IndexDataMap;

public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CreditMetricsBaseCorrelation();
    
    /** public default constructor */
    CreditMetricsBaseCorrelation();

    /** Generate the timeline on which the effective curve will be calculated.
        This overrides the one CreditMetricsModel in order to do some 
        additional validation */
    DateTimeArraySP generateTimeline(
        const DateTime& today,
        const DateTime& lastObservationDate) const;

    /** Overridden to use "Base Correlation" methodology. Note that a
        IFixedTrancheLossCalculator is returned via the conditionalLossCalc 
        parameter. */
    virtual IFixedTrancheLossCalculator* createFixedLossCalculator(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const; // (O)

    /** Analagous to createFixedLossCalculator, but for recovered notional.
        Overridden to use "Base Correlation" methodology. Note that a
        IFixedTrancheLossCalculator is returned via the conditionalLossCalc 
        parameter. */
    virtual IFixedTrancheLossCalculator* createFixedRecoveredNotionalCalculator(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const; // (O)

    /** Invoked by the containing model after the instrument data is fetched
        ie after CInstrument::GetMarket is invoked. This invokes parent method
        and retrieves indexWeightsOverride data */
    void postFetchMarketData(IModel*            model,
                             MarketDataConstSP market);
    
    /** Called immediately after object constructed */
    void validatePop2Object();

    /** Returns all the points on the skew surface (of given name) to which this
        product is sensitive. This implementation returns the points used by the
        Base Correlation methodology */
    virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
        OutputNameConstSP              outputName,
        const DateTime&                lastObservationDate, 
        CreditTrancheLossConfigConstSP tranche,  
        YieldCurveConstSP              discount,
        bool                           tweakAll, 
        const int                      numberOfNeighbours) const;

    /** Compression ratio tweak support */
    virtual string sensName(CompressionRatioTwk* shift) const;
    /** Compression ratio tweak support */
    virtual bool sensShift(CompressionRatioTwk* shift);

    /** Strike mapping override tweak support */
    virtual string sensName(StrikeMappingOverrideTwk* shift) const;
    /** Strike mapping override tweak support */
    virtual bool sensShift(StrikeMappingOverrideTwk* shift);

    /** Base Correlation Strike mapping override tweak support */
    virtual string sensName(BCStrikeMappingOverrideTwk* shift) const;
    /** Base Correlation Strike mapping override tweak support */
    virtual bool sensShift(BCStrikeMappingOverrideTwk* shift);
    
    /** Compression ratio and betaBasis setting support */
    virtual bool sensShift(QuasiContractualBaseCorrelation* shift);

    map<string, double> getDurationWeightedAverage(
        CreditTrancheLossConfigConstSP tranche,
        const DateTime&                maturityDate,
        const YieldCurveConstSP        discount) const;

    /** CAUTION: the spread ratios etc will NOT be computed for the first
        point of the timeline (typically today)
        - includeName: if a NULL SP is passed in, ALL (non-defaulted) names
          will be included
        - indicesToIgnore: if a NULL SP is passed in, NO indices will be 
          ignored */
    void populateDurWeightAvgSpreadData(
        CreditTrancheLossConfigConstSP tranche,
        CBoolArraySP                   includeName,
        const YieldCurveConstSP        discount,
        const DateTimeArray&           reducedTimeLine,
        IndexWeightsMap&               indexWeightsMap,
        StringArrayConstSP             indicesToIgnore,
        IndexDataMap&                  indexDataMap) const;  // (O)
        

protected:
    /** Only build instances of that class using reflection */
    CreditMetricsBaseCorrelation(const CClassConstSP& clazz);
    
private:
    DoubleArrayArraySP calculateBCBetas(
        const DateTimeArray&                timeline,        // includes today
        CreditTrancheLossConfigConstSP      tranche,         // view onto Inst
        OutputRequest*                      betaRequest,     // may be null
        Results*                            results,         //for output result
        double                              strike,          // eg lowStrike
        DoubleArrayConstSP                  bcBetasOverride, //optional override
        SkewSurface::SkewType               skewType,        // fast or full
        double                              pastPortLoss,    // portfolio losses
        double                              sumOutPosNotional,
		double								sumOutNegNotional,
        double                              histAvgBeta,
        double                              avgBeta,         // average beta
        DoubleArrayConstSP                  adjBetas,        // dispersionAdjust
        const DoubleArray&                  timelineAsDouble,// dates->doubles
        const IndexDataMap&                 indexDataMap) const;// from initMaps

    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Returns the indexWeightsMap for all names in the tranche, using
        the override if present */
    void getIndexWeightsMap(
        CreditTrancheLossConfigConstSP tranche,                // (I)
        IndexWeightsMap&               indexWeightsMap) const; // (O)

    void initMaps2(
        CreditTrancheLossConfigConstSP tranche,
        YieldCurveConstSP              discount,
        const DateTimeArray&           timeLine,
        const DoubleArrayArray&        namesSurvivalProb,
        double                         pastPortLoss,
        IndexWeightsMap&               indexWeightsMap,
        IndexDataMap&                  indexDataMap,
        double&                        sumOutPosNotional,
        double&                        sumOutNegNotional,
        double&                        sumOutNetNotional,
        double&                        histAvgBeta) const;


/***** JLH Not used, so remove it - ? */
//     DoubleArrayArraySP ccmBespokeSkew2(
//         const ConvolutionProduct* product,
//         const DateTimeArray&      timeLine,
//         const DoubleArrayArray&   namesSurvivalProb,
//         const DoubleArray&        strikes,
//         SkewSurface::SkewType     skewType) const;

    double calculateAverageBeta2(CreditTrancheLossConfigConstSP tranche) const;

    DoubleArraySP dispersionAdjust2(
        CreditTrancheLossConfigConstSP tranche,
        double                         avgBeta) const;


    /** Calculate the bespoke portfolio skew information for one given strike */
    // NB : we assume lgdFloor = 0.0, lgdCap = 1.0, lgdNotional = 1.0
    // This is checked in CDO when pricing with BaseCorrelation
    void ccmBespokeSkewOneStrike(
        const DateTimeArray&  timeLineEff,
        double                pastPortLoss,       
        double                strike,
        SkewSurface::SkewType skewType,
        const IndexDataMap&   indexDataMap,
        double                sumOutPosNotional,
		double				  sumOutNegNotional,
        double                histAvgBeta,
        DoubleArray&          result) const;

    /** */
    BetaSkewGridPointSet getSensitiveBetaSkewPointsOneStrike(
        const DateTimeArray& timeLine,
        double               pastPortLoss,       
        double               strike,
        IndexDataSP          indexData,
        double               sumOutPosNotional,
        string               surfaceName, 
        const int            numberOfNeighbours) const;

    // ---------
    // CONSTANTS
    // ---------
    
    /** Default value for compression ratio (theta) */
    static double const DEFAULT_THETA;
    
    /** Name of the (optional) beta basis */
    static string const BETA_BASIS_NAME;
    
    /** Default strike used in the calculation of CCC */
    static double const EQUITY_STRIKE_FOR_CCC;
    
protected:
    // -----------------------------------
    // OPTIONAL FIELDS ('MARKET' OVERRIDE)
    // -----------------------------------

    /** 
     * Compression ratio (also called theta) used to rescale the
     * historical betas dispersion.
     * Also called "theta".
     * Range [0,1]
     * Default value is DEFAULT_THETA
     * */
    double compressionRatio;

    /** 
     * Optional global parameter to override strike mapping from market 
     * (global override for all strike mapping parameters in indexSkews).
     * Also called "q".
     * Range [0,1]
     * */
     // use 'double*' instead of 'double' to detect if
     // this field has been populated
    double* strikeMappingOverride;
     
    /** Optional override for lower strike base correlation betas */
    DoubleArraySP lowerBCBetasOverride;
    
    /** Optional override for upper strike base correlation betas */
    DoubleArraySP upperBCBetasOverride;
    
    /**
     * Timeline associated with lowerBCBetasOverride and upperBCBetasOverride.
     * This input is used for checking purpose only (bcBetasTimeline must be
     * exactly the same as the effective curve timeline computed internally)
     * */
    DateTimeArraySP bcBetasTimeline;
    
    /** 
     * Optional global parameter to override the spread ratio used for that trade,
     * typically use 1 for index trades
     * */
     // use 'double*' instead of 'double' to detect if
     // this field has been populated
    double* spreadRatioOverride;
    
    /**
     * Optional override: allow to specify a "local" spread ratio override per
     * region (corresponding regions are defined in "spreadRatiosRegions").
     * Note that "spreadRatioOverride" and "spreadRatiosRegionalOverride" are
     * exclusive options (if both are specified at the same time, we return an
     * error).
     * */
    DoubleArraySP spreadRatiosRegionalOverride;
    
    /**
     * See spreadRatiosRegionalOverride
     * */
    StringArraySP spreadRatiosRegions;
        
    // ---------------------
    // OTHER OPTIONAL FIELDS
    // --------------------- 
    
    /** 
     * Optional parameter used to offset the skew surface.
     * Default is a flat 0 surface.
     * */
    SkewSurfaceSP betaBasis;
    
    /**
     * Historical Beta adjustment.
     * Default value is 0.
     * */
    double historicalBetaAdjustment;
    
    /**
     * Array of names (in the portfolio) for which we want to 
     * override the associated indexWeights
     * */
    StringArraySP namesOverride;
    
    /**
     * Array of overriden indexWeights
     * */
    IndexWeightsArraySP indexWeightsOverride;
    
    /**
     * Flag to authorise negative expected losses (=TRUE).
     * Otherwise negative expected losses will be floored to 0.
     * Default is FALSE (i.e. do NOT authorise negative expected losses).
     * */
    bool authoriseNegativeEL;

    /**
     * Flag to use expected loss ratio (=TRUE).
     * instead of par spread ratio in the spread mapping.
     * Default is FALSE (par spread mapping).
     * */
    bool useExpectedLossRatio;
};

DRLIB_END_NAMESPACE

#endif //QLIB_CREDIT_METRICS_BASE_CORRELATION_HPP
