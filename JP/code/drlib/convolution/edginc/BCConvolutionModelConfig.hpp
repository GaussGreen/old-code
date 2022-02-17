//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 04-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_BCCONVOLUTIONMODELCONFIG_HPP
#define QLIB_BCCONVOLUTIONMODELCONFIG_HPP

#include "edginc/IModelConfigMapper.hpp"
#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IBCCondLossDistributionsGen.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/IndexWeights.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Model config implementation for loss configs from which we
 * can get a loss distribution using a convolution algorithm and
 * the "Base Correlation" methodology
 * */
class CONVOLUTION_DLL BCConvolutionModelConfig:
    public CObject,
    public virtual ICreditLossModelConfig,
    public virtual IModelConfigMapper
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~BCConvolutionModelConfig();

    /** Implementation of ICreditLossModelConfig */
    virtual ICreditLossGenSP lossGenerator(
        ICreditLossConfigConstSP lossConfig,
        IModelConfigMapperConstSP mapper) const;

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP effCurveGenerator(
        ICreditLossConfigConstSP  ntdLossCfg,
        CounterPartyCreditConstSP cpty,
        const bool                recoverNotional,
        IModelConfigMapperConstSP modelConfigMapper) const;

    /** Dummy implementation of IModelConfigMapper (just returns "this") */
    virtual ICreditLossModelConfigConstSP innerModel(
        ICreditLossConfigConstSP lossConfig) const;

    /** Explicit constructor */
    BCConvolutionModelConfig(
        IConditionalDefaultsModelSP condDefaultsModel,
        IConvolutorSP convolutor,
        bool returnBinaryDistribution,
        double compressionRatio,
        CDoubleSP strikeMappingOverride,
        DoubleArraySP lowerBCBetasOverride,
        DoubleArraySP upperBCBetasOverride,
        DateTimeArraySP bcBetasTimeline,
        CDoubleSP spreadRatioOverride,
        DoubleArraySP spreadRatiosRegionalOverride,
        StringArraySP spreadRatiosRegions,
        double historicalBetaAdjustment,
        StringArraySP namesOverride,
        IndexWeightsArraySP indexWeightsOverride,
        bool authoriseNegativeEL,
        bool useExpectedLossRatio);

    // Hack to pass discount curve - note the input discount curve is cloned
    void setDiscountCurve(YieldCurveConstSP discount);

private:
    BCConvolutionModelConfig(const BCConvolutionModelConfig& rhs); // copy constructor - not defined
    BCConvolutionModelConfig& operator=(const BCConvolutionModelConfig& rhs); // operator= method - not defined

    BCConvolutionModelConfig();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    /** Convolution engine */
    IConvolutorSP convolutor;

    // Flag to control if we want to approximate the returned
    // loss distribution by a "binary distribution" (in which case
    // we don't need to do a full convolution, we can just "convolute and
    // integrate").
    bool returnBinaryDistribution;

    /** Conditional default model */
    IConditionalDefaultsModelSP condDefaultsModel;

    // -------------------------------
    // BC specific fields / parameters
    // -------------------------------

    /** Default value for compression ratio (theta) */
    static double const DEFAULT_THETA;
    
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
    CDoubleSP strikeMappingOverride;
     
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
    CDoubleSP spreadRatioOverride;
    
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
    
    // TODO - change ?
    // Hack to have access to this instrument yield curve.
    // An alternative would be to pass the discount curve in the
    // ICondLossDistributionsGen and ILossDistributionsGen methods, but
    // seems heavy since the discountCurve is only used for BC with 
    // wide spreads and "spread ratio" (not used with "expected loss" ratio).
    YieldCurveSP discountCurve;
};

DECLARE(BCConvolutionModelConfig);

DRLIB_END_NAMESPACE

#endif /*BCCONVOLUTIONMODELCONFIG*/
