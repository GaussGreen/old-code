//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 04-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_BCTRANCHECONVOLUTIONLOSSGEN_HPP
#define QLIB_BCTRANCHECONVOLUTIONLOSSGEN_HPP

#include "edginc/IModelConfigMapper.hpp"
#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IBCCondLossDistributionsGen.hpp"
#include "edginc/IIntegrator.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/IMarketFactorModel.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/IndexWeights.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Loss generator implementation for CreditTrancheLossConfig using
 * Base Correlation methodology.
 * Needs all relevant "loss generators" implementations since
 * the "outer model" may use any of them.
 * */
class CONVOLUTION_DLL BCTrancheConvolutionLossGen:
//        public virtual IEffectiveLossCurveGen, TODO
    public virtual ILossDistributionsGen,
    public virtual ICondLossDistributionsGen
{
public:

    virtual ~BCTrancheConvolutionLossGen();

    BCTrancheConvolutionLossGen(
        const ITranchesCombinationPayoff* tranchesCombinationPayoff,
        double lossLevel,
        const ICreditLossConfig* innerLossConfig,
        const DateTime& valueDate,
        const IModelConfigMapper* mapper,
        const IConvolutor* convolutor,
        bool returnBinaryDistribution,
        const IConditionalDefaultsModel* condDefaultsModel,
        double compressionRatio,
        CDoubleConstSP strikeMappingOverride,
        DoubleArrayConstSP lowerBCBetasOverride,
        DoubleArrayConstSP upperBCBetasOverride,
        DateTimeArrayConstSP bcBetasTimeline,
        CDoubleConstSP spreadRatioOverride,
        DoubleArrayConstSP spreadRatiosRegionalOverride,
        StringArrayConstSP spreadRatiosRegions,
        double historicalBetaAdjustment,
        StringArrayConstSP namesOverride,
        IndexWeightsArrayConstSP indexWeightsOverride,
        bool authoriseNegativeEL,
        bool useExpectedLossRatio,
        YieldCurveConstSP discountCurve);

//        /** Implements "IEffectiveLossCurveGen" */ TODO
//        virtual void createEffectiveCurves(
//            // inputs
//            const DateTimeArray&                timeline, // still needed ??
//            CounterPartyCreditConstSP           cpty,
//            Control*                            control,
//            Results*                            results,
//            bool                                recoverNotional,
//            // outputs
//            IDiscountCurveRiskySP& ctgEffCurve,
//            IDiscountCurveRiskySP& feeEffCurve,
//            IDiscountCurveRiskySP& ctgEffCurveCond,
//            IDiscountCurveRiskySP& feeEffCurveCond) const;

    /** [Implements ILossDistributionsGen] */
    virtual IDistribution1DArraySP createLossDistributions(
        const DateTimeArray& timeline) const;
        
    /** [Implements ICondLossDistributionsGen] */
    virtual ICondLossDistributionsGenKeyArrayConstSP initialise(
        const DateTimeArray& timeline) const;

private:
    BCTrancheConvolutionLossGen(const BCTrancheConvolutionLossGen& rhs); // copy constructor - not defined
    BCTrancheConvolutionLossGen& operator=(const BCTrancheConvolutionLossGen& rhs); // operator= method - not defined

    // definition in source file
    struct IndexData;

    typedef refCountPtr<IndexData> IndexDataSP;
    
    /** A (portfolio name, index weights) map*/
    typedef map<string, IndexWeightsConstSP> IndexWeightsMap;

    /** A (index name, IndexData) map */
    typedef map<string, IndexDataSP> IndexDataMap;

    friend class BCSkewCalculator; // to allow access to IndexDataMap structure
    
    /** Initialises indexDataMap and computes sumOutNotional */
    void initIndexDataMap(
        const DateTime& today,
        IBCCondLossDistributionsGenArrayConstSP lossGens,
        YieldCurveConstSP discount,
        const DateTime& timepoint,
        double pastPortLoss,
        IndexDataMap& indexDataMap) const;

    /** Auxiliary function to compute sumOutNotional and histAvgBeta */
    void computeSumNotionalAndBeta(
        IBCCondLossDistributionsGenArrayConstSP lossGens,
        double& sumOutPosNotional,
        double& sumOutNegNotional,
        double& sumOutNetNotional,
        double& histAvgBeta) const;

    // fields
    
    const ITranchesCombinationPayoff* tranchesCombinationPayoff;
    // lossLevel is a unique possible loss level associated to "tranchesCombinationPayoff"
    double lossLevel;
    const ICreditLossConfig* innerLossConfig;
    
    const IModelConfigMapper* mapper;
    const IConvolutor* convolutor;
    bool returnBinaryDistribution;
    const IConditionalDefaultsModel* condDefaultsModel;
    
    DateTime valueDate;

    // BC fields
    double compressionRatio;
    CDoubleConstSP strikeMappingOverride;
    DoubleArrayConstSP lowerBCBetasOverride;
    DoubleArrayConstSP upperBCBetasOverride;
    DateTimeArrayConstSP bcBetasTimeline;
    CDoubleConstSP spreadRatioOverride;
    DoubleArrayConstSP spreadRatiosRegionalOverride;
    StringArrayConstSP spreadRatiosRegions;
    double historicalBetaAdjustment;
    StringArrayConstSP namesOverride;
    IndexWeightsArrayConstSP indexWeightsOverride;
    bool authoriseNegativeEL;
    bool useExpectedLossRatio;
    YieldCurveConstSP discountCurve;
};

DRLIB_END_NAMESPACE

#endif /*BCTRANCHECONVOLUTIONLOSSGEN*/
