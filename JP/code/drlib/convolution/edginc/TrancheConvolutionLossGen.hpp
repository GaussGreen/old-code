//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 04-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHECONVOLUTIONLOSSGEN_HPP
#define QLIB_TRANCHECONVOLUTIONLOSSGEN_HPP

#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"
#include "edginc/IModelConfigMapper.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Loss generator implementation for CreditTrancheLossConfig and more
 * generic payoffs on the loss distribution.
 * Needs all relevant "loss generators" implementations since
 * the "outer model" may use any of them.
 * */
class CONVOLUTION_DLL TrancheConvolutionLossGen:
//        public virtual IEffectiveLossCurveGen, TODO
    public virtual ILossDistributionsGen,
    public virtual ICondLossDistributionsGen
{
public:
    
    virtual ~TrancheConvolutionLossGen();

    TrancheConvolutionLossGen(
        const ICreditLossConfig* lossConfig,
        double lossLevel,
        const ICreditLossConfig* innerLossConfig,
        const IModelConfigMapper* mapper,
        const IConvolutor* convolutor,
        bool returnBinaryDistribution,
        const IConditionalDefaultsModel* condDefaultsModel,
        const DateTime& valueDate);

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
    virtual ICondLossDistributionsGenKeyArrayConstSP initialise(const DateTimeArray& timeline) const;

private:
    TrancheConvolutionLossGen(const TrancheConvolutionLossGen& rhs); // copy constructor - not defined
    TrancheConvolutionLossGen& operator=(const TrancheConvolutionLossGen& rhs); // operator= method - not defined

    // fields
    
    // lossConfig can be a tranche or something more generic (payoff on the loss distribution)
    const ICreditLossConfig* lossConfig;
    
    // lossLevel is a unique possible loss level associated to "lossConfig"
    double lossLevel;
    
    DateTime valueDate;
    
    // innerLossConfig is the "underlying" of lossConfig,
    // typically a CDOPortfolio but not necessary
    const ICreditLossConfig* innerLossConfig;
    
    const IModelConfigMapper* mapper;
    const IConvolutor* convolutor;
    bool returnBinaryDistribution;
    const IConditionalDefaultsModel* condDefaultsModel;
};

DRLIB_END_NAMESPACE

#endif /*TRANCHECONVOLUTIONLOSSGEN*/
