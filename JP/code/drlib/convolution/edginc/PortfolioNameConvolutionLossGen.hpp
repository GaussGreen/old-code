//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 04-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_PORTFOLIONAMECONVOLUTIONLOSSGEN_HPP
#define QLIB_PORTFOLIONAMECONVOLUTIONLOSSGEN_HPP

#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"
#include "edginc/PortfolioName.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Loss generator implementation for PortfolioName.
 * Needs all relevant "loss generators" implementations since
 * the "outer model" may use any of them.
 * */
class CONVOLUTION_DLL PortfolioNameConvolutionLossGen:
//        public virtual IEffectiveLossCurveGen, TODO
    public virtual ILossDistributionsGen,
    public virtual ICondLossDistributionsGen
{
public:

    virtual ~PortfolioNameConvolutionLossGen();

    PortfolioNameConvolutionLossGen(
        const PortfolioName* ptfName,
        const IConditionalDefaultsModel* condDefaultsModel);

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
    PortfolioNameConvolutionLossGen(const PortfolioNameConvolutionLossGen& rhs); // copy constructor - not defined
    PortfolioNameConvolutionLossGen& operator=(const PortfolioNameConvolutionLossGen& rhs); // operator= method - not defined

    // fields
    const PortfolioName* ptfName;
    const IConditionalDefaultsModel* condDefaultsModel;
};

DRLIB_END_NAMESPACE

#endif /*PORTFOLIONAMECONVOLUTIONLOSSGEN*/
