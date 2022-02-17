//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : View onto instrument as required by ConvolutionEngine
//
//   Date        : 18th Nov 2005
//
//   Author      : Mark Robson
//
//----------------------------------------------------------------------------

#ifndef QR_CONVOLUTIONMODEL_HPP
#define QR_CONVOLUTIONMODEL_HPP

#include "edginc/IModel.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ConvolutionProduct);
FORWARD_DECLARE(MarketDataFetcher);
FORWARD_DECLARE(IFixedTrancheLossCalculator);
FORWARD_DECLARE(MarketData);
FORWARD_DECLARE(BetaSkewGridPoint);
FORWARD_DECLARE(OutputName);
FORWARD_DECLARE(CreditTrancheLossConfig);
//FORWARD_DECLARE(CounterPartyCredit);
FORWARD_DECLARE(FlatCDO2LossConfig);


FORWARD_DECLARE_REF_COUNT(ICreditLossGen);
FORWARD_DECLARE(NToDefaultLossConfig);
class IModel;


/** Interface used by ConvolutionEngine (which is a CModel) to drive how
    the convolution is used in order to generated expected losses along a 
    timeline. */
class CONVOLUTION_DLL IConvolutionModel: public virtual IObject {
public:
    static CClassConstSP const TYPE;

    IConvolutionModel();
    virtual ~IConvolutionModel();

    /** Creates a MDF to be used when retrieving market data */
    virtual MarketDataFetcherSP createMDF() const = 0;

    /** Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this process (see IModel::wantsRiskMapping()) */
    virtual IModel::WantsRiskMapping wantsRiskMapping() const = 0;

    /** Invoked by the containing model before the instrument data is fetched
        ie before CInstrument::GetMarket is invoked. Allows this object to
        retrieve any market data it needs (but before the instrument data is
        retrieved). The model parameter will contain the MarketDataFetcher 
        created by the createMDF method above. Also see postFetchMarketData */
    virtual void preFetchMarketData(IModel*           model,
                                    MarketDataConstSP market) = 0;

    /** Invoked by the containing model after the instrument data is fetched
        ie after CInstrument::GetMarket is invoked. Allows this object to
        retrieve any market data it needs (but after the instrument data is
        retrieved). The model parameter will contain the MarketDataFetcher 
        created by the createMDF method above. Also see preFetchMarketData */
    virtual void postFetchMarketData(IModel*           model,
                                     MarketDataConstSP market) = 0;

    /** Generate the timeline on which the effective curve will be calculated */
    virtual DateTimeArraySP generateTimeline(
        const DateTime& today,
        const DateTime& lastObservationDate) const = 0;

	virtual void createLossCalculators(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,                   /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,      /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,              /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const = 0; //(O)

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP createEffCurveGenerator(
        NToDefaultLossConfigConstSP ntdLossCfg,
        CounterPartyCreditConstSP   cpty,
        const bool                  recoverNotional) const = 0;

	virtual void createLossCalculatorsFIX(
        const DateTimeArray&                timeline,           /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,                   /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,      /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,              /* (O) */
		IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc,
		FlatCDO2LossConfigConstSP      ts            /* (I) */) const;

    /** Indicates whether this model supports stochastic recovery rates
        AND there are any engine parameters for any names specifying so. */
    virtual const bool hasStochasticRecoveries(
        CreditTrancheLossConfigConstSP tranche) const = 0;

    /** Recovery of notional from the top of the portfolio requires
        an additional call to the convolution.
        It only has effect if the upper strike is greater than the
        sum of (name notional * name recovery).
        Therefore we can avoid making this additional call with the
        models assistance if the product displays the correct
        characteristics. Typically this means the model must be of
        a fixed recovery type */
    virtual const bool modelRecoveredNotional(
        CreditTrancheLossConfigConstSP tranche) const = 0; /* (I) */

    /** Returns how the effectiveCurve should be interpolated. This should
        really disappear when we return a DefaultRate object say in the 
        calculateEffectiveCurve() method. */
    virtual const string& getLossInterpolation() const = 0;

    /** Returns all the points on the skew surface (of given name) to which this
        product is sensitive. Perhaps this method should go in some derived
        class CONVOLUTION_DLL of IConvolutionModel?
        A null return value is ok - it is interpreted as the greek being 
        NotApplicable. The tweakAll flag, if true, means return all points
        on the skew surface */
    virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
        OutputNameConstSP              outputName,
        const DateTime&                lastObservationDate, 
        CreditTrancheLossConfigConstSP tranche,  
        YieldCurveConstSP              discount,
        bool                           tweakAll,
        const int                      numberOfNeighbours) const = 0;

private:
    static void load(CClassSP& clazz);
};

DECLARE(IConvolutionModel);

DRLIB_END_NAMESPACE

#endif
