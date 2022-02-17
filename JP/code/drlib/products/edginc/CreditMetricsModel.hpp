//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditMetricsModel.hpp
//
//   Description : Simple closed form model for Credit Metrics
//
//   Author      : Antoine Gregoire
//
//   Date        : April 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDIT_METRICS_HPP
#define EDR_CREDIT_METRICS_HPP

#include "edginc/CreditAsset.hpp"
#include "edginc/ConvolutionEngine.hpp"
#include "edginc/ConvolutionModel.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CreditIndexMap.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(NToDefaultLossConfig);

/** CreditMetricsModel has been rewritten to automatically convert itself into
 *  a ConvolutionEngine (so for all intents and purposes it will look like a 
 *  ConvolutionEngine, even in instanceOf calls). Previously the delegation has been 
 *  done by implementing RunMulti and delegating that call to a ConvultionEngine. 
 *  However this did not only delegate the pricing, but also the flow of control, 
 *  which is a bad thing (and it broke MultiModel).
 *  Because of the the automatic conversion CreditMetricsModel is no longer
 *  implementing IModel, but CObject.
 *  The conversion is done in the "convert" method.
 *  Linus Thand, July 2006
*/

class IRGridPointCache;
class CreditMetricsLossCalculatorBase;
class IFixedTrancheLossCalculator;
class ITrancheLossCalculator;
typedef smartPtr<IRGridPointCache> IRGridPointCacheSP;

/** Credit Metrics Model, it also drives how the quanto adjustment is done
    for cds par spreads */

class PRODUCTS_DLL CreditMetricsModel: public CObject,
                                       public virtual IConvolutionModel,
                                       virtual public IHasForwardRatePricer,
                                       virtual public ITypeConvert {
public:
    //// Automatically convert to an instance of ConvolutionEngine.
    virtual void convert(IObjectSP& object, CClassConstSP requiredType) const;

    friend class CreditMetricsModelHelper;
    friend class CDO;
    static CClassConstSP const TYPE;
    
    //values for quickGreeks
    static const string QG_DEFAULT;
    static const string QG_YES;
    static const string QG_NO;

    virtual ~CreditMetricsModel();

    /** Validation */
    virtual void validatePop2Object();                    

    /** Invoked by the containing model before the instrument data is fetched
        ie before CInstrument::GetMarket is invoked. Uses for indexMap data */
    virtual void preFetchMarketData(IModel*            model,
                                    MarketDataConstSP market);

    /** Invoked by the containing model after the instrument data is fetched
        ie after CInstrument::GetMarket is invoked. This implementation 
        sorts out pvToSpot and cfCutOffDate as well as index data */
    virtual void postFetchMarketData(IModel*            model,
                                     MarketDataConstSP market);

    /** Generate the timeline on which the effective curve will be calculated */
    virtual DateTimeArraySP generateTimeline(
        const DateTime& today,
        const DateTime& lastObservationDate) const;

    /** Static method to generate the timeline on which the effective curve 
        will be calculated */
    static DateTimeArraySP generateTimelineStatic(
        const DateTime&  today,
        const DateTime&  lastObservationDate,
        MaturityPeriodSP effSpreadFreqObj);

    virtual void createLossCalculators(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,               /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,  /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,          /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const; //(O)

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP createEffCurveGenerator(
        NToDefaultLossConfigConstSP ntdLossCfg,
        CounterPartyCreditConstSP   cpty,
        const bool                  recoverNotional) const;

    /** Indicates whether this model supports stochastic recovery rates
        AND there are any engine parameters for any names specifying so.
        By default returns false - derived types may override this method */
    virtual const bool hasStochasticRecoveries(
        CreditTrancheLossConfigConstSP tranche) const;

    /** Recovery of notional from the top of the portfolio requires
        an additional call to the convolution.
        It only has effect if the upper strike is greater than the
        sum of (name notional * name recovery).
        Therefore we can avoid making this additional call with the
        models assistance if the product displays the correct
        characteristics. Typically this means the model must be of
        a fixed recovery type */
    virtual const bool modelRecoveredNotional(
        CreditTrancheLossConfigConstSP tranche) const; /* (I) */

    /** Returns how the effectiveCurve should be interpolated. */
    virtual const string& getLossInterpolation() const;

    /** Returns all the points on the skew surface (of given name) to which this
        product is sensitive. This implementation returns null */
    virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
        OutputNameConstSP              outputName,
        const DateTime&                lastObservationDate, 
        CreditTrancheLossConfigConstSP tranche,  
        YieldCurveConstSP              discount,
        bool                           tweakAll, 
        const int                      numberOfNeighbours) const;

    /** Simple constructor */
    CreditMetricsModel(const CClassConstSP& clazz);

    /** Create a MarketDataFetcher which will be used for retrieving
        market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    /** Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model.
     * Returns riskMappingIrrelevant currently but you can maybe imagine RM
     * being useful here?  Note that this implements both
     * IModel::wantsRiskMapping() and IConvolutionModel::wantsRiskMapping(),
     * which is tricksy but does the right thing. */
    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    //this functionality is bypassed for the time being
    // ConvolutionCacheSP getConvolutionCache() const;

    //------------------------------
    // IHasForwardRatePricer methods
    //------------------------------

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;

protected:
    /** Returns true if the 'fast' convolution should be used */
//     bool useFastConvolution(const ConvolutionProduct* product) const;
    bool useFastConvolution(CreditTrancheLossConfigConstSP tranche) const;

    /** Returns an IFixedTrancheLossCalculator which is capable of
        returning expected tranche losses along the specified timeline
        (with the strikes as specified in product). See ConvolutionModel
        for further details. The implementation here uses
        createLossCalculator() fed into 
        ITrancheLossCalculator::createFixedTrancheLossCalculator. The
        conditionalLossCalc parameter is not set */
    virtual IFixedTrancheLossCalculator* createFixedLossCalculator(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const; // (O)

    /** Analogous to the loss calculator creation, but this allows us
        to model recovered notional, by pricing an inverted 1-k2, 1-k1 tranche
        instead. */
    virtual IFixedTrancheLossCalculator* createFixedRecoveredNotionalCalculator(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const; // O

    /** Returns an ITrancheLossCalculator which is capable of
        returning expected tranche losses along the specified timeline
        using the methodology appropriate to this model. The computeCondCurve
        indicates whether losses conditional on the counterparty
        surviving are required. The implementation here just calls 
        createLossCalculatorBase. The indirection/complexity is due to the
        rather dubious inheritance hierarchy employed. */
    virtual ITrancheLossCalculator* createLossCalculator(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** Computes the loss unit (bin size for loss distribution discretization)
        as dictated by this model (NB model can override it) */
    double calculateLossUnit(CreditTrancheLossConfigConstSP tranche) const;

    /** Returns an ITrancheLossCalculator which is capable of
        returning expected tranche losses along the specified timeline
        using the methodology appropriate to this model. The computeCondCurve
        indicates whether losses conditional on the counterparty
        surviving are required. The implementation here just calls 
        createLossCalculatorBase. The indirection/complexity is due to the
        rather dubious inheritance hierarchy employed.  */
    virtual ITrancheLossCalculator* createRecoveredNotionalCalculator(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** Same as above createLossCalculator() method but returns one which
        is derived from CreditMetricsLossCalculatorBase. This method won't
        make sense for certain derived types - or rather its meaning changes.
        The implementation here creates an appropriate CreditMetrics loss
        calculator */
    virtual CreditMetricsLossCalculatorBase* createLossCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** Same as above createLossCalculator() method but returns one which
        is derived from CreditMetricsLossCalculatorBase. This method won't
        make sense for certain derived types - or rather its meaning changes.
        The implementation here creates an appropriate CreditMetrics loss
        calculator */
    virtual CreditMetricsLossCalculatorBase* createRecoveredNotionalCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** A utility method (as inheritance hiearchy is questionable).
        Returns an ITrancheLossCalculator which is capable of returning
        expected tranche losses along the specified timeline using CCM
        algorithm.  The computeCondCurve indicates whether losses
        conditional on the counterparty surviving are required. */
    CreditMetricsLossCalculatorBase* createCCMLossCalculator(
        const DateTimeArray&           timeline,     /* (I) */
        CreditTrancheLossConfigConstSP tranche,      /* (I) */
        CounterPartyCreditConstSP      cpty) const;  /* (I) */
        
    /** A utility method (as inheritance hiearchy is questionable).
        Returns an ITrancheLossCalculator which is capable of returning
        expected tranche losses along the specified timeline using CCM
        algorithm.  The computeCondCurve indicates whether losses
        conditional on the counterparty surviving are required. */
    CreditMetricsLossCalculatorBase* createCCMRecoveredNotionalCalculator(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */


private:
    CreditMetricsModel(const CreditMetricsModel &rhs);
    CreditMetricsModel& operator=(const CreditMetricsModel& rhs);

    // FIELDS
    bool                   pvToSpot;             /* TRUE : return price with pv to spot date on           */
                                                 /*        instrument discount curve                      */
                                                 /* FALSE : return price with pv to today                 */
    DateTime               cfCutOffDate;         /* ignore CF before and on this date                     */
    int                    maxNbSlice;           /* max nb of slices for loss distribution                */
    int                    gcdDivisor;           /* nb by which to divide the GCD of notionals            */
    double                 lossUnitOverride;     /* Loss unit override                                    */
    string                 lossInterpolation;    /* L=linear loss interp, F=flat fwd                      */
    string                 effSpreadFreq;        /* date template frequency for effective spread curve    */
    string                 creditChargeViewType; /* creditCharge view type                                */
    int                    fastConvThreshold;    /* if nbName >= threshold, use fast convolution          */
    string                 quickGreeks;          /* DEFAULT, YES, NO : are quick greeks allowed           */
    ICreditIndexMapWrapper indexMap;             /* allows single names to become index basis adjusted    */
                                                 /* via its corresponding index definition                */
    string                 indexMapTreatment;    /* controls how the index map may be used                */
    CBoolSP                calculateIndexBasis;  /* TRUE  : calculate the basis from the index definition */
                                                 /* FALSE : retrieve named basis from the cache           */

    string calibrationStyle;    // eg CMS - for quanto  - Now deprecated
    string calibrationMaturity; // eg 10Y - for quanto  - Now deprecated
    IForwardRatePricerSP    forwardRateModel;    /* used for determining fee leg cashflows                */

    // transient
    MaturityPeriodSP effSpreadFreqObj; /* date template frequency for effective spread curve */
};

DRLIB_END_NAMESPACE
#endif

