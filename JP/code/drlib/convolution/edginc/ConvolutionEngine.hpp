//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Family of CModels that use convolution to get an 
//                 'effective curve'
//
//   Date        : 18th Nov 2005
//
//   Author      : Mark Robson
//
//----------------------------------------------------------------------------

#ifndef QR_CONVOLUTIONENGINE_HPP
#define QR_CONVOLUTIONENGINE_HPP

#include "edginc/Model.hpp" 
#include "edginc/BetaSkewPointwiseTweak.hpp"
#include "edginc/IRVegaPointwise.hpp"
#include "edginc/QuantoCDSParSpreads.hpp"
#include "edginc/IEffectiveCurveLossModelConfig.hpp"
#include "edginc/IEffectiveCurveLossGen.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/CreditTrancheLossConfig.hpp" // temporarily here, remove
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/ConvolutionModel.hpp" // temporarily here, remove
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IGeneralisedConvolutionProduct);
FORWARD_DECLARE(ConvolutionEngine);
FORWARD_DECLARE(QuantoCDSAlgorithm);
FORWARD_DECLARE(IRGridPointCache);
FORWARD_DECLARE(IRGridPointCache);
FORWARD_DECLARE(ICreditLossConfig);
FORWARD_DECLARE(CreditTrancheLossConfig);
FORWARD_DECLARE_REF_COUNT(ICreditLossGen);
FORWARD_DECLARE(IFixedTrancheLossCalculator);

/** An implementation of CModel that uses 'convolution' [of probabilities] to
    price credit dependent instruments. The approach is to use composition
    rather than inheritance to support the different ways the convolution
    can be driven. In particular, an IConvolutionModel needs to be supplied
    which handles the modelling side of things. This class coordinates the
    activities beteen the IConvolutionModel and the instrument. */
class CONVOLUTION_DLL ConvolutionEngine:
    public CModel,
    public virtual QuantoCDSParSpreads::IAlgorithmBuilder,
    public virtual IRVegaPointwise::ISensitivePoints,
    public virtual BetaSkewPointwiseTweak::ISensitivePoints,
    public virtual IEffectiveCurveLossModelConfig,
    public virtual IHasForwardRatePricer
{
public:
    static CClassConstSP const TYPE;
    virtual ~ConvolutionEngine();

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;

    /** overridden to [reference] copy QuantoCDSAlgorithmSP */
    virtual IObject* clone() const;

    //// called immediately after object constructed
    void validatePop2Object();

    /** interface that the instrument must implement */
    class CONVOLUTION_DLL IIntoProduct: virtual public CModel::IModelIntoProduct{
    public:
        static CClassConstSP const TYPE;

        virtual ~IIntoProduct();

        /** Creates an instance of an IGeneralisedConvolutionProduct */
        virtual IGeneralisedConvolutionProduct* createProduct(
            ConvolutionEngineConstSP model) const = 0;
    private:
        static void load(CClassSP& clazz);
    };
    
    /** Overrides default Model implementation to retrieve today, then
        calls ConvolutionModel::preFetchMarketData before invoking
        parent method, and then finally invoking 
        ConvolutionModel::postFetchMarketData */
    virtual void getInstrumentsAndModelMarket(MarketDataConstSP       market,
                                              IInstrumentCollectionSP insts);

    /** Creates an IEffectiveCurveLossGen for the supplied ICreditLossConfig
        object. The ICreditLossConfig object must implement the IntoLossGen
        interface */
    virtual ICreditLossGenSP lossGenerator(
        ICreditLossConfigConstSP lossConfig,
        IModelConfigMapperConstSP mapper) const;

    /** Creates an IEffectiveCurveGen for the supplied ICreditLossConfig object. 
        The ICreditLossConfig object must implement the IntoEffCurveGen interface */
    // NB: "mapper" is currently NOT USED
    virtual ICreditLossGenSP effCurveGenerator(
        ICreditLossConfigConstSP  lossConfig,
        CounterPartyCreditConstSP cpty,
        const bool                recoverNotional,
        IModelConfigMapperConstSP mapper) const;

    /** Creates an IEffectiveCurveLossGen for a 'tranche'*/
    virtual IEffectiveCurveLossGenSP createLossGenerator(
        CreditTrancheLossConfigConstSP trancheLossCfg) const;

	virtual IEffectiveCurveLossGenSP createLossGenerator(
		FlatCDO2LossConfigConstSP tsLossCfg) const;

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP createEffCurveGenerator(
        NToDefaultLossConfigConstSP ntdLossCfg,
        CounterPartyCreditConstSP   cpty,
        const bool                  recoverNotional) const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument* instrument, 
                       CControl*    control, 
                       CResults*    results);

    /** Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this process.  Delegates to
     * IConvolutionModel::wantsRiskMapping(). See IModel::wantsRiskMapping(). */
    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** Returns all the points on the skew surface to which this
        model for the supplied instrument is sensitive.
        A null return value is ok - it is interpreted as the greek being 
        NotApplicable */
    virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst,
        bool               tweakAll,
        const int          numberOfNeighbours) const;

    /** Returns conservative/aggressive etc. To do: either enums or sort out
        where strings defined */
    const string& getCreditChargeViewType() const;

    /** Kapital backward compatibility mode - to be removed eventually.
        Returns whether we pv to cfCutOffDate() */
    bool pvToSpot() const;
    
    /** Kapital backward compatibility mode - to be removed eventually.
        Ignore cashflows as well as protection etc before this date. This
        method may return an empty DateTime if no cfCutOffDate has been
        specified. */
    const DateTime& cfCutOffDate() const;
    
    /** Constructor - mainly for backward compatibility. Allow models
        to create this object and then run the pricing through this
        class (For use by old style specification of model) */
    ConvolutionEngine(IConvolutionModelSP  convolutionModel,
                      const DateTime&      today,
                      const string&        creditChargeViewType,
                      bool                 pvToSpot,  
                      const DateTime&      cfCutOffDate,
                      const string& calibrationStyle,
                      const string& calibrationMaturity);

    /** Overridden and redirected to ConvolutionModel */
    MarketDataFetcherSP createMDF() const;

    /** when to stop tweaking - to do change infrastructure to route through
        model rather than instrument */
    DateTime endDate(const CInstrument*  inst,
                     const Sensitivity*  sensControl) const;

    /** overrides Models implementation to shift cfCutOffDate */
    virtual bool sensShift(Theta* shift);

    /** Returns an instance of IAlgorithm. Typically the quantoCDSParSpreads
        parameter would be ignored but is there in case you want to switch
        the algorithm dependent upon some property of the quanto'd curve */
    virtual QuantoCDSParSpreads::IAlgorithmSP cdsQuantoAlgorithm(
        const QuantoCDSParSpreads* quantoCDSParSpreads) const;

    /** sets debug state to specified value. If true, then any calls to
        cdsQuantoAlgorithm will return an object that will cache debug data.
        This can be retrieved via getDebugInfo() */
    virtual void selectDebugState(bool switchOn);

    /** Uses quanto calibration parameters to identify relevant points */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;

    /** Returns debug info - may be null eg if the algorithm has not been 
        used */
    virtual IObjectSP getDebugInfo() const;

    // Return convolution model
    IConvolutionModelSP ConvolutionModel() { return convolutionModel; }

private:
    friend class PortfolioSpreadCurvesAddin; // in ConvolutionEngine.cpp
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    ConvolutionEngine();
    IGeneralisedConvolutionProduct* createProduct(const CInstrument* instrument) const;

    void createQuantoAlgorithm(bool debugOn) const;

    static bool isYCSens(const Sensitivity* sensitivity);
    static bool isIRVegaPointwiseSens(const Sensitivity* sensitivity);

    DateTime rhoEndDate(const Sensitivity* sensitivity,
                        const DateTime&    irVegaEndDate) const;

    DateTime endDateForYCTweaks(
        const IGeneralisedConvolutionProduct* product,
        IEffectiveCurveLossGenSP              effCurveLossGen,
        const Sensitivity*                    sensitivity,
        const DateTime&                       creditEndDate) const;

    /// fields //////////////////////////////////////////
    IConvolutionModelSP convolutionModel;
    string              creditChargeViewType;
    bool                doPVToSpot; // transient
    DateTime            theCfCutOffDate; // transient
    string              calibrationStyle;  // eg CMS - for quanto
    string              calibrationMaturity; // eg 10Y - for quanto
    mutable QuantoCDSAlgorithmSP quantoAlgorithm; // transient, not registered $unregistered
    mutable IRGridPointCacheSP   irGridPtsCache;  // transient, not registered $unregistered

    //// this class essentially just wraps the existing IConvolutionModel
    //// and holds the extra data needed to implement the IEffectiveCurveLossGen
    //// interface. Some methods on IConvolutionModel will need changing since
    //// we are no longer using a ConvolutionProduct but a FlatCDO2LossConfig.
    //// The data that the IConvolutionModel requires should be available through
    //// the FlatCDO2LossConfig
    class CONVOLUTION_DLL FlatCDO2LossCurveGen : virtual public IEffectiveCurveLossGen {

    private:
        IConvolutionModelSP convolutionModel;
        ConvolutionEngineConstSP engineForInnerLossConfigs;
        FlatCDO2LossConfigConstSP ts;

    public:
        FlatCDO2LossCurveGen(IConvolutionModelSP convolutionModel,
                            ConvolutionEngineConstSP engineForInnerLossConfigs,
                            FlatCDO2LossConfigConstSP ts) :
            convolutionModel(convolutionModel), 
            engineForInnerLossConfigs(engineForInnerLossConfigs),
            ts(ts)
        {}

        /** Generate the timeline on which the effective curve will be 
            calculated. No customization implemented here, just forward 
            to the convolution model */
        virtual DateTimeArraySP generateTimeline(
            const DateTime& lastObservationDate) const 
        {
            return convolutionModel->generateTimeline(ts->getToday(),
                                                      lastObservationDate);
        }

        /** Returns a set of IFixedFlatCDO2LossCalculator which are capable
            of returning expected TS losses along the specified
            timeline. The cpty indicates whether losses conditional on the
            counterparty surviving are required (null if not). The control
            and results parameters are there for any output requests that
            should be calculated.  */
        virtual void createLossCalculators(
            const DateTimeArray&                timeline,           // (I)
            CounterPartyCreditConstSP           cpty,               // (I)
            const DateTime&                     maturity,           // (I)
            Control*                            control,            // (I)
            Results*                            results,            // (I)
            bool                                recoverNotional,    // (I)
            YieldCurveConstSP                   discount,           // (I)
            IFixedTrancheLossCalculatorConstSP& lossCalculator,               // (O)
            IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,  // (O)
            IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,          // (O)
            IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const // (O)
        {
            return convolutionModel->createLossCalculatorsFIX(
                timeline,
                cpty,
                maturity,
                control,
                results,
                recoverNotional,
                discount,
                lossCalculator,
                recoveredNotionalCalculator,
                conditionalLossCalc,
                conditionalRecNtnlCalc,
				ts); 
        }
        

        /** Returns all the points on the skew surface to which this
            model for the supplied instrument is sensitive.
            A null return value is ok - it is interpreted as the greek being 
            NotApplicable
            The lastObservationDate argument must be instrument's last observation 
            date for both the fee and contingent legs */
        virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
            OutputNameConstSP outputName,
            const DateTime&   lastObservationDate,
            YieldCurveConstSP discount,
            bool              tweakAll,
            const int         numberOfNeighbours) const
        {   
			throw ModelException("Not implemented");
        }

        /** When to stop tweaking the loss config(s) associated to this loss gen */
        virtual DateTime maxLossConfigsTweakableDate(const DateTime& creditEndDate) const {
            DateTime lastExpiry;
            DateTime maxDate(ts->getToday());
            int numNames = ts->numInnerLossConfigs();
            for (int i=0; i < numNames; ++i) {
                SingleCreditAssetConstSP myAsset(ts->nameAsset(i));
                lastExpiry = myAsset->getParSpreadCurve()->stopTweaking(creditEndDate);
                if (lastExpiry.isGreater(maxDate)) {
                    maxDate = lastExpiry;
                }            
            }
            return maxDate;
        }

        /** Recovery of notional from the top of the portfolio requires
            an additional call to the convolution.
            It only has effect if the upper strike is greater than the
            sum of (name notional * name recovery).
            Therefore we can avoid making this additional call with the
            models assistance if the product displays the correct
            characteristics. Typically this means the model must be of
            a fixed recovery type */
        virtual const bool modelRecoveredNotional() const  {
			throw ModelException("Not implemented");
        }

        /** Returns how the effectiveCurve should be interpolated. This should
            really disappear when we return a DefaultRate object say in the 
            calculateEffectiveCurve() method. */
        virtual const string& getLossInterpolation() const {
            return convolutionModel->getLossInterpolation();
        }
    };
    DECLARE_REF_COUNT(FlatCDO2LossCurveGen);

};

DRLIB_END_NAMESPACE

#endif
