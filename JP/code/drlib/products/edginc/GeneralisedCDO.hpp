//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Description : General class that takes an ICreditLossConfig plus a fee
//               and contingent legs
//
// Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_GENERALISEDCDO_HPP
#define QLIB_GENERALISEDCDO_HPP

#include "edginc/CDOPortfolio.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/ConvolutionEngine.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/ICreditLegConvention.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/IProtectionProvider.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/MCPrices.hpp"
#include "edginc/ConditionalLossModel.hpp"
#include "edginc/MCProductClient.hpp"
#include "edginc/GenericSimpleIR.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/IHasMaturityDate.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/KComponent.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_WRAPPER(ICreditLossConfig);
FORWARD_DECLARE_WRAPPER(IForwardRatePricer);
FORWARD_DECLARE(CreditMetricsBaseCorrelation);
FORWARD_DECLARE(FtDProduct);

class CDOPriceDetails;

typedef smartPtr<GeneralisedCDO> GeneralisedCDOSP;
typedef smartConstPtr<GeneralisedCDO> GeneralisedCDOConstSP;

/** Credit Tranche Instrument. */
class PRODUCTS_DLL GeneralisedCDO:
    public KComponent,						// for tree pricing: use KComponent instead of CInstrument
    public virtual LastSensDate,
    public virtual Theta::Shift,
    public virtual ConvolutionEngine::IIntoProduct,
	public virtual FDModel::IIntoProduct,	// for tree pricing
	public virtual ConditionalLossModel::IIntoProduct,
	public virtual MonteCarlo::IIntoProduct,
    public virtual IBadDayAdjuster,
    public virtual IRebateCalculator,
    public virtual IProtectionProvider,
	public virtual IHasMaturityDate
{

protected:

	class MyConvolutionProduct: public ConvolutionProduct {

	public:
        /** Contains data we want to store. */
        class MyOutput : public ConvolutionProduct::Output {
        public:
            /** Constructor */
            MyOutput();

            /** Returns the fair value (aka mark to market) given this
                Output object */
            virtual double price() const;

            /** Returns fee leg price for CCC view (positive cashflows only).
                Requires fLegDebugUnitPrice to be populated */
            double feeLegPriceForCCC(int longOrShort);

            /** Returns fee leg price for CCC view (positive cashflows only).
                Requires cLegDebugUnitPrice to be populated/ */
            double contingentLegPriceForCCC(int longOrShort);

            /** Fudge for Kapital which wants fair value on [yield curve] spot
                date */
            virtual void scaleByForwardFactor(double fwdFact);

            /** Builds CDOPriceDetails object for fee leg */
            CDOPriceDetails* makeFeeLegPriceDetails(int longOrShort) const;

            /** Builds CDOPriceDetails object for contingent leg */
            CDOPriceDetails* makeContingentLegPriceDetails(int longOrShort) const;

            //-------
            // Fields
            //-------
            double          cLegPV;             // fair value of contingent leg
            double          fLegPV;             // fair value of fee leg
            // the remaining fields are for output requests
            double          riskyDurationTotal; // notional weighted risky duration
            double          riskyNotionalsMean; // mean value of notional per period
            double          risklessCFPV;       // fair value of riskless payments
            DoubleArraySP   cLegDebugUnitPrice;
            DoubleArraySP   cLegDebugUnitHistPrice;
            DoubleArraySP   fLegDebugUnitPrice;
            DoubleArraySP   fLegDebugUnitHistPrice;

        private:
            //// just returns sum of unitPrice whose sign is the same as longOrShort
            static double priceCreditChargeView(const DoubleArray& unitPrice,
                                                int                longOrShort);
        };
        typedef smartPtr<MyOutput> MyOutputSP;

        /** Construtor */
        MyConvolutionProduct(ConvolutionEngineConstSP convolutionEngine,
                             const GeneralisedCDO*    cdo);

        /** Returns the object that defines the losses for this product */
        virtual ICreditLossConfigConstSP getLossConfig() const;

        /** Returns the last observation date for credit related data */
        virtual DateTime lastObservationDate() const;

        /** Returns the last 'scheduled' pay date (ie ignoring credit
            events).  Introduced for the output request IND_CDS_PAR_SPREAD
            - not clear if this is right date to use */
        virtual DateTime lastPayDate() const;

        /** Returns the max of maxDate and the last date from when a yield
            curve is used (eg for discounting or rate estimation) by
            payoff method.  This method does not have to worry about what
            dates are used when the clean spread curves are built. Used to
            control when to stop tweaking */
        virtual DateTime lastYCSensDate(
            const DateTime& maxDate) const;

        /** Returns number of names contained within instrument */
        virtual int numNames() const;

        /** Returns the representation of the underlying name corresponding to the
            specified index which must lie in range [0, numNames()-1] */
        virtual SingleCreditAssetConstSP nameAsset(int index) const;

        /** Returns the max of valueDate and protection start date
            for specified name (if any). The specified index must lie in
            range [0, numNames()-1]*/
        virtual const DateTime& nameProtectionStartDate(
            int index) const;

        /** Returns the min of last date parameter and protection end date
            for specified name (if any). The specified index must lie in
            range [0, numNames()-1]. This if for name 'cutoff' */
        virtual const DateTime& nameProtectionEndDate(
            int index,
            const DateTime& lastDate) const;

        /** Returns the recovery for the name corresponding to the
            specified index which must lie in range [0, numNames()-1]. Note that
            the recovery can  be overridden at the trade level hence this method,
            which reflects any trade level overrides, should be used rather than
            the recovery off the par spread curve */
        virtual double nameRecovery(int index) const;

        /** Returns true if the name has defaulted where the name corresponds to the
            specified index which must lie in range [0, numNames()-1] */
        virtual bool nameDefaulted(int index) const;

        /** Returns the notional for the name corresponding to the
            specified index which must lie in range [0, numNames()-1] */
        virtual double nameNotional(int index) const;

        /** Returns the loss given default for the name corresponding to the
            specified index which must lie in range [0, numNames()-1] */
        virtual double nameLossGivenDefault(int index) const;

        /** Returns the 'beta' for the name corresponding to the specified
            index which must lie in range [0, numNames()-1] */
        virtual double nameBeta(int index) const;

        /** Returns the correlation between the loss and correlation market */
        virtual double lossDecretionBeta() const;

        /** Returns the ranges for Loss Given Default for the name
            corresponding to the specified index which must lie in range
            [0, numNames()-1]. In the particular, the LGD of a name is given by
            nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
        virtual void nameLGDRanges(
            int     index,       /* (I) */
            double& lgdNotional, /* (O) */
            double& lgdFloor,    /* (O) */
            double& lgdCap)      /* (O) */ const;

        /** Returns a CounterPartyCredit object containing data required for
            pricing CCC. Note may return null if no counter party information */
        virtual CounterPartyCreditConstSP getCounterParty() const;

        /** Returns the total notional of the portfolio */
        virtual double portfolioNotional() const;

        /** Returns the total long notional of the portfolio */
        virtual double portfolioLongNotional() const;

        /** Returns the total short notional of the portfolio */
        virtual double portfolioShortNotional() const;

        /** Returns the sum of the past losses of the portfolio which eat
            into the total notional returned by portfolioNotional() */
        virtual double portfolioLoss() const;

        /** Returns the sum of the long past losses of the portfolio which eat
            into the total notional returned by portfolioNotional() */
        virtual double portfolioLongLoss() const;

        /** Returns the sum of the short past losses of the portfolio which eat
            into the total notional returned by portfolioNotional() */
        virtual double portfolioShortLoss() const;

        /** Returns whether the instument requires notional to be
            recovered from the top of the portfolio NB may be later
            overridden by the model if it is deemed to be unnecessary */
        virtual bool recoverNotional() const;

        /** Returns the sum of the past recovered notional of the
        * portfolio which eat into the total notional returned by
        * portfolioNotional() */
        virtual double portfolioRecoveredNotional() const;

        /** Returns the sum of the past long recovered notional of the
        * portfolio which eat into the total notional returned by
        * portfolioNotional() */
        virtual double portfolioLongRecoveredNotional() const;

        /** Returns the sum of the past short recovered notional of the
        * portfolio which eat into the total notional returned by
        * portfolioNotional() */
        virtual double portfolioShortRecoveredNotional() const;

        /** Returns the strikes of the tranche. These are relative to the values
            returned by portfolioRanges() */
        virtual void trancheStrikes(double& lowStrike,
                                    double& highStrike) const;

        /** map loss interpolation string to CCMPriceUtil enum */
        static CCMPriceUtil::ExpLossType mapLossInterpolationType(
            const string& lossInterpolation);

        /** Computes payoff for instrument given effectiveCurve */
        virtual OutputSP payoff(
            EffectiveCurveSP     cLegEffectiveCurve,
            EffectiveCurveSP     fLegEffectiveCurve,
            double               forwardFactor,
            IForwardRatePricerSP model) const;

        /** Compute the credit charge using supplied results and credit charge
            view type. To do: sort out rep of creditChargeViewType */
        virtual double computeCreditCharge(
            const EffectiveCurveSP cptyCurve,
            const string&          creditChargeViewType,
            OutputSP               output,
            OutputSP               cccOutput) const;

        /** Populate CResults object with price and any output requests etc. Note
            cccOutput may be null if CCC not being computed */
        virtual void storeExtraResults(CResults* results,
                                       CControl* control,
                                       OutputSP  output,
                                       OutputSP  cccOutput, // may be 0
                                       IForwardRatePricerSP model,
                                       const IConvolutionModel* convolutionModel) const;

    private:
        void computeAndStoreDurWeightAvg(const CreditMetricsBaseCorrelation* bc,
                                         const CreditTrancheLossConfig* tranche,
                                         CResults* results,
                                         OutputRequest* request) const;

        //-------
        // Fields
        //-------
        smartConstPtr<GeneralisedCDO> cdo; // the associated instrument
        CashFlowArraySP    cLegPastNotionalReductions; // uses 1-R for losses
        BoolArraySP        cLegPayPastNotionalReductions; /* Whether to pay for the
                                                           * notional reductions (if
                                                           * the default happens during
                                                           * a protecion period) or not */
        FeeLegReductionPerDefaultArraySP fLegPastReductions;
        double totalLongNotional;
        double totalShortNotional;
        double totalPastLongLoss;
        double totalPastShortLoss;
        double totalPastRecoveredNotional;
        double totalPastLongRecoveredNotional;
        double totalPastShortRecoveredNotional;
        double lowStrike;  // Tranche low strike (absolute value, not a percentage)
        double highStrike; // Tranche high strike (absolute value, not a percentage)
        ICreditLossConfig* lossCfg;  // JLH - temporarily
        CDOPortfolioConstSP portf; // JLH - temporarily, keep a (cast) of the CDOPortfolio
    };

    // ConvolutionProduct needs access
    friend class MyConvolutionProduct;

	// SpreadLossTree product needs access
	friend class CDOTree;

	// N-to-default product needs access
	friend class FtDProduct;

public:

    static CClassConstSP const TYPE;

	virtual ~GeneralisedCDO();

    virtual void validatePop2Object();

    /** called after market data has been retrieved */
    virtual void Validate();

    /** Returns the recover notional flag */
    bool getRecoverNotional() const;

	//----------------------------------------
	// KComponent / CInstrument methods
	//----------------------------------------
	virtual DateTime getValueDate() const;


    /** Retrieve market data from cache for all our components */
	// Note that in the KComponent framework we should not override KComponent::GetMarket
	// Instead should use setup()
	// At the moment we only use the KComponent framework when we are pricing with the SpreadLossTree Model.
	// Could use this more generally once we are happy that the KComponent stuff works.
    virtual void GetMarket(const IModel* model,
                           const CMarketDataSP market);

	/* For LastSensDate interface */
	virtual DateTime endDate(const Sensitivity* sensControl) const;

	/** Setup the component */
	// The derived function must also call KComponent::setup() and the setup() of its underlyings
	// In the KComponent framework this method does everything usually found in GetMarket except
	// for populating market wrappers (this is done in KComponent::GetMarket)
	virtual void setup(const IModel* model, const MarketData* market);

    //----------------------
    // Theta::IShift methods
    //----------------------
    virtual bool sensShift(Theta* shift);

    //------------------------------------------
    //  IBadDayAdjuster methods
    //------------------------------------------
    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;

    //------------------------------------------
    //  IRebateCalculator method
    //------------------------------------------
    /** Returns the rebate amount, i.e., the difference in fee payments under
     * different loss assumptions.
     * Note the amount is not be pv'd in the sense that, if the fees should
     * have been X dollars higher using the "real losses" vs using the "assumed
     * losses", the rebate will be X dollars: the fact that the rebate may be
     * paid on a different date does not impact the rebate computation.
     * Actually this method is not at all aware of when the actual rebate
     * payment will take place. */
    virtual double computeRebate(
        const CashFlowArrayConstSP assumedTrancheReductions,
        const CashFlowArrayConstSP realTrancheReductions,
        const double initialNotional,
        IForwardRatePricerSP model) const;

    //------------------------------------------
    //  IProtectionProvider method
    //------------------------------------------
    /* Checks if the input date is covered for protection */
    virtual bool isDateCoveredForProtection(const DateTime& date) const;

	// Supporting MC model

	/** Forward declaring. Product class when the CDO instrument implements the MC model  */
	FORWARD_DECLARE(GeneralisedCDOMC);

	/** GeneralisedCDOMC access the CDO members */
	friend class GeneralisedCDOMC;

	/** Method of the interface, IMCIntoProduct. Represents the IntoProduct representation of the CDO instrument
		that would be priced using MC
		Caution: we modify the MCPathConfig that is a member of the MonteCarlo. In that sense the const character of the model is
		violated.
		*/
	// can return either the usual MC product
	// or a CondLossMCProduct, depending on a field in *this
	IMCProduct* createProduct(const MonteCarlo* model) const;

    /** Returns the object that defines the losses for this product */
    virtual ICreditLossConfigConstSP getLossConfig() const;

protected:
	GeneralisedCDO(CClassConstSP clazz = TYPE);

    GeneralisedCDO(
        CClassConstSP                clazz,
        bool                         isLong,
        smartPtr<CounterPartyCredit> cptyInfo,
        ICreditContingentLegSP       cLeg,
        ICreditFeeLegSP              fLeg,
        CDOPortfolioSP               portfolio,
        CBoolSP                      recoverNotional,
        const string&                yieldCurveName,
        CIntSP                       triggerDelay = CIntSP(),
        CIntSP                       defaultToCalculationDelay = CIntSP(),
        DateTime                     lastTriggerDate = DateTime(),
        CDoubleSP                    temporaryLossAmount = CDoubleSP(),
        const string&                settlementBDC = "N",
        HolidayWrapper               settlementHols = HolidayWrapper());

    /** Returns the last 'scheduled' pay date (ie ignoring credit
        events).  Introduced for the output request IND_CDS_PAR_SPREAD
        - not clear if this is right date to use */
    DateTime maturityDate() const;

    /** Return the maximum observation end date of the legs */
    DateTime lastObservationDate() const;

    /** Returns the max of maxDate and the last date from when a yield
        curve is used (eg for discounting or rate estimation).
        This method does not have to worry about what dates are used
        when the clean spread curves are built: Just used to
        control when to stop tweaking */
    DateTime lastYCSensDate(const DateTime& maxDate) const;

	//----------------------------------------
	// FDModel::IntoProduct methods
	//----------------------------------------
	virtual FDProductSP createProduct(FDModel * model) const;

    //----------------------------------------
    // ConvolutionEngine::IIntoProduct methods
    //----------------------------------------
    /** Creates an instance of an IGeneralisedConvolutionProduct */
    virtual IGeneralisedConvolutionProduct* createProduct(
        ConvolutionEngineConstSP model) const;

    //------
    //fields
    //------
    DateTime  today;
    bool      isLong;           /* Long or short the tranche deal */

	mutable bool clProduct;

	// from the CL model in the case of a CL pricing
	DateTimeArrayConstSP auxTimeLine;

	// this needs to be here so that it does not change with the tweaking
	mutable FlatCDO2LossConfigSP bucketsCL;

    smartPtr<CounterPartyCredit> cptyInfo;  /* cpty charge info */
    ICreditContingentLegSP       cLeg;      /* contingent leg */
    ICreditFeeLegSP              fLeg;      /* fee leg */
    ICreditLossConfigSP          portfolio; /* portfolio */
	// Note that discount is now in KComponent

    CBoolSP recoverNotional;  /* does the the fee leg recover notional from the top */

    // Default handling parameters
    /** Delay between default and eventDeterminationDate. It is optional.
     * If not present it will be assumed that defaults are triggered
     * on default date. Note that this is not the same as the delay
     * being 0, in which case the name's default will not be considered
     * triggered until a "credit event override" is associated to the
     * defaulted name */
    CIntSP triggerDelay;

    /** Delay between credit event and calculation date.  It is optional.
     * If not present it will be assumed that calculation date happens
     * on trigger date. Note that this is not the same as the delay
     * being 0, in which case the calculation date will keep moving
     * forward until a "credit event override" is associated to the
     * defaulted name */
    CIntSP defaultToCalculationDelay;

    /** Last date when a default occurred during the protection period
     * can be triggered.
     * CAUTION: This is an optional parameter. If not present, it will
     * be set to lastObservationDate() from within Validate - if this
     * method is not called (e.g., if building a GeneralisedCDO "by hand"
     * somewhere) it will remain empty indicating that there is no last
     * trigger date. Note it cannot be set to lastObservationDate() from
     * within validatePop2Object because lastObservationDate requires
     * "today" to be populated, which happens in getMarket. */
    DateTime lastTriggerDate;

    /** Assumption for the loss amount between a credit event and the
     * related calculation date */
    CDoubleSP temporaryLossAmount;

    /** bad day convention for the GeneralisedCDO */
    string settlementBDC;

    /** Holidays */
    HolidayWrapper settlementHols;

    // Transient
    BadDayConventionSP settlementBadDayConv;

private:
    // usual methods
    GeneralisedCDO(const GeneralisedCDO& rhs); // don't use
    GeneralisedCDO& operator=(const GeneralisedCDO& rhs); // don't use

    static IObject* defaultConstructor();
    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    CashFlowArraySP generateKnownCashflows(
        double               fLegOutstandingNotional,
        const CashFlowArray& fLegPastNotionalReductions,
        const CashFlowArray& cLegPastNotionalReductions,
        const BoolArray&     cLegPayPastNotionalReductions,
        IForwardRatePricerSP model) const;

	/** Generates timelines to be used by the MC framework */
	void
	generateMCTimelines(
		DateTimeArraySP& modifiedFeeLegObservationDates,
		DateTimeArraySP& productTimeline,
		DateTimeArraySP& discFactorTimeline,
		IntArraySP& dateToDiscFactorIndex
		) const;

	CMarketDataSP marketPointer;

public:

	DateTimeArray computeTimeLine() const;

	/* Conditional loss product */

	FORWARD_DECLARE(ConditionalLossProduct);
	friend class ConditionalLossProduct;

	FORWARD_DECLARE(ConditionalLossMCProduct);
	friend class ConditionalLossMCProduct;


	ConditionalLossModel::IProduct* createProduct(const ConditionalLossModel* model) const;

};

DECLARE(GeneralisedCDO);

//------------------------------------------------------------------------------
//  Auxiliary classes used for returning fee and contingent leg OutputRequests
//------------------------------------------------------------------------------
class CDOPriceDetails : public CObject {
public:
    static CClassConstSP const TYPE;
    CDOPriceDetails();

    /** Creates CDOPriceDetails from supplied inputs. Note that the DoubleArrays
        must be of the same length. */
    CDOPriceDetails(DoubleArraySP unitPrice,     /* price for each leg unit */
                    DoubleArraySP unitHistPrice, /*price for each leg unit due to historical default*/
                    double        priceCreditChargeView);
    
    // Return price
    double Price() const { return price; }

private:
    static IObject* defaultConstructor();

    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    //-------
    // Fields
    //-------
    double        price;                 // total leg price
    double        histPrice;             // total leg price due to historical default
    double        priceCreditChargeView; // total price of CCC view (positive cf only)
    DoubleArraySP unitPrice;             // price for each leg unit
    DoubleArraySP unitHistPrice;         // price for each leg unit due to historical default
};

typedef smartPtr<CDOPriceDetails> CDOPriceDetailsSP;


/** Leg output. Only used for returning fee and contingent leg OutputRequests*/
class CDOLegOutput : public CObject {
public:
    static CClassConstSP const TYPE;
    //// Simple constructor - just takes references to CDOPriceDetailsSP
    CDOLegOutput(CDOPriceDetailsSP price,
                 CDOPriceDetailsSP priceCond);

    // Return price
    double Price() const { return price->Price(); }

private:
    CDOLegOutput();

    static IObject* defaultConstructor();

    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    //-------
    // Fields
    //-------
    CDOPriceDetailsSP price;
    CDOPriceDetailsSP priceCond;
};

typedef smartPtr<CDOLegOutput> CDOLegOutputSP;
typedef smartPtr<const CDOLegOutput> CDOLegOutputConstSP;

// ################################################################################################################
/** Product class when the CDO instrument implements the MC model  */
class GeneralisedCDO::GeneralisedCDOMC :	public MCProductClient
{
public:
	/** Constructor  */
	GeneralisedCDOMC(
		const GeneralisedCDOConstSP& generalisedCDO,
		const SimSeriesConstSP&   simSeries,
		int numberOfIterations,
		int numberOfSubSamples,
		const DateTimeArrayConstSP& productTimeline,
		const DateTimeArrayConstSP& discFactorTimeline,
		const IntArrayConstSP& dateToDiscFactorIndex,
		const DateTimeArrayConstSP& modifiedFeeLegObservationDates
		); //needed for the MCProductClient

	/** Virtual destructor  */
	virtual ~GeneralisedCDOMC() {}

// ##  Methods of the interface, IMCProduct that MCProductClient inherits from
	/** Populates the svCollector with the SV Gens that are need by this product  */
	virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

	/** Called by the MC framework whenever the path generator is updated  */
	virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

	/** Computes the payoff of the instrument for every simulation path  */
	virtual void payoff(
		const IPathGenerator*  pathGen,
	   	IMCPrices& prices);

	/** stores extra output to be computed for the   */
	virtual void recordExtraOutput(
		Control* control,
		Results* results,
		const IMCPrices& prices) const;

	/** Date to stop tweaking for the supplied tweak  */
	virtual DateTime endDate(const Sensitivity*  sensControl) const;

private:
// ##  Member variables to be used for computation of the contingent leg. We store them since they are reused across
	/** reference to the notionals of the contingent leg period */
	DoubleArray  contNotionals;

	/** array indicating whether the ith contingent leg period is of the type payAsYouGo  */
	BoolArraySP	lpayAsYouGoArray;

	/** array indicating whether for the ith contingent leg period, the payment is after a delay period. Only applicable for payAsYouGo  */
	IntArraySP	lnumDelayDaysArray;

	/** array indicating the start date of the ith contingent leg period  */
	DateTimeArraySP	lstartDateArray;

	/** array indicating the end date of the ith contingent leg period  */
	DateTimeArraySP	lendDateArray;

	/** array indicating the pay date of the ith contingent leg period  */
	DateTimeArraySP	lpaymentDateArray;

	/** array indicating the discount factors at each simulation timeline date  */
	DoubleArray cLegDiscountFactors;

	/** array indicating whether the ith period belongs to any contingent leg period and hence loss should be paid */
	BoolArray payLoss;

	/**
		size of this array would be the product timeline index;
		a loss occuring in a product timeline grid gets paid in a timeline grid that is not necessarily the same;
		this array maps a product timeline index to another product timeline index where it gets paid;

	*/
	IntArray contLegDiscFactorIndex;

	/** array indicating whether the ith period belongs to any contingent leg period and hence loss should be paid */
	DoubleArray contingentLegScalingFactors;

// ##  The following members would be used for the fee leg computations

	/** Number of Fee leg periods that are credit contingent */
	int numRiskyFeeLegPeriods;

	/** Change in notionals per risky fee leg period (due to credit events)
		This array would get refreshed in every call to the payoff(..) method
	*/
	DoubleArray notionalChangePerPeriod;

/** Size of the array would be the length of the product timeline
	The ith element of the array would indicate the risky fee leg period to which the ith product timeline date belongs to.
	An element could be -1 if a product timeline date does not correspond to a risky fee leg period
*/
	IntArray productIndexToRiskyFeeLegPeriod;

// ##  Other member variables
	/** Generalised CDO instrument to be priced using MC  */
	GeneralisedCDOConstSP cdo;

	/** State Variable for the CDO instrument  */
	ICreditLossConfigIndexedSVSP indexedSV;

	/** Generator needed to update the above state variable  */
	ICreditLossConfigIndexedSVGenConstSP indexedSVGen;

	/** State Variable for the Discount factor  */
	SVDiscFactorSP discFactorSV;

	/** State Variable generator for the Discount factor  */
	SVGenDiscFactorConstSP discFactorSVGen;

	/** State Variable for the feeleg of the CDO instrument  */
	ICreditFeeLegSVSP creditFeeLegSV;

	/** Generator needed to update the above state variable  */
	ICreditFeeLegSVGenConstSP creditFeeLegSVGen;

	/** Product timeline */
	DateTimeArrayConstSP productTimeline;

	/** Fee Leg risky observation dates - modified for couponNotionalType */
	DateTimeArrayConstSP modifiedFeeLegObservationDates;

// ##  member variables to store results for each path

	/** Tranche Fee Leg Risk free price  */
	MCPricesSimple trancheFeeLegRFPrice;

	/** Tranche Fee Leg Price  */
	MCPricesSimple trancheFeeLegPrice;

	/** Tranche Fee Leg Price  */
	MCPricesSimple unitCouponTrancheFeeLegPrice;

	/** Tranche Contingent Leg Price  */
	MCPricesSimple trancheContingentLegPrice;

};
// ################################################################################################################

DRLIB_END_NAMESPACE

#endif
