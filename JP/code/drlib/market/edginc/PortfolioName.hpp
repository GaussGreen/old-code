//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Each of the underlying "names" in, say, a CDO
//
//   Author      : Antoine Gregoire
//
//----------------------------------------------------------------------------


#ifndef QLIB_PORTFOLIONAME_HPP
#define QLIB_PORTFOLIONAME_HPP

#include "edginc/Atomic.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/Theta.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/CreditNameNotionalLevel.hpp"
#include "edginc/ISingleDefaultCreditLossConfig.hpp"
#include "edginc/ITrancheCreditEventOverride.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/CreditLossConfigMC.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_WRAPPER(ITrancheCreditEventOverride);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IRebateCalculator);
FORWARD_DECLARE(IProtectionProvider);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);

class MARKET_DLL PortfolioName:
    public CObject,
    virtual public ISingleDefaultCreditLossConfig,
    virtual public CreditNameNotionalLevel::Shift,
    virtual public Theta::Shift,
	virtual public ICreditLossConfigSVGenMC,
	virtual public ICreditLossConfigIndexedSVGenMC
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~PortfolioName();

    /** basic validation */
    virtual void validatePop2Object();

    /** NotionalLevel::Shift implementation */
    bool sensShift(CreditNameNotionalLevel* shift);
    string sensName(CreditNameNotionalLevel* shift) const;

    /** Returns the notional */
    double getNameNotional() const;

    /**
     * Returns the recovery:
     * - field "nameRecovery" if defaultParamOverride=TRUE
     * - market recovery otherwise
     *  */
    double getNameRecovery() const;

    /**
     * Returns the default date:
     * - field "nameDefDate" if specified
     * - market default date (if any) otherwise
     *  */
    DateTime getNameDefaultDate() const;

    /** Indicates whether the name has defaulted before the name-specific
     * "protectionStartDate" - returns false if:
     * - the name has not defaulted, or
     * - there is no name-specific protectionStartDate, or
     * - the name has defaulted on or after the protectionStartDate */
    bool defaultedBeforeProtectionStarts() const;

    /** Returns the underlying credit asset */
    CreditAssetConstSP getCreditAsset() const;

    /** Returns beta */
    double getBeta() const;

    /** Returns protection start date */
    const DateTime& getProtectionStartDate() const;

    /** Returns the min of last date parameter and protection end date
        for specified name (if any). This if for name 'cutoff' */
    const DateTime& getProtectionEndDate(const DateTime& lastDate) const;

    /** Returns the name maturity cut off for the name **/
	const DateTime& getNameMaturityCutOff() const;

    /** Returns the ranges for Loss Given Default for this name.
        In the particular, the LGD of a name is given by
        nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
    void nameLGDRanges(double& lgdNotional, /* (O) */
                       double& lgdFloor,    /* (O) */
                       double& lgdCap)      /* (O) */ const;

    /** Returns
        nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
    double lossGivenDefault() const;

    /** Computes betaTweak = beta + shift * (1 - beta) */
    static double betaTweak(double beta, double shift, bool capAndFloor = true);

    /** Indicates whether this name has a credit event override */
    bool hasEventOverride() const;

    /** Indicates whether this name has a default param override */
    bool defaultParamOverrideExists() const;

    //++++++++ ICreditLossConfig methods
    /** Name of the single curve. */
    virtual string getName() const;

    DateTime getToday() const;

    void impliedParSpreadsAndDurations(
        const YieldCurveConstSP discount,
        const DateTimeArray& dates,
        DoubleArray& impliedSpreads,   /* (Output) */
        DoubleArray& durations) const; /* (Output) */

    /** Returns the historic losses on the contingent leg, due to the
        settlement of this name's default. In the presence of shorts
        the 'losses' can be negative.
        If there are no losses, it may return an empty (0 size)
        CtgLegLossPerDefaultArraySP or a "null" SP.
        Note: the cashflow dates are the CALCULATION DATES (as opposed to
        settlement dates) when the corresponding losses are determined. */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
        CIntConstSP triggerDelay,
        CIntConstSP defaultToCalculationDelay,
        const DateTime& lastTriggerDate,
        IBadDayAdjusterConstSP bda,
        const IProtectionProvider* const protect) const;

    /** Returns the notional reductions on the fee leg (due to losses and/or
        recovered notional) and the corresponding fee rebates if any. In the
        presence of shorts the 'losses' can be negative
        If there are no losses/recovered notional, it may return empty (0 size)
        arrays or "null" SPs.
        The accrualPeriods are potentially required to determine the actual
        date of the notional reductions. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        CIntConstSP               triggerDelay,
        CIntConstSP               defaultToCalculationDelay,
        const double              temporaryLossAmount,
        const DateTime&           lastTriggerDate,
        AccrualPeriodArrayConstSP accrualPeriods,
        IBadDayAdjusterConstSP    bda,
        const bool                recoverNotional) const;

    /** Compute the rebate payments caused by historic fee leg reductions */
    virtual CashFlowArraySP historicRebatePayments(
        const IRebateCalculator* const   rebateCalc,
        FeeLegReductionPerDefaultArraySP reductions,
        IForwardRatePricerSP             model,
        const bool                       recoverNotional) const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the engine parameters for this name (NB Model needs to
        have selected appropriate set during getMarket()). The parameter
        can specify what type of engine parameters are required. An
        exception is thrown if the relevant type is not found */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;

    /** Returns the number of "inner loss configs" contained in this LossConfig.
        In this case returns '1'. */
    virtual int numInnerLossConfigs() const;

    /** Returns the inner loss config number "index".
        In this case throws an exception. */
    virtual ICreditLossConfigConstSP getInnerLossConfig(const int index) const;

    /** Returns lossGivenDefault() */
    virtual double maxPossibleLoss() const;

    /** Returns nameNotional */
    virtual double notional() const;

    /** Computes the loss and recovered notional currently produced by this
        name. Note both can be 0 if this name has not defaulted */
    virtual void currentLossAndRecovered(
        double& loss,               // (O)
        double& recoveredNotional,  // (O)
        const bool recoverNotional) const;

    /** Computes how much of this loss config's original notional has been
        prepaid (amortized) before "prepay date" */
    virtual double getPrepaidNotional(const DateTime& prepayDate) const;

    //-------- ICreditLossConfig methods

    //++++++++++++++++++++++++++++++++++++++++
    // ISingleDefaultCreditLossConfig methods
    //
    /** Returns true if this loss config has defaulted */
    virtual bool hasDefaulted() const;

    /** Returns the default date, or an empty DateTime if the loss config
        has not defaulted */
    virtual DateTime getDefaultDate() const;

    /** Returns the timeline used in this loss config, ie, the timepoints
        defining the regions at which it is safe to assume that the (flat
        forwards) default rate is constant. */
    virtual DateTimeArraySP getTimeLine() const;
    //
    // ISingleDefaultCreditLossConfig methods
    //----------------------------------------

    /** Theta::IShift method */
    virtual bool sensShift(Theta* shift);

	//++++++++ ICreditLossConfigSVGenMC methods
	virtual ICreditLossConfigSVGenConstSP createSVGen(
		const DateTimeLiteVectorConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc,
        const bool recoverNotional
		) const;

	//-------- ICreditLossConfigSVGenMC methods

	//-------- ICreditLossConfigIndexedSVGenMC methods
	virtual ICreditLossConfigIndexedSVGenConstSP createIndexedSVGen(
		const DateTimeArrayConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc,
        const bool recoverNotional
		) const;

	//-------- ICreditLossConfigIndexedSVGenMC methods

// Declare nested classes representing State Variable Generators corresponding to this loss config
	FORWARD_DECLARE(SVGen);
	FORWARD_DECLARE(IndexedSVGen);

private:
    /** Constructor */
    PortfolioName();

    /** Default constructor (for reflection) */
    static IObject* defaultConstructor();

    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    // ------
    // FIELDS
    // ------

    /** Notional of the name */
    double nameNotional;

    /* LGD = Loss Given Default
     * LGD = NameNotional * max(lgdFloor,min(lgdCap, lgdNotional - recovery))
     * Default values are :
     * lgdNotional = 1
     * lgdFloor    = 0
     * lgdCap      = 1 */
    double lgdNotional;
    double lgdFloor;
    double lgdCap;

    /** credit asset */
    CreditAssetWrapper asset;

    // Beta and Decretion beta are now optional: they are MODEL parameters
    // and as such they have been moved to the CreditEngineParameters class.
    // They are still here only for backwards compatibility reasons.
    /** beta */
    CDoubleSP beta;
    /** decretion beta */
    CDoubleSP decretionBeta;

    /*
     * Maturity cut off date of the name.
     * N if no cutoff.
     * */
    DateTime nameMatCutoff;

    /*
     * Say whether we get recovery
     * info and default date from mkt
     * */
    bool defaultParamOverride;

    /** Recovery in % override */
    double nameRecovery;

    /*
     * Default dates of names in
     * portfolio. empty if not defaulted
     * */
    DateTime nameDefDate;

    /**
     * Optional protection start date for the name.
     * If protectionStartDate is empty : uses "today"
     * If protectionStartDate <  today : uses "today"
     * If protectionStartDate >= today : uses "protectionStartDate"
     *
     * Initial motivation for that field is to book tranches on forward
     * starting CDS (called "Forward Tranches")
     * */
     DateTimeSP protectionStartDate;

    /** Override for credit event parameters */
    ITrancheCreditEventOverrideSP creditEventOverride;

    DateTime valueDate;
public:

	virtual DateTime getValueDate() const {return valueDate;};

	virtual double getRecovery() const
	{
		return
			defaultParamOverride?
					nameRecovery:
					asset->recoveryRate();
	};

	void AdjustNotionalAndDates(
		double newNotional,
		DateTime newStart,
		DateTime newCutOff);

	class MARKET_DLL NameView : public CObject {

		string name;
		double notional;
		DateTime protectionStart;
		DateTime protectionEnd;
		double recovery;

	public:
		virtual ~NameView() {};

		static CClassConstSP const TYPE;

		NameView(
			CClassConstSP clazz = TYPE)
			:CObject(clazz) {};

		NameView(
			string name,
			double notional,
			double recovery,
			DateTime protectionStart,
			DateTime protectionEnd,
			CClassConstSP clazz = TYPE)
			:CObject(clazz),
			name(name),
			notional(notional),
			recovery(recovery),
			protectionStart(protectionStart),
			protectionEnd(protectionEnd) {};

		static IObject* defaultNameView()
		{
			return new NameView();
		};

		static void load(CClassSP& clazz);

	};

	DECLARE(NameView);

    //	DECLARE(NameViewArray);

	NameViewSP getNameView() const;

};

DECLARE(PortfolioName);

// ################################################################################################################

/** StateVariable Generator corresponding to the PortfolioName  */
class MARKET_DLL PortfolioName::SVGen:	public virtual ICreditLossConfig::ISVGen,
										public virtual IElemStateVariableGen
{
public:
	/** forward declaration of the inner class representing the actual state variable   */
	FORWARD_DECLARE(SV);

	/** constructor */
	SVGen(){}

	/** another more meaningful constructor */
	SVGen(
		const PortfolioNameConstSP&  portfolioName,
		const DateTimeLiteVectorConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc
		);

	/** virtual destructor */
	virtual ~SVGen()
	{}

	/** access methods for name */
	virtual const string& getName() const;

	/** Actually creates and returns a new instance of the SV */
	virtual ICreditLossConfigSVSP createNewSV(IStateVariableGen::IStateGen* stateGen) const;

private:
	/** same as create but avoids dynamic cast */
	virtual ICreditLossConfigSVSP getLossSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const;


// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const;

// ##  Methods of the parent interface, IStateVariableClient
    virtual void collectStateVars(IStateVariableCollectorSP svdbase) const;

// ##  Methods of the parent interface, IElemStateVariableGen
	virtual void attachSVGen(class IElemStateVariableGenVisitor* sv) const;

    /** Just the declaration of the copy constructor; to avoid compiler creating them   */
    SVGen(const SVGen&);

    /** Just the declaration of the equal to operator; to avoid compiler creating them   */
    SVGen& operator=(const SVGen&);

// ##  Member variables follow

	/** Storage for the product timeline */
	DateTimeLiteVectorConstSP timeline;

	/** name of the portfolio */
	string name;

	/** Portfolio Name  - to be used for calculating Past Default Losses by calling methods on PortfolioName */
	PortfolioNameConstSP portfolioName;

	/** Number of days gap between the Credit Event Date (CED) and the Event Determination Date (EDD)
		To be used for calculating Past Default Losses by calling methods on PortfolioName */
	CIntConstSP triggerDelay;

	/**  Number of days gap between the CED and the Calculation Date (CD)
	*/
	CIntConstSP defaultToCalculationDelay;

	/**  Loss amount to be used when actual recovery amount is not known
	*/
	double temporaryLossAmount;

	/**  Date by which the default must be triggered
	*/
	DateTime lastTriggerDate;

	/**  ???
	*/
	AccrualPeriodArrayConstSP accrualPeriods;

	/**  Typically the CDO */
	IBadDayAdjusterConstSP bda;

	/**  Typically the CDO */
	IProtectionProviderConstSP protect;

	/**  Typically the CDO */
	IRebateCalculatorConstSP rebateCalc;

};

typedef smartPtr<PortfolioName::SVGen> PortfolioNameSVGenSP;
typedef smartConstPtr<PortfolioName::SVGen> PortfolioNameSVGenConstSP;

// ################################################################################################################
/** Interface for state variables representing loss event times  */
class MARKET_DLL PortfolioName::SVGen::SV : 	public virtual ICreditLossConfig::ISVGen::ISV
{
public:
	friend class PortfolioName::SVGen;

	/** Nested sub result - part of the Portfolio that registers with the MC engine */
	FORWARD_DECLARE(SVSubResult);

	/** constructor */
	SV() {}

	/** another, more meaningful constructor */
	SV(	const string& name,
		double notional,
		const CreditAssetConstSP& asset,
		double beta,
		const CDoubleConstSP& recoveryOverride,
		const DateTime& protectionStartDate,
		const DateTime& maturityCutOff,
		const CtgLegLossPerDefaultArrayConstSP& histContLosses,
		const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
		const FeeLegReductionPerDefaultArrayConstSP& feeRecovered,
		const DateTimeLiteVectorConstSP& timeline
		);

	/** virtual destructor */
	virtual ~SV() {};

	/** access methods for name */
	virtual const string& getName() const;

	/** read access method for the sub results */
	virtual const refCountPtr<const vector<SVSubResult> >& getSubResults();

	/** write access method for the sub results */
	virtual void setSubResults(const refCountPtr<const vector<SVSubResult> >&  subResults) const;

	/** uses the subResults and builds the lossEvents */
	virtual void calculateLossEvents() const;

	/** access method for the notional */
	virtual double getNotional() const;

	/** access method for the credit asset  */
	virtual const CreditAssetConstSP& getAsset() const;

	/** access method for the beta  */
	virtual double getBeta() const;

	/** access method for the product timeline */
	virtual const DateTimeLiteVectorConstSP& getTimeline() const;

protected:

	/** read access method for the results */
	virtual const CreditLossConfigSVResultList& getResults() const;

	/** read/write access method for the results */
	virtual CreditLossConfigSVResultList& results();

	/** method of the parent interface, IStateVariable */
	virtual bool doingPast() const;

	/** Just the declaration of the copy constructor; to avoid compiler creating them   */
	SV(const SV&);

    /** Just the declaration of the equal to operator; to avoid compiler creating them   */
    SV& operator=(const SV&);

// ##  Member variables below

	/** storage for product time lines */
	DateTimeLiteVectorConstSP timeline;

	/** name */
	string name;

	/** notional */
	double notional;

	/** credit asset */
    CreditAssetConstSP asset;

	/** beta of this name */
	double beta;

	/** historic losses due to default of this name */
	CtgLegLossPerDefaultArrayConstSP histContLosses;

	/** notional reduction due to loss on account of default of this name */
	FeeLegReductionPerDefaultArrayConstSP feeLosses;

	/** notional reduction due to recovery loss on account of default of this name */
	FeeLegReductionPerDefaultArrayConstSP feeRecovered;

	/** Protection Start Date */
	DateTime protectionStartDate;

	/** Maturity Cutoff date */
	DateTime maturityCutOff;

    /** Recovery of this name if there is a recovery rate override */
	CDoubleConstSP recoveryOverride;

    /** storage for the CreditLossConfigSVResults  */
	mutable CreditLossConfigSVResultList	 creditLossConfigSVResults;

	/** sub result */
	mutable refCountPtr<const vector<SVSubResult> > subResults;

}; // PortfolioName::SVGen::SV

typedef smartPtr<PortfolioName::SVGen::SV> PortfolioNameSVSP;
typedef smartConstPtr<PortfolioName::SVGen::SV> PortfolioNameSVConstSP;

// ################################################################################################################
/** Interface for part of the PorfolioNameSV that can register with the engine
*/
class MARKET_DLL PortfolioName::SVGen::SV::SVSubResult
{
public:
	/** constructor */
	SVSubResult() {}

	/** another, more meaningul constructor */
	SVSubResult(
		const DateTimeLite& defaultTime,
		double recovery):
		defaultTime(defaultTime),
			recovery(recovery) {}

	/** virtual destructor */
	virtual ~SVSubResult() {}

// ## member variables follow
	/** default date returned by the MC engine */
	DateTimeLite defaultTime;

	/** recovery of the porfolio name returned by the MC engine*/
	double recovery;
};

typedef PortfolioName::SVGen::SV::SVSubResult PortfolioNameSVSubResult;
typedef vector<PortfolioNameSVSubResult> PortfolioNameSVSubResultVector;
typedef refCountPtr<PortfolioNameSVSubResultVector> PortfolioNameSVSubResultVectorSP;
typedef refCountPtr<const PortfolioNameSVSubResultVector> PortfolioNameSVSubResultVectorConstSP;



// ################################################################################################################

/** Indexed StateVariable Generator corresponding to the PortfolioName  */
class MARKET_DLL PortfolioName::IndexedSVGen:	public virtual ICreditLossConfig::IIndexedSVGen,
												public virtual IElemStateVariableGen
{
public:
	/** forward declaration of the inner class representing the actual state variable   */
	FORWARD_DECLARE(SV);

	/** constructor */
	IndexedSVGen() {}

	/** another, more meaningful constructor */
	IndexedSVGen(
		const DateTimeArrayConstSP& timeline,
		const PortfolioNameConstSP&  portfolioName,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc
		);

	/** virtual destructor */
	virtual ~IndexedSVGen()
	{}

	/** access methods for name */
	virtual const string& getName() const;

	/** Actually creates and returns a new instance of the SV */
	virtual ICreditLossConfigIndexedSVSP createNewSV(IStateVariableGen::IStateGen* stateGen) const;

private:
	/** same as create but avoids dynamic cast */
	virtual ICreditLossConfigIndexedSVSP getLossSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const;

// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const;

// ##  Methods of the parent interface, IStateVariableClient
    virtual void collectStateVars(IStateVariableCollectorSP svdb) const;

// ##  Methods of the parent interface, IElemStateVariableGen
	virtual void attachSVGen(class IElemStateVariableGenVisitor* sv) const;

// ##  Member variables follow

	/** MC product timeline */
	DateTimeArrayConstSP timeline;

	/** name of the portfolio */
	string name;

	/** Portfolio Name  - to be used for calculating Past Default Losses by calling methods on PortfolioName */
	PortfolioNameConstSP portfolioName;

	/** Number of days gap between the Credit Event Date (CED) and the Event Determination Date (EDD)
		To be used for calculating Past Default Losses by calling methods on PortfolioName */
	CIntConstSP triggerDelay;

	/**  Number of days gap between the CED and the Calculation Date (CD)
	*/
	CIntConstSP defaultToCalculationDelay;

	/**  Loss amount to be used when actual recovery amount is not known
	*/
	double temporaryLossAmount;

	/**  Date by which the default must be triggered
	*/
	DateTime lastTriggerDate;

	/**  ???
	*/
	AccrualPeriodArrayConstSP accrualPeriods;

	/**  Typically the CDO */
	IBadDayAdjusterConstSP bda;

	/**  Typically the CDO */
	IProtectionProviderConstSP protect;

	/**  Typically the CDO */
	IRebateCalculatorConstSP rebateCalc;

};

typedef smartPtr<PortfolioName::IndexedSVGen::SV> PortfolioNameIndexedSVSP;
typedef smartConstPtr<PortfolioName::IndexedSVGen::SV> PortfolioNameIndexedSVConstSP;

// ################################################################################################################
/** Interface for state variables representing loss event times  */
class MARKET_DLL PortfolioName::IndexedSVGen::SV : 	public virtual ICreditLossConfig::IIndexedSVGen::ISV
{
public:
	friend class PortfolioName::IndexedSVGen;

	/** Nested inner class: SVSubResult - part that registers with the MC engine */
	FORWARD_DECLARE(SVSubResult);

	/** constructor */
	SV() {}

	/** another, more meaningful constructor */
	SV(	const DateTimeArrayConstSP& timeline,
		const string& name,
		double notional,
		const CreditAssetConstSP& asset,
		double beta,
		bool hasDefaulted,
		double lossGivenDefault,
		const DateTimeConstSP& defaultSettDate,
		const CDoubleConstSP& recoveryOverride,
		int protectionStartDateIndex,
		int maturityCutOffIndex,
		const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
		const FeeLegReductionPerDefaultArrayConstSP& feeRecovered
		);

	/** virtual destructor */
	virtual ~SV() {};

	/** access method for the defaultSettlement Date */
	virtual const DateTimeConstSP& getDefaultSettDate() const;

	/** This method would be called by the engine.
		Sets the size of the GenericPath that is returned when getResults() is called
		Also, sets the pastDefaultTimeIndex. The SV does not know the engine timeline. So the engine
		first calls getDefaultSettDate() and then converts the date into an index into its timeline and
		sets it to the SV.
	*/
	virtual void initializePath(const CIntConstSP& pastDefaultIndex, int productTimelineIndex,  int pathSize);

	/** access methods for name */
	virtual const string& getName() const;

	/** read access method for the sub result */
	virtual const refCountPtr<const vector<SVSubResult> >& getSubResults();

	/** write access method for the sub result */
	virtual void setSubResults(const refCountPtr<const vector<SVSubResult> >&  subResults) const;

	/** uses the subResults and builds the lossEvents
		will obviously, throw an exception if subResults are not set
	*/
	virtual void calculateLossEvents() const;

	/** access method for the notional */
	virtual double getNotional() const;

	/** access method for the credit asset  */
	virtual const CreditAssetConstSP& getAsset() const;

	/** access method for the beta  */
	virtual double getBeta() const;

protected:

	virtual int getResultPathSize() const;

	/** read access method for the results */
	virtual const CreditLossConfigIndexedSVResultPath& getResults() const;

	/** read/write access method for the results */
	virtual CreditLossConfigIndexedSVResultPath& results();

	/** method of the parent interface, IStateVariable */
	virtual bool doingPast() const;

// ##  Member variables below
	/** storage for the CreditLossConfigSV Results  */
	mutable CreditLossConfigIndexedSVResultVector creditLossConfigIndexedSVResultVector;

	/** GenericPath representing the value of sv results along a monte carlo path */
	mutable CreditLossConfigIndexedSVResultPath	 creditLossConfigIndexedSVResultPath;

	/** indicates whether the member "creditLossConfigIndexedSVResultVector" has been initialized  */
	mutable bool resultPathInitialized;

	/** to ease MC computation. We keep track of the defaultIndices that were updated by the engine in the last iteration  */
	mutable vector<int> defaultIndices;

	/** MC product timeline */
	DateTimeArrayConstSP timeline;

	/** name */
	string name;

	/** notional */
	double notional;

	/** credit asset */
    CreditAssetConstSP asset;

	/** beta of the name */
	double beta;

	/**  Indicates whether the asset has defaulted in the past */
	bool hasDefaulted;

	/**  Default settlement date if the asset has defaulted.  If this member is empty and the asset has defaulted, it means
		that the asset did not default between the protection start and end date or beyond the end of timeline.
	*/
	DateTimeConstSP defaultSettDate;

	/** Loss Amount if this asset has defaulted in the past. */
	double lossGivenDefault;

	/** notional reduction due to loss on account of default of this name */
	FeeLegReductionPerDefaultArrayConstSP feeLosses;

	/** notional reduction due to recovery loss on account of default of this name */
	FeeLegReductionPerDefaultArrayConstSP feeRecovered;

	/** Recovery of this name if there is a recovery rate override */
	CDoubleConstSP recoveryOverride;

	/** Protection Start Date index of the name. This is an optional field */
	int protectionStartDateIndex;

	/** Maturity Cutoff date index of the name. This is an optional field */
	int maturityCutOffIndex;

	/** sub result - this part is registered with the engine*/
	mutable refCountPtr<const vector<SVSubResult> > subResults;

}; // PortfolioNameSVGen::PortfolioNameSV

typedef smartPtr<PortfolioName::IndexedSVGen::SV> PortfolioNameIndexedSVSP;
typedef smartConstPtr<PortfolioName::IndexedSVGen::SV> PortfolioNameIndexedSVConstSP;

// ################################################################################################################

/** Interface for part of the PorfolioNameIndexedSV that can register with the engine
*/
class MARKET_DLL PortfolioName::IndexedSVGen::SV::SVSubResult
{
public:
	/** constructor */
	SVSubResult() {}

	/** constructor */
	SVSubResult(
		int defaultTimeIndex, //in the engine timeline; -1 if the asset does not default
		int productTimelineIndex, // in the engine time, would be empty if the asset does not default
		double recovery
		):
		defaultTimeIndex(defaultTimeIndex),
		productTimelineIndex(productTimelineIndex),
		recovery(recovery)
		{}

	/** virtual destructor */
	virtual ~SVSubResult() {}

	/** Default time (engine) index returned by the MC engine*/
	int defaultTimeIndex;

	/** product timeline index */
	int productTimelineIndex;

	/** recovery of the porfolio name returned by the MC engine*/
	double recovery;
};

typedef PortfolioName::IndexedSVGen::SV::SVSubResult PortfolioNameIndexedSVSubResult;
typedef vector<PortfolioNameIndexedSVSubResult> PortfolioNameIndexedSVSubResultArray;
typedef refCountPtr<PortfolioNameIndexedSVSubResultArray> PortfolioNameIndexedSVSubResultArraySP;
typedef refCountPtr<const PortfolioNameIndexedSVSubResultArray> PortfolioNameIndexedSVSubResultArrayConstSP;

// ################################################################################################################


DRLIB_END_NAMESPACE

#endif
