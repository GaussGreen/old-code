//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : The interface that all types of 'credit loss' representations
//                 implement. Examples include a single name, a 'tranche'  
//                 (which here means a portfolio of names and a pair of strikes)
//                 and a portfolio of names where the loss is given by the first 
//                 to default.
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICREDITLOSSCONFIG_HPP
#define QLIB_ICREDITLOSSCONFIG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/AccrualPeriod.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/MarketWrapper.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/StateVariable.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/DateTimeLite.hpp"
#include "edginc/GenericPath.hpp"
#include <list>

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IRebateCalculator);
FORWARD_DECLARE(IProtectionProvider);
FORWARD_DECLARE(CreditEngineParameters);
FORWARD_DECLARE(IForwardRatePricer);

FORWARD_DECLARE (ICreditLossConfig);

class MARKET_DLL ICreditLossConfig : virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;
    ICreditLossConfig();
    virtual ~ICreditLossConfig();

    /** All instances have a 'name' which can be viewed as the name of
        the portfolio or the name of the single curve. */
    virtual string getName() const = 0;

    /** Returns valuedate */
    virtual DateTime getToday() const = 0;

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
         const IProtectionProvider* const protect) const = 0;

    /** Returns the notional reductions on the fee leg (due to losses and/or
        recovered notional) and the corresponding fee rebates if any. In the
        presence of shorts the 'losses' can be negative
        The accrualPeriods are potentially required to determine the actual  
        date of the notional reductions. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        CIntConstSP               triggerDelay,
        CIntConstSP               defaultToCalculationDelay,
        const double              temporaryLossAmount,
        const DateTime&           lastTriggerDate,    
        AccrualPeriodArrayConstSP accrualPeriods,
        IBadDayAdjusterConstSP    bda,
        const bool                recoverNotional) const = 0;

    /** Compute the rebate payments caused by historic fee leg reductions */
    virtual CashFlowArraySP historicRebatePayments(
        const IRebateCalculator* const   rebateCalc,
        FeeLegReductionPerDefaultArraySP reductions,
        IForwardRatePricerSP             model,
        const bool                       recoverNotional) const = 0;

    /** Returns the engine parameters for this loss cfg (NB Model needs to
        have selected appropriate set during getMarket()). The parameter 
        can specify what type of engine parameters are required. An 
        exception is thrown if the relevant type is not found */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const = 0;

    /** Returns the number of "inner loss configs" contained in this LossConfig.
        The need for this is driven by models that need to know the number of
        underlying "assets", eg, CreditMetricsModel will use slow or fast 
        convolution depending on the number of names.
        Note that composite LossConfigs will typically cascade down the calls
        to this method. LossConfigs of type CDOPortfolio, however, will return
        the number of LossConfigs in the portfolio. */
    virtual int numInnerLossConfigs() const = 0;

    /** Returns the inner loss config number "index". */
    virtual ICreditLossConfigConstSP getInnerLossConfig(const int index) const = 0;

    /** Returns the maximum (realistic) loss that may be produced by this 
        loss config. Equivalent to lossGivenDefault for single names.
        "Realistic" in the sense that this method does NOT assume, eg,
        RR=0, since otherwise this method would be similar to notional() */
    virtual double maxPossibleLoss() const = 0;

    /** Returns the loss config notional, ie, the amount that could be 
        (or could have been) lost in the worst possible scenario, even if 
        this scenario is no longer possible (e.g., if a name has defaulted
        with a RR > 0 it will still be accounted as if RR=0 if the losses
        are greater in this case). */
    virtual double notional() const = 0;

    /** Computes the loss and recovered notional currently produced by this 
        loss config. Note both can be 0 if, say, this is a name that has not
        defaulted */
    virtual void currentLossAndRecovered(
        double& loss,               // (O)
        double& recoveredNotional,  // (O)
        const bool recoverNotional) const = 0;

    /** Computes how much of this loss config's original notional has been
        prepaid (amortized) before "prepay date" */
    virtual double getPrepaidNotional(
        const DateTime& prepayDate) const = 0;

    /** Inner class: declaration of State Variable Gen */
	FORWARD_DECLARE(ISVGen);

    /** Inner class: declaration of State Variable Gen */
	FORWARD_DECLARE(IIndexedSVGen);

private:
    ICreditLossConfig(const ICreditLossConfig&);
    ICreditLossConfig& operator=(const ICreditLossConfig&);
};

typedef MarketWrapper<ICreditLossConfig> ICreditLossConfigWrapper;


// ############################################################################


/** Represents the SV generator for the ICreditLossConfig  */
class MARKET_DLL ICreditLossConfig::ISVGen : public virtual IStateVariableGen,
											public virtual IStateVariableClient,
											public virtual VirtualDestructorBase

{
public:
	/** forward declaration of the inner class representing the actual state variable   */
	FORWARD_DECLARE(ISV);

	/** constructor */
	ISVGen() {}

	/** virtual destructor */
	virtual ~ISVGen() {}

    /** Actually creates and returns a new instance of the SV */
	virtual smartPtr<ICreditLossConfig::ISVGen::ISV> createNewSV(
		IStateVariableGen::IStateGen* stateGen) const = 0;

	/** Same as create but avoids dynamic cast */
	virtual smartPtr<ICreditLossConfig::ISVGen::ISV> getLossSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const = 0;

// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const = 0; //the second argument above is PathGenerator

// ##  Methods of the parent interface, IStateVariableClient
    virtual void collectStateVars(IStateVariableCollectorSP svdbase) const = 0;

private:
    ISVGen(const ISVGen&);
    ISVGen& operator=(const ISVGen&);

}; // ICreditLossConfig::SVGen

typedef smartPtr<ICreditLossConfig::ISVGen> ICreditLossConfigSVGenSP;
typedef smartConstPtr<ICreditLossConfig::ISVGen> ICreditLossConfigSVGenConstSP;
typedef vector<ICreditLossConfigSVGenConstSP> ICreditLossConfigSVGenVector;

// ############################################################################

/** Interface for state variable for the ICreditLossConfig  */
class MARKET_DLL ICreditLossConfig::ISVGen::ISV : 	public virtual IStateVariable
{
public:
	/** Nested State variable result*/
	FORWARD_DECLARE(SVResult);

// ## Methods of ICreditLossConfig::SVGen::SV follow

	/** constructor */
	ISV() {}

	/** virtual destructor */
	virtual ~ISV() {};

	/** access method for the results */
	virtual const list<SVResult>& getResults() const = 0;

	/** read/write access method for the results */
	virtual list<SVResult>& results() = 0;

	/** access method for the product timeline */
	virtual const DateTimeLiteVectorConstSP& getTimeline() const = 0;

	/** access method for the notional */
	virtual double getNotional() const = 0;

	/** method of the parent interface, IStateVariable */
	virtual bool doingPast() const = 0;

private:
    ISV(const ISV&);
    ISV& operator=(const ISV&);

};  // ICreditLossConfig::SVGen::SV

typedef smartPtr<ICreditLossConfig::ISVGen::ISV> ICreditLossConfigSVSP;
typedef smartConstPtr<ICreditLossConfig::ISVGen::ISV> ICreditLossConfigSVConstSP;
typedef vector<ICreditLossConfigSVSP> ICreditLossConfigSVVector;

// ############################################################################

/** Represents state variable result for the ICreditLossConfig
	It is the snap-shot at a time due to a credit loss event.
*/
class MARKET_DLL ICreditLossConfig::ISVGen::ISV::SVResult
{
public:
	/** constructor */
	SVResult()
	:eventTime(DateTimeLite()),
	lossAmount(0),
	notionalAmount(0),
	bucketIndex(0) {}

	/** another, more meaningful constructor */
	SVResult(
		const DateTimeLite& eTime,
		double lAmount,
		double nAmount,
		int    bIndex)
		:eventTime(eTime),
		lossAmount(lAmount),
		notionalAmount(nAmount),
		bucketIndex(bIndex) {}

	/** virtual destructor */
	virtual ~SVResult() {}

    /** less than operator */
	bool operator<(const SVResult& rhs) const
	{
		return (eventTime < rhs.eventTime);
	}

// ## Member variables follow

	/** time in the timeline for which we are fetching the loss status */
	DateTimeLite eventTime;

	/** change in losses at the event time due to the credit event  */
	double lossAmount;

	/** change in notional at the event time due to the credit event */
	double notionalAmount;

	/** defines the interval in a time line grid where the event time belongs   */
	int    bucketIndex;
};

typedef ICreditLossConfig::ISVGen::ISV::SVResult CreditLossConfigSVResult;
typedef list<CreditLossConfigSVResult> CreditLossConfigSVResultList;

// ############################################################################

/** Represents the indexed SV generator for the ICreditLossConfig */
class MARKET_DLL ICreditLossConfig::IIndexedSVGen :	public virtual IStateVariableGen,
													public virtual IStateVariableClient,
													public virtual VirtualDestructorBase
{
public:
	/** Nested State variable*/
	FORWARD_DECLARE(ISV);

	/** constructor */
	IIndexedSVGen() {}

	/** virtual destructor */
	virtual ~IIndexedSVGen() {}

	/** Actually creates and returns a new instance of the SV */
    virtual smartPtr<ICreditLossConfig::IIndexedSVGen::ISV> createNewSV(
		IStateVariableGen::IStateGen* stateGen) const = 0;

	/** Same as create but avoids dynamic cast */
    virtual smartPtr<ICreditLossConfig::IIndexedSVGen::ISV> getLossSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const = 0;

// ##  Methods of the parent interface, IStateVariable
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen *stateGen) const = 0;

// ##  Methods of the parent interface, IStateVariableClient
    virtual void collectStateVars(IStateVariableCollectorSP svdb) const = 0;

private:
    IIndexedSVGen(const IIndexedSVGen&);
    IIndexedSVGen& operator=(const IIndexedSVGen&);

}; // ICreditLossConfig::IndexedSVGen

typedef smartPtr<ICreditLossConfig::IIndexedSVGen> ICreditLossConfigIndexedSVGenSP;
typedef smartConstPtr<ICreditLossConfig::IIndexedSVGen> ICreditLossConfigIndexedSVGenConstSP;
typedef vector<ICreditLossConfigIndexedSVGenConstSP> ICreditLossConfigIndexedSVGenVector;


// ############################################################################

/** Interface for state variables representing indexed loss event times. Indexed means that we bucket
	the credit event time into a bucket
*/
class MARKET_DLL ICreditLossConfig::IIndexedSVGen::ISV : public virtual IStateVariable
{
public:
	/** Nested Indexed State variable result*/
	FORWARD_DECLARE(SVResult);

	/** Constructor  */
	ISV() {};

	/** Virtual Destructor  */
	virtual ~ISV() {};

	/** access method for the results */
	virtual const GenericPath<SVResult>& getResults() const = 0;
	
	/** read/write access method for the results */
	virtual GenericPath<SVResult>& results() = 0;

	/** size of the results path */
	virtual int getResultPathSize() const = 0;

	/** access method for the notional */
	virtual double getNotional() const = 0;

	/** method of the parent interface, IStateVariable */
	virtual bool doingPast() const = 0;

private:
    ISV(const ISV&);
    ISV& operator=(const ISV&);

}; // ICreditLossConfig::IIndexedSVGen::ISV

typedef smartPtr<ICreditLossConfig::IIndexedSVGen::ISV> ICreditLossConfigIndexedSVSP;
typedef smartConstPtr<ICreditLossConfig::IIndexedSVGen::ISV> ICreditLossConfigIndexedSVConstSP;
typedef vector<ICreditLossConfigIndexedSVSP> ICreditLossConfigIndexedSVVector;

// ############################################################################

/** A struct representing the loss and notional change at a time point
	timeline points */
class MARKET_DLL ICreditLossConfig::IIndexedSVGen::ISV::SVResult
{
public:
	SVResult()
		:index(0),
		lossChange(0),
		notionalChange(0) {}

	SVResult(
		int index,
		double lossChange,
		double notionalChange)
		:index(index),
		lossChange(lossChange),
		notionalChange(notionalChange) {}

	virtual ~SVResult() {}

	/** Index on the product timeline  */
	int index;

	/** change in loss at a time point  */
	double lossChange;

	/** change in notional at a time point  */
	double notionalChange;
};

typedef ICreditLossConfig::IIndexedSVGen::ISV::SVResult CreditLossConfigIndexedSVResult;
typedef vector<CreditLossConfigIndexedSVResult> CreditLossConfigIndexedSVResultVector;
typedef GenericPath<CreditLossConfigIndexedSVResult> CreditLossConfigIndexedSVResultPath;

DRLIB_END_NAMESPACE

#endif
