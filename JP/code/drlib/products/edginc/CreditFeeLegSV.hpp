//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditFeeLegSV.hpp
//
//   Description : Base class for ICreditFeeLeg::SVGen::SV
//
//   Date        : Sep 2006
//
//----------------------------------------------------------------------------

#ifndef CREDIT_FEE_LEG_SV_HPP
#define CREDIT_FEE_LEG_SV_HPP

//#include "edginc/IQMCStateVariable.hpp"
#include "edginc/ICreditFeeLeg.hpp"

DRLIB_BEGIN_NAMESPACE

/** The state variable generator corresponding to the CreditFeeLeg
*/
class PRODUCTS_DLL CreditFeeLegSVGen	:	public virtual ICreditFeeLeg::ISVGen
{
public:
//nested inner class representing the SV
	FORWARD_DECLARE(SV)

	/** constructor */
	CreditFeeLegSVGen(
		int productTimelineSize,
		const IntArray& dateToProductTimelineIndexMap,
		const IntArray& dateToDiscFactorIndex,
		double lossConfigNotional,
		const DoubleArray& scalingFactors,
		const DoubleArray& historicalScalingFactors,
		const DateTimeArray& riskyObservationDates,
		const CouponNotionalTypesArray& couponNotionalTypes,
		const AccrualPeriodArray& accrualPeriods,
		const DoubleArray& coupons,
		const DateTimeArray& payDates,
		const CashFlowArraySP& riskFreeCashFlows
		);

	virtual ~CreditFeeLegSVGen() {}

	/** Same as create but avoids dynamic cast */
	virtual ICreditFeeLegSVSP createNewSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const;

// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const; //the second argument above is PathGenerator

	ICreditFeeLegSVSP sv;
};


/** class representing the state variable */
class CreditFeeLegSVGen::SV	:	public virtual ICreditFeeLeg::ISVGen::ISV
{
public:
	/** constructor  */
	SV(
		int productTimelineSize,
		const IntArray& dateToProductTimelineIndexMap,
		const IntArray& dateToDiscFactorIndex,
		double lossConfigNotional,
		const DoubleArray& scalingFactors,
		const DoubleArray& historicalScalingFactors,
		const DateTimeArray& riskyObservationDates,
		const CouponNotionalTypesArray& couponNotionalTypes,
		const AccrualPeriodArray& accrualPeriods,
		const DoubleArray& coupons,
		const DateTimeArray& payDates,
		const CashFlowArraySP& riskFreeCashFlows
		);

	/** virutal destructor */
	virtual ~SV() {}

	/** method of the parent interface, IStateVariable */
	virtual bool doingPast() const;

	/** computes the pv */
	virtual double pv(
		double& riskFreeCashFlow,
		double& unitCouponPV,
		bool computeUnitCouponPV,
		const DateTime& valueDate,
		const CreditLossConfigIndexedSVResultPath& svResult,
		const SVDiscFactor& discountCurve
		) const ;

private:
	/** a map that facilitates quick conversion from any Date to a productTimeline index. */
	IntArray dateToProductTimelineIndexMap;


	/** Notional corresponding to the lossConfig
		This notional would correspond to the notionalChange reported by SVResult
	*/
	double lossConfigNotional;

	/** Pre-calculated for speed
		Coupon Times YearFaction for the accrual periods
		scalingFactor = Coupon X YearFraction /lossConfigNotional * notional
	*/
	DoubleArray scalingFactors;

	/** scaling factors for periods before the value date  */
	DoubleArray historicalScalingFactors;

	/** accrual periods  */
	AccrualPeriodArray accrualPeriods;

	/** coupons  */
	DoubleArray coupons;

	/** risky observation dates  */
	DateTimeArray riskyObservationDates;

	/** wheher OBSERVATION_DATE or AVERAGE  */
	CouponNotionalTypesArray couponNotionalTypes;

	/** Pay Dates of different accrual periods  */
	DateTimeArray payDates;

	/** risk free cashflows  */
	CashFlowArraySP riskFreeCashFlows;

	/** To ease computation
		This corresponds to the dates on the product timeline
		The MC routine would update the notional at each point on the product timeline
	*/
	mutable DoubleArray notionalPerProductTimeline;

	/** To ease computation
		Would store the mapping between a date and a index that can be used to lookup
		on the discountFactorSV
	*/
	IntArray dateToDiscFactorIndex;
};

DRLIB_END_NAMESPACE

#endif //CREDIT_FEE_LEG_SV_BASE_HPP
