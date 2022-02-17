//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditFeeLegSV.cpp
//
//   Description : Base class for ICreditFeeLeg::SVGen::SV
//
//   Date        : Sep 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CreditFeeLegSV.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
CreditFeeLegSVGen::CreditFeeLegSVGen(
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
	):
	sv(
		new CreditFeeLegSVGen::SV(
			productTimelineSize,
			dateToProductTimelineIndexMap,
			dateToDiscFactorIndex,
			lossConfigNotional,
			scalingFactors,
			historicalScalingFactors,
			riskyObservationDates,
			couponNotionalTypes,
			accrualPeriods,
			coupons,
			payDates,
			riskFreeCashFlows))
	{}


ICreditFeeLegSVSP
CreditFeeLegSVGen::createNewSV(
	IStateVariableSP oldStateVar,
	IStateVariableGen::IStateGen* stateGen) const
{
	return ICreditFeeLegSVSP(sv.get());
}

IStateVariableSP
CreditFeeLegSVGen::create(
	IStateVariableSP oldStateVar,
	IStateVariableGen::IStateGen* stateGen) const
{
	return IStateVariableSP(
		this->createNewSV(
			oldStateVar,
			stateGen));
}

CreditFeeLegSVGen::SV::SV(
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
	):
	dateToProductTimelineIndexMap(dateToProductTimelineIndexMap),
	lossConfigNotional(lossConfigNotional),
	scalingFactors(scalingFactors),
	historicalScalingFactors(historicalScalingFactors),
	accrualPeriods(accrualPeriods),
	coupons(coupons),
	riskyObservationDates(riskyObservationDates),
	couponNotionalTypes(couponNotionalTypes),
	payDates(payDates),
	riskFreeCashFlows(riskFreeCashFlows),
	dateToDiscFactorIndex(dateToDiscFactorIndex)
{
	notionalPerProductTimeline.resize(productTimelineSize, 0);
}

bool
CreditFeeLegSVGen::SV::doingPast() const
{
	return false;
}

double
CreditFeeLegSVGen::SV::pv(
	double& riskFreeCashFlow,
	double& unitCouponPV,
	bool computeUnitCoupon,
	const DateTime& valueDate,
	const CreditLossConfigIndexedSVResultPath& svResult,
	const SVDiscFactor& discountCurve) const
{
	static const string routine = "CreditFeeLegSVGen::SV::pv";
//reset
	for (int i=0; i<notionalPerProductTimeline.size(); ++i)
		notionalPerProductTimeline[i] = 0;

//we collect the notional changes at different product timeline periods
	for (int i=svResult.begin(); i <= svResult.end(); ++i)
	{
		int productIndex = svResult[i].index;
		if (productIndex == -1) continue;
		notionalPerProductTimeline[productIndex] -= svResult[i].notionalChange; //what about scaling??
	}

//convert the notional changes to absolute notionals at different product timeline periods
	notionalPerProductTimeline[0] += lossConfigNotional;

	for (int i=1; i< notionalPerProductTimeline.size(); ++i)
		notionalPerProductTimeline[i] += notionalPerProductTimeline[i-1];

//compute the risky PVs
	double feeLegPV = 0;
	unitCouponPV = 0;

	for (int i=0; i < couponNotionalTypes.size(); ++i)
	{
		double pv = 0;
		if (payDates[i] >= valueDate)
		{
			int discFactorIndex = dateToDiscFactorIndex[payDates[i].getDate()];
			if (discFactorIndex == -1)
				throw ModelException(
					routine,
					"unknown payDate - discount factor not calculated for this date ");

			double discountFactor = discountCurve.getDF(discFactorIndex);

			if (!riskyObservationDates[i].empty())
			{
				int productIndex
					= dateToProductTimelineIndexMap[riskyObservationDates[i].getDate()];
				if (productIndex == -1)
					throw ModelException(
						routine,
						"unknown observation Date - product Index not defined for this date ");
				double notional = notionalPerProductTimeline[productIndex];

				pv = notional *
					scalingFactors[i] *
					discountFactor;

				feeLegPV += pv;

				if (computeUnitCoupon)
					unitCouponPV += pv/coupons[i];
			}

//historical accrual
			if (historicalScalingFactors[i] != 0)
			{
				int productIndex
					= dateToProductTimelineIndexMap[accrualPeriods[i]->startDate().getDate()];
				if (productIndex == -1)
					throw ModelException(
						routine,
						"unknown accrual start date - product Index not defined for this date ");

				double notional = notionalPerProductTimeline[productIndex];

				pv = notional *
					historicalScalingFactors[i] *
					discountFactor;

				feeLegPV += pv;

				if (computeUnitCoupon)
					unitCouponPV += pv/coupons[i];
			}
		}
	}
//compute the riskless fee cashflows
	riskFreeCashFlow = 0;
	if (riskFreeCashFlows.get())
	{
		for (int i=0; i< riskFreeCashFlows->size(); ++i)
		{
			int discFactorIndex = dateToDiscFactorIndex[(*riskFreeCashFlows)[i].date.getDate()];
			if (discFactorIndex == -1)
				throw ModelException(
					routine,
					"unknown riskFree cashflow date  - discount factor not calculated for this date ");

			double discountFactor = discountCurve.getDF(discFactorIndex);

			riskFreeCashFlow = (*riskFreeCashFlows)[i].amount *
							discountFactor;

			feeLegPV += riskFreeCashFlow;
		}
	}
	return feeLegPV;
}

DRLIB_END_NAMESPACE
