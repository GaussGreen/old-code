//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FullCreditFeeLeg.cpp
//
//   Description : Implementation of ICreditFeeLeg as a list of cashflows
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FullCreditFeeLeg.hpp"
#include "edginc/AdhocCashFlow.hpp"
#include "edginc/FixedCashFlow.hpp"
#include "edginc/FloatCashFlow.hpp"
#include "edginc/CreditFeeLegSV.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/ActualActual.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

FullCreditFeeLeg::~FullCreditFeeLeg() {}

FullCreditFeeLeg::FullCreditFeeLeg(): CreditFeeLegWithPV(TYPE)
{}

/** Explicit constructor */
FullCreditFeeLeg::FullCreditFeeLeg(AbstractCashFlowArraySP creditCashFlows) :
    CreditFeeLegWithPV(TYPE),
    creditCashFlows(creditCashFlows)
{
    validatePop2Object();
}

IObject* FullCreditFeeLeg::defaultConstructor() {
    return new FullCreditFeeLeg();
}

CClassConstSP const FullCreditFeeLeg::TYPE = CClass::registerClassLoadMethod(
    "FullCreditFeeLeg", typeid(FullCreditFeeLeg), load);

void FullCreditFeeLeg::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(FullCreditFeeLeg, clazz);
    SUPERCLASS(CreditFeeLegWithPV);
	EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(creditCashFlows, "Array of credit cash flows");
}

/** populate from market cache */
void FullCreditFeeLeg::getMarket(const IModel* model, const MarketData* market) {
    for (int i = 0; i < creditCashFlows->size(); ++i) {
       ((*creditCashFlows)[i])->getMarket(model, market);
    }
}

/** Called immediately after object constructed */
void FullCreditFeeLeg::validatePop2Object() {
    static const string method("FullCreditFeeLeg::validatePop2Object");
    try {
        if (creditCashFlows->size() < 1) {
            throw ModelException(method,
                "No cash flow in credit fee leg");
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Return all cash flow dates */
DateTimeArraySP FullCreditFeeLeg::getCashFlowDates() const
{
    //get pay dates from each of the cashflows
    int numFlows = creditCashFlows->size();
    DateTimeArraySP cfDates = DateTimeArraySP(
        new DateTimeArray(numFlows));

    for (int i=0; i<numFlows; i++)
    {
        (*cfDates)[i] = (*creditCashFlows)[i]->getPayDate();
    }

    return cfDates;
}

/** Return risky cash flow dates */
DateTimeArraySP FullCreditFeeLeg::getRiskyCashFlowDates() const
{
    //get pay dates from each of the cashflows
    int numFlows = creditCashFlows->size();
    DateTimeArraySP cfDates = DateTimeArraySP(
        new DateTimeArray(0));

    for (int i=0; i<numFlows; i++)
    {
        if (!((*creditCashFlows)[i]->isRiskFree()))
        {
            cfDates->push_back(((*creditCashFlows)[i])->getPayDate());
        }
    }

    return cfDates;
}

/** Return riskfree cash flow dates */
DateTimeArraySP FullCreditFeeLeg::getRisklessCashFlowDates() const
{
    //get pay dates from each of the cashflows
    int numFlows = creditCashFlows->size();
    DateTimeArraySP cfDates = DateTimeArraySP(
        new DateTimeArray(0));

    for (int i=0; i<numFlows; i++)
    {
        if ((*creditCashFlows)[i]->isRiskFree())
        {
            cfDates->push_back(((*creditCashFlows)[i])->getPayDate());
        }
    }

    return cfDates;
}

/** Returns the accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP FullCreditFeeLeg::getAccrualPeriods() const
{
    int numCreditCashFlows = creditCashFlows->size();

    // create a smart pointer to an array of accrual periods
    AccrualPeriodArraySP accrualPeriods(new AccrualPeriodArray(0));

    // create each accrual period and store it in the array
    for (int i=0; i<numCreditCashFlows; ++i) {
        AccrualPeriodSP newAccrualPeriod =
            (*creditCashFlows)[i]->getAccrualPeriod();

        if (!newAccrualPeriod->startDate().empty()) {
            // an empty accrual start date means this is a riskless
            // cashflow, so there is no associated accrual period
            accrualPeriods->push_back(newAccrualPeriod);
        }
    }
    return accrualPeriods;
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP FullCreditFeeLeg::getRiskyAccrualPeriods() const
{
    int numCreditCashFlows = creditCashFlows->size();

    // create a smart pointer to an array of accrual periods
    AccrualPeriodArraySP accrualPeriods(new AccrualPeriodArray(0));

    // create each accrual period and store it in the array
    for (int i=0; i<numCreditCashFlows; ++i)
    {
        if (!((*creditCashFlows)[i]->isRiskFree()))
        {
            AccrualPeriodSP newAccrualPeriod =
                (*creditCashFlows)[i]->getAccrualPeriod();

            accrualPeriods->push_back(newAccrualPeriod);
        }
    }
    return accrualPeriods;
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP FullCreditFeeLeg::getRisklessAccrualPeriods() const
{
    int numCreditCashFlows = creditCashFlows->size();

    // create a smart pointer to an array of accrual periods
    AccrualPeriodArraySP accrualPeriods(new AccrualPeriodArray(0));

    // create each accrual period and store it in the array
    for (int i=0; i<numCreditCashFlows; ++i)
    {
        if (((*creditCashFlows)[i]->isRiskFree()))
        {
            AccrualPeriodSP newAccrualPeriod =
                (*creditCashFlows)[i]->getAccrualPeriod();

            accrualPeriods->push_back(newAccrualPeriod);
        }
    }
    return accrualPeriods;
}

/** Returns risky coupon notional types */
CouponNotionalTypesArraySP FullCreditFeeLeg::getRiskyCouponNotionalTypes() const
{
	int numCreditCashFlows = creditCashFlows->size();
	CouponNotionalTypesArraySP couponNotionalTypesArray(new CouponNotionalTypesArray(0));

	for (int i=0; i<numCreditCashFlows; ++i)
	{
		try //if it is a creditCashFlow
		{
			CreditCashFlowSP creditCashFlow =
				CreditCashFlowSP::dynamicCast((*creditCashFlows)[i]);
			if (creditCashFlow.get())
				couponNotionalTypesArray->push_back(creditCashFlow->getCouponNotionalType());
		}
		catch(...)
		{
		}
	}
	return couponNotionalTypesArray;
}

/** Returns risky observation dates */
DateTimeArraySP FullCreditFeeLeg::getRiskyObservationDates() const
{
	int numCreditCashFlows = creditCashFlows->size();
	DateTimeArraySP riskObsDatesArray(new DateTimeArray());

	for (int i=0; i < numCreditCashFlows; ++i)
	{
		try //if it is a creditCashFlow
		{
			CreditCashFlowSP creditCashFlow =
				CreditCashFlowSP::dynamicCast((*creditCashFlows)[i]);
			if (creditCashFlow.get())
				riskObsDatesArray->push_back(creditCashFlow->getObservationDate());
		}
		catch(...)
		{
		}
	}
	return riskObsDatesArray;
}

/** Returns all cash flows */
AbstractCashFlowArrayConstSP FullCreditFeeLeg::getCashFlows(IForwardRatePricerSP model) const {
    return creditCashFlows;
}


/** Returns the pay date which is last in terms of time */
DateTime FullCreditFeeLeg::getLastPayDate() const {
    static const string method("FullCreditFeeLeg::getLastPayDate");

    DateTimeArraySP cfRiskyDates = getRiskyCashFlowDates();
    DateTimeArraySP cfRisklessDates = getRisklessCashFlowDates();

    // uses stl max_elements to find the last date (=max) of cfDates
    DateTimeArray::const_iterator maxRiskyDate = max_element(
        cfRiskyDates->begin(), cfRiskyDates->end());

    DateTimeArray::const_iterator maxRisklessDate = max_element(
        cfRisklessDates->begin(), cfRisklessDates->end());

    if (maxRiskyDate != cfRiskyDates->end()) {
        if (maxRisklessDate != cfRisklessDates->end()) {
            return max(*maxRiskyDate,*maxRisklessDate);
        }
        else {
            return *maxRiskyDate;
        }
    }
    else {
        if (maxRisklessDate != cfRisklessDates->end()) {
            return *maxRisklessDate;
        }
        else {
            // cfDates is empty
            throw ModelException(method,
                    "No cash flow in credit fee leg");
        }
    }
}


/** Returns the observation date which is last in terms of time */
DateTime FullCreditFeeLeg::getLastObservationDate() const {
    static const string method("FullCreditFeeLeg::getLastObservationDate");

    //use a closed form model for retrieving cashflows
    //really need to support "get...Dates" methods that do not require a model
    ClosedFormForwardRatePricerSP model = ClosedFormForwardRatePricerSP(
        new ClosedFormForwardRatePricer());

    DateTimeArraySP obsDates = getRiskyNotionalDates(model);

    // uses stl max_elements to find the last date (=max) of obsDates
    DateTimeArray::const_iterator maxDate = max_element(
        obsDates->begin(), obsDates->end());

    if (maxDate != obsDates->end()) {
        return *maxDate;
    } else {
        // obsDates is empty
        // return empty DateTime
        return DateTime();
    }
}

/** When to stop tweaking. Returns max(currentLastDate, this max date)
    where this max date = max(rate maturity, last pay date) */
DateTime FullCreditFeeLeg::lastYCSensDate(const DateTime& currentLastDate) const{
    static const string method("FullCreditFeeLeg::lastYCSensDate");

    DateTime thisMaxDate = currentLastDate;

    int n = creditCashFlows->size();
    for(int i=0; i<n; i++) {
        thisMaxDate = thisMaxDate.max((*creditCashFlows)[i]->lastYCSensDate());
    }

    return thisMaxDate;
}

/** Feedback method for getting information about a fee cashflow */
void FullCreditFeeLeg::getActiveFee(const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
                                    const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start dates, and we are interested in the first fee
                                    IForwardRatePricerSP model,               // (I) for calculating the amount
                                    DateTime&            accrueStartDate,     // (O) when the fee starts accruing
                                    DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
                                    DateTime&            paymentDate,         // (O) when the fee is paid
                                    double&              amount) const        // (O) the cashflow amount
{
    int numFees = creditCashFlows->size();
    int idx;
    for (idx=0; idx<numFees; idx++)
    {
        AccrualPeriodSP ap = (*creditCashFlows)[idx]->getAccrualPeriod();
        accrueStartDate = ap->startDate();
        accrueEndDate = ap->endDate();

        if (accrueStartDate.empty() && accrueEndDate.empty())
        {
            //use payment dates instead
            accrueStartDate = idx==0 ?
                              earliestAccrualDate :
                              (*creditCashFlows)[idx-1]->getPayDate();

            accrueEndDate = (*creditCashFlows)[idx]->getPayDate();
        }
        if (withRespectTo.within(accrueStartDate, accrueEndDate))
        {
            break;
        }
    }

    if (idx == numFees)
    {
        DateTime empty;
        accrueStartDate = empty;
        accrueEndDate   = empty;
        paymentDate     = empty;
        amount = 0;
    }
    else
    {
        paymentDate     = (*creditCashFlows)[idx]->getPayDate();
        amount          = (*creditCashFlows)[idx]->getAmount(model);
    }
}

/** Return known cash flows */
CashFlowArraySP FullCreditFeeLeg::generateKnownCashFlows(
    const DateTime       today,
    const double         initialTrancheSize,
    const DateTimeArray  pastTrancheLossDates,
    const DoubleArray    pastTrancheLosses,
    const double         pastTrancheLoss,
    IForwardRatePricerSP model)
{
    static const string method("FullCreditFeeLeg::generateKnownCashflows");

    try {
        int    i;
        double trancheOutstanding;

        int n = creditCashFlows->size();
        CashFlowArraySP fkcfl(new CashFlowArray(n));

        //manipulate any historic flows
        for (i=0; i<n; i++) {
            (*fkcfl)[i].date   = (*creditCashFlows)[i]->getPayDate();
            (*fkcfl)[i].amount = (*creditCashFlows)[i]->getAmount(model);

            //manipulate any credit type cashflows
            if (CreditCashFlow::TYPE->isInstance((*creditCashFlows)[i])) {
                if ((*fkcfl)[i].date < today) {
                    CreditCashFlow* ccf =
                        dynamic_cast<CreditCashFlow*>(((*creditCashFlows)[i]).get());

                    CashFlow cf = ccf->getKnownCashFlow(initialTrancheSize,
                                                        pastTrancheLossDates,
                                                        pastTrancheLosses,
                                                        model);

                    (*fkcfl)[i].amount = cf.amount;
                }
                else {
                    //future cashflow
                    trancheOutstanding = 1 - (pastTrancheLoss / initialTrancheSize);
                    //scale cashflow
                    (*fkcfl)[i].amount *= trancheOutstanding;
                }
            }
        }
        return fkcfl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns the leg notional */
double FullCreditFeeLeg::getFeeLegNotional() const
{
    static const string method("FullCreditFeeLeg::getFeeLegNotional");
    double ntnl;

    int n = creditCashFlows->size();

    if (n == 0)
    {
        //no components, just return 0
        ntnl = 0.0;
    }
    else
    {
        ntnl = (*creditCashFlows)[0]->getNotional();

        //validate consistency of other notionals
        for (int i=1; i<n; i++)
        {
            if ((*creditCashFlows)[i]->getNotional() != ntnl)
            {
                throw ModelException(method,
                    "Leg has differing notionals across its structure");
            }
        }
    }

    return ntnl;
}

/** Returns the leg notional */
void FullCreditFeeLeg::setFeeLegNotional(double newNotional)
{
    int n = creditCashFlows->size();
    for (int i=0; i<n; i++)
    {
        (*creditCashFlows)[i]->setNotional(newNotional);
    }
}


ICreditFeeLegSVGenSP
FullCreditFeeLeg::createSVGen(
	const DateTime& valDate,
	const DateTimeArray& modifiedFeeLegObservationDates,
	const DateTimeArray& productTimeline, //product timeline
	const IntArray& dateToDiscFactorIndex,
	double lossConfigNotional
) const
{
//risky cashflow parameters
	AccrualPeriodArray riskyAccrualPeriods;
	CouponNotionalTypesArray couponNotionalTypes;
	DateTimeArray observationDates;
	DateTimeArray payDates;
	DayCountConventionArray dayCountConventions;
	DoubleArray coupons;
	DoubleArray notionals;
	IntArray dateToProductTimelineIndexMap;
	DoubleArray scalingFactors;
	DoubleArray historicalScalingFactors;

//riskfree casflow parameters
	CashFlowArraySP riskFreeCashFlows(new CashFlowArray());

	if (this->creditCashFlows.get())
	{
		for (int i=0; i< creditCashFlows->size(); ++i)
		{
			if ( (*creditCashFlows)[i]->isRiskFree())
			{
//we don't handle riskfree floating cash flows yet
				if (FloatCashFlow::TYPE->isInstance((*creditCashFlows)[i]))
					throw ModelException("FullCreditFeeLeg::createSVGen",
										"MC model does not cater for riskfree floating leg cashflows");

//we convert fixed cash flow into an adhoc cash flow for ease of MC computation
				if (FixedCashFlow::TYPE->isInstance((*creditCashFlows)[i]))
				{
					FixedCashFlowConstSP fixedCashFlow =
						FixedCashFlowConstSP::dynamicCast((*creditCashFlows)[i]);
					double amount = fixedCashFlow->getAmount(
														fixedCashFlow->getRate());
					riskFreeCashFlows->push_back(CashFlow(
													fixedCashFlow->getPayDate(),
													amount));
				}
				else //it has to be an adhoc cashflow
				{
					AdhocCashFlowConstSP adhocCashFlow =
						AdhocCashFlowConstSP::dynamicCast((*creditCashFlows)[i]);
					riskFreeCashFlows->push_back(CashFlow(
													adhocCashFlow->getPayDate(),
													adhocCashFlow->getAmount(0)));
				}
			}
			else  //it has to be CreditCashFlow
			{
//too many dynamic casts - needs to be redesigned with each kind of cashflow creating a corresponding SV
				CreditCashFlowConstSP creditCashFlow =
					CreditCashFlowConstSP::dynamicCast((*creditCashFlows)[i]);

//collect observation date related parameters
				couponNotionalTypes.push_back(creditCashFlow->getCouponNotionalType());
				observationDates.push_back(creditCashFlow->getObservationDate());

//collect risky pay dates
				payDates.push_back(creditCashFlow->getPayDate());

//we don't handle risky floating cashflows as yet
				RiskFreeCashFlowConstSP riskFreeCash =
					creditCashFlow->getRiskFreeCashFlow();
				if (FloatCashFlow::TYPE->isInstance(riskFreeCash))
					throw ModelException("FullCreditFeeLeg::createSVGen",
										"MC model does not cater for risky floating leg cashflows");

				if (AdhocCashFlow::TYPE->isInstance(riskFreeCash))
				{
//here we are converting the adhoc cashflow into into a coupon paying cashflow
					AdhocCashFlowConstSP adhocCashFlow =
						AdhocCashFlowConstSP::dynamicCast(riskFreeCash);

//collect all risky accruals
					DateTime payDate = adhocCashFlow->getMaturityDate();
					riskyAccrualPeriods.push_back(AccrualPeriodSP(new AccrualPeriod(valDate, payDate )));

//set the day count convention to ActualActual
					dayCountConventions.push_back(DayCountConventionSP(new ActualActual()));
					double yearFraction = dayCountConventions[i]->years(
																valDate,
																payDate);
					coupons.push_back(1.0/yearFraction);

					notionals.push_back(adhocCashFlow->getAmount(0));
				}
				else
				{
//collect all risky accruals
					AccrualPeriodSP accrualPeriod =
						creditCashFlow->getAccrualPeriod();
					riskyAccrualPeriods.push_back(accrualPeriod);

//here it is guarantedly a fixed leg cashflows
					FixedCashFlowConstSP fixedCashFlow =
						FixedCashFlowConstSP::dynamicCast(riskFreeCash);

//fetch the day count convention
					DayCountConvention* dcc =
						const_cast<DayCountConvention*>(fixedCashFlow->getAccrualDcc().get());
					dayCountConventions.push_back(DayCountConventionSP(dcc));

//fetch the fixed rate coupon
					coupons.push_back(fixedCashFlow->getRate());

//fetch the notionals
					notionals.push_back(fixedCashFlow->getNotional());
				}
			}
		}

//convert productTimeline to a int array
		ICreditFeeLeg::buildDateToProductTimelineIndexMap(
			dateToProductTimelineIndexMap, //output
			productTimeline); //input

//build risky scaling factors
		ICreditFeeLeg::buildScalingFactors(
			scalingFactors,
			historicalScalingFactors,
			valDate,
			riskyAccrualPeriods,
			payDates,
			couponNotionalTypes,
			dayCountConventions,
			coupons,
			notionals,
			lossConfigNotional);
	}
	return ICreditFeeLegSVGenSP(
		new CreditFeeLegSVGen(
			productTimeline.size(),
			dateToProductTimelineIndexMap,
			dateToDiscFactorIndex,
			lossConfigNotional,
			scalingFactors,
			historicalScalingFactors,
			modifiedFeeLegObservationDates,
			couponNotionalTypes,
			riskyAccrualPeriods,
			coupons,
			payDates,
			riskFreeCashFlows
			));
}

//-------------------------------
// IFixedRateCreditFeeLeg methods
//-------------------------------

/**Get the fee rate of the fee leg.*/
double FullCreditFeeLeg::getRate() const
{
    static const string method = "CreditFeeLeg::getRate";

    double rate = 0.0;
    //rely on the abstract cashflows to sanity check
    if (creditCashFlows->size() > 0)
    {
        rate = (*creditCashFlows)[0]->getRate();
        //validate consistency across remaining fees
        for (int i=1; i<creditCashFlows->size(); i++)
        {
            if ((*creditCashFlows)[i]->getRate() != rate)
            {
                throw ModelException(method,
                    "Leg does not have a constant fixed rate");
            }
        }
    }

    return rate;
}

/**Change the fee rate of the fee leg. Note that, if the rate is currently zero,
    this will throw an exception.*/
void FullCreditFeeLeg::setRate(double newRate)
{
    for (int i=0; i<creditCashFlows->size(); i++)
    {
        (*creditCashFlows)[i]->setRate(newRate);
    }
}

DRLIB_END_NAMESPACE
