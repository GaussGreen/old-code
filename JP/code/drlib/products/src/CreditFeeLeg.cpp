//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditFeeLeg.cpp
//
//   Description : Simple implementation of ICreditFeeLeg
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/CreditFeeLeg.hpp"
#include "edginc/CreditFeeLegSV.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FixedCashFlow.hpp"
#include "edginc/FloatCashFlow.hpp"
#include "edginc/AdhocCashFlow.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/IDecretionCurve.hpp"

DRLIB_BEGIN_NAMESPACE

CreditFeeLeg::~CreditFeeLeg() {}

CreditFeeLeg::CreditFeeLeg(): CreditFeeLegWithPV(TYPE), isFixed(false),
                              badDayConvention("NONE"){}

TweakOutcome CreditFeeLeg::sensShift(const PropertyTweak<CreditFeeLegCoupon>& shift) {
    for (int i = 0; i < spreads.size(); ++i) {
        if (accrualEndDates[i] >= valueDate) {
            spreads[i] += shift.coefficient;
        }
    }
    return TweakOutcome(shift.coefficient, false);
}

string CreditFeeLeg::sensName(const CreditFeeLegCoupon* shift) const
{
    return "CreditFeeLeg"; // please don't return "" since that means "don't tweak me"
}

void CreditFeeLeg::sensRestore(const PropertyTweak<CreditFeeLegCoupon>& shift){
    for (int i = 0; i < spreads.size(); ++i) {
        if (accrualEndDates[i] >= valueDate) {
            spreads[i] -= shift.coefficient;
        }
    }
}

CreditFeeLeg::CreditFeeLeg(
        double                    coupon,
        const DateTimeArray       &observationDates, // observation dates
        const DateTimeArray       &accrualStartDates, // accrual dates : can be non contiguous
        const DateTimeArray       &accrualEndDates,
        const DateTimeArray       &payDates,
        const CouponNotionalTypes &couponNotionalType,
        const string              &payDCC,
        const IModel* model, const MarketData* market) :
    CreditFeeLegWithPV(TYPE),
    isFixed(true),
    observationDates(observationDates),
    accrualStartDates(accrualStartDates),
    accrualEndDates(accrualEndDates),
    payDates(payDates),
    couponNotionalType(couponNotionalType),
    payDCC(payDCC),
    badDayConvention("NONE")
{
        for(int i=0; i<payDates.size(); i++) {
                notionals.push_back(1.0);
                fixings.push_back(coupon);
        }
        validatePop2Object();
        getMarket(model, market);
}


IObject* CreditFeeLeg::defaultConstructor() {
    return new CreditFeeLeg();
}


CClassConstSP const CreditFeeLeg::TYPE = CClass::registerClassLoadMethod(
    "CreditFeeLeg", typeid(CreditFeeLeg), load);

void CreditFeeLeg::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CreditFeeLeg, clazz);
    SUPERCLASS(CreditFeeLegWithPV);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(IRestorableWithRespectTo<CreditFeeLegCoupon>);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(isFixed,            "is the whole leg fixed ?");
    FIELD(notionals,          "notionals");
    FIELD(observationDates,   "observation dates");
    FIELD(accrualStartDates,  "accrual start dates");
    FIELD(accrualEndDates,    "accrual end dates");
    FIELD(payDates,           "payment dates");
    FIELD(payDCC,             "The pay day count convention");
    FIELD(fixings,            "fixings");
    FIELD(couponNotionalType, "coupon notional type: "
                                     "OBSERVATION_DATE "
                                     "or AVERAGE "
                                     "or I_TRAXX");
    FIELD       (upfrontFees,        "Upfront fees, if any");
    FIELD_MAKE_OPTIONAL(upfrontFees);

    // optional fields when isFixed = FALSE
    FIELD(refixDates,   "refix dates for IR leg");
    FIELD_MAKE_OPTIONAL(refixDates);
    FIELD(weights,      "weighes");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(refixInterval,"refix interval (e.g. 3M)");
    FIELD_MAKE_OPTIONAL(refixInterval);
    FIELD(rateDCC,      "The rate day count convention");
    FIELD_MAKE_OPTIONAL(rateDCC);
    FIELD(spreads,      "spreads");
    FIELD_MAKE_OPTIONAL(spreads);

    // other optional fields
    FIELD(badDayConvention,  "bad day convention for credit fee leg");
    FIELD_MAKE_OPTIONAL(badDayConvention);
    FIELD(couponCurve,  "Coupon curve");
    FIELD_MAKE_OPTIONAL(couponCurve);
    FIELD(valueDate,    "Value date");
    FIELD_MAKE_OPTIONAL(valueDate);

    // transient fields
    FIELD_NO_DESC(refixIntervalSP);
    FIELD_MAKE_TRANSIENT(refixIntervalSP);
    FIELD_NO_DESC(payDCCSP);
    FIELD_MAKE_TRANSIENT(payDCCSP);
    FIELD_NO_DESC(rateDCCSP);
    FIELD_MAKE_TRANSIENT(rateDCCSP);
    FIELD_NO_DESC(badDayConventionSP);
    FIELD_MAKE_TRANSIENT(badDayConventionSP);
}

double CreditFeeLeg::calculateRate(DateTime             refixDate) const {
    static const string method = "CreditFeeLeg::calculateRate";
    if (isFixed){
        throw ModelException(method, "Leg is fixed not floating");
    }
    //TODO : use the model
    return couponCurve->fwd(refixDate,
                            refixIntervalSP.get(),
                            badDayConventionSP.get(),
                            rateDCCSP.get(),
                            CompoundBasis::SIMPLE);
}

void CreditFeeLeg::setFixingforThetaShift(
    const DateTime&      valueDate,
    const YieldCurve*    discount,
    const DateTime&      rollDate){
    if (!isFixed){
        for (int i=0; i<refixDates.size(); i++)
        {
            if (rollDate.getDate() >= refixDates[i].getDate() &&
                valueDate.getDate() <= refixDates[i].getDate())
                if (!(valueDate == refixDates[i]) || Maths::isZero(fixings[i]))
                    fixings[i] = calculateRate(refixDates[i]);
        }
    }
}

/** populate from market cache */
void CreditFeeLeg::getMarket(const IModel* model, const MarketData* market)
{
    market->GetReferenceDate(valueDate);
    if (!isFixed) {
        couponCurve.getData(model, market);
        couponCurve->setProjectionCurve();
    }
}


/** Implementation of the Theta Shift interface */
bool CreditFeeLeg::sensShift(Theta* shift)
{
    static const string method = "CreditFeeLeg::sensShift(Theta)";
    try  {
        if (!isFixed && couponCurve.get()) {
            const DateTime& newDate = shift->rollDate(valueDate);

            setFixingforThetaShift(valueDate,
                                   couponCurve.get(),
                                   newDate);
            return true; // our components have theta type sensitivity
        }
        else {
            return false;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

static void checkSize(int arraySize, string arrayName, int expectedSize) {
    if (arraySize != expectedSize) {
        throw ModelException("checkArray",
            "Invalid size for " +
            arrayName +
            " : " +
            Format::toString(arraySize) +
            " (" +
            Format::toString(expectedSize) +
            " expected)");
    }
}


void CreditFeeLeg::validatePop2Object() {
    int i;
    static const string method("CreditFeeLeg::validatePop2Object");

    badDayConventionSP.reset(BadDayConventionFactory::make(badDayConvention));

    /* check coupon notional type */
    if (!isCouponNotionalTypeValid())
        throw ModelException(method, "Specified coupon notional type not supported");

    /* check array not empty */
    int notionalsSize = notionals.size();
    if (notionalsSize == 0) {
        throw ModelException(method, "Fee leg should not be empty");
    }

    /* check size of the arrays */
    checkSize(accrualStartDates.size(), "accrualStartDates", notionalsSize);
    checkSize(accrualEndDates.size(), "accrualEndDates", notionalsSize);
    checkSize(observationDates.size(), "observationDates", notionalsSize);
    checkSize(payDates.size(), "payDates", notionalsSize);

    /* check date ordering */
    for (i = 0; i < notionalsSize; i++) {
        if (accrualStartDates[i] >= accrualEndDates[i]) {
            throw ModelException(method, "accrual start date >= "
                                 "accrual end date for #"+
                                 Format::toString(i+1)+" periods");
        }
        if (observationDates[i] > payDates[i]) {
            throw ModelException(method, "observationDate> payDate for #"+
                                 Format::toString(i+1)+" periods");
        }
    }

    if (!isFixed) {
        /* check optional fields are filled */
        string missingFields = "";
        if (refixDates.empty()) {
            missingFields += "refixDates ";
        }
        if (weights.empty()) {
            missingFields += "weights ";
        }
        if (refixInterval.length() == 0) {
            missingFields += "refixInterval ";
        }
        if (rateDCC.length() == 0) {
            missingFields += "rateDCC ";
        }
        if (spreads.empty()) {
            missingFields += "spreads ";
        }
        if (missingFields.length() != 0) {
            throw ModelException(method, "Missing field(s) : " + missingFields);
        }

        /* check refix dates */
        for (i = 0; i < notionalsSize; i++) {
            if (refixDates[i] > payDates[i])
                throw ModelException(method,
                                     "refixDate > payDate for #" +
                                     Format::toString(i+1)+ " periods");
        }
    }

    /* fill transient fields */
    payDCCSP.reset(DayCountConventionFactory::make(payDCC));
    if (!isFixed) {
        rateDCCSP.reset(DayCountConventionFactory::make(rateDCC));
        refixIntervalSP.reset(new MaturityPeriod(refixInterval));
    }

    if (upfrontFees.get()) {
        // Verify that upfrontFees' dates are in increasing order
        CashFlow::ensureDatesIncreasing(*upfrontFees, "upfrontFees", false);
    }
}

bool CreditFeeLeg::isCouponNotionalTypeValid() const {
    return (couponNotionalType == OBSERVATION_DATE
        || couponNotionalType == AVERAGE);
}


/** Returns all cash flows */
AbstractCashFlowArrayConstSP CreditFeeLeg::getCashFlows(IForwardRatePricerSP model) const {
    DateTimeArraySP  upfrontDates;
    CDoubleArraySP   upfrontAmounts;
    int              numCashFlows;
    int              upfrontFeesSize = 0;
    int              payDatesSize = notionals.size();

    if (upfrontFees.get()) {
        // Allocate extra space for the upfront payment cashflow(s)
        upfrontFeesSize = upfrontFees->size();

        // Note that CashFlow::dates returns a DateTimeArray (no SP!)
        upfrontDates =
            DateTimeArraySP(new DateTimeArray(CashFlow::dates(*upfrontFees)));
        upfrontAmounts = CashFlow::amounts(*upfrontFees);
    }
    numCashFlows = payDatesSize + upfrontFeesSize;
    AbstractCashFlowArraySP cashFlows(new AbstractCashFlowArray(numCashFlows));

    // Add all cashflows to the "cashFlows" array
    if (isFixed) {
        FixedCashFlowSP fixedCashFlow;
        AdhocCashFlowSP adhocCashFlow;
        int upfrontIdx = 0;
        int payDatesIdx = 0;
        for(int i=0; i < numCashFlows; ++i) {
            while ((upfrontIdx < upfrontFeesSize) && // only true if upfrontFees!=null, AND
                   ((payDatesIdx >= payDatesSize) || // either no more paydates, or
                    (*upfrontDates)[upfrontIdx] < payDates[payDatesIdx])) // this date goes 1st
            {
                // The upfront payment's date is before the fixed cash flow's
                // date, so add it first
                adhocCashFlow = AdhocCashFlowSP(new AdhocCashFlow(
                    valueDate,
                    (*upfrontDates)[upfrontIdx],
                    (*upfrontAmounts)[upfrontIdx]));

                (*cashFlows)[i] = adhocCashFlow;

                ++i;
                ++upfrontIdx;
            }

            if (payDatesIdx < payDatesSize) {
                fixedCashFlow = FixedCashFlowSP(new FixedCashFlow(
                    valueDate,
                    payDates[payDatesIdx],
                    notionals[payDatesIdx],
                    accrualStartDates[payDatesIdx],
                    accrualEndDates[payDatesIdx],
                    payDCCSP,
                    fixings[payDatesIdx]));

                (*cashFlows)[i] = CreditCashFlowSP(new CreditCashFlow(
                    fixedCashFlow,
                    couponNotionalType,
                    DateTimeSP(new DateTime(observationDates[payDatesIdx]))));

                ++payDatesIdx;
            }
        }
    }
    else {
        FloatCashFlowSP floatCashFlow;
        AdhocCashFlowSP adhocCashFlow;
        int upfrontIdx = 0;
        int payDatesIdx = 0;
        for(int i=0; i < numCashFlows; ++i) {
            while ((upfrontIdx < upfrontFeesSize) && // only true if upfrontFees!=null, AND
                   ((payDatesIdx >= payDatesSize) || // either no more paydates, or
                    (*upfrontDates)[upfrontIdx] < payDates[payDatesIdx])) // this date goes 1st
            {
                // The upfront payment's date is before the floating cash flow's
                // date, so add it first
                adhocCashFlow = AdhocCashFlowSP(new AdhocCashFlow(
                    valueDate,
                    (*upfrontDates)[upfrontIdx],
                    (*upfrontAmounts)[upfrontIdx]));

                (*cashFlows)[i] = adhocCashFlow;

                ++i;
                ++upfrontIdx;
            }

            if (payDatesIdx < payDatesSize) {
                floatCashFlow = FloatCashFlowSP(new FloatCashFlow(
                    valueDate,
                    payDates[payDatesIdx],
                    notionals[payDatesIdx],
                    accrualStartDates[payDatesIdx],
                    accrualEndDates[payDatesIdx],
                    payDCCSP,
                    refixIntervalSP,
                    spreads[payDatesIdx],
                    false, //isCMS
                    refixDates[payDatesIdx],
                    rateDCCSP,
                    BadDayConventionSP(BadDayConventionFactory::clone(badDayConventionSP.get())),
                    couponCurve,
                    weights[payDatesIdx],
                    new double(fixings[payDatesIdx]),
                    0)); //applyAdjustments

                (*cashFlows)[i] = CreditCashFlowSP(new CreditCashFlow(
                    floatCashFlow,
                    couponNotionalType,
                    DateTimeSP(new DateTime(observationDates[payDatesIdx]))));

                ++payDatesIdx;
            }
        }
    }

    return cashFlows;
}

/** Return all cash flow dates */
DateTimeArraySP CreditFeeLeg::getCashFlowDates() const
{
    //make a copy of payment dates
    DateTimeArraySP dates = DateTimeArraySP(dynamic_cast<DateTimeArray*>(payDates.clone()));
    //append upfront dates
    DateTimeArraySP upfrontDates;
    if (upfrontFees.get()) {
        upfrontDates = 
            DateTimeArraySP(new DateTimeArray(CashFlow::dates(*upfrontFees)));
        for (int i=0; i<upfrontDates->size(); i++)
        {
            dates->push_back((*upfrontDates)[i]);
        }
    }
    return dates;   
}

/** Return risky cash flow dates */
DateTimeArraySP CreditFeeLeg::getRiskyCashFlowDates() const
{
    //all fees are risky in this form
    return getCashFlowDates();
}

/** Return riskfree cash flow dates */
DateTimeArraySP CreditFeeLeg::getRisklessCashFlowDates() const
{
    //no fees are riskless in this form
    return DateTimeArraySP();
}

/** Returns the accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP CreditFeeLeg::getAccrualPeriods() const
{
    int numAccrualPeriods = accrualStartDates.size();

    // create a smart pointer to an array of accrual periods, with enough
    // space for numAccrualPeriods
    AccrualPeriodArraySP accrualPeriods(
        new AccrualPeriodArray(numAccrualPeriods));

    // create each accrual period and store it in the array
    for (int i=0; i<numAccrualPeriods; ++i) {
        AccrualPeriodSP newAccrualPeriod =
            AccrualPeriodSP( new AccrualPeriod(accrualStartDates[i],
                                               accrualEndDates[i]));
        (*accrualPeriods)[i] = newAccrualPeriod;
    }
    return accrualPeriods;
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP CreditFeeLeg::getRiskyAccrualPeriods() const
{
    //all fees are risky in this form
    return getAccrualPeriods();
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP CreditFeeLeg::getRisklessAccrualPeriods() const
{
    //no fees are riskless in this form
    return AccrualPeriodArrayConstSP();
}

/** Returns risky coupon notional types */
CouponNotionalTypesArraySP CreditFeeLeg::getRiskyCouponNotionalTypes() const
{
	int size = accrualStartDates.size();
	if (size == 0)
		return CouponNotionalTypesArraySP();
	else
	{
		CouponNotionalTypesArraySP couponNotionalTypesArray(new CouponNotionalTypesArray());
		for (int i=0; i < size; ++i)
			couponNotionalTypesArray->push_back(couponNotionalType);
		return couponNotionalTypesArray;
	}
}

/** Returns risky observation dates */
DateTimeArraySP CreditFeeLeg::getRiskyObservationDates() const
{
	int size = observationDates.size();
	if (size == 0)
		return DateTimeArraySP();
	else
		return DateTimeArraySP(dynamic_cast<DateTimeArray*>(observationDates.clone()));
}


/** Return known cash flows corresponding to a CDO tranche */
CashFlowArraySP CreditFeeLeg::generateKnownCashFlows(
    const DateTime       today,
    const double         initialTrancheSize,
    const DateTimeArray  pastTrancheLossDates,
    const DoubleArray    pastTrancheLosses,
    const double         pastTrancheLoss,
    IForwardRatePricerSP model)
{
    static const string method("CreditFeeLeg::generateKnownCashflows");

    try {
        int    i;
        double trancheOutstanding;

        // get the basic fee leg cashflows (these assume full accrual on the cashflow)
        AbstractCashFlowArrayConstSP cfl = getCashFlows(model);
        int n = cfl->size();
        CashFlowArraySP fkcfl(new CashFlowArray(n));

        for (i=0; i<n; i++) {
            (*fkcfl)[i].date   = (*cfl)[i]->getPayDate();
            (*fkcfl)[i].amount = (*cfl)[i]->getAmount(model);

            //manipulate any credit type cashflows
            if (CreditCashFlow::TYPE->isInstance((*cfl)[i])) {
                //manipulate any historic flows
                if ((*fkcfl)[i].date < today) {
                    CreditCashFlow* ccf = dynamic_cast<CreditCashFlow*>(((*cfl)[i]).get());

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


/** Returns the pay date which is last in terms of time */
DateTime CreditFeeLeg::getLastPayDate() const {
    DateTime lastPay(payDates.back());
    for (int i = payDates.size()-2; i >= 0; i--){
        if (payDates[i] > lastPay){
            lastPay = payDates[i];
        }
    }

    // Search in the upfront payments too, just in case
    if (upfrontFees.get()) {
        for (int j=0; j<upfrontFees->size(); ++j) {
            if ((*upfrontFees)[j].date > lastPay) {
                lastPay = (*upfrontFees)[j].date;
            }
        }
    }
    return lastPay;
}


/** Returns the observation date which is last in terms of time */
DateTime CreditFeeLeg::getLastObservationDate() const {
    DateTime lastObs(observationDates.back());
    for (int i = observationDates.size()-2; i >= 0; i--){
        if (observationDates[i] > lastObs){
            lastObs = observationDates[i];
        }
    }
    return lastObs;
}


/** When to stop tweaking. Returns max(currentLastDate, this max date)
    where this max date = max(rate maturity, last pay date) */
DateTime CreditFeeLeg::lastYCSensDate(const DateTime& currentLastDate) const {
    //NB : payDates are NOT sorted
    DateTime thisMaxDate = getLastPayDate();
    if (!isFixed){
        // find last fixing
        int lastFixing = refixDates.size()-1;
        for (int i = refixDates.size() -2; i >= 0; i--){
            if (refixDates[i] > refixDates[lastFixing]){
                lastFixing = i;
            }
        }
        // find maturity of last fixing
        DateTime start(couponCurve->settles(refixDates[lastFixing]));
        thisMaxDate = thisMaxDate.max(
            couponCurve->rateMaturity(start, refixIntervalSP.get(),
                                      badDayConventionSP.get()));
    }
    return currentLastDate.max(thisMaxDate);
}

/** Feedback method for getting information about a fee cashflow */
void CreditFeeLeg::getActiveFee(const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
                                const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start dates, and we are interested in the first fee
                                IForwardRatePricerSP model,               // (I) for calculating the amount
                                DateTime&            accrueStartDate,     // (O) when the fee starts accruing
                                DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
                                DateTime&            paymentDate,         // (O) when the fee is paid
                                double&              amount) const        // (O) the cashflow amount
{
    int numFees = accrualStartDates.size();
    int idx;
    for (idx=0; idx<numFees; idx++)
    {
        if (withRespectTo.within(accrualStartDates[idx],
                                accrualEndDates[idx]))
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
        amount = 0.0;
    }
    else
    {
        //need to calculate the fee amounts
        AbstractCashFlowArrayConstSP cfls = getCashFlows(model);

        accrueStartDate = accrualStartDates[idx];
        accrueEndDate   = accrualEndDates[idx];
        paymentDate     = payDates[idx];
        amount          = (*cfls)[idx]->getAmount(model);
    }
}

/** Returns the leg notional */
double CreditFeeLeg::getFeeLegNotional() const
{
    static const string method = "CreditFeeLeg::getFeeLegNotional";
    double ntnl;

    //all notionals must be the same for this to succeed
    if (notionals.size() > 0)
    {
        ntnl = notionals[0];

        for (int i=1; i<notionals.size(); i++)
        {
            if (notionals[i] != ntnl)
            {
                throw ModelException(method,
                    "Leg has differing notionals across its structure");
            }
        }
    }
    else
    {
        //no components
        ntnl = 0.0;
    }

    return ntnl;
}

/** Returns the leg notional */
void CreditFeeLeg::setFeeLegNotional(double newNotional)
{
    for (int i=0; i<notionals.size(); i++)
    {
        notionals[i] = newNotional;
    }
}

//Only to handle fixed leg
ICreditFeeLegSVGenSP
CreditFeeLeg::createSVGen(
	const DateTime& valDate,
	const DateTimeArray& modifiedFeeLegObservationDates,
	const DateTimeArray& productTimeline,
	const IntArray& dateToDiscFactorIndex,
	double lossConfigNotional) const
{
//validation for floating leg
	if (!this->isFixed)
		throw ModelException(
						"CreditFeeLeg::createSVGen",
						"MC model cannot not handle floating rate cashflows");

//productTimeline can be null or size = 0
	AccrualPeriodArrayConstSP accrualPeriodArray = this->getRiskyAccrualPeriods();
	DayCountConventionArray dayCountConventions;
	IntArray dateToProductTimelineIndexMap;
	DoubleArray scalingFactors;
	DoubleArray historicalScalingFactors;

	CouponNotionalTypesArraySP couponNotionalTypes =
		this->getRiskyCouponNotionalTypes();

	DateTimeArraySP feePayDates =
		this->getRiskyCashFlowDates();

	if (accrualPeriodArray.get() != 0)
	{
		int size = accrualPeriodArray->size();
		if (size != fixings.size())
			throw ModelException(
							"CreditFeeLeg::createSVGen",
							"Number of accrual periods and fixings do not match");
//convert productTimeline to a int array
		ICreditFeeLeg::buildDateToProductTimelineIndexMap(
			dateToProductTimelineIndexMap, //output
			productTimeline); //input

//build DayCountConvention
		IObject* object = DayCountConvention::createFromString(
								DayCountConvention::TYPE,
								this->payDCC );

		DayCountConventionSP dcc(dynamic_cast<DayCountConvention*>(object));
		for (int i = 0; i< size; ++i)
			dayCountConventions.push_back(dcc);

//build scaling factors
		ICreditFeeLeg::buildScalingFactors(
			scalingFactors,
			historicalScalingFactors,
			valDate,
			*accrualPeriodArray,
			*feePayDates,
			*couponNotionalTypes,
			dayCountConventions,
			fixings,
			notionals,
			lossConfigNotional);
	}

//do the optimization - which is actually applicable to fixed leg pricer
	return ICreditFeeLegSVGenSP(
		new CreditFeeLegSVGen(
			productTimeline.size(),
			dateToProductTimelineIndexMap,
			dateToDiscFactorIndex,
			lossConfigNotional,
			scalingFactors,
			historicalScalingFactors,
			modifiedFeeLegObservationDates,
			*couponNotionalTypes,
			*accrualPeriodArray,
			fixings,
			this->payDates,
			this->upfrontFees
			));

}

//-------------------------------
// IFixedRateCreditFeeLeg methods
//-------------------------------
/**Get the fee rate of the fee leg.*/
double CreditFeeLeg::getRate() const
{
    static const string method = "CreditFeeLeg::getRate";

    if (isFixed)
    {
        double rate = 0.0;
        if (fixings.size() > 0)
        {
            rate = fixings[0];
            //validate the remaining fees
            for (int i=1; i<fixings.size(); i++)
            {
                if (fixings[i] != rate)
                {
                    throw ModelException(method,
                        "Leg does not have a constant fixed rate");
                }
            }
        }
        return rate;
    }
    else
    {
        throw ModelException(method,
            "Fee leg is not fixed");
    }
}

/**Change the fee rate of the fee leg. Note that, if the rate is currently zero,
    this will throw an exception.*/
void CreditFeeLeg::setRate(double newRate)
{
    static const string method = "CreditFeeLeg::setRate";
    if (isFixed)
    {
        for (int i=0; i<fixings.size(); i++)
        {
            fixings[i] = newRate;
        }
    }
    else
    {
        throw ModelException(method,
            "Fee leg is not fixed");
    }
}


DRLIB_END_NAMESPACE
