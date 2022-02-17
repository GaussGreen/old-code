//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditFeeLegWithPV.cpp
//
//   Description : ICreditFeeLeg able to pv itself
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CreditFeeLeg.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/CreditFeeLegWithPV.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


/* Default number of days between calculation date and fee rebate payment date
 * when it is not provided (which may happen when there is no contingent leg in
 * the trade since this data typically comes from the contingent leg). */
#define DEFAULT_DELAY_DAYS 1


/** TYPE for CreditFeeLeg */
CClassConstSP const CreditFeeLegWithPV::TYPE =
    CClass::registerClassLoadMethod(
        "CreditFeeLegWithPV", typeid(CreditFeeLegWithPV), load);

void CreditFeeLegWithPV::load(CClassSP& clazz){
    clazz->setPrivate();
    REGISTER(CreditFeeLegWithPV, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICreditFeeLeg);
}

CreditFeeLegWithPV::CreditFeeLegWithPV(CClassConstSP clazz): CObject(clazz) {}

/**
 * Calculate fee leg price. Unsure about what happens for accrued interest
   for historic defaults. Certainly debugUnitPrice is always zero.
*/
double CreditFeeLegWithPV::price(
    const DateTime&             today,
    const DateTime&             valDateCF,
    const IDiscountCurveRiskySP effectiveCurve,
    const YieldCurveWrapper&    discount,
    const double&               lowStrike,
    const double&               highStrike,
    const double&               outstandingNotional,
    const CashFlowArray&        pastTrancheLosses,
    double&                     riskyDurationTotal, // (O) notional weighted risky duration
    double&                     riskyNotionalsMean, // (O) mean value of notional per period
    double&                     risklessCFPV,       // (O) fair value of riskless payments
    bool                        computeExtra,       // true: populate arrays
    DoubleArray&                debugUnitPrice,     // price for each leg unit
    DoubleArray&                debugUnitHistPrice, // price for each leg unit due to historical default
    CashFlowArraySP             rebatePayments,     
    BoolArrayConstSP            payAsYouGoArray,    
    IntArrayConstSP             numDelayDaysArray,  
    DateTimeArrayConstSP        startDates,         
    DateTimeArrayConstSP        endDates,           
    DateTimeArrayConstSP        paymentDates,
    IBadDayAdjusterConstSP      bda,
    IForwardRatePricerSP        model) const
{
    static const string method("CreditFeeLegWithPV::price");
    try {
        AbstractCashFlowArrayConstSP cashFlows = getCashFlows(model);
        int nbCashFlows = cashFlows->size();
        if (computeExtra){
            debugUnitPrice.resize(nbCashFlows);
            debugUnitHistPrice.resize(nbCashFlows);
        }
        double total = 0.0; // fair value of all periods
        risklessCFPV = 0.0; // fair value of riskless fees
        riskyDurationTotal = 0.0; // notional weighted risky duration
        riskyNotionalsMean = 0.0; // mean value of notional per risky period
        double riskyNotionalsTotal = 0.0; // sum of notional of risky periods
        int numRiskyPeriods = 0; // how many risky periods
        double initialNotional = highStrike - lowStrike;
        double initialNotionalFeeLeg = highStrike - Maths::max(lowStrike, 0.0);

        //loop over cashflows
        for (int i = 0; i < nbCashFlows; i++) {
            const DateTime& cfDate = (*cashFlows)[i]->getPayDate();
         
            if (computeExtra){
                debugUnitHistPrice[i] = 0.0; // what about accrued?
                debugUnitPrice[i] = 0.0; 
            }
            /* skip historical payment w.r.t. value date */
            if (cfDate <= valDateCF){
                continue;
            }
            
            /* discount factor */
            double disc = effectiveCurve->risklessPV(cfDate);

            double outstandingFraction = 1.0;
            if (CreditCashFlow::TYPE->isInstance((*cashFlows)[i])) {
                CreditCashFlowConstSP creditCF = 
                    CreditCashFlowConstSP::dynamicCast((*cashFlows)[i]);

                outstandingFraction = 
                    creditCF->expectedTrancheOutstanding(
                        initialNotionalFeeLeg,
                        outstandingNotional,
                        today,
                        pastTrancheLosses,
                        effectiveCurve,
                        IDiscountCurveRisky::RECOVER_1,
                        0.) / initialNotional;
            }
                    
            // update price
            double cfAmount = (*cashFlows)[i]->getAmount(model);
            double unitPrice = disc * cfAmount * outstandingFraction;
            if (computeExtra){
                debugUnitPrice[i] = unitPrice;
            }
            total += unitPrice; // add to running total
            if (!(*cashFlows)[i]->isAdhoc()) {
                // update risky duration
                riskyDurationTotal += disc * 
                    (*cashFlows)[i]->getAmount(1.0) * outstandingFraction;
                // update sum of fee leg notionals
                riskyNotionalsTotal += (*cashFlows)[i]->getNotional();
                numRiskyPeriods++;
            } else {
                risklessCFPV += unitPrice;
            }
        }
        if (numRiskyPeriods != 0 && riskyNotionalsTotal != 0.0){
            riskyNotionalsMean = riskyNotionalsTotal/numRiskyPeriods;
        }


        // Take the rebate payments into account
        double rebatePrice = 
            priceRebatePayments(effectiveCurve, valDateCF, rebatePayments, 
                                payAsYouGoArray, numDelayDaysArray, startDates, 
                                endDates, paymentDates, bda); 

        total += rebatePrice;

        return total;
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Prices the rebate cashflows by computing the pay date of each of them,
 * according to the "pay as you go" configuration. 
 * If payAsYouGoArray, numDelayDaysArray, startDates, endDates and paymentDates
 * are all null, default values will be used - but if only some of them are null
 * an exception will be thrown */
double CreditFeeLegWithPV::priceRebatePayments(
    const IDiscountCurveRiskySP effectiveCurve,
    const DateTime&             valDateCF,
    CashFlowArraySP             rebatePayments,    
    BoolArrayConstSP            payAsYouGoArray,   
    IntArrayConstSP             numDelayDaysArray, 
    DateTimeArrayConstSP        startDates,        
    DateTimeArrayConstSP        endDates,          
    DateTimeArrayConstSP        paymentDates,
    IBadDayAdjusterConstSP      bda) const
{
    static const string method("CreditFeeLegWithPV::priceRebatePayments");

    int numRebates = !rebatePayments ? 0 : rebatePayments->size();
    if (numRebates == 0) {
        return 0.0; // No rebates
    }
    
    // Need to validate the parameters first
    // If all parameters are null, need to use defaults
    if (!payAsYouGoArray && !numDelayDaysArray && !startDates && 
        !endDates && !paymentDates)
    {
        payAsYouGoArray.reset(new BoolArray(1, true));
        numDelayDaysArray.reset(new IntArray(1, DEFAULT_DELAY_DAYS));
        startDates.reset(new DateTimeArray(1, DateTime()));
        // Use the maximum rebate date as the observation period end date
        DateTime maxDate = (*rebatePayments)[0].date;
        for (int i=1; i < numRebates; ++i) {
            maxDate = maxDate.max((*rebatePayments)[i].date);
        }
        maxDate = bda->addBusinessDays(maxDate, DEFAULT_DELAY_DAYS+1);
        endDates.reset(new DateTimeArray(1, maxDate));
        paymentDates.reset(new DateTimeArray(1, maxDate));       
    }
    else if (!payAsYouGoArray || !numDelayDaysArray || !startDates || 
        !endDates || !paymentDates) 
    {
        throw ModelException(method, 
                             "Internal error: if computing rebate "
                             "payments the payAsYouGo, numDelayDaysArray, "
                             "startDates, endDates and paymentDates arrays "
                             "must be either all null or all non-null");
    }

    int numPeriods = payAsYouGoArray->size();
    if ((numPeriods != numDelayDaysArray->size()) || 
        (numPeriods != startDates->size()) || 
        (numPeriods != endDates->size()) ||
        (numPeriods != paymentDates->size()))
    {
        throw ModelException(method, 
                             "Internal error: if computing rebate "
                             "payments the payAsYouGo, numDelayDaysArray, "
                             "startDates, endDates and paymentDates arrays "
                             "must have the same size.");
    }
    if (!numPeriods) {
        throw ModelException(method, 
                             "Internal error: if computing rebate "
                             "payments the payAsYouGo, startDates and "
                             "endDates arrays must have non-zero size.");
    }

    double total = 0.0;
    for (int i=0; i < numRebates; ++i) {
        const DateTime& rebateDate = (*rebatePayments)[i].date;

		// Find the period this rebate payment falls on
        int periodIdx = 0;
        for (; periodIdx < numPeriods; ++periodIdx) {
            /* In case there are rebates happening between protection periods, 
             * use the information for a period if the rebate falls after the 
             * period's start date or after the previous period's end date, 
             * whichever is earliest */
            const DateTime& minDate = (periodIdx == 0) ?
                (*startDates)[periodIdx] : 
                (*startDates)[periodIdx].min((*endDates)[periodIdx-1]);
            
            if ((rebateDate >= minDate) && 
                (rebateDate <= (*endDates)[periodIdx]))
            {
                break;
            }
        }

        // If periodIdx goes beyond the last period see if using the previous
        // period is acceptable
        if (periodIdx == numPeriods) {
            --periodIdx; // 
            if ((*payAsYouGoArray)[periodIdx]) {
                // The last period is pay as you go, so will use its
                // information
            }
            else {
                // The last period is pay at the end. As long as the rebate 
                // payment falls before the paydate, will use that information
                if (rebateDate > (*paymentDates)[periodIdx]) {
                    throw ModelException(method,
                                         "There is a rebate payment with "
                                         "calculation date on " +
                                         rebateDate.toString() +
                                         " but the last protection period has "
                                         "a pay date of " +
                                         (*paymentDates)[periodIdx].toString() +
                                         " - the pay date of the last period "
                                         "needs to be updated to cover this "
                                         "date.");
                }
            }
        }

        // We are interested in the (updated) observation period 'periodIdx'
        DateTime payDate;
        if ((*payAsYouGoArray)[periodIdx]) {
            int delay = (*numDelayDaysArray)[periodIdx];
            payDate = bda->addBusinessDays(rebateDate, delay);
        }
        else {
            payDate = (*paymentDates)[periodIdx];
        }

        // skip historical payment w.r.t. value date
        if (payDate <= valDateCF) {
            continue;
        }
                
        // PV this rebate payment, and add the result to "total"
        double disc = effectiveCurve->risklessPV(payDate);
        total -= disc * (*rebatePayments)[i].amount;
    }
    return total;
}


/** values all cashflows after cfCutOffDate and before or on endDate*/
/** value of risky cashflows only */
double CreditFeeLegWithPV::pvRisky (const DateTimeArray& timeline, 
                                    const DoubleArray&   effCurve, 
                                    const DateTime&      cfCutOffDate,
                                    const DateTime&      endDate,
                                    YieldCurveConstSP    discCurve,
                                    IForwardRatePricerSP model)
{
	const string method("CreditFeeLegWithPV::pvRisky");
	try {
		int i; //counter
        /** get cash flows	*/      
        AbstractCashFlowArrayConstSP cashFlows = getCashFlows(model);
        int nbCashFlows = cashFlows->size();
        
		double cfAmount = 0;
		double disc;
		double outstandingFraction=0;
		double value = 0;
		DateTime cfDate;
        CreditCashFlowConstSP creditCF;
        
        /**compute outstanding fraction and price for risky cashflows */
        for (i = 0; i < nbCashFlows; ++i) {
           if (CreditCashFlow::TYPE->isInstance((*cashFlows)[i])) {
                CreditCashFlowConstSP creditCF = 
                    CreditCashFlowConstSP::dynamicCast((*cashFlows)[i]);
               
                cfAmount = creditCF->getAmount(model);
                cfDate = creditCF->getPayDate();
                
                /* skip payments on and before cut off date */
                if (cfDate <= cfCutOffDate) continue;
    			/** skip payments after endDate if endDate is given*/
				if(!endDate.empty())
				{
					if (cfDate > endDate) break;
				}
    
                /* discount factor */
                disc = discCurve->pv(cfDate);
                
    			/** need to interpolate effective curve */
                outstandingFraction = CCMPriceUtil::linearLossInterpolate(
                    creditCF->getObservationDate(),
                    timeline,
                    effCurve);
                        
                // update pv
                value += disc * cfAmount * outstandingFraction;
           }
        }
		return value;
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

// overloaded version of pvRisky without end date
double CreditFeeLegWithPV::pvRisky (const DateTimeArray& timeline, 
                                    const DoubleArray&   effCurve, 
                                    const DateTime&      cfCutOffDate,
                                    YieldCurveConstSP    discCurve,
                                    IForwardRatePricerSP model)
{
	return pvRisky(
		timeline,
		effCurve,
		cfCutOffDate,
		DateTime(),
		discCurve,
        model);
}


/** Returns the risky or riskless abstract cash flows, depending on the 
 * value of the "wantRisky" parameter */
AbstractCashFlowArrayConstSP 
CreditFeeLegWithPV::getAbstractCashFlows(bool wantRisky,
                                         IForwardRatePricerSP model) const 
{
    static const string method("CreditFeeLegWithPV::getAbstractCashFlows");
    try {
        AbstractCashFlowArrayConstSP allCashFlows = getCashFlows(model);
        AbstractCashFlowArraySP cashFlows(new AbstractCashFlowArray(0));
        int n = allCashFlows->size();
        
        for (int i=0; i<n; ++i) {
            if (CreditCashFlow::TYPE->isInstance((*allCashFlows)[i]) == wantRisky) {
                // Here isInstance() returns true if the cashflow is risky, so
                // we are picking those cashflows with the same risky status as
                // "wantRisky"
                cashFlows->push_back((*allCashFlows)[i]);
            }
        }
        return cashFlows;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Returns risky cash flows only */
CashFlowArraySP CreditFeeLegWithPV::getRiskyCashFlows(IForwardRatePricerSP model) const {
    static const string method("CreditFeeLegWithPV::getRiskyCashFlows");
    try {
        AbstractCashFlowArrayConstSP creditCashFlows = 
            getAbstractCashFlows(true /*want the risky cashflows*/,
                                 model);

        int n = creditCashFlows->size();
        CashFlowArraySP cashFlows(new CashFlowArray(0));
        cashFlows->reserve(n);
        
        for(int i=0; i<n; ++i) {
            // Create a CashFlow from each (risky) creditCashFlow
            CashFlow cf((*creditCashFlows)[i]->getPayDate(),
                        (*creditCashFlows)[i]->getAmount(model));
            cashFlows->push_back(cf);
        }
        return cashFlows;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Returns risk free cash flows only */
CashFlowArraySP CreditFeeLegWithPV::getRisklessCashFlows(IForwardRatePricerSP model) const {
    static const string method("CreditFeeLegWithPV::getRisklessCashFlows");
    try {
        AbstractCashFlowArrayConstSP creditCashFlows = 
            getAbstractCashFlows(false /*want the riskless cashflows*/,
                                 model);

        int n = creditCashFlows->size();
        CashFlowArraySP cashFlows(new CashFlowArray(0));
        cashFlows->reserve(n);

        
        for(int i=0; i<n; ++i) {
            // Create a CashFlow from each (riskless) creditCashFlow
            CashFlow cf((*creditCashFlows)[i]->getPayDate(),
                        (*creditCashFlows)[i]->getAmount(model));
            cashFlows->push_back(cf);
        }
        return cashFlows;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Returns risky notional dates */
DateTimeArraySP CreditFeeLegWithPV::getRiskyNotionalDates(IForwardRatePricerSP model) const {
    static const string method("CreditFeeLegWithPV::getNotionalDates");
    try {
        AbstractCashFlowArrayConstSP creditCashFlows = getCashFlows(model);
        int n = creditCashFlows->size();
        DateTimeArraySP notionalDates(new DateTimeArray(0));
        
        for(int i=0; i<n; i++) {
            if (CreditCashFlow::TYPE->isInstance((*creditCashFlows)[i]))
            {
                CreditCashFlow* ccf = dynamic_cast<CreditCashFlow*>(((*creditCashFlows)[i]).get());
                notionalDates->push_back(ccf->getObservationDate());
            }
        }
        return notionalDates;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Return known cash flows corresponding to a non-defaulted CDS.
 * CAUTION: 
 * Cashflow amounts are obtained by calling getAmount(). 
 * If the cashflow is not an adhoc cashflow, the amounts are internally 
 * computed by multiplying the coupon by the notional and therefore 
 * both will typically have the same sign (negative coupons are allowed
 * though). For CDSs this is wrong (notional > 0 means long protection, 
 * so fees should typically be negative). 
 * Therefore non-adhoc cashflow amounts will by multiplied 
 * by -1 here, even though this is not consistent with other methods 
 * in this class!  */
CashFlowArraySP CreditFeeLegWithPV::generateKnownCashFlows(IForwardRatePricerSP model) const {
    static const string method = 
        "CreditFeeLegWithPV::generateKnownCashFlows";

    try {
        double amount;
        DateTime payDate;

        // Consider ALL payments, risky and riskless
        AbstractCashFlowArrayConstSP cashFlows = getCashFlows(model);
        int numCashFlows = cashFlows->size();
        
        // Store the cashflows' information
        CashFlowArraySP kcfl(new CashFlowArray(0));
        kcfl->reserve(numCashFlows);
        for (int idx=0; idx < numCashFlows; ++idx)
        {
            payDate = (*cashFlows)[idx]->getPayDate();
            amount = (*cashFlows)[idx]->getAmount(model);

            kcfl->push_back(CashFlow(payDate, amount));
        }
        return kcfl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns the known cash flows corresponding to a defaulted CDS
 * (taking accrued payments into consideration), that happen after a specific
 * date (typically used to exclude cashflows in the past, already paid).
 * If excludePaymentsBeforeDate is empty, all cashflows will be returned.
 *
 * CAUTION:
 * It does not aggregate cashflows happening on the same dates into one.
 * If this is required (e.g., to output these cashflows for MiddleOffice use)
 * you should call CashFlow::agregate(...) on the resulting CashFlowArray -
 * It is not done here for performance reasons */
CashFlowArraySP CreditFeeLegWithPV::generateKnownCashFlowsGivenDefault(
    const DateTime&            valuationDate,
    const DateTime&            defaultDate,
    const DateTime&            excludePaymentsBeforeDate,
    const DateTime&            protectionStartDate,
    const DateTime&            protectionEndDate,
    const bool                 allowIncludingTodaysPayments,
    IForwardRatePricerSP       model,
    ICreditEventOverrideNameSP creditEventOverride,
    IBadDayAdjusterConstSP     badDayAdjuster,
    CIntSP                     triggerDelay,
    CIntSP                     defaultToSettlementDelay,
    const DateTime&            lastTriggerDate) const 
{
    static const string method(
        "CreditFeeLegWithPV::generateKnownCashFlowsGivenDefault");

    try {
        // Consider ALL payments, risky and riskless
        AbstractCashFlowArrayConstSP cashFlows = getCashFlows(model);
        int numCashFlows = cashFlows->size();

        if (numCashFlows == 0) {
            return CashFlowArraySP();
        }
        
        DateTime determinationDate;
        DateTime accrualPayDate;
        getDeterminationAndPayDate(protectionStartDate,
                                   protectionEndDate,
                                   valuationDate,
                                   defaultDate,
                                   creditEventOverride,
                                   lastTriggerDate,
                                   triggerDelay,
                                   defaultToSettlementDelay,
                                   badDayAdjuster,
                                   determinationDate, // (O)
                                   accrualPayDate);   // (O)

        // Flag to indicate whether payments on valuationDate should be included in
        // the valuation or not. In general they should not but for compatibility
        // with the old CredDefSwap implementation, if the triggerDelay and 
        // defaultToSettlementDelay are NOT passed in they will be included 
        bool includeTodaysPayments = !creditEventOverride && 
                                     !triggerDelay &&
                                     allowIncludingTodaysPayments;

        CashFlowArraySP kcfl(new CashFlowArray(0));
        for (int idx=0; idx < numCashFlows; ++idx) {
            AbstractCashFlowSP cf = (*cashFlows)[idx]; // For ease

            const CashFlow& newCashFlow = cf->getKnownCashFlowGivenDefault(
                determinationDate,
                accrualPayDate,
                model);
            
            // only consider it if <amount!=0> and < <payment happens after
            // excludePaymentsBeforeDate> or 
            // < <newCashFlow.date=excludePaymentsBeforeDate> and 
            // <excludePaymentsBeforeDate=valuationDate> and includeTodaysPayments> >
            if (!Maths::isZero(newCashFlow.amount) && 
                ((newCashFlow.date > excludePaymentsBeforeDate) ||
                 ((newCashFlow.date == excludePaymentsBeforeDate) &&
                  (excludePaymentsBeforeDate == valuationDate) && 
                  includeTodaysPayments)))
            {
                kcfl->push_back(newCashFlow);
            }
        }
        return kcfl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Returns the determination and accrual pay date corresponding
    to a defaulted CDS, using the event override if required */
void CreditFeeLegWithPV::getDeterminationAndPayDate(
    const DateTime&            protectionStartDate,
    const DateTime&            protectionEndDate,
    const DateTime&            valuationDate,
    const DateTime&            defaultDate,
    ICreditEventOverrideNameSP creditEventOverride,
    const DateTime&            lastTriggerDate,
    CIntSP                     triggerDelay,
    CIntSP                     defaultToSettlementDelay,
    IBadDayAdjusterConstSP     badDayAdjuster,
    DateTime&                  determinationDate,     // (O)
    DateTime&                  accrualPayDate) const  // (O)
{
    static const string method("CreditFeeLegWithPV::getDeterminationAndPayDate");

    DateTime protectionStart;
    DateTime protectionEnd;

    if (protectionStartDate.empty() || protectionEndDate.empty()) {
        // Assume that protection starts with the first fee's accrue start
        // period and/or ends with the last fee's accrue end.
        AccrualPeriodArrayConstSP accrualPeriods = getAccrualPeriods();

        int numAccrualPeriods = accrualPeriods->size();
        if (numAccrualPeriods == 0) {
            // There is no accrual period information, so there is no
            // meaningful default we can assume
            throw ModelException(method, "There is no contingent leg " 
                                 /*since protection dates are empty */
                                 "and no accrue periods in the fee leg.");
        }

        protectionStart = 
            protectionStartDate.empty() ? accrualPeriods->front()->startDate() :
                                          protectionStartDate;

        protectionEnd = 
            protectionEndDate.empty() ? accrualPeriods->back()->endDate() :
                                        protectionEndDate;
    }
    else {
        protectionStart = protectionStartDate;
        protectionEnd   = protectionEndDate;
    }
    
    double feeLegValue = 0.0;

    if (defaultDate < protectionStart) {
        // If the default happenend before the protection start date, all
        // risky fees will be dropped but riskless fees should still be 
        // included. Set the determination date to the defaultDate.
        determinationDate = defaultDate;
        accrualPayDate = defaultDate;
    }
    else if (defaultDate > protectionEnd) {
        // All fees need to be accrued riskless since the default cannot be 
        // triggered in this contract. No need to look into the override here, 
        // since whatever the eventDeterminationDate etc are, those parameters 
        // will not be used
        // This is flagged by having empty dates
        determinationDate = DateTime();
        accrualPayDate = DateTime();
    }
    else {
        if (!creditEventOverride) {
            // There is no override, so use the values in the instrument
            if (!triggerDelay) {
                // No adjustment is done here.
                determinationDate = defaultDate;
                accrualPayDate = defaultDate;
            }
            else {
                // Note this means both triggerDelay and 
                // defaultToSettlementDelay are not NULL.
                // Obtain the determination and payment dates using this 
                // instruments delay parameters + adjustments
                determinationDate = 
                    ITrancheCreditEventOverride::rollAndAdjustDate(
                        defaultDate,
                        triggerDelay,
                        valuationDate,
                        valuationDate,
                        lastTriggerDate,
                        badDayAdjuster);

                accrualPayDate = 
                    ITrancheCreditEventOverride::rollAndAdjustDate(
                        defaultDate,
                        defaultToSettlementDelay,
                        valuationDate,
                        determinationDate,
                        lastTriggerDate,
                        badDayAdjuster);
            }
        }
        else {
            FeeLegReductionPerDefaultArrayConstSP reductions =
                creditEventOverride->historicFeeLegReductions(
                    1.0, // notional - not used for single names
                    0.0, // recovery rate override - not used for single names
                    lastTriggerDate,
                    defaultDate,
                    badDayAdjuster);
        
            // Sanity-check: verify that all losses have the same determination
            // date and effective date
            int numReductions = reductions->size();
            if (!reductions || (numReductions == 0)) {
                // determinationDate and accrualPayDate are empty 
                determinationDate = DateTime();
                accrualPayDate = DateTime();
            }
            else {
                determinationDate = (*reductions)[0]->determinationDate;
                accrualPayDate = (*reductions)[0]->calculationDate;
                for (int i=1; i < numReductions; ++i) {
                    if (determinationDate != (*reductions)[i]->determinationDate) {
                        throw ModelException(
                            method, "A single name is not expected to "
                            "produce losses with different determination "
                            "dates (" +
                            determinationDate.toString() + " vs " +
                            (*reductions)[i]->determinationDate.toString() +
                            ").");
                    }
                    // The accrualPayDate is the first of all calculationDates
                    if ((*reductions)[i]->calculationDate < accrualPayDate) {
                        accrualPayDate = (*reductions)[i]->calculationDate;
                    }
                }
            }
        }
    }
    return;
}

/** Estimates the known cash flows (so they are not really "known") 
 * corresponding to a CDO tranche - takes into account estimated losses
 * in the future */
CashFlowArraySP CreditFeeLegWithPV::estimateKnownCashFlows( 
    const double         initialTrancheSize,
    const DateTimeArray  pastTrancheLossDates,
    const DoubleArray    pastTrancheLosses,
    IForwardRatePricerSP model)
{
    static const string method("CreditFeeLegWithPV::estimateKnownCashFlows");

    try {
        // get the basic fee leg cashflows (these assume full accrual)
        AbstractCashFlowArrayConstSP cfl = getCashFlows(model);
        int numCashFlows = cfl->size();
        CashFlowArraySP fkcfl(new CashFlowArray(numCashFlows));

        for (int i=0; i < numCashFlows; ++i) {
            (*fkcfl)[i].date   = (*cfl)[i]->getPayDate();
            (*fkcfl)[i].amount = (*cfl)[i]->getAmount(model);
            
            //manipulate any credit type cashflows
            if (CreditCashFlow::TYPE->isInstance((*cfl)[i])) {
                CreditCashFlow* ccf = 
                    static_cast<CreditCashFlow*>(((*cfl)[i]).get());
                
                CashFlow cf = ccf->getKnownCashFlow(initialTrancheSize,
                                                    pastTrancheLossDates,
                                                    pastTrancheLosses,
                                                    model);
                
                (*fkcfl)[i].amount = cf.amount;
            }
        }
        return fkcfl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns the accrued interest of a non-defaulted CDS*/
double CreditFeeLegWithPV::getFeeLegAI(
    const DateTime&              valuationDate, 
    const DateTime&              paymentDate,
    const DateTime&              earliestAccrualStart,
    const DateTime&              latestAccrualEnd,
    const DayCountConventionSP   dcc, //allows an override to be specified
    const IDiscountCurveRisky&   crv,
    const IDiscountCurveConstSP  discount,
    const IDecretionCurveConstSP prepay,
    IForwardRatePricerSP         model) const
{
    static const string method = "FullCreditFeeLeg::getFeeLegAI";

    CashFlowArraySP riskyFees = getRiskyCashFlows(model);
    CashFlowArraySP risklessFees = getRisklessCashFlows(model);

    //protection dates can be estimated
    DateTime protectionStartDate = earliestAccrualStart;
    DateTime protectionEndDate = latestAccrualEnd;

    DefaultRatesSP defRates = crv.getDefaultRates();

    CDSPricer cdsPricer = CDSPricer(riskyFees,
                                    risklessFees,
                                    defRates,
                                    1.0,  //notional,
                                    1.0,  //recovery
                                    true, //payAccruedFee,
                                    dcc,  
                                    paymentDate,
                                    protectionStartDate,
                                    protectionEndDate,
                                    protectionStartDate, //accruedEffectiveDate,
                                    discount,
                                    prepay);

    double accruedInterest = cdsPricer.calculateDefaultPayment(true);

    return accruedInterest;
}

/* Gets the value of the fee leg corresponging to a defaulted CDS 
   under the (optional) credit event */
double CreditFeeLegWithPV::getFeeLegDefaultedPV(
    const DateTime&              valuationDate,
    const DateTime&              defaultDate,
    const DateTime&              protectionStartDate,
    const DateTime&              protectionEndDate,
    const bool                   allowIncludingTodaysPayments,
    const IDiscountCurveConstSP  discount,
    const IDecretionCurveConstSP prepay,
    IForwardRatePricerSP         model,
    IBadDayAdjusterConstSP       badDayAdjuster,
    ICreditEventOverrideNameSP   creditEventOverride,
    CIntSP                       triggerDelay,
    CIntSP                       defaultToSettlementDelay,
    DateTime                     lastTriggerDate) const 
{
    static const string method("CreditFeeLegWithPV::getFeeLegDefaultedPV");

    CashFlowArraySP knownCashFlows = 
        generateKnownCashFlowsGivenDefault(
            valuationDate,
            defaultDate,
            valuationDate, // exclude payments paid before today
            protectionStartDate,
            protectionEndDate,
            allowIncludingTodaysPayments,
            model,
            creditEventOverride,
            badDayAdjuster,
            triggerDelay,
            defaultToSettlementDelay,
            lastTriggerDate);

    // Prepay until default date - sure?
    const double factor = 
        (defaultDate > valuationDate) ? prepay->getFactor(defaultDate) : 1.0;
    const double balance = 
        (defaultDate > valuationDate) ? prepay->pv(defaultDate) : 1.0;

    // PV and sum all cashflows
    double feeLegValue(0.0);
    const int numKnownCashFlows = knownCashFlows->size();
    for (int i=0; i < numKnownCashFlows; ++i) {
        const double pv = discount->pv(valuationDate, (*knownCashFlows)[i].date);
        feeLegValue += (*knownCashFlows)[i].amount * pv * factor * balance;
    }

    return feeLegValue;
}

/* Computes the interest accrued up to valuationDate */
double CreditFeeLegWithPV::getAccruedInterest(const DateTime& valuationDate,
                                              IForwardRatePricerSP model) const
{
    static const string method("CreditFeeLegWithPV::getAccruedInterest");

    try {
        // Consider only risky cashflows
        AbstractCashFlowArrayConstSP cashFlows = 
            getAbstractCashFlows(true, model);

        double ai = 0.0;
        int numCashFlows = cashFlows->size();
        for (int idx=0; idx < numCashFlows; ++idx) {
            AbstractCashFlowSP cf = (*cashFlows)[idx]; // For ease

            const CashFlow& newCashFlow = 
                cf->getAccruedCashFlow(valuationDate, model);
            
            // Compute the accrue including cashflows on valuationDate
            if (!Maths::isZero(newCashFlow.amount) && 
                (newCashFlow.date >= valuationDate)) 
            {
                ai += newCashFlow.amount;
            }
        }
        return ai;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}


/** 'Constructor' for CreditFeeLegWithPV*/
CreditFeeLegWithPV* CreditFeeLegWithPV::buildCreditFeeLeg(
    double                     coupon,
    const DateTimeArray&       observationDates, // observation dates
    const DateTimeArray&       accrualStartDates, // accrual dates : can be non contiguous
    const DateTimeArray&       accrualEndDates,
    const DateTimeArray&       payDates,     
    const CouponNotionalTypes& couponNotionalType,
    const string&              payDCC,
    const IModel*              model, 
    const MarketData*          market)
{
	return new CreditFeeLeg(
					coupon,
					observationDates, // observation dates
					accrualStartDates, // accrual dates : can be non contiguous
					accrualEndDates,
					payDates,     
					couponNotionalType,
					payDCC,
					model, 
					market
					);
}

double CreditFeeLegWithPV::getFeeLegPV(const DateTime&              valuationDate, 
                                       const DateTime&              earliestRiskyDate,
                                       const DateTime&              latestRiskyDate,
                                       const IDiscountCurve&        discount,
                                       const IDiscountCurveRisky&   crv,
                                       const IDecretionCurveConstSP prepay,
                                       const bool                   includeAccrued,
                                       const DayCountConventionSP   dcc,
                                       IForwardRatePricerSP         model) const
{
    return getFeeLegPV(valuationDate, valuationDate,
                       earliestRiskyDate, latestRiskyDate,
                       discount, crv, prepay,
                       includeAccrued, dcc, model);
}

double CreditFeeLegWithPV::getFeeLegPV(const DateTime&              valuationDate, 
                                       const DateTime&              paymentDate, 
                                       const DateTime&              earliestRiskyDate,
                                       const DateTime&              latestRiskyDate,
                                       const IDiscountCurve&        discount,
                                       const IDiscountCurveRisky&   crv,
                                       const IDecretionCurveConstSP prepay,
                                       const bool                   includeAccrued,
                                       const DayCountConventionSP   dcc,
                                       IForwardRatePricerSP         model) const
{
    return getFeeLegPV(valuationDate, paymentDate,
                       earliestRiskyDate, latestRiskyDate,
                       discount, crv, prepay,
                       includeAccrued, dcc, false, model);
}

double CreditFeeLegWithPV::getFeeLegPV(const DateTime&              valuationDate, 
                                       const DateTime&              paymentDate,
                                       const DateTime&              earliestRiskyDate,
                                       const DateTime&              latestRiskyDate,
                                       const IDiscountCurve&        discount,
                                       const IDiscountCurveRisky&   crv,
                                       const IDecretionCurveConstSP prepay,
                                       const bool                   includeAccrued,
                                       const DayCountConventionSP   dcc,
                                       bool                         defaultValueOnly,
                                       IForwardRatePricerSP         model) const
{
   static const string method = "CreditFeeLegWithPV::getFeeLegPV";

   int i;
   double fees = 0.0;

    // Risky fee payments
    CashFlowArrayConstSP riskyFees = getRiskyCashFlows(model);

    //valued at payment date, dependent upon survival from valuation date to payment date
    double sp = crv.survivalProb(valuationDate, paymentDate);

    for (i=0; i<riskyFees->size(); i++)
    {
        const DateTime& cfDate = (*riskyFees)[i].date;
        //for non-ABCDS, settles(cfDate) gives back cfDate since by default there is no holiday
        //no weekends adjustment in prepay::prepaySettle
        DateTime prepaySettleDate = prepay->getSettlement()->settles(cfDate);
        if (prepaySettleDate.isGreater(valuationDate))
        {
            DateTime toRiskyDate = (cfDate > latestRiskyDate) ?
                                   latestRiskyDate :
                                   cfDate;

            // default after protectionEndDate doesn't effect fee payment
            double survPV = 1;
            double prepaypv = 1;
            if(cfDate.isGreater(valuationDate))
            {
                survPV = crv.survivalProb(paymentDate, toRiskyDate);
                prepaypv = prepay->pv(paymentDate, cfDate);
            }

            double discPV = crv.risklessPV(paymentDate, prepaySettleDate);
            fees += sp * (*riskyFees)[i].amount * survPV * discPV * prepaypv;
        } 
    }
    
    double factor = prepay->getFactor(valuationDate);
    fees *= factor;

    // Riskless fee payments
    CashFlowArrayConstSP risklessFees = getRisklessCashFlows(model);

    for (i=0; i<risklessFees->size(); i++)
    {
        const DateTime& cfDate = (*risklessFees)[i].date;
        if (cfDate.isGreater(valuationDate))
        {
            double discPV = crv.risklessPV(valuationDate, cfDate);
            //double prepaypv = prepay->pv(valuationDate, cfDate);

            fees += (*risklessFees)[i].amount * discPV; // * prepaypv;
        } 
    }

    // optionally include accrued interest
    double ai = 0.0;
    if (includeAccrued)
    {
        ai = sp * getFeeLegAI(valuationDate, paymentDate, 
                              earliestRiskyDate, latestRiskyDate,
                              dcc, //allows an override to be specified
                              crv, IDiscountCurveConstSP::attachToRef(&discount), prepay,
                              model);
    }

    return fees + ai;
}

DRLIB_END_NAMESPACE
