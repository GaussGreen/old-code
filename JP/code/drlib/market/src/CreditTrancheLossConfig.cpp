//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : An ICreditLossConfig used for CDO tranche products.
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/IProtectionProvider.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/IEffectiveCurveLossGen.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/DefaultRates.hpp"
#include "edginc/SingleCreditAsset.hpp"

DRLIB_BEGIN_NAMESPACE


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                           Two static methods start
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/** Computes a rebate payment given the original tranche and portfolio
    reductions plus the assumed and real cashflows.
    The reason why it returns a CashFlowArraySP rather than a CashFlow(SP)
    is for better integration with the calling methods, which will typically
    merge the results of this method call with more rebate payments.
    The "trancheReductions" are modified internally in order to incorporate
    the real reductions seen so far. */
CashFlowArraySP computeOneRebate(
    const CashFlowArrayConstSP assumedReductions,
    const CashFlowArrayConstSP realReductions,
    CashFlowArraySP& trancheReductions,        // (I) AND (O)
    const DateTime& rebateCalcDate,
    const IRebateCalculator* const rebateCalc,
    const double initialNotional,
    IForwardRatePricerSP model)
{
    // Clone and change the sign of the assumed reductions
    CashFlowArraySP negativeAssumedReductions(assumedReductions.clone());
    for (int i=0; i < negativeAssumedReductions->size(); ++i) {
        (*negativeAssumedReductions)[i].amount *= -1;
    }

    // Remove these assumed reductions from the list of real reductions
    CashFlowArraySP realTrancheReductions =
        CashFlow::merge(trancheReductions, negativeAssumedReductions);

    // ... and add the real reductions
    realTrancheReductions = CashFlow::merge(realTrancheReductions,
                                            realReductions);

    // Obtain the CUMULATIVE REAL tranche reductions
    CashFlowArraySP cumulativeRealTrancheReductions =
        CashFlow::accumulateCashFlows(realTrancheReductions);

    CashFlowArraySP cumulativeAssumedTrancheReductions =
        CashFlow::accumulateCashFlows(trancheReductions);

    double rebate =
        rebateCalc->computeRebate(cumulativeAssumedTrancheReductions,
                                  cumulativeRealTrancheReductions,
                                  initialNotional,
                                  model);

    CashFlowArraySP rebatePayments;
    if (Maths::isZero(rebate)) {
        rebatePayments.reset(new CashFlowArray(0));
    }
    else {  // Produce the rebate cashflow.
        rebatePayments.reset(new CashFlowArray(
            1, CashFlow(rebateCalcDate, rebate)));
    }

    // Set the trancheReductions to the realTrancheReductions
    trancheReductions.reset(realTrancheReductions.get());

    return rebatePayments;
}


/** Given the fee leg TRANCHE reductions, computes the corresponding rebates,
    ie the fraction of the fees that has been paid when it should not have if
    the final recovery rate of the name had been known at the fee payment date).
    The algorithm is as follows:
    1. Start with the original distribution of tranche reductions,
    2. Find the first rebate date for any of those reductions.
    3. At that point all reductions with that rebate date should be applied
       to their corresponding determination date (ie, as if the information
       had been know from the trigger date).
       The difference in fees (ie, the rebate) is then computed.
    4. Set the new reduction profile as the original reduction profile: since
       that information is now known and the corresponding rebate(s) have
       been paid, it is as if that information had been known at the time.
    5. Go back to 2, looping over all rebate dates */
CashFlowArraySP computeRebatesForTrancheReductions(
    FeeLegReductionPerDefaultArraySP trancheReductions,
    const IRebateCalculator* const rebateCalc,
    const double initialNotional,
    IForwardRatePricerSP model,
    const bool doLosses)
{
    // (1) Set the original tranche reductions cashflows - note they are NOT
    // cumulative at this point
    CashFlowArraySP trancheReductionsCF =
        FeeLegReductionPerDefault::getReductions(trancheReductions, doLosses);

    CashFlowArraySP rebates(new CashFlowArray(0));
    DateTimeArray visitedRebateDates(1, DateTime());
    int numFeeReductions = trancheReductions->size();

    while (true) { // only exit via a break, below
        // (2) Find the smallest rebate date greater than all dates in
        // visitedRebateDates
        bool moreReductionsToVisit = false;
        DateTime minRebateDate;
        for (int i=0; i < numFeeReductions; ++i) {
            const DateTime& reductionRebateDate = (*trancheReductions)[i]->calculationDate;

            if (reductionRebateDate > visitedRebateDates.back()) {
                if ((minRebateDate.empty()) || // first time here, or
                    (reductionRebateDate < minRebateDate)) // smaller than previous min
                {
                    minRebateDate = reductionRebateDate;
                    moreReductionsToVisit = true;
                }
            }
        }
        if (!moreReductionsToVisit) {
            // we are done - get out of the while(true) loop
            break;
        }
        // add this rebateDate to the array of visited dates
        visitedRebateDates.push_back(minRebateDate);

        // (3.1) Get all reductions with that rebate date
        CashFlowArraySP realReduction(new CashFlowArray(0));
        CashFlowArraySP assumedReduction(new CashFlowArray(0));
        for (int i=0; i < numFeeReductions; ++i) {
            const DateTime& reductionRebateDate =
                (*trancheReductions)[i]->calculationDate;
            if (reductionRebateDate == minRebateDate) {
                const DateTime& determinationDate =
                    (*trancheReductions)[i]->determinationDate;

                double amount =
                    (*trancheReductions)[i]->getReductionAmount(doLosses);

                realReduction = CashFlow::merge(
                    realReduction,
                    CashFlowArraySP(new CashFlowArray(
                        1, CashFlow(determinationDate, amount))));

                assumedReduction = CashFlow::merge(
                    assumedReduction,
                    CashFlowArraySP (new CashFlowArray(
                        1, CashFlow((*trancheReductions)[i]->effectiveDate,
                                    amount))));
            }
        }

        // (3.2) and (4): ComputeOneRebate will remove the assumedReductions
        // and add the realReductions to trancheReductionsCF (which will be
        // modified after this call)
        CashFlowArraySP rebateForNamesReduction = computeOneRebate(
                assumedReduction,    // assumed reductions
                realReduction,       // real reductions
                trancheReductionsCF, // note it will be modified here!
                minRebateDate, // rebate calc date
                rebateCalc,
                initialNotional,
                model);

        // Add the delivery rebate to the general rebate payments array
        rebates = CashFlow::merge(rebates, rebateForNamesReduction);
    }
    return rebates;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                          Two static methods end
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Implementation of a parallel shift in strikes
TweakOutcome CreditTrancheLossConfig::sensShift(
    const PropertyTweak<CDOParallelStrike>& shift)
{
    lowStrike  += shift.coefficient;
    highStrike += shift.coefficient;
    return TweakOutcome(shift.coefficient, false); //
}

string CreditTrancheLossConfig::sensName(const CDOParallelStrike* shift) const {
    return "CreditTrancheLossConfig"; // please don't return "" since that means "don't tweak me"
}

void CreditTrancheLossConfig::sensRestore(
    const PropertyTweak<CDOParallelStrike>& shift)
{
    lowStrike  -= shift.coefficient;
    highStrike -= shift.coefficient;
}


// Constructor for tranches without portfolio - i.e. used as pure
// "low strke" - "high strike" payoff
CreditTrancheLossConfig::CreditTrancheLossConfig(
    double lowStrike,
    double highStrike) :
        CObject(TYPE),
        lowStrike(lowStrike),
        highStrike(highStrike),
		isXSub(false),
		sumBufferForXSub(0),
		isOuter(false)
{
	hasCutOff = !(cutOffDate.empty());
}

CreditTrancheLossConfig::CreditTrancheLossConfig() :
    CObject(TYPE),
	isXSub(false),
	sumBufferForXSub(0),
	isOuter(false)
{
	hasCutOff = !(cutOffDate.empty());
}

// Public constructor, with all fields passed in
CreditTrancheLossConfig::CreditTrancheLossConfig(
        const string& name,
        const DateTime& valueDate,
        const double lowStrike,
        const double highStrike,
        ICreditLossConfigSP lossCfg) :
    CObject(TYPE),
    name(name),
    valueDate(valueDate),
    lowStrike(lowStrike),
    highStrike(highStrike),
	portfolio(lossCfg),
	isXSub(false),
	sumBufferForXSub(0),
	isOuter(false)
{
    validatePop2Object();
}


CreditTrancheLossConfig::~CreditTrancheLossConfig()
{}


void CreditTrancheLossConfig::validatePop2Object() {
    static const string method("CreditTrancheLossConfig::validatePop2Object");

    /* check strikes */
    if (lowStrike >= highStrike){
        throw ModelException(method, "LowStrike (" +
                             Format::toString(lowStrike) + ") >= high "
                             "strike (" + Format::toString(highStrike) + ")");
    }

    if (!!portfolio ) {
        // Skip loss configs where no portfolio has been provided (using the
        // constructor that only has the strikes)

        // Convert the poftfolio from ICreditLossConfig into CDOPortfolio, and set it
        // as the underlying portfolio in this CreditTrancheLossConfig
        CDOPortfolio* cdoPortf = dynamic_cast<CDOPortfolio*>(portfolio.get());
        if (!cdoPortf) {
            throw ModelException (method, "Cannot construct a "
                                  "CreditTrancheLossConfig if the supplied "
                                  "loss config is not of type CDOPortfolio.");
        }
        // JLHP - although a hack, still need access to the underlying CDOPortfolio
        // as a CDOPortfolio, so keep a smartpointer to it.
        // This ought to be removed
        portf = CDOPortfolioConstSP::attachToRef(cdoPortf);
    }

	hasCutOff = !(cutOffDate.empty());

	if (isXSub) {
        sumBufferForXSub = 0;
        isOuter = true;
        for (int i=0; i < numInnerLossConfigs(); ++i) {
            const PortfolioName* pname =
                dynamic_cast<const PortfolioName*>(getInnerLossConfig(i).get());

            // components are either single names or CDO's
            if (pname) {
                // do we need to change buffer?
            }
            else {
                const CreditTrancheLossConfig* ctl =
                    dynamic_cast<const CreditTrancheLossConfig*>
					(getInnerLossConfig(i).get());

                if (ctl) {
                    ctl->isXSub = true;
                    ctl->isOuter = false;

                    double k1, k2;

                    // get the original strikes
                    ctl->getTrancheStrikes(k1, k2);

                    sumBufferForXSub += k1;

                    const CDOPortfolio* lc2 =
                        dynamic_cast<const CDOPortfolio*>(ctl->portfolio.get());

                    if (!lc2) {
                        throw ModelException("Odd: CreditTrancheLossConfig "
                                             "should contain a CDOPortfolio");
                    }
                    for (int j = 0; j < lc2->numInnerLossConfigs(); j++) {
                        pname = dynamic_cast<const PortfolioName*>
                            ( lc2->getInnerLossConfig(j).get() );

                        if (!pname) {
                            throw ModelException(method, "Not a CDO^2");
                        }
                    }
                }
                else {
                    throw ModelException(method, "Only a CDO^2 can be XSub");
                }
            }
        }
    }
}


// jlhp comments
IObject* CreditTrancheLossConfig::clone() const {
    CreditTrancheLossConfig* copy = DYNAMIC_CAST(CreditTrancheLossConfig,
                                                 CObject::clone());

    // Convert the lossCfg from ICreditLossConfig into CDOPortfolio, and set it
    // as the underlying portfolio in the copy CreditTrancheLossConfig
    copy->portf.reset(dynamic_cast<CDOPortfolio*>(copy->portfolio.get()));
    return copy;
}

//++++++++ ICreditLossConfig methods
/** Name of the underlying asset */
string CreditTrancheLossConfig::getName() const {
    return name;
}


/** Returns an array with the cashflows corresponding to all "real"
    losses due to the settlement of this name's default.
    If there are no losses, it may return an empty (0 size)
    CtgLegLossPerDefaultArraySP or a "null" SP.
    Note: the cashflow dates are the CALCULATION DATES (as opposed to
    settlement dates) when the corresponding losses are determined. */
CtgLegLossPerDefaultArraySP CreditTrancheLossConfig::historicContingentLegLosses(
    CIntConstSP triggerDelay,
    CIntConstSP defaultToCalculationDelay,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda,
    const IProtectionProvider* const protect) const
{
    static const string method("CreditTrancheLossConfig::historicContingentLegLosses");

    CtgLegLossPerDefaultArraySP portfolioLosses =
        portfolio->historicContingentLegLosses(
            triggerDelay,
            defaultToCalculationDelay,
            lastTriggerDate,
            bda,
            NULL); // Pass null IProtectionProvider - verification is done here

    if (protect) {
        // Validate those losses - will throw an exception if there are issues
        CDOPortfolio::validateCtgLegLosses(portfolioLosses, protect);
    }
    else {
        ; // Nothing to do here - this will be the case, eg, in inner CDOs in a CDO2
    }

    CtgLegLossPerDefaultArraySP cLegPastTrancheLosses =
        getTrancheReductions(portfolioLosses, lowStrike, highStrike);

    return cLegPastTrancheLosses;
}


/** Compute the rebate payments caused by historic fee leg reductions */
CashFlowArraySP CreditTrancheLossConfig::historicRebatePayments(
    const IRebateCalculator* const   rebateCalc,
    FeeLegReductionPerDefaultArraySP trancheReductions,
    IForwardRatePricerSP             model,
    const bool                       recoverNotional) const
{
    static const string method("CreditTrancheLossConfig::historicRebatePayments");

    if (!rebateCalc) {
        // If no rebate calculator is provided, no rebates can be returned.
        // Rather than returning zero rebates, take this is an error
        throw ModelException(method, "Internal error: rebateCalc is null");
    }

    CashFlowArraySP rebates;
    const double initialNotional = notional();
    if (!!trancheReductions) {
        // Compute the rebates for losses
        rebates = computeRebatesForTrancheReductions(trancheReductions,
                                                     rebateCalc,
                                                     initialNotional,
                                                     model,
                                                     true); // do losses

        // Compute the rebates for recovered notionals
        if (recoverNotional) {
            CashFlowArraySP rebateForRecovered =
                computeRebatesForTrancheReductions(trancheReductions,
                                                   rebateCalc,
                                                   initialNotional,
                                                   model,
                                                   false); // do rec notional

            // Add the recovered rebates to the general rebates array
            rebates = CashFlow::merge(rebates, rebateForRecovered);
        }
    }
    return rebates;
}


// In a tranche, if a default settles on an accrue period other than the one
// the default was triggered on, we need to make some assumptions and
// adjustments to the associated fee leg losses: All losses happening on
// later accrue periods need to be aligned to the corresponding period's
// accrue start date (this is the contractual definition). Before that point,
// the "temporaryLossAmount" defined in the contract should be used.
FeeLegReductionPerDefaultArraySP alignReductionsToAccruePeriods(
    FeeLegReductionPerDefaultArrayConstSP portfolioReductions,
    AccrualPeriodArrayConstSP accrualPeriods,
    const double temporaryLossAmount)
{
    static const string method("CreditTrancheLossConfig-"
                               "alignReductionsToAccruePeriods");

    try {
        int numAccrualPeriods = accrualPeriods->size();

        if (numAccrualPeriods == 0) {
            throw ModelException(method, "Internal error: there is a fee leg "
                                 "but there are no accrual periods in it.");
        }

        FeeLegReductionPerDefaultArraySP alignedReductions(
            new FeeLegReductionPerDefaultArray(0));

        const int numPtfReductions = portfolioReductions->size();
        for (int idx=0; idx < numPtfReductions; ++idx) {
            const FeeLegReductionPerDefaultSP reduction = (*portfolioReductions)[idx];

            if (!reduction) {
                // This particular reduction is null (ie, no reduction). Skip
                // and continue with the next one
                continue;
            }

            const DateTime& determinationDate      = reduction->determinationDate;
            const DateTime& calculationDate        = reduction->calculationDate;
            const double defaultedNotionalFraction = reduction->defaultedNotional;
            const double totalNameNotional         = reduction->totalNameNotional;
            const double recRate                   = reduction->recoveryRate;

            // Sanity check - currently we should only "align" the reductions in
            // one place (here) and just once (now) - so make sure the input
            // dates are unaltered
            if (calculationDate != reduction->effectiveDate) {
                throw ModelException(method, "Internal error: one of the inner "
                                     "names produced a fee leg reduction with "
                                     "calculation date (" +
                                     calculationDate.toString() +
                                     ") different to the effective date ("+
                                     reduction->effectiveDate.toString() +").");
            }

            // Obtain which accrual period the determinationDate date falls in
            int determinationDatePeriod = -1;
            bool determinationDateAfterLastPeriod = true;

            for (int i=0; i < numAccrualPeriods; ++i) {
                if (determinationDate <= (*accrualPeriods)[i]->endDate()) {
                    // determinationDate does not fall after the last period
                    determinationDateAfterLastPeriod = false;

                    if ((*accrualPeriods)[i]->startDate() <= determinationDate) {
                        // determinationDate falls in this period
                        // CAUTION: Ignores the fact that any given date may fall on
                        // two or more accrual periods if they overlap
                        determinationDatePeriod = i;
                        break;
                    }
                }
            }

            FeeLegReductionPerDefaultArraySP singleAlignedReductions;
            if (determinationDateAfterLastPeriod) {
                // The determination date falls after the last accrual period.
                // No need to worry about the losses/recovered notional because
                // all fees have been paid (or at least their values have been
                // determined) at this point.
                // Assume no losses/recovered notional (leave as null)
            }
            else {
                // Being here, we have to compute the cashflow
                // corresponding to the real loss as known on calculation date - this
                // cashflow will take effect on "calculationEffectiveDate", computed
                // below
                DateTime calculationEffectiveDate;

                if (determinationDate == calculationDate) {
                    calculationEffectiveDate = calculationDate; // nothing to align here
                }
                else if (determinationDatePeriod == -1) {
                    // The determination date falls between accrual periods.
                    // The actual loss applies from the beginning of the first
                    // period ending after calculationDate, so determine what
                    // this is.
                    DateTime minPeriodEndDate; // empty date
                    int calculationDatePeriod = -1;
                    bool calculationDateAfterLastPeriod = true;
                    for (int i=0; i < accrualPeriods->size(); ++i) {
                        if (calculationDate <= (*accrualPeriods)[i]->endDate()) {
                            // determinationDate does not fall after the last period
                            calculationDateAfterLastPeriod = false;

                            // This accrual period ends after calculation date. We are
                            // interested in the period with the smallest ending date.
                            if (minPeriodEndDate.empty() || // 1st time round
                                ((*accrualPeriods)[i]->endDate()) < minPeriodEndDate)
                            {
                                minPeriodEndDate = (*accrualPeriods)[i]->endDate();
                                calculationDatePeriod = i;
                            }
                        }
                    }

                    if (calculationDateAfterLastPeriod) {
                        // calculationDate falls after the last accrual period, so the
                        // real losses / recovered notional will not impact any future
                        // fees.
                        // Arbitrarily, set calculationEffectiveDate to be 2 days after
                        // the last accrual period end
                        calculationEffectiveDate =
                            accrualPeriods->back()->endDate().rollDate(2);
                    }
                    else {
                        if (calculationDatePeriod < 0) {
                            throw ModelException(method,
                                                 "Internal error: calculationDatePeriod"
                                                 " falls before the last accrual period,"
                                                 " but could not identify the first "
                                                 "period immediately after it");
                        }
                        calculationEffectiveDate =
                            (*accrualPeriods)[calculationDatePeriod]->startDate();
                    }
                }
                else {
                    // The determination date falls in an accrual period
                    // Obtain the accrual period the CALCULATION date falls on
                    int calculationDatePeriod = -1;
                    bool calculationDateAfterLastPeriod = true;
                    for (int i=0; i < accrualPeriods->size(); ++i) {
                        if (calculationDate <= (*accrualPeriods)[i]->endDate()) {
                            // determinationDate does not fall after the last period
                            calculationDateAfterLastPeriod = false;

                            if (calculationDate >= (*accrualPeriods)[i]->startDate()) {
                                // determinationDate falls in this period
                                calculationDatePeriod = i;
                                break;
                            }
                        }
                    }

                    if (calculationDatePeriod == determinationDatePeriod) { // != -1
                        // If calculation date happens in the same accrual period as
                        // the determination date, the real losses/recovered amounts
                        // take effect from the determination date
                        calculationEffectiveDate = determinationDate;
                    }
                    else if (calculationDateAfterLastPeriod) {
                        // calculationDate falls after the last accrual period, so the
                        // real losses/recovered amounts will not impact any future
                        // fees.
                        // Arbitrarily, set calculationEffectiveDate to be 2 days after
                        // the last accrual period end
                        calculationEffectiveDate =
                            accrualPeriods->back()->endDate().rollDate(2);
                    }
                    else if (calculationDatePeriod != -1) {
                        // The calculation date happens on the accrual period number
                        // "calculationDatePeriod", different from the accrual period
                        // when event determination happens. The real losses/recovered
                        // amounts take effect from the beginning of that accrual
                        // period
                        calculationEffectiveDate =
                            (*accrualPeriods)[calculationDatePeriod]->startDate();
                    }
                    else {
                        // The calculation date falls between accrual periods.
                        // The actual loss takes effect from the beginning of the first
                        // accrual period ending after calculationDate, so find out
                        // what that date is.
                        // Note at this point we know that the calculation date falls
                        // before the end of the last period
                        DateTime minPeriodEndDate; // empty date
                        int calculationDatePeriod = -1;
                        for (int i=0; i < accrualPeriods->size(); ++i) {
                            if (calculationDate <= (*accrualPeriods)[i]->endDate()) {
                                if (minPeriodEndDate.empty() || // 1st time round
                                    (minPeriodEndDate > (*accrualPeriods)[i]->endDate()))
                                {
                                    minPeriodEndDate = (*accrualPeriods)[i]->endDate();
                                    calculationDatePeriod = i;
                                }
                            }
                        }

                        if (calculationDatePeriod < 0) {
                            throw ModelException(method,
                                                 "Internal error: calculationDatePeriod"
                                                 " falls before the last accrual period,"
                                                 " but could not identify the first "
                                                 "period immediately after it");
                        }
                        calculationEffectiveDate =
                            (*accrualPeriods)[calculationDatePeriod]->startDate();
                    }
                }


                if (determinationDate == calculationEffectiveDate) {
                    // No need to use the temporaryLossAmount, since we can determine
                    // the real losses and recovered amounts in the same accrual period
                    // as the determination date, ie, here portfolio loss is the REAL
                    // loss etc
                    singleAlignedReductions.reset(new FeeLegReductionPerDefaultArray(
                        1, FeeLegReductionPerDefaultSP(new FeeLegReductionPerDefault(
                            determinationDate,
                            determinationDate, // =calculationEffectiveDate
                            calculationDate, // rebate date (no rebate here)
                            defaultedNotionalFraction,
                            totalNameNotional,
                            recRate))));
                }
                else {
                    // Create the temporary loss assumption first.
                    // Use the totalNameNotional here, assuming at this point that all
                    // the name's notional will be triggered for protection
                    singleAlignedReductions.reset(new FeeLegReductionPerDefaultArray(2));
                    (*singleAlignedReductions)[0] = FeeLegReductionPerDefaultSP(
                         new FeeLegReductionPerDefault(
                             determinationDate,
                             determinationDate, // effective date
                             calculationDate,
                             totalNameNotional, // defaulted notional = all
                             totalNameNotional,
                             1.0 - temporaryLossAmount)); // recovery rate

                    // Now the "real" reduction. This is a bit tricky: need one
                    // FeeLegReductionPerDefault that will produce the following:
                    // - loss = defaultedNotionalFraction * (1.0-recRate) -
                    //          (totalNameNotional * temporaryLossAmount)
                    // - recovered = -loss
                    // This can be done with one FeeLegReductionPerDefault where:
                    // - totalNameNotional = 0 (so that recovery = -loss)
                    // - Loss = N(1-R) =
                    //        = 1 * [defaultedNotionalFraction * (1-recRate) -
                    //          (totalNameNotional * temporaryLossAmount)] =
                    //        = 1 * {1 - [1 + (defaultedNotionalFraction * (recRate-1) +
                    //          totalNameNotional * temporaryLossAmount]}
                    // Note: The recovery rate is no longer a recovery rate here...
                    const double effRecRate =
                        1.0 + defaultedNotionalFraction * (recRate-1) +
                        totalNameNotional * temporaryLossAmount;

                    (*singleAlignedReductions)[1] = FeeLegReductionPerDefaultSP(
                         new FeeLegReductionPerDefault(
                             determinationDate,
                             calculationEffectiveDate, // reduction date
                             calculationDate,
                             1.0, // defaultedNotional (hacky... see above)
                             0.0, // totalNameNotional (hacky... see above)
                             effRecRate));
                }
            }

            // Now merge the singleAlignedReductions into alignedReductions
            if (!!singleAlignedReductions) {
                int numSingleReductions = singleAlignedReductions->size();
                alignedReductions->reserve(alignedReductions->size() +
                                           numSingleReductions);
                for (int i=0; i < numSingleReductions; ++i) {
                    alignedReductions->push_back((*singleAlignedReductions)[i]);
                }
            }
        }
        return alignedReductions;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Computes the notional reductions on the fee leg (due to losses and/or
    recovered notional) and the corresponding fee rebates, if any - based
    on the creditEventOverride in this name.
    If there are no losses, it may return an empty (0 size)
    CtgLegLossPerDefaultArraySP or a "null" SP.
    The accrualPeriods are an argument because they are required to
    determine the actual date of the notional reductions (in the fee leg,
    all notional reductions happen either on the determination date or on
    an accrual start date. */
FeeLegReductionPerDefaultArraySP
CreditTrancheLossConfig::historicFeeLegReductions(
    CIntConstSP triggerDelay,
    CIntConstSP defaultToCalculationDelay,
    const double temporaryLossAmount,
    const DateTime& lastTriggerDate,
    AccrualPeriodArrayConstSP accrualPeriods,
    IBadDayAdjusterConstSP bda,
    const bool recoverNotional) const
{
    // Note the accrual periods do NOT get propagated down: they are used
    // at this (top) level to adjust the reductions coming from the inner
    // names/tranches. This adjustment only needs to take place once, here.
    FeeLegReductionPerDefaultArraySP portfolioReductions =
        portfolio->historicFeeLegReductions(triggerDelay,
                                            defaultToCalculationDelay,
                                            temporaryLossAmount,
                                            lastTriggerDate,
                                            AccrualPeriodArrayConstSP(),
                                            bda,
                                            recoverNotional);

    if (!!accrualPeriods) {
        // Need to adjust the fee leg reductions to the tranche accrue periods
        portfolioReductions = alignReductionsToAccruePeriods(
            portfolioReductions, accrualPeriods, temporaryLossAmount);
    }

    // Collar the portfolio reductions
    FeeLegReductionPerDefaultArraySP trancheReductions =
        getTrancheReductions(portfolioReductions,
                             lowStrike,
                             highStrike,
                             true); // true because these are losses

    // If recovering notional, collar those reductions too
    if (recoverNotional) {
        const double totalLongPortfolioNotional = portfolioLongNotional();
        FeeLegReductionPerDefaultArraySP trancheRecovered =
            getTrancheReductions(portfolioReductions,
                                 totalLongPortfolioNotional - highStrike,
                                 totalLongPortfolioNotional - lowStrike,
                                 false); // false because recovering notional
        // Merge all reductions
        int numRecovered = trancheRecovered->size();
        trancheReductions->reserve(trancheReductions->size() + numRecovered);
        for (int i=0; i < numRecovered; ++i) {
            trancheReductions->push_back((*trancheRecovered)[i]);
        }
    }
    return trancheReductions;
}


/** populate from market cache */
void CreditTrancheLossConfig::getMarket(const IModel* model,
                                        const MarketData* market)
{
    market->GetReferenceDate(valueDate);
    portfolio->getMarket(model, market);
}


/** Returns the engine parameters for this loss cfg (NB Model needs to
    have selected appropriate set during getMarket()). The parameter
    can specify what type of engine parameters are required. An
    exception is thrown if the relevant type is not found */
CreditEngineParametersConstSP CreditTrancheLossConfig::getEngineParams(
    CClassConstSP engineParamsType) const
{
    static const string method("CreditTrancheLossConfig::getEngineParams");
    try {
        if (!engineParamsType){
            throw ModelException(method, "The type of engine parameters "
                                 "requested is Null");
        }
        if (!engineParams){
			return portfolio->getEngineParams(engineParamsType);
			/*
            throw ModelException(method, "Null engine parameters found but "
                                 "require type: " + engineParamsType->getName());
			*/
        }

        return engineParams->getEngineParams(engineParamsType);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed when trying to retrieve credit "
                             "engine parameters of type: " +
                             engineParamsType->getName());
    }
}


/** Returns the number of "inner loss configs" contained in this LossConfig.
    In this case cascade down the call to the portfolio. */
int CreditTrancheLossConfig::numInnerLossConfigs() const {
    return portfolio->numInnerLossConfigs();
}

/** Returns the inner loss config */
ICreditLossConfigConstSP CreditTrancheLossConfig::getInnerLossConfig() const
{
    return portfolio;
}


/** Returns the inner loss config number "index".
    In this case cascade down the call to the portfolio. */
ICreditLossConfigConstSP CreditTrancheLossConfig::getInnerLossConfig(
    const int index) const
{
    return portfolio->getInnerLossConfig(index);
}

/** Returns the maximum (realistic) loss that may be produced by this
    loss config. Equivalent to lossGivenDefault for single names.
    "Realistic" in the sense that this method does NOT assume, eg,
    RR=0, since otherwise this method would be identical to notional() */
double CreditTrancheLossConfig::maxPossibleLoss() const {
    // Need to implement this properly - throw an exception by now
    throw ModelException("CreditTrancheLossConfig::maxPossibleLoss",
                         "Internal error: this method is not yet implemented");
}

/** Returns the loss config notional, ie, the amount that could be
    (or could have been) lost in the worst possible scenario, even if
    this scenario is no longer possible (e.g., if a name has defaulted
    with a RR > 0 it will still be accounted as if RR=0 if the losses
    are greater in this case).
    Note this notional bears no relation to the trade position */
double CreditTrancheLossConfig::notional() const {
    return highStrike - lowStrike;
}


/** Computes the loss and recovered notional currently produced by this
    loss config. Note both can be 0 no names have defaulted */
void CreditTrancheLossConfig::currentLossAndRecovered(
    double& loss,               // (O)
    double& recoveredNotional,  // (O)
    const bool recoverNotional) const
{
    loss = trancheLoss();
    recoveredNotional = (recoverNotional ?
                         trancheRecoveredNotional() : 0.0);
}

/** Computes how much of this loss config's original notional has been
    prepaid (amortized) before "prepay date" */
double CreditTrancheLossConfig::getPrepaidNotional(
    const DateTime& prepayDate) const
{
    // This method is not yet fully implemented!!
    return 0.0;
}
//
//-------- ICreditLossConfig methods


//++++++++++++++++++++++++++++++++++++++++
// Class-specific methods required for convolution purposes (non-virtual)

/** Stores the strikes of the tranche in the method arguments (references) */
void CreditTrancheLossConfig::getTrancheStrikes(double& k1, double& k2) const {
    k1 = lowStrike;
    k2 = highStrike;
}

void CreditTrancheLossConfig::getTrancheXSubAdjustedStrikes(
	double& k1,
	double& k2) const
{
	if (!isXSub) return getTrancheStrikes(k1,k2);

	// it's an inner tranche
	if (!isOuter)
	{
		k1 = 0;
		k2 = highStrike;
		return;
	}

	// outerTranche
	k1 = lowStrike - sumBufferForXSub;
	k2 = highStrike - sumBufferForXSub;
};

/** Utility method: calculates survival probabilities along timeline.
    Resizes survivalProb as required.
    CAUTION: Requires the names in the portfolio to be of type PortfolioName */
void CreditTrancheLossConfig::computeNameSurvivalProb(
    const DateTimeArray& timeline,           // (I)
    DoubleArrayArray&    survivalProb) const // (O) [time index][name index]
{
    static const string method("CreditTrancheLossConfig::computeNameSurvivalProb");

    try {
        if (!portf) {
            throw ModelException(method, "Internal error: The portfolio "
                                 "supplied is not of type CDOPortfolio, so this "
                                 "method cannot be called");
        }
        int numNames = numInnerLossConfigs(); // for ease
        // reserve some space
        int timelineSize = timeline.size();
        survivalProb.resize(timelineSize);

        for (int t=0; t < timelineSize; ++t) {
            survivalProb[t].resize(numNames);
        }

        for (int i=0; i < numNames; ++i) {
            if (nameDefaulted(i)) {
                for (int t=0; t < timelineSize; ++t) {
                    survivalProb[t][i] = 0.0; // override any input
                }
            }
            else {
                const DateTime& protectionStart =
                    portf->getName(i)->getProtectionStartDate();
                const DateTime& protectionEnd =
                    valueDate.max(portf->getName(i)->getProtectionEndDate(timeline.back()));

                // JLHP Here need to get hold of the ICreditLossGen for each of
                // the names - requires access to the model (or ModelMapper),
                // which needs to be passed in
                SingleCreditAssetConstSP asset(portf->nameAsset(i));
                DefaultRatesSP defaultRates(
                    asset->getParSpreadCurve()->defaultRates());

                for (int t=0; t < timelineSize; ++t) {
                    // could be improved by using ICDSParSpreads::Key and also
                    // be stopping once we go beyond protectionEnd
                    // Need to watch out for protection start dates>timeline[t]
                    const DateTime& thisProtectionEnd =
                        timeline[t].min(protectionEnd);

                    // Survival probabilities are UNCONDITIONAL probabilities,
					// from protection start to thisProtection end. This is
                    // needed for fwd tranche.
					// Note that this used to be CONDITIONAL probs in previous
                    // versions
					survivalProb[t][i] = (protectionStart <= thisProtectionEnd) ?
                        1.0 - (defaultRates->calcDefaultPV(valueDate, protectionStart)
                               - defaultRates->calcDefaultPV(valueDate, thisProtectionEnd)) :
                        1.0;
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// JLHP These methods below should NOT be here, they just do not fit
// into the new infrastrucutre...

/** Returns the notional for the name corresponding to the
    specified index which must lie in range [0, numInnerLossConfigs()-1] */
double CreditTrancheLossConfig::nameNotional(int index) const {
    return portfolio->getInnerLossConfig(index)->notional();
}

/** Returns the 'beta' for the name corresponding to the specified
    index which must lie in range [0, numNames()-1] */
double CreditTrancheLossConfig::nameBeta(int index) const {
    CreditEngineParametersConstSP engineParams =
        portfolio->getInnerLossConfig(index)->getEngineParams(CmOnlyParameters::TYPE);

    CmOnlyParametersConstSP cmParams(
        dynamic_cast<const CmOnlyParameters*>(engineParams.get()));

    return cmParams->getBeta();
}

/** Returns the 'decretion beta' for the name corresponding to the specified
    index which must lie in range [0, numNames()-1] */
double CreditTrancheLossConfig::nameDecBeta(int index) const {

    CreditEngineParametersConstSP engineParams =
        portfolio->getInnerLossConfig(index)->getEngineParams(CmOnlyParameters::TYPE);

    CmOnlyParametersConstSP cmParams(
        dynamic_cast<const CmOnlyParameters*>(engineParams.get()));

    return cmParams->getDecretionBeta();
}

/** Returns the recovery for the name corresponding to the
    specified index which must lie in range [0, numNames()-1]. Note that
    the recovery can be overridden at the trade level hence this method,
    which reflects any trade level overrides, should be used rather than
    the recovery off the par spread curve */
double CreditTrancheLossConfig::nameRecovery(int index) const {
    static const string method("CreditTrancheLossConfig::nameRecovery");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->getName(index)->getNameRecovery();
}

/** Returns true if the name has defaulted where the name corresponds to the
    specified index which must lie in range [0, numNames()-1] */
bool CreditTrancheLossConfig::nameDefaulted(int index) const {
    static const string method("CreditTrancheLossConfig::nameDefaulted");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->getName(index)->hasDefaulted();
}

/** Returns the ranges for Loss Given Default for the name
    corresponding to the specified index which must lie in range
    [0, numNames()-1]. In the particular, the LGD of a name is given by
    nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
void CreditTrancheLossConfig::nameLGDRanges(
    int     index,       /* (I) */
    double& lgdNotional, /* (O) */
    double& lgdFloor,    /* (O) */
    double& lgdCap)      /* (O) */ const
{
    static const string method("CreditTrancheLossConfig::nameLGDRanges");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    portf->getName(index)->
        nameLGDRanges(lgdNotional, lgdFloor, lgdCap);
}

/** Returns the max of valueDate and protection start date
    for specified name (if any). The specified index must lie in
    range [0, numNames()-1]*/
const DateTime& CreditTrancheLossConfig::nameProtectionStartDate(
    int index) const
{
    static const string method("CreditTrancheLossConfig::nameProtectionStartDate");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->getName(index)->getProtectionStartDate();
}

/** Returns the min of last date parameter and protection end date
    for specified name (if any). The specified index must lie in
    range [0, numNames()-1]. This if for name 'cutoff' */
const DateTime& CreditTrancheLossConfig::nameProtectionEndDate(
    int index,
    const DateTime& lastDate) const
{
    static const string method("CreditTrancheLossConfig::nameProtectionEndDate");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->getName(index)->getProtectionEndDate(lastDate);
}

/** Returns protection start date */
const DateTime& CreditTrancheLossConfig::getProtectionStartDate() const
{
    if (!protectionStartDate) {
        return valueDate;
    }
    else {
        return protectionStartDate->max(valueDate);
    }
}

/** Returns the 'loss decretion beta' */
double CreditTrancheLossConfig::lossDecretionBeta() const {
    static const string method("CreditTrancheLossConfig::lossDecretionBeta");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->lossDecretionBeta;
}

/** Returns the representation of the underlying name corresponding to the
    specified index which must lie in range [0, numNames()-1] */
SingleCreditAssetConstSP CreditTrancheLossConfig::nameAsset(int index) const {
    static const string method("CreditTrancheLossConfig::nameAsset");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->nameAsset(index);
}

/** Returns the total notional of the portfolio */
double CreditTrancheLossConfig::portfolioNotional() const {
    static const string method("CreditTrancheLossConfig::portfolioNotional");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->portfolioNotional();
}

/** Returns the total short notional of the portfolio */
double CreditTrancheLossConfig::portfolioLongNotional() const {
    static const string method("CreditTrancheLossConfig::portfolioLongNotional");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->portfolioLongNotional();
}

/** Returns the total short notional of the portfolio */
double CreditTrancheLossConfig::portfolioShortNotional() const {
    static const string method("CreditTrancheLossConfig::portfolioShortNotional");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->portfolioShortNotional();
}

/** Returns the TOTAL historical losses in the underlying portfolio,
    ie, before applying the tranche */
double CreditTrancheLossConfig::portfolioLoss() const {
    static const string method("CreditTrancheLossConfig::portfolioLoss");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->portfolioLoss();
}


/** Returns the sum of the past long recovered notional of the
    portfolio which ate into the total notional returned by
    portfolioLongNotional() */
double CreditTrancheLossConfig::portfolioLongRecoveredNotional() const {
    static const string method("CreditTrancheLossConfig::portfolioLongRecoveredNotional");
    if (!portf) {
        throw ModelException(method, "Internal error: The portfolio "
                             "supplied is not of type CDOPortfolio, so this "
                             "method cannot be called");
    }
    return portf->portfolioLongRecoveredNotional();
}


// jlhp This can (should) be private
/** Returns highStrike-lowStrike with highStrike and lowStrike obtained from
    trancheStrikes */
double CreditTrancheLossConfig::trancheNotional() const {
    return (highStrike-lowStrike);
}


/** Returns the outstanding notional of a tranche as seen from today */
double CreditTrancheLossConfig::trancheOutstandingNotional(
    const bool recoverNotional) const
{
    double bottomLoss = Maths::creditCollar(portfolioLoss(), lowStrike, highStrike);
        if (recoverNotional) {
        // need to consider recovered notional separately from losses to the
        // bottom of the tranche. Normally would write
        // collar(recoveredNotional, 1-high strike, 1-low strike)
        // where 1 = portfolio notional. But collar is invariant under
        // a transformation where a constant is added to each parameter.
        // We can use shiftedCollar rather than creditCollar since
        // "1-high strike" etc will be positive provided we enfore
        // high strike <= portfolio notional
        double topLoss = trancheRecoveredNotional();
        return highStrike - Maths::max(lowStrike,0.0) - bottomLoss - topLoss;
    }
    else {
        return highStrike - Maths::max(lowStrike,0.0) - bottomLoss;
    }
}


DateTime CreditTrancheLossConfig::getToday() const {
    return valueDate;
}


const double CreditTrancheLossConfig::trancheLoss() const {
    return Maths::creditCollar(portfolioLoss(), lowStrike, highStrike);
}


/** Normally would write
    collar(recoveredNotional, 1-high strike, 1-low strike)
    where 1 = portfolio notional. But collar is invariant under
    a transformation where a constant is added to each parameter.
    We can use shiftedCollar rather than creditCollar since
    "1-high strike" etc will be positive provided we enfore
    high strike <= portfolio notional */
const double CreditTrancheLossConfig::trancheRecoveredNotional() const {
    return Maths::shiftedCollar(
        -portfolioLongNotional()+ portfolioLongRecoveredNotional(),
        -highStrike,
        -lowStrike);
}
//
//----------------------------------------


// /** Sorts and aggregates the incoming "notional reductions" array, and
//  * returns a new cash flow array, with the actual reductions affecting
//  * this tranche as defined by the k1 and k2 (low and high) strikes */
// CashFlowArraySP CreditTrancheLossConfig::getTrancheReductions(
//     CashFlowArraySP reductions,
//     const double    k1,
//     const double    k2)
// {
//     if (reductions->empty()) {
//         return reductions; // Shortcut
//     }

//     // Sort reductions
//     sort(reductions->begin(), reductions->end(), CashFlow::lessThenForDates);
//     // and make every date unique by merging cashflows on same date
//     CashFlow::aggregate(*reductions);

//     CashFlowArraySP trancheReductions(new CashFlowArray(0));
//     double currentLoss = 0.0; // to avoid inserting the same loss more than once
//     double totalReduction = 0.0;
//     for (int j=0; j < reductions->size(); ++j) {
//         totalReduction += (*reductions)[j].amount;
//         double amount = Maths::creditCollar(totalReduction, k1, k2);

//         if (!Maths::equals(currentLoss, amount)) {
//             // Do not insert a new amount if equal to the previous loss
//             trancheReductions->push_back(CashFlow((*reductions)[j].date,
//                                                   amount));
//             currentLoss = amount;
//         }
//     }
//     return trancheReductions;
// }


/** Sorts and aggregates the incoming "portfolioLosses" array, and
    returns a array of CtgLegLossPerDefault after applying the tranche
    as defined by the k1 and k2 (low and high) strikes.
    CAUTION: the entries in the returned array are non-cumulative */
CtgLegLossPerDefaultArraySP CreditTrancheLossConfig::getTrancheReductions(
    CtgLegLossPerDefaultArraySP portfolioLosses,
    const double                k1,
    const double                k2)
{
    int numPortfolioLosses = portfolioLosses->size();
    if (numPortfolioLosses == 0) { // Shortcut for empty reductions array
        return portfolioLosses;
    }

    sort(portfolioLosses->begin(),
         portfolioLosses->end(),
         CtgLegLossPerDefault::lessThanForLossDates);

    CtgLegLossPerDefaultArraySP trancheLosses(new CtgLegLossPerDefaultArray(0));
    double currentLoss = 0.0; // to avoid inserting the same loss more than once
    double totalReduction = 0.0;
    for (int i=0; i < numPortfolioLosses; ++i) {
        totalReduction += (*portfolioLosses)[i]->loss.amount;
        double amount = Maths::creditCollar(totalReduction, k1, k2);

        if (!Maths::equals(currentLoss, amount)) {
            // Do not insert a new amount if equal to the previous loss
            trancheLosses->push_back(CtgLegLossPerDefaultSP(
                new CtgLegLossPerDefault((*portfolioLosses)[i]->defaultDate,
                                         (*portfolioLosses)[i]->loss.date,
                                         amount - currentLoss)));
            currentLoss = amount;
        }
    }
    return trancheLosses;
}


/** Sorts and aggregates the incoming "portfolioLosses" array, and
    returns a array of FeeLegReductionPerDefault after applying the tranche
    as defined by the k1 and k2 (low and high) strikes.
    CAUTION: the entries in the returned array are non-cumulative */
FeeLegReductionPerDefaultArraySP CreditTrancheLossConfig::getTrancheReductions(
    FeeLegReductionPerDefaultArraySP portfolioLosses,
    const double                     k1,
    const double                     k2,
    const bool                       isLoss)
{
    int numPortfolioLosses = portfolioLosses->size();
    if (numPortfolioLosses == 0) { // Shortcut for empty reductions array
        return portfolioLosses;
    }

    sort(portfolioLosses->begin(),
         portfolioLosses->end(),
         FeeLegReductionPerDefault::lessThanForLossDates);

    FeeLegReductionPerDefaultArraySP trancheLosses(
        new FeeLegReductionPerDefaultArray(0));
    double currentReduction = 0.0; // to avoid inserting the same amount more than once
    double totalReduction = 0.0;
    for (int i=0; i < numPortfolioLosses; ++i) {
        totalReduction += (*portfolioLosses)[i]->getReductionAmount(isLoss);
        double amount = Maths::creditCollar(totalReduction, k1, k2);

        if (!Maths::equals(currentReduction, amount)) {
            // Do not insert a new amount if equal to the previous one
            trancheLosses->push_back(FeeLegReductionPerDefaultSP(
                new FeeLegReductionPerDefault((*(*portfolioLosses)[i]),
                                              isLoss,
                                              amount - currentReduction)));
            currentReduction = amount;
        }
    }
    return trancheLosses;
}


/** Theta::IShift method */
bool CreditTrancheLossConfig::sensShift(Theta* shift) {
    // alter immediate data
    valueDate = shift->rollDate(valueDate);
    return true; // then shift components
}


/** EffectiveCurveLossModelConfig::IIntoLossGen method:
    Create an IEffectiveLossCurveGen as specified by the supplied
    EffectiveLossCurveModelConfig */
IEffectiveCurveLossGenSP CreditTrancheLossConfig::lossGenerator(
    IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig) const
{
    return effCurveLossModelConfig->createLossGenerator(
        CreditTrancheLossConfigConstSP::attachToRef(this));
}


void CreditTrancheLossConfig::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditTrancheLossConfig, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICreditLossConfig);
    IMPLEMENTS(IEffectiveCurveLossModelConfig::IIntoLossGen);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(IRestorableWithRespectTo<CDOParallelStrike>);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name, "Name of the tranche config");
    FIELD(valueDate, "Value date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(lowStrike, "Low strike (absolute value, NOT a percentage)");
    FIELD(highStrike, "High strike (absolute value, NOT a percentage)");
    FIELD(portfolio, "Underlying loss config to apply the tranche to");
    FIELD(engineParams, "Engine parameters");
    FIELD_MAKE_OPTIONAL(engineParams);
    FIELD(isXSub, "Is this a cross-sub CDO^2. Default: False");
    FIELD_MAKE_OPTIONAL(isXSub);
	FIELD_NO_DESC(isOuter);
	FIELD_MAKE_TRANSIENT(isOuter);
	FIELD(cutOffDate, "Risk cut-off for the tranche");
    FIELD_MAKE_OPTIONAL(cutOffDate);
	FIELD(protectionStartDate, "Protection start date");
	FIELD_MAKE_OPTIONAL(protectionStartDate);
}

IObject* CreditTrancheLossConfig::defaultConstructor() {
    return new CreditTrancheLossConfig();
}

CClassConstSP const CreditTrancheLossConfig::TYPE =
    CClass::registerClassLoadMethod("CreditTrancheLossConfig",
                                    typeid(CreditTrancheLossConfig),
                                    load);

DEFINE_TEMPLATE_TYPE(CreditTrancheLossConfigWrapper);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(CreditTrancheLossConfigArray);

bool CreditTrancheLossConfigLinkIn() {
    return CreditTrancheLossConfig::TYPE != NULL;
}

// #########################################################################
/** StateVariable Generator corresponding to the CreditTrancheLossConfig  */
class CreditTrancheLossConfig::SVGen:
	public virtual ICreditLossConfig::ISVGen,
	public virtual IElemStateVariableGen
{
public:
	/** Nested inner class: State variables representing loss event times  */
	class SV : 	public virtual ICreditLossConfig::ISVGen::ISV
	{
	public:
		friend class CreditTrancheLossConfig::SVGen;

		/** constructor */
		SV(
			double lowStrike,
			double highStrike,
			bool recoverNotional,
			double componentNotional,
			const ICreditLossConfigSVSP& componentLoss,
			const DateTimeLiteVectorConstSP& timeline):
		    //maybe we need to pass the pastPathGenerator for the past information
			componentLoss(componentLoss),
			lowStrike(lowStrike),
			highStrike(highStrike),
			recoverNotional(recoverNotional),
			componentNotional(componentNotional),
			timeline(timeline)
		{
			notional = (highStrike - lowStrike);
			if (timeline.get())
				creditLossConfigSVResults.resize(timeline->size());
			else
				throw ModelException("CreditTrancheLossConfig::SVGen::SV::SV",
				"timeline not provided");
		}

		/** virtual destructor */
		virtual ~SV() {}

		/** read access method for the results */
		virtual const CreditLossConfigSVResultList& getResults() const
		{
			this->calculateLossEvents();
			return creditLossConfigSVResults;
		}

		/** read/write access method for the results */
		virtual CreditLossConfigSVResultList& results()
		{
			this->calculateLossEvents();
			return creditLossConfigSVResults;
		}

		/** access method for the product timeline */
		virtual const DateTimeLiteVectorConstSP& getTimeline() const
		{
			return timeline;
		}

		/** access method for the component Notional */
		virtual double getNotional() const
		{
			return notional;
		}

		/** method of the parent interface, IStateVariable */
		virtual bool doingPast() const
		{
			return false;
		}

	private:
		/** Just the declaration of the copy constructor; to avoid compiler creating them   */
		SV(const SV&);

		/** Just the declaration of the equal to operator; to avoid compiler creating them   */
		SV& operator=(const SV&);

		/** aggregates lossEvents of underlying components and populates creditLossConfigSVResults */
		void calculateLossEvents() const
		{
			try
			{
				creditLossConfigSVResults.clear();

//The following assumes that the product timeline for this object and the components are the same. This is ensured by the constructor.
//sub-components would have elements only if losses have occured there
//first collect the losses
				CreditLossConfigSVResultList& ithComponentLosses = componentLoss->results();

//scale down the nominal amounts at this stage
				double scalingFactor = componentNotional/componentLoss->getNotional();

				for (CreditLossConfigSVResultList::iterator iterx = ithComponentLosses.begin();
						iterx != ithComponentLosses.end(); ++iterx)
				{
					iterx->lossAmount *=  scalingFactor;
					iterx->notionalAmount *=  scalingFactor;
				}

//note I am modifying the results of the sub-components and they will get corrupted
//sub-components would have elements only if losses have occured there
				creditLossConfigSVResults.insert(creditLossConfigSVResults.end(), ithComponentLosses.begin(), ithComponentLosses.end());

//sort them on eventTime
				creditLossConfigSVResults.sort();

				if (creditLossConfigSVResults.size() == 0)
					return;

//add and remove duplicates
				CreditLossConfigSVResultList::iterator lastIter =  creditLossConfigSVResults.begin();
				double lastAggLoss = lastIter->lossAmount;
				double lastTranchedLoss = Maths::creditCollar(lastAggLoss, lowStrike, highStrike);
				double lastAggNotional = lastIter->notionalAmount;
				double lastAggPoolRecovery = lastAggNotional - lastAggLoss;
				double lastTranchedRecovery = Maths::max(lastAggPoolRecovery - (componentNotional - highStrike), 0.) -
					Maths::max(lastAggPoolRecovery - (componentNotional - lowStrike), 0.);
				double lastTranchedNotional = recoverNotional ? (lastTranchedLoss + lastTranchedRecovery):lastTranchedLoss;

//update notional later
				CreditLossConfigSVResultList::iterator iter =  creditLossConfigSVResults.begin();
				double aggLoss;
				double tranchedLoss;
				double aggNotional;
				double aggPoolRecovery;
				double tranchedRecovery;
				double tranchedNotional;
				iter++;

				while (iter != creditLossConfigSVResults.end())
				{
//aggregate if the event times are the same
					aggLoss = lastAggLoss + iter->lossAmount;
					tranchedLoss = Maths::creditCollar(aggLoss, lowStrike, highStrike);
					aggNotional = lastAggNotional + iter->notionalAmount;
					aggPoolRecovery = aggNotional - aggLoss;
					tranchedRecovery = Maths::max(aggPoolRecovery - (componentNotional - highStrike), 0.) -
						Maths::max(aggPoolRecovery - (componentNotional - lowStrike), 0.);
					// min(max(R - (N-E),0) , E-A )
					tranchedNotional = recoverNotional ? (tranchedLoss + tranchedRecovery):tranchedLoss;

					if (lastIter->eventTime ==  iter->eventTime)
					{
//delete the object at the iter
						creditLossConfigSVResults.erase(iter);
//update the iter
						iter = lastIter;
					}

					if (iter == creditLossConfigSVResults.begin())
					{
						iter->lossAmount = tranchedLoss;
						iter->notionalAmount = tranchedNotional;
					}
					else
					{
						iter->lossAmount = tranchedLoss - lastTranchedLoss;
						iter->notionalAmount = tranchedNotional - lastTranchedNotional;
					}

					if (iter != creditLossConfigSVResults.end()) //we do this because we may have deleted an item
					{
						iter++;
						lastAggLoss = aggLoss;
						lastTranchedLoss = tranchedLoss;
						lastAggNotional = aggNotional;
						lastAggPoolRecovery = aggPoolRecovery;
						lastTranchedRecovery = tranchedRecovery;
						lastTranchedNotional = tranchedNotional;
					}
				}
			}
			catch(exception& e)
			{
				throw ModelException(e, "CDOPortfolioSVGen::SV", "calculateLossEvents");
			}
		}

// ##  Member variables below
		/** LossEvent state variables for the component   */
		ICreditLossConfigSVSP componentLoss;

		/** lower strike of the tranche */
		double lowStrike;

		/** upper strike of the tranche */
		double highStrike;

		/** recovery notional for the tranche */
		bool recoverNotional;

		/** pool notional for the tranche */
		double componentNotional;

		/** pool notional for the tranche */
		double notional;

		/** storage for the creditLossConfig  */
		mutable CreditLossConfigSVResultList	 creditLossConfigSVResults;

		/** storage for the timeline  */
		DateTimeLiteVectorConstSP timeline;

	}; // CreditTrancheLossConfig::SVGen::SV

	typedef smartPtr<CreditTrancheLossConfig::SVGen::SV> CreditTrancheLossConfigSVSP;

/* ## Methods of CreditTrancheLossConfig::SVGen  */

	/** main constructor */
    SVGen(const CreditTrancheLossConfigConstSP& tranche,
          const DateTimeLiteVectorConstSP& timeline,
          const CIntConstSP& triggerDelay,
          const CIntConstSP& defaultToCalculationDelay,
          double temporaryLossAmount,
          const DateTime& lastTriggerDate,
          const AccrualPeriodArrayConstSP& accrualPeriods,
          const IBadDayAdjusterConstSP& bda,
          const IProtectionProviderConstSP& protect,
          const IRebateCalculatorConstSP& rebateCalc,
          const bool recoverNotional) :
		recoverNotional(recoverNotional),
		componentNotional(tranche->portfolioNotional()),
		timeline(timeline),
        triggerDelay(triggerDelay),
		defaultToCalculationDelay(defaultToCalculationDelay),
		temporaryLossAmount(temporaryLossAmount),
		lastTriggerDate(lastTriggerDate),
		accrualPeriods(accrualPeriods),
		bda(bda),
		protect(protect),
		rebateCalc(rebateCalc)
	{
		tranche->getTrancheXSubAdjustedStrikes(lowStrike, highStrike);
		// tranche->getTrancheStrikes(lowStrike, highStrike); JCP

		//first get the loss config
		ICreditLossConfigConstSP creditLossConfig = tranche->getPortfolio();

		//check whether it can support MC
		const ICreditLossConfigSVGenMC* creditLossConfigSVGenMC =
			dynamic_cast<const ICreditLossConfigSVGenMC*>(creditLossConfig.get());

		if (creditLossConfigSVGenMC == 0)
			throw ModelException("CreditTrancheLossConfig::SVGen", "constructor");

		//create the state variable gen
		componentLossGen = creditLossConfigSVGenMC->createSVGen(
			timeline,
			triggerDelay,
			defaultToCalculationDelay,
			temporaryLossAmount,
			lastTriggerDate,
			accrualPeriods,
			bda,
			protect,
			rebateCalc,
            recoverNotional);
	}

	/** virtual destructor */
	virtual ~SVGen() {}

	/** Actually creates and returns a new instance of the SV */
	virtual ICreditLossConfigSVSP createNewSV(IStateVariableGen::IStateGen* stateGen) const
	{
		try
		{
			ICreditLossConfigSVSP componentSV = componentLossGen->createNewSV(stateGen);
			return ICreditLossConfigSVSP(
				new CreditTrancheLossConfig::SVGen::SV(
					lowStrike,
					highStrike,
					recoverNotional,
					componentNotional,
					componentSV,
					timeline));
		}
		catch(exception& e)
		{
			throw ModelException(e,"CreditTrancheLossConfig::SVGen","createNewSV");
		}
	}

	/** Same as create but avoids dynamic cast */
	virtual ICreditLossConfigSVSP getLossSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const
	{
		try
		{
			//later make PathGenCrEvTimes implement IStateVariableGen::IStateGen
			CreditTrancheLossConfigSVSP creditTrancheLossConfigSV =
				CreditTrancheLossConfigSVSP(dynamic_cast<CreditTrancheLossConfig::SVGen::SV*>(oldStateVar.get()));

			//now we refresh the components inside this with a recursive call
			creditTrancheLossConfigSV->componentLoss =
				this->componentLossGen->getLossSV(creditTrancheLossConfigSV->componentLoss, stateGen);

			return ICreditLossConfigSVSP(creditTrancheLossConfigSV);
		}
		catch(exception& e)
		{
			throw ModelException(e,"CreditTrancheLossConfig::SVGen","getLossSV");
		}
	}


	// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const
	{
		return this->getLossSV(oldStateVar, stateGen);
	}

// ##  Methods of the parent interface, IStateVariableClient
    virtual void collectStateVars(IStateVariableCollectorSP svdb) const
	{
		componentLossGen->collectStateVars(svdb);
	}

// ##  Methods of the parent interface, IElemStateVariableGen
	virtual void attachSVGen(class IElemStateVariableGenVisitor* sv) const
	{
		//sv->processSVGen(this);
	}

private:
	/** Just the declaration of the copy constructor; to avoid compiler creating them   */
	SVGen(const SVGen&);

	/** Just the declaration of the equal to operator; to avoid compiler creating them   */
	SVGen& operator=(const SVGen&);

	// ##  Member variables follow
	/** lower strike of the tranche */
	double lowStrike;

	/** upper strike of the tranche */
	double highStrike;

	/** recovery notional for the tranche */
	bool recoverNotional;

	/** LossEvent state variables generators for the components   */
	ICreditLossConfigSVGenConstSP componentLossGen;

	/** notionals of the component SV inside   */
	double componentNotional;

	/** storage for the timeline  */
	DateTimeLiteVectorConstSP timeline;

	/** Number of days gap between the Credit Event Date (CED) and the Event Determination Date (EDD)
		To be used for calculating Past Default Losses by calling methods on PortfolioName */
	CIntConstSP triggerDelay;

	/**  Number of days gap between the CED and the Calculation Date (CD)
	*/
	CIntConstSP defaultToCalculationDelay;

	/**  Loss amount to be used when actual recovery amount is not known
	*/
	const double temporaryLossAmount;

	/**  Date by which the default must be triggered
	*/
	const DateTime& lastTriggerDate;

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

// ################################################################################################################
/** Indexed StateVariable Generator corresponding to the CreditTrancheLossConfig  */

/** StateVariable Generator corresponding to the CreditTrancheLossConfig  */
class MARKET_DLL CreditTrancheLossConfig::IndexedSVGen:	public virtual ICreditLossConfig::IIndexedSVGen,
														public virtual IElemStateVariableGen
{
public:
	/** Nested inner class: State variables representing indexed loss event times  */
	class SV : 	public virtual ICreditLossConfig::IIndexedSVGen::ISV
	{
	public:
		friend class CreditTrancheLossConfig::IndexedSVGen;

		/** constructor */
		SV(
			double lowStrike,
			double highStrike,
			bool recoverNotional,
			double componentNotional,
			int protectionStartIndex,
			int cutOffIndex,
			const ICreditLossConfigIndexedSVSP& componentLoss): //maybe we need to pass the pastPathGenerator for the past information
			componentLoss(componentLoss),
			lowStrike(lowStrike),
			highStrike(highStrike),
			recoverNotional(recoverNotional),
			componentNotional(componentNotional),
			protectionStartIndex(protectionStartIndex),
			cutOffIndex(cutOffIndex),
			resultPathInitialized(false)
		{
			notional = (highStrike - lowStrike);
		}

		/** virtual destructor */
		virtual ~SV() {};

		/** read access method for the results */
		virtual const CreditLossConfigIndexedSVResultPath& getResults() const
		{
			this->calculateLossEvents();
			return creditLossConfigIndexedSVResultPath;
		}

		/** read/write access method for the results */
		virtual CreditLossConfigIndexedSVResultPath& results()
		{
			this->calculateLossEvents();
			return creditLossConfigIndexedSVResultPath;
		}

		/** access method for the component Notional */
		virtual double getNotional() const
		{
			return notional;
		}

		/** method of the parent interface, IStateVariable */
		virtual bool doingPast() const
		{
			return false;
		}

	private:

		/** returns the size of the resultPath */
		virtual int getResultPathSize() const
		{
			if (resultPathInitialized)
				return creditLossConfigIndexedSVResultVector.size();
			else
				throw ModelException("CreditTrancheLossConfig::getResultPathSize", "ResultPath not initialized");
		}

		/** Just the declaration of the copy constructor; to avoid compiler creating them   */
		SV(const SV&);

		/** Just the declaration of the equal to operator; to avoid compiler creating them   */
		SV& operator=(const SV&);

		/** initializes the ResultPath object. Gets called by calculateLossEvents() */
		void initializeResultPath(int size) const
		{
			creditLossConfigIndexedSVResultVector.resize(size);


			// JCP this was size -1
			creditLossConfigIndexedSVResultPath.initialize(&creditLossConfigIndexedSVResultVector[0], 0, size-1);

			resultPathInitialized = true;
		}


		/** aggregates lossEvents of underlying components and populates  creditLossConfigSVResults */
		void calculateLossEvents() const
		{
			try
			{
				if (componentLoss.get() == 0)
					return;

//The following assumes that the engine timeline for this object and the component are the same. This is ensured by the constructor.
				const CreditLossConfigIndexedSVResultPath& componentResultPath = componentLoss->getResults();

//	we obtain the size from the results of the components.
//  we also assume that, size of the path will not change with the simulation runs
				if (!resultPathInitialized)
					this->initializeResultPath(componentLoss->getResultPathSize());

//first set losses and notionals to zero
				for (unsigned int i=0; i< creditLossConfigIndexedSVResultVector.size(); ++i)
				{
					creditLossConfigIndexedSVResultVector[i].lossChange = 0;
					creditLossConfigIndexedSVResultVector[i].notionalChange = 0;
					creditLossConfigIndexedSVResultVector[i].index = -1;
				}

//scale down the nominal amounts for the subcomponents at this stage
				double scalingFactor = componentNotional/componentLoss->getNotional();

//tranche the results
				double aggregateLoss = 0;
				double lastTranchedLoss = 0;
				double tranchedLoss = 0;

				double aggregateRecovery = 0;
				double tranchRecovery = 0;

				double aggregateNotional = 0;
				double lastTranchedNotional = (highStrike - lowStrike);
				double tranchedNotional = 0;

				for (unsigned int j = 0; j < creditLossConfigIndexedSVResultVector.size(); ++j)
				{
					if (  (componentResultPath[j].index != -1) && (componentResultPath[j].index <= cutOffIndex) &&
						(componentResultPath[j].index >= protectionStartIndex) )
					{
						creditLossConfigIndexedSVResultVector[j].index = componentResultPath[j].index;

						aggregateLoss += componentResultPath[j].lossChange*scalingFactor;
						tranchedLoss = Maths::creditCollar(aggregateLoss, lowStrike, highStrike);
						creditLossConfigIndexedSVResultVector[j].lossChange = tranchedLoss - lastTranchedLoss; //since we store the deltas

						aggregateNotional += componentResultPath[j].notionalChange*scalingFactor;

						aggregateRecovery = aggregateNotional - aggregateLoss;

						tranchRecovery = Maths::shiftedCollar( - componentNotional + aggregateRecovery,
															   - highStrike,
															   - lowStrike);

						tranchedNotional = (highStrike - lowStrike) +
							(recoverNotional ? (- tranchedLoss - tranchRecovery): - tranchedLoss); //tranchedNotional

						creditLossConfigIndexedSVResultVector[j].notionalChange = lastTranchedNotional - tranchedNotional; //since we store the deltas

						lastTranchedLoss = tranchedLoss;
						lastTranchedNotional = tranchedNotional;
					}
				}
			}
			catch(exception& e)
			{
				throw ModelException(e, "CDOPortfolioIndexedSVGen", "calculateLossEvents");
			}
		}

// ##  Member variables below
		/** LossEvent state variables for the component   */
		ICreditLossConfigIndexedSVSP componentLoss;

		/** lower strike of the tranche */
		double lowStrike;

		/** upper strike of the tranche */
		double highStrike;

		/** recovery notional for the tranche */
		bool recoverNotional;

		/** pool notional for the tranche */
		double componentNotional;

		/** pool notional for the tranche */
		double notional;

		/** Index on the productTimeline that corresponds to the protection start date of the tranche.
		*/
		int protectionStartIndex;

		/** Index on the productTimeline that corresponds to the maturity cutoff of the tranche.
			If there is no maturity cutoff, then it is the last point on the product timeline
		*/
		int cutOffIndex;

		/** storage for the CreditLossConfigSV Results  */
		mutable CreditLossConfigIndexedSVResultVector creditLossConfigIndexedSVResultVector;

		/** GenericPath representing the value of sv results along a monte carlo path */
		mutable CreditLossConfigIndexedSVResultPath	 creditLossConfigIndexedSVResultPath;

		/** storage for the timeline  */
		DateTimeLiteVectorConstSP timeline;

		/** indicates whether the member "creditLossConfigIndexedSVResultVector" has been initialized  */
		mutable bool resultPathInitialized;
	};

	typedef smartPtr<IndexedSVGen::SV> CreditTrancheIndexedLossConfigSVSP;

	/* ## Methods of IndexedSVGen  */
	/** main constructor */
	IndexedSVGen(const CreditTrancheLossConfigConstSP& tranche,
                 const DateTimeArrayConstSP& timeline,
                 const CIntConstSP& triggerDelay,
                 const CIntConstSP& defaultToCalculationDelay,
                 double temporaryLossAmount,
                 const DateTime& lastTriggerDate,
                 const AccrualPeriodArrayConstSP& accrualPeriods,
                 const IBadDayAdjusterConstSP& bda,
                 const IProtectionProviderConstSP& protect,
                 const IRebateCalculatorConstSP& rebateCalc,
                 const bool recoverNotional) :
		recoverNotional(recoverNotional),
		componentNotional(tranche->portfolioNotional()),//- ASarma - this call needs revisiting
		triggerDelay(triggerDelay),
		defaultToCalculationDelay(defaultToCalculationDelay),
		temporaryLossAmount(temporaryLossAmount),
		lastTriggerDate(lastTriggerDate),
		accrualPeriods(accrualPeriods),
		bda(bda),
		protect(protect),
		rebateCalc(rebateCalc)
	{
		tranche->getTrancheXSubAdjustedStrikes(lowStrike, highStrike);

	    // first get the loss config
		ICreditLossConfigConstSP creditLossConfig = tranche->getPortfolio();

		//compute the cut-off date
		if (tranche->hasCutOffFlag())
		{
			try
			{
				const DateTime& cutOffDate = tranche->getCutOffDate();
				cutOffIndex = cutOffDate.find(*timeline);
			}
			catch(exception& e)
			{
				throw ModelException(
					e,
					"CreditTrancheLossConfig::IndexedSVGen::constructor",
					"CutOffDate not found on the product timeline");
			}
		}
		else
			cutOffIndex = timeline->size() - 1;

		//compute protection start date index
		try
		{
			const DateTime& protectionStartDate = tranche->getProtectionStartDate();
			protectionStartIndex = protectionStartDate.find(*timeline);
		}
		catch(exception& e)
		{
			throw ModelException(
				e,
				"CreditTrancheLossConfig::IndexedSVGen::constructor",
				"ProtectionStartDate not found on the product timeline");
		}

	    // check whether it can support MC
		const ICreditLossConfigIndexedSVGenMC* creditLossConfigIndexedSVGenMC =
			dynamic_cast<const ICreditLossConfigIndexedSVGenMC*>(creditLossConfig.get());

		if (creditLossConfigIndexedSVGenMC == 0)
			throw ModelException("IndexedSVGen", "constructor");

	    // create the state variable gen for the nested component
		componentLossGen =
            creditLossConfigIndexedSVGenMC->createIndexedSVGen(
                timeline,
                triggerDelay,
                defaultToCalculationDelay,
                temporaryLossAmount,
                lastTriggerDate,
                accrualPeriods,
                bda,
                protect,
                rebateCalc,
                recoverNotional);
	}

	/** virtual destructor */
	virtual ~IndexedSVGen() {}

	/** Actually creates and returns a new instance of the SV */
	virtual ICreditLossConfigIndexedSVSP createNewSV(
		IStateVariableGen::IStateGen* stateGen) const
	{
		try
		{
			ICreditLossConfigIndexedSVSP componentSV = componentLossGen->createNewSV(stateGen);

			return ICreditLossConfigIndexedSVSP(
				new SV(
					lowStrike,
					highStrike,
					recoverNotional,
					componentNotional,
					protectionStartIndex,
					cutOffIndex,
					componentSV));
		}
		catch(exception& e)
		{
			throw ModelException(e,"IndexedSVGen","createNewSV");
		}
	}

	/** Same as create but avoids dynamic cast */
	virtual ICreditLossConfigIndexedSVSP getLossSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const
	{
		try
		{
//later make PathGenCrEvTimes implement IStateVariableGen::IStateGen
			CreditTrancheIndexedLossConfigSVSP creditTrancheIndexedLossConfigSV =
				CreditTrancheIndexedLossConfigSVSP(dynamic_cast<SV*>(oldStateVar.get()));

//now we refresh the components inside this with a recursive call
			creditTrancheIndexedLossConfigSV->componentLoss =
				this->componentLossGen->getLossSV(creditTrancheIndexedLossConfigSV->componentLoss, stateGen);

			return ICreditLossConfigIndexedSVSP(creditTrancheIndexedLossConfigSV);
		}
		catch(exception& e)
		{
			throw ModelException(e,"IndexedSVGen","getLossSV");
		}
	}

// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const
	{
		return this->getLossSV(oldStateVar, stateGen);
	}

// ##  Methods of the parent interface, IStateVariableClient
    virtual void collectStateVars(IStateVariableCollectorSP svdb) const
	{
		componentLossGen->collectStateVars(svdb);
	}

// ##  Methods of the parent interface, IElemStateVariableGen
	virtual void attachSVGen(class IElemStateVariableGenVisitor* sv) const
	{
		//sv->processSVGen(this);
	}
private:
	/** Just the declaration of the copy constructor; to avoid compiler creating them   */
	IndexedSVGen(const IndexedSVGen&);

	/** Just the declaration of the equal to operator; to avoid compiler creating them   */
	IndexedSVGen& operator=(const IndexedSVGen&);

// ##  Member variables follow
	/** lower strike of the tranche */
	double lowStrike;

	/** upper strike of the tranche */
	double highStrike;

	/** recovery notional for the tranche */
	bool recoverNotional;

	/** LossEvent state variables generators for the components   */
	ICreditLossConfigIndexedSVGenConstSP componentLossGen;

	/** notionals of the component SV inside   */
	double componentNotional;

	/** Index on the productTimeline that corresponds to the protection start date of the tranche.
	*/
	int protectionStartIndex;

	/** Index on the productTimeline that corresponds to the maturity cutoff of the tranche.
		If there is no maturity cutoff, then it is the last point on the product timeline
	*/
	int cutOffIndex;

	/** Number of days gap between the Credit Event Date (CED) and the Event Determination Date (EDD)
		To be used for calculating Past Default Losses by calling methods on PortfolioName */
	CIntConstSP triggerDelay;

	/**  Number of days gap between the CED and the Calculation Date (CD)
	*/
	CIntConstSP defaultToCalculationDelay;

	/**  Loss amount to be used when actual recovery amount is not known
	*/
	const double temporaryLossAmount;

	/**  Date by which the default must be triggered
	*/
	const DateTime& lastTriggerDate;

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


// ################################################################################################################

//# Methods of the interface, ICreditTrancheLossConfigSVGenMC follow

/** Creates the corresponding state variable  */
ICreditLossConfigSVGenConstSP
CreditTrancheLossConfig::createSVGen(
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
	) const
{
	CreditTrancheLossConfigConstSP creditTrancheLossConfigSP(this);
	return ICreditLossConfigSVGenConstSP(
		new CreditTrancheLossConfig::SVGen(
			creditTrancheLossConfigSP,
			timeline,
			triggerDelay,
			defaultToCalculationDelay,
			temporaryLossAmount,
			lastTriggerDate,
			accrualPeriods,
			bda,
			protect,
			rebateCalc,
            recoverNotional));
}


//# Methods of the interface, ICreditTrancheLossConfigIndexedSVGenMC follow

ICreditLossConfigIndexedSVGenConstSP
CreditTrancheLossConfig::createIndexedSVGen(
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
	) const
{
	CreditTrancheLossConfigConstSP creditTrancheLossConfigSP(this);
	return ICreditLossConfigIndexedSVGenConstSP(
		new CreditTrancheLossConfig::IndexedSVGen(
			creditTrancheLossConfigSP,
			timeline,
			triggerDelay,
			defaultToCalculationDelay,
			temporaryLossAmount,
			lastTriggerDate,
			accrualPeriods,
			bda,
			protect,
			rebateCalc,
            recoverNotional));
}

/** [Implements ITranchesCombinationPayoff] */
void CreditTrancheLossConfig::linearDecomposition(
    const DateTime& time,
    DoubleArray& baseStrikes,        /* output */
    DoubleArray& baseStrikesWeights, /* output */
    double& expectedLossOffset       /* output */) const
{
    expectedLossOffset = 0.0;
    if (lowStrike == 0.0)
    {
        baseStrikes.resize(1);
        baseStrikesWeights.resize(1);
        baseStrikes[0] = highStrike;
        baseStrikesWeights[0] = 1.0;
    }
    else
    {
        baseStrikes.resize(2);
        baseStrikesWeights.resize(2);
        baseStrikes[0] = lowStrike;
        baseStrikes[1] = highStrike;
        baseStrikesWeights[0] = -1.0;
        baseStrikesWeights[1] = 1.0;
    }
}

DRLIB_END_NAMESPACE
