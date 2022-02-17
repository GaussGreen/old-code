//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Class containing all required information regarding a   
//                 reduction in notional in the fee leg (typically for CDOs).
//                 The reduction can be due to a loss or recovered notional.
//                 The information captured here shall permit the computation
//                 of the rebate payments corresponding to the reductions.
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


FeeLegReductionPerDefault::FeeLegReductionPerDefault(
        const DateTime& determinationDate, 
        const DateTime& effectiveDate,
        const DateTime& calculationDate,
        const double defaultedNotional,
        const double totalNameNotional,
        const double recoveryRate) :
    determinationDate(determinationDate),
    effectiveDate(effectiveDate),
    calculationDate(calculationDate),
    defaultedNotional(defaultedNotional),
    totalNameNotional(totalNameNotional),
    recoveryRate(recoveryRate)
{}

FeeLegReductionPerDefault::FeeLegReductionPerDefault(
        const FeeLegReductionPerDefault& reduction, 
        const bool isLoss,
        const double amount) :
    determinationDate(reduction.determinationDate),
    effectiveDate(reduction.effectiveDate),
    calculationDate(reduction.calculationDate),
    recoveryRate(reduction.recoveryRate)
{
    if (isLoss) {
        defaultedNotional = amount / (1.0 - recoveryRate);
        totalNameNotional = amount; // so that recovered amount is 0
    }
    else {
        defaultedNotional = 0.0; // so that loss is 0;
        totalNameNotional = amount;
    }
}


const double FeeLegReductionPerDefault::getLossAmount() const {
    const double loss = defaultedNotional * (1.0 - recoveryRate);
    return loss;
}

const double FeeLegReductionPerDefault::getRecoveredAmount() const {
    return totalNameNotional - getLossAmount();
}

const double FeeLegReductionPerDefault::getReductionAmount(
    const bool wantLoss) const 
{
    return wantLoss ? getLossAmount() : getRecoveredAmount();
}


/** Returns a CashFlowArray with all the reduction CashFlows (loss or
    recovered notional) inside the argument array of reductions */
CashFlowArraySP FeeLegReductionPerDefault::getReductions(
    const FeeLegReductionPerDefaultArrayConstSP x,
    const bool wantLoss)
{
    CashFlowArraySP reductions(new CashFlowArray(0));

    // Add all reductions in the x array, merging them
    int numReductions = !x ? 0 : x->size();
    for (int i=0; i < numReductions; ++i) {
        CashFlow newCf((*x)[i]->effectiveDate,
                       (*x)[i]->getReductionAmount(wantLoss));

        reductions = CashFlow::merge(
            reductions, CashFlowArraySP(new CashFlowArray(1, newCf)));
    }
    
    // Aggregate all cashflows happening on the same date
    if (!!reductions) {
        CashFlow::aggregate(*reductions);
    }
    return reductions;
}


bool FeeLegReductionPerDefault::lessThanForLossDates(
    const FeeLegReductionPerDefaultSP x, 
    const FeeLegReductionPerDefaultSP y)
{
    if (!x || !y) {
        return true; // One of them is null, so this is arbitrary
    }
    return (x->effectiveDate == y->effectiveDate ?
            x->calculationDate < y->calculationDate : // Equal eff dates so use calc dates
            x->effectiveDate < y->effectiveDate);
}

DRLIB_END_NAMESPACE
