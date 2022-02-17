//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IRebateCalculator.hpp
//
//   Description : Interface to an object capable of computing rebate payments,
//                 i.e., the difference in fee payments under different loss
//                 assumptions
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef IREBATECALCULATOR_HPP
#define IREBATECALCULATOR_HPP

#include "edginc/CashFlow.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE


/** Interface for objects capable of computing rebate payments, i.e., the 
 * difference in fee payments under different loss assumptions */
class MARKET_DLL IRebateCalculator : public virtual IObject {
public:
    static CClassConstSP const TYPE;
    virtual ~IRebateCalculator();

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
        IForwardRatePricerSP model) const = 0;

protected:
    IRebateCalculator();

private:
    static void load(CClassSP& clazz);
};

DECLARE(IRebateCalculator);

DRLIB_END_NAMESPACE

#endif
