//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//
//----------------------------------------------------------------------------

#ifndef QR_TRANCHELOSSCALCULATORLEGACY_HPP
#define QR_TRANCHELOSSCALCULATORLEGACY_HPP
#include "edginc/smartPtr.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tranche Loss Calculator interface.
 *
 * Classes which implement this interface are capable of pricing a tranche
 * with a lower strike and an upper strike.
 *
 * It is not clear where this belongs - it is an ability of any model
 * which can price a CDO tranche and is not limited to the CCM model.
 * So at the point where a new model is introduced for tranche pricing,
 * this interface should be moved out of a CCM-specific file.
 */
class CONVOLUTION_DLL ITrancheLossCalculatorLegacy
{
public:
    ITrancheLossCalculatorLegacy();//in CreditMetricsLossCalculatorBase
    virtual ~ITrancheLossCalculatorLegacy();//in CreditMetricsLossCalculatorBase

    // the loss function - this needs to be defined in the sub-classes
    virtual void loss(
        double        K1,    /* (I) lower strike      */
        double        K2,    /* (I) upper strike      */
        double       &L,     /* (O) tranche loss amt  */
        double       &Lcond) /* (O) tranche loss amt cond on cpty surviving */
        const = 0;
private:
    ITrancheLossCalculatorLegacy(const ITrancheLossCalculatorLegacy& rhs);
    ITrancheLossCalculatorLegacy& operator=(const ITrancheLossCalculatorLegacy& rhs);
};
typedef refCountPtr<
    ITrancheLossCalculatorLegacy> ITrancheLossCalculatorLegacySP;

DRLIB_END_NAMESPACE
#endif
