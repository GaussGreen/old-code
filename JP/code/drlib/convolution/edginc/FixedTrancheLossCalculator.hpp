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

#ifndef QR_FIXEDTRANCHELOSSCALCULATOR_HPP
#define QR_FIXEDTRANCHELOSSCALCULATOR_HPP
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tranche Loss Calculator interface.
 *
 * Classes which implement this interface are capable of returning the
 * expected loss of a tranche at a series of timepoints with a fixed
 * lower strike and an upper strike. The timepoints and the strikes
 * are implicitly assumed to be specified at construction time.
 */
class CONVOLUTION_DLL IFixedTrancheLossCalculator: public virtual VirtualDestructorBase{
public:
    virtual ~IFixedTrancheLossCalculator();// in CreditMetricsLossCalculatorBase

    /** Calculate the expected loss for specified timepoint. */
    virtual void loss(
        int     timePoint,          // (I) do the calculation for this timepoint
        double& loss,                   /* (O) tranche loss amt  */
        double& lossCond) const = 0;    /* (O) tranche loss amt cond on cpty surviving */
    
    /** store calculator specific results */
    virtual void storeResults(Results* result, Control* control) const = 0;
};

typedef smartConstPtr<
    IFixedTrancheLossCalculator> IFixedTrancheLossCalculatorConstSP;

DRLIB_END_NAMESPACE
#endif
