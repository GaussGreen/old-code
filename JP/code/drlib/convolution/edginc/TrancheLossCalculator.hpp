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

#ifndef QR_TRANCHELOSSCALCULATOR_HPP
#define QR_TRANCHELOSSCALCULATOR_HPP
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE
class IFixedTrancheLossCalculator;
class ITrancheLossCalculator;
typedef smartConstPtr<ITrancheLossCalculator> ITrancheLossCalculatorConstSP;

/**
 * Tranche Loss Calculator interface.
 *
 * Classes which implement this interface are capable of pricing a tranche at
 * a series of timepoints with a lower strike and an upper strike. The 
 * timepoints are implicitly assumed to be specified at construction time.
 */
class CONVOLUTION_DLL ITrancheLossCalculator: public virtual VirtualDestructorBase{
public:
    ITrancheLossCalculator();
    virtual ~ITrancheLossCalculator();  // in CreditMetricsLossCalculatorBase

    /** Calculate the expected loss for specified timepoint and strikes. 
        Note if repeated calculations at the same timepoint are required then
        IKey::calc should be used instead */
    virtual void loss(
        int     timePoint,          // (I) do the calculation for this timepoint
        double  k1,                   /* (I) lower strike      */
        double  k2,                   /* (I) upper strike      */
        double& loss,                 /* (O) tranche loss amt  */
        double& lossCond) const = 0;  /* (O) tranche loss amt cond on cpty
                                         surviving */

    virtual void storeResults(Results* result, Control* control) const = 0;

    /** Interface to optimise repeated calculations of losses 
        for different strikes */
    class CONVOLUTION_DLL IKey: public virtual VirtualDestructorBase{
        IKey(const IKey& rhs);
        IKey& operator=(const IKey& rhs);
    public:
        IKey();
        virtual ~IKey();  // in CreditMetricsLossCalculatorBase.cpp

        /** Calculate the expected loss for specified strikes. */
        virtual void loss(
            double  k1,                   /* (I) lower strike      */
            double  k2,                   /* (I) upper strike      */
            double& loss,                 /* (O) tranche loss amt  */
            double& lossCond) const = 0;  /* (O) tranche loss amt cond on cpty
                                             surviving */
    };
    typedef smartPtr<IKey> IKeySP;

    /** Returns a key used to optimise repeated calculations of losses
        at the same timepoint. */
    virtual IKey* lossKey(
        int timePoint) const = 0; // (I) do the calculation for this timepoint

    /** Utility method to build a IFixedTrancheLossCalculator using a
        ITrancheLossCalculator. Implementation in
        CreditMetricsLossCalculatorBase.cpp */
    static IFixedTrancheLossCalculator* createFixedTrancheLossCalculator(
        double                        lowerStrike,  
        double                        upperStrike,  
        ITrancheLossCalculatorConstSP lossCalculator);
private:
    ITrancheLossCalculator(const ITrancheLossCalculator& rhs);
    ITrancheLossCalculator& operator=(const ITrancheLossCalculator& rhs);
    class FixedStrikes;
};

typedef smartConstPtr<ITrancheLossCalculator> ITrancheLossCalculatorConstSP;

DRLIB_END_NAMESPACE
#endif
