//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CCM.hpp
//
//   Description : simple closed form model for CCM
//
//   Author      : Antoine Gregoire
//
//   Date        : 5 Nov 2004
//
//----------------------------------------------------------------------------

#ifndef EDR_CCM_HPP
#define EDR_CCM_HPP
#include "edginc/CreditMetricsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/** Composite Compula Model, it also drives how the quanto adjustment is done
    for cds par spreads */
class PRODUCTS_DLL CCM: public CreditMetricsModel {
public:
    friend class CCMHelper;
    
    static CClassConstSP const TYPE;
    
    virtual ~CCM();

    /** Indicates whether this model supports stochastic recovery rates
        AND there are any engine parameters for any names specifying so. */
    virtual const bool hasStochasticRecoveries(
        CreditTrancheLossConfigConstSP tranche) const;

    /** Static method that indicates whether the tranche contains names with
        engine parameters specifying stochastic recovery rates */
    static const bool ccmHasStochasticRecoveries(
        CreditTrancheLossConfigConstSP tranche);

protected:
    /** Simple constructor */
    CCM(const CClassConstSP& clazz);

    /** Same as CreditMetricsModel::createLossCalculator() method but
        returns one which is derived from
        CreditMetricsLossCalculatorBase. This method won't make sense
        for certain derived types - or rather its meaning changes.
        The implementation here creates an appropriate CCM loss
        calculator */
    virtual CreditMetricsLossCalculatorBase* createLossCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    virtual CreditMetricsLossCalculatorBase* createRecoveredNotionalCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

private:
    /** Don't use copy constructor */
    CCM(const CCM &rhs);
    
    /** Don't use */
    CCM& operator=(const CCM& rhs);
};

DRLIB_END_NAMESPACE
#endif



