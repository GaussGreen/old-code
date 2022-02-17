//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CCMRFL.hpp
//
//   Description : Model Combination of RFL and CCM
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_CCMRFL_HPP
#define EDR_CCMRFL_HPP

#include "edginc/CreditMetricsRFL.hpp"
#include "edginc/CCMLossCalculatorBase.hpp"


DRLIB_BEGIN_NAMESPACE

class ITrancheLossCalculatorLegacy;

/** Composite Compula Model, it also drives how the quanto adjustment is done
    for cds par spreads */
class CCMRFL: public CreditMetricsRFL {

public:
    friend class CCMHelper;
    
    static CClassConstSP const TYPE;

    /** Indicates whether this model supports stochastic recovery rates
        AND there are any engine parameters for any names specifying so. */
    virtual const bool hasStochasticRecoveries(
        CreditTrancheLossConfigConstSP tranche) const;

protected:
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

    virtual CreditMetricsLossCalculatorBase* 
    createRecoveredNotionalCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** Simple constructor */
    CCMRFL(const CClassConstSP& clazz);

    /** Destructor */
    virtual ~CCMRFL();

private:

    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** For CreditMetricsRFL::IIntoProduct */
    static void loadIntoProduct(CClassSP& clazz);

    /** Don't use copy constructor */
    CCMRFL(const CCMRFL &rhs);
    
    /** Don't use */
    CCMRFL& operator=(const CCMRFL& rhs);
};

DRLIB_END_NAMESPACE
#endif



