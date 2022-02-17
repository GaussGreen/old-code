//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CCMBaseCorrelation.hpp
//
//   Description : Composite Copula Model with Base Correlation
//
//   Author      : Antoine Gregoire
//
//   Date        : May 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CCM_BASE_CORRELATION_HPP
#define QLIB_CCM_BASE_CORRELATION_HPP

#include "edginc/CreditMetricsBaseCorrelation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Composite Copula Model with Base Correlation */
class PRODUCTS_DLL CCMBaseCorrelation : public CreditMetricsBaseCorrelation {
public:

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~CCMBaseCorrelation() {}

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

    virtual CreditMetricsLossCalculatorBase* createRecoveredNotionalCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

	/** Only build instances of that class using reflection */
    CCMBaseCorrelation(const CClassConstSP& clazz);
    
    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE

#endif //QLIB_CCM_BASE_CORRELATION_HPP


