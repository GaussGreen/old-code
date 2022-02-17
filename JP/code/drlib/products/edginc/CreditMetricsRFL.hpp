//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CreditMetricsRFL.hpp
//
//   Description : Credit Metrics with RFL
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDIT_METRICS_RFL_HPP
#define EDR_CREDIT_METRICS_RFL_HPP

#include "edginc/CreditMetricsModel.hpp"

DRLIB_BEGIN_NAMESPACE
/** 
 * CreditMetrics model with RFL
 * */
class CreditMetricsRFL : public CreditMetricsModel {
public:

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~CreditMetricsRFL();

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

    /** Called immediately after object constructed */
    //virtual void validatePop2Object();

protected:
    /** Only build instances of that class using reflection */
    CreditMetricsRFL(const CClassConstSP& clazz);
    
private:
    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** For CreditMetricsRFL::IIntoProduct */
    static void loadIntoProduct(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif



