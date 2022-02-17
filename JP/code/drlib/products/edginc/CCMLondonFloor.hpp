//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMLondonFloor.hpp
//
//   Description : London floor adjustment algorithm
//
//   Date        : Dec 2004
//
//----------------------------------------------------------------------------
#ifndef EDR_CCM_LONDON_FLOOR_HPP
#define EDR_CCM_LONDON_FLOOR_HPP

#include "edginc/CreditMetricsLondonFloor.hpp"

DRLIB_BEGIN_NAMESPACE

/** 
 * CCM model with London Floor
 * */
class PRODUCTS_DLL CCMLondonFloor : public CreditMetricsLondonFloor {
public:

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~CCMLondonFloor();

    /** Overridden to apply a 'London Floor' adjustment */
    virtual CreditMetricsLossCalculatorBase* createLossCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** Overridden to apply a 'London Floor' adjustment */
    virtual CreditMetricsLossCalculatorBase* createRecoveredNotionalCalculatorBase(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    /** Indicates whether this model supports stochastic recovery rates
        AND there are any engine parameters for any names specifying so. */
    virtual const bool hasStochasticRecoveries(
        CreditTrancheLossConfigConstSP tranche) const;

private:
    /** Only build instances of that class using reflection */
    CCMLondonFloor(const CClassConstSP& clazz);
    
    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE

#endif
