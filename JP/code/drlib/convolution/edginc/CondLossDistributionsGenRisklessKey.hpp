//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 27-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONDLOSSDISTRIBUTIONSGENRISKLESSKEY_HPP
#define QLIB_CONDLOSSDISTRIBUTIONSGENRISKLESSKEY_HPP

#include "edginc/ICondLossDistributionsGen.hpp"

DRLIB_BEGIN_NAMESPACE

/** 
 * Trivial implementation of ICondLossDistributionsGenKey for
 * riskless assets: always returns a loss distribution with no losses and
 * a survival probability of 1.
 * */
class CondLossDistributionsGenRisklessKey:
    public CObject,
    public virtual ICondLossDistributionsGenKey
{
public:
    
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Constructor */
    CondLossDistributionsGenRisklessKey();
    
    /** [Implements ICondLossDistributionsGenKey] */
    virtual IDistribution1DConstSP conditionalLossDistribution(
        IMarketFactorValueConstSP marketFactorValue) const;

    virtual double conditionalSurvProb(
        IMarketFactorValueConstSP marketFactorValue) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    DiscreteDistributionSP condLossDistribution;
};

DRLIB_END_NAMESPACE

#endif /*QLIB_CONDLOSSDISTRIBUTIONSGENRISKLESSKEY_HPP*/
