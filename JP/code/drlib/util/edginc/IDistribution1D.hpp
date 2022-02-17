//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IDISTRIBUTION1D_HPP
#define QLIB_IDISTRIBUTION1D_HPP

#include "edginc/Object.hpp"
#include "edginc/Function.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(DiscreteDistribution);
/**
 * Interface to represent the distribution of a one-dimensional
 * random variable X.
 * 
 * NB: we do not constrain probabilities to sum up to 1.
 * 
 * */
class UTIL_DLL IDistribution1D: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Probability density function */
    virtual double pdf(double x) const = 0;
    
    /** Cumulative density function */
    virtual double cdf(double x) const = 0;
    
    /** nth derivative of "moment generating function" mgf(z)=E[exp(zX)] */
    virtual double mgf(double z, int n) const = 0;
    
    /** Returns the expectation value of X */
    virtual double expectation() const = 0;

    /** Returns the expectation value of f(X) */
    virtual double expectation(const Function1DDouble& payoff) const = 0;

    /** Returns the variance of X */
    virtual double variance() const = 0;
    
    /** Returns a "discretised" version of itself */
    virtual DiscreteDistributionConstSP discretise() const = 0;
};

DECLARE(IDistribution1D);

DRLIB_END_NAMESPACE

#endif /*QLIB_IDISTRIBUTION1D_HPP*/
