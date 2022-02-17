//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_DISCRETEDISTRIBUTION_HPP
#define QLIB_DISCRETEDISTRIBUTION_HPP

#include "edginc/Object.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/IDistribution1D.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Representation of a discrete distribution i.e.
 * essentially a collection of (value, probability) points sorted by
 * increasing "value".
 * */
class UTIL_DLL DiscreteDistribution:
    public CObject,
    public virtual IDistribution1D
{
public:    
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~DiscreteDistribution();

    /** Internal constructor */
    DiscreteDistribution(
        DoubleArraySP values,
        DoubleArraySP probabilities);
        
    /** Constructor for binary distributions {(0,1-prob); (value, prob)}*/
    DiscreteDistribution(
        double value,
        double probability);
        
    /** Access values */
    DoubleArrayConstSP getValues() const;

    /** Access probabilities */
    DoubleArrayConstSP getProbabilities() const;
    
    /** Update probability */
    void setProbability(int index, double probability);

     /**
     * Bisection search for "x" in field "values", using a "tolerance" such that:
     * - return "true" if found
     *   (in that case output "index" is such that |x - (*values)[index])| < tolerance
     * - return "false" if not found
     *   (in that case output "index" is such that (*values)[index-1] < x - tolerance and
     *    x + tolerance < (*values)[index])
     * 
     * NB: this method internally checks that tolerance < delta / 2 where delta is 
     *     the minimum difference between 2 consecutive values
     *     (this ensures we can find at MOST one x such that |x - (*values)[index])| < tolerance)
     * */
    bool search(double x, int& index, double tolerance = 0.0) const;

    /**
     * Probability density function
     * [implements IDistribution1D]
     * */
    virtual double pdf(double x) const;
    
    /**
     * Cumulative density function
     * [implements IDistribution1D]
     * */
    virtual double cdf(double x) const;
    
    /**
     * nth derivative of "moment generating function" mgf(z)=E[exp(zX)]
     * [implements IDistribution1D]
     * */
    virtual double mgf(double z, int n) const;
    
    /**
     * Returns the expectation value
     * [implements IDistribution1D]
     * */
    virtual double expectation() const;

    /**
     * Returns the expectation value
     * [implements IDistribution1D]
     * */
    virtual double expectation(const Function1DDouble& payoff) const;


    /**
     * Returns the variance
     * [implements IDistribution1D]
     * */
    virtual double variance() const;
    
    /**
     * Returns a "discretised" version of itself
     * [implements IDistribution1D]
     * */
    virtual DiscreteDistributionConstSP discretise() const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

private:
    DiscreteDistribution();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    DoubleArraySP values; // $required
    DoubleArraySP probabilities; // $required
    
    /** delta = minimum difference between 2 consecutive values */ 
    double delta;
};

DECLARE(DiscreteDistribution);

DRLIB_END_NAMESPACE

#endif

