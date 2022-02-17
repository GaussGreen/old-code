//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 09-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_GAUSSIANMARKETFACTORMODEL_HPP
#define QLIB_GAUSSIANMARKETFACTORMODEL_HPP

#include "edginc/IMarketFactorModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * 1-Dimensional gaussian market factor
 * */
class CONVOLUTION_DLL GaussianMarketFactorModel:
    public CObject,
    public virtual IMarketFactorModel
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /**
     * Density function f:R^N->R.
     * N is the dimension of the market factor.
     * We allow time-dependent market factors.
     * */
     // TODO - review returned type, may want to have directly
     //        "double density(const DateTime& time, const DoubleArray& marketFactorValue)"
     //        for performance
    virtual FunctionNDDoubleConstSP density(const DateTime& time) const;
    
    /** Some products can only be priced (closed-form) if the density is 
     * constant over time. Need a method to determine whether this model's
     * density depends on time or not */
    virtual bool isDensityTimeDependent() const;

    /** Dimension N of the market factor */
    // NB: could get N through density(...) method but seems a bit heavy
    virtual int dimension() const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();
    
    /** Overrides clone() method to copy field "gaussianDensity" */
    virtual IObject* clone() const;
    
    /** Public constructor */
    GaussianMarketFactorModel();

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    class GaussianDensity; // define in source file

    // default values for optional fields
    static const double LOWER_BOUND;
    static const double UPPER_BOUND;

    RangeSP range; // $optional
    FunctionNDDoubleSP gaussianDensity; // $unregistered
};

DRLIB_END_NAMESPACE

#endif /*GAUSSIANMARKETFACTORMODEL*/
