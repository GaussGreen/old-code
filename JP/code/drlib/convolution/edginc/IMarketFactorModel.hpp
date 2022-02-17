//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 03-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IMARKETFACTORMODEL_HPP
#define QLIB_IMARKETFACTORMODEL_HPP

#include "edginc/Function.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface to define the market factor "shape" (distribution, random
 * variable generation...)
 * */
class CONVOLUTION_DLL IMarketFactorModel: public virtual IObject {
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
    virtual FunctionNDDoubleConstSP density(const DateTime& time) const = 0;
    
    /** Some products can only be priced (closed-form) if the density is 
     * constant over time. Need a method to determine whether this model's
     * density depends on time or not */
    virtual bool isDensityTimeDependent() const = 0;

    /** Dimension N of the market factor */
    // NB: could get N through density(...) method but seems a bit heavy
    virtual int dimension() const = 0;
    
    // potentially needs other fields / methods for MC
};

DECLARE(IMarketFactorModel);

DRLIB_END_NAMESPACE

#endif /*QLIB_IMARKETFACTORMODEL_HPP*/
