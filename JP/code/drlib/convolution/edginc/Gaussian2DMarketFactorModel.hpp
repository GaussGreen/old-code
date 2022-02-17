//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 06-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_GAUSSIAN2DMARKETFACTORMODEL_HPP
#define QLIB_GAUSSIAN2DMARKETFACTORMODEL_HPP

#include "edginc/IMarketFactorModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * 1-Dimensional gaussian market factor
 * */
class CONVOLUTION_DLL Gaussian2DMarketFactorModel:
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
    
    /** Public empty constructor */
    Gaussian2DMarketFactorModel();

    /** Public constructor */
    Gaussian2DMarketFactorModel(double rho, RangeArray &rangeArray);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    class Gaussian2DDensity; // define in source file

    // default values for optional fields
    static const double LOWER_BOUND_D1;
    static const double UPPER_BOUND_D1;

    static const double LOWER_BOUND_D2;
    static const double UPPER_BOUND_D2;

    RangeArray rangeArray; // $optional
    double rho;            // $optional
    FunctionNDDoubleSP gaussian2DDensity; // $unregistered
};

DRLIB_END_NAMESPACE

#endif /*GAUSSIAN2DMARKETFACTORMODEL*/
