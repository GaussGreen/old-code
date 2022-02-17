#include "edginc/config.hpp"
#include "edginc/GaussianMarketFactorModel.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
#include ext_hash_map
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

/** TYPE (for reflection) */        
CClassConstSP const GaussianMarketFactorModel::TYPE =
CClass::registerClassLoadMethod(
    "GaussianMarketFactorModel",
    typeid(GaussianMarketFactorModel),
    GaussianMarketFactorModel::load);

const double GaussianMarketFactorModel::LOWER_BOUND = -7.575;
const double GaussianMarketFactorModel::UPPER_BOUND = 7.575;

/** Constructor */
GaussianMarketFactorModel::GaussianMarketFactorModel():
    CObject(TYPE),
    range(new Range(LOWER_BOUND, true, UPPER_BOUND, true)){}

IObject* GaussianMarketFactorModel::defaultConstructor()
{
    return new GaussianMarketFactorModel();
}

void GaussianMarketFactorModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(GaussianMarketFactorModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMarketFactorModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(range,
        "Integration bounds for the market factor "
        "(Default is [" +
        Format::toString(LOWER_BOUND) +
        ", " +
        Format::toString(UPPER_BOUND) +
        "]).");
    FIELD_MAKE_OPTIONAL(range);
}

/**
 * Gaussian density - function returned by GaussianMarketFactorModel::density
 * */
class GaussianMarketFactorModel::GaussianDensity: public virtual Function1DDouble {
public:    
    // Constructor
    GaussianDensity(const Range& interval): Function1DDouble(interval){}
    
    // Function
	virtual double operator()(double  x) const
    {
        return N1Density(x);
    }
};

/** Called immediately after object constructed */
void GaussianMarketFactorModel::validatePop2Object()
{
    gaussianDensity.reset(new GaussianDensity(*range));
}

/** Overrides clone() method to copy const field "gaussianDensity" */
IObject* GaussianMarketFactorModel::clone() const{
    IObject* myCopy = CObject::clone();
    GaussianMarketFactorModel& gaussianMF = 
        dynamic_cast<GaussianMarketFactorModel&>(*myCopy);
    gaussianMF.gaussianDensity = gaussianDensity;
    return myCopy;
}

/**
 * Density function f:R^N->R.
 * N is the dimension of the market factor.
 * We allow time-dependent market factors.
 * */
 // TODO - review returned type, may want to have directly
 //        "double density(const DateTime& time, const DoubleArray& marketFactorValue)"
 //        for performance
FunctionNDDoubleConstSP GaussianMarketFactorModel::density(const DateTime& time) const
{
    return gaussianDensity;
}

/** Some products can only be priced (closed-form) if the density is 
 * constant over time. Need a method to determine whether this model's
 * density depends on time or not */
bool GaussianMarketFactorModel::isDensityTimeDependent() const {
    return false;
}

/** Dimension N of the market factor */
// NB: could get N through density(...) method but seems a bit heavy
int GaussianMarketFactorModel::dimension() const
{
    return 1;
}

/* external symbol to allow class to be forced to be linked in */
bool GaussianMarketFactorModelLoad(){
    return (GaussianMarketFactorModel::TYPE != 0);
}

DRLIB_END_NAMESPACE
