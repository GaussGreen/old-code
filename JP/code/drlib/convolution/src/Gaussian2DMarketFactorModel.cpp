#include "edginc/config.hpp"
#include "edginc/Gaussian2DMarketFactorModel.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
//#include ext_hash_map
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

/** TYPE (for reflection) */        
CClassConstSP const Gaussian2DMarketFactorModel::TYPE =
CClass::registerClassLoadMethod(
    "Gaussian2DMarketFactorModel",
    typeid(Gaussian2DMarketFactorModel),
    Gaussian2DMarketFactorModel::load);

const double Gaussian2DMarketFactorModel::LOWER_BOUND_D1 = -7.575;
const double Gaussian2DMarketFactorModel::UPPER_BOUND_D1 = 7.575;
const double Gaussian2DMarketFactorModel::LOWER_BOUND_D2 = -7.575;
const double Gaussian2DMarketFactorModel::UPPER_BOUND_D2 = 7.575;

/** Constructor */
Gaussian2DMarketFactorModel::Gaussian2DMarketFactorModel(double rho, RangeArray &rangeArray):
CObject(TYPE),
rangeArray(rangeArray),
rho(rho)
{}

/** Empty Constructor */
Gaussian2DMarketFactorModel::Gaussian2DMarketFactorModel():
CObject(TYPE),
rho(0),
rangeArray(RangeArray(2))
{
    //RangeArray rangeArray(2);
    rangeArray[0] = RangeSP(new Range(OpenBoundary(LOWER_BOUND_D1), OpenBoundary(UPPER_BOUND_D1)));
    rangeArray[1] = RangeSP(new Range(OpenBoundary(LOWER_BOUND_D2), OpenBoundary(UPPER_BOUND_D2)));
}

IObject* Gaussian2DMarketFactorModel::defaultConstructor()
{
    return new Gaussian2DMarketFactorModel();
}

void Gaussian2DMarketFactorModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(Gaussian2DMarketFactorModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMarketFactorModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(rho, "correlation of the 2 normal market factors");
    FIELD_MAKE_OPTIONAL(rho);
    FIELD(rangeArray,
        "Integration bounds for the market factor "
        "(Default is [" +
        Format::toString(LOWER_BOUND_D1) +
        ", " +
        Format::toString(UPPER_BOUND_D1) +
        "] x [" +
        Format::toString(LOWER_BOUND_D2) +
        ", " +
        Format::toString(UPPER_BOUND_D2) +
        "] ).");
    FIELD_MAKE_OPTIONAL(rangeArray);
}

/**
 * Gaussian density - function returned by GaussianMarketFactorModel::density
 * */
class Gaussian2DMarketFactorModel::Gaussian2DDensity: public virtual FunctionNDDouble {
public:    
    // Constructor
    Gaussian2DDensity(double rho, const RangeArray& intervals):
        rho(rho),
        FunctionNDDouble(2, intervals){}
    
    // Function
	virtual double operator()(const CDoubleArray&  x) const
    {
        return N2Density(x[0], x[1], rho);
    }
private:
    double rho;
};

/** Called immediately after object constructed */
void Gaussian2DMarketFactorModel::validatePop2Object()
{
    gaussian2DDensity.reset(new Gaussian2DDensity(rho, rangeArray));
}

/** Overrides clone() method to copy const field "gaussianDensity" */
IObject* Gaussian2DMarketFactorModel::clone() const{
    IObject* myCopy = CObject::clone();
    Gaussian2DMarketFactorModel& gaussianMF = 
        dynamic_cast<Gaussian2DMarketFactorModel&>(*myCopy);
    gaussianMF.gaussian2DDensity = gaussian2DDensity;
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
FunctionNDDoubleConstSP Gaussian2DMarketFactorModel::density(const DateTime& time) const
{
    return gaussian2DDensity;
}

/** Some products can only be priced (closed-form) if the density is 
 * constant over time. Need a method to determine whether this model's
 * density depends on time or not */
bool Gaussian2DMarketFactorModel::isDensityTimeDependent() const {
    return false;
}

/** Dimension N of the market factor */
// NB: could get N through density(...) method but seems a bit heavy
int Gaussian2DMarketFactorModel::dimension() const
{
    return 2;
}

/* external symbol to allow class to be forced to be linked in */
bool Gaussian2DMarketFactorModelLoad(){
    return (Gaussian2DMarketFactorModel::TYPE != 0);
}

DRLIB_END_NAMESPACE
