#include "edginc/config.hpp"
#include "edginc/IDistribution1D.hpp"
#include "edginc/Class.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

static void IDistribution1DLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IDistribution1D, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IDistribution1D::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IDistribution1D", 
        typeid(IDistribution1D), 
        IDistribution1DLoad);

DEFINE_TEMPLATE_TYPE(IDistribution1DArray);

/**
 * Addin for unit testing of distribution methods.
 * Inputs are the distribution, x, z and n.
 * Result is the following array R:
 * R[0] = pdf(x)
 * R[1] = cdf(x)
 * R[2] = mgf(z,n)
 * R[3] = expectation()
 * R[4] = variance()
 * */
class UTIL_DLL DistributionUnitTest:
    public CObject,
    public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // Addin parameters
    IDistribution1DSP distribution;
    double x;
    double z;
    int n;

    // EdrAction version of addin
    IObjectSP run() {
        DoubleArraySP result(new DoubleArray(5));
        (*result)[0] = distribution->pdf(x);
        (*result)[1] = distribution->cdf(x);
        (*result)[2] = distribution->mgf(z,n);
        (*result)[3] = distribution->expectation();
        (*result)[4] = distribution->variance();

        return result;
    }

    DistributionUnitTest(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DistributionUnitTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(distribution, "The distribution to test");
        FIELD(x, "Value used in pdf(x) and cdf(x)");
        FIELD(z, "Value used in mgf(z,n)");
        FIELD(n, "Value used in mgf(z,n)");
    }
    
    static IObject* defaultConstructor(){
        return new DistributionUnitTest();
    }
};

CClassConstSP const DistributionUnitTest::TYPE =
    CClass::registerClassLoadMethod(
        "DistributionUnitTest",
        typeid(DistributionUnitTest),
        load);

/* external symbol to allow class to be forced to be linked in */
bool IDistribution1DLoad(){
    return (DistributionUnitTest::TYPE != 0 && IDistribution1D::TYPE != 0);
}

DRLIB_END_NAMESPACE
