#include "edginc/config.hpp"
#include "edginc/FunctionIntegratorND.hpp"
#include "edginc/Format.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/FunctionOperations.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(FunctionIntegratorNDArray);

/** TYPE (for reflection) */        
CClassConstSP const FunctionIntegratorND::TYPE =
CClass::registerClassLoadMethod(
    "FunctionIntegratorND",
    typeid(FunctionIntegratorND),
    FunctionIntegratorND::load);

/** Virtual destructor */
FunctionIntegratorND::~FunctionIntegratorND()
{
}

/** Constructor */
FunctionIntegratorND::FunctionIntegratorND():CObject(TYPE) {}

IObject* FunctionIntegratorND::defaultConstructor()
{
    return new FunctionIntegratorND();
}

void FunctionIntegratorND::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(FunctionIntegratorND, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IIntegrator);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(integratorND, "Type of ND integrator to use");
}

/** Integration method */
IObjectSP FunctionIntegratorND::integrate(const IIntegrand* integrand) const
{
    const MFunctionND* f = dynamic_cast<const MFunctionND*>(integrand);
    if (f == 0) {
        throw ModelException(
            "FunctionIntegratorND::integrate",
            "Unable to cast integrand into MFunctionND");
    }
    if (f->getNbFuncs() != 1)
    {
        throw ModelException(
            "FunctionIntegratorND::integrate",
            "Invalid integrand: nb functions = " +
            Format::toString(f->getNbFuncs()) +
            " (1 expected).");
    }

    FunctionNDDoubleConstSP function = DYNAMIC_POINTER_CAST<FunctionNDDouble>(
        FunctionOperations::project(*f, 0));
    
    return CDoubleSP(CDouble::create(integratorND->integrate(*function)));
}

/** Addin for unit testing of ND integrator engines */
class UTIL_DLL IntegratorNDUnitTest:
    public CObject,
    public ClientRunnable {
public:
    
    // inner class for exp function
    class ExpFunction: public virtual FunctionNDDouble {
    public:
        virtual double operator()(const CDoubleArray&  x) const {
            double sum = 0.0;
            for (int i = 0; i < x.size(); ++i) {
				sum +=x[i];
			}
            return exp(sum);
        }
        
        ExpFunction(int N, const RangeArray& intervals):
            FunctionNDDouble(N, intervals){}
	};

    static CClassConstSP const TYPE;

    // EdrAction version of addin
    IObjectSP run() {
        int N = lowerBounds->size();
        RangeArray ranges(N);
        for (int i = 0; i < N; ++i) {
			ranges[i] = RangeSP( new Range((*lowerBounds)[i], true, (*upperBounds)[i], true));
		}
        if (integrand == EXP)
        {
            ExpFunction f(N, ranges);
            return integrator->integrate(&f);
        }
        else
        {
            throw ModelException(
                "IntegratorNDUnitTest::run",
                "Integrand " + integrand + " is not supported");
        }
    }
    
    virtual void validatePop2Object(){
        static const string method("IntegratorNDUnitTest::validatePop2Object");        
        try
        {
            // checks "lowerBounds" and "upperBounds" have the same size
            if (lowerBounds->size() != upperBounds->size())
            {
                throw ModelException(
                    method,
                    "Fields 'lowerBounds' and 'upperBounds' don't have "
                    "the same number of elements.");            
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

private:
    IntegratorNDUnitTest(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IntegratorNDUnitTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(integrator, "The ND integrator to test");
        FIELD(integrand,
            "Function to integrate. "
            "Possible values are: " +
            EXP + " (x->exp(x))");
        FIELD(lowerBounds, "Lower integration bounds");
        FIELD(upperBounds, "Upper integration bounds");
    }
    
    static IObject* defaultConstructor(){
        return new IntegratorNDUnitTest();
    }

    // Addin parameters
    FunctionIntegratorNDSP integrator;
    string integrand;
    DoubleArraySP lowerBounds;
    DoubleArraySP upperBounds;

    // String for basic functions
    static const string EXP;
};

const string IntegratorNDUnitTest::EXP = "EXP";

CClassConstSP const IntegratorNDUnitTest::TYPE =
    CClass::registerClassLoadMethod(
        "IntegratorNDUnitTest",
        typeid(IntegratorNDUnitTest),
        load);

/* external symbol to allow class to be forced to be linked in */
bool FunctionIntegratorNDLoad(){
    return (IntegratorNDUnitTest::TYPE != 0);
}

DRLIB_END_NAMESPACE
