#include "edginc/config.hpp"
#include "edginc/FunctionIntegrator1D.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/FunctionOperations.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(FunctionIntegrator1DArray);

/** TYPE (for reflection) */        
CClassConstSP const FunctionIntegrator1D::TYPE =
CClass::registerClassLoadMethod(
    "FunctionIntegrator1D",
    typeid(FunctionIntegrator1D),
    FunctionIntegrator1D::load);

/** Virtual destructor */
FunctionIntegrator1D::~FunctionIntegrator1D()
{
}

/** Constructor */
FunctionIntegrator1D::FunctionIntegrator1D():CObject(TYPE) {}

/** Explicit constructor */
FunctionIntegrator1D::FunctionIntegrator1D(Integrator1DSP integrator):
    CObject(TYPE), integrator1D(integrator) {}

IObject* FunctionIntegrator1D::defaultConstructor()
{
    return new FunctionIntegrator1D();
}

void FunctionIntegrator1D::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(FunctionIntegrator1D, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IIntegrator);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(integrator1D, "Type of 1D integrator to use");
}

/** Integration method */
IObjectSP FunctionIntegrator1D::integrate(const IIntegrand* integrand) const
{
    const MFunctionND* f = dynamic_cast<const MFunctionND*>(integrand);
    if (f == 0) {
        throw ModelException(
            "FunctionIntegrator1D::integrate",
            "Unable to cast integrand into MFunctionND");
    }
    if (f->getNbFuncs() != 1)
    {
        throw ModelException(
            "FunctionIntegrator1D::integrate",
            "Invalid integrand: nb functions = " +
            Format::toString(f->getNbFuncs()) +
            " (1 expected).");
    }
    if (f->getNbVars() != 1)
    {
        throw ModelException(
            "FunctionIntegrator1D::integrate",
            "Invalid integrand: nb variables = " +
            Format::toString(f->getNbVars()) +
            " (1 expected).");
    }

    Function1DDoubleConstSP function = DYNAMIC_POINTER_CAST<Function1DDouble>(
        FunctionOperations::project(*f, 0));
    
    return CDoubleSP(CDouble::create(integrator1D->integrate(*function)));
}

/** Addin for unit testing of 1D integrator engines */
class UTIL_DLL Integrator1DUnitTest:
    public CObject,
    public ClientRunnable {
public:
    
    // inner class for id function
    class IdentityFunction: public virtual Function1DDouble {
    public:
        virtual double operator()(double  x) const {
            return x;
        }
        
        IdentityFunction(const Range& interval): Function1DDouble(interval){}
	};

    // inner class for exp function
    class ExpFunction: public virtual Function1DDouble {
    public:
        virtual double operator()(double  x) const {
            return exp(x);
        }

        ExpFunction(const Range& interval): Function1DDouble(interval){}
    };

    static CClassConstSP const TYPE;

    // EdrAction version of addin
    IObjectSP run() {
        Range range(lowerBound, true, upperBound, true);
        if (integrand == ID)
        {
            IdentityFunction f(range);
            return integrator->integrate(&f);
        }
        else if (integrand == EXP)
        {
            ExpFunction f(range);
            return integrator->integrate(&f);
        }
        else
        {
            throw ModelException(
                "Integrator1DUnitTest::run",
                "Integrand " + integrand + " is not supported");
        }
    }

private:
    Integrator1DUnitTest(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Integrator1DUnitTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(integrator, "The 1D integrator to test");
        FIELD(integrand,
            "Function to integrate. "
            "Possible values are: " +
            ID + " (x->x), " +
            EXP + " (x->exp(x))");
        FIELD(lowerBound, "Lower integration bound");
        FIELD(upperBound, "Upper integration bound");
    }
    
    static IObject* defaultConstructor(){
        return new Integrator1DUnitTest();
    }

    // Addin parameters
    FunctionIntegrator1DSP integrator;
    string integrand;
    double lowerBound;
    double upperBound;

    // String for basic functions
    static const string ID;
    static const string EXP;
};

const string Integrator1DUnitTest::ID = "ID";
const string Integrator1DUnitTest::EXP = "EXP";

CClassConstSP const Integrator1DUnitTest::TYPE =
    CClass::registerClassLoadMethod(
        "Integrator1DUnitTest",
        typeid(Integrator1DUnitTest),
        load);

/* external symbol to allow class to be forced to be linked in */
bool FunctionIntegrator1DLoad(){
    return (Integrator1DUnitTest::TYPE != 0);
}

DRLIB_END_NAMESPACE
