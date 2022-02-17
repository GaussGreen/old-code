#include "edginc/config.hpp"
#include "edginc/RectangleIntegrator1D.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(RectangleIntegrator1DArray);

/** TYPE (for reflection) */        
CClassConstSP const RectangleIntegrator1D::TYPE =
CClass::registerClassLoadMethod(
    "RectangleIntegrator1D",
    typeid(RectangleIntegrator1D),
    RectangleIntegrator1D::load);

/** Virtual destructor */
RectangleIntegrator1D::~RectangleIntegrator1D()
{
}

/** Constructor */
RectangleIntegrator1D::RectangleIntegrator1D():CObject(TYPE) {}

/** Explicit constructor */
RectangleIntegrator1D::RectangleIntegrator1D(int nbSubdivisions):
    CObject(TYPE), nbSubdivisions(nbSubdivisions) {}

IObject* RectangleIntegrator1D::defaultConstructor()
{
    return new RectangleIntegrator1D();
}

void RectangleIntegrator1D::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(RectangleIntegrator1D, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(Integrator1D);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(nbSubdivisions, "Number of subdivisions");
}

/**
 * Integration method:
 * 
 * We want to compute the integral int(a,b,f(x)dx).
 * 
 * We divide the integration interval [a,b] into n = nbSubdivisions
 * intervals [I_k, I_k+1], k=0..(n-1) where I_k = a+k(b-a)/n
 * 
 * For each interval, we estimate int(I_k,I_k+1,f(x)dx) by
 * E_k = f((I_k+I_k+1)/2)*(I_k+1 - I_k) i.e. f(a+(k+0.5)*(b-a)/n)*(b-a)/n
 * 
 * We estimate int(a,b,f(x)dx) by sum(k=0..(n-1), E_k)
 * 
 * */
double RectangleIntegrator1D::integrate(const Function1DDouble& func) const
{
    static const string method = "RectangleIntegrator1D::integrate";
    
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()){
            throw ModelException(method,
                                 "(Semi-)Open and / or infinite intervals are not supported; got " + interval.toString());
        }

        double a = lower.getValue();
        double b = upper.getValue();

        double dx = (b-a)/nbSubdivisions;
        double x = a + 0.5*dx;
        
        double integral = 0.0;
        for (int k = 0; k < nbSubdivisions; ++k) {
			integral += func(x) * dx;
            x += dx;
		}
        
        return integral;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/* external symbol to allow class to be forced to be linked in */
bool RectangleIntegrator1DLoad(){
    return (RectangleIntegrator1D::TYPE != 0);
}

DRLIB_END_NAMESPACE
