//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 15-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RECTANGLEINTEGRATOR1D_HPP
#define QLIB_RECTANGLEINTEGRATOR1D_HPP

#include "edginc/Integrator.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Basic "rectangle" integrator:

 * We want to compute the integral int(a,b,f(x)dx).
 * 
 * We divide the integration interval [a,b] into n = nbSubdivisions (user input)
 * intervals [I_k, I_k+1], k=0..(n-1) where I_k = a+k(b-a)/n
 * 
 * For each interval, we estimate int(I_k,I_k+1,f(x)dx) by
 * E_k = f((I_k+I_k+1)/2)*(I_k+1 - I_k) i.e. f(a+(k+0.5)*(b-a)/n)*(b-a)/n
 * 
 * We estimate int(a,b,f(x)dx) by sum(k=0..(n-1), E_k)
 * */
class UTIL_DLL RectangleIntegrator1D:
    public CObject,
    public virtual Integrator1D
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~RectangleIntegrator1D();

    /** Integration method */
    virtual double integrate(const Function1DDouble& func) const;
    
    /** Explicit constructor */
    RectangleIntegrator1D(int nbSubdivisions);

private:
    RectangleIntegrator1D();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    int nbSubdivisions; // $required(Number of subdivisions)
};

DECLARE(RectangleIntegrator1D);

DRLIB_END_NAMESPACE

#endif /*RECTANGLEINTEGRATOR1D*/
