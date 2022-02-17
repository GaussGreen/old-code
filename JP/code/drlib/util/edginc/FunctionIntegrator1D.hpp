//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 28-Jul-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FUNCTIONINTEGRATOR1D_HPP
#define QLIB_FUNCTIONINTEGRATOR1D_HPP

#include "edginc/IIntegrator.hpp"
#include "edginc/Integrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** "Classic" integrator for function f: R --> R^M */
class UTIL_DLL FunctionIntegrator1D:
    public CObject,
    public virtual IIntegrator
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~FunctionIntegrator1D();
    
    /** Integration method */
    virtual IObjectSP integrate(const IIntegrand* integrand) const;

    /** Explicit constructor */
    FunctionIntegrator1D(Integrator1DSP integrator);

private:
    FunctionIntegrator1D();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    Integrator1DSP integrator1D; // $required(Type of 1D integrator to use)
};

DECLARE(FunctionIntegrator1D);

DRLIB_END_NAMESPACE

#endif /*FUNCTIONINTEGRATOR1D*/
