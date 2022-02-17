//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 31-Jul-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FUNCTIONINTEGRATORND_HPP
#define QLIB_FUNCTIONINTEGRATORND_HPP

#include "edginc/IIntegrator.hpp"
#include "edginc/IntegratorND.hpp"

DRLIB_BEGIN_NAMESPACE

/** "Classic" integrator for function f: R^N --> R^M */
class UTIL_DLL FunctionIntegratorND:
    public CObject,
    public virtual IIntegrator
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~FunctionIntegratorND();
    
    /** Integration method */
    virtual IObjectSP integrate(const IIntegrand* integrand) const;

private:
    FunctionIntegratorND();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    IntegratorNDSP integratorND; // $required(Type of ND integrator to use)
};

DECLARE(FunctionIntegratorND);

DRLIB_END_NAMESPACE

#endif /*FUNCTIONINTEGRATORND*/
