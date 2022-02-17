//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 28-Jul-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IINTEGRATOR_HPP
#define QLIB_IINTEGRATOR_HPP

#include "edginc/Object.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/IIntegrand.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for integration engines.
 * Types in "integrate" method are deliberately very generic to
 * cope with different possible implementations (eg: "classic" function 
 * integrator, "monte carlo" integrator...). 
 * */
class UTIL_DLL IIntegrator: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /**
     * Generic integration method:
     * - "integrand" could be any type of integrand (a function, a process...)
     * - the integration result could be a double, double array...
     * The range of integration is contained within the integrand.
     * */
    virtual IObjectSP integrate(const IIntegrand* integrand) const = 0; 
};

DECLARE(IIntegrator);

DRLIB_END_NAMESPACE

#endif /*QLIB_IINTEGRATOR_HPP*/
