//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 31-Jul-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_INTEGRATORND_HPP
#define QLIB_INTEGRATORND_HPP

#include "edginc/Function.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for N-dimensional integrators (f:R^N->R) */
class UTIL_DLL IntegratorND: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Integrates func: R^N->R*/
    virtual double integrate(const FunctionNDDouble& func) const = 0;
};

DECLARE(IntegratorND);

DRLIB_END_NAMESPACE

#endif /*QLIB_INTEGRATORND_HPP*/
