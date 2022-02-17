//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 28-Jul-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IINTEGRAND_HPP
#define QLIB_IINTEGRAND_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/VirtualDestructorBase.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_REF_COUNT(IIntegrand);

/** Marker interface for integrands (used in IIntegrator). */
class UTIL_DLL IIntegrand: public virtual VirtualDestructorBase {
public:

    /**
     * Split this integrand into a vector of "1D" integrands
     * (eg: if integrand is a function f:R^N->R^M, to1DIntegrands() will
     * return M integrands that will be functions R^N->R).
     * */
    virtual IIntegrandArrayConstSP to1DIntegrands() const = 0;

};

DRLIB_END_NAMESPACE

#endif /*QLIB_IINTEGRAND_HPP*/
