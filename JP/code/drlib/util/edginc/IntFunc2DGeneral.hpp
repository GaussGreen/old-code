//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IntFunc2DGeneral.hpp
//
//   Description :
//
//   Date        : 4 Sep 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Function.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/IntegratorND.hpp"
#include <math.h>

#ifndef EDR_INTFUNC2DGENERAL_HPP
#define EDR_INTFUNC2DGENERAL_HPP

DRLIB_BEGIN_NAMESPACE

/** Wrapper around imsl's int_fcn  */
class UTIL_DLL IntFunc2DGeneral:public CObject, public IntegratorND {
public:
    static CClassConstSP const TYPE;
    friend class IntFunc2DGeneralHelper;

    virtual double integrate(const FunctionNDDouble& func) const;

    static const double default_err_abs;
    static const double default_err_rel;
    static const int    default_max_evals;
    static const int    default_max_subinter;

	IntFunc2DGeneral(   double err_abs = default_err_abs,
                        double err_rel = default_err_rel,
                        int max_evals = default_max_evals,
                        int max_subinter = default_max_subinter);

    static IntFunc2DGeneral* create(    double err_abs = default_err_abs,
                                        double err_rel = default_err_rel,
                                        int max_evals = default_max_evals,
                                        int max_subinter = default_max_subinter);

private:
    IntFunc2DGeneral(int notUsed);

    double err_abs;
    double err_rel;
    int max_evals;
    int max_subinter;

    // internal
    double (*fc)(double);
    double (*fd)(double);
};

typedef smartPtr<IntFunc2DGeneral> IntFunc2DGeneralSP;
typedef smartConstPtr<IntFunc2DGeneral> IntFunc2DGeneralConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<IntFunc2DGeneral>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<IntFunc2DGeneral>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<IntFunc2DGeneral>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<IntFunc2DGeneral>);
#endif

DRLIB_END_NAMESPACE

#endif
