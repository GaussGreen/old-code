//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IntFuncNDGeneral.hpp
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

#ifndef EDR_INTFUNCNDGENERAL_HPP
#define EDR_INTFUNCNDGENERAL_HPP

DRLIB_BEGIN_NAMESPACE

/** Wrapper around imsl's int_fcn  */
class UTIL_DLL IntFuncNDGeneral:public CObject, public IntegratorND {
public:
    static CClassConstSP const TYPE;
    friend class IntFuncNDGeneralHelper;

    virtual double integrate(const FunctionNDDouble& func) const;

    static const double default_err_abs;
    static const double default_err_rel;
    static const int    default_max_evals;

	IntFuncNDGeneral(   double err_abs = default_err_abs,
                        double err_rel = default_err_rel,
                        int max_evals = default_max_evals);

    static IntFuncNDGeneral* create(    double err_abs = default_err_abs,
                                        double err_rel = default_err_rel,
                                        int max_evals = default_max_evals);

private:
    IntFuncNDGeneral(int notUsed);

    double err_abs;
    double err_rel;
    int max_evals;
};

typedef smartPtr<IntFuncNDGeneral> IntFuncNDGeneralSP;
typedef smartConstPtr<IntFuncNDGeneral> IntFuncNDGeneralConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<IntFuncNDGeneral>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<IntFuncNDGeneral>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<IntFuncNDGeneral>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<IntFuncNDGeneral>);
#endif

DRLIB_END_NAMESPACE

#endif
