//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IntFuncNDGeneral.cpp
//
//   Description : 
//
//   Date        : 4 Sep 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_INTEGRATOR_CPP
#include "edginc/IntFuncNDGeneral.hpp"
#include "edginc/IMSLFunctionNDProvider.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

// IntFuncNDgeneral
const double IntFuncNDGeneral::default_err_abs = sqrt(DBL_EPSILON);
const double IntFuncNDGeneral::default_err_rel = sqrt(DBL_EPSILON);
const int    IntFuncNDGeneral::default_max_evals = 500;

/** imsl's int_FuncND_inf */
double IntFuncNDGeneral::integrate(const FunctionNDDouble& infunc) const {
    static const string method = "IntFuncNDGeneral::integrate";
    
    try{
        int n = infunc.getNbVars();
        const RangeArray& interval = infunc.getIntervals();

        BoundaryArray lower(n);
        BoundaryArray upper(n);

        for (int i=0; i<n; ++i)
        {
            Range::checkIsNonEmpty(*(interval[i]));
            if ((*(interval[i])).isSingleton()){
                return 0.0;
            }
            *(lower[i]) = (*(interval[i])).getLower();
            *(upper[i]) = (*(interval[i])).getUpper();

            if (!(*(lower[i])).isClosedBracket() || (*(lower[i])).isInfinite()
                || !(*(upper[i])).isClosedBracket() || (*(upper[i])).isInfinite()) {
                throw ModelException(method,
                                    "(Semi-)Open and / or infinite intervals are not supported; got " + (*(interval[i])).toString());
            }
        }

        CDoubleArray a = CDoubleArray(n);
        CDoubleArray b = CDoubleArray(n);
        for (int i=0; i<n; ++i)
        {
            a[i] = (*(lower[i])).getValue();
            b[i] = (*(upper[i])).getValue();
        }

        IMSLFunctionNDProvider funcNDProvider(&infunc);
        double err_est;

        double rtn = imsl_d_int_fcn_hyper_rect(funcNDProvider.function(), 
                                    infunc.getNbVars(), &a[0], &b[0],                    
                                    IMSL_ERR_ABS, err_abs,
                                    IMSL_ERR_REL, err_rel,
                                    IMSL_ERR_EST, &err_est,
                                    IMSL_MAX_EVALS, max_evals,
                                    0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

IntFuncNDGeneral* IntFuncNDGeneral::create( double err_abs,
                                            double err_rel,
                                            int max_evals){
    return new IntFuncNDGeneral(err_abs, err_rel, max_evals);
}

IntFuncNDGeneral::IntFuncNDGeneral( double err_abs,
                                    double err_rel,
                                    int max_evals):
CObject(TYPE),
err_abs(err_abs),
err_rel(err_rel),
max_evals(max_evals){}

IntFuncNDGeneral::IntFuncNDGeneral(int notUsed):
CObject(TYPE),
err_abs(default_err_abs),
err_rel(default_err_rel),
max_evals(default_max_evals){}

class IntFuncNDGeneralHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IntFuncNDGeneral, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IntegratorND);
        EMPTY_SHELL_METHOD(defaultIntFuncNDGeneral);
        FIELD(err_abs, "Absolute Error");
        FIELD_MAKE_OPTIONAL(err_abs);
        FIELD(err_rel, "Relative Error");
        FIELD_MAKE_OPTIONAL(err_rel);
        FIELD(max_evals, "Maximum Nb of function evaluation");
        FIELD_MAKE_OPTIONAL(max_evals);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultIntFuncNDGeneral(){
        return new IntFuncNDGeneral(0);
    }
};

CClassConstSP const IntFuncNDGeneral::TYPE = CClass::registerClassLoadMethod(
    "IntFuncNDGeneral", typeid(IntFuncNDGeneral), IntFuncNDGeneralHelper::load);

// force linker to include this file (avoid having header file) */
bool IntFuncNDGeneralLoad() {
    return (IntFuncNDGeneral::TYPE != 0);
}



DRLIB_END_NAMESPACE