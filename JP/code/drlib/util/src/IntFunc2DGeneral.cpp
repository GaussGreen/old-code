//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IntFunc2DGeneral.cpp
//
//   Description : 
//
//   Date        : 4 Sep 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_INTEGRATOR_CPP
#include "edginc/IntFunc2DGeneral.hpp"
#include "edginc/IMSLFunction2DProvider.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

// IntFuncNDgeneral
const double IntFunc2DGeneral::default_err_abs = sqrt(DBL_EPSILON);
const double IntFunc2DGeneral::default_err_rel = sqrt(DBL_EPSILON);
const int    IntFunc2DGeneral::default_max_evals = 500;
const int    IntFunc2DGeneral::default_max_subinter = 500;

/** imsl's int_FuncND_inf */
double IntFunc2DGeneral::integrate(const FunctionNDDouble& infunc) const {
    static const string method = "IntFunc2DGeneral::integrate";
    
    try{
        int n = infunc.getNbVars();
        if (n != 2)
        {
            throw ModelException(method,
                "internal error: function provided is not a 2d function");            
        }

        const RangeArray& interval = infunc.getIntervals();

        Range::checkIsNonEmpty(*(interval[0]));
        if ((*(interval[0])).isSingleton()){
            return 0.0;
        }
        const Boundary& lower = (*(interval[0])).getLower();
        const Boundary& upper = (*(interval[0])).getUpper();

        if ( !lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()) {
            throw ModelException(method,
                                "(Semi-)Open and / or infinite intervals are not supported; got " + (*(interval[0])).toString());
        }
        

        double a, b;
        a = lower.getValue();
        b = upper.getValue();
        
        IMSLFunction2DProvider func2DProvider(&infunc);
        double err_est;
        int n_subinter, n_evals;
        double rtn = imsl_d_int_fcn_2d(func2DProvider.function(),
                                    a, b, func2DProvider.lowerBound(), func2DProvider.upperBound(),
                                    IMSL_ERR_ABS, err_abs,
                                    IMSL_ERR_REL, err_rel,
                                    IMSL_ERR_EST, &err_est,
                                    IMSL_MAX_SUBINTER, max_subinter,
                                    IMSL_N_SUBINTER, &n_subinter,
                                    IMSL_N_EVALS, &n_evals,
                                    0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

IntFunc2DGeneral* IntFunc2DGeneral::create( double err_abs,
                                            double err_rel,
                                            int max_evals,
                                            int max_subinter){
    return new IntFunc2DGeneral(err_abs, err_rel, max_evals);
}

IntFunc2DGeneral::IntFunc2DGeneral( double err_abs,
                                    double err_rel,
                                    int max_evals,
                                    int max_subinter):
CObject(TYPE),
err_abs(err_abs),
err_rel(err_rel),
max_evals(max_evals),
max_subinter(max_subinter){}

IntFunc2DGeneral::IntFunc2DGeneral(int notUsed):
CObject(TYPE),
err_abs(default_err_abs),
err_rel(default_err_rel),
max_evals(default_max_evals),
max_subinter(default_max_subinter){}

class IntFunc2DGeneralHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IntFunc2DGeneral, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IntegratorND);
        EMPTY_SHELL_METHOD(defaultIntFunc2DGeneral);
        FIELD(err_abs, "Absolute Error");
        FIELD_MAKE_OPTIONAL(err_abs);
        FIELD(err_rel, "Relative Error");
        FIELD_MAKE_OPTIONAL(err_rel);
        FIELD(max_evals, "Maximum Nb of function evaluation");
        FIELD_MAKE_OPTIONAL(max_evals);
        FIELD(max_subinter, "Maximum Nb of Sub-intervals");
        FIELD_MAKE_OPTIONAL(max_subinter);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultIntFunc2DGeneral(){
        return new IntFunc2DGeneral(0);
    }
};

CClassConstSP const IntFunc2DGeneral::TYPE = CClass::registerClassLoadMethod(
    "IntFunc2DGeneral", typeid(IntFunc2DGeneral), IntFunc2DGeneralHelper::load);

// force linker to include this file (avoid having header file) */
bool IntFunc2DGeneralLoad() {
    return (IntFunc2DGeneral::TYPE != 0);
}

DRLIB_END_NAMESPACE