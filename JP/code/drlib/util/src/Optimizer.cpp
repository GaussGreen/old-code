//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Optimizer.cpp
//
//   Description : 
//
//   Date        : 20 May 2002
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_OPTIMIZER_CPP
#include "edginc/Optimizer.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/imsl.h"
#include "edginc/Maths.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Nrutil.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CommandLineParams.hpp"

DRLIB_BEGIN_NAMESPACE

// OPTIMIZER ND
CClassConstSP const OptimizerND::TYPE = CClass::registerClassLoadMethod(
    "OptimizerND", typeid(OptimizerND), OptimizerND::load);

OptimizerND::OptimizerND(const CClassConstSP& clazz):
CObject(clazz){}

void OptimizerND::load(CClassSP& clazz){
    REGISTER(OptimizerND, clazz);
    SUPERCLASS(CObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

void OptimizerND::minimize(const MFunctionND&  func,
                           const CDoubleArray& xguess,        // initial guess
                           CDoubleArray&       x,
                           CDoubleArray&       f) const{
    minimize(func,
             xguess,
             x);
    func(x, f);
}

void OptimizerND::minimize(const MFunctionND&  func,
                           const CDoubleArray& xguess,        // initial guess
                           const CStringArray& ids,           // identifiers
                           CDoubleArray&       x,
                           CDoubleArray&       f) const{
    minimize(func,
             xguess,
             ids,
             x);
    func(x, f);
}

CClassConstSP const LeastSquareOptimizer::TYPE = 
CClass::registerClassLoadMethod(
    "LeastSquareOptimizer", typeid(LeastSquareOptimizer), LeastSquareOptimizer::load);

LeastSquareOptimizer::LeastSquareOptimizer(const CClassConstSP& clazz):
OptimizerND(clazz){}

void LeastSquareOptimizer::load(CClassSP& clazz){
    REGISTER(LeastSquareOptimizer, clazz);
    SUPERCLASS(OptimizerND);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

static const MFunctionND* local_func = 0;
static OneToOneMappingConstArray* local_mapping = 0;
static CDoubleArray* local_x = 0;   // no reallocation of memory
static CDoubleArray* local_f = 0;   // no reallocation of memory
static const CStringArray* local_ids = 0;   // no reallocation of memory

static void LM_fcn (int m, int n, double inx[], double outf[]) {
    // if function failed earlier, skip
    if (IMSLError::isError()){
        int j = 0;
        for (; j < m; ++j){
            outf[j] = 0.0;
        }
        return;
    }
    DoubleArray& x = *local_x; // for ease
    DoubleArray& f = *local_f; // for ease
    const StringArray& ids = *local_ids;

    // map vars
    int i = 0;
    for (; i < n; ++i){
        const OneToOneMapping& thisMapping = *(*local_mapping)[i];
        x[i] = thisMapping(inx[i]);  // mapp from R to (a, b)
    }
    try{
        // call objective func
        (*local_func)(x, f);
    }
    catch(exception& e){
        if (ids.size() > 0){
            string message = "LM_fcn: Failed with values";
            int i = 0;
            for (; i < ids.size(); ++i){
                message += "\n\t" + ids[i] + ": " + Format::toString(x[i]);
            }
            IMSLError::appendToStack(ModelException::addTextToException(e, message));
        }
        else{
            IMSLError::appendToStack(e);
        }
        int j = 0;
        for (; j < m; ++j){
            outf[j] = 0.0;
        }
        return;
    }
    // copy func vals
    int j = 0;
    for (; j < m; ++j){
        outf[j] = f[j];
    }
}

void LevenbergMarquardt::minimize(const MFunctionND&  func,
                                  const CDoubleArray& inxguess,      // initial guess
                                  CDoubleArray&       outx) const {  // result
    CStringArray ids(0);
    minimize(func,
             inxguess,
             ids,
             outx);
}

void LevenbergMarquardt::minimize(const MFunctionND&  func,
                                  const CDoubleArray& inxguess,      // initial guess
                                  const CStringArray& ids,           // identifiers
                                  CDoubleArray&       outx) const {  // result
    const string method = "LevenbergMarquardt::minimize";

    try{
        int m = func.getNbFuncs();
        int n = func.getNbVars();

        if (m <= 0){
            throw ModelException(method, 
                                 Format::toString(m) + " number of functions to minimize");
        }

        if (n <= 0){
            throw ModelException(method, 
                                 Format::toString(n) + " number of variables to minimize");
        }

        if (inxguess.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of xguess should be %ld; got %ld",
                                                  n, inxguess.size()));
        }

        if (outx.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of outx should be %ld; got %ld",
                                                  n, outx.size()));
        }

        /* Create mapping functions using ranges */
        const RangeArray& range = func.getIntervals();
        OneToOneMappingConstArray mapping(n);
        int i = 0;
        for(; i < n; ++i){
            mapping[i] = OneToOneMapping::create(*(range[i]));
        }

        /* Map var guesses from (a, b) to R */
        DoubleArray xguess(n);
        for (i = 0; i < n; ++i){
            xguess[i] = mapping[i]->inverse(inxguess[i]);
        }

        /* Create local copies and call imsl */
        local_func = &func;
        DoubleArray x(n); local_x = &x;
        DoubleArray f(m); local_f = &f;
        local_mapping = &mapping;
        local_ids = &ids;

        DoubleArray e(m);   // residuals
        double* imsl_jtj_inv = 0;
        bool ignoreStats;
        int dfe;    // degrees of freedom
        // NB: if m <= n, imsl cores when computing the rank of the system
        // and the inverse of the JTJ, so don't request them in that case
        if (m > n){
            ignoreStats = false;
            int rank;     // rank of Jacobian matrix
            if (!gradTol) // validatePop2Object is such that 
                          // (one optional param is NULL <=> all optional params are NULL)
            {
                imsl_d_nonlin_least_squares (LM_fcn, m, n,
                                         IMSL_XGUESS, &xguess[0],
                                         IMSL_RETURN_USER, &outx[0],
                                         IMSL_FVEC_USER, &e[0],
                                         IMSL_JTJ_INVERSE, &imsl_jtj_inv,
                                         IMSL_RANK, &rank,
                                         0);
            }
            else
            {
                imsl_d_nonlin_least_squares (LM_fcn, m, n,
                                         IMSL_XGUESS, &xguess[0],
                                         IMSL_RETURN_USER, &outx[0],
                                         IMSL_FVEC_USER, &e[0],
                                         IMSL_JTJ_INVERSE, &imsl_jtj_inv,
                                         IMSL_RANK, &rank,
                                         IMSL_GRAD_TOL, gradTol->doubleValue(),
                                         IMSL_STEP_TOL, stepTol->doubleValue(),
                                         IMSL_REL_FCN_TOL, relFcnTol->doubleValue(),
                                         IMSL_ABS_FCN_TOL, absFcnTol->doubleValue(),
                                         IMSL_MAX_STEP, maxStepSize->doubleValue(),
                                         IMSL_MAX_ITN, maxItnNb->intValue(),
                                         IMSL_MAX_FCN, maxFcnEvalNb->intValue(),
                                         IMSL_MAX_JACOBIAN, maxJacobEvalNb->intValue(),
                                         0);
            }
            dfe = m - rank;
        }
        else{
            ignoreStats = true;
            dfe = 0;
            if (!gradTol) // validatePop2Object is such that 
                          // (one optional param is NULL <=> all optional params are NULL)
            {
                imsl_d_nonlin_least_squares (LM_fcn, m, n,
                                         IMSL_XGUESS, &xguess[0],
                                         IMSL_RETURN_USER, &outx[0],
                                         IMSL_FVEC_USER, &e[0],
                                         0);
            }
            else
            {
                 imsl_d_nonlin_least_squares (LM_fcn, m, n,
                                         IMSL_XGUESS, &xguess[0],
                                         IMSL_RETURN_USER, &outx[0],
                                         IMSL_FVEC_USER, &e[0],
                                         IMSL_GRAD_TOL, gradTol->doubleValue(),
                                         IMSL_STEP_TOL, stepTol->doubleValue(),
                                         IMSL_REL_FCN_TOL, relFcnTol->doubleValue(),
                                         IMSL_ABS_FCN_TOL, absFcnTol->doubleValue(),
                                         IMSL_MAX_STEP, maxStepSize->doubleValue(),
                                         IMSL_MAX_ITN, maxItnNb->intValue(),
                                         IMSL_MAX_FCN, maxFcnEvalNb->intValue(),
                                         IMSL_MAX_JACOBIAN, maxJacobEvalNb->intValue(),
                                         0);
           }
        }

        local_ids = 0;
        local_mapping = 0;
        local_f = 0;
        local_x = 0;
        local_func = 0;

        try{
            IMSLError::throwExceptionIfError();
        }
        catch(exception& e){
            // it may be that imsl_d_nonlin_least_squares failed simply because it
            // failed to invert the Jacobian. In that case, don't report stats
            // but do report calibrated values
            if (IMSLError::getErrorCode() == IMSL_REMAINING_ELMNTS_NOT_ZERO){
                ignoreStats = true;
            }
            else{
                free(imsl_jtj_inv);
                if (ids.size() > 0){
                    string message = "imsl_d_nonlin_least_squares : Failed with values";
                    int i = 0;
                    for (; i < ids.size(); ++i){
                        const OneToOneMapping& thisMapping = (*mapping[i]);
                        message += "\n\t" + ids[i] + ": " + Format::toString(thisMapping(outx[i]));
                    }
                    throw ModelException::addTextToException(e, message);
                }
                throw;
            }
        }

        
        // var
        if (haveStats = (!ignoreStats && (dfe > 0))){
            CDoubleArray mapderiv(n);
            for (i = 0; i < n; ++i){
                mapderiv[i] = mapping[i]->derivative(outx[i]);
            }

            // inverse of hessian
            jtj_inv = CDoubleMatrixSP(new DoubleMatrix(n, n));
            for (i = 0; i < n; ++i){
                int j = 0;
                for (; j < n; ++j){
                    (*jtj_inv)[i][j] = mapderiv[i] * imsl_jtj_inv[i * n + j] * mapderiv[j] ;
                }
            }

            sqSigma = 0.0;
            int j = 0;
            for (; j < m; j++){
                sqSigma += e[j] * e[j];
            }
            sqSigma /= dfe;

            // covariance matrix
            covar = CDoubleMatrixSP(new DoubleMatrix(n, n));
            for (i = 0; i < n; ++i){
                int j = 0;
                for (; j < n; ++j){
                    (*covar)[i][j] = sqSigma * (*jtj_inv)[i][j];
                }
            }

            // Lower, upper limits
            lowerLims.resize(n);
            upperLims.resize(n);
            for (i = 0; i < n; ++i){
                double a = imsl_d_t_inverse_cdf (0.975, static_cast<double>(dfe)) 
                            * sqrt(sqSigma * imsl_jtj_inv[i * n + i]);
                const OneToOneMapping& thisMapping = (*mapping[i]);
                lowerLims[i] = thisMapping(outx[i] - a);
                upperLims[i] = thisMapping(outx[i] + a);
            }
 
            // check whether any of the  standard errors is equal to zero
            // if yes then the statistic are not provided
            stdError.resize(n);
            for (i = 0; i < n; i++){
                stdError[i] = sqrt((*covar)[i][i]);
                haveStats = haveStats&&(!Maths::isZero(stdError[i]));
            }
                    
            // correlation matrix       
            correl = CDoubleMatrixSP(new DoubleMatrix(n, n));
            
            for (i = 0; i < n; i++){
                (*correl)[i][i] = 1.0;
                for (j = i+1; j < n; j++){
                    (*correl)[i][j] = (*covar)[i][j]/(stdError[i] * stdError[j]);
                    (*correl)[j][i] = (*correl)[i][j];
                }
            }

        }

        // Map vars from R to (a, b)
        for (i = 0; i < n; ++i){
            const OneToOneMapping& thisMapping = (*mapping[i]);
            outx[i] = thisMapping(outx[i]);
        }

        free(imsl_jtj_inv);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

bool LevenbergMarquardt::statsIsAvailable() const{
    return haveStats;
}

double LevenbergMarquardt::getVariance() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getVariance",
                             "Stats is not available");
    }
    return sqSigma;
}

const DoubleMatrix& LevenbergMarquardt::getCovarianceMatrix() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getCovarianceMatrix",
                             "Stats is not available");
    }
    return *covar;
}

const DoubleMatrix& LevenbergMarquardt::getCorrelationMatrix() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getCorrelationMatrix",
                             "Stats is not available");
    }
    return *correl;
}

const DoubleMatrix& LevenbergMarquardt::getHessianInverse() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getHessianInverse",
                             "Stats is not available");
    }
    return *jtj_inv;
}

const DoubleArray& LevenbergMarquardt::getLowerLimits() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getLowerLimits",
                             "Stats is not available");
    }
    return lowerLims;
}

const DoubleArray& LevenbergMarquardt::getUpperLimits() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getUpperLimits",
                             "Stats is not available");
    }
    return upperLims;
}

const DoubleArray& LevenbergMarquardt::getStandardErrors() const{
    if(!haveStats){
        throw ModelException("LevenbergMarquardt::getStandardError",
                             "Stats is not available");
    }
    return stdError;
}

const double LevenbergMarquardt::DefaultMachinePrecision = 1.0e-16; // just some sensible value, 
                                                                    // might be not far away from the right one
// the numbers below would replicate the default values of IMSP if the machine precision were correct -------
const double LevenbergMarquardt::DEFAULT_gradTol 
    = pow(LevenbergMarquardt::DefaultMachinePrecision, 1.0e0/3.0e0);
const double LevenbergMarquardt::DEFAULT_stepTol 
    = pow(LevenbergMarquardt::DefaultMachinePrecision, 2.0e0/3.0e0);
const double LevenbergMarquardt::DEFAULT_relFcnTol 
    = imsl_f_max(1.0e-20, pow(LevenbergMarquardt::DefaultMachinePrecision, 2.0e0/3.0e0));
const double LevenbergMarquardt::DEFAULT_absFcnTol 
    = imsl_f_max(1.0e-10, pow(LevenbergMarquardt::DefaultMachinePrecision, 2.0e0/3.0e0));
const double LevenbergMarquardt::DEFAULT_maxStepSize = -999.0e0;
const int    LevenbergMarquardt::DEFAULT_maxItnNb = 100;
const int    LevenbergMarquardt::DEFAULT_maxFcnEvalNb = 400;
const int    LevenbergMarquardt::DEFAULT_maxJacobEvalNb = 400;


LevenbergMarquardt::LevenbergMarquardt()
    : LeastSquareOptimizer(TYPE)
    , haveStats(false)
    , sqSigma(0.0)
{}

void LevenbergMarquardt::validatePop2Object()
{
    if (!((!gradTol) &&  (!stepTol) && (!relFcnTol) &&  (!absFcnTol) 
            && (!maxStepSize) &&  (!maxItnNb) && (!maxFcnEvalNb) &&  (!maxJacobEvalNb)))
    // if no optional parameter has been inputed, do nothing.
    // alternatively, assign default values to the parameters that have not been inputed
    {
        if (!gradTol) { gradTol               = CDouble::SP(DEFAULT_gradTol); }
        if (!stepTol) { stepTol               = CDouble::SP(DEFAULT_stepTol); }
        if (!relFcnTol) { relFcnTol           = CDouble::SP(DEFAULT_relFcnTol); }
        if (!absFcnTol) { absFcnTol           = CDouble::SP(DEFAULT_absFcnTol); }
        if (!maxStepSize) { maxStepSize       = CDouble::SP(DEFAULT_maxStepSize); }
        if (!maxItnNb) { maxItnNb             = CInt::SP(DEFAULT_maxItnNb); }
        if (!maxFcnEvalNb) { maxFcnEvalNb     = CInt::SP(DEFAULT_maxFcnEvalNb); }
        if (!maxJacobEvalNb) { maxJacobEvalNb = CInt::SP(DEFAULT_maxJacobEvalNb); }
    }
}

class LevenbergMarquardtHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LevenbergMarquardt, clazz);
        SUPERCLASS(LeastSquareOptimizer);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(gradTol, "Scaled gradient tolerance. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_gradTol));
        FIELD_MAKE_OPTIONAL(gradTol);
        FIELD(stepTol, "Scaled step tolerance. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_stepTol));
        FIELD_MAKE_OPTIONAL(stepTol);
        FIELD(relFcnTol, "Relative function tolerance. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_relFcnTol));
        FIELD_MAKE_OPTIONAL(relFcnTol);
        FIELD(absFcnTol, "Absolute function tolerance. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_absFcnTol));
        FIELD_MAKE_OPTIONAL(absFcnTol);
        FIELD(maxStepSize, string("Maximum allowable step size. Default = ")
            + Format::toString(LevenbergMarquardt::DEFAULT_maxStepSize) + " . If maxStepSize is negative,"
            + "its value be computed as a function of other inputs. If it is = " 
            + Format::toString(LevenbergMarquardt::DEFAULT_maxStepSize) +
            + " the message waring will be suppressed.");
        FIELD_MAKE_OPTIONAL(maxStepSize);
        FIELD(maxItnNb, "Maximum number of iterations. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_maxItnNb));
        FIELD_MAKE_OPTIONAL(maxItnNb);
        FIELD(maxFcnEvalNb, "Maximum number of function evaluations. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_maxFcnEvalNb));
        FIELD_MAKE_OPTIONAL(maxFcnEvalNb);
        FIELD(maxJacobEvalNb, "Maximum number of Jacobian evaluations. Default = " + Format::toString(LevenbergMarquardt::DEFAULT_maxJacobEvalNb));
        FIELD_MAKE_OPTIONAL(maxJacobEvalNb);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new LevenbergMarquardt();
    }
};

CClassConstSP const LevenbergMarquardt::TYPE = CClass::registerClassLoadMethod(
    "LevenbergMarquardt", typeid(LevenbergMarquardt), LevenbergMarquardtHelper::load);


static double QN_fcn (int n, double x[]) {
    double rtn = 0.0;
    LM_fcn(1, n, x, &rtn);
    return rtn;
}

void QuasiNewton::minimize(const MFunctionND&  func,
                           const CDoubleArray& inxguess,      // initial guess
                           CDoubleArray&       outx) const {  // result
    CStringArray ids(0);
    minimize(func,
             inxguess,
             ids,
             outx);
}

void QuasiNewton::minimize(const MFunctionND&  func,
                           const CDoubleArray& inxguess,      // initial guess
                           const CStringArray& ids,           // identifiers
                           CDoubleArray&       outx) const {  // result
    const string method = "QuasiNewton::minimize";

    try{
        if (func.getNbFuncs() != 1){
            throw ModelException(method, 
                                 Format::toString("Nb of funcs should be 1; got %ld",
                                                  func.getNbFuncs()));
        }

        int n = func.getNbVars();

        if (n <= 0){
            throw ModelException(method, 
                                 Format::toString(n) + " number of variables to minimize");
        }

        if (inxguess.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of xguess should be %ld; got %ld",
                                                  n, inxguess.size()));
        }

        if (outx.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of outx should be %ld; got %ld",
                                                  n, outx.size()));
        }

        /* Create mapping functions using ranges */
        const RangeArray& range = func.getIntervals();
        OneToOneMappingConstArray mapping(n);
        int i = 0;
        for(; i < n; ++i){
            mapping[i] = OneToOneMapping::create(*(range[i]));
        }

        /* Map var guesses from (a, b) to R */
        DoubleArray xguess(n);
        for (i = 0; i < n; ++i){
            xguess[i] = mapping[i]->inverse(inxguess[i]);
        }

        /* Create local copies and call imsl */
        local_func = &func;
        DoubleArray x(n); local_x = &x;
        DoubleArray f(1); local_f = &f;
        local_mapping = &mapping;
        local_ids = &ids;

        imsl_d_min_uncon_multivar (QN_fcn, n,
                                   IMSL_XGUESS, &xguess[0],
                                   IMSL_RETURN_USER, &outx[0],
                                   IMSL_INIT_HESSIAN, 1,
                                   0);




        local_ids = 0;
        local_mapping = 0;
        local_f = 0;
        local_x = 0;
        local_func = 0;

        try{
            IMSLError::throwExceptionIfError();
        }
        catch(exception& e){
            if (ids.size() > 0){
                string message = "imsl_d_nonlin_least_squares : Failed with values";
                int i = 0;
                for (; i < ids.size(); ++i){
                    const OneToOneMapping& thisMapping = (*mapping[i]);
                    message += "\n\t" + ids[i] + ": " + Format::toString(thisMapping(outx[i]));
                }
                throw ModelException::addTextToException(e, message);
            }
            throw;
        }

        /* Map vars from R to (a, b) */
        for (i = 0; i < n; ++i){
            const OneToOneMapping& thisMapping = (*mapping[i]);
            outx[i] = thisMapping(outx[i]);
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

QuasiNewton::QuasiNewton():
OptimizerND(TYPE){}

class QuasiNewtonHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(QuasiNewton, clazz);
        SUPERCLASS(OptimizerND);
        EMPTY_SHELL_METHOD(defaultCtor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new QuasiNewton();
    }
};

CClassConstSP const QuasiNewton::TYPE = CClass::registerClassLoadMethod(
    "QuasiNewton", typeid(QuasiNewton), QuasiNewtonHelper::load);

// SIMPLEX
CClassConstSP const Simplex::IBasisMaker::TYPE = 
CClass::registerInterfaceLoadMethod("Simplex::IBasisMaker", typeid(Simplex::IBasisMaker), 0);

CClassConstSP const Simplex::CanonicalBasisMaker::TYPE = CClass::registerClassLoadMethod(
    "Simplex::CanonicalBasisMaker", typeid(Simplex::CanonicalBasisMaker), load);

const string Simplex::RELATIVE  = "relative"; 
const string Simplex::ABSOLUTE  = "absolute"; 
const string Simplex::COMPOSITE = "composite"; 
const string Simplex::USE_AMOEBA2 = "use_amoeba2";

Simplex::CanonicalBasisMaker::CanonicalBasisMaker():
CObject(TYPE){}

/** Invoked when Class is 'loaded' */
void Simplex::CanonicalBasisMaker::load(CClassSP& clazz){
    REGISTER(Simplex::CanonicalBasisMaker, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(Simplex::IBasisMaker);
    EMPTY_SHELL_METHOD(defaultCtor);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

IObject* Simplex::CanonicalBasisMaker::defaultCtor(){
    return new CanonicalBasisMaker();
}

DoubleMatrix Simplex::CanonicalBasisMaker::make(int n) const{
    DoubleMatrix basisMaker(n, n);
    for (int j = 0; j < n; ++j){        
        for (int i = 0; i < n; ++i){
            basisMaker[j][i] = (i == j ? 1 : 0);
        }
    }
    return basisMaker;
}


void Simplex::minimize(const MFunctionND&  func,
                       const CDoubleArray& inxguess,      // initial guess
                       CDoubleArray&       outx) const {  // result
    CStringArray ids(0);
    minimize(func,
             inxguess,
             ids,
             outx);
}


void Simplex::minimize(const MFunctionND&  func,
                       const CDoubleArray& inxguess,      // initial guess
                       const CStringArray& ids,           // identifiers
                       CDoubleArray&       outx) const {  // result
    const string method = "Simplex::minimize";
    try{
        if (func.getNbFuncs() != 1){
            throw ModelException(method, 
                                 Format::toString("Nb of funcs should be 1; got %ld",
                                                  func.getNbFuncs()));
        }
        int n = func.getNbVars();
        if (n <= 0){
            throw ModelException(method, 
                                 Format::toString(n) + " number of variables to minimize");
        }
        me = const_cast<Simplex*>(this);
        me->ndim = n;
        if (inxguess.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of xguess should be %ld; got %ld",
                                                  n, inxguess.size()));
        }
        if (outx.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of outx should be %ld; got %ld",
                                                  n, outx.size()));
        }
        /* Create mapping functions using ranges */
        const RangeArray& range = func.getIntervals();
        me->mappings.resize(n);
        int i = 0;
        for(; i < n; ++i){
            me->mappings[i] = OneToOneMapping::create(*(range[i]));
        }
        /* Map var guesses from (a, b) to R */
        DoubleArray xguess(n);
        for (i = 0; i < n; ++i){
            xguess[i] = mappings[i]->inverse(inxguess[i]);
        }
        /* Create local copies */
        me->func = &func;
        me->ids = &ids;
        me->x.resize(n);
        me->f.resize(1);
        /* create the vertices of the simplex */
        DoubleArrayArray vertices(n + 1);
        // initialize the first vertex to the initial guess
        vertices[0].resize(n);
        for (i = 0; i < n; ++i){
            vertices[0][i] = xguess[i];
        }
        // create 'n' additional vertices along the unit vectors
        DoubleMatrix basis(basisMaker->make(n));
        int j;
        for (j = 1; j < n + 1; ++j){
            vertices[j].resize(n);
            for (i = 0; i < n; ++i){
                if(CString::equalsIgnoreCase(useNewMapping,Simplex::RELATIVE)){
                    if(i == n-1){
                        vertices[j][i] = inxguess[i] * (1. +  0.1 * lengthScale * basis[j-1][i]);
                    }
                    else{
                        vertices[j][i] = inxguess[i] * (1. +  lengthScale * basis[j-1][i]);
                    }
                    vertices[j][i] = mappings[i]->inverse(vertices[j][i]);
                }else{
                    vertices[j][i] = xguess[i] + lengthScale * basis[j-1][i];
                }
                
            }
        }
        // evaluate obj func at each vertex
        vector<double*> p(n + 1);
        DoubleArray y(n + 1);
        for (j = 0; j < n + 1; ++j){
            p[j] = &vertices[j][0] - 1;
            y[j] = funk(p[j]);
        }

        
        // call nr routine
        if(CString::equalsIgnoreCase(isSRM,Simplex::USE_AMOEBA2)){
            amoeba2(&p[0] - 1 , &y[0] - 1 , ndim, ftol,ftolAbs,
                    funk, &me->nfunk, 1);
        }else{
            amoeba(&p[0] - 1 , &y[0] - 1 , ndim, ftol,
                funk, &me->nfunk);
        }
       
        /* Map vars from R to (a, b) */
        for (i = 0; i < n; ++i){
            const OneToOneMapping& thisMapping = (*mappings[i]);
            // the best vertex is slotted in the first slot, so return
            // that one
            outx[i] = thisMapping(vertices[0][i]);
        }

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

                     
const double Simplex::TINY = 1.0e-10; // A small number.
const int Simplex::NMAX_default = 5000;       // Maximum allowed number of function evaluations.

const double Simplex::EXCEPTIONALLY_LARGE_NUMBER = 1.0e+30;     // returned by obj func upon failure
Simplex* Simplex::me = 0;

double Simplex::funk(double inx[]){
    static const string method = "Simplex::funk";
    try{
        DoubleArray& x = me->x;
        DoubleArray& f = me->f;
        OneToOneMappingConstArray& mappings = me->mappings;
        const StringArray& ids = *me->ids;
        int n = me->ndim;
        // map vars
        int i = 0;
        for (; i < n; ++i){
            const OneToOneMapping& thisMapping = *mappings[i];
            x[i] = thisMapping(inx[i+1]);  // map from R to (a, b)
        }
        try{
            // call objective func 
            (*me->func)(x, f);
        }
        catch(exception& e){
            if (me->exitUponFailure){
                if (ids.size() > 0){
                    string message = "Simplex::funk: Failed with values";
                    int i = 0;
                    for (; i < ids.size(); ++i){
                        message += "\n\t" + ids[i] + ": " + Format::toString(x[i]);
                    }
                    throw ModelException::addTextToException(e, message);
                }
                else{
                    throw e;
                }
            }
            else{
                return EXCEPTIONALLY_LARGE_NUMBER;
            }
        }
        return f[0];
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

#ifdef GET_PSUM 
#undef GET_PSUM 
#endif
# ifdef SWAP
#undef SWAP
#endif

#define GET_PSUM \
    for (j=1;j<=ndim;j++) {\
    for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
    psum[j]=sum;}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void Simplex::amoeba(double **p, double y[], int ndim, double ftol,
                      double (*funk)(double []), int *nfunk) const
    /** Multidimensional minimization of the function funk(x) where x[1..ndim] is a dvector in ndim
    dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[1..ndim+1][1..ndim] 
    is input. Its ndim+1 rows are ndim-dimensional dvectors which are the vertices of
    the starting simplex. Also input is the dvector y[1..ndim+1], whose components must be preinitialized
    to the values of funk evaluated at the ndim+1 vertices (rows) of p; and ftol the
    fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and
    y will have been reset to ndim+1 new points all within ftol of a minimum function value, and
    nfunk gives the number of function evaluations taken. */
{
    int i,ihi,ilo,inhi,j,mpts=ndim+1;
    double rtol,sum,swap,ysave,ytry,*psum;
    try{
        psum=dvector(1,ndim);
        *nfunk=0;
        GET_PSUM
            for (;;) {
                ilo=1;
                // First we must determine which point is the highest (worst), next-highest, and lowest
                // (best), by looping over the points in the simplex.
                ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
                for (i=1;i<=mpts;i++) {
                    if (y[i] <= y[ilo]) ilo=i;
                    if (y[i] > y[ihi]) {
                        inhi=ihi;
                        ihi=i;
                    } else if (y[i] > y[inhi] && i != ihi) inhi=i;
                }
                rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
               
                // Compute the fractional range from highest to lowest and return if satisfactory.
                if (rtol < ftol
                    || *nfunk >= Simplex::me->NMAX) {   // don't fail if max nb iters is breached. simply return
                        // If returning, put best point and value in slot 1.
                        SWAP(y[1],y[ilo])
                            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]) 
                                break;
                    }
            
                *nfunk += 2;
                // Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
                // across from the high point, i.e., reflect the simplex from the high point.
                ytry=amotry(p,y,psum,ndim,funk,ihi,ilo,-1.0);

                if (ytry <= y[ilo])
                    // Gives a result better than the best point, so try an additional extrapolation by a
                    // factor 2 
                        ytry=amotry(p,y,psum,ndim,funk,ihi,ilo,2.0);
                else if (ytry >= y[inhi]) {

                    // The reflected point is worse than the second-highest, so look for an intermediate
                    // lower point, i.e., do a one-dimensional contraction.
                    ysave=y[ihi];
                    ytry=amotry(p,y,psum,ndim,funk,ihi,ilo,0.5);
                    if (ytry >= ysave) {    // Can?t seem to get rid of that high point. Better
                        for (i=1;i<=mpts;i++) { // contract around the lowest (best) point.
                            if (i != ilo) {
                                for (j=1;j<=ndim;j++)
                                    p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                                y[i]=(*funk)(psum);
                            }
                        }
                        *nfunk += ndim; // Keep track of function evaluations.
                        GET_PSUM        // Recompute psum.
                    }
                } else --(*nfunk);// Correct the evaluation count.
                                  // Go back for the test of doneness and the next iteration. 
            }                      
            free_dvector(psum,1,ndim);
    }
    catch(exception& e){
        free_dvector(psum,1,ndim);
        throw e;
    }
}

void Simplex::amoeba2(double **p, double y[], int ndim, double ftol,double ftolAbs,
                     double (*funk)(double []), int *nfunk, int idx) const
/** Multidimensional minimization of the function funk(x) where x[1..ndim] is a dvector in ndim
    dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[1..ndim+1][1..ndim] 
    is input. Its ndim+1 rows are ndim-dimensional dvectors which are the vertices of
    the starting simplex. Also input is the dvector y[1..ndim+1], whose components must be preinitialized
    to the values of funk evaluated at the ndim+1 vertices (rows) of p; and ftol the
    fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and
    y will have been reset to ndim+1 new points all within ftol of a minimum function value, and
    nfunk gives the number of function evaluations taken. */
{
    int i,ihi,ilo,inhi,j,mpts=ndim+1;
    double facLow, facUp; 
    double tolGuess = 1.0E-05;
    int maxNbIter = 3;
    int indexsave = 0;
    double maxTol = 4.0E-06;
    double rtolAbs,rtolRel,sum,swap,ysave,ytry,*psum;
    try{
        psum=dvector(1,ndim);
        *nfunk=0;
        GET_PSUM
        for (;;) {
            ilo=1;
            // First we must determine which point is the highest (worst), next-highest, and lowest
            // (best), by looping over the points in the simplex.
            ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
            for (i=1;i<=mpts;i++) {
                if (y[i] <= y[ilo]) ilo=i;
                if (y[i] > y[ihi]) {
                    inhi=ihi;
                    ihi=i;
                } else if (y[i] > y[inhi] && i != ihi) inhi=i;
            }
            rtolRel=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
            rtolAbs=fabs(y[ilo]);
            
            if(idx == 1){
                if(rtolAbs > tolGuess){
                    facLow = -1.0;
                    facUp = 2.0;
                }
                else{
                    facLow = -0.5;
                    facUp = 1.5;
                }
            }
            idx += 1;
            
            if(rtolAbs < maxTol){
                indexsave += 1;
            }

          
            if (CString::equalsIgnoreCase(stoppingCriterionType,Simplex::RELATIVE)){
                 // Compute the fractional range from highest to lowest and return if satisfactory.
                 if (rtolRel < ftol
                     || *nfunk >= Simplex::me->NMAX) {   // don't fail if max nb iters is breached. simply return
                         // If returning, put best point and value in slot 1.
                         SWAP(y[1],y[ilo])
                             for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]) 
                                 break;
                     }
            }
            if (CString::equalsIgnoreCase(stoppingCriterionType,Simplex::ABSOLUTE)){
                // Compute the fractional range from highest to lowest and return if satisfactory.
                if (rtolAbs < ftolAbs  
                    || *nfunk >= Simplex::me->NMAX) {   // don't fail if max nb iters is breached. simply return
                        // If returning, put best point and value in slot 1.
                        SWAP(y[1],y[ilo])
                             for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]) 
                                break;
                    }
                if(indexsave == maxNbIter){
                       SWAP(y[1],y[ilo])
                            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]) 
                                break;
                    }
            }
            if (CString::equalsIgnoreCase(stoppingCriterionType,Simplex::COMPOSITE)){
                // Compute the fractional range from highest to lowest and return if satisfactory.
                if (rtolAbs < ftolAbs || rtolRel < ftol
                    || *nfunk >= Simplex::me->NMAX) {   // don't fail if max nb iters is breached. simply return
                        // If returning, put best point and value in slot 1.
                        SWAP(y[1],y[ilo])
                            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]) 
                                break;
                    }
                if(indexsave == maxNbIter){
                        SWAP(y[1],y[ilo])
                            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i]) 
                                break;
                    }
            }
            *nfunk += 2;
            // Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
            // across from the high point, i.e., reflect the simplex from the high point.
            ytry=amotry(p,y,psum,ndim,funk,ihi,ilo,facLow);

            if (ytry <= y[ilo]){
                // Gives a result better than the best point, so try an additional extrapolation by a
                // factor 2 
                if(fabs(y[ilo]) > ftolAbs){
                    ytry=amotry(p,y,psum,ndim,funk,ihi,ilo,facUp);
                }
            }
            else if (ytry >= y[inhi]) {

                // The reflected point is worse than the second-highest, so look for an intermediate
                // lower point, i.e., do a one-dimensional contraction.
                ysave=y[ihi];
                ytry=amotry(p,y,psum,ndim,funk,ihi,ilo,0.5);
                if (ytry >= ysave) {    // Can?t seem to get rid of that high point. Better
                    for (i=1;i<=mpts;i++) { // contract around the lowest (best) point.
                        if (i != ilo) {
                            for (j=1;j<=ndim;j++)
                                p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                            y[i]=(*funk)(psum);
                        }
                    }
                    *nfunk += ndim; // Keep track of function evaluations.
                    GET_PSUM        // Recompute psum.
                }
            } else --(*nfunk);// Correct the evaluation count.
                              // Go back for the test of doneness and the next iteration. 
        }                      
        free_dvector(psum,1,ndim);
    }
    catch(exception& e){
        free_dvector(psum,1,ndim);
        throw e;
    }
}

#undef GET_PSUM
#undef SWAP

double Simplex::amotry(double **p, double y[], double psum[], int ndim,
                       double (*funk)(double []), int ihi, int ilo, double fac)
/** Extrapolates by a factor fac through the face of the simplex across from the high point, tries
    it, and replaces the high point if the new point is better.*/
{
    int j;
    double fac1,fac2,ytry,*ptry;
  
    try{
        ptry=dvector(1,ndim);
        fac1=(1.0-fac)/ndim;
        fac2=fac1-fac;
        for (j=1;j<=ndim;j++){
            ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
        }
        ytry=(*funk)(ptry);     // Evaluate the function at the trial point.
        if (ytry < y[ihi]) {    // If it?s better than the highest, then replace the highest.
            y[ihi]=ytry;
            for (j=1;j<=ndim;j++) {
                psum[j] += ptry[j]-p[ihi][j];
                p[ihi][j]=ptry[j];
            }
        }
        free_dvector(ptry,1,ndim);
        return ytry;
    }
    catch(exception& e){
        free_dvector(ptry,1,ndim);
        throw e;
    }
}

Simplex::Simplex():
OptimizerND(TYPE),
exitUponFailure(false),
maxNbIter(NMAX_default),
ftol(sqrt(DBL_EPSILON)), // just a bit bigger than sqrt of machine precision by default
ftolAbs(0.000002),
stoppingCriterionType(Simplex::RELATIVE),
lengthScale(1.0),
basisMaker(new CanonicalBasisMaker()){}

void Simplex::validatePop2Object() {
    static const string method = "Simplex::validatePop2Object";
    NMAX = maxNbIter;

    // For overnight intraweek testing allow fewer iterations via
    // command line flag
    bool doingFewIter = 
        CommandLineParams::hasParameter(CommandLineParams::FewIter);
    if(doingFewIter)
    {
        // only do at most 2 function evaluations
        NMAX = 2;
        maxNbIter = 2;
    }
}

class SimplexHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Simplex, clazz);
        SUPERCLASS(OptimizerND);
        EMPTY_SHELL_METHOD(defaultCtor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        FIELD(exitUponFailure, "Whether to terminate if the evaluation of the obj func fails");
        FIELD_MAKE_OPTIONAL(exitUponFailure);
        FIELD(maxNbIter, "Maximum allowed number of function evaluations");
        FIELD_MAKE_OPTIONAL(maxNbIter);
        FIELD(NMAX, "Maximum allowed number of function evaluations");
        FIELD_MAKE_TRANSIENT(NMAX);
        FIELD(ftol, "Function tolerance");
        FIELD_MAKE_OPTIONAL(ftol);
        FIELD(lengthScale, "Length scale of the problem");
        FIELD_MAKE_OPTIONAL(lengthScale);
        FIELD(stoppingCriterionType, "stopping criterion in the simplex");
        FIELD_MAKE_OPTIONAL(stoppingCriterionType);
        FIELD(ftolAbs, "absolute tolerance level");
        FIELD_MAKE_OPTIONAL(ftolAbs);
        FIELD(basisMaker, "Basis used to generate n vertices of the simplex at start");
        FIELD_MAKE_OPTIONAL(basisMaker);
    }

    static IObject* defaultCtor(){
        return new Simplex();
    }
};


CClassConstSP const Simplex::TYPE = CClass::registerClassLoadMethod(
    "Simplex", typeid(Simplex), SimplexHelper::load);

/** alternative to public constructor */
SimplexSP Simplex::create(double _lengthScale)
{
	SimplexSP out = SimplexSP(new Simplex());
	out->lengthScale = _lengthScale;
	out->validatePop2Object();
	return out;
}
double Simplex::getLengthScale(){
    return lengthScale;
}

void Simplex::setAbsTol(double tol)
{
    ftolAbs = tol;
}

void Simplex::setLengthScale(double _lengthScale){
    lengthScale = _lengthScale;
}

void Simplex::useMyMapping(string _useNewMapping){
    useNewMapping = _useNewMapping;
}

void Simplex::setStoppingCriterionType(string _stoppingCriterionType){
    stoppingCriterionType = _stoppingCriterionType;
}

void Simplex::useSRMSimplex(string choose){
    isSRM = choose;
}
// OPTIMIZER 1D
CClassConstSP const Optimizer1D::TYPE = CClass::registerClassLoadMethod(
    "Optimizer1D", typeid(Optimizer1D), Optimizer1D::load);

Optimizer1D::Optimizer1D(const CClassConstSP& clazz):
CObject(clazz){}

void Optimizer1D::load(CClassSP& clazz){
    REGISTER(Optimizer1D, clazz);
    SUPERCLASS(CObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// OPTMIZER 1D BRACKETER
const Function1DDouble* Optimizer1D::Bracketer::func = 0;
OneToOneMappingConstSP Optimizer1D::Bracketer::mapping = OneToOneMappingConstSP( );

double Optimizer1D::Bracketer::local_func(double inx){
    double x = (*mapping)(inx);
    return (*func)(x);
}

/** Wrapper around numerical recipe's mnbrak */
void Optimizer1D::Bracketer::bracket(const Function1DDouble& infunc,
                                     double&                 inxa,
                                     double&                 inxb,
                                     double&                 inxc,
                                     double&                 fa,
                                     double&                 fb,
                                     double&                 fc){
    static const string method("Optimizer1D::Bracketer::bracket");
    try{
        func = &infunc;
        mapping = OneToOneMapping::create(infunc.getInterval());
        double xa = mapping->inverse(inxa);
        double xb = mapping->inverse(inxb);
        double xc = mapping->inverse(inxc);
        mnbrak(&xa, &xb, &xc, &fa, &fb, &fc, local_func);
        inxa = (*mapping)(xa);
        inxb = (*mapping)(xb);
        inxc = (*mapping)(xc);
        func = 0;
    }
    catch(exception& e){
        func = 0;
        throw ModelException(e, method);
    }
}

// OPTIMIZER 1D BRENT
Optimizer1DBrent* Optimizer1DBrent::me = 0;

CClassConstSP const Optimizer1DBrent::TYPE = CClass::registerClassLoadMethod(
    "Optimizer1DBrent", typeid(Optimizer1DBrent), Optimizer1DBrent::load);

Optimizer1DBrent::Optimizer1DBrent():
Optimizer1D(TYPE){}

Optimizer1DBrent::Optimizer1DBrent(double tol):
Optimizer1D(TYPE),
tol(tol){}

void Optimizer1DBrent::load(CClassSP& clazz){
    REGISTER(Optimizer1DBrent, clazz);
    SUPERCLASS(Optimizer1D);
    clazz->setPublic(); // make visible to EAS/spreadsheet
    FIELD(tol, "Relative Tolerance in x")
    FIELD(fmin, "")
    FIELD_MAKE_TRANSIENT(fmin)
}

double Optimizer1DBrent::local_func(double inx){
    double x = (*me->mapping)(inx);
    return (*me->func)(x);
}

double Optimizer1DBrent::minimize(const Function1DDouble& infunc,
                                  double                  inxa,
                                  double                  inxb,
                                  double                  inxc) const{
    static const string method("Optimizer1DBrent::minimize");
    try{
        me = const_cast<Optimizer1DBrent*>(this);
        me->func = &infunc;
        me->mapping = OneToOneMapping::create(infunc.getInterval());
        double xa = me->mapping->inverse(inxa);
        double xb = me->mapping->inverse(inxb);
        double xc = me->mapping->inverse(inxc);
        double mappedtol = tol * inxb / Maths::max(1.0, fabs(xb * mapping->derivative(xb)));
        double xmin = 0.0;
        fmin = brent(xa, xb, xc, local_func, mappedtol, &xmin);
        xmin = (*me->mapping)(xmin);
        me = 0;
        return xmin;
    }
    catch(exception& e){
        me = 0;
        throw ModelException(e, method);
    }
}


DRLIB_END_NAMESPACE
