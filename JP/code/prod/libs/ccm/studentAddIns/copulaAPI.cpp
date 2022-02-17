#include <General/General.h>
#include <Magnet/Magnet.h>
#include <list>
extern "C" {
#include "../delta.h"
//#include "../skew.h"
#include "../skew_util.h"
#include "../beta.h"
#include "../student.h"
#include "../gaussian.h"
#include "../dependence.h"
#include "../independence.h"
#include "../proba_utils.h"
//#include "../payoff_tp.h"
#include "../alpha_stable_random.h"
#include "../gamma_random.h"
#include "../convolution.h"
#include "../morphing.h"
#include "../clayton.h"
#include "../gumbel.h"
}


using namespace CM;
using namespace std;
MAGNET_PREFIX("DR2_")
MAGNET_CATEGORY("Copula Utils tp")
// --------------------------------------------------------------------------
// DRError message handling
static list<String> msgList;
static String errmsg() {
    String out;
    for (list<String>::iterator i=msgList.begin(); i!=msgList.end(); i++)
        out += *i + '\n';
    msgList.clear();
    return out;
}
#define DR_ERROR_CHECK {if (!msgList.empty()) CM_THROW RuntimeError(errmsg());}

extern "C" {
void DR_ErrorCallBack(const char* format, va_list args)
{
    static char msg[1024];
    vsprintf(msg,format,args);
    msgList.push_back(String(msg));
}

typedef void (*ErrFuncPtr)(const char* , va_list);
extern ErrFuncPtr errCallBack;
}
// setup error callback when lib is initialized.
static struct dummy {dummy(){errCallBack = DR_ErrorCallBack;}} s_initLib;

// --------------------------------------------------------------------------
// Magnet to C adaptor code
//

/* 
 * When converting from C struct to C++ object
 * - the C struct contains a pointer to C++ owned data
 * - the C++ class can be constructed from a C struct, but this requires 
 *   a deep copy
 * Conclusion
 * - C++ struct always owns its vectors/matrix and destroys them
 * - C struct own its vector when it is not create from C++ class, in
 *   this case (and this one only) it needs be freed).
 */

/** define an indicator sim object that can be mapped into C struct */
struct IndicatorSim : public CM::Object
{
    Matrix<int> indicator; //[nbPath][nbNames]
    Array<double> weight;
};

inline INDICATOR_SIM IndicatorSim2Struct(const IndicatorSim &is)
{
    INDICATOR_SIM i_s;
    i_s.indicator = const_cast<int *>(&is.indicator[0][0]);
    i_s.nbNames = is.indicator.colSize();
    i_s.nbPaths = is.indicator.rowSize();
    i_s.weight = const_cast<double *>(&is.weight[0]);
    return i_s;
}

inline IndicatorSim Struct2IndicatorSim(const INDICATOR_SIM &i_s)
{
    IndicatorSim is;
    ConstMatrixView<int> m(i_s.indicator,i_s.nbPaths,i_s.nbNames);
    is.indicator = m;
    is.weight = Array<double>(i_s.weight,i_s.weight+i_s.nbPaths);
    return is;
}

/** define a credit portfolio object that can be mapped into C struct */
struct CreditPortfolio : public CM::Object
{
    Array<double> notional;
    Array<double> recovery;
    Array<double> nameMaturity;
    double strike1;
    double strike2;
};

inline CREDIT_PORTFOLIO CreditPortfolio2Struct(const CreditPortfolio &cp)
{
    CREDIT_PORTFOLIO c_p;
    c_p.nbNames = cp.notional.size();
    c_p.notional = const_cast<double *>(&cp.notional[0]);
    c_p.recovery = const_cast<double *>(&cp.recovery[0]);
    c_p.nameMaturity = const_cast<double *>(&cp.nameMaturity[0]);
    c_p.strike1 = cp.strike1;
    c_p.strike2 = cp.strike2;
    return c_p;
}

inline CreditPortfolio Struct2CreditPortfolio(CREDIT_PORTFOLIO &c_p)
{
    CreditPortfolio cp;
    cp.notional = Array<double>(c_p.notional, c_p.notional+c_p.nbNames);
    cp.recovery = Array<double>(c_p.recovery, c_p.recovery+c_p.nbNames);
    cp.nameMaturity = Array<double>(c_p.nameMaturity, c_p.nameMaturity+c_p.nbNames);
    cp.strike1 = c_p.strike1;
    cp.strike2 = c_p.strike2;
    return cp;
}

/* 
 * When converting from C struct to C++ object, a different route was chosen
 * here because some type dependent logic for alloc/free exist in the C code.
 * - the C struct owns some private data, and only CopulaFree can free it.
 * - the C++ class owns the C struct, and calls CopulaFree on
 *   destruction. This scheme forces us to preclude C++ copy and assignment.
 * Conclusion:
 * - C struct always own its param pointer, and frees it with CopulaFree
 * - C++ struct always ends up owning everything,its dtor takes care of
 *   freeing the C struct
 */

/** define a copula object that can be mapped into C struct */
class Copula : public CM::Object
{
public:
    COPULA _copula; 
    // ctor and dtor memory mgt of the C struct.
    Copula(const COPULA &c):_copula(c) {}
    ~Copula() {CopulaFree(&_copula);}
private:
    // disable copy ctor and assgn operator
    Copula operator=(const Copula &a);
    Copula(const Copula &a);
};

// --------------------------------------------------------------------------
// MAGNET code
//
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// CreditPortfolioCreate
//
MAGNET_DESCRIPTION("returns a credit portfolio object")
MAGNET_PARAMETER("",notional,"","Array of the names' notionals")
MAGNET_PARAMETER("",recovery,"","Array of the names' recoveries")
MAGNET_PARAMETER("",nameMaturity,"","Array of the names' maturities")
MAGNET_PARAMETER("",strike1,"","lower strike")
MAGNET_PARAMETER("",strike2,"","upper strike")
MAGNET_FUNCTION5(SharedPointer<CreditPortfolio>, CreditPortfolioCreate,
                   const Array<double> &notional, 
                   const Array<double> &recovery, 
                   const Array<double> &nameMaturity, 
                   double  strike1, 
                   double  strike2)
{
    SharedPointer<CreditPortfolio> cp = Raw2SharedPointer(
        new CreditPortfolio());
    cp->notional = notional;
    cp->recovery = recovery;
    cp->nameMaturity = nameMaturity;
    cp->strike1 = strike1;
    cp->strike2 = strike2;
    return cp;
}

// --------------------------------------------------------------------------
// CopulaCreate
//
MAGNET_DESCRIPTION("returns a copula object")
MAGNET_PARAMETER("",type,"","one of Gaussian, Student, Dependence, Independence, Gumble, Clayton or q")
MAGNET_PARAMETER("",beta,"","vector of weights (or alpha for gumbel)")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom ofthe student variables")
MAGNET_X_FUNCTION14(SharedPointer<Copula>,
                 CopulaCreate,
                 String type, MAGNET_MANDATORY,
                 long nbNames, MAGNET_MANDATORY,
                 const Array<double> &beta, =Array<double>(),
                 long freedomDegree, = 0,
                 double theta, = 1.0,
                 const Array<double> &Msample, =Array<double>(),
                 const Array<double> &chi2Sample, =Array<double>(),
                 long seed, = -7,
                 double qM, = 0.0,
                 double qZ, = 0.0,
                 double xa, = 0.0,
                 double xb, = 0.0,
                 double xc, = 0.0,
                 double delta, = 0.0)
{
    // map types
    SymbolMap<long> theType(
        C_GAUSS, "gaussian",
        C_STUDENT, "student",
        C_DEPENDENCE, "dependence",
        C_INDEPENDENCE, "independence",
        C_GUMBEL, "gumbel",
        C_CLAYTON, "clayton",
        C_SKEW, "q",
        C_DELTA, "delta");
  	long pos = -1;
    int status = FAILURE;
	theType.name2Value(type, pos);
    if (pos < 0 ) 
        CM_THROW RuntimeError("type must be either gaussian, student, dependence, independence, gumbel, clayton, q or delta");

    COPULA c = CopulaCreate((CopulaType)pos,nbNames,
        const_cast<double *>(&beta[0]),
                            freedomDegree,
                            theta,
                            seed,
                            qM,
                            qZ,
                            xa,xb,xc,delta,
                            Msample.size() ? &Msample[0] : (const double *)0,
                            Msample.size(),
                            chi2Sample.size() ? &chi2Sample[0] : (const double *)0,
                            chi2Sample.size(),
                            &status);
    DR_ERROR_CHECK
    // shallow copy
    return Raw2SharedPointer(new Copula(c));
}

// --------------------------------------------------------------------------
// CopulaProductCreate
//
MAGNET_DESCRIPTION("returns a copula object, product of 2 copulas")
MAGNET_PARAMETER("",copula1,"","copula object created by CopulaCreate")
MAGNET_PARAMETER("",copula2,"","copula object created by CopulaCreate")
MAGNET_PARAMETER("",alpha,"","array[nbNames] of the copula weights")
MAGNET_PARAMETER("",nbNames,"","numberof names in the credit portfolio")
MAGNET_FUNCTION4(SharedPointer<Copula>,
                 CopulaProductCreate,
                 const Copula &copula1,
                 const Copula &copula2,
                 const Array<double> &alpha,
                 long nbNames)
{
    COPULA c = CopulaProductCreate(&copula1._copula,&copula2._copula,&alpha[0],nbNames);
    DR_ERROR_CHECK
    // shallow copy
    return Raw2SharedPointer(new Copula(c));
}

// --------------------------------------------------------------------------
// Copula3ProductCreate
//
MAGNET_DESCRIPTION("returns a copula object, product of 3 copulas")
MAGNET_PARAMETER("",copula1,"","copula object created by CopulaCreate")
MAGNET_PARAMETER("",copula2,"","copula object created by CopulaCreate")
MAGNET_PARAMETER("",alpha1,"","array[nbNames] of the 1st copula weights")
MAGNET_PARAMETER("",alpha2,"","array[nbNames] of the last copula weights")
MAGNET_PARAMETER("",nbNames,"","numberof names in the credit portfolio")
MAGNET_FUNCTION6(SharedPointer<Copula>,
                 Copula3ProductCreate,
                 const Copula &copula1,
                 const Copula &copula2,
                 const Copula &copula3,
                 const Array<double> &alpha1,
                 const Array<double> &alpha2,
                 long nbNames)
{
    COPULA c = Copula3ProductCreate(&copula1._copula,
                                    &copula2._copula,
                                    &copula3._copula,
                                    &alpha1[0],
                                    &alpha2[0],
                                    nbNames);
    // shallow copy
    DR_ERROR_CHECK
    return Raw2SharedPointer(new Copula(c));
}

// --------------------------------------------------------------------------
// IndicatorSimCreate
//
MAGNET_DESCRIPTION("returns an indicator object")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_PARAMETER("",survivalProba,"","array[nbNames] of survival probabilities")
MAGNET_PARAMETER("",copula,"","copula object created by CopulaCreate")
MAGNET_X_FUNCTION3(SharedPointer<IndicatorSim>, IndicatorSimCreate,
                 long nbPaths, MAGNET_MANDATORY,
                 const Array<double> &survivalProba, MAGNET_MANDATORY,
                 const Copula &copula, MAGNET_MANDATORY)
{
    long nbNames = survivalProba.size();
    INDICATOR_SIM i_s = 
        IndicatorSimCreate( nbPaths,
                            nbNames,
                            &survivalProba[0],
                            &copula._copula);
    DR_ERROR_CHECK

    SharedPointer<IndicatorSim> is = Raw2SharedPointer(new IndicatorSim());
    // deep copy
    *is = Struct2IndicatorSim(i_s);
    IndicatorSimFree(&i_s);
    DR_ERROR_CHECK
    return is;
}

// --------------------------------------------------------------------------
// IndicatorSimCreateMtp
//
MAGNET_DESCRIPTION("returns a mlti time point indicator object")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_PARAMETER("",survivalProba,"","array[nbNames] of survival probabilities")
MAGNET_PARAMETER("",copula,"","copula object created by CopulaCreate")
MAGNET_X_FUNCTION3(SharedPointer<IndicatorSim>, IndicatorSimCreateMtp,
                 long nbPaths, MAGNET_MANDATORY,
                 const Matrix<double> &survivalProba, MAGNET_MANDATORY,
                 const Copula &copula, MAGNET_MANDATORY)
{
    long nbNames = survivalProba.rowSize();
    long nbTimes = survivalProba.colSize();
    INDICATOR_SIM i_s = 
        IndicatorSimCreate_mtp( nbPaths,
                            nbNames,
                            nbTimes,
                            &survivalProba[0][0],
                            &copula._copula);
    DR_ERROR_CHECK

    SharedPointer<IndicatorSim> is = Raw2SharedPointer(new IndicatorSim());
    // deep copy
    *is = Struct2IndicatorSim(i_s);
    IndicatorSimFree(&i_s);
    DR_ERROR_CHECK
    return is;
}

// --------------------------------------------------------------------------
// ExpectedPayoff
//
MAGNET_DESCRIPTION("returns the expected payoff of a credit portfolio")
MAGNET_PARAMETER("",port,"","credit portfolio structure created by CreditPortfolioCreate")
MAGNET_PARAMETER("",indicatorSim,"","Indicator simulation structure created by IndicatorSimCreate")
MAGNET_FUNCTION2(double, ExpectedPayoff,
                 const CreditPortfolio &creditPortfolio,
                 const IndicatorSim &indicatorSim)
{
    double res;
    CREDIT_PORTFOLIO c_r = CreditPortfolio2Struct(creditPortfolio);
    INDICATOR_SIM i_s = IndicatorSim2Struct(indicatorSim);
    res = ExpectedPayoff_tp(&c_r, &i_s);
    DR_ERROR_CHECK
    return res;
}

// --------------------------------------------------------------------------
// PayoffDistribution
//
MAGNET_FUNCTION3(Array<double>, PayoffDistribution,
                            const CreditPortfolio &creditPortfolio,
                            const IndicatorSim &indicatorSim,
                            const Array<double> &sampleLoss)
{
    long nbSampleLoss = sampleLoss.size();
    int status = FAILURE;
    CREDIT_PORTFOLIO c_r = CreditPortfolio2Struct(creditPortfolio);
    INDICATOR_SIM i_s = IndicatorSim2Struct(indicatorSim);
    Array<double> loss = Array<double>(nbSampleLoss);
    status = PayoffDistribution_tp( &c_r,
                                    &i_s,
                                    &sampleLoss[0],
                                    nbSampleLoss,
                                    &loss[0]);
    DR_ERROR_CHECK;
    return loss;

}

// --------------------------------------------------------------------------
// NbDefaultNameDistribution
//
MAGNET_DESCRIPTION("returns an array[0..nbNames] of discrete probability of k names defaulted")
MAGNET_PARAMETER("",indicatorSim,"","IndicatorSIm object create by IndicatorSimCreate")
MAGNET_FUNCTION1(Array<double>, NbDefaultNameDistribution,
                 const IndicatorSim &indicatorSim)
{
    INDICATOR_SIM i_s = IndicatorSim2Struct(indicatorSim);
    long nbNames = i_s.nbNames;
    Array<double> out = Array<double>(nbNames+1);
    NbDefaultNameDistribution( i_s, &out[0]);
    DR_ERROR_CHECK
    return out;
}

// --------------------------------------------------------------------------
// NbConditionalDefaultNameDistribution_mtp
//
MAGNET_DESCRIPTION("returns an array[0..nbNames] of discrete probability of k names defaulted")
MAGNET_PARAMETER("",indicatorSim,"","IndicatorSIm object create by IndicatorSimCreate")
MAGNET_PARAMETER("",idx_t0,"","index of time points after which we observe the first default")
MAGNET_PARAMETER("",idx_t1,"","index of time points before which we observe the first default")
MAGNET_PARAMETER("",idx_T0,"","index of time points after which we observe the conditional defaults")
MAGNET_PARAMETER("",idx_T1,"","index of time points before which we observe the conditional defaults")
MAGNET_FUNCTION5(Array<double>, NbConditionalDefaultNameDistribution_mtp,
                 const IndicatorSim &indicatorSim,
                 long idx_t0,
                 long idx_t1,
                 long idx_T0,
                 long idx_T1)
{
    INDICATOR_SIM i_s = IndicatorSim2Struct(indicatorSim);
    long nbNames = i_s.nbNames;
    Array<double> out = Array<double>(nbNames+1);
    NbConditionalDefaultNameDistribution_mtp( i_s, &out[0], idx_t0, idx_t1, idx_T0, idx_T1);
    DR_ERROR_CHECK
    return out;
}


// --------------------------------------------------------------------------
// AlphaStableDeviates
//
MAGNET_DESCRIPTION("Create an array of alpha stable deviates")
MAGNET_PARAMETER("",nbPaths,"","number of paths in the simulation, must be an integer >0")
MAGNET_PARAMETER("",c,"","dispersion parameter of the distribution")
MAGNET_PARAMETER("",alpha,"","characteristic exponent of the distribution")
MAGNET_FUNCTION3(Array<double>, AlphaStableDeviates,    
    long nbPaths,
    double c,
    double alpha)
{
    Array<double> out = Array<double>(nbPaths);
    AlphaStableDeviates(    &out[0],
                            nbPaths,
                            c,
                            alpha);
    DR_ERROR_CHECK
    return out;
}

// --------------------------------------------------------------------------
// GammaDeviates
//
MAGNET_DESCRIPTION("Create an array of alpha stable deviates")
MAGNET_PARAMETER("",nbPaths,"","number of paths in the simulation, must be an integer >0")
MAGNET_PARAMETER("",c,"","dispersion parameter of the distribution")
MAGNET_PARAMETER("",alpha,"","characteristic exponent of the distribution")
MAGNET_FUNCTION3(Array<double>, GammaDeviates,    
    long nbPaths,
    double alpha,
    double beta)
{
    Array<double> out = Array<double>(nbPaths);
    GammaDeviates(    &out[0],
                            nbPaths,
                            alpha,
                            beta);
    DR_ERROR_CHECK
    return out;
}

// --------------------------------------------------------------------------
// Convolution
//
MAGNET_DESCRIPTION("Create the convoluted array of two arrays")
MAGNET_PARAMETER("",v1,"","array 1")
MAGNET_PARAMETER("",v2,"","array 2")
MAGNET_FUNCTION2(Array<double>, Convolution,    
    const Array<double> &v1,
    const Array<double> &v2)
{
    long n = v1.size();
    int status;
    Array<double> out(n);
    status = Convolution(&v1[0],&v2[0],&out[0],n);
    DR_ERROR_CHECK
    return out;
}

// --------------------------------------------------------------------------
// Beta
//
MAGNET_FUNCTION2(double, Beta, double a, double b)
{
    return beta(a,b);
}

MAGNET_FUNCTION3(double, Beta_Distribution,
                 double x,
                 double a,
                 double b)
{
    return Beta_Distribution(x,a,b);
}

MAGNET_FUNCTION2(double, Beta_Mean, double a, double b)
{
    return Beta_Mean(a,b);
}

MAGNET_FUNCTION2(double, Beta_Variance, double a, double b)
{
    return Beta_Variance(a,b);
}

MAGNET_FUNCTION2(double, Beta_Skew, double a, double b)
{
    return Beta_Skew(a,b);
}

MAGNET_FUNCTION2(double, Beta_Kurtosis, double a, double b)
{
    return Beta_Kurtosis(a,b);
}

MAGNET_FUNCTION2(double, Beta_Mode, double a, double b)
{
    return Beta_Mode(a,b);
}

// --------------------------------------------------------------------------
// Moments
//
MAGNET_FUNCTION3(Array<double>,
                 Moments,
                 const Array<double> &probas,
                 const Array<double> &values,
                 long nbMoments)
{
    int status;
    long nbValues = probas.size();
    Array<double> out(nbMoments);
    status = Moments(&probas[0], &values[0],nbValues,&out[0],nbMoments);
    DR_ERROR_CHECK
    return out;
}

// --------------------------------------------------------------------------
// ComplexMatrix
//
/*
MAGNET_FUNCTION2(Matrix<double>, ComplexMatrix,
                 Array<double> &realPart,
                 Array<double> &imagPart)
{
    int status;
    long nbValues = realPart.size();
    Matrix<double> out = Matrix<double>(nbValues,nbValues);
    status = ComplexMatrix(&realPart[0], &imagPart[0],nbValues, &out[0][0]);
    return out;
}
*/

// --------------------------------------------------------------------------
// MatrixProduct
//

MAGNET_FUNCTION2(Matrix<double>, MatrixProduct,
                Matrix<double> &M1,
                Matrix<double> &M2)
{
    int status;
    long n = M1.colSize();
    Matrix<double> M(n,n);
    status = MatrixProduct(&M1[0][0], &M2[0][0],n,&M[0][0]);
    DR_ERROR_CHECK
    return M;
}


// --------------------------------------------------------------------------
// F
//
MAGNET_X_FUNCTION6(   double, F,
                    double u, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ, MAGNET_MANDATORY,
                    long nbPoints, = 1000,
                    double eps, = 1e-6)
{
    return F(u,beta,qM,qZ,nbPoints,eps,NULL);
}

// --------------------------------------------------------------------------
// Fq
//
MAGNET_FUNCTION2(double, Fq,
                double x,
                double q)
{
    return Fq(x,q);
}

// --------------------------------------------------------------------------
// Fqinv
//
MAGNET_FUNCTION2(double, Fqinv,
                double y,
                double q)
{
    return Fqinv(y,q);
}

// --------------------------------------------------------------------------
// fq
//
MAGNET_FUNCTION2(double, fq,
                double x,
                double q)
{
    return fq(x,q);
}

// --------------------------------------------------------------------------
// f1q
//
MAGNET_FUNCTION2(double, f1q,
                double x,
                double q)
{
    return f1q(x,q);
}

// --------------------------------------------------------------------------
// Loss
//
MAGNET_X_FUNCTION8(double, Loss,
                 double K, MAGNET_MANDATORY,
                 double u, MAGNET_MANDATORY,
                 double beta, MAGNET_MANDATORY,
                 double qM, MAGNET_MANDATORY,
                 double qZ, MAGNET_MANDATORY,
                 long nbPoints, = 1000,
                 double eps, = 1e-6,
                 double alpha, = 0.0)
{
    return Loss(K,u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// DLoss
//
MAGNET_X_FUNCTION8(double, DLoss,
                 double K, MAGNET_MANDATORY,
                 double u, MAGNET_MANDATORY,
                 double beta, MAGNET_MANDATORY,
                 double qM, MAGNET_MANDATORY,
                 double qZ, MAGNET_MANDATORY,
                 long nbPoints, = 1000,
                 double eps, = 1e-6,
                 double alpha, = 0.0)
{
    return DLoss(K,u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Mode
//
MAGNET_X_FUNCTION7(double, Mode,
                    double u, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ, MAGNET_MANDATORY,
                    long nbPoints, =1000,
                    double eps, =1e-6,
                    double alpha, = 0.0)
{
    return Mode(u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// FMode
//
MAGNET_X_FUNCTION7(double, FMode,
                    double u, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ, MAGNET_MANDATORY,
                    long nbPoints, =1000,
                    double eps, =1e-6,
                    double alpha, =0.)
{
    return FMode(u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Moment
//
MAGNET_X_FUNCTION8(double, Moment,
                    double u, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ, MAGNET_MANDATORY,
                    long   n, = 1,
                    long nbPoints, =1000,
                    double eps, =1e-6,
                    double alpha, = 0.)
{
    return Moment(u,beta,qM,qZ,alpha,n,nbPoints,eps);
}

// --------------------------------------------------------------------------
// IntegerandLn
//
MAGNET_X_FUNCTION6(double, IntegrandLn,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K,MAGNET_MANDATORY,
                double n,MAGNET_MANDATORY)
{
    return IntegrandLn( Finv(u,beta,qM,qZ,10000,1e-10),
                        beta,qM,qZ,K,n);
}

// --------------------------------------------------------------------------
// Kn
//
MAGNET_X_FUNCTION8(double, Kn,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, = 0.)
{
    return Kn(u,beta,qM,qZ,alpha,n,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Sigma
//
MAGNET_X_FUNCTION7(double, Sigma,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, = 0.)
{
    return Sigma(u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Skew
//
MAGNET_X_FUNCTION7(double, Skew,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, = 0.)
{
    return Skew(u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Kurtosis
//
MAGNET_X_FUNCTION7(double, Kurtosis,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, = 0.)
{
    return Kurtosis(u,beta,qM,qZ,alpha,nbPoints,eps);
}

// --------------------------------------------------------------------------
// SeniorProba
//
MAGNET_X_FUNCTION8(double, SeniorProba,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, =0.)
{
    return SeniorProba(u,beta,qM,qZ,alpha,K,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Quantile
//
MAGNET_X_FUNCTION8(double, Quantile,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double proba,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, = 0.)
{
    return Quantile(u,beta,qM,qZ,alpha,proba,nbPoints,eps);
}

// --------------------------------------------------------------------------
// TrancheletPrice
/*
MAGNET_X_FUNCTION8(double, TrancheletPrice,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double alpha, = 0.)
{
    return TrancheletPrice(u,beta,qM,qZ,alpha,K,nbPoints,eps);
}
*/

/***************************************************************************/
/***************************************************************************/
// --------------------------------------------------------------------------
// Derivatives of Kn, Sigma, Mode, FMode, Moment
//
MAGNET_X_FUNCTION8(double, DModeDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DModeDqM(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DModeDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DModeDqZ(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DModeDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DModeDbeta(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DModeDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DModeDalpha(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DFModeDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DFModeDqM(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DFModeDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DFModeDqZ(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DFModeDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DFModeDbeta(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DFModeDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DFModeDalpha(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSigmaDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSigmaDqM(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSigmaDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DSigmaDqZ(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSigmaDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSigmaDbeta(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSigmaDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSigmaDalpha(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSkewDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSkewDqM(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSkewDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DSkewDqZ(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSkewDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DSkewDbeta(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DSkewDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DSkewDalpha(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DKurtosisDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DKurtosisDqM(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DKurtosisDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DKurtosisDqZ(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DKurtosisDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DKurtosisDbeta(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DKurtosisDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DKurtosisDalpha(u,beta,qM,qZ,alpha,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DKnDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DKnDqM(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DKnDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DKnDqZ(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DKnDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DKnDbeta(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DKnDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DKnDalpha(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DMomentDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DMomentDqM(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DMomentDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DMomentDqZ(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DMomentDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DMomentDbeta(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DMomentDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double n, = 1,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DMomentDalpha(u,beta,qM,qZ,alpha,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DSeniorProbaDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DSeniorProbaDqM(u,beta,qM,qZ,alpha,K,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DSeniorProbaDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSeniorProbaDqZ(u,beta,qM,qZ,alpha,K,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DSeniorProbaDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSeniorProbaDbeta(u,beta,qM,qZ,alpha,K,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DSeniorProbaDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DSeniorProbaDalpha(u,beta,qM,qZ,alpha,K,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DQuantileDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double proba, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DQuantileDqM(u,beta,qM,qZ,alpha,proba,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DQuantileDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double proba, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, =0.)
{
    return DQuantileDqZ(u,beta,qM,qZ,alpha,proba,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DQuantileDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double proba, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DQuantileDbeta(u,beta,qM,qZ,alpha,proba,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DQuantileDalpha,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double proba, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01,
                double alpha, = 0.)
{
    return DQuantileDalpha(u,beta,qM,qZ,alpha,proba,nbPoints,eps,h);
}


/*
MAGNET_X_FUNCTION8(double, DTrancheletPriceDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01)
{
    return DTrancheletPriceDqM(u,beta,qM,qZ,K,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DTrancheletPriceDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01)
{
    return DTrancheletPriceDqZ(u,beta,qM,qZ,K,nbPoints,eps,h);
}

MAGNET_X_FUNCTION8(double, DTrancheletPriceDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01)
{
    return DTrancheletPriceDbeta(u,beta,qM,qZ,K,nbPoints,eps,h);
}
*/

MAGNET_X_FUNCTION9(double, DIntegrandLnDqM,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long n, =0,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01)
{
    return DIntegrandLnDqM(u,beta,qM,qZ,K,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DIntegrandLnDqZ,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long n, =0,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01)
{
    return DIntegrandLnDqZ(u,beta,qM,qZ,K,n,nbPoints,eps,h);
}

MAGNET_X_FUNCTION9(double, DIntegrandLnDbeta,
                double u,MAGNET_MANDATORY,
                double beta,MAGNET_MANDATORY,
                double qM,MAGNET_MANDATORY,
                double qZ,MAGNET_MANDATORY,
                double K, MAGNET_MANDATORY,
                long n, =0,
                long nbPoints, = 1000,
                double eps, = 1e-6,
                double h, = 0.01)
{
    return DIntegrandLnDbeta(u,beta,qM,qZ,K,n,nbPoints,eps,h);
}

/***************************************************************************/
/*  COMPONENT IMPLEMENTATION
/***************************************************************************/
MAGNET_FUNCTION10(    double, LossDensity,
                    double K,
                    double u,
                    double beta,
                    const Array<double> &muM,
                    const Array<double> &sigmaM,
                    const Array<double> &weightM,
                    const Array<double> &muZ,
                    const Array<double> &sigmaZ,
                    const Array<double> &weightZ,
                    long var_adjust
                    )
{
    long i;
    long nbCompM = muM.size();
    long nbCompZ = muZ.size();
    DENSITY_COMPONENT *dcompM = (DENSITY_COMPONENT *)
                                malloc(nbCompM*sizeof(DENSITY_COMPONENT));
    DENSITY_COMPONENT *dcompZ = (DENSITY_COMPONENT *)
                                malloc(nbCompZ*sizeof(DENSITY_COMPONENT));
    for(i=0;i<nbCompM;i++)
    {
        dcompM[i].mu = muM[i];
        dcompM[i].sigma = sigmaM[i];
        dcompM[i].weight = weightM[i];
    }

    for(i=0;i<nbCompZ;i++)
    {
        dcompZ[i].mu = muZ[i];
        dcompZ[i].sigma = sigmaZ[i];
        dcompZ[i].weight = weightZ[i];
    }

    FACTOR_DENSITY  fm;
    FACTOR_DENSITY  fz;
    X_DENSITY       fx;

    fm.n = nbCompM;
    fm.t = dcompM;
    fz.n = nbCompZ;
    fz.t = dcompZ;

    fx.var_adjust = var_adjust;
    fx.beta = beta;
    fx.fm = fm;
    fx.fz = fz;

    double result = LossDensity(K, u, &fx);

    free(dcompM);
    free(dcompZ);
    return result;

}

MAGNET_FUNCTION10(    double, TrancheletLoss,
                    double K,
                    double u,
                    double beta,
                    const Array<double> &muM,
                    const Array<double> &sigmaM,
                    const Array<double> &weightM,
                    const Array<double> &muZ,
                    const Array<double> &sigmaZ,
                    const Array<double> &weightZ,
                    long var_adjust
                    )
{
    long i;
    long nbCompM = muM.size();
    long nbCompZ = muZ.size();
    DENSITY_COMPONENT *dcompM = (DENSITY_COMPONENT *)
                                malloc(nbCompM*sizeof(DENSITY_COMPONENT));
    DENSITY_COMPONENT *dcompZ = (DENSITY_COMPONENT *)
                                malloc(nbCompZ*sizeof(DENSITY_COMPONENT));
    for(i=0;i<nbCompM;i++)
    {
        dcompM[i].mu = muM[i];
        dcompM[i].sigma = sigmaM[i];
        dcompM[i].weight = weightM[i];
    }

    for(i=0;i<nbCompZ;i++)
    {
        dcompZ[i].mu = muZ[i];
        dcompZ[i].sigma = sigmaZ[i];
        dcompZ[i].weight = weightZ[i];
    }

    FACTOR_DENSITY  fm;
    FACTOR_DENSITY  fz;
    X_DENSITY       fx;

    fm.n = nbCompM;
    fm.t = dcompM;
    fz.n = nbCompZ;
    fz.t = dcompZ;

    fx.var_adjust = var_adjust;
    fx.beta = beta;
    fx.fm = fm;
    fx.fz = fz;

    double result = TrancheletLoss(K, u, &fx);

    free(dcompM);
    free(dcompZ);

    return result;
}

MAGNET_FUNCTION4(   double, FactorDensity,
                    double x,
                    const Array<double> &mu,
                    const Array<double> &sigma,
                    const Array<double> &weight,
                    )
{
    long i;
    long nbComp = mu.size();
    DENSITY_COMPONENT *dcomp = (DENSITY_COMPONENT *)
                                malloc(nbComp*sizeof(DENSITY_COMPONENT));
    for(i=0;i<nbComp;i++)
    {
        dcomp[i].mu = mu[i];
        dcomp[i].sigma = sigma[i];
        dcomp[i].weight = weight[i];
    }

    FACTOR_DENSITY  f;

    f.n = nbComp;
    f.t = dcomp;

    double result = FactorDensity(x, &f);

    free(dcomp);

    return result;
}

MAGNET_FUNCTION4(   double, FactorCum,
                    double x,
                    const Array<double> &mu,
                    const Array<double> &sigma,
                    const Array<double> &weight,
                    )
{
    long i;
    long nbComp = mu.size();
    DENSITY_COMPONENT *dcomp = (DENSITY_COMPONENT *)
                                malloc(nbComp*sizeof(DENSITY_COMPONENT));
    for(i=0;i<nbComp;i++)
    {
        dcomp[i].mu = mu[i];
        dcomp[i].sigma = sigma[i];
        dcomp[i].weight = weight[i];
    }

    FACTOR_DENSITY  f;

    f.n = nbComp;
    f.t = dcomp;

    double result = FactorCum(x, &f);

    free(dcomp);

    return result;
}

MAGNET_FUNCTION9(   double, F_X,
                    double x,
                    double beta,
                    const Array<double> &muM,
                    const Array<double> &sigmaM,
                    const Array<double> &weightM,
                    const Array<double> &muZ,
                    const Array<double> &sigmaZ,
                    const Array<double> &weightZ,
                    long var_adjust
                    )
{
    long i;
    long nbCompM = muM.size();
    long nbCompZ = muZ.size();
    DENSITY_COMPONENT *dcompM = (DENSITY_COMPONENT *)
                                malloc(nbCompM*sizeof(DENSITY_COMPONENT));
    DENSITY_COMPONENT *dcompZ = (DENSITY_COMPONENT *)
                                malloc(nbCompZ*sizeof(DENSITY_COMPONENT));
    for(i=0;i<nbCompM;i++)
    {
        dcompM[i].mu = muM[i];
        dcompM[i].sigma = sigmaM[i];
        dcompM[i].weight = weightM[i];
    }

    for(i=0;i<nbCompZ;i++)
    {
        dcompZ[i].mu = muZ[i];
        dcompZ[i].sigma = sigmaZ[i];
        dcompZ[i].weight = weightZ[i];
    }

    FACTOR_DENSITY  fm;
    FACTOR_DENSITY  fz;
    X_DENSITY       fx;

    fm.n = nbCompM;
    fm.t = dcompM;
    fz.n = nbCompZ;
    fz.t = dcompZ;

    fx.var_adjust = var_adjust;
    fx.beta = beta;
    fx.fm = fm;
    fx.fz = fz;

    double result = F_X(x, &fx);

    free(dcompM);
    free(dcompZ);

    return result;
}

MAGNET_FUNCTION6(double, DpiDwM,
                double K0,
                double K1,
                double u,
                double beta,
                double muM,
                double sigmaM)
{
    return DpiDwM( K0, K1, u, beta, muM, sigmaM);
}

MAGNET_FUNCTION6(double, DpiDwZ,
                double K0,
                double K1,
                double u,
                double beta,
                double muZ,
                double sigmaZ)
{
    return DpiDwZ( K0, K1, u, beta, muZ, sigmaZ);
}

MAGNET_FUNCTION4(double, DpiDbeta,
                double K0,
                double K1,
                double u,
                double beta)
{
    return DpiDbeta( K0, K1, u, beta);
}


/***************************************************************************/
/***************************************************************************/

//---------------------------------------------------------------------------
// Steps
//
MAGNET_X_FUNCTION6( Array<double>, Steps,
                    double u, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ, MAGNET_MANDATORY,
                    long nbPoints, =1000,
                    double eps, = 1e-8)
{
    DEBUGINFO *debugInfo = NULL;
    Array<double> xp(nbPoints);
    Array<double> yp2(nbPoints);
    long nbp = 0;
    double **yp = NULL;
    debugInfo = (DEBUGINFO*) malloc(sizeof(DEBUGINFO));
    if(debugInfo==NULL) goto RETURN;
    yp = (double**)malloc(sizeof(double**));
    if(yp==NULL) goto RETURN;
    debugInfo->kmax = nbPoints;
    debugInfo->dxsav = 0.0;
    debugInfo->nbPoints = nbp;
    debugInfo->xp = &xp[0];
    *yp = &yp2[0];
    debugInfo->yp = yp;
    F(u,beta,qM,qZ,nbPoints,eps,debugInfo);
RETURN:
    DR_ERROR_CHECK
    if(yp) free(yp);
    if(debugInfo) free(debugInfo);
    return xp;
}



// --------------------------------------------------------------------------
// Fij
/*
MAGNET_FUNCTION7(   double, Fij,
                    double u,
                    double v,
                    double betai,
                    double betaj,
                    double qM,
                    double qZ,
                    long nbPoints)
{
    return Fij(u,v,betai,betaj,qM,qZ,nbPoints);
}
*/

// --------------------------------------------------------------------------
// F2
MAGNET_X_FUNCTION9( double, F2,
                    double u1, MAGNET_MANDATORY,
                    double u2, MAGNET_MANDATORY,
                    double beta1, MAGNET_MANDATORY,
                    double beta2, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ1, MAGNET_MANDATORY,
                    double qZ2, MAGNET_MANDATORY,
                    long nbPoints, MAGNET_MANDATORY,
                    double eps, = 1e-6)
{
    return F2(u1,u2,beta1,beta2,qM,qZ1,qZ2,eps,nbPoints);
}

// --------------------------------------------------------------------------
// Finv
//
MAGNET_X_FUNCTION6(   double, Finv,
                    double p, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double qM, MAGNET_MANDATORY,
                    double qZ, MAGNET_MANDATORY,
                    long nbPoints, = 1000,
                    double eps, = 1e-6)
{
    return Finv(p,beta,qM,qZ,nbPoints,eps);
}

// --------------------------------------------------------------------------
// Integrand
//
MAGNET_FUNCTION5(   double, Integrand,
                    double u,
                    double M,
                    double beta,
                    double qM,
                    double qZ)
{
    return Integrand(u,M,beta,qM,qZ);
}

// --------------------------------------------------------------------------
// QCum
//
MAGNET_FUNCTION2(double, QCum,
                 double x,
                 double q)
{
    return NormalCum(NormMapinv(x,q));
}

// --------------------------------------------------------------------------
// QCuminv
//
MAGNET_FUNCTION2(double, QCuminv,
                 double y,
                 double q)
{
    return NormMap(NormalCumInverse(y),q);
}

// --------------------------------------------------------------------------
// QDensity
//
MAGNET_FUNCTION2(double, QDensity,
                 double x,
                 double q)
{
    return NormalDensity(NormMapinv(x,q))*DNormMapinv(x,q);
}

// --------------------------------------------------------------------------
// Fdelta
//
MAGNET_FUNCTION5(double, Fdelta,
                 double x,
                 double a,
                 double b,
                 double c,
                 double delta)
{
    return Fdelta(x,a,b,c,delta);
}

// --------------------------------------------------------------------------
// FdeltaC
//
MAGNET_FUNCTION4(double, FdeltaC,
                 double x,
                 double a,
                 double b,
                 double delta)
{
    return FdeltaC(x,a,b,delta);
}

// --------------------------------------------------------------------------
// Fdeltainv
//
MAGNET_FUNCTION5(double, Fdeltainv,
                 double y,
                 double a,
                 double b,
                 double c,
                 double delta)
{
    return Fdeltainv(y,a,b,c,delta);
}

// --------------------------------------------------------------------------
// FdeltaCinv
//
MAGNET_FUNCTION4(double, FdeltaCinv,
                 double y,
                 double a,
                 double b,
                 double delta)
{
    return FdeltaCinv(y,a,b,delta);
}


// --------------------------------------------------------------------------
// DeltaMap
//
MAGNET_FUNCTION5(double, DeltaMap,
                 double x,
                 double a,
                 double b,
                 double c,
                 double delta)
{
    return DeltaMap(x, a, b, c, delta);
}
// --------------------------------------------------------------------------
MAGNET_FUNCTION3(double, q, double x,
                            double y,
                            double alpha)
{
    return alpha*alpha/(1.-alpha*alpha)*(x*x-2./alpha*x*y+y*y);
}

// --------------------------------------------------------------------------
MAGNET_FUNCTION1(Array<double>, Integer, long n)
{
    Array<double> out(n);
    Integer(n, &out[0]);
    return out;
}

MAGNET_FUNCTION1(Matrix<double>, MatInteger, long n)
{
    Matrix<double> out(n,n);
    MatInteger(n, &out[0][0]);
    return out;
}

// --------------------------------------------------------------------------
MAGNET_FUNCTION3(Array<double>, calcLossDensity, 
                const Copula &copula,
                const Array<int> &lossAmt,
                const Array<double> &survProbas
                )
{
    size_t maxLoss = 0U;
    for (int i=0; i<lossAmt.size(); ++i)
    {
        if (lossAmt[i] < 0) CM_THROW RuntimeError("all loss amount must be positive");
        maxLoss += lossAmt[i];
    }
    PRODUCT_PARAM *pp = (PRODUCT_PARAM *) (&copula._copula)->param;
    Array<double> density(maxLoss+1U);
    calcLossDensity(
        &lossAmt[0],    
        lossAmt.size(),
        maxLoss,
        &copula._copula, 
        &survProbas[0],         /* (I) array[n] of uncond survival proba */
        &density[0]             /* (O) array[maxLoss+1] of loss density  */
        );
    DR_ERROR_CHECK
    return density;
}

// --------------------------------------------------------------------------
MAGNET_FUNCTION5(double, TrancheletLoss3C,
                        double u,
                        double beta,
                        double a,
                        double c,
                        double K)
{
    return TrancheletLoss3C(u, beta, a, c, K);
}

MAGNET_FUNCTION6(double, TrancheLoss3C,
                        double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1)
{
    return TrancheLoss3C(u, beta, a, c, K0,K1);
}

MAGNET_X_FUNCTION6(double, DTrancheletLoss3CDbeta,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheletLoss3CDbeta(u,beta,a,c,K,h);
}

MAGNET_X_FUNCTION6(double, DTrancheletLoss3CDa,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheletLoss3CDa(u,beta,a,c,K,h);
}

MAGNET_X_FUNCTION6(double, DTrancheletLoss3CDc,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheletLoss3CDc(u,beta,a,c,K,h);
}

MAGNET_X_FUNCTION7(double, DTrancheletLoss3CDcAAdjust,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K, MAGNET_MANDATORY,
                        double Kmatch, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheletLoss3CDcAAdjust(u,beta,a,c,K,Kmatch,h);
}

MAGNET_X_FUNCTION7(double, DTrancheLoss3CDbeta,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K0, MAGNET_MANDATORY,
                        double K1, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheLoss3CDbeta(u,beta,a,c,K0,K1,h);
}

MAGNET_X_FUNCTION7(double, DTrancheLoss3CDa,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K0, MAGNET_MANDATORY,
                        double K1, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheLoss3CDa(u,beta,a,c,K0,K1,h);
}

MAGNET_X_FUNCTION7(double, DTrancheLoss3CDc,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K0, MAGNET_MANDATORY,
                        double K1, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheLoss3CDc(u,beta,a,c,K0,K1,h);
}

MAGNET_X_FUNCTION8(double, DTrancheLoss3CDcAAdjust,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K0, MAGNET_MANDATORY,
                        double K1, MAGNET_MANDATORY,
                        double Kmatch, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheLoss3CDcAAdjust(u,beta,a,c,K0,K1,Kmatch,h);
}

MAGNET_X_FUNCTION8(double, DTrancheletLoss3CDcATAdjust,
                        double u, MAGNET_MANDATORY,
                        double beta, MAGNET_MANDATORY,
                        double a, MAGNET_MANDATORY,
                        double c, MAGNET_MANDATORY,
                        double K, MAGNET_MANDATORY,
                        double Kmatch0, MAGNET_MANDATORY,
                        double Kmatch1, MAGNET_MANDATORY,
                        double h, = 0.01)
{
    return DTrancheletLoss3CDcATAdjust(u,beta,a,c,K,Kmatch0, Kmatch1,h);
}

MAGNET_FUNCTION5(double, Aeq,
                double u,
                double beta_star,
                double beta,
                double c,
                double K)
{
    return Aeq(u,beta_star,beta,c,K);
}


MAGNET_FUNCTION6(double, ATeq,
                double u,
                double beta_star,
                double beta,
                double c,
                double K0,
                double K1)
{
    return ATeq(u,beta_star,beta,c,K0,K1);
}


MAGNET_FUNCTION5(   double, LossDensity3C,
                    double u,
                    double beta,
                    double a,
                    double c,
                    double K)
{
    return LossDensity3C(u,beta,a,c,K);
}


MAGNET_X_FUNCTION7( double, Moment3C,
                    double u, MAGNET_MANDATORY,
                    double beta, MAGNET_MANDATORY,
                    double a, MAGNET_MANDATORY,
                    double c, MAGNET_MANDATORY,
                    long n, = 2,
                    long nbPoint, = 1000,
                    double eps, = 1e-8)
{
    return Moment3C(u,beta,a,c,n,nbPoint,eps);
}

MAGNET_X_FUNCTION6( double, BetaEq,
                    double u, MAGNET_MANDATORY,
                    double beta_star, MAGNET_MANDATORY,
                    double a, MAGNET_MANDATORY,
                    double c, MAGNET_MANDATORY,
                    long nbPoint, =1000,
                    double eps, = 1e-8)
{
    double x = BetaEq(u,beta_star,a,c,nbPoint,eps);
    DR_ERROR_CHECK
    return x;
}

MAGNET_FUNCTION2(   Matrix <double>, GaussianCopulatedUniformDeviates,
                    long nbPaths,
                    const Array<double> &beta)
{
    long nbNames = beta.size();
    Matrix<double> deviates(nbPaths,nbNames);
    GaussianCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths,
                                      &beta[0]);
    return deviates;
}

MAGNET_X_FUNCTION3( Matrix <double>, StudentCopulatedUniformDeviates,
                    long nbPaths, MAGNET_MANDATORY,
                    const Array<double> &beta, MAGNET_MANDATORY,
                    long freedomDegree, = 10)
{
    long nbNames = beta.size();
    Matrix<double> deviates(nbPaths,nbNames);
    StudentCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths,
                                      &beta[0],
                                      freedomDegree);
    return deviates;
}

MAGNET_FUNCTION2(   Matrix <double>, DependenceCopulatedUniformDeviates,
                    long nbPaths,
                    long nbNames)
{
    Matrix<double> deviates(nbPaths,nbNames);
    DependenceCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths);
    return deviates;
}

MAGNET_FUNCTION2(   Matrix <double>, IndependenceCopulatedUniformDeviates,
                    long nbPaths,
                    long nbNames)
{
    Matrix<double> deviates(nbPaths,nbNames);
    IndependenceCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths);
    return deviates;
}


MAGNET_X_FUNCTION3( Matrix <double>, ClaytonCopulatedUniformDeviates,
                    long nbPaths, MAGNET_MANDATORY,
                    long nbNames, MAGNET_MANDATORY,
                    double alpha, = 2.0)
{
    Matrix<double> deviates(nbPaths,nbNames);
    ClaytonCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths,
                                      alpha);
    return deviates;
}

MAGNET_X_FUNCTION3( Matrix <double>, GumbelCopulatedUniformDeviates,
                    long nbPaths, MAGNET_MANDATORY,
                    long nbNames, MAGNET_MANDATORY,
                    double theta, = 2.0)
{
    Matrix<double> deviates(nbPaths,nbNames);
    GumbelCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths,
                                      theta);
    return deviates;
}

MAGNET_X_FUNCTION4( Matrix <double>, QCopulatedUniformDeviates,
                    long nbPaths, MAGNET_MANDATORY,
                    const Array<double> &beta, MAGNET_MANDATORY,
                    double qM, = 0.0,
                    double qZ, = 0.0)
{
    long nbNames = beta.size();
    Matrix<double> deviates(nbPaths,nbNames);
    QCopulatedUniformDeviates(&deviates[0][0],
                                      nbNames,
                                      nbPaths,
                                      &beta[0],
                                      qM,
                                      qZ);
    return deviates;
}
