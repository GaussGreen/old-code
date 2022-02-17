#include <General/General.h>
#include <Magnet/Magnet.h>
extern "C" {
#include "../student.h"
#include "../gaussian.h"
#include "../dependence.h"
#include "../independence.h"
#include "../proba_utils.h"
#include "../payoff_r.h"
#include "../alpha_stable_random.h"
#include "../gamma_random.h"
}

using namespace CM;
MAGNET_PREFIX("DR3_")
MAGNET_CATEGORY("Copula Utils r")
// --------------------------------------------------------------------------
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
struct IndicatorRecoverySim : public CM::Object
{
    Matrix<int> indicator; //[nbPath][nbNames]
    Matrix<double> X;
    Array<double> weight;
};

inline INDICATOR_RECOVERY_SIM IndicatorRecoverySim2Struct(const IndicatorRecoverySim &is)
{
    INDICATOR_RECOVERY_SIM i_s;
    i_s.indicator = const_cast<int *>(&is.indicator[0][0]);
    i_s.X = const_cast<double *>(&is.X[0][0]);
    i_s.nbNames = is.indicator.colSize();
    i_s.nbPaths = is.indicator.rowSize();
    i_s.weight = const_cast<double *>(&is.weight[0]);
    return i_s;
}

inline IndicatorRecoverySim Struct2IndicatorRecoverySim(const INDICATOR_RECOVERY_SIM &i_s)
{
    IndicatorRecoverySim is;
    ConstMatrixView<int> m(i_s.indicator,i_s.nbPaths,i_s.nbNames);
    is.indicator = m;
    ConstMatrixView<double> m2(i_s.X,i_s.nbPaths,i_s.nbNames);
    is.X = m2;
    is.weight = Array<double>(i_s.weight,i_s.weight+i_s.nbPaths);
    return is;
}

/** define a credit portfolio object that can be mapped into C struct */
struct CreditPortfolioR : public CM::Object
{
    Array<double> notional;
    Array<double> recovery;
    Array<double> sigmaRecovery;
    Array<double> nameMaturity;
    double strike1;
    double strike2;
};

inline CREDIT_PORTFOLIO_R CreditPortfolioR2Struct(const CreditPortfolioR &cp)
{
    CREDIT_PORTFOLIO_R c_p;
    c_p.nbNames = cp.notional.size();
    c_p.notional = const_cast<double *>(&cp.notional[0]);
    c_p.recovery = const_cast<double *>(&cp.recovery[0]);
    c_p.sigma_recovery = const_cast<double *>(&cp.sigmaRecovery[0]);
    c_p.nameMaturity = const_cast<double *>(&cp.nameMaturity[0]);
    c_p.strike1 = cp.strike1;
    c_p.strike2 = cp.strike2;
    return c_p;
}

inline CreditPortfolioR Struct2CreditPortfolioR(CREDIT_PORTFOLIO_R &c_p)
{
    CreditPortfolioR cp;
    cp.notional = Array<double>(c_p.notional, c_p.notional+c_p.nbNames);
    cp.recovery = Array<double>(c_p.recovery, c_p.recovery+c_p.nbNames);
    cp.sigmaRecovery = Array<double>(c_p.sigma_recovery, c_p.sigma_recovery+c_p.nbNames);
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
class CopulaR : public CM::Object
{
public:
    COPULA_R _copulaR; 
    // ctor and dtor memory mgt of the C struct.
    CopulaR(const COPULA_R &c):_copulaR(c) {}
    ~CopulaR() {CopulaRFree(&_copulaR);}
private:
    // disable copy ctor and assgn operator
    CopulaR operator=(const CopulaR &a);
    CopulaR(const CopulaR &a);
};

// --------------------------------------------------------------------------
// MAGNET code
//
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// CreditPortfolioRCreate
//
MAGNET_DESCRIPTION("returns a credit portfolio R object")
MAGNET_PARAMETER("",notional,"","Array of the names' notionals")
MAGNET_PARAMETER("",recovery,"","Array of the names' recoveries")
MAGNET_PARAMETER("",nameMaturity,"","Array of the names' maturities")
MAGNET_PARAMETER("",strike1,"","lower strike")
MAGNET_PARAMETER("",strike2,"","upper strike")
MAGNET_FUNCTION6(SharedPointer<CreditPortfolioR>, CreditPortfolioRCreate,
                   const Array<double> &notional, 
                   const Array<double> &recovery, 
                   const Array<double> &sigmaRecovery,
                   const Array<double> &nameMaturity, 
                   double  strike1, 
                   double  strike2)
{
    SharedPointer<CreditPortfolioR> cp = Raw2SharedPointer(
        new CreditPortfolioR());
    cp->notional = notional;
    cp->recovery = recovery;
    cp->sigmaRecovery = sigmaRecovery;
    cp->nameMaturity = nameMaturity;
    cp->strike1 = strike1;
    cp->strike2 = strike2;
    return cp;
}

// --------------------------------------------------------------------------
// CopulaCreate
//
MAGNET_DESCRIPTION("returns a copula R object")
MAGNET_PARAMETER("",type,"","one of Gaussian, Student, Dependence, Independence, Gumble, Clayton")
MAGNET_PARAMETER("",beta,"","vector of weights (or alpha for gumbel)")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom ofthe student variables")
MAGNET_X_FUNCTION4(SharedPointer<CopulaR>,
                 CopulaRCreate,
                 String type, MAGNET_MANDATORY,
                 long nbNames, MAGNET_MANDATORY,
                 const Array<double> &beta, =Array<double>(),
                 const Array<double> &alpha, = Array<double>())
{
    // map types
    SymbolMap<int> theType(
        0, "gaussian",
        1, "student",
        2, "dependence",
        3, "independence",
        4, "gumbel",
        5, "clayton");
  	int pos = -1;
    int status = FAILURE;
	theType.name2Value(type, pos);
    if (pos < 0 ) 
        CM_THROW RuntimeError("type must be either gaussian, student, dependence, independence, gumbel, clayton");

    COPULA_R c = CopulaRCreate(pos,nbNames,
        const_cast<double *>(&beta[0]), const_cast<double *>(&alpha[0]), &status);
    // shallow copy
    return Raw2SharedPointer(new CopulaR(c));
}


// --------------------------------------------------------------------------
// IndicatorRecoverySimCreate
//
MAGNET_DESCRIPTION("returns an indicator recovery object")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_PARAMETER("",survivalProba,"","array[nbNames] of survival probabilities")
MAGNET_PARAMETER("",copula,"","copula object created by CopulaCreate")
MAGNET_X_FUNCTION4(SharedPointer<IndicatorRecoverySim> , IndicatorRecoverySimCreate,
                 long nbPaths, MAGNET_MANDATORY ,
                 const Array<double> &survivalProba, MAGNET_MANDATORY ,
                 const CopulaR &copulaR, MAGNET_MANDATORY ,
                 const Array<double> &Msample , = Array<double> () )
{
    long nbNames = survivalProba.size();
    INDICATOR_RECOVERY_SIM i_s = 
        IndicatorRecoverySimCreate( nbPaths,
                            nbNames,
                            &survivalProba[0],
                            &copulaR._copulaR,
                            Msample.size() ? &Msample[0] : (const double *)0,
                            Msample.size());

    SharedPointer<IndicatorRecoverySim> is = Raw2SharedPointer(new IndicatorRecoverySim());
    // deep copy
    *is = Struct2IndicatorRecoverySim(i_s);
    IndicatorRecoverySimFree(&i_s);
    return is;
}


// --------------------------------------------------------------------------
// ExpectedPayoff
//
MAGNET_DESCRIPTION("returns the expected payoff of a credit portfolio R")
MAGNET_PARAMETER("",port,"","credit portfolio R structure created by CreditPortfolioCreate")
MAGNET_PARAMETER("",indicatorRecoverySim,"","Indicator recovery simulation structure created by IndicatorSimCreate")
MAGNET_FUNCTION2(double, ExpectedPayoff,
                 const CreditPortfolioR &creditPortfolioR,
                 const IndicatorRecoverySim &indicatorRecoverySim)
{
    double res;
    CREDIT_PORTFOLIO_R c_r = CreditPortfolioR2Struct(creditPortfolioR);
    INDICATOR_RECOVERY_SIM i_s = IndicatorRecoverySim2Struct(indicatorRecoverySim);
    res = ExpectedPayoff_R(&c_r, &i_s);
    return res;
}

// --------------------------------------------------------------------------
// PayoffDistribution
//
MAGNET_FUNCTION3(Array<double>, PayoffDistribution,
                            const CreditPortfolioR &creditPortfolioR,
                            const IndicatorRecoverySim &indicatorRecoverySim,
                            const Array<double> &sampleLoss)
{
    long nbSampleLoss = sampleLoss.size();
    int status = FAILURE;
    CREDIT_PORTFOLIO_R c_r = CreditPortfolioR2Struct(creditPortfolioR);
    INDICATOR_RECOVERY_SIM i_s = IndicatorRecoverySim2Struct(indicatorRecoverySim);
    Array<double> loss = Array<double>(nbSampleLoss);
    status = PayoffDistribution_R( &c_r,
                                    &i_s,
                                    &sampleLoss[0],
                                    nbSampleLoss,
                                    &loss[0]);
    return loss;
}

// --------------------------------------------------------------------------
// NbDefaultNameDistribution
//
MAGNET_DESCRIPTION("returns an array[0..nbNames] of discrete probability of k names defaulted")
MAGNET_PARAMETER("",indicatorRecoverySim,"","IndicatorRecoverySim object create by IndicatorSimCreate")
MAGNET_FUNCTION1(Array<double>, NbDefaultNameDistribution,
                 const IndicatorRecoverySim &indicatorRecoverySim)
{
    INDICATOR_RECOVERY_SIM i_s = IndicatorRecoverySim2Struct(indicatorRecoverySim);
    long nbNames = i_s.nbNames;
    Array<double> out = Array<double>(nbNames+1);
    NbDefaultNameDistribution_R( i_s, &out[0]);
    return out;
}
