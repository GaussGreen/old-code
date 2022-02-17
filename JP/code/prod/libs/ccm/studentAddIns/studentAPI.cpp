#include <General/General.h>
#include <Magnet/Magnet.h>
extern "C" {
#include "../student.h"
#include "../gaussian.h"
#include "../dependence.h"
#include "../independence.h"
#include "../proba_utils.h"
#include "../payoff.h"
}

using namespace CM;
MAGNET_PREFIX("DR_")
MAGNET_CATEGORY("Copula Utils")

// --------------------------------------------------------------------------
// DefaultTimeFromUniform
//
MAGNET_DESCRIPTION("Returns a default time given a spread and a trigger level")
MAGNET_PARAMETER("",time,"","vector of times in 'double' format")
MAGNET_PARAMETER("",spread,"","vector of the constant spreads, ith elements is the spread between ith and (i+1)th time")
MAGNET_PARAMETER("",triggerLevel,"","default happens when the 'credit discount factor' hits the trigger level")
MAGNET_FUNCTION3(   double, DefaultTimeFromUniform,
                    const Array<double> &time,
                    const Array<double> &spread,
                    double triggerLevel)
{
    return DefaultTimeFromUniform(  &time[0],
                                    &spread[0],
                                    time.size(),
                                    triggerLevel);
}

// --------------------------------------------------------------------------
// StudentCum
//
MAGNET_DESCRIPTION("Returns the cumulative probability of a Student variable")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom")
MAGNET_FUNCTION2(   double, StudentCum,
                    double x,
                    long freedomDegree)
{
    return StudentCum(  x,
                        freedomDegree);
}

// --------------------------------------------------------------------------
// StudentCumInv
//
MAGNET_DESCRIPTION("Returns the invert cumulative probability of a Student variable")
MAGNET_PARAMETER("",y,"","")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom")
MAGNET_FUNCTION2(   double, StudentCumInv,
                    double y,
                    long freedomDegree)
{
    return StudentCumInverse(y, freedomDegree);
}

// --------------------------------------------------------------------------
// StableCum
//
MAGNET_DESCRIPTION("Returns the cumulative probability of a Stable variable")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",alpha,"","alpha")
MAGNET_PARAMETER("",beta,"","beta")
MAGNET_FUNCTION3(   double, StableCum,
                    double x,
                    double alpha, double beta)
{
    return StableCum(x,alpha,beta);
}

// --------------------------------------------------------------------------
// StableDensity
//
MAGNET_DESCRIPTION("Returns the density of a Stable variable")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",alpha,"","alpha")
MAGNET_PARAMETER("",beta,"","beta")
MAGNET_FUNCTION3(   double, StableDensity,
                    double x,
                    double alpha, double beta)
{
    return StableDensity(x,alpha,beta);
}

// --------------------------------------------------------------------------
// StableCumInv
//
MAGNET_DESCRIPTION("Returns the invert cumulative probability of a Stable variable")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",alpha,"","alpha")
MAGNET_PARAMETER("",beta,"","beta")
MAGNET_FUNCTION3(   double, StableCumInv,
                    double x,
                    double alpha, double beta)
{
    return StableCumInverse(  x,alpha,beta);
}

// --------------------------------------------------------------------------
// Chi2Density
//
MAGNET_DESCRIPTION("returns the chi2 distribution")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom")
MAGNET_FUNCTION2(   double, Chi2Density,
                    double x,
                    long freedomDegree)
{
    return Chi2Density( x,
                        freedomDegree);
}

// --------------------------------------------------------------------------
// Chi2Deviates
//
MAGNET_DESCRIPTION("Create an array of chi2 deviates")
MAGNET_PARAMETER("",nbPaths,"","number of paths in the simulation, must be an integer >0")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom of the chi2 deviates, must be an integer >0")
MAGNET_FUNCTION2(Array<double>, Chi2Deviates,    
    long nbPaths,
    long freedomDegree)
{
    Array<double> out = Array<double>(nbPaths);
    Chi2DeviatesSobol(   &out[0],
                    nbPaths,
                    freedomDegree);
    return out;
}

// --------------------------------------------------------------------------
// GaussianDeviates
//
MAGNET_DESCRIPTION("Create a matrix of student deviates")
MAGNET_PARAMETER("",nbPaths,"","number of paths in the simulation, must be an integer >0")
MAGNET_PARAMETER("",beta,"","vector of weights")
MAGNET_FUNCTION2(Matrix<double>, GaussianMultivariateDeviates, 
    long nbPaths,
    const Array<double> &beta)
{
    long nbNames = beta.size();
    Matrix<double> out = Matrix<double>(nbPaths,nbNames);
    GaussianMultivariateDeviates(   &out[0][0],
                                    nbNames,
                                    nbPaths,
                                    &beta[0]);
    return out;
}

// --------------------------------------------------------------------------
// GaussianCopulatedUniformDeviates
// 
MAGNET_DESCRIPTION("Returns a matrix of Gaussian copulated uniform deviates")
MAGNET_PARAMETER("",nbPaths,"","number of paths in the simulation, must be an integer >0")
MAGNET_PARAMETER("",beta,"","vector of weights")
MAGNET_FUNCTION2(Matrix<double>, GaussianCopulatedUniformDeviates,
                long nbPaths,
                const Array<double> &beta)
{
    long nbNames = beta.size();
    Matrix<double>out = Matrix<double>(nbPaths,nbNames);
    GaussianCopulatedUniformDeviates(   &out[0][0],
                                        nbNames,
                                        nbPaths,
                                        &beta[0]);
    return out;
}

// --------------------------------------------------------------------------
// StudentDeviates
//
MAGNET_DESCRIPTION("Create a matrix of student deviates")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_PARAMETER("",beta,"","vector of weights")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom ofthe student variables")
MAGNET_FUNCTION3(Matrix<double>, StudentMultivariateDeviates, 
    long nbPaths,
    const Array<double> &beta,
    long freedomDegree)
{
    long nbNames = beta.size();
    Matrix<double> out = Matrix<double>(nbPaths,nbNames);
    StudentMultivariateDeviates(&out[0][0],
                                nbNames,
                                nbPaths,
                                &beta[0],
                                freedomDegree);
    return out;
}

// --------------------------------------------------------------------------
// StudentCopulatedUniformDeviates
// 
MAGNET_DESCRIPTION("Returns a matrix of student copulated uniform deviates")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_PARAMETER("",beta,"","vector of weights")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom of the student variables")
MAGNET_FUNCTION3(Matrix<double>, StudentCopulatedUniformDeviates,
                long nbPaths,
                const Array<double> &beta, 
                long freedomDegree)
{
    long nbNames = beta.size();
    Matrix<double>out = Matrix<double>(nbPaths,nbNames);
    StudentCopulatedUniformDeviates(&out[0][0],
                                    nbNames,
                                    nbPaths,
                                    &beta[0],
                                    freedomDegree);
    return out;
}

// --------------------------------------------------------------------------
// DensityBarChart
//
MAGNET_DESCRIPTION("returns a vector of the 'steped density'")
MAGNET_PARAMETER("",sample,"","vector of random numbers")
MAGNET_PARAMETER("",start,"","start valueof the chart")
MAGNET_PARAMETER("",stepSize,"","size of a bar in the chart")
MAGNET_PARAMETER("",nbSteps,"","number of step in the bar chart")
MAGNET_FUNCTION4(   Array<double>, DensityBarChart,
                    const Array<double> &sample,
                    double start, 
                    double stepSize,
                    long nbSteps)
{
    long sampleSize = sample.size();
    Array<double>out = Array<double>(nbSteps);
    DensityBarChart(&out[0],
                    &sample[0],
                    sampleSize,
                    start,
                    stepSize,
                    nbSteps);
    return out;
}

// --------------------------------------------------------------------------
// GaussianCopulaLinearCorrelation
//
MAGNET_DESCRIPTION("returns the linear correlation given a Gaussian copula generating correlation coefficient")
MAGNET_PARAMETER("",beta,"","correlation used to generate the copula")
MAGNET_FUNCTION1(   double, GaussianCopulaLinearCorrelation,
                    double beta)
{
    return GaussianCopulaLinearCorrelation(beta);
}


// --------------------------------------------------------------------------
// StudentCopulaLinearCorrelation
//
MAGNET_DESCRIPTION("returns the linear correlation given a Student copula generating correlation coefficient")
MAGNET_PARAMETER("",beta,"","correlation used to generate the copula")
MAGNET_X_FUNCTION2(double, StudentCopulaLinearCorrelation,
                   double beta, MAGNET_MANDATORY,
                   long freedomDegree, MAGNET_MANDATORY)
{
    return StudentCopulaLinearCorrelation(  beta,
                                            freedomDegree);
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// Magnet to C adaptor code
//

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

/** define a spread object that can be mapped into C struct */
struct SpreadCurve : public CM::Object
{
    Array<double> time;
    Array<double> spread;
};

inline SPREAD_CURVE SpreadCurve2Struct(const SpreadCurve &sc)
{
    SPREAD_CURVE s_c;
    s_c.nbTimes = sc.spread.size();
    s_c.time = const_cast<double *> (&sc.time[0]);
    s_c.spread = const_cast<double *>(&sc.spread[0]);
    return s_c;
}

inline SpreadCurve Struct2SpreadCurve(SPREAD_CURVE &s_c)
{
    SpreadCurve sc;
    sc.time = Array<double>(s_c.time,s_c.time+s_c.nbTimes);
    sc.spread = Array<double>(s_c.spread,s_c.spread+s_c.nbTimes);
    return sc;
}

/** define a uniform sim object that can be mapped into C struct */
struct UniformSim : public CM::Object
{
    Matrix<double> uniformDeviate;
    Array<double> weight;
};

inline UNIFORM_SIM UniformSim2Struct(const UniformSim &us)
{
    UNIFORM_SIM u_s;
    u_s.nbNames = us.uniformDeviate.colSize();
    u_s.nbPaths = us.uniformDeviate.rowSize();
    u_s.uniformDeviate = const_cast<double *>(&us.uniformDeviate[0][0]);
    u_s.weight = const_cast<double *>(&us.weight[0]);
    return u_s;
}

inline UniformSim Struct2UniformSim(const UNIFORM_SIM &u_s)
{
    UniformSim us;
    ConstMatrixView<double> m(u_s.uniformDeviate,u_s.nbPaths, u_s.nbNames);
    us.uniformDeviate = m;
    us.weight = Array<double>(u_s.nbPaths);
    us.weight.fill(*u_s.weight);
    return us;
}

/** define a scenario sim object that can be mapped into C struct */
struct ScenarioSim : public CM::Object
{
    Matrix<double> defaultTime; //[nbPath][nbNames]
    Array<double> weight;
};

inline SCENARIO_SIM ScenarioSim2Struct(const ScenarioSim &s)
{
    SCENARIO_SIM s_s;
    s_s.defaultTime = const_cast<double *>(&s.defaultTime[0][0]);
    s_s.nbNames = s.defaultTime.colSize();
    s_s.nbPaths = s.defaultTime.rowSize();
    s_s.weight = const_cast<double *>(&s.weight[0]);
    return s_s;
}

inline ScenarioSim Struct2ScenarioSim(const SCENARIO_SIM &s_s)
{
    ScenarioSim ss;
    ConstMatrixView<double> m(s_s.defaultTime,s_s.nbPaths,s_s.nbNames);
    ss.defaultTime = m;
    ss.weight = Array<double>(s_s.weight,s_s.weight+s_s.nbPaths);
    return ss;
}
// --------------------------------------------------------------------------
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
// SpreadCurveCreate
//
MAGNET_DESCRIPTION("returns a spread curve object")
MAGNET_PARAMETER("",time,"","vector of times in 'double' format")
MAGNET_PARAMETER("",spread,"","vector of the constant spreads, ith elements is the spread between ith and (i+1)th time")
MAGNET_FUNCTION2(SharedPointer<SpreadCurve>, SpreadCurveCreate,
                const Array<double> &time,
                const Array<double> &spread)
{
    SharedPointer<SpreadCurve> sc = Raw2SharedPointer(
        new SpreadCurve());
    sc->spread = spread;
    sc->time = time;
    return sc;
}

// --------------------------------------------------------------------------
// UniformSimCreate
//
MAGNET_DESCRIPTION("returns a simulation object")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_PARAMETER("",copulaType,"","0 is Gaussian, 1 is Student, 2 is Dependence")
MAGNET_PARAMETER("",beta,"","vector of weights")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom ofthe student variables")
MAGNET_FUNCTION4(SharedPointer<UniformSim>, UniformSimCreate,
                 long nbPaths,
                 int copulaType,
                 const Array<double> &beta,
                 long freedomDegree)
{
    long nbNames = beta.size();
    SharedPointer<UniformSim> us = Raw2SharedPointer(new UniformSim());
    UNIFORM_SIM u_s = 
        UniformSimCreate(   nbPaths,
                            nbNames,
                            copulaType,
                            &beta[0],
                            freedomDegree);
    *us = Struct2UniformSim(u_s);
    UniformSimFree(&u_s);
    return us;
}

// --------------------------------------------------------------------------
// CopulaProduct
//
MAGNET_DESCRIPTION("Returns a matrix of product copulated uniform deviates")
MAGNET_FUNCTION3(SharedPointer<UniformSim>,
                 CopulaProduct,
                 const UniformSim &uniformSim1,
                 const UniformSim &uniformSim2,
                 const Matrix<double> &alpha)
{
    SharedPointer<UniformSim> us = Raw2SharedPointer(
        new UniformSim());
    UNIFORM_SIM u_s1 = UniformSim2Struct(uniformSim1);
    UNIFORM_SIM u_s2 = UniformSim2Struct(uniformSim2);
    UNIFORM_SIM u_s;
    u_s = CopulaProduct(u_s1, u_s2, &alpha[0][0]);
    *us = Struct2UniformSim(u_s);
    return us;
}


// --------------------------------------------------------------------------
// ScenarioSimCreate
//
MAGNET_DESCRIPTION("returns a scenario object")
MAGNET_PARAMETER("",us, "","simulation object")
MAGNET_PARAMETER("",sc, "","array of spread curve objects")
MAGNET_FUNCTION2(SharedPointer<ScenarioSim>, ScenarioSimCreate,
                const UniformSim &us,
                const Array<SharedPointer<SpreadCurve> > &sc)
{
    SharedPointer<ScenarioSim> ss = Raw2SharedPointer(
        new ScenarioSim());
    UNIFORM_SIM u_s = UniformSim2Struct(us);
    SPREAD_CURVE *s_c = new SPREAD_CURVE[sc.size()];
    long i = 0;
    for(i=0;i<sc.size();i++)
    {
        s_c[i] = SpreadCurve2Struct(*sc[i]);
    }
    SCENARIO_SIM s_s = ScenarioSimCreate(&u_s, s_c);
    *ss = Struct2ScenarioSim(s_s);
    ScenarioSimFree(&s_s);
    return ss;
}

// --------------------------------------------------------------------------
// ScenarioSimCreate2
//
MAGNET_DESCRIPTION("returns a scenario object")
MAGNET_PARAMETER("",uniformSim1, "","simulation object created by UniformSimCreate")
MAGNET_PARAMETER("",uniformSim2, "","simulation object created by UniformSimCreate")
MAGNET_PARAMETER("",spreadCurve, "","array of spread curve objects created by SpreadCurveCreate")
MAGNET_FUNCTION4(SharedPointer<ScenarioSim>, ScenarioSimCreate2,
                const UniformSim &uniformSim1,
                const UniformSim &uniformSim2,
                const Array<SharedPointer<SpreadCurve> > &spreadCurve,
                const Array<double> &alpha)
{
    SharedPointer<ScenarioSim> ss = Raw2SharedPointer(
        new ScenarioSim());
    long nbNames = spreadCurve.size();
    long i = 0;
    long k = 0;
    UNIFORM_SIM u_s1;
    UNIFORM_SIM u_s2;

    SPREAD_CURVE *s_c = new SPREAD_CURVE[nbNames];
    for(i=0;i<nbNames;i++)
    {
        s_c[i] = SpreadCurve2Struct(*spreadCurve[i]);
    }

    u_s1 = UniformSim2Struct(uniformSim1);
    u_s2 = UniformSim2Struct(uniformSim2);

    SCENARIO_SIM s_s = ScenarioSimCreate2(&u_s1, &u_s2, s_c, &alpha[0]);
    *ss = Struct2ScenarioSim(s_s);
    ScenarioSimFree(&s_s);
    return ss;
}

// --------------------------------------------------------------------------------
// ExpectedLoss
//
MAGNET_DESCRIPTION("returns a portfolio expected loss against time")
MAGNET_PARAMETER("",creditPortfolio,"","credit portfolio object created by CreditPortfolioCreate")
MAGNET_PARAMETER("",scenarioSim,"","scenario object created by ScenarioSimCreate")
MAGNET_PARAMETER("",samplingTime,"", "array of sampling times")
MAGNET_FUNCTION3(Array<double>, ExpectedLoss,
                   const CreditPortfolio &creditPortfolio,
                   const ScenarioSim     &scenarioSim,
                   const Array<double>   &samplingTime)
{
    CREDIT_PORTFOLIO c_p = CreditPortfolio2Struct(creditPortfolio);
    SCENARIO_SIM     s_s = ScenarioSim2Struct(scenarioSim);
    double *res = ExpectedPayoff(
        &c_p,
        &s_s,
        &samplingTime[0], 
        samplingTime.size());
    return Array<double>(res,res+samplingTime.size());

}


// ---------------------------------------------------------------------------
// CumulatedSpreadToFlatForward
//
MAGNET_FUNCTION2(   Array<double>, CumulatedSpreadToFlatForward,
                    Array<double> &time,
                    Array<double> &cumSpread)
{
    double *res = CumulatedSpreadToFlatForward(&time[0], &cumSpread[0], time.size());
    return Array<double>(res,res+time.size());
}

// --------------------------------------------------------------------------
// DependenceCopulatedUniformDeviates
// 
MAGNET_DESCRIPTION("Returns a matrix of dependence copulated uniform deviates")
MAGNET_PARAMETER("",nbNames,"nbNames must be an integer >0","number of names in the simulation")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_FUNCTION2(Matrix<double>, DependenceCopulatedUniformDeviates,
                 long nbNames,
                 long nbPaths)
{
    Matrix<double>out = Matrix<double>(nbPaths,nbNames);
    DependenceCopulatedUniformDeviates(&out[0][0],
                                    nbNames,
                                    nbPaths);
    return out;
}

// --------------------------------------------------------------------------
// IndependenceCopulatedUniformDeviates
// 
MAGNET_DESCRIPTION("Returns a matrix of independence copulated uniform deviates")
MAGNET_PARAMETER("",nbNames,"nbNames must be an integer >0","number of names in the simulation")
MAGNET_PARAMETER("",nbPaths,"nbPaths must be an integer >0","number of paths in the simulation")
MAGNET_FUNCTION2(Matrix<double>, IndependenceCopulatedUniformDeviates,
                 long nbNames,
                 long nbPaths)
{
    Matrix<double>out = Matrix<double>(nbPaths,nbNames);
    IndependenceCopulatedUniformDeviates(&out[0][0],
                                    nbNames,
                                    nbPaths);
    return out;
}

// --------------------------------------------------------------------------
// Version
//
MAGNET_DESCRIPTION("returns the version of the library")
MAGNET_FUNCTION0(CM::String, Version)
{
    return "Copula Lib 2.1";
}
