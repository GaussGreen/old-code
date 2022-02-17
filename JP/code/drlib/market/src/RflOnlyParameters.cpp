//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Description : Container class for RFL-specific per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/PiecewiseFlatMappingFunction.hpp"
#include "edginc/PiecewiseFlatIncrementalMappingFunction.hpp"
#include "edginc/RFLParameters.hpp"
#include "edginc/PiecewiseLinearMappingFunction.hpp"
#include "edginc/Nrfns.hpp" //for zbrentUseful root finder
#include "edginc/RootFinder.hpp" // for RtSafe root finder
#include "edginc/PortfolioName.hpp" // for the betaTweak function
#define QLIB_RFLONLYPARAMETERS_CPP
#include "edginc/RflOnlyParameters.hpp"
#include "edginc/ClientRunnable.hpp"

/* For the Integ Method */
#define CCM_INTEG_LOW -7.5
#define CCM_INTEG_UP   7.5
#define CCM_INTEG_NB   101

DRLIB_BEGIN_NAMESPACE

/** setupBetas :
 * ugly function to extend the beta vector for piecewise flat betas
 * in case the same nb of betas and thresholds are specified */
static void setupBetas(CDoubleArrayConstSP thresholdValuesInitial, // I
                       CDoubleArrayConstSP betaValuesInitial,      // I
                       CDoubleArray& thresholdValues,              // O
                       CDoubleArray& betaValues)                   // O
{
    int nbBeta = betaValuesInitial->size();
    if (thresholdValuesInitial->size() - betaValuesInitial->size() > 1) {
        throw ModelException("RflOnlyParameters::setupBetas",
                             "inconsistent array size for thresholds and "
                             "beta values");
    }

    // VERY UGLY : for flat  interpolation, if n thresholds
    // and n betas are supplied, we assume the first value in the beta array is valid before the
    // first threshold value and between first and second threshold so
    // we build an extended beta array reflecting that (the first 2 values will be the same).
    // if n+1 betas are supplied, we assume btea[i] it is valid on the intervals over
    //  ] threshold[i-1] , threshold[i] ]
    thresholdValues.resize(thresholdValuesInitial->size());
    if (thresholdValuesInitial->size() == betaValuesInitial->size()) {
        // extend the beta array so that it matches nbThreshold + 1
        nbBeta += 1;
        betaValues.resize(nbBeta);

        betaValues[0] = (*betaValuesInitial)[0];
        for (int i= 0; i<nbBeta-1; ++i) {
            betaValues[i+1]      = (*betaValuesInitial)[i];
            thresholdValues[i]   = (*thresholdValuesInitial)[i];
        }
    }
    else {
        betaValues.resize(nbBeta);
        for (int i = 0; i<nbBeta-1; ++i) {
            betaValues[i] = (*betaValuesInitial)[i];
            thresholdValues[i]   = (*thresholdValuesInitial)[i];
        }
        betaValues[nbBeta-1] = (*betaValuesInitial)[nbBeta-1];
    }
}


/** Computes the mean of Beta(M)*M assuming beta is piecewise flat */
static double MeanCF(int nbBeta, 
                     const CDoubleArray& marketThreshold, 
                     const CDoubleArray& betaValue, 
                     bool isTweak, 
                     double betaHist)
{
    double temp=0;
    double dens1=0,dens2, beta;
    int i;
    
    for (i=0;i<nbBeta-1;++i) {
        dens2=N1Density(marketThreshold[i]);
        if (isTweak) {
            beta = PortfolioName::betaTweak(betaHist, betaValue[i]);
        }
        else {
            beta = betaValue[i];
        }
        temp=temp+beta*(dens2-dens1);
        dens1=dens2;
    }

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValue[nbBeta-1]);
    }
    else {
        beta = betaValue[nbBeta-1];
    }
    return temp - beta*dens1;
}

/** Computes the mean of Beta(M)*M numerically assuming linear interpolation */
static double MeanNum(const RflOnlyParameters *p, double betaHist)
{
    double mean = 0;
    /* set up new market factor containers */
    GaussianIntegrationMethod integMethod = 
        GaussianIntegrationMethod(CCM_INTEG_LOW, CCM_INTEG_UP, CCM_INTEG_NB);
    double M = integMethod.l;

    for (int i=0; i<integMethod.n; i++) {
        double beta;
        if (p->isTweak) {
            beta = PortfolioName::betaTweak(betaHist, (p->betaCurve)->map(M));
        }
        else {
            beta = (p->betaCurve)->map(M);
        }
        mean += M*beta*integMethod.weightCalc(M);
        M += integMethod.step;
    }
    return mean;
}


/** Computes the variance of Beta(M)*M assuming beta is piecewise flat */
static double VarianceCF(int nbBeta, 
                         const CDoubleArray& marketThreshold, 
                         const CDoubleArray& betaValue, 
                         double mean, 
                         bool isTweak, 
                         double betaHist)
{
    double temp=0;
    double dens1,dens2,cum1,cum2, beta;
    int i;
    
    dens1=dens2=N1Density(marketThreshold[0]);
    cum1=cum2=N1(marketThreshold[0]);
    
    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValue[0]);
    }
    else {
        beta = betaValue[0];
    }
    temp = beta*beta*(-marketThreshold[0]*dens2+cum2);

    for (i=1; i<nbBeta-1; ++i) {
        dens2=N1Density(marketThreshold[i]);
        cum2 =N1(marketThreshold[i]);
        if (isTweak) {
            beta = PortfolioName::betaTweak(betaHist, betaValue[i]);
        }
        else {
            beta = betaValue[i];
        }
        temp=temp+beta*beta*(marketThreshold[i-1]*dens1-marketThreshold[i]*dens2+cum2-cum1);
        dens1=dens2;
        cum1=cum2;
    }

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValue[nbBeta-1]);
    }
    else {
        beta = betaValue[nbBeta-1];
    }
    temp = temp + beta*beta*(marketThreshold[nbBeta-2]*dens1+1-cum1) - mean*mean;
    
    return temp;
}

/** Computes the variance of Beta(M)*M numerically */
static double VarianceNum(const RflOnlyParameters *p, double betaHist) {
    double e2 = 0;
    double mean = MeanNum(p, betaHist);
    /* set up new market factor containers */
    GaussianIntegrationMethod integMethod = 
        GaussianIntegrationMethod(CCM_INTEG_LOW,CCM_INTEG_UP,CCM_INTEG_NB);
    double M = integMethod.l;

    for (int i=0; i<integMethod.n; i++) {
        double beta;
        if (p->isTweak) {
            beta = PortfolioName::betaTweak(betaHist, (p->betaCurve)->map(M));
        }
        else {
            beta = (p->betaCurve)->map(M);
        }
  
        e2 += M*beta*M*beta*integMethod.weightCalc(M);
        M += integMethod.step;
    }
    return e2 - mean*mean;
}

static double XCumProbaCFwithVarianceCorrection(double x, int nbBeta, 
                                                const CDoubleArray& thresholdValues, 
                                                const CDoubleArray& betaValues, 
                                                double mean, 
                                                double variance, 
                                                bool isTweak, 
                                                double betaHist)
{
    double temp, a, b, bi1, bi2, beta, aux, aux2;
    int j;

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[0]);
    }
    else {
        beta = betaValues[0];
    }
    
    aux = sqrt(1.-variance);
    a = -beta/aux;
    b = (x-mean)/aux;
    
    aux2 = sqrt(1+a*a);
    bi2 = N2(thresholdValues[0],b/aux2,-a/aux2);
    
    temp = bi2;

    for (j=1;j<nbBeta-1;++j) {
        if (isTweak) {
            beta = PortfolioName::betaTweak(betaHist, betaValues[j]);
        }
        else {
            beta = betaValues[j];
        }

        a = -beta/aux;
        aux2 = sqrt(1+a*a);

        bi1 = N2(thresholdValues[j-1],b/aux2,-a/aux2);
        bi2 = N2(thresholdValues[j],b/aux2,-a/aux2);

        temp = temp + bi2 - bi1;
    }

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[nbBeta-1]);
    }
    else {
        beta = betaValues[nbBeta-1];
    }

    a = -beta/aux;
    aux2 = sqrt(1+a*a);

    bi1 = N2(-thresholdValues[nbBeta-2], b/aux2, a/aux2);
    
    temp = temp + bi1;
    return temp;
}

static double XCumProbaCFwithoutVarianceCorrection(double x, 
                                                   int nbBeta, 
                                                   const CDoubleArray& thresholdValues, 
                                                   const CDoubleArray& betaValues, 
                                                   double mean, 
                                                   bool isTweak, 
                                                   double betaHist)
{
    double temp, bi1, bi2, beta;
    int j;

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[0]);
    }
    else {
        beta = betaValues[0];
    }
    
    bi2 = N2(thresholdValues[0], x-mean, beta);
    
    temp = bi2;

    for (j=1;j<nbBeta-1;++j) {
        if (isTweak) {
            beta = PortfolioName::betaTweak(betaHist, betaValues[j]);
        }
        else {
            beta = betaValues[j];
        }

        bi1 = N2(thresholdValues[j-1], x-mean, beta);
        bi2 = N2(thresholdValues[j], x-mean, beta);

        temp = temp + bi2 - bi1;
    }

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[nbBeta-1]);
    }
    else {
        beta = betaValues[nbBeta-1];
    }

    bi1 = N2(-thresholdValues[nbBeta-2], x-mean, -beta);
    
    temp = temp + bi1;
    return temp;
}


static double XCumProbaCF(double x, 
                          int nbBeta, 
                          const CDoubleArray& thresholdValues, 
                          const CDoubleArray& betaValues, 
                          double mean, 
                          double variance, 
                          bool isTweak, 
                          bool correctMean, 
                          bool correctVariance, 
                          double betaHist)
{
    if (correctVariance) {
        return XCumProbaCFwithVarianceCorrection(
            x, nbBeta, thresholdValues, betaValues, correctMean ? mean : 0.0, variance, isTweak, betaHist);
    }
    else {
        return XCumProbaCFwithoutVarianceCorrection(
            x, nbBeta, thresholdValues, betaValues, correctMean ? mean : 0.0, isTweak, betaHist);
    }
}

static double XCumProbaNum(double x, const RflOnlyParameters *p, double betaHist)
{
    double mean     = p->correctMean?MeanNum(p, betaHist):0.;
    double variance = p->correctVariance?VarianceNum(p, betaHist):1.;
    double fx = 0;
    /* set up new market factor containers */
    GaussianIntegrationMethod integMethod = 
        GaussianIntegrationMethod(CCM_INTEG_LOW,CCM_INTEG_UP,CCM_INTEG_NB);
    double M = integMethod.l;

    for (int i=0; i<integMethod.n; i++) {
        double beta;
        if (p->isTweak) {
            beta = PortfolioName::betaTweak(betaHist, (p->betaCurve)->map(M));
        }
        else {
            beta = (p->betaCurve)->map(M);
        }
        // if variance is not corrected then we take the pure beta*beta value
        variance = p->correctVariance?variance:beta*beta;
        fx += N1((x - beta*M - mean)/sqrt(1.-variance))*integMethod.weightCalc(M);
        M  += integMethod.step;
    }
    return fx;
}



static double XDensProbaCFwithVarianceCorrection(double x, int nbBeta, 
                                                const CDoubleArray& thresholdValues, 
                                                const CDoubleArray& betaValues, 
                                                double mean, 
                                                double variance, 
                                                bool isTweak, 
                                                double betaHist)
{
    double temp, a, b, bi1, bi2, beta, aux, aux2, aux3;
    int j;

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[0]);
    }
    else {
        beta = betaValues[0];
    }
    
    aux = sqrt(1.-variance);
    a = -beta/aux;
    b = (x-mean)/aux;
    aux2 = sqrt(1+a*a);
    aux3 = a*b/aux2;

    temp = N1Density(b/aux2) * N1(aux2*thresholdValues[0]+aux3) / aux2;

    for (j=1;j<nbBeta-1;++j) {
        if (isTweak) {
            beta = PortfolioName::betaTweak(betaHist, betaValues[j]);
        }
        else {
            beta = betaValues[j];
        }

        a = -beta/aux;
        aux2 = sqrt(1+a*a);
        aux3 = a*b/aux2;

        bi1 = N1(aux2 * thresholdValues[j-1] + aux3);
        bi2 = N1(aux2 * thresholdValues[j]   + aux3);

        temp = temp + N1Density(b/aux2)* (bi2 - bi1) / aux2;
    }

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[nbBeta-1]);
    }
    else {
        beta = betaValues[nbBeta-1];
    }

    a = -beta/aux;
    aux2 = sqrt(1+a*a);
    aux3 = a*b/aux2;
    
    temp += N1Density(b/aux2)* (1. - N1(aux2 * thresholdValues[nbBeta-2] + aux3)) / aux2;
    temp /= aux;

    return temp;
}

static double XDensProbaCFwithoutVarianceCorrection(double x, 
                                                   int nbBeta, 
                                                   const CDoubleArray& thresholdValues, 
                                                   const CDoubleArray& betaValues, 
                                                   double mean, 
                                                   bool isTweak, 
                                                   double betaHist)
{
    double temp, bi1, bi2, beta, a, b, aux, aux2, aux3;
    int j;

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[0]);
    }
    else {
        beta = betaValues[0];
    }

    aux = sqrt(1.-beta*beta);
    a = -beta/aux;
    b = (x-mean)/aux;
    aux2 = sqrt(1+a*a);
    aux3 = a*b/aux2;

    temp = 1./(aux*aux2)*N1Density(b/aux2)* N1(aux2*thresholdValues[0]+aux3);

    for (j=1;j<nbBeta-1;++j) {
        if (isTweak) {
            beta = PortfolioName::betaTweak(betaHist, betaValues[j]);
        }
        else {
            beta = betaValues[j];
        }

        aux = sqrt(1.-beta*beta);
        a = -beta/aux;
        b = (x-mean)/aux;
        aux2 = sqrt(1+a*a);
        aux3 = a*b/aux2;

        bi1 = N1(aux2 * thresholdValues[j-1] + aux3);
        bi2 = N1(aux2 * thresholdValues[j]   + aux3);

        temp += 1./(aux*aux2)*N1Density(b/aux2)*(bi2 - bi1);
    }

    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, betaValues[nbBeta-1]);
    }
    else {
        beta = betaValues[nbBeta-1];
    }

    aux = sqrt(1.-beta*beta);
    a = -beta/aux;
    b = (x-mean)/aux;
    aux2 = sqrt(1+a*a);
    aux3 = a*b/aux2;

    temp += 1./(aux*aux2) * N1Density(b/aux2) * (1. - N1(aux2*thresholdValues[nbBeta-2]+aux3));

    return temp;
}


static double XDensProbaCF(double x, 
                          int nbBeta, 
                          const CDoubleArray& thresholdValues, 
                          const CDoubleArray& betaValues, 
                          double mean, 
                          double variance, 
                          bool isTweak, 
                          bool correctMean, 
                          bool correctVariance, 
                          double betaHist)
{
    if (correctVariance) {
        return XDensProbaCFwithVarianceCorrection(
            x, nbBeta, thresholdValues, betaValues, correctMean ? mean : 0.0, variance, isTweak, betaHist);
    }
    else {
        return XDensProbaCFwithoutVarianceCorrection(
            x, nbBeta, thresholdValues, betaValues, correctMean ? mean : 0.0, isTweak, betaHist);
    }
}

static double XDensProbaNum(double x, const RflOnlyParameters *p, double betaHist)
{
    double mean     = p->correctMean?MeanNum(p, betaHist):0.;
    double variance = p->correctVariance?VarianceNum(p, betaHist):1.;
    double fx = 0;
    /* set up new market factor containers */
    GaussianIntegrationMethod integMethod = 
        GaussianIntegrationMethod(CCM_INTEG_LOW,CCM_INTEG_UP,CCM_INTEG_NB);
    double M = integMethod.l;

    for (int i=0; i<integMethod.n; i++) {
        double beta;
        if (p->isTweak) {
            beta = PortfolioName::betaTweak(betaHist, (p->betaCurve)->map(M));
        }
        else {
            beta = (p->betaCurve)->map(M);
        }
        // if variance is not corrected then we take the pure beta*beta value
        variance = p->correctVariance?variance:beta*beta;
        fx += N1Density((x - beta*M - mean)/sqrt(1.-variance))*integMethod.weightCalc(M) / sqrt(1.-variance);
        M  += integMethod.step;
    }
    return fx;
}

////////////////////////// RflOnlyParameters methods //////////////////////////

/** Build an RflOnlyParameters object based on an old-style
    RFLParameters object - which must not be null */
RflOnlyParameters::RflOnlyParameters(RFLParametersConstSP rflParam) :
    RationalisedCreditEngineParameters(TYPE),
    betaCurve(rflParam->betaCurve),
    correctMean(rflParam->correctMean),
    correctVariance(rflParam->correctVariance),
	isTweak(rflParam->isTweak),
	useClosedForm(rflParam->useClosedForm),
    name(rflParam->getName())
{}


/** Method to compute the mean of Beta(M)*M in any case
 * using CF if beta is piecewise flat, otherwise numerical integration */
double RflOnlyParameters::mean(double betaHist) const {
    /* in case of piecewise flat mapping function, there is a close form for the mean and variance */
    if (isFlatBetaCurve() && useClosedForm) {
        PiecewiseMappingFunctionConstSP piecewiseBetaCurve(
            DYNAMIC_CAST(PiecewiseMappingFunction, betaCurve.get()));
 
        CDoubleArrayConstSP thresholdValuesInitial  = piecewiseBetaCurve->getX();
        CDoubleArrayConstSP betaValuesInitial       = piecewiseBetaCurve->getY();        
        CDoubleArray betaValues;
        CDoubleArray thresholdValues;

        setupBetas(thresholdValuesInitial, betaValuesInitial, thresholdValues, betaValues);
        
        int nbBeta = betaValues.size();

        // compute the variance of the beta(M)*M term
        return MeanCF(nbBeta, thresholdValues, betaValues, isTweak, betaHist);
    }
    else {
        return MeanNum(this, betaHist);
    }
}

/** method to compute the variance of Beta(M)*M in any case
 * using CF if beta is piecewise flat, otherwise numerical integration */
double RflOnlyParameters::variance(double betaHist) const {
    double meanValue;
    /* in case of piecewise flat mapping function, there is a close form for the mean and variance */
    if (isFlatBetaCurve() && useClosedForm) {
        PiecewiseMappingFunctionConstSP piecewiseBetaCurve(
            DYNAMIC_CAST(PiecewiseMappingFunction, betaCurve.get()));
 
        CDoubleArrayConstSP thresholdValuesInitial  = piecewiseBetaCurve->getX();
        CDoubleArrayConstSP betaValuesInitial       = piecewiseBetaCurve->getY();        
        CDoubleArray betaValues;
        CDoubleArray thresholdValues;
        setupBetas(thresholdValuesInitial, betaValuesInitial, thresholdValues, betaValues);
        
        int nbBeta = betaValues.size();

        // compute the expected value of beta(M)*M + sqrt(...)*Zi
        meanValue = MeanCF(nbBeta, thresholdValues, betaValues, isTweak, betaHist);

        // compute the variance of the beta(M)*M term
        return VarianceCF(nbBeta, thresholdValues, betaValues, meanValue, isTweak, betaHist);
    }
    else {
        return VarianceNum(this, betaHist);
    }
}

/** compute the conditional survival probability for a name */
void RflOnlyParameters::condSurvivalProba(double pind,
                                          double tgauss,
                                          double betaHist,
                                          double M,
                                          double& pm) const
{
    double meanValue, varianceValue, beta;

    // compute beta(M)
    if (isTweak) {
        beta = PortfolioName::betaTweak(betaHist, (betaCurve)->map(M));
    }
    else {
        beta = (betaCurve)->map(M);
    }

    // compute the expected value of beta(M)*M 
    meanValue = correctMean?mean(betaHist):0.;
    // compute the variance of the beta(M)*M term
    varianceValue = correctVariance?variance(betaHist):beta*beta;
    pm = pind * N1((tgauss - beta*M - meanValue) / sqrt(1.-varianceValue));
}

// structure used only in the threshold solving
// contains a pointer to the class itself
// and the value of the probability p for which
// to find the threshold Fx-1(p)
typedef struct thresholdStruct_ {
    const RflOnlyParameters *t;
    double         p;
    double         betaHist;
} thresholdStruct;

// function used in Brent solving for tgauss
// returns FX(tgauss) - pgauss
static double thresholdFunctionToSolve(double x, void  *s) {
    thresholdStruct *str = (thresholdStruct *) s;
    return str->t->XCumProba(x, str->betaHist) - str->p;
}

// function class used in Brent solving for tgauss
// returns FX(tgauss) - pgauss
class SolverHelper : public Func1D::WtDeriv {
public:
    SolverHelper(const RflOnlyParameters *t, double p, double betaHist) : 
      t(t), p(p), betaHist(betaHist) {};
    void operator()(double x, double& f, double& df) const
    {
        f  = t->XCumProba(x, betaHist) - p;
        df = t->XDensProba(x, betaHist);
    }
private:
    const RflOnlyParameters *t;
    double         p;
    double         betaHist;
};

/** compute the threshold associated to the survival probability for a name */
void RflOnlyParameters::threshold(double pgauss, 
                                  double &tgauss, 
                                  double betaHist) const
{
    static char method[] = "RflOnlyParameters::threshold";

#ifdef ZBRENT
    thresholdStruct str;
    str.t = this;
    str.p = pgauss;
    str.betaHist = betaHist;
    
    // solve for x s.t. XCumProba(x) = pgauss
    double rangeLow  = -10.;
    double rangeHigh =  10.;
    double funcLow;
    double funcHigh;
    ZbracReturn zbrac = zbracUseful(&thresholdFunctionToSolve,
                                    &str,
                                    &rangeLow,
                                    &rangeHigh,
                                    &funcLow,
                                    &funcHigh);

    // zbracUseful returns ZBRAC_SUCCESS if everything went well. Else:
    if (zbrac == ZBRAC_FUNC_ERROR) {
        throw NRException(method,
                          "Zbrac error evaluating the function in the "
                          "initial boundary [" + 
                          Format::toString(rangeLow) + ", " +
                          Format::toString(rangeHigh) + "].");
    }
    else if (zbrac == ZBRAC_BRAC_ERROR) {
        // Throw an NRException
        throw NRException(method,
                          "Zbrac error: Failed to identify the range"
                          " (last range tested: [" +
                          Format::toString("%.12f", rangeLow) + ", " + 
                          Format::toString("%.12f", rangeHigh) + "]).");
    }

    // Call zbrent
    try {
        tgauss = zbrentUsefulBoundary(
            &thresholdFunctionToSolve, // Function to find the root of
            &str,       // cds, for
            rangeLow,  // Low value for x
            rangeHigh, // High value for x
            1.e-10,    // Tolerance
            funcLow,
            funcHigh);
    }
    catch (NRException&) {
        // This is zbrent failing to solve - hide "zbrent" message,
        // which confuses users, and output a more meaningful error
        throw NRException(
             method,
             "Zbrent error: failed to find tgauss");
    }
#else
    SolverHelper helper = SolverHelper(this, pgauss, betaHist);
    
    RtSafe_solve(helper, 
                -7.5, 
                 7.5, 
                 1e-8, // accuracy 
                 true, // throw on Error
                 tgauss, 
                 RtSafe::default_MAXIT); // max iteration
        
#endif
}

/** compute the cumulative function of X = beta(M)*M+ ... */
double RflOnlyParameters::XCumProba(double x, double betaHist) const {
    double meanValue, varianceValue, fx;
    int nbBeta;

    // compute the expected value of beta(M)*M
    meanValue = correctMean?mean(betaHist):0.;
    // compute the variance of beta(M)*M 
    varianceValue = correctVariance?variance(betaHist):1;

    /* in case of piecewise flat mapping function, there is a close form for the theshold */
    if (isFlatBetaCurve() && useClosedForm) {
        PiecewiseMappingFunctionConstSP piecewiseBetaCurve(
            DYNAMIC_CAST(PiecewiseMappingFunction, betaCurve.get()));
 
        CDoubleArrayConstSP thresholdValuesInitial  = piecewiseBetaCurve->getX();
        CDoubleArrayConstSP betaValuesInitial       = piecewiseBetaCurve->getY();        
        CDoubleArray betaValues;
        CDoubleArray thresholdValues;
        setupBetas(thresholdValuesInitial, betaValuesInitial, thresholdValues, betaValues);
        nbBeta = betaValues.size();

        // compute the cumulative of X in closed form:
        fx = XCumProbaCF(x, nbBeta, thresholdValues, betaValues, meanValue, 
                         varianceValue, isTweak, correctMean, correctVariance, 
                         betaHist);
    }
    else { /* in any other case, we use a numerial inegration to compute the threshold */ 
        fx = XCumProbaNum(x, this, betaHist);
    }
    return fx;
}

/** compute the density function of X = beta(M)*M+ ... */
double RflOnlyParameters::XDensProba(double x, double betaHist) const {
    double meanValue, varianceValue, fx;
    int nbBeta;

    // compute the expected value of beta(M)*M
    meanValue = correctMean?mean(betaHist):0.;
    // compute the variance of beta(M)*M 
    varianceValue = correctVariance?variance(betaHist):1;

    /* in case of piecewise flat mapping function, there is a close form for the theshold */
    if (isFlatBetaCurve() && useClosedForm) {
        PiecewiseMappingFunctionConstSP piecewiseBetaCurve(
            DYNAMIC_CAST(PiecewiseMappingFunction, betaCurve.get()));
 
        CDoubleArrayConstSP thresholdValuesInitial  = piecewiseBetaCurve->getX();
        CDoubleArrayConstSP betaValuesInitial       = piecewiseBetaCurve->getY();        
        CDoubleArray betaValues;
        CDoubleArray thresholdValues;
        setupBetas(thresholdValuesInitial, betaValuesInitial, thresholdValues, betaValues);
        nbBeta = betaValues.size();

        // compute the cumulative of X in closed form:
        fx = XDensProbaCF(x, nbBeta, thresholdValues, betaValues, meanValue, 
                         varianceValue, isTweak, correctMean, correctVariance, 
                         betaHist);
    }
    else { /* in any other case, we use a numerial inegration to compute the threshold */ 
        fx = XDensProbaNum(x, this, betaHist);
    }
    return fx;
}


/** Destructor */
RflOnlyParameters::~RflOnlyParameters() 
{}


/** Returns the name of this object */
string RflOnlyParameters::getName() const {
    return name;
}

/** Returns true if underlying "betaCurve" is flat or piecewise flat:
 * in that case, there are closed form formulae for mean, variance and
 * integration of the conditional single-name survival probabilities
 * over the market factor. */
bool RflOnlyParameters::isFlatBetaCurve() const {
    return PiecewiseFlatMappingFunction::TYPE->isInstance(betaCurve) ||
        PiecewiseFlatIncrementalMappingFunction::TYPE->isInstance(betaCurve);
}


/** Invoked when Class is 'loaded' */
void RflOnlyParameters::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(RflOnlyParameters, clazz);
    SUPERCLASS(RationalisedCreditEngineParameters);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name, "Identifier");

    FIELD(correctMean, "Flag to correct or not the mean");
    FIELD(correctVariance, "Flag to correct or not the variance");
    FIELD(isTweak, "Flag set to true if the beta curve is relative "
                 "a tweak from the historical beta");
    FIELD_MAKE_OPTIONAL(isTweak);
    FIELD(useClosedForm, "Boolean set to true if closed form formulae "
                 "should be used when available (e.g. with piecewise flat "
                 "beta curve)");
    FIELD_MAKE_OPTIONAL(useClosedForm);
    FIELD(betaCurve, "Mapping function for beta as a function of market factor");
}

RflOnlyParameters::RflOnlyParameters() :
    RationalisedCreditEngineParameters(TYPE), 
    isTweak(false), 
    useClosedForm(true) 
{}
    
/** Default constructor */
IObject* RflOnlyParameters::defaultConstructor() {
    return new RflOnlyParameters();
}

/** TYPE (for reflection) */        
CClassConstSP const RflOnlyParameters::TYPE =
    CClass::registerClassLoadMethod("RflOnlyParameters",
                                    typeid(RflOnlyParameters),
                                    RflOnlyParameters::load);

DEFINE_TEMPLATE_TYPE(RflOnlyParametersWrapper);


// addin for unit test of the cumulative and its derivative 
class MARKET_DLL RFLProbaUnitTest:
    public CObject,
    public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    double x;
    double betaHist;
    RflOnlyParametersSP params;

    // EdrAction version of addin
    IObjectSP run() {
        DoubleArraySP result = DoubleArraySP(new DoubleArray(2)) ;
        (*result)[0] = params->XDensProba(x,betaHist);
        (*result)[1] = params->XCumProba(x,betaHist);
        return IObjectSP(result);
    }

    RFLProbaUnitTest(): CObject(TYPE) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RFLProbaUnitTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(x, "x value for the function");
        FIELD(betaHist, "historical beta");
        FIELD(params, "rfl parameters");
    }
    
    static IObject* defaultConstructor(){
        return new RFLProbaUnitTest();
    }
};

CClassConstSP const RFLProbaUnitTest::TYPE =
    CClass::registerClassLoadMethod(
        "RFLProbaUnitTest",
        typeid(RFLProbaUnitTest),
        load);

/* external symbol to allow class to be forced to be linked in */
bool RFLOnlyParametersLoad(){
    return (RFLProbaUnitTest::TYPE != 0);
}




DRLIB_END_NAMESPACE
