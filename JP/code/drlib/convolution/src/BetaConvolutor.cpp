//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BetaConvolutor.cpp
//
//   Description : Convolution Algorithm
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/BetaConvolutor.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/Nrfns.hpp"

DRLIB_BEGIN_NAMESPACE

// Class responsible to perform the beta approximation convolution.
// Would probably internally represent the "DiscreteDistribution"s in a different way
// (discretising "values" according to loss unit).

/** Constructor */
BetaConvolutor::BetaConvolutor(  
               double varCutoff,           // maximum variance value in % to be considered to be 0
               int    maxNbIter,           // maximum nb of iterations in the beta function algorithm
               double precision,           // precision in the beta function algorithm
               double flooring):           // flooring value of the internal variables close to 0 
CObject(TYPE),
varCutoff(varCutoff),
maxNbIter(maxNbIter),
precision(precision),
flooring(flooring)
{
}

/** empty Constructor */
BetaConvolutor::BetaConvolutor():
CObject(TYPE),
varCutoff(VAR_CUTOFF),
maxNbIter(MAX_NB_ITER),
precision(PRECISION),
flooring(FLOORING)
{
}

IObject* BetaConvolutor::defaultConstructor()
{
    return new BetaConvolutor();
}


/* Destructor */
BetaConvolutor::~BetaConvolutor()
{}

/** TYPE (for reflection) */        
CClassConstSP const BetaConvolutor::TYPE =
CClass::registerClassLoadMethod(
                                "BetaConvolutor",
                                typeid(BetaConvolutor),
                                BetaConvolutor::load);

const double BetaConvolutor::VAR_CUTOFF = 1e-10;
const int BetaConvolutor::MAX_NB_ITER   = 100;
const double BetaConvolutor::PRECISION  = 3.0e-7;
const double BetaConvolutor::FLOORING   = 1.0e-30;

void BetaConvolutor::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(BetaConvolutor, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IConvolutor);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(varCutoff,
        "Higher value of the variance to be considered as 0 "
        "[default value is " +
        Format::toString(VAR_CUTOFF) +
        "]");
    FIELD_MAKE_OPTIONAL(varCutoff);
    FIELD(maxNbIter,
        "Maximum nb of iteration in the beta function algorithm "
        "[default value is " +
        Format::toString(MAX_NB_ITER) +
        "]");
    FIELD_MAKE_OPTIONAL(maxNbIter);
    FIELD(precision,
        "Precsision in the beta function algorithm "
        "[default value is " +
        Format::toString(PRECISION) +
        "]");
    FIELD_MAKE_OPTIONAL(precision);
    FIELD(flooring,
        "Floor value of close to 0 internally computed variables "
        "[default value is " +
        Format::toString(FLOORING) +
        "]");
    FIELD_MAKE_OPTIONAL(flooring);
}

/**
* "Pure" convolution algorithm that returns the whole convoluted distribution
* [Implements IConvolutor]
* */
// Result could be cached if this method is called more than once.
IDistribution1DConstSP BetaConvolutor::convolute(IDistribution1DArrayConstSP distributions) const
{
    throw ModelException(
        "BetaConvolutor::convolute",
        "Not implemented !");            
}


/************************************************************************
* Fast convolution algorithm
*/
/* ----------------------------------------------------------------------
* Beta incomplete function and closed form integral
*/
/*
*  A better approx of a binomial can be done with an incomplete
*  beta function. Utilities for incomplete beta function NRC p 227
* 
*  \todo use cephes as this seems to have troubles with numerical
*        values like mu=1e-9 s=1e-7... the cephes implementation
*        looks better.
*  log of gamma function
*/
/* incomplete beta function continued fraction expansion */
double BetaConvolutor::betacf(double a, double b, double x) const
{
    static char routine[] = "BetaConvolutor::betacf";

    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;

    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < flooring) d=flooring;
    d=1.0/d;
    h=d;
    for (m=1;m<=maxNbIter;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < flooring) d=flooring;
        c=1.0+aa/c;
        if (fabs(c) < flooring) c=flooring;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < flooring) d=flooring;
        c=1.0+aa/c;
        if (fabs(c) < flooring) c=flooring;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < precision) break;
    }
    if (m > maxNbIter) 
        throw ModelException(routine, "a or b too big, or maxNbIter too small."
                                      "Increase the nb of iterations or the variance cutoff.");

    return h;
}

/* incomplete beta function */
double BetaConvolutor::betai(double a, double b, double x) const
{
    static char routine[] = "BetaConvolutor::betai";

    double bt,y;
    double out;

    if (x < 0.0 || x > 1.0) 
        throw ModelException (routine, "Bad x=" + Format::toString(x));

    if (x == 0.0 || x == 1.0) 
        bt=0.0;
    else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));

    if (x < (a+1.0)/(a+b+2.0)) 
    {
        y = betacf(a,b,x);
        out = bt * y / a;
    }
    else
    {
        y = betacf(b,a,1.0-x);
        out = 1.0 - bt * y / b;
    }
    return out;
}

/**
* Compute the expected loss of a senior tranche
* Assuming that L follows a beta distribution
* with a given expected loss and variance.
*
* we use an "incomplete beta" integ close formula 
*/
double BetaConvolutor::betaTrancheLoss(
                               double K,           /* (I) lower strike   */
                               double min,         /* (I) maximum value  */
                               double max,         /* (I) minimum value  */
                               double mean,        /* (I) expected loss  */
                               double var) const   /* (I) loss variance  */
{
    double a,b,x, y1,y2;
    double e;
    double norm, normMean, normVar, normK;
    norm     = (max - min);
    if (norm == 0.0)
    {
        // in that case expect var to be 0.0 as well !
        if (var == 0.0)
        {
            return Maths::max(mean-K,0.) + Maths::min(K,0.);
        }
        else
        {
            throw ModelException (
                "BetaConvolutor::betaTrancheLoss",
                "Distribution is reduced to a single value (=" +
                Format::toString(min) +
                ") but has non zero variance (=" +
                Format::toString(var) +
                ").");
        }
    }
    
    normMean = (mean - min)   / norm;
    normVar  = var / (norm*norm);
    normK    = (K - min) / norm;

//    if (K > max)
//    {
//        return 0.0;
//    }
//
//    if (K < min)
//    {
//        return mean;
//    }
    
    /* strike is >= 1 */
    if (1.-normK<1e-10)
        return 0.0;
    /* strike is <= 0 */
    if (normK<1e-10)
        return mean;
    

    /* variance is close to 0 */
    if (normVar < varCutoff)
        return Maths::max(mean-K,0.) + Maths::min(K,0.);

    /* general case */
    x = normMean*(1.-normMean)/normVar - 1.; 
    if (x <= 0)
    {
        throw ModelException (
            "BetaConvolutor::betaTrancheLoss",
            "internal error");
    }

    a = normMean*x;
    b = (1.-normMean)*x;
    y1 = betai(b, a+1., 1.-normK);
    y2 = betai(b, a   , 1.-normK);
    e = Maths::min(K,0.) + norm * a/(a+b) * y1 + (min - K) * y2;
    return e;
}

/** [Implements IConvolutor] */
double BetaConvolutor::convoluteAndIntegrate(
    IDistribution1DArrayConstSP distributions,
    ICreditLossConfigConstSP lossConfig,
    const DateTime& timepoint) const
{
    const ITranchesCombinationPayoff* combinationPayoff =
        dynamic_cast<const ITranchesCombinationPayoff*>(lossConfig.get());
    if (combinationPayoff != 0)
    {
        // initialise "payoff" for given timepoint
        DoubleArray baseStrikes;
        DoubleArray baseStrikesWeights;
        double expectedLossOffset;
        combinationPayoff->linearDecomposition(
            timepoint,
            baseStrikes,
            baseStrikesWeights,
            expectedLossOffset);

        // compute the expected value and the variance of the convoluted distribution
        double minValue, maxValue, e, var;
        minValue = 0.;
        maxValue = 0.;
        e   = 0.;
        var = 0.;
        DoubleArrayConstSP tmpValues;
        for (int i = 0; i < (*distributions).size(); ++i)
        {
            tmpValues = (*distributions)[i]->discretise()->getValues();
            minValue += (*tmpValues)[0];
            maxValue += (*tmpValues)[(*tmpValues).size() - 1];
            e   += (*((*distributions)[i])).expectation();
            var += (*((*distributions)[i])).variance();
        }

        // integration
        double expectedLoss = expectedLossOffset;
        for (int k = 0; k < baseStrikes.size(); ++k) {
            // compute the expected loss of a tranche from k1 to maxValue
            double elSenior = betaTrancheLoss(baseStrikes[k], minValue, maxValue, e, var);
            expectedLoss += baseStrikesWeights[k] * (e - elSenior);
        }
        return expectedLoss;
    }
    else
    {
        throw ModelException (
        "BetaConvolutor::convoluteAndIntegrate",
        "Not implemented for payoffs other than linear combination of tranches.");
    }
}


/* external symbol to allow class to be forced to be linked in */
bool BetaConvolutorLoad(){
    return (BetaConvolutor::TYPE != 0);
}

DRLIB_END_NAMESPACE

