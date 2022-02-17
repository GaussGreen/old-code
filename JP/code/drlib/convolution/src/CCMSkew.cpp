//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMSkew.cpp
//
//   Description : 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMSkew.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Nrfns.hpp"

DRLIB_BEGIN_NAMESPACE

/* -------------------------------------------------------------------------
** NormB
** this function is needed in the bounds of M (see below)
** determines a bound of the integral where Integrand is not = 0
*/
static double normB(double q) {
    if(fabs(q)<3e-16) {
        return -1e300;
    } else {
        return (q*q-1.)/sqrt(exp(q*q)*(exp(q*q)-1.));
    }
}

/* -------------------------------------------------------------------------
** normMapinv
** inverse q mapping 
** normalised to keep the variance equal to 1
*/
static double normMapinv(double y, double q){
    if(fabs(q)<3e-16) {
        return y;
    } else if(fabs(q-1.)<3e-16) {
        if(y*sqrt(exp(1.)*(exp(1.)-1.)) >0) {
            return (log(y*sqrt(exp(1.)*(exp(1.)-1.))));
        } else {
            /* we should never go there */
            return -1e300;
        }
    } else {
        if(Maths::sign(q)*y > normB(q)) {
            return (log(y*Maths::sign(q)*sqrt(exp(q*q)*(exp(q*q)-1.)) + 
                        1. - q*q)/q);
        } else {
            /* we should never go there */
            return -1e300;
        }
    }
}

/** for optimization of normMap */
static double normMapOpt(double q) {
    // avoid division by zero
    if (Maths::isPositive(q)){
        double expq2 = exp(q*q);
        return (Maths::sign(q)/sqrt(expq2 * (expq2 - 1.)));
    }
    return 1.0;
}

/* -------------------------------------------------------------------------
** normMapFast
**  q mapping from a normal variable
** we normalise the transformed variable to keep the variance equal to 1
** we keep the possibility to precalculate normMapFactor
*/
static double normMapFast(double x, double q, double normMapFactor){
    if (fabs(q)<3e-16){
        return x;
    }
    return ((q*q+ exp(q*x) - 1.) * normMapFactor);
}

/** root solving implementation */
typedef struct {double y; double beta; double qM; long nbPoint;} ParamStruct;
static double fBrent(double x, void *data) {
    ParamStruct *p = (ParamStruct*)data;
    return (CCMSkew::f(x, p->beta,p->qM,p->nbPoint) - p->y);
}

/**
 Utility that return cond proba given the threshold
 X = b \bar{H}_q(M) + \sqrt{1-b_^2} Z1
 Pr(X<T|M) = N(\{frac{T- b \bar{H}_q(M)}{\sqrt{1-b^2}})
 */
static double condProba(
    double T,        /* (I) Finv(pgauss)         */
    double beta,     /* (I) beta                 */
    double qM,       /* (I) skew of M            */
    double M,        /* (I) M                    */
    double normMapF){/* (I) optimization         */
    double c;
    if (fabs(beta) >.999999) {
        c = T-beta*normMapFast(M,qM,normMapF) > 0. ? 1. : 0.;
    } else {
        double threshold = (T-beta*normMapFast(M,qM,normMapF))/sqrt(1-beta*beta);
        c = threshold<8 ? N1(threshold) : 1.;
    }
    return c;
}


/**
   Inverse of the Cumulative density function F
   The implementation involves root solving for F(u)=Y
*/
double CCMSkew::fInv(
    double p,        /* (I) survival proba (pgauss) */
    double b,        /* (I) beta                    */
    double q,        /* (I) skew of M               */
    int    nbPoint){ /* (I) nb point for integ      */
    static char routine[] = "ccmFinv";

    /* deal with common cases */
    if (fabs(q)<3e-16 || fabs(b) < 3e-16){
        return N1InverseBetter(p);
    }
    if (b > .99999){
        return normMap(N1InverseBetter(p),q);
    }
    /* extreme values of p */
    if (p < 3e-16) {
        return normB(q);
    }
    if (1.-p < 3e-16) {
        return 3e100;
    }
    ParamStruct ps;
    ps.y    = p;
    ps.beta = b;
    ps.qM   = q;
    ps.nbPoint = nbPoint;
    // this zbrent doesn't take all the parameters. Need to review
    try{
        return zbrentUseful(fBrent, // Function to call  
                            &ps,
                            -1e6, 
                            1e6,
                            1e-8);
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/**
   Cumulative density function of the variable X=beta M + srqt{1-beta^2} Z
   where M has the qM distribution and Z the qZ distribution.
   The implementation involves a numerical integral
   X = b \bar{H}_q(M) + \sqrt{1-b^2} Z
   Pr(X<T_1) = \int_M {Pr(X<T_1|M)n(M) dM}
*/
double CCMSkew::f(
    double T,        /* (I) ccmFinv(pgauss)         */
    double b,        /* (I) beta                    */
    double q,        /* (I) skew of M               */
    int    nbPoint){ /* (I) nb point for integ      */
    double low=-7.5, high=+7.5;
    double step = (high-low)/nbPoint;
    double integral = 0.;
    double sumW     = 0.;

    /* deal with common special cases */
    if (fabs(q) < 3e-16 || fabs(b) < 3e-16){
        return N1(T);
    }
    if (b > .99999){
        return N1(normMapinv(T,q));
    }
    double normMapFactor = normMapOpt(q);
    double M = low;
    for (int i=0; i< nbPoint; ++i, M+=step) {
        double c   = condProba(T,b,q,M,normMapFactor);
        double den = N1Density(M) * step;
        integral += c * den;
        sumW     += den;
    }

    ASSERT(fabs(sumW-1) < 1e-10);
    return integral;
}

/**
   binormal cumulative density function for 2 skewed variables X and Y
   such that there is a q-mapping on M.
   X = b_1 \bar{H}_q1(M) + \sqrt{1-b_1^2} Z1
   Y = b_2 \bar{H}_q2(M) + \sqrt{1-b_2^2} Z2
   Pr(X<T_1, Y<T_2) = \int_M {Pr(X<T_1|M)Pr(Y<T_2|M)n(M) dM}
       
   Typically used to get joint survival proba (beware, this is
   not symmetric when qM != 0)
*/
double CCMSkew::biCumF(
    double  T1,      /* (I) ccmFinv(pgauss1)         */
    double  T2,      /* (I) ccmFinv(pgauss2)         */
    double  b1,      /* (I) beta1                    */ 
    double  b2,      /* (I) beta2                    */ 
    double  q1,      /* (I) skew of M for name 1     */ 
    double  q2,      /* (I) skew of M for name 2     */ 
    long    nbPoint){/* (I) nb point for integ       */
    double low=-7.5, high=+7.5;
    double step = (high-low)/nbPoint;
    double normMapFactor1,normMapFactor2; /* optimization */

    if (fabs(q1)<3e-16 && fabs(q2)<3e-16){
        return N2(T1,T2,b1*b2);
    }
    normMapFactor1 = normMapOpt(q1);
    normMapFactor2 = normMapOpt(q2);
    double M = low;
    double sumW = 0.;
    double integral = 0.;
    for (int i = 0; i<nbPoint; ++i, M+=step)
    {
        double c1  = condProba(T1,b1,q1,M,normMapFactor1);
        double c2  = condProba(T2,b2,q2,M,normMapFactor2);
        double den = N1Density(M) * step;
        integral += c1 * c2 * den;
        sumW     += den;
    }

    ASSERT(fabs(sumW-1) < 1e-10);
    return integral;
}

/**
   return conditional survival probability
   Pr(\tau > t) = pgauss pind
   pgauss = Pr(X<T)
       
   X = b \bar{H}_q(M) + \sqrt{1-b_^2} Z1
   Pr(X<T|M) = N(\{frac{T- b \bar{H}_q(M)}{\sqrt{1-b^2}})
*/
double CCMSkew::condSurvivalProba(
    double pind,     /* (I) indep survival proba     */
    double tgauss,   /* (I) gauss survival threshold */
    double beta,     /* (I) beta                     */
    double qM,       /* (I) skew of M                */
    double M){       /* (I) M                        */
    double normMapFactor = (fabs(qM)<3e-16) ? 0. : normMapOpt(qM);
    double c = condProba(tgauss,beta,qM,M,normMapFactor);
    return c*pind;
}
    
/**
   cumulative density mapping function of M=Hqbar(Y)
   where Y is a standard normal variable
   The q map function is defined as H_q(y)=(exp(qy)-1)/q
   The normalized map function is \bar{H}_q(y) = H_q(y) / \sigma(H_q(Y))
*/
double CCMSkew::normMap(double x, double q){
    if(fabs(q)<3e-16){
        return x;
    }
    return normMapFast(x,q, normMapOpt(q));
}


DRLIB_END_NAMESPACE
