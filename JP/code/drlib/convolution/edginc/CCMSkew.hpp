//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMSkew.hpp
//
//   Description : 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_CCMSKEW_HPP
#define EDR_CCMSKEW_HPP

#include "edginc/AtomicArray.hpp"
DRLIB_BEGIN_NAMESPACE

class CONVOLUTION_DLL CCMSkew{
public:

    /**
       Inverse of the Cumulative density function F
       The implementation involves root solving for F(u)=Y
    */
    static double fInv(
        double p,        /* (I) survival proba (pgauss) */
        double beta,     /* (I) beta                    */
        double qM,       /* (I) skew of M               */
        int    nbPoint); /* (I) nb point for integ      */

    /**
       Cumulative density function of the variable X=beta M + srqt{1-beta^2} Z
       where M has the qM distribution and Z the qZ distribution.
       The implementation involves a numerical integral
       X = b \bar{H}_q(M) + \sqrt{1-b^2} Z
       Pr(X<T_1) = \int_M {Pr(X<T_1|M)n(M) dM}
    */
    static double f(
        double T,        /* (I) ccmFinv(pgauss)         */
        double beta,     /* (I) beta                    */
        double qM,       /* (I) skew of M               */
        int    nbPoint); /* (I) nb point for integ      */

    /**
       binormal cumulative density function for 2 skewed variables X and Y
       such that there is a q-mapping on M.
       X = b_1 \bar{H}_q1(M) + \sqrt{1-b_1^2} Z1
       Y = b_2 \bar{H}_q2(M) + \sqrt{1-b_2^2} Z2
       Pr(X<T_1, Y<T_2) = \int_M {Pr(X<T_1|M)Pr(Y<T_2|M)n(M) dM}
       
       Typically used to get joint survival proba (beware, this is
       not symmetric when qM != 0)
    */
    static double biCumF(
        double  T1,      /* (I) ccmFinv(pgauss1)         */
        double  T2,      /* (I) ccmFinv(pgauss2)         */
        double  b1,      /* (I) beta1                    */ 
        double  b2,      /* (I) beta2                    */ 
        double  q1,      /* (I) skew of M for name 1     */ 
        double  q2,      /* (I) skew of M for name 2     */ 
        long    nbPoint);/* (I) nb point for integ       */

    /**
       return conditional survival probability
       Pr(\tau > t) = pgauss pind
       pgauss = Pr(X<T)
       
       X = b \bar{H}_q(M) + \sqrt{1-b_^2} Z1
       Pr(X<T|M) = N(\{frac{T- b \bar{H}_q(M)}{\sqrt{1-b^2}})
    */
    static double condSurvivalProba(
        double pind,     /* (I) indep survival proba     */
        double tgauss,   /* (I) gauss survival threshold */
        double beta,     /* (I) beta                     */
        double qM,       /* (I) skew of M                */
        double M);       /* (I) M                        */
    
    /**
       cumulative density mapping function of M=Hqbar(Y)
       where Y is a standard normal variable
       The q map function is defined as H_q(y)=(exp(qy)-1)/q
       The normalized map function is \bar{H}_q(y) = H_q(y) / \sigma(H_q(Y))
    */
    static double normMap(double x, double q);
};

DRLIB_END_NAMESPACE
#endif
