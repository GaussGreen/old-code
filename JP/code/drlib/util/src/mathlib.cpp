//----------------------------------------------------------------------------
//
//     File           : Mathlib.cpp
//     Author         : Ning Shen
//
//     Description    : maths support functions.
//                     This file does not require any external functions.
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Addin.hpp"  
#include <cstdio>



DRLIB_BEGIN_NAMESPACE

/*
*============================================================================
*
* function definitions
*
*============================================================================
*/

/*------------------------------------------------------------------------------
*
*   Name         :  N1Density
*
*   Description  :      a univariate standard normal distribution density function
*
*   Parameters   :  StandardVariate.
*
*   Returns      :  Probability density of the random variable at the
*                                       given StandardVariate value.
*
*------------------------------------------------------------------------------*/
double N1Density(double  StandardVariate)
{
    double  result;
    
    result = (1.0/Maths::ROOT_TWO_PI)*exp(-(StandardVariate*StandardVariate)/2.0);
    return result;
}

/*------------------------------------------------------------------------------
*
*   Name         :  N1
*
*   Description  :      Approximation to a univariate normal distribution
*                                       cumulative probability function
*
*   Parameters   :  StandardVariate             Number of standard deviations
*                                                                               above the mean. The domain is
*                                                                               all double  numbers.
*
*   Returns      :  Probability of the random variable being below or
*                                       equal to the StandardVariate value.
*
*------------------------------------------------------------------------------*/
// double N1(double  StandardVariate)
// {
//     const double  A1 =   0.2316419;
//     const double  B1 =   0.319381530;
//     const double  B2 =  -0.356563782;
//     const double  B3 =   1.781477937;
//     const double  B4 =  -1.821255978;
//     const double  B5 =      1.330274429;
// 
//     double  absoluteVariate;
//     double  dndd;
//     double  k1;
//     double  result;
// 
//     absoluteVariate = fabs(StandardVariate);
//     if (absoluteVariate > 8.5)
//         result = 1.0;
//     else if (Maths::isZero(absoluteVariate))
//         result = 0.5;
//     else
//     {
//         dndd = exp(- absoluteVariate * absoluteVariate / 2.0) / Maths::ROOT_TWO_PI;
//         k1 = 1 / (1 + (A1 * absoluteVariate));
//         result = 1 - (dndd *    ((((((((B5 * k1) + B4) * k1) + B3) * k1) + B2) * k1) + B1) * k1);
//     }
//     if (StandardVariate < 0.0)
//         result = 1.0 - result;
//     return result;
// }


double N1(double  StandardVariate)
{

   const double SQRT2   = 1.414213562373095049;     /* sqrt(2) */
   const double SQRTPI  = 1.772453850905516027;     /* sqrt(pi) */

   /* Coefficients in expressiion of erf(x) for -0.46875<=x<=0.46875 */
   const double P10 = 3209.377589138469472562;    /* Numerator */
   const double P11 = 377.4852376853020208137;
   const double P12 = 113.8641541510501556495;
   const double P13 = 3.161123743870565596947;
   const double P14 = 0.1857777061846031526730;
   const double Q10 = 2844.236833439170622273;   /* Denominator */
   const double Q11 = 1282.616526077372275645;
   const double Q12 = 244.0246379344441733056;
   const double Q13 = 23.60129095234412093499;
   const double Q14 = 1.0;

   /* Coefficients in expression of erfc(x) for 0.46875<=x<=4.0 */
   const double P20 = 1230.33935479799725272;  /* Numerator */
   const double P21 = 2051.07837782607146532;
   const double P22 = 1712.04761263407058314;
   const double P23 = 881.952221241769090411;
   const double P24 = 298.635138197400131132;
   const double P25 = 66.1191906371416294775;
   const double P26 = 8.88314979438837594118;
   const double P27 = 0.564188496988670089180;
   const double P28 = 2.15311535474403846343e-8;
   const double Q20 = 1230.33935480374942043;  /* Denominator */
   const double Q21 = 3439.36767414372163696;
   const double Q22 = 4362.61909014324715820;
   const double Q23 = 3290.79923573345962678;
   const double Q24 = 1621.38957456669018874;
   const double Q25 = 537.181101862009857509;
   const double Q26 = 117.693950891312499305;
   const double Q27 = 15.7449261107098347253;
   const double Q28 = 1.0;

   /* Coefficients in expression of erfc(x) for x>= 4.0 */
   double P30 = -6.58749161529837803157E-4;    /* Numerator */
   double P31 = -1.60837851487422766278E-2;
   double P32 = -1.25781726111229246204E-1;
   double P33 = -3.60344899949804439429E-1;
   double P34 = -3.05326634961232344035E-1;
   double P35 = -1.63153871373020978498E-2;
   double Q30 =  2.33520497626869185443E-3;   /* Denominator */
   double Q31 =  6.05183413124413191178E-2;
   double Q32 =  5.27905102951428412248E-1;
   double Q33 =  1.87295284992346047209;
   double Q34 =  2.56852019228982242072;
   double Q35 =  1.0;

   double numerator;            /* Numerator of polynomial in expression */
   double denominator;          /* Denominator of polynomial in expression */
   double y;                    /* y = abs(x)/sqrt(2) */
   double y2;                   /* y*y */
   double erf;                  /* Error function value */
   double erfc;                 /* Complimentary Error function value */

   y  = fabs(StandardVariate) / SQRT2;
   y2 = y * y;

   if (y < 0.46875)
   {
      numerator   = P10 + y2*(P11 + y2*(P12 + y2*(P13 +y2*P14)));
      denominator = Q10 + y2*(Q11 + y2*(Q12 + y2*(Q13 +y2*Q14)));
      erf = y * numerator / denominator;
      return (StandardVariate>0.0) ? 0.5 + 0.5*erf : 0.5 - 0.5*erf;
   }
   else if (y < 4.0)
   {
      numerator   = P20 + y*(P21 + y*(P22 + y*(P23 +
                          y*(P24 + y*(P25 + y*(P26 + y*(P27 + y*P28)))))));
      denominator = Q20 + y*(Q21 + y*(Q22 + y*(Q23 +
                          y*(Q24 + y*(Q25 + y*(Q26 + y*(Q27 + y*Q28)))))));
      erfc = exp(-y2) * numerator / denominator;
      return (StandardVariate>0.0) ? 1.0 - 0.5*erfc : 0.5*erfc;
   }
   else /* (y > 4.0) */
   {
      double z2 = 1/y2;
      numerator   = P30 + z2*(P31 + z2*(P32 + z2*(P33 + z2*(P34 +z2*P35))));
      denominator = Q30 + z2*(Q31 + z2*(Q32 + z2*(Q33 + z2*(Q34 +z2*Q35))));
      erfc = (exp(-y2)/y) * (1.0 / SQRTPI + numerator / (denominator * y2));
      return (StandardVariate>0.0) ? 1.0 - 0.5*erfc : 0.5*erfc;
   }
}


// constructor of GaussianIntegrationMethod 
GaussianIntegrationMethod::GaussianIntegrationMethod(double lowerBound, double upperBound, long nbStep)
    :l(lowerBound), u(upperBound), n(nbStep), step((u-l)/(n-1))
{
}

double GaussianIntegrationMethod::weightCalc(double m)
{
    return step * N1Density(m);
}

/*------------------------------------------------------------------------------
*
*   Name         :  N1Inverse
*
*   Description  :      Inverse of the NormalCumulative function
*
*   Parameters   :  Probability         
*
*   Returns      :  Value of the random variable 
*                                       or -999 for invalid probability value
*
*------------------------------------------------------------------------------*/
double N1Inverse(double  Probability)
{
    const double  C0 = 2.53429149812948000;
    const double  C1 = 0.75131810118138500;
    const double  C2 = 0.01032544750876560;
    const double  C3 = 0.02392856341090320;
                
    const double  D1 = 1.41260322359690000;
    const double  D2 = 0.19190495008513600;
    const double  D3 = 0.00130792050519774;
    const double  D4 = 0.00929768260247225;
        
    double  t;
    double  Result;
    int             OverHalf = 0;
                     
    if (Probability > 1.0 || Probability < 0.0) 
        return -999;                                     
    if (Maths::isZero(Probability))
        return -9.0;
    if (Maths::isZero(1.0-Probability))
        return 9.0;
        
    if (Probability > 0.5)
    {
        OverHalf = 1;
        Probability = 1.0 - Probability;
    }
                
    t = sqrt(-2 * log(Probability));
        
    Result = - t + (C0+t*(C1+t*(C2+t*C3)))/(1+t*(D1+t*(D2+t*(D3+t*D4))));
        
    if (OverHalf) Result = - Result;

    return Result;
}
/// *** UPDATED VERSIONS *** ///////////////////////////////////////////////////////

/*****  NormalHP  ************************************************************/
/*                                                                           */
/*      Normal cumulative distribution according to W. J. Cody               */
/*      This function evaluates near-minimax approximations from             */
/*      "Rational Chebyshev approximations for the error function",          */
/*      Math. Comp., 1969, PP. 631-638.                                      */
/*                                                                           */
/*      The original algorithm uses rational functions that theoretically    */
/*      approximate erf(x) to at least 18 significant decimal digits.        */
/*      The result is then mapped into the cumulative normal distribution    */
/*      N(u):=(erfc(-u/sqrt(2))/2;                                           */
/*                                                                           */
double  NormalHP (double u)
{

    const double SQRT2   = 1.414213562373095049;     /* sqrt(2)  */
    const double SQRTPI  = 1.772453850905516027;     /* sqrt(pi) */

    const double a0 = 1.161110663653770e-002;
    const double a1 = 3.951404679838207e-001;
    const double a2 = 2.846603853776254e+001;
    const double a3 = 1.887426188426510e+002;
    const double a4 = 3.209377589138469e+003;

    const double b0 = 1.767766952966369e-001;
    const double b1 = 8.344316438579620e+000;
    const double b2 = 1.725514762600375e+002;
    const double b3 = 1.813893686502485e+003;
    const double b4 = 8.044716608901563e+003;

    const double c0 = 2.15311535474403846e-8;
    const double c1 = 5.64188496988670089e-1;
    const double c2 = 8.88314979438837594e00;
    const double c3 = 6.61191906371416295e01;
    const double c4 = 2.98635138197400131e02;
    const double c5 = 8.81952221241769090e02;
    const double c6 = 1.71204761263407058e03;
    const double c7 = 2.05107837782607147e03;
    const double c8 = 1.23033935479799725e03;

    const double d0 = 1.00000000000000000e00;
    const double d1 = 1.57449261107098347e01;
    const double d2 = 1.17693950891312499e02;
    const double d3 = 5.37181101862009858e02;
    const double d4 = 1.62138957456669019e03;
    const double d5 = 3.29079923573345963e03;
    const double d6 = 4.36261909014324716e03;
    const double d7 = 3.43936767414372164e03;
    const double d8 = 1.23033935480374942e03;

    const double p0 = 1.63153871373020978e-2;
    const double p1 = 3.05326634961232344e-1;
    const double p2 = 3.60344899949804439e-1;
    const double p3 = 1.25781726111229246e-1;
    const double p4 = 1.60837851487422766e-2;
    const double p5 = 6.58749161529837803e-4;

    const double q0 = 1.00000000000000000e00;
    const double q1 = 2.56852019228982242e00;
    const double q2 = 1.87295284992346047e00;
    const double q3 = 5.27905102951428412e-1;
    const double q4 = 6.05183413124413191e-2;
    const double q5 = 2.33520497626869185e-3;

    register double y, z;

    y = fabs(u);
    if (y <= 0.46875 * SQRT2) 
    {
        /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
        z = y * y;
        y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4)
              / ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);

        return (0.5 + y);
    }
 
    z = exp(- y*y/2) / 2;
    if (y <= 4.0 * SQRT2) 
    {
        /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
        y = y / SQRT2;
        y = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8)
          / ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

        y = z * y;
    } 
    else 
    {
        /* evaluate erfc() for |u| > sqrt(2)*4.0 */
        z = z * SQRT2 / y;
        y = 2 / (y * y);
        y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5)
              / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
        
        y = z * (1./SQRTPI - y);
    }

    return (u < 0.0 ? y : 1 - y);
}


/*****  Normal_InvHP  ********************************************************/
/*                                                                           */
/*      Inverse of normal cumulative distribution                            */
/*      This function calculates the normal deviate Z corresponding to a     */
/*      given lower tail area of p. Z is accurate to about 1 part in 10**16. */
/*                                                                           */
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, PP 477-484    */
/*                                                                           */
double  Normal_InvHP (double p)
{
    const double a0 = 3.3871328727963666080e0000;
    const double a1 = 1.3314166789178437745e0002;
    const double a2 = 1.9715909503065514427e0003;
    const double a3 = 1.3731693765509461125e0004;
    const double a4 = 4.5921953931549871457e0004;
    const double a5 = 6.7265770927008700853e0004;
    const double a6 = 3.3430575583588128105e0004;
    const double a7 = 2.5090809287301226727e0003;

    const double b1 = 4.2313330701600911252e0001;
    const double b2 = 6.8718700749205790830e0002;
    const double b3 = 5.3941960214247511077e0003;
    const double b4 = 2.1213794301586595867e0004;
    const double b5 = 3.9307895800092710610e0004;
    const double b6 = 2.8729085735721942674e0004;
    const double b7 = 5.2264952788528545610e0003;

    const double c0 = 1.42343711074968357734e000;
    const double c1 = 4.63033784615654529590e000;
    const double c2 = 5.76949722146069140550e000;
    const double c3 = 3.64784832476320460504e000;
    const double c4 = 1.27045825245236838258e000;
    const double c5 = 2.41780725177450611770e-01;
    const double c6 = 2.27238449892691845833e-02;
    const double c7 = 7.74545014278341407640e-04;

    const double d1 = 2.05319162663775882187e000;
    const double d2 = 1.67638483018380384940e000;
    const double d3 = 6.89767334985100004550e-01;
    const double d4 = 1.48103976427480074590e-01;
    const double d5 = 1.51986665636164571966e-02;
    const double d6 = 5.47593808499534494600e-04;
    const double d7 = 1.05075007164441684324e-09;

    const double e0 = 6.65790464350110377720e000;
    const double e1 = 5.46378491116411436990e000;
    const double e2 = 1.78482653991729133580e000;
    const double e3 = 2.96560571828504891230e-01;
    const double e4 = 2.65321895265761230930e-02;
    const double e5 = 1.24266094738807843860e-03;
    const double e6 = 2.71155556874348757815e-05;
    const double e7 = 2.01033439929228813265e-07;

    const double f1 = 5.99832206555887937690e-01;
    const double f2 = 1.36929880922735805310e-01;
    const double f3 = 1.48753612908506148525e-02;
    const double f4 = 7.86869131145613259100e-04;
    const double f5 = 1.84631831751005468180e-05;
    const double f6 = 1.42151175831644588870e-07;
    const double f7 = 2.04426310338993978564e-15;

    const double const1 = 0.180625;
    const double const2 = 1.6;

    const double split1 = 0.425;
    const double split2 = 5.0;

    register double out, q, r;


    if (p >  1.0 || p <  0.0) { return (-999); }           // not accepted input values //

    q = p - 0.5;
    if (p >= 1.0 || p <= 0.0) { return ( (q > 0.0) ? 9.0 : -9.0); }  // for p=1 and p=0 //  

    if (fabs(q) <= split1)
    {
        r = const1 - q*q;
        out = ( q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0 )
                  / (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + 1.0) );

        return (out);
    }
    else
    {
        r = (q < 0.0) ? p : 1.0 - p;
        r = sqrt(-log(r));
        if (r <= split2)
        {
            r -= const2;
            out = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0 )
                / (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + 1.0);
        }
        else
        {
            r -= split2;
            out = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0 )
                / (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + 1.0);           
        }
    }

    return (q < 0.0 ? -out : out);
}
/// *** END UPDATED VERSIONS *** ///////////////////////////////////////////////////////

/*
*------------------------------------------------------------------------------
*
*   Name         :  N2Density
*
*   Description  :      bivariate normal distribution density function
*
*   Parameters   :  StandardVariate1,StandardVariate2,Correlation
*
*   Returns      :  (double ) density of the random variables at the
*                                       given StandardVariate values.
*
*------------------------------------------------------------------------------*/

double N2Density(double  s1, double  s2, double  rho)
{
    double  result, var;
    var = 0.5*(s1*s1-2.0*rho*s1*s2+s2*s2)/(1-rho*rho);
    result = (1/(2*Maths::PI)/sqrt(1-rho*rho))*exp(-var);
    return(result);
}

#if 0
// new improved version below

/*------------------------------------------------------------------------------
*
*   Name         :  N2
*
*   Description  : Returns Bivariate distribution using approximation in HULL.with some modifications
*
*   Parameters   : Random variables 1 & 2 plus correlation.
*
*   Returns      : Cumulative probability for bivariate normal deviates.
*
*
*------------------------------------------------------------------------------*/
/************ local functions for BiNormalCumulative  */
double f_BiNormalCumulative (double x, double y, double a, double b, double r)
{
    double ans;
    ans = exp(a * (2.0 * x - a) + b * (2.0 * y - b) + 2.0 * r * (x - a) * (y - b));
    return (ans);
}
double m_BiNormalCumulative (double a, double b, double r)
{
    short i, j;
    double c[4]={0.325303, 0.4211071, 0.1334425, 0.006374323};
    double d[4]={0.1337764,0.6243247,1.3425378, 2.2626645};
    double k, kk, sum, aa, bb, ans;

    /* trap for correl =1 or -1     */
    if (Maths::isZero(1.0-r) || Maths::isZero(1.0+r))
    {return 0.0;}

    k = sqrt(1.0 - r * r) / Maths::PI;
    kk = 1.0 / sqrt(2.0 * (1.0 - r * r));
    sum = 0.0;
    aa = a * kk;
    bb = b * kk;
    for (i=0; i<=3; i++)
        for (j=0; j<=3; j++)
            sum += c[i] * c[j] * f_BiNormalCumulative(d[i], d[j], aa, bb, r);
    ans = sum * k;
    return (ans);
}
/********** end of local static functions for BiNormalCumulative  */
double N2 (double  a, double  b, double  r)
{
    double ans, denom, r1, r2, del, ma, mb;

    //  rho = 0
    if (Maths::isZero(r))
    {
        ans = N1(a) * N1(b);
        return ans;
    }
    // rho =1 or -1
    if (Maths::isZero(1.0 - r))
        return N1(Maths::min(a,b));
    if (Maths::isZero(1.0 + r))
    {
        if (a + b <= 0.0)
            return 0.0;
        else
        {
            ans = N1(a) - N1(-b);
            return ans;
        }
    }
    // a=0 and b=0
    if (Maths::isZero(a) && Maths::isZero(b))
    {
        ans = 0.25 + asin(r)/2.0 / Maths::PI;
        return ans;
    }

    if (a * b * r <= 0.0)
    {
        if (a <= 0.0 && b <= 0.0 && r <= 0.0)
            ans = m_BiNormalCumulative(a, b, r);
        else if (a <= 0.0 && b >= 0.0 && r >= 0.0)
            ans = N1(a) - m_BiNormalCumulative(a, -b, -r);
        else if (a >= 0.0 && b <= 0.0 && r >= 0.0)
            ans = N1(b) - m_BiNormalCumulative(-a, b, -r);
        else
            ans = N1(a) + N1(b) - 1.0 + m_BiNormalCumulative(-a, -b, r);
    }
    else
    {
        if (Maths::isZero(a * a - 2.0 * a * b * r + b * b))
            denom = 0.00001;         // this is not a correct solution
        else
            denom = sqrt(a * a - 2.0 * a * b * r + b * b);
        r1 = (r * a - b) / denom * Maths::sign(a);
        r2 = (r * b - a) / denom * Maths::sign(b);
        del = (1.0 - Maths::sign(a) * Maths::sign(b)) / 4.0;

        if (a <= 0.0 && r1 <= 0.0)
            ma = m_BiNormalCumulative(a, 0.0, r1);
        else if (a >= 0.0 && r1 <= 0.0)
            ma = N1(a) - 0.5 + m_BiNormalCumulative(-a, 0.0, r1);
        else if (a >= 0.0 && r1 >= 0.0)
            ma = 0.5 - m_BiNormalCumulative(-a, 0.0, -r1);
        else
            ma = N1(a) - m_BiNormalCumulative(a, 0.0, -r1);

        if (b <= 0.0 && r2 <= 0.0)
            mb = m_BiNormalCumulative(b, 0.0, r2);
        else if (b >= 0.0 && r2 <= 0.0)
            mb = N1(b) - 0.5 + m_BiNormalCumulative(-b, 0.0, r2);
        else if (b >= 0.0 && r2 >= 0.0)
            mb = 0.5 - m_BiNormalCumulative(-b, 0.0, -r2);
        else
            mb = N1(b) - m_BiNormalCumulative(b, 0.0, -r2);

        ans = ma + mb - del;
    }
    return (ans);
}
#endif

/*------------------------------------------------------------------------------
*
*   Name         :  LNDensity
*
*   Description  :      log-normal density function
*
*   Parameters   :  x   log-normal variable at the required density point
*                                       mat time to maturity
*                                       vol     volatility
*                                       fwd (expected) value of x at mat
*
*   Returns      :  density of log-normal function
*
*------------------------------------------------------------------------------*/

double LNDensity(double  x, double mat, double vol, double fwd)
{
    double  result, var;
    var = log(x/fwd)+0.5*vol*vol*mat;
    var = 0.5*var*var/vol/vol/mat;
    result = 1.0/x/vol/sqrt(2.0*Maths::PI*mat)*exp(-var);
    return result;
}


/*------------------------------------------------------------------------------
*
*   Name         :  SolveQuadratic
*
*   Description  :      solve a*x^2 + b*x + c =0. only real roots solved
*       Results:                x1, x2  the two roots (x1 <= x2).
*
*   Returns      :  0   no real roots.
*                                       1       two identical roots
*                                       2       2 distinct real roots.
*
*------------------------------------------------------------------------------*/
int SolveQuadratic(double a, double b, double c, double *x1, double *x2)
{
    double  q;
    q = b*b - 4.0*a*c;
    if (q < 0)
        return 0;
    if (q==0)
    {
        *x1 = *x2 = -b/2/a;
        return 1;
    }
    q = sqrt(q);
    *x1 = (-b - q)/2/a;
    *x2 = (-b + q)/2/a;
    return 2;
}

/*------------------------------------------------------------------------------
*
*   Name         :  SolveCubic
*
*   Description  :      solve x^3 + a*x^2 + b*x + c =0. only real roots returned.
*       Results:                x1, x2, x3      if only one real root x1 is returned
*
*   Returns      :      1       one real root
*                                       2       2 real roots
*                                       3       3 real roots
*
*------------------------------------------------------------------------------*/
int SolveCubic(double a, double b, double c, double *x1, double *x2, double *x3)
{
    double  q, r, q3, s1, s2, mod, theta;
    q = b/3 - a*a/9;
    q3 = q*q*q;
    r = (a*b - 3*c)/6 - a*a*a/27;
    if (q3 + r*r == 0)
    {
        s1 = pow(r,1.0/3.0);
        *x1 = 2*s1 - a/3;
        *x2 = *x3 = -s1 - a/3;
        return 2;
    }
    if (q3 + r*r > 0)
    {
        mod = fabs(r+sqrt(q3+r*r));
        theta = Maths::sign(r+sqrt(q3+r*r));
        s1 = theta*pow(mod,1.0/3.0);
        mod = fabs(r-sqrt(q3+r*r));
        theta = Maths::sign(r-sqrt(q3+r*r));
        s2 = theta*pow(mod,1.0/3.0);
        *x1 = s1 + s2 - a/3;
        return 1;
    }
    mod = sqrt(r*r - (q3+r*r));
    theta = atan(sqrt(-q3-r*r)/r);
    mod = 2*pow(mod,1.0/3.0);
    s1 = -mod*cos(theta/3);
    s2 = mod*sin(theta/3);
        
    *x1 = s1 - a/3;
    *x2 = -s1/2 - a/3 + sqrt(3.0)/2*s2;
    *x3 = -s1/2 - a/3 - sqrt(3.0)/2*s2;
    return 3;
}

/*-----------------------------------------------------------------------------
*       Name:       :   QuadraticInterp
*                                       quadratic interpolation
*
*       Return          :       result
*------------------------------------------------------------------------------*/
double QuadraticInterp(double x, double x1, double x2, double x3, double y1, double y2, double y3)
{
    double result;

    if ((x1 == x2) || (x1 == x3) || (x2 == x3)){
        throw ModelException("QuadraticInterp", 
                    "wrong input in QuadraticInterp!.");           

    }else{
        result = y1 * (x - x2)*(x - x3)/(x1 - x2)/(x1 - x3)
            +y2 * (x - x1)*(x - x3)/(x2 - x1)/(x2 - x3)
            +y3 * (x - x2)*(x - x1)/(x3 - x1)/(x3 - x2);
    }


    return result;
}

/*-----------------------------------------------------------------------------
*       Name:       :   CubicInterp
*                                       cubic interpolation
*
*       Return          :       result
*------------------------------------------------------------------------------*/
double CubicInterp(double x, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4)
{
    double result;

    result = y1 * (x - x2)*(x - x3)*(x - x4)/(x1 - x2)/(x1 - x3)/(x1 - x4)
        +y2 * (x - x1)*(x - x3)*(x - x4)/(x2 - x1)/(x2 - x3)/(x2 - x4)
        +y3 * (x - x2)*(x - x1)*(x - x4)/(x3 - x1)/(x3 - x2)/(x3 - x4)
        +y4 * (x - x2)*(x - x1)*(x - x3)/(x4 - x1)/(x4 - x2)/(x4 - x3);

    return result;
}


/*-------------------------------------------------------------------------------*/

double LinearInterp(double x, double x1, double x2, double y1, double y2){

    double slope;
    double y;

    if (x2 - x1 == 0.0){
        if (y1 ==y2){
            y = y1;
        }else{
            throw ModelException("LinearInterp", 
                        "wrong input in Linear Interp!.");            
        }        
    }else{
        slope = (y2-y1)/(x2-x1);
        y = y1 + slope * (x-x1);
    }
    return y;
}


//linear interp near the boundary, and QuadraticInterp inside 
double interpF2(double x, double y, double* vX, double* vY, double* fxy, int dim_x, int dim_y){

    int ix, iy;
    double f_down, f_mid, f_up;
    double f;


    ix = Neighbour(x, vX, 0, dim_x -1,1);
    iy = Neighbour(y, vY, 0, dim_y -1,1);

    if (ix <= 0 || ix >= dim_x - 1 || iy <= 0 || iy >= dim_y -1){
        int ixx, iyy, ixx1, iyy1;

        if (ix <= 0) {
            ixx = 0;
        }
        else if (ix >= dim_x-1) {
            ixx = dim_x -2;
        }
        else {
            ixx = ix;
        }
        
        if (iy <= 0){
            iyy = 0;
        }
        else if(iy >= dim_y -1){
            iyy = dim_y -2;
        }
        else{
            iyy = iy;
        }

        ixx1 = ixx +1;
        iyy1 = iyy +1;
        
        f_down = LinearInterp(x, vX[ixx], vX[ixx1], fxy[iyy + ixx* dim_y], fxy[iyy + ixx1* dim_y]);
        f_mid = LinearInterp(x, vX[ixx], vX[ixx1], fxy[iyy1 + ixx* dim_y], fxy[iyy1 + ixx1* dim_y]);
        f = LinearInterp(y, vY[iyy], vY[iyy1], f_down, f_mid);

    }else{//inside
        f_down = QuadraticInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 vX[ix+1],
                                 fxy[(iy-1) + (ix-1) * dim_y], 
                                 fxy[(iy-1) + ix * dim_y], 
                                 fxy[(iy-1) + (ix+1) * dim_y]);

        f_mid = QuadraticInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 vX[ix+1],
                                 fxy[iy + (ix-1) * dim_y], 
                                 fxy[iy + ix * dim_y], 
                                 fxy[iy + (ix+1) * dim_y]);

        f_up = QuadraticInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 vX[ix+1],
                                 fxy[(iy+1) + (ix-1) * dim_y], 
                                 fxy[(iy+1) + ix * dim_y], 
                                 fxy[(iy+1) + (ix+1) * dim_y]);

        f = QuadraticInterp(y, 
                                 vY[iy-1], 
                                 vY[iy], 
                                 vY[iy+1],
                                 f_down, 
                                 f_mid, 
                                 f_up);

    }

    return f;

}

//linear interp near the boundary, and QuadraticInterp inside 
double interpF2(double x, double y, vector<double>& vX, vector<double>& vY, double** fxy){

    int ix, iy;
    int dim_x, dim_y;
    double f_down, f_mid, f_up;
    double f;

    dim_x = vX.size();
    dim_y = vY.size();

    ix = Neighbour(x, &vX[0], 0, dim_x -1,1);
    iy = Neighbour(y, &vY[0], 0, dim_y -1,1);

    if (ix <= 0 || ix >= dim_x - 1 || iy <= 0 || iy >= dim_y -1){
        int ixx, iyy, ixx1, iyy1;

        if (ix <= 0) {
            ixx = 0;
            
            if (iy <= 0){
                iyy = 0;
            }else if(iy >= dim_y -1){
                iyy = dim_y -2;
            }else{
                iyy = iy;
            }

        }else if (ix >= dim_x-1){
            ixx = dim_x -2;
            
            if (iy <= 0){
                iyy = 0;
            }else if(iy >= dim_y -1){
                iyy = dim_y -2;
            }else{
                iyy = iy;
            }
            iyy1 = iyy +1;
        }
        ixx1 = ixx +1;
        iyy1 = iyy +1;
        
        f_down = LinearInterp(x, vX[ixx], vX[ixx1], fxy[ixx][iyy], fxy[ixx1][iyy]);
        f_mid = LinearInterp(x, vX[ixx], vX[ixx1], fxy[ixx][iyy1], fxy[ixx1][iyy1]);
        f = LinearInterp(y, vY[iyy], vY[iyy1], f_down, f_mid);

    }else{//inside
        f_down = QuadraticInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 vX[ix+1],
                                 fxy[ix-1][iy-1], 
                                 fxy[ix][iy-1], 
                                 fxy[ix+1][iy-1]);

        f_mid = QuadraticInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 vX[ix+1],
                                 fxy[ix-1][iy], 
                                 fxy[ix][iy], 
                                 fxy[ix+1][iy]);

        f_up = QuadraticInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 vX[ix+1],
                                 fxy[ix-1][iy+1], 
                                 fxy[ix][iy+1], 
                                 fxy[ix+1][iy+1]);

        f = QuadraticInterp(y, 
                                 vY[iy-1], 
                                 vY[iy], 
                                 vY[iy+1],
                                 f_down, 
                                 f_mid, 
                                 f_up);

    }

    return f;

}

//linear interp 
double interpF2_linear(double x, double y, vector<double>& vX, vector<double>& vY, double** fxy){

    int ix, iy;
    int dim_x, dim_y;
    double f_down, f_mid, f_up;
    double f;

    dim_x = vX.size();
    dim_y = vY.size();

    ix = Neighbour(x, &vX[0], 0, dim_x -1,1);
    iy = Neighbour(y, &vY[0], 0, dim_y -1,1);

    if (ix <= 0 || ix >= dim_x - 1 || iy <= 0 || iy >= dim_y -1){
        int ixx, iyy, ixx1, iyy1;

        if (ix <= 0) {
            ixx = 0;
            
            if (iy <= 0){
                iyy = 0;
            }else if(iy >= dim_y -1){
                iyy = dim_y -2;
            }else{
                iyy = iy;
            }

        }else if (ix >= dim_x-1){
            ixx = dim_x -2;
            
            if (iy <= 0){
                iyy = 0;
            }else if(iy >= dim_y -1){
                iyy = dim_y -2;
            }else{
                iyy = iy;
            }
            iyy1 = iyy +1;
        }
        ixx1 = ixx +1;
        iyy1 = iyy +1;
        
        f_down = LinearInterp(x, vX[ixx], vX[ixx1], fxy[ixx][iyy], fxy[ixx1][iyy]);
        f_mid = LinearInterp(x, vX[ixx], vX[ixx1], fxy[ixx][iyy1], fxy[ixx1][iyy1]);
        f = LinearInterp(y, vY[iyy], vY[iyy1], f_down, f_mid);

    }else{//inside
        f_down = LinearInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 fxy[ix-1][iy-1], 
                                 fxy[ix][iy-1]);

        f_mid = LinearInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 fxy[ix-1][iy], 
                                 fxy[ix][iy]);

        f = LinearInterp(y, 
                                 vY[iy-1], 
                                 vY[iy], 
                                 f_down, 
                                 f_mid);

    }

    return f;
}

//linear interp  
double interpF2_linear(double x, double y, double* vX, double* vY, double* fxy, int dim_x, int dim_y){

    int ix, iy;
    double f_down, f_mid, f_up;
    double f;


    ix = Neighbour(x, vX, 0, dim_x -1,1);
    iy = Neighbour(y, vY, 0, dim_y -1,1);

    if (ix <= 0 || ix >= dim_x - 1 || iy <= 0 || iy >= dim_y -1){
        int ixx, iyy, ixx1, iyy1;

        if (ix <= 0) {
            ixx = 0;
        }
        else if (ix >= dim_x-1) {
            ixx = dim_x -2;
        }
        else {
            ixx = ix;
        }
        
        if (iy <= 0){
            iyy = 0;
        }
        else if(iy >= dim_y -1){
            iyy = dim_y -2;
        }
        else{
            iyy = iy;
        }

        ixx1 = ixx +1;
        iyy1 = iyy +1;
        
        f_down = LinearInterp(x, vX[ixx], vX[ixx1], fxy[iyy + ixx* dim_y], fxy[iyy + ixx1* dim_y]);
        f_mid = LinearInterp(x, vX[ixx], vX[ixx1], fxy[iyy1 + ixx* dim_y], fxy[iyy1 + ixx1* dim_y]);
        f = LinearInterp(y, vY[iyy], vY[iyy1], f_down, f_mid);

    }else{//inside
        f_down = LinearInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 fxy[(iy-1) + (ix-1) * dim_y], 
                                 fxy[(iy-1) + ix * dim_y] );

        f_mid = LinearInterp(x, 
                                 vX[ix-1], 
                                 vX[ix], 
                                 fxy[iy + (ix-1) * dim_y], 
                                 fxy[iy + ix * dim_y]);

        f = LinearInterp(y, 
                                 vY[iy-1], 
                                 vY[iy], 
                                 f_down, 
                                 f_mid);
    }

    return f;
}

/*---------------------------------------------------------------------------------*/












double atanh(double x)
{
#define ATANH_HUGE 1.0e10
    if(x < -1.0 + DBL_EPSILON) return - ATANH_HUGE;
    if(x > 1.0 - DBL_EPSILON) return ATANH_HUGE;

    return 0.5 * log( (1.0 + x) / (1.0 - x) );
}

// this is different to Maths::isZero
#define IS_ALMOST_ZERO(x) (fabs(x)<3e-15)?1:0
// and this is different to Maths::PI too
#ifndef PI
#define PI     3.14159265358979323846264338328      /* pi */
#endif

static double BiNormalCumAux (double a, double b, double rho, bool highAccuracy = true);
static double BNC1 (double a, double b, double rho, bool highAccuracy = true);
static double Phi1 (double a, double b, double rho, bool highAccuracy = true);
static double Phi (double a, double b, double rho);

/** Credit Hybrids ....
    Old Description: This function calculates the value for cumulative
    bivariate normal distribution by performing a numeric quadrature using
    transformation of variables.  See Analytics Library Technical Memo for
    mathematical background.
    
    New Description: The new method is based on an approximation that is
    well known and fast (Gauss quadrature). The weights and abscissae were
    computed to 15 significant digits using Mathematica. The accuracy for
    a 10 point quadrature in two dimensions is approximately 10^(-12).
*/
double N2(
    double a,            /* (I) normal point for first variable */
    double b,            /* (I) normal point for second variable */
    double rho,          /* (I) correlation of the two variables */
    bool highAccuracy) /* (I) if true/false, order 10/4 approximation will be used */
{         
    static char routine[]="N2";
    try{
        if (rho < -1. || rho > 1.) {
            throw ModelException(routine, 
                                 "abs(rho) ("+Format::toString(fabs(rho))+
                                 " >= 1.");
        }
        
        double biNormalCum = BiNormalCumAux(a, b, rho, highAccuracy);
        
        /*
         * Perhaps the routine returned something just outside the range 
         * [0,1]. This may be due to instability in GtoNormalCum.
         * We (silently) correct it here.
         */
        
        if ((biNormalCum >= -2.0) && (biNormalCum <= 3.0)) {
            biNormalCum = Maths::min(biNormalCum, 1.0);
            biNormalCum = Maths::max(biNormalCum, 0.0);
        } else {
            /*
             * We do a sanity check for numbers grossly out of range.
             * The routine should be stable enough that this case never 
             * actually occurs.
             */
            throw ModelException(routine, "Probability ("+
                                 Format::toString(biNormalCum)+
                                 ") outide the range [0,1] "
                                 "returned from integration");
        }
        return biNormalCum;
    } catch (exception& e){
        throw ModelException(e, routine,
                             Format::toString("Failed on inputs "
                                              "a=%g, b=%g, rho=%g", a, b, rho));
    }
}


/*
** FUNCTION: BiNormalCum
** AUTHOR:   Peter Taylor (27 July 1999)
**
** This is an algorithm extracted from a book. This is based on
** Drezner's method.
*/
static double BiNormalCumAux (double a, double b, double rho, bool highAccuracy) 
{
    double BNC; /* To be returned */

    if ((IS_ALMOST_ZERO(a)) &&
        (IS_ALMOST_ZERO(b)))
    {
        /*
         * An easy case where a and b are both zero.
         */
        return (0.25 + asin(rho)*(1./(2.*PI)));
    }

    /*
     * We trap the case where rho is either +1 or -1. In this case,
     * the bivariate normal degenerates.
     */
    if (IS_ALMOST_ZERO(fabs(rho) -  1.0))
    {
        /*
        ** Degenerate case - perfectly correlated.
        */
        if (rho > 0.0)
        {
            /*
            ** Perfectly correlated
            */
            BNC = N1(Maths::min (a,b));
        }
        else
        {
            /*
            ** Perfectly negatively correlated
            */
            if (a < -b)
            {
                BNC = 0.0;
            }
            else
            {
                BNC = N1(a) + N1(b) - 1.0;
            }
        }
        return BNC;
    }

    if (a > 15.)
        return N1(b);
    if (b > 15.)
        return N1(a);

    if (a*b*rho <= 0.0)
    {
        BNC = BNC1 (a,b,rho,highAccuracy);
    }
    else
    {
        /*
         * a*b*rho > 0 in this case.
         */
        double rho1;
        double rho2;
        double den2;
        double den;

        if (fabs(a) <= fabs(b))
        {
            double z = a/b;

            den2 = (z-1)*(z-1)+2*(1-rho)*z;
            if ((den2 <= 0) ||
                (IS_ALMOST_ZERO(den2)))
            {
                /*
                 * This case implies that z=1 or z=-1, which implies that
                 * rho=1 or rho=-1.
                 *
                 * This case is theoretically impossible, but we do a reasonable
                 * thing anyway.
                 */
                if (rho > 0.0)
                {
                    rho = 1.0;
                }
                else
                {
                    rho = -1.0;
                }
                return (BiNormalCumAux(a, b, rho, highAccuracy));
            }

            den = sqrt(den2);

            rho1 = (rho*z - 1)/den;
            rho2 = (rho   - z)/den;

            if (((a > 0.0) && (b < 0.0)) ||
                ((a < 0.0) && (b > 0.0)))
            {
                /*
                 * a and b are of opposite sign.
                 */
                rho1 = -rho1;
            }
        }
        else
        {
            double z = b/a;

            den2 = (z-1)*(z-1)+2*(1-rho)*z;
            if ((den2 <= 0) ||
                (IS_ALMOST_ZERO(den2)))
            {
                /*
                 * This case implies that z=1 or z=-1, which implies that
                 * rho=1 or rho=-1.
                 *
                 * This case is theoreticall impossible, but we do a reasonable
                 * thing anyway.
                 */
                if (rho > 0.0)
                {
                    rho = 1.0;
                }
                else
                {
                    rho = -1.0;
                }
                return (BiNormalCumAux(a, b, rho, highAccuracy));
            }

            den = sqrt(den2);

            rho1 = (rho   - z)/den;
            rho2 = (rho*z - 1)/den;

            if (((a > 0.0) && (b < 0.0)) ||
                ((a < 0.0) && (b > 0.0)))
            {
                /*
                 * a and b are of opposite sign.
                 */
                rho2 = -rho2;
            }
        }

        /*
         * Extra safety check to make absolutely sure rho1 and rho2 are in the 
         * correct range [-1, 1].
         */

        rho1 = Maths::min(rho1,  1.0);
        rho1 = Maths::max(rho1, -1.0);
        rho2 = Maths::min(rho2,  1.0);
        rho2 = Maths::max(rho2, -1.0);

        if (((a > 0.0) && (b > 0.0)) ||
            ((a < 0.0) && (b < 0.0)))
        {
            /*
             * a and b are of the same sign.
             */
            
            BNC = BNC1 (a,0.0,rho1,highAccuracy) + BNC1 (b,0.0,rho2,highAccuracy);
        }
        else
        {
            /*
             * a and b are of different sign.
             */
            BNC = BNC1 (a,0.0,rho1,highAccuracy) + BNC1 (b,0.0,rho2,highAccuracy) - 0.5;
        }
    }
    return BNC;
}


const double Smallest = 0.000000001;

/*!
   This function returns an approximate value for the cumulative bivariate
   normal distribution (with unit diagonal elements of the covariance matrix)
   when a, b, and c are all non-positive.
*/
double Phi2(double a, double b, double c)
{
  double ai[4], bi[4];

  // Define "magic" numbers 
  ai[0] = 0.3253030;  ai[1] = 0.4211071;  ai[2] = 0.1334425;  ai[3] = 0.006374323;
  bi[0] = 0.1337764;  bi[1] = 0.6243247;  bi[2] = 1.3425378;  bi[3] = 2.2626645;

  if (c > 1.0 - Smallest)
    c = 1.0 - Smallest;

  if (c < -1.0 + Smallest)
    c = - 1.0 + Smallest;

  double tmp = sqrt(1.0 - c * c);

  double atmp = a / (sqrt(2.0) * tmp);
  double btmp = b / (sqrt(2.0) * tmp);

  double res = 0.0;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      double tmp = a * (2.0 * bi[i] - atmp) + b * (2.0 * bi[j] - btmp) + 2.0 * c * (bi[i] - atmp) * (bi[j] - btmp);
      res += ai[i] * ai[j] * exp(tmp);
    }

  res *= tmp;
  res /= PI;

  return res;
}

static double Phi1(double a1, double b1, double rho, bool highAccuracy)
{
  return highAccuracy ? Phi(a1,b1,rho) : Phi2(a1,b1,rho);
}

static double BNC1 (double a, double b, double rho, bool highAccuracy)
{
    /*
    ** Only use if a * b * rho is zero or negative.
    */

    double BNC; /* To be returned */

    double a1;
    double b1;
    double c = sqrt (2.0 * (1.0 - rho * rho));

    /*
     * NOTE: We must make a check here, because this can be called with a
     * rho (very close to 1) from BiNormalCum that makes this c zero.
     */

    if (IS_ALMOST_ZERO (c))
    {
        /*
        ** Degenerate case - perfectly correlated.
        */
        if (rho > 0.0)
        {
            /*
            ** Perfectly correlated
            */
            BNC = N1(Maths::min (a,b));
        }
        else
        {
            /*
            ** Perfectly negatively correlated
            */
            if (a < -b)
            {
                BNC = 0.0;
            }
            else
            {
                BNC = N1 (a) + N1(b) - 1.0;
            }
        }
        return BNC;
    }

    /*
     * This is now a "safe" division, although the answer might be +Inf or -Inf.
     */
    a1 = a/c;
    b1 = b/c;

    /*
     * Note: It is acceptable for a1 or b1 to be either +Inf or -Inf.
     * The Phi routine handles this case correctly, because it gets called
     * with NEGATIVE arguments.
     */
    if (a <= 0.0 && b <= 0.0 && rho <= 0.0)
    {
        BNC = Phi1 (a1,b1,rho, highAccuracy);
    }
    else if (a <= 0.0 && b >= 0.0 && rho >= 0.0)
    {
        BNC = N1(a) - Phi1 (a1, -b1, -rho, highAccuracy);
    }
    else if (a >= 0.0 && b <= 0.0 && rho >= 0.0)
    {
        BNC = N1(b) - Phi1 (-a1, b1, -rho, highAccuracy);
    }
    else if (a >= 0.0 && b >= 0.0 && rho <= 0.0)
    {
        BNC = N1(a) + N1(b) - 1 + Phi1 (-a1, -b1, rho, highAccuracy);
    }
    else
    {
        /*
         * This case should not occur under any circumstance.
         */
        throw ModelException("mathlib::Phi", "Internal error");
    }

    return BNC;
}

/*
 * This routine performs a 2D numerical integration using a 10 point Gauss 
 * quadrature with weighting function Exp(-x^2) over the interval [0,Infinity].
 *
 * Previous methods integrated over [-Infinity,Infinity] and used the Hermite
 * polynomials's roots as the abscissae, but this leads to an instability at the 
 * origin. It is better to derive quadrature results for the half-infinite 
 * interval.
 *
 * The coefficients were generated by computing the orthogonal polynomials with
 * respect to the inner product <f,g>=Integral(0,Infinity,Exp(-x^2)*f(x)*g(x)).
 *
 * These were computed algebraically exact with Mathematica. Then, the roots of
 * this polynomial were found to 100 decimal places. After this, the weights
 * were easily determined (also to 100 decimal places).
 *
 * The resulting weights and abscissae were checked against the polynomials
 * x^i for i=0, 1, 2, ..., 19.
 *
 * The difference between the theoretical value of the integral and the 
 * quadrature value was on the order of 10^(-90).
 *
 * The values listed here were rounded to exactly 15 digits in scientific 
 * notation.
 *
 * The array ww[i][j] holds the result of multiplying weight i against weight
 * j (We precomputed these in Mathematica for speed).
 *
 * For details, see Z. Drezner, "Computation of the Bivariate Normal Integral"
 * Math. of Computation, Vol 32, Number 141, pp. 277-279 (1978).
 */

static double Phi (double a1, double b1, double rho)
{
    int    i;
    int    j;
    double sum;
    double BNC;

    static double x[10] = 
    {
        3.87385243256994e-2,
        1.98233304012949e-1,
        4.65201111814507e-1,
        8.16861885591907e-1,
        1.23454132402774e+0,
        1.70679814968865e+0,
        2.22994008892444e+0,
        2.80910374689825e+0,
        3.46387241949537e+0,
        4.25536180636561e+0,
    };                   /* Roots of the 10th degree orthogonal polynomial */

    static double ww[10][10] = 
    {
         {
             9.71251592540161e-3,
             2.05656611704430e-2,
             2.48402225754016e-2,
             1.95807584549757e-2,
             9.57910843855530e-3,
             2.66331292395503e-3,
             3.74956200622751e-4,
             2.25572193453660e-5,
             4.28242844331612e-7,
             1.22967113105880e-9,
         },
         {
             2.05656611704430e-2,
             4.35465354832845e-2,
             5.25976590213905e-2,
             4.14610639445273e-2,
             2.02831789389435e-2,
             5.63940297298984e-3,
             7.93946927345222e-4,
             4.77635386718788e-5,
             9.06778151287977e-7,
             2.60375375717955e-9,
         },
         {
             2.48402225754016e-2,
             5.25976590213905e-2,
             6.35300536271685e-2,
             5.00787233660737e-2,
             2.44990265668761e-2,
             6.81154979071513e-3,
             9.58968361136659e-4,
             5.76911640119526e-5,
             1.09525149314805e-6,
             3.14494254883637e-9,
         },
         {
             1.95807584549757e-2,
             4.14610639445273e-2,
             5.00787233660737e-2,
             3.94754669765184e-2,
             1.93118044788806e-2,
             5.36932834442915e-3,
             7.55924299324728e-4,
             4.54761121433375e-5,
             8.63351963521471e-7,
             2.47905831828260e-9,
         },
         {
             9.57910843855530e-3,
             2.02831789389435e-2,
             2.44990265668761e-2,
             1.93118044788806e-2,
             9.44753338706181e-3,
             2.62673065355268e-3,
             3.69805942462381e-4,
             2.22473818155004e-5,
             4.22360660759302e-7,
             1.21278082822666e-9,
         },
         {
             2.66331292395503e-3,
             5.63940297298984e-3,
             6.81154979071513e-3,
             5.36932834442915e-3,
             2.62673065355268e-3,
             7.30319083684035e-4,
             1.02818435790040e-4,
             6.18551714843314e-6,
             1.17430407389781e-7,
             3.37193682946561e-10,
         },
         {
             3.74956200622751e-4,
             7.93946927345222e-4,
             9.58968361136659e-4,
             7.55924299324728e-4,
             3.69805942462381e-4,
             1.02818435790040e-4,
             1.44753587500177e-5,
             8.70831958198587e-7,
             1.65325144471073e-8,
             4.74720266981928e-11,
         },
         {
             2.25572193453660e-5,
             4.77635386718788e-5,
             5.76911640119526e-5,
             4.54761121433375e-5,
             2.22473818155004e-5,
             6.18551714843314e-6,
             8.70831958198587e-7,
             5.23889122553910e-8,
             9.94589645655809e-10,
             2.85589868155718e-12,
         },
         {
             4.28242844331612e-7,
             9.06778151287977e-7,
             1.09525149314805e-6,
             8.63351963521471e-7,
             4.22360660759302e-7,
             1.17430407389781e-7,
             1.65325144471073e-8,
             9.94589645655809e-10,
             1.88820214174986e-11,
             5.42184812670272e-14,
         },
         {
             1.22967113105880e-9,
             2.60375375717955e-9,
             3.14494254883637e-9,
             2.47905831828260e-9,
             1.21278082822666e-9,
             3.37193682946561e-10,
             4.74720266981928e-11,
             2.85589868155718e-12,
             5.42184812670272e-14,
             1.55684799095647e-16,
         },
    };                  /* Array of weights precomputed for 2D integration */

    sum = 0.0;

    for (i = 0; i < 10; ++i)
    {
        for (j = 0; j < 10; ++j)
        {
            double sumand;

            /*
             * Since Phi is called with all NEGATIVE arguments, the exponent 
             * is negative (check that each term is negative). Thus it is 
             * impossible that the exponent can cause an overflow problem. 
             * The worst possible case is that the exponent is -Inf, in which 
             * case exp(-Inf)=0.
             */

            sumand = ww[i][j] * exp (a1 * (x[i]+x[i]-a1) +
                                     b1 * (x[j]+x[j]-b1) +
                                     2.0 * rho * (x[i]-a1) * (x[j]-b1));
            sum += sumand;
        }
    }

    BNC = sqrt (1.0 - rho * rho) / PI * sum;

    return BNC;
}


/* --- --- --- ---  From CCM - from ALIB ...   --- --- --- --- --- ---  --
function: NormalCumInvFast

Created by: David Hait

Description: This function calculates the value of Z given the area
under a cumulative normal distribution.  It uses a Rational approximation due
to Boris Moro.  The article "The Full Monte"  may be found in RISK magazine, 
Feb 1995.

--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  -- */
static double a[4] = {
     2.50662823884,
   -18.61500062529,
    41.39119773534,
   -25.44106049637
   };
   
static double b[4] = {
    -8.47351093090,
    23.08336743743,
   -21.06224101826,
    3.13082909833
   };

static double c[9] = { 
   0.3374754822726147,
   0.9761690190917186,
   0.1607979714918209,
   0.0276438810333863,
   0.0038405729373609,
   0.0003951896511919,
   0.0000321767881768,
   0.0000002888167364,
   0.0000003960315187
  };


static double fastApprox(double prob)
{
   double t,x,r;
   
   t = (prob < 0.5)?(1.0-prob):prob;
   
   x = t-0.5;
   if (fabs(x) < 0.42)
   {
     r=x*x;
     r=x*(((a[3]*r+a[2])*r+a[1])*r+a[0]) /
         ((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);
     return (prob < 0.5)? -r : r;
   }
   else
   {
     r=t;
     if (x>0.0) r = 1.0-t;
     r = log(-log(r));
     r = c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+
              r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));
     if (x<0.0) r = -r;
     return (prob < 0.5)? -r : r;
   }
}

/* Define constants for NormalCumInv  */
#define NORMAL_CUM_INV_RELATIVE_ACCURACY (1.0E-14)
//// this is less than DBL_MAX but does allow the resulting number to be fed
//// back into N1()
#define NORMAL_CUM_INV_DOUBLE_MAX (9.9E99)

/*--------------------------------------------------------------
// NormalCumInverse [0,1] -> R
// this function returns the cumulative normal inverse
// Moro's algorithm
*/
#define MAX_ITERATIONS         (100)
/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --
function: NormalCumInverse

created by: 09/29/92 Krishna Varikooty

description: This function calculates the value of Z given the area
under a cumulative normal distribution.

Modified by: 9/19/97 David Hait
description: Better polynomial for initial guess, accuracy improved to
  1.0E-14 for Newton (to take advantage of better NormalCum).

Modified:    10/3/97 Alexander Ng
Description: The Newton's method has been modified to exit after 
             MAX_ITERATIONS = 100 times and use the resulting answer
             as a best guess for the truth.  If MAX_ITERATIONS is reached,
             a message is written to the error log, but the answer is
             returned.  This was done to correct some pathologies in the
             rootfinding.
             
             It would be highly desirable to add a feature in the
             Newton's method that exits the while loop with status
             FAILURE when MAX_ITERATIONS is reached  AND the absolute
             error exceeds known error bounds.  The known error bounds are
             explicitly stated in Boris Moro's RISK magazine article.
--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  -- */
//// Probably should rename the existing N1Inverse N1InverseFast
double N1InverseBetter(double prob) /* (I) probability */
{
   double prob1;
   double prob_t, prob_d, dif, tzvalue;
   long   idx=0;       /* Number of iterations of Newton's method */

   prob1 = 1.0 - prob;
   if (prob1 <= 0)
   {
      return NORMAL_CUM_INV_DOUBLE_MAX;
   }

   if (prob <= 0.0)
   {
      return -NORMAL_CUM_INV_DOUBLE_MAX;
   }
   /*  use polynomial approximation to get the initial guess  */
   
   tzvalue = fastApprox(prob);
   
   /* use Newton-Raphson to get the final result */
   prob_t = N1(tzvalue);
   prob_d = N1Density(tzvalue);
   dif = prob_t-prob;

   while(fabs(dif) > NORMAL_CUM_INV_RELATIVE_ACCURACY*prob &&
         idx <= MAX_ITERATIONS   )
   {
     idx++;
     tzvalue = tzvalue-dif/prob_d;
     prob_t = N1(tzvalue);
     prob_d = N1Density(tzvalue);
     dif=prob_t-prob;
   }

   return tzvalue;
}

/** Calculates the number of jumps n such that 
    the probability of more than n jumps in time interval tau is less than epsilon. */
int Quantile(double       tau,
            double         jumpRate,
            double         epsilon,
            int            maxJumps,
            double*        quantile)
{
    if (jumpRate <= 0.0 || Maths::isZero(tau))
        return 0;

    double proba    = exp(-tau* jumpRate);
    *quantile = 0;
    int    iStep;

    for( iStep=1; iStep<=maxJumps; iStep++ )
    {
        *quantile += proba;
        proba *= jumpRate * tau / (double)iStep;
        if( *quantile>1.0-epsilon ) break;
    }

    return iStep;
}

/*f----------------------------------------------------------------------------
 * Cumulative error function weighted by exponential factor. Computes
 * exp(b+a^2/2)erfc(x-a)
 *
 * Alib comment:  The routine has a relative accuracy no worse than 1.0E-14, 
 * where relative accuracy is defined as (computed - truth)/truth, and truth
 * comes from a continued fraction calculation.  This is essentially 
 * accurate to the next to last decimal digit of machine accuracy on the Sun.
 */
double ExpCErrFcn (double a, double b, double x) {

    const double SQRT2   = 1.414213562373095049;     /* sqrt(2) */
    const double SQRTPI  = 1.772453850905516027;     /* sqrt(pi) */

    /* Coefficients in expression of erf(x) for -0.46875<=x<=0.46875 */
    const double P10 = 3209.377589138469472562;    /* Numerator */
    const double P11 = 377.4852376853020208137;
    const double P12 = 113.8641541510501556495;
    const double P13 = 3.161123743870565596947;
    const double P14 = 0.1857777061846031526730;
    const double Q10 = 2844.236833439170622273;   /* Denominator */
    const double Q11 = 1282.616526077372275645;
    const double Q12 = 244.0246379344441733056;
    const double Q13 = 23.60129095234412093499;
    const double Q14 = 1.0;

    /* Coefficients in expression of erfc(x) for 0.46875<=x<=4.0 */
    const double P20 = 1230.33935479799725272;  /* Numerator */
    const double P21 = 2051.07837782607146532;
    const double P22 = 1712.04761263407058314;
    const double P23 = 881.952221241769090411;
    const double P24 = 298.635138197400131132;
    const double P25 = 66.1191906371416294775;
    const double P26 = 8.88314979438837594118;
    const double P27 = 0.564188496988670089180;
    const double P28 = 2.15311535474403846343e-8;
    const double Q20 = 1230.33935480374942043;  /* Denominator */
    const double Q21 = 3439.36767414372163696;
    const double Q22 = 4362.61909014324715820;
    const double Q23 = 3290.79923573345962678;
    const double Q24 = 1621.38957456669018874;
    const double Q25 = 537.181101862009857509;
    const double Q26 = 117.693950891312499305;
    const double Q27 = 15.7449261107098347253;
    const double Q28 = 1.0;

    /* Coefficients in expression of erfc(x) for x>= 4.0 */
    const double P30 = -6.58749161529837803157E-4;    /* Numerator */
    const double P31 = -1.60837851487422766278E-2;
    const double P32 = -1.25781726111229246204E-1;
    const double P33 = -3.60344899949804439429E-1;
    const double P34 = -3.05326634961232344035E-1;
    const double P35 = -1.63153871373020978498E-2;
    const double Q30 =  2.33520497626869185443E-3 ;   /* Denominator */
    const double Q31 =  6.05183413124413191178E-2;
    const double Q32 =  5.27905102951428412248E-1;
    const double Q33 =  1.87295284992346047209;
    const double Q34 =  2.56852019228982242072;
    const double Q35 =  1.0;


    double numerator;            /* numerator of polynomial in expression */
    double denominator;          /* denominator of polynomial in expression */
    double y,y2;                 /* y = abs(x)/sqrt(2), y2 = y*y */
    double erfc;                 /* return value */

    double W = b + 0.5 * a * a;

    y  = fabs(x-a) / SQRT2;
    y2 = y * y;

    if (y < 0.46875) 
    {
        numerator   = P10 + y2*(P11 + y2*(P12 + y2*(P13 +y2*P14)));
        denominator = Q10 + y2*(Q11 + y2*(Q12 + y2*(Q13 +y2*Q14)));
        erfc = exp(W) * (1 - y * numerator / denominator);
    }
    else if (y < 4.0) 
    {
        numerator   = P20 + y*(P21 + y*(P22 + y*(P23 +
                                                 y*(P24 + y*(P25 + y*(P26 + y*(P27 + y*P28)))))));
        denominator = Q20 + y*(Q21 + y*(Q22 + y*(Q23 +
                                                 y*(Q24 + y*(Q25 + y*(Q26 + y*(Q27 + y*Q28)))))));
        erfc = exp(-y2 + W) * numerator / denominator;
    }
    else /* (y > 4.0) */ 
    {
        double z2 = 1/y2; 
        numerator   = P30 + z2*(P31 + z2*(P32 + z2*(P33 + z2*(P34 +z2*P35))));
        denominator = Q30 + z2*(Q31 + z2*(Q32 + z2*(Q33 + z2*(Q34 +z2*Q35))));
        erfc = (exp(-y2 + W)/y) * 
            (1.0 / SQRTPI + numerator / (denominator * y2)); 
    }
    
    return erfc;
} /*ExpCerrFcn */



/**Number of permutations of k from n*/
long Permutations(int n, int k) {
    long p = 1;
    if (n<=0) throw ModelException("mathlib::Permutations","n must be strictly greater than zero!");
    if (k<0) throw ModelException("mathlib::Permutations","k must be greater than zero!");
    if (k>n) throw ModelException("mathlib::Permutations","k must be less than or equal to n!");

    // P(n,k) = n (n-1) ... (n-k+1) = n!/(n-k)!
    k = n-k;
    while (n>k) {
        p *= n;
        n--;
    }
    return p;
}

/**Number of combinations of k things from n - binomial coefficient.*/
long Combinations(int n, int k) {
    long p = Permutations(n,k);
    while (k>0) {
        p/=k;
        k--;
    }
    return p;
}

/**Probability of k out of N independent Bernoulli trials with probality p. Returns 0 if (N,p,k) don't make
   sense: (N,p,k) must satisfy N>0, 0<=p<=1, 0<=k<=N */
double BinomialProbability(int N, double p, int k) {
    if (N<=0 || k<0 || k>N || p<0 || p>1) return 0;
    return pow(p,k) * pow((1-p),N-k) * exp(gammln(N+1.0) - gammln(N-k+1.0) - gammln(k+1.0));
}

/**Sets given array to probabilities of outcomes of N trials. resultArray always ends up with N+1 elements,
   and they will always sum to 1.0 exactly: if necessary, the largest element will be adjusted to
   ensure this. If the cumulative probability is > truncationBound, then all remaining probabilities
   are set to zero. */
void BinomialProbabilities(int N, double p, DoubleArray* resultArray, double truncationBound) {

    if (!resultArray || N<=0) return;
    
    if ((N+1)!=resultArray->size()) {
        resultArray->resize(N+1);
    }
    if (p<0 || p>1 || truncationBound<0) return;

    DoubleArray& ary = *resultArray; // to make syntax prettier
    double oneMinusP = 1.0-p;
    double sum = pow(oneMinusP, N);
    int biggestK = 0;
    double biggest = sum;
    ary[0] = sum;
    int k=1;
    while (k<=N) {
        double e = (sum<truncationBound ? ary[k-1]*(N-k)*p/((k+1) * oneMinusP) : 0.0);
        ary[k] = e;
        if (e>biggest) {
            biggestK = k;
            biggest = e;
        }
        sum += e;
        k++;
    }
    ary[biggestK] += 1.0 - sum;
    return;
}

// computes raw material. Called by the constructor ---------------------------
void CoeffsPQ4Pade::resizeSeriesFCoeffs(int gMaxOrder) const
{
    const int _initOrder = seriesF.size();
    if ((gMaxOrder + 1) < _initOrder){
        return;
    }
    const static string function = "CoeffsPQ4Pade::resizeSeriesFCoeffs : ";
    QLIB_VERIFY(gMaxOrder >= 1, function + "maxOrder should be >= 1");
    seriesF.resize(gMaxOrder+1);
    for(int m = 1; m <= gMaxOrder; m++){
        seriesF[m] = seriesF[m-1]/((double) m);
    }
}
// public constructor ---------------------------------------------------------
CoeffsPQ4Pade::CoeffsPQ4Pade(int gMaxOrder)
    : seriesF(1, 1.0) 
{
    this->resizeSeriesFCoeffs(gMaxOrder);
}

// clone constructor ----------------------------------------------------------
CoeffsPQ4Pade::CoeffsPQ4Pade(const CoeffsPQ4Pade& g2Clone)
{
    seriesF = g2Clone.seriesF;
    seriesPQ= g2Clone.seriesPQ;
}

// core computational method --------------------------------------------------
pairDoubleVectorPointers CoeffsPQ4Pade::calsPQ(int gOrderP, int gOrderQ) const
{
    static const string method = "CoeffsPQ4Pade::calsPQ";
    ///////////////////////////////////////////////////////////////////////////
    // check inputs -----------------------------------------------------------
    QLIB_VERIFY((gOrderP >= 0) && (gOrderQ >= 0),
        method + 
        "The approximation order coeffs orderP and OrderQ must both be >= 0.");
    ///////////////////////////////////////////////////////////////////////////
    // check if we already have PQ series -------------------------------------

    // if we do not have an iterator for P, insert a new element  -------------
    map<int, map<int, pair<DoubleVector, DoubleVector> > >::iterator 
                                                 _itP = seriesPQ.find(gOrderP);
    if (_itP == seriesPQ.end())
    {
        // create -------------------------------------------------------------
        seriesPQ.insert(pair<int, map<int, pair<DoubleVector,DoubleVector> > >(
                        gOrderP,map<int, pair<DoubleVector, DoubleVector> >()));
        _itP = seriesPQ.find(gOrderP);
    }
    // if we do not have an iterator for P, insert a new element  -------------
    map<int, pair<DoubleVector, DoubleVector> >::iterator 
                                            _itPQ = _itP->second.find(gOrderQ);
    if (_itPQ != _itP->second.end())
    {
        // return the result straight away ------------------------------------
        return pairDoubleVectorPointers(&(_itPQ->second.first), 
                                        &(_itPQ->second.second));
    }
    ///////////////////////////////////////////////////////////////////////////
    // if we are here, we need to create a new pair of coeffs vectors --------
    _itP->second.insert(pair<int, pair<DoubleVector,DoubleVector> >(
                                  gOrderQ, pair<DoubleVector,DoubleVector>()));
    _itPQ = _itP->second.find(gOrderQ);
    DoubleVector & _resultP = _itPQ->second.first;
    DoubleVector & _resultQ = _itPQ->second.second;
    // ------------------------------------------------------------------------
    // resize seriesF if needed -----------------------------------------------
    this->resizeSeriesFCoeffs(gOrderP + gOrderQ + 1);
    // ------------------------------------------------------------------------
    // do _resultQ ------------------------------------------------------------
    _resultQ.resize(gOrderQ+1, 0.0);
    _resultQ[0] = 1.0; // arbitrary normalization

    if (gOrderQ==1) {
        /* DON'T BOTHER TO DO MATRIX INVERSION FOR Q ORDER 1 */
        _resultQ[1] = -seriesF[gOrderP+1]/seriesF[gOrderP];

    } else { // gOrderQ==1
        /* DETERMINE Q COEFFICIENTS BY MATCHING COEFFICIENTS IN P(x)= Q(x).f(x)
            FROM ORDER  m+1 TO m+n+1: SUM(j=0 to j=k) q[j].f[k-j] = 0
            FOR k=m+1, ..., m+n+1. THIS IS EQUIVALENT TO SOLVING Aq = -*/
        DoubleMatrix A(gOrderQ, gOrderQ); // initialize with zeros
        for (int c=0; c<gOrderQ; c++) {
            for (int r=0; r<gOrderQ; r++) {
                int fidx = gOrderP-c+r;
                // we will scale A to escape numerical issues
                A[c][r] = (fidx<0) ? 0.0 : seriesF[fidx]/seriesF[gOrderP];
            }
        }
        const DoubleMatrix & inverseA = A.computeInverse();
        for (int r=0; r<gOrderQ; r++) {
            _resultQ[r+1] = 0.0;
            for(int s=0; s < gOrderQ; s++) {
                _resultQ[r+1] 
                  -= seriesF[gOrderP+s+1] * inverseA[s][r] / seriesF[gOrderP]; 
            }
        }
    }
    // ------------------------------------------------------------------------
    // do _resultP ------------------------------------------------------------
    /* NOW SET P COEFFICIENTS: p[k] =  SUM(j=0 to k) q[j].f[k-j] for k=0, ..., m */
    _resultP.resize(gOrderP+1);
    for (int k=0; k<=gOrderP; k++) {
        _resultP[k] = 0.0;
        for (int j=0; j<=k; j++) {
            _resultP[k] += _resultQ[j]*seriesF[k-j];
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    return pairDoubleVectorPointers(&_resultP, &_resultQ);
}

 /*=============================================================================
 * Class to provide add-in functions
 *===========================================================================*/
FORWARD_DECLARE(PadeApproximation)
class PadeApproximation : public CObject {
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(PadeApproximation, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPadeApproximation);
        FIELD(fCoefficients, "Series expansion coefficients of the function you want to approximate.");
        FIELD(numeratorOrder, "Order of polynomial P in f ~ P/Q.");
        FIELD(denominatorOrder, "Order of polynomial Q in f ~ P/Q.");
        FIELD(constructed,"Internal flag");
        FIELD(pCoefficients,"Internal array");
        FIELD(qCoefficients,"Internal array");
        FIELD_MAKE_TRANSIENT(constructed);
        FIELD_MAKE_TRANSIENT(pCoefficients);
        FIELD_MAKE_TRANSIENT(qCoefficients);
    }

    PadeApproximation(const DoubleArray& fSeries, int m, int n) : CObject(TYPE), fCoefficients(fSeries),
        numeratorOrder(m), denominatorOrder(m),
        constructed(false), pCoefficients(), qCoefficients() {};

    double operator()(double x) {
        static const char* method = "PadeApproximation::operator()";
        if (!constructed) calculateCoefficients();

        double qVal = 0.0;
        for (int i=qCoefficients.size()-1; i>=0; i--) {
            qVal = qVal*x + qCoefficients[i];
        }
        if (qVal==0.0) throw ModelException(method,"Zero denominator in Pade approximation.");

        double pVal = 0.0;
        for (int i=pCoefficients.size()-1; i>=0; i--) {
            pVal = pVal*x + pCoefficients[i];
        }
        
        return pVal/qVal;

    };

    double getNumeratorCoefficient(int i) {
        if (!constructed) calculateCoefficients();
        if (i>=0 && i<pCoefficients.size()) {
            return pCoefficients[i];
        } else {
            return 0;
        }
    }

    double getDenominatorCoefficient(int i) {
        if (!constructed) calculateCoefficients();
        if (i>=0 && i<qCoefficients.size()) {
            return qCoefficients[i];
        } else {
            return 0;
        }
    }

private:
    DoubleArray fCoefficients;
    int numeratorOrder;
    int denominatorOrder;
    bool constructed;
    DoubleArray pCoefficients;
    DoubleArray qCoefficients;

    void calculateCoefficients() {
        if (constructed) return;

        //ComputePadeCoefficients(fCoefficients, 
        //    numeratorOrder, denominatorOrder, &pCoefficients, &qCoefficients);
        
        constructed = true;
    }

    PadeApproximation(): CObject(TYPE), fCoefficients(), numeratorOrder(0), denominatorOrder(0),
        constructed(false), pCoefficients(), qCoefficients() {};  
    

    static IObject* defaultPadeApproximation() {
        return new PadeApproximation();
    }
};

CClassConstSP const PadeApproximation::TYPE = CClass::registerClassLoadMethod(
    "PadeApproximation", typeid(PadeApproximation), PadeApproximation::load);


    
/*=============================================================================
 * Classes to provide add-in functions
 *===========================================================================*/
class PadeValueAddin : public CObject {
public:
    static CClassConstSP const TYPE;

    PadeApproximationSP padeObject;
    double x;

    PadeValueAddin() : CObject(TYPE), padeObject(), x(0) {};
 
    double getValue() {
       return (*padeObject)(x);
    }
    
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(PadeValueAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPadeValueAddin);
        FIELD(padeObject, "Object you wish to apply function to.");
        FIELD(x, "Value at which to compute Pade approximation P(x)/Q(x).");

        Addin::registerDoubleMethod("PADE_GET_VALUE", 
            Addin::UTILITIES, 
            "Returns the value of a Pade approximation at a point x.", 
            &PadeValueAddin::getValue);

    }

    static IObject* defaultPadeValueAddin() {
        return new PadeValueAddin();
    }
};

CClassConstSP const PadeValueAddin::TYPE = CClass::registerClassLoadMethod(
    "PadeValueAddin", typeid(PadeValueAddin), PadeValueAddin::load);


class PadeCoefficientAddin : public CObject {
public:
    static CClassConstSP const TYPE;

    PadeApproximationSP padeObject;
    bool isNumeratorCoefficient;
    int i;

    PadeCoefficientAddin() : CObject(TYPE), padeObject(), isNumeratorCoefficient(false), i(0) {};
 
    double getValue() {
        return (isNumeratorCoefficient 
            ? padeObject->getNumeratorCoefficient(i) 
            : padeObject->getDenominatorCoefficient(i));
    }
    
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(PadeCoefficientAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPadeCoefficientAddin);
        FIELD(padeObject, "Object you wish to apply function to.");
        FIELD(isNumeratorCoefficient, "If true, return numerator polynomial coefficient, else denominator.");
        FIELD(i, "Order of coefficient to return.");

        Addin::registerDoubleMethod("PADE_GET_COEFFICIENT", 
            Addin::UTILITIES, 
            "Returns the value of the Pade approximation coefficient.", 
            &PadeCoefficientAddin::getValue);

    }

    static IObject* defaultPadeCoefficientAddin() {
        return new PadeCoefficientAddin();
    }
};

CClassConstSP const PadeCoefficientAddin::TYPE = CClass::registerClassLoadMethod(
    "PadeCoefficientAddin", typeid(PadeCoefficientAddin), PadeCoefficientAddin::load);



DRLIB_END_NAMESPACE
