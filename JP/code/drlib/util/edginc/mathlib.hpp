//----------------------------------------------------------------------------
//
//     File           : Mathlib.hpp
//     Author         : Ning Shen
//
//     Description    : maths support functions.
//                     This file does not require any external functions.
//
//
//----------------------------------------------------------------------------

#ifndef _MATHLIB_H_
#define _MATHLIB_H_

/*
*============================================================================
*
* Defines
*
*============================================================================
*/

#include <math.h>
#include <map>
#include "edginc/AtomicArray.hpp"

#ifdef __unix__
#define _DOMAIN DOMAIN
#define _SING SING
#define _OVERFLOW OVERFLOW
#define _UNDERFLOW UNDERFLOW
#define _TLOSS TLOSS
#define _PLOSS PLOSS
#endif // __unix__

DRLIB_BEGIN_NAMESPACE

/*============================================================================
*
* maths function prototypes
*
*============================================================================*/

UTIL_DLL double N1Density(double StandardVariate);
/*
 * utility structure used to compute gaussian weights for discrete integration
 * with constant step
 */
struct GaussianIntegrationMethod{
    double l;    /* lower bound */
    double u;    /* upper bound */
    int    n;    /* nb of steps */
    double step; /* step size   */
    UTIL_DLL GaussianIntegrationMethod(double lowerBound, double upperBound, long nbStep);
    UTIL_DLL double weightCalc(double m);
};

UTIL_DLL double N1(double StandardVariate);
UTIL_DLL double N1Inverse(double Probability);
//// A more accurate version of N1Inverse
UTIL_DLL double N1InverseBetter(double prob); /* (I) probability */
UTIL_DLL double N2Density(double StandardVariate1, double StandardVariate2, double Correlation);
UTIL_DLL double N2(double StandardVariate1, double StandardVariate2, double Correlation, bool highAccuracy = true);
inline UTIL_DLL double N2Std(double StandardVariate1, double StandardVariate2, double Correlation) { return N2(StandardVariate1, StandardVariate2, Correlation, false); }
UTIL_DLL double LNDensity(double  x, double mat, double vol, double fwd);
UTIL_DLL int SolveQuadratic(double a, double b, double c, double *x1, double *x2);
UTIL_DLL int SolveCubic(double a, double b, double c, double *x1, double *x2, double *x3);
UTIL_DLL double QuadraticInterp(double x, double x1, double x2, double x3, double y1, double y2, double y3);
UTIL_DLL double LinearInterp(double x, double x1, double x2, double y1, double y2);
UTIL_DLL double interpF2(double x, double y, vector<double>& vX, vector<double>& vY, double** fxy);
UTIL_DLL double interpF2(double x, double y, double* vX, double* vY, double* fxy, int dim_x, int dim_y);
UTIL_DLL double interpF2_linear(double x, double y, double* vX, double* vY, double* fxy, int dim_x, int dim_y);
UTIL_DLL double interpF2_linear(double x, double y, vector<double>& vX, vector<double>& vY, double** fxy);
UTIL_DLL double CubicInterp(double x, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
UTIL_DLL double atanh(double x);
UTIL_DLL int Quantile(double  tau, double jumpRate, double epsilon, int  maxJumps, double*  quantile);
/**Number of combinations of k things from n - binomial coefficient.*/
UTIL_DLL long Combinations(int n, int k);
/**Number of permutations of k from n*/
UTIL_DLL long Permutations(int n, int k);
/**Probability of k out of N independent Bernoulli trials with probality p. Returns 0 if (N,p,k) don't make
   sense: (N,p,k) must satisfy N>0, 0<=p<=1, 0<=k<=N */
double BinomialProbability(int N, double p, int k);
/**Sets given array to probabilities of outcomes of N trials. resultArray always ends up with N+1 elements,
   and they will always sum to 1.0 exactly. If necessary, the largest element will be adjusted to
   ensure this. If the cumulative probability is > truncationBound, then all remaining probabilities
   are set to zero. */
void BinomialProbabilities(int N, double p, DoubleArray* resultArray, double truncationBound);

/** Cumulative error function weighted by exponential factor. Computes
 * exp(b+a^2/2)erfc(x-a)
 * Lifted from ALIB via IR Q3 library (charles.morcom)
 * Alib comment:  The routine has a relative accuracy no worse than 1.0E-14, 
 * where relative accuracy is defined as (computed - truth)/truth, and truth
 * comes from a continued fraction calculation.  This is essentially 
 * accurate to the next to last decimal digit of machine accuracy on the Sun.
 */
UTIL_DLL double ExpCErrFcn (double a, double b, double x);

/**Computes the (m,n) order Pade approximation to a function, f, with given Taylor series
 * coefficients. At least m+n+1 coefficients of f's series expansion must be provided.
 * The Pade approximation is a rational function P(x)/Q(x), where polynomials P and Q have order
 * m and n respectively. The coefficients of P and Q are determined by matching f(0) and the
 * derivatives D^kf(0) for k = 1,...m+n. The zero-order coefficient of Q is set to 1. */
typedef vector<double>  DoubleVector;
typedef vector<int>     IntVector;
typedef vector<double*> DoublePointerVector;
                        
typedef pair<const DoubleVector *, const DoubleVector *> pairDoubleVectorPointers;
class CoeffsPQ4Pade
{
public:
    // public constructor -----------------------------------------------------
    CoeffsPQ4Pade(int maxOrder); // initial max order does not matter
                        // if it is not sufficient seriesF will be resized
    // core computation method ------------------------------------------------
    pairDoubleVectorPointers calsPQ(int orderP, int orderQ) const;
    // aux
    CoeffsPQ4Pade(const CoeffsPQ4Pade& g2Clone);
private:
    void operator=(const CoeffsPQ4Pade& g2Clone);// not defined do not use
    CoeffsPQ4Pade();                             // not defined do not use
    // computes raw material. Called by the constructor -----------------------
    void resizeSeriesFCoeffs(int gMaxOrder) const;
    // cashed raw material ----------------------------------------------------
    mutable DoubleVector seriesF;
    // we cash the results here -----------------------------------------------
    mutable map<int, map<int, pair<DoubleVector, DoubleVector> > > seriesPQ;
};

//void ComputePadeCoefficients(
//    const DoubleArray&  fCoefficients, /**<series expansion coefficients of taget function       */
//    int                 m,             /**<Order of P                                            */
//    int                 n,             /**<Order of Q                                            */
//    DoubleArray*        pCoefficients, /**<m+1 coefficients of P, if successful                  */
//    DoubleArray*        qCoefficients  /**<n+1 coefficients of Q, if successful.                 */
//    );

DRLIB_END_NAMESPACE

#endif
