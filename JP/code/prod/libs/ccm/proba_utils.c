#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "error2.h"
#include "beta.h"
#include "proba_utils.h"
#include "rootbrent.h"

/* Define constants for NormalCumInv  */
#define NORMAL_CUM_INV_RELATIVE_ACCURACY (1.0E-14)
#define IS_ALMOST_ZERO(x) ((x)<1.0E-3)?1:0
#define LOTUS_DBL_MAX 9.9E99

/*---------------------------------------------------------------
// NormalCum R -> [0,1]
// this function returns the cumulative normal
*/
/*------------------------------------------------------------------------
  Routine: NormalCum

  Coded by: BBroder/NShewmaker  10/14/96

  Uses approximations found in:
      Rational Checbyshev Approximations for the Error Function, by W.J. Cody,
      Mathematics of Computation, published by the AMS, Vol 23, No 107,
      July, 1969, pp.631-637.

  The routine has a relative accuracy no worse than 1.0E-14, where relative
  accuracy is defined as (computed - truth)/truth, and truth comes from
  a continued fraction calculation.  This is essentially accurate to the
  next to last decimal digit of machine accuracy on the Sun.

 *------------------------------------------------------------------------ */
#define SQRT2   1.414213562373095049     /* sqrt(2) */
#define SQRTPI  1.772453850905516027     /* sqrt(pi) */

/* Coefficients in expression of erf(x) for -0.46875<=x<=0.46875 */
#define P10 3209.377589138469472562    /* Numerator */
#define P11 377.4852376853020208137
#define P12 113.8641541510501556495
#define P13 3.161123743870565596947
#define P14 0.1857777061846031526730
#define Q10 2844.236833439170622273   /* Denominator */
#define Q11 1282.616526077372275645
#define Q12 244.0246379344441733056
#define Q13 23.60129095234412093499
#define Q14 1.0

/* Coefficients in expression of erfc(x) for 0.46875<=x<=4.0 */
#define P20 1230.33935479799725272  /* Numerator */
#define P21 2051.07837782607146532
#define P22 1712.04761263407058314
#define P23 881.952221241769090411
#define P24 298.635138197400131132
#define P25 66.1191906371416294775
#define P26 8.88314979438837594118
#define P27 0.564188496988670089180
#define P28 2.15311535474403846343e-8
#define Q20 1230.33935480374942043  /* Denominator */
#define Q21 3439.36767414372163696
#define Q22 4362.61909014324715820
#define Q23 3290.79923573345962678
#define Q24 1621.38957456669018874
#define Q25 537.181101862009857509
#define Q26 117.693950891312499305
#define Q27 15.7449261107098347253
#define Q28 1.0

/* Coefficients in expression of erfc(x) for x>= 4.0 */
#define P30 -6.58749161529837803157E-4    /* Numerator */
#define P31 -1.60837851487422766278E-2
#define P32 -1.25781726111229246204E-1
#define P33 -3.60344899949804439429E-1
#define P34 -3.05326634961232344035E-1
#define P35 -1.63153871373020978498E-2
#define Q30  2.33520497626869185443E-3    /* Denominator */
#define Q31  6.05183413124413191178E-2
#define Q32  5.27905102951428412248E-1
#define Q33  1.87295284992346047209
#define Q34  2.56852019228982242072
#define Q35  1.0

double NormalCum(double x)
{
   double numerator;            /* Numerator of polynomial in expression */
   double denominator;          /* Denominator of polynomial in expression */
   double y;                    /* y = abs(x)/sqrt(2) */
   double y2;                   /* y*y */
   double erf;                  /* Error function value */
   double erfc;                 /* Complimentary Error function value */

   y  = fabs(x) / SQRT2;
   y2 = y * y;

   if (y < 0.46875)
   {
      numerator   = P10 + y2*(P11 + y2*(P12 + y2*(P13 +y2*P14)));
      denominator = Q10 + y2*(Q11 + y2*(Q12 + y2*(Q13 +y2*Q14)));
      erf = y * numerator / denominator;
      return (x>0.0) ? 0.5 + 0.5*erf : 0.5 - 0.5*erf;
   }
   else if (y < 4.0)
   {
      numerator   = P20 + y*(P21 + y*(P22 + y*(P23 +
                          y*(P24 + y*(P25 + y*(P26 + y*(P27 + y*P28)))))));
      denominator = Q20 + y*(Q21 + y*(Q22 + y*(Q23 +
                          y*(Q24 + y*(Q25 + y*(Q26 + y*(Q27 + y*Q28)))))));
      erfc = exp(-y2) * numerator / denominator;
      return (x>0.0) ? 1.0 - 0.5*erfc : 0.5*erfc;
   }
   else /* (y > 4.0) */
   {
      double z2 = 1/y2;
      numerator   = P30 + z2*(P31 + z2*(P32 + z2*(P33 + z2*(P34 +z2*P35))));
      denominator = Q30 + z2*(Q31 + z2*(Q32 + z2*(Q33 + z2*(Q34 +z2*Q35))));
      erfc = (exp(-y2)/y) * (1.0 / SQRTPI + numerator / (denominator * y2));
      return (x>0.0) ? 1.0 - 0.5*erfc : 0.5*erfc;
   }
}

/*---------------------------------------------------------------
// NormalDensity R -> R^+
// this function returns normal density 
*/
/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --
function: NormalDensity

created by: 09/29/92 Krishna Varikooty

description: This program calculates the density of a normal distribution
function

inputs:

outputs:

notes on use: If daysToPayment is zero, we get the values at expiry

modification history
--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  -- */
#define ONE_DIV_SQRT_TWO_PI  0.3989422804014327

double NormalDensity(double x)
{
   return ONE_DIV_SQRT_TWO_PI * exp (- x*x * 0.5);
}

/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --
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


double fastApprox(double prob)
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

double NormalCumInvFast(double prob) /* (I) probability */
{
   double prob1;

   prob1 = 1.0 - prob;
   if (IS_ALMOST_ZERO(prob1) || prob1 < 0)
   {
      return LOTUS_DBL_MAX;
   }

   if (IS_ALMOST_ZERO(prob) || prob < 0)
   {
      return -LOTUS_DBL_MAX;
   }
   /*  use polynomial approximation to get the initial guess  */
   return fastApprox(prob);
}

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

double NormalCumInverse(double prob) /* (I) probability */
{
   double prob1;
   double prob_t, prob_d, dif, tzvalue;
   long   idx=0;       /* Number of iterations of Newton's method */

   prob1 = 1.0 - prob;
   if (prob1 <= 0)
   {
      return LOTUS_DBL_MAX;
   }

   if (prob <= 0.0)
   {
      return -LOTUS_DBL_MAX;
   }
   /*  use polynomial approximation to get the initial guess  */
   
   tzvalue = fastApprox(prob);
   
   /* use Newton-Raphson to get the final result */


   prob_t = NormalCum(tzvalue);
   prob_d = NormalDensity(tzvalue);
   dif = prob_t-prob;

   while(fabs(dif) > NORMAL_CUM_INV_RELATIVE_ACCURACY*prob &&
         idx <= MAX_ITERATIONS   )
   {
     idx++;
     tzvalue = tzvalue-dif/prob_d;
     prob_t = NormalCum(tzvalue);
     prob_d = NormalDensity(tzvalue);
     dif=prob_t-prob;
   }

   return tzvalue;
}

/*---------------------------------------------------------------
// Chi2Density R^+ -> R^+
// this function returns the chi2 density
*/
double Chi2Density(double x,
                   long freedomDegree)
{
    return pow(x,freedomDegree/2.0-1)*exp(-x/2.0 - gammln(freedomDegree/2.0))/
                                            pow(2,freedomDegree/2);
}

/*---------------------------------------------------------------
// Chi2Cum R^+ -> [0,1]
// this function returns the chi2 cumulative
*/
double Chi2Cum(double x,
               long freedomDegree)
{
    return gammp(freedomDegree *0.5, x *0.5);
}

/*--------------------------------------------------------------
// StudentCum R -> [0,1]
// this function returns the cumulative student
*/
double StudentCum(double x,
                  long freedomDegree)
{
    double b;
    if(x>3e100)
    {
        return 1.0;
    }
    else if(x< -3e100)
    {
        return 0.0;
    }

    b = .5 * betai(0.5 * freedomDegree,0.5,freedomDegree/(freedomDegree + x*x));
    return x>0 ? 1.-b : b;
}



/** root solving implementation */
typedef struct {double s; long f;} PairStruct;
static int StudentCumBrent(double x, void *data, double *out)
{
    PairStruct *p = (PairStruct*)data;
    *out = log(StudentCum(x, p->f)) - p->s;
    return 0;
}

/*--------------------------------------------------------------
// StudentCumInverse [0,1] -> R
// this function returns the cumulative student inverse
*/
double StudentCumInverse(double y,
                         long freedomDegree)
{
    double x;
    PairStruct p;
    if (y < 3e-16) return -7e15;
    if (1.-y < 1e-9) return 7e15;
    p.s = log(y);
    p.f = freedomDegree;
    RootFindBrent(&StudentCumBrent, &p, -1e300, 1e300, 
        50, 0., 1., 0., 1e-11, 1e-11, &x);
    return x;
}


/*--------------------------------------------------------------
// StableDensity R -> [0,1]
// this function returns the density stable
*/
extern void stable(int *n, double *y, double *beta, double *alpha, int *npt,
	    double *up, double *eps, int *type, int *err, double *ffy);
double StableDensity(double x, double alpha, double beta)
{
    int n=1, err=0;
    double eps=1e-12,ffy;
    double up = 100.0;
    int npt = 1000;
    int type = 1;

    stable(&n,&x,&beta,&alpha,&npt,&up,&eps,&type,&err,&ffy);
    return ffy;
}


/*--------------------------------------------------------------
// StableCum R -> [0,1]
// this function returns the cumulative stable
*/
extern void pstable(int *n, double *y, double *beta, double *alpha,
	     double *eps, int *err, double *ffy);
double StableCum(double x, double alpha, double beta)
{
    int n=1, err=0;
    double eps=1e-12,ffy;
    pstable(&n,&x,&beta,&alpha,&eps,&err,&ffy);
    return ffy;
}

/** root solving implementation */
typedef struct {double s,a,b;} TripletStruct;
static int StableCumBrent(double x, void *data, double *out)
{
    TripletStruct *p = (TripletStruct*)data;
    *out = log(StableCum(x, p->a,p->b)) - p->s;
    return 0;
}

/*--------------------------------------------------------------
// StableCumInverse R -> [0,1]
// this function returns the cumulative stable
*/
double StableCumInverse(double y, double alpha, double beta)
{
    double x;
    TripletStruct p;
    if (y < 1e-8) return -7e15;
    if (1.-y < 1e-8) return 7e15;
    p.s = log(y);
    p.a = alpha;
    p.b = beta;
    RootFindBrent(&StableCumBrent, &p, -1e300, 1e300, 
        50, 0., 1., 0., 1e-12, 1e-12, &x);
    return x;
}


/*----------------------------------------------------------------
// DensityBarChart
// this function returns the number of sample between two values
// defined by a start value, a step value and a number of steps
*/
void DensityBarChart(double *bar,
                     const double *sample,
                     long sampleSize,
                     double start,
                     double stepSize,
                     long nbSteps)
{
    double lowerBound = start;
    double upperBound = start + stepSize * nbSteps;
    double sample_j;
    int i,j;
    long usedSampleSize = 0;

    for(i=0;i<nbSteps;i++)
    {
        bar[i] = 0;
    }

    for(j=0;j<sampleSize;j++)
    {
        sample_j = sample[j];
        if(sample_j>=lowerBound && sample_j<upperBound)
        {
            bar[(long) floor((sample_j - start)/stepSize)] += 1;
            usedSampleSize += 1;
        }
    }

    for(i=0;i<nbSteps;i++)
    {
        bar[i] /= usedSampleSize * stepSize;
    }
 
}


/* ----------------------------------------------------------------
// Moments
*/
int Moments(const double *probas,
            const double *values,
            long nbValues,
            double *moments,
            long nbMoments)
{
    static char routine[] = "Moments";
    int i,j;
    int status = FAILURE;
    double *product = NULL;
    double m;
    product = malloc(nbValues*sizeof(double));
    if(product==NULL) goto RETURN;
    memcpy(product,probas,nbValues*sizeof(double));
    
    for(i=0;i<nbMoments;i++)
    {
        m = 0;
        for(j=0;j<nbValues;j++)
        {
            product[j] *= values[j];
            m += product[j];
        }
        moments[i] = m;
    }
    status = SUCCESS;
RETURN:
    
    if(product) free(product);
    if(status == FAILURE)
    {
         DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);   
    }
    return status;
}
