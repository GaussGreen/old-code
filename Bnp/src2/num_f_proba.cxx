// S.Galluccio: 21 November 2000
// Added functions to evaluate the cumulative Student distribution
// and the Inverse cumulative Student distribution with n degrees of freedom

/* ==========================================================================
   FILE_NAME:	num_f_proba.cxx

   PURPOSE:     useful probability distributions...
   ========================================================================== */

#include "utallhdr.h"
#include <math.h"
#include <num_h_allhdr.h"
#include <num_h_proba.h"
#include <num_h_simpso.h"

#define INFINITY 1000
#define LOWFLOOR -88

static int Integer_Part(double d) {
  int ret;

  ret = (int)d;
  if (ret > d)
    ret = ret - 1;
  return ret;
}

/* ---------------------------------------------------------------------------
 */
/*                       NORMAL DISTRIBUTION */
/* ---------------------------------------------------------------------------
 */

/* A simple normal distribution calculation accurate at 1.0e-07 */

double norm(double d) {
  double t;
  double z;
  double ans;
  double poly;

  z = fabs(SQRT_TWO_INV * d);
  t = 1.0 / (1.0 + 0.5 * z);
  poly = -1.26551223 +
         t * (1.00002368 +
              t * (0.37409196 +
                   t * (0.09678418 +
                        t * (-0.18628806 +
                             t * (0.27886807 +
                                  t * (-1.13520398 +
                                       t * (1.48851587 +
                                            t * (-0.82215223 +
                                                 t * 0.17087277))))))));
  ans = t * exp(-z * z + poly);

  if (d <= 0.0)
    ans = 0.5 * ans;
  else
    ans = 1.0 - 0.5 * ans;

  return ans;
}

/* A more complex normal distribution calculation accurate at 1.0e-12 */

/* ------------------------------------------------------------------------ */
double norm_accurate(double X) {

  double A[10];
  double B[10];
  double C[10];
  double D[10];
  double P[10];
  double Q[10];
  double d, ans, y, ysq, Xden, Xnum;
  double THRESH, DEL, SQRPI;
  int i;

  A[1] = 3.16112374387056560;
  A[2] = 1.13864154151050156 * 100;
  A[3] = 3.77485237685302021 * 100;
  A[4] = 3.20937758913846947 * 1000;
  A[5] = 1.85777706184603153 * 0.1;

  B[1] = 2.36012909523441209 * 10;
  B[2] = 2.44024637934444173 * 100;
  B[3] = 1.28261652607737228 * 1000;
  B[4] = 2.84423683343917062 * 1000;

  C[1] = 5.64188496988670089 * 0.1;
  C[2] = 8.88314979438837594;
  C[3] = 6.61191906371416295 * 10;
  C[4] = 2.98635138197400131 * 100;
  C[5] = 8.81952221241769090 * 100;
  C[6] = 1.71204761263407058 * 1000;
  C[7] = 2.05107837782607147 * 1000;
  C[8] = 1.23033935479799725 * 1000;
  C[9] = 2.15311535474403846e-8;

  D[1] = 1.57449261107098347 * 10;
  D[2] = 1.17693950891312499 * 100;
  D[3] = 5.37181101862009858 * 100;
  D[4] = 1.62138957456669019 * 1000;
  D[5] = 3.29079923573345963 * 1000;
  D[6] = 4.36261909014324716 * 1000;
  D[7] = 3.43936767414372164 * 1000;
  D[8] = 1.23033935480374942 * 1000;

  P[1] = 3.05326634961232344 * 0.1;
  P[2] = 3.60344899949804439 * 0.1;
  P[3] = 1.25781726111229246 * 0.1;
  P[4] = 1.60837851487422766 * 0.01;
  P[5] = 6.58749161529837803 * 0.0001;
  P[6] = 1.63153871373020978 * 0.01;

  Q[1] = 2.56852019228982242;
  Q[2] = 1.87295284992346047;
  Q[3] = 5.27905102951428412e-1;
  Q[4] = 6.05183413124413191e-2;
  Q[5] = 2.33520497626869185e-3;

  d = X * SQRT_TWO_INV;

  y = fabs(d);

  THRESH = 0.46875;

  SQRPI = 5.6418958354775628695 * 0.1;

  /*	------------------------------------------------------------------
                          Evaluate  erf  for  |X| <= 0.46875
          ------------------------------------------------------------------ */
  if (y <= THRESH) {
    ysq = 0;
    ysq = y * y;
    Xnum = A[5] * ysq;
    Xden = ysq;
    for (i = 1; i < 4; i++) {
      Xnum = (Xnum + A[i]) * ysq;
      Xden = (Xden + B[i]) * ysq;
    }

    ans = d * (Xnum + A[4]) / (Xden + B[4]);

  } /* END if (y<= THRESH) */

  /*  -------------------------------------------------------------------
          Evaluate  erfc  for 0.46875 <= |X| <= 4.0
      ------------------------------------------------------------------- */
  else if (y <= 4) {
    Xnum = C[9] * y;
    Xden = y;
    for (i = 1; i < 8; i++) {
      Xnum = (Xnum + C[i]) * y;
      Xden = (Xden + D[i]) * y;
    }
    ans = (Xnum + C[8]) / (Xden + D[8]);
    ysq = (Integer_Part(y * 16)) / 16;
    DEL = (y - ysq) * (y + ysq);

    ans = exp(-ysq * ysq) * exp(-DEL) * ans;
    /* ------------------------------------------------
       Fix up for negative argument        , erf        , etc.
       ------------------------------------------------ */
    ans = (0.5 - ans) + 0.5;
    if (d < 0)
      ans = -ans;

  } /* END if (y <=4) */
  /*  ------------------------------------------------------------------
          Evaluate  erfc  for |X| > 4.0
          ------------------------------------------------------------------ */
  else {
    ysq = 1 / (y * y);
    Xnum = P[6] * ysq;
    Xden = ysq;
    for (i = 1; i < 5; i++) {
      Xnum = (Xnum + P[i]) * ysq;
      Xden = (Xden + Q[i]) * ysq;
    }
    ans = ysq * (Xnum + P[5]) / (Xden + Q[5]);
    ans = (SQRPI - ans) / y;
    ysq = (Integer_Part(y * 16)) / 16;
    DEL = (y - ysq) * (y + ysq);

    ans = exp(-ysq * ysq) * exp(-DEL) * ans;

    /* -------------------------------------------
       Fix up for negative argument        , erf        , etc.
       ------------------------------------------- */
    ans = (0.5 - ans) + 0.5;
    if (d < 0)
      ans = -ans;

  } /* END else if (x>4.0) */

  /* ------------------------------------------------------------------
     The definitive answer
     ------------------------------------------------------------------ */
  ans = 0.5 * (1 + ans);

  return ans;
}

/* ----------------------------------------------------------------------------------
 */
/*                                 GAUSSIAN */
/* ----------------------------------------------------------------------------------
 */

double gauss(double x) {
  double result;

  result = INV_SQRT_TWO_PI * exp(-0.5 * x * x);

  return (result);
}

/* ========================================================================= */

/* -------------------------------------------------------------------------
    inverse to the cumulative normal        , N^(-1)(p)
    using newton-raphson method to solve the implicit equation
                        N(x)=target
    One can introduce a start number if one knows a number that satisfies the
    following equation:
                        fabs(norm(start))<fabs(target)
    in order to make the iteration faster
   ------------------------------------------------------------------------- */

double inv_cumnorm_newton(double start, double target) {
  double answer, diff;
  double precision = 1e-09;
  int sign;

  if (target < 0.5) {
    target = 1.0 - target;
    start = 1.0 - start;
    sign = -1;
  } else
    sign = 1;
  if (norm(start) <= target)
    answer = start;
  else /*IF start too big        , starts at zero*/
    answer = 0.0;

  diff = target - norm(answer);

  while (fabs(diff) > precision) {
    answer += diff / (exp(-(answer * answer) / 2) * INV_SQRT_TWO_PI);
    diff = target - norm(answer);
  }
  return (sign * answer);
}

/* ========================================================================= */

/* -------------------------------------------------------------------------
    inverse to the cumulative normal        , N^(-1)(p)
        This method uses a polynomial approximation which is more accurate than
    the NAG one (especially for the tails of the distribution)
    It is quicker than the inverse_cum_newton function
   ------------------------------------------------------------------------- */

double inv_cumnorm_fast(double u) {
  double a[4] = {2.50662823884, -18.61500062529, 41.39119773534,
                 -25.44106049637};

  double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826,
                 3.13082909833};

  double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
                 0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
                 0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

  double x, r;

  x = u - 0.5;

  if (fabs(x) < 0.42) {
    r = x * x;
    r = x * (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) /
        ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1);
    return (r);
  }

  r = u;

  if (x > 0.0)
    r = 1.0 - u;

  r = log(-log(r));
  r = c[0] +
      r * (c[1] +
           r * (c[2] +
                r * (c[3] +
                     r * (c[4] +
                          r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));

  if (x < 0.0)
    r = -r;

  return (r);
}

/* --------------------------------------------------------------------------
    inverse to the cumulative normal        , N^(-1)(p)
    this method has been proposed by Nag        , but unfortunately truncates
    too early (at +-2) the normal distribution tales
    Therefore        , the Newton-Raphson or the _fast method is more accurate
    for purposes like balance sampling
   ------------------------------------------------------------------------- */

double normal_dist_inv(double p) {
  double s, x, y;
  double yd, r2pi;
  double pi = 3.14;

  double p1 = 1.0;
  double p2 = -41.633627626163736;
  double p3 = 64.941289794046639;
  double p4 = -36.893553865736869;
  double p5 = 8.9843748879492917;
  double p6 = -0.84862050999166811;
  double p7 = 0.017992274523223912;
  double p8 = 0.00012434847784254829;
  double q1 = 1.0;
  double q2 = -41.800294292830407;
  double q3 = 71.849672176185038;
  double q4 = -46.455347140717776;
  double q5 = 13.576963141528349;
  double q6 = -1.7148871823011679;
  double q7 = 0.069291086876167264;
  double q8 = 1.5e-09;
  double t2t1 = 0.84162123357333769;
  double t2t2 = -6.0020630307113762;
  double t2t3 = -0.47226172781013548;
  double t2t4 = 68.25030402764662;
  double t2t5 = -96.239667846681925;
  double t2t6 = -10.47655255702662;
  double t2b1 = 1.0;
  double t2b2 = -11.37563668160269;
  double t2b3 = 41.338780429985753;
  double t2b4 = -43.59278882675467;
  double t2b5 = -21.255899751727728;
  double t2b6 = 25.0861060307781;
  double t3t1 = -3.0124726079130561;
  double t3t2 = 2.9824192543096469;
  double t3t3 = 22.70676420727861;
  double t3t4 = 7.6707637908225523;
  double t3t5 = 0.55191866290296669;
  double t3t6 = 0.0079854430460765382;
  double t3t7 = 4.1501841515746553e-06;
  double t3b1 = 1.0;
  double t3b2 = 1.2397581258179219;
  double t3b3 = -12.05350980555889;
  double t3b4 = -12.023590582199259;
  double t3b5 = -2.3730357962006429;
  double t3b6 = -0.11934402820315079;
  double t3b7 = -0.001216250189896074;

  if (p > 0.0 && p < 1.0) {
    x = p - 0.5;
    if (fabs(x) <= 0.3) {
      /*
       *    Break up range
       *    For 0.2 <= p <= 0.8 we use a Rational Tchebychev (p/q)
       */
      r2pi = sqrt(pi * 2.0);
      s = x * r2pi;
      x = s * s;
      y = p1 +
          x * (p2 +
               x * (p3 + x * (p4 + x * (p5 + x * (p6 + x * (p7 + x * p8))))));
      yd = q1 +
           x * (q2 +
                x * (q3 + x * (q4 + x * (q5 + x * (q6 + x * (q7 + x * q8))))));
      y /= yd;
      return s * y;
    }
    if (fabs(x) > 0.3) {
      /*
       *    For 0.08 <= p <= 0.2 || 0.8 <= p <= 0.92
       *    We use rational Tchebychev (t2t/t2b)
       */
      if (x > 0)
        s = 1;
      else
        s = -1;

      if (fabs(x) <= 0.42) {
        /*
         *        Break up range some more
         */
        x = fabs(x) - 0.3;
        y = t2t1 + x * (t2t2 + x * (t2t3 + x * (t2t4 + x * (t2t5 + x * t2t6))));
        yd =
            t2b1 + x * (t2b2 + x * (t2b3 + x * (t2b4 + x * (t2b5 + x * t2b6))));
        y /= yd;
        return s * y;
      }
      if (fabs(x) > 0.42) {
        /*
         *        The following case handles asymptotic behaviour.
         *        Unfortunately we must make a transformation.
         */
        if (p > 0.5)
          p = 1.0 - p;
        x = sqrt((double)(-2.0) * log((double)p));
        y = t3t1 +
            x * (t3t2 +
                 x * (t3t3 + x * (t3t4 + x * (t3t5 + x * (t3t6 + x * t3t7)))));
        yd = t3b1 +
             x * (t3b2 +
                  x * (t3b3 + x * (t3b4 + x * (t3b5 + x * (t3b6 + x * t3b7)))));
        y = y / yd + x;
        return s * y;
      }
    }

  } else {
    /*
     *    Error exits
     */
    if (p <= 0.0) {
      return 0;
    } else if (p >= 1.0) {
      return 0;
    }
  }
  return 0.0;
}

/*** bivariate normal cum distribution function        , with correlation rho
 * ***/

double bivar(double x, double y, double rho) {

  int i;

  double c1, c2, cos1, cos2, p, power1, power2, r, r1, root, rpower, rs2,
      rthalf, s1, s1c1, s2, s2c2, sine1, sine2, theta1, theta2, twopi;
  double tni[16];
  double result;

  double zero = 0.0e0;
  double quart = 0.25e0;
  double half = 0.5e0;
  double one = 1.0e0;
  double two = 2.0e0;

  static double b[16] = {0.0e0,
                         1.2533139981133488e+00,
                         -9.9999651036917565e-01,
                         6.2662474713886428e-01,
                         -3.3317729558070630e-01,
                         1.5620747416133978e-01,
                         -6.5782806339063672e-02,
                         2.4915538851209170e-02,
                         -8.3488914684103226e-03,
                         2.3979262116288155e-03,
                         -5.6652772989060614e-04,
                         1.0508139384820295e-04,
                         -1.4495937454221355e-05,
                         1.3832041191137142e-06,
                         -8.0839994686965162e-08,
                         2.1681406621770789e-09};

  /* g01hact.cxx
   *
   * Copyright 1990 Numerical Algorithms Group.
   *
   * NAG C Library
   *
   * Computes the bivariate normal probability;
   * prob(X.le.x        , Y.LE.y)
   * where X and Y are standardised with correlation  RHO.
   * The algorithm given by D. R. Divgi is used with a 15
   * term series. Ref : Divgi        , d.r.        , (1979)        , Calculation
   * of univariate and bivariate normal probability functions. The Annals of
   * Statistics        , 7        ,4        ,903-910
   *
   * Derived from: NAG Fortran Library G01HAF (Last Revision Mark 14).
   *
   * C version: July 1990        , Bode Meduoye        , NAG Ltd        , Oxford
   * , U.K.
   *
   * Mark 1        , 1990
   */

  result = zero;

  rthalf = sqrt(half);
  twopi = SRT_PI * two;

  if (rho == zero)
    result = norm(x) * norm(y);
  else if (x == zero && y == zero)
    result = quart + asin(rho) / twopi;
  else {
    root = sqrt(one - rho * rho);
    if (root != zero) {
      r1 = sqrt(x * x + y * y - two * rho * x * y);
      r = r1 / root;
      rs2 = -half * r * r;
    }
    if (root == zero || rs2 <= LOWFLOOR) /*if rs2 is too small*/
    {                                    /*it returns a result for rho=+-1*/
      if (rho < zero) {
        if (x >= -y)
          result = norm(x) + norm(y) - 1;
        if (x < -y)
          result = zero;
      } else {
        if (x > y) {
          result = norm(y);
        } /*result=norm(MIN(x        ,y))*/
        else {
          result = norm(x);
        }
      }
    } else {
      sine1 = fabs(x) / r;
      sine2 = fabs(y) / r;
      cos1 = (rho * x - y) / r1;
      cos2 = (rho * y - x) / r1;
      s1 = -SIGN(one, x);
      s2 = -SIGN(one, y);
      c1 = SIGN(one, cos1);
      c2 = SIGN(one, cos2);
      s1c1 = s1 * c1;
      s2c2 = s2 * c2;
      cos1 = fabs(cos1);
      cos2 = fabs(cos2);
      if (cos1 > one)
        cos1 = one;
      if (cos2 > one)
        cos2 = one;
      theta1 = acos(cos1);
      theta2 = acos(cos2);
      tni[0] = SIGN(theta1, s1c1) + SIGN(theta2, s2c2);
      tni[1] = SIGN(sine1, s1c1) + SIGN(sine2, s2c2);
      power1 = SIGN(sine1, s1c1);
      power2 = SIGN(sine2, s2c2);
      rpower = r;
      /*
       *            15-term series
       */
      p = tni[0] - b[1] * tni[1] * rpower;
      for (i = 2; i < 16; i++) {
        power1 = power1 * cos1;
        power2 = power2 * cos2;
        rpower = rpower * r;
        tni[i] = ((i - 1) * tni[i - 2] + power1 + power2) / (i);
        p = p - b[i] * tni[i] * rpower;
      }
      p = p * exp(rs2) / twopi;
      if (c1 < zero)
        p = p + s1 * norm(-fabs(x));
      if (c2 < zero)
        p = p + s2 * norm(-fabs(y));
      if (s1 < zero && s2 < zero)
        p = p + one;
      result = p;
    }
  }

  if (result < zero) {
    result = zero;
  }

  return (result);

} /*end program binormale  */

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optintnrm(14)
 *
 * PURPOSE      	: Integration of Normal Distribution Functions
 *
 * DESCRIPTION  	: Useful integration formulae as described in paper
 * 		  "integrating the Normal Distribution Function: Some
 *			Useful Formulae" by S. Rady.
 *
 * CALLS		: bivar()
 *
 *
 * PARAMETERS   	:      	-
 *
 * RETURNS      	: ??          	- ??
 *
 *******************************************************************************/

/*
   Added OVE 04-dec-1995
   If b1 is equal to 0        , we assume that the integral is made on a Dirac
   density function centered on y=a1
   Therefore        , the integral is replaced by the value of the function we
   integrate        , evaluated in y=a1
*/

double srt_f_intnrm_i1(double a0, double b0, double a1, double b1) {
  double b2;
  double x;
  double y;

  double rho;

  double i1;

  if (b1 != 0.00) {
    b2 = sqrt(b0 * b0 + b1 * b1);
    x = a1 / b1;
    y = (a0 + a1) / b2;
    rho = b1 / b2;

    i1 = bivar(x, y, rho);
  } else {
    if (a1 > 0)
      i1 = norm((a0 + a1) / b0);
    else
      i1 = 0; /* Dirac not in [0;+inf[*/
  }

  return i1;

} /* END srt_f_optintnrm_i1() */

/*
   Added OVE 04-dec-1995
   If b1 is equal to 0        , we assume that the integral is made on a Dirac
   density function centered on y=a1
   Therefore        , the integral is replaced by the value of the function we
   integrate        , evaluated in y=a1
*/

double srt_f_intnrm_i2(double a0, double b0, double a1, double b1) {
  double b2;
  double x;
  double y;

  double rho;

  double i2;

  if (b1 != 0.00) {
    b2 = sqrt(b0 * b0 + b1 * b1);
    x = a1 / b1;
    y = (-a0 + a1) / b2;
    rho = b1 / b2;
    i2 = norm(x) - bivar(x, y, rho);
  } else {
    if (a1 > 0)
      i2 = norm((a0 - a1) / b0);
    else
      i2 = 0; /* Dirac not in [0;+inf[*/
  }

  return i2;

} /* END srt_f_optintnrm_i2() */

/*
   Added OVE 04-dec-1995
   If b1 is equal to 0        , we assume that the integral is made on a Dirac
   density function centered on y=a1
   Therefore        , the integral is replaced by the value of the function we
   integrate        , evaluated in y=a1
*/

double srt_f_intnrm_j1(double alpha, double a0, double b0, double a1,
                       double b1) {
  double a3;
  double scalar;

  double i1;
  double j1;

  if (b1 != 0.00) {
    a3 = a1 + (alpha * b1 * b1);
    scalar = exp((alpha * a1) + (0.5 * alpha * alpha * b1 * b1));
    i1 = srt_f_intnrm_i1(a0, b0, a3, b1);

    j1 = scalar * i1;
  } else {
    if (a1 > 0)
      j1 = exp(alpha * a1) * norm((a0 + a1) / b0);
    else
      j1 = 0; /* Dirac not in [0;+inf[*/
  }

  return j1;

} /* END srt_f_optintnrm_j1() */

/*
   Added OVE 04-dec-1995
   If b1 is equal to 0        , we assume that the integral is made on a Dirac
   density function centered on y=a1
   Therefore        , the integral is replaced by the value of the function we
   integrate        , evaluated in y=a1
*/
double srt_f_intnrm_j2(double alpha, double a0, double b0, double a1,
                       double b1) {
  double a3;
  double scalar;

  double i2;
  double j2;

  if (b1 != 0.00) {
    a3 = a1 + (alpha * b1 * b1);
    scalar = exp((alpha * a1) + (0.5 * alpha * alpha * b1 * b1));
    i2 = srt_f_intnrm_i2(a0, b0, a3, b1);
    j2 = scalar * i2;
  } else {
    if (a1 > 0)
      j2 = exp(alpha * a1) * norm((a0 - a1) / b0);
    else
      j2 = 0; /* Dirac not in [0;+inf[*/
  }
  return j2;

} /* END srt_f_optintnrm_j2() */

double srt_f_intnrm_l(double a, double b, double ap, double bp)

{
  double bs;
  bs = sqrt(b * b + bp * bp);

  return (1.0 / bs) * gauss((a - ap) / bs) *
         norm((bp * a / b + b * ap / bp) / bs);
}

double srt_f_intnrm_m1(double a, double b, double ap, double bp)

{

  return bp * norm_accurate(a / b) * gauss(ap / bp) +
         (bp * bp) * srt_f_intnrm_l(-a, b, ap, bp) +
         ap * srt_f_intnrm_i1(a, b, ap, bp);
}

double srt_f_intnrm_m2(double a, double b, double ap, double bp)

{

  return bp * norm_accurate(a / b) * gauss(ap / bp) -
         (bp * bp) * srt_f_intnrm_l(a, b, ap, bp) +
         ap * srt_f_intnrm_i2(a, b, ap, bp);
}

/* ========================================================================
   FUNC: trivar_dens (static)
   DESC: Function under numerical integration for the
         three dimentional normal distr.
   MODIFIES:
   DECLARATION:
   ======================================================================== */

static double trivar_dens(double x, va_list argptr) {
  double a2, a3, r12, r23, r13, ro, b1, b2;
  static double result;

  a2 = va_arg(argptr, double);
  a3 = va_arg(argptr, double);
  r12 = va_arg(argptr, double);
  r23 = va_arg(argptr, double);
  r13 = va_arg(argptr, double);
  ro = va_arg(argptr, double);
  b1 = (a2 - r12 * x) / sqrt(1 - r12 * r12);
  b2 = (a3 - r13 * x) / sqrt(1 - r13 * r13);
  result = exp(-x * x / 2) * bivar(b1, b2, ro);
  return (result);
}

/* ========================================================================
   FUNC: trivar
   DESC: threedimentional normal distribution
   MODIFIES:
   DECLARATION:
   ======================================================================== */

double trivar(double a1, double a2, double a3, double r12, double r23,
              double r13) {
  double ro, result;

  ro = (r23 - r12 * r13) / sqrt((1 - r12 * r12) * (1 - r13 * r13));
  result = sm_qsimp_list(trivar_dens, -INFINITY, a1, a2, a3, r12, r23, r13, ro);
  result /= sqrt(SRT_PI * 2);
  return (result);
}

/* ========================================================================
   FUNC: EvalStudDis
   DESC: Computes the  cumulative Student distribution with n degrees of freedom
   at a given point

                n : n. of degrees of freedom (must be <= 10). 0 corresponds to
   the Gaussian Cumulative distribution x:  abscissa where we want to evaluate
   the Cumulative distribution

   MODIFIES:
   DECLARATION:
   ======================================================================== */

double Student_Dis(const int n, const double x

)

{

  double Value, alfa, beta;

  switch (n) {

  case 0: // this is the gaussian case (i.e. it corresponds effectively to
          // Student degree -> Infinity)
    Value = norm(x);
    break;

  case 1:
    Value = 0.5 + atan(x) / SRT_PI;
    break;

  case 2:
    Value = 0.5 + x / 2 / sqrt(2 + x * x);
    break;

  case 3:
    alfa = sqrt(3.) * x;
    beta = (3. + x * x);
    Value = 0.5 + (alfa + beta * atan(x / sqrt(n))) / SRT_PI / beta;
    break;

  case 4:
    Value = 0.5 + x * (6. + x * x) / 2 /
                      sqrt((4 + x * x) * (4 + x * x) * (4 + x * x));
    break;

  case 5:
    alfa = sqrt(5.) * x * (25. + 3. * x * x);
    beta = 3. * (5. + x * x) * (5. + x * x);
    Value = 0.5 + (alfa + beta * atan(x / sqrt(n))) / SRT_PI / beta;
    break;

  case 6:
    Value = 0.5 + x * (135. + 30. * x * x + 2. * x * x * x * x) / 4. /
                      sqrt((6. + x * x) * (6. + x * x) * (6. + x * x) *
                           (6. + x * x) * (6. + x * x));
    break;

  case 7:
    alfa = sqrt(7.) * x * (1617. + 280. * x * x + 15. * x * x * x * x);
    beta = 15. * (7. + x * x) * (7. + x * x) * (7. + x * x);
    Value = 0.5 + (alfa + beta * atan(x / sqrt(n))) / SRT_PI / beta;
    break;

  case 8:
    Value = 0.5 + x * sqrt(2 + x * x / 4.) *
                      (1120 + 280. * x * x + 28. * x * x * x * x +
                       x * x * x * x * x * x) /
                      (8. + x * x) / (8. + x * x) / (8. + x * x) / (8. + x * x);
    break;

  case 9:
    alfa = 3 * x *
           (67797. + 13797. * x * x + 1155 * x * x * x * x +
            35 * x * x * x * x * x * x);
    beta = 35 * (9. + x * x) * (9. + x * x) * (9. + x * x) * (9. + x * x);
    Value = 0.5 + (alfa + beta * atan(x / sqrt(n))) / SRT_PI / beta;
    break;

  case 10:
    Value = 0.5 + x *
                      (196875. + 52500. * x * x + 6300. * x * x * x * x +
                       360. * x * x * x * x * x * x +
                       8. * x * x * x * x * x * x * x * x) /
                      16. /
                      sqrt((10. + x * x) * (10. + x * x) * (10. + x * x) *
                           (10. + x * x) * (10. + x * x) * (10. + x * x) *
                           (10. + x * x) * (10. + x * x) * (10. + x * x));
    break;

  default:
    break;
  }

  return Value;
}

/* ========================================================================
   FUNC: EvalStudDis
   DESC: Computes the  cumulative Student distribution with n degrees of freedom
   at a given point

                n : n. of degrees of freedom (must be <= 10). 0 corresponds to
   the Gaussian Cumulative distribution x:  abscissa where we want to evaluate
   the Cumulative distribution

   MODIFIES:
   DECLARATION:
   ======================================================================== */
double Student_Integrand(const double x, va_list argptr) {
  double m = va_arg(argptr, double);
  return exp(-0.5 * (m + 1.0) * log(1 + x * x / m));
}

double Student_Distribution(const double x, const double m,
                            const double sigma) {
  double sgn = x > 0.0 ? 1.0 : -1.0;
  double A =
      sgn * exp(gammln(0.5 * (m + 1.0)) + gammln(0.5 * m)) / sqrt(SRT_PI * m);
  return 0.5 + A * sm_qsimp_list(Student_Integrand, 0.0, sgn * x / sigma, m);
}

double Std_Student_Distribution(const double x, const double m) {
  return Student_Distribution(x, m, 1.0);
}

Err Inverse_Student_Diff(double x, double *f, double *df, double *data) {
  *f = Std_Student_Distribution(x, data[0]) - data[1];
  *df = Std_Student_Density(x, data[0]);
  return NULL;
}

// t should be between 0 and 1
Err Inverse_Std_Student_Distribution2(double *x, const double t, const double m,
                                      const double acc) {
  // Local Variables
  double HighGuess, LowGuess, *data;
  unsigned long i = 0;
  Err err = NULL;

  // Return if we are giving a known value
  if (t == 0.0) {
    *x = -INFINITY;
    return err;
  }
  if (t == 1.0) {
    *x = INFINITY;
    return err;
  }
  if (t == 0.5) {
    *x = 0.0;
    return err;
  }

  // Find boundaries on the root.
  if (t < 0.5) {
    HighGuess = 0.0;
    LowGuess = -1.0;
    while (Std_Student_Distribution(LowGuess, m) > t && i < 100) {
      HighGuess = LowGuess;
      LowGuess *= 2.0;
      i++;
    }
  } else {
    LowGuess = 0.0;
    HighGuess = 1.0;
    while (Std_Student_Distribution(HighGuess, m) < t && i < 100) {
      LowGuess = HighGuess;
      HighGuess *= 2.0;
      i++;
    }
  }
  if (i == 100)
    return serror("rtsafe hit max iterations");

  // Put the data for the function in an array
  data = dvector(0, 1);
  data[0] = m;
  data[1] = t;

  err = rtsafe_with_par(&Inverse_Student_Diff, LowGuess, HighGuess, acc, 100, x,
                        data);

  // Release the memory and head home
  free_dvector(data, 0, 1);
  return err;
}

/**************************************************************************************************************************************/
//
// Inverse_Student_Distribution
//
// A wrapper for the functions below
//
// Inputs:
//		t				value of the Student cumulative
// distribution 		m				degree of the
// student distribution 		acc				accuracy of the
// method
//
// Return value:
//		the inverse of the Student cumulative distribution function
//
/**************************************************************************************************************************************/
double Inverse_Student_Distribution(const double t, const unsigned int m,
                                    const double acc) {
  return Student_Dis_Inv(acc, t, m);
}

/* ========================================================================
   FUNC: Student_Dis_Inv
   DESC: Computes the inverse of the cumulative Student distribution
         with a number of degrees of freedom

                 degree: N. of degrees of freedom of the Student distribution
   (must be <=10 ) 0 correspond to the gaussian Cumulative Inverse distribution
                 xacc : accuracy of the numerical root finding algorithm
                 tVal: value of the Student cumulative distribution

   MODIFIES:
   DECLARATION:
   ======================================================================== */
#define MAXIT 100
#define THETA(a, b) ((a) > (b) ? 1 : 0)
double Student_Dis_Inv(const double xacc, const double tVal,
                       const unsigned int degree)

{
  void nrerror(char error_text[]);
  int j;
  double x1, x2;
  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;

  switch (degree) {

  case 0:
    return inv_cumnorm_fast(tVal);
    break;

  case 1:
    return tan(SRT_PI * (tVal - 0.5));
    break;
  case 2:
    return THETA(tVal, 0.5) * sqrt((4 * tVal - 4 * tVal * tVal - 1) /
                                   (2 * tVal * tVal - 2 * tVal)) -
           THETA(0.5, tVal) * sqrt((4 * tVal - 4 * tVal * tVal - 1) /
                                   (2 * tVal * tVal - 2 * tVal));
    break;

  default:

    // determines the boundaries where to find a solution

    if (tVal < 0.5 && tVal > 0.0) {
      x1 = -30;
      x2 = 0.0;
    } else if (tVal >= 0.5 && tVal < 1.0) {
      x1 = 0.0;
      x2 = 30.0;
    } else {
      goto end;
    }

    TryStudent_Dis_Inv(x1, tVal, degree, &fl, &df);
    TryStudent_Dis_Inv(x2, tVal, degree, &fh, &df);
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
      goto end;
    //		nrerror("Root must be bracketed in rtsafe");

    if (fl == 0.0)
      return x1;
    if (fh == 0.0)
      return x2;
    if (fl < 0.0) {
      xl = x1;
      xh = x2;
    } else {
      xh = x1;
      xl = x2;
    }
    rts = 0.5 * (x1 + x2);
    dxold = fabs(x2 - x1);
    dx = dxold;
    TryStudent_Dis_Inv(rts, tVal, degree, &f, &df);

    for (j = 1; j <= MAXIT; j++) {
      if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) ||
          (fabs(2.0 * f) > fabs(dxold * df))) {
        dxold = dx;
        dx = 0.5 * (xh - xl);
        rts = xl + dx;
        if (xl == rts)
          return rts;
      } else {
        dxold = dx;
        dx = f / df;
        temp = rts;
        rts -= dx;
        if (temp == rts)
          return rts;
      }

      if (fabs(dx) < xacc)
        return rts;

      TryStudent_Dis_Inv(rts, tVal, degree, &f, &df);
      if (f < 0.0)
        xl = rts;
      else
        xh = rts;
    }
    //	nrerror("Maximum number of iterations exceeded in rtsafe");
  end:
    return 0.0;

    break;
  }
}
#undef MAXIT
#undef THETA

void TryStudent_Dis_Inv(const double x, const double tVal,
                        const unsigned int degree, double *fVal, double *f1Val)

{

  double alfa, beta, gamma;

  switch (degree) {

  case 3:
    alfa = sqrt(3.) * x;
    beta = (3 + x * x);
    *fVal = x / sqrt(1. * degree) -
            tan((SRT_PI * beta * (tVal - 0.5) - alfa) / beta);
    *f1Val = 6. * sqrt(3.) / SRT_PI / (3. + x * x) / (3. + x * x);
    break;

  case 4:
    alfa = (4. + x * x);
    *fVal = (tVal - 0.5) * (tVal - 0.5) -
            x * x * (6. + x * x) * (6. + x * x) / 4. / alfa / alfa / alfa;
    *f1Val = 12. / sqrt(alfa * alfa * alfa * alfa * alfa);
    break;

  case 5:
    alfa = sqrt(5.) * x * (25. + 3. * x * x);
    beta = 3. * (5. + x * x) * (5. + x * x);
    *fVal = x / sqrt(1. * degree) -
            tan((SRT_PI * beta * (tVal - 0.5) - alfa) / beta);
    *f1Val =
        200. * sqrt(5.) / SRT_PI / (5. + x * x) / (5. + x * x) / (5. + x * x);
    break;

  case 6:
    alfa = (6. + x * x);
    *fVal = (tVal - 0.5) * (tVal - 0.5) -
            x * x * (135. + 30 * x * x + 2. * x * x * x * x) *
                (135. + 30 * x * x + 2. * x * x * x * x) / 16. / alfa / alfa /
                alfa / alfa / alfa;
    *f1Val = 405. / 2. / sqrt(alfa * alfa * alfa * alfa * alfa * alfa * alfa);
    break;

  case 7:
    alfa = sqrt(7.) * x * (1617. + 280 * x * x + 15. * x * x * x * x);
    gamma = (7. + x * x);
    beta = 15. * gamma * gamma * gamma;
    *fVal = x / sqrt(1. * degree) -
            tan((SRT_PI * beta * (tVal - 0.5) - alfa) / beta);
    *f1Val = 5488. * sqrt(7.) / 5. / SRT_PI / gamma / gamma / gamma / gamma;
    break;

  case 8:
    alfa = (8. + x * x);
    *fVal =
        (tVal - 0.5) * alfa * alfa * alfa * alfa -
        x * sqrt(2. + x * x / 4.) *
            (1120 + 280. * x * x + 28. * x * x * x * x + x * x * x * x * x * x);
    *f1Val = 4480. /
             sqrt(alfa * alfa * alfa * alfa * alfa * alfa * alfa * alfa * alfa);
    break;

  case 9:
    alfa = 3 * x *
           (67797. + 13797. * x * x + 1155 * x * x * x * x +
            35. * x * x * x * x * x * x);
    gamma = 9 + x * x;
    beta = 35. * gamma * gamma * gamma * gamma;
    *fVal = x / sqrt(1. * degree) -
            tan((SRT_PI * beta * (tVal - 0.5) - alfa) / beta);
    *f1Val = 2519424 / 35. / SRT_PI / gamma / gamma / gamma / gamma / gamma;
    break;

  case 10:
    alfa = 10. + x * x;
    *fVal =
        (tVal - 0.5) * 16. *
            sqrt(alfa * alfa * alfa * alfa * alfa * alfa * alfa * alfa * alfa) -
        x * (196875. + 52500. * x * x + 6300. * x * x * x * x +
             360. * x * x * x * x * x * x + 8. * x * x * x * x * x * x * x * x);
    *f1Val = 984375. / 8 /
             sqrt(alfa * alfa * alfa * alfa * alfa * alfa * alfa * alfa * alfa *
                  alfa * alfa);
    break;

  default:
    *fVal = 1.0;
    *f1Val = 1.0;
    break;
  }
}

#undef INFINITY
#undef LOWFLOOR
