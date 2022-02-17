/* ==========================================================================
   FILE_NAME:	num_f_gamma.c

   PURPOSE:     computes complete and incomplete gamma functions
                (cf NRC pages 216ff)
   ========================================================================== */

#include "utallhdr.h"
#include <math.h"

#define ITMAX 100
#define GAMMA_EPS 3.0e-7
#define FPMIN 1.0e-30

double gammln(double xx) {
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146,     -86.50532032941677,
                          24.01409824083091,     -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++)
    ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

double betacmp(double alpha, double beta)
/*Returns the value of the Beta function with parameters alpha and beta*/
{
  double gammln(double x);

  return (exp(gammln(alpha) + gammln(beta) - gammln(alpha + beta)));
}

double gammp(double a, double x) {
  void gcf(double *gammcf, double a, double x, double *gln);
  void gser(double *gamser, double a, double x, double *gln);
  double gamser, gammcf, gln;

  if (x < 0.0 || a <= 0.0)
    serror("Invalid arguments in routine gammp");
  if (x < (a + 1.0)) {
    gser(&gamser, a, x, &gln);
    return gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return 1.0 - gammcf;
  }
}

void gcf(double *gammcf, double a, double x, double *gln) {
  double gammln(double xx);
  int i;
  double an, b, c, d, del, h;

  *gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++) {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < FPMIN)
      d = FPMIN;
    c = b + an / c;
    if (fabs(c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < GAMMA_EPS)
      break;
  }
  if (i > ITMAX)
    serror("a too large  , ITMAX too small in gcf");
  *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}

void gser(double *gamser, double a, double x, double *gln) {
  double gammln(double xx);
  int n;
  double sum, del, ap;

  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0)
      serror("x less than 0 in routine gser");
    *gamser = 0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++) {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs(del) < fabs(sum) * GAMMA_EPS) {
        *gamser = sum * exp(-x + a * log(x) - (*gln));
        return;
      }
    }
    serror("a too large  , ITMAX too small in routine gser");
    return;
  }
}

double gammq(double a, double x) {
  void gcf(double *gammcf, double a, double x, double *gln);
  void gser(double *gamser, double a, double x, double *gln);
  double gamser, gammcf, gln;

  if (x < 0.0 || a <= 0.0)
    serror("Invalid arguments in routine gammq");
  if (x < (a + 1.0)) {
    gser(&gamser, a, x, &gln);
    return 1.0 - gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return gammcf;
  }
}

double betacf(double a, double b, double x)
/*Used by betaincmp: Evaluates continued fraction for incomplete beta
function by modified Lentz's method*/
{
  int m, m2;
  double aa, c, d, de1, h, qab, qam, qap;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;
  d = 1.0 - qab * x / qap;
  if (fabs(d) < FPMIN)
    d = FPMIN;
  d = 1.0 / d;
  h = d;
  for (m = 1; m <= ITMAX; m++) {
    m2 = 2 * m;
    aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN)
      d = FPMIN;
    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qam + m2));
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN)
      d = FPMIN;
    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    de1 = d * c;
    h *= de1;
    if (fabs(de1 - 1.0) < GAMMA_EPS)
      break;
  }
  if (m > ITMAX)
    smessage("a or b too big  , or ITMAX too small in betacf");
  return h;
}

double betaincmp(double a, double b, double x)
/*Returns the incomplete Beta function with parameters a and b when
it is computed up to x*/
{
  double betacf(double a, double b, double x);
  double gammln(double xx);
  double bt;

  if (x < 0.0 || x > 1.0)
    smessage("Bad x in routine betaincmp");
  if (x == 0.0 || x == 1.0)
    bt = 0.0;
  else
    bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) +
             b * log(1.0 - x));
  if (x < (a + 1.0) / (a + b + 2.0))
    return bt * betacf(a, b, x) / a;
  else
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

/* Calculates Factorial(n) */
double fact(long n) {
  if (n == 0)
    return (1);
  else
    return (n * fact(n - 1));
}

#undef GAMMA_EPA