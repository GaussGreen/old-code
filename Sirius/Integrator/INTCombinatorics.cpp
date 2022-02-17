#include	"stdafx.h"


#include "MCInt\INTCombinatorics.h"



//----------------------------------------------------------------------------------------------

// ln(\Gamma(x)), x in |R^+
// From NR in C, p. 214
double gammln(double xx)
{
  static double cof[6] = {	 76.18009172947146, -86.50532032941677, 24.01409824083091,  
							-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
  double x, y, tmp, ser;
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5)*log(tmp);
  ser = 1.000000000190015;

  for (j = 0;j <= 5;j++)  
	  ser += cof[j] / ++y;

  return  - tmp + log(2.5066282746310005 * ser / x);
}

//----------------------------------------------------------------------------------------------

// Factorial n! = n*(n-1)*...*1, n in |N
// From NR in C, p. 214
double factr(int n)
{
  if (n < 0) ERROR("In factr: n must be non-negative. ")

  static int ntop = 4;
  static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
  int j;

  if (n > 32) return exp(gammln(n + 1.0));
  while (ntop < n) 
    {
      j = ntop++;
      a[ntop] = a[j] * ntop;
    }
  return a[n];
}

//----------------------------------------------------------------------------------------------

