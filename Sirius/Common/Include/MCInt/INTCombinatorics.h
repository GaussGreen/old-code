
#ifndef __INTCOMBINATORICS_H__
#define __INTCOMBINATORICS_H__


#include "INTSTL.h"


//----------------------------------------------------------------------------------------------

// ln(\Gamma(x)), x in |R^+
 double gammln(double x);


// ln(n!), n in |N
template<class T> double factln(T n)
{
  if (n < 0) ERROR("In factln: n must be non-negative. ")

  static double a[101];

  if (n <= 1) return 0.0;
  if (n <= 100) 
      return a[n] ? a[n] : (a[n] = gammln(n + 1.0)); 
  else 
      return gammln(n + 1.0);
}

// Factorial n! = n*(n-1)*...*1, n in |N
 double factr(int n);



// Binomial coefficient (n choose k) = n! / (k! * (n - k)!), n >= k in |N
template<class T> double bico(T n, T k)
{
  if (n - k < 0) ERROR("In bico: n - k must be non-negative. ")

  return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
}

//----------------------------------------------------------------------------------------------


#endif
