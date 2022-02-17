/* randist/gamma.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include "error2.h"
#include "random_utils.h"
#include "proba_utils.h"
#include "gamma_random.h"

#ifndef M_E
#define M_E        2.71828182845904523536028747135	/* e */
#endif

#ifndef M_PI
#define M_PI	   3.14159265358979323846264338328      /* pi */
#endif

/* The Gamma distribution of order a>0 is defined by:

   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a1 and a2 with the same scale parameter b, then
   X+Y has gamma distribution of order a1+a2.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

int RandomGamma (void *random, double *gammaSequence, long nbPaths, const double a, const double b)
{
    long j = 0;
    for(j=0;j<nbPaths;j++)
    {
        gammaSequence[j] = gsl_ran_gamma(random, a,b);
    }
    return SUCCESS;
}


double gsl_ran_gamma (void *random, const double a, const double b)
{
  /* assume a > 0 */
  unsigned int na = (int) floor (a);

  if (a == na)
    {
      return b * gsl_ran_gamma_int (random, na);
    }
  else if (na == 0)
    {
      return b * gamma_frac (random, a);
    }
  else
    {
      return b * (gsl_ran_gamma_int (random, na) + gamma_frac (random, a - na)) ;
    }
}

double gsl_ran_gamma_int (void *random, const unsigned int a)
{
  if (a < 12)
    {
      unsigned int i;
      double prod = 1;

      for (i = 0; i < a; i++)
	{
	  prod *= RandomGeneratorGet(random);
	}

      /* Note: for 12 iterations we are safe against underflow, since
	 the smallest positive random number is O(2^-32). This means
	 the smallest possible product is 2^(-12*32) = 10^-116 which
	 is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gamma_large (random, (double) a);
    }
}

static double gamma_large (void *random, const double a)
{
  /* Works only if a > 1, and is most efficient if a is large

     This algorithm, reported in Knuth, is attributed to Ahrens.  A
     faster one, we are told, can be found in: J. H. Ahrens and
     U. Dieter, Computing 12 (1974) 223-246.  */

  double sqa, x, y, v;
  sqa = sqrt (2 * a - 1);
  do
    {
      do
	{
	  y = tan (M_PI * RandomGeneratorGet(random));
	  x = sqa * y + a - 1;
	}
      while (x <= 0);
      v = RandomGeneratorGet(random);
    }
  while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

  return x;
}

static double gamma_frac (void *random, const double a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;
  p = M_E / (a + M_E);
  do
    {
      u = RandomGeneratorGet(random);
      v = RandomGeneratorGet(random);

      if (u < p)
	{
	  x = exp ((1 / a) * log (v));
	  q = exp (-x);
	}
      else
	{
	  x = 1 - log (v);
	  q = exp ((a - 1) * log (x));
	}
    }
  while (RandomGeneratorGet(random) >= q);

  return x;
}


int GammaDeviates (   double *alphaSequence,
                            long nbPaths,
                            const double a, 
                            const double b)
{
    static char routine[] = "GammaDeviates";
    int status = FAILURE;
    void *random = CreateRandomGenerator(-7);

    status = RandomGamma(random, alphaSequence,nbPaths,a,b);
    if(status == FAILURE) goto RETURN;

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}
