/* randist/levy.c
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
#include "alpha_stable_random.h"



/* The stable Levy probability distributions have the form

   p(x) dx = (1/(2 pi)) \int dt exp(- it x - |c t|^alpha)

   with 0 < alpha <= 2. 

   For alpha = 1, we get the Cauchy distribution
   For alpha = 2, we get the Gaussian distribution with sigma = sqrt(2) c.

   Fromn Chapter 5 of Bratley, Fox and Schrage "A Guide to
   Simulation". The original reference given there is,

   J.M. Chambers, C.L. Mallows and B. W. Stuck. "A method for
   simulating stable random variates". Journal of the American
   Statistical Association, JASA 71 340-344 (1976).

   */

int RandomAlphaStable (     void *random,
                            double *alphaSequence,
                            long nbPaths,
                            const double c, 
                            const double alpha)
{
    static char routine[] = "RandomAlphaStable";
    int status = FAILURE;
    double u, v, t, s;
    long j = 0;
    
    for(j=0;j<nbPaths;j++)
    {
        u = M_PI * (RandomGeneratorGet(random) - 0.5);

        if (fabs(alpha - 1.)<3e-16)		/* cauchy case */
        {
            t = tan (u);
            alphaSequence[j] = c * t;
            continue;
        }

        v = -log(RandomGeneratorGet(random));

        if (fabs(alpha - 2.)<3e-16)             /* gaussian case */
        {
            t = 2 * sin (u) * sqrt(v);
            alphaSequence[j] = c * t;
            continue;
        }

        /* general case */
        t = sin (alpha * u) / pow (cos (u), 1 / alpha);
        s = pow (cos ((1 - alpha) * u) / v, (1 - alpha) / alpha);

        alphaSequence[j] = c * t * s;
        continue;
    }
    
    status = SUCCESS;
    return status;
}

/* The following routine for the skew-symmetric case was provided by
   Keith Briggs.

   The stable Levy probability distributions have the form

   2*pi* p(x) dx

     = \int dt exp(mu*i*t-|sigma*t|^alpha*(1-i*beta*sign(t)*tan(pi*alpha/2))) for alpha!=1
     = \int dt exp(mu*i*t-|sigma*t|^alpha*(1+i*beta*sign(t)*2/pi*log(|t|)))   for alpha==1

   with 0<alpha<=2, -1<=beta<=1, sigma>0.

   For beta=0, sigma=c, mu=0, we get gsl_ran_levy above.

   For alpha = 1, beta=0, we get the Lorentz distribution
   For alpha = 2, beta=0, we get the Gaussian distribution

   See A. Weron and R. Weron: Computer simulation of Lévy alpha-stable 
   variables and processes, preprint Technical University of Wroclaw.
   http://www.im.pwr.wroc.pl/~hugo/Publications.html

*/

int    RandomAlphaStableSkew (  void *random,
                                double *alphaSequence,
                                long nbPaths,
                                const double c, 
                                const double alpha,
                                const double beta)
{
    static char routine[] = "RandomAlphaStableSkew";
    int status = FAILURE;
    double V, W, X;
    long j = 0;
    
    for(j=0;j<nbPaths;j++)
    {
        if (fabs(beta)<3e-16)  /* symmetric case */
        {
            status = RandomAlphaStable ( random,
                                alphaSequence,
                                nbPaths,
                                c, 
                                alpha);
            if(status == FAILURE) goto RETURN;
            break;
        }

        V = M_PI * (RandomGeneratorGet(random) - 0.5);

        do
        {
            W = -log(RandomGeneratorGet(random));
        }
        while (fabs(W)<3e-16);

        if (fabs(alpha-1)<3e-16)
        {
            X = ((M_PI_2 + beta * V) * tan (V) -
	        beta * log (M_PI_2 * W * cos (V) / (M_PI_2 + beta * V))) / M_PI_2;
            alphaSequence[j] = c * (X + beta * log (c) / M_PI_2);
            continue;
        }
        else
        {
            double t = beta * tan (M_PI_2 * alpha);
            double B = atan (t) / alpha;
            double S = pow (1 + t * t, 1/(2 * alpha));

            X = S * sin (alpha * (V + B)) / pow (cos (V), 1 / alpha)
	    *   pow (cos (V - alpha * (V + B)) / W, (1 - alpha) / alpha);
            alphaSequence[j] = c * X;
            continue;
        }
    }
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);        
    }
    return status;
}



/* -------------------------------------------------------------------------
// AlphaStableDeviates
// This function generates itself the uniform random
// they are not parameters of the function.
*/
int AlphaStableDeviates (   double *alphaSequence,
                            long nbPaths,
                            const double c, 
                            const double alpha)
{
    static char routine[] = "AlphaStableDeviates";
    int status = FAILURE;
    void *random = CreateRandomGenerator(-12);

    status = RandomAlphaStable(random, alphaSequence,nbPaths,c,alpha);
    if(status == FAILURE) goto RETURN;

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

/* -------------------------------------------------------------------------
// AlphaStableSkewDeviates
// This function generates itself the uniform random
// they are not parameters of the function.
*/
int AlphaStableSkewDeviates(double *alphaSequence,
                            long nbPaths,
                            const double c, 
                            const double alpha,
                            const double beta)
{
    static char routine[] = "AlphaStableSkewDeviates";
    int status = FAILURE;
    void *random = CreateRandomGenerator(-12);

    status = RandomAlphaStableSkew(random, alphaSequence,nbPaths,c,alpha,beta);
    if(status == FAILURE) goto RETURN;

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}
