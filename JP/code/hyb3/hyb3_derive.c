/****************************************************************************/
/*      Calculate delta and gamma using 4th interpolation.                  */
/****************************************************************************/
/*      DERIVE.c                                                            */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"



/*****  Delta  **************************************************************/
/*
*       Delta of an option. We have 5 prices and spot level (0, 1, 2, 3, 4)
*  	The central one, indexed 2 corresponds to today's level. We use 
*	them to find the price of the option for an small up and down moves
*	with a 4th degree interpolation. Then we calculate the delta.
*/
double 	Delta (	double 	*Under,                                                 /* Underlying prices */
                double 	*Option)                                                /* Option prices */
{
        double 	
                U,	                                                        /* Central value for the underlying */
                e,	                                                        /* Perturbation of underlying */
                P1,	                                                        /* Value of the option for an up move */
                P2,                                                             /* Value of the option for a down move */
                delta;

        U = Under[2]; e = U * 0.001;

        d4interp (U-e, &P1, Under, Option);
        d4interp (U+e, &P2, Under, Option);
        
        delta = (P2 - P1) / (2. * e);	                                        /*  Approximates df(x)/dx by (f(x+e)-f(x-e))/2e.  */

        return (delta);

}  /* Delta */



/*****  Gamma  **************************************************************/
/*
*       Gamma of an option. 
*/
double 	Gamma (	double 	*Under,                                                 /* Underlying prices */
                double 	*Option)                                                /* Option prices */
{
        double 	
                U,	                                                        /* Central value for the underlying */
                e,	                                                        /* Perturbation of underlying */
                P,	                                                        /* Current value of the option */
                P1,	                                                        /* Value of the option for an up move */
                P2,                                                             /* Value of the option for a down move */
                gamma;

        U = Under[2]; e = U * 0.001; P = Option[2];

        d4interp (U-e, &P1, Under, Option);
        d4interp (U+e, &P2, Under, Option);
        
        gamma = (P1 + P2 - 2. * P) / (e * e);

        return (gamma);

}  /* Gamma */



/*****  Delta1  *************************************************************/
/*
*       Left or right delta of an option. We have 3 prices and spot levels. 
*  	The extreme one, indexed 0 corresponds to today's level. We use 
*	them to find the price of the option for an small move with a 
*	quadratic interpolation. Then we calculate the delta.
*/
double 	Delta1 (double 	*Under,                                                 /* Underlying prices */
                double 	*Option,                                                /* Option prices */
                char	LoR)                                                    /* Left or right derivative */
{
        double 	
                U,	                                                        /* Current value for the underlying */
                e,	                                                        /* Perturbation of underlying */
                P,	                                                        /* Current value of the option */
                P1,	                                                        /* Value of the tweaked option */
                delta;

                
        U = Under[0];
        P = Option[0];
        
        if (LoR == 'L')
        {        
                e = U * 0.001 * SIGN (Under[-1] - Under[0]);	                /* Does not assume underlying is increasing with the index */

                qinterp (       Under - 3,	                                /* qinterp uses index from 1 to 3. Hence the offset */
                                Option - 3,
                                U + e,	                                        
                                &P1,
                                0.);
                
                delta = (P1 - P) / e;	                                        /*  Approximates df(x)/dx by (f(x)-f(x-e))/e */

                return (delta);
        }
        else 
        {        
                e = U * 0.001 * SIGN (Under[1] - Under[0]);

                qinterp (       Under - 1,
                                Option - 1,
                                U + e,
                                &P1,
                                0.);
                
                delta = (P1 - P) / e;	                                        /*  Approximates df(x)/dx by (f(x)-f(x-e))/e */

                return (delta);
                
        }  /* if then else */

}  /* Delta1 */



/*****  Gamma1  *************************************************************/
/*
*       Gamma of an option. 
*/
double 	Gamma1 (double 	*Under,                                                 /* Underlying prices */
                double 	*Option,                                                /* Option prices */
                char	LoR)                                                    /* Left or right derivative */
{
        double 	
                U,	                                                        /* Current value for the underlying */
                e,	                                                        /* Perturbation of underlying */
                P,	                                                        /* Current value of the option */
                P1,	                                                        /* Value of the option for an up move */
                P2,                                                             /* Value of the option for an down move */
                gamma;

        
        U = Under[0];
        P = Option[0];
        
        if (LoR == 'L')
        {        
                e = U * 0.001 * SIGN (Under[-1] - Under[0]);

                qinterp (       Under - 3,	                                /* qinterp uses index from 1 to 3. Hence the offset */
                                Option - 3,
                                U + e,	                                        
                                &P1,
                                0.);
                
                qinterp (       Under - 3,	                                
                                Option - 3,
                                U + 2. * e,	                                        
                                &P2,
                                0.);
                
                gamma = (P2 + P - 2. * P1) / (e * e);

                return (gamma);
        }
        else 
        {        
                e = U * 0.001 * SIGN (Under[1] - Under[0]);

                qinterp (       Under - 1,
                                Option - 1,
                                U + e,
                                &P1,
                                0.);
                
                qinterp (       Under - 1,
                                Option - 1,
                                U + 2. * e,
                                &P2,
                                0.);
                
                gamma = (P2 + P - 2. * P1) / (e * e);

                return (gamma);
                
        }  /* if then else */

}  /* Gamma1 */
