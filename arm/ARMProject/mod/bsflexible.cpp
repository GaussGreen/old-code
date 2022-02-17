/*
 * $Log: bsflexible.cpp,v $
 * Revision 1.2  2003/09/09 09:07:04  mab
 * RCS comments
 *
 */

/*----------------------------------------------------------------------------*
    bsflexible.cpp
 
    Black and Scholes Analytics for Vanilla and Exotic Options
 
    Copyright (c) 2002
*----------------------------------------------------------------------------*/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "armglob.h"
#include "linalg.h"
#include "gaussian.h"
#include "newton.h"
#include "bsflexible.h"

#define kMinSpot 1.0e-4
#define kMaxSpot 1.0e10    // limits for bracketing spot root
#define kSpotTol 1.0e-4
#define kIterMax 20

#define K_VAL_ERROR     -100




/*----------------------------------------------------------------------------*
   NAME    bsOption2

   Compute value of european Vanilla option using Black & Scholes formula

   Output:
   the function returns the Black Scholes price for an equity type option,
   it uses the underlying more general pricing function bsflexible .

*----------------------------------------------------------------------------*/

double bsOption2(double spot,
                double strike,
                double volatility,
                double dividend,
                double discountRate,
                double maturity,
                double CallPut)
{
    double forward,totalvolatility,bondprice,value;

	forward=spot*exp((discountRate-dividend)*maturity);

	totalvolatility=volatility * sqrt(maturity);

	bondprice=exp(-discountRate*maturity);

    value = bsflexible(forward,
                       totalvolatility,
                       bondprice,
                       strike,
                       CallPut);

    return(value);
};

/*----------------------------------------------------------------------------*
   NAME    bsflexible

   Compute value of european Vanilla option using Black & Scholes formula
   in its more general formulation

   Input:

   forward : the expected value of the underlying under the risk neutral 
     measure, assuming that the underlying has a lognormal distribution of
	 at reset time under the given probability..
	 
   totalvolatility: the volatility of underlying at reset time ..
     In cas of a process with a deterministic volatility, this is 
     the square root of the integral of the square of the 
     deterministic volatility of the process  between the start of 
	 the option an the reset time of the option.

   bondprice : the price of a bond with notionnal=1, maturiring at the payment
     time of the option.(payment time >= reset time)

  


   Output:

   the Black Scholes price for an equity type option,
   

*----------------------------------------------------------------------------*/

double bsflexible(double forward,
                       double totalvolatility,
                       double bondprice,
                       double strike,
                       double CallPut)
{
	double  value, d1;

    d1 = (log(forward/strike))/totalvolatility+0.5*totalvolatility ;

    value = CallPut*bondprice*(forward*cdfNormal(CallPut*d1)-
				strike*cdfNormal(CallPut*(d1-totalvolatility)));

    return(value);
};


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
