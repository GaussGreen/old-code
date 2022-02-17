/*
 * $Log: bsxtic.cpp,v $
 * Revision 1.8  2004/03/24 13:43:59  rguillemot
 * SUMMIT Payment Lag
 *
 * Revision 1.7  2004/01/13 09:15:26  jpriaudel
 * added : CallForwardValue and PutForwardValue
 *
 * Revision 1.6  2003/06/30 17:15:09  ebenhamou
 * remove a typo in the code. Previous user had ignored important warning message!
 *
 * Revision 1.5  2003/06/30 09:18:22  ebenhamou
 * remove unused var
 *
 *
 */


/*----------------------------------------------------------------------------*
    bsxtic.cpp
 
    Black and Scholes Analytics for Vanilla and Exotic Options
 
    Copyright (c) 1996
*----------------------------------------------------------------------------*/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "armglob.h"
#include "linalg.h"
#include "gaussian.h"
#include "newton.h"
#include "bsxtic.h"
#include "interpol.h"
#include "fromto.h"
#include "bsmodel.h"
#include "expt.h"

#include "gpbase/timer.h"
using ARM::ARM_Timer;



#define kMinSpot 1.0e-4
#define kMaxSpot 1.0e10    // limits for bracketing spot root
#define kSpotTol 1.0e-4
#define kIterMax 20

#define K_VAL_ERROR     -100


#define ARM_FX_GEN_VOL_PREC 0.0001




/*----------------------------------------------------------------------------*
   NAME    bsOption

   Compute value of european Vanilla option using Black & Scholes formula

   Output:
   the function returns the Black Scholes price

*----------------------------------------------------------------------------*/

double bsOption(double spot,
                double strike,
                double volatility,
                double dividend,
                double discountRate,
                double maturity,
                double CallPut)
{
	// on renvoie la valeur intrinsèque
    if ( maturity <= 1.e-6 )
	{
	   if ( CallPut == K_CALL ) // a CALL option 
		  return(MAX(spot-strike, 0.0));
	   else
		  return(MAX(strike-spot, 0.0));
	}

    if ( volatility <= 0.0 )
    {
       double df = exp(-discountRate*maturity);

 	   if ( CallPut == K_CALL ) // a CALL option 
		  return(MAX(spot-strike, 0.0)*df);
	   else
		  return(MAX(strike-spot, 0.0)*df);      
    }

    double vsq, vsqrtmat, value, d1, d2;

    vsq = SQR(volatility);

    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike)+
            (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

    d2 = d1-vsqrtmat;

    value = CallPut *(
            spot*exp(-dividend*maturity)*cdfNormal(CallPut*d1) -
            strike*exp(-discountRate*maturity)*cdfNormal(CallPut*d2) );

    return(value);
}


/*----------------------------------------------------------------------------*
   NAME    bsOptionSq

   Compute value of the square of an european Vanilla option using Black & 
   Scholes formula

   Output:
   the function returns the Black Scholes square price

*----------------------------------------------------------------------------*/

double bsOptionSq(double spot,
                double strike,
                double volatility,
                double dividend,
                double discountRate,
                double maturity,
                double CallPut)
{
    double vsq, vsqrtmat, value, d1, d2, d3;



    vsq = SQR(volatility);

    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike)+
            (discountRate - dividend + 1.5*vsq)*maturity) / vsqrtmat;

    d2 = d1-vsqrtmat;

    d3 = d2-vsqrtmat;

    value = spot*spot*exp(-2*dividend*maturity)*exp(vsq*maturity)*cdfNormal(CallPut*d1)
            -2*strike*spot*exp(-(dividend+discountRate)*maturity)*cdfNormal(CallPut*d2)
            +strike*strike*exp(-2*discountRate*maturity)*cdfNormal(CallPut*d3);

    return(value);
}


/*----------------------------------------------------------------------------*
    NAMES    bsGreek

        Compute the derivative of european Vanilla option
        price with respect to each argument of bsOption

*----------------------------------------------------------------------------*/

double bsDelta(double spot,
               double strike,
               double volatility,
               double dividend,
               double discountRate,
               double maturity,
               double CallPut)
{
    double   vsq, vsqrtmat, value, d1;


    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

    value = CallPut*exp(-dividend*maturity)*cdfNormal(CallPut*d1);

    return(value);
}



double bsKappa(double spot,
               double strike,
               double volatility,
               double dividend,
               double discountRate,
               double maturity,
               double CallPut)
{
    double   vsq, vsqrtmat, value, d1, d2;



    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;
    d2 = d1 - vsqrtmat;

    value = -CallPut * cdfNormal(CallPut * d2) * exp(-discountRate * maturity);

    return (value);
}



double bsGamma(double spot,
               double strike,
               double volatility,
               double dividend,
               double discountRate,
               double maturity,
               double CallPut)

{
    double   vsq, vsqrtmat, value, d1;

    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

    value = exp(-dividend*maturity)
            /volatility/spot/pow(maturity, 0.5)*dNormal(d1) / 100.0;

    return (value);
}



double bsVega(double spot,
              double strike,
              double volatility,
              double dividend,
              double discountRate,
              double maturity,
              double CallPut)
{
    double   vsq, vsqrtmat, value, d1;


    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);


	if ((fabs(strike ) > K_DOUBLE_TOL) && (fabs(vsqrtmat ) > K_DOUBLE_TOL))
	{
		d1 = (log(spot/strike) +
			(discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

		value = exp(-dividend*maturity)*spot*pow(maturity, 0.5)*dNormal(d1)/100.0;
	}
	else
		value = 0.0;

    return (value);
}



double bsRho(double spot,
             double strike,
             double volatility,
             double dividend,
             double discountRate,
             double maturity,
             double CallPut)

{
    double   vsq, vsqrtmat, value, d1, d2;

    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;
    d2 = d1 - vsqrtmat;
    value = CallPut * maturity * strike * exp(-maturity * discountRate)
             * cdfNormal(CallPut * d2) / 100000.0;

    return (value);
}



double bsTheta(double spot,
               double strike,
               double volatility,
               double dividend,
               double discountRate,
               double maturity,
               double CallPut)

{
    double   vsq, vsqrtmat, value, d1, d2;

    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

    d2 = d1 - vsqrtmat;

	value = (CallPut*dividend*spot*exp(-dividend*maturity)*cdfNormal(CallPut*d1)
			 -CallPut*discountRate*strike*exp(-discountRate*maturity)*cdfNormal(CallPut*d2)
			 -spot*dNormal(d1)*volatility*exp(-dividend*maturity)/2.0/pow(maturity, 0.5) 
			)/365.0;

    return (value);
}



/*----------------------------------------------------------------------------*
   NAME    dekBarrier

   Compute value of european Barrier option using Derman, Ergener & Kani
   static replication method.

    InOut = 1 if the option knocks in at barrier and 0 otherwise
    UpDown = 1 for Up, -1 for Down 
    numOptions is the number of options used to replicate the barrier option
*----------------------------------------------------------------------------*/

double dekBarrier(double spot,
                  double strike,
                  double barrier,
                  double volatility,
                  double dividend,
                  double discountRate,
                  double maturity,
                  double Rebate,
                  double CallPut,
                  double InOut,
                  double UpDown,
                  int numOptions)

{
    ARM_Vector *hdgOptMat, *hdgWeight;
    double brrOption=0.0, hdgOption=0.0, opt;
    int i, j;

    hdgOptMat = new ARM_Vector(numOptions+1, 0.0);
    hdgWeight = new ARM_Vector(numOptions+1, 0.0);

    for (i=1; i<=numOptions; i++)
    {
        (*hdgOptMat)[i] = double(numOptions+1-i)*maturity/numOptions;
    }
   
    for (i=1; i<=numOptions; i++)
    {
        if (InOut == 0.0)
        {
           // For Up and Out, hedge with Call euro options striked 
           // at the barrier For Down and Out, hedge with 
           // Put euro options striked at the barrier 

           hdgOption = bsOption(barrier, strike, volatility, 
                 dividend, discountRate, i*maturity/numOptions, CallPut);
            
           for (j = 1; j < i; j++) 
           {
               hdgOption += (*hdgWeight)[j]*bsOption(barrier, barrier, 
                               volatility, dividend, discountRate, 
                               (i+1-j)*maturity/numOptions, UpDown);
           }
        
           (*hdgWeight)[i] = (Rebate-hdgOption)/bsOption(barrier, barrier, 
                               volatility, dividend, discountRate, 
                               maturity/numOptions, UpDown);
        
           brrOption += (*hdgWeight)[i]*bsOption(spot, barrier, volatility, 
                        dividend, discountRate, 
                        (numOptions+1-i)*maturity/numOptions, UpDown);
        }
        else 
        {    
           // For Up and In,hedge with Call euro options striked at the barrier 
           // For Down and In, hedge with Put euro options 
           // striked at the barrier 

           hdgOption = Rebate*exp(-discountRate*double(i)*maturity/numOptions);
            
           opt = bsOption(barrier, strike, volatility, 
                dividend, discountRate, i*maturity/numOptions, CallPut);
            
           for (j=1; j<i; j++) 
           {
               hdgOption += (*hdgWeight)[j]*bsOption(barrier, barrier, 
                            volatility, 
                            dividend, discountRate, 
                            (i+1-j)*maturity/numOptions, UpDown);
           }
        
           (*hdgWeight)[i] = (opt-hdgOption)/bsOption(barrier, barrier, 
                                  volatility, 
                      dividend, discountRate, maturity/numOptions, UpDown);
                
            brrOption += (*hdgWeight)[i]*bsOption(spot, barrier, 
                         volatility, 
                         dividend, 
                  discountRate, (numOptions+1-i)*maturity/numOptions, UpDown);
        }
    }
    
    // For KnockOut, an European option with the same parameters 
    // as the Barrier option 
    // must be added to the hedge (cf DE&K)
    
    if ( InOut == 0.0 )
    {
       brrOption += bsOption(spot, strike, volatility, 
                     dividend, discountRate, maturity, CallPut);
    } 
    else 
    {
       brrOption += Rebate*exp(-discountRate*maturity);
    }
    
    if (hdgOptMat) 
       delete hdgOptMat;

    if (hdgWeight) 
       delete hdgWeight;

    return(brrOption);
}



/*----------------------------------------------------------------------------*
   NAME    bsBarrier

   Compute value of european Barrier option using Black & Scholes formula
   InOut = 1 if the option knocks in at barrier and 0 otherwise
*----------------------------------------------------------------------------*/

double bsBarrier(double spot,
                  double strike,
                  double barrier,
                  double volatility,
                  double dividend,
                  double discountRate,
                  double maturity,
                  double Rebate,
                  double CallPut,
                  double InOut,
                  double UpDown)

{
    double e,
        Variance,
        mu,
        L,
        d1,
        Integrale1,
        d2,
        Integrale2;

    double d3,
        Integrale3,
        d4,
        Integrale4,
        d5,
        Integrale5,
         Integrale6,
        value;

    double VarAux1,
         VarAux2,
        VarAux3,
        VarAux4;


    e = -(spot <= barrier) + (spot > barrier);

    Variance = volatility * volatility;

    mu = discountRate - dividend - 0.5 * Variance;

    L = 1 + mu / Variance;

    d1 = log(spot / strike) / sqrt(Variance * maturity)
        + L * sqrt(Variance * maturity);

    Integrale1 = CallPut * spot * exp(-dividend * maturity)
                    * cdfNormal(CallPut * d1)
                    - CallPut * strike * exp(-discountRate * maturity)
                    * cdfNormal(CallPut * (d1 - sqrt(Variance * maturity)));

    d2 = log(spot/barrier)/sqrt(Variance*maturity)+L*sqrt(Variance*maturity);

    Integrale2 = CallPut * spot * exp(-dividend * maturity)
                    * cdfNormal(CallPut * d2)
                    - CallPut * strike * exp(-discountRate * maturity)
                    * cdfNormal(CallPut * (d2 - sqrt(Variance * maturity)));

    d3 = log((barrier * barrier) / (spot * strike)) / sqrt(Variance * maturity)
        + L * sqrt(Variance * maturity);

    Integrale3 = CallPut * spot * exp(-dividend * maturity)
                    * pow(barrier / spot, 2 * L) * cdfNormal(e * d3)
                    - CallPut * strike * exp(-discountRate * maturity)
                    * pow(barrier / spot, 2.0 * L - 2)
                    * cdfNormal(e * (d3 - sqrt(Variance * maturity)));

    d4 = log(barrier / spot) / sqrt(Variance * maturity) 
         + L * sqrt(Variance * maturity);

    Integrale4 = CallPut * spot * exp(-dividend * maturity)
                    * pow(barrier / spot, 2 * L) * cdfNormal(e * d4)
                    - CallPut * strike * exp(-discountRate * maturity)
                    * pow(barrier / spot, 2.0 * L - 2)
                    * cdfNormal(e * (d4 - sqrt(Variance * maturity)));

    d5 = log(spot/barrier)/sqrt(Variance*maturity)+L*sqrt(Variance*maturity);

    Integrale5 = (cdfNormal(e * (d5 - volatility * sqrt(maturity)))
                    - pow(barrier / spot, 2.0 * L - 2)
                    * cdfNormal(e * (d4 - volatility * sqrt(maturity))))
                    * Rebate * exp(-discountRate * maturity);

    VarAux1 = mu / Variance;
    VarAux2 = pow(mu, 2.0) + 2 * discountRate * Variance;
    VarAux3 = sqrt(VarAux2) / Variance;
    VarAux4 = log(barrier / spot) / sqrt(Variance * maturity);
    VarAux4 += VarAux3 * sqrt(Variance * maturity);
    Integrale6 = Rebate * pow(barrier / spot, VarAux1 + VarAux3)
                    * cdfNormal(e * VarAux4);

    Integrale6 += Rebate*pow(barrier/spot, VarAux1-VarAux3)
                *cdfNormal(e*(VarAux4-2*VarAux3*sqrt(Variance*maturity)));


    if (InOut == 1.0 )
    {
        if (CallPut * strike >= CallPut * barrier)
        {
            value = Integrale1 * ((1 - CallPut * e) / 2) 
                + Integrale3 * ((1 + CallPut * e) / 2) + Integrale5;
        }
        else
        {
            value = (Integrale1 - Integrale2) * ((1 + CallPut * e) / 2) 
                    +(Integrale2 - Integrale3) * ((1 - CallPut * e) / 2)
                        + Integrale5 + Integrale4;
        }
    }
    else
    {
        if (CallPut * strike >= CallPut * barrier)
        {
            value = Integrale1 * ((1 + CallPut * e) / 2) -
                Integrale3 * ((1 + CallPut * e) / 2) + Integrale6;
        }
        else
        {
                value = (Integrale2 - Integrale4) * ((1 + CallPut * e) / 2) +
                (Integrale1 - Integrale2 + Integrale3 - Integrale4) *
                ((1 - CallPut * e) / 2) + Integrale6;
        }
    }

    return(value);
}



/*----------------------------------------------------------------------------*
   NAME    bsDoubleBarrier

   Compute value of european double Barrier option using Black&Scholes formula
   the option knock out when underlying price touch the up or low barrier
*----------------------------------------------------------------------------*/

double bsDoubleBarrier(double spot,
                        double strike,
                        double UpBarrier,
                        double DownBarrier,
                        double volatility,
                        double dividend,
                        double discountRate,
                        double maturity,
                        double CallPut)
{
    double res, Variance, Probabilite1, Probabilite2,
         PrixVariable1, PrixVariable2, Alpha1, Alpha2, Gamma1, Gamma2;

    double C1, C2, d11, d12,d21, d22, d31, d32, d41, d42;
    int n;

    Variance = pow(volatility, 2.0);

    Probabilite1 = 0.0;
    Probabilite2 = 0.0;

    if (CallPut == 1.0)
    {
        PrixVariable1 = strike;
        PrixVariable2 = UpBarrier;
    }
    else
    {
        PrixVariable1 = DownBarrier;
        PrixVariable2 = strike;
    }

    for (n = 5; n >= -5; n--)
    {
        C1 = ((2 * (discountRate - dividend)) / Variance) + 1;
        C2 = ((2 * (discountRate - dividend)) / Variance) - 1;
        Alpha1 = pow((UpBarrier / DownBarrier), n * C1);
        Alpha2 = pow((UpBarrier / DownBarrier), n * C2);
        Gamma1 = pow((DownBarrier / UpBarrier), n * C1)
              * pow((DownBarrier / spot), C1);

        Gamma2 = pow((DownBarrier / UpBarrier), n * C2)
              * pow((DownBarrier / spot), C2);

        d11 = (log(pow((UpBarrier / DownBarrier), 2.0 * n)
            * (spot / PrixVariable1))+ (discountRate - dividend
                        + 0.5 * Variance) * (maturity))
                              / sqrt(maturity * Variance);

        d12 = (log(pow((UpBarrier / DownBarrier), 2.0 * n)
                  * (spot / PrixVariable1)) 
                                  + (discountRate - dividend - 0.5 *
                         Variance) * (maturity)) / sqrt(maturity * Variance);

        d21 = (log(pow((UpBarrier / DownBarrier), 2.0 * n)
                              * (spot / PrixVariable2))
                  + (discountRate - dividend + 0.5 * Variance) 
                                  * (maturity)) / sqrt(maturity * Variance);

        d22 = (log(pow((UpBarrier / DownBarrier), 2.0 * n)
                               * (spot / PrixVariable2))
                  + (discountRate - dividend - 0.5 * Variance) 
                                * (maturity)) / sqrt(maturity * Variance);

        d31 = (log(pow((DownBarrier / UpBarrier), 2.0 * n)
                  * (DownBarrier * DownBarrier / spot / PrixVariable1))
                  + (discountRate - dividend + 0.5 * Variance) * (maturity))
                        / sqrt(maturity * Variance);

        d32 = (log(pow((DownBarrier / UpBarrier), 2.0 * n)
                 * (DownBarrier * DownBarrier / spot / PrixVariable1))
                  + (discountRate - dividend - 0.5 * Variance)
                  * (maturity)) / sqrt(maturity * Variance);

        d41 = (log(pow((DownBarrier / UpBarrier), 2.0 * n) 
              * (DownBarrier * DownBarrier
                    / spot / PrixVariable2)) 
                            + (discountRate - dividend + 0.5 *
            Variance) * (maturity)) / sqrt(maturity * Variance);

        d42 = (log(pow((DownBarrier / UpBarrier), 2.0 * n) *
              (DownBarrier * DownBarrier / spot / PrixVariable2)) +
              (discountRate - dividend - 0.5 *
            Variance) * (maturity)) / sqrt(maturity * Variance);

        Probabilite1 += Alpha1 * (cdfNormal(d11) - cdfNormal(d21)) - Gamma1 *
            (cdfNormal(d31) - cdfNormal(d41));

        Probabilite2 += Alpha2 * (cdfNormal(d12) - cdfNormal(d22)) - Gamma2 
                                * (cdfNormal(d32) - cdfNormal(d42));
    }

    res = (Probabilite1 * spot * exp(-dividend * maturity) - Probabilite2
             * strike * exp(-discountRate * maturity)) * CallPut;

    return (res);
}




/*----------------------------------------------------------------------------*
    NAME    bsBinary

    Compute value of european Binary option using Black & Scholes formula

    If CashAsset is set to Cash, a binary call(put) pays off 1 unit if the spot
    finishes above(below) the striking price and nothing elsewhere 
    When set to Asset, the call pays off the underlying asset price at maturity  
*----------------------------------------------------------------------------*/

double bsBinary(double spot,
                 double strike,
                 double volatility,
                 double dividend,
                 double discountRate,
                 double maturity,
                 double CallPut,
                 double CashAsset)
{
    double   vsq, vsqrtmat, value, d1, d2;

    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate-dividend+0.5*vsq)*maturity) / vsqrtmat;

        d2 = d1 - vsqrtmat;

    value = CashAsset*
            exp(-discountRate*maturity)*cdfNormal(CallPut*d2)+
            (1.0-CashAsset)*
            spot*exp(-dividend*maturity)*cdfNormal(CallPut*d1);

    return(value);
}




/*----------------------------------------------------------------------------*
    NAME    bsGap

    Compute value of european Gap option using Black & Scholes formula

    A gap Option is like vanilla option except that the payoff is
    computed using a price different from strike: for example, a call pays
    off Max(S-X, 0) if S>K (if the payoff argument is set to X)

*----------------------------------------------------------------------------*/

double bsGap(double spot,
             double strike,
             double payoff,
             double volatility,
             double dividend,
             double discountRate,
             double maturity,
             double CallPut)
{
    double   vsq, vsqrtmat, value, d1, d2;

    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

        d2 = d1 - vsqrtmat;

    value = CallPut *(
            spot*exp(-dividend*maturity)*cdfNormal(CallPut*d1) -
            payoff*exp(-discountRate*maturity)*cdfNormal(CallPut*d2) );

    return(value);
}



/*----------------------------------------------------------------------------*
    NAME    bsSuperShare

    Compute value of european Gap option using Black & Scholes formula

    A superShare pays at expiry spot/infStrike if infStrike<=spot<supStrike
    and expires worthless otherwise
*----------------------------------------------------------------------------*/

double bsSuperShare(double spot,
                     double strikeInf,
                     double strikeSup,
                     double volatility,
                     double dividend,
                     double discountRate,
                     double maturity)
{
    double vsq, vsqrtmat, value, dInf, dSup;

    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    dInf = (log(spot/strikeInf) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;
    dSup = (log(spot/strikeSup) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;

    value = spot*exp(-dividend*maturity)*
            (cdfNormal(dInf) - cdfNormal(dSup))/strikeInf;

    return(value);
}



/*----------------------------------------------------------------------------*
    NAME    bsForwardStart

    Compute value of european Forward Start option using Black&Scholes formula

    After startDate, the buyer will receive an option with time to expiration
    "maturity" and strike equal to IOMoney*Spot price(uncertain) at startDate
    that is, the strike is IOMoney percent of the asset price at startDate 

 !!! Note that "maturity" is the maturity of the option today, not at Startdate
*----------------------------------------------------------------------------*/

double bsForwardStart(double spot,
                         double IOMoney,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double startDate,
                         double maturity,
                         double CallPut)
{
    double   vsq, vsqrtmat, value, d1, d2, optMat;



    optMat = maturity - startDate;

    if (fabs(optMat) < K_DOUBLE_TOL)
    {
       value = exp(-dividend*startDate)*MAX(CallPut*spot*(1-IOMoney), 0.0);

       return(value);
    }
    else
    {
        vsq = SQR(volatility);
        vsqrtmat = volatility * sqrt(optMat);

        d1 = ( -log(IOMoney) +
            (discountRate - dividend + 0.5*vsq)*optMat ) / vsqrtmat;

        d2 = d1 - vsqrtmat;

        value = CallPut*exp(-dividend*startDate)*spot*(
                exp(-dividend*optMat)*cdfNormal(CallPut*d1) -
                IOMoney*exp(-discountRate*optMat)*cdfNormal(CallPut*d2) );

        return(value);
    }
}



/*----------------------------------------------------------------------------*

    NAME    bsPayLater

    Compute value of european Pay Later option using Black & Scholes formula

    The buyer of the option pays a premium which must be paid if the option
    expires in the money at maturity, regardless of whether the option expires
    sufficiently deeply in the money

*----------------------------------------------------------------------------*/

double bsPayLater(double spot,
                     double strike,
                     double volatility,
                     double dividend,
                     double discountRate,
                     double maturity,
                     double CallPut)
{
    double   value, vsq, vsqrtmat, d1, d2;



    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    d1 = (log(spot/strike) +
        (discountRate - dividend + 0.5*vsq)*maturity) / vsqrtmat;
        d2 = d1 - vsqrtmat;

    value = CallPut*(
            spot*exp((discountRate-dividend)*maturity)*
            cdfNormal(CallPut*d1)/cdfNormal(CallPut*d2) - strike );

    return(value);
}



/*----------------------------------------------------------------------------*
    NAME    bsLookBack

    Compute value of european Look Back option using Black & Scholes formula

    The look back call option pays off:
        max(0, Sn-min(S0, S1,..., Sn)) = Sn-min(S0, S1,..., Sn)
    The look back put option pays off:
        max(0, max(S0, S1,..., Sn)-Sn) = max(S0, S1,..., Sn)-Sn)

  currentMinOrMax is the current minimum or maximum that the asset price has
    attained since the inception of the option
*----------------------------------------------------------------------------*/

double bsLookBack(double spot,
                     double currentMinOrMax,
                     double volatility,
                     double dividend,
                     double discountRate,
                     double maturity,
                     double CallPut)
{
    double   value, vsq, vsqrtmat, b, mu, lambda;


    vsq = SQR(volatility);
    vsqrtmat = volatility * sqrt(maturity);

    b = log(spot/currentMinOrMax);
    mu = discountRate - dividend - 0.5*vsq;
    lambda =    0.5*vsq/(discountRate-dividend);

    value = CallPut*(
        spot*exp(-dividend*maturity)-
        currentMinOrMax*exp(-discountRate*maturity)*
            cdfNormal(CallPut*(b+mu*maturity)/vsqrtmat)+
        currentMinOrMax*exp(-discountRate*maturity+b*(1.0-1.0/lambda))*
            lambda*cdfNormal(CallPut*(-b+mu*maturity)/vsqrtmat)-
        spot*exp(-dividend*maturity)*(1.0+lambda)*
            cdfNormal(CallPut*(-b-mu*maturity-vsq*maturity)/vsqrtmat) );

    return(value);
}



/*----------------------------------------------------------------------------*
  NAME    bsAveragePrice

  Compute value of european Average Price option using Black & Scholes formula

        The  Average Price call option pays off:
            max(0, (S1+...+Sn)/n-K)
        The  Average Price put option pays off:
            max(0, K-(S1+...+Sn)/(n)
     the valuation is from M. Curran (Risk-1996) 

*----------------------------------------------------------------------------*/

double bsAveragePrice(double spot,
                         double strike,
                         double numAvrgPoints,
                         double firstAvrgPoints,
                         double dt,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double maturity,
                         double CallPut)
{
     int i;
    double   value, vsq, var, mu, mat_i, vsqmat_i, var_i, mu_i, strikeHat;


    vsq = SQR(volatility);

    mu = log(spot) + (discountRate - dividend - 0.5*vsq)*
        (firstAvrgPoints+0.5*(numAvrgPoints-1)*dt);
    var = vsq*( firstAvrgPoints+
        dt*(numAvrgPoints-1)*(2*numAvrgPoints-1)/(6*numAvrgPoints) );

    strikeHat = 0.0;

    for (i=1; i<=numAvrgPoints; i++)
    {
        mat_i = firstAvrgPoints+(i-1.0)*dt;
        mu_i = log(spot) + (discountRate - dividend - 0.5*vsq)*mat_i;
        vsqmat_i = SQR(volatility)*mat_i;
        var_i = vsq*( firstAvrgPoints+
                    dt*( (i-1)-i*(i-1)/(2*numAvrgPoints) ) );

        strikeHat += exp( mu_i+var_i*(log(strike)-mu)/var +
                  0.5*(vsqmat_i-var_i/var) );
    }
    strikeHat /= numAvrgPoints;
    strikeHat = 2*strike - strikeHat;

    value = 0.0;
    for (i=1; i<=numAvrgPoints;i++)
    {
        mat_i = firstAvrgPoints+(i-1.0)*dt;
        mu_i = log(spot) + (discountRate - dividend - 0.5*vsq)*mat_i;
        vsqmat_i = SQR(volatility)*mat_i;
        var_i = vsq*( firstAvrgPoints+
                    dt*( (i-1)-i*(i-1)/(2*numAvrgPoints) ) );
        value += exp(mu_i+0.5*vsqmat_i)*
            cdfNormal( (mu-log(strikeHat))/sqrt(var)+sqrt(var_i/var) );
    }
    value /= numAvrgPoints;
    value -= strike*cdfNormal((mu-log(strikeHat))/sqrt(var));
    value *= exp(-discountRate*maturity);

    return(value);
}



/*----------------------------------------------------------------------------*
    NAME    bsChooser

    Compute value of european Chooser option using Black & Scholes formula
    A chooser Option gives the buyer to choose the type of option
    (either call or put) at chsCallPutDate.
    The Call strike and maturity may be different from those of the Put
*----------------------------------------------------------------------------*/

double bsChooser(double spot,
                 double callStrike,
                 double putStrike,
                 double volatility,
                 double dividend,
                 double discountRate,
                 double callMaturity,
                 double putMaturity,
                 double chsCallPutDate)
{
    void* fixedParams[7];

    double  spotCP, spotCP1, spotCP2, vsq, vsqrtmat, value, x1, x2,
            tCall, tPut, x1Call, x1Put, x2Call, x2Put;

    if (chsCallPutDate <= K_DOUBLE_TOL)
    {
        x1Call = bsOption(spot, callStrike, volatility, dividend,
            discountRate, callMaturity, 1.0);
        x1Put = bsOption(spot, putStrike, volatility, dividend,
            discountRate, putMaturity, -1.0);

        return(MAX(x1Call, x1Put));
    }
    else if(fabs(callMaturity - chsCallPutDate) <= K_DOUBLE_TOL)
    {
        x1Call = MAX(spot-callStrike, 0.0);

        x1Put = bsOption(spot, putStrike, volatility, dividend,
            discountRate, putMaturity, -1.0);

        return(MAX(x1Call, x1Put));
    }
    else if(fabs(putMaturity - chsCallPutDate) <= K_DOUBLE_TOL)
    {
        x1Call = bsOption(spot, callStrike, volatility, dividend,
            discountRate, callMaturity, 1.0);
        x1Put = MAX(putStrike-spot, 0.0);

        return(MAX(x1Call, x1Put));
    }
    else
    {
        tCall = callMaturity - chsCallPutDate;
        tPut = putMaturity - chsCallPutDate;

        fixedParams[0] = (void *) &callStrike;
        fixedParams[1] = (void *) &putStrike;
        fixedParams[2] = (void *) &tCall;
        fixedParams[3] = (void *) &tPut;
        fixedParams[4] = (void *) &volatility;
        fixedParams[5] = (void *) &dividend;
        fixedParams[6] = (void *) &discountRate;

        // Find the spot price such that Call=Put using <newtonRoot>

        spotCP1 = callStrike;

        while (     bsOption(spotCP1, callStrike, volatility, dividend,
                discountRate, tCall, 1.0) -
                bsOption(spotCP1, putStrike, volatility, dividend,
                discountRate, tPut, -1.0) < 0.0 && spotCP1 > kMinSpot) 
        {
            spotCP1 *= 1.5;
        }

        spotCP2 = putStrike;

        while (    bsOption(spotCP2, callStrike, volatility, dividend,
                discountRate, tCall, 1.0) -
                bsOption(putStrike, spotCP2, volatility, dividend,
                discountRate, tPut, -1.0) > 0.0 && spotCP2 < kMaxSpot)
        {
             spotCP2 /= 1.5;
        }

        // !!! erreur a gerer

        spotCP = newtonRoot( bsSpot_CallPutParity, spotCP1, spotCP2,
                        fixedParams, kSpotTol, kIterMax);

        vsq = SQR(volatility);
        vsqrtmat = volatility * sqrt(chsCallPutDate);

        x1 = (log(spot/spotCP) +
            (discountRate - dividend + 0.5*vsq)*chsCallPutDate) / vsqrtmat;
        x2 = x1 - vsqrtmat;

        vsqrtmat = volatility * sqrt(callMaturity);

        x1Call =     (log(spot/callStrike) +
                (discountRate - dividend + 0.5*vsq)*callMaturity) / vsqrtmat;
        x2Call = x1Call - vsqrtmat;

        vsqrtmat = volatility * sqrt(putMaturity);

        x1Put =     (log(spot/putStrike) +
                (discountRate - dividend + 0.5*vsq)*putMaturity) / vsqrtmat;
        x2Put = x1Put - vsqrtmat;

        tCall = sqrt(chsCallPutDate/callMaturity);
        tPut = sqrt(chsCallPutDate/putMaturity);

        value = spot*exp(-dividend*callMaturity)
                *cdf2Normal(tCall, x1, x1Call)-
                callStrike*exp(-discountRate*callMaturity)
                *cdf2Normal(tCall, x2, x2Call) -

        spot*exp(-dividend*putMaturity)*cdf2Normal(tPut, -x1, -x1Put) +
        putStrike*exp(-discountRate*putMaturity)*cdf2Normal(tPut, -x2, -x2Put);

        return(value);
    }
}



double bsSpot_CallPutParity(double& f, double& df, 
                            double spot, void** fixedParams)
{
    double callStrike, putStrike, callMaturity, putMaturity, volatility,
    dividend, discountRate;

    callStrike = * (double *) fixedParams[0];
    putStrike = * (double *) fixedParams[1];
    callMaturity = * (double *) fixedParams[2];
    putMaturity = * (double *) fixedParams[3];
    volatility = * (double *) fixedParams[4];
    dividend = * (double *) fixedParams[5];
    discountRate = * (double *) fixedParams[6];

    f = bsOption(spot, callStrike, volatility, dividend,
        discountRate, callMaturity, 1.0) -
            bsOption(spot, putStrike, volatility, dividend,
        discountRate, putMaturity, -1.0);
    df = bsDelta(spot, callStrike, volatility, dividend,
        discountRate, callMaturity, 1.0) -
            bsDelta(spot, putStrike, volatility, dividend,
        discountRate, putMaturity, -1.0);

    return(f);
}



/*----------------------------------------------------------------------------*
    NAME    bsCompound

    Compute value of european Compound option using Black & Scholes formula
        
    Undelying option maturity must be less than compound option maturity

*----------------------------------------------------------------------------*/

double bsCompound(double spot,
                     double volatility,
                     double dividend,
                     double discountRate,
                     double uoStrike,
                     double uoMaturity,
                     double uoCallPut,
                     double oStrike,
                     double oMaturity,
                     double oCallPut)
{
    void *fixedParams[7];
    double spot1, spot2, spotStrike, vsq, oVsqrtmat, uoVsqrtmat, value,
    x1, x2, y1, y2, c1, c2, c3, r;


    fixedParams[0] = (void *) &oStrike;
    fixedParams[1] = (void *) &uoStrike;
    fixedParams[2] = (void *) &volatility;
    fixedParams[3] = (void *) &dividend;
    fixedParams[4] = (void *) &discountRate;
    fixedParams[5] = (void *) &uoMaturity;
    fixedParams[6] = (void *) &uoCallPut;

    //    Find the critical spot price
    //    ie, such that option value=oStrike using <newtonRoot>


    spot1 = spot;
    while ( bsOption(spot1, uoStrike, volatility, dividend,
              discountRate, uoMaturity, uoCallPut) 
               - oStrike < 0.0) spot1 += uoCallPut*0.5*spot1;

    spot2 = spot;

    while ( bsOption(spot2, uoStrike, volatility, dividend,
                     discountRate, uoMaturity, uoCallPut) - oStrike > 0.0) 
    {
        spot2 += -uoCallPut*0.5*spot2;
    }
    
    spotStrike = newtonRoot(bsSpot_OptionPrice, spot1, spot2,
        fixedParams, kSpotTol, kIterMax);

    vsq = SQR(volatility);
    oVsqrtmat = volatility * sqrt(oMaturity);
    uoVsqrtmat = volatility * sqrt(uoMaturity);

    x1 = uoCallPut*oCallPut*(
        log(spot/spotStrike) +
        (discountRate - dividend + 0.5*vsq)*oMaturity) / oVsqrtmat;

    x2 = x1 - uoCallPut*oCallPut*oVsqrtmat;

    y1 = uoCallPut*( log(spot/uoStrike) +
        (discountRate - dividend + 0.5*vsq)*uoMaturity) / uoVsqrtmat;
    y2 = y1 - uoCallPut*uoVsqrtmat;

    r = sqrt(oMaturity/uoMaturity);

    c1 = uoCallPut*oCallPut*spot*exp(-dividend*uoMaturity)
         *cdf2Normal(oCallPut*r, x1, y1);

    c2 = uoCallPut*oCallPut*uoStrike*exp(-uoMaturity*discountRate)
          *cdf2Normal(oCallPut*r, x2, y2);

    c3 = oCallPut*oStrike*exp(-oMaturity*discountRate)*cdfNormal(x2);

    value = c1 - c2 - c3;

    return(value);
}



double bsSpot_OptionPrice(double& f, double& df, double spot, 
                          void** fixedParams)
{
    double optionValue, strike, maturity, volatility,
    dividend, discountRate, CallPut;

    optionValue = * (double *) fixedParams[0];
    strike = * (double *) fixedParams[1];
    volatility = * (double *) fixedParams[2];
    dividend = * (double *) fixedParams[3];
    discountRate = * (double *) fixedParams[4];
    maturity = * (double *) fixedParams[5];
    CallPut = * (double *) fixedParams[6];

     f = bsOption(spot, strike, volatility, dividend,
        discountRate, maturity, CallPut) - optionValue;
    df = bsDelta(spot, strike, volatility, dividend,
        discountRate, maturity, CallPut);

    return(f);
}




/*----------------------------------------------------------------------*
          FONCTIONS IMPLICITES
*----------------------------------------------------------------------*/
/* Calcul de la volatilite implicite nouveau algo */

double bsVolImp2(double Price,
                double spot,
                double strike,
                double dividend,
                double discountRate,
                double mat,
                double CallPut,
                double Ytol,
				double InitVal, // Init Value = 0.5
				int algo) 
{
	double vol1 =0.01, vol2 =1., vol, d_Vol, oldVol, bsPrice, d_bsPrice, oldBsPrice;
	try
	{
		double bsPrice1 = Price - bsOption(spot, strike, vol1, dividend, discountRate, mat, CallPut);
		double bsPrice2 = Price - bsOption(spot, strike, vol2, dividend, discountRate, mat, CallPut);

		int sig, nbIter = 0;
		while( (bsPrice1*bsPrice2 > 0.0) && (nbIter < 200))
		{
			if(fabs(bsPrice1) < fabs(bsPrice2))
			{
				vol1 += 1.6*(vol1-vol2);
				bsPrice1 = Price - bsOption(spot, strike, vol1, dividend, discountRate, mat, CallPut);
			}
			else
			{
				vol2 += 1.6*(vol2-vol1);
				bsPrice2 = Price - bsOption(spot, strike, vol2, dividend, discountRate, mat, CallPut);
			}
			nbIter++;
		}
		vol = (vol1+vol2)/2.0;
		double EPS_VOL = 1.0e-15;
		nbIter = 0;
		if(fabs(bsPrice1) < fabs(bsPrice2))
		{
			
			if(bsPrice1 - Price != 0)
			{
				oldBsPrice = bsPrice1;
				oldVol = vol1;
			}
			else
			{
				oldBsPrice = bsPrice2;
				oldVol = vol2;
			}
			
		}
		else
		{
			if(bsPrice2 - Price != 0)
			{
				oldBsPrice = bsPrice2;
				oldVol = vol2;
			}
			else
			{
				oldBsPrice = bsPrice1;
				oldVol = vol1;
			}
		}

		do
		{
			bsPrice = Price - bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut);
			d_bsPrice = bsVega(spot, strike, vol, dividend, discountRate, mat, CallPut);

			if(fabs(d_bsPrice)<EPS_VOL)
			{
				if(fabs(bsPrice)<EPS_VOL)
					break;
				if(bsPrice - Price == 0.0) // correction de l'anomalie de cfdNormal...
				{
					vol = oldVol;
					double pas = (CallPut == 1 ? 10. : 0.1);
					int i = 0;
					while((bsPrice - Price == 0.0)&& (i < 100))
					{

						bsPrice = Price - bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut);
						double h = 0.01*MAX(fabs(vol),0.01);
				
						double h_Price = Price - bsOption(spot, strike, vol-h, dividend, discountRate, mat, CallPut);
						d_bsPrice = (bsPrice - h_Price)/h;
						vol *= pas;
						i++;
					}
					if(i>99)
						return -8888;
					vol /= pas;
					d_Vol = - (bsPrice/d_bsPrice)/100.;
					sig = (d_Vol > 0 ? 1 : -1);
					d_Vol = sig*MIN(fabs(d_Vol),(oldVol == vol ? 1 : fabs(oldVol-vol)));
					oldVol = vol;
					oldBsPrice = bsPrice;
					if(fabs(d_Vol) < 1.0e-6)
						break;
				}
				else
				{
					vol = 0.;
					break;
				}
			}
			else if(fabs(oldBsPrice) < fabs(bsPrice))
				vol = (vol+oldVol)/2.0;
			else
			{
				d_Vol = - (bsPrice/d_bsPrice)/100.;
				sig = (d_Vol > 0 ? 1 : -1);
				d_Vol = sig*MIN(fabs(d_Vol),(oldVol == vol ? 1 : fabs(oldVol-vol)));
				oldVol = vol;
				oldBsPrice = bsPrice;
				if(fabs(d_Vol) < 1.0e-6)
					break;
				else if(vol-d_Vol < 0.75*vol)
					vol = 0.75*vol;
				else if(vol-d_Vol > 1.25*vol)
					vol = 1.25*vol;
				else
					vol -= d_Vol;
			}
			nbIter++;
		}
		while(nbIter < 100);

		if(nbIter > 99)
			vol = -8888;
	}
	catch(Exception& m)
	{
		m.DebugPrint();
		throw m;
	}
	return vol;
}


/* Calcul de la Volatilite implicite   */

double bsVolImp1(double Price,
                double spot,
                double strike,
                double dividend,
                double discountRate,
                double mat,
                double CallPut,
                double Ytol,
				double InitVal, // Init Value = 0.5
				int algo) 
{
    double calcVega;

	double vol = InitVal, fc;

    int n = 0;

    if (( CallPut == 1.0 ) || ( CallPut == -1.0 ))
    {
       double fc0 = bsOption(spot, strike, K_DOUBLE_TOL, dividend, discountRate, mat, CallPut);

       if ( Price - fc0 < Ytol)
       {
          return(0.0);
       }

       fc = bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut);

       // Dichtomie

       if ( Price < fc )
       {
           double sup = vol;
           double inf = K_DOUBLE_TOL;
           double mid = (sup);

           while (( fabs(fc-Price) > 0.05 ) && ( n < 100 ))
           {
               mid = (sup+inf)/2.0;

               fc = bsOption(spot, strike, mid,  dividend, discountRate, mat,CallPut);

               if ( fc > Price )
                  sup = mid;
               else
                  inf = mid;

               n++;
           }

           vol = mid;
       }

       n = 0;

	   double VOL = vol * 100.;
 
	   int stop = 0;

       while (( fabs(fc-Price) > Ytol ) && ( n < 100 ) && ( stop == 0 ))
       {
           calcVega = bsVega(spot, strike, VOL/100., dividend, discountRate,mat, CallPut);

           if ( fabs(calcVega) <= Ytol )
           {				 
			   stop = 1;
           }
           else
		   {
			   if ( algo == 0 )
			   {
				  VOL += (Price - fc) / calcVega;
			   }
			   else
			   {
			      double step = (Price - fc)/calcVega;

			      if (VOL + step < 0.75*VOL)
				     VOL *= 0.75;
			      else if (VOL + step > 1.25*VOL)
				     VOL *= 1.25;
			      else
				     VOL += step;
               }
		   }

           fc = bsOption(spot, strike, VOL/100., dividend, discountRate, mat,CallPut);

           n++;
       }

	   vol = VOL/100.;

       if (( fabs(fc-Price) > Ytol ) || ( stop == 1))
	   {
		  if ( algo == 0 )
		  {
		     algo = 1;

			 // try an other algo

			 vol = bsVolImp(Price,
                            spot,
                            strike,
                            dividend,
                            discountRate,
                            mat,
                            CallPut,
                            Ytol,
							InitVal,
				            1 /* algo 1 */
							);

			 return(vol);
          }
          else
		  {
             vol = K_VAL_ERROR;

			 return(vol);
		  }
	   }
    }
    else
	{
       vol = K_VAL_ERROR;

	   return(vol);
	}

    return(vol);
}

double bsVolImp(double Price,
                double spot,
                double strike,
                double dividend,
                double discountRate,
                double mat,
                double CallPut,
                double Ytol,
				double InitVal, 
				int algo) 
{
	double vol = bsVolImp1(Price, spot, strike, dividend, discountRate, mat, CallPut, Ytol, InitVal, algo);
	double bsPrice = bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut);
	if((vol == 0.0) || (vol == -100.0))
		vol = bsVolImp2(Price, spot, strike, dividend, discountRate, mat, CallPut, Ytol, InitVal, algo);
	return vol;
}

/* Calcul de la Volatilite implicite Normal */

double bsVolImpNor(double Price,
                double spot,
                double strike,
                double dividend,
                double discountRate,
                double mat,
                double CallPut,
                double Ytol,
				ARM_BSNorModel* model)
{
    double calcVega;
	ARM_BSNorModel* mymodel = model;
    double vol = 0.0005, fc;
    int n = 0;
	
	if ( mymodel == NULL )
	   return(K_VAL_ERROR);
 
    if (( CallPut == 1.0 ) || ( CallPut == -1.0 ))
    {
		double fc0 = mymodel->bsOption(spot, strike, 0.00001, dividend, discountRate, mat, CallPut);

		
	   if ( Price < fc0 )
       {
          return(0.0);
       }

       fc = mymodel->bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut);
	
       // Dichtomie

	   double sup = vol;
       double inf = 0.00001;//K_DOUBLE_TOL;
       double mid = (sup);

       if ( Price < fc )
       {
           

           while (( fabs(fc-Price) > 0.001 ) && ( n < 200 ))
           {
               mid = (sup+inf)/2.0;

               fc = mymodel->bsOption(spot, strike, mid,  dividend, discountRate, mat,CallPut);
				   
               if ( fc > Price )
                  sup = mid;
               else
                  inf = mid;

               n++;
           }

           vol = mid;
       }

       n = 0;
	   double tempVol = 0.0;
	  // inf = ( fc > Price ) ?  inf : sup ;

       while (( fabs(fc-Price) > Ytol ) && ( n < 200 ))
       {
			calcVega = mymodel->bsVega(spot, strike, vol, dividend, discountRate,mat, CallPut);
   

			if ( calcVega == 0.0 )
			{

				if (gTrace)
					printf("\n ===> ????? Vega = 0.0 \n");
			}

			tempVol = vol + (Price - fc) / calcVega / 1000000.0;
			if( tempVol < inf )
			{
				vol *= 0.75;
				if( inf > vol )
					vol = inf;
			}
			else if ( tempVol > sup )
			{
				vol *= 1.25;
				if( sup < vol )
					vol = sup;
			}
			else
				vol = tempVol;
       
		   //vol += (Price-fc)/calcVega;

           fc = mymodel->bsOption(spot, strike, vol, dividend, discountRate, mat,CallPut);
			   

		   n++;
       }

       if (fabs(fc - Price) > Ytol)
          vol = K_VAL_ERROR;
    }
    else
       vol = K_VAL_ERROR;

    return (vol);
}

/* Calcul du strike implicite a partir du delta  */


double bsKImpDelta(double delta, 
                   double  spot, 
                   double  vol, 
                   double  dividend, 
                   double  discountRate, 
                   double  mat, 
                   double  CallPut, 
                   double  Ytol)
{
    double strike,
        d1;
 
    if ((CallPut == 1.0) || (CallPut == -1.0))
    {
        d1 = INV_PART_FUNC_NOR(delta / CallPut / exp(-dividend * mat))
                                               / CallPut;

        strike = spot * exp((discountRate - dividend + pow(vol, 2.0) / 2.0)
                                           * mat - d1 * sqrt(mat) * vol);
    }
    else
        strike = K_VAL_ERROR;
 
    return (strike);
}
 


/* Calcul du strike implicite a partir du delta  */


double bsKImpDelta2(double delta, 
                    double spot, 
                    double vol, 
                    double dividend, 
                    double discountRate, 
                    double mat, 
                    double CallPut, 
                    double Ytol)
{
    double strike, dc, d2;
    int n = 0;
 
    if ((CallPut == 1.0) || (CallPut == -1.0))
    {
       strike = bsKImpDelta(delta, spot, vol, dividend, discountRate, 
                                                 mat, CallPut, Ytol);

       dc = bsDelta(spot, strike, vol, dividend, discountRate, mat, CallPut)
            - bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut)
                                             / spot;

       while ((fabs(dc - delta) > Ytol) && (n <= 300))
       {
           d2 = (log(spot/strike)+(discountRate-dividend-vol*vol/2.0)
                          *mat)/vol/pow(mat, 0.5);

           strike += (delta - dc) / (exp(-discountRate * mat) / spot *
                      (cdfNormal(CallPut * d2) - dNormal(CallPut * d2) /
                                                 vol / sqrt(mat)));
 
           dc = bsDelta(spot,strike, vol, dividend,discountRate, mat, CallPut)
                 -bsOption(spot, strike, vol, dividend, discountRate, mat,
                                CallPut) / spot;
 
           n++;
       }
    }
    else
        strike = K_VAL_ERROR;
 
    return (strike);
}
 


/* Calcule du strike implicite  */

 
double bsKImp(double Price, 
                double spot, 
                double vol, 
                double dividend, 
                double discountRate, 
                double mat, 
                double CallPut, 
                double Ytol)
{
    double strike,
        fc;
    int n = 0;


 
    if ((CallPut == 1.0) || (CallPut == -1.0))
    {
       if (CallPut == 1.0)
          strike = 0.01;
 
       if (CallPut == -1.0)
          strike = 100000;
 
 
       fc = bsOption(spot, strike, vol, dividend, discountRate, mat, CallPut);

       while ((fabs(fc - Price) > Ytol) && (n <= 300))
       {
           strike -= (fc - Price) / bsKappa(spot, strike, vol, dividend,
                                  discountRate, mat, CallPut);

           fc = bsOption(spot, strike, vol, dividend, discountRate, 
                         mat, CallPut);
           n++;
       }

       if ( fabs(fc - Price) > Ytol )
          strike = K_VAL_ERROR;
    }
    else
        strike = K_VAL_ERROR;
 
    return (strike);
}
 


/* Forward BS Value of a Call */
double CallForwardValue(double F,
                        double K,
                        double Sigma,
                        double T)
{
    // If null strike, then return the forward:
    if ( K < 1.0e-6 )
       return F;

    // If negative forward, then return NULL:
    if ( F < 0.0 )
       return 0.0;


    double Q;

    // If in the past or present, or if null vol,
    // then return the intrinsic:

    if (( T < 1.0/365.0 ) || ( Sigma < 1.0e-6 ))
    {
       Q = MAX(F-K, 0.0);
    }
    else // Otherwise, price the option:
    {
       double SigRootT = Sigma*sqrt(T);
       double d1 = log(F/K)/SigRootT+0.5*SigRootT;
       double d2 = d1-SigRootT;

       Q = F*cdfNormal(d1)-K*cdfNormal(d2); // call's FV
    }

    return(Q);
}



/* Forward BS Value of a put: */
double PutForwardValue(double F,
                       double K,
                       double Sigma,
                       double T)
{
    // Evaluate the call:
    double Q = CallForwardValue(F, K, Sigma, T);

    // Call / Put parity:
    return(Q-(F-K));
}



double ComputeFwdBSDelta(double fwd,
                         double strike, 
                         double vol,
                         double T,
                         int CallPut)
{
    double delta ;

    if ( T <= 0.0 )
    {
       delta = 0.0;
    }
    else
    {
       // Optimization

       double vsqrt = vol*sqrt(T);

       delta = CallPut*cdfNormal(CallPut*(log(fwd/strike)/vsqrt+0.5*vsqrt));
    }

    return(delta);
} 



extern double FXCallFwd(double F, double K, double Sigma, double T);

extern double FXPutFwd(double F, double K, double Sigma, double T);



double ComputeSplinedSigmaATMF(ARM_Vector* deltas,
                               ARM_Vector* sigmas,
                               double matu,
                               double SigmaZDS,
                               double Precision,
                               double FX_SPOT)
{
    double sigma, prevSigma;

    int nbIterMax = 100; 
    int nbIter;

    double delta;
    double Forward;
    double Strike;
     
    int callPut = K_CALL;

    double sigmaATMF; // Init.

    if ( SigmaZDS < 0.0 )
    {
       sigmaATMF = 50.0;
    }
    else
    {
       sigmaATMF = SigmaZDS;
    }

    // Set the ATM constraint

    Forward = Strike = 100.0;

    sigma   = sigmaATMF/100.0;

    nbIter  = 0;

    do
    {
        prevSigma = sigma;

        // Delta Fwd
        delta = ComputeFwdBSDelta(Forward, Strike, sigma, matu, callPut)*100.0;
   
        if ( FX_SPOT >= 0.0 ) // Case Delta with premium
        {
           double premium =  FXCallFwd(Forward, Strike, sigma, matu)*100.0;
       
           delta = delta-(premium/FX_SPOT);
        }

        sigma = SplineInterpolateFunc(deltas, sigmas,
                                      delta,
                                      NULL, // SecondDerivCalc,
                                      1);   // keep2Der

        sigma /= 100.0;

        nbIter++;
    }
    while (( fabs(sigma-prevSigma) > Precision )
           &&
           ( nbIter < nbIterMax )
          );

    return(sigma*100.0);
}




double CalculateImpliedStrikeFromDeltaWithPremium(ARM_Date& AsOf,
                                                  double matu,
                                                  double sigma,
                                                  double fxSpot,
                                                  double deltaWithPremium,
                                                  ARM_ZeroCurve* domCrv, // JPY
                                                  ARM_ZeroCurve* foreignCrv) // USD
{
    int callPut = K_CALL;
    double Precision = 1.0e-3;
    int nbIter = 0;
    int nbIterMax = 100;

    double EPS_DERIV = 1.0e-3;
    double H;

    double PHI; 
    double DerivPHI;

    ARM_Date matuDate;

    matuDate = AsOf;
    matuDate = matuDate.AddMonths(int(matu*12.0));

    // Calculate the Forward FX

    double F = CalcFwdFXSpot(AsOf,
                             fxSpot,
                             matuDate,
                             domCrv,
                             foreignCrv);

    // Calculate discount factor

    char payCal[30];

    // Get Currency names

    char* domCcy = domCrv->GetCurrencyUnit()->GetCcyName();
    char* forCcy = foreignCrv->GetCurrencyUnit()->GetCcyName();
   
    ARM_Currency MainCurrency(forCcy);
    ARM_Currency MoneyCurrency(domCcy);

    strcpy(payCal, MainCurrency.GetCcyName());

    strcat(payCal, MoneyCurrency.GetCcyName());

    ARM_Date fwdDate(matuDate);

    fwdDate.NextBusinessDay(MoneyCurrency.GetFxSpotDays(), payCal);
                
    double matDelivery = CountYears(KACTUAL_365, AsOf, fwdDate);

    double dfDom = domCrv->DiscountPrice(matDelivery);
    double dfFor = foreignCrv->DiscountPrice(matDelivery);

    double premium;
    
    // Init.

    double strike = fxSpot; // TMP F;
    double prevStrike;

    deltaWithPremium = deltaWithPremium/100.0;

    if ( deltaWithPremium > 0.0 )
    {
       callPut = K_CALL;
    }
    else
    {
       callPut = K_PUT;
    }

    // Solve
    do
    {
        prevStrike = strike;

        // Compute Premium

        if ( callPut == K_CALL )
        {
           premium = FXCallFwd(F, strike, sigma/100.0, matu);
        }
        else
        {
           premium = FXPutFwd(F, strike, sigma/100.0, matu);
        }

        premium *= dfDom;
       
        // TMP premium = premium/100.0;

        // End Premium calculation

        // Compute DELTA Forward

        double deltaFwd = ComputeFwdBSDelta(F,
                                            strike, 
                                            sigma/100.0,
                                            matu,
                                            callPut);


        
        // Calculate Derivative depending on the strike

        PHI = deltaWithPremium-deltaFwd+(premium/fxSpot/dfFor);

        double derivDeltaFwd;

        H = EPS_DERIV*strike;

        derivDeltaFwd = ComputeFwdBSDelta(F,
                                          strike+H, 
                                          sigma/100.0,
                                          matu,
                                          callPut);
                        
        derivDeltaFwd = (derivDeltaFwd-deltaFwd)/H;

        // Calculate Premium derivative
        double derivPremium;

        if ( callPut == K_CALL )
        {
           derivPremium = FXCallFwd(F, strike+H, sigma/100.0, matu);
        }
        else
        {
           derivPremium = FXPutFwd(F, strike+H, sigma/100.0, matu);
        }

        derivPremium *= dfDom;
        // TMP derivPremium = derivPremium/100.0;

        derivPremium = (derivPremium-premium)/H;

        DerivPHI = -derivDeltaFwd+(derivPremium/fxSpot/dfFor);

        strike = prevStrike-(PHI/DerivPHI);

        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax )
          );

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "CalculateImpliedStrikeFromDeltaWithPremium, did not converge");

       return(0.0);
    }

    return(strike);
}




double ComputeDeltaFwdFromDeltaWP(ARM_Date& AsOf,
                                  double matu,
                                  double sigma,
                                  double fxSpot,
                                  double deltaWithPremium,
                                  ARM_ZeroCurve* domCrv, // JPY
                                  ARM_ZeroCurve* foreignCrv) // USD
{
    double deltaFwd;
    int    callPut;


    if ( deltaWithPremium > 0.0 )
    {
       callPut = K_CALL;
    }
    else
    {
       callPut = K_PUT;
    }

    ARM_Date matuDate;

    matuDate = matuDate.AddMonths(int(matu*12.0));


    try
    {
        double strike;

        strike = CalculateImpliedStrikeFromDeltaWithPremium(AsOf,
                                                            matu,
                                                            sigma,
                                                            fxSpot,
                                                            deltaWithPremium,
                                                            domCrv,
                                                            foreignCrv);

        // Calculate the FX Forward

        double F = CalcFwdFXSpot(AsOf,
                                 fxSpot,
                                 matuDate,
                                 domCrv,
                                 foreignCrv);

        // Calculate discount factor

        char payCal[30];

        // Get Currency names

        char* domCcy = domCrv->GetCurrencyUnit()->GetCcyName();
        char* forCcy = foreignCrv->GetCurrencyUnit()->GetCcyName();
   
        ARM_Currency MainCurrency(forCcy);
        ARM_Currency MoneyCurrency(domCcy);

        strcpy(payCal, MainCurrency.GetCcyName());

        strcat(payCal, MoneyCurrency.GetCcyName());

        ARM_Date fwdDate(matuDate);

        fwdDate.NextBusinessDay(MoneyCurrency.GetFxSpotDays(), payCal);
                
        double matDelivery = CountYears(KACTUAL_365, AsOf, fwdDate);
        double dfDom = domCrv->DiscountPrice(matDelivery);
    

        double premium;
    
        // Compute Premium

        if ( callPut == K_CALL )
        {
           premium = FXCallFwd(F, strike, sigma/100.0, matu);
        }
        else
        {
           premium = FXPutFwd(F, strike, sigma/100.0, matu);
        }

        premium *= dfDom;

        // End Premium calculation

        double dfFor = foreignCrv->DiscountPrice(matDelivery);

        deltaFwd = deltaWithPremium+(premium*100.0/fxSpot/dfFor);
    }

    catch(Exception& anExpt)
    {
        throw anExpt;
    }

    catch(...)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeDeltaFwdFromDeltaWP: Unrecognized Failure");

       return(0.0);
    }

    return(deltaFwd);
}




// whatever delta, not necessarily fwd with premium
double ComputeSigmaFromDelta(double curDelta,
                             ARM_Vector* deltaCall, ARM_Vector* volsCall,
                             ARM_Vector* deltaPut, ARM_Vector* volsPut,
                             int interpolType)
{
    double sigmaRes;
    

    if ( curDelta < 0.0 )
    {
       // callPut = K_PUT;

       if ( interpolType == K_LINEAR )
       {
          sigmaRes = linInterpol2(deltaPut, curDelta, volsPut);
       }
       else
       {
          sigmaRes = SplineInterpolateFunc(deltaPut, volsPut,
                                           curDelta,
                                           NULL, // SecondDerivCalc,
                                           1);   // keep2Der 
       }
    }
    else
    {
       // callPut = K_CALL;

       if ( interpolType == K_LINEAR )
       {
          sigmaRes = linInterpol2(deltaCall, curDelta, volsCall);
       }
       else
       {
          sigmaRes = SplineInterpolateFunc(deltaCall, volsCall,
                                           curDelta,
                                           NULL, // SecondDerivCalc,
                                           1);   // keep2Der 
       }
    }
  
    return(sigmaRes);
}                             



double ComputeSigmaFromDeltaSIMPLEX(double curDelta,
                                    ARM_Vector* delta, ARM_Vector* vols,
                                    int interpolType)
{
    double sigmaRes;


    if ( interpolType == K_LINEAR )
    {
       sigmaRes = linInterpol2(delta, curDelta, vols);
    }
    else
    {
       sigmaRes = SplineInterpolateFunc(delta, vols,
                                        curDelta);
    }

    return(sigmaRes);
}                             



/*-----------------------------------------------------------------------------------*/
/*                               New functions For FX VOL                            */
/*-----------------------------------------------------------------------------------*/

// Utility function computing Domestic and Foreign Discount Factor

void ComputeFxDomForDiscountFactorsAtSettlmentDate(ARM_ZeroCurve* zcDom, 
                                                   ARM_ZeroCurve* zcFor,
                                                   double matu,
                                                   // Output
                                                   double& dfDom,
                                                   double& dfFor)
{
    ARM_Date matuDate;


    ARM_Date AsOf = zcDom->GetAsOfDate();

    matuDate      = AsOf;

    matuDate      = matuDate.AddDays(int(matu*365.0));

    // Calculate discount factors

    char payCal[30];

    // Get Currency names

    char* domCcy = zcDom->GetCurrencyUnit()->GetCcyName();
    char* forCcy = zcFor->GetCurrencyUnit()->GetCcyName();

    ARM_Currency MainCurrency(forCcy);
    ARM_Currency MoneyCurrency(domCcy);

    strcpy(payCal, MainCurrency.GetCcyName());

    strcat(payCal, MoneyCurrency.GetCcyName());

    ARM_Date settleDate(matuDate);

    // Get FX spot days

    int fxSpotDom = MoneyCurrency.GetFxSpotDays();
    int fxSpotFor = MainCurrency.GetFxSpotDays();

    int fxSpotDays = MAX(fxSpotDom, fxSpotFor);


    settleDate.NextBusinessDay(fxSpotDays, payCal);
                
    double matDelivery = CountYears(KACTUAL_365, AsOf, settleDate);

    dfDom = zcDom->DiscountPrice(matDelivery);
    dfFor = zcFor->DiscountPrice(matDelivery);
}



// Utility function computing an FX Spot Delta 
// Delta_spot= dPfwd/dS = (dPfwd/dF)*(dF/dS) = (dPfwd/dF)*(BFor/Bdom)

double ComputeSpotFxDeltaWithoutPremium(double volZDS, 
                                        double fxFwd, 
                                        double spot, 
                                        double strike,
                                        double matu,
                                        double dfDom,
                                        double dfFor,
                                        int callPut)
{
    double fwdDelta;

    fwdDelta = ComputeFwdBSDelta(fxFwd, strike, volZDS, matu, callPut);

// TMP: Correction JCH?!    double deltaSpot = fwdDelta*(dfFor/dfDom);

    double deltaSpot = fwdDelta*dfFor;

    return(deltaSpot);
}


/*
double ComputeSpotFxDeltaWithoutPremium(double volZDS, 
                                        double fxFwd, 
                                        double spot, 
                                        double strike,
                                        double matu,
                                        ARM_ZeroCurve* zcDom,
                                        ARM_ZeroCurve* zcFor,
                                        int callPut)
{
    double dfDom = 1.0;
    double dfFor = 1.0;


    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);

    double deltaSpot = ComputeSpotFxDeltaWithoutPremium(volZDS, 
                                                        fxFwd, 
                                                        spot, 
                                                        strike,
                                                        matu,
                                                        dfDom,
                                                        dfFor,
                                                        callPut);

    return(deltaSpot);
}
*/


double ComputeSpotFxDeltaWithPremium(double volZDS, 
                                     double fxFwd, 
                                     double spot, 
                                     double strike,
                                     double matu,
                                     double dfDom,
                                     double dfFor,
                                     int callPut)
{
    double deltaWithoutPrem = ComputeSpotFxDeltaWithoutPremium(volZDS, 
                                                               fxFwd, 
                                                               spot, 
                                                               strike,
                                                               matu,
                                                               dfDom,
                                                               dfFor,
                                                               callPut);

    double premium;

    if ( callPut == K_CALL )
    {
       premium = FXCallFwd(fxFwd, strike, volZDS, matu);
    }
    else
    {
       premium = FXPutFwd(fxFwd, strike, volZDS, matu);
    }

    premium *= dfDom;

    double deltaWithPrem = deltaWithoutPrem-(premium/spot);

    return(deltaWithPrem);
}



double DerivSpotFxDeltaWithPremium(double volZDS, 
                                   double fxFwd, 
                                   double spot, 
                                   double strike,
                                   double matu,
                                   double dfDom,
                                   double dfFor,
                                   int callPut)
{
    double Precision = ARM_FX_GEN_VOL_PREC; // 1.0e-3;
  
    double EPS_DERIV = 1.0e-3;
    double H;

    H = EPS_DERIV*strike;

    double d_Delta_PLUS_H, d_Delta_MOINS_H;

    d_Delta_PLUS_H = ComputeSpotFxDeltaWithPremium(volZDS, 
                                                   fxFwd, 
                                                   spot, 
                                                   strike+H,
                                                   matu,
                                                   dfDom,
                                                   dfFor,
                                                   callPut);

    d_Delta_MOINS_H = ComputeSpotFxDeltaWithPremium(volZDS, 
                                                    fxFwd, 
                                                    spot, 
                                                    strike-H,
                                                    matu,
                                                    dfDom,
                                                    dfFor,
                                                    callPut);

    double deriv = (d_Delta_PLUS_H-d_Delta_MOINS_H)/(2.0*H);

    return(deriv);
}



void ARM_ComputeImpliedStrikeFromSpotDeltaWithPremium(double vol, 
                                                      double fxFwd,
                                                      double spot,
                                                      double matu, 
                                                      double target, 
                                                      int callPut,
                                                      ARM_ZeroCurve* zcDom,
                                                      ARM_ZeroCurve* zcFor,
                                                      double& strike)

{
    // look for implicit strike such that delta spot = target for the call or put
    double Precision = ARM_FX_GEN_VOL_PREC; // 1.0e-3;
    int nbIter = 0;
    int nbIterMax = 100;
    
    double EPS_DERIV = 1.0e-3;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
    double DeltaSpot = 0 ;
    double d_DeltaSpot;
    double dfDom, dfFor;


    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);

    strike = fxFwd;
 
    do
    {
        prevStrike = strike;
    
        // Compute DELTA Spot


        DeltaSpot = ComputeSpotFxDeltaWithPremium(vol/100.0, 
                                                  fxFwd, 
                                                  spot, 
                                                  strike,
                                                  matu,
                                                  dfDom,
                                                  dfFor,
                                                  callPut);
        
        // Calculate Derivative depending on the strike
        PHI =  target-DeltaSpot;
        H = EPS_DERIV*strike;

        d_DeltaSpot = DerivSpotFxDeltaWithPremium(vol/100.0, 
                                                  fxFwd, 
                                                  spot, 
                                                  strike,
                                                  matu,
                                                  dfDom,
                                                  dfFor,
                                                  callPut);

   
        
        // Calculate Premium derivative
        
        DerivPHI = -d_DeltaSpot;
        
        strike = prevStrike-(PHI/DerivPHI);
        
        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax ) );


    if ( nbIter == nbIterMax )
    {
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedStrike, did not converge");
    }
}



void ARM_ComputeImpliedStrikeFromSpotDeltaWithoutPremium(double vol, 
                                                         double fxFwd,
                                                         double spot,
                                                         double matu, 
                                                         double target, 
                                                         int callPut,
                                                         ARM_ZeroCurve* zcDom,
                                                         ARM_ZeroCurve* zcFor,
                                                         double& strike)

{
    // look for implicit strike such that delta spot = target for the call or put
    double Precision = ARM_FX_GEN_VOL_PREC;
    int nbIter = 0;
    int nbIterMax = 100;
    
    double EPS_DERIV = ARM_FX_GEN_VOL_PREC;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
    double DeltaSpot = 0 ;
    double d_DeltaSpot;
    double dfDom, dfFor;


    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);

    strike = spot;
 
    do
    {
        prevStrike = strike;
    
        // Compute DELTA Spot

        DeltaSpot = ComputeSpotFxDeltaWithoutPremium(vol/100.0, 
                                                     fxFwd, 
                                                     spot, 
                                                     strike,
                                                     matu,
                                                     dfDom,
                                                     dfFor,
                                                     callPut);
        
        // Calculate Derivative depending on the strike
        PHI =  target-DeltaSpot;
        H = EPS_DERIV*strike;

        d_DeltaSpot = ComputeSpotFxDeltaWithoutPremium(vol/100.0, 
                                                     fxFwd, 
                                                     spot, 
                                                     strike+H,
                                                     matu,
                                                     dfDom,
                                                     dfFor,
                                                     callPut);

        d_DeltaSpot -= ComputeSpotFxDeltaWithoutPremium(vol/100.0, 
                                                        fxFwd, 
                                                        spot, 
                                                        strike-H,
                                                        matu,
                                                        dfDom,
                                                        dfFor,
                                                        callPut);
        d_DeltaSpot /= 2*H;
        
        // Calculate Premium derivative
        
        DerivPHI = -d_DeltaSpot;
        
        strike = prevStrike-(PHI/DerivPHI);
        
        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax ) );


    if ( nbIter == nbIterMax )
    {
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedStrike, did not converge");
    }
}



// Correspond To: PIVOT_IS_ZDS_SPOT_NOPREM 5
void ARM_ComputeImpliedZDSWithoutPremiumSpot(double volZDS, double fxFwd, double spot, 
                                             double matu,
                                             ARM_ZeroCurve* zcDom,
                                             ARM_ZeroCurve* zcFor,
                                             double& strike, 
                                             double& DeltaCall, double& DeltaPut)
{
    // look for implicit strike such that delta fwd with premium Call+Put=0

    double Precision = ARM_FX_GEN_VOL_PREC; // 1.0e-3;
    int nbIter = 0;
    int nbIterMax = 100;
  
    double EPS_DERIV = 1.0e-3;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
	
    DeltaCall = 0.0;
    DeltaPut  = 0.0;

    double d_DeltaCall;
    double d_DeltaPut;
	

    strike    = fxFwd;	
    	

    double dfDom;
    double dfFor;


    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);
    
    do
    {
        prevStrike = strike;


        // Compute DELTA Spot 

        DeltaCall = ComputeSpotFxDeltaWithoutPremium(volZDS/100.0, 
                                                     fxFwd, 
                                                     spot, 
                                                     strike,
                                                     matu,
                                                     dfDom,
                                                     dfFor,
                                                     K_CALL);

        DeltaPut = ComputeSpotFxDeltaWithoutPremium(volZDS/100.0, 
                                                    fxFwd, 
                                                    spot, 
                                                    strike,
                                                    matu,
                                                    dfDom,
                                                    dfFor,
                                                    K_PUT);
        

        // Calculate Derivative depending on the strike

        PHI = (DeltaCall+DeltaPut);


        H = EPS_DERIV*strike;

        d_DeltaCall = ComputeSpotFxDeltaWithoutPremium(volZDS/100.0, 
                                                       fxFwd, 
                                                       spot, 
                                                       strike+H,
                                                       matu,
                                                       dfDom,
                                                       dfFor,
                                                       K_CALL);

        d_DeltaPut  = ComputeSpotFxDeltaWithoutPremium(volZDS/100.0, 
                                                       fxFwd, 
                                                       spot, 
                                                       strike+H,
                                                       matu,
                                                       dfDom,
                                                       dfFor,
                                                       K_PUT);
       
        d_DeltaCall -= ComputeSpotFxDeltaWithoutPremium(volZDS/100.0, 
                                                        fxFwd, 
                                                        spot, 
                                                        strike-H,
                                                        matu,
                                                        dfDom,
                                                        dfFor,
                                                        K_CALL);

        d_DeltaPut  -= ComputeSpotFxDeltaWithoutPremium(volZDS/100.0, 
                                                        fxFwd, 
                                                        spot, 
                                                        strike-H,
                                                        matu,
                                                        dfDom,
                                                        dfFor,
                                                        K_PUT);

	
        d_DeltaCall /= 2.0*H;
	    d_DeltaPut  /= 2.0*H;


        DerivPHI = (d_DeltaCall+d_DeltaPut);

        strike = prevStrike-(PHI/DerivPHI);

        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax )
          );

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedZDSWithoutPremiumSpot, did not converge");

    }
}



// Correspond To: PIVOT_IS_ZDS_SPOT_WPREM  4
void ARM_ComputeImpliedZDSWithPremiumSpot(double volZDS, double fxFwd, double spot, 
                                          double matu,
                                          ARM_ZeroCurve* zcDom,
                                          ARM_ZeroCurve* zcFor,
                                          double& strike, 
                                          double& DeltaCall, double& DeltaPut)
{
    // look for implicit strike such that delta fwd with premium Call+Put=0

    double Precision = ARM_FX_GEN_VOL_PREC; // 1.0e-3;
    int nbIter = 0;
    int nbIterMax = 100;
  
    double EPS_DERIV = 1.0e-3;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
	
    DeltaCall = 0.0;
    DeltaPut  = 0.0;

    double d_DeltaCall;
    double d_DeltaPut;
	

    strike    = fxFwd;
	
    	

    double dfDom;
    double dfFor;


    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);
    
    do
    {
        prevStrike = strike;


        // Compute DELTA Spot 

        DeltaCall = ComputeSpotFxDeltaWithPremium(volZDS/100.0, 
                                                  fxFwd, 
                                                  spot, 
                                                  strike,
                                                  matu,
                                                  dfDom,
                                                  dfFor,
                                                  K_CALL);

        DeltaPut = ComputeSpotFxDeltaWithPremium(volZDS/100.0, 
                                                 fxFwd, 
                                                 spot, 
                                                 strike,
                                                 matu,
                                                 dfDom,
                                                 dfFor,
                                                 K_PUT);
        

        // Calculate Derivative depending on the strike

        PHI = (DeltaCall+DeltaPut);


        H = EPS_DERIV*strike;

        d_DeltaCall = DerivSpotFxDeltaWithPremium(volZDS/100.0, 
                                                  fxFwd, 
                                                  spot, 
                                                  strike,
                                                  matu,
                                                  dfDom,
                                                  dfFor,
                                                  K_CALL);

        d_DeltaPut  = DerivSpotFxDeltaWithPremium(volZDS/100.0, 
                                                  fxFwd, 
                                                  spot, 
                                                  strike,
                                                  matu,
                                                  dfDom,
                                                  dfFor,
                                                  K_PUT);

        DerivPHI = (d_DeltaCall+d_DeltaPut);

        strike = prevStrike-(PHI/DerivPHI);

        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax )
          );

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedZDSWithPremiumSpot, did not converge");

    }
}



// Correspond To: PIVOT_IS_ZDS_NOPREM      3
void ARM_ComputeImpliedZDSWithoutPremium(double volZDS, double fxFwd, double matu,
                                         double& strike, 
                                         double& deltaCall, double& deltaPut)
{
    // look for implicit strike such that delta fwd with premium Call+Put=0

    double Precision = ARM_FX_GEN_VOL_PREC;
    int nbIter = 0;
    int nbIterMax = 100;
  
    double EPS_DERIV = 1.0e-3;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
	
  
    double fwdDeltaCall = 0.0;
    double fwdDeltaPut  = 0.0;

    double d_fwdDeltaCall;
    double d_fwdDeltaPut;
	
    strike    = fxFwd;
	
    	
    do
    {
        prevStrike = strike;


        // Compute DELTA Forward

        fwdDeltaCall = ComputeFwdBSDelta(fxFwd, strike, volZDS/100.0, matu, K_CALL);
        fwdDeltaPut  = ComputeFwdBSDelta(fxFwd, strike, volZDS/100.0, matu, K_PUT);
        

        // Calculate Derivative depending on the strike

        PHI = (fwdDeltaCall+fwdDeltaPut);


        H = EPS_DERIV*strike;

        d_fwdDeltaCall = ComputeFwdBSDelta(fxFwd, strike+H, volZDS/100.0, matu, K_CALL);
        d_fwdDeltaPut  = ComputeFwdBSDelta(fxFwd, strike+H, volZDS/100.0, matu, K_PUT);
       
        d_fwdDeltaCall -= ComputeFwdBSDelta(fxFwd, strike-H, volZDS/100.0, matu, K_CALL);
        d_fwdDeltaPut  -= ComputeFwdBSDelta(fxFwd, strike-H, volZDS/100.0, matu, K_PUT);
	
        d_fwdDeltaCall /= 2.0*H;
	    d_fwdDeltaPut  /= 2.0*H;


        DerivPHI = (d_fwdDeltaCall+d_fwdDeltaPut);

        strike = prevStrike-(PHI/DerivPHI);

        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax )
          );

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedZDSWithoutPremium, did not converge");

    }

    deltaCall = fwdDeltaCall;
    deltaPut  = fwdDeltaPut;
}



// Correspond To: PIVOT_IS_ATMF_ANDNOTWP   2
void ARM_ComputeImpliedATMF(double volATMF, double fxFwd, double matu,
                            double& strike, double& deltaCall, double& deltaPut)
{
    // look for implicit strike such that delta fwd with premium Call+Put=0

	double fwdDeltaCall = 0 ;
	double fwdDeltaPut = 0 ;
	strike = fxFwd ;
        
    // Compute DELTA Forward

    fwdDeltaCall = ComputeFwdBSDelta(fxFwd, strike, volATMF/100.0, matu, K_CALL);
    fwdDeltaPut  = ComputeFwdBSDelta(fxFwd, strike, volATMF/100.0, matu, K_PUT);
    
    // Calculate Derivative depending on the strike

    deltaCall = fwdDeltaCall;
	deltaPut  = fwdDeltaPut;
}



// Correspond To: PIVOT_IS_ATMF  1
void ARM_ComputeImpliedATMF_FWP(double volATMF, double fxFwd, double matu,
                                double& strike, double& deltaCall, double& deltaPut)
{
    // look for implicit strike such that delta fwd with premium Call+Put=0

    double fwdPremiumCall = 0.0;	
	double fwdDeltaCall = 0.0;
	double fwdPremiumPut = 0.0;	
	double fwdDeltaPut = 0.0;
	
    strike = fxFwd ;
        

    // Compute Premium

    fwdPremiumCall = FXCallFwd(fxFwd, strike, volATMF/100.0, matu);
    fwdPremiumPut  = FXPutFwd(fxFwd, strike, volATMF/100.0, matu);

    // End Premium calculation

    // Compute DELTA Forward

    fwdDeltaCall = ComputeFwdBSDelta(fxFwd, strike, volATMF/100.0, matu, K_CALL);
    fwdDeltaPut  = ComputeFwdBSDelta(fxFwd, strike, volATMF/100.0, matu, K_PUT);
    
    // Calculate Derivative depending on the strike

    deltaCall = fwdDeltaCall-fwdPremiumCall/fxFwd ;
	deltaPut  = fwdDeltaPut-fwdPremiumPut/fxFwd ;
}



// Correspond To: PIVOT_IS_ZDS    0
void ARM_ComputeImpliedZDS_FWP(double volZDS, double fxFwd, double matu,
                               double& strike, double& deltaCall, double& deltaPut)
{
    // look for implicit strike such that delta fwd with premium Call+Put=0

    double Precision = 1.0e-3;
    int nbIter = 0;
    int nbIterMax = 100;
    double EPS_DERIV = 1.0e-3;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;

	double fwdPremiumCall = 0.0;	
	double fwdDeltaCall   = 0.0;
	double fwdPremiumPut  = 0.0;	
	double fwdDeltaPut    = 0.0;
    double d_fwdPremiumCall;
    double d_fwdPremiumPut;
    double d_fwdDeltaCall;
    double d_fwdDeltaPut;
	strike = fxFwd ;
	
    deltaCall = 0.5;	
    	
    do
    {
        prevStrike = strike;

        // Compute Premium

        fwdPremiumCall = FXCallFwd(fxFwd, strike, volZDS/100.0, matu);
        fwdPremiumPut  = FXPutFwd(fxFwd, strike, volZDS/100.0, matu);

        // End Premium calculation


        // Compute DELTA Forward

        fwdDeltaCall = ComputeFwdBSDelta(fxFwd, strike, volZDS/100.0, matu, K_CALL);
        fwdDeltaPut  = ComputeFwdBSDelta(fxFwd, strike, volZDS/100.0, matu, K_PUT);
        

        // Calculate Derivative depending on the strike

        PHI =  0 - (
            (fwdDeltaCall+fwdDeltaPut) - (fwdPremiumCall+fwdPremiumPut)/fxFwd );

        H = EPS_DERIV*strike;

        d_fwdDeltaCall = ComputeFwdBSDelta(fxFwd, strike+H, volZDS/100.0, matu, K_CALL);
        d_fwdDeltaPut  = ComputeFwdBSDelta(fxFwd, strike+H, volZDS/100.0, matu, K_PUT);
       
        d_fwdDeltaCall -= ComputeFwdBSDelta(fxFwd, strike-H, volZDS/100.0, matu, K_CALL);
        d_fwdDeltaPut  -= ComputeFwdBSDelta(fxFwd, strike-H, volZDS/100.0, matu, K_PUT);
	
        d_fwdDeltaCall /= 2.0*H;
	    d_fwdDeltaPut  /= 2.0*H;

        // Calculate Premium derivative

        d_fwdPremiumCall = FXCallFwd(fxFwd, strike+H, volZDS/100.0, matu);
        d_fwdPremiumPut = FXPutFwd(fxFwd, strike+H, volZDS/100.0, matu);
        
        d_fwdPremiumCall -= FXCallFwd(fxFwd, strike-H, volZDS/100.0, matu);
        d_fwdPremiumPut -= FXPutFwd(fxFwd, strike-H, volZDS/100.0, matu);
        
        d_fwdPremiumCall /= 2*H;
        d_fwdPremiumPut  /= 2*H;


        DerivPHI = -(
		(d_fwdDeltaCall+d_fwdDeltaPut)-(d_fwdPremiumCall+d_fwdPremiumPut)/fxFwd);

        strike = prevStrike-(PHI/DerivPHI);

        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax )
          );

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedZDS_FWP, did not converge");

    }

    deltaCall = fwdDeltaCall-fwdPremiumCall/fxFwd ;
    deltaPut  = fwdDeltaPut-fwdPremiumPut/fxFwd ;
}




void ARM_ComputeImpliedStrike_FWP(double vol, double fxFwd, double matu, 
                                  double target, int callPut,
                                  double& strike, double firstStrike, // solution précédente
							      bool isFirstTime)

{

    // look for implicit strike such that delta fwd with :
    // premium = target for the call or put
   
    double Precision = 1.0e-2;
    int nbIter = 0;
    int nbIterMax = 1000;
    double EPS_DERIV = 1e-7;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
    double fwdPremium = 0 ;           
    double fwdDelta = 0 ;
    double d_fwdPremium;
    double d_fwdDelta;
    double move;

    strike = fxFwd;


    do
    {
        prevStrike = strike;
        
        // Compute Premium
        if ( callPut == K_CALL )
           fwdPremium = FXCallFwd(fxFwd, strike, vol/100.0, matu);
        else
           fwdPremium = FXPutFwd(fxFwd, strike, vol/100.0, matu);
        // End Premium calculation

        // Compute DELTA Forward
        fwdDelta = ComputeFwdBSDelta(fxFwd, strike, vol/100.0, matu, callPut);

        // Calculate Derivative depending on the strike
        PHI = target-(fwdDelta-fwdPremium/fxFwd);
        H = EPS_DERIV*strike;
        d_fwdDelta = ComputeFwdBSDelta(fxFwd, strike+H, vol/100.0, matu, callPut);
        d_fwdDelta -= ComputeFwdBSDelta(fxFwd, strike-H, vol/100.0, matu, callPut);
        d_fwdDelta /= 2*H;

        // Calculate Premium derivative
        if ( callPut == K_CALL )
        {
           d_fwdPremium = FXCallFwd(fxFwd, strike+H, vol/100.0, matu);
           d_fwdPremium -= FXCallFwd(fxFwd, strike-H, vol/100.0, matu);
        }
        else
        {
           d_fwdPremium = FXPutFwd(fxFwd, strike+H, vol/100.0, matu);
           d_fwdPremium -= FXPutFwd(fxFwd, strike-H, vol/100.0, matu);
        }

        d_fwdPremium /= 2*H;

        DerivPHI = -( d_fwdDelta - d_fwdPremium/fxFwd);

        move = PHI/DerivPHI;

		if (!(isFirstTime))
		{
			if ( fabs(DerivPHI) < 0.0001)
				move = -move;
		}

		if ( move < -0.1*fxFwd )
           move = -0.1*fxFwd;
        else
        {
           if ( move > 0.1*fxFwd )
              move = 0.1*fxFwd;
        }

        strike = prevStrike-move;

        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision ) && ( nbIter < nbIterMax ) );

	double secondStrike;
    
    if ( nbIter == nbIterMax )
    {
		if (!(isFirstTime))
		{
			if (firstStrike < 0.0)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"ComputeImpliedStrike_FWP, did not converge");
			}
			else
			{
				strike = firstStrike;
				return;
			}
		}

		isFirstTime = 0;
		
		ARM_ComputeImpliedStrike_FWP(vol, fxFwd, matu,
								     target, callPut,
								      secondStrike, firstStrike, isFirstTime);

    }

	if (isFirstTime)
	{
        isFirstTime = 0;

		double ratio = fxFwd / strike;

		if (ratio >= 8.0)
		{
			firstStrike = strike;

			ARM_ComputeImpliedStrike_FWP(vol, fxFwd, matu,
									     target, callPut,
									     secondStrike, firstStrike, isFirstTime);

			double secondRatio = secondStrike / fxFwd;

			if (ratio > secondRatio)
				strike = secondStrike;
		}
	}
}

 

void ARM_ComputeImpliedStrike(double vol, double fxFwd, double matu, double target, 
                              int callPut,
                              double& strike)

{
    // look for implicit strike such that delta fwd = target for the call or put
    double Precision = 1.0e-3;
    int nbIter = 0;
    int nbIterMax = 100;
    double EPS_DERIV = 1.0e-3;
    double H;
    double PHI; 
    double DerivPHI;
    double prevStrike;
    double fwdDelta = 0 ;
    double d_fwdDelta;
    strike = fxFwd ;
 
    do
    {
        prevStrike = strike;
    
        // Compute DELTA Forward
        fwdDelta = ComputeFwdBSDelta(fxFwd, strike, vol/100.0, matu, callPut);
        
        // Calculate Derivative depending on the strike
        PHI =  target-fwdDelta;
        H = EPS_DERIV*strike;
        d_fwdDelta = ComputeFwdBSDelta(fxFwd, strike+H, vol/100.0, matu, callPut);
        d_fwdDelta -= ComputeFwdBSDelta(fxFwd, strike-H, vol/100.0, matu, callPut);
        d_fwdDelta /= 2*H;
        
        // Calculate Premium derivative
        DerivPHI = - d_fwdDelta;
        strike = prevStrike-(PHI/DerivPHI);
        nbIter++;
    }
    while (( fabs(strike-prevStrike) > Precision )
           &&
           ( nbIter < nbIterMax ) );


    if ( nbIter == nbIterMax )
    {
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeImpliedStrike, did not converge");
    }
}




                         /****************************/


double computeDeltaFWP_withStrikeAndVol(double strike, double vol, 
                                        double fxFwd, double matu, int callPut)
{
    double fwdPremium = 0;
    double fwdDelta = 0;

    
    if ( callPut == K_CALL ) 
    {
       fwdPremium = FXCallFwd(fxFwd, strike, vol/100.0, matu);
    } 
    else 
    {
       fwdPremium = FXPutFwd(fxFwd, strike, vol/100.0, matu);
    }

    fwdDelta = ComputeFwdBSDelta(fxFwd, strike, vol/100.0, matu, callPut);
    
    return(fwdDelta-fwdPremium/fxFwd);
}



double computeVol_withDeltaFWPasInput_aux(double strike, double fxFwd, double matu,
                                          ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                          ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                          double pivotStrike, double guessVol, 
                                          int interpolType, int correctSplineWithLinear)
{
    double vol = guessVol;
    double precision = 0.0001;
    double volPrec = guessVol;
    int nbIter = 0; 
    int nbIterMax = 2000;
    double delta = 0.0;
    int callPut;
    
    
    do 
    {
        if ( strike < pivotStrike )
        {
           callPut = K_PUT;
        }
        else
        {
            callPut = K_CALL;
        }

        vol = volPrec;
        delta = computeDeltaFWP_withStrikeAndVol(strike, vol, fxFwd, matu, callPut);
        volPrec = ComputeSigmaFromDelta(delta, deltaCall, volsCall,
                                        deltaPut, volsPut, interpolType);
		nbIter++;

    } while(fabs(vol - volPrec)>precision && nbIter<nbIterMax);

    if ( nbIter == nbIterMax )
    {
        if ( interpolType == K_LINEAR || !correctSplineWithLinear )
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "computeVol_withDeltaFWPasInput_aux, did not converge");
        }
        else  //second chance, usually linear mode converges better
        {
            try
            {
                vol = computeVol_withDeltaFWPasInput_aux(strike, fxFwd, matu,
                                                         deltaCall, volsCall,
                                                         deltaPut, volsPut,
                                                         pivotStrike, guessVol,
                                              K_LINEAR, correctSplineWithLinear);
                return(vol);
            }

            catch(Exception& )
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "computeVol_withDeltaFWPasInput_aux, did not converge (even in linear mode)");
            }
        }
    }

    return(vol);
}



double computeDelta_withStrikeAndVol(double strike, double vol, double fxFwd, 
                                     double matu, int callPut)
{
    double fwdDelta = 0;
 
    fwdDelta = ComputeFwdBSDelta(fxFwd, strike, vol/100.0, matu, callPut);
    
    return fwdDelta;
}



double computeVol_withDeltaAsInput_aux(double strike, double fxFwd, double matu,
                                       ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                       ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                       double pivotStrike,
                                       int interpolType)
{
    double vol = 20;
    double precision = 0.001;
    double volPrec = 20;
    int nbIter = 0; 
    int nbIterMax = 100;
    double delta = 0;
    int callPut;


    do 
    {
        if ( strike < pivotStrike )
        {
           callPut = K_PUT;
        }
        else
        {
           callPut = K_CALL;
        }

        vol = volPrec;
        delta = computeDelta_withStrikeAndVol(strike, vol, fxFwd, matu, callPut);
        
        volPrec = ComputeSigmaFromDelta(delta, deltaCall, volsCall, deltaPut,
                                        volsPut, interpolType);
		
        nbIter++;

    } while(fabs(vol - volPrec)>precision && nbIter<nbIterMax);

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "computeVol_withDeltaAsInput_aux, did not converge");

    }

    return(vol);
}



double ComputeVolWithDeltaSpotAsInputWithPremium(double strike, double fxFwd, 
                                                 double FXSpot,
                                                 double matu,
                                                 ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                                 ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                                 double pivotStrike, double guessVol, 
                                                 int interpolType, 
                                                 int correctSplineWithLinear,
                                                 ARM_ZeroCurve* zcDom,
                                                 ARM_ZeroCurve* zcFor)
{
    double vol = guessVol;
    double precision = ARM_FX_GEN_VOL_PREC;
    double volPrec = guessVol;
    int nbIter = 0; 
    int nbIterMax = 2000;
    double delta = 0.0;
    int callPut;
    
    double dfDom;
    double dfFor;


    // Compute discount factors

    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);    
    do 
    {
        if ( strike < pivotStrike )
        {
           callPut = K_PUT;
        }
        else
        {
           callPut = K_CALL;
        }

        vol = volPrec;

        delta = ComputeSpotFxDeltaWithoutPremium(vol/100.0,
                                                 fxFwd, 
                                                 FXSpot, 
                                                 strike,
                                                 matu,
                                                 dfDom,
                                                 dfFor,
                                                 callPut);
        
        volPrec = ComputeSigmaFromDelta(delta, deltaCall, volsCall,
                                        deltaPut, volsPut, interpolType);
		
        nbIter++;

    } 
    while((fabs(vol-volPrec) > precision ) && ( nbIter < nbIterMax ));

    if ( nbIter == nbIterMax )
    {
        if ( interpolType == K_LINEAR || !correctSplineWithLinear )
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeVolWithDeltaSpotAsInputWithPremium, did not converge");
        }
        else  // second chance, usually linear mode converges better
        {
            try
            {
                vol = ComputeVolWithDeltaSpotAsInputWithPremium(strike, fxFwd, 
                                                                FXSpot,
                                                                matu,
                                                                deltaCall, volsCall,
                                                                deltaPut, volsPut,
                                                                pivotStrike, guessVol,
                                                      K_LINEAR, correctSplineWithLinear,
                                                                zcDom,
                                                                zcFor);
                return(vol);
            }

            catch(Exception& )
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "ComputeVolWithDeltaSpotAsInputWithPremium, did not converge (even in linear mode)");
            }
        }
    }

    return(vol);
}



double ComputeVolWithDeltaSpotAsInputWithoutPremium(double strike, double fxFwd,
                                                    double FXSpot,
                                                    double matu,
                                                    ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                                    ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                                    double pivotStrike,
                                                    int interpolType,
                                                    ARM_ZeroCurve* zcDom,
                                                    ARM_ZeroCurve* zcFor)
{
    double vol = 20.0;
    double precision = ARM_FX_GEN_VOL_PREC;
    double volPrec = 20;
    int nbIter = 0; 
    int nbIterMax = 100;
    double delta = 0;
    int callPut;

    double dfDom;
    double dfFor;


    // Compute discount factors

    ComputeFxDomForDiscountFactorsAtSettlmentDate(zcDom, zcFor,
                                                  matu,
                                                  // Output
                                                  dfDom,
                                                  dfFor);
    do 
    {
        if ( strike < pivotStrike )
        {
           callPut = K_PUT;
        }
        else
        {
           callPut = K_CALL;
        }

        vol = volPrec;

        delta = ComputeSpotFxDeltaWithoutPremium(vol/100.0, 
                                                 fxFwd, 
                                                 FXSpot, 
                                                 strike,
                                                 matu,
                                                 dfDom,
                                                 dfFor,
                                                 callPut);
      
        volPrec = ComputeSigmaFromDelta(delta, deltaCall, volsCall, deltaPut,
                                        volsPut, interpolType);
		
        nbIter++;

    } 
    while (( fabs(vol-volPrec) > precision ) && ( nbIter < nbIterMax ));

    if ( nbIter == nbIterMax )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ComputeVolWithDeltaSpotAsInputWithoutPremium");
    }

    return(vol);
}



void fooFWP(void)
{
    double dCall[100] = {0.1, 0.25, 0.46383155};
    double dPut[100]  = {-0.46383155, -0.25, -0.1};

    double vCall[100] = {11.128743, 10.53413295, 12.25};
    double vPut[100]  = {12.25, 14.93199262, 18.4364718};

    ARM_Vector deltaCall(3, dCall);
    ARM_Vector deltaPut(3, dPut);
    ARM_Vector volsCall(3, vCall);
    ARM_Vector volsPut(3, vPut);

	double strike, fxFwd, maturity ;
	strike   = 100.0;
	fxFwd    = 77.2083711;
	maturity = 10.0054795;
    

  
	int interp = K_LINEAR;
	double volZDS ;
	volZDS=12.25;
	double delta, deltaP;

	ARM_ComputeImpliedZDS_FWP(volZDS, 77.2083711, 10.0054795, strike, delta, deltaP);

	deltaCall.Elt(2)=delta;
	deltaPut.Elt(0)=-delta;

	//res = computeVol_withDeltaFWPasInput_aux(100, fxFwd, maturity, &deltaCall, &volsCall, 
    //                                         &deltaPut, &volsPut, strike, interp); 
}



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/