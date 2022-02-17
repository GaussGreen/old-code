/*
 * $Log: bondmath.cpp,v $
 * Revision 1.3  2003/07/10 07:38:11  ebenhamou
 * remove unused var
 *
 * Revision 1.2  2003/05/06 09:57:51  mab
 * RCS Comments
 *
 */

/*----------------------------------------------------------------------------*     
     bondmath.cpp
 
     This file contains analytic functions for doing bond math

*----------------------------------------------------------------------------*/
 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "expt.h"


#include "bondmath.h"








/*----------------------------------------------------------------------------*     

    NAME    accruedInterest
    
    Compute accrued interest

    Input:
        double coupon    - the coupon in percent
        int    frequency - the frequency of the bond (2 for semi annual, etc...)
        double previousSettl - the year term from previous coupon 
                               date to settlement 

    Output:
        the function returns the accrued interest

*----------------------------------------------------------------------------*/

double accruedInterest(double coupon, int frequency, double previousToSettl)
{
    if ( frequency <= 0 || 12%frequency != 0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                     "Invalid coupon frequency");
    }

    return(coupon*previousToSettl);
}




/*----------------------------------------------------------------------------*     

    NAME    regularPresentValueToYield
    
    Compute yield from present value in the case of a regular coupon period

    Input:
        double    presentValue    -the present value
        double    redemptionValue -the redemption value at maturity (usually 100)
        double    coupon    - the coupon in percent
        int       frequency - the frequency of the bond (2 for semi annual, etc...)
        double    settlToNext  - the year term from settlement to next coupon date
        double    settlToMaturity - the year term from settlement to maturity

    Output:
        the function returns the yield
*----------------------------------------------------------------------------*/

double regularPresentValueToYield(double presentValue, 
                                  double redemptionValue, 
                                  double coupon, 
                                  int    frequency,
                                  double settlToNext, 
                                  double settlToMaturity)
{
    double approxVector[2],
           minYield, maxYield, approxYield,
           f, dfdx, fmin, fmax;

    void* fixedParams[6];

    /*
     * if(presentValue*coupon < 0.0 || presentValue*redemptionValue <= 0.0 )
     * return( K_HUGE_DOUBLE );
     *  Test to see if we are within last coupon period.  In this case
     *  we just want simple interest.
     */

    if ( frequency * settlToMaturity <= 1.0 )
    {
       if ( settlToMaturity < 1.0e-12 ) 
       {
          return(0.0);
       }

       approxYield = redemptionValue+coupon/(double)frequency-presentValue;

       return(100.0*approxYield/(presentValue*settlToMaturity));
    }

    /*
     *  Take care of the zero coupon case.
     */

    if ( coupon == 0.0 )
    {
       approxYield = (double) frequency*(pow((redemptionValue/presentValue),
                          1.0/((double)frequency*settlToMaturity))-1.0);

       return( 100.0 * approxYield );
    }

    /*
     *  Obtain the initial approximation to the yield (in approxVector[0])
     *  and the approximate error (in approxVector[1]).
     */

    if ( approxRegularPvToYield(presentValue, redemptionValue, 
             coupon, frequency,
            settlToNext, settlToMaturity, approxVector ) == (double *)NULL )
    {
       throw Exception(__LINE__, __FILE__, ERR_INITIAL_VALUE_PB,
                             "Initial Yield approximation problem");
    }

    fixedParams[0] = &presentValue;
    fixedParams[1] = &redemptionValue;
    fixedParams[2] = &coupon;
    fixedParams[3] = &frequency;
    fixedParams[4] = &settlToNext;
    fixedParams[5] = &settlToMaturity;

    /*
     *  Try to bracket the root by using the approximation +/- five times
     *  the estimated error.  Five is just pulled out of the hat using the
     *  assumption that any error estimate worth anything should be correct
     *  to within an order of magnitude.
     */

    minYield = approxVector[0]-5.0*approxVector[1];
    maxYield = approxVector[0]+5.0*approxVector[1];
    
    fmin = dummyRegularYieldToPresentValue(f, dfdx, minYield, fixedParams);

    if ( fabs(fmin) < 2.0 * K_DOUBLE_TOL ) 
    {
       return(minYield);
    }

    fmax = dummyRegularYieldToPresentValue(f, dfdx, maxYield, fixedParams);

    if ( fabs(fmax) < 2.0 * K_DOUBLE_TOL ) 
    {
       return(maxYield);
    }

    
    if (fmin * fmax < 0.0)
    {    
        /* root is correctly bracketed */

        try
        {
            approxYield = newtonRoot(dummyRegularYieldToPresentValue, minYield,
                        maxYield, fixedParams, kYieldErr, kIterMax);
        }

        catch(Exception& m)
        {
            m.DebugPrint();
            throw Exception(__LINE__, __FILE__, ERR_NEWTON_ROOT_PB,
                             "Problem in <newtonRoot> function");
        }
    }
    else 
    {
        /*
           * Root is not bracketed. This should happen only for very high or very low
           * (i.e., negative) yields.  In these cases the approximation can
           * be very bad.  However, the approximation will always be negative
           * if the real yield is negative and vice versa.  Thus it is used
           * to set the bracketing boundaries to [-99,0] for negative yields
           * and [0,kMaxYield] for positive ones.
         */

        minYield    = ( approxVector[0] < 0.0 ) ? -99.0: 0.0;
        maxYield    = ( approxVector[0] < 0.0 ) ?   0.0: kMaxYield;
        
        f = dummyRegularYieldToPresentValue(f, dfdx, minYield, fixedParams);

        while (( f >= K_HUGE_DOUBLE && dfdx >= K_HUGE_DOUBLE )
               && minYield < maxYield-1.0 )
        {
            minYield += 1.0;

            f = dummyRegularYieldToPresentValue(f, dfdx, minYield, fixedParams);
        }

        try
        {
            approxYield = newtonRoot(dummyRegularYieldToPresentValue, minYield,
                        maxYield, fixedParams, kYieldErr, kIterMax);
        }

        catch(Exception& m)
        {
            m.DebugPrint();

            throw Exception(__LINE__, __FILE__, ERR_NEWTON_ROOT_PB,
                             "Problem in <newtonRoot> function");
        }
    }

    return(approxYield);
}



/*----------------------------------------------------------------------------*


    NAME    approxRegularPvToYield
    
    Compute approximate yield from present value in the case of a regular
    coupon period

    Input:
        double presentValue    - the present value
        double redemptionValue - the redemption value at maturity (usually 100)
        double coupon     - the coupon in percent
        int    frequency  - the frequency of the bond (2 for semi annual, etc...)
        double settlToNext     - the year term from settlement to next coupon date
        double settlToMaturity - the year term from settlement to maturity

    Output:
        double *approxVector - The first value in this array (approxVector[0]) 
                is the approximate yield in units appropriate to the frequency 
                of the coupon. The second element (approxVector[1]) is an 
                estimate of the error in the approximation.
                                     
        The function returns a pointer to approxVector

*----------------------------------------------------------------------------*/

double* approxRegularPvToYield(
    double    presentValue, 
    double    redemptionValue, 
    double    coupon, 
    int    frequency,
      double    settlToNext, 
      double    settlToMaturity, 
      double    *approxVector )
{
    double    previousToSettl, price, testprice,
            kappa, kappa0, kappa1, kappa2,
            numerator, denominator;
    int        num;

    settlToMaturity  *= (double)frequency;
    settlToNext      *= (double)frequency;
    presentValue    *= 0.01;
    redemptionValue *= 0.01;
    coupon           *= 0.01 / (double)frequency;

    previousToSettl   = 1.0 - settlToNext;

     /*
      *   num = number of payments left.
      *
      *  The 0.00001 is to correctly account for the case where the settlement
      *  date falls on a coupon date.  It will not affect other cases since
      *  0.00001 is much smaller than 1/365 (i.e., a change of one day in a
      *  one year period).
      */

    num = (int)( settlToMaturity - settlToNext + 0.00001 ) + 1;

     /*
      *  Get price from the present value and the zero yield
      *  price (the sum of the undiscounted cash flows).
      */

    price     = presentValue - coupon * previousToSettl;
    testprice = redemptionValue + coupon * ( (double)num - previousToSettl );

    if ( fabs(price) < K_DOUBLE_TOL  ||  fabs(coupon) < K_DOUBLE_TOL )
    {
        return( (double *)NULL );
    }

    if ( price < testprice )
    {
        /*
         *  positive yield
         */

        kappa        = 1.0 + coupon / price;

        kappa0       = pow( kappa, -settlToMaturity );
        kappa1       = pow( kappa, -(double)(num-1) );
        kappa2       = pow( kappa, -settlToNext );

        denominator  = ( kappa - kappa1 ) * ( settlToNext + price / coupon ) -
                       (double)num * kappa1;

        denominator *= price * kappa2;
        denominator += redemptionValue * settlToMaturity * kappa0;

        numerator    = price * kappa2 * ( kappa - kappa1 );
        numerator   += redemptionValue * kappa0 - presentValue;

        if ( fabs(denominator) < K_DOUBLE_TOL )
        {
            return((double *) NULL);
        }

        approxVector[0] = 100.0 * frequency 
                           * ( coupon / price + numerator / denominator);

        approxVector[1] = 100.0 * frequency 
                          * coupon * numerator / ( price * denominator );
    }
    else
    {
        /*
         *  negative yield
         */

        approxVector[0] = 100.0 * frequency * ( testprice - price ) /
                  ( settlToMaturity * redemptionValue 
                  + coupon * (double)num * (settlToNext +0.5 * (double)(num-1)));

        approxVector[1] = 5.0;
    }

    return(approxVector);
}


/*----------------------------------------------------------------------------*

  A dummy function that just calls <regularYieldToPresentValue>
  using the current value of the yield and subtracts the input
  present value (because we want to find the zero of the function).

*----------------------------------------------------------------------------*/

double dummyRegularYieldToPresentValue(double& f, double& dfdx, 
                                       double x1, void **fixedParams )
{

    f = regularYieldToPresentValue( x1, 
      *(double *) fixedParams[1],
      *(double *) fixedParams[2],
      *(int *) fixedParams[3],
      *(double *) fixedParams[4],
      *(double *) fixedParams[5])
      - * (double *) fixedParams[0];
      
    dfdx = regularDpDy( x1,
      *(double *) fixedParams[1],
      *(double *) fixedParams[2],
      *(int *) fixedParams[3],
      *(double *) fixedParams[4],
      *(double *) fixedParams[5]);
      
    return(f);
}



/*----------------------------------------------------------------------------*

    NAME    regularYieldToPresentValue
    
    Compute present value from yield in the case of a regular coupon period

    Input:
        double    yield        - the yield in percent (e.g. 9.50)
        double    redemptionValue - the redemption value at maturity (usually 100)
        double    coupon    - the coupon in percent
        int    frequency    - the frequency of the bond (2 for semi annual, etc...)
        double    settlToNext  - the year term from settlement to next coupon date
        double    settlToMaturity    - the year term from settlement to maturity

    Output:
        the function returns the present value

*----------------------------------------------------------------------------*/

double regularYieldToPresentValue( 
    double    yield, 
    double    redemptionValue, 
    double    coupon, 
    int        frequency,
    double    settlToNext, 
    double    settlToMaturity )
{
    double pv, numerator, idenominator, factor, temp;
    int    num;


    if ( yield < -99.0  || redemptionValue * coupon <  0.0)
    {
        return( K_HUGE_DOUBLE );
    }

    temp              = 0.01 / (double)frequency;
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    redemptionValue *= 0.01;
    coupon           *= temp;
    yield            *= temp;

    /*
     *  num = number of payments left.
     * 
     *  The 0.00001 is to correctly account for the case where the settlement
     *  date falls on a coupon date.  It will not affect other cases since
     *  0.00001 is much smaller than 1/365 (i.e., a change of one day in a
     *  one year period).
     */

    num    = (int)( settlToMaturity - settlToNext + 0.00001 ) + 1;
    factor = 1.0 + yield;

    /*
      *  Check to see:
      *     if    1) we are within the last coupon period (use simple
      *              interest discounting)
      *  or if    2) the yield is very small
      *  or if    3) there is a large discounting factor
      *  or if    4) it is a normal case
      *  else     5) must be huge negative yield and a long term.
     */

    if ( settlToMaturity <= 1.0 )
    {
        if ( settlToMaturity < 1.0e-12 ) 
        {
            pv = redemptionValue;
        }
        else 
        {
            pv = (redemptionValue + coupon) / (1.0 + yield * settlToMaturity);
        }
    }
    else if( fabs(yield) < 1.0e-10 )
    {
        pv = redemptionValue + coupon * (double) num;
    }
    else if ( (numerator = pow(factor, (double) num)) > kTooBig )
    {
        pv = coupon * pow(factor,((double)num-settlToMaturity)) / yield;
    }
    else if ( (idenominator = pow( factor, -settlToMaturity )) < kTooBig )
    {
        pv = idenominator * (redemptionValue + coupon * (numerator-1.0) / yield);
    }
    else
    {
        return( K_HUGE_DOUBLE );
    }

    return( 100.0 * pv );
}



/*----------------------------------------------------------------------------*

    NAME    regularDpDy
    
    Compute derivative of price with respect to yield in the case of a 
    regular coupon period
    
    Input:
        double yield        - the yield in percent
        double redemptionValue  - the redemption value at maturity (usually 100)
        double coupon - the coupon in percent
        int    frequency - the frequency of the bond (2 for semi annual, etc...)
        double settlToNext    - the year term from settlement to next coupon date
        double settlToMaturity    - the year term from settlement to maturity

    Output:
        the function returns the derivative of price with respect to yield

*----------------------------------------------------------------------------*/

double regularDpDy(
    double    yield, 
    double    redemptionValue, 
    double    coupon, 
    int    frequency,
      double    settlToNext, 
      double    settlToMaturity )
{
    double deriv,
    numerator,
    idenominator,
    factor,
    temp;

    int    num;

    if ( yield < -99.0  || redemptionValue * coupon <  0.0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Invalid Yield or Redemption value");
    }

    temp              = 0.01 / (double)frequency;
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    redemptionValue *= 0.01;
    coupon           *= temp;
    yield            *= temp;

    /*
     *  num = number of payments left.
      *  The 0.00001 is to correctly account for the case where the settlement
         *  date falls on a coupon date.  It will not affect other cases since
         *  0.00001 is much smaller than 1/365 (i.e., a change of one day in a
      *  one year period).
     */

    num    = (int)( settlToMaturity - settlToNext + 0.00001 ) + 1;
    factor = 1.0 + yield;

    /*
     *  Check to see:
     *     if    1) we are within the last coupon period (use
     *              derivative of simple interest discounting)
     *  or if    2) the yield is very small
     *  or if    3) there is a large discounting factor
     *  or if    4) it is a normal case
     *  else     5) must be huge negative yield and a long term.
     */

    if ( settlToMaturity <= 1.0 )
    {
        if ( settlToMaturity < 1.0e-12 ) 
        {
            deriv  = 0.0;
        }
        else 
        {
            deriv  = -settlToMaturity * ( redemptionValue + coupon );
            deriv *=  pow( (1.0+yield*settlToMaturity), -2.0 );
        }
    }
    else if ( fabs(yield) < 1.0e-6 )
    {
        deriv = -0.5 * coupon * num * ( num - 1.0 + 2.0 * settlToNext ) -
             redemptionValue * settlToMaturity;
    }
    else if( (numerator = pow( factor, (double)num )) > kTooBig )
    {
        deriv  = -coupon * ( 1.0 / yield + settlToNext ) / yield;
        deriv *=  pow( factor, ((double)(num-1) - settlToMaturity) );
    }
    else if( (idenominator = pow( factor, -(1.0+settlToMaturity) )) < kTooBig )
    {
        temp   =  1.0 / yield;
        deriv  =  coupon * ( temp * factor - numerator * ( temp + settlToNext ) ) +
             settlToMaturity * ( coupon - yield * redemptionValue );
        deriv *=  temp * idenominator;
    }
    else
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Calculation problem in <regularDpDy>");
    }

    return( deriv / (double)frequency );
}



/*----------------------------------------------------------------------------*
    NAME    regularD2pDy2
    
    Compute second derivative of price with respect to yield in the case of 
    a regular coupon period
    
    Input:
        double yield    - the yield in percent
        double redemptionValue    - the redemption value at maturity (usually 100)
        double coupon       - the coupon in percent
        int    frequency    - the frequency of the bond (2 for semi annual, etc...)
        double settlToNext    - the year term from settlement to next coupon date
        double settlToMaturity    - the year term from settlement to maturity

    Output:
        the function returns the second derivative of price with respect to yield
*----------------------------------------------------------------------------*/

double regularD2pDy2( 
    double    yield, 
    double    redemptionValue, 
    double    coupon, 
    int        frequency,
    double    settlToNext, 
    double    settlToMaturity )
{
    double    deriv,numerator, idenominator, factor, temp1, temp;
    int    num, i;

    if( yield < -99.0  || redemptionValue*coupon < 0.0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Invalid Yield or Redemption value");
    }

    temp              = 0.01 / (double)frequency;
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    redemptionValue *= 0.01;
    coupon           *= temp;
    yield            *= temp;

    /*
      *  num = number of payments left.
      *  The 0.00001 is to correctly account for the case where the settlement
      *  date falls on a coupon date.  It will not affect other cases since
      *  0.00001 is much smaller than 1/365 (i.e., a change of one day in a
      *  one year period).
     */

    num    = (int)( settlToMaturity - settlToNext + 0.00001 ) + 1;
    factor = 1.0 + yield;

    /*
       Check to see:
     
         if    1) we are within the last coupon period (use 
                   derivative of simple interest discounting)
     
           or if    2) the yield is very small
     
          or if    3) there is a large discounting factor
    
          or if    4) it is a normal case
    
         else     5) must be huge negative yield and a long term.
    */

    if ( settlToMaturity <= 1.0 )
    {
        if( settlToMaturity < 1.0e-12 )
        {
            deriv  = 0.0;
        }
        else
        {
            deriv  = 2.0 * settlToMaturity * settlToMaturity *
                ( redemptionValue + coupon );
            deriv *= pow( ( 1.0 + yield * settlToMaturity ), -3.0 );
        }
    }
    else if( fabs(yield) < 1.0e-4 )
    {
        deriv = 0.0;

        for ( i = 1; i <= num; i++ ) 
        {
            deriv += (double)i * (double)i;
        }
        deriv += (double)num * ( (double)(num+1) * ( settlToNext - 0.5 ) +
            settlToNext * ( settlToNext - 1.0 ));

        deriv *= coupon;
        deriv += settlToMaturity * ( 1.0 + settlToMaturity ) * redemptionValue;
    }
    else if( (numerator = pow( factor, (double)num )) > kTooBig )
    {
        temp   = 1.0 / yield;
        temp1  = temp * factor;
        deriv  = settlToNext * ( 1.0 + settlToNext ) +
            2.0 * temp * ( temp1 + settlToNext );
        deriv *= temp * coupon *
            pow( factor, ( (double)(num-2) - settlToMaturity ) );
    }
    else if( (idenominator = pow( factor, -(2.0+settlToMaturity) )) < kTooBig )
    {
        temp   = 1.0 / yield;
        temp1  = temp * factor;
        deriv  = numerator * ( settlToNext * ( 1.0 + settlToNext ) +
            2.0 * temp * ( temp1 + settlToNext )
            );
        deriv -= settlToMaturity * ( settlToMaturity + 1.0 + 2.0 * temp1 ) +
            2.0 * temp1 * temp1;
        deriv *= coupon;
        deriv += yield * settlToMaturity * ( 1.0 + settlToMaturity ) *
            redemptionValue;
        deriv *= temp * idenominator;
    }
    else
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Calculation problem in <regularD2pDy2>");
    }

    return( 0.01 * deriv / ( (double)frequency * (double)frequency ) );
}




/*----------------------------------------------------------------------------*

NAME    regularForwardValue
    
    Compute forward value (e.g. forward price + accrued at forward date) of a regular coupon bond

Input:
    double    presentValue    - present value of bond
    double    reinvestRate    - the reinvestment rate in percent (e.g. 9.50). This rate has 
            the same unit than the underlying yield (e.g. semi-annual
                        if semi-annual bond, etc...)
    double    redemptionValue        - the redemption value at maturity (usually 100)
    double    coupon            - the coupon in percent
    int       frequency        - the frequency of the bond (2 for semi annual, etc...)
    double    settlToNext            - the year term from settlement to next coupon date
    double    settlToForward        - the year term from settlement to forward date
    int     mmFlag          - 1 (default) if MM  method is used to compute cost of carry (<1year)

    Output:
        the function returns the forward value.
*----------------------------------------------------------------------------*/

double regularForwardValue( 
    double    presentValue, 
    double    reinvestRate, 
    double    coupon, 
    int       frequency,
    double    settlToNext, 
    double    settlToForward,
    int mmFlag)
{
    double    fvOfCoupons, invCoupon,
        forwardValue, lastToForward, temp;
    int    num, i;
    
    temp            = 0.01 / (double)frequency;
    settlToNext     *= (double)frequency;
    settlToForward  *= (double)frequency;
    coupon          *= temp;
    reinvestRate    *= temp;

     /*
       num = number of payments left.
     
       The 0.00001 is to correctly account for the case where the settlement
       date falls on a coupon date.  It will not affect other cases since
       0.00001 is much smaller than 1/365 (i.e., a change of one day in a
        one year period).
     */

    num = (int) floor( settlToForward - settlToNext + 0.00001 ) + 1;
    lastToForward = settlToForward - settlToNext - floor( settlToForward - settlToNext + 0.00001 );


    if (settlToForward < settlToNext) 
    {
        fvOfCoupons = 0.0;
    }
    else if( fabs(reinvestRate) < 1.0e-10 )
    {
        fvOfCoupons = coupon * (double) num;
    }
    else if (settlToForward <= 1.0) 
    {
        fvOfCoupons = coupon * ((double) num 
          + reinvestRate * lastToForward * (double) num
          + 0.5 * reinvestRate * (double) num * (double) (num-1));
    }
/*
    else if( (numerator = pow( (1.0+reinvestRate), (double)num )) < kTooBig )
    {
        fvOfCoupons =  coupon * ( numerator - 1.0 ) / reinvestRate;
        fvOfCoupons *= pow(1.0+reinvestRate, lastToForward);
        //fvOfCoupons *= 1.0+reinvestRate*lastToForward;
    }
*/
    else if( num>0 )
    {
        fvOfCoupons = 0.0;

        for (i=num; i>0; i--)
        {
            invCoupon = coupon*pow(1.0+reinvestRate, lastToForward);
            fvOfCoupons += invCoupon;
            lastToForward += (double) frequency;
        }
    }
    else
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Calculation problem in <regularForwardValue>");
    }

    fvOfCoupons *= 100.0;


    if (mmFlag) 
    {
        forwardValue = presentValue * (1.0 + reinvestRate * settlToForward) - fvOfCoupons;
    }
    else
    {
        forwardValue = presentValue * pow(1.0 + reinvestRate, settlToForward) - fvOfCoupons;
//        forwardValue = presentValue * (1.0 + reinvestRate * settlToForward) - fvOfCoupons;
    }

    return(forwardValue);
}




/*----------------------------------------------------------------------------*

    NAME    regularDfDr
    
    Compute first derivative of forward value  eith respect to reinvest rate
    in the case of a regular coupon bond

    SYNOPSIS
        double regularDfDr( 
            double    presentValue, 
            double    reinvestRate, 
            double    coupon, 
            int        frequency,
              double    settlToNext, 
              double    settlToForward )

    
    Input:
        double    presentValue    - present value of bond
        double    reinvestRate    - the reinvestment rate in percent (e.g. 9.50). This rate has 
                    the same unit than the underlying yield (e.g. semi-annual
                    if semi-annual bond, etc...)
        double    redemptionValue    - the redemption value at maturity (usually 100)
        double    coupon        - the coupon in percent
        int    frequency    - the frequency of the bond (2 for semi annual, etc...)
        double    settlToNext    - the year term from settlement to next coupon date
        double    settlToForward    - the year term from settlement to forward date

    Output:
        the function returns the first derivatice of the forward value wrt reinvest rate
*----------------------------------------------------------------------------*/

double regularDfDr(
    double    presentValue, 
    double    reinvestRate, 
    double    coupon, 
    int        frequency,
      double    settlToNext, 
      double    settlToForward )
{
    double    lastToForward, numerator, dfvOfCoupons, dforwardValue, temp;
    int    num;
    
    temp             = 0.01 / (double)frequency;
    settlToNext     *= (double)frequency;
    settlToForward  *= (double)frequency;
    coupon          *= temp;
    reinvestRate    *= temp;

     /*
       num = number of payments left.
     
       The 0.00001 is to correctly account for the case where the settlement
       date falls on a coupon date.  It will not affect other cases since
       0.00001 is much smaller than 1/365 (i.e., a change of one day in a
        one year period).
     */

    num = (int) floor( settlToForward - settlToNext + 0.00001 ) + 1;
    lastToForward = settlToForward - settlToNext - floor( settlToForward - settlToNext + 0.00001 );


    if (settlToForward < settlToNext) 
    {
        dfvOfCoupons = 0.0;
    }
    else if( fabs(reinvestRate) < 1.0e-10 )
    {
        dfvOfCoupons = 0.0;
    }
    else if( (numerator = pow( (1.0+reinvestRate), (double)num )) < kTooBig )
    {
        dfvOfCoupons = ((double) num) * pow(1.0+reinvestRate, (double) num -1.0 + lastToForward) 
                  / reinvestRate;

        dfvOfCoupons += lastToForward * (numerator - 1.0) * pow(1.0+reinvestRate, 
                lastToForward - 1.0);

        dfvOfCoupons -= (numerator - 1.0) * pow(1.0+reinvestRate, lastToForward) 
                  / (reinvestRate * reinvestRate);

        dfvOfCoupons *=  coupon;
    }
    else
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Calculation problem in <regularDfDr>");
    }
    
    dfvOfCoupons *= 100.0;


    dforwardValue = presentValue * settlToForward * pow(1.0 + reinvestRate, settlToForward-1.0) 
                  - dfvOfCoupons;
      
    dforwardValue *= temp;
      
    return(dforwardValue);
}



/*----------------------------------------------------------------------------*

    NAME    dummyRegularForwardValue
    
    A dummy function that calls regularForwardValue and regularDfDr, to be used 
    with newtonRoot for computing reinvest rate from forward value

*----------------------------------------------------------------------------*/

double dummyRegularForwardValue(double &f, double &df, double reinvestRate, void **fixedParams)
{
    double    presentValue, coupon, settlToNext, settlToForward, forwardValue;
    int frequency;

    presentValue = * (double *) fixedParams[0];
    coupon = * (double *) fixedParams[1];
    frequency = * (int *) fixedParams[2];
    settlToNext = * (double *) fixedParams[3];
    settlToForward = * (double *) fixedParams[4];
    forwardValue = * (double *) fixedParams[5];


    f = regularForwardValue( presentValue, reinvestRate, coupon, frequency,
       settlToNext, settlToForward ) - forwardValue;
 
     df = regularDfDr( presentValue, reinvestRate, coupon, frequency, settlToNext, settlToForward );
     
     return(f);
}
 


/*----------------------------------------------------------------------------*
    NAME    regularReinvestRate
    
    Compute repo rate from forward value in the case of a regular coupon bond

    Input:
        double    presentValue        - present value of bond
        double    forwardValue        - forward value of bond
        double    redemptionValue        - the redemption value at maturity (usually 100)
        double    coupon            - the coupon in percent
        int    frequeny        - the frequency of the bond (2 for semi annual, etc...)
        double    settlToNext        - the year term from settlement to next coupon date
        double    settlToForward        - the year term from settlement to forward date

    Output:
        the function returns the reinvest rate in the same unit than the bond yield.

*----------------------------------------------------------------------------*/

double regularReinvestRate(double presentValue, 
                            double forwardValue, 
                            double coupon, 
                            int frequency,
                            double settlToNext, 
                            double settlToForward,
                            int mmFlag)
{
    double    r1, r2, reinvestRate, f, df;
    void    *fixedParams[6];

    fixedParams[0] = (void *) &presentValue;
    fixedParams[1] = (void *) &coupon;
    fixedParams[2] = (void *) &frequency;
    fixedParams[3] = (void *) &settlToNext;
    fixedParams[4] = (void *) &settlToForward;
    fixedParams[5] = (void *) &forwardValue;


    /*
     * Bracket root for reinvest rate
     */
    
    r1 = 3.0;
    while (dummyRegularForwardValue(f, df, r1, fixedParams) > 0.0 && r1 > kMinYield) 
    {
        r1 -= 1.0;
    }

    if (r1 <= kMinYield)
    {
        throw Exception(__LINE__, __FILE__, ERR_INITIAL_VALUE_PB,
                             "Cannot bracket root for reinvest rate");
    }

    r2 = 30.0;

    while (dummyRegularForwardValue(f, df, r2, fixedParams) > 0.0 && r2 < kMaxYield) 
    {
        r2 += 1.0;
    }

    if (r1 >= kMaxYield)
    {
        throw Exception(__LINE__, __FILE__, ERR_INITIAL_VALUE_PB,
                             "Cannot bracket root for reinvest rate");
    }

    try
    {
        reinvestRate = newtonRoot(dummyRegularForwardValue, r1, r2,
        fixedParams, kYieldErr, kIterMax);
    }

    catch(Exception& m)
    {
        m.DebugPrint();
        throw Exception(__LINE__, __FILE__, ERR_NEWTON_ROOT_PB,
                             "Problem in <newtonRoot> function");
    }

    return(reinvestRate);
}



/*----------------------------------------------------------------------------*

    NAME    oddFirstPresentValueToYield
    
    Compute yield from present value in the case of an odd first coupon and 
    settlement is within odd period
    
    Input:
        double    yield            - the yield in percent
        double    redemptionValue        - the redemption value at maturity (usually 100)
        double    coupon            - the coupon in percent
        int    frequency        - the frequency of the bond (2 for semi annual, etc...)
        double    issueToNext        - the year term from issue to first coupon
        double    settlToNext        - the year term from settlement to first coupon date
        double    settlToMaturity        - the year term from settlement to maturity

    Output:
        the function returns the yield if ok, K_HUGE_DOUBLE if not.

*----------------------------------------------------------------------------*/

double oddFirstPresentValueToYield(
    double    presentValue, 
    double    redemptionValue, 
    double    coupon, 
    int        frequency,
     double    issueToNext, 
     double    settlToNext, 
     double    settlToMaturity )
{
    double    approxVector[2],
        minYield, maxYield,fmin, fmax, approxYield, f, dfdx;
    void    *fixedParams[7];

    if(presentValue*coupon < 0.0 || presentValue * redemptionValue <= 0.0)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Negative value of PV or coupon or RedValue");
    }

    /*
      *  Test to see if we are within last coupon period.  In this case
      *  we just want simple interest.
     */

    if( frequency * settlToMaturity <= 1.0 )
    {
        if( settlToMaturity < 1.0e-12 ) return( 0.0 );
        approxYield = redemptionValue + coupon/(double)frequency - presentValue;
        return( 100.0 * approxYield / ( presentValue * settlToMaturity ) );
    }

    /*
     *  Obtain the initial approximation to the yield (in approxVector[0])
     *  and the appoximate error (in approxVector[1]).
     */

    if ( approxRegularPvToYield( presentValue, redemptionValue, coupon, frequency,
        settlToNext, settlToMaturity, approxVector ) == (double *)NULL )
    {
        return( K_HUGE_DOUBLE );
    }

    fixedParams[0] = &presentValue;
    fixedParams[1] = &redemptionValue;
    fixedParams[2] = &coupon;
    fixedParams[3] = &frequency;
    fixedParams[4] = &issueToNext;
    fixedParams[5] = &settlToNext;
    fixedParams[6] = &settlToMaturity;

    /*
     *  Try to bracket the root by using the approximation +/- five times
     *  the estimated error.  Five is just pulled out of the hat using the
     *  assumption that any error estimate worth anything should be correct
     *  to within an order of magnitude.
     */

    minYield       = approxVector[0] - 5.0 * approxVector[1];
    maxYield       = approxVector[0] + 5.0 * approxVector[1];

    fmin = dummyOddFirstYieldToPresentValue(f, dfdx, minYield, fixedParams);

    if (fabs(fmin) < 2.0 * K_DOUBLE_TOL) 
    {
        return(minYield);
    }

    fmax = dummyOddFirstYieldToPresentValue(f, dfdx, maxYield, fixedParams);

    if (fabs(fmax) < 2.0 * K_DOUBLE_TOL) 
    {
        return(maxYield);
    }

    
    if (fmin * fmax < 0.0)
    {    /*    root is correctly bracketed    */

        try
        {
            approxYield = newtonRoot(dummyOddFirstYieldToPresentValue, minYield,
                maxYield, fixedParams, kYieldErr, kIterMax);
        }

        catch(Exception& m)
        {
            m.DebugPrint();
            throw Exception(__LINE__, __FILE__, ERR_NEWTON_ROOT_PB,
                             "Problem in <newtonRoot> function");
        }
    }
    else 
    {
           /*
            *  newtonRoot will return the value of K_HUGE_DOUBLE if the root was not
            *  bracketed.  This should happen only for very high or very low
            *  (i.e., negative) yields.  In these cases the approximation can
            *  be very bad.  However, the approximation will always be negative
            *  if the real yield is negative and vice versa.  Thus it is used
            *  to set the bracketing boundaries to [-99,0] for negative yields
            *  and [0,kMaxYield] for positive ones.
            */

        minYield    = ( approxVector[0] < 0.0 ) ? -99.0: 0.0;
        maxYield    = ( approxVector[0] < 0.0 ) ?   0.0: kMaxYield;
        
        
        f = dummyOddFirstYieldToPresentValue(f, dfdx, minYield, fixedParams);

        while ( (f >= K_HUGE_DOUBLE && dfdx >= K_HUGE_DOUBLE) && minYield < maxYield-1.0) 
        {
            minYield += 1.0;
            f = dummyRegularYieldToPresentValue(f, dfdx, minYield, fixedParams);
        }

        try
        {
            approxYield = newtonRoot(dummyOddFirstYieldToPresentValue, minYield,
                maxYield, fixedParams, kYieldErr, kIterMax);
        }
        catch(Exception& m)
        {
            m.DebugPrint();
            throw Exception(__LINE__, __FILE__, ERR_NEWTON_ROOT_PB,
                             "Problem in <newtonRoot> function");
        }
    }
    return( approxYield );
}




/*----------------------------------------------------------------------------*
   A dummy function that just calls <oddFirstYieldToPresentValue>
   using the current value of the yield and subtracts the input
   present value (because we want to find the zero of the function).
*----------------------------------------------------------------------------*/

double dummyOddFirstYieldToPresentValue(double &f, double &dfdx, double x1, void **fixedParams )
{
    f = oddFirstYieldToPresentValue( x1,
        * (double *) fixedParams[1],
        * (double *) fixedParams[2],
        * (int *) fixedParams[3],
        * (double *) fixedParams[4],
        * (double *) fixedParams[5],
        * (double *) fixedParams[6]   ) - * (double *) fixedParams[0];

    dfdx = oddFirstDpDy( x1,
        * (double *) fixedParams[1],
        * (double *) fixedParams[2],
        * (int *) fixedParams[3],
        * (double *) fixedParams[4],
        * (double *) fixedParams[5],
        * (double *) fixedParams[6] );

    return(f);
}



/*----------------------------------------------------------------------------*

    NAME    oddFirstYieldToPresentValue

    Compute present value from yield in the case of an odd first
    coupon and settlement in odd period
    
    Input:
        double    yield    - the yield in percent
        double    redemptionValue    - the redemption value at maturity (usually 100)
        double    coupon        - the coupon in percent
        int    frequency    - the frequency of the bond (2 for semi annual, etc...)
        double    issueToNext    - the year term from issue to first coupon
        double    settlToNext    - the year term from settlement to first coupon date
        double    settlToMaturity    - the year term from settlement to maturity

    Output:
        the function returns the present value
*----------------------------------------------------------------------------*/

double oddFirstYieldToPresentValue(
    double    yield, 
    double    redemptionValue, 
    double    coupon, 
    int    frequency,
     double    issueToNext, 
     double    settlToNext, 
     double    settlToMaturity )
{
    double pv, denominator,factor,temp;

    if ((yield < -99.0)  || (redemptionValue * coupon <  0.0))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Invalid Yield or Redemption value");
    }

    pv = regularYieldToPresentValue(yield, redemptionValue, coupon, frequency,
           settlToNext, settlToMaturity);

    temp              = 0.01 / (double)frequency;
    issueToNext     *= (double) frequency;
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    redemptionValue     *= 0.01;
    coupon           *= temp;
    yield            *= temp;

    factor = 1.0 + yield;
    denominator = pow(factor, settlToNext);
    temp = 100.0 * coupon * (issueToNext - 1.0) / denominator;
    pv += temp;
    

    return(pv);
}




/*----------------------------------------------------------------------------*

    NAME    oddFirstDpDy
    
    Compute derivative of price with respect to yield in the case of 
    an odd first
    coupon and settlement in odd period

    Input:
        double yield           - the yield in percent
        double redemptionValue - the redemption value at maturity (usually 100)
        double coupon    -the coupon in percent
        int    frequency -the frequency of the bond (2 for semi annual, etc...)
        double issueToNext - the year term from issue to first coupon
        double settlToNext - the year term from settlement to first coupon date
        double    settlToMaturity - the year term from settlement to maturity

    Output:
        the function returns the derivative of price with respect to yield

*----------------------------------------------------------------------------*/

double oddFirstDpDy(
     double    yield, 
     double    redemptionValue, 
     double    coupon, 
     int       frequency,
     double    issueToNext, 
     double    settlToNext, 
     double    settlToMaturity )
{
    double deriv, denominator,factor,temp;



    if (( yield < -99.0 ) || ( redemptionValue*coupon <  0.0 ))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Invalid Yield or Redemption value");
    }

    deriv = regularDpDy(yield, redemptionValue, coupon, 
             frequency, settlToNext, settlToMaturity );

    temp              = 0.01 / (double)frequency;
    issueToNext      *= (double)frequency;
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    coupon           *= temp;
    yield            *= temp;

    factor = 1.0 + yield;
    denominator = pow(factor, settlToNext + 1.0);

    temp = -coupon*(issueToNext-1.0)*settlToNext/(denominator
              *(double) frequency);

    deriv += temp;

    return(deriv);

}




/*----------------------------------------------------------------------------*

    NAME    oddFirstD2pDy2
    
    Compute second derivative of price with respect to yield in the case 
    of an odd first coupon and settlement in odd period

    Input:
     double yield            - the yield in percent
     double redemptionValue  - the redemption value at maturity (usually 100)
     double coupon    - the coupon in percent
     int    frequency - the frequency of the bond (2 for semi annual, etc...)
     double issueToNext - the year term from issue to first coupon
     double settlToNext - the year term from settlement to first coupon date
     double settlToMaturity - the year term from settlement to maturity

    Output:
        the function returns the second derivative of price with 
        respect to yield
*----------------------------------------------------------------------------*/


double oddFirstD2pDy2(
     double    yield, 
     double    redemptionValue, 
     double    coupon, 
     int       frequency,
     double    issueToNext, 
     double    settlToNext, 
     double    settlToMaturity )
{
    double deriv, numerator,factor,temp;



    if ((yield < -99.0)  || (redemptionValue * coupon <  0.0))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                             "Invalid Yield or Redemption value");
    }

    deriv = regularD2pDy2(yield, redemptionValue, coupon, 
                          frequency, settlToNext, settlToMaturity );

    temp              = 0.01 / (double)frequency;
    issueToNext      *= (double)frequency;
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    redemptionValue  *= 0.01;
    coupon           *= temp;
    yield            *= temp;
    
    factor = 1.0 + yield;
    numerator = pow(factor, -settlToNext - 2.0);

    temp = coupon*(issueToNext-1.0)*numerator*settlToNext
           *(settlToNext+1.0)*0.01/(double) frequency;

    deriv += temp;

    return(deriv);

}


#undef kMaxYield
#undef kMinYield
#undef kYieldErr
#undef kIterMax
#undef kTooBig





/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
