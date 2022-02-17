/****************************************************************************/
/*      Zero price and zero bank.                                           */
/****************************************************************************/
/*      ZEROS.c                                                             */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esl_zeros.h"
#include "esl_date.h"
#include "esl_util.h"
#include "esl_error.h"

/** 0=Linear Zero Cpn; 1=Flat Fwd */ 
int  ZeroInterpTypeFlag = 1;        
/** 0=Linear Stub; 1=Flat Stub    */ 
int  ZeroInterpTypeFlagStub = 1;    

/*****  ZeroPrice  ****************************************************/
/**
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      (ZeroInterpTypeFlag = 0) Linear Zero Cpn:
 *      -- linear zero cpn interp
 *      -- flat zero cpn extrapolation on both sides
 *
 *      (ZeroInterpTypeFlag = 1) Flat Fwd:
 *      -- flat fwd interp and extrapolation
 */

double   ZeroPrice(long         MatDate    /** (I) Mat date of the zero       */
                  ,long         ValueDate  /** (I) Value date                 */
                  ,int          NbZero     /** (I) Number of zeros in curve   */
                  ,long const*  ZeroDates  /** (I) maturity dates in zero crv */
                  ,double const*ZeroRates  /** (I) Zero rates                 */
		)
{
    double Price = -999.99;
    double t,  ZR, 
           t1, ZR1, Z1, 
           t2, ZR2, Z2;
    int    idx;

    /* basic checks */
    if ((ZeroDates == NULL) || (ZeroRates == NULL)) goto RETURN;
    if (NbZero <= 0) goto RETURN;
    if (MatDate == ValueDate) return(1.0);

    t   = Daysact(ValueDate, MatDate)/365.0;
    idx = GetDLOffset(NbZero, ZeroDates, MatDate, CbkHIGHER);

    /* MatDate <= 1st ZeroDate or crv has only 1 pt */
    /*  Flat Stub                                   */
    if ( ((idx == 0) || (NbZero == 1)) && (ZeroInterpTypeFlagStub == 1) ) 
    {
        t1  = 0.0;
        ZR1 = ZeroRates[0];

        t2  = Daysact(ValueDate, ZeroDates[0])/365.0;
        ZR2 = ZeroRates[0];
    }
    /* MatDate <= 1st ZeroDate or crv has only 1 pt */
    /*  Flat Stub                                   */
    else if ( ((idx == 0) || (NbZero == 1)) && (ZeroInterpTypeFlagStub == 0) ) 
    {
        t1  = Daysact(ValueDate, ZeroDates[0])/365.0;
        ZR1 = ZeroRates[0];

        t2  = Daysact(ValueDate, ZeroDates[1])/365.0;
        ZR2 = ZeroRates[1];
    }
    else if (idx<0) /* i.e. all zero dates are < MatDate, flat fwd extrap */
    {
        switch (ZeroInterpTypeFlag)
        {
        case 0: /* Linear Zero Cpn */
            t1  = Daysact(ValueDate, ZeroDates[NbZero-1])/365.0;
            ZR1 = ZeroRates[NbZero-1];

            t2  = t;
            ZR2 = ZeroRates[NbZero-1];
            break;
        
        case 1: /* Flat Fwd */    
            t1  = Daysact(ValueDate, ZeroDates[NbZero-2])/365.0;
            ZR1 = ZeroRates[NbZero-2];

            t2  = Daysact(ValueDate, ZeroDates[NbZero-1])/365.0;
            ZR2 = ZeroRates[NbZero-1];
            break;        
        
        default:
            goto RETURN;
        }
    }
    else
    {
        t1  = Daysact(ValueDate, ZeroDates[idx-1])/365.0;
        ZR1 = ZeroRates[idx-1];

        t2  = Daysact(ValueDate, ZeroDates[idx])/365.0;
        ZR2 = ZeroRates[idx];
    }
    
    switch (ZeroInterpTypeFlag)
    {
    case 0: /* Linear Zero Cpn */
        if (linterp(t, &ZR,
                    t1,  t2,
                    ZR1, ZR2) == FAILURE) goto RETURN;
                    
        Price = pow(1.0 + ZR, -t);
        break;
        
    case 1: /* Flat Fwd */
        if (IS_EQUAL(t1, t2)) goto RETURN;
        Z1  = pow(1.0 + ZR1, -t1);
        Z2  = pow(1.0 + ZR2, -t2);
        Price = Z1 * pow(Z2/Z1, (t-t1)/(t2-t1));
        break;
        
    default:
        goto RETURN;
    }

RETURN:

    if (Price < 0.0)
    {
        DR_Error("ZeroPrice: failed.");
    }
    return (Price);

} /* ZeroPrice */


/*****	Get_Zero  *************************************************************/
/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *	   curve stored in t_curve style (deterministic).
 *
 */
int  Get_Zero(double   *OutZeroRate   /** (O) Zero rate	                  */
             ,double   *OutZeroPrice  /** (O) Zero price                  */
             ,int       NbZero        /** (I) Number of zeros in the ZC   */
             ,double   *Zero          /** (I) Z rates on an ACT/365 basis */
             ,long     *ZeroDate      /** (I) Zero maturity dates         */
             ,long      CurrentDate   /** (I) Current date                */
             ,long      Maturity      /** (I) Maturity of the zero        */
	     )
{
    double t, Pr;

    t = Daysact (CurrentDate, Maturity) / 365.0;

    Pr = ZeroPrice(Maturity,   
		   CurrentDate, 
		   NbZero,
		   ZeroDate, 
		   Zero);

    *OutZeroPrice = Pr;

    if (fabs(*OutZeroPrice -999.99) < TINY)
	return (FAILURE);

    if (t <= 0.0) 
    {
	    *OutZeroRate = Zero[0];
    }
    else 
    {
	    *OutZeroRate = pow(Pr, -1.0/t ) - 1.0;
    }

    return (SUCCESS);

}  /* Get_Zero */

