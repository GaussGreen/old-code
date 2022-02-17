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
#include "bmx123head.h"


int  ZeroInterpTypeFlag = 1;    /* 0=Linear Zero Cpn; 1=Flat Fwd */

/*****  Zero_t  *************************************************************/
/*
*       Zero coupon
*/
int     Zero_t (double      *Zero,         /* (I/O) Zero prices       */
                long        Reset,         /* (I) Reset flag          */
                int         t,             /* (I) Current time point  */
                int         T,             /* (I) Last time point     */
                int         DCurve,        /* (I) Discount curve      */
                DEV_DATA    *dev_data,     /* (I) Dev data structure  */
                TREE_DATA   *tree_data)    /* (I) Tree data structure */
{

    double  *ZeroL;                         /* Local slice pointer    */

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */

    int     i, j, k;                        /* Node indices           */
    int     status = FAILURE;               /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Reset)
    {
        if (tree_data->NbFactor == 1)
        {
            ZeroL = Zero + Node_Offset(1, 0, 0, t, tree_data);
    
            for (i = Bottom1; i <= Top1; i ++)                                      
            {
                ZeroL[i] = 1.;
            }
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                ZeroL = Zero + Node_Offset(2, i, 0, t, tree_data);
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    ZeroL[j] = 1.;
                }
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    ZeroL = Zero + Node_Offset(3, i, j, t, tree_data);
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        ZeroL[k] = 1.;
                    }
                }  /* for j */  
        }  /* if then else */
    }
    else
    {
        if (Dev (   Zero,
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Zero_t */



/*****  Zero_Bank  **********************************************************/
/*
*       Manipulation of the zero bank.
*/
int     Zero_Bank ( double      **Zero,        /* (I/O) Set of zeros        */
                    long        *ZeroMaturity, /* (I/O) Zero maturities     */
                    int         *NbZero,       /* (I/O) Current nb of zeros */
                    int         TotNbZero,     /* (I) Total number of zeros */
                    long        Reset,         /* (I) Reset flag            */
                    long        CurrentDate,   /* (I) Current date          */
                    int         t,             /* (I) Current time point    */
                    int         T,             /* (I) Last time point       */
                    int         DCurve,        /* (I) Discount curve        */
                    DEV_DATA    *dev_data,     /* (I) Dev data structure    */
                    TREE_DATA   *tree_data)    /* (I) Tree data structure   */
{

    double  *tmpZero;            /* Temporary pointer  */
    long    tmpMaturiy;          /* Temporary maturity */
    int     i;                   /*                    */
    int     status = FAILURE;    /* Error status       */

        
    if (Reset)
    {   
        /* If we have more than one zero, we permute them */
        if (*NbZero > 0)
        {
            tmpZero    = Zero[*NbZero-1];
            tmpMaturiy = ZeroMaturity[*NbZero-1];

            for (i = *NbZero - 1; i > 0; i--)
            {
                Zero[i] = Zero[i-1];
                ZeroMaturity[i] = ZeroMaturity[i-1];
            }

            /* If we don't have all the zeros in the bank we add one */
            if (*NbZero < TotNbZero)
            {
                Zero[0]               = Zero[*NbZero];
                Zero[*NbZero]         = tmpZero;
                ZeroMaturity[*NbZero] = tmpMaturiy;
                (*NbZero)++;
            }
            /* Otherwise we drop the last one */
            else
            {
                Zero[0] = tmpZero;
            }
        }
        else
        {
            (*NbZero)++;
                        
        }  /* if then else */           

        if (Zero_t (Zero[0],
                    TRUE,
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
        }

        ZeroMaturity[0] = CurrentDate;
    } 
    /* Discount the existing zeros */
    else
    {
        for (i = 0; i < *NbZero; i++)                                           
        {                                                                       
            if (Zero_t (Zero[i],                                        
                        0,                                                      
                        t,
                        T,
                        DCurve,
                        dev_data,
                        tree_data) == FAILURE)
            {
                goto RETURN;
            }
        }  /* for i */                        
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Zero_Bank */


/*****  ZeroPrice  ****************************************************/
/*
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

double   ZeroPrice(long     MatDate,   /* (I) Mat date of the zero       */
                   long     ValueDate, /* (I) Value date                 */
                   int      NbZero,    /* (I) Number of zeros in curve   */
                   long    *ZeroDates, /* (I) maturity dates in zero crv */
                   double  *ZeroRates) /* (I) Zero rates                 */
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

    if ((idx == 0) || (NbZero == 1)) /* MatDate <= 1st ZeroDate or crv has only 1 pt */
    {
        t1  = 0;
        ZR1 = ZeroRates[0];

        t2  = Daysact(ValueDate, ZeroDates[0])/365.0;
        ZR2 = ZeroRates[0];
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

