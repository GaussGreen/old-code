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
#include "fix123head.h"


/*****  Fix3_Zero_t  *************************************************************/
/*
*       Zero coupon
*/
int     Fix3_Zero_t (double      *Zero,         /* (I/O) Zero prices       */
                long        Reset,         /* (I) Reset flag          */
                int         t,             /* (I) Current time point  */
                int         T,             /* (I) Last time point     */
                int         DCurve,        /* (I) Discount curve      */
                FIX3_DEV_DATA    *dev_data,     /* (I) Fix3_Dev data structure  */
                FIX3_TREE_DATA   *tree_data)    /* (I) Tree data structure */
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
            ZeroL = Zero + Fix3_Node_Offset(1, 0, 0, t, tree_data);
    
            for (i = Bottom1; i <= Top1; i ++)                                      
            {
                ZeroL[i] = 1.;
            }
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                ZeroL = Zero + Fix3_Node_Offset(2, i, 0, t, tree_data);
    
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
                    ZeroL = Zero + Fix3_Node_Offset(3, i, j, t, tree_data);
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        ZeroL[k] = 1.;
                    }
                }  /* for j */	
        }  /* if then else */
    }
    else
    {
        if (Fix3_Dev (   Zero,
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

}  /* Fix3_Zero_t */



/*****  Fix3_Zero_Bank  **********************************************************/
/*
*       Manipulation of the zero bank.
*/
int     Fix3_Zero_Bank ( double      **Zero,        /* (I/O) Set of zeros        */
                    long        *ZeroMaturity, /* (I/O) Zero maturities     */
                    int         *NbZero,       /* (I/O) Current nb of zeros */
                    int         TotNbZero,     /* (I) Total number of zeros */
                    long        Reset,         /* (I) Reset flag            */
                    long        CurrentDate,   /* (I) Current date          */
                    int         t,             /* (I) Current time point    */
                    int         T,             /* (I) Last time point       */
                    int         DCurve,        /* (I) Discount curve        */
                    FIX3_DEV_DATA    *dev_data,     /* (I) Fix3_Dev data structure    */
                    FIX3_TREE_DATA   *tree_data)    /* (I) Tree data structure   */
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

        if (Fix3_Zero_t (Zero[0],
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
            if (Fix3_Zero_t (Zero[i],                                        
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

}  /* Fix3_Zero_Bank */


