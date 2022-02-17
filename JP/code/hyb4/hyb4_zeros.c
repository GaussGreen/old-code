/****************************************************************************/
/*      Manipulation of the zeros carried along in the interest rate tree.  */
/****************************************************************************/
/*      HYB4_ZEROS.c                                                        */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

/*****  Hyb4_Zero_t  **************************************************************/
/**
 *      'Updates' (i.e. DEV's or resets) a zero coupon bond being carried
 *      backwards in the tree for reference.
 *
 */
int Hyb4_Zero_t (TSLICE      Zero,       /* (I/O) Zero prices                 */
                long        Reset,      /* (I) Reset flag                    */
                int         t,          /* (I) Current time period           */
                int         T,          /* (I) Total number of period        */
                int         DCurve,     /* (I) Discount curve                */
                int         DMode,      /* (I) Dimension of DEV and of zeros */
                HYB4_DEV_DATA    *dev_data,  /* (I) Hyb4_Dev data structure            */
                HYB4_TREE_DATA   *tree_data) /* (I) Tree data structure           */
{


    double	*ZeroL;
    int	                                
              Top1,   Bottom1,                          /* Tree limits (1rst dim) */
             *Top2,  *Bottom2,                          /* Tree limits (2nd dim)  */
            **Top3, **Bottom3,                          /* Tree limits (3rd dim)  */
           ***Top4,***Bottom4,                          /* Tree limits (3rd dim)  */

            i, j, k,                                    /* Node indices           */
            status = FAILURE;                           /* Error status	          */


    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];


    if (Reset)
    {
        if (DMode == DISC_1D_NOCUPS)
        {
            ZeroL = (double *) Zero + tree_data->NodeOffset0[t];

            for (i = Bottom1; i <= Top1; i ++)                                      
            {
                ZeroL[i] = 1.;

            }  /* for i */
            
        }
        else if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
        {
            

            for (i = Bottom1; i <= Top1; i ++)
            {
                ZeroL = (double *) Zero + tree_data->NodeOffset1[t][i];
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    ZeroL[j] = 1.;

                }  /* for j */	
            }

        }
        else if (DMode == DISC_3D_CUPS)
        {

            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    ZeroL = (double *) Zero + tree_data->NodeOffset2[t][i][j];
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        ZeroL[k] = 1.;

                    }  /* for k */	
                }
            }
        }  /* if then else */

    }
    else
    {
        if (Hyb4_Dev(Zero,
                t,
                T,
                DCurve,
                DMode,
                dev_data,
                tree_data) == FAILURE)
        {
            goto RETURN;
                    
        }  /* if */
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb4_Zero_t */

/*****  Hyb4_FXZero_t  **************************************************************/
/**
 *      'Updates' (i.e. DEV's or resets) a (zero coupon bond * FX spot(maturity date))
 *      being carried  backwards in the tree for reference.
 *      This is to be able to use the zero/claim banks to maintain discount slices
 *      for foreign denomiated payments in their domestic equivalent 
 *      (ie. discount 1.0 * FXSpot)
 *      
 *
 */
int     Hyb4_FXZero_t (TSLICE      FXZero,       /* (I/O) FXZero prices      */
                long        Reset,      /* (I) Reset flag                    */
                int         t,          /* (I) Current time period           */
                int         T,          /* (I) Total number of period        */
                int         DCurve,     /* (I) Discount curve                */
                int         DMode,      /* (I) Dimension of DEV and of zeros */
                HYB4_DEV_DATA    *dev_data,  /* (I) Hyb4_Dev data structure            */
                HYB4_TREE_DATA   *tree_data) /* (I) Tree data structure           */
{

    int	                                
              Top1,   Bottom1,                          /* Tree limits (1rst dim) */
             *Top2,  *Bottom2,                          /* Tree limits (2nd dim)  */
            **Top3, **Bottom3,                          /* Tree limits (3rd dim)  */
            status = FAILURE;                           /* Error status	          */


    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];

    if (DMode != DISC_3D_CUPS)
    {
        DR_Error("Discount mode must be DISC_3D_CUPS to discount an FX "
                 "denominated zero bank");
        goto RETURN;
    }

    if (Reset)
    {
        /* set reset slice as current date's spot FX slice
           (instead of 1.0 as in traditional zero bank)*/
        if (Hyb4_CopySlice(FXZero,
                           dev_data->AssetValue2[0],
                           3,  /* FX slice always 3d */
                           t,
                           tree_data) == FAILURE)
        {
            DR_Error("Unable to copy FxSpot slice at time step %d\n", t);
            goto RETURN;
        }
    }
    else
    {
        if (Hyb4_Dev(FXZero,
                t,
                T,
                DCurve,
                DMode,
                dev_data,
                tree_data) == FAILURE)
        {
            goto RETURN;
                    
        }
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb4_FXZero_t */

/*****  Hyb4_Zero_Bank  ***********************************************************/
/*
*       Manipulation of the zero bank.
*/
int   Hyb4_Zero_Bank(TSLICE      *Zero,            /* (I/O) Set of zeros          */
                long        *ZeroMaturity,    /* (I/O) Zero maturities       */
                int         *NbZero,          /* (I/O) Current nb of zeros   */
                int         TotNbZero,        /* (I) Total number of zeros   */
                long        Reset,            /* (I) Reset flag	             */
                long        CurrentDate,      /* (I) Current date            */
                int         t,                /* (I) Current time period     */
                int         T,                /* (I) Total number of period  */
                int         DCurve,           /* (I) Discount curve          */
                int         DMode,            /* (I) Dim of DEV and of zeros */
                HYB4_DEV_DATA    *dev_data,        /* (I) Hyb4_Dev data structure      */
                HYB4_TREE_DATA   *tree_data)       /* (I) Tree data structure     */
{
    double
            *tmpZero;                              /* Temporary pointer      */
    long 
            tmpMaturiy;                            /* Temporary maturity     */
    int	
            i,
            status = FAILURE;	                   /* Error status	         */

        
    if (Reset)
    {	                                
        if (*NbZero > 0)                                                        /* If we have more than one zero, we permute them */
        {
            tmpZero  = Zero[*NbZero-1];
            tmpMaturiy = ZeroMaturity[*NbZero-1];

            for (i = *NbZero - 1; i > 0; i--)                                   /* Shift the zeros */
            {
                Zero[i] = Zero[i-1];
                ZeroMaturity[i] = ZeroMaturity[i-1];

            }  /* for i */

            if (*NbZero < TotNbZero)                                            /* If we don't have all the zeros in the bank we add an extra zero */
            {
                Zero[0]               = Zero[*NbZero];
                Zero[*NbZero]         = tmpZero;
                ZeroMaturity[*NbZero] = tmpMaturiy;
                (*NbZero)++;
            }                                                                   /* Otherwise we drop the last one */
            else
            {
                Zero[0] = tmpZero;

            }  /* if then else */
        }
        else
        {
            (*NbZero)++;
                        
        }  /* if then else */        	

        if (Hyb4_Zero_t (Zero[0],                                                    /* Zero[0] is the zero coupon maturing now. We reset it to 1 */
                    TRUE,
                    t,
                    T,
                    DCurve,
                    DMode,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
                    
        }  /* if */

        ZeroMaturity[0] = CurrentDate;

    } 
    else                                                                        /* Discount the existing zeros */
    {        
        for (i = 0; i < *NbZero; i++)                                           
        {                                                                       
            if (Hyb4_Zero_t (Zero[i],                                        
                        0,                                                      
                        t,
                        T,
                        DCurve,
                        DMode,
                        dev_data,
                        tree_data) == FAILURE)
            {
                goto RETURN;
        
            }  /* if */                                                
        }  /* for i */                        
    }  /* if then else */


    status = SUCCESS;


    RETURN:

    return (status);

}  /* Hyb4_Zero_Bank */

