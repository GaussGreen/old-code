/****************************************************************************/
/*      Building the interest rate tree: calculation of forward rates and   */
/*      interest rate drift at each time step.                              */
/****************************************************************************/
/*      DRIFT.c                                                             */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"

/*****  Drift_1D   **********************************************************/
/*
*       Calculate state prices for the one-factor model
*/
int     Fix3_Drift_1D_Tmx (  
                    MKTVOL_DATA    *mktvol_data,   /* (I) Volatility data */
                    FIX3_TREE_DATA *tree_data)     /* (I/O) Tree data     */
{

    /* Local values of tree data structure */
    double  Beta         = mktvol_data->Beta[0];
    double  **Aweight    = tree_data->Aweight;
    int     *Top         = tree_data->Top1;
    int     *Bottom      = tree_data->Bottom1;
    int     *OutTop      = tree_data->OutTop1;
    int     *OutBottom   = tree_data->OutBottom1;
    int     NbTP         = tree_data->NbTP;
    double  *Length      = tree_data->Length;
    double  *LengthJ     = tree_data->LengthJ;
    int     NbEDevDates  = 0;

    /* Other variables */
    double  *StatePr=NULL;     /* State Prices at t                         */
    double  *StatePr1=NULL;    /* State Prices at t+1                       */

    double  *StatePrL;         /* Local slice pointers                      */
    double  *StatePr1L;     
    double  *EDevPriceL;       /* St prices kept for express DEV'ing        */

    double  Jump;              /* Size of the jump (in log space)           */
    double  PreviousJump;      /* Jump size at previous period              */
    double  JumpCoeff, du;     /* Jump coefficients                         */

    double  BetaLocal;         /* = Beta * Length[t]                        */
    double  Pi;                /* Drift at the current node                 */
    double  d;                 /* Drift due to the change in jump size      */
    double  p;                 /* Total probability                         */
    double  x;                 /*                                           */

    double  pu, pd, p0;        /* Probabilities                             */

    int     EDevIdx=0;         /* Counter of dates for express DEV          */
    int     offset;
    int     i;                 /* Node index                                */
    int     j;                 /* Critical date type and index              */
    int     l, lMin, lMax;     /* Node branching shift                      */
    int     t;                 /* Time step index                           */
    int     status = FAILURE;  /* Error status = FAILURE initially          */


    /* Count critical dates */
    for (t = 0; t <=tree_data->NbTP; t++)
    {
        for (j = 0; j < NBCRITDATE; j++)
        {
            if (tree_data->TPtype[j][t])
            {
                NbEDevDates ++;
                break;
            }
        }
    }
    tree_data->NbEDevDates = NbEDevDates;

    /* Allocate array of pointers to state price slices */
    tree_data->EDevStPrice= (double **)DR_Array(DOUBLE_PTR,0,NbEDevDates-1);  
    tree_data->EDevDate   = (long *)   DR_Array(LONG,      0,NbEDevDates-1);
 
    if ((tree_data->EDevStPrice== NULL) ||
        (tree_data->EDevDate   == NULL))
    {
        goto RETURN;
    }

    /* Store critical dates */
    i = 0;
    for (t = 0; t <= tree_data->NbTP; t++)
    {
        for (j = 0; j < NBCRITDATE; j++)
        {
            if (tree_data->TPtype[j][t])
            {
                tree_data->EDevDate[i] = tree_data->TPDate[t];
                i++;
                break;
            }
        }
    }

    if (i != NbEDevDates)
    { 
        DR_Error ("Not all express DEV dates were processed");
        goto RETURN;
    }

    /* Initialise pointers to NULL for safe freeing */
    for (i=0; i<NbEDevDates; i++)
    {
        tree_data->EDevStPrice[i] = NULL;
    }

    /* Finally allocate the slices themselves in preparation */
    for (i=0; i<NbEDevDates; i++)
    {
        tree_data->EDevStPrice[i]= Fix3_Alloc_Slice(tree_data);
        if (tree_data->EDevStPrice[i] == NULL)
        {
            goto RETURN;
        }
    }


    StatePr  = Fix3_Alloc_Slice (tree_data);
    StatePr1 = Fix3_Alloc_Slice (tree_data);

    if (   (StatePr  == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Drift_1D: could not allocate memory!");
        goto RETURN;        
    }


    /* Develop state prices forward and store on Nmr dates into tree_data */
    offset = Fix3_Node_Offset(1, 0, 0, 0, tree_data);
    StatePrL = StatePr + offset;
    StatePrL[0] = 1.;

    EDevPriceL    = tree_data->EDevStPrice[0] + offset;
    EDevPriceL[0] = 1.;    
    EDevIdx++;

    for (t = 0; t < NbTP; t++)
    {   
        StatePr1L = StatePr1 + Fix3_Node_Offset(1, 0, 0, t+1, tree_data);
        StatePrL  = StatePr  + Fix3_Node_Offset(1, 0, 0, t,   tree_data);
   
        /* 
         *  Precompute jumps and probabilities coefficients.
         */
        du           = sqrt (JUMPCOEFF * LengthJ[t-1]);
        PreviousJump = Aweight[0][t-1] * du;
                                                                            
        du   = sqrt (JUMPCOEFF * LengthJ[t]);
        Jump = Aweight[0][t] * du;
        d    = PreviousJump / Jump - 1.;
    
        BetaLocal = Beta * Length[t] * PreviousJump / Jump;
        JumpCoeff = Length[t]/(LengthJ[t]*JUMPCOEFF);

        /* 
         *   Compute state prices for next period.
         */
        for (i = Bottom[t+1]; i <= Top[t+1]; i++)
        {
            StatePr1L[i] = 0.;                        
        }
            

        /* 
         * Branching has to be inside outer ellipse at next period 
         */
        lMax = OutTop[t+1] - 1;                                             
        lMin = OutBottom[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Drift_1D: problem in building the tree "
                      "(lMin > lMax)!");
            goto RETURN;
        }


        for (i = Bottom[t]; i <= Top[t]; i++)
        {                        
            Pi = (d - BetaLocal) * i;
                                                                            
            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                                                                            
            Pi -= l;
                    
                    
            pu = .5 * (JumpCoeff + Pi + Pi * Pi);
            pd = pu - Pi;
            p0 = 1. - pu - pd;

            x = StatePrL[i];

            StatePr1L[i+1+l] += pu * x;
            StatePr1L[i  +l] += p0 * x;
            StatePr1L[i-1+l] += pd * x;
        }


        /* 
         *   Check that state prices add up to 1.
         */

        p = 0.;

        for (i = Bottom[t+1]; i <= Top[t+1]; i++)
        {
            p += StatePr1L[i];
        }

        if (fabs (p - 1.) > 0.0001)
        {
            DR_Error ("Drift_1D: state prices don't add up to 1: "
                      "increase sigma or Ppy!");
            goto RETURN;
        }

        /* Store ALL state prices */
        if (EDevIdx < NbEDevDates)
        {
            for (j = 0; j < NBCRITDATE; j++)
            {
                if (tree_data->TPtype[j][t+1])
                {
                    offset     = Fix3_Node_Offset(1, 0, 0, t+1, tree_data);
                    EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;

                    for (i = Bottom[t+1]; i <= Top[t+1]; i++)
                    {
                        EDevPriceL[i] = StatePr1L[i];
                    }
                    EDevIdx ++;
                    break;
                }
            }
        }

        /* 
         *   Permute values for the next iteration.
         */
        {
            double *TempPointer = StatePr;

            StatePr     = StatePr1;
            StatePr1    = TempPointer;
        }

    }  /* for t */      
    
    /* Final check on number of Nmr dates processed */
    if (EDevIdx != NbEDevDates)
    {
        DR_Error("NmrTool:  Error in processing  the NmrTool dates.\n"
                 "Please ensure that dates are ordered, that there\n"
                 "is no repetition, and that the NmrTool dates\n"
                 "are also critical dates.\n");
        goto RETURN;
    }

                                                
    status = SUCCESS;

    RETURN:

    Fix3_Free_Slice (StatePr,  tree_data);
    Fix3_Free_Slice (StatePr1, tree_data);
        
    return (status);

}  /* Fix3_Drift_1D_Tmx */

