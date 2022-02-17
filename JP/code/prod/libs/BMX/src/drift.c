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
#include "bmx123head.h"

/*****  Drift_1D   **********************************************************/
/*
*       Calculate state prices for the one-factor model
*/
int     Drift_1D (  MKTVOL_DATA   *mktvol_data,   /* (I) Volatility data */
                    TREE_DATA     *tree_data)     /* (I/O) Tree data     */
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
    double  *XValuesL;

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
    int     j,k;               /* Critical date type and index              */
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
    tree_data->XValues    = (double **)DR_Array(DOUBLE_PTR,0,NbEDevDates-1);
    tree_data->YValues    = (double **)DR_Array(DOUBLE_PTR,0,NbEDevDates-1);
    

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
        tree_data->XValues[i]     = NULL;
        tree_data->YValues[i]     = NULL;
    }

    /* Finally allocate the slices themselves in preparation */
    for (i=0; i<NbEDevDates; i++)
    {
        tree_data->EDevStPrice[i]= Alloc_Slice(tree_data);
        tree_data->XValues[i]    = Alloc_Slice(tree_data);
        tree_data->YValues[i]    = Alloc_Slice(tree_data);      
        
        if (tree_data->EDevStPrice[i] == NULL    ||
            tree_data->XValues[i]     == NULL    ||
            tree_data->YValues[i]     == NULL    )
        {
            goto RETURN;
        }
    }


    StatePr  = Alloc_Slice (tree_data);
    StatePr1 = Alloc_Slice (tree_data);

    if (   (StatePr  == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Drift_1D: could not allocate memory!");
        goto RETURN;        
    }


    /* Develop state prices forward and store on Nmr dates into tree_data */
    offset = Node_Offset(1, 0, 0, 0, tree_data);
    StatePrL = StatePr + offset;
    StatePrL[0] = 1.;

    EDevPriceL    = tree_data->EDevStPrice[0] + offset;
    EDevPriceL[0] = 1.;    
    EDevIdx++;

    for (t = 0; t < NbTP; t++)
    {   
        StatePr1L = StatePr1 + Node_Offset(1, 0, 0, t+1, tree_data);
        StatePrL  = StatePr  + Node_Offset(1, 0, 0, t,   tree_data);
   
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
                    offset     = Node_Offset(1, 0, 0, t+1, tree_data);
                    EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;
                    XValuesL   = tree_data->XValues[EDevIdx] + offset;

                    for (i = Bottom[t+1]; i <= Top[t+1]; i++)
                    {
                        EDevPriceL[i] = StatePr1L[i];
                        XValuesL[i]    = i * Jump;
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

    Free_Slice (StatePr,  tree_data);
    Free_Slice (StatePr1, tree_data);
        
    return (status);

}  /* Drift_1D */




/*****  Drift_2D   **********************************************************/
/*
*       Calculate the interest rate drift for a two factor model.
*/
int     Drift_2D (MKTVOL_DATA     *mktvol_data,   /* (I)   mktvol data   */  
                  TREE_DATA       *tree_data)     /* (I/O) Tree data     */
{

     /* Local values of tree data structure */
    double  *Beta        = mktvol_data->Beta;
    double  **Aweight    = tree_data->Aweight;
    int     *Top1        = tree_data->Top1;
    int     *Bottom1     = tree_data->Bottom1;
    int     **Top2       = tree_data->Top2;
    int     **Bottom2    = tree_data->Bottom2;
    int     *OutTop1     = tree_data->OutTop1;
    int     *OutBottom1  = tree_data->OutBottom1;
    int     **OutTop2    = tree_data->OutTop2;
    int     **OutBottom2 = tree_data->OutBottom2;
    int     NbTP         = tree_data->NbTP;
    double  *Length      = tree_data->Length;
    double  *LengthJ     = tree_data->LengthJ;

    /* Other variables */
    double  *StatePr;            /* State Prices at t                         */
    double  *StatePr1;           /* State Prices at t+1                       */

                                                                              
    double  *StatePrL;           /* Local slice pointers                      */                                                      
    double  *StatePr1L;                                                       
    double  *StatePr2L;                                                       
    double  *StatePr3L;  
    double  *EDevPriceL;         /* St prices kept for express DEV'ing        */
    double  *XValuesL;           /* Local Pointer to x values                 */
    double  *YValuesL;           /* Local pointer to y values                 */



    int     NbEDevDates  = 0;
    long    Idx;                 /* Local variable for tree_data->NbPts       */ 
   
     
    
    double  *TempPointer;                                                 
                                                                              
    double  Jump[3];             /* Jump sizes at current time point          */
    double  PreviousJump[3];     /* Jump sizes at previous time point         */
    double  RateJumpLeft2;       /* Rate space jump size for last dimension   */
    double  RateJumpRight2;      /* Rate space jump size for last dimension   */
    double  JumpCoeff, d[3], du; /* Jump coefficients                         */
    double  Grid;                /* Grid points                               */
                                                                              
    double  Beta1, Beta2, Beta3; /* Modified mean reversion coefficients      */
    double  Pi, Qi, Qij;         /* Total drift                               */
    double  p;                   /* Total probability                         */
    double  x;                   /*                                           */
    
    double  pu, pd, p0;          /* Probabilities in first dimension          */
    double  qu, qd, q0;          /* Probabilities in second dimension         */

    int     EDevIdx = 0;         /* Counter of dates for express DEV          */
                                                                              

    int     offset;
    int     Mid;                 /* Mid of distribution index                 */
    int     i,k, j, j1, j2, j3;  /* Node indices                              */
    int     l, lMin, lMax;       /* Node branching shifts                     */ 
    int     m, mMin, mMax;                                                    
    int     t;                   /* Time step index                           */
    int     status = FAILURE;    /* Error status = FAILURE initially          */


    

    StatePr  = Alloc_Slice (tree_data);
    StatePr1 = Alloc_Slice (tree_data);

    if (   (StatePr  == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Drift_1D: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }


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
    tree_data->XValues    = (double **)DR_Array(DOUBLE_PTR,0,NbEDevDates-1);
    tree_data->YValues    = (double **)DR_Array(DOUBLE_PTR,0,NbEDevDates-1);
    
    
    tree_data->EDevDate   = (long *)   DR_Array(LONG,      0,NbEDevDates-1);    
     
    if ((tree_data->EDevStPrice== NULL) ||
        (tree_data->EDevDate   == NULL))
    {
        goto FREE_MEM_AND_RETURN;
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
        goto FREE_MEM_AND_RETURN;
    }

    /* Initialise pointers to NULL for safe freeing */
    for (i=0; i<NbEDevDates; i++)
    {
        tree_data->EDevStPrice[i] = NULL;
        tree_data->XValues[i]     = NULL;
        tree_data->YValues[i]     = NULL;
    }

    /* Finally allocate the slices themselves in preparation */
    for (i=0; i<NbEDevDates; i++)
    {
        tree_data->EDevStPrice[i]= Alloc_Slice(tree_data);
        tree_data->XValues[i]    = Alloc_Slice(tree_data);
        tree_data->YValues[i]    = Alloc_Slice(tree_data);

        if (tree_data->EDevStPrice[i] == NULL         ||
            tree_data->XValues[i]     == NULL         ||
            tree_data->YValues[i]     == NULL         )
        {
            goto FREE_MEM_AND_RETURN;
        }

    }


    /* Develop state prices forward and store on Nmr dates into tree_data */
    offset = Node_Offset(2, 0, 0, 0, tree_data);
    StatePrL = StatePr + offset;
    StatePrL[0] = 1.;

    EDevPriceL    = tree_data->EDevStPrice[0] + offset;
    EDevPriceL[0] = 1.;    
    EDevIdx++;



    for (t = 0; t < NbTP; t++)
    { 

        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
        PreviousJump[0] = Aweight[0][t-1] * du;
        PreviousJump[1] = Aweight[1][t-1] * du;
        PreviousJump[2] = Aweight[2][t-1] * du;
        
        du = sqrt (JUMPCOEFF * LengthJ[t]);
            
        Jump[0] = Aweight[0][t] * du;
        Jump[1] = Aweight[1][t] * du;
        Jump[2] = Aweight[2][t] * du;
        
        d[0] = (PreviousJump[0] - Jump[0]) / Jump[0];
        d[1] = (PreviousJump[1] - Jump[1]) / Jump[2];
        d[2] = (PreviousJump[2] - Jump[2]) / Jump[2];
            
        Beta1  = Beta[0] * Length[t] * PreviousJump[0] / Jump[0];
        Beta2  = Beta[1] * Length[t] * PreviousJump[1] / Jump[2];
        Beta3  = Beta[1] * Length[t] * PreviousJump[2] / Jump[2];
        
        JumpCoeff = Length[t]/(LengthJ[t]*JUMPCOEFF);

        

        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + Node_Offset(2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr1L[j] = 0.;
            }
        }  /* for i */
                
                        
        lMax = OutTop1[t+1] - 1;
        lMin = OutBottom1[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Drift_2D: problem in building the tree (lMin > lMax)!");
            goto FREE_MEM_AND_RETURN;
        }


        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {                        
            Pi = (d[0] - Beta1) * i;                
            Qi = (d[1] - Beta2) * i - Pi * Jump[1] / Jump[2];

            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                
            Pi -= l;
                                
            pu = .5 * (JumpCoeff + Pi + Pi * Pi);
            pd = pu - Pi;
            p0 = 1. - pu - pd;
                                
                                
            mMax =            OutTop2[t+1][i+l-1];
            mMax = MIN (mMax, OutTop2[t+1][i+l  ]);
            mMax = MIN (mMax, OutTop2[t+1][i+l+1]) - 1;

            mMin =            OutBottom2[t+1][i+l-1];
            mMin = MAX (mMin, OutBottom2[t+1][i+l  ]);
            mMin = MAX (mMin, OutBottom2[t+1][i+l+1]) + 1;

            if (mMin > mMax)
            {
                DR_Error ("Drift_2D: problem in building the tree (mMin > mMax)!");
                goto FREE_MEM_AND_RETURN;
            }


            StatePrL  = StatePr  + Node_Offset (2, i, 0, t, tree_data);
            
                                        
            StatePr1L = StatePr1 + Node_Offset (2, i+1+l, 0, t+1, tree_data);
            StatePr2L = StatePr1 + Node_Offset (2, i  +l, 0, t+1, tree_data);
            StatePr3L = StatePr1 + Node_Offset (2, i-1+l, 0, t+1, tree_data);

            
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                Qij = Qi + (d[2] - Beta3) * j;

                m = NEAR_INT (Qij);
                m = MIN (MAX (mMin - j, m), mMax - j);
                
                Qij -= m;
                                
                qu = .5 * (JumpCoeff + Qij + Qij * Qij);
                qd = qu - Qij;
                q0 = 1. - qu - qd;

                        
                j1 = j + 1 + m;                 
                j2 = j     + m;
                j3 = j - 1 + m;

                x = StatePrL[j] ;

                StatePr1L[j1] += pu * qu * x;
                StatePr1L[j2] += pu * q0 * x;
                StatePr1L[j3] += pu * qd * x;
                StatePr2L[j1] += p0 * qu * x;
                StatePr2L[j2] += p0 * q0 * x;
                StatePr2L[j3] += p0 * qd * x;
                StatePr3L[j1] += pd * qu * x;
                StatePr3L[j2] += pd * q0 * x;
                StatePr3L[j3] += pd * qd * x;
            }
        }  /* for i */

                           
        p = 0.;

        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + Node_Offset(2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {
                p += StatePr1L[j];
            }
        }  /* for i */

        if (fabs (p - 1.0) > 0.0001)
        {
            DR_Error ("Drift_2D: state prices don't add up to 1: "
                        "increase sigma or Ppy!");
            goto FREE_MEM_AND_RETURN;
        }


        /* Store ALL state prices */
        if (EDevIdx < NbEDevDates)
        {
            for (j = 0; j < NBCRITDATE; j++)
            {
                if (tree_data->TPtype[j][t+1])
                {
                    Idx = 0;
                    /* Store corresponding state prices */
                    for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
                    {
                        StatePr1L  = StatePr1 + Node_Offset(2, i, 0, t+1, tree_data);
                        EDevPriceL = tree_data->EDevStPrice[EDevIdx] + Node_Offset(2, i, 0, t+1, tree_data);            
                        XValuesL   = tree_data->XValues[EDevIdx]     + Node_Offset(2, i, 0, t+1, tree_data);            
                        YValuesL   = tree_data->YValues[EDevIdx]     + Node_Offset(2, i, 0, t+1, tree_data);            
                        
                        for (k = Bottom2[t+1][i]; k <= Top2[t+1][i]; k++)
                        {
                            EDevPriceL[k] = StatePr1L[k];
                            XValuesL[k]    = i * Jump[0];
                            YValuesL[k]    = i * Jump[1] + k * Jump[2];                         
                            
                            Idx ++;
                        }
                    }
                    EDevIdx ++;                 
                    break;
                }
            }
        }           
        /* Permute variables in preparation for next time step */          
        TempPointer = StatePr;
        StatePr     = StatePr1;
        StatePr1    = TempPointer;


    }  /* for t */ 
    
    
    /* Final check on number of Nmr dates processed */
    if (EDevIdx != NbEDevDates)
    {
        DR_Error("NmrTool:  Error in processing  the NmrTool dates.\n"
                 "Please ensure that dates are ordered, that there\n"
                 "is no repetition, and that the NmrTool dates\n"
                 "are also critical dates.\n");
        goto FREE_MEM_AND_RETURN;
    }
    
                                        

                                                
    status = SUCCESS;

FREE_MEM_AND_RETURN:

    Free_Slice (StatePr,  tree_data);
    Free_Slice (StatePr1, tree_data);
        
    return (status);


}  /* Drift_2D */


/*****  Drift   **********************************************************/
/*
*       Calculate the interest rate drift.
*/
int     Drift    (MKTVOL_DATA     *mktvol_data,   /* (I)   mktvol data   */  
                  TREE_DATA       *tree_data)     /* (I/O) Tree data     */
{
    int status = FAILURE;

    /* 1 Factor tree */
    if (tree_data->NbFactor == 1)
    {
        if (Drift_1D(mktvol_data,
                     tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }
    
    /* 2 Factor tree */
    else if (tree_data->NbFactor == 2)
    {
        if (Drift_2D(mktvol_data,
                     tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }


    status = SUCCESS;

RETURN:

    return(status);
}

































