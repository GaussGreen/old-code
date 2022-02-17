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



/*****  Fix3_Drift_2D_Smd   *******************************************************/
/*
*       Calculate the interest rate drift for a two factor model.
*/
int     Fix3_Drift_2D_Smd (MKTVOL_DATA     *mktvol_data,   /* (I) Volatility data */
                           FIX3_TREE_DATA  *tree_data)     /* (I/O) Tree data     */
{

    /* Local values of tree data structure */
    /* Use diffused zero and fwd curves    */
    double  *ZCenter     = tree_data->ZCenter;
    double  *ZeroCoupon  = tree_data->ZeroCoupon[tree_data->CvDiff];
    double  *FwdRate     = tree_data->FwdRate[tree_data->CvDiff];
    double  QLeft        = mktvol_data->QLeft;
    double  QRight       = mktvol_data->QRight;
    double  FwdShift     = mktvol_data->FwdShift;
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
    double  *Discount;           /* 1 period discount factor                  */
    double  *StatePr1;           /* State Prices at t+1                       */

                                                                              
    double  *StatePrL;           /* Local slice pointers                      */
    double  *DiscountL;                                                       
    double  *StatePr1L;                                                       
    double  *StatePr2L;                                                       
    double  *StatePr3L;  
    double  *EDevPriceL;         /* St prices kept for express DEV'ing        */
   
     
    
    double  *TempPointer;                                                 
                                                                              
    double  Jump[3];             /* Jump sizes at current time point          */
    double  PreviousJump[3];     /* Jump sizes at previous time point         */
    double  RateJumpLeft2;       /* Rate space jump size for last dimension   */
    double  RateJumpRight2;      /* Rate space jump size for last dimension   */
    double  JumpCoeff, d[3], du; /* Jump coefficients                         */
    double  Grid;                /* Grid points                               */
    double  MeanDrift;
    double  MDJump;
                                                                              
    double  Beta1, Beta2, Beta3; /* Modified mean reversion coefficients      */
    double  Pi, Qi, Qij;         /* Total drift                               */
    double  p;                   /* Total probability                         */
    double  Zt=0;                /* Center of the tree in X-space             */
    double  DelZt;               /* Correction to Zt for current time point   */
    double  Zidx;                /* Zt index i,j,k adjusted                   */
    double  P[3];                /* 2nd degree polynomial coefficients        */
    double  x;                   /*                                           */
    double  FwdRateA;            /* Fwd rate adjusted                         */
    double  MLeft, SLeft;        /* Multiple coeff and shift for grid point   */
    double  MRight, SRight;
    double  VolBbq;              /* Sigma used in bone mapping                */
    double  Bbq;                 /* Backbone parameter                        */
    double  QMid;                /* Average q parameter                       */
    double  QSh;                 /* q parameter for shift adjustment          */
    double  InsRate;             /* Instantaneous rate                        */
    double  Q, U;                /* Q and U weight coeff for Z-polynomial     */ 
    double  D0;                  /* Consequtive power of discount factor      */ 
    double  aL00;                /* Coeff of Z-polynomial across all states   */
    double  aL10, aL11 ;
    double  aL20, aL21, aL22 ;
    double  aR00;               
    double  aR10, aR11 ;
    double  aR20, aR21, aR22 ;


    double  pu, pd, p0;          /* Probabilities in first dimension          */
    double  qu, qd, q0;          /* Probabilities in second dimension         */

    int     EDevOn = FALSE;      /* To indicate whether EDev tool is used     */
    int     EDevIdx = 0;         /* Counter of dates for express DEV          */
                                                                              
    int     offset;
    int     Mid;                 /* Mid of distribution index                 */
    int     i, j, j1, j2, j3;    /* Node indices                              */
    int     l, lMin, lMax;       /* Node branching shifts                     */ 
    int     m, mMin, mMax;                                                    
    int     t;                   /* Time step index                           */
    int     status = FAILURE;    /* Error status = FAILURE initially          */


    StatePr  = Fix3_Alloc_Slice (tree_data);
    Discount = Fix3_Alloc_Slice (tree_data);
    StatePr1 = Fix3_Alloc_Slice (tree_data);

    if (   (StatePr  == NULL)
        || (Discount == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Fix3_Drift_2D: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }


    /* Check if Express DEV tool is being used */
    if (tree_data->NbEDevDates > 0)
    {

        EDevOn = TRUE;
        EDevIdx = 0;

        /* Allocate array of pointers to state price slices */
        tree_data->EDevStPrice = (double **)DR_Array(DOUBLE_PTR, 
                                                     0, 
                                                     tree_data->NbEDevDates-1);
        if (tree_data->EDevStPrice == NULL)
        {
            goto FREE_MEM_AND_RETURN;
        }


        /* Initialise pointers to NULL for safe freeing */
        for (i=0; i<tree_data->NbEDevDates; i++)
        {
            tree_data->EDevStPrice[i] = NULL;
        }

        /* Finally allocate the slices themselves in preparation */
        for (i=0; i<tree_data->NbEDevDates; i++)
        {
            tree_data->EDevStPrice[i] = Fix3_Alloc_Slice(tree_data);
            if (tree_data->EDevStPrice[i] == NULL)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }

    }

    Bbq  = mktvol_data->Bbq;
    QMid = (QLeft + QRight) / 2;

    StatePrL = StatePr + Fix3_Node_Offset(2, 0, 0, 0, tree_data);
    StatePrL[0] = 1.;

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
        MDJump = Beta[1] * Length[t]                   / Jump[2];
        
        JumpCoeff = Length[t]/(LengthJ[t]*JUMPCOEFF);

        FwdRateA  = FwdRate[t] / (1. + FwdShift);
        InsRate   = FwdRate[t] / Length[t];    

        VolBbq  = mktvol_data->VolLogn * Bbq  
                + mktvol_data->VolNorm *(1. - Bbq) / InsRate;
        VolBbq *= (1. + FwdShift) / (1. + QMid * FwdShift);

        if (IS_Q(QLeft))
        {
            MLeft         = FwdRateA / QLeft;
            SLeft         = 1 + FwdRateA - FwdRateA / QLeft;      
            RateJumpLeft2 = exp(QLeft * VolBbq * PreviousJump[2]);
        }
        else
        {
            MLeft         = FwdRateA;
            SLeft         = 1 + FwdRateA;
            RateJumpLeft2 = FwdRateA * VolBbq * PreviousJump[2];
        }
        if (IS_Q(QRight))
        {
            MRight         = FwdRateA / QRight;
            SRight         = 1 + FwdRateA - FwdRateA / QRight;      
            RateJumpRight2 = exp(QRight * VolBbq * PreviousJump[2]);
        }
        else
        {
            MRight         = FwdRateA;
            SRight         = 1 + FwdRateA;
            RateJumpRight2 = FwdRateA * VolBbq * PreviousJump[2];
        }

        if (t == 0)
        {
            QSh = (FwdShift > 0) ? QRight : QLeft;
            if (IS_Q(QSh))
            {
                Zt = log(1. + QSh * FwdShift) / (QSh * VolBbq);
            }
            else
            {
                Zt = FwdShift / VolBbq;
            }
        }

        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            /* LEFT part of distribution */

            Zidx  = Zt 
                  + (PreviousJump[1]) * i;
            Mid   = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);
            Zidx += (PreviousJump[2]) * Bottom2[t][i];

            DiscountL = Discount + Fix3_Node_Offset (2, i, 0, t, tree_data);
    
            if (IS_Q(QLeft))
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid *= RateJumpLeft2;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid += RateJumpLeft2;
                }
            }

            /* Right part of distribution */

            Zidx = Zt 
                 + (PreviousJump[1]) * i
                 + (PreviousJump[2]) * (Mid + 1);

            if (IS_Q(QRight))
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid *= RateJumpRight2;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid += RateJumpRight2;
                }
            }

        }  /* for i */


                
        for (i = 0; i < 3; i++)                     
        {      
            P[i] = 0.; 
        }

        Q = QLeft * VolBbq;
        U = (FwdRateA * (1. - QLeft) - QLeft) * VolBbq; 

        aL00 =  1.;
        aL10 = -Q;          aL11 = -U;
        aL20 =  0.5*Q*Q;    aL21 =  1.5*U*Q;    aL22 =   U*U;

        Q = QRight * VolBbq;
        U = (FwdRateA * (1. - QRight) - QRight) * VolBbq; 

        aR00 =  1.;
        aR10 = -Q;          aR11 = -U;
        aR20 =  0.5*Q*Q;    aR21 =  1.5*U*Q;    aR22 =   U*U;

        
        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {    
            StatePrL  = StatePr  + Fix3_Node_Offset (2, i, 0, t, tree_data);
            DiscountL = Discount + Fix3_Node_Offset (2, i, 0, t, tree_data);
                                                
            Zidx = Zt 
                 + (PreviousJump[1]) * i;
            Mid  = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid  = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);

            for (j = Bottom2[t][i]; j <= Mid ; j++)
            {
                
                x   = DiscountL[j];
                D0  = StatePrL[j] * x;
                P[0] += D0;

                P[1] += (aL10 + aL11 * x) * D0;
            
                P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  

            }

            for (j = Mid + 1; j <= Top2[t][i] ; j++)
            {
                
                x   = DiscountL[j];
                D0  = StatePrL[j] * x;
                P[0] += D0;

                P[1] += (aR10 + aR11 * x) * D0;
            
                P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  

            }
        }  /* for i */
        

        P[0] -= ZeroCoupon[t+1];

        /* Find change in tree centre at current time step*/
        {
            double Root[2] = {0.0, 0.0};

            /* Find root of drift equation (2nd-order expansion in DelZt) */
            if (Quadratic_Solve(Root, P) == FAILURE)
            {
                DR_Error ("Fix3_Drift_2D: no solution exists to drift equation!");
                goto FREE_MEM_AND_RETURN;
            }

            /* Choose root nearest to 0.0 */
            DelZt = ( (fabs(Root[0]) < fabs(Root[1])) ? Root[0] : Root[1] );

            /* Check power series expansion was valid */ 
            if (fabs(DelZt * VolBbq) > 1.5)
            {
                DR_Error ("Fix3_Drift_2D: change in drift is too large");
                goto FREE_MEM_AND_RETURN;
            }
        }


        Zt += DelZt;
        ZCenter[t] = Zt;


        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            /* LEFT part of distribution */

            Zidx  = Zt 
                  + (PreviousJump[1]) * i;
            Mid   = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);
            Zidx += (PreviousJump[2]) * Bottom2[t][i];

            DiscountL = Discount + Fix3_Node_Offset (2, i, 0, t, tree_data);
    
            if (IS_Q(QLeft))
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid *= RateJumpLeft2;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid += RateJumpLeft2;
                }
            }

            /* Right part of distribution */

            Zidx = Zt 
                 + (PreviousJump[1]) * i
                 + (PreviousJump[2]) * (Mid + 1);

            if (IS_Q(QRight))
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid *= RateJumpRight2;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);

                    Grid += RateJumpRight2;
                }
            }

        }  /* for i */


        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + Fix3_Node_Offset(2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr1L[j] = 0.;
            }
        }  /* for i */
                
                        
        lMax = OutTop1[t+1] - 1;
        lMin = OutBottom1[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Fix3_Drift_2D: problem in building the tree (lMin > lMax)!");
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
                DR_Error ("Fix3_Drift_2D: problem in building the tree (mMin > mMax)!");
                goto FREE_MEM_AND_RETURN;
            }


            StatePrL  = StatePr  + Fix3_Node_Offset (2, i, 0, t, tree_data);
            DiscountL = Discount + Fix3_Node_Offset (2, i, 0, t, tree_data);
                                        
            StatePr1L = StatePr1 + Fix3_Node_Offset (2, i+1+l, 0, t+1, tree_data);
            StatePr2L = StatePr1 + Fix3_Node_Offset (2, i  +l, 0, t+1, tree_data);
            StatePr3L = StatePr1 + Fix3_Node_Offset (2, i-1+l, 0, t+1, tree_data);



            /* Mean drift at X1 node i */
            MeanDrift = Smd_MeanDrift(i * PreviousJump[0], 
                                      mktvol_data->Afac,
                                      mktvol_data->Bfac, 
                                      mktvol_data->Cfac);

            
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                Qij = Qi + (d[2] - Beta3) * j + MeanDrift * MDJump;

                m = NEAR_INT (Qij);
                m = MIN (MAX (mMin - j, m), mMax - j);
                
                Qij -= m;
                                
                qu = .5 * (JumpCoeff + Qij + Qij * Qij);
                qd = qu - Qij;
                q0 = 1. - qu - qd;

                        
                j1 = j + 1 + m;                 
                j2 = j     + m;
                j3 = j - 1 + m;

                x = StatePrL[j] * DiscountL[j];

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
            StatePr1L = StatePr1 + Fix3_Node_Offset(2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {
                p += StatePr1L[j];
            }
        }  /* for i */

        if (fabs (p - ZeroCoupon[t+1]) > 0.0001)
        {
            DR_Error ("Fix3_Drift_2D: state prices don't add up to 1: "
                        "increase sigma or Ppy!");
            goto FREE_MEM_AND_RETURN;
        }


        /* Express DEV tool: store state prices if necessary   */
        if ((EDevOn) && (EDevIdx<tree_data->NbEDevDates))
        {
            if (tree_data->TPDate[t+1] == tree_data->EDevDate[EDevIdx])
            {

                /* Store corresponding state prices */
                for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
                {
                    offset =  Fix3_Node_Offset(2, i, 0, t+1, tree_data);
                    StatePr1L = StatePr1 + offset;
                    EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;
            
                    for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
                    {
                        EDevPriceL[j] = StatePr1L[j];
                    }
                }  /* for i */

                /* Increment the EDev date index */
                EDevIdx++;
            }
        }
    


    
        
        /* Permute variables in preparation for next time step */          
        TempPointer = StatePr;
        StatePr     = StatePr1;
        StatePr1    = TempPointer;


    }  /* for t */ 
    
    
    /* Final check on number of EDev dates processed */
    if (EDevOn)
    {
        if (EDevIdx != tree_data->NbEDevDates)
        {
            DR_Error("Express DEV:  Error in processing  the DEV dates.\n"
                     "Please ensure that dates are ordered, that there\n"
                     "is no repetition,that value date is not included\n"
                     "and that the express DEV dates are also critical\n" 
                     "dates .\n");
            goto FREE_MEM_AND_RETURN;
        }
    }
    
                                        

                                                
    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Fix3_Free_Slice (StatePr,  tree_data);
    Fix3_Free_Slice (Discount, tree_data);
    Fix3_Free_Slice (StatePr1, tree_data);
        
    return (status);

}  /* Fix3_Drift_2D_Smd */

