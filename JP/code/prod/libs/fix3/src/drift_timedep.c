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

/*****  Fix3_Drift_1D   **********************************************************/
/*
*       Calculate the interest rate drift for a one factor model.
*/
int     Fix3_Drift_1D_TimeDep (  
                    MKTVOL_DATA     *mktvol_data,        /* (I) Volatility data */
                    FIX3_TREE_DATA       *tree_data)     /* (I/O) Tree data     */
{

    /* Local values of tree data structure */
    /* Use diffused zero and fwd curves    */
    double  *ZCenter     = tree_data->ZCenter;
    double  *ZeroCoupon  = tree_data->ZeroCoupon[tree_data->CvDiff];
    double  *FwdRate     = tree_data->FwdRate[tree_data->CvDiff];
   
    double  **Aweight    = tree_data->Aweight;
    int     *Top         = tree_data->Top1;
    int     *Bottom      = tree_data->Bottom1;
    int     *OutTop      = tree_data->OutTop1;
    int     *OutBottom   = tree_data->OutBottom1;
    int     NbTP         = tree_data->NbTP;
    double  *Length      = tree_data->Length;
    double  *LengthJ     = tree_data->LengthJ;

    /* Other variables */
    double  *StatePr;          /* State Prices at t                    */
    double  *Discount;         /* 1 period discount factor             */
    double  *StatePr1;         /* State Prices at t+1                  */

    double  *StatePrL;         /* Local slice pointers                 */
    double  *DiscountL;
    double  *StatePr1L;
    double  *EDevPriceL;       /* St prices kept for express DEV'ing   */

    double  *TempPointer;

    double  Jump;              /* Size of the jump (in log space)           */
    double  RateJumpLeft;      /* Rate space jump size for last dimension   */
    double  RateJumpRight;     /* Rate space jump size for last dimension   */
    double  PreviousJump;      /* Jump size at previous period              */
    double  JumpCoeff, du;     /* Jump coefficients                         */
    double  Grid;              /* Current grid point                        */

    double  BetaLocal;         /* = Beta * Length[t]                        */
    double  Pi;                /* Drift at the current node                 */
    double  d;                 /* Drift due to the change in jump size      */
    double  Zt=0;              /* Center of the tree in X-space             */
    double  DelZt;             /* Correction to Zt for current time point   */
    double  Zidx;              /* Zt index i,j,k adjusted                   */
    double  P[3];              /* 2rd degree polynomial coefficients        */
    double  p;                 /* Total probability                         */
    double  x;                 /*                                           */
    double  FwdRateA;          /* Fwd rate adjusted                         */
    double  MLeft, SLeft;      /* Multiple coeff and shift for grid point   */
    double  MRight, SRight;
    double  VolBbq;            /* Sigma used in bone mapping                */
    double  Bbq;               /* Backbone parameter                        */
    double  QMid;              /* Average q parameter at current point      */
    double  QLeft;             /* QLeft at current point                    */
    double  QRight;            /* QRight at current point                   */
    double  FwdShift;          /* Forward shift at current point            */
    double  QSh;               /* q parameter for shift adjustment          */
    double  InsRate;           /* Instantaneous rate                        */
    double  Q, U;              /* Q and U weight coeff for Z-polynomial     */ 
    double  D0;                /* Consecutive powers of discount factor     */ 
    double  aL00;              /* Coeff of Z-polynomial across all states   */
    double  aL10, aL11 ;
    double  aL20, aL21, aL22 ;
    double  aR00;               
    double  aR10, aR11 ;
    double  aR20, aR21, aR22 ;

    double  pu, pd, p0;        /* Probabilities                             */

    int     EDevOn = FALSE;    /* To indicate whether EDev tool is used     */
    int     EDevIdx = 0;       /* Counter of dates for express DEV          */

    int     offset;
    int     Mid;               /* Mid of distribution index                 */
    int     i;                 /* Node index                                */
    int     l, lMin, lMax;     /* Node branching shift                      */
    int     t;                 /* Time step index                           */
    int     status = FAILURE;  /* Error status = FAILURE initially          */



    StatePr  = Fix3_Alloc_Slice (tree_data);
    Discount = Fix3_Alloc_Slice (tree_data);
    StatePr1 = Fix3_Alloc_Slice (tree_data);

    if (   (StatePr  == NULL)
        || (Discount == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Fix3_Drift_1D: could not allocate memory!");
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

    Bbq     = mktvol_data->Bbq;
   
            
    StatePrL = StatePr + Fix3_Node_Offset(1, 0, 0, 0, tree_data);
    StatePrL[0] = 1.;

    for (t = 0; t < NbTP; t++)
    {   
        StatePr1L = StatePr1 + Fix3_Node_Offset(1, 0, 0, t+1, tree_data);
        StatePrL  = StatePr  + Fix3_Node_Offset(1, 0, 0, t,   tree_data);
        DiscountL = Discount + Fix3_Node_Offset(1, 0, 0, t,   tree_data);
        QLeft     = tree_data->QLeft[t];
        QRight    = tree_data->QRight[t];
        FwdShift  = tree_data->FwdShift[t];
        QMid      = (QLeft + QRight) / 2;
  

        /* 
        *   Precompute jumps and probabilities coefficients.
        */

        du           = sqrt (JUMPCOEFF * LengthJ[t-1]);
        PreviousJump = Aweight[0][t-1] * du;
                                                                                
        du   = sqrt (JUMPCOEFF * LengthJ[t]);
        Jump = Aweight[0][t] * du;
        d    = PreviousJump / Jump - 1.;
        
        BetaLocal = tree_data->BetaTD[0][t] * Length[t] * PreviousJump / Jump;
        JumpCoeff = Length[t]/(LengthJ[t]*JUMPCOEFF);

        FwdRateA  = FwdRate[t] / (1. + FwdShift);
        InsRate   = FwdRate[t] / Length[t];    

        VolBbq  = mktvol_data->VolLogn * Bbq  
                + mktvol_data->VolNorm *(1. - Bbq) / InsRate;
        VolBbq *= (1. + FwdShift) / (1. + QMid * FwdShift);

        if (IS_Q(QLeft))
        {
            MLeft        = FwdRateA / QLeft;
            SLeft        = 1 + FwdRateA - FwdRateA / QLeft;      
            RateJumpLeft = exp(QLeft * VolBbq * PreviousJump);
        }
        else
        {
            MLeft        = FwdRateA;
            SLeft        = 1 + FwdRateA;
            RateJumpLeft = FwdRateA * VolBbq * PreviousJump;
        }
        if (IS_Q(QRight))
        {
            MRight        = FwdRateA / QRight;
            SRight        = 1 + FwdRateA - FwdRateA / QRight;      
            RateJumpRight = exp(QRight * VolBbq * PreviousJump);
        }
        else
        {
            MRight        = FwdRateA;
            SRight        = 1 + FwdRateA;
            RateJumpRight = FwdRateA * VolBbq * PreviousJump;
        }

        /* Set initial Zt */
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

        /* 
        *   Set up the grid points. 
        */

        /* LEFT part of distribution */

        Zidx  = Zt; 
        Mid   = (int) ceil(-Zidx / PreviousJump) - 1;
        Mid   = MIN ( MAX (Mid, Bottom[t] - 1), Top[t]);
        Zidx += (PreviousJump) * Bottom[t];

        /* printf("(%3d) b=%4d m=%4d t=%4d \n",t ,Bottom[t], Mid, Top[t]); */


        if (IS_Q(QLeft))
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
        
            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;

            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */

        Zidx = Zt 
             + (PreviousJump) * (Mid + 1);

        if (IS_Q(QRight))
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
        
            for (i = Mid + 1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;

            for (i = Mid +1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid += RateJumpRight;
            }
        }


        /*
        *   Calculate polynomial coefficients and drift.
        */

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

        Zidx = Zt; 
        Mid  = (int) ceil(-Zidx / PreviousJump) - 1;
        Mid  = MIN ( MAX (Mid, Bottom[t] - 1), Top[t]);

        for (i = Bottom[t]; i <= Mid; i++)
        {
            x   = DiscountL[i];
            D0  = StatePrL[i] * x;
            P[0] += D0;

            P[1] += (aL10 + aL11 * x) * D0;

            P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  
        }

        for (i = Mid + 1; i <= Top[t]; i++)
        {
            x   = DiscountL[i];
            D0  = StatePrL[i] * x;
            P[0] += D0;

            P[1] += (aR10 + aR11 * x) * D0;

            P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  
        }


        P[0] -= ZeroCoupon[t+1];

        /* Find change in tree centre at current time step*/
        {
            double Root[2] = {0.0, 0.0};

            /* Find root of drift equation (2nd-order expansion in DelZt) */
            if (Quadratic_Solve(Root, P) == FAILURE)
            {
                DR_Error ("Fix3_Drift_1D: no solution exists to drift equation!");
                goto FREE_MEM_AND_RETURN;
            }

            /* Choose root nearest to 0.0 */
            DelZt = ( (fabs(Root[0]) < fabs(Root[1])) ? Root[0] : Root[1] );

            /* Check power series expansion was valid */ 
            if (fabs(DelZt * VolBbq) > 1.)
            {
                DR_Error ("Fix3_Drift_1D: change in drift is too large");
                goto FREE_MEM_AND_RETURN;
            }
        }

        Zt += DelZt;
        ZCenter[t] = Zt;


        /* 
         *  Recalculate the grid using the middle node.
         */

        /* LEFT part of distribution */

        Zidx  = Zt; 
        Mid   = (int) ceil(-Zidx / PreviousJump) - 1;
        Mid   = MIN ( MAX (Mid, Bottom[t] - 1), Top[t]);
        Zidx += (PreviousJump) * Bottom[t];

        if (IS_Q(QLeft))
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
        
            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;

            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */

        Zidx = Zt 
             + (PreviousJump) * (Mid + 1);

        if (IS_Q(QRight))
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
        
            for (i = Mid + 1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;

            for (i = Mid + 1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid += RateJumpRight;
            }
        }


        /* 
        *   Compute state prices for next period.
        */

        for (i = Bottom[t+1]; i <= Top[t+1]; i++)
        {
            StatePr1L[i] = 0.;                        
        }
                

        /* Branching has to be inside outer ellipse at next period */
        lMax = OutTop[t+1] - 1;                                             
        lMin = OutBottom[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Fix3_Drift_1D: problem in building the tree (lMin > lMax)!");
            goto FREE_MEM_AND_RETURN;
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

            x = StatePrL[i] * DiscountL[i];

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

        if (fabs (p - ZeroCoupon[t+1]) > 0.0001)
        {
            DR_Error ("Fix3_Drift_1D: state prices don't add up to 1: "
                        "increase sigma or Ppy!");
            goto FREE_MEM_AND_RETURN;
        }



        /* Express DEV tool: store state prices if necessary   */
        if ((EDevOn) && (EDevIdx<tree_data->NbEDevDates))
        {
            if (tree_data->TPDate[t+1] == tree_data->EDevDate[EDevIdx])
            {

                /* Store corresponding state prices */
                offset = Fix3_Node_Offset(1, 0, 0, t+1, tree_data);
                EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;

                for (i = Bottom[t+1]; i <= Top[t+1]; i++)
                {
                    EDevPriceL[i] = StatePr1L[i];
                }

                /* Increment the EDev date index */
                EDevIdx++;
            }
        }


    

        /* 
        *   Permute values for the next iteration.
        */
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

}  /* Fix3_Drift_1D */



/*****  Fix3_Drift_2D_TimeDep   **********************************************************/
/*
*       Calculate the interest rate drift for a two factor model.
*/
int     Fix3_Drift_2D_TimeDep(  
                    MKTVOL_DATA     *mktvol_data,        /* (I) Volatility data */
                    FIX3_TREE_DATA       *tree_data)     /* (I/O) Tree data     */
{

    /* Local values of tree data structure */
    /* Use diffused zero and fwd curves    */
    double  *ZCenter     = tree_data->ZCenter;
    double  *ZeroCoupon  = tree_data->ZeroCoupon[tree_data->CvDiff];
    double  *FwdRate     = tree_data->FwdRate[tree_data->CvDiff];
    double  QLeft        = mktvol_data->QLeft;
    double  QRight       = mktvol_data->QRight;
    double  FwdShift     = mktvol_data->FwdShift;
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
  
    StatePrL = StatePr + Fix3_Node_Offset(2, 0, 0, 0, tree_data);
    StatePrL[0] = 1.;

    for (t = 0; t < NbTP; t++)
    {                                   
        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
        QLeft    = tree_data->QLeft[t];
        QRight   = tree_data->QRight[t];
        FwdShift = tree_data->FwdShift[t];
        QMid     = (QLeft + QRight) / 2;
        
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
            
        Beta1  = tree_data->BetaTD[0][t] * Length[t] * PreviousJump[0] / Jump[0];
        Beta2  = tree_data->BetaTD[1][t] * Length[t] * PreviousJump[1] / Jump[2];
        Beta3  = tree_data->BetaTD[1][t] * Length[t] * PreviousJump[2] / Jump[2];
        
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
                  + (PreviousJump[0] + PreviousJump[1]) * i;
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
                 + (PreviousJump[0] + PreviousJump[1]) * i
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
                 + (PreviousJump[0] + PreviousJump[1]) * i;
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
            if (fabs(DelZt * VolBbq) > 1.)
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
                  + (PreviousJump[0] + PreviousJump[1]) * i;
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
                 + (PreviousJump[0] + PreviousJump[1]) * i
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

}  /* Fix3_Drift_2D */



/*****  Fix3_Drift_3D_TimeDep   **********************************************************/
/*
*       Calculate the interest rate drift for a three factor model.
*/
int     Fix3_Drift_3D_TimeDep (  MKTVOL_DATA     *mktvol_data,   /* (I) Volatility data */
                    FIX3_TREE_DATA       *tree_data)     /* (I/O) Tree data     */
{

    /* Local values of tree data structure */
    /* Use diffused zero and fwd curves    */
    double  *ZCenter      = tree_data->ZCenter;
    double  *ZeroCoupon   = tree_data->ZeroCoupon[tree_data->CvDiff];
    double  *FwdRate      = tree_data->FwdRate[tree_data->CvDiff];
    double  QLeft         = mktvol_data->QLeft;
    double  QRight        = mktvol_data->QRight;
    double  FwdShift      = mktvol_data->FwdShift;
    double  **Aweight     = tree_data->Aweight;
    int     *Top1         = tree_data->Top1;
    int     *Bottom1      = tree_data->Bottom1;
    int     **Top2        = tree_data->Top2;
    int     **Bottom2     = tree_data->Bottom2;
    int     ***Top3       = tree_data->Top3;
    int     ***Bottom3    = tree_data->Bottom3;
    int     *OutTop1      = tree_data->OutTop1;
    int     *OutBottom1   = tree_data->OutBottom1;
    int     **OutTop2     = tree_data->OutTop2;
    int     **OutBottom2  = tree_data->OutBottom2;
    int     ***OutTop3    = tree_data->OutTop3;
    int     ***OutBottom3 = tree_data->OutBottom3;
    int     NbTP          = tree_data->NbTP;
    double  *Length       = tree_data->Length;
    double  *LengthJ      = tree_data->LengthJ;

    /* Other variables */
    double  *StatePr;            /* State Prices at t                         */
    double  *Discount;           /* 1 period discount factor                  */
    double  *StatePr1;           /* State Prices at t+1                       */    


    double  *StatePrL;           /* Local slice pointers                      */
    double  *DiscountL;
    double  *StatePr11L, *StatePr12L, *StatePr13L;
    double  *StatePr21L, *StatePr22L, *StatePr23L;
    double  *StatePr31L, *StatePr32L, *StatePr33L;
    double  *EDevPriceL;         /* St prices kept for express DEV'ing        */



    double  *TempPointer;

    double  Jump[6];             /* Jump sizes at current time point          */
    double  PreviousJump[6];     /* Jump sizes at previous time point         */
    double  RateJumpLeft5;       /* Rate space jump size for last dimension   */
    double  RateJumpRight5;      /* Rate space jump size for last dimension   */
    double  JumpCoeff, d[6], du; /* Jump coefficients                         */
    double  Grid;                /* Grid points                               */
                                                                              
    double  Beta11, Beta21, Beta22; /* Mean reversion coefficients            */
    double  Beta31, Beta32, Beta33;                                           
    double  Pi, Qi, Qij;            /* Total drift                            */
    double  Ri, Rij, Rijk;                                                    
    double  p;                   /* Total probability                         */
    double  Zt=0;                /* Center of the tree in X-space             */
    double  DelZt;               /* Correction to Zt for current time point   */
    double  Zidx;                /* Zt index i,j,k adjusted                   */
    double  P[3];                /* 2nd degree polynomial coefficients        */
    double  x;                   
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


    double  pu, pd, p0;          /* Probabilities */
    double  qu, qd, q0;
    double  quu, qu0, qud, q0u, q00, q0d, qdu, qd0, qdd; 
    double  ru, rd, r0;

    int     EDevOn = FALSE;      /* To indicate whether EDev tool is used     */
    int     EDevIdx = 0;         /* Counter of dates for express DEV          */


    int     offset;
    int     Mid;                 /* Mid of distribution index                 */
    int     i, j, k, k1, k2, k3; /* Node indices                              */
    int     l, lMin, lMax;       /* Node branching shifts                     */
    int     m, mMin, mMax;       
    int     n, nMin, nMax;
    int     t;                   /* Time step index                           */
    int     status = FAILURE;    /* Error status = FAILURE initially          */


    StatePr  = Fix3_Alloc_Slice (tree_data);
    Discount = Fix3_Alloc_Slice (tree_data);
    StatePr1 = Fix3_Alloc_Slice (tree_data);

    if (   (StatePr  == NULL)
        || (Discount == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Fix3_Drift_3D: could not allocate memory!");
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
   
  
    StatePrL = StatePr + Fix3_Node_Offset(3, 0, 0, 0, tree_data);
    StatePrL[0] = 1.;

    for (t = 0; t < NbTP; t++)               
    {                                   
        du = sqrt (JUMPCOEFF * LengthJ[t-1]);

        QLeft    = tree_data->QLeft[t];
        QRight   = tree_data->QRight[t];
        FwdShift = tree_data->FwdShift[t];
        QMid     = (QLeft + QRight) / 2;
        
        PreviousJump[0] = Aweight[0][t-1] * du;
        PreviousJump[1] = Aweight[1][t-1] * du;
        PreviousJump[2] = Aweight[2][t-1] * du;
        PreviousJump[3] = Aweight[3][t-1] * du;
        PreviousJump[4] = Aweight[4][t-1] * du;
        PreviousJump[5] = Aweight[5][t-1] * du;

        du = sqrt (JUMPCOEFF * LengthJ[t]);
            
        Jump[0] = Aweight[0][t] * du;
        Jump[1] = Aweight[1][t] * du;
        Jump[2] = Aweight[2][t] * du;
        Jump[3] = Aweight[3][t] * du;
        Jump[4] = Aweight[4][t] * du;
        Jump[5] = Aweight[5][t] * du;

        d[0] = (PreviousJump[0] - Jump[0]) / Jump[0];
        d[1] = (PreviousJump[1] - Jump[1]) / Jump[2];
        d[2] = (PreviousJump[2] - Jump[2]) / Jump[2];
        d[3] = (PreviousJump[3] - Jump[3]) / Jump[5];
        d[4] = (PreviousJump[4] - Jump[4]) / Jump[5];
        d[5] = (PreviousJump[5] - Jump[5]) / Jump[5];
            
        Beta11 = tree_data->BetaTD[0][t] * Length[t] * PreviousJump[0] / Jump[0];
        Beta21 = tree_data->BetaTD[1][t] * Length[t] * PreviousJump[1] / Jump[2];
        Beta22 = tree_data->BetaTD[1][t] * Length[t] * PreviousJump[2] / Jump[2];
        Beta31 = tree_data->BetaTD[2][t] * Length[t] * PreviousJump[3] / Jump[5];
        Beta32 = tree_data->BetaTD[2][t] * Length[t] * PreviousJump[4] / Jump[5];
        Beta33 = tree_data->BetaTD[2][t] * Length[t] * PreviousJump[5] / Jump[5];
        
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
            RateJumpLeft5 = exp(QLeft * VolBbq * PreviousJump[5]);
        }
        else
        {
            MLeft         = FwdRateA;
            SLeft         = 1 + FwdRateA;
            RateJumpLeft5 = FwdRateA * VolBbq * PreviousJump[5];
        }
        if (IS_Q(QRight))
        {
            MRight         = FwdRateA / QRight;
            SRight         = 1 + FwdRateA - FwdRateA / QRight;      
            RateJumpRight5 = exp(QRight * VolBbq * PreviousJump[5]);
        }
        else
        {
            MRight         = FwdRateA;
            SRight         = 1 + FwdRateA;
            RateJumpRight5 = FwdRateA * VolBbq * PreviousJump[5];
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
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                /* LEFT part of distribution */

                Zidx  = Zt 
                      + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                      + (PreviousJump[2] + PreviousJump[4]) * j;
                Mid   = (int) ceil(-Zidx / PreviousJump[5]) - 1;
                Mid   = MIN ( MAX (Mid, Bottom3[t][i][j] - 1), Top3[t][i][j]);
                Zidx += (PreviousJump[5]) * Bottom3[t][i][j];
            
                DiscountL = Discount + Fix3_Node_Offset (3, i, j, t, tree_data);
                
                if (IS_Q(QLeft))
                {
                    Grid = MLeft * exp (QLeft * VolBbq * Zidx);
    
                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid *= RateJumpLeft5;
                    }
                }
                else
                {
                    Grid = MLeft * VolBbq * Zidx;

                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid += RateJumpLeft5;
                    }
                }

                /* Right part of distribution */

                Zidx = Zt 
                     + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                     + (PreviousJump[2] + PreviousJump[4]) * j
                     + (PreviousJump[5]) * (Mid + 1);

                if (IS_Q(QRight))
                {
                    Grid = MRight * exp (QRight * VolBbq * Zidx);
    
                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid *= RateJumpRight5;
                    }
                }
                else
                {
                    Grid = MRight * VolBbq * Zidx;

                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid += RateJumpRight5;
                    }
                }

            }  /* for j */     
        
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
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                StatePrL  = StatePr  + Fix3_Node_Offset (3, i, j, t, tree_data);
                DiscountL = Discount + Fix3_Node_Offset (3, i, j, t, tree_data);
                    
                Zidx = Zt 
                     + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                     + (PreviousJump[2] + PreviousJump[4]) * j;
                Mid  = (int) ceil(-Zidx / PreviousJump[5]) - 1;
                Mid  = MIN ( MAX (Mid, Bottom3[t][i][j] - 1), Top3[t][i][j]);

                for (k = Bottom3[t][i][j]; k <= Mid; k++)
                {

                    x   = DiscountL[k];
                    D0  = StatePrL[k] * x;
                    P[0] += D0;

                    P[1] += (aL10 + aL11 * x) * D0;
                              
                    P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  

                }

                for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                {

                    x   = DiscountL[k];
                    D0  = StatePrL[k] * x;
                    P[0] += D0;

                    P[1] += (aR10 + aR11 * x) * D0;
                              
                    P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  

                }

            }  /* for j */
        }  /* for i */
        

        P[0] -= ZeroCoupon[t+1];
                                        
        /* Find change in tree centre at current time step*/
        {
            double Root[2] = {0.0, 0.0};

            /* Find root of drift equation (2nd-order expansion in DelZt) */
            if (Quadratic_Solve(Root, P) == FAILURE)
            {
                DR_Error ("Fix3_Drift_3D: no solution exists to drift equation!");
                goto FREE_MEM_AND_RETURN;
            }

            /* Choose root nearest to 0.0 */
            DelZt = ( (fabs(Root[0]) < fabs(Root[1])) ? Root[0] : Root[1] );

            /* Check power series expansion was valid */ 
            if (fabs(DelZt * VolBbq) > 1.)
            {
                DR_Error ("Fix3_Drift_3D: change in drift is too large");
                goto FREE_MEM_AND_RETURN;
            }
        }


        Zt += DelZt;
        ZCenter[t] = Zt;


        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                /* LEFT part of distribution */

                Zidx  = Zt 
                      + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                      + (PreviousJump[2] + PreviousJump[4]) * j;
                Mid   = (int) ceil(-Zidx / PreviousJump[5]) - 1;
                Mid   = MIN ( MAX (Mid, Bottom3[t][i][j] - 1), Top3[t][i][j]);
                Zidx += (PreviousJump[5]) * Bottom3[t][i][j];
            
                DiscountL = Discount + Fix3_Node_Offset (3, i, j, t, tree_data);
                
                if (IS_Q(QLeft))
                {
                    Grid = MLeft * exp (QLeft * VolBbq * Zidx);
    
                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid *= RateJumpLeft5;
                    }
                }
                else
                {
                    Grid = MLeft * VolBbq * Zidx;

                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid += RateJumpLeft5;
                    }
                }

                /* Right part of distribution */

                Zidx = Zt 
                     + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                     + (PreviousJump[2] + PreviousJump[4]) * j
                     + (PreviousJump[5]) * (Mid + 1);

                if (IS_Q(QRight))
                {
                    Grid = MRight * exp (QRight * VolBbq * Zidx);
    
                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid *= RateJumpRight5;
                    }
                }
                else
                {
                    Grid = MRight * VolBbq * Zidx;

                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid += RateJumpRight5;
                    }
                }

            }  /* for j */     
        
        }  /* for i */

                
        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr11L = StatePr1 + Fix3_Node_Offset (3, i, j, t+1, tree_data);
            
                for (k = Bottom3[t+1][i][j]; k <= Top3[t+1][i][j]; k++)
                {
                    StatePr11L[k] = 0.;
                }
            }  /* for j */
        }  /* for i */
                
                        
        lMax = OutTop1[t+1] - 1;
        lMin = OutBottom1[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Fix3_Drift_3D: problem in building the tree (lMin > lMax)!");
            goto FREE_MEM_AND_RETURN;
        }

        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {                        
            Pi = (d[0] - Beta11) * i;               
            Qi = (d[1] - Beta21) * i - Pi * Jump[1] / Jump[2];
            Ri = (d[3] - Beta31) * i - Pi * Jump[3] / Jump[5];

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
                DR_Error ("Fix3_Drift_3D: problem in building the tree (mMin > mMax)!");
                goto FREE_MEM_AND_RETURN;
            }


            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                Qij = Qi + (d[2] - Beta22) * j;
                Rij = Ri + (d[4] - Beta32) * j - Qij * Jump[4] / Jump[5];

                m = NEAR_INT (Qij);
                m = MIN (MAX (mMin - j, m), mMax - j);
                
                Qij -= m;
                                
                qu = .5 * (JumpCoeff + Qij + Qij * Qij);
                qd = qu - Qij;
                q0 = 1. - qu - qd;

                quu = pu * qu; qu0 = pu * q0; qud = pu * qd;
                q0u = p0 * qu; q00 = p0 * q0; q0d = p0 * qd;
                qdu = pd * qu; qd0 = pd * q0; qdd = pd * qd;

                        
                nMax =            OutTop3[t+1][i+l-1][j+m-1];
                nMax = MIN (nMax, OutTop3[t+1][i+l-1][j+m  ]);
                nMax = MIN (nMax, OutTop3[t+1][i+l-1][j+m+1]);
                nMax = MIN (nMax, OutTop3[t+1][i+l  ][j+m-1]);
                nMax = MIN (nMax, OutTop3[t+1][i+l  ][j+m  ]);
                nMax = MIN (nMax, OutTop3[t+1][i+l  ][j+m+1]);
                nMax = MIN (nMax, OutTop3[t+1][i+l+1][j+m-1]);
                nMax = MIN (nMax, OutTop3[t+1][i+l+1][j+m  ]);
                nMax = MIN (nMax, OutTop3[t+1][i+l+1][j+m+1]) - 1;

                nMin =            OutBottom3[t+1][i+l-1][j+m-1];
                nMin = MAX (nMin, OutBottom3[t+1][i+l-1][j+m  ]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l-1][j+m+1]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l  ][j+m-1]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l  ][j+m  ]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l  ][j+m+1]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l+1][j+m-1]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l+1][j+m  ]);
                nMin = MAX (nMin, OutBottom3[t+1][i+l+1][j+m+1]) + 1;

                if (nMin > nMax)
                {
                    DR_Error ("Fix3_Drift_3D: problem in building the tree (nMin > nMax)!");
                    goto FREE_MEM_AND_RETURN;
                }


                StatePrL  = StatePr  + Fix3_Node_Offset (3, i, j, t, tree_data);
                DiscountL = Discount + Fix3_Node_Offset (3, i, j, t, tree_data);
                    
                StatePr11L = StatePr1 + Fix3_Node_Offset (3, i+1+l, j+1+m, t+1, tree_data);
                StatePr12L = StatePr1 + Fix3_Node_Offset (3, i+1+l, j  +m, t+1, tree_data);
                StatePr13L = StatePr1 + Fix3_Node_Offset (3, i+1+l, j-1+m, t+1, tree_data);
                StatePr21L = StatePr1 + Fix3_Node_Offset (3, i  +l, j+1+m, t+1, tree_data);
                StatePr22L = StatePr1 + Fix3_Node_Offset (3, i  +l, j  +m, t+1, tree_data);
                StatePr23L = StatePr1 + Fix3_Node_Offset (3, i  +l, j-1+m, t+1, tree_data);
                StatePr31L = StatePr1 + Fix3_Node_Offset (3, i-1+l, j+1+m, t+1, tree_data);
                StatePr32L = StatePr1 + Fix3_Node_Offset (3, i-1+l, j  +m, t+1, tree_data);
                StatePr33L = StatePr1 + Fix3_Node_Offset (3, i-1+l, j-1+m, t+1, tree_data);
                
                
                for (k = Bottom3[t][i][j]; k <= Top3[t][i][j]; k++)
                {
                    Rijk = Rij + (d[5] - Beta33) * k;

                    n = NEAR_INT (Rijk);
                    n = MIN (MAX (nMin - k, n), nMax - k);
                
                    Rijk -= n;
                
                    ru = .5 * (JumpCoeff + Rijk + Rijk * Rijk);
                    rd = ru - Rijk;
                    r0 = 1. - ru - rd;


                    k1 = k + 1 + n;
                    k2 = k     + n;
                    k3 = k - 1 + n;

                                
                    x = StatePrL[k] * DiscountL[k];

                    StatePr11L[k1] += quu * ru * x;
                    StatePr11L[k2] += quu * r0 * x;
                    StatePr11L[k3] += quu * rd * x;
                    StatePr12L[k1] += qu0 * ru * x;
                    StatePr12L[k2] += qu0 * r0 * x;
                    StatePr12L[k3] += qu0 * rd * x;
                    StatePr13L[k1] += qud * ru * x;
                    StatePr13L[k2] += qud * r0 * x;
                    StatePr13L[k3] += qud * rd * x;
                    StatePr21L[k1] += q0u * ru * x;
                    StatePr21L[k2] += q0u * r0 * x;
                    StatePr21L[k3] += q0u * rd * x;
                    StatePr22L[k1] += q00 * ru * x;
                    StatePr22L[k2] += q00 * r0 * x;
                    StatePr22L[k3] += q00 * rd * x;
                    StatePr23L[k1] += q0d * ru * x;
                    StatePr23L[k2] += q0d * r0 * x;
                    StatePr23L[k3] += q0d * rd * x;
                    StatePr31L[k1] += qdu * ru * x;
                    StatePr31L[k2] += qdu * r0 * x;
                    StatePr31L[k3] += qdu * rd * x;
                    StatePr32L[k1] += qd0 * ru * x;
                    StatePr32L[k2] += qd0 * r0 * x;
                    StatePr32L[k3] += qd0 * rd * x;
                    StatePr33L[k1] += qdd * ru * x;
                    StatePr33L[k2] += qdd * r0 * x;
                    StatePr33L[k3] += qdd * rd * x;
                }
            }  /* for j */
        }  /* for i */

                             
        p = 0.;

        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr11L = StatePr1 + Fix3_Node_Offset (3, i, j, t+1, tree_data);
            
                for (k = Bottom3[t+1][i][j]; k <= Top3[t+1][i][j]; k++)
                {
                    p += StatePr11L[k];
                }
            }  /* for j */
        }  /* for i */

        if (fabs (p - ZeroCoupon[t+1]) > 0.0001)
        {
            DR_Error ("Fix3_Drift_3D: state prices don't add up to 1: "
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
                    for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
                    {   
                        offset =  Fix3_Node_Offset (3, i, j, t+1, tree_data);                               
                        StatePr11L = StatePr1 + offset;
                        EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;
            
                        for (k = Bottom3[t+1][i][j]; k <= Top3[t+1][i][j]; k++)
                        {
                            EDevPriceL[k] = StatePr11L[k];
                        }
                    }  /* for j */
                }  /* for i */


                /* Increment the EDev date index */
                EDevIdx++;
            }
        }

    
                  
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

}  /* Fix3_Drift_3D */
