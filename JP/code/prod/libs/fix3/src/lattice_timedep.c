/****************************************************************************/
/*      Build the lattice.                                                  */
/****************************************************************************/
/*      LATTICE.c                                                           */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  Fix3_Lattice  **************************************************************/
/*
*       Position the nodes in the lattice: calculate the one period discount
*       factor and the probabilities.
*/
int    Fix3_Lattice_TimeDep (FIX3_DEV_DATA        *dev_data,     /* (O) Fix3_Dev data structure     */
                int             t,             /* (I) Current time point     */
                int             T,             /* (I) Last time point        */
                MKTVOL_DATA     *mktvol_data,  /* (I) Volatility data        */
                FIX3_TREE_DATA       *tree_data)    /* (I) Tree data structure    */
{

    double  *DiscountL[3];              /* Local slice pointers */

    double  *puL,  *pdL,  *p0L;
    double  *quuL, *qu0L, *qudL;
    double  *q0uL, *q00L, *q0dL;
    double  *qduL, *qd0L, *qddL;
    double  *ruL,  *r0L,  *rdL;

    double  pu, pd, p0;                 /* Local probabilities */
    double  qu, qd, q0;

    double  Grid;                       /* Grid points                      */
    double  ZRatio1, ZRatio2;           /* Zero coupon ratios               */
    double  FwdRate;                    /* Fwd rate for current time point  */
    double  Length;                     /* Length of time steps             */
    double  LengthJ;                    /* Length of time steps for jump    */
    double  ZCenter;                    /* Center of the tree in X-space    */
    double  Zidx;                       /* Zt index i,j,k adjusted          */
    double  FwdRateA;                   /* Fwd rate adjusted                */
    double  MLeft, SLeft;               /* Coeff and shift for grid point   */
    double  MRight, SRight;
    double  Jump[6];                    /* Jump sizes at current time step  */
    double  PreviousJump[6];            /* Jump sizes at previous time step */
    double  RateJumpLeft;               /* Rate space jump size for last dim*/
    double  RateJumpRight;              /* Rate space jump size for last dim*/
    double  JumpCoeff, d[6], du;        /* Jump coefficients                */
    double  QLeft;                      /* Left q mapping coefficient       */
    double  QRight;                     /* Right q mapping coefficient      */
    double  FwdShift;                   /* Fwd shift mapping coefficient    */
    double  VolBbq;                     /* Sigma used in bone mapping       */
    double  Bbq;                        /* Backbone parameter               */
    double  QMid;                       /* Average q parameter              */
    double  InsRate;                    /* Instantaneous rate               */
    double  Beta11, Beta21, Beta22; 
    double  Beta31, Beta32, Beta33;     /* Modified mean reversion */
    double  Pi, Qi, Qij, Ri, Rij, Rijk; /* Total drift             */

    int     *Shift1L, *Shift2L, *Shift3L;

    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */
    int     OutTop1, OutBottom1;        /* Outer tree limits      */
    int     *OutTop2, *OutBottom2;
    int     **OutTop3, **OutBottom3;

    int     CvDiff, CvIdx1, CvIdx2;     /* Nb of zero curves      */


    int     i, j, k;                    /* Node indices           */
    int     Mid;                        /* Mid of dist index      */
    int     offset;                     /* Node offset            */
    int     l, lMin, lMax;              /* Node branching shifts  */
    int     m, mMin, mMax;
    int     n, nMin, nMax;

    int     status = FAILURE;           /* Error status = FAILURE initially */


    /* Nothing to do at the back of the tree */

    if (t == T)                                                        
    {
        return (SUCCESS);
    }

    /*
    *   Internal assigments of zero curves.
    */
    CvDiff = tree_data->CvDiff;
    CvIdx1 = tree_data->CvIdx1;
    CvIdx2 = tree_data->CvIdx2;

    /* 
    *	Drifts and coefficients used to calculate probas.
    */

    ZRatio1 = (1.+tree_data->FwdRate[CvDiff][t])/(1.+tree_data->FwdRate[CvIdx1][t]);
    ZRatio2 = (1.+tree_data->FwdRate[CvDiff][t])/(1.+tree_data->FwdRate[CvIdx2][t]);
        
    /* 
    *   Previous period coefficients.
    */

    LengthJ = tree_data->LengthJ[t-1];                    
    du = sqrt (JUMPCOEFF * LengthJ);
    
    PreviousJump[0] = tree_data->Aweight[0][t-1] * du;
    PreviousJump[1] = tree_data->Aweight[1][t-1] * du;
    PreviousJump[2] = tree_data->Aweight[2][t-1] * du;
    PreviousJump[3] = tree_data->Aweight[3][t-1] * du;
    PreviousJump[4] = tree_data->Aweight[4][t-1] * du;
    PreviousJump[5] = tree_data->Aweight[5][t-1] * du;

    /* 
    *   Current period coefficients.
    */

    LengthJ   = tree_data->LengthJ[t];                    
    du        = sqrt (JUMPCOEFF * LengthJ);
    Length    = tree_data->Length[t];
    JumpCoeff = Length/(LengthJ*JUMPCOEFF);
    FwdRate   = tree_data->FwdRate[CvDiff][t];
    ZCenter   = tree_data->ZCenter[t];

    Jump[0] = tree_data->Aweight[0][t] * du;
    Jump[1] = tree_data->Aweight[1][t] * du;
    Jump[2] = tree_data->Aweight[2][t] * du;
    Jump[3] = tree_data->Aweight[3][t] * du;
    Jump[4] = tree_data->Aweight[4][t] * du;
    Jump[5] = tree_data->Aweight[5][t] * du;

    d[0] = (PreviousJump[0] - Jump[0]) / Jump[0];
    d[1] = (PreviousJump[1] - Jump[1]) / Jump[2];
    d[2] = (PreviousJump[2] - Jump[2]) / Jump[2];
    d[3] = (PreviousJump[3] - Jump[3]) / Jump[5];
    d[4] = (PreviousJump[4] - Jump[4]) / Jump[5];
    d[5] = (PreviousJump[5] - Jump[5]) / Jump[5];
    
    Beta11 = tree_data->BetaTD[0][t] * Length * PreviousJump[0] / Jump[0];
    Beta21 = tree_data->BetaTD[1][t] * Length * PreviousJump[1] / Jump[2];
    Beta22 = tree_data->BetaTD[1][t] * Length * PreviousJump[2] / Jump[2];
    Beta31 = tree_data->BetaTD[2][t] * Length * PreviousJump[3] / Jump[5];
    Beta32 = tree_data->BetaTD[2][t] * Length * PreviousJump[4] / Jump[5];
    Beta33 = tree_data->BetaTD[2][t] * Length * PreviousJump[5] / Jump[5];
    
        
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    OutTop1    = tree_data->OutTop1[t+1];
    OutTop2    = tree_data->OutTop2[t+1];
    OutTop3    = tree_data->OutTop3[t+1];
    OutBottom1 = tree_data->OutBottom1[t+1];
    OutBottom2 = tree_data->OutBottom2[t+1];
    OutBottom3 = tree_data->OutBottom3[t+1];

    /* Q distribution parameters */
    QLeft    = tree_data->QLeft[t];
    QRight   = tree_data->QRight[t];
    FwdShift = tree_data->FwdShift[t];
        
    /* Mapping sigma */
    Bbq     = mktvol_data->Bbq;
    InsRate = FwdRate / Length;    
    QMid    = (QLeft + QRight) / 2;
    VolBbq  = mktvol_data->VolLogn * Bbq  
            + mktvol_data->VolNorm *(1. - Bbq) / InsRate;
    VolBbq *= (1. + FwdShift) / (1. + QMid * FwdShift);

    FwdRateA  = FwdRate / (1. + FwdShift);

    if (IS_Q(QLeft))
    {
        MLeft = FwdRateA / QLeft;
        SLeft = 1 + FwdRateA - FwdRateA / QLeft;      
    }
    else
    {
        MLeft = FwdRateA;
        SLeft = 1 + FwdRateA;
    }
    if (IS_Q(QRight))
    {
        MRight = FwdRateA / QRight;
        SRight = 1 + FwdRateA - FwdRateA / QRight;      
    }
    else
    {
        MRight = FwdRateA;
        SRight = 1 + FwdRateA;
    }

                
    /* 
    *   Calculate interest rate grids. Only the index curve is diffused.
    *   Christian's approximation is used for the other two curves.
    */

    if (tree_data->NbFactor == 1)
    {
        /* Grid jumps in both part of distribution */        

        if (IS_Q(QLeft))  RateJumpLeft = exp(QLeft * VolBbq * PreviousJump[0]);
        else              RateJumpLeft =  FwdRateA * VolBbq * PreviousJump[0];

        if (IS_Q(QRight)) RateJumpRight = exp(QRight * VolBbq * PreviousJump[0]);
        else              RateJumpRight =   FwdRateA * VolBbq * PreviousJump[0];    


        /* Discount slice offseting */

        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        DiscountL[0] = dev_data->Discount[0] + offset;
        DiscountL[1] = dev_data->Discount[1] + offset;
        DiscountL[2] = dev_data->Discount[2] + offset;

        /* LEFT part of distribution */

        Zidx  = ZCenter; 
        Mid   = (int) ceil(-Zidx / PreviousJump[0]) - 1;
        Mid   = MIN ( MAX (Mid, Bottom1 - 1), Top1);
        Zidx += (PreviousJump[0]) * Bottom1;

        /* printf("(%3d) b=%4d m=%4d t=%4d \n",t ,Bottom1, Mid, Top1); */

        if (IS_Q(QLeft))
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
        
            for (i = Bottom1; i <= Mid; i++)
            {
                DiscountL[CvDiff][i] = 1. / (SLeft + Grid);

                DiscountL[CvIdx1][i] = ZRatio1 * DiscountL[CvDiff][i];
                DiscountL[CvIdx2][i] = ZRatio2 * DiscountL[CvDiff][i];

                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;

            for (i = Bottom1; i <= Mid; i++)
            {
                DiscountL[CvDiff][i] = 1. / (SLeft + Grid);

                DiscountL[CvIdx1][i] = ZRatio1 * DiscountL[CvDiff][i];
                DiscountL[CvIdx2][i] = ZRatio2 * DiscountL[CvDiff][i];
                
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */

        Zidx = ZCenter 
             + (PreviousJump[0]) * (Mid + 1);

        if (IS_Q(QRight))
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
        
            for (i = Mid + 1; i <= Top1; i++)
            {
                DiscountL[CvDiff][i] = 1. / (SRight + Grid);

                DiscountL[CvIdx1][i] = ZRatio1 * DiscountL[CvDiff][i];
                DiscountL[CvIdx2][i] = ZRatio2 * DiscountL[CvDiff][i];

                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;

            for (i = Mid +1; i <= Top1; i++)
            {
                DiscountL[CvDiff][i] = 1. / (SRight + Grid);

                DiscountL[CvIdx1][i] = ZRatio1 * DiscountL[CvDiff][i];
                DiscountL[CvIdx2][i] = ZRatio2 * DiscountL[CvDiff][i];
                
                Grid += RateJumpRight;
            }
        }

        puL = dev_data->pu + offset;
        p0L = dev_data->p0 + offset;
        pdL = dev_data->pd + offset;

        Shift1L = dev_data->Shift1 + offset;

        lMax = OutTop1 - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Fix3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }

        for (i = Bottom1; i <= Top1; i++)                   
        {
            Pi = (d[0] - Beta11) * i;

            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                                                                                
            Shift1L[i] = l;

            Pi -= l;

            puL[i] = .5 * (JumpCoeff + Pi + Pi * Pi);
            pdL[i] = puL[i] - Pi;
            p0L[i] = 1. - puL[i] - pdL[i];
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        /* Grid jumps in both part of distribution */
        
        if (IS_Q(QLeft))  RateJumpLeft = exp(QLeft * VolBbq * PreviousJump[2]);
        else              RateJumpLeft =  FwdRateA * VolBbq * PreviousJump[2];

        if (IS_Q(QRight)) RateJumpRight = exp(QRight * VolBbq * PreviousJump[2]);
        else              RateJumpRight =   FwdRateA * VolBbq * PreviousJump[2];    

        for (i = Bottom1; i <= Top1; i++)
        {
                
            /* Discount slice offseting */
    
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            DiscountL[0] = dev_data->Discount[0] + offset;
            DiscountL[1] = dev_data->Discount[1] + offset;
            DiscountL[2] = dev_data->Discount[2] + offset;

            /* LEFT part of distribution */

            Zidx  = ZCenter 
                  + (PreviousJump[0] + PreviousJump[1]) * i;
            Mid   = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[i] - 1), Top2[i]);
            Zidx += (PreviousJump[2]) * Bottom2[i];

            if (IS_Q(QLeft))
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[i]; j <= Mid; j++)
                {
                    DiscountL[CvDiff][j] = 1. / (SLeft + Grid);

                    DiscountL[CvIdx1][j] = ZRatio1 * DiscountL[CvDiff][j];
                    DiscountL[CvIdx2][j] = ZRatio2 * DiscountL[CvDiff][j];

                    Grid *= RateJumpLeft;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[i]; j <= Mid; j++)
                {
                    DiscountL[CvDiff][j] = 1. / (SLeft + Grid);

                    DiscountL[CvIdx1][j] = ZRatio1 * DiscountL[CvDiff][j];
                    DiscountL[CvIdx2][j] = ZRatio2 * DiscountL[CvDiff][j];

                    Grid += RateJumpLeft;
                }
            }

            /* Right part of distribution */

            Zidx = ZCenter 
                 + (PreviousJump[0] + PreviousJump[1]) * i
                 + (PreviousJump[2]) * (Mid + 1);

            if (IS_Q(QRight))
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[i]; j++)
                {
                    DiscountL[CvDiff][j] = 1. / (SRight + Grid);

                    DiscountL[CvIdx1][j] = ZRatio1 * DiscountL[CvDiff][j];
                    DiscountL[CvIdx2][j] = ZRatio2 * DiscountL[CvDiff][j];

                    Grid *= RateJumpRight;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[i]; j++)
                {
                    DiscountL[CvDiff][j] = 1. / (SRight + Grid);

                    DiscountL[CvIdx1][j] = ZRatio1 * DiscountL[CvDiff][j];
                    DiscountL[CvIdx2][j] = ZRatio2 * DiscountL[CvDiff][j];

                    Grid += RateJumpRight;
                }
            }

        }  /* for i */


        Shift1L = dev_data->Shift1 + Fix3_Node_Offset(1, 0, 0, t, tree_data);

        lMax = OutTop1 - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Fix3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }

        for (i = Bottom1; i <= Top1; i++)
        {                        
            Pi = (d[0] - Beta11) * i;               
            Qi = (d[1] - Beta21) * i - Pi * Jump[1] / Jump[2];

            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                
            Shift1L[i] = l;

            Pi -= l;
                                
            pu = .5 * (JumpCoeff + Pi + Pi * Pi);
            pd = pu - Pi;
            p0 = 1. - pu - pd;
                                

            mMax =            OutTop2[i+l-1];
            mMax = MIN (mMax, OutTop2[i+l  ]);
            mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

            mMin =            OutBottom2[i+l-1];
            mMin = MAX (mMin, OutBottom2[i+l  ]);
            mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;

            if (mMin > mMax)
            {
                DR_Error ("Fix3_Lattice: problem in building the tree (mMin > mMax)!");
                goto RETURN;
            }


            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            quuL = dev_data->quu + offset;
            qu0L = dev_data->qu0 + offset;
            qudL = dev_data->qud + offset;
            q0uL = dev_data->q0u + offset;
            q00L = dev_data->q00 + offset;
            q0dL = dev_data->q0d + offset;
            qduL = dev_data->qdu + offset;
            qd0L = dev_data->qd0 + offset;
            qddL = dev_data->qdd + offset;

            Shift2L = dev_data->Shift2 + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                Qij = Qi + (d[2] - Beta22) * j;

                m = NEAR_INT (Qij);
                m = MIN (MAX (mMin - j, m), mMax - j);
                
                Shift2L[j] = m;

                Qij -= m;
                                
                qu = .5 * (JumpCoeff + Qij + Qij * Qij);
                qd = qu - Qij;
                q0 = 1. - qu - qd;

                /* 
                 *  Store the 9 probabilities instead of 3 p + 3 q to save 
                 *  time in Fix3_Dev ().
                 */

                quuL[j] = pu * qu; qu0L[j] = pu * q0; qudL[j] = pu * qd;
                q0uL[j] = p0 * qu; q00L[j] = p0 * q0; q0dL[j] = p0 * qd;
                qduL[j] = pd * qu; qd0L[j] = pd * q0; qddL[j] = pd * qd;
            }
        }  /* for i */
    }
    else if (tree_data->NbFactor == 3)
    {
        /* Grid jumps in both part of distribution */
        
        if (IS_Q(QLeft))  RateJumpLeft = exp(QLeft * VolBbq * PreviousJump[5]);
        else              RateJumpLeft =  FwdRateA * VolBbq * PreviousJump[5];

        if (IS_Q(QRight)) RateJumpRight = exp(QRight * VolBbq * PreviousJump[5]);
        else              RateJumpRight =   FwdRateA * VolBbq * PreviousJump[5];    

        for (i = Bottom1; i <= Top1; i++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                /* Discount slice offseting */
    
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                DiscountL[0] = dev_data->Discount[0] + offset;
                DiscountL[1] = dev_data->Discount[1] + offset;
                DiscountL[2] = dev_data->Discount[2] + offset;

                /* LEFT part of distribution */

                Zidx  = ZCenter 
                      + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                      + (PreviousJump[2] + PreviousJump[4]) * j;
                Mid   = (int) ceil(-Zidx / PreviousJump[5]) - 1;
                Mid   = MIN ( MAX (Mid, Bottom3[i][j] - 1), Top3[i][j]);
                Zidx += (PreviousJump[5]) * Bottom3[i][j];
            
                if (IS_Q(QLeft))
                {
                    Grid = MLeft * exp (QLeft * VolBbq * Zidx);
    
                    for (k = Bottom3[i][j]; k <= Mid; k++)
                    {
                        DiscountL[CvDiff][k] = 1. / (SLeft + Grid);
            
                        DiscountL[CvIdx1][k] = ZRatio1 * DiscountL[CvDiff][k];
                        DiscountL[CvIdx2][k] = ZRatio2 * DiscountL[CvDiff][k];
                
                        Grid *= RateJumpLeft;
                    }
                }
                else
                {
                    Grid = MLeft * VolBbq * Zidx;

                    for (k = Bottom3[i][j]; k <= Mid; k++)
                    {
                        DiscountL[CvDiff][k] = 1. / (SLeft + Grid);
            
                        DiscountL[CvIdx1][k] = ZRatio1 * DiscountL[CvDiff][k];
                        DiscountL[CvIdx2][k] = ZRatio2 * DiscountL[CvDiff][k];
                
                        Grid += RateJumpLeft;
                    }
                }

                /* Right part of distribution */

                Zidx = ZCenter 
                     + (PreviousJump[0] + PreviousJump[1] + PreviousJump[3]) * i
                     + (PreviousJump[2] + PreviousJump[4]) * j
                     + (PreviousJump[5]) * (Mid + 1);

                if (IS_Q(QRight))
                {
                    Grid = MRight * exp (QRight * VolBbq * Zidx);
    
                    for (k = Mid + 1; k <= Top3[i][j]; k++)
                   {
                        DiscountL[CvDiff][k] = 1. / (SRight + Grid);
            
                        DiscountL[CvIdx1][k] = ZRatio1 * DiscountL[CvDiff][k];
                        DiscountL[CvIdx2][k] = ZRatio2 * DiscountL[CvDiff][k];
                
                        Grid *= RateJumpRight;
                    }
                }
                else
                {
                    Grid = MRight * VolBbq * Zidx;

                    for (k = Mid + 1; k <= Top3[i][j]; k++)
                    {
                        DiscountL[CvDiff][k] = 1. / (SRight + Grid);
            
                        DiscountL[CvIdx1][k] = ZRatio1 * DiscountL[CvDiff][k];
                        DiscountL[CvIdx2][k] = ZRatio2 * DiscountL[CvDiff][k];
                
                        Grid += RateJumpRight;
                    }
                }

            }  /* for j */     
        
        }  /* for i */

        Shift1L = dev_data->Shift1 + Fix3_Node_Offset(1, 0, 0, t, tree_data);

        lMax = OutTop1 - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Fix3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }

        for (i = Bottom1; i <= Top1; i++)
        {                        
            Pi = (d[0] - Beta11) * i;               
            Qi = (d[1] - Beta21) * i - Pi * Jump[1] / Jump[2];
            Ri = (d[3] - Beta31) * i - Pi * Jump[3] / Jump[5];

            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                
            Shift1L[i] = l;

            Pi -= l;
                                
            pu = .5 * (JumpCoeff + Pi + Pi * Pi);
            pd = pu - Pi;
            p0 = 1. - pu - pd;
                                
                                
            mMax =            OutTop2[i+l-1];
            mMax = MIN (mMax, OutTop2[i+l  ]);
            mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

            mMin =            OutBottom2[i+l-1];
            mMin = MAX (mMin, OutBottom2[i+l  ]);
            mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;

            if (mMin > mMax)
            {
                DR_Error ("Fix3_Lattice: problem in building the tree (mMin > mMax)!");
                goto RETURN;
            }


            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            quuL = dev_data->quu + offset;
            qu0L = dev_data->qu0 + offset;
            qudL = dev_data->qud + offset;
            q0uL = dev_data->q0u + offset;
            q00L = dev_data->q00 + offset;
            q0dL = dev_data->q0d + offset;
            qduL = dev_data->qdu + offset;
            qd0L = dev_data->qd0 + offset;
            qddL = dev_data->qdd + offset;

            Shift2L = dev_data->Shift2 + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                Qij = Qi + (d[2] - Beta22) * j;
                Rij = Ri + (d[4] - Beta32) * j - Qij * Jump[4] / Jump[5];

                m = NEAR_INT (Qij);
                m = MIN (MAX (mMin - j, m), mMax - j);
                
                Shift2L[j] = m;

                Qij -= m;
                                
                qu = .5 * (JumpCoeff + Qij + Qij * Qij);
                qd = qu - Qij;
                q0 = 1. - qu - qd;

                quuL[j] = pu * qu; qu0L[j] = pu * q0; qudL[j] = pu * qd;
                q0uL[j] = p0 * qu; q00L[j] = p0 * q0; q0dL[j] = p0 * qd;
                qduL[j] = pd * qu; qd0L[j] = pd * q0; qddL[j] = pd * qd;


                nMax =            OutTop3[i+l-1][j+m-1];
                nMax = MIN (nMax, OutTop3[i+l-1][j+m  ]);
                nMax = MIN (nMax, OutTop3[i+l-1][j+m+1]);
                nMax = MIN (nMax, OutTop3[i+l  ][j+m-1]);
                nMax = MIN (nMax, OutTop3[i+l  ][j+m  ]);
                nMax = MIN (nMax, OutTop3[i+l  ][j+m+1]);
                nMax = MIN (nMax, OutTop3[i+l+1][j+m-1]);
                nMax = MIN (nMax, OutTop3[i+l+1][j+m  ]);
                nMax = MIN (nMax, OutTop3[i+l+1][j+m+1]) - 1;

                nMin =            OutBottom3[i+l-1][j+m-1];
                nMin = MAX (nMin, OutBottom3[i+l-1][j+m  ]);
                nMin = MAX (nMin, OutBottom3[i+l-1][j+m+1]);
                nMin = MAX (nMin, OutBottom3[i+l  ][j+m-1]);
                nMin = MAX (nMin, OutBottom3[i+l  ][j+m  ]);
                nMin = MAX (nMin, OutBottom3[i+l  ][j+m+1]);
                nMin = MAX (nMin, OutBottom3[i+l+1][j+m-1]);
                nMin = MAX (nMin, OutBottom3[i+l+1][j+m  ]);
                nMin = MAX (nMin, OutBottom3[i+l+1][j+m+1]) + 1;

                if (nMin > nMax)
                {
                    DR_Error ("Fix3_Lattice: problem in building the tree (nMin > nMax)!");
                    goto RETURN;
                }


                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                ruL = dev_data->ru + offset;
                r0L = dev_data->r0 + offset;
                rdL = dev_data->rd + offset;

                Shift3L = dev_data->Shift3 + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    Rijk = Rij + (d[5] - Beta33) * k;

                    n = NEAR_INT (Rijk);
                    n = MIN (MAX (nMin - k, n), nMax - k);
                
                    Shift3L[k] = n;

                    Rijk -= n;
                
                    ruL[k] = .5 * (JumpCoeff + Rijk + Rijk * Rijk);
                    rdL[k] = ruL[k] - Rijk;
                    r0L[k] = 1. - ruL[k] - rdL[k];
                }
            }  /* for j */
        }  /* for i */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Lattice */
