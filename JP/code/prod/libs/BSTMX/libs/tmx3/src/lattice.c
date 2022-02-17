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
#include <string.h>
#include <math.h>
#include "tmx123head.h"



/*****  Lattice  **************************************************************/
/*
*       Position the nodes in the lattice: calculate the numeraire prices
*       and the probabilities.
*/
int    Lattice (DEV_DATA        *dev_data,     /* (O) Dev data structure     */
                int             t,             /* (I) Current time point     */
                int             T,             /* (I) Last time point        */
                MKTVOL_DATA     *mktvol_data,  /* (I) Volatility data        */
                TREE_DATA       *tree_data,    /* (I) Tree data structure    */
                T_CURVE         *t_curve)      /* (I) T_curve structure      */
{

    double  *puL,  *pdL,  *p0L;
//    double  *quuL, *qu0L, *qudL;
//    double  *q0uL, *q00L, *q0dL;
//    double  *qduL, *qd0L, *qddL;
//    double  *ruL,  *r0L,  *rdL;

//    double  pu, pd, p0;                 /* Local probabilities */
//    double  qu, qd, q0;

//    double  Grid;                       /* Grid points                      */
    double  ZRatio1, ZRatio2;           /* Zero coupon ratios               */
    double  FwdRate;                    /* Fwd rate for current time point  */
    double  Length;                     /* Length of time steps             */
    double  LengthJ;                    /* Length of time steps for jump    */
    double  ZCenter;                    /* Center of the tree in X-space    */
//    double  Zidx;                       /* Zt index i,j,k adjusted          */
    double  Jump[6];                    /* Jump sizes at current time step  */
    double  PreviousJump[6];            /* Jump sizes at previous time step */
    double  JumpCoeff, d[6], du;        /* Jump coefficients                */
    double  Bbq;                        /* Backbone parameter               */
    double  InsRate;                    /* Instantaneous rate               */
    double  Beta11, Beta21, Beta22; 
    double  Beta31, Beta32, Beta33;     /* Modified mean reversion          */
    double  Pi;                         /* Total drift                      */

    int     *Shift1L;

    int     Top1, Bottom1;              /* Tree limits (1rst dim)           */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)            */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)            */
    int     OutTop1, OutBottom1;        /* Outer tree limits                */
    int     *OutTop2, *OutBottom2;
    int     **OutTop3, **OutBottom3;

    int     CvDiff, CvIdx1, CvIdx2;     /* Nb of zero curves                */

    int     i;                          /* Node indices                     */
    int     offset;                     /* Node offset                      */
    int     l, lMin, lMax;              /* Node branching shifts            */
//    int     m, mMin, mMax;
//    int     n, nMin, nMax;

    long    CurrentDate;
    int     NmrIdx;

    int     status = FAILURE;           /* Error status = FAILURE initially*/


    CurrentDate = tree_data->TPDate[t];

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /*
     *   Internal assigments of zero curves.
     */
    CvDiff = tree_data->CvDiff;
    CvIdx1 = tree_data->CvIdx1;
    CvIdx2 = tree_data->CvIdx2;

    /* Zero ratios for other curves */
    ZRatio1 = tree_data->TermZero[CvDiff][t] / tree_data->TermZero[CvIdx1][t];
    ZRatio2 = tree_data->TermZero[CvDiff][t] / tree_data->TermZero[CvIdx2][t];

    /* There is no future at the end of the tree */
    if (t == T)                                                        
    {
        NmrIdx = tree_data->NbNmr-1;

        if (dev_data->NmrToCcy)
        {
            if (Copy_Slice (dev_data->NmrInv[CvDiff],
                            tree_data->NmrInv[NmrIdx],
                            t,
                            tree_data) == FAILURE) goto RETURN;

            if (SliceTimesScalar (dev_data->NmrInv[CvIdx1],
                                  dev_data->NmrInv[CvDiff],
                                  ZRatio1,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;

            if (SliceTimesScalar (dev_data->NmrInv[CvIdx2],
                                  dev_data->NmrInv[CvDiff],
                                  ZRatio2,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;
        }
      

        return (SUCCESS);
    }

    /* Store lagged numeraire inverse for slice conversion to 
     * numeraire denominated quantities; used lagged slices */
    if (dev_data->CcyToNmr)
    {
        if (Copy_Slice (dev_data->NmrInvLag[CvDiff],
                        dev_data->NmrInv[CvDiff],
                        t+1,
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }

        if (Copy_Slice (dev_data->NmrInvLag[CvIdx1],
                        dev_data->NmrInv[CvIdx1],
                        t+1,
                        tree_data) == FAILURE) 
        {
            goto RETURN;
        }

        if (Copy_Slice (dev_data->NmrInvLag[CvIdx2],
                        dev_data->NmrInv[CvIdx2],
                        t+1,
                        tree_data) == FAILURE) 
        {
            goto RETURN;
        }
    }

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
    
    Beta11 = mktvol_data->Beta[0] * Length * PreviousJump[0] / Jump[0];
    Beta21 = mktvol_data->Beta[1] * Length * PreviousJump[1] / Jump[2];
    Beta22 = mktvol_data->Beta[1] * Length * PreviousJump[2] / Jump[2];
    Beta31 = mktvol_data->Beta[2] * Length * PreviousJump[3] / Jump[5];
    Beta32 = mktvol_data->Beta[2] * Length * PreviousJump[4] / Jump[5];
    Beta33 = mktvol_data->Beta[2] * Length * PreviousJump[5] / Jump[5];
    
        
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

    /* Mapping sigma */
    Bbq      = mktvol_data->Bbq;
    InsRate  = FwdRate / Length;
    
    /* Calculate probabilities */
    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        puL = dev_data->pu + offset;
        p0L = dev_data->p0 + offset;
        pdL = dev_data->pd + offset;

        Shift1L = dev_data->Shift1 + offset;

        lMax = OutTop1 - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Lattice: problem in building the tree (lMin > lMax)!");
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
        return FAILURE;
    }
    else if (tree_data->NbFactor == 3)
    {
        return FAILURE;
    }  /* if then else */   


    /* If numeraire interpolation is inactive, we are DONE (i.e. in Nmr_Calc) */
    if (!tree_data->NmrInterpOn)
    {
        return (SUCCESS);
    }

    /* Get numeraire index: search for first NmrDate >= CurrentDate */
    NmrIdx = GetDLOffset (mktvol_data->NbNmr,
                          mktvol_data->NmrDate,
                          CurrentDate,
                          CbkHIGHER);

    if (NmrIdx < 0 || NmrIdx >= tree_data->NbNmr-1)
    {
        if (dev_data->NmrToCcy)
        {
            DR_Error ("Interpolated numeraire can only be used before "
                      "%ld (current date: %ld)\n", 
                      mktvol_data->NmrDate[mktvol_data->NbNmr-2],
                      CurrentDate);
            goto RETURN;
        }
        else /* nothing more to do */
        {
            return (SUCCESS);
        }
    }

    /* Special case: on numeraire dates (includes ValueDate for NmrIdx==0) */
    if (CurrentDate == mktvol_data->NmrDate[NmrIdx])
    {
        /* Re-Initialise Next Zero for Nmr Interpolation     */
        if (Copy_Slice (tree_data->LastZero,
                        tree_data->NmrInv[NmrIdx],
                        t,
                        tree_data) == FAILURE) goto RETURN;

        /* Re-Initialise Next Critical date for Nmr Interpolation     */
        tree_data->NxtCritDate = mktvol_data->NmrDate[NmrIdx];



        if (dev_data->NmrToCcy)
        {
            if (Copy_Slice (dev_data->NmrInv[CvDiff],
                            tree_data->NmrInv[NmrIdx],
                            t,
                            tree_data) == FAILURE) goto RETURN;

            if (SliceTimesScalar (dev_data->NmrInv[CvIdx1],
                                  dev_data->NmrInv[CvDiff],
                                  ZRatio1,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;

            if (SliceTimesScalar (dev_data->NmrInv[CvIdx2],
                                  dev_data->NmrInv[CvDiff],
                                  ZRatio2,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;
        }

        return (SUCCESS);
    }
    else /* In-between numeraire dates */
    {
        /* Roll back Zero slice (1^) Using Ev, not Dev */
        if (Ev (tree_data->LastZero,t,T,dev_data,tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Interpolate numeraire if required */
        if (dev_data->NmrToCcy)
        {  
            if (Nmr_Interp (dev_data->NmrInv[CvDiff],
                            NmrIdx,
                            t,
                            mktvol_data,
                            tree_data,
                            t_curve) == FAILURE) goto RETURN;
            
            
            /* Re-Initialise Next Zero for Nmr Interpolation     */
            if (Copy_Slice (tree_data->LastZero,
                            dev_data->NmrInv[CvDiff],
                            t,
                            tree_data) == FAILURE) goto RETURN;
            
            /* Re-Initialise Next Critical date for Nmr Interpolation     */
            tree_data->NxtCritDate = tree_data->TPDate[t];



            if (SliceTimesScalar (dev_data->NmrInv[CvIdx1],
                                  dev_data->NmrInv[CvDiff],
                                  ZRatio1,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;

            if (SliceTimesScalar (dev_data->NmrInv[CvIdx2],
                                  dev_data->NmrInv[CvDiff],
                                  ZRatio2,
                                  t,
                                  tree_data) == FAILURE) goto RETURN;

        } /* if numeraire required */
 
    } /* if not on NmrDate */
  

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
    {    
        DR_Error("Lattice: Failed!");
    }

    return (status);

}  /* Lattice */


