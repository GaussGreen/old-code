/**************************************************************************/
/*      Payoff function for the sticky/adjustable collar                  */
/**************************************************************************/
/*      sticky.c                                                          */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"

/*****  Fix3_Sticky_t  *********************************************************/
/*
*   Calculates the sticky/adjustable price
*   dev's the sticky to this time point (t)
*   if FltResetFlag=TRUE, then add payoff and revise state levels
*   if *PrevStates is NULL, then assume Sticky is uninitialised    
*/
int  Fix3_Sticky_t(
                
        double    **Sticky,                 /* (I/O) Prices for all states */
        long        FltResetFlag,           /* (I)   Reset at this timept  */
        double     *FltIndex,               /* (I)   Index level           */
        double     *ZeroToPmt,              /* (I)   PmtZero 4 this reset  */
        double      Notional,               /* (I)                         */
        double      DCF,                    /* (I)                         */
        double      FltSpread,              /* (I)                         */
        double      StickyFltSpread,        /* (I)                         */
        double      CapSpread,              /* (I)                         */
        double      StickyCapSpread,        /* (I)                         */
        double      FloorSpread,            /* (I)                         */
        double      StickyFloorSpread,      /*(I)                          */
        double      LifeCapLevel,           /* (I)                         */
        double      StickyLifeCapLevel,     /*(I)                          */
        double      LifeFloorLevel,         /* (I)                         */
        double      StickyLifeFloorLevel,   /* (I)                         */
        int         IsEmbedded,             /* (I)   Embedded/Option only  */
        int         IsSimple,               /* (I)   Simple/Compound pmt   */
        int         IsSticky,               /* (I)   sticky/adjustable     */
        int         IsStickyOff,            /* (I)   use cap/flr spread?   */
        int         NbStates,               /* (I)   Nb of states          */
        double      MaxState,               /* (I)   Upper state bound     */
        double      MinState,               /* (I)   Lower state bound     */
        double    **PrevStates,             /* (I)   Prev state levels     */
        int         t,                      /* (I)   Current time point    */
        int         T,                      /* (I)   Last time point       */
        int         DCurve,                 /* (I)   Discount curve        */
        FIX3_DEV_DATA    *dev_data,              /* (I)   Fix3_Dev data structure    */
        FIX3_TREE_DATA   *tree_data)             /* (I)   Tree data structure   */
{
    /* Local slice pointers */

    double  *StickyL[MAXNBSTATES+1];
    double  *FltIndexL;
    double  *ZeroToPmtL;

    /* state-variable variables */

    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs    */
    int     s;                        /* State variable index            */
    int     sj;                       /* Interp state variable index     */
    int     q;                        /* Index for quad, linear intrp    */
    double  *State  = NULL;           /* Current levels of state var     */
    double  *PState = NULL;           /* Previous levels of state var    */
    double  deltaS;                   /* Increament between cons states  */
    double  InterpStateLevel;         /* state level for interp          */

    /* payoff variables */

    double  Payoff[MAXNBSTATES];      /* Intermediate values of sticky   */
    double  tmpPayoff;
    double  CapStrike[MAXNBSTATES];
    double  StickyCapStrike[MAXNBSTATES];
    double  FloorStrike[MAXNBSTATES];
    double  StickyFloorStrike[MAXNBSTATES];
    int     IsLastReset;               /* true=this is last reset date    */
    double  Q, CouponRate;
    double  QSticky, CouponRateSticky; /* variables used to determine     */
                                       /* state level for interp          */
    
    /* Tree variables */

    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    IsLastReset = (*PrevStates == NULL);

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Fix3_Sticky_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* Discount all the state variables */
    
    if (!IsLastReset)
    {
        for (s=0; s<NbStates; s++)
        {
            if (Fix3_Dev(Sticky[s],
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
            {
                goto RETURN;
            }
        }
    }

    /* Add payoff if it's a reset */
    
    if (FltResetFlag)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /*   and the cap/floor strikes              */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("Fix3_Sticky_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Fix3_Sticky_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s = 0; s < NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        for (s = 0; s < NbStates; s++)
        {
            CapStrike[s] = State[s] + CapSpread;
            FloorStrike[s] = State[s] + FloorSpread;
            StickyCapStrike[s] = State[s] + StickyCapSpread;
            StickyFloorStrike[s] = State[s] + StickyFloorSpread;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */


        if (!IsLastReset)
        {
            PState = *PrevStates;

            if (!IS_EQUAL(PState[0], PState[NbStates-1]))
            {
                for (s = 0; s < NbStates-2; s++)
                {
                    D[s][0] = 1. / ((PState[s]-PState[s+1])
                                   *(PState[s]-PState[s+2]));
                    D[s][1] = 1. / ((PState[s+1]-PState[s])
                                   *(PState[s+1]-PState[s+2]));
                    D[s][2] = 1. / ((PState[s+2]-PState[s])
                                   *(PState[s+2]-PState[s+1]));
                }
            }
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            FltIndexL  = FltIndex  + offset;
            ZeroToPmtL = ZeroToPmt + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {
                for (s = 0; s < NbStates; s++)
                {
                    /* ---------------------------- */
                    /*   Add the new coupon amt     */
                    /* ---------------------------- */
        
                    Q = (IsStickyOff) ?                                 \
                         FltIndexL[i] + FltSpread :                     \
                         COLLAR(FltIndexL[i] + FltSpread, CapStrike[s], \
                                                        FloorStrike[s]); 


                    QSticky = (IsStickyOff) ?
                               FltIndexL[i] + StickyFltSpread :
                               COLLAR (FltIndexL[i] + StickyFltSpread,
                                                        StickyCapStrike[s],
                                                        StickyFloorStrike[s]);

                    CouponRate = COLLAR(Q, LifeCapLevel, LifeFloorLevel);
                    CouponRateSticky = COLLAR (QSticky, StickyLifeCapLevel,
                                                        StickyLifeFloorLevel);

                    if (IsEmbedded)
                    {
                        Payoff[s] = Notional * ZeroToPmtL[i] * 
                                    ACC_FN(CouponRate,DCF,IsSimple);
                    }
                    else
                    {
                        Payoff[s] = Notional * ZeroToPmtL[i] * 
                                    ACC_FN_DIFF(FltIndexL[i]+FltSpread,
                                                CouponRate,
                                                DCF,
                                                IsSimple);
                    }

                    /* ------------------------------------- */
                    /*   Add the interpolated Sticky price   */
                    /*   to the payoff if appropriate        */
                    /* ------------------------------------- */

                    if (!IsLastReset)
                    {
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] += StickyL[0][i];
                        }
                        else
                        {
                            InterpStateLevel = (IsSticky) ? CouponRateSticky :
                                                            FltIndexL[i];

                            sj = (int) floor((InterpStateLevel - PState[0]) /
                                             (PState[1] - PState[0]));

                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] += StickyL[NbStates-1][i];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] += StickyL[0][i];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                    
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         StickyL[q][i],StickyL[q+1][i],
                                         StickyL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] += tmpPayoff;
                               
                            } /* if sj */
                        }
                    } /* if !IsLastReset */

                } /* for state idx */

                /* replace sticky prices with new values for each state */

                for (s=0; s<NbStates ; s++)
                {
                    StickyL[s][i] = Payoff[s];
                }

            } /* for i */

        } /* if NbFactor == 1 */

        /************************  2 FACTORS   ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                FltIndexL  = FltIndex  + offset;
                ZeroToPmtL = ZeroToPmt + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ---------------------------- */
                        /*   Add the new coupon amt     */
                        /* ---------------------------- */
            
                        Q = (IsStickyOff) ?                                  \
                             FltIndexL[j]+FltSpread :                        \
                             COLLAR(FltIndexL[j]+FltSpread, CapStrike[s],    \
                                                            FloorStrike[s]);
                        



                        QSticky = (IsStickyOff) ?                           \
                                   FltIndexL[j] + StickyFltSpread :         \
                                   COLLAR(FltIndexL[j] + StickyFltSpread,   \
                                                        StickyCapStrike[s], \
                                                        StickyFloorStrike[s]);


                        CouponRate = COLLAR(Q, LifeCapLevel, LifeFloorLevel);
                        CouponRateSticky = COLLAR (QSticky,
                                                   StickyLifeCapLevel,
                                                   StickyLifeFloorLevel);

                        if (IsEmbedded)
                        {
                            Payoff[s] = Notional * ZeroToPmtL[j] * 
                                        ACC_FN(CouponRate,DCF,IsSimple);
                        }
                        else
                        {
                            Payoff[s] = Notional * ZeroToPmtL[j] * 
                                        ACC_FN_DIFF(FltIndexL[j]+FltSpread,
                                                    CouponRate,
                                                    DCF,
                                                    IsSimple);
                        }
    
                        /* ------------------------------------- */
                        /*   Add the interpolated Sticky price   */
                        /*   to the payoff if appropriate        */
                        /* ------------------------------------- */
    
                        if (!IsLastReset)
                        {
                            if (IS_EQUAL(PState[0], PState[1]))
                            {
                                Payoff[s] += StickyL[0][j];
                            }
                            else
                            {
                                InterpStateLevel = (IsSticky) ?         \
                                                    CouponRateSticky :  \
                                                    FltIndexL[j];

                                sj = (int)floor((InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
    
                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] += StickyL[NbStates-1][j];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += StickyL[0][j];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q], PState[q+1], 
                                             PState[q+2],
                                             StickyL[q][j],StickyL[q+1][j],
                                             StickyL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] += tmpPayoff;
                                   
                                } /* if sj */
                            }
                        } /* if !IsLastReset */
    
                    } /* for state idx */
    
                    /* replace sticky prices with new values for each state */
    
                    for (s=0; s<NbStates ; s++)
                    {
                        StickyL[s][j] = Payoff[s];
                    }
                                    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTORS   ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }

                    FltIndexL  = FltIndex  + offset;
                    ZeroToPmtL = ZeroToPmt + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ---------------------------- */
                            /*   Add the new coupon amt     */
                            /* ---------------------------- */
                
                            Q = (IsStickyOff) ?                         \
                                 FltIndexL[k]+FltSpread :               \
                                 COLLAR(FltIndexL[k]+FltSpread,         \
                                        CapStrike[s],                   \
                                        FloorStrike[s]); 


                            QSticky = (IsStickyOff) ?                        \
                                        FltIndexL[k]+StickyFltSpread :       \
                                        COLLAR(FltIndexL[k]+StickyFltSpread, \
                                        StickyCapStrike[s],                  \
                                        StickyFloorStrike[s]);

        
                            CouponRate = COLLAR(Q, LifeCapLevel, 
                                                LifeFloorLevel);

                            CouponRateSticky = COLLAR (QSticky,             \
                                                       StickyLifeCapLevel,  \
                                                       StickyLifeFloorLevel);
        
                            if (IsEmbedded)
                            {
                                Payoff[s] = Notional * ZeroToPmtL[k] * 
                                            ACC_FN(CouponRate,DCF,IsSimple);
                            }
                            else
                            {
                                Payoff[s] = Notional * ZeroToPmtL[k] * 
                                            ACC_FN_DIFF(FltIndexL[k]+FltSpread,
                                                        CouponRate,
                                                        DCF,
                                                        IsSimple);
                            }
        
                            /* ------------------------------------- */
                            /*   Add the interpolated Sticky price   */
                            /*   to the payoff if appropriate        */
                            /* ------------------------------------- */
        
                            if (!IsLastReset)
                            {
                                if (IS_EQUAL(PState[0], PState[1]))
                                {
                                    Payoff[s] += StickyL[0][k];
                                }
                                else
                                {
                                    InterpStateLevel = (IsSticky) ?           \
                                                          CouponRateSticky :  \
                                                          FltIndexL[k];
                                    sj = (int) floor(
                                                (InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
        
                                    if (sj >= (NbStates - 1)) /* right */
                                    {
                                        Payoff[s] += StickyL[NbStates-1][k];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] += StickyL[0][k];
                                    }
                                    else /* in between */
                                    {
                                        q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(PState[q], PState[q+1], 
                                                 PState[q+2],
                                                 StickyL[q][k],StickyL[q+1][k],
                                                 StickyL[q+2][k],
                                                 D[q][0], D[q][1], D[q][2],
                                                 InterpStateLevel,
                                                 &(tmpPayoff));
                                        Payoff[s] += tmpPayoff;
                                       
                                    } /* if sj */
                                }
                            } /* if !IsLastReset */
        
                        } /* for state idx */
        
                        /* replace sticky prices for each state */
        
                        for (s=0; s<NbStates ; s++)  
                        {
                            StickyL[s][k] = Payoff[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

        /* -------------------------------------- */
        /*   replace PrevStates with new states   */
        /* -------------------------------------- */

        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;

    }  /* if FltResetFlag */

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }

    return (status);

}  /* Fix3_Sticky_t */




/*****  Fix3_StickySwap_t  *****************************************************/
/*
*   Calculates the sticky/adjustable swap price
*   no dev's, just add sticky and floating coupons
*/
int  Fix3_StickySwap_t
             (double    **Sticky,         /* (O) Prices for all states   */
              long        ResetFlagSticky,/* (I) Reset at timept         */
              long        ResetFlagFloat, /* (I) Reset at timept         */
              double     *IndexSticky,    /* (I) Index level             */
              double     *IndexFloat,     /* (I) Index level             */
              double     *ZeroToPmtSticky,/* (I) PmtZero 4 this reset    */
              double     *ZeroToPmtFloat, /* (I) PmtZero 4 this reset    */
              double      WeightSticky,   /* (I) Weight of sticky index  */
              double      FloorSticky,    /* (I) Floor spread            */
              double      CapSticky,      /* (I) Cap spread              */
              double      SpreadSticky,   /* (I) Sticky spread           */
              double      OutsSticky,     /* (I) Outstanding             */
              double      DcfSticky,      /* (I) Day count fract         */
              double      LifeCapLevel,   /* (I)                         */
              double      LifeFloorLevel, /* (I)                         */
              double      FloorFloat,     /* (I) Floor rate              */
              double      CapFloat,       /* (I) Cap rate                */
              double      StepUpFloat,    /* (I) Float stepup            */
              double      OutsFloat,      /* (I) Outstanding             */
              double      DcfFloat,       /* (I) Day count fract         */
              int         IsSticky,       /* (I) Sticky/adjustable       */
              int         IsStickyOff,    /* (I) Use cap/flr spread?     */
              int         IsSimpleSticky, /* (I) Simple/Compound pmt     */
              int         IsSimpleFloat,  /* (I) Simple/Compound pmt     */
              int         NbStates,       /* (I) Nb of states            */
              double      MaxState,       /* (I) Upper state bound       */
              double      MinState,       /* (I) Lower state bound       */
              double    **PrevStates,     /* (I) Prev state levels       */
              int         t,              /* (I) Current time point      */
              FIX3_TREE_DATA   *tree_data)     /* (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *StickyL[MAXNBSTATES+1];
    double  *IndexStickyL;
    double  *IndexFloatL;
    double  *ZeroToPmtStickyL;
    double  *ZeroToPmtFloatL;

    /* state-variable variables */

    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs    */
    int     s;                        /* State variable index            */
    int     sj;                       /* Interp state variable index     */
    int     q;                        /* Index for quad, linear intrp    */
    double  *State  = NULL;           /* Current levels of state var     */
    double  *PState = NULL;           /* Previous levels of state var    */
    double  deltaS;                   /* Increament between cons states  */
    double  InterpStateLevel;         /* state level for interp          */

    /* payoff variables */

    double  Payoff[MAXNBSTATES];      /* Intermediate values of sticky   */
    double  tmpPayoff;
    double  CapStrike[MAXNBSTATES];
    double  FloorStrike[MAXNBSTATES];
    int     IsLastReset;              /* true=this is last reset date    */
    double  Q, CouponRate, CouponAmt;
    
    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    IsLastReset = (*PrevStates == NULL);

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Fix3_Sticky_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* FLOATIN LEG */

    /* Add floating payoff if it's a reset */
    if (ResetFlagFloat)
    {

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            IndexFloatL  = IndexFloat  + offset;
            ZeroToPmtFloatL = ZeroToPmtFloat + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                CouponRate = COLLAR(IndexFloatL[i]+StepUpFloat,  \
                                    CapFloat,                    \
                                    FloorFloat); 

                CouponAmt = OutsFloat * ZeroToPmtFloatL[i] * 
                            ACC_FN(CouponRate,DcfFloat,IsSimpleFloat);

                for (s = 0; s < NbStates ; s++)  
                {
                    StickyL[s][i] -= CouponAmt;
                }

            } /* for i */

        } 

        /************************  2 FACTOR    ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                IndexFloatL  = IndexFloat  + offset;
                ZeroToPmtFloatL = ZeroToPmtFloat + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    CouponRate = COLLAR(IndexFloatL[j]+StepUpFloat,  \
                                        CapFloat,                    \
                                        FloorFloat); 
                
                    CouponAmt = OutsFloat * ZeroToPmtFloatL[j] * 
                                ACC_FN(CouponRate,DcfFloat,IsSimpleFloat);

                    for (s = 0; s < NbStates ; s++)  
                    {
                        StickyL[s][j] -= CouponAmt;
                    }
                                    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTOR    ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }

                    IndexFloatL  = IndexFloat  + offset;
                    ZeroToPmtFloatL = ZeroToPmtFloat + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        CouponRate = COLLAR(IndexFloatL[k]+StepUpFloat,  \
                                            CapFloat,                    \
                                            FloorFloat); 
                
                        CouponAmt = OutsFloat * ZeroToPmtFloatL[k] * 
                                    ACC_FN(CouponRate,DcfFloat,IsSimpleFloat);

                        for (s=0; s<NbStates ; s++)  
                        {
                            StickyL[s][k] -= CouponAmt;
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlagFloat */

    /* STICKY LEG */

    /* Add sticky payoff if it's a reset */
    if (ResetFlagSticky)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /*   and the cap/floor strikes              */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("Fix3_Sticky_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Fix3_Sticky_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        for (s=0; s<NbStates; s++)
        {
            CapStrike[s] = State[s] + CapSticky;
            FloorStrike[s] = State[s] + FloorSticky;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */


        if (!IsLastReset)
        {
            PState = *PrevStates;

            if (!IS_EQUAL(PState[0], PState[NbStates-1]))
            {
                for (s = 0; s < NbStates-2; s++)
                {
                    D[s][0] = 1. / ((PState[s]-PState[s+1])
                                   *(PState[s]-PState[s+2]));
                    D[s][1] = 1. / ((PState[s+1]-PState[s])
                                   *(PState[s+1]-PState[s+2]));
                    D[s][2] = 1. / ((PState[s+2]-PState[s])
                                   *(PState[s+2]-PState[s+1]));
                }
            }
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            IndexStickyL  = IndexSticky  + offset;
            ZeroToPmtStickyL = ZeroToPmtSticky + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates; s++)
                {
                    /* ---------------------------- */
                    /*   Add the new coupon amt     */
                    /* ---------------------------- */
        
                    Q = (IsStickyOff) ?                                       \
                         IndexStickyL[i]*WeightSticky+SpreadSticky :          \
                         COLLAR(IndexStickyL[i]*WeightSticky+SpreadSticky,    \
                                CapStrike[s],                                 \
                                FloorStrike[s]); 

                    CouponRate = COLLAR(Q, LifeCapLevel, LifeFloorLevel);

                    Payoff[s] = OutsSticky * ZeroToPmtStickyL[i] * 
                                ACC_FN(CouponRate,DcfSticky,IsSimpleSticky);

                    /* ------------------------------------- */
                    /*   Add the interpolated Sticky price   */
                    /*   to the payoff if appropriate        */
                    /* ------------------------------------- */

                    if (!IsLastReset)
                    {
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] += StickyL[0][i];
                        }
                        else
                        {
                            InterpStateLevel = (IsSticky) ? CouponRate :  \
                                                            IndexStickyL[i];
                            sj = (int) floor((InterpStateLevel - PState[0]) /
                                             (PState[1] - PState[0]));

                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] += StickyL[NbStates-1][i];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] += StickyL[0][i];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                    
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         StickyL[q][i],StickyL[q+1][i],
                                         StickyL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] += tmpPayoff;
                               
                            } /* if sj */
                        }
                    } 
                    else
                    {
                        Payoff[s] += StickyL[0][i];
                    } /* if !IsLastReset */

                } /* for state idx */

                /* replace sticky prices with new values for each state */

                for (s=0; s<NbStates ; s++)  
                {
                    StickyL[s][i] = Payoff[s];
                }

            } /* for i */

        } /* if NbFactor == 1 */

        /************************  2 FACTORS   ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                IndexStickyL  = IndexSticky  + offset;
                ZeroToPmtStickyL = ZeroToPmtSticky + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ---------------------------- */
                        /*   Add the new coupon amt     */
                        /* ---------------------------- */
            
                        Q = (IsStickyOff) ?                                   \
                             IndexStickyL[j]*WeightSticky+SpreadSticky :      \
                             COLLAR(IndexStickyL[j]*WeightSticky+SpreadSticky,\
                                    CapStrike[s],                             \
                                    FloorStrike[s]); 
    
                        CouponRate = COLLAR(Q, LifeCapLevel, LifeFloorLevel);
    
                        Payoff[s] = OutsSticky * ZeroToPmtStickyL[j] * 
                                    ACC_FN(CouponRate,DcfSticky,IsSimpleSticky);
    
                        /* ------------------------------------- */
                        /*   Add the interpolated Sticky price   */
                        /*   to the payoff if appropriate        */
                        /* ------------------------------------- */
    
                        if (!IsLastReset)
                        {
                            if (IS_EQUAL(PState[0], PState[1]))
                            {
                                Payoff[s] += StickyL[0][j];
                            }
                            else
                            {
                                InterpStateLevel = (IsSticky) ? CouponRate : \
                                                                IndexStickyL[j];
                                sj = (int)floor((InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
    
                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] += StickyL[NbStates-1][j];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += StickyL[0][j];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q], PState[q+1], 
                                             PState[q+2],
                                             StickyL[q][j],StickyL[q+1][j],
                                             StickyL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] += tmpPayoff;
                                   
                                } /* if sj */
                            }
                        } 
                        else
                        {
                            Payoff[s] += StickyL[0][j];
                        } 
    
                    } /* for state idx */
    
                    /* replace sticky prices with new values for each state */
    
                    for (s=0; s<NbStates ; s++)  
                    {
                        StickyL[s][j] = Payoff[s];
                    }
                                    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTORS   ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }

                    IndexStickyL  = IndexSticky  + offset;
                    ZeroToPmtStickyL = ZeroToPmtSticky + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ---------------------------- */
                            /*   Add the new coupon amt     */
                            /* ---------------------------- */
                
                            Q = (IsStickyOff) ?                               \
                                 IndexStickyL[k]*WeightSticky+SpreadSticky :  \
                                 COLLAR(IndexStickyL[k]*WeightSticky          \
                                            +SpreadSticky,                    \
                                        CapStrike[s],                         \
                                        FloorStrike[s]); 
        
                            CouponRate = COLLAR(Q, LifeCapLevel, 
                                                LifeFloorLevel);
        
                            Payoff[s] = OutsSticky * ZeroToPmtStickyL[k] * 
                                        ACC_FN(CouponRate,DcfSticky,IsSimpleSticky);
        
                            /* ------------------------------------- */
                            /*   Add the interpolated Sticky price   */
                            /*   to the payoff if appropriate        */
                            /* ------------------------------------- */
        
                            if (!IsLastReset)
                            {
                                if (IS_EQUAL(PState[0], PState[1]))
                                {
                                    Payoff[s] += StickyL[0][k];
                                }
                                else
                                {
                                    InterpStateLevel = (IsSticky) ?          \
                                                          CouponRate :       \
                                                          IndexStickyL[k];
                                    sj = (int) floor(
                                                (InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
        
                                    if (sj >= (NbStates - 1)) /* right */
                                    {
                                        Payoff[s] += StickyL[NbStates-1][k];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] += StickyL[0][k];
                                    }
                                    else /* in between */
                                    {
                                        q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(PState[q], PState[q+1], 
                                                 PState[q+2],
                                                 StickyL[q][k],StickyL[q+1][k],
                                                 StickyL[q+2][k],
                                                 D[q][0], D[q][1], D[q][2],
                                                 InterpStateLevel,
                                                 &(tmpPayoff));
                                        Payoff[s] += tmpPayoff;
                                       
                                    } /* if sj */
                                }
                            } 
                            else
                            {
                                Payoff[s] += StickyL[0][k];
                            } 
        
                        } /* for state idx */
        
                        /* replace sticky prices for each state */
        
                        for (s=0; s<NbStates ; s++)  
                        {
                            StickyL[s][k] = Payoff[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

        /* -------------------------------------- */
        /*   replace PrevStates with new states   */
        /* -------------------------------------- */

        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;

    }  /* if ResetFlagSticky */

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }

    return (status);

}  /* Fix3_StickySwap_t */


/*****  Fix3_VolBond_t  *****************************************************/
/*
 *   Calculates the vol bond price - NO FUNDING LEG
 *   no dev's, just add complex coupons
 */
int  Fix3_VolBond_t(double    **VolBond,        /* (O) Prices for all states   */
               long        ResetFlagCplx,  /* (I) Reset at timept         */
               double     *IndexCplx,      /* (I) Index level             */
               double     *ZeroToPmtCplx,  /* (I) PmtZero 4 this reset    */
               double      LeverageCplx,   /* (I) Leverage of cplx index  */
               double      FloorCplx,      /* (I) Floor for cplx leg      */
               double      CapCplx,        /* (I) Cap for cplx leg        */
               double      StepUpCplx,     /* (I) StepUp for cplx leg     */
               double      OutsCplx,       /* (I) Outstanding             */
               double      DcfCplx,        /* (I) Day count fract         */
               int         IsCplxPmtOff,   /* (I) TRUE=no cplx leg pmt    */
               int         NbStates,       /* (I) Nb of states            */
               double      MaxState,       /* (I) Upper state bound       */
               double      MinState,       /* (I) Lower state bound       */
               double    **PrevStates,     /* (I) Prev state levels       */
               int         t,              /* (I) Current time point      */
               FIX3_TREE_DATA   *tree_data)     /* (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *VolBondL[MAXNBSTATES+1];
    double  *IndexCplxL;
    double  *ZeroToPmtCplxL;

    /* state-variable variables */

    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs    */
    int     s;                        /* State variable index            */
    int     sj;                       /* Interp state variable index     */
    int     q;                        /* Index for quad, linear intrp    */
    double  *State  = NULL;           /* Current levels of state var     */
    double  *PState = NULL;           /* Previous levels of state var    */
    double  deltaS;                   /* Increament between cons states  */
    double  InterpStateLevel;         /* state level for interp          */

    /* payoff variables */

    double  Payoff[MAXNBSTATES];      /* Intermediate values of volbond  */
    double  tmpPayoff;
    int     IsLastReset;              /* true=this is last reset date    */
    double  CouponRate;

    /* Tree variables */

    int     Top1, Bottom1;            /* Tree limits (1st dim)  */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    IsLastReset = (*PrevStates == NULL);

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Fix3_VolBond_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* COMPLEX LEG */

    /* Add complex payoff if it's a reset */
    if (ResetFlagCplx)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /*   and the cap/floor strikes              */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("Fix3_VolBond_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Fix3_VolBond_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */

        if (!IsLastReset)
        {
            PState = *PrevStates;

            if (!IS_EQUAL(PState[0], PState[NbStates-1]))
            {
                for (s = 0; s < NbStates-2; s++)
                {
                    D[s][0] = 1. / ((PState[s]-PState[s+1])
                                   *(PState[s]-PState[s+2]));
                    D[s][1] = 1. / ((PState[s+1]-PState[s])
                                   *(PState[s+1]-PState[s+2]));
                    D[s][2] = 1. / ((PState[s+2]-PState[s])
                                   *(PState[s+2]-PState[s+1]));
                }
            }
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                VolBondL[s] = VolBond[s] + offset;
            }
            IndexCplxL     = IndexCplx  + offset;
            ZeroToPmtCplxL = ZeroToPmtCplx + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {
                for (s = 0; s < NbStates; s++)
                {
                    /* ---------------------------- */
                    /*   Add the new coupon amt     */
                    /* ---------------------------- */

                    if (IsCplxPmtOff)
                    {
                        Payoff[s] = 0.0;
                    }
                    else
                    {
                        CouponRate = ABS(IndexCplxL[i] - State[s]) * 
                                     LeverageCplx + StepUpCplx;

                        CouponRate = COLLAR(CouponRate, CapCplx, FloorCplx);

                        Payoff[s] = OutsCplx * ZeroToPmtCplxL[i] * 
                                    ACC_FN(CouponRate,DcfCplx,TRUE);
                    }

                    /* ------------------------------------- */
                    /*   Add the interpolated VolBond price  */
                    /*   to the payoff if appropriate        */
                    /* ------------------------------------- */

                    if (!IsLastReset)
                    {
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] += VolBondL[0][i];
                        }
                        else
                        {
                            InterpStateLevel = IndexCplxL[i];

                            sj = (int) floor((InterpStateLevel - PState[0]) /
                                             (PState[1] - PState[0]));

                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] += VolBondL[NbStates-1][i];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] += VolBondL[0][i];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                    
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         VolBondL[q][i],VolBondL[q+1][i],
                                         VolBondL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] += tmpPayoff;
                            } /* if sj */
                        }
                    } 
                    else
                    {
                        Payoff[s] += VolBondL[0][i];
                    } /* if !IsLastReset */
                } /* for state idx */

                /* replace volbond prices with new values for each state */
                for (s=0; s<NbStates ; s++)
                {
                    VolBondL[s][i] = Payoff[s];
                }

            } /* for i */

        } /* if NbFactor == 1 */

        /************************  2 FACTORS   ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    VolBondL[s] = VolBond[s] + offset;
                }
                IndexCplxL     = IndexCplx  + offset;
                ZeroToPmtCplxL = ZeroToPmtCplx + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ---------------------------- */
                        /*   Add the new coupon amt     */
                        /* ---------------------------- */

                        if (IsCplxPmtOff)
                        {
                            Payoff[s] = 0.0;
                        }
                        else
                        {
                            CouponRate = ABS(IndexCplxL[j] - State[s]) * 
                                         LeverageCplx + StepUpCplx;

                            CouponRate = COLLAR(CouponRate, 
                                                CapCplx, 
                                                FloorCplx);

                            Payoff[s] = OutsCplx * ZeroToPmtCplxL[j] * 
                                        ACC_FN(CouponRate,DcfCplx,TRUE);
                        }

                        /* ------------------------------------- */
                        /*   Add the interpolated VolBond price   */
                        /*   to the payoff if appropriate        */
                        /* ------------------------------------- */
    
                        if (!IsLastReset)
                        {
                            if (IS_EQUAL(PState[0], PState[1])) 
                            {
                                Payoff[s] += VolBondL[0][j];
                            }
                            else
                            {
                                InterpStateLevel = IndexCplxL[j];

                                sj = (int)floor((InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
    
                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] += VolBondL[NbStates-1][j];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += VolBondL[0][j];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q], PState[q+1], 
                                             PState[q+2],
                                             VolBondL[q][j],VolBondL[q+1][j],
                                             VolBondL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] += tmpPayoff;
                                   
                                } /* if sj */
                            }
                        } 
                        else
                        {
                            Payoff[s] += VolBondL[0][j];
                        } 
                    } /* for state idx */

                    /* replace volbond prices with new values for each state */
                    for (s=0; s<NbStates ; s++)  
                    {
                        VolBondL[s][j] = Payoff[s];
                    }

                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTORS   ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        VolBondL[s] = VolBond[s] + offset;
                    }

                    IndexCplxL  = IndexCplx  + offset;
                    ZeroToPmtCplxL = ZeroToPmtCplx + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ---------------------------- */
                            /*   Add the new coupon amt     */
                            /* ---------------------------- */

                            if (IsCplxPmtOff)
                            {
                                Payoff[s] = 0.0;
                            }
                            else
                            {
                                CouponRate = ABS(IndexCplxL[k] - State[s]) * 
                                             LeverageCplx + StepUpCplx;

                                CouponRate = COLLAR(CouponRate, 
                                                    CapCplx, 
                                                    FloorCplx);

                                Payoff[s] = OutsCplx * ZeroToPmtCplxL[k] * 
                                            ACC_FN(CouponRate,DcfCplx,TRUE);
                            }

                            /* ------------------------------------- */
                            /*   Add the interpolated VolBond price  */
                            /*   to the payoff if appropriate        */
                            /* ------------------------------------- */
        
                            if (!IsLastReset)
                            {
                                if (IS_EQUAL(PState[0], PState[1])) 
                                {
                                    Payoff[s] += VolBondL[0][k];
                                }
                                else
                                {
                                    InterpStateLevel = IndexCplxL[k];

                                    sj = (int) floor(
                                                (InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
        
                                    if (sj >= (NbStates - 1)) /* right */
                                    {
                                        Payoff[s] += VolBondL[NbStates-1][k];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] += VolBondL[0][k];
                                    }
                                    else /* in between */
                                    {
                                        q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(
                                              PState[q], PState[q+1], 
                                              PState[q+2],
                                              VolBondL[q][k],VolBondL[q+1][k],
                                              VolBondL[q+2][k],
                                              D[q][0], D[q][1], D[q][2],
                                              InterpStateLevel,
                                              &(tmpPayoff));
                                        Payoff[s] += tmpPayoff;
                                       
                                    } /* if sj */
                                }
                            } 
                            else
                            {
                                Payoff[s] += VolBondL[0][k];
                            }
                        } /* for state idx */

                        /* replace volbond prices for each state */
                        for (s=0; s<NbStates ; s++)  
                        {
                            VolBondL[s][k] = Payoff[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

        /* -------------------------------------- */
        /*   replace PrevStates with new states   */
        /* -------------------------------------- */

        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;

    }  /* if ResetFlagCplx */

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }

    return (status);

}  /* Fix3_VolBond_t */



/*****  Fix3_VolSwap_t  *****************************************************/
/*
 *   Calculates the vol swap price
 *   no dev's, just add complex and funding coupons
 */
int  Fix3_VolSwap_t(double    **VolBond,        /* (O) Prices for all states   */
               long        ResetFlagCplx,  /* (I) Reset at timept         */
               long        ResetFlagFund,  /* (I) Reset at timept         */
               double     *IndexCplx,      /* (I) Index level             */
               double     *IndexFund,      /* (I) Index level             */
               double     *ZeroToPmtCplx,  /* (I) PmtZero 4 this reset    */
               double     *ZeroToPmtFund,  /* (I) PmtZero 4 this reset    */
               double      LeverageCplx,   /* (I) Leverage of cplx index  */
               double      FloorCplx,      /* (I) Floor for cplx leg      */
               double      CapCplx,        /* (I) Cap for cplx leg        */
               double      StepUpCplx,     /* (I) StepUp for cplx leg     */
               double      OutsCplx,       /* (I) Outstanding             */
               double      DcfCplx,        /* (I) Day count fract         */
               double      FloorFund,      /* (I) Floor for fund leg      */
               double      CapFund,        /* (I) Cap for fund leg        */
               double      StepUpFund,     /* (I) Stepup for funding leg  */
               double      OutsFund,       /* (I) Outstanding             */
               double      DcfFund,        /* (I) Day count fract         */
               int         IsCplxPmtOff,   /* (I) TRUE=no cplx leg pmt    */
               int         NbStates,       /* (I) Nb of states            */
               double      MaxState,       /* (I) Upper state bound       */
               double      MinState,       /* (I) Lower state bound       */
               double    **PrevStates,     /* (I) Prev state levels       */
               int         t,              /* (I) Current time point      */
               FIX3_TREE_DATA   *tree_data)     /* (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *VolBondL[MAXNBSTATES+1];
    double  *IndexCplxL;
    double  *IndexFundL;
    double  *ZeroToPmtCplxL;
    double  *ZeroToPmtFundL;

    /* state-variable variables */

    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs    */
    int     s;                        /* State variable index            */
    int     sj;                       /* Interp state variable index     */
    int     q;                        /* Index for quad, linear intrp    */
    double  *State  = NULL;           /* Current levels of state var     */
    double  *PState = NULL;           /* Previous levels of state var    */
    double  deltaS;                   /* Increament between cons states  */
    double  InterpStateLevel;         /* state level for interp          */

    /* payoff variables */

    double  Payoff[MAXNBSTATES];      /* Intermediate values of volbond  */
    double  tmpPayoff;
    int     IsLastReset;              /* true=this is last reset date    */
    double  CouponRate, CouponAmt;

    /* Tree variables */

    int     Top1, Bottom1;            /* Tree limits (1st dim)  */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    IsLastReset = (*PrevStates == NULL);

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Fix3_VolSwap_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* FUNDING LEG */

    /* Add funding payoff if it's a reset */
    if (ResetFlagFund)
    {
        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                VolBondL[s] = VolBond[s] + offset;
            }
            IndexFundL     = IndexFund  + offset;
            ZeroToPmtFundL = ZeroToPmtFund + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                CouponRate = COLLAR(IndexFundL[i]+StepUpFund,  \
                                    CapFund,                   \
                                    FloorFund); 

                CouponAmt = OutsFund * ZeroToPmtFundL[i] * 
                            ACC_FN(CouponRate,
                                   DcfFund,
                                   TRUE);  /* simple */

                for (s = 0; s < NbStates ; s++)
                {
                    VolBondL[s][i] -= CouponAmt;
                }

            } /* for i */
        } 

        /************************  2 FACTOR    ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    VolBondL[s] = VolBond[s] + offset;
                }
                IndexFundL     = IndexFund  + offset;
                ZeroToPmtFundL = ZeroToPmtFund + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    CouponRate = COLLAR(IndexFundL[j]+StepUpFund,  \
                                        CapFund,                   \
                                        FloorFund); 
                
                    CouponAmt = OutsFund * ZeroToPmtFundL[j] * 
                                ACC_FN(CouponRate,
                                       DcfFund,
                                       TRUE); /* simple */

                    for (s = 0; s < NbStates ; s++)  
                    {
                        VolBondL[s][j] -= CouponAmt;
                    }
                                    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTOR    ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        VolBondL[s] = VolBond[s] + offset;
                    }

                    IndexFundL     = IndexFund  + offset;
                    ZeroToPmtFundL = ZeroToPmtFund + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        CouponRate = COLLAR(IndexFundL[k]+StepUpFund,  \
                                            CapFund,                   \
                                            FloorFund); 
                
                        CouponAmt = OutsFund * ZeroToPmtFundL[k] * 
                                    ACC_FN(CouponRate,
                                           DcfFund,
                                           TRUE); /* simple */

                        for (s=0; s<NbStates ; s++)  
                        {
                            VolBondL[s][k] -= CouponAmt;
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlagFund */

    /* COMPLEX LEG */

    /* Add complex payoff if it's a reset */
    if (ResetFlagCplx)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /*   and the cap/floor strikes              */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("Fix3_VolSwap_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Fix3_VolSwap_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */

        if (!IsLastReset)
        {
            PState = *PrevStates;

            if (!IS_EQUAL(PState[0], PState[NbStates-1]))
            {
                for (s = 0; s < NbStates-2; s++)
                {
                    D[s][0] = 1. / ((PState[s]-PState[s+1])
                                   *(PState[s]-PState[s+2]));
                    D[s][1] = 1. / ((PState[s+1]-PState[s])
                                   *(PState[s+1]-PState[s+2]));
                    D[s][2] = 1. / ((PState[s+2]-PState[s])
                                   *(PState[s+2]-PState[s+1]));
                }
            }
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                VolBondL[s] = VolBond[s] + offset;
            }
            IndexCplxL     = IndexCplx  + offset;
            ZeroToPmtCplxL = ZeroToPmtCplx + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {
                for (s = 0; s < NbStates; s++)
                {
                    /* ---------------------------- */
                    /*   Add the new coupon amt     */
                    /* ---------------------------- */

                    if (IsCplxPmtOff)
                    {
                        Payoff[s] = 0.0;
                    }
                    else
                    {
                        CouponRate = ABS(IndexCplxL[i] - State[s]) * 
                                     LeverageCplx + StepUpCplx;

                        CouponRate = COLLAR(CouponRate, CapCplx, FloorCplx);

                        Payoff[s] = OutsCplx * ZeroToPmtCplxL[i] * 
                                    ACC_FN(CouponRate,DcfCplx,TRUE);
                    }

                    /* ------------------------------------- */
                    /*   Add the interpolated VolBond price  */
                    /*   to the payoff if appropriate        */
                    /* ------------------------------------- */

                    if (!IsLastReset)
                    {
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] += VolBondL[0][i];
                        }
                        else
                        {
                            InterpStateLevel = IndexCplxL[i];

                            sj = (int) floor((InterpStateLevel - PState[0]) /
                                             (PState[1] - PState[0]));

                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] += VolBondL[NbStates-1][i];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] += VolBondL[0][i];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                    
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         VolBondL[q][i],VolBondL[q+1][i],
                                         VolBondL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] += tmpPayoff;
                            } /* if sj */
                        }
                    } 
                    else
                    {
                        Payoff[s] += VolBondL[0][i];
                    } /* if !IsLastReset */
                } /* for state idx */

                /* replace volbond prices with new values for each state */
                for (s=0; s<NbStates ; s++)
                {
                    VolBondL[s][i] = Payoff[s];
                }

            } /* for i */

        } /* if NbFactor == 1 */

        /************************  2 FACTORS   ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    VolBondL[s] = VolBond[s] + offset;
                }
                IndexCplxL     = IndexCplx  + offset;
                ZeroToPmtCplxL = ZeroToPmtCplx + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ---------------------------- */
                        /*   Add the new coupon amt     */
                        /* ---------------------------- */

                        if (IsCplxPmtOff)
                        {
                            Payoff[s] = 0.0;
                        }
                        else
                        {
                            CouponRate = ABS(IndexCplxL[j] - State[s]) * 
                                         LeverageCplx + StepUpCplx;

                            CouponRate = COLLAR(CouponRate, 
                                                CapCplx, 
                                                FloorCplx);

                            Payoff[s] = OutsCplx * ZeroToPmtCplxL[j] * 
                                        ACC_FN(CouponRate,DcfCplx,TRUE);
                        }

                        /* ------------------------------------- */
                        /*   Add the interpolated VolBond price   */
                        /*   to the payoff if appropriate        */
                        /* ------------------------------------- */
    
                        if (!IsLastReset)
                        {
                            if (IS_EQUAL(PState[0], PState[1])) 
                            {
                                Payoff[s] += VolBondL[0][j];
                            }
                            else
                            {
                                InterpStateLevel = IndexCplxL[j];

                                sj = (int)floor((InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
    
                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] += VolBondL[NbStates-1][j];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += VolBondL[0][j];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q], PState[q+1], 
                                             PState[q+2],
                                             VolBondL[q][j],VolBondL[q+1][j],
                                             VolBondL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] += tmpPayoff;
                                   
                                } /* if sj */
                            }
                        } 
                        else
                        {
                            Payoff[s] += VolBondL[0][j];
                        } 
                    } /* for state idx */

                    /* replace volbond prices with new values for each state */
                    for (s=0; s<NbStates ; s++)  
                    {
                        VolBondL[s][j] = Payoff[s];
                    }

                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTORS   ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        VolBondL[s] = VolBond[s] + offset;
                    }

                    IndexCplxL  = IndexCplx  + offset;
                    ZeroToPmtCplxL = ZeroToPmtCplx + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ---------------------------- */
                            /*   Add the new coupon amt     */
                            /* ---------------------------- */

                            if (IsCplxPmtOff)
                            {
                                Payoff[s] = 0.0;
                            }
                            else
                            {
                                CouponRate = ABS(IndexCplxL[k] - State[s]) * 
                                             LeverageCplx + StepUpCplx;

                                CouponRate = COLLAR(CouponRate, 
                                                    CapCplx, 
                                                    FloorCplx);

                                Payoff[s] = OutsCplx * ZeroToPmtCplxL[k] * 
                                            ACC_FN(CouponRate,DcfCplx,TRUE);
                            }

                            /* ------------------------------------- */
                            /*   Add the interpolated VolBond price  */
                            /*   to the payoff if appropriate        */
                            /* ------------------------------------- */
        
                            if (!IsLastReset)
                            {
                                if (IS_EQUAL(PState[0], PState[1])) 
                                {
                                    Payoff[s] += VolBondL[0][k];
                                }
                                else
                                {
                                    InterpStateLevel = IndexCplxL[k];

                                    sj = (int) floor(
                                                (InterpStateLevel - PState[0])/
                                                (PState[1] - PState[0]));
        
                                    if (sj >= (NbStates - 1)) /* right */
                                    {
                                        Payoff[s] += VolBondL[NbStates-1][k];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] += VolBondL[0][k];
                                    }
                                    else /* in between */
                                    {
                                        q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(
                                              PState[q], PState[q+1], 
                                              PState[q+2],
                                              VolBondL[q][k],VolBondL[q+1][k],
                                              VolBondL[q+2][k],
                                              D[q][0], D[q][1], D[q][2],
                                              InterpStateLevel,
                                              &(tmpPayoff));
                                        Payoff[s] += tmpPayoff;
                                       
                                    } /* if sj */
                                }
                            } 
                            else
                            {
                                Payoff[s] += VolBondL[0][k];
                            }
                        } /* for state idx */

                        /* replace volbond prices for each state */
                        for (s=0; s<NbStates ; s++)  
                        {
                            VolBondL[s][k] = Payoff[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

        /* -------------------------------------- */
        /*   replace PrevStates with new states   */
        /* -------------------------------------- */

        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;

    }  /* if ResetFlagCplx */

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }

    return (status);

}  /* Fix3_VolSwap_t */
