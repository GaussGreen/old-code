/**************************************************************************/
/*      Payoff function for the ladder swap                               */
/**************************************************************************/
/*      ladder.c                                                          */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bmx123head.h"


/*****  LadderSwap_t  *****************************************************/
/*
*   Calculates the sticky/adjustable swap price
*   no dev's, just add sticky and floating coupons
*/
int  LadderSwap_t
             (double    **Sticky,         /* (O) Prices for all states   */
              long        ResetFlagSt,    /* (I) Reset at timept         */
              long        ResetFlagFund,  /* (I) Reset at timept         */
              double     *Step,           /* (I) Step for ladder         */
              double     *Funding,        /* (I) Funding                 */
              double     *ZeroToPmtSt,    /* (I) PmtZero 4 this reset    */
              double      FloorSt,        /* (I) Floor spread            */
              double      CapSt,          /* (I) Cap spread              */
              double      OutsSt,         /* (I) Outstanding             */
              double      DcfSt,          /* (I) Day count fract         */
              char        SoZ,            /* (I) Swap/zero coupon        */
              char        CompSt,         /* (I) Simple/Compound pmt     */
              int         NbStates,       /* (I) Nb of states            */
              double      MaxState,       /* (I) Upper state bound       */
              double      MinState,       /* (I) Lower state bound       */
              double    **PrevStates,     /* (I) Prev state levels       */
              int         t,              /* (I) Current time point      */
              TREE_DATA   *tree_data)     /* (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *StickyL[MAXNBSTATES+1];
    double  *FundingL;
    double  *ZeroToPmtStL;
    double  *StepL;

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
    double  PRate[MAXNBSTATES];
    double  X;
    int     isZeroSt, isSimpleSt;

    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    if (*PrevStates == NULL)
    {
        DR_Error("LadderSwap_t: No previous states !");
        goto RETURN;
    }

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("LadderSwap_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* STICKY LEG */

    isZeroSt   = (SoZ == 'Z');
    isSimpleSt = (CompSt == 'S');
    /* Add sticky payoff if it's a reset */
    if (ResetFlagSt)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("LadderSwap_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("LadderSwap_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */

        PState = *PrevStates;

        if (!IS_EQUAL(PState[0],PState[NbStates-1]))
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

        /* Note: outs is 1 for zero coupon */
        for (s=0; s<NbStates; s++)
        {
            PRate[s] =  ACC_FN(PState[s],DcfSt,isSimpleSt);
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            ZeroToPmtStL = ZeroToPmtSt + offset;
            StepL  = Step  + offset;


            for (i = Bottom1; i <= Top1; i ++)
            {                
                /* ---------------------------- */
                /*   Update current value:      */
                /*   mulitply and add coupon    */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    if (isZeroSt) StickyL[s][i] *= (1.+PRate[s]);
                    StickyL[s][i] += ZeroToPmtStL[i]*PRate[s]*OutsSt;
                }

                /* ---------------------------- */
                /*   Rearrange variables        */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    /* ------------------------------------- */
                    /*   Find next coupon level              */
                    /* ------------------------------------- */
                    
                    X  = State[s];
                    X += StepL[i]; 
                    X  = COLLAR(X,CapSt,FloorSt);

                    /* ------------------------------------- */
                    /*   Payoff = interpolated Sticky price  */
                    /* ------------------------------------- */

                    Payoff[s] = 0.;
                    if (IS_EQUAL(PState[1],PState[0])) /* one state level */
                    {
                        sj = 0;
                        Payoff[s] += StickyL[0][i];
                    }
                    else
                    {
                        InterpStateLevel = X;
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
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                StepL  = Step  + offset;
                ZeroToPmtStL = ZeroToPmtSt + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    /* ---------------------------- */
                    /*   Update current value:      */
                    /*   mulitply and add coupon    */
                    /* ---------------------------- */
                    for (s = 0; s < NbStates; s++)
                    {
                        if (isZeroSt) StickyL[s][j] *= (1.+PRate[s]);
                        StickyL[s][j] += ZeroToPmtStL[j]*PRate[s]*OutsSt;
                    }
                
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ------------------------------------- */
                        /*   Find next coupon level              */
                        /* ------------------------------------- */
                    
                        X  = State[s];
                        X += StepL[j];
                        X  = COLLAR(X,CapSt,FloorSt);

                        /* ------------------------------------- */
                        /*   Payoff = interpolated Sticky price  */
                        /* ------------------------------------- */
        
                        Payoff[s] = 0.;
                        if (IS_EQUAL(PState[1],PState[0])) 
                        {
                            Payoff[s] += StickyL[0][j];
                        }
                        else
                        {
                            InterpStateLevel = X;
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
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }
                    StepL  = Step  + offset;
                    ZeroToPmtStL = ZeroToPmtSt + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        /* ---------------------------- */
                        /*   Update current value:      */
                        /*   mulitply and add coupon    */
                        /* ---------------------------- */
                        for (s = 0; s < NbStates; s++)
                        {
                            if (isZeroSt) StickyL[s][k] *= (1.+PRate[s]);
                            StickyL[s][k] += ZeroToPmtStL[k]*PRate[s]*OutsSt;
                        }
        
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ------------------------------------- */
                            /*   Find next coupon level              */
                            /* ------------------------------------- */
                    
                            X  = State[s];
                            X += StepL[k];
                            X  = COLLAR(X,CapSt,FloorSt);

                            /* ------------------------------------- */
                            /*   Payoff = interpolated Sticky price  */
                            /* ------------------------------------- */
            
                            Payoff[s] = 0.;
                            if (IS_EQUAL(PState[1],PState[0])) 
                            {
                                Payoff[s] += StickyL[0][k];
                            }
                            else
                            {
                                InterpStateLevel = X;
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
        State = NULL;

    }  /* if ResetFlagSticky */

    /* FUNDING LEG */

    /* Add floating payoff if it's a reset BUT ONLY if it's sticky date */
    /* This restriction is achieved by setting flags in calc properly   */
    if (ResetFlagFund)
    {

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            FundingL = Funding + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates ; s++)  
                {
                    StickyL[s][i] -= FundingL[i];
                }
                FundingL[i] = 0.;

            } /* for i */

        } 

        /************************  2 FACTOR    ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                FundingL = Funding + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        StickyL[s][j] -= FundingL[j];
                    }
                    FundingL[j] = 0.;
                                    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTOR    ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }
                    FundingL = Funding + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            StickyL[s][k] -= FundingL[k];
                        }
                        FundingL[k] = 0.;

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlagFloat */

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }

    return (status);

}  /* LadderSwap_t */



/*****  LadderStep_t  ********************************************************/
/*                                                                           
*       Smooth 3-level step function.
*/
int     LadderStep_t (double      *Step,       /* (I/O) Smooth step          */
                      double      *Index,      /* (I) Index prices           */
                      double      LowBarrier,  /* (I) Lower barrier          */
                      double      HighBarrier, /* (I) Higher barrier         */
                      double      DownValue,   /* (I) Down value             */
                      double      MidValue,    /* (I) Mid value              */
                      double      UpValue,     /* (I) Up value               */
                      char        Smoothing,   /* (I) Smoothing ('Y' or 'N') */
                      int         t,           /* (I) Current time point     */
                      TREE_DATA   *tree_data)  /* (I) Structure of tree data */
{

    double  *StepL;             /* Local slice pointers */
    double  *IndexL;

    double  x;                  /* Intermediate smoothed value */
    double  IndexStep;          /* Diff between index values   */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int     i, j, k;            /* Node indices           */
    int     offset;             /* Node offset            */
    int     status = FAILURE;   /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        StepL = Step + offset;
        IndexL = Index + offset;

        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                if (IndexL[i] <= LowBarrier)
                {
                    StepL[i] = DownValue;
                }
                else 
                if (IndexL[i] <= HighBarrier)
                {        
                    StepL[i] = MidValue;
                }
                else
                {
                    StepL[i] = UpValue;
                }
            }  /* for i */  
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                IndexStep = GetIndexStep (Index,
                                          1,
                                          i, 0, 0,
                                          t,
                                          tree_data);
                
                /* 'up' value = mid value    */
                /* 'down' value = down value */
                if (Smooth_Step (   &x,
                                    MidValue, 
                                    DownValue,
                                    IndexL[i],
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                {
                    goto RETURN;
                }
                                
                /* 'up' value = up value         */
                /* 'down' value = smoothed value */
                if (Smooth_Step (   &(StepL[i]),
                                    UpValue,
                                    x,
                                    IndexL[i],
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                {
                    goto RETURN;
                }
            }  /* for i */  
        }  /* if then else */   
    }
    else if (tree_data->NbFactor == 2)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);
  
                StepL = Step + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    if (IndexL[j] <= LowBarrier)
                    {
                        StepL[j] = DownValue;
                    }
                    else 
                    if (IndexL[j] <= HighBarrier)
                    {        
                        StepL[j] = MidValue;
                    }
                    else
                    {
                        StepL[j] = UpValue;
                    }
                }  /* for j */  
            }  /* for i */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                StepL = Step + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    IndexStep = GetIndexStep (Index,
                                              2,
                                              i, j, 0,
                                              t,
                                              tree_data);
                                                  
                    /* 'up' value = mid value    */
                    /* 'down' value = down value */
                    if (Smooth_Step (   &x,
                                        MidValue,
                                        DownValue,
                                        IndexL[j],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                                
                    /* 'up' value = up value         */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step (   &(StepL[j]),
                                        UpValue,
                                        x,
                                        IndexL[j],
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */   
    }
    else if (tree_data->NbFactor == 3)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)          
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    StepL = Step + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        if (IndexL[k] <= LowBarrier)
                        {
                            StepL[k] = DownValue;
                        }
                        else 
                        if (IndexL[k] <= HighBarrier)
                        {        
                            StepL[k] = MidValue;
                        }
                        else
                        {
                            StepL[k] = UpValue;
                        }
                    }  /* for k */  
                }  /* for j */  
        }
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    StepL = Step + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        IndexStep = GetIndexStep (Index,
                                                  3,
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                                                  
                        /* 'up' value = mid value    */
                        /* 'down' value = down value */
                        if (Smooth_Step (   &x,
                                            MidValue,
                                            DownValue,
                                            IndexL[k],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                    
                        /* 'up' value = up value         */
                        /* 'down' value = smoothed value */
                        if (Smooth_Step (   &(StepL[k]),
                                            UpValue,
                                            x,
                                            IndexL[k],
                                            HighBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }  /* for k */  
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* LadderStep_t */

