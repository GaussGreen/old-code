/**************************************************************************/
/*      Payoff function for the cub loan                                  */
/**************************************************************************/
/*      loan.c                                                          */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"
#define  MRATIO 0.000001




/*****  Fix3_CubLoan_t  *********************************************************/
/*
*   Calculates the cub loan price
*   dev's the loan to this time point (t)
*   if FundResetFlag then multiply by princ ratio and revise state levels
*   if AnnResetFlag then cap annuity and revise state levels
*   if *PrevStates is NULL, then assume Loan is uninitialised    
*/
int  Fix3_CubLoan_t (double    **Loan,         /* (I/O) Prices for all states   */
                long      FundResetFlag,  /* (I)   Funding reset at timept */
                long      AnnResetFlag,   /* (I)   Annuity reset at timept */
                double    *FltIndex,      /* (I)   Index level             */
                double    CapRate,        /* (I)                           */
                double    FloorRate,      /* (I)                           */
                double    MaxStateFund,   /* (I)   Upper state bound       */
                double    MinStateFund,   /* (I)   Lower state bound       */
                double    MaxStateAnn,    /* (I)   Upper state bound       */
                double    MinStateAnn,    /* (I)   Lower state bound       */
                double    RemTime,        /* (I)   Time to ann maturity    */
                int       PmtFreq,        /* (I)   Payment frequncy        */
                int       PPF,            /* (I)   Compounding frequency   */
                int       NbStates,       /* (I)   Nb of states            */
                double    **PrevStates,   /* (I)   Prev state levels       */
                int       t,              /* (I)   Current time point      */
                int       T,              /* (I)   Last time point         */
                int       DCurve,         /* (I)   Discount curve          */
                FIX3_DEV_DATA  *dev_data,      /* (I)   Fix3_Dev data structure      */
                FIX3_TREE_DATA *tree_data)     /* (I)   Tree data structure     */
{
    /* Local slice pointers */

    double  *LoanL[MAXNBSTATES+1];
    double  *FltIndexL;

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

    double  Payoff[MAXNBSTATES];      /* Intermediate values of loan   */
    double  tmpPayoff;
    double  CapState[MAXNBSTATES];
    double  FloorState[MAXNBSTATES];
    double  MinState, MaxState;
    double  Ratio, FullRatio;
    double  M;
    double  Time;
    int     IsLastReset;              /* true=this is last reset date    */
    int     c;
    
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
        DR_Error("Loan_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* Discount all the state variables */
    
    for (s=0; s<NbStates; s++)
    {                                  
        if (Fix3_Dev(Loan[s],
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }

    /*  Account for prinicipal change if funding reset */

    if (FundResetFlag)
    {
        /* Prepare state levels for this time pt  */
        /* and the cap/floor ratios               */
        if (AnnResetFlag)
        {
            MinState = MinStateAnn;  MaxState = MaxStateAnn;
        }
        else
        {
            MinState = MinStateFund; MaxState = MaxStateFund;
        }

        if (MaxState < MinState)
        {
            DR_Error("Loan_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Loan_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /* Prepare prev states for interp if appropriate  */
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
                LoanL[s] = Loan[s] + offset;
            }
            FltIndexL  = FltIndex  + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates; s++)
                {
                    /* Calculate next state due to amortization */
                    FullRatio = 1.;
                    M = State[s];
                    for (c=0; c<PPF; c++)
                    {
                        Ratio = 1. + (FltIndexL[i] - M) / PmtFreq;
                        if (Ratio < MRATIO) break;
                        M /= Ratio;
                        FullRatio *= Ratio;
                    }
                
                    /* Check if not fully amortized and interpolate */
                    if (c<PPF)
                    {
                        Payoff[s] = 0.;
                    }
                    else
                    if (IsLastReset)  
                    {
                        /* values for all M are the same */
                        Payoff[s] = LoanL[s][i];
                        
                        /* account for principal ratio */
                        Payoff[s] *= FullRatio;
                    }
                    else
                    {
                        /* pick up next value */
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] = LoanL[0][i];
                        }
                        else
                        {
                            InterpStateLevel = M;
                            sj = (int) floor((InterpStateLevel - PState[0]) /
                                             (PState[1] - PState[0]));

                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] = LoanL[NbStates-1][i];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] = LoanL[0][i];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                    
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         LoanL[q][i],LoanL[q+1][i],
                                         LoanL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] = tmpPayoff;
                               
                            } /* if sj */
                        }
                        /* account for principal ratio */
                        Payoff[s] *= FullRatio;

                    }/* if IsLastReset */
                    
                } /* for state idx */

                /* replace loan prices with new values for each state */

                for (s=0; s<NbStates ; s++)  
                {
                    LoanL[s][i] = Payoff[s];
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
                    LoanL[s] = Loan[s] + offset;
                }
                FltIndexL  = FltIndex  + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* Calculate next state due to amortization */
                        FullRatio = 1.;
                        M = State[s];
                        for (c=0; c<PPF; c++)
                        {
                            Ratio = 1. + (FltIndexL[j] - M) / PmtFreq;
                            if (Ratio < MRATIO) break;
                            M /= Ratio;
                            FullRatio *= Ratio;
                        }
                
                        /* Check if not fully amortized and interpolate */
                        if (c<PPF)
                        {
                            Payoff[s] = 0.;
                        }
                        else
                        if (IsLastReset)  
                        {
                            /* values for all M are the same */
                            Payoff[s] = LoanL[s][j];
                            
                            /* account for principal ratio */
                            Payoff[s] *= FullRatio;
                        }
                        else
                        {
                            /* pick up next value */
                            if (IS_EQUAL(PState[0], PState[1]))
                            {
                                Payoff[s] = LoanL[0][j];
                            }
                            else
                            {
                                InterpStateLevel = M;
                                sj = (int) floor((InterpStateLevel - PState[0]) /
                                                 (PState[1] - PState[0]));
                        
                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] = LoanL[NbStates-1][j];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] = LoanL[0][j];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q], PState[q+1], PState[q+2],
                                             LoanL[q][j],LoanL[q+1][j],
                                             LoanL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] = tmpPayoff;
                                   
                                } /* if sj */
                            }
                            /* account for principal ratio */
                            Payoff[s] *= FullRatio;
                            
                        }/* if IsLastReset */
                        
                    } /* for state idx */

                    /* replace loan prices with new values for each state */

                    for (s=0; s<NbStates ; s++)  
                    {
                        LoanL[s][j] = Payoff[s];
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
                        LoanL[s] = Loan[s] + offset;
                    }
                    FltIndexL  = FltIndex  + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* Calculate next state due to amortization */
                            FullRatio = 1.;
                            M = State[s];
                            for (c=0; c<PPF; c++)
                            {
                                Ratio = 1. + (FltIndexL[k] - M) / PmtFreq;
                                if (Ratio < MRATIO) break;
                                M /= Ratio;
                                FullRatio *= Ratio;
                            }
                    
                            /* Check if not fully amortized and interpolate */
                            if (c<PPF)
                            {
                                Payoff[s] = 0.;
                            }
                            else
                            if (IsLastReset)  
                            {
                                /* values for all M are the same */
                                Payoff[s] = LoanL[s][k];
                                
                                /* account for principal ratio */
                                Payoff[s] *= FullRatio;
                            }
                            else
                            {
                                /* pick up next value */
                                if (IS_EQUAL(PState[0], PState[1]))
                                {
                                    Payoff[s] = LoanL[0][k];
                                }
                                else
                                {
                                    InterpStateLevel = M;
                                    sj = (int) floor((InterpStateLevel - PState[0]) /
                                                     (PState[1] - PState[0]));
                                    
                                    if (sj >= (NbStates - 1)) /* to the right */
                                    {
                                        Payoff[s] = LoanL[NbStates-1][k];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] = LoanL[0][k];
                                    }
                                    else /* in between */
                                    {
                                        q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(PState[q], PState[q+1], PState[q+2],
                                                 LoanL[q][k],LoanL[q+1][k],
                                                 LoanL[q+2][k],
                                                 D[q][0], D[q][1], D[q][2],
                                                 InterpStateLevel,
                                                 &(tmpPayoff));
                                        Payoff[s] = tmpPayoff;
                                       
                                    } /* if sj */
                                }
                                /* account for principal ratio */
                                Payoff[s] *= FullRatio;
                                
                            }/* if IsLastReset */
                        
                        } /* for state idx */
        
                        /* replace loan prices with new values for each state */

                        for (s=0; s<NbStates ; s++)  
                        {
                            LoanL[s][k] = Payoff[s];
                        }
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        
        /* replace PrevStates with new states   */
        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;
        State = NULL;

    }  /* if FundResetFlag */


    /*  Account for annuity coupon change if annuity reset */
    IsLastReset = (*PrevStates == NULL);

    if (AnnResetFlag)
    {
        /* Prepare state levels for this time pt  */
        /* and the cap/floor ratios               */
        MinState = MinStateFund; MaxState = MaxStateFund;

        if (MaxState < MinState)
        {
            DR_Error("Loan_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Loan_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /* only for annuity reset */
        for (s=0; s<NbStates; s++)
        {
            CapState[s]   = State[s] * (1.+CapRate);
            FloorState[s] = State[s] * (1.+FloorRate);
        }

        /* Prepare prev states for interp if appropriate  */
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
        else
        {
            DR_Error("Loan_t: Prev state does not exist for annuity reset!");
            goto RETURN;
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                LoanL[s] = Loan[s] + offset;
            }
            FltIndexL  = FltIndex  + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates; s++)
                {
                    /* Calculate next state due to annuity reset */
                    if (RemTime > 0)
                    {
                        Time = RemTime * PmtFreq / 12.;
                        M = FltIndexL[i] 
                          / (1. - pow(1. + FltIndexL[i] / PmtFreq, -Time));
                    }
                    else
                    {
                        M = State[s];
                    }
                    M = COLLAR(M, CapState[s], FloorState[s]); 
                    
                    if (IsLastReset)
                    {
                        DR_Error("Loan_t: Prev state does not exist for annuity reset!");
                        goto RETURN;
                    }
                    else
                    {
                        /* pick up next value */
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] = LoanL[0][i];
                        }
                        else
                        {
                            InterpStateLevel = M;
                            sj = (int) floor((InterpStateLevel - PState[0]) /
                                             (PState[1] - PState[0]));

                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] = LoanL[NbStates-1][i];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] = LoanL[0][i];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                    
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         LoanL[q][i],LoanL[q+1][i],
                                         LoanL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] = tmpPayoff;
                               
                            } /* if sj */
                        }
                    }/* if IsLastReset */
                    
                } /* for state idx */

                /* replace loan prices with new values for each state */

                for (s=0; s<NbStates ; s++)  
                {
                    LoanL[s][i] = Payoff[s];
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
                    LoanL[s] = Loan[s] + offset;
                }
                FltIndexL  = FltIndex  + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* Calculate next state due to annuity reset */
                        if (RemTime > 0)
                        {
                            Time = RemTime * PmtFreq / 12.;
                            M = FltIndexL[j] 
                              / (1. - pow(1. + FltIndexL[j] / PmtFreq, -Time));
                        }
                        else
                        {
                            M = State[s];
                        }
                        M = COLLAR(M, CapState[s],FloorState[s]); 
                        
                        if (IsLastReset)
                        {
                            DR_Error("Loan_t: Prev state does not exist for annuity reset!");
                            goto RETURN;
                        }
                        else
                        {
                            /* pick up next value */
                            if (IS_EQUAL(PState[0], PState[1]))
                            {
                                Payoff[s] = LoanL[0][j];
                            }
                            else
                            {
                                InterpStateLevel = M;
                                sj = (int) floor((InterpStateLevel - PState[0]) /
                                                 (PState[1] - PState[0]));
                            
                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] = LoanL[NbStates-1][j];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] = LoanL[0][j];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q], PState[q+1], PState[q+2],
                                             LoanL[q][j],LoanL[q+1][j],
                                             LoanL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] = tmpPayoff;
                                   
                                } /* if sj */
                            }
                        }/* if IsLastReset */
                        
                    } /* for state idx */
    
                    /* replace loan prices with new values for each state */

                    for (s=0; s<NbStates ; s++)  
                    {
                        LoanL[s][j] = Payoff[s];
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
                        LoanL[s] = Loan[s] + offset;
                    }
                    FltIndexL  = FltIndex  + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* Calculate next state due to annuity reset */
                            if (RemTime > 0)
                            {
                                Time = RemTime * PmtFreq / 12.;
                                M = FltIndexL[k] 
                                  / (1. - pow(1. + FltIndexL[k] / PmtFreq, -Time));
                            }
                            else
                            {
                                M = State[s];
                            }
                            M = COLLAR(M, CapState[s], FloorState[s]); 
                            
                            if (IsLastReset)
                            {
                                DR_Error("Loan_t: Prev state does not exist for annuity reset!");
                                goto RETURN;
                            }
                            else
                            {
                                /* pick up next value */
                                if (IS_EQUAL(PState[0], PState[1]))
                                {
                                    Payoff[s] = LoanL[0][k];
                                }
                                else
                                {
                                    InterpStateLevel = M;
                                    sj = (int) floor((InterpStateLevel - PState[0]) /
                                                     (PState[1] - PState[0]));
                                    
                                    if (sj >= (NbStates - 1)) /* to the right */
                                    {
                                        Payoff[s] = LoanL[NbStates-1][k];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] = LoanL[0][k];
                                    }
                                    else /* in between */
                                    {
                                        q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(PState[q], PState[q+1], PState[q+2],
                                                 LoanL[q][k],LoanL[q+1][k],
                                                 LoanL[q+2][k],
                                                 D[q][0], D[q][1], D[q][2],
                                                 InterpStateLevel,
                                                 &(tmpPayoff));
                                        Payoff[s] = tmpPayoff;
                                       
                                    } /* if sj */
                                }
                            }/* if IsLastReset */
                            
                        } /*  for state idx */
        
                        /* replace loan prices with new values for each state */
                
                        for (s=0; s<NbStates ; s++)  
                        {
                            LoanL[s][k] = Payoff[s];
                        }       
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        
        /* replace PrevStates with new states   */
        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;

    }  /* if AnnResetFlag */

    status = SUCCESS;
    
RETURN:

    return (status);

}  /* Loan_t */
