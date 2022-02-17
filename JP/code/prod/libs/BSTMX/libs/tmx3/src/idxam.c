/**************************************************************************/
/*      Payoff function for the index amortizing swap                     */
/**************************************************************************/
/*      idxam.c                                                           */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"




/*****  IdxamSwap_t  *****************************************************/
/*
*   Calculates the sticky/adjustable swap price
*   no dev's, just add sticky and floating coupons
*/
int  IdxamSwap_t
             (double    **Idxam,          /* (O) Prices for all states   */
              long        AmortFlag,      /* (I) Princ pmt at timept     */
              long        ResetFlag,      /* (I) Flt reset at timept     */
              long        FixCpnFlag,     /* (I) Fix pmt at timept       */
              double     *AmortIndex,     /* (I) Amort index             */
              double     *FloatCpn,       /* (I) Float cpn amt           */
              double      FixCpn,         /* (I) Fix cpn amt             */
              double      Notional,       /* (I) Deal notional           */
              double      AmortKnown,     /* (I) Known amort rate        */
              double      AmortBase,      /* (I) Base rate               */
              double      AmortFactor,    /* (I) Amort factor            */
              double      CleanUpLevel,   /* (I) Clean up for outs below */
              double      *ATRate,        /* (I) Stoch amort rate values */
              double      *ATSprd,        /* (I) Stoch amort rate spreads*/
              int         ATDim,          /* (I) Stoch amort rate dim    */
              char        PmtType,        /* (I) Pmt type                */
              char        OoR,            /* (I) Outs or remainging amort*/
              char        ArrearsReset,   /* (I) Idxam/adjustable        */
              int         NbStates,       /* (I) Nb of states            */
              double      MaxState,       /* (I) Upper state bound       */
              double      MinState,       /* (I) Lower state bound       */
              double    **PrevStates,     /* (I) Prev state levels       */
              int         t,              /* (I)   Current time point    */
              int         T,              /* (I)   Last time point       */
              int         DCurve,         /* (I)   Discount curve        */
              DEV_DATA    *dev_data,      /* (I)   Dev data structure    */
              TREE_DATA   *tree_data)     /* (I)   Tree data structure   */
{
    /* Local slice pointers */

    double  *IdxamL[MAXNBSTATES+1];
    double  *AmortIndexL;
    double  *FloatCpnL;

    /* state-variable variables */

    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs    */
    int     s;                        /* State variable index            */
    int     sj;                       /* Interp state variable index     */
    int     q;                        /* Index for quad, linear intrp    */
    double  *State  = NULL;           /* Current levels of state var     */
    double  *PState = NULL;           /* Previous levels of state var    */
    double  deltaS;                   /* Increament between cons states  */
    double  IntpState;         /* state level for interp          */

    /* payoff variables */

    double  Payoff[MAXNBSTATES];      /* Intermediate values of index am */
    double  tmpPayoff;
    int     IsLastReset;              /* true=this is last reset date    */
    double  SignNotional;
    double  AmortStoch;    
    double  NextOuts;
    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    if ((PmtType == 'X') || (PmtType =='P'))
    {
        SignNotional = +Notional;
    }
    else
    if (PmtType == 'L')
    {
        SignNotional = -Notional;
    }
    else
    {
        SignNotional = 0.;
    }


    IsLastReset = (*PrevStates == NULL);

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("IdxamSwap_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* Discount all the state variables */
    
    if (!IsLastReset)
    {
        for (s=0; s<NbStates; s++)
        {                                  
            if (Dev(Idxam[s],
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
    
    /* FLOATING LEG IN ADVANCE */

    /* Add floating payoff if it's a reset-in-advance case */
    if (ResetFlag && (ArrearsReset == 'N'))
    {
        PState = *PrevStates;
        if (PState == NULL)
        {
            DR_Error("IdxamSwap_t: No outstanding available !");
            goto RETURN;
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                IdxamL[s] = Idxam[s] + offset;
            }
            FloatCpnL = FloatCpn + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {   
                for (s = 0; s < NbStates ; s++)  
                {
                    IdxamL[s][i] -= FloatCpnL[i] * PState[s];
                }

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
                    IdxamL[s] = Idxam[s] + offset;
                }
                FloatCpnL = FloatCpn + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        IdxamL[s][j] -= FloatCpnL[j] * PState[s];
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
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        IdxamL[s] = Idxam[s] + offset;
                    }
                    FloatCpnL = FloatCpn + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            IdxamL[s][k] -= FloatCpnL[k] * PState[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlag */


    /* AMORTIZATION */

    if (AmortFlag)
    {
        /*   Prepare state levels for this time pt  */
        if (MaxState < MinState)
        {
            DR_Error("Idxam_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Idxam_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /*   Prepare prev states for interp if appropriate  */
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
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                IdxamL[s] = Idxam[s] + offset;
            }
            AmortIndexL  = AmortIndex  + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                /* find stochastic amortization rate */
                tableinterp(AmortIndexL[i]-AmortBase,
                            &AmortStoch,
                            ATSprd,
                            ATRate,
                            ATDim);
                AmortStoch *= AmortFactor;

                for (s = 0; s < NbStates; s++)
                {
                    if (!IsLastReset)
                    {
                        /* find new outstanding  */
                        if (OoR == 'O')
                        {
                            NextOuts = State[s]
                                     - AmortKnown 
                                     - AmortStoch;
                        }
                        else
                        {
                            NextOuts = State[s]
                                     * (1.-AmortKnown) 
                                     * (1.-AmortStoch);
                        }
                    
                        /* make pmts according to next outstanding */
                        if (NextOuts < CleanUpLevel)
                        {
                            Payoff[s] = State[s] * SignNotional;
                        }
                        else
                        {
                            /* possible princ pmt */ 
                            Payoff[s] = (State[s] - NextOuts) * SignNotional;

                            /* add the interpolated Idxam price   */

                            if (IS_EQUAL(PState[0], PState[1]))
                            {
                                Payoff[s] += IdxamL[0][i];
                            }
                            else
                            {
                                IntpState = NextOuts;
                                sj = (int) floor((IntpState - PState[0]) /
                                                 (PState[1] - PState[0]));

                                if (sj >= (NbStates - 1)) /* to the right */
                                {
                                    Payoff[s] += IdxamL[NbStates-1][i];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += IdxamL[0][i];
                                }
                                else /* in between */
                                {
                                   q = MIN(sj, NbStates - 3);
                        
                                    sqinterp(PState[q],PState[q+1],PState[q+2],
                                             IdxamL[q][i],IdxamL[q+1][i],
                                             IdxamL[q+2][i],
                                             D[q][0], D[q][1], D[q][2],
                                             IntpState,
                                             &(tmpPayoff));
                                    Payoff[s] += tmpPayoff;
                                   
                                } /* if sj */
                            } /* if possible to interp */
                        } /* if next oust < clean up level */
                    } 
                    else
                    {
                        Payoff[s] = State[s] * SignNotional;
                    } /* if !IsLastReset */

                } /* for state idx */

                /* replace sticky prices with new values for each state */

                for (s=0; s<NbStates ; s++)  
                {
                    IdxamL[s][i] = Payoff[s];
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
                    IdxamL[s] = Idxam[s] + offset;
                }
                AmortIndexL  = AmortIndex  + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    /* find stochastic amortization rate */
                    tableinterp(AmortIndexL[j]-AmortBase,
                                &AmortStoch,
                                ATSprd,
                                ATRate,
                                ATDim);
                    AmortStoch *= AmortFactor;

                    for (s = 0; s < NbStates; s++)
                    {
                        if (!IsLastReset)
                        {
                            /* find new outstanding  */
                            if (OoR == 'O')
                            {
                                NextOuts = State[s]
                                         - AmortKnown 
                                         - AmortStoch;
                            }
                            else
                            {
                                NextOuts = State[s]
                                         * (1.-AmortKnown) 
                                         * (1.-AmortStoch);
                            }
                    
                            /* make pmts according to next outstanding */
                            if (NextOuts < CleanUpLevel)
                            {
                                Payoff[s] = State[s] * SignNotional;
                            }
                            else
                            {
                                /* possible princ pmt */ 
                                Payoff[s] = (State[s] - NextOuts) * SignNotional;
            
                                /* add the interpolated Idxam price   */
                    
                                if (IS_EQUAL(PState[0], PState[1]))
                                {
                                    Payoff[s] += IdxamL[0][j];
                                }
                                else
                                {
                                    IntpState = NextOuts;
                                    sj = (int) floor((IntpState - PState[0]) /
                                                     (PState[1] - PState[0]));

                                    if (sj >= (NbStates - 1)) /*to the right*/
                                    {
                                        Payoff[s] += IdxamL[NbStates-1][j];
                                    }
                                    else if (sj < 0) /* to the left */
                                    {
                                        Payoff[s] += IdxamL[0][j];
                                    }
                                    else /* in between */
                                    {
                                       q = MIN(sj, NbStates - 3);
                            
                                        sqinterp(PState[q],PState[q+1],PState[q+2],
                                                 IdxamL[q][j],IdxamL[q+1][j],
                                                 IdxamL[q+2][j],
                                                 D[q][0], D[q][1], D[q][2],
                                                 IntpState,
                                                 &(tmpPayoff));
                                        Payoff[s] += tmpPayoff;
                                       
                                    } /* if sj */
                                } /* if possible to interp */
                            } /* if next oust < clean up level */
                        } 
                        else
                        {
                            Payoff[s] = State[s] * SignNotional;
                        } /* if !IsLastReset */

                    } /* for state idx */

                    /* replace sticky prices with new values for each state */
        
                    for (s=0; s<NbStates ; s++)  
                    {
                        IdxamL[s][j] = Payoff[s];
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
                        IdxamL[s] = Idxam[s] + offset;
                    }
                    AmortIndexL  = AmortIndex  + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        /* find stochastic amortization rate */
                        tableinterp(AmortIndexL[k]-AmortBase,
                                    &AmortStoch,
                                    ATSprd,
                                    ATRate,
                                    ATDim);
                        AmortStoch *= AmortFactor;
        
                        for (s = 0; s < NbStates; s++)
                        {
                            if (!IsLastReset)
                            {
                                /* find new outstanding  */
                                if (OoR == 'O')
                                {
                                    NextOuts = State[s]
                                             - AmortKnown 
                                             - AmortStoch;
                                }
                                else
                                {
                                    NextOuts = State[s]
                                             * (1.-AmortKnown) 
                                             * (1.-AmortStoch);
                                }
                        
                                /* make pmts according to next outstanding */
                                if (NextOuts < CleanUpLevel)
                                {
                                    Payoff[s] = State[s] * SignNotional;
                                }
                                else
                                {
                                    /* possible princ pmt */ 
                                    Payoff[s] = (State[s] - NextOuts) 
                                              * SignNotional;
                    
                                    /* add the interpolated Idxam price   */
                            
                                    if (IS_EQUAL(PState[0], PState[1]))
                                    {
                                        Payoff[s] += IdxamL[0][k];
                                    }
                                    else
                                    {
                                        IntpState = NextOuts;
                                        sj = (int) floor((IntpState - PState[0]) /
                                                         (PState[1] - PState[0]));
                                
                                        if (sj >= (NbStates - 1))
                                        {
                                            Payoff[s] += IdxamL[NbStates-1][k];
                                        }
                                        else if (sj < 0) /* to the left */
                                        {
                                            Payoff[s] += IdxamL[0][k];
                                        }
                                        else /* in between */
                                        {
                                           q = MIN(sj, NbStates - 3);
                                    
                                           sqinterp(PState[q],PState[q+1],PState[q+2],
                                                    IdxamL[q][k],IdxamL[q+1][k],
                                                    IdxamL[q+2][k],
                                                    D[q][0], D[q][1], D[q][2],
                                                    IntpState,
                                                    &(tmpPayoff));
                                            Payoff[s] += tmpPayoff;
                                       
                                        } /* if sj */
                                   } /* if possible to interp */
                                } /* if next oust < clean up level */
                            } 
                            else
                            {
                                Payoff[s] = State[s] * SignNotional;
                            } /* if !IsLastReset */
                        
                        } /* for state idx */

                        /* replace sticky prices with new values */
            
                        for (s=0; s<NbStates ; s++)  
                        {
                            IdxamL[s][k] = Payoff[s];
                        }
                    }  /* for k */
                }  /* for j */
        }  /* if then else */

        /*   replace PrevStates with new states   */

        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevStates = State;
        State = NULL;

    }  /* if AmortFlag */


    /* FLOATING LEG IN ARREARS */

    /* Add floating payoff if it's a reset-in-advance case */
    if (ResetFlag && (ArrearsReset == 'Y'))
    {
        /* PState is new outstanding at this moment */
        PState = *PrevStates;
        if (PState == NULL)
        {
            DR_Error("IdxamSwap_t: No outstanding available !");
            goto RETURN;
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                IdxamL[s] = Idxam[s] + offset;
            }
            FloatCpnL = FloatCpn + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {   
                for (s = 0; s < NbStates ; s++)  
                {
                    IdxamL[s][i] -= FloatCpnL[i] * PState[s];
                }

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
                    IdxamL[s] = Idxam[s] + offset;
                }
                FloatCpnL = FloatCpn + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        IdxamL[s][j] -= FloatCpnL[j] * PState[s];
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
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        IdxamL[s] = Idxam[s] + offset;
                    }
                    FloatCpnL = FloatCpn + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            IdxamL[s][k] -= FloatCpnL[k] * PState[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlag */


    /* FIX LEG */

    /* Add floating payoff if it's a reset-in-advance case */
    if (FixCpnFlag)
    {
        PState = *PrevStates;
        if (PState == NULL)
        {
            DR_Error("IdxamSwap_t: No outstanding available !");
            goto RETURN;
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                IdxamL[s] = Idxam[s] + offset;
            }

            for (i = Bottom1; i <= Top1; i ++)
            {   
                for (s = 0; s < NbStates ; s++)  
                {
                    IdxamL[s][i] += FixCpn * PState[s];
                }

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
                    IdxamL[s] = Idxam[s] + offset;
                }
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        IdxamL[s][j] += FixCpn * PState[s];
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
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        IdxamL[s] = Idxam[s] + offset;
                    }
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            IdxamL[s][k] += FixCpn * PState[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if FixCpnFlag */

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }
    return (status);

}  /* IdxamSwap_t */
