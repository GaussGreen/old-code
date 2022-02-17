/**************************************************************************/
/*      Amortloan                                                         */
/**************************************************************************/
/*      AMORTLOAN.c                                                       */
/**************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"





/*****  Fix3_AnnAmortSwap_t  ****************************************************/
/*
*       Amort float price as a function of state variable (outstanding).
*/
int     Fix3_AnnAmortSwap_t 
           (double    **Amortswap,    /* (I/O) Array of amort loans        */
            double     *Collaret,     /* (I) Current collaret              */
            double     *Zero,         /* (I) Zero maturing at next loanpmt */
            double     *Floater,      /* (I) Current funding coupon        */
            char        Leg,          /* (I) Loan, fund or swap leg flag   */
            char        Arrears,      /* (I) 'Y' if reset in arrears       */
            double      FixCoupon,    /* (I) Fixed coupon amount           */
            double      Fee,          /* (I) Fee amount                    */
            double      FinalWeight,  /* (I) Final principal weight        */
            long        CollaretFlag, /* (I) Loan payment day flag         */
            long        FloaterFlag,  /* (I) Funding payment day flat      */
            long        FinalMatFlag, /* (I) Final maturity day flag       */
            int         NbStates,     /* (I) Maximum number of amort loans */
            double      OutsMin,      /* (I) Minimum outstanding           */
            double      OutsMax,      /* (I) Maximum outstanding           */
            double    **PrevState,    /* (I) Ptr to array of prev states   */
            int         t,            /* (I) Current time point            */
            int         T,            /* (I) Last time point               */
            int         DCurve,       /* (I) Discount curve                */
            FIX3_DEV_DATA    *dev_data,    /* (I) Fix3_Dev data structure            */
            FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure           */
{

    double  *AmortswapL[MAXNBSTATES];  /* Local slices pointers             */
    double  *CollaretL;
    double  *ZeroL;
    double  *FloaterL;

    double  D[MAXNBSTATES][3];         /* Precomputed quadratic coeffs      */
    double  *State = NULL;             /* Current levels of state var       */
    double  *PState = NULL;            /* Previous levels of state var      */
    double  IC[MAXNBSTATES];           /* Intermediate values of amortloan  */
    double  deltaS;                    /* Increament between cons states    */
    double  NState;                    /* Next state                        */
    double  IntPmt;                    /* Loan interest payment             */
    double  FltRate;                   /* Capped and floored floating rate  */
    double  PrincPmt;                  /* Principal pmt as part of loan pmt */
               
    int     Top1, Bottom1;             /* Tree limits (1rst dim)            */
    int     *Top2, *Bottom2;           /* Tree limits (2nd dim)             */
    int     **Top3, **Bottom3;         /* Tree limits (3rd dim)             */

    int     i, j, k;                   /* Node indices                      */
    int     offset;                    /* Node offset                       */
    int     s;                         /* State variable index              */
    int     sj;                        /* Interp state variable index       */
    int     q;                         /* Index for quadratic interp        */
    int     status = FAILURE;          /* Error status                      */
            

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Fix3_AnnAmortSwap_t: nb of state vars exceeds limit (200)!");
        goto RETURN;
    }

    if (OutsMax < TINY)
    {
        if (*PrevState != NULL)
        {
            DR_Error("Fix3_AnnAmortSwap_t: outstanding is zero!");
            status = FAILURE;
        }
        else
            status = SUCCESS;

        goto RETURN;
    }

    /* Discount all the state variables */
    
    for (s = 0; s < NbStates; s++)
    {                                  
        if (Fix3_Dev (Amortswap[s],
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
            
        }  /* if */
    }  /* for l */

    /* Amortloan state variable propagation */
    
    if (CollaretFlag)
    {
        if (OutsMax < OutsMin)
        {
            DR_Error("Fix3_AnnAmortSwap_t: max outs smaller than min outs!");
            goto RETURN;
        }
  
        deltaS = (OutsMax - OutsMin) / (NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates - 1);
        if (State == NULL)
        {
            DR_Error("Fix3_AnnAmortSwap_t: could not allocate memory for State!");
            goto RETURN;
        }

        for (s = 0; s < NbStates; s++)
        {
            State[s] = OutsMin + deltaS * s;
        }

        /* Find previous state var array */
        PState = *PrevState;

        if (PState != NULL)
        {
            /* PState exists, it is necessary to interpolate          */
            /* Prepare interpolation quadratic and linear polynomials */
            for (s = 0; s < NbStates-2; s++)
            {
                D[s][0] = 1. / ((PState[s]-PState[s+1])*(PState[s]-PState[s+2]));
                D[s][1] = 1. / ((PState[s+1]-PState[s])*(PState[s+1]-PState[s+2]));
                D[s][2] = 1. / ((PState[s+2]-PState[s])*(PState[s+2]-PState[s+1]));
            }
        }

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                AmortswapL[s] = Amortswap[s] + offset;
            }
            CollaretL = Collaret + offset;
            ZeroL     = Zero     + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {             
                for (s = 0; s < NbStates ; s++)
                {
                    /* The value of the loan with no outstanding is zero */
                    if (State[s] > TINY)
                    {
                        /* Determine new value of the state variable as a   */
                        /* function of index s and collaret payment         */
                        /* Determine indices for interpolation sj           */
                        /* Remember to remove discount factor from collaret */
                        /* if reset-in-advance                              */
                        
                        if (Arrears == 'Y')
                        { 
                            FltRate = CollaretL[i];
                        }
                        else
                        {
                            FltRate =  CollaretL[i] / ZeroL[i];
                        }

                        IntPmt   = MIN(State[s] * FltRate + Fee, FixCoupon);
                        PrincPmt = MIN(FixCoupon - IntPmt, State[s]);

                        NState   = State[s] - PrincPmt;
                        
                        /* Estimate amortswaper at NState */ 
                        if (FinalMatFlag)
                        {
                            if (Arrears == 'Y')
                            { 
                                IC[s] = NState * FinalWeight;
                            }
                            else
                            {
                                IC[s] = NState * ZeroL[i] * FinalWeight;
                            }
                        }
                        else if (PState != NULL)
                        {
                            
                            sj = (int) floor((NState - PState[0]) / (PState[1] - PState[0]));
                            
                            if (sj < 0)
                            {
                                /* New value outside of [0,Level].      */      
                                /* Use cap with amount state[0] */
                                
                                IC[s] = AmortswapL[0][i];
                            }
                            else if (sj >= (NbStates - 1))
                            {
                                /* New value outside of [0,Level].      */      
                                /* Use cap with amount state[NbStates - 1] */
                                
                                IC[s] = AmortswapL[NbStates - 1][i];
                            }
                            else
                            {
                                /* New value inside [0,Level].          */
                                /* Interpolate.                         */
                                
                                q = MIN(sj, NbStates - 3);
                                
                                sqinterp(PState[q], PState[q+1], PState[q+2],
                                         AmortswapL[q][i], AmortswapL[q+1][i], AmortswapL[q+2][i],
                                         D[q][0], D[q][1], D[q][2],
                                         NState,
                                         &(IC[s]));
                            }
                        }
                        else
                        {
                            /* Nothing to interpolate off */
                            IC[s] = 0.;
                            
                        } /* if then else (FinalMatFlag) */
                        
                        /* Loan interest payment */
                        if (Leg != 'F')
                        {
                            if (Arrears == 'Y')
                            {
                                IC[s] += IntPmt;
                            }
                            else
                            {
                                IC[s] += IntPmt * ZeroL[i];
                            }
                        }

                        /* Loan or funding principal payment */
                        if (Leg != 'S')
                        {
                            if (Arrears == 'Y')
                            {
                                IC[s] += PrincPmt;
                            }
                            else
                            {
                                IC[s] += PrincPmt * ZeroL[i];
                            }
                        }
                    }
                    else
                    { 
                        IC[s] = 0.;
                    }

                }  /* for s */
                
                /* Copy from intermadiate variable to amortswap */
                for (s = 0; s < NbStates ; s++)  
                {
                    AmortswapL[s][i] = IC[s];
                }
        
            }  /* for i */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    AmortswapL[s] = Amortswap[s] + offset;
                }
                CollaretL = Collaret + offset;
                ZeroL     = Zero     + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)
                    {
                        /* The value of the loan with no outstanding is zero */
                        if (State[s] > TINY)
                        {
                            /* Determine new value of the state variable as a   */
                            /* function of index s and collaret payment         */
                            /* Determine indices for interpolation sj           */
                            /* Remember to remove discount factor from collaret */
                            /* if reset-in-advance                              */
                            
                            if (Arrears == 'Y')
                            { 
                                FltRate = CollaretL[j];
                            }
                            else
                            {
                                FltRate =  CollaretL[j] / ZeroL[j];
                            }

                            IntPmt   = MIN(State[s] * FltRate + Fee, FixCoupon);
                            PrincPmt = MIN(FixCoupon - IntPmt, State[s]);

                            NState   = State[s] - PrincPmt;

                            /* Estimate amortswaper at NState */ 
                            if (FinalMatFlag)
                            {
                                if (Arrears == 'Y')
                                { 
                                    IC[s] = NState * FinalWeight;
                                }
                                else
                                {
                                    IC[s] = NState * ZeroL[j] * FinalWeight;
                                }
                            }
                            else if (PState != NULL)
                            {
                                
                                sj = (int) floor((NState - PState[0]) / (PState[1] - PState[0]));
                                
                                if (sj < 0)
                                {
                                    /* New value outside of [0,Level].      */      
                                    /* Use cap with amount state[0] */
                                
                                    IC[s] = AmortswapL[0][j];
                                }
                                else if (sj >= (NbStates - 1))
                                {
                                    /* New value outside of [0,Level].      */      
                                    /* Use cap with amount state[NbStates - 1] */
                                    
                                    IC[s] = AmortswapL[NbStates - 1][j];
                                }
                                else
                                {
                                    /* New value inside [0,Level].          */
                                    /* Interpolate.                         */
                                
                                    q = MIN(sj, NbStates - 3);
                                
                                    sqinterp(PState[q], PState[q+1], PState[q+2],
                                             AmortswapL[q][j], AmortswapL[q+1][j], AmortswapL[q+2][j],
                                             D[q][0], D[q][1], D[q][2],
                                             NState,
                                             &(IC[s]));
                                }
                            }
                            else
                            {
                                /* Nothing to interpolate off */
                                IC[s] = 0.;
                                
                            } /* if then else (FinalMatFlag) */
                                                                

                            if (Leg != 'F')
                            {
                                if (Arrears == 'Y')
                                {
                                    IC[s] += IntPmt;
                                }
                                else
                                {
                                    IC[s] += IntPmt * ZeroL[i];
                                }
                            }

                            if (Leg != 'S')
                            {
                                if (Arrears == 'Y')
                                {
                                    IC[s] += PrincPmt;
                                }
                                else
                                {
                                    IC[s] += PrincPmt * ZeroL[j];
                                }
                            }
                        }
                        else
                        { 
                            IC[s] = 0.;
                        }
                        
                    }  /* for s */
                    
                    /* Copy from intermadiate variable to amortswap */
                    for (s = 0; s < NbStates ; s++)  
                    {
                        AmortswapL[s][j] = IC[s];
                    }
                    
                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        AmortswapL[s] = Amortswap[s] + offset;
                    }
                    CollaretL = Collaret + offset;
                    ZeroL     = Zero     + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates ; s++)
                        {
                            /* The value of the loan with no outstanding is zero */
                            if (State[s] > TINY)
                            {
                                /* Determine new value of the state variable as a   */
                                /* function of index s and collaret payment         */
                                /* Determine indices for interpolation sj           */
                                /* Remember to remove discount factor from collaret */
                                /* if reset-in-advance                              */
                                
                                if (Arrears == 'Y')
                                { 
                                    FltRate = CollaretL[k];
                                }
                                else
                                {
                                    FltRate =  CollaretL[k] / ZeroL[k];
                                }

                                IntPmt   = MIN(State[s] * FltRate + Fee, FixCoupon);
                                PrincPmt = MIN(FixCoupon - IntPmt, State[s]);

                                NState   = State[s] - PrincPmt;
                                
                                if (FinalMatFlag)
                                {
                                    if (Arrears == 'Y')
                                    { 
                                        IC[s] = NState * FinalWeight;
                                    }
                                    else
                                    {
                                        IC[s] = NState * ZeroL[k] * FinalWeight;
                                    }
                                }
                                else if (PState != NULL)
                                {
                                    
                                    sj = (int) floor((NState - PState[0]) / (PState[1] - PState[0]));
                                    
                                    if (sj < 0)
                                    {
                                        /* New value outside of [0,Level].      */      
                                        /* Use cap with amount state[0] */
                                        
                                        IC[s] = AmortswapL[0][k];
                                    }
                                    else if (sj >= (NbStates - 1))
                                    {
                                        /* New value outside of [0,Level].      */      
                                        /* Use cap with amount state[NbStates - 1] */
                                        
                                        IC[s] = AmortswapL[NbStates - 1][k];
                                    }
                                    else
                                    {
                                        /* New value inside [0,Level].          */
                                        /* Interpolate.                         */
                                        
                                        q = MIN(sj, NbStates - 3);
                                        
                                        sqinterp(PState[q], PState[q+1], PState[q+2],
                                                 AmortswapL[q][k], AmortswapL[q+1][k], AmortswapL[q+2][k],
                                                 D[q][0], D[q][1], D[q][2],
                                                 NState,
                                                 &(IC[s]));
                                    }
                                }
                                else
                                {
                                    /* Nothing to interpolate off */
                                    IC[s] = 0.;
                                    
                                } /* if then else (FinalMatFlag) */
                            
                                if (Leg != 'F')
                                {
                                    if (Arrears == 'Y')
                                    {
                                        IC[s] += IntPmt;
                                    }
                                    else
                                    {
                                        IC[s] += IntPmt * ZeroL[i];
                                    }
                                }

                                if (Leg != 'S')
                                {
                                    if (Arrears == 'Y')
                                    {
                                        IC[s] += PrincPmt;
                                    }
                                    else
                                    {
                                        IC[s] += PrincPmt * ZeroL[k];
                                    }
                                }
                            }
                            else
                            { 
                                IC[s] = 0.;
                            }
                            
                        }  /* for s */
                        
                        /* Copy from intermadiate variable to amortswap */
                        for (s = 0; s < NbStates ; s++)  
                        {
                            AmortswapL[s][k] = IC[s];
                        }
                        
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        
        /* Free PState and assign State to PrevSate */        
        if (Free_DR_Array (PState, DOUBLE, 0, NbStates - 1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevState = State;
        State = NULL;
    }  /* if CollaretFlag*/

    /* Funding pmt is added on reset (i.e. pmt) days for arrears,  */
    /* but it is added on loan (i.e. collaret) days for advance    */

    if ((Leg != 'L') &&
        (((FloaterFlag)  && (Arrears == 'Y')) ||
         ((CollaretFlag) && (Arrears == 'N')) ))
    {
        /* If funding more frequent than loan, use last outstanding */
        /* info (remembered in PrevState)                           */

        State = *PrevState;

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                AmortswapL[s] = Amortswap[s] + offset;
            }
            FloaterL = Floater + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {             
                for (s = 0; s < NbStates ; s++)  
                {
                    AmortswapL[s][i] += FloaterL[i] * State[s];
                }
        
            }  /* for i */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    AmortswapL[s] = Amortswap[s] + offset;
                }
                FloaterL = Floater + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    /* Copy from intermadiate variable to amortswap */
                    for (s = 0; s < NbStates ; s++)  
                    {
                        AmortswapL[s][j] += FloaterL[j] * State[s];
                    }
                    
                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        AmortswapL[s] = Amortswap[s] + offset;
                    }
                    FloaterL = Floater + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        /* Copy from intermadiate variable to amortswap */
                        for (s = 0; s < NbStates ; s++)  
                        {
                            AmortswapL[s][k] += FloaterL[k] * State[s];
                        }
                        
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
    }  /* if FloaterFlag*/

    status = SUCCESS;
    
    RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates - 1);
    }

    return (status);

}  /* Fix3_AnnAmortSwap_t */
