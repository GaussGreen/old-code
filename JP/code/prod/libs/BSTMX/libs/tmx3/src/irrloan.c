/**************************************************************************/
/*      Payoff function for the irr loan swap                             */
/**************************************************************************/
/*      irrloan.c                                                         */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"



/*****  IRRLoanSwap_t  **************************************************/
/*
*   Calculates the complex swap price
*   no dev's, just add complex and floating coupons
*/
int  IRRLoanSwap_t
             (double    **Complex,        /* (O)Prices for all states     */
              long        ResetCx,        /* (I)Reset loan notionals      */
              long        PayFund,        /* (I)Pay funding               */
              double     *IndexCx,        /* (I)Index level               */
              double     *Funding,        /* (I)Funding                   */
              double     *ZeroToPmtCx,    /* (I)PmtZero 4 this reset      */
              double      Annuity,        /* (I)Annuity rate              */
              double      FloorCx,        /* (I)Floor spread              */
              double      CapCx,          /* (I)Cap spread                */
              double      UpRateCx,       /* (I)Complex up rate           */
              double      DcfCx,          /* (I)Day count fract           */
              double      DcfOut,         /* (I)Day count fract for Outst.*/
              double      Notional,       /* (I)Notional of pmt           */
              char        CompCx,         /* (I)Simple/Compound pmt       */
              char        Freq,           /* (I)Freq of annuity           */
              int         Tenor,          /* (I)Nb of periods in rem ann  */
              int         NbStates,       /* (I)Nb of states              */
              double      MaxState,       /* (I)Upper state bound         */
              double      MinState,       /* (I)Lower state bound         */
              double    **PrevState,      /* (I)Prev state levels         */
              int         t,              /* (I)Current time point        */
              TREE_DATA   *tree_data)     /* (I)Tree data structure       */
{
    /* Local slice pointers */

    double  *ComplexL[MAXNBSTATES+1];
    double  *IndexCxL;
    double  *FundingL;
    double  *ZeroToPmtCxL;

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

    double  Payoff[MAXNBSTATES];      /* Intermediate values of complex   */
    double  tmpPayoff;
    double  IRRStrike[MAXNBSTATES];
    double  PmtRate,PmtDcf,AccFact;
    double  NOuts;
    int     isSimpleCx;
    int     FreqN;
    
    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    if (*PrevState == NULL)
    {
        DR_Error("IRRLoanSwap_t: No previous states !");
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
        DR_Error("IRRLoanSwap_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }


    /* COMPLEX LEG */

    isSimpleCx = (CompCx == 'S');
    FreqN = Conv_Freq(Freq);

    /* Add complex payoff if it's a reset */
    if (ResetCx)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("IRRLoanSwap_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("IRRLoanSwap_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s=0; s<NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */

        PState = *PrevState;

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

        /* Note: outs is 1 for zero coupon */
        for (s=0; s<NbStates; s++)
        {
            if (AnnIRR(State[s],Annuity,Tenor,FreqN,&(IRRStrike[s])) == FAILURE)
            {
                DR_Error("IRRLoanSwap: Strike not found!");
                goto RETURN;
            }
        }



        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                ComplexL[s] = Complex[s] + offset;
            }
            IndexCxL  = IndexCx  + offset;
            ZeroToPmtCxL = ZeroToPmtCx + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                /* ---------------------------- */
                /*   Rearrange variables        */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    /* ------------------------------------- */
                    /*   Find next coupon and outs           */
                    /* ------------------------------------- */
                    
                    PmtRate = COLLAR(IndexCxL[i]+UpRateCx,CapCx,FloorCx);
                    PmtRate = MIN(PmtRate,IRRStrike[s]);
                    AccFact = 1 + PmtRate * DcfOut;
                    NOuts   = MAX(State[s] * AccFact - Annuity/FreqN,0.);

                    /* ------------------------------------- */
                    /*   Payoff += interpolated Complex price */
                    /* ------------------------------------- */
                    PmtDcf    = ACC_FN(PmtRate,DcfCx,isSimpleCx);
                    Payoff[s] = PmtDcf * ZeroToPmtCxL[i] * Notional * State[s];

                    if (IS_EQUAL(PState[0], PState[1]))
                    {
                        Payoff[s] += ComplexL[0][i];
                    }
                    else
                    {
                        InterpStateLevel = NOuts;
                        sj = (int) floor((InterpStateLevel - PState[0]) /
                                         (PState[1] - PState[0]));

                        if (sj >= (NbStates - 1)) /* to the right */
                        {
                            Payoff[s] += ComplexL[NbStates-1][i];
                        }
                        else if (sj < 0) /* to the left */
                        {
                            Payoff[s] += ComplexL[0][i];
                        }
                        else /* in between */
                        {
                            q = MIN(sj, NbStates - 3);
                    
                            sqinterp(PState[q], PState[q+1], PState[q+2],
                                     ComplexL[q][i],ComplexL[q+1][i],
                                     ComplexL[q+2][i],
                                     D[q][0], D[q][1], D[q][2],
                                     InterpStateLevel,
                                     &(tmpPayoff));
                            Payoff[s] += tmpPayoff;
                           
                        } /* if sj */
                    }

                } /* for state idx */

                /* replace complex prices with new values for each state */

                for (s=0; s<NbStates ; s++)  
                {
                    ComplexL[s][i] = Payoff[s];
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
                    ComplexL[s] = Complex[s] + offset;
                }
                IndexCxL  = IndexCx  + offset;
                ZeroToPmtCxL = ZeroToPmtCx + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {                
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ------------------------------------- */
                        /*   Find next coupon level              */
                        /* ------------------------------------- */
                        PmtRate = COLLAR(IndexCxL[j]+UpRateCx,CapCx,FloorCx);
                        PmtRate = MIN(PmtRate,IRRStrike[s]);
                        AccFact = 1 + PmtRate * DcfOut;
                        NOuts   = MAX(State[s] * AccFact - Annuity/FreqN,0.);

                        /* ------------------------------------- */
                        /*   Payoff += interpolated Complex price */
                        /* ------------------------------------- */
                        PmtDcf    = ACC_FN(PmtRate,DcfCx,isSimpleCx);
                        Payoff[s] = PmtDcf*ZeroToPmtCxL[j]*Notional*State[s];
                    
                        if (IS_EQUAL(PState[0], PState[1]))
                        {
                            Payoff[s] += ComplexL[0][j];
                        }
                        else
                        {
                            InterpStateLevel = NOuts;
                            sj = (int)floor((InterpStateLevel - PState[0])/
                                            (PState[1] - PState[0]));
    
                            if (sj >= (NbStates - 1)) /* to the right */
                            {
                                Payoff[s] += ComplexL[NbStates-1][j];
                            }
                            else if (sj < 0) /* to the left */
                            {
                                Payoff[s] += ComplexL[0][j];
                            }
                            else /* in between */
                            {
                                q = MIN(sj, NbStates - 3);
                        
                                sqinterp(PState[q], PState[q+1], 
                                         PState[q+2],
                                         ComplexL[q][j],ComplexL[q+1][j],
                                         ComplexL[q+2][j],
                                         D[q][0], D[q][1], D[q][2],
                                         InterpStateLevel,
                                         &(tmpPayoff));
                                Payoff[s] += tmpPayoff;
                               
                            } /* if sj */
                        }
    
                    } /* for state idx */
    
                    /* replace complex prices with new values for each state */
    
                    for (s=0; s<NbStates ; s++)  
                    {
                        ComplexL[s][j] = Payoff[s];
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
                        ComplexL[s] = Complex[s] + offset;
                    }
                    IndexCxL  = IndexCx  + offset;
                    ZeroToPmtCxL = ZeroToPmtCx + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ------------------------------------- */
                            /*   Find next coupon level              */
                            /* ------------------------------------- */
                            PmtRate = COLLAR(IndexCxL[k]+UpRateCx,CapCx,FloorCx);
                            PmtRate = MIN(PmtRate,IRRStrike[s]);
                            AccFact = 1 + PmtRate * DcfOut;
                            NOuts   = MAX(State[s] * AccFact - Annuity/FreqN,0.);

                            /* ------------------------------------- */
                            /*   Payoff += interpolated Complex price */
                            /* ------------------------------------- */
                            PmtDcf    = ACC_FN(PmtRate,DcfCx,isSimpleCx);
                            Payoff[s] = PmtDcf*ZeroToPmtCxL[k]*Notional*State[s];
                    
                            if (IS_EQUAL(PState[0], PState[1]))
                            {
                                Payoff[s] += ComplexL[0][k];
                            }
                            else
                            {
                                InterpStateLevel = NOuts;
                                sj = (int) floor(
                                            (InterpStateLevel - PState[0])/
                                            (PState[1] - PState[0]));
        
                                if (sj >= (NbStates - 1)) /* right */
                                {
                                    Payoff[s] += ComplexL[NbStates-1][k];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += ComplexL[0][k];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                            
                                    sqinterp(PState[q], PState[q+1], 
                                             PState[q+2],
                                             ComplexL[q][k],ComplexL[q+1][k],
                                             ComplexL[q+2][k],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));
                                    Payoff[s] += tmpPayoff;
                                   
                                } /* if sj */
                            }
        
                        } /* for state idx */
        
                        /* replace complex prices for each state */
        
                        for (s=0; s<NbStates ; s++)  
                        {
                            ComplexL[s][k] = Payoff[s];
                        }

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

        /* -------------------------------------- */
        /*   replace PrevState with new states   */
        /* -------------------------------------- */

        if (Free_DR_Array (PState, DOUBLE, 0, NbStates-1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevState = State;
        State = NULL;

    }  /* if ResetComplex */


    /* FUNDING LEG */

    /* Add floating payoff if it's a reset */
    if (PayFund)
    {

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                ComplexL[s] = Complex[s] + offset;
            }
            FundingL = Funding + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates ; s++)  
                {
                    ComplexL[s][i] -= FundingL[i] * (*PrevState)[s];
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
                    ComplexL[s] = Complex[s] + offset;
                }
                FundingL = Funding + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        ComplexL[s][j] -= FundingL[j] * (*PrevState)[s];
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
                        ComplexL[s] = Complex[s] + offset;
                    }
                    FundingL = Funding + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            ComplexL[s][k] -= FundingL[k] * (*PrevState)[s];
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

}  /* IRRLoanSwap_t */
