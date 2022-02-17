/****************************************************************************/
/*      Calculation routine customised for a sticky turbo swap              */
/****************************************************************************/
/*      STICKY.C                                                             */
/****************************************************************************/


/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

/********** Hyb3_StickyTurboSwap_t ***************************************/
/** Adds currently reset StickyTurbo  and Funding Pmts (if any)       
    and sorts out previous coupon as uncertainty is resolved.         
    NO DEV.                                                           
 ********************************************************************/

int Hyb3_StickyTurboSwap_t(
        TSLICE      *StickyTurbo,           /**< (I/O) prices for all states   */
        /* TURBO DATA */
        int         ArrearsSticky,          /**< (I)ArrearsForAll              */
        int         ResetFlagStickyTurbo,   /**< (I)Reset for sticky           */
        TSLICE      SwapPmt,                /**< (I)Turbo SwapPmt              */
        TSLICE      ZeroToStickyPmt,        /**< (I)domestic zero to sticky pmt*/
        double      FloorSpreadAmt,         /**< (I) floor spread dom. Amount  */
        double      CapSpreadAmt,           /**< (I) cap spread dom Amount     */
        double      LifeCapAmt,             /**< (I) Life Cap dom.amount       */
        double      LifeFloorAmt,           /**< (I) Life Floor dom.amount     */

        /* FUNDING LEG DATA */
        int         FundResetFlag,
        TSLICE      FundPmt,                /**< (I) Funding pmt denom in dom  */
        
        /* STICKY DATA */
        char        StickyCoAoF,            /**< (I) sticky 'C'oupon,'A'djustable,
                                                  or sticky 'F'ormula           */
        int         IsStickyOff,            /**< (I)                            */
        int         NbStates,               /**< (I) Nb of states for state var.*/
        double      MaxState,               /**< (I) U.Bound for state var      */
        double      MinState,               /**< (I) L.Bound "  "  "            */
        double      **PrevStates,           /**< (I/O) array of previous states */

        /* HYB3_TREE_DATA  */
        int         t,
        HYB3_TREE_DATA   *tree_data)
{
    /* Local slice pointers */

    TSLICE  SwapPmtL,FundPmtL;
    TSLICE  StickyTurboL[MAXNBSTATES+1];
    TSLICE  ZeroToStickyPmtL;

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
    double  StickyPmt;
    double  Payoff[MAXNBSTATES];      /* Intermediate values of sticky   */
    double  tmpPayoff;
    double  CapSpdStrike[MAXNBSTATES];
    double  FloorSpdStrike[MAXNBSTATES];
    int     IsLastReset;              /* true=this is last reset date    */
    double  Q;
    
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
        DR_Error("Hyb3_Sticky_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }


    if (StickyCoAoF != 'C' && StickyCoAoF != 'F' &&
        StickyCoAoF != 'A')
    {
        DR_Error ("StickyCoAoF should be either C or A or F!\n");
        goto RETURN;
    }

    /* Funding LEG */

    /* Add floating payoff if it's a reset */
    if (FundResetFlag)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyTurboL[s] = StickyTurbo[s] + offset;
                }

                FundPmtL  = FundPmt  + offset;
                
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {                    
                    for (s = 0; s < NbStates ; s++)  
                    {
                        StickyTurboL[s][k] += FundPmtL[k];  /*ADD Fund Pmt!!!*/
                    }/* for s */
                }  /* for k */
            }  /* for j */
        }/* for i */
 
    }  /* if FundResetFlag */

    /* STICKY LEG */

    /* Add sticky payoff if it's a reset */
    if (ResetFlagStickyTurbo)
    {
        /* ---------------------------------------- */
        /*   Prepare state levels for this time pt  */
        /*   and the cap/floor strikes              */
        /* ---------------------------------------- */

        if (MaxState < MinState)
        {
            DR_Error("Hyb3_Sticky_t: max state smaller than min state!");
            goto RETURN;
        }

        deltaS = (MaxState - MinState)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates-1);
        if (State == NULL)
        {
            DR_Error("Hyb3_Sticky_t: could not allocate memory for State!");
            goto RETURN;
        }
        
        for (s = 0; s < NbStates; s++)
        {
            State[s] = MinState + deltaS * s;
        }

        for (s = 0; s < NbStates; s++)
        {
            CapSpdStrike[s] = State[s] + CapSpreadAmt;
            FloorSpdStrike[s] = State[s] + FloorSpreadAmt;
        }

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */

        if (!IsLastReset)
        {
            PState = *PrevStates;

            if (!IS_EQUAL(PState[0] , PState[NbStates - 1]))
            {
                for (s = 0; s < NbStates - 2; s++)
                {
                    D[s][0] = 1. / ((PState[s]-PState[s+1])
                                   *(PState[s]-PState[s+2]));
                    D[s][1] = 1. / ((PState[s+1]-PState[s])
                                   *(PState[s+1]-PState[s+2]));
                    D[s][2] = 1. / ((PState[s+2]-PState[s])
                                   *(PState[s+2]-PState[s+1]));
                }/* for s */
            }/* if more than 1 state level */
        }/* if !IsLastReset */

        for (i = Bottom1; i <= Top1; i ++)
        {       
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyTurboL[s] = StickyTurbo[s] + offset;
                }

                SwapPmtL  = SwapPmt  + offset;
                ZeroToStickyPmtL = ZeroToStickyPmt + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ---------------------------- */
                        /*   Add the new coupon amt     */
                        /* ---------------------------- */
                        
                        Q = (IsStickyOff) ? SwapPmtL[k] :  
                                            COLLAR( SwapPmtL[k],
                                                    CapSpdStrike[s],
                                                    FloorSpdStrike[s]);
                        StickyPmt = COLLAR(Q,LifeCapAmt,LifeFloorAmt);
                        if (ArrearsSticky)
                        {
                            Payoff[s] = StickyPmt;
                        }
                        else
                        {
                            Payoff[s] = StickyPmt * ZeroToStickyPmtL[k];
                        }
                            
                        /* ------------------------------------- */
                        /*   Add the interpolated Sticky price   */
                        /*   to the payoff if appropriate        */
                        /* ------------------------------------- */
        
                        if (!IsLastReset)
                        {
                            if (IS_EQUAL(PState[1] , PState[0]))
                            {
                                Payoff[s] += StickyTurboL[0][k];
                            }
                            else
                            {   
                                switch (StickyCoAoF)
                                {
                                    case 'C':
                                        InterpStateLevel = StickyPmt;
                                        break;
                                    case 'A':
                                        InterpStateLevel = SwapPmtL[k];
                                        break;
                                    case 'F':
                                        InterpStateLevel = Q;
                                        break;
                                    default:/* does not get there because of earlier check*/
                                            /* it is here to avoid warnings*/
                                    {
                                        DR_Error("StickyCoAoF should be either"
                                                   "C or A or F\n");
                                        goto RETURN;
                                    }
                                                                                              
                                }/* switch */

                                sj = (int) floor(
                                            (InterpStateLevel - PState[0])/
                                            (PState[1] - PState[0]));
                                if (sj >= (NbStates - 1)) /* right */
                                {
                                    Payoff[s] += StickyTurboL[NbStates-1][k];
                                }
                                else if (sj < 0) /* to the left */
                                {
                                    Payoff[s] += StickyTurboL[0][k];
                                }
                                else /* in between */
                                {
                                    q = MIN(sj, NbStates - 3);
                            
                                    sqinterp(PState[q], PState[q+1],
                                             PState[q+2],
                                             StickyTurboL[q][k],
                                             StickyTurboL[q+1][k],
                                             StickyTurboL[q+2][k],
                                             D[q][0], D[q][1], D[q][2],
                                             InterpStateLevel,
                                             &(tmpPayoff));

                                    Payoff[s] += tmpPayoff;
                                } /* if sj */
                            }    
                        }
                        else
                        {
                                Payoff[s] += StickyTurboL[0][k];
                        }/* if !IsLastReset, else */
                    } /* for state idx */

                    /* replace sticky prices for each state */
                    for (s = 0; s < NbStates ; s++)
                    {
                            StickyTurboL[s][k] = Payoff[s];
                    }/* for s*/
                }  /* for k */
            }/* for j*/
        }  /* for i*/
        

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
        Free_DR_Array (State,DOUBLE,0,NbStates - 1);
    }

    return (status);

}  /* Hyb3_StickyTurboSwap_t */

