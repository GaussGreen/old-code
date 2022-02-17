/****************************************************************************/
/*      Callable ladder:                                                    */
/*      Calculation going backward in the tree.                             */
/****************************************************************************/
/*      CALC39.C                                                            */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fix123head.h"
#include "template39.h"
#include "evtstat.h"

/*****  Fix3_Sticky_LadderSwap_t ******************************************/
/**
*   Calculates the ladder swap price
*   no dev's, just add ladder and floating coupons
*   or just revise "coupon rate" states.
*
*   Important: For RIB reset-in-arrears, the coupon rate states are
*              RIBed rates; while 
*              for RIB reset-in-advance, the coupon rate states are
*              unRIBed rates. 
*   New:
*      ResetFlagSt is now "overloaded" to handle the RIB Ladder
*      Reset-in-Advance case:
*      ResetFlagSt = 1 --- RIB Ladder Reset-in-Arrears/nonRIB reset time pt
*                          (Add rib coupons)
*      ResetFlagSt = 2 --- RIB Ladder Reset-in-Advance payment & reset time pt
*                          (Revise "coupon rate" states and Add rib coupons)
*      ResetFlagSt = 3 --- RIB Ladder Reset-in-Advance first reset time pt 
*                          (Revise "coupon rate" states)
*      ResetFlagSt = 4 --- RIB Ladder Reset-in-Advance last payment time pt 
*                          (set the ladder leg to the last payment amount)
*                        
*/
static int  Fix3_Sticky_LadderSwap_t
             (double    **Ladder,         /**< (O) Prices for all states   */
              long        ResetFlagSt,    /**< (I) Reset/Payment at timept */
              long        ResetFlagFund,  /**< (I) Reset at timept         */
              double     *Step,           /**< (I) Step for ladder         */
              double     *Spread,         /**  (I) Spread for ladder       */
              char        AoM,            /**  (I) Additive or mult step   */
              double     *Funding,        /**< (I) Funding                 */
              double     *ZeroToPmtSt,    /**< (I) PmtZero 4 this reset    */
              double      FloorSt,        /**< (I) Floor spread            */
              double      CapSt,          /**< (I) Cap spread              */
              double      StickyCoeff,    /**< (I) Sticky coefficient      */
              double      RibFrac,       /**< (I) Range accrual weight    */
              double      OutsSt,         /**< (I) Outstanding             */
              double      DcfSt,          /**< (I) Day count fract         */
              char        SoZ,            /**< (I) Swap/zero coupon        */
              char        CompSt,         /**< (I) Simple/Compound pmt     */
              int         NbStates,       /**< (I) Nb of states            */
              double     *CurrStates,     /**< (I) Curr state levels at i-1*/
              double     *PrevStates,     /**< (I) Prev state levels at i  */
              int         t,              /**< (I) Current time point      */
              FIX3_TREE_DATA   *tree_data)/**< (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *LadderL[MAXNBSTATES+1];
    double  *FundingL;
    double  *ZeroToPmtStL;
    double  *StepL;
    double  *SpreadL;

    /* state-variable variables */
    int     s;                        /* State variable index            */

    /* Payoff variables */
    double  tmpLadder[MAXNBSTATES];   /* Intermediate values of ladder   */
    double  PRate[MAXNBSTATES];
    double  Cpn[MAXNBSTATES];

    double  X;
    int     isZeroSt, isSimpleSt;

    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int addRibCoupon;                 /* Internal flag */
    int resetRibRate;                 /* Internal flag */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */

    long    CurrDate = tree_data->TPDate[t];


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    addRibCoupon = (ResetFlagSt && (ResetFlagSt != 3));
    resetRibRate = (ResetFlagSt && (ResetFlagSt != 4));

    if (addRibCoupon) {
        /* STICKY LEG */

        isZeroSt   = (SoZ == 'Z');
        isSimpleSt = (CompSt == 'S');
    }

    /* Add ladder payoff or perform state variable revision */
    if (ResetFlagSt)
    {
        if (PrevStates == NULL || CurrStates == NULL)
        {
            DR_Error("Fix3_Sticky_LadderSwap_t: No previous or current "
                     "states on %ld! ResetFlagSt=%d",
                     CurrDate, ResetFlagSt);
            goto RETURN;
        }

        if (ResetFlagSt == 1) 
        {
            for (s=0; s<NbStates; s++)
            {
                PRate[s] =  ACC_FN(PrevStates[s],DcfSt,isSimpleSt);
                Cpn[s]   =  PRate[s]*OutsSt;
            }

        }
        else if (ResetFlagSt == 2 || ResetFlagSt == 4) 
        {
            for (s=0; s<NbStates; s++)
            {
                PRate[s] =  ACC_FN(CurrStates[s]*RibFrac,DcfSt,isSimpleSt);
                Cpn[s]   =  PRate[s]*OutsSt;
            }

        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                LadderL[s] = Ladder[s] + offset;
            }

            if (resetRibRate) 
            {
                StepL  = Step  + offset;
                SpreadL = Spread + offset;
            }

            if (ResetFlagSt == 1)
                ZeroToPmtStL = ZeroToPmtSt + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {
                if (ResetFlagSt == 4)
                {
                    for (s = 0; s < NbStates; s++)
                        LadderL[s][i] = Cpn[s];

                    continue;
                }

                /* ---------------------------- */
                /*   Update current value:      */
                /*   mulitply and add coupon    */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    tmpLadder[s] = LadderL[s][i];

                    if (ResetFlagSt == 1)
                    {
                        if (isZeroSt)
                            tmpLadder[s] = LadderL[s][i] * (1.+PRate[s]);

                        tmpLadder[s] += ZeroToPmtStL[i]*Cpn[s];
                    }
                }

                /* ---------------------------- */
                /*   Rearrange variables        */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    /* ------------------------------------- */
                    /*   Find next coupon level              */
                    /* ------------------------------------- */
                    
                    if (ResetFlagSt == 2)
                        X  = StickyCoeff * CurrStates[s] * RibFrac;
                    else 
                        X  = StickyCoeff * CurrStates[s];

                    if ( AoM == 'A')
                    {
                        X += SpreadL[i];
                        X += StepL[i];
                        X  = COLLAR(X, CapSt, FloorSt);
                    }
                    else if ( AoM == 'B')
                    {
                        X += SpreadL[i] * StepL[i];
                        X  = COLLAR(X, CapSt, FloorSt);
                    }
                    else
                    {
                        X += SpreadL[i];
                        X  = COLLAR(X,CapSt,FloorSt);
                        X *= StepL[i];
                    }

                    if (ResetFlagSt == 1)
                        X *= RibFrac;                
  
                    /* ------------------------------------- */
                    /*   Payoff = interpolated Ladder price  */
                    /* ------------------------------------- */
                    if (Fix3_DoubleQuadraticInterp(PrevStates, 
                                              tmpLadder, 
                                              NbStates,
                                              X, 
                                              &(LadderL[s][i])) == FAILURE)
                        goto RETURN;

                    if (ResetFlagSt == 2) {
                        if (isZeroSt)
                            LadderL[s][i] = LadderL[s][i] * (1.+PRate[s]);

                        LadderL[s][i] += Cpn[s];
                    }

                } /* for state idx */

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
                    LadderL[s] = Ladder[s] + offset;
                }

                if (ResetFlagSt == 1)
                    ZeroToPmtStL = ZeroToPmtSt + offset;

                if (resetRibRate) {
                    StepL  = Step  + offset;
                    SpreadL = Spread + offset;
                }

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {

                    if (ResetFlagSt == 4)
                    {
                        for (s = 0; s < NbStates; s++)
                            LadderL[s][j] = Cpn[s];

                        continue;
                    }

                    /* ---------------------------- */
                    /*   Update current value:      */
                    /*   mulitply and add coupon    */
                    /* ---------------------------- */
                    for (s = 0; s < NbStates; s++)
                    {
                        tmpLadder[s] = LadderL[s][j];

                        if (ResetFlagSt == 1)
                        {
                            if (isZeroSt)
                                tmpLadder[s] = LadderL[s][j] * (1.+PRate[s]);

                            tmpLadder[s] += ZeroToPmtStL[j]*Cpn[s];
                        }
                    }
                
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ------------------------------------- */
                        /*   Find next coupon level              */
                        /* ------------------------------------- */
                    
                        if (ResetFlagSt == 2)
                            X  = StickyCoeff * CurrStates[s] * RibFrac;
                        else 
                            X  = StickyCoeff * CurrStates[s];

                        if ( AoM == 'A')
                        {
                            X += SpreadL[j];
                            X += StepL[j];
                            X  = COLLAR(X,CapSt,FloorSt);
                        }
                        else if ( AoM == 'B')                    
                        {                        
                            X += SpreadL[j] * StepL[j];                        
                            X  = COLLAR(X, CapSt, FloorSt);                    
                        }
                        else
                        {
                            X += SpreadL[j];
                            X  = COLLAR(X,CapSt,FloorSt);
                            X *= StepL[j];
                        }


                        if (ResetFlagSt == 1)
                            X *= RibFrac;                             

                        /* ------------------------------------- */
                        /*   Payoff = interpolated Ladder price  */
                        /* ------------------------------------- */
        
                        if (Fix3_DoubleQuadraticInterp(PrevStates, 
                                                  tmpLadder, 
                                                  NbStates,
                                                  X, 
                                                  &(LadderL[s][j])) == FAILURE)
                            goto RETURN;
    
                        if (ResetFlagSt == 2) {
                            if (isZeroSt)
                                LadderL[s][j] = LadderL[s][j] * (1.+PRate[s]);

                            LadderL[s][j] += Cpn[s];
                        }

                    } /* for state idx */
    
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
                        LadderL[s] = Ladder[s] + offset;
                    }

                    if (ResetFlagSt == 1)
                        ZeroToPmtStL = ZeroToPmtSt + offset;

                    if (resetRibRate) {
                        StepL  = Step  + offset;
                        SpreadL = Spread + offset;
                    }

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        if (ResetFlagSt == 4)
                        {
                            for (s = 0; s < NbStates; s++)
                                LadderL[s][k] = Cpn[s];

                            continue;
                        }

                        /* ---------------------------- */
                        /*   Update current value:      */
                        /*   mulitply and add coupon    */
                        /* ---------------------------- */
                        for (s = 0; s < NbStates; s++)
                        {
                            tmpLadder[s] = LadderL[s][k];


                            if (ResetFlagSt == 1)
                            {
                                if (isZeroSt)
                                    tmpLadder[s] = LadderL[s][k]*(1.+PRate[s]);

                                tmpLadder[s] += ZeroToPmtStL[k]*Cpn[s];
                            }
                        }
        
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ------------------------------------- */
                            /*   Find next coupon level              */
                            /* ------------------------------------- */
                    
                            if (ResetFlagSt == 2)
                                X  = StickyCoeff * CurrStates[s] * RibFrac;
                            else 
                                X  = StickyCoeff * CurrStates[s];

                            if ( AoM == 'A')
                            {
                                X += SpreadL[k];
                                X += StepL[k];
                                X  = COLLAR(X,CapSt,FloorSt);
                            } 
                            else if ( AoM == 'B')                    
                            {                        
                                X += SpreadL[k] * StepL[k];                        
                                X  = COLLAR(X, CapSt, FloorSt);                    
                            }
                            else
                            { 
                                X += SpreadL[k];
                                X  = COLLAR(X,CapSt,FloorSt);
                                X *= StepL[k];
                            }
                           

                            if (ResetFlagSt == 1)
                                X *= RibFrac;                             

                            /* ------------------------------------- */
                            /*   Payoff = interpolated Ladder price  */
                            /* ------------------------------------- */
            
                            Fix3_DoubleQuadraticInterp(PrevStates, 
                                                  tmpLadder, 
                                                  NbStates,
                                                  X, 
                                                  &(LadderL[s][k]));
        
                            if (ResetFlagSt == 2) {
                                if (isZeroSt)
                                    LadderL[s][k] = LadderL[s][k]*(1.+PRate[s]);

                                LadderL[s][k] += Cpn[s];
                            }

                        } /* for state idx */

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

    }  /* if ResetFlagLadder */

    /* FUNDING LEG */

    /* Add floating payoff if it's a reset BUT ONLY if it's ladder date */
    /* This restriction is achieved by setting flags in calc properly   */
    if (ResetFlagFund)
    {

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                LadderL[s] = Ladder[s] + offset;
            }
            FundingL = Funding + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates ; s++)  
                {
                    LadderL[s][i] -= FundingL[i];
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
                    LadderL[s] = Ladder[s] + offset;
                }
                FundingL = Funding + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        LadderL[s][j] -= FundingL[j];
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
                        LadderL[s] = Ladder[s] + offset;
                    }
                    FundingL = Funding + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            LadderL[s][k] -= FundingL[k];
                        }
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlagFloat */

    status = SUCCESS;
    
RETURN:

    return (status);

}  /* Fix3_Sticky_LadderSwap_t */


/*****  Fix3_RibLadderInterpolate ***************************************/
/*
*   Interpolates the Rib Ladder swap in the Rib dimension;
*   no DEV and no coupons are added -- interpolation only
*   Note: 
*         After this interpolation routine, the size of the Rib dimension
*         (NbRibActive) of the slice matrix
*         decreases by at least 1 (since at least 1 new RIB observation is
*         processed and subsequently deactivated), so we only update
*         in the Rib dimension up to (NbRibActive - 1). 
*/
static int  Fix3_RibLadderInterpolate_t
             (double***             RibLadder,  /* (I/O) State variable slices */
              int                   NbStates,   /* (I) Nb of ladder states     */
              int                   NbRibActive,/* (I) Nb of active Rib slices */     
              double const*         ObsIndex,   /* (I) 1 inside, 0 outside     */
              int                   t,          /* (I) Current time point      */
              FIX3_TREE_DATA const* tree_data)  /* (I) Tree data structure     */
{
    
    /* Local slice pointers */
    double        *SliceL = NULL;
    double  const *ObsIndexL = NULL;

    /* Arrays of Rib states */
    double  RibInRangeCount[MAXNBDATE];
    double  tmpInValues    [MAXNBDATE];
    double  tmpOutValues   [MAXNBDATE];

    /* State-variable variables */
    int     s;      /* State variable index */
    int     r;      /* Rib obs index        */

    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */



    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    for (i = 0; i < NbRibActive; i++)
    {
        RibInRangeCount[i] = (double) i;
        tmpInValues    [i] = 0.;
        tmpOutValues   [i] = 0.;
    }


    /************************  1 FACTOR    ****************************/

    if (tree_data->NbFactor == 1)
    {
        offset    = Fix3_Node_Offset(1, 0, 0, t, tree_data);
        ObsIndexL = ObsIndex + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {                
            for (s = 0; s < NbStates; s++)
            {
                /* interpolate in the Rib dimension */
                for (r = 0; r < NbRibActive; r++)
                {
                    SliceL = RibLadder[r][s] + offset;
                    tmpInValues[r] = SliceL[i];
                }
                for (r = 0; r < NbRibActive-1; r++)
                {
                    if (Fix3_DoubleQuadraticInterp (
                            RibInRangeCount,
                            tmpInValues,
                            NbRibActive,
                            r + ObsIndexL[i],
                            &(tmpOutValues[r])) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
                
                /* store results */
                for (r = 0; r < NbRibActive-1; r++)
                {
                    SliceL = RibLadder[r][s] + offset;
                    SliceL[i] = tmpOutValues[r];
                }

            } /* for state idx */

        } /* for i */

    } /* if NbFactor == 1 */

    /************************  2 FACTORS   ****************************/

    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);
            ObsIndexL = ObsIndex + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                for (s = 0; s < NbStates; s++)
                {
                    /* interpolate in the Rib dimension */
                    for (r = 0; r < NbRibActive; r++)
                    {
                        SliceL = RibLadder[r][s] + offset;
                        tmpInValues[r] = SliceL[j];
                    }
                    for (r = 0; r < NbRibActive-1; r++)
                    {
                        if (Fix3_DoubleQuadraticInterp (
                                RibInRangeCount,
                                tmpInValues,
                                NbRibActive,
                                r + ObsIndexL[j],
                                &(tmpOutValues[r])) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }
                
                    /* store results */
                    for (r = 0; r < NbRibActive-1; r++)
                    {
                        SliceL = RibLadder[r][s] + offset;
                        SliceL[j] = tmpOutValues[r];
                    }

                } /* for state idx */

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
                ObsIndexL = ObsIndex + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* interpolate in the Rib dimension */
                        for (r = 0; r < NbRibActive; r++)
                        {
                            SliceL = RibLadder[r][s] + offset;
                            tmpInValues[r] = SliceL[k];
                        }
                        for (r = 0; r < NbRibActive-1; r++)
                        {
                            if (Fix3_DoubleQuadraticInterp (
                                    RibInRangeCount,
                                    tmpInValues,
                                    NbRibActive,
                                    r + ObsIndexL[k],
                                    &(tmpOutValues[r])) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
            
                        /* store results */
                        for (r = 0; r < NbRibActive-1; r++)
                        {
                            SliceL = RibLadder[r][s] + offset;
                            SliceL[k] = tmpOutValues[r];
                        }

                    } /* for state idx */

                }  /* for k */
            }  /* for j */
    }  /* if then else */

    status = SUCCESS;
    
RETURN:

    return (status);

}



/*****  Calc_Ladder  ***********************************************/
/**
*       Main calculation routine: discounted expected value of cash-flows
*       going backward in the tree.
*/
int     Calc_Ladder
             (MKTVOL_DATA         *mktvol_data,      /**< (I) Vol data        */
              LADDER_DATA         *ladder_data,      /**< (I) Deal data       */
              FIX3_TREE_DATA      *tree_data,        /**< (I) Tree data       */
              OPT_OUT_DATA        *opt_out_data)     /**< (I) Output data     */
{

    FIX3_DEV_DATA    dev_data;        /* Fix3_Dev data structure          */
                                                                          
    /* Variables */
    double    ***Ladder = NULL;       /* Callable ladder swap             */
    double      *KnownCpnSt = NULL;   /* Known cpn for ladder             */
    double      *KnownCpnFl = NULL;   /* Known cpn for floater            */
    double      *IndexSt[2] = {NULL, NULL};/* Underlying ladder indices   */
    double      *IndexLadder = NULL;  /* Combination of ladder indices    */
    double      *IndexPaySt = NULL;   /* Payment ladder index             */
    double      *IndexFl = NULL;      /* Underlying floating index        */
    double      *Identity = NULL;     /* Slice set to 1                   */
    double      *Annuity = NULL;      /* Slice for dummy annuity in yield */
    double      *Funding = NULL;      /* Slice for all funding            */
    double      *Binary  = NULL;      /* Slice for ladder binary step     */
    double      *Step    = NULL;      /* Slice for ladder step            */
    double      *Spread  = NULL;      /* Slice for spread                 */
    double      *AuxSlice = NULL;     /* Slice for call proba             */
    double      *Indicator = NULL;    /* Slice for call proba             */

    double      *ZeroToPmtSt = NULL;
    double      *ZeroToPmtFl = NULL;

    /* Flags */
    int         ExerFlag;
    int         KnownCpnFlagSt;
    int         KnownCpnFlagFl;
    int         ResetFlagSt;
    int         ResetFlagFl;
    int         StateVarFlag;
    int         ZbkResetFlag[3];
    int         IsPastFirstResetSt;
    int         IsPastFirstResetFl;
    int         StatFlag;

    int         IdxReset;

    /* Banks and bank variables */
    CLAIM_BANK  ZeroBank[3];          /* zero bank for 3 curves             */
    int         CvToUseSt[2];   
    int         CvToUseFl;
    int         DCurve;    

    /* State variables */
    double     *PrevStates = NULL;    /* previous state levels at i   */
    double     *CurrStates = NULL;    /* current state levels  at i-1 */
    
    double     *RibInitStates = NULL; /* for RIB reset-in-advance     */
    /* Numeric variables */
    double      Strike=0.;            /* Current strike                     */

    double      KCAmtSt=0.;           /* Known coupon amount for ladder     */
    double      KCAmtFl=0.;           /* Known coupon amount for float      */
    double      KCFactorSt=0.;        /* Known coupon factor for ladder     */

    double      CpnFloorSt=0.;        /* St coupon floor                    */
    double      CpnCapSt=0.;          /* St coupon cap                      */
    double      CpnBarrLoSt=0.;       /* St coupon barrier low              */
    double      CpnBarrHiSt=0.;       /* St coupon barrier high             */
    double      CpnOutsSt=0.;         /* St coupon outstanding              */
    double      CpnDcfSt=0.;          /* St coupon day count fraction       */
    double      CpnDownSt=0.;         /* St coupon down rate                */
    double      CpnMidSt=0.;          /* St coupon mid rate                 */
    double      CpnUpSt=0.;           /* St coupon up rate                  */
    double      CpnIdxWtSt[2]={0., 0.};/* St coupon index weight            */
    double      CpnSprdSt=0.;         /* St coupon spread                   */
    double      CpnStickyCoef=0.;     /* St coupon ladder coefficient       */
    double      CpnFloorFl=0.;        /* Fl coupon floor                    */
    double      CpnCapFl=0.;          /* Fl coupon cap                      */
    double      CpnUpFl=0.;           /* Fl coupon up rate                  */
    double      CpnOutsFl=0.;         /* Fl coupon outstanding              */
    double      CpnDcfFl=0.;          /* Fl coupon day count fraction       */
    double      CpnLeverage=0.;       /* St coupon leverage                 */
    double      CpnIdxFloor=0.;       /* St coupon index floor              */
    double      CpnIdxCap=0.;         /* St coupon index cap                */

    long        MatDate;              /* Final maturity date                */
    long        CurrentDate;          /* Current date                       */
    long        ValueDate;            /* Value date                         */
    long        CpnPmtDateSt=0L;
    long        CpnPmtDateFl=0L;
    long        ZbkErDate[3];         /* EarliestUse date for curr zero mat */

    int         isKnownUsedSt=0;      /* 1 if zero for known pmt is used    */
    int         isKnownUsedFl=0;      /* 1 if zero for known pmt is used    */
    int         isArrearsSt=0;
    int         isArrearsFl=0;
    int         NbStates=ladder_data->NbStates;

    int         isRibAdvanceSt=0;     /* 1 if RIB reset-in-advance          */

    char        SoZ    = ladder_data->SoZ;
    char        AoM    = ladder_data->AoM;
    char        CompSt = ladder_data->CompSt;

    /* Rib slices and variables */
    double     *ObsIndex    = NULL;   /* Rib in/out range indicator         */
    double     *IndexRib[2] ={NULL,NULL}; /* Rib obs indices                */

    double      RibFrac= 0.;         /* percentage observations in range   */
    double      RibInRangeWeight = 0.;/* inside  weight this complex period */
    double      RibOutRangeWeight= 0.;/* outside weight this complex period */

    int         RibResetFlag; 
    int         CplxCpnFlag;
    int         CvToUseRib[2];     
    int         RibIdxTotal = 0;      /* Index of this obs in input list    */
    int         RibIdxInPer;          /* Index of this obs in current period*/
    int         NbRibObsInPer = ladder_data->MaxNbRib;    
                                      /* Nb Rib obs this ladder period      */
    int         NbRibActive = NbRibObsInPer + 1;        
                                      /* Nb active Rib slices               */

    int         NbPastRibObsDates = ladder_data->NbPastRibObsDates;
    int         r;

    int         T;                    /* Last time point                    */
    int         t;                    /* Current time point                 */
    int         offset;               /* Node offset at t=0                 */

    int         status = FAILURE;     /* Error status = FAILURE initially   */
    int         i, j, s;

    FIX3_EVENT_STATS_SCHEDULE**  evs = NULL; /* vector of calculators        */
    double***                    xsc = NULL; /* vector of slice pointers     */



    /* allocate event stats calculators */
    if (ladder_data->CalcStats == 'Y')
    {
        evs = (FIX3_EVENT_STATS_SCHEDULE**) 
                malloc (sizeof(FIX3_EVENT_STATS_SCHEDULE*) * NbRibActive);
 
        xsc = (double***) malloc (sizeof(double**) * NbRibActive);               

        if (!evs || !xsc)
        {
            DR_Error ("failed to allocate stats calculator memory");
            goto FREE_MEM_AND_RETURN;
        }

        for (r = 0; r < NbRibActive; r++)
        {
            evs[r] = (FIX3_EVENT_STATS_SCHEDULE*) 
                       malloc (sizeof(FIX3_EVENT_STATS_SCHEDULE) * NbStates);
            xsc[r] = (double **) malloc(sizeof(double*) * NbStates);

            if (!evs[r] || !xsc[r])
            {
                DR_Error ("failed to allocate stats calculator memory");
                goto FREE_MEM_AND_RETURN;
            }
        }

        for (r = 0; r < NbRibActive; r++)
        {
            for (s=0; s<NbStates; ++s)
            {
                if (Fix3EventStatsScheduleInit(
                        &(evs[r][s]), 
                        tree_data, 
                        &dev_data) != SUCCESS)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }
        }
    }

    Fix3_Dev_Init(&dev_data);

    /* Must init claimbanks first to NULL their pointer members */
    for (j=0; j<3; j++) 
    {
        Fix3_CbkInit(&(ZeroBank[j]));
        ZbkErDate[j] = 99999999L; /* init memory: avoid purify's complaint */
    }

    for(i =0 ; i < 2; i++)
        IndexSt[i] = Fix3_Alloc_Slice(tree_data);

    /*  Allocation of variables for tree pricing. */
    KnownCpnSt  = Fix3_Alloc_Slice(tree_data);
    KnownCpnFl  = Fix3_Alloc_Slice(tree_data);
    IndexFl     = Fix3_Alloc_Slice(tree_data);
    IndexLadder = Fix3_Alloc_Slice(tree_data);
    IndexPaySt  = Fix3_Alloc_Slice(tree_data);
    Identity    = Fix3_Alloc_Slice(tree_data);
    Annuity     = Fix3_Alloc_Slice(tree_data);
    Funding     = Fix3_Alloc_Slice(tree_data);
    Binary      = Fix3_Alloc_Slice(tree_data);
    Step        = Fix3_Alloc_Slice(tree_data);
    Spread      = Fix3_Alloc_Slice(tree_data);
    AuxSlice    = Fix3_Alloc_Slice(tree_data);
    Indicator   = Fix3_Alloc_Slice(tree_data);
    ObsIndex    = Fix3_Alloc_Slice (tree_data);
    IndexRib[0] = Fix3_Alloc_Slice (tree_data);
    IndexRib[1] = Fix3_Alloc_Slice (tree_data);

    if (   (KnownCpnSt  == NULL)
        || (KnownCpnFl  == NULL)
        || (IndexSt[0]  == NULL)
        || (IndexSt[1]  == NULL)
        || (IndexLadder == NULL)
        || (IndexPaySt  == NULL)
        || (IndexFl     == NULL)
        || (Identity    == NULL)
        || (Annuity     == NULL)
        || (Funding     == NULL)
        || (Binary      == NULL)
        || (Step        == NULL)
        || (AuxSlice    == NULL)
        || (Indicator   == NULL)
        || (Spread      == NULL)
        || (ObsIndex    == NULL)
        || (IndexRib[0] == NULL)
        || (IndexRib[1] == NULL))
    {
        DR_Error ("Calc_ladder: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;
    }

    Ladder = (double ***) DR_Array (DOUBLE_D_PTR, 0, ladder_data->MaxNbRib);
    if (Ladder == NULL)
    {
        DR_Error ("Calc_ladder: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }

    for (r = 0; r <= ladder_data->MaxNbRib; r++)
    {
        Ladder[r] = (double **) DR_Array(DOUBLE_PTR, 0, ladder_data->NbStates-1);
        if (Ladder[r] == NULL)
        {
            DR_Error ("Calc_ladder: could not allocate memory!");
            goto FREE_MEM_AND_RETURN;        
        }

        for (i = 0; i < ladder_data->NbStates; i++)
        {
            Ladder[r][i] = Fix3_Alloc_Slice (tree_data);
            if (Ladder[r][i] == NULL)
            {
                DR_Error ("Calc_ladder: could not allocate memory!");
                goto FREE_MEM_AND_RETURN;
            }
        }  
    }
        
    if (Fix3_Dev_Alloc (&dev_data, tree_data) == FAILURE)
    {
        DR_Error ("Calc_ladder: could not allocate dev memory!");
        goto FREE_MEM_AND_RETURN;
    }

    /* Allocate zero bank */
    for (j=0; j<3; j++)
    {
        if (tree_data->NbZeros[ZbkEVENT+j] > 0)
        {
            if (Fix3_CbkAlloc(&(ZeroBank[j]),
                         tree_data->NbZeros[ZbkEVENT+j],
                         tree_data) == FAILURE)
            {
                DR_Error ("Calc_ladder: could not allocate ZeroBank!");
                goto FREE_MEM_AND_RETURN;
            }
        }
    }

    T = tree_data->NbTP;        
   
    /* Set a 1.0 slice to be a zero to pmt for in-arrears case */
    if (Fix3_Init_Slice(Identity,
                   1.0,
                   tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    ZeroToPmtSt = Identity;
    ZeroToPmtFl = Identity;

    for(i = 0; i < 2; i++)
        CvToUseSt[i]  = ladder_data->IdxIoDSt[i];
    for (j = 0; j < 2; j++) 
        CvToUseRib[j] = ladder_data->RibIdxIoD[j];

    CvToUseFl   = ladder_data->IdxIoDFl;
    DCurve      = tree_data->CvDisc;
    MatDate     = ladder_data->MatDate;

    CurrentDate = ValueDate = tree_data->TPDate[0];
  
    isArrearsSt = (ladder_data->ArrearsSt == 'Y');
    isArrearsFl = (ladder_data->ArrearsFl == 'Y');

    isRibAdvanceSt = (ladder_data->CplxIsRib && !isArrearsSt);

    /* Main loop going backwards */

    for (t = T; t >= 0; t--)
    {

        CurrentDate = tree_data->TPDate[t];

        /* Test if has passed the first reset date */
        IsPastFirstResetSt = (CurrentDate < ladder_data->FirstResetDateSt);
        IsPastFirstResetFl = (CurrentDate < ladder_data->FirstResetDateFl);

        /* Update event flags... */  
        ExerFlag       = tree_data->TPtype[0][t];
        KnownCpnFlagSt = tree_data->TPtype[1][t];
        KnownCpnFlagFl = tree_data->TPtype[3][t];
        ResetFlagSt    = tree_data->TPtype[2][t];
        ResetFlagFl    = tree_data->TPtype[4][t];
        StateVarFlag   = tree_data->TPtype[5][t];
        RibResetFlag   = tree_data->TPtype[9][t];
        CplxCpnFlag    = tree_data->TPtype[10][t];
        StatFlag       = tree_data->TPtype[11][t];

        if (StateVarFlag != ResetFlagSt) goto FREE_MEM_AND_RETURN;   

        /* Zerobank flags */
        for (j=0; j<3; j++)
        {
            ZbkResetFlag[j] = tree_data->TPtype[ZbkEVENT+j][t];
            if (ZbkResetFlag[j])
            {
                ZbkErDate[j]=(tree_data->CritDate[ZbkEVENT+j][t]).SuppDate[0];
            } 
        }

        /* .. and event amounts stored alongside the timeline */
        if (ExerFlag)
        {
            Strike = (tree_data->CritDate[0][t]).Value[0];
        }

        if (KnownCpnFlagSt && !isRibAdvanceSt) 
        {   
            KCAmtSt    = (tree_data->CritDate[1][t]).Value[0];
            KCFactorSt = (tree_data->CritDate[1][t]).Value[4];
        }

        if (KnownCpnFlagFl) 
        {   
            KCAmtFl  = (tree_data->CritDate[3][t]).Value[0];
        }

        if (ResetFlagSt)
        {   
            CpnFloorSt    = (tree_data->CritDate[2][t]).Value[0];
            CpnCapSt      = (tree_data->CritDate[2][t]).Value[1];
            CpnOutsSt     = (tree_data->CritDate[2][t]).Value[2];
            CpnDcfSt      = (tree_data->CritDate[2][t]).Value[3];
            CpnDownSt     = (tree_data->CritDate[6][t]).Value[0];
            CpnMidSt      = (tree_data->CritDate[6][t]).Value[1];
            CpnUpSt       = (tree_data->CritDate[6][t]).Value[2];
            CpnBarrLoSt   = (tree_data->CritDate[6][t]).Value[3];
            CpnBarrHiSt   = (tree_data->CritDate[6][t]).Value[4];
            CpnIdxWtSt[0] = (tree_data->CritDate[7][t]).Value[0];
            CpnIdxWtSt[1] = (tree_data->CritDate[7][t]).Value[1];
            CpnSprdSt     = (tree_data->CritDate[7][t]).Value[2];
            CpnStickyCoef = (tree_data->CritDate[7][t]).Value[3];
            CpnLeverage   = (tree_data->CritDate[8][t]).Value[0];
            CpnIdxFloor   = (tree_data->CritDate[8][t]).Value[1];
            CpnIdxCap     = (tree_data->CritDate[8][t]).Value[2];
            CpnPmtDateSt  = (tree_data->CritDate[2][t]).SuppDate[0];
        }

        if (ResetFlagFl)
        {   
            CpnFloorFl     = (tree_data->CritDate[4][t]).Value[0];
            CpnCapFl       = (tree_data->CritDate[4][t]).Value[1];
            CpnUpFl        = (tree_data->CritDate[4][t]).Value[2];
            CpnOutsFl      = (tree_data->CritDate[4][t]).Value[3];
            CpnDcfFl       = (tree_data->CritDate[4][t]).Value[4];
            CpnPmtDateFl   = (tree_data->CritDate[4][t]).SuppDate[0];
        }

        if (StateVarFlag)
        {   
           IdxReset        = (int) (tree_data->CritDate[5][t]).Value[0];
           PrevStates      = ladder_data->State[IdxReset];
           CurrStates      = ladder_data->State[IdxReset-1];
        }

        if (RibResetFlag)
        {
            RibIdxTotal    = (int) (tree_data->CritDate[9][t]).Value[0];
            RibIdxInPer    = (int) (tree_data->CritDate[9][t]).Value[1] ; 
        }

        /*  'Update' tree. */
        if (Fix3_Lattice (&dev_data,
                     t,                                      
                     T,
                     mktvol_data,                                   
                     tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;   
        }

        /* Update the zero banks */
        for (j=0; j<3; j++)
        {
            if (Fix3_ZbkUpdate
                  (&(ZeroBank[j]),
                   ZbkResetFlag[j],
                   CurrentDate,
                   ZbkErDate[j], /* EarliestUse date */
                   t,
                   T,
                   j, /* curve id */
                   &dev_data,
                   tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }



        /* Discount all the state variables */
        for (r = 0; r < NbRibActive; r++)
        {
            for (s = 0; s < NbStates; s++)
            {
                if (Fix3_Dev(Ladder[r][s],
                    t,
                    T,
                    DCurve,
                    &dev_data,
                    tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }
        }
        /* Reset Rib index (if present) and update Tarn index; this has to be
         * done BEFORE resetting the complex index */
        if (RibResetFlag)
        {
            for (j = 0; j < 2; j++)
            {
                if (ladder_data->RibIdxOn[j])
                {

                    if (Fix3_ZbkParYield_t (IndexRib[j],
                                            AuxSlice,
                                            &(ZeroBank[CvToUseRib[j]]),
                                            CurrentDate,
                                            CurrentDate,
                                            ladder_data->RibIdxMat[j],
                                            ladder_data->RibIdxDCC[j],
                                            ladder_data->RibIdxFreq[j],
                                            0.,    /* spread */
                                            t,
                                            tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                } /* if RibIdxOn */
            } /* for j */

            /* Combine the two Rib indices */
            if (Fix3_LCombTwoSlices (AuxSlice,
                                     IndexRib[0],
                                     ladder_data->RibIdxWeight[0],
                                     IndexRib[1],
                                     ladder_data->RibIdxWeight[1],
                                     t, 
                                     tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }

            /* Calculate ObsIndex profile: 1 inside, 0 outside */
            if (Fix3_Set_Slice (ObsIndex,
                                1.,
                                t,
                                tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }

            if (Fix3_KoOption_t(ObsIndex,
                                AuxSlice,
                                1,     /* knockout now */
                                ladder_data->RibLoBarrier[RibIdxTotal],
                                ladder_data->RibHiBarrier[RibIdxTotal],
                                0.,
                                'O',                 
                                ladder_data->RibSmoothing,
                                t,
                                tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }

            /* If Statistics required, interpolate calculators */
            if (ladder_data->CalcStats == 'Y')
            {
                for (r = 0; r < NbRibActive; r++)
                    for (s = 0; s < NbStates; s++)
                        xsc[r][s] = evs[r][s].probSlice;
                if (Fix3_RibLadderInterpolate_t (xsc, 
                                                 NbStates, 
                                                 NbRibActive,
                                                 ObsIndex,
                                                 t,
                                                 tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                for (r = 0; r < NbRibActive; r++)
                    for (s = 0; s < NbStates; s++)
                        xsc[r][s] = evs[r][s].timeSlice;
                if (Fix3_RibLadderInterpolate_t (xsc, 
                                                 NbStates, 
                                                 NbRibActive,
                                                 ObsIndex,
                                                 t,
                                                 tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
                   
                for (r = 0; r < NbRibActive; r++)
                    for (s = 0; s < NbStates; s++)
                        xsc[r][s] = evs[r][s].tsqrSlice;
                if (Fix3_RibLadderInterpolate_t (xsc, 
                                                 NbStates, 
                                                 NbRibActive,
                                                 ObsIndex,
                                                 t,
                                                 tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                for (i = 0; i < (int) evs[0][0].size; i++)
                {
                    for (r = 0; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].sched[i].probSlice;
                    if (Fix3_RibLadderInterpolate_t (xsc, 
                                                     NbStates, 
                                                     NbRibActive,
                                                     ObsIndex,
                                                     t,
                                                     tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }
            } /* if Calc Stats */

            /* Revise Rib state variables */
            if (Fix3_RibLadderInterpolate_t (Ladder,
                                             NbStates,
                                             NbRibActive,
                                             ObsIndex,
                                             t,
                                             tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }

            /* Drop one active slice -- fewer Rib observations
             * remaining in the period */
            NbRibActive--;

        } /* if RibResetFlag */

        /* Reset zero for known ladder cpn, or dev it if already used */
        if (KnownCpnFlagSt && !isRibAdvanceSt)
        {
            if (Fix3_Set_Slice(KnownCpnSt,
                          KCAmtSt,
                          t,
                          tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            isKnownUsedSt = 1;
        }
        else
        if (isKnownUsedSt)
        {
            if (Fix3_Dev (KnownCpnSt,
                     t,
                     T,
                     DCurve,  
                     &dev_data,
                     tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
        }

        /* Reset zero for known ladder cpn, or dev it if already used */
        if (KnownCpnFlagFl)
        {
            if (Fix3_Set_Slice(KnownCpnFl,
                          KCAmtFl,
                          t,
                          tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            isKnownUsedFl = 1;
        }
        else
        if (isKnownUsedFl)
        {
            if (Fix3_Dev (KnownCpnFl,
                     t,
                     T,
                     DCurve,  
                     &dev_data,
                     tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
        }

        /* Reset ladder index */
        for(i = 0 ; i < 2; i++)
        {
            if (ResetFlagSt)
            {
                if (ladder_data->IdxOnSt[i])
                {
                    if (Fix3_ZbkParYield_t(IndexSt[i],
                                  Annuity,
                                  &(ZeroBank[CvToUseSt[i]]),
                                  CurrentDate,
                                  CurrentDate,
                                  ladder_data->IdxMatSt[i],
                                  ladder_data->IdxBaseSt[i],
                                  ladder_data->IdxFreqSt[i],
                                  0.0,    /* no spread */
                                  t,
                                  tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    } 
                }
            }
        }
        if (Fix3_LCombTwoSlices(IndexLadder,
                               IndexSt[0],
                               ladder_data->IdxObsWeightSt[0],
                               IndexSt[1],
                               ladder_data->IdxObsWeightSt[1],
                               t,
                               tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
        
        /* Reset floating index, weight it, add spread, cap&floor */
        if (ResetFlagFl)
        {
            if (ladder_data->IdxOnFl)
            {  
                if (Fix3_ZbkParYield_t(IndexFl,
                                  Annuity,
                                  &(ZeroBank[CvToUseFl]),
                                  CurrentDate,
                                  CurrentDate,
                                  ladder_data->IdxMatFl,
                                  ladder_data->IdxBaseFl,
                                  ladder_data->IdxFreqFl,
                                  0.0,    /* no spread */
                                  t,
                                  tree_data) == FAILURE)
                {
                goto FREE_MEM_AND_RETURN;
                }  
            }
            if (Fix3_MultiplyScalar(IndexFl,  
                               ladder_data->IdxWeightFl,
                               t,
                               tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
            if (Fix3_AddScalar(IndexFl,
                          CpnUpFl,
                          t,
                          tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            if (Fix3_MaxMinOnSlice(IndexFl,
                              CpnFloorFl,
                              CpnCapFl,
                              t,
                              tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }

        /* Get flt pmt zero from bank if on reset date AND in advance */
        /* AND not RIBed. Otherwise zeros are simply equal to 1.      */
        if (ResetFlagSt && (ladder_data->ArrearsSt == 'N') &&
            !isRibAdvanceSt)
        {
            ZeroToPmtSt = Fix3_ZbkReadZero(&(ZeroBank[DCurve]),
                                          CpnPmtDateSt,
                                          FALSE,
                                          CurrentDate,
                                          t,
                                          tree_data);
        }

        if (ResetFlagFl && (ladder_data->ArrearsFl == 'N'))
        {
            ZeroToPmtFl = Fix3_ZbkReadZero(&(ZeroBank[DCurve]),
                                          CpnPmtDateFl,
                                          FALSE,
                                          CurrentDate,
                                          t,
                                          tree_data);
        }

        if (CurrentDate <= MatDate)
        {
            /* Update funding for advance reset */
            if (Fix3_Floater_t(Funding,
                          IndexFl,
                          ZeroToPmtFl,
                          ResetFlagFl && !isArrearsFl,
                          CpnDcfFl,
                          0.0,
                          ladder_data->CompFl,
                          ladder_data->ArrearsFl,
                          CpnOutsFl,
                          0.,       
                          0,   
                          t,
                          T,
                          DCurve,
                          &dev_data,
                          tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            
            
            
            /* Prepare step-up rate if on reset date */
            if (ResetFlagSt)
            {
                /* Smooth binarstep */
                if (Fix3_LadderStep_t (Binary,
                                  IndexLadder,
                                  CpnBarrLoSt,
                                  CpnBarrHiSt,
                                  CpnDownSt,
                                  CpnMidSt,
                                  CpnUpSt,
                                  ladder_data->Smoothing,
                                  t,
                                  tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
                
                

                /* Linear combination of indexSt and spread */
                if (Fix3_Set_Slice (Step,
                               CpnSprdSt,
                               t,
                               tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                if (Fix3_LCombTwoSlices (IndexPaySt,
                                         IndexSt[0],
                                         CpnIdxWtSt[0],
                                         IndexSt[1],
                                         CpnIdxWtSt[1],
                                         t, 
                                         tree_data) == FAILURE)
                {
                     goto FREE_MEM_AND_RETURN;
                }

                if (Fix3_AddTwoSlices (Spread,
                                    Step,
                                    IndexPaySt,
                                    t,
                                    tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
              
                /* Collar Spread & multiply by Leverage */
                if (Fix3_MaxMinOnSlice(Spread,
                                       CpnIdxFloor,                                  
                                       CpnIdxCap,
                                       t,
                                       tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                if (Fix3_MultiplyScalar(Spread,
                                        CpnLeverage,
                                        t,
                                        tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

            } /* if ResetFlagSt */

            /* Add ladder and float coupon if reset in advance     */
            /* In case of a RIB, only revise Last Coupon Rate states
               if it is the first reset date                       */
            for (r = 0; r < NbRibActive; r++)
            {
                /* RIB reset-in-advance: if this is the first reset date,
                   revise Coupon Rate states */
                if (isRibAdvanceSt && ResetFlagSt)
                {
                    if (NbRibActive != 1) {
                        DR_Error("RIB reset-in-advance: NbRibActive != 1 "
                                 "at reset");
                        goto FREE_MEM_AND_RETURN;
                    }

                    if (Fix3_Sticky_LadderSwap_t(Ladder[r],
                                 /* RIB Ladder Reset-in-Advance reset
                                    if this is the first reset date */ 
                                 CplxCpnFlag? 0:3,
                                 (ResetFlagSt||IsPastFirstResetSt) && (ResetFlagFl||IsPastFirstResetFl),
                                 Binary,
                                 Spread,
                                 AoM,
                                 Funding,
                                 NULL, /* not used */
                                 CpnFloorSt,
                                 CpnCapSt,
                                 CpnStickyCoef,
                                 1.0, /* RibFrac = 1.0 */
                                 CpnOutsSt,
                                 CpnDcfSt,
                                 SoZ,   
                                 CompSt,
                                 NbStates,
                                 CurrStates,
                                 PrevStates,
                                 t,
                                 tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }
                else {
 
                    if (Fix3_Sticky_LadderSwap_t(Ladder[r],
                                 ResetFlagSt && !isArrearsSt,
                                 (ResetFlagSt||IsPastFirstResetSt) && (ResetFlagFl||IsPastFirstResetFl),
                                 Binary,
                                 Spread,
                                 AoM,
                                 Funding,
                                 ZeroToPmtSt,
                                 CpnFloorSt,
                                 CpnCapSt,
                                 CpnStickyCoef,
                                 1.0, /* RibFrac = 1.0 */
                                 CpnOutsSt,
                                 CpnDcfSt,
                                 SoZ,   
                                 CompSt,
                                 NbStates,
                                 CurrStates,
                                 PrevStates,
                                 t,
                                 tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }   
                }
            }

            /* Clean up Funding -- use the same flag as for FundResetFlag above */
            if ((ResetFlagSt||IsPastFirstResetSt) && (ResetFlagFl||IsPastFirstResetFl))
            {
                if (Fix3_Set_Slice (Funding,
                                    0.,
                                    t,
                                    tree_data) == FAILURE)
                {
                        goto FREE_MEM_AND_RETURN;
                } 
            }

            /* ready to exercise - see below - update call stats */
            if (ladder_data->CalcStats == 'Y') 
            {
                /* update stats calculators */
                for (r = 0; r <= ladder_data->MaxNbRib; r++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        /* calculate indicator slice aux = -(ladder + strike) */
                        if (ExerFlag && 
                                (Fix3_AddScalar2 (AuxSlice, Ladder[r][s], Strike, t, tree_data) == FAILURE ||
                                 Fix3_MultiplyScalar (AuxSlice, -1.0, t, tree_data) == FAILURE ||
                                 Fix3_ExIndicator(Indicator, AuxSlice, 0.0, 'A', ladder_data->Smoothing, t, tree_data) == FAILURE))
                            goto FREE_MEM_AND_RETURN;

                        Fix3EventStatsScheduleUpdate(&evs[r][s], ExerFlag ? Indicator : NULL, StatFlag, t);
                    }
                }

                /* interpolate */
                if (ResetFlagSt && !isArrearsSt)
                {
                    for (r = 0; r < NbRibActive; r++)
                    {
                        if (NbRibActive != 1) {
                            DR_Error("reset-in-advance: NbRibActive != 1 "
                                     "at reset");
                            goto FREE_MEM_AND_RETURN;
                        }

                        RibFrac = 1.;

                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].probSlice;
                        if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                    CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                    NbStates, CurrStates, PrevStates,
                                                    t, tree_data) == FAILURE)
                            goto FREE_MEM_AND_RETURN;

                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].timeSlice;
                        if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                    CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                    NbStates, CurrStates, PrevStates,
                                                    t, tree_data) == FAILURE)
                            goto FREE_MEM_AND_RETURN;

                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].tsqrSlice;
                        if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                    CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                    NbStates, CurrStates, PrevStates,
                                                    t, tree_data) == FAILURE)
                            goto FREE_MEM_AND_RETURN;

                        for (i = 0; i < (int) evs[0][0].size; i++)
                        {
                            for (s = 0; s < NbStates; s++)
                                xsc[r][s] = evs[r][s].sched[i].probSlice;
                            if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                        CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                        NbStates, CurrStates, PrevStates,
                                                        t, tree_data) == FAILURE)
                                goto FREE_MEM_AND_RETURN;
                        }
                    }
                }
            }

            /* Exercise the callable option */
            /* IT WORKS ON FINAL MATURITY AS WELL: STRIKE IS 0 THERE */
            if (ExerFlag)
            {
                for (r = 0; r < NbRibActive; r++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        if (Fix3_MinOnSlice(Ladder[r][s],
                                       -Strike,
                                       t,
                                       tree_data) == FAILURE)
                        {
                            goto FREE_MEM_AND_RETURN;
                        }
                    }
                }
            }

            /* Reset NbRibObsInPer -- waited until now since there may be a Rib
             * Obs on accrual start; from now on, we are processing the ladder
             * reset-in-arrears case as well as reset-in-advance. 
             * Also recall that there are no Rib observations on accrual end
             * dates (=payment dates), hence these settings will work
             * on all upcoming Rib reset observations in the complex period */
            if (CplxCpnFlag)
            {
                NbRibObsInPer  = (int) (tree_data->CritDate[10][t]).Value[0];

                RibInRangeWeight  = tree_data->CritDate[10][t].Value[1];
                RibOutRangeWeight = tree_data->CritDate[10][t].Value[2];

                /* Set number of active Rib slices */
                NbRibActive    = NbRibObsInPer + 1;

                offset = Fix3_Node_Offset(tree_data->NbFactor, 0, 0, t, tree_data);

                /* Copy value of r = 0 slice to all the active ones */
                for (r = 1; r < NbRibActive; r++)
                {
                    for (s = 0; s < NbStates; s++)
                    {
                        if (Fix3_Copy_Slice (Ladder[r][s],
                                             Ladder[0][s],
                                             t,
                                             tree_data) == FAILURE)
                        {
                            goto FREE_MEM_AND_RETURN;
                        }
                    }
                }


                /* If Statistics required, interpolate calculators */
                if (ladder_data->CalcStats == 'Y')
                {
                    for (r = 0; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].probSlice;
                    for (r = 1; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++) 
                            if (Fix3_Copy_Slice (xsc[r][s], 
                                                 xsc[0][s],
                                                 t,
                                                 tree_data) == FAILURE)
                                goto FREE_MEM_AND_RETURN;

                    for (r = 0; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].timeSlice;
                    for (r = 1; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++) 
                            if (Fix3_Copy_Slice (xsc[r][s], 
                                                 xsc[0][s],
                                                 t,
                                                 tree_data) == FAILURE)
                                goto FREE_MEM_AND_RETURN;
                   
                    for (r = 0; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].tsqrSlice;
                    for (r = 1; r < NbRibActive; r++)
                        for (s = 0; s < NbStates; s++) 
                            if (Fix3_Copy_Slice (xsc[r][s], 
                                                 xsc[0][s],
                                                 t,
                                                 tree_data) == FAILURE)
                                goto FREE_MEM_AND_RETURN;

                    for (i = 0; i < (int) evs[0][0].size; i++)
                    {
                        for (r = 0; r < NbRibActive; r++)
                            for (s = 0; s < NbStates; s++)
                                xsc[r][s] = evs[r][s].sched[i].probSlice;
                        for (r = 1; r < NbRibActive; r++)
                            for (s = 0; s < NbStates; s++) 
                                if (Fix3_Copy_Slice (xsc[r][s], 
                                                     xsc[0][s],
                                                     t,
                                                     tree_data) == FAILURE)
                                    goto FREE_MEM_AND_RETURN;
                    }
                } /* if Calc Stats */
            }

            /* For RIB Ladder Reset-in-Advance, we update the ladder leg */
            /* on coupon dates */
            if (CplxCpnFlag && isRibAdvanceSt)
            {

            double  newCpnOutsSt;
            double  newCpnDcfSt;

            /* Check if the reset date of this coupon payment is in the
             * past, which means the unRIBed coupon rate is known      */ 
            if (KnownCpnFlagSt)
            {
                newCpnOutsSt = (tree_data->CritDate[1][t]).Value[3];
                newCpnDcfSt  = (tree_data->CritDate[1][t]).Value[2];
                /* If CurrStates is not set, then this must be the 
                 * last coupon payment date                        */
                if (!CurrStates)
                    CurrStates= ladder_data->State[ladder_data->FirstResetI-1]; 

            }
            else
            {
                /* This is the case on the last payment day */
                /* Search for the corresponding reset date and retrieve
                 *  CpnOutsSt CpnDcfSt and prev/currStates             */
                int     resetIdx;
                for (resetIdx = t; resetIdx >= 0; resetIdx--)
                {
                    if ((tree_data->CritDate[2][resetIdx]).SuppDate[0] ==
                        CurrentDate)   break;
                }

                if (resetIdx < 0)
                {
                    DR_Error("calc_ladder: could not find the reset date from"
                             " payment date\n");
                    goto FREE_MEM_AND_RETURN;
                }

                newCpnOutsSt = (tree_data->CritDate[2][resetIdx]).Value[2];
                newCpnDcfSt  = (tree_data->CritDate[2][resetIdx]).Value[3];

                /* If CurrStates is not set, then this must be the 
                 * last coupon payment date                        */
                if (!CurrStates) 
                {
                    int newIdxReset;

                    if (!tree_data->TPtype[5][resetIdx]) {
                        DR_Error("calc_ladder: StateVarFlag not set on the "
                                "last reset date of a RIB reset-in-advabce\n");
                        goto FREE_MEM_AND_RETURN;
                    }

                    newIdxReset =
                        (int) (tree_data->CritDate[5][resetIdx]).Value[0];
                    CurrStates      = ladder_data->State[newIdxReset];
                } 

            } /* end of if (KnownCpnFlagSt) */

            if (NbRibObsInPer == 0)
            {
                DR_Error("Calc_Ladder: "
                         "No rib observation in a RIB accural period");
                goto FREE_MEM_AND_RETURN;
            }

            for (r = 0; r < NbRibActive; r++)
            {
                /* Set percentage of Rib observations in range if applicable */
                RibFrac = RibInRangeWeight  * r / NbRibObsInPer +
                          RibOutRangeWeight * (1- (double) r / NbRibObsInPer);

                if (!ResetFlagSt)
                {
                    /* This is the last payment date */
                    if (Fix3_Sticky_LadderSwap_t(Ladder[r],
                                 /*--- RIB Ladder Reset-in-Adv last payment */
                                 4,
                                 0, /* No funding added here */
                                 NULL, /* not used */
                                 NULL, /* not used */
                                 0,    /* not used */
                                 NULL, /* not used */
                                 NULL, /* not used */
                                 0.0,  /* not used */
                                 0.0,  /* not used */
                                 0.0,  /* not used */
                                 RibFrac,
                                 newCpnOutsSt,
                                 newCpnDcfSt,
                                 SoZ,
                                 CompSt,
                                 NbStates,
                                 CurrStates, /* unRIBed coupon rate */
                                 CurrStates, /* not used */
                                 t,
                                 tree_data) == FAILURE)

                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }
                else
                {
                    if (Fix3_Sticky_LadderSwap_t(Ladder[r],
                                 /*--- RIB Ladder Reset-in-Adv payment+reset */
                                 2,
                                 0, /* No funding added here */
                                 Binary,
                                 Spread,
                                 AoM,
                                 NULL, /* not used */
                                 NULL, /* not used */
                                 CpnFloorSt,
                                 CpnCapSt,
                                 CpnStickyCoef,
                                 RibFrac,
                                 newCpnOutsSt,
                                 newCpnDcfSt,
                                 SoZ,   
                                 CompSt,
                                 NbStates,
                                 CurrStates,
                                 PrevStates,
                                 t,
                                 tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }

            } // End of for loop
            } // End of if (CplxCpnFlag && isRibAdvanceSt)

            /* Add ladder and float coupon if reset and in-arrears  */
            if (ResetFlagSt && isArrearsSt)
            for (r = 0; r < NbRibActive; r++)
            {
                /* Set percentage of Rib observations in range if applicable */
                if (NbRibObsInPer == 0)
                {
                    RibFrac = 1.;
                }
                else
                {
                    RibFrac = RibInRangeWeight  * r / NbRibObsInPer +
                              RibOutRangeWeight * (1- (double) r / NbRibObsInPer);
                }

                if (Fix3_Sticky_LadderSwap_t(Ladder[r],
                                 1,     /* Add ladder coupon payment */
                                 0,     /* No funding added in arrears */
                                 Binary,
                                 Spread,
                                 AoM,
                                 Funding,
                                 ZeroToPmtSt,
                                 CpnFloorSt,
                                 CpnCapSt,
                                 CpnStickyCoef,
                                 RibFrac,
                                 CpnOutsSt,
                                 CpnDcfSt,
                                 SoZ,   
                                 CompSt,
                                 NbStates,
                                 CurrStates,
                                 PrevStates,
                                 t,
                                 tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }

            /* Update funding for arrears reset */
            if (Fix3_Floater_t(Funding,
                          IndexFl,
                          ZeroToPmtFl,
                          ResetFlagFl&&(isArrearsFl),
                          CpnDcfFl,                                                
                          0.0,
                          ladder_data->CompFl,
                          ladder_data->ArrearsFl,
                          CpnOutsFl,
                          0.,       
                          0,   
                          t,
                          t,   /* NO DEV */
                          DCurve,
                          &dev_data,
                          tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }

            /* revise call stats, ladder reset in arrears */
            if (ladder_data->CalcStats == 'Y') 
            {
                /* interpolate */
                if (ResetFlagSt && isArrearsSt)
                {
                    for (r = 0; r < NbRibActive; r++)
                    {
                        /* Set percentage of Rib observations in range if applicable */
                        if (NbRibObsInPer == 0)
                        {
                            RibFrac = 1.;
                        }
                        else
                        {
                            RibFrac = RibInRangeWeight  * r / NbRibObsInPer +
                                      RibOutRangeWeight * (1- (double) r / NbRibObsInPer);
                        }

                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].probSlice;
                        if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                    CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                    NbStates, CurrStates, PrevStates,
                                                    t, tree_data) == FAILURE)
                            goto FREE_MEM_AND_RETURN;

                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].timeSlice;
                        if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                    CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                    NbStates, CurrStates, PrevStates,
                                                    t, tree_data) == FAILURE)
                            goto FREE_MEM_AND_RETURN;

                        for (s = 0; s < NbStates; s++)
                            xsc[r][s] = evs[r][s].tsqrSlice;
                        if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                    CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                    NbStates, CurrStates, PrevStates,
                                                    t, tree_data) == FAILURE)
                            goto FREE_MEM_AND_RETURN;

                        for (i = 0; i < (int) evs[0][0].size; i++)
                        {
                            for (s = 0; s < NbStates; s++)
                                xsc[r][s] = evs[r][s].sched[i].probSlice;
                            if(Fix3_LadderInterpolate_t(xsc[r], Binary, Spread, CpnStickyCoef,
                                                        CpnFloorSt, CpnCapSt, AoM, RibFrac,
                                                        NbStates, CurrStates, PrevStates,
                                                        t, tree_data) == FAILURE)
                                goto FREE_MEM_AND_RETURN;
                        }
                    }
                }
            }
        } /* if Current Date < Mat Date */
    }  /* END OF MAIN LOOP (for t) */

    offset = Fix3_Node_Offset(tree_data->NbFactor, 0, 0, 0, tree_data);

    /* Revise the Rib state variable once more using percentage of
     * observations in range in current period */
    if (NbPastRibObsDates > 0)
    {
        if (Fix3_Set_Slice (ObsIndex,
                            ladder_data->RibPastObsPerc * NbPastRibObsDates ,
                            t, 
                            tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }


        if (Fix3_RibLadderInterpolate_t (Ladder,
                                         NbStates,
                                         NbRibActive,
                                         ObsIndex,
                                         0,
                                         tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    /* Update notional for known fixing  */
    if (ladder_data->FixingGivenSt && (ladder_data->SoZ == 'Z') &&
        !isRibAdvanceSt)
    {
        (Ladder[0][0] + offset)[0] *= KCFactorSt;
    }

    /* Include known fixing into the ladder */
    (Ladder[0][0]  + offset)[0]  += (KnownCpnSt + offset)[0];

    /* Update funding for known fixing */
    (Funding    + offset)[0]  += (KnownCpnFl + offset)[0];

    opt_out_data->Option  = (Ladder[0][0] + offset)[0] - (Funding + offset)[0];
    
    /* Account for notional */
    opt_out_data->Option *= ladder_data->InitOuts; 

    /* Account for short option */
    if (ladder_data->LoS == 'S')
    {
        opt_out_data->Option *= -1.;
    }
    
    /* fill the rest of the stats */
    if (ladder_data->CalcStats == 'Y')
    {
        size_t cnt = 0;
        opt_out_data->prob_calc[EXER_EVENT]                    = TRUE;
        opt_out_data->prob_out_data[EXER_EVENT].TotalEventProb = Fix3EventStatsScheduleGetProbability(&evs[0][0]);
        opt_out_data->prob_out_data[EXER_EVENT].ExpEventTime   = Fix3EventStatsScheduleGetTimeExp    (&evs[0][0]);
        opt_out_data->prob_out_data[EXER_EVENT].StdEventTime   = Fix3EventStatsScheduleGetTimeStd    (&evs[0][0]);
        opt_out_data->prob_out_data[EXER_EVENT].Fugit          = Fix3EventStatsScheduleGetFugit      (&evs[0][0], 0);

        for (i=0; i<ladder_data->OptNbStats; ++i)
        {
            if (cnt >= OPT_OUT_DATA_SIZE)
                break;
            opt_out_data->prob_out_data[EXER_EVENT].EventDate[cnt] = ladder_data->OptStatDates[i];
            opt_out_data->prob_out_data[EXER_EVENT].EventProb[cnt] = 
                Fix3EventStatsScheduleGetProbabilityAsOfDate(&evs[0][0], ladder_data->OptStatDates[i]);
            ++cnt;
        }
        opt_out_data->prob_out_data[EXER_EVENT].Count = cnt;
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:
        
    if (Ladder != NULL)
    {
        for (r = 0; r <= ladder_data->MaxNbRib; r++)
        {
            if (Ladder[r] != NULL)
            {
                for (i = 0; i < ladder_data->NbStates; i++)
                {
                    if (Ladder[r][i] != NULL)
                    {
                        Fix3_Free_Slice (Ladder[r][i], tree_data);
                    }
                }
                Free_DR_Array (Ladder[r], DOUBLE_PTR, 0,
                               ladder_data->NbStates-1);
            }
        }
        Free_DR_Array (Ladder, DOUBLE_D_PTR, 0,
                       ladder_data->MaxNbRib);
    }

    Fix3_Free_Slice(KnownCpnSt,   tree_data);
    Fix3_Free_Slice(KnownCpnFl,   tree_data);
    for(i = 0; i < 2 ;i++)
        Fix3_Free_Slice(IndexSt[i],tree_data);

    Fix3_Free_Slice(IndexLadder,  tree_data);
    Fix3_Free_Slice(IndexPaySt,   tree_data);
    Fix3_Free_Slice(IndexFl,      tree_data);
    Fix3_Free_Slice(Identity,     tree_data);
    Fix3_Free_Slice(Annuity,      tree_data);
    Fix3_Free_Slice(Funding,      tree_data);
    Fix3_Free_Slice(Binary,       tree_data);
    Fix3_Free_Slice(Step,         tree_data);
    Fix3_Free_Slice(Spread,       tree_data);
    Fix3_Free_Slice(AuxSlice,     tree_data);
    Fix3_Free_Slice(Indicator,    tree_data);
    Fix3_Free_Slice(ObsIndex,     tree_data);
    Fix3_Free_Slice(IndexRib[0],  tree_data);
    Fix3_Free_Slice(IndexRib[1],  tree_data);

    for (j=0; j<3; j++) Fix3_CbkFree(&(ZeroBank[j]), tree_data);

    Fix3_Dev_Free (&dev_data, tree_data);

    if (ladder_data->CalcStats == 'Y')
    {
        for (r = 0; r <= ladder_data->MaxNbRib; r++)
        {
            for (s = 0; s < NbStates; ++s)
                Fix3EventStatsScheduleClear(&(evs[r][s]));
            free(evs[r]);
            free(xsc[r]);
        }

        free(evs);
        free(xsc);
    }
    
    if (RibInitStates)
        Free_DR_Array(RibInitStates,DOUBLE,0,NbStates-1);

    return (status);

} /* Calc_Ladder */

