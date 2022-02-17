/****************************************************************************/
/*      Calculation of a par yield from a zero coupon curve.  	            */
/****************************************************************************/
/*      PARYIELD.c                                                          */
/****************************************************************************/


/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

int Hyb3ZeroBankInterpTypeFlag = 1;

/*****  Hyb3_Par_YieldPlus_t  ****************************************************/
/**
*       Calculation of par yield in the tree. A spread is added on top of it.
*       It also outputs the annuity slice.
*       NOTE: YieldAnnuity = NULL is allowed, and in that case 
*       the annuity is not returned.
*/
int   Hyb3_Par_YieldPlus_t(
           TSLICE       ParYield,      /**< (O) Par yield                */
           TSLICE       YieldAnnuity,  /**< (O) Annuity                  */
           int          YieldDim,      /**< (I) Dimension of yield slice */
           int          NbZero,        /**< (I) Number of zeros          */
           TSLICE*      Zero,          /**< (I) Zero bank                */
           long const*  ZeroMaturity,  /**< (I) Zero maturities          */
           long         Reset,         /**< (I) Reset flag               */
           long         CurrentDate,   /**< (I) Current date             */
           long         StartDate,     /**< (I) Forward start date       */
           int          IndexMat,      /**< (I) Index maturity in months */
           char         DayCount,      /**< (I) Index Day count conv     */
           char         IndexF,        /**< (I) Index payment frequency  */
           double       Spread,        /**< (I) Spread                   */
           int          t,             /**< (I) Current time period      */
           HYB3_TREE_DATA const* tree_data) /**< (I) Tree data structure */
{

    /* This is the cutoff level used to  ensure stability of the par */
    /* yield calculations. It is somewhat arbitrary, but it has been */
    /* proven  adequate for the  "problematic" JPY scenario.  A more */
    /* elegant alternative woud be to calculate the vol  (Vladimir's */
    /* approx) of the yield  in question and use  as a cut off level */
    /* the deterministic level plus the number of standard devs  set */
    /* by the user. The cost of this approach is not justifiable  in */
    /* view of the rarity of such cases.            (LP/LB March'98) */



    #define   DR_YIELD_CUTOFF   9.99


   
    double  *ParYieldL = NULL, *YieldAnnL= NULL;    /* Local slice pointer   */
    double  *ZeroL= NULL, *ZeroR= NULL;

    double  *DayCntFrn=NULL;   /* Accrual periods                       */

    double  Annuity;           /* Annuity price                         */
    double  ZRate, ZPrice;     /* Price and rate of current zero coupon */
    double  ZRateL, ZRateR;    /* Zero rates (ACT/365 basis)            */
    double  tFactor;           /* time factor used in Flat Fwd interp   */

    long    *PmtDate=NULL;     /* Index payment dates                   */
    long    *ZDate=NULL;       /* Zero maturities in days               */

    long    AccStart;          /* Accrual start of current payment      */
    long    AccEnd;            /* Accrual end date of current payment   */

    int     *InterpIndex=NULL; /* Zero index                            */

    int     IndexFInt;         /* Index frequency as integer            */
    int     NbReset;           /* Total number of resets in the index   */

    int     CutOffLevelReached;

    int     Top1, Bottom1;     /* Tree limits (1rst dim)                */
    int     *Top2, *Bottom2;   /* Tree limits (2nd dim)                 */
    int     **Top3, **Bottom3; /* Tree limits (3rd dim)                 */

    int     i, j, k, l, m;
    int     offset;            /* Node offset                           */
    int     status = FAILURE;  /* Error status = FAILURE initially      */




    /*  If no reset, discount par yield and return. */
    if (!Reset)
    {
        return (SUCCESS);
    }


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

        
    IndexFInt = Conv_Freq (IndexF);

    NbReset = IndexMat * IndexFInt / 12;


    DayCntFrn   = (double *) DR_Array (DOUBLE, 0, NbReset);
    PmtDate     = (long *)   DR_Array (LONG,   0, NbReset);
    ZDate       = (long *)   DR_Array (LONG,   0, NbZero);
    InterpIndex = (int *)    DR_Array (INT,    0, NbReset);

    if (   (DayCntFrn   == NULL)
        || (PmtDate     == NULL)
        || (ZDate       == NULL)
        || (InterpIndex == NULL))
    {
        DR_Error("Hyb3_Par_Yield_t: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;
    }


    for (m = 0; m < NbZero; m++)
    {
        ZDate[m] = Daysact (CurrentDate, ZeroMaturity[m]);
    }


    /*
     *  Index payment dates and day count fractions.
     *  Also find index used to interpolate in zero bank.
     */

    for (l = 0; l <= NbReset; l++)
    {
        AccStart = Nxtmth (StartDate, (long) (12 * (l-1) / IndexFInt), 1L);
        AccEnd   = Nxtmth (StartDate, (long) (12 *  l    / IndexFInt), 1L);

        if (DrDayCountFraction( AccStart,
                                AccEnd,
                                DayCount,
                                &(DayCntFrn[l])) == FAILURE)
        {
            DR_Error("Hyb3_Par_Yield_t: could not calculate day count fractions!");
            goto FREE_MEM_AND_RETURN;
        }

        PmtDate[l] = Daysact (CurrentDate, AccEnd);

        m = NbZero - 2;
        while ((PmtDate[l] < ZDate[m]) && (m > 0))
            m--;

        InterpIndex[l] = m;
    }


    if (PmtDate[NbReset] > ZDate[NbZero-1])
    {        
        DR_Error("Hyb3_Par_Yield_t: not enough zeros to calculate the index (%d months). "
                 "CurrentDate = %ld, index startDate = %ld",
                 IndexMat, CurrentDate, StartDate);
        goto FREE_MEM_AND_RETURN;
    }


    switch (Hyb3ZeroBankInterpTypeFlag)
    {
        /***************************/
    case 0: /* Linear Zero Cpn */
        /***************************/

        if (YieldAnnuity != NULL)
        {
            /***   One Factor mode ***/
            if (YieldDim == 1)
            {        
                offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
                ParYieldL = (double *)ParYield     + offset;
                YieldAnnL = (double *)YieldAnnuity + offset;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = ((double **)Zero)[m]  + offset;
                        ZeroR = ((double **)Zero)[m+1] + offset;

                        /* First zero for spot starting swap */
                        if (PmtDate[l] == 0)
                        {
                            ZPrice = 1.;
                        }
                        else if (PmtDate[l] == ZDate[m])
                        {
                            if (ZeroL[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }
                            ZPrice = ZeroL[i];
                        }
                        else
                        {

                            /* The checks include the zero price and (below) */
                            /* the estimated par yield itself.               */
                            if (ZeroL[i] < TINY || ZeroR[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }

                            /* Knowing we are in valid range,perform usual calcs */
                            ZRateL = pow (ZeroL[i], - 365. / ZDate[m]) - 1.;
                            ZRateR = pow (ZeroR[i], - 365. / ZDate[m+1]) - 1.;


                            linterp(PmtDate[l],
                                &ZRate,
                                ZDate[m],
                                ZDate[m+1],
                                ZRateL,
                                ZRateR);

                            ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                        }  /* if then else */

                        if (l == 0)
                            ParYieldL[i] = ZPrice;
                        else
                            Annuity += DayCntFrn[l] * ZPrice;

                    }  /* for l */

                    if (CutOffLevelReached)
                    {
                        ParYieldL[i] = DR_YIELD_CUTOFF; /* See note above */
                        YieldAnnL[i] = 
                            (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                            / DR_YIELD_CUTOFF;
                    }
                    else
                    {
                        ParYieldL[i] -= ZPrice;
                        ParYieldL[i] /= Annuity;
                        ParYieldL[i] += Spread;

                        YieldAnnL[i] = Annuity;

                        if (ParYieldL[i] > DR_YIELD_CUTOFF)
                        {
                            ParYieldL[i] = DR_YIELD_CUTOFF;
                            YieldAnnL[i] = 
                                (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                / DR_YIELD_CUTOFF;

                        }
                    }
                }  /* for i */
            }
            /***   Two Factor mode ***/
            else if (YieldDim == 2)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    ParYieldL = (double *) ParYield     + offset;
                    YieldAnnL = (double *) YieldAnnuity + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        Annuity = 0.;
                        ZPrice  = 1.;
                        CutOffLevelReached = FALSE;

                        for (l = 0; l <= NbReset; l++)
                        {
                            m = InterpIndex[l];

                            ZeroL = ((double **)Zero)[m]   + offset;
                            ZeroR = ((double **)Zero)[m+1] + offset;

                            if (PmtDate[l] == 0)
                            {
                                ZPrice = 1.;
                            }	
                            else if (PmtDate[l] == ZDate[m])
                            {
                                if (ZeroL[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }
                                ZPrice = ZeroL[j];
                            }	
                            else
                            {
                                if (ZeroL[j] < TINY || ZeroR[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }

                                ZRateL = pow (ZeroL[j], - 365. / ZDate[m]) - 1.;
                                ZRateR = pow (ZeroR[j], - 365. / ZDate[m+1]) - 1.;

                                linterp (   PmtDate[l],
                                    &ZRate,
                                    ZDate[m],
                                    ZDate[m+1],
                                    ZRateL,
                                    ZRateR);

                                ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                            }  /* if then else */

                            if (l == 0)
                                ParYieldL[j] = ZPrice;
                            else
                                Annuity += DayCntFrn[l] * ZPrice;

                        }  /* for l */

                        if (CutOffLevelReached)
                        {
                            ParYieldL[j] = DR_YIELD_CUTOFF;
                            YieldAnnL[j] = 
                                (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                / DR_YIELD_CUTOFF;

                        }
                        else
                        {
                            ParYieldL[j] -= ZPrice;
                            ParYieldL[j] /= Annuity;
                            ParYieldL[j] += Spread;

                            YieldAnnL[j] = Annuity;

                            if (ParYieldL[j] > DR_YIELD_CUTOFF)
                            {
                                ParYieldL[j] = DR_YIELD_CUTOFF;
                                YieldAnnL[j] = (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                    / DR_YIELD_CUTOFF;

                            }
                        }
                    }  /* for j */
                }  /* for i */
            }
            else if (YieldDim == 3)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        ParYieldL = (double *) ParYield     + offset;
                        YieldAnnL = (double *) YieldAnnuity + offset;


                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            Annuity = 0.;
                            ZPrice  = 1.;
                            CutOffLevelReached = FALSE;

                            for (l = 0; l <= NbReset; l++)
                            {
                                m = InterpIndex[l];

                                ZeroL = ((double **)Zero)[m]   + offset;
                                ZeroR = ((double **)Zero)[m+1] + offset;

                                if (PmtDate[l] == 0)
                                {
                                    ZPrice = 1.;
                                }	
                                else if (PmtDate[l] == ZDate[m])
                                {
                                    if (ZeroL[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }
                                    ZPrice = ZeroL[k];
                                }	
                                else
                                {
                                    if (ZeroL[k] < TINY || ZeroR[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }

                                    ZRateL = pow (ZeroL[k], -365. / ZDate[m])-1.;
                                    ZRateR = pow (ZeroR[k], -365. / ZDate[m+1])-1.;

                                    linterp (   PmtDate[l],
                                        &ZRate,
                                        ZDate[m],
                                        ZDate[m+1],
                                        ZRateL,
                                        ZRateR);

                                    ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                                }  /* if then else */	

                                if (l == 0)
                                    ParYieldL[k] = ZPrice;
                                else
                                    Annuity += DayCntFrn[l] * ZPrice;

                            }  /* for l */

                            if (CutOffLevelReached)
                            {
                                ParYieldL[k] = DR_YIELD_CUTOFF;
                                YieldAnnL[j] = 
                                    (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                    / DR_YIELD_CUTOFF;
                            }
                            else
                            {
                                ParYieldL[k] -= ZPrice;
                                ParYieldL[k] /= Annuity;
                                ParYieldL[k] += Spread;

                                YieldAnnL[j] = Annuity;

                                if (ParYieldL[k] > DR_YIELD_CUTOFF)
                                {
                                    ParYieldL[k] = DR_YIELD_CUTOFF;
                                    YieldAnnL[j] = (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                        / DR_YIELD_CUTOFF;

                                }
                            }
                        }  /* for k */
                    }  /* for j */
                }  /* for i */
            }  /* if then else */
            else
            {
                DR_Error("ParYieldPlus_t does not support more than 3 factor tree\n");
                goto FREE_MEM_AND_RETURN;
            }
        }/* if YieldAnnuity != NULL */
        else
        {
            /***   One Factor mode ***/
            if (YieldDim == 1)
            {


                offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
                ParYieldL = (double *)ParYield + offset;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = ((double **)Zero)[m]  + offset;
                        ZeroR = ((double **)Zero)[m+1] + offset;

                        /* First zero for spot starting swap */
                        if (PmtDate[l] == 0)                                            
                        {
                            ZPrice = 1.;
                        }	    
                        else if (PmtDate[l] == ZDate[m])
                        {
                            if (ZeroL[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }
                            ZPrice = ZeroL[i];
                        }	
                        else
                        {

                            /* The checks include the zero price and (below) */
                            /* the estimated par yield itself.               */
                            if (ZeroL[i] < TINY || ZeroR[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }

                            /* Knowing we are in valid range,perform usual calcs */
                            ZRateL = pow (ZeroL[i], - 365. / ZDate[m]) - 1.;
                            ZRateR = pow (ZeroR[i], - 365. / ZDate[m+1]) - 1.;


                            linterp(PmtDate[l],
                                &ZRate,
                                ZDate[m],
                                ZDate[m+1],
                                ZRateL,
                                ZRateR);

                            ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                        }  /* if then else */	

                        if (l == 0)
                            ParYieldL[i] = ZPrice;
                        else
                            Annuity += DayCntFrn[l] * ZPrice;

                    }  /* for l */

                    if (CutOffLevelReached)
                    {
                        ParYieldL[i] = DR_YIELD_CUTOFF; /* See note above */
                    }
                    else
                    {
                        ParYieldL[i] -= ZPrice;
                        ParYieldL[i] /= Annuity;
                        ParYieldL[i] += Spread;

                        if (ParYieldL[i] > DR_YIELD_CUTOFF)
                        {
                            ParYieldL[i] = DR_YIELD_CUTOFF;
                        }
                    }

                }  /* for i */
            }
            /***   Two Factor mode ***/
            else if (YieldDim == 2)
            {


                for (i = Bottom1; i <= Top1; i ++)
                {    
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    ParYieldL = (double *)ParYield + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        Annuity = 0.;
                        ZPrice  = 1.;
                        CutOffLevelReached = FALSE;

                        for (l = 0; l <= NbReset; l++)
                        {
                            m = InterpIndex[l];

                            ZeroL = ((double **)Zero)[m]   + offset;
                            ZeroR = ((double **)Zero)[m+1] + offset;

                            if (PmtDate[l] == 0)
                            {
                                ZPrice = 1.;
                            }	
                            else if (PmtDate[l] == ZDate[m])
                            {
                                if (ZeroL[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }
                                ZPrice = ZeroL[j];
                            }	
                            else
                            {
                                if (ZeroL[j] < TINY || ZeroR[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }

                                ZRateL = pow (ZeroL[j], - 365. / ZDate[m]) - 1.;
                                ZRateR = pow (ZeroR[j], - 365. / ZDate[m+1]) - 1.;

                                linterp (   PmtDate[l],
                                    &ZRate,
                                    ZDate[m],
                                    ZDate[m+1],
                                    ZRateL,
                                    ZRateR);

                                ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                            }  /* if then else */	

                            if (l == 0)
                                ParYieldL[j] = ZPrice;
                            else
                                Annuity += DayCntFrn[l] * ZPrice;

                        }  /* for l */

                        if (CutOffLevelReached)
                        {
                            ParYieldL[j] = DR_YIELD_CUTOFF;
                        }
                        else
                        {
                            ParYieldL[j] -= ZPrice;
                            ParYieldL[j] /= Annuity;
                            ParYieldL[j] += Spread;

                            if (ParYieldL[j] > DR_YIELD_CUTOFF)
                            {
                                ParYieldL[j] = DR_YIELD_CUTOFF;
                            }
                        }
                    }  /* for j */
                }  /* for i */
            }/* if 2D */
            else if (YieldDim == 3)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        ParYieldL = (double *) ParYield + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            Annuity = 0.;
                            ZPrice  = 1.;
                            CutOffLevelReached = FALSE;

                            for (l = 0; l <= NbReset; l++)
                            {
                                m = InterpIndex[l];

                                ZeroL = ((double **)Zero)[m]   + offset;
                                ZeroR = ((double **)Zero)[m+1] + offset;

                                if (PmtDate[l] == 0)
                                {
                                    ZPrice = 1.;
                                }	
                                else if (PmtDate[l] == ZDate[m])
                                {
                                    if (ZeroL[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }
                                    ZPrice = ZeroL[k];
                                }	
                                else
                                {
                                    if (ZeroL[k] < TINY || ZeroR[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }

                                    ZRateL = pow (ZeroL[k], -365. / ZDate[m])-1.;
                                    ZRateR = pow (ZeroR[k], -365. / ZDate[m+1])-1.;

                                    linterp (   PmtDate[l],
                                                &ZRate,
                                                ZDate[m],
                                                ZDate[m+1],
                                                ZRateL,
                                                ZRateR);

                                    ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                                }  /* if then else */	

                                if (l == 0)
                                    ParYieldL[k] = ZPrice;
                                else
                                    Annuity += DayCntFrn[l] * ZPrice;

                            }  /* for l */

                            if (CutOffLevelReached)
                            {
                                ParYieldL[k] = DR_YIELD_CUTOFF;
                            }
                            else
                            {
                                ParYieldL[k] -= ZPrice;
                                ParYieldL[k] /= Annuity;
                                ParYieldL[k] += Spread;

                                if (ParYieldL[k] > DR_YIELD_CUTOFF)
                                {
                                    ParYieldL[k] = DR_YIELD_CUTOFF;
                                }
                            }
                        }  /* for k */
                    }  /* for j */
                }  /* for i */
            }  /* if then else */
            else
            {
                DR_Error("PatYieldPlus_t does not support more than 3 factors\n");
                goto FREE_MEM_AND_RETURN;
            }        
        }/* if annuity not required */
        break;

        /********************/
    case 1: /* Flat Fwd */
        /********************/

        if (YieldAnnuity != NULL)
        {
            /***   One Factor mode ***/
            if (YieldDim == 1)
            {        
                offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
                ParYieldL = (double *)ParYield     + offset;
                YieldAnnL = (double *)YieldAnnuity + offset;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = ((double **)Zero)[m]  + offset;
                        ZeroR = ((double **)Zero)[m+1] + offset;

                        /* First zero for spot starting swap */
                        if (PmtDate[l] == 0)
                        {
                            ZPrice = 1.;
                        }
                        else if (PmtDate[l] == ZDate[m])
                        {
                            if (ZeroL[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }
                            ZPrice = ZeroL[i];
                        }	
                        else
                        {

                            /*The checks include the zero price and below*/
                            /*the estimated par yield itself.            */
                            if (ZeroL[i] < TINY || ZeroR[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }

                            /* Knowing we are in valid range,
                            * perform flat fwd interp */
                            tFactor = ((double)(PmtDate[l]-ZDate[m]))
                                /((double)(ZDate[m+1]-ZDate[m]));
                            ZPrice = ZeroL[i] 
                            *  pow(ZeroR[i]/ZeroL[i], tFactor);
                        }  /* if then else */	

                        if (l == 0)
                            ParYieldL[i] = ZPrice;
                        else
                            Annuity += DayCntFrn[l] * ZPrice;

                    }  /* for l */

                    if (CutOffLevelReached)
                    {
                        ParYieldL[i] = DR_YIELD_CUTOFF; /* See note above */
                        YieldAnnL[i] = 
                            (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                            / DR_YIELD_CUTOFF;
                    }
                    else
                    {
                        ParYieldL[i] -= ZPrice;
                        ParYieldL[i] /= Annuity;
                        ParYieldL[i] += Spread;

                        YieldAnnL[i] = Annuity;

                        if (ParYieldL[i] > DR_YIELD_CUTOFF)
                        {
                            ParYieldL[i] = DR_YIELD_CUTOFF;
                            YieldAnnL[i] = 
                                (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                / DR_YIELD_CUTOFF;

                        }
                    }
                }  /* for i */
            }
            /***   Two Factor mode ***/
            else if (YieldDim == 2)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    ParYieldL = (double *) ParYield     + offset;
                    YieldAnnL = (double *) YieldAnnuity + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        Annuity = 0.;
                        ZPrice  = 1.;
                        CutOffLevelReached = FALSE;

                        for (l = 0; l <= NbReset; l++)
                        {
                            m = InterpIndex[l];

                            ZeroL = ((double **)Zero)[m]   + offset;
                            ZeroR = ((double **)Zero)[m+1] + offset;

                            if (PmtDate[l] == 0)
                            {
                                ZPrice = 1.;
                            }	
                            else if (PmtDate[l] == ZDate[m])
                            {
                                if (ZeroL[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }
                                ZPrice = ZeroL[j];
                            }	
                            else
                            {
                                if (ZeroL[j] < TINY || ZeroR[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }

                                /* Knowing we are in the valid range, 
                                * perform flat fwd interp */
                                tFactor = ((double)(PmtDate[l]-ZDate[m]))
                                    /((double)(ZDate[m+1]-ZDate[m]));
                                ZPrice = ZeroL[j] 
                                *  pow(ZeroR[j]/ZeroL[j], tFactor);

                            }  /* if then else */

                            if (l == 0)
                                ParYieldL[j] = ZPrice;
                            else
                                Annuity += DayCntFrn[l] * ZPrice;

                        }  /* for l */

                        if (CutOffLevelReached)
                        {
                            ParYieldL[j] = DR_YIELD_CUTOFF;
                            YieldAnnL[j] = 
                                (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                / DR_YIELD_CUTOFF;

                        }
                        else
                        {
                            ParYieldL[j] -= ZPrice;
                            ParYieldL[j] /= Annuity;
                            ParYieldL[j] += Spread;

                            YieldAnnL[j] = Annuity;

                            if (ParYieldL[j] > DR_YIELD_CUTOFF)
                            {
                                ParYieldL[j] = DR_YIELD_CUTOFF;
                                YieldAnnL[j] = (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                    / DR_YIELD_CUTOFF;

                            }
                        }
                    }  /* for j */
                }  /* for i */
            }
            else if (YieldDim == 3)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        ParYieldL = (double *) ParYield + offset;
                        YieldAnnL = (double *) YieldAnnuity + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            Annuity = 0.;
                            ZPrice  = 1.;
                            CutOffLevelReached = FALSE;

                            for (l = 0; l <= NbReset; l++)
                            {
                                m = InterpIndex[l];

                                ZeroL = ((double **)Zero)[m]   + offset;
                                ZeroR = ((double **)Zero)[m+1] + offset;

                                if (PmtDate[l] == 0)
                                {
                                    ZPrice = 1.;
                                }	
                                else if (PmtDate[l] == ZDate[m])
                                {
                                    if (ZeroL[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }
                                    ZPrice = ZeroL[k];
                                }	
                                else
                                {
                                    if (ZeroL[k] < TINY || ZeroR[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }

                                    /* Knowing we are in the valid range, perform flat fwd interp */
                                    tFactor = ((double)(PmtDate[l]-ZDate[m]))/((double)(ZDate[m+1]-ZDate[m]));
                                    ZPrice = ZeroL[k] *  pow(ZeroR[k]/ZeroL[k], tFactor);

                                }  /* if then else */	

                                if (l == 0)
                                    ParYieldL[k] = ZPrice;
                                else
                                    Annuity += DayCntFrn[l] * ZPrice;

                            }  /* for l */

                            if (CutOffLevelReached)
                            {
                                ParYieldL[k] = DR_YIELD_CUTOFF;
                                YieldAnnL[j] = 
                                    (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                    / DR_YIELD_CUTOFF;
                            }
                            else
                            {
                                ParYieldL[k] -= ZPrice;
                                ParYieldL[k] /= Annuity;
                                ParYieldL[k] += Spread;

                                YieldAnnL[k] = Annuity;

                                if (ParYieldL[k] > DR_YIELD_CUTOFF)
                                {
                                    ParYieldL[k] = DR_YIELD_CUTOFF;
                                    YieldAnnL[j] = (1 - pow(1.+ DR_YIELD_CUTOFF,-PmtDate[NbReset]/365.))
                                        / DR_YIELD_CUTOFF;

                                }
                            }
                        }  /* for k */
                    }  /* for j */
                }  /* for i */
            }  
            else
            {
                DR_Error("ParYieldPlus_t does not support more than 3 factor tree\n");
                goto FREE_MEM_AND_RETURN;
            }
        }/* if YieldAnnuity != NULL */
        else
        {
            /***   One Factor mode ***/
            if (YieldDim == 1)
            {
                offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
                ParYieldL = (double *)ParYield + offset;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = ((double **)Zero)[m]  + offset;
                        ZeroR = ((double **)Zero)[m+1] + offset;

                        /* First zero for spot starting swap */
                        if (PmtDate[l] == 0)
                        {
                            ZPrice = 1.;
                        }	    
                        else if (PmtDate[l] == ZDate[m])
                        {
                            if (ZeroL[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }
                            ZPrice = ZeroL[i];
                        }	
                        else
                        {

                            /*The checks include the zero price and below*/
                            /*the estimated par yield itself.            */
                            if (ZeroL[i] < TINY || ZeroR[i] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;  /* breaks out of for l */
                            }

                            /* Knowing we are in valid range,
                            * perform flat fwd interp */
                            tFactor = ((double)(PmtDate[l]-ZDate[m]))
                                /((double)(ZDate[m+1]-ZDate[m]));
                            ZPrice = ZeroL[i] 
                            *  pow(ZeroR[i]/ZeroL[i], tFactor);

                        }  /* if then else */	

                        if (l == 0)
                            ParYieldL[i] = ZPrice;
                        else
                            Annuity += DayCntFrn[l] * ZPrice;

                    }  /* for l */

                    if (CutOffLevelReached)
                    {
                        ParYieldL[i] = DR_YIELD_CUTOFF; /* See note above */
                    }
                    else
                    {
                        ParYieldL[i] -= ZPrice;
                        ParYieldL[i] /= Annuity;
                        ParYieldL[i] += Spread;

                        if (ParYieldL[i] > DR_YIELD_CUTOFF)
                        {
                            ParYieldL[i] = DR_YIELD_CUTOFF;
                        }
                    }

                }  /* for i */
            }
            /***   Two Factor mode ***/
            else if (YieldDim == 2)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    ParYieldL = (double *)ParYield + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        Annuity = 0.;
                        ZPrice  = 1.;
                        CutOffLevelReached = FALSE;

                        for (l = 0; l <= NbReset; l++)
                        {
                            m = InterpIndex[l];

                            ZeroL = ((double **)Zero)[m]   + offset;
                            ZeroR = ((double **)Zero)[m+1] + offset;

                            if (PmtDate[l] == 0)
                            {
                                ZPrice = 1.;
                            }	
                            else if (PmtDate[l] == ZDate[m])
                            {
                                if (ZeroL[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }
                                ZPrice = ZeroL[j];
                            }	
                            else
                            {
                                if (ZeroL[j] < TINY || ZeroR[j] < TINY)
                                {
                                    CutOffLevelReached = TRUE;
                                    break;
                                }

                                /* Knowing we are in the valid range, 
                                * perform flat fwd interp */
                                tFactor = ((double)(PmtDate[l]-ZDate[m]))
                                    /((double)(ZDate[m+1]-ZDate[m]));
                                ZPrice = ZeroL[j] 
                                * pow(ZeroR[j]/ZeroL[j], tFactor);

                            }  /* if then else */	

                            if (l == 0)
                                ParYieldL[j] = ZPrice;
                            else
                                Annuity += DayCntFrn[l] * ZPrice;

                        }  /* for l */

                        if (CutOffLevelReached)
                        {
                            ParYieldL[j] = DR_YIELD_CUTOFF;
                        }
                        else
                        {
                            ParYieldL[j] -= ZPrice;
                            ParYieldL[j] /= Annuity;
                            ParYieldL[j] += Spread;

                            if (ParYieldL[j] > DR_YIELD_CUTOFF)
                            {
                                ParYieldL[j] = DR_YIELD_CUTOFF;
                            }
                        }
                    }  /* for j */
                }  /* for i */
            }/* if 2D */
            else if (YieldDim == 3)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        ParYieldL = (double *)ParYield + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            Annuity = 0.;
                            ZPrice  = 1.;
                            CutOffLevelReached = FALSE;

                            for (l = 0; l <= NbReset; l++)
                            {
                                m = InterpIndex[l];

                                ZeroL = ((double **)Zero)[m]   + offset;
                                ZeroR = ((double **)Zero)[m+1] + offset;

                                if (PmtDate[l] == 0)
                                {
                                    ZPrice = 1.;
                                }	
                                else if (PmtDate[l] == ZDate[m])
                                {
                                    if (ZeroL[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }
                                    ZPrice = ZeroL[k];
                                }	
                                else
                                {
                                    if (ZeroL[k] < TINY || ZeroR[k] < TINY)
                                    {
                                        CutOffLevelReached = TRUE;
                                        break;
                                    }

                                    /* Knowing we are in the valid range, perform flat fwd interp */
                                    tFactor = ((double)(PmtDate[l]-ZDate[m]))/((double)(ZDate[m+1]-ZDate[m]));
                                    ZPrice = ZeroL[k] *  pow(ZeroR[k]/ZeroL[k], tFactor);

                                }  /* if then else */	

                                if (l == 0)
                                    ParYieldL[k] = ZPrice;
                                else
                                    Annuity += DayCntFrn[l] * ZPrice;

                            }  /* for l */

                            if (CutOffLevelReached)
                            {
                                ParYieldL[k] = DR_YIELD_CUTOFF;
                            }
                            else
                            {
                                ParYieldL[k] -= ZPrice;
                                ParYieldL[k] /= Annuity;
                                ParYieldL[k] += Spread;

                                if (ParYieldL[k] > DR_YIELD_CUTOFF)
                                {
                                    ParYieldL[k] = DR_YIELD_CUTOFF;
                                }
                            }
                        }  /* for k */
                    }  /* for j */
                }  /* for i */
            }
            else
            {
                DR_Error("PatYieldPlus_t does not support more than 3 factors\n");
                goto FREE_MEM_AND_RETURN;
            }        
        }/* if annuity not required */
        break;

        /***************************************/
    default: /* invalid ZeroBankInterpTypeFlag */
        /***************************************/

            goto FREE_MEM_AND_RETURN;

    } /* switch */
    

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (DayCntFrn, DOUBLE, 0, NbReset);
    Free_DR_Array (PmtDate, LONG, 0, NbReset);
    Free_DR_Array (ZDate, LONG, 0, NbZero);
    Free_DR_Array (InterpIndex, INT, 0, NbReset);

    return (status);

}  /* Hyb3_Par_YieldPlus_t */




/*****  Hyb3_Par_Yield_t  ********************************************************/
/**
*       Calculation of par yield in the tree. A spread is added on top of it.
*/
int   Hyb3_Par_Yield_t(
            TSLICE         ParYield,      /**< (O) Par yield                */
            int            YieldDim,      /**< (I) Dimension of yield slice */
            int            NbZero,        /**< (I) Number of zeros          */
            TSLICE*        Zero,          /**< (I) Zero bank                */
            long const*    ZeroMaturity,  /**< (I) Zero maturities          */
            long           Reset,         /**< (I) Reset flag               */
            long           CurrentDate,   /**< (I) Current date             */
            long           StartDate,     /**< (I) Forward start date       */
            int            IndexMat,      /**< (I) Index maturity in months */
            char           DayCount,      /**< (I) Index Day count conv     */
            char           IndexF,        /**< (I) Index payment frequency  */
            double         Spread,        /**< (I) Spread                   */
            int            t,	          /**< (I) Current time period      */
            HYB3_TREE_DATA const* tree_data) /**< (I) Tree data structure   */
{

    /* This is the cutoff level used to  ensure stability of the par */
    /* yield calculations. It is somewhat arbitrary, but it has been */
    /* proven  adequate for the  "problematic" JPY scenario.  A more */
    /* elegant alternative woud be to calculate the vol  (Vladimir's */
    /* approx) of the yield  in question and use  as a cut off level */
    /* the deterministic level plus the number of standard devs  set */
    /* by the user. The cost of this approach is not justifiable  in */
    /* view of the rarity of such cases.            (LP/LB March'98) */
    #define   DR_YIELD_CUTOFF   9.99


   
    double  *ParYieldL;        /* Local slice pointer                   */
    double  *ZeroL, *ZeroR;

    double  *DayCntFrn=NULL;   /* Accrual periods                       */

    double  Annuity;           /* Annuity price                         */
    double  ZRate, ZPrice;     /* Price and rate of current zero coupon */
    double  ZRateL, ZRateR;    /* Zero rates (ACT/365 basis)            */

    long    *PmtDate=NULL;     /* Index payment dates                   */
    long    *ZDate=NULL;       /* Zero maturities in days               */

    long    AccStart;          /* Accrual start of current payment      */
    long    AccEnd;            /* Accrual end date of current payment   */

    int     *InterpIndex=NULL; /* Zero index                            */

    int     IndexFInt;         /* Index frequency as integer            */
    int     NbReset;           /* Total number of resets in the index   */

    int     CutOffLevelReached;

    int     Top1, Bottom1;     /* Tree limits (1rst dim)                */
    int     *Top2, *Bottom2;   /* Tree limits (2nd dim)                 */
    int     **Top3, **Bottom3; /* Tree limits (3rd dim)                 */

    int     i, j, k, l, m;
    int     offset;            /* Node offset                           */
    int     status = FAILURE;  /* Error status = FAILURE initially      */




    /*  If no reset, discount par yield and return. */
    if (!Reset)
    {
        return (SUCCESS);
    }


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

        
    IndexFInt = Conv_Freq (IndexF);

    NbReset = IndexMat * IndexFInt / 12;


    DayCntFrn   = (double *) DR_Array (DOUBLE, 0, NbReset);
    PmtDate     = (long *)   DR_Array (LONG,   0, NbReset);
    ZDate       = (long *)   DR_Array (LONG,   0, NbZero);
    InterpIndex = (int *)    DR_Array (INT,    0, NbReset);

    /* check input pointer fields are valid */
    if (ParYield == NULL)
    {
        DR_Error("Hyb3_Par_Yield_t:  ParYield slice pointer to store result of "
                 "ParYield calculation is NULL");
        goto FREE_MEM_AND_RETURN;
    }
    if (Zero == NULL)
    {
        DR_Error("Hyb3_Par_Yield_t: input Zero bank is invalid - "
                 "NULL pointer for Zero field");
        goto FREE_MEM_AND_RETURN;
    }
    if (ZeroMaturity == NULL)
    {
        DR_Error("Hyb3_Par_Yield_t: input Zero bank is invalid - "
                 "NULL pointer for ZeroMaturity field");
        goto FREE_MEM_AND_RETURN;
    }

    if (   (DayCntFrn   == NULL)
        || (PmtDate     == NULL)
        || (ZDate       == NULL)
        || (InterpIndex == NULL))
    {
        DR_Error("Hyb3_Par_Yield_t: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;
    }


    for (m = 0; m < NbZero; m++)
    {
        ZDate[m] = Daysact (CurrentDate, ZeroMaturity[m]);
    }


    /*
     *  Index payment dates and day count fractions.
     *  Also find index used to interpolate in zero bank.
     */

    for (l = 0; l <= NbReset; l++)
    {
        AccStart = Nxtmth (StartDate, (long) (12 * (l-1) / IndexFInt), 1L);
        AccEnd   = Nxtmth (StartDate, (long) (12 *  l    / IndexFInt), 1L);

        if (DrDayCountFraction( AccStart,     
                                AccEnd,
                                DayCount,
                                &(DayCntFrn[l])) == FAILURE)
        {
            DR_Error("Hyb3_Par_Yield_t: could not calculate day count fractions!");
            goto FREE_MEM_AND_RETURN;
        }

        PmtDate[l] = Daysact (CurrentDate, AccEnd);

        m = NbZero - 2;
        while ((PmtDate[l] < ZDate[m]) && (m > 0))
            m--;

        InterpIndex[l] = m;
    }


    if (PmtDate[NbReset] > ZDate[NbZero-1])
    {        
        DR_Error("Hyb3_Par_Yield_t: not enough zeros to calculate the index "
                 "(tenor %d months).  CurrentDate = %ld, index startDate = %ld",
                 IndexMat, CurrentDate, StartDate);
        goto FREE_MEM_AND_RETURN;
    }

    /***   One Factor mode ***/
    if (YieldDim == 1)
    {
        
        
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        ParYieldL = (double *)ParYield + offset;
    
        for (i = Bottom1; i <= Top1; i ++)
        {
            Annuity = 0.;
            ZPrice  = 1.;
            CutOffLevelReached = FALSE;

            for (l = 0; l <= NbReset; l++)
            {
                m = InterpIndex[l];

                ZeroL = ((double **)Zero)[m]  + offset;
                ZeroR = ((double **)Zero)[m+1] + offset;
    
                /* First zero for spot starting swap */
                if (PmtDate[l] == 0)                                            
                {
                    ZPrice = 1.;
                }	
                else if (PmtDate[l] == ZDate[m])
                {
                    if (ZeroL[i] < TINY)
                    {
                        CutOffLevelReached = TRUE;
                        break;  /* breaks out of for l */
                    }
                    ZPrice = ZeroL[i];
                }	
                else
                {
                    
                    /* The checks include the zero price and (below) */
                    /* the estimated par yield itself.               */
                    if (ZeroL[i] < TINY || ZeroR[i] < TINY)
                    {
                        CutOffLevelReached = TRUE;
                        break;  /* breaks out of for l */
                    }

                    /* Knowing we are in the valid range, perform usual calculations */
                    ZRateL = pow (ZeroL[i], - 365. / ZDate[m]) - 1.;
                    ZRateR = pow (ZeroR[i], - 365. / ZDate[m+1]) - 1.;


                    linterp(PmtDate[l],
                            &ZRate,
                            ZDate[m],
                            ZDate[m+1],
                            ZRateL,
                            ZRateR);

                    ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                }  /* if then else */	

                if (l == 0)
                    ParYieldL[i] = ZPrice;
                else
                    Annuity += DayCntFrn[l] * ZPrice;

            }  /* for l */

            if (CutOffLevelReached)
            {
                ParYieldL[i] = DR_YIELD_CUTOFF; /* See note above */
            }
            else
            {
                ParYieldL[i] -= ZPrice;
                ParYieldL[i] /= Annuity;
                ParYieldL[i] += Spread;

                if (ParYieldL[i] > DR_YIELD_CUTOFF)
                {
                    ParYieldL[i] = DR_YIELD_CUTOFF;
                }
            }

        }  /* for i */
    }
    /***   Two Factor mode ***/
    else if (YieldDim == 2)
    {


        for (i = Bottom1; i <= Top1; i ++)
        {

            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
            ParYieldL = (double *)ParYield + offset;
    
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                Annuity = 0.;
                ZPrice  = 1.;
                CutOffLevelReached = FALSE;

                for (l = 0; l <= NbReset; l++)
                {
                    m = InterpIndex[l];

                    ZeroL = ((double **)Zero)[m]   + offset;
                    ZeroR = ((double **)Zero)[m+1] + offset;
    
                    if (PmtDate[l] == 0)
                    {
                        ZPrice = 1.;
                    }	
                    else if (PmtDate[l] == ZDate[m])
                    {
                        if (ZeroL[j] < TINY)
                        {
                            CutOffLevelReached = TRUE;
                            break;
                        }
                        ZPrice = ZeroL[j];
                    }	
                    else
                    {
                        if (ZeroL[j] < TINY || ZeroR[j] < TINY)
                        {
                            CutOffLevelReached = TRUE;
                            break;
                        }

                        ZRateL = pow (ZeroL[j], - 365. / ZDate[m]) - 1.;
                        ZRateR = pow (ZeroR[j], - 365. / ZDate[m+1]) - 1.;

                        linterp (   PmtDate[l],
                                    &ZRate,
                                    ZDate[m],
                                    ZDate[m+1],
                                    ZRateL,
                                    ZRateR);

                        ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                    }  /* if then else */	

                    if (l == 0)
                        ParYieldL[j] = ZPrice;
                    else
                        Annuity += DayCntFrn[l] * ZPrice;

                }  /* for l */

                if (CutOffLevelReached)
                {
                    ParYieldL[j] = DR_YIELD_CUTOFF;
                }
                else
                {
                    ParYieldL[j] -= ZPrice;
                    ParYieldL[j] /= Annuity;
                    ParYieldL[j] += Spread;

                    if (ParYieldL[j] > DR_YIELD_CUTOFF)
                    {
                        ParYieldL[j] = DR_YIELD_CUTOFF;
                    }
                }
            }  /* for j */
        }  /* for i */
    }
    /*** Three Factor mode ***/
    else if (YieldDim == 3)
    {

        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                ParYieldL = (double *)ParYield + offset;
    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = ((double **)Zero)[m]   + offset;
                        ZeroR = ((double **)Zero)[m+1] + offset;
    
                        if (PmtDate[l] == 0)
                        {
                            ZPrice = 1.;
                        }	
                        else if (PmtDate[l] == ZDate[m])
                        {
                            if (ZeroL[k] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;
                            }

                            ZPrice = ZeroL[k];
                        }	
                        else
                        {
                            if (ZeroL[k] < TINY || ZeroR[k] < TINY)
                            {
                                CutOffLevelReached = TRUE;
                                break;
                            }

                            ZRateL = pow (ZeroL[k], -365. / ZDate[m])-1.;
                            ZRateR = pow (ZeroR[k], -365. / ZDate[m+1])-1.;

                            linterp (   PmtDate[l],
                                        &ZRate,
                                        ZDate[m],
                                        ZDate[m+1],
                                        ZRateL,
                                        ZRateR);

                            ZPrice = pow (1. + ZRate, - PmtDate[l] / 365.);

                        }  /* if then else */	

                        if (l == 0)
                            ParYieldL[k] = ZPrice;
                        else
                            Annuity += DayCntFrn[l] * ZPrice;

                    }  /* for l */

                    if (CutOffLevelReached)
                    {
                        ParYieldL[k] = DR_YIELD_CUTOFF;
                    }
                    else
                    {
                        ParYieldL[k] -= ZPrice;
                        ParYieldL[k] /= Annuity;
                        ParYieldL[k] += Spread;

                        if (ParYieldL[k] > DR_YIELD_CUTOFF)
                        {
                            ParYieldL[k] = DR_YIELD_CUTOFF;
                        }
                    }
                }  /* for k */
            }  /* for j */
        }  /* for i */
    }  /* if then else */


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (DayCntFrn, DOUBLE, 0, NbReset);
    Free_DR_Array (PmtDate, LONG, 0, NbReset);
    Free_DR_Array (ZDate, LONG, 0, NbZero);
    Free_DR_Array (InterpIndex, INT, 0, NbReset);

    return (status);

}  /* Hyb3_Par_Yield_t */


