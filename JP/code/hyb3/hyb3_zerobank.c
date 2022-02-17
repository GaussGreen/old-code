/****************************************************************************/
/*      Utilities to maintain a zero bank structure.                        */
/****************************************************************************/
/*      ZeroBank.c                                                          */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

/*****  Hyb3_ZbkUpdate  *****************************************************/
/*
 *      Hyb3_Dev the active elements in the zerobank
 *      add a new element if ZeroMatFlag is TRUE
 *      Returns SUCCESS or FAILURE
 */

int     Hyb3_ZbkUpdate
            (CLAIM_BANK *ZeroBank,        /* (I/O) payoff eval date list   */
             int         ZeroMatFlag,     /* (I) TRUE => add a new element */
             long        CurrentDate,     /* (I) Current date              */
             long        ErDate,          /* (I) Used for the new element  */
             int         t,               /* (I) Current time point        */
             int         T,               /* (I) Last time point           */
             int         DCurve,          /* (I) Discount curve            */
             int         DMode,           /* (I) Nb of dims of dev required */
             HYB3_DEV_DATA    *dev_data,       /* (I) Hyb3_Dev data structure        */
             HYB3_TREE_DATA   *tree_data)      /* (I) Tree data structure       */
{
    int     status = FAILURE;   /* Error status */
    TSLICE  newZero = NULL;

    /* compact and dev the zero bank */
    if (Hyb3_CbkDev(ZeroBank,
               CurrentDate,
               t,
               T,
               DCurve,
               DMode,
               dev_data,
               tree_data) == FAILURE) goto RETURN;

    if (ZeroMatFlag) /* add a new zero */
    {
        /* get a free slice from the bank */                
        if ((newZero = Hyb3_CbkPopSlice(ZeroBank)) == NULL) goto RETURN;

        if (Hyb3_Zero_t(newZero,
                   TRUE,
                   t,
                   T,
                   DCurve,
                   DMode,
                   dev_data,
                   tree_data) == FAILURE) goto RETURN;

        /* put it back in */
        if (Hyb3_CbkPushSlice(ZeroBank,
                         newZero,
                         CurrentDate,
                         ErDate) == FAILURE) goto RETURN;

    } /* if (ZeroMatFlag) */

    status = SUCCESS;

RETURN:
    
    if (status == FAILURE) DR_Error("Hyb3_ZbkUpdate: Failed!\n");
    return (status);

} /* Hyb3_ZbkUpdate */

/*****  Hyb3_FXZbkUpdate  **************************************************/
/*
 *      Hyb3_Dev the active elements in the (1.0 * FX) zerobank
 *      add a new element if ZeroMatFlag is TRUE
 *      Returns SUCCESS or FAILURE
 *      This bank discounts (1.0 * FXspot(matDate)) to account for the
 *      case where we support discounting foreign payments in domestic 
 *      equivalent curreny by multiplying by 1.0 * FX at the mat date
 *      then discounting in domestic.
 */

int     Hyb3_FXZbkUpdate
            (CLAIM_BANK *FXZeroBank,      /* (I/O) payoff eval date list   */
             int         FXZeroMatFlag,   /* (I) TRUE => add a new element */
             long        CurrentDate,     /* (I) Current date              */
             long        ErDate,          /* (I) Used for the new element  */
             int         t,               /* (I) Current time point        */
             int         T,               /* (I) Last time point           */
             int         DCurve,          /* (I) Discount curve            */
             int         DMode,           /* (I) Nb of dims of dev required */
             HYB3_DEV_DATA    *dev_data,       /* (I) Hyb3_Dev data structure        */
             HYB3_TREE_DATA   *tree_data)      /* (I) Tree data structure       */
{
    int     status = FAILURE;   /* Error status */
    TSLICE  newZero = NULL;

    /* only following modes support this type of Zbk */
    if (tree_data->TreeType != TTYPE_FX2IR)
    {
        DR_Error("Tree mode must be FX2IR to discount an FX denominated "
                 "zero bank");
        goto RETURN;
    }
    if (DMode != DISC_3D_CUPS)
    {
        DR_Error("Discount mode must be DISC_3D_CUPS to discount an FX "
                 "denominated zero bank");
        goto RETURN;
    }

    /* compact and dev the zero bank */
    if (Hyb3_CbkDev(FXZeroBank,
               CurrentDate,
               t,
               T,
               DCurve,
               DMode,
               dev_data,
               tree_data) == FAILURE) goto RETURN;

    if (FXZeroMatFlag) /* add a new zero */
    {
        /* get a free slice from the bank */                
        if ((newZero = Hyb3_CbkPopSlice(FXZeroBank)) == NULL) goto RETURN;

        if (Hyb3_FXZero_t(newZero,
                   TRUE,
                   t,
                   T,
                   DCurve,
                   DMode,
                   dev_data,
                   tree_data) == FAILURE) goto RETURN;

        /* put it back in */
        if (Hyb3_CbkPushSlice(FXZeroBank,
                         newZero,
                         CurrentDate,
                         ErDate) == FAILURE) goto RETURN;

    } /* if (FXZeroMatFlag) */

    status = SUCCESS;

RETURN:
    
    if (status == FAILURE) DR_Error("Hyb3_FXZbkUpdate: Failed!\n");
    return (status);

} /* Hyb3_FXZbkUpdate */

/*****  Hyb3_ZbkReadZero  *****************************************************/
/*
 *      Given the zero bank and a zero maturity date,
 *      returns a pointer to the zero slice. Linear interpolation on zero
 *      rate is used if interpON=TRUE, otherwise no interp is allowed.
 *      Returns a NULL pointer on failure.
 *      Note: the returned slice is for READ ONLY, do not overwrite or free.
 */

double  *Hyb3_ZbkReadZero
            (CLAIM_BANK     *ZBK,           /* (I) the zero bank           */
             long            MatDate,       /* (I) mat date of target zero */
             int             interpON,      /* (I) TRUE => interp allowed  */
             long            CurrentDate,   /* (I) date of current timept  */
             int             t,             /* (I) current timept index    */
             int             dimension,     /* (I) dimension of zero slice */
             HYB3_TREE_DATA      *tree_data)     /* (I) tree data               */
{

    /* --------------------------------------------- */

#undef  Zbk_INTERP_ZERO_FF
#define Zbk_INTERP_ZERO_FF(x)                               \
                                                            \
    if (ZeroL[x]<=TINY || ZeroR[x]<=TINY)                   \
    {                                                       \
        ZeroLocal[x] = TINY;                                \
    }                                                       \
    else                                                    \
    {                                                       \
        ZeroLocal[x] = ZeroL[x] *                           \
                       pow (ZeroR[x]/ZeroL[x], tFactor);    \
    }

#undef  Zbk_INTERP_ZERO_LC
#define Zbk_INTERP_ZERO_LC(x)                               \
                                                            \
    if (ZeroL[x]<=TINY || ZeroR[x]<=TINY)                   \
    {                                                       \
        ZeroLocal[x] = TINY;                                \
    }                                                       \
    else                                                    \
    {                                                       \
        ZRateR = pow(ZeroR[x], -365./ZDaysR) - 1.;          \
        if (ZDaysL==0) {                                    \
            ZRate = ZRateR;                                 \
        }                                                   \
        else {                                              \
            ZRateL = pow(ZeroL[x], -365./ZDaysL) - 1.;      \
            linterp(DaysToMat,                              \
                    &ZRate,                                 \
                    ZDaysL, ZDaysR,                         \
                    ZRateL, ZRateR);                        \
        }                                                   \
        ZeroLocal[x] = pow (1.+ ZRate, -DaysToMat/365.);    \
    }

    /* --------------------------------------------- */

    double *ZeroL = NULL;
    double *ZeroR = NULL;
    double *ZeroLocal = NULL;

    double  ZRate;             /* interp'd zero rate                    */
    double  ZRateL, ZRateR;    /* Zero rates for interp(ACT/365 basis)  */

    double  tFactor;           /* time factor for FlatFwd interp        */

    long    DaysToMat;         /* days from CurrentDate to MatDate      */
    long    ZDaysL, ZDaysR;

    int     Top1, Bottom1;     /* Tree limits (1rst dim)                */
    int     *Top2, *Bottom2;   /* Tree limits (2nd dim)                 */
    int     **Top3, **Bottom3; /* Tree limits (3rd dim)                 */

    int     i, j, k, idx;
    int     offset;            /* Node offset                           */

    double *output  = NULL;    /* pointer to the output zero slice      */

    long ZeroLEvDate = 0;
    long ZeroREvDate = 0;
    double *ZeroLSlice = NULL;
    double *ZeroRSlice = NULL;

    /* basic checks */
    if ((ZBK == NULL) || (tree_data == NULL)) 
        goto RETURN;
    if (ZBK->NbActiveSlices == 0) 
    {
        DR_Error("Hyb3_ZbkReadZero: Number of active claimBank slice == 0 "
                 "(ZBK->NbActiveSlices)");
        goto RETURN;
    }

    /* get the nearest idx higher or equal to MatDate */
    idx = Hyb3_CbkGetOffset(ZBK, MatDate, CbkHIGHER);
    if (idx < 0) 
    {
        DR_Error("Hyb3_ZbkReadZero: Unable to calculate greater or equal claim bank "
                 "idx from zero MatDate %ld\n", MatDate);
        goto RETURN;
    }

    /* if exact, then return the slice pointer */
    if (ZBK->EvDates[idx] == MatDate)
    {
        output = ZBK->Slices[idx];
        goto RETURN;
    }

    /* if interp not allowed or can't interp, then return NULL */
    if  (!interpON)
    {
        DR_Error("Hyb3_ZbkReadZero: No exact match and interpON is false");
        goto RETURN;
    }

    /* We can and we need to interp between idx and (idx-1 or current) */

    ZeroRSlice = ZBK->Slices[idx];
    ZeroREvDate = ZBK->EvDates[idx];

    if (idx>0)
    {
        ZeroLSlice = ZBK->Slices[idx-1];
        ZeroLEvDate = ZBK->EvDates[idx-1];
    }
    else {
        /* the other slice to interpolate with is "now" */
        ZeroLSlice = NULL; /* not needed */
        ZeroLEvDate = CurrentDate;    
    }

    /* We interp between ZeroLSlice and ZeroRSlice */

    /* initialise variables */

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    ZDaysL    = Daysact(CurrentDate, ZeroLEvDate);
    ZDaysR    = Daysact(CurrentDate, ZeroREvDate);

    DaysToMat = Daysact(CurrentDate, MatDate);
    tFactor = ((double)(DaysToMat-ZDaysL))/((double)(ZDaysR-ZDaysL));

    /*************************** 1 Factor *****************************/

    if (dimension == 1)
    {
        offset    = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        ZeroL     = ZeroLSlice   + offset;
        ZeroR     = ZeroRSlice + offset;
        ZeroLocal = ZBK->auxSlice + offset;

        switch (EslGetZeroInterpolation())
        {
        case ESL_INTERP_LINEAR: /* Linear Zero Cpn */

             for (i = Bottom1; i <= Top1; i ++)
             {
                 Zbk_INTERP_ZERO_LC(i);
             }  /* for i */
         break;

        case ESL_INTERP_FLATFWD: /* Flat Fwd */    
            for (i = Bottom1; i <= Top1; i ++)
            {
                Zbk_INTERP_ZERO_FF(i);                
            }  /* for i */
            break;

        default:
            goto RETURN;
        }
         
    }

    /*************************** 2 Factor *****************************/

    else if (dimension == 2)
    {
	switch (EslGetZeroInterpolation())
        {
        case ESL_INTERP_LINEAR: /* Linear Zero Cpn */
             for (i = Bottom1; i <= Top1; i ++)
             {
                 offset    = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                 ZeroL     = ZeroLSlice   + offset;
                 ZeroR     = ZeroRSlice + offset;
                 ZeroLocal = ZBK->auxSlice + offset;

                 for (j = Bottom2[i]; j <= Top2[i]; j++)
                 {
                     Zbk_INTERP_ZERO_LC(j);
                 }  /* for j */
             }  /* for i */
         break;
    
        case ESL_INTERP_FLATFWD: /* Flat Fwd */    
         for (i = Bottom1; i <= Top1; i ++)
             {
                 offset    = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                 ZeroL     = ZeroLSlice   + offset;
                 ZeroR     = ZeroRSlice + offset;
                 ZeroLocal = ZBK->auxSlice + offset;

                 for (j = Bottom2[i]; j <= Top2[i]; j++)
                 {
                     Zbk_INTERP_ZERO_FF(j);
                 }  /* for j */
             }  /* for i */
         break;

    default:
            goto RETURN;
        }
    }


    else if (dimension == 3)
    {
    switch (ZeroInterpTypeFlag)
        {
        case 0: /* Linear Zero Cpn */
             for (i = Bottom1; i <= Top1; i ++)
             {
                 for (j = Bottom2[i]; j <= Top2[i]; j++)
         {
                     offset    = Hyb3_Node_Offset(3, i, j, t, tree_data);
                     ZeroL     = ZBK->Slices[idx]   + offset;
                     ZeroR     = ZBK->Slices[idx+1] + offset;
                     ZeroLocal = ZBK->auxSlice + offset;

                     for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                     {
                         Zbk_INTERP_ZERO_LC(k);
             } /* for k*/
                 }  /* for j */
             }  /* for i */
         break;
    
        case 1: /* Flat Fwd */    
         for (i = Bottom1; i <= Top1; i ++)
             {
         for (j = Bottom2[i]; j <= Top2[i]; j++)
         {
                     offset    = Hyb3_Node_Offset(3, i, j, t, tree_data);
                     ZeroL     = ZBK->Slices[idx]   + offset;
                     ZeroR     = ZBK->Slices[idx+1] + offset;
                     ZeroLocal = ZBK->auxSlice + offset;
                     for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                     {
                         Zbk_INTERP_ZERO_FF(j);
             } /* for k*/
                 }  /* for j */
             }  /* for i */
         break;

    default:
            goto RETURN;
        }
    }
    else
    {
        DR_Error("Hyb3_ZbkReadZero: dimension has to be 1, 2 or 3.\n");
        goto RETURN;
    }

    output = ZBK->auxSlice;

RETURN:
    return (output);

#undef  Zbk_INTERP_ZERO

} /* Hyb3_ZbkReadZero */


/*****  Hyb3_ZbkParYield_t  ******************************************************/
/*
 *      Calculation of par yield in the tree from a zerobank. 
 *      A spread is added on top of it.
 *
 *      NOTE: This routine does NOT DO THE DEV. It's just a payoff function
 */
int   Hyb3_ZbkParYield_t(double     *ParYield,      /* (O) Par yield             */
                    double     *Annuity,       /* (O) Price of the annuity  */
                    int         YieldDim,      /* (I) Dimension of yield slice */
                    CLAIM_BANK *ZBK,           /* (I) the Zero bank         */
                    long        CurrentDate,   /* (I) Current date          */
                    long        SwapStart,     /* (I) Forward swap start    */
                    int         IndexMat,      /* (I) Index maturity in mth */
                    char        DayCount,      /* (I) Index day count       */
                    char        IndexF,        /* (I) Index payment freq    */
                    double      Spread,        /* (I) Spread                */
                    int         t,             /* (I) Current time point    */
                    HYB3_TREE_DATA  *tree_data)     /* (I) Tree data structure   */
{

    /* This is the cutoff level used to ensure stability of the par  */
    /* yield calculations. It is somewhat arbitrary, but it has been */
    /* proven  adequate for the  "problematic" JPY scenario.  A more */
    /* elegant alternative woud be to calculate the vol  (Vladimir's */
    /* approx) of the yield  in question and use  as a cut off level */
    /* the deterministic level plus the number of standard devs  set */
    /* by the user. The cost of this approach is not justifiable  in */
    /* view of the rarity of such cases.            (LP/LB April'98) */
    #define   DR_YIELD_CUTOFF   9.99

    /* --------------------------------------------- */

#undef  Zbk_CALC_PARYIELD
#define Zbk_CALC_PARYIELD(x)                                               \
    if ((AnnuityL[x]<=TINY) || (ParYieldL[x]<=TINY) || (ZPriceL[x]<=TINY)) \
    {                                                                      \
        ParYieldL[x] = DR_YIELD_CUTOFF;                                    \
    }                                                                      \
    else                                                                   \
    {                                                                      \
        ParYieldL[x] -= ZPriceL[x];                                        \
        ParYieldL[x] /= AnnuityL[x];                                       \
        ParYieldL[x] += Spread;                                            \
        if (ParYieldL[x] > DR_YIELD_CUTOFF)                                \
                ParYieldL[x] = DR_YIELD_CUTOFF;                            \
    }

    /* --------------------------------------------- */

    int     status = FAILURE;

    int     Top1, Bottom1;     /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;   /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3; /* Tree limits (3rd dim)  */
    int     i, j, k;
    int     offset;            /* Node offset            */

    long    SwapEnd;           /* end date of forward swap               */
    double  DayCntFrn;         /* day count fraction for each acc period */
    int     NbZMats = 0;       /* Nb of zero mats in the fwd swap        */
    long   *ZMats   = NULL;    /* list of zero mats in the fwd swap      */

    double *ZPrice = NULL;     /* zero slice from zerobank (don't free)  */
    double *AnnuityL;
    double *ParYieldL;
    double *ZPriceL = NULL;

    /* basic checks */
    if (SwapStart < CurrentDate) goto FREE_MEM_AND_RETURN;
    if ((ParYield == NULL) || (ZBK == NULL) || 
        (Annuity  == NULL) || (tree_data == NULL))
        goto FREE_MEM_AND_RETURN;

   

    /* calculate the zero maturity dates of the swap yield */
    SwapEnd = Nxtmth(SwapStart, (long)IndexMat, 1L);

    if (DateListFromFreq(SwapStart,
                         SwapEnd,
                         IndexF,
                         'N',
                         &NbZMats,
                         &ZMats) == FAILURE) goto FREE_MEM_AND_RETURN;

    /* set zero price at swap start */
    if (SwapStart == CurrentDate)
    {
        if (Hyb3_SetSlice(ParYield,
                      YieldDim,
                      1.0,  
                      t,
                      tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;
    }
    else
    {
        ZPrice = Hyb3_ZbkReadZero(ZBK,
                             ZMats[0],
                             TRUE,          /* interp allowed */
                             CurrentDate,
                             t,
                             YieldDim,
                             tree_data);

        if (ZPrice == NULL) goto FREE_MEM_AND_RETURN;

        if (Hyb3_CopySlice(ParYield,
                       ZPrice,
                       YieldDim,
                       t,
                       tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;
    }

    
    /* set annuity to zero  (as caller may not */
    /* have done so) and calculate the annuity */
    if (Hyb3_SetSlice(Annuity,                                       
                  YieldDim,
                  0.0,                                            
                  t,                                              
                  tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;


    for (i=1; i<NbZMats; i++)
    {
        ZPrice = Hyb3_ZbkReadZero(ZBK,
                             ZMats[i],
                             TRUE,          /* interp allowed */
                             CurrentDate,
                             t,
                             YieldDim,
                             tree_data);

        if (ZPrice == NULL) goto FREE_MEM_AND_RETURN;

        if (DrDayCountFraction(ZMats[i-1],
                               ZMats[i],
                               DayCount,
                               &DayCntFrn) == FAILURE) 
                               goto FREE_MEM_AND_RETURN;

        if (Hyb3_LCombTwoSlices(Annuity,
                           YieldDim,
                           Annuity, 1.0,
                           ZPrice,  DayCntFrn,
                           t,
                           tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;

    } /* for each zero mat i */

    /*--------------------------------------------*/
    /* Now, we can calculate the par yield values */
    /*--------------------------------------------*/

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /*************************** 1 Factor *****************************/

    if (YieldDim == 1)
    {
        offset    = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        ParYieldL = ParYield + offset;
        AnnuityL  = Annuity  + offset;
        ZPriceL   = ZPrice   + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            Zbk_CALC_PARYIELD(i);
        }  /* for i */
    }

    /*************************** 2 Factor *****************************/

    else if (YieldDim == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset    = Hyb3_Node_Offset(2, i, 0, t, tree_data);
            ParYieldL = ParYield + offset;
            AnnuityL  = Annuity  + offset;
            ZPriceL   = ZPrice   + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                Zbk_CALC_PARYIELD(j);
            }  /* for j */
        }  /* for i */
    }

    /*************************** 3 Factor *****************************/
    else if (YieldDim == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
        {
                offset    = Hyb3_Node_Offset(3, i, j, t, tree_data);
                ParYieldL = ParYield + offset;
                AnnuityL  = Annuity  + offset;
                ZPriceL   = ZPrice   + offset;
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)  
                {
                    Zbk_CALC_PARYIELD(k);
        } /* for k */
            }  /* for j */
        }  /* for i */
    }
    else 
    {
        DR_Error("Hyb3_ZbkParYield_t: dimension has to be 1, 2 or 3.\n");
        goto FREE_MEM_AND_RETURN;
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    
    Free_DR_Array(ZMats, LONG, 0, NbZMats-1);

    if (status == FAILURE)
    {
        DR_Error("Hyb3_ZbkParYield_t: Failed!\n");
    }

    return (status);

#undef  Zbk_CALC_PARYIELD

} /* Hyb3_ZbkParYield_t */


/*****  Hyb3_ZbkAnnuity_t  ******************************************************/
/*
 *      Calculation of the annuity in the tree from a zerobank. 
 *
 *      NOTE: This routine does NOT DO THE DEV. It's just a payoff function
 */
int   Hyb3_ZbkAnnuity_t(double     *Annuity,       /* (O) Price of the annuity  */
                   int         AnnuityDim,    /* (I) Dimension of annuity slice */
                   CLAIM_BANK *ZBK,           /* (I) the Zero bank         */
                   long        CurrentDate,   /* (I) Current date          */
                   long        SwapStart,     /* (I) Forward swap start    */
                   int         IndexMat,      /* (I) Index maturity in mth */
                   char        DayCount,      /* (I) Index day count       */
                   char        IndexF,        /* (I) Index payment freq    */
                   int         t,             /* (I) Current time point    */
                   HYB3_TREE_DATA  *tree_data)     /* (I) Tree data structure   */
{

    

    int     status = FAILURE;

    
    int     i;
   

    long    SwapEnd;           /* end date of forward swap               */
    double  DayCntFrn;         /* day count fraction for each acc period */
    int     NbZMats = 0;       /* Nb of zero mats in the fwd swap        */
    long   *ZMats   = NULL;    /* list of zero mats in the fwd swap      */

    double *ZPrice = NULL;     /* zero slice from zerobank (don't free)  */

    /* basic checks */
    if (SwapStart < CurrentDate) goto FREE_MEM_AND_RETURN;
    if ((Annuity == NULL) || (ZBK == NULL) || (tree_data == NULL))
        goto FREE_MEM_AND_RETURN;

   

    /* calculate the zero maturity dates of the swap yield */
    SwapEnd = Nxtmth(SwapStart, (long)IndexMat, 1L);

    if (DateListFromFreq(SwapStart,
                         SwapEnd,
                         IndexF,
                         'N',
                         &NbZMats,
                         &ZMats) == FAILURE) goto FREE_MEM_AND_RETURN;

    
    /* set annuity to zero  (as caller may not */
    /* have done so) and calculate the annuity */
    if (Hyb3_SetSlice(Annuity,                                       
                  AnnuityDim,
                  0.0,                                            
                  t,                                              
                  tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;


    for (i=1; i<NbZMats; i++)
    {
        ZPrice = Hyb3_ZbkReadZero(ZBK,
                             ZMats[i],
                             TRUE,          /* interp allowed */
                             CurrentDate,
                             t,
                             AnnuityDim,
                             tree_data);

        if (ZPrice == NULL) goto FREE_MEM_AND_RETURN;

        if (DrDayCountFraction(ZMats[i-1],
                               ZMats[i],
                               DayCount,
                               &DayCntFrn) == FAILURE) 
                               goto FREE_MEM_AND_RETURN;

        if (Hyb3_LCombTwoSlices(Annuity,
                           AnnuityDim,
                           Annuity, 1.0,
                           ZPrice,  DayCntFrn,
                           t,
                           tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;

    } /* for each zero mat i */

   

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    
    Free_DR_Array(ZMats, LONG, 0, NbZMats-1);

    if (status == FAILURE)
    {
        DR_Error("Hyb3_ZbkAnnuity_t: Failed!\n");
    }

    return (status);


} /* Hyb3_ZbkAnnuity_t */

