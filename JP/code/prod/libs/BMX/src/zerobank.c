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
#include "bmx123head.h"


/*****  ZbkUpdate  *****************************************************/
/*
 *      Dev the active elements in the zerobank
 *      add a new element if ZeroMatFlag is TRUE
 *      Returns SUCCESS or FAILURE
 */

int     ZbkUpdate
            (CLAIM_BANK *ZeroBank,        /* (I/O) payoff eval date list   */
             int         ZeroMatFlag,     /* (I) TRUE => add a new element */
             long        CurrentDate,     /* (I) Current date              */
             long        ErDate,          /* (I) Used for the new element  */
             int         t,               /* (I) Current time point        */
             int         T,               /* (I) Last time point           */
             int         DCurve,          /* (I) Discount curve            */
             DEV_DATA    *dev_data,       /* (I) Dev data structure        */
             TREE_DATA   *tree_data)      /* (I) Tree data structure       */
{
    int     status = FAILURE;   /* Error status */
    double *newZero = NULL;

    /* compact and dev the zero bank */
    if (CbkDev(ZeroBank,
               CurrentDate,
               t,
               T,
               DCurve,
               dev_data,
               tree_data) == FAILURE) goto RETURN;

    if (ZeroMatFlag) /* add a new zero */
    {
        /* get a free slice from the bank */                
        if ((newZero = CbkPopSlice(ZeroBank)) == NULL) goto RETURN;

        if (Zero_t(newZero,
                   TRUE,
                   t,
                   T,
                   DCurve,
                   dev_data,
                   tree_data) == FAILURE) goto RETURN;

        /* put it back in */
        if (CbkPushSlice(ZeroBank,
                         newZero,
                         CurrentDate,
                         ErDate) == FAILURE) goto RETURN;

    } /* if (ZeroMatFlag) */

    status = SUCCESS;

RETURN:
    
    if (status == FAILURE) DR_Error("ZbkUpdate: Failed!\n");
    return (status);

} /* ZbkUpdate */


/*****  ZbkDLFromIdx  *****************************************************/
/*
 *      Given a list of reset dates and descriptions of an index
 *      generate a list of zero maturity dates for index estimation and their
 *      corresponding usage dates
 *      Returns SUCCESS or FAILURE
 */

int     ZbkDLFromIdx
            (int        NbResetDates,   /* (I) size of reset date list  */
             long      *ResetDL,        /* (I) reset date list          */
             long      *SwapStDL,       /* (I) SwapRate start date list */
             int        IndexMat,       /* (I) index tenor in months    */
             char       IndexFreq,      /* (I) index freq               */
             int       *NbMatDates,     /* (O) Nb of zero mat dates     */
             long     **MatDL,          /* (O) zero mats                */
             int       *NbUseDates,     /* (O) Nb zero use dates        */
             long     **UseDL)          /* (O) zero use dates           */
{
    int     status = FAILURE;   /* Error status */
    int     i, j;
    long    LastZeroMat;        
    
    long   *tmpDL = NULL;       /* tmp datelist storage */
    int     NbtmpDates = 0;

    for (i=0; i<NbResetDates; i++)
    {
        if (ResetDL[i] > SwapStDL[i]) goto FREE_MEM_AND_RETURN;

        LastZeroMat = Nxtmth(SwapStDL[i], IndexMat, 1L);
            
        if (DateListFromFreq(SwapStDL[i],
                             LastZeroMat,
                             IndexFreq,
                             'N',
                             &NbtmpDates,
                             &tmpDL) == FAILURE) goto FREE_MEM_AND_RETURN;

        /* add them to the mat and use lists */
        j = (ResetDL[i] == SwapStDL[i]) ? 1 : 0;
        for (; j<NbtmpDates; j++)
        {
            if (AddDateToList(NbMatDates, MatDL, tmpDL[j]) == FAILURE) 
                        goto FREE_MEM_AND_RETURN;

            if (AddDateToList(NbUseDates, UseDL, ResetDL[i]) == FAILURE)
                        goto FREE_MEM_AND_RETURN;
        }

        /* free tmp storage for next loop */
        Free_DR_Array(tmpDL,LONG,0,NbtmpDates-1);
        tmpDL = NULL;
        NbtmpDates = 0;

    } /* for i */

    status = SUCCESS;

FREE_MEM_AND_RETURN:
    
    if (status == FAILURE) DR_Error("ZbkDLFromIdx: Failed!\n");
    Free_DR_Array(tmpDL,LONG,0,NbtmpDates-1);
    return (status);

} /* ZbkDLFromIdx */


/*****  ZbkReadZero  *****************************************************/
/*
 *      Given the zero bank and a zero maturity date,
 *      returns a pointer to the zero slice. Linear interpolation on zero
 *      rate is used if interpON=TRUE, otherwise no interp is allowed.
 *      Returns a NULL pointer on failure.
 *      Note: the returned slice is for READ ONLY, do not overwrite or free.
 */

double  *ZbkReadZero
            (CLAIM_BANK     *ZBK,           /* (I) the zero bank           */
             long            MatDate,       /* (I) mat date of target zero */
             int             interpON,      /* (I) TRUE => interp allowed  */
             long            CurrentDate,   /* (I) date of current timept  */
             int             t,             /* (I) current timept index    */
             TREE_DATA      *tree_data)     /* (I) tree data               */
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
        ZRateL = pow(ZeroL[x], -365./ZDaysL) - 1.;          \
        ZRateR = pow(ZeroR[x], -365./ZDaysR) - 1.;          \
        linterp(DaysToMat,                                  \
                &ZRate,                                     \
                ZDaysL, ZDaysR,                             \
                ZRateL, ZRateR);                            \
        ZeroLocal[x] = pow (1.+ ZRate, -DaysToMat/365.);    \
    }

    /* --------------------------------------------- */

    double *ZeroL = NULL;
    double *ZeroR = NULL;
    double *ZeroLocal = NULL;

    double  ZRate;             /* interp'd zero rate                    */
    double  ZRateL, ZRateR;    /* Zero rates for interp(ACT/365 basis)  */

    long    DaysToMat;         /* days from CurrentDate to MatDate      */
    long    ZDaysL, ZDaysR;
    
    double  tFactor;           /* time factor for FlatFwd interp */

    int     Top1, Bottom1;     /* Tree limits (1rst dim)                */
    int     *Top2, *Bottom2;   /* Tree limits (2nd dim)                 */
    int     **Top3, **Bottom3; /* Tree limits (3rd dim)                 */

    int     i, j, k, idx;
    int     offset;            /* Node offset                           */

    double *output  = NULL;    /* pointer to the output zero slice      */

    /* basic checks */
    if ((ZBK == NULL) || (tree_data == NULL)) goto RETURN;
    if (ZBK->NbActiveSlices == 0) goto RETURN;

    /* get the nearest idx lower or equal to MatDate */
    idx = CbkGetOffset(ZBK, MatDate, CbkLOWER);
    if (idx < 0) goto RETURN;

    /* if exact, then return the slice pointer */
    if (ZBK->EvDates[idx] == MatDate)
    {
        output = ZBK->Slices[idx];
        goto RETURN;
    }

    /* if interp not allowed or can't interp, then return NULL */
    if  ((!interpON) || (idx >= ZBK->NbActiveSlices-1)) goto RETURN;

    /* Now, we can and we need to interp between idx and idx+1 */

    /* initialise variables */

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    ZDaysL    = Daysact(CurrentDate, ZBK->EvDates[idx]);
    ZDaysR    = Daysact(CurrentDate, ZBK->EvDates[idx+1]);
    if ((ZDaysL == 0) || (ZDaysR == 0)) 
    {
        DR_Error("ZbkReadZero: One of the zero mats selected coincides "
                 "with current date.\n");
        goto RETURN;
    }


    DaysToMat = Daysact(CurrentDate, MatDate);    
    tFactor = ((double)(DaysToMat-ZDaysL))/((double)(ZDaysR-ZDaysL));

    /*************************** 1 Factor *****************************/

    if (tree_data->NbFactor == 1)
    {
        offset    = Node_Offset(1, 0, 0, t, tree_data);
        ZeroL     = ZBK->Slices[idx]   + offset;
        ZeroR     = ZBK->Slices[idx+1] + offset;
        ZeroLocal = ZBK->auxSlice + offset;

        switch (ZeroInterpTypeFlag)
        {
        case 0: /* Linear Zero Cpn */
            for (i = Bottom1; i <= Top1; i ++)
            {
                Zbk_INTERP_ZERO_LC(i);
            }  /* for i */
            break;
        
        case 1: /* Flat Fwd */    
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

    else if (tree_data->NbFactor == 2)
    {
        switch (ZeroInterpTypeFlag)
        {
        case 0: /* Linear Zero Cpn */
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset    = Node_Offset(2, i, 0, t, tree_data);
                ZeroL     = ZBK->Slices[idx]   + offset;
                ZeroR     = ZBK->Slices[idx+1] + offset;
                ZeroLocal = ZBK->auxSlice + offset;
            
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Zbk_INTERP_ZERO_LC(j);
                }  /* for j */
            }  /* for i */
            break;
    
        case 1: /* Flat Fwd */    
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset    = Node_Offset(2, i, 0, t, tree_data);
                ZeroL     = ZBK->Slices[idx]   + offset;
                ZeroR     = ZBK->Slices[idx+1] + offset;
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

    /*************************** 3 Factor *****************************/

    else if (tree_data->NbFactor == 3)
    {
        switch (ZeroInterpTypeFlag)
        {
        case 0: /* Linear Zero Cpn */
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset    = Node_Offset(3, i, j, t, tree_data);
                    ZeroL     = ZBK->Slices[idx]   + offset;
                    ZeroR     = ZBK->Slices[idx+1] + offset;
                    ZeroLocal = ZBK->auxSlice + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        Zbk_INTERP_ZERO_LC(k);
                    }  /* for k */
                }  /* for j */
            }  /* for i */
            break;
    
        case 1: /* Flat Fwd */    
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset    = Node_Offset(3, i, j, t, tree_data);
                    ZeroL     = ZBK->Slices[idx]   + offset;
                    ZeroR     = ZBK->Slices[idx+1] + offset;
                    ZeroLocal = ZBK->auxSlice + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        Zbk_INTERP_ZERO_FF(k);
                    }  /* for k */
                }  /* for j */
            }  /* for i */
            break;

        default:
            goto RETURN;
        }

    }  /* if then else */

    output = ZBK->auxSlice;

RETURN:

    return (output);

#undef  Zbk_INTERP_ZERO_LC
#undef  Zbk_INTERP_ZERO_FF

} /* ZbkReadZero */


/*****  ZbkBnd  *************************************************************/
/*
 *      Given the start date and the bounds date (the date the Bnds is
 *      calculated), returns a new date satisfying the following:
 *      if mode = +1:
 *          new date - start date = Bnd(bounds date)
 *      if mode = -1:
 *          start date - new date = Bnd(bounds date)
 */

long    ZbkBnd(long  Sdate,     /* (I) date the offset is calculated from */
               long  Bdate,     /* (I) date used to determine the bound   */
               long  ValueDate, /* (I) value date                         */
               int   mode)      /* (I) +1 = Fwd; -1 = Backward            */
{
    long   sign;
    double BdateInYrs;
    long   newDate;

    sign = ((mode > 0) ? +1L : -1L);
    BdateInYrs = ((double)Daysact(ValueDate, Bdate))/365.0;

    if (BdateInYrs <= 2.0)
    {
        newDate = Nxtday(Sdate, sign*31L);   /* 1 mth */
    }
    else
    if (BdateInYrs <= 5.0)
    {
        newDate = Nxtday(Sdate, sign*93L);   /* 3 mth */
    }
    else
    {
        newDate = Nxtday(Sdate, sign*186L);   /* 6 mth */
    }

    return(newDate);

} /* ZbkBnd */


/*****  ZbkG  ***********************************************************/
/*
 *      Performs the first pass of the algorithm. This routine attempts
 *      to span a set of target maturities with the fewest number of 
 *      zero-bank maturity dates.
 *      Returns the zero-bank mats (in descending order) and their
 *      corresponding earliest usage dates, together with the index of 
 *      the smallest processed Mat.
 *      Returns -999 on FAILURE
 *
 *      Only called by ZbkOptDates, do NOT call directly
 */

int     ZbkG(int          NbMat,     /* (I) size of target mats          */
             long        *Mat,       /* (I) target zero mats             */
             long        *LaMat,     /* (I) latest usage for zero mats   */
             long        *ErMat,     /* (I) earliest usage for zero mats */
             long         ValueDate, /* (I) value date                   */
             int         *OutNbZ,    /* (O) size of OutZ and ErOutZ      */
             long       **OutZ,      /* (O) zerobank mat array           */
             long       **OutErZ,    /* (O) zerobank earliest usage list */
             long       **OutFstMat) /* (O) 1st mat in each zero intval  */
{    
    int    outputN  = -999;     /* output: smallest processed Mat[] idx  */

    int    k;

    long  *ZLocal      = NULL;  /* local zerobank mat array              */
    long  *ErZLocal    = NULL;  /* local zerobank earliest usage array   */
    long  *FstMatLocal = NULL;  /* local 1st mat array                   */
    int    NbZ      = 0;        /* size of ZLocal   */
    int    NbErZ    = 0;        /* size of ErZLocal */
    int    NbFstMat = 0;        /* size of FstMat   */

    int    n;                   /* curr smallest processed Mat[] index   */
    long   Zi, ErZi;            /* curr zbank mat and earliest use dates */
    long   Zdash;               /* proposed next zbank mat */

    long   minErMat;            /* min ErMat[k] with M[k] in (Zdash, Zi) */


    /* don't proceed if the output is not empty */
    if (*OutZ   != NULL) goto FREE_MEM_AND_RETURN;
    if (*OutErZ != NULL) goto FREE_MEM_AND_RETURN;

    /* 4.4: process the last target mat */
    n    = NbMat-1;
    Zi   = Mat[n];
    ErZi = ErMat[n];

    if (AddDateToList(&NbZ,   &ZLocal,   Zi)   == FAILURE)
            goto FREE_MEM_AND_RETURN;
    if (AddDateToList(&NbErZ, &ErZLocal, ErZi) == FAILURE)
            goto FREE_MEM_AND_RETURN;
    if (AddDateToList(&NbFstMat, &FstMatLocal, Mat[n]) == FAILURE)
            goto FREE_MEM_AND_RETURN;


    while (n>0) /* there is Mat[0] to Mat[n-1] to process */
    {
        /* 4.5: find max backward jump */
        Zdash   = ZbkBnd(Zi,                            /* start date   */
                         ZbkBnd(Zi, Zi, ValueDate, -1), /* bnd date     */
                         ValueDate, 
                         -1);                           /* go backwards */

        if (Mat[n-1] < Zdash) break; /* 4.6: empty period */

        /* 4.7, 4.8: find next n, next Zi, and next ErZi */

        minErMat = 99999999L;
        for (k=n-1; k>=0; k--)
        {
            /* if Mat[k] is outside [Zdash, Zi)... */
            if (Mat[k] < Zdash) break;

            Zdash = MAX(Zdash, Nxtday(LaMat[k],1L));
            if (Mat[k] > Zdash) minErMat = MIN(minErMat, ErMat[k]);

        } /* for */

        /* prepare for next loop */
        n = k+1;
        Zi = Zdash;
        ErZi = MIN(minErMat, ErMat[n]);

        if (Mat[n] == Zi)
            ErZLocal[NbErZ-1] = MIN(minErMat, ErZLocal[NbErZ-1]);
        else
            ErZLocal[NbErZ-1] = MIN(ErZi,     ErZLocal[NbErZ-1]);

        if (AddDateToList(&NbZ,   &ZLocal,   Zi)   == FAILURE)
                goto FREE_MEM_AND_RETURN;
        if (AddDateToList(&NbErZ, &ErZLocal, ErZi) == FAILURE)
                goto FREE_MEM_AND_RETURN;
        if (AddDateToList(&NbFstMat, &FstMatLocal, Mat[n]) == FAILURE)
                goto FREE_MEM_AND_RETURN;

    } /* while */

    outputN = n;
    *OutNbZ = NbZ;
    *OutZ = ZLocal;
    *OutErZ = ErZLocal;
    *OutFstMat = FstMatLocal;

FREE_MEM_AND_RETURN:

    if (outputN == -999)
    {
        Free_DR_Array(ZLocal,      LONG, 0, NbZ);
        Free_DR_Array(ErZLocal,    LONG, 0, NbErZ);
        Free_DR_Array(FstMatLocal, LONG, 0, NbFstMat);
    }
    return (outputN);

} /* ZbkG */


/*****  ZbkH  ***********************************************************/
/*
 *      Performs the second pass of the algorithm. This routine attempts
 *      to shift the zero-bank mats from the first pass to the right, so 
 *      that they may coincide with the critical dates
 *      Returns the final zero-bank mats (in ascending order) and their
 *      corresponding earliest usage dates.
 *
 *      Only called by ZbkOptDates, do NOT call directly
 */

int     ZbkH(int       NbZ,         /* (I) size of 1st pass zeros          */
             long     *Z,           /* (I) 1st pass zero mats (desc order) */
             long     *ErZ,         /* (I) earliest usage for zero mats    */
             long     *FirstMat,    /* (I) 1st mat date in [Zi, Zi+1)      */
             int       NbC,         /* (I) size of critical dates          */
             long     *C,           /* (I) critical dates                  */
             long      ValueDate,   /* (I) value date                      */
             int      *NbMatDates,  /* (O) Nb of output zerobank mats      */
             long    **ZbkMats,     /* (O) output zerobank mats            */
             int      *NbErDates,   /* (O) Nb of output earliest usage     */
             long    **ZbkErs)      /* (O) zerobank earliest usage list    */
{

    int  status = FAILURE;

    int  i;
    long UpBound;          /* max right shift allowed               */
    int  UseMaxC;          /* TRUE => shift to a critical date      */
    int  MaxT;             /* index of the max critdate to shift to */


    /* 4.11: process the first zerobank mat */
    if (AddDateToList(NbMatDates, ZbkMats, FirstMat[NbZ-1]) == FAILURE)
            goto RETURN;
    if (AddDateToList(NbErDates, ZbkErs, ErZ[NbZ-1]) == FAILURE)
            goto RETURN;

    for (i=NbZ-2; i>=1; i--)
    {
        /* 4.12: find upper bound for right shift */

        UpBound = ZbkBnd((*ZbkMats)[*NbMatDates-1],   /* start date */
                         (*ZbkMats)[*NbMatDates-1],   /* bnd date   */
                         ValueDate, 
                         +1);                         /* go forward */

        UpBound = MIN(UpBound, Z[i-1]);
        UpBound = MIN(UpBound, FirstMat[i]);

        /* find max critical date index <= UpBound */
        MaxT = GetDLOffset(NbC, C, UpBound, CbkLOWER);
        UseMaxC = ((MaxT<0) ? FALSE : (C[MaxT] > Z[i]));

        /* 4.13: do the shift */
        if (UseMaxC)
        {
            /* shift zerobank mat to the critical date */
            if (AddDateToList(NbMatDates, ZbkMats, C[MaxT]) == FAILURE)
                    goto RETURN;
            if (AddDateToList(NbErDates, ZbkErs, ErZ[i]) == FAILURE)
                    goto RETURN;
        }
        else
        {
            /* shift to the max allowed date to allow further right shift */
            if (AddDateToList(NbMatDates, ZbkMats, UpBound) == FAILURE)
                    goto RETURN;
            if (AddDateToList(NbErDates, ZbkErs, ErZ[i]) == FAILURE)
                    goto RETURN;
        }
    } /* for */

    /* 4.15: process the last zerobank mat */
    if (NbZ > 1)
    {
        if (AddDateToList(NbMatDates, ZbkMats, Z[0]) == FAILURE)
                goto RETURN;
        if (AddDateToList(NbErDates, ZbkErs, ErZ[0]) == FAILURE)
                goto RETURN;
    }
    
    status = SUCCESS;

RETURN:

    return (status);

} /* ZbkH */


/*****  ZbkProcessDL  **********************************************/
/*
 *      This routine performs the following to a raw list of zero mat
 *      dates and the corresponding usage dates:
 *      - checks that InpUseDL[i] is strictly < InpMatDL[i]
 *      - sort the dates in MatDates ascending order
 *      - remove duplicate MatDates 
 *      - set the earliest usage date as the minimum usage dates 
 *        of all duplicates
 *      - set the latest usage date as the max usage dates of all
 *        duplicates
 *
 *      Returns SUCCESS or FAILURE
 */

int     ZbkProcessDL(int    NbInpMat,   /* (I) Nb of input mats          */
                     long  *InpMatDL,   /* (I) zero mat datelist         */
                     long  *InpUseDL,   /* (I) zero usage list           */
                     int   *NbOutMat,   /* (O) Nb of processed zero mats */
                     long **Mat,        /* (O) processed zero mats       */
                     long **LaMat,      /* (O) latest usage date list    */ 
                     long **ErMat)      /* (O) earliest usage date list  */ 
{
    int     status = FAILURE;   /* Error status = FAILURE initially */

    char    ErrorMsg[MAXBUFF];
    int     i;

    /* local datelists for the output */
    int     NbMat = 0;
    int     NbLaMat = 0;
    int     NbErMat = 0;
    long   *MatLocal   = NULL;
    long   *ErMatLocal = NULL;
    long   *LaMatLocal = NULL;

    /* ensure usage date < mat date */
    for (i=0; i<NbInpMat; i++)
    {
        if (InpUseDL[i] >= InpMatDL[i])
        {
            sprintf (ErrorMsg, 
                     "ZbkProcessDL: Usage date #%d (%ld) >= "
                     "zero maturity date (%ld)!\n",
                     i+1, InpUseDL[i], InpMatDL[i]);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }

    /* sort the date lists */
    if (SortDateList(NbInpMat, InpMatDL, InpUseDL) == FAILURE) goto RETURN;

    /* compact the list to have unique Mat */

    if (AddDateToList(&NbMat,  &MatLocal,  InpMatDL[0])==FAILURE) goto RETURN;
    if (AddDateToList(&NbLaMat,&LaMatLocal,InpUseDL[0])==FAILURE) goto RETURN;
    if (AddDateToList(&NbErMat,&ErMatLocal,InpUseDL[0])==FAILURE) goto RETURN;
    
    for (i=1; i<NbInpMat; i++)
    {
        if (InpMatDL[i] > MatLocal[NbMat-1])
        {
            /* keep this mat date by adding it to output list */
            if (AddDateToList(&NbMat,  &MatLocal,  InpMatDL[i]) == FAILURE) 
                    goto RETURN;
            if (AddDateToList(&NbLaMat,&LaMatLocal,InpUseDL[i]) == FAILURE)
                    goto RETURN;
            if (AddDateToList(&NbErMat,&ErMatLocal,InpUseDL[i]) == FAILURE)
                    goto RETURN;
        }
        else
        if (InpMatDL[i] == MatLocal[NbMat-1])
        {
            /* ignore this mat, only update the latest usage date */
            /* the earliest usage is already done due to the sort */
            LaMatLocal[NbLaMat-1] = InpUseDL[i];
        }
        else goto RETURN; /* should not get here due to the sort */

    } /* for i */

    *NbOutMat = NbMat;
    *Mat      = MatLocal;
    *LaMat    = LaMatLocal;
    *ErMat    = ErMatLocal;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("ZbkProcessDL: Failed!\n");
        Free_DR_Array(MatLocal,   LONG, 0, NbMat-1);
        Free_DR_Array(LaMatLocal, LONG, 0, NbLaMat-1);
        Free_DR_Array(ErMatLocal, LONG, 0, NbErMat-1);
    }

    return (status);

} /* ZbkProcessDL */


/*****  ZbkOrdCritDates  ****************************************************/
/*
 *      Given a list of critical dates in the CRIT_DATE struct form,
 *      returns a list of unique critical dates sorted in ascending order
 */

int     ZbkOrdCritDates
            (int         NbCritDates,   /* (I) size of CritDate struct     */
             CRIT_DATE  *CritDate,      /* (I) CritDate struct             */
             int        *NbCritDL,      /* (O) Nb of ordered crit dates    */
             long      **CritDL)        /* (O) critical date list          */
{
    int     status = FAILURE;   /* Error status */
    int     i;

    int     NbC_Local = 0;      /* size of local crit datelist            */
    long   *C_Local   = NULL;   /* local crit datelist                    */
    int     NbDummyDL = 0;      /* size of dummy datelist                 */
    long   *DummyDL   = NULL;   /* dummy datelist for use in CbkProcessDL */


    /* basic checks */
    if ((CritDL == NULL) || (NbCritDL == NULL)) goto FREE_MEM_AND_RETURN;
    if ((NbCritDates > 0) && (CritDate == NULL)) goto FREE_MEM_AND_RETURN;

    /* extract critical dates */
    for (i=0; i<NbCritDates; i++)
    {
        if (AddDateToList(&NbC_Local,
                          &C_Local, 
                          CritDate[i].CritDate
                         ) == FAILURE) goto FREE_MEM_AND_RETURN;

        if (AddDateToList(&NbDummyDL,
                          &DummyDL, 
                          CritDate[i].CritDate
                         ) == FAILURE) goto FREE_MEM_AND_RETURN;

    } /* for i */

    /* sort and merge the critical dates */
    if (CbkProcessDL(&NbC_Local,
                     &C_Local,
                     &NbDummyDL,
                     &DummyDL) == FAILURE) goto FREE_MEM_AND_RETURN;

    *NbCritDL = NbC_Local;
    *CritDL   = C_Local;
                    
    status = SUCCESS;

FREE_MEM_AND_RETURN:

    Free_DR_Array(DummyDL, LONG, 0, NbDummyDL-1);
    
    if (status == FAILURE)
    {
        Free_DR_Array(C_Local, LONG, 0, NbC_Local-1);
        DR_Error("ZbkOrdCritDates: Failed!\n");
    }

    return (status);

} /* ZbkOrdCritDates */


/*****  ZbkOptDates  ********************************************************/
/*
 *      Given a set of zero maturities/usage dates that we want to eval,
 *      find an "optimal" set of zero-bank maturities according to the
 *      algorithm. The algorithm minimises the number of zero-bank dates
 *      subject to the Bnd(.) function, and attempts to place the final date
 *      on the critical dates
 */

int     ZbkOptDates
            (int        NbMatDates,      /* (I) Nb of target mat dates      */
             long      *MatDL,           /* (I) target mat date list        */
             int        NbUseDates,      /* (I) Nb of target use dates      */
             long      *UseDL,           /* (I) target use date list        */
             int        NbCrit,          /* (I) Nb of critical dates        */
             CRIT_DATE *Crit,            /* (I) critical date structs       */
             long       ValueDate,       /* (I) value date                  */
             int       *NbZbkMats,       /* (O) Nb of zerobank mat dates    */
             long     **ZbkMats,         /* (O) zerobank mats               */
             int       *NbZbkErs,        /* (O) Nb zbank earliest use dates */
             long     **ZbkErs)          /* (O) zbank earliest use dates    */
{
    int     status = FAILURE;   /* Error status */

    int     NbCritDates = 0;
    long   *CritDates = NULL;

    int     NbMat = 0;
    long   *Mat   = NULL;
    long   *LaMat = NULL;
    long   *ErMat = NULL;

    int     NbZ      = 0;
    long   *Z        = NULL;     /* intermediate zerobank mats */
    long   *ErZ      = NULL;     /* earliest usage for Z's */
    long   *FirstMat = NULL;

    int     N; /* size of unprocessed Mat (Mat[0] .. Mat[N-1]) */
    int     n; /* smallest processed idx (Mat[n] .. Mat[N-1] is processed) */

    /* basic checks */
    if (NbMatDates != NbUseDates) goto FREE_MEM_AND_RETURN;
    if (NbMatDates == 0) return(SUCCESS);
    if ((MatDL == NULL) || (UseDL == NULL)) goto FREE_MEM_AND_RETURN;

    if ((NbCritDates > 0) && (CritDates == NULL)) goto FREE_MEM_AND_RETURN;

    if ((NbZbkMats == NULL) ||
        (NbZbkErs  == NULL) ||
        (ZbkMats == NULL)   ||
        (ZbkErs == NULL)) goto FREE_MEM_AND_RETURN;


    /* convert the critical dates struct to a list */
    if (ZbkOrdCritDates(NbCrit,
                        Crit,
                        &NbCritDates,
                        &CritDates) == FAILURE) goto FREE_MEM_AND_RETURN;
    

    /* sort/check/merge mat datelists, and calc latest usage (LaMat) */
    /* and earliest usage (ErMat) */
    if (ZbkProcessDL(NbMatDates,
                     MatDL,
                     UseDL,
                     &NbMat,
                     &Mat,
                     &LaMat,
                     &ErMat) == FAILURE) goto FREE_MEM_AND_RETURN;

    N = NbMat;

    while (N > 0) /* size of unprocessed Mat (N) is > 0 */
    {
        /* First pass: span the maturities */
        n = ZbkG(N, Mat, LaMat, ErMat, ValueDate,
                 &NbZ, &Z, &ErZ, &FirstMat);
        if (n < 0) goto FREE_MEM_AND_RETURN;

        /* now Z has covered Mat[n] .. Mat[N-1] */

        /* Second pass: shift the dates */
        if (ZbkH(NbZ, Z, ErZ, FirstMat,
                 NbCritDates, CritDates,
                 ValueDate,
                 NbZbkMats, ZbkMats,
                 NbZbkErs,  ZbkErs) == FAILURE) goto FREE_MEM_AND_RETURN;

        N=n;

        /* reset the intermediate zerobank */
        Free_DR_Array(Z,        LONG, 0, NbZ-1);
        Free_DR_Array(ErZ,      LONG, 0, NbZ-1);
        Free_DR_Array(FirstMat, LONG, 0, NbZ-1);
        Z = ErZ = FirstMat = NULL;
        NbZ = 0;

    } /* while N > 0 */

    status = SUCCESS;

FREE_MEM_AND_RETURN:
    
    if (status == FAILURE) DR_Error("ZbkOptDates: Failed!\n");

    Free_DR_Array(CritDates, LONG, 0, NbCritDates-1);

    Free_DR_Array(Mat,   LONG, 0, NbMat-1);
    Free_DR_Array(LaMat, LONG, 0, NbMat-1);
    Free_DR_Array(ErMat, LONG, 0, NbMat-1);

    Free_DR_Array(Z,        LONG, 0, NbZ-1);
    Free_DR_Array(ErZ,      LONG, 0, NbZ-1);
    Free_DR_Array(FirstMat, LONG, 0, NbZ-1);

    return (status);

} /* ZbkOptDates */


/*****  ZbkParYield_t  ******************************************************/
/*
 *      Calculation of par yield in the tree from a zerobank. 
 *      A spread is added on top of it.
 *
 *      NOTE: This routine does NOT DO THE DEV. It's just a payoff function
 */
int   ZbkParYield_t(double     *ParYield,      /* (O) Par yield             */
                    double     *Annuity,       /* (O) Price of the annuity  */
                    CLAIM_BANK *ZBK,           /* (I) the Zero bank         */
                    long        CurrentDate,   /* (I) Current date          */
                    long        SwapStart,     /* (I) Forward swap start    */
                    int         IndexMat,      /* (I) Index maturity in mth */
                    char        DayCount,      /* (I) Index day count       */
                    char        IndexF,        /* (I) Index payment freq    */
                    double      Spread,        /* (I) Spread                */
                    int         t,             /* (I) Current time point    */
                    TREE_DATA  *tree_data)     /* (I) Tree data structure   */
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
                         'F',
                         &NbZMats,
                         &ZMats) == FAILURE) goto FREE_MEM_AND_RETURN;

    /* set zero price at swap start */
    if (SwapStart == CurrentDate)
    {
        if (Set_Slice(ParYield,
                      1.0,
                      t,
                      tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;
    }
    else
    {
        ZPrice = ZbkReadZero(ZBK,
                             ZMats[0],
                             TRUE,          /* interp allowed */
                             CurrentDate,
                             t,
                             tree_data);

        if (ZPrice == NULL) goto FREE_MEM_AND_RETURN;

        if (Copy_Slice(ParYield,
                       ZPrice,
                       t,
                       tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;
    }

    
    /* set annuity to zero  (as caller may not */
    /* have done so) and calculate the annuity */
    if (Set_Slice(Annuity,                                       
                  0.0,                                            
                  t,                                              
                  tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;


    for (i=1; i<NbZMats; i++)
    {
        ZPrice = ZbkReadZero(ZBK,
                             ZMats[i],
                             TRUE,          /* interp allowed */
                             CurrentDate,
                             t,
                             tree_data);

        if (ZPrice == NULL) goto FREE_MEM_AND_RETURN;

        if (DrDayCountFraction(ZMats[i-1],
                               ZMats[i],
                               DayCount,
                               &DayCntFrn) == FAILURE) 
                               goto FREE_MEM_AND_RETURN;

        if (LCombTwoSlices(Annuity,
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

    if (tree_data->NbFactor == 1)
    {
        offset    = Node_Offset(1, 0, 0, t, tree_data);
        ParYieldL = ParYield + offset;
        AnnuityL  = Annuity  + offset;
        ZPriceL   = ZPrice   + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            Zbk_CALC_PARYIELD(i);
        }  /* for i */
    }

    /*************************** 2 Factor *****************************/

    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset    = Node_Offset(2, i, 0, t, tree_data);
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

    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset    = Node_Offset(3, i, j, t, tree_data);
                ParYieldL = ParYield + offset;
                AnnuityL  = Annuity  + offset;
                ZPriceL   = ZPrice   + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    Zbk_CALC_PARYIELD(k);
                }  /* for k */
            }  /* for j */
        }  /* for i */
    }  /* if then else */

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    
    Free_DR_Array(ZMats, LONG, 0, NbZMats-1);

    if (status == FAILURE)
    {
        DR_Error("ZbkParYield_t: Failed!\n");
    }

    return (status);

#undef  Zbk_CALC_PARYIELD

} /* ZbkParYield_t */


/*****  ZbkAnnuity_t  ******************************************************/
/*
 *      Calculation of the annuity in the tree from a zerobank. 
 *
 *      NOTE: This routine does NOT DO THE DEV. It's just a payoff function
 */
int   ZbkAnnuity_t(double     *Annuity,       /* (O) Price of the annuity  */
                   CLAIM_BANK *ZBK,           /* (I) the Zero bank         */
                   long        CurrentDate,   /* (I) Current date          */
                   long        SwapStart,     /* (I) Forward swap start    */
                   int         IndexMat,      /* (I) Index maturity in mth */
                   char        DayCount,      /* (I) Index day count       */
                   char        IndexF,        /* (I) Index payment freq    */
                   int         t,             /* (I) Current time point    */
                   TREE_DATA  *tree_data)     /* (I) Tree data structure   */
{

    

    int     status = FAILURE;

    
    int     i;
   

    long    SwapEnd;           /* end date of forward swap               */
    double  DayCntFrn;         /* day count fraction for each acc period */
    int     NbZMats = 0;       /* Nb of zero mats in the fwd swap        */
    long   *ZMats   = NULL;    /* list of zero mats in the fwd swap      */

    double *ZPrice = NULL;     /* zero slice from zerobank (don't free)  */

    /* basic checks */
    if ((Annuity == NULL) || (ZBK == NULL) || (tree_data == NULL))
        goto FREE_MEM_AND_RETURN;

   

    /* calculate the zero maturity dates of the swap yield */
    SwapEnd = Nxtmth(SwapStart, (long)IndexMat, 1L);

    if (DateListFromFreq(SwapStart,
                         SwapEnd,
                         IndexF,
                         'F',
                         &NbZMats,
                         &ZMats) == FAILURE) goto FREE_MEM_AND_RETURN;

    /* First zero date cannot be earlier than current date */
    if (ZMats[1] < CurrentDate) goto FREE_MEM_AND_RETURN;

    
    /* set annuity to zero  (as caller may not */
    /* have done so) and calculate the annuity */
    if (Set_Slice(Annuity,                                       
                  0.0,                                            
                  t,                                              
                  tree_data) == FAILURE) goto FREE_MEM_AND_RETURN;


    for (i=1; i<NbZMats; i++)
    {
        ZPrice = ZbkReadZero(ZBK,
                             ZMats[i],
                             TRUE,          /* interp allowed */
                             CurrentDate,
                             t,
                             tree_data);

        if (ZPrice == NULL) goto FREE_MEM_AND_RETURN;

        if (DrDayCountFraction(ZMats[i-1],
                               ZMats[i],
                               DayCount,
                               &DayCntFrn) == FAILURE) 
                               goto FREE_MEM_AND_RETURN;

        if (LCombTwoSlices(Annuity,
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
        DR_Error("ZbkAnnuity_t: Failed!\n");
    }

    return (status);


} /* ZbkAnnuity_t */



