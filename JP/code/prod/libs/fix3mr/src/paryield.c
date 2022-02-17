/****************************************************************************/
/*      Calculation of par yield from zero bank.                            */
/****************************************************************************/
/*      PARYIELD.c                                                          */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/paryield.c,v 1.10 2003/09/12 17:31:47 dfung Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  Par_Yield_Plus ******************************************************/
/*
*       Calculate par 4orward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int Par_Yield_Plus(
	double				*ParYield,			/* (O) Par forward yield          */
    double				*Annuity,			/* (O) Annuity                    */
    int					NbZero,				/* (I) Number of zeros            */
    double				*Zero,				/* (I) Zero rates                 */
    long				*ZeroDate,			/* (I) Zero maturity dates        */
    long				CurrentDate,		/* (I) Current date               */
    long				StartDate,			/* (I) Forward start date         */
	char				IndexType,			/* (I) Rate type (cash or swap)	  */
	int					IndexMat,			/* (I) Index maturity             */
	char				DayCount,			/* (I) Index Day count convention */
	char				IndexF)				/* (I) Index payment frequency    */
{

    double				*DayCntFrn=NULL;	/* Accrual periods                         */
    double				*ZeroToPmt=NULL;	/* Zero coupons maturing on index payments */
											
    long				*PmtDate=NULL;		/* Index payment dates                     */
											
    long				AccStart;			/* Accrual start of current payment        */
    long				AccEnd;				/* Accrual end date of current payment     */
											
    int					IndexFInt;			/* Index frequency as integer              */
    int					NbReset=0;			/* Total number of resets in the index     */
    int					j;					/* Loop index */
    int					status = FAILURE;	/* Error status = FAILURE initially        */

	double				DFstart;			/* Discount bond to start date */
	double				DFend;				/* Discount bond to end date */
	long				EndDate;			/* Rate end date */
	double				DCF;				/* Cash rate accrual */

	/* Date ordering */
    if (StartDate < CurrentDate)                            
    {        
        DR_Error("Par_Yield_Plus: start date falls before current date!");
        goto FREE_MEM_AND_RETURN;
    }
	
	/* Compute the rate */
	if (IndexType == 'C')
	{
		/* ------------------------------------------------------- */
		/* CASH RATE	                                           */
		/* ------------------------------------------------------- */

		/* End date */
		EndDate = Nxtmth (StartDate, IndexMat, 1L);

		/* Price the zero coupon bond to StartDate */
        DFstart = ZeroPrice(StartDate,			/* (I) Discount bond maturity     */
                            CurrentDate,	   	/* (I) Value Date                 */
                            NbZero,				/* (I) Number of zeros            */
                            ZeroDate,           /* (I) Zero maturity dates        */
                            Zero);   		   	/* (I) Zero rates                 */
        if (DFstart < 0.0) goto FREE_MEM_AND_RETURN;

		/* Price the zero coupon bond to EndDate */
        DFend   = ZeroPrice(EndDate,			/* (I) Discount bond maturity     */
                            CurrentDate,	   	/* (I) Value Date                 */
                            NbZero,				/* (I) Number of zeros            */
                            ZeroDate,           /* (I) Zero maturity dates        */
                            Zero);   		   	/* (I) Zero rates                 */
        if (DFend < 0.0) goto FREE_MEM_AND_RETURN;

		/* Accrual */
		if (DrDayCountFraction(StartDate, EndDate, DayCount, &DCF) == FAILURE)
		{
			DR_Error("Par_Yield_Plus: Failed to compute cash accrual!");
			goto FREE_MEM_AND_RETURN;
		}

		/* Rate */
		*ParYield = (DFstart/DFend - 1.0) / DCF;
		*Annuity  = DFend * DCF;
	}
	else if (IndexType == 'S')
	{
		/* ------------------------------------------------------- */
		/* SWAP RATE	                                           */
		/* ------------------------------------------------------- */

		/* Number of coupons per year */
		IndexFInt = Conv_Freq (IndexF);

		/* Nb of reset: 1 for Libor. If IndexMat is not an integer */
		/* number of resets it will round it down.                 */
		NbReset = IndexMat * IndexFInt / 12;

		/* Allocate space */
		DayCntFrn = (double *) DR_Array (DOUBLE, 0, NbReset);
		ZeroToPmt = (double *) DR_Array (DOUBLE, 0, NbReset);
		PmtDate   = (long *)   DR_Array (LONG,   0, NbReset);

		/* Check for memory allcoation failure */
		if (   (DayCntFrn == NULL)
			|| (ZeroToPmt == NULL)
			|| (PmtDate   == NULL))
		{
			DR_Error("Par_Yield_Plus: could not allocate memory!");
			goto FREE_MEM_AND_RETURN;
		}

		/*
		 *  Index payment dates and day count fractions.
		 */

		for (j = 0; j <= NbReset; j++)
		{
			AccStart = Nxtmth (StartDate, (long) (12 * (j-1) / IndexFInt), 1L);
			AccEnd   = Nxtmth (StartDate, (long) (12 *  j    / IndexFInt), 1L);
                            
			PmtDate[j] = AccEnd;

			if (DrDayCountFraction( AccStart,
									AccEnd,
									DayCount,
									&(DayCntFrn[j])) == FAILURE)
			{
				DR_Error("Par_Yield_Plus: could not calculate day count fractions!");
				goto FREE_MEM_AND_RETURN;
			}
		}  /* for j */


		/* 
		 *  Interpolate zero bank to find zeros to index payments.
		 */

		for (j = 0; j <= NbReset; j++)
		{
			/* Price the zero coupon bond to that date */
            ZeroToPmt[j] = ZeroPrice(PmtDate[j],	/* (I) Discount bond maturity     */
                                     CurrentDate,	/* (I) Value Date                 */
                                     NbZero,		/* (I) Number of zeros            */
                                     ZeroDate,      /* (I) Zero maturity dates        */
                                     Zero);   		/* (I) Zero rates                 */
            if (ZeroToPmt[j] < 0.0) goto FREE_MEM_AND_RETURN;
		}


		/*
		 *  Calculate annuity price and par yield.
		 */

		*Annuity = 0.;

		for (j = 1; j <= NbReset; j++)
		{
			*Annuity += DayCntFrn[j] * ZeroToPmt[j];

		}  /* for j */
        
		*ParYield = (ZeroToPmt[0] - ZeroToPmt[NbReset]) / *Annuity;
	}
    else
	{
		DR_Error("Par_Yield_Plus: Rate type must be 'C'ash or 'S'wap.");
		goto FREE_MEM_AND_RETURN;        
	}

	/* If we get this far, it worked */
    status = SUCCESS;

FREE_MEM_AND_RETURN:

	/* Release memory */
    Free_DR_Array (DayCntFrn, DOUBLE, 0, NbReset);
    Free_DR_Array (ZeroToPmt, DOUBLE, 0, NbReset);
    Free_DR_Array (PmtDate, LONG, 0, NbReset);

	/* Done */
    return (status);
}  /* Par_Yield_Plus */





/*****  Par_Yield  **********************************************************/
/*
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int    Par_Yield (double      *ParYield,   /* (O) Par forward yield          */
                  double      *Annuity,    /* (O) Annuity                    */
                  int         NbZero,      /* (I) Number of zeros            */
                  double      *Zero,       /* (I) Zero rates                 */
                  long        *ZeroDate,   /* (I) Zero maturity dates        */
                  long        CurrentDate, /* (I) Current date               */
                  long        StartDate,   /* (I) Forward start date         */
                  int         IndexMat,    /* (I) Index maturity             */
                  char        DayCount,    /* (I) Index Day count convention */
                  char        IndexF)      /* (I) Index payment frequency    */
{
	/* Use the new function in swap-rate mode */
	return Par_Yield_Plus(
					ParYield,			/* (O) Par forward yield          */
					Annuity,			/* (O) Annuity                    */
					NbZero,				/* (I) Number of zeros            */
					Zero,				/* (I) Zero rates                 */
					ZeroDate,			/* (I) Zero maturity dates        */
					CurrentDate,		/* (I) Current date               */
					StartDate,			/* (I) Forward start date         */
					'S',				/* (I) Rate type (cash or swap)	  */
					IndexMat,			/* (I) Index maturity             */
					DayCount,			/* (I) Index Day count convention */
					IndexF);			/* (I) Index payment frequency    */
}  /* Par_Yield */





/*****  Par_Yield_From_Dates  ************************************************/
/*
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 */

 int    Par_Yield_From_Dates
                (double  *ParYield,     /* (O) Par forward yield             */
                 double  *Annuity,      /* (O) Annuity                       */
                 long    SwapSt,        /* (I) Underlying swap start         */
                 long    SwapMat,       /* (I) Underlying swap maturity      */
                 char    DCC,           /* (I) Underlying day count conv.    */
                 char    Freq,          /* (I) Underlying frequency          */
                 char    StubConv,      /* (I) F, B or N(one allowed)        */
                 int     NbZero,        /* (I) Number of zeros               */
                 double  *Zero,         /* (I) Zero rates                    */
                 long    *ZeroDate,     /* (I) Zero maturity dates           */
                 long    BaseDate)      /* (I) Zero curve base date          */
{
    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmts */
    long        CpnPmtDate;           /* Current cpn pmt date    */
    long        TempDates[2];         /* To construct cpn list   */

    double  DCCFrac;             /* Current coupon day count fraction  */
    double  AnnuityL=0.0;        /* Local variable for annuity price   */
    double  ParYieldL;           /* Local variable for forward yield   */

    double  ZerotoS;             /* Zero to swap start                 */
    double  ZerotoCpn=0.0;       /* Zero to current cpn                */


    int     j;
    int     status = FAILURE;    /* Error status = FAILURE initially   */

    TempDates[0] = SwapSt;
    TempDates[1] = SwapMat;

    CpnEventList = DrNewEventListFromFreq(2,
                                          TempDates,
                                          Freq,
                                          StubConv,
                                          'N',   /* Dates in not required */
                                          NULL, NULL, NULL, NULL, NULL);
    if (CpnEventList == NULL)
    {
        goto FREE_MEM_AND_RETURN;
    }


    /* Zero to swap start */
    ZerotoS = ZeroPrice(SwapSt,	        /* (I) Discount bond maturity     */
                        BaseDate,	    /* (I) Value Date                 */
                        NbZero,		    /* (I) Number of zeros            */
                        ZeroDate,       /* (I) Zero maturity dates        */
                        Zero);   		/* (I) Zero rates                 */
    if (ZerotoS < 0.0) goto FREE_MEM_AND_RETURN;

    for (j = 1; j < CpnEventList->NbEntries; j++)
    {

        CpnPmtDate = CpnEventList->Dates[j];

        if (DrDayCountFraction (CpnEventList->Dates[j-1],
                                CpnEventList->Dates[j], 
                                DCC,
                                &DCCFrac) == FAILURE)
        {
            DR_Error ("Could not calculate day count fraction "
                      "(Par_Yield_From_Dates)!");
            goto FREE_MEM_AND_RETURN;
        }


        /* Zero to current cpn */
        ZerotoCpn = ZeroPrice(CpnPmtDate,	/* (I) Discount bond maturity     */
                              BaseDate,	    /* (I) Value Date                 */
                              NbZero,		/* (I) Number of zeros            */
                              ZeroDate,     /* (I) Zero maturity dates        */
                              Zero);   		/* (I) Zero rates                 */
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;
        
        AnnuityL += DCCFrac * ZerotoCpn;

    }  /* for j */

    if (AnnuityL < 0.0)
    {
        DR_Error("Annuity is < 0.0. (Par_Yield_From_Dates)");
        goto FREE_MEM_AND_RETURN;
    }
    ParYieldL = (ZerotoS - ZerotoCpn) / AnnuityL;

    /* Transfer values to output */
    *ParYield = ParYieldL;
    *Annuity  = AnnuityL;


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    DrFreeEventList(CpnEventList);

    return (status);

}  /* Par_Yield_From_Dates */





/*****  Par_Yield_t  ********************************************************/
/*
*       Calculation of par yield in the tree. A spread is added on top of it.
*/
int     Par_Yield_t (double      *ParYield,     /* (O) Par yield             */
                     int         NbZero,        /* (I) Number of zeros       */
                     double      **Zero,        /* (I) Zero bank             */
                     long        *ZeroMaturity, /* (I) Zero maturities       */
                     long        Reset,         /* (I) Reset flag            */
                     long        CurrentDate,   /* (I) Current date          */
                     long        StartDate,     /* (I) Forward start date    */
                     int         IndexMat,      /* (I) Index maturity in mth */
                     char        DayCount,      /* (I) Index day count       */
                     char        IndexF,        /* (I) Index payment freq    */
                     double      Spread,        /* (I) Spread                */
                     int         t,             /* (I) Current time point    */
                     int         T,             /* (I) Last time point       */
                     int         DCurve,        /* (I) Discount curve        */
                     DEV_DATA    *dev_data,     /* (I) Dev data structure    */
                     TREE_DATA   *tree_data)    /* (I) Tree data structure   */
{

    /* This is the cutoff level used to  ensure stability of the par */
    /* yield calculations. It is somewhat arbitrary, but it has been */
    /* proven  adequate for the  "problematic" JPY scenario.  A more */
    /* elegant alternative woud be to calculate the vol  (Vladimir's */
    /* approx) of the yield  in question and use  as a cut off level */
    /* the deterministic level plus the number of standard devs  set */
    /* by the user. The cost of this approach is not justifiable  in */
    /* view of the rarity of such cases.            (LP/LB April'98) */
    #define   DR_YIELD_CUTOFF   9.99


   
    double  *ParYieldL;        /* Local slice pointer                   */
    double  *ZeroL, *ZeroR;

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


    /*
     *  In no reset, discount par yield and return.
     */

    if (!Reset)
    {
        if (Dev (   ParYield,
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            return (FAILURE);
        }

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
        DR_Error("Par_Yield_t: could not allocate memory!");
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
            DR_Error("Par_Yield_t: could not calculate day count fractions!");
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
        DR_Error("Par_Yield_t: not enough zeros to calculate the index!");
        goto FREE_MEM_AND_RETURN;
    }


    switch (ZeroInterpTypeFlag)
    {
    /***************************/
    case 0: /* Linear Zero Cpn */
    /***************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            ParYieldL = ParYield + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {
                Annuity = 0.;
                ZPrice  = 1.;
                CutOffLevelReached = FALSE;

                for (l = 0; l <= NbReset; l++)
                {
                    m = InterpIndex[l];

                    ZeroL = Zero[m]   + offset;
                    ZeroR = Zero[m+1] + offset;
    
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

                        linterp (   PmtDate[l],
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
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                ParYieldL = ParYield + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = Zero[m]   + offset;
                        ZeroR = Zero[m+1] + offset;

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
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    ParYieldL = ParYield + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        Annuity = 0.;
                        ZPrice  = 1.;
                        CutOffLevelReached = FALSE;

                        for (l = 0; l <= NbReset; l++)
                        {
                            m = InterpIndex[l];

                            ZeroL = Zero[m]   + offset;
                            ZeroR = Zero[m+1] + offset;
    
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
        break;
    
    /********************/
    case 1: /* Flat Fwd */
    /********************/
    
        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            ParYieldL = ParYield + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {
                Annuity = 0.;
                ZPrice  = 1.;
                CutOffLevelReached = FALSE;

                for (l = 0; l <= NbReset; l++)
                {
                    m = InterpIndex[l];

                    ZeroL = Zero[m]   + offset;
                    ZeroR = Zero[m+1] + offset;
    
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

                        /* Knowing we are in the valid range, perform flat fwd interp */
                        tFactor = ((double)(PmtDate[l]-ZDate[m]))/((double)(ZDate[m+1]-ZDate[m]));
                        ZPrice = ZeroL[i] *  pow(ZeroR[i]/ZeroL[i], tFactor);

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
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                ParYieldL = ParYield + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Annuity = 0.;
                    ZPrice  = 1.;
                    CutOffLevelReached = FALSE;

                    for (l = 0; l <= NbReset; l++)
                    {
                        m = InterpIndex[l];

                        ZeroL = Zero[m]   + offset;
                        ZeroR = Zero[m+1] + offset;

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

                            /* Knowing we are in the valid range, perform flat fwd interp */
                            tFactor = ((double)(PmtDate[l]-ZDate[m]))/((double)(ZDate[m+1]-ZDate[m]));
                            ZPrice = ZeroL[j] *  pow(ZeroR[j]/ZeroL[j], tFactor);

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
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    ParYieldL = ParYield + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        Annuity = 0.;
                        ZPrice  = 1.;
                        CutOffLevelReached = FALSE;

                        for (l = 0; l <= NbReset; l++)
                        {
                            m = InterpIndex[l];

                            ZeroL = Zero[m]   + offset;
                            ZeroR = Zero[m+1] + offset;
    
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
        }  /* if then else */
        break;

    /***************************************/
    default: /* invalid ZeroInterpTypeFlag */
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

}  /* Par_Yield_t */


/*****  ParYieldRatio  ****************************************************/
/*
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */

int  ParYieldRatio(double     *Ratio,      /* (O) Yield1+sprd / Yield2+sprd  */
                   long        StartDate1, /* (I) start date of yield 1      */
                   long        StartDate2, /* (I) start date of yield 2      */
                   double      Spread,     /* (I) added to both yields       */
                   int         NbZero,     /* (I) Number of zeros            */
                   long        *ZeroDates, /* (I) Zero maturity dates        */
                   double      *ZeroRates, /* (I) Zero rates                 */
                   long        ValueDate,  /* (I) Value date of zero curve   */
                   int         IndexMat,   /* (I) Index maturity             */
                   char        DayCount,   /* (I) Index Day count convention */
                   char        IndexF)     /* (I) Index payment frequency    */
{
    int     status = FAILURE;
    double  Yield1,   Yield2;
    double  Annuity1, Annuity2;

    if (Par_Yield(&Yield1,
                  &Annuity1,
                  NbZero,
                  ZeroRates,
                  ZeroDates,
                  ValueDate,
                  StartDate1,
                  IndexMat,
                  DayCount,
                  IndexF) == FAILURE) goto RETURN;

    if (Par_Yield(&Yield2,
                  &Annuity2,
                  NbZero,
                  ZeroRates,
                  ZeroDates,
                  ValueDate,
                  StartDate2,
                  IndexMat,
                  DayCount,
                  IndexF) == FAILURE) goto RETURN;

    *Ratio = (Yield1 + Spread)/(Yield2 + Spread);

    status = SUCCESS;

RETURN:

    if (status == FAILURE) DR_Error("ParYieldRatio: Failed!");

    return(status);

}/* ParYieldRatio */

