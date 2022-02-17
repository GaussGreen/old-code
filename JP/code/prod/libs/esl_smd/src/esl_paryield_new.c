/****************************************************************************/
/*      Calculation of par yield from zero bank.                            */
/****************************************************************************/
/*      ESL_PARYIELD.c                                                      */
/****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esl_paryield.h"
#include "esl_zeros.h"
#include "esl_alloc.h"
#include "esl_util.h"
#include "esl_error.h"


/*****  Par_Yield_Plus ******************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int Par_Yield_Plus(
                   double*      ParYield      /** (O) Par forward yield          */
                  ,double*      Annuity       /** (O) Annuity                    */
                  ,const T_CURVE  *zc   /** (I) Zero curve            */
                  ,long         StartDate     /** (I) Forward start date         */
                  ,char         IndexType     /** (I) Rate type (cash or swap)   */
                  ,int          IndexMat      /** (I) Index maturity             */
                  ,char         DayCount      /** (I) Index Day count convention */
                  ,char         IndexF        /** (I) Index payment frequency    */
		)
{
    double *DayCntFrn=NULL;     /* Accrual periods                         */
    double *ZeroToPmt=NULL;     /* Zero coupons maturing on index payments */
    long   *PmtDate=NULL;       /* Index payment dates                     */
											
    long AccStart;              /* Accrual start of current payment        */
    long AccEnd;                /* Accrual end date of current payment     */
											
    int IndexFInt;              /* Index frequency as integer              */
    int NbReset=0;              /* Total number of resets in the index     */
    int j;                      /* Loop index                              */
    int status = FAILURE;       /* Error status = FAILURE initially        */

    double DFstart;             /* Discount bond to start date */
    double DFend;               /* Discount bond to end date */
    long EndDate;               /* Rate end date */
    double DCF;                 /* Cash rate accrual */

    /* Base date of zc */
    long  CurrentDate = zc->baseDate;


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
        DFstart = GetZeroPrice(StartDate, zc);   	   	/* (I) Zero curve                 */
        if (DFstart < 0.0) goto FREE_MEM_AND_RETURN;

        /* Price the zero coupon bond to EndDate */
        DFend   = GetZeroPrice(EndDate, zc);   	   	/* (I) Zero curve                 */
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
		    ZeroToPmt[j] = GetZeroPrice(PmtDate[j], zc);
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

    return (status);
}  /* Par_Yield_Plus */


/*****  Par_Yield  **********************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int    Par_Yield (double*       ParYield   /** (O) Par forward yield         */
                 ,double*       Annuity    /** (O) Annuity                   */
                 ,const T_CURVE  *zc /** (I) Zero curve                */
                 ,long          StartDate   /** (I) Forward start date        */
                 ,int           IndexMat    /** (I) Index maturity            */
                 ,char          DayCount    /** (I) Index Day count convention*/
                 ,char          IndexF      /** (I) Index payment frequency   */
		)
{
	/* Use the new function in swap-rate mode */
	return Par_Yield_Plus(
			ParYield,			/* (O) Par forward yield          */
			Annuity,			/* (O) Annuity                    */
			zc,				/* (I) Zero curve                 */
			StartDate,			/* (I) Forward start date         */
			'S',				/* (I) Rate type (cash or swap)	  */
			IndexMat,			/* (I) Index maturity             */
			DayCount,			/* (I) Index Day count convention */
			IndexF);			/* (I) Index payment frequency    */

}  /* Par_Yield */

int    ParYld (double*          ParYield   /** (O) Par forward yield         */
              ,double*          Annuity    /** (O) Annuity                   */
              ,T_CURVE const*   crv        /** (I) index zero curve          */
              ,long             StartDate  /** (I) Forward start date        */
              ,int              IndexMat   /** (I) Index maturity            */
              ,char             DayCount   /** (I) Index Day count convention*/
              ,char             IndexF     /** (I) Index payment frequency   */
		)
{
    return Par_Yield(ParYield, Annuity, crv, StartDate, IndexMat, DayCount, IndexF);
}


/*****  Par_Yield_From_Dates  ************************************************/
/**
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 */

 int    Par_Yield_From_Dates
                (double*        ParYield /** (O) Par forward yield          */
                ,double*        Annuity  /** (O) Annuity                    */
                ,long           SwapSt    /** (I) Underlying swap start      */
                ,long           SwapMat   /** (I) Underlying swap maturity   */
                ,char           DCC       /** (I) Underlying day count conv. */
                ,char           Freq      /** (I) Underlying frequency       */
                ,char           StubConv  /** (I) F, B or N(one allowed)     */
                ,const T_CURVE  *zc /** (I) Zero curve                */
		)
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

    if (SwapMat > zc->startDates[zc->numItems-1])
    {        
        DR_Error ("Not enough zeros to calculate par yield/annuity!");
        goto FREE_MEM_AND_RETURN;
    }

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
    ZerotoS = GetZeroPrice(SwapSt, zc); 
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
        ZerotoCpn = GetZeroPrice(CpnPmtDate, zc);
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

/* like one above but the right way */
int    ParYieldFromDates
                (double*                ParYield /** (O) Par forward yield          */
                ,double*                Annuity  /** (O) Annuity                    */
                ,long                   SwapSt   /** (I) Underlying swap start      */
                ,long                   SwapMat  /** (I) Underlying swap maturity   */
                ,char                   DCC      /** (I) Underlying day count conv. */
                ,char                   Freq     /** (I) Underlying frequency       */
                ,char                   StubConv /** (I) F, B or N(one allowed)     */
                ,T_CURVE const*   crv      /** zero curve */
		)
{
    return Par_Yield_From_Dates(
            ParYield,
            Annuity,
            SwapSt,
            SwapMat,
            DCC,
            Freq,
            StubConv,
            crv);
}


/*****  ParYieldRatio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */

int  ParYieldRatio(double*        Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                  ,long           StartDate1 /** (I) start date of yield 1      */
                  ,long           StartDate2 /** (I) start date of yield 2      */
                  ,double         Spread     /** (I) added to both yields       */
                  ,const T_CURVE  *zc  /** (I) Zero curve                 */
                  ,int            IndexMat   /** (I) Index maturity             */
                  ,char           DayCount   /** (I) Index Day count convention */
                  ,char           IndexF     /** (I) Index payment frequency    */
		)
{
    int     status = FAILURE;
    double  Yield1,   Yield2;
    double  Annuity1, Annuity2;
    long    MatDate1, MatDate2;

    MatDate1 = Nxtmth (StartDate1, IndexMat, 1L);
    MatDate2 = Nxtmth (StartDate2, IndexMat, 1L);

    if (Par_Yield_From_Dates
                 (&Yield1,
                  &Annuity1,
                  StartDate1,
                  MatDate1,
                  DayCount,
                  IndexF,
                  'N',
                  zc) == FAILURE) goto RETURN;

    if (Par_Yield_From_Dates
                 (&Yield2,
                  &Annuity2,
                  StartDate2,
                  MatDate2,
                  DayCount,
                  IndexF,
                  'N',
                  zc) == FAILURE) goto RETURN;

    *Ratio = (Yield1 + Spread)/(Yield2 + Spread);

    status = SUCCESS;

RETURN:

    if (status == FAILURE) DR_Error("ParYieldRatio: Failed!");

    return(status);

}/* ParYieldRatio */

/*****  Par_Yield_Ratio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */
int  Par_Yield_Ratio(double*        Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                    ,long           StartDate1 /** (I) start date of yield 1      */
                    ,long           StartDate2 /** (I) start date of yield 2      */
                    ,double         Spread     /** (I) added to both yields       */
                    ,const T_CURVE  *zc  /** (I) Zero curve                 */
                    ,int            IndexMat   /** (I) Index maturity             */
                    ,char           DayCount   /** (I) Index Day count convention */
                    ,char           IndexF     /** (I) Index payment frequency    */
		)
{
        return (ParYieldRatio(Ratio,
                              StartDate1, 
                              StartDate2,
                              Spread,   
                              zc,
                              IndexMat, 
                              DayCount,
                              IndexF));

}

/* frieldly version */
int  ParYldRatio(double *       ratio       /** (O) Yield1+sprd / Yield2+sprd  */
                ,ESL_DATE       startDate1  /** (I) start date of yield 1      */
                ,ESL_DATE       startDate2  /** (I) start date of yield 2      */
                ,double         spread      /** (I) added to both yields       */
                ,const T_CURVE *zc    /** (I) index zero curve           */
                ,size_t         indexMat    /** (I) Index maturity             */
                ,ESL_DCC        dcc         /** (I) Index Day count convention */
                ,ESL_FREQ       freq        /** (I) Index payment frequency    */
		)
{
        return ParYieldRatio(ratio,
                             startDate1, 
                             startDate2,
                             spread,   
                             zc,
                             indexMat, 
                             (char)dcc,
                             (char)freq);
}


