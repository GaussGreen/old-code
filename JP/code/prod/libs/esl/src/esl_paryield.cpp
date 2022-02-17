/****************************************************************************/
/*      Calculation of par yield from zero bank.                            */
/****************************************************************************/
/*      ESL_PARYIELD.c                                                      */
/****************************************************************************/
/****************************************************************************
** IMPORTANT: This file contains conditionally compiled code depending upon
**            whether one is using IRX curves (i.e.
**            -DESL_NEW_DATE -DESL_NEW_CURVE) - IRX dates are a requisite for
**            IRX curves.
**
**            The code is organized so that the conditionally compiled code is
**            placed near the top of the file, and common code (without any
**            conditional compilation) is placed near the end.
**
**            PLEASE respect this convention when modifying this file.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esl_paryield.h"
#include "esl_zeros.h"
#include "esl_alloc.h"
#include "esl_util.h"
#include "esl_error.h"

#ifndef ESL_NEW_CURVE
/* THE FOLLOWING FUNCTIONS IN THIS CONDITIONALLY COMPILED BLOCK ARE DEPRECATED!!! */

/*****  Par_Yield_Plus ******************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
/* DEPRECATED!!! Use ParYieldPlus instead.                                  */
int Par_Yield_Plus(
                   double*      ParYield      /** (O) Par forward yield          */
                  ,double*      Annuity       /** (O) Annuity                    */
                  ,int          NbZero        /** (I) Number of zeros            */
                  ,double const*Zero          /** (I) Zero rates                 */
                  ,long   const*ZeroDate      /** (I) Zero maturity dates        */
                  ,long         CurrentDate   /** (I) Current date               */
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
        DFstart = ZeroPrice(StartDate,		/* (I) Discount bond maturity     */
			    CurrentDate,   	/* (I) Value Date                 */
			    NbZero,		/* (I) Number of zeros            */
			    ZeroDate,           /* (I) Zero maturity dates        */
			    Zero);   	   	/* (I) Zero rates                 */
        if (DFstart < 0.0) goto FREE_MEM_AND_RETURN;

        /* Price the zero coupon bond to EndDate */
        DFend   = ZeroPrice(EndDate,		/* (I) Discount bond maturity     */
			    CurrentDate,   	/* (I) Value Date                 */
		    	    NbZero,		/* (I) Number of zeros            */
	    		    ZeroDate,           /* (I) Zero maturity dates        */
    			    Zero);   	   	/* (I) Zero rates                 */
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

    return (status);
}  /* Par_Yield_Plus */


/*****  Par_Yield  **********************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
/* DEPRECATED!!! Use ParYield instead.                                      */
int    Par_Yield (double*       ParYield   /** (O) Par forward yield         */
                 ,double*       Annuity    /** (O) Annuity                   */
                 ,int           NbZero      /** (I) Number of zeros           */
                 ,double const* Zero       /** (I) Zero rates                */
                 ,long   const* ZeroDate   /** (I) Zero maturity dates       */
                 ,long          CurrentDate /** (I) Current date              */
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
			NbZero,				/* (I) Number of zeros            */
			Zero,				/* (I) Zero rates                 */
			ZeroDate,			/* (I) Zero maturity dates        */
			CurrentDate,		        /* (I) Current date               */
			StartDate,			/* (I) Forward start date         */
			'S',				/* (I) Rate type (cash or swap)	  */
			IndexMat,			/* (I) Index maturity             */
			DayCount,			/* (I) Index Day count convention */
			IndexF);			/* (I) Index payment frequency    */

}  /* Par_Yield */


/*****  Par_Yield_Ratio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */
/* DEPRECATED!!! Use ParYieldRatio instead.                                 */
int  Par_Yield_Ratio(double*        Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                  ,long           StartDate1 /** (I) start date of yield 1      */
                  ,long           StartDate2 /** (I) start date of yield 2      */
                  ,double         Spread     /** (I) added to both yields       */
                  ,int            NbZero     /** (I) Number of zeros            */
                  ,long    const* ZeroDates /** (I) Zero maturity dates        */
                  ,double  const* ZeroRates /** (I) Zero rates                 */
                  ,long           ValueDate  /** (I) Value date of zero curve   */
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
                  NbZero,
                  ZeroRates,
                  ZeroDates,
                  ValueDate) == FAILURE) goto RETURN;

    if (Par_Yield_From_Dates
                 (&Yield2,
                  &Annuity2,
                  StartDate2,
                  MatDate2,
                  DayCount,
                  IndexF,
                  'N',
                  NbZero,
                  ZeroRates,
                  ZeroDates,
                  ValueDate) == FAILURE) goto RETURN;

    *Ratio = (Yield1 + Spread)/(Yield2 + Spread);

    status = SUCCESS;

RETURN:

    if (status == FAILURE) DR_Error("Par_Yield_Ratio: Failed!");

    return(status);

}/* Par_Yield_Ratio */


/*****  Par_Yield_From_Dates  ************************************************/
/**
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 */
/* DEPRECATED!!!  USE ParYieldFromDates, instead!                           */
int    Par_Yield_From_Dates
                (double*        ParYield /** (O) Par forward yield          */
                ,double*        Annuity  /** (O) Annuity                    */
                ,long           SwapSt    /** (I) Underlying swap start      */
                ,long           SwapMat   /** (I) Underlying swap maturity   */
                ,char           DCC       /** (I) Underlying day count conv. */
                ,char           Freq      /** (I) Underlying frequency       */
                ,char           StubConv  /** (I) F, B or N(one allowed)     */
                ,int            NbZero    /** (I) Number of zeros            */
                ,double  const* Zero     /** (I) Zero rates                 */
                ,long    const* ZeroDate /** (I) Zero maturity dates        */
                ,long           BaseDate  /** (I) Zero curve base date       */
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

    /*
    if (SwapMat > ZeroDate[NbZero-1])
    {        
        DR_Error ("Not enough zeros to calculate par yield/annuity!");
        goto FREE_MEM_AND_RETURN;
    }
    */

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
#endif  /* ifndef ESL_NEW_CURVE */


/*****  ParYieldFromDates  ************************************************/
/**
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 *
 */
int ParYieldFromDates
                (double*        ParYield /** (O) Par forward yield          */
                ,double*        Annuity  /** (O) Annuity                    */
                ,IRDate          SwapSt    /** (I) Underlying swap start      */
                ,IRDate          SwapMat   /** (I) Underlying swap maturity   */
                ,char           DCC       /** (I) Underlying day count conv. */
                ,char           Freq      /** (I) Underlying frequency       */
                ,char           StubConv  /** (I) F, B or N(one allowed)     */
                ,const T_CURVE  *zc /** (I) Zero curve                */
		        )
{
    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmts */
    IRDate        CpnPmtDate;           /* Current cpn pmt date    */
    IRDate        TempDates[2];         /* To construct cpn list   */

    double  DCCFrac;             /* Current coupon day count fraction  */
    double  AnnuityL=0.0;        /* Local variable for annuity price   */
    double  ParYieldL;           /* Local variable for forward yield   */

    double  ZerotoS;             /* Zero to swap start                 */
    double  ZerotoCpn=0.0;       /* Zero to current cpn                */

    int     j;
    int     status = FAILURE;    /* Error status = FAILURE initially   */

    const IRDate*    zcStartDates;
    long            zcNumItems;

#ifdef ESL_NEW_CURVE
    zcStartDates = zc->startDates;
    zcNumItems = zc->numItems;
#else
    zcStartDates = zc->ZeroDate;
    zcNumItems = zc->NbZero;
#endif

    /*
    if (SwapMat > zcStartDates[zcNumItems-1])
    {        
        DR_Error ("Not enough zeros to calculate par yield/annuity!");
        goto FREE_MEM_AND_RETURN;
    }
    */

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
                      "(ParYieldFromDates)!");
            goto FREE_MEM_AND_RETURN;
        }

        /* Zero to current cpn */
        ZerotoCpn = GetZeroPrice(CpnPmtDate, zc);
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;
        
        AnnuityL += DCCFrac * ZerotoCpn;

    }  /* for j */

    if (AnnuityL < 0.0)
    {
        DR_Error("Annuity is < 0.0. (ParYieldFromDates)");
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

}  /* ParYieldFromDates */


/***** ParYieldPlus *******************************************************
**
**      Calculate par forward yield for a specific maturity from a zero 
**      coupon curve (deterministic). Does not allow for stub.
*****************************************************************************/
int ParYieldPlus(
                   double*      ParYield      /** (O) Par forward yield          */
                  ,double*      Annuity       /** (O) Annuity                    */
                  ,const T_CURVE  *zc         /** (I) Zero curve                 */
                  ,IRDate        StartDate     /** (I) Forward start date         */
                  ,char         IndexType     /** (I) Rate type (cash or swap)   */
                  ,int          IndexMat      /** (I) Index maturity             */
                  ,char         DayCount      /** (I) Index Day count convention */
                  ,char         IndexF        /** (I) Index payment frequency    */
		)
{
    double *DayCntFrn=NULL;     /* Accrual periods                         */
    double *ZeroToPmt=NULL;     /* Zero coupons maturing on index payments */
    IRDate  *PmtDate=NULL;       /* Index payment dates                     */
											
    IRDate AccStart;              /* Accrual start of current payment        */
    IRDate AccEnd;                /* Accrual end date of current payment     */
											
    int IndexFInt;              /* Index frequency as integer              */
    int NbReset=0;              /* Total number of resets in the index     */
    int j;                      /* Loop index                              */
    int status = FAILURE;       /* Error status = FAILURE initially        */

    double DFstart;             /* Discount bond to start date */
    double DFend;               /* Discount bond to end date */
    IRDate EndDate;               /* Rate end date */
    double DCF;                 /* Cash rate accrual */

    IRDate CurrentDate = GetZeroCurveBaseDate(zc);

	/* Date ordering */
    if (StartDate < CurrentDate)                            
    {        
        DR_Error("ParYieldPlus: start date falls before current date!");
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
		DR_Error("ParYieldPlus: Failed to compute cash accrual!");
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
		PmtDate   = (IRDate *)  DR_Array (IDATE,   0, NbReset);

		/* Check for memory allcoation failure */
		if (   (DayCntFrn == NULL)
			|| (ZeroToPmt == NULL)
			|| (PmtDate   == NULL))
		{
			DR_Error("ParYieldPlus: could not allocate memory!");
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
				DR_Error("ParYieldPlus: could not calculate day count fractions!");
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
		DR_Error("ParYieldPlus: Rate type must be 'C'ash or 'S'wap.");
		goto FREE_MEM_AND_RETURN;        
	}

    /* If we get this far, it worked */
    status = SUCCESS;

FREE_MEM_AND_RETURN:

	/* Release memory */
    Free_DR_Array (DayCntFrn, DOUBLE, 0, NbReset);
    Free_DR_Array (ZeroToPmt, DOUBLE, 0, NbReset);
    Free_DR_Array (PmtDate, IDATE, 0, NbReset);

    return (status);
}  /* ParYieldPlus */


/*****  ParYield  **********************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int    ParYield  (double*       ParYield   /** (O) Par forward yield         */
                 ,double*       Annuity    /** (O) Annuity                   */
                 ,const T_CURVE  *zc       /** (I) Zero curve                */
                 ,IRDate         StartDate  /** (I) Forward start date        */
                 ,int           IndexMat   /** (I) Index maturity            */
                 ,char          DayCount   /** (I) Index Day count convention*/
                 ,char          IndexF     /** (I) Index payment frequency   */
		)
{
	/* Use the new function in swap-rate mode */
	return ParYieldPlus(
			ParYield,			/* (O) Par forward yield          */
			Annuity,			/* (O) Annuity                    */
			zc,				    /* (I) Zero curve                 */
			StartDate,			/* (I) Forward start date         */
			'S',				/* (I) Rate type (cash or swap)	  */
			IndexMat,			/* (I) Index maturity             */
			DayCount,			/* (I) Index Day count convention */
			IndexF);			/* (I) Index payment frequency    */

}  /* ParYield */


/*****  ParYieldRatio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */
int  ParYieldRatio(double*        Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                  ,IRDate           StartDate1 /** (I) start date of yield 1      */
                  ,IRDate           StartDate2 /** (I) start date of yield 2      */
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
    IRDate    MatDate1, MatDate2;

    MatDate1 = Nxtmth (StartDate1, IndexMat, 1L);
    MatDate2 = Nxtmth (StartDate2, IndexMat, 1L);

    if (ParYieldFromDates
                 (&Yield1,
                  &Annuity1,
                  StartDate1,
                  MatDate1,
                  DayCount,
                  IndexF,
                  'N',
                  zc) == FAILURE) goto RETURN;

    if (ParYieldFromDates
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



/*****  Swap_Yield_From_Dates  **********************************************/
/*
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 *       Uses the formula sum(FwdLibor*Z)/A, which is the correct formula to
 *       use in the presence of ccy basis.
 *       Please refer to Par_Yield_From_Dates for comparison
 */

int  Swap_Yield_From_Dates
             (double     *SwapYield,    /* (O) Par forward yield            */
              double     *Annuity,      /* (O) Annuity                      */
              long        SwapSt,       /* (I) Underlying swap start        */
              long        SwapMat,      /* (I) Underlying swap maturity     */
              char        DCC,          /* (I) Underlying day count conv.   */
              char        Freq,         /* (I) Underlying frequency         */
              char        StubConv,     /* (I) F, B or N(one allowed)       */
              const T_CURVE    *IdxZCurve,
              const T_CURVE    *DiscZCurve)
{
    static char  routine[] = "Swap_Yield_From_Dates";
    int          status = FAILURE;
    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmts */
    long         TempDates[2];        /* To construct cpn list   */

    double *DCCFrac = NULL;
    double *FwdZero = NULL;    /* from CurrDate to Ti, in disc crv */
    double *Ratio   = NULL;

    int     j;

    /* temp var for speed */
    double  ZerotoCpn = 0.0;
    double  PrevZero;

    double  SumCoupon = 0.0;    /* sum(Libor*Zero) */
    double  AnnuityL  = 0.0;

    long    ValueDate;

    if ((IdxZCurve  == NULL) || 
        (DiscZCurve == NULL) ||
        (SwapYield  == NULL) ||
        (Annuity    == NULL)) goto FREE_MEM_AND_RETURN;

    ValueDate = IdxZCurve->Today;

    if (ValueDate > SwapSt)
    {
        DR_Error("Swap start is before base date."
                 "(%s)", routine);
        goto FREE_MEM_AND_RETURN;
    }

    /* create the coupon payment date list */
    TempDates[0] = SwapSt;
    TempDates[1] = SwapMat;
    CpnEventList = DrNewEventListFromFreq(2,
                                          TempDates,
                                          Freq,
                                          StubConv,
                                          'N',   /* Dates in not required */
                                          NULL, NULL, NULL, NULL, NULL);
    if (CpnEventList == NULL) goto FREE_MEM_AND_RETURN;

    if (CpnEventList->NbEntries < 2) goto FREE_MEM_AND_RETURN;

    /* allocate memory for temp var */
    FwdZero = (double *) DR_Array (DOUBLE, 0, CpnEventList->NbEntries-1);
    Ratio   = (double *) DR_Array (DOUBLE, 0, CpnEventList->NbEntries-2);
    DCCFrac = (double *) DR_Array (DOUBLE, 1, CpnEventList->NbEntries-1);

    if ((FwdZero == NULL) ||
        (Ratio   == NULL) ||
        (DCCFrac == NULL)) goto FREE_MEM_AND_RETURN;

    /* prepare other intermediate quantities */

    for (j = 0; j < CpnEventList->NbEntries; j++)
    {
        FwdZero[j] = GetZeroPrice(CpnEventList->Dates[j], DiscZCurve);
    }

    for (j = 0; j < CpnEventList->NbEntries-1; j++)
    {
        double  idxZero1, idxZero2;

        idxZero1 = GetZeroPrice(CpnEventList->Dates[j],IdxZCurve);

        idxZero2 = GetZeroPrice(CpnEventList->Dates[j+1],IdxZCurve);

        Ratio[j] = (idxZero1 / idxZero2) * (FwdZero[j+1] / FwdZero[j]);
    }

    for (j = 1; j < CpnEventList->NbEntries; j++)
    {
        if (DrDayCountFraction (CpnEventList->Dates[j-1],
                                CpnEventList->Dates[j],
                                DCC,
                                &(DCCFrac[j])) == FAILURE)
        {
            DR_Error ("Could not calculate day count fraction "
                      "(%s)!", routine);
            goto FREE_MEM_AND_RETURN;
        }
    }

    /* calculate the spot annuity and sum(Libor*Zero) */

    PrevZero = FwdZero[0];

    for (j=1; j<CpnEventList->NbEntries; j++)
    {
        ZerotoCpn = FwdZero[j];

        SumCoupon  += (PrevZero * Ratio[j-1] - ZerotoCpn);
        AnnuityL   += DCCFrac[j] * ZerotoCpn;

        PrevZero = ZerotoCpn;
    }

    /* populate the results in the output fields */ 

    *SwapYield = SumCoupon / AnnuityL;
    *Annuity   = AnnuityL;

    status = SUCCESS;


FREE_MEM_AND_RETURN:

    {
        long NbEntries = (CpnEventList==NULL) ? 1L : CpnEventList->NbEntries;

        Free_DR_Array (FwdZero, DOUBLE, 0, NbEntries-1);
        Free_DR_Array (Ratio,   DOUBLE, 0, NbEntries-2);
        Free_DR_Array (DCCFrac, DOUBLE, 1, NbEntries-1);
    }
    DrFreeEventList(CpnEventList);

    return (status);

}  /* Swap_Yield_From_Dates */

int  SwapYieldFromDates
             (double     *SwapYield,    /* (O) Par forward yield            */
              double     *Annuity,      /* (O) Annuity                      */
              IRDate        SwapSt,       /* (I) Underlying swap start        */
              IRDate        SwapMat,      /* (I) Underlying swap maturity     */
              char        DCC,          /* (I) Underlying day count conv.   */
              char        Freq,         /* (I) Underlying frequency         */
              char        StubConv,     /* (I) F, B or N(one allowed)       */
              const T_CURVE    *IdxZCurve, /* (I) index curve               */
              const T_CURVE    *DiscZCurve /* (I) discount curve            */
              )
{
    return Swap_Yield_From_Dates(SwapYield, Annuity, SwapSt, SwapMat, DCC, Freq, StubConv, IdxZCurve, DiscZCurve);
}
/************************************************************************************
**
** IMPORTANT:  READ COMMENTS AT TOP OF FILE BEFORE ADDING/MODIFYING CODE TO THIS FILE.
**
************************************************************************************/
