/************************************************************************
 * Module:      driface
 * File:        bondsprdopt.c
 * Function:    Spread-struck bond option pricer.
 * Author:      Juan C Pooras Aug. 2004
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "cgeneral.h"
#include "bondcnst.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "ldate.h"              /* GtoDayCountFraction */
#include "busday.h"   
#include "stub.h"               /* GTO_STUB_NONE */
#include "datelist.h"           /* GtoNewDateList */
#include "tcurve.h"             /* GtoDiscountDate */
#include "cerror.h"             /* GtoErrMsg */
#include "date_sup.h"
#include "matswapt.h"           /* GtoSwaptionMatrix2DPrint */
#include "convert.h"
#include "gtonpi.h"
#include "zr2coup.h"
#include "zr2simp.h"
#include "bondyld.h"
#include "normopt.h"		/* GtoBiVariOption */
#include "swappv.h"
#include "check.h"
#include "swapadj.h"
#include "presval.h"
#include "bondflow.h"
#include "yield.h"
#include "gtobf.h"
#include "drlio.h"
#include "drloptio.h" 
#include "drlstr.h"
#include "drlsmat.h"
#include "drltime.h"
#include "dritkwrp.h"		/* TDrWrapperData routines */

#include "eqfwd.h"		
#include "pfmopt.h"		
#include "drlstd.h"

#define xTRI_DEBUG

#define HOLIDAYFILE  "NONE"             /* currently no holiday file available
                                       * in Kapital environment */
#define BUSDAYCONV   GTO_BAD_DAY_MODIFIED
#define DRI_INTERP_TYPE GTO_FLAT_FORWARDS


#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
#undef 	NUMIDX
#define NUMIDX  2       /* 2 indices  */

#define TINY 1e-4

     
PRIVATE int driSpreadStrikeOptionLogInputs(
FILE	       *fpTERM_PRN,	/* (I) TERM.PRN */
TDate          today,		/* (I) spot date */
TDate		   optExpiry,
TDate	       payDate,
TDrWrapperData *drWrap,     /* (I) DR Wrapper Data  */
int            optType,   	/* (I) Max, Min, Diff */
double		   vol,         /* (I) index volatility curve */
TEqStatData    **eqStatData,/* (I) data in equity.sta files	*/
char           *holidayFile,/* (I) "NONE" for weekends only
							 *     "No_Weekends" for no adjustments */
long	       busDayConv,	/* (I) Business day convention  */
double		   coupon1,
long		   freq1,
double		   coupon2,
long		   freq2,
double         strike,      /* (I) strike  */
char           *routine);


DLL_EXPORT(int)
DriOutSpreadStrikeOption( 
		TDate          today,		/* (I) spot date */
		TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
		char           optType,   	/* (I) Max, Min, Diff */
		TDate		   optExpiry,
		TDate	       payDate,
		double		   vol,
		char		   distType,
		double		   dprice1,
		TDate	       matDate1,
		TDateInterval  freq1,
		double         coupon1,
		long	       dcc1,
		double         dprice2,
		TDate	       matDate2,
		TDateInterval  freq2,
		double         coupon2,
		long	       dcc2,
		TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
		double         strike,
		double	       *pv)
{
    static char	routine[] = "DriOutPerformanceIdxOption";
    int	   status = FAILURE; 

    double *idxVol = NULL; 
    double *idxReturn = NULL; 
	double fwd1, fwd2;
	double ytm1, ytm2;
	double nfreq1, nfreq2;
	long lfreq1, lfreq2;

    double pvFactor=0.0;
    double currOption = 0.0;
	double strikePrice;
    int i;
    TDate prevCoupon;
    double accCoupon1, accCoupon2, spot1, spot2;

    double tExp;

	char *holidayFile = NULL;
    TDate  fwdDate[2];


	TCurve *repo1ZC, *repo2ZC;

	double tMat1, tMat2;
	double fwdPrices[2];
	double yieldStrike;

	TOptionProperties optResult, *optResultPtr;

#ifdef TRI_DEBUG
    GtoLoggingSet(1); 
#endif

    /* 
     * Forward dates
     */

    fwdDate[0] = today;
	fwdDate[1] = optExpiry;

	fwdPrices[0] = 0;

    /* 
     * Discount curve:  zero.dat (curve 1)
     * Repo asset 1:	disczero.dat (curve 2)
     * Repo asset 2:	riskzero.dat (curve 3)
     */
    repo1ZC = drWrap->fDiscZcCurve;
    repo2ZC = drWrap->fRiskZcCurve;

    /* 
     * Calculate the forward prices 
     */

    if ((eqStatData[1]->divList->fNumItems < 2) || (eqStatData[2]->divList->fNumItems < 2))
    {
        GtoErrMsg("%s : must have at least two dividend dates per asset for coupon calc.\n");
        goto done;
    }

    i = 0;
    while (today > eqStatData[1]->divList->fArray[i+1].payDivDate)
    {
        i++;
        if (i ==  eqStatData[1]->divList->fNumItems)
        {
            GtoErrMsg("%s : must have dividend dates for all coupons from previous to expiry\n");
            goto done;
        };
    };

    prevCoupon = eqStatData[1]->divList->fArray[i].payDivDate;
    
    GtoDayCountFraction(prevCoupon, today, dcc1, &accCoupon1);
    accCoupon1 *= 100*coupon1;

    i = 0;
    while (today > eqStatData[2]->divList->fArray[i+1].payDivDate)
    {
        i++;
        if (i ==  eqStatData[2]->divList->fNumItems)
        {
            GtoErrMsg("%s : must have dividend dates for all coupons from previous to expiry\n");
            goto done;
        };
    };

    prevCoupon = eqStatData[2]->divList->fArray[i].payDivDate;

    GtoDayCountFraction(prevCoupon, today, dcc2, &accCoupon2);
    accCoupon2 *= 100*coupon2;

    spot1 = dprice1 + accCoupon1;
    spot2 = dprice2 + accCoupon2;

	if (DrForwardPriceGen(  
				spot1,
			    eqStatData[1],
				repo1ZC,
				BUSDAYCONV,
			    HOLIDAYFILE,
				2,
				fwdDate,
				fwdPrices) == FAILURE)
		goto done;
		 
	fwd1 = fwdPrices[1];

	if (DrForwardPriceGen(  
				spot2,
			    eqStatData[2],
				repo2ZC,
				BUSDAYCONV,
			    HOLIDAYFILE,
				2,
				fwdDate,
				fwdPrices) == FAILURE)
		goto done;

	fwd2 = fwdPrices[1];

	GtoDateIntervalToFreq(&freq1, &nfreq1);

	GtoDateIntervalToFreq(&freq2, &nfreq2);

	lfreq1 = (long) nfreq1;
	lfreq2 = (long) nfreq2;

	if (GtoDayCountFraction( optExpiry,
			     matDate1,
			     dcc1,
			     &tMat1) == FAILURE)
		goto done;

	if (GtoDayCountFraction( optExpiry,
			     matDate2,
			     dcc2,
			     &tMat2) == FAILURE)
		goto done;

	if (fwd1 < TINY) {
		GtoErrMsg("Error: Fwd1 =  %12.4f", fwd1);
		goto done;
	}

	if (fwd2 < TINY) {
		GtoErrMsg("Error: Fwd2 =  %12.4f", fwd2);
		goto done;
	}

	if (GtoBondYieldToMaturity(
			coupon1,
			fwd1/100.,
			lfreq1,
			tMat1,
			0.001,         /* yield guess */
			GTO_STUB_NONE,
			&ytm1) == FAILURE)
		goto done;

	if (GtoBondYieldToMaturity(
			coupon2,
			fwd2/100.,
			lfreq2,
			tMat2,
			0.001,
			GTO_STUB_NONE,
			&ytm2) == FAILURE)
		goto done;


	GTO_IF_LOGGING({
		GtoErrMsg( "Fwd1 = %12.4f  : yld1 = %12.6f, mat1 = %s\n", fwd1, ytm1, GtoFormatDate(matDate1));
		GtoErrMsg( "Fwd2 = %12.4f  : yld2 = %12.6f, mat2 = %s\n", fwd2, ytm2, GtoFormatDate(matDate2));
    });

	yieldStrike = ytm1+strike;

	if (fabs(yieldStrike-ytm2) < TINY)

		yieldStrike = ytm2 - TINY;

	GtoBondPVFromYield(coupon2, yieldStrike, lfreq2, tMat2, GTO_STUB_NONE, &strikePrice);

	/* Price option. 
     * Vol applies to the fwd spread, which is the difference of the two ytm's
     */

    if (GtoDayCountFraction( today,
			     optExpiry,
			     GTO_ACT_365F,
			     &tExp) == FAILURE)
	goto done;


     if(GtoDiscountDate( payDate,
                         drWrap->fZcCurve,
                         DRI_INTERP_TYPE,
                         &pvFactor) == FAILURE)
	goto done;


	 if ('C' == optType)
		 optType = 'P';
     else
		optType = 'C';
	 
	 if ('N' == distType) {

		 if ((optResultPtr = 
				GtoNormalOption(
					optType, 
					ytm2-ytm1, 
					strike, 
					tExp, 
					0, 
					vol*fabs(ytm2-ytm1), 
					0, 
					'P')) == NULL)
				goto done;
		 else
			 optResult = *optResultPtr;

	 } else {

		 if (ytm2-ytm1 <= 0e0) {

			 GtoErrMsg("%s: forward spread is zero or negative", routine);
			 goto done;

		 } else {

			if(GtoOptionsAnalytics2(
					optType, 
					ytm2-ytm1, 
					strike, 
					tExp, 
					0, 
					vol, 
					0, 
					0,
					0,
					'P', 
					&optResult) == FAILURE)
				goto done;
		 }
	 }

	*pv = optResult.fPrice  * pvFactor *(fwd2/100.-strikePrice)/(yieldStrike-ytm2);

    status = SUCCESS;
  
  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);

    return(status);
}



/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriTotalRetIdxBiOption}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "dritridx_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriOutSpreadStrikeptionW(char *dataFnam)
{
    static	char	routine[] = "DriOutPerformanceIdxOptionW";
    int	status = FAILURE;

    TDate       today;
    TCurve      **idxVolCurve = NULL;
    TEqStatData **eqStatData = NULL;

    char        *holidayFile = HOLIDAYFILE;
    long        busDayConv  = BUSDAYCONV;

    double      *lastIdxPrice = NULL;
    double      *currIdxPrice = NULL;
    long        numResetDates = 1;   /*  Only one observation date */
    double      *lastIdxPrice1 = NULL;
    double      *lastIdxPrice2 = NULL;
    TDate       *resetDates = NULL;
    TDate       *payDates = NULL;


    char optType;
    double notional, coupon1, coupon2;
	TDate expDate, payDate, mat1Date, mat2Date;
	long dcc1, dcc2;
	TDateInterval freq1, freq2;

    double      strike;
	double correlation1, correlation2;
	double nfreq1, nfreq2;
	long lfreq1, lfreq2;


    double      vol;
	char distType;
	
    FILE		*fp = NULL;
    FILE 	    *fpTERM_PRN = NULL; 

    static	char	defDataFnam[] = "sprdstrike_w.dat";
    static	char	todayFnam[] = "today.dat";
    TDrWrapperData	*drWrap = NULL;
    double		pv;
    
    long i;

    /* Read deal data 
     */

    if (dataFnam == NULL) dataFnam = defDataFnam;
    
    if ((fp = fopen(dataFnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, dataFnam, strerror(errno));
	goto done;
    }

    READ_DATA(DRL_CHAR_T,    &optType,	        "option type");
    READ_DATA(DRL_DOUBLE_T, &notional,	        "notional");
	READ_DATA(DRL_TDATE_T,    &expDate,	        "option expiry");
	READ_DATA(DRL_TDATE_T,    &payDate,	        "payment expiry");
	READ_DATA(DRL_TDATE_T,    &mat1Date,	        "maturity asset #1");
    READ_DATA(DRL_TDATEINTERVAL_T,    &freq1,	        "frequency asset #1");
    READ_DATA(DRL_DOUBLE_T, &coupon1,      	"coupon asset #1");
    READ_DATA(DRL_TDAYCOUNT_T,    &dcc1,	        "dcc asset #1");
	READ_DATA(DRL_TDATE_T,    &mat2Date,	        "maturity asset #2");
    READ_DATA(DRL_TDATEINTERVAL_T,    &freq2,	        "frequency asset #2");
    READ_DATA(DRL_DOUBLE_T, &coupon2,      	"coupon asset #2");
    READ_DATA(DRL_TDAYCOUNT_T,    &dcc2,	        "dcc asset #2");
    READ_DATA(DRL_DOUBLE_T, &strike,      	"spread strike");
    READ_DATA(DRL_CHAR_T,    &distType,	        "distribution type");
    /* check option type */

/*
    switch (optType)
    {
    case 'C':
	optType = GTO_OPTION_CALL;
	break;
    case 'P':
	optType = GTO_OPTION_PUT;
	break;

    default:
	GtoErrMsg("%s: Unknown option type (%d). Allowed values C/P.\n",
                  routine, optType);
        goto done;
    }
*/
    
    /* Close data file */
    if(fp) fclose(fp);


    /* Close data file and read today
     */
    if ((fp = fopen(todayFnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, todayFnam, strerror(errno));
	goto done;
    }
    READ_DATA(DRL_TDATE_T,  &today,	        "today");

    if (fp) fclose(fp);

    /* Check forward date starts after today
     */
    if(expDate < today)
    {	
	GtoErrMsg("%s: forward date (%s) <= today (%s).\n",
		  routine, GtoFormatDate(expDate), GtoFormatDate(today));
	goto done;
    }

    if(expDate > mat1Date)
    {	
	GtoErrMsg("%s: maturity date 1st asset (%s) <= expiry (%s).\n",
		  routine, GtoFormatDate(mat1Date), GtoFormatDate(expDate));
	goto done;
    }

	if(expDate > mat2Date)
    {	
	GtoErrMsg("%s: maturity date 2nd asset (%s) <= expiry (%s).\n",
		  routine, GtoFormatDate(mat2Date), GtoFormatDate(expDate));
	goto done;
    }

    idxVolCurve = NEW_ARRAY(TCurve *, NUMIDX+1);
    if (idxVolCurve == NULL)
        goto done;

    idxVolCurve[0] = (void *)NUMIDX;

    
    currIdxPrice = NEW_ARRAY(double, NUMIDX+1);
    if (currIdxPrice == NULL)
        goto done;

    currIdxPrice[0] = (double)NUMIDX;

    eqStatData = NEW_ARRAY(TEqStatData *, NUMIDX+1);
    if (eqStatData == NULL)
        goto done;
 
    eqStatData[0] = (void *)NUMIDX;


    /* 
     * Read market data. 
     * Discount curve:  zero.dat (curve 1, because it requires vol )
     * Funding curve 1:	disczero.dat (curve 2)
     * Funding curve 2: riskzero.dat (curve 3)
     */
    if (DriTDrWrapperDataGetFull(NULL, 
				 DRI_DRW_TYPE2_3CURVES,
				 &drWrap) != SUCCESS ||
	DriEqDynDataGet( NULL, "equity1.dyn", 
			&currIdxPrice[1],
			&correlation1, 
			&idxVolCurve[1]) != SUCCESS ||
	DriEqDynDataGet( NULL, "equity2.dyn", 
			&currIdxPrice[2],
			&correlation2, 
			&idxVolCurve[2]) != SUCCESS ||
	DriTEqStatGet(NULL, "equity1.sta", 
			busDayConv, 
			holidayFile,
                        &eqStatData[1]) != SUCCESS ||
	DriTEqStatGet(NULL, "equity2.sta", 
			busDayConv, 
			holidayFile,
                        &eqStatData[2]) != SUCCESS)
	goto done;
    

    /* Open TERM.prn
     */
    fpTERM_PRN = fopen("TERM.PRN", "w");
    if (fpTERM_PRN IS NULL)
    {
        GtoErrMsg("%s: Cannot open TERM_PRN.\n",routine);
        goto done;
    }

	GtoDateIntervalToFreq(&freq1, &nfreq1);

	GtoDateIntervalToFreq(&freq2, &nfreq2);

	lfreq1 = (long) nfreq1;
	lfreq2 = (long) nfreq2;

    
	if (GtoInterpRate(expDate,
			   idxVolCurve[2],
			   GTO_LINEAR_INTERP,
			   &vol) != SUCCESS)
	{
	     goto done;

    }	/* idx */

    GtoLoggingSet(1); 

    /* Log inputs in fpTERM_PRN
     */
    GTO_IF_LOGGING({
		driSpreadStrikeOptionLogInputs(
			fpTERM_PRN,
			today, 
			expDate,
			payDate,
			drWrap, 
			optType, 
			vol, 
			eqStatData, 
			HOLIDAYFILE, 
			BUSDAYCONV, 
			coupon1, 
			lfreq1, 
			coupon2, 
			lfreq2, 
			strike, 
			routine);    
	});



    /* Call pricing routine 
     */
    if (DriOutSpreadStrikeOption(today,	   
		  	       drWrap,
			       optType,
				   expDate,
				   payDate,
				   vol,
				   distType,
				   currIdxPrice[1],
				   mat1Date,
				   freq1,
				   coupon1/100.,
				   dcc1,
				   currIdxPrice[2],
				   mat2Date,
				   freq2,
				   coupon2/100.,
				   dcc2,
				   eqStatData,
				   strike/10000.,
				   &pv) == FAILURE)
	goto done;

    pv *= notional;

    GTO_IF_LOGGING({
	DrlFPrintf(fpTERM_PRN, "pv = %12.4f\n", pv);
    });

    /* */
    if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;
    printf("Price:     %f \n", pv);
   
	status = SUCCESS;
  done:
    if (fp) fclose(fp);
    if(fpTERM_PRN != NULL) fclose(fpTERM_PRN);

    DriTDrWrapperDataFree(drWrap);
    
    if(eqStatData != NULL)
    {
	for(i=1; i<=NUMIDX; i++)
	{
    		DrFreeTEqStatData(eqStatData[i]);
	}
    	FREE(eqStatData);
    }
    
    FREE(resetDates);
    FREE(payDates);
    FREE(lastIdxPrice1);
    FREE(lastIdxPrice2);
    FREE(currIdxPrice);
    FREE(lastIdxPrice);


    if (idxVolCurve != NULL)
    {
        for (i=1;i<=NUMIDX;i++) GtoFreeTCurve(idxVolCurve[i]);
        FREE(idxVolCurve);
    }

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

PRIVATE int driSpreadStrikeOptionLogInputs(
FILE	       *fpTERM_PRN,	/* (I) TERM.PRN */
TDate          today,		/* (I) spot date */
TDate		   optExpiry,
TDate	       payDate,
TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
int            optType,   	/* (I) Max, Min, Diff */
double		   vol,         /* (I) index volatility curve */
TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
char           *holidayFile,    /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
long	       busDayConv,	/* (I) Business day convention  */
double		   coupon1,
long		   freq1,
double		   coupon2,
long		   freq2,
double         strike,          /* (I) strike  */
char           *routine)
{
    int	status = FAILURE;
    int timeStamp = GtoErrMsgTimeStamp (0); /* disable temporarily */
    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN);  /* to switch back 
							  to error.log after 
							  logging is done */
    long i;

    int  idx;
  
    /* Dividend and settlement */
    char                settleType;
    char                divStr;
    TDividendList       *divList;
    TEqStmPrivate       *stm;
    long                settleDays;

    GtoErrMsg("\n%s INPUTS:\n", routine);
    GtoErrMsg("\n");
    GtoErrMsg("Today:           %s\n", GtoFormatDate(today));
    GtoErrMsg("\n");

    /* 
     * Discount curve:  zero.dat (curve 1, because it requires vol )
     * Funding curve 1:	disczero.dat (curve 2)
     * Funding curve 2: riskzero.dat (curve 3)
     */
    GtoPrintTCurve(drWrap->fZcCurve, "Discount Curve");
    GtoErrMsg("\n");

    GtoPrintTCurve(drWrap->fDiscZcCurve, "1st Index Curve");
    GtoErrMsg("\n");

    GtoPrintTCurve(drWrap->fRiskZcCurve, "2nd Index Curve");
    GtoErrMsg("\n\n");


    for (idx=1; idx<=2; idx++)
    {
  	GtoErrMsg("Bond%d:\n", idx);


	settleType = eqStatData[idx]->settleType;
	divList    = eqStatData[idx]->divList;
	stm        = (TEqStmPrivate *)eqStatData[idx]->stm;

    	/* Divident List */
    	GtoErrMsg("\n");
    	GtoErrMsg("# Coupon Dates\n\n");
    	GtoErrMsg("NUMBER_OF_POINTS: %d\n", divList->fNumItems);
    	GtoErrMsg("#  =====================================================\n");
        GtoErrMsg("#    exDiv Date	pay Date	Amount        Type\n");
    	GtoErrMsg("#  -----------------------------------------------------\n");
    	for (i=0; i<=divList->fNumItems-1; i++){
		switch (divList->fArray[i].divType)
		{
		case GtoDIVIDENDRATE_AMOUNT:
                	divStr = 'D';
                	break;
        	case GtoDIVIDENDRATE_PERCENT:
                	divStr = 'Y';
                	break;
        	case GtoDIVIDENDRATE_CONTINUOUS:
                	divStr = 'C';
                	break;
        	default:
                	GtoErrMsg("Unknown divident type %d for the %dth"
				  " divident.\n",
                          	  divList->fArray[i].divType, i+1);
                	goto done;
        	}

		GtoErrMsg("     %s\t%s\t%lf\t%c\n", 
		  	  GtoFormatDate(divList->fArray[i].exDivDate),
		  	  GtoFormatDate(divList->fArray[i].payDivDate),
		  	  divList->fArray[i].divAmount,
		  	  divStr);
    	}

    	GtoErrMsg("\n");
    	if (settleType == 'R'){
		settleDays = stm->stmPeriods[0];
		GtoErrMsg("Rolling Settlement: T+%ld\n", settleDays);
    	}
    	else
    	{
		settleDays = stm->lastTradingDates->fNumItems;
    		GtoErrMsg("Fixed Settlement: \n");
    		GtoErrMsg("NUMBER_OF_POINTS: %ld\n", settleDays);
    		GtoErrMsg("#  =============================================\n");
    		GtoErrMsg("#	Trade Date	Settle Date\n");
    		GtoErrMsg("#  ---------------------------------------------\n");
    		for (i=0; i<=settleDays-1; i++){
			GtoErrMsg(" \t%s\t%s\n", 
		  		GtoFormatDate(stm->lastTradingDates->fArray[i]),
		  		GtoFormatDate(stm->settlementDates->fArray[i]));
		}
    	}
	
    	GtoErrMsg("\n\n");

    }  /* idx */

    GtoErrMsg("Holiday file:             %s\n", holidayFile);


    GtoErrMsg("Coupon #1 :                   %f\n", coupon1);
	GtoErrMsg("Freq  #1 :                   %f\n", freq1);

    GtoErrMsg("Coupon #2 :                   %f\n", coupon2);
	GtoErrMsg("Freq  #2 :                   %f\n", freq2);

    GtoErrMsg("Strike:                   %f\n", strike);
    GtoErrMsg("Vol:                   %f\n", vol);

    GtoErrMsg("Expiry Date:         %s\n", GtoFormatDate(optExpiry));
    GtoErrMsg("Payment Date:             %s\n", GtoFormatDate(payDate));

    GtoErrMsg("\n");

    status = SUCCESS;

  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    GtoErrMsgTimeStamp (timeStamp);

    return(status);

}
