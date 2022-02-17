/************************************************************************
 * Module:      driface
 * File:        bsbondopt.c
 * Function:    Simple vanilla option pricer
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
#define NUMIDX  1       /* 2 indices  */

#define TINY 1e-4

     
PRIVATE int driBSBondOptionLogInputs(
FILE	       *fpTERM_PRN,	/* (I) TERM.PRN */
TDate          today,		/* (I) spot date */
TDate		   optExpiry,
TDate	       payDate,
TDrWrapperData *drWrap,     /* (I) DR Wrapper Data  */
int            optType,   	/* (I) Max, Min, Diff */
TEqStatData    **eqStatData,/* (I) data in equity.sta files	*/
char           *holidayFile,/* (I) "NONE" for weekends only
							 *     "No_Weekends" for no adjustments */
long	       busDayConv,	/* (I) Business day convention  */
double         strike,      /* (I) strike  */
char           *routine);


DLL_EXPORT(int)
DriOutBSBondOption( 
		TDate          today,		/* (I) spot date */
		TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
		int            optType,   	/* (I) Max, Min, Diff */
		TDate		   optExpiry,
		TDate	       payDate,
		char		   distType,
		double		   spot,
		TCurve	       *idxVC,		/* (I) vol curves  */
		TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
		double         strike,
		double	       *pv)
{
    static char	routine[] = "DriOutBSBondOption";
    int	   status = FAILURE; 

    double *idxVol = NULL; 
    double *idxReturn = NULL; 
	double fwd1, pvFactor;
	
	double vol;

    double tExp;
	char *holidayFile = NULL;
    TDate  fwdDate[2];


	TCurve *repo1ZC;
	
	double fwdPrices[2];

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
     * Funding asset 1:	disczero.dat (curve 2)
     */
    repo1ZC = drWrap->fDiscZcCurve;

    /* 
     * Calculate the forward prices 
     */

	if (DrForwardPriceGen(  
				spot,
			    eqStatData[1],
				repo1ZC,
				BUSDAYCONV,
			    HOLIDAYFILE,
				2,
				fwdDate,
				fwdPrices) == FAILURE)
		goto done;
		 
	fwd1 = fwdPrices[1];


	/* interpolate vol */

	/* index volatility */ 
	if (GtoInterpRate(optExpiry,
			   idxVC,
			   GTO_LINEAR_INTERP,
			   &vol) != SUCCESS)
	{
	     goto done;

    }	/* idx */



	GTO_IF_LOGGING({
		GtoErrMsg( "Fwd1 = %12.4f  : vol = %12.4f\n", fwd1, vol);
    }); 

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

	 if ('N' == distType) {

		 if ((optResultPtr = 
				GtoNormalOption(
					optType, 
					fwd1, 
					strike, 
					tExp, 
					0, 
					vol*fwd1, 
					0, 
					'P')) == NULL)
				goto done;
		 else
			 optResult = *optResultPtr;

	 } else {

		 if (fwd1 <= 0e0) {

			 GtoErrMsg("%s: forward is zero or negative", routine);
			 goto done;

		 } else {

			if(GtoOptionsAnalytics2(
					optType, 
					fwd1, 
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

	*pv = optResult.fPrice*pvFactor;

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
DriOutBSBondOptionW(char *dataFnam)
{
    static	char	routine[] = "DriOutBSBondOptionW";
    int	status = FAILURE;

    TDate       today;
    TCurve      **idxVolCurve = NULL;
    TEqStatData **eqStatData = NULL;

    char        *holidayFile = HOLIDAYFILE;
    long        busDayConv  = BUSDAYCONV;


    double      *currIdxPrice = NULL;
    long        numResetDates = 1;   /*  Only one observation date */

    
    TDate       *resetDates = NULL;
    TDate       *payDates = NULL;


    char optType;
    double notional;
	TDate expDate, payDate;

    double      strike;
	double correlation1;

	char distType;
	
    FILE		*fp = NULL;
    FILE 	    *fpTERM_PRN = NULL; 

    static	char	defDataFnam[] = "bsbondopt_w.dat";
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
    READ_DATA(DRL_DOUBLE_T, &strike,      	"spread strike");
    READ_DATA(DRL_CHAR_T,    &distType,	        "distribution type");
    /* check option type */


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
        if (DriTDrWrapperDataPutPrice(0) != SUCCESS)
		goto done;

        return;
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
     * Discount curve:  zero.dat 
     * Funding curve 1:	disczero.dat (curve 2)
     */
    if (DriTDrWrapperDataGetFull(NULL, 
				 DRI_DRW_TYPE2_3CURVES,
				 &drWrap) != SUCCESS ||
	DriEqDynDataGet( NULL, "equity.dyn", 
			&currIdxPrice[1],
			&correlation1, 
			&idxVolCurve[1]) != SUCCESS ||
	DriTEqStatGet(NULL, "equity.sta", 
			busDayConv, 
			holidayFile,
                        &eqStatData[1]))
	goto done;
    

    /* Open TERM.prn
     */
    fpTERM_PRN = fopen("TERM.PRN", "w");
    if (fpTERM_PRN IS NULL)
    {
        GtoErrMsg("%s: Cannot open TERM_PRN.\n",routine);
        goto done;
    }


    GtoLoggingSet(1); 

    /* Log inputs in fpTERM_PRN
     */
    GTO_IF_LOGGING({
		driBSBondOptionLogInputs(
			fpTERM_PRN,
			today, 
			expDate,
			payDate,
			drWrap, 
			optType, 
			eqStatData, 
			HOLIDAYFILE, 
			BUSDAYCONV, 
			strike, 
			routine);    
	});


    /* Call pricing routine 
     */
    if (DriOutBSBondOption(today,	   
		  	       drWrap,
			       optType,
				   expDate,
				   payDate,
				   distType,
				   currIdxPrice[1],
				   idxVolCurve[1],
				   eqStatData,
				   strike,
				   &pv) == FAILURE)
	goto done;

    pv *= notional/100.;

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


    FREE(currIdxPrice);
    


    if (idxVolCurve != NULL)
    {
        for (i=1;i<=NUMIDX;i++) GtoFreeTCurve(idxVolCurve[i]);
        FREE(idxVolCurve);
    }

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

PRIVATE int driBSBondOptionLogInputs(
FILE	       *fpTERM_PRN,	/* (I) TERM.PRN */
TDate          today,		/* (I) spot date */
TDate		   optExpiry,
TDate	       payDate,
TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
int            optType,   	/* (I) Max, Min, Diff */
TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
char           *holidayFile,    /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
long	       busDayConv,	/* (I) Business day convention  */
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
     */
    GtoPrintTCurve(drWrap->fZcCurve, "Discount Curve");
    GtoErrMsg("\n");

    GtoPrintTCurve(drWrap->fDiscZcCurve, "1st Index Curve");
    GtoErrMsg("\n");


    for (idx=1; idx<=1; idx++)
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


    GtoErrMsg("Strike:                   %f\n", strike);

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
