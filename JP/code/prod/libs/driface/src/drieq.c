/************************************************************************
 * Module:      driface
 * File:        drieq.c
 * Function:    DR wrapper interface to read in equity.sta 
 * 		and equity.dyn files
 * Author:      David Liu Sept. 1998
 ************************************************************************/
#include "drlstd.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "cgeneral.h"
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
#include "bivari.h"		/* GtoBiVariOption */

#include "drieq.h"

#include "check.h"
#include "swapadj.h"

#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h"
#include "drlsmat.h"
#include "drltime.h"
#include "dritkwrp.h"		/* TDrWrapperData routines */

#define xTRI_DEBUG

#define DRI_INTERP_TYPE GTO_LINEAR_INTERP

#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}

/*f---------------------------------------------------------------------
 * DR wrapper to build equity static oject TEqStatData from equity.sta file
 * Returns SUCCESS/FAILURE.
 * Since dividend pay dates are not given in the equity.sta file, following
 * assumptions are made:  
 * 1) for rolling settlement, dividend is settled after same number 
 *    of settlement days as trade.
 * 2) for fixed settlement, pay date = ex-dividend date.
 */

GTO_EXPORT(int)
DriTEqStatGet(
char 	     *pathdir,		/* (I) tmp direcory name (or NULL) */
char 	     *eqStaFnam,	/* (I) file name (equity.sta)  */
long	     busDayConv,        /* (I) Business Day Conv    */
char	     *holidayFile,      /* (I) Holiday File       */
TEqStatData  **eqStatData)      /* (O) equity static data */
{
    static	char	routine[] = "DriTEqStatGet";
    int	status = FAILURE;

    char   fnam[256];
    FILE   *fp = NULL;

    long   divPaySettleDays = 0;  

    char   buf[256];
    char   *label = "Dividend";

    long   numDivs;
    TDate  *exDivDates = NULL;
    TDate  *payDivDates = NULL;
    double *divAmount = NULL;
    long   *divTypes = NULL;
    char   *divStrs = NULL;

    long     	numSettle, settleDays;
    char 	settleType;
    TDateList   *lastTradingDates = NULL;
    TDateList   *settlementDates = NULL;
  
    TDate       stockDate = 0;
    TDateList   *stmPeriodDates = NULL;
    long        numStmPeriods = 1;
    long	stmPeriod[1];

    long i;
    
    TDividendList     *divList = NULL;
    TEquitySettlement *stm = NULL;


    /* Open equity.sta 
     */
    if (pathdir != NULL) {
	sprintf(fnam, "%s/%s", pathdir, eqStaFnam);
    } else {
	strcpy(fnam, eqStaFnam);
    }

    if ((fp = fopen(fnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, fnam, strerror(errno));
	goto done;
    }
    

    /* 
     * Skip all before divident line start
     */
    if(DrlFAdvanceToToken(fp,
			  label) == FAILURE)
    {
	GtoErrMsg("%s: Cannot find string token %s.\n", routine, label);
	goto done;
    }

    /* 
     * The file pointer stops after the token string.  
     * Need to skip rest of the line 
     */
    if(fgets(buf, sizeof(buf), fp) == NULL) 
    {
	GtoErrMsg("%s: EOF encountered before divident list starts.\n",
                     routine);
	goto done;
    }

    /* 
     * Read the divident list
     */
    READ_DATA(DRL_LONG_T, &numDivs,    "number of dividents");

    if (numDivs <= 0)     /* No dividend. Create a dummy one */ 
    {
	numDivs = 1;
	if((exDivDates = NEW_ARRAY(TDate, numDivs))  == NULL ||
	   (divAmount  = NEW_ARRAY(double, numDivs)) == NULL ||
	   (divStrs    = NEW_ARRAY(char, numDivs))   == NULL)
		goto done;

	exDivDates[0] = 100;
	divAmount[0]  = 0.0;
	divStrs[0]    = 'D';
    }
    else
    {
    	if(DrlLilVectArrayFpReadV(fp, 
			      	  numDivs,
			      	  DRL_TDATE_T,   (void*) &exDivDates,
			      	  DRL_DOUBLE_T,  (void*) &divAmount,
			      	  DRL_CHAR_T,    (void*) &divStrs,
			      	  DRL_NULL_T) == FAILURE)
    	{  
        	GtoErrMsg("%s: Cannot read divident lists.\n", routine);
        	goto done;
    	}
    }

    /*
     * Read settlement date
     */
    READ_DATA(DRL_LONG_T, &numSettle,    "number of settlement points");

    if (numSettle <= 0)    /* Rolling settlement */
    {
	settleType = 'R';
	READ_DATA(DRL_LONG_T, &settleDays,    "settlement days");

    	stmPeriodDates = GtoNewEmptyDateList(0);
    	stmPeriod[0] = settleDays;

	divPaySettleDays = settleDays;
    }
    else	 /* Fixed settlement */
    {
	settleType = 'F';

	if ((lastTradingDates = GtoNewEmptyDateList(0)) == NULL ||
	    (settlementDates  = GtoNewEmptyDateList(0)) == NULL)
		goto done;

	lastTradingDates->fNumItems = numSettle;
	settlementDates->fNumItems = numSettle;

	if(DrlLilVectArrayFpReadV(fp,
                               numSettle,
                               DRL_TDATE_T, (void*) &(lastTradingDates->fArray),
                               DRL_TDATE_T, (void*) &(settlementDates->fArray),
                               DRL_NULL_T) == FAILURE)
	{
	    GtoErrMsg("%s: Cannot read fixed settlement lists.\n", routine);
            goto done;
	}
    }


    /* 
     * create divident pay dates from ex-divident dates and settlement dates
     */
    if((payDivDates = NEW_ARRAY(TDate, numDivs)) == NULL ||
       (divTypes    = NEW_ARRAY(long, numDivs)) == NULL)
    {
	GtoErrMsg("%s: Error allocating memory for payDivDates or divTypes.\n",
		  routine);
	goto done;
    }

    /* 
     * Generate pay dates and dividend types 
     * For continuous dividend, it requires ex-div dates 
     * and payment dates are equal
     */
    for (i=0; i<=numDivs-1; i++){

	switch (toupper(divStrs[i]))
	{
	case 'D':
		divTypes[i] = GtoDIVIDENDRATE_AMOUNT;
		if(GtoDateFromBusDaysOffset(exDivDates[i],
                                    	    divPaySettleDays,
                                    	    holidayFile,
                                    	    &payDivDates[i]) == FAILURE)
	    		goto done;
		break;
	case 'Y':
		divTypes[i] = GtoDIVIDENDRATE_PERCENT;
		if(GtoDateFromBusDaysOffset(exDivDates[i],
                                    	    divPaySettleDays,
                                    	    holidayFile,
                                    	    &payDivDates[i]) == FAILURE)
	    		goto done;
		break;
	case 'C':
		divTypes[i] = GtoDIVIDENDRATE_CONTINUOUS;
		payDivDates[i] = exDivDates[i];
		break;
	default:
		GtoErrMsg("%s: Unknown divident type %c for the %dth divident."
			  "\n",
                  	  routine, divStrs[i], i+1);
        	goto done;
	}
    }


    /* Construct TDividendList and TEquitySettment
     */
    if((divList = GtoNewTDividendList(numDivs,
			              exDivDates,
			              payDivDates,
				      divTypes,
				      divAmount)) == NULL ||
       (stm = GtoEqSettleNew( stockDate,
			      stmPeriodDates,
    			      numStmPeriods,
    			      stmPeriod, 
			      lastTradingDates,
			      settlementDates,
			      holidayFile,
			      busDayConv,
			      holidayFile,
			      busDayConv,
			      GTO_EQ_STM_DIRN_INC)) == NULL)	
	goto done;
			 

    /* 
     * Construct eqStatData
     */

    if((*eqStatData = DrNewTEqStatData(divList,
				      settleType,
				      stm)) == NULL) 
	goto done;


    status = SUCCESS;

  done:

    FREE(exDivDates);
    FREE(payDivDates);
    FREE(divAmount);
    FREE(divStrs);
    FREE(divTypes);
    GtoFreeDateList(stmPeriodDates);
    GtoFreeDateList(lastTradingDates);
    GtoFreeDateList(settlementDates);
    
    if (fp != NULL) fclose(fp);

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);

    return(status);
}


GTO_EXPORT(extern TEqStatData *)
DrNewTEqStatData(TDividendList     *divList,
		 char		   settleType,
                 TEquitySettlement *stm)
{
    static      char    routine[] = "DrNewTEqStatData";
    int status = FAILURE;

    TEqStatData *eqStatData = NULL;
    if ((eqStatData = NEW(TEqStatData)) == NULL)
	goto done;

    eqStatData->divList = divList;
    eqStatData->settleType = settleType;
    eqStatData->stm = stm;

    status = SUCCESS;

done:

    if (status == FAILURE)
    {
        DrFreeTEqStatData(eqStatData);
        GtoErrMsg("%s: Failed\n", routine);
        return NULL;
    }

    return eqStatData;
}
  

GTO_EXPORT(extern void)
DrFreeTEqStatData(TEqStatData *eqStatData)
{
    if (eqStatData == NULL)
	return;

    GtoFreeTDividendList(eqStatData->divList);
    GtoEqSettleFree(eqStatData->stm);

    FREE(eqStatData);
}



/*f---------------------------------------------------------------------
 * Read data from equity.dyn in London format
 * Basis of vol curve is set to 4L (quarterly)
 *
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriEqDynDataGet(   
char    *pathdir,	       	/* (I) tmp direcory name (or NULL) */
char 	*eqDynFnam,		/* (I) file name (equity.dyn)  */
double  *spotValue,            	/* (O)  */
double  *correlation,          	/* (O)  */
TCurve **indxVolCurve)         	/* (O) index volatility curve */
{
    static	char	routine[] = "DriEqDynDataGet";
    int	status = FAILURE;

    char	fnam[256];
    FILE *fp = NULL;

    TDate baseDate;

    long numVolPts;
    long basis = 4L;
    long dcc = GTO_ACT_365F;

    TCurve *vc = *indxVolCurve = NULL;

    long i;

    /* Open equity.sta 
     */
    if (pathdir != NULL) {
	sprintf(fnam, "%s/%s", pathdir, eqDynFnam);
    } else {
	strcpy(fnam, eqDynFnam);
    }

    if ((fp = fopen(fnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, fnam, strerror(errno));
	goto done;
    }
    
    /* Read constants
     */
    READ_DATA(DRL_TDATE_T,  &baseDate,	        "base date");
    READ_DATA(DRL_DOUBLE_T, spotValue,	        "spot value");
    READ_DATA(DRL_DOUBLE_T, correlation,	"correlation");
    READ_DATA(DRL_LONG_T,   &numVolPts,	        "num vol points");

    /* Build vol curve
     */
    if((vc = GtoNewTCurve(baseDate, numVolPts, basis, dcc)) == NULL)
	goto done;

    for(i=0; i<numVolPts; i++)
    {
	if(DrlFScanVType(fp,DRL_TDATE_T, 
			 (void*) &(vc->fArray[i].fDate)) == FAILURE ||
	   DrlFScanVType(fp,DRL_PERCENT_T, 
			 (void*) &(vc->fArray[i].fRate)) == FAILURE)
	{  
	    GtoErrMsg("%s: Cannot read volatility curve data.\n", routine);
	    goto done;
	}
    }
    
    *indxVolCurve = vc;
    
    status = SUCCESS;
    
  done:
    
    if (fp != NULL) fclose(fp);

    if (status != SUCCESS)
    {
	GtoFreeTCurve(vc);
	GtoErrMsg("%s: failed.\n", routine);
    }

    return(status);
}




/*f---------------------------------------------------------------------
 * Given a swap zero curve, build a funding curve for a given index.
 * - read funding spreads from equity DR wrapper environment
 * - recover MM and par swap rates from swap curve
 * - add spreads to swap rates, and generate the corresponding zero curve
 *
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriTFundingCurveGet(   
char *pathdir,			/* (I) tmp direcory name (or NULL) */
TDrWrapperData *drWrapper,	/* (I) dr wrapper data */
TCurve        **indxZC)         /* (O) index zero curve with funding spread */
{
    static	char	routine[] = "DriTFundingCurveGet";
    int	status = FAILURE;

    static	char	eqStaFnam[] = "equity.sta";
    char	fnam[256];
    FILE *fp = NULL;

    TCurve *zc = drWrapper->fZcCurve; 
    TDate valueDate = zc->fBaseDate;
    int mmDenom = (int)(drWrapper->fMMDenom);
    long mmDayCount = mmDenom==365 ? GTO_ACT_365F : GTO_ACT_360;
    int swapFreq = (int)(drWrapper->fCmsSwMat->swapPayFreq);

    TDateInterval payInterval;
    long swapDcc = (long)(drWrapper->fSwDcc);

    long numRates = 0;
    double *rates = NULL;
    double *prices = NULL;
    TDate *dates = NULL;
    char *names = NULL;
    TDateInterval *rateMats = NULL;


    double temp;

    long i;
    
    *indxZC = NULL;

    /* Open equity.sta 
     */
    if (pathdir != NULL) {
	sprintf(fnam, "%s/%s", pathdir, eqStaFnam);
    } else {
	strcpy(fnam, eqStaFnam);
    }

    if ((fp = fopen(fnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, fnam, strerror(errno));
	goto done;
    }
    
    if (DrlFScanVType(fp, DRL_LONG_T, &numRates) != SUCCESS)
    {
	GtoErrMsg("%s: can't read number of benchmarks.\n", routine);
	goto done;
    }

    if((dates=NEW_ARRAY(TDate, numRates)) == NULL ||
       (prices=NEW_ARRAY(double, numRates)) == NULL ||
       (names=NEW_ARRAY(char, numRates)) == NULL)
    {
	GtoErrMsg("%s: Error allocating memory.\n", routine);
	goto done;
    }

     if(DrlLilVectArrayFpReadV(fp, 
			       numRates,
			       DRL_TDATEINTERVAL_T,  (void*) &rateMats,
			       DRL_PERCENT_T,  (void*) &rates,
			       DRL_NULL_T) == FAILURE)
    {  
        GtoErrMsg("%s: Cannot read reset dates.\n", routine);
        goto done;
    }

    /* Fill in arrays for zc generating and compute the funding
     * rates with spreads. Rates <=1y are MM, >1y par rates
     */
    for(i=0; i<numRates; i++)
    {
	if( GtoDtFwdAny(valueDate,
			rateMats+i,
			dates+i) == FAILURE)
	    goto done;

	if(dates[i]-valueDate < 367) /* less than a year = MM */
	{
	    names[i] = 'M';
	    if(GtoZerosToSimplePoint(zc,
				     DRI_INTERP_TYPE,
				     valueDate,
				     dates[i],
				     mmDayCount,
				     &temp) == FAILURE)
	    goto done;
	}
	else
	{
	    names[i] = 'S';
	    if(GtoFreq2TDateInterval(swapFreq,
				     &payInterval) == FAILURE ||
	       GtoZerosToCouponsPoint(zc,
				      DRI_INTERP_TYPE,
				      valueDate,
				      &payInterval,
				      dates[i],
				      swapDcc,
				      GTO_STUB_NONE,
				      FALSE,
				      &temp) == FAILURE)
	    goto done;
	}

	rates[i] += temp;
	prices[i] = 1.;
    }
   
    /* Generate funding curve
     */
    if(((*indxZC) = GtoNPiZC(valueDate,
			     mmDenom,
			     swapFreq,
			     swapDcc,
			     rates,
			     dates,
			     prices,
			     names,
			     GtoSETINSTR,
			     numRates,
			     NULL,
			     NULL,
			     0, /* num futures */
			     0, /* num FRAs */
			     NULL,
			     0.,
			     0,
			     NULL,
			     0,
			     GtoLINEARINTERP, /* for coupons */
			     DRI_INTERP_TYPE, /* for zero rates */
			     0,
			     0)) /* no bus day adjust */
       == NULL) goto done;
			 
    status = SUCCESS;

  done:

    FREE(rates);
    FREE(dates);
    FREE(prices);
    FREE(rateMats);
    FREE(names);
    
    if (fp != NULL) fclose(fp);

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);

    return(status);
}

