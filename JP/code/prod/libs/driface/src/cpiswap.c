/************************************************************************
 * Module:      driface
 * File:        cpiswap.c
 * Function:    CPI linked swap pricer
 * Author:      Julia Chislenko Sep 1999

$ Header:$
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
#include "tcurve.h"             /* GtoDiscountDate */
#include "cerror.h"             /* GtoErrMsg */
#include "date_sup.h"
#include "convert.h"
#include "zr2fwd.h"

#include "check.h"

#include "drlio.h"
#include "drlstr.h"
#include "drltime.h"
#include "dritkwrp.h"		/* TDrWrapperData routines */

#include "cpiswap.h"		/* Prototype consistency */

#define xCPI_DEBUG

#define CPI_NUM_FIXINGS 3

#define CPI_INTERP_TYPE GTO_LINEAR_INTERP

static FILE 	        *fpTERM_PRN = NULL; 

#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 0.000001 : (x)))

#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
    
PRIVATE int cpiSwapLogInputs(
TDate       today,              /* (I) spot date */
char        instrType,   	/* (I) (T)IP/(C)urrentPay/(A)ccretingCurrentPay */
double      origNotl,           /* (I) Original notional */
double      coupon,             /* (I) Coupon */
TDayCount   dayCountConv,       /* (I) For coupon payments */
double     *cpiFixings,         /* (I) [CPI_NUM_FIXINGS] Spot, last, first(for TIPs only) */
TCurve     *nomZC,              /* (I) nominal zero curve */
TCurve     *realZC,             /* (I) real zero curve */
TCurve     *discZC,             /* (I) discount zero curve */
long        numPayDates,        /* (I) Number of payment dates */
TDate      *accStartDates,      /* (I) [numPayDates] accrue start dates */
TDate      *accEndDates,        /* (I) [numPayDates] accrue end dates */
TDate      *payDates,           /* (I) [numPayDates] payment dates */
double      nomVol,             /* (I) nominal rate volatility */
double      realVol,            /* (I) real rate volatility */
double      nomMR,              /* (I) nominal (normal) mean reversion */
double      realMR,             /* (I) real (normal) mean reversion */
double      correlation,        /* (I) between nom and real interest rates */    
char       *routine);

/*f---------------------------------------------------------------------
 * Pricing routine for CPI linked swaps.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCPISwap(   
TDate       today,              /* (I) spot date */
char        instrType,   	/* (I) (T)IP/(C)urrentPay/(A)ccretingCurrentPay */
double      origNotl,           /* (I) original notional */
double      coupon,             /* (I) coupon */
TDayCount   dayCountConv,       /* (I) for coupon payments */
double     *cpiFixings,         /* (I) [CPI_NUM_FIXINGS] Spot, last, first(for TIPs only) */
TCurve     *nomZC,              /* (I) nominal zero curve */
TCurve     *realZC,             /* (I) real zero curve */
TCurve     *discZC,             /* (I) discount zero curve */
long        numPayDates,        /* (I) number of payment dates */
TDate      *accStartDates,      /* (I) [numPayDates] accrue start dates */
TDate      *accEndDates,        /* (I) [numPayDates] accrue end dates */
TDate      *payDates,           /* (I) [numPayDates] payment dates */
double      nomBpVol,           /* (I) nominal rate bp volatility */
double      realBpVol,          /* (I) real rate bp volatility */
double      nomMR,              /* (I) nominal (normal) mean reversion */
double      realMR,             /* (I) real (normal) mean reversion */
double      correlation,        /* (I) between nom and real interest rates */    
double     *pv) 		/* (O) present value */
{
    static	char	routine[] = "DriCPISwap";
    int	status = FAILURE;

    long  payIdx;

    double discFact, realFact, nomFact;
    double pValue = 0.;
    double cpiFactor, cpiReturn;
    double dcFraction;
    TDate  payDate, accStart, accEnd;
    TDate  cpiStart; /* to compute CPI factor */
    TDate  valueDate = discZC->fBaseDate;
    double covAdj = 0.;
    double payment = 0.;
    double spotCPI = cpiFixings[0];
    double lastCPI = cpiFixings[1];
    double firstCPI = cpiFixings[2];
    
    TBoolean usedFixing = FALSE; /* to make sure not more than 1 fixing required */

    TBoolean tipStyle = FALSE;

#ifdef CPI_DEBUG
    GtoLoggingSet(1); 
#endif

    /* Log inputs in fpTERM_PRN
     */
    GTO_IF_LOGGING({
	cpiSwapLogInputs(today,
			 instrType,
			 origNotl,
			 coupon,
			 dayCountConv,
			 cpiFixings,
			 nomZC,
			 realZC,
			 discZC,
/*      		 drWrap->fZcCurve,
			 drWrap->fRiskZcCurve,
			 drWrap->fDiscZcCurve,
*/
			 numPayDates,
			 accStartDates,
			 accEndDates,
			 payDates,
			 nomBpVol,
			 realBpVol,
			 nomMR,
			 realMR,
			 correlation,
			 routine);
    });

    instrType = toupper(instrType);

    /* Structure type
     */
    switch (instrType)
    {
    case 'T':
	tipStyle = TRUE;
	cpiStart = today;
	accStart = accStartDates[0];

	if(accStart>today) { /* fwd starting trade */
		cpiStart = accStart;
		firstCPI = spotCPI; /* not to adjust the orig notional */
	}

	break;

    case 'A':
    case 'C':
	break;

    default:
	GtoErrMsg("%s: Unknown instrument type (%c)."
		  " Allowed values are T/A/C.\n",
		  routine, instrType);
	goto done;
    }

    /* Check correlation
     */
    if(correlation>1. ||
       correlation<-1.)
    {
	GtoErrMsg("%s: Correlation (%f) must be within [-1,1].\n",
		  routine, correlation);
	goto done;
    }

    SHIFT_ZERO(nomMR);
    SHIFT_ZERO(realMR);

    if(numPayDates < 1)
    {
	 GtoErrMsg("%s: At least 2 reset dates must be specified.\n",
		   routine);
	goto done;
    }
    
    SHIFT_ZERO(lastCPI);
    SHIFT_ZERO(firstCPI);

    /* Calculate cash flows
     */
    for (payIdx=0; payIdx<numPayDates; payIdx++) {
	
	payDate = payDates[payIdx];

	/* Skip payment dates before value date
	 */
	if(payDate <= valueDate) 
		continue;
	
	accStart = accStartDates[payIdx];
	accEnd = accEndDates[payIdx];

	/* Compute inflation and disc factor for the corresponding dates
	 */
	if(!tipStyle) {
		cpiStart = MAX(today, accStart);
	}

	if(GtoForwardFromZCurve(nomZC,
				CPI_INTERP_TYPE,
				cpiStart,
				accEnd,
				GTO_ACT_ACT, /* not used for disc fact */
				GTO_DISCOUNT_FACTOR,
				&nomFact) != SUCCESS ||
	   GtoForwardFromZCurve(realZC,
				CPI_INTERP_TYPE,
				cpiStart,
				accEnd, 
				GTO_ACT_ACT, /* not used for disc fact */
				GTO_DISCOUNT_FACTOR,
				&realFact) != SUCCESS ||
	   GtoDiscountDate( payDate,
			    discZC,
			    CPI_INTERP_TYPE,
			    &discFact)  != SUCCESS ||
	   GtoDayCountFraction(accStart,
			       accEnd,
			       dayCountConv,
			       &dcFraction) != SUCCESS)
		goto done;

	cpiFactor = realFact/nomFact;
		
	if(tipStyle) {
		payment = coupon * dcFraction * cpiFactor;

		GTO_IF_LOGGING({
			DrlFPrintf(fpTERM_PRN,"Accrue start %s end %s Pay date %s\n\n"
				   "CPI factor %f PV factor %f Payment %f\n\n",
				   GtoFormatDate(accStart), GtoFormatDate(accEnd), 
				   GtoFormatDate(payDate), 
				   cpiFactor, discFact, 
				   payment);
		});
		
		pValue += payment * discFact;
		continue;
	}

	/* Chesk that if pay > value date, then accrue end is >= today
	 * Otherwise the calculation would require 2 fixings
	 */
	if(accEnd < today) {
		GtoErrMsg("%s: Can not compute payment on %s.\n"
			  "Both accrue start (%s) and accrue end (%s) are < today (%s).\n",
			  routine, GtoFormatDate(payDate),
			  GtoFormatDate(accStart), GtoFormatDate(accEnd), 
			  GtoFormatDate(today));
		goto done;
	}
	/* Process an already started interval, if any
	 */
	if(accStart < today) {
		if(usedFixing) {
			GtoErrMsg("%s: Can not process more than 1 interval "
				  "with accrue start < today < accrue end.\n",
				  routine);
			goto done;
		}
		usedFixing = TRUE;
		
		/* Current index percentage
		 */
		cpiReturn = spotCPI/lastCPI * cpiFactor-1.;
	}
	else { /* Need covariance adjustment */
		double T = (accStart-today)/365.25;
		double dT = (accEnd-accStart)/365.25;
		double expRealMR = exp(-realMR*T);		

		covAdj = realBpVol*(1.-exp(-realMR*dT))/realMR *
			(correlation*nomBpVol/nomMR*
			 ((1.-expRealMR)/realMR-(1.-exp(-(realMR+nomMR)*T))/(realMR+nomMR))
			 -realBpVol*(1.-expRealMR)*(1.-expRealMR)*0.5/(realMR*realMR)) ;

		cpiReturn = cpiFactor -1. + covAdj;
	}
	payment = (coupon * dcFraction * (instrType == 'C'? 1.0 : cpiReturn+1.) +
		   cpiReturn) ;

        /* Floor payment at 0
	 */
	payment = MAX(payment, 0.);

	GTO_IF_LOGGING({
		DrlFPrintf(fpTERM_PRN,"Accrue start %s end %s Pay date %s\n\n"
			"CPI factor %f Covar adj %f CPI return %f \n"
			"PV factor %f Payment %f\n\n",
			GtoFormatDate(accStart), GtoFormatDate(accEnd), 
			GtoFormatDate(payDate), 
			cpiFactor, covAdj, cpiReturn, discFact, 
			payment);
	});
	
	pValue += payment * discFact;

    } /* payIdx */

    /* Add notional payment 
     * For TIPs: notional paym is floored at par; multiply everything 
     * by CPI factor from the beginning of the deal to today
     */
    if(tipStyle) {
	    pValue *= spotCPI/firstCPI;
	    pValue += MAX(1., cpiFactor*spotCPI/firstCPI) * discFact;
    } 
    else {
	    pValue += discFact;
    }

    *pv = pValue * origNotl;

    GTO_IF_LOGGING({
	DrlFPrintf(fpTERM_PRN, "pv=%12.4f\n", *pv);
    });

    status = SUCCESS;
  
  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
    return(status);
}



/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriCPISwap}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "cpiswap_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCpiSwapW(char *dataFnam)
{
    static	char	routine[] = "DriCPISwapW";
    int	status = FAILURE;

    TDate       today;
    char        instrType;
    double      origNotl, coupon;
    TDayCount   dayCountConv;
    double      cpiFixings[CPI_NUM_FIXINGS];
    long        numPayDates;
    TDate      *accStartDates = NULL;
    TDate      *accEndDates = NULL;
    TDate      *payDates = NULL;
    double      nomVol, realVol, nomMR, realMR, correlation;
    TDate       dummyDate; 

    FILE		*fp = NULL;
    static	char	defDataFnam[] = "cpiswap_w.dat";
    static	char	todayFnam[] = "today.dat";
    static	char	eqFnam[] = "equity.dyn";

    TDrWrapperData	*drWrap = NULL;
    double		pv;

    /* Read deal data 
     */
    if (dataFnam == NULL) dataFnam = defDataFnam;
    
    if ((fp = fopen(dataFnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, dataFnam, strerror(errno));
	goto done;
    }

    READ_DATA(DRL_CHAR_T, &instrType,	        "instrument type");
    READ_DATA(DRL_DOUBLE_T, &origNotl,		"notional");
    READ_DATA(DRL_PERCENT_T, &coupon,	        "coupon");
    READ_DATA(DRL_TDAYCOUNT_T, &dayCountConv,	"day count conv");
    READ_DATA(DRL_DOUBLE_T, &cpiFixings[2],	"first CPI fixing");
    READ_DATA(DRL_DOUBLE_T, &cpiFixings[1],	"last CPI fixing");
    READ_DATA(DRL_LONG_T, &numPayDates, 	"num pay dates");

    
    if(numPayDates < 1) 
    {
	GtoErrMsg("%s: need at least 1 reset date.\n", routine);
        goto done;
    }

    if(DrlLilVectArrayFpReadV(fp, 
			      numPayDates,
			      DRL_TDATE_T,  (void*) &accStartDates,
			      DRL_TDATE_T,  (void*) &accEndDates,
			      DRL_TDATE_T,  (void*) &payDates,
			      DRL_NULL_T) == FAILURE)
    {  
        GtoErrMsg("%s: Cannot read dates array.\n", routine);
        goto done;
    }

    /* Read model data 
     */   
    READ_DATA(DRL_PERCENT_T, &nomVol,	        "vols");
    READ_DATA(DRL_PERCENT_T, &realVol,	        "real vol");
    READ_DATA(DRL_DOUBLE_T,  &nomMR,		"nominal MR");
    READ_DATA(DRL_DOUBLE_T,  &realMR,		"real MR");
    READ_DATA(DRL_DOUBLE_T,  &correlation,	"correlation");
    
    /* Close data file */ 
    if(fp) fclose(fp);

    /* Read market data
     */
    if (DriTDrWrapperDataGetFull(NULL, 
				 DRI_DRW_TYPE2_3CURVES, 
				 &drWrap) != SUCCESS)
	goto done;

    /* Read CPI spot value
     */
    if ((fp = fopen(eqFnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, eqFnam, strerror(errno));
	goto done;
    }
    READ_DATA(DRL_TDATE_T,  &dummyDate,	        "base date");
    READ_DATA(DRL_DOUBLE_T, &cpiFixings[0],	"spot value");

    if (fp) fclose(fp);

    /* Read today
     */
    if ((fp = fopen(todayFnam, "r")) == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, todayFnam, strerror(errno));
	goto done;
    }
    READ_DATA(DRL_TDATE_T,  &today,	        "today");

    if (fp) fclose(fp);

    /* Open TERM.prn
     */
    fpTERM_PRN = fopen("TERM.PRN", "w");
    if (fpTERM_PRN IS NULL)
    {
        GtoErrMsg("%s: Cannot open TERM_PRN.\n",routine);
        goto done;
    }

    GtoLoggingSet(1); 

    /* Call pricing routine 
     */
    if (DriCPISwap( today,
		    instrType,
		    origNotl,
		    coupon,
		    dayCountConv,
		    cpiFixings,
		    drWrap->fZcCurve,
		    drWrap->fRiskZcCurve,
		    drWrap->fDiscZcCurve,
		    numPayDates,
		    accStartDates,
		    accEndDates,
		    payDates,
		    nomVol,
		    realVol,
		    nomMR,
		    realMR,
		    correlation,
		    &pv) == FAILURE)
	    goto done;


    if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;
    printf("Price:     %f \n", pv);
   
	status = SUCCESS;
  done:
    if (fp) fclose(fp);
    if(fpTERM_PRN != NULL) fclose(fpTERM_PRN);

    DriTDrWrapperDataFree(drWrap);
    
    FREE(accStartDates);
    FREE(accEndDates);
    FREE(payDates);

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

PRIVATE int cpiSwapLogInputs(
TDate       today,              /* (I) spot date */
char        instrType,   	/* (I) (T)IP/(C)urrentPay/(A)ccretingCurrentPay */
double      origNotl,           /* (I) Original notional */
double      coupon,             /* (I) Coupon */
TDayCount   dayCountConv,       /* (I) For coupon payments */
double     *cpiFixings,         /* (I) [CPI_NUM_FIXINGS] Spot, last, first(for TIPs only) */
TCurve     *nomZC,              /* (I) nominal zero curve */
TCurve     *realZC,             /* (I) real zero curve */
TCurve     *discZC,             /* (I) discount zero curve */
long        numPayDates,        /* (I) Number of payment dates */
TDate      *accStartDates,      /* (I) [numPayDates] accrue start dates */
TDate      *accEndDates,        /* (I) [numPayDates] accrue end dates */
TDate      *payDates,           /* (I) [numPayDates] payment dates */
double      nomVol,             /* (I) nominal rate volatility */
double      realVol,            /* (I) real rate volatility */
double      nomMR,              /* (I) nominal (normal) mean reversion */
double      realMR,             /* (I) real (normal) mean reversion */
double      correlation,        /* (I) between nom and real interest rates */    
char       *routine)
{
    int status = FAILURE;

    int timeStamp = GtoErrMsgTimeStamp (0); /* disable temporarily */
    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN);  /* to switch back 
							  to error.log after 
							  logging is done */
    long 		i;


    GtoErrMsg("\n%s INPUTS:\n", routine);
    GtoErrMsg("\n");
    GtoErrMsg("Today:                %s\n\n", GtoFormatDate(today));
    GtoErrMsg("Instr type:           %c\n\n", instrType);
    GtoErrMsg("Original notional:    %10.4f\n\n", origNotl);
    GtoErrMsg("Coupon:               %10.4f\n\n", coupon);
    GtoErrMsg("Pay day count:        %s\n\n", GtoFormatDayCountConv(dayCountConv));
    GtoErrMsg("Spot CPI fixing:      %10.4f\n\n", cpiFixings[0]);
    GtoErrMsg("Last CPI fixing:      %10.4f\n\n", cpiFixings[1]);
    GtoErrMsg("First CPI fixing:     %10.4f\n\n", cpiFixings[2]);

    GtoErrMsg("\nAccrue start, accrue end, pay dates\n\n");
    for (i=0; i<numPayDates; i++) {
	    GtoErrMsg("%3ld:   %s   %s   %s\n", i+1, 
		      GtoFormatDate(accStartDates[i]),
		      GtoFormatDate(accEndDates[i]),
		      GtoFormatDate(payDates[i]));
    }		     
    GtoErrMsg("\n\n");
    GtoErrMsg("Nominal bp vol:       %10.4f\n\n", nomVol);
    GtoErrMsg("Real bp vol:          %10.4f\n\n", realVol);
    GtoErrMsg("Nominal MR:           %10.4f\n\n", nomMR);
    GtoErrMsg("Real MR:              %10.4f\n\n", realMR);
    GtoErrMsg("Correlation:          %10.4f\n\n", correlation);


    GtoPrintTCurve(nomZC,  "Nominal Curve");
    GtoPrintTCurve(realZC, "Real Curve");
    GtoPrintTCurve(discZC, "Discount Curve");

    GtoErrMsg("\n\n\n");

    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    GtoErrMsgTimeStamp (timeStamp);

    status = SUCCESS;

    return(status);
}
