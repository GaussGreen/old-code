/************************************************************************
 * Module:      driface
 * File:        dritridx.c
 * Function:    Total return index swap and cap/floor. Index forward price.
 * Author:      Julia Chislenko,  David Liu Oct 1998
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

#include "check.h"
#include "swapadj.h"

#include "drieq.h"
#include "eqfwd.h"   

#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h"
#include "drlsmat.h"
#include "drltime.h"
#include "dritkwrp.h"		/* TDrWrapperData routines */

#include "dritridx.h"		/* Prototype consistency */

#define xTRI_DEBUG

#ifndef TINY
#define TINY 1e-4
#endif

#define HOLIDAYFILE  "NONE"		/* currently no holiday file available
                                       * in Kapital environment */
#define BUSDAYCONV   GTO_BAD_DAY_MODIFIED
#define DRI_INTERP_TYPE GTO_LINEAR_INTERP

static	FILE	*fpTERM_PRN = NULL; 

#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
    
PRIVATE int driTotalRetSwapLogInputs(
TDate       today,              /* (I) spot date */
TCurve     *discZC,             /* (I) discount zero curve */
TCurve     *indxZC,             /* (I) index zero curve */
char        instrType,   	/* (I) 'S'wap, 'C'ap, 'F'loor */
TSwaptionMatrix2D *swMat,       /* (I) CMS swaption matrix */
TCurve     *indxVolCurve,       /* (I) index volatility curve */
double      correlation,        /* (I) between indx and interest rate */    
TEqStatData *eqStatData,	/* (I) equity.sta file		*/
char       *holidayFile,        /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
double      spotPrice,     	/* (I) Spot indx price */

long        numIdxFixing,	/* (I) Number of past refixings */
TDate      *refixingDates,	/* (I) Refixing dates		*/
TDate      *refixingEffDates,	/* (I) Refixing effective dates */
double     *indexFixings,	/* (I) Past refixings           */
long        numResetDates,      /* (I) Number of index reset dates */
TDate      *resetStDates,	/* (I) Reset start dates   */
TDate      *resetEndDates,	/* (I) Reset end dates     */
TDate      *payDates,	        /* (I) Payment dates       */
double     *strikes,            /* (I) [numResetDates] for cap/floor
				 *     NULL if swap only */
char       *routine);

PRIVATE
int MatchIdx(TDate *y, long n, TDate x, long *j);

/*f---------------------------------------------------------------------
 * Pricing routine for total return index swap/cap/floor
 * Returns SUCCESS/FAILURE.
 */
DLL_EXPORT(int)
DriTotalRetIndxSwapCapFloor(   
TDate       today,              /* (I) spot date */
TCurve     *discZC,             /* (I) discount zero curve */
TCurve     *indxZC,             /* (I) index zero curve */
char        instrType,   	/* (I) 'F'orward, 'S'wap, 'C'ap, 'F'loor */
TSwaptionMatrix2D *swMat,       /* (I) CMS swaption matrix
				 *     can be NULL if swap only */
TCurve     *indxVolCurve,       /* (I) index volatility curve
				 *     can be NULL if swap only  */
double      correlation,        /* (I) between indx and interest rate */    
long        settleDays,         /* (I) Num bus days btw indx reset and paymt*/
long        busDayConv,         /* (I) Business Day Conv    */
char       *holidayFile,        /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
double      spotPrice,     	/* (I) Spot index price	*/
long        numIdxFixing,	/* (I) Number of past refixings */
TDate      *refixingDates,	/* (I) Refixing dates		*/
TDate      *refixingEffDates,	/* (I) Refixing effective dates */
double     *indexFixings,	/* (I) Past refixings           */
long        numResetDates,      /* (I) Number of index reset dates */
TDate      *resetStDates,	/* (I) Reset start dates > today */
TDate      *resetEndDates,	/* (I) Reset end dates */
TDate      *payDates,	        /* (I) Payment date        */
double     *strikes,            /* (I) [numResetDates] for cap/floor
				 *     NULL if swap only */
double     *pv) 		/* (O) present value */
{
    static	char	routine[] = "DriTotalRetIndxSwapCapFloor";
    int      status = FAILURE;

    long     resetIdx;

    TDate    resetSt, resetEnd, resetPay;
    double   stFactor  = 0., 
             endFactor = 0., 
             pvFactor  = 0.;

    double   tExp = 0., 
             tMat = 0.; /* to interp swaption vol */
    double   bpVol   = 0., 
             indxVol = 0., 
             compVol = 0.;

    long     stFixing, endFixing;

    double   pValue = 0.;
    double   currReturn = 0., currOption = 0.;
    int	     daysInRatePeriod;
    double   daysPerYear;

    TBoolean swapOnly = FALSE;
    TBoolean cpiOnly  = FALSE;

    TBoolean closeFpTerm = FALSE;


#ifdef TRI_DEBUG
    GtoLoggingSet(1); 
#endif

    GTO_IF_LOGGING({ if (fpTERM_PRN == NULL)
        {
            fpTERM_PRN = fopen("TERM.prn", "w");
	    if (fpTERM_PRN == NULL)
	    {
	        GtoErrMsg("%s: Cannot open TERM.prn for output.\n", routine);
	        goto done;
	    }
	    closeFpTerm = TRUE;
        }});

    instrType = toupper(instrType);

    /* Structure type
     */
    switch (instrType)
    {
    case 'S':
	swapOnly = TRUE;
	break;
    case 'C':
    case 'F':
	break;
    case 'H':    /* Call on CPI return */
	cpiOnly = TRUE;
	instrType = 'C';
	break;
    case 'L':    /* Put  on CPI return */
	cpiOnly = TRUE;
	instrType = 'F';
	break;
    default:
	GtoErrMsg("%s: Unknown instrument type (%c)."
		  " Allowed values S/C/F/P/V/H/L.\n",
		  routine, instrType);
	goto done;
    }


    /*
     * To avoid affecting other modes of usage, daysPerYear was
     * only changed for CPI.  While further investigation is required,
     * it seems like it should be 365F (i.e. 365.0) for everything.
     */
    if (cpiOnly)
      daysPerYear = 365.0;
    else
      daysPerYear = 365.25;


    /* Check correlation
     */
    if(correlation>1. ||
       correlation<-1.)
    {
	GtoErrMsg("%s: Correlation (%f) must be within [-1,1].\n",
		  routine, correlation);
	goto done;
    }


    /* Calculate cash flows
     */
    for (resetIdx=0; resetIdx<numResetDates; resetIdx++)
    {
        daysInRatePeriod = resetEndDates[resetIdx] - resetStDates[resetIdx];

	/* 1. Complete past period */
	if (resetStDates[resetIdx]  < today &&
	    resetEndDates[resetIdx] < today &&
	    payDates[resetIdx]      < today)
	{
	    continue;
	}
	/* 2. Both resets are known, but pending payment */
	else 
	if (resetStDates[resetIdx]  <  today &&
	    resetEndDates[resetIdx] <  today &&
	    payDates[resetIdx]      >= today)
	{
	    /* check refixings known on both start/end dates */
	    if (MatchIdx(refixingEffDates, numIdxFixing, 
			 resetStDates[resetIdx],
			 &stFixing) == FAILURE ||
		MatchIdx(refixingEffDates, numIdxFixing, 
			 resetEndDates[resetIdx],
			 &endFixing) == FAILURE)
	    {
		GtoErrMsg("%s: missing refixings for the period with past "
			  "resets(resetSt=%s, resetEnd=%s) and pending "
			  "payment (%s) > today (%s).\n",
			  routine,
			  GtoFormatDate(resetStDates[resetIdx]),
			  GtoFormatDate(resetEndDates[resetIdx]),
			  GtoFormatDate(payDates[resetIdx]),
			  GtoFormatDate(today));
		    goto done;
	    }

	    currReturn = indexFixings[endFixing] 
		       / indexFixings[stFixing] -1.;


	    /* We annualize CPI returns (for the period) to be consistent with annualized input of strikes. */
	    if (cpiOnly)
		currReturn *= daysPerYear / daysInRatePeriod;

	    /* Option intrinsic value */
	    if (instrType == 'C')
		currReturn = MAX(currReturn - strikes[resetIdx], 0e0);
	    else if (instrType == 'F')
		currReturn = MAX(-currReturn + strikes[resetIdx], 0e0);


	    if(GtoDateFromBusDaysOffset(payDates[resetIdx],
					settleDays,
					holidayFile,
					&resetPay) == FAILURE) 
		goto done;
		    
	    /* Today/ValueDate fudge */
	    resetPay = MAX(resetPay, discZC->fBaseDate);

	    if (GtoDiscountDate(resetPay,
				discZC,
				DRI_INTERP_TYPE,
				&pvFactor) == FAILURE)
		goto done;

	    pValue += currReturn * pvFactor;

    
	}

	/* 3. In the middle of reset period, reset start is known*/
	else 
	if (resetStDates[resetIdx]  <  today &&
	    resetEndDates[resetIdx] >= today )
	{
	    /* payDate has to be in the future */
	    if (payDates[resetIdx]  < today)
	    {
		GtoErrMsg("%s: the stub period (today = %s, "
			  "(resetStart = %s, resetEnd = %s) "
			  "has passed the payment date %s!\n", 
			  routine,
			  GtoFormatDate(today),
			  GtoFormatDate(resetStDates[resetIdx]),
			  GtoFormatDate(resetEndDates[resetIdx]),
			  GtoFormatDate(payDates[resetIdx]));
		goto done;
			     
	    }

	    /* check refixings known on start dates */
	    if (MatchIdx(refixingEffDates, numIdxFixing, 
			 resetStDates[resetIdx],
			 &stFixing) == FAILURE) 
	    {
		GtoErrMsg("%s: missing refixing on reset start for period "
			  "with resetSt=%s, resetEnd=%s, and today (%s).\n",
			  routine,
			  GtoFormatDate(resetStDates[resetIdx]),
			  GtoFormatDate(resetEndDates[resetIdx]),
			  GtoFormatDate(today));
		goto done;
	    }

	    /* Settlement date adjustment */
	    if(GtoDateFromBusDaysOffset(resetEndDates[resetIdx],
					settleDays,
					holidayFile,
					&resetEnd) == FAILURE ||
	       GtoDateFromBusDaysOffset(payDates[resetIdx],
					settleDays,
					holidayFile,
					&resetPay) == FAILURE) 
		goto done;

	    /* Today/ValueDate fudge */
	    resetEnd = MAX(resetEnd, indxZC->fBaseDate);
	    resetPay = MAX(resetPay, discZC->fBaseDate);

	    if(GtoDiscountDate( resetEnd,
				indxZC,
				DRI_INTERP_TYPE,
				&endFactor) == FAILURE ||
	       GtoDiscountDate( resetPay,
				discZC,
				DRI_INTERP_TYPE,
				&pvFactor) == FAILURE)
		goto done;


	    /* Return consists resetSt -> today, and from today->resetEnd */
	    currReturn = spotPrice / indexFixings[stFixing] / endFactor - 1.;

	    if (swapOnly)
	    {
		pValue += currReturn * pvFactor;
	    }
	    else /* vol of return is index vol only */
	    {
		tExp = 0.;
		tMat = (resetEndDates[resetIdx] - today) / daysPerYear;

		/* Index price spot vol or CPI return bp vol */
		if (GtoInterpRate( resetEndDates[resetIdx],
				   indxVolCurve,
				   GTO_LINEAR_INTERP,
				   &indxVol) != SUCCESS)
		     goto done;
	
		 compVol = indxVol;
		 bpVol   = 0.;

         if (cpiOnly)
         {
	   /* 
	    * The CPI bp vol input from env is ANNUAL return bp vol,
	    * but we need to annualize the return for the period.
	    */
	     currReturn *= daysPerYear / daysInRatePeriod;
         }
         else
         {
         /* Interp swaption vol between expiry/3 and expiry
		  */
		  if (DrlTSwaptionMatrix2DInterpExpMat( swMat,
							   &bpVol,
							   tMat/3.,
							   tMat - tMat/3.,
							   FALSE) != SUCCESS)
			 goto done;

		     bpVol *= currReturn; /* bp vol of return */  
		     
		     
		     compVol =  (indxVol*indxVol +
		 	       bpVol*bpVol/3. +
		 	       correlation*bpVol*indxVol) ;
		     compVol = sqrt(compVol);
                     

         }
                                 
        /* Price option - if vol = 0 take TINY
		 */
	 	compVol = MAX(TINY, compVol);

		if(DrlNormOption( tExp+tMat,
				  currReturn,
				  compVol,
				  strikes[resetIdx], 
				  instrType=='C'?"C":"P",/* cap=call, 
							 floor=put */
				  "p",     /* price */
				  &currOption) == FAILURE)
		     goto done;
		     
		pValue += currOption * pvFactor;
	    }		 

	}
	/* 4. Forward starting period */
	else 
	{

	    /* Settlement date adjustment */
	    if(GtoDateFromBusDaysOffset(resetStDates[resetIdx],
					settleDays,
					holidayFile,
					&resetSt) == FAILURE || 
	       GtoDateFromBusDaysOffset(resetEndDates[resetIdx],
					settleDays,
					holidayFile,
					&resetEnd) == FAILURE ||
	       GtoDateFromBusDaysOffset(payDates[resetIdx],
					settleDays,
					holidayFile,
					&resetPay) == FAILURE) 
		goto done;

	    /* Today/ValueDate fudge */
	    resetSt  = MAX(resetSt,  indxZC->fBaseDate);
	    resetEnd = MAX(resetEnd, indxZC->fBaseDate);
	    resetPay = MAX(resetPay, discZC->fBaseDate);

	    if(GtoDiscountDate( resetSt,
				indxZC,
				DRI_INTERP_TYPE,
				&stFactor) == FAILURE ||
	       GtoDiscountDate( resetEnd,
				indxZC,
				DRI_INTERP_TYPE,
				&endFactor) == FAILURE ||
	       GtoDiscountDate( resetPay,
				discZC,
				DRI_INTERP_TYPE,
				&pvFactor) == FAILURE)
		goto done;


	    currReturn = stFactor/endFactor-1.;

	    /* Adjust percentage for accruing interval between reset dates
	     */
	    currReturn *= daysInRatePeriod / (double)(resetEnd - resetSt);

	    if (swapOnly)
	    {
		pValue += currReturn * pvFactor;
	    }

	    else /* compute composite vol and price an option */
	    {
		tExp = MAX((resetStDates[resetIdx]-today)/daysPerYear, 0.);
		tMat = daysInRatePeriod / daysPerYear; 
		     
		/* Index price spot vol or CPI return bp vol */ 
		if (GtoInterpRate( resetEndDates[resetIdx],
				   indxVolCurve,
				   GTO_LINEAR_INTERP,
				   &indxVol) != SUCCESS)
		     goto done;

		
		if (cpiOnly)
                {
		/* Input from eqStat is bp vol of CPI return
                 *
                 * The CPI bp vol input from env is ANNUAL return bp vol.
                 * However, we need to annualize the return for the period.
		 *
                 */
		     compVol = indxVol; 
		     currReturn /= tMat;
                }
		else                      
		{
		     /* The total vol consists of libor vol up to tExp, 
		      * plus the composite vol of spot index and return 
		      * in the reset period
		      */

		     /* Interp swaption vol
		      */
		     if (DrlTSwaptionMatrix2DInterpExpMat( swMat,
							   &bpVol,
							   tExp,
							   tMat,
							   FALSE) != SUCCESS)
			 goto done;
		     
		     bpVol *= currReturn; /* bp vol of return */  
		     
		     
		     compVol = (bpVol*bpVol*tExp  + 
			       (indxVol*indxVol +
			       bpVol*bpVol/3. +
			       correlation*bpVol*indxVol)*tMat)
			       /(tExp+tMat);
		     compVol = sqrt(compVol);
	       }		 
		     
	       /* Price option - Take TINY if compVol = 0
		*/
	
	        compVol = MAX(compVol, TINY);	
	       if(DrlNormOption( tExp+tMat,
				 currReturn,
				 compVol,
				 strikes[resetIdx], 
				 instrType=='C'?"C":"P",/* cap=call, 
							 floor=put */
				 "p",     /* price */
				 &currOption) == FAILURE)
		   goto done;
		     
	       pValue += currOption * pvFactor;
	    }
	
	}	/* end if  */

	/* Logging */
	GTO_IF_LOGGING({
		fprintf(fpTERM_PRN,"Reset %s Payment %s\n"
			"Idx factor %f Return %f PV factor %f\n\n",
			GtoFormatDate(resetStDates[resetIdx]), 
			GtoFormatDate(payDates[resetIdx]), 
			stFactor, currReturn, pvFactor);
		if(!swapOnly)
		{    
		    fprintf(fpTERM_PRN,"Vol: swaption %f compound %f\n"
			    "Strike %f Option %f\n\n",
			    bpVol, compVol, strikes[resetIdx], currOption);
		}
	});

    }		/* end for */


    *pv = pValue;

    GTO_IF_LOGGING({
	DrlFPrintf(fpTERM_PRN, "pv=%12.4f\n", *pv);
    });

    status = SUCCESS;
  
  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);

    if (closeFpTerm)
    {
	fclose(fpTERM_PRN);
	fpTERM_PRN = NULL;
    }

    return(status);
}



/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriTotalRetIndxSwapCapFloor}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "dritridx_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriIdxForwardSwapCapFloorW(char *dataFnam)
{
    static	char	routine[] = "DriIdxForwardSwapCapFloorW";
    int	status = FAILURE;

    TDate       today;
    char       instrType;
    TBoolean   forward;

    TCurve     *indxVolCurve = NULL;
    double      correlation;
    long        settleDays;
    char        holidayFile[] = HOLIDAYFILE; 
    long	busDayConv  = BUSDAYCONV;


    long        numIdxFixing;
    TDate      *refixingDates    = NULL;
    TDate      *refixingEffDates = NULL;
    double     *indexFixings     = NULL;

    long        numResetDates;
    TDate      *resetStDates  = NULL;
    TDate      *resetEndDates = NULL;
    TDate      *payDates      = NULL;
    double     *strikes = NULL;

    double     notional;
    double     payFlag;
    long       numStepUpDates;
    TDate     *stepUpDates   = NULL;
    double    *stepUpStrikes = NULL;
   
    FILE	*fp = NULL;
    static	char	defDataFnam[] = "dritridx_w.dat";
    static	char	todayFnam[] = "today.dat";
    TDrWrapperData	*drWrap = NULL;
    double		pv;

    TEqStatData		*eqStatData = NULL;
    TEqStmPrivate       *stm = NULL;
    
    long i, j;


    long    numFwd = 2;
    TDate   fwdDate[2];
    double  fwdDiscZero[2];

    double  fwdDisc;

    double  spotPrice;
    double  *fwdPrice = NULL;

    TBoolean IsIdxVolOverride = FALSE;
    char    StrOverride[512];
    TCurve  *idxVolCurveOverride = NULL;


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
    fp = NULL;


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
    READ_DATA(DRL_DOUBLE_T, &notional,		"notional");
    READ_DATA(DRL_DOUBLE_T, &payFlag,	        "payment flag");
    READ_DATA(DRL_LONG_T,   &numIdxFixing,	"num past fixing dates");

    /* 
     * If numIdxFixing = 0, forward starting trade 
     */
    if(numIdxFixing > 0) 
    {
        if(DrlLilVectArrayFpReadV(fp, 
	                	  numIdxFixing,
			          DRL_TDATE_T,   (void*) &refixingDates,
			          DRL_TDATE_T,   (void*) &refixingEffDates,
			          DRL_DOUBLE_T,  (void*) &indexFixings,
			          DRL_NULL_T) == FAILURE)
        {  
            GtoErrMsg("%s: Cannot read past refixings.\n", routine);
            goto done;
        }
    }
	

    READ_DATA(DRL_LONG_T, &numResetDates,	"num reset dates");

    if(numResetDates < 1) 
    {
	GtoErrMsg("%s: need at least 1 reset date.\n", routine);
        goto done;
    }

    if(DrlLilVectArrayFpReadV(fp, 
			      numResetDates,
			      DRL_TDATE_T,  (void*) &resetStDates,
			      DRL_TDATE_T,  (void*) &resetEndDates,
			      DRL_TDATE_T,  (void*) &payDates,
			      DRL_NULL_T) == FAILURE)
    {  
        GtoErrMsg("%s: Cannot read reset dates.\n", routine);
        goto done;
    }



    instrType = toupper(instrType);
    forward = (instrType == 'V' ||
	       instrType == 'P');

    if(!forward)  /* only used for swap, cap and floor */
    {
    	READ_DATA(DRL_LONG_T, &numStepUpDates,	"num step up dates");

    	if(numStepUpDates<1 && instrType != 'S')
    	{
	     GtoErrMsg("%s: At least one step up date must be given.\n",
		  	routine);
	     goto done;
    	}

    	if(DrlLilVectArrayFpReadV(fp, 
			          numStepUpDates,
  			          DRL_TDATE_T,  (void*) &stepUpDates,
			          DRL_PERCENT_T,  (void*) &stepUpStrikes,
			          DRL_NULL_T) == FAILURE)
    	{  
             GtoErrMsg("%s: Cannot read step up strikes.\n", routine);
             goto done;
    	}

    	/* Flat interpolate the strikes
     	 */
    	if((strikes = NEW_ARRAY(double, numResetDates)) == NULL)
	     goto done;

	 /* Strike date -> Reset End */
    	for(i=0, j=0; i<numResetDates; i++)
    	{
	     while(j<numStepUpDates-1 &&
	   	  stepUpDates[j+1]<=resetEndDates[i])
	    	j++;
	     strikes[i] = stepUpStrikes[j];
    	}

    }


    /* OPTIONAL index vol override section at the end of deal */
    IsIdxVolOverride = FALSE;
    if (fgets (StrOverride, 512, fp) != NULL)     /* EOF */
    {
        /* Space is allowed at the end of file */
        if (!isspace(StrOverride[0]))
        {
            long    numVolDates  = 0;
            TDate   *idxVolDates = NULL;
            double  *idxVolRates = NULL;

            IsIdxVolOverride = TRUE;

    	    READ_DATA(DRL_LONG_T, &numVolDates,	"num idx vol overrid dates");

            if (numVolDates > 0)
            {    
                idxVolCurveOverride = GtoNewTCurve(today, 
                                                   numVolDates, 
                                                   1, 
                                                   GTO_ACT_365F);
                if (idxVolCurveOverride == NULL) goto done;

    	        if(DrlLilVectArrayFpReadV(fp, 
	    		              numVolDates,
  			              DRL_TDATE_T,    (void*) &idxVolDates,
			              DRL_PERCENT_T,  (void*) &idxVolRates,
			              DRL_NULL_T) == FAILURE)
    	        {      
                    GtoErrMsg("%s: Cannot read idx vol override.\n", routine);
                    goto done;
    	        }
            
                /* Construct vol TCurve */
                for (i=0; i<numVolDates; i++)
                {
                    idxVolCurveOverride->fArray[i].fDate = idxVolDates[i];
                    idxVolCurveOverride->fArray[i].fRate = idxVolRates[i];
                }

                /* Free memory */
                FREE(idxVolDates);
                FREE(idxVolRates);

            }    /* end if */ 
        }
    }

 

    /* Close data file and set to NULL */ 
    if(fp) fclose(fp);
    fp = NULL;


    /* Read market data
     */
    if (DriTDrWrapperDataGetFull(NULL, 
				 DRI_DRW_TYPE2_3CURVES, 
				 &drWrap) != SUCCESS ||
	DriEqDynDataGet( NULL, "equity.dyn", &spotPrice,
			&correlation, &indxVolCurve) != SUCCESS ||
	DriTEqStatGet(NULL, "equity.sta", busDayConv, holidayFile,
			&eqStatData) != SUCCESS)
	goto done;


    /* Check for index vol override */
    if (IsIdxVolOverride)
    {
        GtoFreeTCurve(indxVolCurve);
        indxVolCurve = NULL;
        indxVolCurve = idxVolCurveOverride;
        idxVolCurveOverride = NULL;
    }


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
	driTotalRetSwapLogInputs(today,
				 drWrap->fDiscZcCurve,
				 drWrap->fRiskZcCurve,
				 instrType,
				 drWrap->fCmsSwMat,
				 indxVolCurve,
				 correlation,
				 eqStatData,
				 holidayFile,
				 spotPrice,
    				 numIdxFixing,
    				 refixingDates,
    				 refixingEffDates,
    				 indexFixings,
				 numResetDates,
				 resetStDates,
				 resetEndDates,
				 payDates,
				 strikes,
				 routine);
    });

    /* Call pricing routine 
     */
    if (forward)
    {
	if((fwdPrice = NEW_ARRAY(double, numFwd)) == NULL)
                goto done;
	
	fwdDate[0] = today;
 
        fwdDate[1] = resetStDates[0];

        if (fwdDate[1] <= today)
        {
            GtoErrMsg("%s: reset start date (%s) for forward price "
                      "calculation < today (%s).\n",
		      routine,
		      GtoFormatDate(fwdDate[1]),
		      GtoFormatDate(today));
            goto done;
        }


	/* Compute the forward prices at forward dates */
	if (DrForwardPriceGen(spotPrice,
                              eqStatData,
                              drWrap->fRiskZcCurve,     
			      busDayConv,
			      holidayFile,
                              numFwd,
                              fwdDate,
                              fwdPrice)!= SUCCESS)
                goto done;
	
	if (instrType == 'V')	/* pv of fwd price discounted using disczero */
	{
		/* 
		 * discount factors corresponding to 
		 * forward settlement dates 
		 */
		if (DriForwardDiscZero( eqStatData, 
        		        	drWrap->fDiscZcCurve,
        		        	holidayFile, 
        				numFwd, 
        				fwdDate, 
        				fwdDiscZero) != SUCCESS)
			goto done;

		fwdDisc = fwdDiscZero[1];

		pv = fwdPrice[1]*fwdDisc*notional*payFlag; 

		GTO_IF_LOGGING({
		DrlFPrintf(fpTERM_PRN, "pv:  %lf\n", pv);
		});
	}
	else   			/* forward price */
		pv = fwdPrice[1];

		GTO_IF_LOGGING({
		DrlFPrintf(fpTERM_PRN, "Forward Price:  %lf\n", pv);
		});
    }
    else
    {
	
	if (eqStatData->settleType != 'R')
	{
		GtoErrMsg("%s: only rolling settlement allowed"
			  " for total return swaps/cap/floor\n",
			  routine);
		goto done;
	}

	stm = (TEqStmPrivate *)eqStatData->stm;
	settleDays = stm->stmPeriods[0];

    	if (DriTotalRetIndxSwapCapFloor( today,	   
				         drWrap->fDiscZcCurve,
				         drWrap->fRiskZcCurve,
				         instrType,
				         drWrap->fCmsSwMat,
				         indxVolCurve,
				         correlation,
				         settleDays,
					 busDayConv,
				         holidayFile,
				         spotPrice,
    				 	 numIdxFixing,
    				 	 refixingDates,
    				 	 refixingEffDates,
    				 	 indexFixings,
				         numResetDates,
				         resetStDates,
				         resetEndDates,
				         payDates,
				         strikes,
				         &pv) == FAILURE)
	goto done;

    	pv *= notional*payFlag;
    }


    /* */
    if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;
    printf("Price:     %f \n", pv);
   
	status = SUCCESS;
  done:
    if (fp) fclose(fp);
    if(fpTERM_PRN != NULL) fclose(fpTERM_PRN);

    DriTDrWrapperDataFree(drWrap);
    DrFreeTEqStatData(eqStatData);
    
    FREE(refixingDates);
    FREE(refixingEffDates);
    FREE(indexFixings);

    FREE(resetStDates);
    FREE(resetEndDates);
    FREE(payDates);

    if (instrType == 'P')
        FREE(fwdPrice);
    else
    {
    	FREE(stepUpDates);
    	FREE(stepUpStrikes);
    	FREE(strikes);
    }

    GtoFreeTCurve(indxVolCurve);
    if (idxVolCurveOverride != NULL) GtoFreeTCurve(idxVolCurveOverride);
    

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

PRIVATE int driTotalRetSwapLogInputs(
TDate       today,              /* (I) spot date */
TCurve     *discZC,             /* (I) discount zero curve */
TCurve     *indxZC,             /* (I) index zero curve */
char        instrType,   	/* (I) 'S'wap, 'C'ap, 'F'loor */
TSwaptionMatrix2D *swMat,       /* (I) CMS swaption matrix */
TCurve     *indxVolCurve,       /* (I) index volatility curve */
double      correlation,        /* (I) between indx and interest rate */    
TEqStatData *eqStatData,	/* (I) equity.sta file		*/
char       *holidayFile,        /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
double      spotPrice,     	/* (I) Spot indx price */
long        numIdxFixing,	/* (I) Number of past refixings */
TDate      *refixingDates,	/* (I) Refixing dates		*/
TDate      *refixingEffDates,	/* (I) Refixing effective dates */
double     *indexFixings,	/* (I) Past refixings           */
long        numResetDates,      /* (I) Number of index reset dates */
TDate      *resetStDates,	/* (I) Reset start dates > today */
TDate      *resetEndDates,	/* (I) Reset end dates     */
TDate      *payDates,	        /* (I) Payment dates       */
double     *strikes,            /* (I) [numResetDates] for cap/floor
				 *     NULL if swap only */
char       *routine)
{
    int	status = FAILURE;

    int timeStamp = GtoErrMsgTimeStamp (0); /* disable temporarily */
    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN);  /* to switch back 
							  to error.log after 
							  logging is done */
    long 		i;

    char		settleType = eqStatData->settleType;
    char		divStr;
    TDividendList     	*divList   = eqStatData->divList;
    TEqStmPrivate     	*stm       = (TEqStmPrivate *)eqStatData->stm;
    long		settleDays;

    TBoolean capFloor = (toupper(instrType) == 'C' ||
			 toupper(instrType) == 'F');

    TBoolean forward = (toupper(instrType) == 'V' ||
			toupper(instrType) == 'P');

    GtoErrMsg("\n%s INPUTS:\n", routine);
    GtoErrMsg("\n");
    GtoErrMsg("Today:           %s\n", GtoFormatDate(today));
    GtoErrMsg("\n");

    GtoPrintTCurve(discZC, "Discount Curve");
    GtoErrMsg("\n");

    GtoPrintTCurve(indxZC, "Funding Curve");
    GtoErrMsg("\n");

    if(capFloor)
    {
	GtoSwaptionMatrix2DPrint(swMat, "Swaption matrix");

	GtoPrintTCurve(indxVolCurve, "Index Vol Curve");
	GtoErrMsg("\n");

	GtoErrMsg("Correlation:           %f\n", correlation);
    }

    /* Divident List */
    GtoErrMsg("\n");
    GtoErrMsg("# Dividend Dates\n\n");
    GtoErrMsg("NUMBER_OF_POINTS: %d\n", divList->fNumItems);
    GtoErrMsg("#  ==================================================\n");
    GtoErrMsg("#    exDiv Date	pay Date	Amount		Type\n");
    GtoErrMsg("#  --------------------------------------------------\n");
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
                GtoErrMsg("Unknown divident type %d for the %dth divident."
                          "\n",
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
    GtoErrMsg("Holiday file:          %s\n", holidayFile);
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
	
    GtoErrMsg("\n");

    GtoErrMsg("Structure type:        %c\n", instrType);
    GtoErrMsg("Spot Price:     %f\n", spotPrice);

    /* forward/reset dates */
    if (forward)
	GtoErrMsg("Forward Date:   %s", GtoFormatDate(resetStDates[0]));
    else
    {	
	if (numIdxFixing>0)
	{
		GtoErrMsg("Past refixings:\n");
    		for (i=0; i<numIdxFixing; i++)
    		{
			GtoErrMsg("Date [%2ld]:          %s  %s  %lf\n", 
			  	   i, 
               	            	  GtoFormatDate(refixingDates[i]),
			          GtoFormatDate(refixingEffDates[i]),
				  indexFixings[i]);
		}
	}

	GtoErrMsg("Reset schedule:\n");
    	for (i=0; i<numResetDates; i++)
    	{
		GtoErrMsg("Date [%2ld]:          %s  %s  %s", 
		  	   i, 
                           GtoFormatDate(resetStDates[i]),
                           GtoFormatDate(resetEndDates[i]),
                           GtoFormatDate(payDates[i]));
		if(capFloor)
		{
	    		GtoErrMsg("    Strike:          %f", 
		  		  strikes[i]);
		}
		GtoErrMsg("\n");
    	}
    }

    GtoErrMsg("\n");

    status = SUCCESS;

  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    GtoErrMsgTimeStamp (timeStamp);

    return(status);
}



/**
 * Give array y, find an index j such that y[j] = x
 */
PRIVATE
int MatchIdx(TDate *y, long n, TDate x, long *j)
{
        long    i;

        for (i=0; i<n; i++)
        {
                if (y[i] == x)
                {
                        *j = i;
                        return SUCCESS;
                }
        }

        return FAILURE;

}

