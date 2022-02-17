/************************************************************************
 * Module:      driface
 * File:        opfmidx.c
 * Function:    Outperformace options on indices.
 * Author:      David Liu Nov. 1998
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

#include "check.h"
#include "swapadj.h"

#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h"
#include "drlsmat.h"
#include "drltime.h"
#include "dritkwrp.h"		/* TDrWrapperData routines */

#include "eqfwd.h"		
#include "pfmopt.h"		

#define xTRI_DEBUG

#define HOLIDAYFILE  "NONE"             /* currently no holiday file available
                                       * in Kapital environment */
#define BUSDAYCONV   GTO_BAD_DAY_MODIFIED
#define DRI_INTERP_TYPE GTO_LINEAR_INTERP


#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
#undef 	NUMIDX
#define NUMIDX  2       /* 2 indices  */
    
PRIVATE int driOutPerformanceIdxOptionLogInputs(
FILE	       *fpTERM_PRN,	/* (I) TERM.PRN */
TDate          today,		/* (I) spot date */
TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
int            optType,   	/* (I) Max, Min, Diff */
TCurve         **idxVC,         /* (I) index volatility curve */
double	       corr12,		/* (I) Correlation between two indices */
TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
char           *holidayFile,    /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
long	       busDayConv,	/* (I) Business day convention  */
double         *currIdxPrice,   /* (I) Spot idx level */
double         *lastIdxPrice,   /* (I) Initial idx level */
long           numResetDates,   /* (I) Number of index reset dates */
TDate          *resetDates,	/* (I) Observation dates > today */
TDate          *payDates,	/* (I) Payment dates >= Observation dates  */
double         strike,          /* (I) strike  */
double         spread,          /* (I) spread  */
char           *routine);

/*f-@CDOC(catn="Products")---------------------------------------
  FUNCTION:       DriOutPerformanceIdxOption
 
  CREATED BY:     David Liu  Sept, 1998
 
  PURPOSE:        Outperformance Option on two indices.
 
*/
DLL_EXPORT(int)
DriOutPerformanceIdxOption(   
TDate          	today,		/* (I) spot date */
TDrWrapperData 	*drWrap,	/* (I) DR Wrapper data  */
int            	optType, 	/* (I) 1 = Max-Max 
				 *     2 = Min-Max 
				 *     3 = Max-Min   
				 *     4 = Min-Min  
				 *     5 = Spread   */
TCurve	       	**idxVC,	/* (I) Index vol curves  */
double	       	corr12,		/* (I) Correlation between two indices */
TEqStatData    	**eqStatData,	/* (I) data in equity.sta files	*/
char           	*holidayFile,  	/* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
long	       	busDayConv,	/* (I) Business day convention  */
double     	*currIdxPrice,	/* (I) Spot idx level */
double     	*lastIdxPrice,	/* (I) Initial idx level */
long        	numResetDates, 	/* (I) Number of index reset dates */
TDate      	*resetDate,	/* (I) Observation date > today */
TDate      	*payDate,	/* (I) Payment date >= Observation date */
double     	strike,        	/* (I) strike  */
double     	spread,        	/* (I) spread added to fwdIdx1  */
double     	*pv) 		/* (O) present value */
{
    static char	routine[] = "DriOutPerformanceIdxOption";
    int	   status = FAILURE;

    double *idxVol = NULL; 
    double *idxReturn = NULL; 

    double pvFactor=0.0;
    double currOption = 0.0;

    double tExp;

    TDate  fwdDate[2];
    double fwdIdxPrice[2];

    int  idx;
    int	 numIdx = (int)idxVC[0];

    TCurve **idxZC = NULL;

#ifdef TRI_DEBUG
    GtoLoggingSet(1); 
#endif

    if ((idxReturn  = NEW_ARRAY(double, numIdx)) == NULL ||
        (idxVol     = NEW_ARRAY(double, numIdx)) == NULL) 
   	goto done;
	

    /* 
     * Forward dates
     */
    fwdDate[0] = today;
    fwdDate[1] = resetDate[0];

    /* 
     * Index funding curves 
     */
    idxZC = NEW_ARRAY(TCurve *, numIdx+1);
    if (idxZC == NULL)
        goto done;

    /* 
     * Discount curve:  zero.dat (curve 1, because it requires vol )
     * Funding curve 1:	disczero.dat (curve 2)
     * Funding curve 2: riskzero.dat (curve 3)
     */
    idxZC[1] = drWrap->fDiscZcCurve;
    idxZC[2] = drWrap->fRiskZcCurve;


    /* 
     * Calculate the forward price and volatility of each index
     */
    for (idx=1; idx<=numIdx; idx++)
    {
	/* index return */
	if (DrForwardPriceGen(  currIdxPrice[idx],
			        eqStatData[idx],
				idxZC[idx],
				busDayConv,
			      	holidayFile,
				2,
				fwdDate,
				fwdIdxPrice) == FAILURE)
		goto done;

	idxReturn[idx-1] = fwdIdxPrice[1]/lastIdxPrice[idx];

	/* index volatility */ 
	if (GtoInterpRate( resetDate[0],
			   idxVC[idx],
			   GTO_LINEAR_INTERP,
			   &idxVol[idx-1]) != SUCCESS)
	     goto done;

    }	/* idx */
		 
    /* Price option. 
     * The underlying is the ratio between forward and spot, so
     * it will always be positive. The notional "1" is added
     * to the strike, but subtracted from the premium to give correct
     * option price for return indices.
     * For spread option, no changes. 
     */
    if (GtoDayCountFraction( today,
			     resetDate[0],
			     GTO_ACT_365F,
			     &tExp) == FAILURE)
	goto done;

    if(GtoBiVariOption( optType,
		        idxReturn[0],
		        idxReturn[1],
			spread,
			strike + 1.0,
			idxVol[0],	
			idxVol[1],	
			0e0,	
			corr12,
			tExp, 
			tExp, 
			"P",
			&currOption) ==FAILURE)
	goto done;

     if(GtoDiscountDate( payDate[0],
                         drWrap->fZcCurve,
                         DRI_INTERP_TYPE,
                         &pvFactor) == FAILURE)
	goto done;

    if(optType == GTO_SPREAD)
    	*pv = currOption * pvFactor;
    else
    	*pv = (currOption - 1.0) * pvFactor;


    status = SUCCESS;
  
  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
    FREE(idxReturn);
    FREE(idxVol);
    
    FREE(idxZC);

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
DriOutPerformanceIdxOptionW(char *dataFnam)
{
    static	char	routine[] = "DriOutPerformanceIdxOptionW";
    int	status = FAILURE;

    TDate       today;
    TCurve      **idxVolCurve = NULL;
    TEqStatData **eqStatData = NULL;
    int         optType;
    double	corr12;
    double	*correlation = NULL;
    char        *holidayFile = HOLIDAYFILE;
    long        busDayConv  = BUSDAYCONV;

    double      *lastIdxPrice = NULL;
    double      *currIdxPrice = NULL;
    long        numResetDates = 1;   /*  Only one observation date */
    double      *lastIdxPrice1 = NULL;
    double      *lastIdxPrice2 = NULL;
    TDate       *resetDates = NULL;
    TDate       *payDates = NULL;
    double      strike;
    double      spread = 0.0;

    double notional;
    double payFlag;
   
    FILE		*fp = NULL;
    FILE 	        *fpTERM_PRN = NULL; 

    static	char	defDataFnam[] = "pfmopt_w.dat";
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

    READ_DATA(DRL_INT_T,    &optType,	        "option type");
    READ_DATA(DRL_DOUBLE_T, &strike,      	"strike");
    READ_DATA(DRL_DOUBLE_T, &spread,      	"spread");
    READ_DATA(DRL_DOUBLE_T, &notional,	        "notional");
    READ_DATA(DRL_DOUBLE_T, &payFlag,	        "payment flag");

    /* check option type */
    switch (optType)
    {
    case 1:
	optType = GTO_MAX_MAX;
	break;
    case 2:
	optType = GTO_MIN_MAX;
	break;
    case 3:
	optType = GTO_MAX_MIN;
	break;
    case 4:
	optType = GTO_MIN_MIN;
	break;
    case 5:
	optType = GTO_SPREAD;
	break;
    default:
	GtoErrMsg("%s: Unknown option type (%d). Allowed values 1/2/3/4/5.\n",
                  routine, optType);
        goto done;
    }
	
    /* Check spread 
     */
    if(spread < 0 && optType != GTO_SPREAD){
	GtoErrMsg("%s: Spread (%lf) must be non-negative for option type %d.\n",
	  	   routine, spread, optType-4);
	goto done;
    }
	
    READ_DATA(DRL_DOUBLE_T, &corr12,	"correlation between two indices");
    /* Check correlation
     */
    if(corr12>1. || corr12<-1.){
	GtoErrMsg("%s: Correlation (%lf) between the two indices"
		  " must be within [-1,1].\n",
	  	   routine, corr12);
	goto done;
    }

    /* Last index fixing  */
    lastIdxPrice = NEW_ARRAY(double, NUMIDX+1);
    if (lastIdxPrice == NULL)
        goto done;

    lastIdxPrice[0] = (double)NUMIDX;

    if(DrlLilVectArrayFpReadV(fp, 
			      1,
			      DRL_DOUBLE_T,  (void*) &lastIdxPrice1,
			      DRL_DOUBLE_T,  (void*) &lastIdxPrice2,
			      DRL_NULL_T) == FAILURE)
    {  
        GtoErrMsg("%s: Cannot read last indices.\n", routine);
        goto done;
    }
    
    lastIdxPrice[1] = lastIdxPrice1[0];
    lastIdxPrice[2] = lastIdxPrice2[0];

    /* Forward index observation date */
    if(DrlLilVectArrayFpReadV(fp, 
			      1,	/* only one reset date */
			      DRL_TDATE_T,  (void*) &resetDates,
			      DRL_NULL_T) == FAILURE)
    {  
        GtoErrMsg("%s: Cannot read reset dates.\n", routine);
        goto done;
    }

    /* Option payment date */
    if(DrlLilVectArrayFpReadV(fp, 
			      1,	/* only one payment date */
			      DRL_TDATE_T,  (void*) &payDates,
			      DRL_NULL_T) == FAILURE)
    {  
        GtoErrMsg("%s: Cannot read reset dates.\n", routine);
        goto done;
    }


    /* Close data file */
    if(fp) fclose(fp);
	fp = NULL;

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
	fp = NULL;

    /* Check forward date starts after today
     */
    if(resetDates[0] <= today)
    {	
	GtoErrMsg("%s: forward date (%s) <= today (%s).\n",
		  routine, GtoFormatDate(resetDates[0]), GtoFormatDate(today));
	goto done;
    }

    /* Check payment date on or after observation date
     */
    if(resetDates[0] > payDates[0])
    {	
	GtoErrMsg("%s: payment date (%s) <= observation date (%s).\n",
		   routine, 
		   GtoFormatDate(payDates[0]), 
		   GtoFormatDate(resetDates[0]));
	goto done;
    }


    idxVolCurve = NEW_ARRAY(TCurve *, NUMIDX+1);
    if (idxVolCurve == NULL)
        goto done;

    idxVolCurve[0] = (void *)NUMIDX;

    correlation = NEW_ARRAY(double, NUMIDX+1);
    if (correlation == NULL)
        goto done;

    correlation[0] = (double)NUMIDX;

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
			&correlation[1], 
			&idxVolCurve[1]) != SUCCESS ||
	DriEqDynDataGet( NULL, "equity2.dyn", 
			&currIdxPrice[2],
			&correlation[2], 
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
    
/*
    currIdxPrice[1] /= lastIdxPrice[1];
    currIdxPrice[2] /= lastIdxPrice[2];
*/

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
        driOutPerformanceIdxOptionLogInputs(fpTERM_PRN,
					    today,
                                 	    drWrap,
                                 	    optType,
                                 	    idxVolCurve,
                                 	    corr12,
                                 	    eqStatData,
                                 	    holidayFile,
                                 	    busDayConv,
                                 	    currIdxPrice,
                                 	    lastIdxPrice,
                                 	    numResetDates,
                                 	    resetDates,
                                 	    payDates,
                                 	    strike,
                                 	    spread,
                                 	    routine);
    });

    /* Call pricing routine 
     */
    if (DriOutPerformanceIdxOption(today,	   
		  	           drWrap,
			           optType,
				   idxVolCurve,
				   corr12,
				   eqStatData,
				   holidayFile,
				   busDayConv,
				   currIdxPrice,
				   lastIdxPrice,
				   numResetDates,
				   resetDates,
				   payDates,
				   strike,
				   spread,
				   &pv) == FAILURE)
	goto done;

    pv *= notional*payFlag;

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
    FREE(correlation);

    if (idxVolCurve != NULL)
    {
        for (i=1;i<=NUMIDX;i++) GtoFreeTCurve(idxVolCurve[i]);
        FREE(idxVolCurve);
    }

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

PRIVATE int driOutPerformanceIdxOptionLogInputs(
FILE	       *fpTERM_PRN,	/* (I) TERM.PRN */
TDate          today,		/* (I) spot date */
TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
int            optType,   	/* (I) Max, Min, Diff */
TCurve         **idxVC,         /* (I) index volatility curve */
double	       corr12,		/* (I) Correlation between two indices */
TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
char           *holidayFile,    /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
long	       busDayConv,	/* (I) Business day convention  */
double         *currIdxPrice,  /* (I) Spot idx level */
double         *lastIdxPrice,  /* (I) Initial idx level */
long           numResetDates,   /* (I) Number of index reset dates */
TDate          *resetDates,	/* (I) Observation dates > today */
TDate          *payDates,	/* (I) Payment dates >= Observation dates  */
double         strike,          /* (I) strike  */
double         spread,          /* (I) spread  */
char           *routine)
{
    int	status = FAILURE;
    int timeStamp = GtoErrMsgTimeStamp (0); /* disable temporarily */
    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN);  /* to switch back 
							  to error.log after 
							  logging is done */
    long i;

    int  idx;
    int	 numIdx = (int)idxVC[0];

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


    for (idx=1; idx<=numIdx; idx++)
    {
  	GtoErrMsg("Index %d:\n", idx);

    	GtoPrintTCurve(idxVC[idx], "Index Vol Curve");
    	GtoErrMsg("\n");

	settleType = eqStatData[idx]->settleType;
	divList    = eqStatData[idx]->divList;
	stm        = (TEqStmPrivate *)eqStatData[idx]->stm;

    	/* Divident List */
    	GtoErrMsg("\n");
    	GtoErrMsg("# Dividend Dates\n\n");
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
    
    for (idx=1; idx<=numIdx; idx++)
    {
	GtoErrMsg("Initial Index[%d] Setting: %f\n", idx, lastIdxPrice[idx]);
    }

    for (idx=1; idx<=numIdx; idx++)
    {
	GtoErrMsg("Spot Index[%d] Level:      %f\n", idx, currIdxPrice[idx]);
    }

    /* Option Details */
   
    switch (optType)
    {
    case GTO_MAX_MAX:
    	GtoErrMsg("Option Type:              %s\n", "MAX_MAX");
	break;
    case GTO_MIN_MAX:
    	GtoErrMsg("Option Type:              %s\n", "MIN_MAX");
	break;
    case GTO_MAX_MIN:
    	GtoErrMsg("Option Type:              %s\n", "MAX_MIN");
	break;
    case GTO_MIN_MIN:
    	GtoErrMsg("Option Type:              %s\n", "MIN_MIN");
	break;
    case GTO_SPREAD:
    	GtoErrMsg("Option Type:              %s\n", "SPREAD");
	break;
    default:
	GtoErrMsg("%s: Unknown option type (%d). Allowed values 1/2/3/4/5.\n",
                  routine, optType-4);
        goto done;
    }
    GtoErrMsg("Strike:                   %f\n", strike);
    GtoErrMsg("Spread:                   %f\n", spread);
    GtoErrMsg("Correlation:              %f\n", corr12);

    GtoErrMsg("Observation Date:         %s\n", GtoFormatDate(resetDates[0]));
    GtoErrMsg("Payment Date:             %s\n", GtoFormatDate(payDates[0]));

    GtoErrMsg("\n");

    status = SUCCESS;

  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    GtoErrMsgTimeStamp (timeStamp);

    return(status);

}
