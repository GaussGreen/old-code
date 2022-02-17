/************************************************************************
 * Module:      driface
 * File:        cdiopt.c
 * Function:    CDI option pricer
 * Author:      Julia Chislenko July 2001

$Header$
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
#include "cfileio.h"   
#include "matswapt.h"   
#include "interp.h"   
#include "busday.h"   
#include "tcurve.h"             /* GtoDiscountDate */
#include "cerror.h"             /* GtoErrMsg */
#include "date_sup.h"
#include "convert.h"

#include "check.h"

#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h"
#include "drltime.h"
#include "drlsmat.h"
#include "dritkwrp.h"		/* TDrWrapperData routines */
#include "vnfmwrap.h"

#include "cdiopt.h"		/* Prototype consistency */

#define xCDI_DEBUG

#define BUF_SIZE 256L

static FILE 	        *fpTERM_PRN = NULL; 

#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 0.000001 : (x)))

#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
#define	SCAN_DATA(strIn,type,ptr,str)	\
    { if (DrlVTypeScan(strIn, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
#define	SCAN_DATA_2(strIn,ptr1,ptr2,str)  \
    { if (GtoSscanf(strIn, "%lf %lf", ptr1, ptr2) == FAILURE) \
          { GtoErrMsg("%s:  can't read %s.\n", routine, str); \
								 goto done;}}

#define	SCAN_DATA_3(strIn,ptr1,ptr2,ptr3,str)  \
    { if (GtoSscanf(strIn, "%lf %lf %lf", ptr1, ptr2, ptr3) == FAILURE) \
          { GtoErrMsg("%s:  can't read %s.\n", routine, str); \
								 goto done;}}

#define IS_NIL(ptr)    (strstr (ptr,"nil") != NULL)

#define NEW_WRAP_ARRAY(ptr,type,size) \
    { ptr=NEW_ARRAY(type,size+1); \
      if(ptr == NULL) { \
	  GtoErrMsg("%s: Memory allocation error.\n", routine); \
	  goto done;  \
	} \
      ptr[0] = (type)(size); \
}
#define NEW_WRAP_ARRAY_FILL(ptr,type,size,values) \
    { ptr=NEW_ARRAY(type,size+1); \
      if(ptr == NULL) { \
	  GtoErrMsg("%s: Memory allocation error.\n", routine); \
	  goto done;  \
	} else if (values) { \
          long i; \
          for(i=0; i<size; i++) ptr[i+1] = values[i]; \
        } \
      ptr[0] = (type)(size); \
}
    
PRIVATE int cdiOptLogInputs(
char        instrType,   	/* (I) P/R/S/C/F */
double      notional,           /* (I) notional */
TDate       startDate,          /* (I) accrual start */
TDate       endDate,            /* (I) accrual end = expiration for cap/floor */
double      strike,             /* (I) strike rate */
TBoolean    cashSettle,         /* (I) TRUE=cash settle FALSE=physical for swaptions */
long        dcDenom,            /* (I) day count denom for CDI compounding 
				       (act if 360 or 365, bus if 30 or 252) */
long        numCDIFixs,         /* (I) num CDI fixings betw start and value date */
double     *cdiFixings,         /* (I) CDI fixings betw start and value date */
TCurve     *zc,                 /* (I) Brazil zero curve */
long        zeroInterpType,     /* (I) zero interp type */
TCurve     *baseVolCrv,         /* (I) used for caps, can be NULL otherwise */
TSwaptionMatrix2D *swVolMtx,    /* (I) used for swaptions, can be NULL otherwise */
long        numModlPars,        /* (I) number of model params */
double     *modlParams,         /* (I) MRs, weights, corrs */
char       *holidays,           /* (I) for business day compounding */
char       *routine);

/*f---------------------------------------------------------------------
 * Pricing routine for CDI options.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCDIOpt(   
char        instrType,   	/* (I) P/R/S/C/F */
double      notional,           /* (I) notional */
TDate       startDate,          /* (I) accrual start */
TDate       endDate,            /* (I) accrual end = expiration for cap/floor */
double      strike,             /* (I) strike rate */
TBoolean    cashSettle,         /* (I) TRUE=cash settle FALSE=physical for swaptions */
long        dcDenom,            /* (I) day count denom for CDI compounding 
				       (act if 360 or 365, bus if 30 or 252) */
long        numCDIFixs,         /* (I) num CDI fixings betw start and value date */
double     *cdiFixings,         /* (I) CDI fixings betw start and value date */
TCurve     *zc,                 /* (I) Brazil zero curve */
long        zeroInterpType,     /* (I) zero interp type */
TCurve     *baseVolCrv,         /* (I) used for caps, can be NULL otherwise */
TSwaptionMatrix2D *swVolMtx,    /* (I) used for swaptions, can be NULL otherwise */
long        numModlPars,        /* (I) number of model params */
double     *modlParams,         /* (I) MRs, weights, corrs */
char       *holidays,           /* (I) for business day compounding */
double     *pv) 		/* (O) present value */
{
    static	char	routine[] = "DriCDIOption";
    int	status = FAILURE;
    
    TDate valueDate = zc->fBaseDate;
    long numPastFixs = 0, numAccrued = 0, numToStart = 0;
    double floatCFL=1., fixedCFL=1., cdiMat = 1./dcDenom, expBusTime;
    double discStart=1., discEnd=1.;
    char optType[2] = "C";
    double rateVol = 0., zeroVol = 0.;
    double fwdRate;

    long i;

    /* Parameters for average cap vol vnfm
     */
    TDateL *zcDateL = NULL;
    FloatL *zcRateL = NULL;
    FloatL *nfParamsL = NULL;
    TDateL *rateDatesL = NULL;
    FloatL *rateVolL = NULL;
    FloatL *rateMatL = NULL;
    IntL *rateFreqL = NULL;
    
#ifdef CDI_DEBUG
    GtoLoggingSet(1); 
#endif

    /* Log inputs in fpTERM_PRN
     */
    GTO_IF_LOGGING({  cdiOptLogInputs(instrType,
				      notional, 
				      startDate, 
				      endDate, 
				      strike, 
				      cashSettle,
				      dcDenom,  
				      numCDIFixs,
				      cdiFixings,
				      zc,         
				      zeroInterpType,
				      baseVolCrv,   
				      swVolMtx,  
				      numModlPars,
				      modlParams, 
				      holidays,   
				      routine);
			      });

    if(endDate <= valueDate) {
	    GtoErrMsg("%s: swap end date (%s) <= value date (%s).\n",
		      routine, GtoFormatDate(endDate), GtoFormatDate(valueDate));
	    goto done;
    }
    /* Adjust dates for holidays
     */
    if(GtoBusinessDay( startDate,
		       GTO_BAD_DAY_FOLLOW,
		       holidays,
		       &startDate) != SUCCESS ||
       GtoBusinessDay( endDate,
		       GTO_BAD_DAY_FOLLOW,
		       holidays,
		       &endDate) != SUCCESS)
       goto done;

    if(GtoBusinessDaysDiff(startDate,
			   valueDate,
			   holidays,
			   &numToStart) != SUCCESS ||
       GtoBusinessDaysDiff(startDate,
			   endDate,
			   holidays,
			   &numAccrued) != SUCCESS )
	    goto done;

   /* Float accrued before value date
     */
    if(numToStart > 0) {
	    numPastFixs = numToStart;
	    numToStart = 0;
    } else numToStart *=-1;

    if(numPastFixs > numCDIFixs) {
	    GtoErrMsg("%s: num of CDI fixings (%ld) < num of bus days before valueDate (%ld)\n",
		      routine, numCDIFixs, numPastFixs);
	    goto done;
    }
    
    for (i=1; i<=numPastFixs; i++) {
 	    floatCFL *= pow(1.+cdiFixings[numCDIFixs-i],cdiMat);
    }
    /* Full fixed accrued
     */
    fixedCFL *= pow(1.+strike, numAccrued*cdiMat);
    
    if(GtoDiscountDate( endDate,
			zc,
			zeroInterpType,
			&discEnd)  != SUCCESS ||
       (startDate>valueDate &&
	GtoDiscountDate( startDate,
			 zc,
			 zeroInterpType,
			 &discStart)  != SUCCESS))
	    goto done;
    
    /* Discounted cash flows - can use for BS since assume lognormal
     */
    floatCFL *= discStart;
    fixedCFL *= discEnd;

    fwdRate = pow(discStart/discEnd, (double)dcDenom/(numAccrued-numPastFixs)) - 1.;

    instrType = toupper(instrType);

    /* Turn non-cash delivered swaption into swap
     */ 
    if(startDate < valueDate && (instrType == 'P' || instrType == 'R')) {
	    if(cashSettle) {
		    GtoErrMsg("%s: cash settled swaption must have start date > value date.\n",
			      routine);
		    goto done;
	    } else {
		    instrType = 'S'; /* turn physical swaption into a swap */
		    notional *= instrType == 'P'? -1. : 1.; /* flip notional for the payer */
	    }
    }

    /* Structure type
     */
    switch (instrType)
    {
    case 'S': /* swap */
	   *pv = fixedCFL-floatCFL;
	   break;
    case 'R':
	    optType[0] = 'P';
    case 'P':
	    /* Interp swaption vol
	     */
	    if(DrlTSwaptionMatrix2DInterpExpMat (swVolMtx,
						 &rateVol,
						 (startDate-valueDate)/365.,
						 (endDate-startDate)/365.,
						 FALSE) != SUCCESS)
		    goto done;

	    expBusTime = (double)numToStart/dcDenom;
	    zeroVol = fwdRate *rateVol*(1.-rateVol*rateVol/24.*expBusTime)
		    * numAccrued*cdiMat/(1.+fwdRate);

	    /* Call Black-Scholes - use bus day time
	     */
	    if(DrlBlack( expBusTime,
			 floatCFL, /* fwd */
			 zeroVol,   /* atmVol */
			 fixedCFL, /* strike */
			 optType,
			 "p",     /* price */
			 pv) != SUCCESS)
		    goto done;
	    
	    GTO_IF_LOGGING({  
		    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN);  
		    GtoErrMsg("Calling DrlBlack: \n\nExp time = %10.4f\n", expBusTime);
		    GtoErrMsg("Float cfl = %10.4f\n", floatCFL);
		    GtoErrMsg("Rate vol = %10.4f\n", rateVol);
		    GtoErrMsg("Zero vol = %10.4f\n",zeroVol);
		    GtoErrMsg("Fixed cfl (strike) = %10.4f\n",fixedCFL);
		    GtoErrMsg("Opt type = %s\n\n", optType);
		    GtoErrMsg("Num accrued = %10ld num to start = %10ld\n\n",
			      numAccrued, numToStart);
		    GtoErrMsgFilePointer(prevFile); /* return to error.log */
	    });
	    
	    break;
    case 'F':
	    optType[0] = 'P';
    case 'C':
	    /* Call vnfm to compute vol of average rate
	     */
    { 
	    TDateL refDateL[2] = {1,0};
	    FloatL floatScalarsL[3] = {2.,0.5, 1./365.}; /* normal bb, O/N min */
	    TDateL startDateL[2] = {1,0};
	    TDateL endDateL[2] = {1,0};
	    FloatL avgONVolL[2] = {1.,0.};
	    long numZCdates = zc->fNumItems;
	    long numVols = baseVolCrv->fNumItems;
	    NEW_WRAP_ARRAY(zcDateL, TDateL, numZCdates);
	    NEW_WRAP_ARRAY(zcRateL, FloatL, numZCdates);
	    for(i=0; i<numZCdates; i++) {
		    zcDateL[i+1] = zc->fArray[i].fDate;
		    zcRateL[i+1] = zc->fArray[i].fRate;
	    }
	    NEW_WRAP_ARRAY_FILL(nfParamsL, FloatL, numModlPars, modlParams);
	    NEW_WRAP_ARRAY(rateDatesL, TDateL, numVols);
	    NEW_WRAP_ARRAY(rateVolL, FloatL, numVols);
	    NEW_WRAP_ARRAY(rateMatL, FloatL, numVols);
	    NEW_WRAP_ARRAY(rateFreqL, IntL, numVols);
	    for(i=0; i<numVols; i++) {
		    rateDatesL[i+1] = baseVolCrv->fArray[i].fDate;
		    rateVolL[i+1] = baseVolCrv->fArray[i].fRate;
		    rateMatL[i+1] = cdiMat;
		    rateFreqL[i+1] = 1;
	    }
	    
	    expBusTime = (double)(numToStart+numAccrued)/dcDenom;
	    
	    refDateL[1] = valueDate;
	    floatScalarsL[2] = cdiMat;
	    startDateL[1] = MAX(startDate,valueDate);
	    endDateL[1] = endDate;
	    
	    if (VnfmCalib1VNFAverageONL(refDateL,
					zcDateL,
					zcRateL,
					floatScalarsL,
					nfParamsL,
					rateDatesL,
					rateMatL,
					rateFreqL,
					rateVolL,
					startDateL,
					endDateL,
					endDateL,
					avgONVolL) != SUCCESS)
		    goto done;

	    rateVol = avgONVolL[1];
	    zeroVol = rateVol*fwdRate*(1.-rateVol*rateVol/24.*expBusTime)
		    * (double)numAccrued*cdiMat/(1.+fwdRate); 
;

	    /* Call Black-Scholes - use bus day time
	     */
	    if(DrlBlack( expBusTime,
			 floatCFL, /* fwd */
			 zeroVol,   /* atmVol */
			 fixedCFL, /* strike */
			 optType,
			 "p",     /* price */
			 pv) != SUCCESS)
		    goto done;

	    GTO_IF_LOGGING({      
		    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN); 
		    GtoErrMsg("Calling DrlBlack: \n\nExp time = %10.4f\n"
			      "Float cfl = %10.4f\n"
			      "Rate vol = %10.4f\n"
			      "Fwd rate = %10.4f\n"
			      "Zero vol = %10.4f\n"
			      "Fixed cfl (strike) = %10.4f\n"
			      "Opt type = %s\n\n"
			      "Num accrued = %10ld num to start = %10ld\n\n",
			      expBusTime, floatCFL, rateVol,
			      fwdRate, zeroVol,
			      fixedCFL, optType, numAccrued, numToStart); 
                                   
		    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    });
	    
	    break;
    default:
	    GtoErrMsg("%s: Unknown instrument type (%c)."
		      " Allowed values are P/R/S/C/F.\n",
		      routine, instrType);
	    goto done;
    }
    }
    
    *pv *= notional;
    
    GTO_IF_LOGGING({
	    DrlFPrintf(fpTERM_PRN, "pv=%12.4f\n", *pv);
    });
    
    status = SUCCESS;
    
  done:
    
    FREE(zcDateL);
    FREE(zcRateL);
    FREE(nfParamsL);
    FREE(rateDatesL);
    FREE(rateVolL);
    FREE(rateMatL);
    FREE(rateFreqL);
   
    if (status != SUCCESS)
	    GtoErrMsg("%s: failed.\n", routine);
    return(status);
}



/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriCDIOpt}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "cdiopt_w.dat" is used).
 * Returns SUCCESS/FAILURE.
 */

extern DLL_EXPORT(int)
DriCDIOptW(char *dataFnam)
{
    static	char	routine[] = "DriCDIOptW";
    int	status = FAILURE;

    char        instrType, settleType;
    double      notional, strike;
    TDate       startDate, endDate, valueDate, *cdiDates = NULL;
    long        dcDenom, zeroInterpType;
    long        numCDIFixs=0, numModlPars, numFact, numLines, firstCdiIdx = 0; 
    double      *cdiFixings = NULL, modlPars[20];
    char        zeroInterp[BUF_SIZE], holidays[BUF_SIZE], cdiFname[BUF_SIZE];
    char        wgtString[BUF_SIZE], mrString[BUF_SIZE], corrString[BUF_SIZE];
    char        path[BUF_SIZE];

    FILE		*fp = NULL;
    static	char	defDataFnam[] = "cdiopt_w.dat";

    TDrWrapperData	*drWrap = NULL;
    double		pv;

    long i;

    /* Read deal data 
     */
    if (dataFnam == NULL) dataFnam = defDataFnam;

    fp = fopen(dataFnam, "r");
    if (fp == NULL) 
    {
	GtoErrMsg("%s: can't open `%s' (%s).\n",
		  routine, dataFnam, strerror(errno));
	goto done;
    }

    READ_DATA(DRL_CHAR_T, &instrType,	        "instrument type");
    READ_DATA(DRL_DOUBLE_T, &notional,		"notional");
    READ_DATA(DRL_TDATE_T, &startDate,	        "start date");
    READ_DATA(DRL_TDATE_T, &endDate,	        "end date");
    READ_DATA(DRL_PERCENT_T, &strike,	        "strike");
    READ_DATA(DRL_CHAR_T, &settleType,	        "settlement type");
    READ_DATA(DRL_LONG_T, &dcDenom, 	        "day count denom");
    READ_DATA(DRL_CHAR_ARRAY_T, &zeroInterp,	"zeroInterp type");
    READ_DATA(DRL_CHAR_ARRAY_T, holidays,	"holiday fName");
    READ_DATA(DRL_CHAR_ARRAY_T, cdiFname,	"cdi fName");
    READ_DATA(DRL_LONG_T, &numFact, 	        "num factors");


    /* set path for holidays and CDI fixings -- if $LIVERATE_HOME 
     * is defined, append it to file name, otherwise assume the
     * files are in the working directory */
    if (getenv ("LIVERATE_HOME") != NULL)
    {
        strcpy (path, getenv ("LIVERATE_HOME"));
        strcat (path, "/");
        strcat (path, cdiFname);
        strcpy (cdiFname, path);
    
        strcpy (path, getenv ("LIVERATE_HOME"));
        strcat (path, "/");
        strcat (path, holidays);
        strcpy (holidays, path);  
    }

    if(DrlFGetLine(mrString, BUF_SIZE, fp, NULL) == NULL ||
       DrlFGetLine(wgtString, BUF_SIZE, fp, NULL) == NULL ||
       DrlFGetLine(corrString, BUF_SIZE, fp, NULL) == NULL) {
	    GtoErrMsg("%s: Error reading model param overrides.\n", routine);
	    goto done;
    }

    if(GtoStringToInterpType(zeroInterp, &zeroInterpType) != SUCCESS)
      goto done;

    /* Close data file */ 
    if (fp) {
	fclose(fp);
	fp = NULL;
    }

    /* Read market data
     */
    if (DriTDrWrapperDataGetFull(NULL, 
				 DRI_DRW_TYPE2_2CURVES, 
				 &drWrap) != SUCCESS)
	goto done;

    /* Model parameters
     */
    numModlPars = 2*numFact+(long)(numFact*(numFact-1)*0.5);

#ifdef CDI_DEBUG
    GtoErrMsg("%s: num factors: %ld  num modl params %ld.\n\n", routine, numFact, numModlPars);
    GtoErrMsg("MeanRev overrides:  %s \n", mrString);
    GtoErrMsg("Weights overrides:  %s \n", wgtString);
    GtoErrMsg("Correlation overrides:  %s \n\n", corrString);
#endif

    switch (numFact){
    case 1:
	    if(IS_NIL(mrString)) {
		    modlPars[0] = drWrap->f1Beta;
	    } else {  SCAN_DATA(mrString, DRL_DOUBLE_T, &(modlPars[0]), "mean rev");}
	    if(IS_NIL(wgtString)) {
		    modlPars[1] = drWrap->f1Weight;
	    } else {  SCAN_DATA(wgtString, DRL_DOUBLE_T, &(modlPars[1]), "weight");}
	    break;
    case 2:
	    if(IS_NIL(mrString)) {
		    modlPars[0] = drWrap->f2Beta1; 
		    modlPars[1] = drWrap->f2Beta2; 
	    } else {  
		    SCAN_DATA_2(mrString, &(modlPars[0]), &(modlPars[1]), "mean rev");}
	    if(IS_NIL(wgtString)) {
		    modlPars[2] = drWrap->f2Weight1; 
		    modlPars[3] = drWrap->f2Weight2; 
	    } else {  SCAN_DATA_2(wgtString, &(modlPars[2]), &(modlPars[3]), "weight");}
	    if(IS_NIL(corrString)) {
		    modlPars[4] = drWrap->f2Corr12; 
	    } else {  SCAN_DATA(corrString, DRL_DOUBLE_T, &(modlPars[4]), "correlation");}
	    break;
    case 3:
	    if(IS_NIL(mrString)) {
		    modlPars[0] = drWrap->f3Beta1; 
		    modlPars[1] = drWrap->f3Beta2; 
		    modlPars[2] = drWrap->f3Beta3; 
	    } else {  SCAN_DATA_3(mrString, &(modlPars[0]), &(modlPars[1]), &(modlPars[2]), 
				  "mean rev");}
	    if(IS_NIL(wgtString)) {
		    modlPars[3] = drWrap->f3Weight1; 
		    modlPars[4] = drWrap->f3Weight2; 
		    modlPars[5] = drWrap->f3Weight3; 
	    } else {  SCAN_DATA_3(wgtString, &(modlPars[3]), &(modlPars[4]), &(modlPars[5]), 
				  "weight");}
	    if(IS_NIL(corrString)) {
		    modlPars[6] = drWrap->f3Corr12; 
		    modlPars[7] = drWrap->f3Corr13; 
		    modlPars[8] = drWrap->f3Corr23; 
	    } else {  SCAN_DATA_3(corrString, &(modlPars[6]), &(modlPars[7]), &(modlPars[8]), 
				  "correlation");}
	    break;
     default:
	GtoErrMsg("%s: Num of factors (%ld) is not allowed.\n",
		  routine, numFact);
	goto done;
	    
    }
    
    /* Read CDI settings
     */
    valueDate = drWrap->fZcCurve->fBaseDate;
    numCDIFixs = MAX(valueDate - startDate,0);
    if(numCDIFixs>0) {
	    fp = fopen(cdiFname, "r");
	    if (fp IS NULL) {
		    GtoErrMsg("%s: can't open `%s' (%s).\n",
			      routine, cdiFname, strerror(errno));
		    goto done;
	    }

	    READ_DATA(DRL_LONG_T, &numLines,	        "num cdi lines");

	    if(DrlLilVectArrayFpReadV(fp, 
				      numLines,
				      DRL_TDATE_T,  (void*) &cdiDates,
				      DRL_PERCENT_T,  (void*) &cdiFixings,
				      DRL_NULL_T) == FAILURE) {  
	      GtoErrMsg("%s: Cannot read cdi fixings.\n", routine);
	      goto done;
	    }
	    numCDIFixs = 0;
	    for(i=0; i<numLines; i++) {
	      if(firstCdiIdx == 0) firstCdiIdx = i;
	      if(cdiDates[i] >= startDate &&
		       cdiDates[i] < valueDate) numCDIFixs++;
	    }
	    /* Close cdi file */ 
	    if (fp) {
		fclose(fp);
		fp = NULL;
	    }
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

    /* Call pricing routine 
     */
    if (DriCDIOpt( instrType,
		   notional,
		   startDate,
		   endDate,
		   strike,
		   (toupper(settleType) != 'P'),
		   dcDenom,
		   numCDIFixs,
		   cdiFixings+firstCdiIdx,
		   drWrap->fZcCurve,
		   zeroInterpType,
		   drWrap->fBvCurve,
		   drWrap->fCmsSwMat,
		   numModlPars,
		   modlPars,
		   holidays,
		   &pv) !=SUCCESS)
	    goto done;


    if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;
    printf("Price:     %f \n", pv);
   
	status = SUCCESS;
  done:
    if (fp) fclose(fp);
    if(fpTERM_PRN) fclose(fpTERM_PRN);

    DriTDrWrapperDataFree(drWrap);
    
    FREE(cdiFixings);
    FREE(cdiDates);

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

PRIVATE int cdiOptLogInputs(
char        instrType,   	/* (I) P/R/S/C/F */
double      notional,           /* (I) notional */
TDate       startDate,          /* (I) accrual start */
TDate       endDate,            /* (I) accrual end = expiration for cap/floor */
double      strike,             /* (I) strike rate */
TBoolean    cashSettle,         /* (I) TRUE=cash settle FALSE=physical for swaptions */
long        dcDenom,            /* (I) day count denom for CDI compounding 
				       (act if 360 or 365, bus if 30 or 252) */
long        numCDIFixs,         /* (I) num CDI fixings betw start and value date */
double     *cdiFixings,         /* (I) CDI fixings betw start and value date */
TCurve     *zc,                 /* (I) Brazil zero curve */
long        zeroInterpType,     /* (I) zero interp type */
TCurve     *baseVolCrv,         /* (I) used for caps, can be NULL otherwise */
TSwaptionMatrix2D *swVolMtx,    /* (I) used for swaptions, can be NULL otherwise */
long        numModlPars,        /* (I) number of model params */
double     *modlParams,         /* (I) MRs, weights, corrs */
char       *holidays,           /* (I) for business day compounding */
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
    GtoErrMsg("Instr type:           %c\n\n", instrType);
    GtoErrMsg("Notional:             %10.4f\n\n", notional);
    GtoErrMsg("Start date:           %s\n\n", GtoFormatDate(startDate));
    GtoErrMsg("End date:             %s\n\n", GtoFormatDate(endDate));
    GtoErrMsg("Strike:               %10.4f\n\n", strike);
    GtoErrMsg("Cash/pysical:         %s\n\n", cashSettle? "C" : "P");
    GtoErrMsg("Day count denom:      %10ld\n\n", dcDenom);
    GtoErrMsg("Zero curve interp:    %10ld\n\n", zeroInterpType);
    GtoErrMsg("Num CDI fixings:      %10ld\n\n", numCDIFixs);

    GtoErrMsg("\nFixings:\n\n");
    for (i=0; i<numCDIFixs; i++) {
	    GtoErrMsg("%3ld:   %10.4f\n", i+1, cdiFixings[i]);
    }		     
    GtoErrMsg("\n\n");

    GtoErrMsg("\nModel parameters:\n\n");
    for (i=0; i<numModlPars; i++) {
	    GtoErrMsg("    %10.4f", modlParams[i]);
    }		     
    GtoErrMsg("\n\n\n");

    GtoErrMsg("Holiday file:         %s\n\n", holidays);

    GtoPrintTCurve(zc,  "Zero Curve");
    GtoPrintTCurve(baseVolCrv,  "Base Vol Curve");
    GtoSwaptionMatrix2DPrint(swVolMtx,  "Swaption matrix");
    

    GtoErrMsg("\n\n\n");

    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    GtoErrMsgTimeStamp (timeStamp);

    status = SUCCESS;

    return(status);
}
