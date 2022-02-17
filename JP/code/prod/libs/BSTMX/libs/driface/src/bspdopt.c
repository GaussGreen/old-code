/*-----------------------------------------------------------------------
  SOURCE FILE:    bspdopt.c
  
  CREATED BY:     Julia Chislenko  March, 2002
  
  PURPOSE:        Driface: Basis spread option
  
  $Id: bspdopt.c 32629 2002-04-15 15:25:15Z jchislen $
  ---------------------------------------------------------------------- */
#include "drlstd.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "ldate.h"              /* GtoDayCountFraction */
#include "stub.h"               /* GTO_STUB_NONE */
#include "datelist.h"           /* GtoNewDateList */
#include "tcurve.h"             /* GtoDiscountDate */
#include "cerror.h"             /* GtoErrMsg */
#include "date_sup.h"
#include "convert.h"
#include "yearfrac.h"
#include "swaprate.h"
#include "zr2coup.h"
#include "zr2simp.h"
#include "zr2fwd.h"


#include "drlsmat.h"
#include "drltime.h"
#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h" 		/* DrlStrSubsEnv */
#include "dritkwrp.h"


#define	__DEBUG__
#undef	__DEBUG__



#undef	ARGSIZE
#define ARGSIZE(arg)	arg[0]








/*----------------------------------------------------------------------
 * Convenience routine to compute a forward rate - fixed and float have 
 * the same freq and day count
 */

static	int
computeFwdRate(
	TCurve *zcCurve,	        /* (I) zero curve */
	TCurve *discCurve,	        /* (I) discount curve */
	TDate rateEffDate,	        /* (I) rate effective date */
	TDateInterval *fixRateMat,	/* (I) fixed rate maturity interval */
	long fixRateFreq,		/* (I) fixed rate frequency */
	TDayCount fixRateDcc,	        /* (I) fixed rate day count convention */
	double *fwdRate)	        /* (O) forward swap rate */
{
static	char	routine[] = "computeFwdRate";
	TDateInterval	payInterval;
	TDate		maturityDate;
	int		status = FAILURE;


	if (GtoDtFwdAny(rateEffDate, fixRateMat, &maturityDate)
		!= SUCCESS) goto done;

	switch (fixRateFreq) {
	case 1:
	case 2:
	case 4:
	case 12:
		IF_FAILED_DONE (GtoFreq2TDateInterval( fixRateFreq, &payInterval)) ;

		IF_FAILED_DONE(GtoSwapRate2(
			discCurve,
			GTO_LINEAR_INTERP,
			rateEffDate,
			maturityDate,
			&payInterval,
			fixRateDcc,
			TRUE, /* value floating */
			0.,
			zcCurve,
			GTO_LINEAR_INTERP,
			&payInterval, /* floatRateIvl */
			fixRateDcc, /* floatRateDcc */
			FALSE, 0., /* first fixed */
			FALSE, NULL, /* convDelay vol */
			GTO_STUB_SIMPLE, FALSE, /* stub */
			GTO_BAD_DAY_NONE, GTO_BAD_DAY_NONE, GTO_BAD_DAY_NONE, /* acc/pay/reset */
			"NONE", 
			fwdRate));
		break;

	case 0:
		IF_FAILED_DONE(GtoZerosToSimplePoint(zcCurve,
			GTO_LINEAR_INTERP,
			rateEffDate,
			maturityDate,
			fixRateDcc,
			fwdRate));
		break;
	default:
		GtoErrMsg("%s: bad frequency %d.\n", routine, fixRateFreq);
		goto done;
	}

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}

/*----------------------------------------------------------------------
 * Compute a perc spread adjustment 
 */

static	int
percSpreadAdj(
	char *logfile,                  /* (I) for debug data */
	TCurve *discCurve,	        /* (I) discount curve */
	TDate expDate,	                /* (I) option exiration */
	TDate rateEffDate,	        /* (I) rate effective date */
	TDateInterval *rateMat,	        /* (I) rate maturity interval */
	TCurve     *baseVolCrv,         /* (I) used for caps, can be NULL otherwise */
	TSwaptionMatrix2D *swVolMtx,    /* (I) used for swaptions, can be NULL otherwise */
	double spdVol,                  /* (I) spread ln vol */
	double corr,                    /* (I) spread/libor corr */
	double *adj)	                /* (O) adjustment factor */
{
	static	char	routine[] = "percSpreadAdj";
	int		status = FAILURE;
	
	TDate		maturityDate;
	double          maturity, expiry, swMat=0., temp;
	double          zRate, zVol, fwdZero;

	IF_FAILED_DONE (GtoDtFwdAny(rateEffDate, rateMat, &maturityDate)) ;
	
	expiry = (expDate - baseVolCrv->fBaseDate)/365.;
	if(GtoDateIntervalToYears(rateMat,
				  &maturity) != SUCCESS)
	  goto done;

	if(corr < -1. || corr > 1.) {
		GtoErrMsg("%s: correlation (%f) out of bounds.",
			  routine, corr);
		goto done;
	}

	/* Compute zero rate (cont comp)
	 */
	IF_FAILED_DONE (GtoForwardFromZCurve( discCurve,
					      GTO_LINEAR_INTERP,
					      rateEffDate,
					      maturityDate,
					      GTO_ACT_365F,
					      GTO_CONTINUOUS_BASIS,
					      &zRate));
	fwdZero = exp(-zRate*maturity);
				 
	/* Interpolate zero rate vol
	 * if < 1y use base vols
	 * if >= 1y interp for swap maturity that has the same duration as zero
	 */
	if(maturity < 1.1) {
		IF_FAILED_DONE(GtoInterpRate( expDate,
					      baseVolCrv,
					      GTO_LINEAR_INTERP,
					      &zVol)); 
	} else { 
	       /* Skip duration calculation to get the position at rate maturity
	        */
	       /*swMat = -1./zRate * log(MAX(exp(-50.*zRate),1.-maturity*zRate));*/ /* capped at 50.*/

	        swMat = maturity;

		IF_FAILED_DONE(	DrlTSwaptionMatrix2DInterpExpMat (
			swVolMtx,
			&zVol,
			expiry,
			swMat,
			FALSE));
	}
	
	temp = (1.-fwdZero*exp(-corr*spdVol*zVol*zRate*expiry*maturity));
	if(IS_ALMOST_ZERO(temp)) {
		GtoErrMsg("%s: Numerical error while performing  perc adjustment.\n", routine);
		goto done;
	}
	*adj = (1.-fwdZero) / temp;

	/* Debug data
	 */
	DrlFilePrintf(logfile,
		      ":: Adjustment calculation ::\n"
		      "Fwd zero rate:        %f \n"
		      "Fwd pv factor:        %f \n"
		      "Equivalent swap mat:  %f \n"
		      "Interp vol:           %f \n\n",
		      zRate, fwdZero, swMat, zVol);

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}



/*----------------------------------------------------------------------
 */

int
DriBasisSpreadOptW(char *dataFnam)
{
static	char	routine[] = "DriBasisSpreadOptW";
	int	status = FAILURE;

	FILE		*fp = NULL;
	static	char		defDataFnam[] = "bspdopt_w.dat";
	static	char		logfile[] = "TERM.prn";
	TDrWrapperData	*drWrap = NULL;

	char	strBuf[GTO_MAX_STR_LEN];


	double   pvFactor, pvPremium, notional;

	TDate           expDate, optSettle;
	double		strike, expiry;

	TDate   	bsStart1, bsStart2;
	TDateInterval	bsMat1, bsMat2;
	long		bsFreq1, bsFreq2;
	TDayCount	bsDcc1, bsDcc2;
	double          weight1, weight2, corr1, corr2, bsVol1, bsVol2, spdVol;
	char            bsType;
	double		bsFwd1, bsFwd2, swFwd1, swFwd2, spd1, spd2;
	double          adj1=0., adj2=0.;
	TCurve		*bsCurve1, *bsCurve2, *swCurve;
	TBoolean        discOffLibor;

	/* Clear TERM.prn
	 */
	if ((fp = fopen(logfile, "w")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
                        routine, logfile, strerror(errno));
		goto done;
	} else fclose(fp);
	
	/*
	 * Read market data
	 */
	if (DriTDrWrapperDataGetFull(
		NULL,
		DRI_DRW_TYPE2_3CURVES,
		&drWrap) != SUCCESS)
			goto done;

	swCurve = drWrap->fZcCurve;
	bsCurve1 = drWrap->fDiscZcCurve;
	bsCurve2 = drWrap->fRiskZcCurve;

	/*
	 * Read deal data
	 */
#define	READ_DATA(type,ptr,str)	\
		{ if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
		    { GtoErrMsg("%s: can't read %s.\n", routine, str); \
		    goto done;}}

	if (dataFnam == NULL) dataFnam = defDataFnam;

	if ((fp = fopen(dataFnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
                        routine, dataFnam, strerror(errno));
		goto done;
	}

	/* Read option data
	 */
	READ_DATA(DRL_DOUBLE_T,(void*)&notional,    "notional");
	READ_DATA(DRL_TDATE_T, (void*)&expDate,     "option expiry");
	READ_DATA(DRL_TDATE_T, (void*)&optSettle,   "option settlement");
	READ_DATA(DRL_DOUBLE_T,(void*)&strike,      "strike");
	strike *= 1e-4; /* bps */

	/* Read basis rate data
	 */
	READ_DATA(DRL_DOUBLE_T,(void*)&weight1,      "weight1");
	READ_DATA(DRL_DOUBLE_T,(void*)&weight2,      "weight2");
	
	READ_DATA(DRL_TDATE_T, (void*)&bsStart1,     "basis start1");
	READ_DATA(DRL_TDATE_T, (void*)&bsStart2,     "basis start2");

	READ_DATA(DRL_CHAR_ARRAY_T, (void*)strBuf,  "basis mat1");
	IF_FAILED_DONE( GtoStringToDateInterval((void*)strBuf, routine,&bsMat1));
	READ_DATA(DRL_CHAR_ARRAY_T, (void*)strBuf,  "basis mat2");
	IF_FAILED_DONE( GtoStringToDateInterval((void*)strBuf, routine,&bsMat2));

	READ_DATA(DRL_LONG_T,(void*)&bsFreq1,        "bs rate freq1");
	READ_DATA(DRL_LONG_T,(void*)&bsFreq2,        "bs rate freq2");
	
	READ_DATA(DRL_TDAYCOUNT_T,(void*)&bsDcc1,     "bs dcc1");
	READ_DATA(DRL_TDAYCOUNT_T,(void*)&bsDcc2,     "bs dcc2");
	/*	READ_DATA(DRL_TDAYCOUNT_T,(void*)strBuf,     "bs dcc1");
	IF_FAILED_DONE( GtoStringToDayCountConv((void*)strBuf, &bsDcc1));
	READ_DATA(DRL_TDAYCOUNT_T,(void*)strBuf,     "bs dcc2");
	IF_FAILED_DONE( GtoStringToDayCountConv((void*)strBuf, &bsDcc2));
*/
	READ_DATA(DRL_CHAR_ARRAY_T, (void*)strBuf,   "basis type");
	bsType = (char)toupper(strBuf[0]);

	READ_DATA(DRL_DOUBLE_T,(void*)&corr1,        "corr1");
	READ_DATA(DRL_DOUBLE_T,(void*)&corr2,        "corr2");
	
	READ_DATA(DRL_PERCENT_T,(void*)&bsVol1,      "bsVol1");
	READ_DATA(DRL_PERCENT_T,(void*)&bsVol2,      "bsVol2");

	READ_DATA(DRL_DOUBLE_T,(void*)&spdVol,       "spdVol");
	spdVol *= 1e-4; /* bps */

	/* Expiration time
	 */
	expiry = (expDate - drWrap->fToday)/365.;
	if (expiry <=0.) {
		GtoErrMsg("%s: Expiration time (%f) < 0.\n", routine, expiry);
		goto done;
	}
	if (expDate > optSettle) {
		GtoErrMsg("%s: Option settle date (%s) > option expiry (%s).\n", 
			  routine, GtoFormatDate(optSettle), GtoFormatDate(expDate));
		goto done;
	}
	if ( bsStart1 < expDate || bsStart2 < expDate) {
		GtoErrMsg("%s: Basis start date < option expiry (%s).\n", 
			  routine, GtoFormatDate(optSettle), GtoFormatDate(expDate));
		goto done;
	}

	/* Compute fwd rates - for CMT only need par rate off CMT curve
	 */
	switch (bsType) {
	case 'P':
	case 'A':
		discOffLibor = TRUE;
		break;
	case 'T':
		discOffLibor = FALSE;
	break;
	default:
		GtoErrMsg("%s: bad basis type %c.\n", routine, bsType);
		goto done;
	}
	
	IF_FAILED_DONE( computeFwdRate(
		bsCurve1,
		(discOffLibor?swCurve:bsCurve1),
		bsStart1,
		&bsMat1,
		bsFreq1,
		bsDcc1,
		&bsFwd1));
	IF_FAILED_DONE( computeFwdRate(
		swCurve,
		swCurve,
		bsStart1,
		&bsMat1,
		bsFreq1,
		bsDcc1,
		&swFwd1));

	IF_FAILED_DONE( computeFwdRate(
		bsCurve2,
		(discOffLibor?swCurve:bsCurve2),
		bsStart2,
		&bsMat2,
		bsFreq2,
		bsDcc2,
		&bsFwd2));
	IF_FAILED_DONE( computeFwdRate(
		swCurve,
		swCurve,
		bsStart2,
		&bsMat2,
		bsFreq2,
		bsDcc2,
		&swFwd2));

	/* Calculate fwd spreads - perform adjustment for percentage
	 * No adjustments for additive spreads
	 */
	if(bsType == 'P') {
		IF_FAILED_DONE( percSpreadAdj( logfile,
					       swCurve,
					       expDate,
					       bsStart1,
					       &bsMat1,
					       drWrap->fBvCurve,
					       drWrap->fCmsSwMat,
					       bsVol1,
					       corr1,
					       &adj1));
		IF_FAILED_DONE( percSpreadAdj( logfile,
					       swCurve,
					       expDate,
					       bsStart2,
					       &bsMat2,
					       drWrap->fBvCurve,
					       drWrap->fCmsSwMat,
					       bsVol2,
					       corr2,
					       &adj2));
		spd1 = bsFwd1/swFwd1 * adj1;
		spd2 = bsFwd2/swFwd2 * adj2;
	} else {
		spd1 = swFwd1 - bsFwd1;
		spd2 = swFwd2 - bsFwd2;
	}
		
	/* Price the option
	 */
	IF_FAILED_DONE( DrlNormOption( expiry,
				       spd1*weight1+spd2*weight2,
				       spdVol,
				       strike, 
				       "C",
				       "p",     /* price */
				       &pvPremium));
	IF_FAILED_DONE (GtoForwardFromZCurve( swCurve,
					      GTO_LINEAR_INTERP,
					      swCurve->fBaseDate,
					      optSettle,
					      GTO_ACT_365F,
					      GTO_DISCOUNT_FACTOR,
					      &pvFactor));
	
	pvPremium *= pvFactor * notional;

	/* Debug data
	 */
	DrlFilePrintf(logfile, "\n\n" 
		      "Forward swap rates:   %f    %f\n"
		      "Forward basis rates:  %f    %f\n"
		      "Adjustments:          %f    %f\n"
		      "Adjusted spreads:     %f    %f\n\n",
		      swFwd1,swFwd2,  bsFwd1,bsFwd2, adj1,adj2, spd1,spd2);
	DrlFilePrintf(logfile, 
		      "Expiration:  %f \nstrike:  %f \nbpVol:  %f\n\n"
		      "Premium:  %f \npvFactor:  %f\nPV:  %f\n\n",
		      expiry, strike, spdVol, 
		      (pvPremium/pvFactor), pvFactor, pvPremium);
	
	/* Output
	 */
	GtoErrMsg("%s: NPV= %12.8f\n", routine, pvPremium);

	IF_FAILED_DONE( DriTDrWrapperDataPutPrice(pvPremium));

	status = SUCCESS;

done:
	DriTDrWrapperDataFree(drWrap);
    
	if (status != SUCCESS)
		GtoErrMsg("%s: failed (compiled %s %s).\n",
			routine, __DATE__, __TIME__);
	return(status);
}



