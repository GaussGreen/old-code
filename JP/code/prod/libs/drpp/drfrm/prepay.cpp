/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  prepay.c
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Forrest Quinn
 *			   Davis (Shuenn-Tyan) Lee
 *			   (Derivatives Research)
 *	Code version    :  1.30
 *	Extracted	:  5/19/97 at 13:54:17
 *	Last Updated	:  5/19/97 at 13:54:09
 ***************************************************************************
 *      MBS prepay module
 *      Module version:  1.1.30
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/

/* GTO hdrs */
extern "C" {
#include "cdate.h"     /* GtoMakeDateInterval */
#include "ldate.h"     /* DayCountFraction */
#include "cerror.h"    /* GtoErrMsg */
#include "dateconv.h"  /* Date format conversion */
#include "mdydate.h"   /* GtoDateToMDY */
#include "convert.h"   /* GtoFormatDate */
#include "bastypes.h"  /* TCurve, List, Array ... etc */
#include "macros.h"    /* memory allocation/deallocation */
}
/* MBS hdrs */
#include "mbsbase.h"
#include "ppconst.h"
#include "ppbase.h"
#include "mbsconst.h"
#include "prepay.h"
#include "frmpp.h"
#include "drsymboltable.h"

const DRString DAVIS_HIST_REFI30 = "$HIST_REFIRATES_DIR/refifh30.dat";
const DRString DAVIS_HIST_REFI15 = "$HIST_REFIRATES_DIR/refifh15.dat";


static char lastGtoErrStr[ERRSTRLEN] = ""; /* string to hold most recent
					   * GTO error message */

static TDate  timeStamp;
static long numHistFH30Rates = 0;
static long numHistFH15Rates = 0;
static TDate  histFH30Dates[700];
static TDate  histFH15Dates[700];
static double histFH30Rates[700];
static double histFH15Rates[700];
static double histFH30Pts[700];
static double histFH15Pts[700];

PRIVATE 
TBoolean trapGtoError(char *errStr,
		      void *userData);


/*************************************************************************
 *  PUBLIC FUNCS
 *************************************************************************/


/***************************************************************************
 *  MbsAmort()
 *  Computes array of consecutive monthly amortization,
 *  including both MBS prepays and/or MBS scheduled amortization
 *  Need following information contained in mbsDeal, mbsPrepayAssump, and
 *  mbsRateEnv data structure:
 *  <mbsDeal>
 *  mbsType - basic type of MBS (e.g., MBS_MBSTYPE_FIXED or MBS_MBSTYPE_ARM).
 *  mbsSubtype - string with supplementary info to further specify mbs type;
 *               not needed for common MBS types (i.e., set to NULL or "")--
 *               see documentation on prepay model.
 *  agency - issuing agency (e.g., MBS_AGENCY_FNMA, ...).
 *  mbsTerm - Either MBS_PP_MBSTERM_30 or MBS_PP_MBSTERM_15;
 *            determines whether model will treat collateral as 30- or
 *            15-year product.
 *  <mbsPrepayAssump>
 *  prepayModel - prepay model to use (e.g., MBS_PP_MODEL_ARM, etc).
 *  inclSchedAmort - boolean: if true, include sched. amort.; if false,
 *                   exclude scheduled amort.
 *  startDate - first month for which amort needed (day-of-month ignored).
 *  numAmortMons - # months for which to compute amort.
 *  amortForm - form for output amort--one of: MBS_PP_SPD_CPR, MBS_PP_SPD_SMM,
 *              MBS_PP_SPD_PSA.
 *  grossCpn - gross coupon rate (in decimal); should be the gross coupon
 *             accruing in startDate.
 *  wala - WALA (in months) as of startDate.
 *  warm - wgted-avg time to maturity (in months) as of startDate.
 *  <mbsRateEnv>
 *  amortIndexType - type of FHLMC commitment rate used; should be either
 *                   MBS_PP_REFI_TYPE_FH30 or MBS_PP_REFI_TYPE_FH15;
 *                   NOTE: should be consistent with mbsTerm.
 *  amortIndexStartDate - month of first amortIndexRates rate (day-of-month
 *                        ignored).
 *  numAmortIndexRates - number of amortIndexRates[] rates supplied.
 *  amortIndexRates - past & future monthly-average commitment rates
 *                    (in decimal) of FIXED FHLMC MBS OF THE SAME ORIGINAL
 *                    TERM (e.g., for mbsTerm=30yr, use FHLMC-30 rates; for
 *                    15yr, use FHLMC-15 rates);  Should begin at least 6mons
 *                    before startDate (though user can supply earlier rates
 *                    to override the hard-coded historical values stored in
 *                    this module.
 *  fhCommitPts - past & future monthly-average points, corresponding to the
 *                commitment rates in amortIndexRates[]; in decimal (e.g., 2
 *                points would be 0.02); same date range as amortIndexRates[].
 *
 ***************************************************************************/

int EXPORT MbsAmort2
(TMbsDeal     *mbsDeal,                  /* (I) Includes deal information */
 TMbsRateEnv  *mbsRateEnv,               /* (I) Includes rate information */
 TMbsPrepayDefaults  *mbsPrepayDefaults, /* (I) Includes default numbers */
 double *amort)      /* (O) array of monthly prepays (in decimal, in form
		      * determined by amortForm);  note that calling code
		      * must have allocated this array, with at least
		      * numAmortMons rate pts */
{
  static char routine[] = "MbsAmort";
  int    status = FAILURE;
  long   iCnt;
  TDate  today;

  GtoErrMsgOn();
  //    GtoErrMsgAddCallback(trapGtoError, TRUE, NULL);

  /* init mbs base lib to run "silently" */
  MbsBaseInit("JPM_PREPAY",FALSE,FALSE,NULL);
  /* @@ for testing:
     MbsBaseInit("JPM_PREPAY",TRUE,TRUE,"jpm_prepay.log");
     */

  /*
   * allocate memory space
   */
  if ((mbsRateEnv->adjAmortIndexRates = NEW_ARRAY(double,
						  mbsRateEnv->numAmortIndexRates)) IS NULL)
    {
      GtoErrMsg("%s: Failed to alloc adjAmortIndexRates[]\n",routine);
      goto done;
    }

  /*
   * Input checking for all models
   */
  if (((mbsDeal->mbsTerm == MBS_PP_MBSTERM_30) &&
       (mbsRateEnv->amortIndexType ISNT MBS_PP_REFI_TYPE_FH30)) ||
      ((mbsDeal->mbsTerm == MBS_PP_MBSTERM_15) && 
       (mbsRateEnv->amortIndexType ISNT MBS_PP_REFI_TYPE_FH15)))
    {
      GtoErrMsg("%s: mbsTerm supplied(%d) inconsistent"
		" with refi rate type(%d)\n",
		routine, mbsDeal->mbsTerm, mbsRateEnv->amortIndexType);
      goto done;
    }

  /*
   * Regardless of model, adjust input commitment rates to be 2-point
   * FH commitment rates 
   */
  if( PointAdjustCommitRate(mbsRateEnv->numAmortIndexRates,
			    mbsRateEnv->amortIndexType,
			    mbsRateEnv->amortIndexRates,
			    mbsRateEnv->fhCommitPts,
			    mbsPrepayDefaults->desiredRefiPoints,
			    mbsPrepayDefaults->ioMultiplier,
			    mbsRateEnv->adjAmortIndexRates) ISNT SUCCESS )
    {
      GtoErrMsg("%s: Failed to do point-adjust\n",routine);
      goto done;
    }


  /* 
   * Also convert the appropriate historical rate array to be 2-point
   * rates;  First, find "hard-coded" historical rates matching those
   * the user supplied */

  /* @# ask Paul: how to get data files from fimsprod to fimsdev ? */
  /* @# Soon/today, the fimsprod data should contain our files,
   * so this special handling no longer necessary */
  /* @# This is old code for PRIMUS and released on 10/03/96 - DL */
#ifdef HARDCODE_V
  InitHistAmortIndexRates(mbsRateEnv);
#else
  /* historical rafi rates only read in once for the model
   * 12/23/96 - DL
   * Be aware of updating refi. rates when time move on. 1/15/97 - Davis Lee
   */
  GtoToday(&today);
  if (numHistFH30Rates == 0 ||
      numHistFH15Rates == 0 ||
      today != timeStamp)
    {
      timeStamp = today;
      /* Read in FH30 refi rates and save to arrays */
      if (ReadHistRefiRates(MBS_PP_REFI_TYPE_FH30,
			    &mbsRateEnv->numHistAmortIndexRates,
			    &(mbsRateEnv->histAmortIndexDates),
			    &(mbsRateEnv->histAmortIndexRates),
			    &(mbsRateEnv->histFhCommitPts)) ISNT SUCCESS)
        {   
	  GtoErrMsg("%s: Failed to read in hist. refi rates\n",routine);
	  goto done;
        }   

      numHistFH30Rates = mbsRateEnv->numHistAmortIndexRates;
      for (iCnt = 0; iCnt < numHistFH30Rates; iCnt++)
        {   
	  histFH30Dates[iCnt] = mbsRateEnv->histAmortIndexDates[iCnt];
	  histFH30Rates[iCnt] = mbsRateEnv->histAmortIndexRates[iCnt];
	  histFH30Pts[iCnt] = mbsRateEnv->histFhCommitPts[iCnt];
        }   

      /* Free memory */
      FREE_ARRAY(mbsRateEnv->histAmortIndexDates);
      FREE_ARRAY(mbsRateEnv->histAmortIndexRates);
      FREE_ARRAY(mbsRateEnv->histFhCommitPts);

      /* Read in FH15 refi rates and save to arrays */
      if (ReadHistRefiRates(MBS_PP_REFI_TYPE_FH15,
			    &mbsRateEnv->numHistAmortIndexRates,
			    &(mbsRateEnv->histAmortIndexDates),
			    &(mbsRateEnv->histAmortIndexRates),
			    &(mbsRateEnv->histFhCommitPts)) ISNT SUCCESS)
        {   
	  GtoErrMsg("%s: Failed to read in hist. refi rates\n",routine);
	  goto done;
        }   

      numHistFH15Rates = mbsRateEnv->numHistAmortIndexRates;
      for (iCnt = 0; iCnt < numHistFH15Rates; iCnt++)
        {
	  histFH15Dates[iCnt] = mbsRateEnv->histAmortIndexDates[iCnt];
	  histFH15Rates[iCnt] = mbsRateEnv->histAmortIndexRates[iCnt];
	  histFH15Pts[iCnt] = mbsRateEnv->histFhCommitPts[iCnt];
        }
 
      /* Free memory */
      FREE_ARRAY(mbsRateEnv->histAmortIndexDates);
      FREE_ARRAY(mbsRateEnv->histAmortIndexRates);
      FREE_ARRAY(mbsRateEnv->histFhCommitPts);
    }
  if (mbsRateEnv->amortIndexType == MBS_PP_REFI_TYPE_FH30)
    {
      mbsRateEnv->histAmortIndexDates = NEW_ARRAY(TDate,numHistFH30Rates);
      mbsRateEnv->histAmortIndexRates = NEW_ARRAY(double,numHistFH30Rates);
      mbsRateEnv->histFhCommitPts = NEW_ARRAY(double,numHistFH30Rates);
    }
  else if (mbsRateEnv->amortIndexType == MBS_PP_REFI_TYPE_FH15)
    {
      mbsRateEnv->histAmortIndexDates = NEW_ARRAY(TDate,numHistFH15Rates);
      mbsRateEnv->histAmortIndexRates = NEW_ARRAY(double,numHistFH15Rates);
      mbsRateEnv->histFhCommitPts = NEW_ARRAY(double,numHistFH15Rates);
    }
  else if (mbsRateEnv->amortIndexType == MBS_PP_REFI_TYPE_CMT10)
    {
      GtoErrMsg("%s: CMT10 is not available\n",routine);
      goto done;
    }
 
  if (mbsRateEnv->histAmortIndexDates IS NULL)
    {
      GtoErrMsg("%s : Failed to allocate amortIndexDates\n", routine);
      goto done;
    }
  if (mbsRateEnv->histAmortIndexRates IS NULL)
    {
      GtoErrMsg("%s : Failed to allocate amortIndexRates\n", routine);
      goto done;
    }
 
  if (mbsRateEnv->amortIndexType == MBS_PP_REFI_TYPE_FH30)
    {
      mbsRateEnv->numHistAmortIndexRates = numHistFH30Rates;
      for (iCnt = 0; iCnt < mbsRateEnv->numHistAmortIndexRates; iCnt++)
        {
	  mbsRateEnv->histAmortIndexDates[iCnt] = histFH30Dates[iCnt];
	  mbsRateEnv->histAmortIndexRates[iCnt] = histFH30Rates[iCnt];
	  mbsRateEnv->histFhCommitPts[iCnt] = histFH30Pts[iCnt];
        }
    }
  else if (mbsRateEnv->amortIndexType == MBS_PP_REFI_TYPE_FH15)
    {
      mbsRateEnv->numHistAmortIndexRates = numHistFH15Rates;
      for (iCnt = 0; iCnt < mbsRateEnv->numHistAmortIndexRates; iCnt++)
        {
	  mbsRateEnv->histAmortIndexDates[iCnt] = histFH15Dates[iCnt];
	  mbsRateEnv->histAmortIndexRates[iCnt] = histFH15Rates[iCnt];
	  mbsRateEnv->histFhCommitPts[iCnt] = histFH15Pts[iCnt];
        }
    }
 
  /* @@ attached (any) day to start date; will simply this later - DL */
  mbsRateEnv->histAmortIndexStartDate =
    mbsRateEnv->histAmortIndexDates[0] * 100 + 15;
#endif

  /* Allocate a new array for adjusted hist. rates */
  if ((mbsRateEnv->adjHistAmortIndexRates = NEW_ARRAY(double,
						      mbsRateEnv->numHistAmortIndexRates)) IS NULL)
    {
      GtoErrMsg("%s: failed to alloc adjHistAmortIndexRates\n",routine);
      goto done;
    }
  /* Call central function to adjust the rates */
  if (PointAdjustCommitRate(mbsRateEnv->numHistAmortIndexRates,
			    mbsRateEnv->amortIndexType,
			    mbsRateEnv->histAmortIndexRates,
			    mbsRateEnv->histFhCommitPts,
			    mbsPrepayDefaults->desiredRefiPoints,
			    mbsPrepayDefaults->ioMultiplier,
			    mbsRateEnv->adjHistAmortIndexRates) ISNT SUCCESS)
    {
      GtoErrMsg("%s: Failed to point-adjust rates\n",routine);
      goto done;
    }

  /*
   * Branch to appropriate model function
   */

  /* if not amort at all: */
  if ((mbsDeal->mbsPrepayAssump)->prepayModel IS MBS_PP_MODEL_NONE &&
      !(mbsDeal->mbsPrepayAssump)->inclSchedAmort)
    {
      /* @@ Later: set all output amort rates to 0.0 and exit */
      GtoErrMsg("%s Not available\n",routine);
      goto done;
    }
  /* Only do schedule amort: */
  else if ((mbsDeal->mbsPrepayAssump)->prepayModel IS MBS_PP_MODEL_NONE &&
	   (mbsDeal->mbsPrepayAssump)->inclSchedAmort)
    {
      /* @@ Later: depending on whether type = fixed or ARM, compute
	 scheduled amort, and copy into output array, then exit */
      GtoErrMsg("%s : Pure scheduled-amort only not yet available\n",routine);
      goto done;
    }
  /* Branch depending on model type */
  else
    {
      switch((mbsDeal->mbsPrepayAssump)->prepayModel)
        {
        case MBS_PP_MODEL_CONST:
	  /* @@ add later */
	  GtoErrMsg("%s : Not yet available\n",routine);
	  goto done;
        case MBS_PP_MODEL_VECTOR:
	  /* @@ add later */
	  GtoErrMsg("%s : Not yet available\n",routine);
	  goto done;
        case MBS_PP_MODEL_FIX_MGRP:
	  if (mbsDeal->mbsType != MBS_MBSTYPE_FIXED)
            {
	      GtoErrMsg("%s : mbsType should be MBS_MBSTYPE_FIXED\n",
			routine);
	      goto done;
            }
	  if (frm_mgrp_prepays(mbsDeal,
			       mbsRateEnv,
			       mbsPrepayDefaults,
			       amort) IS FAILURE)
            {
	      GtoErrMsg("%s: ERROR returned in frm_mgrp_prepays.\n",
			routine);
	      goto done;
            }
	  break;
        }
    }
  status = SUCCESS;

 done:
  FREE(mbsRateEnv->adjAmortIndexRates);
  FREE(mbsRateEnv->adjHistAmortIndexRates);
  if (status IS FAILURE)
    {
      GtoErrMsg("%s : Failed\n",routine);
    }
  return(status);
}


int EXPORT InitHistAmortIndexRates
   (TMbsRateEnv *mbsRateEnv)
{
    static char routine[] = "InitHistAmortIndexRates";
    int    status = FAILURE;

    switch (mbsRateEnv->amortIndexType) 
    {
      case MBS_PP_REFI_TYPE_FH30:
        mbsRateEnv->numHistAmortIndexRates = num_hist_conv30_commit_rates;
        mbsRateEnv->histAmortIndexRates = hist_conv30_commit_rates;
        mbsRateEnv->histFhCommitPts = hist_conv30_adjust_point;
        mbsRateEnv->histAmortIndexStartDate = hist_conv30_commit_rates_start;
        break;
      case MBS_PP_REFI_TYPE_FH15:
        mbsRateEnv->numHistAmortIndexRates = num_hist_conv15_commit_rates;
        mbsRateEnv->histAmortIndexRates = hist_conv15_commit_rates;
        mbsRateEnv->histFhCommitPts = hist_conv15_adjust_point;
        mbsRateEnv->histAmortIndexStartDate = hist_conv15_commit_rates_start;
        break;
      default:
        GtoErrMsg("%s: Invalid refi type: %d\n",
            routine,mbsRateEnv->amortIndexType);
        goto done;
    }
 
    status = SUCCESS;
 
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}

int EXPORT InitHistAmortAvgIndexRates
   (TMbsRateEnv *mbsRateEnv,
    TMbsDeal    *mbsDeal)
{
    static char routine[] = "InitHistAmortAvgIndexRates";
    int    status = FAILURE;

    switch (mbsRateEnv->amortIndexType) 
    {
    case MBS_PP_REFI_TYPE_FH30:
        (mbsDeal->mbsPrepayAssump)->armAvgHistRefi=DEF_MBS_PP_ARM_HIST_FH30_IDX;
        break;
    case MBS_PP_REFI_TYPE_FH15:
        (mbsDeal->mbsPrepayAssump)->armAvgHistRefi=DEF_MBS_PP_ARM_HIST_FH15_IDX;
        break;
    case MBS_PP_REFI_TYPE_CMT10:
        (mbsDeal->mbsPrepayAssump)->armAvgHistRefi=DEF_MBS_PP_ARM_HIST_CMT10_IDX;
        break;
    default:
        GtoErrMsg("%s: Invalid refi type: %d\n",
            routine,mbsRateEnv->amortIndexType);
        goto done;
    }
 
    status = SUCCESS;
 
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}

/***************************************************************************
 * InitMbsTreeAmort
 * Init function to be called (once) before 
 * using "tree" prepay function MbsTreeAmort() (pre-computes some 
 * parameters, thus speeding up calculation)
 ***************************************************************************/
int EXPORT InitMbsTreeAmort
   (TMbsDeal *mbsDeal,             /* (I) deal to be priced */
    long   refiIndexRateType,      /* (I) type of refi index rate
                                    * e.g., MBS_PP_REFI_TYPE_CMT10 */
    /* @# may not need this */
    double grCpnSprdOffRefiIx)     /* (I) spread (decimal) to add to refi idx*/
{
    F_INIT("InitMbsTreeAmort");
    long iSeas;
    long iLog;
    long iGroup;

    /* Seasonality */
    double mgrp_gnma_seasonality[DEF_MBS_PP_NUM_SEASONALITY];
    double mgrp_fnma_seasonality[DEF_MBS_PP_NUM_SEASONALITY];
    /* S-curve params */
    double mgrp_gnma_scurve_logist[DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST];
    double mgrp_fnma_scurve_logist[DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST];
    /* Seasoning ramp params */
    double mgrp_gnma_seas_logist[DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST];
    double mgrp_fnma_seas_logist[DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST];
    /* Group params */
    double mgrp_gnma_group_ps[DEF_MBS_PP_MGRP_NUM_GROUPPS];
    double mgrp_fnma_group_ps[DEF_MBS_PP_MGRP_NUM_GROUPPS];

    for (iSeas = 0; iSeas < DEF_MBS_PP_NUM_SEASONALITY; iSeas++)
    {
        mgrp_gnma_seasonality[iSeas]=DEF_MBS_PP_MGRP_GNMA_SEASONALITY[iSeas];
        mgrp_fnma_seasonality[iSeas]=DEF_MBS_PP_MGRP_FNMA_SEASONALITY[iSeas];
    }
    for (iLog = 0; iLog < DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST; iLog++)
    {
        mgrp_gnma_scurve_logist[iLog]=DEF_MBS_PP_MGRP_GNMA_SCURVE_LOGIST[iLog];
        mgrp_fnma_scurve_logist[iLog]=DEF_MBS_PP_MGRP_FNMA_SCURVE_LOGIST[iLog];
    }
    for (iSeas = 0; iSeas < DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST; iSeas++)
    {
        mgrp_gnma_seas_logist[iSeas]=DEF_MBS_PP_MGRP_GNMA_SCURVE_LOGIST[iSeas];
        mgrp_fnma_seas_logist[iSeas]=DEF_MBS_PP_MGRP_FNMA_SCURVE_LOGIST[iSeas];
    }
    for (iGroup = 0; iGroup < DEF_MBS_PP_MGRP_NUM_GROUPPS; iGroup++)
    {
        mgrp_gnma_group_ps[iGroup]=DEF_MBS_PP_MGRP_GNMA_SEAS_LOGIST[iGroup];
        mgrp_fnma_group_ps[iGroup]=DEF_MBS_PP_MGRP_FNMA_SEAS_LOGIST[iGroup];
    }

    /* @# maybe toss this */
    mbsDeal->mbsGrCcSprd[refiIndexRateType] = grCpnSprdOffRefiIx;

    /* 
     * Set parameters depending on model to be used 
     */

    if ( (mbsDeal->mbsPrepayAssump)->prepayModel IS MBS_PP_MODEL_ARM )
    {
        /* For ARM model, we can only handle GNMA II, for now */
        XONTRUE(mbsDeal->mbsAgency ISNT MBS_AGENCY_GNMAII,
                "Presently, can only use ARM model for GNMA II" );
        /* Use GNMA II seasonality from ARM model */
        XONFAIL( CopyDouble(DEF_MBS_PP_NUM_SEASONALITY,
                            DEF_MBS_PP_ARM_SEASONALITY,
                            mbsDeal->mbsPrepayAssump->seasonality) );
        (mbsDeal->mbsPrepayAssump)->armAbsRateEffect =
            DEF_MBS_PP_ARM_ABS_RATE_EFF;

        switch (refiIndexRateType)
        {
        case MBS_PP_REFI_TYPE_FH30:
            /* Start with avg historical fh30 rate (ignoring variation
             * in pts from month-to-month */
            (mbsDeal->mbsPrepayAssump)->armAvgHistRefi =
                DEF_MBS_PP_ARM_HIST_FH30_IDX;
            break;
        case MBS_PP_REFI_TYPE_FH15:
            /* Start with avg historical fh15 rate (ignoring variation
             * in pts from month-to-month */
            (mbsDeal->mbsPrepayAssump)->armAvgHistRefi =
                DEF_MBS_PP_ARM_HIST_FH15_IDX;
            break;
        case MBS_PP_REFI_TYPE_CMT10:
            /* Start with avg historical CMT10 rate */
            (mbsDeal->mbsPrepayAssump)->armAvgHistRefi =
                DEF_MBS_PP_ARM_HIST_CMT10_IDX;
            break;
        default:
            GtoErrMsg("%s: Invalid refi index rate type: %ld\n",
                      routine, refiIndexRateType);
            break;
        }
        (mbsDeal->mbsPrepayAssump)->armSmmLogist[0] =
            DEF_MBS_PP_ARM_SMM_LOGIST_RIGHT;
        (mbsDeal->mbsPrepayAssump)->armSmmLogist[1] =
            DEF_MBS_PP_ARM_SMM_LOGIST_LEFT;
        (mbsDeal->mbsPrepayAssump)->armSmmLogist[2] =
            DEF_MBS_PP_ARM_SMM_LOGIST_WIDTH;
        (mbsDeal->mbsPrepayAssump)->armSmmLogist[3] =
            DEF_MBS_PP_ARM_SMM_LOGIST_INFLEC;

        (mbsDeal->mbsPrepayAssump)->armSeasLogist[0] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_RIGHT;
        (mbsDeal->mbsPrepayAssump)->armSeasLogist[1] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_LEFT;
        (mbsDeal->mbsPrepayAssump)->armSeasLogist[2] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_WIDTH;
        (mbsDeal->mbsPrepayAssump)->armSeasLogist[3] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_INFLEC;

        /* Also pre-compute params to use with faster logistic
           function fast_logistic2() */
        (mbsDeal->mbsPrepayAssump)->armRefiIncShift = DEF_MBS_PP_ARM_INC_SHIFT;
        (mbsDeal->mbsPrepayAssump)->armNewWacRefiIncShift = DEF_MBS_PP_ARM_NEW_WAC_INC_SHIFT;

        (mbsDeal->mbsPrepayAssump)->armSmmFastLogist[0] =
            DEF_MBS_PP_ARM_SMM_LOGIST_RIGHT;
        (mbsDeal->mbsPrepayAssump)->armSmmFastLogist[1] =
            DEF_MBS_PP_ARM_SMM_LOGIST_LEFT *
            exp(DEF_MBS_PP_ARM_SMM_LOGIST_INFLEC/
            DEF_MBS_PP_ARM_SMM_LOGIST_WIDTH);
        (mbsDeal->mbsPrepayAssump)->armSmmFastLogist[2] =
            DEF_MBS_PP_ARM_SMM_LOGIST_WIDTH;
        (mbsDeal->mbsPrepayAssump)->armSmmFastLogist[3] = 
            exp(DEF_MBS_PP_ARM_SMM_LOGIST_INFLEC/
            DEF_MBS_PP_ARM_SMM_LOGIST_WIDTH);

        (mbsDeal->mbsPrepayAssump)->armSeasFastLogist[0] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_RIGHT;
        (mbsDeal->mbsPrepayAssump)->armSeasFastLogist[1] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_LEFT *
            exp(DEF_MBS_PP_ARM_SEAS_LOGIST_INFLEC/
            DEF_MBS_PP_ARM_SEAS_LOGIST_WIDTH);
        (mbsDeal->mbsPrepayAssump)->armSeasFastLogist[2] =
            DEF_MBS_PP_ARM_SEAS_LOGIST_WIDTH;
        (mbsDeal->mbsPrepayAssump)->armSeasFastLogist[3] = 
            exp(DEF_MBS_PP_ARM_SEAS_LOGIST_INFLEC/
            DEF_MBS_PP_ARM_SEAS_LOGIST_WIDTH);
    }
    else if ((mbsDeal->mbsPrepayAssump)->prepayModel IS MBS_PP_MODEL_FIX_MGRP)
    {
        switch (mbsDeal->mbsAgency)
        {
        case MBS_AGENCY_GNMAI:
        case MBS_AGENCY_GNMAII:
            XONFAIL( CopyDouble(DEF_MBS_PP_NUM_SEASONALITY,
                                mgrp_gnma_seasonality,
                                mbsDeal->mbsPrepayAssump->seasonality) );
            XONFAIL( CopyDouble(DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST,
                                mgrp_gnma_scurve_logist,
                                mbsDeal->mbsPrepayAssump->mgrpScurveLogist) );
            XONFAIL( CopyDouble(DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST,
                                mgrp_gnma_seas_logist,
                                mbsDeal->mbsPrepayAssump->mgrpSeasLogist) );
            break;
        case MBS_AGENCY_FNMA:
        case MBS_AGENCY_FHLMC:
        case MBS_AGENCY_GOLD:
        case MBS_AGENCY_WHOLE:
            XONFAIL( CopyDouble(DEF_MBS_PP_NUM_SEASONALITY,
                                mgrp_fnma_seasonality,
                                mbsDeal->mbsPrepayAssump->seasonality) );
            XONFAIL( CopyDouble(DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST,
                                mgrp_fnma_scurve_logist,
                                mbsDeal->mbsPrepayAssump->mgrpScurveLogist) );
            XONFAIL( CopyDouble(DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST,
                                mgrp_fnma_seas_logist,
                                mbsDeal->mbsPrepayAssump->mgrpSeasLogist) );
            break;
        default:
            XONTRUE( TRUE, "Invalid agency");
        }
    }
    else
    {
        XONTRUE( TRUE, "Tree prepay function cannot use this prepay model");
    }


    F_END;
}




/************************************************************
 * SchedAmortSMM()
 * Prepay utility to compute scheduled monthly amortization
 * (in decimal SMM) for a level-pay mortgage
 *   Note that here we use the "prepay" time convention--that is,
 * the prepays/amortization of a given month (the "amort. month")
 * describe the drop from the factor of the prev. month to that of 
 * curr. month (that is, change between starting balance of prev. 
 * month and starting balance of curr. month
 ************************************************************/
int EXPORT SchedAmortSMM
   (double      grossCpn,       /* (I) gross coupon (decimal) accruing
                                   in month before amort. month */
    long        warm,           /* (I) # months remaining to maturity,
                                   as of amort. month (e.g., for first
                                   amort. of a mortgage, would be
                                   term-1) */
    double     *smm)            /* (O) scheduled amort., in decimal SMM */
{
    double    mwac;             /*  monthly gross cpn */

    mwac = grossCpn/12.;
    *smm = mwac/(pow(1.+mwac,(double)warm+1.)-1.);

    return(SUCCESS);
}


/***************************************************************************
 *  GetPrepayForPeriod
 *  Prepay utility to compute the prepay rate over a period of several
 *  months, based on array of monthly prepays.
 *  FUNCTION RETURNS: 0=success, nonzero=failure
 *  Note: if output prepay is to be in PSA, the associated age will be age1.
 ***************************************************************************/
int EXPORT GetPrepayForPeriod
   (long     num_prepays,       /* (I) # of prepays in prepays[] array */
    long     amortForm,        /* (I) form of prepays in prepays[]:
                                  MBS_PP_SPD_PSA, MBS_PP_SPD_CPR, or
                                  MBS_PP_SPD_SMM */
    double *prepays,           /* (I) array of monthly prepays, in decimal */
    long     imon1,             /* (I) starting month of period (between 0
                                  and num_prepays-1) */
    long     imon2,             /* (I) ending month of period (>= imon1) */
    long     age1,              /* (I) age as of imon1 (needed if PSA spds
                                  involved) */
    long     final_prepay_form, /* (I) form for output prepay rate:
                                  MBS_PP_SPD_PSA, MBS_PP_SPD_CPR, or
                                  MBS_PP_SPD_SMM */
    double *prepay_over_period)/* (O) prepay rate over period */
{
    static char routine[] = "GetPrepayForPeriod";
    int    status = FAILURE;
    long    imon;
    double factor;
    double smm;

    /* if not already done, init lib */

    if (num_prepays < 1)
    {
        GtoErrMsg("%s: Bad # input prepays supplied\n",routine);
        goto done;
    }
    if (prepays IS NULL)
    {
        GtoErrMsg("%s: No input prepays supplied\n",routine);
        goto done;
    }
    if (!is_valid_prepay_form(amortForm))
    {
        GtoErrMsg("%s: Bad prepay form for input prepays\n",routine);
        goto done;
    }
    if (imon1 > imon2)
    {
        GtoErrMsg("%s: Bad order of start/end month range\n",routine);
        goto done;
    }
    if (imon1 < 0 ||
        imon1 >= num_prepays ||
        imon2 < 0 ||
        imon2 >= num_prepays)
    {
        GtoErrMsg("%s: Invalid month start/end range\n",routine);
        goto done;
    }
    if (!is_valid_prepay_form(final_prepay_form))
    {
        GtoErrMsg("%s: Bad prepay form for output prepays\n",routine);
        goto done;
    }

    /* compute net paydown over period */
    factor = 1.;
    for (imon = imon1; imon <= imon2; imon++)
    {
        /* convert to SMM */
        if (ConvertPrepays(amortForm,
                           MBS_PP_SPD_SMM,
                           prepays[imon],
                           age1+(imon-imon1),
                          &smm) IS FAILURE)
        {
            GtoErrMsg("%s: ERROR returned in ConvertPrepays\n",routine);
            goto done;
        }
        factor *= 1.-smm;
    }
    /* compute smm from this paydown */
    smm = 1. - mbs_pow(factor,1./(double) (imon2-imon1+1));
    /* convert to desired form */
    if (ConvertPrepays(MBS_PP_SPD_SMM,
                       final_prepay_form,
                       smm,age1,
                       prepay_over_period) IS FAILURE)
    {
        GtoErrMsg("%s: ERROR returned in ConvertPrepays\n",routine);
        goto done;
    }
 
    status = SUCCESS;

done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/************************************************************************
 *  ReadHistRefiRates()
 *  This function is to read in hist. refi rates from ascii file and
 *  store dates and rates to array.
 *  Files should locate at $DATA directory for all platforms.
 ************************************************************************/
int EXPORT ReadHistRefiRates
   (long      amortIndexType,     /* (I) type of FHLMC commitment rate used;
                                    * should be either MBS_PP_REFI_TYPE_FH30
                                    * or MBS_PP_REFI_TYPE_FH15;  NOTE: should
                                    * be consistent with mbsTerm */
    long     *numAmortIndexRates, /* (O) # of hoist. amortIndexRates */
    TDate  **amortIndexDates,    /* (O) array of dates of amortIndexRates */
    double **amortIndexRates,    /* (O) hist. monthly-average commitment rates
                                  * (in decimal) of FIXED FHLMC MBS OF THE
                                  * SAME ORIGINAL TERM (e.g., for mbsTerm=30yr,
                                  * use FHLMC-30 rates; for 15yr, use FHLMC-15
                                  * rates) */
    double **fhCommitPts)        /* (O) hist. monthly-average points,
                                  * corresponding to the commitment rates in
                                  * amortIndexRates[]; in decimal (e.g., 2
                                  * points would be 0.02); same date range
                                  * as amortIndexRates[] */
{
    int status = FAILURE;
    static char routine[] = "ReadHistRefiRates";
    FILE  *refiFp;
    long    iNum;
    long    lineNum;
    long   numRefiRates = 0;
    TDate  closingDate = 0;
    char   refiFileName[132];
    char   inString[132];
 
    switch (amortIndexType)
    {
      case MBS_PP_REFI_TYPE_FH30:
        strcpy(refiFileName,CheckPathDelimiters(theSymbolTable.get(DAVIS_HIST_REFI30)).c_str());
        break;
      case MBS_PP_REFI_TYPE_FH15:
        strcpy(refiFileName,CheckPathDelimiters(theSymbolTable.get(DAVIS_HIST_REFI15)).c_str());
        break;
      case MBS_PP_REFI_TYPE_CMT10:
        strcpy(refiFileName,MBS_REFI_CMT10_DAT);
        break;
      default:
        GtoErrMsg("%s: Invalid refi type: %d\n", routine,amortIndexType);
        goto done;
    }
    
    /* Try to open the file */
	refiFp = fopen(refiFileName,"r");
    
    /* If file still not found, complain */
    if( refiFp IS NULL )
    {
        GtoErrMsg("%s: cannot read input file '%s' (open error)",
                  routine,refiFileName);
        goto done;
    }

    /* This is the only to count number of refi rates in the hist.
     * refi data file; using system call does the same job. We need
     * this number to allocate enough memory space for refi date,
     * rate, and point array
     */
    lineNum = 0;
    while ( fgets (inString, sizeof(inString), refiFp) != NULL )
    {
        lineNum++;
    }
    fclose(refiFp);

    refiFp = fopen(refiFileName,"r");
    if( refiFp IS NULL )
    {
        GtoErrMsg("%s: cannot read input file '%s' (open error)",
                  routine,refiFileName);
        goto done;
    }

    read_perl_long(refiFp,"closing_date",1,&closingDate);
/*
    read_perl_long(refiFp,"number",1,&numRefiRates);
*/

    if ((*amortIndexDates = NEW_ARRAY(TDate,lineNum)) IS NULL)
    {
        GtoErrMsg("%s : Failed to allocate amortIndexDates\n", routine);
        goto done;
    }
    if ((*amortIndexRates = NEW_ARRAY(double,lineNum)) IS NULL)
    {
        GtoErrMsg("%s : Failed to allocate amortIndexRates\n", routine);
        goto done;
    }
    if (amortIndexType == MBS_PP_REFI_TYPE_FH30 ||
        amortIndexType == MBS_PP_REFI_TYPE_FH15)
    {
        if ((*fhCommitPts = NEW_ARRAY(double,lineNum)) IS NULL)
        {
            GtoErrMsg("%s : Failed to allocate fhCommitPts\n", routine);
            goto done;
        }
    }

    iNum = 0;
    while ( fgets (inString, sizeof(inString), refiFp) != NULL )
    {
        if (inString[0] == '#' ||
            inString[0] == '\n')
        {
            continue;
        }
 
        if (amortIndexType == MBS_PP_REFI_TYPE_FH30 ||
            amortIndexType == MBS_PP_REFI_TYPE_FH15)
        {
            sscanf(inString,"%ld %lf %lf",
                &((*amortIndexDates)[iNum]),
                &((*amortIndexRates)[iNum]),
                &((*fhCommitPts)[iNum]));
            iNum++;
        }
        else if (amortIndexType == MBS_PP_REFI_TYPE_CMT10)
        {
            sscanf(inString,"%ld %lf",
                amortIndexDates[iNum],amortIndexRates[iNum]);
            iNum++;
        }
    }
    fclose(refiFp);
 
    (*numAmortIndexRates) = iNum;

    status = SUCCESS;
done:
    return status;
}



/*************************************************************************
 *  PRIVATE FUNCS
 *************************************************************************/


/************************************************************
 * trapGtoError()
 * traps error calls in Gto library, so that we can save the
 * most recent error string
 ************************************************************/
PRIVATE
TBoolean trapGtoError(char *errStr,
		      void *userData)
{
    /* copy Gto message to our local string variable */
    strncpy(lastGtoErrStr,errStr,ERRSTRLEN-1);
    /* resume with Gto error handling */
    return TRUE;					     
}
