/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  prepay.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Forrest Quinn
 *			   Davis (Shuenn-Tyan) Lee
 *			   Robert Lenk
 *			   (Derivatives Research)
 *	Code version    :  1.11
 *	Extracted	:  5/19/97 at 13:54:17
 *	Last Updated	:  5/19/97 at 13:54:10
 ***************************************************************************
 *      Public functions of MBS prepay module
 *      Module version:  1.1.11
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __prepay_h
#define __prepay_h

#include "cdate.h"
#include "bastypes.h"

#include "ppconst.h"
#include "mbsconst.h"

/***************************************************************************
 *  DEFAULT PREPAY MODEL CONSTANTS
 ***************************************************************************/

/* These flat files contain historical monthly-average
 * values of the indicated refi index rate. 
 * They are updated nightly by the 
 * /home/mbs/cronjobs/bin/MbsDataGrab.csh script, run from the 
 * cron table of user "mbs", and copied to both the $DRDATA
 * and FIMS $DATA/LATEST data areas. */
/* 
 * 2/97: To allow for future PC compatibility, code now checks
 * first for new DOS-compat. refi name, then (on UNIX only) for older name;
 * eventually, we'll phase out the old filename.
 */
#define MBS_REFI_FH30_DAT    "refifh30.dat"
#define MBS_REFI_FH15_DAT    "refifh15.dat"
#define MBS_REFI_CMT10_DAT   "reficm10.dat"

#ifdef UNIX
#define MBS_REFI_FH30_DAT_OLD    "mbsder_refi_fh30.dat"
#define MBS_REFI_FH15_DAT_OLD    "mbsder_refi_fh15.dat"
#define MBS_REFI_CMT10_DAT_OLD   "mbsder_refi_cmt10.dat"
#endif


/***************************************************************************
 *  PUBLIC FUNCS
 ***************************************************************************/


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
 *  mbsAgency - issuing agency (e.g., MBS_AGENCY_FNMA, ...). 
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
/* AMORTIZATION TO COMPUTE */
   (TMbsDeal     *mbsDeal,                  /* (I) Includes deal information */
    TMbsRateEnv  *mbsRateEnv,               /* (I) Includes rate information */
    TMbsPrepayDefaults  *mbsPrepayDefaults, /* (I) Includes default numbers */
/* OUTPUTS */
    double *amort);     /* (O) array of monthly prepays (in decimal, in form
                         * determined by amortForm);  note that calling code
                         * must have allocated this array, with at least
                         * numAmortMons rate pts */

int EXPORT InitHistAmortIndexRates
   (TMbsRateEnv *mbsRateEnv);

int EXPORT InitHistAmortAvgIndexRates
   (TMbsRateEnv  *mbsRateEnv,
    TMbsDeal     *mbsDeal);

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
    double **fhCommitPts);       /* (O) hist. monthly-average points,
                                  * corresponding to the commitment rates in
                                  * amortIndexRates[]; in decimal (e.g., 2
                                  * points would be 0.02); same date range
                                  * as amortIndexRates[] */

/***************************************************************************
 * InitMbsTreeAmort
 * Init function to be called (once) before 
 * using "tree" prepay function MbsTreeAmort() (pre-computes some 
 * parameters, thus speeding up calculation)
 ***************************************************************************/
int EXPORT InitMbsTreeAmort
   (TMbsDeal *mbsDeal,             /* (I) deal to be priced */
   long    refiIndexRateType,      /* (I) type of refi index rate
                                    * e.g., MBS_PP_REFI_TYPE_CMT10 */
    double grCpnSprdOffRefiIx);    /* (I) spread (decimal) to add to refi idx*/


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
    double     *smm);           /* (O) scheduled amort., in decimal SMM */


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
    double *prepay_over_period); /* (O) prepay rate over period */


#endif
