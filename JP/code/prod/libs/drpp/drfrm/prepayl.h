/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  prepayl.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Forrest Quinn
 *			   Davis (Shuenn-Tyan) Lee
 *			   Robert Lenk
 *			   (Derivatives Research)
 *	Code version    :  1.11
 *	Extracted	:  3/10/97 at 13:04:22
 *	Last Updated	:  3/10/97 at 13:04:17
 ***************************************************************************
 *      Wrapper function for prepay models
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __prepayl_h
#define __prepayl_h

#include "cdate.h"
#include "bastypes.h"
#include "mbsconst.h"

#define MBS_PP_MGRP_VERSION  "MGRP " __MBSMODEL_VER__
#define MBS_PP_FQARM_VERSION "FQ ARM " __MBSMODEL_VER__

/***************************************************************************
 *  PUBLIC FUNCS
 ***************************************************************************/

/***************************************************************************
 *  MbsAmortL
 *  Prepay model wrapper function for PRIMUS, KAPITAL, and Excel 
 *  and Applix spreadsheet add-ins. This wrapper function will compute
 *  amortization vectors using either the FQuinn Multi-group (fixed) or 
 *  ARM prepay models.
 *    FLAT FILES: this module depends on finding several flat files
 *  in the FIMS $DATA/LATEST data directory. Currently, DR (Lee/Lenk) 
 *  updates these files nightly. They are:
 *       mbsder_refi_cmt10.dat, mbsder_refi_fh15.dat, mbsder_refi_fh30.dat
 *    NOTE: any symbolic constants are defined in ppconst.h
 *    NOTE: as per GTO convention, wrapper inputs are pointers 
 *  to arrays of at least two elements; the first element is
 *  always the number of subsequent elements in the array;
 *  e.g., to input two doubles (e.g., 3.3, 4.4), the wrapper array 
 *  should be {2., 3.3, 4.4}
 *    NOTE: since wrapper functions are limited in the number
 *  of inputs, we group together scalar inputs into vectors of
 *  the same type. The components of these "misc scalar" inputs
 *  are detailed below:
 ==========================================================================
 * prepayInfo contains:
 * 0th (# of params):    9 (# of values)
 * 1st (prepayModel):    prepay model to use (e.g., * MBS_PP_MODEL_ARM, etc).
 * 2nd (mbsType):        basic type of MBS (e.g., MBS_MBSTYPE_FIXED or
 *                       MBS_MBSTYPE_ARM)
 * 3rd (mbsAgency);      issuing agency (e.g., * MBS_AGENCY_FNMA, ...)
 * 4th (amortForm):      form for output amort--one of: MBS_PP_SPD_CPR,
 *                       MBS_PP_SPD_SMM, MBS_PP_SPD_PSA
 * 5th (inclSchedAmort): boolean: if true, inc sched. amort.; if false,
 *                       exclude scheduled amort
 * 6th (numAmortMons):   # months for which to compute amort
 * 7th (mbsTerm):        Either MBS_PP_MBSTERM_30 or MBS_PP_MBSTERM_15;
 *                       determines whether model will treat collateral as 30-
 *                       or 15-year product
 * 8th (wala):           WALA (in months) as of startDate
 * 9th (warm):           wgted-avg time to maturity (in months) as of startDate
 ==========================================================================
 * ppTweakInfo contains:
 * 0th (# of params):    4 (# of values)
 * 1st (turnoverTwk):    multiplicative tweak (in CPR) to turnover rate of fully
 *                       seasoned paper (positive for higher turnover prepay
 *                       rate); use 0.0 for no tweak.
 * 2nd (incentPtTwk):    shift in prepay incentive point--the prepay incentive
 *                       ratio at which refi-based prepays begin to kick in is
 *                       shifted deeper in-the-money (lower rates) by this
 *                       amount;  use 0.0 for no tweak.
 * 3rd (totSteepTwk):    multiplicative increase in the rationality of the
 *                       callable region of the prepay function is increased
 *                       by this amount;  note that this tweak applies both
 *                       to burned-out (base) and hazard components of prepay
 *                       function; use 1.0 for no tweak.
 * 4th (baseSteepTwk):   multiplicative increase in the rationality of the
 *                       callable region of the prepay function is increased
 *                       by this amount; note that this tweak ONLY applies to
 *                       the burned-out portion of the prepay function;  use
 *                       1.0 for no tweak.
 ==========================================================================
 * ppArmInfo contains:
 * 0th (# of params):    7 (# of values)
 * 1st (capSpread):      Cpn spread for periodic cap, in decimal (e.g., 0.01)
 * 2nd (floorSpread):    Cpn spread for periodic floor, in decimal(e.g., -0.01)
 * 3rd (lifeCap):        lifetime max on net cpn, in decimal (e.g., 0.12)
 * 4th (lifeFloor):      lifetime min on net cpn, in decimal (e.g., 0.02)
 * 5th (grossMargin):    gross margin for ARM
 * 6th (netMargin):      net margin for ARM
 * 7th (armVsFixedCcSprd): spread between ARM gross current coupon and the
 *                         amortIndexRates[] rate; spread should be decimal;
 *                         e.g., -0.0108);  this spread is added to
 *                         fh_commit[] rates to obtain future values of the
 *                         ARM gross current coupon rate.
 * 8th (useNewWacInfo)   : default is FALSE, otherwise use new wac info.
 * 9th (grossLifeCap)    : gross life cap (assumed the model already takes
 *                         current net & gross coupon, and both net and
 *                         gross margins.
 * 10th (grossLifeFloor) : gross life floor
 ***************************************************************************/

int EXPORT MbsAmortL2
   (TDate  *startDate,             /* (I) first month for which amort needed
                                    * (day-of-month ignored) */
    long   *prepayInfo,            /* (I) several long scalars--see doc. above */
    double *grossCpn,              /* (I) gross coupon rate of MBS, in decimal;
                                    * should that accruing in startDate */
    double *netCpn,                /* (I) net coupon rate of MBS, in decimal;
                                    * should that accruing in startDate */
    double *ppTweakInfo,           /* (I) several twk scalars--see doc. above */
    long   *amortIndexType,        /* (I) type of FHLMC commitment rate
                                    * in amortIndexRates[]; should be either 
                                    * MBS_PP_REFI_TYPE_FH30 or ..._FH15; 
                                    * NOTE: should be consistent with term
                                    * of MBS (30 or 15 yr) */
    TDate  *amortIndexStartDate,   /* (I) month of first amortIndexRates[] rate
                                    * (day-of-month ignored) */
    double *amortIndexRates,       /* (I) past & future monthly-average
                                    * commitment rates (in decimal) of FIXED
                                    * FHLMC MBS OF THE SAME ORIGINAL TERM
                                    * (e.g., for mbsTerm=30yr, use FHLMC-30y
                                    * rates; for 15yr, use FHLMC-15y rates);
                                    * Should begin at least 6mons before
                                    * startDate (though user can supply
                                    * his/her own historical rates to override 
                                    * flat-file values used by this module */
    double *fhCommitPts,           /* (I) past & future monthly-average points,
                                    * corresponding to the commitment rates
                                    * in amortIndexRates[]; in decimal (e.g., 2
                                    * points would be 0.02); same date range
                                    * as amortIndexRates[] */
    TDate  *armResetDates,         /* (I) array of dates for arm reset curve
                                    * NB: These should be "MBS-convention" 
                                    * ARM reset dates (e.g., 4/1/97 for
                                    * "April reset" ARM--that is, the date on
                                    * which accrual of the newly-reset cpn
                                    * begins, not the actual index lookup date*/
    double *armResetRates,         /* (I) array of YEARLY reset rates
                                    * NB: These should be the ARM index rates
                                    * as of the lookup date */
    double *ppArmInfo,             /* (I) several ARM scalars--see doc. above */
    /* Following inputs alter default behavior of both ARM & MGRP models: */
    double *seasonality,           /* (I) seasonality multiplier (12) */
    /* Following inputs alter default behavior of ARM model only: */
    long   *maxLag,                /* (I) maxLag; (months) */ 
    double *absoluteRateEffect,    /* (I) absoluteRateEffect */
    double *histRefiIdx,           /* (I) array of 3: {FH30, FH15, CMT10} */
    double *lagWgts,               /* (I) lag wgts (3) */
    double *armBaseLogistic,       /* (I) armBaseLogistic (4) */
    double *armSeasLogistic,       /* (I) armSeasLogistic (4) */
    /* Following inputs alter default behavior of MGRP model only: */
    long   *numGroups,             /* (I) # of groups */
    double *shiftSmmArray,         /* (I) 1st: shift refinance inc.
                                    *     2nd: Additional smm */
    double *scurveLogistic,        /* (I) scurveLogistic (6) */ 
    double *seasLogistic,          /* (I) seasLogistic (4) */ 
    double *groupPs,               /* (I) groupPs (6) */
    double *ltAges,                /* (I) ltAges (6) */
    double *mgrpLagWgts,           /* (I) mgrpLagWgts (6) */
    /* OUTPUTS */
    double *amort);                /* (O) array of monthly amort. rate (in
                                    * decimal; form determined by amortForm);
                                    * note that calling code must allocate
                                    * this array with at least 1+numAmortMons
                                    * elements */

#endif
