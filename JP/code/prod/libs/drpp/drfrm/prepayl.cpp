/***************************************************************************
 *      SCCS Keyword Information
 *      ------------------------
 *      Module name     :  prepayl.c
 *      Company Name    :  JP Morgan Securities Inc.
 *      Author          :  Davis (Shuenn-Tyan) Lee
 *                         Derivatives Research
 *      Code version    :  1.34
 *      Extracted       :  3/11/97 at 12:21:00
 *      Last Updated    :  3/11/97 at 12:20:59
 ***************************************************************************
 *      Wrapper function for prepay models
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#include <stdio.h>
extern "C" {
#include "cdate.h"
#include "ldate.h"
#include "cerror.h"
#include "dateconv.h"
#include "convert.h"
#include "mdydate.h"
#include "tcurve.h"
#include "datelist.h"
}

#include "mbsbase.h"
#include "ppconst.h"
#include "mbsconst.h"
#include "ppbase.h"
#include "prepay.h"


/***************************************************************************
 *  MbsAmortL2
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
 * 0th (# of params):    5 (# of values)
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
 * 5th (ageRampTwk) :    It is to adjust age ramp to get rid of speed pick
 *                       up from new loan if necessary. Davis Lee 2/10/2000
 ==========================================================================
 * ppArmInfo contains:
 * 0th (# of params):    10 (# of values)
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
int EXPORT
MbsAmortL2
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
    double *amort)                 /* (O) array of monthly amort. rate (in
                                    * decimal; form determined by amortForm);
                                    * note that calling code must allocate
                                    * this array with at least 1+numAmortMons
                                    * elements */
{
    F_INIT("MbsAmortL2");
    long      prepayModel;  /* prepay model to use (e.g., MBS_PP_MODEL_ARM,
                            * etc). */
    long      mbsType;      /* basic type of MBS (e.g., MBS_MBSTYPE_FIXED or
                            * MBS_MBSTYPE_ARM) */
    long      mbsAgency;    /* issuing agency (e.g., MBS_AGENCY_FNMA, ...) */
    long      amortForm;    /* form for output amort--one of: MBS_PP_SPD_CPR,
                            * MBS_PP_SPD_SMM, MBS_PP_SPD_PSA */
    long      inclSchedAmort;/* boolean: if true, include sched. amort.; if
                             * false, exclude scheduled amort */
    long      numAmortMons; /* # months for which to compute amort */
    long      mbsTerm;      /* Either MBS_PP_MBSTERM_30 or MBS_PP_MBSTERM_15;
                            * determines whether model
                            * will treat collateral as 30- or 15-year product */
    long     wala;         /* WALA (in months) as of startDate */
    long     warm;         /* wgted-avg time to maturity (in months) as of
                            * startDate */
    double   armResetFreq = 2;
    double   capSpread;    /* Cpn spread for periodic cap, in decimal */
    double   floorSpread;  /* Cpn spread for periodic floor, in decimal */
    double   lifeCap;      /* lifetime max on net cpn, in decimal */
    double   lifeFloor;    /* lifetime min on net cpn, in decimal */
    double   grossMargin;  /* gross margin for ARM */
    double   netMargin;    /* net margin for ARM */
    double   armVsFixedCcSprd;
    double   fUseNewWacInfo;
    double   grossLifeCap;
    double   grossLifeFloor;
    double   turnoverTwk;
    double   incentPtTwk;
    double   totSteepTwk;
    double   baseSteepTwk;
    double   ageRampTwk;
    double   mbsServFee;
    TBoolean useNewWacInfo;
    long     dayCountConv = GTO_ACT_365;   /* see GtoDayCountFraction */
    long      idx;
    long      i;
    long     numSettleDays[2] = {0,0};
    double grCcSprd[MBS_PP_NUM_REFI_TYPES]; /* (for this mbs type) spread
                                             * between gross current coupon
                                             * and each of all valid refi
                                             * index rate types */
    TMbsDeal     *mbsDeal = NULL;         
    TMbsPrepayAssump *mbsPrepayAssump = NULL; 
    TMbsIoRule   *grossIoRule = NULL;
    TMbsIoRule   *netIoRule = NULL;
    TMbsRateEnv  *mbsRateEnv = NULL; 
    TMbsPrepayDefaults *mbsPrepayDefaults = NULL;


    /* init mbs base lib--see explanation of inputs */
    MbsBaseInit("JPM_PREPAY",
                FALSE,  /* not print msgs to stdout (i.e., no printf()) */
                TRUE,   /* (i.e., GtoErrMsgOn()) turn on error logging */
                "");    /* ...but use default error log file name */

    for (i = 0; i < MBS_PP_NUM_REFI_TYPES; i++)
    {
        grCcSprd[i] = 0.;
    }

    if ( (long)prepayInfo[0] ISNT 9)
    {
        GtoErrMsg("%s: Invalid # [%d] scalars in prepayInfo input vector\n",
            routine,(long)prepayInfo[0]);
        goto done;
    }
    idx = 1;
    prepayModel = (long)prepayInfo[idx++];
    mbsType = (long)prepayInfo[idx++];
    mbsAgency = (long)prepayInfo[idx++];
    amortForm = (long)prepayInfo[idx++];
    inclSchedAmort = (long)prepayInfo[idx++];
    numAmortMons = (long)prepayInfo[idx++];
    mbsTerm = (long)prepayInfo[idx++];
    wala = (long)prepayInfo[idx++];
    warm = (long)prepayInfo[idx++];

    if ( (long) ppTweakInfo[0] ISNT 5)
    {
        GtoErrMsg("%s: Invalid # [%d] scalars in ppTweakInfo input vector\n",
            routine,(long)ppTweakInfo[0]);
        goto done;
    }
    idx = 1;
    turnoverTwk = ppTweakInfo[idx++];
    incentPtTwk = ppTweakInfo[idx++];
    totSteepTwk = ppTweakInfo[idx++];
    baseSteepTwk = ppTweakInfo[idx++];
    ageRampTwk = ppTweakInfo[idx++];

    if ( (long) ppArmInfo[0] ISNT 10)
    {
        GtoErrMsg("%s: Invalid # [%d] scalars in ppArmInfo input vector\n",
            routine,(long)ppArmInfo[0]);
        goto done;
    }
    idx = 1;
    capSpread = ppArmInfo[idx++];
    floorSpread = ppArmInfo[idx++];
    lifeCap = ppArmInfo[idx++];
    lifeFloor = ppArmInfo[idx++];
    grossMargin = ppArmInfo[idx++];
    netMargin = ppArmInfo[idx++];
    armVsFixedCcSprd = ppArmInfo[idx++];
    fUseNewWacInfo = ppArmInfo[idx++];
    grossLifeCap = ppArmInfo[idx++];
    grossLifeFloor = ppArmInfo[idx++];

    if (IS_ALMOST_ZERO(fUseNewWacInfo))
    {
        useNewWacInfo = FALSE;
    }
    else
    {
        useNewWacInfo = TRUE;
    }


    /* 
     * Build prepay assump struc 
     */
    if (MakeMbsPpAssump
        (prepayModel,
        inclSchedAmort,
        (TBoolean)NULL, /* (I) inclBullet */
        useNewWacInfo,  /* (I) use new wac info. */ 
        (TDate)0,       /* (I) bulletDate */
        startDate[1],
        grossCpn[1],
        wala,
        warm,
        amortForm,
        numAmortMons,
        grossMargin,
        (TCurve *)NULL, /* (I) prepayVector */
        (double)1,      /* (I) stdSmmSpdTwk; mult. tweak to prepay (SMM)
                         * speed */
        (double)1,      /* (I) stdSmmFlrTwk; mult. tweak to floor speed (SMM) */
        (double)1,      /* (I) stdSteepTwk; mult. tweak to rationality */
        (double)1,      /* (I) stdIncentPtTwk; mult. tweak to incentive
                         * measure */
        turnoverTwk,
        incentPtTwk,
        totSteepTwk,
        baseSteepTwk,
        ageRampTwk,
         &mbsPrepayAssump) IS FAILURE)
    {
        GtoErrMsg("%s: ERROR returned in MakeMbsPpAssump\n", routine);
        goto done;
    }
#ifdef PP_DEBUG
    PrintMbsPpAssump(TRUE,mbsPrepayAssump);
#endif


    /* 
     * Build IO rule, depending on whether this is ARM or fixed 
     */
    if( netCpn IS NULL )
    {
        mbsServFee = 0.0050;    /* dummy value */
    }
    else
    {
        mbsServFee = grossCpn[1] - netCpn[1];
    }
    if(mbsType IS MBS_MBSTYPE_FIXED)
    {
        if (MakeMbsIoRule(1.,                /* (I) multiplier 1=on */
                          0.,                /* (I) 0=fixed */
                          grossCpn[1],       /* (I) cpnSpread (cpn) */
                          0,                 /* (I) cpnIndexFreq */
                          0,                 /* (I) cpnIndexMatMons */
                          0,                 /* (I) cpnIndexCurve */
                          0,                 /* (I) cpnIndexDC */
                          0,                 /* (I) cpnPayDC */
                          FALSE,             /* (I) ruleIsAmortRule */
                          mbsServFee,        /* (*) curr serv spread */
                          /* no known cpns/info */
                          0,                 /* (I) nKnownCpns */
                          (TDate *)NULL,     /* (I) knownEffResetDates */
                          (double *)NULL,    /* (I) knownCpns */
                          (double *)NULL,    /* (I) knownResetIndexRates */
                          /* no rounding */
                          0.,                /* (I) cpnRounding */
                          0,                 /* (I) effResetLkbkDays */
                          /* reset timing info (not relevant for fixed) */
                          0,                 /* (I) firstEffResetDate */
                          0,                 /* (I) resetIntervalMons */
                          0,                 /* (I) effResetDOM */
                          /* explicit date lists (not used) */
                          0,                 /* numResets */
                          (TDate *)NULL,     /* (I) effResetDates */
                          (TDate *)NULL,     /* (I) 1st acc st dates */
                          FALSE,             /* 1st accrual may be stub */
                          /* caps/floors (not relevant for fixed) */
                          NULL,              /* (I) Cpn spreads of pd. cap for
                                              * each reset */
                          NULL,              /* (I) Cpn spreads of periodic floor
                                              * for each reset */
                          NULL,              /* (I) Normal cap strike for each
                                              * reset (decimal, e.g., 0.11) */
                          NULL,              /* (I) Normal flr strike for each
                                              * reset (decimal, e.g., 0.01) */
                          "NONE",            /* (I) name of GTO holiday file */
                          (long *)NULL,      /* (I) numSettleDays */
                          /* set "today" to be in distant past, to make 
                           * all resets look like future resets */
                          (TDate) 35000,     /* (I) todayDate */
                          &grossIoRule) IS FAILURE)
        {
            GtoErrMsg("%s: Failed to build IO rule for fixed-rate MBS\n", 
                      routine);
            goto done;
        }
    }
    else  /* ARM */
    {
        if (MakeMbsIoRule(1.,                /* (I) multiplier 1=on */
                          1.,                /* (I) cpnIndexMult */
                          (useNewWacInfo ? grossMargin : netMargin + mbsServFee),
                          /* (I) cpnSpread */
                          2,                 /* (I) cpnIndexFreq */
                          12,                /* (I) cpnIndexMatMons */
                          0,                 /* (I) cpnIndexCurve */
                          GTO_B30_360,       /* (I) cpnIndexDC */
                          GTO_B30_360,       /* (I) cpnPayDC */
                          FALSE,             /* (I) ruleIsAmortRule */
                          mbsServFee,        /* (*) curr serv spread */
                          /* no known cpns/info */
                          0,                 /* (I) nKnownCpns */
                          (TDate *)NULL,     /* (I) knownEffResetDates */
                          (double *)NULL,    /* (I) knownCpns */
                          (double *)NULL,    /* (I) knownResetIndexRates */
                          /* no rounding */
                          0.,                /* (I) cpnRounding */
                          0,                 /* (I) effResetLkbkDays */
                          /* reset timing info (doesn't matter for pp model,
                           * though must specify reset date far into future
                           * to keep model from insisting that this is
                           * historical (known) reset) */
                          170000,             /* (I) firstEffResetDate */
                          12,                /* (I) resetIntervalMons */
                          1,                 /* (I) effResetDOM */
                          /* explicit date lists (not used) */
                          0,                 /* numResets */
                          (TDate *)NULL,     /* (I) effResetDates */
                          (TDate *)NULL,     /* (I) 1st acc st dates */
                          FALSE,             /* 1st accrual may be stub */
                          /* caps/floors */
                          &capSpread,        /* (I) Cpn spreads of periodic cap for
                                              * each reset */
                          &floorSpread,      /* (I) Cpn spreads of periodic floor
                                              * for each reset */
                          (useNewWacInfo ? &grossLifeCap : &lifeCap),
                          /* (I) Normal cap strike for each
                           * reset (decimal, e.g., 0.11) */
                          (useNewWacInfo ? &grossLifeFloor : &lifeFloor),
                          /* (I) Normal flr strike for each
                           * reset (decimal, e.g., 0.01) */
                          
                          "NONE",            /* (I) name of GTO holiday file */
                          numSettleDays, /* (I) numSettleDays */
                          /* set "today" to be in distant past, to make 
                           * all resets look like future resets */
                          (TDate) 35000,     /* (I) todayDate */
                          &grossIoRule) IS FAILURE)
        {
            GtoErrMsg("%s: ERROR returned in MakeMbsIoRule\n", routine);
            goto done;
        }
    } /* end of ARM IO rule construction */
#ifdef PP_DEBUG
    PrintMbsIoRule("IO Rule for MBS Gross Coupon",
                   10000,       /* print all date rows */
                   mbsDeal->mbsGrossIoRule);
#endif


    /*
     * To be thorough, create net cpn IO rule by
     * copying gross cpn IO rule 
     */
    XONFAIL( CopyMbsIoRule(grossIoRule, &netIoRule) );


    /*
     * Build deal using lower structures just created
     */
    if (MakeMbsDeal
        (mbsType,
        "",                       /* (I) mbsSubtype */
        mbsAgency,
        mbsTerm,
        startDate[1],             /* warm/wala are valid as-of this date */
        warm,
        wala,
         netIoRule,               /* (I) mbsNetIoRule;  Net cpn behavior
                                   * of underlying */
         grossIoRule,             /* (I) mbsGrossIoRule;  Gross cpn behavior
                                   * of underlying */
        (TMbsPoRule *)NULL,       /* (I) mbsPoRule; PO behavior of underlying */
        (TMbsIoRule *)NULL,       /* (I) mbsAmortRule;  Amortization/prepay
                                   * reset behavior of underlying */
         mbsPrepayAssump,         /* (I) mbsPrepayAssump;  Main block of
                                   * prepay inputs/assump */
        grCcSprd,
        (long)1,                  /* (I) dealPayIntervalMons */
        (long)0,                  /* (I) dealMbsAccDelayDays */
        (TDate)0,                 /* (I) dealFirstAccStDate */
        (TDate)0,                 /* (I) dealRegularAccStDate */
        (long)0,                   /* (I) dealCashflowDelayDays */
        (TDate)0,                 /* (I) dealMatDate */
        (long)0,                  /* (I) numCashflows */
        (TDate *)NULL,            /* (I) dealAccStDates */
        (TDate *)NULL,            /* (I) dealCfDates */
        (TMbsIoRule *)NULL,       /* (I) legAIoRule;  for leg A of swap */
        (TMbsPoRule *)NULL,       /* (I) legAPoRule;  for leg A of swap */
        (TMbsIoRule *)NULL,       /* (I) legBIoRule;  for leg B of swap */
        (TMbsPoRule *)NULL,       /* (I) legBPoRule;  for leg B of swap */
       &mbsDeal) IS FAILURE)
    {
        GtoErrMsg("%s: ERROR returned in MakeMbsDeal\n",routine);
        goto done;
    }
#ifdef PP_DEBUG
    PrintMbsDeal("deal",
                 10000,       /* print all date rows */
                 mbsDeal);
#endif

    if (MakeMbsRateEnv(startDate[1],
                       (TCurve *)NULL,  /* (I) discZeroCurve */
                       (long)0,          /* (I) numIndexZeroCurves */
                       (TCurve **)NULL, /* (I) indexZeroCurves */
                       (long *)NULL,     /* (I) numSettleDays */
                       (TCurve *)NULL,  /* (I) volCurve */
                      &armResetDates[1],
                      &armResetRates[1],
                       armResetDates ? (long)armResetDates[0] : 0,
                       armResetFreq,
                       dayCountConv,
                       (long)amortIndexType[1],
                       amortIndexStartDate[1],
                       (long)amortIndexRates[0], // refi rates
                      &amortIndexRates[1],
                      &fhCommitPts[1],
                       armVsFixedCcSprd,
                      &mbsRateEnv) IS FAILURE)
    {
        GtoErrMsg("%s: ERROR returned in MakeMbsRateEnv\n", routine);
        goto done;
    }
#ifdef PP_DEBUG
    PrintMbsRateEnv(TRUE,mbsRateEnv);
#endif

    if (InitMbsPrepayDefaults((mbsDeal->mbsPrepayAssump)->prepayModel,
                              mbsDeal->mbsAgency,
                              mbsRateEnv->amortIndexType,
                              /* Following inputs are user supplied to alter
                               * global default numbers;  for both ARM model 
                               * and MGRP FIXED model */
                              seasonality ? &seasonality[1] : NULL,
                                               /* (I) seasonality (12) */
                              /* Following inputs are user supplied to alter
                               * global default numbers;  for ARM model only */
                              maxLag ? (long)maxLag[1] : 0,
                                               /* (I) maxLag; (months) */ 
                              absoluteRateEffect ? absoluteRateEffect[1] :
                                           0., /* (I) absoluteRateEffect */
                              histRefiIdx ? histRefiIdx[1] : 0.,
                                               /* (I) histrefi30Idx; (FH30) */
                              histRefiIdx ? histRefiIdx[2] : 0.,
                                               /* (I) histrefi15Idx; (FH15) */
                              histRefiIdx ? histRefiIdx[3] : 0.,
                                               /* (I) histCmt10Idx; (CMT10) */
                              lagWgts ? &lagWgts[1] : NULL,
                                               /* (I) new array of lag wgt */
                              armBaseLogistic ? &armBaseLogistic[1] : NULL,
                                               /* (I) armBaseLogistic (4) */
                              armSeasLogistic ? &armSeasLogistic[1] : NULL,
                                               /* (I) armSeasLogistic (4) */
                              /* Following inputs are user supplied to alter
                               * global default numbers; for MGRP FIXED model
                               * only */
                              numGroups ? (long)numGroups[1] : 0,
                                               /* (I) the number of groups */
                              shiftSmmArray ? shiftSmmArray[1] : 0.,
                                               /* (I) shift refinance inc. */
                              shiftSmmArray ? shiftSmmArray[2] : 0.,
                                               /* (I) Additional smm */
                              scurveLogistic ? &scurveLogistic[1] : NULL,
                                               /* (I) scurveLogistic */ 
                              seasLogistic ? &seasLogistic[1] : NULL,
                                               /* (I) seasLogistic */ 
                              groupPs ? &groupPs[1] : NULL, /* (I) groupPs */
                              ltAges ? &ltAges[1] : NULL,   /* (I) ltAges */    
                              mgrpLagWgts ? &mgrpLagWgts[1] : NULL,
                                               /* (I) mgrpLagWgts */
                             &mbsPrepayDefaults) IS FAILURE)
    {
        GtoErrMsg("%s: ERROR returned in InitMbsPrepayDefaults\n", routine);
        goto done;
    }
#ifdef PP_DEBUG
    PrintMbsPrepayDefaults((mbsDeal->mbsPrepayAssump)->prepayModel,
                           mbsPrepayDefaults);
#endif

    if ( MbsAmort2(mbsDeal,
                  mbsRateEnv,
                  mbsPrepayDefaults,
                 &amort[1]) IS FAILURE )
    {
        GtoErrMsg("%s: ERROR returned in MbsAmort\n",routine);
        goto done;
    }


    status = SUCCESS;
done:
    FreeMbsDeal(&mbsDeal);
    FreeMbsPrepayAssump(&mbsPrepayAssump);
    FreeMbsIoRule(&grossIoRule);
    FreeMbsIoRule(&netIoRule);
    FreeMbsRateEnv(&mbsRateEnv);
    FreeMbsPrepayDefaults(&mbsPrepayDefaults);

    if ( status IS FAILURE )
    {
        GtoErrMsg("%s : ARM/FRM Model Failed\n",routine);
    }
    return(status);
}
