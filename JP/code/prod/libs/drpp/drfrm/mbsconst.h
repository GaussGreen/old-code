/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsconst.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Davis (Shuenn-Tyan) Lee
 *			   Robert Lenk
 *			   (Derivatives Research)
 *	Code version    :  1.35
 *	Extracted	:  2/11/97 at 15:13:03
 *	Last Updated	:  2/11/97 at 15:13:01
 ***************************************************************************
 *      Common constants/structures/methods for prepay/pricing code
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __mbsconst_h
#define __mbsconst_h

#include <stdio.h>
#include <stdlib.h>

#include "timeline.h"
#include "virttree.h"
#include "modltype.h"

#include "mbsbase.h"
#include "ppconst.h"


/************************************************************
 *  CONSTANTS
 ************************************************************/
/*
 * Basic types of MBS collateral
 */
#define	MBS_MBSTYPE_FIXED	1
#define	MBS_MBSTYPE_ARM	        2

/*
 * MBS issuing agency types
 */
#define	MBS_AGENCY_FNMA	        1
#define	MBS_AGENCY_FHLMC	2
#define	MBS_AGENCY_GOLD	        3
#define	MBS_AGENCY_GNMAI	4
#define	MBS_AGENCY_GNMAII	5
#define	MBS_AGENCY_WHOLE	6

/*
 * Interest rate models available via MBS analytics
 */
#define MBS_IRMODEL_TREE1FAC   1       /* 1-factor virtual L-N tree */
#define MBS_IRMODEL_TREE2FAC   2       /* 2-factor virtual L-N tree */
#define MBS_IRMODEL_MCTREE1FAC 3       /* MC generated from older 1-fac. tree */

/*
 * Available prepay models
 */
#define	MBS_PP_MODEL_NONE	0 /* prepays = 0.0 */
#define	MBS_PP_MODEL_CONST	1 /* prepay = const rate */
#define	MBS_PP_MODEL_VECTOR	2 /* prepays from 1-d prepay vector */
#define	MBS_PP_MODEL_FIX_MGRP	3 /* F. Quinn fixed-rate MGRP model */
#define	MBS_PP_MODEL_ARM	4 /* F. Quinn ARM model */
#define	MBS_PP_MODEL_ARMSTATIC  5 /* static vector from F. Quinn ARM model */

/*
 * Possible MBS terms
 */
#define MBS_PP_MBSTERM_30       30
#define MBS_PP_MBSTERM_20       20
#define MBS_PP_MBSTERM_15       15

/*
 * Possible types of refi index rates
 */
#define MBS_PP_REFI_TYPE_NONE   0 /* i.e., for some prepay models, don't
                                   * need ANY refi rate */
#define MBS_PP_REFI_TYPE_FH30   1 /* denotes FHLMC 30yr fixed commitment */
#define MBS_PP_REFI_TYPE_FH15   2 /* denotes FHLMC 15yr fixed commitment */
#define MBS_PP_REFI_TYPE_CMT10  3 /* denotes CMT 10yr yield */
#define MBS_PP_NUM_REFI_TYPES   4

/*
 * Forms of amortization
 */
#define	MBS_PP_SPD_PSA          0 /* PSA, in decimal (e.g., 1.00) */
#define	MBS_PP_SPD_CPR          1 /* CPR, in decimal (e.g., 0.05) */
#define	MBS_PP_SPD_SMM          2 /* SMM, in decimal (e.g., 0.05) */

/* 
 * Misc
 */
#define MBS_NUM_STD_PPTWKS       4 /* number of "standard" prepay twk params */



/* 
 * Macros for convenience in option pricing 
 * (x assumed to be ptr to TMbsPricingAssump) 
 */
#define VALUE_CAP(x) (x->optionsToPrice[0] IS GTO_MODEL_CAP || \
                      x->optionsToPrice[0] IS GTO_MODEL_COLLAR)
#define VALUE_FLR(x) (x->optionsToPrice[0] IS GTO_MODEL_FLOOR || \
                      x->optionsToPrice[0] IS GTO_MODEL_COLLAR)


/************************************************************
 *  TCalibFreeFunc is new in analytic library V8.0; therefore
 *  we need to avoid this if user link to order version of
 *  analytic library.   - Davis Lee  12/18/96
 ************************************************************/
#ifdef PRIMUS_V
typedef void (TCalibFreeFunc)(void *);
#endif

/************************************************************
 *  STRUCTURES
 ************************************************************/

/*
 * TMbsIoRule  (input structure)
 * Defines IO/cpn behavior (current coupon, resets, caps/floors, etc)
 * Can be used for net cpn, gross cpn, prepay resets, etc.
 * (can also think of "IO rule" as a kind of "reset rule")
 */
typedef struct {
    /* CASHFLOW MULTIPLIER */
    double        multiplier;          /* all IO cashflows scaled by this */
    /* TYPE OF COUPON BEHAVIOR */
    TBoolean      cpnIsFloating;       /* TRUE=floating cpns, FALSE=fixed
                                        * NB: IO should be considered fixed
                                        * if the floating reset index wgt = 0 */
    long          cpnPayDC;            /* day count to use in paying cpn
                                        * (usually GTO_B30_360 or GTO_ACT_360) */
    TBoolean      ruleIsAmortRule;     /* TRUE if IO rule is amort reset info;
                                        * FALSE for normal cpn reset info */
    double        constServSpread;     /* @@ cheat for the gross coupon:
                                        * if the MBS is floating, we often
                                        * approximate curr WAC assuming a
                                        * const serv sprd: wac=currnet+sprd */
    /* FOR FIXED CPN */
    double        fixedCoupon;         /* (monthly) coupon rate, in decimal */
    /* FOR FLOATING CPN: INDEX RATE (INCL MULT/SPRD) */
    TFloatRateArray *cpnIndexRate;     /* Reset index for floating cpn;
                                        * includes spread and wgt */
    /* FOR FLOATING CPN: ROUNDING */
    TBoolean      useCpnRounding;      /* TRUE if (floating) cpn is to be
                                        * rounded */
    double        cpnRounding;         /* reset coupons are rounded to
                                        * nearest multiple of this value;
                                        * use zero for no rounding;
                                        * e.g., GNMA ARM resets involve rounding
                                        * cpns to nearest 1/8th of a point, 
                                        * hence would use 0.00125 */
    /* FLOATING CPN LOOKBACK "RULE"
     * (always have this info, even if user supplied reset date list) */
    long          effResetLkbkDays;    /* effective reset date is (approx) this
                                        * many (act) days before start of REGULAR
                                        * accrual period of new cpn */
    /* FLOATING CPN GENERAL RESET PARAMS 
     * (may all be 0, if user supplied reset date list) */
    TDate         firstEffResetDate;   /* Effective reset date of first
                                        * reset for this IO cpn;
                                        * NB: if user supplies reset dates,
                                        * this param will be 0 */
    long          resetIntervalMons;   /* # months between resets;
                                        * NB: if user supplies reset dates,
                                        * this param will be 0 */
    long          effResetDOM;         /* day-of-month of eff. reset date 
                                        * NB: if user supplies reset dates,
                                        * this param will be 0 */
    /* FLOATING CPN RESET DATES/STRIKES, KNOWN CPN */
    long          numResets;           /* # of resets in arrays below */
    TDateList    *effResetDL;          /* list of effective reset dates;
                                        * NB: all dates should be unique,
                                        * no repeated resets */
    TDateList    *firstAccStDL;        /* accrual start date of FIRST accrual
                                        * period corresponding to this reset;
                                        * NB: at present, these only used
                                        * for matching MBS accrual periods;
                                        * NB: if this IO rule describes a
                                        * deal leg, these will be deal, 
                                        * not U/L MBS, accrual pds */
    TDateList    *firstEstRegAccStDL;  /* (approx) reg accrual start date,
                                        * computed from above dates;
                                        * NB: these are the dates used to 
                                        * match deal accrual periods */
    double       *stickyCapSpreads;    /* Cpn spreads of periodic cap for each
                                        * reset, in decimal (e.g., 0.01) */
    double       *stickyFlrSpreads;    /* Cpn spreads of periodic floor for each
                                        * reset, in decimal (e.g., -0.01) */
    double       *capStrikes;          /* Normal cap strike for each reset
                                        * (decimal, e.g., 0.11)*/
    double       *flrStrikes;          /* Normal flr strike for each reset
                                        * (decimal, e.g., 0.01)*/
    TBoolean     *isKnownCpn;          /* array of booleans indicating if
                                        * reset cpn is already known */
    double       *knownCpns;           /* Array of (possibly) known coupon 
                                        * rates (monthly, decimal) at each
                                        * reset; -1.0 if not known */
    double       *knownCpnIndexRates;  /* Array of (possibly) known index rates
                                        * at each reset; -1.0 if not known;
                                        * NB: these are index rates with any
                                        * weight and spread included, to be
                                        * consistent w/defn in cpnIndexRate */
    /* For date calculations */
    TDate         todayDate;
    char          holidayFile[MAXPATHLEN]; /* holiday filename(for date arith)*/
} TMbsIoRule;




/*
 * TMbsPoRule  (input structure)
 * Defines the way the PO cashflows are determined
 */
typedef struct {
    /* CASHFLOW MULTIPLIER */
    double        multiplier;        /* all PO cashflows scaled by this */
} TMbsPoRule;



/*
 * TMbsPrepayAssump  (input structure)
 * Defines the prepay assumptions to use in computing
 * prepays or pricing a deal
 */
typedef struct {
    /* BASIC PREPAY ASSUMPTIONS */
    long           prepayModel;	       /* prepay model to use */
    TBoolean      inclSchedAmort;      /* TRUE=return prepays & sched am,
					* FALSE=return just prepays */
    TBoolean      inclBullet;	       /* TRUE=force total amort on bulletDate;
					* does so for ANY prepay model */
    TBoolean      useNewWacInfo;       /* TRUE=use new wac info.; default is
                                        * FALSE */
    TDate         bulletDate;	       /* date on which to do 100% amort */
    long           amortForm;	       /* form for output amort
					* (e.g., MBS_PP_SPD_PSA, etc.) */
    /* Assumption may changed by pricing model */
    long           numAmortMons;        /* # months for which to compute amort */
    TDate         startDate;           /* first month for which amort needed
                                        * (day-of-month ignored) */
    /* SPECIAL INFO NEEDED FOR SPECIFIC PREPAY MODELS */
    double        grossMargin;         /* gross margin for ARM */
    /* Vector prepay model */
    TCurve       *prepayVector;	       /* (for MBS_PP_MODEL_VECTOR only) 
					* monthly amortization rates (in SMM,
					* decimal), where dates should be 
					* 1st of reported prepay month */
    /* "std" tweaks for prepay models (valid w/any prepay model)
     * Important: should be MBS_NUM_STD_PPTWKS of these
     * Also: defaults (no tweaks) obtained with 1.0 for all */
    double        stdSmmSpdTwk;        /* mult. tweak to prepay (SMM) speed */
    double        stdSmmFlrTwk;        /* mult. tweak to floor spd (in SMM) */
    double        stdSteepTwk;         /* mult. tweak to rationality */
    double        stdIncentPtTwk;      /* mult. tweak to incentive measure */
    /* "FQuinn" tweaks for Mgrp and ARM */
    double        turnoverTwk;	       /* additive tweak (in CPR) to turnover 
					* rateof fully seasoned paper 
					* (positive for higher turnover prepay 
					* rate); use 0.0 for no tweak */
    double        incentPtTwk;         /* shift in prepay incentive point--the
					* prepay incentive ratio at which 
					* refi-based prepays begin to kick in 
					* is shifted deeper in-the-money 
					* (lower rates) by this amount;
					* use 0.0 for no tweak */
    double        totSteepTwk;         /* multiplicative increase in the
					* rationality of the callable region 
					* of the prepay function is increased
					* by this amount; note that this 
					* tweak applies both to burned-out 
					* (base) and hazard components of
					* prepay function; use 1.0 for no 
					* tweak */
    double        baseSteepTwk;        /* multiplicative increase in the
					* rationality of the callable region 
					* of the prepay function is increased 
					* by this amount; note that this tweak
					* ONLY applies to the burned-out 
					* portion of the prepay function;
					* use 1.0 for no tweak */
    double        ageRampTwk;
    /* "INTERMEDIATE STORAGE":
     *    Pre-computed parameters needed by tree prepay function;
     *    set by setupMbsTreeAmort() */
    double        seasonality[DEF_MBS_PP_NUM_SEASONALITY];
    /* To use with logistic2(): */
    double        mgrpScurveLogist[DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST];
    double        mgrpSeasLogist[DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST];
    /* ARM data */
    double        armAbsRateEffect;
    double        armAvgHistRefi;
    /* To use with logistic2(): */
    double        armSmmLogist[DEF_MBS_PP_ARM_NUM_SMM_LOGIST];
    double        armSeasLogist[DEF_MBS_PP_ARM_NUM_SEAS_LOGIST];
    /* To use with fast_logistic2(): */
    double        armSmmFastLogist[DEF_MBS_PP_ARM_NUM_SMM_LOGIST];
    double        armSeasFastLogist[DEF_MBS_PP_ARM_NUM_SEAS_LOGIST];
    /* precompute zero incentive logistic result */
    double        armRefiIncShift;
    double        armNewWacRefiIncShift;
    /* Intermediate value and also leave here for PRIMUS */
    double        grossCpn;            /* (I) gross coupon rate (in decimal);
                                        * should be the gross coupon accruing
                                        * in startDate */
    long          wala;                /* WALA (in months) as of start_date */
    long          warm;                /* wgted-avg time to maturity (in months)
                                        * as of start_date */
} TMbsPrepayAssump;



/*
 *  TMbsDeal  (input structure)
 *  Describes an MBS instrument or prepay-linked swap
 */
typedef struct {
    /* UNDERLYING MBS */
    long          mbsType;	       /* basic type of MBS (fixed/ARM) */
    char          mbsSubtype[LINESTRLEN]; /* supplementary info */
    long          mbsAgency;	       /* issuing agency 
					* (e.g., MBS_AGENCY_FNMA) */
    long          mbsTerm;	       /* Either MBS_PP_MBSTERM_30 or
					* MBS_PP_MBSTERM_15 */ 
    TDate         mbsOrigDate;         /* date(month) of MBS origination */
    TDate         mbsMatDate;          /* maturity date (use this to 
					* compute WARM) */
    long          mbsEffOrigTerm;      /* WARM+WALA */
    long          mbsAccrualStDom;     /* Day-of-month of start of accrual pd
                                        * (for MBS, not deal!); typically 1 */
    long          mbsCashflowDom;      /* Day-of-month of cashflows for this
                                        * type of MBS (not of deal!) */
    long          mbsPayDelayDays;     /* for MBS, not deal */
    /* UNDERLYING MBS CPN BEHAVIOR(S) 
     * Note: the accrual dates stored in these MBS/U-L IO rules 
     * are MBS, not deal, accrual start dates */
    TMbsIoRule   *mbsNetIoRule;	       /* net cpn behavior of underlying */
    TMbsIoRule   *mbsGrossIoRule;      /* gross cpn behavior of underlying */
    TMbsPoRule   *mbsPoRule;	       /* PO behavior of underlying */
    /* PREPAY INFO */
    TMbsIoRule   *mbsAmortRule;        /* amortization/prepay reset behavior
                                        * of underlying;
                                        * (@@ later, perhaps put this in 
                                        * prepay structure...) */
    TMbsPrepayAssump *mbsPrepayAssump; /* main block of prepay inputs/assump */
    double        mbsGrCcSprd[MBS_PP_NUM_REFI_TYPES];
				       /* (for this mbs type) spread 
					* between gross current coupon and 
					* each of all valid refi index 
					* rate types */
    /* BASIC DEAL TIMING INFO 
     * These variables are always defined, 
     * regardless of how cf/accrual dates determined */
    long          dealPayIntervalMons; /* length (in months) of regular deal
                                        * accrual period or, equivalently, 
                                        * months between regular cashflow dates;
                                        * typically 1 */
    long          dealMbsAccDelayDays; /* approx. # days (actual) between 
                                        * start of MBS accrual period and of 
                                        * (later) deal regular accrual period
                                        * (e.g., 19 for ARM classic swaps) */
    /* DEAL TIMING INFO FOR "AUTOMATIC" GENERATION OF DATES:
     * These variables only defined when accrual/cashflow dates 
     * are to be generated by code, not supplied by user */
    TDate         dealFirstAccStDate;  /* start date of 1st accrual
                                        * period of deal (may be stub) */
    TDate         dealRegularAccStDate; /* start date of 1st regular
                                        * accrual period of deal */
    long          dealCashflowDelayDays; /* deal cashflows occur this # days
                                        * after start of NEXT accrual pd */
    TDate         dealMatDate;	       /* date of last cashflow of deal */
    /* DEAL ACCRUAL/CASHFLOW DATES 
     * (always defined) */
    long          numCashflows;        /* # of dates in next two lists */
    TDateList    *dealAccStDL;         /* list of deal accrual start dates
                                        * (NB: real accrual pds, may be stub) */
    TDateList    *dealAccEndDL;        /* list of deal accrual end dates
                                        * (NB: real accrual pds, may be stub) */
    TDateList    *dealCfDL;            /* list of deal cashflow dates
                                        * (NB: if final acc pd is stub, final CF
                                        * will be irregular) */
    /* (APPROX) REGULAR ACC/CASHFLOW DATES
     * (simplifies some date arith; computed from above dates) */
    TDateList    *dealEstRegAccStDL;   /* approx deal reg acc start dates */
    TDateList    *dealEstRegAccEndDL;  /* approx deal reg acc end dates */
    TDateList    *dealEstRegCfDL;      /* approx deal reg cashflow dates */
    /* LEG A & B OF SWAP */
    TBoolean      doingLegAIo;         /* TRUE if leg A IO being priced */
    TBoolean      doingLegAPo;         /* TRUE if leg A PO being priced */
    TBoolean      doingLegBIo;         /* TRUE if leg B IO being priced */
    TBoolean      doingLegBPo;         /* TRUE if leg B PO being priced */
    /* Note: the accrual dates stored in these leg IO rules,
     * are deal, not MBS, accrual start dates,
     * EVEN IF THE RULE IS OTHERWISE IDENTICAL TO MBS/ASSET RULE */
    TMbsIoRule   *legAIoRule;          /* IO rule for leg A of swap */
    TMbsPoRule   *legAPoRule;          /* PO rule for leg A of swap */
    TMbsIoRule   *legBIoRule;          /* IO rule for leg B of swap */
    TMbsPoRule   *legBPoRule;          /* PO rule for leg B of swap */
} TMbsDeal;



/*
 * TMbsRateEnv  (input structure)
 * Defines the interest rate environment to use
 * in pricing the deal
 */
typedef struct {
    TDate    tradeDate;              /* date on which market information */
    TCurve  *discZeroCurve;          /* (I) discounting zero curve */
    long      numIndexZeroCurves;     /* (I) # of index zero curves */
    TCurve **indexZeroCurves;        /* (I) array of index zero curves */
    long    *numSettleDays;          /* (I) for each zero curve, settle day
                                      * convention (0th value for disc curve,
                                      * then values for index curves) */
    TCurve  *volCurve;               /* Base volatility curve; note */
    /* MARKET INFORMATION NEEDED FOR PREPAYS - USER SUPPLIED REFI RATES */
    long     amortIndexType;          /* type of index rate used to drive
                                      * prepays; e.g., MBS_PP_REFI_TYPE_FH30 */
    TDate   amortIndexStartDate;     /* month of first amort index rate */
    long     numAmortIndexRates;      /* number of amortIndexRates[] */
    double *amortIndexRates;         /* past & future amort index rates */
    double *adjAmortIndexRates;      /* points correction - past & future amort
                                      * index rates */
    double *fhCommitPts;             /* (only if amort index is FH15/30 
                                      * commitment rate) past & future
                                      * monthly-average */
    /* MARKET INFORMATION NEEDED FOR PREPAYS - HIST. REFI RATES */
    TDate   histAmortIndexStartDate; /* month of first amort index rate */
    TDate  *histAmortIndexDates;     /* month of first amort index rate */
    long     numHistAmortIndexRates;  /* number of amortIndexRates[] */
    double *histAmortIndexRates;     /* past & future amort index rates */
    double *adjHistAmortIndexRates;  /* point correction - past & future amort
                                      * index rates */
    double *histFhCommitPts;         /* (only if amort index is FH15/30 
                                      * commitment rate) past & future
                                      * monthly-average */
    /* MARKET INFORMATION FOR ARM PREPAYMENT */
    double  armVsFixedCcSprd;        /* spread between ARM gross current
                                      * coupon and the FHLMC commitment rate
                                      * (in fh_commit[]) as of the prepay
                                      * startDate;  spread should be decimal;
                                      * e.g., -0.0108);  this spread is added
                                      * to fh_commit[] rates to obtain future
                                      * values of the ARM gross current
                                      * coupon rate */
    double  adjArmVsFixedCcSprd;     /* point correction */
    TCurve *armResetIndex;           /* future reset index rates (rate from
                                      * which cpn resets are computed;
                                      * typically CMT 1yr), at yearly
                                      * intervals, starting with the next
                                      * reset AFTER start date for prepays
                                      * (start_date);  we define reset date
                                      * as the date (month, actually) in which
                                      * the new (reset) coupon starts to
                                      * accrue */
} TMbsRateEnv;
 
 
/*
 * TMbsPricingAssump  (input structure)
 * Defines the way in which the deal should be priced
 */
typedef struct {
    TDate     valueDate;        /* date for which to compute PV */
    TBoolean  detFirstCfUsingMbsAcc;  /* affects how first cashflow to price
                                       * is determined from valueDate:
                                       * TRUE: use mbs accrual periods
                                       * FALSE: use deal/swap accrual periods */
    TDate     pvBalanceDate;      /* price to be computed rel. to MBS balance
                                   * on this date */
    double    oas;                /* option-adjusted spread */
    TBoolean  computeOasFromPv;   /* FALSE ("normal") for PV from OAS;
                                   * TRUE to compute OAS from PV */
    double    pvForOas;           /* PV to use if computing OAS from PV */
    /* MANNER OF PRICING */
    TBoolean   priceBareOptions;  /* FALSE=compute deal PV w/embed. options; 
                                   * TRUE=compute PV of some/all options only*/
    char       optionsToPrice[2]; /* (when pricing bare options) 2-char string:
                                   * 1st char:
                                   *        GTO_MODEL_CAP      'C'
                                   *        GTO_MODEL_FLOOR    'F'
                                   *        GTO_MODEL_COLLAR   'L'
                                   * 2nd char: dep on context
                                   * (e.g., is 'N','S','B' in ARM model)*/
    /* TYPE OF INT RATE MODEL */
    long       irModel;         /* MBS_IRMODEL_TREE1FAC, 
                                 * MBS_IRMODEL_TREE2FAC,
                                 * MBS_IRMODEL_MCTREE1FAC */
    /* TREE MODELS */
    long       minPpy;
    long       minPpy2;
    TDate      ppySwitchDate;
    long       maxNumSvDims;    /* # values in numStateVarGridPts */
    long      *numStateVarGridPts;
                                  /* # sample pts for each sv dim */
    /* MISC */
    char       holidayFile[MAXPATHLEN]; /* holiday filename */
    char       dataFileDir[MAXPATHLEN]; /* dir for suppl. data files */
} TMbsPricingAssump;




/*
 * TMbsPricingOutput  (output structure)
 * Container for output results from pricing functions
 */
typedef struct {
    double      pv;               /* PV of deal as of value date */
    double      oasFromPv;        /* (if getting OAS from PV) OAS */
    TDate       pvBalanceDate;
    TDate       firstCashflowDate;
} TMbsPricingOutput;



/***************************************************************************
 *  PUBLIC METHODS
 ***************************************************************************/


/***************************************************************************
 * MakeMbsDeal();
 * Constructor for TMbsDeal
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Note: calling code must explicitly construct internal
 * structures (e.g., prepay, mbsio, mbspo) after calling this function
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsDeal
   (/* UNDERLYING MBS */
    long          mbsType,             /* (I) basic type of MBS (fixed/ARM) */
    char         *mbsSubtype,          /* (I) supplementary info */
    long          mbsAgency,           /* (I) issuing agency 
                                        * (e.g., MBS_AGENCY_FNMA) */
    long          mbsTerm,             /* (I) Either MBS_PP_MBSTERM_30 or
                                        * MBS_PP_MBSTERM_15 */

    TDate         mbsWarmWalaAsOfDate, /* (I) date(month) for which WARM and
                                        * WALA (below) are valid */
    long          mbsWarm,             /* (I) WARM of MBS in as-of month */
    long          mbsWala,             /* (I) WALA of MBS in as-of month */
    /* UNDERLYING MBS CPN BEHAVIOR(S) */
    TMbsIoRule   *mbsNetIoRule,	       /* (I) net cpn behavior of underlying */
    TMbsIoRule   *mbsGrossIoRule,      /* (I) gross cpn behavior of underlying */
    TMbsPoRule   *mbsPoRule,	       /* (I) PO behavior of underlying */
    /* PREPAY INFO */
    TMbsIoRule   *mbsAmortRule,        /* amortization/prepay reset behavior
                                        * of underlying */
    TMbsPrepayAssump *mbsPrepayAssump, /* main block of prepay inputs/assump */
    double       *mbsGrCcSprd,         /* array of MBS_PP_NUM_REFI_TYPES
                                        * spreads (for this mbs type) between
					* gross current coupon and each of all
					* valid refi index rate types */
    /* BASIC DEAL TIMING INFO  */
    long          dealPayIntervalMons, /* length (in months) of regular deal
                                        * accrual period or, equivalently, 
                                        * months between regular cashflow dates;
                                        * typically 1 */
    long          dealMbsAccDelayDays, /* approx. # days (act) between start of
                                        * MBS accrual period and of (later) deal
                                        * regular accrual period (e.g., 19 for
                                        * ARM classic swaps) */
    /* DEAL ACCRUAL/CASHFLOW DATES:
     * We allow two methods of defining these dates:
     *   1) "By rule": User can supply these "general" reset 
     *   parameters--constructor will then build internal lists).... */
    TDate         dealFirstAccStDate,  /* start date of 1st accrual
                                        * period of deal (may be stub);
                                        * 0 if using method 2 */
    TDate         dealRegularAccStDate, /* start date of 1st regular
                                        * accrual period of deal;
                                        * 0 if using method 2 */
    long           dealCashflowDelayDays, /* deal cashflows occur this # days
                                        * after start of NEXT accrual pd;
                                        * typ. 0 for swaps;
                                        * 0 if using method 2 */
    TDate         dealMatDate,	       /* date of last cashflow of deal;
                                        * 0 if using method 2 */
    /*   2) "By list": user can supply lists of cf/accrual dates */
    long          numCashflows,        /* # of dates in next two lists;
                                        * 0 if using method 1 */
    TDate        *dealAccStDates,      /* list of deal accrual start dates;
                                        * NULL if using method 1 */
    TDate        *dealCfDates,         /* list of deal cashflow dates;
                                        * NULL if using method 1 */
    /* LEG A & B OF SWAP */
    TMbsIoRule   *legAIoRule,          /* IO rule for leg A of swap */
    TMbsPoRule   *legAPoRule,          /* PO rule for leg A of swap */
    TMbsIoRule   *legBIoRule,          /* IO rule for leg B of swap */
    TMbsPoRule   *legBPoRule,          /* PO rule for leg B of swap */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsDeal    **mbsDeal);         /* (I/O) structure to alloc/modify */


/***************************************************************************
 * FreeMbsDeal()
 * Destructor for TMbsDeal
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsDeal
   (TMbsDeal **mbsDeal);        /* (I/O) structure to free */


/***************************************************************************
 * CheckMbsDeal()
 * Checks contents of deal struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsDeal
   (TMbsDeal *mbsDeal);     /* (I) structure to check */


/***************************************************************************
 * PrintMbsDeal()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsDeal
   (char        *label,         /* (I) label string telling us what
                                 * this deal applies to */
    long         maxNDates,     /* (I) max # date rows to list (saves space) */
    TMbsDeal *mbsDeal);         /* (I) structure contains MBS deal info. */


/***************************************************************************
 * MakeMbsIoRule();
 * Constructor for TMbsIoRule
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsIoRule
   (double        multiplier,          /* (I) all IO cashflows scaled by this */
    /* CPN BEHAVIOR */
    double        cpnIndexMult,        /* (I) multiplier of index rate for cpn;
                                        * use 0.0 for fixed-rate */
    double        cpnRateOrSprd,       /* (I) for fixed cpn: cpn rate;
                                        * for float cpn: sprd to index rate */
    long          cpnIndexFreq,        /* (I) (for flt cpn) pay freq (cpns/yr)
                                        * of index rate (e.g., 2=semiann) */
    long          cpnIndexMatMons,     /* (I) (for flt cpn) maturity (months) of
                                        * index rate (e.g., 12 for 1yr) */
    long          cpnIndexCurve,       /* (I) (for flt cpn) which zero curve 
                                        * to use for cpn (e.g., 1,2) */
    long          cpnIndexDC,          /* (I) (for flt cpn) day count of index
                                        * rate (e.g., GTO_B30_360) */
    long          cpnPayDC,            /* (I) For fixed or flt: day count used 
                                        * in paying cpn (e.g., GTO_B30_360) */
    TBoolean      ruleIsAmortRule,     /* (I) TRUE if IO rule is amort reset info;
                                        * FALSE for normal cpn reset info */
    double        constServSpread,     /* (I) @@ cheat for the gross coupon:
                                        * if the MBS is floating, we often
                                        * approximate curr WAC assuming a
                                        * const serv sprd: wac=currnet+sprd */
    /* FOR FLOATING CPN: KNOWN CPNS */
    long          nKnownCpns,          /* (I) # known coupons 
                                        * (I) (i.e., #values in arrays below) */
    TDate        *knownEffResetDates,  /* (I) Array of effective reset dates
                                        * for each known cpn;
                                        * NB: if this constructor builds list
                                        * of reset dates, these known reset
                                        * dates must be consistent with the
                                        * reset date parameters: first date,
                                        * reset interval, reset d-o-m */
    double       *knownCpns,           /* (I) Array of known coupon
                                        * rates (monthly, decimal) */
    double       *knownResetIndexRates, /* (I) Array of known index rates
                                        * index rates (as of reset date);
                                        * NB: should be raw index rate 
                                        * (e.g., CMT1yr), not including
                                        * cpnIndexMult or cpnRateOrSprd */
    /* FOR FLOATING CPN: ROUNDING */
    double        cpnRounding,         /* (I) reset coupons are rounded to
                                        * nearest multiple of this value;
                                        * use zero for no rounding;
                                        * e.g., GNMA ARM resets involve rounding
                                        * cpns to nearest 1/8th of a point, 
                                        * hence would use 0.00125 */
    /* FLOATING CPN LOOKBACK "RULE" 
     * (always supply this info, even if user supplied reset date list) */
    /* @@ Note that since we now also pass in firstAccStDates[] (below),
     * this parameter is redundant, but keep it for now */
    long          effResetLkbkDays,    /* effective reset date is (approx) this
                                        * many (act) days before start of REGULAR
                                        * accrual period of new cpn 
					* NB: if IO rule describes deal, we mean
					* delay of reset prior to DEAL accrual;
					* if rule is for MBS, then MBS accrual
                                        * NB: if user supplies reset dates
                                        * (in effResetDates[]), this delay
                                        * should be approx. consistent with 
                                        * these reset dates */
    /* FLOATING CPN RESET DATES: 
     * We allow two methods of defining reset dates:
     *   1) "By rule": User can supply these "general" reset 
     *   parameters--constructor will then build internal lists).... */
    TDate         firstEffResetDate,   /* (I) Effective reset date of first
                                        * reset for this IO cpn;
                                        * NB: if user supplies reset dates,
                                        * (meth 2) this param will be 0 */
    long          resetIntervalMons,   /* (I) # months between resets;
                                        * NB: if user supplies reset dates,
                                        * (meth 2) this param will be 0 */
    long          effResetDOM,         /* (I) day-of-month of eff. reset date 
                                        * NB: if user supplies reset dates,
                                        * (meth 2) this param will be 0 */
    /*   2) "By list": user can supply list of reset dates (effResetDates[]) */
    long          numResets,           /* (I) # of resets (& # elements
                                        * in the arrays below);
                                        * use 0 if using meth 1 */
    TDate        *effResetDates,       /* (I) list of effective reset dates;
                                        * NB: all dates should be unique
                                        * (no repeated resets);
                                        * NULL if reset dates are to generated
                                        * automatically (meth 1) */
    TDate        *firstAccStDates,     /* (I) list of first accrual start dates
                                        * corresponding to each reset above;
                                        * NB: if this IO rule describes the
                                        * underlying MBS, these will be
                                        * MBS, not deal, accrual pds */
    TBoolean      firstAccMayBeStub,   /* (I) if TRUE, firstAccStDates[0]
                                        * may be start of stub period */
    /* FLOATING CPN STRIKES (can be NULL):
     *   Depends on method of specifying reset dates:
     *   1) each pointer assumed to point at a single sprd/strike
     *   2) ptr assumed to be an array of numResets sprds/strikes */
    double       *stickyCapSpreads,    /* (I) Cpn spreads of periodic cap for 
                                        * each reset, in decimal (e.g., 0.01) */
    double       *stickyFlrSpreads,    /* (I) Cpn spreads of periodic floor for 
                                        * each reset, in decimal (e.g., -0.01) */
    double       *capStrikes,          /* (I) Normal cap strike for each reset
                                        * (decimal, e.g., 0.11)*/
    double       *flrStrikes,          /* (I) Normal flr strike for each reset
                                        * (decimal, e.g., 0.01)*/
    /* MISC */
    char         *holidayFile,         /* (I) name of GTO holiday file */
    long         *numSettleDays,       /* (I) array of settle days
                                        * for each zero curve */
    TDate         todayDate,           /* (I) (for determining truly historical
                                        * resets) current date */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsIoRule  **mbsIoRule);          /* (I/O) structure to alloc/modify */



/***************************************************************************
 * CopyMbsIoRule();
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CopyMbsIoRule
   (TMbsIoRule   *inIoRule,   /* (I) IO cashflow rule struc */
    TMbsIoRule  **outIoRule); /* (O) copy of IO cashflow rule struc */



/***************************************************************************
 * FreeMbsIoRule()
 * Destructor for TMbsIoRule 
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsIoRule
   (TMbsIoRule  **mbsIoRule);        /* (I/O) structure to free */


/***************************************************************************
 * CheckMbsIoRule()
 * Checks contents of io rule struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsIoRule
   (TMbsIoRule  *mbsIoRule);     /* (I) structure to check */


/***************************************************************************
 * PrintMbsIoRule()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsIoRule
   (char        *label,         /* (I) label string telling us what
                                 * this IO rule applies to */
    long         maxNDates,     /* (I) max # date rows to list (saves space) */
    TMbsIoRule  *mbsIoRule);    /* (I) structure contains IO rule info. */

/***************************************************************************
 * MakeMbsPoRule();
 * Constructor for TMbsPoRule 
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsPoRule
   (double        multiplier,        /* (I) all PO cashflows scaled by this */
    TMbsPoRule  **mbsPoRule);        /* (I/O) structure to alloc/modify */



/***************************************************************************
 * CopyMbsPoRule();
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CopyMbsPoRule
   (TMbsPoRule   *inPoRule,   /* (I) PO cashflow rule struc */
    TMbsPoRule  **outPoRule); /* (O) copy of PO cashflow rule struc */



/***************************************************************************
 * FreeMbsPoRule()
 * Destructor for TMbsPoRule 
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPoRule
   (TMbsPoRule  **mbsPoRule);        /* (I/O) structure to free */


/***************************************************************************
 * CheckMbsPoRule()
 * Checks contents of po rule struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsPoRule
   (TMbsPoRule  *mbsPoRule); /* (I) structure to check */


/***************************************************************************
 * PrintMbsPoRule()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPoRule
   (char        *label,         /* (I) label string telling us what
                                 * this PO rule applies to */
    TMbsPoRule  *mbsPoRule);    /* (I) structure contains PO rule info. */


/***************************************************************************
 * MakeMbsPpAssump();
 * Constructor for TMbsPrepayAssump
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsPpAssump
   (/* BASIC PREPAY ASSUMPTIONS */
    long      prepayModel,       /* (I) prepay model to use */
    TBoolean inclSchedAmort,    /* (I) TRUE=return prepays & sched am,
                                 * FALSE=return just prepays */
    TBoolean inclBullet,        /* (I) TRUE=force total amort on bulletDate;
                                 * does so for ANY prepay model */
    TBoolean useNewWacInfo,     /* (I) TRUE=use new wac info.; default is
                                 * FALSE */
    TDate    bulletDate,        /* (I) date on which to do 100% amort */
    TDate    startDate,         /* (I) first month for which amort needed
                                 * (day-of-month ignored) */
    double   grossCpn,          /* (I) gross coupon rate (in decimal); 
                                 * should be the gross coupon accruing 
                                 * in startDate */ 
    long     wala,              /* (I) WALA (in months) as of start_date */
    long     warm,              /* (I) wgted-avg time to maturity (in months)
                                 * as of start_date */
    long      amortForm,         /* (I) form for output amort
                                 * (e.g., MBS_PP_SPD_PSA, etc.) */
    long      numAmortMons,      /* (I) # months for which to compute amort */
    /* SPECIAL INFO NEEDED FOR SPECIFIC PREPAY MODELS */
    double   grossMargin,       /* (I) gross margin for ARM */
    /* Vector prepay model */
    TCurve  *prepayVector,      /* (I) (for MBS_PP_MODEL_VECTOR only) monthly
                                 * amortization rates (in SMM, decimal), where
                                 * dates should be 1st of reported prepay
                                 * month */
    /* "Std" tweaks for prepay models (valid w/any prepay model);
     * (default) no tweak obtained with 1.0 */
    double   stdSmmSpdTwk,      /* (I) mult. tweak to prepay (SMM) speed */
    double   stdSmmFlrTwk,      /* (I) mult. tweak to floor speed (in SMM)*/
    double   stdSteepTwk,       /* (I) mult. tweak to rationality */
    double   stdIncentPtTwk,    /* (I) mult. tweak to incentive measure */
    /* "FQuinn" tweaks for Mgrp and ARM */
    double   turnoverTwk,       /* (I) additive tweak (in CPR) to turnover rate
                                 * of fully seasoned paper (positive for higher
                                 * turnover prepay rate); use 0.0 for no
                                 * tweak */
    double   incentPtTwk,       /* (I) shift in prepay incentive point--the
                                 * prepay incentive ratio at which refi-based
                                 * prepays begin to kick in is shifted deeper
                                 * in-the-money (lower rates) by this amount;
                                 * use 0.0 for no tweak */
    double   totSteepTwk,       /* (I) multiplicative increase in the
                                 * rationality of the callable region of the
                                 * prepay function is increased by this amount;
                                 * note that this tweak applies both to
                                 * burned-out (base) and hazard components of
                                 * prepay function; use 1.0 for no tweak */
    double   baseSteepTwk,      /* (I) multiplicative increase in the
                                 * rationality of the callable region of the
                                 * prepay function is increased by this amount;
                                 * note that this tweak ONLY applies to the
                                 * burned-out portion of the prepay function;
                                 * use 1.0 for no tweak */
    double   ageRampTwk,
    TMbsPrepayAssump
          **mbsPrepayAssump);   /* (I/O) structure to alloc/modify */


/***************************************************************************
 * CopyMbsPpAssump();
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CopyMbsPpAssump
   (TMbsPrepayAssump  *inPpAssump,   /* (I) input prepay assump struc */
    TMbsPrepayAssump **outPpAssump); /* (O) copy of prepay assump struc */


/***************************************************************************
 * FreeMbsPrepayAssump()
 * Destructor for TMbsPoRule 
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPrepayAssump
   (TMbsPrepayAssump **mbsPrepayAssump);        /* (I/O) structure to free */


/***************************************************************************
 * CheckMbsPrepayAssump()
 * Checks contents of pp assump struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsPrepayAssump
   (TMbsPrepayAssump *mbsPrepayAssump); /* (I) structure to check */


/***************************************************************************
 * PrintMbsPrepayAssump()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPpAssump
   (TBoolean     printPpModelInfo, /* (I) TRUE to also output info specific to
                                    * prepay-model code */
    TMbsPrepayAssump
          *mbsPrepayAssump);   /* (I) structure contains prepay assump. */


/***************************************************************************
 * RefiRateNeeded()
 * From MBS and prepay info, determines type of refi rate needed
 ***************************************************************************/
int EXPORT RefiRateNeeded
   (TMbsDeal         *mbsDeal,           /* (I) deal info */
    TMbsPrepayAssump *mbsPrepayAssump,   /* (I) prepay info */
    long              *refiRateType);     /* (O) type of refi rate needed for
                                          * this MBS and this prepay model
                                          * (e.g., MBS_PP_REFI_TYPE_FH30) */


/***************************************************************************
 * MakeMbsRateEnv();
 * Constructor for TMbsRateEnv
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsRateEnv
   (TDate    tradeDate,          /* (I) date on which market information */
    TCurve  *discZeroCurve,      /* (I) discounting zero curve */
    long      numIndexZeroCurves, /* (I) # of index zero curves */
    TCurve **indexZeroCurves,    /* (I) array of index zero curves */
    long    *numSettleDays,      /* (I) for each zero curve, settle day
                                  * convention (0th value for disc curve,
                                  * then values for index curves) */
    TCurve  *volCurve,           /* (I) Base volatility curve; note */
    TDate   *armResetDates,      /* (I) array of dates for arm reset curve */
    double  *armResetRates,      /* (I) array of rates for dates of curve */
    long      numArmResetPts,     /* (I) # elements in dates & rates array */
    double   armResetFreq,       /* (I) basis. see GtoRateToDiscount */
    long     dayCountConv,       /* (I) see GtoDayCountFraction */
    /* MARKET INFORMATION NEEDED FOR PREPAYS - USER SUPPLIED REFI RATES */
    long     amortIndexType,      /* (I) type of index rate used to drive
                                  * prepays; e.g., MBS_PP_REFI_TYPE_FH30 */
    TDate   amortIndexStartDate, /* (I) month of first amort index rate */
    long     numAmortIndexRates,  /* (I) number of amortIndexRates[] */
    double *amortIndexRates,     /* (I) past & future amort index rates */
    double *fhCommitPts,         /* (I) (only if amort index is FH15/30
                                  * commitment rate) past & future
                                  * monthly-average */
    double  armVsFixedCcSprd,    /* (I) spread between ARM gross current coupon
                                  * and the FHLMC commitment rate (in
                                  * fh_commit[]) as of the prepay startDate;
                                  * spread should be decimal;  e.g., -0.0108);
                                  * this spread is added to fh_commit[] rates
                                  * to obtain future values of the ARM gross
                                  * current coupon rate */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsRateEnv **mbsRateEnv);   /* (I/O) structure to alloc/modify */


/***************************************************************************
 * FreeMbsRateEnv()
 * Destructor for TMbsRateEnv
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsRateEnv
   (TMbsRateEnv **mbsRateEnv); /* (I/O) structure to free */


/***************************************************************************
 * CheckMbsRateEnv()
 * Checks contents of rate env. struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsRateEnv
   (TMbsRateEnv *mbsRateEnv);  /* (I) structure to check */


/***************************************************************************
 * PrintMbsRateEnv()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsRateEnv
   (TBoolean     printPrepayInfo, /* (I) TRUE to also output info specific to
                                   * prepay-model code */
    TMbsRateEnv *mbsRateEnv);    /* (I) structure contains MBS rate Env. info. */


/***************************************************************************
 * MakeMbsPricingAssump();
 * Constructor for TMbsPricingAssump
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsPricingAssump
   (TDate           valueDate,  /* (I) date for which to compute PV */
    TBoolean        detFirstCfUsingMbsAcc, 
                                    /* (I) affects how first cashflow to price
                                     * is determined from valueDate:
                                     * TRUE: use mbs accrual periods
                                     * FALSE: use deal/swap accrual periods */
    TDate           pvBalanceDate, /* (I) price to be computed rel. to MBS 
                                    * balance on this date */
    double          oas,            /* (I) option-adjusted spread */
    TBoolean        computeOasFromPv, /* FALSE ("normal") for PV from OAS;
                                       * TRUE to compute OAS from PV */
    double          pvForOas,       /* PV to use if computing OAS from PV */
    /* MANNER OF PRICING */
    TBoolean        priceBareOptions, /* (I) FALSE=compute deal PV 
                                       * w/embed. options; TRUE=compute PV 
                                       * of some/all options only*/
    char           *optionsToPrice,   /* (I) (when pricing bare options) 
                                       * 2-char string:
                                       * 1st char:
                                       *        GTO_MODEL_CAP      'C'
                                       *        GTO_MODEL_FLOOR    'F'
                                       *        GTO_MODEL_COLLAR   'L'
                                       * 2nd char: dep on context
                                       * (e.g., is 'N','S','B' in ARM model)*/
    /* TYPE OF MODEL */
    long            irModel,        /* (I) int rate model to use 
                                     * (e.g., MBS_IRMODEL_TREE1FAC) */
    /* TREE MODELS */
    long            minPpy,         /* (I) pds per year (early) */
    long            minPpy2,        /* (I) later pds per year */
    TDate           ppySwitchDate,  /* (I) early/later transition date */
    long            maxNumSvDims,   /* (I) # values in numStateVarGridPts */
    long           *numStateVarGridPts, /* (I) # sample pts for each sv dim */
    /* MISC */
    char           *holidayFile,    /* (I) for holiday/date calc */
    char           *dataFileDir,    /* (I) dir for suppl. data files */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsPricingAssump
                  **mbsPricingAssump);/* (I/O) structure to alloc/modify */


/***************************************************************************
 * FreeMbsPricingAssump()
 * Destructor for TMbsPricingAssump
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPricingAssump
   (TMbsPricingAssump **mbsPricingAssump); /* (I/O) structure to free */


/***************************************************************************
 * CheckMbsPricingAssump()
 * Checks contents of pricing assump. struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
    /* @# update this, since we now have explicit scalar params */
CheckMbsPricingAssump
   (TMbsPricingAssump *mbsPricingAssump); /* (I) structure to check */


/***************************************************************************
 * PrintMbsPricingAssump()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
    /* @# update this, since we now have explicit scalar params */
PrintMbsPricingAssump
   (TMbsPricingAssump *mbsPricingAssump); /* (I) structure contains princing
                                           * assumptions  */


/***************************************************************************
 * MakeMbsPricingOutput
 * Constructor for pricing output struc 
 * Note: allocates only; does not fill member data
 ***************************************************************************/
int EXPORT
MakeMbsPricingOutput
   (TMbsPricingOutput    **mbsPricingOutput); /* (I/O) struc to alloc/mod */

/***************************************************************************
 * FreeMbsPricingOutput()
 * Destructor for TMbsPricingOutput
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPricingOutput
   (TMbsPricingOutput **mbsPricingOutput); /* (I/O) structure to free */


#endif
