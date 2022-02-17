/***************************************************************************
 *      SCCS Keyword Information
 *      ------------------------
 *      Module name     :  ppconst.h
 *      Company Name    :  JP Morgan Securities Inc.
 *      Author          :  Robert Lenk
 *                         Davis (Shuenn-Tyan) Lee
 *                         Derivatives Research
 *      Code version    :  1.7
 *      Extracted       :  5/19/97 at 13:54:17
 *      Last Updated    :  5/19/97 at 13:54:08
 ***************************************************************************
 *      Older constructors/destructors used by Primus interface
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __ppconst_h
#define __ppconst_h

#ifdef __cplusplus
extern "C" {
#endif

#include "tcurve.h"
#define	VERSION_3
/***************************************************************************
 *  DEFAULT PREPAY MODEL CONSTANTS
 ***************************************************************************/

/*=========================================================================
 * GLOBAL DEFAULT NUMBERS FOR ARM & MGRP FIXED MODEL
 *========================================================================*/
#define DEF_MBS_PP_NUM_SEASONALITY        12

/* In our models, we generally adjust input commitment rates to be 2-point rates
 */
#define DEF_MBS_PP_DESIRED_REFI_POINTS  0.02

/*=========================================================================
 * GLOBAL DEFAULT NUMBERS FOR ARM MODEL
 *========================================================================*/
/* minimum number of historical refi rates needed(user) */
#define DEF_MBS_PP_ARM_MIN_NUM_REFIS      6     

/* max len (months) of internal prepay arrays */
#define DEF_MBS_PP_ARM_MAX_PREPAYS        401   

/* Magnitude of the "absolute rate effect": 0 means no effect, 1 means that
 * an absolute shift of one basis point has the same impact as a relative
 * shift of 1 basis point. (Orig value=.33). */
#define DEF_MBS_PP_ARM_ABS_RATE_EFF       0.33

/* This constant (originally interpreted as the historical average of the
 * spread of ARM gross current coupon off CMT 10yr) was used as an additive
 * shift in the refi incentive measure in this model; even though the refi
 * index was then changed to be FH 15/30yr commitment rates, this shift must
 * remain the same */
#define DEF_MBS_PP_ARM_INC_SHIFT          0.0033
#define DEF_MBS_PP_ARM_NEW_WAC_INC_SHIFT  0.0010

/* historical avg levels */
#define DEF_MBS_PP_ARM_HIST_FH30_IDX      0.0608
#define DEF_MBS_PP_ARM_HIST_FH15_IDX      0.0713
#define DEF_MBS_PP_ARM_HIST_CMT10_IDX     0.0625

/* Logistic curve params for base prepay func */
#define DEF_MBS_PP_ARM_NUM_SMM_LOGIST     4
#define DEF_MBS_PP_ARM_SMM_LOGIST_RIGHT   5.25
#define DEF_MBS_PP_ARM_SMM_LOGIST_LEFT    37.
#define DEF_MBS_PP_ARM_SMM_LOGIST_WIDTH   110.
#define DEF_MBS_PP_ARM_SMM_LOGIST_INFLEC  -215.

/* Logistic curve params for seasoning ramp */
#define DEF_MBS_PP_ARM_NUM_SEAS_LOGIST    4
#define DEF_MBS_PP_ARM_SEAS_LOGIST_RIGHT  36.
#define DEF_MBS_PP_ARM_SEAS_LOGIST_LEFT   3.
#define DEF_MBS_PP_ARM_SEAS_LOGIST_WIDTH  37.
#define DEF_MBS_PP_ARM_SEAS_LOGIST_INFLEC -93.

/* max # lag months */
#define DEF_MBS_PP_ARM_MAX_NUM_LAGMONS    12	   

/* Maximum lag in months (typically we use 3) */
#define DEF_MBS_PP_ARM_MAX_LAG            3
/* Array of wgts to use w/prev refi rates */
static double DEF_MBS_PP_ARM_LAGWGTS[] = {0.25,0.50,0.25};

/* Seasonal multiplier for turnover component */
static double DEF_MBS_PP_ARM_SEASONALITY[] = {0.941, 0.826, 0.728, 0.865, \
                                              1.011, 1.124, 1.178, 1.196, \
                                              1.183, 1.023, 1.008, 0.917};

/* no longer use these historical avg. spreads, but kept in case needed in
 * future: */
/* arm gr cc off CMT10yr */
#define MBS_PP_ARM_DEF_HISTGR_ARM30_CMT10_SPRD   0.0033 
/* arm gr cc off FH30 commit */
#define MBS_PP_ARM_DEF_HISTGR_ARM30_FH30_SPRD   -0.0103 

/*=========================================================================
 * GLOBAL DEFAULT NUMBERS FOR MGRP FIXED MODEL
 *========================================================================*/
/* shift refinance incentive according to CRITRAT */
#define DEF_MBS_PP_MGRP_CRITRAT            1.1

#define DEF_MBS_PP_MGRP_NUM_GROUPS         99 
#define DEF_MBS_PP_MGRP_MAX_PREPAYS        401  /* max len (months) of 
                                                 * internal prepay arrays */

#define DEF_MBS_PP_MGRP_GNMA_GPADDLSMM     0.1463 
#define DEF_MBS_PP_MGRP_FNMA_GPADDLSMM     0.1463 
static double DEF_MBS_PP_MGRP_GNMA_SEASONALITY[]={0.941, 0.826, 0.728, 0.865, \
                                                  1.011, 1.124, 1.178, 1.196, \
                                                  1.183, 1.023, 1.008, 0.917};
static double DEF_MBS_PP_MGRP_FNMA_SEASONALITY[]={0.941, 0.826, 0.728, 0.865, \
                                                  1.011, 1.124, 1.178, 1.196, \
                                                  1.183, 1.023, 1.008, 0.917};

/* Logistic params: max,min,width,inflec,skew,disc decay */
#ifdef VERSION_3
#define DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST      9
static double DEF_MBS_PP_MGRP_GNMA_SCURVE_LOGIST[] = {0.1960, 0.0099, 0.0114, \
                                                      0.9874, 1.6528, 0.9900, \
                                                      1.1678, 0.0211, 7.2665};
static double DEF_MBS_PP_MGRP_FNMA_SCURVE_LOGIST[] = {0.1880, 0.0011, 0.0171, \
                                                      1.1466, 0.0322, 0.9900, \
                                                      1.1872, 0.0608, 7.9249};
#define DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST        4
/*
static double DEF_MBS_PP_MGRP_GNMA_SEAS_LOGIST[] = {35.34, 35.33, 0.5, 1.0};
static double DEF_MBS_PP_MGRP_FNMA_SEAS_LOGIST[] = {44.30, 44.20, 0.5, 1.0};
*/
static double DEF_MBS_PP_MGRP_GNMA_SEAS_LOGIST[] = {35.33, 3., 0.02, 1.13};
static double DEF_MBS_PP_MGRP_FNMA_SEAS_LOGIST[] = {44.32, 3., 0.02, 1.13};

#define DEF_MBS_PP_MGRP_NUM_GROUPPS             6
static double DEF_MBS_PP_MGRP_GNMA_GROUPPS[] = {0.0001, 1.2830, 0.1550, \
                                                1.1360, 1.4390, 0.1406};
static double DEF_MBS_PP_MGRP_FNMA_GROUPPS[] = {0.0001, 1.2150, 0.1000, \
                                                1.4770, 1.1000, 0.1243};
/* Longer term aging parameters */
#define DEF_MBS_PP_MGRP_NUM_LTAGES              6 

static double DEF_MBS_PP_MGRP_GNMA_LT_AGES[] = {0.06, 1.0625, 0.085, \
                                                0.0100, 0.0100, 0.0938};
static double DEF_MBS_PP_MGRP_FNMA_LT_AGES[] = {0.06, 1.0625, 0.085, \
                                                0.1340, 0.99, 0.0521};
#else
#define DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST      6
static double DEF_MBS_PP_MGRP_GNMA_SCURVE_LOGIST[] = {0.2050, 0.0600, 0.0600, \
                                                      1.1200, 2.2400, 0.7500};
static double DEF_MBS_PP_MGRP_FNMA_SCURVE_LOGIST[] = {0.2050, 0.0650, 0.0600, \
                                                      1.1300, 3.0000, 0.7500};
#define DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST        4
static double DEF_MBS_PP_MGRP_GNMA_SEAS_LOGIST[] = {75, 33., 0.05, 0.84};
static double DEF_MBS_PP_MGRP_FNMA_SEAS_LOGIST[] = {75, 33., 0.05, 0.84};
#define DEF_MBS_PP_MGRP_NUM_GROUPPS             6
static double DEF_MBS_PP_MGRP_GNMA_GROUPPS[] = {0.4912, 1.2970, 0.1000, \
                                                1.0990, 1.3660, 0.2421};
static double DEF_MBS_PP_MGRP_FNMA_GROUPPS[] = {0.0000, 1.2180, 0.1000, \
                                                1.3640, 1.5740, 0.1249};
#define DEF_MBS_PP_MGRP_NUM_LTAGES              6 
#endif

/*************************************************************************
 *  @# For now, put some hard code numbers here and will move to ppconst.h
 *  after the new code is fully tested and accepted. - 12/4/96 Davis Lee
 *************************************************************************/
#define DEF_MBS_PP_MGRP_NUM_WTLAGS             6 
static double DEF_MBS_PP_MGRP_GNMA_WTLAGS[] = {0.0, 0.0787, 0.4136, 0.2337, \
                                               0.1370, 0.1370};
static double DEF_MBS_PP_MGRP_FNMA_WTLAGS[] = {0.0, 0.0787, 0.4136, 0.2337, \
                                               0.1370, 0.1370};


/* Pre-compute yearly average service fee and average gross margin spread;
 * for ARM prepayment model.
 */

static long avgYear[] =            {1984, 1985, 1986, 1987, 1988, 1989, 1990,
                                   1991, 1992, 1993, 1994, 1995, 1996};
static double avgServFee[] =      {51.5, 55.8, 63.5,101.0, 55.1,123.1, 56.7,
                                   56.4, 54.8, 60.0, 59.3, 78.9, 75.0};
static double avgGrMarginSprd[] = {75.0, 75.4, 62.3, 57.7, 56.2, 60.0, 58.1,
                                   57.2, 56.6, 56.9, 57.7, 87.3,110.8};


typedef struct {
    long    maxLag;
    double smmLogistRight;
    double smmLogistLeft;
    double smmLogistWidth;
    double smmLogistInflec;
    double seasLogistRight;
    double seasLogistLeft;
    double seasLogistWidth;
    double seasLogistInflec;
    double histFh30Idx;
    double histFh15Idx;
    double histCmt10Idx;
    double refiIncShift;
    double absRateEffect;
    double lagWgts[DEF_MBS_PP_ARM_MAX_NUM_LAGMONS];
    double seasonality[12];   /* @@ Use constant defined symbol - DL */
} MBS_PP_ARM_DEFAULT;

typedef struct {
    long    numGroups;
    double gpAddlSmm;
    double critRat;
    double seasonality[12];   /* @@ Use constant defined symbol - DL */
    double scurveLogist[DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST];
    double seasLogist[DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST];
    double groupPs[DEF_MBS_PP_MGRP_NUM_GROUPPS];
    double wtLags[DEF_MBS_PP_MGRP_NUM_WTLAGS];
} MBS_PP_MGRP_DEFAULT;

typedef struct {
    /* For ARM & MGRP fixed models */
    double desiredRefiPoints;
    double ioMultiplier;
    double seasonality[DEF_MBS_PP_NUM_SEASONALITY];
    /* For ARM model */
    long    maxLag;
    double smmLogistRight;
    double smmLogistLeft;
    double smmLogistWidth;
    double smmLogistInflec;
    double seasLogistRight;
    double seasLogistLeft;
    double seasLogistWidth;
    double seasLogistInflec;
    double histFh30Idx;
    double histFh15Idx;
    double histCmt10Idx;
    double refiIncShift;
    double newWacRefiIncShift;
    double absoluteRateEffect;
    double armLagWgts[DEF_MBS_PP_ARM_MAX_NUM_LAGMONS];
    /* For MGRP fixed model */
    long    numGroups;
    double gpAddlSmm;
    double critRat;
    double scurveLogist[DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST];
    double seasLogist[DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST];
    double groupPs[DEF_MBS_PP_MGRP_NUM_GROUPPS];
    double ltAges[DEF_MBS_PP_MGRP_NUM_LTAGES];
    double mgrpLagWgts[DEF_MBS_PP_MGRP_NUM_WTLAGS];
} TMbsPrepayDefaults;

/*
 * Historical monthly-average FHLMC commitment rates
 * for 15- and 30-year collateral
 */
#define	hist_conv30_commit_rates_start		19630715
#define	num_hist_conv30_commit_rates		396
/* last date is 199606 */
#define hist_conv30_io_multiple                 4.375
static	double	hist_conv30_commit_rates[num_hist_conv30_commit_rates] =
{0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310,
 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07310, 0.07430, 0.07530,
 0.07600, 0.07700, 0.07690, 0.07630, 0.07550, 0.07480, 0.07440, 0.07320,
 0.07290, 0.07290, 0.07370, 0.07370, 0.07400, 0.07400, 0.07420, 0.07420,
 0.07430, 0.07440, 0.07440, 0.07440, 0.07460, 0.07540, 0.07650, 0.07730,
 0.08050, 0.08500, 0.08820, 0.08770, 0.08580, 0.08540, 0.08540, 0.08460,
 0.08410, 0.08580, 0.08970, 0.09090, 0.09280, 0.09590, 0.09980, 0.09980,
 0.09790, 0.09620, 0.09430, 0.09100, 0.08890, 0.08820, 0.08910, 0.08890,
 0.08890, 0.08940, 0.09120, 0.09220, 0.09150, 0.09100, 0.09020, 0.08810,
 0.08760, 0.08730, 0.08760, 0.08850, 0.08930, 0.09000, 0.08980, 0.08920,
 0.08810, 0.08790, 0.08720, 0.08670, 0.08690, 0.08750, 0.08830, 0.08860,
 0.08940, 0.08940, 0.08900, 0.08920, 0.08920, 0.08960, 0.09010, 0.09140,
 0.09200, 0.09350, 0.09570, 0.09710, 0.09740, 0.09780, 0.09760, 0.09860,
 0.10110, 0.10350, 0.10390, 0.10410, 0.10430, 0.10500, 0.10690, 0.11040,
 0.11090, 0.11090, 0.11300, 0.11640, 0.12830, 0.12900, 0.12880, 0.13040,
 0.15280, 0.16320, 0.14260, 0.12710, 0.12190, 0.12560, 0.13190, 0.13790,
 0.14200, 0.14790, 0.14900, 0.15130, 0.15400, 0.15580, 0.16400, 0.16700,
 0.16830, 0.17280, 0.18160, 0.18450, 0.17820, 0.16950, 0.17480, 0.17600,
 0.17160, 0.16890, 0.16680, 0.16700, 0.16820, 0.16270, 0.15430, 0.14610,
 0.13820, 0.13620, 0.13250, 0.13040, 0.12800, 0.12780, 0.12630, 0.12870,
 0.13430, 0.13810, 0.13730, 0.13540, 0.13440, 0.13420, 0.13370, 0.13230,
 0.13390, 0.13650, 0.13940, 0.14420, 0.14670, 0.14470, 0.14350, 0.14130,
 0.13640, 0.13180, 0.13080, 0.12920, 0.13170, 0.13200, 0.12910, 0.12220,
 0.12030, 0.12190, 0.12190, 0.12140, 0.11780, 0.11260, 0.10890, 0.10710,
 0.10080, 0.09940, 0.10150, 0.10690, 0.10510, 0.10200, 0.10010, 0.09980,
 0.09700, 0.09320, 0.09200, 0.09080, 0.09040, 0.09830, 0.10600, 0.10540,
 0.10280, 0.10330, 0.10890, 0.11260, 0.10650, 0.10640, 0.10380, 0.09890,
 0.09930, 0.10200, 0.10460, 0.10460, 0.10430, 0.10600, 0.10480, 0.10300,
 0.10270, 0.10610, 0.10730, 0.10650, 0.11030, 0.11050, 0.10770, 0.10200,
 0.09880, 0.09990, 0.10130, 0.09950, 0.09770, 0.09740, 0.09900, 0.10200,
 0.10270, 0.10370, 0.10480, 0.10160, 0.10040, 0.10100, 0.10180, 0.10170,
 0.10010, 0.09670, 0.09640, 0.09370, 0.09500, 0.09500, 0.09470, 0.09620,
 0.09580, 0.09240, 0.09010, 0.08860, 0.08710, 0.08500, 0.08430, 0.08760,
 0.08940, 0.08850, 0.08670, 0.08510, 0.08130, 0.07980, 0.07920, 0.08090,
 0.08310, 0.08220, 0.08000, 0.07680, 0.07500, 0.07470, 0.07460, 0.07430,
 0.07210, 0.07080, 0.06915, 0.06840, 0.07110, 0.07172, 0.07060, 0.07153,
 0.07750, 0.08420, 0.08580, 0.08430, 0.08625, 0.08512, 0.08530, 0.08925,
 0.09050, 0.09202, 0.09150, 0.08848, 0.08484, 0.08336, 0.07989, 0.07573,
 0.07586, 0.07845, 0.07629, 0.07492, 0.07380, 0.07204, 0.07027, 0.07099,
 0.07625, 0.07886, 0.08071, 0.08280};

static	double	hist_conv30_adjust_point[num_hist_conv30_commit_rates] =
{0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009,
 0.009, 0.009, 0.009, 0.009, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.012, 0.012, 0.012, 0.012,
 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.011, 0.011,
 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011,
 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012,
 0.012, 0.012, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011,
 0.011, 0.011, 0.011, 0.011, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013,
 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.016, 0.016, 0.016, 0.016,
 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016,
 0.020, 0.019, 0.019, 0.018, 0.018, 0.017, 0.017, 0.017, 0.017, 0.017,
 0.020, 0.020, 0.020, 0.020, 0.021, 0.021, 0.021, 0.021, 0.021, 0.023,
 0.021, 0.021, 0.022, 0.022, 0.022, 0.023, 0.023, 0.022, 0.022, 0.023,
 0.023, 0.022, 0.022, 0.022, 0.022, 0.020, 0.022, 0.021, 0.021, 0.021,
 0.022, 0.022, 0.022, 0.021, 0.021, 0.022, 0.023, 0.024, 0.024, 0.024,
 0.025, 0.025, 0.026, 0.026, 0.026, 0.026, 0.025, 0.025, 0.025, 0.024,
 0.026, 0.026, 0.025, 0.025, 0.026, 0.026, 0.026, 0.025, 0.024, 0.023,
 0.023, 0.023, 0.023, 0.022, 0.023, 0.023, 0.022, 0.021, 0.022, 0.021,
 0.020, 0.021, 0.022, 0.021, 0.021, 0.023, 0.023, 0.022, 0.022, 0.021,
 0.022, 0.022, 0.021, 0.021, 0.020, 0.021, 0.020, 0.021, 0.021, 0.020,
 0.020, 0.022, 0.021, 0.019, 0.021, 0.021, 0.021, 0.022, 0.022, 0.022,
 0.021, 0.021, 0.021, 0.021, 0.020, 0.020, 0.020, 0.020, 0.021, 0.021,
 0.021, 0.021, 0.020, 0.020, 0.020, 0.020, 0.021, 0.022, 0.021, 0.019,
 0.021, 0.020, 0.021, 0.020, 0.020, 0.021, 0.020, 0.019, 0.019, 0.019,
 0.018, 0.018, 0.018, 0.018, 0.019, 0.017, 0.017, 0.017, 0.016, 0.017,
 0.017, 0.018, 0.019, 0.016, 0.016, 0.015, 0.016, 0.017, 0.018, 0.016,
 0.016, 0.015, 0.015, 0.015, 0.016, 0.017, 0.017, 0.018, 0.017, 0.018,
 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.019,
 0.018, 0.019, 0.018, 0.018, 0.018, 0.018, 0.018, 0.019, 0.018, 0.018,
 0.018, 0.017, 0.018, 0.018, 0.018, 0.018};

#define	hist_conv15_commit_rates_start		19790115
#define	num_hist_conv15_commit_rates		210
/* last date is 199606 */
/* @@ For now, use the 30-year multiple: */
#define hist_conv15_io_multiple                 4.375
/* @@ later, update these with the FHLMC data (hardcopy) */
static	double	hist_conv15_commit_rates[num_hist_conv15_commit_rates] =
{0.09890, 0.09910, 0.09930, 0.10000, 0.10190, 0.10540, 0.10590, 0.10590,
 0.10800, 0.11140, 0.12330, 0.12400, 0.12380, 0.12540, 0.14780, 0.15820,
 0.13760, 0.12210, 0.11690, 0.12060, 0.12690, 0.13290, 0.13700, 0.14290,
 0.14400, 0.14630, 0.14900, 0.15080, 0.15900, 0.16200, 0.16330, 0.16780,
 0.17660, 0.17950, 0.17320, 0.16450, 0.16980, 0.17100, 0.16660, 0.16390,
 0.16180, 0.16200, 0.16320, 0.15770, 0.14930, 0.14110, 0.13320, 0.13120,
 0.12750, 0.12540, 0.12300, 0.12280, 0.12130, 0.12370, 0.12930, 0.13310,
 0.13230, 0.13040, 0.12940, 0.12920, 0.12870, 0.12730, 0.12890, 0.13150,
 0.13440, 0.13920, 0.14170, 0.13970, 0.13850, 0.13630, 0.13140, 0.12680,
 0.12580, 0.12420, 0.12670, 0.12700, 0.12410, 0.11720, 0.11530, 0.11690,
 0.11690, 0.11640, 0.11280, 0.10760, 0.10390, 0.10210, 0.09580, 0.09440,
 0.09650, 0.10190, 0.10010, 0.09700, 0.09510, 0.09480, 0.09200, 0.08820,
 0.08700, 0.08580, 0.08540, 0.09330, 0.10100, 0.10040, 0.09780, 0.09830,
 0.10390, 0.10760, 0.10150, 0.10140, 0.09880, 0.09390, 0.09430, 0.09700,
 0.09960, 0.09960, 0.09930, 0.10100, 0.09980, 0.09800, 0.09770, 0.10110,
 0.10230, 0.10150, 0.10530, 0.10550, 0.10270, 0.09700, 0.09380, 0.09490,
 0.09630, 0.09450, 0.09270, 0.09240, 0.09400, 0.09700, 0.09770, 0.09870,
 0.09980, 0.09660, 0.09540, 0.09600, 0.09680, 0.09670, 0.09510, 0.09170,
 0.09140, 0.08870, 0.09000, 0.09000, 0.08970, 0.09120, 0.09080, 0.08740,
 0.08510, 0.08360, 0.08210, 0.08000, 0.07930, 0.08260, 0.08440, 0.08350,
 0.08170, 0.08010, 0.07630, 0.07480, 0.07420, 0.07590, 0.07810, 0.07720,
 0.07500, 0.07180, 0.07000, 0.06970, 0.06960, 0.06930, 0.06710, 0.06580,
 0.06415, 0.06340, 0.06610, 0.06672, 0.06560, 0.06653, 0.07250, 0.07920,
 0.08080, 0.07930, 0.08125, 0.08012, 0.08030, 0.08425, 0.08550, 0.08702,
 0.08650, 0.08348, 0.07984, 0.07836, 0.07489, 0.07073, 0.07086, 0.07345,
 0.07129, 0.06992, 0.06880, 0.06704, 0.06527, 0.06599, 0.07125, 0.07386,
 0.07571, 0.07780};

/* @@ For now, pretend 15yr rates are already 2-pt rates */
static	double	hist_conv15_adjust_point[num_hist_conv15_commit_rates] =
{0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020};


/***************************************************************************
 *  STRUCTURES FOR MODEL-SPECIFIC OUTPUT PARAMETERS FOR TREE PRICING
 ***************************************************************************/
typedef struct {
    long     prepayModel;
    long    useSchAmort;
    double  origTermInMths;
    double  grCpnSprdOffRefiIx;
    double  seasonality[12];
    /* To use with logistic2(): */
    double  mgrpScurveLogist[DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST];
    double  mgrpSeasLogist[DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST];
    /* ARM data */
    double  armAbsRateEffect;
    double  armAvgHistRefi;
    /* To use with logistic2(): */
    double  armSmmLogist[DEF_MBS_PP_ARM_NUM_SMM_LOGIST];
    double  armSeasLogist[DEF_MBS_PP_ARM_NUM_SEAS_LOGIST];
    /* To use with fast_logistic2(): */
    double  armSmmFastLogist[DEF_MBS_PP_ARM_NUM_SMM_LOGIST];
    double  armSeasFastLogist[DEF_MBS_PP_ARM_NUM_SEAS_LOGIST];
    /* precompute zero incentive logistic result */
    double  armRefiIncShift;
    double  armNewWacRefiIncShift;
} MBS_PP_TREE_PARAMS;


/***************************************************************************
 *  STRUCTURES FOR MODEL-SPECIFIC INPUT DATA
 ***************************************************************************/

/* 
 * MBS_PP_MODELINFO_CONST:
 * model-specific inputs for constant-prepay model 
 */
typedef struct {
    long    prepayModel;	        /* always set to MBS_PP_MODEL_CONST */
    double const_prepay;	/* const. prepay rate */
    long    const_prepay_form;	/* form of above */
} MBS_PP_MODELINFO_CONST;

/* 
 * MBS_PP_MODELINFO_VECTOR:
 * model-specific inputs for vector-prepay model 
 */
typedef struct {
    long      prepayModel;         /* always set to MBS_PP_MODEL_VECTOR */
    long      vector_prepay_form;  /* form of prepay speeds in vector table */
    long      vector_num_ages;     /* # age categories (rows) in prepay vector
                                   * table */
    long      vector_num_ratemoves;/* # rate move columns in prepay vector
                                   * table */
    long     *vector_ages;         /* ages (in months) of rows of vector table */
    double  *vector_ratemoves;    /* rate move sizes (decimal; e.g., +100bp
                                   * as 0.01) */
    double **vector_prepays;      /* table of prepay speeds, indexed as:
                                   * vector_prepays[i_age][j_ratemove]  */
    long      vector_prepay_lag;   /* # months lag between rate changes and
                                   * effect on prepays */
} MBS_PP_MODELINFO_VECTOR;


int EXPORT
InitMbsPrepayDefaults
   (long prepayModel,
    long mbsAgency,
    long amortIndexType,
    /* Following inputs are user supplied to alter global default numbers;
     * for both ARM model and MGRP FIXED model
     */
    double *seasonality,        /* (I) ptr to new seasonality array
                                 * (12 values) */
    /* Following inputs are user supplied to alter global default numbers;
     * for ARM model only
     */
    long     maxLag,             /* (I) new value for max lag (months) */ 
    double  absoluteRateEffect, /* (I) absolute rate effect coeff. */
    double  histrefi30Idx,      /* (I) new value for historic value of
                                 * refi_index rate (FH30 commit rate) */
    double  histrefi15Idx,      /* (I) new value for historic value of
                                 * refi_index rate (FH15 commit rate) */
    double  histCmt10Idx,       /* (I) new value for historic value of
                                   refi_index rate (CMT10 rate) */
    double *armLagWgts,         /* (I) ptr to new array of lag wgts (expected
                                 * to contain max_lag values) */
    double *armBaseLogistic,    /* (I) ptr to new array of logistic params
                                 * (currently 4 values) for base prepay
                                 * function logistic */
    double *armSeasLogistic,    /* (I) ptr to new array of logistic params 
                                 * (currently 4 values) for seasoning ramp 
                                 * logistic */
    /* Following inputs are user supplied to alter global default numbers;
     * for MGRP FIXED model only 
     */
    long    numGroups,           /* (I) the number of groups used in
                                 * multi-group */
    double critRat,             /* (I) shift refinance incentive according
                                 * to critRat */
    double gpAddlSmm,           /* (I) Additional smm for a group that is
                                 * triggered */
    double *scurveLogistic,     /* (I) ptr to new array of logistic params 
                                 * (currently 6 values; max, min, width,
                                 * inflec,skew,disc decay) for base prepay
                                 * function logistic */ 
    double *seasLogistic,       /* (I) ptr to new array of logistic params 
                                 * (currently 4 values) for seasoning ramp 
                                 * logistic */ 
    double *groupPs,            /* (I) is an array with the six parameters
                                 * for the group distribution*/    
    double *ltAges,             /* (I) praameters used in adjusting the
                                 * refinance ratio */
    double *mgrpLagWgts,        /* (I) Array of wgts to use w/prev refi
                                 * rates */
    TMbsPrepayDefaults **mbsPrepayDefaults);

void EXPORT
FreeMbsPrepayDefaults
   (TMbsPrepayDefaults **mbsPrepayDefaults);


/***************************************************************************
 * PrintMbsPrepayDefaults()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPrepayDefaults
   (long prepayModel,
    TMbsPrepayDefaults *mbsPrepayDefaults); /* (I) structure contains default
                                             * numbers */


/***************************************************************************
 *  CONSTRUCTORS/DESTRUCTORS FOR MODEL-SPECIFIC STRUCTURES
 ***************************************************************************/
int
make_mbs_pp_modelinfo_const
   (double  const_prepay,       /* (I) constant prepay rate to use */
    long     const_prepay_form,  /* (I) form of constant prepay rate;
                                   MBS_PP_SPD_SMM, etc. */
    /* PTR TO NEW/EXISTING STRUCTURE */
    MBS_PP_MODELINFO_CONST
          **model_data          /* (I/O) structure to alloc/modify */
    );

void
free_mbs_pp_modelinfo_const
   (MBS_PP_MODELINFO_CONST
          **model_data);

int
make_mbs_pp_modelinfo_vector
   (long    vector_prepay_form,   /* (I) form of prepay speeds in vector table;
                                    e.g., MBS_PP_SPD_SMM */
    long    vector_num_ages,      /* (I) # age categories (rows) in prepay
                                    vector table */
    long    vector_num_ratemoves, /* (I) # rate move columns in prepay vector
                                    table */
    long    *vector_ages,         /* (I) ages (in months) of rows of vector
                                    table */
    double *vector_ratemoves,    /* (I) rate move sizes (decimal, e.g., +100bp
                                    as 0.01) */
    double **vector_prepays,     /* (I) table of prepay speeds, indexed as:
                                    vector_prepays[i_age][j_ratemove]  */
    long    vector_prepay_lag,    /* (I) # months lag between rate changes
                                    and effect on prepays */
    /* PTR TO NEW/EXISTING STRUCTURE */
    MBS_PP_MODELINFO_VECTOR
         **model_data            /* (I/O) structure to alloc/modify */
    );

void
free_mbs_pp_modelinfo_vector
   (MBS_PP_MODELINFO_VECTOR
         **model_data);

#ifdef __cplusplus
}
#endif

#endif
