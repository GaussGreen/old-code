/*******************************************************************************
**
**              Copyright (c) 1993 PARIBAS Capital Markets Group
**
********************************************************************************
**
**      SYSTEM:         SRT     SORT      , Fixed Income 2020 Addins
**      SUB_SYSTEM:     SP      Special functions
**
**      MODULE NAME:    SPECIAL_FNC
**
**      PURPOSE:        Import Header
**
**      AUTHOR:         Antoine SAVINE
**
**      DESCRIPTION:    Header for the library of special functions
*******************************************************************************/

/* ========================================================================== */

#ifndef SPFNCTNS_H
#define SPFNCTNS_H

/* ========================================================================== */

/*******************************************************************************
**                      Macros      , Typedefs and Constants
*******************************************************************************/

#include "srt_h_all.h"

#define SRT_IN 1
#define SRT_OUT 0

#define MAX_STRIKE 1000
#define MAX_VOL 10
#define DELTA_STRIKE 5.0e-4
#define TOL 1.0e-5
#define NVOL 5

#define NSOBOL 5000

/******************************************************************************/

/* =============================================================================
  FUNCTION     : srt_f_splocspr(...)
 ============================================================================ */

double srt_f_splocspr(double fwd_cms_quanto_spread, /* forward */
                      double maturity,              /* maturity */
                      double strike,                /* strike */
                      double disc,                  /* discount factor */
                      int num_spots,                /* number of spots */
                      double *spots,                /* array of spots */
                      int num_dates,                /* number of 'maturities' */
                      double *dates,                /* array of 'maturities' */
                      double **vols,                /* array of vols */
                      SrtCallPutType call_put,      /* Call or Put */
                      int num_steps,                /* number of steps */
                      int type,            /* regular = 0      , digital = 1 */
                      int option,          /* euro = 0      , amer = 1 */
                      SrtGreekType greek); /* Greek */

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_sprrf_tok_anal()
 *
 *******************************************************************************/

Err srt_f_sprrf_tok_anal(long today_long, long fixing_long, long pay_long,
                         double CMSfixing, double CMSpay, double cvgpay,
                         double discpay, double margin, double K1, double K2,
                         double slopevol, double fwdvol, long in_out_flag,
                         double *ans);

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_sprrf_tok_mc()
 *
 *******************************************************************************/

Err srt_f_sprrf_tok_mc(long today_long, long fixing_long, long pay_long,
                       double CMSfixing, double CMSpay, double cvgpay,
                       double discpay, double margin, double K1, double K2,
                       double CMSfvol, double CMSpvol, double rhofp,
                       double fwdvol, long in_out_flag, double *ans);

#endif
