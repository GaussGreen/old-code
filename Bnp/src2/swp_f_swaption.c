/* ===================================================================================
   FILENAME:      swp_f_swaption.c

   PURPOSE:       Compute swaption prices  , and implied volatilities
   ===================================================================================
 */
#include "math.h"
#include "num_h_allhdr.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_swap_pricing.h"
#include "swp_h_swap_simple.h"
#include "swp_h_swaption.h"

#define OPTION_EPS 1.0e-15

/* -------------------  SWAPTION PRICE AND DERIVATIVES   -----------------------
 */

/*  ------------------------------------------------------------------------------
 */
/* Compute the price or greeks of a swaption (for use outside SORT) */
Err swp_f_Swaption(long start, long end_nfp, String compStr, String basisStr,
                   double vol, double strike, String recPayStr,
                   String refRateCode, String ycName, String greekStr,
                   String normLogStr, double *result) {
  Err err = NULL;
  SwapDP swapdp;
  SrtReceiverType recPay;
  SrtGreekType greek;
  SrtDiffusionType logNorm;
  SrtCurvePtr yccrv;

  /* Create the SwapDP from start  , end  , comp and basis */
  err = swp_f_initSwapDP(start, end_nfp, compStr, basisStr, &swapdp);
  if (err)
    return err;

  /* Modify the recpayStr into a type */
  err = interp_rec_pay(recPayStr, &recPay);
  if (err)
    return err;

  /* Modify the GreekString into a type */
  greek = PREMIUM;
  if (greekStr) {
    err = interp_greeks(greekStr, &greek);
    if (err)
      return err;
  }

  /* Get the Black Scholes type: normal or lognormal */
  logNorm = SRT_LOGNORMAL;
  if (normLogStr) {
    err = interp_diffusion_type(normLogStr, &logNorm);
    if (err)
      return err;
  }

  /* Get the yield curve */
  yccrv = lookup_curve(ycName);
  if (!yccrv)
    return serror("Could not find yc %s in swp_f_Swaption", ycName);

  /* Main call to the pricing function */
  err = swp_f_Swaption_SwapDP(&swapdp, vol, strike, recPay, refRateCode, ycName,
                              greek, logNorm, result);
  if (err)
    return err;

  /* Return a success message */
  return NULL;

} /* END swp_f_Swaption(...) */

/* -----------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------
 */
/* Compute the price or greeks of a swaption (for use inside SORT)
   Caution: swapdp->spot_lag should be set */

Err swp_f_Swaption_SwapDP(SwapDP *swapdp, double vol, double strike,
                          SrtReceiverType recPay, String refRateCode,
                          String ycName, SrtGreekType greek,
                          SrtDiffusionType logNorm, double *result) {
  Err err = NULL;
  double forward;
  double level;
  double mat;
  double bs;
  Date clcn_date;
  Date fix_date;
  SrtCallPutType call_put;
  SrtCurvePtr yccrv;

  /* Gets the Yield Curve */
  yccrv = lookup_curve(ycName);

  /* Gets calculation date from Yield Curve */
  clcn_date = get_clcndate_from_yldcrv(yccrv);

  /* Get the spot lag from the Yield Curve and sets it into the SwapDP */
  swapdp->spot_lag = get_spotlag_from_curve(yccrv);

  /* Compute the forward swap rate for the relevant swap period */
  err = swp_f_ForwardRate_SwapDP(swapdp, ycName, refRateCode, &forward);
  if (err)
    return err;

  /* Compute the level payment of the swap */
  err = swp_f_Level_SwapDP(swapdp, ycName, &level);
  if (err)
    return err;

  /* Get the maturity of the option: from today to fixing date */
  fix_date = add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  mat = (fix_date - clcn_date) * YEARS_IN_DAY;

  /* Sets the call_put type from receriver or payer (call on rate = payer) */
  call_put = (recPay == SRT_PAYER ? SRT_CALL : SRT_PUT);

  /* Compute the BlackScholes premium or greek */
  if (logNorm == SRT_LOGNORMAL) {
    bs = srt_f_optblksch(forward, strike, vol, mat, 1.0, call_put, greek);
  } else if (logNorm == SRT_NORMAL) {
    bs = srt_f_optblknrm(forward, strike, vol, mat, 1.0, call_put, greek);
  }

  /* Multiply the result by the level on which it is paid */
  *result = bs * level;

  /* Return a success message */
  return NULL;

} /* END swp_f_Swaption_SwapDP(...) */

/* -----------------------------------------------------------------------------------
 */

/* ----------------------  IMPLIED VOLATILITY FOR SWAPTIONS
 * ------------------------ */

/* -----------------------------------------------------------------------------------
 */
/* Compute the implied volatility given the price of a swaption (for use outside
 * SROT) */

Err swp_f_SwaptionImpliedVol(double premium, long start, long end_nfp,
                             String compStr, String basisStr, double strike,
                             String recPayStr, String refRateCode,
                             String ycName, String normLogStr, double *impvol) {
  Err err = NULL;
  SwapDP swapdp;
  SrtReceiverType recPay;
  SrtDiffusionType logNorm;

  /* Create the SwapDP from start  , end  , comp and basis */
  err = swp_f_initSwapDP(start, end_nfp, compStr, basisStr, &swapdp);
  if (err)
    return err;

  /* Modify the recpayStr into a type */
  err = interp_rec_pay(recPayStr, &recPay);
  if (err)
    return err;

  /* Get the Black Scholes type: normal or lognormal */
  logNorm = SRT_LOGNORMAL;
  if (normLogStr) {
    err = interp_diffusion_type(normLogStr, &logNorm);
    if (err)
      return err;
  }

  /* Main call to the pricing function */
  err = swp_f_SwaptionImpliedVol_SwapDP(premium, &swapdp, strike, recPay,
                                        refRateCode, ycName, logNorm, impvol);
  if (err)
    return err;

  /* Return a success message */
  return NULL;
}

/* -----------------------------------------------------------------------------------
 */
/* -----------------------------------------------------------------------------------
 */
/* Compute the implied volatility given the price of a swaption (for use inside
 * SROT) */

Err swp_f_SwaptionImpliedVol_SwapDP(double premium, SwapDP *swapdp,
                                    double strike, SrtReceiverType recPay,
                                    String refRateCode, String ycName,
                                    SrtDiffusionType logNorm, double *impvol) {
  Err err;
  Date clcn_date;
  Date fix_date;
  double forward;
  double level;
  double mat;
  double bs;
  double intrinsic;
  SrtCallPutType call_put;
  SrtCurvePtr yccrv;

  /* Gets the Yield Curve */
  yccrv = lookup_curve(ycName);

  /* Get clcn_date from the yield curve */
  clcn_date = get_clcndate_from_yldcrv(yccrv);

  /* Get the spot lag from the Yield Curve and sets it into the SwapDP */
  swapdp->spot_lag = get_spotlag_from_curve(yccrv);

  /* Compute the intrinsic value of the option  , with a null volatility */
  err = swp_f_Swaption_SwapDP(swapdp, 0.0, strike, recPay, refRateCode, ycName,
                              PREMIUM, logNorm, &intrinsic);
  if (err)
    return err;

  /* Make sure the premium is above the intrinsic */
  if (premium < intrinsic - OPTION_EPS)
    return serror("Premium below intrinsic in swp_f_SwaptionImpliedVol_SwapDP");

  /* Compute the forward swap rate for the relevant swap period */
  err = swp_f_ForwardRate_SwapDP(swapdp, ycName, refRateCode, &forward);
  if (err)
    return err;

  /* Compute the level payment of the swap */
  err = swp_f_Level_SwapDP(swapdp, ycName, &level);
  if (err)
    return err;

  /* Get the maturity of the option: from today to fixing date */
  fix_date = add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  mat = (fix_date - clcn_date) * YEARS_IN_DAY;

  /* Sets the call_put type from receriver or payer (call on rate = payer) */
  call_put = (recPay == SRT_PAYER ? SRT_CALL : SRT_PUT);

  /* Renormalise the BS premium by the level */
  bs = premium / level;

  /* Compute the BS implied vol of the renormalised premium */
  err =
      srt_f_optimpvol(bs, forward, strike, mat, 1.0, call_put, logNorm, impvol);
  if (err)
    return serror("%s in swp_f_SwaptionImpliedVol_SwapDP", err);

  /* Return a success message */
  return NULL;

} /* END swp_f_Swaption_ImpliedVol_SwapDP(...) */

/* ------------------------------------------------------------------------------
 */

/* -------------------------  CASH SETTLED SWAPTIONS
 * -------------------------- */

/* ------------------------------------------------------------------------------
   Static function used to integrate numerically the cash settled swaption
   payoff value under Qlevel (in terms of level)
   ------------------------------------------------------------------------------
 */
static double CS_FWD_SPREAD;
static double CS_COMP;
static double CS_STRIKE;
static double CS_CALL_PUT_FLAG;
static SrtDiffusionType CS_NORM_LOG;
static double CS_STDEV;
static double CS_DRIFT_FWD_SWAP_RATE;
static long CS_NFP;
static double CS_ALPHA;
static double *CS_COVERAGE;

static double cash_swaption_payoff(double x) {
  double fwd_swap_level;
  double fwd_cash_level;
  double l_swap;
  double l_cash;
  double payoff;
  double fwd_swap_rate;
  long i;

  if (CS_NORM_LOG == SRT_LOGNORMAL) {
    fwd_swap_rate = CS_DRIFT_FWD_SWAP_RATE * exp(CS_STDEV * x);
  } else if (CS_NORM_LOG == SRT_NORMAL) {
    fwd_swap_rate = CS_DRIFT_FWD_SWAP_RATE + (CS_STDEV * x);
  }

  /* The cash level is the sum of cf discounted with the swap rate (with broken
   * period) */
  l_cash = 1.0 / (1.0 + fwd_swap_rate / CS_COMP);
  fwd_cash_level = 0.0;
  for (i = CS_NFP; i >= 0; i--) {
    fwd_cash_level += CS_COVERAGE[i];
    fwd_cash_level *= l_cash;
  }
  /*	fwd_cash_level += CS_COVERAGE[1];
          fwd_cash_level *= pow(l_cash  , CS_ALPHA); */

  /* The swap level is computed through the swap df and the zero rate (spread
   * constant) */
  l_swap = 1.0 / (1.0 + (fwd_swap_rate + CS_FWD_SPREAD) / CS_COMP);
  fwd_swap_level = (1.0 - pow(l_swap, CS_NFP + CS_ALPHA)) / fwd_swap_rate;

  payoff = fwd_swap_rate - CS_STRIKE;
  payoff *= CS_CALL_PUT_FLAG;
  if (payoff < 0.0)
    payoff = 0.0;
  payoff *= fwd_cash_level / fwd_swap_level;

  /* Multiply by the Gaussian density */
  payoff *= exp(-0.5 * x * x);

  return payoff;
}

/* -----------------------------------------------------------------------------
 */

Err swp_f_CashSettledSwaption_SwapDP(SwapDP *swapdp, double vol, double strike,
                                     SrtReceiverType rec_pay,
                                     String refRateCode, String ycname,
                                     SrtGreekType greek,
                                     SrtDiffusionType norm_or_log,
                                     double precision, double *answer) {
  double fwd_cash_level, swap_level, fwd_swap_rate, fwd_df, fwd_zero_rate,
      fwd_spread;
  double df_settlement;
  Err err;
  Date fixing_date, settle_date, last_pay_date, today;
  Date *pStart_date;
  Date *pEnd_date;
  Date *pPay_date;
  double *pCvg;
  int num_pay_date;
  int num_date;
  double mat;
  double l;
  double integration;
  double stdev;
  double d;
  double upper_bound, lower_bound;
  double alpha;
  long nfp;
  long i;
  SrtCurvePtr crv;

  /* gets the YC Curve assoicated to the YCName */
  crv = lookup_curve(ycname);

  /* Get today's date from the yield curve */
  today = get_clcndate_from_yldcrv(crv);

  /* Get the swap level to compute expectation under Qlevel*/
  if (err = swp_f_Level_SwapDP(swapdp, ycname, &swap_level))
    return err;

  /* Get the forward swap rate */
  if (err =
          swp_f_ForwardRate_SwapDP(swapdp, ycname, refRateCode, &fwd_swap_rate))
    return err;

  /* Compute the Coverage and then the Coupons */
  swp_f_make_FixedLegDatesAndCoverages(swapdp, today, &pPay_date, &num_pay_date,
                                       &pStart_date, &pEnd_date, &pCvg,
                                       &num_date);

  /* Makes the swap schedule  , and compute the first period fraction (alpha) */
  /*	swap = make_SimpleSwap(swapdp  , 1.0  , 0.0  , 0.0  , today  , SWAP);
          if (swap.dl.type == BROKEN)
          {
                  alpha = coverage(swap.dl.date[0]  , swap.dl.date[1]  ,
     swapdp->basis_code); alpha /= coverage(swap.dl.prev  , swap.dl.date[1]  ,
     swapdp->basis_code); alpha *= (double) swapdp->compd;
          }
          else
          {
                  alpha = 1.0;
          }
  */
  alpha = 1.0;

  /* Makes the swap schedule  , and compute the first period fraction (alpha) */

  /* The number of full periods: forget the first coupon : could be broken */
  nfp = num_date - 1;

  /* Compute the forward cash settled level  , using the forward swap rate to
   * discount */
  l = 1.0 / (1.0 + fwd_swap_rate / (double)swapdp->compd);
  fwd_cash_level = 0.0;
  for (i = nfp; i >= 0; i--) {
    fwd_cash_level += pCvg[i];
    fwd_cash_level *= l;
  }
  /* Add the last (broken) coupon of the swap */
  /*	fwd_cash_level += swap.cpn.d[1]; */
  /* Discount for the remaining (broken) period */
  /*	fwd_cash_level *= pow(l  , alpha); */

  /* The Zero Coupon rate: df = 1/(1+zr/freq) ^ ( alpha + n ) */
  last_pay_date = bus_date_method(swapdp->end, MODIFIED_SUCCEEDING);
  settle_date = bus_date_method(swapdp->start, MODIFIED_SUCCEEDING);
  fwd_df = swp_f_df(settle_date, last_pay_date, ycname);
  fwd_zero_rate =
      (pow(fwd_df, -1.0 / (alpha + nfp)) - 1.0) * (double)swapdp->compd;

  /* The spread between forward zero rate and forward swap rate */
  fwd_spread = fwd_zero_rate - fwd_swap_rate;

  /* Sets a few option parameters */
  fixing_date =
      add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  mat = (fixing_date - today) * YEARS_IN_DAY;
  stdev = vol * sqrt(mat);
  if (stdev != 0.0) {
    if (norm_or_log == SRT_LOGNORMAL) {
      d = log(strike / fwd_swap_rate) / stdev + 0.5 * stdev;
    } else {
      d = (strike - fwd_swap_rate) / stdev;
    }

    /* Sets the static for numerical integration */
    CS_FWD_SPREAD = fwd_spread;
    /*		CS_NFP = (long) swapdp->nfp; */
    CS_COMP = (double)swapdp->compd;
    CS_STRIKE = strike;
    CS_CALL_PUT_FLAG = (rec_pay == SRT_PAYER) ? 1.0 : -1.0;
    CS_NORM_LOG = norm_or_log;
    CS_STDEV = stdev;
    CS_NFP = nfp;
    CS_ALPHA = alpha;
    CS_COVERAGE = pCvg;

    if (norm_or_log == SRT_LOGNORMAL) {
      CS_DRIFT_FWD_SWAP_RATE = fwd_swap_rate * exp(-0.5 * stdev * stdev);
    } else {
      CS_DRIFT_FWD_SWAP_RATE = fwd_swap_rate;
    }

    /* Integrate numerically E[(S-K)+ * fwd_cash_level/fwd_swap_level] */
    if (rec_pay == SRT_PAYER) {
      /* Call on the Rate */
      upper_bound = 5.0;
      lower_bound = d;
      if (lower_bound > upper_bound)
        lower_bound = upper_bound;
    } else if (rec_pay == SRT_RECEIVER) {
      /* Put on the Rate */
      upper_bound = d;
      lower_bound = -5.0;
      if (upper_bound < lower_bound)
        upper_bound = lower_bound;
    }
    integration =
        sm_qsimp(cash_swaption_payoff, lower_bound, upper_bound, precision);
    integration *= INV_SQRT_TWO_PI;

    /* The answer is the discounted expected value (using Level as numeraire) */
    *answer = integration * swap_level;

  } /* END if stdev != 0.0 */
  else {
    *answer = fwd_swap_rate - strike;
    if (rec_pay == SRT_RECEIVER)
      *answer *= -1.0;
    df_settlement = swp_f_df(0, settle_date, ycname);
    *answer *= df_settlement * fwd_cash_level;
    if (*answer < 0.0)
      *answer = 0.0;
  }
  /* Free memory */
  srt_free(pCvg);
  srt_free(pStart_date);
  srt_free(pEnd_date);
  srt_free(pPay_date);
  /*	free_inSimpleSwap(&swap); */

  /* Return a success message */
  return NULL;
}

/* ---------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------
   This function is intended as an interface function for the
   cash_settled_swaptiondp function (compute cash settled swaption prices)
   ---------------------------------------------------------------------------
 */

Err swp_f_CashSettledSwaption(long start, long nfp, String strComp,
                              String strBasis, double strike, double vol,
                              String strRecPay, String refRateCode,
                              String ycname, String strGreek, String strNormLog,
                              double precision, double *answer) {
  Err err = NULL;
  SwapDP swapdp;
  SrtCurvePtr yc_crv;
  SrtReceiverType rec_pay;
  SrtDiffusionType log_or_norm;
  SrtGreekType greek;

  err = swp_f_initSwapDP(start, nfp, strComp, strBasis, &swapdp);
  if (err)
    return err;

  /* Gets the swap yield curve to extract some date information: spot lag... */
  yc_crv = lookup_curve(ycname);
  if (!yc_crv)
    return serror("Could not find %s curve", ycname);
  swapdp.spot_lag = get_spotlag_from_curve(yc_crv);

  err = interp_rec_pay(strRecPay, &rec_pay);
  if (err)
    return err;

  if (strNormLog[0]) {
    err = interp_diffusion_type(strNormLog, &log_or_norm);
    if (err)
      return err;
  } else
    log_or_norm = SRT_LOGNORMAL;

  if (strGreek[0]) {
    err = interp_greeks(strGreek, &greek);
    if (err)
      return err;
  } else
    greek = PREMIUM;

  if (greek != PREMIUM)
    return serror("Cash Settled only available in premium so far");

  if (!precision)
    precision = 1.0e-5;

  err = swp_f_CashSettledSwaption_SwapDP(&swapdp, vol, strike, rec_pay,
                                         refRateCode, ycname, greek,
                                         log_or_norm, precision, answer);
  if (err)
    return err;

  return NULL;
}

/* -----------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------
 */

/* ------------------------- QUANTO CASH SETTLED SWAPTIONS
 * -------------------------- */

/* ------------------------------------------------------------------------------
   Static function used to integrate numerically the cash settled swaption
   payoff value under Qlevel (in terms of level)
   ------------------------------------------------------------------------------
 */
static double QCS_COMP;
static double QCS_STRIKE;
static double QCS_CALL_PUT_FLAG;
static SrtDiffusionType QCS_NORM_LOG;
static double QCS_STDEV;
static double QCS_DRIFT_FWD_SWAP_RATE;
static long QCS_NFP;
static double QCS_ALPHA;
static double *QCS_COVERAGE;

static double quanto_cash_swaption_payoff(double x) {
  double fwd_cash_level;
  double l_cash;
  double payoff;
  double fwd_swap_rate;
  long i;

  if (QCS_NORM_LOG == SRT_LOGNORMAL) {
    fwd_swap_rate = QCS_DRIFT_FWD_SWAP_RATE * exp(QCS_STDEV * x);
  } else if (QCS_NORM_LOG == SRT_NORMAL) {
    fwd_swap_rate = QCS_DRIFT_FWD_SWAP_RATE + (QCS_STDEV * x);
  }

  /* The cash level is the sum of cf discounted with the swap rate (with broken
   * period) */
  l_cash = 1.0 / (1.0 + fwd_swap_rate / QCS_COMP);
  fwd_cash_level = 0.0;
  for (i = QCS_NFP; i >= 0; i--) {
    fwd_cash_level += QCS_COVERAGE[i];
    fwd_cash_level *= l_cash;
  }
  /*	fwd_cash_level += CS_COVERAGE[1];
          fwd_cash_level *= pow(l_cash  , CS_ALPHA); */

  payoff = fwd_swap_rate - QCS_STRIKE;
  payoff *= QCS_CALL_PUT_FLAG;
  if (payoff < 0.0)
    payoff = 0.0;
  payoff *= fwd_cash_level;

  /* Multiply by the Gaussian density */
  payoff *= exp(-0.5 * x * x);

  return payoff;
}

/* -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 */

Err swp_f_QuantoCashSettledSwaption_SwapDP(
    SwapDP *swapdp, double vol, double strike, double correlationfxcms,
    double volfx, SrtReceiverType rec_pay, String refRateCode, String ycname,
    SrtGreekType greek, SrtDiffusionType norm_or_log, double precision,
    double *answer) {
  double fwd_cash_level, swap_level, fwd_swap_rate;
  double df_settlement;
  Err err;
  Date fixing_date, settle_date, last_pay_date, today;
  Date *pStart_date;
  Date *pEnd_date;
  Date *pPay_date;
  double *pCvg;
  int num_pay_date;
  int num_date;
  double mat;
  double l;
  double integration;
  double stdev;
  double d;
  double upper_bound, lower_bound;
  double alpha;
  long nfp;
  long i;
  SrtCurvePtr crv;

  /* gets the YC Curve assoicated to the YCName */
  crv = lookup_curve(ycname);

  /* Get today's date from the yield curve */
  today = get_clcndate_from_yldcrv(crv);

  /* Get the swap level to compute expectation under Qlevel*/
  if (err = swp_f_Level_SwapDP(swapdp, ycname, &swap_level))
    return err;

  /* Compute the Coverage and then the Coupons */
  swp_f_make_FixedLegDatesAndCoverages(swapdp, today, &pPay_date, &num_pay_date,
                                       &pStart_date, &pEnd_date, &pCvg,
                                       &num_date);

  /* The Zero Coupon rate: df = 1/(1+zr/freq) ^ ( alpha + n ) */
  last_pay_date = bus_date_method(swapdp->end, MODIFIED_SUCCEEDING);
  settle_date = bus_date_method(swapdp->start, MODIFIED_SUCCEEDING);

  /* Makes the swap schedule  , and compute the first period fraction (alpha) */
  /*	swap = make_SimpleSwap(swapdp  , 1.0  , 0.0  , 0.0  , today  , SWAP);
          if (swap.dl.type == BROKEN)
          {
                  alpha = coverage(swap.dl.date[0]  , swap.dl.date[1]  ,
     swapdp->basis_code); alpha /= coverage(swap.dl.prev  , swap.dl.date[1]  ,
     swapdp->basis_code); alpha *= (double) swapdp->compd;
          }
          else
          {
                  alpha = 1.0;
          }
  */
  alpha = 1.0;

  /* The number of full periods: forget the first coupon : could be broken */
  swapdp->nfp = nfp = num_date - 1;

  /* Compute the CMSrate */
  if (err = swp_f_cmsratefwd(swapdp, refRateCode, vol, settle_date,
                             DEFAULT_CMS_DELTA, MAX_CMS_SWAPS, ycname,
                             norm_or_log, &fwd_swap_rate))
    return err;

  /* Compute the forward cash settled level  , using the forward swap rate to
   * discount */
  l = 1.0 / (1.0 + fwd_swap_rate / (double)swapdp->compd);
  fwd_cash_level = 0.0;
  for (i = nfp; i >= 0; i--) {
    fwd_cash_level += pCvg[i];
    fwd_cash_level *= l;
  }
  /* Add the last (broken) coupon of the swap */
  /*	fwd_cash_level += swap.cpn.d[1]; */
  /* Discount for the remaining (broken) period */
  /*	fwd_cash_level *= pow(l  , alpha); */

  /* Sets a few option parameters */
  fixing_date =
      add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  mat = (fixing_date - today) * YEARS_IN_DAY;
  stdev = vol * sqrt(mat);
  if (stdev != 0.0) {
    if (norm_or_log == SRT_LOGNORMAL) {
      d = log(strike / fwd_swap_rate) / stdev + 0.5 * stdev;
    } else {
      d = (strike - fwd_swap_rate) / stdev;
    }

    /* Sets the static for numerical integration */
    /*		CS_NFP = (long) swapdp->nfp; */
    QCS_COMP = (double)swapdp->compd;
    QCS_STRIKE = strike;
    QCS_CALL_PUT_FLAG = (rec_pay == SRT_PAYER) ? 1.0 : -1.0;
    QCS_NORM_LOG = norm_or_log;
    QCS_STDEV = stdev;
    QCS_NFP = nfp;
    QCS_ALPHA = alpha;
    QCS_COVERAGE = pCvg;

    if (norm_or_log == SRT_LOGNORMAL) {
      QCS_DRIFT_FWD_SWAP_RATE = fwd_swap_rate * exp(-0.5 * stdev * stdev) *
                                exp(correlationfxcms * volfx * vol * mat);
    } else {
      QCS_DRIFT_FWD_SWAP_RATE =
          fwd_swap_rate + correlationfxcms * volfx * vol * mat;
    }

    /* Integrate numerically E[(S-K)+ * fwd_cash_level/fwd_swap_level] */
    if (rec_pay == SRT_PAYER) {
      /* Call on the Rate */
      upper_bound = 5.0;
      lower_bound = d;
      if (lower_bound > upper_bound)
        lower_bound = upper_bound;
    } else if (rec_pay == SRT_RECEIVER) {
      /* Put on the Rate */
      upper_bound = d;
      lower_bound = -5.0;
      if (upper_bound < lower_bound)
        upper_bound = lower_bound;
    }
    integration = sm_qsimp(quanto_cash_swaption_payoff, lower_bound,
                           upper_bound, precision);
    integration *= INV_SQRT_TWO_PI;

    /* The answer is the discounted expected value (using df as numeraire) */
    df_settlement = swp_f_df(0, settle_date, ycname);
    *answer = integration * df_settlement;

  } /* END if stdev != 0.0 */
  else {
    *answer = fwd_swap_rate - strike;
    if (rec_pay == SRT_RECEIVER)
      *answer *= -1.0;
    df_settlement = swp_f_df(0, settle_date, ycname);
    *answer *= df_settlement * fwd_cash_level;
    if (*answer < 0.0)
      *answer = 0.0;
  }
  /* Free memory */
  srt_free(pCvg);
  srt_free(pStart_date);
  srt_free(pEnd_date);
  srt_free(pPay_date);
  /*	free_inSimpleSwap(&swap); */

  /* Return a success message */
  return NULL;
}

/* ---------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------
   This function is intended as an interface function for the
   cash_settled_swaptiondp function (compute cash settled swaption prices)
   ---------------------------------------------------------------------------
 */

Err swp_f_QuantoCashSettledSwaption(long start, long nfp, String strComp,
                                    String strBasis, double strike, double vol,
                                    double correlationfxcms, double volfx,
                                    String strRecPay, String refRateCode,
                                    String ycname, String strGreek,
                                    String strNormLog, double precision,
                                    double *answer) {
  Err err = NULL;
  SwapDP swapdp;
  SrtCurvePtr yc_crv;
  SrtReceiverType rec_pay;
  SrtDiffusionType log_or_norm;
  SrtGreekType greek;

  err = swp_f_initSwapDP(start, nfp, strComp, strBasis, &swapdp);
  if (err)
    return err;

  /* Gets the swap yield curve to extract some date information: spot lag... */
  yc_crv = lookup_curve(ycname);
  if (!yc_crv)
    return serror("Could not find %s curve", ycname);
  swapdp.spot_lag = get_spotlag_from_curve(yc_crv);

  err = interp_rec_pay(strRecPay, &rec_pay);
  if (err)
    return err;

  if (strNormLog[0]) {
    err = interp_diffusion_type(strNormLog, &log_or_norm);
    if (err)
      return err;
  } else
    log_or_norm = SRT_LOGNORMAL;

  if (strGreek[0]) {
    err = interp_greeks(strGreek, &greek);
    if (err)
      return err;
  } else
    greek = PREMIUM;

  if (greek != PREMIUM)
    return serror("Cash Settled only available in premium so far");

  if (!precision)
    precision = 1.0e-5;

  err = swp_f_QuantoCashSettledSwaption_SwapDP(
      &swapdp, vol, strike, correlationfxcms, volfx, rec_pay, refRateCode,
      ycname, greek, log_or_norm, precision, answer);
  if (err)
    return err;

  return NULL;
}

Err swp_f_QuantoSwaptionPrice(
    double dFwdSwapRate, double dStrike, SrtCallPutType SrtPayRec, Date dToday,
    Date dDomExerciceDate, Date dForFixingDate, long dStartDate,
    Date *pdEndDates, double *pdCvg, int iNumDates, double dCmsNumOfPeriods,
    SrtCompounding SrtForFrequency, SrtDiffusionType SrtDomVolType,
    SrtDiffusionType SrtForVolType, double *pdForCmsVol, double *pdVolForwardFx,
    double *pdCorrelationCmsForwardFx, long lNumOfVolsAndCorrels,
    char *sDomYCName, SRT_Boolean bAdjustCmsVol, double *dQuantoSwaption) {
  Err err = NULL;
  SrtCurvePtr SrtDomYCurve;
  double dOptionMaturity, dCmsMaturity, dAdjRate = 0.0, dCoupon = 0.0,
                                        dCmsVol = 0.0, dCmsOptPremium = 0.0;
  int i, iIndex;

  SrtDomYCurve = lookup_curve(sDomYCName);

  dOptionMaturity = (dDomExerciceDate - dToday) * YEARS_IN_DAY;
  dCmsMaturity = (dForFixingDate - dToday) * YEARS_IN_DAY;

  /* loop on the payment dates */
  for (i = 0; i < iNumDates; i++) {
    /* if there is only one vol value */
    iIndex = min(i, lNumOfVolsAndCorrels);

    /* First Cms rate: using the given vol */
    if (err = swp_f_cmsrate(dFwdSwapRate, dCmsNumOfPeriods, SrtForFrequency,
                            pdForCmsVol[iIndex], dCmsMaturity,
                            (pdEndDates[i] - dStartDate) * YEARS_IN_DAY,
                            DEFAULT_CMS_DELTA, MAX_CMS_SWAPS, SrtForVolType,
                            &dAdjRate))
      return err;

    /* Do we adjust the Vol By using the Cms option */
    if (bAdjustCmsVol == SRT_TRUE) {
      if (SrtPayRec == SRT_PAYER) { /* Foreign Receiver CMS Option  , */
        if (err = swp_f_cmsoption(
                dFwdSwapRate, dCmsNumOfPeriods, SrtForFrequency, dStrike,
                pdForCmsVol[iIndex], dCmsMaturity, SRT_PAYER,
                (pdEndDates[i] - dStartDate) * YEARS_IN_DAY, DEFAULT_CMS_DELTA,
                MAX_CMS_SWAPS, SrtForVolType, &dCmsOptPremium))
          return err;
        /* check the instrinsic value of the result */
        if (dCmsOptPremium < dAdjRate - dStrike)
          dCmsVol = pdForCmsVol[iIndex];
        else
            /* Compute the correct vol for the Cms quanto adjustment */
            if (err = srt_f_optimpvol(dCmsOptPremium, dAdjRate, dStrike,
                                      dCmsMaturity, 1, SrtPayRec, SrtDomVolType,
                                      &dCmsVol))
          return err;
      } else { /* Foreign Payer CMS Option  , strike at Forward */
        if (err = swp_f_cmsoption(
                dFwdSwapRate, dCmsNumOfPeriods, SrtForFrequency, dStrike,
                pdForCmsVol[iIndex], dCmsMaturity, SRT_RECEIVER,
                (pdEndDates[i] - dStartDate) * YEARS_IN_DAY, DEFAULT_CMS_DELTA,
                MAX_CMS_SWAPS, SrtForVolType, &dCmsOptPremium))
          return err;

        if (dCmsOptPremium < dStrike - dAdjRate)
          dCmsVol = pdForCmsVol[iIndex];
        else
            /* Compute the correct vol for the Cms quanto adjustment */
            if (err = srt_f_optimpvol(dCmsOptPremium, dAdjRate, dStrike,
                                      dCmsMaturity, 1, SrtPayRec, SrtDomVolType,
                                      &dCmsVol))
          return err;
      }
    } else
      dCmsVol = pdForCmsVol[iIndex];

    /* adjust the Cms rate using a quanto effect
       and Compute the coupon according to the diffusion type */
    if (SrtDomVolType == SRT_NORMAL) {
      dAdjRate -= pdCorrelationCmsForwardFx[iIndex] * pdVolForwardFx[iIndex] *
                  dCmsVol * dOptionMaturity;
      dCoupon = srt_f_optblknrm(dAdjRate, dStrike, dCmsVol, dOptionMaturity,
                                1.0, SrtPayRec, PREMIUM);
    } else {
      dAdjRate *= exp(-pdCorrelationCmsForwardFx[iIndex] *
                      pdVolForwardFx[iIndex] * dCmsVol * dOptionMaturity);
      dCoupon = srt_f_optblksch(dAdjRate, dStrike, dCmsVol, dOptionMaturity,
                                1.0, SrtPayRec, PREMIUM);
    }

    dCoupon *= swp_f_df(dToday, pdEndDates[i], sDomYCName) * pdCvg[i];
    *dQuantoSwaption += dCoupon;
  }
  return err;
}

Err swp_f_QuantoSwaption(
    double dFwdSwapRate, double dStrike, char *sPayRec, long dStartDate,
    long dEndDate, double dCmsNumOfPeriods, char *sDomCompounding,
    char *sDomBasis, char *sForCompounding, SrtDiffusionType SrtDomVolType,
    SrtDiffusionType SrtForVolType, double *pdForCmsVol, double *pdVolForwardFx,
    double *pdCorrelationCmsForwardFx, long lNumOfVolsAndCorrels,
    char *sDomYCName, char *sForYCName, SRT_Boolean bAdjustCmsVol,
    double *dQuantoSwaption) {
  Err err = NULL;
  SwapDP SrtSwapdp;
  SrtCurvePtr SrtForYCurve, SrtDomYCurve;
  SrtCompounding SrtForFrequency;
  double *pdCvg;
  Date dCalcDate, dDomFixingDate, dForFixingDate, *pdStartDates, *pdEndDates,
      *pdPayDates;
  int iNumPayDates, iNumDates;
  SrtCallPutType SrtPayRec;

  /* Gets the Yield Curve */
  SrtDomYCurve = lookup_curve(sDomYCName);
  SrtForYCurve = lookup_curve(sForYCName);

  if (SrtDomYCurve == NULL || SrtForYCurve == NULL)
    return "Quanto Swaption - No Domestic Yield Curve";

  /* interp all strings */
  if (err = interp_compounding(sForCompounding, &SrtForFrequency))
    return err;
  if (err = interp_rec_pay(sPayRec, &SrtPayRec))
    return err;

  /* Get today's date from the yield curve */
  dCalcDate = get_clcndate_from_yldcrv(SrtDomYCurve);

  /* Create the swap dp */
  if (err = swp_f_initSwapDP(dStartDate, dEndDate, sDomCompounding, sDomBasis,
                             &SrtSwapdp))
    return err;

  /* Get the spot lag from the Yield Curve and sets it into the SwapDP */
  SrtSwapdp.spot_lag = get_spotlag_from_curve(SrtDomYCurve);

  /* Compute swap dates and the Coverage */
  swp_f_make_FixedLegDatesAndCoverages(&SrtSwapdp, dCalcDate, &pdPayDates,
                                       &iNumPayDates, &pdStartDates,
                                       &pdEndDates, &pdCvg, &iNumDates);

  if (lNumOfVolsAndCorrels < iNumDates)
    lNumOfVolsAndCorrels = 0;

  /* Gets calculation date from Yield Curve */
  dCalcDate = get_clcndate_from_yldcrv(SrtDomYCurve);

  /* Get the maturity of the option: from today to fixing date */
  dDomFixingDate =
      add_unit(SrtSwapdp.start, -(SrtSwapdp.spot_lag), SRT_BDAY, SUCCEEDING);

  dForFixingDate =
      add_unit(SrtSwapdp.start, -get_spotlag_from_curve(SrtForYCurve), SRT_BDAY,
               SUCCEEDING);

  err = swp_f_QuantoSwaptionPrice(
      dFwdSwapRate, dStrike, SrtPayRec, dCalcDate, dDomFixingDate,
      dForFixingDate, dStartDate, pdEndDates, pdCvg, iNumDates,
      dCmsNumOfPeriods, SrtForFrequency, SrtDomVolType, SrtForVolType,
      pdForCmsVol, pdVolForwardFx, pdCorrelationCmsForwardFx,
      lNumOfVolsAndCorrels, sDomYCName, bAdjustCmsVol, dQuantoSwaption);

  free(pdPayDates);
  free(pdStartDates);
  free(pdEndDates);
  free(pdCvg);

  return err;
}