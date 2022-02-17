/* ===================================================================================
   FILENAME:      srt_f_closedform.cxx

   PURPOSE:       Compute swaption        , bond option        ,caps/floors
   prices:
                                                - build the cash flow schedule
                                                - call the relevant pricing
   function (LGM or EtaBeta)
   ===================================================================================
 */

#include "srt_h_all.h"
#include "srt_h_betaetaclsdfrm.h"
#include "srt_h_closedform.h"
#include "srt_h_lgmclsdfrm.h"
#include "srt_h_quickclsdfrm.h"
#include "swp_h_swap_generic.h"

static Err clsdfrm_FixedAndSpreadLegs(GenSwapLeg *fixedleg,
                                      GenSwapLeg *spreadleg,
                                      SrtReceiverType rec_pay, StructType type,
                                      SrtUndPtr und, SrtMdlType mdl_type,
                                      SrtMdlDim mdl_dim, double *answer);

/* -------------------------------------------------------------------------
   FUNCTION: srt_f_closedform

   PURPOSE:  LGM model --- closed form prices for swaptions        , bond
   options        , caps and floors. When dealing with a swaption        , a cap
   or a floor        , the bond strike has no impact. Please note that if it is
   an option on a bond        , the bond strike is assumed to be a clean strike.
   The coupons are calculated with the coverage when generating the fixed leg:
   it is a clean bond generated
   ------------------------------------------------------------------------- */

Err srt_f_closed_form(SrtUndPtr und, SrtMdlType mdl_type, SrtMdlDim mdl_dim,
                      SwapDP *sdp, double strike, double bond_strike,
                      SrtReceiverType rec_pay, StructType type,
                      String ref_rate_code, double *answer) {
  double initial;
  double final;
  Date today;
  Err err;
  GenSwapLeg fixedleg;
  GenSwapLeg spreadleg;

  /* Get today from the underlying */
  today = get_today_from_underlying(und);

  /* Sets the initial and final exchanges according to type */
  if (type == SWAPTION) {
    initial = 1.0;
    final = 1.0;
  } else if (type == CAPFLOOR) {
    initial = 0.0;
    final = 0.0;
  } else if (type == BOND_OPTION) {
    initial = bond_strike;
    final = 1.0;
  }

  /* Make the fixed leg (coupon of strike) with initial and final */
  err = swp_f_make_FixedAndNotionalsLeg(sdp, strike, initial, final, today,
                                        &fixedleg);
  if (err)
    return err;

  /* For the moment: assume receive fixed leg */
  fixedleg.rec_pay = SRT_RECEIVER;

  /* Make the spread leg: dates        , times        , cvg        , spreads
   * cash vs Livor...*/
  err = swp_f_make_SpreadLeg(sdp->start, sdp->end, today, ref_rate_code,
                             &spreadleg);
  if (err) {
    swp_f_freein_GenSwapLeg(&fixedleg);
    return err;
  }

  /* For the moment: assume pay floating */
  spreadleg.rec_pay = SRT_PAYER;

  /* Get the price for this configuration */
  err = clsdfrm_FixedAndSpreadLegs(&fixedleg, &spreadleg, rec_pay, type, und,
                                   mdl_type, mdl_dim, answer);

  /* Free memory */
  swp_f_freein_GenSwapLeg(&fixedleg);
  swp_f_freein_GenSwapLeg(&spreadleg);

  /* Return the error message */
  return err;

} /* END srt_f_closed_from(...) */

/* ---------------------------------------------------------------------------
 */

static Err clsdfrm_FixedAndSpreadLegs(GenSwapLeg *fixedleg,
                                      GenSwapLeg *spreadleg,
                                      SrtReceiverType rec_pay, StructType type,
                                      SrtUndPtr und, SrtMdlType mdl_type,
                                      SrtMdlDim mdl_dim, double *answer) {
  double price;
  TermStruct *ts;
  String yc_name = NULL;
  Err err;
  GenSwapLeg bigleg;
  SrtCurvePtr yldcrv;

  /* Get the Term Structure from the underlying */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  /* Merge fixed and floating in one single big leg */
  err = swp_f_merge_SwapLegs(fixedleg, spreadleg, &bigleg);
  if (err)
    return err;

  /* Get the curve and the spot lag */
  yc_name = get_ycname_from_irund(und);
  yldcrv = lookup_curve(yc_name);
  bigleg.sdp.spot_lag = get_spotlag_from_underlying(und);

  /* MAke the times from today to the swap dates */
  err = time_list(bigleg.pay_date, bigleg.today, &(bigleg.time));

  /* MAke the maturities from today to the fixing dates */
  err = maturity_list(bigleg.fixing_date, bigleg.today, &(bigleg.mat));

  /* Compute the discount factors for each cash flow date in the big list */
  err = df_list(bigleg.pay_date, yc_name, &(bigleg.df));

  switch (type) {
  case SWAPTION:
  case BOND_OPTION:
    /* Price (sends the bond clean strike: it is negative in payment.d[0]) */
    if (mdl_type == LGM) {
      price = srt_f_lgm_coupon_bond_option(
          bigleg.leg_length - 1, -bigleg.payment.d[0], ts, bigleg.mat.d[0],
          bigleg.payment.d, bigleg.time.d, bigleg.df.d, rec_pay, mdl_dim);
    } else if (mdl_type == NEWLGM) {
      price = srt_f_newlgm_coupon_bond_option(
          bigleg.leg_length - 1, -bigleg.payment.d[0], ts, bigleg.mat.d[0],
          bigleg.payment.d, bigleg.time.d, bigleg.df.d, rec_pay, mdl_dim);
    } else if (mdl_type == ETABETA) {
      price = srt_f_etabeta_coupon_bond_option(
          bigleg.leg_length - 1, -bigleg.payment.d[0], ts, bigleg.mat.d[0],
          bigleg.payment.d, bigleg.time.d, bigleg.df.d, rec_pay, mdl_dim);
    } else {
      price = srt_f_quick_beta_bond_option(
          bigleg.leg_length - 1, -bigleg.payment.d[0], ts, bigleg.mat.d[0],
          bigleg.payment.d, bigleg.time.d, bigleg.df.d, rec_pay, mdl_type,
          mdl_dim, yc_name);
    }

    break;

  case CAPFLOOR:
    if (mdl_type == LGM) {
      price = srt_f_lgm_capfloor(
          ts, bigleg.mat.d, bigleg.time.d, bigleg.df.d, bigleg.payment.d,
          bigleg.leg_length - 1, /* if n caplets        , n+1 dates */
          rec_pay, mdl_dim);
    } else if (mdl_type == NEWLGM) {
      price = srt_f_lgm_capfloor(
          ts, bigleg.mat.d, bigleg.time.d, bigleg.df.d, bigleg.payment.d,
          bigleg.leg_length - 1, /* if n caplets        , n+1 dates */
          rec_pay, mdl_dim);
    } else if (mdl_type == ETABETA) {
      price = srt_f_etabeta_capfloor(
          ts, bigleg.mat.d, bigleg.time.d, bigleg.df.d, bigleg.payment.d,
          bigleg.leg_length - 1, /* if n caplets        , n+1 dates */
          rec_pay, mdl_dim);
    } else {
      price = srt_f_quick_beta_capfloor(
          ts, bigleg.mat.d, bigleg.time.d, bigleg.df.d, bigleg.payment.d,
          bigleg.leg_length - 1, /* if n caplets        , n+1 dates */
          rec_pay, mdl_type, mdl_dim, yc_name);
    }
    break;
    /*		case RESETCAPFLOOR:
                            fixing_date = srt_calloc (cfs->dl.len + 1        ,
       sizeof(double));
    */			/* No fixing for last pay date */
    /*			for (i = 0; i <df.len; i++)
                            {
                                    fixing_date[i] = add_unit(cfs->dl.date[i] ,
       cfs->sdp.spot_lag        , SRT_BDAY        , SUCCEEDING);
                            }
                            mkt = lookup_und (und_name);
                            price = gen_ir_resetcap(
                                            mkt        ,
                                            fixing_date        ,
                                            cfs->dl.date        ,
                                            df.d        ,
                                            cfs->cpn.d        ,
                                            df.len        ,
                                            rec_pay        ,
                                            strike);
                            srt_free(fixing_date);

                            if (price == ERROR_IN_RESETCAP)
                            {
                                    free_inDlist(df);
                                    return serror ("Closed form: error in
       resetcap pricing");
                            }
                    break;
    */
  default: /* never get here */
    swp_f_freein_GenSwapLeg(&bigleg);
    ;
    return serror("Closed form: unknown structure type %d", type);
    break;
  }

  swp_f_freein_GenSwapLeg(&bigleg);
  ;

  *answer = price;

  return NULL;

} /* END static Err clsdfrm_FixedAndSpreadLegs(...) */

/* ----------------------------------------------------------------------------
 */
