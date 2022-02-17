/* ==============================================================================

   FILENAME:      swp_f_swap_compute.cxx

   PURPOSE:       Provide a few functions for swaps related calculations

   ==============================================================================
 */
#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"
#include "swp_h_irrrng.h"

static SwapType determine_swap_type(Arg_Obj *arg);

Err swp_f_interp_swap_message(String messageStr, SwapOutput *output) {
  Err err = NULL;

  strupper(messageStr);
  strip_white_space(messageStr);

  if (!strcmp(messageStr, "PV")) {
    *output = COMPUTE_PV;
  } else if (!strcmp(messageStr, "FWD_PV")) {
    *output = COMPUTE_FWD_PV;
  } else if (!strcmp(messageStr, "FWD_RATE")) {
    *output = COMPUTE_FWD_RATE;
  } else if (!strcmp(messageStr, "LEVEL")) {
    *output = COMPUTE_LEVEL;
  } else if (!strcmp(messageStr, "MARGIN")) {
    *output = COMPUTE_MARGIN;
  } else if (!strcmp(messageStr, "IRR")) {
    *output = COMPUTE_IRR;
  } else if (!strcmp(messageStr, "FWD_IRR")) {
    *output = COMPUTE_FWD_IRR;
  } else if (!strcmp(messageStr, "DURATION")) {
    *output = COMPUTE_DURATION;
  } else if (!strcmp(messageStr, "CONVEXITY")) {
    *output = COMPUTE_CONVEXITY;
  } else if (!strcmp(messageStr, "MODIFIED_DURATION")) {
    *output = COMPUTE_MODIFIED_DURATION;
  } else if (!strcmp(messageStr, "MODIFIED_CONVEXITY")) {
    *output = COMPUTE_MODIFIED_CONVEXITY;
  } else if ((!strcmp(messageStr, "MATCHING_RATIO")) ||
             (!strcmp(messageStr, "CONVEXITY_DURATION"))) {
    *output = COMPUTE_MODIFIED_CONVEXITY;
  } else {
    return serror("Unknown swap message %s", messageStr);
  }

  return err;
}

/* -------------------------------------------------------------------------------
 */

Err swap_unwind_compute(SrtCurvePtr m, Arg_Obj *arg, SRT_Boolean b, int eod) {
  double ans, *d;
  SwapType swap_type;
  Leg_Obj *leg;
  Swap_Obj *swap;
  int compd, index_start, index_end, i, j, nfp_or_end;
  SrtBasisCode basis;
  int coupon_today, coupon_index;
  double strike, init_not, final_not, factor, cover, leg_payment;
  double value_date;
  Date start, today, prev_date, next_date;
  SRT_Boolean fix_leg_flag, fix_today;

  start = Arg_Dateget(arg, ARG_START);
  today = Arg_Dateget(arg, ARG_TODAY);
  value_date = Arg_Dateget(arg, ARG_VALUE_DATE);
  strike = Arg_get(arg, ARG_STRIKE);
  init_not = Arg_get(arg, ARG_INITIAL_NOT);
  final_not = Arg_get(arg, ARG_FINAL_NOT);
  compd = Arg_Iget(arg, ARG_COMPD);
  basis = Arg_Iget(arg, ARG_BASIS_CODE);
  nfp_or_end = Arg_get(arg, ARG_NUM_FULL_PERIOD);

  fix_leg_flag = ((init_not == 0.0 || final_not == 0.0) ? SRT_YES : SRT_NO);
  if (!fix_leg_flag)
    Arg_set(arg, ARG_STRIKE, 0.0);

  fix_today = ((start == today && eod == SRT_YES) ? SRT_YES : SRT_NO);
  if (start < today || (fix_today == SRT_YES)) {
    Arg_set(arg, ARG_TODAY, today - 366);
    swap_type = determine_swap_type(arg);
    swap = generate_swap(arg, swap_type);

    i = (fix_leg_flag ? 0 : 1);
    leg = Swap_field_get(swap, SWAP_LEG, i);
    index_start = Leg_Iget(leg, LEG_INDEX_START);
    index_end = Leg_Iget(leg, LEG_INDEX_END);

    /* 	  loop through to find first coupon date  */

    j = index_start;
    next_date = start;
    do {
      prev_date = next_date;
      next_date = (Date)Leg_field_get(leg, LEG_DATE, j);
      j += 1;
    } while (next_date < today && j < index_end + 1);

    /*	  add end-of-day funtionality  */

    if (today == next_date) {
      coupon_today = 1;
    } else {
      coupon_today = 0;
    }

    if ((coupon_today == 1) && (eod == 1)) {
      prev_date = next_date;
      next_date = (Date)Leg_field_get(leg, LEG_DATE, j);
      if (value_date == today) {
        value_date += 1;
      }
    }

    Arg_set(arg, ARG_START, start);
    Arg_set(arg, ARG_TODAY, today);
    Arg_set(arg, ARG_VALUE_DATE, value_date);

    swap = generate_swap(arg, swap_type);
    ans = value_swap(m, swap) / DEFAULT_NOTIONAL;
    leg = Swap_field_get(swap, SWAP_LEG, i);
    index_start = Leg_Iget(leg, LEG_INDEX_START);

    if (value_date <= next_date) {
      coupon_index = index_start + 1;
      if (coupon_today == 1) {
        coupon_index = index_start + eod;
      }

      factor = Leg_field_get(leg, LEG_DISC_FACTOR, coupon_index);
      if (b == SRT_YES) {
        cover = coverage(prev_date, today, basis);
      } else {
        cover = 0.0;
      }

      if (fix_leg_flag) {
        ans += strike * cover * factor;
      } else {
        leg_payment = Leg_field_get(leg, LEG_PAYMENT, coupon_index);
        ans -= leg_payment * factor / DEFAULT_NOTIONAL;
        ans += -strike * factor * (cover + coverage(today, next_date, basis));
      }
    }

  } else {
    swap_type = determine_swap_type(arg);
    swap = generate_swap(arg, swap_type);
    ans = value_swap(m, swap) / DEFAULT_NOTIONAL;
  }

  free_Swap(swap);

  d = (double *)(Arg_Dget(arg, ARG_RESULT_TARGET));
  *d = ans;

  return (0);
}

Err SWAP_compute(SrtCurvePtr m, Message mess, Arg_Obj *arg) {
  double ans, floating_leg, level_payment, pv, *d;
  Ddate start;
  SwapType swap_type;
  Leg_Obj *leg;
  Swap_Obj *swap;

  switch (mess) {
  case COMPUTE_PV:
    swap_type = determine_swap_type(arg);
    swap = generate_swap(arg, swap_type);
    ans = value_swap(m, swap) / DEFAULT_NOTIONAL;
    break;
  case COMPUTE_FWD_PV:
    start = Arg_get(arg, ARG_START);
    Arg_set(arg, ARG_VALUE_DATE, start);
    swap_type = determine_swap_type(arg);
    swap = generate_swap(arg, swap_type);
    ans = value_swap(m, swap) / DEFAULT_NOTIONAL;
    break;
  case COMPUTE_FWD_RATE:
    Arg_set(arg, ARG_STRIKE, 1.0);
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL);
    Arg_set(arg, ARG_INITIAL_NOT, 1);
    Arg_set(arg, ARG_FINAL_NOT, 1);
    swap = generate_swap(arg, FIXED_NOTIONALS_SWAP);
    value_swap(m, swap);
    leg = Swap_field_get(swap, SWAP_LEG, 0);
    level_payment = Leg_Dget(leg, LEG_VALUE);
    leg = Swap_field_get(swap, SWAP_LEG, 1);
    floating_leg = Leg_Dget(leg, LEG_VALUE);
    ans = -floating_leg / level_payment * 100;
    break;
  case COMPUTE_LEVEL:
    Arg_set(arg, ARG_STRIKE, 1.0);
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL);
    swap = generate_swap(arg, ONE_LEG_FIXED_SWAP);
    ans = value_swap(m, swap) / DEFAULT_NOTIONAL;
    break;
  case COMPUTE_MARGIN:
    pv = Arg_get(arg, ARG_PV);
    Arg_set(arg, ARG_STRIKE, 1.0);
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL);
    swap = generate_swap(arg, ONE_LEG_FIXED_SWAP);
    ans = value_swap(m, swap) / DEFAULT_NOTIONAL;
    ans = pv / ans * 100.0;
    break;
  case COMPUTE_IRR:
    ans = 0;
    break;
  case COMPUTE_FWD_IRR:
    ans = 0;
    break;
  }
  free_Swap(swap);

  d = (double *)(Arg_get(arg, ARG_RESULT_TARGET));
  *d = ans;

  return (0);
}

static SwapType determine_swap_type(Arg_Obj *arg) {
  double initial_not, final_not;

  initial_not = Arg_get(arg, ARG_INITIAL_NOT);
  final_not = Arg_get(arg, ARG_FINAL_NOT);

  if (final_not == 0 && initial_not == 0) {
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL);
    return (ONE_LEG_FIXED_SWAP);
  } else if (final_not == initial_not) {
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL * final_not);
    return (FIXED_FLOATING_SWAP);
  } else if (final_not == 0)
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL * initial_not);
  else
    Arg_set(arg, ARG_NOTIONAL, DEFAULT_NOTIONAL * final_not);

  Arg_set(arg, ARG_FINAL_NOT, DEFAULT_NOTIONAL * final_not);
  Arg_set(arg, ARG_INITIAL_NOT, DEFAULT_NOTIONAL * initial_not);

  return (ONE_LEG_STANDARD_BOND);
}

double get_greek(struct greek_struct greeks, SrtGreekType greek) {
  if (greek == DELTA_FWD) /** 1 fwd delta underlying **/
    return (greeks.fwd_delta_underlying);
  else if (greek == DELTA) /** 2 spot delta underlying **/
    return (greeks.spot_delta_underlying);
  else if (greek == GAMMA_FWD) /** 5 fwd gamma **/
    return (greeks.fwd_gamma);
  else if (greek == GAMMA) /** 6 spot gammma **/
    return (greeks.spot_gamma);
  else if (greek == VEGA) /** 7 vega **/
    return (greeks.vega);
  else if (greek == THETA) /** 8 theta **/
    return (greeks.theta);
  else
    return (greeks.spot_delta_underlying); /* default */
}

double black_scholes_tmp(double fwd_price, double strike, double vol,
                         double mat, double disc, SrtReceiverType rec_or_pay,
                         struct greek_struct *greeks) {
  double d1, d2, prem;
  double nd1, nd2;

  if (vol != 0 && mat != 0) {
    d1 = (log(fwd_price / strike) + vol * vol / 2 * mat) / (vol * sqrt(mat));
    d2 = d1 - vol * sqrt(mat);
    nd1 = norm(d1);
    nd2 = norm(d2);
  } else {
    if ((fwd_price / strike) > 1) {
      d1 = 1000000; /** Equivalent of +oo **/
      d2 = 1000000; /** Equivalent of +oo **/
      nd1 = 1;
      nd2 = 1;
    } else {
      if ((fwd_price / strike) == 1) {
        d1 = 0;
        d2 = 0;
        nd1 = 0.5;
        nd2 = 0.5;
      } else {
        d1 = -1000000; /** Equivalent of +oo */

        d2 = -1000000; /** Equivalent of +oo **/
        nd1 = 0;
        nd2 = 0;
      }
    }
  }

  if (rec_or_pay == SRT_PAYER) /** if call **/
  {
    prem = disc * (fwd_price * nd1 - strike * nd2);
    greeks->spot_delta_underlying = nd1;
    greeks->spot_delta_strike = -nd2;
    if (mat == 0)
      greeks->theta = 0;
    else
      greeks->theta = -1 / 365.25 *
                      (0.5 * fwd_price * disc * vol / sqrt(mat) * gauss(d1) -
                       log(disc) * disc * strike * nd2 / mat);
  }
  if (rec_or_pay == SRT_RECEIVER) /** if put **/
  {
    prem = disc * (fwd_price * (nd1 - 1) - strike * (nd2 - 1));
    greeks->spot_delta_underlying = (nd1 - 1);
    greeks->spot_delta_strike = -(nd2 - 1);
    if (mat == 0)
      greeks->theta = 0;
    else
      greeks->theta = -1 / 365.25 *
                      (0.5 * fwd_price * disc * vol / sqrt(mat) * gauss(d1) -
                       strike * disc * log(disc) / mat * (nd2 - 1));
    nd1 = 1 - nd1;
    nd2 = 1 - nd2;
  }

  greeks->fwd_delta_underlying = greeks->spot_delta_underlying * disc;
  greeks->fwd_delta_strike = greeks->spot_delta_strike * disc;

  /** Gamma and Vega are the same for calls and puts **/

  if (mat == 0) {
    greeks->vega = 0;
    greeks->spot_gamma = 0;
    greeks->fwd_gamma = 0;
    greeks->spot_gamma_cross = 0;
    greeks->fwd_gamma_strike = 0;
    greeks->spot_gamma_strike = 0;
  } else {
    greeks->vega = 0.01 * disc * fwd_price * sqrt(mat) * gauss(d1);
    greeks->spot_gamma = 1 / (fwd_price * disc * vol * sqrt(mat)) * gauss(d1);
    greeks->fwd_gamma = greeks->spot_gamma * disc * disc;
    greeks->spot_gamma_strike = 1 / (strike * vol * sqrt(mat)) * gauss(d2);
    greeks->fwd_gamma_strike = greeks->spot_gamma_strike * disc * disc;
    greeks->spot_gamma_cross = -1 / (strike * vol * sqrt(mat)) * gauss(d1);
  }

  return (prem);
}
