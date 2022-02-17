/* ===============================================================================

   FILENAME:      swp_f_bond_cox_compute.cxx

   PURPOSE:       Provide a few function for Bond Option Pricing that duplicates
                  Murex mains functions (CoxClean        , CoxFwd        ,
   CoxFwdYield        , CoxFwdClean
   ================================================================================
 */
#include "math.h"
#include "swp_h_all.h"

double cox_clean_fct(SwapDP p, double coupon, double spot, Date today,
                     Date opt_mat, double vol, double strike, double step,
                     SrtCallPutType call_put, SrtGreekType greek, String swap,
                     String repo, double redemption, double first_coupon)

{
  int i, j, cp, dval;
  double u, d, r_swap, df_swap, r_repo, df_repo, proba, mat1, coup;
  double delta, gamma, result;
  double pstart, pend, ti;
  SwapDP p1;
  DateList list;
  double *acc_int, *s, *price;

  acc_int = (double *)calloc((int)step + 1, sizeof(double));
  s = (double *)calloc((int)step + 1, sizeof(double));
  price = (double *)calloc((int)step + 1, sizeof(double));
  if ((acc_int == NULL) || (s == NULL) || (price == NULL)) {
    srt_free(s);
    srt_free(acc_int);
    srt_free(price);
    return MEMORY_ERR;
  }

  /* days of settlement */
  dval = p.start - today;

  /* This is a spot yield */
  if (spot < 0.35) {
    spot = clean_price_fct(p, coupon, spot, redemption, first_coupon);
  }

  /* look if the option maturity is before the first coupon or not */
  /* to get the proper first coupon at maturity of the option */
  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  if (list.date[1] >= (opt_mat + dval))
    coup = first_coupon;
  else
    coup = coupon;
  srt_free(list.date);

  p1 = p;
  /* This is a strike forward yield */
  if (strike < 0.35) {
    p1.start = opt_mat + dval;
    strike = clean_price_fct(p1, coupon, strike, redemption, coup);
  }

  /* Up and Down Amplitude of the clean price */
  mat1 = (double)(opt_mat + dval - p.start) / 365.0;
  u = exp(vol * sqrt(mat1 / step));
  d = 1 / u;

  coup = first_coupon;
  /* Calculate the accrued interest for each discretization step */
  /* with approximation due at the step which are not an integer */
  for (i = 0; i <= step; i++) {
    ti = (p.start + i * mat1 * 365.0 / step);
    p1.start = (int)ti;
    acc_int[i] = acc_int_fct(p1, coup) + (ti - p1.start) * coup / 365.0;

    /* we go from first coupon to coupon */
    if ((i > 0) && (acc_int[i - 1] > acc_int[i]))
      coup = coupon;
  }

  cp = (call_put == SRT_CALL ? 1 : -1);
  /* Initialise the spot for each node at maturity */
  for (j = 0; j <= step; j++) {
    s[j] = spot * pow(u, step - j) * pow(d, j);

    if (cp * (s[j] - strike) > 0) {
      price[j] = cp * (s[j] - strike);
    } else {
      price[j] = 0.0;
    }
  }

  for (i = (int)step - 1; i >= 0; i--) {
    pstart = (p.start + i * mat1 * 365.0 / step);
    pend = (p.start + (i + 1) * mat1 * 365.0 / step);

    /* Get the different market */
    /* Discount factor for the repo and swap */
    df_repo = swp_f_df((Ddate)(pstart), (Ddate)(pend), repo);
    r_repo = 1 / df_repo;

    df_swap = swp_f_df((Ddate)(pstart), (Ddate)(pend), swap);
    r_swap = 1 / df_swap;

    /* if we discounted from one period to another */
    if (acc_int[i + 1] < acc_int[i]) {
      for (j = 0; j <= i; j++) {
        s[j] = s[j] * d;
        proba = (r_repo -
                 (acc_int[i + 1] - r_repo * acc_int[i] +
                  r_swap * coupon / p.cxxompd) /
                     s[j] -
                 d) /
                (u - d);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;
      }
    }
    /* if we are in the same coupon */
    else {
      for (j = 0; j <= i; j++) {
        s[j] = s[j] * d;
        proba = (r_repo - (acc_int[i + 1] - r_repo * acc_int[i]) / s[j] - d) /
                (u - d);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;
      }
    }

    if (i == 2)
      /* this formula is false:
            gamma =
         (d*price[0]-(u+d)*price[1]+d*price[2])/(spot*spot*(u*u-2*u*d+d*d));

         SL 01/97
      */
      gamma = 2. / (spot * spot * (u - d) * (u - d)) *
              ((1. - d * d) / (u * u - d * d) * price[0]     /* price high */
               - price[1]                                    /* price middle */
               + (u * u - 1.) / (u * u - d * d) * price[2]); /* price low */
    else if (i == 1)
      delta = (price[0] - price[1]) / (spot * (u - d));
  }

  result = price[0];

  srt_free(acc_int);
  srt_free(s);
  srt_free(price);

  if (greek == PREMIUM)
    return (result);
  else if (greek == DELTA)
    return (delta);
  else if (greek == GAMMA)
    return (gamma);
  else
    return (result); /* default */
}

double cox_fwd_yield_fct(SwapDP p, double coupon, double forward, Date today,
                         Date opt_mat, double vol, double strike, double step,
                         SrtCallPutType call_put, SrtGreekType greek,
                         String swap, String repo, double redemption,
                         double first_coupon) {
  int i, j, dval, cp;
  double fwd_yield1, dirty_strike, su, sd, coup;
  double u, d, r_swap, df_swap, r_repo, df_repo, proba, mat1;
  double delta, gamma, result;
  SwapDP fwd_p;
  DateList list;
  double *xcc, *s, *price;
  SrtCurvePtr crv;

  xcc = (double *)calloc((int)step + 1, sizeof(double));
  s = (double *)calloc((int)step + 1, sizeof(double));
  price = (double *)calloc((int)step + 1, sizeof(double));
  if ((xcc == NULL) || (s == NULL) || (price == NULL)) {
    srt_free(s);
    srt_free(xcc);
    srt_free(price);
    return MEMORY_ERR;
  }

  /* days of settlement */
  dval = p.start - today;

  /* look if the option maturity is before the first coupon or not */
  /* to get the proper first coupon at maturity of the option */
  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  if (list.date[1] >= (opt_mat + dval))
    coup = first_coupon;
  else
    coup = coupon;
  srt_free(list.date);

  fwd_p = p;
  fwd_p.start = opt_mat + dval;

  /* we may have a forward yield */
  if (forward < 0.35)
    fwd_yield1 = forward;
  else
    fwd_yield1 = yield_fct(fwd_p, coupon, forward, redemption, coup);

  /* This is a strike forward yield */
  if (strike < 0.35)
    strike = clean_price_fct(fwd_p, coupon, strike, redemption, coup);

  /* Discount factor for the repo and swap */
  crv = lookup_curve(repo);
  df_repo = exp(log(swp_f_df((Ddate)(p.start), (Ddate)(opt_mat + dval), repo)) /
                step);

  crv = lookup_curve(swap);
  df_swap = exp(log(swp_f_df((Ddate)(p.start), (Ddate)(opt_mat + dval), swap)) /
                step);

  r_repo = 1 / df_repo;
  r_swap = 1 / df_swap;

  /* Up and Down Amplitude of the clean price */
  mat1 = (double)(opt_mat + dval - p.start) / 365.0;
  u = exp(vol * sqrt(mat1 / step));
  d = 1 / u;

  /* Dirty strike */
  dirty_strike = strike + acc_int_fct(fwd_p, coup);

  cp = (call_put == SRT_CALL ? 1 : -1);

  /* At maturity */
  for (j = 0; j <= step; j++) {
    xcc[j] = fwd_yield1 * pow(u, step - j) * pow(d, j);
    s[j] = dirty_price_fct(fwd_p, coupon, xcc[j], redemption, coup);

    if (cp * (s[j] - dirty_strike) > 0)
      price[j] = cp * (s[j] - dirty_strike);
    else
      price[j] = 0.0;
  }

  /* In the tree */
  for (i = (int)step - 1; i >= 0; i--) {
    for (j = 0; j <= i; j++) {
      xcc[j] = xcc[j] * d;
      sd = s[j];
      su = s[j + 1];
      s[j] = dirty_price_fct(fwd_p, coupon, xcc[j], redemption, coup);

      proba = (su - s[j]) / (su - sd);
      price[j] = ((1 - proba) * price[j + 1] + proba * price[j]) * df_swap;
    }

    if (i == 2)
      gamma = (price[2] - price[1]) / (s[2] - s[1]) -
              (price[1] - price[0]) / (s[1] - s[0]);
    else if (i == 1) {
      delta = (price[1] - price[0]) / (s[1] - s[0]);
      gamma = gamma / (s[1] - s[0]);
    }
  }

  result = price[0];

  srt_free(xcc);
  srt_free(s);
  srt_free(price);

  if (greek == PREMIUM)
    return (result);
  else if (greek == DELTA)
    return (delta);
  else if (greek == GAMMA)
    return (gamma);
  else
    return (result); /* default */
}

double cox_yield_fct(SwapDP p, double coupon, double spot, Date today,
                     Date opt_mat, double vol, double strike, double step,
                     SrtCallPutType call_put, SrtGreekType greek, String swap,
                     String repo, double redemption, double first_coupon)

{
  int i, j, dval, cp;
  double yield, su, sd;
  double u, d, proba, pstart, pend, ti, u1, d1, mat1, r_repo, df_repo, r_swap,
      df_swap;
  double coup, result;
  double delta, gamma;
  double step_first_coupon;
  SwapDP p1;
  DateList calc_list;
  double *xcc, *s, *price, *acc_int;
  SrtCurvePtr crv_repo, crv_swap;

  acc_int = (double *)calloc((int)step + 1, sizeof(double));
  xcc = (double *)calloc((int)step + 1, sizeof(double));
  s = (double *)calloc((int)step + 1, sizeof(double));
  price = (double *)calloc((int)step + 1, sizeof(double));
  if ((xcc == NULL) || (s == NULL) || (price == NULL) || (acc_int == NULL)) {
    srt_free(s);
    srt_free(xcc);
    srt_free(price);
    return MEMORY_ERR;
  }

  /* days of settlement */
  dval = p.start - today;

  /* This is a spot yield */
  if (spot < 0.35) {
    yield = spot;
  } else {
    yield = yield_fct(p, coupon, spot, redemption, first_coupon);
  }

  /* look if the option maturity is before the first coupon or not */
  /* to get the proper first coupon at maturity of the option */
  calc_list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  if (opt_mat + dval >= calc_list.date[1])
    coup = coupon;
  else
    coup = first_coupon;

  p1 = p;
  p1.start = opt_mat + dval;

  /* This is a strike forward yield */
  if (strike < 0.35) {
    strike = clean_price_fct(p1, coupon, strike, redemption, coup);
  }

  /* Up and Down Amplitude of the yield */
  mat1 = (double)(opt_mat + dval - p.start) / 365.0;
  u = exp(vol * sqrt(mat1 / step));
  d = 1 / u;

  coup = first_coupon;
  /* Calculate the accrued interest for each discretization step */
  /* with approximation due at the step wich are not an integer */
  for (i = 0; i <= step; i++) {
    ti = (p.start + i * mat1 * 365.0 / step);

    /* To know the first time we go away the first coupon */
    if (ti <= calc_list.date[1]) {
      step_first_coupon = step;
    }

    p1.start = (int)ti;
    acc_int[i] = acc_int_fct(p1, coup) + (ti - p1.start) * coup / 365.0;

    /* we go from first coupon to coupon */
    if ((i > 0) && (acc_int[i - 1] > acc_int[i]))
      coup = coupon;
  }

  /* Gain the new first coupon */
  if (opt_mat + dval >= calc_list.date[1])
    coup = coupon;
  else
    coup = first_coupon;

  cp = (call_put == SRT_CALL ? 1 : -1);
  /* Initialise the yield and price for each node at maturity */
  for (j = 0; j <= step; j++) {
    xcc[j] = yield * pow(u, step - j) * pow(d, j);
    s[j] = clean_price_fct(p1, coupon, xcc[j], redemption, coup);

    if (cp * (s[j] - strike) > 0)
      price[j] = cp * (s[j] - strike);
    else
      price[j] = 0.0;
  }

  for (i = (int)step - 1; i >= 0; i--) {
    pstart = (p.start + i * mat1 * 365.0 / step);
    pend = (p.start + (i + 1) * mat1 * 365.0 / step);

    /* Get the different market */
    /* Discount factor for the repo and swap */
    crv_repo = lookup_curve(repo);
    df_repo = swp_f_df((Ddate)(pstart), (Ddate)(pend), repo);
    r_repo = 1 / df_repo;

    crv_swap = lookup_curve(swap);
    df_swap = swp_f_df((Ddate)(pstart), (Ddate)(pend), swap);
    r_swap = 1 / df_swap;

    p1.start = (int)pstart;

    /* we go backward */
    if (step == step_first_coupon) {
      coup = first_coupon;
    }

    /* if we discounted from one period to another */
    if (acc_int[i + 1] < acc_int[i]) {
      for (j = 0; j <= i; j++) {
        xcc[j] = xcc[j] * d;
        sd = s[j]; /* xcc[j] > xcc[j+1] */
        su = s[j + 1];
        s[j] = clean_price_fct(p1, coupon, xcc[j], redemption, coup);

        d1 = sd / s[j];
        u1 = su / s[j];
        proba = (r_repo -
                 (acc_int[i + 1] - r_repo * acc_int[i] +
                  r_swap * coupon / p.cxxompd) /
                     s[j] -
                 d1) /
                (u1 - d1);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;
      }
    }
    /* if we are in the same coupon */
    else {
      for (j = 0; j <= i; j++) {
        xcc[j] = xcc[j] * d;
        sd = s[j];
        su = s[j + 1];
        s[j] = clean_price_fct(p1, coupon, xcc[j], redemption, coup);

        d1 = sd / s[j];
        u1 = su / s[j];

        proba = (r_repo - (acc_int[i + 1] - r_repo * acc_int[i]) / s[j] - d1) /
                (u1 - d1);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;
      }
    }

    if (i == 2)
      gamma = (price[2] - price[1]) / (s[2] - s[1]) -
              (price[1] - price[0]) / (s[1] - s[0]);
    else if (i == 1) {
      delta = (price[1] - price[0]) / (s[1] - s[0]);
      gamma = gamma / (s[1] - s[0]);
    }
  }

  result = price[0];

  srt_free(xcc);
  srt_free(s);
  srt_free(price);
  srt_free(acc_int);

  srt_free(calc_list.date);

  if (greek == PREMIUM)
    return (result);
  else if (greek == DELTA)
    return (delta);
  else if (greek == GAMMA)
    return (gamma);
  else
    return (result); /* default */
}

double am_cox_clean_fct(SwapDP p, double coupon, double spot, Date today,
                        Date opt_mat, double vol, double strike, double step,
                        SrtCallPutType call_put, SrtGreekType greek,
                        String swap, String repo, double redemption,
                        double first_coupon)

{
  int i, j, cp, dval;
  double u, d, r_swap, df_swap, r_repo, df_repo, proba, mat1, coup;
  double delta, gamma, ti, result;
  SwapDP p1;
  DateList list;
  double *acc_int, *s, *price;

  acc_int = (double *)calloc((int)step + 1, sizeof(double));
  s = (double *)calloc((int)step + 1, sizeof(double));
  price = (double *)calloc((int)step + 1, sizeof(double));
  if ((acc_int == NULL) || (s == NULL) || (price == NULL)) {
    srt_free(s);
    srt_free(acc_int);
    srt_free(price);
    return MEMORY_ERR;
  }

  /* days of settlement */
  dval = p.start - today;

  /* This is a spot yield */
  if (spot < 0.35) {
    spot = clean_price_fct(p, coupon, spot, redemption, first_coupon);
  }

  /* look if the option maturity is before the first coupon or not */
  /* to get the proper first coupon at maturity of the option */
  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  if (list.date[1] >= (opt_mat + dval))
    coup = first_coupon;
  else
    coup = coupon;
  srt_free(list.date);

  p1 = p;
  /* This is a strike forward yield */
  if (strike < 0.35) {
    p1.start = opt_mat + dval;
    strike = clean_price_fct(p1, coupon, strike, redemption, coup);
  }

  /* Discount factor for the repo and swap */
  df_repo = exp(log(swp_f_df((Ddate)(p.start), (Ddate)(opt_mat + dval), repo)) /
                step);

  df_swap = exp(log(swp_f_df((Ddate)(p.start), (Ddate)(opt_mat + dval), swap)) /
                step);

  r_repo = 1 / df_repo;
  r_swap = 1 / df_swap;

  /* Up and Down Amplitude of the clean price */
  mat1 = (double)(opt_mat + dval - p.start) / 365.0;
  u = exp(vol * sqrt(mat1 / step));
  d = 1 / u;

  coup = first_coupon;
  /* Calculate the accrued interest for each discretization step */
  /* with approximation due at the step wich are not an integer */
  for (i = 0; i <= step; i++) {
    ti = (p.start + i * mat1 * 365.0 / step);
    p1.start = (int)ti;
    acc_int[i] = acc_int_fct(p1, coup) + (ti - p1.start) * coup / 365.0;

    /* we go from first coupon to coupon */
    if ((i > 0) && (acc_int[i - 1] > acc_int[i]))
      coup = coupon;
  }

  cp = (call_put == SRT_CALL ? 1 : -1);

  /* Initialise the spot for each node at maturity */
  for (j = 0; j <= step; j++) {
    s[j] = spot * pow(u, step - j) * pow(d, j);

    if (cp * (s[j] - strike) > 0) {
      price[j] = cp * (s[j] - strike);
    } else {
      price[j] = 0.0;
    }
  }

  for (i = (int)step - 1; i >= 0; i--) {
    /* if we are in the same coupon */
    if (acc_int[i + 1] < acc_int[i]) {
      for (j = 0; j <= i; j++) {
        s[j] = s[j] * d;
        proba = (r_repo -
                 (acc_int[i + 1] - r_repo * acc_int[i] +
                  r_swap * coupon / p.cxxompd) /
                     s[j] -
                 d) /
                (u - d);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;

        if ((price[j] < (s[j] - strike)) && (call_put == SRT_CALL))
          price[j] = s[j] - strike;
        else if ((price[j] < (strike - s[j])) && (call_put == SRT_PUT))
          price[j] = strike - s[j];
      }
    }
    /* if we discounted from one period to another */
    else {
      for (j = 0; j <= i; j++) {
        s[j] = s[j] * d;
        proba = (r_repo - (acc_int[i + 1] - r_repo * acc_int[i]) / s[j] - d) /
                (u - d);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;

        if ((price[j] < (s[j] - strike)) && (call_put == SRT_CALL))
          price[j] = s[j] - strike;
        else if ((price[j] < (strike - s[j])) && (call_put == SRT_PUT))
          price[j] = strike - s[j];
      }
    }

    if (i == 2)
      gamma = (d * price[0] - (u + d) * price[1] + d * price[2]) /
              (spot * spot * (u * u - 2 * u * d + d * d));
    else if (i == 1)
      delta = (price[0] - price[1]) / (spot * (u - d));
  }

  result = price[0];

  srt_free(acc_int);
  srt_free(s);
  srt_free(price);

  if (greek == PREMIUM)
    return (result);
  else if (greek == DELTA)
    return (delta);
  else if (greek == GAMMA)
    return (gamma);
  else
    return (result); /* default */
}

double am_cox_yield_fct(SwapDP p, double coupon, double spot, Date today,
                        Date opt_mat, double vol, double strike, double step,
                        SrtCallPutType call_put, SrtGreekType greek,
                        String swap, String repo, double redemption,
                        double first_coupon)

{
  int i, j, dval, cp;
  double yield, su, sd;
  double u, d, proba, pstart, pend, ti, u1, d1, mat1, r_repo, df_repo, r_swap,
      df_swap;
  double coup, result;
  double delta, gamma;
  double step_first_coupon;
  SwapDP p1;
  DateList calc_list;
  double *xcc, *s, *price, *acc_int;

  acc_int = (double *)calloc((int)step + 1, sizeof(double));
  xcc = (double *)calloc((int)step + 1, sizeof(double));
  s = (double *)calloc((int)step + 1, sizeof(double));
  price = (double *)calloc((int)step + 1, sizeof(double));
  if ((xcc == NULL) || (s == NULL) || (price == NULL) || (acc_int == NULL)) {
    srt_free(s);
    srt_free(xcc);
    srt_free(price);
    return MEMORY_ERR;
  }

  /* days of settlement */
  dval = p.start - today;

  /* This is a spot yield */
  if (spot < 0.35) {
    yield = spot;
  } else {
    yield = yield_fct(p, coupon, spot, redemption, first_coupon);
  }

  /* look if the option maturity is before the first coupon or not */
  /* to get the proper first coupon at maturity of the option */
  calc_list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  if (opt_mat + dval >= calc_list.date[1])
    coup = coupon;
  else
    coup = first_coupon;

  p1 = p;
  p1.start = opt_mat + dval;

  /* This is a strike forward yield */
  if (strike < 0.35) {
    strike = clean_price_fct(p1, coupon, strike, redemption, coup);
  }

  /* Up and Down Amplitude of the yield */
  mat1 = (double)(opt_mat + dval - p.start) / 365.0;
  u = exp(vol * sqrt(mat1 / step));
  d = 1 / u;

  coup = first_coupon;
  /* Calculate the accrued interest for each discretization step */
  /* with approximation due at the step wich are not an integer */
  for (i = 0; i <= step; i++) {
    ti = (p.start + i * mat1 * 365.0 / step);

    /* To know the first time we go away the first coupon */
    if (ti <= calc_list.date[1]) {
      step_first_coupon = step;
    }

    p1.start = (int)ti;
    acc_int[i] = acc_int_fct(p1, coup) + (ti - p1.start) * coup / 365.0;

    /* we go from first coupon to coupon */
    if ((i > 0) && (acc_int[i - 1] > acc_int[i]))
      coup = coupon;
  }

  /* Gain the new first coupon */
  if (opt_mat + dval >= calc_list.date[1])
    coup = coupon;
  else
    coup = first_coupon;

  cp = (call_put == SRT_CALL ? 1 : -1);
  /* Initialise the yield and price for each node at maturity */
  for (j = 0; j <= step; j++) {
    xcc[j] = yield * pow(u, step - j) * pow(d, j);
    s[j] = clean_price_fct(p1, coupon, xcc[j], redemption, coup);

    if (cp * (s[j] - strike) > 0)
      price[j] = cp * (s[j] - strike);
    else
      price[j] = 0.0;
  }

  for (i = (int)step - 1; i >= 0; i--) {
    pstart = (p.start + i * mat1 * 365.0 / step);
    pend = (p.start + (i + 1) * mat1 * 365.0 / step);

    /* Get the different market */
    /* Discount factor for the repo and swap */
    df_repo = swp_f_df((Ddate)(pstart), (Ddate)(pend), repo);
    r_repo = 1 / df_repo;

    df_swap = swp_f_df((Ddate)(pstart), (Ddate)(pend), swap);
    r_swap = 1 / df_swap;

    p1.start = (int)pstart;

    /* we go backward */
    if (step == step_first_coupon) {
      coup = first_coupon;
    }

    /* if we discounted from one period to another */
    if (acc_int[i + 1] < acc_int[i]) {
      for (j = 0; j <= i; j++) {
        xcc[j] = xcc[j] * d;
        sd = s[j]; /* xcc[j] > xcc[j+1] */
        su = s[j + 1];
        s[j] = clean_price_fct(p1, coupon, xcc[j], redemption, coup);

        d1 = sd / s[j];
        u1 = su / s[j];
        proba = (r_repo -
                 (acc_int[i + 1] - r_repo * acc_int[i] +
                  r_swap * coupon / p.cxxompd) /
                     s[j] -
                 d1) /
                (u1 - d1);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;

        if ((price[j] < (s[j] - strike)) && (call_put == SRT_CALL))
          price[j] = s[j] - strike;
        else if ((price[j] < (strike - s[j])) && (call_put == SRT_PUT))
          price[j] = strike - s[j];
      }
    }
    /* if we are in the same coupon */
    else {
      for (j = 0; j <= i; j++) {
        xcc[j] = xcc[j] * d;
        sd = s[j];
        su = s[j + 1];
        s[j] = clean_price_fct(p1, coupon, xcc[j], redemption, coup);

        d1 = sd / s[j];
        u1 = su / s[j];

        proba = (r_repo - (acc_int[i + 1] - r_repo * acc_int[i]) / s[j] - d1) /
                (u1 - d1);

        price[j] = (proba * price[j] + (1 - proba) * price[j + 1]) * df_swap;

        if ((price[j] < (s[j] - strike)) && (call_put == SRT_CALL))
          price[j] = s[j] - strike;
        else if ((price[j] < (strike - s[j])) && (call_put == SRT_PUT))
          price[j] = strike - s[j];
      }
    }

    if (i == 2)
      gamma = (price[2] - price[1]) / (s[2] - s[1]) -
              (price[1] - price[0]) / (s[1] - s[0]);
    else if (i == 1) {
      delta = (price[1] - price[0]) / (s[1] - s[0]);
      gamma = gamma / (s[1] - s[0]);
    }
  }

  result = price[0];

  srt_free(xcc);
  srt_free(s);
  srt_free(price);
  srt_free(acc_int);

  srt_free(calc_list.date);

  if (greek == PREMIUM)
    return (result);
  else if (greek == DELTA)
    return (delta);
  else if (greek == GAMMA)
    return (gamma);
  else
    return (result); /* default */
}
