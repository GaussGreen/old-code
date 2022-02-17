/*************************************************************************************
 * AUTHORS: J.LEBUCHOUX & S.WISNIA
 *
 * FUNCTION: srt_f_approx_CHEYBETA()
 *
 * PURPOSE: compute prices of caplets & swaptions in a CHEYBETA 1F model
 *          Study of the Smile & linkage against the calibration routine
 *
 * METHOD: Small noise expansion of the density probability
 *
 * INPUTS: start_date               start date of the instrument as a long
 *         end_date                 end date of the instrument as a long
 *         frequency                payment frequency as a char*
 *         compounding              compounding convention as a char*
 *         strike                   strike of the instrument as a double
 *         SWAPTIONS/CAPLETS        type of the instrument as a char*
 *         BETA                     exogeneous parameters of the model
 *         und                      underlying
 *         underlying               as a char* to get all the vols & the taus
 *************************************************************************************/

/* NEEDED FUNCTIONS:-reconstruction formula and its first & second derivatives
                    -dichotomy function
                                        -numerical integration procedure */

/*---------------------------------Include
 * files------------------------------------*/
#include "OPFNCTNS.H>
#include "SPFNCTNS.H>
#include "UTALLHDR.H>
#include "math.h"

/*----------------------------------------------------------------------------------*/
/*-----------------------------------Macro
 * definitions------------------------------*/
#define IFR_epsilon 1 /* in days */
#define x1 0.01
#define x2 0.09
#define xacc 0.000001
#define tacc 1.2
#define N 100

typedef Err (*pricing_caplet)(Ddate, Ddate, Ddate, char *, char *, String,
                              String, char *, double, double, double *);

typedef Err (*find_t_star)(Ddate, Ddate, Ddate, double, char *, char *, char *,
                           String, String, char *, double *);

/*----------------------------------------------------------------------------------*/

/*-------------------------COMPUTES THE MID_POINT IN THE
 * APPROX---------------------*/

Err srt_f_mid_point(Ddate fixing_date, Ddate end_date, double strike,
                    char *cap_swaption, String compStr, String refRateCode,
                    char *basis_code, char *yc_name, double *ans) {
  Err err;
  double fra, swap_rate;
  SrtBasisCode val;
  SrtCurvePtr yldcrv;
  StructType type;

  err = interp_basis(basis_code, &val);
  if (err) {
    return err;
  }

  err = interp_struct(cap_swaption, &type);
  if (err) {
    return err;
  }

  yldcrv = lookup_curve(yc_name);

  if (type == CAPFLOOR) {
    fra = swp_f_fwdcash(fixing_date, end_date, val, yc_name);
    *ans = 0.5 * (fra + strike);

  } else if (type == SWAPTION) {
    err = swp_f_ForwardRate((long)fixing_date, (long)end_date, compStr,
                            basis_code, yc_name, refRateCode, &swap_rate);
    if (err) {
      return err;
    }

    *ans = 0.5 * (swap_rate + strike);
  } else {
    return serror("(Fatal):srt_f_mid_point: Unknown dealType");
  }

  return NULL;
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------BETA
 * FUNCTION-------------------------------------*/

Err srt_f_beta_function(Ddate current_time, Ddate maturity, char *und_name,
                        double *beta_t_T) {
  Err err;
  SrtUndPtr und;
  TermStruct *ts;
  Date today;
  double time1, time2;

  und = lookup_und(und_name);

  err = get_underlying_ts(und, &ts);
  if (err) {
    return err;
  }

  today = get_spotdate_from_underlying(und);

  time1 = (double)(maturity - (Ddate)today) * YEARS_IN_DAY;
  time2 = (double)(current_time - (Ddate)today) * YEARS_IN_DAY;

  if (time2 == 0.0) {
    *beta_t_T = Psi_func(time1, ts);
  } else {
    *beta_t_T = (Psi_func(time1, ts) - Psi_func(time2, ts)) / F_func(time2, ts);
  }

  return NULL;
}

/*-------------------------------CAPLET
 * FUNCTION-----------------------------------*/
/* in the function  , current time is the maturity of the IFR  , marutity is set
   14 days after */

Err srt_f_caplet_function(Ddate current_time, Ddate fixing_date, Ddate end_date,
                          char *und_name, char *yc_name, String compStr,
                          String refRateCode, char *basis_code, double phi_t,
                          double r_t, double *ans)

{
  Err err;
  Date today, IFR_maturity;
  SrtCrvPtr yldcrv;
  SrtBasisCode val;
  double fra_coverage, df_fixing, df_end, beta_t_fixing, beta_t_end, IFR,
      spot_date;

  yldcrv = lookup_curve(yc_name);

  today = get_clcndate_from_yldcrv(yldcrv);
  spot_date = (double)get_spotdate_from_yldcrv(yldcrv);

  if (current_time <= (Ddate)spot_date) {
    current_time = (Ddate)spot_date;
  }

  df_fixing = swp_f_df((Ddate)today, (Ddate)fixing_date, yc_name);
  if (df_fixing == SRT_DF_ERROR) {
    return serror("Could not compute df(%d  ,%d)", today, fixing_date);
  }

  df_end = swp_f_df((Ddate)today, (Ddate)end_date, yc_name);
  if (df_end == SRT_DF_ERROR) {
    return serror("Could not compute df(%d  ,%d)", today, end_date);
  }

  err = interp_basis(basis_code, &val);
  if (err) {
    return err;
  }

  IFR_maturity = add_unit((Date)current_time, (long)IFR_epsilon, SRT_DAY,
                          MODIFIED_SUCCEEDING);

  IFR = swp_f_fwdcash(current_time, IFR_maturity, val, yc_name);

  fra_coverage = coverage((Date)fixing_date, (Date)end_date, val);

  err =
      srt_f_beta_function(current_time, fixing_date, und_name, &beta_t_fixing);
  if (err) {
    return err;
  }

  err = srt_f_beta_function(current_time, end_date, und_name, &beta_t_end);
  if (err) {
    return err;
  }

  *ans = (df_fixing / (fra_coverage * df_end)) *
             exp(0.5 * phi_t *
                     (beta_t_end * beta_t_end - beta_t_fixing * beta_t_fixing) -
                 r_t * (beta_t_end - beta_t_fixing) +
                 IFR * (beta_t_end - beta_t_fixing)) -
         1.0 / fra_coverage;

  return NULL;
}

/*-----------------------------------------------------------------------------------*/
/*-------------------FIRST & SECOND DERIVATIVES OF THE CAPLET
 * FUNCTION---------------*/

Err srt_f_caplet_derivative_fct(Ddate current_time, Ddate fixing_date,
                                Ddate end_date, char *und_name, char *yc_name,
                                String compStr, String refRateCode,
                                char *basis_code, double phi_t, double r_t,
                                SrtDerType der_order, double *ans) {
  Err err;
  double beta_t_fixing, beta_t_end, rf_fra, fra_coverage;

  SrtBasisCode val;

  err = interp_basis(basis_code, &val);
  if (err) {
    return err;
  }

  err =
      srt_f_beta_function(current_time, fixing_date, und_name, &beta_t_fixing);
  if (err) {
    return err;
  }

  err = srt_f_beta_function(current_time, end_date, und_name, &beta_t_end);
  if (err) {
    return err;
  }

  fra_coverage = coverage((Date)fixing_date, (Date)end_date, val);

  err = srt_f_caplet_function(current_time, fixing_date, end_date, und_name,
                              yc_name, compStr, refRateCode, basis_code, phi_t,
                              r_t, &rf_fra);
  if (err) {
    return err;
  }

  if (der_order == SRT_FIRSTDER) {
    *ans = (beta_t_fixing - beta_t_end) * (rf_fra + (1.0 / fra_coverage));
    *ans *= 0.84;

  } else if (der_order == SRT_SECONDDER) {
    *ans = (beta_t_fixing - beta_t_end) * (beta_t_fixing - beta_t_end) *
           (rf_fra + (1.0 / fra_coverage));
    *ans *= 0.84;
  } else if (der_order == SRT_TERDER) {
    *ans = (beta_t_fixing - beta_t_end) * (beta_t_fixing - beta_t_end) *
           (beta_t_fixing - beta_t_end) * (rf_fra + (1.0 / fra_coverage));
    *ans *= 0.84;
  }

  else {
    return serror(
        "Unknown derivatives order in srt_f_derivatives_caplet_function");
  }
  return NULL;
}

/*-----------------------------------------------------------------------------------*/
Err srt_f_phi_function_of_r(Ddate fixing_date, Ddate end_date, double strike,
                            char *cap_swaption, String compStr,
                            String refRateCode, char *basis_code, char *yc_name,
                            char *und_name, double r, double *ans_phi) {
  Err err;
  SrtUndPtr und;
  SrtCrvPtr yldcrv;
  TermStruct *ts;
  double beta, time, val_G, val_H;
  Date today;

  und = lookup_und(und_name);

  err = get_underlying_ts(und, &ts);
  if (err) {

    return err;
  }

  beta = find_beta((Ddate)end_date, ts);

  yldcrv = lookup_curve(yc_name);
  today = get_clcndate_from_yldcrv(yldcrv);

  time = (double)(fixing_date - (Ddate)today) * YEARS_IN_DAY;

  G_H_func(time, ts, &val_G, &val_H);

  *ans_phi = 0.5 * pow((r), (2.0 * beta)) * F_func(time, ts) *
             F_func(time, ts) * val_G;

  return NULL;
}
/*---------------------------BISSECTION METHOD USED TO SEARCH
 * MID_STATE_VAR----------*/

Err srt_f_bissection(pricing_caplet srt_f_caplet_function, Ddate fixing_date,
                     Ddate end_date, double strike, char *cap_swaption,
                     String compStr, String refRateCode, char *basis_code,
                     char *yc_name, char *und_name, double *rtb)

{

  Err err;
  SrtCrvPtr yldcrv;
  int i;
  double today, dx, f, fmid, xmid, phi_f_r, x_star;

  yldcrv = lookup_curve(yc_name);
  today = get_spotdate_from_yldcrv(yldcrv);

  err = srt_f_mid_point(fixing_date, end_date, strike, cap_swaption, compStr,
                        refRateCode, basis_code, yc_name, &x_star);
  if (err) {
    return err;
  }

  err = srt_f_phi_function_of_r(fixing_date, end_date, strike, cap_swaption,
                                compStr, refRateCode, basis_code, yc_name,
                                und_name, x1, &phi_f_r);
  if (err) {
    return err;
  }

  err = srt_f_caplet_function((Ddate)today, fixing_date, end_date, und_name,
                              yc_name, compStr, refRateCode, basis_code,
                              phi_f_r, x1, &f);
  if (err) {
    return err;
  }

  f -= x_star;

  err = srt_f_phi_function_of_r(fixing_date, end_date, strike, cap_swaption,
                                compStr, refRateCode, basis_code, yc_name,
                                und_name, x2, &phi_f_r);
  if (err) {
    return err;
  }

  err = srt_f_caplet_function((Ddate)today, fixing_date, end_date, und_name,
                              yc_name, compStr, refRateCode, basis_code,
                              phi_f_r, x2, &fmid);
  if (err) {
    return err;
  }

  fmid -= x_star;

  if (f * fmid >= 0.0) {
    return serror("Root must be bracketed in function srt_f_bissection");
  }

  if (f < 0.0) {
    dx = x2 - x1;
    *rtb = x1;
  } else {
    dx = x1 - x2;
    *rtb = x2;
  }

  i = 1;
  while (i <= MAX_ITER && fabs(dx) >= xacc && fmid != 0.0) {
    dx *= 0.5;
    xmid = *rtb + dx;

    err = srt_f_phi_function_of_r(fixing_date, end_date, strike, cap_swaption,
                                  compStr, refRateCode, basis_code, yc_name,
                                  und_name, xmid, &phi_f_r);

    if (err) {
      return err;
    }

    err = srt_f_caplet_function((Ddate)today, fixing_date, end_date, und_name,
                                yc_name, compStr, refRateCode, basis_code,
                                phi_f_r, xmid, &fmid);

    if (err) {
      return err;
    }

    fmid -= x_star;

    if (fmid <= 0.0) {
      *rtb = xmid;
    }
    i++;
  }
  return NULL;
}
/*-----------------------------------------------------------------------------------*/
/*--------------------------------MID state variable
 * function------------------------*/
/* The flag for CHEY in init_ONE_fac_TermStruct has been removed
   the functions G & H are now calculated in CHEY model */

/*-----------------------------------------------------------------------------------*/
/*-------------------------------------COMPUTES THE MID-STATE
 * VARIABLE---------------*/

Err srt_f_mid_state_variables(Ddate fixing_date, Ddate end_date, double strike,
                              char *cap_swaption, char *und_name, char *yc_name,
                              String compStr, String refRateCode,
                              char *basis_code, double *r_star,
                              double *phi_star) {
  Err err;

  err = srt_f_bissection(srt_f_caplet_function, fixing_date, end_date, strike,
                         cap_swaption, compStr, refRateCode, basis_code,
                         yc_name, und_name, r_star);
  if (err) {
    return err;
  }

  err = srt_f_phi_function_of_r(fixing_date, end_date, strike, cap_swaption,
                                compStr, refRateCode, basis_code, yc_name,
                                und_name, *r_star, phi_star);
  if (err) {
    return err;
  }
  return NULL;
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------COMPUTES THE MID VOLATILITY at current
 * time----------------*/

Err srt_f_mid_volatility(Ddate current_time, Ddate fixing_date, Ddate end_date,
                         double strike, char *cap_swaption, char *und_name,
                         char *yc_name, String compStr, String refRateCode,
                         char *basis_code, double r_star, double phi_star,
                         double *sigma) {
  Err err;
  TermStruct *ts;
  SrtUndPtr und;
  double current_volatility, first_derivative, beta;

  und = lookup_und(und_name);
  err = get_underlying_ts(und, &ts);
  if (err) {
    return err;
  }

  err = srt_f_caplet_derivative_fct(current_time, fixing_date, end_date,
                                    und_name, yc_name, compStr, refRateCode,
                                    basis_code, phi_star, r_star, SRT_FIRSTDER,
                                    &first_derivative);
  if (err) {
    return err;
  }

  beta = find_beta((Ddate)end_date, ts);

  current_volatility = find_sig(current_time, ts);

  *sigma = current_volatility * pow(r_star, beta) * first_derivative;
  return NULL;
}

/*-----------------------------------------------------------------------------------*
/*-----------------------------tau_m:
time_rescaling---------------------------------*/
Err srt_f_time_rescaling(Ddate final_time, Ddate fixing_date, Ddate end_date,
                         double strike, char *cap_swaption, char *und_name,
                         char *yc_name, String compStr, String refRateCode,
                         char *basis_code, double *tau_m) {
  Err err;
  double temp, current_time, sigma, r_star, phi_star;
  Date today;
  SrtUndPtr und;

  int i;

  *tau_m = 0.0;
  temp = 0.0;

  err = srt_f_mid_state_variables(fixing_date, end_date, strike, cap_swaption,
                                  und_name, yc_name, compStr, refRateCode,
                                  basis_code, &r_star, &phi_star);
  if (err) {
    return err;
  }

  und = lookup_und(und_name);
  today = get_spotdate_from_underlying(und);

  for (i = 1; i <= N; i++) {
    temp += (double)(final_time - today) / N;

    current_time =
        add_unit((Date)today, (long)temp, SRT_DAY, MODIFIED_SUCCEEDING);

    err =
        srt_f_mid_volatility((Ddate)current_time, fixing_date, end_date, strike,
                             cap_swaption, und_name, yc_name, compStr,
                             refRateCode, basis_code, r_star, phi_star, &sigma);
    if (err) {
      return err;
    }

    *tau_m += ((double)(final_time - today) / N) * sigma * sigma * YEARS_IN_DAY;
  }
  return NULL;
}
/*-----------------------------------------------------------------------------------*/
/*---------------------------------Bissection method for
 * t_star----------------------*/
Err srt_f_t_star_bissection(find_t_star srt_f_time_rescaling, Ddate fixing_date,
                            Ddate end_date, double strike, char *cap_swaption,
                            char *und_name, char *yc_name, String compStr,
                            String refRateCode, char *basis_code,
                            Ddate *t_star) {
  Err err;
  int i;
  double dt, f, fmid, tmid, r_star, phi_star;
  SrtCurvePtr yldcrv;
  Ddate t1, t2, today;

  yldcrv = lookup_curve(yc_name);
  today = get_clcndate_from_yldcrv(yldcrv);

  t1 = (Ddate)today;
  t2 = fixing_date;

  err = srt_f_mid_state_variables(fixing_date, end_date, strike, cap_swaption,
                                  und_name, yc_name, compStr, refRateCode,
                                  basis_code, &r_star, &phi_star);
  if (err) {
    return err;
  }

  f = 0.0;
  /*err = srt_f_time_rescaling(t1  ,fixing_date  ,end_date  ,strike
  ,cap_swaption  , und_name  ,yc_name  ,compStr  ,refRateCode  ,basis_code ,&f);
  if(err)
  {
          return err;
  }*/
  f -= 0.5 * phi_star;

  err = srt_f_time_rescaling(t2, fixing_date, end_date, strike, cap_swaption,
                             und_name, yc_name, compStr, refRateCode,
                             basis_code, &fmid);
  if (err) {
    return err;
  }
  fmid -= 0.5 * phi_star;

  if (f * fmid >= 0.0) {
    return serror(
        "Root must be bracketed in function srt_f__t_star_bissection");
  }

  if (f < 0.0) {
    dt = t2 - t1;
    *t_star = t1;
  } else {
    dt = t1 - t2;
    *t_star = t2;
  }

  i = 1;
  while (i <= MAX_ITER && fabs(dt) >= 1.0 && fmid != 0.0) {
    dt *= 0.5;
    tmid = *t_star + dt;
    err = srt_f_time_rescaling((Ddate)tmid, fixing_date, end_date, strike,
                               cap_swaption, und_name, yc_name, compStr,
                               refRateCode, basis_code, &fmid);

    if (err) {
      return err;
    }

    fmid -= 0.5 * phi_star;

    if (fmid <= 0.0) {
      *t_star = tmid;
    }
    i++;
  }
  return NULL;
}

/*----------------------------------find
 * t_star--------------------------------------*/
Err srt_f_t_star(Ddate fixing_date, Ddate end_date, double strike,
                 char *cap_swaption, char *und_name, char *yc_name,
                 String compStr, String refRateCode, char *basis_code,
                 Ddate *t_star) {
  Err err;

  err = srt_f_t_star_bissection(srt_f_time_rescaling, fixing_date, end_date,
                                strike, cap_swaption, und_name, yc_name,
                                compStr, refRateCode, basis_code, t_star);
  if (err) {
    return err;
  }
  return NULL;
}
/*------------------------------------------------------------------------------------*/

Err srt_f_analytical_price(Ddate fixing_date, Ddate end_date, double strike,
                           char *cap_swaption, char *und_name, char *yc_name,
                           String compStr, String refRateCode, char *basis_code,
                           double *price) {
  Err err;
  SrtBasisCode val;
  SrtCurvePtr yldcrv;
  SrtUndPtr und;
  TermStruct *ts;

  Ddate t_star;

  double fra, today, fra_coverage, df, beta, tau_m, X0, x_star, r_star,
      phi_star, K_star, first_der, second_der, ter_der, ordre_0, ordre_1,
      ordre_2, maturity, sigma_bs;

  yldcrv = lookup_curve(yc_name);
  today = get_clcndate_from_yldcrv(yldcrv);

  err = interp_basis(basis_code, &val);
  if (err) {
    return err;
  }

  fra = swp_f_fwdcash(fixing_date, end_date, val, yc_name);

  err = srt_f_mid_point(fixing_date, end_date, strike, cap_swaption, compStr,
                        refRateCode, basis_code, yc_name, &x_star);
  if (err) {
    return err;
  }

  err = srt_f_time_rescaling(fixing_date, fixing_date, end_date, strike,
                             cap_swaption, und_name, yc_name, compStr,
                             refRateCode, basis_code, &tau_m);
  if (err) {
    return err;
  }
  K_star = (strike - x_star) / sqrt(tau_m);

  X0 = (fra - x_star) / sqrt(tau_m);

  maturity = (fixing_date - today) * YEARS_IN_DAY;

  df = swp_f_df((Ddate)today, (Ddate)end_date, yc_name);
  fra_coverage = coverage((Date)fixing_date, (Date)end_date, val);

  err = srt_f_mid_state_variables(fixing_date, end_date, strike, cap_swaption,
                                  und_name, yc_name, compStr, refRateCode,
                                  basis_code, &r_star, &phi_star);
  if (err) {
    return err;
  }

  err = srt_f_t_star(fixing_date, end_date, strike, cap_swaption, und_name,
                     yc_name, compStr, refRateCode, basis_code, &t_star);
  if (err) {
    return err;
  }

  err = srt_f_caplet_derivative_fct(t_star, fixing_date, end_date, und_name,
                                    yc_name, compStr, refRateCode, basis_code,
                                    phi_star, r_star, SRT_FIRSTDER, &first_der);
  if (err) {
    return err;
  }

  err = srt_f_caplet_derivative_fct(
      t_star, fixing_date, end_date, und_name, yc_name, compStr, refRateCode,
      basis_code, phi_star, r_star, SRT_SECONDDER, &second_der);
  if (err) {
    return err;
  }

  err = srt_f_caplet_derivative_fct(t_star, fixing_date, end_date, und_name,
                                    yc_name, compStr, refRateCode, basis_code,
                                    phi_star, r_star, SRT_TERDER, &ter_der);
  if (err) {
    return err;
  }

  und = lookup_und(und_name);
  err = get_underlying_ts(und, &ts);
  if (err) {
    return err;
  }
  beta = find_beta((Ddate)end_date, ts);

  ordre_0 = (sqrt(tau_m / maturity)) / x_star;

  ordre_1 =
      (second_der / (first_der * first_der)) + (beta / (r_star * first_der));

  ordre_2 = (ter_der / (first_der * first_der * first_der)) -
            ((second_der * second_der) /
             (first_der * first_der * first_der * first_der)) +
            (second_der * beta / (r_star * first_der * first_der * first_der)) +
            (beta * (beta - 1.0) / (r_star * r_star * first_der * first_der));

  sigma_bs = ordre_0 *
             (1.0 +
              (1.0 / (x_star * x_star) - ordre_1 * ordre_1) *
                  (tau_m + 2.0 * (fra - strike) * (fra - strike)) / 24 +
              ordre_2 * (2.0 * tau_m + (fra - strike) * (fra - strike)) / 24);

  *price =
      fra_coverage * df *
      srt_f_optblksch(fra, strike, sigma_bs, maturity, 1, SRT_CALL, PREMIUM);

  return NULL;
}
