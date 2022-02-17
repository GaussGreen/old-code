/* ==========================================================================
   FILE_NAME:	3FQuantoCalib.c

   PURPOSE:		Modelling of the spot FX vol by taking in consideration:
                                - A foreign LGM2F / domestic LGM1F
                                - A foreign LGM1F / domestic LGM2F
                                - A lognormal dynamics on the Spot FX
                                - A term structure of correlations

                                The calibration is using a 5F approach (2 LGM2F
   , 1 spot FX) and all 3F quanto underlyings must be converted to 2 LGM2F using
   the function convert_3F_Quanto_to_2_LGM2F

                                There are 2 set of functions:

                                A set of functions for a calibration of a 5F
   underlying  , taking a 5x5 correl tensor as input (index details in
   3FQuantoCalib.h)

                                A set of functions for a calibration of a 3F
   quanto underlying  , taking a 4x4 correl tensor as input (index details in
   3FQuantoCalib.h) In this case  , the 4x4 correl tensor is converted to a 5x5
   correl tensor and the 5F calibration routine is called.

   DATE:		26 Sep 2003

   AUTHOR:		J.M.L.
   ========================================================================== */

#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_all3FQuanto.h"
#include "srt_h_allfx3f.h"
#include "srt_h_und_fct.h"

// Functions for exponential integrals
static double Phi_5F(double x, double T, double s, double t) {
  double result;
  result = (exp(-x * (T - t)) - exp(-x * (T - s))) / x;
  return result;
}

static double Ksi_5F(double x, double T, double s, double t) {
  double result;
  result = (t - s - Phi_5F(x, T, s, t)) / x;
  return result;
}

static double Psi_5F(double x, double y, double T, double s, double t) {
  double result;
  result = (t - s - Phi_5F(x, T, s, t) - Phi_5F(y, T, s, t) +
            Phi_5F(x + y, T, s, t)) /
           (x * y);
  return result;
}

// Returns the number of factors of a LGM underlying
Err get_lgm_number_of_factors(char *underlying, int *number_of_factors) {
  Err err = NULL;
  SrtUndPtr und;
  und = lookup_und(underlying);

  // Conducts a few tests
  if (!und) {
    return serror(
        "get_lgm_number_of_factors: Could not find the underlying named %s",
        underlying);
  }
  if (get_underlying_type(und) != INTEREST_RATE_UND) {
    return serror("get_lgm_number_of_factors: Underlying %s is not of type IR",
                  underlying);
  }
  if (get_mdltype_from_irund(und) != LGM) {
    return serror("get_lgm_number_of_factors: Underlying %s is not of type LGM",
                  underlying);
  }

  // Reads the number of factors
  if (get_mdldim_from_irund(und) == ONE_FAC) {
    *number_of_factors = 1;
  } else if (get_mdldim_from_irund(und) == TWO_FAC) {
    *number_of_factors = 2;
  } else {
    return serror("get_lgm_number_of_factors: number of factors of underlying "
                  "%s is undefined",
                  underlying);
  }

  return err;
}

//******************************************************************************************************************
//* *
//*                           3F QUANTO conversion to 5F routines *
//* *
//******************************************************************************************************************

// Filling of the upper right zone of a 5x5 correlation matrix
// See 4x4 and 5x5 index definitions in the 3FQuantoCalib.h
Err convert_4x4_to_5x5_correlation_matrix(double **correl_matrix,
                                          double **correl_matrix_5x5,
                                          int is_lgm2F_dom_for) {

  if (is_lgm2F_dom_for == 1) {
    // LGM2F foreign
    correl_matrix_5x5[W1_D][W2_D] = 1.0;
    correl_matrix_5x5[W1_D][W3_F] = correl_matrix[LGM2F_W1][LGM1F_W3];
    correl_matrix_5x5[W1_D][W4_F] = correl_matrix[LGM2F_W2][LGM1F_W3];
    correl_matrix_5x5[W1_D][FX_5F] = correl_matrix[LGM1F_W3][FX_3F_QUANTO];
    correl_matrix_5x5[W2_D][W3_F] = 0.0;
    correl_matrix_5x5[W2_D][W4_F] = 0.0;
    correl_matrix_5x5[W2_D][FX_5F] = 0.0;
    correl_matrix_5x5[W3_F][W4_F] = correl_matrix[LGM2F_W1][LGM2F_W2];
    correl_matrix_5x5[W3_F][FX_5F] = correl_matrix[LGM2F_W1][FX_3F_QUANTO];
    correl_matrix_5x5[W4_F][FX_5F] = correl_matrix[LGM2F_W2][FX_3F_QUANTO];
  } else if (is_lgm2F_dom_for == 0) {
    // LGM2F domestic
    correl_matrix_5x5[W1_D][W2_D] = correl_matrix[LGM2F_W1][LGM2F_W2];
    correl_matrix_5x5[W1_D][W3_F] = correl_matrix[LGM2F_W1][LGM1F_W3];
    correl_matrix_5x5[W1_D][W4_F] = 0.0;
    ;
    correl_matrix_5x5[W1_D][FX_5F] = correl_matrix[LGM2F_W1][FX_3F_QUANTO];
    correl_matrix_5x5[W2_D][W3_F] = correl_matrix[LGM2F_W2][LGM1F_W3];
    correl_matrix_5x5[W2_D][W4_F] = 0.0;
    correl_matrix_5x5[W2_D][FX_5F] = correl_matrix[LGM2F_W2][FX_3F_QUANTO];
    correl_matrix_5x5[W3_F][W4_F] = 1.0;
    correl_matrix_5x5[W3_F][FX_5F] = correl_matrix[LGM1F_W3][FX_3F_QUANTO];
    correl_matrix_5x5[W4_F][FX_5F] = 0.0;
  } else
    return "lgm2F is neither foreign nor domestic";

  return NULL;
}

// Takes a domestic and a foreign LGM underlying (one being 1F  , the other 2F)
// and returns 2 LGM2F structures
Err convert_3FQuanto_to_2_LGM2F(
    char *dom_und, char *for_und, double **dom_sigma_dates,
    long *n_dom_sigma_dates, double **dom_sigma_values,
    double *dom_fixed_lambda, double *dom_fixed_alpha, double *dom_fixed_gamma,
    double *dom_fixed_rho, double **for_sigma_dates, long *n_for_sigma_dates,
    double **for_sigma_values, double *for_fixed_lambda,
    double *for_fixed_alpha, double *for_fixed_gamma, double *for_fixed_rho,
    int *lgm2F_is_dom_for) {
  Err err = NULL;
  int dom_und_number_of_factors, for_und_number_of_factors;

  double *lgm1F_sigma_dates = NULL, *lgm1F_sigma_values = NULL,
         *lgm1F_tau_dates = NULL, *lgm1F_tau_values = NULL,
         *lgm2F_sigma_dates = NULL, *lgm2F_sigma_values = NULL;

  long n_lgm1F_sigma_dates, n_lgm2F_sigma_dates, n_lgm1F_tau_dates;
  double lgm1F_lambda, lgm2F_tau, alpha, gamma, rho;

  // Checks the LGM underlyings dimensions
  err = get_lgm_number_of_factors(dom_und, &dom_und_number_of_factors);
  if (err)
    goto FREE_RETURN;

  err = get_lgm_number_of_factors(for_und, &for_und_number_of_factors);
  if (err)
    goto FREE_RETURN;

  if (dom_und_number_of_factors == 1 && for_und_number_of_factors == 2) {

    // Reads the DOMESTIC LGM1F
    err = Get_LGM_TermStructure(
        dom_und, &lgm1F_sigma_dates, &lgm1F_sigma_values, &n_lgm1F_sigma_dates,
        &lgm1F_tau_dates, &lgm1F_tau_values, &n_lgm1F_tau_dates);

    if (err)
      goto FREE_RETURN;

    err = get_unique_lambda(lgm1F_tau_values, n_lgm1F_tau_dates, &lgm1F_lambda);

    if (err)
      goto FREE_RETURN;

    // Reads the FOREIGN LGM2F
    *lgm2F_is_dom_for = 1;

    err = Get_LGM2F_TermStructure(for_und, &lgm2F_sigma_dates,
                                  &lgm2F_sigma_values, &n_lgm2F_sigma_dates,
                                  &lgm2F_tau, &alpha, &gamma, &rho);

    if (err)
      goto FREE_RETURN;

    // Setup the proper values for lambda  , alpha  , gamma & correlation matrix
    *dom_fixed_lambda = lgm1F_lambda;
    *dom_fixed_alpha = 0.0;
    *dom_fixed_gamma = 0.0;
    *dom_fixed_rho = 1.0;

    *for_fixed_lambda = 1.0 / lgm2F_tau;
    *for_fixed_alpha = alpha;
    *for_fixed_gamma = gamma;
    *for_fixed_rho = rho;

    // Setup the pointers to the term structures
    *dom_sigma_dates = lgm1F_sigma_dates;
    *n_dom_sigma_dates = n_lgm1F_sigma_dates;
    *dom_sigma_values = lgm1F_sigma_values;

    *for_sigma_dates = lgm2F_sigma_dates;
    *n_for_sigma_dates = n_lgm2F_sigma_dates;
    *for_sigma_values = lgm2F_sigma_values;

  } else if (dom_und_number_of_factors == 2 && for_und_number_of_factors == 1) {

    // Reads the DOMESTIC LGM2F
    *lgm2F_is_dom_for = 0;

    err = Get_LGM2F_TermStructure(dom_und, &lgm2F_sigma_dates,
                                  &lgm2F_sigma_values, &n_lgm2F_sigma_dates,
                                  &lgm2F_tau, &alpha, &gamma, &rho);

    if (err)
      goto FREE_RETURN;

    // Reads the FOREIGN LGM1F
    err = Get_LGM_TermStructure(
        for_und, &lgm1F_sigma_dates, &lgm1F_sigma_values, &n_lgm1F_sigma_dates,
        &lgm1F_tau_dates, &lgm1F_tau_values, &n_lgm1F_tau_dates);

    if (err)
      goto FREE_RETURN;

    err = get_unique_lambda(lgm1F_tau_values, n_lgm1F_tau_dates, &lgm1F_lambda);

    if (err)
      goto FREE_RETURN;

    // Setup the proper values for lambda  , alpha  , gamma & correlation matrix
    *dom_fixed_lambda = 1.0 / lgm2F_tau;
    *dom_fixed_alpha = alpha;
    *dom_fixed_gamma = gamma;
    *dom_fixed_rho = rho;

    *for_fixed_lambda = lgm1F_lambda;
    *for_fixed_alpha = 0.0;
    *for_fixed_gamma = 0.0;
    *for_fixed_rho = 1.0;

    // Setup the pointers to the term structures
    *dom_sigma_dates = lgm2F_sigma_dates;
    *n_dom_sigma_dates = n_lgm2F_sigma_dates;
    *dom_sigma_values = lgm2F_sigma_values;

    *for_sigma_dates = lgm1F_sigma_dates;
    *n_for_sigma_dates = n_lgm1F_sigma_dates;
    *for_sigma_values = lgm1F_sigma_values;
  } else {
    return serror("One LGM should be 2F and the other LGM shoud be 1F");
  }

FREE_RETURN:
  return err;
}

// Takes a 3F Quanto underlying and turns it into one LGM2F and one LGM1F
Err Get_FX_StochRate_TermStructures3FQuanto(
    char *underlying, char **dom_und_name, double **dom_sigma_dates,
    long *n_dom_sigma_dates, double **dom_sigma_values,
    double *dom_fixed_lambda, double *dom_fixed_alpha, double *dom_fixed_gamma,
    double *dom_fixed_rho, char **for_und_name, double **for_sigma_dates,
    long *n_for_sigma_dates, double **for_sigma_values,
    double *for_fixed_lambda, double *for_fixed_alpha, double *for_fixed_gamma,
    double *for_fixed_rho, int *lgm2F_is_dom_for, double **fx_sigma_dates,
    long *n_fx_sigma_dates, double **fx_sigma_values) {
  Err err = NULL;
  SrtUndPtr und;
  char *dom_und, *for_und;
  long today, i, corr_date_n;
  double *corr_date = NULL, *corr = NULL;

  // Prepares the FX curve data
  *fx_sigma_dates = NULL;
  *fx_sigma_values = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Could not find underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", underlying);
    goto FREE_RETURN;
  }

  // Get the type of domestic LGM
  dom_und = get_domname_from_fxund(und);
  strcpy(*dom_und_name, dom_und);

  // Get the type of foreign LGM
  for_und = get_forname_from_fxund(und);
  strcpy(*for_und_name, for_und);

  // Converts 2 underyings (one LGM1F  , the other LGM2F) into 2 LGM2Fs
  err = convert_3FQuanto_to_2_LGM2F(
      dom_und, for_und, dom_sigma_dates, n_dom_sigma_dates, dom_sigma_values,
      dom_fixed_lambda, dom_fixed_alpha, dom_fixed_gamma, dom_fixed_rho,
      for_sigma_dates, n_for_sigma_dates, for_sigma_values, for_fixed_lambda,
      for_fixed_alpha, for_fixed_gamma, for_fixed_rho, lgm2F_is_dom_for);
  if (err)
    goto FREE_RETURN;

  // Now get the FX term structure
  err = srt_f_display_FX_TermStruct(underlying, n_fx_sigma_dates,
                                    fx_sigma_dates, fx_sigma_values,
                                    &corr_date_n, &corr_date, &corr);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(und);
  for (i = 0; i < *n_fx_sigma_dates; i++) {
    (*fx_sigma_dates)[i] = ((*fx_sigma_dates)[i] - today) / 365;
  }

FREE_RETURN:

  if (err) {
    if (*fx_sigma_dates)
      free(*fx_sigma_dates);
    *fx_sigma_dates = NULL;

    if (*fx_sigma_values)
      free(*fx_sigma_values);
    *fx_sigma_values = NULL;
  }

  if (corr_date)
    free(corr_date);
  if (corr)
    free(corr);

  return err;
}

//******************************************************************************************************************
//* *
//*                           Main 5F Calibration core routines *
//* *
//******************************************************************************************************************

// Calculates the coefficients of the quadratic equation
// on an interval [t1  , t2[ where ALL parameters (LGM1F  , LGM2F  , FX sigmas
// , correls and lambdas) are CONSTANT DOMESTIC DYNAMICS: sigma1  , lambda1 +
// sigma2  , lambda2 FOREIGN DYNAMICS: sigma3  , lambda3 + sigma4  , lambda4
// Reminder: in the double** 5x5 correl_matrix  , the indices are:
// 0: DOM W1  , 1: DOM W2  , 3: FOR W3  , 4: FOR W4  , 5: FX
Err Coefs_Partial_Var_5F(double Tk, double t1, double t2, double sigma1,
                         double lambda1, double sigma2, double lambda2,
                         double sigma3, double lambda3, double sigma4,
                         double lambda4, double **correl_matrix, double *a,
                         double *b, double *c) {
  (*a) = t2 - t1;

  (*b) =
      -2 * correl_matrix[W3_F][FX_5F] * sigma3 * Ksi_5F(lambda3, Tk, t1, t2) -
      2 * correl_matrix[W4_F][FX_5F] * sigma4 * Ksi_5F(lambda4, Tk, t1, t2) +
      2 * correl_matrix[W1_D][FX_5F] * sigma1 * Ksi_5F(lambda1, Tk, t1, t2) +
      2 * correl_matrix[W2_D][FX_5F] * sigma2 * Ksi_5F(lambda2, Tk, t1, t2);

  (*c) = sigma3 * sigma3 * Psi_5F(lambda3, lambda3, Tk, t1, t2) +
         sigma4 * sigma4 * Psi_5F(lambda4, lambda4, Tk, t1, t2) +
         sigma1 * sigma1 * Psi_5F(lambda1, lambda1, Tk, t1, t2) +
         sigma2 * sigma2 * Psi_5F(lambda2, lambda2, Tk, t1, t2) +
         2 * correl_matrix[W3_F][W4_F] * sigma3 * sigma4 *
             Psi_5F(lambda3, lambda4, Tk, t1, t2) -
         2 * correl_matrix[W1_D][W3_F] * sigma1 * sigma3 *
             Psi_5F(lambda1, lambda3, Tk, t1, t2) -
         2 * correl_matrix[W2_D][W3_F] * sigma2 * sigma3 *
             Psi_5F(lambda2, lambda3, Tk, t1, t2) -
         2 * correl_matrix[W1_D][W4_F] * sigma1 * sigma4 *
             Psi_5F(lambda1, lambda4, Tk, t1, t2) -
         2 * correl_matrix[W2_D][W4_F] * sigma2 * sigma4 *
             Psi_5F(lambda2, lambda4, Tk, t1, t2) +
         2 * correl_matrix[W1_D][W2_D] * sigma1 * sigma2 *
             Psi_5F(lambda1, lambda2, Tk, t1, t2);

  return NULL;
}

// Calculates the partial volatility of a FX fwd maturing at time
// fx_option_maturity on an interval [t1  , t2[ where ALL parameters (LGM1F  ,
// LGM2F  , FX sigmas  , correls and lambdas) are CONSTANT
Err Partial_Var_5F(double fx_option_maturity, double t1, double t2,
                   double dom_sigma1, double dom_lambda1, double dom_sigma2,
                   double dom_lambda2, double for_sigma3, double for_lambda3,
                   double for_sigma4, double for_lambda4, double fx_sigma,
                   double **correl_matrix, double *var_partial) {
  Err err = NULL;
  double a, b, c;

  // Get coefficients
  err =
      Coefs_Partial_Var_5F(fx_option_maturity, t1, t2, dom_sigma1, dom_lambda1,
                           dom_sigma2, dom_lambda2, for_sigma3, for_lambda3,
                           for_sigma4, for_lambda4, correl_matrix, &a, &b, &c);

  (*var_partial) = fx_sigma * (a * fx_sigma + b) + c;

  return NULL;
}

// Calibration of a FX term structure to a set of European FX options
// given a TERM STRUCTURE OF CORRELATIONS 5x5
// We know the partial vol up to T1 and we bootstrap on the spot fx volatility
// between [T1  , T2] using the FX option maturing at T2
Err Fx5FtsCalibrationCorr(
    double *fx_options_exo_times, double *fx_options_mat_times,
    double *fx_options_vols, long n_fx_options, double *ir_dates,
    long n_ir_dates, double *dom_sigma_values, double dom_fixed_lambda,
    double dom_fixed_alpha, double dom_fixed_gamma, double *for_sigma_values,
    double for_fixed_lambda, double for_fixed_alpha, double for_fixed_gamma,
    double *correl_matrix_dates, long n_correl_matrix_dates,
    double ***correl_matrix_ts_5x5, double **fx_vol_curve) {
  Err err = NULL;
  double a_part, b_part, c_part, a, b, c, c2, T1, T2, t1, t2, tt1, tt2, vol,
      cum_var, Tk, delta;
  double lambda1, lambda2, lambda3, lambda4, sigma1, sigma2, sigma3, sigma4,
      VolMin;
  long i, j, i_fx_options, start_index_ir, end_index_ir, start_index_corr,
      end_index_corr;

  // Initializes the result
  (*fx_vol_curve) = NULL;
  (*fx_vol_curve) = (double *)calloc(n_fx_options, sizeof(double));
  if (!(*fx_vol_curve)) {
    err = "Memory allocation error in FX5FtsCalibration";
    goto FREE_RETURN;
  }

  // Setup the mean reversions
  lambda1 = dom_fixed_lambda;
  lambda2 = lambda1 + dom_fixed_gamma;
  lambda3 = for_fixed_lambda;
  lambda4 = lambda3 + for_fixed_gamma;

  // Perform the bootstrap through all FX Options
  for (i_fx_options = 0; i_fx_options < n_fx_options; i_fx_options++) {
    T2 = fx_options_exo_times[i_fx_options];
    Tk = fx_options_mat_times[i_fx_options];

    // Get the cumulated variance up to T1
    if (i_fx_options > 0) {
      T1 = fx_options_exo_times[i_fx_options - 1];

      err = Fx5FImpliedVolCorr(
          Tk, 0, T1, ir_dates, n_ir_dates, dom_sigma_values, dom_fixed_lambda,
          dom_fixed_alpha, dom_fixed_gamma, for_sigma_values, for_fixed_lambda,
          for_fixed_alpha, for_fixed_gamma, fx_options_exo_times, *fx_vol_curve,
          n_fx_options, correl_matrix_dates, n_correl_matrix_dates,
          correl_matrix_ts_5x5, &vol);

      if (err)
        goto FREE_RETURN;

      cum_var = vol * vol * T1;
    } else {
      T1 = 0;
      cum_var = 0;
    }

    // Calculates the coefficients of the quadratic equation expressing the
    // variance between T1 and T2 We know that fx_sigma will be constant between
    // T1 and T2 but we still need to split by sections [t1  , t2] where LGM
    //sigmas are constants

    a = b = c2 = c = 0;
    start_index_ir = Get_Index(T1, ir_dates, n_ir_dates);
    end_index_ir = Get_Index(T2, ir_dates, n_ir_dates);

    for (i = start_index_ir; i < end_index_ir + 1; i++) {
      if (i > start_index_ir) {
        t1 = ir_dates[i - 1];
      } else {
        t1 = T1;
      }

      if (i == end_index_ir || start_index_ir == end_index_ir) {
        t2 = T2;
      } else {
        t2 = ir_dates[i];
      }

      sigma1 = dom_sigma_values[i];
      sigma2 = sigma1 * dom_fixed_alpha;
      sigma3 = for_sigma_values[i];
      sigma4 = sigma3 * for_fixed_alpha;

      start_index_corr =
          Get_Index(t1, correl_matrix_dates, n_correl_matrix_dates);
      end_index_corr =
          Get_Index(t2, correl_matrix_dates, n_correl_matrix_dates);

      for (j = start_index_corr; j < end_index_corr + 1; j++) {
        if (j > start_index_corr) {
          tt1 = correl_matrix_dates[j - 1];
        } else {
          tt1 = t1;
        }

        if (j == end_index_corr || start_index_corr == end_index_corr) {
          tt2 = t2;
        } else {
          tt2 = correl_matrix_dates[j];
        }

        err = Coefs_Partial_Var_5F(Tk, tt1, tt2, sigma1, lambda1, sigma2,
                                   lambda2, sigma3, lambda3, sigma4, lambda4,
                                   correl_matrix_ts_5x5[j], &a_part, &b_part,
                                   &c_part);

        if (err)
          goto FREE_RETURN;

        // We then sum the partial coefficients calculated on [tt1  , tt2]
        a += a_part;
        b += b_part;
        c2 += c_part;

      } // End loop on correlation

    } // End loop on rates

    c = c2 + cum_var; // Add the cumulative variance calculated on [0  , T1]
    c -= fx_options_vols[i_fx_options] * fx_options_vols[i_fx_options] * T2;

    // Now we solve the quadratic equation aX^2+bX+c
    delta = b * b - 4 * a * c;

    // delta<0 or X^2-SX+P=0 <=> P>0 & S<0 <=> 2 negative solutions
    if ((delta < 0) || ((c > 0) && (b > 0))) {
      VolMin = sqrt((c2 + cum_var - b * b / (4 * a)) / Tk) * 100;
      err = serror("Cannot find solution to match the option %d (its vol "
                   "should be > %.2f %%)",
                   i_fx_options, VolMin);
      goto FREE_RETURN;
    }

    // if 2 positive solutions (P>0  , S>0) take the smallest one
    if ((c > 0) && (b < 0)) {
      (*fx_vol_curve)[i_fx_options] = (-b - sqrt(delta)) / (2 * a);
    } else {
      (*fx_vol_curve)[i_fx_options] = (-b + sqrt(delta)) / (2 * a);
    }
  }

FREE_RETURN:
  if (err) {
    if (*fx_vol_curve) {
      free(*fx_vol_curve);
      *fx_vol_curve = NULL;
    }
  }

  return err;
}

// Calculates the partial implied vol between start_date and end_date of
// the fwd FX of time fx_option_maturity given the term structures of 2 LGM2Fs
// and a TERM STRUCTURE OF CORRELATIONS
// The TS of rates must be merged before and the correl matrix MUST be 5x5
Err Fx5FImpliedVolCorr(double fx_option_maturity, double start_date,
                       double end_date, double *ir_dates, long n_ir_dates,
                       double *dom_sigma_values, double dom_fixed_lambda,
                       double dom_fixed_alpha, double dom_fixed_gamma,
                       double *for_sigma_values, double for_fixed_lambda,
                       double for_fixed_alpha, double for_fixed_gamma,
                       double *fx_sigma_times, double *fx_sigma_values,
                       long n_fx_sigma_times, double *correl_matrix_dates,
                       long n_correl_matrix_dates,
                       double ***correl_matrix_5x5_ts, double *fx_vol) {
  Err err = NULL;
  long start_index_fx, end_index_fx, start_index_ir, end_index_ir,
      start_index_corr, end_index_corr, i, j, k;
  double var, dom_fixed_lambda1, dom_fixed_lambda2, for_fixed_lambda1,
      for_fixed_lambda2;
  double dom_sigma1, dom_sigma2, for_sigma1, for_sigma2, fx_sigma, var_partial;
  double T1, T2, t1, t2, tt1, tt2;

  if (start_date > end_date) {
    err = "end_date before start_date in Fx5FImpliedVolCorr!";
    goto FREE_RETURN;
  }

  if (end_date == 0) {
    (*fx_vol) = 0;
    goto FREE_RETURN;
  }

  // Reads the constant LGM2F lambdas
  dom_fixed_lambda1 = dom_fixed_lambda;
  dom_fixed_lambda2 = dom_fixed_lambda1 + dom_fixed_gamma;
  for_fixed_lambda1 = for_fixed_lambda;
  for_fixed_lambda2 = for_fixed_lambda1 + for_fixed_gamma;

  // Integrating from start_time to end_time is equivalent to integrating
  // the FX vols from start_index to end_index
  start_index_fx = Get_Index(start_date, fx_sigma_times, n_fx_sigma_times);
  end_index_fx = Get_Index(end_date, fx_sigma_times, n_fx_sigma_times);

  var = 0;

  // For all ranges of fx vols on which the spot FX vol is constant  ,
  //(ranges between start_index and end_index)
  for (i = start_index_fx; i < end_index_fx + 1; i++) {
    // We first isolate the [T1  , T2] corresponding to the current range
    // on which the FX vols are constant = fx_sigma_values[i]
    if (i > start_index_fx) {
      T1 = fx_sigma_times[i - 1];
    } else {
      // First bit of the TS
      T1 = start_date;
    }

    if (i == end_index_fx || start_index_fx == end_index_fx) {
      T2 = end_date;
    } else {
      T2 = fx_sigma_times[i];
    }
    fx_sigma = fx_sigma_values[i];

    // We then map the corresponding indices of the LGMs term structures
    // That will correspond to a sum of intervals [t1  , t2] with the leftmost
    // t1 = T1 and the rightmost t2 = T2
    start_index_ir = Get_Index(T1, ir_dates, n_ir_dates);
    end_index_ir = Get_Index(T2, ir_dates, n_ir_dates);

    for (j = start_index_ir; j < end_index_ir + 1; j++) {
      if (j > start_index_ir) {
        t1 = ir_dates[j - 1];
      } else {
        t1 = T1;
      }

      if (j == end_index_ir || start_index_ir == end_index_ir) {
        t2 = T2;
      } else {
        t2 = ir_dates[j];
      }

      // On the interval [t1  , t2]  , the LGM sigmas are constants
      dom_sigma1 = dom_sigma_values[j];
      dom_sigma2 = dom_sigma1 * dom_fixed_alpha;
      for_sigma1 = for_sigma_values[j];
      for_sigma2 = for_sigma1 * for_fixed_alpha;

      // We now slice for correlations
      start_index_corr =
          Get_Index(t1, correl_matrix_dates, n_correl_matrix_dates);
      end_index_corr =
          Get_Index(t2, correl_matrix_dates, n_correl_matrix_dates);

      for (k = start_index_corr; k < end_index_corr + 1; k++) {
        if (k > start_index_corr) {
          tt1 = correl_matrix_dates[k - 1];
        } else {
          tt1 = t1;
        }

        if (k == end_index_corr || start_index_corr == end_index_corr) {
          tt2 = t2;
        } else {
          tt2 = correl_matrix_dates[k];
        }

        err = Partial_Var_5F(fx_option_maturity, tt1, tt2, dom_sigma1,
                             dom_fixed_lambda1, dom_sigma2, dom_fixed_lambda2,
                             for_sigma1, for_fixed_lambda1, for_sigma2,
                             for_fixed_lambda2, fx_sigma,
                             correl_matrix_5x5_ts[k], &var_partial);

        if (err)
          goto FREE_RETURN;

        var += var_partial;

      } // End loop on correlation
    }   // End loop on IR
  }     // End loop on FX Vols

  if (fabs(end_date - start_date) > 1.0e-08) {
    *fx_vol = sqrt(var / (end_date - start_date));
  } else {
    *fx_vol = 0.0;
  }

FREE_RETURN:
  return err;
}

//******************************************************************************************************************
//* *
//*                           5F Entry Functions *
//* *
//******************************************************************************************************************

// FX calibration directly from the underlying 5F (2 LGM2F)
// given a TERM STRUCTURE OF CORRELATIONS
Err Fx5FCalibrationCorr(char *dom_und, char *for_und,
                        double *correl_matrix_dates, double ***correl_matrix_ts,
                        long n_correl_matrix_ts, double *fx_options_exo_times,
                        double *fx_options_mat_times, double *fx_options_vols,
                        long n_fx_options, double **fx_vol_curve) {

  Err err = NULL;
  long n_for_sigma_dates, n_dom_sigma_dates, n_merge_dates;
  double dom_fixed_tau, dom_fixed_lambda, dom_fixed_alpha, dom_fixed_gamma,
      dom_fixed_rho;
  double for_fixed_tau, for_fixed_lambda, for_fixed_alpha, for_fixed_gamma,
      for_fixed_rho;
  int i, lgm2F_is_dom_for = 1; // 0: dom 1: for

  double *dom_sigma_dates = NULL, *dom_sigma_values = NULL,
         *for_sigma_dates = NULL, *for_sigma_values = NULL, *merge_dates = NULL,
         *dom_merge_sigma_values = NULL, *for_merge_sigma_values = NULL;

  // Reads the 2 LGM2F underlyings
  err = Get_LGM2F_TermStructure(
      dom_und, &dom_sigma_dates, &dom_sigma_values, &n_dom_sigma_dates,
      &dom_fixed_tau, &dom_fixed_alpha, &dom_fixed_gamma, &dom_fixed_rho);

  if (err)
    goto FREE_RETURN;

  dom_fixed_lambda = 1.0 / dom_fixed_tau;

  err = Get_LGM2F_TermStructure(
      for_und, &for_sigma_dates, &for_sigma_values, &n_for_sigma_dates,
      &for_fixed_tau, &for_fixed_alpha, &for_fixed_gamma, &for_fixed_rho);

  if (err)
    goto FREE_RETURN;

  for_fixed_lambda = 1.0 / for_fixed_tau;

  if (err)
    goto FREE_RETURN;

  // Merge all term structures to have consistent dates + matrix_5x5!!!
  err = merge_rates_ts(dom_sigma_dates, dom_sigma_values, n_dom_sigma_dates,
                       for_sigma_dates, for_sigma_values, n_for_sigma_dates,
                       &merge_dates, &dom_merge_sigma_values,
                       &for_merge_sigma_values, &n_merge_dates);

  if (err)
    goto FREE_RETURN;

  // Writes the LGM2F correls into the 5x5 correlation matrix
  for (i = 0; i < n_correl_matrix_ts; i++) {
    correl_matrix_ts[i][W1_D][W2_D] = dom_fixed_rho;
    correl_matrix_ts[i][W3_F][W4_F] = for_fixed_rho;
  }

  // Calibrate now....
  err = Fx5FtsCalibrationCorr(
      fx_options_exo_times, fx_options_mat_times, fx_options_vols, n_fx_options,
      merge_dates, n_merge_dates, dom_merge_sigma_values, dom_fixed_lambda,
      dom_fixed_alpha, dom_fixed_gamma, for_merge_sigma_values,
      for_fixed_lambda, for_fixed_alpha, for_fixed_gamma, correl_matrix_dates,
      n_correl_matrix_ts, correl_matrix_ts, fx_vol_curve);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (dom_sigma_dates)
    free(dom_sigma_dates);
  if (dom_sigma_values)
    free(dom_sigma_values);
  if (for_sigma_dates)
    free(for_sigma_dates);
  if (for_sigma_values)
    free(for_sigma_values);
  if (merge_dates)
    free(merge_dates);
  if (dom_merge_sigma_values)
    free(dom_merge_sigma_values);
  if (for_merge_sigma_values)
    free(for_merge_sigma_values);
  return err;
}

// Gets the implied FX fwd vol from a 5F underlying
// given a TERM STRUCTURE OF CORRELATIONS
Err Fx5FImpliedVolFrom5FUnderlying(char *und_5F, double *correl_matrix_dates,
                                   double ***correl_matrix_ts,
                                   long n_correl_matrix_ts, double val_time,
                                   double start_time, double fix_time,
                                   double *fx_vol) {
  Err err = NULL;
  int i;
  double *dom_sigma_dates = NULL, *dom_sigma_values = NULL,
         *for_sigma_dates = NULL, *for_sigma_values = NULL,
         *fx_sigma_dates = NULL, *fx_sigma_values = NULL, *merge_dates = NULL,
         *dom_merge_sigma_values = NULL, *for_merge_sigma_values = NULL;

  long n_dom_sigma_dates, n_for_sigma_dates, n_fx_sigma_dates, n_merge_dates;
  double dom_fixed_tau, dom_fixed_lambda, dom_fixed_alpha, dom_fixed_gamma,
      dom_fixed_rho;
  double for_fixed_tau, for_fixed_lambda, for_fixed_alpha, for_fixed_gamma,
      for_fixed_rho;
  // Dummy correlations:
  double correl_dom_for, correl_dom_fx, correl_for_fx;

  // Get the 2 underlying LGM2Fs
  err = Get_FX_StochRate_TermStructures5F(
      und_5F, &dom_sigma_dates, &dom_sigma_values, &n_dom_sigma_dates,
      &dom_fixed_tau, &dom_fixed_alpha, &dom_fixed_gamma, &dom_fixed_rho,
      &for_sigma_dates, &for_sigma_values, &n_for_sigma_dates, &for_fixed_tau,
      &for_fixed_alpha, &for_fixed_gamma, &for_fixed_rho, &fx_sigma_dates,
      &fx_sigma_values, &n_fx_sigma_dates, &correl_dom_for, &correl_dom_fx,
      &correl_for_fx);

  if (err)
    goto FREE_RETURN;

  dom_fixed_lambda = 1.0 / dom_fixed_tau;
  for_fixed_lambda = 1.0 / for_fixed_tau;

  // Merge all term structures to have consistent dates
  err = merge_rates_ts(dom_sigma_dates, dom_sigma_values, n_dom_sigma_dates,
                       for_sigma_dates, for_sigma_values, n_for_sigma_dates,
                       &merge_dates, &dom_merge_sigma_values,
                       &for_merge_sigma_values, &n_merge_dates);

  if (err)
    goto FREE_RETURN;

  // Writes the LGM2F correls into the 5x5 correlation matrix
  for (i = 0; i < n_correl_matrix_ts; i++) {
    correl_matrix_ts[i][W1_D][W2_D] = dom_fixed_rho;
    correl_matrix_ts[i][W3_F][W4_F] = for_fixed_rho;
  }

  err = Fx5FImpliedVolCorr(
      val_time, start_time, fix_time, merge_dates, n_merge_dates,
      dom_merge_sigma_values, dom_fixed_lambda, dom_fixed_alpha,
      dom_fixed_gamma, for_merge_sigma_values, for_fixed_lambda,
      for_fixed_alpha, for_fixed_gamma, fx_sigma_dates, fx_sigma_values,
      n_fx_sigma_dates, correl_matrix_dates, n_correl_matrix_ts,
      correl_matrix_ts, fx_vol);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (dom_sigma_dates)
    free(dom_sigma_dates);
  if (dom_sigma_values)
    free(dom_sigma_values);
  if (for_sigma_dates)
    free(for_sigma_dates);
  if (for_sigma_values)
    free(for_sigma_values);
  if (fx_sigma_dates)
    free(fx_sigma_dates);
  if (fx_sigma_values)
    free(fx_sigma_values);
  if (merge_dates)
    free(merge_dates);
  if (dom_merge_sigma_values)
    free(dom_merge_sigma_values);
  if (dom_merge_sigma_values)
    free(for_merge_sigma_values);

  return err;
}

//******************************************************************************************************************
//* *
//*                           3F QUANTO Entry Functions *
//* *
//******************************************************************************************************************

// FX calibration directly from the underlying
// given a TERM STRUCTURE OF CORRELATIONS
Err Fx3FQuantoCalibrationCorr(
    char *dom_und, char *for_und, double *correl_matrix_dates,
    double ***correl_matrix_ts, long n_correl_matrix_ts,
    double *fx_options_exo_times, double *fx_options_mat_times,
    double *fx_options_vols, long n_fx_options, double **fx_vol_curve) {

  Err err = NULL;
  long n_for_sigma_dates, n_dom_sigma_dates, n_merge_dates;
  double dom_fixed_lambda, dom_fixed_alpha, dom_fixed_gamma, dom_fixed_rho;
  double for_fixed_lambda, for_fixed_alpha, for_fixed_gamma, for_fixed_rho;
  int lgm2F_is_dom_for = 1, i; // 0: dom 1: for

  double *dom_sigma_dates = NULL, *dom_sigma_values = NULL,
         *for_sigma_dates = NULL, *for_sigma_values = NULL, *merge_dates = NULL,
         *dom_merge_sigma_values = NULL, *for_merge_sigma_values = NULL;

  double ***correl_matrix_ts_5x5 = NULL;

  // Setup the new 5x5 correlation matrix  , note that we put the depth in first
  // so that we can extract a 5x5 slice from it and pass it to our double**
  // functions.
  correl_matrix_ts_5x5 = f3tensor(0, n_correl_matrix_ts - 1, 0, 4, 0, 4);

  // Converts 2 underlyings (one LGM1F  , the other LGM2F) into 2 LGM2Fs
  err = convert_3FQuanto_to_2_LGM2F(
      dom_und, for_und, &dom_sigma_dates, &n_dom_sigma_dates, &dom_sigma_values,
      &dom_fixed_lambda, &dom_fixed_alpha, &dom_fixed_gamma, &dom_fixed_rho,
      &for_sigma_dates, &n_for_sigma_dates, &for_sigma_values,
      &for_fixed_lambda, &for_fixed_alpha, &for_fixed_gamma, &for_fixed_rho,
      &lgm2F_is_dom_for);

  if (err)
    goto FREE_RETURN;

  // Initializes the new correlation matrix 5x5
  for (i = 0; i < n_correl_matrix_ts; i++) {
    if (lgm2F_is_dom_for == 0)
      correl_matrix_ts[i][LGM2F_W1][LGM2F_W2] = dom_fixed_rho;
    else
      correl_matrix_ts[i][LGM2F_W1][LGM2F_W2] = for_fixed_rho;

    err = convert_4x4_to_5x5_correlation_matrix(
        correl_matrix_ts[i], correl_matrix_ts_5x5[i], lgm2F_is_dom_for);

    if (err)
      goto FREE_RETURN;
  }

  // Merge all term structures to have consistent dates + matrix_5x5!!!
  err = merge_rates_ts(dom_sigma_dates, dom_sigma_values, n_dom_sigma_dates,
                       for_sigma_dates, for_sigma_values, n_for_sigma_dates,
                       &merge_dates, &dom_merge_sigma_values,
                       &for_merge_sigma_values, &n_merge_dates);

  if (err)
    goto FREE_RETURN;

  // Calibrate now....
  err = Fx5FtsCalibrationCorr(
      fx_options_exo_times, fx_options_mat_times, fx_options_vols, n_fx_options,
      merge_dates, n_merge_dates, dom_merge_sigma_values, dom_fixed_lambda,
      dom_fixed_alpha, dom_fixed_gamma, for_merge_sigma_values,
      for_fixed_lambda, for_fixed_alpha, for_fixed_gamma, correl_matrix_dates,
      n_correl_matrix_ts, correl_matrix_ts_5x5, fx_vol_curve);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (dom_sigma_dates)
    free(dom_sigma_dates);
  if (dom_sigma_values)
    free(dom_sigma_values);
  if (for_sigma_dates)
    free(for_sigma_dates);
  if (for_sigma_values)
    free(for_sigma_values);
  if (merge_dates)
    free(merge_dates);
  if (dom_merge_sigma_values)
    free(dom_merge_sigma_values);
  if (for_merge_sigma_values)
    free(for_merge_sigma_values);
  if (correl_matrix_ts_5x5)
    free_f3tensor(correl_matrix_ts_5x5, 0, n_correl_matrix_ts - 1, 0, 4, 0, 4);
  return err;
}

// Gets the implied FX fwd vol from a 3F Quanto underlying
Err Fx3FQuantoImpliedVolFrom3FQuantoUnderlying(
    char *und_3F_quanto, double *correl_matrix_dates,
    double ***correl_matrix_ts, long n_correl_matrix_ts, double val_time,
    double start_time, double fix_time, double *fx_vol) {
  Err err = NULL;

  double *dom_sigma_dates = NULL, *dom_sigma_values = NULL,
         *for_sigma_dates = NULL, *for_sigma_values = NULL,
         *fx_sigma_dates = NULL, *fx_sigma_values = NULL, *merge_dates = NULL,
         *dom_merge_sigma_values = NULL, *for_merge_sigma_values = NULL;

  double ***correl_matrix_ts_5x5 = NULL;

  char *dom_name = NULL, *for_name = NULL;

  int i, lgm2F_is_dom_for;
  long n_dom_sigma_dates, n_for_sigma_dates, n_fx_sigma_dates, n_merge_dates;
  double dom_fixed_lambda, dom_fixed_alpha, dom_fixed_gamma, dom_fixed_rho;
  double for_fixed_lambda, for_fixed_alpha, for_fixed_gamma, for_fixed_rho;

  dom_name = (char *)calloc(256, sizeof(char));
  for_name = (char *)calloc(256, sizeof(char));

  err = Get_FX_StochRate_TermStructures3FQuanto(
      und_3F_quanto, &dom_name, &dom_sigma_dates, &n_dom_sigma_dates,
      &dom_sigma_values, &dom_fixed_lambda, &dom_fixed_alpha, &dom_fixed_gamma,
      &dom_fixed_rho, &for_name, &for_sigma_dates, &n_for_sigma_dates,
      &for_sigma_values, &for_fixed_lambda, &for_fixed_alpha, &for_fixed_gamma,
      &for_fixed_rho, &lgm2F_is_dom_for, &fx_sigma_dates, &n_fx_sigma_dates,
      &fx_sigma_values);

  if (err)
    goto FREE_RETURN;

  // Merge all term structures to have consistent dates
  err = merge_rates_ts(dom_sigma_dates, dom_sigma_values, n_dom_sigma_dates,
                       for_sigma_dates, for_sigma_values, n_for_sigma_dates,
                       &merge_dates, &dom_merge_sigma_values,
                       &for_merge_sigma_values, &n_merge_dates);

  if (err)
    goto FREE_RETURN;

  // Setup the new 5x5 correlation matrix  , note that we put the depth in first
  // so that we can extract a 5x5 slice from it and pass it to our double**
  // functions.
  correl_matrix_ts_5x5 = f3tensor(0, n_correl_matrix_ts - 1, 0, 4, 0, 4);

  // Initializes the new correlation matrix 5x5
  for (i = 0; i < n_correl_matrix_ts; i++) {
    if (lgm2F_is_dom_for == 0)
      correl_matrix_ts[i][LGM2F_W1][LGM2F_W2] = dom_fixed_rho;
    else
      correl_matrix_ts[i][LGM2F_W1][LGM2F_W2] = for_fixed_rho;

    err = convert_4x4_to_5x5_correlation_matrix(
        correl_matrix_ts[i], correl_matrix_ts_5x5[i], lgm2F_is_dom_for);
    if (err)
      goto FREE_RETURN;
  }

  err = Fx5FImpliedVolCorr(
      val_time, start_time, fix_time, merge_dates, n_merge_dates,
      dom_merge_sigma_values, dom_fixed_lambda, dom_fixed_alpha,
      dom_fixed_gamma, for_merge_sigma_values, for_fixed_lambda,
      for_fixed_alpha, for_fixed_gamma, fx_sigma_dates, fx_sigma_values,
      n_fx_sigma_dates, correl_matrix_dates, n_correl_matrix_ts,
      correl_matrix_ts_5x5, fx_vol);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (dom_sigma_dates)
    free(dom_sigma_dates);
  if (dom_sigma_values)
    free(dom_sigma_values);
  if (for_sigma_dates)
    free(for_sigma_dates);
  if (for_sigma_values)
    free(for_sigma_values);
  if (fx_sigma_dates)
    free(fx_sigma_dates);
  if (fx_sigma_values)
    free(fx_sigma_values);
  if (merge_dates)
    free(merge_dates);
  if (dom_merge_sigma_values)
    free(dom_merge_sigma_values);
  if (dom_merge_sigma_values)
    free(for_merge_sigma_values);
  if (correl_matrix_ts_5x5)
    free_f3tensor(correl_matrix_ts_5x5, 0, n_correl_matrix_ts - 1, 0, 4, 0, 4);
  if (dom_name)
    free(dom_name);
  if (for_name)
    free(for_name);

  return err;
}
