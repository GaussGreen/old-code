#include "FXLGMSVUnd.h"
#include "Fx3FUtils.h"
#include "LGMSVPDE.h"
#include "LGMSVUtil.h"
#include "opfnctns.h"

#define HESTON_FIRSTBREAK 10
#define HESTON_LASTBREAK 1000000
#define PI2 6.28318530717958647692528676655900576839433879875021164194988918462
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

Err init_LGMSV_str(long today,
                   /* TS Dates and Times */
                   int num_dates, long *dates, double *times,

                   /* Domestic Underlying */
                   int one2F,

                   double *sigma,

                   double *tau,

                   double alpha, double gamma, double rho,

                   double *alphaSV, double *lamSV, double *rhoSV,
                   double *rho2SV,

                   double tstar,

                   irm_sv *irmsv) {
  Err err = NULL;

  irmsv->one2F = one2F;
  irmsv->alpha = alpha;
  irmsv->gamma = gamma;
  irmsv->rho = rho;
  irmsv->tstar = tstar;

  irmsv->sigma_dates = NULL;
  irmsv->sigma_times = NULL;

  irmsv->tau_dates = NULL;
  irmsv->tau_times = NULL;

  irmsv->smile_dates = NULL;
  irmsv->smile_times = NULL;

  irmsv->sigma = NULL;
  irmsv->tau = NULL;

  irmsv->alphaSV = NULL;
  irmsv->lamSV = NULL;
  irmsv->rhoSV = NULL;
  irmsv->rho2SV = NULL;

  irmsv->num_sigma = num_dates;
  irmsv->sigma_dates = (long *)calloc(num_dates, sizeof(long));
  irmsv->sigma_times = (double *)calloc(num_dates, sizeof(double));

  irmsv->num_tau = num_dates;
  irmsv->tau_dates = (long *)calloc(num_dates, sizeof(long));
  irmsv->tau_times = (double *)calloc(num_dates, sizeof(double));

  irmsv->num_smile = num_dates;
  irmsv->smile_dates = (long *)calloc(num_dates, sizeof(long));
  irmsv->smile_times = (double *)calloc(num_dates, sizeof(double));

  irmsv->sigma = (double *)calloc(num_dates, sizeof(double));
  irmsv->tau = (double *)calloc(num_dates, sizeof(double));
  irmsv->alphaSV = (double *)calloc(num_dates, sizeof(double));
  irmsv->lamSV = (double *)calloc(num_dates, sizeof(double));
  irmsv->rhoSV = (double *)calloc(num_dates, sizeof(double));
  irmsv->rho2SV = (double *)calloc(num_dates, sizeof(double));

  if ((!irmsv->sigma_dates) || (!irmsv->sigma_times) || (!irmsv->tau_dates) ||
      (!irmsv->tau_times) || (!irmsv->smile_dates) || (!irmsv->smile_times) ||
      (!irmsv->sigma) || (!irmsv->tau) || (!irmsv->alphaSV) ||
      (!irmsv->lamSV) || (!irmsv->rhoSV) || (!irmsv->rho2SV)) {
    err = "memory allocation failed in init_irmsvUnd";
    irm_sv_free_struct(irmsv);
  }

  memcpy(irmsv->sigma_dates, dates, num_dates * sizeof(long));
  memcpy(irmsv->sigma_times, times, num_dates * sizeof(double));

  memcpy(irmsv->tau_dates, dates, num_dates * sizeof(long));
  memcpy(irmsv->tau_times, times, num_dates * sizeof(double));

  memcpy(irmsv->smile_dates, dates, num_dates * sizeof(long));
  memcpy(irmsv->smile_times, times, num_dates * sizeof(double));

  memcpy(irmsv->sigma, sigma, num_dates * sizeof(double));
  memcpy(irmsv->tau, tau, num_dates * sizeof(double));
  memcpy(irmsv->alphaSV, alphaSV, num_dates * sizeof(double));
  memcpy(irmsv->lamSV, lamSV, num_dates * sizeof(double));
  memcpy(irmsv->rhoSV, rhoSV, num_dates * sizeof(double));
  memcpy(irmsv->rho2SV, rho2SV, num_dates * sizeof(double));

  return err;
}

Err fxlgmsv_free_und(fxlgmsv_str *fxlgmsv) {
  if (fxlgmsv) {
    if (fxlgmsv->dates)
      free(fxlgmsv->dates);
    if (fxlgmsv->times)
      free(fxlgmsv->times);

    if (fxlgmsv->fx_sigma)
      free(fxlgmsv->fx_sigma);

    if (fxlgmsv->correlation)
      free_f3tensor(fxlgmsv->correlation, 0, fxlgmsv->num_dates - 1, 0, 6, 0,
                    6);

    if (fxlgmsv->dom_irmsv) {
      irm_sv_free_struct(fxlgmsv->dom_irmsv);
      free(fxlgmsv->dom_irmsv);
      fxlgmsv->dom_irmsv = NULL;
    }

    if (fxlgmsv->for_irmsv) {
      irm_sv_free_struct(fxlgmsv->for_irmsv);
      free(fxlgmsv->for_irmsv);
      fxlgmsv->for_irmsv = NULL;
    }

    fxlgmsv = NULL;
  }

  return NULL;
}

Err init_FXLGMSVUnd(long today,
                    /* TS Dates and Times */
                    int num_dates, long *dates, double *times,

                    /* Domestic Underlying */
                    int dom_one2F,

                    double *dom_sigma,

                    double *dom_tau,

                    double dom_alpha, double dom_gamma, double dom_rho,

                    double *dom_alphaSV, double *dom_lamSV, double *dom_rhoSV,
                    double *dom_rho2SV,

                    double dom_tstar,

                    /* Foreign Underlying */
                    int for_one2F,

                    double *for_sigma,

                    double *for_tau,

                    double for_alpha, double for_gamma, double for_rho,

                    double *for_alphaSV, double *for_lamSV, double *for_rhoSV,
                    double *for_rho2SV,

                    double for_tstar,

                    /* FX Underlying */
                    double fx_spot,

                    double *fx_sigma,

                    /* Correlations */
                    double ***correlation,

                    fxlgmsv_str *fxlgmsv) {
  Err err = NULL;
  int i, j, k;

  fxlgmsv->dates = NULL;
  fxlgmsv->times = NULL;

  fxlgmsv->fx_spot = fx_spot;

  fxlgmsv->fx_sigma = NULL;

  fxlgmsv->correlation = NULL;

  fxlgmsv->dom_irmsv = NULL;
  fxlgmsv->for_irmsv = NULL;

  fxlgmsv->dom_irmsv = (irm_sv *)calloc(1, sizeof(irm_sv));
  fxlgmsv->for_irmsv = (irm_sv *)calloc(1, sizeof(irm_sv));
  if ((!fxlgmsv->dom_irmsv) || (!fxlgmsv->for_irmsv)) {
    err = "allocation failed in init_FXLGMSVUnd";
    goto FREE_RETURN;
  }

  err = init_LGMSV_str(today, num_dates, dates, times, dom_one2F, dom_sigma,
                       dom_tau, dom_alpha, dom_gamma, dom_rho, dom_alphaSV,
                       dom_lamSV, dom_rhoSV, dom_rho2SV, dom_tstar,
                       fxlgmsv->dom_irmsv);
  if (err) {
    goto FREE_RETURN;
  }

  err = init_LGMSV_str(today, num_dates, dates, times, for_one2F, for_sigma,
                       for_tau, for_alpha, for_gamma, for_rho, for_alphaSV,
                       for_lamSV, for_rhoSV, for_rho2SV, for_tstar,
                       fxlgmsv->for_irmsv);
  if (err) {
    goto FREE_RETURN;
  }

  fxlgmsv->num_dates = num_dates;
  fxlgmsv->dates = (long *)calloc(num_dates, sizeof(long));
  fxlgmsv->times = (double *)calloc(num_dates, sizeof(double));

  fxlgmsv->fx_sigma = (double *)calloc(num_dates, sizeof(double));
  fxlgmsv->correlation = f3tensor(0, num_dates - 1, 0, 6, 0, 6);

  if ((!fxlgmsv->dates) || (!fxlgmsv->times) || (!fxlgmsv->dates) ||
      (!fxlgmsv->fx_sigma) || (!fxlgmsv->correlation)) {
    err = "memory allocation failed in init_FXLGMSVUnd";
    goto FREE_RETURN;
  }

  memcpy(fxlgmsv->dates, dates, num_dates * sizeof(long));
  memcpy(fxlgmsv->times, times, num_dates * sizeof(double));

  memcpy(fxlgmsv->fx_sigma, fx_sigma, num_dates * sizeof(double));

  for (i = 0; i < num_dates; ++i) {
    for (j = 0; j < 7; ++j) {
      for (k = 0; k < 7; ++k) {
        fxlgmsv->correlation[i][j][k] = correlation[i][j][k];
      }
    }
  }

FREE_RETURN:
  if (err) {
    fxlgmsv_free_und(fxlgmsv);
  }

  return err;
}

/*	Merge fx        , rates and corr term structures */
Err FXLGMSV_merge_fx_ir_ir_ts(
    /* Domestic TS */
    double *dom_sigma_time, int dom_n_sigma, double *dom_sigma,

    double *dom_tau_time, int dom_n_tau, double *dom_tau,

    double *dom_SV_time, int dom_n_SV, double *dom_alphaSV, double *dom_lamSV,
    double *dom_rhoSV, double *dom_rho2SV,

    /* Foreign TS */
    double *for_sigma_time, int for_n_sigma, double *for_sigma,

    double *for_tau_time, int for_n_tau, double *for_tau,

    double *for_SV_time, int for_n_SV, double *for_alphaSV, double *for_lamSV,
    double *for_rhoSV, double *for_rho2SV,

    /* FX TS */
    double *fx_sigma_time, int fx_n_sigma, double *fx_sigma,

    /* Correlations TS */
    double *rho_time, int rho_n, double ***correlation,

    /* ------- */
    /* Outputs */
    /* ------- */

    /* Merge TS Times*/
    double **merge_time, int *merge_n,

    /* Domestic TS */
    double **merge_dom_sigma, double **merge_dom_tau,

    double **merge_dom_alphaSV, double **merge_dom_lamSV,
    double **merge_dom_rhoSV, double **merge_dom_rho2SV,

    /* Foreign TS */
    double **merge_for_sigma, double **merge_for_tau,

    double **merge_for_alphaSV, double **merge_for_lamSV,
    double **merge_for_rhoSV, double **merge_for_rho2SV,

    /* FX TS */
    double **merge_fx_sigma,

    /* Correlations TS */
    double ****merge_correlation) {
  int i, j, k, index;
  Err err = NULL;

  *merge_time = NULL;

  *merge_dom_sigma = NULL;
  *merge_dom_tau = NULL;

  *merge_dom_alphaSV = NULL;
  *merge_dom_lamSV = NULL;
  *merge_dom_rhoSV = NULL;
  *merge_dom_rho2SV = NULL;

  *merge_for_sigma = NULL;
  *merge_for_tau = NULL;

  *merge_for_alphaSV = NULL;
  *merge_for_lamSV = NULL;
  *merge_for_rhoSV = NULL;
  *merge_for_rho2SV = NULL;

  *merge_fx_sigma = NULL;

  *merge_correlation = NULL;

  *merge_time = (double *)calloc(dom_n_sigma, sizeof(double));
  if (!(*merge_time)) {
    err = "Memory allocation error in FXLGMSV_merge_fx_ir_ir_ts";
    goto FREE_RETURN;
  }

  memcpy(*merge_time, dom_sigma_time, dom_n_sigma * sizeof(double));
  *merge_n = dom_n_sigma;
  num_f_concat_vector(merge_n, merge_time, dom_n_tau, dom_tau_time);
  num_f_concat_vector(merge_n, merge_time, dom_n_SV, dom_SV_time);

  num_f_concat_vector(merge_n, merge_time, for_n_sigma, for_sigma_time);
  num_f_concat_vector(merge_n, merge_time, for_n_tau, for_tau_time);
  num_f_concat_vector(merge_n, merge_time, for_n_SV, for_SV_time);

  num_f_concat_vector(merge_n, merge_time, fx_n_sigma, fx_sigma_time);
  num_f_concat_vector(merge_n, merge_time, rho_n, rho_time);

  num_f_sort_vector(*merge_n, *merge_time);
  num_f_unique_vector(merge_n, *merge_time);

  *merge_dom_sigma = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_tau = (double *)calloc(*merge_n, sizeof(double));

  *merge_dom_alphaSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_lamSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_rhoSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_rho2SV = (double *)calloc(*merge_n, sizeof(double));

  *merge_for_sigma = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_tau = (double *)calloc(*merge_n, sizeof(double));

  *merge_for_alphaSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_lamSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_rhoSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_rho2SV = (double *)calloc(*merge_n, sizeof(double));

  *merge_fx_sigma = (double *)calloc(*merge_n, sizeof(double));

  *merge_correlation = f3tensor(0, *merge_n - 1, 0, 6, 0, 6);

  if (!(*merge_dom_sigma) || !(*merge_dom_tau) || !(*merge_dom_alphaSV) ||
      !(*merge_dom_lamSV) || !(*merge_dom_rhoSV) || !(*merge_dom_rho2SV) ||
      !(*merge_for_sigma) || !(*merge_for_tau) || !(*merge_for_alphaSV) ||
      !(*merge_for_lamSV) || !(*merge_for_rhoSV) || !(*merge_for_rho2SV) ||
      !(*merge_fx_sigma) || !(*merge_correlation)) {
    err = "Memory allocation error in FXLGMSV_merge_fx_ir_ir_ts";
    goto FREE_RETURN;
  }

  for (i = *merge_n - 1; i >= 0; i--) {
    /* Domestic TS */
    (*merge_dom_sigma)[i] =
        dom_sigma[Get_Index((*merge_time)[i], dom_sigma_time, dom_n_sigma)];
    (*merge_dom_tau)[i] =
        dom_tau[Get_Index((*merge_time)[i], dom_tau_time, dom_n_tau)];

    index = Get_Index((*merge_time)[i], dom_SV_time, dom_n_SV);
    (*merge_dom_alphaSV)[i] = dom_alphaSV[index];
    (*merge_dom_lamSV)[i] = dom_lamSV[index];
    (*merge_dom_rhoSV)[i] = dom_rhoSV[index];
    if (!dom_rho2SV) {
      (*merge_dom_rho2SV)[i] = 0.0;
    } else {
      (*merge_dom_rho2SV)[i] = dom_rho2SV[index];
    }

    /* Foreign TS */
    (*merge_for_sigma)[i] =
        for_sigma[Get_Index((*merge_time)[i], for_sigma_time, for_n_sigma)];
    (*merge_for_tau)[i] =
        for_tau[Get_Index((*merge_time)[i], for_tau_time, for_n_tau)];

    index = Get_Index((*merge_time)[i], for_SV_time, for_n_SV);
    (*merge_for_alphaSV)[i] = for_alphaSV[index];
    (*merge_for_lamSV)[i] = for_lamSV[index];
    (*merge_for_rhoSV)[i] = for_rhoSV[index];
    if (!for_rho2SV) {
      (*merge_for_rho2SV)[i] = 0.0;
    } else {
      (*merge_for_rho2SV)[i] = for_rho2SV[index];
    }

    /* FX TS */
    (*merge_fx_sigma)[i] =
        fx_sigma[Get_Index((*merge_time)[i], fx_sigma_time, fx_n_sigma)];

    /* Correlations TS */
    index = Get_Index((*merge_time)[i], rho_time, rho_n);
    for (j = 0; j < 7; j++) {
      for (k = 0; k < 7; k++) {
        (*merge_correlation)[i][j][k] = correlation[index][j][k];
      }
    }
  }

FREE_RETURN:

  if (err) {
    if (*merge_dom_sigma) {
      free(*merge_dom_sigma);
      *merge_dom_sigma = NULL;
    }

    if (*merge_dom_tau) {
      free(*merge_dom_tau);
      *merge_dom_tau = NULL;
    }

    if (*merge_dom_alphaSV) {
      free(*merge_dom_alphaSV);
      *merge_dom_alphaSV = NULL;
    }

    if (*merge_dom_lamSV) {
      free(*merge_dom_lamSV);
      *merge_dom_lamSV = NULL;
    }

    if (*merge_dom_rhoSV) {
      free(*merge_dom_rhoSV);
      *merge_dom_rhoSV = NULL;
    }

    if (*merge_dom_rho2SV) {
      free(*merge_dom_rho2SV);
      *merge_dom_rho2SV = NULL;
    }

    if (*merge_for_sigma) {
      free(*merge_for_sigma);
      *merge_for_sigma = NULL;
    }

    if (*merge_for_tau) {
      free(*merge_for_tau);
      *merge_for_tau = NULL;
    }

    if (*merge_for_alphaSV) {
      free(*merge_for_alphaSV);
      *merge_for_alphaSV = NULL;
    }

    if (*merge_for_lamSV) {
      free(*merge_for_lamSV);
      *merge_for_lamSV = NULL;
    }

    if (*merge_for_rhoSV) {
      free(*merge_for_rhoSV);
      *merge_for_rhoSV = NULL;
    }

    if (*merge_for_rho2SV) {
      free(*merge_for_rho2SV);
      *merge_for_rho2SV = NULL;
    }

    if (*merge_fx_sigma) {
      free(*merge_fx_sigma);
      *merge_fx_sigma = NULL;
    }

    if (*merge_correlation) {
      free_f3tensor(*merge_correlation, 0, *merge_n - 1, 0, 6, 0, 6);
      *merge_correlation = NULL;
    }
  }

  return err;
}

Err fxlgmsv_free_und_struct(SrtUndPtr pUndDesc) {
  SrtFxDesc *pSrtFxPtr;

  pSrtFxPtr = (SrtFxDesc *)(pUndDesc->spec_desc);
  fxlgmsv_free_und(pSrtFxPtr->spec);
  free(pSrtFxPtr->spec);
  free(pSrtFxPtr);
  free(pUndDesc);
  pUndDesc = NULL;
  return NULL;
}

Err SrtInitFXLGMSVUnd(char *undName, /* und name */

                      char *dom_undName, /* domestic underlying name */
                      char *for_undName, /* foreign underlying name */

                      /*	FX Underlying	*/
                      double fx_spot, int fx_n_sigma, long *fx_sigma_date,
                      double *fx_sigma,

                      /*	Correlations	*/
                      long *rho_date, int rho_n, double ***correlation) {
  Err err = NULL;

  int bCleanUpUndFlag = 1;
  SrtUndPtr pUndDesc;
  SrtFxDesc *pSrtFxPtr;
  SrtUndListPtr und_list;
  SrtCurvePtr pdomYieldCurve;
  SrtCurvePtr pforYieldCurve;
  char *domccy;
  char *forccy;
  char domName[256];
  char forName[256];
  fxlgmsv_str *fxlgmsv = NULL;

  int dom_one2F, for_one2F;
  int dom_n_sigma, for_n_sigma;
  double *dom_sigma_time = NULL, *for_sigma_time = NULL;
  double *dom_sigma = NULL, *for_sigma = NULL;
  int dom_n_tau, for_n_tau;
  double *dom_tau_time = NULL, *for_tau_time = NULL;
  double *dom_tau = NULL, *for_tau = NULL;
  double dom_alpha, dom_gamma, dom_rho;
  double for_alpha, for_gamma, for_rho;
  int dom_SV_n, for_SV_n;
  double *dom_SV_time = NULL, *for_SV_time = NULL;
  double *dom_alphaSV = NULL, *for_alphaSV = NULL;
  double *dom_lamSV = NULL, *for_lamSV = NULL;
  double *dom_rhoSV = NULL, *for_rhoSV = NULL;
  double *dom_rho2SV = NULL, *for_rho2SV = NULL;
  double dom_tstar, for_tstar;

  double *merge_time = NULL;
  int merge_n;

  /* Domestic TS */
  double *merge_dom_sigma = NULL;
  double *merge_dom_tau = NULL;

  double *merge_dom_alphaSV = NULL;
  double *merge_dom_lamSV = NULL;
  double *merge_dom_rhoSV = NULL;
  double *merge_dom_rho2SV = NULL;

  /* Foreign TS */
  double *merge_for_sigma = NULL;
  double *merge_for_tau = NULL;

  double *merge_for_alphaSV = NULL;
  double *merge_for_lamSV = NULL;
  double *merge_for_rhoSV = NULL;
  double *merge_for_rho2SV = NULL;

  /* FX TS */
  double *merge_fx_sigma = NULL;

  /* Correlations TS */
  double ***merge_correlation = NULL;

  SrtUndPtr domund, forund;
  long ltoday;
  int i;
  double *fx_sigma_time = NULL;
  double *rho_time = NULL;
  long *merge_date = NULL;

  irm_sv *dom_irmsv;
  irm_sv *for_irmsv;

  domund = lookup_und(dom_undName);
  if (!domund) {
    return serror("Couldn't find domestic underlying %s", domund);
  }
  forund = lookup_und(for_undName);
  if (!forund) {
    return serror("Couldn't find foreign underlying %s", forund);
  }

  ltoday = get_today_from_underlying(domund);

  fxlgmsv = calloc(1, sizeof(fxlgmsv_str));

  // Get yield curves from the yield curve names
  pdomYieldCurve = lookup_curve(((SrtIrDesc *)(domund->spec_desc))->yc_name);
  domccy = get_curve_ccy(pdomYieldCurve);
  pforYieldCurve = lookup_curve(((SrtIrDesc *)(forund->spec_desc))->yc_name);
  forccy = get_curve_ccy(pforYieldCurve);

  // Create the new FXLGMSV underlying
  und_list = get_underlying_list();
  pUndDesc = (SrtUndPtr)calloc(1, sizeof(SrtUndDesc));
  strcpy(pUndDesc->underl_name, undName);
  strupper(pUndDesc->underl_name);
  strip_white_space(pUndDesc->underl_name);
  strcpy(pUndDesc->underl_lbl, "FOREX_UND");
  pUndDesc->underl_ccy = domccy;
  pUndDesc->underl_type = FOREX_UND;

  err = irm_sv_get_struct_from_und(dom_undName, &dom_irmsv);
  if (dom_irmsv->one2F == 1) {
    err = "Domestic underlying should be of type LGMSV2F";
    return err;
  }

  err = irm_sv_get_struct_from_und(for_undName, &for_irmsv);
  if (for_irmsv->one2F == 1) {
    err = "Foreign underlying should be of type LGMSV2F";
    return err;
  }

  err = irm_sv_get_term_struct(dom_undName, &dom_one2F, &dom_n_sigma,
                               &dom_sigma_time, &dom_sigma, &dom_n_tau,
                               &dom_tau_time, &dom_tau, &dom_alpha, &dom_gamma,
                               &dom_rho, &dom_SV_n, &dom_SV_time, &dom_alphaSV,
                               &dom_lamSV, &dom_rhoSV, &dom_rho2SV, &dom_tstar);
  if (err) {
    return err;
  }

  err = irm_sv_get_term_struct(for_undName, &for_one2F, &for_n_sigma,
                               &for_sigma_time, &for_sigma, &for_n_tau,
                               &for_tau_time, &for_tau, &for_alpha, &for_gamma,
                               &for_rho, &for_SV_n, &for_SV_time, &for_alphaSV,
                               &for_lamSV, &for_rhoSV, &for_rho2SV, &for_tstar);
  if (err) {
    return err;
  }

  fx_sigma_time = (double *)calloc(fx_n_sigma, sizeof(double));
  for (i = 0; i < fx_n_sigma; ++i) {
    fx_sigma_time[i] = (fx_sigma_date[i] - ltoday) / 365.0;
  }

  rho_time = (double *)calloc(rho_n, sizeof(double));
  for (i = 0; i < rho_n; ++i) {
    rho_time[i] = (rho_date[i] - ltoday) / 365.0;
  }

  err = FXLGMSV_merge_fx_ir_ir_ts(
      /* Domestic TS */
      dom_sigma_time, dom_n_sigma, dom_sigma,

      dom_tau_time, dom_n_tau, dom_tau,

      dom_SV_time, dom_SV_n, dom_alphaSV, dom_lamSV, dom_rhoSV, dom_rho2SV,

      /* Foreign TS */
      for_sigma_time, for_n_sigma, for_sigma,

      for_tau_time, for_n_tau, for_tau,

      for_SV_time, for_SV_n, for_alphaSV, for_lamSV, for_rhoSV, for_rho2SV,

      /* FX TS */
      fx_sigma_time, fx_n_sigma, fx_sigma,

      /* Correlations TS */
      rho_time, rho_n, correlation,

      /* ------- */
      /* Outputs */
      /* ------- */

      /* Merge TS Times*/
      &merge_time, &merge_n,

      /* Domestic TS */
      &merge_dom_sigma, &merge_dom_tau,

      &merge_dom_alphaSV, &merge_dom_lamSV, &merge_dom_rhoSV, &merge_dom_rho2SV,

      /* Foreign TS */
      &merge_for_sigma, &merge_for_tau,

      &merge_for_alphaSV, &merge_for_lamSV, &merge_for_rhoSV, &merge_for_rho2SV,

      /* FX TS */
      &merge_fx_sigma,

      /* Correlations TS */
      &merge_correlation);
  if (err) {
    return err;
  }

  merge_date = (long *)calloc(merge_n, sizeof(long));
  for (i = 0; i < merge_n; ++i) {
    merge_date[i] = (long)(ltoday + 365.0 * merge_time[i]);
  }

  err = init_FXLGMSVUnd(
      ltoday,
      /* TS Dates and Times */
      merge_n, merge_date, merge_time,

      /* Domestic Underlying */
      dom_one2F,

      merge_dom_sigma, merge_dom_tau, dom_alpha, dom_gamma, dom_rho,

      merge_dom_alphaSV, merge_dom_lamSV, merge_dom_rhoSV, merge_dom_rho2SV,

      dom_tstar,

      /* Foreign Underlying */
      for_one2F,

      merge_for_sigma, merge_for_tau, for_alpha, for_gamma, for_rho,

      merge_for_alphaSV, merge_for_lamSV, merge_for_rhoSV, merge_for_rho2SV,

      for_tstar,

      /* FX Underlying */
      fx_spot, merge_fx_sigma,

      /* Correlations */
      merge_correlation,

      fxlgmsv);

  pSrtFxPtr = calloc(1, sizeof(SrtFxDesc));

  strcpy(pSrtFxPtr->mdl_lbl, "FXLGMSV_UND");
  pSrtFxPtr->mdl_type = FX_LGMSV;
  pSrtFxPtr->mdl_dim = MULTI_FAC;
  pSrtFxPtr->spot = fx_spot;

  strcpy(domName, dom_undName);
  rem_tick_string(domName, domName);
  strcpy(pSrtFxPtr->dom_name, domName);

  strcpy(forName, for_undName);
  rem_tick_string(forName, forName);
  strcpy(pSrtFxPtr->for_name, forName);
  pSrtFxPtr->spec = fxlgmsv;

  pUndDesc->spec_desc = pSrtFxPtr;

  // Put the underlying into the depot
  err = srt_f_lstins(und_list, pUndDesc->underl_name, 0.0, OBJ_PTR_UND,
                     (void *)pUndDesc, &fxlgmsv_free_und_struct,
                     &(pUndDesc->underl_ticker));
  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (merge_date) {
    free(merge_date);
    merge_date = NULL;
  }

  if (fx_sigma_time) {
    free(fx_sigma_time);
    fx_sigma_time = NULL;
  }

  if (rho_time) {
    free(rho_time);
    rho_time = NULL;
  }

  if (merge_time) {
    free(merge_time);
    merge_time = NULL;
  }

  if (merge_dom_sigma) {
    free(merge_dom_sigma);
    merge_dom_sigma = NULL;
  }

  if (merge_dom_tau) {
    free(merge_dom_tau);
    merge_dom_tau = NULL;
  }

  if (merge_dom_alphaSV) {
    free(merge_dom_alphaSV);
    merge_dom_alphaSV = NULL;
  }

  if (merge_dom_lamSV) {
    free(merge_dom_lamSV);
    merge_dom_lamSV = NULL;
  }

  if (merge_dom_rhoSV) {
    free(merge_dom_rhoSV);
    merge_dom_rhoSV = NULL;
  }

  if (merge_dom_rho2SV) {
    free(merge_dom_rho2SV);
    merge_dom_rho2SV = NULL;
  }

  if (merge_for_sigma) {
    free(merge_for_sigma);
    merge_for_sigma = NULL;
  }

  if (merge_for_tau) {
    free(merge_for_tau);
    merge_for_tau = NULL;
  }

  if (merge_for_alphaSV) {
    free(merge_for_alphaSV);
    merge_for_alphaSV = NULL;
  }

  if (merge_for_lamSV) {
    free(merge_for_lamSV);
    merge_for_lamSV = NULL;
  }

  if (merge_for_rhoSV) {
    free(merge_for_rhoSV);
    merge_for_rhoSV = NULL;
  }

  if (merge_for_rho2SV) {
    free(merge_for_rho2SV);
    merge_for_rho2SV = NULL;
  }

  if (merge_fx_sigma) {
    free(merge_fx_sigma);
    merge_fx_sigma = NULL;
  }

  if (merge_correlation) {
    free_f3tensor(merge_correlation, 0, merge_n - 1, 0, 6, 0, 6);
    merge_correlation = NULL;
  }

  if (dom_sigma_time) {
    free(dom_sigma_time);
    dom_sigma_time = NULL;
  }

  if (dom_sigma) {
    free(dom_sigma);
    dom_sigma = NULL;
  }

  if (dom_tau_time) {
    free(dom_tau_time);
    dom_tau_time = NULL;
  }

  if (dom_tau) {
    free(dom_tau);
    dom_tau = NULL;
  }

  if (dom_SV_time) {
    free(dom_SV_time);
    dom_SV_time = NULL;
  }

  if (dom_alphaSV) {
    free(dom_alphaSV);
    dom_alphaSV = NULL;
  }

  if (dom_lamSV) {
    free(dom_lamSV);
    dom_lamSV = NULL;
  }

  if (dom_rhoSV) {
    free(dom_rhoSV);
    dom_rhoSV = NULL;
  }

  if (dom_rho2SV) {
    free(dom_rho2SV);
    dom_rho2SV = NULL;
  }

  if (for_sigma_time) {
    free(for_sigma_time);
    for_sigma_time = NULL;
  }

  if (for_sigma) {
    free(for_sigma);
    for_sigma = NULL;
  }

  if (for_tau_time) {
    free(for_tau_time);
    for_tau_time = NULL;
  }

  if (for_tau) {
    free(for_tau);
    for_tau = NULL;
  }

  if (for_SV_time) {
    free(for_SV_time);
    for_SV_time = NULL;
  }

  if (for_alphaSV) {
    free(for_alphaSV);
    for_alphaSV = NULL;
  }

  if (for_lamSV) {
    free(for_lamSV);
    for_lamSV = NULL;
  }

  if (for_rhoSV) {
    free(for_rhoSV);
    for_rhoSV = NULL;
  }

  if (for_rho2SV) {
    free(for_rho2SV);
    for_rho2SV = NULL;
  }

  return err;
}

Err fxlgmsv_free_model(fxlgmsv_model *fxlgmsvmodel) {
  if (fxlgmsvmodel) {
    if (fxlgmsvmodel->times)
      free(fxlgmsvmodel->times);

    if (fxlgmsvmodel->fx_sigma)
      free(fxlgmsvmodel->fx_sigma);

    if (fxlgmsvmodel->correlation)
      free_f3tensor(fxlgmsvmodel->correlation, 0, fxlgmsvmodel->num_times - 1,
                    0, 6, 0, 6);

    if (fxlgmsvmodel->dom_lgmsv_model) {
      free_LGMSV_model(fxlgmsvmodel->dom_lgmsv_model);
      fxlgmsvmodel->dom_lgmsv_model = NULL;
    }

    if (fxlgmsvmodel->for_lgmsv_model) {
      free_LGMSV_model(fxlgmsvmodel->for_lgmsv_model);
      fxlgmsvmodel->for_lgmsv_model = NULL;
    }

    fxlgmsvmodel = NULL;
  }

  return NULL;
}

Err fxlgmsv_get_struct_from_und(char *und, fxlgmsv_str **fxlgmsv) {
  Err err = NULL;
  SrtUndPtr pUndDesc;
  SrtFxDesc *pFxDesc;

  // Get the underlying through its name and check it exists
  // Check on the underlying type
  pUndDesc = lookup_und(und);
  if (!pUndDesc) {
    err = "Undefined underlying";
    return err;
  }
  if (!(ISUNDTYPE(pUndDesc, FOREX_UND))) {
    err = "Not a Forex Underlying";
  }
  if (get_mdltype_from_fxund(pUndDesc) != FX_LGMSV) {
    err = serror("Underlying %s is not of type FX_LGMSV", und);
    return err;
  }

  // Extract the information from the underlying
  pFxDesc = (SrtFxDesc *)(pUndDesc->spec_desc);
  *fxlgmsv = (fxlgmsv_str *)(pFxDesc->spec);

  return err;
}

Err qtolgmsv_check_und(char *und, int *for_one2F) {
  Err err = NULL;
  fxlgmsv_str *qtolgmsv;

  err = fxlgmsv_get_struct_from_und(und, &qtolgmsv);
  if (err) {
    return err;
  }

  if ((qtolgmsv->dom_irmsv->alpha > 1.0e-5) ||
      (qtolgmsv->dom_irmsv->alphaSV[0] > 1.0e-5)) {
    err = "domestic underlying should be LGM1F";
    return err;
  }

  *for_one2F = 1;
  if (qtolgmsv->for_irmsv->alpha > 1.0e-5) {
    *for_one2F = 2;
  }

  return err;
}

Err fxlgmsv_get_fx_and_correl_ts(char *und, double *fx_spot, int *num_fx,
                                 double **fx_time, double **fx_vol,
                                 int *num_rho, double **rho_time,
                                 double ****correlation) {
  Err err = NULL;
  fxlgmsv_str *fxlgmsv;
  int i, j, k;

  err = fxlgmsv_get_struct_from_und(und, &fxlgmsv);
  if (err) {
    goto FREE_RETURN;
  }

  *num_fx = fxlgmsv->num_dates;
  *fx_spot = fxlgmsv->fx_spot;
  *fx_time = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *fx_vol = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *num_rho = fxlgmsv->num_dates;
  *rho_time = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *correlation = f3tensor(0, fxlgmsv->num_dates - 1, 0, 6, 0, 6);
  if ((!(*fx_vol)) || (!(*fx_time)) || (!(*rho_time)) || (!(*correlation))) {
    err = "Allocation failed in fxlgmsv_get_fx_ts";
    goto FREE_RETURN;
  }

  for (i = 0; i < fxlgmsv->num_dates; ++i) {
    (*fx_time)[i] = fxlgmsv->times[i];
    (*fx_vol)[i] = fxlgmsv->fx_sigma[i];
    (*rho_time)[i] = fxlgmsv->times[i];
    for (j = 0; j < 7; ++j) {
      for (k = 0; k < 7; ++k) {
        (*correlation)[i][j][k] = fxlgmsv->correlation[i][j][k];
      }
    }
  }

FREE_RETURN:
  if (err) {
    if (*fx_vol) {
      free(*fx_vol);
      *fx_vol = NULL;
    }
    if (*fx_time) {
      free(*fx_time);
      *fx_time = NULL;
    }
    if (*rho_time) {
      free(*rho_time);
      *rho_time = NULL;
    }
    if (*correlation) {
      free_f3tensor(*correlation, 0, fxlgmsv->num_dates - 1, 0, 6, 0, 6);
      *correlation = NULL;
    }
  }

  return err;
}

void init_NULL_FXLGMSV_struct(fxlgmsv_str *fxlgmsv) {
  if (fxlgmsv) {
    fxlgmsv->dates = NULL;
    fxlgmsv->times = NULL;
    fxlgmsv->fx_sigma = NULL;
    fxlgmsv->correlation = NULL;
  }
}

Err fxlgmsv_get_ts_from_und(char *undName, /* fxlgmsv und name */

                            /* TS Dates and Times */
                            int *num_dates, long **dates, double **times,

                            /* Domestic Underlying */
                            int *dom_one2F,

                            double **dom_sigma,

                            double **dom_tau,

                            double *dom_alpha, double *dom_gamma,
                            double *dom_rho,

                            double **dom_alphaSV, double **dom_lamSV,
                            double **dom_rhoSV, double **dom_rho2SV,

                            double *dom_tstar,

                            /* Foreign Underlying */
                            int *for_one2F,

                            double **for_sigma,

                            double **for_tau,

                            double *for_alpha, double *for_gamma,
                            double *for_rho,

                            double **for_alphaSV, double **for_lamSV,
                            double **for_rhoSV, double **for_rho2SV,

                            double *for_tstar,

                            /* FX Underlying */
                            double *fx_spot,

                            double **fx_sigma,

                            /* Correlations */
                            double ****correlation) {
  Err err = NULL;
  fxlgmsv_str *fxlgmsv;
  int i, j, k;

  err = fxlgmsv_get_struct_from_und(undName, &fxlgmsv);
  if (err) {
    goto FREE_RETURN;
  }

  *dom_one2F = fxlgmsv->dom_irmsv->one2F;
  *dom_alpha = fxlgmsv->dom_irmsv->alpha;
  *dom_gamma = fxlgmsv->dom_irmsv->gamma;
  *dom_rho = fxlgmsv->dom_irmsv->rho;
  *dom_tstar = fxlgmsv->dom_irmsv->tstar;

  *for_one2F = fxlgmsv->for_irmsv->one2F;
  *for_alpha = fxlgmsv->for_irmsv->alpha;
  *for_gamma = fxlgmsv->for_irmsv->gamma;
  *for_rho = fxlgmsv->for_irmsv->rho;
  *for_tstar = fxlgmsv->for_irmsv->tstar;

  *fx_spot = fxlgmsv->fx_spot;

  *num_dates = fxlgmsv->num_dates;
  *dates = (long *)calloc(fxlgmsv->num_dates, sizeof(long));
  *times = (double *)calloc(fxlgmsv->num_dates, sizeof(double));

  *dom_sigma = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *dom_tau = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *dom_alphaSV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *dom_lamSV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *dom_rhoSV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *dom_rho2SV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));

  *for_sigma = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *for_tau = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *for_alphaSV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *for_lamSV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *for_rhoSV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));
  *for_rho2SV = (double *)calloc(fxlgmsv->num_dates, sizeof(double));

  *fx_sigma = (double *)calloc(fxlgmsv->num_dates, sizeof(double));

  *correlation = f3tensor(0, fxlgmsv->num_dates - 1, 0, 6, 0, 6);

  if ((!(*dates)) || (!(*times)) || (!(*dates)) || (!(*dom_sigma)) ||
      (!(*dom_tau)) || (!(*dom_alphaSV)) || (!(*dom_lamSV)) ||
      (!(*dom_rhoSV)) || (!(*dom_rho2SV)) || (!(*for_sigma)) || (!(*for_tau)) ||
      (!(*for_alphaSV)) || (!(*for_lamSV)) || (!(*for_rhoSV)) ||
      (!(*for_rho2SV)) || (!(*fx_sigma)) || (!(*correlation))) {
    err = "memory allocation failed in FXLGMSV_get_ts_from_und";
    goto FREE_RETURN;
  }

  memcpy(*dates, fxlgmsv->dates, fxlgmsv->num_dates * sizeof(long));
  memcpy(*times, fxlgmsv->times, fxlgmsv->num_dates * sizeof(double));

  memcpy(*dom_sigma, fxlgmsv->dom_irmsv->sigma,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*dom_tau, fxlgmsv->dom_irmsv->tau,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*dom_alphaSV, fxlgmsv->dom_irmsv->alphaSV,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*dom_lamSV, fxlgmsv->dom_irmsv->lamSV,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*dom_rhoSV, fxlgmsv->dom_irmsv->rhoSV,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*dom_rho2SV, fxlgmsv->dom_irmsv->rho2SV,
         fxlgmsv->num_dates * sizeof(double));

  memcpy(*for_sigma, fxlgmsv->for_irmsv->sigma,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*for_tau, fxlgmsv->for_irmsv->tau,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*for_alphaSV, fxlgmsv->for_irmsv->alphaSV,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*for_lamSV, fxlgmsv->for_irmsv->lamSV,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*for_rhoSV, fxlgmsv->for_irmsv->rhoSV,
         fxlgmsv->num_dates * sizeof(double));
  memcpy(*for_rho2SV, fxlgmsv->for_irmsv->rho2SV,
         fxlgmsv->num_dates * sizeof(double));

  memcpy(*fx_sigma, fxlgmsv->fx_sigma, fxlgmsv->num_dates * sizeof(double));

  for (i = 0; i < fxlgmsv->num_dates; ++i) {
    for (j = 0; j < 7; ++j) {
      for (k = 0; k < 7; ++k) {
        (*correlation)[i][j][k] = fxlgmsv->correlation[i][j][k];
      }
    }
  }

FREE_RETURN:

  if (err) {
    if (*dates) {
      free(*dates);
      *dates = NULL;
    }
    if (*times) {
      free(*times);
      *times = NULL;
    }

    if (*dom_sigma) {
      free(*dom_sigma);
      *dom_sigma = NULL;
    }
    if (*dom_tau) {
      free(*dom_tau);
      *dom_tau = NULL;
    }
    if (*dom_alphaSV) {
      free(*dom_alphaSV);
      *dom_alphaSV = NULL;
    }
    if (*dom_lamSV) {
      free(*dom_lamSV);
      *dom_lamSV = NULL;
    }
    if (*dom_rhoSV) {
      free(*dom_rhoSV);
      *dom_rhoSV = NULL;
    }
    if (*dom_rho2SV) {
      free(*dom_rho2SV);
      *dom_rho2SV = NULL;
    }

    if (*for_sigma) {
      free(*for_sigma);
      *for_sigma = NULL;
    }
    if (*for_tau) {
      free(*for_tau);
      *for_tau = NULL;
    }
    if (*for_alphaSV) {
      free(*for_alphaSV);
      *for_alphaSV = NULL;
    }
    if (*for_lamSV) {
      free(*for_lamSV);
      *for_lamSV = NULL;
    }
    if (*for_rhoSV) {
      free(*for_rhoSV);
      *for_rhoSV = NULL;
    }
    if (*for_rho2SV) {
      free(*for_rho2SV);
      *for_rho2SV = NULL;
    }

    if (*fx_sigma) {
      free(*fx_sigma);
      *fx_sigma = NULL;
    }
    if (*correlation)

    {
      free_f3tensor(*correlation, 0, fxlgmsv->num_dates - 1, 0, 6, 0, 6);
      *correlation = NULL;
    }
  }

  return err;
}

void init_NULL_FXLGMSV_model(fxlgmsv_model *model) {
  if (model) {
    init_NULL_LGMSV_model(model->dom_lgmsv_model);
    init_NULL_LGMSV_model(model->for_lgmsv_model);

    model->times = NULL;
    model->fx_sigma = NULL;
    model->correlation = NULL;
  }
}

Err fxlgmsv_fill_model_from_struct(fxlgmsv_str *str, fxlgmsv_model *model) {
  Err err = NULL;
  int i, j, k;

  model->fx_spot = str->fx_spot;

  model->times = NULL;
  model->fx_sigma = NULL;
  model->correlation = NULL;
  model->dom_lgmsv_model = NULL;
  model->for_lgmsv_model = NULL;

  model->dom_lgmsv_model = (LGMSV_model *)calloc(1, sizeof(LGMSV_model));
  model->for_lgmsv_model = (LGMSV_model *)calloc(1, sizeof(LGMSV_model));
  if ((!model->dom_lgmsv_model) || (!model->for_lgmsv_model)) {
    err = "allocation failed in fxlgmsv_fill_model_from_struct";
    goto FREE_RETURN;
  }

  err = fill_LGMSV_model_from_irm_sv(str->dom_irmsv, model->dom_lgmsv_model);
  if (err) {
    goto FREE_RETURN;
  }

  err = fill_LGMSV_model_from_irm_sv(str->dom_irmsv, model->dom_lgmsv_model);
  if (err) {
    goto FREE_RETURN;
  }

  model->times = (double *)calloc(str->num_dates, sizeof(double));
  model->fx_sigma = (double *)calloc(str->num_dates, sizeof(double));
  model->correlation = f3tensor(0, str->num_dates - 1, 0, 6, 0, 6);
  if ((!model->times) || (!model->fx_sigma) || (!model->fx_sigma)) {
    err = "Allocation failed in fxlgmsv_fill_model_from_struct";
    goto FREE_RETURN;
  }

  memcpy(model->times, str->times, str->num_dates * sizeof(double));
  memcpy(model->fx_sigma, str->fx_sigma, str->num_dates * sizeof(double));

  model->num_times = str->num_dates;
  for (i = 0; i < str->num_dates; ++i) {
    for (j = 0; j < 7; ++j) {
      for (k = 0; k < 7; ++k) {
        model->correlation[i][j][k] = str->correlation[i][j][k];
      }
    }
  }

FREE_RETURN:

  if (err) {
    fxlgmsv_free_model(model);
  }
  return err;
}

Err fxlgmsv_get_model_from_und(char *undName, fxlgmsv_model **fxlgmsv) {
  Err err = NULL;
  fxlgmsv_str *fxlgmsv_struct = NULL;

  *fxlgmsv = (fxlgmsv_model *)calloc(1, sizeof(fxlgmsv_model));
  if (!(*fxlgmsv)) {
    err = "allocation failed in fxlgmsv_get_model_from_und";
    goto FREE_RETURN;
  }

  err = fxlgmsv_get_struct_from_und(undName, &fxlgmsv_struct);
  if (err) {
    goto FREE_RETURN;
  }

  err = fxlgmsv_fill_model_from_struct(fxlgmsv_struct, *fxlgmsv);
  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:
  if (err) {
    if (*fxlgmsv) {
      fxlgmsv_free_model(*fxlgmsv);
    }
  }

  return err;
}

void HestonConstructLegendreGrid(int iNbInt, double *dBreakPoints,
                                 int *iNbPoints, double *dX, double *dW) {
  int i;
  int iStartIndex;
  double dStart, dEnd;

  dStart = 0;
  iStartIndex = 0;

  for (i = 0; i < iNbInt; i++) {
    dEnd = dBreakPoints[i];
    gauleg(dStart, dEnd, &(dX[iStartIndex - 1]), &(dW[iStartIndex - 1]),
           iNbPoints[i]);
    dStart = dEnd;
    iStartIndex += iNbPoints[i];
  }
}

void HestonSolveODEGride(int iNbX, double *X, double Fwd, int endi, double *dt,
                         double *CoefIntT, double *Coef1ReT, double *Coef1ImT,
                         double *Coef2ReT, double *Coef2ImT, double *Coef3ReT,
                         double *Coef3ImT,

                         double dIntegParam,

                         /* Outputs */
                         double *IntRe, double *IntIm) {
  double FreqRe, FreqIm, ResRe, ResIm;
  double coef1, coef2, temp1, temp2, temp3;
  double dExpReal, dImag;
  int i;

  FreqIm = -(dIntegParam + 1);
  coef1 = dIntegParam * (dIntegParam + 1);
  coef2 = 2.0 * dIntegParam + 1;

  for (i = 0; i < iNbX; i++) {
    FreqRe = X[i];

    LGMSVSolveODE(FreqRe, FreqIm, Fwd, endi, dt, CoefIntT, Coef1ReT, Coef1ImT,
                  Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, &ResRe, &ResIm);

    dExpReal = exp(ResRe);
    dImag = fmod(ResIm, PI2);
    ResRe = dExpReal * cos(dImag);
    ResIm = dExpReal * sin(dImag);

    temp1 = coef1 - FreqRe * FreqRe;
    temp2 = FreqRe * coef2;
    temp3 = temp1 * temp1 + temp2 * temp2;

    IntRe[i] = (ResRe * temp1 + ResIm * temp2) / temp3;
    IntIm[i] = (ResIm * temp1 - ResRe * temp2) / temp3;
  }
}

void HestonComputeIntegralLegendre(int iNbX, double *X, double *W,
                                   double *IntRe, double *IntIm, int iNbStrike,
                                   double *dLogStrike, double dIntegParam,
                                   double *result) {
  int i, j;
  double theta, temp, temp1;

  for (j = 0; j < iNbStrike; j++) {
    temp = 0.0;

    for (i = 0; i < iNbX; i++) {
      theta = X[i] * dLogStrike[j];
      temp1 = (cos(theta) * IntRe[i] + sin(theta) * IntIm[i]);
      temp += W[i] * temp1;
    }

    result[j] = exp(-dIntegParam * dLogStrike[j]) * temp / PI;
  }
}

void HestonClosedFormApprox(
    /* Parameter of diffusion */
    int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigma, double *dAlpha, double *dLevelEps,
    double *dLambdaEps, double *dRho,

    /* Product description */
    double dFwd, double dStrike, double dExTime,

    /* Numerical parameters */
    int iNbX, double dIntegParam,

    /* Output */
    double *Price) {
  Err err = NULL;
  double lambdaArray[10];

  double InitX, t1, t2, temp_vol;
  int i, endi;

  double *dt = NULL, *CoefIntT = NULL, *Coef1ReT = NULL, *Coef2ReT = NULL,
         *Coef2ImT = NULL, *Coef3ReT = NULL, *Coef3ImT = NULL;

  double *X = NULL, *W = NULL, *IntRe = NULL, *IntIm = NULL;

  double LogStrike;

  double dBreakPoints[4];
  int iBreakNbX[4];

  endi = Get_Index(dExTime, dPWTime, iNbPWTime);

  /* Memory allocation */

  dt = dvector(0, endi);
  CoefIntT = dvector(0, endi);
  Coef1ReT = dvector(0, endi);
  Coef2ReT = dvector(0, endi);
  Coef2ImT = dvector(0, endi);
  Coef3ReT = dvector(0, endi);
  Coef3ImT = dvector(0, endi);

  X = dvector(0, iNbX - 1);
  IntRe = dvector(0, iNbX - 1);
  IntIm = dvector(0, iNbX - 1);

  if (!dt || !CoefIntT || !Coef1ReT || !Coef2ReT || !Coef2ImT || !Coef3ReT ||
      !Coef3ImT || !X || !IntRe || !IntIm) {
    err = "Memory Allocation failure (1) in LGMSVClosedFormNew";
    goto FREE_RETURN;
  }

  memset(lambdaArray, 0, 10 * sizeof(double));

  /* Precalculations */

  InitX = log(fabs(dFwd));

  t1 = 0.0;
  for (i = 0; i <= endi; i++) {
    /* Precalculation on the option i*/

    /* First the time discretisation */
    if (i < endi) {
      t2 = dPWTime[i];
    } else {
      t2 = dExTime;
    }

    /* Calculate constant values */
    temp_vol = dSigma[i];
    dt[i] = (t2 - t1);
    CoefIntT[i] = dLevelEps[i]; // * dLambdaEps[i];
    Coef1ReT[i] = -0.5 * dAlpha[i] * dAlpha[i];
    Coef2ReT[i] = dLambdaEps[i];
    Coef2ImT[i] = -dAlpha[i] * dRho[i] * temp_vol;
    temp_vol *= temp_vol;
    Coef3ReT[i] = 0.5 * temp_vol;
    Coef3ImT[i] = 0.5 * temp_vol;

    t1 = t2;
  }

  /* construct the grid */

  /* Gauss Legendre */
  W = dvector(0, iNbX - 1);
  if (!W) {
    err = "Memory Allocation failure (2) in LGMSVClosedFormNew";
    goto FREE_RETURN;
  }

  dBreakPoints[0] = HESTON_FIRSTBREAK;
  dBreakPoints[1] = 1000;
  dBreakPoints[2] = 100000;
  dBreakPoints[3] = HESTON_LASTBREAK;

  iBreakNbX[0] = (int)(iNbX / 4);
  iBreakNbX[1] = iBreakNbX[0];
  iBreakNbX[2] = iBreakNbX[0];
  iBreakNbX[3] = iNbX - 3 * iBreakNbX[0];

  HestonConstructLegendreGrid(4, dBreakPoints, iBreakNbX, X, W);

  /* calculate the integrand */
  HestonSolveODEGride(iNbX, X, InitX, endi, dt, CoefIntT, Coef1ReT, Coef1ReT,
                      Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, dIntegParam,
                      IntRe, IntIm);

  /* do the integral */
  LogStrike = log(fabs(dStrike));

  HestonComputeIntegralLegendre(iNbX, X, W, IntRe, IntIm, 1, &LogStrike,
                                dIntegParam, Price);

FREE_RETURN:

  if (dt)
    free_dvector(dt, 0, endi);
  if (CoefIntT)
    free_dvector(CoefIntT, 0, endi);
  if (Coef1ReT)
    free_dvector(Coef1ReT, 0, endi);
  if (Coef2ReT)
    free_dvector(Coef2ReT, 0, endi);
  if (Coef2ImT)
    free_dvector(Coef2ImT, 0, endi);
  if (Coef3ReT)
    free_dvector(Coef3ReT, 0, endi);
  if (Coef3ImT)
    free_dvector(Coef3ImT, 0, endi);

  if (X)
    free_dvector(X, 0, iNbX - 1);
  if (W)
    free_dvector(W, 0, iNbX - 1);
  if (IntRe)
    free_dvector(IntRe, 0, iNbX - 1);
  if (IntIm)
    free_dvector(IntIm, 0, iNbX - 1);
}

void ComputePWCBondVol(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime, double *dSigma,

    double dLambda, double dTStar,

    double *dAlphaLGM, double *dRhoLGM, double dGammaLGM,

    double dFixTime, double dBondMat,

    double *dPWCSigma) {
  int i, FixTimeIndex;
  double ZCVar1, ZCVar2, ZCCoVar;
  double dLambda2;

  FixTimeIndex = Get_Index(dFixTime, dPWTime, iNbPWTime);

  dLambda2 = dLambda + dGammaLGM;

  ZCVar1 = (dSigma[0] / dLambda) * (dSigma[0] / dLambda) *
           (0.5 *
                (exp(2 * dLambda * dTStar) -
                 exp(2 * dLambda * (dTStar - dPWTime[0]))) /
                dLambda -
            2.0 * exp(dLambda * (dTStar - dBondMat)) *
                (exp(dLambda * dTStar) - exp(dLambda * (dTStar - dPWTime[0]))) /
                dLambda +
            exp(2 * dLambda * (dTStar - dBondMat)) * dPWTime[0]);

  ZCVar2 =
      (dAlphaLGM[0] * dSigma[0] / dLambda2) *
      (dAlphaLGM[0] * dSigma[0] / dLambda2) *
      (0.5 *
           (exp(2 * dLambda2 * dTStar) -
            exp(2 * dLambda2 * (dTStar - dPWTime[0]))) /
           dLambda2 -
       2.0 * exp(dLambda2 * (dTStar - dBondMat)) *
           (exp(dLambda2 * dTStar) - exp(dLambda2 * (dTStar - dPWTime[0]))) /
           dLambda2 +
       exp(2 * dLambda2 * (dTStar - dBondMat)) * dPWTime[0]);

  ZCCoVar =
      dAlphaLGM[0] * dRhoLGM[0] * (dSigma[0] / dLambda) *
      (dSigma[0] / dLambda2) *
      ((exp((dLambda + dLambda2) * dTStar) -
        exp((dLambda + dLambda2) * (dTStar - dPWTime[0]))) /
           (dLambda + dLambda2) -
       exp(dLambda * (dTStar - dBondMat)) *
           (exp(dLambda2 * dTStar) - exp(dLambda2 * (dTStar - dPWTime[0]))) /
           dLambda2 -
       exp(dLambda2 * (dTStar - dBondMat)) *
           (exp(dLambda * dTStar) - exp(dLambda * (dTStar - dPWTime[0]))) /
           dLambda +
       exp((dLambda + dLambda2) * (dTStar - dBondMat)) * dPWTime[0]);

  dPWCSigma[0] = sqrt((ZCVar1 + ZCVar2 + 2 * ZCCoVar) / dPWTime[0]);

  for (i = 1; i < FixTimeIndex; ++i) {
    ZCVar1 = (dSigma[i] / dLambda) * (dSigma[i] / dLambda) *
             (0.5 *
                  (exp(2 * dLambda * (dTStar - dPWTime[i - 1])) -
                   exp(2 * dLambda * (dTStar - dPWTime[i]))) /
                  dLambda -
              2.0 * exp(dLambda * (dTStar - dBondMat)) *
                  (exp(dLambda * (dTStar - dPWTime[i - 1])) -
                   exp(dLambda * (dTStar - dPWTime[i]))) /
                  dLambda +
              exp(2 * dLambda * (dTStar - dBondMat)) *
                  (dPWTime[i] - dPWTime[i - 1]));

    ZCVar2 = (dAlphaLGM[i] * dSigma[i] / dLambda2) *
             (dAlphaLGM[i] * dSigma[i] / dLambda2) *
             (0.5 *
                  (exp(2 * dLambda2 * (dTStar - dPWTime[i - 1])) -
                   exp(2 * dLambda2 * (dTStar - dPWTime[i]))) /
                  dLambda2 -
              2.0 * exp(dLambda2 * (dTStar - dBondMat)) *
                  (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
                   exp(dLambda2 * (dTStar - dPWTime[i]))) /
                  dLambda2 +
              exp(2 * dLambda2 * (dTStar - dBondMat)) *
                  (dPWTime[i] - dPWTime[i - 1]));

    ZCCoVar = dAlphaLGM[i] * dRhoLGM[i] * (dSigma[i] / dLambda) *
              (dSigma[i] / dLambda2) *
              ((exp((dLambda + dLambda2) * (dTStar - dPWTime[i - 1])) -
                exp((dLambda + dLambda2) * (dTStar - dPWTime[i]))) /
                   (dLambda + dLambda2) -
               exp(dLambda * (dTStar - dBondMat)) *
                   (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
                    exp(dLambda2 * (dTStar - dPWTime[i]))) /
                   dLambda2 -
               exp(dLambda2 * (dTStar - dBondMat)) *
                   (exp(dLambda * (dTStar - dPWTime[i - 1])) -
                    exp(dLambda * (dTStar - dPWTime[i]))) /
                   dLambda +
               exp((dLambda + dLambda2) * (dTStar - dBondMat)) *
                   (dPWTime[i] - dPWTime[i - 1]));

    dPWCSigma[i] =
        sqrt((ZCVar1 + ZCVar2 + 2 * ZCCoVar) / (dPWTime[i] - dPWTime[i - 1]));
  }

  ZCVar1 =
      (dSigma[i] / dLambda) * (dSigma[i] / dLambda) *
      (0.5 *
           (exp(2 * dLambda * (dTStar - dPWTime[i - 1])) -
            exp(2 * dLambda * (dTStar - dFixTime))) /
           dLambda -
       2.0 * exp(dLambda * (dTStar - dBondMat)) *
           (exp(dLambda * (dTStar - dPWTime[i - 1])) -
            exp(dLambda * (dTStar - dFixTime))) /
           dLambda +
       exp(2 * dLambda * (dTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  ZCVar2 =
      (dAlphaLGM[i] * dSigma[i] / dLambda2) *
      (dAlphaLGM[i] * dSigma[i] / dLambda2) *
      (0.5 *
           (exp(2 * dLambda2 * (dTStar - dPWTime[i - 1])) -
            exp(2 * dLambda2 * (dTStar - dFixTime))) /
           dLambda2 -
       2.0 * exp(dLambda2 * (dTStar - dBondMat)) *
           (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
            exp(dLambda2 * (dTStar - dFixTime))) /
           dLambda2 +
       exp(2 * dLambda2 * (dTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  ZCCoVar = dAlphaLGM[i] * dRhoLGM[i] * (dSigma[i] / dLambda) *
            (dSigma[i] / dLambda2) *
            ((exp((dLambda + dLambda2) * (dTStar - dPWTime[i - 1])) -
              exp((dLambda + dLambda2) * (dTStar - dFixTime))) /
                 (dLambda + dLambda2) -
             exp(dLambda * (dTStar - dBondMat)) *
                 (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
                  exp(dLambda2 * (dTStar - dFixTime))) /
                 dLambda2 -
             exp(dLambda2 * (dTStar - dBondMat)) *
                 (exp(dLambda * (dTStar - dPWTime[i - 1])) -
                  exp(dLambda * (dTStar - dFixTime))) /
                 dLambda +
             exp((dLambda + dLambda2) * (dTStar - dBondMat)) *
                 (dFixTime - dPWTime[i - 1]));

  dPWCSigma[i] =
      sqrt((ZCVar1 + ZCVar2 + 2 * ZCCoVar) / (dFixTime - dPWTime[i - 1]));
}

void FwdAdjustLambdaOnDom(int iNbPWTime, double *dPWTime, double *dSigma,
                          double *dAlpha, double *dLambdaEps, double *dRho,
                          double *dRho2,

                          double dLambda, double *dAlphaLGM, double *dRhoLGM,
                          double dGammaLGM,

                          double dTinitial, double dTfinal,

                          double *dLambdaEpsNew) {
  int i;
  double dLambda2;
  double Coeff1, Coeff2;

  dLambda2 = dLambda + dGammaLGM;
  Coeff1 = (exp(dLambda * (dTinitial - dTfinal)) - 1) / dLambda;
  Coeff2 = (exp(dLambda2 * (dTinitial - dTfinal)) - 1) / dLambda2;

  for (i = 0; i < iNbPWTime; ++i) {
    dLambdaEpsNew[i] =
        dLambdaEps[i] +
        dAlpha[i] * dSigma[i] * (Coeff1 * dRho[i] + Coeff2 * dRho2[i]);
  }
}

void FwdAdjustLevelAndLambdaOnFor(int iNbPWTime, double *dPWTime,
                                  double *dforAlpha, double *dforLevelEps,
                                  double *dforLambdaEps,

                                  double *ddomSigma, double ddomLambda,
                                  double *ddomAlphaLGM, double *ddomRhoLGM,
                                  double ddomGammaLGM,

                                  double ***dCorrelation,

                                  double dTinitial, double dTfinal,

                                  double *dforLevelEpsNew,
                                  double *dforLambdaEpsNew) {
  int i;
  double ddomLambda2;
  double Coeff1, Coeff2;
  double adjustment;

  ddomLambda2 = ddomLambda + ddomGammaLGM;
  Coeff1 = (exp(ddomLambda * (dTinitial - dTfinal)) - 1) / ddomLambda;
  Coeff2 = (exp(ddomLambda2 * (dTinitial - dTfinal)) - 1) / ddomLambda2;

  for (i = 0; i < iNbPWTime; ++i) {
    adjustment =
        0.5 * dforAlpha[i] * ddomSigma[i] *
        (Coeff1 * dCorrelation[i][5][0] + Coeff2 * dCorrelation[i][5][1]);
    dforLambdaEpsNew[i] = dforLambdaEps[i] + adjustment;
    dforLevelEpsNew[i] = dforLevelEps[i] + adjustment;
  }
}

void QtoAdjustLevelAndLambda(int iNbPWTime, double *dPWTime,

                             double *ddomSigma, double ddomLambda,
                             double *ddomAlphaLGM, double *ddomRhoLGM,
                             double ddomGammaLGM,

                             double *dforSigma, double dforLambda,
                             double *dforAlphaLGM, double *dforRhoLGM,
                             double dforGammaLGM, double *dforLevelEps,
                             double *dforLambdaEps, double *dforAlpha,

                             double *dfxSigma,

                             double ***dCorrelation,

                             double dTStar,

                             double *dforLevelEpsNew,
                             double *dforLambdaEpsNew) {
  int i;
  double ddomLambda2;
  double domZCVol1, domZCVol2;
  double dforLambda2;
  double forZCVol1, forZCVol2;

  ddomLambda2 = ddomLambda + ddomGammaLGM;
  dforLambda2 = dforLambda + dforGammaLGM;

  for (i = 0; i < iNbPWTime; ++i) {
    domZCVol1 = -ddomSigma[i] *
                (exp(ddomLambda * (dTStar - dPWTime[i])) - 1.0) / ddomLambda;
    domZCVol2 = -ddomAlphaLGM[i] * ddomSigma[i] *
                (exp(ddomLambda2 * (dTStar - dPWTime[i])) - 1.0) / ddomLambda2;

    dforLevelEpsNew[i] =
        dforLevelEps[i] -
        0.5 * dforAlpha[i] * dfxSigma[i] * dCorrelation[i][5][6] +
        0.5 * dforAlpha[i] * domZCVol1 * dCorrelation[i][5][0] +
        0.5 * dforAlpha[i] * domZCVol2 * dCorrelation[i][5][1];

    forZCVol1 = -dforSigma[i] *
                (exp(dforLambda * (dTStar - dPWTime[i])) - 1.0) / dforLambda;
    forZCVol2 = -dforAlphaLGM[i] * dforSigma[i] *
                (exp(dforLambda2 * (dTStar - dPWTime[i])) - 1.0) / dforLambda2;

    dforLambdaEpsNew[i] =
        dforLambdaEps[i] + dforAlpha[i] * forZCVol1 * dCorrelation[i][5][3] +
        dforAlpha[i] * forZCVol2 * dCorrelation[i][5][4] +
        0.5 * dforAlpha[i] * dfxSigma[i] * dCorrelation[i][5][6] -
        0.5 * dforAlpha[i] * domZCVol1 * dCorrelation[i][5][0] -
        0.5 * dforAlpha[i] * domZCVol2 * dCorrelation[i][5][1];
  }
}

void ComputePWCBondVolAndCorrel(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime, double *dSigma, double *dAlpha, double *dRho,
    double *dRho2,

    int sign,

    double dLambda, double dTStar,

    double *dAlphaLGM, double *dRhoLGM, double dGammaLGM,

    double dFixTime, double dBondMat,

    double *dPWCSigma, double *dPWCRho) {
  int i, FixTimeIndex;
  double ZCVar1, ZCVar2, ZCCoVar;
  double ZCVol, ZCVol1, ZCVol2;
  double dLambda2;

  FixTimeIndex = Get_Index(dFixTime, dPWTime, iNbPWTime);

  dLambda2 = dLambda + dGammaLGM;

  ZCVar1 = (dSigma[0] / dLambda) * (dSigma[0] / dLambda) *
           (0.5 *
                (exp(2 * dLambda * dTStar) -
                 exp(2 * dLambda * (dTStar - dPWTime[0]))) /
                dLambda -
            2.0 * exp(dLambda * (dTStar - dBondMat)) *
                (exp(dLambda * dTStar) - exp(dLambda * (dTStar - dPWTime[0]))) /
                dLambda +
            exp(2 * dLambda * (dTStar - dBondMat)) * dPWTime[0]);

  ZCVar2 =
      (dAlphaLGM[0] * dSigma[0] / dLambda2) *
      (dAlphaLGM[0] * dSigma[0] / dLambda2) *
      (0.5 *
           (exp(2 * dLambda2 * dTStar) -
            exp(2 * dLambda2 * (dTStar - dPWTime[0]))) /
           dLambda2 -
       2.0 * exp(dLambda2 * (dTStar - dBondMat)) *
           (exp(dLambda2 * dTStar) - exp(dLambda2 * (dTStar - dPWTime[0]))) /
           dLambda2 +
       exp(2 * dLambda2 * (dTStar - dBondMat)) * dPWTime[0]);

  ZCCoVar =
      dAlphaLGM[0] * dRhoLGM[0] * (dSigma[0] / dLambda) *
      (dSigma[0] / dLambda2) *
      ((exp((dLambda + dLambda2) * dTStar) -
        exp((dLambda + dLambda2) * (dTStar - dPWTime[0]))) /
           (dLambda + dLambda2) -
       exp(dLambda * (dTStar - dBondMat)) *
           (exp(dLambda2 * dTStar) - exp(dLambda2 * (dTStar - dPWTime[0]))) /
           dLambda2 -
       exp(dLambda2 * (dTStar - dBondMat)) *
           (exp(dLambda * dTStar) - exp(dLambda * (dTStar - dPWTime[0]))) /
           dLambda +
       exp((dLambda + dLambda2) * (dTStar - dBondMat)) * dPWTime[0]);

  dPWCSigma[0] = sqrt((ZCVar1 + ZCVar2 + 2 * ZCCoVar) / dPWTime[0]);

  ZCVol1 = (exp(dLambda * (dTStar - dPWTime[0])) -
            exp(dLambda * (dTStar - dBondMat))) *
           dSigma[0] / dLambda;
  ZCVol1 += (exp(dLambda * dTStar) - exp(dLambda * (dTStar - dBondMat))) *
            dSigma[0] / dLambda;
  ZCVol1 = 0.5 * ZCVol1;

  ZCVol2 = (exp(dLambda2 * (dTStar - dPWTime[0])) -
            exp(dLambda2 * (dTStar - dBondMat))) *
           dAlphaLGM[0] * dSigma[0] / dLambda2;
  ZCVol2 += (exp(dLambda2 * dTStar) - exp(dLambda2 * (dTStar - dBondMat))) *
            dAlphaLGM[0] * dSigma[0] / dLambda2;
  ZCVol2 = 0.5 * ZCVol2;

  ZCVol = sqrt(ZCVol1 * ZCVol1 + ZCVol2 * ZCVol2 +
               2 * dRhoLGM[0] * ZCVol1 * ZCVol2);

  dPWCRho[0] = sign * (ZCVol1 * dRho[0] + ZCVol2 * dRho2[0]) / ZCVol;

  for (i = 1; i < FixTimeIndex; ++i) {
    ZCVar1 = (dSigma[i] / dLambda) * (dSigma[i] / dLambda) *
             (0.5 *
                  (exp(2 * dLambda * (dTStar - dPWTime[i - 1])) -
                   exp(2 * dLambda * (dTStar - dPWTime[i]))) /
                  dLambda -
              2.0 * exp(dLambda * (dTStar - dBondMat)) *
                  (exp(dLambda * (dTStar - dPWTime[i - 1])) -
                   exp(dLambda * (dTStar - dPWTime[i]))) /
                  dLambda +
              exp(2 * dLambda * (dTStar - dBondMat)) *
                  (dPWTime[i] - dPWTime[i - 1]));

    ZCVar2 = (dAlphaLGM[i] * dSigma[i] / dLambda2) *
             (dAlphaLGM[i] * dSigma[i] / dLambda2) *
             (0.5 *
                  (exp(2 * dLambda2 * (dTStar - dPWTime[i - 1])) -
                   exp(2 * dLambda2 * (dTStar - dPWTime[i]))) /
                  dLambda2 -
              2.0 * exp(dLambda2 * (dTStar - dBondMat)) *
                  (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
                   exp(dLambda2 * (dTStar - dPWTime[i]))) /
                  dLambda2 +
              exp(2 * dLambda2 * (dTStar - dBondMat)) *
                  (dPWTime[i] - dPWTime[i - 1]));

    ZCCoVar = dAlphaLGM[i] * dRhoLGM[i] * (dSigma[i] / dLambda) *
              (dSigma[i] / dLambda2) *
              ((exp((dLambda + dLambda2) * (dTStar - dPWTime[i - 1])) -
                exp((dLambda + dLambda2) * (dTStar - dPWTime[i]))) /
                   (dLambda + dLambda2) -
               exp(dLambda * (dTStar - dBondMat)) *
                   (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
                    exp(dLambda2 * (dTStar - dPWTime[i]))) /
                   dLambda2 -
               exp(dLambda2 * (dTStar - dBondMat)) *
                   (exp(dLambda * (dTStar - dPWTime[i - 1])) -
                    exp(dLambda * (dTStar - dPWTime[i]))) /
                   dLambda +
               exp((dLambda + dLambda2) * (dTStar - dBondMat)) *
                   (dPWTime[i] - dPWTime[i - 1]));

    dPWCSigma[i] =
        sqrt((ZCVar1 + ZCVar2 + 2 * ZCCoVar) / (dPWTime[i] - dPWTime[i - 1]));

    ZCVol1 = (exp(dLambda * (dTStar - dPWTime[i])) -
              exp(dLambda * (dTStar - dBondMat))) *
             dSigma[i] / dLambda;
    ZCVol1 += (exp(dLambda * (dTStar - dPWTime[i - 1])) -
               exp(dLambda * (dTStar - dBondMat))) *
              dSigma[i] / dLambda;
    ZCVol1 = 0.5 * ZCVol1;

    ZCVol2 = (exp(dLambda2 * (dTStar - dPWTime[i])) -
              exp(dLambda2 * (dTStar - dBondMat))) *
             dAlphaLGM[i] * dSigma[i] / dLambda2;
    ZCVol2 += (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
               exp(dLambda2 * (dTStar - dBondMat))) *
              dAlphaLGM[i] * dSigma[i] / dLambda2;
    ZCVol2 = 0.5 * ZCVol2;

    ZCVol = sqrt(ZCVol1 * ZCVol1 + ZCVol2 * ZCVol2 +
                 2 * dRhoLGM[i] * ZCVol1 * ZCVol2);

    dPWCRho[i] = sign * (ZCVol1 * dRho[i] + ZCVol2 * dRho2[i]) / ZCVol;
  }

  ZCVar1 =
      (dSigma[i] / dLambda) * (dSigma[i] / dLambda) *
      (0.5 *
           (exp(2 * dLambda * (dTStar - dPWTime[i - 1])) -
            exp(2 * dLambda * (dTStar - dFixTime))) /
           dLambda -
       2.0 * exp(dLambda * (dTStar - dBondMat)) *
           (exp(dLambda * (dTStar - dPWTime[i - 1])) -
            exp(dLambda * (dTStar - dFixTime))) /
           dLambda +
       exp(2 * dLambda * (dTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  ZCVar2 =
      (dAlphaLGM[i] * dSigma[i] / dLambda2) *
      (dAlphaLGM[i] * dSigma[i] / dLambda2) *
      (0.5 *
           (exp(2 * dLambda2 * (dTStar - dPWTime[i - 1])) -
            exp(2 * dLambda2 * (dTStar - dFixTime))) /
           dLambda2 -
       2.0 * exp(dLambda2 * (dTStar - dBondMat)) *
           (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
            exp(dLambda2 * (dTStar - dFixTime))) /
           dLambda2 +
       exp(2 * dLambda2 * (dTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  ZCCoVar = dAlphaLGM[i] * dRhoLGM[i] * (dSigma[i] / dLambda) *
            (dSigma[i] / dLambda2) *
            ((exp((dLambda + dLambda2) * (dTStar - dPWTime[i - 1])) -
              exp((dLambda + dLambda2) * (dTStar - dFixTime))) /
                 (dLambda + dLambda2) -
             exp(dLambda * (dTStar - dBondMat)) *
                 (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
                  exp(dLambda2 * (dTStar - dFixTime))) /
                 dLambda2 -
             exp(dLambda2 * (dTStar - dBondMat)) *
                 (exp(dLambda * (dTStar - dPWTime[i - 1])) -
                  exp(dLambda * (dTStar - dFixTime))) /
                 dLambda +
             exp((dLambda + dLambda2) * (dTStar - dBondMat)) *
                 (dFixTime - dPWTime[i - 1]));

  dPWCSigma[i] =
      sqrt((ZCVar1 + ZCVar2 + 2 * ZCCoVar) / (dFixTime - dPWTime[i - 1]));

  ZCVol1 = (exp(dLambda * (dTStar - dFixTime)) -
            exp(dLambda * (dTStar - dBondMat))) *
           dSigma[i] / dLambda;
  ZCVol1 += (exp(dLambda * (dTStar - dPWTime[i - 1])) -
             exp(dLambda * (dTStar - dBondMat))) *
            dSigma[i] / dLambda;
  ZCVol1 = 0.5 * ZCVol1;

  ZCVol2 = (exp(dLambda2 * (dTStar - dFixTime)) -
            exp(dLambda2 * (dTStar - dBondMat))) *
           dAlphaLGM[i] * dSigma[i] / dLambda2;
  ZCVol2 += (exp(dLambda2 * (dTStar - dPWTime[i - 1])) -
             exp(dLambda2 * (dTStar - dBondMat))) *
            dAlphaLGM[i] * dSigma[i] / dLambda2;
  ZCVol2 = 0.5 * ZCVol2;

  ZCVol = sqrt(ZCVol1 * ZCVol1 + ZCVol2 * ZCVol2 +
               2 * dRhoLGM[i] * ZCVol1 * ZCVol2);

  dPWCRho[i] = sign * (ZCVol1 * dRho[i] + ZCVol2 * dRho2[i]) / ZCVol;
}

Err LGMSVBondVolApprox(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime, double *dSigma, double *dAlpha, double *dLevelEps,
    double *dLambdaEps, double *dRho, double *dRho2,

    double dLambda, double dTStar,

    double *dAlphaLGM, double *dRhoLGM, double dGammaLGM,

    // Product description
    int sign, double dFixTime, double dBondMat,

    /* Numerical parameters */
    int iNbX, double dIntegParam,

    // Outputs
    double *dBondATMVol) {
  Err err = NULL;
  double *dPWCSigma = NULL;
  double *dPWCRho = NULL;
  double dFwdPrice;

  dPWCSigma = (double *)calloc(iNbPWTime, sizeof(double));
  if (!dPWCSigma) {
    err = "allocation failed in LGMSVBondVolApprox";
    return err;
  }

  dPWCRho = (double *)calloc(iNbPWTime, sizeof(double));
  if (!dPWCRho) {
    err = "allocation failed in LGMSVBondVolApprox";
    return err;
  }

  ComputePWCBondVolAndCorrel(iNbPWTime, dPWTime, dSigma, dAlpha, dRho, dRho2,

                             sign,

                             dLambda, dTStar,

                             dAlphaLGM, dRhoLGM, dGammaLGM,

                             dBondMat, dFixTime,

                             dPWCSigma, dPWCRho);

  HestonClosedFormApprox(iNbPWTime, dPWTime, dPWCSigma, dAlpha, dLevelEps,
                         dLambdaEps, dPWCRho,

                         1, 1, dFixTime,

                         iNbX, dIntegParam,

                         &dFwdPrice);

  err = srt_f_optimpvol(dFwdPrice, 1, 1, dFixTime, 1, SRT_CALL, SRT_LOGNORMAL,
                        dBondATMVol);
  if (err) {
    return err;
  }

  if (dPWCSigma)
    free(dPWCSigma);
  if (dPWCRho)
    free(dPWCRho);

  return err;
}

void Compute5FCovariances(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    double *ddomSigma, double ddomLambda, double *ddomAlphaLGM,
    double *ddomRhoLGM, double ddomGammaLGM, double ddomTStar,

    double *dforSigma, double dforLambda, double *dforAlphaLGM,
    double *dforRhoLGM, double dforGammaLGM, double dforTStar,

    double *dfxSigma,

    double ***dCorrelation, // 0 : Dom1
                            // 1 : Dom2
                            // 2 : For1
                            // 3 : For2
                            // 4 : Fx

    double dBondMat, double dFixTime,

    double *ddomVar, double *dforVar, double *dfxVar, double *ddomfxCov,
    double *dforfxCov, double *ddomforCov) {
  int i, FixTimeIndex;
  double domZCVar, domZCVar1, domZCVar2, domZCCoVar;
  double ddomLambda2;
  double forZCVar, forZCVar1, forZCVar2, forZCCoVar;
  double dforLambda2;
  double fxVar;
  double dfxdomCoVar, dfxdomCoVar1, dfxdomCoVar2;
  double dfxforCoVar, dfxforCoVar1, dfxforCoVar2;
  double ddomforCoVar, ddomforCoVar11, ddomforCoVar12, ddomforCoVar21,
      ddomforCoVar22;

  FixTimeIndex = Get_Index(dFixTime, dPWTime, iNbPWTime);

  //	Domestic Bond Variance
  ddomLambda2 = ddomLambda + ddomGammaLGM;

  domZCVar1 = (ddomSigma[0] / ddomLambda) * (ddomSigma[0] / ddomLambda) *
              (0.5 *
                   (exp(2 * ddomLambda * ddomTStar) -
                    exp(2 * ddomLambda * (ddomTStar - dPWTime[0]))) /
                   ddomLambda -
               2.0 * exp(ddomLambda * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda * ddomTStar) -
                    exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
                   ddomLambda +
               exp(2 * ddomLambda * (ddomTStar - dBondMat)) * dPWTime[0]);

  domZCVar2 = (ddomAlphaLGM[0] * ddomSigma[0] / ddomLambda2) *
              (ddomAlphaLGM[0] * ddomSigma[0] / ddomLambda2) *
              (0.5 *
                   (exp(2 * ddomLambda2 * ddomTStar) -
                    exp(2 * ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                   ddomLambda2 -
               2.0 * exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda2 * ddomTStar) -
                    exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                   ddomLambda2 +
               exp(2 * ddomLambda2 * (ddomTStar - dBondMat)) * dPWTime[0]);

  domZCCoVar =
      ddomAlphaLGM[0] * ddomRhoLGM[0] * (ddomSigma[0] / ddomLambda) *
      (ddomSigma[0] / ddomLambda2) *
      ((exp((ddomLambda + ddomLambda2) * ddomTStar) -
        exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[0]))) /
           (ddomLambda + ddomLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(ddomLambda2 * ddomTStar) -
            exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(ddomLambda * ddomTStar) -
            exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
           ddomLambda +
       exp((ddomLambda + ddomLambda2) * (ddomTStar - dBondMat)) * dPWTime[0]);

  domZCVar = domZCVar1 + domZCVar2 + 2 * domZCCoVar;

  //	Foreign Bond Variance
  dforLambda2 = dforLambda + dforGammaLGM;

  forZCVar1 = (dforSigma[0] / dforLambda) * (dforSigma[0] / dforLambda) *
              (0.5 *
                   (exp(2 * dforLambda * dforTStar) -
                    exp(2 * dforLambda * (dforTStar - dPWTime[0]))) /
                   dforLambda -
               2.0 * exp(dforLambda * (dforTStar - dBondMat)) *
                   (exp(dforLambda * dforTStar) -
                    exp(dforLambda * (dforTStar - dPWTime[0]))) /
                   dforLambda +
               exp(2 * dforLambda * (dforTStar - dBondMat)) * dPWTime[0]);

  forZCVar2 = (dforAlphaLGM[0] * dforSigma[0] / dforLambda2) *
              (dforAlphaLGM[0] * dforSigma[0] / dforLambda2) *
              (0.5 *
                   (exp(2 * dforLambda2 * dforTStar) -
                    exp(2 * dforLambda2 * (dforTStar - dPWTime[0]))) /
                   dforLambda2 -
               2.0 * exp(dforLambda2 * (dforTStar - dBondMat)) *
                   (exp(dforLambda2 * dforTStar) -
                    exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
                   dforLambda2 +
               exp(2 * dforLambda2 * (dforTStar - dBondMat)) * dPWTime[0]);

  forZCCoVar =
      dforAlphaLGM[0] * dforRhoLGM[0] * (dforSigma[0] / dforLambda) *
      (dforSigma[0] / dforLambda2) *
      ((exp((dforLambda + dforLambda2) * dforTStar) -
        exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[0]))) /
           (dforLambda + dforLambda2) -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(dforLambda2 * dforTStar) -
            exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(dforLambda * dforTStar) -
            exp(dforLambda * (dforTStar - dPWTime[0]))) /
           dforLambda +
       exp((dforLambda + dforLambda2) * (dforTStar - dBondMat)) * dPWTime[0]);

  forZCVar = forZCVar1 + forZCVar2 + 2 * forZCCoVar;

  // FX Variance
  fxVar = dfxSigma[0] * dfxSigma[0] * dPWTime[0];

  // FX-Dom Covariance
  dfxdomCoVar1 = -dCorrelation[0][0][4] * dfxSigma[0] * ddomSigma[0] /
                 ddomLambda *
                 ((exp(ddomLambda * ddomTStar) -
                   exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
                      ddomLambda -
                  exp(ddomLambda * (ddomTStar - dBondMat)) * dPWTime[0]);

  dfxdomCoVar2 = -dCorrelation[0][1][4] * dfxSigma[0] * ddomAlphaLGM[0] *
                 ddomSigma[0] / ddomLambda2 *
                 ((exp(ddomLambda2 * ddomTStar) -
                   exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                      ddomLambda2 -
                  exp(ddomLambda2 * (ddomTStar - dBondMat)) * dPWTime[0]);

  dfxdomCoVar = dfxdomCoVar1 + dfxdomCoVar2;

  // FX-For Covariance
  dfxforCoVar1 = -dCorrelation[0][2][4] * dfxSigma[0] * dforSigma[0] /
                 dforLambda *
                 ((exp(dforLambda * dforTStar) -
                   exp(dforLambda * (dforTStar - dPWTime[0]))) /
                      dforLambda -
                  exp(dforLambda * (dforTStar - dBondMat)) * dPWTime[0]);

  dfxforCoVar2 = -dCorrelation[0][3][4] * dfxSigma[0] * dforAlphaLGM[0] *
                 dforSigma[0] / dforLambda2 *
                 ((exp(dforLambda2 * dforTStar) -
                   exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
                      dforLambda2 -
                  exp(dforLambda2 * (dforTStar - dBondMat)) * dPWTime[0]);

  dfxforCoVar = dfxforCoVar1 + dfxforCoVar2;

  // Dom-For Covariance
  ddomforCoVar11 =
      dCorrelation[0][0][2] * ddomSigma[0] * dforSigma[0] /
      (dforLambda * ddomLambda) *
      (-(exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda + dforLambda) * ddomTStar)) /
           (ddomLambda + dforLambda) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda * dforTStar) -
            exp(dforLambda * (dforTStar - dPWTime[0]))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda * ddomTStar) -
            exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
           ddomLambda +
       exp((ddomLambda + dforLambda) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar12 =
      dCorrelation[0][0][3] * ddomSigma[0] * dforSigma[0] * dforAlphaLGM[0] /
      (dforLambda2 * ddomLambda) *
      (-(exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda + dforLambda2) * ddomTStar)) /
           (ddomLambda + dforLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * dforTStar) -
            exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda * ddomTStar) -
            exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
           ddomLambda +
       exp((ddomLambda + dforLambda2) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar21 =
      dCorrelation[0][1][2] * ddomSigma[0] * ddomAlphaLGM[0] * dforSigma[0] /
      (dforLambda * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda2 + dforLambda) * ddomTStar)) /
           (ddomLambda2 + dforLambda) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda * dforTStar) -
            exp(dforLambda * (dforTStar - dPWTime[0]))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * ddomTStar) -
            exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar22 =
      dCorrelation[0][1][3] * ddomSigma[0] * ddomAlphaLGM[0] * dforSigma[0] *
      dforAlphaLGM[0] / (dforLambda2 * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda2 + dforLambda2) * ddomTStar)) /
           (ddomLambda2 + dforLambda2) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * dforTStar) -
            exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * ddomTStar) -
            exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda2) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar =
      ddomforCoVar11 + ddomforCoVar12 + ddomforCoVar21 + ddomforCoVar22;

  for (i = 1; i < FixTimeIndex; ++i) {
    //	Domestic Bond Variance
    domZCVar1 = (ddomSigma[i] / ddomLambda) * (ddomSigma[i] / ddomLambda) *
                (0.5 *
                     (exp(2 * ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                      exp(2 * ddomLambda * (ddomTStar - dPWTime[i]))) /
                     ddomLambda -
                 2.0 * exp(ddomLambda * (ddomTStar - dBondMat)) *
                     (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                      exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
                     ddomLambda +
                 exp(2 * ddomLambda * (ddomTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    domZCVar2 = (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
                (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
                (0.5 *
                     (exp(2 * ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                      exp(2 * ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                     ddomLambda2 -
                 2.0 * exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                     (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                      exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                     ddomLambda2 +
                 exp(2 * ddomLambda2 * (ddomTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    domZCCoVar =
        ddomAlphaLGM[i] * ddomRhoLGM[i] * (ddomSigma[i] / ddomLambda) *
        (ddomSigma[i] / ddomLambda2) *
        ((exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[i - 1])) -
          exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[i]))) /
             (ddomLambda + ddomLambda2) -
         exp(ddomLambda * (ddomTStar - dBondMat)) *
             (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
             ddomLambda2 -
         exp(ddomLambda2 * (ddomTStar - dBondMat)) *
             (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
             ddomLambda +
         exp((ddomLambda + ddomLambda2) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    domZCVar += domZCVar1 + domZCVar2 + 2 * domZCCoVar;

    //	Foreign Bond Variance
    forZCVar1 = (dforSigma[i] / dforLambda) * (dforSigma[i] / dforLambda) *
                (0.5 *
                     (exp(2 * dforLambda * (dforTStar - dPWTime[i - 1])) -
                      exp(2 * dforLambda * (dforTStar - dPWTime[i]))) /
                     dforLambda -
                 2.0 * exp(dforLambda * (dforTStar - dBondMat)) *
                     (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                      exp(dforLambda * (dforTStar - dPWTime[i]))) /
                     dforLambda +
                 exp(2 * dforLambda * (dforTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    forZCVar2 = (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
                (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
                (0.5 *
                     (exp(2 * dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                      exp(2 * dforLambda2 * (dforTStar - dPWTime[i]))) /
                     dforLambda2 -
                 2.0 * exp(dforLambda2 * (dforTStar - dBondMat)) *
                     (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                      exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
                     dforLambda2 +
                 exp(2 * dforLambda2 * (dforTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    forZCCoVar =
        dforAlphaLGM[i] * dforRhoLGM[i] * (dforSigma[i] / dforLambda) *
        (dforSigma[i] / dforLambda2) *
        ((exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[i - 1])) -
          exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[i]))) /
             (dforLambda + dforLambda2) -
         exp(dforLambda * (dforTStar - dBondMat)) *
             (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
             dforLambda2 -
         exp(dforLambda2 * (dforTStar - dBondMat)) *
             (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda * (dforTStar - dPWTime[i]))) /
             dforLambda +
         exp((dforLambda + dforLambda2) * (dforTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    forZCVar += forZCVar1 + forZCVar2 + 2 * forZCCoVar;

    // FX Variance
    fxVar += dfxSigma[i] * dfxSigma[i] * (dPWTime[i] - dPWTime[i - 1]);

    // FX-Dom Covariance
    dfxdomCoVar1 = -dCorrelation[i][0][4] * dfxSigma[i] * ddomSigma[i] /
                   ddomLambda *
                   ((exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                     exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
                        ddomLambda -
                    exp(ddomLambda * (ddomTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxdomCoVar2 = -dCorrelation[i][1][4] * dfxSigma[i] * ddomAlphaLGM[i] *
                   ddomSigma[i] / ddomLambda2 *
                   ((exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                     exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                        ddomLambda2 -
                    exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxdomCoVar += dfxdomCoVar1 + dfxdomCoVar2;

    // FX-For Covariance
    dfxforCoVar1 = -dCorrelation[i][2][4] * dfxSigma[i] * dforSigma[i] /
                   dforLambda *
                   ((exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                     exp(dforLambda * (dforTStar - dPWTime[i]))) /
                        dforLambda -
                    exp(dforLambda * (dforTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxforCoVar2 = -dCorrelation[i][3][4] * dfxSigma[i] * dforAlphaLGM[i] *
                   dforSigma[i] / dforLambda2 *
                   ((exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                     exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
                        dforLambda2 -
                    exp(dforLambda2 * (dforTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxforCoVar += dfxforCoVar1 + dfxforCoVar2;

    // Dom-For Covariance
    ddomforCoVar11 =
        dCorrelation[i][0][2] * ddomSigma[i] * dforSigma[i] /
        (dforLambda * ddomLambda) *
        (-(exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda + dforLambda) -
         exp(ddomLambda * (ddomTStar - dBondMat)) *
             (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda * (dforTStar - dPWTime[i]))) /
             dforLambda -
         exp(dforLambda * (dforTStar - dBondMat)) *
             (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
             ddomLambda +
         exp((ddomLambda + dforLambda) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar12 =
        dCorrelation[i][0][3] * ddomSigma[i] * dforSigma[i] * dforAlphaLGM[i] /
        (dforLambda2 * ddomLambda) *
        (-(exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda + dforLambda2) -
         exp(ddomLambda * (ddomTStar - dBondMat)) *
             (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
             dforLambda2 -
         exp(dforLambda2 * (dforTStar - dBondMat)) *
             (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
             ddomLambda +
         exp((ddomLambda + dforLambda2) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar21 =
        dCorrelation[i][1][2] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] /
        (dforLambda * ddomLambda2) *
        (-(exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda2 + dforLambda) -
         exp(ddomLambda2 * (ddomTStar - dBondMat)) *
             (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda * (dforTStar - dPWTime[i]))) /
             dforLambda -
         exp(dforLambda * (dforTStar - dBondMat)) *
             (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
             ddomLambda2 +
         exp((ddomLambda2 + dforLambda) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar22 =
        dCorrelation[i][1][3] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] *
        dforAlphaLGM[i] / (dforLambda2 * ddomLambda2) *
        (-(exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda2 + dforLambda2) -
         exp(ddomLambda2 * (ddomTStar - dBondMat)) *
             (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
             dforLambda2 -
         exp(dforLambda2 * (dforTStar - dBondMat)) *
             (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
             ddomLambda2 +
         exp((ddomLambda2 + dforLambda2) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar +=
        ddomforCoVar11 + ddomforCoVar12 + ddomforCoVar21 + ddomforCoVar22;
  }

  //	Domestic Bond Variance
  domZCVar1 = (ddomSigma[i] / ddomLambda) * (ddomSigma[i] / ddomLambda) *
              (0.5 *
                   (exp(2 * ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                    exp(2 * ddomLambda * (ddomTStar - dFixTime))) /
                   ddomLambda -
               2.0 * exp(ddomLambda * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                    exp(ddomLambda * (ddomTStar - dFixTime))) /
                   ddomLambda +
               exp(2 * ddomLambda * (ddomTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  domZCVar2 = (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
              (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
              (0.5 *
                   (exp(2 * ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                    exp(2 * ddomLambda2 * (ddomTStar - dFixTime))) /
                   ddomLambda2 -
               2.0 * exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                    exp(ddomLambda2 * (ddomTStar - dFixTime))) /
                   ddomLambda2 +
               exp(2 * ddomLambda2 * (ddomTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  domZCCoVar =
      ddomAlphaLGM[i] * ddomRhoLGM[i] * (ddomSigma[i] / ddomLambda) *
      (ddomSigma[i] / ddomLambda2) *
      ((exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[i - 1])) -
        exp((ddomLambda + ddomLambda2) * (ddomTStar - dFixTime))) /
           (ddomLambda + ddomLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda +
       exp((ddomLambda + ddomLambda2) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  domZCVar += domZCVar1 + domZCVar2 + 2 * domZCCoVar;

  //	Foreign Bond Variance
  forZCVar1 = (dforSigma[i] / dforLambda) * (dforSigma[i] / dforLambda) *
              (0.5 *
                   (exp(2 * dforLambda * (dforTStar - dPWTime[i - 1])) -
                    exp(2 * dforLambda * (dforTStar - dFixTime))) /
                   dforLambda -
               2.0 * exp(dforLambda * (dforTStar - dBondMat)) *
                   (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                    exp(dforLambda * (dforTStar - dFixTime))) /
                   dforLambda +
               exp(2 * dforLambda * (dforTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  forZCVar2 = (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
              (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
              (0.5 *
                   (exp(2 * dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                    exp(2 * dforLambda2 * (dforTStar - dFixTime))) /
                   dforLambda2 -
               2.0 * exp(dforLambda2 * (dforTStar - dBondMat)) *
                   (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                    exp(dforLambda2 * (dforTStar - dFixTime))) /
                   dforLambda2 +
               exp(2 * dforLambda2 * (dforTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  forZCCoVar =
      dforAlphaLGM[i] * dforRhoLGM[i] * (dforSigma[i] / dforLambda) *
      (dforSigma[i] / dforLambda2) *
      ((exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[i - 1])) -
        exp((dforLambda + dforLambda2) * (dforTStar - dFixTime))) /
           (dforLambda + dforLambda2) -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda +
       exp((dforLambda + dforLambda2) * (dforTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  forZCVar += forZCVar1 + forZCVar2 + 2 * forZCCoVar;

  // FX Variance
  fxVar += dfxSigma[i] * dfxSigma[i] * (dFixTime - dPWTime[i - 1]);

  // FX-Dom Covariance
  dfxdomCoVar1 =
      -dCorrelation[i][0][4] * dfxSigma[i] * ddomSigma[i] / ddomLambda *
      ((exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
        exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda -
       exp(ddomLambda * (ddomTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxdomCoVar2 =
      -dCorrelation[i][1][4] * dfxSigma[i] * ddomAlphaLGM[i] * ddomSigma[i] /
      ddomLambda2 *
      ((exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
        exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxdomCoVar += dfxdomCoVar1 + dfxdomCoVar2;

  // FX-For Covariance
  dfxforCoVar1 =
      -dCorrelation[i][2][4] * dfxSigma[i] * dforSigma[i] / dforLambda *
      ((exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
        exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxforCoVar2 =
      -dCorrelation[i][3][4] * dfxSigma[i] * dforAlphaLGM[i] * dforSigma[i] /
      dforLambda2 *
      ((exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
        exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxforCoVar += dfxforCoVar1 + dfxforCoVar2;

  // Dom-For Covariance
  ddomforCoVar11 =
      dCorrelation[i][0][2] * ddomSigma[i] * dforSigma[i] /
      (dforLambda * ddomLambda) *
      (-(exp((ddomLambda + dforLambda) * (ddomTStar - dFixTime)) -
         exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda + dforLambda) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda +
       exp((ddomLambda + dforLambda) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar12 =
      dCorrelation[i][0][3] * ddomSigma[i] * dforSigma[i] * dforAlphaLGM[i] /
      (dforLambda2 * ddomLambda) *
      (-(exp((ddomLambda + dforLambda2) * (ddomTStar - dFixTime)) -
         exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda + dforLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda +
       exp((ddomLambda + dforLambda2) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar21 =
      dCorrelation[i][1][2] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] /
      (dforLambda * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda) * (ddomTStar - dFixTime)) -
         exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda2 + dforLambda) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar22 =
      dCorrelation[i][1][3] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] *
      dforAlphaLGM[i] / (dforLambda2 * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda2) * (ddomTStar - dFixTime)) -
         exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda2 + dforLambda2) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda2) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar +=
      ddomforCoVar11 + ddomforCoVar12 + ddomforCoVar21 + ddomforCoVar22;

  *ddomVar = domZCVar;
  *dforVar = forZCVar;
  *dfxVar = fxVar;
  *ddomfxCov = dfxdomCoVar;
  *dforfxCov = dfxforCoVar;
  *ddomforCov = ddomforCoVar;
}

void Compute5FCovariancesVect(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    double *ddomSigma, double ddomLambda, double *ddomAlphaLGM,
    double *ddomRhoLGM, double ddomGammaLGM, double ddomTStar,

    double *dforSigma, double dforLambda, double *dforAlphaLGM,
    double *dforRhoLGM, double dforGammaLGM, double dforTStar,

    double *dfxSigma,

    double *dAlpha,

    double ***dCorrelation, // 0 : Dom1
                            // 1 : Dom2
                            // 2 : For1
                            // 3 : For2
                            // 4 : Fx
                            // 5 : Alpha

    double dBondMat, double dFixTime,

    double *ddomVar, double *dforVar, double *dfxVar, double *dAlphaVariance,
    double *ddomfxCov, double *dforfxCov, double *ddomforCov,
    double *dAlphafxCov, double *dAlphadomCov, double *dAlphaforCov) {
  int i, FixTimeIndex;
  double domZCVar, domZCVar1, domZCVar2, domZCCoVar;
  double ddomLambda2;
  double forZCVar, forZCVar1, forZCVar2, forZCCoVar;
  double dforLambda2;
  double fxVar;
  double dfxdomCoVar, dfxdomCoVar1, dfxdomCoVar2;
  double dfxforCoVar, dfxforCoVar1, dfxforCoVar2;
  double ddomforCoVar, ddomforCoVar11, ddomforCoVar12, ddomforCoVar21,
      ddomforCoVar22;
  double dAlphadomCoVar, dAlphadomCoVar1, dAlphadomCoVar2;
  double dAlphaforCoVar, dAlphaforCoVar1, dAlphaforCoVar2;
  double dAlphafxCoVar;
  double dAlphaVar;

  FixTimeIndex = Get_Index(dFixTime, dPWTime, iNbPWTime);

  //	Domestic Bond Variance
  ddomLambda2 = ddomLambda + ddomGammaLGM;

  domZCVar1 = (ddomSigma[0] / ddomLambda) * (ddomSigma[0] / ddomLambda) *
              (0.5 *
                   (exp(2 * ddomLambda * ddomTStar) -
                    exp(2 * ddomLambda * (ddomTStar - dPWTime[0]))) /
                   ddomLambda -
               2.0 * exp(ddomLambda * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda * ddomTStar) -
                    exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
                   ddomLambda +
               exp(2 * ddomLambda * (ddomTStar - dBondMat)) * dPWTime[0]);

  domZCVar2 = (ddomAlphaLGM[0] * ddomSigma[0] / ddomLambda2) *
              (ddomAlphaLGM[0] * ddomSigma[0] / ddomLambda2) *
              (0.5 *
                   (exp(2 * ddomLambda2 * ddomTStar) -
                    exp(2 * ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                   ddomLambda2 -
               2.0 * exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda2 * ddomTStar) -
                    exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                   ddomLambda2 +
               exp(2 * ddomLambda2 * (ddomTStar - dBondMat)) * dPWTime[0]);

  domZCCoVar =
      ddomAlphaLGM[0] * ddomRhoLGM[0] * (ddomSigma[0] / ddomLambda) *
      (ddomSigma[0] / ddomLambda2) *
      ((exp((ddomLambda + ddomLambda2) * ddomTStar) -
        exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[0]))) /
           (ddomLambda + ddomLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(ddomLambda2 * ddomTStar) -
            exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(ddomLambda * ddomTStar) -
            exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
           ddomLambda +
       exp((ddomLambda + ddomLambda2) * (ddomTStar - dBondMat)) * dPWTime[0]);

  domZCVar = domZCVar1 + domZCVar2 + 2 * domZCCoVar;
  ddomVar[0] = domZCVar;

  //	Foreign Bond Variance
  dforLambda2 = dforLambda + dforGammaLGM;

  forZCVar1 = (dforSigma[0] / dforLambda) * (dforSigma[0] / dforLambda) *
              (0.5 *
                   (exp(2 * dforLambda * dforTStar) -
                    exp(2 * dforLambda * (dforTStar - dPWTime[0]))) /
                   dforLambda -
               2.0 * exp(dforLambda * (dforTStar - dBondMat)) *
                   (exp(dforLambda * dforTStar) -
                    exp(dforLambda * (dforTStar - dPWTime[0]))) /
                   dforLambda +
               exp(2 * dforLambda * (dforTStar - dBondMat)) * dPWTime[0]);

  forZCVar2 = (dforAlphaLGM[0] * dforSigma[0] / dforLambda2) *
              (dforAlphaLGM[0] * dforSigma[0] / dforLambda2) *
              (0.5 *
                   (exp(2 * dforLambda2 * dforTStar) -
                    exp(2 * dforLambda2 * (dforTStar - dPWTime[0]))) /
                   dforLambda2 -
               2.0 * exp(dforLambda2 * (dforTStar - dBondMat)) *
                   (exp(dforLambda2 * dforTStar) -
                    exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
                   dforLambda2 +
               exp(2 * dforLambda2 * (dforTStar - dBondMat)) * dPWTime[0]);

  forZCCoVar =
      dforAlphaLGM[0] * dforRhoLGM[0] * (dforSigma[0] / dforLambda) *
      (dforSigma[0] / dforLambda2) *
      ((exp((dforLambda + dforLambda2) * dforTStar) -
        exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[0]))) /
           (dforLambda + dforLambda2) -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(dforLambda2 * dforTStar) -
            exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(dforLambda * dforTStar) -
            exp(dforLambda * (dforTStar - dPWTime[0]))) /
           dforLambda +
       exp((dforLambda + dforLambda2) * (dforTStar - dBondMat)) * dPWTime[0]);

  forZCVar = forZCVar1 + forZCVar2 + 2 * forZCCoVar;
  dforVar[0] = forZCVar;

  // FX Variance
  fxVar = dfxSigma[0] * dfxSigma[0] * dPWTime[0];
  dfxVar[0] = fxVar;

  // FX-Dom Covariance
  dfxdomCoVar1 = -dCorrelation[0][0][4] * dfxSigma[0] * ddomSigma[0] /
                 ddomLambda *
                 ((exp(ddomLambda * ddomTStar) -
                   exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
                      ddomLambda -
                  exp(ddomLambda * (ddomTStar - dBondMat)) * dPWTime[0]);

  dfxdomCoVar2 = -dCorrelation[0][1][4] * dfxSigma[0] * ddomAlphaLGM[0] *
                 ddomSigma[0] / ddomLambda2 *
                 ((exp(ddomLambda2 * ddomTStar) -
                   exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                      ddomLambda2 -
                  exp(ddomLambda2 * (ddomTStar - dBondMat)) * dPWTime[0]);

  dfxdomCoVar = dfxdomCoVar1 + dfxdomCoVar2;
  ddomfxCov[0] = dfxdomCoVar;

  // FX-For Covariance
  dfxforCoVar1 = -dCorrelation[0][2][4] * dfxSigma[0] * dforSigma[0] /
                 dforLambda *
                 ((exp(dforLambda * dforTStar) -
                   exp(dforLambda * (dforTStar - dPWTime[0]))) /
                      dforLambda -
                  exp(dforLambda * (dforTStar - dBondMat)) * dPWTime[0]);

  dfxforCoVar2 = -dCorrelation[0][3][4] * dfxSigma[0] * dforAlphaLGM[0] *
                 dforSigma[0] / dforLambda2 *
                 ((exp(dforLambda2 * dforTStar) -
                   exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
                      dforLambda2 -
                  exp(dforLambda2 * (dforTStar - dBondMat)) * dPWTime[0]);

  dfxforCoVar = dfxforCoVar1 + dfxforCoVar2;
  dforfxCov[0] = dfxforCoVar;

  // Dom-For Covariance
  ddomforCoVar11 =
      dCorrelation[0][0][2] * ddomSigma[0] * dforSigma[0] /
      (dforLambda * ddomLambda) *
      (-(exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda + dforLambda) * ddomTStar)) /
           (ddomLambda + dforLambda) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda * dforTStar) -
            exp(dforLambda * (dforTStar - dPWTime[0]))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda * ddomTStar) -
            exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
           ddomLambda +
       exp((ddomLambda + dforLambda) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar12 =
      dCorrelation[0][0][3] * ddomSigma[0] * dforSigma[0] * dforAlphaLGM[0] /
      (dforLambda2 * ddomLambda) *
      (-(exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda + dforLambda2) * ddomTStar)) /
           (ddomLambda + dforLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * dforTStar) -
            exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda * ddomTStar) -
            exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
           ddomLambda +
       exp((ddomLambda + dforLambda2) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar21 =
      dCorrelation[0][1][2] * ddomSigma[0] * ddomAlphaLGM[0] * dforSigma[0] /
      (dforLambda * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda2 + dforLambda) * ddomTStar)) /
           (ddomLambda2 + dforLambda) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda * dforTStar) -
            exp(dforLambda * (dforTStar - dPWTime[0]))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * ddomTStar) -
            exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar22 =
      dCorrelation[0][1][3] * ddomSigma[0] * ddomAlphaLGM[0] * dforSigma[0] *
      dforAlphaLGM[0] / (dforLambda2 * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[0])) -
         exp((ddomLambda2 + dforLambda2) * ddomTStar)) /
           (ddomLambda2 + dforLambda2) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * dforTStar) -
            exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * ddomTStar) -
            exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda2) * (ddomTStar - dBondMat)) * dPWTime[0]);

  ddomforCoVar =
      ddomforCoVar11 + ddomforCoVar12 + ddomforCoVar21 + ddomforCoVar22;
  ddomforCov[0] = ddomforCoVar;

  // Alpha-Dom Covariance
  dAlphadomCoVar1 = -dCorrelation[0][0][5] * dAlpha[0] * ddomSigma[0] /
                    ddomLambda *
                    ((exp(ddomLambda * ddomTStar) -
                      exp(ddomLambda * (ddomTStar - dPWTime[0]))) /
                         ddomLambda -
                     exp(ddomLambda * (ddomTStar - dBondMat)) * dPWTime[0]);

  dAlphadomCoVar2 = -dCorrelation[0][1][5] * dAlpha[0] * ddomAlphaLGM[0] *
                    ddomSigma[0] / ddomLambda2 *
                    ((exp(ddomLambda2 * ddomTStar) -
                      exp(ddomLambda2 * (ddomTStar - dPWTime[0]))) /
                         ddomLambda2 -
                     exp(ddomLambda2 * (ddomTStar - dBondMat)) * dPWTime[0]);

  dAlphadomCoVar = dAlphadomCoVar1 + dAlphadomCoVar2;
  dAlphadomCov[0] = dAlphadomCoVar;

  // Alpha-for Covariance
  dAlphaforCoVar1 = -dCorrelation[0][0][5] * dAlpha[0] * dforSigma[0] /
                    dforLambda *
                    ((exp(dforLambda * dforTStar) -
                      exp(dforLambda * (dforTStar - dPWTime[0]))) /
                         dforLambda -
                     exp(dforLambda * (dforTStar - dBondMat)) * dPWTime[0]);

  dAlphaforCoVar2 = -dCorrelation[0][1][5] * dAlpha[0] * dforAlphaLGM[0] *
                    dforSigma[0] / dforLambda2 *
                    ((exp(dforLambda2 * dforTStar) -
                      exp(dforLambda2 * (dforTStar - dPWTime[0]))) /
                         dforLambda2 -
                     exp(dforLambda2 * (dforTStar - dBondMat)) * dPWTime[0]);

  dAlphaforCoVar = dAlphaforCoVar1 + dAlphaforCoVar2;
  dAlphaforCov[0] = dAlphaforCoVar;

  // FX-Alpha Covariance
  dAlphafxCoVar = dfxSigma[0] * dAlpha[0] * dCorrelation[0][4][5] * dPWTime[0];
  dAlphafxCov[0] = dAlphafxCoVar;

  // Alpha Variance
  dAlphaVar = dAlpha[0] * dAlpha[0] * dPWTime[0];
  dAlphaVariance[0] = dAlphaVar;

  for (i = 1; i < FixTimeIndex; ++i) {
    //	Domestic Bond Variance
    domZCVar1 = (ddomSigma[i] / ddomLambda) * (ddomSigma[i] / ddomLambda) *
                (0.5 *
                     (exp(2 * ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                      exp(2 * ddomLambda * (ddomTStar - dPWTime[i]))) /
                     ddomLambda -
                 2.0 * exp(ddomLambda * (ddomTStar - dBondMat)) *
                     (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                      exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
                     ddomLambda +
                 exp(2 * ddomLambda * (ddomTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    domZCVar2 = (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
                (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
                (0.5 *
                     (exp(2 * ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                      exp(2 * ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                     ddomLambda2 -
                 2.0 * exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                     (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                      exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                     ddomLambda2 +
                 exp(2 * ddomLambda2 * (ddomTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    domZCCoVar =
        ddomAlphaLGM[i] * ddomRhoLGM[i] * (ddomSigma[i] / ddomLambda) *
        (ddomSigma[i] / ddomLambda2) *
        ((exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[i - 1])) -
          exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[i]))) /
             (ddomLambda + ddomLambda2) -
         exp(ddomLambda * (ddomTStar - dBondMat)) *
             (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
             ddomLambda2 -
         exp(ddomLambda2 * (ddomTStar - dBondMat)) *
             (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
             ddomLambda +
         exp((ddomLambda + ddomLambda2) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    domZCVar = domZCVar1 + domZCVar2 + 2 * domZCCoVar;
    ddomVar[i] = domZCVar;

    //	Foreign Bond Variance
    forZCVar1 = (dforSigma[i] / dforLambda) * (dforSigma[i] / dforLambda) *
                (0.5 *
                     (exp(2 * dforLambda * (dforTStar - dPWTime[i - 1])) -
                      exp(2 * dforLambda * (dforTStar - dPWTime[i]))) /
                     dforLambda -
                 2.0 * exp(dforLambda * (dforTStar - dBondMat)) *
                     (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                      exp(dforLambda * (dforTStar - dPWTime[i]))) /
                     dforLambda +
                 exp(2 * dforLambda * (dforTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    forZCVar2 = (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
                (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
                (0.5 *
                     (exp(2 * dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                      exp(2 * dforLambda2 * (dforTStar - dPWTime[i]))) /
                     dforLambda2 -
                 2.0 * exp(dforLambda2 * (dforTStar - dBondMat)) *
                     (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                      exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
                     dforLambda2 +
                 exp(2 * dforLambda2 * (dforTStar - dBondMat)) *
                     (dPWTime[i] - dPWTime[i - 1]));

    forZCCoVar =
        dforAlphaLGM[i] * dforRhoLGM[i] * (dforSigma[i] / dforLambda) *
        (dforSigma[i] / dforLambda2) *
        ((exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[i - 1])) -
          exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[i]))) /
             (dforLambda + dforLambda2) -
         exp(dforLambda * (dforTStar - dBondMat)) *
             (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
             dforLambda2 -
         exp(dforLambda2 * (dforTStar - dBondMat)) *
             (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda * (dforTStar - dPWTime[i]))) /
             dforLambda +
         exp((dforLambda + dforLambda2) * (dforTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    forZCVar = forZCVar1 + forZCVar2 + 2 * forZCCoVar;
    dforVar[i] = forZCVar;

    // FX Variance
    fxVar = dfxSigma[i] * dfxSigma[i] * (dPWTime[i] - dPWTime[i - 1]);
    dfxVar[i] = fxVar;

    // FX-Dom Covariance
    dfxdomCoVar1 = -dCorrelation[i][0][4] * dfxSigma[i] * ddomSigma[i] /
                   ddomLambda *
                   ((exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                     exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
                        ddomLambda -
                    exp(ddomLambda * (ddomTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxdomCoVar2 = -dCorrelation[i][1][4] * dfxSigma[i] * ddomAlphaLGM[i] *
                   ddomSigma[i] / ddomLambda2 *
                   ((exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                     exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                        ddomLambda2 -
                    exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxdomCoVar = dfxdomCoVar1 + dfxdomCoVar2;
    ddomfxCov[i] = dfxdomCoVar;

    // FX-For Covariance
    dfxforCoVar1 = -dCorrelation[i][2][4] * dfxSigma[i] * dforSigma[i] /
                   dforLambda *
                   ((exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                     exp(dforLambda * (dforTStar - dPWTime[i]))) /
                        dforLambda -
                    exp(dforLambda * (dforTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxforCoVar2 = -dCorrelation[i][3][4] * dfxSigma[i] * dforAlphaLGM[i] *
                   dforSigma[i] / dforLambda2 *
                   ((exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                     exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
                        dforLambda2 -
                    exp(dforLambda2 * (dforTStar - dBondMat)) *
                        (dPWTime[i] - dPWTime[i - 1]));

    dfxforCoVar = dfxforCoVar1 + dfxforCoVar2;
    dforfxCov[i] = dfxforCoVar;

    // Dom-For Covariance
    ddomforCoVar11 =
        dCorrelation[i][0][2] * ddomSigma[i] * dforSigma[i] /
        (dforLambda * ddomLambda) *
        (-(exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda + dforLambda) -
         exp(ddomLambda * (ddomTStar - dBondMat)) *
             (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda * (dforTStar - dPWTime[i]))) /
             dforLambda -
         exp(dforLambda * (dforTStar - dBondMat)) *
             (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
             ddomLambda +
         exp((ddomLambda + dforLambda) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar12 =
        dCorrelation[i][0][3] * ddomSigma[i] * dforSigma[i] * dforAlphaLGM[i] /
        (dforLambda2 * ddomLambda) *
        (-(exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda + dforLambda2) -
         exp(ddomLambda * (ddomTStar - dBondMat)) *
             (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
             dforLambda2 -
         exp(dforLambda2 * (dforTStar - dBondMat)) *
             (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
             ddomLambda +
         exp((ddomLambda + dforLambda2) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar21 =
        dCorrelation[i][1][2] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] /
        (dforLambda * ddomLambda2) *
        (-(exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda2 + dforLambda) -
         exp(ddomLambda2 * (ddomTStar - dBondMat)) *
             (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda * (dforTStar - dPWTime[i]))) /
             dforLambda -
         exp(dforLambda * (dforTStar - dBondMat)) *
             (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
             ddomLambda2 +
         exp((ddomLambda2 + dforLambda) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar22 =
        dCorrelation[i][1][3] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] *
        dforAlphaLGM[i] / (dforLambda2 * ddomLambda2) *
        (-(exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[i])) -
           exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
             (ddomLambda2 + dforLambda2) -
         exp(ddomLambda2 * (ddomTStar - dBondMat)) *
             (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
              exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
             dforLambda2 -
         exp(dforLambda2 * (dforTStar - dBondMat)) *
             (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
              exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
             ddomLambda2 +
         exp((ddomLambda2 + dforLambda2) * (ddomTStar - dBondMat)) *
             (dPWTime[i] - dPWTime[i - 1]));

    ddomforCoVar =
        ddomforCoVar11 + ddomforCoVar12 + ddomforCoVar21 + ddomforCoVar22;
    ddomforCov[i] = ddomforCoVar;

    // Alpha-Dom Covariance
    dAlphadomCoVar1 = -dCorrelation[i][0][5] * dAlpha[i] * ddomSigma[i] /
                      ddomLambda *
                      ((exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                        exp(ddomLambda * (ddomTStar - dPWTime[i]))) /
                           ddomLambda -
                       exp(ddomLambda * (ddomTStar - dBondMat)) *
                           (dPWTime[i] - dPWTime[i - 1]));

    dAlphadomCoVar2 = -dCorrelation[i][1][5] * dAlpha[i] * ddomAlphaLGM[i] *
                      ddomSigma[i] / ddomLambda2 *
                      ((exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                        exp(ddomLambda2 * (ddomTStar - dPWTime[i]))) /
                           ddomLambda2 -
                       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                           (dPWTime[i] - dPWTime[i - 1]));

    dAlphadomCoVar = dAlphadomCoVar1 + dAlphadomCoVar2;
    dAlphadomCov[i] = dAlphadomCoVar;

    // Alpha-for Covariance
    dAlphaforCoVar1 = -dCorrelation[i][0][5] * dAlpha[i] * dforSigma[i] /
                      dforLambda *
                      ((exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                        exp(dforLambda * (dforTStar - dPWTime[i]))) /
                           dforLambda -
                       exp(dforLambda * (dforTStar - dBondMat)) *
                           (dPWTime[i] - dPWTime[i - 1]));

    dAlphaforCoVar2 = -dCorrelation[i][1][5] * dAlpha[i] * dforAlphaLGM[i] *
                      dforSigma[i] / dforLambda2 *
                      ((exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                        exp(dforLambda2 * (dforTStar - dPWTime[i]))) /
                           dforLambda2 -
                       exp(dforLambda2 * (dforTStar - dBondMat)) *
                           (dPWTime[i] - dPWTime[i - 1]));

    dAlphaforCoVar = dAlphaforCoVar1 + dAlphaforCoVar2;
    dAlphaforCov[i] = dAlphaforCoVar;

    // FX-Alpha Covariance
    dAlphafxCoVar = dfxSigma[i] * dAlpha[i] * dCorrelation[i][4][5] *
                    (dPWTime[i] - dPWTime[i - 1]);
    dAlphafxCov[i] = dAlphafxCoVar;

    // Alpha Variance
    dAlphaVar = dAlpha[i] * dAlpha[i] * (dPWTime[i] - dPWTime[i - 1]);
    dAlphaVariance[i] = dAlphaVar;
  }

  //	Domestic Bond Variance
  domZCVar1 = (ddomSigma[i] / ddomLambda) * (ddomSigma[i] / ddomLambda) *
              (0.5 *
                   (exp(2 * ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                    exp(2 * ddomLambda * (ddomTStar - dFixTime))) /
                   ddomLambda -
               2.0 * exp(ddomLambda * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
                    exp(ddomLambda * (ddomTStar - dFixTime))) /
                   ddomLambda +
               exp(2 * ddomLambda * (ddomTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  domZCVar2 = (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
              (ddomAlphaLGM[i] * ddomSigma[i] / ddomLambda2) *
              (0.5 *
                   (exp(2 * ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                    exp(2 * ddomLambda2 * (ddomTStar - dFixTime))) /
                   ddomLambda2 -
               2.0 * exp(ddomLambda2 * (ddomTStar - dBondMat)) *
                   (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
                    exp(ddomLambda2 * (ddomTStar - dFixTime))) /
                   ddomLambda2 +
               exp(2 * ddomLambda2 * (ddomTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  domZCCoVar =
      ddomAlphaLGM[i] * ddomRhoLGM[i] * (ddomSigma[i] / ddomLambda) *
      (ddomSigma[i] / ddomLambda2) *
      ((exp((ddomLambda + ddomLambda2) * (ddomTStar - dPWTime[i - 1])) -
        exp((ddomLambda + ddomLambda2) * (ddomTStar - dFixTime))) /
           (ddomLambda + ddomLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda +
       exp((ddomLambda + ddomLambda2) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  domZCVar = domZCVar1 + domZCVar2 + 2 * domZCCoVar;
  ddomVar[i] = domZCVar;

  //	Foreign Bond Variance
  forZCVar1 = (dforSigma[i] / dforLambda) * (dforSigma[i] / dforLambda) *
              (0.5 *
                   (exp(2 * dforLambda * (dforTStar - dPWTime[i - 1])) -
                    exp(2 * dforLambda * (dforTStar - dFixTime))) /
                   dforLambda -
               2.0 * exp(dforLambda * (dforTStar - dBondMat)) *
                   (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
                    exp(dforLambda * (dforTStar - dFixTime))) /
                   dforLambda +
               exp(2 * dforLambda * (dforTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  forZCVar2 = (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
              (dforAlphaLGM[i] * dforSigma[i] / dforLambda2) *
              (0.5 *
                   (exp(2 * dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                    exp(2 * dforLambda2 * (dforTStar - dFixTime))) /
                   dforLambda2 -
               2.0 * exp(dforLambda2 * (dforTStar - dBondMat)) *
                   (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
                    exp(dforLambda2 * (dforTStar - dFixTime))) /
                   dforLambda2 +
               exp(2 * dforLambda2 * (dforTStar - dBondMat)) *
                   (dFixTime - dPWTime[i - 1]));

  forZCCoVar =
      dforAlphaLGM[i] * dforRhoLGM[i] * (dforSigma[i] / dforLambda) *
      (dforSigma[i] / dforLambda2) *
      ((exp((dforLambda + dforLambda2) * (dforTStar - dPWTime[i - 1])) -
        exp((dforLambda + dforLambda2) * (dforTStar - dFixTime))) /
           (dforLambda + dforLambda2) -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda +
       exp((dforLambda + dforLambda2) * (dforTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  forZCVar = forZCVar1 + forZCVar2 + 2 * forZCCoVar;
  dforVar[i] = forZCVar;

  // FX Variance
  fxVar = dfxSigma[i] * dfxSigma[i] * (dFixTime - dPWTime[i - 1]);
  dfxVar[i] = fxVar;

  // FX-Dom Covariance
  dfxdomCoVar1 =
      -dCorrelation[i][0][4] * dfxSigma[i] * ddomSigma[i] / ddomLambda *
      ((exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
        exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda -
       exp(ddomLambda * (ddomTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxdomCoVar2 =
      -dCorrelation[i][1][4] * dfxSigma[i] * ddomAlphaLGM[i] * ddomSigma[i] /
      ddomLambda2 *
      ((exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
        exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxdomCoVar = dfxdomCoVar1 + dfxdomCoVar2;
  ddomfxCov[i] = dfxdomCoVar;

  // FX-For Covariance
  dfxforCoVar1 =
      -dCorrelation[i][2][4] * dfxSigma[i] * dforSigma[i] / dforLambda *
      ((exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
        exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxforCoVar2 =
      -dCorrelation[i][3][4] * dfxSigma[i] * dforAlphaLGM[i] * dforSigma[i] /
      dforLambda2 *
      ((exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
        exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dfxforCoVar = dfxforCoVar1 + dfxforCoVar2;
  dforfxCov[i] = dfxforCoVar;

  // Dom-For Covariance
  ddomforCoVar11 =
      dCorrelation[i][0][2] * ddomSigma[i] * dforSigma[i] /
      (dforLambda * ddomLambda) *
      (-(exp((ddomLambda + dforLambda) * (ddomTStar - dFixTime)) -
         exp((ddomLambda + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda + dforLambda) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda +
       exp((ddomLambda + dforLambda) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar12 =
      dCorrelation[i][0][3] * ddomSigma[i] * dforSigma[i] * dforAlphaLGM[i] /
      (dforLambda2 * ddomLambda) *
      (-(exp((ddomLambda + dforLambda2) * (ddomTStar - dFixTime)) -
         exp((ddomLambda + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda + dforLambda2) -
       exp(ddomLambda * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda +
       exp((ddomLambda + dforLambda2) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar21 =
      dCorrelation[i][1][2] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] /
      (dforLambda * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda) * (ddomTStar - dFixTime)) -
         exp((ddomLambda2 + dforLambda) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda2 + dforLambda) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar22 =
      dCorrelation[i][1][3] * ddomSigma[i] * ddomAlphaLGM[i] * dforSigma[i] *
      dforAlphaLGM[i] / (dforLambda2 * ddomLambda2) *
      (-(exp((ddomLambda2 + dforLambda2) * (ddomTStar - dFixTime)) -
         exp((ddomLambda2 + dforLambda2) * (ddomTStar - dPWTime[i - 1]))) /
           (ddomLambda2 + dforLambda2) -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) *
           (exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
            exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) *
           (exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
            exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 +
       exp((ddomLambda2 + dforLambda2) * (ddomTStar - dBondMat)) *
           (dFixTime - dPWTime[i - 1]));

  ddomforCoVar =
      ddomforCoVar11 + ddomforCoVar12 + ddomforCoVar21 + ddomforCoVar22;
  ddomforCov[i] = ddomforCoVar;

  // Alpha-Dom Covariance
  dAlphadomCoVar1 =
      -dCorrelation[i][0][5] * dAlpha[i] * ddomSigma[i] / ddomLambda *
      ((exp(ddomLambda * (ddomTStar - dPWTime[i - 1])) -
        exp(ddomLambda * (ddomTStar - dFixTime))) /
           ddomLambda -
       exp(ddomLambda * (ddomTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dAlphadomCoVar2 =
      -dCorrelation[i][1][5] * dAlpha[i] * ddomAlphaLGM[i] * ddomSigma[i] /
      ddomLambda2 *
      ((exp(ddomLambda2 * (ddomTStar - dPWTime[i - 1])) -
        exp(ddomLambda2 * (ddomTStar - dFixTime))) /
           ddomLambda2 -
       exp(ddomLambda2 * (ddomTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dAlphadomCoVar = dAlphadomCoVar1 + dAlphadomCoVar2;
  dAlphadomCov[i] = dAlphadomCoVar;

  // Alpha-for Covariance
  dAlphaforCoVar1 =
      -dCorrelation[i][0][5] * dAlpha[i] * dforSigma[i] / dforLambda *
      ((exp(dforLambda * (dforTStar - dPWTime[i - 1])) -
        exp(dforLambda * (dforTStar - dFixTime))) /
           dforLambda -
       exp(dforLambda * (dforTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dAlphaforCoVar2 =
      -dCorrelation[i][1][5] * dAlpha[i] * dforAlphaLGM[i] * dforSigma[i] /
      dforLambda2 *
      ((exp(dforLambda2 * (dforTStar - dPWTime[i - 1])) -
        exp(dforLambda2 * (dforTStar - dFixTime))) /
           dforLambda2 -
       exp(dforLambda2 * (dforTStar - dBondMat)) * (dFixTime - dPWTime[i - 1]));

  dAlphaforCoVar = dAlphaforCoVar1 + dAlphaforCoVar2;
  dAlphaforCov[i] = dAlphaforCoVar;

  // FX-Alpha Covariance
  dAlphafxCoVar = dfxSigma[i] * dAlpha[i] * dCorrelation[i][4][5] *
                  (dFixTime - dPWTime[i - 1]);
  dAlphafxCov[i] = dAlphafxCoVar;

  // Alpha Variance
  dAlphaVar = dAlpha[i] * dAlpha[i] * (dFixTime - dPWTime[i - 1]);
  dAlphaVariance[i] = dAlphaVar;
}

void Compute5F_Equiv_Vol_Correls(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    double dFixTime,

    double *ddomVar, double *dforVar, double *dfxVar, double *dAlphaVariance,
    double *ddomfxCov, double *dforfxCov, double *ddomforCov,
    double *dAlphafxCov, double *dAlphadomCov, double *dAlphaforCov,

    double *h1, double *h2, double *h3, double *h4, double ***dCorrelation) {
  int i, j, k, FixTimeIndex;

  FixTimeIndex = Get_Index(dFixTime, dPWTime, iNbPWTime);

  // FX Variance
  h1[0] = sqrt(dfxVar[0] / dPWTime[0]);

  //	Domestic Bond Variance
  h2[0] = sqrt(ddomVar[0] / dPWTime[0]);

  //	Foreign Bond Variance
  h3[0] = sqrt(dforVar[0] / dPWTime[0]);

  //	Alpha
  h4[0] = sqrt(dAlphaVariance[0] / dPWTime[0]);

  // Correlations
  dCorrelation[0][0][1] = ddomfxCov[0] / sqrt(ddomVar[0] * dfxVar[0]);
  dCorrelation[0][0][2] = dforfxCov[0] / sqrt(dforVar[0] * dfxVar[0]);
  dCorrelation[0][0][3] = dAlphafxCov[0] / sqrt(dAlphaVariance[0] * dfxVar[0]);

  dCorrelation[0][1][2] = ddomforCov[0] / sqrt(ddomVar[0] * dforVar[0]);
  dCorrelation[0][1][3] =
      dAlphadomCov[0] / sqrt(dAlphaVariance[0] * ddomVar[0]);

  dCorrelation[0][2][3] =
      dAlphaforCov[0] / sqrt(dAlphaVariance[0] * dforVar[0]);

  for (k = 0; k < 4; ++k) {
    dCorrelation[0][k][k] = 1.0;
    for (j = k + 1; j < 4; ++j) {
      dCorrelation[0][j][k] = dCorrelation[0][k][j];
    }
  }

  for (i = 1; i < FixTimeIndex; ++i) {
    // FX Variance
    h1[i] = sqrt(dfxVar[i] / (dPWTime[i] - dPWTime[i - 1]));

    //	Domestic Bond Variance
    h2[i] = sqrt(ddomVar[i] / (dPWTime[i] - dPWTime[i - 1]));

    //	Foreign Bond Variance
    h3[i] = sqrt(dforVar[i] / (dPWTime[i] - dPWTime[i - 1]));

    //	Alpha
    h4[i] = sqrt(dAlphaVariance[i] / (dPWTime[i] - dPWTime[i - 1]));

    // Correlations
    dCorrelation[i][0][1] = ddomfxCov[i] / sqrt(ddomVar[i] * dfxVar[i]);
    dCorrelation[i][0][2] = dforfxCov[i] / sqrt(dforVar[i] * dfxVar[i]);
    dCorrelation[i][0][3] =
        dAlphafxCov[i] / sqrt(dAlphaVariance[i] * dfxVar[i]);

    dCorrelation[i][1][2] = ddomforCov[i] / sqrt(ddomVar[i] * dforVar[i]);
    dCorrelation[i][1][3] =
        dAlphadomCov[i] / sqrt(dAlphaVariance[i] * ddomVar[i]);

    dCorrelation[i][2][3] =
        dAlphaforCov[i] / sqrt(dAlphaVariance[i] * dforVar[i]);

    for (k = 0; k < 4; ++k) {
      dCorrelation[i][k][k] = 1.0;
      for (j = k + 1; j < 4; ++j) {
        dCorrelation[i][j][k] = dCorrelation[i][k][j];
      }
    }
  }

  // FX Variance
  h1[i] = sqrt(dfxVar[i] / (dFixTime - dPWTime[i - 1]));

  //	Domestic Bond Variance
  h2[i] = sqrt(ddomVar[i] / (dFixTime - dPWTime[i - 1]));

  //	Foreign Bond Variance
  h3[i] = sqrt(dforVar[i] / (dFixTime - dPWTime[i - 1]));

  //	Alpha
  h4[i] = sqrt(dAlphaVariance[i] / (dFixTime - dPWTime[i - 1]));

  // Correlations
  dCorrelation[i][0][1] = ddomfxCov[i] / sqrt(ddomVar[i] * dfxVar[i]);
  dCorrelation[i][0][2] = dforfxCov[i] / sqrt(dforVar[i] * dfxVar[i]);
  dCorrelation[i][0][3] = dAlphafxCov[i] / sqrt(dAlphaVariance[i] * dfxVar[i]);

  dCorrelation[i][1][2] = ddomforCov[i] / sqrt(ddomVar[i] * dforVar[i]);
  dCorrelation[i][1][3] =
      dAlphadomCov[i] / sqrt(dAlphaVariance[i] * ddomVar[i]);

  dCorrelation[i][2][3] =
      dAlphaforCov[i] / sqrt(dAlphaVariance[i] * dforVar[i]);

  for (k = 0; k < 4; ++k) {
    dCorrelation[i][k][k] = 1.0;
    for (j = k + 1; j < 4; ++j) {
      dCorrelation[i][j][k] = dCorrelation[i][j][k];
    }
  }
}

void PerturbationPrecomputations(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    double *ddomSigma, double ddomLambda, double *ddomAlphaLGM,
    double *ddomRhoLGM, double ddomGammaLGM, double ddomTStar,

    double *dforSigma, double dforLambda, double *dforAlphaLGM,
    double *dforRhoLGM, double dforGammaLGM, double dforTStar,

    double *dfxSigma,

    double *dAlpha,

    double ***dCorrelation, // 0 : Dom1
                            // 1 : Dom2
                            // 2 : For1
                            // 3 : For2
                            // 4 : Fx
                            // 5 : ForVol

    double dBondMat, double dFixTime,

    double *h1, double *h2, double *h3, double *h4,
    double ***dhCorrel) // 0 : h1
                        // 1 : h2
                        // 2 : h3
                        // 3 : h4
{
  Err err = NULL;
  double *ddomVar = NULL;
  double *dforVar = NULL;
  double *dfxVar = NULL;
  double *dAlphaVariance = NULL;
  double *ddomfxCov = NULL;
  double *dforfxCov = NULL;
  double *ddomforCov = NULL;
  double *dAlphafxCov = NULL;
  double *dAlphadomCov = NULL;
  double *dAlphaforCov = NULL;

  ddomVar = (double *)calloc(iNbPWTime, sizeof(double));
  dforVar = (double *)calloc(iNbPWTime, sizeof(double));
  dfxVar = (double *)calloc(iNbPWTime, sizeof(double));
  dAlphaVariance = (double *)calloc(iNbPWTime, sizeof(double));
  ddomfxCov = (double *)calloc(iNbPWTime, sizeof(double));
  dforfxCov = (double *)calloc(iNbPWTime, sizeof(double));
  ddomforCov = (double *)calloc(iNbPWTime, sizeof(double));
  dAlphafxCov = (double *)calloc(iNbPWTime, sizeof(double));
  dAlphadomCov = (double *)calloc(iNbPWTime, sizeof(double));
  dAlphaforCov = (double *)calloc(iNbPWTime, sizeof(double));

  if ((!ddomVar) || (!dforVar) || (!dfxVar) || (!dAlphaVariance) ||
      (!ddomfxCov) || (!dforfxCov) || (!ddomforCov) || (!dAlphafxCov) ||
      (!dAlphadomCov) || (!dAlphaforCov)) {
    err = "Allocation failed in PerturbationPrecomputations";
    goto FREE_RETURN;
  }

  Compute5FCovariancesVect(
      iNbPWTime, dPWTime,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM, ddomTStar,

      dforSigma, dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM, dforTStar,

      dfxSigma,

      dAlpha,

      dCorrelation, // 0 : Dom1
                    // 1 : Dom2
                    // 2 : For1
                    // 3 : For2
                    // 4 : Fx
                    // 5 : Alpha

      dBondMat, dFixTime,

      ddomVar, dforVar, dfxVar, dAlphaVariance, ddomfxCov, dforfxCov,
      ddomforCov, dAlphafxCov, dAlphadomCov, dAlphaforCov);

  Compute5F_Equiv_Vol_Correls(iNbPWTime, dPWTime,

                              dFixTime,

                              ddomVar, dforVar, dfxVar, dAlphaVariance,
                              ddomfxCov, dforfxCov, ddomforCov, dAlphafxCov,
                              dAlphadomCov, dAlphaforCov,

                              h1, h2, h3, h4, dhCorrel);

FREE_RETURN:

  if (ddomVar) {
    free(ddomVar);
    ddomVar = NULL;
  }

  if (dforVar) {
    free(dforVar);
    dforVar = NULL;
  }

  if (dfxVar) {
    free(dfxVar);
    dfxVar = NULL;
  }

  if (dAlphaVariance) {
    free(dAlphaVariance);
    dAlphaVariance = NULL;
  }

  if (ddomfxCov) {
    free(ddomfxCov);
    ddomfxCov = NULL;
  }

  if (dforfxCov) {
    free(dforfxCov);
    dforfxCov = NULL;
  }

  if (dforfxCov) {
    free(dforfxCov);
    dforfxCov = NULL;
  }

  if (ddomfxCov) {
    free(ddomfxCov);
    ddomfxCov = NULL;
  }

  if (ddomforCov) {
    free(ddomforCov);
    ddomforCov = NULL;
  }

  if (dAlphafxCov) {
    free(dAlphafxCov);
    dAlphafxCov = NULL;
  }

  if (dAlphadomCov) {
    free(dAlphadomCov);
    dAlphadomCov = NULL;
  }

  if (dAlphaforCov) {
    free(dAlphaforCov);
    dAlphaforCov = NULL;
  }
}

void FXImpliedVol5F(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    double *ddomSigma, double ddomLambda, double *ddomAlphaLGM,
    double *ddomRhoLGM, double ddomGammaLGM, double ddomTStar,

    double *dforSigma, double dforLambda, double *dforAlphaLGM,
    double *dforRhoLGM, double dforGammaLGM, double dforTStar,

    double *dfxSigma,

    double ***dCorrelation, // 0 : Dom1
                            // 1 : Dom2
                            // 2 : For1
                            // 3 : For2
                            // 4 : Fx

    double dBondMat, double dFixTime,

    double *dVol) {
  double ddomVar;
  double dforVar;
  double dfxVar;
  double ddomfxCov;
  double dforfxCov;
  double ddomforCov;

  Compute5FCovariances(
      iNbPWTime, // Piece Wise Constant Term Structures
      dPWTime,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM, ddomTStar,

      dforSigma, dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM, dforTStar,

      dfxSigma,

      dCorrelation, // 0 : Dom1
                    // 1 : Dom2
                    // 2 : For1
                    // 3 : For2
                    // 4 : Fx

      dBondMat, dFixTime,

      &ddomVar, &dforVar, &dfxVar, &ddomfxCov, &dforfxCov, &ddomforCov);

  *dVol = sqrt((ddomVar + dforVar + dfxVar + 2 * ddomfxCov + 2 * dforfxCov +
                2 * ddomfxCov) /
               dFixTime);
}

Err FXLGMSV_FXImpliedVolApprox_FromTS(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    // Domestic Underlying
    double *ddomSigma, double *ddomAlpha, double *ddomLambdaEps,
    double *ddomRho, double *ddomRho2,

    double ddomLambda, double *ddomAlphaLGM, double *ddomRhoLGM,
    double ddomGammaLGM,

    double ddomTStar,

    // Foreign Underlying
    double *dforSigma, double *dforAlpha, double *dforLambdaEps,
    double *dforRho, double *dforRho2,

    double dforLambda, double *dforAlphaLGM, double *dforRhoLGM,
    double dforGammaLGM,

    double dforTStar,

    // FX underlying
    double *dfxSigma,

    // Correlation
    double ***d5FCorrelation, double ***dFullCorrelation,

    // Product description
    double dFixTime, double dBondMat,

    /* Numerical parameters */
    int iNbX, double dIntegParam,

    // Outputs
    double *FXImpliedVol) {
  Err err = NULL;
  double ddomBondATMVol, dforBondATMVol;
  double ddomVar, dforVar, dfxVar, ddomfxCov, dforfxCov, ddomforCov;
  double *dforLevelEpsNew = NULL;
  double *dforLambdaEpsNew = NULL;
  double *ddomLambdaEpsNew = NULL;

  ddomLambdaEpsNew = (double *)calloc(iNbPWTime, sizeof(double));
  if (!ddomLambdaEpsNew) {
    err = "Allocation failed in FXLGMSV_FXImpliedVolApprox_FromTS";
    goto FREE_RETURN;
  }

  FwdAdjustLambdaOnDom(iNbPWTime, dPWTime, ddomSigma, ddomAlpha, ddomLambdaEps,
                       ddomRho, ddomRho2,

                       ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM,

                       ddomTStar, dBondMat,

                       ddomLambdaEpsNew);

  err = LGMSVBondVolApprox(iNbPWTime, dPWTime,

                           ddomSigma, ddomAlpha, ddomLambdaEps,
                           ddomLambdaEpsNew, ddomRho, ddomRho2,

                           ddomLambda, ddomTStar,

                           ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM,

                           -1, dFixTime, dBondMat,

                           iNbX, dIntegParam,

                           &ddomBondATMVol);
  if (err) {
    goto FREE_RETURN;
  }

  dforLevelEpsNew = (double *)calloc(iNbPWTime, sizeof(double));
  dforLambdaEpsNew = (double *)calloc(iNbPWTime, sizeof(double));
  if ((!dforLevelEpsNew) || (!dforLambdaEpsNew)) {
    err = "Allocation failed in FXLGMSV_FXImpliedVolApprox_FromTS";
    goto FREE_RETURN;
  }

  FwdAdjustLevelAndLambdaOnFor(
      iNbPWTime, dPWTime, dforAlpha, dforLambdaEps, dforLambdaEps,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM,

      dFullCorrelation,

      dforTStar, dBondMat,

      dforLevelEpsNew, dforLambdaEpsNew);

  QtoAdjustLevelAndLambda(
      iNbPWTime, dPWTime,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM,

      dforSigma, dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM,
      dforLevelEpsNew, dforLambdaEpsNew, dforAlpha,

      dfxSigma,

      dFullCorrelation,

      dforTStar,

      dforLevelEpsNew, dforLambdaEpsNew);

  err = LGMSVBondVolApprox(iNbPWTime, dPWTime,

                           dforSigma, dforAlpha, dforLevelEpsNew,
                           dforLambdaEpsNew, dforRho, dforRho2,

                           dforLambda, dforTStar,

                           dforAlphaLGM, dforRhoLGM, dforGammaLGM,

                           1, dFixTime, dBondMat,

                           iNbX, dIntegParam,

                           &dforBondATMVol);
  if (err) {
    goto FREE_RETURN;
  }

  Compute5FCovariances(
      iNbPWTime, dPWTime,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM, ddomTStar,

      dforSigma, dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM, dforTStar,

      dfxSigma,

      d5FCorrelation,

      dBondMat, dFixTime,

      &ddomVar, &dforVar, &dfxVar, &ddomfxCov, &dforfxCov, &ddomforCov);

  *FXImpliedVol =
      sqrt((dfxVar + ddomBondATMVol * ddomBondATMVol * dFixTime +
            dforBondATMVol * dforBondATMVol * dFixTime -
            2 * (ddomfxCov / sqrt(ddomVar * dfxVar)) * sqrt(dfxVar / dFixTime) *
                ddomBondATMVol * dFixTime +
            2 * (dforfxCov / sqrt(dforVar * dfxVar)) * sqrt(dfxVar / dFixTime) *
                dforBondATMVol * dFixTime -
            2 * (ddomforCov / sqrt(ddomVar * dforVar)) * ddomBondATMVol *
                dforBondATMVol * dFixTime) /
           dFixTime);

FREE_RETURN:

  if (dforLevelEpsNew) {
    free(dforLevelEpsNew);
    dforLevelEpsNew = NULL;
  }

  if (dforLambdaEpsNew) {
    free(dforLambdaEpsNew);
    dforLambdaEpsNew = NULL;
  }

  if (ddomLambdaEpsNew) {
    free(ddomLambdaEpsNew);
    ddomLambdaEpsNew = NULL;
  }

  return err;
}

/*	Merge fx        , rates and corr model term structures */
Err FXLGMSV_merge_fx_ir_ir_model_ts(
    int target_nstp,
    /* Domestic TS */
    int dom_n, double *dom_time, double *dom_sigma, double *dom_alphaLGM,
    double *dom_rhoLGM, double *dom_alphaSV, double *dom_lamSV,
    double *dom_rhoSV, double *dom_rho2SV,

    /* Foreign TS */
    int for_n, double *for_time, double *for_sigma, double *for_alphaLGM,
    double *for_rhoLGM, double *for_alphaSV, double *for_lamSV,
    double *for_rhoSV, double *for_rho2SV,

    /* FX TS */
    int fx_n, double *fx_time, double *fx_sigma,

    /* Correlations TS */
    int rho_n, double *rho_time, double ***correlation,

    /* ------- */
    /* Outputs */
    /* ------- */

    /* Merge TS Times*/
    double **merge_time, int *merge_n,

    /* Domestic TS */
    double **merge_dom_sigma, double **merge_dom_alphaLGM,
    double **merge_dom_rhoLGM, double **merge_dom_alphaSV,
    double **merge_dom_lamSV, double **merge_dom_rhoSV,
    double **merge_dom_rho2SV,

    /* Foreign TS */
    double **merge_for_sigma, double **merge_for_alphaLGM,
    double **merge_for_rhoLGM, double **merge_for_alphaSV,
    double **merge_for_lamSV, double **merge_for_rhoSV,
    double **merge_for_rho2SV,

    /* FX TS */
    double **merge_fx_sigma,

    /* Correlations TS */
    double ****merge_correlation) {
  int i, j, k, index;
  Err err = NULL;

  *merge_time = NULL;

  *merge_dom_sigma = NULL;
  *merge_dom_alphaLGM = NULL;
  *merge_dom_rhoLGM = NULL;
  *merge_dom_alphaSV = NULL;
  *merge_dom_lamSV = NULL;
  *merge_dom_rhoSV = NULL;
  *merge_dom_rho2SV = NULL;

  *merge_for_sigma = NULL;
  *merge_for_alphaLGM = NULL;
  *merge_for_rhoLGM = NULL;
  *merge_for_alphaSV = NULL;
  *merge_for_lamSV = NULL;
  *merge_for_rhoSV = NULL;
  *merge_for_rho2SV = NULL;

  *merge_fx_sigma = NULL;

  *merge_correlation = NULL;

  *merge_time = (double *)calloc(dom_n, sizeof(double));
  if (!(*merge_time)) {
    err = "Memory allocation error in FXLGMSV_merge_fx_ir_ir_model_ts";
    goto FREE_RETURN;
  }

  memcpy(*merge_time, dom_time, dom_n * sizeof(double));
  *merge_n = dom_n;

  num_f_concat_vector(merge_n, merge_time, for_n, for_time);
  num_f_concat_vector(merge_n, merge_time, fx_n, fx_time);
  num_f_concat_vector(merge_n, merge_time, rho_n, rho_time);

  num_f_sort_vector(*merge_n, *merge_time);
  num_f_unique_vector(merge_n, *merge_time);

  num_f_fill_vector_newalgo(merge_n, merge_time, target_nstp);

  *merge_dom_sigma = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_alphaLGM = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_rhoLGM = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_alphaSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_lamSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_rhoSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_dom_rho2SV = (double *)calloc(*merge_n, sizeof(double));

  *merge_for_sigma = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_alphaLGM = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_rhoLGM = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_alphaSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_lamSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_rhoSV = (double *)calloc(*merge_n, sizeof(double));
  *merge_for_rho2SV = (double *)calloc(*merge_n, sizeof(double));

  *merge_fx_sigma = (double *)calloc(*merge_n, sizeof(double));

  *merge_correlation = f3tensor(0, *merge_n - 1, 0, 6, 0, 6);

  if (!(*merge_dom_sigma) || !(*merge_dom_alphaLGM) || !(*merge_dom_rhoLGM) ||
      !(*merge_dom_alphaSV) || !(*merge_dom_lamSV) || !(*merge_dom_rhoSV) ||
      !(*merge_dom_rho2SV) || !(*merge_for_sigma) || !(*merge_for_alphaLGM) ||
      !(*merge_for_rhoLGM) || !(*merge_for_alphaSV) || !(*merge_for_lamSV) ||
      !(*merge_for_rhoSV) || !(*merge_for_rho2SV) || !(*merge_fx_sigma) ||
      !(*merge_correlation)) {
    err = "Memory allocation error in FXLGMSV_merge_fx_ir_ir_model_ts";
    goto FREE_RETURN;
  }

  for (i = *merge_n - 1; i >= 0; i--) {
    /* Domestic TS */
    index = Get_Index((*merge_time)[i], dom_time, dom_n);
    (*merge_dom_sigma)[i] = dom_sigma[index];
    if (!dom_alphaLGM) {
      (*merge_dom_alphaLGM)[i] = 0.0;
    } else {
      (*merge_dom_alphaLGM)[i] = dom_alphaLGM[index];
    }
    if (!dom_rhoLGM) {
      (*merge_dom_rhoLGM)[i] = 0.0;
    } else {
      (*merge_dom_rhoLGM)[i] = dom_rhoLGM[index];
    }
    (*merge_dom_alphaSV)[i] = dom_alphaSV[index];
    (*merge_dom_lamSV)[i] = dom_lamSV[index];
    (*merge_dom_rhoSV)[i] = dom_rhoSV[index];
    if (!dom_rho2SV) {
      (*merge_dom_rho2SV)[i] = 0.0;
    } else {
      (*merge_dom_rho2SV)[i] = dom_rho2SV[index];
    }

    /* Foreign TS */
    index = Get_Index((*merge_time)[i], for_time, for_n);
    (*merge_for_sigma)[i] = for_sigma[index];
    if (!for_alphaLGM) {
      (*merge_for_alphaLGM)[i] = 0.0;
    } else {
      (*merge_for_alphaLGM)[i] = for_alphaLGM[index];
    }
    if (!for_rhoLGM) {
      (*merge_for_rhoLGM)[i] = 0.0;
    } else {
      (*merge_for_rhoLGM)[i] = for_rhoLGM[index];
    }
    (*merge_for_alphaSV)[i] = for_alphaSV[index];
    (*merge_for_lamSV)[i] = for_lamSV[index];
    (*merge_for_rhoSV)[i] = for_rhoSV[index];
    if (!for_rho2SV) {
      (*merge_for_rho2SV)[i] = 0.0;
    } else {
      (*merge_for_rho2SV)[i] = for_rho2SV[index];
    }

    /* FX TS */
    (*merge_fx_sigma)[i] = fx_sigma[Get_Index((*merge_time)[i], fx_time, fx_n)];

    /* Correlations TS */
    index = Get_Index((*merge_time)[i], rho_time, rho_n);
    for (j = 0; j < 7; j++) {
      for (k = 0; k < 7; k++) {
        (*merge_correlation)[i][j][k] = correlation[index][j][k];
      }
    }
  }

FREE_RETURN:

  if (err) {
    if (*merge_dom_sigma) {
      free(*merge_dom_sigma);
      *merge_dom_sigma = NULL;
    }

    if (*merge_dom_alphaLGM) {
      free(*merge_dom_alphaLGM);
      *merge_dom_alphaLGM = NULL;
    }

    if (*merge_dom_rhoLGM) {
      free(*merge_dom_rhoLGM);
      *merge_dom_rhoLGM = NULL;
    }

    if (*merge_dom_alphaSV) {
      free(*merge_dom_alphaSV);
      *merge_dom_alphaSV = NULL;
    }

    if (*merge_dom_lamSV) {
      free(*merge_dom_lamSV);
      *merge_dom_lamSV = NULL;
    }

    if (*merge_dom_rhoSV) {
      free(*merge_dom_rhoSV);
      *merge_dom_rhoSV = NULL;
    }

    if (*merge_dom_rho2SV) {
      free(*merge_dom_rho2SV);
      *merge_dom_rho2SV = NULL;
    }

    if (*merge_for_sigma) {
      free(*merge_for_sigma);
      *merge_for_sigma = NULL;
    }

    if (*merge_for_alphaLGM) {
      free(*merge_for_alphaLGM);
      *merge_for_alphaLGM = NULL;
    }

    if (*merge_for_rhoLGM) {
      free(*merge_for_rhoLGM);
      *merge_for_rhoLGM = NULL;
    }

    if (*merge_for_alphaSV) {
      free(*merge_for_alphaSV);
      *merge_for_alphaSV = NULL;
    }

    if (*merge_for_lamSV) {
      free(*merge_for_lamSV);
      *merge_for_lamSV = NULL;
    }

    if (*merge_for_rhoSV) {
      free(*merge_for_rhoSV);
      *merge_for_rhoSV = NULL;
    }

    if (*merge_for_rho2SV) {
      free(*merge_for_rho2SV);
      *merge_for_rho2SV = NULL;
    }

    if (*merge_fx_sigma) {
      free(*merge_fx_sigma);
      *merge_fx_sigma = NULL;
    }

    if (*merge_correlation) {
      free_f3tensor(*merge_correlation, 0, *merge_n - 1, 0, 6, 0, 6);
      *merge_correlation = NULL;
    }
  }

  return err;
}

Err FXLGMSV_GetModelTS(
    // Underlying Name
    char *UndName,

    // TS times
    int *iNbPWTime, double **dPWTime,

    // Domestic Underlying
    double **ddomSigma, double **ddomAlpha, double **ddomLambdaEps,
    double **ddomRho, double **ddomRho2,

    double *ddomLambda, double **ddomAlphaLGM, double **ddomRhoLGM,
    double *ddomGammaLGM,

    double *ddomTStar,

    // Foreign Underlying
    double **dforSigma, double **dforAlpha, double **dforLambdaEps,
    double **dforRho, double **dforRho2,

    double *dforLambda, double **dforAlphaLGM, double **dforRhoLGM,
    double *dforGammaLGM,

    double *dforTStar,

    // FX underlying
    double **dfxSigma,

    // Correlation
    double ****d5FCorrelation, double ****dFullCorrelation) {
  Err err = NULL;
  int i, j, k, j2, k2;
  SrtUndPtr fx_und, dom_und, for_und;
  SrtUndPtr *und_ptr = NULL, und = NULL;
  LGMSV_model *dom_model = NULL;
  LGMSV_model *for_model = NULL;
  long today;

  char *domname, *forname;

  double fx_spot;
  int fx_n;
  double *fx_time = NULL;
  double *fx_vol = NULL;

  int rho_n;
  double *rho_time = NULL;
  double ***CorrMatrix = NULL;

  fx_und = lookup_und(UndName);

  if (!fx_und) {
    err = serror("Couldn't find underlying named %s", UndName);
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(fx_und);

  if (get_underlying_type(fx_und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", UndName);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(fx_und) != FX_LGMSV) {
    err = serror("Underlying %s is not of type FX Stoch Rates", UndName);
    goto FREE_RETURN;
  }

  fxlgmsv_get_fx_and_correl_ts(UndName, &fx_spot, &fx_n, &fx_time, &fx_vol,
                               &rho_n, &rho_time, &CorrMatrix);

  //	Initialisation of the LGMSV model
  dom_model = (LGMSV_model *)calloc(1, sizeof(LGMSV_model));
  for_model = (LGMSV_model *)calloc(1, sizeof(LGMSV_model));
  init_NULL_LGMSV_model(dom_model);
  init_NULL_LGMSV_model(for_model);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  if (!dom_und) {
    err = serror("Couldn't find underlying named %s", domname);
    goto FREE_RETURN;
  }
  // Get the model
  err = Get_LGMSV_model(domname, dom_model);

  forname = get_forname_from_fxund(fx_und);
  for_und = lookup_und(forname);
  if (!for_und) {
    err = serror("Couldn't find underlying named %s", forname);
    goto FREE_RETURN;
  }
  // Get the model
  err = Get_LGMSV_model(forname, for_model);

  *ddomLambda = dom_model->dLambdaX;
  *ddomGammaLGM = dom_model->dLGMGamma;
  *ddomTStar = dom_model->dTStar;

  *dforLambda = for_model->dLambdaX;
  *dforGammaLGM = for_model->dLGMGamma;
  *dforTStar = for_model->dTStar;

  err = FXLGMSV_merge_fx_ir_ir_model_ts(
      500,
      // Domestic TS
      dom_model->iNbPWTime, dom_model->dPWTime, dom_model->dSigma,
      dom_model->dLGMAlpha, dom_model->dLGMRho, dom_model->dAlpha,
      dom_model->dLambdaEps, dom_model->dRho, dom_model->dRho2,

      // Foreign TS
      for_model->iNbPWTime, for_model->dPWTime, for_model->dSigma,
      for_model->dLGMAlpha, for_model->dLGMRho, for_model->dAlpha,
      for_model->dLambdaEps, for_model->dRho, for_model->dRho2,

      // FX TS
      fx_n, fx_time, fx_vol,

      // Correlations TS
      rho_n, rho_time, CorrMatrix,

      // -------
      // Outputs
      // -------

      // Merge TS Times
      dPWTime, iNbPWTime,

      // Domestic TS
      ddomSigma, ddomAlphaLGM, ddomRhoLGM, ddomAlpha, ddomLambdaEps, ddomRho,
      ddomRho2,

      dforSigma, dforAlphaLGM, dforRhoLGM, dforAlpha, dforLambdaEps, dforRho,
      dforRho2,

      dfxSigma,

      dFullCorrelation);
  if (err) {
    goto FREE_RETURN;
  }

  *d5FCorrelation = f3tensor(0, *iNbPWTime - 1, 0, 4, 0, 4);
  if (!d5FCorrelation) {
    err = "allocation failed in ";
    goto FREE_RETURN;
  }

  for (i = 0; i < *iNbPWTime; ++i) {
    for (j = 0; j < 5; ++j) {
      j2 = j;
      if ((j == 2) || (j == 3)) {
        j2 = j + 1;
      } else if (j == 4) {
        j2 = 6;
      }
      for (k = 0; k < 5; ++k) {
        k2 = k;
        if ((k == 2) || (k == 3)) {
          k2 = k + 1;
        } else if (k == 4) {
          k2 = 6;
        }

        (*d5FCorrelation)[i][j][k] = (*dFullCorrelation)[i][j2][k2];
      }
    }
  }

FREE_RETURN:

  if (dom_model) {
    free_LGMSV_model(dom_model);
    free(dom_model);
    dom_model = NULL;
  }

  if (for_model) {
    free_LGMSV_model(for_model);
    free(for_model);
    for_model = NULL;
  }

  if (fx_time) {
    free(fx_time);
    fx_time = NULL;
  }

  if (fx_vol) {
    free(fx_time);
    fx_time = NULL;
  }

  if (rho_time) {
    free(rho_time);
    rho_time = NULL;
  }

  if (CorrMatrix) {
    free_f3tensor(CorrMatrix, 0, rho_n - 1, 0, 6, 0, 6);
    CorrMatrix = NULL;
  }

  if (err) {
    if (*dFullCorrelation) {
      free_f3tensor(*dFullCorrelation, 0, *iNbPWTime - 1, 0, 6, 0, 6);
      dFullCorrelation = NULL;
    }

    if (*dPWTime) {
      free(*dPWTime);
      *dPWTime = NULL;
    }

    if (*ddomSigma) {
      free(*ddomSigma);
      *ddomSigma = NULL;
    }

    if (*ddomAlpha) {
      free(*ddomAlpha);
      *ddomAlpha = NULL;
    }

    if (*ddomLambdaEps) {
      free(*ddomLambdaEps);
      *ddomLambdaEps = NULL;
    }

    if (*ddomRho) {
      free(*ddomRho);
      *ddomRho = NULL;
    }

    if (*ddomRho2) {
      free(*ddomRho2);
      *ddomRho2 = NULL;
    }

    if (*ddomAlphaLGM) {
      free(*ddomAlphaLGM);
      *ddomAlphaLGM = NULL;
    }

    if (*ddomRhoLGM) {
      free(*ddomRhoLGM);
      *ddomRhoLGM = NULL;
    }

    if (*dforSigma) {
      free(*dforSigma);
      *dforSigma = NULL;
    }

    if (*dforAlpha) {
      free(*dforAlpha);
      *dforAlpha = NULL;
    }

    if (*dforLambdaEps) {
      free(*dforLambdaEps);
      *dforLambdaEps = NULL;
    }

    if (*dforRho) {
      free(*dforRho);
      *dforRho = NULL;
    }

    if (*dforRho2) {
      free(*dforRho2);
      *dforRho2 = NULL;
    }

    if (*dforAlphaLGM) {
      free(*dforAlphaLGM);
      *dforAlphaLGM = NULL;
    }

    if (*dforRhoLGM) {
      free(*dforRhoLGM);
      *dforRhoLGM = NULL;
    }

    if (*dfxSigma) {
      free(*dfxSigma);
      *dfxSigma = NULL;
    }

    if (*d5FCorrelation) {
      free_f3tensor(*d5FCorrelation, 0, *iNbPWTime - 1, 0, 4, 0, 4);
      *d5FCorrelation = NULL;
    }
  }

  return err;
}

Err FXLGMSV_FXImpliedVolApprox(
    // Underlying Name
    char *UndName,

    // Product description
    double dFixTime, double dBondMat,

    // Numerical parameters
    int iNbX, double dIntegParam,

    // Outputs
    double *FXImpliedVol) {
  Err err = NULL;
  // TS times
  int iNbPWTime;
  double *dPWTime = NULL;

  // Domestic Underlying
  double *ddomSigma = NULL;
  double *ddomAlpha = NULL;
  double *ddomLambdaEps = NULL;
  double *ddomRho = NULL;
  double *ddomRho2 = NULL;

  double ddomLambda;
  double *ddomAlphaLGM = NULL;
  double *ddomRhoLGM = NULL;
  double ddomGammaLGM;

  double ddomTStar;

  // Foreign Underlying
  double *dforSigma = NULL;
  double *dforAlpha = NULL;
  double *dforLambdaEps = NULL;
  double *dforRho = NULL;
  double *dforRho2 = NULL;

  double dforLambda;
  double *dforAlphaLGM = NULL;
  double *dforRhoLGM = NULL;
  double dforGammaLGM;

  double dforTStar;

  // FX underlying
  double *dfxSigma = NULL;

  // Correlation
  double ***d5FCorrelation = NULL;
  double ***dFullCorrelation = NULL;

  err = FXLGMSV_GetModelTS(
      // Underlying Name
      UndName,

      // TS times
      &iNbPWTime, &dPWTime,

      // Domestic Underlying
      &ddomSigma, &ddomAlpha, &ddomLambdaEps, &ddomRho, &ddomRho2,

      &ddomLambda, &ddomAlphaLGM, &ddomRhoLGM, &ddomGammaLGM,

      &ddomTStar,

      // Foreign Underlying
      &dforSigma, &dforAlpha, &dforLambdaEps, &dforRho, &dforRho2,

      &dforLambda, &dforAlphaLGM, &dforRhoLGM, &dforGammaLGM,

      &dforTStar,

      // FX underlying
      &dfxSigma,

      // Correlation
      &d5FCorrelation, &dFullCorrelation);
  if (err) {
    goto FREE_RETURN;
  }

  err = FXLGMSV_FXImpliedVolApprox_FromTS(
      iNbPWTime, // Piece Wise Constant Term Structures
      dPWTime,

      // Domestic Underlying
      ddomSigma, ddomAlpha, ddomLambdaEps, ddomRho, ddomRho2,

      ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM,

      ddomTStar,

      // Foreign Underlying
      dforSigma, dforAlpha, dforLambdaEps, dforRho, dforRho2,

      dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM,

      dforTStar,

      // FX underlying
      dfxSigma,

      // Correlation
      d5FCorrelation, dFullCorrelation,

      // Product description
      dFixTime, dBondMat,

      /* Numerical parameters */
      iNbX, dIntegParam,

      // Outputs
      FXImpliedVol);
  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (dPWTime) {
    free(dPWTime);
    dPWTime = NULL;
  }

  if (ddomSigma) {
    free(ddomSigma);
    ddomSigma = NULL;
  }

  if (ddomAlpha) {
    free(ddomAlpha);
    ddomAlpha = NULL;
  }

  if (ddomLambdaEps) {
    free(ddomLambdaEps);
    ddomLambdaEps = NULL;
  }

  if (ddomRho) {
    free(ddomRho);
    ddomRho = NULL;
  }

  if (ddomRho2) {
    free(ddomRho2);
    ddomRho2 = NULL;
  }

  if (ddomAlphaLGM) {
    free(ddomAlphaLGM);
    ddomAlphaLGM = NULL;
  }

  if (ddomRhoLGM) {
    free(ddomRhoLGM);
    ddomRhoLGM = NULL;
  }

  if (dforSigma) {
    free(dforSigma);
    dforSigma = NULL;
  }

  if (dforAlpha) {
    free(dforAlpha);
    dforAlpha = NULL;
  }

  if (dforLambdaEps) {
    free(dforLambdaEps);
    dforLambdaEps = NULL;
  }

  if (dforRho) {
    free(dforRho);
    dforRho = NULL;
  }

  if (dforRho2) {
    free(dforRho2);
    dforRho2 = NULL;
  }

  if (dforAlphaLGM) {
    free(dforAlphaLGM);
    dforAlphaLGM = NULL;
  }

  if (dforRhoLGM) {
    free(dforRhoLGM);
    dforRhoLGM = NULL;
  }

  if (dfxSigma) {
    free(dfxSigma);
    dfxSigma = NULL;
  }

  if (d5FCorrelation) {
    free_f3tensor(d5FCorrelation, 0, iNbPWTime - 1, 0, 4, 0, 4);
    d5FCorrelation = NULL;
  }

  if (dFullCorrelation) {
    free_f3tensor(dFullCorrelation, 0, iNbPWTime - 1, 0, 6, 0, 6);
    dFullCorrelation = NULL;
  }

  return err;
}

void compute_vol_derivatives( // TS times
    int iNbPWTime, double *dPWTime,

    // Model parameters
    double *h1, double *h2, double *h3, double *h4, double kappa,
    double ***correlation,

    // Product parameters
    double Mat1, double Mat2,

    // Outputs
    double *dvol, double *dd1vol, double *dd2vol) {
  int i, indexMat1, indexMat2;
  double var, d1vol, d2vol;
  double d1Integral, d2Integral;

  //	var = (*dvol) * (*dvol) * Mat1;
  //	d1vol = *dd1vol;
  //	d2vol = *dd2vol;

  var = 0;
  d1vol = 0;
  d2vol = 0;

  d1Integral = 0;
  d2Integral = 0;

  indexMat1 = Get_Index(Mat1, dPWTime, iNbPWTime);
  indexMat2 = Get_Index(Mat2, dPWTime, iNbPWTime);
  i = indexMat1;

  var += (h1[i] * h1[i] + h2[i] * h2[i] + h3[i] * h3[i]

          - 2 * h1[i] * h2[i] * correlation[i][0][1] +
          2 * h1[i] * h3[i] * correlation[i][0][2]

          - 2 * h2[i] * h3[i] * correlation[i][1][2]) *
         (dPWTime[indexMat1] - Mat1);

  d1Integral += (h1[i] * h3[i] * correlation[i][0][2] -
                 h2[i] * h3[i] * correlation[i][1][2]) *
                (exp(-kappa * Mat1) - exp(-kappa * dPWTime[indexMat1])) / kappa;
  ;
  d1Integral += 2 * h3[i] * h3[i] *
                (exp(-kappa * Mat1) - exp(-kappa * dPWTime[indexMat1])) / kappa;

  d2Integral +=
      (h1[i] * h3[i] * correlation[i][0][2] -
       h2[i] * h3[i] * correlation[i][1][2]) *
      (exp(-2 * kappa * Mat1) - exp(-2 * kappa * dPWTime[indexMat1])) /
      (2 * kappa);

  for (i = indexMat1 + 1; i < indexMat2; ++i) {
    var += (h1[i] * h1[i] + h2[i] * h2[i] + h3[i] * h3[i]

            - 2 * h1[i] * h2[i] * correlation[i][0][1] +
            2 * h1[i] * h3[i] * correlation[i][0][2]

            - 2 * h2[i] * h3[i] * correlation[i][1][2]) *
           (dPWTime[i] - dPWTime[i - 1]);

    d1Integral += (h1[i] * h3[i] * correlation[i][0][2] -
                   h2[i] * h3[i] * correlation[i][1][2]) *
                  (exp(-kappa * dPWTime[i - 1]) - exp(-kappa * dPWTime[i])) /
                  kappa;
    ;
    d1Integral += 2 * h3[i] * h3[i] *
                  (exp(-kappa * dPWTime[i - 1]) - exp(-kappa * dPWTime[i])) /
                  kappa;

    d2Integral +=
        (h1[i] * h3[i] * correlation[i][0][2] -
         h2[i] * h3[i] * correlation[i][1][2]) *
        (exp(-2 * kappa * dPWTime[i - 1]) - exp(-2 * kappa * dPWTime[i])) /
        (2 * kappa);
  }

  i = indexMat2;
  var += (h1[i] * h1[i] + h2[i] * h2[i] + h3[i] * h3[i]

          - 2 * h1[i] * h2[i] * correlation[i][0][1] +
          2 * h1[i] * h3[i] * correlation[i][0][2]

          - 2 * h2[i] * h3[i] * correlation[i][1][2]) *
         (Mat2 - dPWTime[max(0, i - 1)]);

  d1Integral += (h1[i] * h3[i] * correlation[i][0][2] -
                 h2[i] * h3[i] * correlation[i][1][2]) *
                (exp(-kappa * dPWTime[max(0, i - 1)]) - exp(-kappa * Mat2)) /
                kappa;
  ;
  d1Integral += 2 * h3[i] * h3[i] *
                (exp(-kappa * dPWTime[max(0, i - 1)]) - exp(-kappa * Mat2)) /
                kappa;

  d2Integral +=
      (h1[i] * h3[i] * correlation[i][0][2] -
       h2[i] * h3[i] * correlation[i][1][2]) *
      (exp(-2 * kappa * dPWTime[max(0, i - 1)]) - exp(-2 * kappa * Mat2)) /
      (2 * kappa);

  *dvol = sqrt(var / (Mat2 - Mat1));
  *dd1vol = 0.25 * d1Integral / ((*dvol) * (Mat2 - Mat1));
  *dd2vol =
      -((*dd1vol) * (*dd1vol) + 0.125 * d2Integral / (Mat2 - Mat1)) / (*dvol);
}

/* The model :
        dS(t)/S(t) = ( h1(t) + h2(t) + sqrt(V(t)) ( h3(t) + h4(t) ) ) * dW(t)
        dV(t) = Kappa (Vinfinity - V(t)) + sqrt(V(t)) h5(t) * dW(t)

        With constant Kappa and Vinfinity = V(0) */

Err Heston_Pert_Approx(
    // TS times
    int iNbPWTime, double *dPWTime,

    // Model parameters
    double *h1, double *h2, double *h3, double *h4, double kappa,
    double ***correlation,

    // Product parameters
    double Maturity,

    // Numerical parameters
    int NHermiteQuad, int NLegendreQuad,

    // Output
    double *dVol) {
  Err err = NULL;
  int i, j, index;
  double *x_her = NULL;
  double *w_her = NULL;
  double *x_leg = NULL;
  double *w_leg = NULL;
  double vol, dvol_full, dvol, dvol2, dd1vol, dd2vol;
  double Mat1, Mat2;
  double price, Fwd, BSVega, BSVolga, BSVanna, coeffVanna, coeffVolga, Volga,
      Vanna;

  x_leg = (double *)calloc(NLegendreQuad + 1, sizeof(double));
  w_leg = (double *)calloc(NLegendreQuad + 1, sizeof(double));

  x_her = (double *)calloc(NHermiteQuad + 1, sizeof(double));
  w_her = (double *)calloc(NHermiteQuad + 1, sizeof(double));

  if ((!x_her) || (!w_her) || (!x_leg) || (!w_leg)) {
    err = "Allocation failed";
    goto FREE_RETURN;
  }

  HermiteStandard(x_her, w_her, NHermiteQuad);
  gauleg(0, Maturity, x_leg, w_leg, NLegendreQuad);

  dvol = 0.0;
  dd1vol = 0.0;
  dd2vol = 0.0;
  compute_vol_derivatives(iNbPWTime, dPWTime, h1, h2, h3, h4, kappa,
                          correlation, 0, Maturity, &dvol, &dd1vol, &dd2vol);
  price = srt_f_optblksch(1.0, 1.0, dvol, Maturity, 1.0, SRT_CALL, PREMIUM);
  dvol_full = dvol;
  Mat1 = 0.0;
  dvol = 0.0;
  dd1vol = 0.0;
  dd2vol = 0.0;
  for (i = 1; i <= NLegendreQuad; ++i) {
    Mat2 = x_leg[i];
    dvol = 0.0;
    dd1vol = 0.0;
    dd2vol = 0.0;
    compute_vol_derivatives(iNbPWTime, dPWTime, h1, h2, h3, h4, kappa,
                            correlation, 0.0, Mat2, &dvol, &dd1vol, &dd2vol);

    dvol2 = 0.0;
    dd1vol = 0.0;
    dd2vol = 0.0;
    compute_vol_derivatives(iNbPWTime, dPWTime, h1, h2, h3, h4, kappa,
                            correlation, Mat2, Maturity, &dvol2, &dd1vol,
                            &dd2vol);

    index = Get_Index(Mat2, dPWTime, iNbPWTime);

    for (j = 1; j <= NHermiteQuad; ++j) {
      Fwd = exp(-0.5 * dvol * dvol * Mat2 + sqrt(Mat2) * dvol * x_her[j]);
      vol = sqrt((dvol_full * dvol_full * Maturity - dvol * dvol * Mat2) /
                 (Maturity - Mat2));
      BSVega =
          srt_f_optblksch(Fwd, 1.0, vol, Maturity - Mat2, 1.0, SRT_CALL, VEGA);
      BSVolga =
          srt_f_optblksch(Fwd, 1.0, vol, Maturity - Mat2, 1.0, SRT_CALL, VOLGA);
      BSVanna =
          srt_f_optblksch(Fwd, 1.0, vol, Maturity - Mat2, 1.0, SRT_CALL, VANNA);
      Volga = dd2vol * BSVega + dd1vol * dd1vol * BSVolga;
      Vanna = dd1vol * BSVanna;
      coeffVanna = h1[index] * h4[index] * correlation[index][0][3] -
                   h2[index] * h4[index] * correlation[index][1][3] +
                   h3[index] * h4[index] * correlation[index][2][3];
      coeffVolga = 0.5 * h4[index] * h4[index];
      price += w_leg[i] * w_her[j] * (coeffVolga * Volga + coeffVanna * Vanna);
    }
    Mat1 = Mat2;
  }

  err =
      srt_f_optimpvol(price, 1, 1, Maturity, 1, SRT_CALL, SRT_LOGNORMAL, dVol);

FREE_RETURN:

  if (x_leg) {
    free(x_leg);
    x_leg = NULL;
  }

  if (w_leg) {
    free(w_leg);
    w_leg = NULL;
  }

  if (x_her) {
    free(x_her);
    x_her = NULL;
  }

  if (w_her) {
    free(w_her);
    w_her = NULL;
  }

  return err;
}

Err FXLGMSV_FXImpliedVol_Pert_From_TS(
    // Parameters of diffusion
    int iNbPWTime, // Piece Wise Constant Term Structures
    double *dPWTime,

    double *ddomSigma, double ddomLambda, double *ddomAlphaLGM,
    double *ddomRhoLGM, double ddomGammaLGM, double ddomTStar,

    double *dforSigma, double dforLambda, double *dforAlphaLGM,
    double *dforRhoLGM, double dforGammaLGM, double dforTStar,

    double *dforAlpha, double *dforLambdaEps,

    double *dfxSigma,

    double ***dCorrelation, // 0 : Dom1
                            // 1 : Dom2
                            // 2 : For1
                            // 3 : For2
                            // 4 : Fx
                            // 5 : ForVol

    double dBondMat, double dFixTime,

    // Numerical parameters
    int NHermiteQuad, int NLegendreQuad,

    double *dVol) {
  Err err = NULL;
  double *h1 = NULL;
  double *h2 = NULL;
  double *h3 = NULL;
  double *h4 = NULL;
  double ***dhCorrel = NULL;

  h1 = (double *)calloc(iNbPWTime, sizeof(double));
  h2 = (double *)calloc(iNbPWTime, sizeof(double));
  h3 = (double *)calloc(iNbPWTime, sizeof(double));
  h4 = (double *)calloc(iNbPWTime, sizeof(double));
  dhCorrel = f3tensor(0, iNbPWTime - 1, 0, 3, 0, 3);
  if ((!h1) || (!h2) || (!h3) || (!h4) || (!dhCorrel)) {
    err = "Allocation failed in FXLGMSV_FXImpliedVol_Pert_From_TS";
    goto FREE_RETURN;
  }

  PerturbationPrecomputations(
      iNbPWTime, dPWTime,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM, ddomTStar,

      dforSigma, dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM, dforTStar,

      dfxSigma,

      dforAlpha,

      dCorrelation,

      dBondMat, dFixTime,

      h1, h2, h3, h4, dhCorrel);

  err = Heston_Pert_Approx(iNbPWTime, dPWTime,

                           h1, h2, h3, h4, dforLambdaEps[0], dhCorrel,

                           dFixTime,

                           NHermiteQuad, NLegendreQuad,

                           dVol);

FREE_RETURN:

  if (h1) {
    free(h1);
    h1 = NULL;
  }

  if (h2) {
    free(h2);
    h2 = NULL;
  }

  if (h3) {
    free(h3);
    h3 = NULL;
  }

  if (h4) {
    free(h4);
    h4 = NULL;
  }

  if (dhCorrel) {
    free_f3tensor(dhCorrel, 0, iNbPWTime - 1, 0, 4, 0, 4);
    dhCorrel = NULL;
  }

  return err;
}

Err FXLGMSV_FXImpliedVolPert(
    // Underlying Name
    char *UndName,

    // Product description
    double dFixTime, double dBondMat,

    // Numerical parameters
    int NHermiteQuad, int NLegendreQuad,

    // Outputs
    double *FXImpliedVol) {
  Err err = NULL;
  // TS times
  int i, j, j2, k, iNbPWTime;
  double *dPWTime = NULL;

  // Domestic Underlying
  double *ddomSigma = NULL;
  double *ddomAlpha = NULL;
  double *ddomLambdaEps = NULL;
  double *ddomRho = NULL;
  double *ddomRho2 = NULL;

  double ddomLambda;
  double *ddomAlphaLGM = NULL;
  double *ddomRhoLGM = NULL;
  double ddomGammaLGM;

  double ddomTStar;

  // Foreign Underlying
  double *dforSigma = NULL;
  double *dforAlpha = NULL;
  double *dforLambdaEps = NULL;
  double *dforRho = NULL;
  double *dforRho2 = NULL;

  double dforLambda;
  double *dforAlphaLGM = NULL;
  double *dforRhoLGM = NULL;
  double dforGammaLGM;

  double dforTStar;

  // FX underlying
  double *dfxSigma = NULL;

  // Correlation
  double ***d5FCorrelation = NULL;
  double ***dFullCorrelation = NULL;
  double ***dCorrelation = NULL;

  err = FXLGMSV_GetModelTS(
      // Underlying Name
      UndName,

      // TS times
      &iNbPWTime, &dPWTime,

      // Domestic Underlying
      &ddomSigma, &ddomAlpha, &ddomLambdaEps, &ddomRho, &ddomRho2,

      &ddomLambda, &ddomAlphaLGM, &ddomRhoLGM, &ddomGammaLGM,

      &ddomTStar,

      // Foreign Underlying
      &dforSigma, &dforAlpha, &dforLambdaEps, &dforRho, &dforRho2,

      &dforLambda, &dforAlphaLGM, &dforRhoLGM, &dforGammaLGM,

      &dforTStar,

      // FX underlying
      &dfxSigma,

      // Correlation
      &d5FCorrelation, &dFullCorrelation);
  if (err) {
    goto FREE_RETURN;
  }

  dCorrelation = f3tensor(0, iNbPWTime - 1, 0, 5, 0, 5);
  for (i = 0; i < iNbPWTime; ++i) {
    for (j = 0; j < 5; ++j) {
      dCorrelation[i][j][j] = 1.0;
      for (k = j + 1; k < 5; ++k) {
        dCorrelation[i][j][k] = d5FCorrelation[i][j][k];
        dCorrelation[i][k][j] = dCorrelation[i][j][k];
      }
      dCorrelation[i][5][5] = 1.0;
      if ((j == 0) || (j == 1)) {
        j2 = j;
      } else if ((j == 2) || (j == 3)) {
        j2 = j + 1;
      } else {
        j2 = 6;
      }

      dCorrelation[i][5][j] = dFullCorrelation[i][5][j2];
      dCorrelation[i][j][5] = dCorrelation[i][5][j];
    }
  }

  err = FXLGMSV_FXImpliedVol_Pert_From_TS(
      iNbPWTime, dPWTime,

      ddomSigma, ddomLambda, ddomAlphaLGM, ddomRhoLGM, ddomGammaLGM, ddomTStar,

      dforSigma, dforLambda, dforAlphaLGM, dforRhoLGM, dforGammaLGM, dforTStar,

      dforAlpha, dforLambdaEps,

      dfxSigma,

      dCorrelation, // 0 : Dom1
                    // 1 : Dom2
                    // 2 : For1
                    // 3 : For2
                    // 4 : Fx
                    // 5 : ForVol
      dBondMat, dFixTime,

      NHermiteQuad, NLegendreQuad,

      FXImpliedVol);
  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (dPWTime) {
    free(dPWTime);
    dPWTime = NULL;
  }

  if (ddomSigma) {
    free(ddomSigma);
    ddomSigma = NULL;
  }

  if (ddomAlpha) {
    free(ddomAlpha);
    ddomAlpha = NULL;
  }

  if (ddomLambdaEps) {
    free(ddomLambdaEps);
    ddomLambdaEps = NULL;
  }

  if (ddomRho) {
    free(ddomRho);
    ddomRho = NULL;
  }

  if (ddomRho2) {
    free(ddomRho2);
    ddomRho2 = NULL;
  }

  if (ddomAlphaLGM) {
    free(ddomAlphaLGM);
    ddomAlphaLGM = NULL;
  }

  if (ddomRhoLGM) {
    free(ddomRhoLGM);
    ddomRhoLGM = NULL;
  }

  if (dforSigma) {
    free(dforSigma);
    dforSigma = NULL;
  }

  if (dforAlpha) {
    free(dforAlpha);
    dforAlpha = NULL;
  }

  if (dforLambdaEps) {
    free(dforLambdaEps);
    dforLambdaEps = NULL;
  }

  if (dforRho) {
    free(dforRho);
    dforRho = NULL;
  }

  if (dforRho2) {
    free(dforRho2);
    dforRho2 = NULL;
  }

  if (dforAlphaLGM) {
    free(dforAlphaLGM);
    dforAlphaLGM = NULL;
  }

  if (dforRhoLGM) {
    free(dforRhoLGM);
    dforRhoLGM = NULL;
  }

  if (dfxSigma) {
    free(dfxSigma);
    dfxSigma = NULL;
  }

  if (d5FCorrelation) {
    free_f3tensor(d5FCorrelation, 0, iNbPWTime - 1, 0, 4, 0, 4);
    d5FCorrelation = NULL;
  }

  if (dFullCorrelation) {
    free_f3tensor(dFullCorrelation, 0, iNbPWTime - 1, 0, 6, 0, 6);
    dFullCorrelation = NULL;
  }

  if (dCorrelation) {
    free_f3tensor(dCorrelation, 0, iNbPWTime - 1, 0, 5, 0, 5);
    dCorrelation = NULL;
  }

  return err;
}
