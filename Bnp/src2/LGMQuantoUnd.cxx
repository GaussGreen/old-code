
#include "LGMQuantoUnd.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define BIG_BIG_NUMBER 1E40

/*	Functions for the underlying */

/*	Fill underlying structure from a predefined underlying */
Err fill_lgmQto_und(char *und3dfx, char *dom_vc, char *for_vc, char *dom_ref,
                    char *dom_swap_freq, char *dom_swap_basis, char *for_ref,
                    char *for_swap_freq, char *for_swap_basis, LGMQTO_UND und) {
  int i;
  long today;
  Err err = NULL;
  SrtUndPtr fx_und, dom_und, for_und;
  TermStruct *fx_ts, *dom_ts, *for_ts;
  SrtCorrLstPtr sCorrlist;
  char *domname, *forname;
  double dom_lam, for_lam;
  char *dom_yc, *for_yc;

  double *domsigtime = NULL;
  double *domsig = NULL;
  long domsig_n;
  double *domtau = NULL;
  double *domtautime = NULL;
  long domtau_n;

  double *forsigtime = NULL;
  double *forsig = NULL;
  long forsig_n;
  double *fortau = NULL;
  double *fortautime = NULL;
  long fortau_n;

  double *fxsigtime = NULL;
  double *fxsig = NULL;
  long fxsig_n;

  double *merge_times = NULL;
  double *merge_dates = NULL;

  double *sigtime = NULL;
  double *sigdom = NULL;
  double *sigfor = NULL;
  double *sigfx = NULL;

  double *correltime = NULL;
  double *quantocorr = NULL;
  double *domforcorr = NULL;
  double *domfxcorr = NULL;
  long correl_n;

  int nb_merge_dates;

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->dom_sigma = NULL;
  und->for_sigma = NULL;
  und->fx_sigma = NULL;
  und->domfor_rho = NULL;
  und->quanto_rho = NULL;

  und->has_inst_data = 0;

  und->has_fwd_iv = 0;
  und->nb_fwdiv = 0;
  und->exercise_date = NULL;
  und->market_fwdiv = NULL;
  und->model_fwdiv = NULL;
  und->extra_fees = NULL;

  //	Now        , lookup underlyings involved and their term structures

  fx_und = lookup_und(und3dfx);

  if (!fx_und) {
    serror("Couldn't find underlying named %s", und3dfx);
    err = "Couldn't find FX3D underlying";
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(fx_und);

  if (get_underlying_type(fx_und) != FOREX_UND) {
    serror("Underlying %s is not of type FX", und3dfx);
    err = "Underlying is not of type FX";
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES) {
    serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
    err = "Underlying is not of type FX Stoch Rates";
    goto FREE_RETURN;
  }

  fx_ts = get_ts_from_fxund(fx_und);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  if (!dom_und) {
    serror("Couldn't find underlying named %s", domname);
    err = "Couldn't find domestic underlying";
    goto FREE_RETURN;
  }
  dom_ts = get_ts_from_irund(dom_und);

  forname = get_forname_from_fxund(fx_und);
  for_und = lookup_und(forname);
  if (!for_und) {
    serror("Couldn't find underlying named %s", forname);
    err = "Couldn't find foreign underlying";
    goto FREE_RETURN;
  }
  for_ts = get_ts_from_irund(for_und);

  //	Fill the model parameters as required by doublelgm1fQuanto_adi2

  sCorrlist = srt_f_GetTheCorrelationList();
  if (!sCorrlist->head->element) {
    err = "correlation list improperly initialised";
    goto FREE_RETURN;
  }

  //	 Get all the term structures
  err = Get_FX_StochRate_TermStructures_corr(
      und3dfx, &domsigtime, &domsig, &domsig_n, &domtautime, &domtau, &domtau_n,
      &forsigtime, &forsig, &forsig_n, &fortautime, &fortau, &fortau_n,
      &fxsigtime, &fxsig, &fxsig_n, &correltime, &domforcorr, &domfxcorr,
      &quantocorr, &correl_n);

  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(domtau, domtau_n, &dom_lam);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(fortau, fortau_n, &for_lam);
  if (err) {
    goto FREE_RETURN;
  }

  //	now merge all the term structure
  merge_dates = (double *)calloc(domsig_n, sizeof(double));
  if (!merge_dates) {
    err = "Memory allocation error (4) in SrtGrfnQuantoPDE";
    goto FREE_RETURN;
  }
  memcpy(merge_dates, domsigtime, domsig_n * sizeof(double));
  nb_merge_dates = domsig_n;
  num_f_concat_vector(&nb_merge_dates, &merge_dates, forsig_n, forsigtime);
  num_f_concat_vector(&nb_merge_dates, &merge_dates, fxsig_n, fxsigtime);
  num_f_concat_vector(&nb_merge_dates, &merge_dates, correl_n, correltime);
  num_f_sort_vector(nb_merge_dates, merge_dates);
  num_f_unique_vector(&nb_merge_dates, merge_dates);

  //	Fill the new term structures

  und->sigma_n = nb_merge_dates;

  und->sigma_date = (double *)calloc(nb_merge_dates, sizeof(double));
  und->sigma_time = (double *)calloc(nb_merge_dates, sizeof(double));
  und->dom_sigma = (double *)calloc(nb_merge_dates, sizeof(double));
  und->for_sigma = (double *)calloc(nb_merge_dates, sizeof(double));
  und->fx_sigma = (double *)calloc(nb_merge_dates, sizeof(double));
  und->domfor_rho = (double *)calloc(nb_merge_dates, sizeof(double));
  und->quanto_rho = (double *)calloc(nb_merge_dates, sizeof(double));

  if (!und->sigma_date || !und->sigma_time || !und->dom_sigma ||
      !und->for_sigma || !und->fx_sigma || !und->domfor_rho ||
      !und->quanto_rho) {
    err = "Memory allocation error (4) in SrtGrfnQuantoPDE";
    goto FREE_RETURN;
  }

  for (i = nb_merge_dates - 1; i >= 0; i--) {
    und->sigma_date[i] = today + merge_dates[i] * YEARS_IN_DAY;
    und->sigma_time[i] = merge_dates[i];
    und->dom_sigma[i] = find_sig(merge_dates[i], dom_ts);
    und->for_sigma[i] = find_sig(merge_dates[i], for_ts);
    und->fx_sigma[i] = find_fx_sig(merge_dates[i], fx_ts);

    err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname,
                                       merge_dates[i], &(und->domfor_rho[i]));
    if (err) {
      goto FREE_RETURN;
    }
    err = srt_f_get_corr_from_CorrList(sCorrlist, forname, und3dfx,
                                       merge_dates[i], &(und->quanto_rho[i]));
    if (err) {
      goto FREE_RETURN;
    }
  }

  //	Get yield curves
  dom_yc = get_ycname_from_irund(dom_und);
  for_yc = get_ycname_from_irund(for_und);

  strcpy(und->name, und3dfx);
  und->today = today;

  strcpy(und->dom_yc, dom_yc);
  strcpy(und->for_yc, for_yc);

  strcpy(und->dom_vc, dom_vc);
  strcpy(und->for_vc, for_vc);

  strcpy(und->dom_ref, dom_ref);
  strcpy(und->dom_swap_freq, dom_swap_freq);
  strcpy(und->dom_swap_basis, dom_swap_basis);

  strcpy(und->for_ref, for_ref);
  strcpy(und->for_swap_freq, for_swap_freq);
  strcpy(und->for_swap_basis, for_swap_basis);

  und->dom_lambda = dom_lam;
  und->for_lambda = for_lam;

  und->has_inst_data = 0;
  //	cpd_init_calib_inst_data (&(und->dom_inst_data));
  //	cpd_init_calib_inst_data (&(und->for_inst_data));

FREE_RETURN:
  if (err) {
    if (und)
      free_lgmQto_und(und);
  }

  return err;
}

void copy_lgmQto_und(LGMQTO_UND src, LGMQTO_UND dest) {
  strcpy(dest->name, src->name);
  dest->today = src->today;

  strcpy(dest->dom_yc, src->dom_yc);
  strcpy(dest->dom_vc, src->dom_vc);
  strcpy(dest->dom_ref, src->dom_ref);
  strcpy(dest->dom_swap_freq, src->dom_swap_freq);
  strcpy(dest->dom_swap_basis, src->dom_swap_basis);

  strcpy(dest->for_yc, src->for_yc);
  strcpy(dest->for_vc, src->for_vc);
  strcpy(dest->for_ref, src->for_ref);
  strcpy(dest->for_swap_freq, src->for_swap_freq);
  strcpy(dest->for_swap_basis, src->for_swap_basis);

  dest->sigma_n = src->sigma_n;
  dest->dom_lambda = src->dom_lambda;
  dest->for_lambda = src->for_lambda;

  if (src->sigma_date) {
    dest->sigma_date = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma_date, src->sigma_date, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma_date = NULL;
  }

  if (src->sigma_time) {
    dest->sigma_time = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma_time, src->sigma_time, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma_time = NULL;
  }

  if (src->dom_sigma) {
    dest->dom_sigma = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->dom_sigma, src->dom_sigma, dest->sigma_n * sizeof(double));
  } else {
    dest->dom_sigma = NULL;
  }

  if (src->for_sigma) {
    dest->for_sigma = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->for_sigma, src->for_sigma, dest->sigma_n * sizeof(double));
  } else {
    dest->for_sigma = NULL;
  }

  if (src->fx_sigma) {
    dest->fx_sigma = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->fx_sigma, src->fx_sigma, dest->sigma_n * sizeof(double));
  } else {
    dest->fx_sigma = NULL;
  }

  if (src->domfor_rho) {
    dest->domfor_rho = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->domfor_rho, src->domfor_rho, dest->sigma_n * sizeof(double));
  } else {
    dest->domfor_rho = NULL;
  }

  if (src->quanto_rho) {
    dest->quanto_rho = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->quanto_rho, src->quanto_rho, dest->sigma_n * sizeof(double));
  } else {
    dest->quanto_rho = NULL;
  }

  if (src->has_inst_data) {
    cpd_copy_calib_inst_data(&(dest->dom_inst_data), &(src->dom_inst_data));
    cpd_copy_calib_inst_data(&(dest->for_inst_data), &(src->for_inst_data));
    dest->has_inst_data = 1;
  } else {
    dest->has_inst_data = 0;
  }

  if (src->has_fwd_iv) {
    dest->has_fwd_iv = 1;
    dest->nb_fwdiv = src->nb_fwdiv;

    dest->exercise_date = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->exercise_date, src->exercise_date,
           dest->nb_fwdiv * sizeof(double));

    dest->market_fwdiv = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->market_fwdiv, src->market_fwdiv,
           dest->nb_fwdiv * sizeof(double));

    dest->model_fwdiv = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->model_fwdiv, src->model_fwdiv,
           dest->nb_fwdiv * sizeof(double));

    dest->extra_fees = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->extra_fees, src->extra_fees, dest->nb_fwdiv * sizeof(double));
  } else {
    dest->has_fwd_iv = 0;
    dest->nb_fwdiv = 0;
    dest->exercise_date = NULL;
    dest->market_fwdiv = NULL;
    dest->model_fwdiv = NULL;
    dest->extra_fees = NULL;
  }
}

Err free_lgmQto_und(LGMQTO_UND und) {
  if (und->sigma_date)
    free(und->sigma_date);
  if (und->sigma_time)
    free(und->sigma_time);
  if (und->dom_sigma)
    free(und->dom_sigma);
  if (und->for_sigma)
    free(und->for_sigma);
  if (und->fx_sigma)
    free(und->fx_sigma);
  if (und->domfor_rho)
    free(und->domfor_rho);
  if (und->quanto_rho)
    free(und->quanto_rho);

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->dom_sigma = NULL;
  und->for_sigma = NULL;
  und->fx_sigma = NULL;
  und->domfor_rho = NULL;
  und->quanto_rho = NULL;

  if (und->has_inst_data) {
    cpd_free_calib_inst_data(&und->dom_inst_data);
    cpd_free_calib_inst_data(&und->for_inst_data);
    und->has_inst_data = 0;
  }

  if (und->has_fwd_iv) {
    if (und->market_fwdiv)
      free(und->market_fwdiv);
    if (und->model_fwdiv)
      free(und->model_fwdiv);
    if (und->extra_fees)
      free(und->extra_fees);
    if (und->exercise_date)
      free(und->exercise_date);

    und->has_fwd_iv = 0;
    und->nb_fwdiv = 0;
  }

  return NULL;
}

//----Warning : (void*) prm is not freed in this function-------------
Err lgmQto_free_adi_arg(LGMQTO_ADI_ARG adi_arg) {
  //	int i;

  if (adi_arg->time)
    free(adi_arg->time);
  if (adi_arg->date)
    free(adi_arg->date);
  if (adi_arg->dom_ifr)
    free(adi_arg->dom_ifr);
  if (adi_arg->for_ifr)
    free(adi_arg->for_ifr);
  if (adi_arg->is_event)
    free(adi_arg->is_event);

  //	for (i=0; i<adi_arg->nstp;++i)
  //	{
  //		if(adi_arg->void_prm[i]) free(adi_arg->void_prm[i]);
  //	}

  if (adi_arg->void_prm)
    free(adi_arg->void_prm);

  adi_arg->time = NULL;
  adi_arg->date = NULL;
  adi_arg->domsig = NULL;
  adi_arg->forsig = NULL;
  adi_arg->fxsig = NULL;
  adi_arg->domforrho = NULL;
  adi_arg->quantorho = NULL;
  adi_arg->sig_time = NULL;
  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;

  return NULL;
}

Err lgmQto_fill_adi_arg(
    LGMQTO_UND und, int nb_event_dates, double *event_time,
    void *prm,          //
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    int req_stp, int req_stpx, LGMQTO_ADI_ARG adi_arg) {
  Err err = NULL;
  int i, j;

  /*	Initialise */
  adi_arg->time = NULL;
  adi_arg->date = NULL;

  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->dom_ifr = NULL;
  adi_arg->for_ifr = NULL;

  /*	Compute time steps */

  /*	Copy event dates */

  adi_arg->nstp = nb_event_dates;
  adi_arg->nstpx = req_stpx;

  adi_arg->time = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->time) {
    err = "Memory allocation error (1) in lgmQto_fill_adi_arg";
    goto FREE_RETURN;
  }
  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->time[i] = event_time[i];
  }

  /*	Fill time vector */

  /*	Add today if required */
  if (adi_arg->time[0] < -EPS) {
    err = "Past event date in lgmQto_fill_adi_arg";
    goto FREE_RETURN;
  }
  if (adi_arg->time[0] > EPS) {
    num_f_add_number(&(adi_arg->nstp), &(adi_arg->time), 0.0);
    num_f_sort_vector(adi_arg->nstp, adi_arg->time);
    num_f_unique_vector(&(adi_arg->nstp), adi_arg->time);
  }

  /*	If only one event today        , add empty event */
  if (adi_arg->nstp == 1) {
    num_f_add_number(&(adi_arg->nstp), (&adi_arg->time), 1.0);
  }

  /*	Fill the vector */
  num_f_fill_vector_newalgo(&(adi_arg->nstp), &(adi_arg->time), req_stp);

  /*	Make dates */
  adi_arg->date = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->date) {
    err = "Memory allocation error (2) in lgmQto_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->date[i] = und->today + DAYS_IN_YEAR * adi_arg->time[i];

    if (i > 0 && adi_arg->date[i] - adi_arg->date[i - 1] >= 1) {
      adi_arg->date[i] = (long)(adi_arg->date[i] + 1.0e-08);
      adi_arg->time[i] = YEARS_IN_DAY * (adi_arg->date[i] - und->today);
    }
  }

  /*	Sigma */

  adi_arg->sig_time = und->sigma_time;
  adi_arg->domsig = und->dom_sigma;
  adi_arg->forsig = und->for_sigma;
  adi_arg->fxsig = und->fx_sigma;
  adi_arg->domforrho = und->domfor_rho;
  adi_arg->quantorho = und->quanto_rho;
  adi_arg->nb_sig = und->sigma_n;

  /*	Lambdas and ratio */

  adi_arg->domlambda = und->dom_lambda;
  adi_arg->forlambda = und->for_lambda;

  /*	Spot fx and yield curves */

  /*	Domestic yield curve */
  strcpy(adi_arg->dom_yc, und->dom_yc);

  /*	Fill distributions */

  adi_arg->dom_ifr = (double *)calloc(adi_arg->nstp, sizeof(double));

  if (!adi_arg->dom_ifr) {
    err = "Memory allocation error (3) in lgmQto_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp - 1; i++) {
    adi_arg->dom_ifr[i] =
        swp_f_zr(adi_arg->date[i], adi_arg->date[i + 1], adi_arg->dom_yc);
  }

  /*	Foreign yield curve */
  strcpy(adi_arg->for_yc, und->for_yc);

  /*	Fill distributions */

  adi_arg->for_ifr = (double *)calloc(adi_arg->nstp, sizeof(double));

  if (!adi_arg->for_ifr) {
    err = "Memory allocation error (3) in lgmQto_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp - 1; i++) {
    adi_arg->for_ifr[i] =
        swp_f_zr(adi_arg->date[i], adi_arg->date[i + 1], adi_arg->for_yc);
  }

  /*	Fill limit conditions (product) */

  adi_arg->is_event = (int *)calloc(adi_arg->nstp, sizeof(int));
  adi_arg->void_prm = (void **)calloc(adi_arg->nstp, sizeof(void *));

  if (!adi_arg->is_event || !adi_arg->void_prm) {
    err = "Memory allocation error (4) in lgmQto_fill_adi_arg";
    goto FREE_RETURN;
  }

  j = nb_event_dates - 1;

  for (i = adi_arg->nstp - 1; i >= 0; i--) {
    if (j >= 0 && fabs(adi_arg->time[i] - event_time[j]) < 1.0e-08) {
      adi_arg->is_event[i] = 1;
      adi_arg->void_prm[i] = (void *)prm;
      j--;
    } else {
      adi_arg->is_event[i] = 0;
      adi_arg->void_prm[i] = NULL;
    }
  }

FREE_RETURN:

  if (err) {
    lgmQto_free_adi_arg(adi_arg);
  }

  return err;
}
