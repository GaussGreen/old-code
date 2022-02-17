/* -------------------------------------------------------------------------

   FILE NAME	: srt_f_grfnfutclsdfrm.cxx

   PURPOSE		: functions to price swaptions        , caps or
   bondoptions without the Grfn interface as of a date in the FUTURE        ,
   taking todays forward curve as input. THese functions take simplified
   descriptions of the products to be priced        , and construct arrays of
   Grfn cells to feed to the Grfn function If LGM can be used        , then
   switch to closed form
   ------------------------------------------------------------------------- */

#include "SrtAccess.h"
#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_grfclsdfrm.h"

#define VERY_SMALL 1.0e-08

Err srt_f_grfn_future_closedform(long future_date, SrtUndPtr und,
                                 SrtGrfnParam *grfnparam, SwapDP *sdp,
                                 double strike, double bondstrike,
                                 SrtReceiverType rec_pay, StructType type,
                                 String ref_rate_code, double *price) {
  Err err = NULL;
  int status = 0;
  String yc_name;
  SrtCurvePtr yldcrv;
  double df;
  double discprice;
  SrtMdlDim mdl_dim;
  SrtMdlType mdl_type;
  TermStruct *initial_ts, *changed_ts;
  Date today;
  double *sigma[6], *tau[3];
  double *sigma_mod[6];
  long new_num_sig;
  long i, next_sigma;
  long num_sig;
  long num_tau;
  long insert_at;

  /* Gets the Number of Factors in the Model */
  err = get_underlying_mdldim(und, &mdl_dim);
  if (err) {
    return err;
  }

  /* Gets the Model type (LGM        , Cheyette        ,...) */
  err = get_underlying_mdltype(und, &mdl_type);
  if (err) {
    return err;
  }

  /* Gets the initial Term Structure from the Underlying */
  err = get_underlying_ts(und, &initial_ts);
  if (err) {
    return err;
  }

  /* Gets th Yield Curve through the name attached to the underlying */
  yc_name = get_discname_from_underlying(und);
  yldcrv = lookup_curve(yc_name);

  /* Gets today's date from the Yield Curve */
  today = get_clcndate_from_yldcrv(yldcrv);

  /* Extract the Term Struct according to the model dimension */
  if (mdl_dim == ONE_FAC) {
    err = srt_f_display_IRM_OneFac_TermStruct(
        initial_ts, &(sigma[0]), &(sigma[1]), &(sigma[2]), &(sigma[3]),
        &(sigma[4]), &(sigma[5]), &num_sig, &(tau[0]), &(tau[1]), &num_tau);
    if (err) {
      return err;
    }
  } else if (mdl_dim == TWO_FAC) {
    err = srt_f_display_IRM_TwoFac_TermStruct(
        initial_ts, &(sigma[0]), &(sigma[1]), &(sigma[2]), &(sigma[3]),
        &(sigma[4]), &(sigma[5]), &num_sig, &(tau[0]), &(tau[1]), &(tau[2]),
        &num_tau);
    if (err) {
      return err;
    }
  } else {
    return err;
  }

  /* Finds the index of the first date in the Sigma Term Struct that is on or
   * after the Future Date */
  next_sigma = 0;
  while ((sigma[0][next_sigma] < future_date) && (next_sigma < num_sig)) {
    next_sigma++;
  }
  if (next_sigma == num_sig) {
    insert_at = num_sig;
    new_num_sig = num_sig + 1;
    next_sigma = num_sig - 1;
  } else if (sigma[0][next_sigma] == future_date) {
    insert_at = next_sigma;
    new_num_sig = num_sig;
  } else {
    insert_at = next_sigma;
    new_num_sig = num_sig + 1;
  }

  /* Allocate memory for a new tableau of sigmas and taus */
  sigma_mod[0] = dvector(0, new_num_sig - 1);
  sigma_mod[1] = dvector(0, new_num_sig - 1);
  sigma_mod[2] = dvector(0, new_num_sig - 1);
  if (mdl_dim == TWO_FAC) {
    sigma_mod[3] = dvector(0, new_num_sig - 1);
    sigma_mod[4] = dvector(0, new_num_sig - 1);
    sigma_mod[5] = dvector(0, new_num_sig - 1);
  }

  /* Puts the new point correspnding to the future date (with the right vol...)
   */
  sigma_mod[0][insert_at] = future_date;
  sigma_mod[1][insert_at] = sigma[1][next_sigma];
  sigma_mod[2][insert_at] = sigma[2][next_sigma];
  if (mdl_dim == TWO_FAC) {
    sigma_mod[3][insert_at] = sigma[3][next_sigma];
    sigma_mod[4][insert_at] = sigma[4][next_sigma];
    sigma_mod[5][insert_at] = sigma[5][next_sigma];
  }

  /* Puts all the vols before this point to 0.0 */
  for (i = 0; i < insert_at; i++) {
    sigma_mod[0][i] = sigma[0][i];
    sigma_mod[1][i] = VERY_SMALL;
    if (mdl_dim == TWO_FAC) {
      sigma_mod[3][i] = VERY_SMALL;
    }
  }

  /* Puts all the vols after as the normal ones */
  for (i = insert_at + 1; i < new_num_sig; i++) {
    sigma_mod[0][i] = sigma[0][i - (new_num_sig - num_sig)];
    sigma_mod[1][i] = sigma[1][i - (new_num_sig - num_sig)];
    sigma_mod[2][i] = sigma[2][i - (new_num_sig - num_sig)];
    if (mdl_dim == TWO_FAC) {
      sigma_mod[3][i] = sigma[3][i - (new_num_sig - num_sig)];
      sigma_mod[4][i] = sigma[4][i - (new_num_sig - num_sig)];
      sigma_mod[5][i] = sigma[5][i - (new_num_sig - num_sig)];
    }
  }

  /* Frees the tableau of initial sigmas and taus */
  srt_free(sigma[0]);
  srt_free(sigma[1]);
  srt_free(sigma[2]);
  srt_free(sigma[3]);
  srt_free(sigma[4]);
  srt_free(sigma[5]);

  /* Builds the term structure depending on the number of factors */
  if (mdl_dim == ONE_FAC) {
    err = srt_f_init_IRM_TermStruct(today, sigma_mod, 2, new_num_sig, tau, 2,
                                    num_tau, mdl_type, mdl_dim, 0.0, /* BEta*/
                                    0.0,                             /* Alpha*/
                                    0.0,                             /* Gamma */
                                    0.0,                             /* Rho */
                                    0.0,                             /*Vovol */
                                    0.0,                             /* Omega */
                                    0.0, /* vasicek parms */
                                    0, 0, NULL, &changed_ts);
    if (err) {
      return err;
    }
  }
  if (mdl_dim == TWO_FAC) {
    err = srt_f_init_IRM_TermStruct(
        today, sigma_mod, 4, new_num_sig, tau, 3, num_tau, mdl_type, mdl_dim,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, NULL, &changed_ts);
    if (err) {
      return err;
    }
  }

  /* Attaches the new Term structure to the Underlying */
  set_irund_ts(und, changed_ts);

  /* Calls Grfn Closed Form to price (as of today) this option with no vol until
   * the future date*/
  err = srt_f_grfn_clsdfrm(und, grfnparam, sdp, strike, bondstrike, rec_pay,
                           type, ref_rate_code, &discprice);
  if (err) {
    return err;
  }

  /* Reset the Initial Term Structure into the Underlying...*/
  set_irund_ts(und, initial_ts);

  /* Frees all allocates memory */
  free_dvector(sigma_mod[0], 0, new_num_sig - 1);
  free_dvector(sigma_mod[1], 0, new_num_sig - 1);
  free_dvector(sigma_mod[2], 0, new_num_sig - 1);
  srt_free(tau[0]);
  srt_free(tau[1]);
  if (mdl_dim == TWO_FAC) {
    free_dvector(sigma_mod[3], 0, new_num_sig - 1);
    free_dvector(sigma_mod[4], 0, new_num_sig - 1);
    free_dvector(sigma_mod[5], 0, new_num_sig - 1);
    srt_free(tau[2]);
  }

  /* Compute the discount factor for future date */
  df = swp_f_df(today, future_date, yc_name);

  /* Return the future price = price / df */
  *price = discprice / df;

  return err;
}
