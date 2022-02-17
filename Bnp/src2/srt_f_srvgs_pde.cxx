#include "SRT_H_ALL.H>
#include "pde_h_struct.h"
#include "srt_h_FeynmannKac.h"
#include "srt_h_srvgs_pde_fct.h"

#include "OPFNCTNS.H>
#define MAX_NUM_ITER 30

static Err make_pde_operators_allocations(SrtBasicPdeInf pde_inf,
                                          double **payoff,
                                          double ****diffusion_mat,
                                          double ****convection_mat,
                                          double ****dividend_tensor,
                                          double **option_pv_vect) {
  (*diffusion_mat) = f3tensor(0, pde_inf.num_pde - 1, 0, pde_inf.num_step, 1,
                              pde_inf.num_mesh);
  (*convection_mat) = f3tensor(0, pde_inf.num_pde - 1, 0, pde_inf.num_step, 1,
                               pde_inf.num_mesh);
  (*dividend_tensor) = f3tensor(0, pde_inf.num_pde - 1, 0, pde_inf.num_step, 1,
                                pde_inf.num_mesh);

  (*option_pv_vect) = dvector(1, pde_inf.num_mesh - 1);
  (*payoff) = dvector(1, pde_inf.num_mesh - 1);

  return NULL;
}

static Err make_pde_operators_desallocations(SrtBasicPdeInf pde_inf,
                                             double **payoff,
                                             double ****diffusion_mat,
                                             double ****convection_mat,
                                             double ****dividend_tensor,
                                             double **option_pv_vect) {
  if ((*payoff))
    free_dvector((*payoff), 1, pde_inf.num_mesh - 1);
  (*payoff) = NULL;
  if ((*option_pv_vect))
    free_dvector((*option_pv_vect), 1, pde_inf.num_mesh - 1);
  (*option_pv_vect) = NULL;

  if ((*diffusion_mat))
    free_f3tensor((*diffusion_mat), 0, pde_inf.num_pde - 1, 0, pde_inf.num_step,
                  1, pde_inf.num_mesh);
  (*diffusion_mat) = NULL;
  if ((*convection_mat))
    free_f3tensor((*convection_mat), 0, pde_inf.num_pde - 1, 0,
                  pde_inf.num_step, 1, pde_inf.num_mesh);
  (*convection_mat) = NULL;
  if ((*dividend_tensor))
    free_f3tensor((*dividend_tensor), 0, pde_inf.num_pde - 1, 0,
                  pde_inf.num_step, 1, pde_inf.num_mesh);
  (*dividend_tensor) = NULL;

  return NULL;
}

static Err free_pde_inf_ptr(SrtBasicPdeInf *pde_inf) {
  long num_mesh, num_time_step;

  num_mesh = (*pde_inf).num_mesh;
  num_time_step = (*pde_inf).num_step;

  free_dvector((*pde_inf).log_spots, 1, num_mesh - 1);
  (*pde_inf).log_spots = NULL;
  free_dvector((*pde_inf).spots, 1, num_mesh - 1);
  (*pde_inf).spots = NULL;

  free_dvector((*pde_inf).time, 0, num_time_step - 1);
  (*pde_inf).time = NULL;

  free_dvector((*pde_inf).basevol, 0, num_time_step - 1);
  (*pde_inf).basevol = NULL;
  free_dvector((*pde_inf).beta, 0, num_time_step - 1);
  (*pde_inf).beta = NULL;

  free_dvector((*pde_inf).gamma, 0, num_time_step - 1);
  (*pde_inf).gamma = NULL;
  free_dvector((*pde_inf).omega, 0, num_time_step - 1);
  (*pde_inf).omega = NULL;

  free_dvector((*pde_inf).voldrift, 0, num_time_step - 1);
  (*pde_inf).voldrift = NULL;
  free_dvector((*pde_inf).vovol, 0, num_time_step - 1);
  (*pde_inf).vovol = NULL;

  free_dvector((*pde_inf).cxxorrelation_und_vol, 0, num_time_step - 1);
  (*pde_inf).cxxorrelation_und_vol = NULL;
  free_dvector((*pde_inf).cxxorrelation_und_rate, 0, num_time_step - 1);
  (*pde_inf).cxxorrelation_und_rate = NULL;

  free_dvector((*pde_inf).forward, 0, num_time_step - 1);
  (*pde_inf).forward = NULL;

  return NULL;
}

static Err fill_pde_payoff(SrtBasicPdeInf pde_inf, double strike,
                           SrtReceiverType rec_pay, long pde_index,
                           double **payoff) {

  long i;
  Err err = NULL;

  for (i = 1; i < pde_inf.num_mesh; i++)
    (*payoff)[i] = 0.0;

  return NULL;
}

/* next step : implement pde fwd */
Err srt_f_stoch_rate_and_vol_gamma_smile_pde(Date event_date,

                                             double max_time, long min_node,
                                             long min_num_mesh,

                                             char *und_name, char *dom_und_name,

                                             char *yc_name,

                                             double strike, char *rec_pay_str,
                                             char *greek_str, double *greeks)

{
  Err err = NULL;
  SrtMdlType dom_mdl_type, und_mdl_type;
  SrtBasicPdeInf pde_inf;
  SrtUndPtr und, dom_und;
  SrtGreekType greek_type;
  SrtReceiverType rec_pay;

  long k;
  double forward, black_vol;
  double ***diffusion_mat, ***convection_mat, ***dividend_tensor, *payoff,
      *option_pv_vect, pde_bounds[2];
  Date date_today;

  err = interp_rec_pay(rec_pay_str, &rec_pay);
  if (err)
    return err;

  /* Get the Underlyings pointer and the model types */
  und = lookup_und(und_name);
  err = get_underlying_mdltype(und, &und_mdl_type);
  if (err)
    return err;

  dom_und = lookup_und(dom_und_name);
  err = get_underlying_mdltype(dom_und, &dom_mdl_type);
  if (err)
    return err;

  /* Builds the pde: localisation and meshing of the state space */

  err = und_pde_lim(event_date, und, dom_und, max_time, min_node, min_num_mesh,
                    &pde_inf);

  if (err)
    return err;

  err = make_pde_operators_allocations(pde_inf, &payoff, &diffusion_mat,
                                       &convection_mat, &dividend_tensor,
                                       &option_pv_vect);
  if (err)
    return err;

  err = srt_f_set_srvgs_pde_operators(event_date, dom_und_name, und_name,
                                      pde_inf, strike, rec_pay, &diffusion_mat,
                                      &convection_mat, &dividend_tensor);
  if (err)
    return err;

  err = interp_greeks(greek_str, &greek_type);
  if (err)
    return err;

  /* Get the forward */
  err = srt_f_eq_forward(event_date, und, yc_name, &forward);
  if (err)
    return err;

  /* Get the date today */
  date_today = get_today_from_underlying(und);

  err = srt_f_eq_implied_vol((event_date - date_today) * YEARS_IN_DAY, und_name,
                             &black_vol);

  if (err)
    return err;

  (*greeks) = srt_f_optblksch(forward, strike, black_vol,
                              (event_date - date_today) * YEARS_IN_DAY, 1.0,
                              rec_pay, greek_type);

  /* solve sucessively the pde */
  for (k = 0; k < pde_inf.num_pde; k++) {
    err = fill_pde_payoff(pde_inf, strike, rec_pay, k, &payoff);
    if (err)
      return err;

    pde_bounds[0] = pde_bounds[1] = 0.0;

    err = srt_f_feynmann_kac_pde(pde_inf.num_mesh, pde_inf.num_step, pde_inf.h,
                                 pde_inf.time_step, pde_inf.log_spots,
                                 convection_mat[k], diffusion_mat[k], LINEAR,
                                 DIRICHLET, pde_bounds, payoff, NULL,
                                 dividend_tensor[k], &option_pv_vect);
    if (err)
      return err;

    if (greek_type == PREMIUM) {
      (*greeks) += option_pv_vect[pde_inf.spot_index];
    } else if (greek_type == DELTA) {
      (*greeks) += (option_pv_vect[pde_inf.spot_index + 1] -
                    option_pv_vect[pde_inf.spot_index]) /
                   pde_inf.h;

    } else if (greek_type == GAMMA) {
      (*greeks) += (option_pv_vect[pde_inf.spot_index + 1] -
                    2 * option_pv_vect[pde_inf.spot_index] +
                    option_pv_vect[pde_inf.spot_index - 1]) /
                   (pde_inf.h * pde_inf.h);

    } else
      return serror("Can not display greeks");
  }

  err = make_pde_operators_desallocations(pde_inf, &payoff, &diffusion_mat,
                                          &convection_mat, &dividend_tensor,
                                          &option_pv_vect);
  if (err)
    return err;

  err = free_pde_inf_ptr(&pde_inf);
  if (err)
    return err;

  return NULL;
}

#undef MAX_NUM_ITER
