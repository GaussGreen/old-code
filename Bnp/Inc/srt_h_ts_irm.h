/* ============================================================================

   FILENAME:   SRT_H_TS_IRM.H

   PURPOSE:    Header file for all the functions used to initialise and handle
               a ONE FACTOR interest rate model Term Structure

   ============================================================================
 */

#ifndef SRT_H_TS_IRM_H
#define SRT_H_TS_IRM_H

#include "srt_h_ts.h"

Err srt_f_get_vasicek_mean_sr(double time, TermStruct *ts,
                              double *vasicek_mean_sr);

Err srt_f_get_vasicek_mean_int_sr(double time, TermStruct *ts,
                                  double *vasicek_mean_int_sr);

Err srt_get_vasicek_risk_neutral_sr_mean(double time, TermStruct *ts,
                                         double *vasicek_mean_sr);

Err srt_get_lgm_int_zc_bond(double bond_mat, double time, TermStruct *ts,
                            double *int_zc_bond);

Err srt_f_get_vasicek_var_sr(double time, TermStruct *ts, double *var_sr);

Err srt_f_vasicek_discount_factor(Date discount_mat_date, SrtUndPtr und,
                                  double *discount_factor);

Err srt_f_vasicek_fra(Date start_date, Date end_date, SrtUndPtr und,
                      double *fra);

Err srt_f_vasicek_swap(Date start_date, Date theo_end_date, char *swap_freq_str,
                       char *swap_basis_str, SrtUndPtr und, double *swap);

Err srt_f_vasicek_fut(double start_date, double pay_date, char *und_name,
                      double *fut_pr);

Err srt_f_vasicek_market_libor_expectation(double start_date, double pay_date,
                                           double risk_premuim, char *und_name,
                                           double *market_libor_expectation);

Err srt_f_vasicek_overnight_interest_rate_swap(
    char *und_name, long start_date, long end_date,

    double *overnight_interest_rate_swap);

Err srt_f_vasicek_cont_fwd_zr(double fixing_date, double pay_date,
                              SrtUndPtr und, double *cont_fwd_zr);

Err srt_f_vasicek_dirty_pr(SrtUndPtr und, long num_dates, double settlt_date,
                           double *pay_dates, double *cpns,
                           SRT_Boolean discount_fwd_dirty_pr,
                           double *dirty_price);

Err srt_f_vasicek_calibrate(double settlt_date, double vasicek_init_cond,
                            long num_bonds, long *num_cpns, double **pay_dates,
                            double **cpns,

                            double *mkt_dirty_prices,
                            SRT_Boolean discount_fwd_dirty_pr,

                            SRT_Boolean price_floor, double *floor_gearing,
                            double *pre_factor_fwd, double *floor_strike,
                            double *floor_implied_vol,

                            char **vasicek_und_name);

Err srt_f_vasicek_overnights_calibrate(long num_overnight_swap_rates,
                                       Date *overnights_swap_end_dates,
                                       double *market_overnight_swap_rates,
                                       double *overnight_swap_rates_weights,
                                       char **und_name);

Err srt_f_vasicek_calibrated_mean_rev_lvl(
    double settlt_date, double vasicek_init_cond, long num_bonds,
    long *num_cpns, double **pay_dates, double **cpns,

    double *bond_dirty_prices, SRT_Boolean discount_fwd_dirty_pr,

    SRT_Boolean price_floor, double *floor_gearing, double *pre_factor_fwd,
    double *floor_strike, double *floor_implied_vol,

    char *vasicek_und_name, double ***calibrated_mean_rev_lvl);

Err srt_f_get_vasicek_init_cond(double time, TermStruct *ts,
                                double *vasicek_init_cond);

Err srt_f_vasicek_money_market_account_vol(Date date, char *und_name,
                                           double *money_market_account_vol);

Err srt_f_inflation_fwd(long start_date, long end_date,
                        char *price_index_und_name, double *inflation_fwd);

Err srt_f_lgm_phi_func(Date date, char *und_name, double *phi_val);

/* -------------------------------------------------------------------------
   The MAIN function to call when a Interest Rate Model TermStructure has
   to be initialised
   ------------------------------------------------------------------------- */

Err srt_f_init_IRM_TermStruct(Date today, double **sig, /* sig[col][row] */
                              int sig_cols, int nsig,
                              double **tau, /* tau[col][row] */
                              int tau_cols, int ntau,

                              SrtMdlType mdl_type, SrtMdlDim mdl_dim,

                              /* SMILE */
                              double beta,

                              /* TWO FACTOR */
                              double alpha, double gamma, double rho,

                              /* STOCH VOL */
                              double vovol,

                              /* ETABETA  or MIXED */
                              double eta_or_omega,

                              /* VASICEK MODEL */
                              double vasicek_init_cond, int num_mean_rev_level,
                              int num_mean_rev_level_cols,
                              double **vasicek_mean_rev_level_data,

                              /* OUTPUT */
                              TermStruct **ts);

/* -------------------------------------------------------------------------- */

/* The MAIN function to FREE an Interest Rate Model Term Structure */

Err srt_f_free_IRM_TermStruct(TermStruct **l);

/* Function required to free a TermStructObject once attached into a linked list
 */
Err srt_f_irmtsvalfree(void *tsvalptr);

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
                  ONE FACTOR SPECIFIC FUNCTIONS
   ------------------------------------------------------------------------- */

/*
   The MAIN function to call when a Interest Rate Model TermStructure has
   to be initialised with a ONE FACTOR MODEL   */

Err srt_f_init_IRM_OneFac_TermStruct(
    TermStruct **ts, Date today, double **sig_data, int sig_cols, int num_sig,
    double **tau_data, int tau_cols, int num_tau, SrtMdlType mdl_type,
    double vovol, double rho, double meanvol, double beta,
    double vasicek_init_cond, int vasicek_mean_rev_level_rows,
    int vasicek_mean_rev_level_cols, double **vasicek_mean_rev_level_vals);

/* -------------------------------------------------------------------------- */

Err srt_f_display_IRM_OneFac_TermStruct(TermStruct *ts, double **sigma_date,
                                        double **sigma, double **beta,
                                        double **vovol, double **rho,
                                        double **meanvol, long *plNumSigmas,
                                        double **tau_date, double **tau,
                                        long *plNumTaus);

/* -------------------------------------------------------------------------- */

Err make_Lambda_vector(TermStruct *l);

Err make_F_Psi_vector(TermStruct *l);

Err make_G_H_vector(TermStruct *l);

/* -------------------------------------------------------------------------- */

double find_tau(double time, TermStruct *l);

/* -> return from TermStruct 'l' the tau value corresponding with 'time' */

/* -------------------------------------------------------------------------- */

double find_sig(double time, TermStruct *l);

/* -> return from TermStruct 'l' the sigma value corresponding with 'time' */

/* -------------------------------------------------------------------------- */

double find_vovol(double time, TermStruct *l);

/* -> return from TermStruct 'l' the vovvol value corresponding with 'time' */

/* -------------------------------------------------------------------------- */

double find_rho(double time, TermStruct *l);

/* -> return from TermStruct 'l' the correlation value corresponding with 'time'
 */

/* -------------------------------------------------------------------------- */

double find_meanvol(double time, TermStruct *l);

/* -> return from TermStruct 'l' the correlation value corresponding with 'time'
 */

/* -------------------------------------------------------------------------- */

double **find_M_eta_beta(double time, TermStruct *l);

/* -> return from 'l' the M matrix corresponding with 'time' for betaeta */

/* -------------------------------------------------------------------------- */

double find_beta(double time, TermStruct *l);
/* -> return from TermStruct 'power' the power value corresponding to the state
 * var */

/* -------------------------------------------------------------------------- */

double find_eta(double time, TermStruct *l);
/* -> return from TermStruct 'power' the power value corresponding to the state
 * var */

/* -------------------------------------------------------------------------- */

double find_vasicek_init_cond(double time, TermStruct *l);
/* -> return from TermStruct 'power' the power value corresponding to the state
 * var */

/* -------------------------------------------------------------------------- */

double find_mean_rev_level(double time, TermStruct *l);
/* -> return from TermStruct 'power' the power value corresponding to the state
 * var */

/* -------------------------------------------------------------------------- */

double find_F(double time, TermStruct *l);
double find_G(double time, TermStruct *l);
double find_H(double time, TermStruct *l);
double find_Psi(double time, TermStruct *l);
double find_I(double time, TermStruct *l);
double find_J(double time, TermStruct *l);
double find_K(double time, TermStruct *l);
double find_L(double time, TermStruct *l);
double find_O(double time, TermStruct *l);
double find_Q(double time, TermStruct *l);
double find_Phi(double time, TermStruct *l);

/*-------------------------------------------------------------------------------*/
double F_func(double time, TermStruct *l);

/* -> return ? */

/* -------------------------------------------------------------------------- */

/*-------------------------------------------------------------------------------*/
double T_func(double time, TermStruct *l);

/* -> return ? */

/* -------------------------------------------------------------------------- */

double Psi_func(double time, TermStruct *l);
/* -> return ? */

/* -------------------------------------------------------------------------- */

double Kappa_func(double time, TermStruct *l);
/* -> return ? */

/* -------------------------------------------------------------------------- */

double J_func(double time, TermStruct *l);
/* -> calculates the exact value of the J function at date using J_vector */

/* -------------------------------------------------------------------------- */

double Zeta_func(double time, TermStruct *l);

/* -> return Zeta for the betaeta model */

/* -------------------------------------------------------------------------- */

double Lambda_func(double starttime, double endtime, TermStruct *l);

/* -> returns the variations of Psi */
/* -------------------------------------------------------------------------- */

double M_eta_beta_func(double s, double theta, double eta, double **M);
/*BigM        *Mstruct);*/

/* -> returns the linearly interpolated M for the betaeta model */

/* -------------------------------------------------------------------------- */

double dln_F_func_dt(double time, TermStruct *l);

/* -> return ? */

/* -------------------------------------------------------------------------- */

double find_sig2_interp(double time1, double time2, TermStruct *l);

/* -> return interpolated value between 'time1' and 'time2'
        We suppose that time1 < time2			    */

/* -------------------------------------------------------------------------- */

void G_H_func(double t, TermStruct *l, double *val_G, double *val_H);

/* -------------------------------------------------------------------------- */

Err srt_f_tsupdate(TermStruct **l);

/* -> update TermStruct 'l': fills blank cells      ,
                             recomputes F      ,Psi      ,G and H */

/* -------------------------------------------------------------------------- */

Err srt_f_tsaddtau(TermStruct *list, Ddate date, double today, double tau);

/* -> add to the TermStruct 'list a date 'date' with tau value 'tau' */

/* -------------------------------------------------------------------------- */

Err srt_f_tsaddsig(TermStruct *list, Ddate date, double today, double sig);

/* -> add to the TermStruct 'list' a date 'date' with sigma value 'sig' */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Computes the LGM cumulative volatility between two dates t1 and t2      ,
   defined by:
                CumVol(t      , T1      , T2) = [Psi(T2)-Psi(T1)]^2 * G(t)
   ------------------------------------------------------------------------- */

double srt_f_lgm_cum_vol(TermStruct *ts, double t, double t1, double t2,
                         double *answer);

/* -------------------------------------------------------------------------- */

/* For Quanto Adjustment in LGM Jumping : computes the M function */

Err srt_f_extend_lgm_jumping_ts_for_quanto(TermStruct *lgm_ts,
                                           TermStruct *fx_ts,
                                           char *szLgmUndName,
                                           char *szFxUndName);

/* FOr Quanto Adjustments in LGM Jumping */

double M_func(double time, TermStruct *lgm_ts);

/* ---------------------------------------------------------------------------
 */
/* For the NEWLGM Model */

typedef enum { G = 0, H = 1 } LabelStruct;

double find_struct_interp(double time, LabelStruct label, TermStruct *l);
double find_struct_interp_der(double time, LabelStruct label, TermStruct *l);

#endif
