
#ifndef LGM2FMCH
#define LGM2FMCH

#include "MCEBOptimisation.h"

typedef struct {
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  long num_dom_df;
  double *dom_df_tms;
  long *dom_df_dts;
  double *dom_dff;
  double *dom_gam1;
  double *dom_gam2;
  double *dom_gam12;
  int do_dom;

} grfn_parm_mc2F, *GRFNPARMMC2F;

/*	Main function. Diffusion can be made either under QTfinal or under
 * Jumping Numeraire */
/*	------------------------------------------------------------------------------------
 */
Err mc_main_lgm2f(
    /*	Time data */
    long npaths, int num_col, double *time, double *date, long nb_dates,

    /*	Model data */
    int jumping_num, double *dom_fwd1, double *dom_fwd2, double *dom_exp1,
    double *dom_exp2, double *dom_phi1, double *dom_phi2, double *dom_phi12,
    double *dom_gam1_fwd, double *dom_gam2_fwd, double *dom_bond_pay,
    double *dom_gam1_pay, double *dom_gam2_pay, double ***covar,

    /*	Product data */
    void **func_parm_tab,

    /* do PECS adjustment */
    int do_pecs,

    /* for Optimisation of exercise boundary */
    int do_optimisation, int *optimise, MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path
                    or NULL if none */
    void (*init_func)(),
    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date, double evt_time, void *func_parm,
        /* Market data */
        double R1D, double R2D,
        /* Results */
        int num_col, double *res, int *stop_path),
    /*	Results */
    double **res);

Err fill_mc_init_lgm2f(
    int do_jump, long pay_date, double pay_time, double *date, double *time,
    long nb_dates, double *sig_dates, long nb_sig_dates, double *sig_curve_dom,
    double lda_dom, double alpha_dom, double gamma_dom, double rho_dom,
    char *dom_yc, double *dom_fwd1, double *dom_fwd2, double *dom_exp1,
    double *dom_exp2, double *dom_phi1, double *dom_phi2, double *dom_phi12,
    double *dom_gam1_fwd, double *dom_gam2_fwd, double *dom_bond_pay,
    double *dom_gam1_pay, double *dom_gam2_pay, double ***covar);

// Merges the term structures of lambda and sigma
Err merge_lambda_sigma_ts(double *sigma, double *sigma_time, int nb_sigma,
                          double *lambda, double *lambda_time, int nb_lambda,
                          double **ts_time, double **lam, double **sig,
                          int *nb_new_time);

// Calculates the fwds      , variances      , phis etc. for a LGM2F with lambda
// term structure
Err fill_mc_init_lgm2f_lambda(
    int do_jump, long pay_date, double pay_time, double *date, double *time,
    long nb_dates, double *sig_dates, long nb_sig_dates, double *sig_curve_dom,
    double *lda_dom, double *lda_dom2, double alpha_dom, double gamma_dom,
    double rho_dom, char *dom_yc, double *dom_fwd1, double *dom_fwd2,
    double *dom_exp1, double *dom_exp2, double *dom_phi1, double *dom_phi2,
    double *dom_phi12, double *dom_gam1_fwd, double *dom_gam2_fwd,
    double *dom_bond_pay, double *dom_gam1_pay, double *dom_gam2_pay,
    double ***covar);

// Calculates gamma_lambda = integral(s      ,T      ,exp(-integral(s      ,u
// ,lambda(w)
//     ,w)      ,du when lambda(w) is piecewise constant
Err gamma_lambda(double s, double T, double *lambda_time, double *lambda,
                 long n_lambda, double *result);

// Calculates exp_lambda = exp(-integral(s      ,T      ,lambda(w)      ,w))
// when lambda(w) is piecewise constant
Err exp_lambda(double s, double T, double *lambda_time, double *lambda,
               long n_lambda, double *result);

// Calculates intermediary results for the calculation of the forwards
Err calculate_Bk_rhoCk(double Tk, double Tk1, long k, double Tn1,
                       double *sigma_dates, long nb_sigma_dates,
                       double *sigma_dom, double *lda_dom, double *lda_dom2,
                       double alpha_dom, double gamma_dom, double rho_dom,
                       double *Bk_rhoCk_fwd1, double *Bk_rhoCk_fwd2);

// Calculates intermediary results for the calculation of the forwards
Err calculate_Dk_rhoEk(double Tk, double Tk1, long k, double Tpay, long iTpay,
                       double Tn1, double *sigma_dates, long nb_sigma_dates,
                       double *sigma_dom, double *lda_dom, double *lda_dom2,
                       double alpha_dom, double gamma_dom, double rho_dom,
                       double *Dk_rhoEk_fwd1, double *Dk_rhoEk_fwd2);

#endif
