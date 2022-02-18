
#ifndef Fx5FMCH
#define Fx5FMCH

typedef struct
{
    GRFNCOMMSTRUCT global;
    FIRSTMktAtT*   local;

    int fx_idx;
    int dom_idx; /* -1 if none */
    int for_idx; /* -1 if none */

    long    num_fx_df;
    double* fx_df_tms;
    long*   fx_df_dts;
    double* fx_dff;
    double* fx_gam1;
    double* fx_gam2;
    double* fx_gam12;
    int     do_fx;

    long    num_dom_df;
    double* dom_df_tms;
    long*   dom_df_dts;
    double* dom_dff;
    double* dom_gam1;
    double* dom_gam2;
    double* dom_gam12;
    int     do_dom;

    long    num_for_df;
    double* for_df_tms;
    long*   for_df_dts;
    double* for_dff;
    double* for_gam1;
    double* for_gam2;
    double* for_gam12;
    int     do_for;

} grfn_parm_mc5F, *GRFNPARMMC5F;

//Err mc_main_5dfx(
//    /*	Time data */
//    long      npaths,
//    int       num_col,
//    double*   time,
//    double*   date,
//    long      nb_dates,
//    double*   dom_ifr,
//    double*   dom_fwd1,
//    double*   dom_fwd2,
//    double*   dom_exp1,
//    double*   dom_exp2,
//    double*   dom_phi1,
//    double*   dom_phi2,
//    double*   dom_phi12,
//    double*   dom_gam1_fwd,
//    double*   dom_gam2_fwd,
//    double*   dom_bond_pay,
//    double*   dom_gam1_pay,
//    double*   dom_gam2_pay,
//    double*   for_ifr,
//    double*   for_fwd1,
//    double*   for_fwd2,
//    double*   for_exp1,
//    double*   for_exp2,
//    double*   for_phi1,
//    double*   for_phi2,
//    double*   for_phi12,
//    double*   for_gam1_fwd,
//    double*   for_gam2_fwd,
//    double*   fx_fwd,
//    double*** covar,
//    /*	Product data */
//    void** func_parm_tab,
//    /*	Model data */
//    double   dom_lam,
//    double   dom_alpha,
//    double   dom_gamma,
//    double   for_lam,
//    double   for_alpha,
//    double   for_gamma,
//    double** correl,
//    /*	Market data */
//    double spot_fx,
//    char*  dom_yc,
//    char*  for_yc,
//    /* do PECS adjustment */
//    int do_pecs,
//    /* for Optimisation of exercise boundary */
//    int        do_optimisation,
//    int*       optimise,
//    MCEBPARAMS params,
//    /*	Initialisation function to be called at the beggining of each path
//                    or NULL if none */
//    void (*init_func)(),
//    /*	Payoff function */
//    Err (*payoff_func)(
//        /* Event */
//        double evt_date,
//        double evt_time,
//        void*  func_parm,
//        /* Market data */
//        double spot_fx,
//        double R1D,
//        double R2D,
//        double R1F,
//        double R2F,
//        double Z,
//        /* Results */
//        int     num_col,
//        double* res,
//        int*    stop_path),
//    /*	Results */
//    double** res);
//
///*	Monte Carlo in the 5F model							*/
///*	Optimise the payoff in column col_pay				*/
///*	against an exercise boundary in column col_bound	*/
///*	------------------------------------------------	*/
//Err mc_exe_bound_5dfx(
//    /*	Time data */
//    long      npaths,
//    int       num_col,
//    int       col_pay,
//    int       col_bound,
//    int*      optimise,
//    int       call_current,
//    int       is_ko,
//    double*   time,
//    double*   date,
//    long      nb_dates,
//    double*   dom_ifr,
//    double*   dom_fwd1,
//    double*   dom_fwd2,
//    double*   dom_exp1,
//    double*   dom_exp2,
//    double*   dom_phi1,
//    double*   dom_phi2,
//    double*   dom_phi12,
//    double*   dom_gam1_fwd,
//    double*   dom_gam2_fwd,
//    double*   dom_bond_pay,
//    double*   dom_gam1_pay,
//    double*   dom_gam2_pay,
//    double*   for_ifr,
//    double*   for_fwd1,
//    double*   for_fwd2,
//    double*   for_exp1,
//    double*   for_exp2,
//    double*   for_phi1,
//    double*   for_phi2,
//    double*   for_phi12,
//    double*   for_gam1_fwd,
//    double*   for_gam2_fwd,
//    double*   fx_fwd,
//    double*** covar,
//    /*	Product data */
//    void** func_parm_tab,
//    /*	Model data */
//    double   dom_lam,
//    double   dom_alpha,
//    double   dom_gamma,
//    double   for_lam,
//    double   for_alpha,
//    double   for_gamma,
//    double** correl,
//    /*	Market data */
//    double spot_fx,
//    char*  dom_yc,
//    char*  for_yc,
//    /* do PECS adjustment */
//    int do_pecs,
//    /*	Initialisation function to be called at the beggining of each path
//                    or NULL if none */
//    void (*init_func)(),
//    /*	Payoff function */
//    Err (*payoff_func)(
//        /* Event */
//        double evt_date,
//        double evt_time,
//        void*  func_parm,
//        /* Market data */
//        double spot_fx,
//        double R1D,
//        double R2D,
//        double R1F,
//        double R2F,
//        double Z,
//        /* Results */
//        int     num_col,
//        double* res,
//        int*    stop_path),
//    /*	Results */
//    double** res);

Err mc_main_5dfx_test(
    /*	Time data */
    long    npaths,
    int     num_col,
    double* time,
    double* date,
    long    nb_dates,
    double* dom_ifr, /*	Distributions */
    double* dom_fwd,
    double* dom_std,
    double* dom_phi1,
    double* dom_phi2,
    double* dom_phi12,
    double* dom_beta,
    double* dom_bond_pay,
    double* dom_beta_pay,
    double* for_ifr,
    double* for_fwd,
    double* for_std,
    double* for_phi1,
    double* for_phi2,
    double* for_phi12,
    double* for_beta,
    double* fx_fwd,
    double* fx_std,
    double* dom_for_cov,
    double* dom_fx_cov,
    double* for_fx_cov,
    /*	Product data */
    void** func_parm_tab,
    /*	Model data */
    double   dom_lam,
    double   dom_alpha,
    double   dom_gamma,
    double   dom_rho,
    double   for_lam,
    double   for_alpha,
    double   for_gamma,
    double   for_rho,
    double** correl,
    /*	Market data */
    double spot_fx,
    char*  dom_yc,
    char*  for_yc,
    /* do PECS adjustment */
    int do_pecs,
    /*	Initialisation function to be called at the beggining of each path
                    or NULL if none */
    void (*init_func)(),
    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date,
        double evt_time,
        void*  func_parm,
        /* Market data */
        double spot_fx,
        double R1D,
        double R2D,
        double R1F,
        double R2F,
        double Z,
        /* Results */
        int     num_col,
        double* res,
        int*    stop_path),
    /*	Results */
    double** res);

Err fill_mc_init5Fsimple(
    long    pay_date,
    double  pay_time,
    double* date,
    double* time,
    long    nb_dates,
    double* sig_dates,
    long    nb_sig_dates,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* sig_curve_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    char*   dom_yc,
    char*   for_yc,
    double* dom_ifr,
    double* dom_fwd,
    double* dom_std,
    double* dom_phi,
    double* dom_beta,
    double* dom_bond_pay,
    double* dom_beta_pay,
    double* for_ifr,
    double* for_fwd,
    double* for_std,
    double* for_phi,
    double* for_beta,
    double* fx_fwd,
    double* fx_std,
    double* dom_for_cov,
    double* dom_fx_cov,
    double* for_fx_cov);

Err fill_mc_init5F(
    long      pay_date,
    double    pay_time,
    double*   date,
    double*   time,
    long      nb_dates,
    double*   sig_dates,
    long      nb_sig_dates,
    double*   sig_curve_dom,
    double    lda_dom,
    double    alpha_dom,
    double    gamma_dom,
    double*   sig_curve_for,
    double    lda_for,
    double    alpha_for,
    double    gamma_for,
    double*   sig_curve_fx,
    double**  correl,
    char*     dom_yc,
    char*     for_yc,
    double*   dom_ifr,
    double*   dom_fwd1,
    double*   dom_fwd2,
    double*   dom_exp1,
    double*   dom_exp2,
    double*   dom_phi1,
    double*   dom_phi2,
    double*   dom_phi12,
    double*   dom_gam1_fwd,
    double*   dom_gam2_fwd,
    double*   dom_bond_pay,
    double*   dom_gam1_pay,
    double*   dom_gam2_pay,
    double*   for_ifr,
    double*   for_fwd1,
    double*   for_fwd2,
    double*   for_exp1,
    double*   for_exp2,
    double*   for_phi1,
    double*   for_phi2,
    double*   for_phi12,
    double*   for_gam1_fwd,
    double*   for_gam2_fwd,
    double*   fx_fwd,
    double*** covar);
#endif
