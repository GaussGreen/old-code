
#ifndef Fx3FMcHelpH
#define Fx3FMcHelpH

Err fill_mc_init(long pay_date, double pay_time, double *date, double *time,
                 long nb_dates, double *sig_dates, long nb_sig_dates,
                 double *sig_curve_dom, double lda_dom, double *sig_curve_for,
                 double lda_for, double *sig_curve_fx, double correl_dom_for,
                 double correl_dom_fx, double correl_for_fx, char *dom_yc,
                 char *for_yc, double *dom_ifr, double *dom_fwd,
                 double *dom_std, double *dom_phi, double *dom_beta,
                 double *dom_bond_pay, double *dom_beta_pay, double *for_ifr,
                 double *for_fwd, double *for_std, double *for_phi,
                 double *for_beta, double *fx_fwd, double *fx_std,
                 double *dom_for_cov, double *dom_fx_cov, double *for_fx_cov);

Err fill_mc_init_corr(long pay_date, double pay_time, double *date,
                      double *time, long nb_dates, double *sig_dates,
                      long nb_sig_dates, double *sig_curve_dom, double lda_dom,
                      double *sig_curve_for, double lda_for,
                      double *sig_curve_fx, double *correl_dom_for_ts,
                      double *correl_dom_fx_ts, double *correl_for_fx_ts,
                      char *dom_yc, char *for_yc, double *dom_ifr,
                      double *dom_fwd, double *dom_std, double *dom_phi,
                      double *dom_beta, double *dom_bond_pay,
                      double *dom_beta_pay, double *for_ifr, double *for_fwd,
                      double *for_std, double *for_phi, double *for_beta,
                      double *fx_fwd, double *fx_std, double *dom_for_cov,
                      double *dom_fx_cov, double *for_fx_cov);

#endif
