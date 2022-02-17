
#ifndef Fx3FBetaMcHelpH
#define Fx3FBetaMcHelpH

Err fill_Betamc_init(double *date, double *time, long nb_dates,
                     double *sig_dates, long nb_sig_dates,
                     double *sig_curve_dom, double lda_dom,
                     double *sig_curve_for, double lda_for,
                     double *sig_curve_fx, char *dom_yc, char *for_yc,
                     double *dom_ifr, double *dom_std, double *dom_phi,
                     double *for_ifr, double *for_std, double *for_phi,
                     double *fx_std);

#endif
