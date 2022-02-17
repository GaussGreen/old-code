
#ifndef Fx3FLocalCalib_h
#define Fx3FLocalCalib_h

Err FxLocal_log_approx(
    long today, double *maturity, long nb_mat, double *sig_dom, double *sig_for,
    double *mat_fx, long nb_mat_fx, double *sig_fx,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double *params),
    double *params, double lam_dom, double lam_for, double corr_dom_for,
    double corr_dom_fx, double corr_for_fx, double spot, char *dom_yc,
    char *for_yc, double *time_fx, long nb_time_fx, double *fx_fwd,
    double *fx_vol, double *fx_fwd0, double max_time);

Err Fx3DLocaltsImpliedVol(
    long today, double opt_maturity, double start_date, double end_date,
    double *maturity, long nbMat, double *sig_curve_dom, double lda_dom,
    double *sig_curve_for, double lda_for, double *mat_fx, long nb_mat_fx,
    double *sig_curve_fx,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double *params),
    double *params, double spot_fx, double corr_dom_for, double corr_dom_fx,
    double corr_for_fx, char *dom_yc, char *for_yc, double *fx_vol,
    double disc_dt, double fx_dt);

Err Fx3DLocalImpliedVol(char *fx_underlying,
                        double (*vol_ln_func)(double t, double Spot, double Fwd,
                                              double *params),
                        double *params, double val_time, double start_time,
                        double end_time, double disc_dt, double fx_dt,
                        double *vol);

Err Fx3DLocaltsCalibration(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, double *maturity, long nbrMat, double *sig_curve_dom,
    double lda_dom, double *sig_curve_for, double lda_for,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double *params),
    double *params, double spot_fx, double correl_dom_for, double correl_dom_fx,
    double correl_for_fx, char *dom_yc, char *for_yc, double **fx_vol_curve,
    double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DLocalCalibration(
    char *dom_underlying, char *for_underlying, double spot_fx,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double *params),
    double *params, double correl_dom_for, double correl_dom_fx,
    double correl_for_fx, double *exercise_opt, double *maturity_opt,
    double *vol_opt, long nbropt, double **fx_vol_curve, double disc_dt,
    double fx_dt, long nbIterMax);

#endif