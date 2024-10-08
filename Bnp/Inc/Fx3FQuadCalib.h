
#ifndef Fx3FQuadCalib_h
#define Fx3FQuadCalib_h

Err FxQuad_log_approx(long today, double *maturity, long nb_mat,
                      double *sig_dom, double *sig_for, double *mat_fx,
                      long nb_mat_fx, double *sig_fx, double alpha, double beta,
                      double gamma, double sig0, double lam_dom, double lam_for,
                      double corr_dom_for, double corr_dom_fx,
                      double corr_for_fx, double spot, char *dom_yc,
                      char *for_yc, double *time_fx, long nb_time_fx,
                      double *fx_fwd, double *fx_vol, double *fx_fwd0,
                      double max_time);

Err Fx3DQuadtsImpliedVol(long today, double opt_maturity, double start_date,
                         double end_date, double *maturity, long nbMat,
                         double *sig_curve_dom, double lda_dom,
                         double *sig_curve_for, double lda_for, double *mat_fx,
                         long nb_mat_fx, double *sig_curve_fx, double alpha,
                         double beta, double gamma, double sig0, double spot_fx,
                         double corr_dom_for, double corr_dom_fx,
                         double corr_for_fx, char *dom_yc, char *for_yc,
                         double *fx_vol, double disc_dt, double fx_dt);

Err Fx3DQuadImpliedVol(char *fx_underlying, double beta, double gamma,
                       double sig0, double val_time, double start_time,
                       double end_time, double disc_dt, double fx_dt,
                       double *vol);

Err Fx3DQuadtsCalibration(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double *params),
    double *param, double spot_fx, double correl_dom_for, double correl_dom_fx,
    double correl_for_fx, char *dom_yc, char *for_yc, double **fx_vol_curve,
    long nbSteps, long nbNewton, double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DQuadCalibration(char *dom_underlying, char *for_underlying,
                        double spot_fx, double beta, double gamma,
                        double correl_dom_for, double correl_dom_fx,
                        double correl_for_fx, double *exercise_opt,
                        double *maturity_opt, double *vol_opt, long nbropt,
                        long nbLong, double **fx_vol_curve, long nbSteps,
                        long nbNewton, double disc_dt, double fx_dt,
                        long nbIterMax);

Err Fx3DQuadtsTreeFxOptions(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double gamma, double sig0,
    double spot_fx, double corr_dom_for, double corr_dom_fx, double corr_for_fx,
    char *dom_yc, char *for_yc, double *vols, double *vols_time, int nbVols,
    double *option_prices, long num_stp);

#endif