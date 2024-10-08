
#ifndef Fx3FBetaCalib_h
#define Fx3FBetaCalib_h

Err Fxbeta_log_approx(long today, double *maturity, long nb_mat,
                      double *sig_dom, double *sig_for, double *mat_fx,
                      long nb_mat_fx, double *sig_fx, double *beta,
                      double lam_dom, double lam_for, double corr_dom_for,
                      double corr_dom_fx, double corr_for_fx, double spot,
                      char *dom_yc, char *for_yc, double *time_fx,
                      long nb_time_fx, double *fx_fwd, double *fx_vol,
                      double max_time);

Err Fx3DBetatsImpliedVol(long today, double opt_maturity, double start_date,
                         double end_date, double *maturity, long nbMat,
                         double *sig_curve_dom, double lda_dom,
                         double *sig_curve_for, double lda_for, double *mat_fx,
                         long nb_mat_fx, double *sig_curve_fx, double *beta,
                         double spot_fx, double corr_dom_for,
                         double corr_dom_fx, double corr_for_fx, char *dom_yc,
                         char *for_yc, double *fx_vol, double disc_dt,
                         double fx_dt);

Err Fx3DBetaImpliedVol(char *fx_underlying, double beta, double val_time,
                       double start_time, double end_time, double disc_dt,
                       double fx_dt, double *vol);

Err Fx3DBetatsCalibration(long today, double *exercise_opt,
                          double *maturity_opt, double *vol_opt, long nbrOpt,
                          double *maturity, long nbrMat, double *sig_curve_dom,
                          double lda_dom, double *sig_curve_for, double lda_for,
                          double *beta, double spot_fx, double correl_dom_for,
                          double correl_dom_fx, double correl_for_fx,
                          char *dom_yc, char *for_yc, double **fx_vol_curve,
                          double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DBetaCalibration(char *dom_underlying, char *for_underlying,
                        double spot_fx, double beta, double correl_dom_for,
                        double correl_dom_fx, double correl_for_fx,
                        double *exercise_opt, double *maturity_opt,
                        double *vol_opt, long nbropt, double **fx_vol_curve,
                        double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DBetatsTreeFxOptions(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double corr_dom_for, double corr_dom_fx, double corr_for_fx, char *dom_yc,
    char *for_yc, double *option_prices, long num_stp);

Err Fx3DAlphaBetatsTreeFxOptions(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double corr_dom_for, double corr_dom_fx, double corr_for_fx, char *dom_yc,
    char *for_yc, double *option_prices, long num_stp);

Err Fx3DBetatsCalibration2(long today, double *exercise_opt,
                           double *maturity_opt, double *vol_opt, long nbrOpt,
                           long nbrLong, double *maturity, long nbrMat,
                           double *sig_curve_dom, double lda_dom,
                           double *sig_curve_for, double lda_for, double beta,
                           double spot_fx, double corr_dom_for,
                           double corr_dom_fx, double corr_for_fx, char *dom_yc,
                           char *for_yc, double **fx_vol_curve, long nbSteps,
                           long nbNewton, double disc_dt, double fx_dt,
                           long nbIterMax);

Err Fx3DBetaCalibration2(char *dom_underlying, char *for_underlying,
                         double spot_fx, double beta, double correl_dom_for,
                         double correl_dom_fx, double correl_for_fx,
                         double *exercise_opt, double *maturity_opt,
                         double *vol_opt, long nbropt, long nbLong,
                         double **fx_vol_curve, long nbSteps, long nbNewton,
                         double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DAlphaBetatsCalibration2(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for, double alpha, double beta, double spot_fx,
    double corr_dom_for, double corr_dom_fx, double corr_for_fx, char *dom_yc,
    char *for_yc, double **fx_vol_curve, long nbSteps, long nbNewton,
    double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DAlphaBetaCalibration2(char *dom_underlying, char *for_underlying,
                              double spot_fx, double alpha, double beta,
                              double correl_dom_for, double correl_dom_fx,
                              double correl_for_fx, double *exercise_opt,
                              double *maturity_opt, double *vol_opt,
                              long nbropt, long nbLong, double **fx_vol_curve,
                              long nbSteps, long nbNewton, double disc_dt,
                              double fx_dt, long nbIterMax);

Err Fxbeta_log_approx_corr(long today, double *maturity, long nb_mat,
                           double *sig_dom, double *sig_for, double *mat_fx,
                           long nb_mat_fx, double *sig_fx, double *beta,
                           double lam_dom, double lam_for, double *mat_corr,
                           double *corr_dom_for_ts, double *corr_dom_fx_ts,
                           double *corr_for_fx_ts, long nb_mat_corr,
                           double spot, char *dom_yc, char *for_yc,
                           double *time_fx, long nb_time_fx, double *fx_fwd,
                           double *fx_vol, double max_time);

/*	Calibration of a fx term structure to a set of fx options Corr TS
 * available */
Err Fx3DBetatsImpliedVol_corr(
    long today, double opt_maturity, double start_date, double end_date,
    double *maturity, long nbMat, double *sig_curve_dom, double lda_dom,
    double *sig_curve_for, double lda_for, double *mat_fx, long nb_mat_fx,
    double *sig_curve_fx, double *beta, double spot_fx, double *mat_corr,
    double *corr_dom_for_ts, double *corr_dom_fx_ts, double *corr_for_fx_ts,
    long nb_mat_corr, char *dom_yc, char *for_yc, double *fx_vol,
    double disc_dt, double fx_dt);

/*	Implied vol direct from underlying */
Err Fx3DBetaImpliedVol_corr(char *fx_underlying, double beta, double val_time,
                            double start_time, double end_time, double disc_dt,
                            double fx_dt, double *vol);

Err Fx3DBetatsCalibration_corr(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, double *maturity, long nbrMat, double *sig_curve_dom,
    double lda_dom, double *sig_curve_for, double lda_for, double *beta,
    double spot_fx, double *mat_corr, double *corr_dom_for_ts,
    double *corr_dom_fx_ts, double *corr_for_fx_ts, long nb_mat_corr,
    char *dom_yc, char *for_yc, double **fx_vol_curve, double disc_dt,
    double fx_dt, long nbIterMax);

Err Fx3DBetaCalibration_corr(char *dom_underlying, char *for_underlying,
                             double spot_fx, double beta, double *mat_corr,
                             double *corr_dom_for_ts, double *corr_dom_fx_ts,
                             double *corr_for_fx_ts, long nb_mat_corr,
                             double *exercise_opt, double *maturity_opt,
                             double *vol_opt, long nbropt,
                             double **fx_vol_curve, double disc_dt,
                             double fx_dt, long nbIterMax);

Err Fx3DBetatsCalibration2_corr(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for, double beta, double spot_fx, double *corr_mat,
    double *corr_dom_for_ts, double *corr_dom_fx_ts, double *corr_for_fx_ts,
    long nb_corr_mat, char *dom_yc, char *for_yc, double **fx_vol_curve,
    long nbSteps, long nbNewton, double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DBetatsTreeFxOptions_corr(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double *corr_mat, double *corr_dom_for_ts, double *corr_dom_fx_ts,
    double *corr_for_fx_ts, long nb_corr_mat, char *dom_yc, char *for_yc,
    double *option_prices, long num_stp);

Err Fx3DBetaCalibration2_corr(char *dom_underlying, char *for_underlying,
                              double spot_fx, double beta, double *correl_mat,
                              double *correl_dom_for, double *correl_dom_fx,
                              double *correl_for_fx, long nb_correl,
                              double *exercise_opt, double *maturity_opt,
                              double *vol_opt, long nbropt, long nbLong,
                              double **fx_vol_curve, long nbSteps,
                              long nbNewton, double disc_dt, double fx_dt,
                              long nbIterMax);

Err Fx3DAlphaBetatsTreeFxOptions_corr(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double *corr_mat, double *corr_dom_for, double *corr_dom_fx,
    double *corr_for_fx, long nb_corr, char *dom_yc, char *for_yc,
    double *option_prices, long num_stp);

Err Fx3DAlphaBetatsCalibration2_corr(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for, double alpha, double beta, double spot_fx, double *corr_mat,
    double *corr_dom_for, double *corr_dom_fx, double *corr_for_fx,
    long nb_corr, char *dom_yc, char *for_yc, double **fx_vol_curve,
    long nbSteps, long nbNewton, double disc_dt, double fx_dt, long nbIterMax);

Err Fx3DAlphaBetaCalibration2_corr(
    char *dom_underlying, char *for_underlying, double spot_fx, double alpha,
    double beta, double *corr_mat, double *correl_dom_for,
    double *correl_dom_fx, double *correl_for_fx, long nb_corr,
    double *exercise_opt, double *maturity_opt, double *vol_opt, long nbropt,
    long nbLong, double **fx_vol_curve, long nbSteps, long nbNewton,
    double disc_dt, double fx_dt, long nbIterMax);
#endif