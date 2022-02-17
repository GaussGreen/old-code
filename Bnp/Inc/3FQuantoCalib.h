#ifndef __3F_QUANTO_CALIB_H
#define __3F_QUANTO_CALIB_H

// Define columns names for handling 4x4 and 5x5 matrices
#define W1_D 0
#define W2_D 1
#define W3_F 2
#define W4_F 3
#define FX_5F 4

#define LGM2F_W1 0
#define LGM2F_W2 1
#define LGM1F_W3 2
#define FX_3F_QUANTO 3

// Returns the number of factors of a LGM underlying
Err get_lgm_number_of_factors(char *underlying, int *number_of_factors);

// Filling of the new 5x5 correlation matrix
// 0=W1F  , 1=W2F  , 2=W1D  , 3=W2D  , 4=FX
Err convert_4x4_to_5x5_correlation_matrix(double **correl_matrix,
                                          double **correl_matrix_5x5,
                                          int is_lgm2F_dom_for);

// Takes a domestic and a foreign LGM underlying (one being 1F  , the other 2F)
// and returns 2 LGM2F structures
Err convert_3FQuanto_to_2_LGM2F(
    char *dom_und, char *for_und, double **dom_sigma_dates,
    long *n_dom_sigma_dates, double **dom_sigma_values,
    double *dom_fixed_lambda, double *dom_fixed_alpha, double *dom_fixed_gamma,
    double *dom_fixed_rho, double **for_sigma_dates, long *n_for_sigma_dates,
    double **for_sigma_values, double *for_fixed_lambda,
    double *for_fixed_alpha, double *for_fixed_gamma, double *for_fixed_rho,
    int *lgm2F_is_dom_for);

// Takes a 3F Quanto underlying and turns it into one LGM2F and one LGM1F
Err Get_FX_StochRate_TermStructures3FQuanto(
    char *underlying, char **dom_und_name, double **dom_sigma_dates,
    long *n_dom_sigma_dates, double **dom_sigma_values,
    double *dom_fixed_lambda, double *dom_fixed_alpha, double *dom_fixed_gamma,
    double *dom_fixed_rho, char **for_und_name, double **for_sigma_dates,
    long *n_for_sigma_dates, double **for_sigma_values,
    double *for_fixed_lambda, double *for_fixed_alpha, double *for_fixed_gamma,
    double *for_fixed_rho, int *lgm2F_is_dom_for, double **fx_sigma_dates,
    long *n_fx_sigma_dates, double **fx_sigma_values);

// Calculates the coefficients of the quadratic equation
// on an interval [t1  , t2[ where ALL parameters (LGM1F  , LGM2F  , FX sigmas
// , correls and lambdas) are CONSTANT DOMESTIC DYNAMICS: sigma1  , lambda1 +
// sigma2  , lambda2 FOREIGN DYNAMICS: sigma3  , lambda3 + sigma4  , lambda4
// Reminder: in the double** 5x5 correl_matrix  , the indices are:
// 0: DOM W1  , 1: DOM W2  , 3: FOR W3  , 4: FOR W4  , 5: FX
Err Coefs_Partial_Var_5F(double Tk, double t1, double t2, double sigma1,
                         double lambda1, double sigma2, double lambda2,
                         double sigma3, double lambda3, double sigma4,
                         double lambda4, double **correl_matrix, double *a,
                         double *b, double *c);

// Calculates the partial volatility of a FX fwd maturing at time
// fx_option_maturity on an interval [t1  , t2[ where ALL parameters (LGM1F  ,
// LGM2F  , FX sigmas  , correls and lambdas) are CONSTANT
Err Partial_Var_5F(double fx_option_maturity, double t1, double t2,
                   double dom_sigma1, double dom_lambda1, double dom_sigma2,
                   double dom_lambda2, double for_sigma3, double for_lambda3,
                   double for_sigma4, double for_lambda4, double fx_sigma,
                   double **correl_matrix, double *var_partial);

// Gets the implied FX fwd vol from a 3F Quanto underlying
Err Fx3FQuantoImpliedVolFrom3FQuantoUnderlying(
    char *und_3F_quanto, double *correl_matrix_dates,
    double ***correl_matrix_ts, long n_correl_matrix_ts, double val_time,
    double start_time, double fix_time, double *fx_vol);

// FX calibration directly from the underlying
// given a TERM STRUCTURE OF CORRELATIONS
Err Fx3FQuantoCalibrationCorr(
    char *dom_und, char *for_und, double *correl_matrix_dates,
    double ***correl_matrix_ts, long n_correl_matrix_ts,
    double *fx_options_exo_times, double *fx_options_mat_times,
    double *fx_options_vols, long n_fx_options, double **fx_vol_curve);

// Calibration of a FX term structure to a set of European FX options
// given a TERM STRUCTURE OF CORRELATIONS 5x5
// We know the partial vol up to T1 and we bootstrap on the spot fx volatility
// between [T1  , T2] using the FX option maturing at T2
Err Fx5FtsCalibrationCorr(
    double *fx_options_exo_times, double *fx_options_mat_times,
    double *fx_options_vols, long n_fx_options, double *ir_dates,
    long n_ir_dates, double *dom_sigma_values, double dom_fixed_lambda,
    double dom_fixed_alpha, double dom_fixed_gamma, double *for_sigma_values,
    double for_fixed_lambda, double for_fixed_alpha, double for_fixed_gamma,
    double *correl_matrix_dates, long n_correl_matrix_dates,
    double ***correl_matrix_ts_5x5, double **fx_vol_curve);

// Calculates the partial implied vol between start_date and end_date of
// the fwd FX of time fx_option_maturity given the term structures of 2 LGM2Fs
// and a TERM STRUCTURE OF CORRELATIONS
// The TS of rates must be merged before and the correl matrix MUST be 5x5
Err Fx5FImpliedVolCorr(double fx_option_maturity, double start_date,
                       double end_date, double *ir_dates, long n_ir_dates,
                       double *dom_sigma_values, double dom_fixed_lambda,
                       double dom_fixed_alpha, double dom_fixed_gamma,
                       double *for_sigma_values, double for_fixed_lambda,
                       double for_fixed_alpha, double for_fixed_gamma,
                       double *fx_sigma_times, double *fx_sigma_values,
                       long n_fx_sigma_times, double *correl_matrix_dates,
                       long n_correl_matrix_dates,
                       double ***correl_matrix_5x5_ts, double *fx_vol);

// FX calibration directly from the underlying 5F (2 LGM2F)
// given a TERM STRUCTURE OF CORRELATIONS
Err Fx5FCalibrationCorr(char *dom_und, char *for_und,
                        double *correl_matrix_dates, double ***correl_matrix_ts,
                        long n_correl_matrix_ts, double *fx_options_exo_times,
                        double *fx_options_mat_times, double *fx_options_vols,
                        long n_fx_options, double **fx_vol_curve);

// Gets the implied FX fwd vol from a 5F underlying
// given a TERM STRUCTURE OF CORRELATIONS
Err Fx5FImpliedVolFrom5FUnderlying(char *und_5F, double *correl_matrix_dates,
                                   double ***correl_matrix_ts,
                                   long n_correl_matrix_ts, double val_time,
                                   double start_time, double fix_time,
                                   double *fx_vol);

#endif