#ifndef __3F_QUANTO_UTILS_H
#define __3F_QUANTO_UTILS_H

// MAXNODE3D is equal to 400 but since the 3F Quanto model needs more nodes
// we have to increase this limit
#define MAXNODE3FQUANTO 500

// Selects the dates in r1_dates  , r3_dates  , fx_dates which are <= last_time
// and put them in a single array. (r1 then r3 then fx)
Err compute_vol_times_3FQuanto(double *r1_dates, double *r1_sigma,
                               int n_r1_dates, double *r3_dates,
                               double *r3_sigma, int n_r3_dates,
                               double *fx_dates, double *fx_sigma,
                               int n_fx_dates, int *num_vol_times,
                               double **vol_times, double last_time);

// Merge rates  , fx and corr term structures and returns all variables in lists
Err merge_rates_fx_corr_ts_3FQuanto(
    double *sig_r1_mat, double *sig_r1, int sig_n_r1, double *sig_r3_mat,
    double *sig_r3, int sig_n_r3, double *sig_fx_mat, double *sig_fx,
    int sig_n_fx, double *corr_4x4_mat, double ***corr_4x4, int corr_4x4_n_mat,

    // Outputs
    double **sig_merge_mat, double **sig_merge_r1, double **sig_merge_r3,
    double **sig_merge_fx, double **corr_merge_r1_r2, double **corr_merge_r1_r3,
    double **corr_merge_r1_fx, double **corr_merge_r2_r3,
    double **corr_merge_r2_fx, double **corr_merge_r3_fx, int *sig_merge_n);

// Fills the expectations (under Q_Beta_Domestic) and variances for a 3FQuanto
// model
void fill_fwd_var_corr_3FQuanto(
    long nstp, double *time, double *date, double *sigma_R1, double *sigma_R2,
    double *sigma_R3, double *sigma_FX, double lambda_R1, double lambda_R2,
    double lambda_R3, double *correl_R1_R2, double *correl_R1_R3,
    double *correl_R1_FX, double *correl_R2_R3, double *correl_R2_FX,
    double *correl_R3_FX, char *lgm2F_yc, char *lgm1F_yc, int lgm2F_is_dom_for,
    /*	To be allocated by caller  , filled by the function */
    double *fwd_R1, double *var_R1, double *fwd_R2, double *var_R2,
    double *covar_R1_R2, double *fwd_R3, double *var_R3, double *lgm2F_ifr,
    double *lgm1F_ifr);

#endif
