#ifndef LGM2FpdeH
#define LGM2FpdeH

#include "srt_h_all.h"

/*	Main function without Tau TS */
/*	---------------------------- */
Err lgm2f_adi(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nsteps,

    /*	Model data		*/
    double lam1, double *sig_time, double *sig1, int nb_sig, double alpha,
    double gamma, double rho,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,

                       /* Model data	*/
                       double *lam, double *ts_time, int nb_ts, double gamma,
                       double rho, double phi1, double phi2, double phi12,

                       /* Gride data	*/
                       int l1, int u1, int l2, int u2, double *x, double **y,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res);

/*	Main function with Tau TS */
/*	------------------------- */
Err lgm2fTau_adi(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int nsteps, int disc_method, /* 0: linear distretisation      , 1: normal
                                    discretisation */

    /*	Model data		*/
    double *sig, double *sig_time, int nb_sig, double *lam, double *lam_time,
    int nb_lam, double alpha, double gamma, double rho,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,

                       /* Model data	*/
                       double *lam, double *ts_time, int nb_ts, double gamma,
                       double rho, double phi1, double phi2, double phi12,

                       /* Gride data	*/
                       int l1, int u1, int l2, int u2, double *r1_dim2,
                       double **r2,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res);

/*	Other function with Tau TS but specially adapted for neg Tau */
/*	------ ----------------------------------------------------- */
Err lgm2fTau2_adi(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int nsteps, int disc_method, /* 0: linear distretisation      , 1: normal
                                    discretisation */

    /*	Model data		*/
    double *sig, double *sig_time, int nb_sig, double *lam, double *lam_time,
    int nb_lam, double alpha, double gamma, double rho,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,

                       /* Model data	*/
                       double *lam, double *ts_time, int nb_ts, double gamma,
                       double rho, double phi1, double phi2, double phi12,

                       /* Gride data	*/
                       int l1, int u1, int l2, int u2, double *r1_dim2,
                       double **r2,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res);

#endif