
#ifndef FxSabrAdiH
#define FxSabrAdiH

#include "srt_h_all.h"

void disc_normal_center(double *x, long nstp, double fwd, double xmin,
                        double xmax, double stdn, double nstd, long *index);
void disc_linleft_logright(double *x, long nstp, double fwd, double xmin,
                           double xmax, double stdln, double nstd, long *index);
void disc_linleft_logright_center(double *x, long nstp, double fwd, double xmin,
                                  double xmax, double stdln, double nstd,
                                  long *index);
void disc_linleft_linright_center(double *x, long nstp, double fwd, double xmin,
                                  double xmax, double stdln, double nstd,
                                  long *index);
void disc_linleft_logright_bar(double **x, long *nstp, double fwd, double xmin,
                               double xmax, double stdln, double nstd,
                               double *barrier, long nb_bar, double *bar_pres,
                               long *index);
void disc_linleft_logright_bar2(double **x, long *nstp, double fwd, double xmin,
                                double xmax, double stdln, double nstd,
                                double *barrier, long nb_bar, double *bar_pres,
                                long *index);
void disc_linleft_logright_strike(double *x, long nstp, double fwd, double xmin,
                                  double xmax, double stdln, double nstd,
                                  double strike, long *index);
void disc_linleft_logright_center_strike(double *x, long nstp, double fwd,
                                         double xmin, double xmax, double stdln,
                                         double nstd, double strike,
                                         long *index);
void disc_SL_center(double *x, long nstp, double fwd, double xmin, double xmax,
                    double a, double b, double std, double nstd, long *index);
void disc_SL_center_linleft_logright(double *x, long nstp, double fwd,
                                     double xmin, double xmax, double a,
                                     double b, double std, double nstd,
                                     long *index);
void disc_SL_center_barrier(double **x, long *nstp, double fwd, double xmin,
                            double xmax, double a, double b, double std,
                            double nstd, double *bar, long nb_bar, long *index);
void disc_SL_strike(double *x, long nstp, double fwd, double xmin, double xmax,
                    double a, double b, double std, double nstd, double strike,
                    int is_center, long *index);
void Interpolate_Payoff(double ***values, double *x, int ux, int nstepx, int lz,
                        int uz, int nprod, int bar_col, double barrier,
                        int interp_bar);

Err FxSabr_adi(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floormu,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int *bar_col,
    int *is_bar,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       long today, double spot_fx, /*	The cash one */
                       void *dom_yc, void *for_yc,

                       /* Grid data	*/
                       int l1, int u1, int l2, int u2, double *x,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res);

Err FxSabrCalibration(
    /*	Underlying	*/
    char *dom_yc, char *for_yc, long today, double spot_fx,

    /*	Model Parameters			*/
    double alpha, double beta, double rho, double lambda,

    double floormu,

    /*	Options Parameters			*/
    double *exercise, double *maturity, double *volatility, int nbOpt,

    /*	Discretisation Parameters	*/
    long nstpt, int nstpfx, int nstpvol, int nbIter, double precision,

    /*	Result						*/
    double **sigma);

Err FxSabr_adi_bar(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floorstd,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int bar_col,
    int *is_bar, int *is_up,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       long today, double spot_fx, /*	The cash one */
                       void *dom_yc, void *for_yc,

                       /* Grid data	*/
                       int l1, int u1, int l2, int u2, double *x,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res);

Err FxSabr_KOOption(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floorstd,

    /*	Product data */
    double strike, int is_call, /* 1 Call  , 0: Put */
    double *bar_lvl, int is_up, /* 1 Up  , 0: Down */
    int is_cvx,                 /* 1 use 1 / Fx  , 0 use Fx */

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Result */
    double *res);

double solve_for_next_coef(double **res_iter, int nb_iter, double premium_tgt,
                           int method);

double find_beta_lim(double forward, double std_beta, double beta,
                     double percent);

Err op_sabrSL_adi(double forward, double *strike, int nb_strike,
                  double maturity, double disc,
                  int call_put_var, /*	0:	(F - K)+
                                                            1:	(K - F)+
                                                            2:	F^2	*/
                  double sigma_beta, double alpha, double beta, double rho,
                  double lambda, double floorstd, int nstp, int nstepx,
                  int nstepz, double *res,
                  /* Additional informations */
                  int calc_greeks,
                  double **greeks, /* array 6 * nprod containing delta  , gamma
                                      , theta  , vega  , volga and vanna */

                  /* For calibration purpose */
                  int calc_at_point, int column, double target, double *vol,
                  double *res_at_point);

Err op_sabrSL_calib_adi(double forward, double strike, double maturity,
                        double tgt_vol, SrtDiffusionType input_vol_type,
                        double alpha, double beta, double rho, double lambda,
                        int nt, int nx, int nz, int nbIter, double precision,
                        double floor_std, double *guess, double *res);

Err op_sabrSL_MC(double forward, double *strike, int nb_strike, double maturity,
                 double sigma_beta, double alpha, double beta, double rho,
                 double lambda, int npaths, int nsteps, double **res);

Err op_sabrSL_MC2(double forward, double *strike, int nb_strike,
                  double maturity, double sigma_beta, double alpha, double beta,
                  double rho, double lambda, int npaths, int nsteps,
                  int do_balsam, double **res);

double op_sabrSL(double F, double K, double T, double sigma, double alpha,
                 double beta, double rho, double lambda);

double op_sabrSLcalib(double F, double K, double T, double sigma, double alpha,
                      double beta, double rho, double lambda);

Err op_sabrQuad_MC(double forward, double *strike, int nb_strike,
                   double maturity, double sigma_beta, double alpha, double a,
                   double b, double c, double rho, double lambda, int npaths,
                   int nsteps, int do_balsam, double **res);

#endif