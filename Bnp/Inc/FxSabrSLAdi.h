
#ifndef FxSabrSLAdiH
#define FxSabrSLAdiH

#include "srt_h_all.h"

double solve_var_diff(double a, double tgt);

Err FxSabrSLPrecalculations(
    /*	Time data		*/
    int nstp, double *time,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double a, double b, double rho,
    double lambda,

    /*	Market data */
    double spot_fx,

    /*	Output			*/
    double *expect_z, double *drift_z, double *std1);

Err FxSabrSL_adi(
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
    int nprod, double *res,

    /* Additional informations */
    int calc_greeks,
    double **greeks, /* array 6 * nprod containing delta  , gamma  , theta  ,
                        vega  , volga and vanna */

    /* For calibration purpose */
    int calc_at_point, int column, double target, double *vol,
    double *res_at_point);

Err FxSabrSLCalibration(
    /*	Underlying	*/
    char *dom_yc, char *for_yc, long today, double spot_fx,

    /*	Model Parameters			*/
    double alpha, double beta, double rho, double lambda,

    double floormu,

    /*	Options Parameters			*/
    double *exercise, double *maturity, double *volatility, int nbOpt,
    int start_opt, int is_straddle,

    /*	Discretisation Parameters	*/
    long nstpt, int nstpfx, int nstpvol, int nbIter, double precision,

    /*	Guess						*/
    double *guess,

    /*	Result						*/
    double *sigma);

Err FxSabrSLCalibSmile(
    /*	Underlying						*/
    char *dom_yc, char *for_yc, long today, double spot_fx,

    /*	Model Parameters				*/
    double *alpha_out, double beta, double *rho_out, double lambda,
    double floormu,

    /*	User input of starting points	*/
    int calib_smile, int use_input, int use_total_ts,

    /*	Options Parameters				*/
    double *exercise, double *maturity, double *volatility, int nbOpt,
    int is_straddle,

    long mat_date, double risk_reversal, double butterfly,

    /*	Discretisation Parameters		*/
    long nstpt_tgt, int nstpfx, int nstpvol, int nbIterATM, double precisionATM,
    int nbIterSmile, double precisionSmile,

    /*	Result							*/
    double *strike1_out, double *vol1_out, double *strike2_out,
    double *vol2_out,

    double *sigma);

Err FxSabrSL_KOOption(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floorstd,

    /*	Product data */
    double strike, int is_call,  /* 1 Call  , 0: Put */
    int is_american, int is_cvx, /* 1 use 1/Fx  , 0: use Fx */
    int is_digital, /* 1: digital payoff  , 0  , regular option payoff */
    double *bar_lvl_up, double *bar_lvl_down, double rebate_up,
    double rebate_down,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Result */
    double *res,

    /* Additional informations */
    int calc_greeks,
    double *greeks); /* array 6 * nprod containing delta  , gamma  , theta  ,
                        vega  , volga and vanna */

#endif