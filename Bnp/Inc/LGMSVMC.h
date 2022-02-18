#ifndef LGMSVMCH
#define LGMSVMCH

#include "MCEBOptimisation.h"
#include "srt_h_all.h"
#include "LGMSVPDE.h"

Err lgmSV_mc_balsam(
    /*	Time Information  */
    int     iNbTime,
    int     iNbEvent,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double  dLambdaX,
    double* dSigma,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,

    /* Parameters for DF(t,T*) reconstruction */
    double* dff_star,
    double* gam_star,
    double* gam2_star,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(
        long   path_index,
        double evt_date,
        double evt_time,
        void*  func_parm,

        double ft,
        double psi,
        double v,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val,
        int*    stop_path),
    /*	Result */
    int      iNbProduct,
    double** res);

Err lgmSV2F_mc_balsam(
    /*	Time Information  */
    int     iNbTime,
    int     iNbEvent,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double dLambdaX1,
    double dLambdaX2,

    double* dSigma,
    double* dAlphaLGM,
    double* dRhoLGM,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,
    double* dRho2,

    /* Parameters for DF(t,T*) reconstruction */
    double* dff_star,
    double* gam1_star,
    double* gam2_star,
    double* gam1_2_star,
    double* gam2_2_star,
    double* gam12_star,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(
        long   path_index,
        double evt_date,
        double evt_time,
        void*  func_parm,

        /* Model data	*/
        double ft1,
        double ft2,
        double phi1,
        double phi2,
        double phi12,
        double v,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val,
        int*    stop_path),
    /*	Result */
    int      iNbProduct,
    double** res);

Err lgmSV_mc_balsam_rev(
    /*	Time Information  */
    int     iNbTime,
    int     iNbEvent,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double  dLambdaX,
    double* dSigma,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,

    /* Parameters for DF(t,T*) reconstruction */
    double* dff_star,
    double* gam_star,
    double* gam2_star,

    /* Parameters */
    LGMSVPARAM Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(
        long   path_index,
        double evt_date,
        double evt_time,
        void*  func_parm,

        double ft,
        double psi,
        double v,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val,
        int*    stop_path),

    Err (*payoff_adjust_function)(
        long       event_index,
        long       npaths,
        long       nprod,
        double***  saved_values,
        void*      func_parm,
        MCEBPARAMS mcebparams),

    /*	Result */
    int      iNbProduct,
    double** res);

Err lgmSV2F_mc_balsam_rev(
    /*	Time Information  */
    int     iNbTime,
    int     iNbEvent,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double dLambdaX1,
    double dLambdaX2,

    double* dSigma,
    double* dAlphaLGM,
    double* dRhoLGM,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,
    double* dRho2,

    /* Parameters for DF(t,T*) reconstruction */
    double* dff_star,
    double* gam1_star,
    double* gam2_star,
    double* gam1_2_star,
    double* gam2_2_star,
    double* gam12_star,

    /* Parameters */
    LGMSVPARAM Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(
        long   path_index,
        double evt_date,
        double evt_time,
        void*  func_parm,

        /* Model data	*/
        double ft1,
        double ft2,
        double phi1,
        double phi2,
        double phi12,
        double v,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val,
        int*    stop_path),

    Err (*payoff_adjust_function)(
        long       event_index,
        long       npaths,
        long       nprod,
        double***  saved_values,
        void*      func_parm,
        MCEBPARAMS mcebparams),

    /*	Result */
    int      iNbProduct,
    double** res);

Err lgmSV2F_mc_balsam_optim_mem(
    /*	Time Information  */
    int     iNbTime,
    int     iNbEvent,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double dLambdaX1,
    double dLambdaX2,

    double* dSigma,
    double* dAlphaLGM,
    double* dRhoLGM,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,
    double* dRho2,

    /* Parameters for DF(t,T*) reconstruction */
    double* dff_star,
    double* gam1_star,
    double* gam2_star,
    double* gam1_2_star,
    double* gam2_2_star,
    double* gam12_star,

    /* Parameters */
    LGMSVPARAM Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(
        long   path_index,
        double evt_date,
        double evt_time,
        void*  func_parm,

        /* Model data	*/
        double ft1,
        double ft2,
        double phi1,
        double phi2,
        double phi12,
        double v,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val,
        int*    stop_path),
    /*	Result */
    int      iNbProduct,
    double** res);

Err lgmSV_mc(
    /*	Time Information  */
    int     iNbTime,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double  dLambdaX,
    double* Sig,
    double* SigPsi,
    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,
    double* dRho2,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /*	Payoff function */
    Err (*payoff_func)(
        double evt_date,
        double evt_time,
        void*  func_parm,

        double ft,
        double psi,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val),
    /*	Result */
    int      iNbProduct,
    double** res);

Err lgmSV_mc_cv(
    /*	Time Information  */
    int     iNbTime,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double  dLambdaX,
    double* Sig,
    double* SigPsi,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,
    double* dRho2,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /*	Payoff function */
    Err (*payoff_func)(
        double evt_date,
        double evt_time,
        void*  func_parm,

        double ft,
        double psi,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val),
    /*	Result */
    int      iNbProduct,
    double** res);

Err LGMSVOptionMC(
    char* und_name,      /*	Name of the underlying */
    char* yc_name,       /*	Name of the yield curve */
    char* ref_rate_name, /*	Name of the reference rate */
    char* swaption_freq, /*	Frequency and basis of underlying swaptions */
    char* swaption_basis,

    long    lExDate,
    long    lEndDate,
    double* dStrike,
    int     nb_strike,
    int     pay_rec, /*	pay:1 rec:-1 */

    long nb_paths,
    long nb_vol,
    int  use_balsam,

    /* Output */
    double* pSwaptionPrice,
    double* pStd);

Err LGMSVOptionMCts(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double  dLambdaX,
    int     iNbSigTime, /* Term Structure of g(t) */
    double* SigTime,
    double* Sig,
    double  dTStar, /* Tstar in years from today */
    double  dAlpha,
    double  dLambdaEps,
    double  dRho,

    /* Product description */
    long     lExDate,   /* Exercice date of the swaption  */
    double   dExTime,   /* Exercice of the swaption in years from today */
    int      iNbCoupon, /* Description of the cashflows */
    double*  CouponTime,
    long*    CouponDate,
    double** Coupon,
    int      nb_strike,
    char*    cYieldCurve, /* Yield Curve */

    /* Parameter of grids */
    long nb_paths,
    long nb_vol,

    /* Outputs */
    double* Price,
    double* pStd);

Err lgmSV2F_mc(
    /*	Time Information  */
    int     iNbTime,
    double* dTime,
    double* dDate,

    int iNumPaths,

    /*	Model data Information	*/
    double  dLambdaX1,
    double  dLambdaX2,
    double* Sig1,
    double* Sig2,
    double* SigPsi1,
    double* SigPsi2,
    double* SigPsi12,

    double* dAlpha,
    double* dLambdaEps,
    double* dLvlEps,
    double* dRho,
    double* dRho2,

    double* Coef_Ortho21,
    double* Coef_Ortho22,
    double* Coef_Ortho31,
    double* Coef_Ortho32,
    double* Coef_Ortho33,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void** func_parm_tab,
    int*   EvalEvent,

    /*	Payoff function */
    Err (*payoff_func)(
        double evt_date,
        double evt_time,
        void*  func_parm,

        double ft1,
        double ft2,
        double psi1,
        double psi2,
        double psi12,

        /* Vector of results to be updated */
        int nprod,
        /* Result	*/
        double* prod_val),
    /*	Result */
    int      iNbProduct,
    double** res);

/* Put points between original points to ensure time between two points does not exceed Max Time */
Err fill_time_vector_max_time(
    int iInitNbTimes, double* dInitTimes, double dMaxTime, int* iNewNbTimes, double** dNewTimes);

#endif