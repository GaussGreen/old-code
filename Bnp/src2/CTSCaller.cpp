#include "CTSCaller.h"

#include "CTSReserve.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srt_h_all_LGMSV.h"
#include "srtaccess.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define CTS_MINALPHA 0.005

/* Set the default values for all the CTS autocal parameters */
Err cts_set_default_params(
    int*    accrue_on_barrier,
    int*    value_zero,
    int*    use_cmsopt,
    double* correl_start,
    double* correl_end,
    int*    float_adjust_type,
    int*    iv_call_sell,
    double* iv_call_spread,
    int*    call_call_sell,
    double* call_call_spread,
    double* numer_call_spread,

    int*    iv_trim_type,
    int*    iv_max_fix,
    double* iv_min_fix_time,
    int*    call_trim_type,
    int*    call_max_fix,
    double* call_min_fix_time,

    int*    iNbX,
    double* iNbSigmaXLeft,
    double* iNbSigmaXRight,
    double* dIntegParam,
    int*    iIntegMethod,
    double* dVolLimit,
    int*    iCalibLGM,
    double* dMinStd,
    double* dMaxStd,
    double* numer_tstar,
    double* precision,
    int*    nbIterMax,
    int*    keep_first,

    int*    pde_or_mc,
    int*    nstpt,
    int*    nstpx,
    int*    nstpvol,
    int*    nstpphi,
    double* mc_mintime,
    long*   npaths,
    double* integ_mintime,

    int*                 nb_factor,
    double*              lgm_alpha,
    double*              lgm_gamma,
    double*              lgm_rho,
    int*                 calib_strategy,
    int*                 fix_lambda,
    int*                 fix_smile,
    int*                 smile_calib_months,
    LGMSV_CalibParams*   lgmsv_calib_params,
    int*                 force_atm,
    int*                 strike_flag_long,
    int*                 strike_flag_short,
    double*              max_std_long,
    double*              max_std_short,
    double*              vol_shift_long,
    DIAGCALIB_VOLTYPE*   vol_type_long,
    DIAGCALIB_SHIFTTYPE* vol_shift_type_long,
    double*              vol_shift_short,
    DIAGCALIB_VOLTYPE*   vol_type_short,
    DIAGCALIB_SHIFTTYPE* vol_shift_type_short,
    double*              lambda_shift,
    double*              min_time,
    int*                 skip_last,
    double*              min_fact,
    double*              max_fact,
    int*                 use_jumps,
    int                  ref_months,
    char*                short_tenor,

    int* calc_fwdiv,
    int* adjust_fee,

    int*    do_one_time,
    int*    one_time_index,
    int*    compute_reserve,
    int*    reserve_method,
    int*    lgm_reserve,
    int*    lgm_nstpt,
    int*    lgm_nstpx,
    int*    midat_reserve,
    double* one_time_vega,
    double* lambda_reserve,
    int*    recalib_european,
    int*    recalc_one_factor,
    int*    euro_nb_iter,

    double* tstar,

    int* compatibility_flag)
{
    char def_short_tenor[256];

    LGMSV_SetDefault_CalibParams(lgmsv_calib_params);

    lgmsv_calib_params->use_sabr_calib      = 1;
    lgmsv_calib_params->sabr_calib_min_time = 0.0;

    *accrue_on_barrier = 1;
    *value_zero        = 0;
    *use_cmsopt        = 0;
    *correl_start      = 1.0;
    *correl_end        = 1.0;
    *float_adjust_type = 3;
    *iv_call_sell      = -1;
    *iv_call_spread    = 10.0 / 10000.0;
    *call_call_sell    = -1;
    *call_call_spread  = 10.0 / 10000.0;
    *numer_call_spread = 35.0 / 10000.0;
    *iv_trim_type      = 2;
    *iv_max_fix        = 100000;
    *iv_min_fix_time   = 1.0 / 16.0;
    *call_trim_type    = 1;
    *call_max_fix      = 2;
    *call_min_fix_time = 0.25;

    *numer_tstar = 0;

    lgmsv_app_set_default_params(
        iNbX,
        iNbSigmaXLeft,
        iNbSigmaXRight,
        dIntegParam,
        iIntegMethod,
        dVolLimit,
        iCalibLGM,
        dMinStd,
        dMaxStd);

    *precision  = 0.00001;
    *nbIterMax  = 10;
    *keep_first = 0;

    *pde_or_mc     = 0;
    *nstpt         = 50;
    *nstpx         = 100;
    *nstpvol       = 30;
    *nstpphi       = 45;
    *mc_mintime    = 1.0 / 24.0;
    *npaths        = 20000;
    *integ_mintime = 0.25;

    *nb_factor = 1;
    *lgm_alpha = 1.2;
    *lgm_gamma = 0.2;
    *lgm_rho   = -0.85;

    *calib_strategy     = 1;
    *fix_lambda         = 1;
    *fix_smile          = 1;
    *smile_calib_months = -1;
    *force_atm          = 0;
    *strike_flag_long   = 0;
    *strike_flag_short  = 0;
    *max_std_long       = 999;
    *max_std_short      = 999;

    *vol_shift_long       = 0.0;
    *vol_type_long        = LGM_VOL;
    *vol_shift_type_long  = MULTIPLICATIVE;
    *vol_shift_short      = 0.0;
    *vol_type_short       = LGM_VOL;
    *vol_shift_type_short = MULTIPLICATIVE;
    *lambda_shift         = 0.0;

    *min_time  = 0.95;
    *skip_last = 1;
    *min_fact  = 0.5 * 0.5;
    *max_fact  = 0.5 * 0.5 / 4.0;
    *use_jumps = 0;

    sprintf(def_short_tenor, "%dM", ref_months);
    strcpy(short_tenor, def_short_tenor);

    *calc_fwdiv = 1;
    *adjust_fee = 1;

    *do_one_time       = 0;
    *one_time_index    = -1000;
    *compute_reserve   = 0;
    *one_time_vega     = -1.0 / 100.0;
    *reserve_method    = 0;
    *lgm_reserve       = 1;
    *lgm_nstpt         = 100;
    *lgm_nstpx         = 200;
    *midat_reserve     = 1;
    *lambda_reserve    = 0.0;
    *recalib_european  = 1;
    *recalc_one_factor = 0;
    *euro_nb_iter      = 5;

    *tstar = LGMSV_Tstar;

    *compatibility_flag = 0;
    return NULL;
}

Err cts_caller(
    /*	Market */
    long  today,
    char* yc,           /*	yc */
    char* vc,           /*	vc */
    char* ref,          /*	ref rate (only if calib) */
    char* swap_freq,    /*	swap freq (only if calib) */
    char* swap_basis,   /*	swap basis (only if calib) */
    char* lib_freq,     /*	libor freq (only if calib) */
    char* lib_basis,    /*	libor basis (only if calib) */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    /*	The underlying */
    double tstar,
    int    use_calib, /*	0: use lgmsvund, 1: calibrate */
    /*		if calib */
    int     nb_factor,
    double  lambda, /*	LGM lambda */
    double  lgm_alpha,
    double  lgm_gamma,
    double  lgm_rho,
    int     nsmilepar,
    double* smilepartime,
    double* alphaeps, /*	alpha */
    double* rhoeps,   /*	rho */
    double* ldaeps,   /*	ldaeps */
    double* rho2eps,  /*	rho2eps */
    /*	End of calib params */
    char* lgmsvund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    int     fund_ccy, /*	0: domestic, 1: other */
    double* fund_not_ts,
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double*
        fund_fix_cpn, /*	Past coupon fixing if relevant,
                                                      includes spr, but not mrg, cvg and notional */
    double* cpn_not_ts,
    /*			number of coupons */
    int ncpn,
    /*			coupon description */
    /*				cpn */
    int* cpn_type, /*	0:	fixed
                                   1:	libor fixed at start
                                   2:	libor fixed at end */
    long*   cpn_start_date,
    long*   cpn_end_date,
    long*   cpn_pay_date,
    double* cpn_coupon,
    char**  cpn_basis,
    /*				pay libor */
    long*   pay_fix_date, /*	fixing date of the Libor */
    int     pay_months,
    char*   pay_freq,
    char*   pay_basis,
    double* pay_fwd_spread,
    double* pay_gearing,
    /*				fix libor */
    /*			general Libor properties */
    int fix_lag_bd,
    /*				fixings */
    int*     nfix,
    double** weights,
    long**   fix_dates,
    char**   fix_tenor,
    char*    fix_freq,
    char*    fix_basis,
    double** fix_fwd_spreads,
    /*				profiles */
    int tstype, /*	0: generic, 1: range */
    /*				profile 0 */
    double*  alpha,
    double*  beta,
    int*     nstr,
    double** str,
    double** nbopt,
    /*				profile 1 */
    int     iv_buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int     call_buy_sell, /*	1: BNPP buys, -1: BNPP sells */
    double  iv_call_spread,
    double  call_call_spread,
    double  numer_call_spread,
    int     accrue_on_barrier,
    int     value_zero, /*	0: Do not value low strike options, 1: do */
    double* lb,         /*	Lower bounds of ranges */
    double* ub,         /*	Upper bounds of ranges */
    double* payoff,     /*	Payoff in ranges */

    /*				fixings	*/
    double** fixed_ref, /*	Past coupon reference fixing if relevant */
    double*  fixed_pay, /*	Past coupon payment libor fixing if relevant */

    /*	Calls */
    int     ncall,
    int     pay_rec, /*	1: rec pd, -1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*		Trimming */
    int iv_trim_type, /*	0: no trim
                                      1: x fixings max
                                      2: x time min between two fixings */
    int    iv_max_fix,
    double iv_min_fix_time,
    int    call_trim_type, /*	0: no trim
                                           1: x fixings max
                                           2: x time min between two fixings */
    int    call_max_fix,
    double call_min_fix_time,
    /*	Extra model parameters*/
    int    use_cmsopt,        /*	Use CmsOption to value fix option, use BS on CmsRate otherwise */
    double correl_start,      /*	Correl between libor fixing and libor paid at start */
    double correl_end,        /*	Correl between libor fixing and libor paid at start */
    int    float_adjust_type, /*	type of adjustment for the floating coupon,
                                                      0: ATM vol, 1: Strike Vol */
    /*	Numerical params */
    /*		CF */
    int    iNbX,
    double iNbSigmaXGridLeft,
    double iNbSigmaXGridRight,
    double dIntegParam,
    int    iIntegMethod,
    double dVolLimit,
    int    iCalibLGM,
    double dMinStd,
    double dMaxStd,
    double numer_tstar,

    /*		ADI */
    int    pde_or_mc,
    int    req_stp,
    int    req_stppsi,
    int    req_stpx,
    int    req_stpz,
    double req_mintime,
    long   req_paths,
    double integ_mintime,

    /*		Calib */
    char*               cal_tenor,
    char*               cal_ref,
    char*               cal_freq,
    char*               cal_basis,
    int                 force_atm, /*	force atm calib */
    double              max_std_long,
    double              max_std_short,
    double              vol_shift_long,
    DIAGCALIB_VOLTYPE   vol_type_long,
    DIAGCALIB_SHIFTTYPE vol_shift_type_long,
    double              vol_shift_short,
    DIAGCALIB_VOLTYPE   vol_type_short,
    DIAGCALIB_SHIFTTYPE vol_shift_type_short,
    double              lambda_shift,
    int                 calib_strategy, /*	-1: autocal, 0: swaptions / cap, 1: cap / swaptions */
    int                 fix_lambda,     /*	0: calib lambda to cap, 1: fix lambda calib
                                                        to diagonal */
    char* short_tenor,
    char* short_refrate,
    char* short_freq,
    char* short_basis,
    int   fix_smile,          /*	0: calib smile parameters to market smile */
    int   smile_calib_months, /* 0: co-terminal swaption, otherwise underlyings with required nb
                                 months */
    LGMSV_CalibParams* lgmsv_calib_params,
    double             min_time,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps, /*	1: we allow jumps on vol, 0: we don't */
    double prec,
    int    maxiter,
    int    keep_first,
    /*	Strike choice */
    int long_strike_flag,  /*	0: ATM
                                                   1: Coupon
                                                   2: Eq (PV/Lvl) */
    int short_strike_flag, /*	0: ATM,
                                                   1: implied digital caplet strike
                                                   2: same number of std */
    /*		IV calculation */
    int calc_fwd_iv,
    int adj_fee,

    /*	Flag for extra calculation / adjustment of one-time callable */
    int    do_one_time,       /*	1: calc the one time */
    int    one_time_index,    /*	0: choose automatically the index, >0: index provided by user */
    int    compute_reserve,   /*	1: adjust the price from the reserve */
    int    reserve_method,    /*	0: old method, 1: new method */
    int    lgm_reserve,       /*	1: calculate the reserve in the one factor */
    int    lgm_nstpt,         /*	number of time steps for LGM midat */
    int    lgm_nstpx,         /*	number of space steps for LGM midat */
    int    midat_reserve,     /*	1: reserve calculated using the midat-replication */
    double one_time_vega,     /*	Vega to be applied one the one time callable if required */
    double lambda_reserve,    /*	Percent of the switch to be reserved */
    int    recalib_european,  /*	If set to 1, then european are recalibrated when changing lambda */
    int    recalc_one_factor, /*	Recalculate MCEB in one factor for reserve */
    int    euro_nb_iter,      /*	Number of Levenberg iteration for calibration of european */

    /*	vol Matrix parameters */
    int              num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double*          strikes_in_vol,
    SrtDiffusionType vol_type, /*	Type of vol in matrix, SRT_NORMAL or SRT_LOGNORMAL */
    int              cash_vol, /*	1: matrix is a cash vol
                                               0: matrix is a swap vol */
    /*	Flags */

    /*		EOD */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    /* compatibility flag for old version of Callable TS */
    int compatibility_flag,

    /*	Results */
    double* fund_val,        /*	Value of the funding leg */
    double* cpn_val,         /*	Value of the Power Dual leg */
    double* call_val,        /*	Value of the callable feature */
    double* onetime_val,     /*	Value of the one-time callable feature */
    double* mostexp_reserve, /*	One Time Callable Reserve */
    double* switch_reserve,  /*	Switch reserve */
    /*	Feedback */
    int                  export_und,
    CTS_UND              und_exp,
    int                  save_inst_data,
    cpd_calib_inst_data* inst_data,
    int                  save_fwdiv,
    CTS_IV               fwd_iv_info,
    int                  save_extra_infos,
    CTS_EXTRA_INFOS      extra_infos)
{
    Err    err = NULL;
    int    i, j, j0;
    int    free_struct, free_struct_iv, free_struct_iv_adj;
    int    for_fund;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;
    int    call_feat, call_feat_iv, call_feat_iv_adj;
    double fund_leg_pv, exo_leg_pv, call, onetime[100], fwd_iv;

    double       cpn_pv, cpn_coupon_temp;
    int          cpntype_temp;
    SrtBasisCode bas;
    cts *        cts_pricing = NULL, *cts_on_iv = NULL, *cts_iv_adj = NULL;
    CTS          cts_for_iv;
    cts_und *    und = NULL, *und_iv = NULL, *und_iv_adj = NULL;
    cts_adi_arg* adi_arg = NULL;

    cts_iv* fwd_iv_secu    = NULL;
    cts_iv* fwd_iv_reserve = NULL;

    double cpn_not, fund_not;

    cts_cpn* cpn2 = NULL;
    cts_mkt* mkt  = NULL;

    LGMSV_NumerParams* NumerParams = NULL;

    double  OtherCcyfundNot;
    char    Copyfund_ccy_yc[255];
    double *fund_mrg2 = NULL, *fund_spr2 = NULL, *fund_not_ts2 = NULL, *fund_fix_cpn2 = NULL;

    /* All Structure initialisations */
    cts_pricing = (cts*)calloc(1, sizeof(cts));
    cts_on_iv   = (cts*)calloc(1, sizeof(cts));
    cts_iv_adj  = (cts*)calloc(1, sizeof(cts));

    und        = (cts_und*)calloc(1, sizeof(cts_und));
    und_iv     = (cts_und*)calloc(1, sizeof(cts_und));
    und_iv_adj = (cts_und*)calloc(1, sizeof(cts_und));

    adi_arg = (cts_adi_arg*)calloc(1, sizeof(cts_adi_arg));

    fwd_iv_secu    = (cts_iv*)calloc(1, sizeof(cts_iv));
    fwd_iv_reserve = (cts_iv*)calloc(1, sizeof(cts_iv));

    cpn2 = (cts_cpn*)calloc(1, sizeof(cts_cpn));
    mkt  = (cts_mkt*)calloc(1, sizeof(cts_mkt));

    NumerParams = (LGMSV_NumerParams*)calloc(1, sizeof(LGMSV_NumerParams));

    if (!cts_pricing || !cts_on_iv || !cts_iv_adj || !und || !und_iv || !und_iv_adj || !adi_arg ||
        !fwd_iv_secu || !fwd_iv_reserve || !cpn2 || !mkt || !NumerParams)
    {
        err = "Memory allocation faillure in cts_caller";
        goto FREE_RETURN;
    }

    /* Added fix for data members not initialized */
    init_NULL_LGMSV_model(&(und->model));
    init_NULL_LGMSV_model(&(und_iv->model));
    init_NULL_LGMSV_model(&(und_iv_adj->model));

    /*	Numer Params initialisation */
    lgmsv_init_numer_params(
        NumerParams,
        iNbX,
        iNbSigmaXGridLeft,
        iNbSigmaXGridRight,
        dIntegParam,
        dVolLimit,
        iCalibLGM,
        iIntegMethod,
        dMinStd,
        dMaxStd);

    /*	Safety check */
    for (i = 0; i < nsmilepar; i++)
    {
        if (fabs(alphaeps[i]) < CTS_MINALPHA)
        {
            alphaeps[i] = CTS_MINALPHA;
            rhoeps[i]   = 0.0;
            ldaeps[i]   = 0.0;
        }
    }

    if (nb_factor == 2)
    {
        pde_or_mc = 1;
    }

    /*	Flags for reserve calculation */
    if (compute_reserve)
    {
        /* for one time callable reserve */
        if (!lgm_reserve)
        {
            do_one_time = 1;
        }

        if (!save_fwdiv)
        {
            save_fwdiv  = 1;
            fwd_iv_info = fwd_iv_secu;
        }

        calc_fwd_iv = 1;
    }

    /* Extra initialisation */
    und_iv->mkt     = NULL;
    und_iv_adj->mkt = NULL;
    und->mkt        = NULL;

    cts_init_cts_iv(fwd_iv_secu);

    /*	Initialise the market */
    err = cts_init_mkt(
        today,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,
        lib_freq,
        lib_basis,
        get_cash_vol,
        GetCashVolAndConvert,
        vol_type,
        cash_vol,
        num_strikes_in_vol,
        strikes_in_vol,
        mkt);

    if (err)
    {
        goto FREE_RETURN;
    }

    free_struct        = 0;
    free_struct_iv     = 0;
    free_struct_iv_adj = 0;

    call_feat        = 0;
    call_feat_iv     = 0;
    call_feat_iv_adj = 0;

    /*	If exercised */
    if (exercised)
    {
        /* Consider only the coupons which start before the ex date */
        i = 0;
        while (i < ncpn && cpn_start_date[i] < ex_date_ex)
        {
            i++;
        }
        ncpn = i;

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        /* Checks */
        if ((ncpn > 0) && (ex_date_set < cpn_pay_date[ncpn - 1]))
        {
            err =
                "Exercised settlement should be > than payment date of the last not called cf "
                "coupon ";
            goto FREE_RETURN;
        }

        if ((fund_ncpn > 0) && (ex_date_set < fund_pay[fund_ncpn - 1]))
        {
            err =
                "Exercised settlement should be > than payment date of the last not called funding "
                "coupon ";
            goto FREE_RETURN;
        }

        /* Gestion of the particular case when the Ex settlement is past *
         * or all the coupons (cf and funding) are called				 */
        if ((ex_date_set < today + eod_pay_flag) || ((ncpn == 0) && (fund_ncpn == 0)))
        {
            /* Fill the output */
            *fund_val = *cpn_val = *call_val = *mostexp_reserve = *switch_reserve = 0.0;

            /* Initial exchange of notional */
            if (start_date >= today + eod_pay_flag)
            {
                if (fund_ccy)
                {
                    *fund_val -= swp_f_df(today, start_date, fund_ccy_yc) * fund_not_ts[0] *
                                 fx_fund_dom * swp_f_df(today, fx_fund_dom_spot_date, yc) /
                                 swp_f_df(today, fx_fund_dom_spot_date, fund_ccy_yc);
                    *cpn_val -= cpn_not_ts[0] * swp_f_df(today, start_date, yc);
                }
            }

            /* Final exchange of notional */
            if (ex_date_set >= today + eod_pay_flag)
            {
                if (fund_ccy)
                {
                    /* Final exchange of notional */
                    *fund_val += swp_f_df(today, ex_date_set, fund_ccy_yc) * fund_not_ts[0] *
                                 fx_fund_dom * swp_f_df(today, fx_fund_dom_spot_date, yc) /
                                 swp_f_df(today, fx_fund_dom_spot_date, fund_ccy_yc);
                    *cpn_val += cpn_not_ts[0] * swp_f_df(today, ex_date_set, yc);
                    *call_val += -ex_fee * swp_f_df(today, ex_date_set, yc);
                }
            }

            return NULL;
        }

        /*  those particular cases should not happen */
        if ((ncpn > 0) && (fund_ncpn == 0))
        {
            err = "All the funding coupons are called but not all the cf coupons : Not allowed ";
            goto FREE_RETURN;
        }
        if ((ncpn == 0) && (fund_ncpn > 0))

        {
            err = "All the cf coupons are called but not all the funding coupons : Not allowed ";
            goto FREE_RETURN;
        }

        /*  Case where (ncpn > 0) && (fund_ncpn > 0) */
        ncall = 0;
    }

    /* save the initial fund margins */
    fund_mrg2     = (double*)calloc(fund_ncpn, sizeof(double));
    fund_spr2     = (double*)calloc(fund_ncpn, sizeof(double));
    fund_fix_cpn2 = (double*)calloc(fund_ncpn, sizeof(double));
    fund_not_ts2  = (double*)calloc(fund_ncpn, sizeof(double));

    if (!fund_mrg2 || !fund_spr2 || !fund_fix_cpn2 || !fund_not_ts2)
    {
        err = "Memory allocation error in cts_caller";
        goto FREE_RETURN;
    }

    memcpy(fund_mrg2, fund_mrg, fund_ncpn * sizeof(double));
    memcpy(fund_spr2, fund_spr, fund_ncpn * sizeof(double));
    memcpy(fund_fix_cpn2, fund_fix_cpn, fund_ncpn * sizeof(double));
    memcpy(fund_not_ts2, fund_not_ts, fund_ncpn * sizeof(double));

    if (fund_ccy == 1)
    {
        /* save the initial fund not in the third ccy
           and the third Ccy yield curve name			*/
        OtherCcyfundNot = fund_not_ts2[fund_ncpn - 1];
        strcpy(Copyfund_ccy_yc, fund_ccy_yc);

        /* Convert into domestic */

        /* check that all notionals are equals */
        if (fund_ncpn > 0)
        {
            fund_not = fund_not_ts2[0];

            for (i = 1; i < fund_ncpn; i++)
            {
                if (fabs(fund_not_ts2[i] - fund_not) > 1.0E-08)
                {
                    err = "Notional must be constant in the case of a foreign funding";
                    return err;
                }
            }
        }

        if (ncpn > 0)
        {
            cpn_not = cpn_not_ts[0];

            for (i = 1; i < ncpn; i++)
            {
                if (fabs(cpn_not_ts[i] - cpn_not) > 1.0E-08)
                {
                    err = "Notional must be constant in the case of a foreign funding";
                    return err;
                }
            }
        }

        fund_ccy = 0;
        for_fund = 1;
        err      = convert_funding_to_domestic(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            cpn_not,
            yc,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            &fund_not,
            fund_spr2,
            fund_mrg2,
            fund_fix_cpn2,
            &fund_start_date,
            &eq_final_ex,
            &eq_init_ex);
        if (err)
        {
            return err;
        }

        for (i = 0; i < fund_ncpn; i++)
        {
            fund_not_ts2[i] = fund_not;
        }
    }
    else
    {
        for_fund = 0;
    }

    /* first calculate the structure for IV */

    err = cts_fill_check_all_struct(
        mkt,
        tstar,
        use_calib,
        lambda,
        nb_factor,
        lgm_alpha,
        lgm_gamma,
        lgm_rho,
        nsmilepar,
        smilepartime,
        alphaeps,
        rhoeps,
        ldaeps,
        rho2eps,
        lgmsvund,
        fund_ncpn,
        fund_not_ts2,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr2,
        fund_mrg2,
        fix_lag_bd,
        ncpn,
        cpn_type,
        cpn_start_date,
        cpn_end_date,
        cpn_pay_date,
        cpn_coupon,
        cpn_basis,
        cpn_not_ts,
        pay_fix_date,
        pay_months,
        pay_freq,
        pay_basis,
        pay_fwd_spread,
        pay_gearing,
        nfix,
        weights,
        fix_dates,
        fix_tenor,
        fix_freq,
        fix_basis,
        fix_fwd_spreads,
        tstype,
        alpha,
        beta,
        nstr,
        str,
        nbopt,
        iv_buy_sell,
        value_zero,
        lb,
        ub,
        payoff,
        iv_call_spread,
        iv_call_spread,
        iv_trim_type,
        iv_max_fix,
        iv_min_fix_time,
        use_cmsopt,
        correl_start,
        correl_end,
        float_adjust_type,
        0,
        pay_rec,
        ex_date,
        set_date,
        fee,
        adj_fee,
        do_one_time,
        one_time_index,
        0,
        NumerParams,
        pde_or_mc,
        req_stp,
        req_stppsi,
        req_stpx,
        req_stpz,
        req_mintime,
        req_paths,
        integ_mintime,
        cal_tenor,
        cal_ref,
        cal_freq,
        cal_basis,
        force_atm,
        max_std_long,
        max_std_short,
        vol_shift_long,
        vol_type_long,
        vol_shift_type_long,
        vol_shift_short,
        vol_type_short,
        vol_shift_type_short,
        lambda_shift,
        calib_strategy,
        fix_lambda,
        short_tenor,
        short_refrate,
        short_freq,
        short_basis,
        fix_smile,
        smile_calib_months,
        lgmsv_calib_params,
        min_time,
        skip_last,
        min_fact,
        max_fact,
        use_jumps,
        numer_tstar,
        prec,
        maxiter,
        keep_first,
        long_strike_flag,
        short_strike_flag,
        NULL,
        1,
        calc_fwd_iv,
        eod_fix_flag,
        eod_ex_flag,
        cts_on_iv,
        und_iv,
        &call_feat_iv,
        adi_arg,
        save_inst_data,
        inst_data,
        save_fwdiv,
        fwd_iv_info);

    if (err)
    {
        goto FREE_RETURN;
    }

    free_struct_iv = 1;

    /* calculate the IV */
    fund_leg_pv = 0.0;

    /*	Initial exchange */
    if (cts_on_iv->fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            fund_leg_pv +=
                swp_f_df(today, cts_on_iv->fund_leg->cpn[0].start_date, yc) * eq_final_ex;
            fund_leg_pv += swp_f_df(
                               today,
                               cts_on_iv->fund_leg->cpn[cts_on_iv->fund_leg->num_cpn - 1].pay_date,
                               mkt->yc) *
                           cts_on_iv->fund_leg->cpn[cts_on_iv->fund_leg->num_cpn - 1].not ;
        }
        else
        {
            fund_leg_pv +=
                swp_f_df(today, cts_on_iv->fund_leg->cpn[0].start_date, yc) * fund_not_ts2[0];
        }
    }
    else
    {
        fin_not_date = fund_pay[fund_ncpn - 1];

        if (fin_not_date >= today + eod_pay_flag)
        {
            if (for_fund)
            {
                fund_leg_pv += swp_f_df(today, fin_not_date, yc) * eq_final_ex;
            }
        }
    }

    /* first the funding leg */

    /* initial notional */

    for (i = 0; i < cts_on_iv->fund_leg->num_cpn; i++)
    {
        fund_leg_pv += cts_on_iv->fund_leg->cpn[i].mkt_val;
    }

    /*	PV of coupons fixed in the past and not yet paid */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
    {
        if (fund_pay[i] >= today + eod_pay_flag)
        {
            err = interp_basis(fund_basis[i], &bas);
            if (err)
            {
                goto FREE_RETURN;
            }

            fund_leg_pv += (fund_fix_cpn2[i] + fund_mrg2[i]) *
                           coverage(fund_start[i], fund_end[i], bas) * fund_not_ts2[i] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	Final Notional exchange */
    if (for_fund)
    {
        if (start_date >= today + eod_pay_flag)
        {
            fund_leg_pv -= swp_f_df(today, start_date, yc) * eq_init_ex;
        }
    }

    /* Exerciced Case */
    /* Force the Notional to be paid at ex_settle_date */
    if ((exercised) && (fund_pay[fund_ncpn - 1] != ex_date_set))
    {
        /* Substract the notional received at fund_pay[fund_ncpn-1] and add notional paid at
         * ex_settle_date */
        if (for_fund)
        {
            /* Other currency case */
            fund_leg_pv += (swp_f_df(today, ex_date_set, Copyfund_ccy_yc) -
                            swp_f_df(today, fund_pay[fund_ncpn - 1], Copyfund_ccy_yc)) *
                           OtherCcyfundNot * fx_fund_dom *
                           swp_f_df(today, fx_fund_dom_spot_date, yc) /
                           swp_f_df(today, fx_fund_dom_spot_date, Copyfund_ccy_yc);
        }
    }

    exo_leg_pv = 0.0;

    if (!compatibility_flag)
    {
        /* then the exotic leg */
        for (i = 0; i < cts_on_iv->exo_leg->num_cpn; i++)
        {
            exo_leg_pv += cts_on_iv->exo_leg->cpn[i].mkt_val;
        }

        i = 0;
        while (i < ncpn && fix_dates[i][0] < today + eod_fix_flag)
        {
            if (cpn_pay_date[i] >= today + eod_pay_flag)
            {
                /* we are in the middle of a coupon */
                j0 = 0;

                while (j0 < nfix[i] && fix_dates[i][j0] < today + eod_fix_flag)
                {
                    j0++;
                }

                err = cts_calc_fixed_cpn_mkt_value(
                    mkt,
                    cpn_start_date[i],
                    cpn_end_date[i],
                    cpn_pay_date[i],
                    cpn_coupon[i],
                    cpn_basis[i],
                    cpn_not_ts[i],
                    fix_lag_bd,
                    cpn_type[i],
                    pay_fix_date[i],
                    pay_months,
                    pay_freq,
                    pay_basis,
                    pay_fwd_spread[i],
                    pay_gearing[i],
                    fixed_pay[i],
                    j0,
                    weights[i],
                    fixed_ref[i],
                    tstype,
                    alpha[i],
                    beta[i],
                    nstr[i],
                    str[i],
                    nbopt[i],
                    lb[i],
                    ub[i],
                    payoff[i],
                    accrue_on_barrier,
                    eod_fix_flag,
                    eod_pay_flag,
                    eod_ex_flag,
                    &cpn_pv);

                exo_leg_pv += cpn_pv;

                if (j0 < nfix[i])
                {
                    /* second part must be evaluated */
                    if (pay_fix_date[i] < today + eod_fix_flag &&
                        (cpn_type[i] == 1 || cpn_type[i] == 2))
                    {
                        cpn_coupon_temp = cpn_coupon[i] + pay_gearing[i] * fixed_pay[i];
                        cpntype_temp    = 0;
                    }
                    else
                    {
                        cpn_coupon_temp = cpn_coupon[i];
                        cpntype_temp    = cpn_type[i];
                    }

                    err = cts_init_cpn(
                        mkt,
                        cpn_start_date[i],
                        cpn_end_date[i],
                        cpn_pay_date[i],
                        cpn_coupon_temp,
                        cpn_basis[i],
                        cpn_not_ts[i],
                        fix_lag_bd,
                        cpntype_temp,
                        pay_fix_date[i],
                        pay_months,
                        pay_freq,
                        pay_basis,
                        pay_fwd_spread[i],
                        pay_gearing[i],
                        nfix[i] - j0,
                        &(weights[i][j0]),
                        &(fix_dates[i][j0]),
                        fix_tenor[i],
                        fix_freq,
                        fix_basis,
                        &(fix_fwd_spreads[i][j0]),
                        tstype,
                        alpha[i],
                        beta[i],
                        nstr[i],
                        str[i],
                        nbopt[i],
                        iv_buy_sell,
                        value_zero,
                        lb[i],
                        ub[i],
                        payoff[i],
                        iv_call_spread,
                        iv_call_spread,
                        iv_trim_type,
                        iv_max_fix,
                        iv_min_fix_time,
                        1,
                        use_cmsopt,
                        correl_start,
                        correl_end,
                        float_adjust_type,
                        cpn2);

                    if (err)
                    {
                        goto FREE_RETURN;
                    }

                    exo_leg_pv += cpn2->mkt_val;

                    if (cpn2)
                    {
                        cts_free_coupon(cpn2);
                    }
                }
            }

            i++;
        }

        /*	Initial and final exchange */
        if (for_fund)
        {
            /*	Final */
            if (cts_on_iv->exo_leg->num_cpn > 0)
            {
                exo_leg_pv +=
                    swp_f_df(
                        today,
                        cts_on_iv->exo_leg->cpn[cts_on_iv->exo_leg->num_cpn - 1].cpn_pay_date,
                        yc) *
                    cts_on_iv->exo_leg->cpn[cts_on_iv->exo_leg->num_cpn - 1].cpn_not;
            }
            else
            {
                fin_not_date = fund_pay[fund_ncpn - 1];

                if (fin_not_date >= today + eod_pay_flag)
                {
                    exo_leg_pv += swp_f_df(today, fin_not_date, yc) * cpn_not_ts[0];
                }
            }

            /*	Initial */
            if (start_date >= today + eod_pay_flag)
            {
                exo_leg_pv -= swp_f_df(today, start_date, yc) * cts_on_iv->exo_leg->cpn[0].cpn_not;
            }
        }
    }
    else
    {
        /* evaluate only the called coupon */

        i = 0;
        while (i < ncall && ex_date[i] < today + eod_ex_flag)
        {
            i++;
        }

        if (i < ncall)
        {
            j = 0;
            while (j < cts_on_iv->exo_leg->num_cpn &&
                   cts_on_iv->exo_leg->cpn[j].cpn_start_date < ex_date[i])
            {
                j++;
            }

            for (j = j; j < cts_on_iv->exo_leg->num_cpn; j++)
            {
                exo_leg_pv += cts_on_iv->exo_leg->cpn[j].mkt_val;
            }
        }
    }

    /* eventuall evaluate the call */
    if (ncall > 0)
    {
        if (iv_buy_sell != call_buy_sell && calc_fwd_iv)
        {
            err = cts_fill_check_all_struct(
                mkt,
                tstar,
                use_calib,
                lambda,
                nb_factor,
                lgm_alpha,
                lgm_gamma,
                lgm_rho,
                nsmilepar,
                smilepartime,
                alphaeps,
                rhoeps,
                ldaeps,
                rho2eps,
                lgmsvund,
                fund_ncpn,
                fund_not_ts2,
                fund_fix,
                fund_start,
                fund_end,
                fund_pay,
                fund_basis,
                fund_spr2,
                fund_mrg2,
                fix_lag_bd,
                ncpn,
                cpn_type,
                cpn_start_date,
                cpn_end_date,
                cpn_pay_date,
                cpn_coupon,
                cpn_basis,
                cpn_not_ts,
                pay_fix_date,
                pay_months,
                pay_freq,
                pay_basis,
                pay_fwd_spread,
                pay_gearing,
                nfix,
                weights,
                fix_dates,
                fix_tenor,
                fix_freq,
                fix_basis,
                fix_fwd_spreads,
                tstype,
                alpha,
                beta,
                nstr,
                str,
                nbopt,
                call_buy_sell,
                value_zero,
                lb,
                ub,
                payoff,
                iv_call_spread,
                iv_call_spread,
                iv_trim_type,
                iv_max_fix,
                iv_min_fix_time,
                use_cmsopt,
                correl_start,
                correl_end,
                float_adjust_type,
                0,
                pay_rec,
                ex_date,
                set_date,
                fee,
                adj_fee,
                do_one_time,
                one_time_index,
                0,
                NumerParams,
                pde_or_mc,
                req_stp,
                req_stppsi,
                req_stpx,
                req_stpz,
                req_mintime,
                req_paths,
                integ_mintime,
                cal_tenor,
                cal_ref,
                cal_freq,
                cal_basis,
                force_atm,
                max_std_long,
                max_std_short,
                vol_shift_long,
                vol_type_long,
                vol_shift_type_long,
                vol_shift_short,
                vol_type_short,
                vol_shift_type_short,
                lambda_shift,
                calib_strategy,
                fix_lambda,
                short_tenor,
                short_refrate,
                short_freq,
                short_basis,
                fix_smile,
                smile_calib_months,
                lgmsv_calib_params,
                min_time,
                skip_last,
                min_fact,
                max_fact,
                use_jumps,
                numer_tstar,
                prec,
                maxiter,
                keep_first,
                long_strike_flag,
                short_strike_flag,
                NULL,
                1,
                calc_fwd_iv,
                eod_fix_flag,
                eod_ex_flag,
                cts_iv_adj,
                und_iv_adj,
                &call_feat_iv_adj,
                adi_arg,
                save_inst_data,
                inst_data,
                save_fwdiv,
                fwd_iv_info);

            if (err)
            {
                goto FREE_RETURN;
            }

            free_struct_iv_adj = 1;

            cts_for_iv = cts_iv_adj;
        }
        else
        {
            cts_for_iv = cts_on_iv;
        }

        err = cts_fill_check_all_struct(
            mkt,
            tstar,
            use_calib,
            lambda,
            nb_factor,
            lgm_alpha,
            lgm_gamma,
            lgm_rho,
            nsmilepar,
            smilepartime,
            alphaeps,
            rhoeps,
            ldaeps,
            rho2eps,
            lgmsvund,
            fund_ncpn,
            fund_not_ts2,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr2,
            fund_mrg2,
            fix_lag_bd,
            ncpn,
            cpn_type,
            cpn_start_date,
            cpn_end_date,
            cpn_pay_date,
            cpn_coupon,
            cpn_basis,
            cpn_not_ts,
            pay_fix_date,
            pay_months,
            pay_freq,
            pay_basis,
            pay_fwd_spread,
            pay_gearing,
            nfix,
            weights,
            fix_dates,
            fix_tenor,
            fix_freq,
            fix_basis,
            fix_fwd_spreads,
            tstype,
            alpha,
            beta,
            nstr,
            str,
            nbopt,
            call_buy_sell,
            value_zero,
            lb,
            ub,
            payoff,
            call_call_spread,
            numer_call_spread,
            call_trim_type,
            call_max_fix,
            call_min_fix_time,
            use_cmsopt,
            correl_start,
            correl_end,
            float_adjust_type,
            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,
            adj_fee,
            do_one_time,
            one_time_index,
            0,
            NumerParams,
            pde_or_mc,
            req_stp,
            req_stppsi,
            req_stpx,
            req_stpz,
            req_mintime,
            req_paths,
            integ_mintime,
            cal_tenor,
            cal_ref,
            cal_freq,
            cal_basis,
            force_atm,
            max_std_long,
            max_std_short,
            vol_shift_long,
            vol_type_long,
            vol_shift_type_long,
            vol_shift_short,
            vol_type_short,
            vol_shift_type_short,
            lambda_shift,
            calib_strategy,
            fix_lambda,
            short_tenor,
            short_refrate,
            short_freq,
            short_basis,
            fix_smile,
            smile_calib_months,
            lgmsv_calib_params,
            min_time,
            skip_last,
            min_fact,
            max_fact,
            use_jumps,
            numer_tstar,
            prec,
            maxiter,
            keep_first,
            long_strike_flag,
            short_strike_flag,
            cts_for_iv,
            0,
            calc_fwd_iv,
            eod_fix_flag,
            eod_ex_flag,
            cts_pricing,
            und,
            &call_feat,
            adi_arg,
            save_inst_data,
            inst_data,
            save_fwdiv,
            fwd_iv_info);

        if (err)
        {
            goto FREE_RETURN;
        }

        free_struct = 1;

        /*	3)	Adjust the fee is needed */
        err = cts_adjust_model_fwd_iv(
            und,
            cts_for_iv,
            cts_pricing,
            for_fund,
            calc_fwd_iv,
            adj_fee,
            pde_or_mc,
            save_fwdiv,
            fwd_iv_info);

        if (err)
        {
            goto FREE_RETURN;
        }
    }

    /*	4)	If there is at least one call after today, value call feature */

    if (call_feat == 1)
    {
        err = cts_launch_algo(cts_pricing, und, adi_arg, &call, &(onetime[0]), &fwd_iv);

        if (err)
            goto FREE_RETURN;

        if (save_extra_infos && extra_infos)
        {
            extra_infos->algo_fwd_iv1 = fwd_iv;
        }
    }
    else
    {
        call = 0.0;
    }

    /* Exerciced Case */
    /* Force the Notional to be paid at ex_settle_date */
    if ((exercised) && (fund_pay[fund_ncpn - 1] != ex_date_set))
    {
        /* Substract the notional received at fund_pay[fund_ncpn-1] and add notional paid at
         * ex_settle_date */
        if (for_fund)
        {
            exo_leg_pv +=
                (swp_f_df(today, ex_date_set, yc) - swp_f_df(today, fund_pay[fund_ncpn - 1], yc)) *
                cpn_not_ts[ncpn - 1];
        }

        /* Add fee */
        if (ex_date_set >= today + eod_pay_flag)
        {
            call += -ex_fee * swp_f_df(today, ex_date_set, yc);
        }
    }

    *fund_val = fund_leg_pv;
    *cpn_val  = exo_leg_pv;
    *call_val = call;

    if (do_one_time)
    {
        if (one_time_index > 0)
        {
            onetime_val[0] = onetime[0];

            for (i = 1; i < cts_pricing->num_calls; i++)
            {
                onetime_val[i] = 0.0;
            }
        }
        else
        {
            for (i = 0; i < cts_pricing->num_calls; i++)
            {
                onetime_val[i] = onetime[i];
            }
        }
    }

    if (export_und && call_feat == 1)
    {
        cts_copy_und(und, und_exp);

        Convert_Tstar_model(&(und_exp->model), und_exp->model.dInitTStar);

        ConvertTS_LGMSV_to_LGM(
            und_exp->model.iNbPWTime,
            und_exp->model.dPWTime,
            und_exp->model.dSigma,
            und_exp->model.dLambdaX,
            und_exp->model.dTStar);

        /* divide the smile parameters by 2 */
        for (i = 0; i < und_exp->model.iNbPWTime; i++)
        {
            und_exp->model.dAlpha[i] /= 2.0;
            und_exp->model.dLambdaEps[i] /= 2.0;
        }
    }

    /* All about reserve calculation now */
    if (compute_reserve && call_feat == 1)
    {
        err = cts_calc_reserve(
            cts_for_iv,
            cts_pricing,
            und,
            adi_arg,
            for_fund,
            call,
            &(onetime[0]),
            nsmilepar,
            smilepartime,
            alphaeps,
            rhoeps,
            ldaeps,
            rho2eps,
            tstar,
            NumerParams,
            numer_tstar,
            cal_tenor,
            cal_ref,
            cal_freq,
            cal_basis,
            force_atm,
            max_std_long,
            max_std_short,
            vol_shift_long,
            vol_type_long,
            vol_shift_type_long,
            vol_shift_short,
            vol_type_short,
            vol_shift_type_short,
            lambda_shift,
            calib_strategy,
            fix_lambda,
            short_tenor,
            short_refrate,
            short_freq,
            short_basis,
            fix_smile,
            smile_calib_months,
            lgmsv_calib_params,
            min_time,
            skip_last,
            min_fact,
            max_fact,
            use_jumps,
            prec,
            maxiter,
            keep_first,
            long_strike_flag,
            short_strike_flag,
            calc_fwd_iv,
            adj_fee,
            lambda_reserve,
            lgm_alpha,
            lgm_gamma,
            lgm_rho,
            one_time_vega,
            vol_shift_type_long,
            reserve_method,
            lgm_reserve,
            midat_reserve,
            recalib_european,
            euro_nb_iter,
            one_time_index,
            recalc_one_factor,
            lgm_nstpt,
            lgm_nstpx,
            req_stppsi,
            req_stpz,
            integ_mintime,
            switch_reserve,
            mostexp_reserve,
            fwd_iv_info,
            save_extra_infos,
            extra_infos);

        if (err)
            goto FREE_RETURN;
    }
    else
    {
        *mostexp_reserve = 0.0;
        *switch_reserve  = 0.0;
    }

FREE_RETURN:

    if (free_struct_iv)
    {
        cts_free_all_struct(cts_on_iv, und_iv, call_feat_iv, NULL);
    }

    if (free_struct_iv_adj)
    {
        cts_free_all_struct(cts_iv_adj, und_iv_adj, call_feat_iv_adj, NULL);
    }

    if (free_struct)
    {
        cts_free_all_struct(cts_pricing, und, call_feat, adi_arg);
    }

    cts_free_cts_iv(fwd_iv_secu);

    if (cts_pricing)
        free(cts_pricing);
    if (cts_on_iv)
        free(cts_on_iv);
    if (cts_iv_adj)
        free(cts_iv_adj);

    if (und)
        free(und);
    if (und_iv)
        free(und_iv);
    if (und_iv_adj)
        free(und_iv_adj);

    if (adi_arg)
        free(adi_arg);

    if (fwd_iv_secu)
        free(fwd_iv_secu);
    if (fwd_iv_reserve)
        free(fwd_iv_reserve);

    if (cpn2)
        free(cpn2);
    if (mkt)
        free(mkt);

    if (fund_mrg2)
        free(fund_mrg2);
    if (fund_spr2)
        free(fund_spr2);
    if (fund_fix_cpn2)
        free(fund_fix_cpn2);
    if (fund_not_ts2)
        free(fund_not_ts2);

    if (NumerParams)
        free(NumerParams);

    return err;
}
