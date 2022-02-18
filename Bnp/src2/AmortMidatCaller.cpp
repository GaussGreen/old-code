
#include "amortmidatcaller.h"

#include "amortmidatprodstruct.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define NUM_HERMITE 10

static double dens(double x)
{
    return INV_SQRT_TWO_PI * exp(-x * x / 2.0);
}

/*
Convert the funding leg into a domestic one
                -	spreads and margins
                -	past fixings
                -	final notional exchange
                -	initial notional exchange
*/

Err convert_funding_to_domestic_amort(
    //	Inputs
    long today,        //	Today
    long not_ex_date,  //	Date at which the
                       //		initial notional
                       //		exchange takes place
                       //		(or has taken place)
    int eod_fix_flag,  //	0: I, 1: E
    int eod_pay_flag,  //	0: I, 1: E

    double fx_fund_dom,            //	Fx fund/dom, 2bd fwd
    long   fx_fund_dom_spot_date,  //	Spot date for Fx

    char*   dom_yc,      //	Domestic discount curve
    int     nb_dom_not,  //	Number of Domestic notional
    double* dom_not,     //	Domestic notional

    int    fund_ncpn,   //	Number of coupon
    long*  fund_fix,    //	Fixing dates
    long*  fund_start,  //	Start dates
    long*  fund_pay,    //	Pay dates
    char** fund_basis,  //	Basis
    //	The following are modified
    char* fund_yc,         //	Funding discount curve,
                           //	changed to domestic
    double* fund_not,      //	Funding notional
                           //		in funding ccy,
                           //	converted to domestic
    double* fund_spr,      //	Spread in funding ccy,
                           //	put to 0
    double* fund_mrg,      //	Margin in funding ccy,
                           //	converted to margin over
                           //	cash libor in domestic currency
    double* fund_fix_cpn,  //	Fixing: contains spread
                           //		but not margin,
                           //	converted to equivalent
                           //	domestic cash-flows
    //	The following are returned
    long*   fund_start_date,  //	Start date of the funding
    double* eq_final_ex,      //	Domestic cash-flow equivalent
                              //		to final exchange
                              //	(to be delivered
                              //	at funding start date)
    double* eq_init_ex)       //	Domestic cash-flow equivalent
                              //		to initial exchange
                              //	(to be delivered
                              //	at initial exchange date)
{
    int    i, i0;
    double ffx;
    double not_ratio;

    //	Find the start date
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
    {
        i++;
    }
    i0 = i;

    if (i0 == fund_ncpn)
    {
        *fund_start_date = fund_pay[i0 - 1];
    }
    else
    {
        *fund_start_date = fund_start[i0];
    }

    not_ratio = dom_not[nb_dom_not - 1] / fund_not[fund_ncpn - 1];

    //	Adjust fixings
    for (i = 0; i < i0; i++)
    {
        if (fund_pay[i] >= today + eod_pay_flag)
        {
            ffx = fx_fund_dom * swp_f_df(fx_fund_dom_spot_date, fund_pay[i], fund_yc) /
                  swp_f_df(fx_fund_dom_spot_date, fund_pay[i], dom_yc);
            fund_fix_cpn[i] = not_ratio * ffx * (fund_fix_cpn[i] + fund_mrg[i]);
        }
        else
        {
            fund_fix_cpn[i] = 0.0;
        }
        fund_mrg[i] = 0.0;
        fund_not[i] = fund_not[i] * not_ratio;
    }

    //	Adjust spreads and margins
    for (i = i0; i < fund_ncpn; i++)
    {
        ffx = fx_fund_dom * swp_f_df(fx_fund_dom_spot_date, fund_pay[i], fund_yc) /
              swp_f_df(fx_fund_dom_spot_date, fund_pay[i], dom_yc);

        fund_mrg[i] = not_ratio * ffx * (fund_spr[i] + fund_mrg[i]);
        fund_spr[i] = 0.0;

        fund_not[i] = fund_not[i] * not_ratio;
    }

    return NULL;
}

Err am_compute_american_fee(
    long today,
    /*	yc */
    char* yc,

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,

    /*		calls */
    int     ncall,
    int     pay_rec,
    long*   ex_date,
    long*   set_date,
    double* fee)
{
    Err          err = NULL;
    int          i;
    SrtBasisCode bas;
    double       x;

    /*	float coupon fixed in the past and not yet paid */
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

            x = (fund_fix_cpn[i] + fund_mrg[i]) * coverage(fund_start[i], fund_pay[i], bas) *
                fund_not[0] * swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	fix coupon of past period and not yet paid */
    i = 0;
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            err = interp_basis(fix_basis[i], &bas);
            if (err)
            {
                goto FREE_RETURN;
            }

            x = fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

FREE_RETURN:

    return err;
}

Err am_caller(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */

    /*		if calib */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char*   fund_ref,
    int     fund_ccy,    /*	0: domestic, 1: other */
    double* fund_not_ts, /*	If different from domestic or foreign (fund_ccy = 1) */
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
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee, /*	Exercise Fee */

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    double mintime,
    double mininterval,
    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    /*	Results */
    double* fund_val,  /*	Value of the funding leg */
    double* fix_val,   /*	Value of the Power Dual leg */
    double* call_val,  /*	Value of the callable feature */
    int     export_ts, /*	1: Export TS, 0: don't */
    AM_UND  und_exp)
{
    am_str*     am      = NULL;
    am_und*     und     = NULL;
    am_adi_arg* adi_arg = NULL;

    AM_FUND_LEG fund_leg;
    AM_FUND_CPN fund_cpn;
    AM_FIX_LEG  fix_leg;
    AM_FIX_CPN  fix_cpn;

    int          call_feat;
    double       fund_leg_pv, fix_leg_pv;
    SrtBasisCode bas;

    double df, coupon;

    double call;
    int    i, j;
    int    free_struct = 0;

    double temp;

    int    for_fund;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;

    double *fund_mrg2 = NULL, *fund_spr2 = NULL, *fund_not_ts2 = NULL, *fund_fix_cpn2 = NULL;

    Err err = NULL;

    /*	If exercised */
    if (exercised)
    {
        i = 0;
        while (i < fix_ncpn && fix_start[i] < ex_date_ex)
        {
            i++;
        }
        fix_ncpn = i;

        /*	Structure is called before start: return 0 */
        if (fix_ncpn == 0)
        {
            *fund_val = *fix_val = *call_val = 0.0;
            return NULL;
        }

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        ncall     = 0;
        exercised = 0;

        err = am_caller(
            today,
            use_calib,
            yc,
            vc,
            ref,
            swap_freq,
            swap_basis,
            lambda,
            alpha,
            gamma,
            rho,
            get_cash_vol,
            lgm2dund,
            start_date,
            theoEndDate,
            fund_ref,
            fund_ccy,
            fund_not_ts,
            fund_ccy_yc,
            fx_fund_dom,
            fx_fund_dom_spot_date,

            fund_ncpn,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,

            fix_not,
            fix_ncpn,
            fix_start,
            fix_end,
            fix_pay,
            fix_basis,
            fix_rate,
            fix_fee,

            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,

            req_stp,
            req_stpx,

            mintime,
            mininterval,
            notperiod,
            one2F,
            use_jump,
            max_var_jump,
            strike_type,
            european_model,

            get_correl,
            CorrelName,

            max_std_short,
            fix_lambda,
            one_f_equi,
            skip_last,

            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,

            exercised,
            ex_date_ex,
            ex_date_set,
            ex_fee,
            fund_val,
            fix_val,
            call_val,
            export_ts,
            und_exp);

        if (err)
        {
            goto FREE_RETURN;
        }

        if (ex_date_set >= today + eod_pay_flag)
        {
            *call_val = -ex_fee * swp_f_df(today, ex_date_set, yc);
        }

        goto FREE_RETURN;
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
        fund_ccy = 0;
        for_fund = 1;
        err      = convert_funding_to_domestic_amort(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            yc,
            fix_ncpn,
            fix_not,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            fund_not_ts2,
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
    }
    else
    {
        for_fund = 0;
    }

    am      = (am_str*)calloc(1, sizeof(am_str));
    und     = (am_und*)calloc(1, sizeof(am_und));
    adi_arg = (am_adi_arg*)calloc(1, sizeof(am_adi_arg));

    if (!am || !und || !adi_arg)
    {
        err = "memory allocation failure in am_caller";
        goto FREE_RETURN;
    }

    /*	Initialise structures */
    free_struct = 0;
    err         = am_fill_check_all_struct(
        today,
        theoEndDate,

        use_calib,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,
        lambda,
        alpha,
        gamma,
        rho,

        lgm2dund,

        fund_ref,
        fund_not_ts2,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr2,
        fund_mrg2,

        fix_not,
        fix_ncpn,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,

        ncall,
        pay_rec,
        ex_date,
        set_date,
        fee,

        req_stp,
        req_stpx,

        get_cash_vol,

        mintime,
        mininterval,
        notperiod,
        one2F,
        use_jump,
        max_var_jump,
        strike_type,
        european_model,

        get_correl,
        CorrelName,
        /*
                        GetVolForBadr,
                        cVolType,
        */
        max_std_short,
        fix_lambda,
        one_f_equi,
        skip_last,

        eod_fix_flag,
        eod_ex_flag,

        am,
        und,

        &call_feat,
        adi_arg);

    if (err)
    {
        goto FREE_RETURN;
    }
    free_struct = 1;

    /*	1)	Value funding leg */

    fund_leg    = am->fund_leg;
    fund_leg_pv = 0.0;

    /*	Cash libor */
    if (fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            for (j = 0; j < fund_leg->num_cpn - 1; ++j)
            {
                if (fund_leg->cpn[j].pay_date >= today + eod_pay_flag)
                {
                    temp = swp_f_df(today, fund_leg->cpn[j].pay_date, yc) *
                           (fund_not_ts2[j] - fund_not_ts2[j + 1]);
                    fund_leg_pv += temp;
                }
            }

            if (fund_leg->cpn[fund_leg->num_cpn - 1].pay_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date, yc) *
                       fund_not_ts2[fund_leg->num_cpn - 1];
                fund_leg_pv += temp;
            }
        }
        else
        {
            fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, yc) * fund_leg->notional[0];
        }
    }

    /*	Coupons: spread + margin */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        fund_cpn = fund_leg->cpn + i;
        temp     = swp_f_df(today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        fund_leg_pv += temp;
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
                           coverage(fund_start[i], fund_pay[i], bas) * fund_not_ts2[0] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	2)	Value fix leg */

    fix_leg    = am->fix_leg;
    fix_leg_pv = 0.0;

    /*	Coupons */
    for (i = 0; i < fix_leg->num_cpn; i++)
    {
        fix_cpn = fix_leg->cpn + i;

        /*	Discount */
        df = swp_f_df(today, fix_cpn->pay_date, yc);

        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        temp = df * coupon;
        fix_leg_pv += temp;
    }

    /*	PV of coupons of past periods and not yet paid */
    i = 0;
    //	while (i < fix_ncpn && fix_fix[i] < today + eod_fix_flag)
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            err = interp_basis(fix_basis[i], &bas);
            if (err)
            {
                goto FREE_RETURN;
            }

            //			fix_leg_pv += fix_fix_cpn[i]
            fix_leg_pv += fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                          swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

    /*	Initial and final exchange */
    if (for_fund)
    {
        /*	Final */
        if (fix_leg->num_cpn > 0)
        {
            for (j = i; j < fix_leg->num_cpn - 1; ++j)
            {
                temp =
                    swp_f_df(today, fix_leg->cpn[j].pay_date, yc) * (fix_not[j] - fix_not[j + 1]);
                fix_leg_pv += temp;
            }

            temp = swp_f_df(today, fix_leg->cpn[fix_leg->num_cpn - 1].pay_date, yc) *
                   fix_not[fix_leg->num_cpn - 1];
            fix_leg_pv += temp;
        }
        else
        {
            fin_not_date = fund_pay[fund_ncpn - 1];
            if (fin_not_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fin_not_date, yc) * fix_not[fund_leg->num_cpn - 1];
                fix_leg_pv += temp;
            }
        }

        /*	Initial */
        if (start_date >= today + eod_pay_flag)
        {
            fix_leg_pv -= swp_f_df(today, start_date, yc) * fix_not[0];
        }
    }

    /*	4)	If there is at least one call after today, value call feature */

    if (call_feat == 1)
    {
        smessage("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);

        err = am_launch_adi(am, und, adi_arg, &call);

        if (err)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        call = 0.0;
    }

    *fund_val = fund_leg_pv;
    *fix_val  = fix_leg_pv;
    *call_val = call;

    if (export_ts)
    {
        am_copy_und(und, und_exp);
    }

FREE_RETURN:

    if (free_struct)
    {
        am_free_all_struct(am, und, call_feat, adi_arg);
    }
    if (am)
        free(am);
    if (und)
        free(und);
    if (adi_arg)
        free(adi_arg);

    if (fund_mrg2)
        free(fund_mrg2);
    if (fund_spr2)
        free(fund_spr2);
    if (fund_fix_cpn2)
        free(fund_fix_cpn2);
    if (fund_not_ts2)
        free(fund_not_ts2);

    return err;
}

Err am_caller_ts(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */

    /*		if calib */
    char*   yc,            /*	yc */
    char*   vc,            /*	vc */
    char*   ref,           /*	ref rate (only if calib) */
    char*   swap_freq,     /*	swap freq (only if calib) */
    char*   swap_basis,    /*	swap basis (only if calib) */
    double* pdLambdaValue, /*	lambda if unique */
    double* pdLambdaTime,
    int     nLambdaSize,
    double  alpha, /*	alpha */
    double  gamma, /*	gamma */
    double  rho,   /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char*   fund_ref,
    int     fund_ccy,    /*	0: domestic, 1: other */
    double* fund_not_ts, /*	If different from domestic or foreign (fund_ccy = 1) */
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
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    //			long		*fix_fix,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee, /*	Exercise Fee */
                     //			double		*fix_fix_cpn,			/*	Past coupon fixing if
                     //relevant
                     //*/

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    double mintime,
    double mininterval,
    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    /*
                            Err (*GetVolForBadr)( Date, Date, double, SRT_Boolean, double *),
                            char *cVolType,
    */
    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    /*	Results */
    double* fund_val,  /*	Value of the funding leg */
    double* fix_val,   /*	Value of the Power Dual leg */
    double* call_val,  /*	Value of the callable feature */
    int     export_ts, /*	1: Export TS, 0: don't */
    AM_UND  und_exp)
{
    am_str*     am      = NULL;
    am_und*     und     = NULL;
    am_adi_arg* adi_arg = NULL;

    AM_FUND_LEG fund_leg;
    AM_FUND_CPN fund_cpn;
    AM_FIX_LEG  fix_leg;
    AM_FIX_CPN  fix_cpn;

    int          call_feat;
    double       fund_leg_pv, fix_leg_pv;
    SrtBasisCode bas;

    double df, coupon;

    double call;
    int    i, j;
    int    free_struct = 0;

    double temp;

    int    for_fund;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;

    double *fund_mrg2 = NULL, *fund_spr2 = NULL, *fund_not_ts2 = NULL, *fund_fix_cpn2 = NULL;

    Err err = NULL;

    /*	If exercised */
    if (exercised)
    {
        i = 0;
        while (i < fix_ncpn && fix_start[i] < ex_date_ex)
        {
            i++;
        }
        fix_ncpn = i;

        /*	Structure is called before start: return 0 */
        if (fix_ncpn == 0)
        {
            *fund_val = *fix_val = *call_val = 0.0;
            return NULL;
        }

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        ncall     = 0;
        exercised = 0;

        err = am_caller_ts(
            today,
            use_calib,
            yc,
            vc,
            ref,
            swap_freq,
            swap_basis,
            pdLambdaValue,
            pdLambdaTime,
            nLambdaSize,

            alpha,
            gamma,
            rho,
            get_cash_vol,
            lgm2dund,
            start_date,
            theoEndDate,
            fund_ref,
            fund_ccy,
            fund_not_ts,
            fund_ccy_yc,
            fx_fund_dom,
            fx_fund_dom_spot_date,

            fund_ncpn,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,

            fix_not,
            fix_ncpn,
            fix_start,
            fix_end,
            fix_pay,
            fix_basis,
            fix_rate,
            fix_fee,

            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,

            req_stp,
            req_stpx,

            mintime,
            mininterval,
            notperiod,
            one2F,
            use_jump,
            max_var_jump,
            strike_type,
            european_model,

            get_correl,
            CorrelName,

            max_std_short,
            fix_lambda,
            one_f_equi,
            skip_last,

            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,

            exercised,
            ex_date_ex,
            ex_date_set,
            ex_fee,
            fund_val,
            fix_val,
            call_val,
            export_ts,
            und_exp);

        if (err)
        {
            goto FREE_RETURN;
        }

        if (ex_date_set >= today + eod_pay_flag)
        {
            *call_val = -ex_fee * swp_f_df(today, ex_date_set, yc);
        }

        goto FREE_RETURN;
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
        fund_ccy = 0;
        for_fund = 1;
        err      = convert_funding_to_domestic_amort(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            yc,
            fix_ncpn,
            fix_not,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            fund_not_ts2,
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
    }
    else
    {
        for_fund = 0;
    }

    am      = (am_str*)calloc(1, sizeof(am_str));
    und     = (am_und*)calloc(1, sizeof(am_und));
    adi_arg = (am_adi_arg*)calloc(1, sizeof(am_adi_arg));

    if (!am || !und || !adi_arg)
    {
        err = "memory allocation failure in am_caller";
        goto FREE_RETURN;
    }

    /*	Initialise structures */
    free_struct = 0;
    err         = am_fill_check_all_struct_ts(
        today,
        theoEndDate,

        use_calib,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,

        nLambdaSize,
        pdLambdaTime,
        pdLambdaValue,

        alpha,
        gamma,
        rho,

        lgm2dund,

        fund_ref,
        fund_not_ts2,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr2,
        fund_mrg2,

        fix_not,
        fix_ncpn,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,

        ncall,
        pay_rec,
        ex_date,
        set_date,
        fee,

        req_stp,
        req_stpx,

        get_cash_vol,

        mintime,
        mininterval,
        notperiod,
        one2F,
        use_jump,
        max_var_jump,
        strike_type,
        european_model,

        get_correl,
        CorrelName,
        /*
                        GetVolForBadr,
                        cVolType,
        */
        max_std_short,
        fix_lambda,
        one_f_equi,
        skip_last,

        eod_fix_flag,
        eod_ex_flag,

        am,
        und,

        &call_feat,
        adi_arg);

    if (err)
    {
        goto FREE_RETURN;
    }
    free_struct = 1;

    /*	1)	Value funding leg */

    fund_leg    = am->fund_leg;
    fund_leg_pv = 0.0;

    /*	Cash libor */
    if (fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            for (j = 0; j < fund_leg->num_cpn - 1; ++j)
            {
                if (fund_leg->cpn[j].pay_date >= today + eod_pay_flag)
                {
                    temp = swp_f_df(today, fund_leg->cpn[j].pay_date, yc) *
                           (fund_not_ts2[j] - fund_not_ts2[j + 1]);
                    fund_leg_pv += temp;
                }
            }

            if (fund_leg->cpn[fund_leg->num_cpn - 1].pay_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date, yc) *
                       fund_not_ts2[fund_leg->num_cpn - 1];
                fund_leg_pv += temp;
            }
        }
        else
        {
            fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, yc) * fund_leg->notional[0];
        }
    }

    /*	Coupons: spread + margin */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        fund_cpn = fund_leg->cpn + i;
        temp     = swp_f_df(today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        fund_leg_pv += temp;
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
                           coverage(fund_start[i], fund_pay[i], bas) * fund_not_ts2[0] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	2)	Value fix leg */

    fix_leg    = am->fix_leg;
    fix_leg_pv = 0.0;

    /*	Coupons */
    for (i = 0; i < fix_leg->num_cpn; i++)
    {
        fix_cpn = fix_leg->cpn + i;

        /*	Discount */
        df = swp_f_df(today, fix_cpn->pay_date, yc);

        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        temp = df * coupon;
        fix_leg_pv += temp;
    }

    /*	PV of coupons of past periods and not yet paid */
    i = 0;
    //	while (i < fix_ncpn && fix_fix[i] < today + eod_fix_flag)
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            err = interp_basis(fix_basis[i], &bas);
            if (err)
            {
                goto FREE_RETURN;
            }

            //			fix_leg_pv += fix_fix_cpn[i]
            fix_leg_pv += fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                          swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

    /*	Initial and final exchange */
    if (for_fund)
    {
        /*	Final */
        if (fix_leg->num_cpn > 0)
        {
            for (j = i; j < fix_leg->num_cpn - 1; ++j)
            {
                temp =
                    swp_f_df(today, fix_leg->cpn[j].pay_date, yc) * (fix_not[j] - fix_not[j + 1]);
                fix_leg_pv += temp;
            }

            temp = swp_f_df(today, fix_leg->cpn[fix_leg->num_cpn - 1].pay_date, yc) *
                   fix_not[fix_leg->num_cpn - 1];
            fix_leg_pv += temp;
        }
        else
        {
            fin_not_date = fund_pay[fund_ncpn - 1];
            if (fin_not_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fin_not_date, yc) * fix_not[fund_leg->num_cpn - 1];
                fix_leg_pv += temp;
            }
        }

        /*	Initial */
        if (start_date >= today + eod_pay_flag)
        {
            fix_leg_pv -= swp_f_df(today, start_date, yc) * fix_not[0];
        }
    }

    /*	4)	If there is at least one call after today, value call feature */

    if (call_feat == 1)
    {
        smessage("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);

        err = am_launch_adi_ts(pdLambdaValue, pdLambdaTime, nLambdaSize, am, und, adi_arg, &call);

        if (err)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        call = 0.0;
    }

    *fund_val = fund_leg_pv;
    *fix_val  = fix_leg_pv;
    *call_val = call;

    if (export_ts)
    {
        am_copy_und(und, und_exp);
    }

FREE_RETURN:

    if (free_struct)
    {
        am_free_all_struct(am, und, call_feat, adi_arg);
    }
    if (am)
        free(am);
    if (und)
        free(und);
    if (adi_arg)
        free(adi_arg);

    if (fund_mrg2)
        free(fund_mrg2);
    if (fund_spr2)
        free(fund_spr2);
    if (fund_fix_cpn2)
        free(fund_fix_cpn2);
    if (fund_not_ts2)
        free(fund_not_ts2);

    return err;
}