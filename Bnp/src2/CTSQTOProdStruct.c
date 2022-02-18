
#include "ctsqtoprodstruct.h"

#include "DoubleLGM1FQuanto_pde.h"
#include "fundlegprodstruct.h"
#include "lgmquantound.h"
#include "math.h"
#include "opfnctns.h"
#include "rangeaccrualprodstruct.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define BIG_BIG_NUMBER 1E40

Err merge_all_term_structure(
    int      num_sig,
    double*  sig_time,
    double*  dom_sig,
    double*  for_sig,
    int      num_fx_vol,
    double*  fx_vol_time,
    double*  fx_vol,
    int      corr_n_times,
    double*  corr_times,
    double*  correl_dom_for,
    double*  correl_dom_fx,
    double*  correl_for_fx,
    int*     sigma_n,
    double** sigma_time,
    double** dom_sigma,
    double** for_sigma,
    double** fx_sigma,
    double** domfor_rho,
    double** quanto_rho)
{
    int i, j;
    int IR_index, FX_index, RHO_index;
    Err err = NULL;
    int nb_merge_time;

    *sigma_time = (double*)calloc(num_sig, sizeof(double));
    if (!(*sigma_time))
    {
        smessage("memory allocation failed in merge_all_term_structure");
        err = "memory allocation failed in merge_all_term_structure";
        goto FREE_RETURN;
    }
    memcpy(*sigma_time, sig_time, num_sig * sizeof(double));
    nb_merge_time = num_sig;
    num_f_concat_vector(&nb_merge_time, sigma_time, num_fx_vol, fx_vol_time);
    num_f_concat_vector(&nb_merge_time, sigma_time, corr_n_times, corr_times);
    num_f_sort_vector(nb_merge_time, *sigma_time);
    num_f_unique_vector(&nb_merge_time, *sigma_time);

    *sigma_n = nb_merge_time;

    *dom_sigma  = (double*)calloc(nb_merge_time, sizeof(double));
    *for_sigma  = (double*)calloc(nb_merge_time, sizeof(double));
    *fx_sigma   = (double*)calloc(nb_merge_time, sizeof(double));
    *domfor_rho = (double*)calloc(nb_merge_time, sizeof(double));
    *quanto_rho = (double*)calloc(nb_merge_time, sizeof(double));
    if ((!dom_sigma) || (!for_sigma) || (!fx_sigma) || (!domfor_rho) || (!quanto_rho))
    {
        smessage("memory allocation failed in merge_all_term_structure");
        err = "Memory Allocation Failed in merge_all_term_structure";
        goto FREE_RETURN;
    }

    IR_index  = 0;
    FX_index  = 0;
    RHO_index = 0;

    smessage("first loop start in merge_all_term_structure");
    i = 0;
    for (IR_index = 0; IR_index < num_sig; ++IR_index)
    {
        while (((*sigma_time)[i] <= sig_time[IR_index]) && (i < nb_merge_time))
        {
            (*dom_sigma)[i] = dom_sig[IR_index];
            (*for_sigma)[i] = for_sig[IR_index];
            i++;
        }
        if (i == nb_merge_time)
        {
            IR_index = num_sig;
        }
    }
    for (j = i; j < nb_merge_time; ++j)
    {
        (*dom_sigma)[j] = dom_sig[num_sig - 1];
        (*for_sigma)[j] = for_sig[num_sig - 1];
    }

    smessage("second loop start in merge_all_term_structure");
    i = 0;
    for (FX_index = 0; FX_index < num_fx_vol; ++FX_index)
    {
        while (((*sigma_time)[i] <= fx_vol_time[FX_index]) && (i < nb_merge_time))
        {
            (*fx_sigma)[i] = fx_vol[FX_index];
            i++;
        }
        if (i == nb_merge_time)
        {
            FX_index = num_fx_vol;
        }
    }
    for (j = i; j < nb_merge_time; ++j)
    {
        (*fx_sigma)[j] = fx_vol[num_fx_vol - 1];
    }

    smessage("third loop start in merge_all_term_structure");
    i = 0;
    for (RHO_index = 0; RHO_index < corr_n_times; ++RHO_index)
    {
        while (((*sigma_time)[i] <= corr_times[RHO_index]) && (i < nb_merge_time))
        {
            (*domfor_rho)[i] = correl_dom_for[RHO_index];
            (*quanto_rho)[i] = correl_for_fx[RHO_index];
            i++;
        }
        if (i == nb_merge_time)
        {
            RHO_index = corr_n_times;
        }
    }
    for (j = i; j < nb_merge_time; ++j)
    {
        (*domfor_rho)[j] = correl_dom_for[corr_n_times - 1];
        (*quanto_rho)[j] = correl_for_fx[corr_n_times - 1];
    }

    smessage("free start in merge_all_term_structure");

FREE_RETURN:

    return err;
}

/*	Fill LGM Quanto underlying structure from CTSQTO calibration instruments */
Err ctsqto_calib_und(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I, 1: E */

    char* dom_yc,         /*	yc */
    char* dom_vc,         /*	vc (only if calib) */
    char* dom_ref,        /*	ref rate (only if calib) */
    char* dom_swap_freq,  /*	swap freq (only if calib) */
    char* dom_swap_basis, /*	swap basis (only if calib) */

    char* for_yc,         /*	yc */
    char* for_vc,         /*	vc (only if calib) */
    char* for_ref,        /*	ref rate (only if calib) */
    char* for_swap_freq,  /*	swap freq (only if calib) */
    char* for_swap_basis, /*	swap basis (only if calib) */

    int forcalib, /*	0 : RA Und, 1 : Diag */

    double dom_lambda, /*	unique lambda */
    double for_lambda, /*	unique lambda */
    /*	Calib params */
    int    dom_force_atm, /*	force atm calib domestic und */
    int    for_force_atm, /*	force atm calib foreign und */
    double max_std_long,
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
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,

    double* corr_times,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    int     corr_n_times,

    /*	Diag Calib Parameters */
    CPD_DIAG_CALIB_PARAM calibparam,

    /*	End of calib params */
    CTSQTO_STR ctsqto,  /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    LGMQTO_UND und)
{
    int   i, nex;
    long* lex    = NULL;
    int*  uselex = NULL;

    long    end_struct, end_swap;
    double *dom_strikes = NULL, *for_strikes = NULL;

    int strike_type;

    RangeAccrualStruct* RA;

    Err err = NULL;

    char** dom_end_tenor = NULL;
    char** for_end_tenor = NULL;

    char* dom_swap_freq_loc;
    char* for_swap_freq_loc;

    /*	Calibration instrument data */
    //	CPD_CALIB_INST_DATA	dom_inst_data;
    //	CPD_CALIB_INST_DATA	for_inst_data;

    int numMonth;

    int     num_sig;
    double* sig_time = NULL;
    double* dom_sig  = NULL;
    double* for_sig  = NULL;
    int     num_fx_vol;
    double* fx_vol_time = NULL;
    double* fx_vol      = NULL;
    double* domfx_rho   = NULL;

    dom_swap_freq_loc = (char*)calloc(256, sizeof(char));
    for_swap_freq_loc = (char*)calloc(256, sizeof(char));
    if ((!dom_swap_freq_loc) || (!for_swap_freq_loc))
    {
        smessage("Memory Allocation failed in ctsqto_calib_und");
        err = "Memory Allocation failed in ctsqto_calib_und";
        goto FREE_RETURN;
    }

    strcpy(dom_swap_freq_loc, dom_swap_freq);
    strcpy(for_swap_freq_loc, for_swap_freq);

    if (fabs(dom_lambda) < 1.0e-08)
    {
        dom_lambda = 1.0e-08;
    }

    if (fabs(for_lambda) < 1.0e-08)
    {
        for_lambda = 1.0e-08;
    }

    /*	Initialise */
    RA = ctsqto->ra_leg;

    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->dom_sigma  = NULL;
    und->for_sigma  = NULL;
    und->fx_sigma   = NULL;
    und->domfor_rho = NULL;
    und->quanto_rho = NULL;

    und->has_fwd_iv    = 0;
    und->nb_fwdiv      = 0;
    und->exercise_date = NULL;
    und->market_fwdiv  = NULL;
    und->model_fwdiv   = NULL;
    und->extra_fees    = NULL;

    und->has_inst_data = 0;
    //	cpd_init_calib_inst_data (&(und->dom_inst_data));
    //	cpd_init_calib_inst_data (&(und->for_inst_data));
    und->today = today;
    strcpy(und->dom_yc, dom_yc);
    strcpy(und->dom_vc, dom_vc);
    strcpy(und->dom_ref, dom_ref);
    strcpy(und->dom_swap_freq, dom_swap_freq);
    strcpy(und->dom_swap_basis, dom_swap_basis);

    strcpy(und->for_yc, for_yc);
    strcpy(und->for_vc, for_vc);
    strcpy(und->for_ref, for_ref);
    strcpy(und->for_swap_freq, for_swap_freq);
    strcpy(und->for_swap_basis, for_swap_basis);

    strcpy(und->name, "CALIB");

    /*	Exercise dates for calibration */

    end_struct = RA->dates[RA->n_periods];
    end_swap   = end_struct;

    if (ctsqto->num_calls > 0 &&
        !(ctsqto->num_calls == 1 && ctsqto->call[0].ex_date <= today + eod_flag))
    {
        /*	If call dates, choose call dates as option expiries for calibration */
        nex = ctsqto->num_calls;
        lex = (long*)calloc(nex, sizeof(long));
        if (!lex)
        {
            err = "Allocation error (2) in ctsqto_calib_und";
            goto FREE_RETURN;
        }

        uselex = (int*)calloc(nex, sizeof(int));
        if (!uselex)
        {
            err = "Allocation error (3) in ctsqto_calib_und";
            goto FREE_RETURN;
        }

        dom_end_tenor = (char**)calloc(nex, sizeof(char*));
        for_end_tenor = (char**)calloc(nex, sizeof(char*));
        if ((!dom_end_tenor) || (!for_end_tenor))
        {
            err = "Memory allocation error in ctsqto calib und";
            goto FREE_RETURN;
        }

        for (i = 0; i < nex; i++)
        {
            dom_end_tenor[i] = (char*)calloc(256, sizeof(char));
            for_end_tenor[i] = (char*)calloc(256, sizeof(char));
            if ((!dom_end_tenor[i]) || (!for_end_tenor[i]))
            {
                err = "Memory allocation error in ctsqto calib und";
                goto FREE_RETURN;
            }
            lex[i]    = ctsqto->call[i].ex_date;
            uselex[i] = 1;
            strcpy(dom_end_tenor[i], "DIAG");

            if (forcalib == 0)
            {
                strcpy(for_end_tenor[i], RA->tenor_f_char);
            }
            else
            {
                strcpy(for_end_tenor[i], "DIAG");
            }

            if ((end_swap - ctsqto->call[i].ex_date) / 365 < 1)
            {
                //				numMonth = (int)((end_swap -
                //ctsqto->call[i].ex_date)/365.0 * 12);
                numMonth = 12;
                sprintf(dom_end_tenor[i], "%dM", numMonth);
                if (forcalib == 0)
                {
                    strcpy(for_end_tenor[i], RA->tenor_f_char);
                }
                else
                {
                    sprintf(for_end_tenor[i], "%dM", numMonth);
                }
            }
        }
    }
    else
    {
        /*	If no call dates, exit */
        ctsqto->num_calls = 0;
        goto FREE_RETURN;
    }

    /*	Implement force atm */
    if (dom_force_atm)
    {
        strike_type = 0;
        dom_strikes = NULL;
        dom_strikes = (double*)calloc(nex, sizeof(double));
        if (!dom_strikes)
        {
            err = "Allocation error in ctsqto calib und";
            goto FREE_RETURN;
        }
        smessage("FORCE ATM flag detected - calibrating ATM");
    }
    else
    {
        /*	Compute strikes */
        strike_type = 0;
        dom_strikes = NULL;
        dom_strikes = (double*)calloc(nex, sizeof(double));
        if (!dom_strikes)
        {
            err = "Allocation error (3) in ctsqto_calib_und";
            goto FREE_RETURN;
        }
    }

    /*	Implement force atm */
    if (for_force_atm)
    {
        strike_type = 0;
        for_strikes = NULL;
        for_strikes = (double*)calloc(nex, sizeof(double));
        if (!for_strikes)
        {
            err = "Allocation error (3) in ctsqto_calib_und";
            goto FREE_RETURN;
        }
        smessage("FORCE ATM flag detected - calibrating ATM");
    }
    else
    {
        /*	Compute strikes */
        for_strikes = (double*)calloc(nex, sizeof(double));
        if (!for_strikes)
        {
            err = "Allocation error (3) in ctsqto_calib_und";
            goto FREE_RETURN;
        }
    }

    /*	Calibrate all */
    err = cpd_calib_all(
        today,
        get_cash_vol,
        dom_yc,
        dom_vc,
        dom_ref,
        dom_swap_freq_loc,
        dom_swap_basis,
        dom_lambda,
        for_yc,
        for_vc,
        for_ref,
        for_swap_freq,
        for_swap_basis,
        for_lambda,
        fx_mkt_vol_date,
        fx_mkt_vol,
        num_fx_mkt_vol,
        corr_times,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        corr_n_times,
        nex,
        lex,
        uselex,
        dom_end_tenor,
        for_end_tenor,
        end_swap,
        dom_strikes,
        for_strikes,
        &num_sig,
        &sig_time,
        &dom_sig,
        &for_sig,
        &num_fx_vol,
        &fx_vol_time,
        &fx_vol,
        calibparam,
        NULL,   //&(und->dom_inst_data),
        NULL);  //&(und->dom_inst_data));

    smessage("Calib all OK");
    if (err)
    {
        smessage("Error in cpd_calib_all");
        goto FREE_RETURN;
    }

    //	err = merge_all_term_structure(
    //		num_sig, sig_time, dom_sig, for_sig,
    //		num_fx_vol, fx_vol_time, fx_vol,
    //		corr_n_times, corr_times, correl_dom_for, correl_dom_fx, correl_for_fx,
    //		&(und->sigma_n), &(und->sigma_time), &(und->dom_sigma), &(und->for_sigma),
    //		&(und->fx_sigma), &(und->domfor_rho), &(und->quanto_rho));

    err = merge_rates_fx_corr_ts(
        sig_time,
        dom_sig,
        num_sig,
        sig_time,
        for_sig,
        num_sig,
        fx_vol_time,
        fx_vol,
        num_fx_vol,
        corr_times,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        corr_n_times,
        &(und->sigma_time),
        &(und->dom_sigma),
        &(und->for_sigma),
        &(und->fx_sigma),
        &(und->domfor_rho),
        &domfx_rho,
        &(und->quanto_rho),
        &(und->sigma_n));

    smessage("merge_all_term_structure OK");
    if (err)
    {
        smessage("Error in merge_all_term_structure");
        goto FREE_RETURN;
    }

    und->has_inst_data = 1;

    und->dom_lambda = dom_lambda;
    und->for_lambda = for_lambda;

FREE_RETURN:

    smessage("Free in ctsqto_calib_und 1");
    if (err)
    {
        smessage("Error in ctsqto calib und");
    }

    smessage("Free in ctsqto_calib_und 2");
    for (i = 0; i < nex; i++)
    {
        if (dom_end_tenor[i])
            free(dom_end_tenor[i]);
        if (for_end_tenor[i])
            free(for_end_tenor[i]);
    }
    if (dom_end_tenor)
        free(dom_end_tenor);
    if (for_end_tenor)
        free(for_end_tenor);

    smessage("Free in ctsqto_calib_und 2");
    if (err)
        free_lgmQto_und(und);

    if (dom_swap_freq_loc)
        free(dom_swap_freq_loc);
    if (for_swap_freq_loc)
        free(for_swap_freq_loc);

    smessage("Free in ctsqto_calib_und 3");
    if (lex)
        free(lex);
    if (uselex)
        free(uselex);
    if (dom_strikes)
        free(dom_strikes);
    if (for_strikes)
        free(for_strikes);

    smessage("Free in ctsqto_calib_und 4");
    if (sig_time)
        free(sig_time);
    if (dom_sig)
        free(dom_sig);
    if (for_sig)
        free(for_sig);
    if (fx_vol_time)
        free(fx_vol_time);
    if (fx_vol)
        free(fx_vol);
    if (domfx_rho)
        free(domfx_rho);

    return err;
}

/*	Functions for the calls */

/*	Check dates consistency */
Err ctsqto_check_calls(CTSQTO_STR ctsqto)
{
    int i;

    /*	Check that ex and set dates are strictly increasing
                    Also check that funding and pd indices are strictly increasing,
                    i.e. there is no redundant calls */
    for (i = 1; i < ctsqto->num_calls; i++)
    {
        if (ctsqto->call[i].ex_date <= ctsqto->call[i - 1].ex_date)
        {
            return "Exercise dates should be increasing";
        }

        if (ctsqto->call[i].set_date <= ctsqto->call[i - 1].set_date)
        {
            return "Settlement dates should be increasing";
        }
    }

    /*	Check that set dates are after ex dates
                    Also check that the call date is before the start, end and fixing dates
                    of the coupons it controls */
    for (i = 0; i < ctsqto->num_calls; i++)
    {
        if (ctsqto->call[i].set_date < ctsqto->call[i].ex_date)
        {
            return "Settlement dates should be after exercise dates";
        }
    }

    /*	OK */
    return NULL;
}

/*	Free */
Err ctsqto_free_calls(CTSQTO_STR ctsqto)
{
    if (ctsqto->call)
    {
        free(ctsqto->call);
        ctsqto->call = NULL;
    }

    return NULL;
}

/*	Fill */
Err ctsqto_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int        eod_flag, /*	0: I, 1: E */
    int        ncall,
    int        pay_rec, /*	0: rec pd, 1: pay pd */
    long*      ex_date,
    long*      set_date,
    double*    fee,
    CTSQTO_STR ctsqto)
{
    int i, j, k;

    RangeAccrualStruct* RA;
    CTSQTO_CALL         call;
    FUNDING_LEG         fund_leg;
    Err                 err = NULL;

    /*	Initialise pointers */
    ctsqto->call = NULL;
    RA           = ctsqto->ra_leg;
    fund_leg     = ctsqto->fund_leg;

    /*	Skip calls to be exercised before today */
    i = 0;
    while (i < ncall && ex_date[i] < today + eod_flag)
    {
        i++;
    }

    /*	Check that at least one call is left */
    if (i == ncall)
    {
        /* err = "All calls are to be exercised before today in ctsqto_fill_calls"; */
        ctsqto->num_calls = 0;
        goto FREE_RETURN;
    }

    /*	Allocate memory */
    ctsqto->num_calls = ncall - i;
    ctsqto->call      = (ctsqto_call*)calloc(ctsqto->num_calls, sizeof(ctsqto_call));
    if (!ctsqto->num_calls)
    {
        err = "Allocation error in ctsqto_fill_calls";
        goto FREE_RETURN;
    }

    /*	Fill coupons information */
    j = 0;

    while (i < ncall)
    {
        call = ctsqto->call + j;

        /*	Dates */
        call->ex_date  = ex_date[i];
        call->set_date = set_date[i];

        /*	Times */
        call->ex_time  = (ctsqto->call[j].ex_date - today) * YEARS_IN_DAY;
        call->set_time = (ctsqto->call[j].set_date - today) * YEARS_IN_DAY;

        /*	Call on funding leg */
        /*	k = index of the first coupon to be called on funding leg,
                        i.e. first coupon with a start date >= ex date */

        k = 0;
        while (k < fund_leg->num_cpn && fund_leg->cpn[k].start_date < call->ex_date)
        {
            k++;
        }
        if (k == fund_leg->num_cpn)
        {
            serror("Call number %d does not control any coupon in funding leg", i);
            err = "One call does not control any coupon in funding leg";
            goto FREE_RETURN;
        }
        call->fund_idx     = k;
        call->num_fund_cpn = fund_leg->num_cpn - k;

        /*	Call on exotic leg */
        /*	k = index of the first coupon to be called on exotic leg,
                        i.e. first coupon with a start date >= ex date */
        k = 0;
        while (k < RA->n_periods && RA->dates[k] < call->ex_date)
        {
            k++;
        }
        if (k == RA->n_periods)
        {
            serror("Call number %d does not control any coupon in exotic leg", i);
            err = "One call does not control any coupon in exotic leg";
            goto FREE_RETURN;
        }
        call->ra_idx     = k;
        call->num_ra_cpn = RA->n_periods - k;

        /*	Payer or receiver */
        call->pay_rec = pay_rec;

        /*	Fee */
        call->fee = fee[i];

        i++;
        j++;
    }

    err = ctsqto_check_calls(ctsqto);

FREE_RETURN:

    if (err)
    {
        ctsqto_free_calls(ctsqto);
    }

    return err;
}

Err ctsqto_payoff_adi(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* dom_yc,
    void* for_yc,

    /* Model data	*/
    double domlam,
    double forlam,
    double domphi,
    double forphi,

    /* Grid data	*/
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,

    /* Vector of results to be updated */
    int       nprod,
    double*** prod_val)
{
    CTSQTO_STR          ctsqtostr;
    RangeAccrualStruct* RA       = NULL;
    FUNDING_LEG         fund_leg = NULL;
    //	ctsqto_call	call;

    double *dom_dff = NULL, *dom_gam = NULL, *dom_gam2 = NULL, *for_dff = NULL, *for_gam = NULL,
           *for_gam2 = NULL, *cpn_float_dff = NULL, *cpn_float_gam = NULL, *cpn_float_gam2 = NULL;

    double *fund_dom_dff = NULL, *fund_dom_gam = NULL, *fund_dom_gam2 = NULL;

    int            i, j, l;
    unsigned short do_dom, do_for;
    Err            err       = NULL;
    int            pn_unds   = 0;
    int*           pn_dfs    = NULL;
    long**         ra_pdates = NULL;
    double**       ra_dfs    = NULL;

    int     fund_ndfs;
    long*   fund_pdates = NULL;
    double* fund_dfs    = NULL;

    double time_i, cpn_float_time_i;

    double ra_payoff;
    double fund_payoff;
    double payoff;

    double fee;

    int    call_idx;
    double fee_dom_df, fee_dom_dff, fee_dom_gam, fee_dom_gam2;
    double time_call;

    /* Get the event */
    ctsqtostr = (CTSQTO_STR)(func_parm);
    RA        = (RangeAccrualStruct*)(ctsqtostr->ra_leg);
    fund_leg  = (FUNDING_LEG)(ctsqtostr->fund_leg);

    i = 0;
    //	call = (ctsqtostr->call[i]);

    while ((evt_date > (ctsqtostr->call[i]).ex_date) && (i < ctsqtostr->num_calls))
    {
        ++i;
    }
    call_idx = i;
    if (call_idx >= ctsqtostr->num_calls)
    {
        err = "All exercise dates are in the past (ctsqto_payoff_adi)";
    }

    err = RA_RequestDfDates(RA, (long)(evt_date), &pn_unds, &pn_dfs, &ra_pdates);

    err = FundLeg_RequestDfDates(fund_leg, (long)(evt_date), &fund_ndfs, &fund_pdates);

    fund_dfs = (double*)calloc(fund_ndfs, sizeof(double));
    if (!fund_dfs)
    {
        err = "Memory allocation error (1) in ctsqto_payoff";
        goto FREE_RETURN;
    }
    ra_dfs = (double**)calloc(pn_unds, sizeof(double*));
    if (!ra_dfs)
    {
        err = "Memory allocation error (1) in ctsqto_payoff";
        goto FREE_RETURN;
    }
    ra_dfs[0] = (double*)calloc(pn_dfs[0], sizeof(double));
    if (!ra_dfs[0])
    {
        err = "Memory allocation error (1) in ctsqto_payoff";
        goto FREE_RETURN;
    }
    ra_dfs[1] = (double*)calloc(pn_dfs[1], sizeof(double));
    if (!ra_dfs[1])
    {
        err = "Memory allocation error (1) in ctsqto_payoff";
        goto FREE_RETURN;
    }

    // In case of floating coupons
    if (RA->cpn_type != 0)
    {
        ra_dfs[2] = (double*)calloc(2 * pn_dfs[2], sizeof(double));
        if (!ra_dfs[2])
        {
            err = "Memory allocation error (1) in ctsqto_payoff";
            goto FREE_RETURN;
        }
    }

    /* FEE : Precalculate DF, gamma and 0.5 * gamma * gamma */
    fee_dom_dff  = swp_f_df(evt_date, (ctsqtostr->call[call_idx]).ex_date, (char*)dom_yc);
    time_call    = ((ctsqtostr->call[call_idx]).ex_date - evt_date) * YEARS_IN_DAY;
    fee_dom_gam  = (1.0 - exp(-domlam * time_call)) / domlam;
    fee_dom_gam2 = 0.5 * fee_dom_gam * fee_dom_gam * domphi;

    /* FUNDING LEG : Precalculate DF, gamma and 0.5 * gamma * gamma */
    fund_dom_dff  = dvector(0, fund_ndfs - 1);
    fund_dom_gam  = dvector(0, fund_ndfs - 1);
    fund_dom_gam2 = dvector(0, fund_ndfs - 1);

    if (!fund_dom_dff || !fund_dom_gam || !fund_dom_gam2)
    {
        err = "Memory allocation error (1) in ctsqto_payoff";
        goto FREE_RETURN;
    }

    for (i = 0; i < fund_ndfs; i++)
    {
        fund_dom_dff[i]  = swp_f_df(evt_date, fund_pdates[i], (char*)dom_yc);
        time_i           = (fund_pdates[i] - evt_date) * YEARS_IN_DAY;
        fund_dom_gam[i]  = (1.0 - exp(-domlam * time_i)) / domlam;
        fund_dom_gam2[i] = 0.5 * fund_dom_gam[i] * fund_dom_gam[i] * domphi;
    }

    /* RA LEG : Precalculate DF, gamma and 0.5 * gamma * gamma */
    if (pn_dfs[0] > 0)
    {
        dom_dff  = dvector(0, pn_dfs[0] - 1);
        dom_gam  = dvector(0, pn_dfs[0] - 1);
        dom_gam2 = dvector(0, pn_dfs[0] - 1);

        if (!dom_dff || !dom_gam || !dom_gam2)
        {
            err = "Memory allocation error (1) in ctsqto_payoff";
            goto FREE_RETURN;
        }

        for (i = 0; i < pn_dfs[0]; i++)
        {
            dom_dff[i]  = swp_f_df(evt_date, ra_pdates[0][i], (char*)dom_yc);
            time_i      = (ra_pdates[0][i] - evt_date) * YEARS_IN_DAY;
            dom_gam[i]  = (1.0 - exp(-domlam * time_i)) / domlam;
            dom_gam2[i] = 0.5 * dom_gam[i] * dom_gam[i] * domphi;
        }

        do_dom = 1;
    }
    else
    {
        do_dom = 0;
    }

    if (pn_dfs[1] > 0)
    {
        for_dff  = dvector(0, pn_dfs[1] - 1);
        for_gam  = dvector(0, pn_dfs[1] - 1);
        for_gam2 = dvector(0, pn_dfs[1] - 1);

        if (!for_dff || !for_gam || !for_gam2)
        {
            err = "Memory allocation error (3) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < pn_dfs[1]; i++)
        {
            for_dff[i]  = swp_f_df(evt_date, ra_pdates[1][i], (char*)for_yc);
            time_i      = (ra_pdates[1][i] - evt_date) * YEARS_IN_DAY;
            for_gam[i]  = (1.0 - exp(-forlam * time_i)) / forlam;
            for_gam2[i] = 0.5 * for_gam[i] * for_gam[i] * forphi;
        }

        do_for = 1;
    }
    else
    {
        do_for = 0;
    }

    // In case of floating coupons
    if (RA->cpn_type != 0)
    {
        cpn_float_dff  = dvector(0, 2 * pn_dfs[2] - 1);
        cpn_float_gam  = dvector(0, 2 * pn_dfs[2] - 1);
        cpn_float_gam2 = dvector(0, 2 * pn_dfs[2] - 1);

        if (!cpn_float_dff || !cpn_float_gam || !cpn_float_gam2)
        {
            err = "Memory allocation error (4) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < 2 * pn_dfs[2]; i++)
        {
            if (RA->float_cpn_is_dom_for == 0)
            {
                cpn_float_dff[i]  = swp_f_df(evt_date, ra_pdates[2][i], (char*)dom_yc);
                cpn_float_time_i  = (ra_pdates[2][i] - evt_date) * YEARS_IN_DAY;
                cpn_float_gam[i]  = (1.0 - exp(-domlam * cpn_float_time_i)) / domlam;
                cpn_float_gam2[i] = 0.5 * cpn_float_gam[i] * cpn_float_gam[i] * domphi;
            }
            else
            {
                cpn_float_dff[i]  = swp_f_df(evt_date, ra_pdates[2][i], (char*)for_yc);
                cpn_float_time_i  = (ra_pdates[2][i] - evt_date) * YEARS_IN_DAY;
                cpn_float_gam[i]  = (1.0 - exp(-forlam * cpn_float_time_i)) / forlam;
                cpn_float_gam2[i] = 0.5 * cpn_float_gam[i] * cpn_float_gam[i] * forphi;
            }
        }
    }

    /* Eval payoff */
    for (i = l1; i < u1; i++)
    {
        for (j = l2; j < u2; j++)
        {
            //-------------- FUND LEG Dfs --------------------------
            for (l = 0; l < fund_ndfs; l++)
            {
                fund_dfs[l] = fund_dom_dff[l] * exp(-fund_dom_gam[l] * r1[i] - fund_dom_gam2[l]);
            }

            FundLeg_Payoff(fund_leg, (long)(evt_date), (long)(evt_date), fund_dfs, &fund_payoff);

            //-------------- RA LEG Dfs ----------------------------
            if (do_dom)
            {
                for (l = 0; l < pn_dfs[0]; l++)
                {
                    ra_dfs[0][l] = dom_dff[l] * exp(-dom_gam[l] * r1[i] - dom_gam2[l]);
                }
            }
            if (do_for)
            {
                for (l = 0; l < pn_dfs[1]; l++)
                {
                    ra_dfs[1][l] = for_dff[l] * exp(-for_gam[l] * r2[i][j] - for_gam2[l]);
                }
            }

            // DFs in order to calculate the floating coupons of the RA leg
            if (RA->cpn_type != 0)
            {
                if (RA->float_cpn_is_dom_for == 0)
                {
                    for (l = 0; l < 2 * pn_dfs[2]; l++)
                    {
                        ra_dfs[2][l] =
                            cpn_float_dff[l] * exp(-cpn_float_gam[l] * r1[i] - cpn_float_gam2[l]);
                    }
                }
                else
                {
                    for (l = 0; l < 2 * pn_dfs[2]; l++)
                    {
                        ra_dfs[2][l] = cpn_float_dff[l] *
                                       exp(-cpn_float_gam[l] * r2[i][j] - cpn_float_gam2[l]);
                    }
                }
            }

            RA_Payoff(RA, (long)(evt_date), (long)(evt_date), ra_dfs, &ra_payoff);

            //-------------- FEE Df ----------------------------
            fee_dom_df = fee_dom_dff * exp(-fee_dom_gam * r1[i] - fee_dom_gam2);
            if (fabs((ctsqtostr->call[call_idx]).fee) > 1.0e-08)
            {
                fee = fee_dom_df * (ctsqtostr->call[call_idx]).fee;
            }
            else
            {
                fee = 0.0;
            }

            //			payoff = ra_payoff + fund_payoff - fee;
            // Was + fee before, changed to - fee for consistency reasons
            if ((ctsqtostr->call[call_idx]).pay_rec == 0)
            {
                payoff = ra_payoff - fund_payoff - fee;
            }
            else
            {
                payoff = -ra_payoff + fund_payoff - fee;
            }

            if (prod_val[i][j][0] < DMAX(0.00, payoff))
            {
                prod_val[i][j][0] = DMAX(0.00, payoff);
            }

            if (err)
            {
                goto FREE_RETURN;
            }
        }
    }

FREE_RETURN:

    if (ra_dfs[0])
        free(ra_dfs[0]);
    if (ra_dfs[1])
        free(ra_dfs[1]);
    if (RA->cpn_type != 0)
    {
        if (ra_dfs[2])
            free(ra_dfs[2]);
    }
    if (ra_dfs)
        free(ra_dfs);

    if (ra_pdates[0])
        free(ra_pdates[0]);
    if (ra_pdates[1])
        free(ra_pdates[1]);
    if (RA->cpn_type != 0)
    {
        if (ra_pdates[2])
            free(ra_pdates[2]);
    }
    if (ra_pdates)
        free(ra_pdates);

    if (dom_dff)
        free_dvector(dom_dff, 0, pn_dfs[0] - 1);
    if (dom_gam)
        free_dvector(dom_gam, 0, pn_dfs[0] - 1);
    if (dom_gam2)
        free_dvector(dom_gam2, 0, pn_dfs[0] - 1);

    if (for_dff)
        free_dvector(for_dff, 0, pn_dfs[1] - 1);
    if (for_gam)
        free_dvector(for_gam, 0, pn_dfs[1] - 1);
    if (for_gam2)
        free_dvector(for_gam2, 0, pn_dfs[1] - 1);

    if (cpn_float_dff)
        free_dvector(cpn_float_dff, 0, 2 * pn_dfs[2] - 1);
    if (cpn_float_gam)
        free_dvector(cpn_float_gam, 0, 2 * pn_dfs[2] - 1);
    if (cpn_float_gam2)
        free_dvector(cpn_float_gam2, 0, 2 * pn_dfs[2] - 1);

    if (fund_dom_dff)
        free_dvector(fund_dom_dff, 0, fund_ndfs - 1);
    if (fund_dom_gam)
        free_dvector(fund_dom_gam, 0, fund_ndfs - 1);
    if (fund_dom_gam2)
        free_dvector(fund_dom_gam2, 0, fund_ndfs - 1);

    if (pn_dfs)
        free(pn_dfs);

    if (fund_dfs)
        free(fund_dfs);
    if (fund_pdates)
        free(fund_pdates);

    return err;
}

Err ctsqto_payoff_adi2(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* dom_yc,
    void* for_yc,

    /* Model data	*/
    double domlam,
    double forlam,
    double domphi,
    double forphi,

    /* Grid data	*/
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,

    /* Vector of results to be updated */
    int       nprod,
    double*** prod_val)
{
    CTSQTO_STR          ctsqtostr;
    RangeAccrualStruct* RA       = NULL;
    FUNDING_LEG         fund_leg = NULL;
    //	ctsqto_call	call;

    double *dom_dff = NULL, *dom_gam = NULL, *dom_gam2 = NULL, *for_dff = NULL, *for_gam = NULL,
           *for_gam2 = NULL, *cpn_float_dff = NULL, *cpn_float_gam = NULL, *cpn_float_gam2 = NULL;

    double *fund_dom_dff = NULL, *fund_dom_gam = NULL, *fund_dom_gam2 = NULL;

    int            i, j, l;
    unsigned short do_dom, do_for;
    Err            err       = NULL;
    int            pn_unds   = 0;
    int*           pn_dfs    = NULL;
    long**         ra_pdates = NULL;
    double**       ra_dfs    = NULL;

    int     fund_ndfs;
    long*   fund_pdates = NULL;
    double* fund_dfs    = NULL;

    double time_i, cpn_float_time_i;

    double ra_payoff;
    double fund_payoff;
    double payoff;

    double fee;

    int    call_idx;
    double fee_dom_df, fee_dom_dff, fee_dom_gam, fee_dom_gam2;
    double time_call;

    /* Get the event */
    ctsqtostr = (CTSQTO_STR)(func_parm);
    RA        = (RangeAccrualStruct*)(ctsqtostr->ra_leg);
    fund_leg  = (FUNDING_LEG)(ctsqtostr->fund_leg);

    i = 0;
    //	call = (ctsqtostr->call[i]);

    while ((evt_date > (ctsqtostr->call[i]).ex_date) && (i < ctsqtostr->num_calls))
    {
        ++i;
    }
    call_idx = i;
    if (call_idx >= ctsqtostr->num_calls)
    {
        err = "All exercise dates are in the past (ctsqto_payoff_adi)";
    }

    err = RA_RequestDfDates(RA, (long)(evt_date), &pn_unds, &pn_dfs, &ra_pdates);

    err = FundLeg_RequestDfDates(fund_leg, (long)(evt_date), &fund_ndfs, &fund_pdates);

    fund_dfs = (double*)calloc(fund_ndfs, sizeof(double));
    if (!fund_dfs)
    {
        err = "Memory allocation error (1) in ctsqto_payoff2";
        goto FREE_RETURN;
    }
    ra_dfs = (double**)calloc(pn_unds, sizeof(double*));
    if (!ra_dfs)
    {
        err = "Memory allocation error (1) in ctsqto_payoff2";
        goto FREE_RETURN;
    }
    ra_dfs[0] = (double*)calloc(pn_dfs[0], sizeof(double));
    if (!ra_dfs[0])
    {
        err = "Memory allocation error (1) in ctsqto_payoff2";
        goto FREE_RETURN;
    }
    ra_dfs[1] = (double*)calloc(pn_dfs[1], sizeof(double));
    if (!ra_dfs[1])
    {
        err = "Memory allocation error (1) in ctsqto_payoff2";
        goto FREE_RETURN;
    }

    // In case of floating coupons
    if (RA->cpn_type != 0)
    {
        ra_dfs[2] = (double*)calloc(2 * pn_dfs[2], sizeof(double));
        if (!ra_dfs[2])
        {
            err = "Memory allocation error (1) in ctsqto_payoff";
            goto FREE_RETURN;
        }
    }

    /* FEE : Precalculate DF, gamma and 0.5 * gamma * gamma */
    fee_dom_dff  = swp_f_df(evt_date, (ctsqtostr->call[call_idx]).ex_date, (char*)dom_yc);
    time_call    = ((ctsqtostr->call[call_idx]).ex_date - evt_date) * YEARS_IN_DAY;
    fee_dom_gam  = (1.0 - exp(-domlam * time_call)) / domlam;
    fee_dom_gam2 = 0.5 * fee_dom_gam * fee_dom_gam * domphi;

    /* FUNDING LEG : Precalculate DF, gamma and 0.5 * gamma * gamma */
    fund_dom_dff  = dvector(0, fund_ndfs - 1);
    fund_dom_gam  = dvector(0, fund_ndfs - 1);
    fund_dom_gam2 = dvector(0, fund_ndfs - 1);

    if (!fund_dom_dff || !fund_dom_gam || !fund_dom_gam2)
    {
        err = "Memory allocation error (1) in ctsqto_payoff";
        goto FREE_RETURN;
    }

    for (i = 0; i < fund_ndfs; i++)
    {
        fund_dom_dff[i]  = swp_f_df(evt_date, fund_pdates[i], (char*)dom_yc);
        time_i           = (fund_pdates[i] - evt_date) * YEARS_IN_DAY;
        fund_dom_gam[i]  = (1.0 - exp(-domlam * time_i)) / domlam;
        fund_dom_gam2[i] = 0.5 * fund_dom_gam[i] * fund_dom_gam[i] * domphi;
    }

    /* RA LEG : Precalculate DF, gamma and 0.5 * gamma * gamma */
    if (pn_dfs[0] > 0)
    {
        dom_dff  = dvector(0, pn_dfs[0] - 1);
        dom_gam  = dvector(0, pn_dfs[0] - 1);
        dom_gam2 = dvector(0, pn_dfs[0] - 1);

        if (!dom_dff || !dom_gam || !dom_gam2)
        {
            err = "Memory allocation error (1) in ctsqto_payoff";
            goto FREE_RETURN;
        }

        for (i = 0; i < pn_dfs[0]; i++)
        {
            dom_dff[i]  = swp_f_df(evt_date, ra_pdates[0][i], (char*)dom_yc);
            time_i      = (ra_pdates[0][i] - evt_date) * YEARS_IN_DAY;
            dom_gam[i]  = (1.0 - exp(-domlam * time_i)) / domlam;
            dom_gam2[i] = 0.5 * dom_gam[i] * dom_gam[i] * domphi;
        }

        do_dom = 1;
    }
    else
    {
        do_dom = 0;
    }

    if (pn_dfs[1] > 0)
    {
        for_dff  = dvector(0, pn_dfs[1] - 1);
        for_gam  = dvector(0, pn_dfs[1] - 1);
        for_gam2 = dvector(0, pn_dfs[1] - 1);

        if (!for_dff || !for_gam || !for_gam2)
        {
            err = "Memory allocation error (3) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < pn_dfs[1]; i++)
        {
            for_dff[i]  = swp_f_df(evt_date, ra_pdates[1][i], (char*)for_yc);
            time_i      = (ra_pdates[1][i] - evt_date) * YEARS_IN_DAY;
            for_gam[i]  = (1.0 - exp(-forlam * time_i)) / forlam;
            for_gam2[i] = 0.5 * for_gam[i] * for_gam[i] * forphi;
        }

        do_for = 1;
    }
    else
    {
        do_for = 0;
    }

    // In case of floating coupons
    if (RA->cpn_type != 0)
    {
        cpn_float_dff  = dvector(0, 2 * pn_dfs[2] - 1);
        cpn_float_gam  = dvector(0, 2 * pn_dfs[2] - 1);
        cpn_float_gam2 = dvector(0, 2 * pn_dfs[2] - 1);

        if (!cpn_float_dff || !cpn_float_gam || !cpn_float_gam2)
        {
            err = "Memory allocation error (4) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < 2 * pn_dfs[2]; i++)
        {
            if (RA->float_cpn_is_dom_for == 0)
            {
                cpn_float_dff[i]  = swp_f_df(evt_date, ra_pdates[2][i], (char*)dom_yc);
                cpn_float_time_i  = (ra_pdates[2][i] - evt_date) * YEARS_IN_DAY;
                cpn_float_gam[i]  = (1.0 - exp(-domlam * cpn_float_time_i)) / domlam;
                cpn_float_gam2[i] = 0.5 * cpn_float_gam[i] * cpn_float_gam[i] * domphi;
            }
            else
            {
                cpn_float_dff[i]  = swp_f_df(evt_date, ra_pdates[2][i], (char*)for_yc);
                cpn_float_time_i  = (ra_pdates[2][i] - evt_date) * YEARS_IN_DAY;
                cpn_float_gam[i]  = (1.0 - exp(-forlam * cpn_float_time_i)) / forlam;
                cpn_float_gam2[i] = 0.5 * cpn_float_gam[i] * cpn_float_gam[i] * forphi;
            }
        }
    }

    /* Eval payoff */
    for (i = l1; i < u1; i++)
    {
        for (j = l2; j < u2; j++)
        {
            //-------------- FUND LEG Dfs --------------------------
            for (l = 0; l < fund_ndfs; l++)
            {
                fund_dfs[l] = fund_dom_dff[l] * exp(-fund_dom_gam[l] * r1[i] - fund_dom_gam2[l]);
            }

            FundLeg_Payoff(fund_leg, (long)(evt_date), (long)(evt_date), fund_dfs, &fund_payoff);

            //-------------- RA LEG Dfs ----------------------------
            if (do_dom)
            {
                for (l = 0; l < pn_dfs[0]; l++)
                {
                    ra_dfs[0][l] = dom_dff[l] * exp(-dom_gam[l] * r1[i] - dom_gam2[l]);
                }
            }
            if (do_for)
            {
                for (l = 0; l < pn_dfs[1]; l++)
                {
                    ra_dfs[1][l] = for_dff[l] * exp(-for_gam[l] * r2[i][j] - for_gam2[l]);
                }
            }

            // DFs in order to calculate the floating coupons of the RA leg
            if (RA->cpn_type != 0)
            {
                if (RA->float_cpn_is_dom_for == 0)
                {
                    for (l = 0; l < 2 * pn_dfs[2]; l++)
                    {
                        ra_dfs[2][l] =
                            cpn_float_dff[l] * exp(-cpn_float_gam[l] * r1[i] - cpn_float_gam2[l]);
                    }
                }
                else
                {
                    for (l = 0; l < 2 * pn_dfs[2]; l++)
                    {
                        ra_dfs[2][l] = cpn_float_dff[l] *
                                       exp(-cpn_float_gam[l] * r2[i][j] - cpn_float_gam2[l]);
                    }
                }
            }

            RA_Payoff(RA, (long)(evt_date), (long)(evt_date), ra_dfs, &ra_payoff);

            //-------------- FEE Df ----------------------------
            fee_dom_df = fee_dom_dff * exp(-fee_dom_gam * r1[i] - fee_dom_gam2);
            if (fabs((ctsqtostr->call[call_idx]).fee) > 1.0e-08)
            {
                fee = fee_dom_df * (ctsqtostr->call[call_idx]).fee;
            }
            else
            {
                fee = 0.0;
            }

            if ((ctsqtostr->call[call_idx]).pay_rec == 0)
            {
                payoff = ra_payoff - fund_payoff - fee;
            }
            else
            {
                payoff = -ra_payoff + fund_payoff - fee;
            }
            //			payoff = ra_payoff + fund_payoff - fee;

            prod_val[i][j][0] = payoff;

            if (err)
            {
                goto FREE_RETURN;
            }
        }
    }

FREE_RETURN:

    if (ra_dfs[0])
        free(ra_dfs[0]);
    if (ra_dfs[1])
        free(ra_dfs[1]);
    if (RA->cpn_type != 0)
    {
        if (ra_dfs[2])
            free(ra_dfs[2]);
    }
    if (ra_dfs)
        free(ra_dfs);

    if (ra_pdates[0])
        free(ra_pdates[0]);
    if (ra_pdates[1])
        free(ra_pdates[1]);
    if (RA->cpn_type != 0)
    {
        if (ra_pdates[2])
            free(ra_pdates[2]);
    }
    if (ra_pdates)
        free(ra_pdates);

    if (dom_dff)
        free_dvector(dom_dff, 0, pn_dfs[0] - 1);
    if (dom_gam)
        free_dvector(dom_gam, 0, pn_dfs[0] - 1);
    if (dom_gam2)
        free_dvector(dom_gam2, 0, pn_dfs[0] - 1);

    if (for_dff)
        free_dvector(for_dff, 0, pn_dfs[1] - 1);
    if (for_gam)
        free_dvector(for_gam, 0, pn_dfs[1] - 1);
    if (for_gam2)
        free_dvector(for_gam2, 0, pn_dfs[1] - 1);

    if (cpn_float_dff)
        free_dvector(cpn_float_dff, 0, 2 * pn_dfs[2] - 1);
    if (cpn_float_gam)
        free_dvector(cpn_float_gam, 0, 2 * pn_dfs[2] - 1);
    if (cpn_float_gam2)
        free_dvector(cpn_float_gam2, 0, 2 * pn_dfs[2] - 1);

    if (fund_dom_dff)
        free_dvector(fund_dom_dff, 0, fund_ndfs - 1);
    if (fund_dom_gam)
        free_dvector(fund_dom_gam, 0, fund_ndfs - 1);
    if (fund_dom_gam2)
        free_dvector(fund_dom_gam2, 0, fund_ndfs - 1);

    if (pn_dfs)
        free(pn_dfs);

    if (fund_dfs)
        free(fund_dfs);
    if (fund_pdates)
        free(fund_pdates);

    return err;
}

/*	Main pricing function */

/*	Launch the pde */
Err ctsqto_launch_adi(
    CTSQTO_STR     cstqto,
    LGMQTO_ADI_ARG adi_arg,
    /*	Result */
    double* prem)
{
    double temp_val[2];
    Err    err;

    /* launch the corresponding ADI */

    err = doublelgm1fQuanto_adi_correl2(
        // Time data
        adi_arg->nstp,
        adi_arg->time,
        adi_arg->date,
        adi_arg->nstpx,
        adi_arg->domlambda,
        adi_arg->forlambda,
        adi_arg->sig_time,
        adi_arg->domsig,
        adi_arg->forsig,
        adi_arg->fxsig,
        adi_arg->nb_sig,
        adi_arg->quantorho,
        adi_arg->domforrho,
        adi_arg->void_prm,
        adi_arg->is_event,
        adi_arg->dom_ifr,
        adi_arg->for_ifr,
        adi_arg->dom_yc,
        adi_arg->for_yc,
        ctsqto_payoff_adi,
        1,
        (double*)temp_val);

    *prem = temp_val[0];

    return err;
}

/*	Free all structures */
Err ctsqto_free_all_struct(CTSQTO_STR ctsqto, LGMQTO_UND und, int call_feat, LGMQTO_ADI_ARG adi_arg)
{
    free_lgmQto_und(und);
    ctsqto_free_calls(ctsqto);

    if (ctsqto->fund_leg)
    {
        free_funding_leg(ctsqto->fund_leg);
        free(ctsqto->fund_leg);
    }

    if (ctsqto->ra_leg)
    {
        ra_free(ctsqto->ra_leg);
        free(ctsqto->ra_leg);
    }

    lgmQto_free_adi_arg(adi_arg);

    return NULL;
}

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

Err ctsqto_fill_check_all_struct(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund, 1: calibrate */
    /*		if calib */
    char*  dom_yc,         /*	yc */
    char*  dom_vc,         /*	vc (only if calib) */
    char*  dom_ref,        /*	ref rate (only if calib) */
    char*  dom_swap_freq,  /*	swap freq (only if calib) */
    char*  dom_swap_basis, /*	swap basis (only if calib) */
    double dom_lambda,     /*	lambda if unique */

    char*  for_yc,         /*	yc */
    char*  for_vc,         /*	vc */
    char*  for_ref,        /*	ref rate */
    char*  for_swap_freq,  /*	swap freq */
    char*  for_swap_basis, /*	swap basis */
    double for_lambda,     /*	lambda if unique */

    int forcalib, /*	0 : RA Und, 1 : Diag */

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
    char* fx3dund,
    /*	The structure */

    /*		funding */
    double  fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    /*		ra */
    double  ra_not,
    int     ra_cpn_type,
    int     n_fxvol_dates,
    long*   fxvol_dates,
    double* fxvol,
    double* qtocorrel,

    int typeVol,

    int      n_periods,
    long*    dates, /* n_periods + 1 */
    double*  cpns,
    char*    basis,
    int*     ra_nfixings,
    long**   ra_fixingdates,
    double** ra_fixings, /*	Past coupon fixing if relevant */

    // RA floating coupons
    int     ra_float_refrate_is_dom_for,
    char*   ra_float_refrate,
    long*   ra_float_startl,
    double* ra_float_past_fixings,
    double* ra_float_gearings,

    double* upper_barr,
    double* lower_barr,
    char*   recpay,
    char*   ra_refrate,
    double  c_spread,
    int     obs_freq,
    double  rho_df,

    // Params for floating coupon
    double correl_start,
    double correl_end,
    int    float_adj_strike,

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
    int    dom_force_atm, /*	force atm calib domestic und */
    int    for_force_atm, /*	force atm calib foreign und */
    double max_std_long,
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
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,

    double* corr_times,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    int     corr_n_times,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    //------Spot Lag-------------
    int fund_spot_lag,
    int ra_spot_lag,

    /*	Results */
    CTSQTO_STR ctsqto,
    LGMQTO_UND und,
    int*       call_feat, /*	0: No callable feature to be valued
                                  1: Callable feature to be valued through adi */
    LGMQTO_ADI_ARG adi_arg)
{
    Err err = NULL;

    int     i;
    double* event_time = NULL;

    CPD_DIAG_CALIB_PARAM calibparam = NULL;

    calibparam = (CPD_DIAG_CALIB_PARAM)malloc(sizeof(cpd_diag_calib_param));
    if (!calibparam)
    {
        err = "Memory allocation failed in ctsqto_fill_check_all_struct";
        goto FREE_RETURN;
    }

    cpd_calib_set_default_param(calibparam);

    /*	Initialisation */
    ctsqto->fund_leg = NULL;
    ctsqto->ra_leg   = NULL;
    ctsqto->call     = NULL;

    und->sigma_n    = 0;
    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->dom_sigma  = NULL;
    und->for_sigma  = NULL;
    und->dom_lambda = dom_lambda;
    und->for_lambda = for_lambda;
    //	cpd_init_calib_inst_data (&(und->dom_inst_data));
    //	cpd_init_calib_inst_data (&(und->for_inst_data));
    und->has_inst_data = 0;
    und->today         = today;

    adi_arg->time      = NULL;
    adi_arg->date      = NULL;
    adi_arg->domsig    = NULL;
    adi_arg->forsig    = NULL;
    adi_arg->fxsig     = NULL;
    adi_arg->domforrho = NULL;
    adi_arg->quantorho = NULL;
    adi_arg->sig_time  = NULL;
    adi_arg->void_prm  = NULL;
    adi_arg->is_event  = NULL;
    adi_arg->dom_ifr   = NULL;
    adi_arg->for_ifr   = NULL;

    /*	Funding leg */

    ctsqto->fund_leg = (FUNDING_LEG)malloc(sizeof(funding_leg));
    if (!ctsqto->fund_leg)
    {
        err = "Memory allocation error (1) in ctsqto_fill_check_all_struct";
        goto FREE_RETURN;
    }

    err = fill_funding_leg(
        und->today,
        eod_fix_flag,
        fund_spot_lag,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,
        ctsqto->fund_leg);
    if (err)
    {
        smessage("Error in fill funding leg");
        goto FREE_RETURN;
    }

    /*	Exotic leg */

    ctsqto->ra_leg = (RangeAccrualStruct*)malloc(sizeof(RangeAccrualStruct));
    if (!ctsqto->ra_leg)
    {
        err = "Memory allocation error (2) in ctsqto_fill_check_all_struct";
        goto FREE_RETURN;
    }

    err = ra_init_struct(
        today,
        dom_yc,
        for_yc,
        dom_vc,
        for_vc,
        dom_ref,
        for_ref,
        ra_not,
        ra_cpn_type,
        n_fxvol_dates,
        fxvol_dates,
        fxvol,
        qtocorrel,
        typeVol,
        n_periods,
        dates,
        cpns,
        basis,
        ra_nfixings,
        ra_fixingdates,
        ra_fixings,

        // RA float coupons
        ra_float_refrate_is_dom_for,
        ra_float_refrate,
        ra_float_startl,
        ra_float_past_fixings,
        ra_float_gearings,

        upper_barr,
        lower_barr,
        recpay,
        ra_refrate,
        c_spread,
        obs_freq,
        rho_df,

        // Params for the floating coupon
        correl_start,
        correl_end,
        float_adj_strike,

        eod_fix_flag,
        ctsqto->ra_leg);
    if (err)
    {
        smessage("Error in ra init");
        goto FREE_RETURN;
    }

    /*	Calls */

    if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag)
    {
        err = ctsqto_fill_calls(
            und->today, eod_ex_flag, ncall, pay_rec, ex_date, set_date, fee, ctsqto);
        if (err)
        {
            smessage("Error in fill calls");
            goto FREE_RETURN;
        }
    }
    else
    {
        ctsqto->num_calls = 0;
        ctsqto->call      = NULL;
    }

    /*	Underlying */
    if (ctsqto->num_calls > 0 && ctsqto->call)
    {
        if (use_calib)
        {
            err = ctsqto_calib_und(
                today,
                eod_ex_flag,
                dom_yc,
                dom_vc,
                dom_ref,
                dom_swap_freq,
                dom_swap_basis,
                for_yc,
                for_vc,
                for_ref,
                for_swap_freq,
                for_swap_basis,
                forcalib,
                dom_lambda,
                for_lambda,
                dom_force_atm,
                for_force_atm,
                max_std_long,
                max_std_short,
                fix_lambda,
                one_f_equi,
                skip_last,
                fx_mkt_vol_date,
                fx_mkt_vol,
                num_fx_mkt_vol,
                corr_times,
                correl_dom_for,
                correl_dom_fx,
                correl_for_fx,
                corr_n_times,
                calibparam,
                ctsqto,
                get_cash_vol,
                und);
            if (err)
            {
                smessage("Error in calibration");
                goto FREE_RETURN;
            }
        }
        else
        {
            err = fill_lgmQto_und(
                fx3dund,
                dom_vc,
                for_vc,
                dom_ref,
                dom_swap_freq,
                dom_swap_basis,
                for_ref,
                for_swap_freq,
                for_swap_basis,
                und);
            if (err)
            {
                goto FREE_RETURN;
            }
        }
    }

    smessage("Calib OK");
    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Tree */

    if (ctsqto->num_calls > 0 && ctsqto->call)
    {
        event_time = dvector(0, ctsqto->num_calls - 1);
        for (i = 0; i < ctsqto->num_calls; ++i)
        {
            event_time[i] = ctsqto->call[i].ex_time;
        }
        err = lgmQto_fill_adi_arg(
            und,
            ctsqto->num_calls,
            event_time,
            (void*)ctsqto,
            get_cash_vol,
            req_stp,
            req_stpx,
            adi_arg);

        smessage("Fill ADI Arg OK");
        if (err)
        {
            smessage("Error in fill adi arg");
            goto FREE_RETURN;
        }

        *call_feat = 1;
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (calibparam)
        free(calibparam);

    if (event_time)
        free_dvector(event_time, 0, ctsqto->num_calls - 1);

    if (err)
    {
        ctsqto_free_all_struct(ctsqto, und, *call_feat, adi_arg);
    }

    return err;
}

static double Ji0_func(double lam, double t1, double t2)
{
    return (exp(lam * t2) - exp(lam * t1)) / lam;
}

static Err LGMQTOPhi(
    int     nstept,
    double* time,
    double  domlam,
    double  forlam,
    double* sig_time,  //	domsig, forsig and fxsig must have
    double* domsig,    //	the same sig_time
    double* forsig,
    double* fxsig,
    int     nb_sig,
    double* quantorho,
    double* domforrho,
    double* fwd1,
    double* fwd2,
    double* std1,
    double* std2,
    double* corr,
    double* domphi,
    double* forphi,
    double* domforphi,
    double* sigmadom,
    double* sigmafor,
    double* quantocorrel,
    double* domforcorrel)
{
    Err    err = NULL;
    double t1, t2, ta, tb;
    int    i, j, nb_sig_minus1;
    double domI1, domI2;
    double forI1, forI2;
    double domforI1, domforI2;
    double domH, forH;
    double doma, fora, quantorhoa, domforrhoa;
    double domb, forb, quantorhob, domforrhob;
    double QAdj;

    nb_sig_minus1 = nb_sig - 1;

    // initialisation
    t1 = 0.0;
    j  = 0;

    fwd1[0] = 0;
    fwd2[0] = 0;

    doma       = domsig[0];
    domb       = 0;
    fora       = forsig[0];
    forb       = 0;
    quantorhoa = quantorho[0];
    quantorhob = 0;
    domforrhoa = domforrho[0];
    domforrhob = 0;

    sigmadom[0]     = domsig[0];
    sigmafor[0]     = forsig[0];
    quantocorrel[0] = quantorho[0];
    domforcorrel[0] = domforrho[0];

    domphi[0]    = 0.0;
    forphi[0]    = 0.0;
    domforphi[0] = 0.0;

    domI1    = 0;
    domI2    = 0;
    forI1    = 0;
    forI2    = 0;
    domforI1 = 0;
    domforI2 = 0;
    domH     = 0;
    forH     = 0;
    QAdj     = 0;

    for (i = 1; i < nstept; i++)
    {
        t2 = time[i];
        ta = t1;
        tb = sig_time[j];

        //		std1 = 0;

        while (tb < t2 && j < nb_sig_minus1)
        {
            doma = domsig[j];
            domb = 0;

            fora = forsig[j];
            forb = 0;

            quantorhoa = quantorho[j];
            quantorhob = 0;

            domforrhoa = domforrho[j];
            domforrhob = 0;

            domI1 += doma * doma * Ji0_func(domlam, ta, tb);

            domI2 += doma * doma * Ji0_func(2 * domlam, ta, tb);

            forI1 += fora * fora * Ji0_func(forlam, ta, tb);

            forI2 += fora * fora * Ji0_func(2 * forlam, ta, tb);

            domforI1 += domforrhoa * doma * fora * Ji0_func(forlam, ta, tb);

            domforI2 += domforrhoa * doma * fora * Ji0_func(domlam + forlam, ta, tb);

            QAdj += fxsig[j] * quantorhoa * fora * Ji0_func(forlam, ta, tb);

            j++;
            ta = tb;
            tb = sig_time[j];
        }

        sigmadom[i]     = domsig[j];
        sigmafor[i]     = forsig[j];
        quantocorrel[i] = quantorho[j];
        domforcorrel[i] = domforrho[j];

        doma = domsig[j];
        domb = 0;

        fora = forsig[j];
        forb = 0;

        quantorhoa = quantorho[j];
        quantorhob = 0;

        domforrhoa = domforrho[j];
        domforrhob = 0;

        domI1 += doma * doma * Ji0_func(domlam, ta, t2);

        domI2 += doma * doma * Ji0_func(2 * domlam, ta, t2);

        forI1 += fora * fora * Ji0_func(forlam, ta, t2);

        forI2 += fora * fora * Ji0_func(2 * forlam, ta, t2);

        domforI1 += domforrhoa * doma * fora * Ji0_func(forlam, ta, t2);

        domforI2 += domforrhoa * doma * fora * Ji0_func(domlam + forlam, ta, t2);

        QAdj += fxsig[j] * quantorhoa * fora * Ji0_func(forlam, ta, t2);

        domphi[i]    = exp(-2 * domlam * t2) * domI2;
        forphi[i]    = exp(-2 * forlam * t2) * forI2;
        domforphi[i] = exp(-(domlam + forlam) * t2) * domforI2;

        fwd1[i] = (exp(-domlam * t2) * domI1 - exp(-2 * domlam * t2) * domI2) / domlam;

        fwd2[i] =
            (exp(-forlam * t2) * forI1 - exp(-2 * forlam * t2) * forI2) / forlam -
            (exp(-forlam * t2) * domforI1 - exp(-(domlam + forlam) * t2) * domforI2) / domlam -
            exp(-forlam * t2) * QAdj;

        std1[i - 1] = sqrt(domphi[i]);
        std2[i - 1] = sqrt(forphi[i]);
        if (std1[i - 1] * std2[i - 1] > 0)
        {
            corr[i - 1] = domforphi[i] / (std1[i - 1] * std2[i - 1]);
        }
        else
        {
            corr[i - 1] = domforrho[0];
        }

        //		std3[i-1] = sqrt(
        //							(domsig[j]/forsig[j])*(domsig[j]/forsig[j]) *
        //forphi[i] / (1-domforrho[j]*domforrho[j]) 							+domforrho[j]*domforrho[j]*domphi[i] /
        //(1-domforrho[j]*domforrho[j]) 							-2*domforrho[j]*(domsig[j]/forsig[j])*domforphi[i] /
        //(1-domforrho[j]*domforrho[j])
        //							);

        t1 = t2;
    }

    return err;
}

Err ctsqto_calc_mdl_iv_fwd(
    CTSQTO_STR     ctsqto,
    LGMQTO_UND     und,
    LGMQTO_ADI_ARG adi_arg,

    int num_hermite,
    /*	Result */
    double* premium)
{
    int     i, j, k, index;
    double *time = NULL, *x = NULL, *w = NULL, *domphi = NULL, *forphi = NULL, *domforphi = NULL,
           *r1 = NULL, **r2 = NULL, ***pv = NULL, *std1vec = NULL, *std2vec = NULL, *corrvec = NULL,
           *fwd1 = NULL, *fwd2 = NULL;

    double* sigmadom     = NULL;
    double* sigmafor     = NULL;
    double* quantocorrel = NULL;
    double* domforcorrel = NULL;

    double std1, std2, corr;
    double sum, sum_part, dom_df, fee;

    double test, test_for_gam, test_for_gam2;

    CTSQTO_CALL call;

    Err err = NULL;

    /* memory allocation */

    time = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));

    x = dvector(1, num_hermite);
    w = dvector(1, num_hermite);

    domphi       = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    forphi       = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    domforphi    = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    sigmadom     = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    sigmafor     = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    quantocorrel = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    domforcorrel = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));

    fwd1 = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));
    fwd2 = (double*)calloc(ctsqto->num_calls + 1, sizeof(double));

    std1vec = (double*)calloc(ctsqto->num_calls, sizeof(double));
    std2vec = (double*)calloc(ctsqto->num_calls, sizeof(double));
    corrvec = (double*)calloc(ctsqto->num_calls, sizeof(double));

    r1 = dvector(0, num_hermite - 1);
    r2 = dmatrix(0, num_hermite - 1, 0, num_hermite - 1);

    pv = f3tensor(0, num_hermite - 1, 0, num_hermite - 1, 0, 0);

    if (!time || !x || !w || !domphi || !forphi || !domforphi || !sigmadom || !sigmafor ||
        !quantocorrel || !domforcorrel || !fwd1 || !fwd2 || !std1vec || !std2vec || !corrvec ||
        !r1 || !r2 || !pv)
    {
        err = "Memory allocation failure in ctsqto_calc_mdl_iv_fwd";
        goto FREE_RETURN;
    }

    /* Hermite calculation */
    err = HermiteStandard(x, w, num_hermite);

    if (err)
    {
        goto FREE_RETURN;
    }

    time[0] = 0.0;
    for (i = 0; i < ctsqto->num_calls; i++)
    {
        call        = ctsqto->call + i;
        time[i + 1] = call->ex_time;
    }

    /* Phi Calculation */
    err = LGMQTOPhi(
        ctsqto->num_calls + 1,
        time,
        und->dom_lambda,
        und->for_lambda,
        und->sigma_time,
        und->dom_sigma,
        und->for_sigma,
        und->fx_sigma,
        und->sigma_n,
        und->quanto_rho,
        und->domfor_rho,
        fwd1,
        fwd2,
        std1vec,
        std2vec,
        corrvec,
        domphi,
        forphi,
        domforphi,
        sigmadom,
        sigmafor,
        quantocorrel,
        domforcorrel);

    if (err)
    {
        goto FREE_RETURN;
    }

    index = 0;

    for (i = 0; i < ctsqto->num_calls; i++)
    {
        /* first remove the fee */
        call      = ctsqto->call + i;
        fee       = call->fee;
        call->fee = 0.0;

        std1 = std1vec[i];
        std2 = std2vec[i];
        corr = corrvec[i];

        // Corresponding dom_r and for_r
        for (j = 0; j < num_hermite; j++)
        {
            r1[j] = std1 * x[j + 1];

            for (k = 0; k < num_hermite; k++)
            {
                r2[j][k] =
                    fwd2[i + 1] + std2 * (corr * x[j + 1] + sqrt(1 - corr * corr) * x[k + 1]);
            }
        }

        /* launch evaluation */

        while (!((adi_arg->void_prm)[index]) && index < 10000)
        {
            index++;
        }

        err = ctsqto_payoff_adi2(
            call->ex_date,
            call->ex_time,
            (adi_arg->void_prm)[index],
            und->dom_yc,
            und->for_yc,
            und->dom_lambda,
            und->for_lambda,
            domphi[i + 1],
            forphi[i + 1],
            0,
            num_hermite - 1,
            0,
            num_hermite - 1,
            r1,
            r2,
            1,
            pv);

        index++;

        /* put the fee back */
        call->fee = fee;

        /* convolution */

        sum      = 0.0;
        sum_part = 0.0;

        test = 0.0;

        for (j = 0; j < num_hermite; j++)
        {
            sum_part = 0.0;
            for (k = 0; k < num_hermite; k++)
            {
                sum_part += w[k + 1] * pv[j][k][0];
                pv[j][k][0] = 0;

                if (i < ctsqto->num_calls - 1)
                {
                    test_for_gam = (1.0 - exp(-und->for_lambda * (ctsqto->call[i + 1].ex_time -
                                                                  ctsqto->call[i].ex_time))) /
                                   und->for_lambda;
                    test_for_gam2 = 0.5 * test_for_gam * test_for_gam * forphi[i + 1];
                    test += w[j + 1] *
                            swp_f_df(call->ex_date, ctsqto->call[i + 1].ex_date, und->for_yc) *
                            exp(-test_for_gam * r2[j][k] - test_for_gam2);
                }
            }
            sum = sum + w[j + 1] * sum_part;
        }

        dom_df = swp_f_df(und->today, call->ex_date, und->dom_yc);

        test *= dom_df;

        sum *= dom_df;

        premium[i] = sum;
    }

FREE_RETURN:

    if (time)
        free(time);

    if (x)
        free_dvector(x, 1, num_hermite);
    if (w)
        free_dvector(w, 1, num_hermite);

    if (domphi)
        free(domphi);
    if (forphi)
        free(forphi);
    if (domforphi)
        free(domforphi);

    if (sigmadom)
        free(sigmadom);
    if (sigmafor)
        free(sigmafor);
    if (quantocorrel)
        free(quantocorrel);
    if (domforcorrel)
        free(domforcorrel);
    if (fwd1)
        free(fwd1);
    if (fwd2)
        free(fwd2);
    if (std1vec)
        free(std1vec);
    if (std2vec)
        free(std2vec);
    if (corrvec)
        free(corrvec);

    if (r1)
        free_dvector(r1, 0, num_hermite - 1);
    if (r2)
        free_dmatrix(r2, 0, num_hermite - 1, 0, num_hermite - 1);

    if (pv)
        free_f3tensor(pv, 0, num_hermite - 1, 0, num_hermite - 1, 0, 0);

    return err;
}
