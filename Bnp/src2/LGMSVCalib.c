#include "lgmsvcalib.h"

#include "cpdcalib.h"
#include "lgmsvgrfn.h"
#include "lgmsvpde.h"
#include "lgmsvutil.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

#define MAX_CPN 600

static Err get_end_date(
    long ex_date, long struct_end_date, char* tenor, int theo_act, long* end_date)
{
    long start_date;
    Err  err;

    start_date = add_unit(ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    strupper(tenor);
    strip_white_space(tenor);
    if (*tenor == 'D')
    {
        *end_date = struct_end_date;
    }
    else
    {
        err = add_tenor(start_date, tenor, NO_BUSDAY_CONVENTION, end_date);
        if (err)
        {
            return err;
        }
    }

    if (theo_act)
    {
        *end_date = bus_date_method(*end_date, MODIFIED_SUCCEEDING);
    }

    return NULL;
}

/*	Setup G function from lambda */
static void static_lgmsetupG(
    double lambda,
    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G[],    /*	Output: G at cash-flow dates
                                                       G(T) = (1.0 - exp (- lambda * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G[])     /*	Output: G at exercise dates */
{
    int i;

    for (i = 0; i < ncpn; i++)
    {
        cpn_G[i] = (1.0 - exp(-lambda * cpn_time[i])) / lambda;
    }

    for (i = 0; i < nex; i++)
    {
        ex_G[i] = (1.0 - exp(-lambda * ex_time[i])) / lambda;
    }
}

/*	Calibrate and price cap given lambda */
Err lgmSVprcapgivenlambda_num(
    int    ncpn,         /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int    nex,          /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int    ex_cpn[],     /*	Index of the first cash-flow to be exercised */
    int    ex_endcpn[],  /*	Index of the last cash-flow to be exercised */
    int    ex_sncpn[],   /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    double ex_sstrike[], /*	Strikes for cap */
    double vol[],        /*	Output: volatility structure */
    double lambda,       /*	Lambda */

    /*	Market */
    long  today,
    char* yc,

    /*	Lambda Alpha and Rho */
    double alpha,
    double lambdaEps,
    double rho,

    int skip_last,     /*	If 1, the last option is disregarded
                                               and the forward volatility is flat from option
                                               n-1 */
    int     price_cap, /*	0: just calibrate */
    double* ex_sprice, /*	Cap price as output */

    /*	Space Discretisation	*/
    int nbT,
    int nbPhi,
    int nbX,
    int nbEps,

    /*	Newton Parameters		*/
    int    nbIterMax,
    double precision,

    /*	Model Params */
    LGMSVParam Params)
{
    int i, j, nbIter;

    static double ex_zeta[MAX_CPN], cpn_G[MAX_CPN], ex_G[MAX_CPN];

    double *time = NULL, *date = NULL, *ifr = NULL;

    int* is_event = NULL;

    static lgmsv_swap_param param;

    void** func_parm = NULL;

    double dt, df;

    double Var, newVar, vol1, vol2, price1, price2;
    double error, save_vol, vega;

    double level, fwd, alphaeq, bsvol, bsvol1, sigbeta, sigbeta1;

    time_t t1, t2;

    Err err = NULL;

    t1 = clock();
    smessage("Starting Calibration of LGM SV");

    if (nex == 1)
    {
        skip_last = 0;
    }

    /* First calibrate a regular 1F model */

    /*	Setup G function */
    static_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

    err = lgmcalibzeta1F_2(
        ncpn,
        cpn_time,
        cpn_df,
        cpn_cvg,
        cpn_G,
        nex,
        ex_time,
        ex_cpn,
        ex_endcpn,
        ex_G,
        ex_lstrike,
        ex_lprice,
        ex_zeta,
        lambda,
        skip_last);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Transform zeta into sigma */
    vol[0] = sqrt(ex_zeta[0] * 2 * lambda / (exp(2 * lambda * ex_time[0]) - 1.0));

    for (i = 1; i < nex; i++)
    {
        if (ex_zeta[i] > ex_zeta[i - 1])
        {
            vol[i] = sqrt(
                (ex_zeta[i] - ex_zeta[i - 1]) * 2 * lambda /
                (exp(2 * lambda * ex_time[i]) - exp(2 * lambda * ex_time[i - 1])));
        }
        else
        {
            smessage(
                "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                ex_time[i]);
            for (j = i; j < nex; j++)
            {
                vol[j] = vol[i - 1];
            }
            i = nex;
        }
    }

    /*	Setup the time vector */
    date = dvector(0, nbT - 1);
    time = dvector(0, nbT - 1);
    ifr  = dvector(0, nbT - 1);

    is_event  = (int*)calloc(nbT, sizeof(int));
    func_parm = (void**)malloc(nbT * sizeof(void*));

    if (!date || !time || !ifr || !is_event || !func_parm)
    {
        err = "Memory Allocation error (1) in lgmSVprcapgivenlambda_num";
        goto FREE_RETURN;
    }

    is_event[nbT - 1] = 1;
    memcpy(param.cpn_cvg, cpn_cvg, ncpn * sizeof(double));

    func_parm[nbT - 1] = (void*)(&param);

    Var = 0.0;

    for (i = 0; i < nex - skip_last; i++)
    {
        nbIter = 0;

        /*	Fill time information */
        dt = ex_time[i] / (nbT - 1);

        time[nbT - 1] = ex_time[i];
        date[nbT - 1] = today + time[nbT - 1] * DAYS_IN_YEAR + 1.0E-08;
        ifr[nbT - 1]  = swp_f_zr(date[nbT - 1], date[nbT - 1] + 1, yc);

        for (j = nbT - 2; j >= 0; j--)
        {
            time[j] = time[j + 1] - dt;
            date[j] = today + time[j] * DAYS_IN_YEAR + 1.0E-08;

            if (date[j + 1] > date[j])
            {
                ifr[j] = swp_f_zr(date[j], date[j + 1], yc);
            }
            else
            {
                ifr[j] = swp_f_zr(date[j], date[j] + 1, yc);
            }
        }

        time[0] = 0.0;
        date[0] = today;

        /*	Fill the payoff parameters */
        param.start_index = ex_cpn[i];
        param.strike      = ex_lstrike[i];
        param.ncpn        = ex_endcpn[i] + 1;

        df = swp_f_df(today, today + ex_time[i] * DAYS_IN_YEAR + 1.0E-08, yc);

        for (j = param.start_index; j < param.ncpn; j++)
        {
            param.cpn_df[j]   = cpn_df[j] / df;
            param.cpn_time[j] = cpn_time[j] - ex_time[i];
            param.gam[j]      = (1.0 - exp(-lambda * param.cpn_time[j])) / lambda;
            param.gam2[j]     = 0.5 * param.gam[j] * param.gam[j];
        }

        level = 0.0;

        for (j = param.start_index + 1; j < param.ncpn; j++)
        {
            level += cpn_df[j] * cpn_cvg[j];
        }

        fwd = (cpn_df[param.start_index] - cpn_df[param.ncpn - 1]) / level;

        err = srt_f_optimpvol(
            ex_lprice[i], fwd, ex_lstrike[i], ex_time[i], level, SRT_PUT, SRT_LOGNORMAL, &bsvol);

        if (fabs(lambdaEps) > 1.0E-10)
        {
            alphaeq =
                alpha *
                sqrt((1.0 - exp(-2.0 * lambdaEps * ex_time[i])) / (2.0 * lambdaEps * ex_time[i]));
        }
        else
        {
            alphaeq = alpha;
        }

        err = vol_conv(
            bsvol,
            SABR_ATM_LOG,
            &sigbeta,
            SABR_ATM_BETA,
            fwd,
            ex_lstrike[i],
            ex_time[i],
            alphaeq,
            0.0,
            rho);

        /*	First guess and first evaluation */

        if (i > 0)
        {
            newVar = vol[i - 1] * vol[i - 1] * ex_time[i - 1] + vol[i] * vol[i] * ex_time[i] -
                     save_vol * save_vol * ex_time[i - 1];
            save_vol = vol[i];

            if (newVar > 0.0)
            {
                vol[i] = sqrt(newVar / ex_time[i]);
            }
            else
            {
                vol[i] = vol[i - 1];
            }

            newVar = Var + vol[i] * vol[i] * (ex_time[i] - ex_time[i - 1]);
        }
        else
        {
            save_vol = vol[0];
            newVar   = vol[0] * vol[0] * ex_time[0];
        }

        vol1 = vol[i];

        err = lgmSV_adi(
            nbT,
            time,
            date,
            nbPhi,
            nbX,
            nbEps,
            lambda,
            1.0,
            ex_time,
            vol,
            i + 1,
            lambdaEps,
            alpha,
            rho,
            Params,
            func_parm,
            is_event,
            ifr,
            yc,
            payoff_lgmsv_rec_swaption,
            1,
            &price1);

        if (err)
        {
            goto FREE_RETURN;
        }

        error = fabs(price1 - ex_lprice[i]);

        if (error > precision)
        {
            /* Need to adjust */

            /* Second guess */

            /*
            if (i > 0)
            {
                    vol2 = (ex_lprice[i] * ex_lprice[i] / price1 / price1 * newVar - Var) /
            (ex_time[i] - ex_time[i-1]);
            }
            else
            {
                    vol2 = ex_lprice[i] * ex_lprice[i] / price1 / price1 * newVar / ex_time[i];
            }

            if (vol2 > 0.0)
            {
                    vol2 = sqrt(vol2);
            }
            else
            {
                    vol2 = vol1 * 0.75;
            }
            */

            err = srt_f_optimpvol(
                price1, fwd, ex_lstrike[i], ex_time[i], level, SRT_PUT, SRT_LOGNORMAL, &bsvol1);

            err = vol_conv(
                bsvol1,
                SABR_ATM_LOG,
                &sigbeta1,
                SABR_ATM_BETA,
                fwd,
                ex_lstrike[i],
                ex_time[i],
                alphaeq,
                0.0,
                rho);

            newVar *= sigbeta / sigbeta1 * sigbeta / sigbeta1;

            if (i > 0)
            {
                vol2 = (newVar - Var) / (ex_time[i] - ex_time[i - 1]);
            }
            else
            {
                vol2 = (newVar - Var) / ex_time[i];
            }

            if (vol2 > 0.0)
            {
                vol2 = sqrt(vol2);
            }
            else
            {
                vol2 = vol1 * 0.75;
            }

            vol[i] = vol2;

            /*	Simple Newtonalgorithm */

            nbIter++;

            while (error > precision && nbIter < nbIterMax)
            {
                err = lgmSV_adi(
                    nbT,
                    time,
                    date,
                    nbPhi,
                    nbX,
                    nbEps,
                    lambda,
                    1.0,
                    ex_time,
                    vol,
                    i + 1,
                    lambdaEps,
                    alpha,
                    rho,
                    /* Hard coded to be modified */
                    Params,
                    func_parm,
                    is_event,
                    ifr,
                    yc,
                    payoff_lgmsv_rec_swaption,
                    1,
                    &price2);

                if (err)
                {
                    goto FREE_RETURN;
                }

                error = fabs(price2 - ex_lprice[i]);

                vega   = (price2 - price1) / (vol2 - vol1);
                vol1   = vol2;
                price1 = price2;

                if (vega < 0.0)
                {
                    /* strange, isn't it ? */
                    vol2 = vol1 * 0.75;
                }
                else
                {
                    vol2 = vol1 + (ex_lprice[i] - price1) / vega;

                    if (vol2 < 0.0)
                    {
                        vol2 = vol1 * 0.75;
                    }
                }

                if (i > 0 && (vol2 < 0.5 * vol[i - 1]))
                {
                    if (fabs(vol1 - 0.5 * vol[i - 1]) < 1.0E-08)
                    {
                        /* Failed */
                        vol2   = vol1;
                        nbIter = nbIterMax + 1;
                    }
                    else
                    {
                        vol2 = 0.5 * vol[i - 1];
                    }

                    error = 1.0E9;
                }

                vol[i] = vol2;

                nbIter++;
            }
        }

        if (error < precision)
        {
            smessage("Option %d: success at iteration %d", i + 1, nbIter);
        }
        else
        {
            smessage("Option %d: may have failed", i + 1, nbIter - 1);
        }

        if (i > 0)
        {
            Var += vol[i] * vol[i] * (ex_time[i] - ex_time[i - 1]);
        }
        else
        {
            Var = vol[0] * vol[0] * ex_time[0];
        }
    }

    if (skip_last)
    {
        vol[i] = vol[i - 1];
    }

    t2 = clock();
    smessage("Calibration time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

    if (date)
        free_dvector(date, 0, nbT - 1);
    if (time)
        free_dvector(time, 0, nbT - 1);
    if (ifr)
        free_dvector(ifr, 0, nbT - 1);
    if (is_event)
        free(is_event);
    if (func_parm)
        free(func_parm);

    return err;
}

/*	Calibrate lgm: main function */
Err lgmsv_calib_diagonal(
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    char* ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    double vol_shift,
    int    shift_type,    /*	0:	Additive
                                          1:	Multiplicative */
                          /*	If ex_date is NULL,
                          exercise dates will be generated 2bd before start */
    int num_ex_dates,     /*	Exercise dates,
                                                          all supposed to be on or after today */
    long* ex_date,        /*	Supposed to be sorted
                                                          NULL = 2bd before each coupon */
    char** end_tenor,     /*	Tenors of the underlying instruments
                                                          or "DIAG" */
    long    end_date,     /*	End date for diagonal */
    double* long_strike,  /*	Diagonal swaption strikes
                                                                  NULL = ATM */
    double* short_strike, /*	Short swaption strikes
                                                                  NULL = ATM */
    int strike_type,      /*	0: ATM
                                                  1: CASH
                                                  2: SWAP
                                                  3: STD */
    double max_std_long,
    double max_std_short,
    char*  swaption_freq, /*	Frequency and basis of underlying swaptions */
    char*  swaption_basis,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                               if set to 1, then 2F lambda will calibrate
                                               to the cap priced within calibrated 1F
                                               with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                               and the forward volatility is flat from option
                                               n-1 */
    double* lambda,    /*	Lambda: may be changed in the process */

    /*	Alpha, Gamma, Rho (2F only) */
    double   alpha,
    double   lambdaEps,
    double   rho,
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,

    /*	Space Discretisation	*/
    int nbT,
    int nbPhi,
    int nbX,
    int nbEps,

    /*	Newton Parameters		*/
    int    nbIterMax,
    double precision,

    /*	Model Params */
    LGMSVParam ModelParams,

    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data) /*	NULL = don't save calibration instrument data */
{
    int            i, j, k, l, nex, ncpn;
    SrtCompounding ifreq;
    SrtBasisCode   ibasis;
    double         ex_time[MAX_CPN], ex_lfwd[MAX_CPN], ex_llvl[MAX_CPN], ex_lstrike[MAX_CPN],
        ex_lvol[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
        ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sprice[MAX_CPN], ex_sweight[MAX_CPN],
        ex_zeta[MAX_CPN], ex_G[MAX_CPN];
    int         ex_cpn[MAX_CPN], ex_sncpn[MAX_CPN], ex_endcpn[MAX_CPN];
    long        cpn_date[MAX_CPN];
    double      cpn_time[MAX_CPN], cpn_cvg[MAX_CPN], cpn_df[MAX_CPN], cpn_G[MAX_CPN];
    long        tmplng1[MAX_CPN], tmplng2[MAX_CPN];
    long        theo_date, act_date, temp_date, temp_date2;
    long *      theo_end_dates, *act_end_dates;
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      cap_price;
    double      std;
    SrtCurvePtr yc_ptr;
    Err         err = NULL;

    theo_end_dates = &(tmplng1[0]);
    act_end_dates  = &(tmplng2[0]);

    *sig_time = NULL;
    *sig      = NULL;

    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    /*	1.)	Setup the bond schedule and its coupons */

    /*	Coupons */

    err = interp_compounding(swaption_freq, &ifreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = interp_basis(swaption_basis, &ibasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Find the end date as the longest total maturity */
    theo_date = end_date;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    for (i = 0; i < num_ex_dates; i++)
    {
        err = get_end_date(ex_date[i], end_date, end_tenor[i], 0, &(theo_end_dates[i]));
        if (err)
        {
            goto FREE_RETURN;
        }
        act_end_dates[i] = bus_date_method(theo_end_dates[i], MODIFIED_SUCCEEDING);
    }
    for (i = 0; i < num_ex_dates; i++)
    {
        if (theo_end_dates[i] > theo_date || act_end_dates[i] > act_date)
        {
            theo_date = theo_end_dates[i];
            act_date  = act_end_dates[i];
        }
    }
    ncpn       = 1;
    temp_date  = theo_date;
    temp_date2 = act_date;

    while (act_date > today)
    {
        theo_date = add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
        act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
        ncpn++;
    }
    ncpn--;

    if (ncpn < 2)
    {
        err = "Not enough coupons in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    theo_date = temp_date;
    act_date  = temp_date2;
    i         = ncpn - 1;

    while (i >= 0)
    {
        cpn_time[i] = (act_date - today) * YEARS_IN_DAY;
        cpn_date[i] = act_date;

        theo_date = add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

        temp_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
        cpn_cvg[i] = coverage(temp_date, act_date, ibasis);
        cpn_df[i]  = swp_f_df(today, act_date, yc_name);
        act_date   = temp_date;

        i--;
    }
    cpn_cvg[0] = 0.0;

    /*	Exercise */

    /*	Remove past dates */
    while (ex_date[0] <= today)
    {
        ex_date++;
        end_tenor++;
        theo_end_dates++;
        act_end_dates++;
        if (long_strike)
            long_strike++;
        if (short_strike)
            short_strike++;
        num_ex_dates--;
        if (num_ex_dates == 0)
        {
            "All exercise dates are past in cpd_calib_diagonal";
            goto FREE_RETURN;
        }
    }

    /*	Remove redundant dates */
    j = ncpn - 1;
    l = ncpn + 1;
    for (i = num_ex_dates - 1; i >= 0; i--)
    {
        while (j > 0 && cpn_date[j] > ex_date[i])
        {
            j--;
        }
        if (cpn_date[j] < ex_date[i])
        {
            j++;
        }

        if (j >= ncpn - 1 || j == l)
        {
            for (k = i - 1; k >= 0; k--)
            {
                ex_date[k + 1] = ex_date[k];
                strcpy(end_tenor[k + 1], end_tenor[k]);
                theo_end_dates[k + 1] = theo_end_dates[k];
                act_end_dates[k + 1]  = act_end_dates[k];
                if (long_strike)
                {
                    long_strike[k + 1] = long_strike[k];
                }
                if (short_strike)
                {
                    short_strike[k + 1] = short_strike[k];
                }
            }

            ex_date++;
            end_tenor++;
            theo_end_dates++;
            act_end_dates++;
            if (long_strike)
                long_strike++;
            if (short_strike)
                short_strike++;
            num_ex_dates--;
            if (num_ex_dates < 1)
            {
                err = "All exercise dates are past in cpd_calib_diagonal";
                goto FREE_RETURN;
            }
        }
        else
        {
            l = j;
        }
    }

    if (num_ex_dates < 1)
    {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    for (i = 0; i < num_ex_dates; i++)
    {
        ex_sweight[i] = 1.0;
    }

    nex = num_ex_dates;
    j   = 0;
    for (i = 0; i < nex; i++)
    {
        while (cpn_date[j] < ex_date[i])
        {
            j++;
        }

        ex_cpn[i]  = j;
        ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;

        k = j;
        while (cpn_date[k] < act_end_dates[i])
        {
            k++;
        }
        if (k > 0 && cpn_date[k] - act_end_dates[i] > act_end_dates[i] - cpn_date[k - 1])
        {
            k--;
        }

        if (k <= j)
        {
            k = j + 1;
        }
        ex_endcpn[i] = k;

        if (j >= ncpn || k >= ncpn)
        {
            err = "Coupon date bug in cpd_calib_diagonal";
            goto FREE_RETURN;
        }
    }

    /*	Underlyings */

    /*	Long */

    for (i = 0; i < nex; i++)
    {
        j = ex_cpn[i];
        l = ex_endcpn[i];

        lvl = 0.0;
        for (k = j + 1; k <= l; k++)
        {
            lvl += cpn_cvg[k] * cpn_df[k];
        }

        dfi = swp_f_df(today, cpn_date[j], yc_name);
        dff = swp_f_df(today, cpn_date[l], yc_name);

        ex_llvl[i] = lvl;
        ex_lfwd[i] = (dfi - dff) / lvl;

        /*	ATM std */
        err = get_cash_vol(
            vol_curve_name,
            add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
            theo_end_dates[i],
            ex_lfwd[i],
            0,
            ref_rate_name,
            &std,
            &power);

        if (err)
        {
            goto FREE_RETURN;
        }
        std += (shift_type == 1 ? std * vol_shift : vol_shift);
        if (power > 0.5)
        {
            power =
                srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, ex_time[i], 1.0, SRT_CALL, PREMIUM);
            err = srt_f_optimpvol(
                power, ex_lfwd[i], ex_lfwd[i], ex_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
        }
        std *= sqrt(ex_time[i]);

        /*	Strike */
        if ((!long_strike) || (!strike_type))
        {
            ex_lstrike[i] = ex_lfwd[i];
        }
        else if (strike_type == 1)
        {
            ex_lstrike[i] = long_strike[i];
        }
        else if (strike_type == 2)
        {
            if (err = swp_f_ForwardRate(
                    cpn_date[j],
                    theo_end_dates[i],
                    swaption_freq,
                    swaption_basis,
                    yc_name,
                    ref_rate_name,
                    &swp_rte))
            {
                goto FREE_RETURN;
            }

            spr = swp_rte - ex_lfwd[i];

            ex_lstrike[i] = long_strike[i] - spr;
        }
        else if (strike_type == 3)
        {
            ex_lstrike[i] = ex_lfwd[i] + long_strike[i] * std;
        }

        /*	Apply max std */
        if (ex_lstrike[i] > ex_lfwd[i] + max_std_long * std)
        {
            ex_lstrike[i] = ex_lfwd[i] + max_std_long * std;
        }
        else if (ex_lstrike[i] < ex_lfwd[i] - max_std_long * std)
        {
            ex_lstrike[i] = ex_lfwd[i] - max_std_long * std;
        }

        /*	Make sure strikes are positive (actually more than 1bp)
                        otherwise use ATM	*/
        if (ex_lstrike[i] < 1.0e-04)
        {
            ex_lstrike[i] = ex_lfwd[i];
        }

        err = get_cash_vol(
            vol_curve_name,
            add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
            theo_end_dates[i],
            ex_lstrike[i],
            0,
            ref_rate_name,
            &(ex_lvol[i]),
            &power);
        if (err)
        {
            goto FREE_RETURN;
        }
        ex_lvol[i] += (shift_type == 1 ? ex_lvol[i] * vol_shift : vol_shift);

        if (power > 0.5)
        {
            ex_lprice[i] = srt_f_optblksch(
                ex_lfwd[i], ex_lstrike[i], ex_lvol[i], ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
        }
        else
        {
            ex_lprice[i] = srt_f_optblknrm(
                ex_lfwd[i], ex_lstrike[i], ex_lvol[i], ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
        }
    }

    /*	Short */

    if (!fix_lambda)
    {
        cap_price = 0.0;
        for (i = 0; i < nex; i++)
        {
            if (i < nex - 1)
            {
                ex_sncpn[i] = ex_cpn[i + 1] - ex_cpn[i] + 1;
            }
            else
            {
                ex_sncpn[i] = ncpn - ex_cpn[i];
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            lvl = 0.0;
            for (k = ex_cpn[i] + 1; k < ex_cpn[i] + ex_sncpn[i]; k++)
            {
                lvl += cpn_cvg[k] * cpn_df[k];
            }
            dfi = swp_f_df(today, cpn_date[ex_cpn[i]], yc_name);
            dff = swp_f_df(today, cpn_date[ex_cpn[i] + ex_sncpn[i] - 1], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            /*	ATM std */
            err = get_cash_vol(
                vol_curve_name,
                add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                cpn_date[ex_cpn[i] + ex_sncpn[i] - 1],
                ex_sfwd[i],
                0,
                ref_rate_name,
                &std,
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            std += (shift_type == 1 ? std * vol_shift : vol_shift);
            if (power > 0.5)
            {
                power = srt_f_optblksch(
                    ex_sfwd[i], ex_sfwd[i], std, ex_time[i], 1.0, SRT_CALL, PREMIUM);
                err = srt_f_optimpvol(
                    power, ex_sfwd[i], ex_sfwd[i], ex_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
            }
            std *= sqrt(ex_time[i]);

            /*	Strike */
            if ((!short_strike) || (!strike_type))
            {
                ex_sstrike[i] = ex_sfwd[i];
            }
            else if (strike_type == 1)
            {
                ex_sstrike[i] = short_strike[i];
            }
            else if (strike_type == 2)
            {
                if (err = swp_f_ForwardRate(
                        cpn_date[ex_cpn[i]],
                        cpn_date[ex_cpn[i] + ex_sncpn[i] - 1],
                        swaption_freq,
                        swaption_basis,
                        yc_name,
                        ref_rate_name,
                        &swp_rte))
                {
                    goto FREE_RETURN;
                }

                spr = swp_rte - ex_sfwd[i];

                ex_sstrike[i] = short_strike[i] - spr;
            }
            else if (strike_type == 3)
            {
                ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
            }

            /*	Apply max std */
            if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
            }
            else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
            }

            /*	Make sure strikes are positive (actually more than 1bp)
                            otherwise use ATM	*/
            if (ex_sstrike[i] < 1.0e-04)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            err = get_cash_vol(
                vol_curve_name,
                add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                cpn_date[ex_cpn[i] + ex_sncpn[i] - 1],
                ex_sstrike[i],
                0,
                ref_rate_name,
                &(ex_svol[i]),
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

            if (power > 0.5)
            {
                ex_sprice[i] = srt_f_optblksch(
                    ex_sfwd[i],
                    ex_sstrike[i],
                    ex_svol[i],
                    ex_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }
            else
            {
                ex_sprice[i] = srt_f_optblknrm(
                    ex_sfwd[i],
                    ex_sstrike[i],
                    ex_svol[i],
                    ex_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }

            cap_price += ex_sprice[i];
        }
    }

    /*	The 1F equivalent case */

    if (fix_lambda && one_f_equi)
    {
        cap_price = 0.0;
        for (i = 0; i < nex; i++)
        {
            if (i < nex - 1)
            {
                ex_sncpn[i] = ex_cpn[i + 1] - ex_cpn[i] + 1;
            }
            else
            {
                ex_sncpn[i] = ncpn - ex_cpn[i];
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            lvl = 0.0;
            for (k = ex_cpn[i] + 1; k < ex_cpn[i] + ex_sncpn[i]; k++)
            {
                lvl += cpn_cvg[k] * cpn_df[k];
            }
            dfi = swp_f_df(today, cpn_date[ex_cpn[i]], yc_name);
            dff = swp_f_df(today, cpn_date[ex_cpn[i] + ex_sncpn[i] - 1], yc_name);

            ex_slvl[i]    = lvl;
            ex_sfwd[i]    = (dfi - dff) / lvl;
            ex_sstrike[i] = ex_sfwd[i];
        }

        err = lgmcalibzetalambda(
            ncpn,
            cpn_time,
            cpn_df,
            cpn_cvg,
            nex,
            ex_time,
            ex_cpn,
            ex_sncpn,
            ex_lstrike,
            ex_lprice,
            ex_sstrike,
            ex_sweight,
            0.0,
            ex_zeta,
            1,
            lambda,
            1,
            0.0,
            0.0,
            0.0,
            skip_last,
            CALPRES,
            CALPRES,
            0.25,
            0.25 / 4.0,
            0);

        if (err)
        {
            goto FREE_RETURN;
        }

        static_lgmsetupG(*lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

        cap_price = lgmcapval1F(
            ncpn,
            cpn_df,
            cpn_cvg,
            cpn_G,
            nex,
            ex_cpn,
            ex_sncpn,
            ex_sweight,
            ex_zeta,
            ex_G,
            ex_sstrike);

        fix_lambda = 0;
    }

    /*	2.)	Calibrate lambda and zeta */

    err = lgmSVprcapgivenlambda_num(
        ncpn,
        cpn_time,
        cpn_df,
        cpn_cvg,
        nex,
        ex_time,
        ex_cpn,
        ex_endcpn,
        ex_sncpn,
        ex_lstrike,
        ex_lprice,
        ex_sstrike,
        ex_zeta,
        *lambda,
        today,
        yc_name,
        2.0 * alpha,
        2.0 * lambdaEps,
        rho,
        skip_last,
        0,
        NULL,
        nbT,
        nbPhi,
        nbX,
        nbEps,
        nbIterMax,
        precision,
        ModelParams);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	3.)	Transform into sigma */

    *num_sig  = nex;
    *sig_time = (double*)calloc(nex, sizeof(double));
    *sig      = (double*)calloc(nex, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in lgmsv_calib_diagonal";
        goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++)
    {
        (*sig_time)[i] = ex_time[i];
        (*sig)[i]      = ex_zeta[i];
    }

    /*	4.) Convert the TS */

    /*
    ConvertTS_LGMSV_to_LGM(	nex,
                                                    ex_time,
                                                    *sig,
                                                    *lambda,
                                                    ModelParams.Tstar);
    */

    /*	5.)	Save instrument data if required */
    if (inst_data)
    {
        inst_data->num_inst      = nex;
        inst_data->start_dates   = (long*)calloc(nex, sizeof(long));
        inst_data->end_dates     = (long*)calloc(nex, sizeof(long));
        inst_data->short_strikes = (double*)calloc(nex, sizeof(double));
        inst_data->long_strikes  = (double*)calloc(nex, sizeof(double));

        if (!inst_data->start_dates || !inst_data->end_dates || !inst_data->short_strikes ||
            !inst_data->long_strikes)
        {
            err = "Allocation error (4) in cpd_calib_diagonal";
            goto FREE_RETURN;
        }

        for (i = 0; i < nex; i++)
        {
            inst_data->start_dates[i] = cpn_date[ex_cpn[i]];
            inst_data->end_dates[i]   = cpn_date[ncpn - 1];
            if (!fix_lambda)
            {
                inst_data->short_strikes[i] = ex_sstrike[i];
            }
            else
            {
                inst_data->short_strikes[i] = 0.0;
            }
            inst_data->long_strikes[i] = ex_lstrike[i];
        }
    }

FREE_RETURN:

    if (err)
    {
        if (*sig_time)
            free(*sig_time);
        *sig_time = NULL;

        if (*sig)
            free(*sig);
        *sig = NULL;

        if (inst_data)
        {
            cpd_free_calib_inst_data(inst_data);
        }
    }

    return err;
}
