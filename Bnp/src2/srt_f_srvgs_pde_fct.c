#include "opfnctns.h"
#include "math.h"
#include "pde_h_struct.h"
#include "srt_h_all.h"
#include "srt_h_fwdcurve.h"
#include "srt_h_ts_eq.h"

#define LIM_IN_STDEV 6

static Err compute_diffusion_generators(
    Date           ex_date,
    SrtUndPtr      und,
    SrtUndPtr      dom_und,
    SrtBasicPdeInf pde_inf,
    double****     diff_mat,
    double****     conv_mat);

static Err compute_black_scholes_greeks(
    double          fwd,
    double          strike,
    double          volatility,
    double          time,
    double          df,
    SrtReceiverType rec_pay,

    double* gamma,
    double* vega,
    double* volga,
    double* vanna);

static Err make_greeks_allocations(
    SrtBasicPdeInf pde_inf,
    double***      gamma_matrix,
    double***      vega_matrix,
    double***      volga_matrix,
    double***      vanna_matrix)
{
    (*gamma_matrix) = dmatrix(0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*vega_matrix)  = dmatrix(0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*volga_matrix) = dmatrix(0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*vanna_matrix) = dmatrix(0, pde_inf.num_step, 1, pde_inf.num_mesh);

    if (((*gamma_matrix) == NULL) || ((*vega_matrix) == NULL) || ((*volga_matrix) == NULL) ||
        ((*vanna_matrix) == NULL))

        return serror("Allocation Memory Failure");

    else
        return NULL;
}

static Err make_greeks_desallocations(
    SrtBasicPdeInf pde_inf,
    double***      gamma_matrix,
    double***      vega_matrix,
    double***      volga_matrix,
    double***      vanna_matrix)
{
    if ((*gamma_matrix))
        free_dmatrix((*gamma_matrix), 0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*gamma_matrix) = NULL;
    if ((*vega_matrix))
        free_dmatrix((*vega_matrix), 0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*vega_matrix) = NULL;
    if ((*vanna_matrix))
        free_dmatrix((*vanna_matrix), 0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*vanna_matrix) = NULL;
    if ((*volga_matrix))
        free_dmatrix((*volga_matrix), 0, pde_inf.num_step, 1, pde_inf.num_mesh);
    (*volga_matrix) = NULL;

    return NULL;
}

static Err allocate_pde_inf_ptr(SrtBasicPdeInf* pde_inf)
{
    long num_mesh, num_time_step;

    num_mesh      = (*pde_inf).num_mesh;
    num_time_step = (*pde_inf).num_step;

    (*pde_inf).log_spots = dvector(1, num_mesh - 1);
    (*pde_inf).spots     = dvector(1, num_mesh - 1);

    (*pde_inf).time = dvector(0, num_time_step);

    (*pde_inf).basevol = dvector(0, num_time_step);
    (*pde_inf).beta    = dvector(0, num_time_step);
    (*pde_inf).gamma   = dvector(0, num_time_step);
    (*pde_inf).omega   = dvector(0, num_time_step);

    (*pde_inf).voldrift = dvector(0, num_time_step);
    (*pde_inf).vovol    = dvector(0, num_time_step);

    (*pde_inf).correlation_und_vol  = dvector(0, num_time_step);
    (*pde_inf).correlation_und_rate = dvector(0, num_time_step);

    (*pde_inf).df      = dvector(0, num_time_step);
    (*pde_inf).forward = dvector(0, num_time_step);

    return NULL;
}

Err srt_f_set_srvgs_pde_operators(
    Date           ex_date,
    String         dom_und_name,
    String         und_name,
    SrtBasicPdeInf pde_inf,

    double          strike,
    SrtReceiverType rec_pay,
    double****      diff_mat,
    double****      conv_mat,
    double****      div_tensor)
{
    Err           err = NULL;
    TermStruct *  dom_ts, *ts;
    SrtUndPtr     dom_und, und;
    SrtCorrLstPtr corr_list;

    Date   date_today;
    long   mesh_index, step_index;
    String discount_crv_name;

    double **gamma_matrix, **vega_matrix, **volga_matrix, **vanna_matrix;

    double basevol, beta, gamma, omega, voldrift, vovol, correlation_und_vol, correlation_und_rate,
        time, stdev, tanh_stdev, skew, wing;

    /* Get the underlying & its term structure */
    und = lookup_und(und_name);

    err = get_underlying_ts(und, &ts);
    if (err)
        return err;

    /* Get the correlation list */
    corr_list = srt_f_GetTheCorrelationList();
    if (!corr_list)
        return serror("Cannot get the correlation structure");

    /* Get the domestic underlying term structure */
    dom_und = lookup_und(dom_und_name);

    err = get_underlying_ts(dom_und, &dom_ts);
    if (err)
        return err;

    /* Get today from und */
    date_today = get_today_from_underlying(und);

    /* Get the discount curve for the domestic underlying */
    discount_crv_name = get_discname_from_underlying(dom_und);

    /* make the greeks allocations */
    err =
        make_greeks_allocations(pde_inf, &gamma_matrix, &vega_matrix, &volga_matrix, &vanna_matrix);
    if (err)
        return err;

    err = compute_diffusion_generators(ex_date, und, dom_und, pde_inf, diff_mat, conv_mat);
    if (err)
        return err;

    for (step_index = 0; step_index <= pde_inf.num_step; step_index++)
    {
        /* Get the (only) time dependent variables */
        time = pde_inf.time[step_index];

        beta  = pde_inf.beta[step_index];
        gamma = pde_inf.gamma[step_index];
        omega = pde_inf.omega[step_index];

        basevol  = pde_inf.basevol[step_index];
        voldrift = pde_inf.voldrift[step_index];
        vovol    = pde_inf.vovol[step_index];

        correlation_und_vol  = pde_inf.correlation_und_vol[step_index];
        correlation_und_rate = pde_inf.correlation_und_rate[step_index];

        for (mesh_index = 1; mesh_index < pde_inf.num_mesh;
             mesh_index++) /* loop on each mesh point */
        {
            /* compute the greeks */
            err = compute_black_scholes_greeks(
                pde_inf.spots[mesh_index] * (pde_inf.forward[0] / pde_inf.forward[step_index]),
                strike,
                basevol,
                pde_inf.time[0] - pde_inf.time[step_index],
                pde_inf.df[0] / pde_inf.df[step_index],

                rec_pay,
                &gamma_matrix[step_index][mesh_index],
                &vega_matrix[step_index][mesh_index],
                &volga_matrix[step_index][mesh_index],
                &vanna_matrix[step_index][mesh_index]);

            if (err)
                return err;

            stdev      = log(pde_inf.spots[mesh_index] / pde_inf.spot) / (basevol * sqrt(time));
            tanh_stdev = omega * tanh(stdev / omega);

            skew = beta * tanh_stdev;
            wing = gamma * tanh_stdev * tanh_stdev;

            gamma_matrix[step_index][mesh_index] *=
                pow(pde_inf.forward[0] / pde_inf.forward[step_index], 2);
            vega_matrix[step_index][mesh_index] *= 100;
            vanna_matrix[step_index][mesh_index] *=
                (pde_inf.forward[0] / pde_inf.forward[step_index]);

            /* cost of hedge in presence of stochastic volatility */
            (*div_tensor)[0][step_index][mesh_index] =
                -voldrift * vega_matrix[step_index][mesh_index] -
                gamma_matrix[step_index][mesh_index] * pde_inf.spots[mesh_index] *
                    pde_inf.spots[mesh_index] * basevol * voldrift * time -
                0.5 * vovol * vovol * volga_matrix[step_index][mesh_index] -
                vovol * correlation_und_vol * basevol * pde_inf.spots[mesh_index] *
                    vanna_matrix[step_index][mesh_index];

            (*div_tensor)[0][step_index][mesh_index] *= (pde_inf.df[step_index] / pde_inf.df[0]);

            /* cost of hedge in presence of smile  */
            (*div_tensor)[1][step_index][mesh_index] =
                -gamma_matrix[step_index][mesh_index] * pde_inf.spots[mesh_index] *
                pde_inf.spots[mesh_index] * basevol * basevol * (skew + wing);
            (*div_tensor)[1][step_index][mesh_index] *= (pde_inf.df[step_index] / pde_inf.df[0]);

        } /* End of loop on mesh_index */
    }

    err = make_greeks_desallocations(
        pde_inf, &gamma_matrix, &vega_matrix, &volga_matrix, &vanna_matrix);
    if (err)
        return err;

    return NULL;
}

Err compute_black_scholes_greeks(
    double          fwd,
    double          strike,
    double          volatility,
    double          time_to_exp,
    double          df,
    SrtReceiverType rec_pay,

    double* gamma,
    double* vega,
    double* volga,
    double* vanna)

{
    *gamma = srt_f_optblksch(fwd, strike, volatility, time_to_exp, df, rec_pay, GAMMA);
    *vega  = srt_f_optblksch(fwd, strike, volatility, time_to_exp, df, rec_pay, VEGA);
    *volga = srt_f_optblksch(fwd, strike, volatility, time_to_exp, df, rec_pay, VOLGA);
    *vanna = srt_f_optblksch(fwd, strike, volatility, time_to_exp, df, rec_pay, VANNA);

    return NULL;
}

Err compute_diffusion_generators(
    Date           ex_date,
    SrtUndPtr      und,
    SrtUndPtr      dom_und,
    SrtBasicPdeInf pde_inf,
    double****     diff_mat,
    double****     conv_mat)
{
    TermStruct* ts;
    Date        date_today;
    String      yc_name;
    double      drift, stdev, tanh_stdev, spot, time;
    double      skew, wing;
    long        mesh_index, step_index;

    Date date;
    Err  err = NULL;

    err = get_underlying_ts(und, &ts);
    if (err)
        return err;

    yc_name = get_discname_from_underlying(dom_und);

    spot = get_spot_from_eqund(und);

    date_today = get_today_from_underlying(und);

    for (step_index = 0; step_index <= pde_inf.num_step; step_index++)
    {
        if (step_index == 0)
            time = (ex_date - date_today) * YEARS_IN_DAY;
        else
            time = max(pde_inf.time_step, pde_inf.time[step_index - 1] - pde_inf.time_step);

        pde_inf.time[step_index] = time;
        date                     = (Date)(date_today + time * DAYS_IN_YEAR);

        err = srt_f_eq_forward(date, und, yc_name, &pde_inf.forward[step_index]);

        if (err)
            return err;

        pde_inf.df[step_index] = swp_f_df(date_today, date, yc_name);

        pde_inf.beta[step_index]                = find_eq_beta(time, ts);
        pde_inf.gamma[step_index]               = find_eq_gamma(time, ts);
        pde_inf.omega[step_index]               = find_eq_omega(time, ts);
        pde_inf.basevol[step_index]             = find_eq_basevol(time, ts);
        pde_inf.voldrift[step_index]            = find_eq_voldrift(time, ts);
        pde_inf.vovol[step_index]               = find_eq_vovol(time, ts);
        pde_inf.correlation_und_vol[step_index] = find_eq_rho(time, ts);

        err = srt_f_get_corr_from_CorrList(
            srt_f_GetTheCorrelationList(),
            und->underl_name,
            dom_und->underl_name,
            time,
            &pde_inf.correlation_und_rate[step_index]);
        if (err)
            return err;
    }

    for (step_index = 0; step_index <= pde_inf.num_step; step_index++)
    {
        for (mesh_index = 1; mesh_index < pde_inf.num_mesh;
             mesh_index++) /* loop on each mesh point */
        {
            time = pde_inf.time[step_index];

            stdev =
                log(pde_inf.spots[mesh_index] / spot) / (pde_inf.basevol[step_index] * sqrt(time));
            tanh_stdev = tanh(stdev / pde_inf.omega[step_index]);

            skew = pde_inf.beta[step_index] * pde_inf.omega[step_index] * tanh_stdev;
            wing = pde_inf.gamma[step_index] * pde_inf.omega[step_index] *
                   pde_inf.omega[step_index] * tanh_stdev * tanh_stdev;

            (*diff_mat)[0][step_index][mesh_index] = pde_inf.basevol[step_index];
            (*diff_mat)[1][step_index][mesh_index] = pde_inf.basevol[step_index];

            if (step_index == pde_inf.num_step)
                drift = log(pde_inf.forward[step_index] / spot) / pde_inf.time_step;
            else
                drift = log(pde_inf.forward[step_index] / pde_inf.forward[step_index + 1]) /
                        pde_inf.time_step;

            (*conv_mat)[0][step_index][mesh_index] =
                drift - 0.5 * pde_inf.basevol[step_index] * pde_inf.basevol[step_index];
            (*conv_mat)[1][step_index][mesh_index] = (*conv_mat)[0][step_index][mesh_index];

        } /* End of loop on mesh_index */

    } /* End of loop on step index */

    return NULL;
}

Err und_pde_lim(
    Date            ex_date,
    SrtUndPtr       und,
    SrtUndPtr       dom_und,
    double          max_time,
    long            min_node,
    long            min_num_mesh,
    SrtBasicPdeInf* pde_inf)
{
    Err         err = NULL;
    SrtMdlType  dom_mdl_type, und_mdl_type;
    TermStruct* ts;
    String      yc_name;
    Date        date_today;
    long        num_mesh, k, i;
    double      log_spot_min, log_spot_max;
    double      black_vol, forward, spot, exp_mat;
    double      max_h, time_step, h;

    /* Check the interest rate and the equity underlying type */
    err = get_underlying_mdltype(dom_und, &dom_mdl_type);
    if (err)
        return err;

    err = get_underlying_mdltype(und, &und_mdl_type);
    if (err)
        return err;

    if ((dom_mdl_type != LGM) || (und_mdl_type != EQ_STOCH_RATES_SRVGS))
        return serror("Wrong type of model ");

    /* Get the discount curve and the date today*/
    yc_name = get_discname_from_underlying(dom_und);
    if (!yc_name)
        return serror("can not find yield curve");

    date_today = get_today_from_underlying(und);

    /* Get the equity underlying term structure */
    err = get_underlying_ts(und, &ts);
    if (err)
        return err;

    /* Compute the equity forward and its volatility  */

    spot = spot = get_spot_from_eqund(und);
    err         = srt_f_eq_forward(ex_date, und, yc_name, &forward);
    if (err)
        return err;

    exp_mat = (ex_date - date_today) * YEARS_IN_DAY;
    err     = srt_f_eq_implied_vol(exp_mat, und->underl_name, &black_vol);
    if (err)
        return err;

    log_spot_min =
        log(forward) - 0.5 * black_vol * black_vol * exp_mat - LIM_IN_STDEV * black_vol * exp_mat;

    log_spot_max =
        log(forward) - 0.5 * black_vol * black_vol * exp_mat + LIM_IN_STDEV * black_vol * exp_mat;

    /* Build the discretisation scheme */
    max_h = black_vol * black_vol /
            fabs(swp_f_zr(ex_date, ex_date + 14.0, yc_name) - 0.5 * black_vol * black_vol);

    h = (log_spot_max - log_spot_min) / min_num_mesh;
    h = min(max_h, h);

    /* adjust h so that log_spot is a mesh point */
    k        = (long)((log(spot) - log_spot_min) / h) + 1;
    h        = (log(spot) - log_spot_min) / k;
    num_mesh = (long)((log_spot_max - log_spot_min) / h);

    /* get the time step */
    time_step = min(max_time, (ex_date - date_today) * YEARS_IN_DAY / min_node);

    /* fill the pde inf structure */
    (*pde_inf).h         = h;
    (*pde_inf).time_step = time_step;

    (*pde_inf).num_step = (long)((ex_date - date_today) * YEARS_IN_DAY / time_step);

    (*pde_inf).num_mesh = num_mesh;

    (*pde_inf).log_spot_min = log_spot_min;
    (*pde_inf).spot_index   = k;

    err = allocate_pde_inf_ptr(pde_inf);
    if (err)
        return err;

    for (i = 1; i < num_mesh; i++)
    {
        (*pde_inf).log_spots[i] = log_spot_min + i * h;
        (*pde_inf).spots[i]     = exp((*pde_inf).log_spots[i]);
    }

    (*pde_inf).num_pde = 2;

    (*pde_inf).spot     = spot;
    (*pde_inf).spot_min = exp(log_spot_min);
    (*pde_inf).spot_max = exp(log_spot_max);

    return (NULL);
}

#undef LIM_IN_STDEV
