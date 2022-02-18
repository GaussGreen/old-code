/**********************************************************************
 *      Name: SrtGrfnMain.c                                           *
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 18/10/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 18/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "BGMEval.h"
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

#define MAX_STP 3000

static double H_func(double lam, double t1, double t2)
{
    return exp(lam * t2) - exp(lam * t1);
}

static void LGM2FExpectations(
    int     nstept,
    double* time,
    double  lam1,
    double* sig_time,
    double* sig1,
    int     nb_sig,
    double  alpha,
    double  gamma,
    double  rho,
    double* phi1,
    double* phi2,
    double* phi12)
{
    double lam2, vol2, lam21, lam22, lam12;
    double t1, t2, ta, tb;
    double coef11, coef12, coef13, coef21, coef22, coef23, fact1, fact2, rhoalpha;

    double H1, H2, H21, H22, H12;
    double st1;

    int i, j, nb_sig_minus1;

    lam2          = lam1 + gamma;
    lam21         = 2.0 * lam1;
    lam22         = 2.0 * lam2;
    lam12         = (lam1 + lam2);
    rhoalpha      = rho * alpha;
    nb_sig_minus1 = nb_sig - 1;

    fact1 = 1.0 / alpha / sqrt(1.0 - rho * rho);
    fact2 = rho / sqrt(1.0 - rho * rho);

    /* constant for expectation r1_dim1 */
    coef11 = (lam2 + rho * alpha * lam1) / (lam1 * lam1 * lam2);
    coef12 = -1.0 / (2.0 * lam1 * lam1);
    coef13 = -rhoalpha / (lam2 * lam12);

    /* constant for expectation r3 */
    coef21 = alpha * (rho * lam2 + alpha * lam1) / (lam1 * lam2 * lam2);
    coef22 = -alpha * alpha / (2.0 * lam2 * lam2);
    coef23 = -rhoalpha / (lam1 * lam12);

    /* initialisation */
    H1 = H2 = H21 = H22 = H12 = 0.0;
    t1                        = 0.0;
    j                         = 0;
    vol2                      = sig1[0] * sig1[0];

    phi1[0]  = 0.0;
    phi2[0]  = 0.0;
    phi12[0] = 0.0;

    for (i = 1; i < nstept; i++)
    {
        t2 = time[i];
        ta = t1;
        tb = sig_time[j];

        st1 = 0.0;

        while (tb < t2 && j < nb_sig_minus1)
        {
            H1 += vol2 * H_func(lam1, ta, tb);
            H2 += vol2 * H_func(lam2, ta, tb);
            H21 += vol2 * H_func(lam21, ta, tb);
            H22 += vol2 * H_func(lam22, ta, tb);
            H12 += vol2 * H_func(lam12, ta, tb);
            st1 += vol2 * (tb - ta);

            j++;
            vol2 = sig1[j] * sig1[j];
            ta   = tb;
            tb   = sig_time[j];
        }

        H1 += vol2 * H_func(lam1, ta, t2);
        H2 += vol2 * H_func(lam2, ta, t2);
        H21 += vol2 * H_func(lam21, ta, t2);
        H22 += vol2 * H_func(lam22, ta, t2);
        H12 += vol2 * H_func(lam12, ta, t2);
        st1 += vol2 * (t2 - ta);

        phi1[i]  = exp(-lam21 * t2) * H21;
        phi2[i]  = exp(-lam22 * t2) * H22;
        phi12[i] = exp(-lam12 * t2) * H12;

        phi1[i] /= lam21;
        phi2[i] *= alpha * alpha / lam22;
        phi12[i] *= rho * alpha / lam12;

        t1 = t2;
    }
}

char* SrtGrfn5DFXMc(
    char*     und3dfx,
    double**  correl,
    int       numeventdates,
    long*     eventdates,
    long      tableauRows,
    long*     tableauCols,
    char***   tableauStrings,
    int**     tableauMask,
    long      auxWidth,
    long*     auxLen,
    double**  aux,
    long      num_paths,
    int       do_pecs,
    double*** prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;

    double *time = NULL, *date = NULL;

    double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL;

    double *dom_ifr = NULL, *dom_fwd1 = NULL, *dom_fwd2 = NULL, *dom_exp1 = NULL, *dom_exp2 = NULL,
           *dom_phi1 = NULL, *dom_phi2 = NULL, *dom_phi12 = NULL, *dom_gam1_fwd = NULL,
           *dom_gam2_fwd = NULL, *dom_bond_pay = NULL, *dom_gam1_pay = NULL, *dom_gam2_pay = NULL,
           *for_ifr = NULL, *for_fwd1 = NULL, *for_fwd2 = NULL, *for_exp1 = NULL, *for_exp2 = NULL,
           *for_phi1 = NULL, *for_phi2 = NULL, *for_phi12 = NULL, *for_gam1_fwd = NULL,
           *for_gam2_fwd = NULL, *fx_fwd = NULL, ***covar = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC5F grfn_prm;

    long        today, spot_date;
    int         i, j, k, num_col;
    SrtUndPtr   fx_und, dom_und, for_und;
    TermStruct *fx_ts, *dom_ts, *for_ts;
    char *      domname, *forname;
    double lam_dom, lam_dom2, tau_dom, alpha_dom, gamma_dom, rho_dom, lam_for, lam_for2, tau_for,
        alpha_for, gamma_for, rho_for, correl_dom_for, correl_dom_fx, correl_for_fx;
    double spot_fx;
    char * dom_yc, *for_yc;
    int    fx_idx, dom_idx, for_idx;

    double *sigma_date_dom = NULL, *sigma_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
           *sigma_date_fx = NULL, *sigma_fx = NULL, *merge_dates = NULL;

    double *eigen_val = NULL, **eigen_vec = NULL, **corr_temp = NULL;

    int nrot;

    double pay_time, df;

    long nb_merge_dates, pay_date;

    long sigma_n_dom, sigma_n_for, sigma_n_fx;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /* First we check that the correlation matrix has positive eigen values */

    eigen_val = dvector(0, 4);
    eigen_vec = dmatrix(0, 4, 0, 4);
    corr_temp = dmatrix(0, 4, 0, 4);

    if (!eigen_val || !eigen_vec || !corr_temp)
    {
        err = "Memory allocation error in SrtGrfn5DFXMcEB";
        goto FREE_RETURN;
    }

    for (i = 0; i < 5; i++)
    {
        corr_temp[i][i] = 1.0;

        for (j = i + 1; j < 5; j++)
        {
            corr_temp[i][j] = corr_temp[j][i] = correl[i][j];
        }
    }

    err = jacobi_diagonalisation(corr_temp, 5, eigen_val, eigen_vec, &nrot);

    if (err || eigen_val[0] < 0.0 || eigen_val[1] < 0.0 || eigen_val[2] < 0.0 ||
        eigen_val[3] < 0.0 || eigen_val[4] < 0.0)
    {
        if (!err)
        {
            err = "Correlation matrix is not a positive matrix";
        }

        goto FREE_RETURN;
    }

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        und3dfx,
        &defParm,
        &forback,
        &xStr);

    if (err)
    {
        goto FREE_RETURN;
    }

    free_str = 1;

    /*	Now, lookup underlyings involved and their term structures */

    fx_und = lookup_und(und3dfx);

    if (!fx_und)
    {
        err = serror("Couldn't find underlying named %s", und3dfx);
        goto FREE_RETURN;
    }

    today = get_today_from_underlying(fx_und);

    if (get_underlying_type(fx_und) != FOREX_UND)
    {
        err = serror("Underlying %s is not of type FX", und3dfx);
        goto FREE_RETURN;
    }

    if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES)
    {
        err = serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
        goto FREE_RETURN;
    }

    fx_ts = get_ts_from_fxund(fx_und);

    domname = get_domname_from_fxund(fx_und);
    dom_und = lookup_und(domname);
    if (!dom_und)
    {
        err = serror("Couldn't find underlying named %s", domname);
        goto FREE_RETURN;
    }
    dom_ts = get_ts_from_irund(dom_und);

    forname = get_forname_from_fxund(fx_und);
    for_und = lookup_und(forname);
    if (!for_und)
    {
        err = serror("Couldn't find underlying named %s", forname);
        goto FREE_RETURN;
    }
    for_ts = get_ts_from_irund(for_und);

    num_col = xStr.num_cols;
    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }

    time = (double*)calloc(nstp, sizeof(double));
    date = (double*)calloc(nstp, sizeof(double));

    if (!time || !date)
    {
        err = "Memory allocation error (1) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    memcpy(time, xStr.tms, nstp * sizeof(double));

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    pay_date = (long)(date[nstp - 1] + 1.0E-8);
    pay_time = time[nstp - 1];

    /* Get all the term structures */
    err = Get_FX_StochRate_TermStructures5F(
        und3dfx,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_for,
        &alpha_for,
        &gamma_for,
        &rho_for,
        &sigma_date_fx,
        &sigma_fx,
        &sigma_n_fx,
        &correl_dom_for,
        &correl_dom_fx,
        &correl_for_fx);
    if (err)
    {
        goto FREE_RETURN;
    }

    rho_dom = correl[0][1];
    rho_for = correl[2][3];

    lam_dom  = 1.0 / tau_dom;
    lam_dom2 = lam_dom + gamma_dom;
    lam_for  = 1.0 / tau_for;
    lam_for2 = lam_for + gamma_for;

    /* now merge all the term structure */
    merge_dates = (double*)calloc(sigma_n_dom, sizeof(double));
    memcpy(merge_dates, sigma_date_dom, sigma_n_dom * sizeof(double));
    nb_merge_dates = sigma_n_dom;
    num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_for, sigma_date_for);
    num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
    num_f_sort_vector(nb_merge_dates, merge_dates);
    num_f_unique_vector(&nb_merge_dates, merge_dates);

    /*	Fill the new term structures */
    sig_dom = (double*)calloc(nb_merge_dates, sizeof(double));
    sig_for = (double*)calloc(nb_merge_dates, sizeof(double));
    sig_fx  = (double*)calloc(nb_merge_dates, sizeof(double));

    if (!sig_dom || !sig_for || !sig_fx)
    {
        err = "Memory allocation error (2) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    for (i = nb_merge_dates - 1; i >= 0; i--)
    {
        sig_dom[i] = sigma_dom[Get_Index(merge_dates[i], sigma_date_dom, sigma_n_dom)];
        sig_for[i] = sigma_for[Get_Index(merge_dates[i], sigma_date_for, sigma_n_for)];
        sig_fx[i]  = sigma_fx[Get_Index(merge_dates[i], sigma_date_fx, sigma_n_fx)];
    }

    /*	Get Fx spot and yield curves */
    dom_yc = get_ycname_from_irund(dom_und);
    for_yc = get_ycname_from_irund(for_und);

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    spot_fx = get_spot_from_fxund(fx_und) * swp_f_df(today, spot_date, dom_yc) /
              swp_f_df(today, spot_date, for_yc);

    /*	Get distributions */

    dom_ifr      = (double*)calloc(nstp, sizeof(double));
    dom_fwd1     = (double*)calloc(nstp, sizeof(double));
    dom_fwd2     = (double*)calloc(nstp, sizeof(double));
    dom_exp1     = (double*)calloc(nstp, sizeof(double));
    dom_exp2     = (double*)calloc(nstp, sizeof(double));
    dom_phi1     = (double*)calloc(nstp, sizeof(double));
    dom_phi2     = (double*)calloc(nstp, sizeof(double));
    dom_phi12    = (double*)calloc(nstp, sizeof(double));
    dom_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    dom_gam2_fwd = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam1_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam2_pay = (double*)calloc(nstp, sizeof(double));

    for_ifr      = (double*)calloc(nstp, sizeof(double));
    for_fwd1     = (double*)calloc(nstp, sizeof(double));
    for_fwd2     = (double*)calloc(nstp, sizeof(double));
    for_exp1     = (double*)calloc(nstp, sizeof(double));
    for_exp2     = (double*)calloc(nstp, sizeof(double));
    for_phi1     = (double*)calloc(nstp, sizeof(double));
    for_phi2     = (double*)calloc(nstp, sizeof(double));
    for_phi12    = (double*)calloc(nstp, sizeof(double));
    for_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    for_gam2_fwd = (double*)calloc(nstp, sizeof(double));

    fx_fwd = (double*)calloc(nstp, sizeof(double));

    covar = f3tensor(0, nstp - 1, 0, 4, 0, 4);

    if (!dom_ifr || !dom_fwd1 || !dom_fwd2 || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_exp1 ||
        !dom_exp2 || !dom_gam1_fwd || !dom_gam2_fwd || !dom_bond_pay || !dom_gam1_pay ||
        !dom_gam2_pay || !for_ifr || !for_fwd1 || !for_fwd2 || !for_phi1 || !for_phi2 ||
        !for_phi12 || !for_exp1 || !for_exp2 || !for_gam1_fwd || !for_gam2_fwd || !fx_fwd || !covar)
    {
        err = "Memory allocation error (3) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    fill_mc_init5F(
        pay_date,
        pay_time,
        date,
        time,
        nstp,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        lam_dom,
        alpha_dom,
        gamma_dom,
        sig_for,
        lam_for,
        alpha_for,
        gamma_for,
        sig_fx,
        correl,
        dom_yc,
        for_yc,
        dom_ifr,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        for_ifr,
        for_fwd1,
        for_fwd2,
        for_exp1,
        for_exp2,
        for_phi1,
        for_phi2,
        for_phi12,
        for_gam1_fwd,
        for_gam2_fwd,
        fx_fwd,
        covar);

    /*	Fill product structure */

    strupper(und3dfx);
    strip_white_space(und3dfx);
    strupper(domname);
    strip_white_space(domname);
    strupper(forname);
    strip_white_space(forname);
    for (i = 0; i < xStr.num_und; i++)
    {
        strupper(xStr.und_data[i].und_name);
        strip_white_space(xStr.und_data[i].und_name);
    }

    fx_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, und3dfx))
        {
            fx_idx = i;
        }
    }
    if (fx_idx == -1)
    {
        err = "The Fx underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    dom_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, domname))
        {
            dom_idx = i;
        }
    }
    if (dom_idx == -1)
    {
        err = "The domestic underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    for_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, forname))
        {
            for_idx = i;
        }
    }
    if (for_idx == -1)
    {
        err = "The foreign underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (4) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    for (i = xStr.num_evt - 1; i >= 0; i--)
    {
        if (xStr.evt[i].evt)
        {
            grfn_prm          = malloc(sizeof(grfn_parm_mc5F));
            grfn_prm->global  = &xStr;
            grfn_prm->local   = xStr.evt + i;
            grfn_prm->fx_idx  = fx_idx;
            grfn_prm->dom_idx = dom_idx;
            grfn_prm->for_idx = for_idx;

            grfn_prm->num_fx_df = xStr.evt[i].evt->dflen[fx_idx];
            grfn_prm->fx_df_tms = xStr.evt[i].evt->dft[fx_idx];
            grfn_prm->fx_df_dts = xStr.evt[i].evt->dfd[fx_idx];

            if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
            {
                grfn_prm->fx_dff   = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam1  = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam2  = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam12 = dvector(0, grfn_prm->num_fx_df - 1);

                if (!grfn_prm->fx_dff || !grfn_prm->fx_gam1 || !grfn_prm->fx_gam2 ||
                    !grfn_prm->fx_gam12)
                {
                    err = "Memory allocation error (7) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_fx_df; k++)
                {
                    grfn_prm->fx_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->fx_df_dts[k], (char*)dom_yc);
                    grfn_prm->fx_gam1[k] = (1.0 - exp(-lam_dom * grfn_prm->fx_df_tms[k])) / lam_dom;
                    grfn_prm->fx_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->fx_df_tms[k])) / lam_dom2;
                    grfn_prm->fx_gam12[k] =
                        -0.5 * (grfn_prm->fx_gam1[k] * grfn_prm->fx_gam1[k] * dom_phi1[i + flag] +
                                grfn_prm->fx_gam2[k] * grfn_prm->fx_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->fx_gam1[k] * grfn_prm->fx_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_fx = 1;
            }
            else
            {
                grfn_prm->do_fx = 0;
            }

            grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    grfn_prm->dom_gam1[k] =
                        (1.0 - exp(-lam_dom * grfn_prm->dom_df_tms[k])) / lam_dom;
                    grfn_prm->dom_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->dom_df_tms[k])) / lam_dom2;
                    grfn_prm->dom_gam12[k] =
                        -0.5 *
                            (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i + flag] +
                             grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            grfn_prm->num_for_df = xStr.evt[i].evt->dflen[for_idx];
            grfn_prm->for_df_tms = xStr.evt[i].evt->dft[for_idx];
            grfn_prm->for_df_dts = xStr.evt[i].evt->dfd[for_idx];

            if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
            {
                grfn_prm->for_dff   = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam1  = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam2  = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam12 = dvector(0, grfn_prm->num_for_df - 1);

                if (!grfn_prm->for_dff || !grfn_prm->for_gam1 || !grfn_prm->for_gam2 ||
                    !grfn_prm->for_gam12)
                {
                    err = "Memory allocation error (9) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_for_df; k++)
                {
                    grfn_prm->for_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->for_df_dts[k], (char*)for_yc);
                    grfn_prm->for_gam1[k] =
                        (1.0 - exp(-lam_for * grfn_prm->for_df_tms[k])) / lam_for;
                    grfn_prm->for_gam2[k] =
                        (1.0 - exp(-lam_for2 * grfn_prm->for_df_tms[k])) / lam_for2;
                    grfn_prm->for_gam12[k] =
                        -0.5 *
                            (grfn_prm->for_gam1[k] * grfn_prm->for_gam1[k] * for_phi1[i + flag] +
                             grfn_prm->for_gam2[k] * grfn_prm->for_gam2[k] * for_phi2[i + flag]) -
                        grfn_prm->for_gam1[k] * grfn_prm->for_gam2[k] * for_phi12[i + flag];
                }

                grfn_prm->do_for = 1;
            }
            else
            {
                grfn_prm->do_for = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *prod_val = dmatrix(0, num_col - 1, 0, 1);

    if (!(*prod_val))
    {
        err = "Memory allocation error";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_5dfx(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        dom_ifr,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        for_ifr,
        for_fwd1,
        for_fwd2,
        for_exp1,
        for_exp2,
        for_phi1,
        for_phi2,
        for_phi12,
        for_gam1_fwd,
        for_gam2_fwd,
        fx_fwd,
        covar,
        void_prm,
        /*	Model data */
        lam_dom,
        alpha_dom,
        gamma_dom,
        lam_for,
        alpha_for,
        gamma_for,
        correl,
        /*	Market data */
        spot_fx,
        dom_yc,
        for_yc,
        /* Do PECS adjustment */
        do_pecs,
        0,
        NULL,
        NULL,
        NULL,
        /*	Payoff function */
        grfn_payoff_4_5dfx_mc, /*	Result */
        *prod_val);

    df = swp_f_zr(today, pay_date, dom_yc);
    df = exp(-df * pay_time);

    *tableauCols = num_col;
    for (i = 0; i < num_col; i++)
    {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
    }

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (eigen_val)
        free_dvector(eigen_val, 0, 4);
    if (eigen_vec)
        free_dmatrix(eigen_vec, 0, 4, 0, 4);
    if (corr_temp)
        free_dmatrix(corr_temp, 0, 4, 0, 4);

    if (time)
        free(time);
    if (date)
        free(date);
    if (merge_dates)
        free(merge_dates);
    if (sig_dom)
        free(sig_dom);
    if (sig_for)
        free(sig_for);
    if (sig_fx)
        free(sig_fx);
    if (dom_ifr)
        free(dom_ifr);
    if (dom_fwd1)
        free(dom_fwd1);
    if (dom_fwd2)
        free(dom_fwd2);
    if (dom_exp1)
        free(dom_exp1);
    if (dom_exp2)
        free(dom_exp2);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_gam1_fwd)
        free(dom_gam1_fwd);
    if (dom_gam2_fwd)
        free(dom_gam2_fwd);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_gam1_pay)
        free(dom_gam1_pay);
    if (dom_gam2_pay)
        free(dom_gam2_pay);

    if (for_ifr)
        free(for_ifr);
    if (for_fwd1)
        free(for_fwd1);
    if (for_fwd2)
        free(for_fwd2);
    if (for_exp1)
        free(for_exp1);
    if (for_exp2)
        free(for_exp2);
    if (for_phi1)
        free(for_phi1);
    if (for_phi2)
        free(for_phi2);
    if (for_phi12)
        free(for_phi12);
    if (for_gam1_fwd)
        free(for_gam1_fwd);
    if (for_gam2_fwd)
        free(for_gam2_fwd);

    if (fx_fwd)
        free(fx_fwd);

    if (covar)
        free_f3tensor(covar, 0, nstp - 1, 0, 4, 0, 4);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);

    if (sigma_date_for)
        free(sigma_date_for);
    if (sigma_for)
        free(sigma_for);

    if (sigma_date_fx)
        free(sigma_date_fx);
    if (sigma_fx)
        free(sigma_fx);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC5F)void_prm[i];

                if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
                {
                    if (grfn_prm->fx_dff)
                        free_dvector(grfn_prm->fx_dff, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam1)
                        free_dvector(grfn_prm->fx_gam1, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam2)
                        free_dvector(grfn_prm->fx_gam2, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam12)
                        free_dvector(grfn_prm->fx_gam12, 0, grfn_prm->num_fx_df - 1);
                }

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                if (grfn_prm->do_for && grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
                {
                    if (grfn_prm->for_dff)
                        free_dvector(grfn_prm->for_dff, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam1)
                        free_dvector(grfn_prm->for_gam1, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam2)
                        free_dvector(grfn_prm->for_gam2, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam12)
                        free_dvector(grfn_prm->for_gam12, 0, grfn_prm->num_for_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}

char* SrtGrfn5DFXMcEB(
    char*      und3dfx,
    double**   correl,
    int        numeventdates,
    long*      eventdates,
    int*       nUsedEventDates,
    int*       optimise,
    double*    fwd_iv,
    MCEBPARAMS params,
    long*      resRows,
    long       tableauRows,
    long*      tableauCols,
    char***    tableauStrings,
    int**      tableauMask,
    long       auxWidth,
    long*      auxLen,
    double**   aux,
    long       num_paths,
    int        do_pecs,
    double***  prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;
    int*         ivOptimise = 0;

    double *time = NULL, *date = NULL;

    double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL;

    double *dom_ifr = NULL, *dom_fwd1 = NULL, *dom_fwd2 = NULL, *dom_exp1 = NULL, *dom_exp2 = NULL,
           *dom_phi1 = NULL, *dom_phi2 = NULL, *dom_phi12 = NULL, *dom_gam1_fwd = NULL,
           *dom_gam2_fwd = NULL, *dom_bond_pay = NULL, *dom_gam1_pay = NULL, *dom_gam2_pay = NULL,
           *for_ifr = NULL, *for_fwd1 = NULL, *for_fwd2 = NULL, *for_exp1 = NULL, *for_exp2 = NULL,
           *for_phi1 = NULL, *for_phi2 = NULL, *for_phi12 = NULL, *for_gam1_fwd = NULL,
           *for_gam2_fwd = NULL, *fx_fwd = NULL, ***covar = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC5F grfn_prm;

    long        today, spot_date;
    int         i, j, k, num_col;
    SrtUndPtr   fx_und, dom_und, for_und;
    TermStruct *fx_ts, *dom_ts, *for_ts;
    char *      domname, *forname;
    double lam_dom, lam_dom2, tau_dom, alpha_dom, gamma_dom, rho_dom, lam_for, lam_for2, tau_for,
        alpha_for, gamma_for, rho_for, correl_dom_for, correl_dom_fx, correl_for_fx;
    double spot_fx;
    char * dom_yc, *for_yc;
    int    fx_idx, dom_idx, for_idx;

    double *sigma_date_dom = NULL, *sigma_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
           *sigma_date_fx = NULL, *sigma_fx = NULL, *merge_dates = NULL;

    double *eigen_val = NULL, **eigen_vec = NULL, **corr_temp = NULL;

    int nrot;

    double pay_time, df;

    long nb_merge_dates, pay_date;

    long sigma_n_dom, sigma_n_for, sigma_n_fx;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /* First we check that the correlation matrix has positive eigen values */

    eigen_val = dvector(0, 4);
    eigen_vec = dmatrix(0, 4, 0, 4);
    corr_temp = dmatrix(0, 4, 0, 4);

    if (!eigen_val || !eigen_vec || !corr_temp)
    {
        err = "Memory allocation error in SrtGrfn5DFXMcEB";
        goto FREE_RETURN;
    }

    for (i = 0; i < 5; i++)
    {
        corr_temp[i][i] = 1.0;

        for (j = i + 1; j < 5; j++)
        {
            corr_temp[i][j] = corr_temp[j][i] = correl[i][j];
        }
    }

    err = jacobi_diagonalisation(corr_temp, 5, eigen_val, eigen_vec, &nrot);

    if (err || eigen_val[0] < 0.0 || eigen_val[1] < 0.0 || eigen_val[2] < 0.0 ||
        eigen_val[3] < 0.0 || eigen_val[4] < 0.0)
    {
        if (!err)
        {
            err = "Correlation matrix is not a positive matrix";
        }

        goto FREE_RETURN;
    }

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        und3dfx,
        &defParm,
        &forback,
        &xStr);

    if (err)
    {
        goto FREE_RETURN;
    }

    free_str = 1;

    /*	Now, lookup underlyings involved and their term structures */

    fx_und = lookup_und(und3dfx);

    if (!fx_und)
    {
        err = serror("Couldn't find underlying named %s", und3dfx);
        goto FREE_RETURN;
    }

    today = get_today_from_underlying(fx_und);

    if (get_underlying_type(fx_und) != FOREX_UND)
    {
        err = serror("Underlying %s is not of type FX", und3dfx);
        goto FREE_RETURN;
    }

    if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES)
    {
        err = serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
        goto FREE_RETURN;
    }

    fx_ts = get_ts_from_fxund(fx_und);

    domname = get_domname_from_fxund(fx_und);
    dom_und = lookup_und(domname);
    if (!dom_und)
    {
        err = serror("Couldn't find underlying named %s", domname);
        goto FREE_RETURN;
    }
    dom_ts = get_ts_from_irund(dom_und);

    forname = get_forname_from_fxund(fx_und);
    for_und = lookup_und(forname);
    if (!for_und)
    {
        err = serror("Couldn't find underlying named %s", forname);
        goto FREE_RETURN;
    }
    for_ts = get_ts_from_irund(for_und);

    num_col = xStr.num_cols;
    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }

    time = (double*)calloc(nstp, sizeof(double));
    date = (double*)calloc(nstp, sizeof(double));

    if (!time || !date)
    {
        err = "Memory allocation error (1) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    memcpy(time, xStr.tms, nstp * sizeof(double));

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    pay_date = (long)(date[nstp - 1] + 1.0E-8);
    pay_time = time[nstp - 1];

    /* Get all the term structures */
    err = Get_FX_StochRate_TermStructures5F(
        und3dfx,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_for,
        &alpha_for,
        &gamma_for,
        &rho_for,
        &sigma_date_fx,
        &sigma_fx,
        &sigma_n_fx,
        &correl_dom_for,
        &correl_dom_fx,
        &correl_for_fx);
    if (err)
    {
        goto FREE_RETURN;
    }

    rho_dom = correl[0][1];
    rho_for = correl[2][3];

    lam_dom  = 1.0 / tau_dom;
    lam_dom2 = lam_dom + gamma_dom;
    lam_for  = 1.0 / tau_for;
    lam_for2 = lam_for + gamma_for;

    /* now merge all the term structure */
    merge_dates = (double*)calloc(sigma_n_dom, sizeof(double));
    memcpy(merge_dates, sigma_date_dom, sigma_n_dom * sizeof(double));
    nb_merge_dates = sigma_n_dom;
    num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_for, sigma_date_for);
    num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
    num_f_sort_vector(nb_merge_dates, merge_dates);
    num_f_unique_vector(&nb_merge_dates, merge_dates);

    /*	Fill the new term structures */
    sig_dom = (double*)calloc(nb_merge_dates, sizeof(double));
    sig_for = (double*)calloc(nb_merge_dates, sizeof(double));
    sig_fx  = (double*)calloc(nb_merge_dates, sizeof(double));

    if (!sig_dom || !sig_for || !sig_fx)
    {
        err = "Memory allocation error (2) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    for (i = nb_merge_dates - 1; i >= 0; i--)
    {
        sig_dom[i] = sigma_dom[Get_Index(merge_dates[i], sigma_date_dom, sigma_n_dom)];
        sig_for[i] = sigma_for[Get_Index(merge_dates[i], sigma_date_for, sigma_n_for)];
        sig_fx[i]  = sigma_fx[Get_Index(merge_dates[i], sigma_date_fx, sigma_n_fx)];
    }

    /*	Get Fx spot and yield curves */
    dom_yc = get_ycname_from_irund(dom_und);
    for_yc = get_ycname_from_irund(for_und);

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    spot_fx = get_spot_from_fxund(fx_und) * swp_f_df(today, spot_date, dom_yc) /
              swp_f_df(today, spot_date, for_yc);

    /*	Get distributions */

    dom_ifr      = (double*)calloc(nstp, sizeof(double));
    dom_fwd1     = (double*)calloc(nstp, sizeof(double));
    dom_fwd2     = (double*)calloc(nstp, sizeof(double));
    dom_exp1     = (double*)calloc(nstp, sizeof(double));
    dom_exp2     = (double*)calloc(nstp, sizeof(double));
    dom_phi1     = (double*)calloc(nstp, sizeof(double));
    dom_phi2     = (double*)calloc(nstp, sizeof(double));
    dom_phi12    = (double*)calloc(nstp, sizeof(double));
    dom_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    dom_gam2_fwd = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam1_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam2_pay = (double*)calloc(nstp, sizeof(double));

    for_ifr      = (double*)calloc(nstp, sizeof(double));
    for_fwd1     = (double*)calloc(nstp, sizeof(double));
    for_fwd2     = (double*)calloc(nstp, sizeof(double));
    for_exp1     = (double*)calloc(nstp, sizeof(double));
    for_exp2     = (double*)calloc(nstp, sizeof(double));
    for_phi1     = (double*)calloc(nstp, sizeof(double));
    for_phi2     = (double*)calloc(nstp, sizeof(double));
    for_phi12    = (double*)calloc(nstp, sizeof(double));
    for_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    for_gam2_fwd = (double*)calloc(nstp, sizeof(double));

    fx_fwd = (double*)calloc(nstp, sizeof(double));

    covar = f3tensor(0, nstp - 1, 0, 4, 0, 4);

    if (!dom_ifr || !dom_fwd1 || !dom_fwd2 || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_exp1 ||
        !dom_exp2 || !dom_gam1_fwd || !dom_gam2_fwd || !dom_bond_pay || !dom_gam1_pay ||
        !dom_gam2_pay || !for_ifr || !for_fwd1 || !for_fwd2 || !for_phi1 || !for_phi2 ||
        !for_phi12 || !for_exp1 || !for_exp2 || !for_gam1_fwd || !for_gam2_fwd || !fx_fwd || !covar)
    {
        err = "Memory allocation error (3) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    fill_mc_init5F(
        pay_date,
        pay_time,
        date,
        time,
        nstp,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        lam_dom,
        alpha_dom,
        gamma_dom,
        sig_for,
        lam_for,
        alpha_for,
        gamma_for,
        sig_fx,
        correl,
        dom_yc,
        for_yc,
        dom_ifr,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        for_ifr,
        for_fwd1,
        for_fwd2,
        for_exp1,
        for_exp2,
        for_phi1,
        for_phi2,
        for_phi12,
        for_gam1_fwd,
        for_gam2_fwd,
        fx_fwd,
        covar);

    /*	Fill product structure */

    strupper(und3dfx);
    strip_white_space(und3dfx);
    strupper(domname);
    strip_white_space(domname);
    strupper(forname);
    strip_white_space(forname);
    for (i = 0; i < xStr.num_und; i++)
    {
        strupper(xStr.und_data[i].und_name);
        strip_white_space(xStr.und_data[i].und_name);
    }

    fx_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, und3dfx))
        {
            fx_idx = i;
        }
    }
    if (fx_idx == -1)
    {
        err = "The Fx underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    dom_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, domname))
        {
            dom_idx = i;
        }
    }
    if (dom_idx == -1)
    {
        err = "The domestic underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    for_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, forname))
        {
            for_idx = i;
        }
    }
    if (for_idx == -1)
    {
        err = "The foreign underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (4) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    for (i = xStr.num_evt - 1; i >= 0; i--)
    {
        if (xStr.evt[i].evt)
        {
            grfn_prm          = malloc(sizeof(grfn_parm_mc5F));
            grfn_prm->global  = &xStr;
            grfn_prm->local   = xStr.evt + i;
            grfn_prm->fx_idx  = fx_idx;
            grfn_prm->dom_idx = dom_idx;
            grfn_prm->for_idx = for_idx;

            grfn_prm->num_fx_df = xStr.evt[i].evt->dflen[fx_idx];
            grfn_prm->fx_df_tms = xStr.evt[i].evt->dft[fx_idx];
            grfn_prm->fx_df_dts = xStr.evt[i].evt->dfd[fx_idx];

            if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
            {
                grfn_prm->fx_dff   = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam1  = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam2  = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam12 = dvector(0, grfn_prm->num_fx_df - 1);

                if (!grfn_prm->fx_dff || !grfn_prm->fx_gam1 || !grfn_prm->fx_gam2 ||
                    !grfn_prm->fx_gam12)
                {
                    err = "Memory allocation error (7) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_fx_df; k++)
                {
                    grfn_prm->fx_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->fx_df_dts[k], (char*)dom_yc);
                    grfn_prm->fx_gam1[k] = (1.0 - exp(-lam_dom * grfn_prm->fx_df_tms[k])) / lam_dom;
                    grfn_prm->fx_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->fx_df_tms[k])) / lam_dom2;
                    grfn_prm->fx_gam12[k] =
                        -0.5 * (grfn_prm->fx_gam1[k] * grfn_prm->fx_gam1[k] * dom_phi1[i + flag] +
                                grfn_prm->fx_gam2[k] * grfn_prm->fx_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->fx_gam1[k] * grfn_prm->fx_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_fx = 1;
            }
            else
            {
                grfn_prm->do_fx = 0;
            }

            grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    grfn_prm->dom_gam1[k] =
                        (1.0 - exp(-lam_dom * grfn_prm->dom_df_tms[k])) / lam_dom;
                    grfn_prm->dom_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->dom_df_tms[k])) / lam_dom2;
                    grfn_prm->dom_gam12[k] =
                        -0.5 *
                            (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i + flag] +
                             grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            grfn_prm->num_for_df = xStr.evt[i].evt->dflen[for_idx];
            grfn_prm->for_df_tms = xStr.evt[i].evt->dft[for_idx];
            grfn_prm->for_df_dts = xStr.evt[i].evt->dfd[for_idx];

            if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
            {
                grfn_prm->for_dff   = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam1  = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam2  = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam12 = dvector(0, grfn_prm->num_for_df - 1);

                if (!grfn_prm->for_dff || !grfn_prm->for_gam1 || !grfn_prm->for_gam2 ||
                    !grfn_prm->for_gam12)
                {
                    err = "Memory allocation error (9) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_for_df; k++)
                {
                    grfn_prm->for_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->for_df_dts[k], (char*)for_yc);
                    grfn_prm->for_gam1[k] =
                        (1.0 - exp(-lam_for * grfn_prm->for_df_tms[k])) / lam_for;
                    grfn_prm->for_gam2[k] =
                        (1.0 - exp(-lam_for2 * grfn_prm->for_df_tms[k])) / lam_for2;
                    grfn_prm->for_gam12[k] =
                        -0.5 *
                            (grfn_prm->for_gam1[k] * grfn_prm->for_gam1[k] * for_phi1[i + flag] +
                             grfn_prm->for_gam2[k] * grfn_prm->for_gam2[k] * for_phi2[i + flag]) -
                        grfn_prm->for_gam1[k] * grfn_prm->for_gam2[k] * for_phi12[i + flag];
                }

                grfn_prm->do_for = 1;
            }
            else
            {
                grfn_prm->do_for = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *tableauCols     = num_col;
    *resRows         = max(num_col + 1, nstp);
    *nUsedEventDates = xStr.num_evt;
    df               = swp_f_df(today, pay_date, dom_yc);

    /* create an optimisation vector from the input */
    ivOptimise = ivector(0, nstp - 1);
    if (flag)
        ivOptimise[0] = 0;
    for (i = flag; i < nstp; i++)
        ivOptimise[i] = optimise[i + numeventdates - xStr.num_evt - flag];

    if (params->iMultiIndex)
    {
        params->iNbIndex = params->iMultiIndex;
    }
    else
    {
        params->iNbIndex = 1;
    }

    err = mceb_allocate_params(params, nstp);

    if (err)
        goto FREE_RETURN;

    if (params->iAdjustIV)
    {
        if (flag)
        {
            params->dMarketFwdIV[0] = 0.0;
        }

        for (i = flag; i < nstp; i++)
        {
            params->dMarketFwdIV[i] = fwd_iv[i + numeventdates - xStr.num_evt - flag] / df;
        }
    }

    *prod_val = dmatrix(0, *resRows - 1, 0, 2 + params->iNbIndex);

    if (!(*prod_val))
    {
        err = "Memory allocation error";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_5dfx(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        dom_ifr,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        for_ifr,
        for_fwd1,
        for_fwd2,
        for_exp1,
        for_exp2,
        for_phi1,
        for_phi2,
        for_phi12,
        for_gam1_fwd,
        for_gam2_fwd,
        fx_fwd,
        covar,
        void_prm,
        /*	Model data */
        lam_dom,
        alpha_dom,
        gamma_dom,
        lam_for,
        alpha_for,
        gamma_for,
        correl,
        /*	Market data */
        spot_fx,
        dom_yc,
        for_yc,
        /* Do PECS adjustment */
        do_pecs,
        1,
        ivOptimise,
        params,
        NULL,
        /*	Payoff function */
        grfn_payoff_4_5dfx_mc, /*	Result */
        *prod_val);

    /* Recopy Barrier / CoefLin for the moment */
    for (i = 0; i < nstp; i++)
    {
        (*prod_val)[i][2] = params->dBarrier[i];

        for (j = 0; j < params->iNbIndex; j++)
        {
            (*prod_val)[i][3 + j] = params->dCoefLin[i][j + 1];
        }
    }

    df = swp_f_zr(today, pay_date, dom_yc);
    df = exp(-df * pay_time);

    for (i = 0; i < num_col + 1; i++)
    {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
    }

    for (i = 0; i < nstp; i++)
    {
        if (params->iCalcOneTime)
            params->dOneTimeCall[i] *= df;
        if (params->iCalcOneTimePartial)
            params->dOneTimePartial[i] *= df;
        if (params->iCalcIV || params->iAdjustIV)
            params->dModelFwdIV[i] *= df;
    }

    if (flag)
    {
        for (i = 0; i < nstp - 1; i++)
        {
            (*prod_val)[i][2] = (*prod_val)[i + 1][2];

            for (k = 0; k < params->iNbIndex; k++)
            {
                (*prod_val)[i][3 + k] = (*prod_val)[i + 1][3 + k];
            }
        }

        mceb_shift_extrainfos(params);
    }

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;
    (*prod_val)[num_col][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (eigen_val)
        free_dvector(eigen_val, 0, 4);
    if (eigen_vec)
        free_dmatrix(eigen_vec, 0, 4, 0, 4);
    if (corr_temp)
        free_dmatrix(corr_temp, 0, 4, 0, 4);

    if (time)
        free(time);
    if (date)
        free(date);
    if (merge_dates)
        free(merge_dates);
    if (sig_dom)
        free(sig_dom);
    if (sig_for)
        free(sig_for);
    if (sig_fx)
        free(sig_fx);
    if (dom_ifr)
        free(dom_ifr);
    if (dom_fwd1)
        free(dom_fwd1);
    if (dom_fwd2)
        free(dom_fwd2);
    if (dom_exp1)
        free(dom_exp1);
    if (dom_exp2)
        free(dom_exp2);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_gam1_fwd)
        free(dom_gam1_fwd);
    if (dom_gam2_fwd)
        free(dom_gam2_fwd);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_gam1_pay)
        free(dom_gam1_pay);
    if (dom_gam2_pay)
        free(dom_gam2_pay);

    if (for_ifr)
        free(for_ifr);
    if (for_fwd1)
        free(for_fwd1);
    if (for_fwd2)
        free(for_fwd2);
    if (for_exp1)
        free(for_exp1);
    if (for_exp2)
        free(for_exp2);
    if (for_phi1)
        free(for_phi1);
    if (for_phi2)
        free(for_phi2);
    if (for_phi12)
        free(for_phi12);
    if (for_gam1_fwd)
        free(for_gam1_fwd);
    if (for_gam2_fwd)
        free(for_gam2_fwd);

    if (fx_fwd)
        free(fx_fwd);

    if (covar)
        free_f3tensor(covar, 0, nstp - 1, 0, 4, 0, 4);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);

    if (sigma_date_for)
        free(sigma_date_for);
    if (sigma_for)
        free(sigma_for);

    if (sigma_date_fx)
        free(sigma_date_fx);
    if (sigma_fx)
        free(sigma_fx);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC5F)void_prm[i];

                if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
                {
                    if (grfn_prm->fx_dff)
                        free_dvector(grfn_prm->fx_dff, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam1)
                        free_dvector(grfn_prm->fx_gam1, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam2)
                        free_dvector(grfn_prm->fx_gam2, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam12)
                        free_dvector(grfn_prm->fx_gam12, 0, grfn_prm->num_fx_df - 1);
                }

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                if (grfn_prm->do_for && grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
                {
                    if (grfn_prm->for_dff)
                        free_dvector(grfn_prm->for_dff, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam1)
                        free_dvector(grfn_prm->for_gam1, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam2)
                        free_dvector(grfn_prm->for_gam2, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam12)
                        free_dvector(grfn_prm->for_gam12, 0, grfn_prm->num_for_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}

char* SrtGrfn5DFXMcTest(
    char*     und3dfx,
    double**  correl,
    int       numeventdates,
    long*     eventdates,
    long      tableauRows,
    long*     tableauCols,
    char***   tableauStrings,
    int**     tableauMask,
    long      auxWidth,
    long*     auxLen,
    double**  aux,
    long      num_dates,
    long      num_paths,
    int       do_pecs,
    double*** prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;

    double *time = NULL, *date = NULL;

    double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL;

    double *dom_ifr = NULL, *dom_fwd = NULL, *dom_std = NULL, *dom_phi1 = NULL, *dom_phi2 = NULL,
           *dom_phi12 = NULL, *dom_beta = NULL, *dom_bond_pay = NULL, *dom_beta_pay = NULL,

           *for_ifr = NULL, *for_fwd = NULL, *for_std = NULL, *for_phi1 = NULL, *for_phi2 = NULL,
           *for_phi12 = NULL, *for_beta = NULL,

           *fx_fwd = NULL, *fx_std = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC5F grfn_prm;

    long        today, spot_date;
    int         i, j, k, num_col;
    SrtUndPtr   fx_und, dom_und, for_und;
    TermStruct *fx_ts, *dom_ts, *for_ts;
    char *      domname, *forname;
    double lam_dom, lam_dom2, tau_dom, alpha_dom, gamma_dom, rho_dom, lam_for, lam_for2, tau_for,
        alpha_for, gamma_for, rho_for, correl_dom_for, correl_dom_fx, correl_for_fx;
    double spot_fx;
    char * dom_yc, *for_yc;
    int    fx_idx, dom_idx, for_idx;

    double *sigma_date_dom = NULL, *sigma_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
           *sigma_date_fx = NULL, *sigma_fx = NULL, *dom_for_cov = NULL, *dom_fx_cov = NULL,
           *for_fx_cov = NULL, *merge_dates = NULL;

    double *eigen_val = NULL, **eigen_vec = NULL, **corr_temp = NULL;

    int nrot;

    double pay_time, df;

    long nb_merge_dates, pay_date;

    long sigma_n_dom, sigma_n_for, sigma_n_fx;

    double sqdt;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /* First we check that the correlation matrix has positive eigen values */

    eigen_val = dvector(0, 4);
    eigen_vec = dmatrix(0, 4, 0, 4);
    corr_temp = dmatrix(0, 4, 0, 4);

    if (!eigen_val || !eigen_vec || !corr_temp)
    {
        err = "Memory allocation error in SrtGrfn5DFXMcEB";
        goto FREE_RETURN;
    }

    for (i = 0; i < 5; i++)
    {
        corr_temp[i][i] = 1.0;

        for (j = i + 1; j < 5; j++)
        {
            corr_temp[i][j] = corr_temp[j][i] = correl[i][j];
        }
    }

    err = jacobi_diagonalisation(corr_temp, 5, eigen_val, eigen_vec, &nrot);

    if (err || eigen_val[0] < 0.0 || eigen_val[1] < 0.0 || eigen_val[2] < 0.0 ||
        eigen_val[3] < 0.0 || eigen_val[4] < 0.0)
    {
        if (!err)
        {
            err = "Correlation matrix is not a positive matrix";
        }

        goto FREE_RETURN;
    }

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        und3dfx,
        &defParm,
        &forback,
        &xStr);

    if (err)
    {
        goto FREE_RETURN;
    }

    free_str = 1;

    /*	Now, lookup underlyings involved and their term structures */

    fx_und = lookup_und(und3dfx);

    if (!fx_und)
    {
        err = serror("Couldn't find underlying named %s", und3dfx);
        goto FREE_RETURN;
    }

    today = get_today_from_underlying(fx_und);

    if (get_underlying_type(fx_und) != FOREX_UND)
    {
        err = serror("Underlying %s is not of type FX", und3dfx);
        goto FREE_RETURN;
    }

    if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES)
    {
        err = serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
        goto FREE_RETURN;
    }

    fx_ts = get_ts_from_fxund(fx_und);

    domname = get_domname_from_fxund(fx_und);
    dom_und = lookup_und(domname);
    if (!dom_und)
    {
        err = serror("Couldn't find underlying named %s", domname);
        goto FREE_RETURN;
    }
    dom_ts = get_ts_from_irund(dom_und);

    forname = get_forname_from_fxund(fx_und);
    for_und = lookup_und(forname);
    if (!for_und)
    {
        err = serror("Couldn't find underlying named %s", forname);
        goto FREE_RETURN;
    }
    for_ts = get_ts_from_irund(for_und);

    num_col = xStr.num_cols;
    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }
    time = (double*)calloc(nstp, sizeof(double));
    if (!time)
    {
        err = "Memory allocation error (1) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
    }
    memcpy(time, xStr.tms, nstp * sizeof(double));

    /*	Fill the time vector */

    err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, num_dates);

    date = (double*)calloc(nstp, sizeof(double));
    if (!date)
    {
        err = "Memory allocation error (3) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
    }

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    pay_date = (long)(date[nstp - 1] + 1.0E-8);
    pay_time = time[nstp - 1];

    /* Get all the term structures */
    err = Get_FX_StochRate_TermStructures5F(
        und3dfx,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_for,
        &alpha_for,
        &gamma_for,
        &rho_for,
        &sigma_date_fx,
        &sigma_fx,
        &sigma_n_fx,
        &correl_dom_for,
        &correl_dom_fx,
        &correl_for_fx);
    if (err)
    {
        goto FREE_RETURN;
    }

    rho_dom = correl[0][1];
    rho_for = correl[2][3];

    lam_dom  = 1.0 / tau_dom;
    lam_dom2 = lam_dom + gamma_dom;
    lam_for  = 1.0 / tau_for;
    lam_for2 = lam_for + gamma_for;

    /* now merge all the term structure */
    merge_dates = (double*)calloc(sigma_n_dom, sizeof(double));
    memcpy(merge_dates, sigma_date_dom, sigma_n_dom * sizeof(double));
    nb_merge_dates = sigma_n_dom;
    num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_for, sigma_date_for);
    num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
    num_f_sort_vector(nb_merge_dates, merge_dates);
    num_f_unique_vector(&nb_merge_dates, merge_dates);

    /*	Fill the new term structures */

    sig_dom = (double*)calloc(nb_merge_dates, sizeof(double));
    sig_for = (double*)calloc(nb_merge_dates, sizeof(double));
    sig_fx  = (double*)calloc(nb_merge_dates, sizeof(double));

    if (!sig_dom || !sig_for || !sig_fx)
    {
        err = "Memory allocation error (4) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    for (i = nb_merge_dates - 1; i >= 0; i--)
    {
        sig_dom[i] = sigma_dom[Get_Index(merge_dates[i], sigma_date_dom, sigma_n_dom)];
        sig_for[i] = sigma_for[Get_Index(merge_dates[i], sigma_date_for, sigma_n_for)];
        sig_fx[i]  = sigma_fx[Get_Index(merge_dates[i], sigma_date_fx, sigma_n_fx)];
    }

    /*	Get Fx spot and yield curves */

    dom_yc = get_ycname_from_irund(dom_und);
    for_yc = get_ycname_from_irund(for_und);

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    spot_fx = get_spot_from_fxund(fx_und) * swp_f_df(today, spot_date, dom_yc) /
              swp_f_df(today, spot_date, for_yc);

    /*	Get distributions */

    dom_ifr   = (double*)calloc(nstp, sizeof(double));
    dom_fwd   = (double*)calloc(nstp, sizeof(double));
    dom_std   = (double*)calloc(nstp, sizeof(double));
    dom_phi1  = (double*)calloc(nstp, sizeof(double));
    dom_phi2  = (double*)calloc(nstp, sizeof(double));
    dom_phi12 = (double*)calloc(nstp, sizeof(double));

    dom_beta     = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_beta_pay = (double*)calloc(nstp, sizeof(double));

    for_ifr = (double*)calloc(nstp, sizeof(double));
    for_fwd = (double*)calloc(nstp, sizeof(double));
    for_std = (double*)calloc(nstp, sizeof(double));

    for_phi1  = (double*)calloc(nstp, sizeof(double));
    for_phi2  = (double*)calloc(nstp, sizeof(double));
    for_phi12 = (double*)calloc(nstp, sizeof(double));

    for_beta    = (double*)calloc(nstp, sizeof(double));
    fx_fwd      = (double*)calloc(nstp, sizeof(double));
    fx_std      = (double*)calloc(nstp, sizeof(double));
    dom_for_cov = (double*)calloc(nstp, sizeof(double));
    dom_fx_cov  = (double*)calloc(nstp, sizeof(double));
    for_fx_cov  = (double*)calloc(nstp, sizeof(double));

    if (!dom_ifr || !dom_fwd || !dom_std || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_beta ||
        !dom_bond_pay || !dom_beta_pay || !for_ifr || !for_fwd || !for_std || !for_phi1 ||
        !for_phi2 || !for_phi12 || !for_beta || !fx_fwd || !fx_std)
    {
        err = "Memory allocation error (5) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    LGM2FExpectations(
        nstp,
        time,
        lam_dom,
        merge_dates,
        sig_dom,
        nb_merge_dates,
        alpha_dom,
        gamma_dom,
        rho_dom,
        dom_phi1,
        dom_phi2,
        dom_phi12);

    LGM2FExpectations(
        nstp,
        time,
        lam_for,
        merge_dates,
        sig_for,
        nb_merge_dates,
        alpha_for,
        gamma_for,
        rho_for,
        for_phi1,
        for_phi2,
        for_phi12);

    for (i = 1; i < nstp; i++)
    {
        sqdt       = sqrt(time[i] - time[i - 1]);
        dom_std[i] = sigma_dom[Get_Index(time[i], sigma_date_dom, sigma_n_dom)] * sqdt;
        for_std[i] = sigma_for[Get_Index(time[i], sigma_date_for, sigma_n_for)] * sqdt;
        fx_std[i]  = sigma_fx[Get_Index(time[i], sigma_date_fx, sigma_n_fx)] * sqdt;
        fx_fwd[i]  = -0.5 * fx_std[i] * fx_std[i];

        dom_ifr[i - 1] = swp_f_zr(date[i - 1], date[i], dom_yc);
        for_ifr[i - 1] = swp_f_zr(date[i - 1], date[i], for_yc);
    }

    /*	Fill product structure */

    strupper(und3dfx);
    strip_white_space(und3dfx);
    strupper(domname);
    strip_white_space(domname);
    strupper(forname);
    strip_white_space(forname);
    for (i = 0; i < xStr.num_und; i++)
    {
        strupper(xStr.und_data[i].und_name);
        strip_white_space(xStr.und_data[i].und_name);
    }

    fx_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, und3dfx))
        {
            fx_idx = i;
        }
    }
    if (fx_idx == -1)
    {
        err = "The Fx underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    dom_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, domname))
        {
            dom_idx = i;
        }
    }
    if (dom_idx == -1)
    {
        err = "The domestic underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    for_idx = -1;
    for (i = 0; i < xStr.num_und; i++)
    {
        if (!strcmp(xStr.und_data[i].und_name, forname))
        {
            for_idx = i;
        }
    }
    if (for_idx == -1)
    {
        err = "The foreign underlying is not present in the mdlcomm structure";
        goto FREE_RETURN;
    }

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (6) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    j = xStr.num_evt - 1;

    for (i = nstp - 1; i >= 0; i--)
    {
        if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04)
        {
            grfn_prm          = malloc(sizeof(grfn_parm_mc5F));
            grfn_prm->global  = &xStr;
            grfn_prm->local   = xStr.evt + j;
            grfn_prm->fx_idx  = fx_idx;
            grfn_prm->dom_idx = dom_idx;
            grfn_prm->for_idx = for_idx;

            grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
            grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
            grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

            if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
            {
                grfn_prm->fx_dff   = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam1  = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam2  = dvector(0, grfn_prm->num_fx_df - 1);
                grfn_prm->fx_gam12 = dvector(0, grfn_prm->num_fx_df - 1);

                if (!grfn_prm->fx_dff || !grfn_prm->fx_gam1 || !grfn_prm->fx_gam2 ||
                    !grfn_prm->fx_gam12)
                {
                    err = "Memory allocation error (7) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_fx_df; k++)
                {
                    grfn_prm->fx_dff[k] =
                        swp_f_df(xStr.dts[j], grfn_prm->fx_df_dts[k], (char*)dom_yc);
                    grfn_prm->fx_gam1[k] = (1.0 - exp(-lam_dom * grfn_prm->fx_df_tms[k])) / lam_dom;
                    grfn_prm->fx_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->fx_df_tms[k])) / lam_dom2;
                    grfn_prm->fx_gam12[k] =
                        -0.5 * (grfn_prm->fx_gam1[k] * grfn_prm->fx_gam1[k] * dom_phi1[i] +
                                grfn_prm->fx_gam2[k] * grfn_prm->fx_gam2[k] * dom_phi2[i]) -
                        grfn_prm->fx_gam1[k] * grfn_prm->fx_gam2[k] * dom_phi12[i];
                }

                grfn_prm->do_fx = 1;
            }
            else
            {
                grfn_prm->do_fx = 0;
            }

            grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[j], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    grfn_prm->dom_gam1[k] =
                        (1.0 - exp(-lam_dom * grfn_prm->dom_df_tms[k])) / lam_dom;
                    grfn_prm->dom_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->dom_df_tms[k])) / lam_dom2;
                    grfn_prm->dom_gam12[k] =
                        -0.5 * (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i] +
                                grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i];
                }

                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
            grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
            grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

            if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
            {
                grfn_prm->for_dff   = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam1  = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam2  = dvector(0, grfn_prm->num_for_df - 1);
                grfn_prm->for_gam12 = dvector(0, grfn_prm->num_for_df - 1);

                if (!grfn_prm->for_dff || !grfn_prm->for_gam1 || !grfn_prm->for_gam2 ||
                    !grfn_prm->for_gam12)
                {
                    err = "Memory allocation error (9) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_for_df; k++)
                {
                    grfn_prm->for_dff[k] =
                        swp_f_df(xStr.dts[j], grfn_prm->for_df_dts[k], (char*)for_yc);
                    grfn_prm->for_gam1[k] =
                        (1.0 - exp(-lam_for * grfn_prm->for_df_tms[k])) / lam_for;
                    grfn_prm->for_gam2[k] =
                        (1.0 - exp(-lam_for2 * grfn_prm->for_df_tms[k])) / lam_for2;
                    grfn_prm->for_gam12[k] =
                        -0.5 * (grfn_prm->for_gam1[k] * grfn_prm->for_gam1[k] * for_phi1[i] +
                                grfn_prm->for_gam2[k] * grfn_prm->for_gam2[k] * for_phi2[i]) -
                        grfn_prm->for_gam1[k] * grfn_prm->for_gam2[k] * for_phi12[i];
                }

                grfn_prm->do_for = 1;
            }
            else
            {
                grfn_prm->do_for = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;

            j--;
            while (j >= 0 && xStr.evt[j].evt == NULL)
            {
                j--;
            }
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *prod_val = dmatrix(0, num_col - 1, 0, 1);

    if (!(*prod_val))
    {
        err = "Memory allocation error";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_5dfx_test(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        dom_ifr, /*	Distributions */
        dom_fwd,
        dom_std,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_beta,
        dom_bond_pay,
        dom_beta_pay,
        for_ifr,
        for_fwd,
        for_std,
        for_phi1,
        for_phi2,
        for_phi12,
        for_beta,
        fx_fwd,
        fx_std,
        dom_for_cov,
        dom_fx_cov,
        for_fx_cov,
        void_prm,
        /*	Model data */
        lam_dom,
        alpha_dom,
        gamma_dom,
        rho_dom,
        lam_for,
        alpha_for,
        gamma_for,
        rho_for,
        correl,
        /*	Market data */
        spot_fx,
        dom_yc,
        for_yc,
        /* Do PECS adjustment */
        do_pecs,
        NULL,
        /*	Payoff function */
        grfn_payoff_4_5dfx_mc, /*	Result */
        *prod_val);

    df = swp_f_zr(today, pay_date, dom_yc);
    df = exp(-df * pay_time);

    *tableauCols = num_col;

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (eigen_val)
        free_dvector(eigen_val, 0, 4);
    if (eigen_vec)
        free_dmatrix(eigen_vec, 0, 4, 0, 4);
    if (corr_temp)
        free_dmatrix(corr_temp, 0, 4, 0, 4);

    if (time)
        free(time);
    if (date)
        free(date);
    if (merge_dates)
        free(merge_dates);
    if (sig_dom)
        free(sig_dom);
    if (sig_for)
        free(sig_for);
    if (sig_fx)
        free(sig_fx);
    if (dom_ifr)
        free(dom_ifr);
    if (dom_fwd)
        free(dom_fwd);

    if (dom_std)
        free(dom_std);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_beta)
        free(dom_beta);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_beta_pay)
        free(dom_beta_pay);

    if (for_ifr)
        free(for_ifr);
    if (for_fwd)
        free(for_fwd);
    if (for_std)
        free(for_std);
    if (for_phi1)
        free(for_phi1);
    if (for_phi2)
        free(for_phi2);
    if (for_phi12)
        free(for_phi12);
    if (for_beta)
        free(for_beta);

    if (fx_fwd)
        free(fx_fwd);
    if (fx_std)
        free(fx_std);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);

    if (sigma_date_for)
        free(sigma_date_for);
    if (sigma_for)
        free(sigma_for);

    if (sigma_date_fx)
        free(sigma_date_fx);
    if (sigma_fx)
        free(sigma_fx);

    if (dom_for_cov)
        free(dom_for_cov);
    if (dom_fx_cov)
        free(dom_fx_cov);
    if (for_fx_cov)
        free(for_fx_cov);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC5F)void_prm[i];

                if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
                {
                    if (grfn_prm->fx_dff)
                        free_dvector(grfn_prm->fx_dff, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam1)
                        free_dvector(grfn_prm->fx_gam1, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam2)
                        free_dvector(grfn_prm->fx_gam2, 0, grfn_prm->num_fx_df - 1);
                    if (grfn_prm->fx_gam12)
                        free_dvector(grfn_prm->fx_gam12, 0, grfn_prm->num_fx_df - 1);
                }

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                if (grfn_prm->do_for && grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
                {
                    if (grfn_prm->for_dff)
                        free_dvector(grfn_prm->for_dff, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam1)
                        free_dvector(grfn_prm->for_gam1, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam2)
                        free_dvector(grfn_prm->for_gam2, 0, grfn_prm->num_for_df - 1);
                    if (grfn_prm->for_gam12)
                        free_dvector(grfn_prm->for_gam12, 0, grfn_prm->num_for_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}

static double Z_Func(double x, double t1, double t2)
{
    return ((exp(x * t2) - exp(x * t1)) / x);
}

static double Phi_Func(double x, double T, double s, double t)
{
    double result;

    result = (exp(-x * (T - t)) - exp(-x * (T - s))) / x;

    return result;
}

static double Etha_Func(double x, double T, double s, double t)
{
    double result;

    result = (t - s - Phi_Func(x, T, s, t)) / x;

    return result;
}

static double B_Func(double x, double T, double s, double t)
{
    double result;

    result = -(t * t - s * s) / 2.0 +
             1.0 / x * ((t - 1.0 / x) * exp(-x * (T - t)) - (s - 1.0 / x) * exp(-x * (T - s)));

    return result;
}

static double Psi_Func(double x, double y, double T, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, T, s, t) - Phi_Func(y, T, s, t) + Phi_Func(x + y, T, s, t));

    return result;
}

static double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, Tx, s, t) - Phi_Func(y, Ty, s, t) +
              exp(-x * (Tx - Ty)) * Phi_Func(x + y, Ty, s, t));

    return result;
}

Err fill_mc_init5Fsimple(
    long    pay_date,
    double  pay_time,
    double* date,
    double* time,
    long    nb_dates,
    double* sig_dates,
    long    nb_sig_dates,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* sig_curve_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    char*   dom_yc,
    char*   for_yc,
    double* dom_ifr,
    double* dom_fwd,
    double* dom_std,
    double* dom_phi,
    double* dom_beta,
    double* dom_bond_pay,
    double* dom_beta_pay,
    double* for_ifr,
    double* for_fwd,
    double* for_std,
    double* for_phi,
    double* for_beta,
    double* fx_fwd,
    double* fx_std,
    double* dom_for_cov,
    double* dom_fx_cov,
    double* for_fx_cov)
{
    double sig_dom, sig_dom2;
    double sig_for, sig_for2;
    double sig_fx, sig_dom_for;
    double T1, T2, start_date, end_date, start_mat, end_mat;
    double var_dom, expect_dom;
    double var_for, expect_for;
    double expect_fx, QTexpect, fx_vol;
    double adj_fx_pay, adj_quanto, adj_pay_dom, adj_pay_for;
    double x_dom, y_dom, x_domfor;
    double x_for, y_for;
    double phi_dom, phi_for, zc_dom, zc_for, zc_pay;
    double mat, pay_mat, mat_pay;
    double correl12, correl13_1, correl13_2, correl13_3, correl23_1, correl23_2, correl23_3;
    int    i, k;
    long   StartIndex, EndIndex;
    Err    err = NULL;

    dom_fwd[0] = 0;
    dom_std[0] = 0;
    dom_phi[0] = 0;
    phi_dom    = 0;

    for_fwd[0] = 0;
    for_std[0] = 0;
    for_phi[0] = 0;
    phi_for    = 0;

    mat_pay = (pay_date - date[0]) / 365.0;

    start_date = date[0];
    start_mat  = time[0];
    StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);

    for (k = 0; k < nb_dates - 1; k++)
    {
        end_date = date[k + 1];
        end_mat  = time[k + 1];
        EndIndex = Get_Index(end_mat, sig_dates, nb_sig_dates);
        mat      = (end_mat - start_mat);

        var_dom     = 0;
        expect_dom  = 0;
        adj_pay_dom = 0;

        var_for     = 0;
        expect_for  = 0;
        adj_quanto  = 0;
        adj_pay_for = 0;

        dom_ifr[k] = swp_f_zr(start_date, start_date + 1, dom_yc);
        for_ifr[k] = swp_f_zr(start_date, start_date + 1, for_yc);

        zc_dom = swp_f_zr(start_date, end_date, dom_yc);
        zc_for = swp_f_zr(start_date, end_date, for_yc);

        dom_beta[k] = -1.0 / lda_dom * (1.0 - exp(-lda_dom * mat));
        for_beta[k] = -1.0 / lda_for * (1.0 - exp(-lda_for * mat));

        /* QTexpect of the log Fx */
        QTexpect = -mat * (zc_for - zc_dom) - 0.5 * (for_beta[k] * for_beta[k] * for_phi[k] -
                                                     dom_beta[k] * dom_beta[k] * dom_phi[k]);

        zc_pay          = swp_f_zr(start_date, pay_date, dom_yc);
        pay_mat         = (pay_date - start_date) / 365.0;
        dom_bond_pay[k] = exp(-zc_pay * pay_mat);

        /*	Implied Fx Vol */
        err = Fx3DtsImpliedVol(
            end_mat,
            start_mat,
            end_mat,
            sig_dates,
            nb_sig_dates,
            sig_curve_dom,
            lda_dom,
            sig_curve_for,
            lda_for,
            sig_dates,
            sig_curve_fx,
            nb_sig_dates,
            correl_dom_for,
            correl_dom_fx,
            correl_for_fx,
            &fx_vol);

        /*	Calculate expectation of the log Fx under Q-Tfix */
        expect_fx = QTexpect - 0.5 * fx_vol * fx_vol * mat;

        adj_fx_pay = 0;

        correl12   = 0;
        correl13_1 = correl13_2 = correl13_3 = 0;
        correl23_1 = correl23_2 = correl23_3 = 0;

        for (i = StartIndex; i < EndIndex + 1; i++)
        {
            if (i > StartIndex)
            {
                T1 = sig_dates[i - 1];
            }
            else
            {
                /* First part */
                T1 = start_mat;
            }

            if (i == EndIndex || StartIndex == EndIndex)
            {
                /* Last part */
                T2 = end_mat;
            }
            else
            {
                T2 = sig_dates[i];
            }

            sig_dom = sig_dom2 = sig_dom_for = sig_curve_dom[i];
            sig_dom2 *= sig_dom2;
            sig_for = sig_for2 = sig_curve_for[i];
            sig_for2 *= sig_for2;
            sig_dom_for *= sig_for;
            sig_fx = sig_curve_fx[i];

            x_dom = Z_Func(lda_dom, T1, T2);
            y_dom = Z_Func(2 * lda_dom, T1, T2);

            x_for = Z_Func(lda_for, T1, T2);
            y_for = Z_Func(2 * lda_for, T1, T2);

            x_domfor = Z_Func(lda_dom + lda_for, T1, T2);

            expect_dom += sig_dom2 * (Phi_Func(lda_dom, end_mat, T1, T2) -
                                      Phi_Func(2 * lda_dom, end_mat, T1, T2));

            /* Q-Tpay adjustment */
            adj_pay_dom += sig_dom2 * (x_dom - exp(-lda_dom * mat_pay) * y_dom);

            var_dom += sig_dom2 * y_dom;

            expect_for += sig_for2 * (Phi_Func(lda_for, end_mat, T1, T2) -
                                      Phi_Func(2 * lda_for, end_mat, T1, T2));

            /* Quanto and Q-T pay adjustment */
            adj_quanto += sig_for * sig_fx * x_for;
            adj_pay_for += sig_dom_for * (x_for - exp(-lda_dom * mat_pay) * x_domfor);

            var_for += sig_for2 * y_for;

            /* Fx adjustment to go to Qpay */
            adj_fx_pay +=
                correl_dom_fx * sig_dom * sig_fx *
                    (Etha_Func(lda_dom, end_mat, T1, T2) - Etha_Func(lda_dom, mat_pay, T1, T2)) -
                correl_dom_for * sig_dom_for *
                    (Psi_Func(lda_for, lda_dom, end_mat, T1, T2) -
                     Psi2_Func(lda_for, lda_dom, end_mat, mat_pay, T1, T2)) +
                sig_dom2 * (Psi_Func(lda_dom, lda_dom, end_mat, T1, T2) -
                            Psi2_Func(lda_dom, lda_dom, end_mat, mat_pay, T1, T2));

            /* Correlations */

            /* domestic / foreign */
            correl12 += sig_dom * sig_for * x_domfor;

            /* domestic / fx */
            correl13_1 += sig_dom2 * (x_dom - exp(-lda_dom * end_mat) * y_dom);
            correl13_2 += sig_dom_for * (x_dom - exp(-lda_for * end_mat) * x_domfor);
            correl13_3 += sig_dom * sig_fx * x_dom;

            /* foreign fx */
            correl23_1 += sig_for2 * (x_for - exp(-lda_for * end_mat) * y_for);
            correl23_2 += sig_dom_for * (x_for - exp(-lda_dom * end_mat) * x_domfor);
            correl23_3 += sig_for * sig_fx * x_for;
        }

        phi_dom += var_dom;
        phi_for += var_for;

        /* Forward and Standard deviation */
        dom_fwd[k + 1] =
            1.0 / lda_dom * (expect_dom - exp(-lda_dom * end_mat) * adj_pay_dom) +
            dom_phi[k] * exp(-lda_dom * mat) * Phi_Func(-lda_dom, start_mat, start_mat, end_mat);

        dom_std[k + 1] = exp(-lda_dom * end_mat) * sqrt(var_dom);
        for_fwd[k + 1] =
            1.0 / lda_for * expect_for +
            exp(-lda_for * end_mat) *
                (-correl_dom_for * adj_pay_for / lda_dom - correl_for_fx * adj_quanto) +
            for_phi[k] * exp(-lda_for * mat) * Phi_Func(-lda_for, start_mat, start_mat, end_mat);

        for_std[k + 1] = exp(-lda_for * end_mat) * sqrt(var_for);
        fx_fwd[k + 1]  = expect_fx + adj_fx_pay;
        fx_std[k + 1]  = fx_vol * sqrt(mat);

        /* Covariance*/
        dom_for_cov[k + 1] = correl_dom_for * correl12 * exp(-(lda_dom + lda_for) * end_mat);

        dom_fx_cov[k + 1] = exp(-lda_dom * end_mat) *
                            (correl13_1 / lda_dom - correl_dom_for * correl13_2 / lda_for +
                             correl_dom_fx * correl13_3);

        for_fx_cov[k + 1] = exp(-lda_for * end_mat) *
                            (-correl23_1 / lda_for + correl_dom_for * correl23_2 / lda_dom +
                             correl_for_fx * correl23_3);

        /* Phi */
        dom_phi[k + 1] = phi_dom * exp(-2 * lda_dom * end_mat);
        for_phi[k + 1] = phi_for * exp(-2 * lda_for * end_mat);

        start_date = end_date;
        start_mat  = end_mat;
        StartIndex = EndIndex;
    }

    dom_ifr[k]      = swp_f_zr(date[k], date[k] + 1, dom_yc);
    for_ifr[k]      = swp_f_zr(date[k], date[k] + 1, for_yc);
    zc_pay          = swp_f_zr(date[k], pay_date, dom_yc);
    pay_mat         = (pay_date - date[k]) / 365.0;
    dom_bond_pay[k] = exp(-zc_pay * pay_mat);
    dom_beta_pay[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));

    return err;
}

Err fill_mc_init5F(
    long      pay_date,
    double    pay_time,
    double*   date,
    double*   time,
    long      nb_dates,
    double*   sig_dates,
    long      nb_sig_dates,
    double*   sig_curve_dom,
    double    lda_dom,
    double    alpha_dom,
    double    gamma_dom,
    double*   sig_curve_for,
    double    lda_for,
    double    alpha_for,
    double    gamma_for,
    double*   sig_curve_fx,
    double**  correl,
    char*     dom_yc,
    char*     for_yc,
    double*   dom_ifr,
    double*   dom_fwd1,
    double*   dom_fwd2,
    double*   dom_exp1,
    double*   dom_exp2,
    double*   dom_phi1,
    double*   dom_phi2,
    double*   dom_phi12,
    double*   dom_gam1_fwd,
    double*   dom_gam2_fwd,
    double*   dom_bond_pay,
    double*   dom_gam1_pay,
    double*   dom_gam2_pay,
    double*   for_ifr,
    double*   for_fwd1,
    double*   for_fwd2,
    double*   for_exp1,
    double*   for_exp2,
    double*   for_phi1,
    double*   for_phi2,
    double*   for_phi12,
    double*   for_gam1_fwd,
    double*   for_gam2_fwd,
    double*   fx_fwd,
    double*** covar)
{
    double alpha_dom2, alpha_for2, alpha_domfor;
    double lda_dom2, sig_dom, sig_dom2;
    double lda_for2, sig_for, sig_for2;
    double sig_fx, sig_dom_fx, sig_for_fx, sig_dom_for;
    double T1, T2, start_date, end_date, start_mat, end_mat;

    double phif_dom1, phif_dom2, phif_dom11, phif_dom12, phif_dom21, phif_dom22;
    double phi_dom1, phi_dom2, phi_dom12;
    double phi_domfor11, phi_domfor12, phi_domfor21, phi_domfor22;
    double etha_dom1_Tend, etha_dom1_Tpay, etha_dom2_Tend, etha_dom2_Tpay;
    double expect_dom1, expect_dom2, expect_dom12, expect_dom21;
    double adj_pay_dom1, adj_pay_dom12, adj_pay_dom2, adj_pay_dom21;

    double phif_for1, phif_for2, phif_for11, phif_for12, phif_for21, phif_for22;
    double phi_for1, phi_for2, phi_for12;
    double expect_for1, expect_for2, expect_for12, expect_for21;
    double adj_pay_for1, adj_pay_for12, adj_pay_for2, adj_pay_for21;
    double adj_quanto_for1, adj_quanto_for2;

    double phif_domfor11, phif_domfor21, phif_domfor12, phif_domfor22;
    double psif_domfor12_Tend, psif_domfor12_Tpay, psif_domfor21_Tend, psif_domfor21_Tpay,
        psif_domfor11_Tend, psif_domfor11_Tpay, psif_domfor22_Tend, psif_domfor22_Tpay,
        psif_domdom11_Tend, psif_domdom11_Tpay, psif_domdom12_Tend, psif_domdom12_Tpay,
        psif_domdom21_Tend, psif_domdom21_Tpay, psif_domdom22_Tend, psif_domdom22_Tpay;

    double adj_pay_fx1, adj_pay_fx2, adj_pay_fx3, adj_pay_fx4, adj_pay_fx5, adj_pay_fx6,
        adj_pay_fx7, adj_pay_fx8, adj_pay_fx9, adj_pay_fx10;

    double covar_dom1_fx1, covar_dom1_fx2, covar_dom1_fx3, covar_dom1_fx4, covar_dom1_fx5,
        covar_dom2_fx1, covar_dom2_fx2, covar_dom2_fx3, covar_dom2_fx4, covar_dom2_fx5,
        covar_for1_fx1, covar_for1_fx2, covar_for1_fx3, covar_for1_fx4, covar_for1_fx5,
        covar_for2_fx1, covar_for2_fx2, covar_for2_fx3, covar_for2_fx4, covar_for2_fx5;

    double mat_pay, mat, pay_mat, pay_mat2;
    double zc_for, zc_dom, zc_pay;
    double FxQTexpect, expect_fx, fx_vol;
    double exp_dom_mat, exp_dom_mat2, exp_dom_pay, exp_dom_pay2, exp_for_mat, exp_for_mat2,
        exp_for_pay, exp_for_pay2;
    int  i, k;
    long StartIndex, EndIndex;

    double** cov;

    Err err = NULL;

    alpha_dom2 = alpha_dom * alpha_dom;
    lda_dom2   = lda_dom + gamma_dom;

    dom_fwd1[0] = 0.0;
    dom_fwd2[0] = 0.0;

    dom_phi1[0]  = 0.0;
    dom_phi2[0]  = 0.0;
    dom_phi12[0] = 0.0;

    alpha_for2 = alpha_for * alpha_for;
    lda_for2   = lda_for + gamma_for;

    for_fwd1[0] = 0;
    for_fwd2[0] = 0;

    for_phi1[0]  = 0;
    for_phi2[0]  = 0;
    for_phi12[0] = 0;

    alpha_domfor = alpha_dom * alpha_for;

    mat_pay = (pay_date - date[0]) / 365.0;

    start_date = date[0];
    start_mat  = time[0];
    StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);

    for (k = 0; k < nb_dates - 1; k++)
    {
        end_date = date[k + 1];
        end_mat  = time[k + 1];
        EndIndex = Get_Index(end_mat, sig_dates, nb_sig_dates);
        mat      = (end_mat - start_mat);
        pay_mat  = (pay_date - start_date) / 365.0;
        pay_mat2 = (pay_date - end_date) / 365.0;

        /* calculation of the IFR and ZC */
        dom_ifr[k] = swp_f_zr(start_date, start_date + 1, dom_yc);
        for_ifr[k] = swp_f_zr(start_date, start_date + 1, for_yc);

        zc_dom = swp_f_zr(start_date, end_date, dom_yc);
        zc_for = swp_f_zr(start_date, end_date, for_yc);

        /* QTexpect of the log Fx */

        dom_gam1_fwd[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * mat));
        dom_gam2_fwd[k] = 1.0 / lda_dom2 * (1 - exp(-lda_dom2 * mat));

        for_gam1_fwd[k] = 1.0 / lda_for * (1 - exp(-lda_for * mat));
        for_gam2_fwd[k] = 1.0 / lda_for2 * (1 - exp(-lda_for2 * mat));

        FxQTexpect = -mat * (zc_for - zc_dom) -
                     0.5 * (for_gam1_fwd[k] * for_gam1_fwd[k] * for_phi1[k] +
                            for_gam2_fwd[k] * for_gam2_fwd[k] * for_phi2[k]) -
                     for_gam1_fwd[k] * for_gam2_fwd[k] * for_phi12[k] +
                     0.5 * (dom_gam1_fwd[k] * dom_gam1_fwd[k] * dom_phi1[k] +
                            dom_gam2_fwd[k] * dom_gam2_fwd[k] * dom_phi2[k]) +
                     dom_gam1_fwd[k] * dom_gam2_fwd[k] * dom_phi12[k];

        /* For Df to pay date */
        zc_pay = swp_f_zr(start_date, pay_date, dom_yc);

        dom_gam1_pay[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));
        dom_gam2_pay[k] = 1.0 / lda_dom2 * (1 - exp(-lda_dom2 * pay_mat));
        dom_bond_pay[k] =
            exp(-zc_pay * pay_mat -
                0.5 * (dom_gam1_pay[k] * dom_gam1_pay[k] * dom_phi1[k] +
                       dom_gam2_pay[k] * dom_gam2_pay[k] * dom_phi2[k]) -
                dom_gam1_pay[k] * dom_gam2_pay[k] * dom_phi12[k]);

        /*	Implied Fx Vol */
        err = Fx5DtsImpliedVol(
            end_mat,
            start_mat,
            end_mat,
            sig_dates,
            nb_sig_dates,
            sig_curve_dom,
            lda_dom,
            alpha_dom,
            gamma_dom,
            sig_curve_for,
            lda_for,
            alpha_for,
            gamma_for,
            sig_dates,
            sig_curve_fx,
            nb_sig_dates,
            correl,
            &fx_vol);

        /*	Calculate expectation of the log Fx under Q-Tfix */
        expect_fx = FxQTexpect - 0.5 * fx_vol * fx_vol * mat;

        /*	Variables initialisation */
        exp_dom_mat  = exp(-lda_dom * mat);
        exp_dom_mat2 = exp(-lda_dom2 * mat);
        exp_dom_pay  = exp(-lda_dom * pay_mat2);
        exp_dom_pay2 = exp(-lda_dom2 * pay_mat2);

        exp_for_mat  = exp(-lda_for * mat);
        exp_for_mat2 = exp(-lda_for2 * mat);
        exp_for_pay  = exp(-lda_for * pay_mat2);
        exp_for_pay2 = exp(-lda_for2 * pay_mat2);

        phi_dom1 = phi_dom2 = phi_dom12 = 0.0;
        expect_dom1 = expect_dom12 = expect_dom2 = expect_dom21 = 0.0;
        adj_pay_dom1 = adj_pay_dom12 = adj_pay_dom2 = adj_pay_dom21 = 0.0;

        phi_for1 = phi_for2 = phi_for12 = 0.0;
        expect_for1 = expect_for12 = expect_for2 = expect_for21 = 0.0;
        adj_pay_for1 = adj_pay_for12 = adj_pay_for2 = adj_pay_for21 = 0.0;
        adj_quanto_for1 = adj_quanto_for2 = 0.0;

        adj_pay_fx1 = adj_pay_fx2 = adj_pay_fx3 = adj_pay_fx4 = adj_pay_fx5 = 0.0;
        adj_pay_fx6 = adj_pay_fx7 = adj_pay_fx8 = adj_pay_fx9 = adj_pay_fx10 = 0.0;

        phi_domfor11 = phi_domfor12 = phi_domfor21 = phi_domfor22 = 0.0;
        covar_dom1_fx1 = covar_dom1_fx2 = covar_dom1_fx3 = covar_dom1_fx4 = covar_dom1_fx5 = 0.0;
        covar_dom2_fx1 = covar_dom2_fx2 = covar_dom2_fx3 = covar_dom2_fx4 = covar_dom2_fx5 = 0.0;
        covar_for1_fx1 = covar_for1_fx2 = covar_for1_fx3 = covar_for1_fx4 = covar_for1_fx5 = 0.0;
        covar_for2_fx1 = covar_for2_fx2 = covar_for2_fx3 = covar_for2_fx4 = covar_for2_fx5 = 0.0;

        for (i = StartIndex; i < EndIndex + 1; i++)
        {
            if (i > StartIndex)
            {
                T1 = sig_dates[i - 1];
            }
            else
            {
                /* First part */
                T1 = start_mat;
            }

            if (i == EndIndex || StartIndex == EndIndex)
            {
                /* Last part */
                T2 = end_mat;
            }
            else
            {
                T2 = sig_dates[i];
            }

            sig_dom = sig_dom2 = sig_dom_for = sig_dom_fx = sig_curve_dom[i];
            sig_dom2 *= sig_dom2;
            sig_for = sig_for2 = sig_for_fx = sig_curve_for[i];
            sig_for2 *= sig_for2;
            sig_dom_for *= sig_for;
            sig_fx = sig_curve_fx[i];
            sig_dom_fx *= sig_fx;
            sig_for_fx *= sig_fx;

            phif_dom1  = Phi_Func(lda_dom, end_mat, T1, T2);
            phif_dom2  = Phi_Func(lda_dom2, end_mat, T1, T2);
            phif_dom11 = Phi_Func(2.0 * lda_dom, end_mat, T1, T2);
            phif_dom22 = Phi_Func(2.0 * lda_dom2, end_mat, T1, T2);
            phif_dom12 = Phi_Func(lda_dom + lda_dom2, end_mat, T1, T2);
            phif_dom21 = phif_dom12;

            phif_for1  = Phi_Func(lda_for, end_mat, T1, T2);
            phif_for2  = Phi_Func(lda_for2, end_mat, T1, T2);
            phif_for11 = Phi_Func(2.0 * lda_for, end_mat, T1, T2);
            phif_for22 = Phi_Func(2.0 * lda_for2, end_mat, T1, T2);
            phif_for12 = Phi_Func(lda_for + lda_for2, end_mat, T1, T2);
            phif_for21 = phif_for12;

            phif_domfor11 = Phi_Func(lda_dom + lda_for, end_mat, T1, T2);
            phif_domfor21 = Phi_Func(lda_dom2 + lda_for, end_mat, T1, T2);
            phif_domfor22 = Phi_Func(lda_dom2 + lda_for2, end_mat, T1, T2);
            phif_domfor12 = Phi_Func(lda_dom + lda_for2, end_mat, T1, T2);

            etha_dom1_Tend = Etha_Func(lda_dom, end_mat, T1, T2);
            etha_dom1_Tpay = Etha_Func(lda_dom, mat_pay, T1, T2);
            etha_dom2_Tend = Etha_Func(lda_dom2, end_mat, T1, T2);
            etha_dom2_Tpay = Etha_Func(lda_dom2, mat_pay, T1, T2);

            psif_domfor11_Tend = Psi2_Func(lda_dom, lda_for, end_mat, end_mat, T1, T2);
            psif_domfor11_Tpay = Psi2_Func(lda_dom, lda_for, mat_pay, end_mat, T1, T2);

            psif_domfor21_Tend = Psi2_Func(lda_dom2, lda_for, end_mat, end_mat, T1, T2);
            psif_domfor21_Tpay = Psi2_Func(lda_dom2, lda_for, mat_pay, end_mat, T1, T2);

            psif_domfor12_Tend = Psi2_Func(lda_dom, lda_for2, end_mat, end_mat, T1, T2);
            psif_domfor12_Tpay = Psi2_Func(lda_dom, lda_for2, mat_pay, end_mat, T1, T2);

            psif_domfor22_Tend = Psi2_Func(lda_dom2, lda_for2, end_mat, end_mat, T1, T2);
            psif_domfor22_Tpay = Psi2_Func(lda_dom2, lda_for2, mat_pay, end_mat, T1, T2);

            psif_domdom11_Tend = Psi2_Func(lda_dom, lda_dom, end_mat, end_mat, T1, T2);
            psif_domdom11_Tpay = Psi2_Func(lda_dom, lda_dom, mat_pay, end_mat, T1, T2);

            psif_domdom21_Tend = Psi2_Func(lda_dom2, lda_dom, end_mat, end_mat, T1, T2);
            psif_domdom21_Tpay = Psi2_Func(lda_dom2, lda_dom, mat_pay, end_mat, T1, T2);

            psif_domdom12_Tend = Psi2_Func(lda_dom, lda_dom2, end_mat, end_mat, T1, T2);
            psif_domdom12_Tpay = Psi2_Func(lda_dom, lda_dom2, mat_pay, end_mat, T1, T2);

            psif_domdom22_Tend = Psi2_Func(lda_dom2, lda_dom2, end_mat, end_mat, T1, T2);
            psif_domdom22_Tpay = Psi2_Func(lda_dom2, lda_dom2, mat_pay, end_mat, T1, T2);

            /* Domestic phi and expectations */
            phi_dom1 += sig_dom2 * phif_dom11;
            phi_dom2 += sig_dom2 * phif_dom22;
            phi_dom12 += sig_dom2 * phif_dom12;

            expect_dom1 += sig_dom2 * (phif_dom1 - phif_dom11);
            expect_dom12 += sig_dom2 * (phif_dom1 - phif_dom12);

            adj_pay_dom1 += sig_dom2 * (phif_dom1 - exp_dom_pay * phif_dom11);
            adj_pay_dom12 += sig_dom2 * (phif_dom1 - exp_dom_pay2 * phif_dom12);

            expect_dom2 += sig_dom2 * (phif_dom2 - phif_dom22);
            expect_dom21 += sig_dom2 * (phif_dom2 - phif_dom21);

            adj_pay_dom2 += sig_dom2 * (phif_dom2 - exp_dom_pay2 * phif_dom22);
            adj_pay_dom21 += sig_dom2 * (phif_dom2 - exp_dom_pay * phif_dom21);

            /* Foreign phi and expectations */
            phi_for1 += sig_for2 * phif_for11;
            phi_for2 += sig_for2 * phif_for22;
            phi_for12 += sig_for2 * phif_for12;

            expect_for1 += sig_for2 * (phif_for1 - phif_for11);
            expect_for12 += sig_for2 * (phif_for1 - phif_for12);

            adj_pay_for1 += sig_dom_for * (phif_for1 - exp_dom_pay * phif_domfor11);
            adj_pay_for12 += sig_dom_for * (phif_for1 - exp_dom_pay2 * phif_domfor21);

            adj_quanto_for1 += sig_for_fx * phif_for1;

            expect_for2 += sig_for2 * (phif_for2 - phif_for22);
            expect_for21 += sig_for2 * (phif_for2 - phif_for12);

            adj_pay_for2 += sig_dom_for * (phif_for2 - exp_dom_pay * phif_domfor12);
            adj_pay_for21 += sig_dom_for * (phif_for2 - exp_dom_pay2 * phif_domfor22);

            adj_quanto_for2 += sig_for_fx * phif_for2;

            /* Fx pay adjustemt*/
            adj_pay_fx1 += sig_dom_fx * (etha_dom1_Tend - etha_dom1_Tpay);
            adj_pay_fx2 += sig_dom_fx * (etha_dom2_Tend - etha_dom2_Tpay);

            adj_pay_fx3 += sig_dom_for * (psif_domfor11_Tend - psif_domfor11_Tpay);
            adj_pay_fx4 += sig_dom_for * (psif_domfor21_Tend - psif_domfor21_Tpay);

            adj_pay_fx5 += sig_dom_for * (psif_domfor12_Tend - psif_domfor12_Tpay);
            adj_pay_fx6 += sig_dom_for * (psif_domfor22_Tend - psif_domfor22_Tpay);

            adj_pay_fx7 += sig_dom2 * (psif_domdom11_Tend - psif_domdom11_Tpay);
            adj_pay_fx8 += sig_dom2 * (psif_domdom21_Tend - psif_domdom21_Tpay);

            adj_pay_fx9 += sig_dom2 * (psif_domdom12_Tend - psif_domdom12_Tpay);
            adj_pay_fx10 += sig_dom2 * (psif_domdom22_Tend - psif_domdom22_Tpay);

            /* Covariances */
            phi_domfor11 += sig_dom_for * phif_domfor11;
            phi_domfor12 += sig_dom_for * phif_domfor12;
            phi_domfor21 += sig_dom_for * phif_domfor21;
            phi_domfor22 += sig_dom_for * phif_domfor22;

            /* dom1 / fx */
            covar_dom1_fx1 += sig_dom_fx * phif_dom1;
            covar_dom1_fx2 += sig_dom_for * (phif_dom1 - phif_domfor11);
            covar_dom1_fx3 += sig_dom_for * (phif_dom1 - phif_domfor12);
            covar_dom1_fx4 += sig_dom2 * (phif_dom1 - phif_dom11);
            covar_dom1_fx5 += sig_dom2 * (phif_dom1 - phif_dom12);

            /* dom2 / fx */
            covar_dom2_fx1 += sig_dom_fx * phif_dom2;
            covar_dom2_fx2 += sig_dom_for * (phif_dom2 - phif_domfor21);
            covar_dom2_fx3 += sig_dom_for * (phif_dom2 - phif_domfor22);
            covar_dom2_fx4 += sig_dom2 * (phif_dom2 - phif_dom21);
            covar_dom2_fx5 += sig_dom2 * (phif_dom2 - phif_dom22);

            /* for1 / fx */
            covar_for1_fx1 += sig_for_fx * phif_for1;
            covar_for1_fx2 += sig_for2 * (phif_for1 - phif_for11);
            covar_for1_fx3 += sig_for2 * (phif_for1 - phif_for12);
            covar_for1_fx4 += sig_dom_for * (phif_for1 - phif_domfor11);
            covar_for1_fx5 += sig_dom_for * (phif_for1 - phif_domfor12);

            /* for2 / fx */
            covar_for2_fx1 += sig_for_fx * phif_for2;
            covar_for2_fx2 += sig_for2 * (phif_for2 - phif_for21);
            covar_for2_fx3 += sig_for2 * (phif_for2 - phif_for22);
            covar_for2_fx4 += sig_dom_for * (phif_for2 - phif_domfor21);
            covar_for2_fx5 += sig_dom_for * (phif_for2 - phif_domfor22);
        }

        dom_exp1[k + 1] = exp_dom_mat;
        dom_exp2[k + 1] = exp_dom_mat2;

        dom_phi1[k + 1] = dom_phi1[k] * exp_dom_mat * exp_dom_mat + phi_dom1;
        dom_phi2[k + 1] = dom_phi2[k] * exp_dom_mat2 * exp_dom_mat2 + alpha_dom2 * phi_dom2;
        dom_phi12[k + 1] =
            dom_phi12[k] * exp_dom_mat * exp_dom_mat2 + correl[0][1] * alpha_dom * phi_dom12;

        dom_fwd1[k + 1] =
            dom_phi1[k] * exp_dom_mat * Phi_Func(-lda_dom, start_mat, start_mat, end_mat) +
            expect_dom1 / lda_dom +
            dom_phi12[k] * exp_dom_mat * Phi_Func(-lda_dom2, start_mat, start_mat, end_mat) +
            correl[0][1] * alpha_dom * expect_dom12 / lda_dom2 - adj_pay_dom1 / lda_dom -
            correl[0][1] * alpha_dom * adj_pay_dom12 / lda_dom2;

        dom_fwd2[k + 1] =
            dom_phi2[k] * exp_dom_mat2 * Phi_Func(-lda_dom2, start_mat, start_mat, end_mat) +
            alpha_dom2 * expect_dom2 / lda_dom2 +
            dom_phi12[k] * exp_dom_mat2 * Phi_Func(-lda_dom, start_mat, start_mat, end_mat) +
            correl[0][1] * alpha_dom * expect_dom21 / lda_dom -
            adj_pay_dom2 * alpha_dom2 / lda_dom2 -
            correl[0][1] * alpha_dom * adj_pay_dom21 / lda_dom;

        for_exp1[k + 1] = exp_for_mat;
        for_exp2[k + 1] = exp_for_mat2;

        for_phi1[k + 1] = for_phi1[k] * exp_for_mat * exp_for_mat + phi_for1;
        for_phi2[k + 1] = for_phi2[k] * exp_for_mat2 * exp_for_mat2 + alpha_for2 * phi_for2;
        for_phi12[k + 1] =
            for_phi12[k] * exp_for_mat * exp_for_mat2 + correl[2][3] * alpha_for * phi_for12;

        for_fwd1[k + 1] =
            for_phi1[k] * exp_for_mat * Phi_Func(-lda_for, start_mat, start_mat, end_mat) +
            expect_for1 / lda_for +
            for_phi12[k] * exp_for_mat * Phi_Func(-lda_for2, start_mat, start_mat, end_mat) +
            correl[2][3] * alpha_for * expect_for12 / lda_for2 -
            correl[0][2] * adj_pay_for1 / lda_dom -
            correl[1][2] * alpha_dom * adj_pay_for12 / lda_dom2 - correl[2][4] * adj_quanto_for1;

        for_fwd2[k + 1] =
            for_phi2[k] * exp_for_mat2 * Phi_Func(-lda_for2, start_mat, start_mat, end_mat) +
            alpha_for2 * expect_for2 / lda_for2 +
            for_phi12[k] * exp_for_mat2 * Phi_Func(-lda_for, start_mat, start_mat, end_mat) +
            correl[2][3] * alpha_for * expect_for21 / lda_for -
            correl[0][3] * alpha_for * adj_pay_for2 / lda_dom -
            correl[1][3] * alpha_dom * alpha_for * adj_pay_for21 / lda_dom2 -
            correl[3][4] * alpha_for * adj_quanto_for2;

        fx_fwd[k + 1] =
            expect_fx + correl[0][4] * adj_pay_fx1 + correl[1][4] * alpha_dom * adj_pay_fx2 -
            correl[0][2] * adj_pay_fx3 - correl[1][2] * alpha_dom * adj_pay_fx4 -
            correl[0][3] * alpha_for * adj_pay_fx5 - correl[1][3] * alpha_domfor * adj_pay_fx6 +
            adj_pay_fx7 + correl[0][1] * alpha_dom * adj_pay_fx8 +
            correl[0][1] * alpha_dom * adj_pay_fx9 + alpha_dom2 * adj_pay_fx10;

        cov = covar[k + 1];

        cov[0][0] = phi_dom1;
        cov[0][1] = cov[1][0] = correl[0][1] * alpha_dom * phi_dom12;
        cov[0][2] = cov[2][0] = correl[0][2] * phi_domfor11;
        cov[0][3] = cov[3][0] = correl[0][3] * alpha_for * phi_domfor12;
        cov[0][4]             = cov[4][0] =
            correl[0][4] * covar_dom1_fx1 - correl[0][2] * covar_dom1_fx2 / lda_for -
            correl[0][3] * alpha_for * covar_dom1_fx3 / lda_for2 + covar_dom1_fx4 / lda_dom +
            correl[0][1] * alpha_dom * covar_dom1_fx5 / lda_dom2;

        cov[1][1] = alpha_dom2 * phi_dom2;
        cov[1][2] = cov[2][1] = correl[1][2] * alpha_dom * phi_domfor21;
        cov[1][3] = cov[3][1] = correl[1][3] * alpha_domfor * phi_domfor22;
        cov[1][4] = cov[4][1] = correl[1][4] * alpha_dom * covar_dom2_fx1 -
                                correl[1][2] * alpha_dom * covar_dom2_fx2 / lda_for -
                                correl[1][3] * alpha_domfor * covar_dom2_fx3 / lda_for2 +
                                correl[0][1] * alpha_dom * covar_dom2_fx4 / lda_dom +
                                alpha_dom2 * covar_dom2_fx5 / lda_dom2;

        cov[2][2] = phi_for1;
        cov[2][3] = cov[3][2] = correl[2][3] * alpha_for * phi_for12;
        cov[2][4] = cov[4][2] = correl[2][4] * covar_for1_fx1 - covar_for1_fx2 / lda_for -
                                correl[2][3] * alpha_for * covar_for1_fx3 / lda_for2 +
                                correl[0][2] * covar_for1_fx4 / lda_dom +
                                correl[1][2] * alpha_dom * covar_for1_fx5 / lda_dom2;

        cov[3][3] = alpha_for2 * phi_for2;
        cov[3][4] = cov[4][3] = correl[3][4] * alpha_for * covar_for2_fx1 -
                                correl[2][3] * alpha_for * covar_for2_fx2 / lda_for -
                                alpha_for2 * covar_for2_fx3 / lda_for2 +
                                correl[0][3] * alpha_for * covar_for2_fx4 / lda_dom +
                                correl[1][3] * alpha_domfor * covar_for2_fx5 / lda_dom2;

        cov[4][4] = fx_vol * fx_vol * mat;

        start_date = end_date;
        start_mat  = end_mat;
        StartIndex = EndIndex;
    }

    dom_ifr[nb_dates - 1]      = swp_f_zr(date[k], date[k] + 1, dom_yc);
    for_ifr[nb_dates - 1]      = swp_f_zr(date[k], date[k] + 1, for_yc);
    dom_gam1_pay[nb_dates - 1] = 0.0;
    dom_gam2_pay[nb_dates - 1] = 0.0;
    dom_bond_pay[nb_dates - 1] = 1.0;

    return err;
}
