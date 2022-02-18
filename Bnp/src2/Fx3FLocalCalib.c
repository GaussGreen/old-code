/* ==========================================================================
   FILE_NAME:	Fx3FCalib.c

   PURPOSE:		Modelling of the spot FX vol by taking in consideration a LGM 1 factor
                                on the domestic and foreign market and a lognormal model on the Spot
   Fx

   DATE:		05/25/00

   AUTHOR:		L.C.
   ========================================================================== */

#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

#define MAX_TIME 0.02
#define SHIFT_VOL 0.005
#define MAX_ERROR 0.0001

Err FxLocal_log_approx(
    long    today,
    double* maturity,
    long    nb_mat,
    double* sig_dom,
    double* sig_for,
    double* mat_fx,
    long    nb_mat_fx,
    double* sig_fx,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double* params),
    double* params,
    double  lam_dom,
    double  lam_for,
    double  corr_dom_for,
    double  corr_dom_fx,
    double  corr_for_fx,
    double  spot,
    char*   dom_yc,
    char*   for_yc,
    double* time_fx,
    long    nb_time_fx,
    double* fx_fwd,
    double* fx_vol,
    double* fx_fwd0,
    double  max_time)
{
    long   i, j, k, k2;
    long   index, nb_time, last_index;
    double t, prev_t, dt, T, last_mat;
    double logSpot;
    double volfx;
    double F, sigdom, sigfor, sigfx, x;

    double *time2 = NULL, *S = NULL, *Sfwd = NULL;

    Err err = NULL;

    logSpot = log(spot);

    /* discretise time */
    last_mat = time_fx[nb_time_fx - 1];
    nb_time  = (long)(last_mat / max_time - 1.0E-08) + 1;
    if (nb_time < 2)
    {
        nb_time = 2;
    }

    time2 = calloc(nb_time, sizeof(double));
    if (!time2)
    {
        err = "Memory allocation failure (1) in fwd_approx";
        goto FREE_RETURN;
    }

    time2[0] = 0.0;

    for (i = 1; i < nb_time - 1; i++)
    {
        time2[i] = time2[i - 1] + max_time;
    }
    time2[i] = last_mat;

    num_f_concat_vector(&nb_time, &time2, nb_mat, maturity);
    num_f_concat_vector(&nb_time, &time2, nb_mat_fx, mat_fx);
    num_f_concat_vector(&nb_time, &time2, nb_time_fx, time_fx);
    num_f_sort_vector(nb_time, time2);
    num_f_unique_vector(&nb_time, time2);

    /* find the index of the last maturity */
    last_index = Get_Index(last_mat, time2, nb_time);

    S    = (double*)calloc(last_index + 1, sizeof(double));
    Sfwd = (double*)calloc(last_index + 1, sizeof(double));

    if (!S || !Sfwd)
    {
        err = "Memory allocation failure (2) in fwd_approx";
        goto FREE_RETURN;
    }

    /* diffuse */

    fx_fwd[0]  = logSpot;
    fx_fwd0[0] = spot;

    S[0]    = spot;
    Sfwd[0] = spot;

    for (i = 0; i < last_index; i++)
    {
        t = 0;
        T = time2[i + 1];
        /* initialisation at the fwd */
        F = logSpot + (swp_f_zr(today, (long)(today + T * 365.0 + 1.0E-08), dom_yc) -
                       swp_f_zr(today, (long)(today + T * 365.0 + 1.0E-08), for_yc)) *
                          T;

        Sfwd[i + 1] = exp(F);
        k           = 0;
        k2          = 0;
        sigdom      = sig_dom[0];
        sigfor      = sig_for[0];
        sigfx       = sig_fx[0];

        for (j = 0; j < i + 1; j++)
        {
            prev_t = t;
            t      = time2[j + 1];

            if (t > maturity[k])
            {
                if (k < nb_mat - 1)
                {
                    k += 1;
                }
                sigdom = sig_dom[k];
                sigfor = sig_for[k];
            }
            if (t > mat_fx[k2])
            {
                if (k2 < nb_mat_fx - 1)
                {
                    k2 += 1;
                }
                sigfx = sig_fx[k2];
            }

            dt = t - prev_t;
            x  = sigdom * (1.0 - exp(-lam_dom * (T - prev_t))) / lam_dom;

            F += x * dt *
                 (corr_dom_fx * sigfx * vol_ln_func(prev_t, S[j], Sfwd[j], params) -
                  corr_dom_for * sigfor * (1.0 - exp(-lam_for * (T - prev_t))) / lam_for + x);
        }
        S[i + 1] = exp(F);
    }

    j      = 0;
    prev_t = 0;
    /* fill the out structure */
    for (i = 1; i < nb_time_fx; i++)
    {
        t              = time_fx[i];
        index          = Get_Index(t, mat_fx, nb_mat_fx);
        volfx          = sig_fx[index];
        index          = Get_Index(prev_t, time2, last_index + 1);
        fx_fwd[i - 1]  = S[index];
        fx_fwd0[i - 1] = Sfwd[index];
        fx_vol[i]      = volfx * vol_ln_func(t, fx_fwd[i - 1], fx_fwd0[i - 1], params);
        prev_t         = t;
    }
    index          = Get_Index(t, time2, last_index + 1);
    fx_fwd[i - 1]  = S[index];
    fx_fwd0[i - 1] = Sfwd[index];

FREE_RETURN:

    if (time2)
    {
        free(time2);
    }

    if (S)
    {
        free(S);
    }

    if (Sfwd)
    {
        free(Sfwd);
    }

    return err;
}

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx3DLocaltsImpliedVol(
    long    today,
    double  opt_maturity,
    double  start_date,
    double  end_date,
    double* maturity,
    long    nbMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* mat_fx,
    long    nb_mat_fx,
    double* sig_curve_fx,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double* params),
    double* params,
    double  spot_fx,
    double  corr_dom_for,
    double  corr_dom_fx,
    double  corr_for_fx,
    char*   dom_yc,
    char*   for_yc,
    double* fx_vol,
    double  disc_dt,
    double  fx_dt)

{
    double *fx_vol_curve = NULL, *fx_fwd = NULL, *fx_time = NULL;

    long i;
    long nb_fx_time, nb_fx_time2;

    Err err = NULL;

    /* simple translation to make everything begin at 0 */

    nb_fx_time = (long)(end_date / fx_dt - 1.0E-08) + 1;
    if (nb_fx_time < 2)
    {
        nb_fx_time = 2;
    }

    fx_time = calloc(nb_fx_time, sizeof(double));

    if (!fx_time)
    {
        err = "Memory allocation error in Fx3DLocaltsImpliedVol (1)";
        goto FREE_RETURN;
    }

    fx_time[0] = 0.0;

    /* fill the time where we want the volatility structure */
    for (i = 1; i < nb_fx_time - 1; i++)
    {
        fx_time[i] = fx_time[i - 1] + fx_dt;
    }

    fx_time[i] = end_date;

    num_f_concat_vector(&nb_fx_time, &fx_time, nbMat, maturity);
    num_f_concat_vector(&nb_fx_time, &fx_time, nb_mat_fx, mat_fx);
    num_f_sort_vector(nb_fx_time, fx_time);
    num_f_unique_vector(&nb_fx_time, fx_time);

    nb_fx_time2 = Get_Index(end_date, fx_time, nb_fx_time) + 1;

    fx_vol_curve = calloc(nb_fx_time2, sizeof(double));
    fx_fwd       = calloc(nb_fx_time2, sizeof(double));

    if (!fx_vol || !fx_fwd)
    {
        err = "Memory allocation error in Fx3DLocaltsImpliedVol (2)";
        goto FREE_RETURN;
    }

    err = FxLocal_log_approx(
        today,
        maturity,
        nbMat,
        sig_curve_dom,
        sig_curve_for,
        mat_fx,
        nb_mat_fx,
        sig_curve_fx,
        vol_ln_func,
        params,
        lda_dom,
        lda_for,
        corr_dom_for,
        corr_dom_fx,
        corr_for_fx,
        spot_fx,
        dom_yc,
        for_yc,
        fx_time,
        nb_fx_time2,
        fx_fwd,
        fx_vol_curve,
        fx_fwd,
        disc_dt);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = Fx3DtsImpliedVol(
        opt_maturity,
        start_date,
        end_date,
        maturity,
        nbMat,
        sig_curve_dom,
        lda_dom,
        sig_curve_for,
        lda_for,
        fx_time,
        fx_vol_curve,
        nb_fx_time2,
        corr_dom_for,
        corr_dom_fx,
        corr_for_fx,
        fx_vol);

FREE_RETURN:

    if (fx_vol_curve)
        free(fx_vol_curve);
    if (fx_time)
        free(fx_time);
    if (fx_fwd)
        free(fx_fwd);

    return err;
}

/*	Implied vol direct from underlying */
Err Fx3DLocalImpliedVol(
    char* fx_underlying,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double* params),
    double* params,
    double  val_time,
    double  start_time,
    double  end_time,
    double  disc_dt,
    double  fx_dt,
    double* vol)
{
    long    sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx, nb_merge_dates;
    long    i;
    double *sigma_date_dom = NULL, *sigma_dom = NULL;
    double *tau_date_dom = NULL, *tau_dom = NULL;
    double *sigma_date_for = NULL, *sigma_for = NULL;
    double *tau_date_for = NULL, *tau_for = NULL;
    double *sigma_date_fx = NULL, *sigma_fx = NULL;
    double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL, *merge_dates = NULL;

    double correl_dom_for, correl_dom_fx, correl_for_fx;
    double lda_dom, lda_for;

    SrtUndPtr fx_und, dom_und, for_und;
    char *    domname, *forname;

    long   today;
    double spot_fx;
    char * dom_yc, *for_yc;

    Err err = NULL;

    err = Get_FX_StochRate_TermStructures(
        fx_underlying,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_date_dom,
        &tau_dom,
        &tau_n_dom,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_date_for,
        &tau_for,
        &tau_n_for,
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

    err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
    if (err)
    {
        goto FREE_RETURN;
    }

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
        err = "Memory allocation error (1) in Fx3DLocalImpliedVol";
        goto FREE_RETURN;
    }

    for (i = nb_merge_dates - 1; i >= 0; i--)
    {
        sig_dom[i] = sigma_dom[Get_Index(merge_dates[i], sigma_date_dom, sigma_n_dom)];
        sig_for[i] = sigma_for[Get_Index(merge_dates[i], sigma_date_for, sigma_n_for)];
        sig_fx[i]  = sigma_fx[Get_Index(merge_dates[i], sigma_date_fx, sigma_n_fx)];
    }

    /* get today, spot, dom_yc, for_yc */
    fx_und  = lookup_und(fx_underlying);
    domname = get_domname_from_fxund(fx_und);
    forname = get_forname_from_fxund(fx_und);
    dom_und = lookup_und(domname);
    for_und = lookup_und(forname);

    err = get_underlying_discname(dom_und, &dom_yc);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = get_underlying_discname(for_und, &for_yc);
    if (err)
    {
        goto FREE_RETURN;
    }

    today = get_today_from_underlying(dom_und);

    err = Fx3DtsFwdFx(dom_und, for_und, fx_und, today, &spot_fx);

    err = Fx3DLocaltsImpliedVol(
        today,
        val_time,
        start_time,
        end_time,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        lda_dom,
        sig_for,
        lda_for,
        merge_dates,
        nb_merge_dates,
        sig_fx,
        vol_ln_func,
        params,
        spot_fx,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        dom_yc,
        for_yc,
        vol,
        disc_dt,
        fx_dt);

FREE_RETURN:

    if (sigma_date_dom)
    {
        free(sigma_date_dom);
    }

    if (sigma_dom)
    {
        free(sigma_dom);
    }

    if (tau_date_dom)
    {
        free(tau_date_dom);
    }

    if (tau_dom)
    {
        free(tau_dom);
    }

    if (sigma_date_for)
    {
        free(sigma_date_for);
    }

    if (sigma_for)
    {
        free(sigma_for);
    }

    if (tau_date_for)
    {
        free(tau_date_for);
    }

    if (tau_for)
    {
        free(tau_for);
    }

    if (sigma_date_fx)
    {
        free(sigma_date_fx);
    }

    if (sigma_fx)
    {
        free(sigma_fx);
    }

    if (sig_dom)
    {
        free(sig_dom);
    }

    if (sig_for)
    {
        free(sig_for);
    }

    if (sig_fx)
    {
        free(sig_fx);
    }

    if (merge_dates)
    {
        free(merge_dates);
    }

    return err;
}

Err Fx3DLocaltsCalibration(
    long    today,
    double* exercise_opt,
    double* maturity_opt,
    double* vol_opt,
    long    nbrOpt,
    double* maturity,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double* params),
    double*  params,
    double   spot_fx,
    double   correl_dom_for,
    double   correl_dom_fx,
    double   correl_for_fx,
    char*    dom_yc,
    char*    for_yc,
    double** fx_vol_curve,
    double   disc_dt,
    double   fx_dt,
    long     nbIterMax)

{
    long    i, nIter, spot_date;
    double  T, prev_t, mat, Fwd;
    double* fx_vol_log = NULL;

    double vol_imp1, vol_imp2;
    double vol_tgt;
    double error;

    Err err = NULL;

    (*fx_vol_curve) = NULL;
    (*fx_vol_curve) = (double*)calloc(nbrOpt, sizeof(double));

    if (!(*fx_vol_curve))
    {
        err = "Memory allocation error in Fx3DLocaltsCalibration";
        goto FREE_RETURN;
    }

    /* First we calibrate a lognormal model */
    err = Fx3DtsCalibration(
        exercise_opt,
        maturity_opt,
        vol_opt,
        nbrOpt,
        maturity,
        nbrMat,
        sig_curve_dom,
        lda_dom,
        sig_curve_for,
        lda_for,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        &fx_vol_log);
    if (err)
    {
        goto FREE_RETURN;
    }

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    spot_fx   = spot_fx / swp_f_df(today, spot_date, for_yc) * swp_f_df(today, spot_date, dom_yc);

    prev_t = 0;
    for (i = 0; i < nbrOpt; i++)
    {
        /* Initialisation: */
        mat = maturity_opt[i];
        T   = exercise_opt[i];
        Fwd = spot_fx * swp_f_df(today, today + prev_t * 365.0, for_yc) /
              swp_f_df(today, today + prev_t * 365.0, dom_yc);
        vol_tgt = vol_opt[i];

        (*fx_vol_curve)[i] = fx_vol_log[i] / vol_ln_func(prev_t, Fwd, Fwd, params);

        err = Fx3DLocaltsImpliedVol(
            today,
            mat,
            0.0,
            T,
            maturity,
            nbrMat,
            sig_curve_dom,
            lda_dom,
            sig_curve_for,
            lda_for,
            exercise_opt,
            nbrOpt,
            (*fx_vol_curve),
            vol_ln_func,
            params,
            spot_fx,
            correl_dom_for,
            correl_dom_fx,
            correl_for_fx,
            dom_yc,
            for_yc,
            &vol_imp1,
            disc_dt,
            fx_dt);
        if (err)
        {
            goto FREE_RETURN;
        }

        error = fabs(vol_imp1 - vol_tgt);

        nIter = 0;

        while (nIter <= nbIterMax && error > MAX_ERROR)
        {
            (*fx_vol_curve)[i] += SHIFT_VOL;
            err = Fx3DLocaltsImpliedVol(
                today,
                mat,
                0.0,
                T,
                maturity,
                nbrMat,
                sig_curve_dom,
                lda_dom,
                sig_curve_for,
                lda_for,
                exercise_opt,
                nbrOpt,
                (*fx_vol_curve),
                vol_ln_func,
                params,
                spot_fx,
                correl_dom_for,
                correl_dom_fx,
                correl_for_fx,
                dom_yc,
                for_yc,
                &vol_imp2,
                disc_dt,
                fx_dt);

            if (err)
            {
                goto FREE_RETURN;
            }

            (*fx_vol_curve)[i] -= SHIFT_VOL;
            (*fx_vol_curve)[i] += (vol_tgt - vol_imp1) * SHIFT_VOL / (vol_imp2 - vol_imp1);

            if ((*fx_vol_curve)[i] < 0)
            {
                err = serror("Cannot calibrate option %d", i);
                goto FREE_RETURN;
            }

            err = Fx3DLocaltsImpliedVol(
                today,
                mat,
                0.0,
                T,
                maturity,
                nbrMat,
                sig_curve_dom,
                lda_dom,
                sig_curve_for,
                lda_for,
                exercise_opt,
                nbrOpt,
                (*fx_vol_curve),
                vol_ln_func,
                params,
                spot_fx,
                correl_dom_for,
                correl_dom_fx,
                correl_for_fx,
                dom_yc,
                for_yc,
                &vol_imp1,
                disc_dt,
                fx_dt);

            if (err)
            {
                goto FREE_RETURN;
            }

            error = fabs(vol_imp1 - vol_tgt);

            nIter += 1;
        }

        prev_t = T;
    }

FREE_RETURN:

    if (fx_vol_log)
        free(fx_vol_log);

    return err;
}

Err Fx3DLocalCalibration(
    char*  dom_underlying,
    char*  for_underlying,
    double spot_fx,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double* params),
    double*  params,
    double   correl_dom_for,
    double   correl_dom_fx,
    double   correl_for_fx,
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbropt,
    double** fx_vol_curve,
    double   disc_dt,
    double   fx_dt,
    long     nbIterMax)

{
    long    sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
    long    nb_merge_dates;
    double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL, *tau_dom = NULL,
           *sigma_date_for = NULL, *sigma_for = NULL, *tau_date_for = NULL, *tau_for = NULL,
           *merge_dates = NULL, *sig_dom = NULL, *sig_for = NULL;

    double lda_dom, lda_for;

    SrtUndPtr dom_und, for_und;

    long  today;
    char *dom_yc, *for_yc;

    Err err = NULL;

    err = Get_LGM_TermStructure(
        dom_underlying,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_date_dom,
        &tau_dom,
        &tau_n_dom);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = Get_LGM_TermStructure(
        for_underlying,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_date_for,
        &tau_for,
        &tau_n_for);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
    if (err)
    {
        goto FREE_RETURN;
    }

    dom_und = lookup_und(dom_underlying);
    for_und = lookup_und(for_underlying);

    err = get_underlying_discname(dom_und, &dom_yc);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = get_underlying_discname(for_und, &for_yc);
    if (err)
    {
        goto FREE_RETURN;
    }

    today = get_today_from_underlying(dom_und);

    err = merge_rates_ts(
        sigma_date_dom,
        sigma_dom,
        sigma_n_dom,
        sigma_date_for,
        sigma_for,
        sigma_n_for,
        &merge_dates,
        &sig_dom,
        &sig_for,
        &nb_merge_dates);

    if (err)
    {
        goto FREE_RETURN;
    }

    err = Fx3DLocaltsCalibration(
        today,
        exercise_opt,
        maturity_opt,
        vol_opt,
        nbropt,
        merge_dates,
        nb_merge_dates,
        sigma_dom,
        lda_dom,
        sigma_for,
        lda_for,
        vol_ln_func,
        params,
        spot_fx,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        dom_yc,
        for_yc,
        fx_vol_curve,
        disc_dt,
        fx_dt,
        nbIterMax);

FREE_RETURN:

    if (sigma_date_dom)
    {
        free(sigma_date_dom);
    }

    if (sigma_dom)
    {
        free(sigma_dom);
    }

    if (tau_date_dom)
    {
        free(tau_date_dom);
    }

    if (tau_dom)
    {
        free(tau_dom);
    }

    if (sigma_date_for)
    {
        free(sigma_date_for);
    }

    if (sigma_for)
    {
        free(sigma_for);
    }

    if (tau_date_for)
    {
        free(tau_date_for);
    }

    if (tau_for)
    {
        free(tau_for);
    }

    if (sig_dom)
    {
        free(sig_dom);
    }

    if (sig_for)
    {
        free(sig_for);
    }

    if (merge_dates)
    {
        free(merge_dates);
    }

    return err;
}
