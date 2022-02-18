
#include "LGM1FQuantoFwd.h"

#include "FX3FCalib.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

static double Ji0_func(double lam, double t1, double t2)
{
    return (exp(lam * t2) - exp(lam * t1)) / lam;
}

static double LGM1FQuantoFwdAdjustment(
    double  mattime,
    double  starttime,
    double  endtime,
    double* sigtime,
    long    nbsig,
    double* domsig,
    double* forsig,
    double* fxsig,
    double  domlam,
    double  forlam,
    double  domforrho,
    double  forfxrho)
{
    long   i;
    double I1, I2, I3;
    double I2_a, I2_b;
    double I3_a, I3_b;
    double QAdj;
    double coverage;

    coverage = endtime - starttime;

    i    = 0;
    I1   = 0.0;
    I2_a = 0.0;
    I2_b = 0.0;
    I3_a = 0.0;
    I3_b = 0.0;

    if (mattime <= sigtime[0])
    {
        I1 += forfxrho * fxsig[0] * forsig[0] * Ji0_func(forlam, 0, mattime);

        I2_a += forsig[0] * forsig[0] * Ji0_func(forlam, 0, mattime);
        I2_b += forsig[0] * forsig[0] * Ji0_func(2 * forlam, 0, mattime);

        I3_a += domforrho * forsig[0] * domsig[0] * Ji0_func(forlam, 0, mattime);
        I3_b += domforrho * forsig[0] * domsig[0] * Ji0_func(domlam + forlam, 0, mattime);

        I1 = -I1 * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / forlam;

        I2_a = I2_a * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / (forlam * forlam);
        I2_b = I2_b * exp(-forlam * starttime) * exp(-forlam * endtime) *
               (1 - exp(-forlam * coverage)) / (forlam * forlam);
        I2 = I2_a - I2_b;

        I3_a = I3_a * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / (forlam * domlam);
        I3_b = I3_b * exp(-forlam * starttime) * exp(-domlam * endtime) *
               (1 - exp(-forlam * coverage)) / (forlam * domlam);
        I3 = -I3_a + I3_b;

        QAdj = exp(I1 + I2 + I3);
    }
    else
    {
        I1 += forfxrho * fxsig[0] * forsig[0] * Ji0_func(forlam, 0, sigtime[0]);

        I2_a += forsig[0] * forsig[0] * Ji0_func(forlam, 0, sigtime[0]);
        I2_b += forsig[0] * forsig[0] * Ji0_func(2 * forlam, 0, sigtime[0]);

        I3_a += domforrho * forsig[0] * domsig[0] * Ji0_func(forlam, 0, sigtime[0]);
        I3_b += domforrho * forsig[0] * domsig[0] * Ji0_func(domlam + forlam, 0, sigtime[0]);

        while (sigtime[i + 1] < mattime && i < nbsig - 1)
        {
            I1 += forfxrho * fxsig[i + 1] * forsig[i + 1] *
                  Ji0_func(forlam, sigtime[i], sigtime[i + 1]);

            I2_a += forsig[i + 1] * forsig[i + 1] * Ji0_func(forlam, sigtime[i], sigtime[i + 1]);
            I2_b +=
                forsig[i + 1] * forsig[i + 1] * Ji0_func(2 * forlam, sigtime[i], sigtime[i + 1]);

            I3_a += domforrho * forsig[i + 1] * domsig[i + 1] *
                    Ji0_func(forlam, sigtime[i], sigtime[i + 1]);
            I3_b += domforrho * forsig[i + 1] * domsig[i + 1] *
                    Ji0_func(domlam + forlam, sigtime[i], sigtime[i + 1]);

            i++;
        }

        if (i == 0)
        {
            I1 += forfxrho * fxsig[1] * forsig[1] * Ji0_func(forlam, sigtime[0], mattime);

            I2_a += forsig[1] * forsig[1] * Ji0_func(forlam, sigtime[0], mattime);
            I2_b += forsig[1] * forsig[1] * Ji0_func(2 * forlam, sigtime[0], mattime);

            I3_a += domforrho * forsig[1] * domsig[1] * Ji0_func(forlam, sigtime[0], mattime);
            I3_b +=
                domforrho * forsig[1] * domsig[1] * Ji0_func(domlam + forlam, sigtime[0], mattime);

            I1 = -I1 * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / forlam;

            I2_a =
                I2_a * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / (forlam * forlam);
            I2_b = I2_b * exp(-forlam * starttime) * exp(-forlam * (endtime)) *
                   (1 - exp(-forlam * coverage)) / (forlam * forlam);
            I2 = I2_a - I2_b;

            I3_a =
                I3_a * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / (forlam * domlam);
            I3_b = I3_b * exp(-forlam * starttime) * exp(-domlam * (endtime)) *
                   (1 - exp(-forlam * coverage)) / (forlam * domlam);
            I3 = -I3_a + I3_b;

            QAdj = exp(I1 + I2 + I3);
        }
        else
        {
            I1 += forfxrho * fxsig[i + 1] * forsig[i + 1] * Ji0_func(forlam, sigtime[i], mattime);

            I2_a += forsig[i + 1] * forsig[i + 1] * Ji0_func(forlam, sigtime[i], mattime);
            I2_b += forsig[i + 1] * forsig[i + 1] * Ji0_func(2 * forlam, sigtime[i], mattime);

            I3_a +=
                domforrho * forsig[i + 1] * domsig[i + 1] * Ji0_func(forlam, sigtime[i], mattime);
            I3_b += domforrho * forsig[i + 1] * domsig[i + 1] *
                    Ji0_func(domlam + forlam, sigtime[i], mattime);

            I1 = -I1 * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / forlam;

            I2_a =
                I2_a * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / (forlam * forlam);
            I2_b = I2_b * exp(-forlam * starttime) * exp(-forlam * (endtime)) *
                   (1 - exp(-forlam * coverage)) / (forlam * forlam);
            I2 = I2_a - I2_b;

            I3_a =
                I3_a * exp(-forlam * starttime) * (1 - exp(-forlam * coverage)) / (forlam * domlam);
            I3_b = I3_b * exp(-forlam * starttime) * exp(-domlam * (endtime)) *
                   (1 - exp(-forlam * coverage)) / (forlam * domlam);
            I3 = -I3_a + I3_b;

            QAdj = exp(I1 + I2 + I3);
        }
    }

    return QAdj;
}

Err LGM1FQuantoFwd(
    char*    und3dfx,
    long     start,
    long     end,
    char*    freq,
    char*    basis,
    char*    refrate,
    long*    nbFwd,
    double** QuantoFwd)
{
    double QAdj;

    long today, spot_date;

    int           i;
    SrtUndPtr     fx_und, dom_und, for_und;
    TermStruct *  fx_ts, *dom_ts, *for_ts;
    SrtCorrLstPtr sCorrlist;
    char *        domname, *forname;
    double        dom_lam, for_lam;
    char *        dom_yc, *for_yc;

    double* domsigtime = NULL;
    double* domsig     = NULL;
    long    domsig_n;
    double* domtau     = NULL;
    double* domtautime = NULL;
    long    domtau_n;

    double* forsigtime = NULL;
    double* forsig     = NULL;
    long    forsig_n;
    double* fortau     = NULL;
    double* fortautime = NULL;
    long    fortau_n;

    double* fxsigtime = NULL;
    double* fxsig     = NULL;
    long    fxsig_n;

    double* merge_dates = NULL;

    double* sigtime = NULL;
    double* sigdom  = NULL;
    double* sigfor  = NULL;
    double* sigfx   = NULL;

    double quantocorr;
    double domforcorr;
    double domfxcorr;

    int nb_merge_dates;

    Err err = NULL;

    SwapDP  swapdp;
    long    float_nb_dates, float_nb_pay_dates;
    long *  float_fixing_dates = NULL, *float_pay_dates = NULL;
    long *  float_start_dates = NULL, *float_end_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL;

    //---------------------------------------------------------------------------
    //------------------------Get FX Underlying----------------------------------
    //---------------------------------------------------------------------------
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

    //---------------------------------------------------------------------------
    //------------------------Get Domestic Underlying----------------------------
    //---------------------------------------------------------------------------
    domname = get_domname_from_fxund(fx_und);
    dom_und = lookup_und(domname);
    if (!dom_und)
    {
        err = serror("Couldn't find underlying named %s", domname);
        goto FREE_RETURN;
    }
    dom_ts = get_ts_from_irund(dom_und);

    //---------------------------------------------------------------------------
    //------------------------Get Domestic Underlying----------------------------
    //---------------------------------------------------------------------------
    forname = get_forname_from_fxund(fx_und);
    for_und = lookup_und(forname);
    if (!for_und)
    {
        err = serror("Couldn't find underlying named %s", forname);
        goto FREE_RETURN;
    }
    for_ts = get_ts_from_irund(for_und);

    //---------------------------------------------------------------------------
    //------------------------Get Correlation Matrix-----------------------------
    //---------------------------------------------------------------------------
    sCorrlist = srt_f_GetTheCorrelationList();
    if (!sCorrlist->head->element)
    {
        err = "correlation list improperly initialised";
        goto FREE_RETURN;
    }

    //---------------------------------------------------------------------------
    //------------------------Get All Term Structures----------------------------
    //---------------------------------------------------------------------------
    err = Get_FX_StochRate_TermStructures(
        und3dfx,
        &domsigtime,
        &domsig,
        &domsig_n,
        &domtautime,
        &domtau,
        &domtau_n,
        &forsigtime,
        &forsig,
        &forsig_n,
        &fortautime,
        &fortau,
        &fortau_n,
        &fxsigtime,
        &fxsig,
        &fxsig_n,
        &domforcorr,
        &domfxcorr,
        &quantocorr);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = get_unique_lambda(domtau, domtau_n, &dom_lam);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = get_unique_lambda(fortau, fortau_n, &for_lam);
    if (err)
    {
        goto FREE_RETURN;
    }

    //---------------------------------------------------------------------------
    //----------------------Merge All Term Structures----------------------------
    //---------------------------------------------------------------------------
    merge_dates = (double*)calloc(domsig_n, sizeof(double));
    memcpy(merge_dates, domsigtime, domsig_n * sizeof(double));
    nb_merge_dates = domsig_n;
    num_f_concat_vector(&nb_merge_dates, &merge_dates, forsig_n, forsigtime);
    num_f_concat_vector(&nb_merge_dates, &merge_dates, fxsig_n, fxsigtime);
    num_f_sort_vector(nb_merge_dates, merge_dates);
    num_f_unique_vector(&nb_merge_dates, merge_dates);

    //----------------------Fill the new term structures-------------------------
    sigdom = (double*)calloc(nb_merge_dates, sizeof(double));
    sigfor = (double*)calloc(nb_merge_dates, sizeof(double));
    sigfx  = (double*)calloc(nb_merge_dates, sizeof(double));

    if (!sigdom || !sigfor || !sigfx)
    {
        err = "Memory allocation error in LGM1FQuantoFwd";
        goto FREE_RETURN;
    }

    for (i = nb_merge_dates - 1; i >= 0; i--)
    {
        sigdom[i] = find_sig(merge_dates[i], dom_ts);
        sigfor[i] = find_sig(merge_dates[i], for_ts);
        sigfx[i]  = find_fx_sig(merge_dates[i], fx_ts);
    }

    //---------------------------------------------------------------------------
    //-------------------Get yield curves and spot date--------------------------
    //---------------------------------------------------------------------------
    dom_yc = get_ycname_from_irund(dom_und);
    for_yc = get_ycname_from_irund(for_und);

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    //---------------------------------------------------------------------------
    //------------------------Generate Schedule----------------------------------
    //---------------------------------------------------------------------------
    err = swp_f_initSwapDP(start, end, freq, basis, &swapdp);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp,
        today,
        refrate,
        &float_pay_dates,
        &float_nb_pay_dates,
        &float_fixing_dates,
        &float_start_dates,
        &float_end_dates,
        &float_cvgs,
        &float_spreads,
        &float_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }

    *nbFwd       = float_nb_dates;
    (*QuantoFwd) = (double*)calloc(float_nb_dates, sizeof(double));
    if (!(*QuantoFwd))
    {
        err = "Memory allocation error in LGM1FQuantoFwd";
        goto FREE_RETURN;
    }

    //---------------------------------------------------------------------------
    //------------------------Compute Quanto Fwds--------------------------------
    //---------------------------------------------------------------------------

    for (i = 0; i < float_nb_dates; ++i)
    {
        QAdj = LGM1FQuantoFwdAdjustment(
            (float_fixing_dates[i] - today) / 365.0,
            (float_start_dates[i] - today) / 365.0,
            (float_end_dates[i] - today) / 365.0,
            merge_dates,
            nb_merge_dates,
            sigdom,
            sigfor,
            sigfx,
            dom_lam,
            for_lam,
            domforcorr,
            quantocorr);

        (*QuantoFwd)[i] = (QAdj * swp_f_df(today, float_start_dates[i], for_yc) /
                               swp_f_df(today, float_end_dates[i], for_yc) -
                           1.0) /
                              float_cvgs[i] +
                          float_spreads[i];
    }

FREE_RETURN:

    if (domsigtime)
        free(domsigtime);
    if (domsig)
        free(domsig);
    if (domtau)
        free(domtau);
    if (domtautime)
        free(domtautime);

    if (forsigtime)
        free(forsigtime);
    if (forsig)
        free(forsig);
    if (fortau)
        free(fortau);
    if (fortautime)
        free(fortautime);

    if (fxsigtime)
        free(fxsigtime);
    if (fxsig)
        free(fxsig);

    if (merge_dates)
        free(merge_dates);

    if (sigtime)
        free(sigtime);
    if (sigdom)
        free(sigdom);
    if (sigfor)
        free(sigfor);
    if (sigfx)
        free(sigfx);

    if (sigtime)
        free(sigtime);

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_pay_dates)
        free(float_pay_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);

    return err;
}
