/*	==========================================================================
        FILE_NAME:	Fx3FBetaDLM.c

        PURPOSE:	Modelling of the FFX(t,T*) as a log-quadratic function of a gaussian
                                variable:

                                FFX(t, T*) = FFX(0, T*) * Exp[A(t) + B(t)*Xt + C(T)*Xt^2]

                                The IR follow classical LGM 1 Factor dynamics

        DATE:		10/15/03

        AUTHOR:		L.C.
        ========================================================================== */

#include "Fx3FBetaDLMUtil.h"

#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "math.h"

#define MIN_FXBETADLM_HERMITE 20
#define INV_SQRT_TWO_PI 0.3989422804014327 /**  1/sqrt(2*Pi)  **/

#define DELTA_RHO 0.05

Err FxBetaDLM_get_C0_from_beta(double TStar, double beta, double* C0)
{
    *C0 = beta;

    return NULL;
}

Err FxBetaDLM_get_beta_from_C0(double TStar, double C0, double* beta)
{
    *beta = C0;

    return NULL;
}

Err FxBetaDLM_get_C0_from_EquiBeta(double TStar, double B0, double beta, double C0[3])
{
    double Guess;
    double Function, Derivative, Error;
    double var;
    int    counter, Limit;
    Err    err;

    err = FxBetaDLM_get_limit_beta(TStar, B0, &C0[2]);

    if (beta < C0[2])
    {
        err =
            "FxBetaDLM_get_C0_from_EquiBeta: the beta required is too low for the TStar and B0 "
            "entered";
        C0[0] = 0;
        C0[1] = 0;
    }
    var   = 0.01 * TStar;
    Error = 10E-4;
    Limit = 30;

    if (beta == 1.0)
    {
        C0[0] = 0.0;
        C0[1] = 0.0;
    }
    else if (TStar == 0.0)
    {
        C0[0] = B0 * B0 * (beta - 1.0) / 2.0;
        C0[1] = C0[0];
    }
    else
    {
        /* First Solution */
        counter = 0;
        Guess   = (1.0 - sqrt(1.0 - 4.0 * var * B0 * B0 * (beta - 1.0) * (beta - 1.0))) /
                (4.0 * var * (beta - 1.0));
        do
        {
            if (1.0 + 2.0 * Guess * var < 0.0)
            {
                Guess = -1.0 / (2.0 * var) + 0.001;
            }
            Function =
                (B0 * B0 + 2.0 * Guess * log(1.0 + 2.0 * Guess * var)) * (beta - 1.0) - 2.0 * Guess;
            Derivative = (beta - 1.0) * (2.0 * log(1.0 + 2.0 * Guess * var) +
                                         4.0 * Guess * var / (1.0 + 2.0 * Guess * var)) -
                         2.0;
            Guess = Guess - Function / Derivative;
            counter++;
        } while (fabs(Function) > Error && counter < Limit);

        C0[0] = Guess;
        if (counter >= Limit)
            err = "FxBetaDLM_get_C0_from_EquiBeta: did not converge";

        /* Second Solution */
        counter = 0;
        Guess   = (1.0 + sqrt(1.0 - 4.0 * var * B0 * B0 * (beta - 1.0) * (beta - 1.0))) /
                (4.0 * var * (beta - 1.0));

        do
        {
            if (1.0 + 2.0 * Guess * var < 0.0)
            {
                Guess = -1.0 / (2.0 * var) + 10E-10;
            }
            Function =
                (B0 * B0 + 2.0 * Guess * log(1.0 + 2.0 * Guess * var)) * (beta - 1.0) - 2.0 * Guess;
            Derivative = (beta - 1.0) * (2.0 * log(1.0 + 2.0 * Guess * var) +
                                         4.0 * Guess * var / (1.0 + 2.0 * Guess * var)) -
                         2.0;
            Guess = Guess - Function / Derivative;
            counter++;
        } while (fabs(Function) > Error && counter < Limit);

        C0[1] = Guess;
        if (counter >= Limit)
            err = "FxBetaDLM_get_C0_from_EquiBeta: did not converge";
    }
    return err;
}

Err FxBetaDLM_get_EquiBeta_from_C0(double TStar, double B0, double C0, double* beta)
{
    double var;
    var   = 0.01 * TStar;
    *beta = 1.0 + 2.0 * C0 / (B0 * B0 + 2.0 * C0 * log(1.0 + 2.0 * C0 * var));
    return NULL;
}

Err FxBetaDLM_get_limit_beta(double TStar, double B0, double* BetaLimit)
{
    double var;
    double C0;
    Err    err;

    var = 0.01 * TStar;

    if (TStar == 0.0)
        *BetaLimit = -10E6;
    else
    {
        C0 = 0.25 * (B0 * B0 * var - sqrt(B0 * B0 * var * B0 * B0 * var + 4 * B0 * B0 * var)) / var;
        err = FxBetaDLM_get_EquiBeta_from_C0(TStar, B0, C0, BetaLimit);
    }
    return NULL;
}

void FxBetaDLM_free_str(FXBETADLM_STR str)
{
    if (str)
    {
        if (str->time_dom)
            free(str->time_dom);
        str->time_dom = NULL;

        if (str->sig_dom)
            free(str->sig_dom);
        str->sig_dom = NULL;

        if (str->time_for)
            free(str->time_for);
        str->time_for = NULL;

        if (str->sig_for)
            free(str->sig_for);
        str->sig_for = NULL;

        if (str->time_fx)
            free(str->time_fx);
        str->time_fx = NULL;

        if (str->sig_fx)
            free(str->sig_fx);
        str->sig_fx = NULL;

        if (str->time_3F_corr)
            free(str->time_3F_corr);
        str->time_3F_corr = NULL;

        if (str->dom_for_3F_corr)
            free(str->dom_for_3F_corr);
        str->dom_for_3F_corr = NULL;

        if (str->dom_fx_3F_corr)
            free(str->dom_fx_3F_corr);
        str->dom_fx_3F_corr = NULL;

        if (str->for_fx_3F_corr)
            free(str->for_fx_3F_corr);
        str->for_fx_3F_corr = NULL;

        if (str->time_DLM_corr)
            free(str->time_DLM_corr);
        str->time_DLM_corr = NULL;

        if (str->dom_for_DLM_corr)
            free(str->dom_for_DLM_corr);
        str->dom_for_DLM_corr = NULL;

        if (str->dom_fx_DLM_corr)
            free(str->dom_fx_DLM_corr);
        str->dom_fx_DLM_corr = NULL;

        if (str->for_fx_DLM_corr)
            free(str->for_fx_DLM_corr);
        str->for_fx_DLM_corr = NULL;

        if (str->dom_fx_DLM_corr_down)
            free(str->dom_fx_DLM_corr_down);
        str->dom_fx_DLM_corr_down = NULL;

        if (str->for_fx_DLM_corr_down)
            free(str->for_fx_DLM_corr_down);
        str->for_fx_DLM_corr_down = NULL;

        if (str->yc_dom)
            free(str->yc_dom);
        str->yc_dom = NULL;

        if (str->yc_for)
            free(str->yc_for);
        str->yc_for = NULL;
    }
}

Err FxBetaDLM_fill_str(
    char* undname_fx,  /* und name */
    char* undname_dom, /* domestic underlying name */
    char* undname_for, /* foreign underlying name */

    /*	FX Underlying	*/
    double  tstar,
    double  B0,
    double  C0,
    double  alpha,
    double  lambda,
    double  spot_fx,
    int     nb_fx,
    long*   date_fx,
    double* sigma_fx,

    /*	Correlations	*/
    /*	3 Factor correlations */
    int     compute_3F_corr,
    int     nb_3F_corr,
    long*   date_3F_corr,
    double* dom_for_3F_corr,
    double* dom_fx_3F_corr,
    double* for_fx_3F_corr,

    /*	Model correlations */
    int    compute_DLM_corr,
    double max_time_correl,
    double min_time_space_correl,

    int     nb_DLM_corr,
    long*   date_DLM_corr,
    double* dom_for_DLM_corr,
    double* dom_fx_DLM_corr,
    double* for_fx_DLM_corr,

    /* For Alpha case */
    double* dom_fx_DLM_corr_down,
    double* for_fx_DLM_corr_down,

    FXBETADLM_STR             str,
    FxBetaDLM_OptNumerParams* NumParams)
{
    SrtUndPtr      und_dom, und_for;
    SrtCorrLstPtr  cls;
    SrtLst*        ls;
    SrtCorrLstVal* corrval;
    char *         name_dom, *name_for;
    long           today_dom, today_for;
    char*          temp;
    int            i;
    Err            err = NULL;

    /* Get and check the today */

    und_dom = lookup_und(undname_dom);
    if (!und_dom)
    {
        err = serror("Couldn't fin underlying named %s", undname_dom);
        goto FREE_RETURN;
    }

    today_dom = get_today_from_underlying(und_dom);

    und_for = lookup_und(undname_for);
    if (!und_for)
    {
        err = serror("Couldn't fin underlying named %s", undname_for);
        goto FREE_RETURN;
    }

    today_for = get_today_from_underlying(und_for);

    if (today_dom != today_for)
    {
        err = serror("Today don't match in %s and %s", undname_dom, undname_for);
        goto FREE_RETURN;
    }

    str->today = today_dom;

    /* Get and save domestic underlying */
    err = Get_LGM_TermStructure2(
        undname_dom, &(str->time_dom), &(str->sig_dom), &(str->nb_dom), &(str->tau_dom));

    if (err)
        goto FREE_RETURN;

    str->lam_dom = 1.0 / str->tau_dom;

    temp        = get_ycname_from_underlying(und_dom);
    str->yc_dom = calloc(strlen(temp) + 1, sizeof(char));

    if (!str->yc_dom)
    {
        err = "Memory allocation faillure in FxBetaDLM_fill_str";
        goto FREE_RETURN;
    }

    strcpy(str->yc_dom, temp);

    /* Get and save foreign underlying */
    err = Get_LGM_TermStructure2(
        undname_for, &(str->time_for), &(str->sig_for), &(str->nb_for), &(str->tau_for));

    if (err)
        goto FREE_RETURN;

    str->lam_for = 1.0 / str->tau_for;

    temp        = get_ycname_from_underlying(und_for);
    str->yc_for = calloc(strlen(temp) + 1, sizeof(char));

    if (!str->yc_for)
    {
        err = "Memory allocation faillure in FxBetaDLM_fill_str";
        goto FREE_RETURN;
    }

    strcpy(str->yc_for, temp);

    /* Save Fx volatility */

    str->spot_fx   = spot_fx;
    str->nb_fx     = nb_fx;
    str->equi_beta = C0;
    str->tstar     = tstar;

    str->b0     = B0;
    str->c0     = C0;
    str->alpha  = alpha;
    str->lambda = lambda;
    // err = FxBetaDLM_get_C0_from_beta(tstar, beta, &str->c0);

    str->time_fx = calloc(nb_fx, sizeof(double));
    str->sig_fx  = calloc(nb_fx, sizeof(double));

    if (!str->time_fx || !str->sig_fx)
    {
        err = "Memory allocation faillure (1) in FxBetaDLM_fill_str";
        goto FREE_RETURN;
    }

    for (i = 0; i < str->nb_fx; i++)
    {
        str->time_fx[i] = (date_fx[i] - str->today) * YEARS_IN_DAY;
        str->sig_fx[i]  = sigma_fx[i];
    }

    /* Save the correlations */

    if (compute_3F_corr && compute_DLM_corr)
    {
        err = "Need to provide at least one correlation set (3F or DLM) in FxBetaDLM_fill_str";
        goto FREE_RETURN;
    }

    /* 3 Factor Correlations */

    if (!compute_3F_corr)
    {
        if (nb_3F_corr == 0)
        {
            name_dom = get_underlying_name(und_dom);
            name_for = get_underlying_name(und_for);

            cls = srt_f_GetTheCorrelationList();
            if (!cls->head->element)
            {
                err = "Correlation list improperly initialised";
                goto FREE_RETURN;
            }

            /* Take the correlations from the static correlation TS of SORT ! */
            ls = cls->head;
            while (ls != NULL)
            {
                nb_3F_corr++;
                ls = ls->next;
            }

            str->nb_3F_corr      = nb_3F_corr;
            str->time_3F_corr    = calloc(nb_3F_corr, sizeof(double));
            str->dom_for_3F_corr = calloc(nb_3F_corr, sizeof(double));
            str->dom_fx_3F_corr  = calloc(nb_3F_corr, sizeof(double));
            str->for_fx_3F_corr  = calloc(nb_3F_corr, sizeof(double));

            if (!str->time_3F_corr || !str->dom_for_3F_corr || !str->dom_fx_3F_corr ||
                !str->for_fx_3F_corr)
            {
                err = "Memory allocation faillure (2) in FxBetaDLM_fill_str";
                goto FREE_RETURN;
            }

            ls = cls->head;

            for (i = 0; i < str->nb_3F_corr; i++)
            {
                corrval              = (SrtCorrLstVal*)ls->element->val.pval;
                str->time_3F_corr[i] = corrval->time;

                err = srt_f_get_corr_from_CorrList(
                    cls, name_dom, name_for, str->time_3F_corr[i], &(str->dom_for_3F_corr[i]));

                if (err)
                    goto FREE_RETURN;

                err = srt_f_get_corr_from_CorrList(
                    cls, name_dom, undname_fx, str->time_3F_corr[i], &(str->dom_fx_3F_corr[i]));

                if (err)
                    goto FREE_RETURN;

                err = srt_f_get_corr_from_CorrList(
                    cls, name_for, undname_fx, str->time_3F_corr[i], &(str->for_fx_3F_corr[i]));

                if (err)
                    goto FREE_RETURN;
            }
        }
        else
        {
            str->nb_3F_corr      = nb_3F_corr;
            str->time_3F_corr    = calloc(nb_3F_corr, sizeof(double));
            str->dom_for_3F_corr = calloc(nb_3F_corr, sizeof(double));
            str->dom_fx_3F_corr  = calloc(nb_3F_corr, sizeof(double));
            str->for_fx_3F_corr  = calloc(nb_3F_corr, sizeof(double));

            if (!str->time_3F_corr || !str->dom_for_3F_corr || !str->dom_fx_3F_corr ||
                !str->for_fx_3F_corr)
            {
                err = "Memory allocation faillure (2) in FxBetaDLM_fill_str";
                goto FREE_RETURN;
            }

            for (i = 0; i < str->nb_3F_corr; i++)
            {
                str->time_3F_corr[i]    = (date_3F_corr[i] - str->today) * YEARS_IN_DAY;
                str->dom_for_3F_corr[i] = dom_for_3F_corr[i];
                str->dom_fx_3F_corr[i]  = dom_fx_3F_corr[i];
                str->for_fx_3F_corr[i]  = for_fx_3F_corr[i];
            }
        }
    }

    /* DLM Correlations */
    if (!compute_DLM_corr)
    {
        if (nb_DLM_corr > 0)
        {
            /* Fill the structure */
            str->nb_DLM_corr          = nb_DLM_corr;
            str->time_DLM_corr        = calloc(nb_DLM_corr, sizeof(double));
            str->dom_for_DLM_corr     = calloc(nb_DLM_corr, sizeof(double));
            str->dom_fx_DLM_corr      = calloc(nb_DLM_corr, sizeof(double));
            str->for_fx_DLM_corr      = calloc(nb_DLM_corr, sizeof(double));
            str->dom_fx_DLM_corr_down = calloc(nb_DLM_corr, sizeof(double));
            str->for_fx_DLM_corr_down = calloc(nb_DLM_corr, sizeof(double));

            if (!str->time_DLM_corr || !str->dom_for_DLM_corr || !str->dom_fx_DLM_corr ||
                !str->for_fx_DLM_corr || !str->dom_fx_DLM_corr_down || !str->for_fx_DLM_corr_down)
            {
                err = "Memory allocation faillure (2) in FxBetaDLM_fill_str";
                goto FREE_RETURN;
            }

            for (i = 0; i < str->nb_DLM_corr; i++)
            {
                str->time_DLM_corr[i]        = (date_DLM_corr[i] - str->today) * YEARS_IN_DAY;
                str->dom_for_DLM_corr[i]     = dom_for_DLM_corr[i];
                str->dom_fx_DLM_corr[i]      = dom_fx_DLM_corr[i];
                str->for_fx_DLM_corr[i]      = for_fx_DLM_corr[i];
                str->dom_fx_DLM_corr_down[i] = dom_fx_DLM_corr_down[i];
                str->for_fx_DLM_corr_down[i] = for_fx_DLM_corr_down[i];
            }
        }
        else
        {
            /* Allocate */
            str->nb_DLM_corr          = nb_DLM_corr;
            str->time_DLM_corr        = NULL;
            str->dom_for_DLM_corr     = NULL;
            str->dom_fx_DLM_corr      = NULL;
            str->for_fx_DLM_corr      = NULL;
            str->dom_fx_DLM_corr_down = NULL;
            str->for_fx_DLM_corr_down = NULL;
        }
    }

    if (compute_3F_corr)
    {
        FxBetaDLM_get3FCorrelFromDLM(str);
    }

    if (compute_DLM_corr)
    {
        FxBetaDLM_getDLMCorrelFrom3F(str, NumParams);
    }

FREE_RETURN:

    if (err)
    {
        FxBetaDLM_free_str(str);
    }

    return err;
}

Err FxBetaDLM_free_und(SrtUndPtr pUndDesc)
{
    SrtFxDesc* pSrtFxPtr;

    pSrtFxPtr = (SrtFxDesc*)(pUndDesc->spec_desc);

    if (pSrtFxPtr->spec)
    {
        FxBetaDLM_free_str(pSrtFxPtr->spec);
        free(pSrtFxPtr->spec);
    }

    free(pUndDesc->spec_desc);
    free(pUndDesc);
    pUndDesc = NULL;
    return NULL;
}

/* Initialisation of the underlying in memory */
Err SrtInitFXBetaDLMUnd(
    char* undname_fx,  /* und name */
    char* undname_dom, /* domestic underlying name */
    char* undname_for, /* foreign underlying name */

    /*	FX DLM Underlying	*/
    double  tstar,
    double  B0,
    double  C0,
    double  alpha,
    double  lambda,
    double  spot_fx,
    int     nb_fx,
    long*   date_fx,
    double* sigma_fx,

    /*	Correlations	*/
    /*	3 Factor correlations */
    int     compute_3F_corr,
    int     nb_3F_corr,
    long*   date_3F_corr,
    double* dom_for_3F_corr,
    double* dom_fx_3F_corr,
    double* for_fx_3F_corr,

    /*	Model correlations */
    int                       compute_DLM_corr,
    int                       nb_DLM_corr,
    long*                     date_DLM_corr,
    double*                   dom_for_DLM_corr,
    double*                   dom_fx_DLM_corr,
    double*                   for_fx_DLM_corr,
    double*                   dom_fx_DLM_corr_down,
    double*                   for_fx_DLM_corr_down,
    FxBetaDLM_OptNumerParams* NumParams)
{
    FXBETADLM_STR    str;
    FxBetaDLM_model* model;
    SrtUndListPtr    und_list;
    SrtUndPtr        pUndDesc;
    SrtCurvePtr      pdomYieldCurve;
    SrtFxDesc*       und_fx;

    Err err = NULL;

    /* Allocation of the structure */
    str      = calloc(1, sizeof(FxBetaDLM_str));
    und_fx   = calloc(1, sizeof(SrtFxDesc));
    pUndDesc = (SrtUndPtr)calloc(1, sizeof(SrtUndDesc));
    model    = calloc(1, sizeof(FxBetaDLM_model));

    if (!str || !und_fx || !pUndDesc || !model)
    {
        err = "Memory allocation faillure (1) in SrtInitFXBETADLMUnd";
        goto FREE_RETURN;
    }

    /* Fill the structure */
    err = FxBetaDLM_fill_str(
        undname_fx,  /* und name */
        undname_dom, /* domestic underlying name */
        undname_for, /* foreign underlying name */

        /*	FX Underlying	*/
        tstar,
        B0,
        C0,
        alpha,
        lambda,
        spot_fx,
        nb_fx,
        date_fx,
        sigma_fx,

        /*	Correlations	*/
        /*	3 Factor correlations */
        compute_3F_corr,
        nb_3F_corr,
        date_3F_corr,
        dom_for_3F_corr,
        dom_fx_3F_corr,
        for_fx_3F_corr,

        /*	Model correlations */
        compute_DLM_corr,
        30.0, /* max_time_correl */
        1.0,  /* min_time_space_correl */
        nb_DLM_corr,
        date_DLM_corr,
        dom_for_DLM_corr,
        dom_fx_DLM_corr,
        for_fx_DLM_corr,

        /* For Alpha case */
        dom_fx_DLM_corr_down,
        for_fx_DLM_corr_down,

        str,
        NumParams);

    if (err)
        goto FREE_RETURN;

    strcpy(pUndDesc->underl_name, undname_fx);
    strupper(pUndDesc->underl_name);
    strip_white_space(pUndDesc->underl_name);
    strcpy(pUndDesc->underl_lbl, "FOREX_UND");

    pdomYieldCurve        = lookup_curve(str->yc_dom);
    pUndDesc->underl_ccy  = get_curve_ccy(pdomYieldCurve);
    pUndDesc->underl_type = FOREX_UND;

    strcpy(und_fx->mdl_lbl, "FXBETADLM_UND");
    und_fx->mdl_type = FX_BETADLM;
    und_fx->mdl_dim  = 3;
    und_fx->spot     = spot_fx;

    rem_tick_string(undname_dom, und_fx->dom_name);
    rem_tick_string(undname_for, und_fx->for_name);

    und_fx->spec        = str;
    pUndDesc->spec_desc = und_fx;

    und_list = get_underlying_list();

    err = srt_f_lstins(
        und_list,
        pUndDesc->underl_name,
        0.0,
        OBJ_PTR_UND,
        (void*)pUndDesc,
        &FxBetaDLM_free_und,
        &(pUndDesc->underl_ticker));

    if (err)
        goto FREE_RETURN;

FREE_RETURN:

    if (err)
    {
        if (str)
        {
            FxBetaDLM_free_str(str);
            free(str);
            str = NULL;
        }
    }

    if (model)
    {
        free_FxBetaDLM_model(model);
        free(model);
    }

    return err;
}

Err FxBetaDLM_get_struct_from_und(char* und, FXBETADLM_STR* str)
{
    Err        err = NULL;
    SrtUndPtr  pUndDesc;
    SrtFxDesc* pFxDesc;

    pUndDesc = lookup_und(und);
    if (!pUndDesc)
    {
        err = serror("Underlying %s not defined", und);
        return err;
    }

    if (!(ISUNDTYPE(pUndDesc, FOREX_UND)))
    {
        err = serror("Underlying %s is not a Forex Underlying", und);
    }

    if (get_mdltype_from_fxund(pUndDesc) != FX_BETADLM)
    {
        err = serror("Underlying %s is not of type FX_BETADLM", und);
        return err;
    }

    // Extract the information from the underlying
    pFxDesc = (SrtFxDesc*)(pUndDesc->spec_desc);
    *str    = (FxBetaDLM_str*)(pFxDesc->spec);

    return err;
}

Err FxBetaDLM_Get_TermStruct(
    char*    und,
    long*    today,
    double*  tau_dom,
    int*     nb_dom,
    double** time_dom,
    double** sig_dom,
    double*  tau_for,
    int*     nb_for,
    double** time_for,
    double** sig_for,
    double*  tstar,
    int*     nb_fx,
    double** time_fx,
    double** sig_fx,
    int*     nb_3F_corr,
    double** time_3F_corr,
    double** dom_for_3F_corr,
    double** dom_fx_3F_corr,
    double** for_fx_3F_corr,
    int*     nb_DLM_corr,
    double** time_DLM_corr,
    double** dom_for_DLM_corr,
    double** dom_fx_DLM_corr,
    double** for_fx_DLM_corr,
    double** dom_fx_DLM_corr_down,
    double** for_fx_DLM_corr_down,
    double*  equi_beta)
{
    FXBETADLM_STR str = NULL;
    Err           err = NULL;

    err = FxBetaDLM_get_struct_from_und(und, &str);

    if (err)
        goto FREE_RETURN;

    *today = str->today;

    *tau_dom     = str->tau_dom;
    *nb_dom      = str->nb_dom;
    *tau_for     = str->tau_for;
    *nb_for      = str->nb_for;
    *nb_fx       = str->nb_fx;
    *nb_3F_corr  = str->nb_3F_corr;
    *nb_DLM_corr = str->nb_DLM_corr;

    *equi_beta = str->equi_beta;
    *tstar     = str->tstar;

    /* Memory allocations */
    *time_dom             = calloc(*nb_dom, sizeof(double));
    *sig_dom              = calloc(*nb_dom, sizeof(double));
    *time_for             = calloc(*nb_for, sizeof(double));
    *sig_for              = calloc(*nb_for, sizeof(double));
    *time_fx              = calloc(*nb_fx, sizeof(double));
    *sig_fx               = calloc(*nb_fx, sizeof(double));
    *time_3F_corr         = calloc(*nb_3F_corr, sizeof(double));
    *dom_for_3F_corr      = calloc(*nb_3F_corr, sizeof(double));
    *dom_fx_3F_corr       = calloc(*nb_3F_corr, sizeof(double));
    *for_fx_3F_corr       = calloc(*nb_3F_corr, sizeof(double));
    *time_DLM_corr        = calloc(*nb_DLM_corr, sizeof(double));
    *dom_for_DLM_corr     = calloc(*nb_DLM_corr, sizeof(double));
    *dom_fx_DLM_corr      = calloc(*nb_DLM_corr, sizeof(double));
    *for_fx_DLM_corr      = calloc(*nb_DLM_corr, sizeof(double));
    *dom_fx_DLM_corr_down = calloc(*nb_DLM_corr, sizeof(double));
    *for_fx_DLM_corr_down = calloc(*nb_DLM_corr, sizeof(double));

    if (!time_dom || !sig_dom || !time_for || !sig_for || !time_fx || !sig_fx || !time_3F_corr ||
        !time_DLM_corr || !dom_for_3F_corr || !dom_fx_3F_corr || !for_fx_3F_corr ||
        !dom_for_DLM_corr || !dom_fx_DLM_corr || !for_fx_DLM_corr || !dom_fx_DLM_corr_down ||
        !for_fx_DLM_corr_down)
    {
        err = "Memory allocation faillure in FxBetaDLM_Get_TermStruct";
        goto FREE_RETURN;
    }

    memcpy(*time_dom, str->time_dom, *nb_dom * sizeof(double));
    memcpy(*sig_dom, str->sig_dom, *nb_dom * sizeof(double));
    memcpy(*time_for, str->time_for, *nb_for * sizeof(double));
    memcpy(*sig_for, str->sig_for, *nb_for * sizeof(double));
    memcpy(*time_fx, str->time_fx, *nb_fx * sizeof(double));
    memcpy(*sig_fx, str->sig_fx, *nb_fx * sizeof(double));
    memcpy(*time_3F_corr, str->time_3F_corr, *nb_3F_corr * sizeof(double));
    memcpy(*dom_for_3F_corr, str->dom_for_3F_corr, *nb_3F_corr * sizeof(double));
    memcpy(*dom_fx_3F_corr, str->dom_fx_3F_corr, *nb_3F_corr * sizeof(double));
    memcpy(*for_fx_3F_corr, str->for_fx_3F_corr, *nb_3F_corr * sizeof(double));
    memcpy(*time_DLM_corr, str->time_DLM_corr, *nb_DLM_corr * sizeof(double));
    memcpy(*dom_for_DLM_corr, str->dom_for_DLM_corr, *nb_DLM_corr * sizeof(double));
    memcpy(*dom_fx_DLM_corr, str->dom_fx_DLM_corr, *nb_DLM_corr * sizeof(double));
    memcpy(*for_fx_DLM_corr, str->for_fx_DLM_corr, *nb_DLM_corr * sizeof(double));
    memcpy(*dom_fx_DLM_corr_down, str->dom_fx_DLM_corr_down, *nb_DLM_corr * sizeof(double));
    memcpy(*for_fx_DLM_corr_down, str->for_fx_DLM_corr_down, *nb_DLM_corr * sizeof(double));

FREE_RETURN:

    return err;
}

Err init_FxBetaDLM_model(
    double           MinTimeCorrel,
    double           MaxTimeCorrel,
    long             today,
    char*            yc_dom,
    double           tau_dom,
    int              nb_dom,
    double*          time_dom,
    double*          sig_dom,
    char*            yc_for,
    double           tau_for,
    int              nb_for,
    double*          time_for,
    double*          sig_for,
    double           spot_fx,
    double           tstar,
    double           B0,
    double           C0,
    double           alpha,
    double           lambda,
    int              nb_fx,
    double*          time_fx,
    double*          sig_fx3F,
    double*          sig_fx,
    int              nb_corr,
    double*          time_corr,
    double*          dom_for_corr,
    double*          dom_fx_corr,
    double*          for_fx_corr,
    FxBetaDLM_model* model)
{
    Err    err = NULL;
    long   spot_date;
    int    i, j, index, NbTime, n;
    double t, dt;

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    model->lToday  = today;
    model->cYcDom  = yc_dom;
    model->cYcFor  = yc_for;
    model->dB0     = B0;
    model->dC0     = C0;
    model->dAlpha  = alpha;
    model->dLambda = lambda;

    model->dSpotFx = spot_fx;
    model->dCashFx =
        spot_fx * swp_f_df(today, spot_date, yc_dom) / swp_f_df(today, spot_date, yc_for);

    err = FxBetaDLM_get_beta_from_C0(tstar, C0, &model->dEquiBeta);

    if (err)
        goto FREE_RETURN;

    model->dTStar     = tstar;
    model->dInitTStar = tstar;
    model->lTStarDate = (long)(today + tstar * DAYS_IN_YEAR + 0.5);

    model->dTauDom = tau_dom;
    if (fabs(tau_dom) < 1.0E-08)
    {
        err = "FxBetaDLM: tau cannot be null";
        goto FREE_RETURN;
    }

    model->dLambdaDom = 1.0 / tau_dom;

    model->dTauFor = tau_for;
    if (fabs(tau_for) < 1.0E-08)
    {
        err = "FxBetaDLM: tau cannot be null";
        goto FREE_RETURN;
    }

    model->dLambdaFor = 1.0 / tau_for;

    /* Merge the dates */

    model->iNbPWTime = nb_fx;
    model->dPWTime   = calloc(nb_fx, sizeof(double));

    if (!model->dPWTime)
    {
        err = "Memory allcation faillure in init_FxBetaDLM_model";
        goto FREE_RETURN;
    }

    memcpy(model->dPWTime, time_fx, nb_fx * sizeof(double));

    num_f_concat_vector(&model->iNbPWTime, &model->dPWTime, nb_dom, time_dom);
    num_f_concat_vector(&model->iNbPWTime, &model->dPWTime, nb_for, time_for);
    num_f_concat_vector(&model->iNbPWTime, &model->dPWTime, nb_corr, time_corr);
    num_f_sort_vector(model->iNbPWTime, model->dPWTime);
    num_f_unique_vector(&model->iNbPWTime, model->dPWTime);

    /* add the intermediary dates for the correlation TS */
    /* MinTimeCorrel = 1.0; */
    /* MaxTimeCorrel = 20.0; */
    if (model->dPWTime[model->iNbPWTime - 1] < MaxTimeCorrel)
    {
        num_f_add_number(&model->iNbPWTime, &model->dPWTime, MaxTimeCorrel);
    }

    NbTime = model->iNbPWTime;
    for (i = 0; i < model->iNbPWTime; i++)
    {
        if (!i)
        {
            if (model->dPWTime[i] > MinTimeCorrel)
            {
                n  = (int)((model->dPWTime[i]) / MinTimeCorrel);
                dt = (model->dPWTime[i]) / (n + 1);
                t  = dt;
                for (j = 0; j < n; j++)
                {
                    num_f_add_number(
                        &NbTime,
                        &model->dPWTime,
                        (double)(((long)(t * DAYS_IN_YEAR)) * YEARS_IN_DAY));
                    t += dt;
                }
            }
        }
        else
        {
            if (model->dPWTime[i] - model->dPWTime[i - 1] > MinTimeCorrel)
            {
                n  = (int)((model->dPWTime[i] - model->dPWTime[i - 1]) / MinTimeCorrel);
                dt = (model->dPWTime[i] - model->dPWTime[i - 1]) / (n + 1);
                t  = model->dPWTime[i - 1] + dt;
                for (j = 0; j < n; j++)
                {
                    num_f_add_number(
                        &NbTime,
                        &model->dPWTime,
                        (double)(((long)(t * DAYS_IN_YEAR)) * YEARS_IN_DAY));
                    t += dt;
                }
            }
        }
    }
    model->iNbPWTime = NbTime;
    num_f_sort_vector(model->iNbPWTime, model->dPWTime);
    num_f_unique_vector(&model->iNbPWTime, model->dPWTime);

    /* Allocate memory */
    model->dSigmaDom       = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFor       = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFx3F      = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFx        = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFxUp      = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFxDown    = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFor     = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFx      = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrForFx      = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFxDown  = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrForFxDown  = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFxInput = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrForFxInput = calloc(model->iNbPWTime, sizeof(double));

    if (!model->dSigmaDom || !model->dSigmaFor || !model->dSigmaFx3F || !model->dSigmaFx ||
        !model->dSigmaFxUp || !model->dSigmaFxDown || !model->dCorrDomFor || !model->dCorrDomFx ||
        !model->dCorrForFx || !model->dCorrDomFxDown || !model->dCorrForFxDown ||
        !model->dCorrDomFxInput || !model->dCorrForFxInput)
    {
        err = "Memory allcation faillure in init_FxBetaDLM_model";
        goto FREE_RETURN;
    }

    /* Fill the vectors */

    /* Dom */
    for (i = 0; i < model->iNbPWTime; i++)
    {
        index               = Get_Index(model->dPWTime[i], time_dom, nb_dom);
        model->dSigmaDom[i] = sig_dom[index];

        index               = Get_Index(model->dPWTime[i], time_for, nb_for);
        model->dSigmaFor[i] = sig_for[index];

        index                = Get_Index(model->dPWTime[i], time_fx, nb_fx);
        model->dSigmaFx[i]   = sig_fx[index];
        model->dSigmaFx3F[i] = sig_fx3F[index];
        if (fabs(model->dLambda) < TINY)
        {
            model->dSigmaFxUp[i]   = sig_fx[index] * exp(model->dAlpha * sqrt(time_fx[index]));
            model->dSigmaFxDown[i] = sig_fx[index] * exp(-model->dAlpha * sqrt(time_fx[index]));
        }
        else
        {
            model->dSigmaFxUp[i] =
                sig_fx[index] *
                exp(model->dAlpha * sqrt(
                                        (1.0 - exp(-2.0 * model->dLambda * time_fx[index])) /
                                        (2.0 * model->dLambda)));
            model->dSigmaFxDown[i] =
                sig_fx[index] *
                exp(-model->dAlpha * sqrt(
                                         (1.0 - exp(-2.0 * model->dLambda * time_fx[index])) /
                                         (2.0 * model->dLambda)));
        }

        index                 = Get_Index(model->dPWTime[i], time_corr, nb_corr);
        model->dCorrDomFor[i] = dom_for_corr[index];
        model->dCorrDomFx[i]  = dom_fx_corr[index];
        model->dCorrForFx[i]  = for_fx_corr[index];
    }

    memcpy(model->dCorrForFxDown, model->dCorrForFx, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrDomFxDown, model->dCorrDomFx, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrForFxInput, model->dCorrForFx, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrDomFxInput, model->dCorrDomFx, model->iNbPWTime * sizeof(double));

    /* No structure associated */
    model->str = NULL;

FREE_RETURN:

    if (err)
    {
        free_FxBetaDLM_model(model);
    }

    return err;
}

void free_FxBetaDLM_model(FxBetaDLM_model* model)
{
    if (model)
    {
        if (model->dPWTime)
            free(model->dPWTime);
        model->dPWTime = NULL;

        if (model->dSigmaDom)
            free(model->dSigmaDom);
        model->dSigmaDom = NULL;

        if (model->dSigmaFor)
            free(model->dSigmaFor);
        model->dSigmaFor = NULL;

        if (model->dSigmaFx3F)
            free(model->dSigmaFx3F);
        model->dSigmaFx3F = NULL;

        if (model->dSigmaFx)
            free(model->dSigmaFx);
        model->dSigmaFx = NULL;

        if (model->dSigmaFxUp)
            free(model->dSigmaFxUp);
        model->dSigmaFxUp = NULL;

        if (model->dSigmaFxDown)
            free(model->dSigmaFxDown);
        model->dSigmaFxDown = NULL;

        if (model->dCorrDomFor)
            free(model->dCorrDomFor);
        model->dCorrDomFor = NULL;

        if (model->dCorrDomFx)
            free(model->dCorrDomFx);
        model->dCorrDomFx = NULL;

        if (model->dCorrForFx)
            free(model->dCorrForFx);
        model->dCorrForFx = NULL;

        if (model->dCorrDomFxDown)
            free(model->dCorrDomFxDown);
        model->dCorrDomFxDown = NULL;

        if (model->dCorrForFxDown)
            free(model->dCorrForFxDown);
        model->dCorrForFxDown = NULL;

        if (model->dCorrDomFxInput)
            free(model->dCorrDomFxInput);
        model->dCorrDomFxInput = NULL;

        if (model->dCorrForFxInput)
            free(model->dCorrForFxInput);
        model->dCorrForFxInput = NULL;
    }
}

Err FxBetaDLM_Get_Model(char* und, FxBetaDLM_model* model)
{
    Err           err = NULL;
    FXBETADLM_STR str;

    err = FxBetaDLM_get_struct_from_und(und, &str);

    if (err)
        goto FREE_RETURN;

    err = FxBetaDLM_GetModelFromStr(str, model);

    if (err)
        goto FREE_RETURN;

    model->str = str;

FREE_RETURN:

    if (err)
    {
        free_FxBetaDLM_model(model);
    }

    return err;
}

void FxBetaDLM_SetDefaultOptNumerParams(FxBetaDLM_OptNumerParams* NumParams)
{
    NumParams->dLowerBound = -5.0;
    NumParams->dUpperBound = 5.0;
    NumParams->iNbPoints   = 25;
    NumParams->iMethod     = 2;
    NumParams->dAlpha      = 0.0;
    NumParams->dMinTime    = 0.1;
    NumParams->dX0         = 0.0;

    /* Const for CPD */
    NumParams->iPrecalcSmile  = 1;
    NumParams->iPriceVol      = 1;
    NumParams->iNbX           = 15;
    NumParams->dNbStdX        = 3.0;
    NumParams->iNbFwd         = 30;
    NumParams->dNbStdFwd      = 3.0;
    NumParams->iUseCxtyAdj    = 0;
    NumParams->iEquiSpacedFwd = 0;

    /* Const for Correl Mapping */
    NumParams->dMinTimeSpaceCorrel = 1.0;
    NumParams->dMaxTimeCorrel      = 0.0;
    NumParams->iMappingMethod      = 2;
    NumParams->dFFxCorrelTolerance = 0.1;

    /* For The Alpha Fudge */
    NumParams->iMidCorrelation = 0;

    /* Const for Init */
    NumParams->iInitWith3FCorrel = 1;
}

void free_FxBetaDLM_Hermite(FxBetaDLM_Hermite* hermite)
{
    if (hermite)
    {
        if (hermite->W)
            free(hermite->W);
        hermite->W = NULL;

        if (hermite->X)
            free(hermite->X);
        hermite->X = NULL;
    }
}

Err initialise_FxBetaDLM_Hermite(FxBetaDLM_OptNumerParams* NumParams, FxBetaDLM_Hermite* hermite)
{
    Err err = NULL;
    int i;

    hermite->iNbPoints = NumParams->iNbPoints;
    hermite->W         = calloc(hermite->iNbPoints + 1, sizeof(double));
    hermite->X         = calloc(hermite->iNbPoints + 1, sizeof(double));

    if (!hermite->W || !hermite->X)
    {
        err = "Memory allocation faillure in initialise_FxBetaDLM_Hermite";
        goto FREE_RETURN;
    }

    /* hermite standard
            err = HermiteStandard(hermite->X, hermite->W, nb_points);
    */
    gauleg(
        NumParams->dLowerBound, NumParams->dUpperBound, hermite->X, hermite->W, hermite->iNbPoints);

    for (i = 1; i <= hermite->iNbPoints; i++)
    {
        hermite->W[i] *= exp(-0.5 * hermite->X[i] * hermite->X[i]) * INV_SQRT_TWO_PI;
    }

FREE_RETURN:

    if (err)
    {
        free_FxBetaDLM_Hermite(hermite);
    }

    return err;
}

Err FxBetaDLM_Allocation_Precalculations(int iNbPWTime, FxBetaDLM_ModelPrecalc* Precalc)
{
    Err err = NULL;

    Precalc->iNbDone   = 0;
    Precalc->iNbPWTime = iNbPWTime;

    Precalc->dVarRatesDom      = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarRatesFor      = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarRatesFor3D    = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRates       = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRates3D     = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarFFx           = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxDom      = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxFor      = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxFor3D    = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarX             = calloc(iNbPWTime, sizeof(double));
    Precalc->dExpectFor        = calloc(iNbPWTime, sizeof(double));
    Precalc->dIntegral         = calloc(iNbPWTime, sizeof(double));
    Precalc->dIntegral2        = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRatesAdjust = calloc(iNbPWTime, sizeof(double));

    Precalc->dVarRatesFor_down      = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarRatesFor3D_down    = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRates_down       = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRates3D_down     = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarFFx_down           = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxDom_down      = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxFor_down      = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxFor3D_down    = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarX_down             = calloc(iNbPWTime, sizeof(double));
    Precalc->dExpectFor_down        = calloc(iNbPWTime, sizeof(double));
    Precalc->dIntegral_down         = calloc(iNbPWTime, sizeof(double));
    Precalc->dIntegral2_down        = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRatesAdjust_down = calloc(iNbPWTime, sizeof(double));

    Precalc->dVarRatesFor_mid      = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarRatesFor3D_mid    = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRates_mid       = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRates3D_mid     = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarFFx_mid           = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxDom_mid      = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxFor_mid      = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarFFxFor3D_mid    = calloc(iNbPWTime, sizeof(double));
    Precalc->dVarX_mid             = calloc(iNbPWTime, sizeof(double));
    Precalc->dExpectFor_mid        = calloc(iNbPWTime, sizeof(double));
    Precalc->dIntegral_mid         = calloc(iNbPWTime, sizeof(double));
    Precalc->dIntegral2_mid        = calloc(iNbPWTime, sizeof(double));
    Precalc->dCovarRatesAdjust_mid = calloc(iNbPWTime, sizeof(double));

    if (!Precalc->dVarRatesDom || !Precalc->dVarRatesFor || !Precalc->dVarRatesFor3D ||
        !Precalc->dCovarRates || !Precalc->dCovarRates3D || !Precalc->dVarFFx ||
        !Precalc->dCovarFFxDom || !Precalc->dCovarFFxFor || !Precalc->dCovarFFxFor3D ||
        !Precalc->dVarX || !Precalc->dExpectFor || !Precalc->dIntegral || !Precalc->dIntegral2 ||
        !Precalc->dCovarRatesAdjust || !Precalc->dVarRatesFor_down ||
        !Precalc->dVarRatesFor3D_down || !Precalc->dCovarRates_down ||
        !Precalc->dCovarRates3D_down || !Precalc->dVarFFx_down || !Precalc->dCovarFFxDom_down ||
        !Precalc->dCovarFFxFor_down || !Precalc->dCovarFFxFor3D_down || !Precalc->dVarX_down ||
        !Precalc->dExpectFor_down || !Precalc->dIntegral_down || !Precalc->dIntegral2_down ||
        !Precalc->dCovarRatesAdjust_down || !Precalc->dVarRatesFor_mid ||
        !Precalc->dVarRatesFor3D_mid || !Precalc->dCovarRates_mid || !Precalc->dCovarRates3D_mid ||
        !Precalc->dVarFFx_mid || !Precalc->dCovarFFxDom_mid || !Precalc->dCovarFFxFor_mid ||
        !Precalc->dCovarFFxFor3D_mid || !Precalc->dVarX_mid || !Precalc->dExpectFor_mid ||
        !Precalc->dIntegral_mid || !Precalc->dIntegral2_mid || !Precalc->dCovarRatesAdjust_mid)
    {
        err = "Memory allocation faillure in FxBetaDLM_Allocation_Precalculations";
        goto FREE_RETURN;
    }

FREE_RETURN:

    if (err)
    {
        FxBetaDLM_Free_Precalculations(Precalc);
    }

    return err;
}

void FxBetaDLM_Free_Precalculations(FxBetaDLM_ModelPrecalc* Precalc)
{
    if (Precalc)
    {
        if (Precalc->dVarRatesDom)
        {
            free(Precalc->dVarRatesDom);
            Precalc->dVarRatesDom = NULL;
        }

        if (Precalc->dVarRatesFor)
        {
            free(Precalc->dVarRatesFor);
            Precalc->dVarRatesFor = NULL;
        }

        if (Precalc->dVarRatesFor3D)
        {
            free(Precalc->dVarRatesFor3D);
            Precalc->dVarRatesFor3D = NULL;
        }

        if (Precalc->dCovarRates)
        {
            free(Precalc->dCovarRates);
            Precalc->dCovarRates = NULL;
        }

        if (Precalc->dCovarRates3D)
        {
            free(Precalc->dCovarRates3D);
            Precalc->dCovarRates3D = NULL;
        }

        if (Precalc->dVarFFx)
        {
            free(Precalc->dVarFFx);
            Precalc->dVarFFx = NULL;
        }

        if (Precalc->dCovarFFxDom)
        {
            free(Precalc->dCovarFFxDom);
            Precalc->dCovarFFxDom = NULL;
        }

        if (Precalc->dCovarFFxFor)
        {
            free(Precalc->dCovarFFxFor);
            Precalc->dCovarFFxFor = NULL;
        }

        if (Precalc->dCovarFFxFor3D)
        {
            free(Precalc->dCovarFFxFor3D);
            Precalc->dCovarFFxFor3D = NULL;
        }

        if (Precalc->dVarX)
        {
            free(Precalc->dVarX);
            Precalc->dVarX = NULL;
        }

        if (Precalc->dExpectFor)
        {
            free(Precalc->dExpectFor);
            Precalc->dExpectFor = NULL;
        }

        if (Precalc->dIntegral)
        {
            free(Precalc->dIntegral);
            Precalc->dIntegral = NULL;
        }

        if (Precalc->dIntegral2)
        {
            free(Precalc->dIntegral2);
            Precalc->dIntegral2 = NULL;
        }

        if (Precalc->dCovarRatesAdjust)
        {
            free(Precalc->dCovarRatesAdjust);
            Precalc->dCovarRatesAdjust = NULL;
        }

        if (Precalc->dVarRatesFor_down)
        {
            free(Precalc->dVarRatesFor_down);
            Precalc->dVarRatesFor_down = NULL;
        }

        if (Precalc->dVarRatesFor3D_down)
        {
            free(Precalc->dVarRatesFor3D_down);
            Precalc->dVarRatesFor3D_down = NULL;
        }

        if (Precalc->dCovarRates_down)
        {
            free(Precalc->dCovarRates_down);
            Precalc->dCovarRates_down = NULL;
        }

        if (Precalc->dCovarRates3D_down)
        {
            free(Precalc->dCovarRates3D_down);
            Precalc->dCovarRates3D_down = NULL;
        }

        if (Precalc->dVarFFx_down)
        {
            free(Precalc->dVarFFx_down);
            Precalc->dVarFFx_down = NULL;
        }

        if (Precalc->dCovarFFxDom_down)
        {
            free(Precalc->dCovarFFxDom_down);
            Precalc->dCovarFFxDom_down = NULL;
        }

        if (Precalc->dCovarFFxFor_down)
        {
            free(Precalc->dCovarFFxFor_down);
            Precalc->dCovarFFxFor_down = NULL;
        }

        if (Precalc->dCovarFFxFor3D_down)
        {
            free(Precalc->dCovarFFxFor3D_down);
            Precalc->dCovarFFxFor3D_down = NULL;
        }

        if (Precalc->dVarX_down)
        {
            free(Precalc->dVarX_down);
            Precalc->dVarX_down = NULL;
        }

        if (Precalc->dExpectFor_down)
        {
            free(Precalc->dExpectFor_down);
            Precalc->dExpectFor_down = NULL;
        }

        if (Precalc->dIntegral_down)
        {
            free(Precalc->dIntegral_down);
            Precalc->dIntegral_down = NULL;
        }

        if (Precalc->dIntegral2_down)
        {
            free(Precalc->dIntegral2_down);
            Precalc->dIntegral2_down = NULL;
        }

        if (Precalc->dCovarRatesAdjust_down)
        {
            free(Precalc->dCovarRatesAdjust_down);
            Precalc->dCovarRatesAdjust_down = NULL;
        }

        if (Precalc->dVarRatesFor_mid)
        {
            free(Precalc->dVarRatesFor_mid);
            Precalc->dVarRatesFor_mid = NULL;
        }

        if (Precalc->dVarRatesFor3D_mid)
        {
            free(Precalc->dVarRatesFor3D_mid);
            Precalc->dVarRatesFor3D_mid = NULL;
        }

        if (Precalc->dCovarRates_mid)
        {
            free(Precalc->dCovarRates_mid);
            Precalc->dCovarRates_mid = NULL;
        }

        if (Precalc->dCovarRates3D_mid)
        {
            free(Precalc->dCovarRates3D_mid);
            Precalc->dCovarRates3D_mid = NULL;
        }

        if (Precalc->dVarFFx_mid)
        {
            free(Precalc->dVarFFx_mid);
            Precalc->dVarFFx_mid = NULL;
        }

        if (Precalc->dCovarFFxDom_mid)
        {
            free(Precalc->dCovarFFxDom_mid);
            Precalc->dCovarFFxDom_mid = NULL;
        }

        if (Precalc->dCovarFFxFor_mid)
        {
            free(Precalc->dCovarFFxFor_mid);
            Precalc->dCovarFFxFor_mid = NULL;
        }

        if (Precalc->dCovarFFxFor3D_mid)
        {
            free(Precalc->dCovarFFxFor3D_mid);
            Precalc->dCovarFFxFor3D_mid = NULL;
        }

        if (Precalc->dVarX_mid)
        {
            free(Precalc->dVarX_mid);
            Precalc->dVarX_mid = NULL;
        }

        if (Precalc->dExpectFor_mid)
        {
            free(Precalc->dExpectFor_mid);
            Precalc->dExpectFor_mid = NULL;
        }

        if (Precalc->dIntegral_mid)
        {
            free(Precalc->dIntegral_mid);
            Precalc->dIntegral_mid = NULL;
        }

        if (Precalc->dIntegral2_mid)
        {
            free(Precalc->dIntegral2_mid);
            Precalc->dIntegral2_mid = NULL;
        }

        if (Precalc->dCovarRatesAdjust_mid)
        {
            free(Precalc->dCovarRatesAdjust_mid);
            Precalc->dCovarRatesAdjust_mid = NULL;
        }
    }
}

/*
Calculates E[(exp(constant_coef + linear_coef X + quadratic_coef X^2) - strike)+]
where X has a normal distribution with mean and variance
*/

Err quadratic_bs(
    double         constant_coef,
    double         linear_coef,
    double         quadratic_coef,
    double         mean,
    double         variance,
    double         strike,
    SrtCallPutType callput,
    double*        Premium)
{
    double root_1, root_2, swap_root, k1, k2, J1, J2, delta;
    Err    err = NULL;

    if (variance == 0.0)
    {
        if (callput == SRT_CALL)
        {
            *Premium =
                exp(constant_coef + linear_coef * mean + quadratic_coef * mean * mean) - strike;
            if (*Premium < 0)
            {
                *Premium = 0;
            }
        }
        else
        {
            *Premium =
                strike - exp(constant_coef + linear_coef * mean + quadratic_coef * mean * mean);
            if (*Premium < 0)
            {
                *Premium = 0;
            }
        }
    }
    else
    {
        if (quadratic_coef >= 0.5 / variance)
        {
            err = serror(
                "quadratic_bs  The process has infinite variance, C = %f > %f",
                quadratic_coef,
                0.5 / variance);
            goto FREE_RETURN;
        }
        else if (quadratic_coef == 0)
        {
            root_1 = (log(strike) - constant_coef) / linear_coef;
            if (callput == SRT_CALL)
            {
                J1 =
                    exp(constant_coef - mean * mean / (2.0 * variance) +
                        (linear_coef * variance + mean) * (linear_coef * variance + mean) /
                            (2.0 * variance));
                J1 *= norm((mean + linear_coef * variance - root_1) / sqrt(variance));

                J2 = strike * norm((mean - root_1) / sqrt(variance));

                *Premium = J1 - J2;
            }
            else
            {
                J1 =
                    exp(constant_coef - mean * mean / (2.0 * variance) +
                        (linear_coef * variance + mean) * (linear_coef * variance + mean) /
                            (2.0 * variance));
                J1 *= norm((root_1 - mean - linear_coef * variance) / sqrt(variance));

                J2 = strike * norm((root_1 - mean) / sqrt(variance));

                *Premium = J2 - J1;
            }
        }
        else
        {
            delta = linear_coef * linear_coef - 4 * quadratic_coef * (constant_coef - log(strike));
            if (delta < 0)
            {
                if (callput == SRT_CALL)
                {
                    if (quadratic_coef < 0)
                    {
                        *Premium = 0;
                    }
                    else
                    {
                        k1       = linear_coef + mean / variance;
                        k2       = sqrt(1 / variance - 2 * quadratic_coef);
                        *Premium = exp(constant_coef - mean * mean / (2.0 * variance) +
                                       k1 * k1 / (2 * k2 * k2)) /
                                   sqrt(variance) / k2;
                        *Premium += -strike;
                    }
                }
                else
                {
                    if (quadratic_coef > 0)
                    {
                        *Premium = 0;
                    }
                    else
                    {
                        k1       = linear_coef + mean / variance;
                        k2       = sqrt(1 / variance - 2 * quadratic_coef);
                        *Premium = strike;
                        *Premium += -exp(
                                        constant_coef - mean * mean / (2.0 * variance) +
                                        k1 * k1 / (2 * k2 * k2)) /
                                    sqrt(variance) / k2;
                    }
                }
            }

            else
            {
                root_1 = (-linear_coef + sqrt(delta)) / (2 * quadratic_coef);
                root_2 = (-linear_coef - sqrt(delta)) / (2 * quadratic_coef);

                /* put the roots in the right order */
                if (root_1 > root_2)
                {
                    swap_root = root_1;
                    root_1    = root_2;
                    root_2    = swap_root;
                }
                k1 = linear_coef + mean / variance;
                k2 = sqrt(1 / variance - 2 * quadratic_coef);

                if (quadratic_coef > 0)
                {
                    if (callput == SRT_CALL)
                    {
                        J1 = norm(k2 * root_1 - k1 / k2);
                        J1 += norm(k1 / k2 - k2 * root_2);

                        J2 = norm((root_1 - mean) / sqrt(variance));
                        J2 += norm((mean - root_2) / sqrt(variance));
                    }
                    else
                    {
                        J1 = norm(k2 * root_2 - k1 / k2);
                        J1 += -norm(k2 * root_1 - k1 / k2);

                        J2 = norm((root_2 - mean) / sqrt(variance));
                        J2 += -norm((root_1 - mean) / sqrt(variance));
                    }
                }
                else
                {
                    if (callput == SRT_CALL)
                    {
                        J1 = norm(k2 * root_2 - k1 / k2);
                        J1 += -norm(k2 * root_1 - k1 / k2);

                        J2 = norm((root_2 - mean) / sqrt(variance));
                        J2 += -norm((root_1 - mean) / sqrt(variance));
                    }
                    else
                    {
                        J1 = norm(k2 * root_1 - k1 / k2);
                        J1 += norm(k1 / k2 - k2 * root_2);

                        J2 = norm((root_1 - mean) / sqrt(variance));
                        J2 += norm((mean - root_2) / sqrt(variance));
                    }
                }

                J1 *=
                    exp(constant_coef - mean * mean / (2.0 * variance) + k1 * k1 / (2 * k2 * k2)) /
                    sqrt(variance) / k2;
                J2 *= strike;
                if (callput == SRT_CALL)
                {
                    *Premium = J1 - J2;
                }
                else
                {
                    *Premium = J2 - J1;
                }
            }
        }
    }
FREE_RETURN:

    return err;
}

void FxBetaDLM_Init_GRFNNumerParams(FxBetaDLM_GRFNNumerParams* NumParams)
{
    /* Default Parmas */
    NumParams->min_time    = 0.1;
    NumParams->numeraire   = 0;
    NumParams->do_pecs     = 0;
    NumParams->do_discount = 1;
}

/* Finds the correlation using F(t,T*) and checking with S */
Err FxBetaDLM_correl_mapping_forward(
    double  CorrelTolerance,
    double  t,
    double  previous_t,
    double  B0,
    double  C0,
    double  varFFx,
    double  TStar,
    double  SigmaDom,
    double  LambdaDom,
    double  SigmaFor,
    double  LambdaFor,
    double  Sigma3F,
    double  SigmaDLM,
    double  RhoDF,
    double  RhoDS,
    double  RhoFS,
    double* RhoDX,
    double* RhoFX)
{
    double gammaDom, gammaFor;
    double sigFFx, RhoDomFFx, RhoForFFx;
    double b, c;
    double delta;
    double RhoDX1, RhoFX1, RhoDomFFx1, RhoForFFx1, sigFFx1, varFFx1, B1, C1, volS1, volX2_1, RhoDS1,
        RhoFS1;
    double RhoDX2, RhoFX2, RhoDomFFx2, RhoForFFx2, sigFFx2, varFFx2, B2, C2, volS2, volX2_2, RhoDS2,
        RhoFS2;
    int admissible1 = 1, admissible2 = 1;

    Err err = NULL;

    gammaDom = -SigmaDom * (1.0 - exp(-LambdaDom * (TStar - t))) / LambdaDom;
    gammaFor = -SigmaFor * (1.0 - exp(-LambdaFor * (TStar - t))) / LambdaFor;

    sigFFx = Sigma3F * Sigma3F + gammaDom * gammaDom + gammaFor * gammaFor +
             2.0 * RhoFS * Sigma3F * gammaFor - 2.0 * RhoDS * Sigma3F * gammaDom -
             2.0 * RhoDF * gammaDom * gammaFor;

    if (sigFFx < 1.0E-16)
    {
        err = "FxBetaDLM_correl_mapping: the FFx at Tstar has negative variance";
        return err;
    }

    sigFFx = sqrt(sigFFx);

    RhoDomFFx = (RhoDS * Sigma3F + RhoDF * gammaFor - gammaDom) / sigFFx;
    RhoForFFx = (RhoFS * Sigma3F + gammaFor - RhoDF * gammaDom) / sigFFx;

    b = RhoDF * gammaFor - gammaDom + RhoDomFFx * (RhoDomFFx * gammaDom - RhoForFFx * gammaFor);
    c = -(SigmaDLM * SigmaDLM + gammaDom * gammaDom - gammaFor * gammaFor);
    c += -2.0 * RhoForFFx / RhoDomFFx * (RhoDF * gammaFor - gammaDom) * gammaFor;
    c *= RhoDomFFx * RhoDomFFx;
    c += pow(RhoDF * gammaFor - gammaDom, 2);

    delta = b * b - c;

    if (delta < 0.0)
    {
        delta = 0.0;

        /*
        err = "FxBetaDLM_correl_mapping: delta < 0 there is no solution";
        return err;
        */
    }

    /* First solution */
    RhoDX1 = -b + sqrt(delta);
    RhoFX1 = RhoForFFx / RhoDomFFx * (RhoDX1 + RhoDF * gammaFor - gammaDom) -
             (gammaFor - RhoDF * gammaDom);

    sigFFx1 = SigmaDLM * SigmaDLM + gammaDom * gammaDom + gammaFor * gammaFor +
              2.0 * RhoFX1 * gammaFor - 2.0 * RhoDX1 * gammaDom - 2.0 * RhoDF * gammaDom * gammaFor;

    if (sigFFx1 > 0.0)
    {
        sigFFx1 = sqrt(sigFFx1);

        RhoDomFFx1 = RhoDX1 + RhoDF * gammaFor - gammaDom;
        RhoForFFx1 = RhoFX1 + gammaFor - RhoDF * gammaDom;

        RhoDomFFx1 /= sigFFx1;
        RhoForFFx1 /= sigFFx1;
    }
    else
        admissible1 = 0;

    /* Second solution */
    RhoDX2 = -b - sqrt(delta);
    RhoFX2 = RhoForFFx / RhoDomFFx * (RhoDX2 + RhoDF * gammaFor - gammaDom) -
             (gammaFor - RhoDF * gammaDom);

    sigFFx2 = SigmaDLM * SigmaDLM + gammaDom * gammaDom + gammaFor * gammaFor +
              2.0 * RhoFX2 * gammaFor - 2.0 * RhoDX2 * gammaDom - 2.0 * RhoDF * gammaDom * gammaFor;

    if (sigFFx2 > 0.0)
    {
        sigFFx2 = sqrt(sigFFx2);

        RhoDomFFx2 = RhoDX2 + RhoDF * gammaFor - gammaDom;
        RhoForFFx2 = RhoFX2 + gammaFor - RhoDF * gammaDom;

        RhoDomFFx2 /= sigFFx2;
        RhoForFFx2 /= sigFFx2;
    }
    else
        admissible2 = 0;

    /* Check if the solutions are possible */
    if (admissible1 && (fabs(RhoDomFFx1 - RhoDomFFx) > CorrelTolerance ||
                        fabs(RhoForFFx1 - RhoForFFx) > CorrelTolerance))
        admissible1 = 0;

    if (admissible2 && (fabs(RhoDomFFx2 - RhoDomFFx) > CorrelTolerance ||
                        fabs(RhoForFFx2 - RhoForFFx) > CorrelTolerance))
        admissible2 = 0;

    if (admissible1 && admissible2)
    {
        /* First solution */
        varFFx1 = varFFx;
        varFFx1 += SigmaDLM * SigmaDLM * (t - previous_t);
        varFFx1 += -2.0 * RhoFX1 * SigmaFor * Etha_Func(LambdaFor, TStar, previous_t, t) +
                   2.0 * RhoDX1 * SigmaDom * Etha_Func(LambdaDom, TStar, previous_t, t);
        varFFx1 += SigmaFor * SigmaFor * Psi_Func(LambdaFor, LambdaFor, TStar, previous_t, t) +
                   SigmaDom * SigmaDom * Psi_Func(LambdaDom, LambdaDom, TStar, previous_t, t) -
                   2.0 * RhoDF * SigmaDom * SigmaFor *
                       Psi_Func(LambdaDom, LambdaFor, TStar, previous_t, t);

        C1 = 1.0 + 2.0 * C0 * varFFx1;
        B1 = B0 / C1;
        C1 = C0 / C1;

        volX2_1 = SigmaDLM * SigmaDLM;
        volX2_1 += 2.0 * (RhoFX1 * gammaFor - RhoDX1 * gammaDom);
        volX2_1 += gammaDom * gammaDom + gammaFor * gammaFor - 2.0 * RhoDF * gammaDom * gammaFor;

        volS1 = B1 * B1 * SigmaDLM * SigmaDLM +
                (B1 - 1.0) * (B1 - 1.0) * (gammaFor * gammaFor + gammaDom * gammaDom);
        volS1 += 2.0 * (B1 - 1.0) *
                 (B1 * (RhoFX1 * gammaFor - RhoDX1 * gammaDom) -
                  (B1 - 1.0) * RhoDF * gammaDom * gammaFor);
        volS1 += 4.0 * C1 * C1 * varFFx1 * volX2_1;

        if (volS1 < 0.0)
        {
            // err = "FxBetaDLM_correl_mapping: negative vol of Spot";
            // goto FREE_RETURN;
            admissible1 = 0;
        }
        else
        {
            volS1 = sqrt(volS1);

            RhoDS1 = B1 * RhoDX1 + (B1 - 1.0) * (RhoDF * gammaFor - gammaDom);
            RhoFS1 = B1 * RhoFX1 + (B1 - 1.0) * (gammaFor - RhoDF * gammaDom);

            RhoDS1 /= volS1;
            RhoFS1 /= volS1;
        }

        /* Second solution */
        varFFx2 = varFFx;
        varFFx2 += SigmaDLM * SigmaDLM * (t - previous_t);
        varFFx2 += -2.0 * RhoFX2 * SigmaFor * Etha_Func(LambdaFor, TStar, previous_t, t) +
                   2.0 * RhoDX2 * SigmaDom * Etha_Func(LambdaDom, TStar, previous_t, t);
        varFFx2 += SigmaFor * SigmaFor * Psi_Func(LambdaFor, LambdaFor, TStar, previous_t, t) +
                   SigmaDom * SigmaDom * Psi_Func(LambdaDom, LambdaDom, TStar, previous_t, t) -
                   2.0 * RhoDF * SigmaDom * SigmaFor *
                       Psi_Func(LambdaDom, LambdaFor, TStar, previous_t, t);

        C2 = 1.0 + 2.0 * C0 * varFFx2;
        B2 = B0 / C2;
        C2 = C0 / C2;

        volX2_2 = SigmaDLM * SigmaDLM;
        volX2_2 += 2.0 * (RhoFX2 * gammaFor - RhoDX2 * gammaDom);
        volX2_2 += gammaDom * gammaDom + gammaFor * gammaFor - 2.0 * RhoDF * gammaDom * gammaFor;

        volS2 = B2 * B2 * SigmaDLM * SigmaDLM +
                (B2 - 1.0) * (B2 - 1.0) * (gammaFor * gammaFor + gammaDom * gammaDom);
        volS2 += 2.0 * (B2 - 1.0) *
                 (B2 * (RhoFX2 * gammaFor - RhoDX2 * gammaDom) -
                  (B2 - 1.0) * RhoDF * gammaDom * gammaFor);
        volS2 += 4.0 * C2 * C2 * varFFx2 * volX2_2;

        if (volS2 < 0.0)
        {
            // err = "FxBetaDLM_correl_mapping: negative vol of Spot";
            // goto FREE_RETURN;
            admissible2 = 0;
        }
        else
        {
            volS2 = sqrt(volS2);

            RhoDS2 = B2 * RhoDX2 + (B2 - 1.0) * (RhoDF * gammaFor - gammaDom);
            RhoFS2 = B2 * RhoFX2 + (B2 - 1.0) * (gammaFor - RhoDF * gammaDom);

            RhoDS2 /= volS2;
            RhoFS2 /= volS2;
        }

        /* take the solution with the closest mapping */
        if (admissible1 && admissible2)
        {
            if (fabs(RhoDS - RhoDS1) < fabs(RhoDS - RhoDS2))
            {
                admissible2 = 0;
            }
            else
            {
                admissible1 = 0;
            }
        }
    }

    if (admissible1)
    {
        *RhoDX = RhoDX1 / SigmaDLM;
        *RhoFX = RhoFX1 / SigmaDLM;
    }
    else if (admissible2)
    {
        *RhoDX = RhoDX2 / SigmaDLM;
        *RhoFX = RhoFX2 / SigmaDLM;
    }
    else
    {
        err = "FxBetaDLM_correl_mapping: No solution for correlation mapping";
        goto FREE_RETURN;
    }

    RhoDX1 = *RhoDX;
    RhoFX1 = *RhoFX;

    err = FxBetaDLM_CheckOrMakePosMatrix(RhoDF, RhoDX, RhoFX);

    if (err)
        goto FREE_RETURN;

    if (RhoDX1 != *RhoDX || RhoFX1 != *RhoFX)
    {
        sigFFx1 = SigmaDLM * SigmaDLM + gammaDom * gammaDom + gammaFor * gammaFor +
                  2.0 * *RhoFX * SigmaDLM * gammaFor - 2.0 * *RhoDX * SigmaDLM * gammaDom -
                  2.0 * RhoDF * gammaDom * gammaFor;

        if (sigFFx1 > 0.0)
        {
            sigFFx1 = sqrt(sigFFx1);

            RhoDomFFx1 = *RhoDX * SigmaDLM + RhoDF * gammaFor - gammaDom;
            RhoForFFx1 = *RhoFX * SigmaDLM + gammaFor - RhoDF * gammaDom;

            RhoDomFFx1 /= sigFFx1;
            RhoForFFx1 /= sigFFx1;
        }

        if (sigFFx < 0.0 || fabs(RhoDomFFx1 - RhoDomFFx) > CorrelTolerance ||
            fabs(RhoForFFx1 - RhoForFFx) > CorrelTolerance)
        {
            err =
                "FxBetaDLM_correl_mapping: No solution for correlation mapping after adapting the "
                "correlations";
            goto FREE_RETURN;
        }
    }

FREE_RETURN:
    return err;
}

/* Finds the correlation using F(t,T*) */
Err FxBetaDLM_correl_mapping8(
    double  t,
    double  TStar,
    double  SigmaDom,
    double  LambdaDom,
    double  SigmaFor,
    double  LambdaFor,
    double  Sigma3F,
    double  SigmaDLM,
    double  RhoDF,
    double  RhoDS,
    double  RhoFS,
    double* RhoDX,
    double* RhoFX)
{
    double gammaDom, gammaFor;
    double sigFFx, RhoDomFFx, RhoForFFx;
    double RhoDomFFxTest, sigFFxTest;
    double b, c;
    double delta;
    Err    err = NULL;

    gammaDom = -SigmaDom * (1.0 - exp(-LambdaDom * (TStar - t))) / LambdaDom;
    gammaFor = -SigmaFor * (1.0 - exp(-LambdaFor * (TStar - t))) / LambdaFor;

    sigFFx = Sigma3F * Sigma3F + gammaDom * gammaDom + gammaFor * gammaFor +
             2.0 * RhoFS * Sigma3F * gammaFor - 2.0 * RhoDS * Sigma3F * gammaDom -
             2.0 * RhoDF * gammaDom * gammaFor;

    if (sigFFx < 1.0E-16)
    {
        err = "FxBetaDLM_correl_mapping: the FFx at Tstar has negative variance";
        return err;
    }

    sigFFx = sqrt(sigFFx);

    RhoDomFFx = (RhoDS * Sigma3F + RhoDF * gammaFor - gammaDom) / sigFFx;
    RhoForFFx = (RhoFS * Sigma3F + gammaFor - RhoDF * gammaDom) / sigFFx;

    b = RhoDF * gammaFor - gammaDom + RhoDomFFx * (RhoDomFFx * gammaDom - RhoForFFx * gammaFor);
    c = -(SigmaDLM * SigmaDLM + gammaDom * gammaDom - gammaFor * gammaFor);
    c += -2.0 * RhoForFFx / RhoDomFFx * (RhoDF * gammaFor - gammaDom) * gammaFor;
    c *= RhoDomFFx * RhoDomFFx;
    c += pow(RhoDF * gammaFor - gammaDom, 2);

    delta = b * b - c;

    if (delta < 0.0)
    {
        err = "FxBetaDLM_correl_mapping: delta < 0 there is no solution";
        return err;
    }

    *RhoDX = -b + sqrt(delta);
    *RhoFX = RhoForFFx / RhoDomFFx * (*RhoDX + RhoDF * gammaFor - gammaDom) -
             (gammaFor - RhoDF * gammaDom);

    sigFFxTest = SigmaDLM * SigmaDLM + gammaDom * gammaDom + gammaFor * gammaFor +
                 2.0 * *RhoFX * gammaFor - 2.0 * *RhoDX * gammaDom -
                 2.0 * RhoDF * gammaDom * gammaFor;

    RhoDomFFxTest = *RhoDX + RhoDF * gammaFor - gammaDom;

    if (RhoDomFFx * RhoDomFFxTest > 0.0)
    {
        *RhoDX /= SigmaDLM;
        *RhoFX /= SigmaDLM;
    }
    else
    {
        *RhoDX = -b - sqrt(delta);
        *RhoFX = RhoForFFx / RhoDomFFx * (*RhoDX + RhoDF * gammaFor - gammaDom) -
                 (gammaFor - RhoDF * gammaDom);
        *RhoDX /= SigmaDLM;
        *RhoFX /= SigmaDLM;
    }

    if (fabs(*RhoDX) > 1.0)
        err = "FxBetaDLM_correl_mapping: RhoDX not in [-1;1]";

    if (fabs(*RhoFX) > 1.0)
        err = "FxBetaDLM_correl_mapping: RhoFX not in [-1;1]";

    return err;
}

/* Finds the Fix Point Solution */
Err FxBetaDLM_correl_mapping(
    double  volX,
    double  B0,
    double  C0,
    double  TStar,
    double  SigmaDom,
    double  LambdaDom,
    double  SigmaFor,
    double  LambdaFor,
    double  RhoDF,
    double  RhoDS,
    double  RhoFS,
    double  previous_t,
    double  next_t,
    double  VarFFx,
    double* RhoDX,
    double* RhoFX)
{
    Err err = NULL;

    double RhoDX1, RhoFX1;
    double VarFFx_previous;
    double B, C, vol_dom, vol_for, volDom, volFor, vol_X, Const1, Const2;
    double t, middle_t;
    double a, b, c, delta;
    double tolerance = 10e-6;
    int    iteration, iteration_limit;

    /* initialisation */
    iteration       = 0;
    iteration_limit = 20;

    t        = previous_t;
    middle_t = 0.5 * (previous_t + next_t);

    VarFFx_previous = VarFFx;

    RhoDX1 = RhoDS;
    RhoFX1 = RhoFS;

    do
    {
        *RhoDX = RhoDX1;
        *RhoFX = RhoFX1;

        iteration++;

        if ((1.0 + 2.0 * C0 * VarFFx) <= 0.0)
        {
            err = "FxBetaDLM_correl_mapping: The cumulative variance is too high";
            return err;
        }
        B      = B0 / (1.0 + 2.0 * C0 * VarFFx);
        C      = C0 / (1.0 + 2.0 * C0 * VarFFx);
        Const1 = 4.0 * C * C * VarFFx;
        Const2 = Const1 / B;

        volDom = -SigmaDom * (1.0 - exp(-LambdaDom * (TStar - t))) / LambdaDom;
        volFor = -SigmaFor * (1.0 - exp(-LambdaFor * (TStar - t))) / LambdaFor;

        /* add the coefficients */
        vol_dom = -volDom * (B - 1.0);
        vol_for = volFor * (B - 1.0);
        vol_X   = volX * B;

        if (fabs(RhoDS) > 1.0E-08)
        {
            /* solve polynomial */
            a = 1.0;
            b = vol_dom + RhoDF * vol_for +
                RhoDS * (B - 1.0 + Const2) * (RhoDS * volDom - RhoFS * volFor);
            c = pow(RhoDF * vol_for + vol_dom, 2) - pow(RhoDS * volX, 2) * (B * B + Const1) -
                RhoDS * RhoDS * ((B - 1.0) * (B - 1.0) + Const1) *
                    (volFor * volFor - 2.0 * RhoDF * volDom * volFor + volDom * volDom) -
                2.0 * RhoDS * vol_for * (B - 1.0 + Const2) *
                    (RhoFS * (volFor * RhoDF - volDom) - RhoDS * (volFor - volDom * RhoDF));

            delta = b * b - a * c;

            if (delta < 0.0)
            {
                err = serror("No solution for correlation mapping, delta < 0 for t = %f", t);
                return err;
            }

            delta = sqrt(delta);

            if (RhoDS > 0.0)
            {
                /* we choose the positive solution */
                RhoDX1 = (-b + delta) / a;
            }
            else
            {
                RhoDX1 = (-b - delta) / a;
            }

            /*			if (fabs(RhoDX1) > vol_X * 0.999)
                                    {
                                            err = "No solution for correlation mapping, correlation
               dom X not in [-1;1]"; return err;
                                    }*/

            RhoFX1 =
                RhoFS / RhoDS * (vol_dom + RhoDF * vol_for + RhoDX1) - RhoDF * vol_dom - vol_for;

            /*			if (fabs(RhoFX1) > vol_X * 0.999)
                                    {
                                            err = "No solution for correlation mapping, correlation
               for X not in [-1;1]"; return err;
                                    }*/
        }
        else
        {
            RhoDX1 = -vol_dom - RhoDF * vol_for;

            /*			if (fabs(RhoDX1) > vol_X * 0.999)
                                    {
                                            err = "No solution for correlation mapping, correlation
               dom X not in [-1;1]"; return err;
                                    }*/

            a = 1.0;
            b = vol_for + RhoDF * vol_dom - RhoFS * RhoFS * vol_for;
            c = pow(RhoDF * vol_dom, 2) + vol_for * vol_for +
                RhoFS * RhoFS * (vol_dom * vol_dom - vol_for * vol_for - vol_X * vol_X) +
                2.0 * RhoDF * vol_dom * vol_for;

            delta = b * b - a * c;

            if (delta < 0.0)
            {
                err = "No solution for correlation mapping, delta < 0";
                return err;
            }

            delta = sqrt(delta);

            if (RhoFS > 0.0)
            {
                RhoFX1 = (-b + delta) / a;
            }
            else
            {
                RhoFX1 = (-b - delta) / a;
            }

            /*			if (fabs(RhoFX1) > vol_X * 0.999)
                                    {
                                            err = "No solution for correlation mapping, correlation
               For X not in [-1;1]"; return err;
                                    }*/
        }

        RhoDX1 /= vol_X;
        RhoFX1 /= vol_X;

        VarFFx = VarFFx_previous;

        VarFFx += volX * volX * (middle_t - previous_t);
        VarFFx +=
            -2.0 * RhoFX1 * SigmaFor * volX * Etha_Func(LambdaFor, TStar, previous_t, middle_t) +
            2.0 * RhoDX1 * SigmaDom * volX * Etha_Func(LambdaDom, TStar, previous_t, middle_t);
        VarFFx +=
            SigmaFor * SigmaFor * Psi_Func(LambdaFor, LambdaFor, TStar, previous_t, middle_t) +
            SigmaDom * SigmaDom * Psi_Func(LambdaDom, LambdaDom, TStar, previous_t, middle_t) -
            2.0 * RhoDF * SigmaDom * SigmaFor *
                Psi_Func(LambdaDom, LambdaFor, TStar, previous_t, middle_t);

        t = middle_t;

        if (iteration > iteration_limit)
        {
            err = "No Fix Point Solution for correlation mapping";
            return err;
        }
    } while ((fabs(*RhoDX - RhoDX1) > tolerance && fabs(*RhoFX - RhoFX1) > tolerance) ||
             iteration < 2);

    *RhoDX = RhoDX1;
    *RhoFX = RhoFX1;

    if (err)
        goto FREE_RETURN;

    err = FxBetaDLM_CheckOrMakePosMatrix(RhoDF, RhoDX, RhoFX);

    if (err)
        goto FREE_RETURN;

FREE_RETURN:
    return err;
}

/* Finds the Fix Point Solution */
Err FxBetaDLM_correl_mapping4(
    int              index,
    double           VarFFx,
    double           volX,
    double           RhoDS,
    double           RhoFS,
    FxBetaDLM_model* model,
    double*          RhoDX,
    double*          RhoFX)
{
    Err err = NULL;

    double RhoDF, RhoDX1, RhoFX1;
    double VarFFx_previous;
    double B, vol_dom, vol_for, vol_X;
    double t, next_t, previous_t;
    double a, b, c, delta;
    double tolerance = 10e-6;
    double iteration, iteration_limit;

    /* initialisation */
    iteration       = 0;
    iteration_limit = 20;

    if (index > 0)
    {
        t = model->dPWTime[index - 1];
    }
    else
    {
        t = 0.0;
    }

    previous_t      = t;
    next_t          = model->dPWTime[index];
    VarFFx_previous = VarFFx;

    RhoDF  = model->dCorrDomFor[index];
    RhoDX1 = 0;
    RhoFX1 = 0;

    do
    {
        *RhoDX = RhoDX1;
        *RhoFX = RhoFX1;

        iteration++;

        if ((1.0 + 2.0 * model->dC0 * VarFFx) <= 0.0)
        {
            err = "FxBetaDLM_correl_mapping: The cumulative variance is too high";
            return err;
        }
        B       = model->dB0 / (1.0 + 2.0 * model->dC0 * VarFFx);
        vol_dom = -model->dSigmaDom[index] * (1.0 - exp(-model->dLambdaDom * (model->dTStar - t))) /
                  model->dLambdaDom;
        vol_for = -model->dSigmaFor[index] * (1.0 - exp(-model->dLambdaFor * (model->dTStar - t))) /
                  model->dLambdaFor;

        /* add the coefficients */
        vol_dom *= -(B - 1.0);
        vol_for *= (B - 1.0);
        vol_X = volX * B;

        if (fabs(RhoDS) > 1.0E-08)
        {
            /* solve polynomial */
            a = 1.0;
            b = vol_dom + RhoDF * vol_for - RhoDS * RhoDS * vol_dom - RhoDS * RhoFS * vol_for;
            c = vol_dom * vol_dom + pow(RhoDF * vol_for, 2) + 2.0 * RhoDF * vol_dom * vol_for -
                pow(RhoDS * vol_dom, 2) + pow(RhoDS * vol_for, 2) - pow(RhoDS * vol_X, 2) -
                2.0 * RhoDS * RhoFS * vol_for * (vol_dom + RhoDF * vol_for);

            delta = b * b - a * c;

            if (delta < 0.0)
            {
                err = serror("No solution for correlation mapping, delta < 0 for t = %f", t);
                return err;
            }

            delta = sqrt(delta);

            if (RhoDS > 0.0)
            {
                /* we choose the positive solution */
                RhoDX1 = (-b + delta) / a;
            }
            else
            {
                RhoDX1 = (-b - delta) / a;
            }

            if (fabs(RhoDX1) > vol_X * 0.999)
            {
                err = "No solution for correlation mapping, correlation dom X not in [-1;1]";
                return err;
            }

            RhoFX1 =
                RhoFS / RhoDS * (vol_dom + RhoDF * vol_for + RhoDX1) - RhoDF * vol_dom - vol_for;
        }
        else
        {
            RhoDX1 = -vol_dom - RhoDF * vol_for;

            a = 1.0;
            b = vol_for + RhoDF * vol_dom - RhoFS * RhoFS * vol_for;
            c = pow(RhoDF * vol_dom, 2) + vol_for * vol_for +
                RhoFS * RhoFS * (vol_dom * vol_dom - vol_for * vol_for - vol_X * vol_X) +
                2.0 * RhoDF * vol_dom * vol_for;

            delta = b * b - a * c;

            if (delta < 0.0)
            {
                err = "No solution for correlation mapping, delta < 0";
                return err;
            }

            delta = sqrt(delta);

            if (RhoFS > 0.0)
            {
                RhoFX1 = (-b + delta) / a;
            }
            else
            {
                RhoFX1 = (-b - delta) / a;
            }

            if (fabs(RhoFX1) > vol_X * 0.999)
            {
                err = "No solution for correlation mapping, correlation For|X not in [-1;1]";
                return err;
            }
        }

        RhoDX1 /= vol_X;
        RhoFX1 /= vol_X;

        /* update var FFX */
        VarFFx = VarFFx_previous;

        VarFFx += volX * volX * (next_t - previous_t);
        VarFFx += -2.0 * RhoFX1 * model->dSigmaFor[index] * volX *
                      Etha_Func(model->dLambdaFor, model->dTStar, previous_t, next_t) +
                  2.0 * RhoDX1 * model->dSigmaDom[index] * volX *
                      Etha_Func(model->dLambdaDom, model->dTStar, previous_t, next_t);
        VarFFx +=
            model->dSigmaFor[index] * model->dSigmaFor[index] *
                Psi_Func(model->dLambdaFor, model->dLambdaFor, model->dTStar, previous_t, next_t) +
            model->dSigmaDom[index] * model->dSigmaDom[index] *
                Psi_Func(model->dLambdaDom, model->dLambdaDom, model->dTStar, previous_t, next_t) -
            2.0 * RhoDF * model->dSigmaDom[index] * model->dSigmaFor[index] *
                Psi_Func(model->dLambdaDom, model->dLambdaFor, model->dTStar, previous_t, next_t);

        t = next_t;

        if (iteration > iteration_limit)
        {
            err = "No Fix Point Solution for correlation mapping";
            return err;
        }
    } while (fabs(*RhoDX - RhoDX1) > tolerance && fabs(*RhoFX - RhoFX1) > tolerance);

    *RhoDX = RhoDX1;
    *RhoFX = RhoFX1;

    return err;
}

/* Solution with 2 iterations */
Err FxBetaDLM_correl_mapping3(
    int              index,
    double           VarFFx,
    double           volX,
    double           RhoDS,
    double           RhoFS,
    FxBetaDLM_model* model,
    double*          RhoDX,
    double*          RhoFX)
{
    Err err = NULL;

    double RhoDF, RhoDX1, RhoFX1, RhoDX2, RhoFX2;
    double TestRhoDX, TestRhoFX, VarFFxTest;
    double B, vol_dom, vol_for, vol_X;
    double t, next_t;
    double a, b, c, delta;
    int    step;

    VarFFxTest = VarFFx;

    if (index > 0)
    {
        t = model->dPWTime[index - 1];
    }
    else
    {
        t = 0.0;
    }

    next_t = model->dPWTime[index];

    RhoDF = model->dCorrDomFor[index];

    for (step = 0; step < 2; step++)
    {
        if ((1.0 + 2.0 * model->dC0 * VarFFx) <= 0.0)
        {
            err = serror(
                "FxBetaDLM_correl_mapping: The cumulative variance is too high for t = %f", t);
            return err;
        }
        B       = model->dB0 / (1.0 + 2.0 * model->dC0 * VarFFx);
        vol_dom = -model->dSigmaDom[index] * (1.0 - exp(-model->dLambdaDom * (model->dTStar - t))) /
                  model->dLambdaDom;
        vol_for = -model->dSigmaFor[index] * (1.0 - exp(-model->dLambdaFor * (model->dTStar - t))) /
                  model->dLambdaFor;

        /* add the coefficients */
        vol_dom *= -(B - 1.0);
        vol_for *= (B - 1.0);
        vol_X = volX * B;

        if (fabs(RhoDS) > 1.0E-08)
        {
            /* solve polynomial */
            a = 1.0;
            b = vol_dom + RhoDF * vol_for - RhoDS * RhoDS * vol_dom - RhoDS * RhoFS * vol_for;
            c = vol_dom * vol_dom + pow(RhoDF * vol_for, 2) + 2.0 * RhoDF * vol_dom * vol_for -
                pow(RhoDS * vol_dom, 2) + pow(RhoDS * vol_for, 2) - pow(RhoDS * vol_X, 2) -
                2.0 * RhoDS * RhoFS * vol_for * (vol_dom + RhoDF * vol_for);

            delta = b * b - a * c;

            if (delta < 0.0)
            {
                err = serror("No solution for correlation mapping, delta < 0 for t = %f", t);
                return err;
            }

            delta = sqrt(delta);

            if (RhoDS > 0.0)
            {
                /* we choose the positive solution */
                RhoDX1 = (-b + delta) / a;
            }
            else
            {
                RhoDX1 = (-b - delta) / a;
            }

            if (fabs(RhoDX1) > vol_X * 0.999)
            {
                err = serror(
                    "No solution for correlation mapping, correlation dom_X not in [-1;1] for t = "
                    "%f",
                    t);
                return err;
            }

            RhoFX1 =
                RhoFS / RhoDS * (vol_dom + RhoDF * vol_for + RhoDX1) - RhoDF * vol_dom - vol_for;
        }
        else
        {
            RhoDX1 = -vol_dom - RhoDF * vol_for;

            a = 1.0;
            b = vol_for + RhoDF * vol_dom - RhoFS * RhoFS * vol_for;
            c = pow(RhoDF * vol_dom, 2) + vol_for * vol_for +
                RhoFS * RhoFS * (vol_dom * vol_dom - vol_for * vol_for - vol_X * vol_X) +
                2.0 * RhoDF * vol_dom * vol_for;

            delta = b * b - a * c;

            if (delta < 0.0)
            {
                err = "No solution for correlation mapping, delta < 0";
                return err;
            }

            delta = sqrt(delta);

            if (RhoFS > 0.0)
            {
                RhoFX1 = (-b + delta) / a;
            }
            else
            {
                RhoFX1 = (-b - delta) / a;
            }

            if (fabs(RhoFX1) > vol_X * 0.999)
            {
                err = serror(
                    "No solution for correlation mapping, correlation for_X not in [-1;1] for t = "
                    "%f",
                    t);
                return err;
            }
        }

        RhoDX1 /= vol_X;
        RhoFX1 /= vol_X;

        /* update var FFX */

        if (step == 0)
        {
            VarFFx += volX * volX * (next_t - t);
            VarFFx += -2.0 * RhoFX1 * model->dSigmaFor[index] * volX *
                          Etha_Func(model->dLambdaFor, model->dTStar, t, next_t) +
                      2.0 * RhoDX1 * model->dSigmaDom[index] * volX *
                          Etha_Func(model->dLambdaDom, model->dTStar, t, next_t);
            VarFFx += model->dSigmaFor[index] * model->dSigmaFor[index] *
                          Psi_Func(model->dLambdaFor, model->dLambdaFor, model->dTStar, t, next_t) +
                      model->dSigmaDom[index] * model->dSigmaDom[index] *
                          Psi_Func(model->dLambdaDom, model->dLambdaDom, model->dTStar, t, next_t) -
                      2.0 * RhoDF * model->dSigmaDom[index] * model->dSigmaFor[index] *
                          Psi_Func(model->dLambdaDom, model->dLambdaFor, model->dTStar, t, next_t);

            t      = next_t;
            RhoDX2 = RhoDX1;
            RhoFX2 = RhoFX1;
        }
    }

    *RhoDX = 0.5 * (RhoDX1 + RhoDX2);
    *RhoFX = 0.5 * (RhoFX1 + RhoFX2);

    if (index > 0)
    {
        t = model->dPWTime[index - 1];
    }
    else
    {
        t = 0.0;
    }

    VarFFxTest += volX * volX * (next_t - t);
    VarFFxTest += -2.0 * *RhoFX * model->dSigmaFor[index] * volX *
                      Etha_Func(model->dLambdaFor, model->dTStar, t, next_t) +
                  2.0 * *RhoDX * model->dSigmaDom[index] * volX *
                      Etha_Func(model->dLambdaDom, model->dTStar, t, next_t);
    VarFFxTest += model->dSigmaFor[index] * model->dSigmaFor[index] *
                      Psi_Func(model->dLambdaFor, model->dLambdaFor, model->dTStar, t, next_t) +
                  model->dSigmaDom[index] * model->dSigmaDom[index] *
                      Psi_Func(model->dLambdaDom, model->dLambdaDom, model->dTStar, t, next_t) -
                  2.0 * RhoDF * model->dSigmaDom[index] * model->dSigmaFor[index] *
                      Psi_Func(model->dLambdaDom, model->dLambdaFor, model->dTStar, t, next_t);

    /* check */
    err = Get_Correl_FromModel(index, VarFFx, volX, *RhoDX, *RhoFX, model, &TestRhoDX, &TestRhoFX);

    return err;
}

/* Newton in two dimensions */
Err FxBetaDLM_correl_mapping2(
    int              index,
    double           VarFFx,
    double           VolX,
    double           RhoSD,
    double           RhoSF,
    FxBetaDLM_model* model,
    double*          RhoXD,
    double*          RhoXF)
{
    Err    err = NULL;
    double RhoSD1, RhoSD2, RhoSF1, RhoSF2, RhoXD1, RhoXD2, RhoXD3, RhoXF1, RhoXF2, RhoXF3;

    double sens11, sens12, sens21, sens22;
    double delta;

    double error, error1, error2;
    double coef1, coef2, coef3;
    int    nb_iter;

    /* First Guess */
    RhoXD1 = *RhoXD;
    RhoXF1 = *RhoXF;

    /* First Pricing */
    err = Get_Correl_FromModel(index, VarFFx, VolX, RhoXD1, RhoXF1, model, &RhoSD1, &RhoSF1);

    if (err)
        return err;

    error1 = RhoSD1 - RhoSD;
    error2 = RhoSF1 - RhoSF;
    error  = sqrt(pow(error1, 2) + pow(error2, 2));

    nb_iter = 0;

    while (nb_iter < 10 && error > 1.0E-03)
    {
        /* Calculate the gradient against RhoXD */
        RhoXD2 = RhoXD1 + DELTA_RHO;
        if (RhoXD2 > 0.999)
        {
            RhoXD2 = RhoXD1 - DELTA_RHO;
        }

        RhoXF2 = RhoSF1;

        err = Get_Correl_FromModel(index, VarFFx, VolX, RhoXD2, RhoXF2, model, &RhoSD2, &RhoSF2);

        if (err)
            return err;

        sens11 = (RhoSD2 - RhoSD1) / (RhoXD2 - RhoXD1);
        sens21 = (RhoSF2 - RhoSF1) / (RhoXD2 - RhoXD1);

        /* Calculate the gradient against RhoXF */
        RhoXD2 = RhoXD1;
        RhoXF2 = RhoSF1 + DELTA_RHO;

        if (RhoXF2 > 0.999)
        {
            RhoXF2 = RhoXF1 - DELTA_RHO;
        }

        err = Get_Correl_FromModel(index, VarFFx, VolX, RhoXD2, RhoXF2, model, &RhoSD2, &RhoSF2);

        if (err)
            return err;

        sens12 = (RhoSD2 - RhoSD1) / (RhoXF2 - RhoXF1);
        sens22 = (RhoSF2 - RhoSF1) / (RhoXF2 - RhoXF1);

        delta = sens11 * sens22 - sens12 * sens21;

        if (fabs(delta) < 1.0E-10)
        {
            err = "cannot calibrate the correlations";
            return err;
        }

        RhoXD3 = RhoXD1 - (sens22 * error1 - sens21 * error2) / delta;
        RhoXF3 = RhoXF1 - (-sens12 * error1 + sens11 * error2) / delta;

        if (fabs(RhoXD3) < 0 || fabs(RhoXF3) > 0.9999)
        {
            if (RhoXD3 < -0.9999)
            {
                coef1 = (RhoXD1 + 0.99) / ((sens22 * error1 - sens21 * error2) / delta);
            }
            else if (RhoXD3 > 0.9999)
            {
                coef1 = (RhoXD1 - 0.99) / ((sens22 * error1 - sens21 * error2) / delta);
            }

            if (RhoXF3 < -0.9999)
            {
                coef2 = (RhoXF1 + 0.99) / ((-sens12 * error1 + sens11 * error2) / delta);
            }
            else if (RhoXF3 > 0.9999)
            {
                coef2 = (RhoXF1 - 0.99) / ((-sens12 * error1 + sens11 * error2) / delta);
            }

            coef3 = min(coef1, coef2);

            RhoXD3 = RhoXD1 - coef3 * (sens22 * error1 - sens21 * error2) / delta;
            RhoXF3 = RhoXF1 - coef3 * (-sens12 * error1 + sens11 * error2) / delta;
        }

        RhoXD1 = RhoXD3;
        RhoXF1 = RhoXF3;

        err = Get_Correl_FromModel(index, VarFFx, VolX, RhoXD1, RhoXF1, model, &RhoSD1, &RhoSF1);

        if (err)
            return err;

        error1 = RhoSD1 - RhoSD;
        error2 = RhoSF1 - RhoSF;
        error  = sqrt(pow(error1, 2) + pow(error2, 2));

        nb_iter++;
    }

    *RhoXD = RhoXD1;
    *RhoXF = RhoXF1;

    return err;
}

Err Get_Correl_FromModel(
    int              index,
    double           VarFFx,
    double           VolX,
    double           RhoXD,
    double           RhoXF,
    FxBetaDLM_model* model,
    double*          RhoSD,
    double*          RhoSF)
{
    double t1, t2;
    double var_ffx;
    double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor, sig_domfx, sig_forfx;
    double dom_lam, for_lam, tstar;
    double RhoDF;
    double vol_X2;

    double B, C;
    double volS;
    double gamma_for;
    double gamma_dom;
    Err    err = NULL;

    dom_lam = model->dLambdaDom;
    for_lam = model->dLambdaFor;
    tstar   = model->dTStar;

    if (index == 0)
    {
        t1      = 0.0;
        var_ffx = 0.0;
    }
    else
    {
        t1      = model->dPWTime[index - 1];
        var_ffx = VarFFx;
    }
    t2 = model->dPWTime[index];

    sig_dom  = model->dSigmaDom[index];
    sig_for  = model->dSigmaFor[index];
    sig_dom2 = sig_dom * sig_dom;
    sig_for2 = sig_for * sig_for;

    RhoDF      = model->dCorrDomFor[index];
    sig_domfor = RhoDF * sig_dom * sig_for;
    sig_fx     = VolX;
    sig_fx2    = sig_fx * sig_fx;
    sig_domfx  = RhoXD * sig_dom * sig_fx;
    sig_forfx  = RhoXF * sig_for * sig_fx;

    B = model->dB0 / (1.0 + 2.0 * model->dC0 * var_ffx);
    C = model->dC0 / (1.0 + 2.0 * model->dC0 * var_ffx);

    gamma_dom = -sig_dom * (1.0 - exp(-dom_lam * (tstar - t2))) / dom_lam;
    gamma_for = -sig_for * (1.0 - exp(-for_lam * (tstar - t2))) / for_lam;

    vol_X2 = sig_fx2;
    vol_X2 += 2.0 * RhoXF * gamma_for * sig_fx - 2.0 * RhoXD * gamma_dom * sig_fx;
    vol_X2 += gamma_dom * gamma_dom + gamma_for * gamma_for - 2.0 * RhoDF * gamma_dom * gamma_for;

    if (vol_X2 < 0.0)
    {
        err = "Get_Correl_FromModel negative var of X";
        goto FREE_RETURN;
    }

    volS = B * B * VolX * VolX +
           (B - 1.0) * (B - 1.0) * (gamma_for * gamma_for + gamma_dom * gamma_dom);
    volS += 2.0 * (B - 1.0) *
            (B * (RhoXF * gamma_for - RhoXD * gamma_dom) * VolX -
             (B - 1.0) * RhoDF * gamma_dom * gamma_for);
    volS += 4.0 * C * C * var_ffx * vol_X2;

    if (volS < 0.0)
    {
        err = "Get_Correl_FromModel negative var of Spot";
        goto FREE_RETURN;
    }
    volS = sqrt(volS);

    *RhoSD = B * VolX * RhoXD + (B - 1.0) * (gamma_for * RhoDF - gamma_dom);
    *RhoSD /= volS;
    *RhoSF = B * VolX * RhoXF + (B - 1.0) * (gamma_for - gamma_dom * RhoDF);
    *RhoSF /= volS;

FREE_RETURN:

    return err;
}

Err FxBetaDLM_GetFirstGuessFromB0(
    double  B0,
    double  Tstar,
    double* exercise_opt,
    double* maturity_opt,
    long    nbrOpt,
    double* maturity_rates_corr,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    double* cal_vol_3f,
    double* fx_vol_curve)
{
    int    i, j;
    int    start_index, end_index;
    double t1, t2, start_mat, end_mat;
    double Integral, B0minus1;
    Err    err = NULL;

    if (!B0)
        err = "FxBetaDLM_GetFirstGuessFromB0:  B0 = 0 not supported ";

    if (err)
        goto FREE_RETURN;

    B0minus1 = B0 - 1.0;

    start_mat   = 0.0;
    start_index = Get_Index(start_mat, maturity_rates_corr, nbrMat);

    for (i = 0; i < nbrOpt; i++)
    {
        Integral  = 0.0;
        end_mat   = exercise_opt[i];
        end_index = Get_Index(end_mat, maturity_rates_corr, nbrMat);

        for (j = start_index; j < end_index + 1; j++)
        {
            if (j > start_index)
                t1 = maturity_rates_corr[j - 1];
            else
                t1 = start_mat; /* First part */

            if (j == end_index || start_index == end_index)
                t2 = end_mat; /* Last part */
            else
                t2 = maturity_rates_corr[j];

            Integral += B0minus1 * sig_curve_for[j] * sig_curve_for[j] *
                        Psi_Func(lda_for, lda_for, Tstar, t1, t2);
            Integral += B0minus1 * sig_curve_dom[j] * sig_curve_dom[j] *
                        Psi_Func(lda_dom, lda_dom, Tstar, t1, t2);
            Integral += 2.0 * correl_for_fx[j] * cal_vol_3f[i] * sig_curve_for[j] *
                        Etha_Func(lda_for, Tstar, t1, t2);
            Integral -= 2.0 * correl_dom_fx[j] * cal_vol_3f[i] * sig_curve_dom[j] *
                        Etha_Func(lda_dom, Tstar, t1, t2);
            Integral -= B0minus1 * 2.0 * correl_dom_for[j] * sig_curve_dom[j] * sig_curve_for[j] *
                        Psi_Func(lda_dom, lda_for, Tstar, t1, t2);
        }

        fx_vol_curve[i] =
            (cal_vol_3f[i] * cal_vol_3f[i] + B0minus1 * Integral / (end_mat - start_mat));
        fx_vol_curve[i] /= B0 * B0;
        fx_vol_curve[i] = sqrt(fx_vol_curve[i]);

        start_index = end_index;
    }

FREE_RETURN:
    return err;
}

Err FxBetaDLM_GetModelFromStr(FXBETADLM_STR str, FxBetaDLM_model* model)
{
    Err  err = NULL;
    long spot_date;
    int  i, index;

    spot_date = add_unit(str->today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    model->lToday  = str->today;
    model->cYcDom  = str->yc_dom;
    model->cYcFor  = str->yc_for;
    model->dB0     = str->b0;
    model->dC0     = str->c0;
    model->dAlpha  = str->alpha;
    model->dLambda = str->lambda;

    model->dSpotFx = str->spot_fx;
    model->dCashFx = str->spot_fx * swp_f_df(str->today, spot_date, str->yc_dom) /
                     swp_f_df(str->today, spot_date, str->yc_for);

    err = FxBetaDLM_get_beta_from_C0(str->tstar, str->c0, &model->dEquiBeta);

    if (err)
        goto FREE_RETURN;

    model->dTStar     = str->tstar;
    model->dInitTStar = str->tstar;
    model->lTStarDate = (long)(str->today + str->tstar * DAYS_IN_YEAR + 0.5);

    model->dTauDom = str->tau_dom;
    if (fabs(str->tau_dom) < 1.0E-08)
    {
        err = "FxBetaDLM: tau cannot be null";
        goto FREE_RETURN;
    }

    model->dLambdaDom = 1.0 / str->tau_dom;

    model->dTauFor = str->tau_for;
    if (fabs(str->tau_for) < 1.0E-08)
    {
        err = "FxBetaDLM: tau cannot be null";
        goto FREE_RETURN;
    }

    model->dLambdaFor = 1.0 / str->tau_for;

    /* Get the merged dates */

    model->iNbPWTime = str->nb_DLM_corr;

    /* Allocate memory */
    model->dPWTime         = calloc(str->nb_DLM_corr, sizeof(double));
    model->dSigmaDom       = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFor       = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFx        = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFxUp      = calloc(model->iNbPWTime, sizeof(double));
    model->dSigmaFxDown    = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFor     = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFx      = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrForFx      = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFxDown  = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrForFxDown  = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrDomFxInput = calloc(model->iNbPWTime, sizeof(double));
    model->dCorrForFxInput = calloc(model->iNbPWTime, sizeof(double));

    if (!model->dPWTime || !model->dSigmaDom || !model->dSigmaFor || !model->dSigmaFx ||
        !model->dSigmaFxUp || !model->dSigmaFxDown || !model->dCorrDomFor || !model->dCorrDomFx ||
        !model->dCorrForFx || !model->dCorrDomFxDown || !model->dCorrForFxDown ||
        !model->dCorrDomFxInput || !model->dCorrForFxInput)
    {
        err = "Memory allcation faillure in init_FxBetaDLM_model";
        goto FREE_RETURN;
    }

    /* Fill the vectors */
    memcpy(model->dPWTime, str->time_DLM_corr, str->nb_DLM_corr * sizeof(double));

    /* Dom */
    for (i = 0; i < model->iNbPWTime; i++)
    {
        index               = Get_Index(model->dPWTime[i], str->time_dom, str->nb_dom);
        model->dSigmaDom[i] = str->sig_dom[index];

        index               = Get_Index(model->dPWTime[i], str->time_for, str->nb_for);
        model->dSigmaFor[i] = str->sig_for[index];

        index              = Get_Index(model->dPWTime[i], str->time_fx, str->nb_fx);
        model->dSigmaFx[i] = str->sig_fx[index];
        if (fabs(model->dLambda) < TINY)
        {
            model->dSigmaFxUp[i] =
                str->sig_fx[index] * exp(model->dAlpha * sqrt(str->time_fx[index]));
            model->dSigmaFxDown[i] =
                str->sig_fx[index] * exp(-model->dAlpha * sqrt(str->time_fx[index]));
        }
        else
        {
            model->dSigmaFxUp[i] =
                str->sig_fx[index] *
                exp(model->dAlpha * sqrt(
                                        (1.0 - exp(-2.0 * model->dLambda * str->time_fx[index])) /
                                        (2.0 * model->dLambda)));
            model->dSigmaFxDown[i] =
                str->sig_fx[index] *
                exp(-model->dAlpha * sqrt(
                                         (1.0 - exp(-2.0 * model->dLambda * str->time_fx[index])) /
                                         (2.0 * model->dLambda)));
        }

        index = Get_Index(model->dPWTime[i], str->time_3F_corr, str->nb_3F_corr);
        model->dCorrDomFxInput[i] = str->dom_fx_3F_corr[index];
        model->dCorrForFxInput[i] = str->for_fx_3F_corr[index];
    }

    memcpy(model->dCorrForFx, str->for_fx_DLM_corr, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrDomFx, str->dom_fx_DLM_corr, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrForFxDown, str->for_fx_DLM_corr_down, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrDomFxDown, str->dom_fx_DLM_corr_down, model->iNbPWTime * sizeof(double));
    memcpy(model->dCorrDomFor, str->dom_for_DLM_corr, model->iNbPWTime * sizeof(double));

    /* No structure associated */
    model->str = NULL;

FREE_RETURN:

    if (err)
    {
        free_FxBetaDLM_model(model);
    }

    return err;
}

Err FxBetaDLM_get3FCorrelFromDLM(FXBETADLM_STR str)
{
    int    i, index;
    double t1, t2, middle_t;
    double vol_X2, var_ffx;
    double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor, sig_domfx, sig_forfx;
    double dom_lam, for_lam, tstar;
    double RhoDF, RhoXD, RhoXF;

    double B, C;
    double volS;
    double gamma_for;
    double gamma_dom;
    Err    err = NULL;

    str->nb_3F_corr      = str->nb_DLM_corr;
    str->dom_for_3F_corr = calloc(str->nb_DLM_corr, sizeof(double));
    str->time_3F_corr    = calloc(str->nb_DLM_corr, sizeof(double));
    str->dom_fx_3F_corr  = calloc(str->nb_DLM_corr, sizeof(double));
    ;
    str->for_fx_3F_corr = calloc(str->nb_DLM_corr, sizeof(double));
    ;

    if (!str->dom_for_3F_corr || !str->time_3F_corr || !str->dom_fx_3F_corr || !str->for_fx_3F_corr)
    {
        err = "Memory allocation faillure in FxBetaDLM_get3FCorrelFromDLM";
        goto FREE_RETURN;
    }

    memcpy(str->dom_for_3F_corr, str->dom_for_DLM_corr, str->nb_DLM_corr * sizeof(double));
    memcpy(str->time_3F_corr, str->time_DLM_corr, str->nb_DLM_corr * sizeof(double));

    dom_lam = str->lam_dom;
    for_lam = str->lam_for;
    tstar   = str->tstar;

    t1       = 0.0;
    t2       = str->time_DLM_corr[0];
    middle_t = 0.5 * (t1 + t2);
    var_ffx  = 0.0;

    for (i = 0; i < str->nb_DLM_corr; i++)
    {
        index   = Get_Index(t2, str->time_dom, str->nb_dom);
        sig_dom = str->sig_dom[index];

        index   = Get_Index(t2, str->time_for, str->nb_for);
        sig_for = str->sig_for[index];

        index  = Get_Index(t2, str->time_fx, str->nb_fx);
        sig_fx = str->sig_fx[index];

        sig_dom2 = sig_dom * sig_dom;
        sig_for2 = sig_for * sig_for;
        sig_fx2  = sig_fx * sig_fx;

        RhoDF = str->dom_for_3F_corr[i];
        RhoXD = str->dom_fx_DLM_corr[i];
        RhoXF = str->for_fx_DLM_corr[i];

        sig_domfor = RhoDF * sig_dom * sig_for;
        sig_domfx  = RhoXD * sig_dom * sig_fx;
        sig_forfx  = RhoXF * sig_for * sig_fx;

        var_ffx += sig_fx2 * (middle_t - t1);
        var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, middle_t) +
                   2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, middle_t);
        var_ffx += sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, middle_t) +
                   sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, middle_t) -
                   2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, middle_t);

        B = str->b0 / (1.0 + 2.0 * str->c0 * var_ffx);
        C = str->c0 / (1.0 + 2.0 * str->c0 * var_ffx);

        gamma_dom = -sig_dom * (1.0 - exp(-dom_lam * (tstar - middle_t))) / dom_lam;
        gamma_for = -sig_for * (1.0 - exp(-for_lam * (tstar - middle_t))) / for_lam;

        vol_X2 = sig_fx2;
        vol_X2 += 2.0 * RhoXF * gamma_for * sig_fx - 2.0 * RhoXD * gamma_dom * sig_fx;
        vol_X2 +=
            gamma_dom * gamma_dom + gamma_for * gamma_for - 2.0 * RhoDF * gamma_dom * gamma_for;

        if (vol_X2 < 0.0)
        {
            err = "Get_Correl_FromModel negative var of X";
            goto FREE_RETURN;
        }

        volS = B * B * sig_fx * sig_fx +
               (B - 1.0) * (B - 1.0) * (gamma_for * gamma_for + gamma_dom * gamma_dom);
        volS += 2.0 * (B - 1.0) *
                (B * (RhoXF * gamma_for - RhoXD * gamma_dom) * sig_fx -
                 (B - 1.0) * RhoDF * gamma_dom * gamma_for);
        volS += 4.0 * C * C * var_ffx * vol_X2;

        if (volS < 0.0)
        {
            err = "Get_Correl_FromModel negative var of Spot";
            goto FREE_RETURN;
        }

        volS = sqrt(volS);

        str->dom_fx_3F_corr[i] = B * sig_fx * RhoXD + (B - 1.0) * (gamma_for * RhoDF - gamma_dom);
        str->dom_fx_3F_corr[i] /= volS;
        str->for_fx_3F_corr[i] = B * sig_fx * RhoXF + (B - 1.0) * (gamma_for - gamma_dom * RhoDF);
        str->for_fx_3F_corr[i] /= volS;

        var_ffx += sig_fx2 * (t2 - middle_t);
        var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, middle_t, t2) +
                   2.0 * sig_domfx * Etha_Func(dom_lam, tstar, middle_t, t2);
        var_ffx += sig_for2 * Psi_Func(for_lam, for_lam, tstar, middle_t, t2) +
                   sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, middle_t, t2) -
                   2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, middle_t, t2);

        t1       = t2;
        t2       = str->time_DLM_corr[i + 1];
        middle_t = 0.5 * (t1 + t2);
    }

FREE_RETURN:

    return err;
}

Err FxBetaDLM_getDLMCorrelFrom3F(FXBETADLM_STR str, FxBetaDLM_OptNumerParams* NumParams)
{
    int    n, i, j, NbTime;
    int    indexSigFx, indexSigDom, indexSigFor, indexCorr;
    double t, dt, t1, t2;
    double varFFx;
    Err    err = NULL;

    str->nb_DLM_corr   = str->nb_fx;
    str->time_DLM_corr = calloc(str->nb_fx, sizeof(double));

    if (!str->time_DLM_corr)
    {
        err = "Memory allcation faillure in init_FxBetaDLM_model";
        goto FREE_RETURN;
    }

    memcpy(str->time_DLM_corr, str->time_fx, str->nb_fx * sizeof(double));

    num_f_concat_vector(&str->nb_DLM_corr, &str->time_DLM_corr, str->nb_dom, str->time_dom);
    num_f_concat_vector(&str->nb_DLM_corr, &str->time_DLM_corr, str->nb_for, str->time_for);
    num_f_concat_vector(&str->nb_DLM_corr, &str->time_DLM_corr, str->nb_3F_corr, str->time_3F_corr);
    num_f_sort_vector(str->nb_DLM_corr, str->time_DLM_corr);
    num_f_unique_vector(&str->nb_DLM_corr, str->time_DLM_corr);

    /* add the intermediary dates for the correlation TS */
    if (str->time_DLM_corr[str->nb_DLM_corr - 1] < NumParams->dMaxTimeCorrel)
    {
        num_f_add_number(&str->nb_DLM_corr, &str->time_DLM_corr, NumParams->dMaxTimeCorrel);
    }

    NbTime = str->nb_DLM_corr;
    for (i = 0; i < str->nb_DLM_corr; i++)
    {
        if (!i)
        {
            if (str->time_DLM_corr[i] > NumParams->dMinTimeSpaceCorrel)
            {
                n  = (int)((str->time_DLM_corr[i]) / NumParams->dMinTimeSpaceCorrel);
                dt = (str->time_DLM_corr[i]) / (n + 1);
                t  = dt;
                for (j = 0; j < n; j++)
                {
                    num_f_add_number(
                        &NbTime,
                        &str->time_DLM_corr,
                        (double)(((long)(t * DAYS_IN_YEAR)) * YEARS_IN_DAY));
                    t += dt;
                }
            }
        }
        else
        {
            if (str->time_DLM_corr[i] - str->time_DLM_corr[i - 1] > NumParams->dMinTimeSpaceCorrel)
            {
                n = (int)((str->time_DLM_corr[i] - str->time_DLM_corr[i - 1]) / NumParams->dMinTimeSpaceCorrel);
                dt = (str->time_DLM_corr[i] - str->time_DLM_corr[i - 1]) / (n + 1);
                t  = str->time_DLM_corr[i - 1] + dt;
                for (j = 0; j < n; j++)
                {
                    num_f_add_number(
                        &NbTime,
                        &str->time_DLM_corr,
                        (double)(((long)(t * DAYS_IN_YEAR)) * YEARS_IN_DAY));
                    t += dt;
                }
            }
        }
    }
    str->nb_DLM_corr = NbTime;
    num_f_sort_vector(str->nb_DLM_corr, str->time_DLM_corr);
    num_f_unique_vector(&str->nb_DLM_corr, str->time_DLM_corr);

    /* We need now to convert the correlations */
    varFFx = 0.0;
    t1     = 0.0;

    str->dom_for_DLM_corr     = calloc(str->nb_DLM_corr, sizeof(double));
    str->dom_fx_DLM_corr      = calloc(str->nb_DLM_corr, sizeof(double));
    str->for_fx_DLM_corr      = calloc(str->nb_DLM_corr, sizeof(double));
    str->dom_fx_DLM_corr_down = calloc(str->nb_DLM_corr, sizeof(double));
    str->for_fx_DLM_corr_down = calloc(str->nb_DLM_corr, sizeof(double));

    if (!str->dom_for_DLM_corr || !str->dom_fx_DLM_corr || !str->for_fx_DLM_corr ||
        !str->dom_fx_DLM_corr_down || !str->for_fx_DLM_corr_down)
    {
        err = "Memory allocation faillure in FxBetaDLM_getDLMCorrelFrom3F";
        goto FREE_RETURN;
    }

    for (i = 0; i < str->nb_DLM_corr; i++)
    {
        t2 = str->time_DLM_corr[i];

        indexSigFx  = Get_Index(t2, str->time_fx, str->nb_fx);
        indexSigDom = Get_Index(t2, str->time_dom, str->nb_dom);
        indexSigFor = Get_Index(t2, str->time_for, str->nb_for);
        indexCorr   = Get_Index(t2, str->time_3F_corr, str->nb_3F_corr);

        str->dom_for_DLM_corr[i] = str->dom_for_3F_corr[indexCorr];

        /* transform correlations */
        err = FxBetaDLM_correl_mapping(
            str->sig_fx[indexSigFx],
            str->b0,
            str->c0,
            str->tstar,
            str->sig_dom[indexSigDom],
            str->lam_dom,
            str->sig_for[indexSigFor],
            str->lam_for,
            str->dom_for_3F_corr[indexCorr],
            str->dom_fx_3F_corr[indexCorr],
            str->for_fx_3F_corr[indexCorr],
            t1,
            t2,
            varFFx,
            &(str->dom_fx_DLM_corr[i]),
            &(str->for_fx_DLM_corr[i]));
        if (err)
            goto FREE_RETURN;

        str->dom_fx_DLM_corr_down[i] = str->dom_fx_DLM_corr[i];
        str->for_fx_DLM_corr_down[i] = str->for_fx_DLM_corr[i];

        /* update var ffx */
        varFFx += str->sig_fx[indexSigFx] * str->sig_fx[indexSigFx] * (t2 - t1);
        varFFx += -2.0 * str->for_fx_DLM_corr[i] * str->sig_for[indexSigFor] *
                      str->sig_fx[indexSigFx] * Etha_Func(str->lam_for, str->tstar, t1, t2) +
                  2.0 * str->dom_fx_DLM_corr[i] * str->sig_dom[indexSigDom] *
                      str->sig_fx[indexSigFx] * Etha_Func(str->lam_dom, str->tstar, t1, t2);
        varFFx += str->sig_for[indexSigFor] * str->sig_for[indexSigFor] *
                      Psi_Func(str->lam_for, str->lam_for, str->tstar, t1, t2) +
                  str->sig_dom[indexSigDom] * str->sig_dom[indexSigDom] *
                      Psi_Func(str->lam_dom, str->lam_dom, str->tstar, t1, t2) -
                  2.0 * str->dom_for_3F_corr[indexCorr] * str->sig_dom[indexSigDom] *
                      str->sig_for[indexSigFor] *
                      Psi_Func(str->lam_dom, str->lam_for, str->tstar, t1, t2);

        t1 = t2;
    }

FREE_RETURN:
    return err;
}

Err FxBetaDLM_CheckOrMakePosMatrix2(double RhoDF, double* RhoDX, double* RhoFX)
{
    Err      err = NULL;
    double   check;
    double   vp = 1.0;
    int      i, j;
    double** correl_matrix = NULL;
    double   eigen_val[3], eigen_val_mult = 0.05, eigen_val_step = 0.05;
    double** eigen_vector    = NULL;
    double   det_tol         = 1e-4;
    int      iteration_limit = 10, iteration = 0;

    check = 1 - RhoDF * RhoDF - *RhoDX * *RhoDX - *RhoFX * *RhoFX + 2.0 * RhoDF * *RhoDX * *RhoFX;

    if (check > det_tol || fabs(RhoDF - 1.0) < 1.0e-6)
        return err;
    else
    {
        correl_matrix = dmatrix(0, 2, 0, 2);
        eigen_vector  = dmatrix(0, 2, 0, 2);

        if (!correl_matrix || !eigen_vector)
        {
            err = "Memory allocation error in FxBetaDLM_CheckOrMakePosMatrix";
            goto FREE_RETURN;
        }

        /*	Construct correliance Matrix */
        correl_matrix[0][1] = RhoDF;
        correl_matrix[0][2] = *RhoDX;
        correl_matrix[1][2] = *RhoFX;

        for (i = 0; i < 3; i++)
        {
            correl_matrix[i][i] = 1.0;
            for (j = i + 1; j < 3; j++)
            {
                correl_matrix[j][i] = correl_matrix[i][j];
            }
        }

        /*	Diagonalise correliance Matrix */
        err = diagonalise_symmetric_matrix(correl_matrix, 3, eigen_val, eigen_vector);

        if (err)
            goto FREE_RETURN;

        do
        {
            eigen_val[2] = eigen_val_mult * eigen_val[1];

            *RhoDX = eigen_val[0] * eigen_vector[2][0] * eigen_vector[0][0] +
                     eigen_val[1] * eigen_vector[2][1] * eigen_vector[0][1] +
                     eigen_val[2] * eigen_vector[2][2] * eigen_vector[0][2];
            *RhoFX = eigen_val[0] * eigen_vector[2][0] * eigen_vector[1][0] +
                     eigen_val[1] * eigen_vector[2][1] * eigen_vector[1][1] +
                     eigen_val[2] * eigen_vector[2][2] * eigen_vector[1][2];

            check = 1 - RhoDF * RhoDF - *RhoDX * *RhoDX - *RhoFX * *RhoFX +
                    2.0 * RhoDF * *RhoDX * *RhoFX;

            eigen_val_mult += eigen_val_step;
            iteration++;

        } while (check < 0.0 && iteration < iteration_limit);
    }

    if (fabs(*RhoDX) > 1.0 || fabs(*RhoFX) > 1.0 || check < det_tol)
        err = "Matrix not positive definite";

FREE_RETURN:

    if (correl_matrix)
        free_dmatrix(correl_matrix, 0, 2, 0, 2);
    if (eigen_vector)
        free_dmatrix(eigen_vector, 0, 2, 0, 2);

    return err;
}

Err FxBetaDLM_CheckOrMakePosMatrix(double RhoDF, double* RhoDX, double* RhoFX)
{
    Err    err = NULL;
    double check;
    double delta;
    double Epsilon1, Epsilon2;
    double shift;
    double maxShift = 0.0001;
    int    sign     = 1;

    check = 1 - RhoDF * RhoDF - *RhoDX * *RhoDX - *RhoFX * *RhoFX + 2.0 * RhoDF * *RhoDX * *RhoFX;

    if (check > 0.0001 || fabs(RhoDF - 1.0) < 1.0e-6)
        return err;
    else
    {
        if (*RhoDX * *RhoFX > 0.0)
        {
            if (*RhoDX > 0.0)
                sign = -1;

            if (fabs(*RhoDX + *RhoFX) > 2 * maxShift)
                shift = sign * maxShift;
            else
                shift = sign * 0.5 * fabs(*RhoDX + *RhoFX);

            delta = (*RhoDX + *RhoFX) * (*RhoDX + *RhoFX) - 2.0 * check / (RhoDF - 1.0);

            if (delta < 0.0)
            {
                err = "Matrix not positive definite";
                return err;
            }

            delta = sqrt(delta);

            Epsilon1 = 0.5 * (-*RhoDX - *RhoFX + delta);
            Epsilon2 = 0.5 * (-*RhoDX - *RhoFX - delta);

            if (fabs(Epsilon1) < fabs(Epsilon2))
            {
                *RhoDX += Epsilon1 + shift;
                *RhoFX += Epsilon1 + shift;
            }
            else
            {
                *RhoDX += Epsilon2 + shift;
                *RhoFX += Epsilon2 + shift;
            }
        }
        else
        {
            if (*RhoDX < 0.0)
                sign = -1;

            if (fabs(*RhoDX - *RhoFX) > 2 * maxShift)
                shift = maxShift;
            else
                shift = 0.5 * fabs(*RhoDX - *RhoFX);

            delta = (*RhoDX - *RhoFX) * (*RhoDX - *RhoFX) + 2.0 * check / (1.0 + RhoDF);

            if (delta < 0.0)
            {
                err = "Matrix not positive definite";
                return err;
            }

            delta = sqrt(delta);

            Epsilon1 = 0.5 * (*RhoDX - *RhoFX + delta);
            Epsilon2 = 0.5 * (*RhoDX - *RhoFX - delta);

            if (fabs(Epsilon1) < fabs(Epsilon2))
            {
                *RhoDX -= sign * (fabs(Epsilon1) + shift);
                *RhoFX += sign * (fabs(Epsilon1) + shift);
            }
            else
            {
                *RhoDX -= sign * (fabs(Epsilon2) + shift);
                *RhoFX += sign * (fabs(Epsilon2) + shift);
            }
        }
    }

    if (fabs(*RhoDX) > 1.0 || fabs(*RhoFX) > 1.0)
        err = "Matrix not positive definite";
    return err;
}

Err copy_FxBetaDLM_model(FxBetaDLM_model* model_source, FxBetaDLM_model* model_dest)
{
    Err err = NULL;

    model_dest->lToday     = model_source->lToday;
    model_dest->cYcDom     = model_source->cYcDom;
    model_dest->cYcFor     = model_source->cYcFor;
    model_dest->dB0        = model_source->dB0;
    model_dest->dC0        = model_source->dC0;
    model_dest->dAlpha     = model_source->dAlpha;
    model_dest->dLambda    = model_source->dLambda;
    model_dest->dSpotFx    = model_source->dSpotFx;
    model_dest->dCashFx    = model_source->dCashFx;
    model_dest->dEquiBeta  = model_source->dEquiBeta;
    model_dest->dTStar     = model_source->dTStar;
    model_dest->dInitTStar = model_source->dInitTStar;
    model_dest->lTStarDate = model_source->lTStarDate;
    model_dest->dTauDom    = model_source->dTauDom;
    model_dest->dLambdaDom = model_source->dLambdaDom;
    model_dest->dTauFor    = model_source->dTauFor;
    model_dest->dLambdaFor = model_source->dLambdaFor;
    model_dest->iNbPWTime  = model_source->iNbPWTime;

    /* Allocate memory */
    model_dest->dPWTime         = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dSigmaDom       = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dSigmaFor       = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dSigmaFx3F      = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dSigmaFx        = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dSigmaFxUp      = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dSigmaFxDown    = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrDomFor     = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrDomFx      = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrForFx      = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrDomFxDown  = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrForFxDown  = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrDomFxInput = calloc(model_dest->iNbPWTime, sizeof(double));
    model_dest->dCorrForFxInput = calloc(model_dest->iNbPWTime, sizeof(double));

    if (!model_dest->dPWTime || !model_dest->dSigmaDom || !model_dest->dSigmaFor ||
        !model_dest->dSigmaFx3F || !model_dest->dSigmaFx || !model_dest->dSigmaFxUp ||
        !model_dest->dSigmaFxDown || !model_dest->dCorrDomFor || !model_dest->dCorrDomFx ||
        !model_dest->dCorrForFx || !model_dest->dCorrDomFxDown || !model_dest->dCorrForFxDown ||
        !model_dest->dCorrDomFxInput || !model_dest->dCorrForFxInput)
    {
        err = "Memory allcation faillure in init_FxBetaDLM_model_dest";
        goto FREE_RETURN;
    }

    /* Copy the vectors */
    memcpy(model_dest->dPWTime, model_source->dPWTime, model_dest->iNbPWTime * sizeof(double));
    memcpy(model_dest->dSigmaDom, model_source->dSigmaDom, model_dest->iNbPWTime * sizeof(double));
    memcpy(model_dest->dSigmaFor, model_source->dSigmaFor, model_dest->iNbPWTime * sizeof(double));
    memcpy(model_dest->dSigmaFx, model_source->dSigmaFx, model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dSigmaFx3F, model_source->dSigmaFx3F, model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dSigmaFxUp, model_source->dSigmaFxUp, model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dSigmaFxDown,
        model_source->dSigmaFxDown,
        model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrDomFor, model_source->dCorrDomFor, model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrDomFx, model_source->dCorrDomFx, model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrForFx, model_source->dCorrForFx, model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrForFxDown,
        model_source->dCorrForFxDown,
        model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrDomFxDown,
        model_source->dCorrDomFxDown,
        model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrForFxInput,
        model_source->dCorrForFxInput,
        model_dest->iNbPWTime * sizeof(double));
    memcpy(
        model_dest->dCorrDomFxInput,
        model_source->dCorrDomFxInput,
        model_dest->iNbPWTime * sizeof(double));

FREE_RETURN:
    return (err);
}