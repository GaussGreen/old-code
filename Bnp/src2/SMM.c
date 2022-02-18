// Swap Market Model

#include "CopulaGaussian.h"
#include "Product.h"
#include "RandomGen.h"
#include "math.h"
#include "opfnctns.h"
#include "opsabrgenericinterp.h"
#include "srt_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_vol.h"
#undef SIGN
#include "SMM.h"
#include "nag.h"
#include "nag_stdlib.h"
#include "nagf02.h"

static double level_cash(double x, int n, double* cvg)
{
    int    i, j;
    double tmp, tmpi, lvlc = 0.0;

    for (i = 0; i < n; i++)
    {
        tmp = 1.0 + cvg[i] * x;
        for (j = 1, tmpi = tmp; j < i + 1; j++)
            tmpi *= tmp;
        lvlc += cvg[i] / tmpi;
    }
    return lvlc;
}

static Err gf_calc_deriv(double x, smm_swap* swp, double shift, double* g, double* f)
{
    Err    err = NULL;
    int    i, j;
    double tmp, tmpi, lvlc0 = 0.0, lvlc1 = 0.0, lvlc2 = 0.0, p[2];

    // g derivs:

    for (i = 0; i < swp->nfix; i++)
    {
        tmp = 1.0 + swp->fix_cvg[i] * x;
        for (j = 1, tmpi = tmp; j < i + 1; j++)
            tmpi *= tmp;
        lvlc0 += swp->fix_cvg[i] / tmpi;
        lvlc1 += swp->fix_cvg[i] * swp->fix_cvg[i] * (i + 1) / (tmpi * tmp);
        lvlc2 += swp->fix_cvg[i] * swp->fix_cvg[i] * swp->fix_cvg[i] * (i + 1) * (i + 2) /
                 (tmpi * tmp * tmp);
    }
    g[0] = 1.0 / lvlc0;
    g[1] = g[0] * g[0] * lvlc1;
    g[2] = 2.0 * g[0] * g[1] * lvlc1 - g[0] * g[0] * lvlc2;

    // test:
    g[0] = 1.0 / swp->lvlc0;
    g[1] = g[2] = 0.0;

    // f derivs:

    err = BMM_Option_Price_From_States(
        x,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_CALL,
        &f[0]);
    if (err)
        return err;

    err = BMM_Option_Price_From_States(
        x - shift,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_CALL,
        &p[0]);
    if (err)
        return err;

    err = BMM_Option_Price_From_States(
        x + shift,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_CALL,
        &p[1]);
    if (err)
        return err;

    f[1] = (p[1] - p[0]) / (2.0 * shift);

    return NULL;
}

// Compute the Cumulative of a BMM distribution !!!
// Recurcive version with control of interpolation

static Err smm_BMM_linterp_cumulative(
    smm_swap* swp,

    double dMaxError,
    double dNbStdLeft,
    double dNbStdRight,

    long*    lNumPoints,
    double** dPoints,
    double** dCumulative)
{
    Err    err = NULL;
    double dShift, dShiftCumul;
    double dMinPoint, dMaxPoint;
    double dNbStdDevLeft, dNbStdDevRight, dLogStd, dLogCvx;
    double dLastPoint, dNewPoint, dInterpPoint;
    double dLastCumul, dNewCumul, dInterpCumul, dRealCumul;
    double dLastF[2], dNewF[2], dRealF[2];
    double dLastG[3], dNewG[3], dRealG[3];
    long   i, lAllocTotPoints;
    double dError;

    // First Guess
    srt_f_optbmmvolfromstates(
        swp->cms_fwd,  //		 swp->swap0,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_LOGNORMAL,
        &dLogStd);

    dLogStd *= sqrt(swp->fxg_time);
    dLogCvx = -0.5 * dLogStd * dLogStd;

    dNbStdDevLeft  = min(max(dNbStdLeft * sqrt(swp->fxg_time), 5.0), 10);
    dNbStdDevRight = dNbStdRight + 1.0;

    dMinPoint = /*swp->swap0 */ swp->cms_fwd * exp(dLogCvx - dNbStdDevLeft * dLogStd);

    do
    {
        dNbStdDevRight -= 1.0;
        dMaxPoint = /*swp->swap0 */ swp->cms_fwd * exp(dLogCvx + dNbStdDevRight * dLogStd);

        err = BMM_Option_Price_From_States(
            dMaxPoint,
            swp->fxg_time,
            swp->beta,
            swp->fwd1,
            swp->fwd2,
            swp->sigbeta1,
            swp->sigbeta2,
            swp->pi,
            SRT_CALL,
            &dNewF[0]);
        if (err)
            return err;

    } while (dNewF[0] < 1e-16);

    lAllocTotPoints = *lNumPoints;

    (*dPoints)[0] = dMinPoint;

    dShiftCumul = min(0.99 * (*dPoints)[0], /*swp->swap0 */ swp->cms_fwd / 100.0);

    err = gf_calc_deriv(0.0, swp, dShiftCumul, dLastG, dLastF);
    if (err)
        return err;

    dLastF[1] = -1.0;

    err = gf_calc_deriv((*dPoints)[0], swp, dShiftCumul, dNewG, dNewF);
    if (err)
        return err;

    (*dCumulative)[0] =
        (dNewG[0] * dNewF[1] - dLastG[0] * dLastF[1] - dNewG[1] * dNewF[0] + dLastG[1] * dLastF[0] +
         0.5 * (dNewG[2] * dNewF[0] + dLastG[2] * dLastF[0]) * (*dPoints)[0]) *
        swp->lvlc0;

    // Do the iterations
    dShift = exp((dNbStdDevLeft + dNbStdDevRight) * dLogStd / (*lNumPoints - 1));

    memcpy(dLastF, dNewF, 2 * sizeof(double));
    memcpy(dLastG, dNewG, 3 * sizeof(double));
    dLastCumul = (*dCumulative)[0];
    dLastPoint = (*dPoints)[0];
    dNewPoint  = dMinPoint * dShift;
    i          = 1;

    while (dNewPoint < dMaxPoint)
    {
        // Check the interpolation error
        dError = 1000.0;

        // dShiftCumul = 0.5 * (dNewPoint - dLastPoint);

        err = gf_calc_deriv(dNewPoint, swp, dShiftCumul, dNewG, dNewF);
        if (err)
            return err;

        dNewCumul = dLastCumul + (dNewG[0] * dNewF[1] - dLastG[0] * dLastF[1] -
                                  dNewG[1] * dNewF[0] + dLastG[1] * dLastF[0] +
                                  0.5 * (dNewG[2] * dNewF[0] + dLastG[2] * dLastF[0]) *
                                      (dNewPoint - dLastPoint)) *
                                     swp->lvlc0;

        dInterpPoint = dNewPoint;
        dRealCumul   = dNewCumul;
        memcpy(dRealF, dNewF, 2 * sizeof(double));
        memcpy(dRealG, dNewG, 3 * sizeof(double));

        while (dError > dMaxError)
        {
            dNewPoint = dInterpPoint;
            dNewCumul = dRealCumul;
            memcpy(dNewF, dRealF, 2 * sizeof(double));
            memcpy(dNewG, dRealG, 3 * sizeof(double));

            dInterpPoint = 0.5 * (dLastPoint + dNewPoint);

            // dShiftCumul = 0.5 * (dInterpPoint - dLastPoint);

            err = gf_calc_deriv(dInterpPoint, swp, dShiftCumul, dRealG, dRealF);
            if (err)
                return err;

            dRealCumul = dLastCumul + (dRealG[0] * dRealF[1] - dLastG[0] * dLastF[1] -
                                       dRealG[1] * dRealF[0] + dLastG[1] * dLastF[0] +
                                       0.5 * (dRealG[2] * dRealF[0] + dLastG[2] * dLastF[0]) *
                                           (dInterpPoint - dLastPoint)) *
                                          swp->lvlc0;

            dInterpCumul = dLastCumul + (dNewCumul - dLastCumul) / (dNewPoint - dLastPoint) *
                                            (dInterpPoint - dLastPoint);
            dError = fabs(dInterpCumul - dRealCumul);
        }

        (*dPoints)[i]     = dNewPoint;
        (*dCumulative)[i] = dNewCumul;
        dLastPoint        = dNewPoint;
        dLastCumul        = dNewCumul;
        memcpy(dLastF, dNewF, 2 * sizeof(double));
        memcpy(dLastG, dNewG, 3 * sizeof(double));
        dNewPoint *= dShift;
        i++;

        if (i >= lAllocTotPoints)
        {
            lAllocTotPoints += *lNumPoints;
            *dPoints     = realloc(*dPoints, lAllocTotPoints * sizeof(double));
            *dCumulative = realloc(*dCumulative, lAllocTotPoints * sizeof(double));
        }
    }

    *lNumPoints = i;

    return err;
}

static Err smm_swap_density(smm_swap* swp, double x, double shift, double* density)
{
    Err    err = NULL;
    double p[3];

    err = BMM_Option_Price_From_States(
        x - shift,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_CALL,
        &p[0]);
    if (err)
        return err;

    err = BMM_Option_Price_From_States(
        x,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_CALL,
        &p[1]);
    if (err)
        return err;

    err = BMM_Option_Price_From_States(
        x + shift,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_CALL,
        &p[2]);
    if (err)
        return err;

    *density = (p[0] + p[2] - 2.0 * p[1]) / (shift * shift) *
               (swp->lvlc0 / level_cash(x, swp->nfix, swp->fix_cvg));

    return NULL;
}

// Compute the Cumulative of a BMM distribution !!!
// Recurcive version with control of interpolation

static Err smm_BMM_linterp_cumulative2(
    smm_swap* swp,

    double dMaxError,
    double dNbStdLeft,
    double dNbStdRight,

    long*    lNumPoints,
    double** dPoints,
    double** dCumulative)
{
    Err    err = NULL;
    double dShift, dShiftCumul;
    double dMinPoint, dMaxPoint;
    double dNbStdDevLeft, dNbStdDevRight, dLogStd, dLogCvx;
    double dLastPoint, dNewPoint, dInterpPoint, dLastCumul, dNewCumul, dInterpCumul, dRealCumul,
        dLastDensity, dNewDensity, dRealDensity;
    long   i, lAllocTotPoints;
    double dError;

    // First Guess
    srt_f_optbmmvolfromstates(
        swp->swap0,
        swp->fxg_time,
        swp->beta,
        swp->fwd1,
        swp->fwd2,
        swp->sigbeta1,
        swp->sigbeta2,
        swp->pi,
        SRT_LOGNORMAL,
        &dLogStd);

    dLogStd *= sqrt(swp->fxg_time);
    dLogCvx = -0.5 * dLogStd * dLogStd;

    dNbStdDevLeft  = min(max(dNbStdLeft * sqrt(swp->fxg_time), 5.0), 10);
    dNbStdDevRight = dNbStdRight;

    dMinPoint = swp->swap0 * exp(dLogCvx - dNbStdDevLeft * dLogStd);
    dMaxPoint = swp->swap0 * exp(dLogCvx + dNbStdDevRight * dLogStd);

    lAllocTotPoints = *lNumPoints;

    (*dPoints)[0] = dMinPoint;

    dShiftCumul = min(0.99 * (*dPoints)[0], swp->swap0 / 100.0);

    err = smm_swap_density(swp, (*dPoints)[0], dShiftCumul, &dNewDensity);
    if (err)
        return err;

    (*dCumulative)[0] = 0.5 * dNewDensity * (*dPoints)[0];

    // Do the iterations
    dShift = exp((dNbStdDevLeft + dNbStdDevRight) * dLogStd / (*lNumPoints - 1));

    dLastDensity = dNewDensity;
    dLastCumul   = (*dCumulative)[0];
    dLastPoint   = (*dPoints)[0];
    dNewPoint    = dMinPoint * dShift;
    i            = 1;

    while (dNewPoint < dMaxPoint)
    {
        // Check the interpolation error
        dError = 1000.0;

        // dShiftCumul = 0.5 * (dNewPoint - dLastPoint);

        err = smm_swap_density(swp, dNewPoint, dShiftCumul, &dNewDensity);
        if (err)
            return err;

        dNewCumul = dLastCumul + 0.5 * (dLastDensity + dNewDensity) * (dNewPoint - dLastPoint);

        dInterpPoint = dNewPoint;
        dRealCumul   = dNewCumul;
        dRealDensity = dNewDensity;

        while (dError > dMaxError)
        {
            dNewPoint   = dInterpPoint;
            dNewCumul   = dRealCumul;
            dNewDensity = dRealDensity;

            dInterpPoint = 0.5 * (dLastPoint + dNewPoint);

            // dShiftCumul = 0.5 * (dInterpPoint - dLastPoint);

            err = smm_swap_density(swp, dInterpPoint, dShiftCumul, &dRealDensity);
            if (err)
                return err;

            dRealCumul =
                dLastCumul + 0.5 * (dLastDensity + dRealDensity) * (dInterpPoint - dLastPoint);

            dInterpCumul = dLastCumul + (dNewCumul - dLastCumul) / (dNewPoint - dLastPoint) *
                                            (dInterpPoint - dLastPoint);
            dError = fabs(dInterpCumul - dRealCumul);
        }

        (*dPoints)[i]     = dNewPoint;
        (*dCumulative)[i] = dNewCumul;
        dLastPoint        = dNewPoint;
        dLastCumul        = dNewCumul;
        dLastDensity      = dNewDensity;
        dNewPoint *= dShift;
        i++;

        if (i >= lAllocTotPoints)
        {
            lAllocTotPoints += *lNumPoints;
            *dPoints     = realloc(*dPoints, lAllocTotPoints * sizeof(double));
            *dCumulative = realloc(*dCumulative, lAllocTotPoints * sizeof(double));
        }
    }

    *lNumPoints = i;

    return err;
}

static Err smm_simulate(
    int nswaps, smm_swap* swaps, double** corr_matrix, smm_params* params, double** samples)
{
    Err        err = NULL;
    SRandomGen rg;
    double **  chol = NULL, **gg = NULL;
    long       i;
    int        j, k;
    double     t = swaps[0].fxg_time, sqrtt = sqrt(t), tmp;
    const long seed = -12345678;

    memset(&rg, 0, sizeof(SRandomGen));

    chol = dmatrix(0, nswaps - 1, 0, nswaps - 1);
    gg   = dmatrix(0, params->npaths - 1, 0, nswaps - 1);
    if (!chol || !gg)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    err = choldc(nswaps, corr_matrix, chol);
    if (err)
        goto FREE_RETURN;

    // Generate uncorrelated gaussians:

    err = ABS_Init(&rg, seed, params->npaths, nswaps, 0);
    if (err)
        goto FREE_RETURN;

    for (i = 0; i < params->npaths; i++)
    {
        for (j = 0; j < nswaps; j++)
        {
            err = rg.Gauss(&rg, &gg[i][j]);
            if (err)
                goto FREE_RETURN;
        }
    }

    // Correlate gauusians & calculate log-normals:

    for (i = 0; i < params->npaths; i++)
    {
        for (j = 0; j < nswaps; j++)
        {
            tmp = 0.0;
            for (k = 0; k < nswaps; k++)
                tmp += chol[j][k] * gg[i][k];

            samples[i][j] = swaps[j].cms_fwd /*swaps[j].swap0 */ *
                            exp(swaps[j].sigbeta1 * sqrtt * tmp -
                                0.5 * swaps[j].sigbeta1 * swaps[j].sigbeta1 * t);
        }
    }

FREE_RETURN:
    ABS_Free(&rg);
    if (chol)
        free_dmatrix(chol, 0, nswaps - 1, 0, nswaps - 1);
    if (gg)
        free_dmatrix(gg, 0, params->npaths - 1, 0, nswaps - 1);

    return err;
}

static Err smm_simulate2(
    int nswaps, smm_swap* swaps, double** corr_matrix, smm_params* params, double** samples)
{
    Err                   err = NULL;
    int                   i;
    long*                 npts = NULL;
    double **             pts = NULL, **cumul = NULL;
    COPULAGAUSSIAN_Params copula_prm;

    // Allocate memory:
    npts  = (long*)calloc(nswaps, sizeof(long));
    pts   = (double**)calloc(nswaps, sizeof(double*));
    cumul = (double**)calloc(nswaps, sizeof(double*));

    if (!npts || !pts || !cumul)
    {
        err = serror("Memory failure in smm_simulate");
        goto FREE_RETURN;
    }
    for (i = 0; i < nswaps; i++)
    {
        npts[i]  = params->cum_npts;
        pts[i]   = (double*)calloc(npts[i], sizeof(double));
        cumul[i] = (double*)calloc(npts[i], sizeof(double));
        if (!pts[i] || !cumul[i])
        {
            err = serror("Memory failure in smm_simulate");
            goto FREE_RETURN;
        }
    }

    // Calculate marginal cumulatives for all swaps:
    for (i = 0; i < nswaps; i++)
    {
        err = smm_BMM_linterp_cumulative(
            &swaps[i], params->maxerror, params->nstd, params->nstd, &npts[i], &pts[i], &cumul[i]);
        if (err)
            goto FREE_RETURN;
    }

    // Get samples from the copula:

    copula_gaussian_set_default_num_params(&copula_prm);

    copula_prm.iHasSavedGaussian = 0;
    copula_prm.dSavedGaussian    = NULL;
    copula_prm.eMCType           = params->mc_type;
    copula_prm.lNumPaths         = params->npaths;
    copula_prm.iNbFwd            = nswaps;
    copula_prm.iInterpMethod     = 5;  // quadratic

    err = copula_gaussian_numer(nswaps, pts, cumul, npts, 0, corr_matrix, &copula_prm, samples);
    if (err)
        goto FREE_RETURN;

FREE_RETURN:

    if (pts)
        for (i = 0; i < nswaps; i++)
            free(pts[i]);
    free(pts);
    if (cumul)
        for (i = 0; i < nswaps; i++)
            free(cumul[i]);
    free(cumul);
    free(npts);

    return err;
}

#define MAXNDATES 3000

static Err smm_setup_swaps(
    char*       yc,
    char*       vc,
    smm_params* params,
    char*       freq_str,
    char*       bas_str,
    char*       ref_rate,
    long        ex_date,
    long        last_date,
    smm_swap**  swaps,
    int*        nswaps)
{
    Err            err = NULL;
    long           dates[MAXNDATES], start_date, theo_date, act_date, temp_date, today;
    SrtCompounding fix_freq, flt_freq;
    SrtBasisCode   fix_basis, flt_basis;
    int            i, j, k, spotlag;
    double         flt_leg;
    smm_swap*      swp;
    double         sabr_sigbeta, sabr_alpha, sabr_beta, sabr_rho, power, calib_pres;
    double         imp_vols[3], std;
    const double   BMM_pi = 0.3;
    SrtCrvPtr      crv    = lookup_curve(yc);

    // For CMS forward calculation:
    CMSVol    cms_vol;
    CMSParams cms_params;
    double    cms_opt_atm;

    memset(&cms_vol, 0, sizeof(CMSVol));
    memset(&cms_params, 0, sizeof(CMSParams));

    if (!crv)
        return serror("Yield curve not found");
    today   = get_today_from_curve(crv);
    spotlag = get_spotlag_from_curve(crv);

    start_date = add_unit(ex_date, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);

    err = interp_compounding(freq_str, &fix_freq);
    if (err)
        return err;

    err = interp_basis(bas_str, &fix_basis);
    if (err)
        return err;

    err = swp_f_get_ref_rate_details(ref_rate, &flt_basis, &flt_freq);
    if (err)
        return err;

    theo_date = act_date = start_date;
    i                    = 0;
    do
    {
        theo_date  = add_unit(theo_date, 12 / fix_freq, SRT_MONTH, NO_BUSDAY_CONVENTION);
        dates[i++] = act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    } while (act_date < last_date);

    *nswaps = i;
    *swaps  = (smm_swap*)calloc(*nswaps, sizeof(smm_swap));
    if (!*swaps)
        return serror("Memory failure");

    for (i = 0; i < *nswaps; i++)
    {
        swp             = *swaps + i;
        swp->fxg_date   = ex_date;
        swp->start_date = start_date;
        swp->end_date   = dates[i];
        swp->fxg_time   = (swp->fxg_date - today) * YEARS_IN_DAY;
        swp->start_time = (swp->start_date - today) * YEARS_IN_DAY;
        swp->end_time   = (swp->end_date - today) * YEARS_IN_DAY;

        swp->df_fix   = swp_f_df(today, swp->fxg_date, yc);
        swp->df_start = swp_f_df(today, swp->start_date, yc);
        swp->df_end   = swp_f_df(today, swp->end_date, yc);

        swp->nfix = swp->nflt = i + 1;
        swp->fix_d            = (long*)calloc(swp->nfix, sizeof(long));
        swp->fix_t            = (double*)calloc(swp->nfix, sizeof(double));
        swp->fix_cvg          = (double*)calloc(swp->nfix, sizeof(double));
        swp->fix_df           = (double*)calloc(swp->nfix, sizeof(double));
        swp->flt_d            = (long*)calloc(swp->nflt, sizeof(long));
        swp->flt_t            = (double*)calloc(swp->nflt, sizeof(double));
        swp->flt_cvg          = (double*)calloc(swp->nflt, sizeof(double));
        swp->flt_df           = (double*)calloc(swp->nflt, sizeof(double));
        swp->flt_spr          = (double*)calloc(swp->nflt, sizeof(double));

        if (!swp->fix_d || !swp->fix_t || !swp->fix_cvg || !swp->fix_df || !swp->flt_d ||
            !swp->flt_t || !swp->flt_cvg || !swp->flt_df || !swp->flt_spr)
            return serror("Memory failure");

        memcpy(swp->fix_d, dates, swp->nfix * sizeof(long));
        for (j = 0; j < swp->nfix; j++)
        {
            swp->fix_t[j] = (swp->fix_d[j] - today) * YEARS_IN_DAY;
            swp->fix_cvg[j] =
                coverage((j == 0 ? swp->start_date : swp->fix_d[j - 1]), swp->fix_d[j], fix_basis);
            swp->fix_df[j] = swp_f_df(today, swp->fix_d[j], yc);
        }
        memcpy(swp->flt_d, dates, swp->nflt * sizeof(long));
        for (j = 0; j < swp->nflt; j++)
        {
            swp->flt_t[j] = (swp->flt_d[j] - today) * YEARS_IN_DAY;
            swp->flt_cvg[j] =
                coverage((j == 0 ? swp->start_date : swp->flt_d[j - 1]), swp->flt_d[j], flt_basis);
            swp->flt_df[j] = swp_f_df(today, swp->flt_d[j], yc);

            // Calculate spread:

            theo_date = act_date = (j == 0 ? swp->start_date : swp->flt_d[j - 1]);
            flt_leg              = 0.0;
            for (k = flt_freq / fix_freq - 1; k >= 0; k--)
            {
                theo_date = add_unit(theo_date, 12 / flt_freq, SRT_MONTH, NO_BUSDAY_CONVENTION);
                temp_date =
                    (k == 0 ? swp->flt_d[j] : bus_date_method(theo_date, MODIFIED_SUCCEEDING));
                flt_leg += swp_f_spread(act_date, temp_date, ref_rate) *
                           coverage(act_date, temp_date, flt_basis) *
                           swp_f_df(today, temp_date, yc);
                act_date = temp_date;
            }
            swp->flt_spr[j] = flt_leg / (swp->flt_cvg[j] * swp->flt_df[j]);
        }
        // Calculate swap0 & lvl0:

        swp->lvl0 = flt_leg = 0.0;
        for (j = 0; j < swp->nfix; j++)
            swp->lvl0 += swp->fix_cvg[j] * swp->fix_df[j];
        for (j = 0; j < swp->nflt; j++)
            flt_leg += swp->flt_spr[j] * swp->flt_cvg[j] * swp->flt_df[j];
        flt_leg += swp->df_start - swp->df_end;

        swp->swap0 = flt_leg / swp->lvl0;
        swp->lvlc0 = level_cash(swp->swap0, swp->nfix, swp->fix_cvg);

        // Now calibrate BMM parameters for the swap:

        // Get the SABR params from the vol curve first:
        err = swp_f_SABRvol(
            vc, swp->start_date, swp->end_date, 0.0, &sabr_sigbeta, &power, SABR_BETAVOL);
        if (err)
            return err;

        err =
            swp_f_SABRvol(vc, swp->start_date, swp->end_date, 0.0, &sabr_alpha, &power, SABR_ALPHA);
        if (err)
            return err;

        err = swp_f_SABRvol(vc, swp->start_date, swp->end_date, 0.0, &sabr_beta, &power, SABR_BETA);
        if (err)
            return err;

        err = swp_f_SABRvol(vc, swp->start_date, swp->end_date, 0.0, &sabr_rho, &power, SABR_RHO);
        if (err)
            return err;

        err = srt_f_optsarbvol(
            swp->swap0,
            swp->swap0,
            swp->fxg_time,
            sabr_sigbeta,
            sabr_alpha,
            sabr_beta,
            sabr_rho,
            SRT_BETAVOL,
            SRT_LOGNORMAL,
            &imp_vols[0]);
        if (err)
            return err;

        imp_vols[0] *= params->vol_mult;
        std = imp_vols[0] * sqrt(swp->fxg_time);

        swp->swpn_strikes[0] = swp->swap0;
        swp->swpn_strikes[1] = swp->swap0 * exp(-params->smile_nstd * std);
        swp->swpn_strikes[2] = swp->swap0 * exp(params->smile_nstd * std);

        sabr_alpha = 0.0;
        sabr_beta  = 1.0;
        sabr_rho   = 0.0;

        err = srt_f_optsarbvol(
            swp->swap0,
            swp->swpn_strikes[1],
            swp->fxg_time,
            imp_vols[0],
            sabr_alpha,
            sabr_beta,
            sabr_rho,
            SRT_LOGNORMAL,
            SRT_LOGNORMAL,
            &imp_vols[1]);
        if (err)
            return err;

        err = srt_f_optsarbvol(
            swp->swap0,
            swp->swpn_strikes[2],
            swp->fxg_time,
            imp_vols[0],
            sabr_alpha,
            sabr_beta,
            sabr_rho,
            SRT_LOGNORMAL,
            SRT_LOGNORMAL,
            &imp_vols[2]);
        if (err)
            return err;

        for (j = 0; j < 3; j++)
            swp->swpn_values[j] = srt_f_optblksch(
                swp->swap0,
                swp->swpn_strikes[j],
                imp_vols[j],
                swp->fxg_time,
                swp->lvl0,
                SRT_CALL,
                PREMIUM);

        swp->pi = BMM_pi;

        // Calculate CMS forward:

        CMSTECinitParams(&cms_params);
        CMSTECinitvol(&cms_vol);

        cms_vol.comp_method  = 0;
        cms_vol.voltype_used = SRT_BETAVOL;
        cms_vol.input_method = 0;

        err = srt_f_optsarbvol(
            swp->swap0,
            swp->swap0,
            swp->fxg_time,
            imp_vols[0],
            sabr_alpha,
            sabr_beta,
            sabr_rho,
            SRT_LOGNORMAL,
            SRT_BETAVOL,
            &cms_vol.vol);
        if (err)
            return err;

        cms_vol.beta  = sabr_beta;
        cms_vol.alpha = sabr_alpha;
        cms_vol.rho   = sabr_rho;

        //		cms_params.Nx = 201;

        //		err = swp_f_cmsrate(swp->swap0, swp->nfix, fix_freq, sabr_sigbeta,
        //swp->fxg_time, 			0.0, 0.0, 0, SRT_LOGNORMAL, &swp->cms_fwd);

        err = swp_f_CMSTEC_rate(
            swp->swap0,  // forward
            swp->nfix,  // nfp - (we suppose that we're always taking standard swaps as underlyings)
            fix_freq,   // standard swap freq
            1.0,        // so, ratio is always = 1
            swp->fxg_time,  // maturity
            1.0,            // CMS
            0.0,            // margin
            0.0,            // delay
            1.0,            // df_pay / df_start
            &cms_vol,
            NULL,
            &cms_params,
            &swp->cms_fwd);
        if (err)
            return err;

        err = swp_f_CMSTEC_Option(
            swp->swap0,
            swp->nfix,
            fix_freq,
            1.0,
            swp->fxg_time,
            swp->cms_fwd,
            SRT_RECEIVER,
            1.0,
            0.0,
            0.0,
            1.0,
            &cms_vol,
            NULL,
            &cms_params,
            &cms_opt_atm);
        if (err)
            return err;

        err = srt_f_optimpvol(
            cms_opt_atm,
            swp->cms_fwd,
            swp->cms_fwd,
            swp->fxg_time,
            1.0,
            SRT_PUT,
            SRT_LOGNORMAL,
            &imp_vols[0]);
        if (err)
            return err;

        err = srt_f_optsarbvol(
            swp->cms_fwd,
            swp->cms_fwd,
            swp->fxg_time,
            imp_vols[0],
            sabr_alpha,
            sabr_beta,
            sabr_rho,
            SRT_LOGNORMAL,
            SRT_BETAVOL,
            &sabr_sigbeta);
        if (err)
            return err;

        //		swp->cms_fwd *= exp(-0.0010 * (swp->end_time - swp->start_time));

        err = BMMCalibOnSabrStates(
            swp->cms_fwd,  // swp->swap0,
            swp->fxg_time,
            sabr_sigbeta,
            sabr_beta,
            sabr_alpha,
            sabr_rho,
            swp->pi,
            params->smile_nstd,
            &swp->fwd1,
            &swp->fwd2,
            &swp->sigbeta1,
            &swp->sigbeta2,
            &swp->pi,
            SRT_BETAVOL,
            &calib_pres);

        if (err)
            return err;
        if (calib_pres > 1e-5)
            return serror("BMM calibration failed");

        swp->beta = sabr_beta;
    }
    return NULL;
}

static void smm_free_swap(smm_swap* swp)
{
    free(swp->fix_d);
    free(swp->fix_t);
    free(swp->fix_cvg);
    free(swp->fix_df);
    free(swp->flt_d);
    free(swp->flt_t);
    free(swp->flt_cvg);
    free(swp->flt_df);
    free(swp->flt_spr);
}

#define NPASS 0

static void smm_strip(
    int       nswaps,
    smm_swap* swaps,
    double*   orig_samples,
    double*   swp_times,
    double*   swp_dfs,
    double*   zc_rate)
{
    double samples[MAXNDATES];
    int    i, swp_ndfs = swaps[nswaps - 1].nfix + 1;
    double fix_leg, flt_leg, fix_tmp, flt_tmp;
    int    ipass = 0;

    memcpy(samples, orig_samples, nswaps * sizeof(double));

PASS_START:
    fix_leg = flt_leg = 0.0;
    swp_dfs[0]        = exp(-samples[0] * (swp_times[0] - swaps[0].fxg_time));
    zc_rate[0]        = samples[0];

    for (i = 0; i < nswaps; i++)
    {
        fix_tmp = swaps[i].fix_cvg[swaps[i].nfix - 1];
        flt_tmp = swaps[i].flt_spr[swaps[i].nflt - 1] * swaps[i].flt_cvg[swaps[i].nflt - 1];
        swp_dfs[i + 1] =
            (swp_dfs[0] - samples[i] * fix_leg + flt_leg) / (samples[i] * fix_tmp + 1.0 - flt_tmp);
        //		(swp_dfs[0] - samples[i] * level_cash(samples[i], swaps[i].nfix,
        //swaps[i].fix_cvg) + flt_leg) / 		(1.0 - flt_tmp);

        if (swp_dfs[i + 1] > 0)
            zc_rate[i + 1] = -log(swp_dfs[i + 1]) / (swp_times[i + 1] - swaps[0].fxg_time);
        else
            zc_rate[i + 1] = 999.99;
        fix_leg += fix_tmp * swp_dfs[i + 1];
        flt_leg += flt_tmp * swp_dfs[i + 1];

        if (ipass < NPASS)
            samples[i] = orig_samples[i] * swaps[i].lvl0 / swaps[i].df_start / fix_leg;
    }

    if (ipass++ < NPASS)
        goto PASS_START;
}

static void smm_interp_dfs(
    int     swp_ndfs,
    double* swp_times,
    double* swp_dfs,
    double* zc_rate,
    double  mat,
    int     ndfs,
    double* times,
    double* dfs)
{
    int    i, idx = 0;
    double zcr, coef;

    for (i = 0; i < ndfs; i++)
    {
        if (times[i] <= swp_times[0])
            dfs[i] = swp_dfs[0];
        else if (times[i] >= swp_times[swp_ndfs - 1])
            dfs[i] = swp_dfs[swp_ndfs - 1];
        else
        {
            while (times[i] >= swp_times[idx + 1])
                idx++;
            coef = (times[i] - swp_times[idx]) / (swp_times[idx + 1] - swp_times[idx]);
            if (swp_dfs[idx + 1] > 0 && swp_dfs[idx] > 0)
            {
                zcr    = coef * zc_rate[idx + 1] + (1.0 - coef) * zc_rate[idx];
                dfs[i] = exp(-zcr * (times[i] - mat));
            }
            else
                dfs[i] = coef * swp_dfs[idx + 1] + (1.0 - coef) * swp_dfs[idx];
        }
    }
}

static Err smm_match(
    int nswaps, smm_swap* swaps, smm_params* params, double** swp_dfs, double* weights)
{
    Err      err = NULL;
    long     j, iter, failinfo;
    NagError fail;
    int      i, k, n, swp_ndfs = swaps[nswaps - 1].nfix + 1;
    long     N = params->npaths;
    int      M = swp_ndfs + nswaps * 4;  // total number of things to match
    double **A = NULL, **AA = NULL;
    double * b = NULL, *sv = NULL, *e = NULL;
    double   svmin, sum, sum_s[3], fix_leg, flt_leg;

    memset(&fail, 0, sizeof(NagError));

    A  = dmatrix(0, M - 1, 0, N - 1);
    AA = dmatrix(0, M - 1, 0, M - 1);
    b  = (double*)calloc(M, sizeof(double));
    sv = (double*)calloc(M, sizeof(double));
    e  = (double*)calloc(M, sizeof(double));

    if (!A || !AA || !b || !sv || !e)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    // Calculate A and b = Ac - d:

    for (i = 0; i < swp_ndfs - 1; i++)
    {
        for (j = 0, sum = 0.0; j < N; j++)
        {
            A[i][j] = swp_dfs[j][i + 1] / swp_dfs[j][0];
            sum += A[i][j] * weights[j];
        }
        b[i] = sum - swaps[nswaps - 1].fix_df[i] / swaps[nswaps - 1].df_start;
    }

    for (j = 0; j < N; j++)
        A[i][j] = 1.0;
    b[i++] = 0.0;

    for (k = 0; k < nswaps; k++)
    {
        sum = sum_s[0] = sum_s[1] = sum_s[2] = 0.0;
        for (j = 0; j < N; j++)
        {
            // Calculate level and CMS at the simulated point:

            fix_leg = flt_leg = 0.0;
            for (n = 0; n < swaps[k].nfix; n++)
            {
                fix_leg += swaps[k].fix_cvg[n] * swp_dfs[j][n + 1];
                flt_leg += swaps[k].flt_cvg[n] * swp_dfs[j][n + 1] * swaps[k].flt_spr[n];
            }
            flt_leg += swp_dfs[j][0] - swp_dfs[j][k + 1];

            A[i][j] = flt_leg / fix_leg / swp_dfs[j][0];  // CMS
            sum += A[i][j] * weights[j];

            for (n = 0; n < 3; n++)
            {
                A[i + n + 1][j] = (flt_leg - swaps[k].swpn_strikes[n] * fix_leg) / swp_dfs[j][0];
                if (A[i + n + 1][j] < 0.0)
                    A[i + n + 1][j] = 0.0;  // Swaptions
                sum_s[n] += A[i + n + 1][j] * weights[j];
            }
        }
        b[i++] = sum - swaps[k].cms_fwd;
        for (n = 0; n < 3; n++)
            b[i++] = sum_s[n] - swaps[k].swpn_values[n] / swaps[k].df_start;
    }

    // Calculate AA':

    for (i = 0; i < M; i++)
    {
        for (k = i; k < M; k++)
        {
            for (j = 0, sum = 0.0; j < N; j++)
                sum += A[i][j] * A[k][j];
            AA[i][k] = AA[k][i] = sum;
        }
    }

    // SVD (AA is replaced by P' and b - by Q'b):

    nag_real_svd(
        M, M, &AA[0][0], M, 1, b, 1, 0, NULL, 0, sv, 1, NULL, 0, &iter, e, &failinfo, &fail);
    if (fail.code != NE_NOERROR)
    {
        err = serror(fail.message);
        goto FREE_RETURN;
    }

    // Zero small singular values and calculate diag{1/SVi} * Q'b:

    svmin = sv[0] * params->sv_min;
    for (i = 0; i < M; i++)
    {
        if (sv[i] < svmin)
            b[i] = 0.0;
        else
            b[i] /= sv[i];
    }
    // Calculate P * diag{1/SVi} * Q'b = (AA')^(-1) * b:

    for (i = 0; i < M; i++)
    {
        for (k = 0, sum = 0.0; k < M; k++)
            sum += AA[k][i] * b[k];
        e[i] = sum;
    }

    // Finally: adjust weights: x = c - A'(AA')^(-1)(Ac - d)

    for (j = 0; j < N; j++)
    {
        for (i = 0, sum = 0.0; i < M; i++)
            sum += A[i][j] * e[i];
        weights[j] -= sum;
    }

FREE_RETURN:
    if (A)
        free_dmatrix(A, 0, M - 1, 0, N - 1);
    if (AA)
        free_dmatrix(AA, 0, M - 1, 0, M - 1);
    free(b);
    free(sv);
    free(e);

    return err;
}

Err smm_otc_price(
    char*         yc,
    char*         vc,
    char*         freq_str,
    char*         bas_str,
    char*         ref_rate,
    char*         corr_cube,
    smm_params*   params,
    SProductDesc* g,
    int           iex,
    double*       pv,
    double*       stddev)
{
    Err       err   = NULL;
    smm_swap* swaps = NULL;
    int       i, j, nswaps, swp_ndfs;
    double ** corr_matrix = NULL, **samples = NULL, **pvs = NULL;
    double ** swp_dfs = NULL, **zc_rate = NULL, *swp_times = NULL, *weights = NULL;
    double    tmp, dfs[MAXNDATES], *dfs_ = dfs;
    double    sabr_sigbeta[MAXNDATES], sabr_alpha[MAXNDATES], sabr_beta[MAXNDATES],
        sabr_rho[MAXNDATES];
    double *     sabr_sigbeta_ = sabr_sigbeta, *sabr_alpha_ = sabr_alpha;
    double *     sabr_beta_ = sabr_beta, *sabr_rho_ = sabr_rho;
    SSabrMktDesc mkt_sabr_desc = {&sabr_sigbeta_, &sabr_alpha_, &sabr_beta_, &sabr_rho_, 0, NULL};
    SMktData     mkt_data      = {&dfs_, MKT_SABR, &mkt_sabr_desc};
    double       power, *iv = NULL, *fee = NULL;
    long         today, start_corr, end_corr, slide, ipath;
    SrtCrvPtr    crv = lookup_curve(yc);

    if (!crv)
        return serror("Yield curve not found");
    today = get_today_from_curve(crv);

    // Calculate today's iv first:

    iv  = (double*)calloc(g->ninst, sizeof(double));
    fee = (double*)calloc(g->ninst, sizeof(double));
    if (!iv || !fee)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    if (g->nvol)
        for (i = 0; i < g->nvol[0][iex]; i++)
        {
            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i],
                g->vol_end[0][iex][i],
                0.0,
                &sabr_sigbeta[i],
                &power,
                SABR_BETAVOL);
            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i],
                g->vol_end[0][iex][i],
                0.0,
                &sabr_alpha[i],
                &power,
                SABR_ALPHA);
            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i],
                g->vol_end[0][iex][i],
                0.0,
                &sabr_beta[i],
                &power,
                SABR_BETA);
            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i],
                g->vol_end[0][iex][i],
                0.0,
                &sabr_rho[i],
                &power,
                SABR_RHO);
            if (err)
                goto FREE_RETURN;
        }

    if (params->use_fee)
    {
        for (i = 0; i < g->nmat[0][iex]; i++)
            dfs[i] = swp_f_df(today, g->mat_d[0][iex][i], yc);

        err = g->Payoff(g, iex, g->ex[iex], &mkt_data, iv);
        if (err)
            goto FREE_RETURN;
    }

    slide = g->ex_d[iex] - today;

    if (g->nvol && (params->vol_hyp & 1))  // sigma-beta sliding
    {
        for (i = 0; i < g->nvol[0][iex]; i++)
        {
            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i] - slide,
                g->vol_end[0][iex][i] - slide,
                0.0,
                &sabr_sigbeta[i],
                &power,
                SABR_BETAVOL);
            if (err)
                goto FREE_RETURN;
        }
    }
    if (g->nvol && (params->vol_hyp & 2))  // SABR params sliding
    {
        for (i = 0; i < g->nvol[0][iex]; i++)
        {
            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i] - slide,
                g->vol_end[0][iex][i] - slide,
                0.0,
                &sabr_alpha[i],
                &power,
                SABR_ALPHA);
            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i] - slide,
                g->vol_end[0][iex][i] - slide,
                0.0,
                &sabr_beta[i],
                &power,
                SABR_BETA);
            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                vc,
                g->vol_start[0][iex][i] - slide,
                g->vol_end[0][iex][i] - slide,
                0.0,
                &sabr_rho[i],
                &power,
                SABR_RHO);
            if (err)
                goto FREE_RETURN;
        }
    }

    // Setup swaps:

    err = smm_setup_swaps(
        yc,
        vc,
        params,
        freq_str,
        bas_str,
        ref_rate,
        g->ex_d[iex],
        g->mat_d[0][iex][g->nmat[0][iex] - 1],
        &swaps,
        &nswaps);
    if (err)
        goto FREE_RETURN;

    // Fill corr matrix:

    corr_matrix = dmatrix(0, nswaps - 1, 0, nswaps - 1);
    if (!corr_matrix)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    for (i = 0; i < nswaps; i++)
    {
        for (j = 0; j < i; j++)
        {
            start_corr = today + swaps[j].end_date - swaps[j].start_date;
            end_corr   = start_corr + swaps[i].end_date - swaps[j].end_date;
            err = swp_f_vol(corr_cube, start_corr, end_corr, g->ex[iex], &corr_matrix[i][j], &tmp);
            if (err)
                goto FREE_RETURN;
            corr_matrix[j][i] = corr_matrix[i][j];
        }
        corr_matrix[i][i] = 1.0;
    }

    err = PositiveMatrix(corr_matrix, nswaps);
    if (err)
        goto FREE_RETURN;

    // Simulate samples:

    samples = dmatrix(0, params->npaths - 1, 0, nswaps - 1);
    if (!samples)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    err = smm_simulate(nswaps, swaps, corr_matrix, params, samples);
    if (err)
        goto FREE_RETURN;

    // Strip dfs:

    swp_ndfs  = swaps[nswaps - 1].nfix + 1;
    swp_times = (double*)calloc(swp_ndfs, sizeof(double));
    swp_dfs   = dmatrix(0, params->npaths - 1, 0, swp_ndfs - 1);
    zc_rate   = dmatrix(0, params->npaths - 1, 0, swp_ndfs - 1);

    if (!swp_times || !swp_dfs || !zc_rate)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    swp_times[0] = swaps[0].start_time;
    memcpy(swp_times + 1, swaps[nswaps - 1].fix_t, swaps[nswaps - 1].nfix * sizeof(double));

    for (ipath = 0; ipath < params->npaths; ipath++)
        smm_strip(nswaps, swaps, samples[ipath], swp_times, swp_dfs[ipath], zc_rate[ipath]);

    // Match stuff:

    weights = (double*)calloc(params->npaths, sizeof(double));
    if (!weights)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    for (ipath = 0; ipath < params->npaths; ipath++)
        weights[ipath] = 1.0 / params->npaths;  // Initial point

    if (params->do_match)
    {
        err = smm_match(nswaps, swaps, params, swp_dfs, weights);
        if (err)
            goto FREE_RETURN;
    }

    // Pricing:

    pvs = dmatrix(0, params->npaths - 1, 0, g->ninst - 1);
    if (!pvs)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    memset(pv, 0, g->ninst * sizeof(double));
    memset(stddev, 0, g->ninst * sizeof(double));

    for (ipath = 0; ipath < params->npaths; ipath++)
    {
        smm_interp_dfs(
            swp_ndfs,
            swp_times,
            swp_dfs[ipath],
            zc_rate[ipath],
            swaps[0].fxg_time,
            g->nmat[0][iex],
            g->mat[0][iex],
            dfs);

        err = g->Payoff(g, iex, g->ex[iex], &mkt_data, pvs[ipath]);
        if (err)
            goto FREE_RETURN;

        for (i = 0; i < g->ninst; i++)
        {
            pvs[ipath][i] /= swp_dfs[ipath][0];
            pv[i] += pvs[ipath][i] * weights[ipath] * swaps[0].df_start;
        }
    }

    // Calculate fees:
    for (i = 0; i < g->ninst; i++)
        fee[i] = (params->use_fee ? (iv[i] - pv[i]) / swaps[0].df_start : 0.0);

    // Apply params->otc_type:

    memset(pv, 0, g->ninst * sizeof(double));

    for (ipath = 0; ipath < params->npaths; ipath++)
    {
        for (i = 0; i < g->ninst; i++)
        {
            pvs[ipath][i] += fee[i];
            if (params->otc_type == 1)
                pvs[ipath][i] = max(pvs[ipath][i], 0.0);  // call
            if (params->otc_type == 2)
                pvs[ipath][i] = max(-pvs[ipath][i], 0.0);  // put
            pv[i] += pvs[ipath][i] * weights[ipath] * swaps[0].df_start;
        }
    }

    // Calculate stddev:

    for (i = 0; i < g->ninst; i++)
    {
        for (ipath = 0; ipath < params->npaths; ipath++)
        {
            tmp = pvs[ipath][i] * swaps[0].df_start - pv[i];
            stddev[i] += tmp * tmp / (params->npaths - 1);
        }
        stddev[i] = sqrt(stddev[i] / params->npaths);
    }

FREE_RETURN:
    if (pvs)
        free_dmatrix(pvs, 0, params->npaths - 1, 0, g->ninst - 1);

    free(swp_times);
    if (swp_dfs)
        free_dmatrix(swp_dfs, 0, params->npaths - 1, 0, swp_ndfs - 1);
    if (zc_rate)
        free_dmatrix(zc_rate, 0, params->npaths - 1, 0, swp_ndfs - 1);
    free(weights);

    if (samples)
        free_dmatrix(samples, 0, params->npaths - 1, 0, nswaps - 1);
    if (corr_matrix)
        free_dmatrix(corr_matrix, 0, nswaps - 1, 0, nswaps - 1);
    if (swaps)
        for (i = 0; i < nswaps; i++)
            smm_free_swap(&swaps[i]);
    free(swaps);

    free(iv);
    free(fee);

    return err;
}

static Err IndexDates(int ndates, long* dates, int* new_ndates, long* new_dates, long* index)
{
    Err  err = NULL;
    long i, j, lastd = -1, idx[10 * MAXNDATES];

    if (ndates > 10 * MAXNDATES)
        return serror("Too many dates in IndexDates");

    err = indexx_ll(dates, idx, ndates);
    if (err)
        return err;

    for (i = j = 0; i < ndates; i++)
    {
        if (dates[idx[i]] != lastd)
            new_dates[j++] = lastd = dates[idx[i]];
        index[idx[i]] = j - 1;
    }
    *new_ndates = j;
    return NULL;
}

static Err smm_swaptions_payoff(
    SProductDesc* g,
    int           idx,  // must be = 0
    double        time,
    SMktData*     mkt_data,  // market data (dfs)
    double*       pv)              // ninst sized vector with PVs (output)
{
    double*             dfs      = mkt_data->dfs[0];
    smm_swaptions_desc* swp_desc = (smm_swaptions_desc*)g->spec_desc;
    int                 i, j, cf_index;

    // test CMS:
    double lvl;

    for (i = 0; i < g->ninst; i++)
    {
        pv[i] = 0.0;
        for (j = 0; j < swp_desc->ncf[i]; j++)
        {
            cf_index = swp_desc->cf_idx[i] + j;
            pv[i] += swp_desc->cf[cf_index] * dfs[swp_desc->dfs_idx[cf_index]];
        }
        // test CMS:
        lvl = 0.0;
        /*		for (j = swp_desc->cffix_idx[i] - swp_desc->cf_idx[i]; j < swp_desc->ncf[i];
           j++)
                        {
                                cf_index = swp_desc->cf_idx[i] + j;
                                lvl += -swp_desc->cf[cf_index] / swp_desc->K[i] *
           dfs[swp_desc->dfs_idx[cf_index]];
                        }
                        pv[i] /= lvl;*/

        //		if (swp_desc->pay_rec == SRT_RECEIVER) pv[i] = -pv[i];
        //		if (pv[i] < 0.0) pv[i] = 0.0;
    }

    return NULL;
}

static Err smm_cms_payoff(
    SProductDesc* g,
    int           idx,  // must be = 0
    double        time,
    SMktData*     mkt_data,  // market data (dfs)
    double*       pv)              // ninst sized vector with PVs (output)
{
    double*             dfs      = mkt_data->dfs[0];
    smm_swaptions_desc* swp_desc = (smm_swaptions_desc*)g->spec_desc;
    int                 i, j, cf_index;

    // test CMS:
    double lvl;

    for (i = 0; i < g->ninst; i++)
    {
        pv[i] = 0.0;
        for (j = 0; j < swp_desc->ncf[i]; j++)
        {
            cf_index = swp_desc->cf_idx[i] + j;
            pv[i] += swp_desc->cf[cf_index] * dfs[swp_desc->dfs_idx[cf_index]];
        }
        // test CMS:
        lvl = 0.0;
        for (j = swp_desc->cffix_idx[i] - swp_desc->cf_idx[i]; j < swp_desc->ncf[i]; j++)
        {
            cf_index = swp_desc->cf_idx[i] + j;
            lvl += -swp_desc->cf[cf_index] / swp_desc->K[i] * dfs[swp_desc->dfs_idx[cf_index]];
        }
        pv[i] /= lvl;

        if (swp_desc->pay_rec == SRT_RECEIVER)
            pv[i] = -pv[i];
        if (pv[i] < 0.0)
            pv[i] = 0.0;
    }

    return NULL;
}

static Err smm_levelcash_payoff(
    SProductDesc* g,
    int           idx,  // must be = 0
    double        time,
    SMktData*     mkt_data,  // market data (dfs)
    double*       pv)              // ninst sized vector with PVs (output)
{
    double*             dfs      = mkt_data->dfs[0];
    smm_swaptions_desc* swp_desc = (smm_swaptions_desc*)g->spec_desc;
    int                 i, j, cf_index;
    double              cvg[1000];
    int                 ncvg;

    // test CMS:
    double lvl;

    for (i = 0; i < g->ninst; i++)
    {
        pv[i] = 0.0;
        for (j = 0; j < swp_desc->ncf[i]; j++)
        {
            cf_index = swp_desc->cf_idx[i] + j;
            pv[i] += swp_desc->cf[cf_index] * dfs[swp_desc->dfs_idx[cf_index]];
        }
        // test CMS:
        lvl  = 0.0;
        ncvg = 0;
        for (j = swp_desc->cffix_idx[i] - swp_desc->cf_idx[i]; j < swp_desc->ncf[i]; j++)
        {
            cf_index  = swp_desc->cf_idx[i] + j;
            cvg[ncvg] = -swp_desc->cf[cf_index] / swp_desc->K[i];
            lvl += cvg[ncvg] * dfs[swp_desc->dfs_idx[cf_index]];
            ncvg++;
        }
        pv[i] /= lvl;

        pv[i] = level_cash(pv[i], ncvg, cvg);

        //		if (swp_desc->pay_rec == SRT_RECEIVER) pv[i] = -pv[i];
        //		if (pv[i] < 0.0) pv[i] = 0.0;
    }

    return NULL;
}

Err smm_cif_init_product(
    SProductDesc* g,
    char*         yc,

    // Coupons:
    char*   freq_str,  // CMS fix freq
    char*   bas_str,   // CMS fix basis
    char*   ref_rate,  // CMS ref rate
    int     ncpn,
    long*   starts,
    long*   ends,
    double* alpha,
    double* beta,
    int*    floored,
    double* floor,
    int*    capped,
    double* cap,
    long*   cpn_pay,
    double* cpn_cvg,

    // Funding:
    int     fund_ncpn,
    long*   fund_start,  //	Start dates
    long*   fund_pay,    //	Pay dates
    char**  fund_basis,  //	Basis
    double* fund_spr,    //	Forward spreads
    double* fund_mrg,    //	Margins

    // Calls:
    int   ncall,    //	Number of calls
    long* ex_date,  //	Call dates
    long* set_date  //	Call settlement dates
)
{
    Err            err = NULL;
    int            i, j, n_all_dates;
    long           today, spotlag, theo_date, act_date, temp_date;
    long           all_dates[10 * MAXNDATES], new_dates[MAXNDATES];
    double         all_cash_flows[10 * MAXNDATES];
    SrtCrvPtr      crv = lookup_curve(yc);
    cif_desc*      cif;
    SrtCompounding fix_freq, flt_freq;
    SrtBasisCode   fix_basis, flt_basis, fund_bas_code;
    double         cvg;

    if (!crv)
        return serror("Yield curve not found");
    today   = get_today_from_curve(crv);
    spotlag = get_spotlag_from_curve(crv);

    err = interp_compounding(freq_str, &fix_freq);
    if (err)
        return err;

    err = interp_basis(bas_str, &fix_basis);
    if (err)
        return err;

    err = swp_f_get_ref_rate_details(ref_rate, &flt_basis, &flt_freq);
    if (err)
        return err;

    g->type = PRODUCT_CIF;
    g->nccy = g->ninst = 1;
    g->nex             = ncall;

    g->ex   = (double*)calloc(g->nex, sizeof(double));
    g->ex_d = (long*)calloc(g->nex, sizeof(long));

    g->nmat = imatrix(0, 0, 0, g->nex - 1);
    g->nvol = imatrix(0, 0, 0, g->nex - 1);

    g->spec_desc = calloc(1, sizeof(cif_desc));

    if (!g->ex || !g->ex_d || !g->nmat || !g->nvol || !g->spec_desc)
        return serror("Memory failure");

    memcpy(g->ex_d, ex_date, ncall * sizeof(long));
    for (i = 0; i < ncall; i++)
        g->ex[i] = (g->ex_d[i] - today) * YEARS_IN_DAY;

    cif         = (cif_desc*)g->spec_desc;
    n_all_dates = 0;

    // Init funding cpns:

    cif->fund_ncpn = fund_ncpn;
    cif->fund_cpn  = (cif_fund_cpn*)calloc(cif->fund_ncpn, sizeof(cif_fund_cpn));
    if (!cif->fund_cpn)
        return serror("Memory failure");

    for (i = 0; i < fund_ncpn; i++)
    {
        err = interp_basis(fund_basis[i], &fund_bas_code);
        if (err)
            return err;

        cvg                     = coverage(fund_start[i], fund_pay[i], fund_bas_code);
        cif->fund_cpn[i].cpn    = cvg * (fund_spr[i] + fund_mrg[i]);
        cif->fund_cpn[i].df_idx = n_all_dates;
        all_dates[n_all_dates]  = fund_pay[i];
        n_all_dates++;
    }

    return err;
}

Err smm_swaptions_init_product(
    SProductDesc* g,
    char*         yc,
    char*         freq_str,
    char*         bas_str,
    char*         ref_rate,
    char*         pay_rec,
    int           is_cms,
    long          start,
    int           nends,
    long*         ends,
    int*          nK,
    double**      K)
{
    Err                 err = NULL;
    int                 i, j, inst_idx, n_all_dates;
    long                today, spotlag, theo_date, act_date, temp_date;
    long                all_dates[10 * MAXNDATES], new_dates[MAXNDATES];
    double              all_cash_flows[10 * MAXNDATES];
    SrtCrvPtr           crv = lookup_curve(yc);
    smm_swaptions_desc* swp_desc;
    SrtCompounding      fix_freq, flt_freq;
    SrtBasisCode        fix_basis, flt_basis;

    if (!crv)
        return serror("Yield curve not found");
    today   = get_today_from_curve(crv);
    spotlag = get_spotlag_from_curve(crv);

    err = interp_compounding(freq_str, &fix_freq);
    if (err)
        return err;

    err = interp_basis(bas_str, &fix_basis);
    if (err)
        return err;

    err = swp_f_get_ref_rate_details(ref_rate, &flt_basis, &flt_freq);
    if (err)
        return err;

    g->type = PRODUCT_SMM_SWAPTIONS;
    g->nex = g->nccy = 1;
    g->ninst         = 0;
    for (i = 0; i < nends; i++)
        g->ninst += nK[i];

    g->ex        = (double*)calloc(1, sizeof(double));
    g->ex_d      = (long*)calloc(1, sizeof(long));
    g->nmat      = imatrix(0, 0, 0, 0);
    g->spec_desc = calloc(1, sizeof(smm_swaptions_desc));

    if (!g->ex || !g->ex_d || !g->nmat || !g->spec_desc)
        return serror("Memory failure");

    swp_desc         = (smm_swaptions_desc*)g->spec_desc;
    swp_desc->cf_idx = (long*)calloc(g->ninst, sizeof(long));
    swp_desc->ncf    = (int*)calloc(g->ninst, sizeof(long));

    if (!swp_desc->cf_idx || !swp_desc->ncf)
        return serror("Memory failure");

    err = interp_rec_pay(pay_rec, &swp_desc->pay_rec);
    if (err)
        return err;

    // test CMS:
    swp_desc->K         = (double*)calloc(g->ninst, sizeof(double));
    swp_desc->cffix_idx = (long*)calloc(g->ninst, sizeof(long));

    g->ex_d[0] = add_unit(start, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);
    g->ex[0]   = (g->ex_d[0] - today) * YEARS_IN_DAY;

    // Fill all_dates:

    n_all_dates = 0;
    inst_idx    = 0;
    for (i = 0; i < nends; i++)
    {
        for (j = 0; j < nK[i]; j++)
        {
            swp_desc->cf_idx[inst_idx] = n_all_dates;
            swp_desc->ncf[inst_idx]    = 0;

            // test CMS:
            swp_desc->K[inst_idx] = K[i][j];

            // Notionals:
            all_dates[n_all_dates]          = start;
            all_dates[n_all_dates + 1]      = ends[i];
            all_cash_flows[n_all_dates]     = 1.0;
            all_cash_flows[n_all_dates + 1] = -1.0;
            n_all_dates += 2;
            swp_desc->ncf[inst_idx] += 2;

            // Spreads:
            theo_date = act_date = temp_date = start;

            while (act_date < ends[i])
            {
                theo_date = add_unit(theo_date, 12 / flt_freq, SRT_MONTH, NO_BUSDAY_CONVENTION);
                act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
                if (act_date > ends[i])
                    act_date = ends[i];
                all_dates[n_all_dates]      = act_date;
                all_cash_flows[n_all_dates] = swp_f_spread(temp_date, act_date, ref_rate) *
                                              coverage(temp_date, act_date, flt_basis);
                n_all_dates++;
                swp_desc->ncf[inst_idx]++;
                temp_date = act_date;
            }

            // test CMS:
            swp_desc->cffix_idx[inst_idx] = n_all_dates;

            // Fixed leg:
            theo_date = act_date = temp_date = start;

            while (act_date < ends[i])
            {
                theo_date = add_unit(theo_date, 12 / fix_freq, SRT_MONTH, NO_BUSDAY_CONVENTION);
                act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
                if (act_date > ends[i])
                    act_date = ends[i];
                all_dates[n_all_dates]      = act_date;
                all_cash_flows[n_all_dates] = -K[i][j] * coverage(temp_date, act_date, fix_basis);
                n_all_dates++;
                swp_desc->ncf[inst_idx]++;
                temp_date = act_date;
            }
            inst_idx++;
        }
    }
    swp_desc->cf      = (double*)calloc(n_all_dates, sizeof(double));
    swp_desc->dfs_idx = (long*)calloc(n_all_dates, sizeof(long));

    if (!swp_desc->cf || !swp_desc->dfs_idx)
        return serror("Memory failure");
    memcpy(swp_desc->cf, all_cash_flows, n_all_dates * sizeof(double));

    err = IndexDates(n_all_dates, all_dates, &g->nmat[0][0], new_dates, swp_desc->dfs_idx);
    if (err)
        return err;

    g->mat   = f3tensor(0, 0, 0, 0, 0, g->nmat[0][0] - 1);
    g->mat_d = l3tensor(0, 0, 0, 0, 0, g->nmat[0][0] - 1);

    if (!g->mat || !g->mat_d)
        return serror("Memory failure");

    memcpy(g->mat_d[0][0], new_dates, g->nmat[0][0] * sizeof(long));

    for (i = 0; i < g->nmat[0][0]; i++)
        g->mat[0][0][i] = (g->mat_d[0][0][i] - today) * YEARS_IN_DAY;

    g->Payoff = (is_cms ? smm_cms_payoff : smm_swaptions_payoff);

    return NULL;
}

void smm_swaptions_free_product(SProductDesc* g)
{
    smm_swaptions_desc* swp_desc;

    if (g->mat)
        free_f3tensor(g->mat, 0, 0, 0, 0, 0, g->nmat[0][0] - 1);
    if (g->mat_d)
        free_l3tensor(g->mat_d, 0, 0, 0, 0, 0, g->nmat[0][0] - 1);

    if (g->spec_desc)
    {
        swp_desc = (smm_swaptions_desc*)g->spec_desc;
        free(swp_desc->cf_idx);
        free(swp_desc->ncf);
        free(swp_desc->cf);
        free(swp_desc->dfs_idx);

        // test CMS:
        free(swp_desc->K);
        free(swp_desc->cffix_idx);
    }
    free(g->spec_desc);

    free(g->ex);
    free(g->ex_d);
    if (g->nmat)
        free_imatrix(g->nmat, 0, 0, 0, 0);
}

static Err smm_dfs_payoff(
    SProductDesc* g,
    int           idx,  // must be = 0
    double        time,
    SMktData*     mkt_data,  // market data (dfs)
    double*       pv)              // ninst sized vector with PVs (output)
{
    double* dfs = mkt_data->dfs[0];
    int     i;

    for (i = 0; i < g->ninst; i++)
        pv[i] = dfs[i];

    return NULL;
}

Err smm_dfs_init_product(SProductDesc* g, char* yc, long ex_date, int ndates, long* dates)
{
    Err       err = NULL;
    long      today;
    int       i;
    SrtCrvPtr crv = lookup_curve(yc);

    if (!crv)
        return serror("Yield curve not found");
    today = get_today_from_curve(crv);

    g->nex = g->nccy = 1;
    g->ninst         = ndates;

    g->ex   = (double*)calloc(1, sizeof(double));
    g->ex_d = (long*)calloc(1, sizeof(long));

    g->ex_d[0] = ex_date;
    g->ex[0]   = (ex_date - today) * YEARS_IN_DAY;

    g->nmat       = imatrix(0, 0, 0, 0);
    g->nmat[0][0] = ndates;

    g->mat   = f3tensor(0, 0, 0, 0, 0, ndates - 1);
    g->mat_d = l3tensor(0, 0, 0, 0, 0, ndates - 1);

    memcpy(g->mat_d[0][0], dates, ndates * sizeof(long));
    for (i = 0; i < ndates; i++)
        g->mat[0][0][i] = (g->mat_d[0][0][i] - today) * YEARS_IN_DAY;

    g->Payoff = smm_dfs_payoff;

    return NULL;
}

void smm_dfs_free_product(SProductDesc* g)
{
    if (g->mat_d)
        free_l3tensor(g->mat_d, 0, 0, 0, 0, 0, g->ninst - 1);
    if (g->mat)
        free_f3tensor(g->mat, 0, 0, 0, 0, 0, g->ninst - 1);
    if (g->nmat)
        free_imatrix(g->nmat, 0, 0, 0, 0);
    free(g->ex_d);
    free(g->ex);
}
