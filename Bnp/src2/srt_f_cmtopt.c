#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_cmtopt.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "uttypes.h"

#define PROP_SHIFT 1.0e-03
#define ZERO_SHIFT 1.0e-06
#define MIN_SPREAD_MEAN_REVERSION -1.0
#define MAX_SPREAD_MEAN_REVERSION 1.0
#define MIN_SPREAD_FWD_LEVEL 0
#define MAX_SPREAD_FWD_LEVEL 1
#define MIN_SPREAD_LOC_VOL 1e-4
#define MAX_SPREAD_LOC_VOL 5e-1
#define MAX_RHO 1
#define MIN_RHO -MAX_RHO
#define MAX_EXP_ARG 25
#define MIN_EXP_ARG -25

#define SRT_BIG 1.0e+20

static double trunc_exp(double X)
{
    if (X > MAX_EXP_ARG)
        X = MAX_EXP_ARG;
    else if (X < MIN_EXP_ARG)
        X = MIN_EXP_ARG;

    return exp(X);
}

typedef struct
{
    String           swp_yc_name;
    SrtDiffusionType swp_vol_type;

    double yr_to_exp, cms_rate, cmt_rate, atm_sigma_beta, sabr_alpha, sabr_beta, sabr_rho,
        spot_spread, fwd_spread, *cmt_calibration_strikes;

    double** ppdParamBounds;

    long num_strike;
} static_params_struc;

static static_params_struc static_params;

static Err set_static_params_for_calib(
    String           swp_yc_name,
    SrtDiffusionType swp_vol_type,
    double           yr_to_exp,
    double           cms_rate,
    double           cmt_rate,
    double           atm_sigma_beta,
    double           sabr_alpha,
    double           sabr_beta,
    double           sabr_rho,
    double           spot_spread,
    double           fwd_spread,
    double*          cmt_calibration_strikes,
    long             num_strike)

{
    static_params.swp_yc_name             = swp_yc_name;
    static_params.swp_vol_type            = swp_vol_type;
    static_params.yr_to_exp               = yr_to_exp;
    static_params.cms_rate                = cms_rate;
    static_params.cmt_rate                = cmt_rate;
    static_params.atm_sigma_beta          = atm_sigma_beta;
    static_params.sabr_alpha              = sabr_alpha;
    static_params.sabr_beta               = sabr_beta;
    static_params.sabr_rho                = sabr_rho;
    static_params.spot_spread             = spot_spread;
    static_params.fwd_spread              = fwd_spread;
    static_params.cmt_calibration_strikes = cmt_calibration_strikes;
    static_params.num_strike              = num_strike;

    static_params.ppdParamBounds = dmatrix(1, 2, 1, 3);

    return NULL;
}

static Err set_CMT_VOL_parameters_bounds(void)
{
    static_params.ppdParamBounds[1][1] = MIN_SPREAD_MEAN_REVERSION;
    static_params.ppdParamBounds[2][1] = MAX_SPREAD_MEAN_REVERSION;

    static_params.ppdParamBounds[1][2] = MIN_SPREAD_LOC_VOL;
    static_params.ppdParamBounds[2][2] = MAX_SPREAD_LOC_VOL;

    static_params.ppdParamBounds[1][3] = MIN_RHO;
    static_params.ppdParamBounds[2][3] = MAX_RHO;

    return NULL;
}

static Err from_optparam_to_CMT_VOL_parameters(
    double* pdOptimParams, double* spread_mean_reversion, double* spread_loc_vol, double* rho)
{
    (*spread_mean_reversion) = pdOptimParams[1];
    (*spread_loc_vol)        = pdOptimParams[2];
    (*rho)                   = pdOptimParams[3];

    return NULL;
}

static Err levenberg_CMT_VOL_funcs(
    double instr_index, double* pdOptimParams, double* value, int lNumParams)
{
    Err              err          = NULL;
    String           swp_yc_name  = static_params.swp_yc_name;
    SrtDiffusionType swp_vol_type = static_params.swp_vol_type;
    double           yr_to_exp    = static_params.yr_to_exp;
    double atm_sigma_beta = static_params.atm_sigma_beta, sabr_alpha = static_params.sabr_alpha,
           sabr_beta = static_params.sabr_beta, sabr_rho = static_params.sabr_rho,
           spot_spread              = static_params.spot_spread,
           *cmt_calibration_strikes = static_params.cmt_calibration_strikes,
           cms_rate = static_params.cms_rate, cmt_rate = static_params.cmt_rate;

    double spread_mean_reversion, spread_tau;
    double spread_fwd_level = 0.0;
    double spread_loc_vol;
    double rho;
    long   num_strike = static_params.num_strike;
    double cmt_vol;

    if (are_calib_parameters_within_band(pdOptimParams, lNumParams, static_params.ppdParamBounds) ==
        SRT_NO)
    {
        *value = SRT_BIG;
        return NULL;
    }

    err = from_optparam_to_CMT_VOL_parameters(
        pdOptimParams, &spread_mean_reversion, &spread_loc_vol, &rho);
    if (err)
        return err;

    spread_tau = 1.0 / spread_mean_reversion;

    err = srt_f_get_cmt_vol(
        swp_yc_name,
        swp_vol_type,
        yr_to_exp,
        cms_rate,
        cmt_rate,
        cmt_calibration_strikes[(long)instr_index],

        atm_sigma_beta,
        sabr_alpha,
        sabr_beta,
        sabr_rho,

        spot_spread,
        spread_tau,
        spread_fwd_level,
        spread_loc_vol,
        rho,
        SRT_TRUE,
        &cmt_vol);

    (*value) = cmt_vol * cmt_vol;

    if (err)
        return err;

    return NULL;
}

static Err levenberg_calib_funcs(
    double instr_index, double pdOptParams[], double* value, double deriv[], int lNumParams)
{
    Err    err = NULL;
    double shift;
    long   i;

    err = levenberg_CMT_VOL_funcs(instr_index, pdOptParams, value, lNumParams);
    if (err)
        return err;

    for (i = 1; i <= lNumParams; i++)
    {
        if (pdOptParams[i] == 0.0)
            shift = ZERO_SHIFT;
        else
            shift = PROP_SHIFT * pdOptParams[i];
        pdOptParams[i] += shift;

        err = levenberg_CMT_VOL_funcs(instr_index, pdOptParams, &(deriv[i]), lNumParams);
        if (err)
            return err;

        deriv[i] -= *value;
        deriv[i] /= shift;

        pdOptParams[i] -= shift;
    }

    return NULL;
}

Err srt_f_calibrate_CMT_CMS_spread_parameters(
    double yr_to_exp,
    String swp_yc_name,
    String swp_vol_type_str,
    double atm_sigma_beta,
    double sabr_alpha,
    double sabr_beta,
    double sabr_rho,

    double cmt_rate,
    double spot_spread,
    double fwd_spread,

    double spread_mean_reversion_guess,
    double spread_fwd_level_guess,
    double spread_loc_vol_guess,
    double rho_guess,

    long    num_strike,
    double* cmt_calibration_strikes,
    double* cmt_vols,

    double** param)
{
    Err              err = NULL;
    long             ndata, i;
    double           cms_rate;
    SrtDiffusionType swp_vol_type;
    double *         data, *target, *weight, *opt_params;
    double           chisq;

    err = interp_diffusion_type(swp_vol_type_str, &swp_vol_type);
    if (err)
        return err;

    cms_rate = (cmt_rate + spot_spread);

    err = set_static_params_for_calib(
        swp_yc_name,
        swp_vol_type,
        yr_to_exp,
        cms_rate,
        cmt_rate,
        atm_sigma_beta,
        sabr_alpha,
        sabr_beta,
        sabr_rho,
        spot_spread,
        fwd_spread,
        cmt_calibration_strikes,
        num_strike); /* NBR OF PARAMETER */
    if (err)
        return err;

    err = set_CMT_VOL_parameters_bounds();
    if (err)
        return err;

    ndata = num_strike;

    data       = dvector(1, ndata);
    target     = dvector(1, ndata);
    weight     = dvector(1, ndata);
    opt_params = dvector(1, 3);

    for (i = 1; i <= ndata; i++)
        data[i] = (double)(i);
    for (i = 1; i <= ndata; i++)
        target[i] = (cmt_vols[i] * cmt_vols[i]);
    for (i = 1; i <= ndata; i++)
        weight[i] = target[i];

    opt_params[1] = spread_mean_reversion_guess;
    opt_params[2] = spread_loc_vol_guess;
    opt_params[3] = rho_guess;

    err = levenberg_marquardt(
        data, target, weight, ndata, opt_params, 3, 100, &levenberg_calib_funcs, &chisq);
    if (err)
        return err;

    (*param)[1] = opt_params[1];
    (*param)[2] = (fwd_spread - spot_spread * trunc_exp(-opt_params[1] * yr_to_exp)) /
                  (1 - trunc_exp(-opt_params[1] * yr_to_exp));
    (*param)[3] = opt_params[2];
    (*param)[4] = opt_params[3];

    if (data)
        free_dvector(data, 1, ndata);
    data = NULL;
    if (target)
        free_dvector(target, 1, ndata);
    target = NULL;
    if (weight)
        free_dvector(weight, 1, ndata);
    weight = NULL;
    if (opt_params)
        free_dvector(opt_params, 1, 3);
    opt_params = NULL;
    if (static_params.ppdParamBounds)
        free_dmatrix(static_params.ppdParamBounds, 1, 2, 1, 3);
    static_params.ppdParamBounds = NULL;

    return NULL;
}

/* COMPUTE THE NORMAL CMS VOLATILTY FROM THE CMS SABR PARAMATERS AND ADJUST IT BY THE CONVEXITY OF
 * THE SPREAD OVER CMT */
Err srt_f_get_cmt_vol(
    String           swp_yc_name,
    SrtDiffusionType swp_vol_type,

    double yr_to_exp,
    double cms_rate,
    double cmt_rate,
    double cmt_strike,

    /* SABR CMS PARAMETERS */
    double atm_sigma_beta,
    double sabr_alpha,
    double sabr_beta,
    double sabr_rho,

    /* SPREAD DYNAMICS PARAMETERS */
    double      spot_spread,      /* THE SPREAD INITIAL CONDITION */
    double      spread_tau,       /* ONE OVER THE MEAN REVERSION */
    double      spread_fwd_level, /* THE MEAN REVERSION LEVEL */
    double      spread_loc_vol,   /* THE SPREAD VOLATILITY */
    double      rho,              /* CORRELATION BETWEEN CMS AND SPREAD */
    SRT_Boolean use_mkt_fwd_spread,

    double* cmt_vol)
{
    Err    err = NULL;
    double cms_adjusted_strike, fwd_spread, cms_norm_vol, cms_market_vol;
    double spread_vol, cmt_norm_vol;

    /* COMPUTE THE FWD SPREAD AND ADJUST THE CMT STRIKE BY THIS SPREAD */

    if (use_mkt_fwd_spread == SRT_FALSE)
        fwd_spread = spread_fwd_level +
                     (spot_spread - spread_fwd_level) * trunc_exp(-yr_to_exp / spread_tau);
    else
        fwd_spread = static_params.fwd_spread;

    cms_adjusted_strike = cmt_strike + fwd_spread;

    /* GET THE CMS MARKET VOL :SABR IMPLIED BLACK AND SCHOLES VOLATILITY */
    err = srt_f_optsarbvol(
        cms_rate,
        cms_adjusted_strike,
        yr_to_exp,
        atm_sigma_beta,
        sabr_alpha,
        sabr_beta,
        sabr_rho,
        swp_vol_type,
        swp_vol_type,
        &cms_market_vol);
    if (err)
        return err;

    /* CONVERT THE CMS VOLATILITY  INTO NORMAL */
    err = srt_f_optsarbvol(
        cms_rate,
        cms_adjusted_strike,
        yr_to_exp,
        cms_market_vol,
        sabr_alpha,
        sabr_beta,
        sabr_rho,
        swp_vol_type,
        SRT_NORMAL,
        &cms_norm_vol);
    if (err)
        return err;

    /* COMPUTE THE NORMAL VOLATILITY OF THE CMT */
    spread_vol =
        spread_loc_vol *
        sqrt(0.5 * (1 - trunc_exp(-2 * yr_to_exp / spread_tau)) * (spread_tau / yr_to_exp));
    cmt_norm_vol = sqrt(
        cms_norm_vol * cms_norm_vol - 2 * rho * cms_norm_vol * spread_vol +
        spread_vol * spread_vol);

    /* CONVERT THE CMT NORMAL VOLATILTY INTO swp_vol_type TYPE */
    if (swp_vol_type == SRT_NORMAL)
        (*cmt_vol) = cmt_norm_vol;
    else
        (*cmt_vol) = (cmt_norm_vol / cmt_rate);

    return err;
}

#undef PROP_SHIFT
#undef ZERO_SHIFT
#undef MIN_SPREAD_MEAN_REVERSION
#undef MAX_SPREAD_MEAN_REVERSION
#undef MIN_SPREAD_FWD_LEVEL
#undef MAX_SPREAD_FWD_LEVEL
#undef MIN_SPREAD_LOC_VOL
#undef MAX_SPREAD_LOC_VOL
#undef MAX_RHO
#undef MIN_RHO
#undef SRT_BIG
#undef MAX_EXP_ARG
#undef MIN_EXP_ARG
