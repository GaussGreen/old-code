
#include "LGMSVCalibApprox.h"

#include "CTSProdStruct.h"
#include "DiagCalibDLM.h"
#include "DiagCalibDLMSV.h"
#include "DiagCalibGen.h"
#include "FXSABRAdi.h"
#include "Fx3FUtils.h"
#include "LGMSVCalibGen.h"
#include "LGMSVClosedForm.h"
#include "LGMSVClosedFormApprox.h"
#include "LGMSVPDE.h"
#include "math.h"
#include "opfnctns.h"
#include "swp_h_vol.h"

#define MAX_CPN 600
#define ONE_MONTH 0.083333333
#define MAXJUMP_LGMSV 0.5

#define EPS_MAT_TO_LONG_DATE_CONVERSION 1.0e-8

#define LIMIT_LAM_SV -0.30

#define LGMSV_DEFAULT_LAMEPS 0.05
#define LGMSV_DEFAULT_ALPHAEPS 0.30
#define LGMSV_DEFAULT_RHOEPS 0.30
#define LGMSV_ALPHA_MAXITER 20
#define LGMSV_SABRBETA 0.0

/* Set default Calib Params */
void LGMSV_SetDefault_CalibParams(LGMSV_CALIBPARAMS params)
{
    params->fix_lambda       = 1;
    params->use_lgm_lambda   = 0;
    params->calib_flat_smile = 1;

    params->use_sabr_calib            = 1;
    params->use_sabr_levenberg        = 0;
    params->sabr_calib_min_time       = 3.0;
    params->sabr_calib_max_lambda     = 0.25;
    params->sabr_calib_default_lambda = 0.05;

    params->max_calib_alpha = 1.0;

    params->calib_rr_bt  = 1;
    params->calib_alpha  = 0;
    params->calib_lameps = 0;
    params->calib_rho    = 0;
    params->calib_rho2   = 0;

    params->calib_smile_on_prim       = 1;
    params->novolcalib_on_smile_calib = 0;
    params->pricerho_on_alpha         = 0;

    params->onefac_rho = 1;
    params->rho_mat1   = 0.5;
    params->rho_mat2   = 3.0;

    params->use_lgm_first_guess = 1;

    params->recalib_at_end_lambda = 1;
    params->recalib_at_end_alpha  = 1;
    params->recalib_at_end_lameps = 1;
    params->recalib_at_end_rho    = 1;

    params->alpha_sv_shift = 0.0;
    params->lam_sv_shift   = 0.0;
    params->rho_sv_shift   = 0.0;
    params->rho2_sv_shift  = 0.0;
}

void LGMSV_Copy_CalibParams(LGMSV_CALIBPARAMS source, LGMSV_CALIBPARAMS target)
{
    target->fix_lambda                = source->fix_lambda;
    target->use_lgm_lambda            = source->use_lgm_lambda;
    target->calib_flat_smile          = source->calib_flat_smile;
    target->use_sabr_calib            = source->use_sabr_calib;
    target->use_sabr_levenberg        = source->use_sabr_levenberg;
    target->sabr_calib_min_time       = source->sabr_calib_min_time;
    target->sabr_calib_max_lambda     = source->sabr_calib_max_lambda;
    target->sabr_calib_default_lambda = source->sabr_calib_default_lambda;
    target->max_calib_alpha           = source->max_calib_alpha;

    target->calib_rr_bt  = source->calib_rr_bt;
    target->calib_alpha  = source->calib_alpha;
    target->calib_lameps = source->calib_lameps;
    target->calib_rho    = source->calib_rho;
    target->calib_rho2   = source->calib_rho2;

    target->calib_smile_on_prim       = source->calib_smile_on_prim;
    target->novolcalib_on_smile_calib = source->novolcalib_on_smile_calib;
    target->pricerho_on_alpha         = source->pricerho_on_alpha;

    target->onefac_rho          = source->onefac_rho;
    target->rho_mat1            = source->rho_mat1;
    target->rho_mat2            = source->rho_mat2;
    target->use_lgm_first_guess = source->use_lgm_first_guess;

    target->recalib_at_end_lambda = source->recalib_at_end_lambda;
    target->recalib_at_end_alpha  = source->recalib_at_end_alpha;
    target->recalib_at_end_lameps = source->recalib_at_end_lameps;
    target->recalib_at_end_rho    = source->recalib_at_end_rho;

    target->alpha_sv_shift = source->alpha_sv_shift;
    target->lam_sv_shift   = source->lam_sv_shift;
    target->rho_sv_shift   = source->rho_sv_shift;
    target->rho2_sv_shift  = source->rho2_sv_shift;
}

Err LGMSV_Get_Rho_From_RhoTarget(
    double  rho_target,
    double  rho_mat1,
    double  rho_mat2,
    double  LGM_alpha,
    double  LGM_gamma,
    double  LGM_rho,
    double* rho1,
    double* rho2)
{
    Err err = NULL;

    double s1, s2;

    s1 = 1.0 + LGM_alpha * LGM_alpha * exp(-2.0 * LGM_gamma * rho_mat1) +
         2.0 * LGM_rho * LGM_alpha * exp(-LGM_gamma * rho_mat1);

    if (s1 < 0.0)
    {
        err = "Cannot convert rho Target in LGMSV 2 Factor";
        return err;
    }

    s1 = sqrt(s1);

    s2 = 1.0 + LGM_alpha * LGM_alpha * exp(-2.0 * LGM_gamma * rho_mat2) +
         2.0 * LGM_rho * LGM_alpha * exp(-LGM_gamma * rho_mat2);

    if (s2 < 0.0)
    {
        err = "Cannot convert rho Target in LGMSV 2 Factor";
        return err;
    }

    s2 = sqrt(s2);

    *rho2 = (s1 - s2) * rho_target / LGM_alpha /
            (exp(-LGM_gamma * rho_mat1) - exp(-LGM_gamma * rho_mat2));
    *rho1 = s2 * rho_target - LGM_alpha * (*rho2) * exp(-LGM_gamma * rho_mat2);

    return err;
}

Err LGMSV_Get_RhoTarget_From_Rho(
    double  rho1,
    double  rho2,
    double  rho_mat1,
    double  rho_mat2,
    double  LGM_alpha,
    double  LGM_gamma,
    double  LGM_rho,
    double* rho_target)
{
    Err    err = NULL;
    double rho_target1, rho_target2;

    double s1, s2;

    s1 = 1.0 + LGM_alpha * LGM_alpha * exp(-2.0 * LGM_gamma * rho_mat1) +
         2.0 * LGM_rho * LGM_alpha * exp(-LGM_gamma * rho_mat1);

    if (s1 < 0.0)
    {
        err = "Cannot convert rho Target in LGMSV 2 Factor";
        return err;
    }

    s1 = sqrt(s1);

    s2 = 1.0 + LGM_alpha * LGM_alpha * exp(-2.0 * LGM_gamma * rho_mat2) +
         2.0 * LGM_rho * LGM_alpha * exp(-LGM_gamma * rho_mat2);

    if (s2 < 0.0)
    {
        err = "Cannot convert rho Target in LGMSV 2 Factor";
        return err;
    }

    s2 = sqrt(s2);

    rho_target1 = rho1 * (exp(LGM_gamma * rho_mat1) - exp(LGM_gamma * rho_mat2)) /
                  (exp(LGM_gamma * rho_mat1) * s1 - exp(LGM_gamma * rho_mat2) * s2);
    rho_target2 =
        rho2 * LGM_alpha * (exp(-LGM_gamma * rho_mat1) - exp(-LGM_gamma * rho_mat2)) / (s1 - s2);

    *rho_target = 0.5 * (rho_target1 + rho_target2);

    return err;
}

Err LGMSV_Get_RhoTargetSensitivity_From_Rho(
    double  rho1,
    double  rho2,
    double  rho_mat1,
    double  rho_mat2,
    double  LGM_alpha,
    double  LGM_gamma,
    double  LGM_rho,
    double* sensi1,
    double* sensi2)
{
    Err err = NULL;

    double s1, s2;

    s1 = 1.0 + LGM_alpha * LGM_alpha * exp(-2.0 * LGM_gamma * rho_mat1) +
         2.0 * LGM_rho * LGM_alpha * exp(-LGM_gamma * rho_mat1);

    if (s1 < 0.0)
    {
        err = "Cannot convert rho Target in LGMSV 2 Factor";
        return err;
    }

    s1 = sqrt(s1);

    s2 = 1.0 + LGM_alpha * LGM_alpha * exp(-2.0 * LGM_gamma * rho_mat2) +
         2.0 * LGM_rho * LGM_alpha * exp(-LGM_gamma * rho_mat2);

    if (s2 < 0.0)
    {
        err = "Cannot convert rho Target in LGMSV 2 Factor";
        return err;
    }

    s2 = sqrt(s2);

    *sensi1 = (exp(LGM_gamma * rho_mat1) - exp(LGM_gamma * rho_mat2)) /
              (exp(LGM_gamma * rho_mat1) * s1 - exp(LGM_gamma * rho_mat2) * s2);
    *sensi2 = LGM_alpha * (exp(-LGM_gamma * rho_mat1) - exp(-LGM_gamma * rho_mat2)) / (s1 - s2);

    return err;
}

Err LGMSVCalibApprox(
    /* Instrument informations	*/
    int    nex, /*	Total number of exercise dates */
    double ex_time[],
    double ex_lfwd[],
    double ex_llvl[],
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double shift[],
    double coef_vol[],
    double coef_meanrev[],

    /* First guess */
    double sig[],

    /* Model informations */
    double  dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,
    long*   lSigIndex,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Newton parameters */
    double Precision,
    int    NbIterMax,

    /* Output */
    double* dSigmaTS)
{
    double** res_iter = NULL;
    long     new_index, last_index;
    double   last_time, cum_var, cum_var_temp;
    double   price1, price2, sig1, sig2;
    int      i, j, k, l;
    double   dShiftStrike, dShiftFwd;
    Err      err = NULL;

    last_index = -1;
    last_time  = 0.0;
    cum_var    = 0.0;
    cum_var_temp;

    res_iter = dmatrix(0, NbIterMax, 0, 1);

    if (!res_iter)
    {
        err = "Memory allocation faillure in LGMSVCalibApprox";
        goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++)
    {
        new_index    = lSigIndex[i];
        sig1         = (sig)[i];
        dShiftFwd    = ex_lfwd[i] + shift[i];
        dShiftStrike = ex_lstrike[i] + shift[i];

        /* update sigma */
        for (j = last_index + 1; j <= new_index; j++)
        {
            dSigmaTS[j] = sig1;
        }

        cum_var_temp = cum_var + sig1 * sig1 * (ex_time[i] - last_time);

        /* first pricing */
        LGMSVClosedFormApprox(
            dLambdaX,
            iNbPWTime,
            dPWTime,
            dSigmaTS,
            dAlphaTS,
            dLambdaEpsTS,
            dRhoTS,
            dTStar,
            dShiftFwd,
            1,
            &dShiftStrike,
            ex_time[i],
            coef_vol[i],
            1000,
            coef_meanrev[i],
            0,
            0,
            0,
            NumerParams,
            &price1);

        if (NumerParams->iIntegMethod == 2)
            if (1.0 > dShiftStrike)
                price1 *= (ex_lfwd[i] + shift[i]);
            else
                price1 = (price1 + dShiftStrike - 1.0) * (ex_lfwd[i] + shift[i]);

        price1 *= ex_llvl[i];

        if (fabs(ex_lprice[i] - price1) > Precision)
        {
            /* We do a Newton */
            res_iter[0][0] = sig1;
            res_iter[0][1] = price1;

            sig2 = sqrt(
                (ex_lprice[i] * ex_lprice[i] / price1 / price1 * cum_var_temp - cum_var) /
                (ex_time[i] - last_time));

            k = 1;

            while (fabs(ex_lprice[i] - price1) > Precision && k < NbIterMax)
            {
                k++;

                /* update sigma */
                for (j = last_index + 1; j <= new_index; j++)
                {
                    dSigmaTS[j] = sig2;
                }

                LGMSVClosedFormApprox(
                    dLambdaX,
                    iNbPWTime,
                    dPWTime,
                    dSigmaTS,
                    dAlphaTS,
                    dLambdaEpsTS,
                    dRhoTS,
                    dTStar,
                    dShiftFwd,
                    1,
                    &dShiftStrike,
                    ex_time[i],
                    coef_vol[i],
                    1000,
                    coef_meanrev[i],
                    0,
                    0,
                    0,
                    NumerParams,
                    &price2);

                if (NumerParams->iIntegMethod == 2)
                    if (1.0 > dShiftStrike)
                        price2 *= (ex_lfwd[i] + shift[i]);
                    else
                        price2 = (price2 + dShiftStrike - 1.0) * (ex_lfwd[i] + shift[i]);

                price2 *= ex_llvl[i];

                /* Save Res */
                l = 0;
                while (l < k - 1 && res_iter[l][0] < sig2)
                {
                    l++;
                }

                if (l < k - 1)
                {
                    for (j = k - 2; j >= l; j--)
                    {
                        res_iter[j + 1][0] = res_iter[j][0];
                        res_iter[j + 1][1] = res_iter[j][1];
                    }
                }

                res_iter[l][0] = sig2;
                res_iter[l][1] = price2;

                sig1   = sig2;
                price1 = price2;

                sig2 = solve_for_next_coef(res_iter, k, ex_lprice[i], 2);
            }
        }
        else
        {
            sig2 = sig1;
        }

        cum_var += sig2 * sig2 * (ex_time[i] - last_time);
        last_time = ex_time[i];

        /* update sigma */
        for (j = last_index + 1; j <= new_index; j++)
        {
            dSigmaTS[j] = sig2;
        }

        (sig)[i] = sig2;

        last_index = new_index;
    }

    for (j = last_index + 1; j < iNbPWTime; j++)
    {
        dSigmaTS[j] = dSigmaTS[j - 1];
    }

FREE_RETURN:

    if (res_iter)
        free_dmatrix(res_iter, 0, NbIterMax, 0, 1);

    return err;
}

Err cpd_calibSV_approx(
    /*	Market */
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
    char* instr_freq, /*	Frequency and basis of instruments */
    char* instr_basis,
    /*	If ex_date is NULL,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates, /*	Exercise dates,
                                                      all supposed to be on or after today */
    long*  ex_date,   /*	Supposed to be sorted */
    char** end_tenor, /*	Tenors of the underlying instruments
                                                              or "DIAG" */
    long    end_date, /*	End date for diagonal */
    double* strike,   /*	Strikes
                                              0: ATM */
    /*	Model Parameters */
    double  dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Newton parameters */
    double Precision,
    int    NbIterMax,
    /*	Output */
    int*      numres,
    double**  restime,
    double*** result,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data) /*	NULL = don't save calibration instrument data */
{
    int            i, j, k, l, nex, ncpn;
    SrtCompounding ifreq;
    SrtBasisCode   ibasis;
    double         ex_time[MAX_CPN], ex_lfwd[MAX_CPN], ex_lswp[MAX_CPN], ex_llvl[MAX_CPN],
        ex_lstrike[MAX_CPN], ex_lvol[MAX_CPN], ex_lprice[MAX_CPN], atm_price[MAX_CPN],
        ex_zeta[MAX_CPN];
    int         ex_cpn[MAX_CPN], ex_endcpn[MAX_CPN];
    long        cpn_date[MAX_CPN];
    double      cpn_time[MAX_CPN], cpn_cvg[MAX_CPN], cpn_df[MAX_CPN];
    long        tmplng1[MAX_CPN], tmplng2[MAX_CPN];
    long *      theo_end_dates, *act_end_dates;
    long        theo_date, act_date, temp_date, temp_date2;
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      std;
    SrtCurvePtr yc_ptr;
    double      atm_vol, fact;

    double coef_vol[MAX_CPN], coef_meanrev[MAX_CPN], shift[MAX_CPN];

    int     iNbPWTimeNew, num_sig;
    double *sig = NULL, *sig_time = NULL, *dPWTimeNew = NULL, *dSigmaTS = NULL, *dAlphaTSNew = NULL,
           *dRhoTSNew = NULL, *dLambdaEpsTSNew = NULL;

    long* lSigIndex = NULL;

    long last_index;

    double sumAlpha2, sumLambda, sumRho, AlphaEq, LambdaEq, RhoEq, temp;

    Err err = NULL;

    theo_end_dates = &(tmplng1[0]);
    act_end_dates  = &(tmplng2[0]);

    sig_time = NULL;
    sig      = NULL;

    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    /*	1.)	Setup the bond schedule and its coupons */

    /*	Coupons */

    err = interp_compounding(instr_freq, &ifreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = interp_basis(instr_basis, &ibasis);
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
        strike++;
        num_ex_dates--;
        if (num_ex_dates == 0)
        {
            err = "All exercise dates are past in cpd_calib_diagonal";
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
                strike[k + 1]         = strike[k];
            }

            ex_date++;
            end_tenor++;
            theo_end_dates++;
            act_end_dates++;
            strike++;
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

    /*	Remove close dates */
    j = num_ex_dates - 1;
    for (i = num_ex_dates - 2; i >= 0; i--)
    {
        if ((ex_date[j] - ex_date[i]) * YEARS_IN_DAY < param->min_time - ONE_MONTH)
        {
            for (k = i - 1; k >= 0; k--)
            {
                ex_date[k + 1] = ex_date[k];
                strcpy(end_tenor[k + 1], end_tenor[k]);
                theo_end_dates[k + 1] = theo_end_dates[k];
                act_end_dates[k + 1]  = act_end_dates[k];
                strike[k + 1]         = strike[k];
            }

            ex_date++;
            end_tenor++;
            theo_end_dates++;
            act_end_dates++;
            strike++;
            num_ex_dates--;
            j--;
            if (num_ex_dates < 1)
            {
                err = "All exercise dates are past in cpd_calib_diagonal";
                goto FREE_RETURN;
            }
        }
        else
        {
            j = i;
        }
    }

    /*	Remove last? */
    if (param->skip_last && num_ex_dates > 1)
    {
        num_ex_dates--;
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
        std += (param->shift_type == 1 ? std * param->vol_shift : param->vol_shift);
        if (power > 0.5)
        {
            power = srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, ex_time[i], 1.0, SRT_PUT, PREMIUM);
            err   = srt_f_optimpvol(
                power, ex_lfwd[i], ex_lfwd[i], ex_time[i], 1.0, SRT_PUT, SRT_NORMAL, &std);
        }
        std *= sqrt(ex_time[i]);

        /*	Strike */
        if (param->strike_type == 0 || strike[i] < 1.0e-04)
        {
            ex_lstrike[i] = ex_lfwd[i];
        }
        else if (param->strike_type == 1)
        {
            ex_lstrike[i] = strike[i];
        }
        else if (param->strike_type == 2)
        {
            if (err = swp_f_ForwardRate(
                    cpn_date[j],
                    theo_end_dates[i],
                    instr_freq,
                    instr_basis,
                    yc_name,
                    ref_rate_name,
                    &swp_rte))
            {
                goto FREE_RETURN;
            }

            spr = swp_rte - ex_lfwd[i];

            ex_lstrike[i] = strike[i] - spr;
        }
        else if (param->strike_type == 3)
        {
            ex_lstrike[i] = ex_lfwd[i] + strike[i] * std;
        }

        /*	Apply max std */
        if (ex_lstrike[i] > ex_lfwd[i] + param->max_std * std)
        {
            ex_lstrike[i] = ex_lfwd[i] + param->max_std * std;
        }
        else if (ex_lstrike[i] < ex_lfwd[i] - param->max_std * std)
        {
            ex_lstrike[i] = ex_lfwd[i] - param->max_std * std;
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
        ex_lvol[i] += (param->shift_type == 1 ? ex_lvol[i] * param->vol_shift : param->vol_shift);

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

    num_sig  = nex;
    sig_time = (double*)calloc(nex, sizeof(double));
    sig      = (double*)calloc(nex, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++)
    {
        (sig_time)[i] = ex_time[i];
    }

    /*	Merge the term structures */

    dPWTimeNew   = (double*)calloc(num_ex_dates, sizeof(double));
    iNbPWTimeNew = num_ex_dates;

    if (!(sig_time) || !(sig) || !dPWTime)
    {
        err = "Memory allocation faillure in cpd_calibSV_approx";
        goto FREE_RETURN;
    }

    for (i = 0; i < num_ex_dates; i++)
    {
        dPWTimeNew[i] = (sig_time)[i];
    }

    num_f_concat_vector(&iNbPWTimeNew, &dPWTimeNew, iNbPWTime, dPWTime);
    num_f_sort_vector(iNbPWTimeNew, dPWTimeNew);
    num_f_unique_vector(&iNbPWTimeNew, dPWTimeNew);

    dSigmaTS        = dvector(0, iNbPWTimeNew - 1);
    dAlphaTSNew     = dvector(0, iNbPWTimeNew - 1);
    dLambdaEpsTSNew = dvector(0, iNbPWTimeNew - 1);
    dRhoTSNew       = dvector(0, iNbPWTimeNew - 1);
    lSigIndex       = lvector(0, num_ex_dates - 1);

    if (!dSigmaTS || !dAlphaTSNew || !dLambdaEpsTSNew || !dRhoTSNew || !lSigIndex)
    {
        err = "Memory allocation faillure (2) in cpd_calibSV_approx";
        goto FREE_RETURN;
    }

    for (i = 0; i < iNbPWTimeNew; i++)
    {
        j                  = Get_Index(dPWTimeNew[i], dPWTime, iNbPWTime);
        dAlphaTSNew[i]     = 2.0 * dAlphaTS[j];
        dLambdaEpsTSNew[i] = 2.0 * dLambdaEpsTS[j];
        dRhoTSNew[i]       = dRhoTS[j];
    }

    for (i = 0; i < num_ex_dates; i++)
    {
        lSigIndex[i] = Get_Index((sig_time)[i], dPWTimeNew, iNbPWTimeNew);
    }

    /* 3 ) Computation of the first guess */

    /* First find the equivalent ATM volatilities */

    last_index = -1;
    sumAlpha2  = 0.0;
    sumLambda  = 0.0;
    sumRho     = 0.0;

    if (power > 0.5)
    {
        for (i = 0; i < nex; i++)
        {
            for (j = last_index + 1; j <= lSigIndex[i]; j++)
            {
                if (j > 0)
                {
                    sumAlpha2 +=
                        dAlphaTSNew[j] * dAlphaTSNew[j] * (dPWTimeNew[j] - dPWTimeNew[j - 1]);
                    sumLambda += dLambdaEpsTSNew[j] * (dPWTimeNew[j] - dPWTimeNew[j - 1]);
                    sumRho += dRhoTSNew[j] * (dPWTimeNew[j] - dPWTimeNew[j - 1]);
                }
                else
                {
                    sumAlpha2 = dAlphaTSNew[0] * dAlphaTSNew[0] * dPWTimeNew[0];
                    sumLambda = dLambdaEpsTSNew[0] * dPWTimeNew[0];
                    sumRho    = dRhoTSNew[0] * dPWTimeNew[0];
                }
            }

            last_index = lSigIndex[i];

            AlphaEq  = sqrt(sumAlpha2 / (sig_time)[i]) / 2.0;
            LambdaEq = sumLambda / (sig_time)[i] / 2.0;
            RhoEq    = sumRho / (sig_time)[i];

            if (LambdaEq > 1.0E-16)
            {
                AlphaEq *=
                    sqrt((1.0 - exp(-2.0 * LambdaEq * ex_time[i])) / (2.0 * LambdaEq * ex_time[i]));
            }

            vol_conv(
                ex_lvol[i],
                SABR_STR_LOG,
                &atm_vol,
                SABR_ATM_LOG,
                ex_lfwd[i],
                ex_lstrike[i],
                ex_time[i],
                AlphaEq,
                0.0,
                RhoEq);

            atm_price[i] = srt_f_optblksch(
                ex_lfwd[i], ex_lfwd[i], atm_vol, ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
        }
    }
    else
    {
        for (i = 0; i < nex; i++)
        {
            for (j = last_index + 1; j <= lSigIndex[i]; j++)
            {
                if (j > 0)
                {
                    sumAlpha2 +=
                        dAlphaTSNew[j] * dAlphaTSNew[j] * (dPWTimeNew[j] - dPWTimeNew[j - 1]);
                    sumLambda += dLambdaEpsTSNew[j] * (dPWTimeNew[j] - dPWTimeNew[j - 1]);
                    sumRho += dRhoTSNew[j] * (dPWTimeNew[j] - dPWTimeNew[j - 1]);
                }
                else
                {
                    sumAlpha2 = dAlphaTSNew[0] * dAlphaTSNew[0] * dPWTimeNew[0];
                    sumLambda = dLambdaEpsTSNew[0] * dPWTimeNew[0];
                    sumRho    = dRhoTSNew[0] * dPWTimeNew[0];
                }
            }

            last_index = lSigIndex[i];

            AlphaEq  = sqrt(sumAlpha2 / (sig_time)[i]) / 2.0;
            LambdaEq = sumLambda / (sig_time)[i] / 2.0;
            RhoEq    = sumRho / (sig_time)[i];

            if (LambdaEq > 1.0E-16)
            {
                AlphaEq *=
                    sqrt((1.0 - exp(-2.0 * LambdaEq * ex_time[i])) / (2.0 * LambdaEq * ex_time[i]));
            }

            vol_conv(
                ex_lvol[i],
                SABR_STR_NORM,
                &atm_vol,
                SABR_ATM_NORM,
                ex_lfwd[i],
                ex_lstrike[i],
                ex_time[i],
                AlphaEq,
                0.001,
                RhoEq);

            atm_price[i] = srt_f_optblknrm(
                ex_lfwd[i], ex_lfwd[i], atm_vol, ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
        }
    }

    err = lgmprcapgivenlambda_2(
        ncpn,
        cpn_time,
        cpn_df,
        cpn_cvg,
        nex,
        ex_time,
        ex_cpn,
        ex_endcpn,
        ex_lfwd,
        atm_price,
        ex_zeta,
        dLambdaX,
        1,
        0.0,
        0.0,
        0.0,
        0);
    if (err)
    {
        goto FREE_RETURN;
    }

    fact = exp(dLambdaX * dTStar);

    (sig)[0] = sqrt(ex_zeta[0]) / fact;

    for (i = 1; i < nex; i++)
    {
        (sig_time)[i] = ex_time[i];
        if (ex_zeta[i] > ex_zeta[i - 1])
        {
            (sig)[i] = sqrt(ex_zeta[i] - ex_zeta[i - 1]) / fact;
        }
        else
        {
            smessage(
                "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                ex_time[i]);
            for (j = i; j < nex; j++)
            {
                (sig)[j] = (sig)[i - 1];
            }
            i = nex;
        }
    }

    /*	Compute the coefs */
    for (i = 0; i < nex; i++)
    {
        j = ex_cpn[i];
        l = ex_endcpn[i];

        err = Calculate_HestonEquivalent_LGMSV(
            l - j + 1,
            &(cpn_date[j]),
            &(cpn_time[j]),
            &(cpn_cvg[j]),
            &(cpn_df[j]),
            0,
            0,
            today,
            yc_name,
            1.0 / dLambdaX,
            dTStar,
            &swp_rte,
            &lvl,
            &(shift[i]),
            &(coef_vol[i]),
            &(coef_meanrev[i]),
            &temp);

        if (err = swp_f_ForwardRate(
                cpn_date[j],
                theo_end_dates[i],
                instr_freq,
                instr_basis,
                yc_name,
                ref_rate_name,
                &(ex_lswp[i])))
        {
            goto FREE_RETURN;
        }

        spr = ex_lswp[i] - ex_lfwd[i];
        shift[i] -= spr;
    }

    err = LGMSVCalibApprox(
        nex,
        ex_time,
        ex_lswp,
        ex_llvl,
        ex_lstrike,
        ex_lprice,
        shift,
        coef_vol,
        coef_meanrev,
        sig,
        dLambdaX,
        iNbPWTimeNew,
        dPWTimeNew,
        dAlphaTSNew,
        dLambdaEpsTSNew,
        dRhoTSNew,
        dTStar,
        lSigIndex,
        NumerParams,
        Precision,
        NbIterMax,
        dSigmaTS);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	5.) Convert the TS */

    ConvertTS_LGMSV_to_LGM(iNbPWTimeNew, dPWTimeNew, dSigmaTS, dLambdaX, dTStar);

    /*	6.)	Create the result */

    *numres  = iNbPWTimeNew;
    *restime = dvector(0, *numres - 1);
    *result  = dmatrix(0, *numres - 1, 0, 3);

    if (!(*restime) || !(*result))
    {
        err = "Memory allocation faillure in LGMSVCalibApprox";
        goto FREE_RETURN;
    }

    for (i = 0; i < *numres; i++)
    {
        (*restime)[i]   = dPWTimeNew[i];
        (*result)[i][0] = dSigmaTS[i];
        (*result)[i][1] = dAlphaTSNew[i] / 2.0;
        (*result)[i][2] = dLambdaEpsTSNew[i] / 2.0;
        (*result)[i][3] = dRhoTSNew[i];
    }

    /*	6.)	Save instrument data if required */
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
            inst_data->start_dates[i]   = cpn_date[ex_cpn[i]];
            inst_data->end_dates[i]     = cpn_date[ex_endcpn[i]];
            inst_data->short_strikes[i] = 0.0;
            inst_data->long_strikes[i]  = ex_lstrike[i];
        }
    }

FREE_RETURN:

    if (err)
    {
        if (sig_time)
            free(sig_time);
        sig_time = NULL;

        if (sig)
            free(sig);
        sig = NULL;

        if (inst_data)
        {
            cpd_free_calib_inst_data(inst_data);
        }
    }

    if (dPWTimeNew)
        free(dPWTimeNew);

    if (sig)
        free(sig);
    if (sig_time)
        free(sig_time);

    if (dSigmaTS)
        free_dvector(dSigmaTS, 0, iNbPWTimeNew - 1);
    if (dAlphaTSNew)
        free_dvector(dAlphaTSNew, 0, iNbPWTimeNew - 1);
    if (dLambdaEpsTSNew)
        free_dvector(dLambdaEpsTSNew, 0, iNbPWTimeNew - 1);
    if (dRhoTSNew)
        free_dvector(dRhoTSNew, 0, iNbPWTimeNew - 1);
    if (lSigIndex)
        free_lvector(lSigIndex, 0, num_ex_dates - 1);

    return err;
}

void Fill_Equivalent_ATM_Prices(
    LGMSV_MODEL            model,
    ALLCALININSTRUMENTSDLM AllCalibInst,
    long*                  lSigIndex,
    int                    iIndexStrike,
    CALIBEXESCHEDULEDLM    CalibExeSchedule)

{
    int                i, j, last_index;
    double             atm_vol, sumAlpha2, sumLambda, sumRho, AlphaEq, LambdaEq, RhoEq;
    CALIBINSTRUMENTDLM CalibInst;

    last_index = -1;
    sumAlpha2  = 0.0;
    sumLambda  = 0.0;
    sumRho     = 0.0;

    for (i = 0; i < CalibExeSchedule->iNExe; i++)
    {
        CalibInst = &(AllCalibInst->sCalibInst[i]);

        for (j = last_index + 1; j <= lSigIndex[i]; j++)
        {
            if (j > 0)
            {
                sumAlpha2 += model->dAlpha[j] * model->dAlpha[j] *
                             (model->dPWTime[j] - model->dPWTime[j - 1]);
                sumLambda += model->dLambdaEps[j] * (model->dPWTime[j] - model->dPWTime[j - 1]);
                sumRho += model->dRho[j] * (model->dPWTime[j] - model->dPWTime[j - 1]);
            }
            else
            {
                sumAlpha2 = model->dAlpha[0] * model->dAlpha[0] * model->dPWTime[0];
                sumLambda = model->dLambdaEps[0] * model->dPWTime[0];
                sumRho    = model->dRho[0] * model->dPWTime[0];
            }
        }

        last_index = lSigIndex[i];

        AlphaEq  = sqrt(sumAlpha2 / CalibInst->dExeTime) / 2.0;
        LambdaEq = sumLambda / CalibInst->dExeTime / 2.0;
        RhoEq    = sumRho / CalibInst->dExeTime;

        if (LambdaEq > 1.0E-16)
        {
            AlphaEq *= sqrt(
                (1.0 - exp(-2.0 * LambdaEq * CalibInst->dExeTime)) /
                (2.0 * LambdaEq * CalibInst->dExeTime));
        }

        if (fabs(CalibInst->dFwdCash - CalibInst->dStrike[iIndexStrike]) > 1.0E-08)
        {
            vol_conv(
                CalibInst->dNormalVol[iIndexStrike],
                SABR_STR_NORM,
                &atm_vol,
                SABR_ATM_NORM,
                CalibInst->dFwdCash,
                CalibInst->dStrike[iIndexStrike],
                CalibInst->dExeTime,
                AlphaEq,
                0.001,
                RhoEq);
        }
        else
        {
            atm_vol = CalibInst->dNormalVol[iIndexStrike];
        }

        CalibExeSchedule->dATMPrice[i] = srt_f_optblknrm(
            CalibInst->dFwdCash,
            CalibInst->dFwdCash,
            atm_vol,
            CalibInst->dExeTime,
            CalibInst->dLevel,
            SRT_PUT,
            PREMIUM);

        CalibExeSchedule->dATMVega[i] = srt_f_optblknrm(
            CalibInst->dFwdCash,
            CalibInst->dFwdCash,
            atm_vol,
            CalibInst->dExeTime,
            CalibInst->dLevel,
            SRT_PUT,
            VEGA);
    }
}

Err LGMSV_calib_const_alpha_and_rho_from_SABR(
    AllCalibInstrumentsDLM* Instruments,
    CTS_MKT                 mkt,
    double*                 alpha,
    double*                 lambda,
    double*                 rho,
    LGMSV_CALIBPARAMS       sParams,
    CPD_CALIB_INST_DATA     inst_data)
{
    double *mkt_alpha = NULL, *mkt_rho = NULL, *exercise_times = NULL, *weights = NULL;

    double fitting_error;
    double param[3];
    int    use_param[3];
    double temp;
    double dATMLog, dAlpha, dBeta, dRho;

    int iIsSABRMarket;
    int iIsSABRAFMarket;

    int i, iIndexStart;

    Err err = NULL;

    mkt_alpha      = calloc(Instruments->iNbInst, sizeof(double));
    mkt_rho        = calloc(Instruments->iNbInst, sizeof(double));
    exercise_times = calloc(Instruments->iNbInst, sizeof(double));
    weights        = calloc(Instruments->iNbInst, sizeof(double));

    if (!mkt_alpha || !mkt_rho || !exercise_times || !weights)
    {
        err = "Memory allocation faillure (1) in LGMSV_calib_const_alpha_and_rho_from_SABR";
        goto FREE_RETURN;
    }

    /* Check if it is a SABR Market */
    /* old method
            err = swp_f_SABRvol(mkt->vc,
                                                    Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[0].iStartCpn],
                                                    Instruments->sCalibExeSchedule->lTheoEndDates[0],
                                                    0.05,
                                                    &temp,
                                                    &power,
                                                    SABR_ISMARKETSABR);
                                                    */

    err = swp_f_IsSmileVol(mkt->vc, &temp);

    if (err || (fabs(temp) < 1e-8))
    {
        iIsSABRMarket = 0;
        err           = NULL;
    }
    else
    {
        err = swp_f_SABRvol(
            mkt->vc,
            Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[0].iStartCpn],
            Instruments->sCalibExeSchedule->lTheoEndDates[0],
            0.05,
            &temp,
            &temp,
            SABR_ZETA);
        if (err)
        {
            iIsSABRMarket   = 1;
            iIsSABRAFMarket = 0;
        }
        else
        {
            iIsSABRMarket   = 0;
            iIsSABRAFMarket = 1;
        }
    }

    iIndexStart = 0;

    while (iIndexStart < Instruments->iNbInst - 1 &&
           Instruments->sCalibExeSchedule->dExeTimes[iIndexStart] < sParams->sabr_calib_min_time)
    {
        iIndexStart++;
    }

    for (i = iIndexStart; i < Instruments->iNbInst; i++)
    {
        exercise_times[i] = Instruments->sCalibInst[i].dExeTime;

        /* First guess for Alpha and Rho */
        if (iIsSABRMarket)
        {
            /* the JoeB formula */
            err = swp_f_SABRvol(
                mkt->vc,
                Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[i].iStartCpn],
                Instruments->sCalibExeSchedule->lTheoEndDates[i],
                0.05,
                &dATMLog,
                &temp,
                SABR_ATMLOG);

            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                mkt->vc,
                Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[i].iStartCpn],
                Instruments->sCalibExeSchedule->lTheoEndDates[i],
                0.05,
                &dAlpha,
                &temp,
                SABR_ALPHA);

            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                mkt->vc,
                Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[i].iStartCpn],
                Instruments->sCalibExeSchedule->lTheoEndDates[i],
                0.05,
                &dBeta,
                &temp,
                SABR_BETA);

            if (err)
                goto FREE_RETURN;

            err = swp_f_SABRvol(
                mkt->vc,
                Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[i].iStartCpn],
                Instruments->sCalibExeSchedule->lTheoEndDates[i],
                0.05,
                &dRho,
                &temp,
                SABR_RHO);

            if (err)
                goto FREE_RETURN;

            mkt_alpha[i] = dAlpha;
            mkt_rho[i]   = dRho + (dBeta - 0.0) * dATMLog / dAlpha;

            mkt_rho[i] = DMAX(mkt_rho[i], -0.99);
            mkt_rho[i] = DMIN(mkt_rho[i], 0.99);
        }
        else
        {
            if (i == iIndexStart)
            {
                mkt_alpha[i] = LGMSV_DEFAULT_ALPHAEPS;
                mkt_rho[i]   = LGMSV_DEFAULT_RHOEPS;
            }
            else
            {
                mkt_alpha[i] = mkt_alpha[i - 1];
                mkt_rho[i]   = mkt_rho[i - 1];
            }
        }

        err = cts_get_mkt_sabr_beta_param(
            mkt,
            Instruments->sCalibExeSchedule->lExeDates[i],
            Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[i].iStartCpn],
            Instruments->sCalibExeSchedule->lTheoEndDates[i],
            LGMSV_SABRBETA,
            sParams->use_sabr_levenberg,
            &mkt_alpha[i],
            &mkt_rho[i],
            &fitting_error);

        if (err)
            goto FREE_RETURN;
    }

    if (Instruments->iNbInst - iIndexStart > 1)
    {
        /* we take rho as the minimum of all the rhos */

        temp = 0.0;

        for (i = iIndexStart; i < Instruments->iNbInst; i++)
        {
            temp += mkt_rho[i];
        }

        *rho = temp / Instruments->iNbInst;

        /* we fit alpha and mean-reversion to the mkt */

        /* rescale the target */
        for (i = iIndexStart; i < Instruments->iNbInst; i++)
        {
            mkt_alpha[i] *= mkt_alpha[i];
            weights[i] = 1.0;
        }

        /* first guess */
        param[1]     = mkt_alpha[iIndexStart];
        param[2]     = 2.0 * 0.1;
        use_param[1] = 1;
        use_param[2] = 1;

        err = levenberg_marquardt_select(
            exercise_times - 1 + iIndexStart,
            mkt_alpha - 1 + iIndexStart,
            weights - 1 + iIndexStart,
            Instruments->iNbInst - iIndexStart,
            param,
            use_param,
            2,
            LGMSV_ALPHA_MAXITER,
            cts_implied_alpha_approx,
            &fitting_error);

        if (err)
            goto FREE_RETURN;

        *alpha  = sqrt(param[1]);
        *lambda = param[2] / 2.0;

        *lambda = DMIN(*lambda, sParams->sabr_calib_max_lambda);
    }
    else
    {
        if (Instruments->iNbInst - iIndexStart == 1)
        {
            *lambda = sParams->sabr_calib_default_lambda;
            *alpha  = sqrt(
                mkt_alpha[iIndexStart] * 2.0 * (*lambda) * exercise_times[iIndexStart] /
                (1.0 - exp(-2.0 * (*lambda) * exercise_times[iIndexStart])));
            *rho = mkt_rho[iIndexStart];
        }
        else
        {
            *alpha  = LGMSV_DEFAULT_ALPHAEPS;
            *lambda = sParams->sabr_calib_default_lambda;
            *rho    = LGMSV_DEFAULT_RHOEPS;
        }
    }

    if (inst_data && Instruments->iNbInst)
    {
        /* Save Smile Information */
        inst_data->num_inst_smile    = Instruments->iNbInst - iIndexStart;
        inst_data->start_dates_smile = calloc(inst_data->num_inst_smile, sizeof(long));
        inst_data->end_dates_smile   = calloc(inst_data->num_inst_smile, sizeof(long));
        inst_data->alpha             = calloc(inst_data->num_inst_smile, sizeof(double));
        inst_data->rho               = calloc(inst_data->num_inst_smile, sizeof(double));

        if (!inst_data->start_dates_smile || !inst_data->end_dates_smile || !inst_data->alpha ||
            !inst_data->rho)
        {
            err = "Memory allocation faillure in LGMSV_calib_const_alpha_and_rho_from_SABR";
            goto FREE_RETURN;
        }

        for (i = iIndexStart; i < Instruments->iNbInst; i++)
        {
            inst_data->start_dates_smile[i - iIndexStart] =
                Instruments->sCalibCpnSchedule->lCpnDate[Instruments->sCalibInst[i].iStartCpn];
            inst_data->end_dates_smile[i - iIndexStart] =
                Instruments->sCalibExeSchedule->lTheoEndDates[i];
            inst_data->alpha[i - iIndexStart] = sqrt(mkt_alpha[i]);
            inst_data->rho[i - iIndexStart]   = mkt_rho[i];
        }
    }

FREE_RETURN:

    if (mkt_alpha)
        free(mkt_alpha);
    if (mkt_rho)
        free(mkt_rho);
    if (exercise_times)
        free(exercise_times);
    if (weights)
        free(weights);

    return err;
}

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to diagonal
                with lambda calibration */
Err cpd_calib_diagonal_LGMSV_new_dlm(
    /*	Market */
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /* Get Cash Vol ref */
    char* vol_ref_rate_name,

    /* Long Instruments */
    char* instr_long_freq, /*	Frequency and basis of instruments */
    char* instr_long_basis,
    char* long_ref_rate_name, /*	Name of the reference rate */

    int     num_ex_datesl, /*	Long Exercise dates */
    long*   ex_datel_,     /*	Supposed to be sorted */
    int*    cal_datel,     /*	1: use ex_date as calibration date, 0: don't */
    char**  end_tenorl,    /*	Tenors of the underlying instruments */
    long    end_datel,     /*	End date for diagonal */
    double* strikel_,      /*	Strikes */
    double* strikeS1l_,
    double* strikeS2l_,

    CPD_DIAG_CALIB_PARAM paraml,

    /* Short Instruments */
    char* instr_short_freq, /*	Frequency and basis of instruments */
    char* instr_short_basis,
    char* short_ref_rate_name, /*	Name of the reference rate */

    int     num_ex_datess, /*	Short Exercise dates */
    long*   ex_dates_,     /*	Supposed to be sorted */
    int*    cal_dates,     /*	1: use ex_date as calibration date, 0: don't */
    char**  end_tenors,    /*	Tenors of the underlying instruments */
    long    end_dates,     /*	End date for diagonal */
    double* strikes_,      /*	Strikes */
    double* strikeS1s_,    /*	Strike1 for smile calib */
    double* strikeS2s_,    /*	Strike2 for smile calib */
    double* weights_,      /*	Weights on secondary instruments */

    CPD_DIAG_CALIB_PARAM params,

    /*	Model */
    LGMSV_CALIBPARAMS calib_params,

    int     iOne2F,
    double* dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,
    double  dLGMAlpha,
    double  dLGMGamma,
    double  dLGMRho,
    double* dRho2TS,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /*	Output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    double** alpha,
    double** lambdaeps,
    double** rho,
    double** rho2,

    /*	Parameters */
    DIAG_CALIB_LM_PARAMS lm_params,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data) /*	NULL = don't save calibration instrument data */
{
    long today;

    double temp, fact;

    SrtCurvePtr yc_ptr;

    int i, j, last_index;

    Err err = NULL;

    long ex_datel[MAX_INST], start_datel[MAX_INST], theo_end_datel[MAX_INST],
        act_end_datel[MAX_INST], ex_dates[MAX_INST], start_dates[MAX_INST],
        theo_end_dates[MAX_INST], act_end_dates[MAX_INST];

    double strikel[MAX_INST], strikeS1l[MAX_INST], strikeS2l[MAX_INST], strikeS1s[MAX_INST],
        strikeS2s[MAX_INST], ex_timel[MAX_INST], strikes[MAX_INST], weightl[MAX_INST],
        weights[MAX_INST], ex_times[MAX_INST];

    double ex_zeta[MAX_INST];

    CalibCpnScheduleDLM *CalibCpnLongSchedule = NULL, *CalibCpnShortSchedule = NULL;

    CalibExeScheduleDLM *CalibExeLongSchedule = NULL, *CalibExeShortSchedule = NULL;

    AllCalibInstrumentsDLM *AllCalibLongInst = NULL, *AllCalibShortInst = NULL,
                           *AllCalibShortInstOne = NULL;

    void **AllLongInstruments = NULL, **AllShortInstruments = NULL;

    LGMSV_AllParams* AllParams           = NULL;
    CALIBGEN_Params* CalibParamsForSmile = NULL;

    LGMSV_model* model = NULL;

    LGMSV_CalibParams*  local_calib_params = NULL;
    CalibInstrumentDLM* NextInst;

    int     iNbPWTimeNew;
    double *dPWTimeNew = NULL, *dSigmaTS = NULL, *dAlphaTSNew = NULL, *dRhoTSNew = NULL,
           *dRho2TSNew = NULL, *dLambdaEpsTSNew = NULL, *dLvlEpsTSNew = NULL;

    long * lSigLongIndex = NULL, *lSigShortIndex = NULL;
    double time0 = 100;

    int     nb_long_smile, nb_short_smile;
    double  flat_alpha, flat_lambda, flat_rho;
    CTS_MKT mkt = NULL;
    double  null_strikes_in_vol;
    long    total_und_long, total_und_short;
    double  lam_sens;
    int     shift_lam;
    clock_t time1, time2;

    time1 = clock();

    /*	Copy data so as not to change the original */
    *sig_time  = NULL;
    *sig       = NULL;
    *alpha     = NULL;
    *lambdaeps = NULL;
    *rho       = NULL;
    *rho2      = NULL;

    /* Initialisations */
    CalibCpnLongSchedule  = calloc(1, sizeof(CalibCpnScheduleDLM));
    CalibCpnShortSchedule = calloc(1, sizeof(CalibCpnScheduleDLM));
    CalibExeLongSchedule  = calloc(1, sizeof(CalibExeScheduleDLM));
    CalibExeShortSchedule = calloc(1, sizeof(CalibExeScheduleDLM));
    AllCalibLongInst      = calloc(1, sizeof(AllCalibInstrumentsDLM));
    AllCalibShortInst     = calloc(1, sizeof(AllCalibInstrumentsDLM));
    AllCalibShortInstOne  = calloc(1, sizeof(AllCalibInstrumentsDLM));
    AllParams             = calloc(1, sizeof(LGMSV_AllParams));
    model                 = calloc(1, sizeof(LGMSV_model));
    local_calib_params    = calloc(1, sizeof(LGMSV_CalibParams));

    if (!CalibCpnLongSchedule || !CalibCpnShortSchedule || !CalibExeLongSchedule ||
        !CalibExeShortSchedule || !AllCalibLongInst || !AllCalibShortInst || !AllParams || !model ||
        !AllCalibShortInstOne)
    {
        err = "Memory allocation faillure in cpd_calib_diagonal_LGMSV_dlm";
        goto FREE_RETURN;
    }

    /* Copy the calib params */
    LGMSV_Copy_CalibParams(calib_params, local_calib_params);

    /* basic checks */
    if (iOne2F == 1)
    {
        local_calib_params->calib_rho2 = 0;
    }

    if (local_calib_params->calib_rho2 && local_calib_params->use_sabr_calib == 0)
    {
        local_calib_params->onefac_rho = 0;
    }

    if (local_calib_params->fix_lambda && !local_calib_params->calib_smile_on_prim)
    {
        err = "Cannot calibrate Smile on secondary if lambda not calibrated";
        goto FREE_RETURN;
    }

    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }

    today = get_today_from_curve(yc_ptr);

    /* Long instruments */

    /* allocation */
    memset(start_datel, 0, MAX_INST * sizeof(long));
    memset(start_dates, 0, MAX_INST * sizeof(long));

    AllocateCalibExeSchedule(
        &(ex_datel[0]),
        &(ex_timel[0]),
        &(start_datel[0]),
        &(theo_end_datel[0]),
        &(act_end_datel[0]),
        &(strikel[0]),
        &(strikeS1l[0]),
        &(strikeS2l[0]),
        &(weightl[0]),
        CalibExeLongSchedule);

    err = Construct_CalibSchedule(
        yc_name,
        today,
        instr_long_freq,
        instr_long_basis,
        num_ex_datesl,
        ex_datel_,
        cal_datel,
        end_tenorl,
        end_datel,
        strikel_,
        strikeS1l_,
        strikeS2l_,
        NULL,
        CalibCpnLongSchedule,
        CalibExeLongSchedule);

    if (err)
        goto FREE_RETURN;

    err = Reduct_ExeSchedule(CalibCpnLongSchedule, CalibExeLongSchedule, cal_datel, paraml);

    if (err)
        goto FREE_RETURN;

    if (local_calib_params->calib_alpha || local_calib_params->calib_rho)
    {
        if (local_calib_params->calib_rr_bt)
        {
            nb_long_smile = 3;
        }
        else
        {
            nb_long_smile = 1 + local_calib_params->calib_alpha + local_calib_params->calib_rho;
        }
    }
    else
    {
        nb_long_smile = 1;
    }

    err = AllocateAllCalibInst(CalibExeLongSchedule->iNExe, nb_long_smile, AllCalibLongInst);

    if (err)
        goto FREE_RETURN;

    AllLongInstruments = (void**)calloc(CalibExeLongSchedule->iNExe, sizeof(void*));

    if (!AllLongInstruments)
    {
        err = "Memory allocation faillure in cpd_calib_diagonal_LGMSV_new_dlm";
        goto FREE_RETURN;
    }

    for (i = 0; i < CalibExeLongSchedule->iNExe; i++)
    {
        AllLongInstruments[i] = (void*)(&AllCalibLongInst->sCalibInst[i]);
    }

    err = Calculate_CalibInst(
        today,
        yc_name,
        vol_curve_name,
        get_cash_vol,
        vol_ref_rate_name,
        instr_long_freq,
        instr_long_basis,
        long_ref_rate_name,
        1,
        paraml,
        CalibCpnLongSchedule,
        CalibExeLongSchedule,
        AllCalibLongInst,
        AllCalibLongInst);

    if (err)
        goto FREE_RETURN;

    /* Short instruments */
    if (!local_calib_params->fix_lambda)
    {
        /* allocation */
        AllocateCalibExeSchedule(
            &(ex_dates[0]),
            &(ex_times[0]),
            &(start_dates[0]),
            &(theo_end_dates[0]),
            &(act_end_dates[0]),
            &(strikes[0]),
            &(strikeS1s[0]),
            &(strikeS2s[0]),
            &(weights[0]),
            CalibExeShortSchedule);

        err = Construct_CalibSchedule(
            yc_name,
            today,
            instr_short_freq,
            instr_short_basis,
            num_ex_datess,
            ex_dates_,
            cal_dates,
            end_tenors,
            end_dates,
            strikes_,
            strikeS1s_,
            strikeS2s_,
            weights_,
            CalibCpnShortSchedule,
            CalibExeShortSchedule);

        if (err)
        {
            goto FREE_RETURN;
        }

        err = Reduct_ExeSchedule(CalibCpnShortSchedule, CalibExeShortSchedule, cal_dates, params);

        if (err)
            goto FREE_RETURN;

        if (!local_calib_params->calib_smile_on_prim)
        {
            if (local_calib_params->calib_alpha || local_calib_params->calib_rho)
            {
                if (local_calib_params->calib_rr_bt)
                {
                    nb_short_smile = 3;
                }
                else
                {
                    nb_short_smile =
                        1 + local_calib_params->calib_alpha + local_calib_params->calib_rho;
                }
            }
            else
            {
                nb_short_smile = 1;
            }
        }
        else
        {
            nb_short_smile = 1 + local_calib_params->calib_rho2;

            if (strikeS2s_ && strikeS2s_[0] > -100)
            {
                nb_short_smile += 1;
            }
        }

        err = AllocateAllCalibInst(CalibExeShortSchedule->iNExe, nb_short_smile, AllCalibShortInst);

        if (err)
            goto FREE_RETURN;

        AllShortInstruments = (void**)calloc(CalibExeShortSchedule->iNExe, sizeof(void*));

        if (!AllShortInstruments)
        {
            err = "Memory allocation faillure in cpd_calib_diagonal_LGMSV_new_dlm";
            goto FREE_RETURN;
        }

        for (i = 0; i < CalibExeShortSchedule->iNExe; i++)
        {
            AllShortInstruments[i] = (void*)(&AllCalibShortInst->sCalibInst[i]);
        }

        err = Calculate_CalibInst(
            today,
            yc_name,
            vol_curve_name,
            get_cash_vol,
            vol_ref_rate_name,
            instr_short_freq,
            instr_short_basis,
            short_ref_rate_name,
            1,
            params,
            CalibCpnShortSchedule,
            CalibExeShortSchedule,
            AllCalibShortInst,
            AllCalibLongInst);

        if (err)
            goto FREE_RETURN;
    }

    if ((local_calib_params->use_sabr_calib &&
         (local_calib_params->calib_alpha || local_calib_params->calib_rho ||
          local_calib_params->calib_lameps)) ||
        local_calib_params->calib_lameps)
    {
        /*	Initialise the market */
        mkt = calloc(1, sizeof(cts_mkt));

        if (!mkt)
        {
            err = "Memory allocation faillure in lgmsv_calib_dlm_new";
            goto FREE_RETURN;
        }

        null_strikes_in_vol = 1.0;

        err = cts_init_mkt(
            today,
            yc_name,
            vol_curve_name,
            vol_ref_rate_name,
            instr_long_freq,
            instr_long_basis,
            instr_short_freq,
            instr_short_basis,
            get_cash_vol,
            GetCashVolAndConvert,
            0,
            0,
            1,
            &null_strikes_in_vol,
            mkt);

        if (err)
            goto FREE_RETURN;

        if (local_calib_params->calib_smile_on_prim)
        {
            err = LGMSV_calib_const_alpha_and_rho_from_SABR(
                AllCalibLongInst,
                mkt,
                &flat_alpha,
                &flat_lambda,
                &flat_rho,
                local_calib_params,
                inst_data);
        }
        else
        {
            err = LGMSV_calib_const_alpha_and_rho_from_SABR(
                AllCalibShortInst,
                mkt,
                &flat_alpha,
                &flat_lambda,
                &flat_rho,
                local_calib_params,
                inst_data);
        }

        if (err)
            goto FREE_RETURN;

        if (local_calib_params->use_sabr_calib)
        {
            /* Override all the smile parameters */
            for (i = 0; i < iNbPWTime; i++)
            {
                if (local_calib_params->calib_alpha)
                {
                    dAlphaTS[i] = flat_alpha + local_calib_params->alpha_sv_shift;
                }

                if (local_calib_params->calib_lameps)
                {
                    dLambdaEpsTS[i] = flat_lambda + local_calib_params->lam_sv_shift;
                    ;
                }

                if (local_calib_params->calib_rho)
                {
                    dRhoTS[i] = flat_rho + local_calib_params->rho_sv_shift;
                    ;
                }
            }

            local_calib_params->calib_alpha = 0;
            local_calib_params->calib_rho   = 0;
            local_calib_params->calib_rho2  = 0;

            local_calib_params->alpha_sv_shift = 0.0;
            local_calib_params->lam_sv_shift   = 0.0;
            local_calib_params->rho_sv_shift   = 0.0;
            local_calib_params->rho2_sv_shift  = 0.0;
        }
        else
        {
            /* Just Update the Lambda */
            for (i = 0; i < iNbPWTime; i++)
            {
                dLambdaEpsTS[i] = flat_lambda;
            }
        }
    }

    /* Rescaling of rho */
    if (iOne2F == 2)
    {
        if (local_calib_params->onefac_rho)
        {
            model->dFlatOneFactorRho = dRhoTS[0];

            for (i = 0; i < iNbPWTime; i++)
            {
                err = LGMSV_Get_Rho_From_RhoTarget(
                    dRhoTS[i],
                    local_calib_params->rho_mat1,
                    local_calib_params->rho_mat2,
                    dLGMAlpha,
                    dLGMGamma,
                    dLGMRho,
                    &(dRhoTS[i]),
                    &(dRho2TS[i]));

                if (err)
                    goto FREE_RETURN;
            }
        }
        else
        {
            model->dFlatOneFactorRho = 0.0;

            for (i = 0; i < iNbPWTime; i++)
            {
                err = LGMSV_Get_RhoTarget_From_Rho(
                    dRhoTS[i],
                    dRho2TS[i],
                    local_calib_params->rho_mat1,
                    local_calib_params->rho_mat2,
                    dLGMAlpha,
                    dLGMGamma,
                    dLGMRho,
                    &temp);

                if (err)
                    goto FREE_RETURN;

                model->dFlatOneFactorRho += temp;
            }

            model->dFlatOneFactorRho /= iNbPWTime;
        }
    }

    /*	2.)	Merge the term structures */
    dPWTimeNew   = (double*)calloc(CalibExeLongSchedule->iNExe, sizeof(double));
    iNbPWTimeNew = CalibExeLongSchedule->iNExe;

    /* Init with Primary Exe Times */
    memcpy(
        dPWTimeNew, CalibExeLongSchedule->dExeTimes, CalibExeLongSchedule->iNExe * sizeof(double));

    if (iNbPWTime > 1)
    {
        /* Add the TS of Smile Times */
        num_f_concat_vector(&iNbPWTimeNew, &dPWTimeNew, iNbPWTime, dPWTime);
    }

    if (!local_calib_params->fix_lambda)
    {
        /* Add the Secondary Exe Times */
        num_f_concat_vector(
            &iNbPWTimeNew,
            &dPWTimeNew,
            CalibExeShortSchedule->iNExe,
            CalibExeShortSchedule->dExeTimes);
    }

    num_f_sort_vector(iNbPWTimeNew, dPWTimeNew);
    num_f_unique_vector(&iNbPWTimeNew, dPWTimeNew);

    dSigmaTS        = calloc(iNbPWTimeNew, sizeof(double));
    dAlphaTSNew     = calloc(iNbPWTimeNew, sizeof(double));
    dLambdaEpsTSNew = calloc(iNbPWTimeNew, sizeof(double));
    dLvlEpsTSNew    = calloc(iNbPWTimeNew, sizeof(double));
    dRhoTSNew       = calloc(iNbPWTimeNew, sizeof(double));
    dRho2TSNew      = calloc(iNbPWTimeNew, sizeof(double));

    lSigLongIndex = lvector(0, CalibExeLongSchedule->iNExe - 1);

    if (!local_calib_params->fix_lambda)
    {
        lSigShortIndex = lvector(0, CalibExeShortSchedule->iNExe - 1);
    }

    if (!dSigmaTS || !dAlphaTSNew || !dLambdaEpsTSNew || !dLvlEpsTSNew || !dRhoTSNew ||
        !dRho2TSNew || !lSigLongIndex || (!lSigShortIndex && !local_calib_params->fix_lambda))
    {
        err = "Memory allocation faillure (2) in cpd_calibSV_approx";
        goto FREE_RETURN;
    }

    for (i = 0; i < iNbPWTimeNew; i++)
    {
        j                  = Get_Index(dPWTimeNew[i], dPWTime, iNbPWTime);
        dAlphaTSNew[i]     = 2.0 * dAlphaTS[j];
        dLambdaEpsTSNew[i] = 2.0 * dLambdaEpsTS[j];
        dLvlEpsTSNew[i]    = dLambdaEpsTSNew[i];
        dRhoTSNew[i]       = dRhoTS[j];
        if (iOne2F == 2)
            dRho2TSNew[i] = dRho2TS[j];
    }

    for (i = 0; i < CalibExeLongSchedule->iNExe; i++)
    {
        lSigLongIndex[i] = Get_Index(CalibExeLongSchedule->dExeTimes[i], dPWTimeNew, iNbPWTimeNew);
    }

    if (!local_calib_params->fix_lambda)
    {
        for (i = 0; i < CalibExeShortSchedule->iNExe; i++)
        {
            lSigShortIndex[i] =
                Get_Index(CalibExeShortSchedule->dExeTimes[i], dPWTimeNew, iNbPWTimeNew);
        }
    }

    err = init_LGMSV_model(
        model,
        today,
        iOne2F,
        iNbPWTimeNew,
        *dLambdaX,
        dPWTimeNew,
        dSigmaTS,
        dAlphaTSNew,
        dLambdaEpsTSNew,
        dLvlEpsTSNew,
        dRhoTSNew,
        dTStar,
        dLGMAlpha,
        dLGMGamma,
        dLGMRho,
        dRho2TSNew);

    if (err)
        goto FREE_RETURN;

    Fill_Equivalent_ATM_Prices(model, AllCalibLongInst, lSigLongIndex, 0, CalibExeLongSchedule);

    if (!local_calib_params->fix_lambda)
    {
        Fill_Equivalent_ATM_Prices(
            model, AllCalibShortInst, lSigShortIndex, 0, CalibExeShortSchedule);
    }

    /* Check Lambda Shift */
    if (params && fabs(params->lambda_shift) > 1.0E-08)
    {
        shift_lam = 1;

        if (local_calib_params->fix_lambda)
        {
            /* We shift before calibration */
            model->dLambdaX += params->lambda_shift;
            shift_lam = 0;
        }
    }
    else
    {
        shift_lam = 0;
    }

    /* 3) Compute the first guess */
    cpd_init_hermite_for_calib(model->iOne2F);

    err = lgmcalibzetalambda_tauts_dlm(

        CalibCpnLongSchedule->iNCpn,
        CalibCpnLongSchedule->dCpnTime,
        CalibCpnLongSchedule->dCpnDf,
        CalibCpnLongSchedule->dCpnCvg,

        CalibExeLongSchedule->iNExe,
        CalibExeLongSchedule->dExeTimes,
        CalibExeLongSchedule->iStartCpn,
        CalibExeLongSchedule->iEndCpn,
        CalibExeLongSchedule->dATMStrike,
        CalibExeLongSchedule->dATMPrice,
        CalibExeLongSchedule->dATMVega,

        CalibCpnShortSchedule->iNCpn,
        CalibCpnShortSchedule->dCpnTime,
        CalibCpnShortSchedule->dCpnDf,
        CalibCpnShortSchedule->dCpnCvg,

        CalibExeShortSchedule->iNExe,
        CalibExeShortSchedule->dExeTimes,
        CalibExeShortSchedule->iStartCpn,
        CalibExeShortSchedule->iEndCpn,
        CalibExeShortSchedule->dATMStrike,
        CalibExeShortSchedule->dWeights,
        CalibExeShortSchedule->dATMPrice,
        CalibExeShortSchedule->dATMPrice,

        ex_zeta,
        local_calib_params->fix_lambda,
        1,
        &time0,
        &model->dLambdaX,
        iOne2F,
        dLGMAlpha,
        dLGMGamma,
        dLGMRho,
        paraml,
        params,
        lm_params,
        &lam_sens);

    if (err)
        goto FREE_RETURN;

    /* Check if we recalibrate the lambda */
    if (local_calib_params->use_lgm_lambda)
    {
        local_calib_params->fix_lambda = 1;

        if (shift_lam)
        {
            /* shift and recalibrate */
            model->dLambdaX += params->lambda_shift;

            err = lgmcalibzetalambda_tauts_dlm(

                CalibCpnLongSchedule->iNCpn,
                CalibCpnLongSchedule->dCpnTime,
                CalibCpnLongSchedule->dCpnDf,
                CalibCpnLongSchedule->dCpnCvg,

                CalibExeLongSchedule->iNExe,
                CalibExeLongSchedule->dExeTimes,
                CalibExeLongSchedule->iStartCpn,
                CalibExeLongSchedule->iEndCpn,
                CalibExeLongSchedule->dATMStrike,
                CalibExeLongSchedule->dATMPrice,
                CalibExeLongSchedule->dATMVega,

                CalibCpnShortSchedule->iNCpn,
                CalibCpnShortSchedule->dCpnTime,
                CalibCpnShortSchedule->dCpnDf,
                CalibCpnShortSchedule->dCpnCvg,

                CalibExeShortSchedule->iNExe,
                CalibExeShortSchedule->dExeTimes,
                CalibExeShortSchedule->iStartCpn,
                CalibExeShortSchedule->iEndCpn,
                CalibExeShortSchedule->dATMStrike,
                CalibExeShortSchedule->dWeights,
                CalibExeShortSchedule->dATMPrice,
                CalibExeShortSchedule->dATMPrice,

                ex_zeta,
                local_calib_params->fix_lambda,
                1,
                &time0,
                &model->dLambdaX,
                iOne2F,
                dLGMAlpha,
                dLGMGamma,
                dLGMRho,
                paraml,
                params,
                lm_params,
                &lam_sens);

            if (err)
                goto FREE_RETURN;
        }
    }

    /* Initialise the model with calibrated parameters */
    model->dTau = 1.0 / model->dLambdaX;

    if (model->iOne2F == 2)
    {
        model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;
        model->dTau2     = 1.0 / model->dLambdaX2;

        ConvertAlphaRho_LGM_to_LGMSV(
            model->iNbPWTime,
            model->dPWTime,
            model->dInitTStar,
            model->dLambdaX,
            model->dInitLGMAlpha,
            model->dLGMGamma,
            model->dInitLGMRho,
            model->dLGMAlpha,
            model->dLGMRho);
    }

    fact = exp(model->dLambdaX * dTStar);
    temp = sqrt(ex_zeta[0] / CalibExeLongSchedule->dExeTimes[0]) / fact;

    for (j = 0; j <= lSigLongIndex[0]; j++)
    {
        model->dSigma[j] = temp;
    }

    last_index = lSigLongIndex[0];

    for (i = 1; i < CalibExeLongSchedule->iNExe; i++)
    {
        if (ex_zeta[i] > ex_zeta[i - 1])
        {
            temp = sqrt(
                       (ex_zeta[i] - ex_zeta[i - 1]) / (CalibExeLongSchedule->dExeTimes[i] -
                                                        CalibExeLongSchedule->dExeTimes[i - 1])) /
                   fact;

            for (j = last_index + 1; j <= lSigLongIndex[i]; j++)
            {
                model->dSigma[j] = temp;
            }

            last_index = lSigLongIndex[i];
        }
        else
        {
            smessage(
                "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                CalibExeLongSchedule->dExeTimes[i]);

            for (i = i; i < CalibExeLongSchedule->iNExe; i++)
            {
                for (j = last_index + 1; j <= lSigLongIndex[i]; j++)
                {
                    model->dSigma[j] = temp;
                }

                last_index = lSigLongIndex[i];
            }
        }
    }

    /* we will price put */
    NumerParams->dIntegParam = -1.0 - NumerParams->dIntegParam;

    err = Initialise_AllParams(
        model,
        lSigLongIndex,
        CalibCpnLongSchedule,
        AllCalibLongInst,
        paraml,
        lSigShortIndex,
        CalibCpnShortSchedule,
        AllCalibShortInst,
        params,
        local_calib_params,
        NumerParams,
        model->iNbPWTime,
        AllParams);

    if (err)
        goto FREE_RETURN;

    /* Save Lam Sensi */
    AllParams->CalibParamsLambda->sensi = lam_sens;

    /* Set the instruments */
    AllParams->LongInstParams = calloc(AllCalibLongInst->iNbInst, sizeof(INSTPARAMS));

    if (!AllParams->LongInstParams)
    {
        err = "Memory allocation faillure (3) in cpd_calibSV_approx";
        goto FREE_RETURN;
    }

    for (i = 0; i < AllCalibLongInst->iNbInst; i++)
    {
        AllParams->LongInstParams[i]                    = calloc(1, sizeof(InstParams));
        AllParams->LongInstParams[i]->HestonInst        = &(AllParams->HestonLongInst[i]);
        AllParams->LongInstParams[i]->NumerInst         = &(AllParams->NumerLongInst[i]);
        AllParams->LongInstParams[i]->sCalibCpnSchedule = CalibCpnLongSchedule;
    }

    for (i = 0; i < AllCalibLongInst->iNbInst - 1; i++)
    {
        AllCalibLongInst->sCalibInst[i].NextInst     = &(AllCalibLongInst->sCalibInst[i + 1]);
        AllParams->LongInstParams[i]->NextInstParams = AllParams->LongInstParams[i + 1];
    }

    AllCalibLongInst->sCalibInst[i].NextInst     = NULL;
    AllParams->LongInstParams[i]->NextInstParams = NULL;

    if (!local_calib_params->fix_lambda)
    {
        AllParams->ShortInstParams = calloc(AllCalibShortInst->iNbInst, sizeof(INSTPARAMS));

        if (!AllParams->ShortInstParams)
        {
            err = "Memory allocation faillure (3) in cpd_calibSV_approx";
            goto FREE_RETURN;
        }

        for (i = 0; i < AllCalibShortInst->iNbInst; i++)
        {
            AllParams->ShortInstParams[i]                    = calloc(1, sizeof(InstParams));
            AllParams->ShortInstParams[i]->HestonInst        = &(AllParams->HestonShortInst[i]);
            AllParams->ShortInstParams[i]->NumerInst         = &(AllParams->NumerShortInst[i]);
            AllParams->ShortInstParams[i]->sCalibCpnSchedule = CalibCpnShortSchedule;
        }

        for (i = 0; i < AllCalibShortInst->iNbInst - 1; i++)
        {
            AllCalibShortInst->sCalibInst[i].NextInst     = &(AllCalibShortInst->sCalibInst[i + 1]);
            AllParams->ShortInstParams[i]->NextInstParams = AllParams->ShortInstParams[i + 1];
        }

        AllCalibShortInst->sCalibInst[i].NextInst     = NULL;
        AllParams->ShortInstParams[i]->NextInstParams = NULL;

        for (i = 0; i < AllCalibShortInst->iNbInst; i++)
        {
            AllCalibShortInst->sCalibInst[i].dWeight = CalibExeShortSchedule->dWeights[i];
        }

        /* Set the Lambda Sensitivity */
        total_und_long = 0;

        for (i = 0; i < AllCalibLongInst->iNbInst; i++)
        {
            total_und_long +=
                CalibCpnLongSchedule->lCpnDate[AllCalibLongInst->sCalibInst[i].iEndCpn] -
                CalibCpnLongSchedule->lCpnDate[AllCalibLongInst->sCalibInst[i].iStartCpn];
        }

        total_und_long /= AllCalibLongInst->iNbInst;

        total_und_short = 0;

        for (i = 0; i < AllCalibShortInst->iNbInst; i++)
        {
            total_und_short +=
                CalibCpnShortSchedule->lCpnDate[AllCalibShortInst->sCalibInst[i].iEndCpn] -
                CalibCpnShortSchedule->lCpnDate[AllCalibShortInst->sCalibInst[i].iStartCpn];
        }

        total_und_short /= AllCalibShortInst->iNbInst;

        if (total_und_long > total_und_short)
        {
            AllParams->sens_lambda = 1.0;
        }
        else
        {
            AllParams->sens_lambda = -1.0;
        }
    }

    AllCalibShortInstOne->iNbInst           = 1;
    AllCalibShortInstOne->sCalibCpnSchedule = AllCalibShortInst->sCalibCpnSchedule;
    AllCalibShortInstOne->sCalibExeSchedule = AllCalibShortInst->sCalibExeSchedule;
    AllCalibShortInstOne->sCalibInst        = AllCalibShortInst->sCalibInst;

    AllParams->lSmileLongIndex  = AllParams->lSigLongIndex;
    AllParams->lSmileShortIndex = AllParams->lSigShortIndex;

    AllParams->CalibCpnLongSchedule  = CalibCpnLongSchedule;
    AllParams->CalibExeLongSchedule  = CalibExeLongSchedule;
    AllParams->CalibCpnShortSchedule = CalibCpnShortSchedule;
    AllParams->CalibExeShortSchedule = CalibExeShortSchedule;

    AllParams->AllLongInst  = AllLongInstruments;
    AllParams->AllShortInst = AllShortInstruments;

    /* Set the Calibration Vol Function */
    AllParams->CalibFunctionsForVol.GetTarget              = GetTargetVol_LGMSV;
    AllParams->CalibFunctionsForVol.BumpParam              = BumpVol_LGMSV;
    AllParams->CalibFunctionsForVol.ExtrapolParam          = ExtrapolVol_LGMSV;
    AllParams->CalibFunctionsForVol.GetFirstGuess          = GetFirstGuessVol_LGMSV;
    AllParams->CalibFunctionsForVol.GetLimitAndLastParam   = GetLimitAndLastVol_LGMSV;
    AllParams->CalibFunctionsForVol.GetSecondGuess         = GetSecondGuessVol_LGMSV;
    AllParams->CalibFunctionsForVol.PriceInst              = PriceInstVol_LGMSV;
    AllParams->CalibFunctionsForVol.SetParam               = SetVol_LGMSV;
    AllParams->CalibFunctionsForVol.UpdateConstsAfterParam = UpdateParamsAfterVol_LGMSV;

    if (local_calib_params->calib_flat_smile)
    {
        AllParams->CalibFunctionsForRho.GetTarget              = GetTargetFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.BumpParam              = BumpFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.ExtrapolParam          = ExtrapolFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.GetFirstGuess          = GetFirstGuessFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.GetLimitAndLastParam   = GetLimitAndLastFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.GetSecondGuess         = GetSecondGuessFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.PriceInst              = PriceInstFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.SetParam               = SetFlatRho_LGMSV;
        AllParams->CalibFunctionsForRho.UpdateConstsAfterParam = UpdateParamsAfterFlatRho_LGMSV;

        AllParams->CalibFunctionsForAlpha.GetTarget              = GetTargetFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.BumpParam              = BumpFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.ExtrapolParam          = ExtrapolFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.GetFirstGuess          = GetFirstGuessFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.GetLimitAndLastParam   = GetLimitAndLastFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.GetSecondGuess         = GetSecondGuessFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.PriceInst              = PriceInstFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.SetParam               = SetFlatAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.UpdateConstsAfterParam = UpdateParamsAfterFlatAlpha_LGMSV;

        /* set the number of LongInstruments to 1 */
        AllCalibLongInst->iNbInst = 1;
    }
    else
    {
        AllParams->CalibFunctionsForRho.GetTarget              = GetTargetRho_LGMSV;
        AllParams->CalibFunctionsForRho.BumpParam              = BumpRho_LGMSV;
        AllParams->CalibFunctionsForRho.ExtrapolParam          = ExtrapolRho_LGMSV;
        AllParams->CalibFunctionsForRho.GetFirstGuess          = GetFirstGuessRho_LGMSV;
        AllParams->CalibFunctionsForRho.GetLimitAndLastParam   = GetLimitAndLastRho_LGMSV;
        AllParams->CalibFunctionsForRho.GetSecondGuess         = GetSecondGuessRho_LGMSV;
        AllParams->CalibFunctionsForRho.PriceInst              = PriceInstRho_LGMSV;
        AllParams->CalibFunctionsForRho.SetParam               = SetRho_LGMSV;
        AllParams->CalibFunctionsForRho.UpdateConstsAfterParam = UpdateParamsAfterRho_LGMSV;

        AllParams->CalibFunctionsForAlpha.GetTarget              = GetTargetAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.BumpParam              = BumpAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.ExtrapolParam          = ExtrapolAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.GetFirstGuess          = GetFirstGuessAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.GetLimitAndLastParam   = GetLimitAndLastAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.GetSecondGuess         = GetSecondGuessAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.PriceInst              = PriceInstAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.SetParam               = SetAlpha_LGMSV;
        AllParams->CalibFunctionsForAlpha.UpdateConstsAfterParam = UpdateParamsAfterAlpha_LGMSV;
    }

    AllParams->CalibFunctionsForLambda.GetTarget              = GetTargetLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.BumpParam              = BumpLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.ExtrapolParam          = ExtrapolLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.GetFirstGuess          = GetFirstGuessLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.GetLimitAndLastParam   = GetLimitAndLastLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.GetSecondGuess         = GetSecondGuessLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.PriceInst              = PriceInstLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.SetParam               = SetLambda_LGMSV;
    AllParams->CalibFunctionsForLambda.UpdateConstsAfterParam = UpdateParamsAfterLambda_LGMSV;

    AllParams->CalibFunctionsForRho2.GetTarget              = GetTargetRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.BumpParam              = BumpRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.ExtrapolParam          = ExtrapolRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.GetFirstGuess          = GetFirstGuessRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.GetLimitAndLastParam   = GetLimitAndLastRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.GetSecondGuess         = GetSecondGuessRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.PriceInst              = PriceInstRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.SetParam               = SetRho2_LGMSV;
    AllParams->CalibFunctionsForRho2.UpdateConstsAfterParam = UpdateParamsAfterRho2_LGMSV;

    if (local_calib_params->fix_lambda)
    {
        err = CalibrateParamTS(
            0,
            AllCalibLongInst->iNbInst - 1,
            AllLongInstruments,
            AllParams->LongInstParams,
            AllParams,
            model,
            AllParams->CalibParamsRho,
            &(AllParams->CalibFunctionsForRho));

        if (err)
            goto FREE_RETURN;
    }
    else
    {
        if (iOne2F == 1)
        {
            err = CalibrateParamTS(
                0,
                0,
                AllShortInstruments,
                AllParams->ShortInstParams,
                AllParams,
                model,
                AllParams->CalibParamsLambda,
                &(AllParams->CalibFunctionsForLambda));
        }
        else
        {
            err = CalibrateParamTS(
                0,
                0,
                AllShortInstruments,
                AllParams->ShortInstParams,
                AllParams,
                model,
                AllParams->CalibParamsRho2,
                &(AllParams->CalibFunctionsForRho2));
        }

        if (err)
            goto FREE_RETURN;
    }

    time2 = clock();
    smessage("Calibration time in sec: %.2f", (double)(time2 - time1) / CLOCKS_PER_SEC);

    smessage(
        "Total Iterations: %d, Rho2: %d, Lambda: %d, Rho: %d, Alpha: %d",
        AllParams->CalibParamsRho2->total_computed + AllParams->CalibParamsLambda->total_computed +
            AllParams->CalibParamsRho->total_computed + AllParams->CalibParamsAlpha->total_computed,
        AllParams->CalibParamsRho2->total_computed,
        AllParams->CalibParamsLambda->total_computed,
        AllParams->CalibParamsRho->total_computed,
        AllParams->CalibParamsAlpha->total_computed);

    /* if there are some shifts, then apply them and reclibrate */
    if (AllParams->CalibParams->do_calib == 0 || shift_lam ||
        fabs(local_calib_params->alpha_sv_shift) > 1.0E-08 ||
        fabs(local_calib_params->lam_sv_shift) > 1.0E-08 ||
        fabs(local_calib_params->rho_sv_shift) > 1.0E-08 ||
        (model->iOne2F == 2 && fabs(local_calib_params->rho2_sv_shift) > 1.0E-08))
    {
        /* Update Lambda */
        if (shift_lam)
        {
            model->dLambdaX += params->lambda_shift;

            if (fabs(model->dLambdaX) < 1.0E-08)
            {
                model->dLambdaX = 1.0E-08;
            }

            model->dTau = 1.0 / model->dLambdaX;

            if (model->iOne2F == 2)
            {
                model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;

                if (fabs(model->dLambdaX2) < 1.0E-08)
                {
                    model->dLambdaX2 = 1.0E-08;
                }

                model->dTau2 = 1.0 / model->dLambdaX2;

                ConvertAlphaRho_LGM_to_LGMSV(
                    model->iNbPWTime,
                    model->dPWTime,
                    model->dInitTStar,
                    model->dLambdaX,
                    model->dInitLGMAlpha,
                    model->dLGMGamma,
                    model->dInitLGMRho,
                    model->dLGMAlpha,
                    model->dLGMRho);
            }

            /* Update Long Intruments */
            for (i = 0; i < AllParams->iNbLongInst; i++)
            {
                NextInst = (CALIBINSTRUMENTDLM)(AllParams->AllLongInst[i]);

                err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
                    NextInst,
                    AllParams->CalibCpnLongSchedule,
                    model,
                    AllParams->LongInstParams[i]->HestonInst);

                if (err)
                    goto FREE_RETURN;
            }

            /* Update the Equivalent LGM */
            err = Update_EquiLGM(model, AllParams->EquiLGM);

            if (err)
                goto FREE_RETURN;
        }

        local_calib_params->fix_lambda = 1;

        /* Update Smile Term Structure */
        for (i = 0; i < model->iNbPWTime; i++)
        {
            model->dAlpha[i] += 2.0 * local_calib_params->alpha_sv_shift;
            model->dLambdaEps[i] += 2.0 * local_calib_params->lam_sv_shift;

            if (!local_calib_params->onefac_rho || model->iOne2F == 1)
            {
                model->dRho[i] += local_calib_params->rho_sv_shift;
                if (model->iOne2F == 2)
                    model->dRho2[i] += local_calib_params->rho2_sv_shift;
            }

            model->dAlpha[i]     = DMAX(model->dAlpha[i], 0.01);
            model->dLambdaEps[i] = DMAX(model->dLambdaEps[i], 0.01);
            model->dLvlEps[i]    = model->dLambdaEps[i];
            model->dRho[i]       = DMAX(model->dRho[i], -0.9999);
            model->dRho[i]       = DMIN(model->dRho[i], 0.9999);

            if (model->iOne2F == 2)
            {
                model->dRho2[i] = DMAX(model->dRho2[i], -0.9999);
                model->dRho2[i] = DMIN(model->dRho2[i], 0.9999);
            }
        }

        if (local_calib_params->onefac_rho && model->iOne2F == 2 &&
            fabs(local_calib_params->rho_sv_shift) > 1.0E-08)
        {
            for (i = 0; i < model->iNbPWTime; i++)
            {
                err = LGMSV_Get_RhoTarget_From_Rho(
                    model->dRho[i],
                    model->dRho2[i],
                    local_calib_params->rho_mat1,
                    local_calib_params->rho_mat2,
                    model->dInitLGMAlpha,
                    model->dLGMGamma,
                    model->dInitLGMRho,
                    &model->dFlatOneFactorRho);

                if (err)
                    goto FREE_RETURN;

                model->dFlatOneFactorRho += local_calib_params->rho_sv_shift;
                model->dFlatOneFactorRho = DMAX(model->dFlatOneFactorRho, -0.9999);
                model->dFlatOneFactorRho = DMIN(model->dFlatOneFactorRho, 0.9999);

                err = LGMSV_Get_Rho_From_RhoTarget(
                    model->dFlatOneFactorRho,
                    local_calib_params->rho_mat1,
                    local_calib_params->rho_mat2,
                    model->dInitLGMAlpha,
                    model->dLGMGamma,
                    model->dInitLGMRho,
                    &(model->dRho[i]),
                    &(model->dRho2[i]));

                if (err)
                    goto FREE_RETURN;
            }
        }

        err = CalibrateParamTS(
            0,
            AllParams->iNbLongInst - 1,
            AllParams->AllLongInst,
            AllParams->LongInstParams,
            AllParams,
            model,
            AllParams->CalibParams,
            &(AllParams->CalibFunctionsForVol));

        if (err)
            goto FREE_RETURN;
    }

    /* set back the number of instruments */
    if (local_calib_params->calib_flat_smile)
    {
        AllCalibLongInst->iNbInst = AllParams->iNbLongInst;
    }

    if (err)
        goto FREE_RETURN;

    NumerParams->dIntegParam = -1.0 - NumerParams->dIntegParam;

    /*	3.)	Transform into LGM sigma */
    ConvertTS_LGMSV_to_LGM(
        model->iNbPWTime, model->dPWTime, model->dSigma, model->dLambdaX, model->dTStar);

    /*	4.) Shift the Calibrated Vol */
    if (paraml->vol_type == LGM_VOL && fabs(paraml->vol_shift) > 1.0E-08)
    {
        for (i = 0; i < model->iNbPWTime; i++)
        {
            if (paraml->shift_type == MULTIPLICATIVE)
            {
                model->dSigma[i] += model->dSigma[i] * paraml->vol_shift;
            }
            else
            {
                model->dSigma[i] += paraml->vol_shift;
            }
        }
    }

    /*	5.)	Create the result */
    *dLambdaX  = model->dLambdaX;
    *num_sig   = iNbPWTimeNew;
    *sig_time  = calloc(*num_sig, sizeof(double));
    *sig       = calloc(*num_sig, sizeof(double));
    *alpha     = calloc(*num_sig, sizeof(double));
    *lambdaeps = calloc(*num_sig, sizeof(double));
    *rho       = calloc(*num_sig, sizeof(double));
    if (iOne2F == 2)
        *rho2 = calloc(*num_sig, sizeof(double));

    if (!(*sig_time) || !(*sig) || !(*alpha) || !(*lambdaeps) || !(*rho) ||
        (iOne2F == 2 && !(*rho2)))
    {
        err = "Memory allocation faillure in LGMSVCalibApprox";
        goto FREE_RETURN;
    }

    for (i = 0; i < *num_sig; i++)
    {
        (*sig_time)[i]  = model->dPWTime[i];
        (*sig)[i]       = model->dSigma[i];
        (*alpha)[i]     = model->dAlpha[i] / 2.0;
        (*lambdaeps)[i] = model->dLambdaEps[i] / 2.0;
        (*rho)[i]       = model->dRho[i];
    }

    if (iOne2F == 2)
    {
        for (i = 0; i < *num_sig; i++)
        {
            (*rho2)[i] = model->dRho2[i];
        }
    }

    /*	6.)	Save instrument data if required */
    if (inst_data)
    {
        num_ex_datesl = CalibExeLongSchedule->iNExe;

        inst_data->num_inst = num_ex_datesl;

        inst_data->exer_dates_long    = (long*)calloc(num_ex_datesl, sizeof(long));
        inst_data->start_dates        = (long*)calloc(num_ex_datesl, sizeof(long));
        inst_data->end_dates          = (long*)calloc(num_ex_datesl, sizeof(long));
        inst_data->long_strikes       = (double*)calloc(num_ex_datesl, sizeof(double));
        inst_data->market_prices_long = (double*)calloc(num_ex_datesl, sizeof(double));

        if (!calib_params->fix_lambda)
        {
            num_ex_datess = CalibExeShortSchedule->iNExe;

            inst_data->num_insts           = num_ex_datess;
            inst_data->exer_dates_short    = (long*)calloc(num_ex_datess, sizeof(long));
            inst_data->start_datess        = (long*)calloc(num_ex_datess, sizeof(long));
            inst_data->end_datess          = (long*)calloc(num_ex_datess, sizeof(long));
            inst_data->short_strikes       = (double*)calloc(num_ex_datess, sizeof(double));
            inst_data->short_weights       = (double*)calloc(num_ex_datess, sizeof(double));
            inst_data->market_prices_short = (double*)calloc(num_ex_datess, sizeof(double));

            if (!inst_data->exer_dates_long || !inst_data->start_dates || !inst_data->end_dates ||
                !inst_data->long_strikes || !inst_data->market_prices_long ||
                !inst_data->exer_dates_short || !inst_data->start_datess ||
                !inst_data->end_datess || !inst_data->short_weights || !inst_data->short_strikes ||
                !inst_data->market_prices_short)
            {
                err = "Allocation error (4) in cpd_calib_diagonal";
                goto FREE_RETURN;
            }
        }
        else
        {
            inst_data->num_insts = 0;
            if (!inst_data->exer_dates_long || !inst_data->start_dates || !inst_data->end_dates ||
                !inst_data->long_strikes || !inst_data->market_prices_long)
            {
                err = "Allocation error (4) in cpd_calib_diagonal";
                goto FREE_RETURN;
            }
        }

        for (i = 0; i < num_ex_datesl; i++)
        {
            inst_data->exer_dates_long[i] =
                (long)(CalibExeLongSchedule->dExeTimes[i] * DAYS_IN_YEAR + EPS_MAT_TO_LONG_DATE_CONVERSION) +
                today;
            inst_data->start_dates[i] =
                CalibCpnLongSchedule->lCpnDate[CalibExeLongSchedule->iStartCpn[i]];
            inst_data->end_dates[i] =
                CalibCpnLongSchedule->lCpnDate[CalibExeLongSchedule->iEndCpn[i]];
            inst_data->long_strikes[i]       = AllCalibLongInst->sCalibInst[i].dStrike[0];
            inst_data->market_prices_long[i] = CalibExeLongSchedule->dPrices[i];
        }

        if (!calib_params->fix_lambda)
        {
            for (i = 0; i < num_ex_datess; i++)
            {
                inst_data->exer_dates_short[i] =
                    (long)(CalibExeShortSchedule->dExeTimes[i] * DAYS_IN_YEAR + EPS_MAT_TO_LONG_DATE_CONVERSION) +
                    today;
                inst_data->start_datess[i] =
                    CalibCpnShortSchedule->lCpnDate[CalibExeShortSchedule->iStartCpn[i]];
                inst_data->end_datess[i] =
                    CalibCpnShortSchedule->lCpnDate[CalibExeShortSchedule->iEndCpn[i]];
                inst_data->short_strikes[i]       = AllCalibShortInst->sCalibInst[i].dStrike[0];
                inst_data->short_weights[i]       = 1.0;
                inst_data->market_prices_short[i] = CalibExeShortSchedule->dPrices[i];
            }
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

    if (dPWTimeNew)
        free(dPWTimeNew);
    if (dSigmaTS)
        free(dSigmaTS);
    if (dAlphaTSNew)
        free(dAlphaTSNew);
    if (dLambdaEpsTSNew)
        free(dLambdaEpsTSNew);
    if (dLvlEpsTSNew)
        free(dLvlEpsTSNew);
    if (dRhoTSNew)
        free(dRhoTSNew);
    if (dRho2TSNew)
        free(dRho2TSNew);
    if (lSigLongIndex)
        free_lvector(lSigLongIndex, 0, CalibExeLongSchedule->iNExe - 1);
    if (lSigShortIndex)
        free_lvector(lSigShortIndex, 0, CalibExeShortSchedule->iNExe - 1);

    if (CalibCpnLongSchedule)
    {
        free(CalibCpnLongSchedule);
    }

    if (CalibCpnShortSchedule)
    {
        free(CalibCpnShortSchedule);
    }

    if (CalibExeLongSchedule)
    {
        free(CalibExeLongSchedule);
    }

    if (CalibExeShortSchedule)
    {
        free(CalibExeShortSchedule);
    }

    if (AllCalibLongInst)
    {
        FreeAllCalibInst(AllCalibLongInst);
        free(AllCalibLongInst);
    }

    if (AllLongInstruments)
    {
        free(AllLongInstruments);
    }

    if (AllCalibShortInst)
    {
        FreeAllCalibInst(AllCalibShortInst);
        free(AllCalibShortInst);
    }

    if (AllShortInstruments)
    {
        free(AllShortInstruments);
    }

    if (AllCalibShortInstOne)
    {
        free(AllCalibShortInstOne);
    }

    if (AllParams)
    {
        Free_AllParams(AllParams);
        free(AllParams);
    }

    if (model)
    {
        free_LGMSV_model(model);
        free(model);
    }

    Free_CalibParams(CalibParamsForSmile);

    if (local_calib_params)
        free(local_calib_params);

    if (mkt)
    {
        cts_free_mkt(mkt);
        free(mkt);
    }

    return err;
}

Err Initialise_PricingConst(
    LGMSV_MODEL model, LGMSV_NUMERPARAMS NumerParams, LGMSV_PRICINGCONST PricingConst)
{
    int    i;
    double dBreakPoints[4];
    int    iBreakNbX[4];
    Err    err = NULL;

    /* allocate memory + 2 for intermediate exercise and switch time */
    PricingConst->iNbCoef = model->iNbPWTime + 2;

    PricingConst->CoefIntT = dvector(0, PricingConst->iNbCoef - 1);
    PricingConst->Coef1ReT = dvector(0, PricingConst->iNbCoef - 1);
    PricingConst->Coef2ReT = dvector(0, PricingConst->iNbCoef - 1);
    PricingConst->Coef2ImT = dvector(0, PricingConst->iNbCoef - 1);
    PricingConst->Coef3ReT = dvector(0, PricingConst->iNbCoef - 1);
    PricingConst->Coef3ImT = dvector(0, PricingConst->iNbCoef - 1);

    PricingConst->dt = dvector(0, PricingConst->iNbCoef - 1);

    PricingConst->X     = dvector(0, NumerParams->iNbX - 1);
    PricingConst->W     = dvector(0, NumerParams->iNbX - 1);
    PricingConst->IntRe = dvector(0, NumerParams->iNbX - 1);
    PricingConst->IntIm = dvector(0, NumerParams->iNbX - 1);

    if (!PricingConst->CoefIntT || !PricingConst->Coef1ReT || !PricingConst->Coef2ReT ||
        !PricingConst->Coef2ImT || !PricingConst->Coef3ReT || !PricingConst->Coef3ImT ||
        !PricingConst->dt || !PricingConst->X || !PricingConst->W || !PricingConst->IntRe ||
        !PricingConst->IntIm)
    {
        err = "Memory allocation failure in Initialise_PricingConst";
        return err;
    }

    if (NumerParams->iIntegMethod == 3 || NumerParams->iIntegMethod == 4)
    {
        PricingConst->XX = dvector(0, NumerParams->iNbX - 1);
        PricingConst->WW = dvector(0, NumerParams->iNbX - 1);

        if (!PricingConst->XX || !PricingConst->WW)
        {
            err = "Memory allocation failure (2) in Initialise_PricingConst";
            return err;
        }

        if (NumerParams->iIntegMethod == 3)
        {
            if (PricingConst->X)
                free_dvector(PricingConst->X, 0, PricingConst->iNbX - 1);
            PricingConst->X = NULL;

            if (PricingConst->W)
                free_dvector(PricingConst->W, 0, PricingConst->iNbX - 1);
            PricingConst->W = NULL;
        }
    }
    else
    {
        PricingConst->XX = NULL;
        PricingConst->WW = NULL;
    }

    PricingConst->iNbX        = NumerParams->iNbX;
    PricingConst->NumerParams = NumerParams;

    if (NumerParams->iIntegMethod == 1 || NumerParams->iIntegMethod == 4)
    {
        /* Gauss Laguerre */
        gaulag(PricingConst->X - 1, PricingConst->W - 1, NumerParams->iNbX, 0.0);
    }
    else if (NumerParams->iIntegMethod == 2)
    {
        /* Gauss Laguerre */
        gauss_hermite(PricingConst->X - 1, PricingConst->W - 1, NumerParams->iNbX);
    }

    if (NumerParams->iIntegMethod == 3 || NumerParams->iIntegMethod == 4)
    {
        /* Gauss Legendre */

        dBreakPoints[0] = LGMSV_FIRSTBREAK;
        dBreakPoints[1] = PricingConst->NumerParams->dParam1;
        dBreakPoints[2] = PricingConst->NumerParams->dParam2;
        dBreakPoints[3] = LGMSV_LASTBREAK;

        iBreakNbX[0] = (int)(NumerParams->iNbX / 4);
        iBreakNbX[1] = iBreakNbX[0];
        iBreakNbX[2] = iBreakNbX[0];
        iBreakNbX[3] = NumerParams->iNbX - 3 * iBreakNbX[0];

        LGMSVConstructLegendreGrid(4, dBreakPoints, iBreakNbX, PricingConst->XX, PricingConst->WW);
    }

    PricingConst->dt[0] = model->dPWTime[0];

    for (i = 1; i < model->iNbPWTime; i++)
    {
        PricingConst->dt[i] = model->dPWTime[i] - model->dPWTime[i - 1];
    }

    return err;
}

void Free_PricingConst(LGMSV_PRICINGCONST PricingConst)
{
    if (PricingConst)
    {
        if (PricingConst->CoefIntT)
            free_dvector(PricingConst->CoefIntT, 0, PricingConst->iNbCoef - 1);
        PricingConst->CoefIntT = NULL;

        if (PricingConst->Coef1ReT)
            free_dvector(PricingConst->Coef1ReT, 0, PricingConst->iNbCoef - 1);
        PricingConst->Coef1ReT = NULL;

        if (PricingConst->Coef2ReT)
            free_dvector(PricingConst->Coef2ReT, 0, PricingConst->iNbCoef - 1);
        PricingConst->Coef2ReT = NULL;

        if (PricingConst->Coef2ImT)
            free_dvector(PricingConst->Coef2ImT, 0, PricingConst->iNbCoef - 1);
        PricingConst->Coef2ImT = NULL;

        if (PricingConst->Coef3ReT)
            free_dvector(PricingConst->Coef3ReT, 0, PricingConst->iNbCoef - 1);
        PricingConst->Coef3ReT = NULL;

        if (PricingConst->Coef3ImT)
            free_dvector(PricingConst->Coef3ImT, 0, PricingConst->iNbCoef - 1);
        PricingConst->Coef3ImT = NULL;

        if (PricingConst->dt)
            free_dvector(PricingConst->dt, 0, PricingConst->iNbCoef - 1);
        PricingConst->dt = NULL;

        if (PricingConst->X)
            free_dvector(PricingConst->X, 0, PricingConst->iNbX - 1);
        PricingConst->X = NULL;

        if (PricingConst->W)
            free_dvector(PricingConst->W, 0, PricingConst->iNbX - 1);
        PricingConst->W = NULL;

        if (PricingConst->XX)
            free_dvector(PricingConst->XX, 0, PricingConst->iNbX - 1);
        PricingConst->XX = NULL;

        if (PricingConst->WW)
            free_dvector(PricingConst->WW, 0, PricingConst->iNbX - 1);
        PricingConst->WW = NULL;

        if (PricingConst->IntRe)
            free_dvector(PricingConst->IntRe, 0, PricingConst->iNbX - 1);
        PricingConst->IntRe = NULL;

        if (PricingConst->IntIm)
            free_dvector(PricingConst->IntIm, 0, PricingConst->iNbX - 1);
        PricingConst->IntIm = NULL;
    }
}

Err Initialise_NumerInst(
    LGMSV_MODEL        model,
    CALIBINSTRUMENTDLM CalibInst,
    LGMSV_HESTONINST   HestonParam,
    LGMSV_NUMERINST    NumerInst)
{
    Err err = NULL;

    NumerInst->iNbStrike = CalibInst->iNbStrike;

    if (CalibInst->iNbStrike > 0)
    {
        NumerInst->dLogShiftStrike = calloc(CalibInst->iNbStrike, sizeof(double));

        if (!NumerInst->dLogShiftStrike)
        {
            err = "Memory allocation faillure in Initialise_NumerInst";
            return err;
        }
    }
    else
    {
        NumerInst->dLogShiftStrike = NULL;
    }

    NumerInst->iEndIndex = Get_Index(CalibInst->dExeTime, model->dPWTime, model->iNbPWTime);

    if (CalibInst->iIsLiborOption)
    {
        NumerInst->iLiborIndex =
            Get_Index(CalibInst->dLiborFixTime, model->dPWTime, model->iNbPWTime);
        NumerInst->iSwitchIndex = NumerInst->iLiborIndex;

        if (fabs(CalibInst->dLiborFixTime - model->dPWTime[NumerInst->iSwitchIndex]) > 1.0E-08 &&
            CalibInst->dLiborFixTime < CalibInst->dExeTime)
        {
            NumerInst->iDoSwitch      = 1;
            NumerInst->dNewSwitchTime = CalibInst->dLiborFixTime;
        }
        else
        {
            NumerInst->iDoSwitch = 0;

            if (CalibInst->dLiborFixTime > CalibInst->dExeTime)
            {
                NumerInst->dNewSwitchTime = CalibInst->dExeTime;
                NumerInst->iSwitchIndex   = NumerInst->iEndIndex;
            }
        }
    }
    else
    {
        NumerInst->iDoSwitch      = 0;
        NumerInst->dNewSwitchTime = CalibInst->dExeTime;
        NumerInst->iSwitchIndex   = NumerInst->iEndIndex;
    }

    return err;
}

void Update_NumerInst(
    CALIBINSTRUMENTDLM CalibInst, LGMSV_HESTONINST HestonParam, LGMSV_NUMERINST NumerInst)
{
    int i;

    for (i = 0; i < CalibInst->iNbStrike; i++)
    {
        NumerInst->dLogShiftStrike[i] = log(fabs(CalibInst->dStrike[i] + HestonParam->dShift));
    }

    NumerInst->dLogShiftFwd =
        log(fabs(CalibInst->dFwdCash + CalibInst->dSpread + HestonParam->dShift));
}

void Free_NumerInst(LGMSV_NUMERINST NumerInst)
{
    if (NumerInst)
    {
        if (NumerInst->dLogShiftStrike)
            free(NumerInst->dLogShiftStrike);
        NumerInst->dLogShiftStrike = NULL;
    }
}

Err Initialise_EquiLGM(LGMSV_MODEL model, LGMSV_EQUILGM equi_lgm)
{
    int nb_lgm_vol;
    Err err = NULL;

    nb_lgm_vol = model->iNbPWTime;

    equi_lgm->cum_var_sv  = (double*)calloc(nb_lgm_vol, sizeof(double));
    equi_lgm->cum_var_lgm = (double*)calloc(nb_lgm_vol, sizeof(double));
    equi_lgm->exp_fact    = (double*)calloc(nb_lgm_vol + 1, sizeof(double));
    equi_lgm->lgm_vol     = (double*)calloc(nb_lgm_vol, sizeof(double));

    if (!equi_lgm->cum_var_sv || !equi_lgm->cum_var_lgm || !equi_lgm->exp_fact ||
        !equi_lgm->lgm_vol)
    {
        return "Memory allocation faillure in Initialise_EquiLGM";
    }

    err = Update_EquiLGM(model, equi_lgm);

    return err;
}

Err Update_EquiLGM(LGMSV_MODEL model, LGMSV_EQUILGM equi_lgm)
{
    int i, nb_lgm_vol;

    equi_lgm->exp_fact[0] = exp(-2.0 * model->dLambdaX * model->dInitTStar);

    nb_lgm_vol = model->iNbPWTime;
    i          = 0;

    for (i = 0; i < nb_lgm_vol; i++)
    {
        equi_lgm->exp_fact[i + 1] =
            exp(-2.0 * model->dLambdaX * (model->dInitTStar - model->dPWTime[i]));
    }

    return NULL;
}

void Free_EquiLGM(LGMSV_EQUILGM equi_lgm)
{
    if (equi_lgm)
    {
        if (equi_lgm->cum_var_sv)
            free(equi_lgm->cum_var_sv);
        equi_lgm->cum_var_sv = NULL;

        if (equi_lgm->cum_var_lgm)
            free(equi_lgm->cum_var_lgm);
        equi_lgm->cum_var_lgm = NULL;

        if (equi_lgm->exp_fact)
            free(equi_lgm->exp_fact);
        equi_lgm->exp_fact = NULL;

        if (equi_lgm->lgm_vol)
            free(equi_lgm->lgm_vol);
        equi_lgm->lgm_vol = NULL;
    }
}

Err Calculate_SV_Sensi(
    LGMSV_MODEL            model,
    ALLCALININSTRUMENTSDLM AllCalibInst,
    CPD_DIAG_CALIB_PARAM   CPDParams,
    LGMSV_CALIBPARAMS      CalibParams,
    double*                sensi_alpha,
    double*                sensi_rho)
{
    int                i, j, k;
    double             alpha, rho, beta_vol;
    double             alpha_shift, rho_shift;
    double             init_price[3], alpha_price[3], rho_price[3];
    double             price;
    CALIBINSTRUMENTDLM Inst;
    Err                err = NULL;

    /* Initialisation */
    for (i = 0; i < 3; i++)
    {
        init_price[i]  = 0.0;
        alpha_price[i] = 0.0;
        rho_price[i]   = 0.0;
    }

    alpha_shift = 0.2;
    rho_shift   = -0.1;

    for (i = 0; i < AllCalibInst->iNbInst; i++)
    {
        Inst = &(AllCalibInst->sCalibInst[i]);

        for (k = 0; k < 3; k++)
        {
            switch (k)
            {
            case 0:
            {
                /* Equi Alpha and Rho */
                alpha = 0.5 * model->dAlpha[0] *
                        sqrt(
                            (1.0 - exp(-model->dLambdaEps[0] * Inst->dExeTime)) /
                            (model->dLambdaEps[0] * Inst->dExeTime));
                rho = model->dFlatOneFactorRho;
                break;
            }
            case 1:
            {
                /* Equi Alpha and Rho */
                alpha = 0.5 * (model->dAlpha[0] + alpha_shift) *
                        sqrt(
                            (1.0 - exp(-model->dLambdaEps[0] * Inst->dExeTime)) /
                            (model->dLambdaEps[0] * Inst->dExeTime));
                break;
            }
            case 2:
            {
                /* Equi Alpha and Rho */
                alpha = 0.5 * model->dAlpha[0] *
                        sqrt(
                            (1.0 - exp(-model->dLambdaEps[0] * Inst->dExeTime)) /
                            (model->dLambdaEps[0] * Inst->dExeTime));
                rho += rho_shift;
                break;
            }
            }

            err = srt_f_optsarbvol(
                Inst->dFwdCash,
                Inst->dFwdCash,
                Inst->dExeTime,
                Inst->dATMVol,
                alpha,
                0.0,
                rho,
                SRT_NORMAL,
                SRT_BETAVOL,
                &beta_vol);

            if (err)
                return err;

            /*  Init price */
            for (j = 0; j < Inst->iNbStrike; j++)
            {
                err = srt_f_optsarbvol(
                    Inst->dFwdCash,
                    Inst->dStrike[j],
                    Inst->dExeTime,
                    beta_vol,
                    alpha,
                    0.0,
                    rho,
                    SRT_BETAVOL,
                    SRT_NORMAL,
                    &price);

                if (err)
                    return err;

                price = srt_f_optblknrm(
                    Inst->dFwdCash,
                    Inst->dStrike[j],
                    price,
                    Inst->dExeTime,
                    Inst->dLevel,
                    SRT_PUT,
                    SRT_PREMIUM);

                if (CPDParams->vega_prec)
                {
                    price /= Inst->dVega[0];
                }

                switch (k)
                {
                case 0:
                {
                    init_price[j] += price;
                    break;
                }

                case 1:
                {
                    alpha_price[j] += price;
                    break;
                }

                case 2:
                {
                    rho_price[j] += price;
                    break;
                }
                }
            }
        }
    }

    /* Final Result */
    *sensi_alpha = 0.0;
    *sensi_rho   = 0.0;

    if (CalibParams->calib_rr_bt)
    {
        if (CalibParams->calib_alpha)
        {
            *sensi_alpha = ((alpha_price[1] + alpha_price[2] - 2.0 * alpha_price[0]) -
                            (init_price[1] + init_price[2] - 2.0 * init_price[0])) /
                           alpha_shift;
        }

        if (CalibParams->calib_rho || CalibParams->calib_rho2)
        {
            *sensi_rho =
                ((rho_price[1] - rho_price[2]) - (init_price[1] - init_price[2])) / rho_shift;
        }
    }
    else
    {
        if (CalibParams->calib_alpha)
        {
            *sensi_alpha = (alpha_price[1] - init_price[1]) / 0.1;
        }

        if (CalibParams->calib_rho || CalibParams->calib_rho2)
        {
            *sensi_rho = (rho_price[1] - init_price[1]) / (-0.1);
        }
    }

    return err;
}

Err Initialise_AllParams(
    LGMSV_MODEL model,

    long*                  lSigLongIndex,
    CALIBCPNSCHEDULEDLM    CalibLongCpnSchedule,
    ALLCALININSTRUMENTSDLM AllCalibLongInst,
    CPD_DIAG_CALIB_PARAM   paraml,

    long*                  lSigShortIndex,
    CALIBCPNSCHEDULEDLM    CalibShortCpnSchedule,
    ALLCALININSTRUMENTSDLM AllCalibShortInst,
    CPD_DIAG_CALIB_PARAM   params,

    LGMSV_CalibParams* CalibParams,

    LGMSV_NUMERPARAMS NumerParams,

    int NbLGMVol,

    LGMSV_ALLPARAMS AllParams)
{
    int    i, calib_vol;
    double sensi1, sensi2;
    Err    err = NULL;

    AllParams->HestonLongInst = calloc(AllCalibLongInst->iNbInst, sizeof(LGMSV_HestonInst));
    AllParams->NumerLongInst  = calloc(AllCalibLongInst->iNbInst, sizeof(LGMSV_NumerInst));
    AllParams->lSigLongIndex  = calloc(AllCalibLongInst->iNbInst, sizeof(long));

    if (!CalibParams->fix_lambda)
    {
        AllParams->HestonShortInst = calloc(AllCalibShortInst->iNbInst, sizeof(LGMSV_HestonInst));
        AllParams->NumerShortInst  = calloc(AllCalibShortInst->iNbInst, sizeof(LGMSV_NumerInst));
        AllParams->lSigShortIndex  = calloc(AllCalibShortInst->iNbInst, sizeof(long));
    }

    AllParams->PricingConst = calloc(1, sizeof(LGMSV_PricingConst));

    AllParams->CalibParams       = calloc(1, sizeof(CALIBGEN_Params));
    AllParams->CalibParamsLambda = calloc(1, sizeof(CALIBGEN_Params));
    AllParams->CalibParamsAlpha  = calloc(1, sizeof(CALIBGEN_Params));
    AllParams->CalibParamsLamEps = calloc(1, sizeof(CALIBGEN_Params));
    AllParams->CalibParamsRho    = calloc(1, sizeof(CALIBGEN_Params));
    AllParams->CalibParamsRho2   = calloc(1, sizeof(CALIBGEN_Params));

    AllParams->EquiLGM = calloc(1, sizeof(LGMSV_EquiLGM));

    if (!AllParams->HestonLongInst || !AllParams->NumerLongInst || !AllParams->lSigLongIndex ||
        (!AllParams->HestonShortInst && !CalibParams->fix_lambda) ||
        (!AllParams->NumerShortInst && !CalibParams->fix_lambda) ||
        (!AllParams->lSigShortIndex && !CalibParams->fix_lambda) || !AllParams->PricingConst ||
        !AllParams->CalibParams || !AllParams->EquiLGM || !AllParams->CalibParamsLambda ||
        !AllParams->CalibParamsRho || !AllParams->CalibParamsAlpha || !AllParams->CalibParamsRho2)
    {
        err = "Memory allocation faillure in Initialise_AllParams";
        return err;
    }

    AllParams->calib_params = CalibParams;
    AllParams->long_param   = paraml;
    AllParams->short_param  = params;

    err = Initialise_PricingConst(model, NumerParams, AllParams->PricingConst);

    if (err)
        return err;

    AllParams->iIndexStrike = 0;
    AllParams->SkipLast     = paraml->skip_last;
    AllParams->MinFact      = paraml->min_fact;
    AllParams->MaxFact      = paraml->max_fact;

    if (!CalibParams->fix_lambda || CalibParams->calib_alpha || CalibParams->calib_lameps ||
        CalibParams->calib_rho || CalibParams->calib_rho2)
    {
        if (CalibParams->novolcalib_on_smile_calib)
        {
            calib_vol = 0;
        }
        else
        {
            calib_vol = 1;
        }
    }

    err = Initialise_CalibParams(
        paraml->use_jumps,
        paraml->precision,
        paraml->nb_iter_max,
        0,
        calib_vol,
        0,
        0,
        0.0,
        AllParams->CalibParams);

    if (err)
        return err;

    err = Initialise_CalibParams(
        params->use_jumps,
        params->precision,
        params->nb_iter_max,
        CalibParams->recalib_at_end_lambda,
        !CalibParams->fix_lambda,
        0,
        1,
        0.001,
        AllParams->CalibParamsLambda);

    if (err)
        return err;

    err = Initialise_CalibParams(
        paraml->smile_use_jumps,
        paraml->smile_precision,
        paraml->smile_nb_iter_max,
        CalibParams->recalib_at_end_alpha,
        CalibParams->calib_alpha,
        0,
        1,
        0.005,
        AllParams->CalibParamsAlpha);

    if (err)
        return err;

    err = Initialise_CalibParams(
        paraml->smile_use_jumps,
        paraml->smile_precision,
        paraml->smile_nb_iter_max,
        CalibParams->recalib_at_end_lameps,
        CalibParams->calib_lameps,
        0,
        1,
        0.001,
        AllParams->CalibParamsLamEps);

    if (err)
        return err;

    err = Initialise_CalibParams(
        paraml->smile_use_jumps,
        paraml->smile_precision,
        paraml->smile_nb_iter_max,
        CalibParams->recalib_at_end_rho,
        CalibParams->calib_rho,
        0,
        1,
        0.005,
        AllParams->CalibParamsRho);

    if (err)
        return err;

    err = Initialise_CalibParams(
        params->smile_use_jumps,
        params->smile_precision,
        params->smile_nb_iter_max,
        CalibParams->recalib_at_end_rho,
        CalibParams->calib_rho2,
        0,
        1,
        0.005,
        AllParams->CalibParamsRho2);

    if (err)
        return err;

    memcpy(AllParams->lSigLongIndex, lSigLongIndex, AllCalibLongInst->iNbInst * sizeof(long));
    AllParams->iNbLongInst = AllCalibLongInst->iNbInst;

    for (i = 0; i < AllCalibLongInst->iNbInst; i++)
    {
        Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
            &(AllCalibLongInst->sCalibInst[i]),
            CalibLongCpnSchedule,
            model,
            &(AllParams->HestonLongInst[i]));

        Calculate_ShiftAndVol_LGMSV_FromLGM_struct(
            &(AllCalibLongInst->sCalibInst[i]),
            CalibLongCpnSchedule,
            model,
            &(AllParams->HestonLongInst[i]));

        err = Initialise_NumerInst(
            model,
            &(AllCalibLongInst->sCalibInst[i]),
            &(AllParams->HestonLongInst[i]),
            &(AllParams->NumerLongInst[i]));

        if (err)
            return err;
    }

    if (!CalibParams->fix_lambda)
    {
        memcpy(
            AllParams->lSigShortIndex, lSigShortIndex, AllCalibShortInst->iNbInst * sizeof(long));
        AllParams->iNbShortInst = AllCalibShortInst->iNbInst;

        for (i = 0; i < AllCalibShortInst->iNbInst; i++)
        {
            Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
                &(AllCalibShortInst->sCalibInst[i]),
                CalibShortCpnSchedule,
                model,
                &(AllParams->HestonShortInst[i]));

            Calculate_ShiftAndVol_LGMSV_FromLGM_struct(
                &(AllCalibShortInst->sCalibInst[i]),
                CalibShortCpnSchedule,
                model,
                &(AllParams->HestonShortInst[i]));

            err = Initialise_NumerInst(
                model,
                &(AllCalibShortInst->sCalibInst[i]),
                &(AllParams->HestonShortInst[i]),
                &(AllParams->NumerShortInst[i]));

            if (err)
            {
                return err;
            }
        }
    }

    err = Initialise_EquiLGM(model, AllParams->EquiLGM);

    if (err)
        return err;

    /* Sensitivities */
    if (CalibParams->calib_alpha || CalibParams->calib_rho || CalibParams->calib_rho2)
    {
        err = Calculate_SV_Sensi(
            model,
            AllCalibLongInst,
            paraml,
            CalibParams,
            &AllParams->CalibParamsAlpha->sensi,
            &AllParams->CalibParamsRho->sensi);

        if (CalibParams->onefac_rho)
        {
            AllParams->CalibParamsRho2->sensi = AllParams->CalibParamsRho->sensi;
        }
        else
        {
            if (model->iOne2F == 2)
            {
                err = LGMSV_Get_RhoTargetSensitivity_From_Rho(
                    model->dRho[0],
                    model->dRho2[0],
                    CalibParams->rho_mat1,
                    CalibParams->rho_mat2,
                    model->dInitLGMAlpha,
                    model->dLGMGamma,
                    model->dInitLGMRho,
                    &sensi1,
                    &sensi2);

                if (err)
                    return err;
            }

            /* doesn't work */
            /*
            AllParams->CalibParamsRho2->sensi = AllParams->CalibParamsRho->sensi * sensi2;
            AllParams->CalibParamsRho->sensi *= sensi1;
            */

            AllParams->CalibParamsRho->sensi  = 0.0;
            AllParams->CalibParamsRho2->sensi = 0.0;
        }
    }

    if (err)
        return err;

    return err;
}

void Free_AllParams(LGMSV_ALLPARAMS AllParams)
{
    int i;
    Err err = NULL;

    if (AllParams)
    {
        if (AllParams->lSigLongIndex)
            free(AllParams->lSigLongIndex);
        AllParams->lSigLongIndex = NULL;

        for (i = 0; i < AllParams->iNbLongInst; i++)
        {
            Free_NumerInst(&(AllParams->NumerLongInst[i]));
        }

        if (AllParams->LongInstParams)
        {
            for (i = 0; i < AllParams->iNbLongInst; i++)
            {
                if (AllParams->LongInstParams[i])
                    free(AllParams->LongInstParams[i]);
            }

            free(AllParams->LongInstParams);
            AllParams->LongInstParams = NULL;
        }

        if (AllParams->HestonLongInst)
            free(AllParams->HestonLongInst);
        AllParams->HestonLongInst = NULL;

        if (AllParams->NumerLongInst)
            free(AllParams->NumerLongInst);
        AllParams->NumerLongInst = NULL;

        if (AllParams->CalibParams && !AllParams->calib_params->fix_lambda)
        {
            if (AllParams->lSigShortIndex)
                free(AllParams->lSigShortIndex);
            AllParams->lSigShortIndex = NULL;

            for (i = 0; i < AllParams->iNbShortInst; i++)
            {
                Free_NumerInst(&(AllParams->NumerShortInst[i]));
            }

            if (AllParams->ShortInstParams)
            {
                for (i = 0; i < AllParams->iNbShortInst; i++)
                {
                    if (AllParams->ShortInstParams[i])
                        free(AllParams->ShortInstParams[i]);
                }

                free(AllParams->ShortInstParams);
                AllParams->ShortInstParams = NULL;
            }

            if (AllParams->HestonShortInst)
                free(AllParams->HestonShortInst);
            AllParams->HestonShortInst = NULL;

            if (AllParams->NumerShortInst)
                free(AllParams->NumerShortInst);
            AllParams->NumerShortInst = NULL;
        }

        Free_PricingConst(AllParams->PricingConst);
        if (AllParams->PricingConst)
            free(AllParams->PricingConst);
        AllParams->PricingConst = NULL;

        if (AllParams->CalibParams)
        {
            Free_CalibParams(AllParams->CalibParams);
            free(AllParams->CalibParams);
        }

        if (AllParams->CalibParamsLambda)
        {
            Free_CalibParams(AllParams->CalibParamsLambda);
            free(AllParams->CalibParamsLambda);
        }

        if (AllParams->CalibParamsAlpha)
        {
            Free_CalibParams(AllParams->CalibParamsAlpha);
            free(AllParams->CalibParamsAlpha);
        }

        if (AllParams->CalibParamsRho)
        {
            Free_CalibParams(AllParams->CalibParamsRho);
            free(AllParams->CalibParamsRho);
        }

        if (AllParams->CalibParamsLamEps)
        {
            Free_CalibParams(AllParams->CalibParamsLamEps);
            free(AllParams->CalibParamsLamEps);
        }

        if (AllParams->CalibParamsRho2)
        {
            Free_CalibParams(AllParams->CalibParamsRho2);
            free(AllParams->CalibParamsRho2);
        }

        if (AllParams->EquiLGM)
        {
            Free_EquiLGM(AllParams->EquiLGM);
            free(AllParams->EquiLGM);
        }
    }
}

Err LGMSVBumpLambdaModel(void* model_, int iIndex, double dNewValue)
{
    LGMSV_MODEL model = (LGMSV_MODEL)(model_);

    model->dLambdaX = dNewValue;

    if (fabs(dNewValue) > 1.0E-08)
    {
        model->dTau = 1.0 / dNewValue;
    }
    else
    {
        return "LambdaX cannot be NULL";
    }

    if (model->iOne2F == 2)
    {
        model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;

        if (fabs(model->dLambdaX2) > 1.0E-08)
        {
            model->dTau2 = 1.0 / model->dLambdaX2;
        }
        else
        {
            return "Invalid Gamma in init_LGMSV_model";
        }

        /* adjust correlations */
        ConvertAlphaRho_LGM_to_LGMSV(
            model->iNbPWTime,
            model->dPWTime,
            model->dTStar,
            model->dLambdaX,
            model->dInitLGMAlpha,
            model->dLGMGamma,
            model->dInitLGMRho,
            model->dLGMAlpha,
            model->dLGMRho);
    }

    return NULL;
}

Err LGMSVPricingShortInstruments(/* Instrument informations	*/
                                 void* AllCalibInst_,

                                 /* Model informations */
                                 void* model_,

                                 /* Parameter of grids */
                                 void* AllParams_,

                                 int     iIndex,
                                 double* dPrice)
{
    ALLCALININSTRUMENTSDLM AllCalibInst = (ALLCALININSTRUMENTSDLM)(AllCalibInst_);
    LGMSV_MODEL            model        = (LGMSV_MODEL)(model_);
    LGMSV_ALLPARAMS        AllParams    = (LGMSV_ALLPARAMS)(AllParams_);

    CALIBINSTRUMENTDLM CalibInst;

    Err err = NULL;

    CalibInst = &(AllCalibInst->sCalibInst[iIndex]);

    Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
        &(AllCalibInst->sCalibInst[iIndex]),
        AllCalibInst->sCalibCpnSchedule,
        model,
        &(AllParams->HestonShortInst[iIndex]));

    if (err)
    {
        return err;
    }

    err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(
        CalibInst, AllCalibInst->sCalibCpnSchedule, model, &(AllParams->HestonShortInst[iIndex]));

    if (err)
    {
        return err;
    }

    LGMSVClosedFormApprox_struct(
        model,
        CalibInst,
        &(AllParams->HestonShortInst[iIndex]),
        AllParams->PricingConst,
        &(AllParams->NumerShortInst[iIndex]),
        dPrice);

    return err;
}
