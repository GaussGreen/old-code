
#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMMC.h"
#include "Fx3FBetaDLMTree.h"
#include "Fx3FBetaDLMUtil.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"

#define CPD_MINCALDAYS 25

/* DLM parameters */
void cpd_init_beta_dlm_params(cpd_beta_dlm_params* cpd_dlm_params)
{
    cpd_dlm_params->use_beta_dlm         = 0;
    cpd_dlm_params->use_beta_dlm_for_iv  = 0;
    cpd_dlm_params->use_tree_iv_for_call = 0;
    cpd_dlm_params->calib_smile          = 0;
    cpd_dlm_params->calib_smile_std_up   = 0.5;
    cpd_dlm_params->calib_smile_std_down = 0.5;
    cpd_dlm_params->calib_smile_maturity = 3.0;
    cpd_dlm_params->calib_smile_shift    = 0.0;

    cpd_dlm_params->tstar  = 30.0;
    cpd_dlm_params->B0     = 1.0;
    cpd_dlm_params->C0     = 0.0;
    cpd_dlm_params->alpha  = 0.0;
    cpd_dlm_params->lambda = 1.0;

    FxBetaDLM_Init_GRFNNumerParams(cpd_dlm_params->GrfnParams);
    FxBetaDLM_SetDefaultOptNumerParams(cpd_dlm_params->NumParams);
}

Err cpd_free_beta_dlm_params(cpd_beta_dlm_params* cpd_dlm_params)
{
    Err err = NULL;

    if (cpd_dlm_params)
    {
        if (cpd_dlm_params->NumParams)
        {
            free(cpd_dlm_params->NumParams);
            cpd_dlm_params->NumParams = NULL;
        }

        if (cpd_dlm_params->GrfnParams)
        {
            free(cpd_dlm_params->GrfnParams);
            cpd_dlm_params->GrfnParams = NULL;
        }
    }

    return err;
}

/*	Functions for the funding leg */

Err cpd_fill_fund_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int         eod_flag, /*	0: I, 1: E */
    double      fund_not,
    int         fund_ccy, /*	0: domestic, 1: foreign */
    int         fund_ncpn,
    long*       fund_fix,
    long*       fund_start,
    long*       fund_pay,
    char**      fund_basis,
    double*     fund_spr,
    double*     fund_mrg,
    PD_FUND_LEG fund_leg)
{
    int          i, j;
    SrtBasisCode bas;
    double       cvg;
    Err          err = NULL;

    /* Initialise notional */
    fund_leg->dom_for  = fund_ccy;
    fund_leg->notional = fund_not;

    /*	Initialise pointers to NULL */
    fund_leg->cpn = NULL;

    /*	Skip coupons fixed before today */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_flag)
        i++;

    /*	Check that at least one coupon is left */
    if (i == fund_ncpn)
    {
        /* err = "All funding coupons start before today in cpd_fill_fund_leg"; */
        fund_leg->num_cpn = 0;
        goto FREE_RETURN;
    }

    /*	Allocate memory */
    fund_leg->num_cpn = fund_ncpn - i;
    fund_leg->cpn     = (pd_fund_cpn*)calloc(fund_leg->num_cpn, sizeof(pd_fund_cpn));
    if (!fund_leg->cpn)
    {
        err = "Allocation error in cpd_fill_fund_leg";
        goto FREE_RETURN;
    }

    /*	Fill coupons information */
    j = 0;
    while (i < fund_ncpn)
    {
        /*	Dates */
        fund_leg->cpn[j].start_date = fund_start[i];
        fund_leg->cpn[j].pay_date   = fund_pay[i];

        /*	Times */
        fund_leg->cpn[j].start_time = (fund_leg->cpn[j].start_date - today) * YEARS_IN_DAY;
        fund_leg->cpn[j].pay_time   = (fund_leg->cpn[j].pay_date - today) * YEARS_IN_DAY;

        /*	Coupon */
        err = interp_basis(fund_basis[i], &bas);
        if (err)
            goto FREE_RETURN;

        cvg                  = coverage(fund_start[i], fund_pay[i], bas);
        fund_leg->cpn[j].cpn = fund_not * cvg * (fund_spr[i] + fund_mrg[i]);

        i++;
        j++;
    }

    err = cpd_check_fund_leg(fund_leg);

FREE_RETURN:

    if (err)
        cpd_free_fund_leg(fund_leg);
    return err;
}

/*	Check dates consistency */
Err cpd_check_fund_leg(PD_FUND_LEG fund_leg)
{
    int i;

    /*	Check that start and pay dates are increasing */
    for (i = 1; i < fund_leg->num_cpn; i++)
    {
        if (fund_leg->cpn[i].start_date < fund_leg->cpn[i - 1].start_date)
            return "Start dates should be increasing in funding leg";

        if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i - 1].pay_date)
            return "Pay dates should be increasing in funding leg";
    }

    /*	Check that pay dates are after start dates */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i].start_date)
            return "Pay dates should be after start dates in funding leg";
    }

    /*	OK */
    return NULL;
}

/*	Free */
Err cpd_free_fund_leg(PD_FUND_LEG fund_leg)
{
    if (fund_leg->cpn)
    {
        free(fund_leg->cpn);
        fund_leg->cpn = NULL;
    }

    return NULL;
}

/*	Functions for the exotic leg */

Err cpd_fill_exo_leg(
    /*	Coupons that fixed before today are disregarded */
    long today,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    double  pd_not,
    int     pd_ncpn,
    long*   pd_fix,
    long*   pd_start,
    long*   pd_pay,
    char**  pd_basis,
    double* pd_alpha,
    double* pd_beta,
    int*    pd_floored,
    double* pd_floor,
    int*    pd_capped,
    double* pd_cap,

    //		pd interp coupon:
    int      use_cpn_opt_str,
    int*     pd_num_strikes,
    double*  pd_wcst,
    double*  pd_wspot,
    double** pd_strikes,
    double** pd_weights,

    /*		pd not refund */
    long   pd_not_ref_fix,
    double pd_not_ref_alpha,
    double pd_not_ref_beta,
    int    pd_not_ref_floored,
    double pd_not_ref_floor,
    int    pd_not_ref_capped,
    double pd_not_ref_cap,

    //		pd interp notional:
    int     use_not_opt_str,
    int     pd_not_num_strikes,
    double  pd_not_wcst,
    double  pd_not_wspot,
    double* pd_not_strikes,
    double* pd_not_weights,

    PD_EXO_LEG exo_leg)
{
    int          i, j, k;
    SrtBasisCode bas;
    double       cvg;
    Err          err = NULL;

    /*	Initialise pointers to NULL */
    exo_leg->cpn = NULL;

    /*	Skip coupons fixed before today */
    i = 0;
    while (i < pd_ncpn && pd_fix[i] < today + eod_flag)
        i++;

    /*	Check that at least one coupon is left */
    if (i == pd_ncpn)
    {
        /* err = "All funding coupons start before today in cpd_fill_exo_leg"; */
        exo_leg->num_cpn = 0;
        goto DO_NOTIONAL_REF;
    }

    /*	Allocate memory */
    exo_leg->num_cpn = pd_ncpn - i;
    exo_leg->cpn     = (pd_exo_cpn*)calloc(exo_leg->num_cpn, sizeof(pd_exo_cpn));
    if (!exo_leg->cpn)
    {
        err = "Allocation error in cpd_fill_exo_leg";
        goto FREE_RETURN;
    }

    /*	Fill coupons information */
    j = 0;
    while (i < pd_ncpn)
    {
        /*	Dates */
        exo_leg->cpn[j].start_date  = pd_start[i];
        exo_leg->cpn[j].pay_date    = pd_pay[i];
        exo_leg->cpn[j].fx_fix_date = pd_fix[i];
        exo_leg->cpn[j].fx_val_date =
            add_unit(exo_leg->cpn[j].fx_fix_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

        /*	Times */
        exo_leg->cpn[j].start_time  = (exo_leg->cpn[j].start_date - today) * YEARS_IN_DAY;
        exo_leg->cpn[j].pay_time    = (exo_leg->cpn[j].pay_date - today) * YEARS_IN_DAY;
        exo_leg->cpn[j].fx_fix_time = (exo_leg->cpn[j].fx_fix_date - today) * YEARS_IN_DAY;
        exo_leg->cpn[j].fx_val_time = (exo_leg->cpn[j].fx_val_date - today) * YEARS_IN_DAY;

        /*	Coupon */
        err = interp_basis(pd_basis[i], &bas);
        if (err)
            goto FREE_RETURN;

        cvg = coverage(pd_start[i], pd_pay[i], bas) * pd_not;

        if (!use_cpn_opt_str)
        {
            exo_leg->cpn[j].alpha = cvg * pd_alpha[i];
            exo_leg->cpn[j].beta  = cvg * pd_beta[i];

            exo_leg->cpn[j].floored = pd_floored[i];
            exo_leg->cpn[j].capped  = pd_capped[i];

            exo_leg->cpn[j].floor = cvg * pd_floor[i];
            exo_leg->cpn[j].cap   = cvg * pd_cap[i];
        }
        else
        {
            // Interp coupon:
            exo_leg->cpn[j].use_opt_str = 1;
            exo_leg->cpn[j].nstrikes    = pd_num_strikes[i];
            exo_leg->cpn[j].wcst        = cvg * pd_wcst[i];
            exo_leg->cpn[j].wspot       = cvg * pd_wspot[i];
            if (pd_num_strikes[i] > 0)
            {
                exo_leg->cpn[j].strikes = (double*)calloc(pd_num_strikes[i], sizeof(double));
                exo_leg->cpn[j].weights = (double*)calloc(pd_num_strikes[i], sizeof(double));
                if (!exo_leg->cpn[j].strikes || !exo_leg->cpn[j].weights)
                {
                    err = serror("Memory failure in cpd_fill_exo_leg");
                    goto FREE_RETURN;
                }
                memcpy(exo_leg->cpn[j].strikes, pd_strikes[i], pd_num_strikes[i] * sizeof(double));
                for (k = 0; k < pd_num_strikes[i]; k++)
                    exo_leg->cpn[j].weights[k] = cvg * pd_weights[i][k];
            }
        }
        i++;
        j++;
    }

DO_NOTIONAL_REF:
    /*	Notional refund */

    /*	Dates */
    exo_leg->not_ref.start_date  = pd_start[pd_ncpn - 1];
    exo_leg->not_ref.pay_date    = pd_pay[pd_ncpn - 1];
    exo_leg->not_ref.fx_fix_date = pd_not_ref_fix;
    exo_leg->not_ref.fx_val_date =
        add_unit(exo_leg->not_ref.fx_fix_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    /*	Times */
    exo_leg->not_ref.start_time  = (exo_leg->not_ref.start_date - today) * YEARS_IN_DAY;
    exo_leg->not_ref.pay_time    = (exo_leg->not_ref.pay_date - today) * YEARS_IN_DAY;
    exo_leg->not_ref.fx_fix_time = (exo_leg->not_ref.fx_fix_date - today) * YEARS_IN_DAY;
    exo_leg->not_ref.fx_val_time = (exo_leg->not_ref.fx_val_date - today) * YEARS_IN_DAY;

    /*	Coupon */
    if (!use_not_opt_str)
    {
        exo_leg->not_ref.alpha = pd_not * pd_not_ref_alpha;
        exo_leg->not_ref.beta  = pd_not * pd_not_ref_beta;

        exo_leg->not_ref.floored = pd_not_ref_floored;
        exo_leg->not_ref.capped  = pd_not_ref_capped;

        exo_leg->not_ref.floor = pd_not * pd_not_ref_floor;
        exo_leg->not_ref.cap   = pd_not * pd_not_ref_cap;
    }
    else
    // Interp redemption:
    {
        exo_leg->not_ref.use_opt_str = 1;
        exo_leg->not_ref.nstrikes    = pd_not_num_strikes;
        exo_leg->not_ref.wcst        = pd_not * pd_not_wcst;
        exo_leg->not_ref.wspot       = pd_not * pd_not_wspot;
        if (pd_not_num_strikes > 0)
        {
            exo_leg->not_ref.strikes = (double*)calloc(pd_not_num_strikes, sizeof(double));
            exo_leg->not_ref.weights = (double*)calloc(pd_not_num_strikes, sizeof(double));
            if (!exo_leg->not_ref.strikes || !exo_leg->not_ref.weights)
            {
                err = serror("Memory failure in cpd_not_fill_exo_leg");
                goto FREE_RETURN;
            }
            memcpy(exo_leg->not_ref.strikes, pd_not_strikes, pd_not_num_strikes * sizeof(double));
            for (k = 0; k < pd_not_num_strikes; k++)
                exo_leg->not_ref.weights[k] = pd_not * pd_not_weights[k];
        }
    }

    err = cpd_check_exo_leg(exo_leg);

FREE_RETURN:

    if (err)
        cpd_free_exo_leg(exo_leg);
    return err;
}

/*	Check dates consistency */
Err cpd_check_exo_leg(PD_EXO_LEG exo_leg)
{
    int i;

    /*	Check that start, pay, fix and val dates are increasing */
    for (i = 1; i < exo_leg->num_cpn; i++)
    {
        if (exo_leg->cpn[i].start_date < exo_leg->cpn[i - 1].start_date)
            return "Start dates should be increasing in exotic leg";

        if (exo_leg->cpn[i].pay_date < exo_leg->cpn[i - 1].pay_date)
            return "Pay dates should be increasing in exotic leg";

        if (exo_leg->cpn[i].fx_fix_date < exo_leg->cpn[i - 1].fx_fix_date)
            return "Fixing dates should be increasing in exotic leg";

        if (exo_leg->cpn[i].fx_val_date < exo_leg->cpn[i - 1].fx_val_date)
            return "Value dates should be increasing in exotic leg";
    }

    /*	Check that pay dates are after start dates and fix dates
                    and that val dates are after fix dates */
    for (i = 0; i < exo_leg->num_cpn; i++)
    {
        if (exo_leg->cpn[i].pay_date < exo_leg->cpn[i].start_date)
            return "Pay dates should be after start dates in exotic leg";

        if (exo_leg->cpn[i].pay_date < exo_leg->cpn[i].fx_fix_date)
            return "Pay dates should be after fixing dates in exotic leg";

        if (exo_leg->cpn[i].fx_val_date < exo_leg->cpn[i].fx_fix_date)
            return "Value dates should be after fixing dates in exotic leg";
    }
    /*	Check that floor is less than cap */
    for (i = 0; i < exo_leg->num_cpn; i++)
    {
        if (!exo_leg->cpn[i].use_opt_str)
        {
            if (exo_leg->cpn[i].floored && exo_leg->cpn[i].capped &&
                exo_leg->cpn[i].floor > exo_leg->cpn[i].cap)
                return "Floor should be lower than cap";
        }
    }

    /*	Same for notional */
    if (exo_leg->not_ref.pay_date < exo_leg->not_ref.start_date)
        return "Pay dates should be after start dates in exotic leg notional refund";

    if (exo_leg->not_ref.pay_date < exo_leg->not_ref.fx_fix_date)
        return "Pay dates should be after fixing dates in exotic leg notional refund";

    if (exo_leg->not_ref.fx_val_date < exo_leg->not_ref.fx_fix_date)
        return "Value dates should be after fixing dates in exotic leg notional refund";

    if (!exo_leg->not_ref.use_opt_str)
    {
        if (exo_leg->not_ref.floored && exo_leg->not_ref.capped &&
            exo_leg->not_ref.floor > exo_leg->not_ref.cap)
            return "Floor should be lower than cap";
    }

    /*	Check that notional refund starts, fixes and pays after the last coupon */
    if (exo_leg->num_cpn > 0)
    {
        if (exo_leg->not_ref.start_date < exo_leg->cpn[exo_leg->num_cpn - 1].start_date)
            return "Notional refund should start after last coupon in exotic leg";

        if (exo_leg->not_ref.pay_date < exo_leg->cpn[exo_leg->num_cpn - 1].pay_date)
            return "Notional refund should pay after last coupon in exotic leg";

        if (exo_leg->not_ref.fx_fix_date < exo_leg->cpn[exo_leg->num_cpn - 1].fx_fix_date)
            return "Notional refund should fix after last coupon in exotic leg";

        if (exo_leg->not_ref.fx_val_date < exo_leg->cpn[exo_leg->num_cpn - 1].fx_val_date)
            return "Notional refund value date should be after last coupon value date in exotic "
                   "leg";
    }

    /*	OK */
    return NULL;
}

/*	Free */
Err cpd_free_exo_leg(PD_EXO_LEG exo_leg)
{
    int i;
    if (exo_leg->cpn)
    {
        for (i = 0; i < exo_leg->num_cpn; i++)
        {
            if (exo_leg->cpn[i].use_opt_str && exo_leg->cpn[i].nstrikes > 0)
            {
                free(exo_leg->cpn[i].strikes);
                free(exo_leg->cpn[i].weights);
            }
        }
        free(exo_leg->cpn);
    }
    if (exo_leg->not_ref.use_opt_str && exo_leg->not_ref.nstrikes > 0)
    {
        free(exo_leg->not_ref.strikes);
        free(exo_leg->not_ref.weights);
    }
    memset(exo_leg, 0, sizeof(pd_exo_leg));

    return NULL;
}

/*	Functions for the TARN */
Err cpd_fill_exo_and_fund_leg_TARN(long today, CPD_STR cpd)
{
    int          j;
    pd_exo_leg*  pd_leg   = NULL;
    pd_fund_leg* fund_leg = NULL;

    Err err = NULL;

    pd_leg   = (PD_EXO_LEG)calloc(1, sizeof(pd_exo_leg));
    fund_leg = (PD_FUND_LEG)calloc(1, sizeof(pd_fund_leg));

    //	EXOTIC LEG
    pd_leg->cpn     = NULL;
    pd_leg->num_cpn = cpd->pd_leg->num_cpn + 1;
    pd_leg->cpn     = (pd_exo_cpn*)calloc(cpd->pd_leg->num_cpn + 1, sizeof(pd_exo_cpn));

    for (j = 0; j < cpd->pd_leg->num_cpn; j++)
    {
        memcpy(&(pd_leg->cpn[j]), &(cpd->pd_leg->cpn[j]), sizeof(pd_exo_cpn));

        if (cpd->pd_leg->cpn[j].use_opt_str)
        {
            // Interp coupon:
            pd_leg->cpn[j].use_opt_str = 1;
            if (pd_leg->cpn[j].nstrikes > 0)
            {
                pd_leg->cpn[j].strikes = (double*)calloc(pd_leg->cpn[j].nstrikes, sizeof(double));
                pd_leg->cpn[j].weights = (double*)calloc(pd_leg->cpn[j].nstrikes, sizeof(double));
                if (!pd_leg->cpn[j].strikes || !pd_leg->cpn[j].weights)
                {
                    err = serror("Memory failure in cpd_fill_exo_leg");
                    goto FREE_RETURN;
                }
                memcpy(
                    pd_leg->cpn[j].strikes,
                    cpd->pd_leg->cpn[j].strikes,
                    pd_leg->cpn[j].nstrikes * sizeof(double));
                memcpy(
                    pd_leg->cpn[j].weights,
                    cpd->pd_leg->cpn[j].weights,
                    pd_leg->cpn[j].nstrikes * sizeof(double));
            }
        }
    }

    /*	Dates */
    pd_leg->cpn[cpd->pd_leg->num_cpn].start_date =
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date;
    pd_leg->cpn[cpd->pd_leg->num_cpn].pay_date =
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date;
    pd_leg->cpn[cpd->pd_leg->num_cpn].fx_fix_date =
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date;
    pd_leg->cpn[cpd->pd_leg->num_cpn].fx_val_date = add_unit(
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    /*	Times */
    pd_leg->cpn[cpd->pd_leg->num_cpn].start_time =
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_time;
    pd_leg->cpn[cpd->pd_leg->num_cpn].pay_time =
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_time;
    pd_leg->cpn[cpd->pd_leg->num_cpn].fx_fix_time =
        cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_time;
    pd_leg->cpn[cpd->pd_leg->num_cpn].fx_val_time =
        (pd_leg->cpn[cpd->pd_leg->num_cpn - 1].fx_val_date - today) * YEARS_IN_DAY;

    if (!cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].use_opt_str)
    {
        pd_leg->cpn[cpd->pd_leg->num_cpn].use_opt_str = 0;
        pd_leg->cpn[cpd->pd_leg->num_cpn].alpha       = 0.0;
        pd_leg->cpn[cpd->pd_leg->num_cpn].beta        = 0.0;

        pd_leg->cpn[cpd->pd_leg->num_cpn].floored = 0;
        pd_leg->cpn[cpd->pd_leg->num_cpn].capped  = 0;

        pd_leg->cpn[cpd->pd_leg->num_cpn].floor = 0.0;
        pd_leg->cpn[cpd->pd_leg->num_cpn].cap   = 0.0;
    }
    else
    {
        pd_leg->cpn[cpd->pd_leg->num_cpn].use_opt_str = 1;
        pd_leg->cpn[cpd->pd_leg->num_cpn].nstrikes    = 0;
        pd_leg->cpn[cpd->pd_leg->num_cpn].wcst        = 0.0;
        pd_leg->cpn[cpd->pd_leg->num_cpn].wspot       = 0.0;
    }

    /*	Notional refund */
    memcpy(&(pd_leg->not_ref), &(cpd->pd_leg->not_ref), sizeof(pd_exo_cpn));

    if (cpd->pd_leg->not_ref.use_opt_str)
    {
        pd_leg->not_ref.use_opt_str = 1;
        if (pd_leg->not_ref.nstrikes > 0)
        {
            pd_leg->not_ref.strikes = (double*)calloc(pd_leg->not_ref.nstrikes, sizeof(double));
            pd_leg->not_ref.weights = (double*)calloc(pd_leg->not_ref.nstrikes, sizeof(double));
            if (!pd_leg->not_ref.strikes || !pd_leg->not_ref.weights)
            {
                err = serror("Memory failure in cpd_not_fill_exo_leg");
                goto FREE_RETURN;
            }
            memcpy(
                pd_leg->not_ref.strikes,
                cpd->pd_leg->not_ref.strikes,
                pd_leg->not_ref.nstrikes * sizeof(double));
            memcpy(
                pd_leg->not_ref.weights,
                cpd->pd_leg->not_ref.weights,
                pd_leg->not_ref.nstrikes * sizeof(double));
        }
    }

    //	FUNDING LEG
    fund_leg->dom_for  = cpd->fund_leg->dom_for;
    fund_leg->notional = cpd->fund_leg->notional;
    fund_leg->cpn      = NULL;
    fund_leg->num_cpn  = cpd->fund_leg->num_cpn + 1;
    fund_leg->cpn      = (pd_fund_cpn*)calloc(cpd->fund_leg->num_cpn + 1, sizeof(pd_fund_cpn));

    if (!fund_leg->cpn)
    {
        err = "Allocation error in cpd_fill_fund_leg";
        goto FREE_RETURN;
    }

    /*	Fill coupons information */
    for (j = 0; j < fund_leg->num_cpn; j++)
        memcpy(&(fund_leg->cpn[j]), &(cpd->fund_leg->cpn[j]), sizeof(pd_fund_cpn));

    /*	Dates */
    fund_leg->cpn[fund_leg->num_cpn].start_date =
        cpd->fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;
    fund_leg->cpn[fund_leg->num_cpn].pay_date = cpd->fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;

    /*	Times */
    fund_leg->cpn[fund_leg->num_cpn].start_time =
        cpd->fund_leg->cpn[fund_leg->num_cpn - 1].pay_time;
    fund_leg->cpn[fund_leg->num_cpn].pay_time = cpd->fund_leg->cpn[fund_leg->num_cpn - 1].pay_time;

    /*	Coupon */
    fund_leg->cpn[fund_leg->num_cpn].cpn = 0.0;

    // UPDATE THE CPD STRUCTURE
    cpd_free_exo_leg(cpd->pd_leg);
    cpd_free_fund_leg(cpd->fund_leg);
    cpd->pd_leg   = pd_leg;
    cpd->fund_leg = fund_leg;

FREE_RETURN:
    if (err)
        cpd_free_exo_leg(pd_leg);
    if (err)
        cpd_free_fund_leg(fund_leg);
    return err;
}

/*	Functions for the calls */

Err cpd_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int   eod_flag, /*	0: I, 1: E */
    int   ncall,
    int*  type,    /*	0: call, 1: KO */
    int   pay_rec, /*	0: rec pd, 1: pay pd */
    long* ex_date,
    long* set_date,
    long* pd_not_ref_fix,  // Fx fixing dates in case of exotic redemption at call date
    /*	Notionals on both legs in their own currencies,
                    to be refunded upon call
                    in order to determine the Bond Strike */
    double* fund_not,
    double* pd_not,
    //		pd interp notional:
    int      use_not_opt_str,
    int*     pd_not_num_strikes,
    double*  pd_not_wcst,
    double*  pd_not_wspot,
    double** pd_not_strikes,
    double** pd_not_weights,

    double* barrier,  /*	KO only */
    int*    bar_type, /*	0: up and in, 1: down and in */
    double* fees,     /*  fees if deal is called in domestic currency */
    int     TARN_Do,
    double  smooth, /*	Smoothing factor */
    int     fx_bound,
    int     use_bound,
    int     force_optim,
    CPD_STR cpd)
{
    int i, j, k;
    Err err = NULL;

    /*	Initialise pointers to NULL */
    cpd->call = NULL;

    /*	Skip calls to be exercised before today */
    i = 0;
    while (i < ncall && ex_date[i] < today + eod_flag)
        i++;

    /*	Check that at least one call is left */
    if (i == ncall)
    {
        err = "All calls are to be exercised before today in cpd_fill_calls";
        goto FREE_RETURN;
    }

    /*	Allocate memory */
    cpd->num_calls = ncall - i;
    cpd->call      = (pd_call*)calloc(cpd->num_calls, sizeof(pd_call));
    if (!cpd->num_calls)
    {
        err = "Allocation error in cpd_fill_calls";
        goto FREE_RETURN;
    }

    /*	Fill calls information */
    j         = 0;
    cpd->type = type[i];

    while (i < ncall)
    {
        /*	Dates */
        cpd->call[j].ex_date    = ex_date[i];
        cpd->call[j].ex_date2bd = add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING);
        cpd->call[j].set_date   = set_date[i];

        /*	Times */
        cpd->call[j].ex_time    = (cpd->call[j].ex_date - today) * YEARS_IN_DAY;
        cpd->call[j].ex_time2bd = (cpd->call[j].ex_date2bd - today) * YEARS_IN_DAY;
        cpd->call[j].set_time   = (cpd->call[j].set_date - today) * YEARS_IN_DAY;

        /*	Call on funding leg */
        /*	k = index of the first coupon to be called on funding leg,
                        i.e. first coupon with a start date >= ex date */
        k = 0;
        while (k < cpd->fund_leg->num_cpn &&
               cpd->fund_leg->cpn[k].start_date < cpd->call[j].ex_date)
            k++;
        if (k == cpd->fund_leg->num_cpn)
        {
            err = serror("Call number %d does not control any coupon in funding leg", i);
            goto FREE_RETURN;
        }
        cpd->call[j].fund_idx     = k;
        cpd->call[j].num_fund_cpn = cpd->fund_leg->num_cpn - k;

        /*	Call on exotic leg */
        /*	k = index of the first coupon to be called on exotic leg,
                        i.e. first coupon with a start date >= ex date */
        k = 0;
        while (k < cpd->pd_leg->num_cpn && cpd->pd_leg->cpn[k].start_date < cpd->call[j].ex_date)
            k++;
        if (k == cpd->pd_leg->num_cpn)
        {
            err = serror("Call number %d does not control any coupon in exotic leg", i);
            goto FREE_RETURN;
        }
        cpd->call[j].pd_idx     = k;
        cpd->call[j].num_pd_cpn = cpd->pd_leg->num_cpn - k;

        /*	Payer or receiver */
        cpd->call[j].pay_rec  = pay_rec;
        cpd->call[j].fee      = fees[i];
        cpd->call[j].orig_fee = fees[i];

        if (type[i] == 1 && !force_optim)
        /*	KO */
        {
            /* TARN */
            cpd->call[j].TARN_Do = TARN_Do;

            cpd->call[j].call_type    = 1;
            cpd->call[j].barrier      = barrier[i];
            cpd->call[j].orig_barrier = barrier[i];
            cpd->call[j].smooth       = smooth;
            cpd->call[j].bar_type     = bar_type[i];
            cpd->call[j].fx_bound     = 1;
        }
        else
        {
            cpd->call[j].fx_bound = fx_bound; /*	Index for optimisation */

            if (use_bound)
            {
                /* if use bound is equal to 1 we turn the call features into UO features */
                cpd->call[j].call_type    = 1; /* turn to call */
                cpd->call[j].bar_type     = 0; /* turn to UO */
                cpd->call[j].barrier      = barrier[i];
                cpd->call[j].orig_barrier = barrier[i];
            }
            else
            {
                cpd->call[j].call_type = 0;
            }
        }

        if (type[i] != cpd->type || force_optim)
        {
            cpd->type = -1;
        }

        /*	Notionals */
        cpd->call[j].fund_not_amt = fund_not[i];
        cpd->call[j].pd_not_amt   = pd_not[i];

        // Exotic redemption:
        if (use_not_opt_str)
        {
            cpd->call[j].fx_fix_date = pd_not_ref_fix[j];
            cpd->call[j].fx_fix_time = (cpd->call[j].fx_fix_date - today) * YEARS_IN_DAY;
            cpd->call[j].fx_val_date =
                add_unit(cpd->call[j].fx_fix_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
            cpd->call[j].fx_val_time = (cpd->call[j].fx_val_date - today) * YEARS_IN_DAY;
            cpd->call[j].use_opt_str = 1;
            cpd->call[j].nstrikes    = pd_not_num_strikes[i];
            cpd->call[j].wcst        = pd_not[i] * pd_not_wcst[i];
            cpd->call[j].wspot       = pd_not[i] * pd_not_wspot[i];
            if (pd_not_num_strikes[i] > 0)
            {
                cpd->call[j].strikes = (double*)calloc(pd_not_num_strikes[i], sizeof(double));
                cpd->call[j].weights = (double*)calloc(pd_not_num_strikes[i], sizeof(double));
                if (!cpd->call[j].strikes || !cpd->call[j].weights)
                {
                    err = serror("Memory failure in cpd_fill_calls");
                    goto FREE_RETURN;
                }
                memcpy(
                    cpd->call[j].strikes,
                    pd_not_strikes[i],
                    pd_not_num_strikes[i] * sizeof(double));
                for (k = 0; k < pd_not_num_strikes[i]; k++)
                    cpd->call[j].weights[k] = pd_not[i] * pd_not_weights[i][k];
            }
        }
        i++;
        j++;
    }

    if (use_bound)
        cpd->type = 1;

    err = cpd_check_calls(cpd);

FREE_RETURN:

    if (err)
        cpd_free_calls(cpd);

    return err;
}

/*	Check dates consistency */
Err cpd_check_calls(CPD_STR cpd)
{
    int i;

    /*	Check that ex and set dates are strictly increasing
                    Also check that funding and pd indices are strictly increasing,
                    i.e. there is no redundant calls */
    for (i = 1; i < cpd->num_calls; i++)
    {
        if (cpd->call[i].ex_date <= cpd->call[i - 1].ex_date)
            return "Exercise dates should be increasing";

        if (cpd->call[i].set_date <= cpd->call[i - 1].set_date)
            return "Settlement dates should be increasing";

        if (cpd->call[i].fund_idx < cpd->call[i - 1].fund_idx)
            return "Number of funding coupons controlled by calls should be decreasing";

        if (cpd->call[i].pd_idx < cpd->call[i - 1].pd_idx)
            return "Number of exotic coupons controlled by calls should be decreasing";

        if (cpd->call[i].fund_idx <= cpd->call[i - 1].fund_idx &&
            cpd->call[i].pd_idx <= cpd->call[i - 1].pd_idx)
            return serror("Calls %d and %d -indexed after today- are redundant", i - 1, i);
    }

    /*	Check that set dates are after ex dates
                    Also check that the call date is before the start, end and fixing dates
                    of the coupons it controls */
    for (i = 0; i < cpd->num_calls; i++)
    {
        if (cpd->call[i].set_date < cpd->call[i].ex_date)
            return "Settlement dates should be after exercise dates";

        if (cpd->call[i].use_opt_str && cpd->call[i].fx_fix_date < cpd->call[i].ex_date)
            return serror("Fx fixing date cannot be before the corresponding exercise date");

        if (cpd->fund_leg->cpn[cpd->call[i].fund_idx].start_date < cpd->call[i].ex_date ||
            cpd->fund_leg->cpn[cpd->call[i].fund_idx].pay_date < cpd->call[i].ex_date)
            return "A funding coupon starts before its exercise date";

        if (cpd->pd_leg->cpn[cpd->call[i].pd_idx].start_date < cpd->call[i].ex_date ||
            cpd->pd_leg->cpn[cpd->call[i].pd_idx].pay_date < cpd->call[i].ex_date ||
            cpd->pd_leg->cpn[cpd->call[i].pd_idx].fx_fix_date < cpd->call[i].ex_date)
            return "An exotic coupon starts or fixes before its exercise date";
    }

    /*	Check that the last call date is before the exotic notional
                    start, pay, and fixing dates */
    if (cpd->num_calls > 0)
    {
        if (cpd->call[cpd->num_calls - 1].ex_date > cpd->pd_leg->not_ref.start_date ||
            cpd->call[cpd->num_calls - 1].ex_date > cpd->pd_leg->not_ref.pay_date ||
            cpd->call[cpd->num_calls - 1].ex_date > cpd->pd_leg->not_ref.fx_fix_date)
            return serror("Last call date should be before the exotic notional refund date");
    }
    /*	OK */
    return NULL;
}

/*	Free */
Err cpd_free_calls(CPD_STR cpd)
{
    int i;
    if (cpd->call)
    {
        for (i = 0; i < cpd->num_calls; i++)
        {
            if (cpd->call[i].use_opt_str && cpd->call[i].nstrikes > 0)
            {
                free(cpd->call[i].strikes);
                free(cpd->call[i].weights);
            }
        }
        free(cpd->call);
        cpd->call = NULL;
    }

    return NULL;
}

/*	Functions for the underlying */

/*	Fill underlying structure from a predefined underlying */
Err cpd_fill_und(
    char*            fx3dund,
    CPD_UND          und,
    CPDBETADLMPARAMS cpd_dlm_params,
    double           dom_vol_shift,
    double           for_vol_shift,
    double           fx_vol_shift)
{
    double *sig_time_dom = NULL, *sig_dom = NULL, *tau_time_dom = NULL, *tau_dom = NULL,
           *sig_time_for = NULL, *sig_for = NULL, *tau_time_for = NULL, *tau_for = NULL;
    long      tau_n_dom, tau_n_for, sig_n_dom, sig_n_for;
    SrtUndPtr srtund;
    char *    dom_und_name, *for_und_name;
    SrtUndPtr dom_und, for_und;
    int       i, index;
    Err       err = NULL;

    /*	Initialise */
    und->sigma_date_rates = NULL;
    und->sigma_time_rates = NULL;
    und->sigma_dom        = NULL;
    und->sigma_for        = NULL;
    und->sigma_date_fx    = NULL;
    und->sigma_time_fx    = NULL;
    und->sigma_fx         = NULL;
    und->corr_times       = NULL;
    und->correl_dom_for   = NULL;
    und->correl_dom_fx    = NULL;
    und->correl_for_fx    = NULL;
    und->model            = NULL;
    und->fees             = NULL;
    und->fees_dates       = NULL;
    und->nb_fees          = 0;

    strcpy(und->name, fx3dund);

    if (cpd_dlm_params && cpd_dlm_params->use_beta_dlm)
    {
        und->use_beta_dlm = 1;

        und->model   = (FxBetaDLM_model*)calloc(1, sizeof(FxBetaDLM_model));
        und->hermite = (FxBetaDLM_Hermite*)calloc(1, sizeof(FxBetaDLM_Hermite));

        if (!und->model || !und->hermite)
        {
            err = "Allocation error (1) in cpd_fill_und";
            goto FREE_RETURN;
        }

        err = FxBetaDLM_Get_Model(fx3dund, und->model);

        if (err)
            goto FREE_RETURN;

        und->lda_dom = und->model->dLambdaDom;
        und->lda_for = und->model->dLambdaFor;
        und->today   = und->model->lToday;
        und->spot_fx = und->model->dSpotFx;
        strcpy(und->dom_yc, und->model->cYcDom);
        strcpy(und->for_yc, und->model->cYcFor);

        und->num_params  = cpd_dlm_params->NumParams;
        und->grfn_params = cpd_dlm_params->GrfnParams;

        err = initialise_FxBetaDLM_Hermite(und->num_params, und->hermite);

        if (err)
            goto FREE_RETURN;

        /* Fill the equivalent 3 Factor for the moment */
        und->sigma_n_fx    = und->model->str->nb_fx;
        und->sigma_time_fx = (double*)calloc(und->sigma_n_fx, sizeof(double));

        if (!und->sigma_time_fx)
        {
            err = "Allocation error (1.5) in cpd_fill_und";
            goto FREE_RETURN;
        }

        memcpy(und->sigma_time_fx, und->model->str->time_fx, und->sigma_n_fx * sizeof(double));

        err = FxBetaDLM_GetEqui3FactorTS(
            und->model,
            und->num_params,
            und->hermite,
            und->sigma_n_fx,
            und->sigma_time_fx,
            &(und->sigma_fx));

        if (err)
            goto FREE_RETURN;

        und->model->dSigmaFx3F = (double*)calloc(und->model->iNbPWTime, sizeof(double));
        for (i = 0; i < und->model->iNbPWTime; i++)
        {
            index = Get_Index(und->model->dPWTime[i], und->sigma_time_fx, und->sigma_n_fx);
            und->model->dSigmaFx3F[i] = und->sigma_time_fx[index];
        }

        err = merge_rates_ts(
            und->model->str->time_dom,
            und->model->str->sig_dom,
            und->model->str->nb_dom,
            und->model->str->time_for,
            und->model->str->sig_for,
            und->model->str->nb_for,
            &(und->sigma_time_rates),
            &(und->sigma_dom),
            &(und->sigma_for),
            (int*)&(und->sigma_n_rates));

        if (err)
            goto FREE_RETURN;

        und->sigma_date_rates = (double*)calloc(und->sigma_n_rates, sizeof(double));
        und->sigma_date_fx    = (double*)calloc(und->sigma_n_fx, sizeof(double));

        und->corr_n_times   = und->model->str->nb_3F_corr;
        und->corr_times     = (double*)calloc(und->model->str->nb_3F_corr, sizeof(double));
        und->correl_dom_for = (double*)calloc(und->model->str->nb_3F_corr, sizeof(double));
        und->correl_dom_fx  = (double*)calloc(und->model->str->nb_3F_corr, sizeof(double));
        und->correl_for_fx  = (double*)calloc(und->model->str->nb_3F_corr, sizeof(double));

        if (!und->sigma_date_rates || !und->sigma_date_fx || !und->corr_times ||
            !und->correl_dom_for || !und->correl_dom_fx || !und->correl_for_fx)
        {
            err = "Allocation error (2) in cpd_fill_und";
            goto FREE_RETURN;
        }

        memcpy(
            und->corr_times,
            und->model->str->time_3F_corr,
            und->model->str->nb_3F_corr * sizeof(double));
        memcpy(
            und->correl_dom_for,
            und->model->str->dom_for_3F_corr,
            und->model->str->nb_3F_corr * sizeof(double));
        memcpy(
            und->correl_dom_fx,
            und->model->str->dom_fx_3F_corr,
            und->model->str->nb_3F_corr * sizeof(double));
        memcpy(
            und->correl_for_fx,
            und->model->str->for_fx_3F_corr,
            und->model->str->nb_3F_corr * sizeof(double));

        for (i = 0; i < und->sigma_n_rates; i++)
        {
            und->sigma_date_rates[i] = und->today + und->sigma_time_rates[i] * DAYS_IN_YEAR;
        }

        for (i = 0; i < und->sigma_n_fx; i++)
        {
            und->sigma_date_fx[i] = und->today + und->sigma_time_fx[i] * DAYS_IN_YEAR;
        }

        /* Shift the values */
        for (i = 0; i < und->model->iNbPWTime; i++)
        {
            und->model->dSigmaDom[i] += dom_vol_shift;
            und->model->dSigmaFor[i] += for_vol_shift;
            und->model->dSigmaFx[i] += fx_vol_shift;
        }

        for (i = 0; i < und->sigma_n_rates; i++)
        {
            und->sigma_dom[i] += dom_vol_shift;
            und->sigma_for[i] += for_vol_shift;
        }

        for (i = 0; i < und->sigma_n_fx; i++)
        {
            und->sigma_fx[i] += fx_vol_shift;
        }
    }
    else
    {
        und->use_beta_dlm = 0;

        /*	Get term structures */
        err = Get_FX_StochRate_TermStructures_corr(
            fx3dund,
            &sig_time_dom,
            &sig_dom,
            &sig_n_dom,
            &tau_time_dom,
            &tau_dom,
            &tau_n_dom,
            &sig_time_for,
            &sig_for,
            &sig_n_for,
            &tau_time_for,
            &tau_for,
            &tau_n_for,
            &(und->sigma_time_fx),
            &(und->sigma_fx),
            &(und->sigma_n_fx),
            &(und->corr_times),
            &(und->correl_dom_for),
            &(und->correl_dom_fx),
            &(und->correl_for_fx),
            &(und->corr_n_times));

        if (err)
        {
            goto FREE_RETURN;
        }

        err = get_unique_lambda(tau_dom, tau_n_dom, &(und->lda_dom));
        if (err)
        {
            goto FREE_RETURN;
        }

        err = get_unique_lambda(tau_for, tau_n_for, &(und->lda_for));
        if (err)
        {
            goto FREE_RETURN;
        }

        err = merge_rates_ts(
            sig_time_dom,
            sig_dom,
            sig_n_dom,
            sig_time_for,
            sig_for,
            sig_n_for,
            &(und->sigma_time_rates),
            &(und->sigma_dom),
            &(und->sigma_for),
            (int*)&(und->sigma_n_rates));

        if (err)
        {
            goto FREE_RETURN;
        }

        /*	Fill dates */
        srtund                = lookup_und(fx3dund);
        und->today            = get_today_from_underlying(srtund);
        und->sigma_date_rates = (double*)calloc(und->sigma_n_rates, sizeof(double));
        und->sigma_date_fx    = (double*)calloc(und->sigma_n_fx, sizeof(double));

        if (!und->sigma_date_rates || !und->sigma_date_fx)
        {
            err = "Allocation error in cpd_fill_und";
            goto FREE_RETURN;
        }

        for (i = 0; i < und->sigma_n_rates; i++)
        {
            und->sigma_date_rates[i] =
                und->today + und->sigma_time_rates[i] * DAYS_IN_YEAR + 1.0e-08;
        }

        for (i = 0; i < und->sigma_n_fx; i++)
        {
            und->sigma_date_fx[i] = und->today + und->sigma_time_fx[i] * DAYS_IN_YEAR + 1.0e-08;
        }

        /*	Spot fx, term structures and yield curves */

        dom_und_name = get_domname_from_fxund(srtund);
        dom_und      = lookup_und(dom_und_name);
        strcpy(und->dom_yc, get_ycname_from_irund(dom_und));

        for_und_name = get_forname_from_fxund(srtund);
        for_und      = lookup_und(for_und_name);
        strcpy(und->for_yc, get_ycname_from_irund(for_und));

        err = Fx3DFwdFx(fx3dund, und->today, &(und->spot_fx));

        /*	Shifts */

        for (i = 0; i < und->sigma_n_rates; i++)
        {
            und->sigma_dom[i] += dom_vol_shift;
            und->sigma_for[i] += for_vol_shift;
        }

        for (i = 0; i < und->sigma_n_fx; i++)
        {
            und->sigma_fx[i] += fx_vol_shift;
        }
    }

FREE_RETURN:

    if (err)
    {
        cpd_free_und(und);
    }

    if (tau_time_dom)
    {
        free(tau_time_dom);
    }

    if (tau_dom)
    {
        free(tau_dom);
    }

    if (tau_time_for)
    {
        free(tau_time_for);
    }

    if (tau_for)
    {
        free(tau_for);
    }

    if (sig_time_dom)
    {
        free(sig_time_dom);
    }

    if (sig_dom)
    {
        free(sig_dom);
    }

    if (sig_time_for)
    {
        free(sig_time_for);
    }

    if (sig_for)
    {
        free(sig_for);
    }

    return err;
}

/*	Fill underlying structure from calibration instruments */
Err cpd_calib_und(
    long today,
    /*	EOD Flag */
    int    eod_flag, /*	0: I, 1: E */
    double fx_spot,
    long   fx_spot_date,
    int    dom_calib,      /*	Calibrate domestic underlying */
    char*  dom_und,        /*	If no, domestic underlying to be used */
    char*  dom_yc,         /*	Domestic yc */
    char*  dom_vc,         /*	Domestic vc (only if calib) */
    char*  dom_ref,        /*	Domestic ref rate (only if calib) */
    char*  dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char*  dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,        /*	Domestic lambda */
    int    for_calib,      /*	Same for foreign */
    char*  for_und,
    char*  for_yc,
    char*  for_vc,
    char*  for_ref,
    char*  for_swap_freq,
    char*  for_swap_basis,
    double for_lam,

    double min_fact,  /*	Maximum down jump on variance */
    double max_fact,  /*	Maximum up jump on variance */
    int    use_jumps, /*	Allow vol term structure to jump */

    double*          corr_times,
    double*          correl_dom_for, /*	Correlations */
    double*          correl_dom_fx,
    double*          correl_for_fx,
    long             corr_n_times,
    CPDBETADLMPARAMS cpd_dlm_params,
    CPD_STR          cpd,  /*	Structure */
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char*   vol_curve_name,
                           double  start_date,
                           double  end_date,
                           double  cash_strike,
                           int     zero,
                           char*   ref_rate_name,
                           double* vol,
                           double* power),
    /*	Fx vol from the market */
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,
    CPD_UND und,
    double  dom_vol_shift,
    double  for_vol_shift,
    double  fx_vol_shift)
{
    int   i, i0, nex;
    long *lex = NULL, last;

    double* fx_mkt_vol_time = NULL;

    double *temp_sig_time1 = NULL, *temp_sig_time2 = NULL, *temp_sig1 = NULL, *temp_sig2 = NULL;

    double *temp_tau_time = NULL, *temp_tau = NULL;

    int temp_tau_n, temp_sig_n1, temp_sig_n2;

    double* opt_strikes     = NULL;
    long*   fx_mkt_fix_date = NULL;

    Err err = NULL;

    /*	Eliminate zero/negative lambdas */

    if (dom_lam < 1.0e-08)
        dom_lam = 1.0e-08;
    if (for_lam < 1.0e-08)
        for_lam = 1.0e-08;

    /*	Initialise */

    memset(und, 0, sizeof(cpd_und));
    und->today = today;

    und->spot_fx =
        fx_spot * swp_f_df(today, fx_spot_date, dom_yc) / swp_f_df(today, fx_spot_date, for_yc);

    und->corr_n_times   = corr_n_times;
    und->corr_times     = (double*)calloc(corr_n_times, sizeof(double));
    und->correl_dom_for = (double*)calloc(corr_n_times, sizeof(double));
    und->correl_dom_fx  = (double*)calloc(corr_n_times, sizeof(double));
    und->correl_for_fx  = (double*)calloc(corr_n_times, sizeof(double));

    if (!und->corr_times || !und->correl_dom_for || !und->correl_dom_fx || !und->correl_for_fx)
    {
        err = "Allocation error (1) in cpd_calib_und";
        goto FREE_RETURN;
    }

    memcpy(und->corr_times, corr_times, corr_n_times * sizeof(double));
    memcpy(und->correl_dom_for, correl_dom_for, corr_n_times * sizeof(double));
    memcpy(und->correl_dom_fx, correl_dom_fx, corr_n_times * sizeof(double));
    memcpy(und->correl_for_fx, correl_for_fx, corr_n_times * sizeof(double));

    strcpy(und->dom_yc, dom_yc);
    strcpy(und->for_yc, for_yc);
    strcpy(und->name, "CALIB");

    /*	Find exercise dates for the calibration of interest rates */

    if (cpd->pd_leg->num_cpn > 0)
        last = cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date;
    else
        last = 0;

    if (last < today + 740)
        last = today + 740;

    if (dom_calib || for_calib)
    {
        /* check that we don't calibrate too short options */
        i0 = 0;
        while (i0 < cpd->num_calls && (cpd->call[i0].ex_date - today) < CPD_MINCALDAYS)
            i0++;

        if (cpd->num_calls - i0 > 0 &&
            !(cpd->num_calls - i0 == 1 && cpd->call[i0].ex_date <= today + eod_flag))
        {
            /*	If call dates, choose call dates as option expiries for calibration */
            nex = cpd->num_calls - i0;
            lex = (long*)calloc(nex, sizeof(long));
            if (!lex)
            {
                err = "Allocation error (1) in cpd_calib_und";
                goto FREE_RETURN;
            }
            for (i = 0; i < nex; i++)
                lex[i] = cpd->call[i + i0].ex_date;
        }
        else
        {
            /* Calibration every year up to last */
            nex = (last - today) / 365;

            /* check that we don't calibrate too short options */
            if (((last - nex * 365 - today) < CPD_MINCALDAYS) && (nex > 1))
                nex--;

            lex = (long*)calloc(nex, sizeof(long));
            if (!lex)
            {
                err = "Allocation error (1') in cpd_calib_und";
                goto FREE_RETURN;
            }

            for (i = 0; i < nex; i++)
                lex[i] = last - (nex - i) * 365;
        }
    }

    /*	Domestic underlying */

    if (dom_calib)
    {
        und->lda_dom = dom_lam;

        err = cpd_calib_diagonal(
            dom_yc,
            dom_vc,
            dom_ref,
            get_ir_cash_vol,
            dom_vol_shift,
            0,
            nex,
            lex,
            last,
            NULL,
            NULL,
            0,
            1.0,
            1.0,
            dom_swap_freq,
            dom_swap_basis,
            1,
            0,
            1,
            CALPRES,
            CALPRES,
            min_fact,
            max_fact,
            use_jumps,
            0,
            NULL,
            &dom_lam,
            1,
            0.0,
            0.0,
            0.0,
            &temp_sig_n1,
            &temp_sig_time1,
            &temp_sig1,
            NULL);

        if (err && cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date < today + 740)
        {
            err = cpd_calib_diagonal(
                dom_yc,
                dom_vc,
                dom_ref,
                get_ir_cash_vol,
                dom_vol_shift,
                0,
                nex,
                lex,
                last,
                NULL,
                NULL,
                0,
                1.0,
                1.0,
                "Q",
                dom_swap_basis,
                1,
                0,
                1,
                CALPRES,
                CALPRES,
                min_fact,
                max_fact,
                use_jumps,
                0,
                NULL,
                &dom_lam,
                1,
                0.0,
                0.0,
                0.0,
                &temp_sig_n1,
                &temp_sig_time1,
                &temp_sig1,
                NULL);

            if (err && strcmp(err, "All exercise dates are past in cpd_calib_diagonal"))
            {
                temp_sig_n1 = 1;

                temp_sig_time1 = (double*)calloc(temp_sig_n1, sizeof(double));
                temp_sig1      = (double*)calloc(temp_sig_n1, sizeof(double));

                if (!temp_sig_time1 || !temp_sig1)
                {
                    err = "cpd_calib_und: memory allocation failure";
                    goto FREE_RETURN;
                }

                temp_sig_time1[0] = 1.0;
                temp_sig1[0]      = 0.01;

                err = NULL;
            }
        }

        if (err)
            goto FREE_RETURN;
    }
    else
    {
        err = Get_LGM_TermStructure(
            dom_und,
            &temp_sig_time1,
            &temp_sig1,
            (long*)&temp_sig_n1,
            &temp_tau_time,
            &temp_tau,
            (long*)&temp_tau_n);
        if (err)
            goto FREE_RETURN;

        err = get_unique_lambda(temp_tau, temp_tau_n, &(und->lda_dom));
        if (err)
            goto FREE_RETURN;
        free(temp_tau_time);
        temp_tau_time = NULL;
        free(temp_tau);
        temp_tau = NULL;

        for (i = 0; i < temp_sig_n1; i++)
            temp_sig1[i] += dom_vol_shift;
    }

    /*	Foreign underlying */

    if (for_calib)
    {
        und->lda_for = for_lam;

        err = cpd_calib_diagonal(
            for_yc,
            for_vc,
            for_ref,
            get_ir_cash_vol,
            for_vol_shift,
            0,
            nex,
            lex,
            last,
            NULL,
            NULL,
            0,
            1.0,
            1.0,
            for_swap_freq,
            for_swap_basis,
            1,
            0,
            1,
            CALPRES,
            CALPRES,
            min_fact,
            max_fact,
            use_jumps,
            0,
            NULL,
            &for_lam,
            1,
            0.0,
            0.0,
            0.0,
            &temp_sig_n2,
            &temp_sig_time2,
            &temp_sig2,
            NULL);

        if (err && cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date < today + 740)
        {
            err = cpd_calib_diagonal(
                for_yc,
                for_vc,
                for_ref,
                get_ir_cash_vol,
                for_vol_shift,
                0,
                nex,
                lex,
                last,
                NULL,
                NULL,
                0,
                1.0,
                1.0,
                "Q",
                for_swap_basis,
                1,
                0,
                1,
                CALPRES,
                CALPRES,
                min_fact,
                max_fact,
                use_jumps,
                0,
                NULL,
                &for_lam,
                1,
                0.0,
                0.0,
                0.0,
                &temp_sig_n2,
                &temp_sig_time2,
                &temp_sig2,
                NULL);

            if (err && strcmp(err, "All exercise dates are past in cpd_calib_diagonal"))
            {
                temp_sig_n2 = 1;

                temp_sig_time2 = (double*)calloc(temp_sig_n2, sizeof(double));
                temp_sig2      = (double*)calloc(temp_sig_n2, sizeof(double));

                if (!temp_sig_time2 || !temp_sig2)
                {
                    err = "cpd_calib_und: memory allocation failure";
                    goto FREE_RETURN;
                }

                temp_sig_time2[0] = 1.0;
                temp_sig2[0]      = 0.01;

                err = NULL;
            }
        }

        if (err)
            goto FREE_RETURN;
    }
    else
    {
        err = Get_LGM_TermStructure(
            for_und,
            &temp_sig_time2,
            &temp_sig2,
            (long*)&temp_sig_n2,
            &temp_tau_time,
            &temp_tau,
            (long*)&temp_tau_n);
        if (err)
            goto FREE_RETURN;

        err = get_unique_lambda(temp_tau, temp_tau_n, &(und->lda_for));
        if (err)
            goto FREE_RETURN;
        free(temp_tau_time);
        temp_tau_time = NULL;
        free(temp_tau);
        temp_tau = NULL;

        for (i = 0; i < temp_sig_n2; i++)
            temp_sig2[i] += for_vol_shift;
    }

    /*	Merge rates structures */

    err = merge_rates_ts(
        temp_sig_time1,
        temp_sig1,
        temp_sig_n1,
        temp_sig_time2,
        temp_sig2,
        temp_sig_n2,
        &(und->sigma_time_rates),
        &(und->sigma_dom),
        &(und->sigma_for),
        (int*)&(und->sigma_n_rates));

    if (err)
        goto FREE_RETURN;

    und->sigma_date_rates = (double*)calloc(und->sigma_n_rates, sizeof(double));
    if (!und->sigma_date_rates)
    {
        err = "Allocation error (2) in cpd_calib_und";
        goto FREE_RETURN;
    }

    for (i = 0; i < und->sigma_n_rates; i++)
        und->sigma_date_rates[i] = today + und->sigma_time_rates[i] * DAYS_IN_YEAR + 1.0e-16;

    /*	Calibrate Fx */

    /*	Remove unused fx vol dates */
    i = num_fx_mkt_vol - 1;
    while (i > 0 && fx_mkt_vol_date[i] > last)
        i--;

    if (i < num_fx_mkt_vol - 1 && fx_mkt_vol_date[i] < last)
        i++;
    num_fx_mkt_vol = i + 1;

    /*	Remove unused correl dates */
    i = und->corr_n_times - 1;
    while (i > 0 && und->corr_times[i] > (last - today) * YEARS_IN_DAY)
        i--;

    if (i < und->corr_n_times - 1 && und->corr_times[i] < (last - today) * YEARS_IN_DAY)
        i++;
    und->corr_n_times = i + 1;

    und->sigma_n_fx    = num_fx_mkt_vol;
    fx_mkt_vol_time    = (double*)calloc(num_fx_mkt_vol, sizeof(double));
    fx_mkt_fix_date    = (long*)calloc(num_fx_mkt_vol, sizeof(long));
    und->sigma_date_fx = (double*)calloc(num_fx_mkt_vol, sizeof(double));
    und->sigma_time_fx = (double*)calloc(num_fx_mkt_vol, sizeof(double));

    if (!fx_mkt_vol_time || !fx_mkt_fix_date || !und->sigma_date_fx || !und->sigma_time_fx)
    {
        err = "Allocation error (2) in cpd_calib_und";
        goto FREE_RETURN;
    }

    for (i = 0; i < num_fx_mkt_vol; i++)
    {
        fx_mkt_vol_time[i] = (fx_mkt_vol_date[i] - today) * YEARS_IN_DAY;
        fx_mkt_vol[i] += fx_vol_shift;
        fx_mkt_fix_date[i]    = add_unit(fx_mkt_vol_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
        und->sigma_date_fx[i] = fx_mkt_fix_date[i] * 1.0;
        und->sigma_time_fx[i] = (und->sigma_date_fx[i] - today) * YEARS_IN_DAY;
    }

    err = Fx3DtsCalibration_corr(
        und->sigma_time_fx,
        fx_mkt_vol_time,
        fx_mkt_vol,
        num_fx_mkt_vol,
        und->sigma_time_rates,
        und->sigma_n_rates,
        und->sigma_dom,
        und->lda_dom,
        und->sigma_for,
        und->lda_for,
        und->corr_times,
        und->correl_dom_for,
        und->correl_dom_fx,
        und->correl_for_fx,
        und->corr_n_times,
        &(und->sigma_fx));

    if (err)
        goto FREE_RETURN;

    if (cpd_dlm_params->use_beta_dlm)
    {
        und->use_beta_dlm = 1;

        /* Create the model */
        und->model   = (FxBetaDLM_model*)calloc(1, sizeof(FxBetaDLM_model));
        und->hermite = (FxBetaDLM_Hermite*)calloc(1, sizeof(FxBetaDLM_Hermite));
        opt_strikes  = (double*)calloc(num_fx_mkt_vol, sizeof(double));

        if (!und->model || !und->hermite || !opt_strikes)
        {
            err = err = "Allocation error (3) in cpd_calib_und";
            goto FREE_RETURN;
        }

        und->num_params  = cpd_dlm_params->NumParams;
        und->grfn_params = cpd_dlm_params->GrfnParams;

        /* Initialisation of the model */
        err = init_FxBetaDLM_model(
            1.0,                            // MinTimeCorrel
            (last - today) * YEARS_IN_DAY,  // MaxTimeCorrel
            today,
            dom_yc,
            1.0 / und->lda_dom,
            und->sigma_n_rates,
            und->sigma_time_rates,
            und->sigma_dom,
            for_yc,
            1.0 / und->lda_for,
            und->sigma_n_rates,
            und->sigma_time_rates,
            und->sigma_for,
            fx_spot,
            cpd_dlm_params->tstar,
            cpd_dlm_params->B0,
            cpd_dlm_params->C0,
            cpd_dlm_params->alpha,
            cpd_dlm_params->lambda,
            num_fx_mkt_vol,
            und->sigma_time_fx,
            und->sigma_fx,
            und->sigma_fx,
            und->corr_n_times,
            und->corr_times,
            und->correl_dom_for,
            und->correl_dom_fx,
            und->correl_for_fx,
            und->model);

        if (err)
            goto FREE_RETURN;

        err = initialise_FxBetaDLM_Hermite(und->num_params, und->hermite);

        if (err)
            goto FREE_RETURN;

        for (i = 0; i < num_fx_mkt_vol; i++)
            opt_strikes[i] = 0.0;

        err = FxBetaDLM_CalibrationModel(
            num_fx_mkt_vol,
            fx_mkt_vol_date,
            fx_mkt_fix_date,
            opt_strikes,
            fx_mkt_vol,
            cpd_dlm_params->calib_smile,
            cpd_dlm_params->smile_settlmt_date,
            cpd_dlm_params->smile_strikes,
            cpd_dlm_params->smile_bs_vols,
            und->model,
            und->num_params,
            und->hermite);

        if (err)
            goto FREE_RETURN;
    }
    else
    {
        und->use_beta_dlm = 0;
    }

FREE_RETURN:

    if (err)
    {
        cpd_free_und(und);
        und->model   = NULL;
        und->hermite = NULL;
    }

    if (fx_mkt_vol_time)
        free(fx_mkt_vol_time);

    if (lex)
        free(lex);
    if (temp_sig_time1)
        free(temp_sig_time1);
    if (temp_sig1)
        free(temp_sig1);
    if (temp_sig_time2)
        free(temp_sig_time2);
    if (temp_sig2)
        free(temp_sig2);
    if (temp_tau_time)
        free(temp_tau_time);
    if (temp_tau)
        free(temp_tau);

    if (opt_strikes)
        free(opt_strikes);
    if (fx_mkt_fix_date)
        free(fx_mkt_fix_date);

    return err;
}

void cpd_copy_und(CPD_UND src, CPD_UND dest)
{
    Err err = NULL;

    dest->today   = src->today;
    dest->spot_fx = src->spot_fx;

    strcpy(dest->name, src->name);
    strcpy(dest->dom_yc, src->dom_yc);
    strcpy(dest->for_yc, src->for_yc);

    dest->sigma_n_rates    = src->sigma_n_rates;
    dest->sigma_date_rates = (double*)calloc(dest->sigma_n_rates, sizeof(double));
    memcpy(dest->sigma_date_rates, src->sigma_date_rates, dest->sigma_n_rates * sizeof(double));
    dest->sigma_time_rates = (double*)calloc(dest->sigma_n_rates, sizeof(double));
    memcpy(dest->sigma_time_rates, src->sigma_time_rates, dest->sigma_n_rates * sizeof(double));
    dest->sigma_dom = (double*)calloc(dest->sigma_n_rates, sizeof(double));
    memcpy(dest->sigma_dom, src->sigma_dom, dest->sigma_n_rates * sizeof(double));
    dest->sigma_for = (double*)calloc(dest->sigma_n_rates, sizeof(double));
    memcpy(dest->sigma_for, src->sigma_for, dest->sigma_n_rates * sizeof(double));

    dest->lda_dom = src->lda_dom;
    dest->lda_for = src->lda_for;

    dest->sigma_n_fx    = src->sigma_n_fx;
    dest->sigma_date_fx = (double*)calloc(dest->sigma_n_fx, sizeof(double));
    memcpy(dest->sigma_date_fx, src->sigma_date_fx, dest->sigma_n_fx * sizeof(double));
    dest->sigma_time_fx = (double*)calloc(dest->sigma_n_fx, sizeof(double));
    memcpy(dest->sigma_time_fx, src->sigma_time_fx, dest->sigma_n_fx * sizeof(double));
    dest->sigma_fx = (double*)calloc(dest->sigma_n_fx, sizeof(double));
    memcpy(dest->sigma_fx, src->sigma_fx, dest->sigma_n_fx * sizeof(double));

    dest->corr_n_times = src->corr_n_times;
    dest->corr_times   = (double*)calloc(dest->corr_n_times, sizeof(double));
    memcpy(dest->corr_times, src->corr_times, dest->corr_n_times * sizeof(double));
    dest->correl_dom_for = (double*)calloc(dest->corr_n_times, sizeof(double));
    memcpy(dest->correl_dom_for, src->correl_dom_for, dest->corr_n_times * sizeof(double));
    dest->correl_dom_fx = (double*)calloc(dest->corr_n_times, sizeof(double));
    memcpy(dest->correl_dom_fx, src->correl_dom_fx, dest->corr_n_times * sizeof(double));
    dest->correl_for_fx = (double*)calloc(dest->corr_n_times, sizeof(double));
    memcpy(dest->correl_for_fx, src->correl_for_fx, dest->corr_n_times * sizeof(double));

    if (src->use_beta_dlm)
    {
        dest->use_beta_dlm = 1;
        dest->model        = (FxBetaDLM_model*)calloc(1, sizeof(FxBetaDLM_model));
        dest->hermite      = (FxBetaDLM_Hermite*)calloc(1, sizeof(FxBetaDLM_Hermite));

        err = copy_FxBetaDLM_model(src->model, dest->model);

        err = initialise_FxBetaDLM_Hermite(src->num_params, dest->hermite);
    }
    else
    {
        dest->use_beta_dlm = 0;
    }

    /* Copy the extra fees informations */

    if (src->nb_fees > 0 && src->fees && src->fees_dates)
    {
        dest->nb_fees    = src->nb_fees;
        dest->fees       = (double*)calloc(dest->nb_fees, sizeof(double));
        dest->fees_dates = (double*)calloc(dest->nb_fees, sizeof(double));

        memcpy(dest->fees, src->fees, dest->nb_fees * sizeof(double));
        memcpy(dest->fees_dates, src->fees_dates, dest->nb_fees * sizeof(double));
    }
    else
    {
        dest->nb_fees    = 0;
        dest->fees       = NULL;
        dest->fees_dates = NULL;
    }
}

Err cpd_AlphaFudgeUpdateUnd(CPD_UND und, int UpOrDown) /* Up = 0 and Down = 1 */
{
    int i;
    Err err = NULL;

    und->model->dAlpha = 0.0;

    for (i = 0; i < und->model->iNbPWTime; i++)
    {
        if (UpOrDown)
        {
            und->model->dSigmaFx[i]   = und->model->dSigmaFxDown[i];
            und->model->dCorrDomFx[i] = und->model->dCorrDomFxDown[i];
            und->model->dCorrForFx[i] = und->model->dCorrForFxDown[i];
        }
        else
        {
            und->model->dSigmaFx[i] = und->model->dSigmaFxUp[i];
        }
    }

    if (und->sigma_fx)
    {
        free(und->sigma_fx);
        und->sigma_fx = NULL;
    }

    /* not needed anymore !!! */
    /*
    err = FxBetaDLM_GetEqui3FactorTS(	und->model,
                                                                            und->num_params,
                                                                            und->hermite,
                                                                            und->sigma_n_fx,
                                                                            und->sigma_time_fx,
                                                                            &(und->sigma_fx));
    */

    return err;
}

Err cpd_free_und(CPD_UND und)
{
    if (und->sigma_date_rates)
        free(und->sigma_date_rates);
    if (und->sigma_time_rates)
        free(und->sigma_time_rates);
    if (und->sigma_dom)
        free(und->sigma_dom);
    if (und->sigma_for)
        free(und->sigma_for);
    if (und->sigma_date_fx)
        free(und->sigma_date_fx);
    if (und->sigma_time_fx)
        free(und->sigma_time_fx);
    if (und->sigma_fx)
        free(und->sigma_fx);
    if (und->corr_times)
        free(und->corr_times);
    if (und->correl_dom_for)
        free(und->correl_dom_for);
    if (und->correl_dom_fx)
        free(und->correl_dom_fx);
    if (und->correl_for_fx)
        free(und->correl_for_fx);

    und->sigma_date_rates = NULL;
    und->sigma_time_rates = NULL;
    und->sigma_dom        = NULL;
    und->sigma_for        = NULL;
    und->sigma_date_fx    = NULL;
    und->sigma_time_fx    = NULL;
    und->sigma_fx         = NULL;
    und->corr_times       = NULL;
    und->correl_dom_for   = NULL;
    und->correl_dom_fx    = NULL;
    und->correl_for_fx    = NULL;

    if (und->use_beta_dlm)
    {
        if (und->model)
        {
            free_FxBetaDLM_model(und->model);
            free(und->model);
        }

        if (und->hermite)
        {
            free_FxBetaDLM_Hermite(und->hermite);
            free(und->hermite);
        }
    }

    und->grfn_params = NULL;
    und->num_params  = NULL;

    /* Free the extra fees informations */
    if (und->fees)
        free(und->fees);
    und->fees = NULL;
    if (und->fees_dates)
        free(und->fees_dates);
    und->fees_dates = NULL;
    und->nb_fees    = 0;

    return NULL;
}

Err cpd_fill_eval_const(
    CPD_UND und,
    CPD_STR cpd,
    /*	Index of the current call/KO */
    int            call_idx,
    CPD_EVAL_CONST eval_const)
{
    PD_CALL    call, next_call;
    PD_EXO_CPN cpn;
    long       evt_date;
    double     evt_time;
    char*      fund_yc;
    double     dom_lam, for_lam, dom_phi, for_phi;
    double     fund_lam, fund_phi;
    long       today, temp_date;
    double     temp_time;
    char *     dom_yc, *for_yc;
    int        i;
    double     temp;
    Err        err = NULL;

    eval_const->call_idx = call_idx;

    /*	Extract info from call */
    call = cpd->call + call_idx;
    if (call_idx < cpd->num_calls - 1)
    {
        next_call = cpd->call + (call_idx + 1);
    }
    evt_date = call->ex_date;
    evt_time = call->ex_time;

    /*	Extract info from und */
    today   = und->today;
    dom_yc  = (char*)(und->dom_yc);
    for_yc  = (char*)(und->for_yc);
    dom_lam = und->lda_dom;
    for_lam = und->lda_for;
    err     = Fx3DtsPhi(
        evt_time,
        und->sigma_time_rates,
        und->sigma_n_rates,
        und->sigma_dom,
        und->lda_dom,
        &dom_phi);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = Fx3DtsPhi(
        evt_time,
        und->sigma_time_rates,
        und->sigma_n_rates,
        und->sigma_for,
        und->lda_for,
        &for_phi);
    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Extract funding info */
    if (cpd->fund_leg->dom_for == 0)
    {
        fund_yc              = dom_yc;
        fund_lam             = dom_lam;
        fund_phi             = dom_phi;
        eval_const->fund_var = 0;
    }
    else
    {
        fund_yc              = for_yc;
        fund_lam             = for_lam;
        fund_phi             = for_phi;
        eval_const->fund_var = 1;
    }

    /*	Precalculate DFF, gamma, 0.5 * gamma * gamma, and forward vols
            --------------------------------------------------------------	*/

    /*	Spot fx */

    eval_const->fx_spot_date = add_unit(evt_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    eval_const->fx_spot_time = (eval_const->fx_spot_date - today) * YEARS_IN_DAY;

    eval_const->spot_fx_dff_log_ratio =
        log(swp_f_df(evt_date, eval_const->fx_spot_date, for_yc) /
            swp_f_df(evt_date, eval_const->fx_spot_date, dom_yc));

    eval_const->spot_fx_dom_gam =
        (1.0 - exp(-dom_lam * (eval_const->fx_spot_time - evt_time))) / dom_lam;

    eval_const->spot_fx_for_gam =
        (1.0 - exp(-for_lam * (eval_const->fx_spot_time - evt_time))) / for_lam;

    eval_const->spot_fx_cvx =
        -0.5 * eval_const->spot_fx_for_gam * eval_const->spot_fx_for_gam * for_phi +
        0.5 * eval_const->spot_fx_dom_gam * eval_const->spot_fx_dom_gam * dom_phi;

    /*	Discounting funding leg */

    if (call->num_fund_cpn > 0)
    {
        eval_const->do_fund      = 1;
        eval_const->num_fund_cpn = call->num_fund_cpn;

        for (i = call->fund_idx; i < cpd->fund_leg->num_cpn; i++)
        {
            eval_const->fund_dff[i - call->fund_idx] =
                swp_f_df(evt_date, cpd->fund_leg->cpn[i].pay_date, fund_yc);
            eval_const->fund_gam[i - call->fund_idx] =
                (1.0 - exp(-fund_lam * (cpd->fund_leg->cpn[i].pay_time - evt_time))) / fund_lam;
            eval_const->fund_gam2[i - call->fund_idx] =
                0.5 * eval_const->fund_gam[i - call->fund_idx] *
                eval_const->fund_gam[i - call->fund_idx] * fund_phi;
        }
    }
    else
    {
        eval_const->do_fund      = 0;
        eval_const->num_fund_cpn = 0;
        err                      = "No funding coupons left at call date";
        goto FREE_RETURN;
    }

    eval_const->start_dff =
        swp_f_df(evt_date, cpd->fund_leg->cpn[call->fund_idx].start_date, fund_yc);
    eval_const->start_gam =
        (1.0 - exp(-fund_lam * (cpd->fund_leg->cpn[call->fund_idx].start_time - evt_time))) /
        fund_lam;
    eval_const->start_gam2 = 0.5 * eval_const->start_gam * eval_const->start_gam * fund_phi;

    /*	Discounting pd leg */

    if (call->num_pd_cpn > 0)
    {
        eval_const->do_pd_disc = 1;
        eval_const->num_pd_cpn = call->num_pd_cpn;

        for (i = call->pd_idx; i < cpd->pd_leg->num_cpn; i++)
        {
            eval_const->pd_disc_dff[i - call->pd_idx] =
                swp_f_df(evt_date, cpd->pd_leg->cpn[i].pay_date, dom_yc);
            eval_const->pd_disc_gam[i - call->pd_idx] =
                (1.0 - exp(-dom_lam * (cpd->pd_leg->cpn[i].pay_time - evt_time))) / dom_lam;
            eval_const->pd_disc_gam2[i - call->pd_idx] =
                0.5 * eval_const->pd_disc_gam[i - call->pd_idx] *
                eval_const->pd_disc_gam[i - call->pd_idx] * dom_phi;
        }
    }
    else
    {
        eval_const->do_pd_disc = 0;
        eval_const->num_pd_cpn = 0;
        err                    = "No exotic coupons left at call date";
        goto FREE_RETURN;
    }

    /*	Forwards for pd leg: domestic df, foreign df and adjustments,
            also volatilities for pd leg */

    if (call->num_pd_cpn > 0)
    {
        eval_const->do_pd_fwd = 1;

        for (i = call->pd_idx; i < cpd->pd_leg->num_cpn; i++)
        {
            cpn = cpd->pd_leg->cpn + i;

            eval_const->pd_fwd_dff_ratio[i - call->pd_idx] =
                swp_f_df(evt_date, cpn->fx_val_date, for_yc) /
                swp_f_df(evt_date, cpn->fx_val_date, dom_yc);

            eval_const->pd_fwd_dom_gam[i - call->pd_idx] =
                (1.0 - exp(-dom_lam * (cpn->fx_val_time - evt_time))) / dom_lam;

            eval_const->pd_fwd_for_gam[i - call->pd_idx] =
                (1.0 - exp(-for_lam * (cpn->fx_val_time - evt_time))) / for_lam;

            eval_const->pd_fwd_cvx[i - call->pd_idx] =
                -0.5 * eval_const->pd_fwd_for_gam[i - call->pd_idx] *
                    eval_const->pd_fwd_for_gam[i - call->pd_idx] * for_phi +
                0.5 * eval_const->pd_fwd_dom_gam[i - call->pd_idx] *
                    eval_const->pd_fwd_dom_gam[i - call->pd_idx] * dom_phi;

            err = Fx3DtsFwdPayAdjustment_corr(
                evt_time,
                cpn->fx_val_time,
                cpn->fx_val_time,
                cpn->pay_time,
                cpn->fx_fix_time,
                und->sigma_time_rates,
                und->sigma_n_rates,
                und->sigma_dom,
                und->lda_dom,
                und->sigma_for,
                und->lda_for,
                und->sigma_time_fx,
                und->sigma_fx,
                und->sigma_n_fx,
                und->corr_times,
                und->correl_dom_for,
                und->correl_dom_fx,
                und->correl_for_fx,
                und->corr_n_times,
                &temp);

            if (err)
            {
                goto FREE_RETURN;
            }

            eval_const->pd_fwd_cvx[i - call->pd_idx] += temp;

            err = Fx3DtsImpliedVol_corr(
                cpn->fx_val_time,
                evt_time,
                cpn->fx_fix_time,
                und->sigma_time_rates,
                und->sigma_n_rates,
                und->sigma_dom,
                und->lda_dom,
                und->sigma_for,
                und->lda_for,
                und->sigma_time_fx,
                und->sigma_fx,
                und->sigma_n_fx,
                und->corr_times,
                und->correl_dom_for,
                und->correl_dom_fx,
                und->correl_for_fx,
                und->corr_n_times,
                &(eval_const->pd_std[i - call->pd_idx]));

            if (err)
            {
                goto FREE_RETURN;
            }

            eval_const->pd_std[i - call->pd_idx] *= sqrt(cpn->fx_fix_time - evt_time);
            eval_const->pd_half_std[i - call->pd_idx] = 0.5 * eval_const->pd_std[i - call->pd_idx];

            // If option string coupon specification - skip strike calculation:
            if (!cpn->use_opt_str)
            {
                eval_const->pd_abs_beta[i - call->pd_idx] = fabs(cpn->beta);

                if (cpn->floored && fabs(cpn->beta) > 1.0e-16 && cpn->floor > -1.0e-16)
                {
                    eval_const->pd_floor_str[i - call->pd_idx] =
                        (cpn->floor - cpn->alpha) / cpn->beta;

                    if (cpn->beta > 0.0)
                    {
                        if (eval_const->pd_std[i - call->pd_idx] > 1.0e-16 &&
                            eval_const->pd_floor_str[i - call->pd_idx] > 1.0e-16)
                            eval_const->pd_floor_type[i - call->pd_idx] = 4; /*	Put */
                        else
                            eval_const->pd_floor_type[i - call->pd_idx] = 2; /*	Put IV */
                    }
                    else
                    {
                        if (eval_const->pd_std[i - call->pd_idx] > 1.0e-16 &&
                            eval_const->pd_floor_str[i - call->pd_idx] > 1.0e-16)
                            eval_const->pd_floor_type[i - call->pd_idx] = 3; /*	Call */
                        else
                            eval_const->pd_floor_type[i - call->pd_idx] = 1; /*	Call IV */
                    }
                }
                else
                {
                    eval_const->pd_floor_type[i - call->pd_idx] = 0; /*	No floor */
                }

                if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16)
                {
                    eval_const->pd_cap_str[i - call->pd_idx] = (cpn->cap - cpn->alpha) / cpn->beta;

                    if (cpn->beta > 0.0)
                    {
                        if (eval_const->pd_std[i - call->pd_idx] > 1.0e-16 &&
                            eval_const->pd_cap_str[i - call->pd_idx] > 1.0e-16)
                            eval_const->pd_cap_type[i - call->pd_idx] = 3; /*	Call */
                        else
                            eval_const->pd_cap_type[i - call->pd_idx] = 1; /*	Call IV */
                    }
                    else
                    {
                        if (eval_const->pd_std[i - call->pd_idx] > 1.0e-16 &&
                            eval_const->pd_cap_str[i - call->pd_idx] > 1.0e-16)
                            eval_const->pd_cap_type[i - call->pd_idx] = 4; /*	Put */
                        else
                            eval_const->pd_cap_type[i - call->pd_idx] = 2; /*	Put IV */
                    }
                }
                else
                {
                    eval_const->pd_cap_type[i - call->pd_idx] = 0; /*	No cap */
                }
            }  // if (!cpn->use_opt_str)
        }
    }
    else
    {
        eval_const->do_pd_fwd = 0;
        err                   = "No exotic coupons left at call date";
        goto FREE_RETURN;
    }

    /*	Discount to initial notional exchange */

    eval_const->in_not_fund_dff = swp_f_df(evt_date, call->set_date, fund_yc);
    eval_const->in_not_fund_gam = (1.0 - exp(-fund_lam * (call->set_time - evt_time))) / fund_lam;
    eval_const->in_not_fund_gam2 =
        0.5 * eval_const->in_not_fund_gam * eval_const->in_not_fund_gam * fund_phi;

    eval_const->in_not_pd_dff = swp_f_df(evt_date, call->set_date, dom_yc);
    eval_const->in_not_pd_gam = (1.0 - exp(-dom_lam * (call->set_time - evt_time))) / dom_lam;
    eval_const->in_not_pd_gam2 =
        0.5 * eval_const->in_not_pd_gam * eval_const->in_not_pd_gam * dom_phi;

    //	Option string specification case:
    if (call->use_opt_str)
    {
        eval_const->in_not_pd_fwd_dff_ratio = swp_f_df(evt_date, call->fx_val_date, for_yc) /
                                              swp_f_df(evt_date, call->fx_val_date, dom_yc);

        eval_const->in_not_pd_fwd_dom_gam =
            (1.0 - exp(-dom_lam * (call->fx_val_time - evt_time))) / dom_lam;

        eval_const->in_not_pd_fwd_for_gam =
            (1.0 - exp(-for_lam * (call->fx_val_time - evt_time))) / for_lam;

        eval_const->in_not_pd_fwd_cvx =
            -0.5 * eval_const->in_not_pd_fwd_for_gam * eval_const->in_not_pd_fwd_for_gam * for_phi +
            0.5 * eval_const->in_not_pd_fwd_dom_gam * eval_const->in_not_pd_fwd_dom_gam * dom_phi;

        err = Fx3DtsFwdPayAdjustment_corr(
            evt_time,
            call->fx_val_time,
            call->fx_val_time,
            call->set_time,
            call->fx_fix_time,
            und->sigma_time_rates,
            und->sigma_n_rates,
            und->sigma_dom,
            und->lda_dom,
            und->sigma_for,
            und->lda_for,
            und->sigma_time_fx,
            und->sigma_fx,
            und->sigma_n_fx,
            und->corr_times,
            und->correl_dom_for,
            und->correl_dom_fx,
            und->correl_for_fx,
            und->corr_n_times,
            &temp);

        if (err)
        {
            goto FREE_RETURN;
        }

        eval_const->in_not_pd_fwd_cvx += temp;

        err = Fx3DtsImpliedVol_corr(
            call->fx_val_time,
            evt_time,
            call->fx_fix_time,
            und->sigma_time_rates,
            und->sigma_n_rates,
            und->sigma_dom,
            und->lda_dom,
            und->sigma_for,
            und->lda_for,
            und->sigma_time_fx,
            und->sigma_fx,
            und->sigma_n_fx,
            und->corr_times,
            und->correl_dom_for,
            und->correl_dom_fx,
            und->correl_for_fx,
            und->corr_n_times,
            &(eval_const->in_not_pd_std));

        if (err)
        {
            goto FREE_RETURN;
        }

        eval_const->in_not_pd_std *= sqrt(call->fx_fix_time - evt_time);
        eval_const->in_not_pd_half_std = 0.5 * eval_const->in_not_pd_std;
    }

    if (call->pay_rec == 0)
    {
        eval_const->fee = call->fee;
    }
    else
    {
        eval_const->fee = -call->fee;
    }

    /*	Discount to next call date initial notional exchange */
    if (call_idx < cpd->num_calls - 1)
    {
        eval_const->next_fund_dff = swp_f_df(evt_date, next_call->set_date, fund_yc);
        eval_const->next_fund_gam =
            (1.0 - exp(-fund_lam * (next_call->set_time - evt_time))) / fund_lam;
        eval_const->next_fund_gam2 =
            0.5 * eval_const->next_fund_gam * eval_const->next_fund_gam * fund_phi;

        eval_const->next_pd_dff = swp_f_df(evt_date, next_call->set_date, dom_yc);
        eval_const->next_pd_gam =
            (1.0 - exp(-dom_lam * (next_call->set_time - evt_time))) / dom_lam;
        eval_const->next_pd_gam2 =
            0.5 * eval_const->next_pd_gam * eval_const->next_pd_gam * dom_phi;

        temp_date                       = cpd->fund_leg->cpn[next_call->fund_idx].start_date;
        temp_time                       = (temp_date - today) * YEARS_IN_DAY;
        eval_const->next_start_fund_dff = swp_f_df(evt_date, temp_date, fund_yc);
        eval_const->next_start_fund_gam =
            (1.0 - exp(-fund_lam * (temp_time - evt_time))) / fund_lam;
        eval_const->next_start_fund_gam2 =
            0.5 * eval_const->next_start_fund_gam * eval_const->next_start_fund_gam * fund_phi;
    }

    /*	Final notional for pd leg: discount, forward and volatility */

    cpn = &(cpd->pd_leg->not_ref);

    eval_const->fin_not_disc_dff = swp_f_df(evt_date, cpn->pay_date, dom_yc);
    eval_const->fin_not_disc_gam = (1.0 - exp(-dom_lam * (cpn->pay_time - evt_time))) / dom_lam;
    eval_const->fin_not_disc_gam2 =
        0.5 * eval_const->fin_not_disc_gam * eval_const->fin_not_disc_gam * dom_phi;

    eval_const->fin_not_fwd_dff_ratio =
        swp_f_df(evt_date, cpn->fx_val_date, for_yc) / swp_f_df(evt_date, cpn->fx_val_date, dom_yc);

    eval_const->fin_not_fwd_dom_gam =
        (1.0 - exp(-dom_lam * (cpn->fx_val_time - evt_time))) / dom_lam;

    eval_const->fin_not_fwd_for_gam =
        (1.0 - exp(-for_lam * (cpn->fx_val_time - evt_time))) / for_lam;

    eval_const->fin_not_fwd_cvx =
        -0.5 * eval_const->fin_not_fwd_for_gam * eval_const->fin_not_fwd_for_gam * for_phi +
        0.5 * eval_const->fin_not_fwd_dom_gam * eval_const->fin_not_fwd_dom_gam * dom_phi;

    err = Fx3DtsFwdPayAdjustment_corr(
        evt_time,
        cpn->fx_val_time,
        cpn->fx_val_time,
        cpn->pay_time,
        cpn->fx_fix_time,
        und->sigma_time_rates,
        und->sigma_n_rates,
        und->sigma_dom,
        und->lda_dom,
        und->sigma_for,
        und->lda_for,
        und->sigma_time_fx,
        und->sigma_fx,
        und->sigma_n_fx,
        und->corr_times,
        und->correl_dom_for,
        und->correl_dom_fx,
        und->correl_for_fx,
        und->corr_n_times,
        &temp);

    if (err)
    {
        goto FREE_RETURN;
    }

    eval_const->fin_not_fwd_cvx += temp;

    err = Fx3DtsImpliedVol_corr(
        cpn->fx_val_time,
        evt_time,
        cpn->fx_fix_time,
        und->sigma_time_rates,
        und->sigma_n_rates,
        und->sigma_dom,
        und->lda_dom,
        und->sigma_for,
        und->lda_for,
        und->sigma_time_fx,
        und->sigma_fx,
        und->sigma_n_fx,
        und->corr_times,
        und->correl_dom_for,
        und->correl_dom_fx,
        und->correl_for_fx,
        und->corr_n_times,
        &(eval_const->fin_not_std));

    if (err)
    {
        goto FREE_RETURN;
    }

    eval_const->fin_not_std *= sqrt(cpn->fx_fix_time - evt_time);
    eval_const->fin_not_half_std = 0.5 * eval_const->fin_not_std;

    // If option string notional refund specification - skip strike calculation:
    if (!cpn->use_opt_str)
    {
        eval_const->fin_not_abs_beta = fabs(cpn->beta);

        if (cpn->floored && fabs(cpn->beta) > 1.0e-16 && cpn->floor > -1.0e-16)
        {
            eval_const->fin_not_floor_str = (cpn->floor - cpn->alpha) / cpn->beta;
            if (cpn->beta > 0.0)
            {
                if (eval_const->fin_not_std > 1.0e-16 && eval_const->fin_not_floor_str > 1.0e-16)
                    eval_const->fin_not_floor_type = 4; /*	Put */
                else
                    eval_const->fin_not_floor_type = 2; /*	Put IV */
            }
            else
            {
                if (eval_const->fin_not_std > 1.0e-16 && eval_const->fin_not_floor_str > 1.0e-16)
                    eval_const->fin_not_floor_type = 3; /*	Call */
                else
                    eval_const->fin_not_floor_type = 1; /*	Call IV */
            }
        }
        else
        {
            eval_const->fin_not_floor_type = 0; /*	No floor */
        }

        if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16)
        {
            eval_const->fin_not_cap_str = (cpn->cap - cpn->alpha) / cpn->beta;

            if (cpn->beta > 0.0)
            {
                if (eval_const->fin_not_std > 1.0e-16 && eval_const->fin_not_cap_str > 1.0e-16)
                    eval_const->fin_not_cap_type = 3; /*	Call */
                else
                    eval_const->fin_not_cap_type = 1; /*	Call IV */
            }
            else
            {
                if (eval_const->fin_not_std > 1.0e-16 && eval_const->fin_not_cap_str > 1.0e-16)
                    eval_const->fin_not_cap_type = 4; /*	Put */
                else
                    eval_const->fin_not_cap_type = 2; /*	Put IV */
            }
        }
        else
        {
            eval_const->fin_not_cap_type = 0; /*	No cap */
        }
    }

    /*	End of precalculations
            ---------------------- */

FREE_RETURN:

    return err;
}

Err cpd_fill_FxInstrument_dlm(
    CPD_UND        und,
    CPD_STR        cpd,
    PD_CALL        call,
    int            index_cpn,
    PD_EXO_CPN     cpn,
    CPD_EVAL_CONST eval_const,
    int            nb_strikes,
    double*        strikes)
{
    FxBetaDLM_FxOptInst*   Inst;
    FxBetaDLM_FxOptInst*   PrecalcInst = NULL;
    FxBetaDLM_InstPrecalc* InstConst;

    FxBetaDLM_FxOptInst*   PrecalcInstForwardATM = NULL;
    FxBetaDLM_InstPrecalc* InstConstForwardATM   = NULL;

    Err     err = NULL;
    double *used_strikes, *premium, *alpha, *fwd_strike;
    int     i, j;
    double  step_Fwd, step_X, coef_Fwd, coef_X, fwd_premium;
    double  std_X, std_Fwd, std_Fwd_adj;
    double  X0, fwd, vol;
    double  coef, log_coef;

    if (index_cpn >= 0)
    {
        Inst      = &(eval_const->pd_inst[index_cpn]);
        InstConst = &(eval_const->pd_beta_precalc[index_cpn]);
    }
    else
    {
        Inst      = &(eval_const->fin_not_inst);
        InstConst = &(eval_const->fin_not_beta_precalc);
    }

    /* Instrument allocation */
    err = FxBetaDLM_Allocate_FxOptInst(0, nb_strikes, Inst);

    if (err)
        goto FREE_RETURN;

    /* Calculation of the strikes */
    if (und->num_params->iPrecalcSmile)
    {
        used_strikes = &(InstConst->dPrecalcFwdFloor[0]);

        /* First Calculate the STD */
        std_X = sqrt(eval_const->varFFX);

        /* For FFX calculation, we use the equivalent 3 Factor */
        fwd_strike            = (double*)calloc(1, sizeof(double));
        PrecalcInstForwardATM = (FxBetaDLM_FxOptInst*)calloc(1, sizeof(FxBetaDLM_FxOptInst));
        InstConstForwardATM   = (FxBetaDLM_InstPrecalc*)calloc(1, sizeof(FxBetaDLM_InstPrecalc));

        if (!PrecalcInstForwardATM || !InstConstForwardATM || !fwd_strike)
        {
            err = "Memory allocation faillure in cpd_fill_FxInstrument_dlm";
            goto FREE_RETURN;
        }

        fwd_strike[0] = und->model->dCashFx *
                        swp_f_df(und->model->lToday, cpn->fx_val_date, und->model->cYcFor) /
                        swp_f_df(und->model->lToday, cpn->fx_val_date, und->model->cYcDom);

        /* Instrument allocation */
        err = FxBetaDLM_Allocate_FxOptInst(0, 1, PrecalcInstForwardATM);

        if (err)
            goto FREE_RETURN;

        err = FxBetaDLM_Setup_FxOptInst(
            cpn->fx_val_date,
            cpn->pay_date,
            und->model->lToday,
            call->ex_date,
            1,
            fwd_strike,
            NULL,
            SRT_CALL,
            und->model,
            PrecalcInstForwardATM);

        if (err)
            goto FREE_RETURN;

        /* Update the precalculations */
        err = FxBetaDLM_Calculate_AllConst(
            PrecalcInstForwardATM, und->model, und->num_params, InstConstForwardATM);

        if (err)
            goto FREE_RETURN;

        err = FxBetaDLM_Update_InstPrecalc_FromMoment(
            PrecalcInstForwardATM, und->model, und->num_params, InstConstForwardATM);

        if (err)
            goto FREE_RETURN;

        err = FxBetaDLM_Price_FxOptInst(
            PrecalcInstForwardATM,
            und->model,
            InstConstForwardATM,
            und->hermite,
            und->num_params,
            &(fwd_premium));

        err = srt_f_optimpvol(
            fwd_premium,
            fwd_strike[0],
            fwd_strike[0],
            PrecalcInstForwardATM->dVolTime,
            PrecalcInstForwardATM->dDfPayDom,
            SRT_CALL,
            SRT_LOGNORMAL,
            &std_Fwd);

        /*
        err = Fx3DtsImpliedVol_corr(cpn->fx_val_time,
                                                                0.0,
                                                                call->ex_time,
                                                                und->sigma_time_rates,
                                                                und->sigma_n_rates,
                                                                und->sigma_dom,
                                                                und->lda_dom,
                                                                und->sigma_for,
                                                                und->lda_for,
                                                                und->sigma_time_fx,
                                                                und->sigma_fx,
                                                                und->sigma_n_fx,
                                                                und->corr_times,
                                                                und->correl_dom_for,
                                                                und->correl_dom_fx,
                                                                und->correl_for_fx,
                                                                und->corr_n_times,
                                                                &std_Fwd);
        */

        if (err)
            goto FREE_RETURN;

        std_Fwd *= sqrt(call->ex_time);

        if (und->num_params->iUseCxtyAdj)
        {
            std_Fwd_adj = -0.5 * std_Fwd * std_Fwd;
        }
        else
        {
            std_Fwd_adj = 0.0;
        }

        if (und->num_params->iEquiSpacedFwd)
        {
            step_Fwd = 2.0 * und->num_params->dNbStdFwd * std_Fwd / (und->num_params->iNbFwd - 1.0);
            coef_Fwd = -std_Fwd_adj + und->num_params->dNbStdFwd * std_Fwd;
            for (i = 0; i < und->num_params->iNbFwd; i++)
            {
                used_strikes[i] = strikes[0] * exp(coef_Fwd);
                coef_Fwd -= step_Fwd;
            }
        }
        else
        {
            step_Fwd = (norm(-std_Fwd_adj / std_Fwd + und->num_params->dNbStdFwd) -
                        norm(-std_Fwd_adj / std_Fwd - und->num_params->dNbStdFwd)) /
                       (und->num_params->iNbFwd - 1.0);
            coef_Fwd = norm(-std_Fwd_adj / std_Fwd + und->num_params->dNbStdFwd);
            for (i = 0; i < und->num_params->iNbFwd; i++)
            {
                used_strikes[i] = strikes[0] * exp(std_Fwd * inv_cumnorm_fast(coef_Fwd));
                coef_Fwd -= step_Fwd;
            }
        }

        /* Setup the X */
        step_X = 2.0 * und->num_params->dNbStdX * std_X / (und->num_params->iNbX - 1.0);
        coef_X = -und->num_params->dNbStdX * std_X;

        for (i = 0; i < und->num_params->iNbX; i++)
        {
            InstConst->dPrecalcX[i] = coef_X;
            coef_X += step_X;
        }
    }
    else
    {
        used_strikes = &(strikes[0]);
    }

    err = FxBetaDLM_Setup_FxOptInst(
        cpn->fx_val_date,
        cpn->pay_date,
        call->ex_date,
        cpn->fx_fix_date,
        nb_strikes,
        strikes,
        NULL,
        SRT_CALL,
        und->model,
        Inst);

    if (err)
        goto FREE_RETURN;

    /* Remove the discounting */
    Inst->dDfPayDom = 1.0;

    /* Update the precalculations */
    err = FxBetaDLM_Calculate_AllConst(Inst, und->model, und->num_params, InstConst);

    if (err)
        goto FREE_RETURN;

    InstConst->iIsCPD = 1;

    err = FxBetaDLM_Update_InstPrecalc_FromMoment(Inst, und->model, und->num_params, InstConst);

    if (err)
        goto FREE_RETURN;

    for (i = 0; i < und->hermite->iNbPoints; i++)
    {
        InstConst->coupon_hermite[i] = und->hermite->X[i] * InstConst->dQuadStdY;
    }

    /* Do all the precalculations */
    if (und->num_params->iPrecalcSmile)
    {
        /* Create the instrument */
        PrecalcInst = (FxBetaDLM_FxOptInst*)calloc(1, sizeof(FxBetaDLM_FxOptInst));

        if (!PrecalcInst)
        {
            err = "Memory allocation faillure in cpd_fill_FxInstrument_dlm";
            goto FREE_RETURN;
        }

        err = FxBetaDLM_Allocate_FxOptInst(0, und->num_params->iNbFwd, PrecalcInst);

        if (err)
            goto FREE_RETURN;

        err = FxBetaDLM_Setup_FxOptInst(
            cpn->fx_val_date,
            cpn->pay_date,
            call->ex_date,
            cpn->fx_fix_date,
            und->num_params->iNbFwd,
            used_strikes,
            NULL,
            SRT_CALL,
            und->model,
            PrecalcInst);

        if (err)
            goto FREE_RETURN;

        /* Remove the discounting */
        PrecalcInst->dDfPayDom = 1.0;

        /* Transform strike into fwd */
        fwd = PrecalcInst->dFwdFx;

        for (i = 0; i < PrecalcInst->iNbStrike; i++)
        {
            InstConst->dPrecalcFwdFloor[i] = fwd * strikes[0] / PrecalcInst->dStrike[i];
        }

        for (i = 0; i < und->num_params->iNbX; i++)
        {
            X0 = InstConst->dPrecalcX[i];

            err = FxBetaDLM_Update_InstPrecalc_FromX0(
                X0, PrecalcInst, und->model, und->num_params, InstConst);

            if (err)
                goto FREE_RETURN;

            premium = &(InstConst->dPrecalcFloorPriceStd[i][0]);
            alpha   = &(InstConst->dPrecalcFloorPriceStdAlpha[i][0]);

            /* Price the instrument */
            err = FxBetaDLM_Price_FxOptInst(
                PrecalcInst, und->model, InstConst, und->hermite, und->num_params, premium);

            if (err)
                goto FREE_RETURN;

            for (j = 0; j < PrecalcInst->iNbStrike; j++)
            {
                /* Calculate premium on new fwd */
                premium[j] *= strikes[0] / PrecalcInst->dStrike[j];

                if (und->num_params->iPriceVol)
                {
                    /* Transform price into std */
                    err = srt_f_optimpvol(
                        premium[j],
                        InstConst->dPrecalcFwdFloor[j],
                        strikes[0],
                        PrecalcInst->dVolTime,
                        1.0,
                        SRT_CALL,
                        SRT_LOGNORMAL,
                        &vol);

                    if (err)
                    {
                        vol = 1.0E-10;
                        err = NULL;
                    }

                    premium[j] = vol * PrecalcInst->dSqVolTime;
                }

                if (j > 0)
                {
                    alpha[j] = (premium[j] - premium[j - 1]) / (InstConst->dPrecalcFwdFloor[j] -
                                                                InstConst->dPrecalcFwdFloor[j - 1]);
                }
            }
        }

        /* If needed we do the same for the cap */
        if (nb_strikes > 2)
        {
            used_strikes = &(InstConst->dPrecalcFwdCap[0]);
            coef         = strikes[1] / strikes[0];
            log_coef     = log(coef);

            for (i = 0; i < und->num_params->iNbFwd; i++)
            {
                used_strikes[i]         = PrecalcInst->dStrike[i] * coef;
                PrecalcInst->dStrike[i] = used_strikes[i];
                PrecalcInst->dLogStrike[i] += log_coef;
            }

            /* Transform strike into fwd */
            fwd = PrecalcInst->dFwdFx;

            for (i = 0; i < PrecalcInst->iNbStrike; i++)
            {
                InstConst->dPrecalcFwdCap[i] = fwd * strikes[1] / PrecalcInst->dStrike[i];
            }

            for (i = 0; i < und->num_params->iNbX; i++)
            {
                X0 = InstConst->dPrecalcX[i];

                err = FxBetaDLM_Update_InstPrecalc_FromX0(
                    X0, PrecalcInst, und->model, und->num_params, InstConst);

                if (err)
                    goto FREE_RETURN;

                premium = &(InstConst->dPrecalcCapPriceStd[i][0]);
                alpha   = &(InstConst->dPrecalcCapPriceStdAlpha[i][0]);

                /* Price the instrument */
                err = FxBetaDLM_Price_FxOptInst(
                    PrecalcInst, und->model, InstConst, und->hermite, und->num_params, premium);

                if (err)
                    goto FREE_RETURN;

                for (j = 0; j < PrecalcInst->iNbStrike; j++)
                {
                    /* Calculate premium on new fwd */
                    premium[j] *= strikes[1] / PrecalcInst->dStrike[j];

                    if (und->num_params->iPriceVol)
                    {
                        /* Transform price into std */
                        err = srt_f_optimpvol(
                            premium[j],
                            InstConst->dPrecalcFwdCap[j],
                            strikes[1],
                            PrecalcInst->dVolTime,
                            1.0,
                            SRT_CALL,
                            SRT_LOGNORMAL,
                            &vol);

                        if (err)
                        {
                            vol = 1.0E-10;
                            err = NULL;
                        }

                        premium[j] = vol * PrecalcInst->dSqVolTime;
                    }

                    if (j > 0)
                    {
                        alpha[j] =
                            (premium[j] - premium[j - 1]) /
                            (InstConst->dPrecalcFwdCap[j] - InstConst->dPrecalcFwdCap[j - 1]);
                    }
                }
            }
        }
    }

FREE_RETURN:

    if (PrecalcInst)
    {
        FxBetaDLM_Free_FxOptInst(PrecalcInst);
        free(PrecalcInst);
    }

    if (PrecalcInstForwardATM)
    {
        FxBetaDLM_Free_FxOptInst(PrecalcInstForwardATM);
        free(PrecalcInstForwardATM);
    }

    if (InstConstForwardATM)
        free(InstConstForwardATM);
    return err;
}

Err cpd_fill_eval_const_dlm(
    CPD_UND und,
    CPD_STR cpd,
    /*	Index of the current call/KO */
    int            call_idx,
    int            is_mc,
    CPD_EVAL_CONST eval_const)
{
    PD_CALL    call, next_call;
    PD_EXO_CPN cpn;
    long       evt_date;
    double     evt_time;
    char*      fund_yc;
    double     dom_lam, for_lam, dom_phi, for_phi, dom_beta, for_beta;
    double     fund_lam, fund_phi;
    long       today, temp_date;
    double     temp_time;
    double     beta_temp;
    char *     dom_yc, *for_yc;
    int        last_used_pd, last_used_fund;
    int        eval_redemption;
    double     fwd_fx;
    int        nb_strikes;
    double     strikes[3];
    int        i;
    Err        err = NULL;

    /*	Extract info from call */
    call = cpd->call + call_idx;
    if (call_idx < cpd->num_calls - 1)
    {
        next_call = cpd->call + (call_idx + 1);
    }

    evt_date = call->ex_date;
    evt_time = call->ex_time;

    /*	Extract info from und */
    today   = und->today;
    dom_yc  = (char*)(und->dom_yc);
    for_yc  = (char*)(und->for_yc);
    dom_lam = und->lda_dom;
    for_lam = und->lda_for;

    dom_phi  = eval_const->phi_dom;
    for_phi  = eval_const->phi_for;
    dom_beta = eval_const->beta_dom;
    for_beta = eval_const->beta_for;

    /*	Extract funding info */
    if (cpd->fund_leg->dom_for == 0)
    {
        fund_yc              = dom_yc;
        fund_lam             = dom_lam;
        fund_phi             = dom_phi;
        eval_const->fund_var = 0;
    }
    else
    {
        fund_yc              = for_yc;
        fund_lam             = for_lam;
        fund_phi             = for_phi;
        eval_const->fund_var = 1;
    }

    /* Calculate the number of needed coupons */
    if (call_idx == cpd->num_calls - 1)
    {
        last_used_pd    = cpd->pd_leg->num_cpn - 1;
        last_used_fund  = cpd->fund_leg->num_cpn - 1;
        eval_redemption = 1;
    }
    else
    {
        last_used_pd    = (cpd->call + (call_idx + 1))->pd_idx - 1;
        last_used_fund  = (cpd->call + (call_idx + 1))->fund_idx - 1;
        eval_redemption = 0;
    }

    eval_const->nb_pd_inst = last_used_pd - cpd->call->pd_idx + 1;

    if (eval_const->nb_pd_inst > 0)
    {
        eval_const->pd_inst =
            (FxBetaDLM_FxOptInst*)calloc(eval_const->nb_pd_inst, sizeof(FxBetaDLM_FxOptInst));
        eval_const->pd_beta_precalc =
            (FxBetaDLM_InstPrecalc*)calloc(eval_const->nb_pd_inst, sizeof(FxBetaDLM_InstPrecalc));

        if (!eval_const->pd_inst || !eval_const->pd_beta_precalc)
        {
            err = "Memory allocation faillure in cpd_fill_eval_const_dlm";
            goto FREE_RETURN;
        }
    }

    /*	Precalculate DFF, gamma, 0.5 * gamma * gamma, and forward vols
            --------------------------------------------------------------	*/

    /*	Spot fx */

    eval_const->fx_spot_date = add_unit(evt_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    eval_const->fx_spot_time = (eval_const->fx_spot_date - today) * YEARS_IN_DAY;

    eval_const->spot_fx_dff_log_ratio =
        log(swp_f_df(evt_date, eval_const->fx_spot_date, for_yc) /
            swp_f_df(evt_date, eval_const->fx_spot_date, dom_yc));

    if (is_mc)
    {
        beta_temp =
            (1.0 - exp(-dom_lam * (eval_const->fx_spot_time - und->model->dTStar))) / dom_lam;
        eval_const->spot_fx_dom_gam = beta_temp - dom_beta;
        eval_const->spot_fx_cvx     = 0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;

        beta_temp =
            (1.0 - exp(-for_lam * (eval_const->fx_spot_time - und->model->dTStar))) / for_lam;
        eval_const->spot_fx_for_gam = beta_temp - for_beta;
        eval_const->spot_fx_cvx -= 0.5 * (beta_temp * beta_temp - for_beta * for_beta) * for_phi;
    }
    else
    {
        eval_const->spot_fx_dom_gam =
            (1.0 - exp(-dom_lam * (eval_const->fx_spot_time - evt_time))) / dom_lam;

        eval_const->spot_fx_for_gam =
            (1.0 - exp(-for_lam * (eval_const->fx_spot_time - evt_time))) / for_lam;

        eval_const->spot_fx_cvx =
            -0.5 * eval_const->spot_fx_for_gam * eval_const->spot_fx_for_gam * for_phi +
            0.5 * eval_const->spot_fx_dom_gam * eval_const->spot_fx_dom_gam * dom_phi;
    }

    /*	Discounting funding leg */

    if (call->num_fund_cpn > 0)
    {
        eval_const->do_fund      = 1;
        eval_const->num_fund_cpn = call->num_fund_cpn;

        for (i = call->fund_idx; i <= last_used_fund; i++)
        {
            eval_const->fund_dff[i - call->fund_idx] =
                swp_f_df(evt_date, cpd->fund_leg->cpn[i].pay_date, fund_yc);

            if (is_mc)
            {
                beta_temp =
                    (1.0 - exp(-fund_lam * (cpd->fund_leg->cpn[i].pay_time - und->model->dTStar))) /
                    fund_lam;
                eval_const->fund_gam[i - call->fund_idx] = beta_temp - dom_beta;
                eval_const->fund_gam2[i - call->fund_idx] =
                    0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * fund_phi;
            }
            else
            {
                eval_const->fund_gam[i - call->fund_idx] =
                    (1.0 - exp(-fund_lam * (cpd->fund_leg->cpn[i].pay_time - evt_time))) / fund_lam;
                eval_const->fund_gam2[i - call->fund_idx] =
                    0.5 * eval_const->fund_gam[i - call->fund_idx] *
                    eval_const->fund_gam[i - call->fund_idx] * fund_phi;
            }
        }
    }
    else
    {
        eval_const->do_fund      = 0;
        eval_const->num_fund_cpn = 0;
        err                      = "No funding coupons left at call date";
        goto FREE_RETURN;
    }

    eval_const->start_dff =
        swp_f_df(evt_date, cpd->fund_leg->cpn[call->fund_idx].start_date, fund_yc);

    if (is_mc)
    {
        beta_temp = (1.0 - exp(-fund_lam * (cpd->fund_leg->cpn[call->fund_idx].start_time -
                                            und->model->dTStar))) /
                    fund_lam;
        eval_const->start_gam  = beta_temp - dom_beta;
        eval_const->start_gam2 = 0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * fund_phi;
    }
    else
    {
        eval_const->start_gam =
            (1.0 - exp(-fund_lam * (cpd->fund_leg->cpn[call->fund_idx].start_time - evt_time))) /
            fund_lam;
        eval_const->start_gam2 = 0.5 * eval_const->start_gam * eval_const->start_gam * fund_phi;
    }

    /*	Discounting pd leg */

    if (call->num_pd_cpn > 0)
    {
        eval_const->do_pd_disc = 1;
        eval_const->num_pd_cpn = call->num_pd_cpn;

        for (i = call->pd_idx; i <= last_used_pd; i++)
        {
            eval_const->pd_disc_dff[i - call->pd_idx] =
                swp_f_df(evt_date, cpd->pd_leg->cpn[i].pay_date, dom_yc);

            if (is_mc)
            {
                beta_temp =
                    (1.0 - exp(-dom_lam * (cpd->pd_leg->cpn[i].pay_time - und->model->dTStar))) /
                    dom_lam;
                eval_const->pd_disc_gam[i - call->pd_idx] = beta_temp - dom_beta;
                eval_const->pd_disc_gam2[i - call->pd_idx] =
                    0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;
            }
            else
            {
                eval_const->pd_disc_gam[i - call->pd_idx] =
                    (1.0 - exp(-dom_lam * (cpd->pd_leg->cpn[i].pay_time - evt_time))) / dom_lam;
                eval_const->pd_disc_gam2[i - call->pd_idx] =
                    0.5 * eval_const->pd_disc_gam[i - call->pd_idx] *
                    eval_const->pd_disc_gam[i - call->pd_idx] * dom_phi;
            }
        }
    }
    else
    {
        eval_const->do_pd_disc = 0;
        eval_const->num_pd_cpn = 0;
        err                    = "No exotic coupons left at call date";
        goto FREE_RETURN;
    }

    /* X and spot Fx reconstruction */

    fwd_fx =
        und->model->dCashFx * swp_f_df(today, evt_date, for_yc) / swp_f_df(today, evt_date, dom_yc);

    if (is_mc)
    {
        eval_const->pd_S_lin = 1.0 + 2.0 * und->model->dC0 * eval_const->varFFX;

        eval_const->pd_S_const = log(fwd_fx);
        eval_const->pd_S_const +=
            -0.5 * (und->model->dB0 * und->model->dB0 * eval_const->varFFX / eval_const->pd_S_lin +
                    log(eval_const->pd_S_lin));
        eval_const->pd_S_const +=
            0.5 * (dom_beta * dom_beta * dom_phi - for_beta * for_beta * for_phi);

        eval_const->pd_S_quad = und->model->dC0 / eval_const->pd_S_lin;
        eval_const->pd_S_lin  = und->model->dB0 / eval_const->pd_S_lin;

        eval_const->pd_S_dom = dom_beta;
        eval_const->pd_S_for = -for_beta;
    }
    else
    {
        eval_const->pd_S_lin  = 1.0 / (1.0 + 2.0 * und->model->dC0 * eval_const->varFFX);
        eval_const->pd_S_quad = und->model->dC0 * eval_const->pd_S_lin;
        eval_const->pd_S_const =
            -0.5 * (eval_const->varFFX * eval_const->pd_S_lin * und->model->dB0 * und->model->dB0 -
                    log(eval_const->pd_S_lin));
        eval_const->pd_S_lin *= und->model->dB0;

        eval_const->pd_X_const =
            0.5 * (eval_const->varFFX +
                   eval_const->beta_dom * eval_const->beta_dom * eval_const->phi_dom -
                   eval_const->beta_for * eval_const->beta_for * eval_const->phi_for);
        eval_const->pd_X_dom = eval_const->beta_dom;
        eval_const->pd_X_for = -eval_const->beta_for;

        eval_const->pd_S_const += log(fwd_fx);
        eval_const->pd_S_const +=
            -0.5 * (eval_const->beta_dom * eval_const->beta_dom * eval_const->phi_dom -
                    eval_const->beta_for * eval_const->beta_for * eval_const->phi_for);

        eval_const->pd_S_dom = -eval_const->pd_X_dom;
        eval_const->pd_S_for = -eval_const->pd_X_for;
    }

    /*	Forwards for pd leg: domestic df, foreign df and adjustments,
            also volatilities for pd leg */
    if (call->num_pd_cpn > 0)
    {
        eval_const->do_pd_fwd = 1;

        for (i = call->pd_idx; i <= last_used_pd; i++)
        {
            cpn = cpd->pd_leg->cpn + i;

            nb_strikes = 0;

            eval_const->pd_abs_beta[i - call->pd_idx] = fabs(cpn->beta);

            if (cpn->floored && fabs(cpn->beta) > 1.0e-16 && cpn->floor > -1.0e-16)
            {
                eval_const->pd_floor_str[i - call->pd_idx] = (cpn->floor - cpn->alpha) / cpn->beta;

                if (cpn->beta > 0.0)
                {
                    if (eval_const->pd_floor_str[i - call->pd_idx] > 1.0e-16)
                    {
                        eval_const->pd_floor_type[i - call->pd_idx] = 4; /*	Put */
                        strikes[nb_strikes] = eval_const->pd_floor_str[i - call->pd_idx];
                        nb_strikes++;
                    }
                    else
                    {
                        eval_const->pd_floor_type[i - call->pd_idx] = 2; /*	Put IV */
                    }
                }
                else
                {
                    if (eval_const->pd_floor_str[i - call->pd_idx] > 1.0e-16)
                    {
                        eval_const->pd_floor_type[i - call->pd_idx] = 3; /*	Call */
                        strikes[nb_strikes] = eval_const->pd_floor_str[i - call->pd_idx];
                        nb_strikes++;
                    }
                    else
                    {
                        eval_const->pd_floor_type[i - call->pd_idx] = 1; /*	Call IV */
                    }
                }
            }
            else
            {
                eval_const->pd_floor_type[i - call->pd_idx] = 0; /*	No floor */
            }

            if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16)
            {
                eval_const->pd_cap_str[i - call->pd_idx] = (cpn->cap - cpn->alpha) / cpn->beta;

                if (cpn->beta > 0.0)
                {
                    if (eval_const->pd_cap_str[i - call->pd_idx] > 1.0e-16)
                    {
                        eval_const->pd_cap_type[i - call->pd_idx] = 3; /*	Call */
                        strikes[nb_strikes] = eval_const->pd_cap_str[i - call->pd_idx];
                        eval_const->index_cap[i - call->pd_idx] = nb_strikes;
                        nb_strikes++;
                    }
                    else
                    {
                        eval_const->pd_cap_type[i - call->pd_idx] = 1; /*	Call IV */
                    }
                }
                else
                {
                    if (eval_const->pd_cap_str[i - call->pd_idx] > 1.0e-16)
                    {
                        eval_const->pd_cap_type[i - call->pd_idx] = 4; /*	Put */
                        strikes[nb_strikes] = eval_const->pd_cap_str[i - call->pd_idx];
                        eval_const->index_cap[i - call->pd_idx] = nb_strikes;
                        nb_strikes++;
                    }
                    else
                    {
                        eval_const->pd_cap_type[i - call->pd_idx] = 2; /*	Put IV */
                    }
                }
            }
            else
            {
                eval_const->pd_cap_type[i - call->pd_idx] = 0; /*	No cap */
            }

            if (eval_const->pd_abs_beta[i - call->pd_idx] > 1.0E-16 ||
                eval_const->pd_floor_type[i - call->pd_idx] == 1 ||
                eval_const->pd_floor_type[i - call->pd_idx] == 2 ||
                eval_const->pd_floor_type[i - call->pd_idx] == 4 ||
                eval_const->pd_cap_type[i - call->pd_idx] > 0)
            {
                /* calculation of the fwd */
                strikes[nb_strikes]                     = fwd_fx / 10000.0;
                eval_const->index_fwd[i - call->pd_idx] = nb_strikes;
                nb_strikes++;
            }
            else
            {
                eval_const->index_fwd[i - call->pd_idx] = -1;
            }

            /* Setup the instrument */
            if (nb_strikes > 0)
            {
                err = cpd_fill_FxInstrument_dlm(
                    und, cpd, call, i - call->pd_idx, cpn, eval_const, nb_strikes, strikes);
                if (err)
                    goto FREE_RETURN;
            }
            else
            {
                eval_const->pd_inst[i - call->pd_idx].iNbStrike = 0;
            }

            /* Fwd Fx reconstruction */

            eval_const->pd_fwd_dff_ratio[i - call->pd_idx] =
                swp_f_df(evt_date, cpn->fx_val_date, for_yc) /
                swp_f_df(evt_date, cpn->fx_val_date, dom_yc);

            if (is_mc)
            {
                beta_temp =
                    (1.0 - exp(-dom_lam * (cpn->fx_val_time - und->model->dTStar))) / dom_lam;
                eval_const->pd_fwd_dom_gam[i - call->pd_idx] = beta_temp - dom_beta;
                eval_const->pd_fwd_cvx[i - call->pd_idx] =
                    0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;

                beta_temp =
                    (1.0 - exp(-for_lam * (cpn->fx_val_time - und->model->dTStar))) / for_lam;
                eval_const->pd_fwd_for_gam[i - call->pd_idx] = beta_temp - for_beta;
                eval_const->pd_fwd_cvx[i - call->pd_idx] -=
                    0.5 * (beta_temp * beta_temp - for_beta * for_beta) * for_phi;
            }
            else
            {
                eval_const->pd_fwd_dom_gam[i - call->pd_idx] =
                    (1.0 - exp(-dom_lam * (cpn->fx_val_time - evt_time))) / dom_lam;

                eval_const->pd_fwd_for_gam[i - call->pd_idx] =
                    (1.0 - exp(-for_lam * (cpn->fx_val_time - evt_time))) / for_lam;

                eval_const->pd_fwd_cvx[i - call->pd_idx] =
                    -0.5 * eval_const->pd_fwd_for_gam[i - call->pd_idx] *
                        eval_const->pd_fwd_for_gam[i - call->pd_idx] * eval_const->phi_for +
                    0.5 * eval_const->pd_fwd_dom_gam[i - call->pd_idx] *
                        eval_const->pd_fwd_dom_gam[i - call->pd_idx] * eval_const->phi_dom;
            }
        }
    }
    else
    {
        eval_const->do_pd_fwd = 0;
        err                   = "No exotic coupons left at call date";
        goto FREE_RETURN;
    }

    /*	Discount to initial notional exchange */

    eval_const->in_not_fund_dff = swp_f_df(evt_date, call->set_date, fund_yc);

    if (is_mc)
    {
        beta_temp = (1.0 - exp(-fund_lam * (call->set_time - und->model->dTStar))) / fund_lam;
        eval_const->in_not_fund_gam = beta_temp - dom_beta;
        eval_const->in_not_fund_gam2 =
            0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * fund_phi;
    }
    else
    {
        eval_const->in_not_fund_gam =
            (1.0 - exp(-fund_lam * (call->set_time - evt_time))) / fund_lam;
        eval_const->in_not_fund_gam2 =
            0.5 * eval_const->in_not_fund_gam * eval_const->in_not_fund_gam * fund_phi;
    }

    eval_const->in_not_pd_dff = swp_f_df(evt_date, call->set_date, dom_yc);

    if (is_mc)
    {
        beta_temp = (1.0 - exp(-dom_lam * (call->set_time - und->model->dTStar))) / dom_lam;
        eval_const->in_not_pd_gam  = beta_temp - dom_beta;
        eval_const->in_not_pd_gam2 = 0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;
    }
    else
    {
        eval_const->in_not_pd_gam = (1.0 - exp(-dom_lam * (call->set_time - evt_time))) / dom_lam;
        eval_const->in_not_pd_gam2 =
            0.5 * eval_const->in_not_pd_gam * eval_const->in_not_pd_gam * dom_phi;
    }

    if (call->pay_rec == 0)
    {
        eval_const->fee = call->fee;
    }
    else
    {
        eval_const->fee = -call->fee;
    }

    /*	Discount to next call date initial notional exchange */
    if (!eval_redemption)
    {
        eval_const->next_fund_dff = swp_f_df(evt_date, next_call->set_date, fund_yc);

        if (is_mc)
        {
            beta_temp =
                (1.0 - exp(-fund_lam * (next_call->set_time - und->model->dTStar))) / fund_lam;
            eval_const->next_fund_gam = beta_temp - dom_beta;
            eval_const->next_fund_gam2 =
                0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * fund_phi;
        }
        else
        {
            eval_const->next_fund_gam =
                (1.0 - exp(-fund_lam * (next_call->set_time - evt_time))) / fund_lam;
            eval_const->next_fund_gam2 =
                0.5 * eval_const->next_fund_gam * eval_const->next_fund_gam * fund_phi;
        }

        eval_const->next_pd_dff = swp_f_df(evt_date, next_call->set_date, dom_yc);

        if (is_mc)
        {
            beta_temp =
                (1.0 - exp(-dom_lam * (next_call->set_time - und->model->dTStar))) / dom_lam;
            eval_const->next_pd_gam = beta_temp - dom_beta;
            eval_const->next_pd_gam2 =
                0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;
        }
        else
        {
            eval_const->next_pd_gam =
                (1.0 - exp(-dom_lam * (next_call->set_time - evt_time))) / dom_lam;
            eval_const->next_pd_gam2 =
                0.5 * eval_const->next_pd_gam * eval_const->next_pd_gam * dom_phi;
        }

        temp_date                       = cpd->fund_leg->cpn[next_call->fund_idx].start_date;
        temp_time                       = (temp_date - today) * YEARS_IN_DAY;
        eval_const->next_start_fund_dff = swp_f_df(evt_date, temp_date, fund_yc);

        if (is_mc)
        {
            beta_temp = (1.0 - exp(-fund_lam * (temp_time - und->model->dTStar))) / fund_lam;
            eval_const->next_start_fund_gam = beta_temp - dom_beta;
            eval_const->next_start_fund_gam2 =
                0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * fund_phi;
        }
        else
        {
            eval_const->next_start_fund_gam =
                (1.0 - exp(-fund_lam * (temp_time - evt_time))) / fund_lam;
            eval_const->next_start_fund_gam2 =
                0.5 * eval_const->next_start_fund_gam * eval_const->next_start_fund_gam * fund_phi;
        }
    }
    else
    {
        /*	Final notional for pd leg: discount, forward and volatility */

        cpn = &(cpd->pd_leg->not_ref);

        eval_const->fin_not_disc_dff = swp_f_df(evt_date, cpn->pay_date, dom_yc);

        if (is_mc)
        {
            beta_temp = (1.0 - exp(-dom_lam * (cpn->pay_time - und->model->dTStar))) / dom_lam;
            eval_const->fin_not_disc_gam = beta_temp - dom_beta;
            eval_const->fin_not_disc_gam2 =
                0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;
        }
        else
        {
            eval_const->fin_not_disc_gam =
                (1.0 - exp(-dom_lam * (cpn->pay_time - evt_time))) / dom_lam;
            eval_const->fin_not_disc_gam2 =
                0.5 * eval_const->fin_not_disc_gam * eval_const->fin_not_disc_gam * dom_phi;
        }

        eval_const->fin_not_fwd_dff_ratio = swp_f_df(evt_date, cpn->fx_val_date, for_yc) /
                                            swp_f_df(evt_date, cpn->fx_val_date, dom_yc);

        if (is_mc)
        {
            beta_temp = (1.0 - exp(-dom_lam * (cpn->fx_val_time - und->model->dTStar))) / dom_lam;
            eval_const->fin_not_fwd_dom_gam = beta_temp - dom_beta;
            eval_const->fin_not_fwd_cvx =
                0.5 * (beta_temp * beta_temp - dom_beta * dom_beta) * dom_phi;

            beta_temp = (1.0 - exp(-for_lam * (cpn->fx_val_time - und->model->dTStar))) / for_lam;
            eval_const->fin_not_fwd_for_gam = beta_temp - for_beta;
            eval_const->fin_not_fwd_cvx -=
                0.5 * (beta_temp * beta_temp - for_beta * for_beta) * for_phi;
        }
        else
        {
            eval_const->fin_not_fwd_dom_gam =
                (1.0 - exp(-dom_lam * (cpn->fx_val_time - evt_time))) / dom_lam;

            eval_const->fin_not_fwd_for_gam =
                (1.0 - exp(-for_lam * (cpn->fx_val_time - evt_time))) / for_lam;

            eval_const->fin_not_fwd_cvx =
                -0.5 * eval_const->fin_not_fwd_for_gam * eval_const->fin_not_fwd_for_gam * for_phi +
                0.5 * eval_const->fin_not_fwd_dom_gam * eval_const->fin_not_fwd_dom_gam * dom_phi;
        }

        eval_const->fin_not_abs_beta = fabs(cpn->beta);

        nb_strikes = 0;

        if (cpn->floored && fabs(cpn->beta) > 1.0e-16 && cpn->floor > -1.0e-16)
        {
            eval_const->fin_not_floor_str = (cpn->floor - cpn->alpha) / cpn->beta;

            if (cpn->beta > 0.0)
            {
                if (eval_const->fin_not_floor_str > 1.0e-16)
                {
                    eval_const->fin_not_floor_type = 4; /*	Put */
                    strikes[nb_strikes]            = eval_const->fin_not_floor_str;
                    nb_strikes++;
                }
                else
                {
                    eval_const->fin_not_floor_type = 2; /*	Put IV */
                }
            }
            else
            {
                if (eval_const->fin_not_floor_str > 1.0e-16)
                {
                    eval_const->fin_not_floor_type = 3; /*	Call */
                    strikes[nb_strikes]            = eval_const->fin_not_floor_str;
                    nb_strikes++;
                }
                else
                {
                    eval_const->fin_not_floor_type = 1; /*	Call IV */
                }
            }
        }
        else
        {
            eval_const->fin_not_floor_type = 0; /*	No floor */
        }

        if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16)
        {
            eval_const->fin_not_cap_str = (cpn->cap - cpn->alpha) / cpn->beta;

            if (cpn->beta > 0.0)
            {
                if (eval_const->fin_not_cap_str > 1.0e-16)
                {
                    eval_const->fin_not_cap_type  = 3; /*	Call */
                    strikes[nb_strikes]           = eval_const->fin_not_cap_str;
                    eval_const->fin_not_index_cap = nb_strikes;
                    nb_strikes++;
                }
                else
                {
                    eval_const->fin_not_cap_type = 1; /*	Call IV */
                }
            }
            else
            {
                if (eval_const->fin_not_cap_str > 1.0e-16)
                {
                    eval_const->fin_not_cap_type  = 4; /*	Put */
                    strikes[nb_strikes]           = eval_const->fin_not_cap_str;
                    eval_const->fin_not_index_cap = nb_strikes;
                    nb_strikes++;
                }
                else
                {
                    eval_const->fin_not_cap_type = 2; /*	Put IV */
                }
            }
        }
        else
        {
            eval_const->fin_not_cap_type = 0; /*	No cap */
        }

        if (eval_const->fin_not_abs_beta > 1.0E-16 || eval_const->fin_not_floor_type == 1 ||
            eval_const->fin_not_floor_type == 2 || eval_const->fin_not_floor_type == 4 ||
            eval_const->fin_not_cap_type > 0)
        {
            /* calculation of the fwd */
            strikes[nb_strikes]           = fwd_fx / 10000.0;
            eval_const->fin_not_index_fwd = nb_strikes;
            nb_strikes++;
        }
        else
        {
            eval_const->fin_not_index_fwd = -1;
        }

        /* Setup the instrument */
        if (nb_strikes > 0)
        {
            err =
                cpd_fill_FxInstrument_dlm(und, cpd, call, -1, cpn, eval_const, nb_strikes, strikes);
            if (err)
                goto FREE_RETURN;
        }
        else
        {
            eval_const->fin_not_inst.iNbStrike = 0;
        }
    }

    /*	End of precalculations
            ---------------------- */

FREE_RETURN:

    return err;
}

void cpd_free_eval_const(CPD_UND und, CPD_STR cpd, CPD_EVAL_CONST eval_const)
{
    int i;

    if (und->use_beta_dlm)
    {
        if (eval_const->pd_inst)
        {
            for (i = 0; i < eval_const->nb_pd_inst; i++)
            {
                FxBetaDLM_Free_FxOptInst(&(eval_const->pd_inst[i]));
            }

            free(eval_const->pd_inst);
        }

        if (eval_const->pd_beta_precalc)
        {
            free(eval_const->pd_beta_precalc);
        }
    }
}

/*	Function for the arguments of the tree function */

static void static_find_sig(
    CPD_UND und,
    double  time,
    double* sig_dom,
    double* sig_for,
    double* sig_fx,
    double* corr_dom_for,
    double* corr_dom_fx,
    double* corr_for_fx)
{
    int i = 0, j = 0, k = 0;

    while (i < und->sigma_n_rates - 1 && und->sigma_time_rates[i] < time - 1.0e-08)
    {
        i++;
    }

    while (j < und->sigma_n_fx - 1 && und->sigma_time_fx[j] < time - 1.0e-08)
    {
        j++;
    }

    while (k < und->corr_n_times - 1 && und->corr_times[k] < time - 1.0e-08)
    {
        k++;
    }

    *sig_dom      = und->sigma_dom[i];
    *sig_for      = und->sigma_for[i];
    *sig_fx       = und->sigma_fx[j];
    *corr_dom_for = und->correl_dom_for[k];
    *corr_dom_fx  = und->correl_dom_fx[k];
    *corr_for_fx  = und->correl_for_fx[k];
}

Err cpd_fill_tree_arg(
    CPD_UND und,
    CPD_STR cpd,
    /*	Required number of time steps */
    long         req_stp,
    CPD_TREE_ARG tree_arg)
{
    Err         err = NULL;
    CPD_PAY_ARG cpd_prm;
    double      temp, next_bar;
    int         i, j;

    double *dom_phi = NULL, *for_phi = NULL, *dom_beta = NULL, *for_beta = NULL, *var_ffx = NULL;

    /*	Initialise */
    tree_arg->time              = NULL;
    tree_arg->date              = NULL;
    tree_arg->vol_change        = NULL;
    tree_arg->sig_dom           = NULL;
    tree_arg->sig_for           = NULL;
    tree_arg->sig_fx            = NULL;
    tree_arg->corr_dom_for      = NULL;
    tree_arg->corr_dom_fx       = NULL;
    tree_arg->corr_for_fx       = NULL;
    tree_arg->dom_ifr           = NULL;
    tree_arg->dom_fwd           = NULL;
    tree_arg->dom_var           = NULL;
    tree_arg->for_ifr           = NULL;
    tree_arg->for_fwd           = NULL;
    tree_arg->for_var           = NULL;
    tree_arg->fx_fwd            = NULL;
    tree_arg->fx_var            = NULL;
    tree_arg->void_prm          = NULL;
    tree_arg->is_event          = NULL;
    tree_arg->bar_lvl           = NULL;
    tree_arg->bar_cl            = NULL;
    tree_arg->is_bar            = NULL;
    tree_arg->dr_const_dom      = NULL;
    tree_arg->dr_coef_dom       = NULL;
    tree_arg->dr_const_for      = NULL;
    tree_arg->dr_coef1_for      = NULL;
    tree_arg->dr_coef2_for      = NULL;
    tree_arg->dr_coef3_for      = NULL;
    tree_arg->dr_const_fx       = NULL;
    tree_arg->dr_coef1_fx       = NULL;
    tree_arg->dr_coef2_fx       = NULL;
    tree_arg->dr_coef3_fx       = NULL;
    tree_arg->glob_corr_dom_for = NULL;
    tree_arg->glob_corr_dom_fx  = NULL;
    tree_arg->glob_corr_for_fx  = NULL;

    /*	Compute time steps */

    /*	Copy event dates */

    tree_arg->nstp = cpd->num_calls;
    tree_arg->time = (double*)calloc(tree_arg->nstp, sizeof(double));
    if (!tree_arg->time)
    {
        err = "Memory allocation error (1) in cpd_fill_tree_arg";
        goto FREE_RETURN;
    }
    for (i = 0; i < tree_arg->nstp; i++)
    {
        tree_arg->time[i] = cpd->call[i].ex_time;
    }

    /*	Fill time vector */

    /*	Add today if required */
    if (tree_arg->time[0] < -EPS)
    {
        err = "Past event date in cpd_fill_tree_arg";
        goto FREE_RETURN;
    }
    if (tree_arg->time[0] > EPS)
    {
        num_f_add_number((int*)&(tree_arg->nstp), &(tree_arg->time), 0.0);
        num_f_sort_vector(tree_arg->nstp, tree_arg->time);
        num_f_unique_vector((int*)&(tree_arg->nstp), tree_arg->time);
    }

    /*	If only one event today, add empty event */
    if (tree_arg->nstp == 1)
    {
        num_f_add_number((int*)&(tree_arg->nstp), (&tree_arg->time), 1.0);
    }

    /*	Fill the vector */
    num_f_fill_vector_newalgo((int*)&(tree_arg->nstp), &(tree_arg->time), req_stp);

    /*	Make dates */
    tree_arg->date = (double*)calloc(tree_arg->nstp, sizeof(double));
    if (!tree_arg->date)
    {
        err = "Memory allocation error (2) in cpd_fill_tree_arg";
        goto FREE_RETURN;
    }

    for (i = 0; i < tree_arg->nstp; i++)
    {
        tree_arg->date[i] = und->today + DAYS_IN_YEAR * tree_arg->time[i];

        if (i > 0 && tree_arg->date[i] - tree_arg->date[i - 1] >= 1)
        {
            tree_arg->date[i] = (long)(tree_arg->date[i] + 1.0e-08);
            tree_arg->time[i] = YEARS_IN_DAY * (tree_arg->date[i] - und->today);
        }
    }

    /*	Lambdas and correls */

    tree_arg->dom_lam = und->lda_dom;
    tree_arg->for_lam = und->lda_for;

    /*	Spot fx and yield curves */

    tree_arg->spot_fx = und->spot_fx;
    strcpy(tree_arg->dom_yc, und->dom_yc);
    strcpy(tree_arg->for_yc, und->for_yc);

    /*	Fill time info */
    /*	Sigmas */

    if (und->use_beta_dlm)
    {
        tree_arg->vol_change = (int*)calloc(tree_arg->nstp, sizeof(int));
        tree_arg->sig_dom    = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->sig_for    = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->sig_fx     = (double*)calloc(tree_arg->nstp, sizeof(double));

        tree_arg->corr_dom_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->corr_dom_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->corr_for_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));

        tree_arg->dr_const_dom = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef_dom  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_const_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef1_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef2_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef3_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_const_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef1_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef2_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dr_coef3_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));

        /*	Get distributions */
        tree_arg->dom_ifr = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dom_fwd = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dom_var = (double*)calloc(tree_arg->nstp, sizeof(double));
        dom_phi           = (double*)calloc(tree_arg->nstp, sizeof(double));
        dom_beta          = (double*)calloc(tree_arg->nstp, sizeof(double));

        tree_arg->for_ifr = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->for_fwd = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->for_var = (double*)calloc(tree_arg->nstp, sizeof(double));
        for_phi           = (double*)calloc(tree_arg->nstp, sizeof(double));
        for_beta          = (double*)calloc(tree_arg->nstp, sizeof(double));

        tree_arg->fx_fwd = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->fx_var = (double*)calloc(tree_arg->nstp, sizeof(double));
        var_ffx          = (double*)calloc(tree_arg->nstp, sizeof(double));

        tree_arg->glob_corr_dom_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->glob_corr_dom_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->glob_corr_for_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));

        if (!tree_arg->vol_change || !tree_arg->sig_dom || !tree_arg->sig_for ||
            !tree_arg->sig_fx || !tree_arg->corr_dom_for || !tree_arg->corr_dom_fx ||
            !tree_arg->corr_for_fx || !tree_arg->dr_const_dom || !tree_arg->dr_coef_dom ||
            !tree_arg->dr_const_for || !tree_arg->dr_coef1_for || !tree_arg->dr_coef2_for ||
            !tree_arg->dr_coef3_for || !tree_arg->dr_const_fx || !tree_arg->dr_coef1_fx ||
            !tree_arg->dr_coef2_fx || !tree_arg->dr_coef3_fx || !tree_arg->dom_ifr ||
            !tree_arg->dom_fwd || !tree_arg->dom_var || !dom_phi || !dom_beta ||
            !tree_arg->for_ifr || !tree_arg->for_fwd || !tree_arg->for_var || !for_phi ||
            !for_beta || !tree_arg->fx_fwd || !tree_arg->fx_var || !var_ffx ||
            !tree_arg->glob_corr_dom_for || !tree_arg->glob_corr_dom_fx ||
            !tree_arg->glob_corr_for_fx)
        {
            err = "Memory allocation error (5) in SrtGrfn3DFXTree";
            goto FREE_RETURN;
        }

        /* Fill the glob expectations and variances */
        Fx3DBetaDLM_PrecalcGRFNTreeQBeta(
            tree_arg->nstp,
            tree_arg->time,
            und->model,
            tree_arg->dom_fwd,
            tree_arg->dom_var,
            dom_phi,
            dom_beta,
            tree_arg->for_fwd,
            tree_arg->for_var,
            for_phi,
            for_beta,
            tree_arg->fx_fwd,
            tree_arg->fx_var,
            var_ffx,
            tree_arg->glob_corr_dom_for,
            tree_arg->glob_corr_dom_fx,
            tree_arg->glob_corr_for_fx,
            und->grfn_params->min_time);

        /* Fill the local drifts, vols and correlations */
        Fx3DBetaDLM_DriftVolAndCorr(
            tree_arg->nstp,
            tree_arg->time,
            und->model,
            dom_phi,
            dom_beta,
            for_phi,
            for_beta,
            var_ffx,

            tree_arg->dr_const_dom,
            tree_arg->dr_coef_dom,
            tree_arg->dr_const_for,
            tree_arg->dr_coef1_for,
            tree_arg->dr_coef2_for,
            tree_arg->dr_coef3_for,
            tree_arg->dr_const_fx,
            tree_arg->dr_coef1_fx,
            tree_arg->dr_coef2_fx,
            tree_arg->dr_coef3_fx,

            tree_arg->vol_change,
            tree_arg->sig_dom,
            tree_arg->sig_for,
            tree_arg->sig_fx,
            tree_arg->corr_dom_for,
            tree_arg->corr_dom_fx,
            tree_arg->corr_for_fx);

        /* fill the ifr */
        for (i = 0; i < tree_arg->nstp; i++)
        {
            if (i < tree_arg->nstp - 1)
            {
                tree_arg->dom_ifr[i] =
                    swp_f_zr(tree_arg->date[i], tree_arg->date[i + 1], tree_arg->dom_yc);
            }
            else
            {
                tree_arg->dom_ifr[i] =
                    swp_f_zr(tree_arg->date[i], tree_arg->date[i] + 1, tree_arg->dom_yc);
            }
        }
    }
    else
    {
        tree_arg->vol_change = (int*)calloc(tree_arg->nstp, sizeof(int));
        tree_arg->sig_dom    = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->sig_for    = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->sig_fx     = (double*)calloc(tree_arg->nstp, sizeof(double));

        tree_arg->corr_dom_for = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->corr_dom_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->corr_for_fx  = (double*)calloc(tree_arg->nstp, sizeof(double));

        if (!tree_arg->vol_change || !tree_arg->sig_dom || !tree_arg->sig_for ||
            !tree_arg->sig_fx || !tree_arg->corr_dom_for || !tree_arg->corr_dom_fx ||
            !tree_arg->corr_for_fx)
        {
            err = "Memory allocation error (3) in cpd_fill_tree_arg";
            goto FREE_RETURN;
        }

        tree_arg->vol_change[tree_arg->nstp - 1] = 1;

        static_find_sig(
            und,
            tree_arg->time[tree_arg->nstp - 1],
            &(tree_arg->sig_dom[tree_arg->nstp - 1]),
            &(tree_arg->sig_for[tree_arg->nstp - 1]),
            &(tree_arg->sig_fx[tree_arg->nstp - 1]),
            &(tree_arg->corr_dom_for[tree_arg->nstp - 1]),
            &(tree_arg->corr_dom_fx[tree_arg->nstp - 1]),
            &(tree_arg->corr_for_fx[tree_arg->nstp - 1]));

        for (i = tree_arg->nstp - 2; i >= 0; i--)
        {
            static_find_sig(
                und,
                tree_arg->time[i],
                &(tree_arg->sig_dom[i]),
                &(tree_arg->sig_for[i]),
                &(tree_arg->sig_fx[i]),
                &(tree_arg->corr_dom_for[i]),
                &(tree_arg->corr_dom_fx[i]),
                &(tree_arg->corr_for_fx[i]));

            if (fabs(tree_arg->sig_dom[i] - tree_arg->sig_dom[i + 1]) +
                    fabs(tree_arg->sig_for[i] - tree_arg->sig_for[i + 1]) +
                    fabs(tree_arg->sig_fx[i] - tree_arg->sig_fx[i + 1]) +
                    fabs(tree_arg->corr_dom_for[i] - tree_arg->corr_dom_for[i + 1]) +
                    fabs(tree_arg->corr_dom_fx[i] - tree_arg->corr_dom_fx[i + 1]) +
                    fabs(tree_arg->corr_for_fx[i] - tree_arg->corr_for_fx[i + 1]) >
                EPS)
            {
                tree_arg->vol_change[i] = 1;
            }
            else
            {
                tree_arg->vol_change[i] = 0;
            }
        }

        /*	Fill distributions */

        tree_arg->dom_ifr = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dom_fwd = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->dom_var = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->for_ifr = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->for_fwd = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->for_var = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->fx_fwd  = (double*)calloc(tree_arg->nstp, sizeof(double));
        tree_arg->fx_var  = (double*)calloc(tree_arg->nstp, sizeof(double));

        if (!tree_arg->dom_ifr || !tree_arg->dom_fwd || !tree_arg->dom_var || !tree_arg->for_ifr ||
            !tree_arg->for_fwd || !tree_arg->for_var || !tree_arg->fx_fwd || !tree_arg->fx_var)
        {
            err = "Memory allocation error (4) in cpd_fill_tree_arg";
            goto FREE_RETURN;
        }

        fill_fwd_var_corr(
            tree_arg->nstp,
            tree_arg->time,
            tree_arg->date,
            tree_arg->sig_dom,
            tree_arg->sig_for,
            tree_arg->sig_fx,
            tree_arg->dom_lam,
            tree_arg->for_lam,
            tree_arg->corr_dom_for,
            tree_arg->corr_dom_fx,
            tree_arg->corr_for_fx,
            tree_arg->dom_yc,
            tree_arg->for_yc,
            tree_arg->dom_ifr,
            tree_arg->dom_fwd,
            tree_arg->dom_var,
            tree_arg->for_ifr,
            tree_arg->for_fwd,
            tree_arg->for_var,
            tree_arg->fx_fwd,
            tree_arg->fx_var);

        /*	Overwrite beta-forward fx by a more accurate calculation */
        for (i = 0; i < tree_arg->nstp; i++)
        {
            tree_arg->fx_fwd[i] = log(swp_f_df(und->today, tree_arg->date[i], und->for_yc) /
                                      swp_f_df(und->today, tree_arg->date[i], und->dom_yc)) -
                                  0.5 * tree_arg->fx_var[i];

            err = Fx3DtsFwdBetaAdjustment_corr(
                0.0,
                tree_arg->time[i],
                tree_arg->time[i],
                tree_arg->time[i],
                und->sigma_time_rates,
                und->sigma_n_rates,
                und->sigma_dom,
                und->lda_dom,
                und->sigma_for,
                und->lda_for,
                und->sigma_date_fx,
                und->sigma_fx,
                und->sigma_n_fx,
                und->corr_times,
                und->correl_dom_for,
                und->correl_dom_fx,
                und->correl_for_fx,
                und->corr_n_times,
                &temp);

            if (err)
            {
                goto FREE_RETURN;
            }

            tree_arg->fx_fwd[i] += temp;
        }
    }

    /*	Fill limit conditions (product) */

    tree_arg->is_event = (int*)calloc(tree_arg->nstp, sizeof(int));
    tree_arg->void_prm = (void**)calloc(tree_arg->nstp, sizeof(void*));

    tree_arg->is_bar  = (int*)calloc(tree_arg->nstp, sizeof(int));
    tree_arg->bar_lvl = (double*)calloc(tree_arg->nstp, sizeof(double));
    tree_arg->bar_cl  = (int*)calloc(tree_arg->nstp, sizeof(int));

    if (!tree_arg->is_event || !tree_arg->void_prm || !tree_arg->is_bar || !tree_arg->bar_lvl ||
        !tree_arg->bar_cl)
    {
        err = "Memory allocation error (5) in cpd_fill_tree_arg";
        goto FREE_RETURN;
    }

    j        = cpd->num_calls - 1;
    next_bar = 0.0;

    for (i = tree_arg->nstp - 1; i >= 0; i--)
    {
        if (j >= 0 && fabs(tree_arg->date[i] - cpd->call[j].ex_date) < 1.0e-04)
        {
            cpd_prm                     = (cpd_pay_arg*)calloc(1, sizeof(cpd_pay_arg));
            cpd_prm->und                = und;
            cpd_prm->cpd                = cpd;
            cpd_prm->eval_const.pd_inst = NULL;

            cpd_prm->call_idx = j;

            if (und->use_beta_dlm)
            {
                cpd_prm->eval_const.beta_dom = dom_beta[i];
                cpd_prm->eval_const.beta_for = for_beta[i];
                cpd_prm->eval_const.phi_dom  = dom_phi[i];
                cpd_prm->eval_const.phi_for  = for_phi[i];
                cpd_prm->eval_const.varFFX   = var_ffx[i];

                err = cpd_fill_eval_const_dlm(und, cpd, j, 0, &(cpd_prm->eval_const));
            }
            else
            {
                err = cpd_fill_eval_const(und, cpd, j, &(cpd_prm->eval_const));
            }

            if (err)
                goto FREE_RETURN;

            tree_arg->is_event[i] = 1;

            if (cpd->call[j].call_type == 1)
            {
                /*	Knock Out in the Tree	*/
                next_bar = cpd->call[j].barrier *
                           swp_f_df(
                               cpd->call[j].ex_date,
                               add_unit(cpd->call[j].ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                               tree_arg->dom_yc)

                           / swp_f_df(
                                 cpd->call[j].ex_date,
                                 add_unit(cpd->call[j].ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                                 tree_arg->for_yc);

                next_bar             = log(next_bar / tree_arg->spot_fx);
                cpd->call[j].barrier = next_bar;

                tree_arg->is_bar[i] = 1;
                tree_arg->bar_cl[i] = 0;
            }
            else
            {
                tree_arg->is_bar[i] = 0;
                tree_arg->bar_cl[i] = -99999;
            }

            tree_arg->void_prm[i] = (void*)cpd_prm;

            j--;
        }
        else
        {
            tree_arg->is_event[i] = 0;
            tree_arg->void_prm[i] = NULL;
            tree_arg->is_bar[i]   = 0;
            tree_arg->bar_cl[i]   = -99999;
        }

        if (fabs(next_bar) > 1.0e-08)
        {
            tree_arg->bar_lvl[i] = next_bar;
        }
        else
        {
            tree_arg->bar_lvl[i] = -1.0e08;
        }

        if ((i > 0) && (fabs(next_bar) > 1.0e-08))
        {
            next_bar += (tree_arg->for_ifr[i - 1] + tree_arg->for_fwd[i - 1] -
                         tree_arg->dom_ifr[i - 1] - tree_arg->dom_fwd[i - 1]) *
                        (tree_arg->time[i] - tree_arg->time[i - 1]);
        }
    }

FREE_RETURN:

    if (err)
    {
        cpd_free_tree_arg(tree_arg);
    }

    if (dom_phi)
        free(dom_phi);
    if (for_phi)
        free(for_phi);
    if (dom_beta)
        free(dom_beta);
    if (for_beta)
        free(for_beta);
    if (var_ffx)
        free(var_ffx);

    return err;
}

Err cpd_free_tree_arg(CPD_TREE_ARG tree_arg)
{
    int         i;
    CPD_PAY_ARG cpd_prm;

    if (tree_arg->time)
        free(tree_arg->time);
    if (tree_arg->date)
        free(tree_arg->date);
    if (tree_arg->vol_change)
        free(tree_arg->vol_change);
    if (tree_arg->sig_dom)
        free(tree_arg->sig_dom);
    if (tree_arg->sig_for)
        free(tree_arg->sig_for);
    if (tree_arg->sig_fx)
        free(tree_arg->sig_fx);
    if (tree_arg->corr_dom_for)
        free(tree_arg->corr_dom_for);
    if (tree_arg->corr_dom_fx)
        free(tree_arg->corr_dom_fx);
    if (tree_arg->corr_for_fx)
        free(tree_arg->corr_for_fx);
    if (tree_arg->dom_ifr)
        free(tree_arg->dom_ifr);
    if (tree_arg->dom_fwd)
        free(tree_arg->dom_fwd);
    if (tree_arg->dom_var)
        free(tree_arg->dom_var);
    if (tree_arg->for_ifr)
        free(tree_arg->for_ifr);
    if (tree_arg->for_fwd)
        free(tree_arg->for_fwd);
    if (tree_arg->for_var)
        free(tree_arg->for_var);
    if (tree_arg->fx_fwd)
        free(tree_arg->fx_fwd);
    if (tree_arg->fx_var)
        free(tree_arg->fx_var);
    if (tree_arg->is_event)
        free(tree_arg->is_event);
    if (tree_arg->bar_lvl)
        free(tree_arg->bar_lvl);
    if (tree_arg->bar_cl)
        free(tree_arg->bar_cl);
    if (tree_arg->is_bar)
        free(tree_arg->is_bar);

    if (tree_arg->dr_const_dom)
        free(tree_arg->dr_const_dom);
    if (tree_arg->dr_coef_dom)
        free(tree_arg->dr_coef_dom);
    if (tree_arg->dr_const_for)
        free(tree_arg->dr_const_for);
    if (tree_arg->dr_coef1_for)
        free(tree_arg->dr_coef1_for);
    if (tree_arg->dr_coef2_for)
        free(tree_arg->dr_coef2_for);
    if (tree_arg->dr_coef3_for)
        free(tree_arg->dr_coef3_for);
    if (tree_arg->dr_const_fx)
        free(tree_arg->dr_const_fx);
    if (tree_arg->dr_coef1_fx)
        free(tree_arg->dr_coef1_fx);
    if (tree_arg->dr_coef2_fx)
        free(tree_arg->dr_coef2_fx);
    if (tree_arg->dr_coef3_fx)
        free(tree_arg->dr_coef3_fx);

    if (tree_arg->glob_corr_dom_for)
        free(tree_arg->glob_corr_dom_for);
    if (tree_arg->glob_corr_dom_fx)
        free(tree_arg->glob_corr_dom_fx);
    if (tree_arg->glob_corr_for_fx)
        free(tree_arg->glob_corr_for_fx);

    if (tree_arg->void_prm)
    {
        for (i = 0; i < tree_arg->nstp; i++)
        {
            cpd_prm = (CPD_PAY_ARG)(tree_arg->void_prm[i]);

            if (cpd_prm)
            {
                cpd_free_eval_const(cpd_prm->und, cpd_prm->cpd, &(cpd_prm->eval_const));

                free(cpd_prm);
            }
        }

        free(tree_arg->void_prm);
    }

    tree_arg->time              = NULL;
    tree_arg->date              = NULL;
    tree_arg->vol_change        = NULL;
    tree_arg->sig_dom           = NULL;
    tree_arg->sig_for           = NULL;
    tree_arg->sig_fx            = NULL;
    tree_arg->corr_dom_for      = NULL;
    tree_arg->corr_dom_fx       = NULL;
    tree_arg->corr_for_fx       = NULL;
    tree_arg->dom_ifr           = NULL;
    tree_arg->dom_fwd           = NULL;
    tree_arg->dom_var           = NULL;
    tree_arg->for_ifr           = NULL;
    tree_arg->for_fwd           = NULL;
    tree_arg->for_var           = NULL;
    tree_arg->fx_fwd            = NULL;
    tree_arg->fx_var            = NULL;
    tree_arg->void_prm          = NULL;
    tree_arg->is_event          = NULL;
    tree_arg->bar_lvl           = NULL;
    tree_arg->bar_cl            = NULL;
    tree_arg->is_bar            = NULL;
    tree_arg->dr_const_dom      = NULL;
    tree_arg->dr_coef_dom       = NULL;
    tree_arg->dr_const_for      = NULL;
    tree_arg->dr_coef1_for      = NULL;
    tree_arg->dr_coef2_for      = NULL;
    tree_arg->dr_coef3_for      = NULL;
    tree_arg->dr_const_fx       = NULL;
    tree_arg->dr_coef1_fx       = NULL;
    tree_arg->dr_coef2_fx       = NULL;
    tree_arg->dr_coef3_fx       = NULL;
    tree_arg->glob_corr_dom_for = NULL;
    tree_arg->glob_corr_dom_fx  = NULL;
    tree_arg->glob_corr_for_fx  = NULL;

    return NULL;
}

/*	Function for the arguments of the MC function */

Err cpd_fill_mc_arg(
    CPD_UND und,
    CPD_STR cpd,
    /*	Required number of paths */
    long req_pth,
    /*	Do PECS */
    int        do_pecs,
    CPD_MC_ARG mc_arg)
{
    Err         err = NULL;
    CPD_PAY_ARG cpd_prm;
    double *    merge_times = NULL, *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL,
           *corr_dom_for = NULL, *corr_dom_fx = NULL, *corr_for_fx = NULL;

    int nb_merge_times, i, j;

    /*	Initialise */
    mc_arg->time         = NULL;
    mc_arg->date         = NULL;
    mc_arg->dom_ifr      = NULL;
    mc_arg->dom_fwd      = NULL;
    mc_arg->dom_std      = NULL;
    mc_arg->dom_phi      = NULL;
    mc_arg->dom_beta     = NULL;
    mc_arg->dom_bond_pay = NULL;
    mc_arg->dom_beta_pay = NULL;
    mc_arg->for_ifr      = NULL;
    mc_arg->for_fwd      = NULL;
    mc_arg->for_std      = NULL;
    mc_arg->for_phi      = NULL;
    mc_arg->for_beta     = NULL;
    mc_arg->fx_fwd       = NULL;
    mc_arg->fx_std       = NULL;
    mc_arg->dom_for_cov  = NULL;
    mc_arg->dom_fx_cov   = NULL;
    mc_arg->for_fx_cov   = NULL;
    mc_arg->void_prm     = NULL;

    mc_arg->for_fwd_const = NULL;
    mc_arg->for_fwd_lin   = NULL;
    mc_arg->ffx_var       = NULL;

    mc_arg->npaths  = req_pth;
    mc_arg->num_col = 1;
    mc_arg->do_pecs = do_pecs;

    mc_arg->dom_lam      = und->lda_dom;
    mc_arg->for_lam      = und->lda_for;
    mc_arg->corr_times   = und->corr_times;
    mc_arg->corr_dom_for = und->correl_dom_for;
    mc_arg->corr_dom_fx  = und->correl_dom_fx;
    mc_arg->corr_for_fx  = und->correl_for_fx;
    mc_arg->corr_n_times = und->corr_n_times;

    mc_arg->spot_fx = und->spot_fx;
    strcpy(mc_arg->dom_yc, und->dom_yc);
    strcpy(mc_arg->for_yc, und->for_yc);

    /*	Compute time steps */

    /*	Copy event dates */

    mc_arg->nb_dates = cpd->num_calls;
    mc_arg->time     = (double*)calloc(mc_arg->nb_dates, sizeof(double));
    if (!mc_arg->time)
    {
        err = "Memory allocation error (1) in cpd_fill_mc_arg";
        goto FREE_RETURN;
    }
    for (i = 0; i < mc_arg->nb_dates; i++)
    {
        mc_arg->time[i] = cpd->call[i].ex_time;
    }

    /*	Fill time vector */

    /*	Add today if required */
    if (mc_arg->time[0] < -EPS)
    {
        err = "Past event date in cpd_fill_mc_arg";
        goto FREE_RETURN;
    }
    if (mc_arg->time[0] > EPS)
    {
        num_f_add_number((int*)&(mc_arg->nb_dates), &(mc_arg->time), 0.0);
        num_f_sort_vector(mc_arg->nb_dates, mc_arg->time);
        num_f_unique_vector((int*)&(mc_arg->nb_dates), mc_arg->time);
    }

    /*	If only one event today, add empty event */
    if (mc_arg->nb_dates == 1)
    {
        num_f_add_number((int*)&(mc_arg->nb_dates), (&mc_arg->time), 1.0);
    }

    /*	Make dates */
    mc_arg->date = (double*)calloc(mc_arg->nb_dates, sizeof(double));
    if (!mc_arg->date)
    {
        err = "Memory allocation error (2) in cpd_fill_mc_arg";
        goto FREE_RETURN;
    }

    for (i = 0; i < mc_arg->nb_dates; i++)
    {
        mc_arg->date[i] = und->today + DAYS_IN_YEAR * mc_arg->time[i];

        if (i > 0 && mc_arg->date[i] - mc_arg->date[i - 1] >= 1)
        {
            mc_arg->date[i] = (long)(mc_arg->date[i] + 1.0e-08);
            mc_arg->time[i] = YEARS_IN_DAY * (mc_arg->date[i] - und->today);
        }
    }

    /*	Fill distributions */

    if (und->use_beta_dlm)
    {
        mc_arg->dom_fwd       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_std       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_phi       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_beta      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_bond_pay  = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_beta_pay  = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_fwd_const = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_fwd_lin   = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_std       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_phi       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_beta      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->fx_fwd        = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->fx_std        = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->ffx_var       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_for_cov   = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_fx_cov    = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_fx_cov    = (double*)calloc(mc_arg->nb_dates, sizeof(double));

        if (!mc_arg->dom_fwd || !mc_arg->dom_std || !mc_arg->dom_phi || !mc_arg->dom_beta ||
            !mc_arg->dom_bond_pay || !mc_arg->dom_beta_pay || !mc_arg->for_fwd_const ||
            !mc_arg->for_fwd_lin || !mc_arg->for_std || !mc_arg->for_phi || !mc_arg->for_beta ||
            !mc_arg->fx_fwd || !mc_arg->fx_std || !mc_arg->ffx_var || !mc_arg->dom_for_cov ||
            !mc_arg->dom_fx_cov || !mc_arg->for_fx_cov)
        {
            err = "Memory allocation error (4) in cpd_fill_mc_arg";
            goto FREE_RETURN;
        }

        Fx3DBetaDLM_ExpectAndVarGrfn(
            mc_arg->nb_dates,
            mc_arg->time,
            und->model,
            mc_arg->dom_fwd,
            mc_arg->dom_std,
            mc_arg->dom_phi,
            mc_arg->dom_beta,
            mc_arg->for_fwd_const,
            mc_arg->for_fwd_lin,
            mc_arg->for_std,
            mc_arg->for_phi,
            mc_arg->for_beta,
            mc_arg->fx_fwd,
            mc_arg->fx_std,
            mc_arg->ffx_var,
            0,
            mc_arg->dom_for_cov,
            mc_arg->dom_fx_cov,
            mc_arg->for_fx_cov,
            und->grfn_params->min_time);
    }
    else
    {
        mc_arg->dom_ifr      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_fwd      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_std      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_phi      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_beta     = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_bond_pay = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_beta_pay = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_ifr      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_fwd      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_std      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_phi      = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_beta     = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->fx_fwd       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->fx_std       = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_for_cov  = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->dom_fx_cov   = (double*)calloc(mc_arg->nb_dates, sizeof(double));
        mc_arg->for_fx_cov   = (double*)calloc(mc_arg->nb_dates, sizeof(double));

        if (!mc_arg->dom_ifr || !mc_arg->dom_fwd || !mc_arg->dom_std || !mc_arg->dom_phi ||
            !mc_arg->dom_beta || !mc_arg->dom_bond_pay || !mc_arg->dom_beta_pay ||
            !mc_arg->for_ifr || !mc_arg->for_fwd || !mc_arg->for_std || !mc_arg->for_phi ||
            !mc_arg->for_beta || !mc_arg->fx_fwd || !mc_arg->fx_std || !mc_arg->dom_for_cov ||
            !mc_arg->dom_fx_cov || !mc_arg->for_fx_cov)
        {
            err = "Memory allocation error (4) in cpd_fill_mc_arg";
            goto FREE_RETURN;
        }

        merge_times = (double*)calloc(und->sigma_n_rates, sizeof(double));
        if (!merge_times)
        {
            err = "Memory allocation error (5) in cpd_fill_mc_arg";
            goto FREE_RETURN;
        }
        memcpy(merge_times, und->sigma_time_rates, und->sigma_n_rates * sizeof(double));
        nb_merge_times = und->sigma_n_rates;
        num_f_concat_vector(&nb_merge_times, &merge_times, und->sigma_n_fx, und->sigma_time_fx);
        num_f_concat_vector(&nb_merge_times, &merge_times, und->corr_n_times, und->corr_times);
        num_f_sort_vector(nb_merge_times, merge_times);
        num_f_unique_vector(&nb_merge_times, merge_times);

        sig_dom      = (double*)calloc(nb_merge_times, sizeof(double));
        sig_for      = (double*)calloc(nb_merge_times, sizeof(double));
        sig_fx       = (double*)calloc(nb_merge_times, sizeof(double));
        corr_dom_for = (double*)calloc(nb_merge_times, sizeof(double));
        corr_dom_fx  = (double*)calloc(nb_merge_times, sizeof(double));
        corr_for_fx  = (double*)calloc(nb_merge_times, sizeof(double));

        if (!sig_dom || !sig_for || !sig_fx || !corr_dom_for || !corr_dom_fx || !corr_for_fx)
        {
            err = "Memory allocation error (6) in cpd_fill_mc_arg";
            goto FREE_RETURN;
        }

        for (i = nb_merge_times - 1; i >= 0; i--)
        {
            static_find_sig(
                und,
                merge_times[i],
                &(sig_dom[i]),
                &(sig_for[i]),
                &(sig_fx[i]),
                &(corr_dom_for[i]),
                &(corr_dom_fx[i]),
                &(corr_for_fx[i]));
        }

        err = fill_mc_init_corr(
            (long)(mc_arg->date[mc_arg->nb_dates - 1] + 1.0e-08),
            mc_arg->time[mc_arg->nb_dates - 1],
            mc_arg->date,
            mc_arg->time,
            mc_arg->nb_dates,
            merge_times,
            nb_merge_times,
            sig_dom,
            und->lda_dom,
            sig_for,
            und->lda_for,
            sig_fx,
            corr_dom_for,
            corr_dom_fx,
            corr_for_fx,
            und->dom_yc,
            und->for_yc,
            mc_arg->dom_ifr,
            mc_arg->dom_fwd,
            mc_arg->dom_std,
            mc_arg->dom_phi,
            mc_arg->dom_beta,
            mc_arg->dom_bond_pay,
            mc_arg->dom_beta_pay,
            mc_arg->for_ifr,
            mc_arg->for_fwd,
            mc_arg->for_std,
            mc_arg->for_phi,
            mc_arg->for_beta,
            mc_arg->fx_fwd,
            mc_arg->fx_std,
            mc_arg->dom_for_cov,
            mc_arg->dom_fx_cov,
            mc_arg->for_fx_cov);

        free(merge_times);
        free(sig_dom);
        free(sig_for);
        free(sig_fx);
        merge_times = NULL;
        sig_dom     = NULL;
        sig_for     = NULL;
        sig_fx      = NULL;

        if (err)
        {
            goto FREE_RETURN;
        }
    }

    /*	Fill limit conditions (product) */

    mc_arg->void_prm = (void**)calloc(mc_arg->nb_dates, sizeof(void*));

    if (!mc_arg->void_prm)
    {
        err = "Memory allocation error (7) in cpd_fill_mc_arg";
        goto FREE_RETURN;
    }

    j = cpd->num_calls - 1;
    for (i = mc_arg->nb_dates - 1; i >= 0; i--)
    {
        if (j >= 0 && fabs(mc_arg->date[i] - cpd->call[j].ex_date) < 1.0e-04)
        {
            cpd_prm                     = (cpd_pay_arg*)calloc(1, sizeof(cpd_pay_arg));
            cpd_prm->und                = und;
            cpd_prm->cpd                = cpd;
            cpd_prm->eval_const.pd_inst = NULL;

            cpd_prm->call_idx = j;

            if (und->use_beta_dlm)
            {
                cpd_prm->eval_const.beta_dom = mc_arg->dom_beta[i];
                cpd_prm->eval_const.beta_for = mc_arg->for_beta[i];
                cpd_prm->eval_const.phi_dom  = mc_arg->dom_phi[i];
                cpd_prm->eval_const.phi_for  = mc_arg->for_phi[i];
                cpd_prm->eval_const.varFFX   = mc_arg->ffx_var[i];

                err = cpd_fill_eval_const_dlm(und, cpd, j, 1, &(cpd_prm->eval_const));

                mc_arg->dom_bond_pay[i] =
                    log(swp_f_df(cpd->call[j].ex_date, und->model->lTStarDate, und->dom_yc));
                mc_arg->dom_bond_pay[i] += 0.5 * cpd_prm->eval_const.beta_dom *
                                           cpd_prm->eval_const.beta_dom *
                                           cpd_prm->eval_const.phi_dom;
                mc_arg->dom_beta_pay[i] = cpd_prm->eval_const.beta_dom;
            }
            else
            {
                err = cpd_fill_eval_const(und, cpd, j, &(cpd_prm->eval_const));
            }

            if (err)
                goto FREE_RETURN;

            mc_arg->void_prm[i] = (void*)cpd_prm;

            j--;
        }
        else
        {
            mc_arg->void_prm[i] = NULL;
        }
    }

FREE_RETURN:

    if (err)
    {
        cpd_free_mc_arg(mc_arg);
    }

    if (sig_dom)
        free(sig_dom);
    if (sig_for)
        free(sig_for);
    if (sig_fx)
        free(sig_fx);
    if (corr_dom_for)
        free(corr_dom_for);
    if (corr_dom_fx)
        free(corr_dom_fx);
    if (corr_for_fx)
        free(corr_for_fx);

    return err;
}

Err cpd_free_mc_arg(CPD_MC_ARG mc_arg)
{
    int         i;
    CPD_PAY_ARG cpd_prm;

    if (mc_arg->time)
        free(mc_arg->time);
    if (mc_arg->date)
        free(mc_arg->date);
    if (mc_arg->dom_ifr)
        free(mc_arg->dom_ifr);
    if (mc_arg->dom_fwd)
        free(mc_arg->dom_fwd);
    if (mc_arg->dom_std)
        free(mc_arg->dom_std);
    if (mc_arg->dom_phi)
        free(mc_arg->dom_phi);
    if (mc_arg->dom_beta)
        free(mc_arg->dom_beta);
    if (mc_arg->dom_bond_pay)
        free(mc_arg->dom_bond_pay);
    if (mc_arg->dom_beta_pay)
        free(mc_arg->dom_beta_pay);
    if (mc_arg->for_ifr)
        free(mc_arg->for_ifr);
    if (mc_arg->for_fwd)
        free(mc_arg->for_fwd);
    if (mc_arg->for_std)
        free(mc_arg->for_std);
    if (mc_arg->for_phi)
        free(mc_arg->for_phi);
    if (mc_arg->for_beta)
        free(mc_arg->for_beta);
    if (mc_arg->fx_fwd)
        free(mc_arg->fx_fwd);
    if (mc_arg->fx_std)
        free(mc_arg->fx_std);
    if (mc_arg->dom_for_cov)
        free(mc_arg->dom_for_cov);
    if (mc_arg->dom_fx_cov)
        free(mc_arg->dom_fx_cov);
    if (mc_arg->for_fx_cov)
        free(mc_arg->for_fx_cov);

    if (mc_arg->for_fwd_const)
        free(mc_arg->for_fwd_const);
    if (mc_arg->for_fwd_lin)
        free(mc_arg->for_fwd_lin);
    if (mc_arg->ffx_var)
        free(mc_arg->ffx_var);

    if (mc_arg->void_prm)
    {
        for (i = 0; i < mc_arg->nb_dates; i++)
        {
            cpd_prm = (CPD_PAY_ARG)(mc_arg->void_prm[i]);

            if (cpd_prm)
            {
                cpd_free_eval_const(cpd_prm->und, cpd_prm->cpd, &(cpd_prm->eval_const));

                free(cpd_prm);
            }
        }

        free(mc_arg->void_prm);
    }

    mc_arg->time          = NULL;
    mc_arg->date          = NULL;
    mc_arg->dom_ifr       = NULL;
    mc_arg->dom_fwd       = NULL;
    mc_arg->dom_std       = NULL;
    mc_arg->dom_phi       = NULL;
    mc_arg->dom_beta      = NULL;
    mc_arg->dom_bond_pay  = NULL;
    mc_arg->dom_beta_pay  = NULL;
    mc_arg->for_ifr       = NULL;
    mc_arg->for_fwd       = NULL;
    mc_arg->for_std       = NULL;
    mc_arg->for_phi       = NULL;
    mc_arg->for_beta      = NULL;
    mc_arg->fx_fwd        = NULL;
    mc_arg->fx_std        = NULL;
    mc_arg->dom_for_cov   = NULL;
    mc_arg->dom_fx_cov    = NULL;
    mc_arg->for_fx_cov    = NULL;
    mc_arg->for_fwd_const = NULL;
    mc_arg->for_fwd_lin   = NULL;
    mc_arg->for_fx_cov    = NULL;
    mc_arg->ffx_var       = NULL;

    return NULL;
}

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

Err cpd_fill_check_all_struct(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund, 1: calibrate */
    /*		if calib */
    double fx_spot,
    long   fx_spot_date,
    int    dom_calib,      /*	Calibrate domestic underlying */
    char*  dom_und,        /*	If no, domestic underlying to be used */
    char*  dom_yc,         /*	Domestic yc */
    char*  dom_vc,         /*	Domestic vc (only if calib) */
    char*  dom_ref,        /*	Domestic ref rate (only if calib) */
    char*  dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char*  dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,        /*	Domestic lambda */
    int    for_calib,      /*	Same for foreign */
    char*  for_und,
    char*  for_yc,
    char*  for_vc,
    char*  for_ref,
    char*  for_swap_freq,
    char*  for_swap_basis,
    double for_lam,

    double min_fact,  /*	Maximum down jump on variance */
    double max_fact,  /*	Maximum up jump on variance */
    int    use_jumps, /*	Allow vol term structure to jump */

    double*          corr_times,
    double*          correl_dom_for, /*	Correlations */
    double*          correl_dom_fx,
    double*          correl_for_fx,
    long             corr_n_times,
    CPDBETADLMPARAMS cpd_dlm_params,
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char*   vol_curve_name,
                           double  start_date,
                           double  end_date,
                           double  cash_strike,
                           int     zero,
                           char*   ref_rate_name,
                           double* vol,
                           double* power),
    /*	Fx vol from the market */
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,
    /*		if no calilb */
    char* fx3dund,
    /*	The structure */
    /*		funding */
    double  fund_not,
    int     fund_ccy, /*	0: domestic, 1: foreign */
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    /*		pd */
    double  pd_not,
    int     pd_ncpn,
    long*   pd_fix,
    long*   pd_start,
    long*   pd_pay,
    char**  pd_basis,
    double* pd_alpha,
    double* pd_beta,
    int*    pd_floored,
    double* pd_floor,
    int*    pd_capped,
    double* pd_cap,

    //		pd interp coupon:
    int      use_cpn_opt_str,
    int*     pd_num_strikes,
    double*  pd_wcst,
    double*  pd_wspot,
    double** pd_strikes,
    double** pd_weights,

    /*		pd not refund */
    long*  pd_not_ref_fix,
    double pd_not_ref_alpha,
    double pd_not_ref_beta,
    int    pd_not_ref_floored,
    double pd_not_ref_floor,
    int    pd_not_ref_capped,
    double pd_not_ref_cap,

    //		pd interp notional:
    int      use_not_opt_str,
    int*     pd_not_num_strikes,
    double*  pd_not_wcst,
    double*  pd_not_wspot,
    double** pd_not_strikes,
    double** pd_not_weights,

    /*		calls */
    int*    call_type, /*	0: Call, 1: KO */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* barrier,  /*	KO only */
    int*    bar_type, /*	0: up and in, 1: down and in */
    double* fees,     /*  fees if deal is called in domestic currency */
    int     TARN_Do,
    /*	Numerical params */
    long   req_stp,
    long   req_pth,
    int    do_pecs,
    int    forcetree,
    int    do_optim,    /*	If equal to 1 then the call are replaced by optimal KO	*/
    int    force_optim, /*	If equal to 1 then all call will be replaced by optimal KO	*/
    int    fx_bound,    /*	If equal to 1 then optimisation on the Fx, on the IV otherwise	*/
    int    use_bound,
    double smooth,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	Results */
    CPD_STR cpd,
    CPD_UND und,
    int*    call_feat, /*	0: No callable feature to be valued
                                       1: Callable feature to be valued through tree
                                       2: Callable feature (KO) to be valued through MC */
    CPD_TREE_ARG tree_arg,
    CPD_MC_ARG   mc_arg,
    double       dom_vol_shift,
    double       for_vol_shift,
    double       fx_vol_shift,
    int          skip_fill)
{
    Err     err        = NULL;
    double *fund_not_s = NULL, *pd_not_s = NULL;
    int     i;

    memset(tree_arg, 0, sizeof(cpd_tree_arg));
    memset(mc_arg, 0, sizeof(cpd_mc_arg));

    if (skip_fill == 0)
    {
        /*	Initialisation */
        memset(cpd, 0, sizeof(cpd_str));
        memset(und, 0, sizeof(cpd_und));

        und->today = today;

        /*	Funding leg */

        cpd->fund_leg = (PD_FUND_LEG)calloc(1, sizeof(pd_fund_leg));
        if (!cpd->fund_leg)
        {
            err = "Memory allocation error (1) in cpd_fill_check_all_struct";
            goto FREE_RETURN;
        }

        err = cpd_fill_fund_leg(
            und->today,
            eod_fix_flag,
            fund_not,
            fund_ccy,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            cpd->fund_leg);
        if (err)
        {
            goto FREE_RETURN;
        }

        /*	Exotic leg */

        cpd->pd_leg = (PD_EXO_LEG)calloc(1, sizeof(pd_exo_leg));
        if (!cpd->pd_leg)
        {
            err = "Memory allocation error (2) in cpd_fill_check_all_struct";
            goto FREE_RETURN;
        }

        err = cpd_fill_exo_leg(
            und->today,
            eod_fix_flag,
            pd_not,
            pd_ncpn,
            pd_fix,
            pd_start,
            pd_pay,
            pd_basis,
            pd_alpha,
            pd_beta,
            pd_floored,
            pd_floor,
            pd_capped,
            pd_cap,
            use_cpn_opt_str,
            pd_num_strikes,
            pd_wcst,
            pd_wspot,
            pd_strikes,
            pd_weights,
            pd_not_ref_fix[ncall],
            pd_not_ref_alpha,
            pd_not_ref_beta,
            pd_not_ref_floored,
            pd_not_ref_floor,
            pd_not_ref_capped,
            pd_not_ref_cap,
            use_not_opt_str,
            use_not_opt_str ? pd_not_num_strikes[ncall] : 0,
            use_not_opt_str ? pd_not_wcst[ncall] : 0.0,
            use_not_opt_str ? pd_not_wspot[ncall] : 0.0,
            use_not_opt_str ? pd_not_strikes[ncall] : NULL,
            use_not_opt_str ? pd_not_weights[ncall] : NULL,
            cpd->pd_leg);

        if (err)
            goto FREE_RETURN;

        if (1 == 2 && TARN_Do)
        {
            // We have to add one coupon for both funding and exotic leg
            cpd_fill_exo_and_fund_leg_TARN(und->today, cpd);
        }

        /*	Calls */

        if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag)
        {
            fund_not_s = (double*)calloc(ncall, sizeof(double));
            pd_not_s   = (double*)calloc(ncall, sizeof(double));

            if (!fund_not_s || !pd_not_s)
            {
                err = serror("Memory failure");
                goto FREE_RETURN;
            }

            for (i = 0; i < ncall; i++)
            {
                fund_not_s[i] = fund_not;
                pd_not_s[i]   = pd_not;
            }

            err = cpd_fill_calls(
                und->today,
                eod_ex_flag,
                ncall,
                call_type,
                pay_rec,
                ex_date,
                set_date,
                pd_not_ref_fix,
                fund_not_s,
                pd_not_s,
                use_not_opt_str,
                pd_not_num_strikes,
                pd_not_wcst,
                pd_not_wspot,
                pd_not_strikes,
                pd_not_weights,
                barrier,
                bar_type,
                fees,
                TARN_Do,
                smooth,
                fx_bound,
                use_bound,
                force_optim,
                cpd);

            if (err)
                goto FREE_RETURN;
        }
        else
        {
            cpd->num_calls = 0;
            cpd->call      = NULL;
        }

        /*	Underlying */
        if (use_calib)
        {
            err = cpd_calib_und(
                today,
                eod_ex_flag,
                fx_spot,
                fx_spot_date,
                dom_calib,
                dom_und,
                dom_yc,
                dom_vc,
                dom_ref,
                dom_swap_freq,
                dom_swap_basis,
                dom_lam,
                for_calib,
                for_und,
                for_yc,
                for_vc,
                for_ref,
                for_swap_freq,
                for_swap_basis,
                for_lam,
                min_fact,
                max_fact,
                use_jumps,
                corr_times,
                correl_dom_for,
                correl_dom_fx,
                correl_for_fx,
                corr_n_times,
                cpd_dlm_params,
                cpd,
                get_ir_cash_vol,
                fx_mkt_vol_date,
                fx_mkt_vol,
                num_fx_mkt_vol,
                und,
                dom_vol_shift,
                for_vol_shift,
                fx_vol_shift);
        }
        else
        {
            err = cpd_fill_und(
                fx3dund, und, cpd_dlm_params, dom_vol_shift, for_vol_shift, fx_vol_shift);
        }

        if (err)
            goto FREE_RETURN;
    }

    /*	Tree */
    if (cpd->num_calls > 0 && cpd->call)
    {
        /* First we check if we want to force the tree */
        if (cpd->type == 1 && forcetree)
            cpd->type = -1;

        switch (cpd->type)
        {
        case -1:
            /*	Callable KO	in the Tree		*/
            if (do_optim)
            {
                if (!cpd_dlm_params->use_beta_dlm || fabs(und->model->dAlpha) < TINY)
                    err = cpd_fill_mc_arg(und, cpd, req_pth, do_pecs, mc_arg);

                /*	Barrier handling			*/
                for (i = 0; i < cpd->num_calls; i++)
                {
                    if (cpd->call[i].call_type == 1)
                    {
                        cpd->call[i].smooth  = cpd->call[i].smooth / cpd->call[i].barrier;
                        cpd->call[i].barrier = log(cpd->call[i].barrier / und->spot_fx);
                    }
                }
            }
            else
            {
                if (!cpd_dlm_params->use_beta_dlm || fabs(und->model->dAlpha) < TINY)
                    err = cpd_fill_tree_arg(und, cpd, req_stp, tree_arg);
            }

            if (err)
            {
                goto FREE_RETURN;
            }

            *call_feat = 1;
            break;

        case 0:
            /*	Callable only in the Tree	*/
            if (!cpd_dlm_params->use_beta_dlm || fabs(und->model->dAlpha) < TINY)
                err = cpd_fill_tree_arg(und, cpd, req_stp, tree_arg);
            if (err)
            {
                goto FREE_RETURN;
            }
            *call_feat = 1;
            break;

        case 1:
            /*	KO in Monte Carlo			*/
            if (!cpd_dlm_params->use_beta_dlm || fabs(und->model->dAlpha) < TINY)
                err = cpd_fill_mc_arg(und, cpd, req_pth, do_pecs, mc_arg);
            if (err)
            {
                goto FREE_RETURN;
            }

            /*	Barrier handling			*/
            if (skip_fill == 0)
            {
                for (i = 0; i < cpd->num_calls; i++)
                {
                    if (cpd->call[i].fx_bound)
                    {
                        cpd->call[i].smooth  = cpd->call[i].smooth / cpd->call[i].barrier;
                        cpd->call[i].barrier = log(cpd->call[i].barrier / und->spot_fx);
                    }
                }
            }

            *call_feat = 2;
        }
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (err)
    {
        cpd_free_all_struct(cpd, und, tree_arg, mc_arg);
    }

    if (fund_not_s)
    {
        free(fund_not_s);
    }

    if (pd_not_s)
    {
        free(pd_not_s);
    }

    return err;
}

/*	Free all structures */
Err cpd_free_all_struct(CPD_STR cpd, CPD_UND und, CPD_TREE_ARG tree_arg, CPD_MC_ARG mc_arg)
{
    cpd_free_und(und);
    cpd_free_calls(cpd);

    if (cpd->fund_leg)
    {
        cpd_free_fund_leg(cpd->fund_leg);
        free(cpd->fund_leg);
    }

    if (cpd->pd_leg)
    {
        cpd_free_exo_leg(cpd->pd_leg);
        free(cpd->pd_leg);
    }

    cpd_free_tree_arg(tree_arg);
    cpd_free_mc_arg(mc_arg);

    return NULL;
}

/*	Payoff function */
/*	---------------	*/

#define POS_VAL(X) ((X) > 0 ? (X) : 0)

#define CALL_VAL(FWD, STRIKE, D, S) ((FWD)*norm((D) + (S)) - (STRIKE)*norm((D)))

#define PUT_VAL(FWD, STRIKE, D, S) (-(FWD)*norm(-(D) - (S)) + (STRIKE)*norm(-(D)))

#define OPT_VAL_MACRO(TYPE, FWD, STRIKE, STD, HALF_STD)                                         \
    ((TYPE) == 0                                                                                \
         ? 0.0                                                                                  \
         : ((TYPE) == 1                                                                         \
                ? POS_VAL((FWD) - (STRIKE))                                                     \
                : ((TYPE) == 2 ? POS_VAL((STRIKE) - (FWD))                                      \
                               : ((TYPE) == 3 ? CALL_VAL(                                       \
                                                    (FWD),                                      \
                                                    (STRIKE),                                   \
                                                    log((FWD) / (STRIKE)) / (STD) - (HALF_STD), \
                                                    (STD))                                      \
                                              : PUT_VAL(                                        \
                                                    (FWD),                                      \
                                                    (STRIKE),                                   \
                                                    log((FWD) / (STRIKE)) / (STD) - (HALF_STD), \
                                                    (STD))))))

#define OPT_LOGVAL_MACRO(TYPE, FWD, LOGFWD, STRIKE, LOGSTRIKE, STD, HALF_STD)                      \
    ((TYPE) == 0                                                                                   \
         ? 0.0                                                                                     \
         : ((TYPE) == 1                                                                            \
                ? POS_VAL((FWD) - (STRIKE))                                                        \
                : ((TYPE) == 2 ? POS_VAL((STRIKE) - (FWD))                                         \
                               : ((TYPE) == 3 ? CALL_VAL(                                          \
                                                    (FWD),                                         \
                                                    (STRIKE),                                      \
                                                    ((LOGFWD) - (LOGSTRIKE)) / (STD) - (HALF_STD), \
                                                    (STD))                                         \
                                              : PUT_VAL(                                           \
                                                    (FWD),                                         \
                                                    (STRIKE),                                      \
                                                    ((LOGFWD) - (LOGSTRIKE)) / (STD) - (HALF_STD), \
                                                    (STD))))))

Err cpd_payoff_4_3dfx_tree(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    /* Nodes data */
    long n1,
    long n2,
    long n3,
    /* i: d1, j: d2, k: d3, l = {0: xDom, 1: xFor, 2: log (Fx/Fx0)} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call, next_call;
    CPD_UND        und;
    PD_EXO_CPN     cpn;
    CPD_EVAL_CONST eval_const;
    int            call_idx;

    int i, j, k, l, str_idx, opt_type;

    /*
    double			***svi,
                                    **svij,
                                    *svijk, */
    double fx;

    double fund_leg, df, fwd, floor, cap, opt_string, pd_leg, iv, fee;

    int num_fund_cpn, num_pd_cpn, eval_redemption;

    Err err = NULL;

    /* For Optimisation */
    double fx_I, fx_J, fx_K, fx_DI, fx_DJ, fx_DK;

    double start_I, start_J, start_K, start_DI, start_DJ, start_DK;

    double in_not_fund_I, in_not_fund_J, in_not_fund_K, in_not_fund_DI, in_not_fund_DJ,
        in_not_fund_DK;

    double in_not_pd_I, in_not_pd_J, in_not_pd_K, in_not_pd_DI, in_not_pd_DJ, in_not_pd_DK;

    double in_not_pd_fwd_I, in_not_pd_fwd_J, in_not_pd_fwd_K, in_not_pd_fwd_DI, in_not_pd_fwd_DJ,
        in_not_pd_fwd_DK;

    double next_start_fund_I, next_start_fund_J, next_start_fund_K, next_start_fund_DI,
        next_start_fund_DJ, next_start_fund_DK;

    double next_fund_I, next_fund_J, next_fund_K, next_fund_DI, next_fund_DJ, next_fund_DK;

    double next_pd_I, next_pd_J, next_pd_K, next_pd_DI, next_pd_DJ, next_pd_DK;

    double fin_not_disc_I, fin_not_disc_J, fin_not_disc_K, fin_not_disc_DI, fin_not_disc_DJ,
        fin_not_disc_DK;

    double fin_not_fwd_I, fin_not_fwd_J, fin_not_fwd_K, fin_not_fwd_DI, fin_not_fwd_DJ,
        fin_not_fwd_DK;

    /*	Get the event */
    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    /*	Calculate number of coupons to be computed */
    if (call_idx == cpd->num_calls - 1)
    {
        /*	Last call: eval all + redemption */
        num_fund_cpn    = call->num_fund_cpn;
        num_pd_cpn      = call->num_pd_cpn;
        eval_redemption = 1;
    }
    else
    {
        /*	Not the last call: only eval coupons up to next call date */
        next_call       = cpd->call + (call_idx + 1);
        num_fund_cpn    = next_call->fund_idx - call->fund_idx;
        num_pd_cpn      = next_call->pd_idx - call->pd_idx;
        eval_redemption = 0;
    }

    /*	Eval payoff
            ----------- */

    /* ---------------------------------------------------------------------------------------------------------------------------
            For Optimisation
       ---------------------------------------------------------------------------------------------------------------------------
     */

    /* For Fx */
    fx_I = spot_fx * exp(sv[0][0][0][2]);

    fx_J = fx_I;
    fx_K = fx_I;

    if (n1 > 0)
        fx_DI = exp(sv[1][0][0][2] - sv[0][0][0][2]);
    else
        fx_DI = 1.0;

    if (n2 > 0)
        fx_DJ = exp(sv[0][1][0][2] - sv[0][0][0][2]);
    else
        fx_DJ = 1.0;

    if (n3 > 0)
        fx_DK = exp(sv[0][0][1][2] - sv[0][0][0][2]);
    else
        fx_DK = 1.0;

    /* for funding start */
    start_I =
        eval_const->start_dff *
        exp(-eval_const->start_gam * sv[0][0][0][eval_const->fund_var] - eval_const->start_gam2) *
        cpd->fund_leg->notional;

    start_J = start_I;
    start_K = start_I;

    if (n1 > 0)
        start_DI =
            exp(-eval_const->start_gam *
                (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    else
        start_DI = 1.0;

    if (n2 > 0)
        start_DJ =
            exp(-eval_const->start_gam *
                (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    else
        start_DJ = 1.0;

    if (n3 > 0)
        start_DK =
            exp(-eval_const->start_gam *
                (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    else
        start_DK = 1.0;

    /* For initial notional */
    in_not_fund_I = eval_const->in_not_fund_dff *
                    exp(-eval_const->in_not_fund_gam * sv[0][0][0][eval_const->fund_var] -
                        eval_const->in_not_fund_gam2) *
                    call->fund_not_amt;

    in_not_fund_J = in_not_fund_I;
    in_not_fund_K = in_not_fund_I;

    if (n1 > 0)
        in_not_fund_DI =
            exp(-eval_const->in_not_fund_gam *
                (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    else
        in_not_fund_DI = 1.0;

    if (n2 > 0)
        in_not_fund_DJ =
            exp(-eval_const->in_not_fund_gam *
                (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    else
        in_not_fund_DJ = 1.0;

    if (n3 > 0)
        in_not_fund_DK =
            exp(-eval_const->in_not_fund_gam *
                (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    else
        in_not_fund_DK = 1.0;

    /*	for fee */
    in_not_pd_I = eval_const->in_not_pd_dff *
                  exp(-eval_const->in_not_pd_gam * sv[0][0][0][0] - eval_const->in_not_pd_gam2);

    in_not_pd_J = in_not_pd_I;
    in_not_pd_K = in_not_pd_I;

    if (n1 > 0)
        in_not_pd_DI = exp(-eval_const->in_not_pd_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
    else
        in_not_pd_DI = 1.0;

    if (n2 > 0)
        in_not_pd_DJ = exp(-eval_const->in_not_pd_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
    else
        in_not_pd_DJ = 1.0;

    if (n3 > 0)
        in_not_pd_DK = exp(-eval_const->in_not_pd_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
    else
        in_not_pd_DK = 1.0;

    // In case of exotic redemption at call date:
    if (call->use_opt_str)
    {
        /*	For Fwd fx */
        in_not_pd_fwd_I =
            eval_const->in_not_pd_fwd_dff_ratio *
            exp(-eval_const->in_not_pd_fwd_for_gam * sv[0][0][0][1] +
                eval_const->in_not_pd_fwd_dom_gam * sv[0][0][0][0] + eval_const->in_not_pd_fwd_cvx);

        in_not_pd_fwd_J = in_not_pd_fwd_I;
        in_not_pd_fwd_K = in_not_pd_fwd_I;

        if (n1 > 0)
            in_not_pd_fwd_DI =
                exp(-eval_const->in_not_pd_fwd_for_gam * (sv[1][0][0][1] - sv[0][0][0][1]) +
                    eval_const->in_not_pd_fwd_dom_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        else
            in_not_pd_fwd_DI = 1.0;

        if (n2 > 0)
            in_not_pd_fwd_DJ =
                exp(-eval_const->in_not_pd_fwd_for_gam * (sv[0][1][0][1] - sv[0][0][0][1]) +
                    eval_const->in_not_pd_fwd_dom_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        else
            in_not_pd_fwd_DJ = 1.0;

        if (n3 > 0)
            in_not_pd_fwd_DK =
                exp(-eval_const->in_not_pd_fwd_for_gam * (sv[0][0][1][1] - sv[0][0][0][1]) +
                    eval_const->in_not_pd_fwd_dom_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        else
            in_not_pd_fwd_DK = 1.0;
    }

    /*	Next exchange if relevant */
    if (!eval_redemption)
    {
        next_start_fund_I =
            eval_const->next_start_fund_dff *
            exp(-eval_const->next_start_fund_gam * sv[0][0][0][eval_const->fund_var] -
                eval_const->next_start_fund_gam2) *
            cpd->fund_leg->notional;

        next_start_fund_J = next_start_fund_I;
        next_start_fund_K = next_start_fund_I;

        if (n1 > 0)
            next_start_fund_DI =
                exp(-eval_const->next_start_fund_gam *
                    (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            next_start_fund_DI = 1.0;

        if (n2 > 0)
            next_start_fund_DJ =
                exp(-eval_const->next_start_fund_gam *
                    (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            next_start_fund_DJ = 1.0;

        if (n3 > 0)
            next_start_fund_DK =
                exp(-eval_const->next_start_fund_gam *
                    (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            next_start_fund_DK = 1.0;

        next_fund_I = eval_const->next_fund_dff *
                      exp(-eval_const->next_fund_gam * sv[0][0][0][eval_const->fund_var] -
                          eval_const->next_fund_gam2) *
                      next_call->fund_not_amt;

        next_fund_J = next_fund_I;
        next_fund_K = next_fund_I;

        if (n1 > 0)
            next_fund_DI =
                exp(-eval_const->next_fund_gam *
                    (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            next_fund_DI = 1.0;

        if (n2 > 0)
            next_fund_DJ =
                exp(-eval_const->next_fund_gam *
                    (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            next_fund_DJ = 1.0;

        if (n3 > 0)
            next_fund_DK =
                exp(-eval_const->next_fund_gam *
                    (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            next_fund_DK = 1.0;

        next_pd_I = eval_const->next_pd_dff *
                    exp(-eval_const->next_pd_gam * sv[0][0][0][0] - eval_const->next_pd_gam2);

        next_pd_J = next_pd_I;
        next_pd_K = next_pd_I;

        if (n1 > 0)
            next_pd_DI = exp(-eval_const->next_pd_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        else
            next_pd_DI = 1.0;

        if (n2 > 0)
            next_pd_DJ = exp(-eval_const->next_pd_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        else
            next_pd_DJ = 1.0;

        if (n3 > 0)
            next_pd_DK = exp(-eval_const->next_pd_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        else
            next_pd_DK = 1.0;
    }
    else
    /*	Final notional */
    {
        /* For df */
        fin_not_disc_I =
            eval_const->fin_not_disc_dff *
            exp(-eval_const->fin_not_disc_gam * sv[0][0][0][0] - eval_const->fin_not_disc_gam2);

        fin_not_disc_J = fin_not_disc_I;
        fin_not_disc_K = fin_not_disc_I;

        if (n1 > 0)
            fin_not_disc_DI =
                exp(-eval_const->fin_not_disc_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        else
            fin_not_disc_DI = 1.0;

        if (n2 > 0)
            fin_not_disc_DJ =
                exp(-eval_const->fin_not_disc_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        else
            fin_not_disc_DJ = 1.0;

        if (n3 > 0)
            fin_not_disc_DK =
                exp(-eval_const->fin_not_disc_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        else
            fin_not_disc_DK = 1.0;

        /*	For Fwd fx */
        fin_not_fwd_I =
            eval_const->fin_not_fwd_dff_ratio *
            exp(-eval_const->fin_not_fwd_for_gam * sv[0][0][0][1] +
                eval_const->fin_not_fwd_dom_gam * sv[0][0][0][0] + eval_const->fin_not_fwd_cvx);

        fin_not_fwd_J = fin_not_fwd_I;
        fin_not_fwd_K = fin_not_fwd_I;

        if (n1 > 0)
            fin_not_fwd_DI =
                exp(-eval_const->fin_not_fwd_for_gam * (sv[1][0][0][1] - sv[0][0][0][1]) +
                    eval_const->fin_not_fwd_dom_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        else
            fin_not_fwd_DI = 1.0;

        if (n2 > 0)
            fin_not_fwd_DJ =
                exp(-eval_const->fin_not_fwd_for_gam * (sv[0][1][0][1] - sv[0][0][0][1]) +
                    eval_const->fin_not_fwd_dom_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        else
            fin_not_fwd_DJ = 1.0;

        if (n3 > 0)
            fin_not_fwd_DK =
                exp(-eval_const->fin_not_fwd_for_gam * (sv[0][0][1][1] - sv[0][0][0][1]) +
                    eval_const->fin_not_fwd_dom_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        else
            fin_not_fwd_DK = 1.0;
    }

    for (l = 0; l < num_fund_cpn; l++)
    {
        /* for funding coupon */
        eval_const->fund_I[l] = eval_const->fund_dff[l] *
                                exp(-eval_const->fund_gam[l] * sv[0][0][0][eval_const->fund_var] -
                                    eval_const->fund_gam2[l]) *
                                cpd->fund_leg->cpn[call->fund_idx + l].cpn;

        eval_const->fund_J[l] = eval_const->fund_I[l];

        eval_const->fund_K[l] = eval_const->fund_I[l];

        if (n1 > 0)
            eval_const->fund_DI[l] =
                exp(-eval_const->fund_gam[l] *
                    (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            eval_const->fund_DI[l] = 1.0;

        if (n2 > 0)
            eval_const->fund_DJ[l] =
                exp(-eval_const->fund_gam[l] *
                    (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            eval_const->fund_DJ[l] = 1.0;

        if (n3 > 0)
            eval_const->fund_DK[l] =
                exp(-eval_const->fund_gam[l] *
                    (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        else
            eval_const->fund_DK[l] = 1.0;
    }

    for (l = 0; l < num_pd_cpn; l++)
    {
        /*	Coupons */
        eval_const->pd_disc_I[l] =
            eval_const->pd_disc_dff[l] *
            exp(-eval_const->pd_disc_gam[l] * sv[0][0][0][0] - eval_const->pd_disc_gam2[l]);

        eval_const->pd_disc_J[l] = eval_const->pd_disc_I[l];

        eval_const->pd_disc_K[l] = eval_const->pd_disc_I[l];

        if (n1 > 0)
            eval_const->pd_disc_DI[l] =
                exp(-eval_const->pd_disc_gam[l] * (sv[1][0][0][0] - sv[0][0][0][0]));
        else
            eval_const->pd_disc_DI[l] = 1.0;

        if (n2 > 0)
            eval_const->pd_disc_DJ[l] =
                exp(-eval_const->pd_disc_gam[l] * (sv[0][1][0][0] - sv[0][0][0][0]));
        else
            eval_const->pd_disc_DJ[l] = 1.0;

        if (n3 > 0)
            eval_const->pd_disc_DK[l] =
                exp(-eval_const->pd_disc_gam[l] * (sv[0][0][1][0] - sv[0][0][0][0]));
        else
            eval_const->pd_disc_DK[l] = 1.0;

        /* for the fwd */

        eval_const->pd_fwd_I[l] =
            eval_const->pd_fwd_dff_ratio[l] *
            exp(-eval_const->pd_fwd_for_gam[l] * sv[0][0][0][1] +
                eval_const->pd_fwd_dom_gam[l] * sv[0][0][0][0] + eval_const->pd_fwd_cvx[l]);

        eval_const->pd_fwd_J[l] = eval_const->pd_fwd_I[l];

        eval_const->pd_fwd_K[l] = eval_const->pd_fwd_I[l];

        if (n1 > 0)
            eval_const->pd_fwd_DI[l] =
                exp(-eval_const->pd_fwd_for_gam[l] * (sv[1][0][0][1] - sv[0][0][0][1]) +
                    eval_const->pd_fwd_dom_gam[l] * (sv[1][0][0][0] - sv[0][0][0][0]));
        else
            eval_const->pd_fwd_DI[l] = 1.0;

        if (n2 > 0)
            eval_const->pd_fwd_DJ[l] =
                exp(-eval_const->pd_fwd_for_gam[l] * (sv[0][1][0][1] - sv[0][0][0][1]) +
                    eval_const->pd_fwd_dom_gam[l] * (sv[0][1][0][0] - sv[0][0][0][0]));
        else
            eval_const->pd_fwd_DJ[l] = 1.0;

        if (n3 > 0)
            eval_const->pd_fwd_DK[l] =
                exp(-eval_const->pd_fwd_for_gam[l] * (sv[0][0][1][1] - sv[0][0][0][1]) +
                    eval_const->pd_fwd_dom_gam[l] * (sv[0][0][1][0] - sv[0][0][0][0]));
        else
            eval_const->pd_fwd_DK[l] = 1.0;
    }

    /* ---------------------------------------------------------------------------------------------------------------------------
            End of init of the Optimisation
       ---------------------------------------------------------------------------------------------------------------------------
     */

    for (i = 0; i < n1; i++)
    {
        /* svi = sv[i]; */

        for (j = 0; j < n2; j++)
        {
            /* svij = svi[j]; */
            for (k = 0; k < n3; k++)
            {
                /* svijk = svij[k]; */

                /* fx = spot_fx * exp (svijk[2]); */
                fx = fx_K;

                /*	PV of funding leg */

                fund_leg = 0.0;

                /*	Libor */

                /*fund_leg += eval_const->start_dff
                 * exp (- eval_const->start_gam * svijk[eval_const->fund_var] -
                 * eval_const->start_gam2) cpd->fund_leg->notional; */

                fund_leg += start_K;

                /*	Coupons */

                for (l = 0; l < num_fund_cpn; l++)
                {
                    /*fund_leg +=
                            eval_const->fund_dff[l]
            * exp (- eval_const->fund_gam[l] * svijk[eval_const->fund_var] - eval_const->fund_gam2[l])
                            * cpd->fund_leg->cpn[call->fund_idx+l].cpn; */
                    fund_leg += eval_const->fund_K[l];
                }

                /*	Initial notional */

                /*fund_leg -= eval_const->in_not_fund_dff
                 * exp (- eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                 * eval_const->in_not_fund_gam2) call->fund_not_amt; */

                fund_leg -= in_not_fund_K;

                /*	Next exchange if relevant */
                if (!eval_redemption)
                {
                    /* fund_leg -= eval_const->next_start_fund_dff
                     * exp (- eval_const->next_start_fund_gam * svijk[eval_const->fund_var] -
                     * eval_const->next_start_fund_gam2) cpd->fund_leg->notional; */
                    fund_leg -= next_start_fund_K;

                    /*fund_leg += eval_const->next_fund_dff
                     * exp (- eval_const->next_fund_gam * svijk[eval_const->fund_var] -
                     * eval_const->next_fund_gam2) next_call->fund_not_amt; */
                    fund_leg += next_fund_K;
                }

                /*	Fx */

                if (cpd->fund_leg->dom_for != 0)
                {
                    fund_leg *= fx;
                }

                /*	PV of pd leg */

                pd_leg = 0.0;

                /*	Coupons */
                for (l = 0; l < num_pd_cpn; l++)
                {
                    /*	Coupon access */
                    cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                    /*	Discount */
                    /*df =
                            eval_const->pd_disc_dff[l]
                            * exp (- eval_const->pd_disc_gam[l] * svijk[0] - eval_const->pd_disc_gam2[l]);
                    */
                    df = eval_const->pd_disc_K[l];

                    /*	Fwd fx */
                    /*fwd = fx
                            * eval_const->pd_fwd_dff_ratio[l]
                            * exp (	-	eval_const->pd_fwd_for_gam[l] * svijk[1]
                                            +	eval_const->pd_fwd_dom_gam[l] * svijk[0]
                                            +	eval_const->pd_fwd_cvx[l]); */
                    fwd = fx * eval_const->pd_fwd_K[l];

                    // interp coupon case:
                    if (cpn->use_opt_str)
                    {
                        opt_string = 0.0;
                        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++)
                        {
                            if (fabs(cpn->weights[str_idx]) > 1e-16)
                            {
                                if (eval_const->pd_std[l] > 1e-16 && cpn->strikes[str_idx] > 1e-16)
                                    opt_type = 3;
                                else
                                    opt_type = 1;

                                opt_string +=
                                    cpn->weights[str_idx] * OPT_VAL_MACRO(
                                                                opt_type,
                                                                fwd,
                                                                cpn->strikes[str_idx],
                                                                eval_const->pd_std[l],
                                                                eval_const->pd_half_std[l]);
                            }
                        }
                        pd_leg += df * (cpn->wcst + cpn->wspot * fwd + opt_string);
                    }
                    else
                    {
                        /*	Floor */
                        floor = OPT_VAL_MACRO(
                            eval_const->pd_floor_type[l],
                            fwd,
                            eval_const->pd_floor_str[l],
                            eval_const->pd_std[l],
                            eval_const->pd_half_std[l]);

                        /*	Cap */
                        cap = OPT_VAL_MACRO(
                            eval_const->pd_cap_type[l],
                            fwd,
                            eval_const->pd_cap_str[l],
                            eval_const->pd_std[l],
                            eval_const->pd_half_std[l]);

                        /*	Coupon pv */

                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                        eval_const->pd_abs_beta[l] * (floor - cap));
                    }
                }
                /*	Final notional */

                if (eval_redemption)
                {
                    /*	Coupon access */
                    cpn = &(cpd->pd_leg->not_ref);

                    /*	Discount */
                    /*df =
                            eval_const->fin_not_disc_dff
                            * exp (- eval_const->fin_not_disc_gam * svijk[0] - eval_const->fin_not_disc_gam2); */
                    df = fin_not_disc_K;

                    /*	Fwd fx */
                    /* fwd = fx
                            * eval_const->fin_not_fwd_dff_ratio
                            * exp (	-	eval_const->fin_not_fwd_for_gam * svijk[1]
                                            +	eval_const->fin_not_fwd_dom_gam * svijk[0]
                                            +	eval_const->fin_not_fwd_cvx); */
                    fwd = fx * fin_not_fwd_K;

                    // interp redemption case:
                    if (cpn->use_opt_str)
                    {
                        opt_string = 0.0;
                        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++)
                        {
                            if (fabs(cpn->weights[str_idx]) > 1e-16)
                            {
                                if (eval_const->fin_not_std > 1e-16 &&
                                    cpn->strikes[str_idx] > 1e-16)
                                    opt_type = 3;
                                else
                                    opt_type = 1;

                                opt_string +=
                                    cpn->weights[str_idx] * OPT_VAL_MACRO(
                                                                opt_type,
                                                                fwd,
                                                                cpn->strikes[str_idx],
                                                                eval_const->fin_not_std,
                                                                eval_const->fin_not_half_std);
                            }
                        }
                        pd_leg += df * (cpn->wcst + cpn->wspot * fwd + opt_string);
                    }
                    else
                    {
                        /*	Floor */
                        floor = OPT_VAL_MACRO(
                            eval_const->fin_not_floor_type,
                            fwd,
                            eval_const->fin_not_floor_str,
                            eval_const->fin_not_std,
                            eval_const->fin_not_half_std);

                        /*	Cap */
                        cap = OPT_VAL_MACRO(
                            eval_const->fin_not_cap_type,
                            fwd,
                            eval_const->fin_not_cap_str,
                            eval_const->fin_not_std,
                            eval_const->fin_not_half_std);

                        /*	Coupon pv */

                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                        eval_const->fin_not_abs_beta * (floor - cap));
                    }
                }
                else
                {
                    /* pd_leg += eval_const->next_pd_dff
                     * exp (- eval_const->next_pd_gam * svijk[0] - eval_const->next_pd_gam2)
                     * next_call->pd_not_amt; */

                    df = next_pd_K;
                    pd_leg += df * next_call->pd_not_amt;
                }

                /*	Initial notional */

                /* fee = eval_const->in_not_pd_dff
                 * exp (- eval_const->in_not_pd_gam * svijk[0] - eval_const->in_not_pd_gam2); */
                df = in_not_pd_K;

                if (call->use_opt_str)
                {
                    fwd = fx * in_not_pd_fwd_K;

                    opt_string = 0.0;
                    for (str_idx = 0; str_idx < call->nstrikes; str_idx++)
                    {
                        if (fabs(call->weights[str_idx]) > 1e-16)
                        {
                            if (eval_const->in_not_pd_std > 1e-16 && call->strikes[str_idx] > 1e-16)
                                opt_type = 3;
                            else
                                opt_type = 1;

                            opt_string +=
                                call->weights[str_idx] * OPT_VAL_MACRO(
                                                             opt_type,
                                                             fwd,
                                                             call->strikes[str_idx],
                                                             eval_const->in_not_pd_std,
                                                             eval_const->in_not_pd_half_std);
                        }
                    }
                    fee = call->wcst + call->wspot * fwd + opt_string - call->pd_not_amt;
                    fee = df * (call->fee + (call->pay_rec == 0 ? fee : -fee));
                }
                else
                {
                    fee = df * call->fee;
                }

                pd_leg -= df * call->pd_not_amt;

                /*	Finally, intrinsic value */

                if (call->pay_rec == 0)
                {
                    iv = pd_leg - fund_leg;
                }
                else
                {
                    iv = fund_leg - pd_leg;
                }

                /*	Process max */

                if (prod_val[i][j][k][0] - iv + fee > 0)
                {
                    prod_val[i][j][k][0] -= iv;
                }
                else
                {
                    /* we pay the fee and call */
                    prod_val[i][j][k][0] = -fee;
                }

                /* Get the IV from the tree */
                if (nprod == 2)
                {
                    prod_val[i][j][k][1] += iv;
                }

                /* For Optimisation */
                fx_K *= fx_DK;
                start_K *= start_DK;
                in_not_fund_K *= in_not_fund_DK;
                in_not_pd_K *= in_not_pd_DK;
                if (call->use_opt_str)
                {
                    in_not_pd_fwd_K *= in_not_pd_fwd_DK;
                }

                if (!eval_redemption)
                {
                    next_start_fund_K *= next_start_fund_DK;
                    next_fund_K *= next_fund_DK;
                    next_pd_K *= next_pd_DK;
                }
                else
                {
                    fin_not_disc_K *= fin_not_disc_DK;
                    fin_not_fwd_K *= fin_not_fwd_DK;
                }

                for (l = 0; l < num_fund_cpn; l++)
                {
                    eval_const->fund_K[l] *= eval_const->fund_DK[l];
                }

                for (l = 0; l < num_pd_cpn; l++)
                {
                    eval_const->pd_disc_K[l] *= eval_const->pd_disc_DK[l];

                    eval_const->pd_fwd_K[l] *= eval_const->pd_fwd_DK[l];
                }
            }

            /* For Optimisation */
            fx_J *= fx_DJ;
            fx_K = fx_J;

            start_J *= start_DJ;
            start_K = start_J;

            in_not_fund_J *= in_not_fund_DJ;
            in_not_fund_K = in_not_fund_J;

            in_not_pd_J *= in_not_pd_DJ;
            in_not_pd_K = in_not_pd_J;

            if (call->use_opt_str)
            {
                in_not_pd_fwd_J *= in_not_pd_fwd_DJ;
                in_not_pd_fwd_K = in_not_pd_fwd_J;
            }

            if (!eval_redemption)
            {
                next_start_fund_J *= next_start_fund_DJ;
                next_start_fund_K = next_start_fund_J;

                next_fund_J *= next_fund_DJ;
                next_fund_K = next_fund_J;

                next_pd_J *= next_pd_DJ;
                next_pd_K = next_pd_J;
            }
            else
            {
                fin_not_disc_J *= fin_not_disc_DJ;
                fin_not_disc_K = fin_not_disc_J;

                fin_not_fwd_J *= fin_not_fwd_DJ;
                fin_not_fwd_K = fin_not_fwd_J;
            }

            for (l = 0; l < num_fund_cpn; l++)
            {
                eval_const->fund_J[l] *= eval_const->fund_DJ[l];
                eval_const->fund_K[l] = eval_const->fund_J[l];
            }

            for (l = 0; l < num_pd_cpn; l++)
            {
                eval_const->pd_disc_J[l] *= eval_const->pd_disc_DJ[l];
                eval_const->pd_disc_K[l] = eval_const->pd_disc_J[l];

                eval_const->pd_fwd_J[l] *= eval_const->pd_fwd_DJ[l];
                eval_const->pd_fwd_K[l] = eval_const->pd_fwd_J[l];
            }
        }

        /* For Optimisation */
        fx_I *= fx_DI;
        fx_J = fx_I;
        fx_K = fx_I;

        start_I *= start_DI;
        start_J = start_I;
        start_K = start_I;

        in_not_fund_I *= in_not_fund_DI;
        in_not_fund_J = in_not_fund_I;
        in_not_fund_K = in_not_fund_I;

        in_not_pd_I *= in_not_pd_DI;
        in_not_pd_J = in_not_pd_I;
        in_not_pd_K = in_not_pd_I;

        if (call->use_opt_str)
        {
            in_not_pd_fwd_I *= in_not_pd_fwd_DI;
            in_not_pd_fwd_J = in_not_pd_fwd_I;
            in_not_pd_fwd_K = in_not_pd_fwd_I;
        }

        if (!eval_redemption)
        {
            next_start_fund_I *= next_start_fund_DI;
            next_start_fund_J = next_start_fund_I;
            next_start_fund_K = next_start_fund_I;

            next_fund_I *= next_fund_DI;
            next_fund_J = next_fund_I;
            next_fund_K = next_fund_I;

            next_pd_I *= next_pd_DI;
            next_pd_J = next_pd_I;
            next_pd_K = next_pd_I;
        }
        else
        {
            fin_not_disc_I *= fin_not_disc_DI;
            fin_not_disc_J = fin_not_disc_I;
            fin_not_disc_K = fin_not_disc_I;

            fin_not_fwd_I *= fin_not_fwd_DI;
            fin_not_fwd_J = fin_not_fwd_I;
            fin_not_fwd_K = fin_not_fwd_I;
        }

        for (l = 0; l < num_fund_cpn; l++)
        {
            eval_const->fund_I[l] *= eval_const->fund_DI[l];
            eval_const->fund_J[l] = eval_const->fund_I[l];
            eval_const->fund_K[l] = eval_const->fund_I[l];
        }

        for (l = 0; l < num_pd_cpn; l++)
        {
            eval_const->pd_disc_I[l] *= eval_const->pd_disc_DI[l];
            eval_const->pd_disc_J[l] = eval_const->pd_disc_I[l];
            eval_const->pd_disc_K[l] = eval_const->pd_disc_I[l];

            eval_const->pd_fwd_I[l] *= eval_const->pd_fwd_DI[l];
            eval_const->pd_fwd_J[l] = eval_const->pd_fwd_I[l];
            eval_const->pd_fwd_K[l] = eval_const->pd_fwd_I[l];
        }
    }

    /*	End of payoff valuation
            ----------------------- */

    return err;
}

Err cpd_payoff_4_3dfx_tree_ko(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    /* Nodes data */
    long n1,
    long n2,
    long n3,
    /* i: d1, j: d2, k: d3, l = {0: xDom, 1: xFor, 2: log (Fx/Fx0)} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call, next_call;
    CPD_UND        und;
    PD_EXO_CPN     cpn;
    CPD_EVAL_CONST eval_const;

    int call_idx, start_idx;

    int i, j, k, l;

    double ***svi, **svij, *svijk, fx;

    double fund_leg, df, fwd, floor, cap, pd_leg, iv, fee, jump;

    int num_fund_cpn, num_pd_cpn, eval_redemption;

    Err err = NULL;

    /*	Get the event */
    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    /*	Calculate number of coupons to be computed */
    if (call_idx == cpd->num_calls - 1)
    {
        /*	Last call: eval all + redemption */
        num_fund_cpn    = call->num_fund_cpn;
        num_pd_cpn      = call->num_pd_cpn;
        eval_redemption = 1;
    }
    else
    {
        /*	Not the last call: only eval coupons up to next call date */
        next_call       = cpd->call + (call_idx + 1);
        num_fund_cpn    = next_call->fund_idx - call->fund_idx;
        num_pd_cpn      = next_call->pd_idx - call->pd_idx;
        eval_redemption = 0;
    }

    /*	Eval payoff
            ----------- */

    if (cpd->call[call_idx].call_type == 0)
    {
        /*	Callable Case	*/
        for (i = 0; i < n1; i++)
        {
            svi = sv[i];
            for (j = 0; j < n2; j++)
            {
                svij = svi[j];
                for (k = 0; k < n3; k++)
                {
                    svijk = svij[k];

                    fx = spot_fx * exp(svijk[2]);

                    /*	PV of funding leg */

                    fund_leg = 0.0;

                    /*	Libor */

                    fund_leg += eval_const->start_dff *
                                exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                    eval_const->start_gam2) *
                                cpd->fund_leg->notional;

                    /*	Coupons */

                    for (l = 0; l < num_fund_cpn; l++)
                    {
                        fund_leg += eval_const->fund_dff[l] *
                                    exp(-eval_const->fund_gam[l] * svijk[eval_const->fund_var] -
                                        eval_const->fund_gam2[l]) *
                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                    }

                    /*	Initial notional */

                    fund_leg -= eval_const->in_not_fund_dff *
                                exp(-eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                                    eval_const->in_not_fund_gam2) *
                                call->fund_not_amt;

                    /*	Next exchange if relevant */
                    if (!eval_redemption)
                    {
                        fund_leg -=
                            eval_const->next_start_fund_dff *
                            exp(-eval_const->next_start_fund_gam * svijk[eval_const->fund_var] -
                                eval_const->next_start_fund_gam2) *
                            cpd->fund_leg->notional;

                        fund_leg += eval_const->next_fund_dff *
                                    exp(-eval_const->next_fund_gam * svijk[eval_const->fund_var] -
                                        eval_const->next_fund_gam2) *
                                    next_call->fund_not_amt;
                    }

                    /*	Fx */

                    if (cpd->fund_leg->dom_for != 0)
                    {
                        fund_leg *= fx;
                    }

                    /*	PV of pd leg */

                    pd_leg = 0.0;

                    /*	Coupons */
                    for (l = 0; l < num_pd_cpn; l++)
                    {
                        /*	Coupon access */
                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                        /*	Discount */
                        df = eval_const->pd_disc_dff[l] *
                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                 eval_const->pd_disc_gam2[l]);

                        /*	Fwd fx */
                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                  eval_const->pd_fwd_cvx[l]);

                        /*	Floor */
                        floor = OPT_VAL_MACRO(
                            eval_const->pd_floor_type[l],
                            fwd,
                            eval_const->pd_floor_str[l],
                            eval_const->pd_std[l],
                            eval_const->pd_half_std[l]);

                        /*	Cap */
                        cap = OPT_VAL_MACRO(
                            eval_const->pd_cap_type[l],
                            fwd,
                            eval_const->pd_cap_str[l],
                            eval_const->pd_std[l],
                            eval_const->pd_half_std[l]);

                        /*	Coupon pv */

                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                        eval_const->pd_abs_beta[l] * (floor - cap));
                    }
                    /*	Final notional */

                    if (eval_redemption)
                    {
                        /*	Coupon access */
                        cpn = &(cpd->pd_leg->not_ref);

                        /*	Discount */
                        df = eval_const->fin_not_disc_dff *
                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                 eval_const->fin_not_disc_gam2);

                        /*	Fwd fx */
                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                  eval_const->fin_not_fwd_cvx);

                        /*	Floor */
                        floor = OPT_VAL_MACRO(
                            eval_const->fin_not_floor_type,
                            fwd,
                            eval_const->fin_not_floor_str,
                            eval_const->fin_not_std,
                            eval_const->fin_not_half_std);

                        /*	Cap */
                        cap = OPT_VAL_MACRO(
                            eval_const->fin_not_cap_type,
                            fwd,
                            eval_const->fin_not_cap_str,
                            eval_const->fin_not_std,
                            eval_const->fin_not_half_std);

                        /*	Coupon pv */

                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                        eval_const->fin_not_abs_beta * (floor - cap));
                    }
                    else
                    {
                        pd_leg +=
                            eval_const->next_pd_dff *
                            exp(-eval_const->next_pd_gam * svijk[0] - eval_const->next_pd_gam2) *
                            next_call->pd_not_amt;
                    }

                    /*	Initial notional */
                    fee = eval_const->in_not_pd_dff *
                          exp(-eval_const->in_not_pd_gam * svijk[0] - eval_const->in_not_pd_gam2);

                    pd_leg -= fee * call->pd_not_amt;

                    fee *= call->fee;

                    /*	Finally, intrinsic value */

                    if (call->pay_rec == 0)
                    {
                        iv = pd_leg - fund_leg;
                    }
                    else
                    {
                        iv = fund_leg - pd_leg;
                    }

                    /*	Process max */

                    if (prod_val[i][j][k][0] - iv + fee > 0)
                    {
                        prod_val[i][j][k][0] -= iv;
                    }
                    else
                    {
                        /* we pay the fee and call */
                        prod_val[i][j][k][0] = -fee;
                    }
                }
            }
        }
    }
    else
    {
        /* Knock Out case	*/

        if (is_bar == 1)
        {
            if (evt_time == 0)
            {
                for (i = 0; i < n1; i++)
                {
                    svi = sv[i];
                    for (j = 0; j < n2; j++)
                    {
                        svij = svi[j];
                        for (k = 0; k < n3; k++)
                        {
                            svijk = svij[k];

                            if ((call->bar_type == 0 && svijk[2] > bar_lvl) ||
                                (call->bar_type == 1 && svijk[2] < bar_lvl))
                            {
                                prod_val[i][j][k][0] = 0.0;
                            }
                        }
                    }
                }
            }
            else
            {
                switch (bar_k)
                {
                case 0:

                    for (j = 0; j < n2; j++)
                    {
                        for (k = 0; k < n3; k++)
                        {
                            /*	First calculate the Jump	*/
                            i     = bar_idx[j][k];
                            svijk = sv[i][j][k];

                            if (i >= 0)
                            {
                                /*	Evaluation of the IV when Fx = Barrier	*/

                                fx = spot_fx * exp(bar_lvl);

                                /*	PV of funding leg */

                                fund_leg = 0.0;

                                /*	Libor */

                                fund_leg +=
                                    eval_const->start_dff *
                                    exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                        eval_const->start_gam2) *
                                    cpd->fund_leg->notional;

                                /*	Coupons */

                                for (l = 0; l < num_fund_cpn; l++)
                                {
                                    fund_leg +=
                                        eval_const->fund_dff[l] *
                                        exp(-eval_const->fund_gam[l] * svijk[eval_const->fund_var] -
                                            eval_const->fund_gam2[l]) *
                                        cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                }

                                /*	Initial notional */

                                fund_leg -=
                                    eval_const->in_not_fund_dff *
                                    exp(-eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                                        eval_const->in_not_fund_gam2) *
                                    call->fund_not_amt;

                                /*	Next exchange if relevant */
                                if (!eval_redemption)
                                {
                                    fund_leg -= eval_const->next_start_fund_dff *
                                                exp(-eval_const->next_start_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->next_start_fund_gam2) *
                                                cpd->fund_leg->notional;

                                    fund_leg += eval_const->next_fund_dff *
                                                exp(-eval_const->next_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->next_fund_gam2) *
                                                next_call->fund_not_amt;
                                }

                                if (cpd->fund_leg->dom_for != 0)
                                {
                                    fund_leg *= fx;
                                }

                                /*	PV of pd leg */

                                pd_leg = 0.0;

                                /*	Coupons */
                                for (l = 0; l < num_pd_cpn; l++)
                                {
                                    /*	Coupon access */
                                    cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                    /*	Discount */
                                    df = eval_const->pd_disc_dff[l] *
                                         exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                             eval_const->pd_disc_gam2[l]);

                                    /*	Fwd fx */
                                    fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                          exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                              eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                              eval_const->pd_fwd_cvx[l]);

                                    /*	Floor */
                                    floor = OPT_VAL_MACRO(
                                        eval_const->pd_floor_type[l],
                                        fwd,
                                        eval_const->pd_floor_str[l],
                                        eval_const->pd_std[l],
                                        eval_const->pd_half_std[l]);

                                    /*	Cap */
                                    cap = OPT_VAL_MACRO(
                                        eval_const->pd_cap_type[l],
                                        fwd,
                                        eval_const->pd_cap_str[l],
                                        eval_const->pd_std[l],
                                        eval_const->pd_half_std[l]);

                                    /*	Coupon pv */

                                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                    eval_const->pd_abs_beta[l] * (floor - cap));
                                }
                                /*	Final notional */

                                if (eval_redemption)
                                {
                                    /*	Coupon access */
                                    cpn = &(cpd->pd_leg->not_ref);

                                    /*	Discount */
                                    df = eval_const->fin_not_disc_dff *
                                         exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                             eval_const->fin_not_disc_gam2);

                                    /*	Fwd fx */
                                    fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                          exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                              eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                              eval_const->fin_not_fwd_cvx);

                                    /*	Floor */
                                    floor = OPT_VAL_MACRO(
                                        eval_const->fin_not_floor_type,
                                        fwd,
                                        eval_const->fin_not_floor_str,
                                        eval_const->fin_not_std,
                                        eval_const->fin_not_half_std);

                                    /*	Cap */
                                    cap = OPT_VAL_MACRO(
                                        eval_const->fin_not_cap_type,
                                        fwd,
                                        eval_const->fin_not_cap_str,
                                        eval_const->fin_not_std,
                                        eval_const->fin_not_half_std);

                                    /*	Coupon pv */

                                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                    eval_const->fin_not_abs_beta * (floor - cap));
                                }
                                else
                                {
                                    pd_leg += eval_const->next_pd_dff *
                                              exp(-eval_const->next_pd_gam * svijk[0] -
                                                  eval_const->next_pd_gam2) *
                                              next_call->pd_not_amt;
                                }

                                /*	Initial notional */

                                fee = eval_const->in_not_pd_dff *
                                      exp(-eval_const->in_not_pd_gam * svijk[0] -
                                          eval_const->in_not_pd_gam2);

                                pd_leg -= fee * call->pd_not_amt;
                                fee *= call->fee;

                                /*	Finally, intrinsic value */

                                if (call->pay_rec == 0)
                                {
                                    iv = pd_leg - fund_leg;
                                }
                                else
                                {
                                    iv = fund_leg - pd_leg;
                                }

                                jump = -(prod_val[i][j][k][0] - iv + fee);

                                /*
                                iv = fx;
                                jump = -fx;
                                */

                                if (call->bar_type == 1)
                                {
                                    jump *= -1.0;
                                }

                                for (l = 0; l < n1; l++)
                                {
                                    prod_val[l][j][k][nprod] = jump;
                                }
                            }
                            else
                            {
                                jump = 0.0;
                            }

                            start_idx = max(bar_idx[j][k], 0);

                            /*	First Part: before the barrier	*/
                            for (i = 0; i < start_idx; i++)
                            {
                                if (call->bar_type == 0)
                                {
                                    /* We get pv - iv */

                                    svijk = sv[i][j][k];

                                    fx = spot_fx * exp(svijk[2]);

                                    /*	PV of funding leg */

                                    fund_leg = 0.0;

                                    /*	Libor */

                                    fund_leg +=
                                        eval_const->start_dff *
                                        exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                            eval_const->start_gam2) *
                                        cpd->fund_leg->notional;

                                    /*	Coupons */

                                    for (l = 0; l < num_fund_cpn; l++)
                                    {
                                        fund_leg += eval_const->fund_dff[l] *
                                                    exp(-eval_const->fund_gam[l] *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->fund_gam2[l]) *
                                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                    }

                                    /*	Initial notional */

                                    fund_leg -= eval_const->in_not_fund_dff *
                                                exp(-eval_const->in_not_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->in_not_fund_gam2) *
                                                call->fund_not_amt;

                                    /*	Next exchange if relevant */
                                    if (!eval_redemption)
                                    {
                                        fund_leg -= eval_const->next_start_fund_dff *
                                                    exp(-eval_const->next_start_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_start_fund_gam2) *
                                                    cpd->fund_leg->notional;

                                        fund_leg += eval_const->next_fund_dff *
                                                    exp(-eval_const->next_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_fund_gam2) *
                                                    next_call->fund_not_amt;
                                    }

                                    if (cpd->fund_leg->dom_for != 0)
                                    {
                                        fund_leg *= fx;
                                    }

                                    /*	PV of pd leg */

                                    pd_leg = 0.0;

                                    /*	Coupons */
                                    for (l = 0; l < num_pd_cpn; l++)
                                    {
                                        /*	Coupon access */
                                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                        /*	Discount */
                                        df = eval_const->pd_disc_dff[l] *
                                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                                 eval_const->pd_disc_gam2[l]);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                                  eval_const->pd_fwd_cvx[l]);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->pd_floor_type[l],
                                            fwd,
                                            eval_const->pd_floor_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->pd_cap_type[l],
                                            fwd,
                                            eval_const->pd_cap_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Coupon pv */

                                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                        eval_const->pd_abs_beta[l] * (floor - cap));
                                    }
                                    /*	Final notional */

                                    if (eval_redemption)
                                    {
                                        /*	Coupon access */
                                        cpn = &(cpd->pd_leg->not_ref);

                                        /*	Discount */
                                        df = eval_const->fin_not_disc_dff *
                                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                                 eval_const->fin_not_disc_gam2);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                                  eval_const->fin_not_fwd_cvx);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->fin_not_floor_type,
                                            fwd,
                                            eval_const->fin_not_floor_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->fin_not_cap_type,
                                            fwd,
                                            eval_const->fin_not_cap_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Coupon pv */

                                        pd_leg +=
                                            df * (cpn->alpha + cpn->beta * fwd +
                                                  eval_const->fin_not_abs_beta * (floor - cap));
                                    }
                                    else
                                    {
                                        pd_leg += eval_const->next_pd_dff *
                                                  exp(-eval_const->next_pd_gam * svijk[0] -
                                                      eval_const->next_pd_gam2) *
                                                  next_call->pd_not_amt;
                                    }

                                    /*	Initial notional */

                                    fee = eval_const->in_not_pd_dff *
                                          exp(-eval_const->in_not_pd_gam * svijk[0] -
                                              eval_const->in_not_pd_gam2);

                                    pd_leg -= fee * call->pd_not_amt;
                                    fee *= call->fee;

                                    /*	Finally, intrinsic value */

                                    if (call->pay_rec == 0)
                                    {
                                        iv = pd_leg - fund_leg;
                                    }
                                    else
                                    {
                                        iv = fund_leg - pd_leg;
                                    }

                                    /*
                                    iv = -fx;
                                    */

                                    /*	Process max: we are under the barrier */
                                    prod_val[i][j][k][0] -= iv;
                                }
                                else
                                {
                                    /* We get 0 */
                                    prod_val[i][j][k][0] = 0.0;
                                }
                            }

                            /*	Second part: after the barrierand we remove the jump	*/

                            for (i = start_idx; i < n1; i++)
                            {
                                if (call->bar_type == 0)
                                {
                                    /* We get 0 and remove the jump */
                                    prod_val[i][j][k][0] = -jump;
                                }
                                else
                                {
                                    /* We get pv - iv and remove the jump */

                                    svijk = sv[i][j][k];

                                    fx = spot_fx * exp(svijk[2]);

                                    /*	PV of funding leg */

                                    fund_leg = 0.0;

                                    /*	Libor */

                                    fund_leg +=
                                        eval_const->start_dff *
                                        exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                            eval_const->start_gam2) *
                                        cpd->fund_leg->notional;

                                    /*	Coupons */

                                    for (l = 0; l < num_fund_cpn; l++)
                                    {
                                        fund_leg += eval_const->fund_dff[l] *
                                                    exp(-eval_const->fund_gam[l] *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->fund_gam2[l]) *
                                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                    }

                                    /*	Initial notional */

                                    fund_leg -= eval_const->in_not_fund_dff *
                                                exp(-eval_const->in_not_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->in_not_fund_gam2) *
                                                call->fund_not_amt;

                                    /*	Next exchange if relevant */
                                    if (!eval_redemption)
                                    {
                                        fund_leg -= eval_const->next_start_fund_dff *
                                                    exp(-eval_const->next_start_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_start_fund_gam2) *
                                                    cpd->fund_leg->notional;

                                        fund_leg += eval_const->next_fund_dff *
                                                    exp(-eval_const->next_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_fund_gam2) *
                                                    next_call->fund_not_amt;
                                    }

                                    if (cpd->fund_leg->dom_for != 0)
                                    {
                                        fund_leg *= fx;
                                    }

                                    /*	PV of pd leg */

                                    pd_leg = 0.0;

                                    /*	Coupons */
                                    for (l = 0; l < num_pd_cpn; l++)
                                    {
                                        /*	Coupon access */
                                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                        /*	Discount */
                                        df = eval_const->pd_disc_dff[l] *
                                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                                 eval_const->pd_disc_gam2[l]);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                                  eval_const->pd_fwd_cvx[l]);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->pd_floor_type[l],
                                            fwd,
                                            eval_const->pd_floor_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->pd_cap_type[l],
                                            fwd,
                                            eval_const->pd_cap_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Coupon pv */

                                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                        eval_const->pd_abs_beta[l] * (floor - cap));
                                    }
                                    /*	Final notional */

                                    if (eval_redemption)
                                    {
                                        /*	Coupon access */
                                        cpn = &(cpd->pd_leg->not_ref);

                                        /*	Discount */
                                        df = eval_const->fin_not_disc_dff *
                                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                                 eval_const->fin_not_disc_gam2);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                                  eval_const->fin_not_fwd_cvx);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->fin_not_floor_type,
                                            fwd,
                                            eval_const->fin_not_floor_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->fin_not_cap_type,
                                            fwd,
                                            eval_const->fin_not_cap_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Coupon pv */

                                        pd_leg +=
                                            df * (cpn->alpha + cpn->beta * fwd +
                                                  eval_const->fin_not_abs_beta * (floor - cap));
                                    }
                                    else
                                    {
                                        pd_leg += eval_const->next_pd_dff *
                                                  exp(-eval_const->next_pd_gam * svijk[0] -
                                                      eval_const->next_pd_gam2) *
                                                  next_call->pd_not_amt;
                                    }

                                    /*	Initial notional */

                                    fee = eval_const->in_not_pd_dff *
                                          exp(-eval_const->in_not_pd_gam * svijk[0] -
                                              eval_const->in_not_pd_gam2);

                                    pd_leg -= fee * call->pd_not_amt;
                                    fee *= call->fee;

                                    /*	Finally, intrinsic value */

                                    if (call->pay_rec == 0)
                                    {
                                        iv = pd_leg - fund_leg;
                                    }
                                    else
                                    {
                                        iv = fund_leg - pd_leg;
                                    }

                                    /*
                                    iv = -fx;
                                    */

                                    /*	Process max: we are under the barrier */
                                    prod_val[i][j][k][0] -= iv + jump;
                                }
                            }
                        }
                    }
                    break;

                case 1:

                    for (i = 0; i < n1; i++)
                    {
                        for (k = 0; k < n3; k++)
                        {
                            /*	First calculate the Jump	*/
                            j     = bar_idx[i][k];
                            svijk = sv[i][j][k];

                            if (j >= 0)
                            {
                                /*	Evaluation of the IV when Fx = Barrier	*/

                                fx = spot_fx * exp(bar_lvl);

                                /*	PV of funding leg */

                                fund_leg = 0.0;

                                /*	Libor */

                                fund_leg +=
                                    eval_const->start_dff *
                                    exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                        eval_const->start_gam2) *
                                    cpd->fund_leg->notional;

                                /*	Coupons */

                                for (l = 0; l < num_fund_cpn; l++)
                                {
                                    fund_leg +=
                                        eval_const->fund_dff[l] *
                                        exp(-eval_const->fund_gam[l] * svijk[eval_const->fund_var] -
                                            eval_const->fund_gam2[l]) *
                                        cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                }

                                /*	Initial notional */

                                fund_leg -=
                                    eval_const->in_not_fund_dff *
                                    exp(-eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                                        eval_const->in_not_fund_gam2) *
                                    call->fund_not_amt;

                                /*	Next exchange if relevant */
                                if (!eval_redemption)
                                {
                                    fund_leg -= eval_const->next_start_fund_dff *
                                                exp(-eval_const->next_start_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->next_start_fund_gam2) *
                                                cpd->fund_leg->notional;

                                    fund_leg += eval_const->next_fund_dff *
                                                exp(-eval_const->next_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->next_fund_gam2) *
                                                next_call->fund_not_amt;
                                }

                                if (cpd->fund_leg->dom_for != 0)
                                {
                                    fund_leg *= fx;
                                }

                                /*	PV of pd leg */

                                pd_leg = 0.0;

                                /*	Coupons */
                                for (l = 0; l < num_pd_cpn; l++)
                                {
                                    /*	Coupon access */
                                    cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                    /*	Discount */
                                    df = eval_const->pd_disc_dff[l] *
                                         exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                             eval_const->pd_disc_gam2[l]);

                                    /*	Fwd fx */
                                    fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                          exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                              eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                              eval_const->pd_fwd_cvx[l]);

                                    /*	Floor */
                                    floor = OPT_VAL_MACRO(
                                        eval_const->pd_floor_type[l],
                                        fwd,
                                        eval_const->pd_floor_str[l],
                                        eval_const->pd_std[l],
                                        eval_const->pd_half_std[l]);

                                    /*	Cap */
                                    cap = OPT_VAL_MACRO(
                                        eval_const->pd_cap_type[l],
                                        fwd,
                                        eval_const->pd_cap_str[l],
                                        eval_const->pd_std[l],
                                        eval_const->pd_half_std[l]);

                                    /*	Coupon pv */

                                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                    eval_const->pd_abs_beta[l] * (floor - cap));
                                }
                                /*	Final notional */

                                if (eval_redemption)
                                {
                                    /*	Coupon access */
                                    cpn = &(cpd->pd_leg->not_ref);

                                    /*	Discount */
                                    df = eval_const->fin_not_disc_dff *
                                         exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                             eval_const->fin_not_disc_gam2);

                                    /*	Fwd fx */
                                    fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                          exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                              eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                              eval_const->fin_not_fwd_cvx);

                                    /*	Floor */
                                    floor = OPT_VAL_MACRO(
                                        eval_const->fin_not_floor_type,
                                        fwd,
                                        eval_const->fin_not_floor_str,
                                        eval_const->fin_not_std,
                                        eval_const->fin_not_half_std);

                                    /*	Cap */
                                    cap = OPT_VAL_MACRO(
                                        eval_const->fin_not_cap_type,
                                        fwd,
                                        eval_const->fin_not_cap_str,
                                        eval_const->fin_not_std,
                                        eval_const->fin_not_half_std);

                                    /*	Coupon pv */

                                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                    eval_const->fin_not_abs_beta * (floor - cap));
                                }
                                else
                                {
                                    pd_leg += eval_const->next_pd_dff *
                                              exp(-eval_const->next_pd_gam * svijk[0] -
                                                  eval_const->next_pd_gam2) *
                                              next_call->pd_not_amt;
                                }

                                /*	Initial notional */

                                fee = eval_const->in_not_pd_dff *
                                      exp(-eval_const->in_not_pd_gam * svijk[0] -
                                          eval_const->in_not_pd_gam2);

                                pd_leg -= fee * call->pd_not_amt;
                                fee *= call->fee;

                                /*	Finally, intrinsic value */

                                if (call->pay_rec == 0)
                                {
                                    iv = pd_leg - fund_leg;
                                }
                                else
                                {
                                    iv = fund_leg - pd_leg;
                                }

                                jump = -(prod_val[i][j][k][0] - iv + fee);

                                /*
                                iv = fx;
                                jump = -fx;
                                */

                                if (call->bar_type == 1)
                                {
                                    jump *= -1.0;
                                }

                                for (l = 0; l < n2; l++)
                                {
                                    prod_val[i][l][k][nprod] = jump;
                                }
                            }
                            else
                            {
                                jump = 0.0;
                            }

                            start_idx = max(bar_idx[i][k], 0);

                            /*	First Part: before the barrier	*/
                            for (j = 0; j < start_idx; j++)
                            {
                                if (call->bar_type == 0)
                                {
                                    /* We get pv - iv */

                                    svijk = sv[i][j][k];

                                    fx = spot_fx * exp(svijk[2]);

                                    /*	PV of funding leg */

                                    fund_leg = 0.0;

                                    /*	Libor */

                                    fund_leg +=
                                        eval_const->start_dff *
                                        exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                            eval_const->start_gam2) *
                                        cpd->fund_leg->notional;

                                    /*	Coupons */

                                    for (l = 0; l < num_fund_cpn; l++)
                                    {
                                        fund_leg += eval_const->fund_dff[l] *
                                                    exp(-eval_const->fund_gam[l] *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->fund_gam2[l]) *
                                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                    }

                                    /*	Initial notional */

                                    fund_leg -= eval_const->in_not_fund_dff *
                                                exp(-eval_const->in_not_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->in_not_fund_gam2) *
                                                call->fund_not_amt;

                                    /*	Next exchange if relevant */
                                    if (!eval_redemption)
                                    {
                                        fund_leg -= eval_const->next_start_fund_dff *
                                                    exp(-eval_const->next_start_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_start_fund_gam2) *
                                                    cpd->fund_leg->notional;

                                        fund_leg += eval_const->next_fund_dff *
                                                    exp(-eval_const->next_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_fund_gam2) *
                                                    next_call->fund_not_amt;
                                    }

                                    if (cpd->fund_leg->dom_for != 0)
                                    {
                                        fund_leg *= fx;
                                    }

                                    /*	PV of pd leg */

                                    pd_leg = 0.0;

                                    /*	Coupons */
                                    for (l = 0; l < num_pd_cpn; l++)
                                    {
                                        /*	Coupon access */
                                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                        /*	Discount */
                                        df = eval_const->pd_disc_dff[l] *
                                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                                 eval_const->pd_disc_gam2[l]);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                                  eval_const->pd_fwd_cvx[l]);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->pd_floor_type[l],
                                            fwd,
                                            eval_const->pd_floor_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->pd_cap_type[l],
                                            fwd,
                                            eval_const->pd_cap_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Coupon pv */

                                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                        eval_const->pd_abs_beta[l] * (floor - cap));
                                    }
                                    /*	Final notional */

                                    if (eval_redemption)
                                    {
                                        /*	Coupon access */
                                        cpn = &(cpd->pd_leg->not_ref);

                                        /*	Discount */
                                        df = eval_const->fin_not_disc_dff *
                                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                                 eval_const->fin_not_disc_gam2);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                                  eval_const->fin_not_fwd_cvx);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->fin_not_floor_type,
                                            fwd,
                                            eval_const->fin_not_floor_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->fin_not_cap_type,
                                            fwd,
                                            eval_const->fin_not_cap_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Coupon pv */

                                        pd_leg +=
                                            df * (cpn->alpha + cpn->beta * fwd +
                                                  eval_const->fin_not_abs_beta * (floor - cap));
                                    }
                                    else
                                    {
                                        pd_leg += eval_const->next_pd_dff *
                                                  exp(-eval_const->next_pd_gam * svijk[0] -
                                                      eval_const->next_pd_gam2) *
                                                  next_call->pd_not_amt;
                                    }

                                    /*	Initial notional */

                                    pd_leg -= eval_const->in_not_pd_dff *
                                              exp(-eval_const->in_not_pd_gam * svijk[0] -
                                                  eval_const->in_not_pd_gam2) *
                                              (call->pd_not_amt + call->fee);

                                    /*	Finally, intrinsic value */

                                    if (call->pay_rec == 0)
                                    {
                                        iv = pd_leg - fund_leg;
                                    }
                                    else
                                    {
                                        iv = fund_leg - pd_leg;
                                    }

                                    /*
                                    iv = -fx;
                                    */

                                    /*	Process max: we are under the barrier */
                                    prod_val[i][j][k][0] -= iv;
                                }
                                else
                                {
                                    /* We get 0 */
                                    prod_val[i][j][k][0] = 0;
                                }
                            }

                            /*	Second part: after the barrierand we remove the jump	*/

                            for (j = start_idx; j < n2; j++)
                            {
                                if (call->bar_type == 0)
                                {
                                    /* We get 0 and remove the jump */
                                    prod_val[i][j][k][0] = -jump;
                                }
                                else
                                {
                                    /* We get pv - iv and remove the jump */

                                    svijk = sv[i][j][k];

                                    fx = spot_fx * exp(svijk[2]);

                                    /*	PV of funding leg */

                                    fund_leg = 0.0;

                                    /*	Libor */

                                    fund_leg +=
                                        eval_const->start_dff *
                                        exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                            eval_const->start_gam2) *
                                        cpd->fund_leg->notional;

                                    /*	Coupons */

                                    for (l = 0; l < num_fund_cpn; l++)
                                    {
                                        fund_leg += eval_const->fund_dff[l] *
                                                    exp(-eval_const->fund_gam[l] *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->fund_gam2[l]) *
                                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                    }

                                    /*	Initial notional */

                                    fund_leg -= eval_const->in_not_fund_dff *
                                                exp(-eval_const->in_not_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->in_not_fund_gam2) *
                                                call->fund_not_amt;

                                    /*	Next exchange if relevant */
                                    if (!eval_redemption)
                                    {
                                        fund_leg -= eval_const->next_start_fund_dff *
                                                    exp(-eval_const->next_start_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_start_fund_gam2) *
                                                    cpd->fund_leg->notional;

                                        fund_leg += eval_const->next_fund_dff *
                                                    exp(-eval_const->next_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_fund_gam2) *
                                                    next_call->fund_not_amt;
                                    }

                                    if (cpd->fund_leg->dom_for != 0)
                                    {
                                        fund_leg *= fx;
                                    }

                                    /*	PV of pd leg */

                                    pd_leg = 0.0;

                                    /*	Coupons */
                                    for (l = 0; l < num_pd_cpn; l++)
                                    {
                                        /*	Coupon access */
                                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                        /*	Discount */
                                        df = eval_const->pd_disc_dff[l] *
                                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                                 eval_const->pd_disc_gam2[l]);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                                  eval_const->pd_fwd_cvx[l]);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->pd_floor_type[l],
                                            fwd,
                                            eval_const->pd_floor_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->pd_cap_type[l],
                                            fwd,
                                            eval_const->pd_cap_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Coupon pv */

                                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                        eval_const->pd_abs_beta[l] * (floor - cap));
                                    }
                                    /*	Final notional */

                                    if (eval_redemption)
                                    {
                                        /*	Coupon access */
                                        cpn = &(cpd->pd_leg->not_ref);

                                        /*	Discount */
                                        df = eval_const->fin_not_disc_dff *
                                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                                 eval_const->fin_not_disc_gam2);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                                  eval_const->fin_not_fwd_cvx);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->fin_not_floor_type,
                                            fwd,
                                            eval_const->fin_not_floor_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->fin_not_cap_type,
                                            fwd,
                                            eval_const->fin_not_cap_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Coupon pv */

                                        pd_leg +=
                                            df * (cpn->alpha + cpn->beta * fwd +
                                                  eval_const->fin_not_abs_beta * (floor - cap));
                                    }
                                    else
                                    {
                                        pd_leg += eval_const->next_pd_dff *
                                                  exp(-eval_const->next_pd_gam * svijk[0] -
                                                      eval_const->next_pd_gam2) *
                                                  next_call->pd_not_amt;
                                    }

                                    /*	Initial notional */

                                    pd_leg -= eval_const->in_not_pd_dff *
                                              exp(-eval_const->in_not_pd_gam * svijk[0] -
                                                  eval_const->in_not_pd_gam2) *
                                              (call->pd_not_amt + call->fee);

                                    /*	Finally, intrinsic value */

                                    if (call->pay_rec == 0)
                                    {
                                        iv = pd_leg - fund_leg;
                                    }
                                    else
                                    {
                                        iv = fund_leg - pd_leg;
                                    }

                                    /*
                                    iv = -fx;
                                    */

                                    /*	Process max: we are under the barrier */
                                    prod_val[i][j][k][0] -= iv + jump;
                                }
                            }
                        }
                    }
                    break;

                case 2:

                    for (i = 0; i < n1; i++)
                    {
                        for (j = 0; j < n2; j++)
                        {
                            /*	First calculate the Jump	*/
                            k     = bar_idx[i][j];
                            svijk = sv[i][j][k];

                            if (k >= 0)
                            {
                                /*	Evaluation of the IV when Fx = Barrier	*/

                                fx = spot_fx * exp(bar_lvl);

                                /*	PV of funding leg */

                                fund_leg = 0.0;

                                /*	Libor */

                                fund_leg +=
                                    eval_const->start_dff *
                                    exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                        eval_const->start_gam2) *
                                    cpd->fund_leg->notional;

                                /*	Coupons */

                                for (l = 0; l < num_fund_cpn; l++)
                                {
                                    fund_leg +=
                                        eval_const->fund_dff[l] *
                                        exp(-eval_const->fund_gam[l] * svijk[eval_const->fund_var] -
                                            eval_const->fund_gam2[l]) *
                                        cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                }

                                /*	Initial notional */

                                fund_leg -=
                                    eval_const->in_not_fund_dff *
                                    exp(-eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                                        eval_const->in_not_fund_gam2) *
                                    call->fund_not_amt;

                                /*	Next exchange if relevant */
                                if (!eval_redemption)
                                {
                                    fund_leg -= eval_const->next_start_fund_dff *
                                                exp(-eval_const->next_start_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->next_start_fund_gam2) *
                                                cpd->fund_leg->notional;

                                    fund_leg += eval_const->next_fund_dff *
                                                exp(-eval_const->next_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->next_fund_gam2) *
                                                next_call->fund_not_amt;
                                }

                                if (cpd->fund_leg->dom_for != 0)
                                {
                                    fund_leg *= fx;
                                }

                                /*	PV of pd leg */

                                pd_leg = 0.0;

                                /*	Coupons */
                                for (l = 0; l < num_pd_cpn; l++)
                                {
                                    /*	Coupon access */
                                    cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                    /*	Discount */
                                    df = eval_const->pd_disc_dff[l] *
                                         exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                             eval_const->pd_disc_gam2[l]);

                                    /*	Fwd fx */
                                    fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                          exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                              eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                              eval_const->pd_fwd_cvx[l]);

                                    /*	Floor */
                                    floor = OPT_VAL_MACRO(
                                        eval_const->pd_floor_type[l],
                                        fwd,
                                        eval_const->pd_floor_str[l],
                                        eval_const->pd_std[l],
                                        eval_const->pd_half_std[l]);

                                    /*	Cap */
                                    cap = OPT_VAL_MACRO(
                                        eval_const->pd_cap_type[l],
                                        fwd,
                                        eval_const->pd_cap_str[l],
                                        eval_const->pd_std[l],
                                        eval_const->pd_half_std[l]);

                                    /*	Coupon pv */

                                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                    eval_const->pd_abs_beta[l] * (floor - cap));
                                }
                                /*	Final notional */

                                if (eval_redemption)
                                {
                                    /*	Coupon access */
                                    cpn = &(cpd->pd_leg->not_ref);

                                    /*	Discount */
                                    df = eval_const->fin_not_disc_dff *
                                         exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                             eval_const->fin_not_disc_gam2);

                                    /*	Fwd fx */
                                    fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                          exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                              eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                              eval_const->fin_not_fwd_cvx);

                                    /*	Floor */
                                    floor = OPT_VAL_MACRO(
                                        eval_const->fin_not_floor_type,
                                        fwd,
                                        eval_const->fin_not_floor_str,
                                        eval_const->fin_not_std,
                                        eval_const->fin_not_half_std);

                                    /*	Cap */
                                    cap = OPT_VAL_MACRO(
                                        eval_const->fin_not_cap_type,
                                        fwd,
                                        eval_const->fin_not_cap_str,
                                        eval_const->fin_not_std,
                                        eval_const->fin_not_half_std);

                                    /*	Coupon pv */

                                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                    eval_const->fin_not_abs_beta * (floor - cap));
                                }
                                else
                                {
                                    pd_leg += eval_const->next_pd_dff *
                                              exp(-eval_const->next_pd_gam * svijk[0] -
                                                  eval_const->next_pd_gam2) *
                                              next_call->pd_not_amt;
                                }

                                /*	Initial notional */

                                fee = eval_const->in_not_pd_dff *
                                      exp(-eval_const->in_not_pd_gam * svijk[0] -
                                          eval_const->in_not_pd_gam2);

                                pd_leg -= fee * call->pd_not_amt;
                                fee *= call->fee;

                                /*	Finally, intrinsic value */

                                if (call->pay_rec == 0)
                                {
                                    iv = pd_leg - fund_leg;
                                }
                                else
                                {
                                    iv = fund_leg - pd_leg;
                                }

                                jump = -(prod_val[i][j][k][0] - iv + fee);

                                /*
                                iv = fx;
                                jump = -fx;
                                */

                                if (call->bar_type == 1)
                                {
                                    jump *= -1.0;
                                }

                                for (l = 0; l < n3; l++)
                                {
                                    prod_val[i][j][l][nprod] = jump;
                                }
                            }
                            else
                            {
                                jump = 0.0;
                            }

                            start_idx = max(bar_idx[i][j], 0);

                            /*	First Part: before the barrier	*/
                            for (k = 0; k < start_idx; k++)
                            {
                                if (call->bar_type == 0)
                                {
                                    /* We get pv - iv */

                                    svijk = sv[i][j][k];

                                    fx = spot_fx * exp(svijk[2]);

                                    /*	PV of funding leg */

                                    fund_leg = 0.0;

                                    /*	Libor */

                                    fund_leg +=
                                        eval_const->start_dff *
                                        exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                            eval_const->start_gam2) *
                                        cpd->fund_leg->notional;

                                    /*	Coupons */

                                    for (l = 0; l < num_fund_cpn; l++)
                                    {
                                        fund_leg += eval_const->fund_dff[l] *
                                                    exp(-eval_const->fund_gam[l] *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->fund_gam2[l]) *
                                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                    }

                                    /*	Initial notional */

                                    fund_leg -= eval_const->in_not_fund_dff *
                                                exp(-eval_const->in_not_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->in_not_fund_gam2) *
                                                call->fund_not_amt;

                                    /*	Next exchange if relevant */
                                    if (!eval_redemption)
                                    {
                                        fund_leg -= eval_const->next_start_fund_dff *
                                                    exp(-eval_const->next_start_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_start_fund_gam2) *
                                                    cpd->fund_leg->notional;

                                        fund_leg += eval_const->next_fund_dff *
                                                    exp(-eval_const->next_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_fund_gam2) *
                                                    next_call->fund_not_amt;
                                    }

                                    if (cpd->fund_leg->dom_for != 0)
                                    {
                                        fund_leg *= fx;
                                    }

                                    /*	PV of pd leg */

                                    pd_leg = 0.0;

                                    /*	Coupons */
                                    for (l = 0; l < num_pd_cpn; l++)
                                    {
                                        /*	Coupon access */
                                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                        /*	Discount */
                                        df = eval_const->pd_disc_dff[l] *
                                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                                 eval_const->pd_disc_gam2[l]);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                                  eval_const->pd_fwd_cvx[l]);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->pd_floor_type[l],
                                            fwd,
                                            eval_const->pd_floor_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->pd_cap_type[l],
                                            fwd,
                                            eval_const->pd_cap_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Coupon pv */

                                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                        eval_const->pd_abs_beta[l] * (floor - cap));
                                    }
                                    /*	Final notional */

                                    if (eval_redemption)
                                    {
                                        /*	Coupon access */
                                        cpn = &(cpd->pd_leg->not_ref);

                                        /*	Discount */
                                        df = eval_const->fin_not_disc_dff *
                                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                                 eval_const->fin_not_disc_gam2);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                                  eval_const->fin_not_fwd_cvx);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->fin_not_floor_type,
                                            fwd,
                                            eval_const->fin_not_floor_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->fin_not_cap_type,
                                            fwd,
                                            eval_const->fin_not_cap_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Coupon pv */

                                        pd_leg +=
                                            df * (cpn->alpha + cpn->beta * fwd +
                                                  eval_const->fin_not_abs_beta * (floor - cap));
                                    }
                                    else
                                    {
                                        pd_leg += eval_const->next_pd_dff *
                                                  exp(-eval_const->next_pd_gam * svijk[0] -
                                                      eval_const->next_pd_gam2) *
                                                  next_call->pd_not_amt;
                                    }

                                    /*	Initial notional */

                                    pd_leg -= eval_const->in_not_pd_dff *
                                              exp(-eval_const->in_not_pd_gam * svijk[0] -
                                                  eval_const->in_not_pd_gam2) *
                                              (call->pd_not_amt + call->fee);

                                    /*	Finally, intrinsic value */

                                    if (call->pay_rec == 0)
                                    {
                                        iv = pd_leg - fund_leg;
                                    }
                                    else
                                    {
                                        iv = fund_leg - pd_leg;
                                    }

                                    /*
                                    iv = -fx;
                                    */

                                    /*	Process max: we are under the barrier */
                                    prod_val[i][j][k][0] -= iv;
                                }
                                else
                                {
                                    /* We get 0 */
                                    prod_val[i][j][k][0] = 0;
                                }
                            }

                            /*	Second part: after the barrierand we remove the jump	*/

                            for (k = start_idx; k < n3; k++)
                            {
                                if (call->bar_type == 0)
                                {
                                    /* We get 0 and remove the jump */
                                    prod_val[i][j][k][0] = -jump;
                                }
                                else
                                {
                                    /* We get pv - iv and remove the jump */

                                    svijk = sv[i][j][k];

                                    fx = spot_fx * exp(svijk[2]);

                                    /*	PV of funding leg */

                                    fund_leg = 0.0;

                                    /*	Libor */

                                    fund_leg +=
                                        eval_const->start_dff *
                                        exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                            eval_const->start_gam2) *
                                        cpd->fund_leg->notional;

                                    /*	Coupons */

                                    for (l = 0; l < num_fund_cpn; l++)
                                    {
                                        fund_leg += eval_const->fund_dff[l] *
                                                    exp(-eval_const->fund_gam[l] *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->fund_gam2[l]) *
                                                    cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                                    }

                                    /*	Initial notional */

                                    fund_leg -= eval_const->in_not_fund_dff *
                                                exp(-eval_const->in_not_fund_gam *
                                                        svijk[eval_const->fund_var] -
                                                    eval_const->in_not_fund_gam2) *
                                                call->fund_not_amt;

                                    /*	Next exchange if relevant */
                                    if (!eval_redemption)
                                    {
                                        fund_leg -= eval_const->next_start_fund_dff *
                                                    exp(-eval_const->next_start_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_start_fund_gam2) *
                                                    cpd->fund_leg->notional;

                                        fund_leg += eval_const->next_fund_dff *
                                                    exp(-eval_const->next_fund_gam *
                                                            svijk[eval_const->fund_var] -
                                                        eval_const->next_fund_gam2) *
                                                    next_call->fund_not_amt;
                                    }

                                    if (cpd->fund_leg->dom_for != 0)
                                    {
                                        fund_leg *= fx;
                                    }

                                    /*	PV of pd leg */

                                    pd_leg = 0.0;

                                    /*	Coupons */
                                    for (l = 0; l < num_pd_cpn; l++)
                                    {
                                        /*	Coupon access */
                                        cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                                        /*	Discount */
                                        df = eval_const->pd_disc_dff[l] *
                                             exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                                 eval_const->pd_disc_gam2[l]);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                              exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                                  eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                                  eval_const->pd_fwd_cvx[l]);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->pd_floor_type[l],
                                            fwd,
                                            eval_const->pd_floor_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->pd_cap_type[l],
                                            fwd,
                                            eval_const->pd_cap_str[l],
                                            eval_const->pd_std[l],
                                            eval_const->pd_half_std[l]);

                                        /*	Coupon pv */

                                        pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                                        eval_const->pd_abs_beta[l] * (floor - cap));
                                    }
                                    /*	Final notional */

                                    if (eval_redemption)
                                    {
                                        /*	Coupon access */
                                        cpn = &(cpd->pd_leg->not_ref);

                                        /*	Discount */
                                        df = eval_const->fin_not_disc_dff *
                                             exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                                 eval_const->fin_not_disc_gam2);

                                        /*	Fwd fx */
                                        fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                              exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                                  eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                                  eval_const->fin_not_fwd_cvx);

                                        /*	Floor */
                                        floor = OPT_VAL_MACRO(
                                            eval_const->fin_not_floor_type,
                                            fwd,
                                            eval_const->fin_not_floor_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Cap */
                                        cap = OPT_VAL_MACRO(
                                            eval_const->fin_not_cap_type,
                                            fwd,
                                            eval_const->fin_not_cap_str,
                                            eval_const->fin_not_std,
                                            eval_const->fin_not_half_std);

                                        /*	Coupon pv */

                                        pd_leg +=
                                            df * (cpn->alpha + cpn->beta * fwd +
                                                  eval_const->fin_not_abs_beta * (floor - cap));
                                    }
                                    else
                                    {
                                        pd_leg += eval_const->next_pd_dff *
                                                  exp(-eval_const->next_pd_gam * svijk[0] -
                                                      eval_const->next_pd_gam2) *
                                                  next_call->pd_not_amt;
                                    }

                                    /*	Initial notional */

                                    pd_leg -= eval_const->in_not_pd_dff *
                                              exp(-eval_const->in_not_pd_gam * svijk[0] -
                                                  eval_const->in_not_pd_gam2) *
                                              (call->pd_not_amt + call->fee);

                                    /*	Finally, intrinsic value */

                                    if (call->pay_rec == 0)
                                    {
                                        iv = pd_leg - fund_leg;
                                    }
                                    else
                                    {
                                        iv = fund_leg - pd_leg;
                                    }

                                    /*
                                    iv = -fx;
                                    */

                                    /*	Process max: we are under the barrier */
                                    prod_val[i][j][k][0] -= iv + jump;
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;
                }
            }
        }
        else if ((call->bar_type == 0 && is_bar == -1) || (call->bar_type == 1 && is_bar == -2))
        {
            for (i = 0; i < n1; i++)
            {
                for (j = 0; j < n2; j++)
                {
                    for (k = 0; k < n3; k++)
                    {
                        prod_val[i][j][k][0] = 0.0;
                    }
                }
            }
        }
        else
        {
            /* we pay the coupon */
            for (i = 0; i < n1; i++)
            {
                svi = sv[i];
                for (j = 0; j < n2; j++)
                {
                    svij = svi[j];
                    for (k = 0; k < n3; k++)
                    {
                        svijk = svij[k];

                        fx = spot_fx * exp(svijk[2]);

                        /*	PV of funding leg */

                        fund_leg = 0.0;

                        /*	Libor */

                        fund_leg += eval_const->start_dff *
                                    exp(-eval_const->start_gam * svijk[eval_const->fund_var] -
                                        eval_const->start_gam2) *
                                    cpd->fund_leg->notional;

                        /*	Coupons */

                        for (l = 0; l < num_fund_cpn; l++)
                        {
                            fund_leg += eval_const->fund_dff[l] *
                                        exp(-eval_const->fund_gam[l] * svijk[eval_const->fund_var] -
                                            eval_const->fund_gam2[l]) *
                                        cpd->fund_leg->cpn[call->fund_idx + l].cpn;
                        }

                        /*	Initial notional */

                        fund_leg -= eval_const->in_not_fund_dff *
                                    exp(-eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                                        eval_const->in_not_fund_gam2) *
                                    call->fund_not_amt;

                        /*	Next exchange if relevant */
                        if (!eval_redemption)
                        {
                            fund_leg -=
                                eval_const->next_start_fund_dff *
                                exp(-eval_const->next_start_fund_gam * svijk[eval_const->fund_var] -
                                    eval_const->next_start_fund_gam2) *
                                cpd->fund_leg->notional;

                            fund_leg +=
                                eval_const->next_fund_dff *
                                exp(-eval_const->next_fund_gam * svijk[eval_const->fund_var] -
                                    eval_const->next_fund_gam2) *
                                next_call->fund_not_amt;
                        }

                        /*	Fx */

                        if (cpd->fund_leg->dom_for != 0)
                        {
                            fund_leg *= fx;
                        }

                        /*	PV of pd leg */

                        pd_leg = 0.0;

                        /*	Coupons */
                        for (l = 0; l < num_pd_cpn; l++)
                        {
                            /*	Coupon access */
                            cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                            /*	Discount */
                            df = eval_const->pd_disc_dff[l] *
                                 exp(-eval_const->pd_disc_gam[l] * svijk[0] -
                                     eval_const->pd_disc_gam2[l]);

                            /*	Fwd fx */
                            fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                                  exp(-eval_const->pd_fwd_for_gam[l] * svijk[1] +
                                      eval_const->pd_fwd_dom_gam[l] * svijk[0] +
                                      eval_const->pd_fwd_cvx[l]);

                            /*	Floor */
                            floor = OPT_VAL_MACRO(
                                eval_const->pd_floor_type[l],
                                fwd,
                                eval_const->pd_floor_str[l],
                                eval_const->pd_std[l],
                                eval_const->pd_half_std[l]);

                            /*	Cap */
                            cap = OPT_VAL_MACRO(
                                eval_const->pd_cap_type[l],
                                fwd,
                                eval_const->pd_cap_str[l],
                                eval_const->pd_std[l],
                                eval_const->pd_half_std[l]);

                            /*	Coupon pv */

                            pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                            eval_const->pd_abs_beta[l] * (floor - cap));
                        }
                        /*	Final notional */

                        if (eval_redemption)
                        {
                            /*	Coupon access */
                            cpn = &(cpd->pd_leg->not_ref);

                            /*	Discount */
                            df = eval_const->fin_not_disc_dff *
                                 exp(-eval_const->fin_not_disc_gam * svijk[0] -
                                     eval_const->fin_not_disc_gam2);

                            /*	Fwd fx */
                            fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                                  exp(-eval_const->fin_not_fwd_for_gam * svijk[1] +
                                      eval_const->fin_not_fwd_dom_gam * svijk[0] +
                                      eval_const->fin_not_fwd_cvx);

                            /*	Floor */
                            floor = OPT_VAL_MACRO(
                                eval_const->fin_not_floor_type,
                                fwd,
                                eval_const->fin_not_floor_str,
                                eval_const->fin_not_std,
                                eval_const->fin_not_half_std);

                            /*	Cap */
                            cap = OPT_VAL_MACRO(
                                eval_const->fin_not_cap_type,
                                fwd,
                                eval_const->fin_not_cap_str,
                                eval_const->fin_not_std,
                                eval_const->fin_not_half_std);

                            /*	Coupon pv */

                            pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                            eval_const->fin_not_abs_beta * (floor - cap));
                        }
                        else
                        {
                            pd_leg += eval_const->next_pd_dff *
                                      exp(-eval_const->next_pd_gam * svijk[0] -
                                          eval_const->next_pd_gam2) *
                                      next_call->pd_not_amt;
                        }

                        /*	Initial notional */
                        fee =
                            eval_const->in_not_pd_dff *
                            exp(-eval_const->in_not_pd_gam * svijk[0] - eval_const->in_not_pd_gam2);

                        pd_leg -= fee * call->pd_not_amt;

                        fee *= call->fee;

                        /*	Finally, intrinsic value */

                        if (call->pay_rec == 0)
                        {
                            iv = pd_leg - fund_leg;
                        }
                        else
                        {
                            iv = fund_leg - pd_leg;
                        }

                        prod_val[i][j][k][0] -= iv;
                    }
                }
            }
        }
    }

    /*	End of payoff valuation
            ----------------------- */

    return err;
}

/*	Main pricing function */

/*	Launch the tree */
Err cpd_launch_tree(
    CPD_STR      cpd,
    CPD_UND      und,
    CPD_TREE_ARG tree_arg,
    /*	Result */
    double* prem,
    double* IV)
{
    double* temp_val;
    int     nprod;
    Err     err = NULL;

    nprod    = (int)(*IV + 1e-8);
    temp_val = (double*)calloc(nprod + 1, sizeof(double));
    if (!temp_val)
    {
        err = "cpd_launch_tree : Memory allocation error";
        return err;
    }

    if (und->use_beta_dlm)
    {
        err = tree_main_3dBetaDLM_QBeta(
            tree_arg->nstp,
            tree_arg->time,
            tree_arg->date,
            tree_arg->vol_change,
            tree_arg->dom_lam,
            tree_arg->for_lam,
            tree_arg->sig_dom,
            tree_arg->sig_for,
            tree_arg->sig_fx,

            tree_arg->corr_dom_for,
            tree_arg->corr_dom_fx,
            tree_arg->corr_for_fx,

            tree_arg->dr_const_dom,
            tree_arg->dr_coef_dom,
            tree_arg->dr_const_for,
            tree_arg->dr_coef1_for,
            tree_arg->dr_coef2_for,
            tree_arg->dr_coef3_for,
            tree_arg->dr_const_fx,
            tree_arg->dr_coef1_fx,
            tree_arg->dr_coef2_fx,
            tree_arg->dr_coef3_fx,

            tree_arg->dom_fwd,
            tree_arg->dom_var,
            tree_arg->for_fwd,
            tree_arg->for_var,
            tree_arg->fx_fwd,
            tree_arg->fx_var,

            tree_arg->glob_corr_dom_for,
            tree_arg->glob_corr_dom_fx,
            tree_arg->glob_corr_for_fx,

            tree_arg->spot_fx,
            tree_arg->dom_yc,
            tree_arg->for_yc,
            tree_arg->dom_ifr,

            tree_arg->void_prm,
            tree_arg->is_event,
            tree_arg->bar_lvl,
            tree_arg->bar_cl,
            tree_arg->is_bar,

            cpd_payoff_4_3dfx_BetaDLM_tree,
            nprod,
            1,
            (double*)temp_val);
    }
    else
    {
        err = tree_main_3dfx(
            tree_arg->nstp,
            tree_arg->time,
            tree_arg->date,
            tree_arg->vol_change,
            tree_arg->sig_dom,
            tree_arg->sig_for,
            tree_arg->sig_fx,
            tree_arg->dom_ifr,
            tree_arg->dom_fwd,
            tree_arg->dom_var,
            tree_arg->for_ifr,
            tree_arg->for_fwd,
            tree_arg->for_var,
            tree_arg->fx_fwd,
            tree_arg->fx_var,
            tree_arg->void_prm,
            tree_arg->is_event,
            tree_arg->bar_lvl,
            tree_arg->bar_cl,
            tree_arg->is_bar,
            tree_arg->dom_lam,
            tree_arg->for_lam,
            tree_arg->corr_dom_for,
            tree_arg->corr_dom_fx,
            tree_arg->corr_for_fx,
            tree_arg->spot_fx,
            tree_arg->dom_yc,
            tree_arg->for_yc,
            cpd_payoff_4_3dfx_tree,
            nprod,
            1,
            (double*)temp_val);
    }

    if (nprod == 2)
        *IV = temp_val[1];

    *prem = temp_val[0];

    if (temp_val)
        free(temp_val);
    return err;
}

/*	Launch the tree in case of a Callable KO */
Err cpd_launch_tree_ko(
    CPD_STR      cpd,
    CPD_UND      und,
    CPD_TREE_ARG tree_arg,
    /*	Result */
    double* prem)
{
    double temp_val[2];
    Err    err = NULL;

    err = tree_main_3dfx(
        tree_arg->nstp,
        tree_arg->time,
        tree_arg->date,
        tree_arg->vol_change,
        tree_arg->sig_dom,
        tree_arg->sig_for,
        tree_arg->sig_fx,
        tree_arg->dom_ifr,
        tree_arg->dom_fwd,
        tree_arg->dom_var,
        tree_arg->for_ifr,
        tree_arg->for_fwd,
        tree_arg->for_var,
        tree_arg->fx_fwd,
        tree_arg->fx_var,
        tree_arg->void_prm,
        tree_arg->is_event,
        tree_arg->bar_lvl,
        tree_arg->bar_cl,
        tree_arg->is_bar,
        tree_arg->dom_lam,
        tree_arg->for_lam,
        tree_arg->corr_dom_for,
        tree_arg->corr_dom_fx,
        tree_arg->corr_for_fx,
        tree_arg->spot_fx,
        tree_arg->dom_yc,
        tree_arg->for_yc,
        cpd_payoff_4_3dfx_tree_ko,
        1,
        1,
        (double*)temp_val);

    *prem = temp_val[0];

    return err;
}

/*	Payoff function for mc (KO) */
/*	---------------------------	*/

Err cpd_fut_intr_val(
    /* State variables */
    double Xdomfor[], /*	[0]: Xdom, [1]: Xfor */
    double Zfx,
    double spot_fx,
    /* Structure */
    CPD_STR cpd,
    /*	Call */
    int    fund_idx,
    int    num_fund_cpn,
    double fund_not_amt,
    int    pd_idx,
    int    num_pd_cpn,
    double pd_not_amt,
    /* Precalculated constants */
    CPD_EVAL_CONST eval_const,
    /* Results */
    double* fund_leg,
    double* pd_leg)
{
    double     fx, df, fwd, cap, floor, opt_string;
    PD_EXO_CPN cpn;
    PD_CALL    call;
    int        l, str_idx, opt_type;

    fx = spot_fx * exp(Zfx);

    /*	PV of funding leg */

    *fund_leg = 0.0;

    /*	Libor */

    *fund_leg +=
        eval_const->start_dff *
        exp(-eval_const->start_gam * Xdomfor[eval_const->fund_var] - eval_const->start_gam2) *
        cpd->fund_leg->notional;

    /*	Coupons */

    for (l = 0; l < num_fund_cpn; l++)
    {
        *fund_leg += eval_const->fund_dff[l] *
                     exp(-eval_const->fund_gam[l] * Xdomfor[eval_const->fund_var] -
                         eval_const->fund_gam2[l]) *
                     cpd->fund_leg->cpn[fund_idx + l].cpn;
    }

    /*	Initial notional */

    *fund_leg -= eval_const->in_not_fund_dff *
                 exp(-eval_const->in_not_fund_gam * Xdomfor[eval_const->fund_var] -
                     eval_const->in_not_fund_gam2) *
                 fund_not_amt;

    /*	Fx */

    if (cpd->fund_leg->dom_for != 0)
    {
        *fund_leg *= fx;
    }

    /*	PV of pd leg */

    *pd_leg = 0.0;

    /*	Coupons */
    for (l = 0; l < num_pd_cpn; l++)
    {
        /*	Coupon access */
        cpn = cpd->pd_leg->cpn + (pd_idx + l);

        /*	Discount */
        df = eval_const->pd_disc_dff[l] *
             exp(-eval_const->pd_disc_gam[l] * Xdomfor[0] - eval_const->pd_disc_gam2[l]);

        /*	Fwd fx */
        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
              exp(-eval_const->pd_fwd_for_gam[l] * Xdomfor[1] +
                  eval_const->pd_fwd_dom_gam[l] * Xdomfor[0] + eval_const->pd_fwd_cvx[l]);

        // string of options case:
        if (cpn->use_opt_str)
        {
            opt_string = 0.0;
            for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++)
            {
                if (fabs(cpn->weights[str_idx]) > 1e-16)
                {
                    if (eval_const->pd_std[l] > 1e-16 && cpn->strikes[str_idx] > 1e-16)
                        opt_type = 3;
                    else
                        opt_type = 1;

                    opt_string += cpn->weights[str_idx] * OPT_VAL_MACRO(
                                                              opt_type,
                                                              fwd,
                                                              cpn->strikes[str_idx],
                                                              eval_const->pd_std[l],
                                                              eval_const->pd_half_std[l]);
                }
            }
            *pd_leg += df * (cpn->wcst + cpn->wspot * fwd + opt_string);
        }
        else
        {
            /*	Floor */
            floor = OPT_VAL_MACRO(
                eval_const->pd_floor_type[l],
                fwd,
                eval_const->pd_floor_str[l],
                eval_const->pd_std[l],
                eval_const->pd_half_std[l]);

            /*	Cap */
            cap = OPT_VAL_MACRO(
                eval_const->pd_cap_type[l],
                fwd,
                eval_const->pd_cap_str[l],
                eval_const->pd_std[l],
                eval_const->pd_half_std[l]);

            /*	Coupon pv */

            *pd_leg +=
                df * (cpn->alpha + cpn->beta * fwd + eval_const->pd_abs_beta[l] * (floor - cap));
        }
    }
    /*	Final notional */

    /*	Coupon access */
    cpn = &(cpd->pd_leg->not_ref);

    /*	Discount */
    df = eval_const->fin_not_disc_dff *
         exp(-eval_const->fin_not_disc_gam * Xdomfor[0] - eval_const->fin_not_disc_gam2);

    /*	Fwd fx */
    fwd = fx * eval_const->fin_not_fwd_dff_ratio *
          exp(-eval_const->fin_not_fwd_for_gam * Xdomfor[1] +
              eval_const->fin_not_fwd_dom_gam * Xdomfor[0] + eval_const->fin_not_fwd_cvx);

    // interp redemption case:
    if (cpn->use_opt_str)
    {
        opt_string = 0.0;
        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++)
        {
            if (fabs(cpn->weights[str_idx]) > 1e-16)
            {
                if (eval_const->fin_not_std > 1e-16 && cpn->strikes[str_idx] > 1e-16)
                    opt_type = 3;
                else
                    opt_type = 1;

                opt_string += cpn->weights[str_idx] * OPT_VAL_MACRO(
                                                          opt_type,
                                                          fwd,
                                                          cpn->strikes[str_idx],
                                                          eval_const->fin_not_std,
                                                          eval_const->fin_not_half_std);
            }
        }
        *pd_leg += df * (cpn->wcst + cpn->wspot * fwd + opt_string);
    }
    else
    {
        /*	Floor */
        floor = OPT_VAL_MACRO(
            eval_const->fin_not_floor_type,
            fwd,
            eval_const->fin_not_floor_str,
            eval_const->fin_not_std,
            eval_const->fin_not_half_std);

        /*	Cap */
        cap = OPT_VAL_MACRO(
            eval_const->fin_not_cap_type,
            fwd,
            eval_const->fin_not_cap_str,
            eval_const->fin_not_std,
            eval_const->fin_not_half_std);

        /*	Coupon pv */

        *pd_leg +=
            df * (cpn->alpha + cpn->beta * fwd + eval_const->fin_not_abs_beta * (floor - cap));
    }

    /*	Initial notional */

    call = cpd->call + eval_const->call_idx;

    df = eval_const->in_not_pd_dff *
         exp(-eval_const->in_not_pd_gam * Xdomfor[0] - eval_const->in_not_pd_gam2);

    if (call->use_opt_str)
    {
        fwd = fx * eval_const->in_not_pd_fwd_dff_ratio *
              exp(-eval_const->in_not_pd_fwd_for_gam * Xdomfor[1] +
                  eval_const->in_not_pd_fwd_dom_gam * Xdomfor[0] + eval_const->in_not_pd_fwd_cvx);

        opt_string = 0.0;
        for (str_idx = 0; str_idx < call->nstrikes; str_idx++)
        {
            if (fabs(call->weights[str_idx]) > 1e-16)
            {
                if (eval_const->in_not_pd_std > 1e-16 && call->strikes[str_idx] > 1e-16)
                    opt_type = 3;
                else
                    opt_type = 1;

                opt_string += call->weights[str_idx] * OPT_VAL_MACRO(
                                                           opt_type,
                                                           fwd,
                                                           call->strikes[str_idx],
                                                           eval_const->in_not_pd_std,
                                                           eval_const->in_not_pd_half_std);
            }
        }
        *pd_leg -= df * (call->wcst + call->wspot * fwd + opt_string + eval_const->fee);
    }
    else
    {
        *pd_leg -= df * (pd_not_amt + eval_const->fee);
    }

    return NULL;
}

Err cpd_current_pd_cpn_val(
    /* State variables */
    double Xdomfor[], /*	[0]: Xdom, [1]: Xfor */
    double Zfx,
    double spot_fx,
    /* Structure */
    CPD_STR cpd,
    /*	Call */
    int    pd_idx,
    int    num_pd_cpn,
    double pd_not_amt,
    /* Precalculated constants */
    CPD_EVAL_CONST eval_const,
    /* Results */
    double* pd_leg)
{
    double     fx, df, cap, floor, opt_string, call;
    PD_EXO_CPN cpn;
    int        str_idx;

    fx = spot_fx * exp(Zfx);

    /*	PV of pd leg */

    *pd_leg = 0.0;

    /*	Coupon access */
    cpn = cpd->pd_leg->cpn + pd_idx;

    /*	Discount */  // To do
    df = 1.0;

    // string of options case:
    if (cpn->use_opt_str)
    {
        opt_string = 0.0;
        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++)
        {
            if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not_amt)
            {
                call = fx - cpn->strikes[str_idx];
                if (call > 0.0)
                    opt_string += cpn->weights[str_idx] * call;
            }
        }
        *pd_leg += df * (cpn->wcst + cpn->wspot * fx + opt_string);
    }
    else
    {
        if (cpn->floored && fabs(cpn->beta) > 1.0e-10)
        {
            if (cpn->beta > 0.0)
                floor = (cpn->floor - cpn->alpha) / cpn->beta - fx; /*	Put IV */
            else
                floor = fx - (cpn->floor - cpn->alpha) / cpn->beta; /*	Call IV */

            if (floor < 0.0)
                floor = 0.0;
        }
        else
            floor = 0.0;

        if (cpn->capped && fabs(cpn->beta) > 1.0e-10)
        {
            if (cpn->beta < 0.0)
                cap = (cpn->cap - cpn->alpha) / cpn->beta - fx; /*	Put IV */
            else
                cap = fx - (cpn->cap - cpn->alpha) / cpn->beta; /*	Call IV */

            if (cap < 0.0)
                cap = 0.0;
        }
        else
            cap = 0.0;

        /*	Coupon pv */

        *pd_leg += df * (cpn->alpha + cpn->beta * fx + fabs(cpn->beta) * (floor - cap));
    }

    return NULL;
}

Err cpd_price_FxInstrument_dlm(
    double         Fwd,
    double         LogFwd,
    double         X0,
    CPD_UND        und,
    CPD_STR        cpd,
    int            index_cpn,
    PD_EXO_CPN     cpn,
    CPD_EVAL_CONST eval_const,
    double*        prices)
{
    FxBetaDLM_FxOptInst*   Inst;
    FxBetaDLM_InstPrecalc* InstConst;
    Err                    err = NULL;
    int                    index_X1, index_X2, index_Fwd1, index_Fwd2;
    double                 price1, price2, vol1, vol2, vol;
    double*                dPrecalcPriceStdInX;
    double*                dPrecalcPriceStdAlphaInX;

    if (index_cpn >= 0)
    {
        Inst      = &(eval_const->pd_inst[index_cpn]);
        InstConst = &(eval_const->pd_beta_precalc[index_cpn]);
    }
    else
    {
        Inst      = &(eval_const->fin_not_inst);
        InstConst = &(eval_const->fin_not_beta_precalc);
    }

    if (und->num_params->iPrecalcSmile)
    {
        /* find the volatility */
        index_X1   = Get_Index(X0, InstConst->dPrecalcX, und->num_params->iNbX);
        index_Fwd1 = Get_Index(Fwd, InstConst->dPrecalcFwdFloor, und->num_params->iNbFwd);

        if (und->num_params->iPriceVol == 0)  // We do an interpolation of the prices
        {
            if (index_X1 == 0)
            {
                index_X1 = 1;
            }

            if (index_Fwd1 == 0)
            {
                index_Fwd1 = 1;
            }
            index_Fwd2 = index_Fwd1 - 1;
            index_X2   = index_X1 - 1;

            /* Compute price 1 */
            dPrecalcPriceStdInX      = InstConst->dPrecalcFloorPriceStd[index_X1];
            dPrecalcPriceStdAlphaInX = InstConst->dPrecalcFloorPriceStdAlpha[index_X1];
            price1                   = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                         (Fwd - InstConst->dPrecalcFwdFloor[index_Fwd2]) +
                     dPrecalcPriceStdInX[index_Fwd2];

            /* Compute price2 */
            dPrecalcPriceStdInX      = InstConst->dPrecalcFloorPriceStd[index_X2];
            dPrecalcPriceStdAlphaInX = InstConst->dPrecalcFloorPriceStdAlpha[index_X2];
            price2                   = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                         (Fwd - InstConst->dPrecalcFwdFloor[index_Fwd2]) +
                     dPrecalcPriceStdInX[index_Fwd2];

            if ((price1 <= 0.0) && (price2 <= 0.0))
            {
                prices[0] = 0.0;
            }
            else if (price2 <= 0.0)
            {
                prices[0] = price1 /
                            (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                            (X0 - InstConst->dPrecalcX[index_X2]);
            }
            else if (price1 <= 0.0)
            {
                prices[0] = -price2 /
                                (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                                (X0 - InstConst->dPrecalcX[index_X2]) +
                            price2;
            }
            else
            {
                prices[0] = (price1 - price2) /
                                (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                                (X0 - InstConst->dPrecalcX[index_X2]) +
                            price2;
            }

            if (index_Fwd1 > 10)
            {
                index_Fwd1 = index_Fwd1;
            }

            if (Inst->iNbStrike > 2)
            {
                index_Fwd1 = Get_Index(Fwd, InstConst->dPrecalcFwdCap, und->num_params->iNbFwd);
                if (index_Fwd1 == 0)
                {
                    index_Fwd1 = 1;
                }
                index_Fwd2 = index_Fwd1 - 1;

                /* we need to value the cap as well */
                /* Compute price 1 */
                dPrecalcPriceStdInX      = InstConst->dPrecalcCapPriceStd[index_X1];
                dPrecalcPriceStdAlphaInX = InstConst->dPrecalcCapPriceStdAlpha[index_X1];
                price1                   = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                             (Fwd - InstConst->dPrecalcFwdCap[index_Fwd2]) +
                         dPrecalcPriceStdInX[index_Fwd2];

                /* Compute price2 */
                dPrecalcPriceStdInX      = InstConst->dPrecalcCapPriceStd[index_X2];
                dPrecalcPriceStdAlphaInX = InstConst->dPrecalcCapPriceStdAlpha[index_X2];
                price2                   = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                             (Fwd - InstConst->dPrecalcFwdCap[index_Fwd2]) +
                         dPrecalcPriceStdInX[index_Fwd2];

                if ((price1 <= 0.0) && (price2 <= 0.0))
                {
                    prices[1] = 0.0;
                }
                else if (price2 <= 0.0)
                {
                    prices[1] = price1 /
                                (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                                (X0 - InstConst->dPrecalcX[index_X2]);
                }
                else if (price1 <= 0.0)
                {
                    prices[1] =
                        -price2 /
                            (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                            (X0 - InstConst->dPrecalcX[index_X2]) +
                        price2;
                }
                else
                {
                    prices[1] =
                        (price1 - price2) /
                            (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                            (X0 - InstConst->dPrecalcX[index_X2]) +
                        price2;
                }
                prices[2] = Fwd;
            }
            else
            {
                prices[1] = Fwd;
            }
        }
        else  // We do an interpolation of the volatilities
        {
            dPrecalcPriceStdInX      = InstConst->dPrecalcFloorPriceStd[index_X1];
            dPrecalcPriceStdAlphaInX = InstConst->dPrecalcFloorPriceStdAlpha[index_X1];

            /* Compute vol 1 */
            if (Fwd < InstConst->dPrecalcFwdFloor[0])
            {
                vol1 = dPrecalcPriceStdInX[0];
            }
            else if (Fwd > InstConst->dPrecalcFwdFloor[und->num_params->iNbFwd - 1])
            {
                vol1 = dPrecalcPriceStdInX[und->num_params->iNbFwd - 1];
            }
            else
            {
                /* interpolation */
                index_Fwd2 = index_Fwd1 - 1;

                vol1 = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                           (Fwd - InstConst->dPrecalcFwdFloor[index_Fwd2]) +
                       dPrecalcPriceStdInX[index_Fwd2];

                // vol1 = (dPrecalcPriceStdInX[index_Fwd1] - dPrecalcPriceStdInX[index_Fwd2])
                //	/ (InstConst->dPrecalcFwdFloor[index_Fwd1] -
                // InstConst->dPrecalcFwdFloor[index_Fwd2])
                //	* (Fwd - InstConst->dPrecalcFwdFloor[index_Fwd2]) +
                // dPrecalcPriceStdInX[index_Fwd2];
            }

            /* Compute Vol2 if needed */
            if (X0 < InstConst->dPrecalcX[0])
            {
                vol = vol1;
            }
            else if (X0 > InstConst->dPrecalcX[und->num_params->iNbX - 1])
            {
                vol = vol1;
            }
            else
            {
                /* interpolation */
                index_X2                 = index_X1 - 1;
                dPrecalcPriceStdInX      = InstConst->dPrecalcFloorPriceStd[index_X2];
                dPrecalcPriceStdAlphaInX = InstConst->dPrecalcFloorPriceStdAlpha[index_X2];

                /* Compute vol 1 */
                if (Fwd < InstConst->dPrecalcFwdFloor[0])
                {
                    vol2 = dPrecalcPriceStdInX[0];
                }
                else if (Fwd > InstConst->dPrecalcFwdFloor[und->num_params->iNbFwd - 1])
                {
                    vol2 = dPrecalcPriceStdInX[und->num_params->iNbFwd - 1];
                }
                else
                {
                    /* interpolation */
                    index_Fwd2 = index_Fwd1 - 1;

                    vol2 = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                               (Fwd - InstConst->dPrecalcFwdFloor[index_Fwd2]) +
                           dPrecalcPriceStdInX[index_Fwd2];

                    // vol2 = (dPrecalcPriceStdInX[index_Fwd1] - dPrecalcPriceStdInX[index_Fwd2])
                    //	/ (InstConst->dPrecalcFwdFloor[index_Fwd1] -
                    // InstConst->dPrecalcFwdFloor[index_Fwd2])
                    //	* (Fwd - InstConst->dPrecalcFwdFloor[index_Fwd2]) +
                    // dPrecalcPriceStdInX[index_Fwd2];
                }

                vol = (vol1 - vol2) /
                          (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                          (X0 - InstConst->dPrecalcX[index_X2]) +
                      vol2;
            }

            prices[0] = OPT_LOGVAL_MACRO(
                3, Fwd, LogFwd, Inst->dStrike[0], Inst->dLogStrike[0], vol, 0.5 * vol);

            /*prices[0] = OPT_VAL_MACRO(
                                                            3,
                                                            Fwd,
                                                            Inst->dStrike[0],
                                                            vol,
                                                            0.5 * vol);
            */
            if (Inst->iNbStrike > 2)
            {
                /* we need to value the cap as well */
                dPrecalcPriceStdInX      = InstConst->dPrecalcCapPriceStd[index_X1];
                dPrecalcPriceStdAlphaInX = InstConst->dPrecalcCapPriceStdAlpha[index_X1];

                index_Fwd1 = Get_Index(Fwd, InstConst->dPrecalcFwdCap, und->num_params->iNbFwd);

                /* Compute vol 1 */
                if (Fwd < InstConst->dPrecalcFwdCap[0])
                {
                    vol1 = dPrecalcPriceStdInX[0];
                }
                else if (Fwd > InstConst->dPrecalcFwdCap[und->num_params->iNbFwd - 1])
                {
                    vol1 = dPrecalcPriceStdInX[und->num_params->iNbFwd - 1];
                }
                else
                {
                    /* interpolation */
                    index_Fwd2 = index_Fwd1 - 1;

                    vol1 = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                               (Fwd - InstConst->dPrecalcFwdCap[index_Fwd2]) +
                           dPrecalcPriceStdInX[index_Fwd2];

                    //*vol1 = (dPrecalcPriceStdInX[index_Fwd1] - dPrecalcPriceStdInX[index_Fwd2])
                    //	/ (InstConst->dPrecalcFwdCap[index_Fwd1] -
                    // InstConst->dPrecalcFwdCap[index_Fwd2])
                    //	* (Fwd - InstConst->dPrecalcFwdCap[index_Fwd2]) +
                    // dPrecalcPriceStdInX[index_Fwd2];
                }

                /* Compute Vol2 if needed */
                if (X0 < InstConst->dPrecalcX[0])
                {
                    vol = vol1;
                }
                else if (X0 > InstConst->dPrecalcX[und->num_params->iNbX - 1])
                {
                    vol = vol1;
                }
                else
                {
                    /* interpolation */
                    dPrecalcPriceStdInX      = InstConst->dPrecalcCapPriceStd[index_X2];
                    dPrecalcPriceStdAlphaInX = InstConst->dPrecalcCapPriceStdAlpha[index_X2];

                    /* Compute vol 1 */
                    if (Fwd < InstConst->dPrecalcFwdCap[0])
                    {
                        vol2 = dPrecalcPriceStdInX[0];
                    }
                    else if (Fwd > InstConst->dPrecalcFwdCap[und->num_params->iNbFwd - 1])
                    {
                        vol2 = dPrecalcPriceStdInX[und->num_params->iNbFwd - 1];
                    }
                    else
                    {
                        /* interpolation */
                        index_Fwd2 = index_Fwd1 - 1;

                        vol2 = dPrecalcPriceStdAlphaInX[index_Fwd1] *
                                   (Fwd - InstConst->dPrecalcFwdCap[index_Fwd2]) +
                               dPrecalcPriceStdInX[index_Fwd2];

                        // vol2 = (dPrecalcPriceStdInX[index_Fwd1] -
                        // dPrecalcPriceStdInX[index_Fwd2]) 	/
                        // (InstConst->dPrecalcFwdCap[index_Fwd1]
                        //- InstConst->dPrecalcFwdCap[index_Fwd2])
                        //	* (Fwd - InstConst->dPrecalcFwdCap[index_Fwd2]) +
                        // dPrecalcPriceStdInX[index_Fwd2];
                    }

                    vol = (vol1 - vol2) /
                              (InstConst->dPrecalcX[index_X1] - InstConst->dPrecalcX[index_X2]) *
                              (X0 - InstConst->dPrecalcX[index_X2]) +
                          vol2;
                }
                prices[1] = OPT_LOGVAL_MACRO(
                    3, Fwd, LogFwd, Inst->dStrike[1], Inst->dLogStrike[1], vol, 0.5 * vol);
                /*prices[1] = OPT_VAL_MACRO(
                                                                3,
                                                                Fwd,
                                                                Inst->dStrike[1],
                                                                vol,
                                                                0.5 * vol);*/

                prices[2] = Fwd;
            }
            else
            {
                prices[1] = Fwd;
            }
        }
    }
    else
    {
        Inst->dLnFwdFx = LogFwd;

        err = FxBetaDLM_Update_InstPrecalc_FromX0(X0, Inst, und->model, und->num_params, InstConst);

        if (err)
            return err;

        err = FxBetaDLM_Price_FxOptInst(
            Inst, und->model, InstConst, und->hermite, und->num_params, prices);

        if (err)
            return err;
    }

    return err;
}

Err cpd_fut_intr_val_dlm(
    /* State variables */
    double Xdomfor[], /*	[0]: Xdom, [1]: Xfor */
    double X,
    double fx,
    /* Structure */
    CPD_STR cpd,
    CPD_UND und,
    /*	Call */
    int    fund_idx,
    int    num_fund_cpn,
    double fund_not_amt,
    int    pd_idx,
    int    num_pd_cpn,
    double pd_not_amt,
    /* Precalculated constants */
    CPD_EVAL_CONST eval_const,
    /* Results */
    double* fund_leg,
    double* pd_leg)
{
    double     df, fwd, fwd_adj, cap, floor;
    PD_EXO_CPN cpn;
    int        l;
    double     prices[3];
    Err        err = NULL;

    /*
    FILE		*stream = fopen("C:\\MCOption.txt", "a");
    */

    /*	PV of funding leg */

    *fund_leg = 0.0;

    /*	Libor */

    *fund_leg +=
        eval_const->start_dff *
        exp(-eval_const->start_gam * Xdomfor[eval_const->fund_var] - eval_const->start_gam2) *
        cpd->fund_leg->notional;

    /*	Coupons */

    for (l = 0; l < num_fund_cpn; l++)
    {
        *fund_leg += eval_const->fund_dff[l] *
                     exp(-eval_const->fund_gam[l] * Xdomfor[eval_const->fund_var] -
                         eval_const->fund_gam2[l]) *
                     cpd->fund_leg->cpn[fund_idx + l].cpn;
    }

    /*	Initial notional */

    *fund_leg -= eval_const->in_not_fund_dff *
                 exp(-eval_const->in_not_fund_gam * Xdomfor[eval_const->fund_var] -
                     eval_const->in_not_fund_gam2) *
                 fund_not_amt;

    /*	Fx */

    if (cpd->fund_leg->dom_for != 0)
    {
        *fund_leg *= fx;
    }

    /*	PV of pd leg */

    *pd_leg = 0.0;

    /*	Coupons */
    for (l = 0; l < num_pd_cpn; l++)
    {
        /*	Coupon access */
        cpn = cpd->pd_leg->cpn + (pd_idx + l);

        /*	Discount */
        df = eval_const->pd_disc_dff[l] *
             exp(-eval_const->pd_disc_gam[l] * Xdomfor[0] - eval_const->pd_disc_gam2[l]);

        /*	Fwd fx */
        fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
              exp(-eval_const->pd_fwd_for_gam[l] * Xdomfor[1] +
                  eval_const->pd_fwd_dom_gam[l] * Xdomfor[0] + eval_const->pd_fwd_cvx[l]);

        /*	Floor */

        if (eval_const->pd_inst[l].iNbStrike > 0)
        {
            err =
                cpd_price_FxInstrument_dlm(fwd, log(fwd), X, und, cpd, l, cpn, eval_const, prices);

            if (err)
                return err;
        }

        /* get the fwd */
        if (eval_const->index_fwd[l] >= 0)
        {
            fwd_adj = prices[eval_const->index_fwd[l]] +
                      eval_const->pd_inst[l].dStrike[eval_const->index_fwd[l]];
        }
        else
        {
            fwd_adj = fwd;
        }

        switch (eval_const->pd_floor_type[l])
        {
        case 0:
            floor = 0.0;
            break;
        case 1:
        case 2:
            floor = OPT_VAL_MACRO(
                eval_const->pd_floor_type[l], fwd_adj, eval_const->pd_floor_str[l], 1.0, 1.0);
            break;
        case 3:
            floor = prices[0];
            break;
        case 4:
            floor = prices[0] + eval_const->pd_floor_str[l] - fwd_adj;
            break;
        }

        switch (eval_const->pd_cap_type[l])
        {
        case 0:
            cap = 0.0;
            break;
        case 1:
        case 2:
            cap = OPT_VAL_MACRO(
                eval_const->pd_cap_type[l], fwd_adj, eval_const->pd_cap_str[l], 1.0, 1.0);
            break;
        case 3:
            cap = prices[eval_const->index_cap[l]];
            break;
        case 4:
            cap = prices[eval_const->index_cap[l]] + eval_const->pd_cap_str[l] - fwd_adj;
            break;
        }

        /*	Coupon pv */

        *pd_leg +=
            df * (cpn->alpha + cpn->beta * fwd_adj + eval_const->pd_abs_beta[l] * (floor - cap));
    }

    /*	Final notional */

    /*	Coupon access */
    cpn = &(cpd->pd_leg->not_ref);

    /*	Discount */
    df = eval_const->fin_not_disc_dff *
         exp(-eval_const->fin_not_disc_gam * Xdomfor[0] - eval_const->fin_not_disc_gam2);

    /*	Fwd fx */
    fwd = fx * eval_const->fin_not_fwd_dff_ratio *
          exp(-eval_const->fin_not_fwd_for_gam * Xdomfor[1] +
              eval_const->fin_not_fwd_dom_gam * Xdomfor[0] + eval_const->fin_not_fwd_cvx);

    if (eval_const->fin_not_inst.iNbStrike > 0)
    {
        err = cpd_price_FxInstrument_dlm(fwd, log(fwd), X, und, cpd, -1, cpn, eval_const, prices);

        if (err)
            return err;
    }

    /* get the fwd */
    if (eval_const->fin_not_index_fwd >= 0)
    {
        fwd_adj = prices[eval_const->fin_not_index_fwd];
    }
    else
    {
        fwd_adj = fwd;
    }

    switch (eval_const->fin_not_floor_type)
    {
    case 0:
        floor = 0.0;
        break;
    case 1:
    case 2:
        floor = OPT_VAL_MACRO(
            eval_const->fin_not_floor_type, fwd_adj, eval_const->fin_not_floor_str, 1.0, 1.0);
        break;
    case 3:
        floor = prices[0];
        break;
    case 4:
        floor = prices[0] + eval_const->fin_not_floor_str - fwd_adj;
        break;
    }

    switch (eval_const->fin_not_cap_type)
    {
    case 0:
        cap = 0.0;
        break;
    case 1:
    case 2:
        cap = OPT_VAL_MACRO(
            eval_const->fin_not_cap_type, fwd_adj, eval_const->fin_not_cap_str, 1.0, 1.0);
        break;
    case 3:
        cap = prices[eval_const->fin_not_index_cap];
        break;
    case 4:
        cap = prices[eval_const->fin_not_index_cap] + eval_const->fin_not_cap_str - fwd_adj;
        break;
    }

    /*	Coupon pv */

    *pd_leg += df * (cpn->alpha + cpn->beta * fwd + eval_const->fin_not_abs_beta * (floor - cap));

    /*	Initial notional */
    *pd_leg -= eval_const->in_not_pd_dff *
               exp(-eval_const->in_not_pd_gam * Xdomfor[0] - eval_const->in_not_pd_gam2) *
               (pd_not_amt + eval_const->fee);

    /*
    fclose(stream);
    */

    return NULL;
}

Err cpd_fut_total_val_dlm(
    /* State variables */
    double Xdomfor[], /*	[0]: Xdom, [1]: Xfor */
    double X,
    double fx,
    /* Structure */
    CPD_STR cpd,
    CPD_UND und,
    /*	Call */
    int    fund_idx,
    int    num_fund_cpn,
    double fund_not_amt,
    double next_fund_not_amt,
    int    pd_idx,
    int    num_pd_cpn,
    double pd_not_amt,
    double next_pd_not_amt,
    int    eval_redemption,
    /* Precalculated constants */
    CPD_EVAL_CONST eval_const,
    /* Results */
    int     calc_coupon,
    double* fund_leg,
    double* pd_leg,
    double* fee)
{
    double     df, fwd, fwd_adj, cap, floor;
    PD_EXO_CPN cpn;
    int        l;
    double     prices[3];
    Err        err = NULL;

    if (calc_coupon)
    {
        /* funding leg */
        *fund_leg = 0.0;

        *fund_leg +=
            eval_const->start_dff *
            exp(-eval_const->start_gam * Xdomfor[eval_const->fund_var] - eval_const->start_gam2) *
            cpd->fund_leg->notional;

        for (l = 0; l < num_fund_cpn; l++)
        {
            *fund_leg += eval_const->fund_dff[l] *
                         exp(-eval_const->fund_gam[l] * Xdomfor[eval_const->fund_var] -
                             eval_const->fund_gam2[l]) *
                         cpd->fund_leg->cpn[fund_idx + l].cpn;
        }

        *fund_leg -= eval_const->in_not_fund_dff *
                     exp(-eval_const->in_not_fund_gam * Xdomfor[eval_const->fund_var] -
                         eval_const->in_not_fund_gam2) *
                     fund_not_amt;

        if (!eval_redemption)
        {
            *fund_leg -= eval_const->next_start_fund_dff *
                         exp(-eval_const->next_start_fund_gam * Xdomfor[eval_const->fund_var] -
                             eval_const->next_start_fund_gam2) *
                         fund_not_amt;

            *fund_leg += eval_const->next_fund_dff *
                         exp(-eval_const->next_fund_gam * Xdomfor[eval_const->fund_var] -
                             eval_const->next_fund_gam2) *
                         next_fund_not_amt;
        }

        /*	Fx */
        if (cpd->fund_leg->dom_for != 0)
        {
            *fund_leg *= fx;
        }

        /* PD leg */

        *pd_leg = 0.0;

        for (l = 0; l < num_pd_cpn; l++)
        {
            /*	Coupon access */
            cpn = cpd->pd_leg->cpn + (pd_idx + l);

            /*	Discount */
            df = eval_const->pd_disc_dff[l] *
                 exp(-eval_const->pd_disc_gam[l] * Xdomfor[0] - eval_const->pd_disc_gam2[l]);

            /*	Fwd fx */
            fwd = fx * eval_const->pd_fwd_dff_ratio[l] *
                  exp(-eval_const->pd_fwd_for_gam[l] * Xdomfor[1] +
                      eval_const->pd_fwd_dom_gam[l] * Xdomfor[0] + eval_const->pd_fwd_cvx[l]);

            /*	Floor */

            if (eval_const->pd_inst[l].iNbStrike > 0)
            {
                err = cpd_price_FxInstrument_dlm(
                    fwd, log(fwd), X, und, cpd, l, cpn, eval_const, prices);

                if (err)
                    return err;
            }

            /* get the fwd */
            if (eval_const->index_fwd[l] >= 0)
            {
                fwd_adj = prices[eval_const->index_fwd[l]] +
                          eval_const->pd_inst[l].dStrike[eval_const->index_fwd[l]];
            }
            else
            {
                fwd_adj = fwd;
            }

            switch (eval_const->pd_floor_type[l])
            {
            case 0:
                floor = 0.0;
                break;
            case 1:
            case 2:
                floor = OPT_VAL_MACRO(
                    eval_const->pd_floor_type[l], fwd_adj, eval_const->pd_floor_str[l], 1.0, 1.0);
                break;
            case 3:
                floor = prices[0];
                break;
            case 4:
                floor = prices[0] + eval_const->pd_floor_str[l] - fwd_adj;
                break;
            }

            switch (eval_const->pd_cap_type[l])
            {
            case 0:
                cap = 0.0;
                break;
            case 1:
            case 2:
                cap = OPT_VAL_MACRO(
                    eval_const->pd_cap_type[l], fwd_adj, eval_const->pd_cap_str[l], 1.0, 1.0);
                break;
            case 3:
                cap = prices[eval_const->index_cap[l]];
                break;
            case 4:
                cap = prices[eval_const->index_cap[l]] + eval_const->pd_cap_str[l] - fwd_adj;
                break;
            }

            /*	Coupon pv */

            *pd_leg += df * (cpn->alpha + cpn->beta * fwd_adj +
                             eval_const->pd_abs_beta[l] * (floor - cap));
        }

        if (eval_redemption)
        {
            /*	Coupon access */
            cpn = &(cpd->pd_leg->not_ref);

            /*	Discount */
            df = eval_const->fin_not_disc_dff *
                 exp(-eval_const->fin_not_disc_gam * Xdomfor[0] - eval_const->fin_not_disc_gam2);

            /*	Fwd fx */
            fwd = fx * eval_const->fin_not_fwd_dff_ratio *
                  exp(-eval_const->fin_not_fwd_for_gam * Xdomfor[1] +
                      eval_const->fin_not_fwd_dom_gam * Xdomfor[0] + eval_const->fin_not_fwd_cvx);

            if (eval_const->fin_not_inst.iNbStrike > 0)
            {
                err = cpd_price_FxInstrument_dlm(
                    fwd, log(fwd), X, und, cpd, -1, cpn, eval_const, prices);

                if (err)
                    return err;
            }

            /* get the fwd */
            if (eval_const->fin_not_index_fwd >= 0)
            {
                fwd_adj = prices[eval_const->fin_not_index_fwd];
            }
            else
            {
                fwd_adj = fwd;
            }

            switch (eval_const->fin_not_floor_type)
            {
            case 0:
                floor = 0.0;
                break;
            case 1:
            case 2:
                floor = OPT_VAL_MACRO(
                    eval_const->fin_not_floor_type,
                    fwd_adj,
                    eval_const->fin_not_floor_str,
                    1.0,
                    1.0);
                break;
            case 3:
                floor = prices[0];
                break;
            case 4:
                floor = prices[0] + eval_const->fin_not_floor_str - fwd_adj;
                break;
            }

            switch (eval_const->fin_not_cap_type)
            {
            case 0:
                cap = 0.0;
                break;
            case 1:
            case 2:
                cap = OPT_VAL_MACRO(
                    eval_const->fin_not_cap_type, fwd_adj, eval_const->fin_not_cap_str, 1.0, 1.0);
                break;
            case 3:
                cap = prices[eval_const->fin_not_index_cap];
                break;
            case 4:
                cap = prices[eval_const->fin_not_index_cap] + eval_const->fin_not_cap_str - fwd_adj;
                break;
            }

            /*	Coupon pv */

            *pd_leg +=
                df * (cpn->alpha + cpn->beta * fwd + eval_const->fin_not_abs_beta * (floor - cap));
        }
        else
        {
            *pd_leg += eval_const->next_pd_dff *
                       exp(-eval_const->next_pd_gam * Xdomfor[0] - eval_const->next_pd_gam2) *
                       next_pd_not_amt;
        }
    }

    *fee = eval_const->in_not_pd_dff *
           exp(-eval_const->in_not_pd_gam * Xdomfor[0] - eval_const->in_not_pd_gam2);

    if (calc_coupon)
    {
        *pd_leg -= (*fee) * pd_not_amt;
    }

    *fee *= eval_const->fee;

    return NULL;
}

/*	Remaining notional of the structure */
static double cpd_rem_not_4_3dfx_mc;
#define REMNO cpd_rem_not_4_3dfx_mc

static double cpd_sum_coupon_paid_4_3dfx_tarn;
#define CPNPAID cpd_sum_coupon_paid_4_3dfx_tarn

/*	Reset remaining notional */
void cpd_init_4_3dfx_mc(void)
{
    REMNO   = 1.0;
    CPNPAID = 0.0;
}

Err cpd_payoff_4_3dfx_mc_dlm(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call;
    CPD_UND        und;
    CPD_EVAL_CONST eval_const;
    int            call_idx;
    double         fx, fx_now;
    double         Xdomfor[2];

    double fund_leg, pd_leg, iv;
    double SMO, part;

    Err err = NULL;

    *res       = 0.0;
    *stop_path = 0;

    if (!func_parm)
    {
        return NULL;
    }

    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    SMO = call->smooth;

    Xdomfor[0] = Xdom;
    Xdomfor[1] = Yfor;

    fx =
        exp(eval_const->pd_S_const + Zfx * (eval_const->pd_S_lin + Zfx * eval_const->pd_S_quad) +
            eval_const->pd_S_dom * Xdom + eval_const->pd_S_for * Yfor);

    fx_now = log(fx / und->model->dCashFx) + eval_const->spot_fx_dff_log_ratio -
             eval_const->spot_fx_for_gam * Yfor + eval_const->spot_fx_dom_gam * Xdom +
             eval_const->spot_fx_cvx;

    if (call->fx_bound && ((call->bar_type == 0 && fx_now < call->barrier - SMO) ||
                           (call->bar_type == 1 && fx_now > call->barrier + SMO)))
    {
        return NULL;
    }

    err = cpd_fut_intr_val_dlm(
        Xdomfor,
        Zfx,
        fx,
        cpd,
        und,
        call->fund_idx,
        call->num_fund_cpn,
        call->fund_not_amt,
        call->pd_idx,
        call->num_pd_cpn,
        call->pd_not_amt,
        eval_const,
        &fund_leg,
        &pd_leg);

    if (err)
        return err;

    if (call->pay_rec == 0)
    {
        iv = REMNO * (pd_leg - fund_leg);
    }
    else
    {
        iv = REMNO * (fund_leg - pd_leg);
    }

    if (call->fx_bound)
    {
        if (call->bar_type == 0)
        {
            if (iv > 0.0)
            {
                if (fx_now <= call->barrier)
                {
                    return NULL;
                }
                else if (fx_now >= call->barrier + SMO)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now >= call->barrier)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier + SMO) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
        else
        {
            if (iv > 0.0)
            {
                if (fx_now >= call->barrier)
                {
                    return NULL;
                }
                else if (fx_now <= call->barrier - SMO)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier - fx_now) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now <= call->barrier)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier + SMO - fx_now) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
    }
    else
    {
        /* we knock on the IV !!!! */
        /* no smoothing */
        if (call->bar_type == 0)
        {
            if (iv <= call->barrier)
            {
                return NULL;
            }
            else
            {
                *res       = iv;
                *stop_path = 1;
                return NULL;
            }
        }
        else
        {
            if (iv >= call->barrier)
            {
                return NULL;
            }
            else
            {
                *res       = iv;
                *stop_path = 1;
                return NULL;
            }
        }
    }

    return NULL;
}

Err cpd_payoff_total_4_3dfx_mc_dlm(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call, next_call;
    CPD_UND        und;
    CPD_EVAL_CONST eval_const;
    int            call_idx;
    double         fx, fx_now;
    double         Xdomfor[2];

    double fund_leg, pd_leg, fee, iv;
    double SMO, part;

    int num_fund_cpn, num_pd_cpn, eval_redemption, calc_coupon;

    Err err = NULL;

    *res       = 0.0;
    *stop_path = 0;

    if (!func_parm)
    {
        return NULL;
    }

    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    /*	Calculate number of coupons to be computed */
    if (call_idx == cpd->num_calls - 1)
    {
        /*	Last call: eval all + redemption */
        next_call       = cpd->call;
        num_fund_cpn    = call->num_fund_cpn;
        num_pd_cpn      = call->num_pd_cpn;
        eval_redemption = 1;
    }
    else
    {
        /*	Not the last call: only eval coupons up to next call date */
        next_call       = cpd->call + (call_idx + 1);
        num_fund_cpn    = next_call->fund_idx - call->fund_idx;
        num_pd_cpn      = next_call->pd_idx - call->pd_idx;
        eval_redemption = 0;
    }

    SMO = call->smooth;

    Xdomfor[0] = Xdom;
    Xdomfor[1] = Yfor;

    fx =
        exp(eval_const->pd_S_const + Zfx * (eval_const->pd_S_lin + Zfx * eval_const->pd_S_quad) +
            eval_const->pd_S_dom * Xdom + eval_const->pd_S_for * Yfor);

    fx_now = log(fx / und->model->dCashFx) + eval_const->spot_fx_dff_log_ratio -
             eval_const->spot_fx_for_gam * Yfor + eval_const->spot_fx_dom_gam * Xdom +
             eval_const->spot_fx_cvx;

    if (call->fx_bound && ((call->bar_type == 0 && fx_now < call->barrier + SMO) ||
                           (call->bar_type == 1 && fx_now > call->barrier - SMO)))
    {
        calc_coupon = 1;
    }
    else
    {
        calc_coupon = 0;
    }

    err = cpd_fut_total_val_dlm(
        Xdomfor,
        Zfx,
        fx,
        cpd,
        und,
        call->fund_idx,
        num_fund_cpn,
        call->fund_not_amt,
        next_call->fund_not_amt,
        call->pd_idx,
        num_pd_cpn,
        call->pd_not_amt,
        next_call->pd_not_amt,
        eval_redemption,
        eval_const,
        calc_coupon,
        &fund_leg,
        &pd_leg,
        &fee);

    if (err)
        return err;

    if (call->pay_rec == 0)
    {
        iv = REMNO * (pd_leg - fund_leg);
    }
    else
    {
        iv = REMNO * (fund_leg - pd_leg);
    }

    fee *= REMNO;

    if (call->fx_bound)
    {
        if (call->bar_type == 0)
        {
            if (-fee - (-iv) > 0.0)
            {
                if (fx_now <= call->barrier)
                {
                    *res = -iv;
                    return NULL;
                }
                else if (fx_now >= call->barrier + SMO)
                {
                    *res       = -fee;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier) / SMO;
                    *res = -(part * fee + (1.0 - part) * iv);
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now >= call->barrier)
                {
                    *res       = -fee;
                    *stop_path = 1;
                    return NULL;
                }
                else if (fx_now <= call->barrier - SMO)
                {
                    *res = -iv;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier + SMO) / SMO;
                    *res = -(part * fee + (1.0 - part) * iv);
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
        else
        {
            if (-fee - (-iv) > 0.0)
            {
                if (fx_now >= call->barrier)
                {
                    *res = -iv;
                    return NULL;
                }
                else if (fx_now <= call->barrier - SMO)
                {
                    *res       = -fee;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier - fx_now) / SMO;
                    *res = -(part * fee + (1.0 - part) * iv);
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now <= call->barrier)
                {
                    *res       = -fee;
                    *stop_path = 1;
                    return NULL;
                }
                else if (fx_now >= call->barrier + SMO)
                {
                    *res = -iv;
                    return NULL;
                }
                else
                {
                    part = (call->barrier + SMO - fx_now) / SMO;
                    *res = -(part * fee + (1.0 - part) * iv);
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
    }
    else
    {
        /* we knock on the IV !!!! */
        /* no smoothing */
        if (call->bar_type == 0)
        {
            if (iv <= call->barrier)
            {
                *res = -iv;
                return NULL;
            }
            else
            {
                *res       = -fee;
                *stop_path = 1;
                return NULL;
            }
        }
        else
        {
            if (iv >= call->barrier)
            {
                *res = -iv;
                return NULL;
            }
            else
            {
                *res       = -fee;
                *stop_path = 1;
                return NULL;
            }
        }
    }

    return NULL;
}

Err cpd_payoff_4_3dfx_tarn(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call;
    CPD_UND        und;
    CPD_EVAL_CONST eval_const;
    int            call_idx;
    double         Xdomfor[2];

    double partial_pd_leg;
    double fund_leg, pd_leg, iv;
    double SMO, part;

    Err err = NULL;

    *res       = 0.0;
    *stop_path = 0;

    if (!func_parm)
    {
        return NULL;
    }

    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    SMO = call->smooth;

    Xdomfor[0] = Xdom;
    Xdomfor[1] = Yfor;

    err = cpd_current_pd_cpn_val(
        Xdomfor, Zfx, spot_fx, cpd, call_idx, 1, call->pd_not_amt, eval_const, &partial_pd_leg);

    if (err)
        return err;

    CPNPAID += REMNO * partial_pd_leg;  // we add the current coupon to the coupons paid

    if (call->bar_type == 0 && CPNPAID < call->orig_barrier - SMO)  // bar_type 0: up and in
        return NULL;

    err = cpd_fut_intr_val(
        Xdomfor,
        Zfx,
        spot_fx,
        cpd,
        call->fund_idx,
        call->num_fund_cpn,
        call->fund_not_amt,
        call->pd_idx,
        call->num_pd_cpn,
        call->pd_not_amt,
        eval_const,
        &fund_leg,
        &pd_leg);

    if (err)
        return err;

    if (call->pay_rec == 0)
    {
        iv = REMNO * (pd_leg - fund_leg);
    }
    else
    {
        iv = REMNO * (fund_leg - pd_leg);
    }

    if (call->bar_type == 0)
    {
        if (iv > 0.0)
        {
            if (CPNPAID <= call->orig_barrier)
            {
                return NULL;
            }
            else if (CPNPAID >= call->orig_barrier + SMO)
            {
                *res       = iv;
                *stop_path = 1;
                return NULL;
            }
            else
            {
                part = (CPNPAID - call->orig_barrier) / SMO;
                *res = part * iv;
                REMNO *= (1.0 - part);
                return NULL;
            }
        }
        else
        {
            if (CPNPAID >= call->orig_barrier)
            {
                *res       = iv;
                *stop_path = 1;
                return NULL;
            }
            else
            {
                part = (CPNPAID - call->orig_barrier + SMO) / SMO;
                *res = part * iv;
                REMNO *= (1.0 - part);
                return NULL;
            }
        }
    }
    else
    {
        err = "TARN is an UO not DO";

        return NULL;
    }

    return NULL;
}

Err cpd_payoff_4_3dfx_mc(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call;
    CPD_UND        und;
    CPD_EVAL_CONST eval_const;
    int            call_idx;
    double         fx_now;
    double         Xdomfor[2];

    double fund_leg, pd_leg, iv;
    double SMO, part;

    Err err = NULL;

    *res       = 0.0;
    *stop_path = 0;

    if (!func_parm)
    {
        return NULL;
    }

    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    SMO = call->smooth;

    Xdomfor[0] = Xdom;
    Xdomfor[1] = Yfor;

    fx_now = Zfx + eval_const->spot_fx_dff_log_ratio - eval_const->spot_fx_for_gam * Yfor +
             eval_const->spot_fx_dom_gam * Xdom + eval_const->spot_fx_cvx;

    if (call->fx_bound && ((call->bar_type == 0 && fx_now < call->barrier - SMO) ||
                           (call->bar_type == 1 && fx_now > call->barrier + SMO)))
    {
        return NULL;
    }

    err = cpd_fut_intr_val(
        Xdomfor,
        Zfx,
        spot_fx,
        cpd,
        call->fund_idx,
        call->num_fund_cpn,
        call->fund_not_amt,
        call->pd_idx,
        call->num_pd_cpn,
        call->pd_not_amt,
        eval_const,
        &fund_leg,
        &pd_leg);

    if (err)
        return err;

    if (call->pay_rec == 0)
    {
        iv = REMNO * (pd_leg - fund_leg);
    }
    else
    {
        iv = REMNO * (fund_leg - pd_leg);
    }

    if (call->fx_bound)
    {
        if (call->bar_type == 0)
        {
            if (iv > 0.0)
            {
                if (fx_now <= call->barrier)
                {
                    return NULL;
                }
                else if (fx_now >= call->barrier + SMO)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now >= call->barrier)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier + SMO) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
        else
        {
            if (iv > 0.0)
            {
                if (fx_now >= call->barrier)
                {
                    return NULL;
                }
                else if (fx_now <= call->barrier - SMO)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier - fx_now) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now <= call->barrier)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier + SMO - fx_now) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
    }
    else
    {
        /* we knock on the IV !!!! */
        /* no smoothing */
        if (call->bar_type == 0)
        {
            if (iv <= call->barrier)
            {
                return NULL;
            }
            else
            {
                *res       = iv;
                *stop_path = 1;
                return NULL;
            }
        }
        else
        {
            if (iv >= call->barrier)
            {
                return NULL;
            }
            else
            {
                *res       = iv;
                *stop_path = 1;
                return NULL;
            }
        }
    }

    return NULL;
}

Err cpd_payoff_4_3dfx_optiFx_mc(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call;
    CPD_UND        und;
    CPD_EVAL_CONST eval_const;
    int            call_idx;
    double         fx_now;
    double         Xdomfor[2];

    double fund_leg, pd_leg, iv;
    double part, SMO;

    Err err = NULL;

    res[0]     = 0.0;
    res[1]     = 0.0;
    *stop_path = 0;

    if (!func_parm)
    {
        return NULL;
    }

    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    SMO = call->smooth;

    fx_now = Zfx + eval_const->spot_fx_dff_log_ratio - eval_const->spot_fx_for_gam * Yfor +
             eval_const->spot_fx_dom_gam * Xdom + eval_const->spot_fx_cvx;

    if (cpd->call[call_idx].call_type == 1 &&
        ((call->bar_type == 0 && fx_now < call->barrier - SMO) ||
         (call->bar_type == 1 && fx_now > call->barrier + SMO)))
    {
        return NULL;
    }

    Xdomfor[0] = Xdom;
    Xdomfor[1] = Yfor;

    err = cpd_fut_intr_val(
        Xdomfor,
        Zfx,
        spot_fx,
        cpd,
        call->fund_idx,
        call->num_fund_cpn,
        call->fund_not_amt,
        call->pd_idx,
        call->num_pd_cpn,
        call->pd_not_amt,
        eval_const,
        &fund_leg,
        &pd_leg);
    if (err)
    {
        return err;
    }

    if (call->pay_rec == 0)
    {
        iv = REMNO * (pd_leg - fund_leg);
    }
    else
    {
        iv = REMNO * (fund_leg - pd_leg);
    }

    if (cpd->call[call_idx].call_type == 1)
    {
        if (call->bar_type == 0)
        {
            if (iv > 0.0)
            {
                if (fx_now <= call->barrier)
                {
                    return NULL;
                }
                else if (fx_now >= call->barrier + SMO)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now >= call->barrier)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (fx_now - call->barrier + SMO) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
        else
        {
            if (iv > 0.0)
            {
                if (fx_now >= call->barrier)
                {
                    return NULL;
                }
                else if (fx_now <= call->barrier - SMO)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier - fx_now) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
            else
            {
                if (fx_now <= call->barrier)
                {
                    *res       = iv;
                    *stop_path = 1;
                    return NULL;
                }
                else
                {
                    part = (call->barrier + SMO - fx_now) / SMO;
                    *res = part * iv;
                    REMNO *= (1.0 - part);
                    return NULL;
                }
            }
        }
    }
    else
    {
        res[0] = iv;
        if (call->fx_bound)
        {
            res[1] = spot_fx * exp(fx_now);
        }
        else
        {
            res[1] = iv;
        }
        *stop_path = 0;
    }

    return NULL;
}

/*	Main pricing function for mc */

/*	Launch the mc */
Err cpd_launch_mc(
    CPD_STR    cpd,
    CPD_UND    und,
    CPD_MC_ARG mc_arg,
    /*	Result */
    double* prem,
    double* std)
{
    double*  temp_val_ = (double*)calloc(2, sizeof(double));
    double** temp_val  = &temp_val_;
    double   df;
    Err      err = NULL;

    if (!cpd->call[0].TARN_Do)
    {
        if (und->use_beta_dlm)
        {
            err = mc_main_3dBetaDLMfx(
                mc_arg->npaths,
                mc_arg->num_col,
                mc_arg->time,
                mc_arg->date,
                mc_arg->nb_dates,
                mc_arg->dom_fwd,
                mc_arg->dom_std,
                mc_arg->for_fwd_const,
                mc_arg->for_fwd_lin,
                mc_arg->for_std,
                mc_arg->fx_fwd,
                mc_arg->fx_std,
                mc_arg->dom_bond_pay,
                mc_arg->dom_beta_pay,
                mc_arg->dom_for_cov,
                mc_arg->dom_fx_cov,
                mc_arg->for_fx_cov,
                mc_arg->void_prm,
                mc_arg->do_pecs,
                0,
                NULL,
                NULL,
                cpd_init_4_3dfx_mc,
                cpd_payoff_total_4_3dfx_mc_dlm,
                temp_val);

            df = swp_f_df(und->today, und->model->lTStarDate, und->dom_yc);
        }
        else
        {
            err = mc_main_3dfx(
                mc_arg->npaths,
                mc_arg->num_col,
                mc_arg->time,
                mc_arg->date,
                mc_arg->nb_dates,
                mc_arg->dom_ifr,
                mc_arg->dom_fwd,
                mc_arg->dom_std,
                mc_arg->dom_phi,
                mc_arg->dom_beta,
                mc_arg->dom_bond_pay,
                mc_arg->dom_beta_pay,
                mc_arg->for_ifr,
                mc_arg->for_fwd,
                mc_arg->for_std,
                mc_arg->for_phi,
                mc_arg->for_beta,
                mc_arg->fx_fwd,
                mc_arg->fx_std,
                mc_arg->dom_for_cov,
                mc_arg->dom_fx_cov,
                mc_arg->for_fx_cov,
                mc_arg->void_prm,
                mc_arg->dom_lam,
                mc_arg->for_lam,
                mc_arg->corr_dom_for[0],
                mc_arg->corr_dom_fx[0],
                mc_arg->corr_for_fx[0],
                mc_arg->spot_fx,
                mc_arg->dom_yc,
                mc_arg->for_yc,
                mc_arg->do_pecs,
                0,
                NULL,
                NULL,
                cpd_init_4_3dfx_mc,
                cpd_payoff_4_3dfx_mc,
                temp_val,
                NULL);

            df = swp_f_df(und->today, mc_arg->date[mc_arg->nb_dates - 1] + 1.0e-16, und->dom_yc);
        }
    }
    else
    {
        if (und->use_beta_dlm)
        {
            err = "TARN not yet supported in the DLM";
            return err;

            err = mc_main_3dBetaDLMfx(
                mc_arg->npaths,
                mc_arg->num_col,
                mc_arg->time,
                mc_arg->date,
                mc_arg->nb_dates,
                mc_arg->dom_fwd,
                mc_arg->dom_std,
                mc_arg->for_fwd_const,
                mc_arg->for_fwd_lin,
                mc_arg->for_std,
                mc_arg->fx_fwd,
                mc_arg->fx_std,
                mc_arg->dom_bond_pay,
                mc_arg->dom_beta_pay,
                mc_arg->dom_for_cov,
                mc_arg->dom_fx_cov,
                mc_arg->for_fx_cov,
                mc_arg->void_prm,
                mc_arg->do_pecs,
                0,
                NULL,
                NULL,
                cpd_init_4_3dfx_mc,
                cpd_payoff_total_4_3dfx_mc_dlm,
                temp_val);

            df = swp_f_df(und->today, und->model->lTStarDate, und->dom_yc);
        }
        else
        {
            err = mc_main_3dfx(
                mc_arg->npaths,
                mc_arg->num_col,
                mc_arg->time,
                mc_arg->date,
                mc_arg->nb_dates,
                mc_arg->dom_ifr,
                mc_arg->dom_fwd,
                mc_arg->dom_std,
                mc_arg->dom_phi,
                mc_arg->dom_beta,
                mc_arg->dom_bond_pay,
                mc_arg->dom_beta_pay,
                mc_arg->for_ifr,
                mc_arg->for_fwd,
                mc_arg->for_std,
                mc_arg->for_phi,
                mc_arg->for_beta,
                mc_arg->fx_fwd,
                mc_arg->fx_std,
                mc_arg->dom_for_cov,
                mc_arg->dom_fx_cov,
                mc_arg->for_fx_cov,
                mc_arg->void_prm,
                mc_arg->dom_lam,
                mc_arg->for_lam,
                mc_arg->corr_dom_for[0],
                mc_arg->corr_dom_fx[0],
                mc_arg->corr_for_fx[0],
                mc_arg->spot_fx,
                mc_arg->dom_yc,
                mc_arg->for_yc,
                mc_arg->do_pecs,
                0,
                NULL,
                NULL,
                cpd_init_4_3dfx_mc,
                cpd_payoff_4_3dfx_tarn,
                temp_val,
                NULL);

            df = swp_f_df(und->today, mc_arg->date[mc_arg->nb_dates - 1] + 1.0e-16, und->dom_yc);
        }
    }
    *prem = temp_val[0][0] * df;
    *std  = temp_val[0][1] * df;

    free(temp_val_);

    return err;
}

/*	Main pricing function for mc */

/*	Launch the mc */
Err cpd_launch_opti_mc(
    CPD_STR    cpd,
    CPD_UND    und,
    CPD_MC_ARG mc_arg,
    /*	Result */
    double*   prem,
    double*   std,
    int       do_infos,
    double*** optim_bar)
{
    double**   temp_val = dmatrix(0, 2, 0, 2);
    int*       optimise = (int*)calloc(mc_arg->nb_dates, sizeof(int));
    MCEBParams params_, *params = &params_;

    double df;
    int    i, j;
    Err    err = NULL;

    if (cpd->call[0].ex_time < 1.0E-10)
    {
        for (i = 0; i < mc_arg->nb_dates; i++)
        {
            optimise[i] = 1 - cpd->call[i].call_type;
        }
    }
    else
    {
        optimise[0] = 0;
        for (i = 1; i < mc_arg->nb_dates; i++)
        {
            optimise[i] = 1 - cpd->call[i - 1].call_type;
        }
    }

    mceb_set_default_params(params);

    params->iCallCurrent          = 1;
    params->iIsKO                 = 0;
    params->iDoInfos              = do_infos;
    params->iMultiIndex           = 0;
    params->iColPay               = 0;
    params->iColBound             = 1;
    params->iNbIndex              = 1;
    params->iKnockInCol           = 0;
    params->iAddNonOptimisedForKI = 1;

    *optim_bar = (double**)dmatrix(
        0, mc_arg->nb_dates - 1, 0, do_infos * (3 + 2 * mc_arg->nb_dates) + (1 - do_infos) * 1);

    err = mceb_allocate_params(params, mc_arg->nb_dates);
    if (err)
        goto FREE_RETURN;

    err = mc_main_3dfx(
        mc_arg->npaths,
        2,
        mc_arg->time,
        mc_arg->date,
        mc_arg->nb_dates,
        mc_arg->dom_ifr,
        mc_arg->dom_fwd,
        mc_arg->dom_std,
        mc_arg->dom_phi,
        mc_arg->dom_beta,
        mc_arg->dom_bond_pay,
        mc_arg->dom_beta_pay,
        mc_arg->for_ifr,
        mc_arg->for_fwd,
        mc_arg->for_std,
        mc_arg->for_phi,
        mc_arg->for_beta,
        mc_arg->fx_fwd,
        mc_arg->fx_std,
        mc_arg->dom_for_cov,
        mc_arg->dom_fx_cov,
        mc_arg->for_fx_cov,
        mc_arg->void_prm,
        mc_arg->dom_lam,
        mc_arg->for_lam,
        mc_arg->corr_dom_for[0],
        mc_arg->corr_dom_fx[0],
        mc_arg->corr_for_fx[0],
        mc_arg->spot_fx,
        mc_arg->dom_yc,
        mc_arg->for_yc,
        mc_arg->do_pecs,
        1,
        optimise,
        params,
        cpd_init_4_3dfx_mc,
        cpd_payoff_4_3dfx_optiFx_mc,
        temp_val,
        *optim_bar);

    if (err)
        goto FREE_RETURN;

    /* Recopy Barrier / CoefLin for the moment */
    for (i = 0; i < mc_arg->nb_dates; i++)
    {
        (*optim_bar)[i][0] = params->dBarrier[i];

        for (j = 0; j < params->iNbIndex; j++)
        {
            (*optim_bar)[i][1 + j] = params->dCoefLin[i][j + 1];
        }
    }

    df = swp_f_df(und->today, mc_arg->date[mc_arg->nb_dates - 1] + 1.0e-16, und->dom_yc);

    *prem = temp_val[2][0] * df;
    *std  = temp_val[2][1] * df;

    if (do_infos)
    {
        for (i = 0; i < mc_arg->nb_dates; i++)
        {
            (*optim_bar)[i][1] *= df;
            (*optim_bar)[i][3] *= df;

            for (j = 0; j < mc_arg->nb_dates; j++)
            {
                (*optim_bar)[i][4 + mc_arg->nb_dates + j] *= df;
            }
        }
    }

FREE_RETURN:

    free_dmatrix(temp_val, 0, max(mc_arg->nb_dates - 1, 2), 0, 2);
    free(optimise);

    mceb_free_params(params);

    return err;
}

Err cpd_payoff_4_3dfx_BetaDLM_tree(
    double evt_date,
    double evt_time,
    void*  func_parm,
    double spot_fx,
    void*  dom_yc,
    void*  for_yc,
    long   n1,
    long   n2,
    long   n3,
    /* i: d1, j: d2, k: d3, l = {0: xDom, 1: xFor, 2: log (Fx/Fx0)} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl)
{
    CPD_PAY_ARG    cpd_arg;
    CPD_STR        cpd;
    PD_CALL        call, next_call;
    CPD_UND        und;
    PD_EXO_CPN     cpn;
    CPD_EVAL_CONST eval_const;
    int            call_idx;

    int i, j, k, l;

    double ***svi, **svij, *svijk;
    double    X, fx;

    double fund_leg, df, log_fwd, fwd, fwd_adj, floor, cap, pd_leg, iv, fee;

    int num_fund_cpn, num_pd_cpn, eval_redemption;

    double prices[3];

    Err err = NULL;

    /* For Optimisation */
    double X_I, X_J, X_K, X_DI, X_DJ, X_DK;

    double fx_I, fx_J, fx_K, fx_DI_Const, fx_DJ_Const, fx_DK_Const, fx_DI_Lin, fx_DJ_Lin, fx_DK_Lin;

    double start_I, start_J, start_K, start_DI, start_DJ, start_DK;

    double in_not_fund_I, in_not_fund_J, in_not_fund_K, in_not_fund_DI, in_not_fund_DJ,
        in_not_fund_DK;

    double in_not_pd_I, in_not_pd_J, in_not_pd_K, in_not_pd_DI, in_not_pd_DJ, in_not_pd_DK;

    double next_start_fund_I, next_start_fund_J, next_start_fund_K, next_start_fund_DI,
        next_start_fund_DJ, next_start_fund_DK;

    double next_fund_I, next_fund_J, next_fund_K, next_fund_DI, next_fund_DJ, next_fund_DK;

    double next_pd_I, next_pd_J, next_pd_K, next_pd_DI, next_pd_DJ, next_pd_DK;

    double fin_not_disc_I, fin_not_disc_J, fin_not_disc_K, fin_not_disc_DI, fin_not_disc_DJ,
        fin_not_disc_DK;

    double fin_not_log_fwd_I, fin_not_log_fwd_J, fin_not_log_fwd_K, fin_not_log_fwd_DI,
        fin_not_log_fwd_DJ, fin_not_log_fwd_DK;

    double fin_not_fwd_I, fin_not_fwd_J, fin_not_fwd_K, fin_not_fwd_DI, fin_not_fwd_DJ,
        fin_not_fwd_DK;

    double maxIV = 0, maxProdVal = 0;

    /*	Get the event */
    cpd_arg    = (CPD_PAY_ARG)func_parm;
    eval_const = (CPD_EVAL_CONST)(&(cpd_arg->eval_const));
    cpd        = cpd_arg->cpd;
    call_idx   = cpd_arg->call_idx;
    call       = cpd->call + call_idx;
    und        = (CPD_UND)(cpd_arg->und);

    /*	Calculate number of coupons to be computed */
    if (call_idx == cpd->num_calls - 1)
    {
        /*	Last call: eval all + redemption */
        num_fund_cpn    = call->num_fund_cpn;
        num_pd_cpn      = call->num_pd_cpn;
        eval_redemption = 1;
    }
    else
    {
        /*	Not the last call: only eval coupons up to next call date */
        next_call       = cpd->call + (call_idx + 1);
        num_fund_cpn    = next_call->fund_idx - call->fund_idx;
        num_pd_cpn      = next_call->pd_idx - call->pd_idx;
        eval_redemption = 0;
    }

    /*	Eval payoff
            ----------- */

    /* ---------------------------------------------------------------------------------------------------------------------------
            For Optimisation
       ---------------------------------------------------------------------------------------------------------------------------
     */

    /* For X */
    X_I = eval_const->pd_X_const + eval_const->pd_X_dom * sv[0][0][0][0] +
          eval_const->pd_X_for * sv[0][0][0][1] + sv[0][0][0][2];
    X_J = X_I;
    X_K = X_I;
    if (n1 > 0)
    {
        X_DI = eval_const->pd_X_dom * (sv[1][0][0][0] - sv[0][0][0][0]) +
               eval_const->pd_X_for * (sv[1][0][0][1] - sv[0][0][0][1]) +
               (sv[1][0][0][2] - sv[0][0][0][2]);
    }
    else
    {
        X_DI = 0.0;
    }

    if (n2 > 0)
    {
        X_DJ = eval_const->pd_X_dom * (sv[0][1][0][0] - sv[0][0][0][0]) +
               eval_const->pd_X_for * (sv[0][1][0][1] - sv[0][0][0][1]) +
               (sv[0][1][0][2] - sv[0][0][0][2]);
    }
    else
    {
        X_DJ = 0.0;
    }

    if (n3 > 0)
    {
        X_DK = eval_const->pd_X_dom * (sv[0][0][1][0] - sv[0][0][0][0]) +
               eval_const->pd_X_for * (sv[0][0][1][1] - sv[0][0][0][1]) +
               (sv[0][0][1][2] - sv[0][0][0][2]);
    }
    else
    {
        X_DK = 0.0;
    }

    /* For Fx */
    fx_I = eval_const->pd_S_const + X_I * (eval_const->pd_S_lin + X_I * eval_const->pd_S_quad) +
           eval_const->pd_S_dom * sv[0][0][0][0] + eval_const->pd_S_for * sv[0][0][0][1];
    fx_J = fx_I;
    fx_K = fx_I;

    if (n1 > 0)
    {
        fx_DI_Const = X_DI * (eval_const->pd_S_lin + X_DI * eval_const->pd_S_quad) +
                      eval_const->pd_S_dom * (sv[1][0][0][0] - sv[0][0][0][0]) +
                      eval_const->pd_S_for * (sv[1][0][0][1] - sv[0][0][0][1]);
        fx_DI_Lin = 2.0 * X_DI * eval_const->pd_S_quad;
    }
    else
    {
        fx_DI_Const = 0.0;
    }

    if (n2 > 0)
    {
        fx_DJ_Const = X_DJ * (eval_const->pd_S_lin + X_DJ * eval_const->pd_S_quad) +
                      eval_const->pd_S_dom * (sv[0][1][0][0] - sv[0][0][0][0]) +
                      eval_const->pd_S_for * (sv[0][1][0][1] - sv[0][0][0][1]);
        fx_DJ_Lin = 2.0 * X_DJ * eval_const->pd_S_quad;
    }
    else
    {
        fx_DJ_Const = 0.0;
    }

    if (n3 > 0)
    {
        fx_DK_Const = X_DK * (eval_const->pd_S_lin + X_DK * eval_const->pd_S_quad) +
                      eval_const->pd_S_dom * (sv[0][0][1][0] - sv[0][0][0][0]) +
                      eval_const->pd_S_for * (sv[0][0][1][1] - sv[0][0][0][1]);
        fx_DK_Lin = 2.0 * X_DK * eval_const->pd_S_quad;
    }
    else
    {
        fx_DK_Const = 0.0;
    }

    /* for funding start */
    start_I =
        eval_const->start_dff *
        exp(-eval_const->start_gam * sv[0][0][0][eval_const->fund_var] - eval_const->start_gam2) *
        cpd->fund_leg->notional;

    start_J = start_I;
    start_K = start_I;

    if (n1 > 0)
    {
        start_DI =
            exp(-eval_const->start_gam *
                (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    }
    else
    {
        start_DI = 1.0;
    }

    if (n2 > 0)
    {
        start_DJ =
            exp(-eval_const->start_gam *
                (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    }
    else
    {
        start_DJ = 1.0;
    }

    if (n3 > 0)
    {
        start_DK =
            exp(-eval_const->start_gam *
                (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    }
    else
    {
        start_DK = 1.0;
    }

    /* For initial notional */
    in_not_fund_I = eval_const->in_not_fund_dff *
                    exp(-eval_const->in_not_fund_gam * sv[0][0][0][eval_const->fund_var] -
                        eval_const->in_not_fund_gam2) *
                    call->fund_not_amt;

    in_not_fund_J = in_not_fund_I;
    in_not_fund_K = in_not_fund_I;

    if (n1 > 0)
    {
        in_not_fund_DI =
            exp(-eval_const->in_not_fund_gam *
                (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    }
    else
    {
        in_not_fund_DI = 1.0;
    }

    if (n2 > 0)
    {
        in_not_fund_DJ =
            exp(-eval_const->in_not_fund_gam *
                (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    }
    else
    {
        in_not_fund_DJ = 1.0;
    }

    if (n3 > 0)
    {
        in_not_fund_DK =
            exp(-eval_const->in_not_fund_gam *
                (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
    }
    else
    {
        in_not_fund_DK = 1.0;
    }

    /*	for fee */
    in_not_pd_I = eval_const->in_not_pd_dff *
                  exp(-eval_const->in_not_pd_gam * sv[0][0][0][0] - eval_const->in_not_pd_gam2);

    in_not_pd_J = in_not_pd_I;
    in_not_pd_K = in_not_pd_I;

    if (n1 > 0)
    {
        in_not_pd_DI = exp(-eval_const->in_not_pd_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
    }
    else
    {
        in_not_pd_DI = 1.0;
    }

    if (n2 > 0)
    {
        in_not_pd_DJ = exp(-eval_const->in_not_pd_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
    }
    else
    {
        in_not_pd_DJ = 1.0;
    }

    if (n3 > 0)
    {
        in_not_pd_DK = exp(-eval_const->in_not_pd_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
    }
    else
    {
        in_not_pd_DK = 1.0;
    }

    /*	Next exchange if relevant */
    if (!eval_redemption)
    {
        next_start_fund_I =
            eval_const->next_start_fund_dff *
            exp(-eval_const->next_start_fund_gam * sv[0][0][0][eval_const->fund_var] -
                eval_const->next_start_fund_gam2) *
            cpd->fund_leg->notional;

        next_start_fund_J = next_start_fund_I;
        next_start_fund_K = next_start_fund_I;

        if (n1 > 0)
        {
            next_start_fund_DI =
                exp(-eval_const->next_start_fund_gam *
                    (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            next_start_fund_DI = 1.0;
        }

        if (n2 > 0)
        {
            next_start_fund_DJ =
                exp(-eval_const->next_start_fund_gam *
                    (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            next_start_fund_DJ = 1.0;
        }

        if (n3 > 0)
        {
            next_start_fund_DK =
                exp(-eval_const->next_start_fund_gam *
                    (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            next_start_fund_DK = 1.0;
        }

        next_fund_I = eval_const->next_fund_dff *
                      exp(-eval_const->next_fund_gam * sv[0][0][0][eval_const->fund_var] -
                          eval_const->next_fund_gam2) *
                      next_call->fund_not_amt;

        next_fund_J = next_fund_I;
        next_fund_K = next_fund_I;

        if (n1 > 0)
        {
            next_fund_DI =
                exp(-eval_const->next_fund_gam *
                    (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            next_fund_DI = 1.0;
        }

        if (n2 > 0)
        {
            next_fund_DJ =
                exp(-eval_const->next_fund_gam *
                    (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            next_fund_DJ = 1.0;
        }

        if (n3 > 0)
        {
            next_fund_DK =
                exp(-eval_const->next_fund_gam *
                    (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            next_fund_DK = 1.0;
        }

        next_pd_I = eval_const->next_pd_dff *
                    exp(-eval_const->next_pd_gam * sv[0][0][0][0] - eval_const->next_pd_gam2) *
                    next_call->pd_not_amt;

        next_pd_J = next_pd_I;
        next_pd_K = next_pd_I;

        if (n1 > 0)
        {
            next_pd_DI = exp(-eval_const->next_pd_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        }
        else
        {
            next_pd_DI = 1.0;
        }

        if (n2 > 0)
        {
            next_pd_DJ = exp(-eval_const->next_pd_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        }
        else
        {
            next_pd_DJ = 1.0;
        }

        if (n3 > 0)
        {
            next_pd_DK = exp(-eval_const->next_pd_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        }
        else
        {
            next_pd_DK = 1.0;
        }
    }
    else
    /*	Final notional */
    {
        /* For df */
        fin_not_disc_I =
            eval_const->fin_not_disc_dff *
            exp(-eval_const->fin_not_disc_gam * sv[0][0][0][0] - eval_const->fin_not_disc_gam2);

        fin_not_disc_J = fin_not_disc_I;
        fin_not_disc_K = fin_not_disc_I;

        if (n1 > 0)
        {
            fin_not_disc_DI =
                exp(-eval_const->fin_not_disc_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        }
        else
        {
            fin_not_disc_DI = 1.0;
        }

        if (n2 > 0)
        {
            fin_not_disc_DJ =
                exp(-eval_const->fin_not_disc_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        }
        else
        {
            fin_not_disc_DJ = 1.0;
        }

        if (n3 > 0)
        {
            fin_not_disc_DK =
                exp(-eval_const->fin_not_disc_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        }
        else
        {
            fin_not_disc_DK = 1.0;
        }

        /*	For log Fwd fx */
        fin_not_log_fwd_I = log(eval_const->fin_not_fwd_dff_ratio) -
                            eval_const->fin_not_fwd_for_gam * sv[0][0][0][1] +
                            eval_const->fin_not_fwd_dom_gam * sv[0][0][0][0] +
                            eval_const->fin_not_fwd_cvx;

        fin_not_log_fwd_J = fin_not_log_fwd_I;
        fin_not_log_fwd_K = fin_not_log_fwd_I;

        if (n1 > 0)
        {
            fin_not_log_fwd_DI =
                -eval_const->fin_not_fwd_for_gam * (sv[1][0][0][1] - sv[0][0][0][1]) +
                eval_const->fin_not_fwd_dom_gam * (sv[1][0][0][0] - sv[0][0][0][0]);
        }
        else
        {
            fin_not_log_fwd_DI = 0.0;
        }

        if (n2 > 0)
        {
            fin_not_log_fwd_DJ =
                -eval_const->fin_not_fwd_for_gam * (sv[0][1][0][1] - sv[0][0][0][1]) +
                eval_const->fin_not_fwd_dom_gam * (sv[0][1][0][0] - sv[0][0][0][0]);
        }
        else
        {
            fin_not_log_fwd_DJ = 0.0;
        }

        if (n3 > 0)
        {
            fin_not_log_fwd_DK =
                -eval_const->fin_not_fwd_for_gam * (sv[0][0][1][1] - sv[0][0][0][1]) +
                eval_const->fin_not_fwd_dom_gam * (sv[0][0][1][0] - sv[0][0][0][0]);
        }
        else
        {
            fin_not_log_fwd_DK = 0.0;
        }

        /*	For Fwd fx */
        fin_not_fwd_I =
            eval_const->fin_not_fwd_dff_ratio *
            exp(-eval_const->fin_not_fwd_for_gam * sv[0][0][0][1] +
                eval_const->fin_not_fwd_dom_gam * sv[0][0][0][0] + eval_const->fin_not_fwd_cvx);

        fin_not_fwd_J = fin_not_fwd_I;
        fin_not_fwd_K = fin_not_fwd_I;

        if (n1 > 0)
        {
            fin_not_fwd_DI =
                exp(-eval_const->fin_not_fwd_for_gam * (sv[1][0][0][1] - sv[0][0][0][1]) +
                    eval_const->fin_not_fwd_dom_gam * (sv[1][0][0][0] - sv[0][0][0][0]));
        }
        else
        {
            fin_not_fwd_DI = 1.0;
        }

        if (n2 > 0)
        {
            fin_not_fwd_DJ =
                exp(-eval_const->fin_not_fwd_for_gam * (sv[0][1][0][1] - sv[0][0][0][1]) +
                    eval_const->fin_not_fwd_dom_gam * (sv[0][1][0][0] - sv[0][0][0][0]));
        }
        else
        {
            fin_not_fwd_DJ = 1.0;
        }

        if (n3 > 0)
        {
            fin_not_fwd_DK =
                exp(-eval_const->fin_not_fwd_for_gam * (sv[0][0][1][1] - sv[0][0][0][1]) +
                    eval_const->fin_not_fwd_dom_gam * (sv[0][0][1][0] - sv[0][0][0][0]));
        }
        else
        {
            fin_not_fwd_DK = 1.0;
        }
    }

    for (l = 0; l < num_fund_cpn; l++)
    {
        /* for funding coupon */
        eval_const->fund_I[l] = eval_const->fund_dff[l] *
                                exp(-eval_const->fund_gam[l] * sv[0][0][0][eval_const->fund_var] -
                                    eval_const->fund_gam2[l]) *
                                cpd->fund_leg->cpn[call->fund_idx + l].cpn;

        eval_const->fund_J[l] = eval_const->fund_I[l];

        eval_const->fund_K[l] = eval_const->fund_I[l];

        if (n1 > 0)
        {
            eval_const->fund_DI[l] =
                exp(-eval_const->fund_gam[l] *
                    (sv[1][0][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            eval_const->fund_DI[l] = 1.0;
        }

        if (n2 > 0)
        {
            eval_const->fund_DJ[l] =
                exp(-eval_const->fund_gam[l] *
                    (sv[0][1][0][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            eval_const->fund_DJ[l] = 1.0;
        }

        if (n3 > 0)
        {
            eval_const->fund_DK[l] =
                exp(-eval_const->fund_gam[l] *
                    (sv[0][0][1][eval_const->fund_var] - sv[0][0][0][eval_const->fund_var]));
        }
        else
        {
            eval_const->fund_DK[l] = 1.0;
        }
    }

    for (l = 0; l < num_pd_cpn; l++)
    {
        /*	Coupons */
        eval_const->pd_disc_I[l] =
            eval_const->pd_disc_dff[l] *
            exp(-eval_const->pd_disc_gam[l] * sv[0][0][0][0] - eval_const->pd_disc_gam2[l]);

        eval_const->pd_disc_J[l] = eval_const->pd_disc_I[l];

        eval_const->pd_disc_K[l] = eval_const->pd_disc_I[l];

        if (n1 > 0)
        {
            eval_const->pd_disc_DI[l] =
                exp(-eval_const->pd_disc_gam[l] * (sv[1][0][0][0] - sv[0][0][0][0]));
        }
        else
        {
            eval_const->pd_disc_DI[l] = 1.0;
        }

        if (n2 > 0)
        {
            eval_const->pd_disc_DJ[l] =
                exp(-eval_const->pd_disc_gam[l] * (sv[0][1][0][0] - sv[0][0][0][0]));
        }
        else
        {
            eval_const->pd_disc_DJ[l] = 1.0;
        }

        if (n3 > 0)
        {
            eval_const->pd_disc_DK[l] =
                exp(-eval_const->pd_disc_gam[l] * (sv[0][0][1][0] - sv[0][0][0][0]));
        }
        else
        {
            eval_const->pd_disc_DK[l] = 1.0;
        }

        /* for the log fwd */
        eval_const->pd_log_fwd_I[l] =
            log(eval_const->pd_fwd_dff_ratio[l]) - eval_const->pd_fwd_for_gam[l] * sv[0][0][0][1] +
            eval_const->pd_fwd_dom_gam[l] * sv[0][0][0][0] + eval_const->pd_fwd_cvx[l];

        eval_const->pd_log_fwd_J[l] = eval_const->pd_log_fwd_I[l];

        eval_const->pd_log_fwd_K[l] = eval_const->pd_log_fwd_I[l];

        if (n1 > 0)
        {
            eval_const->pd_log_fwd_DI[l] =
                -eval_const->pd_fwd_for_gam[l] * (sv[1][0][0][1] - sv[0][0][0][1]) +
                eval_const->pd_fwd_dom_gam[l] * (sv[1][0][0][0] - sv[0][0][0][0]);
        }
        else
        {
            eval_const->pd_log_fwd_DI[l] = 0.0;
        }

        if (n2 > 0)
        {
            eval_const->pd_log_fwd_DJ[l] =
                -eval_const->pd_fwd_for_gam[l] * (sv[0][1][0][1] - sv[0][0][0][1]) +
                eval_const->pd_fwd_dom_gam[l] * (sv[0][1][0][0] - sv[0][0][0][0]);
        }
        else
        {
            eval_const->pd_log_fwd_DJ[l] = 0.0;
        }

        if (n3 > 0)
        {
            eval_const->pd_log_fwd_DK[l] =
                -eval_const->pd_fwd_for_gam[l] * (sv[0][0][1][1] - sv[0][0][0][1]) +
                eval_const->pd_fwd_dom_gam[l] * (sv[0][0][1][0] - sv[0][0][0][0]);
        }
        else
        {
            eval_const->pd_log_fwd_DK[l] = 0.0;
        }

        /* for the fwd */

        eval_const->pd_fwd_I[l] =
            eval_const->pd_fwd_dff_ratio[l] *
            exp(-eval_const->pd_fwd_for_gam[l] * sv[0][0][0][1] +
                eval_const->pd_fwd_dom_gam[l] * sv[0][0][0][0] + eval_const->pd_fwd_cvx[l]);

        eval_const->pd_fwd_J[l] = eval_const->pd_fwd_I[l];

        eval_const->pd_fwd_K[l] = eval_const->pd_fwd_I[l];

        if (n1 > 0)
        {
            eval_const->pd_fwd_DI[l] =
                exp(-eval_const->pd_fwd_for_gam[l] * (sv[1][0][0][1] - sv[0][0][0][1]) +
                    eval_const->pd_fwd_dom_gam[l] * (sv[1][0][0][0] - sv[0][0][0][0]));
        }
        else
        {
            eval_const->pd_fwd_DI[l] = 1.0;
        }

        if (n2 > 0)
        {
            eval_const->pd_fwd_DJ[l] =
                exp(-eval_const->pd_fwd_for_gam[l] * (sv[0][1][0][1] - sv[0][0][0][1]) +
                    eval_const->pd_fwd_dom_gam[l] * (sv[0][1][0][0] - sv[0][0][0][0]));
        }
        else
        {
            eval_const->pd_fwd_DJ[l] = 1.0;
        }

        if (n3 > 0)
        {
            eval_const->pd_fwd_DK[l] =
                exp(-eval_const->pd_fwd_for_gam[l] * (sv[0][0][1][1] - sv[0][0][0][1]) +
                    eval_const->pd_fwd_dom_gam[l] * (sv[0][0][1][0] - sv[0][0][0][0]));
        }
        else
        {
            eval_const->pd_fwd_DK[l] = 1.0;
        }
    }

    /* ---------------------------------------------------------------------------------------------------------------------------
            End of init of the Optimisation
       ---------------------------------------------------------------------------------------------------------------------------
     */

    for (i = 0; i < n1; i++)
    {
        svi = sv[i];

        for (j = 0; j < n2; j++)
        {
            svij = svi[j];

            for (k = 0; k < n3; k++)
            {
                svijk = svij[k];

                // X = eval_const->pd_X_const + eval_const->pd_X_dom * svijk[0] +
                // eval_const->pd_X_for * svijk[1] + svijk[2];
                X = X_K;

                /*fx = exp (eval_const->pd_S_const + X * (eval_const->pd_S_lin + X *
                   eval_const->pd_S_quad)
                        + eval_const->pd_S_dom * svijk[0] + eval_const->pd_S_for * svijk[1]);*/
                fx = exp(fx_K);

                /*	PV of funding leg */
                fund_leg = 0.0;

                /*	Libor */

                /*fund_leg += eval_const->start_dff
                 * exp (- eval_const->start_gam * svijk[eval_const->fund_var] -
                 * eval_const->start_gam2) cpd->fund_leg->notional; */

                fund_leg += start_K;

                /*	Coupons */

                for (l = 0; l < num_fund_cpn; l++)
                {
                    /*fund_leg +=
                            eval_const->fund_dff[l]
            * exp (- eval_const->fund_gam[l] * svijk[eval_const->fund_var] - eval_const->fund_gam2[l])
                            * cpd->fund_leg->cpn[call->fund_idx+l].cpn; */
                    fund_leg += eval_const->fund_K[l];
                }

                /*	Initial notional */

                /*fund_leg -= eval_const->in_not_fund_dff
                 * exp (- eval_const->in_not_fund_gam * svijk[eval_const->fund_var] -
                 * eval_const->in_not_fund_gam2) call->fund_not_amt; */

                fund_leg -= in_not_fund_K;

                /*	Next exchange if relevant */
                if (!eval_redemption)
                {
                    /* fund_leg -= eval_const->next_start_fund_dff
                     * exp (- eval_const->next_start_fund_gam * svijk[eval_const->fund_var] -
                     * eval_const->next_start_fund_gam2) cpd->fund_leg->notional; */
                    fund_leg -= next_start_fund_K;

                    /*fund_leg += eval_const->next_fund_dff
                     * exp (- eval_const->next_fund_gam * svijk[eval_const->fund_var] -
                     * eval_const->next_fund_gam2) next_call->fund_not_amt; */
                    fund_leg += next_fund_K;
                }

                /*	Fx */

                if (cpd->fund_leg->dom_for != 0)
                {
                    fund_leg *= fx;
                }

                /*	PV of pd leg */

                pd_leg = 0.0;

                /*	Coupons */
                for (l = 0; l < num_pd_cpn; l++)
                {
                    /*	Coupon access */
                    cpn = cpd->pd_leg->cpn + (call->pd_idx + l);

                    /*	Discount */
                    /*df =
                            eval_const->pd_disc_dff[l]
                            * exp (- eval_const->pd_disc_gam[l] * svijk[0] - eval_const->pd_disc_gam2[l]);
                    */
                    df = eval_const->pd_disc_K[l];

                    /*	Fwd fx */
                    /*fwd = fx
                            * eval_const->pd_fwd_dff_ratio[l]
                            * exp (	-	eval_const->pd_fwd_for_gam[l] * svijk[1]
                                            +	eval_const->pd_fwd_dom_gam[l] * svijk[0]
                                            +	eval_const->pd_fwd_cvx[l]); */
                    log_fwd = fx_K + eval_const->pd_log_fwd_K[l];
                    fwd     = fx * eval_const->pd_fwd_K[l];

                    /*	Floor and Cap */

                    if (eval_const->pd_inst[l].iNbStrike > 0)
                    {
                        err = cpd_price_FxInstrument_dlm(
                            fwd, log_fwd, X, und, cpd, l, cpn, eval_const, prices);

                        if (err)
                            return err;
                    }

                    /* get the fwd */
                    if (eval_const->index_fwd[l] >= 0)
                    {
                        fwd_adj = prices[eval_const->index_fwd[l]];
                    }
                    else
                    {
                        fwd_adj = fwd;
                    }

                    switch (eval_const->pd_floor_type[l])
                    {
                    case 0:
                        floor = 0.0;
                        break;
                    case 1:
                    case 2:
                        floor = OPT_VAL_MACRO(
                            eval_const->pd_floor_type[l],
                            fwd_adj,
                            eval_const->pd_floor_str[l],
                            1.0,
                            1.0);
                        break;
                    case 3:
                        floor = prices[0];
                        break;
                    case 4:
                        floor = prices[0] + eval_const->pd_floor_str[l] - fwd_adj;
                        break;
                    }

                    switch (eval_const->pd_cap_type[l])
                    {
                    case 0:
                        cap = 0.0;
                        break;
                    case 1:
                    case 2:
                        cap = OPT_VAL_MACRO(
                            eval_const->pd_cap_type[l],
                            fwd_adj,
                            eval_const->pd_cap_str[l],
                            1.0,
                            1.0);
                        break;
                    case 3:
                        cap = prices[eval_const->index_cap[l]];
                        break;
                    case 4:
                        cap =
                            prices[eval_const->index_cap[l]] + eval_const->pd_cap_str[l] - fwd_adj;
                        break;
                    }

                    /*	Coupon pv */

                    pd_leg += df * (cpn->alpha + cpn->beta * fwd_adj +
                                    eval_const->pd_abs_beta[l] * (floor - cap));
                }
                /*	Final notional */

                if (eval_redemption)
                {
                    /*	Coupon access */
                    cpn = &(cpd->pd_leg->not_ref);

                    /*	Discount */
                    /*df =
                            eval_const->fin_not_disc_dff
                            * exp (- eval_const->fin_not_disc_gam * svijk[0] - eval_const->fin_not_disc_gam2); */
                    df = fin_not_disc_K;

                    /*	Fwd fx */
                    /* fwd = fx
                            * eval_const->fin_not_fwd_dff_ratio
                            * exp (	-	eval_const->fin_not_fwd_for_gam * svijk[1]
                                            +	eval_const->fin_not_fwd_dom_gam * svijk[0]
                                            +	eval_const->fin_not_fwd_cvx); */
                    log_fwd = fx_K + fin_not_log_fwd_K;
                    fwd     = fx * fin_not_fwd_K;

                    if (eval_const->fin_not_inst.iNbStrike > 0)
                    {
                        err = cpd_price_FxInstrument_dlm(
                            fwd, log_fwd, X, und, cpd, -1, cpn, eval_const, prices);

                        if (err)
                            return err;
                    }

                    /* get the fwd */
                    if (eval_const->fin_not_index_fwd >= 0)
                    {
                        fwd_adj = prices[eval_const->fin_not_index_fwd];
                    }
                    else
                    {
                        fwd_adj = fwd;
                    }

                    switch (eval_const->fin_not_floor_type)
                    {
                    case 0:
                        floor = 0.0;
                        break;
                    case 1:
                    case 2:
                        floor = OPT_VAL_MACRO(
                            eval_const->fin_not_floor_type,
                            fwd_adj,
                            eval_const->fin_not_floor_str,
                            1.0,
                            1.0);
                        break;
                    case 3:
                        floor = prices[0];
                        break;
                    case 4:
                        floor = prices[0] + eval_const->fin_not_floor_str - fwd_adj;
                        break;
                    }

                    switch (eval_const->fin_not_cap_type)
                    {
                    case 0:
                        cap = 0.0;
                        break;
                    case 1:
                    case 2:
                        cap = OPT_VAL_MACRO(
                            eval_const->fin_not_cap_type,
                            fwd_adj,
                            eval_const->fin_not_cap_str,
                            1.0,
                            1.0);
                        break;
                    case 3:
                        cap = prices[eval_const->fin_not_index_cap];
                        break;
                    case 4:
                        cap = prices[eval_const->fin_not_index_cap] + eval_const->fin_not_cap_str -
                              fwd_adj;
                        break;
                    }

                    /*	Coupon pv */

                    pd_leg += df * (cpn->alpha + cpn->beta * fwd +
                                    eval_const->fin_not_abs_beta * (floor - cap));
                }
                else
                {
                    /* pd_leg += eval_const->next_pd_dff
                     * exp (- eval_const->next_pd_gam * svijk[0] - eval_const->next_pd_gam2)
                     * next_call->pd_not_amt; */
                    pd_leg += next_pd_K;
                }

                /*	Initial notional */

                /* fee = eval_const->in_not_pd_dff
                 * exp (- eval_const->in_not_pd_gam * svijk[0] - eval_const->in_not_pd_gam2); */
                fee = in_not_pd_K;

                pd_leg -= fee * call->pd_not_amt;

                fee *= call->fee;

                /*	Finally, intrinsic value */

                if (call->pay_rec == 0)
                {
                    iv = pd_leg - fund_leg;
                }
                else
                {
                    iv = fund_leg - pd_leg;
                }

                /*	Process max */

                if (prod_val[i][j][k][0] - iv + fee > 0)
                {
                    prod_val[i][j][k][0] -= iv;
                }
                else
                {
                    /* we pay the fee and call */
                    prod_val[i][j][k][0] = -fee;
                }

                /* Get the IV priced in the tree */
                if (nprod == 2)
                {
                    if (iv > maxIV)
                        maxIV = iv;
                    if (fabs(prod_val[i][j][k][1]) > maxProdVal)
                        maxProdVal = fabs(prod_val[i][j][k][1]);

                    prod_val[i][j][k][1] += iv;
                }

                /* For Optimisation */
                fx_K += fx_DK_Const + X_K * fx_DK_Lin;
                X_K += X_DK;
                start_K *= start_DK;
                in_not_fund_K *= in_not_fund_DK;
                in_not_pd_K *= in_not_pd_DK;

                if (!eval_redemption)
                {
                    next_start_fund_K *= next_start_fund_DK;
                    next_fund_K *= next_fund_DK;
                    next_pd_K *= next_pd_DK;
                }
                else
                {
                    fin_not_disc_K *= fin_not_disc_DK;
                    fin_not_log_fwd_K += fin_not_log_fwd_DK;
                    fin_not_fwd_K *= fin_not_fwd_DK;
                }

                for (l = 0; l < num_fund_cpn; l++)
                {
                    eval_const->fund_K[l] *= eval_const->fund_DK[l];
                }

                for (l = 0; l < num_pd_cpn; l++)
                {
                    eval_const->pd_disc_K[l] *= eval_const->pd_disc_DK[l];

                    eval_const->pd_log_fwd_K[l] += eval_const->pd_log_fwd_DK[l];
                    eval_const->pd_fwd_K[l] *= eval_const->pd_fwd_DK[l];
                }
            }

            /* For Optimisation */
            fx_J += fx_DJ_Const + X_J * fx_DJ_Lin;
            fx_K = fx_J;

            X_J += X_DJ;
            X_K = X_J;

            start_J *= start_DJ;
            start_K = start_J;

            in_not_fund_J *= in_not_fund_DJ;
            in_not_fund_K = in_not_fund_J;

            in_not_pd_J *= in_not_pd_DJ;
            in_not_pd_K = in_not_pd_J;

            if (!eval_redemption)
            {
                next_start_fund_J *= next_start_fund_DJ;
                next_start_fund_K = next_start_fund_J;

                next_fund_J *= next_fund_DJ;
                next_fund_K = next_fund_J;

                next_pd_J *= next_pd_DJ;
                next_pd_K = next_pd_J;
            }
            else
            {
                fin_not_disc_J *= fin_not_disc_DJ;
                fin_not_disc_K = fin_not_disc_J;

                fin_not_log_fwd_J += fin_not_log_fwd_DJ;
                fin_not_log_fwd_K = fin_not_log_fwd_J;
                fin_not_fwd_J *= fin_not_fwd_DJ;
                fin_not_fwd_K = fin_not_fwd_J;
            }

            for (l = 0; l < num_fund_cpn; l++)
            {
                eval_const->fund_J[l] *= eval_const->fund_DJ[l];
                eval_const->fund_K[l] = eval_const->fund_J[l];
            }

            for (l = 0; l < num_pd_cpn; l++)
            {
                eval_const->pd_disc_J[l] *= eval_const->pd_disc_DJ[l];
                eval_const->pd_disc_K[l] = eval_const->pd_disc_J[l];

                eval_const->pd_log_fwd_J[l] += eval_const->pd_log_fwd_DJ[l];
                eval_const->pd_log_fwd_K[l] = eval_const->pd_log_fwd_J[l];
                eval_const->pd_fwd_J[l] *= eval_const->pd_fwd_DJ[l];
                eval_const->pd_fwd_K[l] = eval_const->pd_fwd_J[l];
            }
        }

        /* For Optimisation */
        fx_I += fx_DI_Const + X_I * fx_DI_Lin;
        fx_J = fx_I;
        fx_K = fx_I;

        X_I += X_DI;
        X_J = X_I;
        X_K = X_I;

        start_I *= start_DI;
        start_J = start_I;
        start_K = start_I;

        in_not_fund_I *= in_not_fund_DI;
        in_not_fund_J = in_not_fund_I;
        in_not_fund_K = in_not_fund_I;

        in_not_pd_I *= in_not_pd_DI;
        in_not_pd_J = in_not_pd_I;
        in_not_pd_K = in_not_pd_I;

        if (!eval_redemption)
        {
            next_start_fund_I *= next_start_fund_DI;
            next_start_fund_J = next_start_fund_I;
            next_start_fund_K = next_start_fund_I;

            next_fund_I *= next_fund_DI;
            next_fund_J = next_fund_I;
            next_fund_K = next_fund_I;

            next_pd_I *= next_pd_DI;
            next_pd_J = next_pd_I;
            next_pd_K = next_pd_I;
        }
        else
        {
            fin_not_disc_I *= fin_not_disc_DI;
            fin_not_disc_J = fin_not_disc_I;
            fin_not_disc_K = fin_not_disc_I;

            fin_not_log_fwd_I += fin_not_log_fwd_DI;
            fin_not_log_fwd_J = fin_not_log_fwd_I;
            fin_not_log_fwd_K = fin_not_log_fwd_I;
            fin_not_fwd_I *= fin_not_fwd_DI;
            fin_not_fwd_J = fin_not_fwd_I;
            fin_not_fwd_K = fin_not_fwd_I;
        }

        for (l = 0; l < num_fund_cpn; l++)
        {
            eval_const->fund_I[l] *= eval_const->fund_DI[l];
            eval_const->fund_J[l] = eval_const->fund_I[l];
            eval_const->fund_K[l] = eval_const->fund_I[l];
        }

        for (l = 0; l < num_pd_cpn; l++)
        {
            eval_const->pd_disc_I[l] *= eval_const->pd_disc_DI[l];
            eval_const->pd_disc_J[l] = eval_const->pd_disc_I[l];
            eval_const->pd_disc_K[l] = eval_const->pd_disc_I[l];

            eval_const->pd_log_fwd_I[l] += eval_const->pd_log_fwd_DI[l];
            eval_const->pd_log_fwd_J[l] = eval_const->pd_log_fwd_I[l];
            eval_const->pd_log_fwd_K[l] = eval_const->pd_log_fwd_I[l];
            eval_const->pd_fwd_I[l] *= eval_const->pd_fwd_DI[l];
            eval_const->pd_fwd_J[l] = eval_const->pd_fwd_I[l];
            eval_const->pd_fwd_K[l] = eval_const->pd_fwd_I[l];
        }
    }

    /*	End of payoff valuation
            ----------------------- */

    return err;
}

#undef REMNO
#undef SMO
