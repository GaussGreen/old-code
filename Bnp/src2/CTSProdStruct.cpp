
#include "CTSProdStruct.h"

#include "AmortMidatCalib.h"
#include "Fx3FUtils.h"
#include "FxSABRADI.h"
#include "LGMSVCalib.h"
#include "LGMSVClosedForm.h"
#include "LGMSVClosedFormApprox.h"
#include "LGMSVGrfn.h"
#include "LGMSVMC.h"
#include "LGMSVPDE.h"
#include "LGMSVUtil.h"
#include "SABRCalibRRBT.h"
#include "math.h"
#include "opfnctns.h"
#include "opsabrcalib.h"
#include "srt_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define CTS_CALLSPREADMAXSTRIKE 10
#define CTS_CALLSPREADMAXITER 5
#define CTS_CALLSPREADPRECISION 0.00001

#define ONE_MONTH 0.083333333
#define MAX_INST 600
#define MIN_CALIB_TIME 0.03

#define LGMSV_DEFAULT_LAMEPS 0.05
#define LGMSV_DEFAULT_ALPHAEPS 0.30
#define LGMSV_DEFAULT_RHOEPS 0.30
#define LGMSV_ALPHA_MAXITER 20
#define LGMSV_SABRBETA 0.0

/*	Functions for the market structure */
/*	Init */
Err cts_init_mkt(
    long                      today,
    char*                     yc,
    char*                     vc,
    char*                     ref,
    char*                     swap_freq,
    char*                     swap_basis,
    char*                     lib_freq,
    char*                     lib_basis,
    GET_CASH_VOL_FUNC         get_cash_vol,
    GET_CASH_VOL_CONVERT_FUNC get_cash_vol_convert,
    SrtDiffusionType          vol_type,
    int                       cash_vol,
    int                       num_strikes_in_vol,
    double*                   strikes_in_vol,
    CTS_MKT                   mkt)
{
    Err err = NULL;

    mkt->strikes_in_vol = NULL;

    mkt->today = today;
    strcpy(mkt->yc, yc);
    strcpy(mkt->vc, vc);
    strcpy(mkt->ref, ref);
    strcpy(mkt->swap_freq, swap_freq);
    strcpy(mkt->swap_basis, swap_basis);

    err = interp_compounding(swap_freq, &mkt->swap_srt_freq);
    if (err)
        return err;

    err = interp_basis(swap_basis, &mkt->swap_srt_basis);
    if (err)
        return err;

    strcpy(mkt->lib_freq, lib_freq);
    strcpy(mkt->lib_basis, lib_basis);

    err = interp_compounding(lib_freq, &mkt->lib_srt_freq);
    if (err)
        return err;

    err = interp_basis(lib_basis, &mkt->lib_srt_basis);
    if (err)
        return err;

    mkt->get_cash_vol         = get_cash_vol;
    mkt->get_cash_vol_convert = get_cash_vol_convert;
    mkt->vol_type             = vol_type;
    mkt->cash_vol             = cash_vol;

    mkt->num_strikes_in_vol = num_strikes_in_vol;
    mkt->strikes_in_vol     = (double*)calloc(mkt->num_strikes_in_vol, sizeof(double));
    if (!mkt->strikes_in_vol)
        return "Allocation error in cts_init_mkt";
    memcpy(mkt->strikes_in_vol, strikes_in_vol, mkt->num_strikes_in_vol * sizeof(double));

    return NULL;
}

/*	Copy */
Err cts_copy_mkt(CTS_MKT src, CTS_MKT dest)
{
    dest->strikes_in_vol = NULL;

    dest->today = src->today;
    strcpy(dest->yc, src->yc);
    strcpy(dest->vc, src->vc);
    strcpy(dest->ref, src->ref);
    strcpy(dest->swap_freq, src->swap_freq);
    strcpy(dest->swap_basis, src->swap_basis);
    dest->swap_srt_freq  = src->swap_srt_freq;
    dest->swap_srt_basis = src->swap_srt_basis;

    strcpy(dest->lib_freq, src->lib_freq);
    strcpy(dest->lib_basis, src->lib_basis);
    dest->lib_srt_freq  = src->lib_srt_freq;
    dest->lib_srt_basis = src->lib_srt_basis;

    dest->get_cash_vol = src->get_cash_vol;

    dest->num_strikes_in_vol = src->num_strikes_in_vol;
    dest->strikes_in_vol     = (double*)calloc(dest->num_strikes_in_vol, sizeof(double));
    if (!dest->strikes_in_vol)
        return "Allocation error in cts_copy_mkt";
    memcpy(dest->strikes_in_vol, src->strikes_in_vol, dest->num_strikes_in_vol * sizeof(double));

    return NULL;
}

/*	Free */
void cts_free_mkt(CTS_MKT mkt)
{
    if (mkt)
    {
        if (mkt->strikes_in_vol)
            free(mkt->strikes_in_vol);
        mkt->strikes_in_vol = NULL;
    }
}

/*	Functions for the funding leg */

Err cts_fill_fund_leg(
    CTS_MKT mkt,
    /*	Coupons that started before today are disregarded */
    /*	EOD Flag */
    int          eod_flag, /*	0: I, 1: E */
    int          fund_ncpn,
    double*      fund_not,
    long*        fund_fix,
    long*        fund_start,
    long*        fund_end,
    long*        fund_pay,
    char**       fund_basis,
    double*      fund_spr,
    double*      fund_mrg,
    CTS_FUND_LEG fund_leg)
{
    int          i, i0, j;
    SrtBasisCode bas;
    CTS_FUND_CPN cpn;
    long         today;
    Err          err = NULL;

    /*	Get today */
    today = mkt->today;

    /*	Initialise pointers to NULL */
    fund_leg->cpn = NULL;

    /*	Skip coupons fixed before today */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_flag)
    {
        i++;
    }

    /*	Check that at least one coupon is left */
    if (i == fund_ncpn)
    {
        /*	err = "All funding coupons start before today in cts_fill_fund_leg"; */
        fund_leg->num_cpn = 0;
        goto FREE_RETURN;
    }

    /*	Allocate memory */
    fund_leg->num_cpn = fund_ncpn - i;
    fund_leg->cpn     = (cts_fund_cpn*)calloc(fund_leg->num_cpn, sizeof(cts_fund_cpn));
    if (!fund_leg->cpn)
    {
        err = "Allocation error in cts_fill_fund_leg";
        goto FREE_RETURN;
    }

    /*	Fill coupons information */
    j  = 0;
    i0 = i;

    while (i < fund_ncpn)
    {
        cpn = fund_leg->cpn + j;

        /*	Notional */
        cpn->not = fund_not[i];

        /*	Dates */
        cpn->start_date = fund_start[i];
        cpn->end_date   = fund_end[i];
        cpn->pay_date   = fund_pay[i];

        /*	Times */
        cpn->start_time = (cpn->start_date - today) * YEARS_IN_DAY;
        cpn->end_time   = (cpn->end_date - today) * YEARS_IN_DAY;
        cpn->pay_time   = (cpn->pay_date - today) * YEARS_IN_DAY;

        /*	Coupon */
        err = interp_basis(fund_basis[i], &bas);

        if (err)
        {
            goto FREE_RETURN;
        }

        cpn->basis = bas;

        cpn->cvg = coverage(fund_start[i], fund_end[i], bas);
        cpn->cpn = cpn->not *cpn->cvg * (fund_spr[i] + fund_mrg[i]);

        /* calculate the notional refund/increase */
        cpn->cpn_plus_ex         = cpn->cpn - cpn->not ;
        cpn->cpn_plus_ex_partial = cpn->cpn_plus_ex;

        if (i < fund_ncpn - 1)
        {
            cpn->cpn_plus_ex += fund_not[i + 1];
        }

        cpn->mkt_val = cpn->cpn_plus_ex * swp_f_df(today, cpn->pay_date, mkt->yc);

        i++;
        j++;
    }

    err = cts_check_fund_leg(fund_leg);

FREE_RETURN:

    if (err)
    {
        cts_free_fund_leg(fund_leg);
    }

    return err;
}

/*	Check dates consistency */
Err cts_check_fund_leg(CTS_FUND_LEG fund_leg)
{
    int i;

    /*	Check that start and pay dates are increasing */
    for (i = 1; i < fund_leg->num_cpn; i++)
    {
        if (fund_leg->cpn[i].start_date < fund_leg->cpn[i - 1].start_date)
        {
            return "Start dates should be increasing in funding leg";
        }

        if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i - 1].pay_date)
        {
            return "Pay dates should be increasing in funding leg";
        }
    }

    /*	Check that pay dates are after start dates */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i].start_date)
        {
            return "Pay dates should be after start dates in funding leg";
        }
    }

    /*	OK */
    return NULL;
}

/*	Free */
void cts_free_fund_leg(CTS_FUND_LEG fund_leg)
{
    if (fund_leg->cpn)
    {
        free(fund_leg->cpn);
        fund_leg->cpn = NULL;
    }
}

/*	Structures and functions for the exotic leg */

/*	Fixing */

/*	Schedule construction */
static Err Construct_Schedule(
    long           lTodayDate,
    char*          cYcName,
    long           lStartDate,
    long           lTheoEndDate,
    SrtCompounding iFreq,
    SrtBasisCode   iBasis,
    int*           iNbCoupon,
    long*          lCouponDate,
    double*        dCouponTime,
    double*        dCouponCvg,
    double*        dCouponDf,
    double*        dSwpCash,
    double*        dSumDf,
    double*        dLevel)
{
    long   theo_date, act_date, temp_date;
    double level, sum_df;
    int    i, ncpn;
    Err    err = NULL;

    /*	1) First calculate the number of coupon */
    theo_date = lTheoEndDate;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn      = 1;

    while (act_date >= lStartDate - 5)
    {
        theo_date = add_unit(theo_date, -12 / iFreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
        act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
        ncpn++;
    }
    ncpn--;

    if (ncpn < 2)
    {
        err = "Not enough coupons to construct a schedule: check RefRate frequency";
        goto FREE_RETURN;
    }

    /*	3) Fill the schedule informations */

    theo_date = lTheoEndDate;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    i         = ncpn - 1;

    level  = 0.0;
    sum_df = 0.0;

    while (i >= 0)
    {
        dCouponTime[i] = (act_date - lTodayDate) * YEARS_IN_DAY;
        lCouponDate[i] = act_date;

        dCouponDf[i] = swp_f_df(lTodayDate, act_date, cYcName);

        theo_date = add_unit(theo_date, -12 / iFreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

        temp_date     = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
        dCouponCvg[i] = coverage(temp_date, act_date, iBasis);

        act_date = temp_date;

        if (i > 0)
        {
            level += dCouponCvg[i] * dCouponDf[i];
            sum_df += dCouponDf[i];
        }

        i--;

        if (i == 0)
        {
            /* force the first coupon to be on start date */
            act_date = lStartDate;
        }
    }

    *dSwpCash = (dCouponDf[0] - dCouponDf[ncpn - 1]) / level;
    *dSumDf   = sum_df;
    *dLevel   = level;

    *iNbCoupon = ncpn;

FREE_RETURN:

    return err;
}

/*	General initialisation */
Err cts_init_fixing(
    long           today,
    char*          yc,
    int            ref_fix_lag,
    long           ref_start_date,
    long           ref_theo_end_date,
    SrtCompounding srt_freq,
    SrtBasisCode   srt_basis,
    double         ref_fwd_spread,
    CTS_FIX        fixing)
{
    Err err = NULL;

    /*	Fill data structure */

    fixing->ref_fix_date = add_unit(ref_start_date, -ref_fix_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

    if (fixing->ref_fix_date >= today)
    {
        /*	CMS params */
        if (srt_basis == BASIS_ACT_360)
        {
            fixing->ref_conv = 365.0 / 360.0;
        }
        else
        {
            fixing->ref_conv = 1.0;
        }

        fixing->ref_fix_time = (fixing->ref_fix_date - today) * YEARS_IN_DAY;

        err = Construct_Schedule(
            today,
            yc,
            ref_start_date,
            ref_theo_end_date,
            srt_freq,
            srt_basis,
            (int*)&fixing->schedule.iNCpn,
            &(fixing->schedule.lCpnDate[0]),
            &(fixing->schedule.dCpnTime[0]),
            &(fixing->schedule.dCpnCvg[0]),
            &(fixing->schedule.dCpnDf[0]),
            &(fixing->ref_swp_cash),
            &(fixing->ref_sum_df),
            &(fixing->ref_level));

        if (err)
            goto FREE_RETURN;

        fixing->ref_cpnd = 365.0 / (fixing->schedule.lCpnDate[1] - fixing->schedule.lCpnDate[0]);
        fixing->ref_cpnd =
            (fixing->ref_cpnd - floor(fixing->ref_cpnd) < 0.50 ? floor(fixing->ref_cpnd)
                                                               : floor(fixing->ref_cpnd) + 1);

        fixing->schedule.lToday = today;

        fixing->ref_theo_end_date = ref_theo_end_date;
        fixing->ref_fwd_spread    = ref_fwd_spread;
    }
    else
    {
        fixing->schedule.iNCpn = 0;
    }

FREE_RETURN:

    return err;
}

/*	Check */
Err cts_check_fixing(CTS_FIX fixing)
{
    if (fixing->schedule.iNCpn > 0)
    {
        if (fixing->schedule.lCpnDate[0] >= fixing->schedule.lCpnDate[fixing->schedule.iNCpn - 1])
            return "Date inconsistency in cts_check_fixing";
        if (fixing->ref_fix_date > fixing->schedule.lCpnDate[0])
            return "Date inconsistency in cts_check_fixing";
    }

    return NULL;
}

Err cts_init_cpn_payoff(
    double alpha, double beta, int nstr, double* str, double* nbopt, CTS_CPN cpn)
{
    int i;
    Err err = NULL;

    /*	Init pointers */

    cpn->str   = NULL;
    cpn->nbopt = NULL;

    /*	Fill option structure */

    cpn->alpha = alpha;
    cpn->beta  = beta;
    cpn->nstr  = nstr;

    if (nstr)
    {
        cpn->str   = (double*)calloc(cpn->nstr, sizeof(double));
        cpn->nbopt = (double*)calloc(cpn->nstr, sizeof(double));

        if (!cpn->str || !cpn->nbopt)
        {
            err = "Allocation error in cts_init_cpn_payoff";
            goto FREE_RETURN;
        }

        for (i = 0; i < cpn->nstr; i++)
        {
            cpn->str[i]   = str[i];
            cpn->nbopt[i] = nbopt[i];
        }
    }

FREE_RETURN:

    return err;
}

Err cts_init_cpn_payoff_for_range(
    int     buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int     value_zero, /*	0: Do not value low strike options, 1: do */
    double  min_barrier,
    double  max_barrier,
    double  lb,     /*	Lower bound of range */
    double  ub,     /*	Upper bound of range */
    double  payoff, /*	Payoff in range */
    double  call_spread,
    CTS_CPN cpn)
{
    Err     err = NULL;
    double  alpha, beta;
    int     nstr;
    double *str = NULL, *nbopt = NULL;

    /*	Quick check on range */
    if (lb >= ub - 2.0 * call_spread)
        return "Range inconsistency in cts_init_fixing_for_range";
    if (ub <= call_spread && !value_zero)
        return "Range inconsistency in cts_init_fixing_for_range";
    if (call_spread < 1.0e-08)
        return "Zero call spread in cts_init_fixing_for_range";

    /*	Negative strike means no strike */
    if (lb < min_barrier)
        value_zero = 0;

    /*	Calculation of strikes and number of options */
    if (value_zero || lb >= call_spread)
    {
        if (ub < max_barrier)
        {
            /*	Case 1: lower and upper strikes, 4 options */

            alpha = 0.0;
            beta  = 0.0;
            nstr  = 4;

            str   = (double*)calloc(nstr, sizeof(double));
            nbopt = (double*)calloc(nstr, sizeof(double));
            if (!str || !nbopt)
            {
                err = "Allocation error in cts_init_cpn_payoff_for_range";
                goto FREE_RETURN;
            }

            if (buy_sell == 1)
            /*	We buy: subreplicate */
            {
                str[0] = lb;
                str[1] = lb + call_spread;
                str[2] = ub - call_spread;
                str[3] = ub;
            }
            else
            /*	We sell: overreplicate */
            {
                str[0] = lb - call_spread;
                str[1] = lb;
                str[2] = ub;
                str[3] = ub + call_spread;
            }
            /*	Calculate option structure */
            nbopt[0] = payoff / call_spread;
            nbopt[1] = -nbopt[0];
            nbopt[2] = nbopt[1];
            nbopt[3] = nbopt[0];
        }
        else
        {
            /* Case 2: no higer strike: 2 options */
            alpha = 0.0;
            beta  = 0;
            nstr  = 2;

            str   = (double*)calloc(nstr, sizeof(double));
            nbopt = (double*)calloc(nstr, sizeof(double));
            if (!str || !nbopt)
            {
                err = "Allocation error in cts_init_cpn_payoff_for_range";
                goto FREE_RETURN;
            }

            if (buy_sell == 1)
            /*	We buy: subreplicate */
            {
                str[0] = lb;
                str[1] = lb + call_spread;
            }
            else
            /*	We sell: overreplicate */
            {
                str[0] = lb - call_spread;
                str[1] = lb;
            }
            /*	Calculate option structure */
            nbopt[0] = payoff / call_spread;
            nbopt[1] = -nbopt[0];
        }
    }
    else if (ub < max_barrier)
    /*	Case 3: no lower strike, 2 options */
    {
        alpha = payoff;
        beta  = 0;
        nstr  = 2;

        str   = (double*)calloc(nstr, sizeof(double));
        nbopt = (double*)calloc(nstr, sizeof(double));
        if (!str || !nbopt)
        {
            err = "Allocation error in cts_init_cpn_payoff_for_range";
            goto FREE_RETURN;
        }

        if (buy_sell == 1)
        /*	We buy: subreplicate */
        {
            str[0] = ub - call_spread;
            str[1] = ub;
        }
        else
        /*	We sell: overreplicate */
        {
            str[0] = ub;
            str[1] = ub + call_spread;
        }
        /*	Calculate option structure */
        nbopt[0] = -payoff / call_spread;
        nbopt[1] = -nbopt[0];
    }
    else
    {
        /* Case 4: no higher, no lower strikes */
        alpha = payoff;
        beta  = 0;
        nstr  = 0;
    }

    /*	Init fixing with option structure */
    err = cts_init_cpn_payoff(alpha, beta, nstr, str, nbopt, cpn);

FREE_RETURN:

    if (str)
        free(str);
    if (nbopt)
        free(nbopt);

    return err;
}

/* payoff = (gearing * PayLibor + margin) floored and capped */
Err cts_init_cpn_payoff_for_cif(
    int     value_zero,
    double  min_barrier,
    double  max_barrier,
    double  gearing,
    double  margin,
    double  floor,
    double  cap,
    CTS_CPN cpn)
{
    double  strike1, strike2;
    Err     err = NULL;
    double  alpha, beta;
    int     nstr;
    double *str = NULL, *nbopt = NULL;

    /*	Quick check on range */
    if (floor > cap)
        return "Cap / Floor inconsistency in cts_init_cpn_payoff_for_cif";

    if (fabs(gearing) < 1.0E-08)
    {
        /* simple midat */
        alpha = margin;
        beta  = 0;
        nstr  = 0;
    }
    else
    {
        if (gearing > 0)
        {
            strike1 = (floor - margin) / gearing;
            strike2 = (cap - margin) / gearing;
        }
        else
        {
            strike2 = (floor - margin) / gearing;
            strike1 = (cap - margin) / gearing;
        }

        if (strike1 > min_barrier)
        {
            if (strike2 < max_barrier)
            {
                /* Both floor and cap are valid */
                nstr = 2;

                str   = (double*)calloc(nstr, sizeof(double));
                nbopt = (double*)calloc(nstr, sizeof(double));

                if (!str || !nbopt)
                {
                    err = "Allocation error in cts_init_cpn_payoff_for_cif";
                    goto FREE_RETURN;
                }

                alpha = margin + gearing * strike1;
                beta  = 0.0;

                nbopt[0] = gearing;
                str[0]   = strike1;
                nbopt[1] = -gearing;
                str[1]   = strike2;
            }
            else
            {
                nstr = 1;

                str   = (double*)calloc(nstr, sizeof(double));
                nbopt = (double*)calloc(nstr, sizeof(double));

                if (!str || !nbopt)
                {
                    err = "Allocation error in cts_init_cpn_payoff_for_cif";
                    goto FREE_RETURN;
                }

                alpha = margin + gearing * strike1;
                beta  = 0.0;

                nbopt[0] = gearing;
                str[0]   = strike1;
            }
        }
        else
        {
            if (strike2 > min_barrier && strike2 < max_barrier)
            {
                nstr = 1;

                str   = (double*)calloc(nstr, sizeof(double));
                nbopt = (double*)calloc(nstr, sizeof(double));

                if (!str || !nbopt)
                {
                    err = "Allocation error in cts_init_cpn_payoff_for_cif";
                    goto FREE_RETURN;
                }

                alpha = margin;
                beta  = gearing;

                nbopt[0] = -gearing;
                str[0]   = strike2;
            }
            else if (strike2 > max_barrier)
            {
                nstr  = 0;
                alpha = margin;
                beta  = gearing;
            }
            else
            {
                nstr  = 0;
                alpha = margin + gearing * strike2;
                beta  = 0.0;
            }
        }
    }

    /*	Init fixing with option structure */
    err = cts_init_cpn_payoff(alpha, beta, nstr, str, nbopt, cpn);

FREE_RETURN:

    if (str)
        free(str);
    if (nbopt)
        free(nbopt);

    return err;
}

/* Coupon simple valuation */
Err cts_calc_shifted_range_pv(
    int     buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int     value_zero, /*	0: Do not value low strike options, 1: do */
    double  min_barrier,
    double  max_barrier,
    double  lb_init, /*	Lower bound of range */
    double  ub_init, /*	Upper bound of range */
    double  payoff,  /*	Payoff in range */
    double  shift_strike,
    double  call_spread,
    double  numer_spread,
    CTS_MKT mkt,
    CTS_CPN cpn,
    double* alpha,
    double* beta,
    int*    nstr,
    double* str,
    double* nbopt,
    double* price)
{
    Err err = NULL;

    double ref_swp_cash, ref_swp_fwd, vol, opt;

    double lb, ub;
    double cpn_value;
    double power;

    int i, j;

    CTS_FIX fix;

    /*	Shift the barriers */
    lb = lb_init - shift_strike;
    ub = ub_init + shift_strike;

    /*	Negative strike means no strike */
    if (lb_init < -min_barrier)
        value_zero = 0;

    /*	Calculation of strikes and number of options */
    if (value_zero || lb_init >= call_spread)
    {
        if (ub_init < max_barrier)
        {
            /*	Case 1: lower and upper strikes, 4 options */

            *alpha = 0.0;
            *beta  = 0.0;
            *nstr  = 4;

            if (buy_sell == 1)
            /*	We buy: subreplicate */
            {
                str[0] = lb;
                str[1] = lb + numer_spread;
                str[2] = ub - numer_spread;
                str[3] = ub;
            }
            else
            /*	We sell: overreplicate */
            {
                str[0] = lb - numer_spread;
                str[1] = lb;
                str[2] = ub;
                str[3] = ub + numer_spread;
            }
            /*	Calculate option structure */
            nbopt[0] = payoff / numer_spread;
            nbopt[1] = -nbopt[0];
            nbopt[2] = nbopt[1];
            nbopt[3] = nbopt[0];
        }
        else
        {
            /* Case 2: no higer strike: 2 options */
            *alpha = 0.0;
            *beta  = 0;
            *nstr  = 2;

            if (buy_sell == 1)
            /*	We buy: subreplicate */
            {
                str[0] = lb;
                str[1] = lb + numer_spread;
            }
            else
            /*	We sell: overreplicate */
            {
                str[0] = lb - numer_spread;
                str[1] = lb;
            }
            /*	Calculate option structure */
            nbopt[0] = payoff / numer_spread;
            nbopt[1] = -nbopt[0];
        }
    }
    else if (ub_init < max_barrier)
    /*	Case 3: no lower strike, 2 options */
    {
        *alpha = payoff;
        *beta  = 0;
        *nstr  = 2;

        if (buy_sell == 1)
        /*	We buy: subreplicate */
        {
            str[0] = ub - numer_spread;
            str[1] = ub;
        }
        else
        /*	We sell: overreplicate */
        {
            str[0] = ub;
            str[1] = ub + numer_spread;
        }
        /*	Calculate option structure */
        nbopt[0] = -payoff / numer_spread;
        nbopt[1] = -nbopt[0];
    }
    else
    {
        /* Case 4: no higher, no lower strikes */
        *alpha = payoff;
        *beta  = 0;
        *nstr  = 0;
    }

    /* Check the strikes */
    if (*nstr > 0 && str[0] < 0.0)
    {
        err = "Negative implied strike, decrease the numerical call spread";
        goto FREE_RETURN;
    }

    /*	Do the valuation */
    cpn_value = 0.0;

    for (i = 0; i < cpn->used_nfix; i++)
    {
        fix = &(cpn->fix[i]);

        /* Calculate the forward */
        ref_swp_cash = fix->ref_swp_cash;
        ref_swp_fwd  = ref_swp_cash + fix->ref_fwd_spread;

        cpn_value += *alpha;

        for (j = 0; j < *nstr; j++)
        {
            err = mkt->get_cash_vol(
                mkt->vc,
                fix->schedule.lCpnDate[0],
                fix->ref_theo_end_date,
                str[j],
                0,
                mkt->ref,
                &vol,
                &power);

            if (err)
                goto FREE_RETURN;

            if (power > 0.5)
            {
                opt = srt_f_optblksch(
                    ref_swp_fwd,
                    str[j],
                    vol,
                    fix->ref_fix_time,
                    1.0,
                    SRT_CALL,
                    (SrtGreekType)SRT_PREMIUM);
            }
            else
            {
                opt = srt_f_optblknrm(
                    ref_swp_fwd,
                    str[j],
                    vol,
                    fix->ref_fix_time,
                    1.0,
                    SRT_CALL,
                    (SrtGreekType)SRT_PREMIUM);
            }

            cpn_value += nbopt[j] * opt;
        }
    }

    *price = cpn_value;

FREE_RETURN:

    return err;
}

/*	Adjust call spread values */
Err cts_adjust_call_spread(
    int     buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int     value_zero, /*	0: Do not value low strike options, 1: do */
    double  min_barrier,
    double  max_barrier,
    double  lb,     /*	Lower bound of range */
    double  ub,     /*	Upper bound of range */
    double  payoff, /*	Payoff in range */
    CTS_CPN cpn,
    CTS_MKT mkt,
    double  call_spread,
    double  numer_spread,
    double  precision,
    int     maxiter)
{
    Err err = NULL;

    double   price_tgt, price1, price2;
    double** res_iter;
    double   shift1, shift2;
    int      i, j, l;

    double alpha, beta;
    int    nstr;
    double nbopt[CTS_CALLSPREADMAXSTRIKE], str[CTS_CALLSPREADMAXSTRIKE];

    res_iter = dmatrix(0, maxiter, 0, 1);

    if (!res_iter)
    {
        err = "Memory allocation faillure in cts_adjust_call_spread";
        goto FREE_RETURN;
    }

    err = cts_calc_shifted_range_pv(
        buy_sell,
        value_zero,
        min_barrier,
        max_barrier,
        lb,
        ub,
        payoff,
        0,
        call_spread,
        call_spread,
        mkt,
        cpn,
        &alpha,
        &beta,
        &nstr,
        str,
        nbopt,
        &price_tgt);

    if (err)
        goto FREE_RETURN;

    /* Up shift */
    shift1 = -numer_spread / 2.0;

    err = cts_calc_shifted_range_pv(
        buy_sell,
        value_zero,
        min_barrier,
        max_barrier,
        lb,
        ub,
        payoff,
        shift1,
        call_spread,
        numer_spread,
        mkt,
        cpn,
        &alpha,
        &beta,
        &nstr,
        str,
        nbopt,
        &price1);

    if (err)
        goto FREE_RETURN;

    res_iter[0][0] = shift1;
    res_iter[0][1] = price1;

    i = 1;

    /* Second shift, down */
    if (lb > 0.0)
    {
        /* make sure to not get a negative strike */
        if (buy_sell == 1)
        {
            shift2 = 0.5 * lb;
        }
        else
        {
            shift2 = 0.5 * (lb - numer_spread);
        }

        shift2 = min(shift2, -shift1);
    }
    else
    {
        shift2 = -shift1;
    }

    while (fabs(price_tgt - price1) > precision && i < maxiter)
    {
        i++;

        err = cts_calc_shifted_range_pv(
            buy_sell,
            value_zero,
            min_barrier,
            max_barrier,
            lb,
            ub,
            payoff,
            shift2,
            call_spread,
            numer_spread,
            mkt,
            cpn,
            &alpha,
            &beta,
            &nstr,
            str,
            nbopt,
            &price2);

        if (err)
            goto FREE_RETURN;

        /* Save Res */
        l = 0;
        while (l < i - 1 && res_iter[l][1] < price2)
        {
            l++;
        }

        if (l < i - 1)
        {
            for (j = i - 2; j >= l; j--)
            {
                res_iter[j + 1][0] = res_iter[j][0];
                res_iter[j + 1][1] = res_iter[j][1];
            }
        }

        res_iter[l][0] = shift2;
        res_iter[l][1] = price2;

        shift1 = shift2;
        price1 = price2;

        shift2 = solve_for_next_coef_dlm(res_iter, i, price_tgt, 2);
    }

    /* Now adjust the cpn values */

    err = cts_calc_shifted_range_pv(
        buy_sell,
        value_zero,
        min_barrier,
        max_barrier,
        lb,
        ub,
        payoff,
        shift2,
        call_spread,
        numer_spread,
        mkt,
        cpn,
        &alpha,
        &beta,
        &nstr,
        str,
        nbopt,
        &price2);

    if (err)
        goto FREE_RETURN;

    cpn->alpha = alpha;
    cpn->beta  = beta;
    cpn->nstr  = nstr;

    for (i = 0; i < nstr; i++)
    {
        cpn->str[i]   = str[i];
        cpn->nbopt[i] = nbopt[i];
    }

FREE_RETURN:

    if (res_iter)
        free_dmatrix(res_iter, 0, maxiter, 0, 1);

    return err;
}

/*	Coupon */
/*	Init */
Err cts_init_cpn(
    /*	Market */
    CTS_MKT mkt,
    /*	Coupon details */
    long   cpn_start_date,
    long   cpn_end_date,
    long   cpn_pay_date,
    double cpn_coupon,
    char*  cpn_basis,
    double cpn_not,
    /*	General libor fixing and basis properties */
    int fix_lag_bd,
    /*	Coupon type */
    int cpn_type, /*	0:	fixed
                                  1:	libor fixed at start
                                  2:	libor fixed at end */
    /*	Paid libor details */
    long   pay_fix_date, /*	fixing date of the Libor */
    int    pay_months,   /*	Length of the paid libor in months */
    char*  pay_freq,
    char*  pay_basis,
    double pay_fwd_spread,
    double pay_gearing,
    /*	Fixing details */
    int     nfix,
    double* weights,
    /*		ref libor */
    long*   ref_fix_dates, /*	Fixing dates */
    char*   fix_tenor,     /*	Fixing tenors */
    char*   fix_freq,
    char*   fix_basis,
    double* ref_fwd_spreads, /*	Reference libor spreads corresponding to fixing dates */
    /*	Profile details */
    int tstype, /*	0: generic, 1: range */
    /*		Profile 0 */
    double  alpha,
    double  beta,
    int     nstr,
    double* str,
    double* nbopt,
    /*		Profile 1 */
    int    buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int    value_zero, /*	0: Do not value low strike options, 1: do */
    double lb,         /*	Lower bound of range */
    double ub,         /*	Upper bound of range */
    double payoff,     /*	Payoff in range */
    double call_spread,
    double numer_spread,
    /*		Trimming */
    int trim_type, /*	0: no trim
                                   1: x fixings max
                                   2: x time min between two fixings */
    int    max_fix,
    double min_fix_time,
    /*		Pricing params */
    int    calc_mkt_iv,
    int    use_cmsopt,   /*	Use CmsOption to value fix option, use BS on CmsRate otherwise */
    double correl_start, /*	Correl between libor fixing and libor paid at start */
    double correl_end,   /*	Correl between libor fixing and libor paid at start */

    int float_adjust_type, /*	type of adjustment for the floating coupon,
                                   0: ATM vol, 1: Strike Vol */
    CTS_CPN cpn)
{
    SrtBasisCode   srt_basis, srt_used_basis;
    SrtCompounding srt_used_freq;
    long           fix, start, end;
    int            i;
    long           today;
    Err            err = NULL;

    /*	Get today */
    today = mkt->today;

    /*	Init arrays */
    cpn->weights         = NULL;
    cpn->fix             = NULL;
    cpn->used_fixidx     = NULL;
    cpn->used_fixweights = NULL;

    /*	Libor for floating coupon */
    cpn->type = cpn_type;
    if (cpn->type)
    /*	Floating coupon */
    {
        /*	Paid libor start and end dates */
        cpn->pay_fix_date = pay_fix_date;
        cpn->pay_start_date =
            add_unit(cpn->pay_fix_date, fix_lag_bd, SRT_BDAY, MODIFIED_SUCCEEDING);

        if (cpn->type == 3 || cpn->type == 4)
        {
            err = add_tenor(
                cpn->pay_start_date, fix_tenor, NO_BUSDAY_CONVENTION, &cpn->pay_theo_end_date);
        }
        else
        {
            cpn->pay_theo_end_date =
                add_unit(cpn->pay_start_date, pay_months, SRT_MONTH, NO_BUSDAY_CONVENTION);
        }

        cpn->pay_end_date = bus_date_method(cpn->pay_theo_end_date, MODIFIED_SUCCEEDING);

        /*	Times */
        cpn->pay_fix_time   = (cpn->pay_fix_date - today) * YEARS_IN_DAY;
        cpn->pay_start_time = (cpn->pay_start_date - today) * YEARS_IN_DAY;
        cpn->pay_end_time   = (cpn->pay_end_date - today) * YEARS_IN_DAY;

        /*	Coverage */
        /*	Read basis */
        err = interp_basis(pay_basis, &srt_basis);
        if (err)
            goto FREE_RETURN;

        cpn->pay_cvg = coverage(cpn->pay_start_date, cpn->pay_end_date, srt_basis);

        /*	DRS informations */
        if (srt_basis == BASIS_ACT_360)
        {
            cpn->pay_conv = 365.0 / 360.0;
        }
        else
        {
            cpn->pay_conv = 1.0;
        }

        cpn->pay_cpnd = 12 / pay_months;

        /*	Forward spread */
        cpn->pay_fwd_spread = pay_fwd_spread;
    }
    cpn->pay_gearing = pay_gearing;

    /*	Coupon */

    /*	Dates */
    cpn->cpn_start_date = cpn_start_date;
    cpn->cpn_end_date   = cpn_end_date;
    cpn->cpn_pay_date   = cpn_pay_date;

    /*	Times */
    cpn->cpn_start_time = (cpn->cpn_start_date - today) * YEARS_IN_DAY;
    cpn->cpn_end_time   = (cpn->cpn_end_date - today) * YEARS_IN_DAY;
    cpn->cpn_pay_time   = (cpn->cpn_pay_date - today) * YEARS_IN_DAY;

    /*	Paid Coupon  */
    cpn->cpn_coupon = cpn_coupon;

    /*	Coverage */
    /*	Read basis */
    err = interp_basis(cpn_basis, &srt_basis);
    if (err)
        goto FREE_RETURN;

    cpn->cpn_basis = srt_basis;

    cpn->cpn_cvg = coverage(cpn->cpn_start_date, cpn->cpn_end_date, srt_basis);

    /*	Notional */
    cpn->cpn_not = cpn_not;

    /*	Fixings */
    cpn->nfix    = nfix;
    cpn->weights = (double*)calloc(cpn->nfix, sizeof(double));

    if (!cpn->weights)
    {
        err = "Allocation error in cts_init_cpn";
        goto FREE_RETURN;
    }

    memcpy(cpn->weights, weights, cpn->nfix * sizeof(double));

    /*	Select the fixings we will use */
    err = cts_trim_fixings(cpn, trim_type, max_fix, min_fix_time);
    if (err)
        goto FREE_RETURN;

    cpn->fix = (cts_fix*)calloc(cpn->used_nfix, sizeof(cts_fix));

    if (!cpn->fix)
    {
        err = "Allocation error (2) in cts_init_cpn";
        goto FREE_RETURN;
    }

    if (cpn_type < 3)
    {
        /*	Interpretation of basis */
        err = interp_basis(fix_basis, &srt_used_basis);
        if (err)
            goto FREE_RETURN;

        err = interp_compounding(fix_freq, &srt_used_freq);
        if (err)
            goto FREE_RETURN;
    }
    else
    {
        /*	Interpretation of basis */
        err = interp_basis(pay_basis, &srt_used_basis);
        if (err)
            goto FREE_RETURN;

        err = interp_compounding(pay_freq, &srt_used_freq);
        if (err)
            goto FREE_RETURN;
    }

    /*	Initialise the used fixings */
    for (i = 0; i < cpn->used_nfix; i++)
    {
        fix   = ref_fix_dates[cpn->used_fixidx[i]];
        start = add_unit(fix, fix_lag_bd, SRT_BDAY, MODIFIED_SUCCEEDING);
        err   = add_tenor(start, fix_tenor, NO_BUSDAY_CONVENTION, &end);

        err = cts_init_fixing(
            today,
            mkt->yc,
            fix_lag_bd,
            start,
            end,
            srt_used_freq,
            srt_used_basis,
            ref_fwd_spreads[cpn->used_fixidx[i]],
            &(cpn->fix[i]));

        if (err)
        /*	Error: free the fixings already initialised */
        {
            free(cpn->fix);
            cpn->fix = NULL;

            goto FREE_RETURN;
        }
    }

    /*  Options */

    if (tstype == 0)
    /*	General */
    {
        err = cts_init_cpn_payoff(alpha, beta, nstr, str, nbopt, cpn);
    }
    else
    /*	Range */
    {
        if (cpn->type < 3)
        {
            /* classical Time Swap */
            err = cts_init_cpn_payoff_for_range(
                buy_sell,
                value_zero,
                CTS_MINBARRIER,
                CTS_MAXBARRIER,
                lb,
                ub,
                payoff,
                call_spread,
                cpn);

            numer_spread = max(numer_spread, call_spread);

            if (numer_spread - call_spread > 1.0E-08)
            {
                /* adjust the call_spread */
                err = cts_adjust_call_spread(
                    buy_sell,
                    value_zero,
                    CTS_MINBARRIER,
                    CTS_MAXBARRIER,
                    lb,
                    ub,
                    payoff,
                    cpn,
                    mkt,
                    call_spread,
                    numer_spread,
                    CTS_CALLSPREADPRECISION,
                    CTS_CALLSPREADMAXITER);

                if (err)
                    goto FREE_RETURN;
            }
        }
        else
        {
            /* CIF type */
            err = cts_init_cpn_payoff_for_cif(
                value_zero, 1.0E-08, CTS_MAXBARRIER, pay_gearing, cpn_coupon, lb, ub, cpn);
        }
    }

    if (err)
        goto FREE_RETURN;

    /*	Store caplet strike as ub only in case of a range */
    if (tstype)
    {
        if (cpn->type < 3)
        {
            /* Time Swap case: strike on barrier */

            if (lb < -CTS_MINBARRIER)
                value_zero = 0;

            /*	Calculation of strikes and number of options */
            if (value_zero || lb >= numer_spread)
            {
                if (ub < CTS_MAXBARRIER)
                {
                    cpn->caplet_strike = CTS_NA;
                }
                else
                {
                    if (lb > numer_spread)
                    {
                        cpn->caplet_strike = lb;
                    }
                    else
                    {
                        cpn->caplet_strike = CTS_NA;
                    }
                }
            }
            else if (ub < CTS_MAXBARRIER)
            /*	Case 3: no lower strike, 2 options */
            {
                cpn->caplet_strike = ub;
            }
            else
            {
                cpn->caplet_strike = CTS_NA;
            }
        }
        else
        {
            /* CIF case strike on IF */
            if (cpn->nstr > 0)
            {
                /* we select the first one ! */
                cpn->caplet_strike = cpn->str[0];
            }
            else
            {
                cpn->caplet_strike = CTS_NA;
            }
        }
    }
    else
    {
        cpn->caplet_strike = CTS_NA;
    }

    /*	Calc PV */
    if (calc_mkt_iv)
    {
        err = cts_calc_cpn_mkt_value(
            cpn, mkt, use_cmsopt, correl_start, correl_end, float_adjust_type);
        if (err)
            goto FREE_RETURN;
    }

    /*	Check */
    err = cts_check_coupon(cpn);
    if (err)
        goto FREE_RETURN;

FREE_RETURN:

    if (err)
    {
        cts_free_coupon(cpn);
    }

    return err;
}

/*	Trim fixings */
Err cts_trim_fixings(
    CTS_CPN cpn,
    int     trim_type, /*	0: no trim
                                       1: x fixings max
                                       2: x time min between two fixings */
    int    max_fix,
    double min_fix_time)
{
    int    i, step, nb_points;
    double new_weight;
    Err    err = NULL;

    if (trim_type == 2 && fabs(min_fix_time) < YEARS_IN_DAY)
    {
        trim_type = 0;
    }

    if (trim_type == 0)
    {
        /* just recopy all of them */
        cpn->used_nfix = cpn->nfix;

        cpn->used_fixidx     = (int*)calloc(cpn->used_nfix, sizeof(int));
        cpn->used_fixweights = (double*)calloc(cpn->used_nfix, sizeof(double));

        if (!cpn->used_fixidx || !cpn->used_fixweights)
        {
            err = "Memory allocation faillure in cts_trim_fixings";
            goto FREE_RETURN;
        }

        for (i = 0; i < cpn->used_nfix; i++)
        {
            cpn->used_fixidx[i]     = i;
            cpn->used_fixweights[i] = cpn->weights[i];
        }
    }
    else if (trim_type == 2)
    {
        max_fix   = (int)((cpn->cpn_pay_time - cpn->cpn_start_time) / min_fix_time + 0.5) - 1;
        max_fix   = max(max_fix, 1);
        trim_type = 1;
    }

    if (trim_type == 1)
    {
        /* we choose max_fix fixings */

        /* choose the actual number of points */
        max_fix = min(cpn->nfix, max_fix);

        nb_points = 10000;
        max_fix   = max_fix + 1;

        while (nb_points > cpn->nfix && max_fix > 1)
        {
            max_fix--;
            step      = (int)((cpn->nfix - 1) / (max_fix + 1));
            nb_points = (max_fix - 1) * step + max_fix;
        }

        cpn->used_nfix = max_fix;

        cpn->used_fixidx     = (int*)calloc(cpn->used_nfix, sizeof(int));
        cpn->used_fixweights = (double*)calloc(cpn->used_nfix, sizeof(double));

        if (!cpn->used_fixidx || !cpn->used_fixweights)
        {
            err = "Memory allocation faillure in cts_trim_fixings";
            goto FREE_RETURN;
        }

        /*  choose the used fixings */

        cpn->used_fixidx[0] = (int)((cpn->nfix - nb_points) / 2.0);

        new_weight = 0.0;

        for (i = 0; i < cpn->nfix; i++)
        {
            new_weight += cpn->weights[i];
        }

        new_weight = new_weight / ((double)(cpn->used_nfix));

        cpn->used_fixweights[0] = new_weight;

        for (i = 1; i < cpn->used_nfix; i++)
        {
            cpn->used_fixidx[i]     = cpn->used_fixidx[i - 1] + step + 1;
            cpn->used_fixweights[i] = new_weight;
        }
    }

FREE_RETURN:

    return err;
}

Err cts_calc_fixed_cpn_mkt_value(
    /*	Market */
    CTS_MKT mkt,
    /*	Coupon details */
    long   cpn_start_date,
    long   cpn_end_date,
    long   cpn_pay_date,
    double cpn_coupon,
    char*  cpn_basis,
    double cpn_not,
    /*	General libor fixing and basis properties */
    int fix_lag_bd,
    /*	Coupon type */
    int cpn_type, /*	0:	fixed
                                  1:	libor fixed at start
                                  2:	libor fixed at end
                                  3:	midat */
    /*	Paid libor details */
    long   pay_fix_date, /*	fixing date of the Libor */
    int    pay_months,   /*	Length of the paid libor in months */
    char*  pay_freq,
    char*  pay_basis,
    double pay_fwd_spread,
    double pay_gearing,
    double fixed_pay,
    /*	Fixing details */
    int     nfix,
    double* weights,
    /*		ref libor */
    double* fixed_ref, /*	Fixing dates */
    /*	Profile details */
    int tstype, /*	0: generic, 1: range */
    /*		Profile 0 */
    double  alpha,
    double  beta,
    int     nstr,
    double* str,
    double* nbopt,
    /*		Profile 1 */
    double  lb,     /*	Lower bound of range */
    double  ub,     /*	Upper bound of range */
    double  payoff, /*	Payoff in range */
    double  accrue_on_barrier,
    int     eod_fix_flag, /*	0: I, 1: E */
    int     eod_pay_flag, /*	0: I, 1: E */
    int     eod_ex_flag,  /*	0: I, 1: E */
    double* mkt_value)
{
    Err          err = NULL;
    double       df_pay, cpn_cvg, cpn_pv, pay_cvg, libor;
    int          j, l;
    double       temp, adj_bar;
    SrtBasisCode bas;
    long         pay_start_date, pay_end_date;

    if (accrue_on_barrier)
    {
        adj_bar = 1.0E-08;
    }
    else
    {
        adj_bar = 0.0;
    }

    df_pay  = swp_f_df(mkt->today, cpn_pay_date, mkt->yc);
    err     = interp_basis(cpn_basis, &bas);
    cpn_cvg = coverage(cpn_start_date, cpn_end_date, bas);

    /* first part has fixed fixing */

    cpn_pv = 0.0;

    if (tstype == 1)
    {
        if (cpn_type < 3)
        {
            /* classical time swap range */
            for (j = 0; j < nfix; j++)
            {
                if (fixed_ref[j] < ub + adj_bar && fixed_ref[j] > lb - adj_bar)
                {
                    cpn_pv += weights[j];
                }
            }

            switch (cpn_type)
            {
            case 0:
                cpn_pv *= cpn_coupon * cpn_cvg * cpn_not * df_pay;
                break;
            case 1:
            case 2:

                /*
                err = add_tenor (cpn->pay_start_date, fix_tenor, NO_BUSDAY_CONVENTION,
                &cpn->pay_end_date);
                */

                if (pay_fix_date < mkt->today + eod_fix_flag)
                {
                    libor = fixed_pay;
                }
                else
                {
                    pay_start_date =
                        add_unit(pay_fix_date, fix_lag_bd, SRT_BDAY, MODIFIED_SUCCEEDING);
                    pay_end_date =
                        add_unit(pay_start_date, pay_months, SRT_MONTH, MODIFIED_SUCCEEDING);

                    err     = interp_basis(pay_basis, &bas);
                    pay_cvg = coverage(pay_start_date, pay_end_date, bas);

                    libor = (swp_f_df(mkt->today, pay_start_date, mkt->yc) /
                                 swp_f_df(mkt->today, pay_end_date, mkt->yc) -
                             1.0) /
                                pay_cvg +
                            pay_fwd_spread;
                }

                cpn_pv *= (pay_gearing * libor + cpn_coupon) * cpn_cvg * cpn_not * df_pay;

                break;
            }
        }
        else
        {
            cpn_pv = pay_gearing * fixed_pay + cpn_coupon;

            if (cpn_pv < lb)
            {
                /* floor */
                cpn_pv = lb;
            }
            if (cpn_pv > ub)
            {
                /* cap */
                cpn_pv = ub;
            }

            cpn_pv *= cpn_cvg * cpn_not * df_pay;
        }
    }
    else
    {
        for (j = 0; j < nfix; j++)
        {
            temp = alpha + beta * fixed_ref[j];

            for (l = 0; l < nstr; l++)
            {
                if (fixed_ref[j] > str[l])
                {
                    temp += fixed_ref[j] - str[l];
                }
            }

            cpn_pv += temp * weights[j];
        }

        cpn_pv *= cpn_cvg * cpn_not * df_pay;
    }

    *mkt_value = cpn_pv;

    return err;
}

/*	Calc market values */
Err cts_calc_cpn_mkt_value(
    CTS_CPN cpn,
    CTS_MKT mkt,
    int     use_cmsopt,
    double  correl_start,  /*	Correl between libor fixing and libor paid at start */
    double  correl_end,    /*	Correl between libor fixing and libor paid at start */
    int     float_adjust_type) /*	type of adjustment for the floating coupon,
0: ATM vol, 1: Strike Vol */
{
    Err     err = NULL;
    CTS_FIX fix;
    double  fixing_value, fixing_value_mrg, cpn_value, cpn_value_mrg;
    double  df_start, df_end, df_pay, opt, delta;
    double  ref_swp_cash, ref_swp_fwd, ref_swp_drs, ref_delay, pay_adjlogvol, ref_adjlogvol,
        ref_atmlogvol, ref_atmlogstd, ref_strikelogvol, ref_adj;
    double ref_strikestd, pay_new_strike;
    double pay_swp_cash, pay_swp_fwd, pay_swp_drs, pay_delay, pay_atmlogvol, pay_atmlogstd;
    double pay_regpaydelay, ref_correl;
    int    i, j;
    int    do_margin;
    double temp;

    cpn_value     = 0.0;
    cpn_value_mrg = 0.0;

    if ((cpn->type == 1 || cpn->type == 2) && fabs(cpn->cpn_coupon) < 1.0E-10)
    {
        do_margin = 0;
    }
    else
    {
        do_margin = 1;
    }

    if (cpn->type == 1 || cpn->type == 2)
    {
        /* Calculation on the payment Libor */
        df_start     = swp_f_df(mkt->today, cpn->pay_start_date, mkt->yc);
        df_end       = swp_f_df(mkt->today, cpn->pay_end_date, mkt->yc);
        pay_swp_cash = (df_start - df_end) / (df_end * cpn->pay_cvg);
        pay_swp_fwd  = pay_swp_cash + cpn->pay_fwd_spread;

        pay_delay       = cpn->cpn_pay_time - cpn->pay_start_time;
        pay_regpaydelay = cpn->cpn_pay_time - cpn->cpn_start_time;

        err = swp_f_Cms_Rate(
            pay_swp_fwd,
            cpn->pay_fix_time,
            1,
            cpn->pay_cpnd,
            pay_delay,
            cpn->pay_conv,
            mkt->vol_type,
            0.0,
            1,
            cpn->pay_start_date,
            cpn->pay_theo_end_date,
            (SRT_Boolean)mkt->cash_vol,
            cpn->pay_fwd_spread,
            mkt->vc,
            mkt->num_strikes_in_vol,
            mkt->strikes_in_vol,
            &pay_swp_drs);

        if (err)
            goto FREE_RETURN;

        if (cpn->pay_fix_time > 0.0)
        {
            /* Get the lognormal ATM volatility */
            err = mkt->get_cash_vol_convert(
                mkt->get_cash_vol,
                mkt->vc,
                cpn->pay_start_date,
                cpn->pay_theo_end_date,
                cpn->pay_fix_time,
                pay_swp_cash,
                pay_swp_cash,
                mkt->ref,
                SRT_LOGNORMAL,
                &pay_atmlogvol);

            if (err)
                goto FREE_RETURN;

            if (float_adjust_type == 2)
            {
                pay_atmlogstd = pay_atmlogvol * sqrt(cpn->pay_fix_time);
            }
        }
        else
        {
            pay_atmlogvol = 0.0;
            pay_atmlogstd = 0.0;
        }
    }
    else
    {
        pay_swp_drs   = 0.0;
        pay_swp_fwd   = 0.0;
        pay_atmlogvol = 0.0;
        pay_atmlogstd = 0.0;
    }

    /* value each fixing */
    for (i = 0; i < cpn->used_nfix; i++)
    {
        fix = &(cpn->fix[i]);

        if (cpn->type == 3)
        {
            fixing_value = 0.0;
        }
        else
        {
            fixing_value = cpn->alpha;
        }

        fixing_value_mrg = cpn->alpha;

        ref_delay = cpn->cpn_pay_time - fix->schedule.dCpnTime[0];

        /* Calculate the forward */
        ref_swp_cash = fix->ref_swp_cash;
        ref_swp_fwd  = ref_swp_cash + fix->ref_fwd_spread;

        /* Calculate the DRS adjusted forward */

        if (cpn->beta || cpn->type == 1 || cpn->type == 2 ||
            ((cpn->type == 0 || cpn->type == 3 || cpn->type == 4) && cpn->nstr &&
             (!use_cmsopt || cpn->type == 4)))
        {
            err = swp_f_Cms_Rate(
                ref_swp_fwd,
                fix->ref_fix_time,
                fix->schedule.iNCpn - 1,
                fix->ref_cpnd,
                ref_delay,
                fix->ref_conv,
                mkt->vol_type,
                0.0,
                1,
                fix->schedule.lCpnDate[0],
                fix->ref_theo_end_date,
                (SRT_Boolean)mkt->cash_vol,
                fix->ref_fwd_spread,
                mkt->vc,
                mkt->num_strikes_in_vol,
                mkt->strikes_in_vol,
                &ref_swp_drs);

            if (err)
                goto FREE_RETURN;

            fixing_value += cpn->beta * ref_swp_drs;
        }

        /* Fix coupon valuation */
        if (cpn->type == 0 || do_margin)
        {
            for (j = 0; j < cpn->nstr; j++)
            {
                if (use_cmsopt)
                {
                    /* use CmsOption function */
                    err = swp_f_Cms_Option(
                        ref_swp_fwd,
                        fix->ref_fix_time,
                        fix->schedule.iNCpn - 1,
                        cpn->str[j],
                        fix->ref_cpnd,
                        SRT_PAYER,
                        ref_delay,
                        fix->ref_conv,
                        mkt->vol_type,
                        0.0,
                        1,
                        fix->schedule.lCpnDate[0],
                        fix->ref_theo_end_date,
                        (SRT_Boolean)mkt->cash_vol,
                        fix->ref_fwd_spread,
                        mkt->vc,
                        mkt->num_strikes_in_vol,
                        mkt->strikes_in_vol,
                        &opt);

                    if (err)
                        goto FREE_RETURN;
                }
                else
                {
                    if (cpn->str[j] - fix->ref_fwd_spread > 1.0E-08)
                    {
                        /* Get the lognormal vol of the ref rate */
                        err = mkt->get_cash_vol_convert(
                            mkt->get_cash_vol,
                            mkt->vc,
                            fix->schedule.lCpnDate[0],
                            fix->ref_theo_end_date,
                            fix->ref_fix_time,
                            ref_swp_cash,
                            cpn->str[j] - fix->ref_fwd_spread,
                            mkt->ref,
                            SRT_LOGNORMAL,
                            &ref_strikelogvol);

                        if (err)
                            goto FREE_RETURN;

                        /* Price the option with BS */
                        opt = srt_f_optblksch(
                            ref_swp_drs - fix->ref_fwd_spread,
                            cpn->str[j] - fix->ref_fwd_spread,
                            ref_strikelogvol,
                            fix->ref_fix_time,
                            1.0,
                            SRT_CALL,
                            (SrtGreekType)SRT_PREMIUM);
                    }
                    else
                    {
                        ref_strikelogvol = 1.0E-08;
                        opt              = ref_swp_drs - cpn->str[j];
                    }
                }

                if (cpn->type == 0 || cpn->type == 3 || cpn->type == 4)
                {
                    if (cpn->type == 3)
                    {
                        /* We seperate the value into two parts */
                        if (use_cmsopt)
                        {
                            /* use CmsOption function */
                            err = swp_f_Cms_Option(
                                ref_swp_fwd,
                                fix->ref_fix_time,
                                fix->schedule.iNCpn - 1,
                                cpn->str[j] * 1.01,
                                fix->ref_cpnd,
                                SRT_PAYER,
                                ref_delay,
                                fix->ref_conv,
                                mkt->vol_type,
                                0.0,
                                1,
                                fix->schedule.lCpnDate[0],
                                fix->ref_theo_end_date,
                                (SRT_Boolean)mkt->cash_vol,
                                fix->ref_fwd_spread,
                                mkt->vc,
                                mkt->num_strikes_in_vol,
                                mkt->strikes_in_vol,
                                &delta);

                            if (err)
                                goto FREE_RETURN;

                            delta = 1.0 + (opt - delta) / (0.01 * cpn->str[j]);
                        }
                        else
                        {
                            if (cpn->str[j] - fix->ref_fwd_spread > 1.0E-08)
                            {
                                /* Price the option with BS */
                                delta = srt_f_optblksch(
                                    ref_swp_drs - fix->ref_fwd_spread,
                                    cpn->str[j] - fix->ref_fwd_spread,
                                    ref_strikelogvol,
                                    fix->ref_fix_time,
                                    1.0,
                                    SRT_CALL,
                                    DELTA);
                            }
                            else
                            {
                                delta = 1.0;
                            }
                        }

                        fixing_value += cpn->nbopt[j] * delta * ref_swp_drs;
                        fixing_value_mrg += cpn->nbopt[j] * (opt - delta * ref_swp_drs);
                    }
                    else
                    {
                        fixing_value += cpn->nbopt[j] * opt;
                    }
                }
                else
                {
                    fixing_value_mrg += cpn->nbopt[j] * opt;
                }
            }
        }

        /* Floating coupon valuation */
        if (cpn->type == 1 || cpn->type == 2)
        {
            /* Precalculation of the ATM Ref log vol */
            if (float_adjust_type == 0 || float_adjust_type == 2)
            {
                /* we need the ATM lognormal volatility */
                err = mkt->get_cash_vol_convert(
                    mkt->get_cash_vol,
                    mkt->vc,
                    fix->schedule.lCpnDate[0],
                    fix->ref_theo_end_date,
                    fix->ref_fix_time,
                    ref_swp_cash,
                    ref_swp_cash,
                    mkt->ref,
                    SRT_LOGNORMAL,
                    &ref_atmlogvol);

                if (err)
                    goto FREE_RETURN;

                ref_atmlogstd = ref_atmlogvol * sqrt(fix->ref_fix_time);
            }

            /* Calculate the correlation */
            ref_correl = (correl_end - correl_start) / pay_regpaydelay *
                             (fix->schedule.dCpnTime[0] - cpn->cpn_start_time) +
                         correl_start;

            /* Valuation of each fixing */
            for (j = 0; j < cpn->nstr; j++)
            {
                /* Get the lognormal vol of the ref rate */
                err = mkt->get_cash_vol_convert(
                    mkt->get_cash_vol,
                    mkt->vc,
                    fix->schedule.lCpnDate[0],
                    fix->ref_theo_end_date,
                    fix->ref_fix_time,
                    ref_swp_cash,
                    cpn->str[j] - fix->ref_fwd_spread,
                    mkt->ref,
                    SRT_LOGNORMAL,
                    &ref_strikelogvol);

                if (err)
                    goto FREE_RETURN;

                if (j % 2 == 0)
                {
                    switch (float_adjust_type)
                    {
                    case 0:
                    {
                        /* adjustment ATM / ATM */
                        pay_adjlogvol = pay_atmlogvol;
                        ref_adjlogvol = ref_atmlogvol;

                        break;
                    }
                    case 1:
                    {
                        /* adjustment ATM / ATS */
                        pay_adjlogvol = pay_atmlogvol;
                        ref_adjlogvol = ref_strikelogvol;
                        break;
                    }
                    case 2:
                    {
                        /* adjustment STD / ATS */
                        ref_adjlogvol  = ref_strikelogvol;
                        ref_strikestd  = log(cpn->str[j] / ref_swp_fwd) / ref_atmlogstd;
                        pay_new_strike = pay_swp_fwd * exp(ref_strikestd * pay_atmlogstd);

                        /* Get the lognormal vol of the ref rate */
                        err = mkt->get_cash_vol_convert(
                            mkt->get_cash_vol,
                            mkt->vc,
                            cpn->pay_start_date,
                            cpn->pay_theo_end_date,
                            cpn->pay_fix_time,
                            pay_swp_cash,
                            pay_new_strike - cpn->pay_fwd_spread,
                            mkt->ref,
                            SRT_LOGNORMAL,
                            &pay_adjlogvol);

                        if (err)
                            goto FREE_RETURN;

                        break;
                    }
                    case 3:
                    {
                        /* adjustment ATS / ATS */
                        ref_adjlogvol = ref_strikelogvol;

                        /* Get the lognormal vol of the ref rate */
                        err = mkt->get_cash_vol_convert(
                            mkt->get_cash_vol,
                            mkt->vc,
                            cpn->pay_start_date,
                            cpn->pay_theo_end_date,
                            cpn->pay_fix_time,
                            pay_swp_cash,
                            cpn->str[j] - cpn->pay_fwd_spread,
                            mkt->ref,
                            SRT_LOGNORMAL,
                            &pay_adjlogvol);

                        if (err)
                            goto FREE_RETURN;

                        break;
                    }
                    }

                    /* Calculate the adjustment in the ref numeraire */
                    if (cpn->type == 1)
                    {
                        ref_adj =
                            exp(ref_correl * ref_adjlogvol * pay_adjlogvol * cpn->pay_fix_time);
                    }
                    else
                    {
                        ref_adj =
                            exp(ref_correl * ref_adjlogvol * pay_adjlogvol * fix->ref_fix_time);
                    }
                }

                /* Price the option */
                if (cpn->str[j] - fix->ref_fwd_spread > 1.0E-08)
                {
                    opt = srt_f_optblksch(
                        (ref_swp_drs - fix->ref_fwd_spread) * ref_adj,
                        cpn->str[j] - fix->ref_fwd_spread,
                        ref_strikelogvol,
                        fix->ref_fix_time,
                        1.0,
                        SRT_CALL,
                        (SrtGreekType)SRT_PREMIUM);
                }
                else
                {
                    opt = (ref_swp_drs - fix->ref_fwd_spread) * ref_adj -
                          (cpn->str[j] - fix->ref_fwd_spread);
                }

                fixing_value += cpn->nbopt[j] * opt;
            }
        }

        cpn_value += cpn->used_fixweights[i] * fixing_value;
        cpn_value_mrg += cpn->used_fixweights[i] * fixing_value_mrg;
    }

    if (cpn->type == 1 || cpn->type == 2)
    {
        cpn_value *= pay_swp_drs * cpn->pay_gearing;
        cpn_value_mrg *= cpn->cpn_coupon;
    }
    else if (cpn->type == 0)
    {
        cpn_value *= cpn->cpn_coupon;
        cpn_value_mrg = cpn_value;
        cpn_value     = 0.0;
    }

    df_pay = swp_f_df(mkt->today, cpn->cpn_pay_date, mkt->yc);

    temp                = cpn->cpn_cvg * cpn->cpn_not * df_pay;
    cpn->mkt_fixed_part = cpn_value_mrg * temp;
    cpn->mkt_float_part = cpn_value * temp;

    cpn->mkt_val = cpn->mkt_fixed_part + cpn->mkt_float_part;

FREE_RETURN:

    return err;
}

/* Find equivalent cash strikes */
Err cts_calc_equivalent_strikes(CTS cts_iv, CTS cts_call, CTS_MKT mkt)
{
    Err    err = NULL;
    int    i, j, last_exo, last_fund;
    double exo_pv, fund_pv, level_exo, swap_pv;
    double last_df;
    long   end_exo;

    last_exo  = cts_iv->exo_leg->num_cpn;
    last_fund = cts_iv->fund_leg->num_cpn;

    last_df = swp_f_df(mkt->today, cts_iv->fund_leg->cpn[last_fund - 1].pay_date, mkt->yc);
    end_exo = cts_iv->exo_leg->cpn[last_exo - 1].cpn_pay_date;

    exo_pv  = 0.0;
    fund_pv = last_df * cts_iv->fund_leg->cpn[last_fund - 1].not ;

    for (i = cts_call->num_calls - 1; i >= 00; i--)
    {
        for (j = cts_call->call[i].exo_idx; j < last_exo; j++)
        {
            exo_pv += cts_iv->exo_leg->cpn[j].mkt_val;
        }

        last_exo = cts_call->call[i].exo_idx;

        for (j = cts_call->call[i].fund_idx; j < last_fund; j++)
        {
            fund_pv += cts_iv->fund_leg->cpn[j].mkt_val;
        }

        last_fund = cts_call->call[i].fund_idx;

        swap_pv = exo_pv / cts_iv->exo_leg->cpn[last_exo].cpn_not -
                  fund_pv / cts_iv->fund_leg->cpn[last_fund].not ;

        err = swp_f_LevelPayment(
            cts_iv->exo_leg->cpn[last_exo].cpn_start_date,
            end_exo,
            mkt->swap_freq,
            mkt->swap_basis,
            mkt->yc,
            mkt->ref,
            &level_exo);

        if (err)
            goto FREE_RETURN;

        cts_call->call[i].implied_strike = swap_pv / level_exo;
    }

FREE_RETURN:

    return err;
}

/*	Check */
Err cts_check_coupon(CTS_CPN cpn)
{
    Err err = NULL;
    int i;

    if (cpn->type)
    {
        if (cpn->pay_fix_date > cpn->pay_start_date)
            return "Date inconsistency in cts_check_coupon";
        if (cpn->pay_start_date >= cpn->pay_end_date)
            return "Date inconsistency in cts_check_coupon";
    }

    if (cpn->cpn_start_date >= cpn->cpn_pay_date)
        return "Date inconsistency in cts_check_coupon";

    if (cpn->fix[cpn->used_nfix - 1].ref_fix_date > cpn->cpn_pay_date)
    {
        return "Fixing date inconsistency in cts_check_coupon";
    }

    /* Verify that strikes are sorted */
    for (i = 0; i < cpn->nstr - 1; i++)
    {
        if (cpn->str[i] >= cpn->str[i + 1])
        {
            return "Option strikes are not sorted in cts_check_coupon";
        }
    }

    if (cpn->fix)
    {
        for (i = 0; i < cpn->used_nfix; i++)
        {
            err = cts_check_fixing(&(cpn->fix[i]));
            if (err)
                return err;
        }
    }

    return NULL;
}

/*	Free */
void cts_free_coupon(CTS_CPN cpn)
{
    if (cpn)
    {
        if (cpn->weights)
            free(cpn->weights);
        cpn->weights = NULL;

        if (cpn->str)
            free(cpn->str);
        cpn->str = NULL;
        if (cpn->nbopt)
            free(cpn->nbopt);
        cpn->nbopt = NULL;

        if (cpn->fix)
        {
            free(cpn->fix);
            cpn->fix = NULL;
        }

        if (cpn->used_fixidx)
            free(cpn->used_fixidx);
        cpn->used_fixidx = NULL;
        if (cpn->used_fixweights)
            free(cpn->used_fixweights);
        cpn->used_fixweights = NULL;
    }
}

/*	Exotic leg */

/*	Init */
Err cts_fill_cpn_leg(
    CTS_MKT mkt,
    /*	Coupons that fixed before today are disregarded */
    /*	EOD Flag */
    int eod_flag, /*	0: I, 1: E */
    /*	General Libor properties */
    int fix_lag_bd,
    /*	Number of coupons */
    int ncpn,
    /*	Coupon description */
    /*		cpn */
    int* cpn_type, /*	0:	fixed
                                   1:	libor fixed at start
                                   2:	libor fixed at end */
    long*   cpn_start_date,
    long*   cpn_end_date,
    long*   cpn_pay_date,
    double* cpn_coupon,
    char**  cpn_basis,
    double* cpn_not,
    /*		pay libor */
    long*   pay_fix_date, /*	fixing date of the Libor */
    int     pay_months,
    char*   pay_freq,
    char*   pay_basis,
    double* pay_fwd_spread,
    double* pay_gearing,
    /*		fix libor */
    /*	fixings */
    int*     nfix,
    double** weights,
    long**   ref_fix_dates,
    char**   fix_tenor, /*	Fixing tenors */
    char*    fix_freq,
    char*    fix_basis,
    double** ref_fwd_spreads,
    /*	Profiles */
    int tstype, /*	0: generic, 1: range */
    /*	profile 0 */
    double*  alpha,
    double*  beta,
    int*     nstr,
    double** str,
    double** nbopt,
    /*	profile 1 */
    int     buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int     value_zero, /*	0: Do not value low strike options, 1: do */
    double* lb,         /*	Lower bounds of ranges */
    double* ub,         /*	Upper bounds of ranges */
    double* payoff,     /*	Payoff in ranges */
    double  call_spread,
    double  numer_spread,
    /*		Trimming */
    int trim_type, /*	0: no trim
                                   1: x fixings max
                                   2: x time min between two fixings */
    int    max_fix,
    double min_fix_time,
    int    calc_mkt_iv,
    int    use_cmsopt,        /*	Use CmsOption to value fix option, use BS on CmsRate otherwise */
    double correl_start,      /*	Correl between libor fixing and libor paid at start */
    double correl_end,        /*	Correl between libor fixing and libor paid at start */
    int    float_adjust_type, /*	type of adjustment for the floating coupon,
                                                      0: ATM vol, 1: Strike Vol */
    CTS_EXO_LEG exo_leg)
{
    Err  err = NULL;
    long today;
    int  i, j;

    /*	Get today */
    today = mkt->today;

    /*	Init array */
    exo_leg->cpn = NULL;

    /*	Skip past coupons */
    j = 0;
    while (j < ncpn && ref_fix_dates[j][0] < today + eod_flag)
    {
        j++;
    }
    /*	j: index of the first coupon to be considered */

    exo_leg->num_cpn = ncpn - j;
    if (exo_leg->num_cpn < 1)
    {
        /* err = "All funding coupons start before today in cts_fill_cpn_leg"; */
        exo_leg->num_cpn = 0;
        goto FREE_RETURN;
    }

    /*	Allocate */
    exo_leg->cpn = (cts_cpn*)calloc(exo_leg->num_cpn, sizeof(cts_cpn));
    if (!exo_leg->cpn)
    {
        err = "Allocation error in cts_fill_cpn_leg";
        goto FREE_RETURN;
    }

    /*	Fill coupons */
    for (i = 0; i < exo_leg->num_cpn; i++)
    {
        err = cts_init_cpn(
            mkt,
            cpn_start_date[j + i],
            cpn_end_date[j + i],
            cpn_pay_date[j + i],
            cpn_coupon[j + i],
            cpn_basis[j + i],
            cpn_not[j + i],
            fix_lag_bd,
            cpn_type[j + i],
            pay_fix_date[j + i],
            pay_months,
            pay_freq,
            pay_basis,
            pay_fwd_spread[j + i],
            pay_gearing[j + i],
            nfix[j + i],
            weights[j + i],
            ref_fix_dates[j + i],
            fix_tenor[j + i],
            fix_freq,
            fix_basis,
            ref_fwd_spreads[j + i],
            tstype,
            alpha[j + i],
            beta[j + i],
            nstr[j + i],
            str[j + i],
            nbopt[j + i],
            buy_sell,
            value_zero,
            lb[j + i],
            ub[j + i],
            payoff[j + i],
            call_spread,
            numer_spread,
            trim_type,
            max_fix,
            min_fix_time,
            calc_mkt_iv,
            use_cmsopt,
            correl_start,
            correl_end,
            float_adjust_type,
            &(exo_leg->cpn[i]));

        /*	Error: free the fixings already initialised */
        if (err)
        {
            for (j = 0; j < i; j++)
            {
                cts_free_coupon(&(exo_leg->cpn[j]));
            }
            free(exo_leg->cpn);
            exo_leg->cpn = NULL;

            goto FREE_RETURN;
        }
    }

    err = cts_check_exo_leg(exo_leg);

FREE_RETURN:

    if (err)
    {
        cts_free_exo_leg(exo_leg);
    }

    return err;
}

Err cts_check_exo_leg(CTS_EXO_LEG exo_leg)
{
    int i;
    Err err = NULL;

    if (exo_leg->cpn)
    {
        for (i = 0; i < exo_leg->num_cpn; i++)
        {
            err = cts_check_coupon(&(exo_leg->cpn[i]));
            if (err)
                return err;
        }
    }
    else
    {
        return "Coupons not allocated in exo leg in cts_check_exo_leg";
    }

    return err;
}

void cts_free_exo_leg(CTS_EXO_LEG exo_leg)
{
    int i;

    if (exo_leg->cpn)
    {
        for (i = 0; i < exo_leg->num_cpn; i++)
        {
            cts_free_coupon(&(exo_leg->cpn[i]));
        }
        free(exo_leg->cpn);
        exo_leg->cpn = NULL;
    }
}

/*	Functions for the calls */

/*	Init */
Err cts_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    int     ncall,
    int     pay_rec, /*	1: rec pd, -1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,
    CTS     cts)
{
    CTS_FUND_LEG fund_leg;
    CTS_EXO_LEG  exo_leg;
    CTS_CALL     call;
    int          i, j, k, l, is_midat;
    Err          err = NULL;

    cts->call = NULL;
    fund_leg  = cts->fund_leg;
    exo_leg   = cts->exo_leg;

    /*	Skip calls to be exercised before today */
    i = 0;
    while (i < ncall && ex_date[i] < today + eod_flag)
    {
        i++;
    }

    /*	Check that at least one call is left */
    if (i == ncall)
    {
        /* err = "All calls are to be exercised before today in cts_fill_calls"; */
        cts->num_calls = 0;
        goto FREE_RETURN;
    }

    /*	Allocate memory */
    cts->num_calls = ncall - i;
    cts->call      = (cts_call*)calloc(cts->num_calls, sizeof(cts_call));
    if (!cts->num_calls)
    {
        err = "Allocation error in cts_fill_calls";
        goto FREE_RETURN;
    }

    /*	Midat Flag */
    is_midat = 1;

    for (j = cts->call[0].exo_idx; j < exo_leg->num_cpn; j++)
    {
        if (exo_leg->cpn[j].nstr > 0 || fabs(exo_leg->cpn[j].beta) > 1.0E-08 ||
            exo_leg->cpn[j].type > 0)
        {
            is_midat = 0;
            break;
        }
    }

    /*	Fill call information */
    j = 0;
    k = 0;
    l = 0;

    while (i < ncall)
    {
        call        = cts->call + j;
        call->index = j;

        /*	Dates */
        call->ex_date  = ex_date[i];
        call->set_date = set_date[i];

        /*	Times */
        call->ex_time  = (call->ex_date - today) * YEARS_IN_DAY;
        call->set_time = (call->set_date - today) * YEARS_IN_DAY;

        /*	Call on funding leg */
        /*	k = index of the first coupon to be called on funding leg,
                        i.e. first coupon with a start date >= ex date */

        while (k < fund_leg->num_cpn && fund_leg->cpn[k].start_date < call->ex_date)
        {
            k++;
        }
        if (k == fund_leg->num_cpn)
        {
            err = serror("Call number %d does not control any coupon in funding leg", i);
            goto FREE_RETURN;
        }
        call->fund_idx             = k;
        call->num_fund_cpn         = fund_leg->num_cpn - k;
        call->num_partial_fund_cpn = call->num_fund_cpn;

        if (j > 0)
        {
            (call - 1)->num_partial_fund_cpn = (call - 1)->num_fund_cpn - call->num_fund_cpn;
        }

        /*	Call on exotic leg */
        /*	k = index of the first coupon to be called on exotic leg,
                        i.e. first coupon with a start date >= ex date */
        while (l < exo_leg->num_cpn && exo_leg->cpn[l].cpn_start_date < call->ex_date)
        {
            l++;
        }
        if (l == exo_leg->num_cpn)
        {
            err = serror("Call number %d does not control any coupon in exotic leg", i);
            goto FREE_RETURN;
        }

        call->exo_idx             = l;
        call->num_exo_cpn         = exo_leg->num_cpn - l;
        call->num_partial_exo_cpn = call->num_exo_cpn;

        if (j > 0)
        {
            (call - 1)->num_partial_exo_cpn = (call - 1)->num_exo_cpn - call->num_exo_cpn;
        }

        /*	Payer or receiver */
        call->pay_rec = pay_rec;

        call->is_midat = is_midat;

        /*	Fee */
        call->fee = fee[i];

        i++;
        j++;
    }

    err = cts_check_calls(cts);

    if (err)
        goto FREE_RETURN;

FREE_RETURN:

    if (err)
    {
        cts_free_calls(cts);
    }

    return err;
}

/* Adjustment for fixing in floating case */
Err cts_adjust_floating_fixing(CTS cts, int pde_or_mc)
{
    int i;

    /*	Check that ex and set dates are strictly increasing
                    Also check that funding and pd indices are strictly increasing,
                    i.e. there is no redundant calls */
    if (pde_or_mc == 1)
    {
        for (i = 0; i < cts->num_calls; i++)
        {
            if (cts->call[i].exo_idx > 0)
            {
                /* look at the fixing of the floating rate before */
                if (cts->exo_leg->cpn[cts->call[i].exo_idx - 1].type == 2)
                {
                    if (cts->exo_leg->cpn[cts->call[i].exo_idx - 1].pay_fix_date >=
                        cts->call[i].ex_date)
                    {
                        cts->exo_leg->cpn[cts->call[i].exo_idx - 1].pay_fix_date =
                            cts->call[i].ex_date - 1;
                        cts->exo_leg->cpn[cts->call[i].exo_idx - 1].pay_fix_time =
                            cts->call[i].ex_time - YEARS_IN_DAY;
                    }
                }
            }
        }
    }

    /*	OK */
    return NULL;
}

/*	Check dates consistency */
Err cts_check_calls(CTS cts)
{
    int i;

    /*	Check that ex and set dates are strictly increasing
                    Also check that funding and pd indices are strictly increasing,
                    i.e. there is no redundant calls */
    for (i = 1; i < cts->num_calls; i++)
    {
        if (cts->call[i].ex_date <= cts->call[i - 1].ex_date)
        {
            return "Exercise dates should be increasing";
        }

        if (cts->call[i].set_date <= cts->call[i - 1].set_date)
        {
            return "Settlement dates should be increasing";
        }

        if (cts->call[i].fund_idx < cts->call[i - 1].fund_idx)
        {
            return "Number of funding coupons controlled by calls should be decreasing";
        }

        if (cts->call[i].exo_idx < cts->call[i - 1].exo_idx)
        {
            return "Number of exotic coupons controlled by calls should be decreasing";
        }

        if (cts->call[i].fund_idx <= cts->call[i - 1].fund_idx &&
            cts->call[i].exo_idx <= cts->call[i - 1].exo_idx)
        {
            return serror("Calls %d and %d -indexed after today- are redundant", i - 1, i);
        }

        /* Check that we are not going to call coupon from the previous period */
        /* It returns an error for now by safety, later it will be handled properly */
        if (cts->call[i].exo_idx > 0)
        {
            if (cts->exo_leg->cpn[cts->call[i].exo_idx - 1].used_nfix > 0 &&
                cts->exo_leg->cpn[cts->call[i].exo_idx - 1]
                        .fix[cts->exo_leg->cpn[cts->call[i].exo_idx - 1].used_nfix - 1]
                        .ref_fix_date > cts->call[i].ex_date)
            {
                return serror(
                    "Either too many fixings requested, or Floating Rate fixe after call date");
            }
        }
    }

    /*	Check that set dates are after ex dates
                    Also check that the call date is before the start, end and fixing dates
                    of the coupons it controls */
    for (i = 0; i < cts->num_calls; i++)
    {
        if (cts->call[i].set_date < cts->call[i].ex_date)
        {
            return "Settlement dates should be after exercise dates";
        }

        if (cts->fund_leg->cpn[cts->call[i].fund_idx].start_date < cts->call[i].ex_date ||
            cts->fund_leg->cpn[cts->call[i].fund_idx].pay_date < cts->call[i].ex_date)
        {
            return "A funding coupon starts before its exercise date";
        }

        if (cts->exo_leg->cpn[cts->call[i].exo_idx].cpn_start_date < cts->call[i].ex_date ||
            cts->exo_leg->cpn[cts->call[i].exo_idx].cpn_pay_date < cts->call[i].ex_date ||
            cts->exo_leg->cpn[cts->call[i].exo_idx].fix[0].ref_fix_date < cts->call[i].ex_date)
        {
            return "An exotic coupon starts or fixes before its exercise date";
        }
    }

    /*	OK */
    return NULL;
}

/*	Free */
void cts_free_calls(CTS cts)
{
    if (cts->call)
    {
        free(cts->call);
        cts->call = NULL;
    }
}

/*	Functions for the underlying */

/*	Fill from existing */
Err cts_fill_und(CTS_MKT mkt, char* lgmsvund, double tstar, CTS_UND und)
{
    Err err = NULL;

    /* Initialisation */
    init_NULL_LGMSV_model(&(und->model));

    und->mkt = mkt;
    strcpy(und->name, lgmsvund);

    /* Get the underlying */
    err = Get_LGMSV_model(lgmsvund, &(und->model));

    if (err)
    {
        free_LGMSV_model(&(und->model));
        return err;
    }

    und->model.dInitLambdaEps = und->model.dLambdaEps[0];

    return err;
}

Err cts_get_mkt_sabr_beta_param(
    CTS_MKT mkt,
    long    exe_date,
    long    start_date,
    long    theo_end_date,
    double  fixed_beta,
    int     use_levenberg,
    double* alpha,
    double* rho,
    double* fitting_error)
{
    Err           err = NULL;
    double        atm_log_vol;
    static double Strikes[6], Vols[6];
    int           i;
    long          today = mkt->today;

    double swap_cash, maturity, atm_nor_std, atm_log_std, sigma_beta;
    long   theo_date, act_date, last_act_date;
    double level, df_start, df_end;
    int    nb_months;

    /* first calculate the forward */
    nb_months = 12 / mkt->swap_srt_freq;

    theo_date     = theo_end_date;
    act_date      = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    last_act_date = act_date;

    df_end   = swp_f_df(today, act_date, mkt->yc);
    df_start = df_end;

    level = 0.0;

    theo_date = add_unit(theo_date, -nb_months, SRT_MONTH, NO_BUSDAY_CONVENTION);
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);

    i = 1;

    while (act_date >= start_date - 5)
    {
        level += df_start * coverage(act_date, last_act_date, mkt->swap_srt_basis);
        df_start      = swp_f_df(today, act_date, mkt->yc);
        last_act_date = act_date;

        theo_date = add_unit(theo_date, -nb_months, SRT_MONTH, NO_BUSDAY_CONVENTION);
        act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);

        i++;
    }

    if (i == 1)
    {
        /* this is a short option */
        act_date = start_date;
        level += df_start * coverage(act_date, last_act_date, mkt->swap_srt_basis);
        df_start = swp_f_df(today, act_date, mkt->yc);
    }

    swap_cash = (df_start - df_end) / level;

    /* get volatilities */

    maturity = (exe_date - today) * YEARS_IN_DAY;

    err = mkt->get_cash_vol_convert(
        mkt->get_cash_vol,
        mkt->vc,
        start_date,
        theo_end_date,
        maturity,
        swap_cash,
        swap_cash,
        mkt->ref,
        SRT_LOGNORMAL,
        &atm_log_vol);

    if (err)
        return err;

    atm_log_std = atm_log_vol * sqrt(maturity);

    err = mkt->get_cash_vol_convert(
        mkt->get_cash_vol,
        mkt->vc,
        start_date,
        theo_end_date,
        maturity,
        swap_cash,
        swap_cash,
        mkt->ref,
        SRT_NORMAL,
        &atm_nor_std);

    if (err)
        return err;

    atm_nor_std *= sqrt(maturity);

    /* Get the market smile */

    Strikes[1] = swap_cash * exp(-1.0 * atm_log_std);
    Strikes[2] = swap_cash * exp(-0.5 * atm_log_std);
    Strikes[3] = swap_cash;
    Strikes[4] = swap_cash + 0.5 * atm_nor_std;
    Strikes[5] = swap_cash + 1.0 * atm_nor_std;

    for (i = 1; i < 6; i++)
    {
        err = mkt->get_cash_vol_convert(
            mkt->get_cash_vol,
            mkt->vc,
            start_date,
            theo_end_date,
            maturity,
            swap_cash,
            Strikes[i],
            mkt->ref,
            SRT_LOGNORMAL,
            &Vols[i]);

        if (err)
            return err;
    }

    /* Get the SABR parameters */
    if (use_levenberg)
    {
        err = opsabrcalib(
            swap_cash,
            maturity,
            5,
            Strikes,
            Vols,
            &atm_log_vol,
            alpha,
            1,
            &fixed_beta,
            0,
            rho,
            1,
            fitting_error);
    }
    else
    {
        err = calib_sabr_rr_bt_given_beta(
            swap_cash,
            maturity,
            atm_log_vol,
            Strikes[1],
            Vols[1],
            Strikes[5],
            Vols[5],
            &sigma_beta,
            alpha,
            fixed_beta,
            rho,
            SRT_LOGNORMAL,
            1,
            0.00001,
            7,
            fitting_error);
    }

    if (err)
        return err;

    return err;
}

/* Implied Alpha^2 = Alpha0^2 * (1.0 - exp(-2.0 * LamEps * T)) / (2.0 * LamEps * T) */
Err cts_implied_alpha_approx(
    double  maturity,
    double  param[], /* first param = Alpha^2, second param = 2.0 * LamEps */
    double* price,
    double* gradient,
    int     nb_param)
{
    double expcoef;
    double meancoef;

    if (param[2] < -1.0E-08 || param[1] < -1.0E-08)
    {
        /* out of bounds */
        *price      = -100000;
        gradient[1] = 100000;
        gradient[2] = 100000;
        return NULL;
    }

    if (param[2] > 1.0E-08)
    {
        meancoef    = param[2] * maturity;
        expcoef     = (1.0 - exp(-meancoef)) / meancoef;
        *price      = param[1] * expcoef;
        gradient[1] = expcoef;
        gradient[2] = (param[1] * exp(-meancoef) - *price) / param[2];
    }
    else
    {
        *price      = param[1];
        gradient[1] = 1.0;
        gradient[2] = 0.0;
    }

    return NULL;
}

Err cts_calib_const_alpha_and_rho(
    CTS_MKT mkt,
    CTS     cts,
    int calib_months, /* 0: co-terminal swaption, otherwise underlyings with required nb months */
    double               fixed_beta,
    double               min_time,
    double*              alpha,
    double*              rho,
    double*              lameps,
    int                  save_inst_data,
    cpd_calib_inst_data* inst_data)
{
    Err     err = NULL;
    int     i, i0, j, k;
    double *mkt_alpha = NULL, *mkt_rho = NULL;

    long exercise_dates_[MAX_INST];

    long *  exercise_dates = NULL, *start_dates = NULL, *end_dates = NULL;
    double *exercise_times = NULL, *weights = NULL;

    long   start_date, cts_end_date, theo_end_date;
    int    num_ex_date;
    double fitting_error;
    double param[3];
    int    use_param[3];
    double temp;

    /* first select the calibration dates */

    if (calib_months == 0)
    {
        /* Remove short marturities */
        i0 = 0;

        while (i0 < cts->num_calls && cts->call[i0].ex_time < MIN_CALIB_TIME)
        {
            i0++;
        }

        num_ex_date = cts->num_calls - i0;

        for (i = 0; i < num_ex_date; i++)
        {
            exercise_dates_[i] = cts->call[i + i0].ex_date;
        }
    }
    else
    {
        /* Remove short marturities */
        i0 = cts->call[0].exo_idx;

        while (i0 < cts->exo_leg->num_cpn &&
               cts->exo_leg->cpn[cts->call[0].exo_idx + i0].cpn_start_time < MIN_CALIB_TIME)
        {
            i0++;
        }

        num_ex_date = cts->exo_leg->num_cpn - i0;

        for (i = 0; i < num_ex_date; i++)
        {
            exercise_dates_[i] = cts->exo_leg->cpn[i0 + i].cpn_start_date;
            exercise_dates_[i] = add_unit(exercise_dates_[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
    }

    j              = num_ex_date - 1;
    exercise_dates = &(exercise_dates_[0]);

    for (i = num_ex_date - 2; i >= 0; i--)
    {
        if ((exercise_dates[j] - exercise_dates[i]) * YEARS_IN_DAY < min_time - ONE_MONTH)
        {
            for (k = i - 1; k >= 0; k--)
            {
                exercise_dates[k + 1] = exercise_dates[k];
            }

            exercise_dates++;
            num_ex_date--;
            j--;
            if (num_ex_date < 1)
            {
                err = "All exercise dates are past in cts_calib_const_alpha_and_rho";
                goto FREE_RETURN;
            }
        }
        else
        {
            j = i;
        }
    }

    mkt_alpha      = (double*)calloc(num_ex_date, sizeof(double));
    mkt_rho        = (double*)calloc(num_ex_date, sizeof(double));
    exercise_times = (double*)calloc(num_ex_date, sizeof(double));
    weights        = (double*)calloc(num_ex_date, sizeof(double));
    start_dates    = (long*)calloc(num_ex_date, sizeof(long));
    end_dates      = (long*)calloc(num_ex_date, sizeof(long));

    if (!mkt_alpha || !mkt_rho || !exercise_times || !weights || !start_dates || !end_dates)
    {
        err = "Memory allocation faillure (2) in cts_calib_const_alpha_and_rho";
        goto FREE_RETURN;
    }

    /* compute the mkt implied SABR params */

    cts_end_date = cts->fund_leg->cpn[cts->fund_leg->num_cpn - 1].pay_date;

    for (i = 0; i < num_ex_date; i++)
    {
        exercise_times[i] = (exercise_dates[i] - mkt->today) * YEARS_IN_DAY;
        start_date        = add_unit(exercise_dates[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING);

        if (calib_months == 0)
        {
            theo_end_date = cts_end_date;
        }
        else
        {
            theo_end_date = add_unit(start_date, calib_months, SRT_MONTH, NO_BUSDAY_CONVENTION);
        }

        start_dates[i] = start_date;
        end_dates[i]   = theo_end_date;

        /* First guess for Alpha and Rho */
        if (i == 0)
        {
            mkt_alpha[i] = LGMSV_DEFAULT_ALPHAEPS;
            mkt_rho[i]   = LGMSV_DEFAULT_RHOEPS;
        }
        else
        {
            mkt_alpha[i] = mkt_alpha[i - 1];
            mkt_rho[i]   = mkt_rho[i - 1];
        }

        err = cts_get_mkt_sabr_beta_param(
            mkt,
            exercise_dates[i],
            start_date,
            theo_end_date,
            fixed_beta,
            1,
            &mkt_alpha[i],
            &mkt_rho[i],
            &fitting_error);

        if (err)
            goto FREE_RETURN;
    }

    if (num_ex_date > 1)
    {
        /* we take rho as the minimum of all the rhos */

        temp = 0.0;

        for (i = 0; i < num_ex_date; i++)
        {
            temp += mkt_rho[i];
        }

        *rho = temp / num_ex_date;

        /* we fit alpha and mean-reversion to the mkt */

        /* rescale the target */
        for (i = 0; i < num_ex_date; i++)
        {
            mkt_alpha[i] *= mkt_alpha[i];
            weights[i] = 1.0;
        }

        /* first guess */
        param[1]     = mkt_alpha[0];
        param[2]     = 2.0 * 0.1;
        use_param[1] = 1;
        use_param[2] = 1;

        err = levenberg_marquardt_select(
            exercise_times - 1,
            mkt_alpha - 1,
            weights - 1,
            num_ex_date,
            param,
            (long*)use_param,
            2,
            LGMSV_ALPHA_MAXITER,
            cts_implied_alpha_approx,
            &fitting_error);

        if (err)
            goto FREE_RETURN;

        *alpha  = sqrt(param[1]);
        *lameps = param[2] / 2.0;
    }
    else
    {
        if (num_ex_date == 1)
        {
            *alpha  = mkt_alpha[0];
            *lameps = LGMSV_DEFAULT_LAMEPS;
            *rho    = mkt_rho[0];
        }
        else
        {
            *alpha  = LGMSV_DEFAULT_ALPHAEPS;
            *lameps = LGMSV_DEFAULT_LAMEPS;
            *rho    = LGMSV_DEFAULT_RHOEPS;
        }
    }

    if (save_inst_data)
    {
        /* Memory allocation */
        inst_data->num_inst_smile = num_ex_date;

        if (num_ex_date > 0)
        {
            inst_data->start_dates_smile = (long*)calloc(inst_data->num_inst_smile, sizeof(long));
            inst_data->end_dates_smile   = (long*)calloc(inst_data->num_inst_smile, sizeof(long));
            inst_data->alpha = (double*)calloc(inst_data->num_inst_smile, sizeof(double));
            inst_data->rho   = (double*)calloc(inst_data->num_inst_smile, sizeof(double));

            if (!inst_data->start_dates_smile || !inst_data->end_dates_smile || !inst_data->alpha ||
                !inst_data->rho)
            {
                err = "Memory allocation faillure in cts_calib_const_alpha_and_rho";
                goto FREE_RETURN;
            }

            for (i = 0; i < num_ex_date; i++)
            {
                inst_data->start_dates_smile[i] = start_dates[i];
                inst_data->end_dates_smile[i]   = end_dates[i];
                inst_data->alpha[i]             = sqrt(mkt_alpha[i]);
                inst_data->rho[i]               = mkt_rho[i];
            }
        }
    }

FREE_RETURN:

    if (exercise_times)
        free(exercise_times);
    if (weights)
        free(weights);
    if (start_dates)
        free(start_dates);
    if (end_dates)
        free(end_dates);
    if (mkt_alpha)
        free(mkt_alpha);
    if (mkt_rho)
        free(mkt_rho);

    return err;
}

/*	Calibrate */
/*	Fill underlying structure from calibration instruments */
Err cts_calib_und(
    CTS_MKT mkt,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    double  lambda,   /*	lambda */
    int     nb_factor,
    double  lgm_alpha,
    double  lgm_gamma,
    double  lgm_rho,
    int     nsmilepar,
    double* smilepartime,
    double* alpha,  /*	alpha */
    double* rho,    /*	rho */
    double* ldaeps, /*	ldaeps */
    double* rho2,   /*	rho2eps */
    double  tstar,  /*	TStar of the model */

    /*	Numerical CF params */
    LGMSV_NUMERPARAMS NumerParams,
    /*	Calib params */
    char* cal_tenor,
    char* cal_ref,
    char* cal_freq,
    char* cal_basis,

    int                 force_atm, /*	force atm calib */
    double              max_std_long,
    double              max_std_short,
    double              vol_shift_long,
    DIAGCALIB_VOLTYPE   vol_type_long,
    DIAGCALIB_SHIFTTYPE vol_shift_type_long,
    double              vol_shift_short,
    DIAGCALIB_VOLTYPE   vol_type_short,
    DIAGCALIB_SHIFTTYPE vol_shift_type_short,
    double              lambda_shift,

    int   calib_strategy, /*	0: autocal, 1: swaptions / cap, 2: cap / swaptions */
    int   fix_lambda,     /*	0: calib lambda to cap, 1: fix lambda calib to diagonal */
    char* short_tenor,
    char* short_refrate,
    char* short_freq,
    char* short_basis,
    int   fix_smile,          /*	0: calib smile parameters to market smile */
    int   smile_calib_months, /*	0: co-terminal swaption, otherwise underlyings with required nb
                                 months */
    LGMSV_CalibParams* lgmsv_calib_params,
    double             min_time,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps, /*	1: we allow jumps on vol, 0: we don't */
    double prec,
    int    maxiter,
    int    keep_first,
    /*	Strike choice */
    int long_strike_flag,  /*	0: ATM
                                                   1: Coupon
                                                   2: Eq (PV/Lvl) */
    int short_strike_flag, /*	0: ATM,
                                                   1: implied digital caplet strike
                                                   2: same number of std */
    /*	End of calib params */
    CTS                  cts, /*	structure */
    CTS_UND              und,
    int                  save_inst_data,
    cpd_calib_inst_data* inst_data)
{
    int*    cal_dates    = NULL;
    long*   ex_dates     = NULL;
    double *prim_strikes = NULL, *prim_strikesS1 = NULL, *prim_strikesS2 = NULL,
           *sec_strikes = NULL, *sec_strikesS1 = NULL, *sec_strikesS2 = NULL;

    char **prim_tenor = NULL, **sec_tenor = NULL;

    int     used_nsmilepar;
    double *used_smilepartime = NULL, *used_alpha = NULL, *used_rho = NULL, *used_ldaeps = NULL,
           *used_rho2 = NULL;

    long   temp_date, end_date;
    double cal_alpha, cal_rho, cal_lameps;

    char Diag_tenor[20] = "DIAG", Short_tenor[20] = "12M";

    char *prim_ten, *prim_ref, *prim_freq, *prim_basis, *sec_ten, *sec_ref, *sec_freq, *sec_basis;

    cpd_diag_calib_param *prim_param = NULL, *sec_param = NULL;

    diag_calib_lm_params* lm_params = NULL;

    LGMSV_CalibParams* calib_params = NULL;

    int    i;
    long   today;
    double step_calib;
    int    nb_new_dates;
    int    skip_dates, nb_cal_dates;
    int    one_gearing;

    Err err = NULL;

    /* get today */
    today = mkt->today;

    /* Initialisation */
    prim_param   = (cpd_diag_calib_param*)calloc(1, sizeof(cpd_diag_calib_param));
    sec_param    = (cpd_diag_calib_param*)calloc(1, sizeof(cpd_diag_calib_param));
    lm_params    = (diag_calib_lm_params*)calloc(1, sizeof(diag_calib_lm_params));
    calib_params = (LGMSV_CalibParams*)calloc(1, sizeof(LGMSV_CalibParams));

    used_nsmilepar    = nsmilepar;
    used_smilepartime = (double*)calloc(used_nsmilepar, sizeof(double));
    used_alpha        = (double*)calloc(used_nsmilepar, sizeof(double));
    used_rho          = (double*)calloc(used_nsmilepar, sizeof(double));
    used_ldaeps       = (double*)calloc(used_nsmilepar, sizeof(double));
    used_rho2         = (double*)calloc(used_nsmilepar, sizeof(double));

    if (!prim_param || !sec_param || !lm_params || !calib_params || !used_smilepartime ||
        !used_alpha || !used_rho || !used_ldaeps || !used_rho2)
    {
        err = "Memory allocation faillure in cts_calib_und";
        goto FREE_RETURN;
    }

    /* Recopy */
    if (smilepartime)
        memcpy(used_smilepartime, smilepartime, nsmilepar * sizeof(double));
    if (alpha)
        memcpy(used_alpha, alpha, nsmilepar * sizeof(double));
    if (rho)
        memcpy(used_rho, rho, nsmilepar * sizeof(double));
    if (ldaeps)
        memcpy(used_ldaeps, ldaeps, nsmilepar * sizeof(double));
    if (rho2)
        memcpy(used_rho2, rho2, nsmilepar * sizeof(double));

    und->model.iNbPWTime  = 0;
    und->model.dPWTime    = NULL;
    und->model.dSigma     = NULL;
    und->model.dAlpha     = NULL;
    und->model.dLambdaEps = NULL;
    und->model.dLvlEps    = NULL;
    und->model.dRho       = NULL;
    und->model.dRho2      = NULL;

    und->mkt = mkt;
    strcpy(und->name, "CALIB");

    cpd_init_calib_inst_data(inst_data);

    if (cts->num_calls == 0 || (cts->num_calls == 1 && cts->call[0].ex_date <= today + eod_flag))
    {
        /*	If no call dates, exit */
        cts->num_calls = 0;
        goto FREE_RETURN;
    }

    /* Eliminate dates too close to today */

    skip_dates = 0;
    while (skip_dates < cts->num_calls && cts->call[skip_dates].ex_time < MIN_CALIB_TIME)
    {
        skip_dates++;
    }

    if (skip_dates == cts->num_calls)
    {
        /* all exercise dates are too close: we keep the last one */
        skip_dates = cts->num_calls - 1;
    }

    nb_cal_dates = cts->num_calls - skip_dates;

    /* Memory Allocation */

    cal_dates      = (int*)calloc(nb_cal_dates, sizeof(int));
    ex_dates       = (long*)calloc(nb_cal_dates, sizeof(long));
    prim_strikes   = (double*)calloc(nb_cal_dates, sizeof(double));
    prim_strikesS1 = (double*)calloc(nb_cal_dates, sizeof(double));
    prim_strikesS2 = (double*)calloc(nb_cal_dates, sizeof(double));
    sec_strikes    = (double*)calloc(nb_cal_dates, sizeof(double));
    sec_strikesS1  = (double*)calloc(nb_cal_dates, sizeof(double));
    sec_strikesS2  = (double*)calloc(nb_cal_dates, sizeof(double));
    prim_tenor     = (char**)calloc(nb_cal_dates, sizeof(char*));
    sec_tenor      = (char**)calloc(nb_cal_dates, sizeof(char*));

    if (!cal_dates || !ex_dates || !prim_strikes || !sec_strikes || !prim_tenor || !sec_tenor ||
        !prim_strikesS1 || !prim_strikesS2 || !sec_strikesS1 || !sec_strikesS2)
    {
        err = "Memory allocation faillure in cts_calib_und";
        goto FREE_RETURN;
    }

    /* Fill calibration informations */

    /* Default strategy */
    if (calib_strategy == 0)
    {
        one_gearing = 1;

        for (i = cts->call[0].exo_idx; i < cts->exo_leg->num_cpn; i++)
        {
            if (cts->exo_leg->cpn[i].type == 0 ||
                fabs(cts->exo_leg->cpn[i].pay_gearing - 1.0) > 1.0E-08)
            {
                one_gearing = 0;
                break;
            }
        }

        if (one_gearing)
        {
            calib_strategy = 2;
        }
        else
        {
            calib_strategy = 1;
        }
    }

    switch (calib_strategy)
    {
    case 1:
    {
        /* Swaption / Cap */
        prim_ten   = &(Diag_tenor[0]);
        prim_freq  = mkt->swap_freq;
        prim_basis = mkt->swap_basis;
        prim_ref   = mkt->ref;

        sec_ten   = short_tenor;
        sec_freq  = mkt->lib_freq;
        sec_basis = mkt->lib_basis;
        sec_ref   = short_refrate;

        if (smile_calib_months < 0)
        {
            smile_calib_months = 0;
        }

        break;
    }

    case 2:
    {
        /* Cap / Swaption */
        prim_ten   = short_tenor;
        prim_freq  = mkt->lib_freq;
        prim_basis = mkt->lib_basis;
        prim_ref   = short_refrate;

        sec_ten   = &(Diag_tenor[0]);
        sec_freq  = mkt->swap_freq;
        sec_basis = mkt->swap_basis;
        sec_ref   = mkt->ref;

        err = add_tenor(mkt->today, short_tenor, NO_BUSDAY_CONVENTION, &temp_date);

        if (err)
            goto FREE_RETURN;

        if (smile_calib_months < 0)
        {
            smile_calib_months = (int)((temp_date - mkt->today) * YEARS_IN_DAY * 12.0 + 0.5);
        }

        break;
    }

    case 3:
    {
        /* User + Default 0 */
        /* Primary */
        if (!cal_tenor || strlen(cal_tenor) == 0)
        {
            prim_ten = &(Diag_tenor[0]);
        }
        else
        {
            prim_ten = cal_tenor;
        }

        if (!cal_ref || strlen(cal_ref) == 0)
        {
            prim_ref = mkt->ref;
        }
        else
        {
            prim_ref = cal_ref;
        }

        if (!cal_freq || strlen(cal_freq) == 0)
        {
            prim_freq = mkt->swap_freq;
        }
        else
        {
            prim_freq = cal_freq;
        }

        if (!cal_basis || strlen(cal_basis) == 0)
        {
            prim_basis = mkt->swap_basis;
        }
        else
        {
            prim_basis = cal_basis;
        }

        /* Secondary */
        if (!short_tenor || strlen(short_tenor) == 0)
        {
            sec_ten = &(Short_tenor[0]);
        }
        else
        {
            sec_ten = short_tenor;
        }

        if (!short_refrate || strlen(short_refrate) == 0)
        {
            sec_ref = mkt->ref;
        }
        else
        {
            sec_ref = short_refrate;
        }

        if (!short_freq || strlen(short_freq) == 0)
        {
            sec_freq = mkt->lib_freq;
        }
        else
        {
            sec_freq = short_freq;
        }

        if (!short_basis || strlen(short_basis) == 0)
        {
            sec_basis = mkt->lib_basis;
        }
        else
        {
            sec_basis = short_basis;
        }

        if (smile_calib_months < 0)
        {
            /* Autocalibration to SWAPTION smile */
            smile_calib_months = 0;
        }

        break;
    }

    default:
    {
        err = "Calibration strategy of CTS not recognised";
        goto FREE_RETURN;
    }
    }

    for (i = 0; i < nb_cal_dates; i++)
    {
        cal_dates[i]      = 1;
        ex_dates[i]       = cts->call[i + skip_dates].ex_date;
        prim_strikes[i]   = 0.0;
        prim_strikesS1[i] = 1.0;
        prim_strikesS2[i] = -1.0;
        sec_strikes[i]    = 0.0;
        sec_strikesS1[i]  = 1.0;
        sec_strikesS2[i]  = -1.0;
        prim_tenor[i]     = prim_ten;
        sec_tenor[i]      = sec_ten;
    }

    if (long_strike_flag == 1 && force_atm == 0)
    {
        for (i = 0; i < nb_cal_dates; i++)
        {
            prim_strikes[i] = cts->call[i + skip_dates].implied_strike;
        }
    }

    if (short_strike_flag == 1 && force_atm == 0)
    {
        for (i = 0; i < nb_cal_dates; i++)
        {
            sec_strikes[i] = cts->exo_leg->cpn[cts->call[i + skip_dates].exo_idx].caplet_strike;
        }
    }

    if (cts->exo_leg->num_cpn > 0)
    {
        end_date = cts->exo_leg->cpn[cts->exo_leg->num_cpn - 1].cpn_pay_date;
    }
    else
    {
        end_date = 0;
    }

    if (end_date < mkt->today + 370)
    {
        end_date = mkt->today + 370;
    }

    /* Fill calibration parameters */
    cpd_calib_set_default_param(prim_param);
    cpd_calib_set_default_param(sec_param);

    prim_param->precision   = prec;
    prim_param->vol_shift   = vol_shift_long;
    prim_param->shift_type  = vol_shift_type_long;
    prim_param->vol_type    = vol_type_long;
    prim_param->strike_type = 1;
    prim_param->max_std     = max_std_long;
    prim_param->min_time    = min_time;
    prim_param->skip_last   = skip_last;
    prim_param->min_fact    = min_fact;
    prim_param->max_fact    = max_fact;
    prim_param->use_jumps   = use_jumps;
    prim_param->keep_first  = keep_first;

    sec_param->precision    = prec;
    sec_param->vol_shift    = vol_shift_short;
    sec_param->shift_type   = vol_shift_type_short;
    sec_param->vol_type     = vol_type_short;
    sec_param->lambda_shift = lambda_shift;
    sec_param->strike_type  = 1;
    sec_param->max_std      = max_std_short;
    sec_param->min_time     = min_time;
    sec_param->skip_last    = skip_last;
    sec_param->min_fact     = min_fact;
    sec_param->max_fact     = max_fact;
    sec_param->use_jumps    = 1;
    sec_param->keep_first   = keep_first;

    switch (short_strike_flag)
    {
    case 0:
        sec_param->strike_type = 0;
        break;
    case 1:
        sec_param->strike_type = 2;
        break;
    case 2:
        sec_param->strike_type = 4;
        break;
    }

    sec_param->max_std    = max_std_short;
    sec_param->min_time   = min_time;
    sec_param->skip_last  = skip_last;
    sec_param->keep_first = keep_first;

    lm_params->nb_iter      = maxiter;
    lm_params->use_moment   = 1;
    lm_params->vega_weight  = 0;
    lm_params->freq_short   = 1;
    lm_params->shift_freq   = 0;
    lm_params->nb_moment    = 1;
    lm_params->break_moment = 0;
    lm_params->precision    = prec;

    /* fill the calib params */
    LGMSV_Copy_CalibParams(lgmsv_calib_params, calib_params);

    /* First calibrate the smile if needed */
    if (!fix_smile && calib_params->use_sabr_calib)
    {
        err = cts_calib_const_alpha_and_rho(
            mkt,
            cts,
            smile_calib_months,
            LGMSV_SABRBETA,
            min_time,
            &cal_alpha,
            &cal_rho,
            &cal_lameps,
            save_inst_data,
            inst_data);

        if (err)
            goto FREE_RETURN;

        for (i = 0; i < nsmilepar; i++)
        {
            used_alpha[i]  = cal_alpha + calib_params->alpha_sv_shift;
            used_ldaeps[i] = cal_lameps + calib_params->lam_sv_shift;
            used_rho[i]    = cal_rho + calib_params->rho_sv_shift;
        }

        calib_params->alpha_sv_shift = 0.0;
        calib_params->lam_sv_shift   = 0.0;
        calib_params->rho_sv_shift   = 0.0;
    }
    else
    {
        for (i = 0; i < nsmilepar; i++)
        {
            used_alpha[i]  = alpha[i];
            used_ldaeps[i] = ldaeps[i];
            used_rho[i]    = rho[i];
            if (rho2)
                used_rho2[i] = rho2[i];
        }
    }

    /* Check Limits */
    for (i = 0; i < nsmilepar; i++)
    {
        used_alpha[i]  = DMAX(used_alpha[i], 0.01);
        used_rho[i]    = DMIN(DMAX(used_rho[i], -0.99), 0.99);
        used_ldaeps[i] = DMAX(used_ldaeps[i], 0.01);

        if (rho2)
            used_rho2[i] = DMIN(DMAX(used_rho2[i], -0.99), 0.99);
    }

    calib_params->fix_lambda = fix_lambda;

    if (calib_params->use_sabr_calib || fix_smile)
    {
        calib_params->calib_alpha  = 0;
        calib_params->calib_lameps = 0;
        calib_params->calib_rho    = 0;
        calib_params->calib_rho2   = 0;
    }
    else
    {
        calib_params->calib_alpha  = 1;
        calib_params->calib_lameps = 1;
        calib_params->calib_rho    = 1;
        calib_params->calib_rho2   = 0;
    }

    /* Launch the calibration */
    err = cpd_calib_diagonal_LGMSV_new_dlm(
        mkt->yc,
        mkt->vc,
        mkt->get_cash_vol,
        mkt->ref,

        prim_freq,
        prim_basis,
        prim_ref,
        nb_cal_dates,
        ex_dates,
        cal_dates,
        prim_tenor,
        end_date,
        prim_strikes,
        prim_strikesS1,
        prim_strikesS2,
        prim_param,

        sec_freq,
        sec_basis,
        sec_ref,
        nb_cal_dates,
        ex_dates,
        cal_dates,
        sec_tenor,
        end_date,
        sec_strikes,
        sec_strikesS1,
        sec_strikesS2,
        NULL,
        sec_param,

        calib_params,

        nb_factor,
        &lambda,
        used_nsmilepar,
        used_smilepartime,
        used_alpha,
        used_ldaeps,
        used_rho,
        tstar,
        lgm_alpha,
        lgm_gamma,
        lgm_rho,
        used_rho2,

        NumerParams,

        &(und->model.iNbPWTime),
        &(und->model.dPWTime),
        &(und->model.dSigma),
        &(und->model.dAlpha),
        &(und->model.dLambdaEps),
        &(und->model.dRho),
        &(und->model.dRho2),

        lm_params,
        inst_data);

    if (err)
        goto FREE_RETURN;

    if (fabs(lambda) < 1.0E-10)
    {
        lambda = 0.0001;
    }

    und->model.dLambdaX   = lambda;
    und->model.dTau       = 1.0 / lambda;
    und->model.dTStar     = tstar;
    und->model.dInitTStar = tstar;
    und->model.lToday     = mkt->today;
    und->model.lTStarDate = (long)(und->model.lToday + und->model.dTStar * DAYS_IN_YEAR + 1.0E-08);
    und->model.iOne2F     = nb_factor;

    /* mutliply the smile parameters by 2 */
    for (i = 0; i < und->model.iNbPWTime; i++)
    {
        und->model.dAlpha[i] *= 2.0;
        und->model.dLambdaEps[i] *= 2.0;
    }

    /* Save the Lambda Eps for Reserve */
    und->model.dInitLambdaEps = und->model.dLambdaEps[0];

    /*	add discretisation dates till the end of the deal */

    if (und->model.iNbPWTime > 1)
    {
        step_calib = und->model.dPWTime[1] - und->model.dPWTime[0];
    }
    else
    {
        step_calib = 1.0;
    }

    nb_new_dates = (int) ((cts->exo_leg->cpn[cts->exo_leg->num_cpn-1].cpn_pay_time - und->model.dPWTime[und->model.iNbPWTime - 1]) / step_calib) + 1;

    und->model.iNbPWTime += nb_new_dates;

    und->model.dPWTime =
        (double*)realloc(und->model.dPWTime, und->model.iNbPWTime * sizeof(double));
    und->model.dSigma = (double*)realloc(und->model.dSigma, und->model.iNbPWTime * sizeof(double));
    und->model.dAlpha = (double*)realloc(und->model.dAlpha, und->model.iNbPWTime * sizeof(double));
    und->model.dRho   = (double*)realloc(und->model.dRho, und->model.iNbPWTime * sizeof(double));
    und->model.dLambdaEps =
        (double*)realloc(und->model.dLambdaEps, und->model.iNbPWTime * sizeof(double));
    und->model.dLvlEps = (double*)calloc(und->model.iNbPWTime, sizeof(double));

    if (!und->model.dPWTime || !und->model.dSigma || !und->model.dAlpha || !und->model.dRho ||
        !und->model.dLambdaEps || !und->model.dLvlEps)
    {
        err = "Memory allocation error in cts_calib_und";
        goto FREE_RETURN;
    }

    if (und->model.iOne2F == 2)
    {
        und->model.dRho2 =
            (double*)realloc(und->model.dRho2, und->model.iNbPWTime * sizeof(double));
        und->model.dLGMAlpha = (double*)calloc(und->model.iNbPWTime, sizeof(double));
        und->model.dLGMRho   = (double*)calloc(und->model.iNbPWTime, sizeof(double));

        if (!und->model.dRho2 || !und->model.dLGMAlpha || !und->model.dLGMRho)
        {
            err = "Memory allocation error in cts_calib_und";
            goto FREE_RETURN;
        }

        und->model.dInitLGMAlpha = lgm_alpha;
        und->model.dLGMGamma     = lgm_gamma;
        und->model.dInitLGMRho   = lgm_rho;

        und->model.dLambdaX2 = und->model.dLambdaX + und->model.dLGMGamma;

        if (fabs(und->model.dLambdaX2) < 1.0E-08)
        {
            und->model.dLambdaX2 = 0.0001;
        }

        und->model.dTau2 = 1.0 / und->model.dLambdaX2;
    }

    for (i = 0; i < nb_new_dates; i++)
    {
        und->model.dPWTime[und->model.iNbPWTime - nb_new_dates + i] =
            und->model.dPWTime[und->model.iNbPWTime - nb_new_dates + i - 1] + step_calib;
        und->model.dSigma[und->model.iNbPWTime - nb_new_dates + i] =
            und->model.dSigma[und->model.iNbPWTime - nb_new_dates + i - 1];
        und->model.dAlpha[und->model.iNbPWTime - nb_new_dates + i] =
            und->model.dAlpha[und->model.iNbPWTime - nb_new_dates + i - 1];
        und->model.dRho[und->model.iNbPWTime - nb_new_dates + i] =
            und->model.dRho[und->model.iNbPWTime - nb_new_dates + i - 1];
        und->model.dLambdaEps[und->model.iNbPWTime - nb_new_dates + i] =
            und->model.dLambdaEps[und->model.iNbPWTime - nb_new_dates + i - 1];

        if (und->model.iOne2F == 2)
        {
            und->model.dRho2[und->model.iNbPWTime - nb_new_dates + i] =
                und->model.dRho2[und->model.iNbPWTime - nb_new_dates + i - 1];
        }
    }

    memcpy(und->model.dLvlEps, und->model.dLambdaEps, und->model.iNbPWTime * sizeof(double));

    ConvertTS_LGM_to_LGMSV(
        und->model.iNbPWTime,
        und->model.dPWTime,
        und->model.dSigma,
        und->model.dLambdaX,
        und->model.dTStar,
        und->model.iOne2F,
        und->model.dInitLGMAlpha,
        und->model.dLGMGamma,
        und->model.dInitLGMRho,
        und->model.dLGMAlpha,
        und->model.dLGMRho);

FREE_RETURN:

    if (prim_param)
        free(prim_param);
    if (sec_param)
        free(sec_param);
    if (lm_params)
        free(lm_params);
    if (calib_params)
        free(calib_params);

    if (cal_dates)
        free(cal_dates);
    if (ex_dates)
        free(ex_dates);
    if (prim_strikes)
        free(prim_strikes);
    if (prim_strikesS1)
        free(prim_strikesS1);
    if (prim_strikesS2)
        free(prim_strikesS2);

    if (sec_strikes)
        free(sec_strikes);
    if (sec_strikesS1)
        free(sec_strikesS1);
    if (sec_strikesS2)
        free(sec_strikesS2);

    if (prim_tenor)
        free(prim_tenor);
    if (sec_tenor)
        free(sec_tenor);

    if (used_smilepartime)
        free(used_smilepartime);
    if (used_alpha)
        free(used_alpha);
    if (used_rho)
        free(used_rho);
    if (used_ldaeps)
        free(used_ldaeps);
    if (used_rho2)
        free(used_rho2);

    if (err)
        cts_free_und(und);

    return err;
}

/*	Copy */
Err cts_copy_und(CTS_UND src, CTS_UND dest)
{
    Err err = NULL;

    if (dest && src)
    {
        err = cts_copy_mkt(src->mkt, dest->mkt);
        if (err)
            return err;

        strcpy(dest->name, src->name);

        err = init_LGMSV_model(
            &(dest->model),
            src->model.lToday,
            src->model.iOne2F,
            src->model.iNbPWTime,
            src->model.dLambdaX,
            src->model.dPWTime,
            src->model.dSigma,
            src->model.dAlpha,
            src->model.dLambdaEps,
            src->model.dLvlEps,
            src->model.dRho,
            src->model.dTStar,
            src->model.dInitLGMAlpha,
            src->model.dLGMGamma,
            src->model.dInitLGMRho,
            src->model.dRho2);

        dest->model.dInitTStar = src->model.dInitTStar;
    }
    else
    {
        err = "unable to copy the underlying model";
    }

    return err;
}

/*	Free */
void cts_free_und(CTS_UND und)
{
    if (und)
    {
        free_LGMSV_model(&(und->model));
        cts_free_mkt(und->mkt);
    }

    und = NULL;
}

/*	Calculate forward and market IVs */

/*	Initialise the structure */
void cts_init_cts_iv(CTS_IV iv_infos)
{
    iv_infos->call_date    = NULL;
    iv_infos->extra_fee    = NULL;
    iv_infos->market_fwdiv = NULL;
    iv_infos->model_fwdiv  = NULL;

    iv_infos->one_time_call     = NULL;
    iv_infos->one_time_coupon   = NULL;
    iv_infos->one_time_notional = NULL;
    iv_infos->one_time_log_vol  = NULL;
    iv_infos->one_time_norm_vol = NULL;
    iv_infos->one_time_strike   = NULL;

    iv_infos->ncall        = 0;
    iv_infos->has_one_time = 0;
}

/*	Calculate forward and market IVs */
void cts_free_cts_iv(CTS_IV iv_infos)
{
    if (iv_infos)
    {
        if (iv_infos->call_date)
            free(iv_infos->call_date);
        iv_infos->call_date = NULL;

        if (iv_infos->market_fwdiv)
            free(iv_infos->market_fwdiv);
        iv_infos->market_fwdiv = NULL;

        if (iv_infos->model_fwdiv)
            free(iv_infos->model_fwdiv);
        iv_infos->model_fwdiv = NULL;

        if (iv_infos->extra_fee)
            free(iv_infos->extra_fee);
        iv_infos->extra_fee = NULL;

        if (iv_infos->has_one_time)
        {
            if (iv_infos->one_time_call)
                free(iv_infos->one_time_call);
            iv_infos->one_time_call = NULL;

            if (iv_infos->one_time_strike)
                free(iv_infos->one_time_strike);
            iv_infos->one_time_strike = NULL;

            if (iv_infos->one_time_norm_vol)
                free(iv_infos->one_time_norm_vol);
            iv_infos->one_time_norm_vol = NULL;

            if (iv_infos->one_time_log_vol)
                free(iv_infos->one_time_log_vol);
            iv_infos->one_time_log_vol = NULL;

            if (iv_infos->one_time_coupon)
                free(iv_infos->one_time_coupon);
            iv_infos->one_time_coupon = NULL;

            if (iv_infos->one_time_notional)
                free(iv_infos->one_time_notional);
            iv_infos->one_time_notional = NULL;

            iv_infos->has_one_time = 0;
        }

        iv_infos->ncall = 0;
    }
}

void cts_free_extra_infos(CTS_EXTRA_INFOS extra_infos)
{
    if (extra_infos)
    {
        if (extra_infos->extra_infos && extra_infos->nrow > 0 && extra_infos->ncol > 0)
        {
            free_dmatrix(
                extra_infos->extra_infos, 0, extra_infos->nrow - 1, 0, extra_infos->ncol - 1);
            extra_infos->extra_infos = NULL;
            extra_infos->nrow        = 0;
            extra_infos->ncol        = 0;
        }

        extra_infos = NULL;
    }
}

void cts_switch_coef_libor_option(LGMSV_HESTONINST cHestonInst)
{
    double temp;

    temp                       = cHestonInst->dCoefMeanRev;
    cHestonInst->dCoefMeanRev  = cHestonInst->dCoefMeanRev2;
    cHestonInst->dCoefMeanRev2 = temp;

    temp                   = cHestonInst->dCoefCMS;
    cHestonInst->dCoefCMS  = cHestonInst->dCoefCMS2;
    cHestonInst->dCoefCMS2 = temp;

    temp                      = cHestonInst->dNewCoefCMS;
    cHestonInst->dNewCoefCMS  = cHestonInst->dNewCoefCMS2;
    cHestonInst->dNewCoefCMS2 = temp;

    /* For 2Factor */
    temp                          = cHestonInst->dCoefMeanRev_2F;
    cHestonInst->dCoefMeanRev_2F  = cHestonInst->dCoefMeanRev2_2F;
    cHestonInst->dCoefMeanRev2_2F = temp;

    temp                      = cHestonInst->dCoefCMS_2F;
    cHestonInst->dCoefCMS_2F  = cHestonInst->dCoefCMS2_2F;
    cHestonInst->dCoefCMS2_2F = temp;

    temp                         = cHestonInst->dNewCoefCMS_2F;
    cHestonInst->dNewCoefCMS_2F  = cHestonInst->dNewCoefCMS2_2F;
    cHestonInst->dNewCoefCMS2_2F = temp;

    temp                               = cHestonInst->dNewCoefCMS_cross_2F;
    cHestonInst->dNewCoefCMS_cross_2F  = cHestonInst->dNewCoefCMS2_cross_2F;
    cHestonInst->dNewCoefCMS2_cross_2F = temp;
}

/*	Calculate forward and market IVs */
Err cts_calc_model_fwd_iv(
    CTS_UND und,
    CTS     cts,
    int     calc_fwd_iv,
    int     adj_fee,
    /*	Numerical CF params */
    LGMSV_NUMERPARAMS NumerParams)
{
    int     i, j, k;
    int     first_coupon;
    double  df_pay, dLiborFwd;
    CTS_CPN exo_cpn;
    CTS_FIX cts_fix;
    CTS_MKT mkt;

    LGMSV_PricingConst* sPricingConst = NULL;
    LGMSV_NumerInst*    sNumerParams  = NULL;

    CalibInstrumentDLM* cCalibInst  = NULL;
    LGMSV_HestonInst*   cHestonInst = NULL;

    CalibCpnScheduleDLM *cCpnSchedule = NULL, *cCpnSchedule2;

    double InstPrices_[10], InstPrices2_[10], InstPrices_temp_[10], InstStrikes_[10], InstStds_[10],
        *InstPrices = &(InstPrices_[0]), *InstPrices2 = &(InstPrices2_[0]),
        *InstPrices_temp = &(InstPrices_temp_[0]), *InstStrikes = &(InstStrikes_[0]),
        *InstStds = &(InstStds_[0]);

    long   today;
    double fixing_iv, fixing_mrg_iv, coupon_iv, coupon_mrg_iv;
    int    first_strike, last_strike, nb_strike, do_min_strike, do_max_strike;
    double min_vol, max_vol, dLiborRefFwd;

    Err err = NULL;

    /* Get today */
    mkt   = und->mkt;
    today = mkt->today;

    /* Memory allocations */
    sPricingConst = (LGMSV_PricingConst*)calloc(1, sizeof(LGMSV_PricingConst));
    sNumerParams  = (LGMSV_NumerInst*)calloc(1, sizeof(LGMSV_NumerInst));
    cCalibInst    = (CalibInstrumentDLM*)calloc(1, sizeof(CalibInstrumentDLM));
    cHestonInst   = (LGMSV_HestonInst*)calloc(1, sizeof(LGMSV_HestonInst));
    cCpnSchedule  = (CalibCpnScheduleDLM*)calloc(1, sizeof(CalibCpnScheduleDLM));

    if (!sPricingConst || !sNumerParams || !cCalibInst || !cHestonInst || !cCpnSchedule)
    {
        err = "Memory allocation faillure in cts_calc_model_fwd_iv";
        goto FREE_RETURN;
    }

    if (calc_fwd_iv)
    {
        first_coupon = cts->call[0].exo_idx;

        /* Initialise the structures */
        err = Initialise_PricingConst(&(und->model), NumerParams, sPricingConst);

        cCpnSchedule->iNCpn  = 2;
        cCpnSchedule->lToday = today;

        cCalibInst->iIsCMS    = 1;
        cCalibInst->iStartCpn = 0;
        cCalibInst->dStrike   = InstStrikes;
        cCalibInst->dNbStd    = InstStds;

        cHestonInst->dMinStd = NumerParams->dMinStd;
        cHestonInst->dMaxStd = NumerParams->dMaxStd;

        /* calculate the exotic IV for each coupon */

        for (i = first_coupon; i < cts->exo_leg->num_cpn; i++)
        {
            exo_cpn              = &(cts->exo_leg->cpn[i]);
            cCalibInst->lPayDate = exo_cpn->cpn_pay_date;
            cCalibInst->dPayTime = exo_cpn->cpn_pay_time;
            df_pay               = swp_f_df(today, cCalibInst->lPayDate, mkt->yc);

            coupon_iv     = 0.0;
            coupon_mrg_iv = 0.0;

            if (exo_cpn->type == 0 || exo_cpn->type == 3 || exo_cpn->type == 4)
            {
                cCalibInst->iIsLiborOption = 0;
                cCalibInst->dLiborMargin   = 0.0;
            }
            else
            {
                /* Fill Libor informations */
                cCalibInst->iNbCoupon = 2;
                cCalibInst->iEndCpn   = 1;

                cCalibInst->iIsLiborOption  = 1;
                cCalibInst->lLiborFixDate   = exo_cpn->pay_fix_date;
                cCalibInst->dLiborFixTime   = exo_cpn->pay_fix_time;
                cCalibInst->lLiborStartDate = exo_cpn->pay_start_date;
                cCalibInst->dLiborStartTime = exo_cpn->pay_start_time;
                cCalibInst->lLiborEndDate   = exo_cpn->pay_end_date;
                cCalibInst->dLiborEndTime   = exo_cpn->pay_end_time;
                cCalibInst->dLiborCvg       = exo_cpn->pay_cvg;
                cCalibInst->dLiborSpread    = exo_cpn->pay_fwd_spread;
                cCalibInst->dLiborMargin    = exo_cpn->cpn_coupon;

                /* Calculate its expectation */

                /* Setup the Schedule */
                cCpnSchedule->lCpnDate[0] = exo_cpn->pay_start_date;
                cCpnSchedule->lCpnDate[1] = exo_cpn->pay_end_date;
                cCpnSchedule->dCpnTime[0] = (cCpnSchedule->lCpnDate[0] - today) * YEARS_IN_DAY;
                cCpnSchedule->dCpnTime[1] = (cCpnSchedule->lCpnDate[1] - today) * YEARS_IN_DAY;
                cCpnSchedule->dCpnCvg[1]  = exo_cpn->pay_cvg;
                cCpnSchedule->dCpnDf[0]   = swp_f_df(today, cCpnSchedule->lCpnDate[0], mkt->yc);
                cCpnSchedule->dCpnDf[1]   = swp_f_df(today, cCpnSchedule->lCpnDate[1], mkt->yc);

                /* Setup the instrument */
                cCalibInst->lExeDate = exo_cpn->pay_fix_date;
                cCalibInst->dExeTime = exo_cpn->pay_fix_time;
                cCalibInst->dLevel   = cCpnSchedule->dCpnCvg[1] * cCpnSchedule->dCpnDf[1];
                cCalibInst->dSumDf   = cCpnSchedule->dCpnDf[1];
                cCalibInst->dFwdCash =
                    (cCpnSchedule->dCpnDf[0] - cCpnSchedule->dCpnDf[1]) / cCalibInst->dLevel;
                cCalibInst->dSpread = exo_cpn->pay_fwd_spread;

                cCalibInst->iNbStrike = 0;

                /* Calculate the constant */
                err = Initialise_NumerInst(&(und->model), cCalibInst, cHestonInst, sNumerParams);

                if (err)
                {
                    goto FREE_RETURN;
                }

                /* For the forward calculation we remove the flag*/
                cCalibInst->iIsLiborOption = 0;

                err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
                    cCalibInst, cCpnSchedule, &(und->model), cHestonInst);

                if (err)
                {
                    goto FREE_RETURN;
                }

                err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(
                    cCalibInst, cCpnSchedule, &(und->model), cHestonInst);

                if (err)
                {
                    goto FREE_RETURN;
                }

                Update_NumerInst(cCalibInst, cHestonInst, sNumerParams);

                /* Call the Fwd Pricer */

                LGMSVFwdClosedFormApprox_struct(
                    &(und->model),
                    cCalibInst,
                    cHestonInst,
                    sPricingConst,
                    sNumerParams,
                    &dLiborFwd);

                cCalibInst->iIsLiborOption = 1;
            }

            /* Evaluate Exotic coupon */

            for (j = 0; j < exo_cpn->used_nfix; j++)
            {
                fixing_iv     = 0.0;
                fixing_mrg_iv = 0.0;

                cts_fix = &(exo_cpn->fix[j]);

                /* Setup the schedule */
                cCpnSchedule2 = &(cts_fix->schedule);

                /* Setup the instrument */
                cCalibInst->iNbCoupon = cCpnSchedule2->iNCpn;
                cCalibInst->iEndCpn   = cCpnSchedule2->iNCpn - 1;

                cCalibInst->lExeDate = cts_fix->ref_fix_date;
                cCalibInst->dExeTime = cts_fix->ref_fix_time;
                cCalibInst->dLevel   = cts_fix->ref_level;
                cCalibInst->dSumDf   = cts_fix->ref_sum_df;
                cCalibInst->dFwdCash = cts_fix->ref_swp_cash;
                cCalibInst->dSpread  = cts_fix->ref_fwd_spread;

                cCalibInst->iNbStrike = 10;

                /* Calculate the constant */
                err = Initialise_NumerInst(&(und->model), cCalibInst, cHestonInst, sNumerParams);

                if (err)
                {
                    goto FREE_RETURN;
                }

                err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
                    cCalibInst, cCpnSchedule2, &(und->model), cHestonInst);

                if (err)
                {
                    goto FREE_RETURN;
                }

                err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(
                    cCalibInst, cCpnSchedule2, &(und->model), cHestonInst);

                if (err)
                {
                    goto FREE_RETURN;
                }

                /* calculate the needed strikes */
                first_strike = 0;
                nb_strike    = 0;

                while (first_strike < exo_cpn->nstr &&
                       exo_cpn->str[first_strike] < cHestonInst->dMinStrike)
                {
                    first_strike++;
                }

                if (first_strike > 0 && first_strike < exo_cpn->nstr)
                {
                    do_min_strike                  = 1;
                    cCalibInst->dStrike[nb_strike] = cHestonInst->dMinStrike;
                    nb_strike++;
                }
                else
                {
                    do_min_strike = 0;
                }

                last_strike = first_strike;

                while (last_strike < exo_cpn->nstr &&
                       exo_cpn->str[last_strike] < cHestonInst->dMaxStrike)
                {
                    cCalibInst->dStrike[nb_strike] = exo_cpn->str[last_strike];
                    last_strike++;
                    nb_strike++;
                }

                if (last_strike < exo_cpn->nstr)
                {
                    do_max_strike                  = 1;
                    cCalibInst->dStrike[nb_strike] = cHestonInst->dMaxStrike;
                    nb_strike++;
                }
                else
                {
                    do_max_strike = 0;
                }

                cCalibInst->iNbStrike = nb_strike;

                InstPrices[exo_cpn->nstr] = 0.0;

                if (fabs(exo_cpn->beta) > 1.0E-08 || do_min_strike || do_max_strike)
                {
                    /* we add another strike */
                    cCalibInst->iNbStrike += 1;
                    cCalibInst->dStrike[nb_strike] = 0.0;
                }

                sNumerParams->iNbStrike = cCalibInst->iNbStrike;

                Update_NumerInst(cCalibInst, cHestonInst, sNumerParams);

                /* adjustment for Libor Option */
                if (cCalibInst->iIsLiborOption)
                {
                    sNumerParams->dNewSwitchTime = exo_cpn->pay_fix_time;
                    sNumerParams->iSwitchIndex   = Get_Index(
                        sNumerParams->dNewSwitchTime, und->model.dPWTime, und->model.iNbPWTime);

                    if (cCalibInst->dExeTime < sNumerParams->dNewSwitchTime)
                    {
                        /* no special switch */
                        sNumerParams->iDoSwitch      = 0;
                        sNumerParams->dNewSwitchTime = cCalibInst->dExeTime;
                        sNumerParams->iSwitchIndex   = sNumerParams->iEndIndex;
                    }
                    else
                    {
                        if (fabs(
                                sNumerParams->dNewSwitchTime -
                                und->model.dPWTime[sNumerParams->iSwitchIndex]) > 1.0E-08 ||
                            sNumerParams->iSwitchIndex == sNumerParams->iEndIndex)
                        {
                            /* we add the switch */
                            sNumerParams->iDoSwitch = 1;
                        }
                        else
                        {
                            sNumerParams->iDoSwitch = 0;
                        }
                    }

                    sNumerParams->iEndIndex += sNumerParams->iDoSwitch;
                }

                /* Call the pricer */
                if (cCalibInst->iNbStrike > 0)
                {
                    LGMSVClosedFormApprox_struct(
                        &(und->model),
                        cCalibInst,
                        cHestonInst,
                        sPricingConst,
                        sNumerParams,
                        InstPrices_temp);
                }

                /* Reprice the instruments */
                if (do_min_strike || do_max_strike)
                {
                    dLiborRefFwd = InstPrices_temp[cCalibInst->iNbStrike - 1];
                }

                if (do_min_strike)
                {
                    err = srt_f_optimpvol(
                        InstPrices_temp[0],
                        dLiborRefFwd + 100,
                        cCalibInst->dStrike[0] + 100,
                        cCalibInst->dExeTime,
                        df_pay,
                        SRT_CALL,
                        SRT_NORMAL,
                        &min_vol);

                    if (err)
                    {
                        min_vol = 1.0E-08;
                        err     = NULL;
                    }

                    for (k = 0; k < first_strike; k++)
                    {
                        InstPrices[k] = srt_f_optblknrm(
                            dLiborRefFwd,
                            exo_cpn->str[k],
                            min_vol,
                            cCalibInst->dExeTime,
                            df_pay,
                            SRT_CALL,
                            (SrtGreekType)SRT_PREMIUM);
                    }
                }

                for (k = first_strike; k < last_strike; k++)
                {
                    InstPrices[k] = InstPrices_temp[k - first_strike + do_min_strike];
                }

                if (do_max_strike)
                {
                    err = srt_f_optimpvol(
                        InstPrices_temp[nb_strike - 1],
                        dLiborRefFwd + 100,
                        cCalibInst->dStrike[nb_strike - 1] + 100,
                        cCalibInst->dExeTime,
                        df_pay,
                        SRT_CALL,
                        SRT_NORMAL,
                        &max_vol);

                    if (err)
                    {
                        max_vol = 1.0E-08;
                        err     = NULL;
                    }

                    for (k = last_strike; k < exo_cpn->nstr; k++)
                    {
                        InstPrices[k] = srt_f_optblknrm(
                            dLiborRefFwd,
                            exo_cpn->str[k],
                            max_vol,
                            cCalibInst->dExeTime,
                            df_pay,
                            SRT_CALL,
                            (SrtGreekType)SRT_PREMIUM);
                    }
                }

                if (cCalibInst->iIsLiborOption)
                {
                    /* roll back to the initial setting */

                    sNumerParams->iEndIndex -= sNumerParams->iDoSwitch;
                    cts_switch_coef_libor_option(cHestonInst);

                    /* calculate the non adjusted option as well */

                    sNumerParams->dNewSwitchTime = cCalibInst->dExeTime;
                    sNumerParams->iSwitchIndex   = sNumerParams->iEndIndex;
                    sNumerParams->iDoSwitch      = 0;

                    if (cCalibInst->iNbStrike > 0)
                    {
                        LGMSVClosedFormApprox_struct(
                            &(und->model),
                            cCalibInst,
                            cHestonInst,
                            sPricingConst,
                            sNumerParams,
                            InstPrices_temp);
                    }

                    if (do_min_strike || do_max_strike)
                    {
                        dLiborRefFwd = InstPrices_temp[cCalibInst->iNbStrike - 1];
                    }

                    if (do_min_strike)
                    {
                        err = srt_f_optimpvol(
                            InstPrices_temp[0],
                            dLiborRefFwd + 100,
                            cCalibInst->dStrike[0] + 100,
                            cCalibInst->dExeTime,
                            df_pay,
                            SRT_CALL,
                            SRT_NORMAL,
                            &min_vol);

                        if (err)
                        {
                            min_vol = 1.0E-08;
                            err     = NULL;
                        }

                        for (k = 0; k < first_strike; k++)
                        {
                            InstPrices2[k] = srt_f_optblknrm(
                                dLiborRefFwd,
                                exo_cpn->str[k],
                                min_vol,
                                cCalibInst->dExeTime,
                                df_pay,
                                SRT_CALL,
                                (SrtGreekType)SRT_PREMIUM);
                        }
                    }

                    for (k = first_strike; k < last_strike; k++)
                    {
                        InstPrices2[k] = InstPrices_temp[k - first_strike + do_min_strike];
                    }

                    if (do_max_strike)
                    {
                        err = srt_f_optimpvol(
                            InstPrices_temp[nb_strike - 1],
                            dLiborRefFwd + 100,
                            cCalibInst->dStrike[nb_strike - 1] + 100,
                            cCalibInst->dExeTime,
                            df_pay,
                            SRT_CALL,
                            SRT_NORMAL,
                            &max_vol);

                        if (err)
                        {
                            max_vol = 1.0E-08;
                            err     = NULL;
                        }

                        for (k = last_strike; k < exo_cpn->nstr; k++)
                        {
                            InstPrices2[k] = srt_f_optblknrm(
                                dLiborRefFwd,
                                exo_cpn->str[k],
                                max_vol,
                                cCalibInst->dExeTime,
                                df_pay,
                                SRT_CALL,
                                (SrtGreekType)SRT_PREMIUM);
                        }
                    }

                    /* update the price */
                    for (k = 0; k < cCalibInst->iNbStrike; k++)
                    {
                        InstPrices[k] = 1.0 / exo_cpn->pay_cvg *
                                            (InstPrices[k] * (dLiborFwd * exo_cpn->pay_cvg + 1.0) -
                                             InstPrices2[k]) +
                                        exo_cpn->cpn_coupon;
                    }

                    fixing_iv +=
                        exo_cpn->alpha * dLiborFwd + exo_cpn->beta * InstPrices[exo_cpn->nstr];
                    fixing_mrg_iv += exo_cpn->alpha + exo_cpn->beta * dLiborFwd;

                    for (k = 0; k < exo_cpn->nstr; k++)
                    {
                        fixing_iv += exo_cpn->nbopt[k] * InstPrices[k];
                        fixing_mrg_iv += exo_cpn->nbopt[k] * InstPrices2[k];
                    }
                }
                else
                {
                    if (fabs(exo_cpn->beta) > 1.0E-08)
                    {
                        fixing_iv +=
                            exo_cpn->alpha + exo_cpn->beta * InstPrices_temp[exo_cpn->nstr];
                    }
                    else
                    {
                        fixing_iv += exo_cpn->alpha;
                    }

                    for (k = 0; k < exo_cpn->nstr; k++)
                    {
                        fixing_iv += exo_cpn->nbopt[k] * InstPrices[k];
                    }
                }

                coupon_iv += fixing_iv * exo_cpn->used_fixweights[j];
                coupon_mrg_iv += fixing_mrg_iv * exo_cpn->used_fixweights[j];

                /* Free local Memory */
                Free_NumerInst(sNumerParams);
            }

            if (cCalibInst->iIsLiborOption)
            {
                coupon_iv *= exo_cpn->pay_gearing;
                coupon_iv += coupon_mrg_iv * exo_cpn->cpn_coupon;
            }
            else if (exo_cpn->type == 0)
            {
                coupon_iv *= exo_cpn->cpn_coupon;
            }

            coupon_iv *= exo_cpn->cpn_cvg * exo_cpn->cpn_not * df_pay;
            exo_cpn->mdl_val = coupon_iv;
        }
    }

FREE_RETURN:

    if (calc_fwd_iv)
    {
        if (sPricingConst)
        {
            Free_PricingConst(sPricingConst);
            free(sPricingConst);
        }
    }

    if (sNumerParams)
        free(sNumerParams);
    if (cCalibInst)
        free(cCalibInst);
    if (cHestonInst)
        free(cHestonInst);
    if (cCpnSchedule)
        free(cCpnSchedule);

    return err;
}

/*	Adjust model to match market IVs */
Err cts_adjust_model_fwd_iv(
    CTS_UND und,
    CTS     market_cts,
    CTS     model_cts,
    int     for_fund,
    int     calc_fwd_iv,
    int     adj_fee,
    int     pde_or_mc,
    /*	Feedback */
    int    save_fwdiv,
    CTS_IV fwd_iv_info)
{
    int    i, j;
    double fund_iv, mdl_exo_iv, mkt_exo_iv;
    long   today;
    Err    err = NULL;

    /* Get today */
    today = und->mkt->today;

    if (calc_fwd_iv)
    {
        /* allocate memory */
        if (save_fwdiv)
        {
            fwd_iv_info->call_date    = NULL;
            fwd_iv_info->market_fwdiv = NULL;
            fwd_iv_info->model_fwdiv  = NULL;
            fwd_iv_info->extra_fee    = NULL;

            fwd_iv_info->call_date    = (long*)calloc(model_cts->num_calls, sizeof(long));
            fwd_iv_info->market_fwdiv = (double*)calloc(model_cts->num_calls, sizeof(double));
            fwd_iv_info->model_fwdiv  = (double*)calloc(model_cts->num_calls, sizeof(double));
            fwd_iv_info->extra_fee    = (double*)calloc(model_cts->num_calls, sizeof(double));

            if (!fwd_iv_info->call_date || !fwd_iv_info->market_fwdiv ||
                !fwd_iv_info->model_fwdiv || !fwd_iv_info->extra_fee)
            {
                err = "Memory allocation faillure in cts_adjust_model_fwd_iv";
                goto FREE_RETURN;
            }

            fwd_iv_info->ncall = model_cts->num_calls;
        }

        for (i = 0; i < model_cts->num_calls; i++)
        {
            mdl_exo_iv = 0.0;
            mkt_exo_iv = 0.0;

            for (j = model_cts->call[i].exo_idx; j < model_cts->exo_leg->num_cpn; j++)
            {
                mdl_exo_iv += model_cts->exo_leg->cpn[j].mdl_val;
                mkt_exo_iv += market_cts->exo_leg->cpn[j].mkt_val;
            }

            fund_iv = swp_f_df(
                today,
                model_cts->fund_leg->cpn[model_cts->call[i].fund_idx].start_date,
                und->mkt->yc);
            fund_iv *= model_cts->fund_leg->cpn[model_cts->call[i].fund_idx].not ;

            for (j = model_cts->call[i].fund_idx; j < model_cts->fund_leg->num_cpn; j++)
            {
                fund_iv += model_cts->fund_leg->cpn[j].mkt_val;
            }

            model_cts->call[i].extra_fee =
                model_cts->call[i].pay_rec * (mdl_exo_iv - mkt_exo_iv) /
                swp_f_df(today, model_cts->call[i].set_date, und->mkt->yc);

            if (adj_fee)
            {
                model_cts->call[i].total_fee =
                    model_cts->call[i].fee + model_cts->call[i].extra_fee;
            }
            else
            {
                model_cts->call[i].total_fee = model_cts->call[i].fee;
            }

            if (save_fwdiv)
            {
                fwd_iv_info->call_date[i]    = model_cts->call[i].ex_date;
                fwd_iv_info->market_fwdiv[i] = model_cts->call[i].pay_rec * (mkt_exo_iv - fund_iv);
                fwd_iv_info->model_fwdiv[i]  = model_cts->call[i].pay_rec * (mdl_exo_iv - fund_iv);
                fwd_iv_info->extra_fee[i]    = model_cts->call[i].extra_fee;
            }
        }

        if (pde_or_mc == 1)
        {
            /* readjust the fees to be coupon by coupon */
            for (i = 0; i < model_cts->num_calls - 1; i++)
            {
                model_cts->call[i].total_fee -=
                    model_cts->call[i + 1].total_fee *
                    swp_f_df(
                        model_cts->call[i].set_date, model_cts->call[i + 1].set_date, und->mkt->yc);
            }
        }
    }
    else
    {
        for (i = 0; i < model_cts->num_calls; i++)
        {
            model_cts->call[i].extra_fee = 0.0;
            model_cts->call[i].total_fee = model_cts->call[i].fee;
        }
    }

FREE_RETURN:

    return err;
}

/*	Constants used for reconstruction and evaluation  */
/*	------------------------------------------------- */

/*	Fill fixing structure */
Err cts_fill_eval_const_fix(
    CTS_UND und,
    CTS     cts,
    /*	Index of the current cpn */
    int cpn_idx,
    /*	Index of the current fixing */
    int                fix_idx,
    CTS_EVAL_CONST_FIX eval_const_fix)
{
    CTS_FIX fix;
    CTS_CPN cpn;
    double  temp, beta1, beta2, beta1_2, beta2_2;
    int     i;
    long    today;

    /* Get today */
    today = und->mkt->today;

    cpn = &(cts->exo_leg->cpn[cpn_idx]);
    fix = &(cpn->fix[fix_idx]);

    /* reconstruction of DF(t,Tpay) / DF(t,Tstar) */
    eval_const_fix->tpay_tstar_alpha = swp_f_df(today, cpn->cpn_pay_date, und->mkt->yc) /
                                       swp_f_df(today, und->model.lTStarDate, und->mkt->yc);
    eval_const_fix->tpay_tstar_alpha = -log(eval_const_fix->tpay_tstar_alpha);
    eval_const_fix->tpay_tstar_beta =
        (1.0 - exp(-und->model.dLambdaX * (cpn->cpn_pay_time - und->model.dTStar))) /
        und->model.dLambdaX;
    eval_const_fix->tpay_tstar_gamma =
        0.5 * eval_const_fix->tpay_tstar_beta * eval_const_fix->tpay_tstar_beta;

    /* reconstruction of DF(t,Ts) / DF(t,Te) */
    eval_const_fix->ts_te_alpha = swp_f_df(today, fix->schedule.lCpnDate[0], und->mkt->yc) /
                                  swp_f_df(today, fix->schedule.lCpnDate[1], und->mkt->yc);
    eval_const_fix->ts_te_alpha = -log(eval_const_fix->ts_te_alpha);

    beta1 = (1.0 - exp(-und->model.dLambdaX * (fix->schedule.dCpnTime[0] - und->model.dTStar))) /
            und->model.dLambdaX;
    beta2 = (1.0 - exp(-und->model.dLambdaX * (fix->schedule.dCpnTime[1] - und->model.dTStar))) /
            und->model.dLambdaX;

    eval_const_fix->ts_te_beta  = beta1 - beta2;
    eval_const_fix->ts_te_gamma = 0.5 * eval_const_fix->ts_te_beta * (beta1 + beta2);

    if (und->model.iOne2F == 2)
    {
        /* Fill the Extras */
        eval_const_fix->tpay_tstar_beta2 =
            (1.0 - exp(-und->model.dLambdaX2 * (cpn->cpn_pay_time - und->model.dTStar))) /
            und->model.dLambdaX2;
        eval_const_fix->tpay_tstar_gamma2 =
            0.5 * eval_const_fix->tpay_tstar_beta2 * eval_const_fix->tpay_tstar_beta2;
        eval_const_fix->tpay_tstar_gamma12 =
            eval_const_fix->tpay_tstar_beta * eval_const_fix->tpay_tstar_beta2;

        beta1_2 =
            (1.0 - exp(-und->model.dLambdaX2 * (fix->schedule.dCpnTime[0] - und->model.dTStar))) /
            und->model.dLambdaX2;
        beta2_2 =
            (1.0 - exp(-und->model.dLambdaX2 * (fix->schedule.dCpnTime[1] - und->model.dTStar))) /
            und->model.dLambdaX2;

        eval_const_fix->ts_te_beta2   = beta1_2 - beta2_2;
        eval_const_fix->ts_te_gamma2  = 0.5 * eval_const_fix->ts_te_beta2 * (beta1_2 + beta2_2);
        eval_const_fix->ts_te_gamma12 = beta1 * beta1_2 - beta2 * beta2_2;
    }

    /* special flag for initialisation of the payoff column 2 */
    if (cpn->type == 1 && fix_idx == cpn->used_nfix - 1)
    {
        eval_const_fix->do_init_col2 = 1;
    }
    else
    {
        eval_const_fix->do_init_col2 = 0;
    }

    /*	Precalculations */
    eval_const_fix->a =
        cpn->alpha + cpn->beta * (fix->ref_fwd_spread - 1.0 / fix->schedule.dCpnCvg[1]);
    eval_const_fix->b = cpn->beta / fix->schedule.dCpnCvg[1];

    if (cpn->type == 0)
    {
        /* we include the coupon in the fixings */
        temp = cpn->used_fixweights[fix_idx] * cpn->cpn_coupon * cpn->cpn_cvg * cpn->cpn_not;
    }
    else
    {
        /* we keep the gearing outside the fixing */
        temp = cpn->used_fixweights[fix_idx] * cpn->cpn_cvg * cpn->cpn_not;
    }

    eval_const_fix->a *= temp;
    eval_const_fix->b *= temp;

    for (i = 0; i < cpn->nstr; i++)
    {
        eval_const_fix->s[i] = (cpn->str[i] - fix->ref_fwd_spread) * fix->schedule.dCpnCvg[1] + 1.0;
        eval_const_fix->n[i] = cpn->nbopt[i] / fix->schedule.dCpnCvg[1];
        eval_const_fix->n[i] *= temp;
    }

    return NULL;
}

/*	Fill fixing structure */
Err cts_fill_eval_const_fix_cms(
    CTS_UND und,
    CTS     cts,
    /*	Index of the current cpn */
    int cpn_idx,
    /*	Index of the current fixing */
    int                    fix_idx,
    CTS_EVAL_CONST_FIX_CMS eval_const_fix_cms)
{
    CTS_FIX fix;
    CTS_CPN cpn;
    double  temp, df_tstar;
    int     i;
    long    today;

    /* Get today */
    today = und->mkt->today;

    cpn = &(cts->exo_leg->cpn[cpn_idx]);
    fix = &(cpn->fix[fix_idx]);

    df_tstar = swp_f_df(today, und->model.lTStarDate, und->mkt->yc);

    /* reconstruction of DF(t,Tpay) / DF(t,Tstar) */
    eval_const_fix_cms->tpay_tstar_alpha =
        swp_f_df(today, cpn->cpn_pay_date, und->mkt->yc) / df_tstar;
    eval_const_fix_cms->tpay_tstar_alpha = -log(eval_const_fix_cms->tpay_tstar_alpha);
    eval_const_fix_cms->tpay_tstar_beta =
        (1.0 - exp(-und->model.dLambdaX * (cpn->cpn_pay_time - und->model.dTStar))) /
        und->model.dLambdaX;
    eval_const_fix_cms->tpay_tstar_gamma =
        0.5 * eval_const_fix_cms->tpay_tstar_beta * eval_const_fix_cms->tpay_tstar_beta;

    /* reconstruction of DF(t,Ti) / DF(t,Tstar) * coverage ! */
    for (i = 1; i < fix->schedule.iNCpn - 1; i++)
    {
        eval_const_fix_cms->ti_tstar_alpha[i] =
            fix->schedule.dCpnCvg[i] * swp_f_df(today, fix->schedule.lCpnDate[i], und->mkt->yc) /
            df_tstar;
        eval_const_fix_cms->ti_tstar_alpha[i] = -log(eval_const_fix_cms->ti_tstar_alpha[i]);
        eval_const_fix_cms->ti_tstar_beta[i] =
            (1.0 - exp(-und->model.dLambdaX * (fix->schedule.dCpnTime[i] - und->model.dTStar))) /
            und->model.dLambdaX;
        eval_const_fix_cms->ti_tstar_gamma[i] =
            0.5 * eval_const_fix_cms->ti_tstar_beta[i] * eval_const_fix_cms->ti_tstar_beta[i];
    }

    if (fix->schedule.iNCpn > 0)
    {
        /* reconstruction of the first DF without the coverave */
        eval_const_fix_cms->ti_tstar_alpha[0] =
            swp_f_df(today, fix->schedule.lCpnDate[0], und->mkt->yc) / df_tstar;
        eval_const_fix_cms->ti_tstar_alpha[0] = -log(eval_const_fix_cms->ti_tstar_alpha[0]);
        eval_const_fix_cms->ti_tstar_beta[0] =
            (1.0 - exp(-und->model.dLambdaX * (fix->schedule.dCpnTime[0] - und->model.dTStar))) /
            und->model.dLambdaX;
        eval_const_fix_cms->ti_tstar_gamma[0] =
            0.5 * eval_const_fix_cms->ti_tstar_beta[0] * eval_const_fix_cms->ti_tstar_beta[0];

        /* reconstruction of the last DF without the coverave */
        eval_const_fix_cms->ti_tstar_alpha[fix->schedule.iNCpn - 1] =
            swp_f_df(today, fix->schedule.lCpnDate[fix->schedule.iNCpn - 1], und->mkt->yc) /
            df_tstar;
        eval_const_fix_cms->ti_tstar_alpha[fix->schedule.iNCpn - 1] =
            -log(eval_const_fix_cms->ti_tstar_alpha[fix->schedule.iNCpn - 1]);
        eval_const_fix_cms->ti_tstar_beta[fix->schedule.iNCpn - 1] =
            (1.0 - exp(-und->model.dLambdaX *
                       (fix->schedule.dCpnTime[fix->schedule.iNCpn - 1] - und->model.dTStar))) /
            und->model.dLambdaX;
        eval_const_fix_cms->ti_tstar_gamma[fix->schedule.iNCpn - 1] =
            0.5 * eval_const_fix_cms->ti_tstar_beta[fix->schedule.iNCpn - 1] *
            eval_const_fix_cms->ti_tstar_beta[fix->schedule.iNCpn - 1];
    }

    if (und->model.iOne2F == 2)
    {
        /* Fill the Extras */
        eval_const_fix_cms->tpay_tstar_beta2 =
            (1.0 - exp(-und->model.dLambdaX2 * (cpn->cpn_pay_time - und->model.dTStar))) /
            und->model.dLambdaX2;
        eval_const_fix_cms->tpay_tstar_gamma2 =
            0.5 * eval_const_fix_cms->tpay_tstar_beta2 * eval_const_fix_cms->tpay_tstar_beta2;
        eval_const_fix_cms->tpay_tstar_gamma12 =
            eval_const_fix_cms->tpay_tstar_beta * eval_const_fix_cms->tpay_tstar_beta2;

        for (i = 1; i < fix->schedule.iNCpn - 1; i++)
        {
            eval_const_fix_cms->ti_tstar_beta2[i] =
                (1.0 -
                 exp(-und->model.dLambdaX2 * (fix->schedule.dCpnTime[i] - und->model.dTStar))) /
                und->model.dLambdaX2;
            eval_const_fix_cms->ti_tstar_gamma2[i] =
                0.5 * eval_const_fix_cms->ti_tstar_beta2[i] * eval_const_fix_cms->ti_tstar_beta2[i];
            eval_const_fix_cms->ti_tstar_gamma12[i] =
                eval_const_fix_cms->ti_tstar_beta[i] * eval_const_fix_cms->ti_tstar_beta2[i];
        }

        if (fix->schedule.iNCpn > 0)
        {
            /* reconstruction of the first DF without the coverave */
            eval_const_fix_cms->ti_tstar_beta2[0] =
                (1.0 -
                 exp(-und->model.dLambdaX2 * (fix->schedule.dCpnTime[0] - und->model.dTStar))) /
                und->model.dLambdaX2;
            eval_const_fix_cms->ti_tstar_gamma2[0] =
                0.5 * eval_const_fix_cms->ti_tstar_beta2[0] * eval_const_fix_cms->ti_tstar_beta2[0];
            eval_const_fix_cms->ti_tstar_gamma12[0] =
                eval_const_fix_cms->ti_tstar_beta[0] * eval_const_fix_cms->ti_tstar_beta2[0];

            /* reconstruction of the last DF without the coverave */
            eval_const_fix_cms->ti_tstar_beta2[fix->schedule.iNCpn - 1] =
                (1.0 - exp(-und->model.dLambdaX2 *
                           (fix->schedule.dCpnTime[fix->schedule.iNCpn - 1] - und->model.dTStar))) /
                und->model.dLambdaX2;
            eval_const_fix_cms->ti_tstar_gamma2[fix->schedule.iNCpn - 1] =
                0.5 * eval_const_fix_cms->ti_tstar_beta2[fix->schedule.iNCpn - 1] *
                eval_const_fix_cms->ti_tstar_beta2[fix->schedule.iNCpn - 1];
            eval_const_fix_cms->ti_tstar_gamma12[fix->schedule.iNCpn - 1] =
                eval_const_fix_cms->ti_tstar_beta[fix->schedule.iNCpn - 1] *
                eval_const_fix_cms->ti_tstar_beta2[fix->schedule.iNCpn - 1];
        }
    }

    /* special flag for initialisation of the payoff column 2 */
    if (cpn->type == 1 && fix_idx == cpn->used_nfix - 1)
    {
        eval_const_fix_cms->do_init_col2 = 1;
    }
    else
    {
        eval_const_fix_cms->do_init_col2 = 0;
    }

    /*	Precalculations */
    eval_const_fix_cms->a = cpn->alpha + cpn->beta * fix->ref_fwd_spread;
    eval_const_fix_cms->b = cpn->beta;

    if (cpn->type == 0)
    {
        /* we include the coupon in the fixings */
        temp = cpn->used_fixweights[fix_idx] * cpn->cpn_coupon * cpn->cpn_cvg * cpn->cpn_not;
    }
    else
    {
        /* we keep the gearing outside the fixing */
        temp = cpn->used_fixweights[fix_idx] * cpn->cpn_cvg * cpn->cpn_not;
    }

    eval_const_fix_cms->a *= temp;
    eval_const_fix_cms->b *= temp;

    for (i = 0; i < cpn->nstr; i++)
    {
        eval_const_fix_cms->s[i] = (cpn->str[i] - fix->ref_fwd_spread);
        eval_const_fix_cms->n[i] = cpn->nbopt[i];
        eval_const_fix_cms->n[i] *= temp;
    }

    return NULL;
}

/*	Fill coupon structure */
Err cts_fill_eval_const_cpn(
    CTS_UND und,
    CTS     cts,
    /*	Index of the current cpn */
    int                cpn_idx,
    CTS_EVAL_CONST_CPN eval_const_cpn)
{
    CTS_CPN cpn;
    double  beta1, beta2, beta1_2, beta2_2;
    long    today;

    /* Get today */
    today = und->mkt->today;

    cpn = &(cts->exo_leg->cpn[cpn_idx]);

    /* reconstruction of DF(t,Tpay) / DF(t,Tstar) */
    eval_const_cpn->tpay_tstar_alpha = swp_f_df(today, cpn->cpn_pay_date, und->mkt->yc) /
                                       swp_f_df(today, und->model.lTStarDate, und->mkt->yc);
    eval_const_cpn->tpay_tstar_alpha = -log(eval_const_cpn->tpay_tstar_alpha);
    eval_const_cpn->tpay_tstar_beta =
        (1.0 - exp(-und->model.dLambdaX * (cpn->cpn_pay_time - und->model.dTStar))) /
        und->model.dLambdaX;
    eval_const_cpn->tpay_tstar_gamma =
        0.5 * eval_const_cpn->tpay_tstar_beta * eval_const_cpn->tpay_tstar_beta;

    /* reconstruction of DF(t,Ts) / DF(t,Te) */
    eval_const_cpn->ts_te_alpha = swp_f_df(today, cpn->pay_start_date, und->mkt->yc) /
                                  swp_f_df(today, cpn->pay_end_date, und->mkt->yc);
    eval_const_cpn->ts_te_alpha = -log(eval_const_cpn->ts_te_alpha);

    beta1 = (1.0 - exp(-und->model.dLambdaX * (cpn->pay_start_time - und->model.dTStar))) /
            und->model.dLambdaX;
    beta2 = (1.0 - exp(-und->model.dLambdaX * (cpn->pay_end_time - und->model.dTStar))) /
            und->model.dLambdaX;

    eval_const_cpn->ts_te_beta  = beta1 - beta2;
    eval_const_cpn->ts_te_gamma = 0.5 * eval_const_cpn->ts_te_beta * (beta1 + beta2);

    if (und->model.iOne2F == 2)
    {
        /* Fill the Extras */
        eval_const_cpn->tpay_tstar_beta2 =
            (1.0 - exp(-und->model.dLambdaX2 * (cpn->cpn_pay_time - und->model.dTStar))) /
            und->model.dLambdaX2;
        eval_const_cpn->tpay_tstar_gamma2 =
            0.5 * eval_const_cpn->tpay_tstar_beta2 * eval_const_cpn->tpay_tstar_beta2;
        eval_const_cpn->tpay_tstar_gamma12 =
            eval_const_cpn->tpay_tstar_beta * eval_const_cpn->tpay_tstar_beta2;

        beta1_2 = (1.0 - exp(-und->model.dLambdaX2 * (cpn->pay_start_time - und->model.dTStar))) /
                  und->model.dLambdaX2;
        beta2_2 = (1.0 - exp(-und->model.dLambdaX2 * (cpn->pay_end_time - und->model.dTStar))) /
                  und->model.dLambdaX2;

        eval_const_cpn->ts_te_beta2   = beta1_2 - beta2_2;
        eval_const_cpn->ts_te_gamma2  = 0.5 * eval_const_cpn->ts_te_beta2 * (beta1_2 + beta2_2);
        eval_const_cpn->ts_te_gamma12 = beta1 * beta1_2 - beta2 * beta2_2;
    }

    /* For the geared Libor reconstruction : Libor = (1/cvg * DF(s)/DF(e) - 1/cvg) */
    eval_const_cpn->shift =
        cpn->pay_gearing * (-1.0 / cpn->pay_cvg + cpn->pay_fwd_spread) + cpn->cpn_coupon;
    eval_const_cpn->multip = cpn->pay_gearing / cpn->pay_cvg;

    return NULL;
}

/*	Fill call structure */
Err cts_fill_eval_const_call(
    CTS_UND und,
    CTS     cts,
    /*	Index of the current call */
    int                 call_idx,
    CTS_EVAL_CONST_CALL eval_const_call)
{
    long         today;
    CTS_CALL     call;
    CTS_FUND_LEG fund;
    CTS_FUND_CPN cpn;
    CTS_CPN      exo_cpn;
    int          i;

    /* Get today */
    today = und->mkt->today;

    call = &(cts->call[call_idx]);
    fund = cts->fund_leg;

    /* reconstruction of DF(t, fund) / DF(t, Tstar) */

    for (i = 0; i < call->num_fund_cpn; i++)
    {
        cpn = &(fund->cpn[call->fund_idx + i]);

        eval_const_call->tpay_tstar_alpha[i] = swp_f_df(today, cpn->pay_date, und->mkt->yc) /
                                               swp_f_df(today, und->model.lTStarDate, und->mkt->yc);
        eval_const_call->tpay_tstar_alpha[i] = -log(eval_const_call->tpay_tstar_alpha[i]);
        eval_const_call->tpay_tstar_beta[i] =
            (1.0 - exp(-und->model.dLambdaX * (cpn->pay_time - und->model.dTStar))) /
            und->model.dLambdaX;
        eval_const_call->tpay_tstar_gamma[i] =
            0.5 * eval_const_call->tpay_tstar_beta[i] * eval_const_call->tpay_tstar_beta[i];
    }

    /* then add the first discount factor at start date */
    cpn = &(fund->cpn[call->fund_idx]);

    i                                    = call->num_fund_cpn;
    eval_const_call->tpay_tstar_alpha[i] = swp_f_df(today, cpn->start_date, und->mkt->yc) /
                                           swp_f_df(today, und->model.lTStarDate, und->mkt->yc);
    eval_const_call->tpay_tstar_alpha[i] = -log(eval_const_call->tpay_tstar_alpha[i]);
    eval_const_call->tpay_tstar_beta[i] =
        (1.0 - exp(-und->model.dLambdaX * (cpn->start_time - und->model.dTStar))) /
        und->model.dLambdaX;
    eval_const_call->tpay_tstar_gamma[i] =
        0.5 * eval_const_call->tpay_tstar_beta[i] * eval_const_call->tpay_tstar_beta[i];

    /* reconstruction of DF(t,Tset) / DF(t,Tstar) */
    eval_const_call->tset_tstar_alpha = swp_f_df(today, call->set_date, und->mkt->yc) /
                                        swp_f_df(today, und->model.lTStarDate, und->mkt->yc);
    eval_const_call->tset_tstar_alpha = -log(eval_const_call->tset_tstar_alpha);
    eval_const_call->tset_tstar_beta =
        (1.0 - exp(-und->model.dLambdaX * (call->set_time - und->model.dTStar))) /
        und->model.dLambdaX;
    eval_const_call->tset_tstar_gamma =
        0.5 * eval_const_call->tset_tstar_beta * eval_const_call->tset_tstar_beta;

    /* in the midat case, evaluation of the fixed coupon leg */
    if (call->is_midat)
    {
        for (i = 0; i < call->num_exo_cpn; i++)
        {
            exo_cpn = &(cts->exo_leg->cpn[call->exo_idx + i]);

            eval_const_call->tpay_tstar_alpha2[i] =
                swp_f_df(today, exo_cpn->cpn_pay_date, und->mkt->yc) /
                swp_f_df(today, und->model.lTStarDate, und->mkt->yc);
            eval_const_call->tpay_tstar_alpha2[i] *=
                exo_cpn->cpn_cvg * exo_cpn->cpn_coupon * exo_cpn->cpn_not;
            eval_const_call->tpay_tstar_alpha2[i] = -log(eval_const_call->tpay_tstar_alpha2[i]);
            eval_const_call->tpay_tstar_beta2[i] =
                (1.0 - exp(-und->model.dLambdaX * (exo_cpn->cpn_pay_time - und->model.dTStar))) /
                und->model.dLambdaX;
            eval_const_call->tpay_tstar_gamma2[i] =
                0.5 * eval_const_call->tpay_tstar_beta2[i] * eval_const_call->tpay_tstar_beta2[i];
        }
    }

    if (und->model.iOne2F == 2)
    {
        /* Fill the Extras */
        for (i = 0; i < call->num_fund_cpn; i++)
        {
            cpn = &(fund->cpn[call->fund_idx + i]);

            eval_const_call->tpay_tstar_beta_2[i] =
                (1.0 - exp(-und->model.dLambdaX2 * (cpn->pay_time - und->model.dTStar))) /
                und->model.dLambdaX2;
            eval_const_call->tpay_tstar_gamma_2[i] =
                0.5 * eval_const_call->tpay_tstar_beta_2[i] * eval_const_call->tpay_tstar_beta_2[i];
            eval_const_call->tpay_tstar_gamma_12[i] =
                eval_const_call->tpay_tstar_beta[i] * eval_const_call->tpay_tstar_beta_2[i];
        }

        /* then add the first discount factor at start date */
        cpn = &(fund->cpn[call->fund_idx]);

        i = call->num_fund_cpn;
        eval_const_call->tpay_tstar_beta_2[i] =
            (1.0 - exp(-und->model.dLambdaX2 * (cpn->start_time - und->model.dTStar))) /
            und->model.dLambdaX2;
        eval_const_call->tpay_tstar_gamma_2[i] =
            0.5 * eval_const_call->tpay_tstar_beta_2[i] * eval_const_call->tpay_tstar_beta_2[i];
        eval_const_call->tpay_tstar_gamma_12[i] =
            eval_const_call->tpay_tstar_beta[i] * eval_const_call->tpay_tstar_beta_2[i];

        /* reconstruction of DF(t,Tset) / DF(t,Tstar) */
        eval_const_call->tset_tstar_beta2 =
            (1.0 - exp(-und->model.dLambdaX2 * (call->set_time - und->model.dTStar))) /
            und->model.dLambdaX2;
        eval_const_call->tset_tstar_gamma2 =
            0.5 * eval_const_call->tset_tstar_beta2 * eval_const_call->tset_tstar_beta2;
        eval_const_call->tset_tstar_gamma12 =
            eval_const_call->tset_tstar_beta * eval_const_call->tset_tstar_beta2;

        if (call->is_midat)
        {
            for (i = 0; i < call->num_exo_cpn; i++)
            {
                exo_cpn = &(cts->exo_leg->cpn[call->exo_idx + i]);

                eval_const_call->tpay_tstar_beta2_2[i] =
                    (1.0 -
                     exp(-und->model.dLambdaX2 * (exo_cpn->cpn_pay_time - und->model.dTStar))) /
                    und->model.dLambdaX2;
                eval_const_call->tpay_tstar_gamma2_2[i] = 0.5 *
                                                          eval_const_call->tpay_tstar_beta2_2[i] *
                                                          eval_const_call->tpay_tstar_beta2_2[i];
                eval_const_call->tpay_tstar_gamma2_12[i] =
                    eval_const_call->tpay_tstar_beta2[i] * eval_const_call->tpay_tstar_beta2_2[i];
            }
        }
    }

    return NULL;
}

/*	function to guess what is the most expensive one-time callable */
Err cts_find_one_time_index(CTS_UND und, CTS cts, int* one_time_index)
{
    double exe_mat;
    int    i;
    Err    err = NULL;

    if (*one_time_index == 0)
    {
        exe_mat = cts->exo_leg->cpn[cts->exo_leg->num_cpn - 1].cpn_pay_time / 4.0;

        i = 0;
        while (i < cts->num_calls && cts->call[i].ex_time < exe_mat)
        {
            i++;
        }

        if (i == 0)
            i++;

        *one_time_index = i;
    }

    return err;
}

/*	Arguments to all payoff evaluation functions */
/*	-------------------------------------------- */

Err cts_fill_algo_arg(
    CTS_UND und,
    CTS     cts,
    int     pde_or_mc, /* 0: PDE, 1: MC */
    /*	Required number of steps for PDE */
    int req_stp,
    int req_stppsi,
    int req_stpx,
    int req_stpz,
    /*	Required number of steps for MC */
    double req_mintime,
    long   req_paths,
    /*	Extra numerical parameter */
    double integ_mintime,
    /*	Flag for extra calculation / adjustment of one-time callable */
    int do_one_time,    /*	1: calc the one time */
    int one_time_index, /*	0: choose automatically the index, >0: index provided by user */
    /*	T-star */
    CTS_ADI_ARG adi_arg)
{
    Err         err = NULL;
    CTS_PAY_ARG cts_prm;
    CTS_CALL    call;
    CTS_CPN     cpn, cpn_for_fix;
    long        today;
    double      times[MAXTS];
    int         i, j, k;
    int         index, index_call, index_cpn, index_fix, index_cpn_for_fix, index_event;
    int         is_midat, is_floating;
    int         init_nstp, nstp;
    double      local_min_time;

    /*	Initialise */
    adi_arg->time = NULL;
    adi_arg->date = NULL;

    adi_arg->void_prm = NULL;
    adi_arg->is_event = NULL;
    adi_arg->ifr      = NULL;

    adi_arg->pde_or_mc = pde_or_mc;

    adi_arg->initnstp = req_stp;
    adi_arg->nstpx    = req_stpx;
    adi_arg->nstppsi  = req_stppsi;
    adi_arg->nstpz    = req_stpz;

    adi_arg->npaths     = req_paths;
    adi_arg->mc_mintime = req_mintime;

    /*	Get today */
    today = und->mkt->today;

    /*	Midat flag */
    is_midat = cts->call[0].is_midat;

    /*	Times */

    /*	Get fixing times and coupon times */
    k = 0;

    if (is_midat == 0)
    {
        /* include coupon and fixing time */
        for (i = cts->call[0].exo_idx; i < cts->exo_leg->num_cpn; i++)
        {
            for (j = 0; j < cts->exo_leg->cpn[i].used_nfix; j++)
            {
                times[k] = cts->exo_leg->cpn[i].fix[j].ref_fix_time;
                k++;
            }
            if (cts->exo_leg->cpn[i].type)
            {
                times[k] = cts->exo_leg->cpn[i].pay_fix_time;
                k++;
            }
        }
    }

    /*	Get call times */
    for (i = 0; i < cts->num_calls; i++)
    {
        times[k] = cts->call[i].ex_time;
        k++;
    }

    adi_arg->nstp = k;
    adi_arg->time = (double*)calloc(adi_arg->nstp, sizeof(double));
    if (!adi_arg->time)
    {
        err = "Memory allocation error (1) in cts_fill_algo_arg";
        goto FREE_RETURN;
    }

    /*	Copy and sort */
    memcpy(adi_arg->time, times, adi_arg->nstp * sizeof(double));
    num_f_sort_vector(adi_arg->nstp, adi_arg->time);
    num_f_unique_vector(&(adi_arg->nstp), adi_arg->time);

    adi_arg->nb_event = adi_arg->nstp;

    /*	Fill time vector */

    /*	Add today if required */
    if (adi_arg->time[0] < -EPS)
    {
        err = "Past event date in cts_fill_algo_arg";
        goto FREE_RETURN;
    }
    if (adi_arg->time[0] > EPS)
    {
        num_f_add_number(&(adi_arg->nstp), &(adi_arg->time), 0.0);
        num_f_sort_vector(adi_arg->nstp, adi_arg->time);
        num_f_unique_vector(&(adi_arg->nstp), adi_arg->time);
    }

    /*	If only one event today, add empty event */
    if (adi_arg->nstp == 1)
    {
        num_f_add_number(&(adi_arg->nstp), (&adi_arg->time), 1.0);
    }

    /*	Fill the vector */
    if (pde_or_mc == 0)
    {
        num_f_fill_vector_newalgo(&(adi_arg->nstp), &(adi_arg->time), req_stp);
    }
    else
    {
        memcpy(times, adi_arg->time, adi_arg->nstp * sizeof(double));

        init_nstp = adi_arg->nstp;
        nstp      = 1;

        for (i = 0; i < adi_arg->nstp - 1; i++)
        {
            /* Fill between times[i] and times[i+1] */
            nstp += max((int)((times[i + 1] - times[i]) / adi_arg->mc_mintime + 0.5), 1);
        }

        free(adi_arg->time);

        adi_arg->nstp       = nstp;
        adi_arg->time       = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->optimise   = (int*)calloc(adi_arg->nb_event, sizeof(int));
        adi_arg->mcebparams = (MCEBParams*)calloc(1, sizeof(MCEBParams));

        adi_arg->dSigma     = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dAlpha     = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dLambdaEps = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dLvlEps    = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dRho       = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dRho2      = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dLGMAlpha  = (double*)calloc(adi_arg->nstp, sizeof(double));
        adi_arg->dLGMRho    = (double*)calloc(adi_arg->nstp, sizeof(double));

        if (!adi_arg->time || !adi_arg->optimise || !adi_arg->mcebparams || !adi_arg->dSigma ||
            !adi_arg->dAlpha || !adi_arg->dLambdaEps || !adi_arg->dLvlEps || !adi_arg->dRho ||
            !adi_arg->dRho2 || !adi_arg->dLGMAlpha || !adi_arg->dLGMRho)
        {
            err = "Memory allocation error (2) in cts_fill_algo_arg";
            goto FREE_RETURN;
        }

        /* MCEBPArams */
        mceb_set_default_params(adi_arg->mcebparams);

        adi_arg->mcebparams->iIsKO             = 1;
        adi_arg->mcebparams->iColPay           = 0;
        adi_arg->mcebparams->iColBound         = 1;
        adi_arg->mcebparams->iMultiIndex       = und->model.iOne2F;
        adi_arg->mcebparams->iNbIndex          = und->model.iOne2F + 1;
        adi_arg->mcebparams->iRemoveLastOnLast = und->model.iOne2F - 1;

        adi_arg->time[0] = 0.0;
        j                = 1;

        for (i = 0; i < init_nstp - 1; i++)
        {
            /* Fill between times[i] and times[i+1] */
            nstp           = max((int)((times[i + 1] - times[i]) / adi_arg->mc_mintime + 0.5), 1);
            local_min_time = (times[i + 1] - times[i]) / (nstp * 1.0);

            for (k = 0; k < nstp - 1; k++)
            {
                adi_arg->time[j] = adi_arg->time[j - 1] + local_min_time;
                j++;
            }

            adi_arg->time[j] = times[i + 1];
            j++;
        }

        /* Fill the und information */
        for (i = 0; i < adi_arg->nstp; i++)
        {
            index = Get_Index(adi_arg->time[i], und->model.dPWTime, und->model.iNbPWTime);
            adi_arg->dSigma[i]     = und->model.dSigma[index];
            adi_arg->dAlpha[i]     = und->model.dAlpha[index];
            adi_arg->dLambdaEps[i] = und->model.dLambdaEps[index];
            adi_arg->dLvlEps[i]    = und->model.dLvlEps[index];
            adi_arg->dRho[i]       = und->model.dRho[index];

            if (und->model.iOne2F == 2)
            {
                adi_arg->dRho2[i]     = und->model.dRho2[index];
                adi_arg->dLGMAlpha[i] = und->model.dLGMAlpha[index];
                adi_arg->dLGMRho[i]   = und->model.dLGMRho[index];
            }
        }
    }

    /*	Make dates */
    adi_arg->date = (double*)calloc(adi_arg->nstp, sizeof(double));
    if (!adi_arg->date)
    {
        err = "Memory allocation error (2) in cts_fill_algo_arg";
        goto FREE_RETURN;
    }

    for (i = 0; i < adi_arg->nstp; i++)
    {
        adi_arg->date[i] = today + DAYS_IN_YEAR * adi_arg->time[i];

        /*
        if (i > 0 && adi_arg->date[i] - adi_arg->date[i-1] >= 1)
        {
                adi_arg->date[i] = (long) (adi_arg->date[i] + 1.0e-08);
                adi_arg->time[i] = YEARS_IN_DAY * (adi_arg->date[i] - today);
        }
        */
    }

    /*	IFR */

    adi_arg->ifr = (double*)calloc(adi_arg->nstp, sizeof(double));
    if (!adi_arg->ifr)
    {
        err = "Memory allocation error (3) in cts_fill_algo_arg";
        goto FREE_RETURN;
    }

    for (i = 0; i < adi_arg->nstp - 1; i++)
    {
        adi_arg->ifr[i] = swp_f_zr(adi_arg->date[i], adi_arg->date[i + 1], und->mkt->yc);
    }

    /* Find the type of the CTS */
    is_floating = 0;
    for (i = 0; i < cts->call[0].num_exo_cpn; i++)
    {
        if (cts->exo_leg->cpn[cts->call[0].exo_idx + i].type == 1 ||
            cts->exo_leg->cpn[cts->call[0].exo_idx + i].type == 2)
        {
            is_floating = 1;
            break;
        }
    }

    if (is_floating && pde_or_mc == 1)
    {
        /* Extra allocation for path dependent informations */
        adi_arg->dPathInfos = (double*)calloc(req_paths + 1, sizeof(double));

        if (!adi_arg->dPathInfos)
        {
            err = "Memory allocation faillure in cts_fill_algo_arg";
            goto FREE_RETURN;
        }
    }
    else
    {
        adi_arg->dPathInfos = NULL;
    }

    /*	Fill limit conditions (product) */

    adi_arg->is_event = (int*)calloc(adi_arg->nstp, sizeof(int));
    adi_arg->void_prm = (void**)calloc(adi_arg->nstp, sizeof(void*));

    if (!adi_arg->is_event || !adi_arg->void_prm)
    {
        err = "Memory allocation error (4) in cts_fill_algo_arg";
        goto FREE_RETURN;
    }

    index_call = 0;
    call       = cts->call;

    index_cpn = call->exo_idx;
    cpn       = &(cts->exo_leg->cpn[index_cpn]);

    index_fix         = 0;
    index_cpn_for_fix = index_cpn;
    cpn_for_fix       = cpn;
    index_event       = 0;

    adi_arg->has_one_time = 0;

    for (i = 0; i < adi_arg->nstp; i++)
    {
        adi_arg->is_event[i] = 0;
        adi_arg->void_prm[i] = NULL;

        cts_prm = (cts_pay_arg*)calloc(1, sizeof(cts_pay_arg));

        if (!cts_prm)
        {
            err = "memory allocation faillure in cts_fill_algo_arg";
            goto FREE_RETURN;
        }

        cts_prm->is_call       = 0;
        cts_prm->is_cpn        = 0;
        cts_prm->is_fixing     = 0;
        cts_prm->is_fixing_cms = 0;

        if (index_call < cts->num_calls && fabs(adi_arg->date[i] - call->ex_date) < 1.0e-08)
        {
            /* this is a call date */
            cts_prm->und     = und;
            cts_prm->cts     = cts;
            cts_prm->call    = call;
            cts_prm->is_call = 1;

            err = cts_fill_eval_const_call(und, cts, index_call, &(cts_prm->eval_const_call));

            if (err)
                goto FREE_RETURN;

            if (do_one_time && (index_call + 1 == one_time_index ||
                                (one_time_index < 0 && index_call + 1 <= -one_time_index)))
            {
                cts_prm->eval_const_call.eval_one_time = 1;
                adi_arg->has_one_time += 1;
            }
            else
            {
                cts_prm->eval_const_call.eval_one_time = 0;
            }

            adi_arg->is_event[i] = 1;
            adi_arg->void_prm[i] = (void*)cts_prm;

            if (pde_or_mc == 1 && adi_arg->optimise)
            {
                adi_arg->optimise[index_event] = 1;
            }

            call++;
            index_call++;
        }

        if (index_cpn < cts->exo_leg->num_cpn &&
            fabs(adi_arg->date[i] - cts->exo_leg->cpn[index_cpn].pay_fix_date) < 1.0e-08 &&
            is_midat == 0)
        {
            /* this is a cpn date */
            /* this is an event date only in the floating coupon case */
            if (cpn->type == 1 || cpn->type == 2)
            {
                cts_prm->und    = und;
                cts_prm->cts    = cts;
                cts_prm->cpn    = cpn;
                cts_prm->is_cpn = 1;

                err = cts_fill_eval_const_cpn(und, cts, index_cpn, &(cts_prm->eval_const_cpn));

                if (err)
                    goto FREE_RETURN;

                adi_arg->is_event[i] = 1;
                adi_arg->void_prm[i] = (void*)cts_prm;
            }
            else
            {
                adi_arg->is_event[i] = 0;
                adi_arg->void_prm[i] = NULL;
            }

            cpn++;
            index_cpn++;
        }

        if (index_cpn_for_fix < cts->exo_leg->num_cpn && index_fix < cpn_for_fix->used_nfix &&
            fabs(adi_arg->date[i] - cpn_for_fix->fix[index_fix].ref_fix_date) < 1.0e-08 &&
            is_midat == 0)
        {
            /* this is a fixing date */
            cts_prm->und = und;
            cts_prm->cts = cts;
            cts_prm->cpn = cpn_for_fix;
            cts_prm->fix = &(cpn_for_fix->fix[index_fix]);

            if (cts_prm->fix->schedule.iNCpn == 2)
            {
                cts_prm->is_fixing = 1;

                err = cts_fill_eval_const_fix(
                    und, cts, index_cpn_for_fix, index_fix, &(cts_prm->eval_const_fix));
            }
            else
            {
                cts_prm->is_fixing_cms = 1;

                err = cts_fill_eval_const_fix_cms(
                    und, cts, index_cpn_for_fix, index_fix, &(cts_prm->eval_const_fix_cms));
            }

            if (err)
                goto FREE_RETURN;

            adi_arg->is_event[i] = 1;
            adi_arg->void_prm[i] = (void*)cts_prm;

            index_fix++;

            if (index_fix == cpn_for_fix->used_nfix)
            {
                /* switch to the next coupon */

                if (cpn_for_fix->type == 0 || cpn_for_fix->type == 3)
                {
                    index_cpn++;
                    cpn++;
                }

                index_cpn_for_fix++;
                cpn_for_fix++;

                index_fix = 0;
            }
        }

        cts_prm->dPathInfos = adi_arg->dPathInfos;

        if (!adi_arg->is_event[i])
        {
            if (cts_prm)
                free(cts_prm);
        }

        cts_prm = NULL;
        index_event += adi_arg->is_event[i];
    }

    /* check that all the dates have been caught */
    if (index_call < cts->num_calls || (index_cpn_for_fix < cts->exo_leg->num_cpn && is_midat == 0))
    {
        err = "All event dates are not in the time discretisation";
        goto FREE_RETURN;
    }

    /* fill the PDE parameters */
    adi_arg->params = (LGMSVParam*)calloc(1, sizeof(LGMSVParam));

    if (!adi_arg->params)
    {
        err = "memory allocation error in cts_fill_algo_arg";
        goto FREE_RETURN;
    }

    err = Fill_lgmSV_defaultParam(adi_arg->params);

    adi_arg->params->Tstar              = und->model.dTStar;
    adi_arg->params->dMultiIntegMinTime = integ_mintime;

FREE_RETURN:

    if (err)
    {
        if (cts_prm)
            free(cts_prm);
        cts_prm = NULL;

        cts_free_adi_arg(adi_arg);
    }

    return err;
}

void cts_free_adi_arg(CTS_ADI_ARG adi_arg)
{
    int         i;
    CTS_PAY_ARG cts_prm;

    if (adi_arg)
    {
        if (adi_arg->time)
            free(adi_arg->time);
        if (adi_arg->date)
            free(adi_arg->date);
        if (adi_arg->ifr)
            free(adi_arg->ifr);
        if (adi_arg->is_event)
            free(adi_arg->is_event);

        if (adi_arg->void_prm)
        {
            for (i = 0; i < adi_arg->nstp; i++)
            {
                cts_prm = (CTS_PAY_ARG)(adi_arg->void_prm[i]);
                free(cts_prm);
            }

            free(adi_arg->void_prm);
        }

        if (adi_arg->params)
            free(adi_arg->params);

        if (adi_arg->mcebparams)
        {
            mceb_free_params(adi_arg->mcebparams);
            free(adi_arg->mcebparams);
        }

        if (adi_arg->optimise)
            free(adi_arg->optimise);

        if (adi_arg->dSigma)
            free(adi_arg->dSigma);
        if (adi_arg->dAlpha)
            free(adi_arg->dAlpha);
        if (adi_arg->dLambdaEps)
            free(adi_arg->dLambdaEps);
        if (adi_arg->dLvlEps)
            free(adi_arg->dLvlEps);
        if (adi_arg->dRho)
            free(adi_arg->dRho);
        if (adi_arg->dRho2)
            free(adi_arg->dRho2);
        if (adi_arg->dLGMAlpha)
            free(adi_arg->dLGMAlpha);
        if (adi_arg->dLGMRho)
            free(adi_arg->dLGMRho);

        if (adi_arg->dPathInfos)
            free(adi_arg->dPathInfos);

        adi_arg->time     = NULL;
        adi_arg->date     = NULL;
        adi_arg->ifr      = NULL;
        adi_arg->is_event = NULL;

        adi_arg->void_prm = NULL;
        adi_arg->params   = NULL;

        adi_arg->mcebparams = NULL;
        adi_arg->optimise   = NULL;

        adi_arg->dSigma     = NULL;
        adi_arg->dAlpha     = NULL;
        adi_arg->dLambdaEps = NULL;
        adi_arg->dLvlEps    = NULL;
        adi_arg->dRho       = NULL;
        adi_arg->dRho2      = NULL;
        adi_arg->dLGMAlpha  = NULL;
        adi_arg->dLGMRho    = NULL;

        adi_arg->dPathInfos = NULL;
    }
}

/*	Find Numer TStar */
void cts_get_numer_tstar(long today, CTS cts, double* numer_tstar)
{
    long   temp_date;
    double temp_numer_tstar;

    /* Choose the new Tstar in the middle of the structure */
    if (fabs(*numer_tstar) < 1.0E-08)
    {
        temp_numer_tstar =
            0.5 * (cts->exo_leg->cpn[cts->call[0].exo_idx].cpn_start_time +
                   cts->exo_leg->cpn[cts->call[cts->num_calls - 1].exo_idx].cpn_pay_time);
        temp_date    = (long)(today + DAYS_IN_YEAR * temp_numer_tstar + 0.5);
        *numer_tstar = (temp_date - today) * YEARS_IN_DAY;
    }
}

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

/*	Fill and check all the relevant structures */
Err cts_fill_check_all_struct(
    /*	Market */
    CTS_MKT mkt,
    /*	The underlying */
    double tstar,
    int    use_calib, /*	0: use lgmsvund, 1: calibrate */
    /*		if calib */
    double  lambda, /*	LGM lambda */
    int     nb_factor,
    double  lgm_alpha,
    double  lgm_gamma,
    double  lgm_rho,
    int     nsmilepar,
    double* smilepartime,
    double* alphaeps, /*	alpha */
    double* rhoeps,   /*	rho */
    double* ldaeps,   /*	ldaeps */
    double* rho2eps,
    /*	End of calib params */
    char* lgmsvund,
    /*	The structure */
    /*		funding */
    int     fund_ncpn,
    double* fund_not,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    /*		ts */
    /*			general Libor properties */
    int fix_lag_bd,
    /*			number of coupons */
    int ncpn,
    /*			coupon description */
    /*				cpn */
    int* cpn_type, /*	0:	fixed
                                   1:	libor fixed at start
                                   2:	libor fixed at end */
    long*   cpn_start_date,
    long*   cpn_end_date,
    long*   cpn_pay_date,
    double* cpn_coupon,
    char**  cpn_basis,
    double* cpn_not,
    /*				pay libor */
    long*   pay_fix_date, /*	fixing date of the Libor */
    int     pay_months,
    char*   pay_freq,
    char*   pay_basis,
    double* pay_fwd_spread,
    double* pay_gearing,
    /*				fix libor */
    /*				fixings */
    int*     nfix,
    double** weights,
    long**   fix_dates,
    char**   fix_tenor, /*	Fixing tenors */
    char*    fix_freq,
    char*    fix_basis,
    double** fix_fwd_spreads,
    /*				profiles */
    int tstype, /*	0: generic, 1: range */
    /*				profile 0 */
    double*  alpha,
    double*  beta,
    int*     nstr,
    double** str,
    double** nbopt,
    /*				profile 1 */
    int     buy_sell,   /*	1: BNPP buys, -1: BNPP sells */
    int     value_zero, /*	0: Do not value low strike options, 1: do */
    double* lb,         /*	Lower bounds of ranges */
    double* ub,         /*	Upper bounds of ranges */
    double* payoff,     /*	Payoff in ranges */
    double  call_spread,
    double  numer_spread,
    /*		Trimming */
    int trim_type, /*	0: no trim
                                   1: x fixings max
                                   2: x time min between two fixings */
    int    max_fix,
    double min_fix_time,
    /*	Extra model parameters*/
    int    use_cmsopt,        /*	Use CmsOption to value fix option, use BS on CmsRate otherwise */
    double correl_start,      /*	Correl between libor fixing and libor paid at start */
    double correl_end,        /*	Correl between libor fixing and libor paid at start */
    int    float_adjust_type, /*	type of adjustment for the floating coupon,
                                                      0: ATM vol, 1: Strike Vol */
    /*	Calls */
    int     ncall,
    int     pay_rec, /*	1: rec pd, -1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,
    int     adj_fee,

    /*	Flag for extra calculation / adjustment of one-time callable */
    int do_one_time,    /*	1: calc the one time */
    int one_time_index, /*	0: choose automatically the index, >0: index provided by user */
    int just_recalib,   /*	Just recalibrate the underlying and skip all the other parts */

    /*	Numerical params */
    /*		CF */
    LGMSV_NUMERPARAMS NumerParams,
    /*		ADI / MC */
    int    pde_or_mc,
    int    req_stp,
    int    req_stppsi,
    int    req_stpx,
    int    req_stpz,
    double req_mintime,
    long   req_paths,
    double integ_mintime,
    /*		Calib */
    char*               cal_tenor,
    char*               cal_ref,
    char*               cal_freq,
    char*               cal_basis,
    int                 force_atm, /*	force atm calib */
    double              max_std_long,
    double              max_std_short,
    double              vol_shift_long,
    DIAGCALIB_VOLTYPE   vol_type_long,
    DIAGCALIB_SHIFTTYPE vol_shift_type_long,
    double              vol_shift_short,
    DIAGCALIB_VOLTYPE   vol_type_short,
    DIAGCALIB_SHIFTTYPE vol_shift_type_short,
    double              lambda_shift,
    int                 calib_strategy, /*	-1: autocal, 0: swaptions / cap, 1: cap / swaptions */
    int                 fix_lambda,     /*	0: calib lambda to cap, 1: fix lambda calib
                                                        to diagonal */
    char* short_tenor,
    char* short_refrate,
    char* short_freq,
    char* short_basis,
    int   fix_smile,          /*	0: calib smile parameters to market smile */
    int   smile_calib_months, /* 0: co-terminal swaption, otherwise underlyings with required nb
                                 months */
    LGMSV_CalibParams* lgmsv_calib_params,
    double             min_time,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps, /*	1: we allow jumps on vol, 0: we don't */
    double numer_tstar,
    double prec,
    int    maxiter,
    int    keep_first,
    /*	Strike choice */
    int long_strike_flag,  /*	0: ATM
                                                   1: Coupon
                                                   2: Eq (PV/Lvl) */
    int short_strike_flag, /*	0: ATM,
                                                   1: implied digital caplet strike
                                                   2: same number of std */
    CTS cts_iv,            /*	Needed to compute implied strike */
    /*	Flags */
    /*		IV calculation */
    int calc_mkt_iv,
    int calc_fwd_iv,
    /*		EOD */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	Results */
    CTS     cts,
    CTS_UND und,
    int*    call_feat, /*	0: No callable feature to be valued
                                       1: Callable feature to be valued through adi */
    CTS_ADI_ARG adi_arg,
    /*	Feedback */
    int                  save_inst_data,
    cpd_calib_inst_data* inst_data,
    int                  save_fwdiv,
    CTS_IV               fwd_iv_info)
{
    Err err = NULL;

    /*	Initialisation */
    und->model.iNbPWTime  = 0;
    und->model.dPWTime    = NULL;
    und->model.dSigma     = NULL;
    und->model.dAlpha     = NULL;
    und->model.dLambdaEps = NULL;
    und->model.dLvlEps    = NULL;
    und->model.dRho       = NULL;
    und->model.dRho2      = NULL;

    adi_arg->time     = NULL;
    adi_arg->date     = NULL;
    adi_arg->void_prm = NULL;
    adi_arg->is_event = NULL;
    adi_arg->ifr      = NULL;
    adi_arg->params   = NULL;

    if (!just_recalib)
    {
        cts->fund_leg = NULL;
        cts->exo_leg  = NULL;
        cts->call     = NULL;

        cpd_init_calib_inst_data(inst_data);

        /*	Funding leg */

        cts->fund_leg = (CTS_FUND_LEG)malloc(sizeof(cts_fund_leg));
        if (!cts->fund_leg)
        {
            err = "Memory allocation error (1) in cts_fill_check_all_struct";
            goto FREE_RETURN;
        }

        err = cts_fill_fund_leg(
            mkt,
            eod_fix_flag,
            fund_ncpn,
            fund_not,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            cts->fund_leg);

        if (err)
        {
            goto FREE_RETURN;
        }

        /*	Exotic leg */

        cts->exo_leg = (CTS_EXO_LEG)malloc(sizeof(cts_exo_leg));
        if (!cts->exo_leg)
        {
            err = "Memory allocation error (2) in cts_fill_check_all_struct";
            goto FREE_RETURN;
        }

        err = cts_fill_cpn_leg(
            mkt,
            eod_fix_flag,
            fix_lag_bd,
            ncpn,
            cpn_type,
            cpn_start_date,
            cpn_end_date,
            cpn_pay_date,
            cpn_coupon,
            cpn_basis,
            cpn_not,
            pay_fix_date,
            pay_months,
            pay_freq,
            pay_basis,
            pay_fwd_spread,
            pay_gearing,
            nfix,
            weights,
            fix_dates,
            fix_tenor,
            fix_freq,
            fix_basis,
            fix_fwd_spreads,
            tstype,
            alpha,
            beta,
            nstr,
            str,
            nbopt,
            buy_sell,
            value_zero,
            lb,
            ub,
            payoff,
            call_spread,
            numer_spread,
            trim_type,
            max_fix,
            min_fix_time,
            calc_mkt_iv,
            use_cmsopt,
            correl_start,
            correl_end,
            float_adjust_type,
            cts->exo_leg);

        if (err)
        {
            goto FREE_RETURN;
        }

        /*	Calls */

        if (ncall > 0 && ex_date[ncall - 1] >= mkt->today + eod_ex_flag)
        {
            err = cts_fill_calls(
                mkt->today, eod_ex_flag, ncall, pay_rec, ex_date, set_date, fee, cts);

            if (err)
            {
                goto FREE_RETURN;
            }

            err = cts_calc_equivalent_strikes(cts_iv, cts, mkt);

            if (err)
            {
                goto FREE_RETURN;
            }

            /* Adjustment of floating dates */
            err = cts_adjust_floating_fixing(cts, pde_or_mc);

            if (err)
                goto FREE_RETURN;
        }
        else
        {
            cts->num_calls = 0;
            cts->call      = NULL;
        }
    }

    /*	Underlying */
    if (cts->num_calls > 0 && cts->call)
    {
        if (use_calib)
        {
            err = cts_calib_und(
                mkt,
                eod_ex_flag,
                lambda,
                nb_factor,
                lgm_alpha,
                lgm_gamma,
                lgm_rho,
                nsmilepar,
                smilepartime,
                alphaeps,
                rhoeps,
                ldaeps,
                rho2eps,
                tstar,
                NumerParams,
                cal_tenor,
                cal_ref,
                cal_freq,
                cal_basis,
                force_atm,
                max_std_long,
                max_std_short,
                vol_shift_long,
                vol_type_long,
                vol_shift_type_long,
                vol_shift_short,
                vol_type_short,
                vol_shift_type_short,
                lambda_shift,
                calib_strategy,
                fix_lambda,
                short_tenor,
                short_refrate,
                short_freq,
                short_basis,
                fix_smile,
                smile_calib_months,
                lgmsv_calib_params,
                min_time,
                skip_last,
                min_fact,
                max_fact,
                use_jumps,
                prec,
                maxiter,
                keep_first,
                long_strike_flag,
                short_strike_flag,
                cts,
                und,
                save_inst_data,
                inst_data);
        }
        else
        {
            /* for now all the external underlyings are using the default TStar */
            err = cts_fill_und(mkt, lgmsvund, LGMSV_Tstar, und);
        }

        if (err)
            goto FREE_RETURN;

        /* choose the numeraire on which we diffuse */
        cts_get_numer_tstar(und->mkt->today, cts, &numer_tstar);

        Convert_Tstar_model(&(und->model), numer_tstar);
    }
    else
    {
        und->mkt = NULL;
    }

    /*	Fwd IVS for call adjustment */

    if (cts->num_calls > 0)
    {
        err = cts_calc_model_fwd_iv(und, cts, calc_fwd_iv, adj_fee, NumerParams);

        if (err)
        {
            goto FREE_RETURN;
        }
    }

    /*	ADI */
    if (cts->num_calls > 0 && cts->call)
    {
        if (do_one_time && one_time_index == 0)
        {
            /* First Calculate the index of the one time Call */
            err = cts_find_one_time_index(und, cts, &one_time_index);
        }

        err = cts_fill_algo_arg(
            und,
            cts,
            pde_or_mc,
            req_stp,
            req_stppsi,
            req_stpx,
            req_stpz,
            req_mintime,
            req_paths,
            integ_mintime,
            do_one_time,
            one_time_index,
            adi_arg);

        if (err)
        {
            goto FREE_RETURN;
        }

        *call_feat = 1;
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (err)
    {
        cts_free_all_struct(cts, und, *call_feat, adi_arg);
    }

    return err;
}

/*	Free all structures */
void cts_free_all_struct(CTS cts, CTS_UND und, int call_feat, CTS_ADI_ARG adi_arg)
{
    cts_free_und(und);
    cts_free_calls(cts);

    if (cts->fund_leg)
    {
        cts_free_fund_leg(cts->fund_leg);
        free(cts->fund_leg);
    }

    if (cts->exo_leg)
    {
        cts_free_exo_leg(cts->exo_leg);
        free(cts->exo_leg);
    }

    cts_free_adi_arg(adi_arg);
}

/*	Payoff function for adi (callable) */
/*	-----------------------------------	*/

Err cts_payoff_4_lgmsv_adi(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*  yc,
    double lam,
    double tstar,
    double alpha,
    double rho,
    double sigma,
    /* Nodes data */
    int     lpsi,
    int     upsi,
    int     lx,
    int     ux,
    int     lz,
    int     uz,
    double* psi,
    double* x,
    double* z,
    int*    nprod,
    /* Vector of results to be updated */
    double**** prod_val)
{
    CTS_PAY_ARG            cts_arg;
    CTS                    cts;
    CTS_CALL               call;
    CTS_UND                und;
    CTS_FUND_LEG           fund_leg;
    CTS_FUND_CPN           fund_cpn;
    CTS_EXO_LEG            exo_leg;
    CTS_CPN                exo_cpn;
    CTS_FIX                fix;
    CTS_EVAL_CONST_FIX     eval_const_fix;
    CTS_EVAL_CONST_FIX_CMS eval_const_fix_cms;
    CTS_EVAL_CONST_CPN     eval_const_cpn;
    CTS_EVAL_CONST_CALL    eval_const_call;
    CALIBCPNSCHEDULEDLM    schedule;

    int    i, j, k, l;
    double dUtRhoDivAlpha, tmp, tmpf;
    double DFR_start, DFR_end, DFR_fix, DFR_pay, DFR_set, coupon, libor, funding, exotic, iv, fee;
    double cms;
    int    eval_one_time, one_time_col;
    int    call_at_begining;

    Err err = NULL;

    /*	Get the event */
    cts_arg  = (CTS_PAY_ARG)func_parm;
    cts      = cts_arg->cts;
    und      = (CTS_UND)(cts_arg->und);
    fund_leg = cts->fund_leg;
    exo_leg  = cts->exo_leg;

    /*	Constant for orthogonalization */
    dUtRhoDivAlpha = sigma * rho / alpha;

    call_at_begining = 1;

    if (cts_arg->is_call)
    {
        exo_cpn = &(cts->exo_leg->cpn[cts_arg->call->exo_idx]);

        /* check */
        switch (exo_cpn->type)
        {
        case 0:
        case 3:
        case 4:

            if (exo_cpn->fix[0].ref_fix_date > evt_date)
            {
                call_at_begining = 1;
            }
            else
            {
                call_at_begining = 0;
            }

            break;

        case 1:
        case 2:

            if (exo_cpn->fix[0].ref_fix_date > evt_date && exo_cpn->pay_fix_date > evt_date)
            {
                call_at_begining = 1;
            }
            else
            {
                call_at_begining = 0;
            }

            break;
        }
    }

    /* Call Before fixings in the In arrears Time Swap case */
    if (cts_arg->is_call && call_at_begining)
    {
        eval_const_call = &(cts_arg->eval_const_call);
        call            = cts_arg->call;

        eval_one_time = eval_const_call->eval_one_time;
        *nprod += eval_one_time;

        one_time_col = *nprod - 1;

        if (call->is_midat)
        {
            for (i = lz; i <= uz; i++)
            {
                tmp = dUtRhoDivAlpha * (z[i] - 1.0);

                for (j = lx; j <= ux; j++)
                {
                    tmpf = x[j] + tmp;

                    for (k = lpsi; k <= upsi; k++)
                    {
                        DFR_set =
                            exp(-eval_const_call->tset_tstar_alpha -
                                eval_const_call->tset_tstar_beta * tmpf -
                                eval_const_call->tset_tstar_gamma * psi[k]);

                        /* first evaluate the funding */
                        funding  = 0.0;
                        fund_cpn = &(fund_leg->cpn[call->fund_idx]);

                        for (l = 0; l < call->num_fund_cpn; l++)
                        {
                            DFR_pay =
                                exp(-eval_const_call->tpay_tstar_alpha[l] -
                                    eval_const_call->tpay_tstar_beta[l] * tmpf -
                                    eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                            funding += fund_cpn->cpn_plus_ex * DFR_pay;
                            fund_cpn++;
                        }

                        /* first coupon */
                        DFR_pay =
                            exp(-eval_const_call->tpay_tstar_alpha[l] -
                                eval_const_call->tpay_tstar_beta[l] * tmpf -
                                eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                        funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

                        /* then evaluate the fixed leg */
                        exotic = 0.0;

                        for (l = 0; l < call->num_exo_cpn; l++)
                        {
                            DFR_pay =
                                exp(-eval_const_call->tpay_tstar_alpha2[l] -
                                    eval_const_call->tpay_tstar_beta2[l] * tmpf -
                                    eval_const_call->tpay_tstar_gamma2[l] * psi[k]);
                            exotic += DFR_pay;
                        }

                        iv = exotic - funding;

                        /* then calculate the fee */
                        fee = DFR_set * call->total_fee;

                        prod_val[i][j][k][0] = max(prod_val[i][j][k][0], call->pay_rec * iv - fee);

                        /* Value the one time call if required */
                        if (eval_one_time)
                        {
                            prod_val[i][j][k][one_time_col] = max(0.0, call->pay_rec * iv - fee);
                        }
                    }
                }
            }
        }
        else
        {
            for (i = lz; i <= uz; i++)
            {
                tmp = dUtRhoDivAlpha * (z[i] - 1.0);

                for (j = lx; j <= ux; j++)
                {
                    tmpf = x[j] + tmp;

                    for (k = lpsi; k <= upsi; k++)
                    {
                        DFR_set =
                            exp(-eval_const_call->tset_tstar_alpha -
                                eval_const_call->tset_tstar_beta * tmpf -
                                eval_const_call->tset_tstar_gamma * psi[k]);

                        /* first evaluate the funding */
                        funding  = 0.0;
                        fund_cpn = &(fund_leg->cpn[call->fund_idx]);

                        for (l = 0; l < call->num_fund_cpn; l++)
                        {
                            DFR_pay =
                                exp(-eval_const_call->tpay_tstar_alpha[l] -
                                    eval_const_call->tpay_tstar_beta[l] * tmpf -
                                    eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                            funding += fund_cpn->cpn_plus_ex * DFR_pay;
                            fund_cpn++;
                        }

                        /* first coupon */
                        DFR_pay =
                            exp(-eval_const_call->tpay_tstar_alpha[l] -
                                eval_const_call->tpay_tstar_beta[l] * tmpf -
                                eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                        funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

                        /* then calculate the fee */
                        fee = DFR_set * call->total_fee;

                        prod_val[i][j][k][0] =
                            max(prod_val[i][j][k][0],
                                call->pay_rec * (prod_val[i][j][k][1] - funding) - fee);

                        /* Value the one time call if required */
                        if (eval_one_time)
                        {
                            prod_val[i][j][k][one_time_col] =
                                max(0.0, call->pay_rec * (prod_val[i][j][k][1] - funding) - fee);
                        }
                    }
                }
            }
        }
    }

    /* first the fixing case */
    if (cts_arg->is_fixing)
    {
        eval_const_fix = &(cts_arg->eval_const_fix);
        fix            = cts_arg->fix;
        exo_cpn        = cts_arg->cpn;

        for (i = lz; i <= uz; i++)
        {
            tmp = dUtRhoDivAlpha * (z[i] - 1.0);

            for (j = lx; j <= ux; j++)
            {
                tmpf = x[j] + tmp;

                for (k = lpsi; k <= upsi; k++)
                {
                    DFR_fix =
                        exp(-eval_const_fix->ts_te_alpha - eval_const_fix->ts_te_beta * tmpf -
                            eval_const_fix->ts_te_gamma * psi[k]);
                    DFR_pay = exp(
                        -eval_const_fix->tpay_tstar_alpha - eval_const_fix->tpay_tstar_beta * tmpf -
                        eval_const_fix->tpay_tstar_gamma * psi[k]);

                    coupon = eval_const_fix->a + eval_const_fix->b * DFR_fix;

                    l = 0;

                    while (l < exo_cpn->nstr && DFR_fix > eval_const_fix->s[l])
                    {
                        coupon += eval_const_fix->n[l] * (DFR_fix - eval_const_fix->s[l]);
                        l++;
                    }

                    switch (exo_cpn->type)
                    {
                    case 0:
                        prod_val[i][j][k][1] += coupon * DFR_pay;
                        break;
                    case 1:
                        prod_val[i][j][k][2] += coupon * DFR_pay;
                        break;
                    case 2:
                        prod_val[i][j][k][1] += coupon * prod_val[i][j][k][2];
                        break;
                    case 3:
                        prod_val[i][j][k][1] += coupon * DFR_pay;
                        break;
                    case 4:
                        prod_val[i][j][k][1] += coupon * DFR_pay;
                        break;
                    }
                }
            }
        }
    }

    /* or the CMS fixing case */
    if (cts_arg->is_fixing_cms)
    {
        eval_const_fix_cms = &(cts_arg->eval_const_fix_cms);
        fix                = cts_arg->fix;
        exo_cpn            = cts_arg->cpn;
        schedule           = &(fix->schedule);

        for (i = lz; i <= uz; i++)
        {
            tmp = dUtRhoDivAlpha * (z[i] - 1.0);

            for (j = lx; j <= ux; j++)
            {
                tmpf = x[j] + tmp;

                for (k = lpsi; k <= upsi; k++)
                {
                    /* First DF */
                    if (schedule->iNCpn - 1 > 0)
                    {
                        DFR_start =
                            exp(-eval_const_fix_cms->ti_tstar_alpha[0] -
                                eval_const_fix_cms->ti_tstar_beta[0] * tmpf -
                                eval_const_fix_cms->ti_tstar_gamma[0] * psi[k]);
                        cms = 0.0;

                        /* Intermediary DF, coverage already included in alpha */
                        for (l = 1; l < schedule->iNCpn - 1; l++)
                        {
                            DFR_fix =
                                exp(-eval_const_fix_cms->ti_tstar_alpha[l] -
                                    eval_const_fix_cms->ti_tstar_beta[l] * tmpf -
                                    eval_const_fix_cms->ti_tstar_gamma[l] * psi[k]);
                            cms += DFR_fix;
                        }

                        /* Last DF */
                        DFR_end =
                            exp(-eval_const_fix_cms->ti_tstar_alpha[schedule->iNCpn - 1] -
                                eval_const_fix_cms->ti_tstar_beta[schedule->iNCpn - 1] * tmpf -
                                eval_const_fix_cms->ti_tstar_gamma[schedule->iNCpn - 1] * psi[k]);
                        cms += DFR_end * schedule->dCpnCvg[schedule->iNCpn - 1];

                        cms = (DFR_start - DFR_end) / cms;

                        coupon = eval_const_fix_cms->a + eval_const_fix_cms->b * cms;

                        l = 0;

                        while (l < exo_cpn->nstr && cms > eval_const_fix_cms->s[l])
                        {
                            coupon += eval_const_fix_cms->n[l] * (cms - eval_const_fix_cms->s[l]);
                            l++;
                        }
                    }
                    else
                    {
                        coupon = eval_const_fix_cms->a;
                    }

                    DFR_pay =
                        exp(-eval_const_fix_cms->tpay_tstar_alpha -
                            eval_const_fix_cms->tpay_tstar_beta * tmpf -
                            eval_const_fix_cms->tpay_tstar_gamma * psi[k]);

                    switch (exo_cpn->type)
                    {
                    case 0:
                        prod_val[i][j][k][1] += coupon * DFR_pay;
                        break;
                    case 1:
                        prod_val[i][j][k][2] += coupon * DFR_pay;
                        break;
                    case 2:
                        prod_val[i][j][k][1] += coupon * prod_val[i][j][k][2];
                        break;
                    case 3:
                        prod_val[i][j][k][1] += coupon * DFR_pay;
                        break;
                    case 4:
                        prod_val[i][j][k][1] += coupon * DFR_pay;
                        break;
                    }
                }
            }
        }
    }

    /* the coupon */
    if (cts_arg->is_cpn)
    {
        eval_const_cpn = &(cts_arg->eval_const_cpn);
        exo_cpn        = cts_arg->cpn;

        for (i = lz; i <= uz; i++)
        {
            tmp = dUtRhoDivAlpha * (z[i] - 1.0);

            for (j = lx; j <= ux; j++)
            {
                tmpf = x[j] + tmp;

                for (k = lpsi; k <= upsi; k++)
                {
                    DFR_fix =
                        exp(-eval_const_cpn->ts_te_alpha - eval_const_cpn->ts_te_beta * tmpf -
                            eval_const_cpn->ts_te_gamma * psi[k]);

                    libor = eval_const_cpn->multip * DFR_fix + eval_const_cpn->shift;

                    switch (exo_cpn->type)
                    {
                    case 1:
                        prod_val[i][j][k][1] += prod_val[i][j][k][2] * libor;
                        prod_val[i][j][k][2] = 0.0;
                        break;
                    case 2:
                        DFR_pay =
                            exp(-eval_const_cpn->tpay_tstar_alpha -
                                eval_const_cpn->tpay_tstar_beta * tmpf -
                                eval_const_cpn->tpay_tstar_gamma * psi[k]);
                        prod_val[i][j][k][2] = libor * DFR_pay;
                        break;
                    }
                }
            }
        }
    }

    /* Eventually call the coupon except in the CTS Floating in arrears */
    if (cts_arg->is_call && !call_at_begining)
    {
        eval_const_call = &(cts_arg->eval_const_call);
        call            = cts_arg->call;

        eval_one_time = eval_const_call->eval_one_time;
        *nprod += eval_one_time;

        one_time_col = *nprod - 1;

        if (call->is_midat)
        {
            for (i = lz; i <= uz; i++)
            {
                tmp = dUtRhoDivAlpha * (z[i] - 1.0);

                for (j = lx; j <= ux; j++)
                {
                    tmpf = x[j] + tmp;

                    for (k = lpsi; k <= upsi; k++)
                    {
                        DFR_set =
                            exp(-eval_const_call->tset_tstar_alpha -
                                eval_const_call->tset_tstar_beta * tmpf -
                                eval_const_call->tset_tstar_gamma * psi[k]);

                        /* first evaluate the funding */
                        funding  = 0.0;
                        fund_cpn = &(fund_leg->cpn[call->fund_idx]);

                        for (l = 0; l < call->num_fund_cpn; l++)
                        {
                            DFR_pay =
                                exp(-eval_const_call->tpay_tstar_alpha[l] -
                                    eval_const_call->tpay_tstar_beta[l] * tmpf -
                                    eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                            funding += fund_cpn->cpn_plus_ex * DFR_pay;
                            fund_cpn++;
                        }

                        /* first coupon */
                        DFR_pay =
                            exp(-eval_const_call->tpay_tstar_alpha[l] -
                                eval_const_call->tpay_tstar_beta[l] * tmpf -
                                eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                        funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

                        /* then evaluate the fixed leg */
                        exotic = 0.0;

                        for (l = 0; l < call->num_exo_cpn; l++)
                        {
                            DFR_pay =
                                exp(-eval_const_call->tpay_tstar_alpha2[l] -
                                    eval_const_call->tpay_tstar_beta2[l] * tmpf -
                                    eval_const_call->tpay_tstar_gamma2[l] * psi[k]);
                            exotic += DFR_pay;
                        }

                        iv = exotic - funding;

                        /* then calculate the fee */
                        fee = DFR_set * call->total_fee;

                        prod_val[i][j][k][0] = max(prod_val[i][j][k][0], call->pay_rec * iv - fee);

                        /* Value the one time call if required */
                        if (eval_one_time)
                        {
                            prod_val[i][j][k][one_time_col] = max(0.0, call->pay_rec * iv - fee);
                        }
                    }
                }
            }
        }
        else
        {
            for (i = lz; i <= uz; i++)
            {
                tmp = dUtRhoDivAlpha * (z[i] - 1.0);

                for (j = lx; j <= ux; j++)
                {
                    tmpf = x[j] + tmp;

                    for (k = lpsi; k <= upsi; k++)
                    {
                        DFR_set =
                            exp(-eval_const_call->tset_tstar_alpha -
                                eval_const_call->tset_tstar_beta * tmpf -
                                eval_const_call->tset_tstar_gamma * psi[k]);

                        /* first evaluate the funding */
                        funding  = 0.0;
                        fund_cpn = &(fund_leg->cpn[call->fund_idx]);

                        for (l = 0; l < call->num_fund_cpn; l++)
                        {
                            DFR_pay =
                                exp(-eval_const_call->tpay_tstar_alpha[l] -
                                    eval_const_call->tpay_tstar_beta[l] * tmpf -
                                    eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                            funding += fund_cpn->cpn_plus_ex * DFR_pay;
                            fund_cpn++;
                        }

                        /* first coupon */
                        DFR_pay =
                            exp(-eval_const_call->tpay_tstar_alpha[l] -
                                eval_const_call->tpay_tstar_beta[l] * tmpf -
                                eval_const_call->tpay_tstar_gamma[l] * psi[k]);
                        funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

                        /* then calculate the fee */
                        fee = DFR_set * call->total_fee;

                        prod_val[i][j][k][0] =
                            max(prod_val[i][j][k][0],
                                call->pay_rec * (prod_val[i][j][k][1] - funding) - fee);

                        /* Value the one time call if required */
                        if (eval_one_time)
                        {
                            prod_val[i][j][k][one_time_col] =
                                max(0.0, call->pay_rec * (prod_val[i][j][k][1] - funding) - fee);
                        }
                    }
                }
            }
        }
    }

    /*	End of payoff valuation
    ----------------------- */

    return err;
}

/*	Payoff function for MC (callable with optimise KO) */
/*	-------------------------------------------------- */

Err cts_payoff_4_lgmsv_mc(
    /* Event */
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,
    double ft,
    double psi,
    double v,
    int    nprod,
    /* Vector of results to be updated */
    double* prod_val,
    int*    stop_path)
{
    CTS_PAY_ARG            cts_arg;
    CTS                    cts;
    CTS_CALL               call;
    CTS_UND                und;
    CTS_FUND_LEG           fund_leg;
    CTS_FUND_CPN           fund_cpn;
    CTS_EXO_LEG            exo_leg;
    CTS_CPN                exo_cpn;
    CTS_FIX                fix;
    CTS_EVAL_CONST_FIX     eval_const_fix;
    CTS_EVAL_CONST_FIX_CMS eval_const_fix_cms;
    CTS_EVAL_CONST_CPN     eval_const_cpn;
    CTS_EVAL_CONST_CALL    eval_const_call;
    CALIBCPNSCHEDULEDLM    schedule;

    int    l;
    double DFR_start, DFR_end, DFR_fix, DFR_pay, DFR_set, coupon, libor, partial_funding, funding,
        exotic, iv, fee;
    double cms;

    Err err = NULL;

    /*	Get the event */
    cts_arg  = (CTS_PAY_ARG)func_parm;
    cts      = cts_arg->cts;
    und      = (CTS_UND)(cts_arg->und);
    fund_leg = cts->fund_leg;
    exo_leg  = cts->exo_leg;

    memset(prod_val, 0, nprod * sizeof(double));

    /* Call Before fixings in the In arrears Time Swap case */
    if (cts_arg->is_call)
    {
        eval_const_call = &(cts_arg->eval_const_call);
        call            = cts_arg->call;

        if (call->is_midat)
        {
            DFR_set =
                exp(-eval_const_call->tset_tstar_alpha - eval_const_call->tset_tstar_beta * ft -
                    eval_const_call->tset_tstar_gamma * psi);

            /* first evaluate the funding */
            funding         = 0.0;
            partial_funding = 0.0;
            fund_cpn        = &(fund_leg->cpn[call->fund_idx]);

            for (l = 0; l < call->num_fund_cpn; l++)
            {
                DFR_pay =
                    exp(-eval_const_call->tpay_tstar_alpha[l] -
                        eval_const_call->tpay_tstar_beta[l] * ft -
                        eval_const_call->tpay_tstar_gamma[l] * psi);
                funding += fund_cpn->cpn_plus_ex * DFR_pay;

                if (l < call->num_partial_fund_cpn)
                {
                    if (l < call->num_partial_fund_cpn - 1)
                    {
                        partial_funding += fund_cpn->cpn_plus_ex * DFR_pay;
                    }
                    else
                    {
                        partial_funding += fund_cpn->cpn_plus_ex_partial * DFR_pay;
                    }
                }

                fund_cpn++;
            }

            /* first coupon */
            DFR_pay = exp(
                -eval_const_call->tpay_tstar_alpha[l] - eval_const_call->tpay_tstar_beta[l] * ft -
                eval_const_call->tpay_tstar_gamma[l] * psi);
            funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;
            partial_funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

            /* then evaluate the fixed leg */
            exotic = 0.0;

            for (l = 0; l < call->num_partial_exo_cpn; l++)
            {
                DFR_pay =
                    exp(-eval_const_call->tpay_tstar_alpha2[l] -
                        eval_const_call->tpay_tstar_beta2[l] * ft -
                        eval_const_call->tpay_tstar_gamma2[l] * psi);
                exotic += DFR_pay;
            }

            iv = exotic - partial_funding;

            /* then calculate the fee */
            fee = DFR_set * call->total_fee;

            prod_val[0] += -call->pay_rec * iv - fee;
            prod_val[1] = -call->pay_rec * (-funding);
            prod_val[2] = v;
        }
        else
        {
            DFR_set =
                exp(-eval_const_call->tset_tstar_alpha - eval_const_call->tset_tstar_beta * ft -
                    eval_const_call->tset_tstar_gamma * psi);

            /* first evaluate the funding */
            funding         = 0.0;
            partial_funding = 0.0;
            fund_cpn        = &(fund_leg->cpn[call->fund_idx]);

            for (l = 0; l < call->num_fund_cpn; l++)
            {
                DFR_pay =
                    exp(-eval_const_call->tpay_tstar_alpha[l] -
                        eval_const_call->tpay_tstar_beta[l] * ft -
                        eval_const_call->tpay_tstar_gamma[l] * psi);
                funding += fund_cpn->cpn_plus_ex * DFR_pay;

                if (l < call->num_partial_fund_cpn)
                {
                    if (l < call->num_partial_fund_cpn - 1)
                    {
                        partial_funding += fund_cpn->cpn_plus_ex * DFR_pay;
                    }
                    else
                    {
                        partial_funding += fund_cpn->cpn_plus_ex_partial * DFR_pay;
                    }
                }

                fund_cpn++;
            }

            /* first coupon */
            DFR_pay = exp(
                -eval_const_call->tpay_tstar_alpha[l] - eval_const_call->tpay_tstar_beta[l] * ft -
                eval_const_call->tpay_tstar_gamma[l] * psi);
            funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;
            partial_funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

            /* then calculate the fee */
            fee = DFR_set * call->total_fee;

            prod_val[0] += -call->pay_rec * (-partial_funding) + fee;
            prod_val[1] = -call->pay_rec * (-funding);
            prod_val[2] = v;
        }
    }

    /* the coupon */
    if (cts_arg->is_cpn)
    {
        eval_const_cpn = &(cts_arg->eval_const_cpn);
        exo_cpn        = cts_arg->cpn;

        DFR_fix =
            exp(-eval_const_cpn->ts_te_alpha - eval_const_cpn->ts_te_beta * ft -
                eval_const_cpn->ts_te_gamma * psi);

        libor = eval_const_cpn->multip * DFR_fix + eval_const_cpn->shift;

        switch (exo_cpn->type)
        {
        case 1:
            cts_arg->dPathInfos[path_index] = libor;
            prod_val[0] += 0.0;
            break;
        case 2:
            DFR_pay =
                exp(-eval_const_cpn->tpay_tstar_alpha - eval_const_cpn->tpay_tstar_beta * ft -
                    eval_const_cpn->tpay_tstar_gamma * psi);
            prod_val[0] += -cts->call->pay_rec * libor * DFR_pay * cts_arg->dPathInfos[path_index];
            cts_arg->dPathInfos[path_index] = 0.0;
            break;
        }
    }

    /* first the fixing case */
    if (cts_arg->is_fixing)
    {
        eval_const_fix = &(cts_arg->eval_const_fix);
        fix            = cts_arg->fix;
        exo_cpn        = cts_arg->cpn;

        DFR_fix =
            exp(-eval_const_fix->ts_te_alpha - eval_const_fix->ts_te_beta * ft -
                eval_const_fix->ts_te_gamma * psi);
        DFR_pay =
            exp(-eval_const_fix->tpay_tstar_alpha - eval_const_fix->tpay_tstar_beta * ft -
                eval_const_fix->tpay_tstar_gamma * psi);

        coupon = eval_const_fix->a + eval_const_fix->b * DFR_fix;

        l = 0;

        while (l < exo_cpn->nstr && DFR_fix > eval_const_fix->s[l])
        {
            coupon += eval_const_fix->n[l] * (DFR_fix - eval_const_fix->s[l]);
            l++;
        }

        switch (exo_cpn->type)
        {
        case 0:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 1:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay * cts_arg->dPathInfos[path_index];
            break;
        case 2:
            cts_arg->dPathInfos[path_index] += coupon;
            prod_val[0] += 0.0;
            break;
        case 3:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 4:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        }
    }

    /* or the CMS fixing case */
    if (cts_arg->is_fixing_cms)
    {
        eval_const_fix_cms = &(cts_arg->eval_const_fix_cms);
        fix                = cts_arg->fix;
        exo_cpn            = cts_arg->cpn;
        schedule           = &(fix->schedule);

        /* First DF */
        if (schedule->iNCpn - 1 > 0)
        {
            DFR_start = exp(
                -eval_const_fix_cms->ti_tstar_alpha[0] - eval_const_fix_cms->ti_tstar_beta[0] * ft -
                eval_const_fix_cms->ti_tstar_gamma[0] * psi);
            cms = 0.0;

            /* Intermediary DF, coverage already included in alpha */
            for (l = 1; l < schedule->iNCpn - 1; l++)
            {
                DFR_fix =
                    exp(-eval_const_fix_cms->ti_tstar_alpha[l] -
                        eval_const_fix_cms->ti_tstar_beta[l] * ft -
                        eval_const_fix_cms->ti_tstar_gamma[l] * psi);
                cms += DFR_fix;
            }

            /* Last DF */
            DFR_end =
                exp(-eval_const_fix_cms->ti_tstar_alpha[schedule->iNCpn - 1] -
                    eval_const_fix_cms->ti_tstar_beta[schedule->iNCpn - 1] * ft -
                    eval_const_fix_cms->ti_tstar_gamma[schedule->iNCpn - 1] * psi);
            cms += DFR_end * schedule->dCpnCvg[schedule->iNCpn - 1];

            cms = (DFR_start - DFR_end) / cms;

            coupon = eval_const_fix_cms->a + eval_const_fix_cms->b * cms;

            l = 0;

            while (l < exo_cpn->nstr && cms > eval_const_fix_cms->s[l])
            {
                coupon += eval_const_fix_cms->n[l] * (cms - eval_const_fix_cms->s[l]);
                l++;
            }
        }
        else
        {
            coupon = eval_const_fix_cms->a;
        }

        DFR_pay =
            exp(-eval_const_fix_cms->tpay_tstar_alpha - eval_const_fix_cms->tpay_tstar_beta * ft -
                eval_const_fix_cms->tpay_tstar_gamma * psi);

        switch (exo_cpn->type)
        {
        case 0:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 1:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay * cts_arg->dPathInfos[path_index];
            break;
        case 2:
            cts_arg->dPathInfos[path_index] += coupon;
            prod_val[0] += 0.0;
            break;
        case 3:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 4:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        }
    }

    /*	End of payoff valuation
    ----------------------- */

    return err;
}

Err cts_payoff_4_lgmsv2F_mc(
    /* Event */
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,
    double ft1,
    double ft2,
    double psi1,
    double psi2,
    double psi12,
    double v,
    int    nprod,
    /* Vector of results to be updated */
    double* prod_val,
    int*    stop_path)
{
    CTS_PAY_ARG            cts_arg;
    CTS                    cts;
    CTS_CALL               call;
    CTS_UND                und;
    CTS_FUND_LEG           fund_leg;
    CTS_FUND_CPN           fund_cpn;
    CTS_EXO_LEG            exo_leg;
    CTS_CPN                exo_cpn;
    CTS_FIX                fix;
    CTS_EVAL_CONST_FIX     eval_const_fix;
    CTS_EVAL_CONST_FIX_CMS eval_const_fix_cms;
    CTS_EVAL_CONST_CPN     eval_const_cpn;
    CTS_EVAL_CONST_CALL    eval_const_call;
    CALIBCPNSCHEDULEDLM    schedule;

    int    l;
    double DFR_start, DFR_end, DFR_fix, DFR_pay, DFR_set, coupon, libor, partial_funding, funding,
        exotic, iv, fee;
    double cms;

    Err err = NULL;

    /*	Get the event */
    cts_arg  = (CTS_PAY_ARG)func_parm;
    cts      = cts_arg->cts;
    und      = (CTS_UND)(cts_arg->und);
    fund_leg = cts->fund_leg;
    exo_leg  = cts->exo_leg;

    memset(prod_val, 0, nprod * sizeof(double));

    /* Call Before fixings in the In arrears Time Swap case */
    if (cts_arg->is_call)
    {
        eval_const_call = &(cts_arg->eval_const_call);
        call            = cts_arg->call;

        if (call->is_midat)
        {
            DFR_set = exp(
                -eval_const_call->tset_tstar_alpha - eval_const_call->tset_tstar_beta * ft1 -
                eval_const_call->tset_tstar_beta2 * ft2 - eval_const_call->tset_tstar_gamma * psi1 -
                eval_const_call->tset_tstar_gamma2 * psi2 -
                eval_const_call->tset_tstar_gamma12 * psi12);

            /* first evaluate the funding */
            funding         = 0.0;
            partial_funding = 0.0;
            fund_cpn        = &(fund_leg->cpn[call->fund_idx]);

            for (l = 0; l < call->num_fund_cpn; l++)
            {
                DFR_pay =
                    exp(-eval_const_call->tpay_tstar_alpha[l] -
                        eval_const_call->tpay_tstar_beta[l] * ft1 -
                        eval_const_call->tpay_tstar_beta_2[l] * ft2 -
                        eval_const_call->tpay_tstar_gamma[l] * psi1 -
                        eval_const_call->tpay_tstar_gamma_2[l] * psi2 -
                        eval_const_call->tpay_tstar_gamma_12[l] * psi12);
                funding += fund_cpn->cpn_plus_ex * DFR_pay;

                if (l < call->num_partial_fund_cpn)
                {
                    if (l < call->num_partial_fund_cpn - 1)
                    {
                        partial_funding += fund_cpn->cpn_plus_ex * DFR_pay;
                    }
                    else
                    {
                        partial_funding += fund_cpn->cpn_plus_ex_partial * DFR_pay;
                    }
                }

                fund_cpn++;
            }

            /* first coupon */
            DFR_pay = exp(
                -eval_const_call->tpay_tstar_alpha[l] - eval_const_call->tpay_tstar_beta[l] * ft1 -
                eval_const_call->tpay_tstar_beta_2[l] * ft2 -
                eval_const_call->tpay_tstar_gamma[l] * psi1 -
                eval_const_call->tpay_tstar_gamma_2[l] * psi2 -
                eval_const_call->tpay_tstar_gamma_12[l] * psi12);
            funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;
            partial_funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

            /* then evaluate the fixed leg */
            exotic = 0.0;

            for (l = 0; l < call->num_partial_exo_cpn; l++)
            {
                DFR_pay =
                    exp(-eval_const_call->tpay_tstar_alpha2[l] -
                        eval_const_call->tpay_tstar_beta2[l] * ft1 -
                        eval_const_call->tpay_tstar_beta2_2[l] * ft2 -
                        eval_const_call->tpay_tstar_gamma2[l] * psi1 -
                        eval_const_call->tpay_tstar_gamma2_2[l] * psi2 -
                        eval_const_call->tpay_tstar_gamma2_12[l] * psi12);
                exotic += DFR_pay;
            }

            iv = exotic - partial_funding;

            /* then calculate the fee */
            fee = DFR_set * call->total_fee;

            prod_val[0] += -call->pay_rec * iv - fee;
            prod_val[1] = -call->pay_rec * (-funding);
            prod_val[2] = v;
            prod_val[3] = -call->pay_rec * (-partial_funding);
        }
        else
        {
            DFR_set = exp(
                -eval_const_call->tset_tstar_alpha - eval_const_call->tset_tstar_beta * ft1 -
                eval_const_call->tset_tstar_beta2 * ft2 - eval_const_call->tset_tstar_gamma * psi1 -
                eval_const_call->tset_tstar_gamma2 * psi2 -
                eval_const_call->tset_tstar_gamma12 * psi12);

            /* first evaluate the funding */
            funding         = 0.0;
            partial_funding = 0.0;
            fund_cpn        = &(fund_leg->cpn[call->fund_idx]);

            for (l = 0; l < call->num_fund_cpn; l++)
            {
                DFR_pay =
                    exp(-eval_const_call->tpay_tstar_alpha[l] -
                        eval_const_call->tpay_tstar_beta[l] * ft1 -
                        eval_const_call->tpay_tstar_beta_2[l] * ft2 -
                        eval_const_call->tpay_tstar_gamma[l] * psi1 -
                        eval_const_call->tpay_tstar_gamma_2[l] * psi2 -
                        eval_const_call->tpay_tstar_gamma_12[l] * psi12);
                funding += fund_cpn->cpn_plus_ex * DFR_pay;

                if (l < call->num_partial_fund_cpn)
                {
                    if (l < call->num_partial_fund_cpn - 1)
                    {
                        partial_funding += fund_cpn->cpn_plus_ex * DFR_pay;
                    }
                    else
                    {
                        partial_funding += fund_cpn->cpn_plus_ex_partial * DFR_pay;
                    }
                }

                fund_cpn++;
            }

            /* first coupon */
            DFR_pay = exp(
                -eval_const_call->tpay_tstar_alpha[l] - eval_const_call->tpay_tstar_beta[l] * ft1 -
                eval_const_call->tpay_tstar_beta_2[l] * ft2 -
                eval_const_call->tpay_tstar_gamma[l] * psi1 -
                eval_const_call->tpay_tstar_gamma_2[l] * psi2 -
                eval_const_call->tpay_tstar_gamma_12[l] * psi12);
            funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;
            partial_funding += fund_leg->cpn[call->fund_idx].not *DFR_pay;

            /* then calculate the fee */
            fee = DFR_set * call->total_fee;

            prod_val[0] += -call->pay_rec * (-partial_funding) + fee;
            prod_val[1] = -call->pay_rec * (-funding);
            prod_val[2] = v;
            prod_val[3] = -call->pay_rec * (-partial_funding);
        }
    }

    /* the coupon */
    if (cts_arg->is_cpn)
    {
        eval_const_cpn = &(cts_arg->eval_const_cpn);
        exo_cpn        = cts_arg->cpn;

        DFR_fix =
            exp(-eval_const_cpn->ts_te_alpha - eval_const_cpn->ts_te_beta * ft1 -
                eval_const_cpn->ts_te_beta2 * ft2 - eval_const_cpn->ts_te_gamma * psi1 -
                eval_const_cpn->ts_te_gamma2 * psi2 - eval_const_cpn->ts_te_gamma12 * psi12);

        libor = eval_const_cpn->multip * DFR_fix + eval_const_cpn->shift;

        switch (exo_cpn->type)
        {
        case 1:
            cts_arg->dPathInfos[path_index] = libor;
            prod_val[0] += 0.0;
            break;
        case 2:
            DFR_pay = exp(
                -eval_const_cpn->tpay_tstar_alpha - eval_const_cpn->tpay_tstar_beta * ft1 -
                eval_const_cpn->tpay_tstar_beta2 * ft2 - eval_const_cpn->tpay_tstar_gamma * psi1 -
                eval_const_cpn->tpay_tstar_gamma2 * psi2 -
                eval_const_cpn->tpay_tstar_gamma12 * psi12);
            prod_val[0] += -cts->call->pay_rec * libor * DFR_pay * cts_arg->dPathInfos[path_index];
            cts_arg->dPathInfos[path_index] = 0.0;
            break;
        }
    }

    /* first the fixing case */
    if (cts_arg->is_fixing)
    {
        eval_const_fix = &(cts_arg->eval_const_fix);
        fix            = cts_arg->fix;
        exo_cpn        = cts_arg->cpn;

        DFR_fix =
            exp(-eval_const_fix->ts_te_alpha - eval_const_fix->ts_te_beta * ft1 -
                eval_const_fix->ts_te_beta2 * ft2 - eval_const_fix->ts_te_gamma * psi1 -
                eval_const_fix->ts_te_gamma2 * psi2 - eval_const_fix->ts_te_gamma12 * psi12);
        DFR_pay = exp(
            -eval_const_fix->tpay_tstar_alpha - eval_const_fix->tpay_tstar_beta * ft1 -
            eval_const_fix->tpay_tstar_beta2 * ft2 - eval_const_fix->tpay_tstar_gamma * psi1 -
            eval_const_fix->tpay_tstar_gamma2 * psi2 - eval_const_fix->tpay_tstar_gamma12 * psi12);

        coupon = eval_const_fix->a + eval_const_fix->b * DFR_fix;

        l = 0;

        while (l < exo_cpn->nstr && DFR_fix > eval_const_fix->s[l])
        {
            coupon += eval_const_fix->n[l] * (DFR_fix - eval_const_fix->s[l]);
            l++;
        }

        switch (exo_cpn->type)
        {
        case 0:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 1:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay * cts_arg->dPathInfos[path_index];
            break;
        case 2:
            cts_arg->dPathInfos[path_index] += coupon;
            prod_val[0] += 0.0;
            break;
        case 3:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 4:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        }
    }

    /* or the CMS fixing case */
    if (cts_arg->is_fixing_cms)
    {
        eval_const_fix_cms = &(cts_arg->eval_const_fix_cms);
        fix                = cts_arg->fix;
        exo_cpn            = cts_arg->cpn;
        schedule           = &(fix->schedule);

        /* First DF */
        if (schedule->iNCpn - 1 > 0)
        {
            DFR_start =
                exp(-eval_const_fix_cms->ti_tstar_alpha[0] -
                    eval_const_fix_cms->ti_tstar_beta[0] * ft1 -
                    eval_const_fix_cms->ti_tstar_beta2[0] * ft2 -
                    eval_const_fix_cms->ti_tstar_gamma[0] * psi1 -
                    eval_const_fix_cms->ti_tstar_gamma2[0] * psi2 -
                    eval_const_fix_cms->ti_tstar_gamma12[0] * psi12);
            cms = 0.0;

            /* Intermediary DF, coverage already included in alpha */
            for (l = 1; l < schedule->iNCpn - 1; l++)
            {
                DFR_fix =
                    exp(-eval_const_fix_cms->ti_tstar_alpha[l] -
                        eval_const_fix_cms->ti_tstar_beta[l] * ft1 -
                        eval_const_fix_cms->ti_tstar_beta2[l] * ft2 -
                        eval_const_fix_cms->ti_tstar_gamma[l] * psi1 -
                        eval_const_fix_cms->ti_tstar_gamma2[l] * psi2 -
                        eval_const_fix_cms->ti_tstar_gamma12[l] * psi12);
                cms += DFR_fix;
            }

            /* Last DF */
            DFR_end =
                exp(-eval_const_fix_cms->ti_tstar_alpha[schedule->iNCpn - 1] -
                    eval_const_fix_cms->ti_tstar_beta[schedule->iNCpn - 1] * ft1 -
                    eval_const_fix_cms->ti_tstar_beta2[schedule->iNCpn - 1] * ft2 -
                    eval_const_fix_cms->ti_tstar_gamma[schedule->iNCpn - 1] * psi1 -
                    eval_const_fix_cms->ti_tstar_gamma2[schedule->iNCpn - 1] * psi2 -
                    eval_const_fix_cms->ti_tstar_gamma12[schedule->iNCpn - 1] * psi12);
            cms += DFR_end * schedule->dCpnCvg[schedule->iNCpn - 1];

            cms = (DFR_start - DFR_end) / cms;

            coupon = eval_const_fix_cms->a + eval_const_fix_cms->b * cms;

            l = 0;

            while (l < exo_cpn->nstr && cms > eval_const_fix_cms->s[l])
            {
                coupon += eval_const_fix_cms->n[l] * (cms - eval_const_fix_cms->s[l]);
                l++;
            }
        }
        else
        {
            coupon = eval_const_fix_cms->a;
        }

        DFR_pay =
            exp(-eval_const_fix_cms->tpay_tstar_alpha - eval_const_fix_cms->tpay_tstar_beta * ft1 -
                eval_const_fix_cms->tpay_tstar_beta2 * ft2 -
                eval_const_fix_cms->tpay_tstar_gamma * psi1 -
                eval_const_fix_cms->tpay_tstar_gamma2 * psi2 -
                eval_const_fix_cms->tpay_tstar_gamma12 * psi12);

        switch (exo_cpn->type)
        {
        case 0:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 1:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay * cts_arg->dPathInfos[path_index];
            break;
        case 2:
            cts_arg->dPathInfos[path_index] += coupon;
            prod_val[0] += 0.0;
            break;
        case 3:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        case 4:
            prod_val[0] += -cts->call->pay_rec * coupon * DFR_pay;
            break;
        }
    }

    /*	End of payoff valuation
    ----------------------- */

    return err;
}

/*	Main pricing function */

/*	Launch the PDE */
Err cts_launch_algo(
    CTS         cts,
    CTS_UND     und,
    CTS_ADI_ARG adi_arg,

    /*	Result */
    double* multi_pv,
    double* onetime_pv,
    double* fwd_iv)
{
    double *temp_val = NULL, **temp_val_mc = NULL;

    int    i, ncol_init, ncol_tot;
    double dfstar;
    Err    err = NULL;

    if (adi_arg->pde_or_mc == 0)
    {
        smessage(
            "Launching adi, time steps requested: %d, actual: %d",
            adi_arg->initnstp,
            adi_arg->nstp);
    }
    else
    {
        smessage("Launching MC, Nb Steps: %d, Nb Paths: %d", adi_arg->nstp, adi_arg->npaths);
    }

    if (adi_arg->pde_or_mc == 0)
    {
        if (cts->call[0].is_midat)
        {
            ncol_init = 1;
        }
        else
        {
            ncol_init = 2;

            for (i = cts->call[0].exo_idx; i < cts->exo_leg->num_cpn; i++)
            {
                if (cts->exo_leg->cpn[i].type == 1 || cts->exo_leg->cpn[i].type == 2)
                {
                    ncol_init = 3;
                    break;
                }
            }
        }

        ncol_tot = ncol_init + adi_arg->has_one_time;

        temp_val = (double*)calloc(ncol_tot, sizeof(double));

        if (!temp_val)
        {
            err = "Memory allocation faillure in cts_launch_algo";
            goto FREE_RETURN;
        }

        /* launch the corresponding ADI */
        err = lgmSV_adi_UtPieceWise(
            adi_arg->nstp,
            adi_arg->time,
            adi_arg->date,
            adi_arg->nstppsi,
            adi_arg->nstpx,
            adi_arg->nstpz,
            und->model.dLambdaX,
            0.0,
            und->model.dPWTime,
            und->model.dSigma,
            und->model.iNbPWTime,
            und->model.dLambdaEps,
            und->model.dLvlEps,
            und->model.dAlpha,
            und->model.dRho,
            adi_arg->params,
            adi_arg->void_prm,
            adi_arg->is_event,
            adi_arg->ifr,
            und->mkt->yc,
            cts_payoff_4_lgmsv_adi,
            ncol_init,
            ncol_tot,
            &(temp_val[0]));

        if (err)
            return err;

        *multi_pv = temp_val[0];

        if (adi_arg->has_one_time)
        {
            for (i = 0; i < ncol_tot - ncol_init; i++)
            {
                onetime_pv[i] = temp_val[ncol_tot - 1 - i];
            }
            for (i = i; i < cts->num_calls; i++)
            {
                onetime_pv[i] = 0.0;
            }
        }
        else
        {
            if (onetime_pv)
            {
                onetime_pv[0] = 0.0;
            }
        }

        if (fwd_iv)
        {
            *fwd_iv = temp_val[1];
        }
    }
    else
    {
        /* Allocate MCEBParams */
        err = mceb_allocate_params(adi_arg->mcebparams, adi_arg->nb_event);
        if (err)
            goto FREE_RETURN;

        temp_val_mc = dmatrix(
            0,
            adi_arg->nb_event + 4,
            0,
            2 + adi_arg->mcebparams->iNbIndex + adi_arg->mcebparams->iCalcOneTime);

        if (!temp_val_mc)
        {
            err = "Memory allocation faillure in cts_launch_algo";
            goto FREE_RETURN;
        }

        dfstar = swp_f_df(und->mkt->today, und->model.lTStarDate, und->mkt->yc);

        if (und->model.iOne2F == 1)
        {
            err = lgmSV_mc_balsam_rev(
                adi_arg->nstp,
                adi_arg->nb_event,
                adi_arg->time,
                adi_arg->date,
                adi_arg->npaths,
                und->model.dLambdaX,
                adi_arg->dSigma,
                adi_arg->dAlpha,
                adi_arg->dLambdaEps,
                adi_arg->dLvlEps,
                adi_arg->dRho,
                NULL,
                NULL,
                NULL,
                adi_arg->params,
                adi_arg->void_prm,
                adi_arg->is_event,
                1,
                adi_arg->optimise,
                adi_arg->mcebparams,
                NULL,
                cts_payoff_4_lgmsv_mc,
                NULL,
                3,
                temp_val_mc);

            if (err)
                goto FREE_RETURN;

            *multi_pv = (temp_val_mc[3][0] - temp_val_mc[0][0]) * dfstar;
        }
        else
        {
            err = lgmSV2F_mc_balsam_rev(
                adi_arg->nstp,
                adi_arg->nb_event,
                adi_arg->time,
                adi_arg->date,
                adi_arg->npaths,
                und->model.dLambdaX,
                und->model.dLambdaX2,
                adi_arg->dSigma,
                adi_arg->dLGMAlpha,
                adi_arg->dLGMRho,
                adi_arg->dAlpha,
                adi_arg->dLambdaEps,
                adi_arg->dLvlEps,
                adi_arg->dRho,
                adi_arg->dRho2,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                adi_arg->params,
                adi_arg->void_prm,
                adi_arg->is_event,
                1,
                adi_arg->optimise,
                adi_arg->mcebparams,
                NULL,
                cts_payoff_4_lgmsv2F_mc,
                NULL,
                4,
                temp_val_mc);

            if (err)
                goto FREE_RETURN;

            *multi_pv = (temp_val_mc[4][0] - temp_val_mc[0][0]) * dfstar;
        }

        if (fwd_iv)
        {
            *fwd_iv = temp_val_mc[0][0] * dfstar;
        }
    }

FREE_RETURN:

    if (temp_val)
        free(temp_val);
    if (temp_val_mc)
        free_dmatrix(
            temp_val_mc,
            0,
            adi_arg->nb_event + 4,
            0,
            2 + adi_arg->mcebparams->iNbIndex + adi_arg->mcebparams->iCalcOneTime);

    return err;
}
