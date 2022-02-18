/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_GEN_DEAL                                        */
/*                                                                            */
/*      PURPOSE:                                                              */
/*                                                                            */
/*      AUTHORS:        Guillaume Amblard              			      */
/*                                                                            */
/*      DATE:           19th September 1993                                   */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:    XX                                                    */
/*                                                                            */
/*      FUNCTIONS USED: XXX_X_XXXXXXXXX                                       */
/*              Must include all imported function call made by the module    */
/*                                                                            */
/*      PARAMETERS:     <not applicable>                                      */
/*                                                                            */
/*      RETURNS:                                                              */
/*                                                                            */
/*      DATA ACCESSED:                                                        */
/*              Tables, files, and global variables accessed.                 */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           18th March, 1994                                      */
/*                                                                            */
/*      VERSION:        <not applicable>                                      */
/*                                                                            */
/*      REASON:         Restructuring for re-use                              */
/*                                                                            */
/*      REQUEST NO:     <not applicable>                                      */
/*                                                                            */
/*      DESCRIPTION:    <not applicable>                                      */
/*                                                                            */
/******************************************************************************/

/**
        REQUIRED:	Other utility routines developed by E.Auld & A. Litke
**/

#include "math.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_vol_interpol.h"

static void generate_leg_dates(Arg_Obj* arg, Leg_Obj* leg);
static int  set_date(Leg_Obj* target, int index, double value);
static int  set_time_in_years(Leg_Obj* target, int index, double value);

/* For a standard swap, crv->current has to point a yc_crv
   For a CMT or a CMS swap, crv->current has to point to a cmt_crv
*/
double value_swap(SrtCurvePtr crv, Swap_Obj* swap)
{
    Leg_Obj* leg;
    double   act_payment, disc_factor, swap_value = 0, leg_value = 0, value_date;
    int      i, j, num_leg, index_start, index_end;
    SwapType mess;
    String   yc_name, crv_name;

    num_leg    = Swap_get(swap, SWAP_NUM_LEG);
    value_date = Swap_get(swap, SWAP_VALUE_DATE);
    mess       = Swap_Iget(swap, SWAP_TYPE);

    if ((mess == CMS_MARGIN_FLOATING_SWAP) || (mess == ONE_LEG_CMS_MARGIN) ||
        (mess == CMT_MARGIN_FLOATING_SWAP) || (mess == ONE_LEG_CMT_MARGIN) ||
        (mess == CMS_FOR_TEC_MARGIN_FLOATING_SWAP) || (mess == ONE_LEG_CMS_FOR_TEC_MARGIN) ||
        (mess == TEC_MARGIN_FLOATING_SWAP) || (mess == ONE_LEG_TEC_MARGIN))
    {
        yc_name  = get_ycname_from_cmtcrv(crv);
        crv_name = get_fwdTname_from_cmtcrv(crv);
    }
    else
    {
        yc_name = crv_name = get_fwdTname_from_cmtcrv(crv);
    }
    for (i = 0; i < num_leg; i++)
    {
        leg  = Swap_field_get(swap, SWAP_LEG, i);
        mess = Leg_Iget(leg, LEG_TYPE);

        index_start = Leg_Iget(leg, LEG_INDEX_START);
        index_end   = Leg_Iget(leg, LEG_INDEX_END);

        if ((mess == CMS_MARGIN_LEG) || (mess == CMT_MARGIN_LEG) || (mess == TEC_MARGIN_LEG) ||
            (mess = CMS_FOR_TEC_MARGIN_LEG))
        {
            /* Evaluate the payments in the leg */
            leg_eval_payment(leg, crv_name);
        }
        else
        {
            /* Evaluate the payments in the leg */
            leg_eval_payment(leg, yc_name);
        }

        /* Compute the sum of the pv of the leg payments */
        crv = lookup_curve(yc_name);
        leg_eval_disc_factor(leg, yc_name);
        crv = lookup_curve(crv_name);

        for (j = index_start; j <= index_end; j++)
        {
            if (Leg_field_Dateget(leg, LEG_DATE, j) >= DTOL(value_date))
            {
                act_payment = Leg_field_get(leg, LEG_PAYMENT, j);
                disc_factor = Leg_field_get(leg, LEG_DISC_FACTOR, j);
                leg_value += act_payment * disc_factor;
            }
        }
        Leg_set(leg, LEG_VALUE, leg_value);
        swap_value += leg_value;
        leg_value = 0;
    }
    Swap_set(swap, SWAP_VALUE, swap_value);
    return (swap_value);
}

Swap_Obj* generate_swap(Arg_Obj* arg, SwapType mess)
{
    int       rec_pay;
    int       i;
    int       num_leg;
    double    value_date;
    Swap_Obj* swap;
    Leg_Obj*  leg[2]; /** Array of pointers **/

    swap       = new_Swap();
    rec_pay    = Arg_get(arg, ARG_REC_PAY);
    value_date = Arg_Dget(arg, ARG_VALUE_DATE);

    rec_pay = 1 - 2 * rec_pay;

    switch (mess)
    {
    case ONE_LEG_FIXED_AND_NOTIONALS:
        leg[0]  = generate_leg(arg, FIXED_AND_NOTIONALS_LEG, rec_pay);
        num_leg = 1;
        break;
    case FIXED_NOTIONALS_SWAP:
        leg[0]  = generate_leg(arg, FIXED_LEG, rec_pay);
        leg[1]  = generate_leg(arg, NOTIONALS_LEG, rec_pay);
        num_leg = 2;
        break;
    case FIXED_FLOATING_SWAP:
        leg[0]  = generate_leg(arg, FIXED_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FLOATING_LEG, -rec_pay);
        num_leg = 2;
        break;
    case ONE_LEG_FIXED_SWAP:
        leg[0]  = generate_leg(arg, FIXED_LEG, rec_pay);
        num_leg = 1;
        break;
    case DRS_MARGIN_FLOATING_SWAP:
        leg[0]  = generate_leg(arg, DRS_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FIXED_AND_NOTIONALS_LEG, -rec_pay);
        num_leg = 2;
        break;
    case DRS_MARGIN_FIXED_SWAP:
        leg[0]  = generate_leg(arg, DRS_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FIXED_LEG, -rec_pay);
        num_leg = 2;
        break;
    case ONE_LEG_DRS:
        leg[0]  = generate_leg(arg, DRS_LEG, rec_pay);
        num_leg = 1;
        break;
    case CMS_MARGIN_FLOATING_SWAP:
        leg[0]  = generate_leg(arg, CMS_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, NOTIONALS_LEG, rec_pay);
        num_leg = 2;
        break;
    case CMS_MARGIN_FIXED_SWAP:
        leg[0]  = generate_leg(arg, CMS_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FIXED_LEG, -rec_pay);
        num_leg = 2;
        break;
    case ONE_LEG_CMS_MARGIN:
        leg[0]  = generate_leg(arg, CMS_MARGIN_LEG, rec_pay);
        num_leg = 1;
        break;
    case CMT_MARGIN_FLOATING_SWAP:
        leg[0]  = generate_leg(arg, CMT_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, NOTIONALS_LEG, rec_pay);
        num_leg = 2;
        break;
    case CMT_MARGIN_FIXED_SWAP:
        leg[0]  = generate_leg(arg, CMT_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FIXED_LEG, rec_pay);
        num_leg = 2;
        break;
    case ONE_LEG_CMT_MARGIN:
        leg[0]  = generate_leg(arg, CMT_MARGIN_LEG, rec_pay);
        num_leg = 1;
        break;
    case CMS_FOR_TEC_MARGIN_FLOATING_SWAP:
        leg[0]  = generate_leg(arg, CMS_FOR_TEC_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, NOTIONALS_LEG, rec_pay);
        num_leg = 2;
        break;
    case CMS_FOR_TEC_MARGIN_FIXED_SWAP:
        leg[0]  = generate_leg(arg, CMS_FOR_TEC_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FIXED_LEG, rec_pay);
        num_leg = 2;
        break;
    case ONE_LEG_CMS_FOR_TEC_MARGIN:
        leg[0]  = generate_leg(arg, TEC_MARGIN_LEG, rec_pay);
        num_leg = 1;
        break;
    case TEC_MARGIN_FLOATING_SWAP:
        leg[0]  = generate_leg(arg, TEC_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, NOTIONALS_LEG, rec_pay);
        num_leg = 2;
        break;
    case TEC_MARGIN_FIXED_SWAP:
        leg[0]  = generate_leg(arg, TEC_MARGIN_LEG, rec_pay);
        leg[1]  = generate_leg(arg, FIXED_LEG, rec_pay);
        num_leg = 2;
        break;
    case ONE_LEG_TEC_MARGIN:
        leg[0]  = generate_leg(arg, TEC_MARGIN_LEG, rec_pay);
        num_leg = 1;
        break;
    case ONE_LEG_STANDARD_BOND:
        leg[0]  = generate_leg(arg, FIXED_AND_NOTIONALS_LEG, rec_pay);
        num_leg = 1;
        break;
    default:
        leg[0]  = generate_leg(arg, FIXED_AND_NOTIONALS_LEG, rec_pay);
        num_leg = 1;
        break;
    }
    Swap_set(swap, SWAP_TYPE, mess);
    Swap_set(swap, SWAP_NUM_LEG, num_leg);
    Swap_set(swap, SWAP_REC_PAY, rec_pay);
    Swap_set(swap, SWAP_VALUE_DATE, value_date);

    for (i = 0; i < num_leg; i++)
    {
        Swap_field_set(swap, SWAP_LEG, i, leg[i]);
    }
    return (swap);
}

/* ======================================================================== */

/* Used only for CMS-CMT swaps : stores the forward swap rate at each fixing
   date
   cmt_crv refers to a treasury forward curve market*/
static void leg_eval_fwd_swaps(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    double      fwd;
    Date        today;
    Ddate       period_start;
    int         i, index_start, index_end;
    SwapDP      swapdp;
    String      yc_name, crv_name;
    SrtCurvePtr yc_crv;
    CMT_Param*  cmt;
    String      ccy_str;

    cmt      = get_cmtparam_from_cmtcrv(cmt_crv);
    crv_name = get_fwdTname_from_cmtcrv(cmt_crv);
    yc_name  = get_ycname_from_cmtcrv(cmt_crv);
    yc_crv   = lookup_curve(yc_name);
    ccy_str  = get_curve_ccy(yc_crv);

    today = DTOL(Leg_get(leg, LEG_TODAY));

    index_start = Leg_get(leg, LEG_INDEX_START);
    index_end   = Leg_get(leg, LEG_INDEX_END);

    swapdp.nfp        = (int)cmt->cmt_mat * cmt->swap_compd;
    swapdp.end        = swapdp.nfp;
    swapdp.direction  = FWD;
    swapdp.compd      = cmt->swap_compd;
    swapdp.basis_code = cmt->swap_basis_code;
    swapdp.spot_lag   = cmt->spot_lag;

    for (i = index_start; i <= index_end; i++)
    {
        period_start             = Leg_field_get(leg, LEG_DATE, i);
        swapdp.start             = (Date)period_start;
        swapdp.first_full_fixing = swapdp.start;
        swp_f_ForwardRate_SwapDP(&swapdp, yc_name, ccy_str, &fwd);
        Leg_field_set(leg, LEG_SWAP_FORWARD, i, fwd);
        Leg_field_set(leg, LEG_FORWARD, i, fwd);
    }
    /* Makes sure that crv_list->current points to the CMT crv */
    cmt_crv = lookup_curve(crv_name);
}

/* Used only for CMS-CMT swaps : stores the forward treasury  rate at each
   fixing date
   cmt_crv->current has to point to a treasury forward curve market*/
static void leg_eval_fwd_treas(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    double     fwd_treas;
    Date       today;
    String     crv_name;
    Ddate      period_start;
    int        i, index_start, index_end;
    CMT_Param* cmt;
    Err        err = NULL;

    cmt      = get_cmtparam_from_cmtcrv(cmt_crv);
    crv_name = get_fwdTname_from_cmtcrv(cmt_crv);

    today = DTOL(Leg_get(leg, LEG_TODAY));

    index_start = Leg_get(leg, LEG_INDEX_START);
    index_end   = Leg_get(leg, LEG_INDEX_END);

    for (i = index_start; i < index_end; i++)
    {
        period_start = Leg_field_get(leg, LEG_DATE, i);
        err          = fwd_treas_rate(period_start, cmt_crv, &fwd_treas);
        Leg_field_set(leg, LEG_FORWARD, i, fwd_treas);
    }

    cmt_crv = lookup_curve(crv_name);
}

/* Used only for CMS-CMT swaps : stores the convexity adjusted rate at each
        fixing date
   cmt_crv->current has to point to a treasury forward curve market*/
static void leg_eval_cm_rate(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    Err            err;
    double         fwd, vol, cm_rate;
    Date           today;
    Ddate          period_start, pay_date;
    long           fixing_date;
    double         delay, maturity;
    int            i, index_start, index_end, nfp;
    BasisCode      leg_basis_code;
    SrtCompounding cm_compounding;
    BasisCode      cm_basiscode;
    LegType        mess;
    CMT_Param*     cmt_param;

    cmt_param = get_cmtparam_from_cmtcrv(cmt_crv);

    today = DTOL(Leg_get(leg, LEG_TODAY));

    index_start = Leg_get(leg, LEG_INDEX_START);
    index_end   = Leg_get(leg, LEG_INDEX_END);

    leg_basis_code = Leg_get(leg, LEG_BASIS_CODE);
    mess           = Leg_get(leg, LEG_TYPE);

    /* If CMS leg, use swap parameters */
    if (mess == CMS_MARGIN_LEG)
    {
        cm_compounding = cmt_param->swap_compd;
        cm_basiscode   = cmt_param->swap_basis_code;
    }
    else
        /* If CMT leg, use bond parameters */
        if ((mess == CMT_MARGIN_LEG) || (mess == TEC_MARGIN_LEG))
    {
        cm_compounding = cmt_param->bond_compd;
        cm_basiscode   = cmt_param->bond_basis_code;
    }

    nfp = (int)(cmt_param->cmt_mat) * cm_compounding;

    for (i = index_start; i < index_end; i++)
    {
        period_start = Leg_field_get(leg, LEG_DATE, i);
        /* OVE fix: a swap is described by theoretical payment dates: spot lag after fixing,
           vol is effective only until fixing date */
        fixing_date = add_unit((long)period_start, -cmt_param->spot_lag, SRT_BDAY, SUCCEEDING);
        pay_date    = Leg_field_get(leg, LEG_DATE, i + 1);
        delay       = (double)coverage((Date)period_start, (Date)pay_date, cm_basiscode);
        /* OVE fix: maturity is the maturity of the swaption : vol starts from today */
        maturity = ((double)fixing_date - (double)today) * YEARS_IN_DAY;

        fwd = Leg_field_get(leg, LEG_FORWARD, i);
        vol = Leg_field_get(leg, LEG_VOLATILITY, i);

        /* This is needed for the cms function considers a 365 basis */
        if (cm_basiscode == BASIS_ACT_360)
            fwd = fwd * 365.0 / 360.0;

        if (err = swp_f_cmsrate(
                fwd,
                nfp,
                cm_compounding,
                vol,
                maturity,
                delay,
                DEFAULT_CMS_DELTA,
                MAX_CMS_SWAPS,
                SRT_LOGNORMAL,
                &cm_rate))
            cm_rate = -1;

        if (cm_basiscode == BASIS_ACT_360)
            cm_rate = cm_rate * 360.0 / 365.0;

        Leg_field_set(leg, LEG_CM_RATE, i + 1, cm_rate);
    }
}

/* Used only for CMS-CMT swaps
   cmt_crv->current has to point to a treasury forward curve market*/
static Err leg_eval_vols(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    Err    err = NULL;
    double vol;
    Date   today;

    Ddate  period_start;
    Ddate  period_end;
    double strike;

    CMT_Param* cmt_param;
    int        i, index_start, index_end;
    int        method = 0;

    today = DTOL(Leg_get(leg, LEG_TODAY));

    cmt_param = get_cmtparam_from_cmtcrv(cmt_crv);

    index_start = Leg_get(leg, LEG_INDEX_START);
    index_end   = Leg_get(leg, LEG_INDEX_END);

    for (i = index_start; i < index_end; i++)
    {
        period_start = Leg_field_get(leg, LEG_DATE, i);
        period_end   = add_unit((long)period_start, cmt_param->cmt_mat, SRT_YEAR, SUCCEEDING);
        strike       = Leg_field_get(leg, LEG_FORWARD, i);
        /* OVE fix: a swap is described by theoretical payment dates: spot lag after fixing,
           vol is effective only until fixing date
                        tgt_exp_date = add_unit((long)period_start, -cmt_param->spot_lag, SRT_BDAY,
           SUCCEEDING);
        */
        /* use the getvol function */
        if (cmt_param->flatvol != 0.)
        {
            vol = cmt_param->flatvol;
        }
        else
        {
            err = cmt_param->cmt_getvol((Date)period_start, (Date)period_end, strike, &vol);
            if (err)
                return err;
        }

        Leg_field_set(leg, LEG_SWAP_VOLATILITY, i, vol);
        /* By default, VOLATILITY is SWAP_VOLATILITY */
        Leg_field_set(leg, LEG_VOLATILITY, i, vol);
    }

    return err;
}

/* Used only for CMS-CMT swaps
   cmt_crv->current has to point to a treasury forward curve market*/
static void leg_eval_fwd_adj_vols(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    double vol;
    Date   today;

    CMT_Param* cmt_param;

    int         i, index_start, index_end;
    SRT_Boolean use_prop_vol_flg;
    double      fwd_R_rate, fwd_T_rate;

    today = DTOL(Leg_get(leg, LEG_TODAY));

    cmt_param        = get_cmtparam_from_cmtcrv(cmt_crv);
    use_prop_vol_flg = cmt_param->use_prop_vol_flg;

    index_start = Leg_get(leg, LEG_INDEX_START);
    index_end   = Leg_get(leg, LEG_INDEX_END);

    for (i = index_start; i < index_end; i++)
    {
        vol = Leg_field_get(leg, LEG_SWAP_VOLATILITY, i);
        if (use_prop_vol_flg == SRT_YES)
        {
            fwd_R_rate = Leg_field_get(leg, LEG_SWAP_FORWARD, i);
            fwd_T_rate = Leg_field_get(leg, LEG_FORWARD, i);
            vol        = vol * fwd_R_rate / fwd_T_rate;
        }

        Leg_field_set(leg, LEG_VOLATILITY, i, vol);
    }
}

/* ======================================================================== */

Leg_Obj* generate_leg(Arg_Obj* arg, LegType mess, int rec_pay)
{
    Leg_Obj* leg;
    Date     start_init, end_init, today_date;
    Date     start, end;

    int            i, index_start, index_end, toalloc;
    SrtBasisCode   basis_code;
    SrtCompounding compounding;
    double         notional, initial_not, final_not, cvg, strike;
    double         spread;
    Ddate          first_full_fixing;
    Date           value_date;

    if (Arg_Iget(arg, ARG_DATE_DIR) == FWD)
    {
        toalloc = Arg_Iget(arg, ARG_NUM_FULL_PERIOD) + 10;
    }
    else
    {
        toalloc = (year(Arg_Dateget(arg, ARG_END)) - year(Arg_Dateget(arg, ARG_START)) + 1) *
                      Arg_Iget(arg, ARG_COMPD) +
                  5;
    }

    leg_allocate(&leg, toalloc);

    /* Generates the leg dates, and sets the arg index start and end */
    generate_leg_dates(arg, leg);

    /* Gets arguments from the Arg */
    today_date = DTOL(Arg_get(arg, ARG_TODAY));
    notional   = Arg_get(arg, ARG_NOTIONAL);
    value_date = (Date)Arg_get(arg, ARG_VALUE_DATE);

    start_init        = DTOL(Arg_get(arg, ARG_START));
    end_init          = DTOL(Arg_get(arg, ARG_END));
    index_start       = Arg_Iget(arg, ARG_INDEX_START);
    index_end         = Arg_Iget(arg, ARG_INDEX_END);
    basis_code        = Arg_Iget(arg, ARG_BASIS_CODE);
    compounding       = Arg_Iget(arg, ARG_COMPD);
    first_full_fixing = Arg_Dget(arg, ARG_FIRST_FULL_FIXING);
    initial_not       = Arg_Dget(arg, ARG_INITIAL_NOT);
    final_not         = Arg_Dget(arg, ARG_FINAL_NOT);

    /* Sets arguments in the Leg */
    Leg_Iset(leg, LEG_TYPE, mess);
    Leg_set(leg, LEG_TODAY, (double)today_date);
    Leg_set(leg, LEG_VALUE_DATE, (double)value_date);
    Leg_set(leg, LEG_NOTIONAL, notional);
    Leg_Iset(leg, LEG_REC_PAY, rec_pay);
    Leg_set(leg, LEG_START, (double)start_init);
    Leg_set(leg, LEG_END, (double)end_init);
    Leg_Iset(leg, LEG_INDEX_START, index_start);
    Leg_Iset(leg, LEG_INDEX_END, index_end);
    Leg_Iset(leg, LEG_BASIS_CODE, basis_code);
    Leg_Iset(leg, LEG_COMPOUNDING, compounding);
    Leg_Dset(leg, LEG_FIRST_FULL_FIXING, first_full_fixing);

    switch (mess)
    {
    case NOTIONALS_LEG:
        Leg_set(leg, LEG_STRIKE, 0.0);
        Leg_Dset(leg, LEG_INITIAL_NOT, initial_not);
        Leg_Dset(leg, LEG_FINAL_NOT, final_not);
        break;
    case FIXED_LEG:
        strike = Arg_get(arg, ARG_STRIKE);
        Leg_set(leg, LEG_STRIKE, strike);
        /*	start = start_init;   changed ea for busday bug 26/11/93 5:10  */
        start = (Date)Leg_field_get(leg, LEG_DATE, index_start);
        for (i = index_start + 1; i <= index_end; i++)
        {
            end = (Date)Leg_field_get(leg, LEG_DATE, i);
            cvg = coverage(start, end, basis_code);
            Leg_field_set(leg, LEG_COVERAGE, i, cvg);
            start = end;
        }
        break;
    case FLOATING_LEG:
    case DRS_LEG:
        Leg_set(leg, LEG_STRIKE, 0);
        /*	start = start_init;   changed ea for busday bug 26/11/93 5:10  */
        start = (Date)Leg_field_get(leg, LEG_DATE, index_start);
        for (i = index_start + 1; i <= index_end; i++)
        {
            end = (Date)Leg_field_get(leg, LEG_DATE, i);
            cvg = coverage(start, end, basis_code);
            Leg_field_set(leg, LEG_COVERAGE, i, cvg);
            start = end;
        }
        break;
    case CMS_MARGIN_LEG:
        Leg_set(leg, LEG_STRIKE, 0);
        spread = Arg_get(arg, ARG_SPREAD);
        Leg_set(leg, LEG_SPREAD, spread);
        start = (Date)Leg_field_get(leg, LEG_DATE, index_start);
        for (i = index_start + 1; i <= index_end; i++)
        {
            end = (Date)Leg_field_get(leg, LEG_DATE, i);
            cvg = coverage(start, end, basis_code);
            Leg_field_set(leg, LEG_COVERAGE, i, cvg);
            start = end;
        }
        break;
    case CMT_MARGIN_LEG:
        Leg_set(leg, LEG_STRIKE, 0);
        spread = Arg_get(arg, ARG_SPREAD);
        Leg_set(leg, LEG_SPREAD, spread);
        start = (Date)Leg_field_get(leg, LEG_DATE, index_start);
        for (i = index_start + 1; i <= index_end; i++)
        {
            end = (Date)Leg_field_get(leg, LEG_DATE, i);
            cvg = coverage(start, end, basis_code);
            Leg_field_set(leg, LEG_COVERAGE, i, cvg);
            start = end;
        }
        break;
    case TEC_MARGIN_LEG:
        Leg_set(leg, LEG_STRIKE, 0);
        spread = Arg_get(arg, ARG_SPREAD);
        Leg_set(leg, LEG_SPREAD, spread);
        start = (Date)Leg_field_get(leg, LEG_DATE, index_start);
        for (i = index_start + 1; i <= index_end; i++)
        {
            end = (Date)Leg_field_get(leg, LEG_DATE, i);
            cvg = 0.0; /* coverage(start, end, basis_code); */
            Leg_field_set(leg, LEG_COVERAGE, i, cvg);
            start = end;
        }
        break;

    case FIXED_AND_NOTIONALS_LEG:
    default:
        strike = Arg_get(arg, ARG_STRIKE);
        Leg_set(leg, LEG_STRIKE, strike);
        /*	start = start_init;   changed ea for busday bug 26/11/93 5:10  */
        start = (Date)Leg_field_get(leg, LEG_DATE, index_start);
        for (i = index_start + 1; i <= index_end; i++)
        {
            end = (Date)Leg_field_get(leg, LEG_DATE, i);
            cvg = coverage(start, end, basis_code);
            Leg_field_set(leg, LEG_COVERAGE, i, cvg);
            start = end;
        }
        break;
    }

    return leg;
}

/* For a standard leg, crv->current has to point to a yc_crv
   For a CMT or a CMS swap, crv->current has to point to a cmt_crv
*/
void leg_eval_payment(Leg_Obj* leg, String curvename)
{
    Message        mess;
    Date           start, end;
    SrtBasisCode   basis_code;
    SrtCompounding compd;
    int            i, index_end, index_start, rec_pay;
    double         coupon, notional, initial_not, final_not, spread, strike, cvg, payment;
    SrtCurvePtr    crv = lookup_curve(curvename);

    notional    = Leg_Dget(leg, LEG_NOTIONAL);
    rec_pay     = Leg_Iget(leg, LEG_REC_PAY);
    index_start = Leg_Iget(leg, LEG_INDEX_START);
    index_end   = Leg_Iget(leg, LEG_INDEX_END);
    start       = DTOL(Leg_field_get(leg, LEG_DATE, index_start));
    end         = DTOL(Leg_get(leg, LEG_END));
    basis_code  = Leg_Iget(leg, LEG_BASIS_CODE);
    compd       = Leg_Iget(leg, LEG_COMPOUNDING);
    mess        = Leg_Iget(leg, LEG_TYPE);

    switch (mess)
    {
    case NOTIONALS_LEG:
        initial_not = Leg_get(leg, LEG_INITIAL_NOT) * notional;
        final_not   = Leg_get(leg, LEG_FINAL_NOT) * notional;
        Leg_field_set(leg, LEG_PAYMENT, index_start, -rec_pay * initial_not);
        Leg_field_set(leg, LEG_PAYMENT, index_end, rec_pay * final_not);
        break;
    case FIXED_LEG:
        strike = Leg_get(leg, LEG_STRIKE);
        for (i = index_start + 1; i <= index_end; i++)
        {
            cvg     = Leg_field_get(leg, LEG_COVERAGE, i);
            payment = notional * rec_pay * strike * cvg;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
        }
        break;
    case FLOATING_LEG:
        for (i = index_start + 1; i <= index_end; i++)
        {
            end     = DTOL(Leg_field_get(leg, LEG_DATE, i));
            coupon  = swp_f_fwdcash((double)start, (double)end, (BasisCode)basis_code, curvename);
            cvg     = Leg_field_get(leg, LEG_COVERAGE, i);
            payment = coupon * cvg * notional * rec_pay;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
            start = end;
        }
        break;
    case DRS_LEG:
        break;
    case CMS_MARGIN_LEG:
        spread = Leg_Dget(leg, LEG_SPREAD);
        Leg_field_set(leg, LEG_PAYMENT, index_start, 0.0);
        for (i = index_start + 1; i <= index_end; i++)
        {
            cvg     = Leg_field_get(leg, LEG_COVERAGE, i);
            coupon  = Leg_field_get(leg, LEG_CM_RATE, i);
            payment = (coupon + spread) * cvg * notional;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
        }
        break;

    case CMS_FOR_TEC_MARGIN_LEG:
        spread = Leg_Dget(leg, LEG_SPREAD);
        Leg_field_set(leg, LEG_PAYMENT, index_start, 0.0);
        for (i = index_start + 1; i <= index_end; i++)
        {
            coupon  = Leg_field_get(leg, LEG_CM_RATE, i);
            payment = (exp(1.0 / (double)compd * log(1 + coupon + spread)) - 1.0) * notional;
            ;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
        }
        break;

    case CMT_MARGIN_LEG:
        spread = Leg_Dget(leg, LEG_SPREAD);
        /* This is computed here for the vols...
           ... depend on the iteration on the fwd_T curve   */
        /* everything is generated from index_start to index_end */
        /* however cash flows are paid at the end of the period so we add 1 */
        /* when considering the cms rate */
        leg_eval_fwd_treas(leg, crv);
        leg_eval_vols(leg, crv);
        leg_eval_fwd_adj_vols(leg, crv);
        leg_eval_cm_rate(leg, crv);

        Leg_field_set(leg, LEG_PAYMENT, index_start, 0.0);
        for (i = index_start + 1; i <= index_end; i++)
        {
            cvg     = Leg_field_get(leg, LEG_COVERAGE, i);
            coupon  = Leg_field_get(leg, LEG_CM_RATE, i);
            payment = (coupon + spread) * cvg * notional;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
        }
        break;

    case TEC_MARGIN_LEG:
        spread = Leg_Dget(leg, LEG_SPREAD);
        /* This is computed here for the vols...
           ... depend on the iteration on the fwd_T curve   */
        leg_eval_fwd_treas(leg, crv);
        leg_eval_vols(leg, crv);
        leg_eval_fwd_adj_vols(leg, crv);
        leg_eval_cm_rate(leg, crv);

        Leg_field_set(leg, LEG_PAYMENT, index_start, 0.0);
        for (i = index_start + 1; i <= index_end; i++)
        {
            coupon  = Leg_field_get(leg, LEG_CM_RATE, i);
            payment = (exp(1.0 / (double)compd * log(1 + coupon + spread)) - 1.0) * notional;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
        }
        break;
    case FIXED_AND_NOTIONALS_LEG:
        strike = Leg_get(leg, LEG_STRIKE);
        Leg_field_set(leg, LEG_PAYMENT, index_start, -rec_pay * notional);
        for (i = index_start + 1; i <= index_end; i++)
        {
            cvg     = Leg_field_get(leg, LEG_COVERAGE, i);
            payment = notional * rec_pay * strike * cvg;
            Leg_field_set(leg, LEG_PAYMENT, i, payment);
        }
        payment += rec_pay * notional;
        Leg_field_set(leg, LEG_PAYMENT, index_end, payment);
        break;
    default: /** In the case of FIXED CASHFLOWS, no need to
                        redo any calculations **/
        break;
    }
}

void leg_eval_disc_factor(Leg_Obj* leg, String ycname)
{
    int    i, index_start, index_end;
    Date   end;
    double df;

    index_start = Leg_get(leg, LEG_INDEX_START);
    index_end   = Leg_get(leg, LEG_INDEX_END);

    for (i = index_start; i <= index_end; i++)
    {
        end = DTOL(Leg_field_get(leg, LEG_DATE, i));
        df  = swp_f_df(0.0, end, ycname);
        Leg_field_set(leg, LEG_DISC_FACTOR, i, df);
    }
}

static void generate_leg_dates(Arg_Obj* arg, Leg_Obj* leg)
{
    SwapDP     swap;
    BusDayConv conv;
    DateList   list;
    Date       today;
    int        i, k;

    swap.start             = Arg_Dateget(arg, ARG_START);
    swap.end               = Arg_Dateget(arg, ARG_END);
    swap.basis_code        = Arg_Iget(arg, ARG_BASIS_CODE);
    swap.first_full_fixing = (Date)Arg_Dget(arg, ARG_FIRST_FULL_FIXING);
    swap.nfp               = Arg_Iget(arg, ARG_NUM_FULL_PERIOD);
    swap.direction         = Arg_Iget(arg, ARG_DATE_DIR);
    swap.compd             = Arg_Iget(arg, ARG_COMPD);
    conv                   = MODIFIED_SUCCEEDING;
    today                  = Arg_Dateget(arg, ARG_TODAY);

    list = SwapDP_to_DateList(&swap, conv);

    k = 0;
    while (today >= list.date[k] && k < list.len)
        k++;

    Leg_field_Dateset(leg, LEG_DATE, 0, today);
    Leg_field_Dset(leg, LEG_TIME_IN_YEARS, 0, 0.0);

    for (i = k; i < list.len; i++)
    {
        Leg_field_Dateset(leg, LEG_DATE, i + 1 - k, list.date[i]);
        Leg_field_Dateset(leg, LEG_TIME_IN_YEARS, i + 1 - k, (list.date[i] - today) / 365.0);
    }

    if (k == 0)
        k++;
    Arg_Iset(arg, ARG_INDEX_START, 1 - (today >= list.date[k - 1]));
    Arg_Iset(arg, ARG_INDEX_END, list.len - k + 1 - (today >= list.date[k - 1]));

    srt_free(list.date);
}

static int set_date(Leg_Obj* target, int index, double value)
{
    Leg_field_Dset(target, LEG_DATE, index, value);

    return (0);
}

static int set_time_in_years(Leg_Obj* target, int index, double value)
{
    Leg_field_Dset(target, LEG_TIME_IN_YEARS, index, value);

    return (0);
}

/* ======================================================================== */

Err leg_cms_populate(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    Err err = NULL;

    leg_eval_fwd_swaps(leg, cmt_crv);
    err = leg_eval_vols(leg, cmt_crv);
    if (err)
        return err;
    leg_eval_cm_rate(leg, cmt_crv);

    return NULL;
}

Err leg_cmt_populate(Leg_Obj* leg, SrtCurvePtr cmt_crv)
{
    Err err = NULL;

    leg_eval_fwd_swaps(leg, cmt_crv);
    err = leg_eval_vols(leg, cmt_crv);
    return err;
}
