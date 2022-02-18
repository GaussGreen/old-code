/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_CMT_SWAP	                                      */
/*                                                                            */
/*      PURPOSE:        Compute a few values on CMT swaps... 				*/
/*                                                                            */
/*      AUTHORS:        Olivier VAN EYSEREN                    		      */
/*                                                                            */
/*      DATE:           24th June 1996                                     */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:     ........
 */
/*                                                                            */
/*      FUNCTIONS USED: XXX_X_XXXXXXXXX                                       */
/*              Must include all imported function call made by the module    */
/*                                                                            */
/*      PARAMETERS:     <not applicable>                                      */
/*                                                                            */
/*      RETURNS:                                                              */
/*                                                                            */
/*      DATA ACCESSED:                                                        */
/*      	<none>                                                        */
/*                                                                            */
/******************************************************************************/

#include "num_h_allhdr.h"
#include "swp_h_all.h"
#include "swp_h_cmt_swap.h"
#include "swp_h_swap_simple.h"

#define MAXITER 3

Err swp_f_CMT_swap(
    long        start,
    long        end_or_nfp,
    String      cmt_freq_str,
    String      cmt_basis_str,
    double      cmt_libor_spread, /* absolute (ex: -5bp or -0.0005) */
    double      strike,           /* absolute (ex: 5% or 0.05)*/
    double      initial_exchange,
    double      final_exchange,
    String      info_str,
    SrtCurvePtr crv,
    double*     value)
{
    Err         err;
    SwapDP      cmt_swapdp;
    SwapDP      swp_swapdp;
    Date        today;
    Date        spot;
    SwapOutput  swap_info;
    SrtCurvePtr yc_crv;
    SrtCurvePtr cmt_crv;
    String      yc_name;
    String      cmt_name;
    Arg_Obj*    cmt_arg;

    Swap_Obj* cmt_swap;
    Leg_Obj*  cmt_leg;

    CMT_Param* cmt_param;
    CMTCode    cmt_code;

    double swap_price;
    double fixed_leg_price;
    double level;
    double cmt_val;
    double cmt_price;
    double cmt_spread_price;
    double my_notional;

    int    niter;
    double yans;
    int    count;
    double nstop;
    double a[3], b[3];
    double spread_value;

    /* Checks the underlying type */
    if (!ISCURVETYPE(crv, CMT_CURVE))
        return serror("Wrong underlying type: need a CMT_CURVE");

    /* Sets conventions and dates for the cmt leg*/
    if (err = swp_f_initSwapDP(start, end_or_nfp, cmt_freq_str, cmt_basis_str, &cmt_swapdp))
    {
        return err;
    }

    /* Interprets the result wanted */
    if (err = interp_swap_message(info_str, &swap_info))
    {
        return err;
    }

    /* Gets yield curve from cmt_und */
    cmt_name = get_curve_name(crv);
    yc_name  = get_ycname_from_cmtcrv(crv);
    yc_crv   = lookup_curve(yc_name);
    if (!yc_crv)
        return serror("Could not find the %s yield curve", yc_name);

    /* Get today */
    today = get_clcndate_from_yldcrv(yc_crv);
    spot  = get_spotdate_from_yldcrv(yc_crv);

    /* Checks that start is at least equal to spot */
    if (start < spot)
        return serror("Start of swap should be after spot...");

    /* Sets the notional for the calculation */
    my_notional     = 555555555;
    fixed_leg_price = 0.00;

    /* Initialises the cmt_arg */
    crv = lookup_curve(cmt_name);
    init_Arg(&cmt_arg, (int)GENERATE_SWAP);
    Arg_Dateset(cmt_arg, ARG_TODAY, today);
    Arg_Dateset(cmt_arg, ARG_VALUE_DATE, today);
    Arg_Dateset(cmt_arg, ARG_START, start);
    Arg_Dateset(cmt_arg, ARG_FIRST_FULL_FIXING, start);
    Arg_set(cmt_arg, ARG_NUM_FULL_PERIOD, cmt_swapdp.nfp);
    Arg_set(cmt_arg, ARG_END, cmt_swapdp.end);
    Arg_set(cmt_arg, ARG_DATE_DIR, cmt_swapdp.direction);
    Arg_set(cmt_arg, ARG_COMPD, cmt_swapdp.compd);
    Arg_set(cmt_arg, ARG_BASIS_CODE, cmt_swapdp.basis_code);
    Arg_Dset(cmt_arg, ARG_STRIKE, 0.00);
    Arg_set(cmt_arg, ARG_NOTIONAL, my_notional);
    Arg_Dset(cmt_arg, ARG_INITIAL_NOT, 0.00);
    Arg_Dset(cmt_arg, ARG_FINAL_NOT, 0.00);
    Arg_set(cmt_arg, ARG_INFO_TYPE, swap_info);

    /* Gets the swap compounding and basis from the CMT_Param */
    cmt_param = get_cmtparam_from_cmtcrv(crv);
    cmt_code  = cmt_param->cmt_code;

    /* Initialises the swp_swpdp */
    err = swp_f_setSwapDP(
        start, end_or_nfp, cmt_param->swap_compd, cmt_param->swap_basis_code, &swp_swapdp);
    if (err)
        return err;

    /* Compute result wanted */
    switch (swap_info)
    {
    case COMPUTE_PV:
    case COMPUTE_FWD_PV:
        /* Check does not initilise the notionals and the fixed rate */
        if ((initial_exchange != 0.00 || final_exchange != 0.00) && (strike != 0))
            return serror("Cannot create a LIBOR leg AND a FIXED leg");
        /* Sets Libor leg */
        if (initial_exchange != 0.00 || final_exchange != 0.00)
        {
            Arg_Dset(cmt_arg, ARG_INITIAL_NOT, initial_exchange);
            Arg_Dset(cmt_arg, ARG_FINAL_NOT, final_exchange);
        }
        else
        /* Or compute separately the PV of the Fixed leg (with swap)*/
        {
            Arg_Dset(cmt_arg, ARG_STRIKE, strike);
            err =
                zcswap(&swp_swapdp, strike, 0.0, 0.0, COMPUTE_PV, yc_name, today, &fixed_leg_price);
            fixed_leg_price *= my_notional;
        }
        cmt_crv = lookup_curve(cmt_name);
        /* Sets the cmt_libor spread */
        Arg_Dset(cmt_arg, ARG_SPREAD, cmt_libor_spread);
        /* Build the CMT+spread / LIBOR swap (the vols depend on the fwdT rate */
        if (cmt_code == TEC10)
            cmt_swap = generate_swap(cmt_arg, TEC_MARGIN_FLOATING_SWAP);
        else
            cmt_swap = generate_swap(cmt_arg, CMT_MARGIN_FLOATING_SWAP);
        /* Extracts the CMT+spread Leg of the CMS - Libor swap */
        cmt_leg = Swap_field_get(cmt_swap, SWAP_LEG, 0);
        /* Computes and store fwd swap rates, and swap vols */
        leg_cmt_populate(cmt_leg, cmt_crv);
        /* Price of the swap */
        swap_price = value_swap(cmt_crv, cmt_swap);
        swap_price -= fixed_leg_price;
        *value = swap_price / my_notional;
        crv    = lookup_curve(cmt_name);
        break;

    case COMPUTE_FWD_RATE:
        cmt_crv = lookup_curve(cmt_name);
        /* Sets the cmt_libor spread to its value  */
        Arg_Dset(cmt_arg, ARG_SPREAD, cmt_libor_spread);
        /* Sets Libor leg to 0 */
        Arg_Dset(cmt_arg, ARG_INITIAL_NOT, 0);
        Arg_Dset(cmt_arg, ARG_FINAL_NOT, 0);
        /* Sets the fixed rate to 0 (though it is not used) */
        Arg_set(cmt_arg, ARG_STRIKE, 0.0);
        /* Build the CMT-spread leg (the vols depend on the fwdT rate */
        if (cmt_code == TEC10)
            cmt_swap = generate_swap(cmt_arg, ONE_LEG_TEC_MARGIN);
        else
            cmt_swap = generate_swap(cmt_arg, ONE_LEG_CMT_MARGIN);
        /* Extracts the CMT-spread Leg of the CMS - Libor swap */
        cmt_leg = Swap_field_get(cmt_swap, SWAP_LEG, 0);
        /* Computes and store fwd swap rates, and swap vols */
        leg_cmt_populate(cmt_leg, cmt_crv);
        /* Prices the swap */
        cmt_val = value_swap(cmt_crv, cmt_swap);
        /* Gets the equivalent swap level */
        err = zcswap(&swp_swapdp, 1.0, 0.0, 0.0, COMPUTE_PV, yc_name, today, &level);
        /* The forward rate is given by the ratio of the values */
        *value = cmt_val / level / my_notional;
        crv    = lookup_curve(cmt_name);
        break;

    case COMPUTE_FWD_SPREAD:
    case COMPUTE_SPREAD:
        cmt_crv = lookup_curve(cmt_name);
        /* Sets the cmt_libor spread to 0.0 */
        Arg_Dset(cmt_arg, ARG_SPREAD, 0.0);
        /* Sets the fixed rate to 0 (though it is not used) */
        Arg_set(cmt_arg, ARG_STRIKE, 0.0);
        /* Sets Libor leg */
        Arg_Dset(cmt_arg, ARG_INITIAL_NOT, 1);
        Arg_Dset(cmt_arg, ARG_FINAL_NOT, 1);
        /* Build the CMT-spread / LIBOR swap (the vols depend on the fwdT rate */
        if (cmt_code == TEC10)
            cmt_swap = generate_swap(cmt_arg, TEC_MARGIN_FLOATING_SWAP);
        else
            cmt_swap = generate_swap(cmt_arg, CMT_MARGIN_FLOATING_SWAP);
        /* Extracts the CMT-spread Leg of the CMS - Libor swap */
        cmt_leg = Swap_field_get(cmt_swap, SWAP_LEG, 0);
        /* Computes and store fwd swap rates, and swap vols */
        leg_cmt_populate(cmt_leg, cmt_crv);
        /* Prices the swap without spread*/
        cmt_price = value_swap(cmt_crv, cmt_swap);
        /* Sets the cmt_libor spread to -0.01 */
        cmt_leg = Swap_field_get(cmt_swap, SWAP_LEG, 0);
        Leg_set(cmt_leg, LEG_SPREAD, -0.01);
        /* Prices the swap with a spread of -0.01 (vol have not changed) */
        cmt_spread_price = value_swap(cmt_crv, cmt_swap);
        /* Guess the value of the spread for a linear dependence */
        spread_value = cmt_price * 0.01 / (cmt_spread_price - cmt_price);

        if (cmt_code != TEC10)
        /* FOr NON TEC compounding, price linear in spread */
        {
            *value = spread_value;
        }
        else
        /* FOr TEC compounding, solve using Newton from this point  */
        {
            /* Sets a few useful parameters for Newton*/
            niter = 5;
            yans  = 0.0;
            count = 0;
            nstop = 0.0;

            a[0] = spread_value;
            Leg_set(cmt_leg, LEG_SPREAD, a[0]);
            b[0] = value_swap(cmt_crv, cmt_swap);

            a[1] = a[0] + 0.0001;
            Leg_set(cmt_leg, LEG_SPREAD, a[1]);
            b[1] = value_swap(cmt_crv, cmt_swap);

            a[2] = a[1] + 0.0001;
            while ((count < MAXITER) && (nstop < 1.0))
            {
                Leg_set(cmt_leg, LEG_SPREAD, a[2]);
                b[2] = value_swap(cmt_crv, cmt_swap);
                newton(yans, niter, a, b, &nstop);
                count += 1;
            }

            *value = a[2];
        } /* END if (cmt_code == TEC10) */

        break;
    }
    free_Swap(cmt_swap);

    return NULL;
}

/* ====================================================================== */

#undef MAXITER
