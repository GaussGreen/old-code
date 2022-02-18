/**********************************************************************
 *      Name: SwpSwapUnwind.c                                         *
 *  Function: Swap_unwind EXCEL addin                                 *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 12/12/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 12/12/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "SwpAccess.h"
#include "swp_h_all.h"

char* SwpSwapUnwind(
    long    value_date,
    long    start,
    long    nfp_or_end,
    long    today,
    char*   cpdStr,
    char*   basisStr,
    double  strike,
    char*   first_short_full,
    double  initial_not,
    double  final_not,
    char*   liborCpdStr,
    char*   liborBasStr,
    double  libor_fix,
    char*   dcntName,
    double* pv)
{
    Err         err;
    double      swap_res_target = 0.0;
    double      tmp_res;
    Arg_Obj*    arg;
    int         dir;
    int         compd;
    int         basis_code;
    int         endofday;
    int         libor_term;
    int         libor_basis;
    SRT_Boolean b;
    SrtCurvePtr crv;

    /* pass end of day as optional argument ? */
    /* It can be 1 or 0 */

    endofday = 1;

    if ((0 < nfp_or_end) && (nfp_or_end < 500))
    {
        dir = FWD;
    }
    else if (nfp_or_end < 30317)
    {
        return serror("end date too long ago");
    }
    else
    {
        dir = BKWD;
    }

    /* check start date against today */

    if ((dir == BKWD) && (start >= nfp_or_end))
    {
        return serror("end date must be after start date");
    }

    /* interpret fixed leg basis and compounding */

    strupper(basisStr);
    strupper(cpdStr);
    strupper(liborBasStr);
    strupper(liborCpdStr);

    err = interp_basis(basisStr, &basis_code);
    if (err)
    {
        return err;
    }
    err = interp_compounding(cpdStr, &compd);
    if (err)
    {
        return err;
    }

    /* interpret floating leg basis and compounding */

    err = interp_basis(liborBasStr, &libor_basis);
    if (err)
    {
        return err;
    }
    err = interp_compounding(liborCpdStr, &libor_term);
    if (err)
    {
        return err;
    }

    strupper(first_short_full);
    strip_white_space(first_short_full);

    if (!strcmp(first_short_full, "FULL"))
    {
        b = SRT_YES;
    }
    else if (!strcmp(first_short_full, "F"))
    {
        b = SRT_YES;
    }
    else if (!strcmp(first_short_full, "SHORT"))
    {
        b = SRT_NO;
    }
    else if (!strcmp(first_short_full, "S"))
    {
        b = SRT_NO;
    }
    else
    {
        return serror("invalid first_short_full");
    }

    if (value_date < today)
    {
        return serror("Value date must be after today");
    }

    if (dcntName)
    {
        crv = lookup_curve(dcntName);
    }
    else
    {
        return serror("Heinous error: No Underlying passed to SwpSwapUnwind");
    }

    if (!crv)
    {
        return serror("Could not find discount curve in market list");
    }

    /* calculate fixed leg */
    if (strike != 0.0)
    {
        /**  Initialize arg calls  **/
        init_Arg(&arg, (int)COMPUTE_DISC_FACTOR);
        Arg_Dateset(arg, ARG_TODAY, today);
        Arg_Dateset(arg, ARG_VALUE_DATE, value_date);
        Arg_Dateset(arg, ARG_START, start);
        Arg_Dateset(arg, ARG_FIRST_FULL_FIXING, start);
        Arg_set(arg, ARG_NUM_FULL_PERIOD, nfp_or_end);
        Arg_set(arg, ARG_END, nfp_or_end);
        Arg_set(arg, ARG_DATE_DIR, dir);
        Arg_set(arg, ARG_COMPD, compd);
        Arg_set(arg, ARG_BASIS_CODE, basis_code);
        Arg_Dset(arg, ARG_STRIKE, strike / 100.0);
        Arg_Dset(arg, ARG_INITIAL_NOT, 0.0);
        Arg_Dset(arg, ARG_FINAL_NOT, 0.0);
        Arg_set(arg, ARG_INFO_TYPE, COMPUTE_PV);
        Arg_set(arg, ARG_RESULT_TARGET, &tmp_res);

        if (err = swap_unwind_compute(crv, arg, b, endofday))
        {
            return err;
        }

        swap_res_target = tmp_res;
        free_Arg(arg);
    }

    /* calculate floating leg */
    if (initial_not != 0.0 || final_not != 0.0)
    {
        if (dir == FWD)
        {
            nfp_or_end = nfp_or_end * libor_term / compd;
        }

        /**  Initialize arg calls  **/
        init_Arg(&arg, (int)COMPUTE_DISC_FACTOR);
        Arg_Dateset(arg, ARG_TODAY, today);
        Arg_Dateset(arg, ARG_VALUE_DATE, value_date);
        Arg_Dateset(arg, ARG_START, start);
        Arg_Dateset(arg, ARG_FIRST_FULL_FIXING, start);
        Arg_set(arg, ARG_NUM_FULL_PERIOD, nfp_or_end);
        Arg_set(arg, ARG_END, nfp_or_end);
        Arg_set(arg, ARG_DATE_DIR, dir);
        Arg_set(arg, ARG_COMPD, libor_term);
        Arg_set(arg, ARG_BASIS_CODE, libor_basis);
        Arg_Dset(arg, ARG_STRIKE, libor_fix / 100.0);
        Arg_Dset(arg, ARG_INITIAL_NOT, initial_not);
        Arg_Dset(arg, ARG_FINAL_NOT, final_not);
        Arg_set(arg, ARG_INFO_TYPE, COMPUTE_PV);
        Arg_set(arg, ARG_RESULT_TARGET, &tmp_res);

        if (err = swap_unwind_compute(crv, arg, b, endofday))
        {
            return err;
        }
        swap_res_target += tmp_res;
        free_Arg(arg);
    }

    *pv = swap_res_target;

    return 0;
}
