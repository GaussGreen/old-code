/**********************************************************************
 *      Name: SwpSwapCalc.c                                         *
 *  Function: Performs bootstrap calibration                          *
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
#include "Swp_h_all.h"

char* SwpSwapCalc(
    long    value_date,
    long    start,
    long    nfp_or_end,
    long    today,
    char*   cpdStr,
    char*   basisStr,
    double  strike,
    double  initial_not,
    double  final_not,
    char*   info_message,
    char*   dcntName,
    double* answer)
{
    Err         err;
    double      swap_res_target = 0.0;
    Arg_Obj*    arg;
    int         dir;
    int         compd;
    int         basis_code;
    SrtCurvePtr crv;
    Message     info_type;

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

    strupper(cpdStr);
    strupper(basisStr);

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

    if (value_date < today || start < today)
    {
        return serror("Value date and start date must be after today");
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

    /* Transfer strings into Messages */

    if (err = interp_swap_message(info_message, &info_type))
    {
        return err;
    }

    /*  Initialize arg calls  */

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
    Arg_Dset(arg, ARG_INITIAL_NOT, initial_not);
    Arg_Dset(arg, ARG_FINAL_NOT, final_not);
    Arg_set(arg, ARG_INFO_TYPE, info_type);
    Arg_set(arg, ARG_RESULT_TARGET, &swap_res_target);

    if (err = SWAP_compute(crv, info_type, arg))
    {
        return err;
    }

    free_Arg(arg);

    *answer = swap_res_target;

    return 0;
}
