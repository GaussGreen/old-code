/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_CMT_FWDTREAS.C                                  */
/*                                                                            */
/*      PURPOSE:        Return a fwd treas rate from a CMT strip    	      */
/*                                                                            */
/*      AUTHORS:        Olivier VAN EYSEREN                    		      */
/*                                                                            */
/*      DATE:           23rd January 1996                                     */
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
/*      	<none>                                                        */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     	                                              */
/*                                                                            */
/*      DATE:                          			                      */
/*                                                                            */
/*                                                                            */
/*      REASON:                                                               */
/*                                                                            */
/*                                                                            */
/******************************************************************************/

/**  Include Statements  ******************************************************/

#include "float.h"
#include "math.h"
#include "swp_h_all.h"
#include "swp_h_cmt_fwdtreas.h"
#include "swp_h_cmt_init.h"
#include "swp_h_cmt_param.h"
#include "swp_h_cmt_types.h"

/* ======================================================================== */
/* When calling this function crv_list->current->element has to be
   a CMT market */

Err fwd_treas_rate(Ddate start, SrtCurvePtr tc_crv, double* fwd_T_rate)
{
    int         interp_method, num_rate;
    int         i, i1, i2;
    double      spread, swap_rate;
    double      fwd_treas_high, fwd_treas_low;
    double      spread_high, spread_low;
    double      fwd_treas;
    double      time_start, time_low, time_high;
    Date        today;
    String      crv_name, yc_name;
    FwdT_Obj*   tc_obj;
    SRT_Boolean use_spread_interpolation;
    CMT_Param*  cmt;
    SwapDP      swapdp;
    Err         err = NULL;

    tc_obj = get_fwdTobj_from_cmtcrv(tc_crv);
    if (!tc_obj)
        return ("Can't find specified treasury curve object!");

    crv_name = get_curve_name(tc_crv);

    cmt = get_cmtparam_from_cmtcrv(tc_crv);
    if (!cmt)
        return ("Can't find specified CMT parameters!");

    use_spread_interpolation = cmt->interp_on_spread_flg;

    yc_name  = get_ycname_from_cmtcrv(tc_crv);
    today    = (Date)FwdT_Dget(tc_obj, FwdT_TODAY);
    num_rate = FwdT_Iget(tc_obj, FwdT_NUMSPREAD);

    if (use_spread_interpolation == SRT_YES)
    {
        swapdp.start             = (Date)start;
        swapdp.first_full_fixing = swapdp.start;
        swapdp.nfp               = (int)(cmt->cmt_mat) * cmt->swap_compd;
        swapdp.end               = swapdp.nfp;
        swapdp.direction         = FWD;
        swapdp.compd             = cmt->swap_compd;
        swapdp.basis_code        = cmt->swap_basis_code;
        swapdp.spot_lag          = cmt->spot_lag;

        /* Need the currency as reference rate */
        err = swp_f_ForwardRate_SwapDP(&swapdp, yc_name, cmt->swap_ref_rate, &swap_rate);
        if (err)
            return err;
    }

    interp_method = FwdT_Iget(tc_obj, FwdT_INTERP_METHOD);

    if (num_rate == 1)
    {
        if (use_spread_interpolation == SRT_YES)
        {
            spread    = FwdT_field_Dget(tc_obj, FwdT_SPREAD, 0);
            fwd_treas = swap_rate + spread;
        }
        else if (use_spread_interpolation == SRT_NO)
        {
            fwd_treas = FwdT_field_Dget(tc_obj, FwdT_RATE, 0);
        }
    } /* END if (num_rate == 1) */
    else
    {
        /* COMPUTE THE ENDPOINTS TO INTERPOLATE THE FORWARD RATES   */
        i = 0;
        while ((DTOL(FwdT_field_Dget(tc_obj, FwdT_DATE, i)) < DTOL(start)) && (i < num_rate))
            i++;

        i1 = i - 1;
        i2 = i;

        if (i == num_rate)
        {
            if (use_spread_interpolation == SRT_YES)
            {
                spread    = FwdT_field_Dget(tc_obj, FwdT_SPREAD, i1);
                fwd_treas = swap_rate + spread;
            }
            else if (use_spread_interpolation == SRT_NO)
            {
                fwd_treas = FwdT_field_Dget(tc_obj, FwdT_RATE, i1);
            }
        }
        else if (i2 == 0)
        {
            if (use_spread_interpolation == SRT_YES)
            {
                spread    = FwdT_field_Dget(tc_obj, FwdT_SPREAD, 0);
                fwd_treas = swap_rate + spread;
            }
            else if (use_spread_interpolation == SRT_NO)
            {
                fwd_treas = FwdT_field_Dget(tc_obj, FwdT_RATE, 0);
            }
        }
        else
        {
            time_start = (start - today) * YEARS_IN_DAY;
            time_low   = (FwdT_field_Dget(tc_obj, FwdT_DATE, i1) - today) * YEARS_IN_DAY;
            time_high  = (FwdT_field_Dget(tc_obj, FwdT_DATE, i2) - today) * YEARS_IN_DAY;

            if (use_spread_interpolation == SRT_YES)
            {
                spread_low  = FwdT_field_Dget(tc_obj, FwdT_SPREAD, i1);
                spread_high = FwdT_field_Dget(tc_obj, FwdT_SPREAD, i2);
            }
            else if (use_spread_interpolation == SRT_NO)
            {
                fwd_treas_low  = FwdT_field_Dget(tc_obj, FwdT_RATE, i1);
                fwd_treas_high = FwdT_field_Dget(tc_obj, FwdT_RATE, i2);
            }

            /* THIS IS LINEAR INTERPOLATION IN R*/
            if (interp_method == LIN_R)
            {
                if (use_spread_interpolation == SRT_YES)
                {
                    spread = spread_low + (spread_high - spread_low) * (time_start - time_low) /
                                              (time_high - time_low);

                    fwd_treas = swap_rate + spread;
                }
                else if (use_spread_interpolation == SRT_NO)
                {
                    fwd_treas = fwd_treas_low + (fwd_treas_high - fwd_treas_low) *
                                                    (time_start - time_low) /
                                                    (time_high - time_low);
                }
            }
            else
                /* THIS IS LINEAR INTERPOLATION IN R_T*/
                if (interp_method == LIN_RT)
            {
                if (use_spread_interpolation == SRT_YES)
                {
                    spread =
                        spread_low * time_low + ((spread_high * time_high - spread_low * time_low) *
                                                 (time_start - time_low) / (time_high - time_low));
                    spread /= time_start;

                    fwd_treas = swap_rate + spread;
                }
                else if (use_spread_interpolation == SRT_NO)
                {
                    fwd_treas = fwd_treas_low * time_low +
                                ((fwd_treas_high * time_high - fwd_treas_low * time_low) *
                                 (time_start - time_low) / (time_high - time_low));
                    fwd_treas /= time_start;
                }
            } /* END if (interp_method == LIN_RT) */
        }     /* END if (i2 != 0) */
    }         /* END if (num_rate != 1) */

    /* Make sure that crv_list->current points to the CMT makret*/
    tc_crv = lookup_curve(crv_name);

    *fwd_T_rate = fwd_treas;

    return (NULL);
}