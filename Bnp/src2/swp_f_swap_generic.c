/* =============================================================================

   FILENAME:    swp_f_swap_generic.c

   PURPOSE:     Provide a few function for generic swaps, including:
                                - Fixed leg
                                - Floating leg
                                - Spreads
                                - CMS, CMT, TEC
                                - ...
                                THE FOLLOWING STRUCTURES AND FUNCTIONS SEEMS MORE ROBUST
                                (OR AT LEAST EASIER TO WORK WITH) AND SHOULD BE USED IN THE FUTURE)

   ============================================================================= */
#include "math.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_curve_fct.h"
#include "swp_h_swap_generic.h"
#include "swp_h_swap_structures.h"

/* ------------------------------------------------------------------------------ */
/* Creates the swap leg and sets everything that is stripper independent */
/* Initial and final exchange are assumed positive:
      initial will be accounted for as a +,
          final   will be accounted for as a -  ( we receive the floating cash) */
/* crv is needed only for CMS and CMT */

static Err create_GenSwapLeg(
    long            today,
    SwapDP*         leg_sdp,
    double          fixed,
    double          initial_exchange,
    double          final_exchange,
    LegType         leg_type,
    SrtReceiverType rec_pay,
    String          curve_name,
    String          ref_rate_code,
    GenSwapLeg*     leg)
{
    Err         err = NULL;
    SrtCurvePtr crv;
    int         i;

    /* Sets everything inside leg to NULL*/
    memset(leg, 0, sizeof(GenSwapLeg));

    /* Sets leg type */
    leg->leg_type = leg_type;

    /* Sets leg today */
    leg->today = today;

    /* Sets start, end, comp, basis ... */
    leg->sdp = *leg_sdp;

    /* Allocate space for and generate the dates of the swap */
    leg->pay_date   = SwapDP_to_DateList(leg_sdp, MODIFIED_SUCCEEDING);
    leg->start_date = new_DateList(leg->pay_date.len - 1);
    leg->end_date   = new_DateList(leg->pay_date.len - 1);
    for (i = 0; i < leg->pay_date.len - 1; i++)
    {
        leg->start_date.date[i] = leg->pay_date.date[i];
        leg->end_date.date[i]   = leg->pay_date.date[i + 1];
    }

    /* Allocate space for and generate the fixing dates of the swap */
    err = fixing_list(leg->start_date, leg_sdp->spot_lag, &(leg->fixing_date));
    if (err)
        return (err);

    /* Sets leg_length : number of dates generated */
    leg->leg_length = leg->pay_date.len;

    /* Allocate space for and generate the times (from today to swap dates) of the swap */
    err = time_list(leg->pay_date, today, &(leg->time));
    if (err)
        return (err);

    /* Allocate space for and generate the coverages of the swap
       (the coverage is stroed on the last date of the period) */
    err = cvg_list(leg->start_date, leg->end_date, leg_sdp->basis_code, &(leg->cvg));
    if (err)
        return (err);

    switch (leg_type)
    {
    case NOTIONALS_LEG:
        leg->initial_exchange = initial_exchange;
        leg->final_exchange   = final_exchange;
        break;
    case FIXED_LEG:
        err = const_list(leg->end_date, fixed, &(leg->cpn));
        if (err)
            return (err);

        break;
    case FIXED_AND_NOTIONALS_LEG:
        leg->initial_exchange = initial_exchange;
        leg->final_exchange   = final_exchange;
        err                   = const_list(leg->end_date, fixed, &(leg->cpn));
        if (err)
            return (err);
        break;
    case SPREAD_AND_NOTIONALS_LEG:
        /* Sets the spreads RefRate (Libor,...) vs Cash + exchange of notionals */
        leg->initial_exchange = initial_exchange;
        leg->final_exchange   = final_exchange;
        err = spread_list(leg->start_date, leg->end_date, ref_rate_code, &(leg->spread));
        if (err)
            return (err);

        break;
    case SPREAD_LEG:
        /* Sets the spreads RefRate (Libor,...) vs Cash */
        err = spread_list(leg->start_date, leg->end_date, ref_rate_code, &(leg->spread));
        if (err)
            return (err);
        if ((leg->spread).len == 0)
            return serror("Error in computing spreads in swap leg");
        break;
    case FLOATING_LEG:
        /* Sets the spreads RefRate (Libor,...) vs Cash + cash_fra + exchange of notionals */
        err = spread_list(leg->start_date, leg->end_date, ref_rate_code, &(leg->spread));
        if (err)
            return (err);
        if ((leg->spread).len == 0)
            return serror("Error in computing spreads in swap leg");
        err = fwd_cash_list(
            leg->start_date, leg->end_date, leg_sdp->basis_code, curve_name, &(leg->fwd));
        break;
    case DRS_LEG:
        break;
    case CMS_MARGIN_LEG:
        /* Get the curve */
        crv = lookup_curve(curve_name);

        /* Sets the swaption volatilities and the forward swap rates */
        err = fwd_swap_list(leg->pay_date, crv, &(leg->fwd));
        if (err)
            return (err);
        err = swap_vol_list(leg->pay_date, leg->fwd, crv, &(leg->vol));
        if (err)
            return (err);
        break;
    case CMT_MARGIN_LEG:
        /* Get the curve */
        crv = lookup_curve(curve_name);

        /* Sets the swap data: swaption vols and forward swap rates */
        err = fwd_swap_list(leg->pay_date, crv, &(leg->temp_fwd));
        if (err)
            return (err);
        err = swap_vol_list(leg->pay_date, leg->temp_fwd, crv, &(leg->vol));
        if (err)
            return (err);
        break;
    case TEC_MARGIN_LEG:
        /* Get the curve */
        crv = lookup_curve(curve_name);

        /* Sets the swap data : swaption vols and forward swap rate */
        err = fwd_swap_list(leg->pay_date, crv, &(leg->temp_fwd));
        if (err)
            return (err);
        err = swap_vol_list(leg->pay_date, leg->temp_fwd, crv, &(leg->vol));
        if (err)
            return (err);
        break;
    }

    /* Receive or pay the leg */
    leg->rec_pay = rec_pay;

    /* Return a success message */
    return err;

} /* END GenSwapLeg create_GenSwapLeg(...) */

/* ---------------------------------------------------------------------------------- */

Err swp_f_freein_GenSwapLeg(GenSwapLeg* leg)
{
    Err err = NULL;

    if (leg->fixing_date.date)
        srt_free(leg->fixing_date.date);
    if (leg->start_date.date)
        srt_free(leg->start_date.date);
    if (leg->end_date.date)
        srt_free(leg->end_date.date);
    if (leg->pay_date.date)
        srt_free(leg->pay_date.date);
    if (leg->time.d)
        free_inDlist(&leg->time);
    if (leg->mat.d)
        free_inDlist(&leg->mat);
    if (leg->fwd.d)
        free_inDlist(&leg->fwd);
    if (leg->temp_fwd.d)
        free_inDlist(&leg->temp_fwd);
    if (leg->vol.d)
        free_inDlist(&leg->vol);
    if (leg->cpn.d)
        free_inDlist(&leg->cpn);
    if (leg->cvg.d)
        free_inDlist(&leg->cvg);
    if (leg->payment.d)
        free_inDlist(&leg->payment);
    if (leg->df.d)
        free_inDlist(&leg->df);
    if (leg->spread.d)
        free_inDlist(&leg->spread);

    leg = NULL;

    return err;

} /* END freein_GenSwapLeg(...) */

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */

Err swp_f_create_GenFullSwap(
    SwapDP*      rec_leg_sdp,
    SwapDP*      pay_leg_sdp,
    double       fixed,
    double       initial_exchange,
    double       final_exchange,
    SwapType     swap_type,
    String       curve_name,
    GenFullSwap* swap)
{
    Err         err;
    GenFullSwap full_swap;
    long        today;
    long        spot;
    SrtCurvePtr crv;
    String      yc_name;

    /* Gets the curve */
    if ((swap_type == CMS_MARGIN_FLOATING_SWAP) || (swap_type == CMS_FOR_TEC_MARGIN_FIXED_SWAP) ||
        (swap_type == CMT_MARGIN_FLOATING_SWAP) || (swap_type == CMS_MARGIN_FIXED_SWAP) ||
        (swap_type == CMS_FOR_TEC_MARGIN_FLOATING_SWAP) || (swap_type == CMT_MARGIN_FIXED_SWAP) ||
        (swap_type == TEC_MARGIN_FLOATING_SWAP) || (swap_type == TEC_MARGIN_FIXED_SWAP))
    {
        crv     = lookup_curve(curve_name);
        yc_name = get_ycname_from_cmtcrv(crv);
    }
    else
    {
        yc_name = curve_name;
    }
    crv = lookup_curve(yc_name);

    /* Sets everything inside swap to NULL*/
    memset(&full_swap, 0, sizeof(GenSwapLeg));

    /* Gets today and spotdate from the underlying */
    today              = get_clcndate_from_yldcrv(crv);
    spot               = get_spotdate_from_yldcrv(crv);
    full_swap.today    = today;
    full_swap.spotdate = spot;

    /* Sets notional */
    full_swap.notional = DEFAULT_NOTIONAL;

    /* Sets swap type */
    full_swap.swap_type = swap_type;

    /* According to the swap type, sets the legs accordingly */
    switch (swap_type)
    {
        /* Fixed leg vs exchange of notionals (initial & final) leg*/
    case FIXED_NOTIONALS_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            fixed,
            0.0,
            0.0,
            FIXED_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            NOTIONALS_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* Fixed leg vs floating (with cash spreads) leg */
    case FIXED_FLOATING_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            fixed,
            0.0,
            0.0,
            FIXED_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            FLOATING_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* DRS leg minus margin vs floating (with cash spreads) */
    case DRS_MARGIN_FLOATING_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            DRS_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            FLOATING_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* DRS leg minus margin vs fixed leg */
    case DRS_MARGIN_FIXED_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            DRS_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            FIXED_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* Incomplete swaps: only one leg
                /* Fixed leg and (both) notionals in one leg
                        ONE_LEG_FIXED_AND_NOTIONALS,
                /* Fixed leg and final notional in one leg
                        ONE_LEG_STANDARD_BOND,
                /* Fixed leg only in one leg (for level payment)
                        ONE_LEG_FIXED_SWAP,

                /* CMT, CMS legs...
                        ONE_LEG_CMS_MARGIN,
                        ONE_LEG_CMT_MARGIN,
                        ONE_LEG_CMS_FOR_TEC_MARGIN,
                        ONE_LEG_TEC_MARGIN,

                        ONE_LEG_DRS
        */

        /* CMS minus margin leg vs Libor leg(with cash spreads) */
    case CMS_MARGIN_FLOATING_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            CMS_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            NOTIONALS_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* CMS minus margin (the TEC way) leg vs fixed leg */
    case CMS_FOR_TEC_MARGIN_FLOATING_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            CMS_FOR_TEC_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            NOTIONALS_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* CMS minus margin (the TEC way) leg vs fixed leg */
    case CMS_FOR_TEC_MARGIN_FIXED_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            CMS_FOR_TEC_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            fixed,
            1.0,
            1.0,
            FIXED_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* CMT minus margin leg vs Libor leg (with cash spreads) */
    case CMT_MARGIN_FLOATING_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            CMT_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            NOTIONALS_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* TEC minus margin (the TEC way) leg vs Libor leg (with cash spreads) */
    case TEC_MARGIN_FLOATING_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            TEC_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            0.0,
            1.0,
            1.0,
            NOTIONALS_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* CMS minus margin leg vs fixed leg */
    case CMS_MARGIN_FIXED_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            CMS_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            fixed,
            0.0,
            0.0,
            FIXED_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* CMT minus margin leg vs fixed leg */
    case CMT_MARGIN_FIXED_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            CMT_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            fixed,
            0.0,
            0.0,
            FIXED_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;

        /* TEC minus margin (the TEC way) leg vs fixed leg */
    case TEC_MARGIN_FIXED_SWAP:
        err = create_GenSwapLeg(
            today,
            rec_leg_sdp,
            0.0,
            0.0,
            0.0,
            TEC_MARGIN_LEG,
            SRT_RECEIVER,
            curve_name,
            NULL,
            &full_swap.leg[0]);
        if (err)
            return (err);
        err = create_GenSwapLeg(
            today,
            pay_leg_sdp,
            fixed,
            0.0,
            0.0,
            FIXED_LEG,
            SRT_PAYER,
            curve_name,
            NULL,
            &full_swap.leg[1]);
        if (err)
            return (err);
        break;
    default:
        memset(&full_swap, 0, sizeof(GenFullSwap));
    }

    /* Sets the initialised structure*/
    *swap = full_swap;

    /* return a success message*/
    return err;

} /* END swp_f_create_GenFullSwap(...) */

/* ------------------------------------------------------------------------ */

Err swp_f_freein_GenFullSwap(GenFullSwap* swap)
{
    Err err = NULL;

    err = swp_f_freein_GenSwapLeg(&swap->leg[0]);
    if (err)
        return err;
    err = swp_f_freein_GenSwapLeg(&swap->leg[1]);
    if (err)
        return err;

    swap = NULL;

    return err;
} /* END swp_f_freein_GenFullSwap(...) */

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------- */
/* Here, everything that could be SrtCurvePtr dependent (for a stripper)
   is set up (sets payments)
   We assume we receive the leg
   The SrtCurvePtr is used only for CMS, CMT...
*/
static Err populate_GenSwapLeg(double margin, SrtCurvePtr crv, GenSwapLeg* leg)
{
    int i;
    Err err = NULL;

    /* Allocate space for leg->payments*/
    err = new_Dlist(leg->pay_date.len, &(leg->payment));
    if (err)
        return (err);

    /* Sets all the payments on the leg depending on the type of leg */
    switch (leg->leg_type)
    {
    case NOTIONALS_LEG:
        leg->payment.d[0]                     = -leg->initial_exchange;
        leg->payment.d[leg->pay_date.len - 1] = leg->final_exchange;
        break;
    case FIXED_AND_NOTIONALS_LEG:
        leg->payment.d[0]                     = -leg->initial_exchange;
        leg->payment.d[leg->pay_date.len - 1] = leg->final_exchange;
        for (i = 0; i < leg->end_date.len; i++)
        {
            leg->payment.d[i + 1] += leg->cpn.d[i] * leg->cvg.d[i];
        }
        break;
    case FIXED_LEG:
        for (i = 0; i < leg->end_date.len; i++)
        {
            leg->payment.d[i + 1] = leg->cpn.d[i] * leg->cvg.d[i];
        }
        break;
    case SPREAD_LEG:
        for (i = 0; i < leg->end_date.len; i++)
        {
            leg->payment.d[i + 1] += leg->spread.d[i] * leg->cvg.d[i];
        }
        break;
    case SPREAD_AND_NOTIONALS_LEG:
        leg->payment.d[0]                     = leg->initial_exchange;
        leg->payment.d[leg->pay_date.len - 1] = -leg->final_exchange;
        for (i = 0; i < leg->end_date.len; i++)
        {
            leg->payment.d[i + 1] += leg->spread.d[i] * leg->cvg.d[i];
        }
        break;
    case FLOATING_LEG:
        for (i = 0; i < leg->end_date.len; i++)
        {
            leg->payment.d[i + 1] += (leg->spread.d[i] + leg->fwd.d[i]) * leg->cvg.d[i];
        }
        break;
    case CMS_MARGIN_LEG:
        leg->margin = margin;
        err         = conv_fwd_list(
            leg->today, leg->pay_date, leg->fwd, leg->vol, crv, SWAP_RATE, &(leg->cpn));
        if (err)
            return (err);

        for (i = 1; i < leg->pay_date.len; i++)
        {
            leg->payment.d[i] = (leg->cpn.d[i] + leg->margin) * leg->cvg.d[i - 1];
        }
        break;
    case CMT_MARGIN_LEG:
        err = fwd_treas_list(leg->pay_date, crv, &(leg->fwd));
        if (err)
            return (err);
        err = treas_vol_list(leg->pay_date, leg->temp_fwd, leg->fwd, crv, &(leg->vol));
        if (err)
            return (err);
        err = conv_fwd_list(
            leg->today, leg->pay_date, leg->fwd, leg->vol, crv, TREASURY_RATE, &(leg->cpn));
        if (err)
            return (err);
        for (i = 1; i < leg->pay_date.len; i++)
        {
            leg->payment.d[i] = leg->cpn.d[i] * leg->cvg.d[i - 1];
        }
        break;
    case CMS_FOR_TEC_MARGIN_LEG:
        leg->margin = margin;
        err         = conv_fwd_list(
            leg->today, leg->pay_date, leg->fwd, leg->vol, crv, SWAP_RATE, &(leg->cpn));
        if (err)
            return (err);
        for (i = 1; i < leg->pay_date.len; i++)
        {
            leg->payment.d[i] =
                exp(1.0 / (double)(leg->sdp.compd) * log(1.0 + leg->cpn.d[i])) - 1.0;
        }
        break;
    case TEC_MARGIN_LEG:
        err = fwd_treas_list(leg->pay_date, crv, &(leg->fwd));
        if (err)
            return (err);
        err = treas_vol_list(leg->pay_date, leg->temp_fwd, leg->fwd, crv, &(leg->vol));
        if (err)
            return (err);
        err = conv_fwd_list(
            leg->today, leg->pay_date, leg->fwd, leg->vol, crv, TREASURY_RATE, &(leg->cpn));
        if (err)
            return (err);
        for (i = 1; i < leg->pay_date.len; i++)
        {
            leg->payment.d[i] =
                exp(1.0 / (double)(leg->sdp.compd) * log(1.0 + leg->cpn.d[i])) - 1.0;
        }
        break;
    }

    /* Return a success message */
    return NULL;

} /* END populate_GenSwapLeg(...) */

/* ------------------------------------------------------------------------ */
/* According to the swap type, sets what is loop/stripper/crv dependent */

Err swp_f_populate_GenFullSwap(double margin, SrtCurvePtr crv, GenFullSwap* swap)
{
    Err err = NULL;
    int i;

    for (i = 0; i < 2; i++)
    {
        err = populate_GenSwapLeg(margin, crv, &swap->leg[i]);
        if (err)
            return err;
    }

    /* Return a success message */
    return NULL;

} /* END swp_f_populate_GenFullSwap(...) */

/* ------------------------------------------------------------------------- */
/* Sets discount factors and values the leg */

static Err value_GenSwapLeg(GenSwapLeg* leg, SrtCurvePtr crv)
{
    int    i;
    String crv_name;
    String yc_name;
    double sign;
    double pv;
    Err    err = NULL;

    /* Set the leg value to zero by default */
    leg->leg_value = 0.0;

    /* If there is no leg: the value is zero */
    if (leg->leg_length == 0)
        return NULL;

    /* Sets signs depending on receiver or payer */
    if (leg->rec_pay == SRT_RECEIVER)
        sign = +1.0;
    else if (leg->rec_pay == SRT_PAYER)
        sign = -1.0;

    /* Gets market and yc names */
    crv_name = get_curve_name(crv);
    yc_name  = get_ycname_from_curve(crv);
    crv      = lookup_curve(yc_name);

    /* Sets discount factors (needs a yc_crv)*/
    err = df_list(leg->pay_date, yc_name, &(leg->df));
    if (err)
        return (err);

    /* Goes through the payments (from first date to last one included)*/
    pv = 0.0;
    for (i = 0; i < leg->pay_date.len; i++)
    {
        pv += leg->df.d[i] * leg->payment.d[i];
    }

    /* Sets the market to its initial value */
    crv = lookup_curve(crv_name);

    /* Sets the value of the leg with the right sign*/
    pv *= sign;
    leg->leg_value = pv;

    /* Return a success message */
    return NULL;

} /* END value_GenSwapLeg(...) */

/* ------------------------------------------------------------------------- */

Err swp_f_value_GenFullSwap(GenFullSwap* swap, SrtCurvePtr crv)
{
    int i;
    Err err;

    /* Sets swap value at 0 initially */
    swap->swap_value = 0.0;

    /* Compute the value of each leg */
    for (i = 0; i < 2; i++)
    {
        err = value_GenSwapLeg(&swap->leg[i], crv);
        if (err)
            return err;
        swap->swap_value += swap->leg[i].leg_value;
    }

    /* Multiply the value by the notional */
    swap->swap_value *= swap->notional;

    /* Returns a success message */
    return NULL;

} /* END swp_f_value_GenFullSwap(...) */

/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/* THESE FUNCTIONS ARE VERY USEFUL, AND EASIER TO USE IN SIMPLE CASES        */
/* ------------------------------------------------------------------------- */

/* Make the dates and coupons for a fixed leg + populates the leg */
Err swp_f_make_FixedLeg(SwapDP* fixedsdp, double fixed, Date today, GenSwapLeg* fixedleg)
{
    Err err = NULL;

    /* Create the fixed leg: dates, cvg, coupon */
    err = create_GenSwapLeg(
        today, fixedsdp, fixed, 0.0, 0.0, FIXED_LEG, SRT_RECEIVER, NULL, NULL, fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(fixedleg);
        return err;
    }

    /* Sets the payments in the fixed leg (no margin): payment = coupon * cvg */
    err = populate_GenSwapLeg(0.0, NULL, fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(fixedleg);
        return err;
    }

    /* Return a success message */
    return NULL;

} /* END make_FixedLeg(...) */

/* ------------------------------------------------------------------------- */

/* Make the dates and coupons for a fixed leg + notional (cash flow at spot) */
Err swp_f_make_FixedAndNotionalsLeg(
    SwapDP* fixedsdp, double fixed, double initial, double final, Date today, GenSwapLeg* fixedleg)
{
    Err err = NULL;

    /* Create the fixed and notionals leg: dates, cvg, coupon, initial, final */
    err = create_GenSwapLeg(
        today,
        fixedsdp,
        fixed,
        initial,
        final,
        FIXED_AND_NOTIONALS_LEG,
        SRT_RECEIVER,
        NULL,
        NULL,
        fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(fixedleg);
        return err;
    }

    /* Sets the payments in the fixed & notionals leg (no margin): payment = coupon*cvg */
    err = populate_GenSwapLeg(0.0, NULL, fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(fixedleg);
        return err;
    }

    /* Return a success message */
    return NULL;

} /* END make_FixedAndNotionalsLeg(...) */

Err swp_f_make_NotionalsLeg(
    SwapDP* fixedsdp, double initial, double final, Date today, GenSwapLeg* fixedleg)
{
    Err err = NULL;

    /* Create the notionals leg: initial, final */
    err = create_GenSwapLeg(
        today, fixedsdp, 0.0, initial, final, NOTIONALS_LEG, SRT_RECEIVER, NULL, NULL, fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(fixedleg);
        return err;
    }

    /* Sets the payments in the fixed & notionals leg (no margin): payment = coupon*cvg */
    err = populate_GenSwapLeg(0.0, NULL, fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(fixedleg);
        return err;
    }

    /* Return a success message */
    return NULL;
}

/* ------------------------------------------------------------------------- */

/* Make the dates and spreads + fra's for a floating leg + sets payments */
Err swp_f_make_FloatingLeg(
    long start, long end, Date today, String ref_rate_code, SrtCurvePtr crv, GenSwapLeg* floatleg)
{
    Err            err = NULL;
    SwapDP         floatsdp;
    SrtCompounding comp;
    BasisCode      basis;

    /* Get the details of the floating reference rate: compounding and basis */
    err = swp_f_get_ref_rate_details(ref_rate_code, &basis, &comp);
    if (err)
    {
        return err;
    }

    /* Set up the corresponding SwapDP: start, end,... generated BACKWARD */
    err                = swp_f_setSwapDP(start, end, comp, basis, &floatsdp);
    floatsdp.direction = BKWD;
    if (err)
        return err;

    /* Get the spot lag from the Yield Curve and sets it into the SwapDP */
    floatsdp.spot_lag = get_spotlag_from_curve(crv);

    /* Create the floating leg: dates, cvg, spreads fra's */
    err = create_GenSwapLeg(
        today,
        &floatsdp,
        0.0,
        0.0,
        0.0,
        FLOATING_LEG,
        SRT_RECEIVER,
        get_curve_name(crv),
        ref_rate_code,
        floatleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(floatleg);
        return err;
    }

    /* Sets the payments in the floating leg (no margin): payment = (spread + fra)*cvg */
    err = populate_GenSwapLeg(0.0, crv, floatleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(floatleg);
        return err;
    }

    /* Return a success message */
    return NULL;

} /* END swp_f_make_FloatingLeg(...) */

/* ------------------------------------------------------------------------- */
/* Make the dates and spreads (no notionals) for a floating leg (no fra's) */
Err swp_f_make_SpreadLeg(
    long start, long end, Date today, String ref_rate_code, GenSwapLeg* spreadleg)
{
    Err            err = NULL;
    SwapDP         spreadsdp;
    int            spotlag;
    SrtCompounding comp;
    BasisCode      basis;

    /* Get the details of the floating reference rate: compounding and basis */
    err = swp_f_get_ref_rate_details(ref_rate_code, &basis, &comp);
    if (err)
        return err;

    /* Set up the corresponding SwapDP: start, end,... generated BACKWARD for theo end date */
    err                 = swp_f_setSwapDP(start, end, comp, basis, &spreadsdp);
    spreadsdp.direction = BKWD;
    if (err)
        return err;

    /* Get the spot lag from the ref rate and attach it to the SwapDp */
    err = srt_f_get_spot_lag_from_refrate(ref_rate_code, &spotlag);
    if (err)
        return err;

    spreadsdp.spot_lag = spotlag;

    /* Create the floating leg: dates, cvg, spreads  */
    err = create_GenSwapLeg(
        today, &spreadsdp, 0.0, 0.0, 0.0, SPREAD_LEG, SRT_RECEIVER, NULL, ref_rate_code, spreadleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(spreadleg);
        return err;
    }

    /* Sets the payments in the spread leg (no margin): payment = spread*cvg */
    err = populate_GenSwapLeg(0.0, NULL, spreadleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(spreadleg);
        return err;
    }

    /* Return a success message */
    return err;

} /* END swp_f_make_SpreadLeg(...) */

Err swp_f_make_SpreadAndNotionalsLeg(
    long        start,
    long        end,
    double      initial,
    double      final,
    Date        today,
    String      ref_rate_code,
    GenSwapLeg* floatleg)
{
    Err            err = NULL;
    SwapDP         floatsdp;
    int            spotlag;
    SrtCompounding comp;
    BasisCode      basis;

    /* Get the details of the floating reference rate: compounding and basis */
    err = swp_f_get_ref_rate_details(ref_rate_code, &basis, &comp);
    if (err)
        return err;

    /* Set up the corresponding SwapDP: start, end,... generated BACKWARD for theo end date */
    err                = swp_f_setSwapDP(start, end, comp, basis, &floatsdp);
    floatsdp.direction = BKWD;
    if (err)
        return err;

    /* Get the spot lag from the ref rate and attach it to the SwapDp */
    err = srt_f_get_spot_lag_from_refrate(ref_rate_code, &spotlag);
    if (err)
        return err;

    floatsdp.spot_lag = spotlag;

    /* Create the floating leg: dates, cvg, spreads  */
    err = create_GenSwapLeg(
        today,
        &floatsdp,
        0.0,
        initial,
        final,
        SPREAD_AND_NOTIONALS_LEG,
        SRT_RECEIVER,
        NULL,
        ref_rate_code,
        floatleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(floatleg);
        return err;
    }

    /* Sets the payments in the spread leg (no margin): payment = spread*cvg */
    err = populate_GenSwapLeg(0.0, NULL, floatleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(floatleg);
        return err;
    }

    /* Return a success message */
    return err;

} /* END swp_f_make_SpreadAndNotionalsLeg(...) */

/* ------------------------------------------------------------------------- */

/* A useful function that creates (and returns) the pay dates and coverages of
   a fixed leg specified by a SwapDP (memory allocated inside) */
Err swp_f_make_FixedLegDatesAndCoverages(
    SwapDP*  swapdp,
    Date     today,
    Date**   pay_dates,
    int*     num_pay_dates,
    Date**   start_dates,
    Date**   end_dates,
    double** coverages,
    int*     num_dates)
{
    GenSwapLeg fixedleg;
    Err        err = NULL;

    /* Create a fixed leg with dates and cvg  and a coupon of 1.0 */
    err = create_GenSwapLeg(
        today, swapdp, 1.0, 0.0, 0.0, FIXED_LEG, SRT_RECEIVER, NULL, NULL, &fixedleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        return err;
    }

    /* Transfer required information */
    *pay_dates     = (Date*)fixedleg.pay_date.date;
    *num_pay_dates = fixedleg.pay_date.len;
    *start_dates   = fixedleg.start_date.date;
    *end_dates     = fixedleg.end_date.date;
    *coverages     = fixedleg.cvg.d;
    *num_dates     = fixedleg.end_date.len;

    /* Free what is not needed */
    free_inDlist(&(fixedleg.time));
    free_inDlist(&(fixedleg.cpn));
    swp_f_free_in_DateList(fixedleg.fixing_date);

    /* Return a success message */
    return NULL;

} /* END Err swp_f_make_FixedLegDatesAndCoverages(...) */

/* ------------------------------------------------------------------------- */

/* A useful function that creates (and returns) the pay dates and coverages
  of a floating leg specified by a SwapDP (memory allocated inside) */
Err swp_f_make_FloatLegDatesAndCoverages(
    SwapDP*  swapdp,
    Date     today,
    Date**   pay_dates,
    int*     num_pay_dates,
    Date**   fixing_dates,
    Date**   start_dates,
    Date**   end_dates,
    double** coverages,
    int*     num_dates)
{
    GenSwapLeg floatleg;
    Err        err = NULL;

    /* Make a CASH floating leg (without computing the FRA's): dates, times, cvg */
    err = create_GenSwapLeg(
        today, swapdp, 1.0, 0.0, 0.0, SPREAD_LEG, SRT_RECEIVER, NULL, "CASH", &floatleg);
    if (err)
    {
        swp_f_freein_GenSwapLeg(&floatleg);
        return err;
    }

    /* Transfer required information */
    *pay_dates     = (Date*)floatleg.pay_date.date;
    *num_pay_dates = floatleg.pay_date.len;
    *fixing_dates  = (Date*)floatleg.fixing_date.date;
    *start_dates   = (Date*)floatleg.start_date.date;
    *end_dates     = (Date*)floatleg.end_date.date;
    *coverages     = floatleg.cvg.d;
    *num_dates     = floatleg.start_date.len;

    /* Free what is not needed */
    free_inDlist(&(floatleg.time));
    free_inDlist(&(floatleg.payment));
    free_inDlist(&(floatleg.spread));

    /* Return a success message */
    return NULL;

} /* END swp_f_make_FloatLegDatesAndCoverages(...) */

/* ---------------------------------------------------------------------------- */

/* A useful function that creates (and returns) the pay dates, coverages and
   spreads of a floating leg specified by a SwapDP (memory allocated inside) */
Err swp_f_make_FloatLegDatesCoveragesAndSpreads(
    SwapDP*  swapdp,
    Date     today,
    String   refrate,
    Date**   pay_dates,
    int*     num_pay_dates,
    Date**   fixing_dates,
    Date**   start_dates,
    Date**   end_dates,
    double** coverages,
    double** spreads,
    int*     num_dates)
{
    GenSwapLeg floatleg;
    Err        err = NULL;

    /* Make the floating leg: dates, times, cvg, spreads ...*/
    err = swp_f_make_SpreadLeg(swapdp->start, swapdp->end, today, refrate, &floatleg);
    if (err)
        return err;

    /* Transfer required information */
    *pay_dates     = (Date*)floatleg.pay_date.date;
    *num_pay_dates = floatleg.pay_date.len;
    *fixing_dates  = (Date*)floatleg.fixing_date.date;
    *start_dates   = (Date*)floatleg.start_date.date;
    *end_dates     = (Date*)floatleg.end_date.date;
    *coverages     = floatleg.cvg.d;
    *spreads       = floatleg.spread.d;
    *num_dates     = floatleg.start_date.len;

    /* Free what is not needed */
    free_inDlist(&(floatleg.time));
    free_inDlist(&(floatleg.payment));

    /* Return a success message */
    return NULL;

} /* END swp_f_make_FloatLegDatesCoveragesAndSpreads(...) */

/* ---------------------------------------------------------------------------- */

/* A weird function that can be needed to merge on one leg two different ones:
   dates are merged, and payment are computed as : + for recleg, - for payleg  */
Err swp_f_merge_SwapLegs(GenSwapLeg* recleg, GenSwapLeg* payleg, GenSwapLeg* bigleg)
{
    Err err = NULL;
    int i, j, k, i_over, j_over;

    /* Set everything in the new Leg to 0 */
    memset(bigleg, 0, sizeof(GenSwapLeg));

    /* Check if at least start date and end date match (otherwise no sense) */
    if (recleg->sdp.start != payleg->sdp.start)
        return serror("Start dates do not match in swp_f_merge_SwapLegs");
    if (recleg->sdp.end != payleg->sdp.end)
        return serror("End dates do not match in swp_f_merge_SwapLegs");
    if (recleg->today != payleg->today)
        return serror("Today's dates do not match in swp_f_merge_SwapLegs");

    /* Sets today in the leg */
    bigleg->today = recleg->today;

    /* See how many effective date they are in both legs */
    bigleg->leg_length = 0;
    i                  = 0;
    j                  = 0;
    i_over             = 0;
    j_over             = 0;
    while ((i_over == 0) || (j_over == 0))
    {
        if (recleg->pay_date.date[i] < payleg->pay_date.date[j])
            i++;
        else if (recleg->pay_date.date[i] > payleg->pay_date.date[j])
            j++;
        else if (recleg->pay_date.date[i] == payleg->pay_date.date[j])
        {
            j++;
            i++;
        }

        /* Make sure we do not go beyond the last date of the leg : if we do, set the over flag to 1
         */
        if (i == recleg->leg_length)
        {
            i      = recleg->leg_length - 1;
            i_over = 1;
        }
        if (j == payleg->leg_length)
        {
            j      = payleg->leg_length - 1;
            j_over = 1;
        }

        /* This date will have to be added in the big list */
        bigleg->leg_length++;
    } /* END while (...) loop on all leg dates */

    /* Allocate space for the list of dates of the big leg */
    bigleg->pay_date = new_DateList(bigleg->leg_length);

    /* Allocate space for the list of dates of the big leg */
    bigleg->fixing_date = new_DateList(bigleg->leg_length - 1);

    /* Allocate space for bigleg->payment */
    err = new_Dlist(bigleg->leg_length, &(bigleg->payment));
    if (err)
        return (err);

    /* Duplicate the dates and sets the payment information */
    i      = 0;
    j      = 0;
    k      = 0;
    i_over = 0;
    j_over = 0;
    while ((i_over == 0) || (j_over == 0))
    {
        if (recleg->pay_date.date[i] < payleg->pay_date.date[j])
        {
            bigleg->pay_date.date[k] = recleg->pay_date.date[i];
            if (i < recleg->leg_length - 1)
                bigleg->fixing_date.date[k] = recleg->fixing_date.date[i];
            bigleg->payment.d[k] = recleg->payment.d[i];
            i++;
        }
        else if (recleg->pay_date.date[i] > payleg->pay_date.date[j])
        {
            bigleg->pay_date.date[k] = payleg->pay_date.date[j];
            if (j < payleg->leg_length - 1)
                bigleg->fixing_date.date[k] = payleg->fixing_date.date[j];
            bigleg->payment.d[k] = -payleg->payment.d[j];
            j++;
        }
        else if (recleg->pay_date.date[i] == payleg->pay_date.date[j])
        {
            bigleg->pay_date.date[k] = recleg->pay_date.date[i];
            if (i < recleg->leg_length - 1)
                bigleg->fixing_date.date[k] = recleg->fixing_date.date[i];
            bigleg->payment.d[k] = recleg->payment.d[i] - payleg->payment.d[j];
            j++;
            i++;
        }

        /* Make sure we do not go beyond the last date of the leg: if we do, set the over flag to 1
         */
        if (i == recleg->leg_length)
        {
            i      = recleg->leg_length - 1;
            i_over = 1;
        }
        if (j == payleg->leg_length)
        {
            j      = payleg->leg_length - 1;
            j_over = 1;
        }

        /* Move on to the next date */
        k++;
    } /* END while (...) loop on all leg dates */

    /* Return a success message */
    return NULL;

} /* END Err swp_f_merge_SwapLegs(...) */

/* ====================================================================== */
