/* ===========================================================================

         FILENAME:     swp_h_swap_generic.h

     PURPOSE:      A few useful functions to deal with a generic swap structure

   =========================================================================== */
#ifndef SWP_H_SWAP_GENERIC_H
#define SWP_H_SWAP_GENERIC_H

#include "swp_h_datelist.h"
#include "swp_h_swap_structures.h"
#include "swp_h_swap_utils.h"

#define DEFAULT_NOTIONAL 100000000.0

/* ------------------------------------------------------------------------ */
/* Creates the swap and sets everything that is stripper independent */

Err swp_f_create_GenFullSwap(
    SwapDP*      main_leg_sdp,
    SwapDP*      oth_leg_sdp,
    double       fixed,
    double       initial_exchange,
    double       final_exchange,
    SwapType     swap_type,
    String       curve_name,
    GenFullSwap* swap);

/* ------------------------------------------------------------------------ */
/* Frees everything in a swap  or the leg of a swap */
Err swp_f_freein_GenSwapLeg(GenSwapLeg* leg);

Err swp_f_freein_GenFullSwap(GenFullSwap* swap);

/* ------------------------------------------------------------------------ */
/* According to the swap type, sets what is loop/stripper/mkt dependent */

Err swp_f_populate_GenFullSwap(double margin, SrtCurvePtr crv, GenFullSwap* swap);

/* ------------------------------------------------------------------------- */
/* Sets discount factors and values the swap, or the swap leg                */
Err swp_f_value_GenFullSwap(GenFullSwap* swap, SrtCurvePtr crv);

/* ------------------------------------------------------------------------- */
/* THESE FUNCTIONS ARE VERY USEFUL, AND EASIER TO USE IN SIMPLE CASES        */
/* ------------------------------------------------------------------------- */

/* Make the dates and coupons for a fixed leg (+ sets payments) */
Err swp_f_make_FixedLeg(SwapDP* fixedsdp, double fixed, Date today, GenSwapLeg* fixedleg);

/* Make the dates and coupons for a fixed leg + notional (cash flow at spot) */
Err swp_f_make_FixedAndNotionalsLeg(
    SwapDP* fixedsdp, double fixed, double initial, double final, Date today, GenSwapLeg* fixedleg);

Err swp_f_make_NotionalsLeg(
    SwapDP* fixedsdp, double initial, double final, Date today, GenSwapLeg* fixedleg);

/* Make the dates and spreads + initial + final for a floating leg + payments */
Err swp_f_make_FloatingLeg(
    long start, long end, Date today, String ref_rate_code, SrtCurvePtr crv, GenSwapLeg* floatleg);

/* Make the dates and spreads + initial + final  for a floating leg (no fra's) */
Err swp_f_make_SpreadLeg(
    long start, long end, Date today, String ref_rate_code, GenSwapLeg* spreadleg);

Err swp_f_make_SpreadAndNotionalsLeg(
    long        start,
    long        end,
    double      initial,
    double      final,
    Date        today,
    String      ref_rate_code,
    GenSwapLeg* floatleg);

/* ------------------------------------------------------------------------- */

/* A useful function that creates (and returns) the pay dates and coverages of
   a leg specified by a SwapDP (memory allocated inside) */
Err swp_f_make_FixedLegDatesAndCoverages(
    SwapDP*  swapdp,
    Date     today,
    Date**   pay_dates,
    int*     num_pay_dates,
    Date**   start_dates,
    Date**   end_dates,
    double** coverages,
    int*     num_dates);

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
    int*     num_dates);

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
    int*     num_dates);

/* ------------------------------------------------------------------------- */

/* A weird function that can be needed to merge on one leg two different ones:
dates are merged, and payment are computed as : + for recleg, - for payleg  */

Err swp_f_merge_SwapLegs(GenSwapLeg* recleg, GenSwapLeg* payleg, GenSwapLeg* bigleg);

#endif
