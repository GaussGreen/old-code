/* ===========================================================================

         FILENAME:     swp_h_swap_structures.h

     PURPOSE:      A few useful structures to represent a generic swap

   =========================================================================== */
#ifndef SWP_H_SWAP_STRUCTURES_H
#define SWP_H_SWAP_STRUCTURES_H

#include "swp_h_datelist.h"
#include "swp_h_swap_types.h"
#include "swp_h_swap_utils.h"

/* ------------------------------------------------------------------------ */
/* THE structure to represent a swap leg: type, dates, df, cvg, ...         */
typedef struct gen_swap_leg
{
    LegType         leg_type;
    CcyCode         ccy;
    SrtReceiverType rec_pay;
    SwapDP          sdp;

    DateList fixing_date;
    DateList start_date;
    DateList end_date;
    DateList pay_date; /* Includes notional payment dates */

    int leg_length; /* Num of payment dates */

    Dlist  time; /* Time from today to start date */
    Dlist  mat;  /* Time from today to fixing date (spot lag before) */
    Dlist  temp_fwd;
    Dlist  fwd;
    Dlist  vol;
    Dlist  cpn; /* == convexity adjusted forward or fixed or fra...*/
    Dlist  spread;
    Dlist  cvg;
    Dlist  payment;
    Dlist  df;
    double margin;
    double initial_exchange;
    double final_exchange;
    double leg_value;
    long   today;
} GenSwapLeg;

/* ------------------------------------------------------------------------ */
/* THE structure to represent a full swap, made of two legs                 */
typedef struct gen_full_swap
{
    SwapType   swap_type;
    GenSwapLeg leg[2];
    long       today;
    long       spotdate;
    long       enddate;
    double     notional;
    double     swap_value;
} GenFullSwap;

/* ------------------------------------------------------------------------ */

#endif
