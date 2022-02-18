/**
        FILE_NAME:	swp_h_gen_deal.h
        Include file for swaps.c

        AUTHOR	:	G. Amblard

        First Version:	September 10, 1993
        Last Modified:	September 11, 1993

**/
#ifndef SWP_H_GEN_DEAL_H
#define SWP_H_GEN_DEAL_H

#include "swp_h_swap_structures.h"

/** Functions called from other files **/

double value_swap(SrtCurvePtr crv, Swap_Obj* swap);

Swap_Obj* generate_swap(Arg_Obj* arg, SwapType mess);

/** Support functions which could potentially be called directly **/

Leg_Obj* generate_leg(Arg_Obj* arg, LegType mess, int rec_pay);

void leg_eval_payment(Leg_Obj* leg, String curvename);
void leg_eval_disc_factor(Leg_Obj* leg, String curvename);

/* For CMT stripper */
Err leg_cms_populate(Leg_Obj* leg, SrtCurvePtr crv);
Err leg_cmt_populate(Leg_Obj* leg, SrtCurvePtr crv);

#endif
