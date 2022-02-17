/* ===========================================================================

         FILENAME:     swp_h_swap_types.h

     PURPOSE:      All the possible types of swaps and swap legs.

   ===========================================================================
 */
#ifndef SWP_H_SWAP_TYPES_H
#define SWP_H_SWAP_TYPES_H

/* ------------------------------------------------------------------------ */
/* The different types of legs that one can find in a swap                  */
typedef enum LegType {
  NOTIONALS_LEG,
  FIXED_LEG,
  FIXED_AND_NOTIONALS_LEG,
  FLOATING_LEG,
  SPREAD_AND_NOTIONALS_LEG,
  SPREAD_LEG,
  FLOATING_MARGIN_LEG,
  DRS_MARGIN_LEG,
  CMS_MARGIN_LEG,
  CMT_MARGIN_LEG,
  CMS_FOR_TEC_MARGIN_LEG,
  TEC_MARGIN_LEG,
  DRS_LEG,
  STANDARD_BOND_LEG
} LegType;

/* ------------------------------------------------------------------------ */
/* The different types of swaps that can be created                         */
typedef enum SwapType {

  /* Full swaps: two different legs */
  /* Fixed leg vs exchange of notionals (initial & final) leg*/
  FIXED_NOTIONALS_SWAP,
  /* Fixed leg vs floating (with cash spreads) leg */
  FIXED_FLOATING_SWAP,

  /* CMS minus margin leg vs Libor leg(with cash spreads) */
  CMS_MARGIN_FLOATING_SWAP,
  /* CMS minus margin leg vs fixed leg */
  CMS_MARGIN_FIXED_SWAP,
  /* CMT minus margin leg vs Libor leg (with cash spreads) */
  CMT_MARGIN_FLOATING_SWAP,
  /* CMT minus margin leg vs fixed leg */
  CMT_MARGIN_FIXED_SWAP,
  /* CMS minus margin (the TEC way) leg vs fixed leg */
  CMS_FOR_TEC_MARGIN_FLOATING_SWAP,
  /* CMS minus margin (the TEC way) leg vs fixed leg */
  CMS_FOR_TEC_MARGIN_FIXED_SWAP,
  /* TEC minus margin (the TEC way) leg vs Libor leg (with cash spreads) */
  TEC_MARGIN_FLOATING_SWAP,
  /* TEC minus margin (the TEC way) leg vs fixed leg */
  TEC_MARGIN_FIXED_SWAP,

  /* DRS leg minus margin vs floating (with cash spreads) */
  DRS_MARGIN_FLOATING_SWAP,
  /* DRS leg minus margin vs fixed leg */
  DRS_MARGIN_FIXED_SWAP,

  /* Incomplete swaps: only one leg */
  /* Fixed leg and (both) notionals in one leg */
  ONE_LEG_FIXED_AND_NOTIONALS,
  /* Fixed leg and final notional in one leg */
  ONE_LEG_STANDARD_BOND,
  /* Fixed leg only in one leg (for level payment) */
  ONE_LEG_FIXED_SWAP,

  /* CMT      , CMS legs... */
  ONE_LEG_CMS_MARGIN,
  ONE_LEG_CMT_MARGIN,
  ONE_LEG_CMS_FOR_TEC_MARGIN,
  ONE_LEG_TEC_MARGIN,

  ONE_LEG_DRS
} SwapType;

/* ------------------------------------------------------------------------ */

#endif
