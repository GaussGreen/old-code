/* =============================================================================

   FILENAME:    swp_f_swap_pricing.c

   PURPOSE:     Provide a few function for generic swap calculations:
                                - PV
                                - Fwd Rate
                                - Level
                                - Margin...
   =============================================================================
 */

#include "swp_h_all.h"
#include "swp_h_swap_pricing.h"

/* ---------------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------------
   This function is intended as an interface function
   ---------------------------------------------------------------------------------
 */
Err swp_f_GenericSwap(long start, long end_nfp, String compStr, String basisStr,
                      double coupon, double initial, double final,
                      String ycname, String strMessage, long value_date,
                      double *answer) {
  Err err = NULL;

  return err;
}

/* ------------------------------------------------------------------------------
 */

/* -----------------------------  FORWARD (SWAP)
 * -------------------------------- */

/* Computes a forward Swap rate  , taking spreads into account */
Err swp_f_ForwardRate(long start, long end_nfp, String compStr, String basisStr,
                      String ycName, String refRateCode, double *forward) {
  SwapDP swapdp;
  Err err = NULL;

  /* Sets the dates parameters for the swap */
  err = swp_f_initSwapDP(start, end_nfp, compStr, basisStr, &swapdp);
  if (err)
    return err;

  /* Compute the forward rate */
  err = swp_f_ForwardRate_SwapDP(&swapdp, ycName, refRateCode, forward);
  if (err)
    return err;

  /* Return a success string */
  return NULL;
}

/* ------------------------------------------------------------------------------
 */
/* Computes a forward Swap rate  , taking spreads into account with a SwapDP */

Err swp_f_ForwardRate_SwapDP(SwapDP *swapdp, String ycName, String refRateCode,
                             double *forward) {
  Err err = NULL;
  GenFullSwap fullswap;
  double float_value;
  double level_value;
  Date today;
  SrtCurvePtr yldcrv;

  /* Get the yield curve */
  yldcrv = lookup_curve(ycName);
  if (!yldcrv)
    return serror("Could not find  %s yc in swp_f_ForwardRate", ycName);

  /* Gets today from the YC */
  today = get_clcndate_from_yldcrv(yldcrv);

  /* Sets all the elements in the Swap to null to prevent memory explosions */
  memset(&fullswap, 0, sizeof(GenFullSwap));

  /* Make the fixed leg : a level payment */
  err = swp_f_make_FixedLeg(swapdp, 1.0, today, &fullswap.leg[0]);
  if (err)
    return err;

  /* Construct the floating leg  , keeping the same end date as the fixed leg */
  err = swp_f_make_FloatingLeg(swapdp->start, swapdp->end, today, refRateCode,
                               yldcrv, &fullswap.leg[1]);
  if (err) {
    swp_f_freein_GenFullSwap(&fullswap);
    return err;
  }

  /* Value the swap: each leg independently will matter */
  err = swp_f_value_GenFullSwap(&fullswap, yldcrv);
  if (err) {
    swp_f_freein_GenFullSwap(&fullswap);
    return err;
  }

  /* The forward rate is such that  rate * level = pv(float) */
  level_value = fullswap.leg[0].leg_value;
  float_value = fullswap.leg[1].leg_value;

  *forward = float_value / level_value;

  /* Free memory */
  swp_f_freein_GenFullSwap(&fullswap);

  /* Return a success string */

  return NULL;
}

/* ---------------------------------------------------------------------------------
 */

/* ----------------------------------- LEVEL
 * --------------------------------------- */

Err swp_f_LevelPayment(long start, long end_nfp, String compStr,
                       String basisStr, String ycName, String refRateCode,
                       double *forward) {
  SwapDP swapdp;
  Err err = NULL;

  /* Sets the dates parameters for the swap */
  err = swp_f_initSwapDP(start, end_nfp, compStr, basisStr, &swapdp);
  if (err)
    return err;

  /* Compute the level payment */
  err = swp_f_Level_SwapDP(&swapdp, ycName, forward);
  if (err)
    return err;

  /* Return a success string */
  return NULL;
}

/* ---------------------------------------------------------------------------------
 */

Err swp_f_Level_SwapDP(SwapDP *swapdp, String ycName, double *level) {
  Err err = NULL;
  GenFullSwap fullswap;
  Date today;
  SrtCurvePtr yccrv;

  /* Get the Yield Curve crv */
  yccrv = lookup_curve(ycName);
  if (!yccrv)
    return serror("Could not find  %s yc in swp_f_ForwardRate", ycName);

  /* Gets today from the YC */
  today = get_clcndate_from_yldcrv(yccrv);

  /* Sets everything to zero in the swap */
  memset(&fullswap, 0, sizeof(GenFullSwap));

  /* Make the fixed leg : a level payment (coupon of 1.0) and sets payments */
  err = swp_f_make_FixedLeg(swapdp, 1.0, today, &fullswap.leg[0]);
  if (err)
    return err;

  /* Value the one legged swap */
  err = swp_f_value_GenFullSwap(&fullswap, yccrv);
  if (err) {
    swp_f_freein_GenFullSwap(&fullswap);
    return err;
  }

  /* The level is the pv of the leg */
  *level = fullswap.leg[0].leg_value;

  /* Free memory */
  swp_f_freein_GenFullSwap(&fullswap);

  /* Return a success string */

  return NULL;

} /* END Err swp_f_Level_SwapDP(...) */

/* ------------------------------------------------------------------------------------
 */

/* -----------------------------  MARGIN (ON ONE LEG)
 * -------------------------------- */

/* Computes a margin on the leg of a swap to reach a given PV */
Err swp_f_SwapMargin(double pv, long start, long end_nfp, String compStr,
                     String basisStr, String ycName, double *margin) {
  SwapDP swapdp;
  Err err = NULL;

  /* Sets the dates parameters for the swap */
  err = swp_f_initSwapDP(start, end_nfp, compStr, basisStr, &swapdp);
  if (err)
    return err;

  /* Compute the margin */
  err = swp_f_SwapMargin_SwapDP(pv, &swapdp, ycName, margin);
  if (err)
    return err;

  /* Return a success string */
  return NULL;

} /* END Err swp_f_SwapMargin(...) */

/* ------------------------------------------------------------------------------
 */

/* Computes a margin on the leg of a swap to reach a given PV with a SwapDP */

Err swp_f_SwapMargin_SwapDP(double pv, SwapDP *swapdp, String ycName,
                            double *margin) {
  Err err = NULL;
  double level;
  SrtCurvePtr yccrv;

  /* Get the Yield Curve crv */
  yccrv = lookup_curve(ycName);
  if (!yccrv)
    return serror("Could not find  %s yc in swp_f_ForwardRate", ycName);

  /* Compute the level payment */
  err = swp_f_Level_SwapDP(swapdp, ycName, &level);
  if (err)
    return err;

  /* The Margin is the PV divided by the level */
  *margin = pv / level;

  /* Return a success string */

  return NULL;
}

/* ---------------------------------------------------------------------------------
 */

/* -----------------------------  PV (SWAP) -------------------------------- */

/* Computes the PV of a Swap  , taking spreads into account */
Err swp_f_SwapPv(double coupon, long start, long end_nfp, String compStr,
                 String basisStr, String ycName, String refRateCode,
                 double *pv) {
  SwapDP swapdp;
  Err err = NULL;

  /* Sets the dates parameters for the swap */
  err = swp_f_initSwapDP(start, end_nfp, compStr, basisStr, &swapdp);
  if (err)
    return err;

  /* Compute the PV */
  err = swp_f_SwapPv_SwapDP(coupon, &swapdp, ycName, refRateCode, pv);
  if (err)
    return err;

  /* Return a success string */
  return NULL;
}

/* ------------------------------------------------------------------------------
 */
/* Computes the PV of a Swap  , taking spreads into account with a SwapDP */

Err swp_f_SwapPv_SwapDP(double coupon, SwapDP *swapdp, String ycName,
                        String refRateCode, double *pv) {
  Err err = NULL;
  GenFullSwap fullswap;
  double float_value;
  double fixed_value;
  double swap_value;
  Date today;
  SrtCurvePtr yccrv;

  /* Get the Yield Curve crv */
  yccrv = lookup_curve(ycName);
  if (!yccrv)
    return serror("Could not find  %s yc in swp_f_ForwardRate", ycName);

  /* Gets today from the YC */
  today = get_clcndate_from_yldcrv(yccrv);

  /* Sets all the elements in the Swap to null to prevent memory explosions */
  memset(&fullswap, 0, sizeof(GenFullSwap));

  /* Make the fixed leg with coupon */
  err = swp_f_make_FixedLeg(swapdp, coupon, today, &fullswap.leg[0]);
  if (err)
    return err;

  /* Construct the floating leg  , keeping the same end date as the fixed leg */
  err = swp_f_make_FloatingLeg(swapdp->start, swapdp->end, today, refRateCode,
                               yccrv, &fullswap.leg[1]);
  if (err) {
    swp_f_freein_GenFullSwap(&fullswap);
    return err;
  }

  /* Value the swap (both legs) */
  err = swp_f_value_GenFullSwap(&fullswap, yccrv);
  if (err) {
    swp_f_freein_GenFullSwap(&fullswap);
    return err;
  }

  /* Swap PV is the fixed leg PV minus the floating leg PV */
  fixed_value = fullswap.leg[0].leg_value;
  float_value = fullswap.leg[1].leg_value;
  swap_value = fixed_value - float_value;

  /* Free memory */
  swp_f_freein_GenFullSwap(&fullswap);

  /* Return a success string */
  return NULL;
}

/* ---------------------------------------------------------------------------------
 */
