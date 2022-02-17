/* =============================================================================

   FILENAME:    swp_h_swap_pricing.h

  PURPOSE:     Provide a few function for generic swap calculations:
                                - PV
                                - Fwd Rate
                                - Level
                                - Margin...
   =============================================================================
*/

#ifndef SWP_H_SWAP_PRICING_H
#define SWP_H_SWAP_PRICING_H

Err swp_f_GenericSwap(long start, long end_nfp, String compStr, String basisStr,
                      double coupon, double initial, double final,
                      String ycname, String strMessage, long value_date,
                      double *answer);

/* -----------------------------  FORWARD (SWAP)
 * -------------------------------- */

/* Compute a forward swap rate (for use outside SORT) */
Err swp_f_ForwardRate(long start, long end_nfp, String compStr, String basisStr,
                      String ycName, String refRateCode, double *forward);

/* Compute a forward swap rate (for use inside SORT) */
Err swp_f_ForwardRate_SwapDP(SwapDP *swapdp, String ycName, String refRateCode,
                             double *forward);

/* ----------------------------------- LEVEL
 * --------------------------------------- */

/* Compute a level payment (for use outside SORT) */
Err swp_f_LevelPayment(long start, long end_nfp, String compStr,
                       String basisStr, String ycName, String refRateCode,
                       double *forward);

/* Compute a level payment (for use inside SORT) */
Err swp_f_Level_SwapDP(SwapDP *swapdp, String ycName, double *level);

/* -----------------------------  MARGIN (ON ONE LEG)
 * -------------------------------- */
/* Computes a margin on the leg of a swap to reach a given PV (for use outside
 * SORT) */
Err swp_f_SwapMargin(double pv, long start, long end_nfp, String compStr,
                     String basisStr, String ycName, double *margin);

/* Computes a margin on the leg of a swap to reach a given PV with a SwapDP */
Err swp_f_SwapMargin_SwapDP(double pv, SwapDP *swapdp, String ycName,
                            double *margin);

/* -----------------------------  PV (SWAP) -------------------------------- */

/* Computes the PV of a Swap  , taking spreads into account (for use outside
 * SORT) */
Err swp_f_SwapPv(double coupon, long start, long end_nfp, String compStr,
                 String basisStr, String ycName, String refRateCode,
                 double *pv);

/* Computes the PV of a Swap  , taking spreads into account */
Err swp_f_SwapPv_SwapDP(double coupon, SwapDP *swapdp, String ycName,
                        String refRateCode, double *pv);

#endif
