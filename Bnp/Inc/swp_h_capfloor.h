/* ===================================================================================
   FILENAME:      swp_h_capfloor.h

   PURPOSE:       Compute cap/floors prices  , and implied volatilities
   ===================================================================================
 */

#ifndef SWP_H_CAPFLOOR_H
#define SWP_H_CAPFLOOR_H

/* -----------------------------------------------------------------------------
 */

/* -------------------  CAP/FLOOR PRICE AND DERIVATIVES   ----------------------
 */

/* -----------------------------------------------------------------------------
 */
/* Compute the price or Greeks of a cap/floor on a Reference Rate (use outside
 * SORT)*/
Err swp_f_CapFloor(long start, long end_nfp, double strike,
                   Err (*GetBSVol)(double start, double end, double strike,
                                   double dForward, double dSpread,
                                   double *vol),
                   String capFloorStr, String refRateCode, String ycName,
                   String greekStr, String logNormStr, double *result);

/* -----------------------------------------------------------------------------
 */
/* Compute the price or Greeks of a cap/floor on a Reference Rate (use inside
 * SORT)*/
Err swp_f_CapFloor_SwapDP(SwapDP *swapdp, double strike,
                          Err (*GetBSVol)(double start, double end,
                                          double strike, double dForward,
                                          double dSpread, double *vol),
                          SrtReceiverType capFloor, String refRateCode,
                          String ycName, SrtGreekType greek,
                          SrtDiffusionType logNorm, double *result);

/* -----------------------------------------------------------------------------
 */
/* Call to the relevant BlackScholes formula for the caplet pricing */
Err swp_f_Caplet_Struct(double forward, double strike, double vol, double mat,
                        SrtReceiverType capFloor, SrtGreekType greek,
                        SrtDiffusionType logNorm, double *premium);

/* -----------------------------------------------------------------------------
 */

/* -----------------  IMPLIED VOLATILITY FOR CAP/FLOOR   -----------------------
 */

/* -----------------------------------------------------------------------------
 */
/* Computes the implied vol of a cap/floor using Newton (for use outside SORT)
 */
Err swp_f_CapFloorImpliedVol(double premium, long start, long end_nfp,
                             double strike, String capFloorStr,
                             String refRateCode, String ycName,
                             String logNormStr, double *impvol);

Err swp_f_CapFloorImpVol_SwapDP(double premium, SwapDP *swapdp, double strike,
                                SrtReceiverType capFloor, String refRateCode,
                                String ycName, SrtDiffusionType logNorm,
                                double *impvol);

#endif