/* ------------------------------------------------------------------------------
   FILENAME:      swp_h_bond_swaption.h

   PURPOSE:       Given a bond option, finds the swaption that best matches it
   ------------------------------------------------------------------------------ */
#ifndef SWP_H_BOND_SWAPTION_H
#define SWP_H_BOND_SWAPTION_H

Err swp_f_match_bond_swaption(
    long   option_expiry,
    long   settle_bond_date,
    long   maturity_date,
    String strBondComp,
    String strBondBasis,
    double coupon,          /* as 0.05 for 5% */
    double fwd_clean_price, /* as 1.00 for 100 */
    double clean_strike,
    String strCallPut,
    double bond_option_value,
    String yc_name,
    String refRateCode,
    String strSwapComp,
    String strSwapBasis,
    String strLogorNormal,

    long*   swap_settle,
    long*   swap_maturity,
    double* swap_coupon,
    double* implied_vol,
    double* vega_ratio,
    double* hedge_ratio);

#endif