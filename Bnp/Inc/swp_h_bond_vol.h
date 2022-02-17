/* ------------------------------------------------------------------------------
   FILENAME:		swp_h_bond_vol.h

   PURPOSE:		    Transform a beta bond option vol into a BS vol
   ------------------------------------------------------------------------------
 */
#ifndef SWP_H_BOND_VOL_H
#define SWP_H_BOND_VOL_H

Err swp_f_bond_vol(Date fut, SwapDP p, double coupon, double clean_price,
                   double strike, double volatility, Date today, double alpha,
                   double beta, double rho, String swap, String repo,
                   double redemption, double first_coupon, double *implied_vol);

Err swp_f_bond_vol_yldbeta_to_pricebs(double fwd_clean_price,
                                      double clean_strike, double yield,
                                      double yieldmean, double yieldstar,
                                      double PsiY, double PsiYY, double PsiYYY,
                                      double volatility, double optmat,
                                      double alpha, double beta, double rho,
                                      double *implied_vol);

Err swp_f_bond_vol_yldbeta_to_yldbs(double fwd_clean_price, double clean_strike,
                                    double yield, double yieldmean,
                                    double yieldstar, double PsiY, double PsiYY,
                                    double PsiYYY, double volatility,
                                    double optmat, double alpha, double beta,
                                    double rho, double *implied_vol);

#endif