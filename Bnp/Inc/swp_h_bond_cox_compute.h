/* BOND_OPTIONS.H version 1.00			SORT	*/
/* declariation of all functions used by BOND_OPTIONS.C */
/* Last version June 17th				*/
/* Author Alexandre BENECH				*/
/* Last modified: Dec 2  , 1994 */
/* Author: K L Chau */
/* Reason: use new SrtMkt structure */

#ifndef SRT_H_BOND_COX_COMPUTE_H
#define SRT_H_BOND_COX_COMPUTE_H

double cox_clean_fct(SwapDP p, double coupon, double spot, Date value_date,
                     Date opt_mat, double vol, double strike, double step,
                     SrtCallPutType call_put, SrtGreekType greek, String swap,
                     String repo, double redemption, double first_coupon);

double cox_fwd_yield_fct(SwapDP p, double coupon, double spot, Date value_date,
                         Date opt_mat, double vol, double strike, double step,
                         SrtCallPutType call_put, SrtGreekType greek,
                         String swap, String repo, double redemption,
                         double first_coupon);

double cox_yield_fct(SwapDP p, double coupon, double spot, Date value_date,
                     Date opt_mat, double vol, double strike, double step,
                     SrtCallPutType call_put, SrtGreekType greek, String swap,
                     String repo, double redemption, double first_coupon);

double am_cox_clean_fct(SwapDP p, double coupon, double spot, Date value_date,
                        Date opt_mat, double vol, double strike, double step,
                        SrtCallPutType call_put, SrtGreekType greek,
                        String swap, String repo, double redemption,
                        double first_coupon);

double am_cox_yield_fct(SwapDP p, double coupon, double spot, Date value_date,
                        Date opt_mat, double vol, double strike, double step,
                        SrtCallPutType call_put, SrtGreekType greek,
                        String swap, String repo, double redemption,
                        double first_coupon);

#endif