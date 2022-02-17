/*  SWP_H_BOND_COMPUTE.H   version 1.0      SWAPS and OPTIONS RESEARCH TEAM

Declaration of all functions and constants used by srt_f_bond_compute.c

Last version : June 23th, 1994.             
Author : Matthias Lennkh.
Mods   : Julia Matsumoto, Jasbir Malhi, Alexandre Benech.            

Last modified: Dec 16, 1994, K L Chau, use new SrtMktPtr structure
*/
#ifndef SWP_H_BOND_COMPUTE_H
#define SWP_H_BOND_COMPUTE_H

#define 	EPSILON 	(double)0.00000001
double dirty_price_fct(SwapDP p,double coupon, double irr, 
			double redemption, double first_coupon);   
double acc_int_fct(SwapDP p, double coupon);
double clean_price_fct(SwapDP p,double coupon, double irr, 
			double redemption, double first_coupon);
double new_set_price_fct(Date acc_date,Date new_date,SwapDP p,
			double coupon, double clean_price,
			String swap, String repo, 
			double redemption, double first_coupon);  
double yield_fct(SwapDP p,double coupon, double clean_price, 
			double redemption, double first_coupon);
double bond_sens_fct(SwapDP p,double coupon, double irr, 
			double redemption, double first_coupon);
double fwd_clean_price_fct( Date fut, SwapDP p, double coupon, 
			double clean_price, String swap, 
			String repo, double first_coupon);
double tdy_price_fct( Date fut, SwapDP p, double coupon, double fwd_price,
                           String swap, String repo, double first_coupon);
double fwd_irr_fct( Date fut, SwapDP p, double coupon, double clean_price,
                    String swap, String repo, double redemption,
			double first_coupon); 
double fwd_irr_fct_smp( Date fut, SwapDP p, double coupon, double clean_price, 
			double volatility, String swap, String repo, 
			double redemption, double first_coupon); 


#endif