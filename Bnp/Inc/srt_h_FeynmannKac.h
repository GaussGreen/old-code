

Err srt_f_feynmann_kac_pde(long			    num_mesh,
						 long			    num_time_step,
						 double             dX,
						 double				time_step,
						 double				*X,
						 double				**conv_mat, /* CONVECTION MATRIX */
						 double				**diff_mat, /* DIFFUSION MATRIX */
						 SrtPdePDEType		pde_type,
						 SrtPdeBoundaryCond	pde_bound_cond,
						 double				*pde_bound,
						 double				*term_payoff, /* FINAL MATURITY PAYOFF */
						 double				**amer_payoff, /* AMERICAN PAYOFF */

						 double				**div_mat,
					  
						 double				**U);

Err srt_f_black_scholes_pde(Date		eval_date,
							Date		exp_date,
							double		spot,
							double      forward,
							double      equity_strike,
							double      equity_vol,
							double		min_node,
							double		min_mesh,
							
							/*OUTPUT */
							double		*bs_pv);
