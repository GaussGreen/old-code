
	
Err srt_f_cpi_forward(char		*cpi_und_name,
					  Date		fixing_date,
					  Date		pay_date,
					  double	*fwd_cpi);

Err srt_f_vasicek_calibrated_mean_rev_lvl_and_initial_cond(Date		   spot_date,
														   Date			ipc_date_lm,
														   Date			ipc_date_hm,
														   double		interp_coeff,
														   double		daily_ref_index_target,

														   long			num_bonds,
														   long			*num_cpns,
														   double		**pay_dates,
														   double		**cpns,
														   double		*bond_dirty_prices,
														   SRT_Boolean	discount_fwd_dirty_pr, 

														   SRT_Boolean	price_floor,
														   double		*floor_gearing,
														   double		*pre_factor_fwd,
														   double		*floor_strike,
														   double		*floor_implied_vol,

														   char			*vasicek_und_name,
														   char			*price_und_name,
														   char			*nom_swap_market,
														   SRT_Boolean	price_floor_boolean,

														   double		***calibrated_mean_rev_lvl,
														   double		*calibrated_vasicek_init_cond);