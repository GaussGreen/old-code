#ifndef SRT_H_TS_CPI_FCT_H
#define SRT_H_TS_CPI_FCT_H
#include "math.h"
#include "srt_h_all.h"
#define MAX_ITER 10

Err srt_f_cpi_forward(char				*price_index_und_name,
					  Date				fixing_date,
					  Date				pay_date,
					  double			*fwd_cpi)
{
	Err			err = NULL;
	SrtUndPtr	price_index_und,nom_rate_und,real_rate_und;
	SrtMdlType  mdl_type;
	char		*nom_rate_und_name,*real_rate_und_name;
	char		*nom_rate_yc;
	double      df_real_fixingdt,df_real_settdt,df_real,df_nom;
	double		fixing_t,pay_t,adj,nom_H,nom_G,date_today;
	TermStruct	*cpi_ts,*nom_ts,*real_ts;
	Date		settlt_date;

	/* get the price index unerlying */
	price_index_und = lookup_und(price_index_und_name);
	if(!price_index_und) 
	{
		return serror("can not get the price index underlying");
	}

	settlt_date = get_spotdate_from_underlying(price_index_und);

	/* get the real rate underlying and its type */

	real_rate_und_name = get_forname_from_fxund(price_index_und);
	real_rate_und = lookup_und(real_rate_und_name);
	if(!real_rate_und) 
	{
		return serror("can not get the real. rate underlying");
	}

	err = get_underlying_mdltype(real_rate_und,&mdl_type);
	if(err) 
	{
		return err;
	}

	if(mdl_type != VASICEK) return serror("real rate und must be of VASICEK type");

	/* get the nominal rate underlying */
	nom_rate_und_name = get_domname_from_fxund(price_index_und);
	nom_rate_und = lookup_und(nom_rate_und_name);
	if(!nom_rate_und) 
	{
		return serror("can not get the nom. rate underlying");
	}
	/* get the nominal rate yield curve */
	err = get_underlying_discname(nom_rate_und, &nom_rate_yc);
	if (err) 
	{
		return err;
	}

	/* get the forward nom. rate discount factor */
	df_nom = swp_f_df(settlt_date,fixing_date,nom_rate_yc);

	/* get the forward real rate discount factor */
	err = srt_f_vasicek_discount_factor(settlt_date,
										real_rate_und,
										&df_real_settdt);
	if(err) return err;

	err = srt_f_vasicek_discount_factor(fixing_date,
										real_rate_und,
										&df_real_fixingdt);
	if(err) return err;

	df_real = df_real_fixingdt/df_real_settdt;
											
	(*fwd_cpi) = get_spot_from_fxund(price_index_und)*df_real/df_nom;

	/* compute the adjustment */
	date_today = get_today_from_underlying (price_index_und);

	fixing_t = (fixing_date-date_today)*YEARS_IN_DAY;
	pay_t = (pay_date-date_today)*YEARS_IN_DAY;

	err = get_underlying_ts(price_index_und,&cpi_ts);
	if(err) return err;

	err = get_underlying_ts(nom_rate_und,&nom_ts);
	if(err) return err;

	err = get_underlying_ts(real_rate_und,&real_ts);
	if(err) return err;

	G_H_func(fixing_t,nom_ts,&nom_G,&nom_H);

	adj = 0.0;
	adj += (Psi_func(fixing_t,nom_ts)*V_dx_func(fixing_t,cpi_ts)-W_dx_func(fixing_t,cpi_ts));
	adj -= (Psi_func(pay_t,nom_ts)*V_dx_func(fixing_t,cpi_ts)-W_dx_func(fixing_t,cpi_ts));

	adj += O_fd_func(fixing_t,cpi_ts)*Psi_func(fixing_t,real_ts)*(Psi_func(pay_t,nom_ts)-Psi_func(fixing_t,nom_ts)); 
	adj -= P_fd_func(fixing_t,cpi_ts)*(Psi_func(pay_t,nom_ts)-Psi_func(fixing_t,nom_ts));

	adj += 2*Psi_func(fixing_t,nom_ts)*L_func(fixing_t,nom_ts)
		 - Psi_func(fixing_t,nom_ts)*Psi_func(pay_t,nom_ts)*nom_G
		 +(Psi_func(fixing_t,nom_ts) + Psi_func(pay_t,nom_ts))*(Psi_func(fixing_t,nom_ts)*nom_G-L_func(fixing_t,nom_ts))
		 - nom_G*Psi_func(fixing_t,nom_ts)*Psi_func(fixing_t,nom_ts);
		   		   
	(*fwd_cpi) *= exp(adj);

	return NULL;

}

Err srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func(double		settlt_date,
														   double		vasicek_init_cond,

														   Date			ipc_date_hm,
														   Date			ipc_date_lm,
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

														   char			*real_rate_und_name,
														   char			*price_und_name,

														   double		***calibrated_mean_rev_lvl,
														   double		*model_daily_ref_index)
{
	Err			err = NULL;
	SrtUndPtr	price_index_und,nom_rate_und,real_rate_und;	
	char	    *nom_rate_und_name,*nom_rate_yc;
	double		df_nom_ipc_date_hm,df_nom_ipc_date_lm,df_real_settlt,df_real_ipc_date_lm,df_real_ipc_date_hm;
	double		ipc_lm,ipc_hm;

	/* Calibrate the mean reversion level given the initial condition */
	err = srt_f_vasicek_calibrated_mean_rev_lvl(settlt_date,
												vasicek_init_cond,
												num_bonds,
												num_cpns,
												pay_dates,
												cpns,

												bond_dirty_prices,
												discount_fwd_dirty_pr,

												price_floor,
												floor_gearing,
												pre_factor_fwd,
												floor_strike,
												floor_implied_vol,

												real_rate_und_name,
												calibrated_mean_rev_lvl);
	if(err) return err;


	/* Get the vasicek sort pointer */
	real_rate_und = lookup_und(real_rate_und_name);
	if(!real_rate_und)
	{
		return serror("can not find vasicek underlying in srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func");

	}

	/* Ge the price index underlying pointer */
	price_index_und = lookup_und(price_und_name);
	if(!price_index_und)
	{
		return serror("can not find price index und in srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func");
	}

	nom_rate_und_name = get_domname_from_fxund(price_index_und);
	nom_rate_und = lookup_und(nom_rate_und_name);
	if(!nom_rate_und) 
	{
		return serror("can not get the nom. rate underlying");
	}
	/* get the nominal rate yield curve */
	err = get_underlying_discname(nom_rate_und, &nom_rate_yc);
	if (err) 
	{
		return err;
	}

	df_nom_ipc_date_lm = swp_f_df(settlt_date,ipc_date_lm,nom_rate_yc);
	df_nom_ipc_date_hm = swp_f_df(settlt_date,ipc_date_hm,nom_rate_yc);

	/* get the forward real rate discount factor */
	err = srt_f_vasicek_discount_factor((long)settlt_date,
										real_rate_und,
										&df_real_settlt);
	if(err) return err;

	err = srt_f_vasicek_discount_factor(ipc_date_lm,
										real_rate_und,
										&df_real_ipc_date_lm);
	if(err) return err;

	err = srt_f_vasicek_discount_factor(ipc_date_hm,
										real_rate_und,
										&df_real_ipc_date_hm);
	if(err) return err;

	ipc_lm = (df_real_ipc_date_lm/df_real_settlt)/df_nom_ipc_date_lm;
	ipc_hm = (df_real_ipc_date_hm/df_real_settlt)/df_nom_ipc_date_hm;
					
	(*model_daily_ref_index) = get_spot_from_fxund(price_index_und)*(ipc_lm + interp_coeff*(ipc_hm-ipc_lm));

	return NULL;

}

Err srt_f_vasicek_calibrated_mean_rev_lvl_and_initial_cond(Date			settlt_date,
														   Date			ipc_date_hm,
														   Date			ipc_date_lm,
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
														   double		*calibrated_vasicek_init_cond)
{
	Err			err = NULL;
	SrtUndPtr	und;
	TermStruct	*ts;
	double	    vasicek_init_cond;
	double		a[3],b[3],nstop;
	long		k;
	long		date_today;


	und = lookup_und(vasicek_und_name);
	if(!und) 
		return serror("can not get underlying in srt_f_vasicek_calibrated_mean_rev_lvl_and_initial_cond");

	err = get_underlying_ts(und,&ts);
	if(err) return err;

	date_today = get_today_from_underlying(und);

	err = srt_f_get_vasicek_init_cond((settlt_date-date_today)*YEARS_IN_DAY,
									  ts,
									  &vasicek_init_cond);


	a[0] = vasicek_init_cond;
	err = srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func(settlt_date,
																 a[0],	
																 
																 ipc_date_hm,
																 ipc_date_lm,
																 interp_coeff,
																 daily_ref_index_target,

																 num_bonds,
																 num_cpns,
																 pay_dates,
																 cpns,
																 bond_dirty_prices,
																 discount_fwd_dirty_pr,

																 price_floor,
															     floor_gearing,
															     pre_factor_fwd,
															     floor_strike,
															     floor_implied_vol,
														
																 vasicek_und_name,
																 price_und_name,
																 calibrated_mean_rev_lvl,
																 &b[0]);
	if(err) return err;

	a[1] = a[0] + 1.0/100;
	err = srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func(settlt_date,
																 a[1],	
																 
																 ipc_date_hm,
																 ipc_date_lm,
																 interp_coeff,
																 daily_ref_index_target,

																 num_bonds,
																 num_cpns,
																 pay_dates,
																 cpns,
																 bond_dirty_prices,
																 discount_fwd_dirty_pr,

																 price_floor,
															     floor_gearing,
															     pre_factor_fwd,
															     floor_strike,
															     floor_implied_vol,
														
																 vasicek_und_name,
																 price_und_name,
																 calibrated_mean_rev_lvl,
																 &b[1]);
	if(err) return err;

	a[2] = a[1] + 1.0/100;
	err = srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func(settlt_date,
																 a[2],	
																 
																 ipc_date_hm,
																 ipc_date_lm,
																 interp_coeff,
																 daily_ref_index_target,

																 num_bonds,
																 num_cpns,
																 pay_dates,
																 cpns,
																 bond_dirty_prices,
																 discount_fwd_dirty_pr,

																 price_floor,
															     floor_gearing,
															     pre_factor_fwd,
															     floor_strike,
															     floor_implied_vol,
														
																 vasicek_und_name,
																 price_und_name,
																 calibrated_mean_rev_lvl,
																 &b[2]);
	if(err) return err;

	nstop = 0.0;
	k = 0;
	while(nstop < 1 && k <MAX_ITER)
	{
		newton(daily_ref_index_target,2.0,a,b,&nstop);

		err = srt_f_vasicek_mean_rev_lvl_and_initial_cond_calib_func(settlt_date,
																 a[2],	
																 
																 ipc_date_hm,
																 ipc_date_lm,
																 interp_coeff,
																 daily_ref_index_target,

																 num_bonds,
																 num_cpns,
																 pay_dates,
																 cpns,
																 bond_dirty_prices,
																 discount_fwd_dirty_pr,

																 price_floor,
															     floor_gearing,
															     pre_factor_fwd,
															     floor_strike,
															     floor_implied_vol,
														
																 vasicek_und_name,
																 price_und_name,
																 calibrated_mean_rev_lvl,
																 &b[2]);
		if(err) return err;

		k++;
	}

	(*calibrated_vasicek_init_cond) = a[2];
	return NULL;
}


#endif