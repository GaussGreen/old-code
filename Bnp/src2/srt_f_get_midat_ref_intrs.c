#include "math.h"
#include "OPFNCTNS.H>
#include "srt_h_get_midat_ref_instrs.h"

#define NBR_REIMANN_POINTS 10
#define BUTTERFLY_SHIFT 1e-4

Err srt_f_get_likely_most_expensive_european(Date				 *midat_start_dates,
											 Date				 midat_theo_end_date,
											 String				 ref_rate_code,
											 SrtCompounding		 fixed_freq, 
											 BasisCode			 fixed_basis,
											 long				 nEx,
											 Err				(*get_vol)(Date, Date, double, SRT_Boolean, double*),
											 String				 yc_name,
											 SrtDiffusionType	 swp_vol_type,
											 char				 *pay_rec_str,
											 double				 *swp_strikes,

											 double				 **swp_fwd_vol,
											 
											 long				 **likely_most_expensive)
{
	Err				 err = NULL;
	double			 dx,*swp_spread,*cvg,ref_swap_rate,swap_rate,yr_to_exp_from_exer_date,yr_to_midat_exer;
	long			 i,j,k,l,jmax,num_dates,nfp,num_pay_dates;
	SwapDP			 ref_swp_dp,swp_dp;
	Date			 today,*pay_dates,*start_dates,*end_dates,*midat_exer_dates;	
	SrtCurvePtr		 yld_crv;
	long			 spot_lag;
	double           **swp_premium,*tex_ref_swp_rate,zero_cpn_cash,level_cash,prob_dens_func,swp_vol;
	SrtReceiverType	 srt_pay_rec;
	double			 *ex_frontier_bounds;

	yld_crv = lookup_curve(yc_name);

	today  = get_clcndate_from_yldcrv(yld_crv);
	spot_lag  = get_spotlag_from_curve(yld_crv);

	err = interp_rec_pay (pay_rec_str, &srt_pay_rec);	
	if (err) return err;

	swp_spread = dvector(1, nEx-1);
	swp_premium = dmatrix(1,nEx,1,nEx-1);
	midat_exer_dates = lngvector(1,nEx);
	ex_frontier_bounds = dvector(1,nEx);
	tex_ref_swp_rate = dvector(1,NBR_REIMANN_POINTS);

	for( i = 1; i <= nEx; i++) midat_exer_dates[i] = add_unit(midat_start_dates[i],-spot_lag,SRT_BDAY, SUCCEEDING);

	for( i = 1; i <= nEx; i++)
	{
		if(srt_pay_rec == SRT_RECEIVER)
		{
			ex_frontier_bounds[i] = 0.0;
			dx = (swp_strikes[i]-ex_frontier_bounds[i])/NBR_REIMANN_POINTS; 
		}
		else 
		{
			ex_frontier_bounds[i] = swp_strikes[i] + swp_strikes[i];
			dx = (ex_frontier_bounds[i]-swp_strikes[i])/NBR_REIMANN_POINTS; 
		}
		yr_to_midat_exer = (double)(midat_exer_dates[i]-today)*YEARS_IN_DAY;

		err = swp_f_setSwapDP(midat_start_dates[i], 
						  midat_theo_end_date, 
						  fixed_freq, 
						  fixed_basis, 
						  &ref_swp_dp);
		if(err) return err;

		err =  swp_f_ForwardRate_SwapDP(&ref_swp_dp,
										yc_name,
										ref_rate_code,
										&ref_swap_rate);
		if(err) return err;

		for( j = 1; j <= (nEx-i); j++)
		{
			err = swp_f_setSwapDP(midat_start_dates[i+j], 
							  midat_theo_end_date, 
							  fixed_freq, 
							  fixed_basis, 
							  &swp_dp);
			
			if(err) return err;

			err =  swp_f_ForwardRate_SwapDP(&swp_dp,
											yc_name,
											ref_rate_code,
											&swap_rate);
			if(err) return err;

			swp_spread[j] = (swap_rate-ref_swap_rate);
							
			err = swp_f_make_FixedLegDatesAndCoverages(&swp_dp,
													   today,
													   &pay_dates,
													   &num_pay_dates,
													   &start_dates,
													   &end_dates,
													   &cvg,
													   &num_dates);
			if(err) return err;

			nfp = num_dates-1;
			yr_to_exp_from_exer_date = (double)(midat_exer_dates[j+i]-midat_exer_dates[i])*YEARS_IN_DAY;

			swp_premium[i][j] = 0.;
			for(k = 1; k <= NBR_REIMANN_POINTS; k++) /* LOOP ON THE NUMBER OF DISCRETISATION POINTS */
			{
				if(srt_pay_rec == SRT_RECEIVER) tex_ref_swp_rate[k] = ex_frontier_bounds[i] + k*dx;
				else tex_ref_swp_rate[k] = swp_strikes[i] + k*dx;
				
				zero_cpn_cash = 1.0/(1 + tex_ref_swp_rate[k]/(double)fixed_freq);

				level_cash = 0.0;
				for(l = nfp; l >=0; l--) level_cash += cvg[l]*pow(zero_cpn_cash,l); /* has to be changed*/

				prob_dens_func = 0.0;
				if(swp_vol_type == SRT_LOGNORMAL)
				{
					err = get_vol(midat_start_dates[i],
								  midat_theo_end_date,
								  (tex_ref_swp_rate[k]+BUTTERFLY_SHIFT),
								  SRT_FALSE,
								  &swp_vol);	

					if(err) return err;

					prob_dens_func = srt_f_optblksch(ref_swap_rate,
													(tex_ref_swp_rate[k]+BUTTERFLY_SHIFT),
													swp_vol,
													yr_to_midat_exer,
													1.0,
													srt_pay_rec,
													SRT_PREMIUM);

					err = get_vol(midat_start_dates[i],
								  midat_theo_end_date,
								  tex_ref_swp_rate[k],
								  SRT_FALSE,
								  &swp_vol);	

					if(err) return err;

					prob_dens_func -= 2*srt_f_optblksch(ref_swap_rate,
													tex_ref_swp_rate[k],
													swp_vol,
													yr_to_midat_exer,
													1.0,
													srt_pay_rec,
													SRT_PREMIUM);

					err = get_vol(midat_start_dates[i],
								  midat_theo_end_date,
								  (tex_ref_swp_rate[k]-BUTTERFLY_SHIFT),
								  SRT_FALSE,
								  &swp_vol);	

					if(err) return err;

					prob_dens_func += srt_f_optblksch(ref_swap_rate,
													(tex_ref_swp_rate[k]-BUTTERFLY_SHIFT),
													swp_vol,
													yr_to_midat_exer,
													1.0,
													srt_pay_rec,
													SRT_PREMIUM);

					prob_dens_func /= (BUTTERFLY_SHIFT*BUTTERFLY_SHIFT);

					prob_dens_func *= dx;

					swp_premium[i][j] +=  prob_dens_func*level_cash*srt_f_optblksch(k*dx + swp_spread[j],
																	swp_strikes[i+j],
																	swp_fwd_vol[i][j],
																	yr_to_exp_from_exer_date,
																	1,
																	srt_pay_rec,
																	SRT_PREMIUM);
				}
				else
				{

					err = get_vol(midat_start_dates[i],
								  midat_theo_end_date,
								  (tex_ref_swp_rate[k]+BUTTERFLY_SHIFT),
								  SRT_FALSE,
								  &swp_vol);	

					if(err) return err;

					prob_dens_func = srt_f_optblknrm(ref_swap_rate,
													(tex_ref_swp_rate[k]+BUTTERFLY_SHIFT),
													swp_vol,
													yr_to_midat_exer,
													1.0,
													srt_pay_rec,
													SRT_PREMIUM);

					err = get_vol(midat_start_dates[i],
								  midat_theo_end_date,
								  tex_ref_swp_rate[k],
								  SRT_FALSE,
								  &swp_vol);	

					if(err) return err;

					prob_dens_func -= 2*srt_f_optblknrm(ref_swap_rate,
													tex_ref_swp_rate[k],
													swp_vol,
													yr_to_midat_exer,
													1.0,
													srt_pay_rec,
													SRT_PREMIUM);

					err = get_vol(midat_start_dates[i],
								  midat_theo_end_date,
								  (tex_ref_swp_rate[k]-BUTTERFLY_SHIFT),
								  SRT_FALSE,
								  &swp_vol);	

					if(err) return err;

					prob_dens_func += srt_f_optblknrm(ref_swap_rate,
													(tex_ref_swp_rate[k]-BUTTERFLY_SHIFT),
													swp_vol,
													yr_to_midat_exer,
													1.0,
													srt_pay_rec,
													SRT_PREMIUM);

					prob_dens_func /= (BUTTERFLY_SHIFT*BUTTERFLY_SHIFT);
					
					prob_dens_func *= dx;


					swp_premium[i][j] +=  prob_dens_func*level_cash*srt_f_optblknrm(k*dx + swp_spread[j],
																	swp_strikes[i+j],
																	swp_fwd_vol[i][j],
																	yr_to_exp_from_exer_date,
																	1,
																	srt_pay_rec,
																	SRT_PREMIUM);
				}

			}/* END OF LOOP ON k */
		}/* END OF LOOP ON j */	

		/* GET THE INDEX j SUCH THAT swp_premium[i][j] IS MAXIMUM FOR EACH i */
		jmax = 1;
		for( j = 1; j <= (nEx-i); j++)
		{
			if(swp_premium[i][j] > swp_premium[i][jmax]) jmax = j;
			
		}
		(*likely_most_expensive)[i] = jmax;


	} /* END OF LOOP ON i */

	if(swp_spread) free_dvector(swp_spread,1, nEx-1); swp_spread = NULL;
	if(swp_premium) free_dmatrix(swp_premium,1,nEx,1,nEx-1);swp_premium = NULL;
	if(midat_exer_dates) free_lngvector(midat_exer_dates,1,nEx);midat_exer_dates = NULL;
	if(tex_ref_swp_rate) free_dvector(tex_ref_swp_rate,1,NBR_REIMANN_POINTS); tex_ref_swp_rate = NULL;
	if(ex_frontier_bounds) free_dvector(ex_frontier_bounds,1,nEx); ex_frontier_bounds = NULL;

	return err;

}

#undef NBR_REIMANN_POINTS
#undef BUTTERFLY_SHIFT 

