/*	
	author:		Antoine Savine
	file:		srt_f_lgmresetcap.c
	purpose:	pricing of reset caps
	models:		all
*/
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmresetcap.h"
#include "opfnctns.h"

/*	
	part 1
	pricing in LGM 1F and 2F
	procedure:
	1:	find joint distribution of rt (ie. -log(df)) at strike and spot fixing
		under appropriate probability measure (cms effect)
	2:	denoting Y rt at strike and X rt at spot, write the payoff of reset caplet as
		(a*exp(X) - b*exp(Y) - k)+
	3:	price caplets using spread numeraire (srt_f_optsprnum)
	4:	for full caps: strip dates and call caplet pricing function
*/

/*	static utility functions used to find distribution of rt */

static double		temp_swap;
#define				dswap(a,b)		temp_swap=a;a=b;b=temp_swap;
#define				START			0
#define				END				1
#define				FIRST			0
#define				SECOND			1

/*	compute phi function(s) at fixing date
	model: LGM 1F and 2F
	undptr must point to ir underlying of type LGM */
Err		srt_f_lgm_phifunc	(SrtUndPtr		undptr,
							 Ddate			fixing,
							 double			phi[2][2])
{
	long		today_date;
	double		time;
	double		f, g, temp;
	SrtTFTSVec	tf_f;
	SrtTFTSMat	tf_g;
	
	SrtMdlDim	mdldim;
	TermStruct	*ts;
	Err			err;

	/* convert fixing date to time from today */
	today_date = get_today_from_underlying (undptr);
	time = (fixing - today_date) * YEARS_IN_DAY;
	
	err = get_underlying_ts (undptr, &ts);
	err = get_underlying_mdldim (undptr, &mdldim);
	if (mdldim == ONE_FAC)
	{
		f = F_func (time, ts);
		G_H_func (time, ts, &g, &temp);
	
		phi[FIRST][FIRST] = f * f * g;
	}
	else if (mdldim == TWO_FAC)
	{
		err = get_2f_F_funcs (time, ts, &tf_f);
		err = get_2f_G_funcs (time, ts, &tf_g);
		phi[FIRST][FIRST] = tf_f[FIRST] * tf_f[FIRST] * tf_g[FIRST][FIRST];
		phi[SECOND][SECOND] = tf_f[SECOND] * tf_f[SECOND] * tf_g[SECOND][SECOND];
		phi[FIRST][SECOND] = tf_f[FIRST] * tf_f[SECOND] * tf_g[FIRST][SECOND];
		phi[SECOND][FIRST] = phi[FIRST][SECOND];
	}

	return NULL;
}

/*	compute moments of spot rate at fixing under Qpay measure 
	model: LGM 1F and 2F
	undptr must point to ir underlying of type LGM */
static Err		srt_f_lgm_spotratemoments(
							SrtUndPtr		undptr,
							Ddate			fixing,
							Ddate			pay,
							double			exp[2],
							double			var[2][2])
{
	long		today_date;
	double		time_fix, time_pay;
	double		f_fix, g_fix, psi_fix, psi_pay, temp;
	SrtTFTSVec	tf_f_fix, tf_psi_fix, tf_psi_pay;
	SrtTFTSMat	tf_g_fix;

	SrtMdlDim	mdldim;
	TermStruct	*ts;
	Err			err;

	/* convert dates to time from today */
	today_date = get_today_from_underlying (undptr);
	time_fix = (fixing - today_date) * YEARS_IN_DAY;
	time_pay = (pay - today_date) * YEARS_IN_DAY;
	
	err = get_underlying_ts (undptr, &ts);	
	err = get_underlying_mdldim (undptr, &mdldim);
	if (mdldim == ONE_FAC)
	{
		f_fix = F_func (time_fix, ts);
		G_H_func (time_fix, ts, &g_fix, &temp);
		psi_fix = Psi_func (time_fix, ts);
		psi_pay = Psi_func (time_pay, ts);

		exp[FIRST] = - (psi_pay - psi_fix) * f_fix * g_fix;
		var[FIRST][FIRST] = f_fix * f_fix * g_fix;
	}
	else if (mdldim == TWO_FAC)
	{
		err = get_2f_F_funcs (time_fix, ts, &tf_f_fix);
		err = get_2f_G_funcs (time_fix, ts, &tf_g_fix);
		err = get_2f_Psi_funcs (time_fix, ts, &tf_psi_fix);
		err = get_2f_Psi_funcs (time_pay, ts, &tf_psi_pay);

		exp[FIRST] = - (tf_psi_pay[FIRST] - tf_psi_fix[FIRST]) 
					 * tf_f_fix[FIRST] 
					 * tf_g_fix[FIRST][FIRST]
					 - (tf_psi_pay[SECOND] - tf_psi_fix[SECOND])
					 * tf_f_fix[FIRST] 
					 * tf_g_fix[FIRST][SECOND];
		exp[SECOND] = - (tf_psi_pay[SECOND] - tf_psi_fix[SECOND]) 
					 * tf_f_fix[SECOND] 
					 * tf_g_fix[SECOND][SECOND]
					 - (tf_psi_pay[FIRST] - tf_psi_fix[FIRST])
					 * tf_f_fix[SECOND] 
					 * tf_g_fix[SECOND][FIRST];
		var[FIRST][FIRST] =   tf_f_fix[FIRST] 
			                * tf_f_fix[FIRST] 
						    * tf_g_fix[FIRST][FIRST];
		var[SECOND][SECOND] =   tf_f_fix[SECOND] 
							  * tf_f_fix[SECOND] 
							  * tf_g_fix[SECOND][SECOND];
		var[FIRST][SECOND] =   tf_f_fix[FIRST] 
							 * tf_f_fix[SECOND] 
							 * tf_g_fix[FIRST][SECOND];
		var[SECOND][FIRST] = var[FIRST][SECOND];
	}

	return NULL;
}

/*	compute moments of R.T (zero rate times time period)
	fixed for period start_end under Qpay measure  
	model: LGM 1F and 2F
	undptr must point to ir underlying of type LGM */
static Err		srt_f_lgm_zeroratemoments(
					SrtUndPtr		undptr,
					Ddate			fixing,
					Ddate			start_end[2],
					Ddate			pay,
					double			*exp,
					double			*var)
{
	YTatt_param YT[2];
	double		spot_exp[2], spot_var[2][2], phi[2][2];
	SrtMdlDim	mdldim;
	Err			err;

	/* get spot rate moments and phi */
	err = srt_f_lgm_spotratemoments (undptr, fixing, pay, spot_exp, spot_var); 
	if (err)
	{
		return serror (err);
	}
	err = srt_f_lgm_phifunc	(undptr, fixing, phi);
	if (err)
	{
		return serror (err);
	}

	/* get reconstruction coefficients */
	err = Y_T_at_t_param (fixing, start_end, 2, undptr, YT);
	if (err)
	{
		return serror (err);
	}
	
	err = get_underlying_mdldim (undptr, &mdldim);
	if (mdldim == ONE_FAC)
	{
		*exp =   YT[END].fwd_zero_rate 
			   - YT[START].fwd_zero_rate 
			   + (YT[END].phi_coeff[FIRST][FIRST] - YT[START].phi_coeff[FIRST][FIRST])
			   * phi[FIRST][FIRST]
			   + (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST])
			   * spot_exp[FIRST];
		*var =   spot_var[FIRST][FIRST] 
			   * (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST])		  		
			   * (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST]);
	}
	else if (mdldim == TWO_FAC)
	{
				
		*exp =   YT[END].fwd_zero_rate 
			   - YT[START].fwd_zero_rate 
			   
			   + (YT[END].phi_coeff[FIRST][FIRST] - YT[START].phi_coeff[FIRST][FIRST])
			   * phi[FIRST][FIRST]
			   
			   + (YT[END].phi_coeff[SECOND][SECOND] - YT[START].phi_coeff[SECOND][SECOND])
			   * phi[SECOND][SECOND]
			   
			   + (YT[END].phi_coeff[FIRST][SECOND] - YT[START].phi_coeff[FIRST][SECOND])
			   * phi[FIRST][SECOND]
			   
			   + (YT[END].phi_coeff[SECOND][FIRST] - YT[START].phi_coeff[SECOND][FIRST])
			   * phi[SECOND][FIRST]
			   
			   + (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST])
			   * spot_exp[FIRST] 
			   
			   + (YT[END].x_coeff[SECOND] - YT[START].x_coeff[SECOND])
			   * spot_exp[SECOND];
		
		*var =   spot_var[FIRST][FIRST] 
			   * (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST])		  
			   * (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST])
		
			   + spot_var[SECOND][SECOND] 
			   * (YT[END].x_coeff[SECOND] - YT[START].x_coeff[SECOND])		  
			   * (YT[END].x_coeff[SECOND] - YT[START].x_coeff[SECOND])
			   
			   + spot_var[FIRST][SECOND] 
			   * (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST])		  
			   * (YT[END].x_coeff[SECOND] - YT[START].x_coeff[SECOND])
			   
			   + spot_var[SECOND][FIRST] 
			   * (YT[END].x_coeff[SECOND] - YT[START].x_coeff[SECOND])		  
			   * (YT[END].x_coeff[FIRST] - YT[START].x_coeff[FIRST]);
	}

	return NULL;
}

/*	compute covariance of spot rate(s) at different times
	model: LGM 1F and 2F
	undptr must point to ir underlying of type LGM */
static Err		srt_f_lgm_spotratecovar			 (
					SrtUndPtr		undptr,
					Ddate			fix1,
					Ddate			fix2,
					double			covar[2][2])
{
	long		today_date;
	double		time1, time2;
	double		f1, f2, g1, temp;
	SrtTFTSVec	tf_f1, tf_f2;
	SrtTFTSMat	tf_g1;
	
	SrtMdlDim	mdldim;
	TermStruct	*ts;
	Err			err;

	/* make sure that fix1 <= fix2 */
	if (fix1 > fix2)
	{
		dswap (fix1, fix2);
	}

	/* convert dates to time from today */
	today_date = get_today_from_underlying (undptr);
	time1 = (fix1 - today_date) * YEARS_IN_DAY;
	time2 = (fix2 - today_date) * YEARS_IN_DAY;

	err = get_underlying_ts (undptr, &ts);
	err = get_underlying_mdldim (undptr, &mdldim);
	if (mdldim == ONE_FAC)
	{
		f1 = F_func (time1, ts);
		f2 = F_func (time2, ts);
		G_H_func (time1, ts, &g1, &temp);

		covar[FIRST][FIRST] = f1 * f2 * g1;
	}
	else if (mdldim == TWO_FAC)
	{
		err = get_2f_F_funcs (time1, ts, &tf_f1);
		err = get_2f_F_funcs (time2, ts, &tf_f2);
		err = get_2f_G_funcs (time1, ts, &tf_g1);

		covar[FIRST][FIRST]   = tf_f1[FIRST]  * tf_f2[FIRST]  * tf_g1[FIRST][FIRST];
		covar[SECOND][SECOND] = tf_f1[SECOND] * tf_f2[SECOND] * tf_g1[SECOND][SECOND];
		covar[FIRST][SECOND]  = tf_f1[FIRST]  * tf_f2[SECOND] * tf_g1[FIRST][SECOND];
		covar[SECOND][FIRST]  = tf_f1[SECOND] * tf_f2[FIRST]  * tf_g1[SECOND][FIRST];
	}

	return NULL;
}

/*	compute covariance of two different R.T   
	model: LGM 1F and 2F
	undptr must point to ir underlying of type LGM */
static Err		srt_f_lgm_zeroratecovar			 (
												  SrtUndPtr		undptr,
												 Ddate			fix1,
												 Ddate			fix2,
												 Ddate			start_end1[2],
												 Ddate			start_end2[2],
												 double			*covar)
{
	YTatt_param YT[2][2];
	double		spot_covar[2][2];
	SrtMdlDim	mdldim;
	Err			err;

	/* make sure that fix1 <= fix2 */
	if (fix1 > fix2)
	{
		dswap (fix1, fix2);
		dswap (start_end1[0], start_end2[0]);
		dswap (start_end1[1], start_end2[1]);
	}
	
	/* get spot rate covariance */
	err = srt_f_lgm_spotratecovar (undptr, fix1, fix2, spot_covar);	
	if (err)
	{
		return serror (err);
	}

	/* get reconstruction coefficients */
	err = Y_T_at_t_param (fix1, start_end1, 2, undptr, YT[0]);
	if (err)
	{
		return serror (err);
	}
	err = Y_T_at_t_param (fix2, start_end2, 2, undptr, YT[1]);
	if (err)
	{
		return serror (err);
	}

	err = get_underlying_mdldim (undptr, &mdldim);
	if (mdldim == ONE_FAC)
	{
		*covar =   spot_covar[FIRST][FIRST] 
			     * (YT[0][END].x_coeff[FIRST] - YT[0][START].x_coeff[FIRST])		  
			     * (YT[1][END].x_coeff[FIRST] - YT[1][START].x_coeff[FIRST]);		  
	}
	else if (mdldim == TWO_FAC)
	{
		*covar =   spot_covar[FIRST][FIRST] 
					 * (YT[0][END].x_coeff[FIRST] - YT[0][START].x_coeff[FIRST])		  
					 * (YT[1][END].x_coeff[FIRST] - YT[1][START].x_coeff[FIRST])

				 + spot_covar[SECOND][SECOND] 
					 * (YT[0][END].x_coeff[SECOND] - YT[0][START].x_coeff[SECOND])		  
					 * (YT[1][END].x_coeff[SECOND] - YT[1][START].x_coeff[SECOND])
				 
				 + spot_covar[FIRST][SECOND] 
					 * (YT[0][END].x_coeff[FIRST] - YT[0][START].x_coeff[FIRST])		  
					 * (YT[1][END].x_coeff[SECOND] - YT[1][START].x_coeff[SECOND])		  
				 
				 + spot_covar[SECOND][FIRST] 
					 * (YT[0][END].x_coeff[SECOND] - YT[0][START].x_coeff[SECOND])		  
					 * (YT[1][END].x_coeff[FIRST] - YT[1][START].x_coeff[FIRST]);		  
	}

	return NULL;
}

/*	main function for reset caplet pricing in LGM */

/*	price reset caplet in LGM using srt_f_optsprnum
	dates, coverages and level are to be set up before use */
Err		srt_f_lgm_resetcaplet (SrtUndPtr		undptr, 
									   Ddate			str_fix,
									   Ddate			str_start_end[2],
									   double			str_cvg,
									   Ddate			spot_fix,
									   Ddate			spot_start_end[2],
									   double			spot_cvg,
									   Ddate			pay_date,	
									   double			pay_level,
									   SrtReceiverType	rec_pay,
									   double			*premium,
									   double			fixed_strike)
{
	double			str_frt_exp, str_frt_var, spot_frt_exp, spot_frt_var, str_spot_frt_covar;
	double			a, b, k, Fx, Fy, sx, sy, rho;
	Err			err;

	/* moments of forward rates */
	/* strike */
	err = srt_f_lgm_zeroratemoments (undptr, str_fix, str_start_end, pay_date,
									 &str_frt_exp, &str_frt_var);
	if (err)
	{
		return serror (err);
	}
	/* spot */
	err = srt_f_lgm_zeroratemoments (undptr, spot_fix, spot_start_end, pay_date,
									 &spot_frt_exp, &spot_frt_var);
	if (err)
	{
		return serror (err);
	}
	/* cross */
	err = srt_f_lgm_zeroratecovar (undptr, str_fix, spot_fix, 
								   str_start_end, spot_start_end,
								   &str_spot_frt_covar);
	if (err)
	{
		return serror (err);
	}

	/* arguments for srt_f_optsprnum */
	a = 1.0 / spot_cvg;
	b = 1.0 / str_cvg;
	k = 1.0 / spot_cvg - 1.0 / str_cvg;
	Fx = exp (spot_frt_exp + 0.5 * spot_frt_var);
	Fy = exp (str_frt_exp + 0.5 * str_frt_var);
	sx = sqrt (spot_frt_var);
	sy = sqrt (str_frt_var);

	if( (sx == 0) || (sy == 0) )
		return serror("Fatal error in srt_f_lgmresetcap() ");

	rho = str_spot_frt_covar / sx / sy;

	/* fixed strike case */
	if (fixed_strike > EPS)
	{
		sy = 0.0;
		rho = 0.0;
		Fy = str_cvg * fixed_strike + 1;
	}
		
	/* main call */ 
	if (sy > EPS)
	/* general case: strike Libor is not deterministic, call srt_f_optsprnum */
	{
		if (rec_pay == SRT_STRADDLE)
		{
			*premium = srt_f_optsprnum (Fx, Fy, sx, sy, rho, 
									   1.0, 
									   a, b, k, 
									   1.0, 
									   SRT_CALL, PREMIUM)
					+ srt_f_optsprnum (Fx, Fy, sx, sy, rho, 
									   1.0, 
									   a, b, k, 
									   1.0, 
									   SRT_PUT, PREMIUM);
		}
		else if (rec_pay == SRT_RECEIVER)
		{
			*premium = srt_f_optsprnum (Fx, Fy, sx, sy, rho, 
									   1.0, 
									   a, b, k, 
									   1.0, 
									   SRT_PUT, PREMIUM);
		}
		else if (rec_pay == SRT_PAYER)
		{
			*premium = srt_f_optsprnum (Fx, Fy, sx, sy, rho, 
									   1.0, 
									   a, b, k, 
									   1.0, 
									   SRT_CALL, PREMIUM);
		}
		else
		{
			return serror (" couldn't interpret rec_pay ");
		}
	}
	else
	/* first caplet: strike is fixed now and therefore is deterministic, 
	   call srt_f_optblksch */	
	{
		if (rec_pay == SRT_STRADDLE)
		{
			*premium = srt_f_optblksch (a * Fx, b * Fy + k, sx, 
									   1.0, 
									   1.0, 
									   SRT_CALL, PREMIUM)
					 + srt_f_optblksch (a * Fx, b * Fy + k, sx,  
									   1.0, 
									   1.0, 
									   SRT_PUT, PREMIUM);
		}
		else if (rec_pay == SRT_RECEIVER)
		{
			*premium = srt_f_optblksch (a * Fx, b * Fy + k, sx, 
									   1.0, 
									   1.0, 
									   SRT_PUT, PREMIUM);
		}
		else if (rec_pay == SRT_PAYER)
		{
			*premium = srt_f_optblksch (a * Fx, b * Fy + k, sx, 
									   1.0, 
									   1.0, 
									   SRT_CALL, PREMIUM);
		}
		else
		{
			return serror (" couldn't interpret rec_pay ");
		}
	}

	/* discount */
	*premium *= pay_level;

	return NULL;
}

/*	main function for reset cap pricing in LGM */

/*	price reset cap in LGM using srt_f_optsprnum */
double	gen_ir_resetcap		(SrtUndPtr			undptr, 
							 long				*fixing_date, 
							 long				*start_end_pay_dts,
							 double				*df, 
							 double				*cvg, 
							 long				nfp,
							 SrtReceiverType	rec_pay,
							 double				first_fixing)
{
	long		i;
	double		prem, temp;
	Date		today_date;
	Ddate		str_start_end[2], spot_start_end[2], pay;
	Err			err;

	today_date = get_today_from_underlying (undptr);	

	prem = 0.0;

	/* first caplet */
	if (fixing_date[0] >= today_date )
	{
		str_start_end[START] = (Ddate) start_end_pay_dts[0];
		str_start_end[END] = (Ddate) start_end_pay_dts[1];
		spot_start_end[START] = (Ddate) start_end_pay_dts[1];
		spot_start_end[END] = (Ddate) start_end_pay_dts[2];
		pay = (Ddate) start_end_pay_dts[1];

		err = srt_f_lgm_resetcaplet	   (undptr, 
									   (Ddate) fixing_date[0], 
									   str_start_end,
									   cvg[1],
									   fixing_date[1],
									   spot_start_end,
									   cvg[2],
									   pay,	
									   df[1] * cvg[1],
									   rec_pay,
									   &temp,
									   first_fixing);

		if (err)
		{
			smessage (err);
			return ERROR_IN_RESETCAP;
		}

		prem += temp;
	}
	/* next caplets */
	for (i=1; i<nfp-2; i++)
	{
		if (fixing_date[i] >= today_date )
		{
			str_start_end[START] = (Ddate) start_end_pay_dts[i];
			str_start_end[END] = (Ddate) start_end_pay_dts[i+1];
			spot_start_end[START] = (Ddate) start_end_pay_dts[i+1];
			spot_start_end[END] = (Ddate) start_end_pay_dts[i+2];
			pay = (Ddate) start_end_pay_dts[i+1];

			err = srt_f_lgm_resetcaplet	   (undptr, 
										   (Ddate) fixing_date[i], 
										   str_start_end,
										   cvg[i+1],
										   fixing_date[i+1],
										   spot_start_end,
										   cvg[i+2],
										   pay,	
										   df[i+1] * cvg[i+1],
										   rec_pay,
										   &temp,
										   0.0);

			if (err)
			{
				smessage (err);
				return ERROR_IN_RESETCAP;
			}

			prem += temp;
		}
	}

	return prem;
}

/*	
	part 2
	pricing in Cheyette 1F and 2F
	procedure: make grfn tableau (grfn_SwapDP_to_GrfnCells)
*/

/*
	part 3
	pricing in Black-Scholes
	to be implemented in srt_f_bsresetcap.c
*/
