/* ==========================================================================
   FILE_NAME:	Fx3FCalib.c

   PURPOSE:		Modelling of the spot FX vol by taking in consideration a LGM 1 factor
				on the domestic and foreign market and a lognormal model on the Spot Fx

   DATE:		05/25/00
   
   AUTHOR:		L.C.
   ========================================================================== */

#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "opfnctns.h"
#include "math.h"

#define	MAX_TIME	0.02
#define	SHIFT_VOL	0.005
#define	MAX_ERROR	0.0001
#define	MAX_ERROR2	0.0005
#define	SHIFT_VOL2	0.05

#define MAX_STP 3000

#define	MAX_Z	2.0

/*	Fill the time vector */
static Err fill_time_vector(
double				**time, 
int					*nstp, 
int					num_bar_times,
double				*bar_times, 
int					num_vol_times, 
double				*vol_times, 
int					target_nstp)
{
	Err				err = NULL;

	/*	Add today if required */
	if ((*time)[0] < -EPS)
	{
		err = "Past event date in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	if ((*time)[0] > EPS)
	{
		num_f_add_number (nstp, time, 0.0);
		num_f_sort_vector (*nstp, *time);
		num_f_unique_vector (nstp, *time);
	}

	/*	If only one event today, add empty event */
	if (*nstp == 1)
	{
		num_f_add_number (nstp, time, 1.0);		
	}

	/*	Fill the vector */		
	/*	New algorithm */
	num_f_fill_vector_newalgo (nstp, time, target_nstp);

	/*	Old algorithm */
/*	num_f_fill_vector (*nstp, time, target_nstp, 0);
	*nstp = (*nstp > target_nstp? *nstp: target_nstp);
	num_f_fill_vector_maxtime (nstp, time, 1.1 * (*time)[*nstp-1]/(*nstp-1)); */

FREE_RETURN:

	return err;
}

Err FxQuad_log_approx(	
						long		today,						
						double		*maturity,
						long		nb_mat,
						double		*sig_dom,		
						double		*sig_for, 
						double		*mat_fx,
						long		nb_mat_fx,
						double		*sig_fx,
						double		alpha,
						double		beta,
						double		gamma,
						double		sig0,
						double		lam_dom, 
						double		lam_for, 
						double		corr_dom_for,
						double		corr_dom_fx,
						double		corr_for_fx,
						double		spot,
						char		*dom_yc,
						char		*for_yc,
						double		*time_fx,
						long		nb_time_fx,
						double		*fx_fwd,
						double		*fx_vol,
						double		*fx_fwd0,
						double		max_time)
{
long	i, j, k, k2;
long	index;
double	t, prev_t, dt, T;	
double	logSpot;
double	*time2	= NULL;
double	*S		= NULL,
		*Sfwd	= NULL;
double	volfx;
double	F, sigdom, sigfor, sigfx, x, bet;
double	y, stdev;
long	nb_time;
double	last_mat;
long	last_index;
	
Err		err = NULL;

	logSpot = log(spot);

	/* discretise time */
	last_mat = time_fx[nb_time_fx - 1];
	nb_time = (long) (last_mat / max_time -1.0E-08) + 1;
	if (nb_time < 2)
	{
		nb_time = 2;
	}
	
	time2 = calloc(nb_time, sizeof(double));
	if (!time2)
	{
		err = "Memory allocation failure (1) in fwd_approx";
		goto FREE_RETURN;
	}

	time2[0] = 0.0;

	for (i=1; i<nb_time-1; i++)
	{
		time2[i] = time2[i-1] + max_time;
	}
	time2[i] = last_mat;

	num_f_concat_vector (&nb_time, &time2, nb_mat, maturity);
	num_f_concat_vector (&nb_time, &time2, nb_mat_fx, mat_fx);
	num_f_concat_vector (&nb_time, &time2, nb_time_fx, time_fx);
	num_f_sort_vector (nb_time, time2);	
	num_f_unique_vector(&nb_time, time2);

	/* find the index of the last maturity */
	last_index = Get_Index(last_mat, time2, nb_time);
	
	S = (double *) calloc(last_index+1, sizeof(double));	
	Sfwd = (double *) calloc(last_index+1, sizeof(double));	

	if (!S || !Sfwd)
	{
		err = "Memory allocation failure (2) in fwd_approx";
		goto FREE_RETURN;
	}

	/* diffuse */
	
	
	fx_fwd[0] = logSpot;
	fx_fwd0[0] = spot;

	S[0] = spot;
	Sfwd[0] = spot;
	bet = 1.0 - beta;
			
	for (i=0; i<last_index; i++)
	{		
		t = 0;
		T = time2[i+1];
		/* initialisation at the fwd */
		F = logSpot + 
			(swp_f_zr (today, (long) (today + T * 365.0 + 1.0E-08), dom_yc) - swp_f_zr (today, (long) (today + T * 365.0 + 1.0E-08), for_yc)) * T;		

		Sfwd[i+1] = exp(F);
		k = 0;
		k2 = 0;
		sigdom = sig_dom[0];
		sigfor = sig_for[0];
		sigfx = sig_fx[0];

					
		for (j=0; j<i+1; j++)
		{	
			prev_t = t;
			t = time2[j+1];

			if (t > maturity[k])
			{
				if (k < nb_mat-1)
				{
					k += 1;
				}
				sigdom = sig_dom[k];
				sigfor = sig_for[k];							
				
			}	
			if (t > mat_fx[k2])
			{	
				if (k2 < nb_mat_fx-1)
				{
					k2 += 1;
				}
				sigfx = sig_fx[k2];										
			}	

			
			dt = t - prev_t;
			x = sigdom * (1.0 - exp(-lam_dom * (T - prev_t))) / lam_dom;
			y = Sfwd[j] / S[j];
			stdev = sig0 * sqrt(prev_t);
			F += x * dt * (
							 corr_dom_fx * sigfx * (alpha * stdev * y + (1.0 - y) * (beta + gamma / stdev * (1.0 - y)))
						   - corr_dom_for * sigfor * (1.0 - exp(-lam_for * (T - prev_t))) / lam_for
						   + x
						  );							
		}
		S[i+1] = exp(F);		
	}

	j = 0;
	/* fill the out structure */
	for (i=0; i<nb_time_fx; i++)
	{
		t = time_fx[i];
		index = Get_Index(t, mat_fx, nb_mat_fx);
		volfx = sig_fx[index];
		index = Get_Index(t, time2, last_index+1);
		fx_fwd[i] = S[index];		
		fx_fwd0[i] = Sfwd[index];
		y = fx_fwd0[i] / fx_fwd[i];
		if (t>0)
		{
			stdev = sig0 * sqrt(time_fx[i-1]);
			fx_vol[i] = alpha * stdev * y + (1.0 - y) * (beta + gamma / stdev * (1.0 - y));
		}
		else
		{
		}
	}
		

FREE_RETURN:

	if (time2)
	{
		free (time2);
	}
	
	if (S)
	{
		free (S);
	}

	if (Sfwd)
	{
		free (Sfwd);
	}

	return err;
}


/*	Calibration of a fx term structure to a set of fx options  */
Err Fx3DQuadtsImpliedVol(
						long	today,
						double	opt_maturity,
						double	start_date,
						double	end_date,
						double	*maturity,
						long	nbMat,
						double	*sig_curve_dom,
						double	lda_dom,
						double	*sig_curve_for,
						double	lda_for,						
						double	*mat_fx,
						long	nb_mat_fx,
						double	*sig_curve_fx,
						double	alpha,
						double	beta,
						double	gamma,
						double	sig0,
						double	spot_fx,
						double	corr_dom_for,
						double	corr_dom_fx,
						double	corr_for_fx,
						char	*dom_yc,
						char	*for_yc,						
						double	*fx_vol,
						double	disc_dt,
						double	fx_dt)

{
double	*fx_vol_curve	= NULL,
		*fx_fwd			= NULL,		
		*fx_time		= NULL;


long	i;		
long	nb_fx_time, nb_fx_time2;

Err err = NULL;

	/* simple translation to make everything begin at 0 */
 
	nb_fx_time = (long) (end_date / fx_dt - 1.0E-08) + 1;
	if (nb_fx_time  < 2)
	{
		nb_fx_time = 2;
	}

	fx_time = calloc(nb_fx_time, sizeof(double));
	
	if (!fx_time)
	{
		err = "Memory allocation error in Fx3DBetatsImpliedVol (1)";
		goto FREE_RETURN;
	}

	fx_time[0] = 0.0;
	
	/* fill the time where we want the volatility structure */
	for (i=1; i<nb_fx_time-1; i++)
	{
		fx_time[i] = fx_time[i-1] + fx_dt;
	}

	fx_time[i] = end_date;

	num_f_concat_vector (&nb_fx_time, &fx_time, nbMat, maturity);
	num_f_concat_vector (&nb_fx_time, &fx_time, nb_mat_fx, mat_fx);	
	num_f_sort_vector (nb_fx_time, fx_time);	
	num_f_unique_vector(&nb_fx_time, fx_time);

	nb_fx_time2 = Get_Index(end_date, fx_time, nb_fx_time) + 1;

	fx_vol_curve = calloc(nb_fx_time2, sizeof(double));
	fx_fwd = calloc(nb_fx_time2, sizeof(double));

	if (!fx_vol || !fx_fwd)
	{
		err = "Memory allocation error in Fx3DBetatsImpliedVol (2)";
		goto FREE_RETURN;
	}

	err = FxQuad_log_approx(	
						today,						
						maturity,
						nbMat,
						sig_curve_dom,		
						sig_curve_for, 
						mat_fx,
						nb_mat_fx,
						sig_curve_fx,
						alpha,
						beta,
						gamma,
						sig0,
						lda_dom,
						lda_for,
						corr_dom_for,
						corr_dom_fx,
						corr_for_fx,
						spot_fx,
						dom_yc,
						for_yc,
						fx_time,
						nb_fx_time2,
						fx_fwd,						
						fx_vol_curve,
						fx_fwd,
						disc_dt);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = Fx3DtsImpliedVol(	opt_maturity,
							start_date,
							end_date,
							maturity,
							nbMat,
							sig_curve_dom,
							lda_dom,
							sig_curve_for,
							lda_for,
							fx_time,
							fx_vol_curve,
							nb_fx_time2,							
							corr_dom_for,
							corr_dom_fx,
							corr_for_fx,
							fx_vol);

FREE_RETURN:

	if (fx_vol_curve) free (fx_vol_curve);
	if (fx_time) free (fx_time);
	if (fx_fwd) free (fx_fwd);

	return err;
}


/*	Implied vol direct from underlying */
Err Fx3DQuadImpliedVol(
					   char		*fx_underlying,
					   double	beta,
					   double	gamma,
					   double	sig0,
					   double	val_time,
					   double	start_time,
					   double	end_time,
					   double	disc_dt,
					   double	fx_dt,
					   double	*vol
					   )
{
long		sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx, nb_merge_dates;
long		i;
double		*sigma_date_dom	= NULL, 
			*sigma_dom		= NULL;
double		*tau_date_dom	= NULL, 
			*tau_dom		= NULL;
double		*sigma_date_for	= NULL, 
			*sigma_for		= NULL;
double		*tau_date_for	= NULL,  
			*tau_for		= NULL;
double		*sigma_date_fx	= NULL,  
			*sigma_fx		= NULL;
double		correl_dom_for, correl_dom_fx, correl_for_fx;
double		lda_dom, lda_for;
double		*sig_dom		= NULL,
			*sig_for		= NULL,
			*sig_fx			= NULL,
			*merge_dates	= NULL;

SrtUndPtr	fx_und, dom_und, for_und;
char		*domname, *forname;

long		today;
double		spot_fx;
char		*dom_yc, *for_yc;

double		alpha2, beta2, gamma2;

Err			err = NULL;


	err = Get_FX_StochRate_TermStructures(fx_underlying,
									 &sigma_date_dom,  &sigma_dom,  &sigma_n_dom,
									 &tau_date_dom,  &tau_dom,  &tau_n_dom,
									 &sigma_date_for,  &sigma_for,  &sigma_n_for,
									 &tau_date_for,  &tau_for,  &tau_n_for,
									 &sigma_date_fx,  &sigma_fx,  &sigma_n_fx,									
									 &correl_dom_for,  &correl_dom_fx,  &correl_for_fx);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (tau_dom, tau_n_dom, &lda_dom);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (tau_for, tau_n_for, &lda_for);
	if (err)
	{
		goto FREE_RETURN;
	}

	/* now merge all the term structure */
	merge_dates = (double*) calloc (sigma_n_dom, sizeof (double));
	memcpy (merge_dates, sigma_date_dom, sigma_n_dom * sizeof (double));
	nb_merge_dates = sigma_n_dom;
	num_f_concat_vector (&nb_merge_dates, &merge_dates, sigma_n_for, sigma_date_for);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
	num_f_sort_vector (nb_merge_dates, merge_dates);	
	num_f_unique_vector (&nb_merge_dates, merge_dates);

	/*	Fill the new term structures */	
	sig_dom = (double*) calloc (nb_merge_dates, sizeof (double));
	sig_for = (double*) calloc (nb_merge_dates, sizeof (double));
	sig_fx = (double*) calloc (nb_merge_dates, sizeof (double));

	if (!sig_dom || !sig_for || !sig_fx)
	{
		err = "Memory allocation error (1) in Fx3DBetaImpliedVol";
		goto FREE_RETURN;
	}
	
	for (i=nb_merge_dates-1; i>=0; i--)
	{
		sig_dom[i] = sigma_dom[Get_Index(merge_dates[i], sigma_date_dom, sigma_n_dom)];
		sig_for[i] = sigma_for[Get_Index(merge_dates[i], sigma_date_for, sigma_n_for)];
		sig_fx[i] = sigma_fx[Get_Index(merge_dates[i], sigma_date_fx, sigma_n_fx)];
	}
		
	/* get today, spot, dom_yc, for_yc */
	fx_und = lookup_und (fx_underlying);
	domname = get_domname_from_fxund (fx_und);
	forname = get_forname_from_fxund (fx_und);
	dom_und =lookup_und(domname);
	for_und =lookup_und(forname);
	
	err = get_underlying_discname(dom_und, &dom_yc);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = get_underlying_discname(for_und, &for_yc);
	if (err)
	{
		goto FREE_RETURN;
	}

	today = get_today_from_underlying (dom_und);


	err = Fx3DtsFwdFx(dom_und, for_und, fx_und, today, &spot_fx);


	/* reparametrisation of the model */
	//alpha2 = (1.0 - beta) * spot_fx + gamma * spot_fx * spot_fx;
	alpha2 = 1.0 - beta + gamma;
	//beta2 = beta - 2.0 * gamma * spot_fx;
	beta2 = beta - 2.0 * gamma;
	gamma2 = gamma;

	err = Fx3DQuadtsImpliedVol(
							today,
							val_time,
							start_time,
							end_time,
							merge_dates,
							nb_merge_dates,
							sig_dom,
							lda_dom,
							sig_for,
							lda_for,						
							merge_dates,
							nb_merge_dates,
							sig_fx,
							alpha2,
							beta2,
							gamma2,
							sig0,
							spot_fx,
							correl_dom_for,
							correl_dom_fx,
							correl_for_fx,
							dom_yc,
							for_yc,						
							vol,
							disc_dt,
							fx_dt);


FREE_RETURN:

	if (sigma_date_dom)
	{
		free (sigma_date_dom);
	}

	if (sigma_dom)
	{
		free (sigma_dom);
	}

	if (tau_date_dom)
	{
		free (tau_date_dom);
	}

	if (tau_dom)
	{
		free (tau_dom);
	}

	if (sigma_date_for)
	{
		free (sigma_date_for);
	}

	if (sigma_for)
	{
		free (sigma_for);
	}

	if (tau_date_for)
	{
		free (tau_date_for);
	}

	if (tau_for)
	{
		free (tau_for);
	}

	if (sigma_date_fx)
	{
		free (sigma_date_fx);
	}

	if (sigma_fx)
	{
		free (sigma_fx);
	}

	if (sig_dom)
	{
		free (sig_dom);
	}

	if (sig_for)
	{
		free (sig_for);
	}

	if (sig_fx)
	{
		free (sig_fx);
	}

	if (merge_dates)
	{
		free (merge_dates);
	}

	return err;
}

Err Fx3DQuadtsCalibration(	
							long	today,
							double	*exercise_opt,
							double	*maturity_opt,
							double	*vol_opt,
							long	nbrOpt,
							long	nbrLong,
							double	*maturity,
							long	nbrMat,
							double	*sig_curve_dom,
							double	lda_dom,
							double	*sig_curve_for,
							double	lda_for,
							double	(*vol_ln_func)(	double	t,
													double	Spot,
													double	Fwd,
													double *params),
							double	*param,
							double	spot_fx,
							double	corr_dom_for,
							double	corr_dom_fx,
							double	corr_for_fx,
							char	*dom_yc,
							char	*for_yc,												
							double	**fx_vol_curve,			
							long	nbSteps,
							long	nbNewton,
							double	disc_dt,
							double	fx_dt,
							long	nbIterMax)

{
long	i, j, k;
long	nbrShort;
long	maturity_date, exercise_date;
double	Fwd;
double	strikes[1];
double	df;
double  price_tgt, price1[2], price2[2], error;
double	shift;



Err		err = NULL;

	nbrShort = nbrOpt - nbrLong;

	err = Fx3DLocaltsCalibration(	
								today,
								exercise_opt,
								maturity_opt,
								vol_opt,
								nbrOpt,
								maturity,
								nbrMat,
								sig_curve_dom,
								lda_dom,
								sig_curve_for,
								lda_for,
								vol_ln_func,
								param,
								spot_fx,
								corr_dom_for,
								corr_dom_fx,
								corr_for_fx,
								dom_yc,
								for_yc,
								fx_vol_curve,
								disc_dt,
								fx_dt,
								nbIterMax);

	
	if (err)
	{
		goto FREE_RETURN;
	}

	for (i=nbrShort; i<nbrOpt; i++)
	{
		maturity_date = (long) (today + maturity_opt[i] * 365.0000000001);
		exercise_date = (long) (today + exercise_opt[i] * 365.0000000001);		
		df = swp_f_df(today, exercise_date, dom_yc);
		Fwd = spot_fx * swp_f_df(today, exercise_date, for_yc) / df;
		err = OptBlkSch(Fwd, Fwd, vol_opt[i], exercise_opt[i], df, "CALL", "PREMIUM", &price_tgt);

		/* prices in the tree */
		strikes[0] = Fwd;
	
		error = 1.0E10;
		j = 0;

		while (j<nbNewton && (fabs(error) / price_tgt > MAX_ERROR2))
		{
			err = Fx3DQuadtsTreeFxOptions(	
											today,
											exercise_date,
											strikes,
											1,
											maturity,
											nbrMat,
											sig_curve_dom,
											lda_dom,
											sig_curve_for,
											lda_for,
											exercise_opt,
											nbrOpt,
											*fx_vol_curve,
											param[0],
											param[1],
											param[2],
											vol_opt[i],
											spot_fx,
											corr_dom_for,
											corr_dom_fx,
											corr_for_fx,
											dom_yc,
											for_yc,
											vol_opt,
											exercise_opt,
											nbrOpt,
											&(price1[0]),
											nbSteps
											);
			
			error = price_tgt - price1[0];

			if (fabs(error) / price_tgt > MAX_ERROR2)
			{
				/* calculates the first derivative */
				shift = SHIFT_VOL2 * (*fx_vol_curve)[i];
				if (error < 0)
				{
					shift *= -1.0;
				}
				

				(*fx_vol_curve)[i] += shift;
				err = Fx3DQuadtsTreeFxOptions(	
											today,
											exercise_date,
											strikes,
											1,
											maturity,
											nbrMat,
											sig_curve_dom,
											lda_dom,
											sig_curve_for,
											lda_for,
											exercise_opt,
											nbrOpt,
											*fx_vol_curve,
											param[0],
											param[1],
											param[2],
											vol_opt[i],
											spot_fx,
											corr_dom_for,
											corr_dom_fx,
											corr_for_fx,
											dom_yc,
											for_yc,
											vol_opt,
											exercise_opt,
											nbrOpt,												
											&(price2[0]),
											nbSteps
											);
				
				(*fx_vol_curve)[i] -= shift;
				shift *= error / (price2[0] - price1[0]);

				if ((error * (price2[0] - price1[0])) < 0)
				{
					/* nnot good !!!, add more steps */
					shift = 0.0;
					nbSteps = (long) (1.2 * nbSteps);
				}

				for (k=i; k<nbrOpt; k++)
				{
					(*fx_vol_curve)[k] += shift ;
				}
			}
			j += 1;
		}
	}


FREE_RETURN:

	
	return err;
}

Err Fx3DQuadCalibration(
						char	*dom_underlying,
						char	*for_underlying,		
						double	spot_fx,
						double	beta,
						double	gamma,
						double	correl_dom_for,
						double	correl_dom_fx,
						double	correl_for_fx,
						double	*exercise_opt,
						double	*maturity_opt,
						double	*vol_opt,
						long	nbropt,
						long	nbLong,
						double	**fx_vol_curve,
						long	nbSteps,
						long	nbNewton,
						double	disc_dt,
						double	fx_dt,
						long	nbIterMax)

{
long		sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
long		nb_merge_dates;
double		*sigma_date_dom = NULL,
			*sigma_dom		= NULL,
			*tau_date_dom	= NULL,
			*tau_dom		= NULL,
			*sigma_date_for	= NULL,
			*sigma_for		= NULL,
			*tau_date_for	= NULL,
			*tau_for		= NULL,
			*merge_dates	= NULL,
			*sig_dom		= NULL,
			*sig_for		= NULL;

double		lda_dom, lda_for;
double		Param[3];

SrtUndPtr	dom_und, for_und;

long		today, spot_date;
char		*dom_yc, *for_yc;


Err			err = NULL;

	err = Get_LGM_TermStructure (dom_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);	
	if (err)
	{
		goto FREE_RETURN;
	}
	
	err = Get_LGM_TermStructure (for_underlying, &sigma_date_for, &sigma_for, &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);	
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (tau_dom, tau_n_dom, &lda_dom);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (tau_for, tau_n_for, &lda_for);
	if (err)
	{
		goto FREE_RETURN;
	}

	dom_und =lookup_und(dom_underlying);
	for_und =lookup_und(for_underlying);
	
	err = get_underlying_discname(dom_und, &dom_yc);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = get_underlying_discname(for_und, &for_yc);
	if (err)
	{
		goto FREE_RETURN;
	}

	today = get_today_from_underlying (dom_und);

	spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	spot_fx = spot_fx / swp_f_df (today, spot_date, for_yc) * swp_f_df (today, spot_date, dom_yc);

	err = merge_rates_ts (
		sigma_date_dom, 
		sigma_dom,
		sigma_n_dom,
		sigma_date_for, 
		sigma_for,
		sigma_n_for,
		&merge_dates,
		&sig_dom,
		&sig_for,
		&nb_merge_dates);

	if (err)
	{
		goto FREE_RETURN;
	}	

	Param[0] = 1.0 - beta + gamma;
	Param[1] = beta - 2.0 * gamma;
	Param[2] = gamma;

	err = Fx3DQuadtsCalibration(	
								today,
								exercise_opt,
								maturity_opt,
								vol_opt,
								nbropt,
								nbLong,
								merge_dates,
								nb_merge_dates,
								sigma_dom,
								lda_dom,
								sigma_for,
								lda_for,
								vol_lnATM_func_Quad,
								&(Param[0]),								
								spot_fx,
								correl_dom_for,
								correl_dom_fx,
								correl_for_fx,
								dom_yc,
								for_yc,												
								fx_vol_curve,
								nbSteps,
								nbNewton,
								disc_dt,
								fx_dt,
								nbIterMax);



FREE_RETURN:

	if (sigma_date_dom)
	{
		free (sigma_date_dom);
	}

	if (sigma_dom)
	{
		free (sigma_dom);
	}

	if (tau_date_dom)
	{
		free (tau_date_dom);
	}

	if (tau_dom)
	{
		free (tau_dom);
	}

	if (sigma_date_for)
	{
		free (sigma_date_for);
	}

	if (sigma_for)
	{
		free (sigma_for);
	}

	if (tau_date_for)
	{
		free (tau_date_for);
	}

	if (tau_for)
	{
		free (tau_for);
	}

	if (sig_dom)
	{
		free (sig_dom);
	}

	if (sig_for)
	{
		free (sig_for);
	}

	if (merge_dates)
	{
		free (merge_dates);
	}

	return err;
}

Err Fx3DQuadtsTreeFxOptions(	
							long	today,
							long	maturity_date,
							double	*strikes,
							long	nbrOpt,
							double	*maturity,
							long	nbrMat,
							double	*sig_curve_dom,
							double	dom_lam,
							double	*sig_curve_for,
							double	for_lam,
							double	*maturity_fx,
							long	nbrMat_fx,
							double	*sig_curve_fx,
							double  alpha,
							double	beta,
							double	gamma,
							double	sig0,
							double	spot_fx,
							double	corr_dom_for,
							double	corr_dom_fx,
							double	corr_for_fx,
							char	*dom_yc,
							char	*for_yc,	
							double	*vols,
							double	*vols_time,
							int		nbVols,
							double	*option_prices,
							long	num_stp
							)
{
double	maturity_opt;
double	*time		= NULL,
		*date		= NULL;
long	nstp = 1;
long	i, j;

int		*vol_change = NULL;
double	*dom_vol	= NULL,
		*for_vol	= NULL,
		*fx_vol		= NULL;

double	*dom_ifr	= NULL,		
		*dom_fwd	= NULL,
		*dom_var	= NULL,
		*for_ifr	= NULL,
		*for_fwd	= NULL,
		*for_var	= NULL,
		*fx_fwd		= NULL,
		*fx_fwd0	= NULL,
		*fx_var		= NULL;

double	*sig_fx_approx	= NULL,
		*fx_fwd_approx	= NULL;

double	**func_parm		= NULL;
int		*has_evt		= NULL;

int		*is_bar = NULL;
double	*bar_lvl = NULL;
int		*bar_cl = NULL;

double	temp;

double	param[3];

LOCAL_MODEL_FUNC	*model_func = NULL;
LOCAL_MODEL_PARAM	*model_param = NULL;

Err		err = NULL;

	maturity_opt = (maturity_date -  today) / DAYS_IN_YEAR;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in Fx3DBetatsTreeFxOptions";
		goto FREE_RETURN;
	}
	
	time[0] = maturity_opt;

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&nstp, 
								0,
								NULL,
								0, 
								NULL, 
								num_stp);
	
	if (err)
	{
		goto FREE_RETURN;
	}

	date = (double*) calloc (nstp, sizeof (double));
	has_evt = (int*) calloc (nstp, sizeof (int));
	func_parm = dmatrix(0, nstp-1, 0, nbrOpt-1);

	if (!date  || !func_parm || !has_evt)
	{
		err = "Memory allocation error (2) in Fx3DBetatsTreeFxOptions";
		goto FREE_RETURN;
	}
	
	for (j=0; j<nbrOpt; j++)
	{
		func_parm[nstp-1][j] = strikes[j];
	}
	

	for (i=0; i<nstp; i++)
	{		
		date[i] = today + DAYS_IN_YEAR * time[i];
		has_evt[i] = 0;

		if (i > 0 && date[i] - date[i-1] >= 1)
		{
			date[i] = (long) (date[i] + 1.0e-08);
			time[i] = YEARS_IN_DAY * (date[i] - today);	
		}
	}
	has_evt[nstp-1] = 1;
	
	/*	Fill the model parameters as required by tree_main_3dfx */

	vol_change = (int*) calloc (nstp, sizeof (int));
	dom_vol = (double*) calloc (nstp, sizeof (double));
	for_vol = (double*) calloc (nstp, sizeof (double));
	fx_vol = (double*) calloc (nstp, sizeof (double));	

	is_bar = (int*) calloc (nstp, sizeof (int));
	bar_lvl = (double*) calloc (nstp, sizeof (double));
	bar_cl = (int*) calloc (nstp, sizeof (int));

	if (!vol_change || !dom_vol || !for_vol || !fx_vol || !is_bar || !bar_lvl || !bar_cl)
	{
		err = "Memory allocation error (3) in Fx3DQuadtsTreeFxOptions";
		goto FREE_RETURN;
	}

	vol_change[nstp-1] = 1;

	dom_vol[nstp-1] = sig_curve_dom[Get_Index(time[nstp-1], maturity, nbrMat)];
	for_vol[nstp-1] = sig_curve_for[Get_Index(time[nstp-1], maturity, nbrMat)];
	fx_vol[nstp-1] = sig_curve_fx[Get_Index(time[nstp-1], maturity_fx, nbrMat_fx)];	

	for (i=nstp-2; i>=0; i--)
	{
		dom_vol[i] = sig_curve_dom[Get_Index(time[i], maturity, nbrMat)];
		for_vol[i] = sig_curve_for[Get_Index(time[i], maturity, nbrMat)];
		fx_vol[i] = sig_curve_fx[Get_Index(time[i], maturity_fx, nbrMat_fx)];

		if (fabs (dom_vol[i] - dom_vol[i+1]) 
			+ fabs (for_vol[i] - for_vol[i+1]) 
			+ fabs (fx_vol[i] - fx_vol[i+1]) > EPS)
		{
			vol_change[i] = 1;
		}	
		else
		{
			vol_change[i] = 0;
		}
	}

	/*	Get distributions */
	dom_ifr = (double*) calloc (nstp, sizeof (double));
	dom_fwd = (double*) calloc (nstp, sizeof (double));
	dom_var = (double*) calloc (nstp, sizeof (double));
	for_ifr = (double*) calloc (nstp, sizeof (double));
	for_fwd = (double*) calloc (nstp, sizeof (double));
	for_var = (double*) calloc (nstp, sizeof (double));
	fx_fwd = (double*) calloc (nstp, sizeof (double));
	fx_fwd0 = (double*) calloc (nstp, sizeof (double));
	fx_var = (double*) calloc (nstp, sizeof (double));

	sig_fx_approx = (double*) calloc (nstp, sizeof (double));
	fx_fwd_approx = (double*) calloc (nstp, sizeof (double));

	if (!dom_ifr || !dom_fwd || !dom_var 
		|| !for_ifr || !for_fwd || !for_var
		|| !fx_fwd || !fx_fwd0 || !fx_var || !sig_fx_approx || !fx_fwd_approx)
	{
		err = "Memory allocation error (3) in Fx3DQuadtsTreeFxOptions";
		goto FREE_RETURN;
	}

	/* Compute all the fwd */
	fx_fwd0[0] = spot_fx;
	for (i=1; i<nstp; i++)
	{
		fx_fwd0[i] = spot_fx 
						/ swp_f_df (today, date[i], dom_yc)
						* swp_f_df (today, date[i], for_yc);		
	}



	param[0] = alpha;
	param[1] = beta;
	param[2] = gamma;

	/* computation of the model parameters */
	err = allocate_local_model(nstp, &model_func, &model_param);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = fill_local_model(model_func, model_param, nstp, time, alpha, beta, gamma, sig0, fx_fwd0);
	if (err)
	{
		goto FREE_RETURN;
	}

	/* first get the coresponding lognormal volatilities */
	err = FxLocal_log_approx(
							today,
							time,
							nstp,							
							dom_vol,		
							for_vol, 
							time,
							nstp,
							fx_vol,
							vol_lnATM_func_Quad,
							param,
							dom_lam, 
							for_lam, 
							corr_dom_for,
							corr_dom_fx,
							corr_for_fx,
							spot_fx,
							dom_yc,
							for_yc,	
							time,
							nstp,
							fx_fwd_approx,
							sig_fx_approx,
							fx_fwd0,
							MAX_TIME);


	fill_fwd_var (	nstp,
					time,
					date,
					dom_vol,		
					for_vol, 
					sig_fx_approx,
					dom_lam, 
					for_lam, 
					corr_dom_for,
					corr_dom_fx,
					corr_for_fx,
					dom_yc,
					for_yc,
					dom_ifr,
					dom_fwd,
					dom_var,
					for_ifr,
					for_fwd,
					for_var,
					fx_fwd,
					fx_var);		

	if (nbVols > 0)
	{
		/* we do not use the calculation but the vols given in input */
		for (i=1; i<nstp; i++)
		{
			fx_var[i] = interp(vols_time, vols, nbVols, time[i], 5, &temp);
			fx_var[i] *= fx_var[i] * time[i];
		}
	}


	err = treeQuad_main_3dfx(
						nstp,
						time,
						date,
						vol_change,	
						dom_vol,	
						for_vol, 
						fx_vol,						
						dom_ifr,	
						dom_fwd,
						dom_var,
						for_ifr,
						for_fwd,
						for_var,
						fx_fwd,
						fx_fwd0,
						fx_var,
						func_parm, 
						has_evt,
						bar_lvl,
						bar_cl,
						is_bar,
						dom_lam, 
						for_lam, 
						corr_dom_for,
						corr_dom_fx,
						corr_for_fx,
						spot_fx,
						dom_yc,
						for_yc,
						model_func,
						model_param,
						FxCall_payoff_4_3dfxQuad_tree,
						nbrOpt,
						1,
						option_prices);

FREE_RETURN:

	if (time) free (time);
	if (date) free (date);
	if (vol_change) free (vol_change);
	if (dom_vol) free (dom_vol);
	if (for_vol) free (for_vol);
	if (fx_vol) free (fx_vol);
	if (dom_ifr) free (dom_ifr);
	if (dom_fwd) free (dom_fwd);
	if (dom_var) free (dom_var);
	if (for_ifr) free (for_ifr);
	if (for_fwd) free (for_fwd);
	if (for_var) free (for_var);
	if (fx_fwd) free (fx_fwd);
	if (fx_fwd0) free (fx_fwd0);
	if (fx_var) free (fx_var);

	if (fx_fwd_approx) free (fx_fwd_approx);
	if (sig_fx_approx) free (sig_fx_approx);

	if (bar_lvl) free (bar_lvl);
	if (is_bar) free (is_bar);
	if (bar_cl) free (bar_cl);

	free_local_model_param(model_func, model_param);

	if (func_parm)
	{
		free_dmatrix(func_parm, 0, nstp-1, 0, nbrOpt-1);
	}

	return err;
}
