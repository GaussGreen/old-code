/**********************************************************************
 *      Name: SrtGrfnMainLgm2F.c                                           * 
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: L.C.				                                      *
 *      Date: 03/01/01                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 18/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "math.h"
#include "srt_h_all.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "SrtAccess.h"
#include "BGMEval.h"
#include "LGMSVpde.h"
#include "LGMSVMC.h"
#include "LGMSVGrfn.h"
#include "LGMSVUtil.h"
#include "LgmSVClosedForm.h"

char *SrtGrfnLGMSVpde(
				  char		*underlying,				  
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,

				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  
				  LGMSVParam Params,

				  int		nstept,
				  int		nstepx,
				  int		nstepeps,
				  int		nstepphi,

				  int		*nb_prod,
				  double	**prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMLGMSV		grfn_prm;
	int					forback;
	int					flag = 0;
	long				nstp;

	double				next_d;
						
	double				*evt_tms	= NULL,						
						*time		= NULL,
						*ifr		= NULL,
						*date		= NULL;		

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;

	double				*dff		= NULL,
						*gam		= NULL,
						*gam2		= NULL;
	
	double				final_mat, new_tstar;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domestic_name,
						*yc;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	int					i, j;

	LGMSV_model			model_, *model = &model_;
	
	Err					err = NULL;
	
	/*	Initialisation of the LGMSV model */
	init_NULL_LGMSV_model(model);

	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	/* End of Day Fixing */
	if (is_end_of_day_fixing)
	{
		defParm.end_of_day_fixings = SRT_YES;
		defParm.end_of_day_flg = SRT_YES;
	}
	else
	{
		defParm.end_of_day_fixings = SRT_NO;
		defParm.end_of_day_flg = SRT_NO;
	}

	/* End of Day Payment */
	if (is_end_of_day_payment)
	{
		defParm.end_of_day_payment = SRT_YES;
	}
	else
	{
		defParm.end_of_day_payment = SRT_NO;
	}

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved */	
	err = FIRSTGetUndFromDeal(
								&xStr,
								&num_und,
								&und_ptr);
					
	if (err)
	{
		goto FREE_RETURN;
	}

	if (num_und != 1)
	{
		err = "Product should involve only one underlying";
		goto FREE_RETURN;
	}

	und = und_ptr[0];
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	domestic_name = und -> underl_name;

	if (strcmp(domestic_name, und -> underl_name))
	{
		err = "Tableau uses different underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	yc = (char *) get_ycname_from_irund (und);
	
	/* Get number of columns */
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get the maximum number of dfs required	*/
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	/*	Next, get the time steps */
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);

	if (err) goto FREE_RETURN;

	/* Get the model */
	err = Get_LGMSV_model(	underlying,
							model);

	if (err) goto FREE_RETURN;

	if (model->iOne2F != 1)
	{
		err = "PDE pricer not available for LGMSV 2 Factor or more";
		goto FREE_RETURN;
	}		

	/* change of TStar */
	new_tstar = (double) ((long) (0.5 * (evt_dts[num_evt-1] - today)) * YEARS_IN_DAY);
	
	if (fabs(new_tstar) < 1.0E-08)
	{
		new_tstar = model->dTStar;
	} 

	Convert_Tstar_model(model, new_tstar);
	Params.Tstar = model->dTStar;

	/* discretise in time			*/

	nstp = num_evt;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));
	
	final_mat = xStr.tms[nstp-1];
	i = 0;
	while (i < model->iNbPWTime && model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, model->dPWTime[i]);
		i++;
	}
	
	num_f_sort_vector (nstp, time);
	num_f_unique_vector (&nstp, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&nstp,
								0,
								NULL,
								0,
								NULL, 
								nstept);				
	if (err)
	{
		goto FREE_RETURN;
	}

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));
	ifr = (double *) calloc (nstp, sizeof (double));

	if (max_num_df > 0)
	{
		dff = dvector(0, max_num_df-1);
		gam = dvector(0, max_num_df-1);
		gam2 = dvector(0, max_num_df-1);
	}

	if (!void_prm || !is_event || !ifr || !date ||  ((!dff || !gam || !gam2) && (max_num_df > 0)))
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;

	next_d = evt_dts[j] + 1;

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;

		if (next_d > date[i])
		{
			ifr[i] = swp_f_zr (date[i], next_d, yc);
		}
		else
		{
			ifr[i] = swp_f_zr (date[i], date[i] + 1, yc);
		}

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc(sizeof(grfn_parm_lgm_sv));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
			grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
			grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

			grfn_prm->dff = dff;
			grfn_prm->gam = gam;
			grfn_prm->gam2 = gam2;
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc(sizeof(grfn_parm_lgm_sv));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
			grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
			grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

			grfn_prm->dff = dff;
			grfn_prm->gam = gam;
			grfn_prm->gam2 = gam2;			

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}
		next_d = date[i];
	}

	/*	Eventually! call to function */

	(*prod_val) = (double *) dvector (0, num_col-1);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}
	if ( eventdates[numeventdates-1] >= today ) 
	{
		err = lgmSV_adi_UtPieceWise(nstp,
							time,
									date,
									nstepphi,
									nstepx,
									nstepeps,
									model->dLambdaX,
									1.0,
									model->dPWTime,
									model->dSigma,
									model->iNbPWTime,
									model->dLambdaEps,
									model->dLvlEps,
									model->dAlpha,
									model->dRho,
									&Params,
									void_prm,
									is_event,
									ifr,
									yc,
									payoff_lgmsv_pde_UtPieceWise,
									num_col,
									num_col,
									*prod_val);
	}	
	if (err)
	{
		goto FREE_RETURN;
	}
	
	*nb_prod = num_col;
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1] += xStr.gd->pv_of_past;

FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										nstp,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMLGMSV) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (max_num_df > 0)
	{
		if (dff) free_dvector(dff, 0, max_num_df-1);
		if (gam) free_dvector(gam, 0, max_num_df-1);
		if (gam2) free_dvector(gam2, 0, max_num_df-1);
	}

	if (is_event) free (is_event);
	if (date) free (date);
	if (ifr) free (ifr);
	if (time) free (time);
	
	free_LGMSV_model(model);

	return err;
}

char *SrtGrfnLGMSVMC(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					int			*nUsedEventDates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,
					
					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					double		*fwd_iv,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	Err			err = NULL;
	SrtUndPtr	pUndPtr;
	SrtMdlDim   mdl_dim;

	// Get the underlying through its name and check it exists
	// Check on the underlying type 
	pUndPtr  = lookup_und(underlying);
	if (!pUndPtr ) 
	{
		return serror("Underlying %s is not defined", underlying);				
	}

	if (get_underlying_type (pUndPtr) != INTEREST_RATE_UND)
	{
		return serror("Underlying %s is not of type IR", underlying);
	}

	if (get_mdltype_from_irund (pUndPtr) != LGM_STOCH_VOL)
	{
		return serror("Underlying %s is not of type LGMSV", underlying);
	}

	mdl_dim = get_mdldim_from_irund(pUndPtr);

	if (mdl_dim == ONE_FAC)
	{
		err = SrtGrfnLGMSVMC_1F(underlying,
								numeventdates, 
								eventdates,
								nUsedEventDates,
								tableauRows,
								tableauCols,
								tableauStrings,
								tableauMask,
								auxWidth,
								auxLen,
								aux,
								is_end_of_day_fixing,
								is_end_of_day_payment,
								Params,
								do_optimisation,
								optimise,
								fwd_iv,
								params,
								resRows,
								nstept,
								numpaths,		
								nb_prod,
								prod_val);
	}
	else if (mdl_dim == TWO_FAC)
	{
		err = SrtGrfnLGMSVMC_2F(underlying,
								numeventdates, 
								eventdates,
								nUsedEventDates,
								tableauRows,
								tableauCols,
								tableauStrings,
								tableauMask,
								auxWidth,
								auxLen,
								aux,
								is_end_of_day_fixing,
								is_end_of_day_payment,
								Params,
								do_optimisation,
								optimise,
								fwd_iv,
								params,
								resRows,
								nstept,
								numpaths,		
								nb_prod,
								prod_val);
	}
	else
	{
		return serror("Underlying %s has more than 2 dimensions", underlying);
	}

	return err;
}

char *SrtGrfnLGMSVMC_1F(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					int			*nUsedEventDates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					double		*fwd_iv,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMLGMSV		grfn_prm;
	int					forback;
	int					flag = 0;
	long				nstp;
						
	double				*evt_tms	= NULL,						
						*time		= NULL,
						*date		= NULL,
						*sigma		= NULL,
						*alpha		= NULL,
						*rho		= NULL,
						*lameps		= NULL,
						*lvleps		= NULL,

						*dff_star	= NULL,
						*gam_star	= NULL,
						*gam2_star	= NULL;

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;	

	long				index;
	double				final_mat, new_tstar;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domestic_name,
						*yc;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	LGMSV_model			model_, *model = &model_;

	int					i, j, k;

	double				df, temp;
	
	Err					err = NULL;
	
	/*	Initialisation of the LGMSV model */
	init_NULL_LGMSV_model(model);

	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	/* End of Day Fixing */
	if (is_end_of_day_fixing)
	{
		defParm.end_of_day_fixings = SRT_YES;
		defParm.end_of_day_flg = SRT_YES;
	}
	else
	{
		defParm.end_of_day_fixings = SRT_NO;
		defParm.end_of_day_flg = SRT_NO;
	}

	/* End of Day Payment */
	if (is_end_of_day_payment)
	{
		defParm.end_of_day_payment = SRT_YES;
	}
	else
	{
		defParm.end_of_day_payment = SRT_NO;
	}

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							*tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved */	
	err = FIRSTGetUndFromDeal(
								&xStr,
								&num_und,
								&und_ptr);
					
	if (err)
	{
		goto FREE_RETURN;
	}

	if (num_und != 1)
	{
		err = "Product should involve only one underlying";
		goto FREE_RETURN;
	}

	und = und_ptr[0];
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	domestic_name = und -> underl_name;

	if (strcmp(domestic_name, und -> underl_name))
	{
		err = "Tableau uses different underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	yc = (char *) get_ycname_from_irund (und);
	
	/* Get number of columns */
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get the maximum number of dfs required	*/
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	/*	Next, get the time steps */
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	/* Get the model */
	err = Get_LGMSV_model(	underlying,
							model);

	if (err) goto FREE_RETURN;

	if (Params.UseNewTStarMC)
	{
		/* change of TStar */
		new_tstar = xStr.tms[num_evt-1];
		
		if (fabs(new_tstar) < 1.0E-08)
		{
			new_tstar = model->dTStar;
		} 

		Convert_Tstar_model(model, new_tstar);		
	}

	Params.Tstar = model->dTStar;
	
	/* discretise in time			*/

	nstp = num_evt;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	final_mat = xStr.tms[nstp-1];
	i = 0;
	while (i < model->iNbPWTime && model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, model->dPWTime[i]);
		i++;
	}
	
	num_f_sort_vector (nstp, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&nstp,
								0,
								NULL,
								0,
								NULL, 
								nstept);

	if (err) goto FREE_RETURN;	

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));

	sigma = (double *) calloc (nstp, sizeof (double));
	alpha = (double *) calloc (nstp, sizeof (double));
	rho = (double *) calloc (nstp, sizeof (double));	
	lameps = (double *) calloc (nstp, sizeof (double));
	lvleps = (double *) calloc (nstp, sizeof (double));

	dff_star = (double *) calloc(num_evt, sizeof(double));
	gam_star = (double *) calloc(num_evt, sizeof(double));
	gam2_star = (double *) calloc(num_evt, sizeof(double));
	
	if (!void_prm || !is_event || !date || !sigma ||
		!alpha || !rho || !lameps || !lvleps ||
		!dff_star || !gam_star || !gam2_star)

	{
		err = "Memory allocation failure in SrtGrfnMainLGMSV_1F";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;
	
	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;
		
		/* For diffusion */
		if (i < nstp-1)
		{
			index = Get_Index(time[i+1], model->dPWTime, model->iNbPWTime);

			sigma[i+1] = model->dSigma[index];
			alpha[i+1] = model->dAlpha[index];
			rho[i+1] = model->dRho[index];			
			lameps[i+1] = model->dLambdaEps[index];
			lvleps[i+1] = model->dLvlEps[index];
		}

		/* For reconstruction */
		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc(sizeof(grfn_parm_lgm_sv));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
			grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
			grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];
			 
			/* DF(t, T*) reconstruction */
			dff_star[j] = swp_f_df (xStr.dts[j], model->lTStarDate, yc);
			dff_star[j] = log(dff_star[j]);
			gam_star[j] = (1.0 - exp(model->dLambdaX * (model->dTStar - evt_tms[j]))) / model->dLambdaX;
			gam2_star[j] = 0.5 * gam_star[j] * gam_star[j];

			if (grfn_prm->num_df > 0)
			{
				grfn_prm->dff = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam2 = dvector (0, grfn_prm->num_df-1);

				if (!grfn_prm->dff || !grfn_prm->gam || !grfn_prm->gam2)
				{
					err = "Memory allocation error in SrtGrfnMainLGMSV_1F";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_df; k++)
				{
					grfn_prm->dff[k] = swp_f_df (xStr.dts[j], grfn_prm->df_dts[k], yc);
					grfn_prm->dff[k] = log(grfn_prm->dff[k]);
					temp = (1.0 - exp(model->dLambdaX * (model->dTStar - evt_tms[j] - grfn_prm->df_tms[k]))) / model->dLambdaX;
					grfn_prm->gam[k] = temp - gam_star[j];
					grfn_prm->gam2[k] = 0.5 * temp * temp - gam2_star[j];
				}
			}

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc(sizeof(grfn_parm_lgm_sv));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
			grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
			grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}		
	}

	/*	Eventually! call to function */
	
	*nb_prod = num_col;	
	*tableauCols = num_col;	
	*resRows = max(num_col + 1, num_evt);
	*nUsedEventDates = num_evt;

	df = swp_f_df (today, model->lTStarDate, yc);

	if (do_optimisation)
	{
		if (params->iMultiIndex)
		{
			params->iNbIndex = params->iMultiIndex;
		}
		else
		{
			params->iNbIndex = 1;
		}

		mceb_allocate_params(	params,
								num_evt);

		if (params->iAdjustIV)
		{
			for (i=0; i<num_evt; i++)
			{
				params->dMarketFwdIV[i] = fwd_iv[i+numeventdates-xStr.num_evt] / df;
			}
		}

		*prod_val = dmatrix (0, *resRows - 1, 0, 2 + params->iNbIndex);
	}
	else
	{
		*prod_val = dmatrix (0, *nb_prod - 1, 0, 2);
	}

	if (!*prod_val)
	{
		err = "Memory allocation failure in SrtGrfnMainLGMSVMC";
		goto FREE_RETURN;
	}	
	
	switch (Params.UseReverseMC)
	{
		case 0:
		{
			err = lgmSV_mc_balsam(	nstp,
									num_evt,
									time,
									date,
									numpaths,
									model->dLambdaX,
									sigma,
									alpha,
									lameps,
									lvleps,
									rho,
									dff_star,
									gam_star,
									gam2_star,
									Params,
									void_prm,
									is_event,
									do_optimisation,
									optimise,
									params,
									NULL,
									payoff_lgmsv_mc,
									num_col,
									*prod_val);
			break;
		}

		case 1:
		{
			err = lgmSV_mc_balsam_rev(	nstp,
										num_evt,
										time,
										date,
										numpaths,
										model->dLambdaX,
										sigma,
										alpha,
										lameps,
										lvleps,
										rho,
										dff_star,
										gam_star,
										gam2_star,
										&Params,
										void_prm,
										is_event,
										do_optimisation,
										optimise,
										params,
										NULL,
										payoff_lgmsv_mc,
										NULL,
										num_col,
										*prod_val);
		}

		case 2:
		{
			err = lgmSV_mc_balsam(	nstp,
									num_evt,
									time,
									date,
									numpaths,
									model->dLambdaX,
									sigma,
									alpha,
									lameps,
									lvleps,
									rho,
									dff_star,
									gam_star,
									gam2_star,
									Params,
									void_prm,
									is_event,
									do_optimisation,
									optimise,
									params,
									NULL,
									payoff_lgmsv_mc,
									num_col,
									*prod_val);
			break;
		}
	}		

	if (err) goto FREE_RETURN;	
	
	if (do_optimisation)
	{
		/* Recopy Barrier / CoefLin for the moment */
		for (i=0; i<num_evt; i++)
		{
			(*prod_val)[i][2] = params->dBarrier[i];

			for (k=0; k<params->iNbIndex; k++)
			{
				(*prod_val)[i][3+k] = params->dCoefLin[i][k+1];
			}
		}		
	}

	for (i=0; i<num_col; i++)
	{
		(*prod_val)[i][0] *= df;
		(*prod_val)[i][1] *= df;
	}

	if (do_optimisation)
	{
		(*prod_val)[num_col][0] *= df;
		(*prod_val)[num_col][1] *= df;

		for (i=0; i<num_evt; i++)
		{
			if (params->iCalcOneTime) params->dOneTimeCall[i] *= df;
			if (params->iCalcOneTimePartial) params->dOneTimePartial[i] *= df;
			if (params->iCalcIV || params->iAdjustIV) params->dModelFwdIV[i] *= df;
		}
	}

	if (flag)
	{		
		for (i=0; i<num_evt-1; i++)
		{			
			(*prod_val)[i][2] = (*prod_val)[i+1][2];			
			
			if (do_optimisation)
			{
				for (k=0; k<params->iNbIndex; k++)
				{
					(*prod_val)[i][3+k] = (*prod_val)[i+1][3+k];
				}				
			}
		}

		if (do_optimisation)
		{
			mceb_shift_extrainfos(params);
		}
	}
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;

	if (do_optimisation)
	{
		(*prod_val)[num_col][0] += xStr.gd->pv_of_past;
	}
		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										nstp,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMLGMSV) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}


	if (is_event) free (is_event);
	if (date) free (date);
	if (time) free (time);
	
	if (sigma) free (sigma);
	if (alpha) free (alpha);	
	if (lameps) free (lameps);
	if (lvleps) free (lvleps);
	if (rho) free (rho);

	if (dff_star) free(dff_star);
	if (gam_star) free(gam_star);
	if (gam2_star) free(gam2_star);

	free_LGMSV_model(model);

	return err;
}

char *SrtGrfnLGMSVMC_2F(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					int			*nUsedEventDates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					double		*fwd_iv,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMLGMSV2F		grfn_prm;
	int					forback;
	int					flag = 0;
	long				nstp;	
						
	double				*evt_tms		= NULL,						
						*time			= NULL,
						*date			= NULL,

						*sigma			= NULL,
						*alpha			= NULL,
						*rho			= NULL,
						*rho2			= NULL,						
						*lameps			= NULL,
						*lvleps			= NULL,
						*newalphaLGM	= NULL,
						*newrhoLGM		= NULL,
						
						*dff_star		= NULL,
						*gam1_star		= NULL,
						*gam2_star		= NULL,
						*gam1_2_star	= NULL,
						*gam2_2_star	= NULL,
						*gam12_star		= NULL;

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;
		
	long				index;
	double				final_mat;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domestic_name,
						*yc;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	LGMSV_model			model_, *model = &model_;

	int					i, j, k;

	double				df;
	double				temp1, temp2;
	
	Err					err = NULL;
	
	/*	Initialisation of the LGMSV model */
	init_NULL_LGMSV_model(model);

	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	/* End of Day Fixing */
	if (is_end_of_day_fixing)
	{
		defParm.end_of_day_fixings = SRT_YES;
		defParm.end_of_day_flg = SRT_YES;
	}
	else
	{
		defParm.end_of_day_fixings = SRT_NO;
		defParm.end_of_day_flg = SRT_NO;
	}

	/* End of Day Payment */
	if (is_end_of_day_payment)
	{
		defParm.end_of_day_payment = SRT_YES;
	}
	else
	{
		defParm.end_of_day_payment = SRT_NO;
	}

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							*tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved */	
	err = FIRSTGetUndFromDeal(
								&xStr,
								&num_und,
								&und_ptr);
					
	if (err)
	{
		goto FREE_RETURN;
	}

	if (num_und != 1)
	{
		err = "Product should involve only one underlying";
		goto FREE_RETURN;
	}

	und = und_ptr[0];
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	domestic_name = und -> underl_name;

	if (strcmp(domestic_name, und -> underl_name))
	{
		err = "Tableau uses different underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	yc = (char *) get_ycname_from_irund (und);
	
	/* Get number of columns */
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get the maximum number of dfs required	*/
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	/*	Next, get the time steps */
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	/* Get the model */
	err = Get_LGMSV_model(	underlying,
							model);

	if (err) goto FREE_RETURN;

	if (Params.UseNewTStarMC)
	{
		/* change of TStar */

		/* NOT READY !!! */
		/*
		new_tstar = xStr.tms[num_evt-1];
		
		if (fabs(new_tstar) < 1.0E-08)
		{
			new_tstar = model->dTStar;
		} 

		Convert_Tstar_model(model, new_tstar);
		*/
	}

	Params.Tstar = model->dTStar;
																					
	/* discretise in time			*/

	nstp = num_evt;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	final_mat = xStr.tms[nstp-1];
	i = 0;
	while (i < model->iNbPWTime && model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, model->dPWTime[i]);
		i++;
	}
	
	num_f_sort_vector (nstp, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&nstp,
								0,
								NULL,
								0,
								NULL, 
								nstept);				
	if (err)
	{
		goto FREE_RETURN;
	}

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));

	sigma = (double *) calloc (nstp, sizeof (double));	
	alpha = (double *) calloc (nstp, sizeof (double));
	rho = (double *) calloc (nstp, sizeof (double));	
	rho2 = (double *) calloc (nstp, sizeof (double));
	lameps = (double *) calloc (nstp, sizeof (double));
	lvleps = (double *) calloc (nstp, sizeof (double));
	newalphaLGM = (double *) calloc (nstp, sizeof (double));
	newrhoLGM = (double *) calloc (nstp, sizeof (double));
	
	dff_star = (double *) calloc(num_evt, sizeof(double));
	gam1_star = (double *) calloc(num_evt, sizeof(double));
	gam2_star = (double *) calloc(num_evt, sizeof(double));
	gam1_2_star = (double *) calloc(num_evt, sizeof(double));
	gam2_2_star = (double *) calloc(num_evt, sizeof(double));
	gam12_star = (double *) calloc(num_evt, sizeof(double));
	
	if (!void_prm || !is_event || !date ||
		!sigma || !alpha || !rho || !rho2 || !lameps || !lvleps || !newalphaLGM || !newrhoLGM ||
		!dff_star || !gam1_star || !gam2_star || !gam1_2_star || !gam2_2_star || !gam12_star)
	{
		err = "Memory allocation failure in LGMSV2F MC";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;	
		
		if (i < nstp-1)
		{
			index = Get_Index(time[i+1], model->dPWTime, model->iNbPWTime);

			sigma[i+1] = model->dSigma[index];
			newalphaLGM[i+1] = model->dLGMAlpha[index];
			newrhoLGM[i+1] = model->dLGMRho[index];									
			
			alpha[i+1] = model->dAlpha[index];
			rho[i+1] = model->dRho[index];
			rho2[i+1] = model->dRho2[index];
			lameps[i+1] = model->dLambdaEps[index];
			lvleps[i+1] = model->dLvlEps[index];
		}

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc(sizeof(grfn_parm_lgm_sv_2F));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
			grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
			grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

			/* DF(t, T*) reconstruction */
			dff_star[j] = swp_f_df (xStr.dts[j], model->lTStarDate, yc);
			dff_star[j] = log(dff_star[j]);
			gam1_star[j] = (1.0 - exp(model->dLambdaX * (model->dTStar - evt_tms[j]))) / model->dLambdaX;			
			gam2_star[j] = (1.0 - exp(model->dLambdaX2 * (model->dTStar - evt_tms[j]))) / model->dLambdaX2;
			gam12_star[j] = gam1_star[j] * gam2_star[j];
			gam1_2_star[j] = 0.5 * gam1_star[j] * gam1_star[j];
			gam2_2_star[j] = 0.5 * gam2_star[j] * gam2_star[j];

			if (grfn_prm->num_df > 0)
			{
				grfn_prm->dff = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam1 = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam1_2 = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam2 = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam2_2 = dvector (0, grfn_prm->num_df-1);
				grfn_prm->gam12 = dvector (0, grfn_prm->num_df-1);				

				if (!grfn_prm->dff || !grfn_prm->gam1 || !grfn_prm->gam1_2 || !grfn_prm->gam2 || !grfn_prm->gam2_2 || !grfn_prm->gam12)
				{
					err = "Memory allocation error in SrtGrfnMainLGMSV_2F";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_df; k++)
				{
					grfn_prm->dff[k] = swp_f_df (xStr.dts[j], grfn_prm->df_dts[k], yc);
					grfn_prm->dff[k] = log(grfn_prm->dff[k]);
					temp1 = (1.0 - exp(model->dLambdaX * (model->dTStar - evt_tms[j] - grfn_prm->df_tms[k]))) / model->dLambdaX;
					temp2 = (1.0 - exp(model->dLambdaX2 * (model->dTStar - evt_tms[j] - grfn_prm->df_tms[k]))) / model->dLambdaX2;
					grfn_prm->gam1[k] = temp1 - gam1_star[j];
					grfn_prm->gam2[k] = temp2 - gam2_star[j];					
					grfn_prm->gam1_2[k] = 0.5 * temp1 * temp1 - gam1_2_star[j];
					grfn_prm->gam2_2[k] = 0.5 * temp2 * temp2 - gam2_2_star[j];
					grfn_prm->gam12[k] = temp1 * temp2 - gam12_star[j];
				}
			}
						
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc(sizeof(grfn_parm_lgm_sv));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
			grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
			grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}
	}

	/*	Eventually! call to function */
	
	*nb_prod = num_col;	
	*tableauCols = num_col;
	*resRows = max(num_col + 1, num_evt);
	*nUsedEventDates = num_evt;

	df = swp_f_df (today, model->lTStarDate, yc);

	if (do_optimisation)
	{
		if (params->iMultiIndex)
		{
			params->iNbIndex = params->iMultiIndex;
		}
		else
		{
			params->iNbIndex = 1;
		}		

		mceb_allocate_params(	params,
								num_evt);

		if (params->iAdjustIV)
		{
			for (i=0; i<num_evt; i++)
			{
				params->dMarketFwdIV[i] = fwd_iv[i+numeventdates-xStr.num_evt] / df;
			}
		}

		*prod_val = dmatrix (0, *resRows - 1, 0, 2 + params->iNbIndex);
	}
	else
	{
		*prod_val = dmatrix (0, *nb_prod - 1, 0, 2);
	}
	
/* Check to see if there are any event dates 
	if ( eventdates[numeventdates-1] >= today ) 
	{*/
		
		switch (Params.UseReverseMC)
		{
			case 0:
			{
				err = lgmSV2F_mc_balsam(nstp,
										num_evt,
										time,
										date,
										numpaths,
										model->dLambdaX,
										model->dLambdaX2,					
										sigma,
										newalphaLGM,
										newrhoLGM,
										alpha,
										lameps,
										lvleps,
										rho,
										rho2,
										dff_star,
										gam1_star,
										gam2_star,
										gam1_2_star,
										gam2_2_star,
										gam12_star,
										Params,
										void_prm,
										is_event,
										do_optimisation,
										&optimise[numeventdates-xStr.num_evt],
										params,
										NULL,
										payoff_lgmsv2F_mc,
										num_col,
										*prod_val);

				break;
			}

			case 1:
			{
				err = lgmSV2F_mc_balsam_rev(nstp,
											num_evt,
											time,
											date,
											numpaths,
											model->dLambdaX,
											model->dLambdaX2,					
											sigma,
											newalphaLGM,
											newrhoLGM,
											alpha,
											lameps,
											lvleps,
											rho,
											rho2,
											dff_star,
											gam1_star,
											gam2_star,
											gam1_2_star,
											gam2_2_star,
											gam12_star,
											&Params,
											void_prm,
											is_event,
											do_optimisation,
											&optimise[numeventdates-xStr.num_evt],
											params,
											NULL,
											payoff_lgmsv2F_mc,
											NULL,
											num_col,
											*prod_val);

				break;
			}

			case 2:
			{
				err = lgmSV2F_mc_balsam_optim_mem(	nstp,
													num_evt,
													time,
													date,
													numpaths,
													model->dLambdaX,
													model->dLambdaX2,					
													sigma,
													newalphaLGM,
													newrhoLGM,
													alpha,
													lameps,
													lvleps,
													rho,
													rho2,
													dff_star,
													gam1_star,
													gam2_star,
													gam1_2_star,
													gam2_2_star,
													gam12_star,
													&Params,
													void_prm,
													is_event,
													do_optimisation,
													&optimise[numeventdates-xStr.num_evt],
													params,
													NULL,
													payoff_lgmsv2F_mc,
													num_col,
													*prod_val);
				break;
			}				
		}

		if (err) goto FREE_RETURN;

		if (do_optimisation)
		{
			/* Recopy Barrier / CoefLin for the moment */
			for (i=0; i<num_evt; i++)
			{
				(*prod_val)[i][2] = params->dBarrier[i];

				for (k=0; k<params->iNbIndex; k++)
				{
					(*prod_val)[i][3+k] = params->dCoefLin[i][k+1];
				}
			}
		}

		for (i=0; i<num_col; i++)
		{
			(*prod_val)[i][0] *= df;
			(*prod_val)[i][1] *= df;
		}

		if (do_optimisation)
		{
			(*prod_val)[num_col][0] *= df;
			(*prod_val)[num_col][1] *= df;

			for (i=0; i<num_evt; i++)
			{
				if (params->iCalcOneTime) params->dOneTimeCall[i] *= df;
				if (params->iCalcOneTimePartial) params->dOneTimePartial[i] *= df;
				if (params->iCalcIV || params->iAdjustIV) params->dModelFwdIV[i] *= df;
			}
		}

		if (flag)
		{		
			for (i=0; i<num_evt-1; i++)
			{			
				(*prod_val)[i][2] = (*prod_val)[i+1][2];			
				
				if (do_optimisation)
				{
					for (k=0; k<params->iNbIndex; k++)
					{
						(*prod_val)[i][3+k] = (*prod_val)[i+1][3+k];
					}					
				}
			}

			if (do_optimisation)
			{
				mceb_shift_extrainfos(params);
			}
		}
//	}
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;
	if (do_optimisation)
		(*prod_val)[num_col][0] += xStr.gd->pv_of_past;

		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										nstp,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMLGMSV2F) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}	

	if (time) free (time);
	if (is_event) free (is_event);
	if (date) free (date);	

	if (sigma) free (sigma);
	if (alpha) free (alpha);
	if (rho) free (rho);
	if (rho2) free (rho2);
	if (lameps) free (lameps);
	if (lvleps) free (lvleps);
	if (newalphaLGM) free (newalphaLGM);
	if (newrhoLGM) free (newrhoLGM);

	if (dff_star) free(dff_star);
	if (gam1_star) free(gam1_star);
	if (gam2_star) free(gam2_star);
	if (gam1_2_star) free(gam1_2_star);
	if (gam2_2_star) free(gam2_2_star);
	if (gam12_star) free(gam12_star);

	free_LGMSV_model(model);
	
	return err;
}

char *SrtGrfnLGMSVFFT(
					/* Model */
					char		*underlying,

					/* Product */
					int			numeventdates, 
					long		*eventdates,
					long		tableauRows,
					long		tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,
					
					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,
				  
					/* Parameter of grids */
					int			iNbPhi,			/* Number of Phi : Should be a power of two */
					int			iNbft,			/* Number of ft : Should be a power of two */	
					double		iNbSigmaPhiGridLeft,
					double		iNbSigmaPhiGridRight,
					double		iNbSigmaftLeft,
					double		iNbSigmaftRight,
									
					double		iRatioPhi,
					double		iRatioFt,
					int			iPriorityFreqPhi,
					int			iPriorityFreqFt,

					/* outputs */
					int			*nb_prod,
					double		**prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMLGMSV		grfn_prm = NULL;
	int					forback;
	int					flag = 0;	
						
	double				*evt_tms	= NULL,						
						*time		= NULL,
						*SigTime	= NULL,
						*Sig		= NULL,
						*alphaTS	= NULL,
						*rhoTS		= NULL,
						*lamepsTS	= NULL,
						*ifr		= NULL,
						*date		= NULL;

	double				dAlpha, dRho, dLambdaEps;

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;

	double				*dff		= NULL,
						*gam		= NULL,
						*gam2		= NULL;

	long				iNbSigTime;
	double				tau;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domestic_name,
						*yc;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	int					i, endi;

	double				dLambdaX, dTStar;
	double				dExTime;
	long				dExDate;
	double				Coef1Re, Coef2Re, Coef2ImInit, Coef3ReInit, Coef3ImInit;
	double				t1, t2, sig1, AlphaEq;
	double				dPhitMean, dPhitStd, dPhiStep;
	double				dftStd, dftStep;
	double				dPhiFreqStep, dftFreqStep;
	int					iIndexPhiMean, iIndexft0;

	double				LimitPhi, Limitft;
	double				dB0Tstar;
	int					one2F;

	double				*dt			= NULL,
						*Coef2ImT	= NULL,
						*Coef3ReT	= NULL,
						*Coef3ImT	= NULL;

	double			lambdaArray[10];

	/* For moments calculation */
	LGMSVSolFunc	FuncPhi_, *FuncPhi = &FuncPhi_;
	LGMSVSolFunc	FuncPhi2_, *FuncPhi2 = &FuncPhi2_;
	LGMSVSolFunc	FuncV2_, *FuncV2 = &FuncV2_;
	LGMSVSolFunc	FuncPhiV_, *FuncPhiV = &FuncPhiV_;
	double			ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV;
	
	Err					err = NULL;
	
	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	/* End of Day Fixing */
	if (is_end_of_day_fixing)
	{
		defParm.end_of_day_fixings = SRT_YES;
		defParm.end_of_day_flg = SRT_YES;
	}
	else
	{
		defParm.end_of_day_fixings = SRT_NO;
		defParm.end_of_day_flg = SRT_NO;
	}

	/* End of Day Payment */
	if (is_end_of_day_payment)
	{
		defParm.end_of_day_payment = SRT_YES;
	}
	else
	{
		defParm.end_of_day_payment = SRT_NO;
	}

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	if (numeventdates > 1)
	{
		err = "Algo doesn't accepte more than one row in the GRFN tableau";
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved */	
	err = FIRSTGetUndFromDeal(
								&xStr,
								&num_und,
								&und_ptr);
					
	if (err)
	{
		goto FREE_RETURN;
	}

	if (num_und != 1)
	{
		err = "Product should involve only one underlying";
		goto FREE_RETURN;
	}

	und = und_ptr[0];
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	domestic_name = und -> underl_name;

	if (strcmp(domestic_name, und -> underl_name))
	{
		err = "Tableau uses different underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	yc = (char *) get_ycname_from_irund (und);
	
	/* Get number of columns */
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get the maximum number of dfs required	*/
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	/*	Next, get the time steps */
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = Get_LGMSV_TermStructure(	domestic_name,
									&SigTime,
									&Sig,
									&alphaTS,
									&rhoTS,
									&lamepsTS,
									&dTStar,
									&iNbSigTime,
									&tau,
									&one2F,
									NULL,
									NULL,
									NULL);

	if (err)
	{
		goto FREE_RETURN;
	}

	/* for now no TS */

	dAlpha = alphaTS[0];
	dRho = rhoTS[0],
	dLambdaEps = lamepsTS[0];

	for (i=1; i<iNbSigTime; i++)
	{
		if ((fabs(alphaTS[i] - dAlpha) > 1.0E-08) || (fabs(rhoTS[i] - dRho) > 1.0E-08) || (fabs(lamepsTS[i] - dLambdaEps) > 1.0E-08))
		{
			err = "no alpha, rho, lameps TS allowed for now";
			goto FREE_RETURN;
		}
	}

	dLambdaX = 1.0 / tau;	

	/* discretise in time			*/

	dExTime = evt_tms[0];
	dExDate = evt_dts[0]; 

	endi = Get_Index(dExTime, SigTime, iNbSigTime);

	dt = dvector(0, endi);
	Coef2ImT = dvector(0, endi);
	Coef3ReT = dvector(0, endi);
	Coef3ImT = dvector(0, endi);			
	
	if (max_num_df > 0)
	{
		dff = dvector(0, max_num_df-1);
		gam = dvector(0, max_num_df-1);
		gam2 = dvector(0, max_num_df-1);
	}

	if (!dt || !Coef2ImT || !Coef3ReT || !Coef3ImT || ((!dff || !gam || !gam2) && (max_num_df > 0)))
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	/* Constant independent on time */	
	Coef1Re	= -0.5 * dAlpha * dAlpha;
	Coef2Re = dLambdaEps;
	Coef2ImInit = -dAlpha * dRho;
	Coef3ReInit = 0.5;
	Coef3ImInit = exp(-2.0 * dLambdaX * (dExTime - dTStar));

	Limitft = log(iRatioFt);
	LimitPhi = log(iRatioPhi);

	LGMSVMomentInit2(0.0, dAlpha, dLambdaEps, dRho, FuncPhi, FuncPhi2, FuncV2, FuncPhiV);
	ExpectPhi = 0.0;
	ExpectPhi2 = 0.0;
	ExpevtV2 = 1.0;
	ExpectPhiV = 0.0;
	
	t1 = 0.0;
	for (i=0; i<=endi; i++)
	{
		/* Precalculation on the option i*/

		/* First the time discretisation */			
		if (i < endi)
		{
			t2 = SigTime[i];
		}
		else
		{
			t2 = dExTime;
		}
				
		AlphaEq = dAlpha / 2.0;
		if (dLambdaEps > 1.0E-16)
		{
			AlphaEq *= sqrt((1.0 - exp(-dLambdaEps * t2)) / (dLambdaEps * t2));
		}		

		/* Calculate constant values */
		sig1 = Sig[i];
		dt[i] = (t2 - t1);
		Coef2ImT[i] = Coef2ImInit * sig1;
		Coef3ReT[i] = Coef3ReInit * sig1 * sig1;
		Coef3ImT[i] = -Coef3ImInit * sig1 * sig1;		

		LGMSVMomentCalculation2(sig1 * sig1, t2 - t1, ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV,
			FuncPhi, FuncPhi2, FuncV2, FuncPhiV, lambdaArray, &ExpectPhi, &ExpectPhi2, &ExpevtV2, &ExpectPhiV);

		t1 = t2;
	}

	dPhitMean = ExpectPhi * Coef3ImInit;
	dPhitStd = sqrt(ExpectPhi2 - ExpectPhi * ExpectPhi) * Coef3ImInit;
	dftStd = sqrt(ExpectPhi);

	/* Find Frequence */		
	LGMSVFindFreqNew(iNbPhi, iNbft, AlphaEq, dRho, endi, dt,
					Coef1Re, Coef2Re, Coef2ImT, Coef3ReT, Coef3ImT,
					dPhitMean, dPhitStd, dftStd, iNbSigmaPhiGridLeft, iNbSigmaPhiGridRight, iNbSigmaftLeft, iNbSigmaftRight,
					LimitPhi, Limitft, iPriorityFreqPhi, iPriorityFreqFt,
					&dPhiFreqStep, &dPhiStep, &iIndexPhiMean, &dftFreqStep, &dftStep, &iIndexft0);


	
	grfn_prm = malloc(sizeof(grfn_parm_lgm_sv));
	grfn_prm->global = &xStr;
	grfn_prm->local = xStr.evt;

	grfn_prm->num_df = xStr.evt[0].evt->dflen[0];
	grfn_prm->df_tms = xStr.evt[0].evt->dft[0];
	grfn_prm->df_dts = xStr.evt[0].evt->dfd[0];

	grfn_prm->dff = dff;
	grfn_prm->gam = gam;
	grfn_prm->gam2 = gam2;
						
	void_prm = &((void*) grfn_prm);		

	/*	Eventually! call to function */

	*prod_val = dvector (0, num_col-1);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	err = LGMSVOptionGRFNPrice(	
								iNbPhi,
								iNbft,
								dLambdaX,
								dTStar,
								endi,
								dt,
								Coef1Re,
								Coef2Re,
								Coef2ImT,
								Coef3ReT,
								Coef3ImT,
								dPhitMean,									
								0,
								dExTime,
								dExDate,													
								dPhiFreqStep,
								dPhiStep,
								iIndexPhiMean,
								dftFreqStep,
								dftStep,
								iIndexft0,
								yc,
								void_prm, 
								payoff_lgmsv_FFT,
								num_col, 							
								*prod_val);
	
	if (err)
	{
		goto FREE_RETURN;
	}

	dB0Tstar  = swp_f_df (today, today + dTStar * 365.000000000000001, yc);

	for (i=0; i<num_col; i++)
	{
		(*prod_val)[i] *= dB0Tstar;
	}

	*nb_prod = num_col;
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1] += xStr.gd->pv_of_past;

		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										numeventdates,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (grfn_prm)
	{
		grfn_prm = (GRFNPARMLGMSV) (*void_prm);
		free (grfn_prm);
	}
			

	if (max_num_df > 0)
	{
		if (dff) free_dvector(dff, 0, max_num_df-1);
		if (gam) free_dvector(gam, 0, max_num_df-1);
		if (gam2) free_dvector(gam2, 0, max_num_df-1);
	}	

	if (SigTime) free (SigTime);
	if (Sig) free (Sig);
	if (alphaTS) free (alphaTS);
	if (rhoTS) free (rhoTS);
	if (lamepsTS) free (lamepsTS);

	return err;
}