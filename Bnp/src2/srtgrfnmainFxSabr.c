/**********************************************************************
 *      Name: SrtGrfnMainFxSabr.c                                     * 
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

#include "srt_h_all.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "SrtAccess.h"
#include "BGMEval.h"
#include "FxSabrAdi.h"
#include "FxSabrSLAdi.h"
#include "FxSabrQuadAdi.h"
#include "FxSabrGrfn.h"
#include "opfnctns.h"
#include "math.h"

char *SrtGrfnFxSabrAdi(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMFXSABR		grfn_prm;

	double				next_d;
	
	long				*evt_dts	= NULL;
	double				*evt_tms	= NULL;

	double				*time		= NULL,
						*sig		= NULL,
						*drift	= NULL,
						*date		= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;

	int					ns;

	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;

	int					fx_idx,
						dom_idx,
						for_idx;

	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL;

	long				sigma_n_fx, index;

	double				last_mat, coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	TermStruct			*fx_ts;

	int					*is_bar = NULL;
	double				*bar_lvl = NULL;
	int					*bar_cl = NULL;
	int					has_bar;

	double				next_bar;

	int					i, j, forback;
	
	Err					err = NULL;
	
	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

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

	if (num_und > 3)
	{
		err = "Product should involve only maximum 3 underlyings";
		goto FREE_RETURN;
	}

	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	fx_idx = dom_idx = for_idx = -1;
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);

		if (num_und > 1 || strcmp (und_ptr[0]->underl_name, underlying))
		{
			err = "Error, the GRFN tableau has too many underlyings involved";
			goto FREE_RETURN;
		}

		fx_idx = 0;
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);

		for (i=0; i<num_und; i++)
		{
			if (!strcmp (und_ptr[i]->underl_name, underlying))
			{
				fx_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, domname))
			{
				dom_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, forname))
			{
				for_idx = i;
			}
		}

		if (fx_idx == -1)
		{
			err = "The Fx underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (dom_idx == -1)
		{
			err = "The domestic underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (for_idx == -1)
		{
			err = "The foreign underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

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

	/* discretise in time			*/

	ns = num_evt;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, ns * sizeof (double));

	last_mat = xStr.tms[ns-1];

	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < last_mat)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	is_bar = (int*) calloc (nstp, sizeof (int));
	bar_lvl = (double*) calloc (nstp, sizeof (double));
	bar_cl = (int*) calloc (nstp, sizeof (int));

	if (!void_prm || !is_event || !sig || !date || !drift || !is_bar || !bar_lvl || !bar_cl)
	{
		err = "Memory allocation failure in SrtGrfnMainAdi2";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}
	
	j = xStr.num_evt - 1;
	next_bar = 0.0;
	next_d = evt_dts[j] + 1;
	has_bar = 0;

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		//sig[i] = find_fx_sig (time[i], fx_ts);

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx != -1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			if (barrier[j] > 1.0e-08)
			{
				/* mat_date = add_unit ((long) (date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);*/
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					has_bar = 1;
					bar_cl[i] = bar_col[j];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}				
			}
			else
			{
				is_bar[i] = 0;
			}			

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx !=-1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}

			if (barrier[j] > 1.0e-08)
			{
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					bar_cl[i] = bar_col[j];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}
			}
			else
			{
				is_bar[i] = 0;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			if (is_cont_bar[j+1])
			{
				is_bar[i] = 1;
				bar_lvl[i] = barrier[j+1] * swp_f_df(spot_date, date[i], dom_yc) / swp_f_df(spot_date, date[i], for_yc);
				bar_cl[i] = bar_col[j+1];
			}
			else
			{
				is_bar[i] = 0;
			}

			is_event[i] = 0;
			void_prm[i] = NULL;
		}
		next_d = date[i];
	}

	/* Has bar */
	if (!has_bar)
	{
		free (bar_lvl);
		bar_lvl = NULL;
	}
	
			
	/*	Eventually! call to function */

	*prod_val = dvector (0, num_col-1);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}
		
	err = FxSabr_adi(	nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						void_prm,
						is_event,
						bar_lvl,
						bar_cl,
						is_bar,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi,
						num_col,
						*prod_val);					

	if (err)
	{
		goto FREE_RETURN;
	}
	
	*nprod = num_col;
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1] += xStr.gd->pv_of_past;
		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(	num_und,
								&und_ptr);

		FIRSTFreeEvtDatesFromDeal(	nstp,
									&evt_dts,
									&evt_tms);

		FIRSTFreeMktStruct (&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMFXSABR) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);
	if (is_event) free (is_event);

	if (bar_lvl) free (bar_lvl);
	if (is_bar) free (is_bar);
	if (bar_cl) free (bar_cl);
		
	return err;
}

char *SrtGrfnFxSabrAdi2(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		*is_up_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMFXSABR		grfn_prm;

	double				next_d;
	
	long				*evt_dts	= NULL;
	double				*evt_tms	= NULL;

	double				*time		= NULL,
						*sig		= NULL,
						*drift	= NULL,
						*date		= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;

	int					ns;

	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;

	int					fx_idx,
						dom_idx,
						for_idx;

	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL;

	long				sigma_n_fx, index;

	double				last_mat, coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	TermStruct			*fx_ts;

	int					*is_bar = NULL,
						*is_up	= NULL;
	double				*bar_lvl = NULL;
	int					bar_cl;

	double				next_bar;

	int					i, j, forback;
	
	Err					err = NULL;
	
	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

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

	if (num_und > 3)
	{
		err = "Product should involve only maximum 3 underlyings";
		goto FREE_RETURN;
	}

	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	fx_idx = dom_idx = for_idx = -1;
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);

		if (num_und > 1 || strcmp (und_ptr[0]->underl_name, underlying))
		{
			err = "Error, the GRFN tableau has too many underlyings involved";
			goto FREE_RETURN;
		}

		fx_idx = 0;
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);

		for (i=0; i<num_und; i++)
		{
			if (!strcmp (und_ptr[i]->underl_name, underlying))
			{
				fx_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, domname))
			{
				dom_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, forname))
			{
				for_idx = i;
			}
		}

		if (fx_idx == -1)
		{
			err = "The Fx underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (dom_idx == -1)
		{
			err = "The domestic underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (for_idx == -1)
		{
			err = "The foreign underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

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

	/* discretise in time			*/

	ns = num_evt;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, ns * sizeof (double));

	last_mat = xStr.tms[ns-1];

	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < last_mat)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	is_bar = (int*) calloc (nstp, sizeof (int));
	is_up = (int*) calloc (nstp, sizeof (int));
	bar_lvl = (double*) calloc (nstp, sizeof (double));
	/*
	bar_cl = (int*) calloc (nstp, sizeof (int));
	*/

	if (!void_prm || !is_event || !sig || !date || !drift || !is_bar || !bar_lvl || !is_up)
	{
		err = "Memory allocation failure in SrtGrfnMainFXSAbrAdi";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}
	
	j = xStr.num_evt - 1;
	next_bar = 0.0;
	next_d = evt_dts[j] + 1;

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		//sig[i] = find_fx_sig (time[i], fx_ts);

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx != -1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			if (barrier[j] > 1.0e-08)
			{
				/* mat_date = add_unit ((long) (date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);*/
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					is_up[i] = is_up_bar[j];
					/*
					bar_cl[i] = bar_col[j];
					*/
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}
			}
			else
			{
				is_bar[i] = 0;
			}			

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>0 && xStr.am[j-1])
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j - 1;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j-1].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j-1].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j-1].evt->dfd[fx_idx];
			
			if (dom_idx !=-1)
			{
				grfn_prm->num_dom_df = xStr.evt[j-1].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j-1].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j-1].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j-1].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j-1].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j-1].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}

			if (barrier[j-1] > 1.0e-08)
			{
				next_bar = barrier[j-1] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j-1] >= 0 && bar_col[j-1] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					/*
					bar_cl[i] = bar_col[j-1];
					*/
					is_up[i] = is_up_bar[j-1];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}
			}
			else
			{
				is_bar[i] = 0;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			if (is_cont_bar[j+1])
			{
				is_bar[i] = 1;
				bar_lvl[i] = barrier[j+1] * swp_f_df(spot_date, date[i], dom_yc) / swp_f_df(spot_date, date[i], for_yc);
				/*
				bar_cl[i] = bar_col[j+1];
				*/
				is_up[i] = is_up_bar[j+1];
			}
			else
			{
				is_bar[i] = 0;
			}

			is_event[i] = 0;
			void_prm[i] = NULL;
		}
		next_d = date[i];
	}

	/* find the barrier column */
	
	j = xStr.num_evt - 1;
	
	bar_cl = bar_col[j];

	for (i=0; i<j-1 ; i++)
	{
		if (bar_col[i] != bar_cl)
		{
			err = "The Barrier column should be constant in SrtGrfnMainFXSAbrAdi";
		}
	}
	
			
	/*	Eventually! call to function */

	*prod_val = dvector (0, num_col-1);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}
		
	err = FxSabr_adi_bar(nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						void_prm,
						is_event,
						bar_lvl,
						bar_cl,
						is_bar,
						is_up,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi,
						num_col,
						*prod_val);					

	if (err)
	{
		goto FREE_RETURN;
	}
	
	*nprod = num_col;
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1] += xStr.gd->pv_of_past;
		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(	num_und,
								&und_ptr);

		FIRSTFreeEvtDatesFromDeal(	nstp,
									&evt_dts,
									&evt_tms);

		FIRSTFreeMktStruct (&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMFXSABR) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);
	if (is_event) free (is_event);

	if (bar_lvl) free (bar_lvl);
	if (is_bar) free (is_bar);
	/*if (bar_cl) free (bar_cl);*/
	if (is_up) free (is_up);
		
	return err;
}

Err FxSabrAdiKO(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  double	maturity,
   				  double	strike,
				  int		is_call,	/* 1 Call, 0: Put */
				  double	barrier,
				  int		is_up,		/* 1 Up, 0: Down */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res)
{	
	double				*time		= NULL,
						*sig		= NULL,
						*drift		= NULL,
						*date		= NULL;
	
	long				today, spot_date;
						
	int					num_col,
						ns;
	
	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;
	
	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL,
						*prod_val		= NULL;

	long				sigma_n_fx, index;

	double				coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	
	TermStruct			*fx_ts;
	
	double				*bar_lvl	= NULL,
						**func_parm	= NULL;

	int					*is_event	= NULL;

	double				premium;

	int					i;
	
	Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und);

	if ((is_up && spot_fx > barrier) || (!is_up && spot_fx < barrier))
	{
		*res = 0.0;
		goto FREE_RETURN;
	}
	
	spot_fx *= swp_f_df (today, spot_date, dom_yc)
				/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

	
	/* discretise in time			*/
	
	ns = 2;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = maturity;
	
	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	bar_lvl = (double*) calloc (nstp, sizeof (double));	

	if (!sig || !date || !drift || !bar_lvl)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		bar_lvl[i] = barrier * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);		
	}

	
			
	/*	Eventually! call to function */

	num_col = 1;

	prod_val = dvector (0, num_col-1);

	if (!prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}
		
	err = FxSabr_KOOption(	nstp,
							time,
							date,
							nstpfx,
							nstpvol,
							sigma_fx[0],
							drift,
							alpha,
							beta,
							rho,
							lambda,
							floormu,
							strike,
							is_call,
							bar_lvl,
							is_up,
							is_cvx,
							spot_fx,
							dom_yc,
							for_yc,							
							prod_val);					

	if (err)
	{
		goto FREE_RETURN;
	}

	if (!is_ko)
	{
		/* now we need to price the regular option */

		func_parm = dmatrix(0, nstp-1, 0, 0);
		is_event = calloc(nstp, sizeof(int));

		if (!func_parm || !is_event)
		{
			err = "Memory allocation faillure in FxSabrAdiKO";
			goto FREE_RETURN;
		}

		func_parm[nstp-1][0] = strike;
		is_event[nstp-1] = 1;

		err = FxSabr_adi(	nstp,
							time,
							date,
							nstpfx,
							nstpvol,
							sigma_fx[0],
							drift,
							alpha,
							beta,
							rho,
							lambda,
							floormu,
							func_parm,
							is_event,
							NULL,
							NULL,
							NULL,
							spot_fx,
							dom_yc,
							for_yc,
							payoff_fx_sabr_adi_opt,
							1,
							&premium
							);

		if (err)
		{
			goto FREE_RETURN;
		}

		*res = premium - prod_val[0];
	}
	else
	{
		*res = prod_val[0];
	}
		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);	

	if (bar_lvl) free (bar_lvl);
	if (prod_val) free_dvector (prod_val, 0, num_col-1);

	if (func_parm) free_dmatrix(func_parm, 0, nstp-1, 0, 0);
	if (is_event) free (is_event);
		
	return err;
}

Err FxSabrSmile(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  long		mat_date,
   				  double	*strike,
				  int		nb_strike,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	**res)
{	
	double				*time		= NULL,
						*sig		= NULL,
						*drift		= NULL,
						*date		= NULL;
	
	long				today, spot_date, exe_date;
						
	int					num_col,
						ns;
	
	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;
	
	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL,
						*prod_val		= NULL;

	long				sigma_n_fx, index;

	double				coef, fwd_fx, df, maturity, bs_vol;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	
	TermStruct			*fx_ts;
	
	double				**func_parm	= NULL;

	int					*is_event	= NULL;	

	int					i;
	
	Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) * swp_f_df (today, spot_date, dom_yc)
										/ swp_f_df (today, spot_date, for_yc);

	exe_date = add_unit (mat_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	maturity = (exe_date - today) * YEARS_IN_DAY;

	fwd_fx = get_spot_from_fxund (und) * swp_f_df (spot_date, mat_date, for_yc)
										/ swp_f_df (spot_date, mat_date, dom_yc);
	
	df = swp_f_df(today, exe_date, dom_yc);

	/* discretise in time			*/
	
	ns = 2;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = maturity;
	
	/*	Add revelant vol times	*/
	
	fx_ts = get_ts_from_fxund (und);

	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));
		
	if (!sig || !date || !drift)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
	}

	
			
	/*	Eventually! call to function */

	num_col = nb_strike;	
	
	/* now we need to price the regular option */

	func_parm = dmatrix(0, nstp-1, 0, num_col-1);	
	*res = dvector(0, nb_strike-1);
	is_event = calloc(nstp, sizeof(int));

	if (!func_parm || !is_event || !(*res))
	{
		err = "Memory allocation faillure in FxSabrAdiKO";
		goto FREE_RETURN;
	}

	for (i=0; i<num_col; i++)
	{
		func_parm[nstp-1][i] = strike[i];
	}
	
	is_event[nstp-1] = 1;

	err = FxSabr_adi(	nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						func_parm,
						is_event,
						NULL,
						NULL,
						NULL,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi_opt,
						num_col,
						*res
						);

	if (err)
	{
		goto FREE_RETURN;
	}

	for (i=0; i<nb_strike; i++)
	{
		err = srt_f_optimpvol(		
								(*res)[i],
								fwd_fx,
								strike[i],
								maturity,								
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								&bs_vol);

		(*res)[i] = bs_vol;
	}	

		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);	

	if (func_parm) free_dmatrix(func_parm, 0, nstp-1, 0, num_col-1);
	if (is_event) free (is_event);
		
	return err;
}

char *SrtGrfnFxSabrSLAdi(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMFXSABR		grfn_prm;

	double				next_d;
	
	long				*evt_dts	= NULL;
	double				*evt_tms	= NULL;

	double				*time		= NULL,
						*sig		= NULL,
						*drift	= NULL,
						*date		= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;

	int					ns;

	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;

	int					fx_idx,
						dom_idx,
						for_idx;

	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL;

	long				sigma_n_fx, index;

	double				last_mat, coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	TermStruct			*fx_ts;

	int					*is_bar = NULL;
	double				*bar_lvl = NULL;
	int					*bar_cl = NULL;
	int					has_bar;

	double				next_bar;

	int					i, j, forback;
	
	Err					err = NULL;
	
	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

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

	if (num_und > 3)
	{
		err = "Product should involve only maximum 3 underlyings";
		goto FREE_RETURN;
	}

	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	fx_idx = dom_idx = for_idx = -1;
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);

		if (num_und > 1 || strcmp (und_ptr[0]->underl_name, underlying))
		{
			err = "Error, the GRFN tableau has too many underlyings involved";
			goto FREE_RETURN;
		}

		fx_idx = 0;
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);

		for (i=0; i<num_und; i++)
		{
			if (!strcmp (und_ptr[i]->underl_name, underlying))
			{
				fx_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, domname))
			{
				dom_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, forname))
			{
				for_idx = i;
			}
		}

		if (fx_idx == -1)
		{
			err = "The Fx underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (dom_idx == -1)
		{
			err = "The domestic underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (for_idx == -1)
		{
			err = "The foreign underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

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

	/* discretise in time			*/

	ns = num_evt;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, ns * sizeof (double));

	last_mat = xStr.tms[ns-1];

	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < last_mat)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	is_bar = (int*) calloc (nstp, sizeof (int));
	bar_lvl = (double*) calloc (nstp, sizeof (double));
	bar_cl = (int*) calloc (nstp, sizeof (int));

	if (!void_prm || !is_event || !sig || !date || !drift || !is_bar || !bar_lvl || !bar_cl)
	{
		err = "Memory allocation failure in SrtGrfnMainAdi2";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}
	
	j = xStr.num_evt - 1;
	next_bar = 0.0;
	next_d = evt_dts[j] + 1;
	has_bar = 0;

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		//sig[i] = find_fx_sig (time[i], fx_ts);

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx != -1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			if (barrier[j] > 1.0e-08)
			{
				/* mat_date = add_unit ((long) (date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);*/
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					has_bar = 1;
					bar_cl[i] = bar_col[j];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}				
			}
			else
			{
				is_bar[i] = 0;
			}			

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx !=-1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}

			if (barrier[j] > 1.0e-08)
			{
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					bar_cl[i] = bar_col[j];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}
			}
			else
			{
				is_bar[i] = 0;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			if (is_cont_bar[j+1])
			{
				is_bar[i] = 1;
				bar_lvl[i] = barrier[j+1] * swp_f_df(spot_date, date[i], dom_yc) / swp_f_df(spot_date, date[i], for_yc);
				bar_cl[i] = bar_col[j+1];
			}
			else
			{
				is_bar[i] = 0;
			}

			is_event[i] = 0;
			void_prm[i] = NULL;
		}
		next_d = date[i];
	}

	/* Has bar */
	if (!has_bar)
	{
		free (bar_lvl);
		bar_lvl = NULL;
	}
	
			
	/*	Eventually! call to function */

	*prod_val = dvector (0, num_col-1);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}	
		
	err = FxSabrSL_adi(	nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						void_prm,
						is_event,
						bar_lvl,
						bar_cl,
						is_bar,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi,
						num_col,
						*prod_val,
						0,
						NULL,
						0,
						0,
						0,
						NULL,
						NULL
						);					

	if (err)
	{
		goto FREE_RETURN;
	}
	
	*nprod = num_col;
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1] += xStr.gd->pv_of_past;
		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(	num_und,
								&und_ptr);

		FIRSTFreeEvtDatesFromDeal(	nstp,
									&evt_dts,
									&evt_tms);

		FIRSTFreeMktStruct (&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMFXSABR) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);
	if (is_event) free (is_event);

	if (bar_lvl) free (bar_lvl);
	if (is_bar) free (is_bar);
	if (bar_cl) free (bar_cl);
		
	return err;
}

Err FxSabrSLSmile(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  long		mat_date,
   				  double	*strike,
				  int		nb_strike,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	**res)
{	
	double				*time		= NULL,
						*sig		= NULL,
						*drift		= NULL,
						*date		= NULL;
	
	long				today, spot_date, exe_date;
						
	int					num_col,
						ns;
	
	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;
	
	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL,
						*prod_val		= NULL;

	long				sigma_n_fx, index;

	double				coef, fwd_fx, df, maturity, bs_vol;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	
	TermStruct			*fx_ts;
	
	double				**func_parm	= NULL;

	int					*is_event	= NULL;	

	int					i;
	
	Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) * swp_f_df (today, spot_date, dom_yc)
										/ swp_f_df (today, spot_date, for_yc);

	exe_date = add_unit (mat_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	maturity = (exe_date - today) * YEARS_IN_DAY;

	fwd_fx = get_spot_from_fxund (und) * swp_f_df (spot_date, mat_date, for_yc)
										/ swp_f_df (spot_date, mat_date, dom_yc);
	
	df = swp_f_df(today, exe_date, dom_yc);

	/* discretise in time			*/
	
	ns = 2;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = maturity;
	
	/*	Add revelant vol times	*/
	
	fx_ts = get_ts_from_fxund (und);

	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));
		
	if (!sig || !date || !drift)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
	}

	
			
	/*	Eventually! call to function */

	num_col = nb_strike;	
	
	/* now we need to price the regular option */

	func_parm = dmatrix(0, nstp-1, 0, num_col-1);
	*res = dvector(0, nb_strike-1);
	is_event = calloc(nstp, sizeof(int));

	if (!func_parm || !is_event || !(*res))
	{
		err = "Memory allocation faillure in FxSabrAdiKO";
		goto FREE_RETURN;
	}

	for (i=0; i<num_col; i++)
	{
		func_parm[nstp-1][i] = strike[i];
	}
	
	is_event[nstp-1] = 1;

	err = FxSabrSL_adi(	nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						func_parm,
						is_event,
						NULL,
						NULL,
						NULL,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi_opt,
						num_col,
						*res,
						0,
						NULL,
						0,
						0,
						0,
						NULL,
						NULL
						);

	if (err)
	{
		goto FREE_RETURN;
	}

	for (i=0; i<nb_strike; i++)
	{
		err = srt_f_optimpvol(		
								(*res)[i],
								fwd_fx,
								strike[i],
								maturity,								
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								&bs_vol);

		(*res)[i] = bs_vol;
	}	

		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);	

	if (func_parm) free_dmatrix(func_parm, 0, nstp-1, 0, num_col-1);
	if (is_event) free (is_event);
		
	return err;
}

Err FxSabrSLAdiKO(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  double	maturity,
   				  double	strike,
				  int		is_call,	/* 1 Call, 0: Put */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		is_digital,	/* 1: digital payoff, 0, regular option payoff */
				  int		is_american,
				  double	barrier_up,
				  double	barrier_down,
				  double	rebate_up,
				  double	rebate_down,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res,
				  double	**greeks)
{	
	double				*time		= NULL,
						*sig		= NULL,
						*drift		= NULL,
						*date		= NULL;
	
	long				today, spot_date;
						
	int					num_col,
						ns;
	
	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;
	
	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL,
						*prod_val		= NULL;

	long				sigma_n_fx, index;

	double				coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	
	TermStruct			*fx_ts;
	
	double				*bar_lvl_up		= NULL,
						*bar_lvl_down	= NULL,
						*prod_val_ki	= NULL,
						**greeks_ki		= NULL;	

	int					i;
	
	Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und);	
	
	spot_fx *= swp_f_df (today, spot_date, dom_yc)
				/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

	
	/* discretise in time			*/
	
	ns = 2;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = maturity;
	
	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	bar_lvl_up = (double*) calloc (nstp, sizeof (double));
	bar_lvl_down = (double*) calloc (nstp, sizeof (double));

	if (!sig || !date || !drift || !bar_lvl_up || !bar_lvl_down)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrSLAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		bar_lvl_up[i] = barrier_up * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
		bar_lvl_down[i] = barrier_down * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
	}
			
	/*	Eventually! call to function */

	num_col = 1;

	prod_val = dvector (0, num_col-1);
	*greeks = NULL;
	*greeks = calloc(6, sizeof(double));

	if (!prod_val || !greeks)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	if ((spot_fx > barrier_up) || (spot_fx < barrier_down))
	{
		prod_val[0] = 0.0;		
	}
	else
	{		
		err = FxSabrSL_KOOption(nstp,
								time,
								date,
								nstpfx,
								nstpvol,
								sigma_fx[0],
								drift,
								alpha,
								beta,
								rho,
								lambda,
								floormu,
								strike,
								is_call,
								is_american,
								is_cvx,
								is_digital,
								bar_lvl_up,
								bar_lvl_down,
								rebate_up,
								rebate_down,
								spot_fx,
								dom_yc,
								for_yc,							
								prod_val,
								1,
								*greeks);
	}

	if (err)
	{
		goto FREE_RETURN;
	}

	if (!is_ko)
	{
		/* now we need to price the regular option */		

		greeks_ki = dmatrix(0, 5, 0, 0);
		prod_val_ki = dvector (0, num_col-1);

		for (i=nstp-1; i>=0; i--)
		{
			bar_lvl_up[i] = 1.0E09;
			bar_lvl_down[i] = -1.0E09;
		}

		err = FxSabrSL_KOOption(nstp,
								time,
								date,
								nstpfx,
								nstpvol,
								sigma_fx[0],
								drift,
								alpha,
								beta,
								rho,
								lambda,
								floormu,
								strike,
								is_call,
								is_american,
								is_cvx,
								is_digital,
								bar_lvl_up,
								bar_lvl_down,
								0.0,
								0.0,
								spot_fx,
								dom_yc,
								for_yc,							
								prod_val_ki,
								1,
								*greeks_ki);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		*res = prod_val_ki[0] - prod_val[0];

		(*greeks)[0] = greeks_ki[0][0] - (*greeks)[0];
		(*greeks)[1] = greeks_ki[1][0] - (*greeks)[1];
		(*greeks)[2] = greeks_ki[2][0] - (*greeks)[2];
		(*greeks)[3] = greeks_ki[3][0] - (*greeks)[3];
		(*greeks)[4] = greeks_ki[4][0] - (*greeks)[4];
		(*greeks)[5] = greeks_ki[5][0] - (*greeks)[5];
	}
	else
	{
		*res = prod_val[0];
	}
		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);	

	if (bar_lvl_up) free (bar_lvl_up);
	if (bar_lvl_down) free (bar_lvl_down);
	if (prod_val) free_dvector (prod_val, 0, num_col-1);	

	if (greeks_ki) free_dmatrix(greeks_ki, 0, 5, 0, 0);
	if (prod_val_ki) free_dvector (prod_val_ki, 0, num_col-1);
		
	return err;
}

char *SrtGrfnFxSabrQuadAdi(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMFXSABR		grfn_prm;

	double				next_d;
	
	long				*evt_dts	= NULL;
	double				*evt_tms	= NULL;

	double				*time		= NULL,
						*sig		= NULL,
						*drift	= NULL,
						*date		= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;

	int					ns;

	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;

	int					fx_idx,
						dom_idx,
						for_idx;

	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL;

	long				sigma_n_fx, index;

	double				last_mat, coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	TermStruct			*fx_ts;

	int					*is_bar = NULL;
	double				*bar_lvl = NULL;
	int					*bar_cl = NULL;
	int					has_bar;

	double				next_bar;

	int					i, j, forback;

	double				a, b, c;
	
	Err					err = NULL;
	
	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

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

	if (num_und > 3)
	{
		err = "Product should involve only maximum 3 underlyings";
		goto FREE_RETURN;
	}

	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	fx_idx = dom_idx = for_idx = -1;
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);

		if (num_und > 1 || strcmp (und_ptr[0]->underl_name, underlying))
		{
			err = "Error, the GRFN tableau has too many underlyings involved";
			goto FREE_RETURN;
		}

		fx_idx = 0;
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);

		for (i=0; i<num_und; i++)
		{
			if (!strcmp (und_ptr[i]->underl_name, underlying))
			{
				fx_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, domname))
			{
				dom_idx = i;
			}

			if (!strcmp (und_ptr[i]->underl_name, forname))
			{
				for_idx = i;
			}
		}

		if (fx_idx == -1)
		{
			err = "The Fx underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (dom_idx == -1)
		{
			err = "The domestic underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}

		if (for_idx == -1)
		{
			err = "The foreign underlying is not present in the mdlcomm structure";
			goto FREE_RETURN;
		}
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);	
	
	fx_ts = get_ts_from_fxund (und);

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

	/* discretise in time			*/

	ns = num_evt;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, ns * sizeof (double));

	last_mat = xStr.tms[ns-1];

	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < last_mat)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	is_bar = (int*) calloc (nstp, sizeof (int));
	bar_lvl = (double*) calloc (nstp, sizeof (double));
	bar_cl = (int*) calloc (nstp, sizeof (int));

	if (!void_prm || !is_event || !sig || !date || !drift || !is_bar || !bar_lvl || !bar_cl)
	{
		err = "Memory allocation failure in SrtGrfnMainAdi2";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}
	
	j = xStr.num_evt - 1;
	next_bar = 0.0;
	next_d = evt_dts[j] + 1;
	has_bar = 0;

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		//sig[i] = find_fx_sig (time[i], fx_ts);

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx != -1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			if (barrier[j] > 1.0e-08)
			{
				/* mat_date = add_unit ((long) (date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);*/
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					has_bar = 1;
					bar_cl[i] = bar_col[j];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}				
			}
			else
			{
				is_bar[i] = 0;
			}			

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc (sizeof (grfn_parm_fx_sabr));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;

			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];
			
			if (dom_idx !=-1)
			{
				grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
				grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
				grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];
			}
			else
			{
				grfn_prm->num_dom_df = 0;
				grfn_prm->dom_df_tms = NULL;
				grfn_prm->dom_df_dts = NULL;
			}

			if (for_idx != -1)
			{
				grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
				grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
				grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];
			}
			else
			{
				grfn_prm->num_for_df = 0;
				grfn_prm->for_df_tms = NULL;
				grfn_prm->for_df_dts = NULL;
			}

			if (barrier[j] > 1.0e-08)
			{
				next_bar = barrier[j] * swp_f_df(today, date[i], dom_yc) / swp_f_df(today, date[i], for_yc);
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					bar_cl[i] = bar_col[j];
					bar_lvl[i] = next_bar;
				}
				else
				{
					is_bar[i] = 0;
				}
			}
			else
			{
				is_bar[i] = 0;
			}
			
			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			if (is_cont_bar[j+1])
			{
				is_bar[i] = 1;
				bar_lvl[i] = barrier[j+1] * swp_f_df(spot_date, date[i], dom_yc) / swp_f_df(spot_date, date[i], for_yc);
				bar_cl[i] = bar_col[j+1];
			}
			else
			{
				is_bar[i] = 0;
			}

			is_event[i] = 0;
			void_prm[i] = NULL;
		}
		next_d = date[i];
	}

	/* Has bar */
	if (!has_bar)
	{
		free (bar_lvl);
		bar_lvl = NULL;
	}
	
			
	/*	Eventually! call to function */

	*prod_val = dvector (0, num_col-1);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}	

	/* Transformation from beta and convexity to b and c */
	SabrQuadGetParams(spot_fx, gamma, beta, &a, &b, &c);
		
	err = FxSabrQuad_adi(nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						a,
						b,
						c,
						rho,
						lambda,
						floormu,
						void_prm,
						is_event,
						bar_lvl,
						bar_cl,
						is_bar,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi,
						num_col,
						*prod_val,
						0,
						NULL,
						0,
						0,
						0,
						NULL,
						NULL
						);					

	if (err)
	{
		goto FREE_RETURN;
	}
	
	*nprod = num_col;
	
	/*	Add PV of Past */
	(*prod_val)[num_col-1] += xStr.gd->pv_of_past;
		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(	num_und,
								&und_ptr);

		FIRSTFreeEvtDatesFromDeal(	nstp,
									&evt_dts,
									&evt_tms);

		FIRSTFreeMktStruct (&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMFXSABR) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);
	if (is_event) free (is_event);

	if (bar_lvl) free (bar_lvl);
	if (is_bar) free (is_bar);
	if (bar_cl) free (bar_cl);
		
	return err;
}

Err FxSabrQuadSmile(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  long		mat_date,
   				  double	*strike,
				  int		nb_strike,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	**res)
{	
double				*time		= NULL,
					*sig		= NULL,
					*drift		= NULL,
					*date		= NULL;

long				today, spot_date, exe_date;
					
int					num_col,
					ns;

double				spot_fx;
char				*domname,
					*forname,
					*dom_yc,
					*for_yc;

double				*sigma_fx		= NULL,
					*sigma_time_fx	= NULL,
					*prod_val		= NULL;

long				sigma_n_fx, index;

double				coef, fwd_fx, df, maturity, bs_vol;

SrtUndPtr			*und_ptr	= NULL,
					und			= NULL,
					dom_und		= NULL,
					for_und		= NULL;

TermStruct			*fx_ts;

double				**func_parm	= NULL;

int					*is_event	= NULL;	

int					i;

double				a, b, c;

Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und) * swp_f_df (today, spot_date, dom_yc)
										/ swp_f_df (today, spot_date, for_yc);
	
	exe_date = add_unit (mat_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	maturity = (exe_date - today) * YEARS_IN_DAY;

	fwd_fx = get_spot_from_fxund (und) * swp_f_df (spot_date, mat_date, for_yc)
										/ swp_f_df (spot_date, mat_date, dom_yc);
	
	df = swp_f_df(today, exe_date, dom_yc);

	/* discretise in time			*/
	
	ns = 2;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = maturity;
	
	/*	Add revelant vol times	*/
	
	fx_ts = get_ts_from_fxund (und);

	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));
		
	if (!sig || !date || !drift)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
	}

	
			
	/*	Eventually! call to function */

	num_col = nb_strike;	
	
	/* now we need to price the regular option */

	func_parm = dmatrix(0, nstp-1, 0, num_col-1);
	*res = dvector(0, nb_strike-1);
	is_event = calloc(nstp, sizeof(int));

	if (!func_parm || !is_event || !(*res))
	{
		err = "Memory allocation faillure in FxSabrAdiKO";
		goto FREE_RETURN;
	}

	for (i=0; i<num_col; i++)
	{
		func_parm[nstp-1][i] = strike[i];
	}
	
	is_event[nstp-1] = 1;

	/* Transformation from beta and convexity to b and c */
	SabrQuadGetParams(spot_fx, gamma, beta, &a, &b, &c);

	err = FxSabrQuad_adi(nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						a,
						b,
						c,
						rho,
						lambda,
						floormu,
						func_parm,
						is_event,
						NULL,
						NULL,
						NULL,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi_opt,
						num_col,
						*res,
						0,
						NULL,
						0,
						0,
						0,
						NULL,
						NULL
						);

	if (err)
	{
		goto FREE_RETURN;
	}

	for (i=0; i<nb_strike; i++)
	{
		err = srt_f_optimpvol(		
								(*res)[i],
								fwd_fx,
								strike[i],
								maturity,								
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								&bs_vol);

		(*res)[i] = bs_vol;
	}	

		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);	

	if (func_parm) free_dmatrix(func_parm, 0, nstp-1, 0, num_col-1);
	if (is_event) free (is_event);
		
	return err;
}

Err FxSabrQuadAdiKO(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  double	maturity,
   				  double	strike,
				  int		is_call,	/* 1 Call, 0: Put */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		is_digital,	/* 1: digital payoff, 0, regular option payoff */
				  int		is_american,
				  double	barrier_up,
				  double	barrier_down,
				  double	rebate_up,
				  double	rebate_down,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res,
				  double	**greeks)
{	
	double				*time		= NULL,
						*sig		= NULL,
						*drift		= NULL,
						*date		= NULL;
	
	long				today, spot_date;
						
	int					num_col,
						ns;
	
	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;
	
	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL,
						*prod_val		= NULL;

	long				sigma_n_fx, index;

	double				coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	
	TermStruct			*fx_ts;
	
	double				*bar_lvl_up		= NULL,
						*bar_lvl_down	= NULL,
						*prod_val_ki	= NULL,
						**greeks_ki		= NULL;

	double				a, b, c;
	long				date_temp, settlmt_date;
	int					i;
	
	Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und);	
	
	spot_fx *= swp_f_df (today, spot_date, dom_yc)
				/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

	
	/* discretise in time			*/
	
	ns = 2;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = maturity;
	
	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity)
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	bar_lvl_up = (double*) calloc (nstp, sizeof (double));
	bar_lvl_down = (double*) calloc (nstp, sizeof (double));

	if (!sig || !date || !drift || !bar_lvl_up || !bar_lvl_down)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrSLAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		date_temp = add_unit ((long) (date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);
		bar_lvl_up[i] = barrier_up * swp_f_df(today, date_temp, dom_yc) / swp_f_df(today, date_temp, for_yc);
		bar_lvl_down[i] = barrier_down * swp_f_df(today, date_temp, dom_yc) / swp_f_df(today, date_temp, for_yc);
	}

	settlmt_date = add_unit ((long) (today + maturity * DAYS_IN_YEAR), 2, SRT_BDAY, MODIFIED_SUCCEEDING);
			
	/*	Eventually! call to function */

	/* Transformation from beta and convexity to b and c */
	SabrQuadGetParams(spot_fx, gamma, beta, &a, &b, &c);

	num_col = 1;

	prod_val = dvector (0, num_col-1);
	*greeks = NULL;
	*greeks = calloc(6, sizeof(double));

	if (!prod_val || !greeks)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	if ((spot_fx > barrier_up) || (spot_fx < barrier_down))
	{
		prod_val[0] = 0.0;		
	}
	else
	{		
		err = FxSabrQuad_KOOption(nstp,
								time,
								date,
								nstpfx,
								nstpvol,
								sigma_fx[0],
								drift,
								alpha,
								a,
								b,
								c,
								rho,
								lambda,
								floormu,
								settlmt_date,
								strike,
								is_call,
								is_american,
								is_cvx,
								is_digital,
								bar_lvl_up,
								bar_lvl_down,
								rebate_up,
								rebate_down,
								spot_fx,
								dom_yc,
								for_yc,	
								0,
								0,
								0,
								prod_val,
								1,
								*greeks);
	}

	if (err)
	{
		goto FREE_RETURN;
	}

	if (!is_ko)
	{
		/* now we need to price the regular option */		

		greeks_ki = dmatrix(0, 5, 0, 0);
		prod_val_ki = dvector (0, num_col-1);

		for (i=nstp-1; i>=0; i--)
		{
			bar_lvl_up[i] = 1.0E09;
			bar_lvl_down[i] = -1.0E09;
		}

		err = FxSabrQuad_KOOption(nstp,
								time,
								date,
								nstpfx,
								nstpvol,
								sigma_fx[0],
								drift,
								alpha,
								a,
								b,
								c,
								rho,
								lambda,
								floormu,
								settlmt_date,
								strike,
								is_call,
								is_american,
								is_cvx,
								is_digital,
								bar_lvl_up,
								bar_lvl_down,
								0.0,
								0.0,
								spot_fx,
								dom_yc,
								for_yc,
								0,
								0,
								0,
								prod_val_ki,
								1,
								*greeks_ki);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		*res = prod_val_ki[0] - prod_val[0];

		(*greeks)[0] = greeks_ki[0][0] - (*greeks)[0];
		(*greeks)[1] = greeks_ki[1][0] - (*greeks)[1];
		(*greeks)[2] = greeks_ki[2][0] - (*greeks)[2];
		(*greeks)[3] = greeks_ki[3][0] - (*greeks)[3];
		(*greeks)[4] = greeks_ki[4][0] - (*greeks)[4];
		(*greeks)[5] = greeks_ki[5][0] - (*greeks)[5];
	}
	else
	{
		*res = prod_val[0];
	}
		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);	

	if (bar_lvl_up) free (bar_lvl_up);
	if (bar_lvl_down) free (bar_lvl_down);
	if (prod_val) free_dvector (prod_val, 0, num_col-1);	

	if (sigma_fx) free(sigma_fx);
	if (sigma_time_fx) free(sigma_time_fx);

	if (greeks_ki) free_dmatrix(greeks_ki, 0, 5, 0, 0);
	if (prod_val_ki) free_dvector (prod_val_ki, 0, num_col-1);
		
	return err;
}

Err FxSabrQuadAdi_MultiKO(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */	
				  int		nb_product,
				  double	*notional,
				  long		*exe_date,
				  long		*settlmt_date,
   				  double	*strike,
				  int		*is_call,	/* 1 Call, 0: Put */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		*is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		*is_digital,	/* 1: digital payoff, 0, regular option payoff */				  
				  double	barrier_up,
				  double	barrier_down,
				  double	*rebate_up,
				  double	*rebate_down,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res,
				  double	**greeks)
{	
	double				*time			= NULL,
						*sig			= NULL,
						*drift			= NULL,
						*date			= NULL;

	double				*maturity		= NULL,
						*settlment		= NULL;

	int					*eval_evt		= NULL;
	
	long				today, spot_date;
						
	int					num_col,
						ns;
	
	double				spot_fx;
	char				*domname,
						*forname,
						*dom_yc,
						*for_yc;
	
	double				*sigma_fx		= NULL,
						*sigma_time_fx	= NULL,
						*prod_val		= NULL;

	long				sigma_n_fx, index;

	double				coef;

	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL,
						dom_und		= NULL,
						for_und		= NULL;
	
	TermStruct			*fx_ts;
	
	double				*bar_lvl_up		= NULL,
						*bar_lvl_down	= NULL,
						*new_reb_down	= NULL,
						*new_reb_up		= NULL,
						*prod_val_ki	= NULL,
						**greeks_ki		= NULL;

	double				a, b, c;
	long				date_temp;
	int					i;
	double				sum_reb_down, sum_reb_up;
	
	Err					err = NULL;
	
	/* look for the underlying name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	/* look for the today date */
	today = get_today_from_underlying (und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	if (get_underlying_type (und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}
	
	if (get_mdltype_from_fxund (und) == BLACK_SCHOLES)
	{
		dom_yc = (char *) get_domname_from_fxund (und);
		for_yc = (char *) get_forname_from_fxund (und);
	}
	else if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domname = get_domname_from_fxund (und);
		dom_und = lookup_und (domname);
		if (!dom_und)
		{
			err = serror ("Couldn't find underlying named %s", domname);
			goto FREE_RETURN;
		}
		dom_yc = get_ycname_from_irund (dom_und);

		forname = get_forname_from_fxund (und);
		for_und = lookup_und (forname);
		if (!for_und)
		{
			err = serror ("Couldn't find underlying named %s", forname);
			goto FREE_RETURN;
		}
		for_yc = get_ycname_from_irund (for_und);
	}
	else
	{
		err = "Model must be BLACK_SCHOLES or FX_STOCH_RATES";
		goto FREE_RETURN;
	}	

	/* Get the term structure */

	err = display_FX1F_TermStruct(	underlying,
									&sigma_n_fx,
									&sigma_time_fx,
									&sigma_fx);

	if (err)
	{
		goto FREE_RETURN;
	}		


	spot_fx = get_spot_from_fxund (und);	
	
	spot_fx *= swp_f_df (today, spot_date, dom_yc)
				/ swp_f_df (today, spot_date, for_yc);
	
	fx_ts = get_ts_from_fxund (und);

	
	/* Select relevant products */	
	while (nb_product > 0 && exe_date[0] < today)
	{
		nb_product--;
		notional++;
		exe_date++;
		settlmt_date++;
		strike++;
		is_call++;
		is_cvx++;
		is_digital++;
		rebate_up++;
		rebate_down++;
	}

	if (nb_product == 0)
	{
		/* all products are past */
		prod_val[0] = 0.0;
		goto FREE_RETURN;
	}	

	/* remove 0 notional */
	i = nb_product - 1;

	/* MODIF */
	/* while (i >= 0 && fabs(notional[i]) < 1.0E-08)
	{
		i--;
	} */

	nb_product = i + 1;

	if (nb_product == 0)
	{
		/* all products are past */
		prod_val[0] = 0.0;
		goto FREE_RETURN;
	}

	
	/* skip first options */


	/* discretise in time			*/
	maturity = (double*) calloc (nb_product, sizeof(double));
	settlment = (double*) calloc (nb_product, sizeof(double));

	if (!maturity || !settlment)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrSLAdiBar";
		goto FREE_RETURN;
	}

	for (i=0; i<nb_product; i++)
	{
		maturity[i] = (exe_date[i] - today) * YEARS_IN_DAY;
		settlment[i] = (settlmt_date[i] - today) * YEARS_IN_DAY;
	}
	
	ns = 1;

	time = (double*) calloc (ns, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	
	/*	Add exercise dates	*/	
	for (i=0; i<nb_product; i++)
	{		
		num_f_add_number (&ns, &time, maturity[i]);	
	}
	
	/*	Add revelant vol times	*/
	
	for (i=0; i<sigma_n_fx; i++)
	{
		if (sigma_time_fx[i] < maturity[nb_product-1])
		{
			num_f_add_number (&ns, &time, sigma_time_fx[i]);
		}
		else
		{
			break;
		}
	}

	num_f_sort_vector(ns, time);
	num_f_unique_vector(&ns, time);

	/*	Fill the time vector */
	err = fill_time_vector (	&time, 
								&ns,
								0,
								NULL,
								0,
								NULL, 
								nstp);				
	if (err)
	{
		goto FREE_RETURN;
	}
	
	nstp = ns;

	date = (double *) calloc (nstp, sizeof (double));
	sig = (double *) calloc (nstp, sizeof (double));
	drift = (double *) calloc (nstp, sizeof (double));

	bar_lvl_up = (double*) calloc (nstp, sizeof (double));
	bar_lvl_down = (double*) calloc (nstp, sizeof (double));	

	new_reb_down = (double*) calloc (nb_product, sizeof(double));
	new_reb_up = (double*) calloc (nb_product, sizeof(double));

	eval_evt = (int*) calloc (nstp, sizeof(int));

	if (!sig || !date || !drift || !bar_lvl_up || !bar_lvl_down || !eval_evt ||
		!new_reb_down || !new_reb_up)
	{
		err = "Memory allocation failure in SrtGrfnFxSabrSLAdiBar";
		goto FREE_RETURN;
	}

	/* first part */

	i = 0;
	index = 0;

	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}

	/* middle and end part */
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{
			/* middle part */
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);

			index++;

			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}

		}
		else
		{
			/* end part */
			coef = 0.0;

			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	index = nb_product - 1;

	sum_reb_down = 0.0;
	sum_reb_up = 0.0;

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = today + time[i] * DAYS_IN_YEAR;
		date_temp = add_unit ((long) (date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);
		bar_lvl_up[i] = barrier_up * swp_f_df(today, date_temp, dom_yc) / swp_f_df(today, date_temp, for_yc);
		bar_lvl_down[i] = barrier_down * swp_f_df(today, date_temp, dom_yc) / swp_f_df(today, date_temp, for_yc);

		eval_evt[i] = 0;

		while (index >= 0 && fabs(time[i] - maturity[index]) < 1.0E-08)
		{
			/* if (fabs(notional[index]) > 1.0E-08)
			{ */
				eval_evt[i] += 1;
			
				sum_reb_down += rebate_down[index] * notional[index];
				sum_reb_up += rebate_up[index] * notional[index];
			/*} */

			new_reb_down[index] = sum_reb_down;
			new_reb_up[index] = sum_reb_up;			

			index--;
		}
	}	
			
	/*	Eventually! call to function */

	/* Transformation from beta and convexity to b and c */
	SabrQuadGetParams(spot_fx, gamma, beta, &a, &b, &c);

	num_col = 1;

	prod_val = dvector (0, num_col-1);
	*greeks = NULL;
	*greeks = calloc(6, sizeof(double));

	if (!prod_val || !greeks)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	if ((spot_fx > barrier_up) || (spot_fx < barrier_down))
	{
		prod_val[0] = 0.0;		
	}
	else
	{		
		err = FxSabrQuad_KO_MultiOption(nstp,
										time,
										date,
										eval_evt,
										nstpfx,
										nstpvol,
										sigma_fx[0],
										drift,
										alpha,
										a,
										b,
										c,
										rho,
										lambda,
										floormu,
										nb_product,
										notional,
										exe_date,
										settlmt_date,
										strike,
										is_call,
										is_cvx,
										is_digital,
										bar_lvl_up,
										bar_lvl_down,
										new_reb_up,
										new_reb_down,
										spot_fx,
										dom_yc,
										for_yc,	
										0,
										0,
										0,
										prod_val,
										1,
										*greeks);
	}

	if (err)
	{
		goto FREE_RETURN;
	}

	if (!is_ko)
	{
		/* now we need to price the regular option */		

		greeks_ki = dmatrix(0, 5, 0, 0);
		prod_val_ki = dvector (0, num_col-1);

		for (i=nstp-1; i>=0; i--)
		{
			bar_lvl_up[i] = 1.0E09;
			bar_lvl_down[i] = -1.0E09;
		}

		err = FxSabrQuad_KO_MultiOption(nstp,
										time,
										date,
										eval_evt,
										nstpfx,
										nstpvol,
										sigma_fx[0],
										drift,
										alpha,
										a,
										b,
										c,
										rho,
										lambda,
										floormu,
										nb_product,
										notional,
										exe_date,
										settlmt_date,
										strike,
										is_call,
										is_cvx,
										is_digital,
										bar_lvl_up,
										bar_lvl_down,
										new_reb_up,
										new_reb_down,
										spot_fx,
										dom_yc,
										for_yc,	
										0,
										0,
										0,
										prod_val_ki,
										1,
										*greeks);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		*res = prod_val_ki[0] - prod_val[0];

		(*greeks)[0] = greeks_ki[0][0] - (*greeks)[0];
		(*greeks)[1] = greeks_ki[1][0] - (*greeks)[1];
		(*greeks)[2] = greeks_ki[2][0] - (*greeks)[2];
		(*greeks)[3] = greeks_ki[3][0] - (*greeks)[3];
		(*greeks)[4] = greeks_ki[4][0] - (*greeks)[4];
		(*greeks)[5] = greeks_ki[5][0] - (*greeks)[5];
	}
	else
	{
		*res = prod_val[0];
	}
		
FREE_RETURN:
		
	if (date) free (date);
	if (time) free (time);
	if (sig) free (sig);
	if (drift) free (drift);

	if (sigma_fx) free(sigma_fx);
	if (sigma_time_fx) free(sigma_time_fx);

	if (bar_lvl_up) free (bar_lvl_up);
	if (bar_lvl_down) free (bar_lvl_down);
	if (prod_val) free_dvector (prod_val, 0, num_col-1);	

	if (greeks_ki) free_dmatrix(greeks_ki, 0, 5, 0, 0);
	if (prod_val_ki) free_dvector (prod_val_ki, 0, num_col-1);

	if (maturity) free(maturity);
	if (settlment) free(settlment);

	if (new_reb_down) free(new_reb_down);
	if (new_reb_up) free(new_reb_up);

	if (eval_evt) free(eval_evt);
		
	return err;
}