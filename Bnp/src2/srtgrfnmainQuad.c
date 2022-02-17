/**********************************************************************
 *      Name: SrtGrfnMainQuad.c                                       * 
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: L.C				                                      *
 *      Date: 18/12/00                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------* 
 **********************************************************************/

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "SrtAccess.h"
#include "BGMEval.h"

#define MAX_TIME	0.02
#define	MAX_Z		2.0
#define NB_POINTS	200
#define MAX_STP		3000

char *SrtGrfn3DFXQuadTree(
				  char		*und3dfx,
				  double	beta,
				  double	gamma,
				  double	sig0,
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
				  long		num_stp, 
				  int		*num_prod, 
				  int		discount,
				  double	*vols,
				  double	*vols_time,
				  int		nbVols, 
				  double	**prod_val)
{
	int				free_str = 0;
	FIRSTAllMkts	xStr;
	SrtGrfnParam	defParm;
	int				forback;
	long			nstp;
	
	double			*time = NULL,
					*date = NULL;
	int				*vol_change = NULL;
	double			*dom_vol = NULL,
					*for_vol = NULL,
					*fx_vol = NULL,
					*sig_fx_approx = NULL,
					*fx_fwd_approx = NULL;
	
	double			*dom_ifr = NULL,		
					*dom_fwd = NULL,
					*dom_var = NULL,
					*for_ifr = NULL,
					*for_fwd = NULL,
					*for_var = NULL,
					*fx_fwd = NULL,
					*fx_fwd0 = NULL,
					*fx_var = NULL;
					
	void			**void_prm = NULL;
	GRFNPARMTREE	grfn_prm;
	int				*is_event = NULL;

	long			today,
					spot_date;
	int				i, j;
	SrtUndPtr		fx_und, dom_und, for_und;
	TermStruct		*fx_ts, *dom_ts, *for_ts;
	SrtCorrLstPtr	sCorrlist;
	char			*domname,
					*forname;
	double			dom_lam, 
					for_lam, 
					corr_dom_for,
					corr_dom_fx,
					corr_for_fx;
	double			spot_fx;
	char			*dom_yc,
					*for_yc;
	int				fx_idx,
					dom_idx,
					for_idx;

	int				num_vol_times;
	double			*vol_times = NULL;

	int				num_bar_times;
	double			*bar_times = NULL;

	int				*is_bar = NULL;
	double			*bar_lvl = NULL;
	int				*bar_cl = NULL;
	double			next_bar;
	double			alpha2, beta2, gamma2;
	double			param[3];

	LOCAL_MODEL_FUNC	*model_func = NULL;
	LOCAL_MODEL_PARAM	*model_param = NULL;

	double			temp;

	clock_t			t1, t2;

	Err				err = NULL;


	
	t1 = clock();

	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */
	
	err = srt_f_set_default_GrfnParams (&defParm);
	defParm.min_nodes_per_path = num_stp;

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
		und3dfx,
		&defParm,
		&forback,
		&xStr);
	
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved and their term structures */

	fx_und = lookup_und (und3dfx);
	
	if (!fx_und)
	{
		err = serror ("Couldn't find underlying named %s", und3dfx);
		goto FREE_RETURN;
	}
	
	today = get_today_from_underlying (fx_und);

	if (get_underlying_type (fx_und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", und3dfx);
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (fx_und) != FX_STOCH_RATES)
	{
		err = serror ("Underlying %s is not of type FX Stoch Rates", und3dfx);
		goto FREE_RETURN;
	}

	fx_ts = get_ts_from_fxund (fx_und);

	domname = get_domname_from_fxund (fx_und);
	dom_und = lookup_und (domname);
	if (!dom_und)
	{
		err = serror ("Couldn't find underlying named %s", domname);
		goto FREE_RETURN;
	}
	dom_ts = get_ts_from_irund (dom_und);

	forname = get_forname_from_fxund (fx_und);
	for_und = lookup_und (forname);
	if (!for_und)
	{
		err = serror ("Couldn't find underlying named %s", forname);
		goto FREE_RETURN;
	}
	for_ts = get_ts_from_irund (for_und);

	/*	Next, get the time steps */

	/*	Copy event dates */
	nstp = xStr.num_evt;
	while (nstp >= 1 && xStr.evt[nstp-1].evt == NULL)
	{
		nstp--;
	}
	if (nstp < 1)
	{
		err = "No event in Tableau";
		goto FREE_RETURN;
	}
	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	/*	Get vol dates */
	err = compute_vol_times (und3dfx, &num_vol_times, &vol_times, time[nstp-1]);
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get barrier dates */
	num_bar_times = 0;
	bar_times = (double*) calloc (nstp, sizeof (double));
	if (!bar_times)
	{
		err = "Memory allocation error (2) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	for (i=0; i<xStr.num_evt; i++)
	{
		if (eventdates[i] > today 
			&& eventdates[i] <= xStr.dts[nstp-1]
			&& barrier[i] > 1.0e-08)
		{
			bar_times[num_bar_times] = (eventdates[i] - today) * YEARS_IN_DAY;
			num_bar_times++;
		}
	}
	
	/*	Fill the time vector */

	err = fill_time_vector (	&time, 
								&nstp, 
								num_bar_times,
								bar_times,
								num_vol_times, 
								vol_times, 
								num_stp);
	
	if (err)
	{
		goto FREE_RETURN;
	}

	date = (double*) calloc (nstp, sizeof (double));
	if (!date)
	{
		err = "Memory allocation error (3) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	
	for (i=0; i<nstp; i++)
	{		
		date[i] = today + DAYS_IN_YEAR * time[i];

		if (i > 0 && date[i] - date[i-1] >= 1)
		{
			date[i] = (long) (date[i] + 1.0e-08);
			time[i] = YEARS_IN_DAY * (date[i] - today);	
		}
	}	
						  
	/*	Fill the model parameters as required by tree_main_3dfx */

	vol_change = (int*) calloc (nstp, sizeof (int));
	dom_vol = (double*) calloc (nstp, sizeof (double));
	for_vol = (double*) calloc (nstp, sizeof (double));
	fx_vol = (double*) calloc (nstp, sizeof (double));
	
	if (!vol_change || !dom_vol || !for_vol || !fx_vol)
	{
		err = "Memory allocation error (4) in SrtGrfn3DFXQuadTree";
		goto FREE_RETURN;
	}

	vol_change[nstp-1] = 1;
	dom_vol[nstp-1] = find_sig (time[nstp-1], dom_ts);
	for_vol[nstp-1] = find_sig (time[nstp-1], for_ts);
	fx_vol[nstp-1] = find_fx_sig (time[nstp-1], fx_ts);

	for (i=nstp-2; i>=0; i--)
	{
		dom_vol[i] = find_sig (time[i], dom_ts);
		for_vol[i] = find_sig (time[i], for_ts);
		fx_vol[i] = find_fx_sig (time[i], fx_ts);

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

	/*	Get lambdas and correls */

	err = get_lambda_from_ir_ts (dom_ts, &dom_lam);
	err = get_lambda_from_ir_ts (for_ts, &for_lam);
	
	if (err)
	{
		goto FREE_RETURN;
	}

	sCorrlist = srt_f_GetTheCorrelationList();
	if (!sCorrlist->head->element)
	{
		err = "correlation list improperly initialised";
		goto FREE_RETURN;
	}
		
	/*	It is assumed that correlation has no term structure */
	err = srt_f_get_corr_from_CorrList(
					sCorrlist,		
					domname,									    
					forname,									 
					1.0,
					&corr_dom_for);

	err = srt_f_get_corr_from_CorrList(
					sCorrlist,		
					domname,									    
					und3dfx,									 
					1.0,
					&corr_dom_fx);

	err = srt_f_get_corr_from_CorrList(
					sCorrlist,		
					forname,									    
					und3dfx,									 
					1.0,
					&corr_for_fx);

	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get Fx spot and yield curves */
	
	dom_yc = get_ycname_from_irund (dom_und);
	for_yc = get_ycname_from_irund (for_und);

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	spot_fx = get_spot_from_fxund (fx_und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);

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
		err = "Memory allocation error (5) in SrtGrfn3DFXQuadTree";
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

	/* reparametrisation of the model */
	alpha2 = 1.0 - beta + gamma;	
	beta2 = beta - 2.0 * gamma;
	gamma2 = gamma;

	param[0] = alpha2;
	param[1] = gamma2;
	param[2] = beta2;

	/* computation of the model parameters */

	err = allocate_local_model(nstp, &model_func, &model_param);
	//err = allocate_local_model_cust(NB_POINTS, nstp, &model_func, &model_param);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = fill_local_model(model_func, model_param, nstp, time, alpha2, beta2, gamma2, sig0, fx_fwd0);
	//err = fill_local_model_cust(model_func, model_param, NB_POINTS, time, nstp, alpha2, beta2, gamma2, sig0, fx_fwd0);
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

	/* calculates the fx_var of U as if it was a Martingale 
	t = 0;
	fx_var[0] = 0;
	for (i=1; i<nstp; i++)
	{
		prev_t = t;
		t = time[i];
		fx_var[i] = fx_var[i-1] + fx_vol[i] * fx_vol[i] * (t - prev_t);		 		
	}
	*/

	
	if (nbVols > 0)
	{		
		for (i=1; i<nstp; i++)
		{
			fx_var[i] = interp(vols_time, vols, nbVols, time[i], 5, &temp);
			fx_var[i] *= fx_var[i] * time[i];
		}
	}
	

	/*
	for (i=1; i<nstp; i++)
	{
		fx_vol[i] /= fx_fwd0[i-1];
	}
	*/

	/*	Fill product structure */
	strupper (und3dfx);
	strip_white_space (und3dfx);
	strupper (domname);
	strip_white_space (domname);
	strupper (forname);
	strip_white_space (forname);
	for (i=0; i<xStr.num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	fx_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, und3dfx))
		{
			fx_idx = i;
		}
	}
	if (fx_idx == -1)
	{
		err = "The Fx underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	dom_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, domname))
		{
			dom_idx = i;
		}
	}
	if (dom_idx == -1)
	{
		err = "The domestic underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	for_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, forname))
		{
			for_idx = i;
		}
	}
	if (for_idx == -1)
	{
		err = "The foreign underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	is_event = (int*) calloc (nstp, sizeof (int));
	void_prm = (void**) calloc (nstp, sizeof (void*));

	is_bar = (int*) calloc (nstp, sizeof (int));
	bar_lvl = (double*) calloc (nstp, sizeof (double));
	bar_cl = (int*) calloc (nstp, sizeof (int));

	if (!is_event || !void_prm || !is_bar || !bar_lvl || !bar_cl)
	{
		err = "Memory allocation error (6) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;
	next_bar = 0.0;
	for (i=nstp-1; i>=0; i--)
	{
		if (j >= 0 && fabs (date[i] - xStr.dts[j]) < 1.0e-04)
		{
			grfn_prm = malloc (sizeof (grfn_parm_tree));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
			
			if (barrier[j] > 1.0e-08)
			{
				next_bar = model_func -> Z_to_U(i, model_func -> S_to_Z(i, barrier[j], fx_fwd0[i], model_param),
												model_param);
			
				if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					bar_cl[i] = bar_col[j];
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
		else
		if (j>0 && xStr.am[j-1])
		{
			grfn_prm = malloc (sizeof (grfn_parm_tree));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j - 1;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[j-1].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j-1].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j-1].evt->dfd[fx_idx];

			grfn_prm->num_dom_df = xStr.evt[j-1].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j-1].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j-1].evt->dfd[dom_idx];

			grfn_prm->num_for_df = xStr.evt[j-1].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j-1].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j-1].evt->dfd[for_idx];

			if (barrier[j-1] > 1.0e-08)
			{
				next_bar = model_func -> Z_to_U(i, model_func -> S_to_Z(i, barrier[j-1], fx_fwd0[i], model_param),
												model_param);
				
				if (bar_col[j-1] >= 0 && bar_col[j-1] <= xStr.num_cols-1)
				{
					is_bar[i] = 1;
					bar_cl[i] = bar_col[j-1];
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
			is_event[i] = 0;
			void_prm[i] = NULL;
			is_bar[i] = 0;
		}

		if (fabs (next_bar) > 1.0e-08)
		{
			bar_lvl[i] = next_bar;
		}
		else
		{
			bar_lvl[i] = model_func -> Z_to_U(i, model_func -> S_to_Z(i, fx_fwd0[i], fx_fwd0[i], model_param), model_param);
		}

		if (i > 0)
		{
			/* not done, should be the backward expected barrier */
			next_bar = 0;
		}
	}

	is_bar[0] = 0;	
	bar_lvl[0] = model_func -> Z_to_U(0, 0.0, model_param);
	bar_cl[0] = -1;

	is_bar[1] = 0;
	bar_lvl[1] = model_func -> Z_to_U(1, model_func -> S_to_Z(i, fx_fwd0[1], fx_fwd0[1], model_param), model_param);
	bar_cl[1] = -1;
	
	/*	Eventually! call to function */

	*num_prod = xStr.num_cols;
	*prod_val = (double*) calloc (*num_prod+1, sizeof (double));

	if (!(*prod_val))
	{
		err = "Memory allocation error (7) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Number of times steps, required: %d, actual: %d", num_stp, nstp);
	
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
						void_prm, 
						is_event,
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
						grfn_payoff_4_3dfxQuad_tree,
						*num_prod,
						discount,
						*prod_val);

	/*	Add PV of Past */
	(*prod_val)[*num_prod-1] += xStr.gd->pv_of_past;
	

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

	if (vol_times) free (vol_times);
	if (bar_times) free (bar_times);

	if (is_event) free (is_event);

	if (bar_lvl) free (bar_lvl);
	if (is_bar) free (is_bar);
	if (bar_cl) free (bar_cl);

	free_local_model_param(model_func, model_param);

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMTREE) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (free_str)
	{
		FIRSTFreeMktStruct (&xStr);
	}

	return err;
}


char *SrtGrfn3DFXQuadMc(
				  char		*und3dfx,
				  double	beta,
				  double	gamma,
				  double	sig0,
				  int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,				  
				  long		num_paths,  
				  long		num_stp,
				  int		do_pecs,
				  double	***prod_val)
{
	int				free_str = 0;
	FIRSTAllMkts	xStr;
	SrtGrfnParam	defParm;
	int				forback;
	int				flag = 0;
	long			nstp;
	
	double			*time = NULL,
					*date = NULL;
		
	double			*dom_ifr = NULL,		
					*dom_vol = NULL,
					*dom_phi = NULL,
					
					*for_ifr = NULL,
					*for_vol = NULL,
					*for_phi = NULL,
										
					*fx_vol = NULL,
					*fx_fwd	= NULL,
					
					*corr_dom_for	= NULL,
					*corr_dom_fx	= NULL,
					*corr_for_fx	= NULL;
					
	double			dom_phi_at_t, for_phi_at_t;

	void			**void_prm = NULL;
	GRFNPARMMC		grfn_prm;

	long			today,
					spot_date;
	int				i, j, num_col;
	SrtUndPtr		fx_und, dom_und, for_und;
	TermStruct		*fx_ts, *dom_ts, *for_ts;	
	char			*domname,
					*forname;
	double			dom_lam, 
					for_lam;
	double			spot_fx;
	char			*dom_yc,
					*for_yc;
	int				fx_idx,
					dom_idx,
					for_idx;

	double			*sigma_date_dom = NULL,
					*sigma_dom		= NULL,
					*tau_date_dom	= NULL,
					*tau_dom		= NULL,
					*sigma_date_for = NULL,
					*sigma_for		= NULL,
					*tau_date_for	= NULL,
					*tau_for		= NULL,
					*sigma_date_fx	= NULL,
					*sigma_fx		= NULL,
					*dom_for_cov	= NULL,
					*dom_fx_cov		= NULL,
					*for_fx_cov		= NULL,
					*correl_mat		= NULL,
					*correl_dom_for	= NULL,
					*correl_dom_fx	= NULL,
					*correl_for_fx	= NULL;	
	
	double			*merge_dates	= NULL,
					*sig_dom		= NULL,
					*sig_for		= NULL,
					*sig_fx			= NULL;

	long			nb_merge_dates, nb_correl;

	long			sigma_n_dom,
					tau_n_dom,
					sigma_n_for,
					tau_n_for,
					sigma_n_fx;
		
	int				*has_evt = NULL;

	double			param[3];

	double			alpha2, beta2, gamma2;

	LOCAL_MODEL_FUNC	*model_func = NULL;
	LOCAL_MODEL_PARAM	*model_param = NULL;

	clock_t			t1, t2;
	
	Err				err = NULL;

	
	t1 = clock();

	/*	Initialise the GRFN tableau */

	/*	First, initialise the param struct */
	
	err = srt_f_set_default_GrfnParams (&defParm);
	defParm.num_MCarlo_paths = num_paths;
	defParm.max_time_per_slice = 1000;
	defParm.min_nodes_per_path = 1;
	defParm.force_mc = 1;
	defParm.jumping = 1;

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
							und3dfx,
							&defParm,
							&forback,
							&xStr);
	
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved and their term structures */

	fx_und = lookup_und (und3dfx);
	
	if (!fx_und)
	{
		err = serror ("Couldn't find underlying named %s", und3dfx);
		goto FREE_RETURN;
	}
	
	today = get_today_from_underlying (fx_und);

	if (get_underlying_type (fx_und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", und3dfx);
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (fx_und) != FX_STOCH_RATES)
	{
		err = serror ("Underlying %s is not of type FX Stoch Rates", und3dfx);
		goto FREE_RETURN;
	}

	fx_ts = get_ts_from_fxund (fx_und);

	domname = get_domname_from_fxund (fx_und);
	dom_und = lookup_und (domname);
	if (!dom_und)
	{
		err = serror ("Couldn't find underlying named %s", domname);
		goto FREE_RETURN;
	}
	dom_ts = get_ts_from_irund (dom_und);

	forname = get_forname_from_fxund (fx_und);
	for_und = lookup_und (forname);
	if (!for_und)
	{
		err = serror ("Couldn't find underlying named %s", forname);
		goto FREE_RETURN;
	}
	for_ts = get_ts_from_irund (for_und);

	num_col = xStr.num_cols;
	/*	Next, get the time steps */

	/*	Copy event dates */
	nstp = xStr.num_evt;
	while (nstp >= 1 && xStr.evt[nstp-1].evt == NULL)
	{
		nstp--;
	}
	if (nstp < 1)
	{
		err = "No event in Tableau";
		goto FREE_RETURN;
	}

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	/* fill time vector */
	err = fill_time_vector (	
							&time, 
							&nstp, 
							0,
							NULL,
							0, 
							NULL, 
							num_stp);

	date = (double*) calloc (nstp, sizeof (double));
	has_evt = (int*) calloc (nstp, sizeof (int));

	if (!date || !has_evt)
	{
		err = "Memory allocation error (3) in SrtGrfn3DFXBetaTree";
		goto FREE_RETURN;
	}
	
	for (i=0; i<nstp; i++)
	{		
		date[i] = today + DAYS_IN_YEAR * time[i];

		if (i > 0 && date[i] - date[i-1] >= 1)
		{
			date[i] = (long) (date[i] + 1.0e-08);
			time[i] = YEARS_IN_DAY * (date[i] - today);	
		}		
	}	

	if (time[0] > 0)
	{
		/* add the zero time */
		num_f_add_number(&nstp, &time, 0);
		num_f_sort_vector (nstp, time);	
		nstp -= 1;
		num_f_add_number(&nstp, &date, today);
		num_f_sort_vector (nstp, date);	
		flag = 1;
	}
			

	/* Get all the term structures */
	err = Get_FX_StochRate_TermStructures_corr(und3dfx,
										 &sigma_date_dom,  &sigma_dom,  &sigma_n_dom,
										 &tau_date_dom,  &tau_dom,  &tau_n_dom,
										 &sigma_date_for,  &sigma_for,  &sigma_n_for,
										 &tau_date_for,  &tau_for,  &tau_n_for,
										 &sigma_date_fx,  &sigma_fx,  &sigma_n_fx,									
										 &correl_mat, &correl_dom_for,  &correl_dom_fx,  &correl_for_fx, &nb_correl);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (tau_dom, tau_n_dom, &dom_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (tau_for, tau_n_for, &for_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Get Fx spot and yield curves */	
	dom_yc = get_ycname_from_irund (dom_und);
	for_yc = get_ycname_from_irund (for_und);

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	spot_fx = get_spot_from_fxund (fx_und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);

	/* now merge all the term structure */	
	merge_dates = (double*) calloc (sigma_n_dom, sizeof (double));

	if (!merge_dates)
	{
		err = "Memory allocation error (3) in SrtGrfn3DFXBetaMc";
		goto FREE_RETURN;
	}

	memcpy (merge_dates, sigma_date_dom, sigma_n_dom * sizeof (double));
	nb_merge_dates = sigma_n_dom;
	num_f_concat_vector (&nb_merge_dates, &merge_dates, sigma_n_for, sigma_date_for);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, nb_correl, correl_mat);
	num_f_sort_vector (nb_merge_dates, merge_dates);	
	num_f_unique_vector (&nb_merge_dates, merge_dates);
	
	/*	Fill the new term structures */
	
	sig_dom = (double*) calloc (nb_merge_dates, sizeof (double));
	sig_for = (double*) calloc (nb_merge_dates, sizeof (double));
	sig_fx = (double*) calloc (nb_merge_dates, sizeof (double));

	dom_vol = (double*) calloc (nstp, sizeof (double));
	for_vol = (double*) calloc (nstp, sizeof (double));
	fx_vol = (double*) calloc (nstp, sizeof (double));
	fx_fwd = (double*) calloc (nstp, sizeof (double));

	corr_dom_for = (double*) calloc (nstp, sizeof (double));
	corr_dom_fx = (double*) calloc (nstp, sizeof (double));
	corr_for_fx = (double*) calloc (nstp, sizeof (double));

	dom_ifr = (double*) calloc (nstp, sizeof (double));
	dom_phi= (double*) calloc (nstp, sizeof (double));
	for_ifr = (double*) calloc (nstp, sizeof (double));
	for_phi = (double*) calloc (nstp, sizeof (double));

	if (!sig_dom || !sig_for || !sig_fx || !dom_vol || !for_vol || !fx_vol || !fx_fwd
		|| !corr_dom_for || !corr_dom_fx || !corr_for_fx 
		|| !dom_ifr || !dom_phi || !for_ifr || !for_phi)
	{
		err = "Memory allocation error (4) in SrtGrfn3DFXBetaMc";
		goto FREE_RETURN;
	}
	
	for (i=nb_merge_dates-1; i>=0; i--)
	{
		sig_dom[i] = find_sig (merge_dates[i], dom_ts);
		sig_for[i] = find_sig (merge_dates[i], for_ts);
		sig_fx[i] = find_fx_sig (merge_dates[i], fx_ts);
	}	

	/*	Get phi, std... */
	err = fill_Betamc_init(
							date,
							time,
							nstp,
							merge_dates, 
							nb_merge_dates,
							sig_dom,
							dom_lam,
							sig_for,
							for_lam,
							sig_fx,
							dom_yc,
							for_yc,
							dom_ifr,
							dom_vol,
							dom_phi,					
							for_ifr,					
							for_vol,
							for_phi,
							fx_vol
							);

	if (err)
	{
		goto FREE_RETURN;
	}

	for (i=nstp-1; i>=0; i--)
	{
		corr_dom_for[i] = correl_dom_for[Get_Index(time[i], correl_mat, nb_correl)];
		corr_dom_fx[i] = correl_dom_fx[Get_Index(time[i], correl_mat, nb_correl)];
		corr_for_fx[i] = correl_for_fx[Get_Index(time[i], correl_mat, nb_correl)];
	}

	/* compute the fwd */
	for (i=0; i<nstp; i++)
	{		
		fx_fwd[i] =	  (swp_f_zr (today, date[i], dom_yc)
					- swp_f_zr (today, date[i], for_yc)) * time[i];
		fx_fwd[i] = spot_fx * exp(fx_fwd[i]);
	}

	/* compute the parameters */	
	alpha2 = 1.0 - beta + gamma;	
	beta2 = beta - 2.0 * gamma;
	gamma2 = gamma;
	
	err = allocate_local_model(nstp, &model_func, &model_param);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = fill_local_model(model_func, model_param, nstp, time, alpha2, beta2, gamma2, sig0, fx_fwd);
	if (err)
	{
		goto FREE_RETURN;
	}

	/* Check if there is an event */
	j = 0;
	for (i=0; (i<nstp && j<xStr.num_evt); i++)
	{
		if (time[i] == xStr.tms[j])
		{
			has_evt[i] = 1;
			j += 1;
		}
	}


	/*	Fill product structure */
	strupper (und3dfx);
	strip_white_space (und3dfx);
	strupper (domname);
	strip_white_space (domname);
	strupper (forname);
	strip_white_space (forname);
	for (i=0; i<xStr.num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	/* get the index of each underlying */
	fx_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, und3dfx))
		{
			fx_idx = i;
		}
	}
	if (fx_idx == -1)
	{
		err = "The Fx underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	dom_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, domname))
		{
			dom_idx = i;
		}
	}
	if (dom_idx == -1)
	{
		err = "The domestic underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	for_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, forname))
		{
			for_idx = i;
		}
	}
	if (for_idx == -1)
	{
		err = "The foreign underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}
	
	void_prm = (void**) calloc (nstp, sizeof (void*));

	if (!void_prm)
	{
		err = "Memory allocation error (6) in SrtGrfn3DFXMc";
		goto FREE_RETURN;
	}
	
	for (i=xStr.num_evt-1; i>=0; i--)
	{
		if (xStr.evt[i].evt)
		{
			err = find_phi (date, dom_phi, nstp, xStr.dts[i], &dom_phi_at_t);
			err = find_phi (date, for_phi, nstp, xStr.dts[i], &for_phi_at_t);
			
			if (err)
			{
				goto FREE_RETURN;
			}

			grfn_prm = malloc (sizeof (grfn_parm_mc));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + i;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_fx_df = xStr.evt[i].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[i].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[i].evt->dfd[fx_idx];

			if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
			{
				grfn_prm->fx_dff = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2 = dvector (0, grfn_prm->num_fx_df-1);
				
				if (!grfn_prm->fx_dff || !grfn_prm->fx_gam || !grfn_prm->fx_gam2)
				{
					err = "Memory allocation error (7) in SrtGrfn3DFXBetaMc";
					goto FREE_RETURN;
				}

				for (j=0; j<grfn_prm->num_fx_df; j++)
				{
					grfn_prm->fx_dff[j] = swp_f_df (xStr.dts[i], grfn_prm->fx_df_dts[j], (char*) dom_yc);
					grfn_prm->fx_gam[j] = (1.0 - exp ( - dom_lam * grfn_prm->fx_df_tms[j] )) / dom_lam;
					grfn_prm->fx_gam2[j] = 0.5 * grfn_prm->fx_gam[j] * grfn_prm->fx_gam[j] * dom_phi_at_t;
				}

				grfn_prm->do_fx = 1;
			}
			else
			{
				grfn_prm->do_fx = 0;
			}

			grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

			if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
			{
				grfn_prm->dom_dff = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2 = dvector (0, grfn_prm->num_dom_df-1);
				
				if (!grfn_prm->dom_dff || !grfn_prm->dom_gam || !grfn_prm->dom_gam2)
				{
					err = "Memory allocation error (8) in SrtGrfn3DFXMc";
					goto FREE_RETURN;
				}

				for (j=0; j<grfn_prm->num_dom_df; j++)
				{
					grfn_prm->dom_dff[j] = swp_f_df (xStr.dts[i], grfn_prm->dom_df_dts[j], (char*) dom_yc);
					grfn_prm->dom_gam[j] = (1.0 - exp ( - dom_lam * grfn_prm->dom_df_tms[j] )) / dom_lam;
					grfn_prm->dom_gam2[j] = 0.5 * grfn_prm->dom_gam[j] * grfn_prm->dom_gam[j] * dom_phi_at_t;
				}

				grfn_prm->do_dom = 1;
			}
			else
			{
				grfn_prm->do_dom = 0;
			}

			grfn_prm->num_for_df = xStr.evt[i].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[i].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[i].evt->dfd[for_idx];
			
			if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
			{
				grfn_prm->for_dff = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2 = dvector (0, grfn_prm->num_for_df-1);
				
				if (!grfn_prm->for_dff || !grfn_prm->for_gam || !grfn_prm->for_gam2)
				{
					err = "Memory allocation error (9) in SrtGrfn3DFXMc";
					goto FREE_RETURN;
				}

				for (j=0; j<grfn_prm->num_for_df; j++)
				{
					grfn_prm->for_dff[j] = swp_f_df (xStr.dts[i], grfn_prm->for_df_dts[j], (char*) for_yc);
					grfn_prm->for_gam[j] = (1.0 - exp ( - for_lam * grfn_prm->for_df_tms[j] )) / for_lam;
					grfn_prm->for_gam2[j] = 0.5 * grfn_prm->for_gam[j] * grfn_prm->for_gam[j] * for_phi_at_t;
				}

				grfn_prm->do_for = 1;
			}
			else
			{
				grfn_prm->do_for = 0;
			}
					
			void_prm[i+flag] = (void*) grfn_prm;
		}
		else
		{
			void_prm[i+flag] = NULL;
		}
	}
	
	/*	Eventually! call to function */

	*prod_val = dmatrix (0, num_col-1, 0, 1);

	if (!(*prod_val))
	{
		err = "Memory allocation error";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);

	param[0] = 1.0;
	param[1] = beta;
	param[2] = gamma;

	err = mcLocal_main_3dfx(	
					/*	Time data */
					num_paths,
					nstp,
					num_col,
					time,
					date,
					dom_ifr,		/*	Distributions */
					dom_vol,
					dom_phi,
					for_ifr,
					for_vol,
					for_phi,
					fx_fwd,
					fx_vol,
					model_func,
					model_param,
					void_prm, 
					has_evt,
					/*	Model data */
					dom_lam, 
					for_lam, 
					correl_dom_for,
					correl_dom_fx,
					correl_for_fx,
					/*	Market data */
					spot_fx,
					dom_yc,
					for_yc,
					/* Do PECS adjustment */
					do_pecs,
					/*	Payoff function */
					grfn_payoff_4_3dfx_Quadmc,/*	Result */					
					*prod_val);


	*tableauCols = num_col;

	/*	Add PV of Past */
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;
	
FREE_RETURN:

	if (time) free (time);
	if (date) free (date);
	
	if (dom_ifr) free (dom_ifr);
	if (dom_vol) free (dom_vol);
	if (dom_phi) free (dom_phi);
	
	if (for_ifr) free (for_ifr);
	if (for_vol) free (for_vol);
	if (for_phi) free (for_phi);
	
	if (fx_fwd) free (fx_fwd);
	if (fx_vol) free (fx_vol);

	if (sigma_date_dom)	free (sigma_date_dom);
	if (sigma_dom)	free (sigma_dom);
	if (tau_date_dom)	free (tau_date_dom);
	if (tau_dom)	free (tau_dom);

	if (sigma_date_for)	free (sigma_date_for);
	if (sigma_for)	free (sigma_for);
	if (tau_date_for)	free (tau_date_for);
	if (tau_for)	free (tau_for);
	
	if (sigma_date_fx)	free (sigma_date_fx);
	if (sigma_fx)	free (sigma_fx);

	if (correl_mat) free (correl_mat);
	if (correl_dom_for) free(correl_dom_for);
	if (correl_dom_fx) free(correl_dom_fx);
	if (correl_for_fx) free(correl_for_fx);

	if (corr_dom_for) free(corr_dom_for);
	if (corr_dom_fx) free(corr_dom_fx);
	if (corr_for_fx) free(corr_for_fx);

	if (dom_for_cov) free (dom_for_cov);
	if (dom_fx_cov)	free (dom_fx_cov);
	if (for_fx_cov)	free (for_fx_cov);

	free_local_model_param(model_func, model_param);

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMMC) void_prm[i];
				
				if (	grfn_prm->do_fx 
					&& grfn_prm->num_fx_df > 0 
					&& grfn_prm->fx_idx != -1)
				{
					if (grfn_prm->fx_dff) free_dvector (grfn_prm->fx_dff, 0, grfn_prm->num_fx_df-1);
					if (grfn_prm->fx_gam) free_dvector (grfn_prm->fx_gam, 0, grfn_prm->num_fx_df-1);
					if (grfn_prm->fx_gam2) free_dvector (grfn_prm->fx_gam2, 0, grfn_prm->num_fx_df-1);
				}

				if (	grfn_prm->do_dom
					&& grfn_prm->num_dom_df > 0 
					&& grfn_prm->dom_idx != -1)
				{
					if (grfn_prm->dom_dff) free_dvector (grfn_prm->dom_dff, 0, grfn_prm->num_dom_df-1);
					if (grfn_prm->dom_gam) free_dvector (grfn_prm->dom_gam, 0, grfn_prm->num_dom_df-1);
					if (grfn_prm->dom_gam2) free_dvector (grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df-1);
				}

				if (	grfn_prm->do_for
					&& grfn_prm->num_for_df > 0 
					&& grfn_prm->for_idx != -1)
				{
					if (grfn_prm->for_dff) free_dvector (grfn_prm->for_dff, 0, grfn_prm->num_for_df-1);
					if (grfn_prm->for_gam) free_dvector (grfn_prm->for_gam, 0, grfn_prm->num_for_df-1);
					if (grfn_prm->for_gam2) free_dvector (grfn_prm->for_gam2, 0, grfn_prm->num_for_df-1);
				}

				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (free_str)
	{
		FIRSTFreeMktStruct (&xStr);
	}

	return err;
}



