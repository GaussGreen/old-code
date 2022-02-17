/**********************************************************************
 *      Name: SrtGrfnMainQuanto.c                                     * 
 *  Function: Entry point to GRFN with raw data                       *
 *            in the case of the 3Factor model						  *
 *			  degenerated in Quanto	2 Factor model					  *
 * Copyright: (C) BNP Paribas Capital Markets Ltd.                    *
 *--------------------------------------------------------------------*
 *    Author: Bertrand Baraduc	                                      *
 *      Date: 02/09/01                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 05/01/00 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/

#include "math.h"
#include "srt_h_all.h"
#include "SrtAccess.h"
#include "srt_h_allFx3F.h"
#include "DoubleLGM1FQuanto_pde.h"
#include "DoubleLGM1FQuantoGrfn.h"

Err Get_FX_TS(char *underlying, double **sigma_time, double **sigma, long *sigma_n)
{
	SrtUndPtr und;
	TermStruct *ts;
	long corr_date_n;
	double *corr_date = NULL;
	double *corr = NULL;
	long today;
	int i;
	Err err;

	und = lookup_und (underlying);
	if (!und)
	{
		return serror ("Couldn't find underlying named %s", underlying);
	}
	if (get_underlying_type (und) != FOREX_UND)
	{
		return serror ("Underlying %s is not of type FX", underlying);
	}
	if (get_mdltype_from_fxund (und) != FX_STOCH_RATES)
	{
		return serror ("Underlying %s is not of type FX Stoch Rates", underlying);
	}

	ts = get_ts_from_irund (und);

	today = get_today_from_underlying (und);

	err = srt_f_display_FX_TermStruct(underlying,
										sigma_n,
										sigma_time,
										sigma, 
										&corr_date_n,
										&corr_date,
										&corr);

	for (i = 0; i < *sigma_n; i++)
	{
		(*sigma_time)[i] = ((*sigma_time)[i] - today) / 365.0;
	}

	free (corr_date);
	free (corr);

	return err;
}


void fill_fwd(	
					long		nstp,
					double		*time,
					double		*date,
					char		*dom_yc,
					char		*for_yc,
					/*	To be allocated by caller, filled by the function */
					double		*dom_ifr,
					double		*for_ifr)
{
	int i;

	dom_ifr[0] = swp_f_zr (date[0], date[1], dom_yc);
	for_ifr[0] = swp_f_zr (date[0], date[1], for_yc);
	
	for (i=1; i<nstp ; i++)
	{
		if (i < nstp - 1)
		{
			dom_ifr[i] = swp_f_zr (date[i], date[i+1], dom_yc);
			for_ifr[i] = swp_f_zr (date[i], date[i+1], for_yc);
		}
		else
		{
			dom_ifr[i] = swp_f_zr (date[i], date[i]+1, dom_yc);
			for_ifr[i] = swp_f_zr (date[i], date[i]+1, for_yc);
		}
	}

}


Err fill_time_vector_without_barrier(
					double				**time, 
					int					*nstp, 
					int					num_vol_times, 
					double				*vol_times, 
					int					target_nstp)
{
	Err				err = NULL;

	/*	Add today if required */
	if ((*time)[0] < -EPS)
	{
		err = "Past event date in SrtGrfnQuantoPDE";
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


FREE_RETURN:

	return err;
}


char *SrtGrfnQuantoPDE(
				  char		*und3dfx,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  long		num_stp, 
				  long		num_stpx, 
				  int		*num_prod, 
				  double	**prod_val)
{
	int				free_str = 0;
	FIRSTAllMkts	xStr;
	SrtGrfnParam	defParm;
	int				forback;
	long			nstp;
	
	double			*time			= NULL,
					*date			= NULL;
	
	double			*dom_ifr	= NULL,		
					*for_ifr	= NULL;
					
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
					for_lam;
	double			spot_fx;
	char			*dom_yc,
					*for_yc;
	int				fx_idx,
					dom_idx,
					for_idx;

	int				num_vol_times;
	double			*vol_times = NULL;

	double *domsigtime = NULL;
	double *domsig = NULL;
	long domsig_n;
	double *domtau = NULL;
	double *domtautime = NULL;
	long domtau_n;

	double *forsigtime = NULL;
	double *forsig = NULL;
	long forsig_n;
	double *fortau = NULL;
	double *fortautime = NULL;
	long fortau_n;

	double *fxsigtime = NULL;
	double *fxsig = NULL;
	long fxsig_n;

	double *merge_dates = NULL;

	double *sigtime	= NULL;
	double *sigdom = NULL;
	double *sigfor = NULL;
	double *sigfx = NULL;

	double *sigtime_temp	= NULL;
	double *sigdom_temp		= NULL;
	double *sigfor_temp		= NULL;
	double *sigfx_temp		= NULL;

	double quantocorr;
	double domforcorr;
	double domfxcorr;

	int nb_merge_dates;

	clock_t			t1, t2;

	Err				err = NULL;
	
	t1 = clock();

	//	Initialise the GRFN tableau 

	//	First, initialise the param struct 
	
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

	//	Now, lookup underlyings involved and their term structures 

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

	//	Next, get the time steps 

	//	Copy event dates 
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
		err = "Memory allocation error (1) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	//	Get vol dates 
	err = compute_vol_times (und3dfx, &num_vol_times, &vol_times, time[nstp-1]);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = fill_time_vector_without_barrier(	&time, 
								&nstp, 
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
		err = "Memory allocation error (3) in SrtGrfnQuantoPDE";
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
						  
	//	Fill the model parameters as required by doublelgm1fQuanto_adi2 

	sCorrlist = srt_f_GetTheCorrelationList();
	if (!sCorrlist->head->element)
	{
		err = "correlation list improperly initialised";
		goto FREE_RETURN;
	}


//	 Get all the term structures 
	err = Get_FX_StochRate_TermStructures(und3dfx,
										 &domsigtime,  &domsig,  &domsig_n,
										 &domtautime,  &domtau,  &domtau_n,
										 &forsigtime,  &forsig,  &forsig_n,
										 &fortautime,  &fortau,  &fortau_n,
										 &fxsigtime,  &fxsig,  &fxsig_n,									
										 &domforcorr,  &domfxcorr,  &quantocorr);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (domtau, domtau_n, &dom_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (fortau, fortau_n, &for_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

//	/* now merge all the term structure 
	merge_dates = (double*) calloc (domsig_n, sizeof (double));
	memcpy (merge_dates, domsigtime, domsig_n * sizeof (double));
	nb_merge_dates = domsig_n;
	num_f_concat_vector (&nb_merge_dates, &merge_dates, forsig_n, forsigtime);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, fxsig_n, fxsigtime);
	num_f_sort_vector (nb_merge_dates, merge_dates);	
	num_f_unique_vector (&nb_merge_dates, merge_dates);

//	/*	Fill the new term structures 
	
	sigdom = (double*) calloc (nb_merge_dates, sizeof (double));
	sigfor = (double*) calloc (nb_merge_dates, sizeof (double));
	sigfx = (double*) calloc (nb_merge_dates, sizeof (double));

	if (!sigdom || !sigfor || !sigfx)
	{
		err = "Memory allocation error (4) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}
	
	for (i=nb_merge_dates-1; i>=0; i--)
	{
		sigdom[i] = find_sig (merge_dates[i], dom_ts);
		sigfor[i] = find_sig (merge_dates[i], for_ts);
		sigfx[i] = find_fx_sig (merge_dates[i], fx_ts);
	}


	//	Get Fx spot and yield curves 
	
	dom_yc = get_ycname_from_irund (dom_und);
	for_yc = get_ycname_from_irund (for_und);

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	spot_fx = get_spot_from_fxund (fx_und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);

	//	Get ifr

	dom_ifr = (double*) calloc (nstp, sizeof (double));
	for_ifr = (double*) calloc (nstp, sizeof (double));

	if (!dom_ifr || !for_ifr)
	{
		err = "Memory allocation error (5) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	fill_fwd( nstp, time, date, dom_yc, for_yc, dom_ifr, for_ifr);


	//	Fill product structure 

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

	if (!is_event || !void_prm)
	{
		err = "Memory allocation error (6) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;
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
			
			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else
		if (j>=0 && xStr.am[j])
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
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}

	}

	//	Eventually! call to function 

	*num_prod = xStr.num_cols;
	*prod_val = (double*) calloc (*num_prod+1, sizeof (double));

	if (!(*prod_val))
	{
		err = "Memory allocation error (7) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Number of times steps, required: %d, actual: %d", num_stp, nstp);

	err = doublelgm1fQuanto_adi(
					nstp,
					time,
					date,
					num_stpx,
					dom_lam, 
					for_lam, 
					merge_dates,	
					sigdom,	
					sigfor, 
					sigfx,
					nb_merge_dates,
					quantocorr,
					domforcorr,
					void_prm, 
					is_event,
					dom_ifr,	
					for_ifr,
					dom_yc,
					for_yc,
					payoff_doublelgm1fquanto_pde,
					*num_prod,
					*prod_val);

	//	Add PV of Past 
	(*prod_val)[*num_prod-1] += xStr.gd->pv_of_past;


FREE_RETURN:

	if (time) free (time);

	if (date) free (date);

	if (sigdom) free (sigdom);

	if (sigfor) free (sigfor);

	if (sigfx) free (sigfx);

	if (merge_dates) free (merge_dates);

	if (domsigtime) free (domsigtime);
	if (domsig) free (domsig);

	if (domtautime) free (domtautime);
	if (domtau) free (domtau);

	if (forsigtime) free (forsigtime);
	if (forsig) free (forsig);

	if (fortautime) free (fortautime);
	if (fortau) free (fortau);

	if (fxsigtime) free (fxsigtime);
	if (fxsig) free (fxsig);

	if (dom_ifr) free (dom_ifr);
	if (for_ifr) free (for_ifr);

	if (vol_times) free (vol_times);

	if (is_event) free (is_event);

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


char *SrtGrfnQuantoPDEWithCorrelTS(
				  char		*und3dfx,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  long		num_stp, 
				  long		num_stpx, 
				  int		*num_prod, 
				  double	**prod_val)
{
	int				free_str = 0;
	FIRSTAllMkts	xStr;
	SrtGrfnParam	defParm;
	int				forback;
	long			nstp;
	
	double			*time			= NULL,
					*date			= NULL;
	
	double			*dom_ifr	= NULL,		
					*for_ifr	= NULL;
					
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
					for_lam;
	double			spot_fx;
	char			*dom_yc,
					*for_yc;
	int				fx_idx,
					dom_idx,
					for_idx;

	int				num_vol_times;
	double			*vol_times = NULL;

	double *domsigtime = NULL;
	double *domsig = NULL;
	long domsig_n;
	double *domtau = NULL;
	double *domtautime = NULL;
	long domtau_n;

	double *forsigtime = NULL;
	double *forsig = NULL;
	long forsig_n;
	double *fortau = NULL;
	double *fortautime = NULL;
	long fortau_n;

	double *fxsigtime = NULL;
	double *fxsig = NULL;
	long fxsig_n;

	double *merge_dates = NULL;

	double *sigtime	= NULL;
	double *sigdom = NULL;
	double *sigfor = NULL;
	double *sigfx = NULL;

	double *correltime	= NULL;
	double *quantocorr	= NULL;
	double *domforcorr	= NULL;
	double *domfxcorr	= NULL;
	long correl_n;

	int nb_merge_dates;

	clock_t			t1, t2;

	Err				err = NULL;
	
	t1 = clock();

	//	Initialise the GRFN tableau 

	//	First, initialise the param struct 
	
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

	//	Now, lookup underlyings involved and their term structures 

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

	//	Next, get the time steps 

	//	Copy event dates 
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
		err = "Memory allocation error (1) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));


	//	Get vol dates 
	err = compute_vol_times (und3dfx, &num_vol_times, &vol_times, time[nstp-1]);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = fill_time_vector_without_barrier(	&time, 
								&nstp, 
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
		err = "Memory allocation error (3) in SrtGrfnQuantoPDE";
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
						  
	//	Fill the model parameters as required by doublelgm1fQuanto_adi2 

	sCorrlist = srt_f_GetTheCorrelationList();
	if (!sCorrlist->head->element)
	{
		err = "correlation list improperly initialised";
		goto FREE_RETURN;
	}


//	 Get all the term structures 
	err = Get_FX_StochRate_TermStructures_corr(und3dfx,
										 &domsigtime,  &domsig,  &domsig_n,
										 &domtautime,  &domtau,  &domtau_n,
										 &forsigtime,  &forsig,  &forsig_n,
										 &fortautime,  &fortau,  &fortau_n,
										 &fxsigtime,  &fxsig,  &fxsig_n,
										 &correltime,  
										 &domforcorr,  &domfxcorr,  &quantocorr,
										 &correl_n);

	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (domtau, domtau_n, &dom_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (fortau, fortau_n, &for_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

//	/* now merge all the term structure 
	merge_dates = (double*) calloc (domsig_n, sizeof (double));
	memcpy (merge_dates, domsigtime, domsig_n * sizeof (double));
	nb_merge_dates = domsig_n;
	num_f_concat_vector (&nb_merge_dates, &merge_dates, forsig_n, forsigtime);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, fxsig_n, fxsigtime);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, correl_n, correltime);
	num_f_sort_vector (nb_merge_dates, merge_dates);	
	num_f_unique_vector (&nb_merge_dates, merge_dates);

//	/*	Fill the new term structures 
	
	sigdom = (double*) calloc (nb_merge_dates, sizeof (double));
	sigfor = (double*) calloc (nb_merge_dates, sizeof (double));
	sigfx = (double*) calloc (nb_merge_dates, sizeof (double));
	domforcorr = (double*) calloc (nb_merge_dates, sizeof (double));
	quantocorr = (double*) calloc (nb_merge_dates, sizeof (double));

	if (!sigdom || !sigfor || !sigfx || !domforcorr || !quantocorr)
	{
		err = "Memory allocation error (4) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}
	
	for (i=nb_merge_dates-1; i>=0; i--)
	{
		sigdom[i] = find_sig (merge_dates[i], dom_ts);
		sigfor[i] = find_sig (merge_dates[i], for_ts);
		sigfx[i] = find_fx_sig (merge_dates[i], fx_ts);
		err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname, merge_dates[i], &(domforcorr[i]));
		err = srt_f_get_corr_from_CorrList(sCorrlist, forname, und3dfx, merge_dates[i], &(quantocorr[i]));
	}

	//	Get Fx spot and yield curves 
	
	dom_yc = get_ycname_from_irund (dom_und);
	for_yc = get_ycname_from_irund (for_und);

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	spot_fx = get_spot_from_fxund (fx_und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);

	//	Get ifr

	dom_ifr = (double*) calloc (nstp, sizeof (double));
	for_ifr = (double*) calloc (nstp, sizeof (double));

	if (!dom_ifr || !for_ifr)
	{
		err = "Memory allocation error (5) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	fill_fwd( nstp, time, date, dom_yc, for_yc, dom_ifr, for_ifr);


	//	Fill product structure 

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

	if (!is_event || !void_prm)
	{
		err = "Memory allocation error (6) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;
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
			
			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else
		if (j>=0 && xStr.am[j])
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
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}

	}

	//	Eventually! call to function 

	*num_prod = xStr.num_cols;
	*prod_val = (double*) calloc (*num_prod+1, sizeof (double));

	if (!(*prod_val))
	{
		err = "Memory allocation error (7) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Number of times steps, required: %d, actual: %d", num_stp, nstp);

	err = doublelgm1fQuanto_adi_correl(
					nstp,
					time,
					date,
					num_stpx,
					dom_lam, 
					for_lam, 
					merge_dates,	
					sigdom,	
					sigfor, 
					sigfx,
					nb_merge_dates,
					quantocorr,
					domforcorr,
					void_prm, 
					is_event,
					dom_ifr,	
					for_ifr,
					dom_yc,
					for_yc,
					payoff_doublelgm1fquanto_pde,
					*num_prod,
					*prod_val);

	//	Add PV of Past 
	(*prod_val)[*num_prod-1] += xStr.gd->pv_of_past;


FREE_RETURN:

	if (time) free (time);

	if (date) free (date);

	if (sigdom) free (sigdom);

	if (sigfor) free (sigfor);

	if (sigfx) free (sigfx);

	if (merge_dates) free (merge_dates);

	if (domsigtime) free (domsigtime);
	if (domsig) free (domsig);

	if (domtautime) free (domtautime);
	if (domtau) free (domtau);

	if (forsigtime) free (forsigtime);
	if (forsig) free (forsig);

	if (fortautime) free (fortautime);
	if (fortau) free (fortau);

	if (fxsigtime) free (fxsigtime);
	if (fxsig) free (fxsig);

	if (dom_ifr) free (dom_ifr);
	if (for_ifr) free (for_ifr);

	if (vol_times) free (vol_times);

	if (is_event) free (is_event);

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



char *SrtGrfnQuantoPDEWithCorrelTS2(
				  char		*und3dfx,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  long		num_stp, 
				  long		num_stpx, 
				  int		*num_prod, 
				  double	**prod_val)
{
	int				free_str = 0;
	FIRSTAllMkts	xStr;
	SrtGrfnParam	defParm;
	int				forback;
	long			nstp;
	long			nb_time;
	
	double			*time			= NULL,
					*date			= NULL;
	
	double			*dom_ifr	= NULL,		
					*for_ifr	= NULL;
					
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
					for_lam;
	double			spot_fx;
	char			*dom_yc,
					*for_yc;
	int				fx_idx,
					dom_idx,
					for_idx;

	int				num_vol_times;
	double			*vol_times = NULL;

	double *domsigtime = NULL;
	double *domsig = NULL;
	long domsig_n;
	double *domtau = NULL;
	double *domtautime = NULL;
	long domtau_n;

	double *forsigtime = NULL;
	double *forsig = NULL;
	long forsig_n;
	double *fortau = NULL;
	double *fortautime = NULL;
	long fortau_n;

	double *fxsigtime = NULL;
	double *fxsig = NULL;
	long fxsig_n;

	double *merge_dates = NULL;

	double *sigtime	= NULL;
	double *sigdom = NULL;
	double *sigfor = NULL;
	double *sigfx = NULL;

	double *correltime	= NULL;
	double *quantocorr	= NULL;
	double *domforcorr	= NULL;
	double *domfxcorr	= NULL;
	long correl_n;

	int nb_merge_dates;

	double last_time;

	clock_t			t1, t2;

	Err				err = NULL;
	
	t1 = clock();

	//	Initialise the GRFN tableau 

	//	First, initialise the param struct 
	
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

	//	Now, lookup underlyings involved and their term structures 

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

	//	Next, get the time steps 

	//	Copy event dates 
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
		err = "Memory allocation error (1) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));


	//	Get vol dates 
	err = compute_vol_times (und3dfx, &num_vol_times, &vol_times, time[nstp-1]);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = fill_time_vector_without_barrier(	&time, 
								&nstp, 
								num_vol_times, 
								vol_times, 
								num_stp);
	if (err)
	{
		goto FREE_RETURN;
	}

	//	Fill the model parameters as required by doublelgm1fQuanto_adi2 

	sCorrlist = srt_f_GetTheCorrelationList();
	if (!sCorrlist->head->element)
	{
		err = "correlation list improperly initialised";
		goto FREE_RETURN;
	}


//	 Get all the term structures 
	err = Get_FX_StochRate_TermStructures_corr(und3dfx,
										 &domsigtime,  &domsig,  &domsig_n,
										 &domtautime,  &domtau,  &domtau_n,
										 &forsigtime,  &forsig,  &forsig_n,
										 &fortautime,  &fortau,  &fortau_n,
										 &fxsigtime,  &fxsig,  &fxsig_n,
										 &correltime,  
										 &domforcorr,  &domfxcorr,  &quantocorr,
										 &correl_n);

	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (domtau, domtau_n, &dom_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = get_unique_lambda (fortau, fortau_n, &for_lam);
	if (err)
	{
		goto FREE_RETURN;
	}

//	/* now merge all the term structure 
	merge_dates = (double*) calloc (domsig_n, sizeof (double));
	memcpy (merge_dates, domsigtime, domsig_n * sizeof (double));
	nb_merge_dates = domsig_n;
	num_f_concat_vector (&nb_merge_dates, &merge_dates, forsig_n, forsigtime);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, fxsig_n, fxsigtime);
	num_f_concat_vector (&nb_merge_dates, &merge_dates, correl_n, correltime);
	num_f_sort_vector (nb_merge_dates, merge_dates);	
	num_f_unique_vector (&nb_merge_dates, merge_dates);

//	/*	Fill the new term structures 
	
	sigdom = (double*) calloc (nb_merge_dates, sizeof (double));
	sigfor = (double*) calloc (nb_merge_dates, sizeof (double));
	sigfx = (double*) calloc (nb_merge_dates, sizeof (double));
	domforcorr = (double*) calloc (nb_merge_dates, sizeof (double));
	quantocorr = (double*) calloc (nb_merge_dates, sizeof (double));

	if (!sigdom || !sigfor || !sigfx || !domforcorr || !quantocorr)
	{
		err = "Memory allocation error (4) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}
	
	for (i=nb_merge_dates-1; i>=0; i--)
	{
		sigdom[i] = find_sig (merge_dates[i], dom_ts);
		sigfor[i] = find_sig (merge_dates[i], for_ts);
		sigfx[i] = find_fx_sig (merge_dates[i], fx_ts);
		err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname, merge_dates[i], &(domforcorr[i]));
		err = srt_f_get_corr_from_CorrList(sCorrlist, forname, und3dfx, merge_dates[i], &(quantocorr[i]));
	}


//----------------------------------------------------------------------
//---------------Merge Vol times and discretisation times---------------
//----------------------------------------------------------------------
	nb_time = nstp;
	last_time = time[nstp-1];
	num_f_concat_vector (&nb_time, &time, nb_merge_dates, merge_dates);
	num_f_sort_vector (nb_time, time);	
	num_f_unique_vector (&nb_time, time);

	i=0;
	while(time[i] != last_time)
	{
		i++;
	}
	nstp = i+1;


	date = (double*) calloc (nstp, sizeof (double));
	if (!date)
	{
		err = "Memory allocation error (3) in SrtGrfnQuantoPDE";
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


	//	Get Fx spot and yield curves 
	
	dom_yc = get_ycname_from_irund (dom_und);
	for_yc = get_ycname_from_irund (for_und);

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	spot_fx = get_spot_from_fxund (fx_und) 
		* swp_f_df (today, spot_date, dom_yc)
		/ swp_f_df (today, spot_date, for_yc);

	//	Get ifr

	dom_ifr = (double*) calloc (nstp, sizeof (double));
	for_ifr = (double*) calloc (nstp, sizeof (double));

	if (!dom_ifr || !for_ifr)
	{
		err = "Memory allocation error (5) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	fill_fwd( nstp, time, date, dom_yc, for_yc, dom_ifr, for_ifr);


	//	Fill product structure 

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

	if (!is_event || !void_prm)
	{
		err = "Memory allocation error (6) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;
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
			
			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else
		if (j>=0 && xStr.am[j])
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
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}

	}

	//	Eventually! call to function 

	*num_prod = xStr.num_cols;
	*prod_val = (double*) calloc (*num_prod+1, sizeof (double));

	if (!(*prod_val))
	{
		err = "Memory allocation error (7) in SrtGrfnQuantoPDE";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Number of times steps, required: %d, actual: %d", num_stp, nstp);

	err = doublelgm1fQuanto_adi_correl2(
					nstp,
					time,
					date,
					num_stpx,
					dom_lam, 
					for_lam, 
					merge_dates,	
					sigdom,	
					sigfor, 
					sigfx,
					nb_merge_dates,
					quantocorr,
					domforcorr,
					void_prm, 
					is_event,
					dom_ifr,	
					for_ifr,
					dom_yc,
					for_yc,
					payoff_doublelgm1fquanto_pde,
					*num_prod,
					*prod_val);

	//	Add PV of Past 
	(*prod_val)[*num_prod-1] += xStr.gd->pv_of_past;


FREE_RETURN:

	if (time) free (time);

	if (date) free (date);

	if (sigdom) free (sigdom);

	if (sigfor) free (sigfor);

	if (sigfx) free (sigfx);

	if (merge_dates) free (merge_dates);

	if (domsigtime) free (domsigtime);
	if (domsig) free (domsig);

	if (domtautime) free (domtautime);
	if (domtau) free (domtau);

	if (forsigtime) free (forsigtime);
	if (forsig) free (forsig);

	if (fortautime) free (fortautime);
	if (fortau) free (fortau);

	if (fxsigtime) free (fxsigtime);
	if (fxsig) free (fxsig);

	if (dom_ifr) free (dom_ifr);
	if (for_ifr) free (for_ifr);

	if (vol_times) free (vol_times);

	if (is_event) free (is_event);

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

