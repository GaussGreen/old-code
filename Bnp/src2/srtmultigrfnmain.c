/**********************************************************************
 *      Name: SrtGrfnMain.c                                           * 
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 18/10/95                                                *
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
#include "srt_h_allFx3F.h"
#include "SrtAccess.h"
#include "BGMEval.h"

/* returns the first index whose correspinding maturity is greater or equal to the maturity T */
static long Get_Index(double T, double *Maturity, long nbrMat)
{
static int i;

	for (i=0 ; (i < nbrMat) && (Maturity[i] < T) ; i++);

	if (i >= nbrMat)
		i = i-1;

	return i;
}

static long Get_Index_Und(
						  SrtUndPtr	*und_ptr,
						  int		num_und,
						  char*		und_name)
{
static int i;
	
	for (i=0; i<num_und; i++)
	{
		if (!strcmp (und_ptr[i] -> underl_name, und_name))
			return i;
	}

	return -1;
}


static long Get_Index_FxUnd(
						  SrtUndPtr	*und_ptr,
						  int		num_und,
						  char*		domname,
						  char*		forname)
{
static int i;
static char		*temp_domname,
				*temp_forname;	

	for (i=0; i<num_und; i++)
	{
		if (get_mdltype_from_fxund (und_ptr[i]) == FX_STOCH_RATES)
		{
			temp_domname = get_domname_from_fxund (und_ptr[i]);
			temp_forname = get_forname_from_fxund (und_ptr[i]);

			if ((!strcmp (temp_domname, domname)) && (!strcmp (temp_forname, forname)))
				return i;
		}
	}

	return -1;
}

static Err  CompleteUndPtrList(
								  SrtUndPtr **und_ptr,
								  int		*num_und,
								  char		*domestic_name,
								  char		*domestic_ccy,
								  int		**type,
								  int		**dom_forex,
								  int		**for_forex,
								  int		**fx_index,
								  int		*dom_index
								  )

{
int					i, index, num_und_temp;
char				*domname, *forname;
SrtUndListPtr		all_und;
Err					err = NULL;

	all_und = get_underlying_list();
	num_und_temp = *num_und;	

	for (i=0; i<num_und_temp; i++)
	{
		switch get_mdltype_from_fxund (((*und_ptr)[i]))			
		{
			case FX_STOCH_RATES:

				/*	we have to check that the corresponding domestic and foreign
					are in the structure
				*/
				(*type)[i] = 2;
				domname = get_domname_from_fxund (((*und_ptr)[i]));
				forname = get_forname_from_fxund (((*und_ptr)[i]));

				index = Get_Index_Und (*und_ptr, num_und_temp, domname);

				if (index < 0)
				{
					/* we have to add this underlying */
					num_und_temp += 1;
					*und_ptr = (SrtUndPtr*) realloc (*und_ptr, num_und_temp * sizeof(SrtUndPtr));
					*type = (int *) realloc (*type, num_und_temp * sizeof(int));
					*dom_forex = (int *) realloc (*dom_forex, num_und_temp * sizeof(int));
					*for_forex = (int *) realloc (*for_forex, num_und_temp * sizeof(int));
					*fx_index  = (int *) realloc (*fx_index, num_und_temp * sizeof(int));

					if (!*und_ptr || !*type || !*dom_forex || !*for_forex || !*fx_index)
					{
						err = "Memory reallocation failure";
						goto FREE_RETURN;
					}

					err =  srt_f_getundinlist (all_und, domname, &((*und_ptr)[num_und_temp-1])); 

					if (err)
					{
						goto FREE_RETURN;
					}
					index = num_und_temp-1;
				}
				(*dom_forex)[i] = index;

				index = Get_Index_Und (*und_ptr, num_und_temp, forname);
				if (index < 0)
				{
					/* we have to add this underlying */
					num_und_temp += 1;
					*und_ptr = (SrtUndPtr*) realloc (*und_ptr, num_und_temp * sizeof(SrtUndPtr));
					*type = (int *) realloc (*type, num_und_temp * sizeof(int));
					*dom_forex = (int *) realloc (*dom_forex, num_und_temp * sizeof(int));
					*for_forex = (int *) realloc (*for_forex, num_und_temp * sizeof(int));
					*fx_index  = (int *) realloc (*fx_index, num_und_temp * sizeof(int));

					if (!*und_ptr || !*type || !*dom_forex || !*for_forex || !*fx_index)
					{
						err = "Memory reallocation failure";
						goto FREE_RETURN;
					}

					err =  srt_f_getundinlist (all_und, forname, &((*und_ptr)[num_und_temp-1])); 

					if (err)
					{
						goto FREE_RETURN;
					}
					index = num_und_temp-1;
				}
				(*for_forex)[i] = index;

				break;
			
			case  LGM:

				if (!strcmp (((*und_ptr)[i]) -> underl_name, domestic_name))
				{
					*dom_index = i;
					(*type)[i] = 0;
					break;
				}

				if (!strcmp (((*und_ptr)[i]) -> underl_ccy, domestic_ccy))
				{
					/* This in another domestic LGM */
					(*type)[i] = 0;
					break;
				}

				/*	it is a foreign LGM */
				/*	we have to check that the corresponding fxd rate is in the structure */
				
				(*type)[i] = 1;
				index = Get_Index_FxUnd (*und_ptr, num_und_temp, domestic_name, ((*und_ptr)[i]) ->underl_name);
				if (index < 0)
				{
					/* we have to add this underlying */
					num_und_temp += 1;
					*und_ptr = (SrtUndPtr*) realloc (*und_ptr, num_und_temp * sizeof(SrtUndPtr));
					*type = (int *) realloc (*type, num_und_temp * sizeof(int));
					*dom_forex = (int *) realloc (*dom_forex, num_und_temp * sizeof(int));
					*for_forex = (int *) realloc (*for_forex, num_und_temp * sizeof(int));
					*fx_index  = (int *) realloc (*fx_index, num_und_temp * sizeof(int));

					if (!*und_ptr || !*type || !*dom_forex || !*for_forex || !*fx_index)
					{
						err = "Memory reallocation failure";
						goto FREE_RETURN;
					}

					(*und_ptr)[num_und_temp-1] =  srt_f_lookup_fxund (domestic_ccy, ((*und_ptr)[i]) -> underl_ccy); 
					if (!((*und_ptr)[num_und_temp-1]))
					{
						goto FREE_RETURN;
					}
					index = num_und_temp-1;
				}
				(*fx_index)[i] = index;
				break;

			default:

				err = "Only LGM and FX_STOCH_RATE are allowed";	
		}			
	}
FREE_RETURN:

	*num_und = num_und_temp;
	return err;
}

char *SrtMultiGrfn3DFXMc(
				  char		*underlying,
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
				  int		do_pecs,
				  double	***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	int					forback;
	int					flag = 0;
	long				nstp;
						
	double				*evt_tms = NULL,
						*sig_dates = NULL,
						**und_sig_dates = NULL,
						**und_sig_curve = NULL,
						**sig_curves = NULL,
						pay_time;

	int					*type = NULL,
						*dom_forex = NULL,
						*for_forex = NULL,
						*fx_index = NULL,
						dom_index;

	double				*corr_date = NULL,
						*tau_date = NULL,
						*tau = NULL,
						tau_temp,
						df;
						
	long				nb_tau,
						*nb_und_sig = NULL,
						nb_corr_date,
						nb_sig_dates,
						pay_date;

	double				*corr = NULL;

	long				*evt_dts = NULL,
						today, spot_date;
						
	int					num_und,
						num_col,
						max_num_df,
						num_evt;
						
	char				*domestic_name,
						*domestic_ccy;

	SrtUndPtr			*und_ptr = NULL,
						dom_und,						
						und;

	void				**void_prm = NULL;
	GRFNPARM_MULTIMC	grfn_prm;
	LINK_UND			und_link = NULL;
	SrtCorrLstPtr		corr_list = NULL;	 
	long				i, j, k, l;
	clock_t				t1, t2;

	int					sig_dates_alloc = 0;

	Err					err = NULL,
						err2 = NULL;

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

	/* look for the today date */
	today = get_today_from_underlying (und_ptr[0]);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	/* look for the main domestic name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domestic_name = get_domname_from_fxund (und);
	}
	else
	{		
		domestic_name = und -> underl_name;
	}

	domestic_ccy = get_underlying_ccy(und);

	type = (int *) calloc(num_und, sizeof(int));
	dom_forex = (int *) calloc(num_und, sizeof(int));
	for_forex = (int *) calloc(num_und, sizeof(int));
	fx_index  = (int *) calloc(num_und, sizeof(int));

	if (!type || !dom_forex || !for_forex || !fx_index)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	/* Important complete the underlying list */

	err = CompleteUndPtrList(
							  &und_ptr,
							  &num_und,
							  domestic_name,
							  domestic_ccy,
							  &type,
							  &dom_forex,
							  &for_forex,
							  &fx_index,
							  &dom_index
							  );
	
	if (err)
	{
		goto FREE_RETURN;
	}
	
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
	nstp = num_evt;
	
	/* add today to the event dates */
	if (evt_tms[0] > 0)
	{
		/* add the zero time */
		num_f_add_number (&nstp, &evt_tms, 0);
		num_f_sort_vector (nstp, evt_tms);	
		nstp -= 1;
		
		evt_dts = (long*) realloc (evt_dts, (nstp + 1) * sizeof (double));
		
		if (!evt_dts)
		{
			err = "Memory reallocation failure";
			goto FREE_RETURN;
		}
				
		for (i=nstp; i>=1; i--)
		{
			evt_dts[i] = evt_dts[i-1];
		}
		evt_dts[0] = today;
		nstp += 1;
		flag = 1;
	}
		
	/* pay_date is the last event_date */
	pay_date = (long) (evt_dts[nstp-1] + 1.0E-8); 
	pay_time = evt_tms[nstp-1];
	
	/* Fill the link and grfn_prm structure */

	und_link = malloc (sizeof (link_und));
	if (!und_link)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	und_link -> num_und = num_und;
	und_link -> nb_dates = nstp;

	err = allocate_link_und (
							und_link
							);
	if (err)
	{ 
		goto FREE_RETURN;
	}

	und_link -> type = type;
	und_link -> dom_forex = dom_forex;
	und_link -> for_forex = for_forex;
	und_link -> fx_index = fx_index;
	und_link -> dom_index = dom_index;

	type = NULL;
	dom_forex = NULL;
	for_forex = NULL;
	fx_index = NULL;

	/* clean all the underlyings name */
	for (i=0; i<num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	/* get the correlation list */
	corr_list = srt_f_GetTheCorrelationList();

	if (!corr_list)
	{
		err = "Cannot get the correlation structure";
		goto FREE_RETURN;
	}

	/* fill the link structure */
	und_sig_curve = (double **) calloc (num_und, sizeof(double *));
	und_sig_dates = (double **) calloc (num_und, sizeof(double *));
	nb_und_sig = (long *) calloc (num_und, sizeof(long));

	if (!und_sig_curve || !und_sig_dates || !nb_und_sig)
	{
		err = "Memory allocation failure in SrtMultiGrfn3DFXMc";
		goto FREE_RETURN;
	}

	sig_dates_alloc = 0;
	for (i=0; i<num_und; i++)
	{		
		/* first do the LGN Underlyings */
		if (und_link -> type[i] < 2)		
		{
			dom_und = lookup_und (und_ptr[i] -> underl_name);
			if (!dom_und)
			{
				err = serror ("Couldn't find underlying named %s", und_ptr[i] -> underl_name);
				goto FREE_RETURN;
			}
			und_link -> yc[i] = (char *) get_ycname_from_irund (dom_und);
			err = Get_LGM_TermStructure (und_ptr[i] -> underl_name,
										&(und_sig_dates[i]), &(und_sig_curve[i]), &(nb_und_sig[i]), &tau_date, &tau, &nb_tau);
			if (err)
			{				
				goto FREE_RETURN;
			}
			
			if (nb_tau > 1)
			{
				tau_temp = tau[0];
				for (j=1; j<nb_tau; j++)
				{
					if (tau[j] != tau_temp)
					{
						err = "No Tau term structure is allowed on LGM";
						goto FREE_RETURN;
					}
				}					
			}
			und_link -> lambda[i] = 1.0 / tau[0];

			if (!sig_dates_alloc)
			{
				nb_sig_dates = nb_und_sig[i];
				sig_dates = (double*) calloc (nb_sig_dates, sizeof (double));
				if (!sig_dates)
				{
					err = "Memory allocation failure";
					goto FREE_RETURN;
				}
				
				memcpy (sig_dates, und_sig_dates[i], nb_sig_dates * sizeof (double));

				sig_dates_alloc = 1;
			}
			else
			{
				num_f_concat_vector (&nb_sig_dates, &sig_dates, nb_und_sig[i], und_sig_dates[i]);
			}

			if (tau_date)
			{
				free (tau_date);
				tau_date = NULL;
			}
			if (tau)
			{
				free (tau);
				tau = NULL;
			}
		}
	}

	for (i=0; i<num_und; i++)
	{		
		if (und_link -> type[i]  == 2)
		{						
			err = srt_f_display_FX_TermStruct (und_ptr[i] -> underl_name, &(nb_und_sig[i]),
				&(und_sig_dates[i]), &(und_sig_curve[i]), &nb_corr_date, &corr_date, &corr);
			if (err)
			{				
				goto FREE_RETURN;
			}
			
			for (l=0; l<nb_und_sig[i]; l++)
			{
				und_sig_dates[i][l] = (und_sig_dates[i][l] - today) / 365.0;
			}

			und_link -> spot_fx[i] = get_spot_from_fxund (und_ptr[i]) 
						* swp_f_df (today, spot_date, und_link -> yc[und_link -> dom_forex[i]])
						/ swp_f_df (today, spot_date, und_link -> yc[und_link -> for_forex[i]]);

			und_link -> lambda[i] = und_link -> lambda[und_link -> dom_forex[i]];

			for (l=0; l < nb_corr_date; l++) corr_date[l] = (corr_date[l] - today) / 365.0;

			num_f_concat_vector (&nb_sig_dates, &sig_dates, nb_corr_date, corr_date);
			num_f_concat_vector (&nb_sig_dates, &sig_dates, nb_und_sig[i], und_sig_dates[i]);

			if (corr_date) 
			{
				free (corr_date);
				corr_date = NULL;
			}
			if (corr) 
			{
				free (corr);			
				corr = NULL;
			}
		}				
	}

	num_f_sort_vector (nb_sig_dates, sig_dates);
	num_f_unique_vector (&nb_sig_dates, sig_dates);

	und_link->correlations = f3tensor(0, num_und-1, 0, num_und-1, 0, nb_sig_dates-1);

	/* get the correlations at each date */
	for (l=0; l < nb_sig_dates; l++)
	{
		for (i=0; i<num_und; i++)
		{		
			und_link -> correlations[i][i][l] = 1.0;
			for (j=i+1; j<num_und; j++)
			{		
				err = srt_f_get_corr_from_CorrList(corr_list, und_ptr[i] -> underl_name,
					und_ptr[j] -> underl_name, sig_dates[l], &(und_link -> correlations[i][j][l]));
				if (err)
				{
					goto FREE_RETURN;
				}
				und_link -> correlations[j][i][l] = und_link -> correlations[i][j][l];
			}		
		}
	}

	/* fill all the term structures */
	sig_curves = dmatrix(0, num_und-1, 0, nb_sig_dates-1);

	for (i=0; i<num_und; i++)
	{
		for (j=0; j<nb_sig_dates; j++)
		{
			sig_curves[i][j] = und_sig_curve[i][Get_Index(sig_dates[j], und_sig_dates[i], nb_und_sig[i])];
		}
		if (und_sig_dates[i])
		{
			free (und_sig_dates[i]);
			und_sig_dates[i] = NULL;
		}
		if (und_sig_curve[i])
		{
			free (und_sig_curve[i]);
			und_sig_curve[i] = NULL;
		}
	}

	if (und_sig_curve)
	{
		free (und_sig_curve);
		und_sig_curve = NULL;
	}
	if (und_sig_dates)
	{
		free (und_sig_dates);
		und_sig_dates = NULL;
	}
	if (nb_und_sig)
	{
		free (nb_und_sig);
		nb_und_sig = NULL;
	}

	und_link -> nb_sig_dates = nb_sig_dates;
	und_link -> sig_dates = sig_dates;
	und_link -> sig_curve = sig_curves;
	
	sig_dates = NULL;
	sig_curves = NULL;

	/* now calculate all the fwd, std, covariances of the underlyings */	
	err = fill_mc_multi_init_corr(
							pay_date,
							pay_time,
							evt_dts,
							evt_tms,
							nstp,
							und_link
							);

	if (err)
	{
		goto FREE_RETURN;
	}

	/* fill the void_parm structure */
	void_prm = (void**) calloc (nstp, sizeof (void*));

	if (!void_prm)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	for (i=num_evt-1; i>=0; i--)
	{		
		if (xStr.evt[i].evt)
		{
			grfn_prm = malloc (sizeof (grfn_parm_multi_mc));
			if (!grfn_prm)
			{
				err = "Memory allocation failure";
				goto FREE_RETURN;
			}

			grfn_prm -> num_df = lvector (0, num_und-1);
			grfn_prm -> df_tms = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> df_dts = (long **) calloc (num_und, sizeof(long *));

			grfn_prm -> dff = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> gam = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> gam2 = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> do_und = (int *) calloc (num_und, sizeof(int));

			if (!grfn_prm -> num_df || !grfn_prm -> df_tms || !grfn_prm -> df_dts
				|| !grfn_prm -> dff || !grfn_prm -> gam || !grfn_prm -> gam2
				|| !grfn_prm -> do_und)
			{
				err = "Memory allocation failure";
				goto FREE_RETURN;
			}

			grfn_prm -> global = &xStr;
			grfn_prm -> local = xStr.evt + i;
			grfn_prm -> link = und_link;

			for (j=0; j<num_und; j++)
			{
				grfn_prm->num_df[j] = xStr.evt[i].evt->dflen[j];
				grfn_prm->df_tms[j] = xStr.evt[i].evt->dft[j];
				grfn_prm->df_dts[j] = xStr.evt[i].evt->dfd[j];
			
				if (grfn_prm->num_df[j] > 0)
				{
					grfn_prm -> dff[j] = dvector (0, grfn_prm -> num_df[j] - 1);
					grfn_prm -> gam[j] = dvector (0, grfn_prm -> num_df[j] - 1);
					grfn_prm -> gam2[j] = dvector (0, grfn_prm -> num_df[j] - 1);
				
					if (!grfn_prm -> dff[j] || !grfn_prm -> gam[j] || !grfn_prm -> gam2[j])
					{
						err = "Memory allocation failure";
						goto FREE_RETURN;
					}
		
					if (und_link -> type[j] < 2)
					{
						for (k=0; k<grfn_prm->num_df[j]; k++)
						{
							grfn_prm -> dff[j][k] = swp_f_df (evt_dts[i+flag], grfn_prm->df_dts[j][k], und_link -> yc[j]);
							grfn_prm -> gam[j][k] = (1.0 - exp ( - und_link -> lambda[j] * grfn_prm->df_tms[j][k] )) / und_link -> lambda[j];
							grfn_prm -> gam2[j][k] = 0.5 * grfn_prm->gam[j][k] * grfn_prm->gam[j][k] * und_link -> phi[j][i+flag];
						}
					}
					else
					{
						/* in the fx case we are calculating domestic df */
						for (k=0; k<grfn_prm->num_df[j]; k++)
						{
							grfn_prm -> dff[j][k] = swp_f_df (evt_dts[i+flag], grfn_prm->df_dts[j][k], und_link -> yc[j]);
							grfn_prm -> gam[j][k] = (1.0 - exp ( - und_link -> lambda[und_link -> dom_forex[j]] * grfn_prm->df_tms[j][k] )) / und_link -> lambda[und_link -> dom_forex[j]];
							grfn_prm -> gam2[j][k] = 0.5 * grfn_prm->gam[und_link -> dom_forex[j]][k] * grfn_prm->gam[und_link -> dom_forex[j]][k] * und_link -> phi[und_link -> dom_forex[j]][i+flag];
						}
					}

					grfn_prm->do_und[j] = 1;
				}
				else
				{
					grfn_prm->do_und[j] = 0;
				}
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

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);

	err = mc_main_multi_3dfx(	
							/*	Time data */
							num_paths,
							num_col,
							evt_tms,
							evt_dts,
							nstp,
							do_pecs,
							0,
							NULL,
							NULL,
							und_link,
							void_prm, 
							grfn_payoff_multi_3dfx_mc,
							*prod_val);

	if (err)
	{
		goto FREE_RETURN;
	}

	df = swp_f_zr (today, pay_date, und_link -> yc[dom_index]);
	df = exp (-df * pay_time);

	*tableauCols = num_col;
	for (i=0; i<num_col; i++)
	{
		(*prod_val)[i][0] *= df;
		(*prod_val)[i][1] *= df;
	}

	/*	Add PV of Past */
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;

		
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

	if (corr)
	{ 
		free (corr);
	}

	if (corr_date)
	{
		free (corr_date);
	}

	if (tau_date)
	{
		free (tau_date);
	}

	if (tau)
	{
		free (tau);
	}

	if (void_prm)
	{
		for (i=0; i<num_evt; i++)
		{
			if (void_prm[i+flag])
			{
				grfn_prm = (GRFNPARM_MULTIMC) void_prm[i+flag];
				
				for (l=0; l<num_und; l++)
				{
					if (	grfn_prm->do_und[l] 
						&& grfn_prm->num_df[l] > 0 )
						
					{
						if (grfn_prm->dff[l]) free_dvector (grfn_prm->dff[l], 0, grfn_prm->num_df[l]-1);
						if (grfn_prm->gam[l]) free_dvector (grfn_prm->gam[l], 0, grfn_prm->num_df[l]-1);
						if (grfn_prm->gam2[l]) free_dvector (grfn_prm->gam2[l], 0, grfn_prm->num_df[l]-1);
					}								
				}
				if (grfn_prm->dff)
					free (grfn_prm->dff);
				if (grfn_prm->gam)
					free (grfn_prm->gam);
				if (grfn_prm->gam2)
					free (grfn_prm->gam2);
				if (grfn_prm->do_und)
					free (grfn_prm->do_und);

				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (sig_dates)
	{
		free (sig_dates);
	}

	if (sig_curves)
	{
		free_dmatrix (sig_curves, 0, num_und - 1, 0, nb_sig_dates - 1);
	}

	if (type)
	{
		free (type);
	}

	if (dom_forex)
	{
		free (dom_forex);
	}
	
	if (for_forex)
	{
		free (for_forex);
	}
	
	if (fx_index)
	{
		free (fx_index);
	}

	if (und_sig_dates)
	{
		for (i=0; i<num_und; i++)
		{
			if (und_sig_dates[i])
			{
				free (und_sig_dates[i]);
			}
		}
		free (und_sig_dates);
	}

	if (und_sig_curve)
	{
		for (i=0; i<num_und; i++)
		{
			if (und_sig_curve[i])
			{
				free (und_sig_curve[i]);
			}
		}
		free (und_sig_curve);
	}

	if (nb_und_sig)
	{
		free (nb_und_sig);
	}

	free_link_und (und_link);

	return err;
}

char *SrtMultiGrfn3DFXMcEB(
				  char			*underlying,
	              int			numeventdates, 
				  long			*eventdates,
				  int			*optimise,
				  MCEBPARAMS	params,
				  long			*resRows,
                  long			tableauRows,
				  long			*tableauCols,
				  char			***tableauStrings,
				  int			**tableauMask,
				  long			auxWidth,
				  long			*auxLen,
				  double		**aux,				  
				  long			num_paths,  
				  int			do_pecs,
				  double		***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	int					forback;
	int					flag = 0;
	long				nstp;
						
	double				*evt_tms = NULL,
						*sig_dates = NULL,
						**und_sig_dates = NULL,
						**und_sig_curve = NULL,
						**sig_curves = NULL,
						pay_time;

	int					*type = NULL,
						*dom_forex = NULL,
						*for_forex = NULL,
						*fx_index = NULL,
						dom_index;

	double				*corr_date = NULL,
						*tau_date = NULL,
						*tau = NULL,
						tau_temp,
						df;
						
	long				nb_tau,
						*nb_und_sig = NULL,
						nb_corr_date,
						nb_sig_dates,
						pay_date;

	double				*corr = NULL;

	long				*evt_dts = NULL,
						today, spot_date;
						
	int					num_und,
						num_col,
						max_num_df,
						num_evt;
						
	char				*domestic_name,
						*domestic_ccy;

	SrtUndPtr			*und_ptr = NULL,
						dom_und,						
						und;

	void				**void_prm = NULL;
	GRFNPARM_MULTIMC	grfn_prm;
	LINK_UND			und_link = NULL;
	SrtCorrLstPtr		corr_list = NULL;	 
	long				i, j, k, l;
	clock_t				t1, t2;

	int					sig_dates_alloc = 0;

	Err					err = NULL,
						err2 = NULL;

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

	/* look for the today date */
	today = get_today_from_underlying (und_ptr[0]);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	/* look for the main domestic name */
	und = lookup_und (underlying);
	if (!und)
	{
		err = "cannot find the underlying";
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (und) == FX_STOCH_RATES)
	{
		domestic_name = get_domname_from_fxund (und);
	}
	else
	{		
		domestic_name = und -> underl_name;
	}

	domestic_ccy = get_underlying_ccy(und);

	type = (int *) calloc(num_und, sizeof(int));
	dom_forex = (int *) calloc(num_und, sizeof(int));
	for_forex = (int *) calloc(num_und, sizeof(int));
	fx_index  = (int *) calloc(num_und, sizeof(int));

	if (!type || !dom_forex || !for_forex || !fx_index)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	/* Important complete the underlying list */

	err = CompleteUndPtrList(
							  &und_ptr,
							  &num_und,
							  domestic_name,
							  domestic_ccy,
							  &type,
							  &dom_forex,
							  &for_forex,
							  &fx_index,
							  &dom_index
							  );
	
	if (err)
	{
		goto FREE_RETURN;
	}
	
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
	nstp = num_evt;
	
	/* add today to the event dates */
	if (evt_tms[0] > 0)
	{
		/* add the zero time */
		num_f_add_number (&nstp, &evt_tms, 0);
		num_f_sort_vector (nstp, evt_tms);	
		nstp -= 1;
		
		evt_dts = (long*) realloc (evt_dts, (nstp + 1) * sizeof (double));
		
		if (!evt_dts)
		{
			err = "Memory reallocation failure";
			goto FREE_RETURN;
		}
				
		for (i=nstp; i>=1; i--)
		{
			evt_dts[i] = evt_dts[i-1];
		}
		evt_dts[0] = today;
		nstp += 1;
		flag = 1;
	}
		
	/* pay_date is the last event_date */
	pay_date = (long) (evt_dts[nstp-1] + 1.0E-8); 
	pay_time = evt_tms[nstp-1];
	
	/* Fill the link and grfn_prm structure */

	und_link = malloc (sizeof (link_und));
	if (!und_link)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	und_link -> num_und = num_und;
	und_link -> nb_dates = nstp;

	err = allocate_link_und (
							und_link
							);
	if (err)
	{ 
		goto FREE_RETURN;
	}

	und_link -> type = type;
	und_link -> dom_forex = dom_forex;
	und_link -> for_forex = for_forex;
	und_link -> fx_index = fx_index;
	und_link -> dom_index = dom_index;

	type = NULL;
	dom_forex = NULL;
	for_forex = NULL;
	fx_index = NULL;

	/* clean all the underlyings name */
	for (i=0; i<num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	/* get the correlation list */
	corr_list = srt_f_GetTheCorrelationList();

	if (!corr_list)
	{
		err = "Cannot get the correlation structure";
		goto FREE_RETURN;
	}

	/* fill the link structure */
	und_sig_curve = (double **) calloc (num_und, sizeof(double *));
	und_sig_dates = (double **) calloc (num_und, sizeof(double *));
	nb_und_sig = (long *) calloc (num_und, sizeof(long));

	if (!und_sig_curve || !und_sig_dates || !nb_und_sig)
	{
		err = "Memory allocation failure in SrtMultiGrfn3DFXMc";
		goto FREE_RETURN;
	}

	sig_dates_alloc = 0;
	for (i=0; i<num_und; i++)
	{		
		/* first do the LGN Underlyings */
		if (und_link -> type[i] < 2)		
		{
			dom_und = lookup_und (und_ptr[i] -> underl_name);
			if (!dom_und)
			{
				err = serror ("Couldn't find underlying named %s", und_ptr[i] -> underl_name);
				goto FREE_RETURN;
			}
			und_link -> yc[i] = (char *) get_ycname_from_irund (dom_und);
			err = Get_LGM_TermStructure (und_ptr[i] -> underl_name,
										&(und_sig_dates[i]), &(und_sig_curve[i]), &(nb_und_sig[i]), &tau_date, &tau, &nb_tau);
			if (err)
			{				
				goto FREE_RETURN;
			}
			
			if (nb_tau > 1)
			{
				tau_temp = tau[0];
				for (j=1; j<nb_tau; j++)
				{
					if (fabs(tau[j] - tau_temp) > 1.0E-08)
					{
						err = "No Tau term structure is allowed on LGM";
						goto FREE_RETURN;
					}
				}					
			}
			und_link -> lambda[i] = 1.0 / tau[0];

			if (!sig_dates_alloc)
			{
				nb_sig_dates = nb_und_sig[i];
				sig_dates = (double*) calloc (nb_sig_dates, sizeof (double));
				if (!sig_dates)
				{
					err = "Memory allocation failure";
					goto FREE_RETURN;
				}
				
				memcpy (sig_dates, und_sig_dates[i], nb_sig_dates * sizeof (double));

				sig_dates_alloc = 1;
			}
			else
			{
				num_f_concat_vector (&nb_sig_dates, &sig_dates, nb_und_sig[i], und_sig_dates[i]);
			}

			if (tau_date)
			{
				free (tau_date);
				tau_date = NULL;
			}
			if (tau)
			{
				free (tau);
				tau = NULL;
			}
		}
	}

	for (i=0; i<num_und; i++)
	{		
		if (und_link -> type[i]  == 2)
		{						
			err = srt_f_display_FX_TermStruct (und_ptr[i] -> underl_name, &(nb_und_sig[i]),
				&(und_sig_dates[i]), &(und_sig_curve[i]), &nb_corr_date, &corr_date, &corr);
			if (err)
			{				
				goto FREE_RETURN;
			}
			
			for (l=0; l<nb_und_sig[i]; l++)
			{
				und_sig_dates[i][l] = (und_sig_dates[i][l] - today) / 365.0;
			}

			und_link -> spot_fx[i] = get_spot_from_fxund (und_ptr[i]) 
						* swp_f_df (today, spot_date, und_link -> yc[und_link -> dom_forex[i]])
						/ swp_f_df (today, spot_date, und_link -> yc[und_link -> for_forex[i]]);

			und_link -> lambda[i] = und_link -> lambda[und_link -> dom_forex[i]];
			
			for (l=0; l < nb_corr_date; l++) corr_date[l] = (corr_date[l] - today) / 365.0;

			num_f_concat_vector (&nb_sig_dates, &sig_dates, nb_corr_date, corr_date);
			num_f_concat_vector (&nb_sig_dates, &sig_dates, nb_und_sig[i], und_sig_dates[i]);

			if (corr_date) 
			{
				free (corr_date);
				corr_date = NULL;
			}
			if (corr) 
			{
				free (corr);			
				corr = NULL;
			}
		}				
	}

	num_f_sort_vector (nb_sig_dates, sig_dates);
	num_f_unique_vector (&nb_sig_dates, sig_dates);

	und_link->correlations = f3tensor(0, num_und-1, 0, num_und-1, 0, nb_sig_dates-1);

	/* get the correlations at each date */
	for (l=0; l < nb_sig_dates; l++)
	{
		for (i=0; i<num_und; i++)
		{		
			und_link -> correlations[i][i][l] = 1.0;
			for (j=i+1; j<num_und; j++)
			{		
				err = srt_f_get_corr_from_CorrList(corr_list, und_ptr[i] -> underl_name,
					und_ptr[j] -> underl_name, sig_dates[l], &(und_link -> correlations[i][j][l]));
				if (err)
				{
					goto FREE_RETURN;
				}
				und_link -> correlations[j][i][l] = und_link -> correlations[i][j][l];
			}		
		}
	}


	/* fill all the term structures */
	sig_curves = dmatrix(0, num_und-1, 0, nb_sig_dates-1);

	for (i=0; i<num_und; i++)
	{
		for (j=0; j<nb_sig_dates; j++)
		{
			sig_curves[i][j] = und_sig_curve[i][Get_Index(sig_dates[j], und_sig_dates[i], nb_und_sig[i])];
		}
		if (und_sig_dates[i])
		{
			free (und_sig_dates[i]);
			und_sig_dates[i] = NULL;
		}
		if (und_sig_curve[i])
		{
			free (und_sig_curve[i]);
			und_sig_curve[i] = NULL;
		}
	}

	if (und_sig_curve)
	{
		free (und_sig_curve);
		und_sig_curve = NULL;
	}
	if (und_sig_dates)
	{
		free (und_sig_dates);
		und_sig_dates = NULL;
	}
	if (nb_und_sig)
	{
		free (nb_und_sig);
		nb_und_sig = NULL;
	}

	und_link -> nb_sig_dates = nb_sig_dates;
	und_link -> sig_dates = sig_dates;
	und_link -> sig_curve = sig_curves;
	
	sig_dates = NULL;
	sig_curves = NULL;

	/* now calculate all the fwd, std, covariances of the underlyings */	
	err = fill_mc_multi_init_corr(
							pay_date,
							pay_time,
							evt_dts,
							evt_tms,
							nstp,
							und_link
							);

	if (err)
	{
		goto FREE_RETURN;
	}

	/* fill the void_parm structure */
	void_prm = (void**) calloc (nstp, sizeof (void*));

	if (!void_prm)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	for (i=num_evt-1; i>=0; i--)
	{		
		if (xStr.evt[i].evt)
		{
			grfn_prm = malloc (sizeof (grfn_parm_multi_mc));
			if (!grfn_prm)
			{
				err = "Memory allocation failure";
				goto FREE_RETURN;
			}

			grfn_prm -> num_df = lvector (0, num_und-1);
			grfn_prm -> df_tms = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> df_dts = (long **) calloc (num_und, sizeof(long *));

			grfn_prm -> dff = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> gam = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> gam2 = (double **) calloc (num_und, sizeof(double *));
			grfn_prm -> do_und = (int *) calloc (num_und, sizeof(int));

			if (!grfn_prm -> num_df || !grfn_prm -> df_tms || !grfn_prm -> df_dts
				|| !grfn_prm -> dff || !grfn_prm -> gam || !grfn_prm -> gam2
				|| !grfn_prm -> do_und)
			{
				err = "Memory allocation failure";
				goto FREE_RETURN;
			}

			grfn_prm -> global = &xStr;
			grfn_prm -> local = xStr.evt + i;
			grfn_prm -> link = und_link;

			for (j=0; j<num_und; j++)
			{
				grfn_prm->num_df[j] = xStr.evt[i].evt->dflen[j];
				grfn_prm->df_tms[j] = xStr.evt[i].evt->dft[j];
				grfn_prm->df_dts[j] = xStr.evt[i].evt->dfd[j];
			
				if (grfn_prm->num_df[j] > 0)
				{
					grfn_prm -> dff[j] = dvector (0, grfn_prm -> num_df[j] - 1);
					grfn_prm -> gam[j] = dvector (0, grfn_prm -> num_df[j] - 1);
					grfn_prm -> gam2[j] = dvector (0, grfn_prm -> num_df[j] - 1);
				
					if (!grfn_prm -> dff[j] || !grfn_prm -> gam[j] || !grfn_prm -> gam2[j])
					{
						err = "Memory allocation failure";
						goto FREE_RETURN;
					}
		
					if (und_link -> type[j] < 2)
					{
						for (k=0; k<grfn_prm->num_df[j]; k++)
						{
							grfn_prm -> dff[j][k] = swp_f_df (evt_dts[i+flag], grfn_prm->df_dts[j][k], und_link -> yc[j]);
							grfn_prm -> gam[j][k] = (1.0 - exp ( - und_link -> lambda[j] * grfn_prm->df_tms[j][k] )) / und_link -> lambda[j];
							grfn_prm -> gam2[j][k] = 0.5 * grfn_prm->gam[j][k] * grfn_prm->gam[j][k] * und_link -> phi[j][i+flag];
						}
					}
					else
					{
						/* in the fx case we are calculating domestic df */
						for (k=0; k<grfn_prm->num_df[j]; k++)
						{
							grfn_prm -> dff[j][k] = swp_f_df (evt_dts[i+flag], grfn_prm->df_dts[j][k], und_link -> yc[j]);
							grfn_prm -> gam[j][k] = (1.0 - exp ( - und_link -> lambda[und_link -> dom_forex[j]] * grfn_prm->df_tms[j][k] )) / und_link -> lambda[und_link -> dom_forex[j]];
							grfn_prm -> gam2[j][k] = 0.5 * grfn_prm->gam[und_link -> dom_forex[j]][k] * grfn_prm->gam[und_link -> dom_forex[j]][k] * und_link -> phi[und_link -> dom_forex[j]][i+flag];
						}
					}

					grfn_prm->do_und[j] = 1;
				}
				else
				{
					grfn_prm->do_und[j] = 0;
				}
			}
			void_prm[i+flag] = (void*) grfn_prm;
		}
		else
		{
			void_prm[i+flag] = NULL;
		}
	}	

	/*	Eventually! call to function */

	*tableauCols = num_col;
	*resRows = max(num_col + 1, nstp);

	if (params->iMultiIndex)
	{
		params->iNbIndex = params->iMultiIndex;
	}
	else
	{
		params->iNbIndex = 1;
	}

	err = mceb_allocate_params(params, nstp);
	if (err) goto FREE_RETURN;

	*prod_val = dmatrix (0, *resRows - 1, 0, 2 + params->iNbIndex);

	if (!*prod_val)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);

	err = mc_main_multi_3dfx(
							/*	Time data */
							num_paths,
							num_col,
							evt_tms,
							evt_dts,
							nstp,
							do_pecs,
							1,
							optimise,
							params,
							und_link,
							void_prm, 
							grfn_payoff_multi_3dfx_mc,
							*prod_val);

	if (err) goto FREE_RETURN;

	/* Recopy Barrier / CoefLin for the moment */
	for (i=0; i<nstp; i++)
	{
		(*prod_val)[i][2] = params->dBarrier[i];

		for (j=0; j<params->iNbIndex; j++)
		{
			(*prod_val)[i][3+j] = params->dCoefLin[i][j+1];
		}
	}

	df = swp_f_zr (today, pay_date, und_link -> yc[dom_index]);
	df = exp (-df * pay_time);

	
	for (i=0; i<num_col+1; i++)
	{
		(*prod_val)[i][0] *= df;
		(*prod_val)[i][1] *= df;
	}

	if (flag)
	{
		for (i=0; i<nstp-1; i++)
		{			
			(*prod_val)[i][2] = (*prod_val)[i+1][2];			
			
			for (k=0; k<params->iNbIndex; k++)
			{
				(*prod_val)[i][3+k] = (*prod_val)[i+1][3+k];
			}
		}
	}

	/*	Add PV of Past */
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;

		
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

	if (corr)
	{ 
		free (corr);
	}

	if (corr_date)
	{
		free (corr_date);
	}

	if (tau_date)
	{
		free (tau_date);
	}

	if (tau)
	{
		free (tau);
	}

	if (void_prm)
	{
		for (i=0; i<num_evt; i++)
		{
			if (void_prm[i+flag])
			{
				grfn_prm = (GRFNPARM_MULTIMC) void_prm[i+flag];
				
				for (l=0; l<num_und; l++)
				{
					if (	grfn_prm->do_und[l] 
						&& grfn_prm->num_df[l] > 0 )
						
					{
						if (grfn_prm->dff[l]) free_dvector (grfn_prm->dff[l], 0, grfn_prm->num_df[l]-1);
						if (grfn_prm->gam[l]) free_dvector (grfn_prm->gam[l], 0, grfn_prm->num_df[l]-1);
						if (grfn_prm->gam2[l]) free_dvector (grfn_prm->gam2[l], 0, grfn_prm->num_df[l]-1);
					}								
				}
				if (grfn_prm->dff)
					free (grfn_prm->dff);
				if (grfn_prm->gam)
					free (grfn_prm->gam);
				if (grfn_prm->gam2)
					free (grfn_prm->gam2);
				if (grfn_prm->do_und)
					free (grfn_prm->do_und);

				free (grfn_prm);
			}
		}

		free (void_prm);
	}

	if (sig_dates)
	{
		free (sig_dates);
	}

	if (sig_curves)
	{
		free_dmatrix (sig_curves, 0, num_und - 1, 0, nb_sig_dates - 1);
	}

	if (type)
	{
		free (type);
	}

	if (dom_forex)
	{
		free (dom_forex);
	}
	
	if (for_forex)
	{
		free (for_forex);
	}
	
	if (fx_index)
	{
		free (fx_index);
	}

	if (und_sig_dates)
	{
		for (i=0; i<num_und; i++)
		{
			if (und_sig_dates[i])
			{
				free (und_sig_dates[i]);
			}
		}
		free (und_sig_dates);
	}

	if (und_sig_curve)
	{
		for (i=0; i<num_und; i++)
		{
			if (und_sig_curve[i])
			{
				free (und_sig_curve[i]);
			}
		}
		free (und_sig_curve);
	}

	if (nb_und_sig)
	{
		free (nb_und_sig);
	}

	free_link_und (und_link);

	return err;
}