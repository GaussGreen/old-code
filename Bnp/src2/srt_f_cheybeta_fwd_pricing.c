
/* ==============================================================================
   
   FILE NAME:      srt_f_cheybeta_fwd_pricing.c 
   
   OBJECT:         use the local vol Cheyette forward PDE to price stuff

  =============================================================================== */

#include "srt_h_all.h"
#include "math.h"

static double Local_Lambda_func(double pt, double GT, 
						 CHEYBETA_MDL *pmdl)
{
	double lambda = pmdl->lambda[0];
	return ( 1.0 - exp(-lambda*(GT-pt)) ) / lambda;
}

	
/*	Init swaption data	*/

/*	Insert one expiry in the global expiry list so as to leave it sorted,
	and links it with a given swaption index; if expiry already exists,
	just adds the link to the swaption	*/
static Err insert_expiry_inst(
double							expiry_time, 
int								swapt_idx,
SWAPTIONS_DATA					*g)
{
	int i, j;

	/*	Case 1: first expiry to be added	*/
	if (g->nex == 0)
	{
		g->nex = 1;
		g->ex[0] = expiry_time;
		g->ndata_e[0] = 1;
		g->data_e[0][0] = swapt_idx;
		return NULL;
	}

	/*	Case 2: expiry already exists	*/
	for (i=0; i<g->nex; i++)
	{
		if (g->ex[i] == expiry_time)
		{
			(g->ndata_e[i])++;
			g->data_e[i][g->ndata_e[i]-1] = swapt_idx;
			return NULL;
		}
	}

	/*	Case 3: Expiry is after the last one	*/
	if (expiry_time > g->ex[g->nex-1])
	{
		g->nex++;
		g->ex[g->nex-1] = expiry_time;
		g->ndata_e[g->nex-1] = 1;
		g->data_e[g->nex-1][0] = swapt_idx;
		return NULL;
	}
	
	/*	Case 4: Expiry is to be inserted	*/

	i = 0;
	while (g->ex[i] < expiry_time)
	{
		i++;
	}

	for (j=g->nex; j>i; j--)
	{
		g->ex[j] = g->ex[j-1];
		g->ndata_e[j] = g->ndata_e[j-1];
		memcpy (&(g->data_e[j][0]), &(g->data_e[j-1][0]), g->ndata_e[j-1] * sizeof (int));
	}

	g->nex++;
	g->ex[i] = expiry_time;
	g->ndata_e[i] = 1;
	g->data_e[i][0] = swapt_idx;
		
	return NULL;
}

/*	Fill one given SWAPTION_DATA structure	*/
static Err fill_swaption_data (	
SWAPTION_DATA					*data,
/*	Underlying info	*/
CHEYBETA_MDL					*pmdl,
char							*yc_name,
long							today,
int								spot_lag,
/*	Swaption info	*/
long							expiry,
double							expiry_time,
long							start,
long							end,
char							*freq,
char							*basis,
double							strike,
double							bond_strike,
char							*pay_rec_str,
char 							*refrate)
{
	SwapDP fix_sdp;
	GenSwapLeg fix_leg, float_leg, big_leg;
	SrtReceiverType pay_rec;
	Err err;
	int i;
	double df_to_exp;

	if (err = interp_rec_pay (pay_rec_str, &pay_rec))
	{
		return err;
	}

	/*	Make fix leg	*/
	fix_sdp.spot_lag = spot_lag;
	if (err = swp_f_initSwapDP (start, end, freq, basis, &fix_sdp))
	{
		return err;
	}
	if (err = swp_f_make_FixedAndNotionalsLeg (	&fix_sdp, 
												strike, 
												1.0, 
												1.0, 
												today, 
												&fix_leg))
	{
		return err;
	}
	fix_leg.rec_pay = pay_rec;

	/*	Make float leg	*/
	if (err = swp_f_make_SpreadLeg (start, end, today, refrate, &float_leg))
	{
		swp_f_freein_GenSwapLeg (&fix_leg);
		return err;
	}
	float_leg.rec_pay = (pay_rec == SRT_RECEIVER? SRT_PAYER: SRT_RECEIVER);

	/*	Merge legs	*/
	if (pay_rec == SRT_RECEIVER)
	{
		err = swp_f_merge_SwapLegs (&fix_leg, &float_leg, &big_leg);
	}
	else
	{
		err = swp_f_merge_SwapLegs (&float_leg, &fix_leg, &big_leg);
	}
	swp_f_freein_GenSwapLeg (&fix_leg);
	swp_f_freein_GenSwapLeg (&float_leg);
	if (err)
	{
		return err;
	}
	big_leg.sdp.spot_lag = spot_lag;


	/*	Compute pay times	*/
	if (err = time_list (big_leg.pay_date, big_leg.today, &(big_leg.time)))
	{
		swp_f_freein_GenSwapLeg (&big_leg);
		return err;
	}

	/*	Compute forward DFs	*/
	if (err = df_list (big_leg.pay_date, yc_name, &(big_leg.df)))
	{
		swp_f_freein_GenSwapLeg (&big_leg);
		return err;
	}

	/*	Copy data	*/

	/*	A few checks first	*/

	if (big_leg.leg_length < 1)
	{
		swp_f_freein_GenSwapLeg (&big_leg);
		return serror ("Instrument has no cash-flows after expiry");
	}

	if (expiry_time < 0)
	{
		swp_f_freein_GenSwapLeg (&big_leg);
		return serror ("Instrument has already expired");
	}

	if (big_leg.time.d[0] < expiry_time)
	{
		swp_f_freein_GenSwapLeg (&big_leg);
		return serror ("Instrument delivers cash-flows before expiration");
	}

	if (fabs (bond_strike - 1.0) > 1.0e-08)
	{
		swp_f_freein_GenSwapLeg (&big_leg);
		return serror (	"Forward Chey Beta PDE can only price "
						"Options with bond strike = 100 so far");
	}

	/*	Number of coupons	*/
	data->ncpn = big_leg.leg_length;
	
	/*	Expiry time	*/
	data->exp_time = expiry_time;
	
	/*	Pay times	*/
	data->pay_time = calloc (data->ncpn, sizeof (double));
	memcpy (data->pay_time, big_leg.time.d, data->ncpn * sizeof (double));

	/*	Coupons	*/
	data->pay_cpn = calloc (data->ncpn, sizeof (double));
	memcpy (data->pay_cpn, big_leg.payment.d, data->ncpn * sizeof (double));

	/*	Forward discount factors	*/
	df_to_exp = swp_f_df (today, expiry, yc_name);
	data->pay_dff = calloc (data->ncpn, sizeof (double));
	for (i=0; i<data->ncpn; i++)
	{
		data->pay_dff[i] = big_leg.df.d[i] / df_to_exp;
	}

	/*	Reconstruction formula parameters	*/
	data->pay_gamma = calloc (data->ncpn, sizeof (double));
	for (i=0; i<data->ncpn; i++)
	{
		data->pay_gamma[i] = Local_Lambda_func ( expiry_time,
											data->pay_time[i],
											pmdl);
	}

	data->pay_0_50_gamma_sqr = calloc (data->ncpn, sizeof (double));
	for (i=0; i<data->ncpn; i++)
	{
		data->pay_0_50_gamma_sqr[i] = 0.50 * data->pay_gamma[i] * data->pay_gamma[i];
	}

	/*	Free leg data	*/
	swp_f_freein_GenSwapLeg (&big_leg);
	return NULL;
}

/*	Adds one option in the swaptions data	*/
Err add_new_swaption (
/*	Underlyin info	*/
CHEYBETA_MDL					*pmdl,
char							*yc_name,
long							today,
int								spot_lag,
/*	Index of the instrument the option belongs to	*/
int								inst_idx,
/*	Option info	*/
long							expiry,
long							start,
long							end,
char							*freq,
char							*basis,
double							strike,
double							bond_strike,
char							*pay_rec,
char							*refrate,
SWAPTIONS_DATA					*g)
{
	double expiry_time = (expiry - today) * YEARS_IN_DAY;
	Err err;

	/*	Increment number of swaptions	*/
	g->nswapt++;

	/*	Fill data	*/
	if (err = fill_swaption_data (	g->data + g->nswapt - 1,
									pmdl,
									yc_name,
									today,
									spot_lag,
									expiry,
									expiry_time,
									start,
									end,
									freq,
									basis,
									strike,
									bond_strike,
									pay_rec,
									refrate))
	{
		return err;
	}

	/*	Adds swaption to the instrument it belongs to	*/
	(g->ndata_i[inst_idx])++;
	g->data_i[inst_idx][g->ndata_i[inst_idx]-1] = g->nswapt-1;

	/*	Adds swaption to the expiry it belongs to	*/
	insert_expiry_inst (expiry_time, g->nswapt-1, g);
	return NULL;
}

/*	Maximum number of caplets in a cap	*/
#define MAX_CAPLET			200
/*	Maximum number of options per expiry	*/
#define MAX_OPT_PER_EXP		200

Err chey_beta_initmem_swaption_data(
int								ninst,
SWAPTIONS_DATA					*g )
{
	g->ninst = ninst;
	g->data = (SWAPTION_DATA*) calloc (MAX_CAPLET * ninst, sizeof (SWAPTION_DATA));
	g->nswapt = 0;
	g->ndata_i = (int*) calloc (ninst, sizeof (int));
	g->data_i = imatrix (0, ninst, 0, MAX_CAPLET-1);
	g->nex = 0;
	g->ex = (double*) calloc (MAX_CAPLET * ninst, sizeof (double));
	g->ndata_e = (int*) calloc (MAX_CAPLET * ninst, sizeof (int));
	g->data_e = imatrix (0, MAX_CAPLET * ninst, 0, MAX_OPT_PER_EXP-1);
	if ( !g->data || !g->ndata_i || !g->data_i ||
		 !g->ex || !g->ndata_e || !g->data_e )
		return serror("Memory failure in chey_beta_initmem_swaption_data");
	return NULL;
}

Err chey_beta_fwd_init_swaption_data(
CHEYBETA_MDL					*pmdl,
int								ninst,
long							*start,
long							*end,
char							**freq_str,
char							**basis_str,
double							*strike,
double							*bond_strike,
char							**type_str,
char							**pay_rec_str,
char							**refrate_str,
SWAPTIONS_DATA					*g)
{
	long today = (long)pmdl->today_date;
	char *yc_name = pmdl->yc_name;
	SrtCurvePtr yccrv = lookup_curve(yc_name);
	int spot_lag = get_spotlag_from_curve(yccrv);
	long expiry;

	int i, j;
	StructType type;
	GenSwapLeg float_leg;
	SrtCompounding rate_freq;
	BasisCode rate_basis;
	char *rate_freq_str, *rate_basis_str;
	Err err;

	memset(g, 0, sizeof(SWAPTIONS_DATA));
	err = chey_beta_initmem_swaption_data(ninst, g);
	if (err)
	{
	   	chey_beta_fwd_free_swaption_data (g);
		return err;
	}
	
	/*	Loop on instrument i	*/
	for (i=0; i<ninst; i++)
	{
		if (err = interp_struct (type_str[i], &type))
		{
		   	chey_beta_fwd_free_swaption_data (g);
			return err;
		}

		if (type == SWAPTION)
		{
			expiry = add_unit (start[i], -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
			if (err = add_new_swaption (	pmdl,
											yc_name,
											today,
											spot_lag,
											i,
											expiry,
											start[i],
											end[i],
											freq_str[i],
											basis_str[i],
											strike[i],
											bond_strike[i],
											pay_rec_str[i],
											refrate_str[i],
											g))
			{
				chey_beta_fwd_free_swaption_data (g);
				return err;
			}	
		}
		else
		if (type == CAPFLOOR)
		{

			if (err = swp_f_get_ref_rate_details (	refrate_str[i], 
													&rate_basis, 
													&rate_freq))
			{
				chey_beta_fwd_free_swaption_data (g);
				return err;
			}

			translate_basis	(&rate_basis_str, rate_basis);
			translate_compounding (&rate_freq_str, rate_freq);

			if (err = swp_f_make_SpreadLeg (	start[i], 
												end[i], 
												today, 
												refrate_str[i], 
												&float_leg))
			{
				chey_beta_fwd_free_swaption_data (g);
				return err;
			}
			
			/*	Loop on caplets	*/
			for (j=0; j<float_leg.leg_length-1; j++)
			{
				expiry = add_unit (float_leg.start_date.date[j],
					-spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
				if (err = add_new_swaption (	pmdl,
												yc_name,
												today,
												spot_lag,
												i,
												expiry,
												float_leg.start_date.date[j],
												float_leg.end_date.date[j],
												rate_freq_str,
												rate_basis_str,
												strike[i],
												bond_strike[i],
												pay_rec_str[i],
												refrate_str[i],
												g))
				{
					chey_beta_fwd_free_swaption_data (g);
					return err;
				}
			}

			swp_f_freein_GenSwapLeg (&float_leg);
		}
		else
		{
			return serror (	"Forward Chey Beta PDE can only price "
							"SWAPTIONs or CAPFLOORs so far");
		}
	}

	return NULL;
}

/*	Update gammas in swaption data as model parameters change	*/

Err chey_beta_fwd_update_swaption_data_for_model(
/*	new model	*/
CHEYBETA_MDL					*pmdl,
SWAPTIONS_DATA					*g)
{
	int i, j;

	for (i=0; i<g->nswapt; i++)
	{

		for (j=0; j<g->data[i].ncpn; j++)
		{
			g->data[i].pay_gamma[j] = Local_Lambda_func (	g->data[i].exp_time, 
													g->data[i].pay_time[j], 
													pmdl);
		
			g->data[i].pay_0_50_gamma_sqr[j] =		0.50 
												*	g->data[i].pay_gamma[j] 
												*	g->data[i].pay_gamma[j];
		}
	}
	
	return NULL;	
}

/*	Free globals	*/

Err chey_beta_fwd_free_swaption_data (
SWAPTIONS_DATA					*g)
{
	int i;

	if (g->data)
	{
		for (i=0; i<g->nswapt; i++)
		{
			free (g->data[i].pay_time);
			free (g->data[i].pay_cpn);
			free (g->data[i].pay_dff);
			free (g->data[i].pay_gamma);
			free (g->data[i].pay_0_50_gamma_sqr);
		}
		free (g->data);
	}

	free (g->ndata_i);
	if (g->data_i) free_imatrix (g->data_i, 0, g->ninst, 0, MAX_CAPLET-1);
	free (g->ex);
	free (g->ndata_e);
	if (g->data_e) free_imatrix (g->data_e, 0, g->ninst, 0, MAX_OPT_PER_EXP-1);

	return NULL;
}

/*	Do price	*/

Err chey_beta_fwd_price(
CHEYBETA_MDL					*pmdl,
int								nt,
int								nx,
int								nphi,
SWAPTIONS_DATA					*g,
/*	Allocated inside, to be freed by caller	*/
int								*ninst,
double							**inst_prem)
{
	CHEYBETA_GRID grid;
	CHEYBETA_FWD_PDE_TEMP temp;

	double *x, *phi, **ad;
	int stp_idx, exp_idx;
	int n;
	int i, j, k, l;
	double temp_df, iv;

	/*	Init chey beta fwd pode grids and structures	*/
	chey_beta_grid_init (&grid);
	chey_beta_grid_time (pmdl, g->nex, g->ex, nt, &grid);
	chey_beta_grid_phi_x (pmdl, &grid, nx, nphi);
	chey_beta_grid_ad (&grid);
	srt_f_chey_beta_fwd_pde_alloc_temp (&temp, nx, nphi);
	srt_f_chey_beta_fwd_pde_first_step (pmdl, &grid);
	
	/*	Start @ step 1	*/
	stp_idx = 1;
	/*	Start @ expiry 0	*/
	exp_idx = 0;

	/*	First, treat instruments that expire today	*/
	if (g->ex[exp_idx] < 1.0e-08)
	{
		/*	Loop on swaptions	*/
		for (i=0; i<g->ndata_e[exp_idx]; i++)
		{
			/*	Initialise premium to 0	*/
			g->data[g->data_e[exp_idx][i]].prem = 0.0;

			/*	Calc intrinsic	*/
			iv = 0.0;

			for (l=0; l<g->data[g->data_e[exp_idx][i]].ncpn; l++)
			{
				temp_df
					= g->data[g->data_e[exp_idx][i]].pay_dff[l];

				iv += temp_df
					* g->data[g->data_e[exp_idx][i]].pay_cpn[l];
			}
			
			if (iv > 0.0)
			{
				g->data[g->data_e[exp_idx][i]].prem = iv;
			}
		}
	
		/*	Next expiry	*/
		exp_idx++;
	}

	do
	{
		n = 0;
		/*	Find step corresponding to next expiry	*/
		while (fabs (grid.t[stp_idx+n] - g->ex[exp_idx]) > 1.0e-08)
		{
			n++;
		}

		/*	Go there	*/
		srt_f_chey_beta_fwd_pde_n_steps (&temp, pmdl, &grid, &stp_idx, n);
		
		/*	Copy density	*/
		chey_beta_density_at_t (&grid, stp_idx, &nx, &x, &nphi, &phi, &ad);

		/*	Loop on swaptions	*/
		for (i=0; i<g->ndata_e[exp_idx]; i++)
		{
			/*	Initialise premium to 0	*/
			g->data[g->data_e[exp_idx][i]].prem = 0.0;

			/*	Loop on x	*/
			for (j=0; j<nx; j++)
			{
				/*	Loop on phi	*/
				for (k=0; k<nphi; k++)
				{
					/*	Calc intrinsic	*/
					iv = 0.0;

					for (l=0; l<g->data[g->data_e[exp_idx][i]].ncpn; l++)
					{
						temp_df
							= g->data[g->data_e[exp_idx][i]].pay_dff[l]
							* exp ( - g->data[g->data_e[exp_idx][i]].pay_gamma[l] * x[j]
						- g->data[g->data_e[exp_idx][i]].pay_0_50_gamma_sqr[l] * phi[k]);

						iv += temp_df
							* g->data[g->data_e[exp_idx][i]].pay_cpn[l];
					}
					
					if (iv > 0.0)
					{
						g->data[g->data_e[exp_idx][i]].prem += ad[k][j] * iv;
					}
				}
			}
		}

		chey_beta_free_density_at_t (&nx, &x, &nphi, &phi, &ad);
		exp_idx++;

	} while (stp_idx < grid.nt - 1);

	/*	Copy results	*/

	*ninst = g->ninst;
	*inst_prem = (double*) calloc (*ninst, sizeof (double));

	for (i=0; i<g->ninst; i++)
	{
		for (j=0; j<g->ndata_i[i]; j++)
		{
			(*inst_prem)[i] += g->data[g->data_i[i][j]].prem;
		}
	}

	srt_f_chey_beta_fwd_pde_free_temp (&temp);
	chey_beta_grid_free (&grid);

	return NULL;
}
