#include "srt_h_all.h"
#include "LGM2Fgrfn.h"
#include "LGM2FMC.h"
#include "math.h"

Err payoff_lgm2f_pde(
					double	evt_date,
					double	evt_time,
					void	*func_parm, 
					
					/* Market data	*/										
					void	*yc,
					
					/* Model data	*/
					double	lam1,
					double	lam2,
					double	rho,
					double	phi1,
					double	phi2,
					double	phi12,
					
					/* Gride data	*/
					int		l1,
					int		u1,
					int		l2,
					int		u2,
					double	*r1,
					double	**r2,
											
					/* Vector of results to be updated */
					int		nprod,
					double	***prod_val
					)
{
GRFNPARMLGM2F	total;
GRFNCOMMSTRUCT	global;
FIRSTMktAtT		*local;

static double	temp;
static int		i, j, k;

static Err		err	= NULL;

	/*	Get the event		*/
	total = (GRFNPARMLGM2F) func_parm;
	global = total->global;
	local = total->local;

	/*	Pre calculations		*/
	if (total->num_df > 0)
	{		
		for (k=0; k<total->num_df; k++)
		{
			total -> dff[k] = swp_f_df (evt_date, total->df_dts[k], (char*) yc);
			total -> dff[k] = log(total -> dff[k]);
			total -> gam1[k] = (1.0 - exp(-lam1 * total->df_tms[k])) / lam1;
			total -> gam2[k] = (1.0 - exp(-lam2 * total->df_tms[k])) / lam2;			
			
			total -> gam12[k] = -0.5 * (total ->gam1[k] * total ->gam1[k] * phi1 
										+ total ->gam2[k] * total ->gam2[k] * phi2)
								- rho * total -> gam1[k] * total -> gam2[k] * phi12;			
		}
	}

	/*	Evaluation			*/
	for (i=l1; i<=u1; i++)
	{
		for (j=l2; j<=u2; j++)
		{
			for (k=0; k<total->num_df; k++)
			{
				local->evt->df[0][k] = exp(total -> dff[k] + total -> gam12[k]	- total -> gam1[k] * r1[i]+
																				- total -> gam2[k] * r2[i][j]);
			}

			err = FIRSTEvalEvent(
								global,
								local,
								nprod,
								2,
								NULL,
								NULL,
								prod_val[i][j],
								&temp);
			if (err)
			{
				return err;
			}			
		}
	}
	
	return err;
}

Err payoff_lgm2fTau_pde(
					double	evt_date,
					double	evt_time,
					void	*func_parm, 
					
					/* Market data	*/										
					void	*yc,
					
					/* Model data	*/
					double	*lam,
					double	*ts_time,
					int		nb_ts,
					double	gamma,
					double	rho,
					double	phi1,
					double	phi2,
					double	phi12,
					
					/* Gride data	*/
					int		l1,
					int		u1,
					int		l2,
					int		u2,
					double	*r1,
					double	**r2,
											
					/* Vector of results to be updated */
					int		nprod,
					double	***prod_val
					)
{
GRFNPARMLGM2F	total;
GRFNCOMMSTRUCT	global;
FIRSTMktAtT		*local;

static double	temp, lam1, lam2, gam1, gam2, gam12, fact1, fact1temp, fact2, fact2temp;
static double	ta, tb, t1, t2, dt;
static int		i, j, k, nb_ts_minus1;

static Err		err	= NULL;

	/*	Get the event		*/
	total = (GRFNPARMLGM2F) func_parm;
	global = total->global;
	local = total->local;

	nb_ts_minus1 = nb_ts - 1;

	/*	Pre calculations		*/
	if (total->num_df > 0)
	{
		/* Find first used lambda */
		j = 0;
		while (j < nb_ts && ts_time[j] < evt_time)
		{
			j++;
		}
		j = min(j, nb_ts - 1);

		gam1 = gam2 = gam12 = 0.0;
		lam1 = lam[j];
		lam2 = lam1 + gamma;

		fact1 = 0.0;
		fact2 = 0.0;

		t1 = evt_time;

		for (k=0; k<total->num_df; k++)
		{
			t2 = evt_time + total->df_tms[k];
			total -> dff[k] = swp_f_df (evt_date, total->df_dts[k], (char*) yc);

			ta = t1;
			tb = ts_time[j];

			gam1 = 0.0;
			gam2 = 0.0;
			fact1temp = 0.0;
			fact2temp = 0.0;

			while (tb < t2 && j < nb_ts_minus1)
			{
				dt = (tb - ta);
				
				gam1 += exp(-fact1temp) * (1.0 - exp(-lam1 * dt)) / lam1;				
				gam2 += exp(-fact2temp) * (1.0 - exp(-lam2 * dt)) / lam2;

				fact1temp += lam1 * dt;
				fact2temp += lam2 * dt;

				j++;
				ta = tb;
				tb = ts_time[j];
				lam1 = lam[j];
				lam2 = lam1 + gamma;
			}

			dt = (t2 - ta);
			
			gam1 += exp(-fact1temp) * (1.0 - exp(-lam1 * dt)) / lam1;			
			gam2 += exp(-fact2temp) * (1.0 - exp(-lam2 * dt)) / lam2;

			fact1temp += lam1 * dt;
			fact2temp += lam2 * dt;

			if (k > 0)
			{
				total -> gam1[k] = total -> gam1[k-1] + exp(-fact1) * gam1;
				total -> gam2[k] = total -> gam2[k-1] + exp(-fact2) * gam2;

			}
			else
			{
				total -> gam1[k] = gam1;
				total -> gam2[k] = gam2;
			}

			total -> gam12[k] = -0.5 * (total ->gam1[k] * total ->gam1[k] * phi1 
										+ total ->gam2[k] * total ->gam2[k] * phi2)
								- rho * total -> gam1[k] * total -> gam2[k] * phi12;

			fact1 += fact1temp;
			fact2 += fact2temp;			

			t1 = t2;
		}
	}

	/*	Evaluation			*/
	for (i=l1; i<=u1; i++)
	{
		for (j=l2; j<=u2; j++)
		{
			for (k=0; k<total->num_df; k++)
			{
				local->evt->df[0][k] = total -> dff[k] * exp(	total -> gam12[k]	- total -> gam1[k] * r1[i]+
																					- total -> gam2[k] * r2[i][j]);
			}

			err = FIRSTEvalEvent(
								global,
								local,
								nprod,
								2,
								NULL,
								NULL,
								prod_val[i][j],
								&temp);
			if (err)
			{
				return err;
			}			
		}
	}
	
	return err;
}

Err grfn_payoff_lgm2f_mc(
						/* Event */
						double	evt_date,
						double	evt_time,
						void	*func_parm,
						/* Market data */
						double	R1D,
						double	R2D,
						/* Results */
						int		num_col,
						double	*res,
						int		*stop_path)
{
	GRFNPARMMC2F	total;
	GRFNCOMMSTRUCT	global;
	FIRSTMktAtT		*local;
	int				l;
	double			temp;
	Err				err;

	if (func_parm == NULL)
	{
		memset (res, 0, num_col * sizeof (double));
		*stop_path = 0;
		return NULL;
	}
	else
	{
		err = NULL;
	}

	memset (res, 0, num_col * sizeof (double));

	/* Get the event */
	total = (GRFNPARMMC2F) func_parm;
	global = total->global;
	local = total->local;
		
	/* Calc market data */	
	if (total->do_dom)
	{
		for (l=0; l<total->num_dom_df; l++)
		{			
			local->evt->df[0][l] = total->dom_dff[l] * exp(	total->dom_gam12[l]	- total->dom_gam1[l] * R1D +
																				- total->dom_gam2[l] * R2D);
		}
	}	
	
	err = FIRSTEvalEvent(
		global,
		local,
		num_col,
		2,
		NULL,
		NULL,
		res,
		&temp);

	*stop_path = 0;
	return err;
}
