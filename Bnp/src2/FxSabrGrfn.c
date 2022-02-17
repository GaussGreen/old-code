#include "srt_h_all.h"
#include "FxSabrGrfn.h"
#include "math.h"

Err payoff_fx_sabr_adi(
					double	evt_date,
					double	evt_time,
					void	*func_parm, 
															
					/* Market data	*/										
					long	today,
					double	spot_fx,	/*	The cash one */
					void	*dom_yc,
					void	*for_yc,
															
					/* Grid data	*/
					int		l1,
					int		u1,
					int		l2,
					int		u2,
					double	*x,
											
					/* Vector of results to be updated */
					int		nprod,
					double	***prod_val
					)
{
GRFNPARMFXSABR	total;
GRFNCOMMSTRUCT	global;
FIRSTMktAtT		*local;

static double	df_ratio;
static double	temp;
static int		i, j, k;

static Err		err	= NULL;

	/*	Get the event		*/
	total = (GRFNPARMFXSABR) func_parm;
	global = total->global;
	local = total->local;

	/*	Pre calculations		*/
		
	df_ratio = swp_f_df (today, evt_date, (char*) for_yc) / swp_f_df (today, evt_date, (char*) dom_yc);						

	if (total->num_fx_df > 0)
	{		
		for (k=0; k<total->num_fx_df; k++)
		{
			local->evt->df[total->fx_idx][k] = swp_f_df (evt_date, total->fx_df_dts[k], (char*) dom_yc);
		}
	}

	if (total->num_dom_df > 0)
	{		
		for (k=0; k<total->num_dom_df; k++)
		{
			local->evt->df[total->dom_idx][k] = swp_f_df (evt_date, total->dom_df_dts[k], (char*) dom_yc);
		}
	}

	if (total->num_for_df > 0)
	{		
		for (k=0; k<total->num_for_df; k++)
		{
			local->evt->df[total->for_idx][k] = swp_f_df (evt_date, total->for_df_dts[k], (char*) for_yc);
		}
	}

	/*	Evaluation			*/
	for (i=l1; i<=u1; i++)
	{
		local->smp.und[total->fx_idx].sv[SPOT] = x[i] * df_ratio;

		for (j=l2; j<=u2; j++)
		{
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

Err payoff_fx_sabr_adi_ATM_opt(
								double	evt_date,
								double	evt_time,
								void	*func_parm, 
																		
								/* Market data	*/										
								long	today,
								double	spot_fx,	/*	The cash one */
								void	*dom_yc,
								void	*for_yc,
																		
								/* Grid data	*/
								int		l1,
								int		u1,
								int		l2,
								int		u2,
								double	*x,
														
								/* Vector of results to be updated */
								int		nprod,
								double	***prod_val
								)
{
int		i, j;
double	pay, strike;

	i = l1;

	strike = ((double*) (func_parm))[0];

	while (x[i] < strike && i <= u1)
	{
		for (j=l2; j<=u2; j++)
		{									
			prod_val[i][j][0] = 0.0;
		}

		i++;
	}

	
	for (i; i<=u1; i++)
	{		
		pay = x[i] - strike;
		for (j=l2; j<=u2; j++)
		{									
			prod_val[i][j][0] = pay;
		}
	}		

	return NULL;
}

Err payoff_fx_sabr_adi_opt(
							double	evt_date,
							double	evt_time,
							void	*func_parm, 
																	
							/* Market data	*/										
							long	today,
							double	spot_fx,	/*	The cash one */
							void	*dom_yc,
							void	*for_yc,
																	
							/* Grid data	*/
							int		l1,
							int		u1,
							int		l2,
							int		u2,
							double	*x,
													
							/* Vector of results to be updated */
							int		nprod,
							double	***prod_val
							)
{
int		i, j, k;
double	pay, df_ratio, strike2;
double	*strike;		
long	mat_date;

	mat_date = add_unit ((long) (evt_date + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	strike = (double *) func_parm;	
	df_ratio = swp_f_df (today, mat_date, (char*) for_yc) / swp_f_df (today, mat_date, (char*) dom_yc);	
	
	for (k=0; k<nprod; k++)
	{
		strike2 = strike[k] / df_ratio;

		i = l1;

		while (x[i] < strike2 && i <= u1)
		{
			for (j=l2; j<=u2; j++)
			{									
				prod_val[i][j][k] = 0.0;
			}

			i++;
		}		

		for (i; i<=u1; i++)
		{
			pay = (x[i] - strike2) * df_ratio;
			for (j=l2; j<=u2; j++)
			{									
				prod_val[i][j][k] = pay;
			}
		}
	}

	return NULL;
}