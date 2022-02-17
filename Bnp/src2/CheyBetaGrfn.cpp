#include "srt_h_all.h"
#include "CheyBetaGrfn.h"
#include "math.h"

Err payoff_cheybeta_pde(
					double	evt_date,
					double	evt_time,
					void	*func_parm, 
					
					/* Market data	*/										
					void	*yc,
					double  ifr,					
					/* Model data	*/
					TermStruct	*ts,					
					/* Gride data	*/
					int		lx,
					int		ux,
					int		lphi,
					int		uphi,
					double	*x,
					double	*phi,
											
					/* Vector of results to be updated */
					int		nprod,
					double	***prod_val
					)
{
GRFNPARMCHEYBETA	total;
GRFNCOMMSTRUCT	global;
FIRSTMktAtT		*local;

static double	temp;
static int		i, j, k;

static Err		err	= NULL;

	/*	Get the event		*/
	total = (GRFNPARMCHEYBETA) func_parm;
	global = total->global;
	local = total->local;

	/*	Pre calculations		*/
	if (total->num_df > 0)
	{		
		for (k=0; k<total->num_df; k++)
		{
			if (total->df_tms[k] >= YEARS_IN_DAY)
			{
				total -> dff[k] = swp_f_df (evt_date, total->df_dts[k], (char*) yc);
				total -> gam[k] = Lambda_func(evt_time,evt_time+total->df_tms[k],ts);
				total -> gam_sqr[k] = 0.5 *total -> gam[k]*total -> gam[k];
			}
			else
			{
				total -> dff[k] = 1.0;
				total -> gam[k] = 0.0;
				total -> gam_sqr[k] = 0.0;
			}
			
		}
	}

	/*	Evaluation			*/
	for (i=lx; i<=ux; i++)
	{
		for (j=lphi; j<=uphi; j++)
		{
			for (k=0; k<total->num_df; k++)
			{
				local->evt->df[0][k] = total -> dff[k] * exp(- total -> gam[k] * x[i] - total -> gam_sqr[k] * phi[j]);
			}

			local->smp.und[0].sv[0] = x[i] + ifr;
			local->smp.und[0].sv[1] = phi[j];

			err = FIRSTEvalEvent(
								global,
								local,
								nprod,
								2,
								NULL,
								NULL,
								prod_val[j][i],
								&temp);
			if (err)
			{
				return err;
			}			
		}
	}
	
	return err;
}

Err payoff_cheybeta_mc(
					double	evt_date,
					double	evt_time,
					void	*func_parm, 
					
					/* Market data	*/										
					void	*yc,
					/* var of model	*/
					double  r,
					double	x,
					double	phi,				
					/* Vector of results to be updated */
					int		ncols,
					double	*cols_val
					)
{
GRFNPARMCHEYBETA	total;
GRFNCOMMSTRUCT	global;
FIRSTMktAtT		*local;
static double	temp;
static int		k;

Err	err;

	memset (cols_val, 0, ncols * sizeof (double));

	if (func_parm == NULL)
	{
		return NULL;
	}
	else
	{
		err = NULL;
	}

	total = (GRFNPARMCHEYBETA) func_parm;
	global = total->global;
	local = total->local;

	for (k=0; k<total->num_df; k++)
	{
		local->evt->df[0][k] = total -> dff[k] * exp(- total -> gam[k] * x 
			- total -> gam_sqr[k] * phi);
	}
	
	/* set the values of the state variable */
	local->smp.und[0].sv[0] = r;
	local->smp.und[0].sv[1] = phi;

	
	err = FIRSTEvalEvent(
		global,
		local,
		ncols,
		2,
		NULL,
		NULL,
		cols_val,
		&temp);

	return err;
}