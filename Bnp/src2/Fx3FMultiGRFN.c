
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "math.h"

Err allocate_link_und(						
						LINK_UND	link
					 )
{
int nb_dates, num_und;
Err err = NULL;

	if (link)
	{
		link -> type = NULL;
		link -> dom_forex = NULL;
		link -> for_forex = NULL;
		link -> fx_index = NULL;
		link -> yc = NULL;
		link -> lambda = NULL;
		link -> spot_fx = NULL;
		link -> sig_dates = NULL;
		link -> sig_curve = NULL;
		link -> phi = NULL;
		link -> fwd = NULL;
		link -> std = NULL;
		link -> beta = NULL;
		link -> dom_bond_pay = NULL;
		link -> dom_beta_pay = NULL;
		link -> correlations = NULL;
		link -> covariance = NULL;

		nb_dates = link -> nb_dates;
		num_und = link -> num_und;
	
		link -> yc = (char **) calloc (num_und, sizeof(char *));
		link -> lambda = (double *) calloc (num_und, sizeof(double));
		link -> spot_fx = (double *) calloc (num_und, sizeof(double));

		link -> phi = dmatrix(0, num_und-1, 0, nb_dates-1);
		link -> fwd = dmatrix(0, num_und-1, 0, nb_dates-1);
		link -> std = dmatrix(0, num_und-1, 0, nb_dates-1);
		link -> beta = dmatrix(0, num_und-1, 0, nb_dates-1);

		link -> dom_bond_pay = dvector(0, nb_dates-1);
		link -> dom_beta_pay = dvector(0, nb_dates-1);

		link -> covariance = f3tensor (0, nb_dates-1, 0, num_und-1, 0, num_und-1);

		if (!link -> yc || !link -> lambda || !link -> spot_fx
			|| !link -> phi || !link -> fwd || !link -> std || !link -> beta
			|| !link -> dom_bond_pay || !link -> dom_beta_pay || !link -> covariance
			
			)	
			err = "Memory allocation erro in allocate_link_und";				
	}
	else
	{
		err = "Link_Und not allocated ";
	}
					
	return err;
}

Err free_link_und(
					LINK_UND link
					)
{
int nb_dates, num_und;
Err err = NULL;

	if (link)
	{
		nb_dates = link -> nb_dates;
		num_und = link -> num_und;
		
		if (link -> phi)
			free_dmatrix(link -> phi, 0, num_und-1, 0, nb_dates-1); 

		if (link -> fwd)
			free_dmatrix(link -> fwd, 0, num_und-1, 0, nb_dates-1); 

		if (link -> std)
			free_dmatrix(link -> std, 0, num_und-1, 0, nb_dates-1); 

		if (link -> beta)
			free_dmatrix(link -> beta, 0, num_und-1, 0, nb_dates-1); 

		if (link -> dom_bond_pay)
			free_dvector(link -> dom_bond_pay, 0, nb_dates-1);

		if (link -> dom_beta_pay)
			free_dvector(link -> dom_beta_pay, 0, nb_dates-1);

		if (link -> correlations)
			free_f3tensor(link -> correlations, 0, num_und-1, 0, num_und-1, 0, link -> nb_sig_dates-1);

		if (link -> covariance)
			free_f3tensor(link -> covariance, 0, nb_dates-1, 0, num_und-1, 0, num_und-1);

		if (link -> sig_curve)
			free_dmatrix(link -> sig_curve, 0, num_und-1, 0, link -> nb_sig_dates-1);

		if (link -> sig_dates)
			free (link -> sig_dates);

		if (link -> type)
			free (link -> type);

		if (link -> dom_forex)
			free (link -> dom_forex);

		if (link -> for_forex)
			free (link -> for_forex);

		if (link -> fx_index)
			free (link -> fx_index);

		if (link -> yc)
			free (link -> yc);

		if (link -> lambda)
			free (link -> lambda);

		if (link -> spot_fx)
			free (link -> spot_fx);
				
		free (link);

		link = NULL;
	}
	
	return err;
}

Err grfn_payoff_multi_3dfx_mc (
								/* Event */
								double		evt_date,
								double		evt_time,
								void		*func_parm, 
								/* Market data */
								LINK_UND	link,
								double		*sv,
								/* Results */
								int			num_col,
								double		*res)
{
	GRFNPARM_MULTIMC		total;
	GRFNCOMMSTRUCT			global;
	FIRSTMktAtT				*local;
	int						l, und;
	double					temp;
	Err						err = NULL;

	memset (res, 0, num_col * sizeof (double));
	
	if (func_parm == NULL)
	{	
		return NULL;
	}

	/* Get the event */
	total = (GRFNPARM_MULTIMC) func_parm;
	global = total->global;
	local = total->local;
		
	/* Calc market data */
	for (und=0; und<link -> num_und; und++)
	{
		if (link -> type[und] == 2)
		{
			FIRSTSetSVValue(local, und, SPOT, link -> spot_fx[und] * exp (sv[und]));
			
			if (total->do_und[und])
			{
				for (l=0; l<total->num_df[und]; l++)
				{
					FIRSTSetDFValue(local, und, l,					
						 total->dff[und][l] * exp (-total->gam[und][l] * sv[link -> dom_forex[und]] - total->gam2[und][l]));
				}
			}
		}
		else
		{
			if (total->do_und[und])
			{
				for (l=0; l<total->num_df[und]; l++)
				{
					FIRSTSetDFValue(local, und, l,
						total->dff[und][l] * exp (-total->gam[und][l] * sv[und] - total->gam2[und][l]));					
				}
			}
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
	
	return err;
}
