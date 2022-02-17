#ifndef CHEYBETAGRFNH
#define	CHEYBETAGRFNH

#include "grf_h_mdlcomm.h"

typedef struct
{
	GRFNCOMMSTRUCT			global;
	FIRSTMktAtT				*local;
	
	long					num_df;
	double					*df_tms;
	long					*df_dts;
	double					*dff;

	double					*gam;
	double                  *gam_sqr;

} grfn_parm_cheybeta, *GRFNPARMCHEYBETA;


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
					);

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
					);
#endif