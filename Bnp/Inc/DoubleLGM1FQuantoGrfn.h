#ifndef DoubleLGM1FQuantoGrfnH
#define	DoubleLGM1FQuantoGrfnH

#include "grf_h_mdlcomm.h"

typedef struct
{
	GRFNCOMMSTRUCT			global;
	FIRSTMktAtT				*local;
	
	long					num_df;
	double					*df_tms;
	long					*df_dts;

	double					*dff;
	double					*gam1;
	double					*gam2;
	double					*gam12;

} grfn_parm_lgm, *GRFNPARMDBLELGM1FQTO;


Err payoff_doublelgm1fquanto_pde(
					double	evt_date,
					double	evt_time,
					void	*func_parm, 
					
					/* Market data	*/										
					void	*dom_yc,
					void	*for_yc,
					
					/* Model data	*/
					double	domlam,
					double	forlam,
					double	domphi,
					double	forphi,
					
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
					);

#endif