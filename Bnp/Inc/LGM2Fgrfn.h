#ifndef LGM2FGrfnH
#define	LGM2FGrfnH

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

} grfn_parm_lgm, *GRFNPARMLGM2F;


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
					);

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
					);

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
						int		*stop_path);

#endif