
#ifndef Fx3FMCBETADLMH
#define	Fx3FMCBETADLMH

#include "MCEBOptimisation.h"

/*	Main function */
/*	------------- */
Err	 mc_main_3dBetaDLMfx(	
					/*	Time data */
					long		nb_paths,
					int			nb_col,
					double		*time,
					double		*date,					
					long		nb_dates,

					/*	Model Data */
					double		*dom_fwd,
					double		*dom_std,
					double		*for_fwd_const,
					double		*for_fwd_lin,
					double		*for_std,
					double		*fx_fwd,
					double		*fx_std,
					
					double		*bond_pay_const,
					double		*bond_pay_lin,

					double		*dom_for_cov,
					double		*dom_fx_cov,
					double		*for_fx_cov,
					
					/*	Product data */
					void		**func_parm_tab, 
															
					/* do PECS adjustment */
					int			do_pecs,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,

					/*	Initialisation function to be called at the beggining of each path 
							or NULL if none */
					void (*init_func)(),

					/*	Payoff function */
					Err (*payoff_func)(	/* Event */
										double	evt_date,
										double	evt_time,
										void	*func_parm, 
										/* Market data */										
										double	Xdom,
										double	Yfor,
										double	Zfx,
										/* Results */
										int		nb_col,
										double	*res,
										int		*stop_path),
					/*	Results */
					double		**res);

#endif
