#ifndef DoubleLGM1FQuanto_pdeH
#define	DoubleLGM1FQuanto_pdeH

#include "srt_h_all.h"


//void mergeVol(int nb_sig1, double *sig1_time, double *sig1, 
//	  int nb_sig2, double *sig2_time, double *sig2,
//			  int *nb_sig, double **sig_time, double **sig);

Err mergeSchedule(int nb1, double *time1, 
				   int nb2, double *time2, 
				   int *nb, double **time);

Err enlargeMesh(int nb_sig, double *time_sig, double *sig, 
				   int nb_mesh, double *time_mesh, 
				   int *nb, double **time, double **longsig);

/*	---------------------------- */
/*	Main function without Tau TS */
/*	---------------------------- */

Err	 doublelgm1fQuanto_adi(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation	
					int			nsteps,
										
						//Model data		
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		quantorho,
					double		domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
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
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res);


Err	 doublelgm1fQuanto_adi2(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation	
					int			nsteps,
										
						//Model data		
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		quantorho,
					double		domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
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
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res);


Err	 doublelgm1fQuanto_adi_correl(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation	
					int			nsteps,
										
						//Model data		
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig, fxsig,
					double		*domsig,	//	quantorho and domforrho must have 	
					double		*forsig,	//	the same sig_time
					double		*fxsig,
					int			nb_sig,
					double		*quantorho,
					double		*domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
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
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res);


Err	 doublelgm1fQuanto_adi_correl2(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation
					int			nsteps,
										
						//Model data
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		*quantorho,
					double		*domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
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
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res);


#endif