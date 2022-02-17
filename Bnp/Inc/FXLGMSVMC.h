#ifndef FXLGMSVMC
#define	FXLGMSVMC

#include "srt_h_all.h"
#include "math.h"
#include "LGMSVUtil.h"


Err	 fxlgmsv_mc_balsam(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,

					//Domestic
					/*	Model data Information	*/
					double		ddomLambdaX1,
					double		ddomLambdaX2,

					double		*ddomSigma,
					double		*ddomAlphaLGM,
					double		*ddomRhoLGM,
										
					double		*ddomAlpha,
					double		*ddomLambdaEps,
					double		*ddomLvlEps,
					double		*ddomRho,
					double		*ddomRho2,

					double		*domzcvol1_star,
					double		*domzcvol2_star,

					/* Parameters for DF(t,T*) reconstruction */
					double		*ddomff_star,
					double		*domgam1_star,
					double		*domgam2_star,
					double		*domgam1_2_star,
					double		*domgam2_2_star,
					double		*domgam12_star,

					//Foreign
					/*	Model data Information	*/
					double		dforLambdaX1,
					double		dforLambdaX2,

					double		*dforSigma,
					double		*dforAlphaLGM,
					double		*dforRhoLGM,
										
					double		*dforAlpha,
					double		*dforLambdaEps,
					double		*dforLvlEps,
					double		*dforRho,
					double		*dforRho2,

					double		*forzcvol1_star,
					double		*forzcvol2_star,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dforff_star,
					double		*forgam1_star,
					double		*forgam2_star,
					double		*forgam1_2_star,
					double		*forgam2_2_star,
					double		*forgam12_star,

					//FX
					double		fwd_fx_TStar,
					double		*fx_vol,

					//	Correlation
					double		***CorrMatrix,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										//Domestic
										/* Model data	*/
										double	domft1,
										double	domft2,
										double	domphi1,
										double	domphi2,
										double	domphi12,
										double	domv,
																
										//Foreign
										/* Model data	*/
										double	forft1,
										double	forft2,
										double	forphi1,
										double	forphi2,
										double	forphi12,
										double	forv,

										//FX
										double	fx_spot,

										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),
					/*	Result */
					int			iNbProduct, 
					double		**res);


Err	 qtolgmsv2f_mc_balsam(	

					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,

					//Domestic

					double		ddomLambdaX,

					double		*ddomSigma,

					double		*domzcvol_star,


					double		*ddomff_star,
					double		*domgam1_star,
					double		*domgam1_2_star,

					//Foreign

					double		dforLambdaX1,
					double		dforLambdaX2,

					double		*dforSigma,
					double		*dforAlphaLGM,
					double		*dforRhoLGM,
										
					double		*dforAlpha,
					double		*dforLambdaEps,
					double		*dforLvlEps,
					double		*dforRho,
					double		*dforRho2,

					double		*forzcvol1_star,
					double		*forzcvol2_star,


					double		*dforff_star,
					double		*forgam1_star,
					double		*forgam2_star,
					double		*forgam1_2_star,
					double		*forgam2_2_star,
					double		*forgam12_star,

					//FX Vol
					double		*fx_vol,

					//	Correlation
					double		***CorrMatrix,
					/*	0 : Dom
						1 : For1
						2 : For2
						3 : ForSV
						4 : FX	*/


					void		**func_parm_tab, 
					int			*EvalEvent,


					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					
					void		(*init_func)(),
									
					// Payoff function of QUANTOLGMSV2F
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										//Domestic
										double	domft,
										double	domphi,
																
										//Foreign										
										double	forft1,
										double	forft2,
										double	forphi1,
										double	forphi2,
										double	forphi12,
										double	forv,

										
										int		nprod,
										
										double	*prod_val,
										int		*stop_path),
					
					int			iNbProduct, 
					double		**res);

/* ORIGINAL CALL
Err	 qtolgmsv2f_mc_balsam(	
					//	Time Information
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,

					//Domestic
					//	Model data Information
					double		ddomLambdaX,

					double		*ddomSigma,

					double		*domzcvol_star,

					// Parameters for DF(t,T*) reconstruction
					double		*ddomff_star,
					double		*domgam1_star,
					double		*domgam1_2_star,

					//Foreign
					//	Model data Information
					double		dforLambdaX1,
					double		dforLambdaX2,

					double		*dforSigma,
					double		*dforAlphaLGM,
					double		*dforRhoLGM,
										
					double		*dforAlpha,
					double		*dforLambdaEps,
					double		*dforLvlEps,
					double		*dforRho,
					double		*dforRho2,

					double		*forzcvol1_star,
					double		*forzcvol2_star,

					// Parameters for DF(t,T*) reconstruction
					double		*dforff_star,
					double		*forgam1_star,
					double		*forgam2_star,
					double		*forgam1_2_star,
					double		*forgam2_2_star,
					double		*forgam12_star,

					//FX Vol
					double		*fx_vol,

					//	Correlation
					double		***CorrMatrix,
					//	0 : Dom
					//	1 : For1
					//	2 : For2
					//	3 : ForSV
					//	4 : FX

					//	Product data
					void		**func_parm_tab, 
					int			*EvalEvent,

					// for Optimisation of exercise boundary
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					//	Initialisation function to be called at the beggining of each path or NULL if none
					void		(*init_func)(),
					
					//	Payoff function
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										//Domestic
										// Model data
										double	domft,
										double	domphi,
																
										//Foreign
										//Model data
										double	forft1,
										double	forft2,
										double	forphi1,
										double	forphi2,
										double	forphi12,
										double	forv,

										//Vector of results to be updated
										int		nprod,
										//Result
										double	*prod_val,
										int		*stop_path),
					//Result
					int			iNbProduct, 
					double		**res);
					*/


Err	 qtolgmsv1f_mc_balsam(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,

					//Domestic
					/*	Model data Information	*/
					double		ddomLambdaX,

					double		*ddomSigma,

					double		*domzcvol_star,

					/* Parameters for DF(t,T*) reconstruction */
					double		*ddomff_star,
					double		*domgam1_star,
					double		*domgam1_2_star,

					//Foreign
					/*	Model data Information	*/
					double		dforLambdaX1,

					double		*dforSigma,
										
					double		*dforAlpha,
					double		*dforLambdaEps,
					double		*dforLvlEps,
					double		*dforRho,

					double		*forzcvol_star,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dforff_star,
					double		*forgam1_star,
					double		*forgam1_2_star,

					//FX Vol
					double		*fx_vol,

					//	Correlation
					double		***CorrMatrix,
					/*	0 : Dom
						1 : For1
						2 : ForSV
						3 : Fx	*/

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										//Domestic
										/* Model data	*/
										double	domft,
										double	domphi,
																
										//Foreign
										/* Model data	*/
										double	forft,
										double	forphi,
										double	forv,

										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),
					/*	Result */
					int			iNbProduct, 
					double		**res);


#endif
