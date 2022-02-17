#ifndef DIAGCALIBGENH
#define	DIAGCALIBGENH

#include "srt_h_all.h"
#include "CPDCalib.h"
#include "DiagCalibDLM.h"

double solve_for_next_coef_gen(	double	**res_iter,
								int		nb_iter,
								double	premium_tgt,
								double	precision,
								int		method);			/* 0: linear 1: quadratic */

void solve_for_next_sensi(	double	**res_iter,
							int		nb_computed,
							double	value,
							double	delta,
							double	*sensi);

typedef struct
{	
	int		UseJumps;
	double	Precision;
	int		NbIterMax;
	int		RecalibAtEnd;

	int		do_calib;
	int		retry_after_fail;

	int		nb_computed;
	double	**res_iter;

	int		compute_sensi;
	double	delta;
	double	sensi;

	int		total_computed;

} CALIBGEN_Params, *CALIBGEN_PARAMS;

Err	Initialise_CalibParams(	int						UseJumps,
							double					Precision,
							int						NbIterMax,
							int						RecalibAtEnd,
							int						do_calib,
							int						retry_after_fail,
							int						compute_sensi,
							double					sensi_minshift,
							CALIBGEN_PARAMS			CalibParams);

Free_CalibParams(CALIBGEN_PARAMS		CalibParams);

typedef struct
{
	double	dVolPrecision;
	int		iNbVolIterMax;

	double	dLamPrecision;
	int		iNbLamIterMax;

	double	dSmilePrecision;
	int		iNbSmileIterMax;

} diag_calib_gen_params, DIAG_CALIB_GEN_PARAMS;

typedef struct
{
	/* All functions for fitting of one parameter */

	Err	(*GetTarget)(		void				*Inst,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target);

	/* For Newton first iterations */
	Err	(*GetFirstGuess)(	void			*Model,
							void			*GlobalConst,
							int				index_param,
							double			target,
							double			*param1);

	Err	(*GetSecondGuess)(	void			*Model,
							void			*GlobalConst,
							int				index_param,
							double			param1,
							double			price1,
							double			target,
							double			*param2);

	/* Set the bounds and the default values */
	Err	(*GetLimitAndLastParam)(void				*Model,
								CALIBGEN_PARAMS		CalibConsts,
								void				*GlobalConst,
								int					index_param,
								double				*last_param,
								double				*limit_down,
								double				*limit_up);

	/* Bump the parameter */
	Err	(*BumpParam)(		void			*Model,
							void			*GlobalConst,
							int				index_param,
							double			param);

	/* Set the parameter at its final value */
	Err	(*SetParam)(		void				*Model,
							CALIBGEN_PARAMS		CalibConsts,
							void				*GlobalConst,
							int					index_param,
							double				param);

	/* At the end of the calibration, fill the remaining parameters */
	Err	(*ExtrapolParam)(	void			*Model,
							void			*GlobalConst,
							int				LastParamIndex);

	/* Update of pricing consts after the parameter change */																			
	Err	(*UpdateConstsAfterParam) (	void				*Inst,
									void				*InstConst,
									void				*GlobalConst,
									void				*Model,
									CALIBGEN_PARAMS	CalibConsts);

	/* Pricing of the instrument */
	Err	(*PriceInst)(	void				*Inst,
						void				*InstConst,
						void				*GlobalConst,
						void				*Model,						
						double				*InstPrice);

} CalibFunctions, *CALIBFUNCTIONS;


Err	CalibrateNextParam(	void				*Inst,
						void				*InstConst,
						void				*GlobalConst,
						void				*Model,							
											
						CALIBGEN_PARAMS	CalibConsts,
											
						int					index_param,
						double				last_param,							
						double				limit_down,
						double				limit_up,

						CALIBFUNCTIONS		AllFunctions,
											
						double				*res_param,
						int					*success);


Err	CalibrateParamTS(	int						StartIndexInst,
						int						EndIndexInst,
						void					**AllInst,
						void					**InstConst,
						void					*GlobalConst,
						void					*Model,							
											
						CALIBGEN_PARAMS		CalibConsts,

						CALIBFUNCTIONS			AllFunctions);


Err	CalibrateSmileAndParam(	ALLCALININSTRUMENTSDLM	AllInst,
									int						*DoSmileCalib,
									void					*InstConst,
									void					*GlobalConst,
									void					*Model,							
														
									CALIBGEN_PARAMS		CalibConsts,

									Err					(*GetLimitLastParam)(		void				*Model,
																					CALIBGEN_PARAMS	CalibConsts,
																					void				*GlobalConst,

																					int					index_param,
																					double				*last_param,
																					double				*limit_down,
																					double				*limit_up),

									Err					(*GetFirstGuess)(			void			*Model,
																					void			*GlobalConst,
																					int				index_param,
																					double			target,
																					double			*param1),

									Err					(*GetSecondGuess)(			void			*Model,
																					void			*GlobalConst,
																					int				index_param,
																					double			param1,
																					double			price1,
																					double			target,
																					double			*param2),
									
									Err					(*BumpParam)(					void			*Model,
																					void			*GlobalConst,
																					int				index_param,
																					double			param),
														
									Err					(*UpdateConstsAfterParam) (	void			*Inst,
																					void			*InstConst,
																					void			*GlobalConst,
																					void			*Model),
														
														
									Err					(*PriceInst)(				void			*Inst,
																					void			*InstConst,
																					void			*GlobalConst,
																					void			*Model,	
																					Err				(*UpdateConstsAfterParam) (	void			*Inst,
																																void			*InstConst,
																																void			*GlobalConst,
																																void			*Model),
																					double			*InstPrice),

									Err					(*SetParam)(				void			*Model,
																					void			*GlobalConst,
																					int				index_param,
																					double			param),

									Err					(*ExtrapolParam)(			void			*Model,
																					void			*GlobalConst,
																					int				LastParamIndex));

#endif