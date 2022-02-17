 
#ifndef FX3FBETADLMCALIBRATION_H
#define FX3FBETADLMCALIBRATION_H

#include "srt_h_all.h"
#include "DiagCalibGen.h"

typedef struct
{
	int							iNbOpt;
	void						**cAllInst;
	FxBetaDLM_InstPrecalc		**cAllInstPrecalc;
								
	int							iCalibSmile;
	int							iSmileIndex;
								
	long						*lSigIndex;
								
	FxBetaDLM_ModelPrecalc		*cPrecalc;
	FxBetaDLM_Hermite			*cHermite;
								
	CALIBGEN_PARAMS			CalibParams;
	CALIBGEN_PARAMS			CalibParamsForBeta;
								
	CalibFunctions				*CalibFunctionsForVol;
	CalibFunctions				*CalibFunctionsForBeta;

	FxBetaDLM_OptNumerParams	*cNumParams;

	double						dMinFact;
	double						dMaxFact;

	double						first_guess;
	double						first_guess_for_smile;

	Err	(*ExtrapolAfterSetBeta)(	void			*Model,
									void			*GlobalConst,
									int				LastParamIndex);


} FxBetaDLM_CalibConst, FXBETADLM_CALIBCONST;

void FxBetaDLM_Free_CalibConst(FxBetaDLM_CalibConst	*CalibConst);

Err FxBetaDLM_Fill_CalibConst(	int							iNbOpt,
								void						**cAllInst,						
								int							iCalibSmile,
								int							iSmileIndex,
								FxBetaDLM_model				*cModel,								
								FxBetaDLM_OptNumerParams	*cNumParams,
								FxBetaDLM_Hermite			*cHermite,
								FxBetaDLM_CalibConst		*CalibConst);

Err FxBetaDLM_Fill_ModelPrecalc(FxBetaDLM_model			*model,
								FxBetaDLM_CalibConst	*CalibConst);

Err FxBetaDLM_Fill_InstPrecalc(	FxBetaDLM_FxOptInst		*Inst,
								FxBetaDLM_model			*model,
								FxBetaDLM_InstPrecalc	*CalibConst);

Err	FxBetaDLM_Calculate_AllConst(FxBetaDLM_FxOptInst		*Inst,
								 FxBetaDLM_model			*model,
								 FxBetaDLM_OptNumerParams	*NumParams,
								 FxBetaDLM_InstPrecalc		*InstConst);

Err FxBetaDLM_Update_InstPrecalc_FromMoment(FxBetaDLM_FxOptInst			*Inst,
											FxBetaDLM_model				*model,
											FxBetaDLM_OptNumerParams	*NumParams,
											FxBetaDLM_InstPrecalc		*InstConst);

Err FxBetaDLM_Update_InstPrecalc_FromX0(double						X0,
										FxBetaDLM_FxOptInst			*Inst,
										FxBetaDLM_model				*model,
										FxBetaDLM_OptNumerParams	*NumParams,
										FxBetaDLM_InstPrecalc		*InstConst);

Err FxBetaDLM_Update_ModelPrecalc(	int							i,
									FxBetaDLM_model				*model,
									FxBetaDLM_CalibConst		*CalibConst);

Err FxBetaDLM_Update_InstPrecalc(	FxBetaDLM_FxOptInst		*Inst,
									FxBetaDLM_model			*model,
									FxBetaDLM_CalibConst	*CalibConst,
									FxBetaDLM_InstPrecalc	*InstConst);

Err FxBetaDLM_Calibration(	/* Model informations */
							char	*undname_fx,
							char	*undname_dom,
							char	*undname_for,
							double	tstar,
							double	B0,
							double	*C0,
							double	alpha,
							double	lambda,
							double	spot_fx,

							/*	Correlations	*/
							int		nb_3F_corr,
							long	*date_3F_corr,					
							double	*dom_for_3F_corr,
							double	*dom_fx_3F_corr,
							double	*for_fx_3F_corr,
							
							/* Volatility informations */
							int		nb_opt,
							long	*opt_settlmt_date,
							long	*opt_fixing_date,							
							double	*opt_strikes,
							double	*opt_bs_vols,

							/* Smile informations */
							int		calib_smile,
							long	smile_settlmt_date,								
							double	smile_strikes[2],
							double	smile_bs_vols[2],

							/* Numerical parameters */
							FxBetaDLM_OptNumerParams	*NumParams,
								
							/* Output */
							double		**fx_vols,
							int			*nb_DLM_corr,
							long		**date_DLM_corr,
							double		**dom_for_DLM_corr,
							double		**dom_fx_DLM_corr,
							double		**for_fx_DLM_corr,
							double		**dom_fx_DLM_corr_down,
							double		**for_fx_DLM_corr_down);

Err FxBetaDLM_CalibrationModel(	/* Volatility informations */
								int							nb_opt,
								long						*opt_settlmt_date,
								long						*opt_fixing_date,							
								double						*opt_strikes,
								double						*opt_bs_vols,

								/* Smile informations */
								int							calib_smile,
								long						smile_settlmt_date,								
								double						smile_strikes[2],
								double						smile_bs_vols[2],								

								/* Model informations */
								FxBetaDLM_model				*model,
								FxBetaDLM_OptNumerParams	*NumParams,
								FxBetaDLM_Hermite			*hermite);

Err FxBetaDLM_GetTargetVol(	void				*Inst_,
							void				*GlobalConst,
							void				*Model,
							CALIBGEN_PARAMS	CalibConsts,
							double				*target);

Err FxBetaDLM_GetFirstGuessVol(	void			*Model_,
								void			*GlobalConst,
								int				index_param,
								double			target,
								double			*vol1);

Err FxBetaDLM_GetSecondGuessVol(	void			*Model_,
								void			*GlobalConst,
								int				vol_index,
								double			vol1,
								double			price1,
								double			target,
								double			*vol2);

Err	FxBetaDLM_GetLimitAndLastVol(	void				*Model_,
										CALIBGEN_PARAMS	CalibParams,
										void				*GlobalConst,

										int					vol_index,
										double				*last_vol,
										double				*limit_down,
										double				*limit_up);

Err FxBetaDLM_BumpVol(	void			*Model_,
						void			*GlobalConst,
						int				vol_index,
						double			vol);

Err FxBetaDLM_SetVol(	void				*Model_,
						CALIBGEN_PARAMS	CalibConsts,
						void				*GlobalConst,
						int					vol_index,
						double				vol);

Err	FxBetaDLM_ExtrapolVol(	void	*Model_,						
							void	*GlobalConst,
							int		last_vol_index);

Err	FxBetaDLM_ExtrapolVol_WhenBetaCalib(	void	*Model_,						
											void	*GlobalConst,
											int		last_vol_index);

Err	FxBetaDLM_UpdateConstsAfterVol(	void				*Inst_,
									void				*InstConst_,
									void				*GlobalConst,
									void				*Model_,
									CALIBGEN_PARAMS	CalibConsts);

Err	FxBetaDLM_PriceInstVol(	void				*Inst_,
							void				*InstConst_,
							void				*GlobalConst,
							void				*Model_,
							double				*InstPrice);

Err FxBetaDLM_GetTargetBeta(	void				*Inst_,
							void				*GlobalConst,
							void				*Model,
							CALIBGEN_PARAMS	CalibConsts,
							double				*target);

Err FxBetaDLM_GetFirstGuessBeta(	void			*Model_,
								void			*GlobalConst,
								int				index_param,
								double			target,
								double			*beta1);

Err FxBetaDLM_GetSecondGuessBeta(	void			*Model_,
								void			*GlobalConst,
								int				beta_index,
								double			beta1,
								double			price1,
								double			target,
								double			*beta2);

Err	FxBetaDLM_GetLimitAndLastBeta(	void				*Model_,
										CALIBGEN_PARAMS	CalibParams,
										void				*GlobalConst,

										int					beta_index,
										double				*last_beta,
										double				*limit_down,
										double				*limit_up);

Err FxBetaDLM_BumpBeta(	void			*Model_,
						void			*GlobalConst,
						int				beta_index,
						double			beta);

Err FxBetaDLM_SetBeta(	void				*Model_,
						CALIBGEN_PARAMS	CalibConsts,
						void				*GlobalConst,
						int					beta_index,
						double				beta);

Err	FxBetaDLM_ExtrapolBeta(	void	*Model_,						
							void	*GlobalConst,
							int		last_beta_index);

Err	FxBetaDLM_UpdateConstsAfterBeta(void				*Inst_,
									void				*InstConst_,
									void				*GlobalConst,
									void				*Model_,
									CALIBGEN_PARAMS	CalibConsts);

Err	FxBetaDLM_PriceInstBeta(void				*Inst_,
							void				*InstConst_,
							void				*GlobalConst,
							void				*Model_,
							double				*InstPrice);

#endif