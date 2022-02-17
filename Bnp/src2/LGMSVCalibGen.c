
#include "math.h"
#include "opfnctns.h"
#include "LGMSVClosedForm.h"
#include "LGMSVClosedFormApprox.h"
#include "Fx3FUtils.h"
#include "LGMSVPDE.h"
#include "LGMSVCalibApprox.h"
#include "DiagCalibDLM.h"
#include "DiagCalibDLMSV.h"
#include "DiagCalibGen.h"
#include "LGMSVcalibGen.h"

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to diagonal
		with lambda calibration */


Err	GetTargetVol_LGMSV(	void				*Inst_,
						void				*GlobalParam,
						void				*Model,
												
						CALIBGEN_PARAMS	CalibConsts,
						double				*target)
{
LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

	*target = Inst->dPrice[0];

	if (AllParams->long_param->vega_prec)
	{
		*target /= Inst->dVega[0];
	}

	return NULL;
}


Err	GetFirstGuessVol_LGMSV(void	*Model,
						void	*GlobalParam,
						int		vol_index,
						double	target,
						double	*vol1)
{
LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	*vol1 = pow(model->dSigma[AllParams->lSigLongIndex[vol_index]], 2);

	return NULL;
}

Err	GetSecondGuessVol_LGMSV(	void	*Model,
							void	*GlobalParam,
							int		vol_index,
							double	var1,
							double	price1,
							double	target,
							double	*var2)
{
LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

double	t1, t2, dt;
double	cum_var_temp, cum_var_sv;

	if (vol_index == 0)
	{		
		t1 = 0.0;
		cum_var_temp = 0.0;
		cum_var_sv = 0.0;
	}
	else
	{	
		t1 = model->dPWTime[AllParams->lSigLongIndex[vol_index-1]];
		cum_var_sv = AllParams->EquiLGM->cum_var_sv[AllParams->lSigLongIndex[vol_index-1]];
	}
		
	t2 = model->dPWTime[AllParams->lSigLongIndex[vol_index]];

	dt = t2 - t1;

	cum_var_temp = cum_var_sv + var1 * dt;

	*var2 = (target * target / price1 / price1 * cum_var_temp - cum_var_sv) / dt;

	if (*var2 < 0.0)
	{
		*var2 = 1.0E-08;
	}

	return NULL;
}

Err	BumpVol_LGMSV(	void			*Model,
					void			*GlobalParam,
					int				vol_index,
					double			var)
{
Err	err = NULL;

int		i;
int		last_index, new_index;
double	vol;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (vol_index == 0)
	{
		last_index = 0;
	}
	else
	{
		last_index = AllParams->lSigLongIndex[vol_index-1] + 1;
	}
	
	new_index = AllParams->lSigLongIndex[vol_index];

	vol = sqrt(var);

	for (i=last_index; i<=new_index; i++)
	{
		model->dSigma[i] = vol;
	}

	return err;
}

Err	SetVol_LGMSV(	void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					vol_index,
					double				var)
{
Err	err = NULL;

int		i;
int		last_index, new_index;
double	t1, t2, dt;
double	fact_exp, sig_lgm2;
double	vol, cum_var_sv, cum_var_lgm;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (vol_index == 0)
	{
		last_index = -1;
		t1 = 0.0;
		cum_var_sv = 0.0;
		cum_var_lgm = 0.0;
	}
	else
	{
		last_index = AllParams->lSigLongIndex[vol_index-1];
		t1 = model->dPWTime[last_index];
		cum_var_sv = AllParams->EquiLGM->cum_var_sv[last_index];
		cum_var_lgm = AllParams->EquiLGM->cum_var_lgm[last_index];
	}
	
	new_index = AllParams->lSigLongIndex[vol_index];
	t2 = model->dPWTime[new_index];

	dt = t2 - t1;

	vol = sqrt(var);

	for (i=last_index + 1; i<=new_index; i++)
	{
		model->dSigma[i] = vol;
	}

	AllParams->EquiLGM->cum_var_sv[new_index] = cum_var_sv + var * dt;		

	fact_exp = sqrt((AllParams->EquiLGM->exp_fact[new_index+1] - AllParams->EquiLGM->exp_fact[last_index+1]) 
		/ (2.0 * model->dLambdaX) / dt);

	sig_lgm2 = var / fact_exp / fact_exp;
	AllParams->EquiLGM->cum_var_lgm[new_index] = cum_var_lgm + sig_lgm2 * dt;
	AllParams->EquiLGM->lgm_vol[new_index] = sqrt(sig_lgm2);
	
	return err;
}

Err	GetLimitAndLastVol_LGMSV(	void				*Model,
								CALIBGEN_PARAMS	CalibParams,
								void				*GlobalParam,

								int					vol_index,
								double				*last_vol,
								double				*limit_down,
								double				*limit_up)
{
Err	err = NULL;

double	t1, t2, dt;
int		last_index, new_index;
double	fact_exp;
double	sig_limit_down_lgm, sig_limit_up_lgm;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (vol_index == 0)
	{
		*last_vol = 0.00001;
		*limit_down = 0.000000001;
		*limit_up = 10000.0;
		return err;
	}
			
	last_index = AllParams->lSigLongIndex[vol_index-1];
	new_index = AllParams->lSigLongIndex[vol_index];

	t1 = model->dPWTime[last_index];	
	t2 = model->dPWTime[new_index];
	dt = t2 - t1;

	sig_limit_down_lgm = sqrt(AllParams->EquiLGM->cum_var_lgm[last_index] * AllParams->MinFact / t1);
	sig_limit_up_lgm = sqrt(AllParams->EquiLGM->cum_var_lgm[last_index] / AllParams->MaxFact / t1);

	fact_exp = sqrt((AllParams->EquiLGM->exp_fact[new_index+1] - AllParams->EquiLGM->exp_fact[last_index+1]) 
		/ (2.0 * model->dLambdaX) / dt);
	
	*last_vol = pow(AllParams->EquiLGM->lgm_vol[last_index] * fact_exp, 2);
	*limit_down = pow(sig_limit_down_lgm * fact_exp, 2);
	*limit_up = pow(sig_limit_up_lgm * fact_exp, 2);

	return err;
}

Err	ExtrapolVol_LGMSV(	void				*Model,						
						void				*GlobalParam,
						int					last_vol_index)
{
Err	err = NULL;

double	dt;
int		last_index;
double	fact_exp;
int		i;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	last_index = AllParams->lSigLongIndex[last_vol_index];
	
	for (i=last_index + 1; i<model->iNbPWTime; i++)
	{
		dt = model->dPWTime[i] - model->dPWTime[i-1];

		fact_exp = sqrt((AllParams->EquiLGM->exp_fact[i+1] - AllParams->EquiLGM->exp_fact[i]) 
						/ (2.0 * model->dLambdaX) / dt);
	
		model->dSigma[i] = AllParams->EquiLGM->lgm_vol[last_index] * fact_exp;
	}

	return err;
}

Err	UpdateParamsAfterVol_LGMSV (void				*Inst_,
								void				*InstParam,
								void				*GlobalParam,
								void				*Model,
								CALIBGEN_PARAMS	CalibConsts)
{
Err	err = NULL;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS		HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM) Inst_;

	err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(	Inst,
														HestonParam->sCalibCpnSchedule,
														model,
														HestonParam->HestonInst);

	return err;
}

Err	PriceInstVol_LGMSV(	void					*Inst_,
						void					*InstParam,
						void					*GlobalParam,
						void					*Model,
						double					*InstPrice)
{
Err	err = NULL;

int	iNbStrike;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS		HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM) Inst_;
	
	iNbStrike = HestonParam->NumerInst->iNbStrike;
	HestonParam->NumerInst->iNbStrike = 1;

	LGMSVClosedFormApprox_struct(	model,
									Inst,
									HestonParam->HestonInst,
									AllParams->PricingConst,
									HestonParam->NumerInst,
									Inst->dModelPrice);

	HestonParam->NumerInst->iNbStrike = iNbStrike;	

	if (AllParams->long_param->vega_prec)
	{
		InstPrice[0] = Inst->dModelPrice[0] / Inst->dVega[0];
	}
	else
	{
		InstPrice[0] = Inst->dModelPrice[0];
	}

	return err;
}


Err	GetTargetRho_LGMSV(		void				*Inst_,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM) Inst_;

	if (AllParams->calib_params->calib_alpha || AllParams->calib_params->calib_rr_bt)
	{
		*target = Inst->dPrice[1] - Inst->dPrice[2];
	}
	else
	{
		*target = Inst->dPrice[1];
	}

	return NULL;
}

Err	GetFirstGuessRho_LGMSV(	void			*Model,
							void			*GlobalConst,
							int				index_smile,
							double			target,
							double			*rho1)
{
	LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
	LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);

	if (index_smile == 0)
	{
		*rho1 = model->dRho[AllParams->lSigLongIndex[index_smile]];
	}
	else
	{
		*rho1 = model->dRho[AllParams->lSigLongIndex[index_smile - 1]];
	}

	return NULL;	
}

Err	GetSecondGuessRho_LGMSV(void			*Model,
							void			*GlobalParam,
							int				index_smile,
							double			rho1,
							double			price1,							
							double			target,
							double			*rho2)
{
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (fabs(target - price1) < AllParams->CalibParamsRho->Precision)
	{
		*rho2 = rho1 * 1.001;
	}
	else
	if (target > price1)
	{
		*rho2 = rho1 + 0.2;
	}
	else
	{
		*rho2 = rho1 - 0.2;
	}

	return NULL;
}

Err	GetLimitAndLastRho_LGMSV(	void				*Model,
								CALIBGEN_PARAMS	CalibConsts,
								void				*GlobalParam,
								int					index_smile,
								double				*last_rho,
								double				*limit_down,
								double				*limit_up)
{
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (index_smile == 0)
	{
		*last_rho = 0.0;
	}
	else
	{
		*last_rho = model->dRho[AllParams->lSmileLongIndex[index_smile-1]];
	}

	*limit_down = -0.9;
	*limit_up = 0.9;

	return NULL;
}

Err	BumpRho_LGMSV(	void			*Model,
					void			*GlobalParam,
					int				smile_param,
					double			rho)
{
Err	err = NULL;

int		i;
int		last_index, new_index;
double	rho1, rho2;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);


	if (smile_param == 0)
	{
		last_index = 0;
	}
	else
	{
		last_index = AllParams->lSmileLongIndex[smile_param-1] + 1;
	}
	
	new_index = AllParams->lSmileLongIndex[smile_param];

	if (model->iOne2F == 2 && AllParams->calib_params->onefac_rho)
	{
		err = LGMSV_Get_Rho_From_RhoTarget(	rho,
											AllParams->calib_params->rho_mat1,
											AllParams->calib_params->rho_mat2,
											model->dInitLGMAlpha,
											model->dLGMGamma,
											model->dInitLGMRho,
											&rho1,
											&rho2);

		if (err) return err;

		for (i=last_index; i<=new_index; i++)
		{
			model->dRho[i] = rho1;
			model->dRho2[i] = rho2;
		}
	}
	else
	{
		for (i=last_index; i<=new_index; i++)
		{
			model->dRho[i] = rho;
		}
	}		

	return err;
}

Err	SetRho_LGMSV(	void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					smile_param,
					double				rho)
{
Err	err = NULL;

int		i;
int		last_index, new_index;
double	rho1, rho2;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
	
	if (CalibConsts->do_calib)
	{
		if (smile_param == 0)
		{
			last_index = 0;
		}
		else
		{
			last_index = AllParams->lSmileLongIndex[smile_param-1] + 1;
		}
		
		new_index = AllParams->lSmileLongIndex[smile_param];
		
		if (model->iOne2F == 2 && AllParams->calib_params->onefac_rho)
		{
			err = LGMSV_Get_Rho_From_RhoTarget(	rho,
												AllParams->calib_params->rho_mat1,
												AllParams->calib_params->rho_mat2,
												model->dInitLGMAlpha,
												model->dLGMGamma,
												model->dInitLGMRho,
												&rho1,
												&rho2);

			if (err) return err;

			for (i=last_index; i<=new_index; i++)
			{
				model->dRho[i] = rho1;
				model->dRho2[i] = rho2;
			}
		}
		else
		{
			for (i=last_index; i<=new_index; i++)
			{
				model->dRho[i] = rho;
			}
		}
	}
	
	AllParams->index_vol += 1;

	return err;
}

Err	PriceInstRho_LGMSV(	void			*Inst_,
						void			*InstParam,
						void			*GlobalParam,
						void			*Model,
						double			*InstPrice)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
					
	LGMSVClosedFormApprox_struct(	model,
									Inst,
									HestonParam->HestonInst,
									AllParams->PricingConst,
									HestonParam->NumerInst,
									Inst->dModelPrice);

	if (AllParams->calib_params->calib_alpha || AllParams->calib_params->calib_rr_bt)
	{
		*InstPrice = Inst->dModelPrice[1] - Inst->dModelPrice[2];
	}
	else
	{
		*InstPrice = Inst->dModelPrice[1];
	}

	return err;
}

Err	ExtrapolRho_LGMSV(	void				*Model,						
						void				*GlobalParam,
						int					last_smile_index)
{
Err	err = NULL;

int		last_index;
int		i;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	last_index = AllParams->lSmileLongIndex[last_smile_index];
	
	if (model->iOne2F == 2 && AllParams->calib_params->onefac_rho)
	{
		for (i=last_index + 1; i<model->iNbPWTime; i++)
		{		
			model->dRho[i] = model->dRho[last_index];
			model->dRho2[i] = model->dRho2[last_index];
		}
	}
	else
	{
		for (i=last_index + 1; i<model->iNbPWTime; i++)
		{		
			model->dRho[i] = model->dRho[last_index];
		}
	}

	return err;
}

Err	UpdateParamsAfterRho_LGMSV (	void					*Inst_,
									void					*InstParam,
									void					*GlobalParam,
									void					*Model,
									CALIBGEN_PARAMS		CalibConsts)
{
Err	err = NULL;
double	var, last_var, alpha, last_alpha, limit_down, limit_up;
int		success;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

	if (AllParams->calib_params->calib_alpha)
	{
		err = AllParams->CalibFunctionsForAlpha
			.GetLimitAndLastParam(	Model,
									AllParams->CalibParamsAlpha,
									AllParams,
									AllParams->index_vol,
									&last_alpha,
									&limit_down,
									&limit_up);

		if (err) goto FREE_RETURN;

		err = CalibrateNextParam(	Inst,
									HestonParam,
									AllParams,
									Model,
									AllParams->CalibParamsAlpha,
									AllParams->index_vol,
									last_alpha,
									limit_down,
									limit_up,
									&(AllParams->CalibFunctionsForAlpha),
									&alpha,
									&success);

		if (err) goto FREE_RETURN;

		err = AllParams->CalibFunctionsForAlpha
								.SetParam(	Model,
											AllParams->CalibParamsAlpha,
											AllParams,
											AllParams->index_vol,
											alpha);

		if (err) goto FREE_RETURN;
	}
	else	
	{
		err = AllParams->CalibFunctionsForVol
			.GetLimitAndLastParam(	Model,
									AllParams->CalibParams,
									AllParams,
									AllParams->index_vol,
									&last_var,
									&limit_down,
									&limit_up);

		if (err) goto FREE_RETURN;

		err = CalibrateNextParam(	Inst,
									HestonParam,
									AllParams,
									Model,
									AllParams->CalibParams,
									AllParams->index_vol,
									last_var,
									limit_down,
									limit_up,
									&(AllParams->CalibFunctionsForVol),
									&var,
									&success);

		if (err) goto FREE_RETURN;

		err = AllParams->CalibFunctionsForVol
								.SetParam(	Model,
											AllParams->CalibParams,
											AllParams,
											AllParams->index_vol,
											var);

		if (err) goto FREE_RETURN;		
	}

FREE_RETURN:

	return err;
}

Err	GetTargetAlpha_LGMSV(	void				*Inst_,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM) Inst_;

	if (AllParams->calib_params->calib_rho || AllParams->calib_params->calib_rr_bt)
	{
		*target = Inst->dPrice[2] + Inst->dPrice[1] - 2.0 * Inst->dPrice[0];
	}
	else
	{
		*target = Inst->dPrice[1];
	}

	return NULL;
}

Err	GetFirstGuessAlpha_LGMSV(	void			*Model,
								void			*GlobalConst,
								int				index_smile,
								double			target,
								double			*alpha1)
{
	LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
	LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);

	if (index_smile == 0)
	{
		*alpha1 = model->dAlpha[AllParams->lSigLongIndex[index_smile]];
	}
	else
	{
		*alpha1 = model->dAlpha[AllParams->lSigLongIndex[index_smile - 1]];
	}

	return NULL;
}

Err	GetSecondGuessAlpha_LGMSV(void			*Model,
							void			*GlobalParam,
							int				index_smile,
							double			alpha1,
							double			price1,							
							double			target,
							double			*alpha2)
{
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (fabs(target - price1) < AllParams->CalibParamsAlpha->Precision)
	{
		*alpha2 = alpha1 * 1.001;
	}
	else
	if (target > price1)
	{
		*alpha2 = alpha1 + 0.2;
	}
	else
	{
		*alpha2 = alpha1 - 0.2;
	}

	return NULL;
}

Err	UpdateParamsAfterAlpha_LGMSV (	void					*Inst_,
									void					*InstParam,
									void					*GlobalParam,
									void					*Model,
									CALIBGEN_PARAMS		CalibConsts)
{
Err	err = NULL;
double	var, last_var, limit_down, limit_up;
int		success;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

	err = AllParams->CalibFunctionsForVol
		.GetLimitAndLastParam(	Model,
								AllParams->CalibParams,
								AllParams,
								AllParams->index_vol,
								&last_var,
								&limit_down,
								&limit_up);

	if (err) goto FREE_RETURN;

	err = CalibrateNextParam(	Inst,
								HestonParam,
								AllParams,
								Model,
								AllParams->CalibParams,
								AllParams->index_vol,
								last_var,
								limit_down,
								limit_up,
								&(AllParams->CalibFunctionsForVol),
								&var,
								&success);

	if (err) goto FREE_RETURN;

	err = AllParams->CalibFunctionsForVol
							.SetParam(	Model,
										AllParams->CalibParams,
										AllParams,
										AllParams->index_vol,
										var);

	if (err) goto FREE_RETURN;

FREE_RETURN:

	return err;
}

Err	PriceInstAlpha_LGMSV(	void			*Inst_,
							void			*InstParam,
							void			*GlobalParam,
							void			*Model,
							double			*InstPrice)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
					
	LGMSVClosedFormApprox_struct(	model,
									Inst,
									HestonParam->HestonInst,
									AllParams->PricingConst,
									HestonParam->NumerInst,
									Inst->dModelPrice);

	if (AllParams->calib_params->calib_rho || AllParams->calib_params->calib_rr_bt)
	{
		*InstPrice = Inst->dModelPrice[2] + Inst->dModelPrice[1] - 2.0 * Inst->dModelPrice[0];
	}
	else
	{
		*InstPrice = Inst->dModelPrice[1];
	}

	return err;
}

Err	GetLimitAndLastAlpha_LGMSV(	void				*Model,
								CALIBGEN_PARAMS	CalibConsts,
								void				*GlobalParam,
								int					index_smile,
								double				*last_alpha,
								double				*limit_down,
								double				*limit_up)
{
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (index_smile == 0)
	{
		*last_alpha = 0.0;
	}
	else
	{
		*last_alpha = model->dAlpha[AllParams->lSmileLongIndex[index_smile-1]];
	}

	*limit_down = 0.1;
	*limit_up = 3.0;

	return NULL;
}

Err	BumpAlpha_LGMSV(	void			*Model,
					void			*GlobalParam,
					int				smile_param,
					double			alpha)
{
Err	err = NULL;

int		i;
int		last_index, new_index;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (smile_param == 0)
	{
		last_index = 0;
	}
	else
	{
		last_index = AllParams->lSmileLongIndex[smile_param-1] + 1;
	}
	
	new_index = AllParams->lSmileLongIndex[smile_param];
	

	for (i=last_index; i<=new_index; i++)
	{
		model->dAlpha[i] = alpha;
	}

	return err;
}

Err	SetAlpha_LGMSV(	void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					smile_param,
					double				alpha)
{
Err	err = NULL;

int		i;
int		last_index, new_index;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (CalibConsts->do_calib)
	{
		if (smile_param == 0)
		{
			last_index = 0;
		}
		else
		{
			last_index = AllParams->lSmileLongIndex[smile_param-1] + 1;
		}
		
		new_index = AllParams->lSmileLongIndex[smile_param];
		

		for (i=last_index; i<=new_index; i++)
		{
			model->dAlpha[i] = alpha;
		}
	}
	
	return err;
}

Err	ExtrapolAlpha_LGMSV(	void				*Model,						
						void				*GlobalParam,
						int					last_smile_index)
{
Err	err = NULL;

int		last_index;
int		i;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	last_index = AllParams->lSmileLongIndex[last_smile_index];
	
	for (i=last_index + 1; i<model->iNbPWTime; i++)
	{		
		model->dAlpha[i] = model->dAlpha[last_index];
	}

	return err;
}

Err	GetTargetLambda_LGMSV(	void				*Inst_,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
CalibInstrumentDLM	*NextInst;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
double			total;

	NextInst = Inst;
	total = 0.0;

	while (NextInst != NULL)
	{
		total += NextInst->dPrice[0] * NextInst->dWeight;
		NextInst = NextInst->NextInst;
	}

	*target = total;
	
	return NULL;
}

Err	GetFirstGuessLambda_LGMSV(	void			*Model,
							void			*GlobalConst,
							int				index_lambda,
							double			target,
							double			*lambda1)
{
	LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
	LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
	
	*lambda1 = model->dLambdaX;	

	return NULL;	
}

Err	GetSecondGuessLambda_LGMSV(void			*Model,
							void			*GlobalParam,
							int				index_lambda,
							double			lambda1,
							double			price1,							
							double			target,
							double			*lambda2)
{
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (fabs(target - price1) < AllParams->CalibParamsLambda->Precision)
	{
		*lambda2 = lambda1 * 1.001;
	}
	else
	{
		if (target > price1)
		{
			*lambda2 = lambda1 + AllParams->sens_lambda * 0.005;
		}
		else
		{
			*lambda2 = lambda1 - AllParams->sens_lambda * 0.005;
		}
	}

	return NULL;
}

Err	GetLimitAndLastLambda_LGMSV(void				*Model,
								CALIBGEN_PARAMS	CalibConsts,
								void				*GlobalParam,
								int					index_lambda,
								double				*last_lambda,
								double				*limit_down,
								double				*limit_up)
{
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
	
	*last_lambda = model->dLambdaX;
	*limit_down = AllParams->short_param->lambda_min;
	*limit_up = AllParams->short_param->lambda_max;

	return NULL;
}

Err	BumpLambda_LGMSV(void			*Model,
					void			*GlobalParam,
					int				index_lambda,
					double			lambda)
{
Err	err = NULL;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
	
	model->dLambdaX = lambda;
	model->dTau = 1.0 / model->dLambdaX;
	
	if (model->iOne2F == 2)
	{
		model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;
		model->dTau2 = 1.0 / model->dLambdaX2;

		ConvertAlphaRho_LGM_to_LGMSV(	model->iNbPWTime,
										model->dPWTime,
										model->dInitTStar,
										model->dLambdaX,
										model->dInitLGMAlpha,
										model->dLGMGamma,
										model->dInitLGMRho,
										model->dLGMAlpha,
										model->dLGMRho);
	}


	
	return err;
}

Err	SetLambda_LGMSV(void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					index_lambda,
					double				lambda)
{
Err	err = NULL;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
	
	if (CalibConsts->do_calib)
	{
		model->dLambdaX = lambda;
		model->dTau = 1.0 / model->dLambdaX;
	
		if (model->iOne2F == 2)
		{
			model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;
			model->dTau2 = 1.0 / model->dLambdaX2;		

			ConvertAlphaRho_LGM_to_LGMSV(	model->iNbPWTime,
											model->dPWTime,
											model->dInitTStar,
											model->dLambdaX,
											model->dInitLGMAlpha,
											model->dLGMGamma,
											model->dInitLGMRho,
											model->dLGMAlpha,
											model->dLGMRho);
		}
	}
	
	AllParams->index_vol = 0;

	return err;
}

Err	PriceInstLambda_LGMSV(	void				*Inst_,
							void				*InstParam,
							void				*GlobalParam,
							void				*Model,
							double				*InstPrice)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
InstParams			*NextHestonParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

CalibInstrumentDLM	*NextInst;
double				total;

	total = 0.0;
	NextInst = Inst;
	NextHestonParam = HestonParam;

	while (NextInst != NULL)
	{					
		LGMSVClosedFormApprox_struct(	model,
										NextInst,
										NextHestonParam->HestonInst,
										AllParams->PricingConst,
										NextHestonParam->NumerInst,
										NextInst->dModelPrice);

		total += NextInst->dModelPrice[0] * NextInst->dWeight;
		NextInst = NextInst->NextInst;
		NextHestonParam = NextHestonParam->NextInstParams;
	}

	*InstPrice = total;	

	return err;
}

Err	ExtrapolLambda_LGMSV(	void				*Model,						
							void				*GlobalParam,
							int					last_smile_index)
{
	return NULL;
}

Err	UpdateParamsAfterLambda_LGMSV (	void					*Inst_,
									void					*InstParam,
									void					*GlobalParam,
									void					*Model,
									CALIBGEN_PARAMS		CalibConsts)
{
Err	err = NULL;
int		i, j, last_index, nb_inst;
double	time0 = 100;
double	ex_zeta[MAX_CPN];
double	temp, fact;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CalibInstrumentDLM	*NextInst;
InstParams			*NextHestonParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
	
	/* All the precalculations !!! */
	if (AllParams->calib_params->use_lgm_first_guess)
	{
		err = lgmcalibzetalambda_tauts_dlm(
						AllParams->CalibCpnLongSchedule->iNCpn,						
						AllParams->CalibCpnLongSchedule->dCpnTime,
						AllParams->CalibCpnLongSchedule->dCpnDf,
						AllParams->CalibCpnLongSchedule->dCpnCvg,
						
						AllParams->CalibExeLongSchedule->iNExe,
						AllParams->CalibExeLongSchedule->dExeTimes,
						AllParams->CalibExeLongSchedule->iStartCpn,
						AllParams->CalibExeLongSchedule->iEndCpn,
						AllParams->CalibExeLongSchedule->dATMStrike,
						AllParams->CalibExeLongSchedule->dATMPrice,
						AllParams->CalibExeLongSchedule->dATMVega,

						AllParams->CalibCpnShortSchedule->iNCpn,						
						AllParams->CalibCpnShortSchedule->dCpnTime,
						AllParams->CalibCpnShortSchedule->dCpnDf,
						AllParams->CalibCpnShortSchedule->dCpnCvg,
						
						AllParams->CalibExeShortSchedule->iNExe,
						AllParams->CalibExeShortSchedule->dExeTimes,
						AllParams->CalibExeShortSchedule->iStartCpn,
						AllParams->CalibExeShortSchedule->iEndCpn,
						AllParams->CalibExeShortSchedule->dATMStrike,
						NULL,
						AllParams->CalibExeShortSchedule->dATMPrice,
						AllParams->CalibExeShortSchedule->dATMPrice,

						ex_zeta,
						1,
						1,
						&time0,
						&model->dLambdaX,
						model->iOne2F,
						model->dInitLGMAlpha,
						model->dLGMGamma,
						model->dInitLGMRho,
						AllParams->long_param,
						AllParams->short_param,
						NULL,
						&temp);
		 
		fact = exp(model->dLambdaX * model->dInitTStar);
		temp = sqrt( ex_zeta[0] / AllParams->CalibExeLongSchedule->dExeTimes[0]) / fact;

		for (j=0; j<=AllParams->lSigLongIndex[0]; j++)
		{
			model->dSigma[j] = temp;
		}

		last_index = AllParams->lSigLongIndex[0];

		for (i=1; i<AllParams->CalibExeLongSchedule->iNExe; i++)
		{		
			if (ex_zeta[i] > ex_zeta[i-1])
			{
				temp = sqrt ((ex_zeta[i] - ex_zeta[i-1]) / (AllParams->CalibExeLongSchedule->dExeTimes[i] - AllParams->CalibExeLongSchedule->dExeTimes[i-1])) / fact;

				for (j=last_index+1; j<=AllParams->lSigLongIndex[i]; j++)
				{
					model->dSigma[j] = temp;
				}
				
				last_index = AllParams->lSigLongIndex[i];
			}
			else
			{
				smessage ("Diagonal calibration failed at exercise year %.2f - Calibration stopped", AllParams->CalibExeLongSchedule->dExeTimes[i]);
				
				for (i=i; i<AllParams->CalibExeLongSchedule->iNExe; i++)
				{
					for (j=last_index+1; j<=AllParams->lSigLongIndex[i]; j++)
					{
						model->dSigma[j] = temp;
					}
				
					last_index = AllParams->lSigLongIndex[i];
				}
			}
		}
	}

	/* update the long instruments */
	
	for (i=0; i<AllParams->iNbLongInst; i++)
	{
		NextInst = (CALIBINSTRUMENTDLM) (AllParams->AllLongInst[i]);

		err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(	NextInst,
																AllParams->CalibCpnLongSchedule,
																model,
																AllParams->LongInstParams[i]->HestonInst);

		if (err) goto FREE_RETURN;		
	}

	/* update the EQ LGM */
	err = Update_EquiLGM(	model,
							AllParams->EquiLGM);

	AllParams->index_vol = 0;

	if (CalibConsts->do_calib)
	{	
		/* update the short instruments */
		
		NextInst = Inst;
		NextHestonParam = HestonParam;

		while (NextInst != NULL)
		{
			err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(	NextInst,
																	AllParams->CalibCpnShortSchedule,
																	model,
																	NextHestonParam->HestonInst);

			if (err) goto FREE_RETURN;

			NextInst = NextInst->NextInst;
			NextHestonParam = NextHestonParam->NextInstParams;
		}
	}

	nb_inst = AllParams->iNbLongInst;

	if (AllParams->calib_params->calib_flat_smile)
	{
		nb_inst = 1;
	}

	if (AllParams->calib_params->calib_smile_on_prim || nb_inst > 1)
	{
		err = CalibrateParamTS(	0,
								nb_inst-1,
								AllParams->AllLongInst,
								AllParams->LongInstParams,
								AllParams,
								model,
								AllParams->CalibParamsRho,
								&(AllParams->CalibFunctionsForRho));
	}
	else
	{
		err = CalibrateParamTS(	0,
								nb_inst-1,
								AllParams->AllShortInst,
								AllParams->ShortInstParams,
								AllParams,
								model,
								AllParams->CalibParamsRho,
								&(AllParams->CalibFunctionsForRho));
	}

	if (err) goto FREE_RETURN;

	NextInst = Inst;
	NextHestonParam = HestonParam;

	while (NextInst != NULL)
	{
		err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(	NextInst,
															NextHestonParam->sCalibCpnSchedule,
															model,
															NextHestonParam->HestonInst);

		if (err) goto FREE_RETURN;

		NextInst = NextInst->NextInst;
		NextHestonParam = NextHestonParam->NextInstParams;
	}
	
FREE_RETURN:

	return err;
}

Err	GetTargetRho2_LGMSV(	void				*Inst_,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
CalibInstrumentDLM	*NextInst;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
double			total;

	NextInst = Inst;
	total = 0.0;

	while (NextInst != NULL)
	{
		total += NextInst->dPrice[1];

		if (NextInst->iNbStrike > 2)
		{
			total -= NextInst->dPrice[2];
		}

		NextInst = NextInst->NextInst;
	}

	*target = total;
	
	return NULL;
}

Err	GetFirstGuessRho2_LGMSV(void			*Model,
							void			*GlobalConst,
							int				index_rho2,
							double			target,
							double			*rho21)
{
	LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
	LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
	
	*rho21 = model->dRho2[0];

	return NULL;	
}

Err	GetSecondGuessRho2_LGMSV(void			*Model,
							void			*GlobalParam,
							int				index_rho2,
							double			rho21,
							double			price1,							
							double			target,
							double			*rho22)
{
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (fabs(target - price1) < AllParams->CalibParamsRho2->Precision)
	{
		*rho22 = rho21 * 1.0001;
	}
	else
	{
		if (rho21 < 0.9)
		{
			*rho22 = rho21 + 0.05;
		}
		else
		{
			*rho22 = rho21 - 0.05;
		}
	}

	return NULL;
}

Err	GetLimitAndLastRho2_LGMSV(void				*Model,
								CALIBGEN_PARAMS	CalibConsts,
								void				*GlobalParam,
								int					index_rho2,
								double				*last_rho2,
								double				*limit_down,
								double				*limit_up)
{
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	
	*last_rho2 = model->dRho2[0];
	*limit_down = -0.9;
	*limit_up = 0.9;

	return NULL;
}

Err	BumpRho2_LGMSV(	void			*Model,
					void			*GlobalParam,
					int				index_rho2,
					double			rho2)
{
Err	err = NULL;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
int	i;
	
	for (i=0; i<model->iNbPWTime; i++)
	{
		model->dRho2[i] = rho2;
	}
		
	return err;
}

Err	SetRho2_LGMSV(	void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					index_rho2,
					double				rho2)
{
Err	err = NULL;

int	i;
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
	
	if (CalibConsts->do_calib)
	{
		for (i=0; i<model->iNbPWTime; i++)
		{
			model->dRho2[i] = rho2;
		}
	}
	
	AllParams->index_vol = 0;

	return err;
}

Err	PriceInstRho2_LGMSV(	void				*Inst_,
							void				*InstParam,
							void				*GlobalParam,
							void				*Model,
							double				*InstPrice)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
InstParams			*NextHestonParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

CalibInstrumentDLM	*NextInst;
double				total;

	total = 0.0;
	NextInst = Inst;
	NextHestonParam = HestonParam;

	while (NextInst != NULL)
	{					
		LGMSVClosedFormApprox_struct(	model,
										NextInst,
										NextHestonParam->HestonInst,
										AllParams->PricingConst,
										NextHestonParam->NumerInst,
										NextInst->dModelPrice);

		total += NextInst->dModelPrice[1];

		if (NextInst->iNbStrike > 2)
		{
			total -= NextInst->dModelPrice[2];
		}

		NextInst = NextInst->NextInst;
		NextHestonParam = NextHestonParam->NextInstParams;
	}

	*InstPrice = total;	

	return err;
}

Err	ExtrapolRho2_LGMSV(	void				*Model,						
						void				*GlobalParam,
						int					last_smile_index)
{
	return NULL;
}

Err	UpdateParamsAfterRho2_LGMSV (	void					*Inst_,
									void					*InstParam,
									void					*GlobalParam,
									void					*Model,
									CALIBGEN_PARAMS		CalibConsts)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
	
	err = CalibrateParamTS(	0,
							0,
							AllParams->AllShortInst,
							&InstParam,
							AllParams,
							model,
							AllParams->CalibParamsLambda,
							&(AllParams->CalibFunctionsForLambda));

	if (err) goto FREE_RETURN;
	
	
FREE_RETURN:

	return err;
}

Err	GetTargetFlatAlpha_LGMSV(	void				*Inst_,
								void				*GlobalConst,
								void				*Model,
												
								CALIBGEN_PARAMS	CalibConsts,
								double				*target)
{
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
CalibInstrumentDLM	*NextInst;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
double			total;

	NextInst = Inst;
	total = 0.0;

	while (NextInst != NULL)
	{
		if (AllParams->calib_params->calib_rho || AllParams->calib_params->calib_rr_bt)
		{		
			total += NextInst->dPrice[2] + NextInst->dPrice[1] - 2.0 * NextInst->dPrice[0];
		}
		else
		{			
			total += NextInst->dPrice[1];
		}
		
		NextInst = NextInst->NextInst;
	}

	*target = total;
	
	return NULL;
}

Err	GetFirstGuessFlatAlpha_LGMSV(	void			*Model,
							void			*GlobalConst,
							int				index_flat_alpha,
							double			target,
							double			*flat_alpha1)
{
LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
	
	*flat_alpha1 = model->dAlpha[0];	

	return NULL;	
}

Err	GetSecondGuessFlatAlpha_LGMSV(	void			*Model,
									void			*GlobalParam,
									int				index_flat_alpha,
									double			flat_alpha1,
									double			price1,							
									double			target,
									double			*flat_alpha2)
{
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (fabs(target - price1) < AllParams->CalibParamsAlpha->Precision)
	{
		*flat_alpha2 = flat_alpha1 * 1.001;
	}
	else
	if (target > price1)
	{
		*flat_alpha2 = flat_alpha1 + 0.2;
	}
	else
	{
		*flat_alpha2 = flat_alpha1 - 0.2;
	}

	return NULL;
}

Err	GetLimitAndLastFlatAlpha_LGMSV(void				*Model,
								CALIBGEN_PARAMS		CalibConsts,
								void				*GlobalParam,
								int					index_flat_alpha,
								double				*last_flat_alpha,
								double				*limit_down,
								double				*limit_up)
{
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	
	*last_flat_alpha = model->dAlpha[0];
	*limit_down = 0.1;
	*limit_up = AllParams->calib_params->max_calib_alpha * 2.0;

	return NULL;
}

Err	BumpFlatAlpha_LGMSV(void			*Model,
					void			*GlobalParam,
					int				index_flat_alpha,
					double			flat_alpha)
{
Err	err = NULL;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
int	i;
	
	for (i=0; i<model->iNbPWTime; i++)
	{
		model->dAlpha[i] = flat_alpha;
	}
	
	return err;
}

Err	SetFlatAlpha_LGMSV(void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					index_flat_alpha,
					double				flat_alpha)
{
Err	err = NULL;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
int				i;
	
	if (CalibConsts->do_calib)
	{
		for (i=0; i<model->iNbPWTime; i++)
		{
			model->dAlpha[i] = flat_alpha;
		}		
	}
	
	AllParams->index_vol = 0;

	return err;
}

Err	PriceInstFlatAlpha_LGMSV(	void				*Inst_,
								void				*InstParam,
								void				*GlobalParam,
								void				*Model,
								double				*InstPrice)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
InstParams			*NextHestonParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

CalibInstrumentDLM	*NextInst;
double				total;

	total = 0.0;
	NextInst = Inst;
	NextHestonParam = HestonParam;

	while (NextInst != NULL)
	{
		if (!AllParams->calib_params->calib_smile_on_prim)
		{
			/* First Recalculate the Shifted Log params */
			err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(	NextInst,
																NextHestonParam->sCalibCpnSchedule,
																model,
																NextHestonParam->HestonInst);

			if (err) return err;
		}

		LGMSVClosedFormApprox_struct(	model,
										NextInst,
										NextHestonParam->HestonInst,
										AllParams->PricingConst,
										NextHestonParam->NumerInst,
										NextInst->dModelPrice);

		if (AllParams->calib_params->calib_rho || AllParams->calib_params->calib_rr_bt)
		{		
			total += NextInst->dModelPrice[2] + NextInst->dModelPrice[1] - 2.0 * NextInst->dModelPrice[0];
		}
		else
		{			
			total += NextInst->dModelPrice[1];
		}

		NextInst = NextInst->NextInst;
		NextHestonParam = NextHestonParam->NextInstParams;
	}

	*InstPrice = total;	

	return err;
}

Err	ExtrapolFlatAlpha_LGMSV(void				*Model,						
							void				*GlobalParam,
							int					last_smile_index)
{
	return NULL;
}

Err	UpdateParamsAfterFlatAlpha_LGMSV (	void					*Inst_,
										void					*InstParam,
										void					*GlobalParam,
										void					*Model,
										CALIBGEN_PARAMS		CalibConsts)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
		
	if (AllParams->CalibParams->do_calib)
	{
		err = CalibrateParamTS(	0,
								AllParams->iNbLongInst-1,
								AllParams->AllLongInst,
								AllParams->LongInstParams,
								GlobalParam,
								Model,
								AllParams->CalibParams,
								&(AllParams->CalibFunctionsForVol));
	}

	AllParams->index_vol = 0;
		
	return err;
}

Err	GetTargetFlatRho_LGMSV(	void				*Inst_,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
CalibInstrumentDLM	*NextInst;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
double			total;

	NextInst = Inst;
	total = 0.0;

	while (NextInst != NULL)
	{
		if (AllParams->calib_params->calib_alpha || AllParams->calib_params->calib_rr_bt)
		{		
			total += NextInst->dPrice[1] - NextInst->dPrice[2];
		}
		else
		{			
			total += NextInst->dPrice[1];
		}
		
		NextInst = NextInst->NextInst;
	}

	*target = total;
	
	return NULL;
}

Err	GetFirstGuessFlatRho_LGMSV(	void			*Model,
								void			*GlobalConst,
								int				index_flat_rho,
								double			target,
								double			*flat_rho1)
{
	LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
	LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalConst);
	
	if (model->iOne2F == 2 && AllParams->calib_params->onefac_rho)
	{
		*flat_rho1 = model->dFlatOneFactorRho;
	}
	else
	{
		*flat_rho1 = model->dRho[0];
	}

	return NULL;	
}

Err	GetSecondGuessFlatRho_LGMSV(	void			*Model,
									void			*GlobalParam,
									int				index_flat_rho,
									double			flat_rho1,
									double			price1,							
									double			target,
									double			*flat_rho2)
{
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	if (fabs(target - price1) < AllParams->CalibParamsRho->Precision)
	{
		*flat_rho2 = flat_rho1 * 1.001;
	}
	else
	if (target > price1)
	{
		*flat_rho2 = flat_rho1 + 0.2;
	}
	else
	{
		*flat_rho2 = flat_rho1 - 0.2;
	}

	return NULL;
}

Err	GetLimitAndLastFlatRho_LGMSV(void				*Model,
								CALIBGEN_PARAMS	CalibConsts,
								void				*GlobalParam,
								int					index_flat_rho,
								double				*last_flat_rho,
								double				*limit_down,
								double				*limit_up)
{
LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);

	
	*last_flat_rho = model->dRho[0];
	*limit_down = -0.9;
	*limit_up = 0.9;

	return NULL;
}

Err	BumpFlatRho_LGMSV(void			*Model,
					void			*GlobalParam,
					int				index_flat_rho,
					double			flat_rho)
{
Err	err = NULL;

LGMSV_MODEL					model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS				AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
int	i;
double	rho1, rho2;
	
	if (model->iOne2F == 2 && AllParams->calib_params->onefac_rho)
	{
		model->dFlatOneFactorRho = flat_rho;

		err = LGMSV_Get_Rho_From_RhoTarget(	flat_rho,
											AllParams->calib_params->rho_mat1,
											AllParams->calib_params->rho_mat2,
											model->dInitLGMAlpha,
											model->dLGMGamma,
											model->dInitLGMRho,
											&rho1,
											&rho2);

		if (err) return err;

		for (i=0; i<model->iNbPWTime; i++)
		{
			model->dRho[i] = rho1;
			model->dRho2[i] = rho2;
		}
	}
	else
	{
		for (i=0; i<model->iNbPWTime; i++)
		{
			model->dRho[i] = flat_rho;
		}
	}
	
	return err;
}

Err	SetFlatRho_LGMSV(void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					index_flat_rho,
					double				flat_rho)
{
Err	err = NULL;

LGMSV_MODEL		model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS	AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
int				i;
double	rho1, rho2;
	
	if (CalibConsts->do_calib)
	{
		if (model->iOne2F == 2 && AllParams->calib_params->onefac_rho)
		{
			model->dFlatOneFactorRho = flat_rho;

			err = LGMSV_Get_Rho_From_RhoTarget(	flat_rho,
												AllParams->calib_params->rho_mat1,
												AllParams->calib_params->rho_mat2,
												model->dInitLGMAlpha,
												model->dLGMGamma,
												model->dInitLGMRho,
												&rho1,
												&rho2);

			if (err) return err;

			for (i=0; i<model->iNbPWTime; i++)
			{
				model->dRho[i] = rho1;
				model->dRho2[i] = rho2;
			}
		}
		else
		{
			for (i=0; i<model->iNbPWTime; i++)
			{
				model->dRho[i] = flat_rho;
			}
		}		
	}
	
	AllParams->index_vol = 0;

	return err;
}

Err	PriceInstFlatRho_LGMSV(	void				*Inst_,
							void				*InstParam,
							void				*GlobalParam,
							void				*Model,
							double				*InstPrice)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
InstParams			*NextHestonParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;

CalibInstrumentDLM	*NextInst;
double				total;

	total = 0.0;
	NextInst = Inst;
	NextHestonParam = HestonParam;

	while (NextInst != NULL)
	{
		if (!AllParams->calib_params->pricerho_on_alpha || !AllParams->calib_params->calib_alpha)
		{			
			if (!AllParams->calib_params->calib_alpha && !AllParams->calib_params->calib_smile_on_prim)
			{
				/* First Recalculate the Shifted Log params */
				err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(	NextInst,
																	NextHestonParam->sCalibCpnSchedule,
																	model,
																	NextHestonParam->HestonInst);

				if (err) return err;
			}

			LGMSVClosedFormApprox_struct(	model,
											NextInst,
											NextHestonParam->HestonInst,
											AllParams->PricingConst,
											NextHestonParam->NumerInst,
											NextInst->dModelPrice);
		}

		if (AllParams->calib_params->calib_alpha || AllParams->calib_params->calib_rr_bt)
		{		
			total += NextInst->dModelPrice[1] - NextInst->dModelPrice[2];
		}
		else
		{			
			total += NextInst->dModelPrice[1];
		}

		NextInst = NextInst->NextInst;
		NextHestonParam = NextHestonParam->NextInstParams;
	}

	*InstPrice = total;	

	return err;
}

Err	ExtrapolFlatRho_LGMSV(	void				*Model,						
							void				*GlobalParam,
							int					last_smile_index)
{
	return NULL;
}

Err	UpdateParamsAfterFlatRho_LGMSV (	void					*Inst_,
										void					*InstParam,
										void					*GlobalParam,
										void					*Model,
										CALIBGEN_PARAMS		CalibConsts)
{
Err	err = NULL;

LGMSV_MODEL			model = (LGMSV_MODEL) (Model);
LGMSV_ALLPARAMS		AllParams = (LGMSV_ALLPARAMS) (GlobalParam);
INSTPARAMS			HestonParam = (INSTPARAMS) InstParam;
CALIBINSTRUMENTDLM	Inst = (CALIBINSTRUMENTDLM) Inst_;
	
	if (AllParams->calib_params->calib_smile_on_prim)
	{
		err = CalibrateParamTS(	0,
								0,
								AllParams->AllLongInst,
								&InstParam,
								GlobalParam,
								Model,
								AllParams->CalibParamsAlpha,
								&(AllParams->CalibFunctionsForAlpha));
	}
	else
	{
		err = CalibrateParamTS(	0,
								0,
								AllParams->AllShortInst,
								&InstParam,
								GlobalParam,
								Model,
								AllParams->CalibParamsAlpha,
								&(AllParams->CalibFunctionsForAlpha));
	}
				
	AllParams->index_vol = 0;
	
	return err;
}