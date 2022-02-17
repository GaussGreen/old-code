
#include "SABRfwdCalib.h"
#include "math.h"
#include "opfnctns.h"

#define	SABR_fwd_DEFAULT_ALPHA	0.35

/* Free the models											 */
/* ********************************************************* */
Free_SABRfwd_Model(SABR_FWD_MODEL		Model)
{
	if (Model)
	{
		if (Model->matrix) free_dmatrix(Model->matrix, 0, Model->nSimul - 1, 0, 2);
		Model->matrix = NULL;
		
		if (Model->mergeTimes) free(Model->mergeTimes);
		Model->mergeTimes = NULL;

		if (Model->corDF) free(Model->corDF);
		Model->corDF = NULL;

		if (Model->corDFx) free(Model->corDFx);
		Model->corDFx = NULL;

		if (Model->corFFx) free(Model->corFFx);
		Model->corFFx = NULL;

		if (Model->sigDom) free(Model->sigDom);
		Model->sigDom = NULL;

		if (Model->sigFor) free(Model->sigFor);
		Model->sigFor = NULL;

		if (Model->sigFx) free(Model->sigFx);
		Model->sigFx = NULL;
	}
}

/* All the Functions for calibration of RHO on Risk Reversal */
/* ********************************************************* */

Err	SABR_fwd_GetTarget_RR(	void				*Inst_,
							void				*Params_,
							void				*Model_,
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
SABR_FWD_INST		Inst	= (SABR_FWD_INST) Inst_;
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;

	if(Params->solve_on_vol)
		*target = Inst->target_vol_up - Inst->target_vol_down;
	else
		*target = Inst->target_price_up - Inst->target_price_down;

	return NULL;
}

Err	SABR_fwd_GetFirstGuess_RR(	void				*Model_,
								void				*Params_,
								int					index_param,
								double				target,
								double				*param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	*param1 = Model->RhoAlphaFwd;

	return NULL;
}

Err	SABR_fwd_GetSecondGuess_RR(void				*Model_,
								void				*Params_,
								int					index_param,
								double				param1,
								double				price1,
								double				target,
								double				*param2)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;

	if (price1 > target)
	{
		*param2 = param1 - (param1 + 0.99 * Model->AlphaFwd) / 4.0;
	}
	else
	{
		*param2 = param1 + (0.99 * Model->AlphaFwd - param1) / 4.0;
	}	

	return NULL;
}

Err	SABR_fwd_GetLimitAndLastParam_RR(	void				*Model_,
										CALIBGEN_PARAMS	CalibConsts,
										void				*Params_,
										int					index_param,
										double				*last_param,
										double				*limit_down,
										double				*limit_up)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;

	*limit_down = -0.99 * 2.0;
	*limit_up = 0.99 * 2.0;

	return NULL;
}

Err	SABR_fwd_BumpParam_RR(	void				*Model_,
							void				*Params_,
							int					index_param,
							double				param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	Model->RhoAlphaFwd = param1;
	Model->RhoFwd = param1 / Model->AlphaFwd;

	if (Model->RhoFwd > 0.99)
	{
		Model->RhoFwd = 0.99;
		Model->AlphaFwd = Model->RhoAlphaFwd / 0.99;
	}
	else if (Model->RhoFwd < -0.99)
	{
		Model->RhoFwd = -0.99;
		Model->AlphaFwd = -Model->RhoAlphaFwd / 0.99;
	}

	return NULL;
}

Err	SABR_fwd_SetParam_RR(	void				*Model_,
							CALIBGEN_PARAMS	CalibConsts,
							void				*Params_,
							int					index_param,
							double				param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	Model->RhoAlphaFwd = param1;
	Model->RhoFwd = param1 / Model->AlphaFwd;

	if (Model->RhoFwd > 0.99)
	{
		Model->RhoFwd = 0.99;
		Model->AlphaFwd = Model->RhoAlphaFwd / 0.99;
	}
	else if (Model->RhoFwd < -0.99)
	{
		Model->RhoFwd = -0.99;
		Model->AlphaFwd = -Model->RhoAlphaFwd / 0.99;
	}

	return NULL;
}

Err	SABR_fwd_ExtrapolParam_RR(	void				*Model_,
								void				*Params_,
								int					index_param)
{
	return NULL;
}

Err	SABR_fwd_UpdateConstsAfterParam_RR(void				*Inst_,
										void				*InstConst_,
										void				*Params_,
										void				*Model_,
										CALIBGEN_PARAMS	CalibConsts)
{
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;
Err	err = NULL;

	err = CalibrateParamTS(	0,
							0,
							&Inst_,
							&Params_,
							Params_,
							Model_,
							Params->CalibParams_ForAlpha,
							Params->AllFunctions_ForAlpha);

	return err;
}

Err	SABR_fwd_PriceInst_RR(	void				*Inst_,
							void				*InstConst_,
							void				*Params_,
							void				*Model_,
							double				*InstPrice)
{
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;
SABR_FWD_INST		Inst	= (SABR_FWD_INST) Inst_;

	if (Params->solve_on_vol)
		*InstPrice = Inst->vol_up - Inst->vol_down;
	else
		*InstPrice = Inst->price_up - Inst->price_down;

	return NULL;
}

/* All the Functions for calibration of ALPHA on Butterfly */
/* ********************************************************* */

Err	SABR_fwd_GetTarget_BT(	void				*Inst_,
							void				*Params_,
							void				*Model_,
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
SABR_FWD_INST		Inst	= (SABR_FWD_INST) Inst_;
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;

	if (Params->solve_on_vol)
		*target = Inst->target_vol_up + Inst->target_vol_down - 2.0 * Inst->target_vol_atm;
	else
		*target = Inst->target_price_up + Inst->target_price_down - 2.0 * Inst->target_price_atm;

	return NULL;
}

Err	SABR_fwd_GetFirstGuess_BT(	void				*Model_,
								void				*Params_,
								int					index_param,
								double				target,
								double				*param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	*param1 = Model->AlphaFwd * Model->AlphaFwd;

	return NULL;
}

Err	SABR_fwd_GetSecondGuess_BT(void				*Model_,
								void				*Params_,
								int					index_param,
								double				param1,
								double				price1,
								double				target,
								double				*param2)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;

	if (price1 > target)
	{
		*param2 = param1 - 0.2 * 0.2;
	}
	else
	{
		*param2 = param1 + 0.2 * 0.2;
	}	

	return NULL;
}

Err	SABR_fwd_GetLimitAndLastParam_BT(	void				*Model_,
										CALIBGEN_PARAMS	CalibConsts,
										void				*Params_,
										int					index_param,
										double				*last_param,
										double				*limit_down,
										double				*limit_up)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;
		
	*limit_down = max(pow(fabs(Model->RhoAlphaFwd) / 0.99, 2.0), 0.01 * 0.01);
	*limit_up = 2.00 * 2.00;	
	
	return NULL;
}

Err	SABR_fwd_BumpParam_BT(	void				*Model_,
							void				*Params_,
							int					index_param,
							double				param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	Model->AlphaFwd = sqrt(param1);
	Model->RhoFwd = Model->RhoAlphaFwd / Model->AlphaFwd;

	return NULL;
}

Err	SABR_fwd_SetParam_BT(	void				*Model_,
							CALIBGEN_PARAMS	CalibConsts,
							void				*Params_,
							int					index_param,
							double				param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	Model->AlphaFwd = sqrt(param1);
	Model->RhoFwd = Model->RhoAlphaFwd / Model->AlphaFwd;

	return NULL;
}

Err	SABR_fwd_ExtrapolParam_BT(void				*Model_,
								void				*Params_,
								int					index_param)
{
	return NULL;
}

Err	SABR_fwd_UpdateConstsAfterParam_BT(void				*Inst_,
										void				*InstConst_,
										void				*Params_,
										void				*Model_,
										CALIBGEN_PARAMS	CalibConsts)
{
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;
Err	err = NULL;

	err = CalibrateParamTS(	0,
							0,
							&Inst_,
							&Params_,
							Params_,
							Model_,
							Params->CalibParams_ForVol,
							Params->AllFunctions_ForVol);

	return err;
}

Err	SABR_fwd_PriceInst_BT(	void				*Inst_,
							void				*InstConst_,
							void				*Params_,
							void				*Model_,
							double				*InstPrice)
{
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;
SABR_FWD_INST		Inst	= (SABR_FWD_INST) Inst_;
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

double	*strikes = NULL, *vols = NULL;

Err	err = NULL;

	strikes = (double*)calloc(2,sizeof(double));
	vols = (double*)calloc(2,sizeof(double));
	if(!strikes || !vols)
	{
		err = "SABR_fwd_PriceInst_BT: memory allocation failure";
		goto FREE_RETURN;
	}

	strikes[0] = Inst->strike_down;
	strikes[1] = Inst->strike_up;

	/* First Get the Prices */
	if(!Model->SVfudge)
	{
		err = fwdSABRpayoff(	Model->today,
								Model->forwardFixDate,
								Model->forwardFixTime,
								Model->valDate,
								Model->fixTime,
								Model->valTime,
								Model->fwd,
								Model->SigmaFwd,
								Model->AlphaFwd,
								Model->BetaFwd,
								Model->RhoFwd,
								Model->matrix,
								Model->nSimul,
								Model->tstarTime,
								Model->tstarDate,
								Model->dom_yc,
								Model->for_yc,
								Model->mergeTimes,
								Model->mergeNtimes,
								Model->sigDom,
								Model->lda_dom,
								Model->sigFor,
								Model->lda_for,
								Model->mergeTimes,
								Model->sigFx,
								Model->mergeNtimes,
								Model->mergeTimes,
								Model->corDF,
								Model->corDFx,
								Model->corFFx,
								Model->mergeNtimes,
								2,
								strikes,
								Params->solve_on_vol,
								Params->solve_on_vol,
								vols);

		if(err)	goto FREE_RETURN;
	}
	else
	{
		err = fwdSABRpayoffAlphaFudge(	Model->today,
										Model->forwardFixDate,
										Model->forwardFixTime,
										Model->valDate,
										Model->fixTime,
										Model->valTime,
										Model->fwd,
										Model->SigmaFwd,
										Model->AlphaFwd,
										Model->BetaFwd,
										Model->RhoFwd,
										Model->matrix,
										Model->nSimul,
										Model->tstarTime,
										Model->tstarDate,
										Model->dom_yc,
										Model->for_yc,
										Model->mergeTimes,
										Model->mergeNtimes,
										Model->sigDom,
										Model->lda_dom,
										Model->sigFor,
										Model->lda_for,
										Model->mergeTimes,
										Model->sigFx,
										Model->mergeNtimes,
										Model->mergeTimes,
										Model->corDF,
										Model->corDFx,
										Model->corFFx,
										Model->mergeNtimes,
										2,
										strikes,
										Params->solve_on_vol,
										Params->solve_on_vol,
										Model->SVfudgeSwitchLevel,
										Model->SVfudgeProba,
										Model->SVfudgeAlpha,
										vols);

		if(err)	goto FREE_RETURN;
	}

	if (Params->solve_on_vol)
	{
		Inst->vol_down = vols[0];
		Inst->vol_up = vols[1];
		*InstPrice = Inst->vol_up + Inst->vol_down - 2.0 * Inst->vol_atm;
	}
	else
	{
		Inst->price_down = vols[0];
		Inst->price_up = vols[1];
		*InstPrice = Inst->price_up + Inst->price_down - 2.0 * Inst->price_atm;
	}

FREE_RETURN:
	if(strikes) free(strikes);
	if(vols) free(vols);

	return err;
}

/* All the Functions for calibration of VOL on ATM VOL */
/* ********************************************************* */

Err	SABR_fwd_GetTarget_VOL(	void				*Inst_,
							void				*Params_,
							void				*Model_,
							CALIBGEN_PARAMS	CalibConsts,
							double				*target)
{
SABR_FWD_INST		Inst	= (SABR_FWD_INST) Inst_;
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;

	
	if (Params->solve_on_vol)
	{
		*target = Inst->target_vol_atm;
	}
	else
	{
		*target = Inst->target_price_atm;
	}

	return NULL;
}

Err	SABR_fwd_GetFirstGuess_VOL(	void				*Model_,
								void				*Params_,
								int					index_param,
								double				target,
								double				*param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;
double	logVol;

Err err = NULL;

	err = Fx3DtsImpliedVol_corr(	Model->valTime,
									Model->forwardFixTime,
									Model->fixTime,
									Model->mergeTimes,
									(long) Model->mergeNtimes,
									Model->sigDom,
									Model->lda_dom,
									Model->sigFor,
									Model->lda_for,
									Model->mergeTimes,
									Model->sigFx,
									(long) Model->mergeNtimes,
									Model->mergeTimes,
									Model->corDF,
									Model->corDFx,
									Model->corFFx,
									(long) Model->mergeNtimes,
									&logVol); 	
	
	err = srt_f_optsarbvol(	Model->fwd,
							Model->fwd,
							Model->valTime,
							logVol,
							Model->AlphaFwd,
							Model->BetaFwd,
							Model->RhoFwd,
							SRT_LOGNORMAL,
							SRT_BETAVOL,
							param1);

	return NULL;
}

Err	SABR_fwd_GetSecondGuess_VOL(void				*Model_,
								void				*Params_,
								int					index_param,
								double				param1,
								double				price1,
								double				target,
								double				*param2)
{
	if (price1 > target)
	{
		*param2 = param1 * 0.98;
	}
	else
	{
		*param2 = param1 * 1.02;
	}	

	return NULL;
}

Err	SABR_fwd_GetLimitAndLastParam_VOL(	void				*Model_,
										CALIBGEN_PARAMS	CalibConsts,
										void				*Params_,
										int					index_param,
										double				*last_param,
										double				*limit_down,
										double				*limit_up)
{
	*limit_down = 1e-5;
	*limit_up = 1e5;	

	return NULL;
}

Err	SABR_fwd_BumpParam_VOL(	void				*Model_,
							void				*Params_,
							int					index_param,
							double				param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	Model->SigmaFwd = param1;

	return NULL;
}

Err	SABR_fwd_SetParam_VOL(	void				*Model_,
							CALIBGEN_PARAMS	CalibConsts,
							void				*Params_,
							int					index_param,
							double				param1)
{
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

	Model->SigmaFwd = param1;

	return NULL;
}

Err	SABR_fwd_ExtrapolParam_VOL(void				*Model_,
								void				*Params_,
								int					index_param)
{
	return NULL;
}

Err	SABR_fwd_UpdateConstsAfterParam_VOL(void				*Inst_,
										void				*InstConst_,
										void				*Params_,
										void				*Model_,
										CALIBGEN_PARAMS	CalibConsts)
{
	return NULL;
}

Err	SABR_fwd_PriceInst_VOL(	void				*Inst_,
							void				*InstConst_,
							void				*Params_,
							void				*Model_,
							double				*InstPrice)
{
SABR_FWD_PARAMS	Params	= (SABR_FWD_PARAMS) Params_;
SABR_FWD_INST		Inst	= (SABR_FWD_INST) Inst_;
SABR_FWD_MODEL		Model	= (SABR_FWD_MODEL) Model_;

double	*strikes = NULL, *vols = NULL;

Err	err = NULL;

	strikes = (double*)calloc(1,sizeof(double));
	vols = (double*)calloc(1,sizeof(double));
	if(!strikes || !vols)
	{
		err = "SABR_fwd_PriceInst_BT: memory allocation failure";
		goto FREE_RETURN;
	}

	strikes[0] = Model->fwd;

	/* First Get the Prices */
	if(!Model->SVfudge)
	{
		err = fwdSABRpayoff(	Model->today,
								Model->forwardFixDate,
								Model->forwardFixTime,
								Model->valDate,
								Model->fixTime,
								Model->valTime,
								Model->fwd,
								Model->SigmaFwd,
								Model->AlphaFwd,
								Model->BetaFwd,
								Model->RhoFwd,
								Model->matrix,
								Model->nSimul,
								Model->tstarTime,
								Model->tstarDate,
								Model->dom_yc,
								Model->for_yc,
								Model->mergeTimes,
								Model->mergeNtimes,
								Model->sigDom,
								Model->lda_dom,
								Model->sigFor,
								Model->lda_for,
								Model->mergeTimes,
								Model->sigFx,
								Model->mergeNtimes,
								Model->mergeTimes,
								Model->corDF,
								Model->corDFx,
								Model->corFFx,
								Model->mergeNtimes,
								1,
								strikes,
								Params->solve_on_vol,
								Params->solve_on_vol,
								vols);

		if(err) goto FREE_RETURN;
	}
	else
	{
		err = fwdSABRpayoffAlphaFudge(	Model->today,
										Model->forwardFixDate,
										Model->forwardFixTime,
										Model->valDate,
										Model->fixTime,
										Model->valTime,
										Model->fwd,
										Model->SigmaFwd,
										Model->AlphaFwd,
										Model->BetaFwd,
										Model->RhoFwd,
										Model->matrix,
										Model->nSimul,
										Model->tstarTime,
										Model->tstarDate,
										Model->dom_yc,
										Model->for_yc,
										Model->mergeTimes,
										Model->mergeNtimes,
										Model->sigDom,
										Model->lda_dom,
										Model->sigFor,
										Model->lda_for,
										Model->mergeTimes,
										Model->sigFx,
										Model->mergeNtimes,
										Model->mergeTimes,
										Model->corDF,
										Model->corDFx,
										Model->corFFx,
										Model->mergeNtimes,
										1,
										strikes,
										Params->solve_on_vol,
										Params->solve_on_vol,
										Model->SVfudgeSwitchLevel,
										Model->SVfudgeProba,
										Model->SVfudgeAlpha,
										vols);

		if(err) goto FREE_RETURN;
	}

	if (Params->solve_on_vol)
	{
		Inst->vol_atm = vols[0];
		*InstPrice = Inst->vol_atm;
	}
	else
	{
		Inst->price_atm = vols[0];
		*InstPrice = Inst->price_atm;
	}

FREE_RETURN:
	if(strikes) free(strikes);
	if(vols) free(vols);

	return err;
}
/*=====================================
PRECALCULATIONS
=====================================*/
Err calib_sabr_fwd_rr_bt_given_beta_precalc(	long	today,
												long	spot_date,
												double	spot_fx,
												long	forwardFixDate,
												double	forwardFixTime,
												long	forwardValDate,
												double	forwardValTime,
												long	fixDate,
												double	fixTime,
												long	valDate,
												double	valTime,
												char	*dom_yc,
												char	*for_yc,
												double	*sigma_time_dom, 
												int		sigma_n_dom,
												double	*sigma_dom,
												double	lda_dom,
												double	*sigma_time_for, 
												int		sigma_n_for,
												double	*sigma_for, 
												double	lda_for,
												double	*sigma_time_fx, 
												double	*sigma_fx, 
												int		sigma_n_fx,
												double	*corr_times, 
												double	*correl_dom_for, 
												double	*correl_dom_fx, 
												double	*correl_for_fx,
												int		corr_n_times,
												double	SABRvol,
												double	SABRalpha,
												double	SABRbeta,
												double	SABRrho,
												double	SABRvol2,
												double	SABRalpha2,
												double	SABRbeta2,
												double	SABRrho2,
												double	FirstGuessSABRvol,
												double	FirstGuessSABRalpha,
												double	FirstGuessSABRbeta,
												double	FirstGuessSABRrho,
												int		nSimul,
												int		nPoints,
												int		nIter,
												int		nStd,
												double	std,
												int		smileModel,
												int		outputVol,
												int		SVfudge,
												double	SVfudgeSwitchLevel,
												double	SVfudgeProba,
												double	SVfudgeAlpha,
												double  *Results)
{
double		forward, discount;
double		*mergeTimes = NULL, 
			*sigDom = NULL, 
			*sigFor = NULL, 
			*sigFx = NULL, 
			*corDF = NULL, 
			*corDFx = NULL, 
			*corFFx = NULL;
int			mergeNtimes;

double		**matrix = NULL;

/* For the measure */
double		tstarTime; 
int			tstarDate;

/* For the calibration */
double		strike_down, strike_up,
			vol_down, vol_atm, vol_up,
			price_down, price_atm, price_up;
double		precision = 1e-4;
int			nb_iter_max = 10;
double		n_std = 1.0;


Err err = NULL;


	/* ==================================
	get the end time and date
	================================== */
	tstarDate = forwardFixDate;
	tstarTime = forwardFixTime;

	forward = spot_fx 
			* swp_f_df( spot_date, tstarDate, for_yc)
			/ swp_f_df( spot_date, tstarDate, dom_yc);

	matrix = dmatrix (0, nSimul-1, 0, 3-1);

	/* ==================================
	Merge the TS
	================================== */
	err = merge_rates_fx_corr_ts(	sigma_time_dom, sigma_dom, sigma_n_dom,
									sigma_time_for, sigma_for, sigma_n_for,
									sigma_time_fx, sigma_fx, sigma_n_fx,
									corr_times, correl_dom_for, correl_dom_fx, 
									correl_for_fx,corr_n_times,
									&mergeTimes, &sigDom, &sigFor, &sigFx, &corDF, &corDFx, &corFFx, &mergeNtimes);

	if(err) goto FREE_RETURN;

	err  = OTCgetCopula(	forward,
							forwardFixTime,
							tstarTime,
							mergeTimes, 
							mergeNtimes,
							sigDom,
							lda_dom,
							sigFor, 
							lda_for,
							sigFx, 
							corDF, 
							corDFx, 
							corFFx,
							SABRvol,
							SABRalpha,
							SABRbeta,
							SABRrho,
							nSimul,
							nPoints,
							nIter,
							nStd,
							std,
							smileModel,
							matrix);

	if (err)	goto FREE_RETURN;

	/* ==================================
	Payoff 
	================================== */
	forward = spot_fx 
				* swp_f_df( spot_date, valDate, for_yc)
				/ swp_f_df( spot_date, valDate, dom_yc);

	Results[0] = FirstGuessSABRvol;
	Results[1] = FirstGuessSABRalpha;
	Results[2] = FirstGuessSABRbeta;
	Results[3] = FirstGuessSABRrho;
	Results[4] = 0.0;

	strike_down = forward * exp(- n_std * SABRvol2 * sqrt(fixTime));
	strike_up = forward * exp(n_std * SABRvol2 * sqrt(fixTime));

	vol_atm = SABRvol2;

	err = srt_f_optsarbvol(	forward,
							strike_down,
							valTime,
							SABRvol2,
							SABRalpha2,
							SABRbeta2,
							SABRrho2,
							SRT_LOGNORMAL,
							SRT_LOGNORMAL,
							&vol_down);

	if(err) goto FREE_RETURN;

	err = srt_f_optsarbvol(	forward,
							strike_up,
							valTime,
							SABRvol2,
							SABRalpha2,
							SABRbeta2,
							SABRrho2,
							SRT_LOGNORMAL,
							SRT_LOGNORMAL,
							&vol_up);

	if(err) goto FREE_RETURN;

	discount = swp_f_df( today, valDate, dom_yc);

	price_atm = srt_f_optblksch(	forward,
									forward,
									vol_atm,
									valTime,
									discount,
									outputVol? SRT_CALL:SRT_PUT,
									PREMIUM);

	price_down = srt_f_optblksch(	forward,
									strike_down,
									vol_down,
									valTime,
									discount,
									outputVol? SRT_CALL:SRT_PUT,
									PREMIUM);

	price_up = srt_f_optblksch(	forward,
								strike_up,
								vol_up,
								valTime,
								discount,
								outputVol? SRT_CALL:SRT_PUT,
								PREMIUM);

	err =  calib_sabr_fwd_rr_bt_given_beta(	forward,
											today,
											forwardFixDate,
											valDate,
											tstarDate,
											forwardFixTime,
											fixTime,
											valTime,
											tstarTime,
											matrix,
											nSimul,
											outputVol,
											dom_yc,
											for_yc,
											mergeTimes,
											mergeNtimes,	
											sigDom,
											lda_dom,
											sigFor,
											lda_for,
											sigFx,
											corDF,
											corDFx,
											corFFx,
											vol_atm,
											price_atm,
											strike_down,
											vol_down,
											price_down,
											strike_up,
											vol_up,
											price_up,
											SVfudge,
											SVfudgeSwitchLevel,
											SVfudgeProba,
											SVfudgeAlpha,
											&(Results[0]),
											&(Results[1]),
											Results[2],
											&(Results[3]),
											precision,
											nb_iter_max,
											&(Results[4]));

	if(err) goto FREE_RETURN;
	
FREE_RETURN:

	return err;
}

/*=====================================
CALIBRATION
=====================================*/
Err calib_sabr_fwd_rr_bt_given_beta(	double				fwd,
										long				today,
										long				forwardFixDate,
										long				valDate,
										long				tstarDate,
										double				forwardFixTime,
										double				fixTime,
										double				valTime,
										double				tstarTime,
										double				**matrix,
										int					nSimul,
										int					outputVol,
										char				*dom_yc,
										char				*for_yc,
										double				*mergeTimes,
										int					mergeNtimes,	
										double				*sigDom,
										double				lda_dom,
										double				*sigFor,
										double				lda_for,
										double				*sigFx,
										double				*corDF,
										double				*corDFx,
										double				*corFFx,
										double				vol_atm,
										double				price_atm,
										double				strike_down,
										double				vol_down,
										double				price_down,
										double				strike_up,
										double				vol_up,
										double				price_up,
										int					SVfudge,
										double				SVfudgeSwitchLevel,
										double				SVfudgeProba,
										double				SVfudgeAlpha,
										double				*sigma,
										double				*alpha,
										double				fixed_beta,
										double				*rho,
										double				precision,
										int					nb_iter_max,
										double				*calib_err)
{
CALIBGEN_PARAMS	CalibParamsVol		= NULL;
CALIBGEN_PARAMS	CalibParamsRho		= NULL;
CALIBGEN_PARAMS	CalibParamsAlpha	= NULL;
CALIBFUNCTIONS		AllFunctionsVol		= NULL;
CALIBFUNCTIONS		AllFunctionsRho		= NULL;
CALIBFUNCTIONS		AllFunctionsAlpha	= NULL;
SABR_FWD_INST		Inst				= NULL;
SABR_FWD_MODEL		Model				= NULL;
SABR_FWD_PARAMS	Params				= NULL;

Err	err = NULL;

	CalibParamsVol = calloc(1, sizeof(CALIBGEN_Params));
	CalibParamsAlpha = calloc(1, sizeof(CALIBGEN_Params));
	CalibParamsRho = calloc(1, sizeof(CALIBGEN_Params));

	AllFunctionsVol = calloc(1, sizeof(CalibFunctions));
	AllFunctionsAlpha = calloc(1, sizeof(CalibFunctions));
	AllFunctionsRho = calloc(1, sizeof(CalibFunctions));

	Inst = calloc(1, sizeof(sabr_fwd_inst));
	Model = calloc(1, sizeof(sabr_fwd_model));
	Params = calloc(1, sizeof(sabr_fwd_params));

	if (!CalibParamsVol || !CalibParamsAlpha || !CalibParamsRho || 
		!AllFunctionsVol || !AllFunctionsAlpha || !AllFunctionsRho
		|| !Inst || !Model || !Params)
	{
		err = "Memory allocation faillure in calib_sabr_rr_bt_given_beta";
		goto FREE_RETURN;
	}
	
	/* Initialise Newton params */
	err = Initialise_CalibParams(	1,
									precision,
									nb_iter_max,
									1,
									1,
									1,
									0,
									0.0,
									CalibParamsVol);

	if (err) goto FREE_RETURN;

	err = Initialise_CalibParams(	1,
									precision,
									nb_iter_max,
									1,
									1,
									1,
									0,
									0.0,
									CalibParamsAlpha);

	if (err) goto FREE_RETURN;

	err = Initialise_CalibParams(	1,
									precision,
									nb_iter_max,
									1,
									1,
									1,
									0,
									0.0,
									CalibParamsRho);

	if (err) goto FREE_RETURN;

	/* Initialise Model */
	Model->fwd = fwd;

	Model->today = today;
	Model->forwardFixDate = forwardFixDate;
	Model->valDate = valDate;
	Model->tstarDate = tstarDate;

	Model->forwardFixTime = forwardFixTime;
	Model->fixTime = fixTime;
	Model->valTime = valTime;
	Model->tstarTime = tstarTime;
	
	Model->matrix = matrix;
	Model->nSimul = nSimul;

	Model->dom_yc = dom_yc;
	Model->for_yc = for_yc;

	Model->mergeTimes = mergeTimes;
	Model->mergeNtimes = mergeNtimes;	

	Model->sigDom = sigDom;
	Model->lda_dom = lda_dom;
	Model->sigFor = sigFor;
	Model->lda_for = lda_for;
	Model->sigFx = sigFx;
	Model->corDF = corDF;
	Model->corDFx = corDFx;
	Model->corFFx = corFFx;

	Model->fwd = fwd;
	Model->SigmaFwd = *sigma;

	if (fabs(*alpha) > 1.0E-04)
	{
		Model->AlphaFwd = *alpha;
	}
	else
	{
		Model->AlphaFwd = SABR_fwd_DEFAULT_ALPHA;
	}

	Model->BetaFwd = fixed_beta;
	Model->RhoFwd = *rho;
	Model->RhoAlphaFwd = Model->RhoFwd * Model->AlphaFwd;

	/* For the SV Fudge */
	Model->SVfudge = SVfudge;
	Model->SVfudgeSwitchLevel = SVfudgeSwitchLevel;
	Model->SVfudgeProba = SVfudgeProba;
	Model->SVfudgeAlpha = SVfudgeAlpha;


	/* Initialise Instrument */
	Inst->target_vol_atm = vol_atm;
	Inst->target_price_atm = price_atm;
	Inst->strike_down = strike_down;
	Inst->target_vol_down = vol_down;
	Inst->target_price_down = price_down;
	Inst->strike_up = strike_up;
	Inst->target_vol_up = vol_up;
	Inst->target_price_up = price_up;
	

	/* Initialise Functions */
	AllFunctionsVol->GetTarget = SABR_fwd_GetTarget_VOL;
	AllFunctionsVol->BumpParam = SABR_fwd_BumpParam_VOL;
	AllFunctionsVol->ExtrapolParam = SABR_fwd_ExtrapolParam_VOL;
	AllFunctionsVol->GetFirstGuess = SABR_fwd_GetFirstGuess_VOL;
	AllFunctionsVol->GetLimitAndLastParam = SABR_fwd_GetLimitAndLastParam_VOL;
	AllFunctionsVol->GetSecondGuess = SABR_fwd_GetSecondGuess_VOL;
	AllFunctionsVol->PriceInst = SABR_fwd_PriceInst_VOL;
	AllFunctionsVol->SetParam = SABR_fwd_SetParam_VOL;
	AllFunctionsVol->UpdateConstsAfterParam = SABR_fwd_UpdateConstsAfterParam_VOL;

	AllFunctionsAlpha->GetTarget = SABR_fwd_GetTarget_BT;
	AllFunctionsAlpha->BumpParam = SABR_fwd_BumpParam_BT;
	AllFunctionsAlpha->ExtrapolParam = SABR_fwd_ExtrapolParam_BT;
	AllFunctionsAlpha->GetFirstGuess = SABR_fwd_GetFirstGuess_BT;
	AllFunctionsAlpha->GetLimitAndLastParam = SABR_fwd_GetLimitAndLastParam_BT;
	AllFunctionsAlpha->GetSecondGuess = SABR_fwd_GetSecondGuess_BT;
	AllFunctionsAlpha->PriceInst = SABR_fwd_PriceInst_BT;
	AllFunctionsAlpha->SetParam = SABR_fwd_SetParam_BT;
	AllFunctionsAlpha->UpdateConstsAfterParam = SABR_fwd_UpdateConstsAfterParam_BT;	

	AllFunctionsRho->GetTarget = SABR_fwd_GetTarget_RR;
	AllFunctionsRho->BumpParam = SABR_fwd_BumpParam_RR;
	AllFunctionsRho->ExtrapolParam = SABR_fwd_ExtrapolParam_RR;
	AllFunctionsRho->GetFirstGuess = SABR_fwd_GetFirstGuess_RR;
	AllFunctionsRho->GetLimitAndLastParam = SABR_fwd_GetLimitAndLastParam_RR;
	AllFunctionsRho->GetSecondGuess = SABR_fwd_GetSecondGuess_RR;
	AllFunctionsRho->PriceInst = SABR_fwd_PriceInst_RR;
	AllFunctionsRho->SetParam = SABR_fwd_SetParam_RR;
	AllFunctionsRho->UpdateConstsAfterParam = SABR_fwd_UpdateConstsAfterParam_RR;

	Params->solve_on_vol = outputVol;
	Params->AllFunctions_ForVol = AllFunctionsVol;
	Params->CalibParams_ForVol = CalibParamsVol;
	Params->AllFunctions_ForRho = AllFunctionsRho;
	Params->CalibParams_ForRho = CalibParamsRho;
	Params->AllFunctions_ForAlpha = AllFunctionsAlpha;
	Params->CalibParams_ForAlpha = CalibParamsAlpha;

	/* Solve !!! */
	err = CalibrateParamTS(	0,
							0,
							&Inst,
							&Params,
							Params,
							Model,
							CalibParamsRho,
							AllFunctionsRho);

	if (err) goto FREE_RETURN;

	/* Outputs */
	*sigma = Model->SigmaFwd;
	*alpha = Model->AlphaFwd;
	*rho = Model->RhoFwd;

	if(outputVol)
	{
		*calib_err = pow(Inst->vol_down - Inst->target_vol_down, 2.0)
					+ pow(Inst->vol_up - Inst->target_vol_up, 2.0)
					+ pow(Inst->vol_atm - Inst->target_vol_atm, 2.0);
	}
	else
	{
		*calib_err = pow(Inst->price_down - Inst->target_price_down, 2.0)
					+ pow(Inst->price_up - Inst->target_price_up, 2.0)
					+ pow(Inst->price_atm - Inst->target_price_atm, 2.0);
	}
	*calib_err = sqrt(*calib_err);

FREE_RETURN:

	if (CalibParamsRho)
	{
		Free_CalibParams(CalibParamsRho);
		free(CalibParamsRho);
	}

	if (CalibParamsAlpha)
	{
		Free_CalibParams(CalibParamsAlpha);
		free(CalibParamsAlpha);
	}

	if (CalibParamsVol)
	{
		Free_CalibParams(CalibParamsVol);
		free(CalibParamsVol);
	}

	if (AllFunctionsRho) free(AllFunctionsRho);
	if (AllFunctionsAlpha) free(AllFunctionsAlpha);
	if (AllFunctionsVol) free(AllFunctionsVol);
	if (Inst) free (Inst);

	if (Model)
	{
		Free_SABRfwd_Model(Model);
		free (Model);
	}
	if (Params)
	{
		free (Params);
	}

	return err;
}









