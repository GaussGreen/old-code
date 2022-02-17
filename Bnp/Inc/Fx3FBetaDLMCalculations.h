 
#ifndef FX3FBETADLMCALCULATIONS_H
#define FX3FBETADLMCALCULATIONS_H
#define MAX_HERMITE	100
#define	MAX_FWD_DLM	50
#define	MAX_X_DLM	50

#include "srt_h_all.h"
#include "Fx3FBetaDLMUtil.h"

Err	FxBetaDLM_FxOptionTS_struct(	long						opt_settlmt_date,
									long						fwd_start_date,
									long						fixing_date,
									int							nb_strike,
									double						*strikes,
									SrtCallPutType				callput,
									FxBetaDLM_model				*model,
									FxBetaDLM_OptNumerParams	*NumParams,
									FxBetaDLM_Hermite			*hermite,
									double						*premium);

Err	FxBetaDLM_FxOptionUnd(	char						*und,
							long						opt_settlmt_date,
							long						fwd_start_date,
							long						fixing_date,
							int							nb_strike,
							double						*strikes,
							SrtCallPutType				callput,
							FxBetaDLM_OptNumerParams	*NumParams,
							double						*premium);

typedef struct
{
	int				iInstIndex;

	long			lFwdStartDate;
	long			lExeDate;	
	long			lSettlementDate;
	long			lPayDate;
					
	double			dFwdStartTime;
	double			dExeTime;
	double			dSettlementTime;
	double			dPayTime;
	double			dVolTime;
	double			dSqVolTime;
					
	double			dFwdFx;
	double			dLnFwdFx;
	double			dDfDom;
	double			dDfFor;
	double			dDfPayDom;

	SrtCallPutType	sCallPutType;
	int				iNbStrike;
	double			*dStrike;
	double			*dLogStrike;
	double			*dPrice;
	double			*dLogVol;
	double			*dVega;

} FxBetaDLM_FxOptInst, *FXBETADLM_FXOPTINST;

Err	FxBetaDLM_Allocate_FxOptInst(	int					iIndex,
									int					iNbStrike,
									FxBetaDLM_FxOptInst	*Inst);

void FxBetaDLM_Free_FxOptInst(FxBetaDLM_FxOptInst	*Inst);

Err	FxBetaDLM_Setup_FxOptInst(	long				opt_settlmt_date,
								long				opt_pay_date,
								long				fwd_start_date,
								long				fixing_date,
								int					nb_strike,
								double				*strikes,
								double				*bs_vols,
								SrtCallPutType		callput,
								FxBetaDLM_model		*model,
								FxBetaDLM_FxOptInst	*Inst);

								

typedef struct
{
	double	dBetaDom;
	double	dBetaFor;

	double	dVarRatesDom;
	double	dVarRatesFor;
	double	dCovarRates;
	double	dStdRates;

	double	dTotalVarFFx;
	double	dVarFFx;
	double	dStdFFX;
	double	dCovarFFxDom;
	double	dCovarFFxFor;

	double	dVarRatesFor_down;
	double	dCovarRates_down;
	double	dTotalVarFFx_down;
	double	dVarFFx_down;
	double	dStdFFX_down;
	double	dCovarFFxDom_down;
	double	dCovarFFxFor_down;

	double	dVarRatesFor_mid;
	double	dCovarRates_mid;
	double	dTotalVarFFx_mid;
	double	dVarFFx_mid;
	double	dStdFFX_mid;
	double	dCovarFFxDom_mid;
	double	dCovarFFxFor_mid;

	double	dConstCoef;
	double	dLinCoef;
	double	dQuadCoef;

	double	dQuadConst;
	double	dQuadLin;
	double	dQuadQuad;
	double	dQuadMean;
	double	dQuadVar;
	double	dQuadStd;
	double	dQuadStdY;

	double	dQuadConst_down;
	double	dQuadLin_down;
	double	dQuadQuad_down;
	double	dQuadMean_down;
	double	dQuadVar_down;
	double	dQuadStd_down;
	double	dQuadStdY_down;

	double	dFwdAlpha,
			dFwdBeta,
			dFwdQuad2;

	double	dFwdAlpha_down,
			dFwdBeta_down,
			dFwdQuad2_down;


	double	coupon_hermite[MAX_HERMITE];

	/* Constant for CPD */
	int		iIsCPD;
	double	coef_b_const;
	double	coef_b_lin;

	double	fixed_adj_const;
	double	fixed_adj_lin;
	double	fixed_adj_lin_const;

	double	coef_b_const_down;
	double	coef_b_lin_down;

	double	fixed_adj_const_down;
	double	fixed_adj_lin_down;
	double	fixed_adj_lin_const_down;

	/* Smile Precalculation !!! */
	double	dPrecalcFwdFloor[MAX_FWD_DLM];
	double	dPrecalcFwdCap[MAX_FWD_DLM];
	double	dPrecalcX[MAX_X_DLM];
	double	dPrecalcFloorPriceStd[MAX_X_DLM][MAX_FWD_DLM];
	double	dPrecalcCapPriceStd[MAX_X_DLM][MAX_FWD_DLM];
	double	dPrecalcFloorPriceStdAlpha[MAX_X_DLM][MAX_FWD_DLM];
	double	dPrecalcCapPriceStdAlpha[MAX_X_DLM][MAX_FWD_DLM];

} FxBetaDLM_InstPrecalc, FXBETADLM_INSTPRECALC;

Err FxBetaDLM_Price_FxOptInst(	FxBetaDLM_FxOptInst			*Inst,
								FxBetaDLM_model				*Model,
								FxBetaDLM_InstPrecalc		*Precalc,
								FxBetaDLM_Hermite			*Hermite,
								FxBetaDLM_OptNumerParams	*cNumParams,
								double						*Premium);

Fx3DBetaDLM_ExpectAndVarGrfn(	long			nstp,
								double			*time,
								FxBetaDLM_model	*model,

								double			*dom_fwd,
								double			*dom_var_std,
								double			*dom_phi,
								double			*dom_beta,
								double			*for_fwd_const,
								double			*for_fwd_lin,
								double			*for_var_std,
								double			*for_phi,
								double			*for_beta,
								double			*fx_fwd,
								double			*fx_var_std,
								double			*ffx_var,
												
								int				mc_tree, /*0: MC, 1: Tree */
								double			*dom_for_cov_corr,
								double			*dom_fx_cov_corr,
								double			*for_fx_cov_corr,
								
								double			min_time);

Fx3DBetaDLM_PrecalcGRFNTreeQBeta(	long			nstp,
									double			*time,
									FxBetaDLM_model	*model,

									double			*dom_fwd,
									double			*dom_var,
									double			*dom_phi,
									double			*dom_beta,
									double			*for_fwd,
									double			*for_var,
									double			*for_phi,
									double			*for_beta,
									double			*fx_fwd,
									double			*fx_var,
									double			*ffx_var,
																					
									double			*dom_for_corr,
									double			*dom_fx_corr,
									double			*for_fx_corr,
									
									double			min_time);

Fx3DBetaDLM_DriftVolAndCorr(long			nstp,
							double			*time,
							FxBetaDLM_model	*model,

							double			*dom_phi,
							double			*dom_beta,
							double			*for_phi,
							double			*for_beta,
							double			*var_ffx,
									
							double			*dr_const_dom,
							double			*dr_coef_dom,
							double			*dr_const_for,
							double			*dr_coef1_for,
							double			*dr_coef2_for,
							double			*dr_coef3_for,
							double			*dr_const_fx,
							double			*dr_coef1_fx,
							double			*dr_coef2_fx,
							double			*dr_coef3_fx,
											
							int				*vol_change,
							double			*sig_dom,
							double			*sig_for,
							double			*sig_fx,
							double			*corr_dom_for,
							double			*corr_dom_fx,
							double			*corr_for_fx);

Err quadratic_bs_CPD(	FxBetaDLM_InstPrecalc	*Precalc,
						int						i,
						double					strike,
						double					LogStrike,
						SrtCallPutType			callput,
						double					*Premium);

Err FxBetaDLM_Price_ForwardFx(	FxBetaDLM_InstPrecalc	*Precalc,
								double					*Premium);

Err FxBetaDLM_GetEqui3FactorTS(	FxBetaDLM_model				*Model,
								FxBetaDLM_OptNumerParams	*NumParams,
								FxBetaDLM_Hermite			*hermite,
								int							nb_fx,
								double						*fx_maturities,
								double						**fx_cal_vols);

#endif