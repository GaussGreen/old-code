#ifndef SABRFWDCALIB
#define	SABRFWDCALIB

#include "srt_h_all.h"
#include "DiagcalibGen.h"
#include "OTCutils.h"
#include "Fx3FCalib.h"

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
												double  *Results);

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
										double				*calib_err);

/* All the structures needed */
typedef struct
{
	double				target_vol_atm;
	double				target_price_atm;
	double				vol_atm;
	double				price_atm;

	double				strike_up;
	double				target_vol_up;
	double				target_price_up;
	double				vol_up;
	double				price_up;
						
	double				strike_down;
	double				target_vol_down;
	double				target_price_down;
	double				vol_down;
	double				price_down;

	SrtDiffusionType	vol_type;

}	sabr_fwd_inst, *SABR_FWD_INST;

typedef struct
{
	double				fwd;

	long				today;
	long				forwardFixDate;
	long				valDate;
	long				tstarDate;

	double				forwardFixTime;
	double				fixTime;
	double				valTime;
	double				tstarTime;
	
	double				SigmaFwd;
	double				AlphaFwd;
	double				BetaFwd;
	double				RhoFwd;
	double				RhoAlphaFwd;
	
	double				**matrix;
	int					nSimul;

	char				*dom_yc;
	char				*for_yc;

	double				*mergeTimes;
	int					mergeNtimes;	

	double				*sigDom;
	double				lda_dom;
	double				*sigFor;
	double				lda_for;
	double				*sigFx;
	double				*corDF;
	double				*corDFx;
	double				*corFFx;

	/* For the Stoch Vol Fudge */
	int					SVfudge;
	double				SVfudgeSwitchLevel;
	double				SVfudgeProba;
	double				SVfudgeAlpha;
}	sabr_fwd_model, *SABR_FWD_MODEL;

typedef struct
{
	int					solve_on_vol;
	int					solve_on_rhoalpha;
	CALIBFUNCTIONS		AllFunctions_ForVol;	
	CALIBGEN_PARAMS		CalibParams_ForVol;
	CALIBFUNCTIONS		AllFunctions_ForAlpha;	
	CALIBGEN_PARAMS		CalibParams_ForAlpha;
	CALIBFUNCTIONS		AllFunctions_ForRho;
	CALIBGEN_PARAMS		CalibParams_ForRho;

}	sabr_fwd_params, *SABR_FWD_PARAMS;

Free_SABRfwd_Model(SABR_FWD_MODEL		Model);


#endif