
#ifndef __DIAGCALIBDLM_H
#define __DIAGCALIBDLM_H

#define MAX_CPN			600
#define MAX_INST		600


void constant_or_linear_interpolation_dlm(	double	*dates,
											double	*values,
											int		n_dates,
											double	date_wanted,
											double	*value_wanted);

double solve_for_next_coef_dlm(	double	**res_iter,
								int		nb_iter,
								double	premium_tgt,
								int		method);			/* 0: linear 1: quadratic */

Err cpd_calib_diagonal_dlm(
	/*	Market */
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */	
	Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),

	/* Get vol ref */
	char			*vol_ref_rate_name,
	
	/* Long Instruments */
	char			*instr_long_freq,					/*	Frequency and basis of instruments */
	char			*instr_long_basis,
	char			*long_ref_rate_name,				/*	Name of the reference rate */

	int				num_ex_datesl,						/*	Long Exercise dates */
	long			*ex_datel_,							/*	Supposed to be sorted */
	int				*cal_datel,							/*	1: use ex_date as calibration date, 0: don't */
	char			**end_tenorl,						/*	Tenors of the underlying instruments */
	long			end_datel,							/*	End date for diagonal */
	double			*strikel_,							/*	Strikes */	

	CPD_DIAG_CALIB_PARAM paraml,

	/* Short Instruments */
	char			*instr_short_freq,					/*	Frequency and basis of instruments */
	char			*instr_short_basis,
	char			*short_ref_rate_name,				/*	Name of the reference rate */

	int				num_ex_datess,						/*	Short Exercise dates */
	long			*ex_dates_,							/*	Supposed to be sorted */
	int				*cal_dates,							/*	1: use ex_date as calibration date, 0: don't */
	char			**end_tenors,						/*	Tenors of the underlying instruments */
	long			end_dates,							/*	End date for diagonal */
	double			*strikes_,							/*	Strikes */	
	double			*weights_,

	CPD_DIAG_CALIB_PARAM params,
						
	/*	Model */
	int				fix_lambda,
	int				one_factor_equi,				/*	Calibrate 2 Factor to 1 Factor price */
	int				nlam,							/*	Lambda TS: may NOT be changed in the process */
	double			lam_time[],
	double			lam[],
	double			lam_shift[],
	int				one2F,							/*	Number of factors */
	double			alpha,							/*	Alpha, Gamma, Rho (2F only) */
	double			gamma,
	double			rho,	

	/*	Shift Parameters */
	int				num_long_vol_shift,
	double			*long_vol_shift_time,
	double			*long_vol_shift,

	int				num_short_vol_shift,
	double			*short_vol_shift_time,
	double			*short_vol_shift,

	/*	Output */
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	/*	Parameters */
	
	DIAG_CALIB_LM_PARAMS lm_params,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data);					/*	NULL = don't save calibration instrument data */

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err lgmcalibzetalambda_tauts_dlm(
int				ncpnl,								/*	Total number of cash-flow dates */
double			cpn_timel[],							/*	Cash-Flow times */
double			cpn_dfl[],							/*	Df to cash-flow dates */
double			cpn_cvgl[],							/*	cvg from i-1 to i */
int				nexl,								/*	Total number of exercise dates */
double			ex_timel[],							/*	Exercise times */
int				ex_cpnl[],							/*	Index of the first cash-flow to be exercised */
int				ex_ncpnl[],							/*	Number of coupons in each caplet */
double			ex_strikel[],						/*	Strikes for diagonal */
double			ex_pricel[],						/*	Market prices for diagonal */
double			ex_vegal[],							/*	Market vegas for diagonal */

int				ncpns,								/*	Total number of cash-flow dates */
double			cpn_times[],							/*	Cash-Flow times */
double			cpn_dfs[],							/*	Df to cash-flow dates */
double			cpn_cvgs[],							/*	cvg from i-1 to i */
int				nexs,								/*	Total number of exercise dates */
double			ex_times[],							/*	Exercise times */
int				ex_cpns[],							/*	Index of the first cash-flow to be exercised */
int				ex_ncpns[],							/*	Number of coupons in each caplet */
double			ex_strikes[],						/*	Strikes for diagonal */
double			ex_prices[],						/*	Market prices for short */
double			ex_weights[],						/*	Weights for short */
double			ex_vegas[],							/*	Market vega for short */

double			ex_zeta[],							/*	Output: zetas */
int				fix_lambda,							/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
int				nlam,								/*	Lambda TS: may NOT be changed in the process */
double			lam_time[],
double			lam[],
int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
CPD_DIAG_CALIB_PARAM paraml,
CPD_DIAG_CALIB_PARAM params,
DIAG_CALIB_LM_PARAMS lm_params,
double			*lam_sens);

typedef struct
{
	long	lToday;

	long	iNCpn;
	long	lCpnDate[MAX_CPN];
	long	lCpnTheoDate[MAX_CPN];
	double	dCpnTime[MAX_CPN];	
	double	dCpnCvg[MAX_CPN];
	double	dCpnDf[MAX_CPN];	

} CalibCpnScheduleDLM, *CALIBCPNSCHEDULEDLM;

typedef struct
{
	int		iNExe;
	long	*lExeDates;
	double	*dExeTimes;
	long	*lStartDate;
	long	*lTheoEndDates;
	long	*lActEndDates;
	int		iStartCpn[MAX_INST];
	int		iEndCpn[MAX_INST];	
	double	*dStrikes;
	double	*dStrikesS1;
	double	*dStrikesS2;
	double	*dWeights;
	double	dPrices[MAX_INST];
	double	dATMStrike[MAX_INST];
	double	dATMPrice[MAX_INST];
	double	dATMVega[MAX_INST];

} CalibExeScheduleDLM, *CALIBEXESCHEDULEDLM;

void	AllocateCalibExeSchedule(	long				*lExeDates,
									double				*dExeTimes,
									long				*lStartDates,
									long				*lTheoEndDates,
									long				*lActEndDates,
									double				*dStrikes,
									double				*dStrikesS1,
									double				*dStrikesS2,
									double				*dWeights,
									CALIBEXESCHEDULEDLM CalibExeSchedule);

struct CalibInstrumentDLM;
typedef struct CalibInstrumentDLM CalibInstrumentDLM;
typedef CalibInstrumentDLM *CALIBINSTRUMENTDLM;

struct CalibInstrumentDLM
{			
	int					iNbCoupon;
	int					iStartCpn;
	int					iEndCpn;

	double				dCpnValues[MAX_CPN];

	char				*cFreq;
	char				*cBasis;
	char				*cRefRate;

	SrtCompounding		iFreq;
	SrtBasisCode		iBasis;

	long				*lCouponDate;
	double				*dCouponTime;
	double				*dCouponCvg;
	double				*dCouponDf;

	int					iIsCMS;
	long				lPayDate;
	double				dPayTime;

	int					iIsLiborOption;
	long				lLiborFixDate;
	double				dLiborFixTime;
	long				lLiborStartDate;
	double				dLiborStartTime;
	long				lLiborEndDate;
	double				dLiborEndTime;
	double				dLiborCvg;
	double				dLiborFwd;
	double				dLiborSpread;
	double				dLiborMargin;

	long				lExeDate;
	double				dExeTime;
	long				lStartDate;	
	long				lTheoEndDate;
	long				lEndDate;

	double				dLevel;
	double				dSumDf;
	double				dFwdCash;
	double				dSpread;
	double				dATMVol;
						
	int					iNbStrike;
	double				*dStrike;	
	double				*dPrice;
	double				*dNormalVol;
	double				*dVega;
	double				*dNbStd;
	double				dWeight;

	double				*dModelPrice;	

	CalibInstrumentDLM	*NextInst;

};

Err	AllocateCalibInst(	int					iNbStrike,
						CALIBINSTRUMENTDLM	CalibInst);

void FreeCalibInst(	CALIBINSTRUMENTDLM	CalibInst);

typedef struct
{
	CALIBCPNSCHEDULEDLM	sCalibCpnSchedule;
	CALIBEXESCHEDULEDLM	sCalibExeSchedule;

	int					iNbInst;
	CALIBINSTRUMENTDLM	sCalibInst;

} AllCalibInstrumentsDLM, *ALLCALININSTRUMENTSDLM;

Err	AllocateAllCalibInst(int					iNbInst,
						 int					iNbStrike,
						 ALLCALININSTRUMENTSDLM	AllCalibInst);

void FreeAllCalibInst(ALLCALININSTRUMENTSDLM	AllCalibInst);

Err	Construct_CalibSchedule(
							char				*yc_name,
							long				today,
							char				*instr_freq,
							char				*instr_basis,
							int					num_ex_dates,					/*	Exercise dates */
							long				*ex_date,						/*	Supposed to be sorted */							
							int					*cal_date,						/*	1: use ex_date as calibration date, 0: don't */
							char				**end_tenor,					/*	Tenors of the underlying instruments */
							long				end_date,
							double				*strikes,
							double				*strikesS1,
							double				*strikesS2,
							double				*weights,

							CALIBCPNSCHEDULEDLM	CalibCpnSchedule,
							CALIBEXESCHEDULEDLM	CalibExeSchedule);

Err	Reduct_ExeSchedule(	
						CALIBCPNSCHEDULEDLM		CalibCpnSchedule,
						CALIBEXESCHEDULEDLM		CalibExeSchedule,
						int						*cal_date,
						CPD_DIAG_CALIB_PARAM	param);

Err	Calculate_CalibInst(long					today,
						char					*yc_name,						/*	Name of the yield curve */
						char					*vol_curve_name,				/*	Name of the market vol curve */	
						Err						(*get_cash_vol)(				/*	Function to get cash vol from the market */
															char	*vol_curve_name,	
															double	start_date, 
															double	end_date,
															double	cash_strike,
															int		zero,
															char	*ref_rate_name,
															double	*vol,
															double	*power),
						char					*vol_ref_rate_name,
						char					*instr_freq,
						char					*instr_basis,
						char					*ref_rate_name,

						int						is_for_calibration,
						
						CPD_DIAG_CALIB_PARAM	param,

						CALIBCPNSCHEDULEDLM		CalibCpnSchedule,
						CALIBEXESCHEDULEDLM		CalibExeSchedule,

						ALLCALININSTRUMENTSDLM	AllCalibInst,

						/* information in the case of the short instruments */
						ALLCALININSTRUMENTSDLM	AllCalibShortInst);

typedef Err (*GET_CASH_VOL_FUNC)(
				char	*vol_curve_name,	
				double	start_date, 
				double	end_date,
				double	cash_strike,
				int		zero,
				char	*ref_rate_name,
				double	*vol,
				double	*power);

typedef Err (*GET_CASH_VOL_CONVERT_FUNC)(
				GET_CASH_VOL_FUNC	get_cash_vol,
				char				*vol_curve_name,	
				double				start_date, 
				double				end_date,
				double				mat,
				double				cash_fwd,
				double				cash_strike,
				char				*ref_rate_name,
				SrtDiffusionType	vol_type,
				double				*vol);

/*	Function to get and convert the volatility */
Err	GetCashVolAndConvert(
				GET_CASH_VOL_FUNC	get_cash_vol,
				char				*vol_curve_name,	
				double				start_date, 
				double				end_date,
				double				mat,
				double				cash_fwd,
				double				cash_strike,
				char				*ref_rate_name,
				SrtDiffusionType	vol_type,
				double				*vol);

#endif