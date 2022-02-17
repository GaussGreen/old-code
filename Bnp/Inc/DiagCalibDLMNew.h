
#ifndef __DIAGCALIBDLMNEW_H
#define __DIAGCALIBDLMNEW_H

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to diagonal
		with lambda calibration */
Err cpd_calib_diagonal_dlm_new(
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

	/* Get Cash Vol ref */
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
	double			*weights_,							/*	Weights on secondary instruments */

	CPD_DIAG_CALIB_PARAM params,
							
	/*	Model */
	int				fix_lambda,
	int				one_factor_equi,				/*	Calibrate 2 Factor to 1 Factor price */
	int				nlam,							/*	Lambda TS: may NOT be changed in the process */
	double			lam_time[],
	double			lam[],
	double			lam_shift[],
	int				nfactor,						/*	Number of factors */
	double			alpha,							/*	Alpha, Gamma, Rho (2F only) */
	double			gamma,
	double			rho,		
		
	/*	Output */
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	
	/*	Parameters */	
	DIAG_CALIB_LM_PARAMS lm_params,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data);					/*	NULL = don't save calibration instrument data */

typedef struct
{
	/* Today Date */
	long	lToday;

	/* Number of Factors */
	int		iNbFactor;
	
	/* Merged Term Structures */
	int		iNbPWTime;
	double	*dPWTime;		
	double	*dSigma;					
	double	*dLambda;

	/* 2 Factor Parameters */		
	double	dAlpha;
	double	dGamma;
	double	dRho;
	
	/* Initial Lambda TS */
	int		iNbLam;
	double	*dLambdaTime;
	double	*dInitLambda;
	
	double	*dNewLambda;
	double	*dNewLambda2;

} LGM_model, *LGM_MODEL;

typedef struct
{
	double				*dCpn_G1;
	double				*dCpn_G2;

	double				*dEx_G1;
	double				*dEx_G2;

	double				*dExpFact1;
	double				*dExpFact2;
	double				*dExpFact12;	

	double				*dZeta1;
	double				*dZeta2;
	double				*dZeta12;

} LGM_ModelFactors, *LGM_MODELFACTORS;

typedef struct
{
	int					iNbLongInst;
	void				**AllLongInst;
	long				*lSigLongIndex;
	CALIBCPNSCHEDULEDLM	CalibCpnLongSchedule;
	CALIBEXESCHEDULEDLM	CalibExeLongSchedule;

	int					iNbShortInst;
	void				**AllShortInst;
	long				*lSigShortIndex;
	CALIBCPNSCHEDULEDLM	CalibCpnShortSchedule;
	CALIBEXESCHEDULEDLM	CalibExeShortSchedule;
	
	CALIBGEN_PARAMS		CalibParams;
	CALIBGEN_PARAMS		CalibParamsLambda;	
	
	CalibFunctions		CalibFunctionsForVol;
	CalibFunctions		CalibFunctionsForLambda;	
			
	CPD_DIAG_CALIB_PARAM	long_param;
	CPD_DIAG_CALIB_PARAM	short_param;

	double				sens_lambda;
	double				cum_var_lgm;

	/* All Factors for Long and Short */
	LGM_MODELFACTORS	sLongFactors;
	LGM_MODELFACTORS	sShortFactors;	

	int					index_vol;

} LGM_AllParams, *LGM_ALLPARAMS;

typedef struct
{
	long	iNCpn;
	double	*dCpn;		
		

} LGM_InstNumerParams, LGM_INSTNUMERPARAMS;

Err	GetTargetVol_LGM(	void				*Inst,
						void				*GlobalConst,
						void				*Model,
												
						CALIBGEN_PARAMS	CalibConsts,
						double				*target);

Err	GetFirstGuessVol_LGM(	void	*Model,
							void	*GlobalParam,
							int		vol_index,
							double	target,
							double	*vol1);

Err	GetSecondGuessVol_LGM(	void	*Model,
								void	*GlobalParam,
								int		vol_index,
								double	vol1,
								double	price1,
								double	target,
								double	*vol2);

Err	BumpVol_LGM(	void			*Model,
					void			*GlobalParam,
					int				vol_index,
					double			vol);

Err	SetVol_LGM(	void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					vol_index,
					double				vol);

Err	GetLimitAndLastVol_LGM(	void				*Model,
								CALIBGEN_PARAMS	CalibParams,
								void				*GlobalParam,

								int					vol_index,
								double				*last_vol,
								double				*limit_down,
								double				*limit_up);

Err	ExtrapolVol_LGM(	void				*Model,						
						void				*GlobalParam,
						int					last_vol_index);

Err	UpdateParamsAfterVol_LGM (void				*Inst_,
								void				*InstParam,
								void				*GlobalParam,
								void				*Model,
								CALIBGEN_PARAMS	CalibConsts);

Err	PriceInstVol_LGM(void				*Inst_,
						void			*InstParam,
						void			*GlobalParam,
						void			*Model,
						double			*InstPrice);

Err	GetTargetLambda_LGM(		void				*Inst_,
							void				*GlobalConst,
							void				*Model,
												
							CALIBGEN_PARAMS	CalibConsts,
							double				*target);

Err	GetFirstGuessLambda_LGM(	void			*Model,
							void			*GlobalConst,
							int				index_lambda,
							double			target,
							double			*param1);

Err	GetSecondGuessLambda_LGM(void			*Model,
							void			*GlobalConst,
							int				index_lambda,
							double			lambda1,
							double			price1,							
							double			target,
							double			*lambda2);

Err	GetLimitAndLastLambda_LGM(	void				*Model,
									CALIBGEN_PARAMS	CalibConsts,
									void				*GlobalConst,
									int					index_lambda,
									double				*last_lambda,
									double				*limit_down,
									double				*limit_up);

Err	BumpLambda_LGM(	void			*Model,
					void			*GlobalParam,
					int				index_lambda,
					double			lambda);

Err	SetLambda_LGM(	void				*Model,
					CALIBGEN_PARAMS	CalibConsts,
					void				*GlobalParam,
					int					index_lambda,
					double				lambda);

Err	ExtrapolLambda_LGM(	void				*Model,						
						void				*GlobalParam,
						int					last_smile_index);

Err	UpdateParamsAfterLambda_LGM (	void				*Inst_,
									void					*InstParam,
									void					*GlobalParam,
									void					*Model,
									CALIBGEN_PARAMS		CalibConsts);

Err	PriceInstLambda_LGM(void				*Inst_,
						void			*InstParam,
						void			*GlobalParam,
						void			*Model,
						double			*InstPrice);
#endif