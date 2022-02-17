#ifndef _SRTGRFNLGM2FAUTOCAL_H_
#define _SRTGRFNLGM2FAUTOCAL_H_


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Struct s_grfn_pricing_params to hold Grfn pricing Params for a variety of algorithms  //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
typedef struct 
{
	long  num_steps_time;
	long  num_steps_space;
	long  num_steps_time_mc;
	long  num_paths;
	int   do_pecs;
	int   do_jump;
	int	  use_stabler_pde_with_tau_term_structure;

} s_grfn_pricing_params;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Sort Static level function to calibrate an LGM 2 factor model and then price a Grfn tableau using	//
//	either PDE, MC or MC with optimal ex boundary.		PMc		29th Sept 03							//
//  Note the memory management policy in this routine. In particular the fact that we free the Product  //
//  Values when there has been an error but leave the calling routine to free in the absence of errors. //
//																										//
//	Modified to return primary and secondary model prices.		PMc		1stDec03						//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
Err  GrfnLGM2FAutocalCaller( char			*YieldCurveName,
							 char			*VolCubeName,
							 Err		   (*GetCashVol)(char*,
														 double,
														 double,
														 double,
														 int,
														 char*,
														 double*,
														 double*),
							 char			*DefaultRefRate,
							 char			*PrimSwapFreq,
							 char			*PrimSwapBasis,
							 char			*PrimRefRate,
							 int			 NumPrimCalibDates,
							 long			*ExerDatesPrim,
							 int			*DateSelecSpecifiersPrim,
							 char		   **CalibTenorsPrim,
							 long			 EndDatePrim,
							 double			*PrimStrikes,

							 cpd_diag_calib_param  *PrimParams,

							 char			*SecSwapFreq,
							 char			*SecSwapBasis,
							 char			*SecRefRate,
							 int			 NumSecCalibDates,
							 long			*ExerDatesSec,
							 int			*DateSelecSpecifiersSec,
							 char		   **CalibTenorsSec,
							 long			 EndDateSec,
							 double			*SecStrikes,

							 cpd_diag_calib_param  *SecParams,

							 short			 FixLambda,
							 short			 DoOneFactorEquivalent,
							 int			 NumLambdas,
							 double			*LambdaTimes,
							 double			*LambdaVals,
							 int			 OneOrTwoFactor,
							 double			 Alpha,
							 double			 Gamma,
							 double			 Rho,
							 int			*NumSigmas,						
							 double		   **SigmaTimes,
							 double		   **SigmaVals,

							 diag_calib_lm_params  *LmParams,

							 long			 Today,
							 char			*UndName,
							 int			 NumEventDates,
							 long			*EventDates,
							 long			 TableauRows,
							 long			 TableauCols,
							 char		  ***TableauStrings,
							 int		   **TableauMask,
							 long			 AuxWidth,
							 long			*AuxLen,
							 double		   **Aux,
							 /* End of Day Flags */
							 int			is_end_of_day_fixing,
							 int			is_end_of_day_payment,

							 int			 TreeOrMCOrMCoebAsInteger,

							 s_grfn_pricing_params  *PricingParams,
						
							 int			*ExerBoundaryOptimDateSpecifiers,
							 MCEBParams		*ExerBoundaryOptimParams,
							 int			*NumProducts,
							 double		  ***ProductVals,
							 long			*NumRowsExerBoundary,

							 cpd_calib_inst_data  *CalibrationInstrumentData,
							 
							 double		   **PrimCalibratedModelPrices,
							 double		   **SecCalibratedModelPrices );



#endif