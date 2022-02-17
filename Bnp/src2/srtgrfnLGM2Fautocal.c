//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Sort Static level function to calibrate an LGM 2 factor model and then price a Grfn tableau using	//
//	either PDE, MC or MC with optimal ex boundary.														//
//																										//
//						Author:  Paul McCallum				Date:  29th Sept 2003						//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


#include "SrtAccess.h"
#include "CPDCalib.h"
#include "DiagCalibDLM.h"
#include "AmortMidatCalib.h"
#include "SrtGrfnLGM2FAutocal.h"



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//  Note the memory management policy in this routine. In particular the fact that we free the Product  //
//  Values when there has been an error but leave the calling routine to free in the absence of errors. //
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
							 double		   **SecCalibratedModelPrices )
{
	Err	 err = NULL;
	int  i, sum;
	int nUsedEventDates;
	double *vol_termstructure[2];
	double *tau_termstructure[2];

	double *temp_prod_vals = NULL;

	SrtUndPtr  und_ptr;

	int		notice_period;
	long	temp_date;


	//	Here we force the secondary calibration to be off when there are no intruments passed.
	//
	sum = 0;

	for (i=0; i<NumSecCalibDates; i++)
	{
		sum += DateSelecSpecifiersSec[i];
	}

	if (sum == 0) FixLambda = 1;

	
	//	Straight into calling the calibration function  //
	//
	err = cpd_calib_diagonal_dlm( YieldCurveName,
								  VolCubeName,
								  GetCashVol,
								  DefaultRefRate,
								  PrimSwapFreq,
								  PrimSwapBasis,
								  PrimRefRate,
								  NumPrimCalibDates,
								  ExerDatesPrim,
								  DateSelecSpecifiersPrim,
								  CalibTenorsPrim,
								  EndDatePrim,
								  PrimStrikes,
								  PrimParams,
								  SecSwapFreq,
								  SecSwapBasis,
								  SecRefRate,
								  NumSecCalibDates,
								  ExerDatesSec,
								  DateSelecSpecifiersSec,
								  CalibTenorsSec,
								  EndDateSec,
								  SecStrikes,
								  NULL,
								  SecParams,
								  FixLambda,
								  DoOneFactorEquivalent,
								  NumLambdas,
								  LambdaTimes,
								  LambdaVals,
								  NULL,
								  OneOrTwoFactor,
								  Alpha,
								  Gamma,
								  Rho,

								  0,
								  NULL,
								  NULL,
								  0,
								  NULL,
								  NULL,

								  NumSigmas,
								  SigmaTimes,
								  SigmaVals,
								  LmParams,
								  CalibrationInstrumentData );
	if (err)  return err;


		
	//  Initialise the Underlying and get the pointer  //
	//
	vol_termstructure[0] = *SigmaTimes;
	vol_termstructure[1] = *SigmaVals;
	
	for (i=0; i<*NumSigmas; i++)
	{
		(*SigmaTimes)[i] = (long) (Today + DAYS_IN_YEAR * (*SigmaTimes)[i] + 1.0e-08);
	}


	tau_termstructure[0] = LambdaTimes;
	tau_termstructure[1] = LambdaVals;

	// note how we are now storing tau values in LambdaVals
	for (i=0; i<NumLambdas; i++)
	{
		LambdaVals[i]  = 1.0 / LambdaVals[i];
		LambdaTimes[i] = (long) (Today + DAYS_IN_YEAR * LambdaTimes[i] + 1.0e-08);
	}


	err = SrtInitIRUnd( UndName,
						YieldCurveName,
						OneOrTwoFactor == 1 ? "LGM": "LGM2F", 
						*NumSigmas, 2, vol_termstructure, 
						NumLambdas, 2, tau_termstructure, 
						0.0, 
						Alpha, Gamma, Rho, 
						0,  0, 0, 0, 0, 0);

	if (err) 
	{
		if (*SigmaTimes)  free (*SigmaTimes);
		if (*SigmaVals)  free (*SigmaVals);

		return err;	
	}

	
	und_ptr = lookup_und (UndName);
	
	if (!und_ptr)
	{
		if (*SigmaTimes)  free (*SigmaTimes);
		if (*SigmaVals)  free (*SigmaVals);

		return "Couldn't initialise the underlying";	
	}



	//	Now for the sake of the output at West level and in order to determine how good the calibration is,	   //
	//	we now record the model prices. The market prices are already contained in CalibrationInstrumentData   // 
	//	returned by the calibration. Note the fact that we hard code the option type in the calls below. This  //
	//	is because it is nigh impossible to extract the types from the original code, which are hard coded to  //
	//	be "REC". If this should ever change, the code below will become out of sync. PMc  1stDec03			   //									   //
	//
	if (CalibrationInstrumentData)
	{
		*PrimCalibratedModelPrices = (double*) calloc (CalibrationInstrumentData->num_inst, sizeof (double));

		if (NumLambdas == 1)
		{
			for (i=0; i<CalibrationInstrumentData->num_inst; i++)
			{
				/* find the notice period */		
				notice_period = 0;
				temp_date = CalibrationInstrumentData->start_dates[i];

				while (temp_date >= CalibrationInstrumentData->exer_dates_long[i])
				{
					notice_period++;
					temp_date = add_unit (CalibrationInstrumentData->start_dates[i], -notice_period, SRT_BDAY, MODIFIED_SUCCEEDING);
				}

				notice_period--;

				err = europSwaption_clsdfrm( YieldCurveName,
											 und_ptr,
											 notice_period,
											 CalibrationInstrumentData->start_dates[i],
											 CalibrationInstrumentData->end_dates[i],
											 PrimSwapFreq,
											 PrimSwapBasis,
											 PrimRefRate,
											 0.,
											 CalibrationInstrumentData->long_strikes[i],
											 0.,
											 "REC",
											&(*PrimCalibratedModelPrices)[i] );
				if (err) 
				{
					if (*SigmaTimes)  free (*SigmaTimes);
					if (*SigmaVals)  free (*SigmaVals);

					srt_f_destroy_und (UndName);

					if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
					*PrimCalibratedModelPrices = NULL;			

					return err;
				}
			}
		}

		//	do same for secondary prices  //
		//
		*SecCalibratedModelPrices = (double*) calloc (CalibrationInstrumentData->num_insts, sizeof (double));

		if (NumLambdas == 1)
		{
			for (i=0; i<CalibrationInstrumentData->num_insts; i++)
			{
				/* find the notice period */		
				notice_period = 0;
				temp_date = CalibrationInstrumentData->start_datess[i];

				while (temp_date >= CalibrationInstrumentData->exer_dates_short[i])
				{
					notice_period++;
					temp_date = add_unit (CalibrationInstrumentData->start_datess[i], -notice_period, SRT_BDAY, MODIFIED_SUCCEEDING);
				}

				notice_period--;

				err = europSwaption_clsdfrm( YieldCurveName,
											 und_ptr,
											 notice_period,
											 CalibrationInstrumentData->start_datess[i],
											 CalibrationInstrumentData->end_datess[i],
											 SecSwapFreq,
											 SecSwapBasis,
											 SecRefRate,
											 0.,
											 CalibrationInstrumentData->short_strikes[i],
											 0.,
											 "REC",
											&(*SecCalibratedModelPrices)[i] );
				if (err) 
				{
					if (*SigmaTimes)  free (*SigmaTimes);
					if (*SigmaVals)  free (*SigmaVals);

					srt_f_destroy_und (UndName);

					if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
					*PrimCalibratedModelPrices = NULL;
					if (*SecCalibratedModelPrices)  free (*SecCalibratedModelPrices);
					*SecCalibratedModelPrices = NULL;

					return err;
				}
			}
		}
	}

	//  Price now using Tree/MC/MC optimal exercise boundary. We now have the provision to calculate	//
	//	with a Lambda term structure (12th Nov 2003) (case when NumLambdas != 1). To be sure, we still  //
	//	call old routine when there is only one value of Tau/Lambda. Also modified so that we call two  //
	//	different pde tau routines depending on status of the stability parameter which is now part		//
	//	of the pricing parameter structure (13th Nov 2003),												//
	//
	switch (TreeOrMCOrMCoebAsInteger)
	{
	case 0:
		if (NumLambdas == 1)
		{
			err = SrtGrfnLGM2Fpde( UndName,
								   NumEventDates,
								   EventDates,
								   TableauRows,
								   TableauCols,
								   TableauStrings,
								   TableauMask,
								   AuxWidth,
								   AuxLen,
								   Aux,
								   is_end_of_day_fixing,
								   is_end_of_day_payment,
								   PricingParams->num_steps_time,
								   PricingParams->num_steps_space,
								   NumProducts,
								  &temp_prod_vals );
		}
		else if (PricingParams->use_stabler_pde_with_tau_term_structure == 0)
		{
			err = SrtGrfnLGM2FTaupde( UndName,
									  NumEventDates,
									  EventDates,
									  TableauRows,
									  TableauCols,
									  TableauStrings,
									  TableauMask,
									  AuxWidth,
									  AuxLen,
									  Aux,
									  is_end_of_day_fixing,
									  is_end_of_day_payment,
									  PricingParams->num_steps_time,
									  PricingParams->num_steps_space,
									  NumProducts,
									 &temp_prod_vals );
		}
		else
		{
			err = SrtGrfnLGM2FTaupde2( UndName,
									   NumEventDates,
									   EventDates,
									   TableauRows,
									   TableauCols,
									   TableauStrings,
									   TableauMask,
									   AuxWidth,
									   AuxLen,
									   Aux,
									   is_end_of_day_fixing,
									   is_end_of_day_payment,
									   PricingParams->num_steps_time,
									   PricingParams->num_steps_space,
									   NumProducts,
									  &temp_prod_vals );
		}

		if (err) 
		{
			if (*SigmaTimes)  free (*SigmaTimes);
			if (*SigmaVals)  free (*SigmaVals);
			
			srt_f_destroy_und (UndName);

			if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
			*PrimCalibratedModelPrices = NULL;
			if (*SecCalibratedModelPrices)  free (*SecCalibratedModelPrices);
			*SecCalibratedModelPrices = NULL;

			if (temp_prod_vals)
			{
				free_dvector (temp_prod_vals, 0, *NumProducts - 1);
				temp_prod_vals = NULL;
			}
			return err;	
		}

		*ProductVals = dmatrix (0, *NumProducts-1, 0, 1);

		if (!*ProductVals)
		{
			if (*SigmaTimes)  free (*SigmaTimes);
			if (*SigmaVals)  free (*SigmaVals);
			
			srt_f_destroy_und (UndName);

			if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
			*PrimCalibratedModelPrices = NULL;
			if (*SecCalibratedModelPrices)  free (*SecCalibratedModelPrices);
			*SecCalibratedModelPrices = NULL;

			if (temp_prod_vals)
			{
				free_dvector (temp_prod_vals, 0, *NumProducts - 1);
				temp_prod_vals = NULL;
			}

			return "Not enough memory to return Product Values";	
		}

		// now transfer data
		for (i=0; i<*NumProducts; i++)
		{
			(*ProductVals)[i][0] = temp_prod_vals[i];
		}

		free_dvector (temp_prod_vals, 0, *NumProducts - 1);
		temp_prod_vals = NULL;
		break;

	case 1:
		if (NumLambdas == 1)
		{
			err = SrtGrfnLGM2FMC( UndName,
								  NumEventDates,
								  EventDates,
								  TableauRows,
								 &TableauCols,
								  TableauStrings,
								  TableauMask,
								  AuxWidth,
								  AuxLen,
								  Aux,
								  is_end_of_day_fixing,
								  is_end_of_day_payment,
								  PricingParams->num_paths,
								  PricingParams->do_pecs,
								  PricingParams->do_jump,
								  ProductVals );
		}
		else
		{
			err = SrtGrfnLGM2FMClambda( UndName,
										NumEventDates,
										EventDates,
										TableauRows,
									   &TableauCols,
										TableauStrings,
										TableauMask,
										AuxWidth,
										AuxLen,
										Aux,
										is_end_of_day_fixing,
										is_end_of_day_payment,
										PricingParams->num_paths,
										PricingParams->do_pecs,
										PricingParams->do_jump,
										ProductVals );
		}

		if (err) 
		{
			if (*SigmaTimes)  free (*SigmaTimes);
			if (*SigmaVals)  free (*SigmaVals);
			
			srt_f_destroy_und (UndName);

			if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
			*PrimCalibratedModelPrices = NULL;
			if (*SecCalibratedModelPrices)  free (*SecCalibratedModelPrices);
			*SecCalibratedModelPrices = NULL;

			if (*ProductVals)
			{
				free_dmatrix (*ProductVals, 0, TableauCols - 1, 0, 1);
				*ProductVals = NULL;
			}
			return err;	
		}
		
		// no need to fill out the output array
		*NumProducts = TableauCols;
		break;

	case 2:
		if (NumLambdas == 1)
		{
			err = SrtGrfnLGM2FMCEB( UndName,
									NumEventDates,
									EventDates,
								   &nUsedEventDates,
									ExerBoundaryOptimDateSpecifiers,
									NULL,
									ExerBoundaryOptimParams,
									NumRowsExerBoundary,
									TableauRows,
								   &TableauCols,
									TableauStrings,
									TableauMask,
									AuxWidth,
									AuxLen,
									Aux,
									is_end_of_day_fixing,
									is_end_of_day_payment,
									PricingParams->num_paths,
									PricingParams->do_pecs,
									PricingParams->do_jump,
									ProductVals );
		}
		else
		{
			err = SrtGrfnLGM2FMCEBlambda( UndName,
										  NumEventDates,
										  EventDates,
										  &nUsedEventDates,
										  ExerBoundaryOptimDateSpecifiers,
										  NULL,
										  ExerBoundaryOptimParams,
										  NumRowsExerBoundary,
										  TableauRows,
										 &TableauCols,
										  TableauStrings,
										  TableauMask,
										  AuxWidth,
										  AuxLen,
										  Aux,
										  is_end_of_day_fixing,
										  is_end_of_day_payment,
										  PricingParams->num_paths,
										  PricingParams->do_pecs,
										  PricingParams->do_jump,
										  ProductVals );
		}

		if (err) 
		{
			if (*SigmaTimes)  free (*SigmaTimes);
			if (*SigmaVals)  free (*SigmaVals);
			
			srt_f_destroy_und (UndName);

			if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
			*PrimCalibratedModelPrices = NULL;
			if (*SecCalibratedModelPrices)  free (*SecCalibratedModelPrices);
			*SecCalibratedModelPrices = NULL;

			if (*ProductVals)
			{
				free_dmatrix (*ProductVals, 0, *NumRowsExerBoundary - 1, 0, 2 + ExerBoundaryOptimParams->iNbIndex);
				*ProductVals = NULL;
			}
			return err;	
		}
		
		// no need to fill out the output array
		*NumProducts = TableauCols;
		break;

	default:
		if (*SigmaTimes)  free (*SigmaTimes);
		if (*SigmaVals)  free (*SigmaVals);
		
		srt_f_destroy_und (UndName);

		if (*PrimCalibratedModelPrices)  free (*PrimCalibratedModelPrices);
		*PrimCalibratedModelPrices = NULL;
		if (*SecCalibratedModelPrices)  free (*SecCalibratedModelPrices);
		*SecCalibratedModelPrices = NULL;

		return "Wrong value of TreeOrMCOrMCoeb";
	}


	return err;
}
