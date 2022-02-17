//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Sort Static level function to calibrate an LGM SV model having either 1 or 2 LGM factors and 1 SV	//
//	factor, initialise a SV underlying then price a Grfn tableau using either PDE, MC or MC with		//
//	optimal ex boundary.																				//
//																										//
//						Author:  Paul McCallum				Date:  13th Oct 2003						//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


#include "SrtAccess.h"
#include "LGMSVCalibApprox.h"		// this includes CPDcalib.h, DiagCalibDLM.h	and LGMSVClosedformApprox.h...		
#include "DiagCalibDLMSV.h"

#include "SrtGrfnLGM2FAutocal.h"
#include "SrtGrfnLGMSVAutocal.h"



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//  Note the memory management policy in this routine. In particular the fact that we free the Product  //
//  Values when there has been an error but leave the calling routine to free in the absence of errors. //
//	Also how we have the capability to reassign smile values in the same arrays as the trial values by	//
//	repointing the pointers. These must be freed here since any possible client callnig routine may		//
//	not pass preassigned pointers and therefore expect all memory management to happen here.			//																			//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
Err  GrfnLGMSVAutocalCaller( char			*YieldCurveName,
							 char			*VolCubeName,
							 Err		   (*pGetCashVol)(char*,
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
							 double		   **PrimStrikes,

							 cpd_diag_calib_param  *pPrimParams,

							 char			*SecSwapFreq,
							 char			*SecSwapBasis,
							 char			*SecRefRate,
							 int			 NumSecCalibDates,
							 long			*ExerDatesSec,
							 int			*DateSelecSpecifiersSec,
							 char		   **CalibTenorsSec,
							 long			 EndDateSec,
							 double		   **SecStrikes,

							 cpd_diag_calib_param  *pSecParams,
							 
							 LGMSV_CalibParams *pCalibParams,
							 
							 int			 OneOrTwoFactor,
							 double			*pLambda,							 
							 int			*pNumSetsSmileParams,
							 double		   **pSmileTimes,
							 double		   **pAlphaSmileVals,
							 double		   **pLambdaSmileVals,
							 double		   **pRhoSmileVals,
							 double		   **pRho2SmileVals,
							 double			 TStar,
							 double			 LGMAlpha,
							 double			 LGMGamma,
							 double			 LGMRho,

							 LGMSV_NumerParams	*pNumericalParams,

							 int			*pNumSigmas,						
							 double		   **pSigmaTimes,
							 double		   **pSigmaVals,

							 diag_calib_lm_params  *pLmParams,

							 long			 Today,
							 char			*UndName,
							 LGMSVParam		*pAlgorithmSpecification,
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

							 s_grfn_SV_pricing_params  *pPricingParams,
						
							 int			*ExerBoundaryOptimDateSpecifiers,
							 MCEBParams		*pExerBoundaryOptimParams,
							 int			*pNumProducts,
							 double		  ***pProductVals,
							 long			*pNumRowsExerBoundary,

							 cpd_calib_inst_data  *pCalibrationInstrumentData,
							 
							 double		   **pPrimCalibratedModelPrices,
							 double		   **pSecCalibratedModelPrices )
{
	Err	 err = NULL;
	int  i, sum;
	int nUsedEventDates;
	
	double **vol_termstructure;
	double  *p_single_tau_value;
	double **tau_termstructure;
	double **smile_termstructure;
	double   tau;

	double  *calib_alpha_smile_vals = NULL, *calib_lambda_smile_vals = NULL,
		    *calib_rho_smile_vals = NULL, *calib_rho2_smile_vals = NULL;

	double  *temp_prod_vals = NULL;

	SrtUndPtr  und_ptr;
	
	long	 und_ticker;
	char	 full_und_name[256];

	short	 do_optimisation = TreeOrMCOrMCoebAsInteger == 1 ? 0 : 1;
	

	memset (full_und_name, 0, 256);

	//	Here we force the secondary calibration to be off when there are no instruments passed. //
	//
	sum = 0;

	for (i=0; i<NumSecCalibDates; i++)
	{
		sum += DateSelecSpecifiersSec[i];
	}

	if (sum == 0) pCalibParams->fix_lambda = 1;


	//	Straight into calling the calibration function  //
	//
	err = cpd_calib_diagonal_LGMSV_new_dlm( YieldCurveName,
											VolCubeName,
											pGetCashVol,
											DefaultRefRate,

											PrimSwapFreq,
											PrimSwapBasis,
											PrimRefRate,
											NumPrimCalibDates,
											ExerDatesPrim,
											DateSelecSpecifiersPrim,
											CalibTenorsPrim,
											EndDatePrim,
											PrimStrikes[0],
											PrimStrikes[1],
											PrimStrikes[2],
											pPrimParams,

											SecSwapFreq,
											SecSwapBasis,
											SecRefRate,
											NumSecCalibDates,
											ExerDatesSec,
											DateSelecSpecifiersSec,
											CalibTenorsSec,
											EndDateSec,
											SecStrikes[0],
											SecStrikes[1],
											SecStrikes[2],
											NULL,
											pSecParams,

											pCalibParams,

											OneOrTwoFactor,
											pLambda,
										   *pNumSetsSmileParams,
										   *pSmileTimes,
										   *pAlphaSmileVals,
										   *pLambdaSmileVals,
										   *pRhoSmileVals,
											TStar,
											LGMAlpha,
											LGMGamma,
											LGMRho,
										   *pRho2SmileVals,
											pNumericalParams,

											pNumSigmas,
											pSigmaTimes,
											pSigmaVals,
										   &calib_alpha_smile_vals,
										   &calib_lambda_smile_vals,
										   &calib_rho_smile_vals,
										   &calib_rho2_smile_vals,
											pLmParams,
											pCalibrationInstrumentData );
		
	if (err)
	{	
		// we would like to think that this memory is freed in routine above, however just in case... 
		if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;
		if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

		if (calib_alpha_smile_vals)  free (calib_alpha_smile_vals);
		if (calib_lambda_smile_vals)  free (calib_lambda_smile_vals);
		if (calib_rho_smile_vals)  free (calib_rho_smile_vals);
		if (calib_rho2_smile_vals)  free (calib_rho2_smile_vals);

		return err;
	}

	
	//	Now set up arrays we need for initialising the underlying. Note that we could have done what we did		   //
	//	in the previous non-SV routine. ie had 2 x 2 row arrays of double pointers and point them to the		   //
	//	Lambda and Sigma data, albeit this time callocing the smile row size because of runtime 1 or 2 factor	   //
	//	decision. However it is safer practice and is not much more time consuming, for small arrays at least,	   //
	//	to copy values before passing them to a C routine (which of course could modify the values unbeknownigly)  //
	//
	vol_termstructure	= dmatrix (0, 1, 0, *pNumSigmas - 1);
	smile_termstructure = dmatrix (0, 3 + OneOrTwoFactor - 1, 0, *pNumSigmas - 1);
		
	for (i=0; i<*pNumSigmas; i++)
	{
		(*pSigmaTimes)[i] = (long) (Today + DAYS_IN_YEAR * (*pSigmaTimes)[i] + 1.0e-08);

		vol_termstructure[0][i] = (*pSigmaTimes)[i];
		vol_termstructure[1][i] = (*pSigmaVals)[i];

		smile_termstructure[0][i] = (*pSigmaTimes)[i];
		smile_termstructure[1][i] = calib_alpha_smile_vals[i];
		smile_termstructure[2][i] = calib_lambda_smile_vals[i];
		smile_termstructure[3][i] = calib_rho_smile_vals[i];
	}

	if (OneOrTwoFactor == 2)
	{
		for (i=0; i<*pNumSigmas; i++)
		{
			smile_termstructure[4][i] = calib_rho2_smile_vals[i];
		}
	}

	tau = 1.0 / *pLambda;
	p_single_tau_value = &tau;
	tau_termstructure  = &p_single_tau_value;
	
	
	//  Before setting up underlying, take this opportunity to repoint the external smile value arrays and free     //
	//	the local ones. Believe that where we have originally passed a single set of values and not recalibrated,	//
	//	or where we	calibrate but do not have enough information, there may well be duplicate sets of smile values. //		      //
	//	Note that we are now pointing pSmileTimes to *pSigmaTimes so there we must be careful to NULL after we		//
	//	free one of them so as not to try to free for a second time! Eg the error handling we have below is a 		//
	//	perfect example - the first one frees, the second freeing attempt is academic.								//
	//
	*pNumSetsSmileParams = *pNumSigmas;
									   
	if (*pSmileTimes)  free (*pSmileTimes);
	*pSmileTimes = *pSigmaTimes;

	if (*pAlphaSmileVals)  free (*pAlphaSmileVals);
	*pAlphaSmileVals = calib_alpha_smile_vals;

	if (*pLambdaSmileVals)  free (*pLambdaSmileVals);
	*pLambdaSmileVals = calib_lambda_smile_vals;

	if (*pRhoSmileVals)  free (*pRhoSmileVals);
	*pRhoSmileVals = calib_rho_smile_vals;

	if (*pRho2SmileVals)  free (*pRho2SmileVals);
	*pRho2SmileVals = calib_rho2_smile_vals;

	
	
	//  Set up the underlying and check validity  //
	//
	err = SrtInitIRMSVUnd( UndName,
						   YieldCurveName,
						   OneOrTwoFactor,
						   vol_termstructure,
						  *pNumSigmas,
						   2,
						   tau_termstructure,
						   1,
						   1,
						   LGMAlpha,
						   LGMGamma,
						   LGMRho,
						   smile_termstructure,
						  *pNumSigmas,
						   3 + OneOrTwoFactor,
						   TStar );

	// want to free these arrays now that we are done with them (the values are copied into the underlying) 

	free_dmatrix (vol_termstructure, 0, 1, 0, *pNumSigmas - 1);
	free_dmatrix (smile_termstructure, 0, 3 + OneOrTwoFactor - 1, 0, *pNumSigmas - 1);

	if (err)
	{
		if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
		if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

		if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
		if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
		if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
		if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

		return err;
	}

	
	und_ptr = lookup_und (UndName);
	
	if (!und_ptr)
	{
		if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
		if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

		if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
		if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
		if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
		if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

		return "Couldn't initialise the underlying";	
	}

	und_ticker = get_underlying_ticker (und_ptr);
	sprintf (full_und_name, "%s.%d", UndName, und_ticker);


	//	Now for the sake of the output at West level and in order to determine how good the calibration is,	   //
	//	we now record the model prices. The market prices are already contained in CalibrationInstrumentData   // 
	//	returned by the calibration. Note the fact that we hard code the option type in the calls below. This  //
	//	is because it is nigh impossible to extract the types from the original code, which are hard coded to  //
	//	be "REC". If this should ever change, the code below will become out of sync. PMc  1stDec03			   //									   //
	//

	if (pCalibrationInstrumentData)
	{
		*pPrimCalibratedModelPrices = (double*) calloc (pCalibrationInstrumentData->num_inst, sizeof (double));

		for (i=0; i<pCalibrationInstrumentData->num_inst; i++)
		{
			err = LGMSVOptionApprox_new( full_und_name, //und_ptr,
										 YieldCurveName,
										 PrimRefRate,
										 PrimSwapFreq,
										 PrimSwapBasis,
										 pCalibrationInstrumentData->exer_dates_long[i],
										 pCalibrationInstrumentData->start_dates[i],
										 pCalibrationInstrumentData->end_dates[i],
										 0, 0,
										 0, 0, 0,
										 1,
										&pCalibrationInstrumentData->long_strikes[i],
										-1,													// equiv to "REC"
										 pNumericalParams,
										&(*pPrimCalibratedModelPrices)[i] );
			if (err) 
			{
				if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
				if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

				if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
				if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
				if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
				if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

				if (*pPrimCalibratedModelPrices)  free (*pPrimCalibratedModelPrices);
				*pPrimCalibratedModelPrices = NULL;			

				return err;
			}
		}	

		//	do same for secondary prices  //
		//
		*pSecCalibratedModelPrices = (double*) calloc (pCalibrationInstrumentData->num_insts, sizeof (double));

		for (i=0; i<pCalibrationInstrumentData->num_insts; i++)
		{
			err = LGMSVOptionApprox_new( full_und_name, //und_ptr,
										 YieldCurveName,
										 SecRefRate,
										 SecSwapFreq,
										 SecSwapBasis,
										 pCalibrationInstrumentData->exer_dates_short[i],
										 pCalibrationInstrumentData->start_datess[i],
										 pCalibrationInstrumentData->end_datess[i],
										 0, 0,
										 0, 0, 0,
										 1,
										&pCalibrationInstrumentData->short_strikes[i],
										-1,													// equiv to "REC"
										 pNumericalParams,
										&(*pSecCalibratedModelPrices)[i] );
			if (err) 
			{
				if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
				if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

				if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
				if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
				if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
				if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

				if (*pPrimCalibratedModelPrices)  free (*pPrimCalibratedModelPrices);
				*pPrimCalibratedModelPrices = NULL;
				if (*pSecCalibratedModelPrices)  free (*pSecCalibratedModelPrices);
				*pSecCalibratedModelPrices = NULL;

				return err;
			}
		}
	}

	//  Price now using Tree/MC/MC optimal exercise boundary  //
	//
	switch (TreeOrMCOrMCoebAsInteger)
	{
	case 0:
		err = SrtGrfnLGMSVpde( UndName,
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
							  *pAlgorithmSpecification,
							   pPricingParams->standard_params.num_steps_time,
							   pPricingParams->standard_params.num_steps_space,
						       pPricingParams->num_steps_vol,
							   pPricingParams->num_steps_phi,

							   pNumProducts,
							  &temp_prod_vals );
		if (err) 
		{
			if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
			if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

			if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
			if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
			if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
			if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

			if (*pPrimCalibratedModelPrices)  free (*pPrimCalibratedModelPrices);
			*pPrimCalibratedModelPrices = NULL;
			if (*pSecCalibratedModelPrices)  free (*pSecCalibratedModelPrices);
			*pSecCalibratedModelPrices = NULL;

			if (temp_prod_vals)
			{
				free_dvector (temp_prod_vals, 0, *pNumProducts - 1);
				temp_prod_vals = NULL;
			}
			return err;	
		}


		// now transfer data
		*pProductVals = dmatrix (0, *pNumProducts-1, 0, 1);

		for (i=0; i<*pNumProducts; i++)
		{
			(*pProductVals)[i][0] = temp_prod_vals[i];
		}

		free_dvector (temp_prod_vals, 0, *pNumProducts - 1);
		temp_prod_vals = NULL;
		break;

	case 1:
	case 2:
		// nb we now do Optimisation case in same routine
		err = SrtGrfnLGMSVMC( UndName,
							  NumEventDates,
							  EventDates,
							 &nUsedEventDates,
							  TableauRows,
							 &TableauCols,
							  TableauStrings,
							  TableauMask,
							  AuxWidth,
							  AuxLen,
							  Aux,
							  is_end_of_day_fixing,
							  is_end_of_day_payment,

							 *pAlgorithmSpecification,
							  do_optimisation,
							  ExerBoundaryOptimDateSpecifiers,
							  NULL,
							  pExerBoundaryOptimParams,
							  pNumRowsExerBoundary,
							  pPricingParams->standard_params.num_steps_time_mc,
							  pPricingParams->standard_params.num_paths,

							  pNumProducts,
							  pProductVals );
		if (err) 
		{
			if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
			if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

			if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
			if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
			if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
			if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

			if (*pPrimCalibratedModelPrices)  free (*pPrimCalibratedModelPrices);
			*pPrimCalibratedModelPrices = NULL;
			if (*pSecCalibratedModelPrices)  free (*pSecCalibratedModelPrices);
			*pSecCalibratedModelPrices = NULL;

			if (*pProductVals)
			{
				if (do_optimisation)
				{
					free_dmatrix (*pProductVals, 0, *pNumRowsExerBoundary - 1, 0, 2 + pExerBoundaryOptimParams->iNbIndex);
				}
				else
					free_dmatrix (*pProductVals, 0, *pNumProducts - 1, 0, 2);

				*pProductVals = NULL;
			}
			return err;	
		}
		break;

	default:
		if (*pSigmaTimes)  free (*pSigmaTimes);  *pSigmaTimes = NULL;  *pSmileTimes = NULL;
		if (*pSigmaVals)  free (*pSigmaVals);  *pSigmaVals = NULL;

		if (*pAlphaSmileVals)  free (*pAlphaSmileVals);  *pAlphaSmileVals = NULL;
		if (*pLambdaSmileVals)  free (*pLambdaSmileVals);  *pLambdaSmileVals = NULL;
		if (*pRhoSmileVals)  free (*pRhoSmileVals);  *pRhoSmileVals = NULL;
		if (*pRho2SmileVals)  free (*pRho2SmileVals);  *pRho2SmileVals = NULL;

		if (*pPrimCalibratedModelPrices)  free (*pPrimCalibratedModelPrices);
		*pPrimCalibratedModelPrices = NULL;
		if (*pSecCalibratedModelPrices)  free (*pSecCalibratedModelPrices);
		*pSecCalibratedModelPrices = NULL;

		return "Wrong value of TreeOrMCOrMCoeb";
	}


	return err;
}
