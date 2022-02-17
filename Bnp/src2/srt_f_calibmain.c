/* ===============================================================================
   
   FILENAME:  srt_f_calibmain                                     
  
   PURPOSE	: Performs the calibration of interest rate model.

			  If the FixedPoint algorithm is not chosen, the calibration is done
			  immediately with the real model (call to srt_f_calibcore)
			  Otherwise, the calibration is done using a FIXED POINT ALGORITHM principle.

			  The complex model price is approximated by a rough LGM price (including
			  a rough Beta correction)
			  The Algorithm works as follows:
				- the approximate parameters LGM model is calibrated,
				  using the required method (LM,...)
				- the exact error in prices is computed with the real model
				- the price target is corrected by the error (symmetry)
				- the approximate LGM model is re-calibrated on those new prices
				- ...
			  Please note that the correlation is calibrated using LGM (exact)
			  and is not modified after


   =============================================================================== */

 
#include "srt_h_all.h"
#include "grf_h_public.h"
#include "srt_h_grfclsdfrm.h"
#include "utallhdr.h"
#include "srt_h_calib.h" 

#ifdef PVMPI
#include "parallel.h"
#endif

Err srt_f_calib_main(
		SrtGrfnParam     *psGrfnParams,
		SrtMdlType        eModelType,
		SrtMdlDim		  eModelDim,
		SwapDP            *psSwapDp,                    
		double            *pdStrike,
		double            *pdBondStrike,
		StructType 		  *peOptionType,
		SrtReceiverType   *peRecPay,
		double            *pdPrice,
		double            *pdVega,
		String            *pszRefRateCode,
		long               lNumInstruments,
		double            *dFraMaturities,
		double            **ppdCorrelationMatrix,
		long              lNumTenors,
		SrtUndPtr         sUndPtr,
		SrtCalibParam     *psCalibParams,
		double            **ppdSigmaCurve,
		long              lNumSigmas,
		long              lNumSigmaCols,
		double            **ppdTauCurve,
		long              lNumTaus,
		long              lNumTauCols,
		double            *pdChiSquare)
{
SrtMdlType            eInitialModelType;
String               szYieldCurveName;
long                  lInitialNumIter;
SrtCalibAlgoType      eInitialAlgoType;
double            **ppdLgmSigmaCurve = NULL;
Err                    szErr = NULL;
#ifdef PVMPI
int nb_cpus;
int* cpus;
#endif



/* Send a message warning the user: calibration starts */
#ifdef PVMPI
	GetCPUs(&nb_cpus, &cpus);
#endif
	smessage(" ");
	smessage("Starts Calibration");
#ifdef PVMPI
	smessage("with %d processors", nb_cpus);
#endif
	smessage(" ");
	
/* If FIXED POINT algorithm is required, LGM Starting point should be there */
	if (psCalibParams->eAlgoType == FIXED_POINT)
	{
		psCalibParams->bLgmStartPoint = SRT_YES;
	}

/* If Model is LGM, no need to do LGM starting point */
	if (eModelType == LGM)
	{
		psCalibParams->bLgmStartPoint = SRT_NO;
	}

/* LGM starting point */
	if (psCalibParams->bLgmStartPoint == SRT_YES )
	{
	/* Send message */
		smessage ("Calibrating LGM first");
		smessage ("");
		
	/* Replace the model type in the underlying */
		eInitialModelType = eModelType;
		eModelType = LGM;
		set_irund_mdltype(sUndPtr,eModelType);

	/* Allocate Memory for a LGM Sigma Curve */
		ppdLgmSigmaCurve = dmatrix (0, lNumSigmaCols-1, 0, lNumSigmas-1);
		
	/* Gets the yield curve name from the underlying */
		szYieldCurveName = get_discname_from_underlying(sUndPtr);
	
	/* Rebuilds a TS by a quick Normal approximation of the Volatility */
		szErr =  from_model_ts_to_LGM_ts(
					eModelDim,
					ppdSigmaCurve,
					lNumSigmas,
					szYieldCurveName,
					ppdLgmSigmaCurve);
		if (szErr)
		{
			free_dmatrix (ppdLgmSigmaCurve , 0, lNumSigmaCols-1, 0, lNumSigmas-1);
			return szErr;
		}

	/* Replace the number of iterations and the Algorithm */
		eInitialAlgoType = psCalibParams->eAlgoType;
		psCalibParams->eAlgoType = psCalibParams->eLgmAlgoType;
		lInitialNumIter = psCalibParams->lNumIter;
		psCalibParams->lNumIter = psCalibParams->lLgmNumIter;
		
	/* Calibrate LGM */
		szErr = srt_f_calib_core(
					psGrfnParams,    
					eModelType,
					eModelDim,
					psSwapDp,                    
					pdStrike,
					pdBondStrike,
					peOptionType,
					peRecPay,
					pdPrice,
					pdVega,
					pszRefRateCode,
					lNumInstruments,
					dFraMaturities,
					ppdCorrelationMatrix,	
					lNumTenors,
					sUndPtr,
					psCalibParams,
					ppdLgmSigmaCurve,
					lNumSigmas,
					lNumSigmaCols,
					ppdTauCurve,
					lNumTaus,
					lNumTauCols,
					pdChiSquare);

		if (szErr)
		{
			free_dmatrix (ppdLgmSigmaCurve , 0, lNumSigmaCols-1, 0, lNumSigmas-1);
			return szErr;
		}
	
	/* Reinstore the initial model type */
		set_irund_mdltype(sUndPtr,eInitialModelType);
		eModelType = eInitialModelType;

	/* Force FIXED correlation: it has already been calibrated */
		if (eModelDim == TWO_FAC)
		{
			psCalibParams->eCalibType = FIXED_CALIB;
		}
		psCalibParams->dOptionsWeight = 1.00;
		lNumTenors = 0;

	/* Rebuilds a pseudo-equivalent TS from the calibrated LGM (inverse of approximation) */
		szErr =  from_LGM_ts_to_model_ts(
					eModelDim,
					ppdLgmSigmaCurve,
					lNumSigmas,
					szYieldCurveName,
					ppdSigmaCurve);
		if (szErr)
		{
			free_dmatrix (ppdLgmSigmaCurve , 0, lNumSigmaCols-1, 0, lNumSigmas-1);
			return szErr;
		}

	/* Reinstates the number of iterations and the Algorithm */
		psCalibParams->eAlgoType = eInitialAlgoType;
		psCalibParams->lNumIter = lInitialNumIter;
		
	/* Frees the memory allocated for the LGM Term Struct */
		free_dmatrix (ppdLgmSigmaCurve , 0, lNumSigmaCols-1, 0, lNumSigmas-1);
		
	} /* END of LGM starting point */
	
	
	if (psCalibParams->eAlgoType != FIXED_POINT)
	{

	/* Warn the user: the real calibration is going to start */
		smessage ("Calibrating the Real Model Now");
		smessage ("");

	/* Main calibration, with the real model type */
		szErr = srt_f_calib_core(
					psGrfnParams,    
					eModelType,
					eModelDim,
					psSwapDp,                    
					pdStrike,
					pdBondStrike,
					peOptionType,
					peRecPay,
					pdPrice,
					pdVega,
					pszRefRateCode,
					lNumInstruments,
					dFraMaturities,
					ppdCorrelationMatrix,	
					lNumTenors,
					sUndPtr,
					psCalibParams,
					ppdSigmaCurve,
					lNumSigmas,
					lNumSigmaCols,
					ppdTauCurve,
					lNumTaus,
					lNumTauCols,
					pdChiSquare);
		if (szErr)
		{
			return szErr;
		}
	}
	else
	{
/* Treats the FIXED POINT algorithm separately */
	}


	return NULL;	
	
} /* END Err srt_f_calib_fixed_point(...)  */
