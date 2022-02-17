/* ========================================================================

  FILENAME:  srt_h_calibutils.h	

  AUTHORS:   Olivier Van Eyseren, Antoine Savine

  PURPOSE:   Provide a few utility functions that could help during 
             calibration:
               - automatic set up of dates
               - move from optimisation parameters to sig/tau
               - ...

  ========================================================================== */
#ifndef SRT_H_CALIBUTILS_H
#define SRT_H_CALIBUTILS_H

#ifdef __cplusplus
extern "C" {
#endif
/* -----------------------------------------------------------------------
              FUNCTIONS TO FIND SIGMA OR TAU DATES THAT ARE
                    RELEVANT FOR A SET OF INSTRUMENTS
   ----------------------------------------------------------------------- */

Err find_dates_from_deals(
		SwapDP        *sdparray, 
		StructType    *opt_type, 
		long          numInstruments, 
		SrtCcyParam   *ccy_param, 
		double        **sigDates, 
		long          *lNumSigmas,
		double        **tauDates,
		long          *lNumTaus);

Err find_sigdates_from_deals(
		SwapDP         *sdparray, 
		StructType     *opt_type, 
		long           numInstruments, 
		SrtCcyParam    *ccy_param, 
		double         **sigDates, 
		long           *numSigs);

Err find_taudates_from_deals(
		SwapDP         *sdparray,  
		long           numInstruments, 
		SrtCcyParam    *ccy_param, 
		double         **tauDates, 
		long           *numTaus);

/* ------------------------------------------------------------------------- */


/* -----------------------------------------------------------------------
              FUNCTIONS TO FIND A GOOD STARTING POINT FOR THE
              CALIBRATION ALGORITHM USING A SOBOL PRESAMPLING
   ----------------------------------------------------------------------- */

Err sobol_startingpoint(
	double        *param,
	double        **param_bound,
	long          numParams,
	long          npath,
	double        (*funcs)(double *),
	double        *criteria);

/* ------------------------------------------------------------------------- */


/* -----------------------------------------------------------------------
              FUNCTIONS TO MOVE FROMSIG-TAU TO PARAM
                         (AND VICE-VERSA) 
   ----------------------------------------------------------------------- */

Err from_sigtau_to_optparam( 
		double         **ppdSigmaValues,
		long               lNumSigmas,
		double         **ppdTauValues,
		long               lNumUsedTaus,
		long               lNumUsedBetas,
		long               lNumUsedOmegas,
		SrtCalibType       eCalibType,
		SrtMdlDim		   eModelDim,
		double           *pdOptParams			/* From [1] to [param_number] */
		);            

Err from_optparam_to_sigtau( 
		double           *pdOptParams,			/* From [1] to [param_number] */
		double         **ppdSigmaValues,
		long               lNumSigmas,
		double         **ppdTauValues,
		long               lNumTaus,
		SRT_Boolean            bFreezeTau,
		SRT_Boolean            bOneTau,
		double             *dFixedTau,
		SRT_Boolean            bFreezeBeta,
		SRT_Boolean            bOneBeta,
		double             dFixedBeta,
		SrtMdlDim	       eModelDim,
		SRT_Boolean            bFreezeOmega,
		SRT_Boolean            bOneOmega,
		double             dFixedOmega,
		SrtCalibType       eCalibType,
		double             dFixedAlpha,
		double             dFixedGamma,
		double		       dFixedRho);



/* -----------------------------------------------------------------------
              FUNCTIONS TO SET UP PARAMETERS NEEDED BY CALIBRATION
   ----------------------------------------------------------------------- */
/* -----------------------------------------------------------------------
	Check if parameters (pdOptimParams goes from [1] to [lNumParams]) are
	within the min/max ppdParamBounds ([0]: min; [1] : max)
   ----------------------------------------------------------------------- */

SRT_Boolean are_calib_parameters_within_band(
		double         *pdOptimParams, 
		long           lNumParams,
		double         **ppdParamBounds);

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
	From the calibration parameters (default or input on spreadsheet),
	gives the ppdParamBounds of the parameters in an array (for Sobol):
		[1][1..n] is the minimum value of the parameter
		[2][1..n] is the maximum value of the parameter
	----------------------------------------------------------------------- */


Err set_calib_parameters_bounds(
		SrtCalibParam    *psCalibParams,
		double         **ppdParamBounds,
		long               lNumSigmas,
		long               lNumUsedTaus,
		long               lNumUsedBetas,
		long               lNumUsedOmegas,
		SrtMdlDim	       eModelDim);

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
	Sets pdMarketWeights for minimisation criteria as in LM for SIMPLEX/ANNEALING 
				criteria = Sum(0<=i<n) [(y_i-y0_i)/pdMarketWeights_i]^2 
	For FRD prices, pdMarketWeights = pdMktVega * sqrt( n * dOptionsWeight ).
	For correlation, pdMarketWeights = n / sqrt( 1- dOptionsWeight )
	The global criteria is therefore:
				   dOptionsWeight^2 * 1/n * Sum(i) [(y_i-y0_i)/vega_i]^2
				+ (1-frd_weght)^2 * (1/N)^2 * Sum(j) [corr_i - mktcorr_i]^2
    ----------------------------------------------------------------------- */
Err define_criteria_weights( 
		long           lNumInstruments,
		long		   lNumTenors,
		SRT_Boolean		   bSmoothSigma,
		double         dOptionsWeight,
		double		   dCorrelWeight,
		double		   bSmoothSigmaWeight,
		double         *pdMktVega,                 /* From [0] to [lNumInstruments-1] */
		double         *pdMarketWeights);               /* From [1] to [lNumData] */


/* -----------------------------------------------------------------------
	Sets a vector of market target for minimisation criteria:
		if (index)<=lNumInstruments:
			set market pdMktPrice of caliration instrument 
		if (index)>lNumInstruments:
			set market correlation 
    ----------------------------------------------------------------------- */

Err define_market_targets( 
		long           lNumInstruments,
		SRT_Boolean		   bSmoothSigma,		   	
		double         *mkt_price,            /* From [0] to [lNumInstruments-1] */ 
		double         **ppdCorrelationMatrix,         /* From [0] to [lNumData-1] */
		long           lNumTenors,
		long           lNumData,
		double         *pdMarketTargets);           /* From [1] to [lNumData] */


/* ======================================================================= */
#ifdef PVMPI
#include "parallel.h"
#endif

/* ----------------------------------------------------------------------------
    
	  FOR QUICK TRANSFORMATIONS FORM A MODEL TS TO AN LGM ONE 

   ---------------------------------------------------------------------------- */

/* From a FULL model SigmaCurve to a FULL LGM one */
Err from_model_ts_to_LGM_ts(
		SrtMdlDim         eModelDim,
		double            **ppdModelSigmaCurve,
		long              lNumSigmas,
		String            szYieldCurveName,
		double            **ppdLgmSigmaCurve
		);

/* From a FULL (rough) LGM SigmaCurve to a FULL Model one (inverse of previous) */
Err from_LGM_ts_to_model_ts(
		SrtMdlDim         eModelDim,
		double            **ppdLgmSigmaCurve,
		long              lNumSigmas,
		String            szYieldCurveName,
		double            **ppdModelSigmaCurve);

/* -------------------------------------------------------------------

	CRITERIA FOR THE CALIBRATION TO THE HISTORICAL CORRELATION

  --------------------------------------------------------------------- */

Err srt_f_SquareRootSymMatrix(double **A, long dim, double **result);
Err srt_f_HellingerDist(double **HistSquare, int flag1, double **ModelCorr, int flag2, long dim, double *result);




#ifdef __cplusplus
}
#endif

#endif
