/* ===============================================================================
   
   FILENAME:       srt_h_grfnknockinswap.h

   PURPOSE:        Pricing of a KnockIn Swap with Grfn for a given model
					-> Builds a tableau and prices the deal
					-> AutoCalibrate on Caps/Floors for that Deal

   =============================================================================== */

#ifndef SRT_H_GRFNKNOCKINSWAP_H
#define SRT_H_GRFNKNOCKINSWAP_H


/* ---------------------------------------------------------------------------
		Prices a KnockIn Swap in Grfn, building the Tableau internally
   --------------------------------------------------------------------------- */


Err SrtGrfnKnockInSwapPrice(
	long            *plKnockInFixingDates,
	long            *plKnockInIndexStartDates,
	long            *plKnockInSwapStartDates,
	double          *pdKnockInLevels,
	long              lNumKnockInDates,
	
	String           szIndexRefRate,
	String           szKnockUpOrDown,
	double            dCallSpreadWidth,
	SRT_Boolean           bStubPaidAtStart,

	long            *plFixedStartDates,
	long            *plFixedPayDates,
	double          *pdFixedCoupons,
	double          *pdNotionals,
	long              lNumFixedDates,
	String           szFixedBasis,
	String           szRecPay,
		
	long            *plFloatStartDates,
	long            *plFloatPayDates,
	long              lNumFloatDates,
	double            dFloatMargin,
	String           szFloatRefRate,
	
	String           szUndName,

	String         *pszGrfnParamNames,
	String         *pszGrfnParamValues,
	int               iNumGrfnParams, 
	
/* Outputs from Grfn */	
	double 	        *pdKnockInSwapPrice,
	long            *plNumRows,
	long            *plNumCols,
	GrfnCell 	***pppsGrfnTableau,
	double        ***pppdLastPath,
	long            *plNumAuxColumns,
	long          **pplAuxRangesLength,
	double      ***pppdAuxRanges
	)				;			


/* ---------------------------------------------------------------------------
	Calibrates on a set of Caps/Floors for a Flexi Cap
   --------------------------------------------------------------------------- */
Err SrtGrfnKnockInSwapCalibrate(
	long            *plIndexFixingDates,
	long            *plIndexStartDates,
	double          *pdKnockInLevels,
	long              lNumDates,
	String           szIndexRefRateCode,
	String           szKnockUpOrDown,
	
	String           szModelName,
	double            dTau,
	SRT_Boolean           bUseTwoFactor,
	double            dAlpha,
	double            dGamma,
	double            dRho,
	String         *pszGrfnParamNames,
	String         *pszGrfnParamValues,
	int               iNumGrfnParams, 
	String           szYieldCurveName,
	Err              (*pfGetVol)(double dStart, double dEnd, double dStrike, 
						double dForward, double dSpread, double *pdBsVol),
	String           szVolType,
/* Outputs from this calibration */	
	String           szUndName,
	double 	      **ppdCapletPrices,
	long            *plNumCaplets
	);

#endif