/* ===============================================================================
   
   FILENAME:       srt_h_grfnflexicap.h

   PURPOSE:        Pricing of a FlexiCap with Grfn for a given model
					->Builds a tableau and prices the deal
					-> AutoCalibrate on Caps for that Deal

   =============================================================================== */

#ifndef SRT_H_GRFNFLEXICAP_H
#define SRT_H_GRFNFLEXICAP_H


void GrfnFlexiCapFraCell(
			String	    szCashBasis,
			String      szUndName,
			String      szRefRateCode,
			String      szFraCell);

Err NewGrfnFlexiCapMakeTableau(
			long              lNumEventDates,
			double            dStrike,
			double			  dEps,
			long              lMaxNumExercise,
			long              lNumAlreadyExercised,
			SrtReceiverType   eCapFloor,
			String	         szCashBasis,
			String           szUndName,
			String           szRefRateCode,
			long            *plNumRows,
			long            *plNumCols,
			GrfnCell    ***pppsTableau);

Err GrfnFlexiCapMakeTableau(
		long              lNumEventDates,
		double            dStrike,
		long              lMaxNumExercise,
		long              lNumAlreadyExercised,
		SrtReceiverType   eCapFloor,
		String	         szCashBasis,
		String           szUndName,
		String           szRefRateCode,
		long            *plNumRows,
		long            *plNumCols,
		GrfnCell    ***pppsTableau
		);

/* ---------------------------------------------------------------------------
	Prices a Flexi Cap in Grfn, building the Schedule and the Tableau
   --------------------------------------------------------------------------- */

Err SrtGrfnFlexiCapPrice(
	long            *plFixingDates,
	long            *plStartDates,
	long            *plPayDates,
	double          *pdCoverages,
	long              lNumDates,
	String           szRefRateCode,
	double 	          dStrike,
	double			  dEps,
	long              lMaxNumExercise,
	long              lNumAlreadyExercised,
	String           szCapFloor,
	String           szUndName,
	String         *pszGrfnParamNames,
	String         *pszGrfnParamValues,
	int               iNumGrfnParams, 
	double 	          dFullCapRealPrice,
	
/* Outputs from Grfn */	
	double 	        *pdFlexiCapPrice,
	double			*pdGrfnFullCapFloorPrice,
	long            *plNumRows,
	long            *plNumCols,
	GrfnCell 	 ***ppsGrfnTableau,
	double        ***pppdLastPath,
	long            *plNumAuxColumns,
	long          **pplAuxRangesLength,
	double      ***pppdAuxRanges
	)			
		;			


/* ---------------------------------------------------------------------------
	Calibrates on a set of Caps/Floors for a Flexi Cap
   --------------------------------------------------------------------------- */
Err SrtGrfnFlexiCapCalibrate(
	long            *plFixingDates,
	long            *plStartDates,
	long              lNumDates,
	String           szRefRateCode,
	double 	          dStrike,
	String           szCapFloor,
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
	double 	      **ppdCapletRealPrices,
	long            *plNumCaplets
);

#endif