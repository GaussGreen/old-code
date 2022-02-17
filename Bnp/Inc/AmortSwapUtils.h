#ifndef __AMORTSWAP_UTILS_H
#define __AMORTSWAP_UTILS_H


Err ComputeAmortSwapDiagonalIRRs(char *cYCname, char *cVCname,
					   char *cRefRname,
					   long StartDate,
					   long EndDate,
					   SrtCompounding srtFreq,
					   SrtBasisCode srtBasis,
					   long	lNFixNot,
					   double *dFixNotionals,
					   double *dFixRates,
					   double *dExerciseFees,
					   long	lNFloatNot,
					   double *dFloatNotionals,
					   double *dMargins,
					   double *dIRRs,
					   int UseVol);

Err ComputeAmortSwapDiagonalIRRsNew(char *cYCname, 
								 char *cVCname,

								 char *cRefRname,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,

								 long	lNFix,
								 long	*lFixStartDates,
								 long	*lFixEndDates,
								 double *dFixCoverages,
								 double *dFixRates,
								 double *dFixNotionals,
								 double *dExerciseFees,
								 
								 long	lNFloat,
								 long	*lFloatStartDates,
								 long	*lFloatEndDates,
								 double *dFloatCoverages,
								 double *dFloatMargins,
								 double *dFloatSpreads,
								 double *dFloatNotionals,

								 double *dIRRs,
								 int UseVol);

Err ComputeAmortSwapDiagonalIRRsNew2(char *cYCname, 
								 char *cVCname,

								 char *cRefRname,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,

								 long	lNFix,
								 long	*lFixStartDates,
								 long	*lFixEndDates,
								 double *dFixCoverages,
								 double *dFixRates,
								 double *dFixNotionals,
								 double *dExerciseFees,
								 
								 long	lNFloat,
								 long	*lFloatStartDates,
								 long	*lFloatEndDates,
								 double *dFloatCoverages,
								 double *dFloatMargins,
								 double *dFloatSpreads,
								 double *dFloatNotionals,

								 double *dIRRs,
								 int UseVol);

Err		ComputeSwapRateFromIRR(char *cYCname, char *cVCname, long lToday,
						char *cRefRname,
						long StartDate,
						long EndDate,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						double IRR,
						int UseVol,
						double *dSwapRate);

Err		ComputeSwapRateFromIRR2(char *cYCname, char *cVCname, long lToday,
					    char *cFreq, char *cBasis, char *cRefRname,
						int shortstub,
						int floatshortstub,
						int iNFixPayDates,
						long *lFixPayDates,
						double *dFixCoverages,
						long iNFloatPayDates,
						long *lFloatPayDates,
						double *dFloatCoverages,
						double *dFloatSpreads,
						double IRR,
						int UseVol,
						double *dSwapRate);

Err		ComputeSwapRateFromIRR3(char *cYCname, char *cVCname, long lToday,
					    char *cFreq, char *cBasis, char *cRefRname,
						int iNFixPayDates,
						long *lFixPayDates,
						double *dFixCoverages,
						long *lFloatPayDates,
						double *dFloatCoverages,
						double *dFloatSpreads,
						double IRR,
						int UseVol,
						double *dSwapRate);

Err		ConvertAmortSwapWithMarginsInCoInitSwapPortfolio(
								 char *cYCname,
								 char *cVCname,
								 char *cRefRname,
								 SrtCallPutType srtCallPut,
								 long StartDate,
								 long EndDate,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,
								 double exer_fee,
								 long	lNFixNot,
								 double *dFixNotionals,
								 double *dFixRates,
								 long	lNFloatNot,
								 double *dFloatNotionals,
								 double *dMargins,
								 double  *dStrikes,
								 double  *dAmounts,
								 int UseVol);

Err		ConvertAmortSwapWithMarginsInCoInitSwapPortfolio2(
								 char *cYCname,
								 char *cVCname,
								 char *cRefRname,
								 SrtCallPutType srtCallPut,
								 long StartDate,
								 long EndDate,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,
								 double exer_fee,
								 long	lNFixNot,
								 double *dFixNotionals,
								 double *dFixRates,
								 long	lNFloatNot,
								 double *dFloatNotionals,
								 double *dMargins,
								 double  *dStrikes,
								 double  *dAmounts,
								 double  *dLvls,
								 double  *dFwdSwaps,
								 int UseVol);

Err		ConvertAmortSwapWithMarginsInCoInitSwapPortfolioForMAD(
								 char *cYCname,
								 char *cVCname,
								 char *cRefRname,

								 SrtCallPutType srtCallPut,
								 double			exer_fee,

								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,
								 
								 long	lNFixDates,
								 long	*lFixStartDates,
								 long	*lFixEndDates,
								 double *dFixCoverages,
								 double *dFixNotionals,
								 double *dFixRates,

								 long	lNFloatDates,
								 long	*lFloatStartDates,
								 long	*lFloatEndDates,
								 double *dFloatCoverages,
								 double *dFloatNotionals,
								 double *dFloatSpreads,
								 double *dMargins,

								 double  *dStrikes,
								 double  *dAmounts,

								 double  *dLvls,
								 double  *dFwdSwaps,
								 
								 int UseVol);


Err		ConvertAmortSwapWithMarginsInCoInitSwapPortfolioForMAD2(
								 char *cYCname,
								 char *cVCname,
								 char *cRefRname,

								 SrtCallPutType srtCallPut,
								 double			exer_fee,

								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,
								 
								 long	lNFixDates,
								 long	*lFixStartDates,
								 long	*lFixEndDates,
								 double *dFixCoverages,
								 double *dFixNotionals,
								 double *dFixRates,

								 long	lNFloatDates,
								 long	*lFloatStartDates,
								 long	*lFloatEndDates,
								 double *dFloatCoverages,
								 double *dFloatNotionals,
								 double *dFloatSpreads,
								 double *dMargins,

								 double  *dStrikes,
								 double  *dAmounts,

								 double  *dLvls,
								 double  *dFwdSwaps,
								 
								 int UseVol);


#endif