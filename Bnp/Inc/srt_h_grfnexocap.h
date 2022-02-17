/* ===============================================================================

   FILENAME:       srt_h_grfnexocap.h

   PURPOSE:        Pricing of Exotic Caps with Grfn for a given model
                                        -> Builds tableaux and prices of deals
                                        -> AutoCalibrate on Caps for those Deals

   ===============================================================================
 */

#ifndef SRT_H_GRFNEXOCAP_H
#define SRT_H_GRFNEXOCAP_H

Err SrtGRFNQuickExoticCapPrice(
    long StartDate, long EndDate, long lNumExercices, String szRefRateCode,
    double dStrike, double dEps, String szCapFloor, String szYieldCurveName,
    Err (*pfGetVol)(double dStart, double dEnd, double dStrike, double dForward,
                    double dSpread, double *pdBsVol),
    String szVolType, double dTau, String szReturnUnd, double dAlpha,
    double dGamma, double dRho, SRT_Boolean bIsChooser, String szUndName,
    double *pdGrfnFlexiCapPrice, double *pdGrfnFullCapFloorPrice,
    SRT_Boolean *pbReturnUnd);

Err SrtGRFNQuickBarrierCapPrice(
    long StartDate, long EndDate, String szRefRateCode, double dStrike,
    double dBarrier, String szCapFloor, String szUpDown, String szInOut,
    String szYieldCurveName,
    Err (*pfGetVol)(double dStart, double dEnd, double dStrike, double dForward,
                    double dSpread, double *pdBsVol),
    String szVolType, double dTau, String szReturnUnd, double dAlpha,
    double dGamma, double dRho, String szUndName, double *pdGrfnBarrierPrice,
    double *pdGrfnFullCapFloorPrice, double *pdGrfnCapFloorAtBarrier,
    SRT_Boolean *pbReturnUnd);

Err GrfnQuickDisplayTableau(String szProductName, long StartDate, long EndDate,
                            long lNumExercises, String szRefRateCode,
                            double dStrike, double dBarrier, String szCapFloor,
                            String szUpDown, String szInOut,
                            String szYieldCurveName, long *plNumRows,
                            long *plNumCols, GrfnCell ***pppsGrfnTableau);

Err SrtGrfnChooserPrice(
    long *plFixingDates, long *plStartDates, long *plPayDates,
    double *pdCoverages, long lNumDates, String szRefRateCode, double dStrike,
    long lMaxNumExercise, long lNumAlreadyExercised, String szCapFloor,
    String szUndName, String *pszGrfnParamNames, String *pszGrfnParamValues,
    int iNumGrfnParams, double dFullCapRealPrice, double *pdChooserPrice,
    double *pdGrfnFullCapFloorPrice, long *plNumRows, long *plNumCols,
    GrfnCell ***pppsGrfnTableau, double ***pppdLastPath, long *plNumAuxColumns,
    long **pplAuxRangesLength, double ***pppdAuxRanges);
#endif
