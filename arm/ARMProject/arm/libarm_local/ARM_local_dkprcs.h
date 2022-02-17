/*---------------------------------------------------------------------------------*/

#ifndef _ARM_LOCAL_DKPRCS
#define _ARM_LOCAL_DKPRCS







extern long ARMLOCAL_PRCS3F_Lattice_HWVFDK_Pricing(VECTOR<double>& dLatticeGeometryDataIn,
		 										   double dNumTimeLinesBeforeFirstNotice,
												   double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
												   double dNumTimeLines, 
												   double evalDate, 
												   long dBaseYieldCurveId,
												   long dForeignYieldCurveId,
												   long dBaseRatesNoBasisCurveId,
											       long dForeignRatesNoBasisCurveId,
												   long volSwopBaseId,
												   long volSwopForeignId,
												   long volFxId,
												   VECTOR<double>& dNoticeDatesIn,
												   double dStrike,
												   double dType,
												   double dOptionExpiry, 
												   double dMeanReversionBase,
												   double dMeanReversionForeign,  
												   double dSpotFX,
												   long dBaseForeignCorrelationId,
												   long dBaseSpotFXCorrelationId,
												   long dForeignSpotFXCorrelationId,
												   double dProductModelCodeId,
												   double dX1Limit,
												   double dX2Limit,
												   double dX3Limit,
												   double dI1Limit,
												   double dI2Limit,
												   double dI3Limit,
												   double dADPLimit,
												   double dOptimal,
												   double dTimeBoost,
												   double dDeltaFlag,
												   double dStringModel,
												   long nbCol,
												   VECTOR<double>& dBoosterDataIn,
												   ARM_result& result);


extern long ARMLOCAL_TREE3FACT (double asof,
								long xbxfxId,
								long volSwopBaseId,
								long volSwopForeignId,
								long dBaseForeignCorrelationId,
								VECTOR<double>& dLatticeGeometryDataIn,
								double dNumTimeLinesBeforeFirstNotice,
								double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
								double dNumTimeLines,
								double dMeanReversionBase,
								double dMeanReversionForeign,
								double dX1Limit,
								double dX2Limit,
								double dX3Limit,
								double dI1Limit,
								double dI2Limit,
								double dI3Limit,
								double dADPLimit,
								double dOptimal,
								double dTimeBoost,
								double dDeltaFlag,
								double dSmoothingValue,
								long calcProbaSurvOrNotId,
								double QBaseSmile,
				                double QForeignSmile,
							    long spotVolOrFwdVol,
								double histoVolLongTerm,
								long convInVolFwd,
								long calibBasisIncluded,
								ARM_result& result,
								long objId = -1);

extern long  ARMLOCAL_Bootstrapping_VFDK_HW1To3F (long volId,
								                  long zcId,
								                  VECTOR<double>& noticeDates,
								                  VECTOR<double>& swapStartDates,
                                                  VECTOR<double>& swapEndDates,
                                                  VECTOR<double>& HW3FParameters,
                                                  double observationDate,
                                                  ARM_result& result);

extern long  ARMLOCAL_SwaptionPrice_VFDK_HW1To3F ( double dSwaptionExpiryInYears,
                                            double dSwaptionTenorInYears,
                                            double dNoticePeriodInDays,
                                            double dStrike,
                                            double dCallPut,
                                            long zcId,
                                            VECTOR<double>& dNoticeDates,
                                            VECTOR<double>& dSigma,
                                            VECTOR<double>& HW3FParameters,
                                            double observationDate,
                                            ARM_result& result);

long  ARMLOCAL_ImpliedFwdCorrelation_VFDK_HW1To3F ( double dSwaptionExpiryInYears,
                                                    double dSwaptionTenorInYears,
                                                    double dSwaptionTenor2InYears,
                                                    double dNoticePeriodInDays,
                                                    long zcId,
                                                    VECTOR<double>& dNoticeDates,
                                                    VECTOR<double>& dSigma,
                                                    VECTOR<double>& HW3FParameters,
                                                    double observationDate,
                                                    ARM_result& result);


long  ARMLOCAL_HW3F_CALIBRATION (   double dSpotDate,
                                    long zcId,
                                    VECTOR<double>& HW3FParametersIn,
                                    long volId,
                                    long correlationId,
                                    long volWeightsId,
                                    long correlationWeightsId,
                                    VECTOR<double>& dNoticeDates,
                                    VECTOR<double>& dSwapStartDates,
                                    VECTOR<double>& dSwapEndDates,
                                    ARM_result& result);




#endif
/*---------------------------------------------------------------------------------*/
/*---- End Of File ----*/