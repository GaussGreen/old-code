/*
 * $Log: hw_vfdk_analytics.h,v $
 * Revision 1.11  2004/06/28 18:02:56  spannetier
 * Integration du Q Model
 *
 * Revision 1.10  2004/05/10 12:45:33  mcampet
 * MC from SP several changes
 *
 * Revision 1.9  2004/04/22 12:02:56  mcampet
 * MC for SP GetPriceFromAbsVol
 * 
 *
 */


#include "DKMaille.h"
#include "DKMaille2D.h"

double BasisConversion(double dSpotDate,
											 DKMaille<double> &dFundingDiscountCurveDates,
											 DKMaille<double> &dFundingDiscountCurve,
											 DKMaille<double> &dFundingDiscountCurveAdjusted,
											 DKMaille<double> &dDomesticDiscountCurveDates,
											 DKMaille<double> &dDomesticDiscountCurve,
											 DKMaille<double> &dDomesticDiscountCurveAdjusted,
											 DKMaille<double> &dFundingStartDates,
											 DKMaille<double> &dFundingEndDates,
											 DKMaille<double> &dDomesticStartDates,
											 DKMaille<double> &dDomesticEndDates,
											 double dActualisationFunding,
											 double dActualisationDomestic,
											 double dFundingNotional,
											 double dDomesticNotional,
											 double dFundingSpread,
											 double dFXSpotRate);

double Analytical_QModel(double dSpotDate,
                         DKMaille<double> dDiscountDates,
                         DKMaille<double> dDiscountRates,
                         DKMaille<double> dSwaptionDates,
                         double dNoticePeriod,
												 double dStrike,
												 double dLognormalNormalSmileModelType,
                         double dQVol,
												 double dQParameter,
												 double dCallPut);  

void TransformVolatilitiesSingle( DKMaille<double> dStripDates,
                 DKMaille<double> &dNewStripDates,
                 DKMaille<double> &dVolStrip);

double Delta_QVol(double dSpotDate,
								  DKMaille<double> dDiscountDates,
								  DKMaille<double> dDiscountRates,
								  DKMaille<double> dSwaptionDates,
								  double dNoticePeriod,
								  double dLognormalNormalSmileModelType,
								  double dQVol,
								  double dQParameter,
									double dShift); 
 
double  Analytics_DK_3F_(double dSpotDate, 
												 DKMaille<double> dBaseDates,
												 DKMaille<double> dBaseRates,
												 DKMaille<double> dForeignDates,
												 DKMaille<double> dForeignRates,
												 DKMaille<double> dStdDevBaseX,
												 DKMaille<double> dStdDevBaseZ,
												 double dMeanReversionBase,  
												 DKMaille<double> dStdDevForeignX,
												 DKMaille<double> dStdDevForeignZ,
												 double dMeanReversionForeign,  
												 double dSpotFX,  
												 DKMaille<double> dSpotFXVolDatesTD,
												 DKMaille<double> dSpotFXVolTD,
												 double dBaseForeignCorrelation,
												 double dBaseSpotFXCorrelation,
												 double dForeignSpotFXCorrelation,
												 DKMaille2D<double> dBoosterData,
												 DKMaille<double> dRedemptionData);

double  Analytics_HWVFDK_3F_(double dSpotDate, 
														 DKMaille<double> dBaseDates,
														 DKMaille<double> dBaseRates,
														 DKMaille<double> dForeignDates,
														 DKMaille<double> dForeignRates,
														 DKMaille<double> dStdDevBaseX,
														 DKMaille<double> dStdDevBaseZ,
														 double dMeanReversionBase,  
														 DKMaille<double> dStdDevForeignX,
														 DKMaille<double> dStdDevForeignZ,
														 double dMeanReversionForeign,  
														 double dSpotFX,  
														 DKMaille<double> dSpotFXVolDatesTD,
														 DKMaille<double> dSpotFXVolTD,
														 double dBaseForeignCorrelation,
														 double dBaseSpotFXCorrelation,
														 double dForeignSpotFXCorrelation,
														 DKMaille2D<double> dBoosterData,
														 DKMaille<double> dRedemptionData);

double Analytical_VFDK_HW1F(double dSpotDate,
                            DKMaille<double> dDiscountDates,
                            DKMaille<double> dDiscountRates,
                            DKMaille<double> dSwaptionDates,
                            double dNoticePeriod,
                            double dMeanReversion,
                            double dBlackScholesVol,
							double dOutput);  



DKMaille<double> Bootstrapping_VFDK_HW1To3F(DKMaille2D<double> dAbsVol,
																						DKMaille<double> dSwaptionExpiry,
																						DKMaille<double> dSwaptionTenor,
																						DKMaille<double> dDiscountDates,
																						DKMaille<double> dDiscountRates,
																						DKMaille<double> dAdjustedDiscountRates,
																						DKMaille<double> dNoticeDates,
																						DKMaille<double> dSwapStartDates,
																						DKMaille<double> dSwapEndDates,
																						DKMaille<double> dModelParameters,
																						double dJulianObservationDate);

DKMaille<double> BootstrappingSpotFXVolatility3FTD(double dSpotDate,
																									 DKMaille<double> dStartDate,
																									 DKMaille<double> dEndDate,
																									 DKMaille<double> dForwardMaturityDate,
																									 DKMaille<double> dForwardVol,
																									 DKMaille<double> dDomesticStrip,
																									 DKMaille<double> dForeignStrip,
																									 DKMaille<double> dDomesticStdDev,
																									 DKMaille<double> dForeignStdDev,
																									 double dDomesticMeanReversion,
																									 double dForeignMeanReversion,
																									 double dCorrelationSpotFXDomestic,
																									 double dCorrelationSpotFXForeign,
																									 double dCorrelationForeignDomestic);

DKMaille<double> BootstrappingSpotFXVolatility3FTD_DK(double dSpotDate,
																											DKMaille<double> dStartDate,
																											DKMaille<double> dEndDate,
																											DKMaille<double> dForwardMaturityDate,
																											DKMaille<double> dForwardVol,
																											DKMaille<double> dDomesticStrip,
																											DKMaille<double> dForeignStrip,
																											DKMaille<double> dDomesticStdDev,
																											DKMaille<double> dForeignStdDev,
																											double dDomesticMeanReversion,
																											double dForeignMeanReversion,
																											double dCorrelationSpotFXDomestic,
																											double dCorrelationSpotFXForeign,
																											double dCorrelationForeignDomestic,
																											DKMaille<double> dCurveDates,
																											DKMaille<double> dCurve,
																											DKMaille<double> dCurveForeignDates,
																											DKMaille<double> dCurveForeign);




double SwaptionPrice_VFDK_HW1To3F(double dSwaptionExpiryInYears,
																	double dSwaptionTenorInYears,
																	double dNoticePeriodInDays,
																	double dStrike,
																	double dCallPut,
																	DKMaille<double> dDiscountCurveDates,
																	DKMaille<double> dDiscountCurveRates,
																	DKMaille<double> dAdjustedDiscountCurveRates,
																	DKMaille<double> dNoticeDates,
																	DKMaille<double> dTimeDependentStandardDeviation,
																	DKMaille<double> dModelParameters,
																	double dJulianObservationdate); 

double ImpliedFwdCorrelation_VFDK_HW1To3F(double dSwaptionExpiryInYears,
																					double dSwaptionTenorInYears,
																					double dSwaptionTenor2InYears,
																					double dNoticePeriodInDays,
																					DKMaille<double> dDiscountCurveDates,
																					DKMaille<double> dDiscountCurveRates,
																					DKMaille<double> dAdjustedDiscountCurveRates,
																					DKMaille<double> dNoticeDates,
																					DKMaille<double> dTimeDependentStandardDeviation,
																					DKMaille<double> dModelParameters,																	
																					double dJulianObservationdate); 

double numerical_derivative(double dDate,
														DKMaille<double> dDates,
														DKMaille<double> dFunction);

double numerical_function(double dDate,
													DKMaille<double> dDates,
													DKMaille<double> dFunction);

double  ForwardFXVolatility3FTD(double dSpotDate,
																double dStartDate,
																double dEndDate,
																double dForwardMaturityDate,
																DKMaille<double> dDomesticStrip,
																DKMaille<double> dForeignStrip,
																DKMaille<double> dSpotFXStrip,
																DKMaille<double> dDomesticStdDev,
																DKMaille<double> dForeignStdDev,
																DKMaille<double> dSpotFXVol,
																double dDomesticMeanReversion,
																double dForeignMeanReversion,
																double dCorrelationSpotFXDomestic,
												 				double dCorrelationSpotFXForeign,
																double dCorrelationForeignDomestic);

void CalcShortRates(double dSpotDate,
										DKMaille<double> &dDates,
										DKMaille<double> &dRates,
										DKMaille<double> &dAdjRates,
										DKMaille<double> &dShortRateDates,
										DKMaille<double> &dShortRates);

void CalcXCCYShortRates(double dSpotDate,
												DKMaille<double> &dDates,
												DKMaille<double> &dRates,
												DKMaille<double> &dAdjRates,
												DKMaille<double> &dForeignDates,
												DKMaille<double> &dForeignRates,
												DKMaille<double> &dForeignAdjRates,
												DKMaille<double> &dShortRateDates,
												DKMaille<double> &dShortRates,
												DKMaille<double> &dShortRatesForeign);

double ForwardFXVolatility3FTD_DK(double dSpotDate,
															    double dStartDate,
																	double dEndDate,
																	double dForwardMaturityDate,
																	DKMaille<double> dDomesticStrip,
																	DKMaille<double> dForeignStrip,
																	DKMaille<double> dSpotFXStrip,
																	DKMaille<double> dDomesticStdDev,
																	DKMaille<double> dForeignStdDev,
																	DKMaille<double> dSpotFXVol,
																	double dDomesticMeanReversion,
																	double dForeignMeanReversion,
																	double dCorrelationSpotFXDomestic,
												 					double dCorrelationSpotFXForeign,
																	double dCorrelationForeignDomestic,
																	DKMaille<double> shortRateDates,
																	DKMaille<double> shortRates,
																	DKMaille<double> shortRatesForeign);

void dIntermediateCalcSwap(double dSpotDate,
													 const DKMaille<double> &dDiscountDates,
													 const DKMaille<double> &dDiscountRates,
													 const DKMaille<double> &dAdjustedDiscountRates,
													 const DKMaille<double> &dSwaptionDates,
													 const DKMaille<double> &dBasisSwaptionDates,
                           const DKMaille<double> &dAccrualPeriods,
                           const DKMaille<double> &dBasis,
                           const DKMaille<double> &dBasisAccrualPeriods,
													 DKMaille<double> &dDiscountFactors,
													 DKMaille<double> &dDiscountFactorsOnBasisDates,
													 double *dSwapRate,
													 double *dAnnuity,
													 double *dBasisAnnuity,
													 double dJulianObservationDate);

double GetSwapRate(double dSpotDate,
                   const DKMaille<double> &dSwaptionDates,
                   const DKMaille<double> &dAccrualPeriods,
									 const DKMaille<double> &dDiscountDates,
									 const DKMaille<double> &dDiscountRates,
									 const DKMaille<double> &dAdjustedDiscountRates);




void TransformVolatilities(DKMaille<double> dDomesticVolStripDates,
													 DKMaille<double> dForeignVolStripDates,
													 DKMaille<double> dSpotFXStrip,
													 DKMaille<double> &dNewStripDates,
													 DKMaille<double> &dDomesticVolStrip,
													 DKMaille<double> &dForeignVolStrip,
													 DKMaille<double> &dSpotFXVol);


void TransformFXVolatilities(DKMaille<double> dSpotFXStrip,
													   DKMaille<double> &dNewStripDates,
													   DKMaille<double> &dSpotFXVol);

void TransformSwaptionVolatilities(DKMaille<double> dDomesticVolStripDates,DKMaille<double> dForeignVolStripDates,
																	 DKMaille<double> dSpotFXStrip,DKMaille<double> &dNewStripDates,DKMaille<double> &dDomesticVolStrip,
																	 DKMaille<double> &dForeignVolStrip);




DKMaille<double>  &FromSpotToForwardVol3F( double dSpotDate,
                      DKMaille<double> dStartDate,
                      DKMaille<double> dEndDate,
                      DKMaille<double> dForwardMaturityDate,
                      DKMaille<double> dBaseDiscountCurveDates,
                      DKMaille<double> dBaseDiscountCurveRates,
                      DKMaille<double> dBaseAdjustedDiscountCurveRates,
                      DKMaille<double> dForeignDiscountCurveDates,
                      DKMaille<double> dForeignDiscountCurveRates,
                      DKMaille<double> dForeignAdjustedDiscountCurveRates,
                      DKMaille<double> dDomesticVolX,
                      DKMaille<double> dDomesticVolY,
                      DKMaille2D<double> dDomesticVolZ,
                      DKMaille<double> dForeignVolX,
                      DKMaille<double> dForeignVolY,
                      DKMaille2D<double> dForeignVolZ,
                      DKMaille<double> dSpotFXStrip,
                      DKMaille<double> dSpotFXVol,
                      double dDomesticMeanReversion,
                      double dForeignMeanReversion,
                      double dCorrelationSpotFXDomestic,
                      double dCorrelationSpotFXForeign,
                      double dCorrelationForeignDomestic,
                      double dCurrencyPair1, // 0=jpy 1=usd 2=aud 3=eur
                      double dCurrencyPair2, // 0=jpy 1=usd 2=aud 3=eur
                      DKMaille<double> dNoticeDates,
                      DKMaille<double> dFXCouponResetDates,
                      DKMaille<double> dFXCouponPaymentDates,
											double dIsSwaptionCalibrationWithBasis);

DKMaille<double>  &FromSpotToForwardVolDK3F(double dSpotDate,
																						DKMaille<double> dStartDate,
																						DKMaille<double> dEndDate,
																						DKMaille<double> dForwardMaturityDate,
																						DKMaille<double> dBaseDiscountCurveDates,
																						DKMaille<double> dBaseDiscountCurveRates,
																						DKMaille<double> dBaseAdjustedDiscountCurveRates,
																						DKMaille<double> dForeignDiscountCurveDates,
																						DKMaille<double> dForeignDiscountCurveRates,
																						DKMaille<double> dForeignAdjustedDiscountCurveRates,
																						DKMaille<double> dDomesticVolX,
																						DKMaille<double> dDomesticVolY,
																						DKMaille2D<double> dDomesticVolZ,
																						DKMaille<double> dForeignVolX,
																						DKMaille<double> dForeignVolY,
																						DKMaille2D<double> dForeignVolZ,
																						DKMaille<double> dSpotFXStrip,
																						DKMaille<double> dSpotFXVol,
																						double dDomesticMeanReversion,
																						double dForeignMeanReversion,
																						double dCorrelationSpotFXDomestic,
																						double dCorrelationSpotFXForeign,
																						double dCorrelationForeignDomestic,
																						double dCurrencyPair1, // 0=jpy 1=usd 2=aud 3=eur
																						double dCurrencyPair2, // 0=jpy 1=usd 2=aud 3=eur
																						DKMaille<double> dNoticeDates,
																						DKMaille<double> dFXCouponResetDates,
																						DKMaille<double> dFXCouponPaymentDates,
																						double dIsSwaptionCalibrationWithBasis);  

double ForwardFXVolatility3FConstant(double dTimeToReset,
																		 double dDomesticStdDev,
																		 double dForeignStdDev,
																		 double dSpotFXVol,
																		 double dDomesticMeanReversion,
																		 double dForeignMeanReversion,
																		 double dCorrelationSpotFXDomestic,
																		 double dCorrelationSpotFXForeign,
																		 double dCorrelationForeignDomestic);


double BondFromVFDKAnalytics_HW1F(double dNorm,
                                  double dDecay,
                                  double dBrownian);

double BondFromVFDKAnalytics_HW2F(double dNorm,
                                  double dDecay1,
                                  double dBrownian1,
                                  double dDecay2,
                                  double dBrownian2);

double BondFromVFDKAnalytics_HW3F(double dNorm,
                                  double dDecay1,
                                  double dBrownian1,
                                  double dDecay2,
                                  double dBrownian2,
                                  double dDecay3,
                                  double dBrownian3);

double ki_cap_analytics_price(double dSpotDate,
                              double dDateOne,
                              double dDateTwo,
                              double dStdDev1,
                              double dStdDev2,
                              double dRate1,
                              double dRate2,
                              double dStrike1,
                              double dStrike2,
                              double dCorrelation,
                              double dOutput);

double GetOptionPrice(double dCallPut,
                      double dSwapRate,
                      double dAnnuity,
                      double dStrike,
                      double dForwardStdDev,
                      double dTimeToExpiry);

void Processing_3F_Data_DK(DKMaille<double> &dStripDatesVolFX,
                         DKMaille<double> &dStripDatesVolBase,
                         DKMaille<double> &dStripDatesVolForeign,
													 DKMaille<double> &dVolFX,
													 DKMaille<double> &dVolBase,
													 DKMaille<double> &dVolForeign,
													 double dNumTimeLinesBeforeFirstNotice,
													 double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
													 double dNumTimeLinesPerYear, 
													 double dSpotDate, 
													 DKMaille<double> &dBaseDates,
													 DKMaille<double> &dBaseRates,
													 DKMaille<double> &dForeignDates,
													 DKMaille<double> &dForeignRates,
													 DKMaille<double> &dStdDevBaseX,
													 DKMaille<double> &dStdDevBaseY,
													 DKMaille2D<double> &dStdDevBaseZ,
													 double dMeanReversionBase,  
													 DKMaille<double> &dStdDevForeignX,
													 DKMaille<double> &dStdDevForeignY,
													 DKMaille2D<double> &dStdDevForeignZ,
													 double dMeanReversionForeign,  
													 DKMaille<double> dSpotFXVolDatesTD,
													 DKMaille<double> dSpotFXVolTD,
													 double dBaseForeignCorrelation,
													 double dBaseSpotFXCorrelation,
													 double dForeignSpotFXCorrelation,           
													 double dOptimal,
													 DKMaille2D<double> &dBoosterData,
													 DKMaille<double> &dBaseRatesNoBasis,
													 DKMaille<double> &dForeignRatesNoBasis,
													 double dIsSwaptionCalibrationWithBasis);


void Processing_3F_Data( DKMaille<double> &dStripDatesVolFX,
                         DKMaille<double> &dStripDatesVolBase,
                         DKMaille<double> &dStripDatesVolForeign,
                         DKMaille<double> &dVolFX,
                         DKMaille<double> &dVolBase,
                         DKMaille<double> &dVolForeign,
                         double dNumTimeLinesBeforeFirstNotice,
                         double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                         double dNumTimeLinesPerYear,
                         double dSpotDate,
                         DKMaille<double> &dBaseDates,
                         DKMaille<double> &dBaseRates,
                         DKMaille<double> &dForeignDates,
                         DKMaille<double> &dForeignRates,
                         DKMaille<double> &dStdDevBaseX,
                         DKMaille<double> &dStdDevBaseY,
                         DKMaille2D<double> &dStdDevBaseZ,
                         double dMeanReversionBase,
                         DKMaille<double> &dStdDevForeignX,
                         DKMaille<double> &dStdDevForeignY,
                         DKMaille2D<double> &dStdDevForeignZ,
                         double dMeanReversionForeign,
                         DKMaille<double> dSpotFXVolDatesTD,
                         DKMaille<double> dSpotFXVolTD,
                         double dBaseForeignCorrelation,
                         double dBaseSpotFXCorrelation,
                         double dForeignSpotFXCorrelation,
                         double dOptimal,
                         DKMaille2D<double> &dBoosterData,
                         DKMaille<double> &dBaseRatesNoBasis,
                         DKMaille<double> &dForeignRatesNoBasis,
                         double dIsSwaptionCalibrationWithBasis);



void SetSpotFXVolAnalytics( DKMaille<double> &dFXVolDates,
              DKMaille<double> &dFXVol,
              double dSpotDate,
              DKMaille<double> &dStdDevBaseX,
              DKMaille<double> &dStdDevForeignX,
              DKMaille<double> &dStdDevBaseZ,
              DKMaille<double> &dStdDevForeignZ,
              double dMeanReversionBase, 
              double dMeanReversionForeign,  
              double dBaseForeignCorrelation,
              double dBaseSpotFXCorrelation,
              double dForeignSpotFXCorrelation,
              DKMaille2D<double> &dBoosterData); 

void SetSpotFXVolAnalytics_DK(DKMaille<double> &dFXVolDates,
															DKMaille<double> &dFXVol,
															double dSpotDate,
															DKMaille<double> &dStdDevBaseX,
															DKMaille<double> &dStdDevForeignX,
															DKMaille<double> &dStdDevBaseZ,
															DKMaille<double> &dStdDevForeignZ,
															double dMeanReversionBase, 
															double dMeanReversionForeign,  
															double dBaseForeignCorrelation,
															double dBaseSpotFXCorrelation,
															double dForeignSpotFXCorrelation,
															DKMaille2D<double> &dBoosterData,
															DKMaille<double> &dCurveDates,
															DKMaille<double> &dCurve,
															DKMaille<double> &dCurveForeignDates,
															DKMaille<double> &dCurveForeign);
				 
DKMaille<double> Create_Strip_Analytics( DKMaille<double> NoticeDates, 
                     int dNumTimeLinesBeforeFirstNotice,
                     int dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                     int dNumTimeLinesPerYear,
                     double dOptimalDate, 
                     double LastDate); 
 
void BootstrappedVolsSingleCurrencyFixedIncome1F( DKMaille<double> &dNewVolStripDates,
                         DKMaille<double> &dOutputStrip,
                         DKMaille<double> &dNoticeDatesBis,
                         DKMaille<double> &dDates,
                         DKMaille<double> &dRates,
                         DKMaille<double> &dRatesNoBasis,
                         DKMaille<double> &dStdDevX,
                         DKMaille<double> &dStdDevY,
                          DKMaille2D<double> &dStdDevZ,
                         double dSpotDate,
                         DKMaille2D<double> &dBoosterData,
                         double dMeanReversion);

 void BootstrappedVolsSingleCurrencyFixedIncome1F_DK(DKMaille<double> &dNewVolStripDates,
													DKMaille<double> &dOutputStrip,
													DKMaille<double> &dNoticeDatesBis,
													DKMaille<double> &dDates,
													DKMaille<double> &dRates,
													DKMaille<double> &dRatesNoBasis,
													DKMaille<double> &dStdDevX,
													DKMaille<double> &dStdDevY,
													DKMaille2D<double> &dStdDevZ,
													double dSpotDate,
													DKMaille2D<double> &dBoosterData,
													double dMeanReversion);

double SwaptionAbsVol_VFDK_HW1To3F(double dSwaptionExpiryInYears,
								 double dSwaptionTenorInYears,
								 double dNoticePeriodInDays,
								 DKMaille<double> dDiscountCurveDates,
								 DKMaille<double> dDiscountCurveRates,
								 DKMaille<double> dAdjustedDiscountCurveRates,
								 DKMaille<double> dNoticeDates,
								 DKMaille<double> dTimeDependentStandardDeviation,
								 DKMaille<double> dModelParameters,
								 double dJulianObservationDate);

double GetAdjSwapRate(double dSpotDate,
                      const DKMaille<double> &dSwaptionDates,
                      const DKMaille<double> &dAccrualPeriods,
                      const DKMaille<double> &dBasisSwaptionDates,
                      const DKMaille<double> &dBasis,
                      const DKMaille<double> &dBasisAccrualPeriods,
											const DKMaille<double> &dDiscountDates,
											const DKMaille<double> &dDiscountRates,
											const DKMaille<double> &dAdjustedDiscountRates);

double GetPriceFromAbsVol(double dSpotDate,
													double dJulianObservationDate,
                          const DKMaille<double> &dSwaptionDates,
                          const DKMaille<double> &dAccrualPeriods,
                          double dNoticePeriod,
                          double dAbsoluteVol,
													const DKMaille<double> &dDiscountDates,
													const DKMaille<double> &dDiscountRates,
													const DKMaille<double> &dAdjustedDiscountRates);

double GetAdjAnnuity(double dSpotDate,
                     const DKMaille<double> &dSwaptionDates,
                     const DKMaille<double> &dAccrualPeriods,
									   const DKMaille<double> &dDiscountDates,
									   const DKMaille<double> &dDiscountRates,
									   const DKMaille<double> &dAdjustedDiscountRates,
										 double dJulianObservationDate);

void FromForwardToSpotVol3F(double dSpotDate,
                            DKMaille<double> dBaseDiscountCurveDates,
                            DKMaille<double> dBaseDiscountCurveRates,
                            DKMaille<double> dBaseAdjustedDiscountCurveRates,
                            DKMaille<double> dForeignDiscountCurveDates,
                            DKMaille<double> dForeignDiscountCurveRates,
                            DKMaille<double> dForeignAdjustedDiscountCurveRates,
                            DKMaille<double> dDomesticVolX,
                            DKMaille<double> dDomesticVolY,
                            DKMaille2D<double> dDomesticVolZ,
                            DKMaille<double> dForeignVolX,
                            DKMaille<double> dForeignVolY,
                            DKMaille2D<double> dForeignVolZ,
                            DKMaille<double> dForwardFXStrip,
                            DKMaille<double> dForwardFXVol,
                            double dDomesticMeanReversion,
                            double dForeignMeanReversion,
                            double dCorrelationSpotFXDomestic,
                            double dCorrelationSpotFXForeign,
                            double dCorrelationForeignDomestic,
                            double dCurrencyPair1, // 0=jpy 1=usd 2=aud 3=eur
                            double dCurrencyPair2, // 0=jpy 1=usd 2=aud 3=eur
                            DKMaille<double> dNoticeDates,
                            DKMaille<double> dFXCouponResetDates,
                            DKMaille<double> dFXCouponPaymentDates,
                            double dIsSwaptionCalibrationWithBasis, 
                            DKMaille<double> &dSpotFXStrip,
                            DKMaille<double> &dSpotFXVol);
 
