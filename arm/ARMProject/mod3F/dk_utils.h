
#include "DKMaille.h"
#include "DKMaille2D.h"

int normalDK(double *res,double xx);

void Sort(DKMaille<double> &d,
          int n);

double dBlackScholes(double dTimeToExpiry,
										 double dVol,
										 double dStrike,
										 double dForward,
										 double dCallPut);

void AggregateVolStrip(DKMaille<double> dStrip,DKMaille<double> dInputStrip,DKMaille<double> dVolStrip,DKMaille<double> &dStripOutput);
double LinearInterpolation(double target_point,const DKMaille<double> &x,const DKMaille<double> &y);

DKMaille<double> Reverse(DKMaille<double> &Input);
double isConstant(DKMaille2D<double> &dAbsVol);
double dItoIntegral_HW1F_TD(double dMeanFieldDecay, 
                            double dTerminalDate, 
                            const DKMaille<double> &dDatesStrip, 
                            const DKMaille<double> &dSigmaStrip);

double SpotFXVolIntegral(DKMaille<double> &dStrip,
												 DKMaille<double> &dSpotFXVol,
												 double dStartDate,
												 double dEndDate);

double SpotFXVolIRStdDevQuantoIntegral(DKMaille<double> &dStrip,
																			 DKMaille<double> &dSpotFXVol,
																			 DKMaille<double> &dStripIRStdDev,
																			 double dMeanReversion,
																			 double dStartDate,
																			 double dEndDate);

double CumulativeNormalBivariateIntegral(double da, double db, double drho);

double GetForwardFXVol(double dSpotDate,
											 double dStartDate,
											 double dEndDate,
											 double dForwardMaturityDate,
											 DKMaille<double> dStrip,
											 DKMaille<double> &dStripDomesticStdDev,
											 DKMaille<double> &dStripForeignStdDev,
											 DKMaille<double> &dStripSpotFXVol,
											 double dDomesticMeanReversion,
											 double dForeignMeanReversion,
											 double dCorrelationSpotFXDomestic,
											 double dCorrelationSpotFXForeign,
											 double dCorrelationForeignDomestic);

void AppendDates(const DKMaille<double> &dNoticeDates,
								 const DKMaille<double> &dSwaptionDates,
								 const DKMaille<double> &dBasisSwaptionDates,
								 DKMaille<double> &dSliceDates);

void AppendFXDates(const DKMaille<double> &dNoticeDates,
								   double dOptionMaturityDate,
								   DKMaille<double> &dSliceDates);

void GetSlices(const DKMaille<double> &shortRateDates,
							 const DKMaille<double> &shortRates,
							 const DKMaille<double> &dInputDates,
							 DKMaille<double> &dSlices,
							 DKMaille<double> &dShortRateOnSlice,
							 double dt);

void GetShortRateOnSlice(const DKMaille<double> &shortRateDates,
												 const DKMaille<double> &shortRates,
												 const DKMaille<double> &dSlices,
												 DKMaille<double> &dShortRateOnSlices,
												 DKMaille<double> &dShortRateKOnSlices);

void GetXCCYSlices(const DKMaille<double> &shortRateDates,
									 const DKMaille<double> &shortRates,
									 const DKMaille<double> &shortRatesForeign,
									 const DKMaille<double> &dInputDates,
									 DKMaille<double> &dSlices,
									 DKMaille<double> &dShortRateOnSlice,
									 DKMaille<double> &dShortRateForeignOnSlice,
									 double dt);

double GetForwardFXVol_DK(double dSpotDate,
											    double dStartDate,
											    double dEndDate,
											    double dForwardMaturityDate,
											    DKMaille<double> dStrip,
											    DKMaille<double> &dStripDomesticStdDev,
											    DKMaille<double> &dStripForeignStdDev,
											    DKMaille<double> &dStripSpotFXVol,
											    double dDomesticMeanReversion,
											    double dForeignMeanReversion,
											    double dCorrelationSpotFXDomestic,
											    double dCorrelationSpotFXForeign,
											    double dCorrelationForeignDomestic,
											    const DKMaille<double> &shortRateDates,
													const DKMaille<double> &shortRates,
													const DKMaille<double> &shortRatesForeign);

double IRStdDevIntegral(DKMaille<double> &dStrip,
												DKMaille<double> &dStripIR1StdDev,
												DKMaille<double> &dStripIR2StdDev,
												double dMeanReversion1,
												double dMeanReversion2,
												double dStartDate,
												double dEndDate,
												double dForwardMaturityDate);

double IRStdDevIntegral_DK(DKMaille<double> &dStrip,
												   DKMaille<double> &dStripIR1StdDev,
												   DKMaille<double> &dStripIR2StdDev,
												   double dMeanReversion1,
												   double dMeanReversion2,
												   double dStartDate,
												   double dEndDate,
												   double dForwardMaturityDate,
													 const DKMaille<double> &shortRateDates,
													 const DKMaille<double> &shortRates,
													 const DKMaille<double> &shortRatesForeign);

double dKacDeterminant2F_TD_DK(double dMeanFieldDecay, 
															 double dMeanFieldDecay2, 
                               double dObservationDate, 
															 double dTerminalDate, 
															 double dDate1, 
                               double dDate2,
															 const DKMaille<double> &dSlices,
															 const DKMaille<double> &dShortRateOnSlice,
															 const DKMaille<double> &dShortRateOnSlice2);

DKMaille<double> &dCheckDates(DKMaille<double> &dStrip,
														 double dStartDate,
														 double dEndDate);

DKMaille<double> &dNewVols(DKMaille<double> &dStrip,
													 DKMaille<double> &dVols,
													 DKMaille<double> &dNewStrip);

double dItoCrossIntegral_HW1F_TD(double dMeanFieldDecay, 
                            double dTerminalDate, 
                            const DKMaille<double> &dDatesStrip, 
                            const DKMaille<double> &dSigmaStrip,
                            const DKMaille<double> &dFXSigmaStrip);

double dKacDeterminant2F_TD(double dMeanFieldDecay, double dMeanFieldDecay2, 
                            double dObservationDate, double dTerminalDate, 
                            double dSigma1, double dSigma2, double dDate1, 
                            double dDate2);

double dItoIntegral_HW2F_TD(double dMeanFieldDecay, 
                            double dMeanFieldDecay2,
                            double dObservationDate, 
                            double dTerminalDate, 
                            const DKMaille<double> &dDatesStrip, 
                            const DKMaille<double> &dSigmaStrip, 
                            const DKMaille<double> &dSigmaStrip2);

DKMaille<double> CreateDateStrip(double dSpotDate,
																 DKMaille<double> &dDomesticStrip,
																 DKMaille<double> &dForeignStrip,
																 DKMaille<double> &dSpotFXStrip);

double interpolation2(int method, 
				              double target_point, 
							        const DKMaille<double> &time_array, 
											const DKMaille<double> &rate_array);

double rateinterpolation_dk_maille(int method, 
					                         double target_point, 
										               const DKMaille<double> &time_array, 
                                   const DKMaille<double> &rate_array, 
                                   int array_length);

DKMaille<double> dKacVector(unsigned int uiSize, 
														double dSwapRate, 
                            DKMaille<double> dPeriods, 
														DKMaille<double> dDiscountFactors);

double dKacDeterminant(double dMeanFieldDecay, 
											 double dObservationDate, 
											 double dTerminalDate, 
                       double dSigmaDate1, 
											 double dSigmaDate2, 
											 double dDate1, 
											 double dDate2);

double ATMpricevol_DK(double tte, double inputvol);

double InterpolateMatrix(const DKMaille2D<double> &matrix,
			                   double dX, 
							           double dY,
										     DKMaille<double>  &xAxis,
                         DKMaille<double>  &yAxis);

void CreateVanillaSwap(double dSpotDate,
											 double dXX, 
											 double dYY,
                       DKMaille<double> &dSwaptionDates,
											 DKMaille<double> &dAccrualPeriods,
                       DKMaille<double> &dBasisSwaptionDates,
											 DKMaille<double> &dBasisAccrualPeriods,
                       DKMaille<double> &dBasis,
                       double dNoticePeriod,
											 DKMaille<double> &dCurveDates,
											 DKMaille<double> &dCurve,
                       DKMaille<double> &dCurveAdjusted);

void LatticeJacobiTransformation(DKMaille2D<double> &a,
                                 int n, 
                                 DKMaille<double> &d,
                                 DKMaille2D<double> &v,
                                 int *nrot);

void EigenSort(DKMaille<double> &d,
               DKMaille2D<double> &v,
               int n);

void multiple(double dCoupon, DKMaille<double> *Payoff);

void FromDiscountToZero(double dSpotDate,
											  DKMaille<double> &dZCRates,
											  DKMaille<double> &dZCAdjRates,
											  const DKMaille<double> &dDiscountDates,
											  const DKMaille<double> &dDiscountRates,
											  const DKMaille<double> &dAdjustedDiscountRates);