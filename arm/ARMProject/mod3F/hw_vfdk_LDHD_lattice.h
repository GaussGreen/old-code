

#include "DKMaille.h"
#include "DKMaille2D.h"

class ARM_Matrix;
class ARM_Vector;

extern ARM_Matrix DomVolGlobalExtraction;
extern ARM_Matrix ForVolGlobalExtraction;
extern ARM_Matrix FxVolGlobalExtraction;
extern ARM_Matrix VolAnalyticsGlobalExtraction;
extern ARM_Matrix LatticeVolsGlobalExtraction;
extern ARM_Vector NbNodesPerSliceGlobalExtraction;
extern ARM_Matrix AnalyticFxOptionGlobalExtraction;
extern int FxOptIdxGlobal;

void SelectBootstrappingDates(DKMaille<double> &dNewNoticeDates,DKMaille<double> &dNoticeDates,DKMaille<double> &dSurfaceX)	;

void InterpolationFXVol(DKMaille<double> &dNewVolStripDates, 
												DKMaille<double> &dFXVolDatesPrime,
												DKMaille<double> &dFXVolPrime,
												DKMaille<double> &dFXVolDates, 
												DKMaille<double> &dFXVol);

DKMaille<double>  Lattice_HWVFDK_3F_(DKMaille<double> dLatticeGeometryData,
																		 double dNumTimeLinesBeforeFirstNotice,
																		 double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
																		 double dNumTimeLinesPerYear, 
																		 double dSpotDate, 
																		 DKMaille<double> dBaseDates,
																		 DKMaille<double> dBaseRates,
																		 DKMaille<double> dForeignDates,
																		 DKMaille<double> dForeignRates,
																		 DKMaille<double> dRedemptionData,
																		 double dStrike,
																		 double dType,
																		 double dOptionExpiry,
																		 DKMaille<double> dStdDevBaseX,
																		 DKMaille<double> dStdDevBaseY,
																		 DKMaille2D<double> dStdDevBaseZ,
																		 double dMeanReversionBase,  
																		 DKMaille<double> dStdDevForeignX,
																		 DKMaille<double> dStdDevForeignY,
																		 DKMaille2D<double> dStdDevForeignZ,
																		 double dMeanReversionForeign,  
																		 double dSpotFX,
																		 DKMaille<double> dSpotFXVolDatesTD,																		 
																		 DKMaille<double> dSpotFXVolTD,
																		 double dBaseForeignCorrelation,
																		 double dBaseSpotFXCorrelation,
																		 double dForeignSpotFXCorrelation,
																		 double dProductModelCode,
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
																		 double dSmoothing,
																		 double dSurvivalCalculation,
																		 DKMaille2D<double> dBoosterData,
																		 DKMaille<double> dBaseRatesNoBasis,
																		 DKMaille<double> dForeignRatesNoBasis,
																		 double dSmileParameterBase,
																		 double dSmileParameterForeign,
																		 double dIsSwaptionCalibrationWithBasis);



