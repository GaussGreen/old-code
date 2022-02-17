#ifndef ARM_LOCAL_OPTION_H
#define ARM_LOCAL_OPTION_H


#include <ARM\libarm\ARM_result.h>



extern long ARMLOCAL_EXOPTION (long underId,
							   long optionType,
							   long styleId,
							   long kRefValId,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_OPTION (long underId,
							 double maturityDate,
							 double strike,
							 long optionType,
							 long exerciseType,
							 long strikeType,
							 double FstXDate,
							 double PayDate,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_VolImp (long secId,
							 long modId,
							 double price,
							 ARM_result& result);

extern long ARMLOCAL_bsOption (double spot,
							   double strike,
							   double volatility,
							   double dividend,
							   double discountRate,
							   double maturity,
							   double CallPut,
							   ARM_result& result);

extern long ARMLOCAL_bsDelta (double spot,
							  double strike,
							  double volatility,
							  double dividend,
							  double discountRate,
							  double maturity,
							  double CallPut,
							  ARM_result& result);

extern long ARMLOCAL_bsVega (double spot,
							 double strike,
							 double volatility,
							 double dividend,
							 double discountRate,
							 double maturity,
							 double CallPut,
							 ARM_result& result);

extern long ARMLOCAL_bsGamma (double spot,
							  double strike,
							  double volatility,
							  double dividend,
							  double discountRate,
							  double maturity,
							  double CallPut,
							  ARM_result& result);

extern long ARMLOCAL_bsTheta (double spot,
							  double strike,
							  double volatility,
							  double dividend,
							  double discountRate,
							  double maturity,
							  double CallPut,
							  ARM_result& result);

extern long ARMLOCAL_ARM_OPTIONPORTFOLIO (long portfolioId,
										  long styleId,
										  long kRefValId,
										  long optionType,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (long calibratorId,
									                 CCString portfolioType,
									                 ARM_result& result,
									                 long objId = -1);

extern long ARMLOCAL_SpreadOptionFormula(
	double fwd1, 
	double fwd2, 
	double vol1, 
	double vol2, 
	double Correl, 
	double strike, 
	double optMat, 
	int optType, 
	int modelType, 
	double spreadVol,
	ARM_result& result);

long ARMLOCAL_SumOption(
		const long& DateStripId,
		const double& strike,
		const long& capFloor,
		const double& coeff,
		const long& coeffCurveId,
		const double& payDate,
		const long& dayCount,
		ARM_result& result,
		long objId = -1);

long ARMLOCAL_SmiledSwaption(
		const long& BaseSwaptionId,
		const VECTOR<double>& Data,
		ARM_result& result,
		long objId);

extern long ARMLOCAL_STRIPOPTION (	double asOfDate,
								    // Option data
									long underId,
									long optionType,
									long strikesId,
									long schedId,
									CCString PorS,
									long fxFixingsId,
									long leverageId,
									double leverageValue,
									// Results
									ARM_result& result,
									long objId);

long ARMLOCAL_STRIPDIGITALOPTION (	double asOfDate,
								    // Option data
									long underId,
									long optionType,
									long strikesId,
									long schedId,
									double correl,
									CCString PorS,
									long fxFixingsId,
									int callSpreadFlag,
									double epsilon,
									long payoffCurveId,
									long leverageId,
									double leverageValue,
									// Results
									ARM_result& result,
									long objId);


long ARMLOCAL_FxOptionStrip( double asOfDate,
							 long C_underlyingId,
							 long C_strikesCurveId,
							 long optionType,
							 double C_startDate,
							 double C_endDate,
							 long C_notionalId,
							 CCString C_paymentCcy,
							 long C_resetFreq,
							 long C_dayCount,
							 CCString C_resetCalendar,
							 long C_fwdRule,
							 long C_intRule,
							 long C_stubRule,
							 long C_resetGap,
							 long C_payFreq, 
							 long C_payGap, 
							 CCString C_payCalendar, 
							 long C_resetTiming, 
							 long C_payTiming,
							 CCString C_PorS,
							 long C_fxFixingsId,
							 bool isDigital,
							 int C_callSpreadFlag,
							 double C_epsilon,
							 long C_payoffCurveId,
							 long C_leverageId,
							 double leverageValue,
							 ARM_result& result,
							 long objId = -1 );

#endif /* ARM_LOCAL_OPTION_H */