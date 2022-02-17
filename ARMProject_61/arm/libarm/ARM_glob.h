#ifndef ARM_GLOB_H
#define ARM_GLOB_H



#include "ARM_result.h"


extern double ARMDATE2XLDATE (const CCString& armdate);

extern CCString XLDATE2ARMDATE (double xldate);

extern long ARM_ARM_Price (long secId,
						   long modId,
						   ARM_result& result);

extern long ARM_ARM_View (long instId,
						  const CCString& sockId,
						  ARM_result& result);

extern long ARM_GetFRMShortRateVols(long modId,
									const CCString& outFile, 
                                    ARM_result& result);

extern long ARM_FreeObject (long secId,
							ARM_result& result);

extern long ARM_FreeAllObjects (ARM_result& result);

extern long	ARM_ExitArm ();

extern long ARM_ARM_GetDefaultCurrency (ARM_result& result);

extern long ARM_ARM_SetDefaultCurrency (const CCString& isoCCy,
										ARM_result& result);

extern long ARM_NextBusinessDay (double date,
								 const CCString& cal,
								 long days,
								 ARM_result& result);

extern long ARM_IsBusinessDay (double date,
							   const CCString& isoccyname,
							   ARM_result& result);

extern long ARM_ADJUSTTOBUSDATE (double date,
								 const CCString& currency,
								 long ruleId,
								 ARM_result& result);

extern long ARM_ARM_Accrued (long secId,
							 double fwdDate,
							 long modId,
							 ARM_result& result);

extern long ARM_Sensitivity (long secId,
							 long modId,
							 long paramId,
							 ARM_result& result);

extern long ARM_CvSensitivity (long secId,
							   long modId,
							   long paramId,
							   long viewFlagId,
							   const CCString& id,
							   ARM_result& result);

extern long ARM_SetMarketPrice (long id,
								double price,
								ARM_result& result);

extern long ARM_IMPLIEDVOL (long instId,
							long modelId,
							double price,
							ARM_result& result);

extern long ARM_ParallelShift (long secId,
							   double value,
							   ARM_result& result,
							   long objId = -1);

extern long ARM_SetNotional (long secId,
							 long rId,
							 double percentRemainder,
							 ARM_result& result);

extern long ARM_NextBusinessDayDt (double date,
								   const CCString& cal,
								   long days,
								   ARM_result& result);

extern long ARM_GetExpiry (long secId,
						   ARM_result& result);

extern long ARM_FwdPrice (long secId,
						  long modId,
						  double fwdDate,
						  ARM_result& result);

extern long ARM_BSSpot (long secId,
						long modId,
						double date,
						ARM_result& result);

extern long ARM_SummitSwapInfo (const CCString& tradeId,
								const CCString& curveId,
								double valoDate,
								const CCString& whatTo,
								ARM_result& result);

extern long ARM_FxConvert (const CCString& ccy1,
					       const CCString& ccy2,
					       double asOfDate,
					       double amount,
					       const CCString& cvname,
					       ARM_result& result);

extern long ARM_XCccyAdjustment( long startDate,
								long enDate,
								long payFreq,
								const CCString& domCcy,
								long forIndexTypeId,
								const CCString& forCcy,
								long spreadsId,
								long zcDomId,
								long discDomId,
								long zcForId,
								long discForId,
								double FX,
								long couponId, 
								long domDc,
								long forDc,
								ARM_result& result,
								long objId = -1);


extern long ARM_EvalSummitAsset (const CCString& tradeId,
								 const CCString& tradeType,
								 const CCString& curveId,
								 double asOfDate,
								 long assetId,
								 long isNet,
								 ARM_result &result);

extern long ARM_KImp (long secId,
					  long modId,
					  double price,
					  long param,
					  ARM_result &result);

extern long ARM_SetPrice (long instId,
						  double price,
						  ARM_result& result);

extern long ARM_ExProba (long secId,
						 long modId,
						 long numEx,
						 ARM_result& result);

extern long ARM_SummitValueGreeks (const CCString& tradeId,
								   const CCString& cvId,
								   const CCString& tradeType,
								   double valoDate,
								   const CCString& greek,
								   const CCString& analytic,
								   ARM_result& result);

extern long ARM_SummitValueTrade (const CCString& tradeId,
								  const CCString& tradeType,
								  const CCString& curveId,
								  double date,
								  long isNetVal,
								  long version,
								  ARM_result& result);
							
extern long ARM_ARM_GetPID (ARM_result& result);

extern long ARM_ARM_SCredit (const CCString& SummitFilter,
							 const CCString& SummitCurveId,
							 const CCString& FileName,
							 ARM_result& result);

extern long ARM_ARM_BetweenDates (long date1,
								  long date2,
								  long daycountId,
								  long isYearFrac,
								  ARM_result& result);

extern long ARM_ARM_ADDYEARS  (long date,
							   long nb,
							   long ruleId,
							   const CCString& Ccy,
							   ARM_result& result);

extern long ARM_ARM_ADDMONTHS  (long date,
								long nb,
								long ruleId,
								const CCString& Ccy,
								ARM_result& result);

extern long ARM_ARM_ADDPERIOD  (long date,
								long freq,
								const CCString& ccy,
								long nbPeriods,
								long adjRuleId,
								ARM_result& result);

extern long ARM_ARM_INTERPOL (const VECTOR<double>& vecX,
							  const VECTOR<double>& vecY,
							  double X,
							  long interpId,
							  ARM_result& result);

extern long ARM_DiscountPriceRefvalue(long zcId,
									  long refvalId,
									  long ccyId,
									  double startDate,
									  double enddate,
									  ARM_result& result);

extern long ARM_ARM_Price_OptUnder (long secId,
									long modId,
									ARM_result& result);

extern long ARM_ARM_ClonedAndSetNotional (long secId,
										  long rId,
										  double percentRemainder,
										  ARM_result& result,
										  long objId = -1);

extern long ARM_ARM_DisplayScheduleValues (long instId,
										   long valuesType,
										   long recId,
										   long modelId,
										   const CCString& sockId,
										   ARM_result& result);

extern long ARM_ARM_DisplayScheduleDates (long instId,
										  long datesType,
										  long recId,
										  const CCString& sockId,
										  ARM_result& result);

extern long ARM_ARM_TRIANGULARINTERPOL(const VECTOR<double>& vecX,
									   const VECTOR<double>& vecY,
									   const VECTOR<double>& matZ,
									   double X,
									   double Y,
									   ARM_result& result);
#endif	// ARM_GLOB_H

// EOF %M%