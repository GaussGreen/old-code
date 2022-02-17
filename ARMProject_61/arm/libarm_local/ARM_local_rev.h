#ifndef ARM_LOCAL_REV_H
#define ARM_LOCAL_REV_H



#include "ARM_result.h"



extern long ARMLOCAL_REVERSE (long structSwapLegId,
							  long classSwapLegId,
							  long ReceiveOrPay,
							  long couponId,
							  long exeId,
							  long redempId,
							  long classRedempId,
							  double dualDate,
							  double dualStrike,
							  const CCString& dualFlag,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_STRUCTREVERSECOUPON(const VECTOR<double>& date,
										 const VECTOR<double>& strike,
										 const VECTOR<double>& power,
										 const VECTOR<double>& callput,
										 double xo,
										 const VECTOR<double>& floor,
										 const VECTOR<double>& cap,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_REVERSE_CALENDAR (long structSwapLegId,
									   long classSwapLegId,
									   long ReceiveOrPay,
									   long couponId,
									   long exeId,
									   long redempId,
									   long classRedempId,
									   double dualDate,
									   double dualStrike,
									   const CCString& dualFlag,
									   const VECTOR<double>& dStartDates, 
									   const VECTOR<double>& dEndDates, 
									   const VECTOR<double>& dFixingDates, 
									   const VECTOR<double>& dPaymentDates, 
									   const VECTOR<double>& fStartDates, 
									   const VECTOR<double>& fEndDates, 
									   const VECTOR<double>& fFixingDates, 
									   const VECTOR<double>& fPaymentDates, 
									   ARM_result& result,
									   long objId = -1);

extern long ARMLOCAL_REVERSENOTIONAL_CALENDAR (long structSwapLegId,
											   long classSwapLegId,
											   long ReceiveOrPay,
											   long couponId,
											   long exeId,
											   long redempId,
											   long classRedempId,
											   long notExchFlag,
											   double dualDate,
											   double dualStrike,
											   const CCString& dualFlag,
											   const VECTOR<double>& dStartDates,
											   const VECTOR<double>& dEndDates,
											   const VECTOR<double>& dFixingDates,
											   const VECTOR<double>& dPaymentDates,
											   const VECTOR<double>& fStartDates,
											   const VECTOR<double>& fEndDates,
											   const VECTOR<double>& fFixingDates,
											   const VECTOR<double>& fPaymentDates,
											   ARM_result& result,
											   long objId = -1);

extern long ARMLOCAL_POWERREVERSE (long initPeriodLegId,
								   long realFundLegId,
								   long fxUnderLegId,
								   long fxNumLegId,
								   const VECTOR<double>& dNoticeDates,
								   const VECTOR<double>& dCancelDates,
								   long FX0type,
								   double FX0,
								   double fxStep,
								   long capValueType,
								   double capValue,
								   double capStepValue,
								   long floorValueType,
								   double floorValue,
								   double floorStepValue,
								   long dualOptionFlag,
								   double dualOptionStrike,
								   double redempNoticeDate,
								   double lastLiborFixing,
								   double lastFxSpotFixing,
                                   ARM_result& result,
								   long objId = -1);
extern long ARMLOCAL_DatePowerReverseGetData(
	long prcsId,
	long dataType,
	VECTOR<double>& Data,
	ARM_result& result );
extern long ARMLOCAL_GETPRCSDATA (long prcsId,
								  ARM_result& result);

extern long ARMLOCAL_DELETENEXTCALLFROMPRCS (long prcsId,
											 double asOf,
											 ARM_result& result,
											 long objId = -1);

extern long ARMLOCAL_GETOPTIONDATES(long prcsId,
									const CCString& dateType,
									ARM_result& result);

extern long ARMLOCAL_GETDUALOPTIONSTRIKE(long prcsId,
										 ARM_result& result);

#endif /* ARM_LOCAL_REV_H */