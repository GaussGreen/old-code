#ifndef ARM_REV_H
#define ARM_REV_H



#include "ARM_result.h"



extern long ARM_REVERSE (long structSwapLegId,
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

extern long ARM_REVERSE_CALENDAR (long structSwapLegId,
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

extern long ARM_REVERSECOUPON (double strike,
							   double power,
							   double callput,
							   long size,
							   ARM_result& result,
							   long objId = -1);

extern long ARM_STRUCTREVERSECOUPON(const VECTOR<double>& date,
        							const VECTOR<double>& strike,
							        const VECTOR<double>& power,
							        const VECTOR<double>& callput,
							        double xo,
                                    const VECTOR<double>& floor,
                                    const VECTOR<double>& cap,
							        ARM_result& result,
							        long objId = -1);


extern long ARM_STRUCTREVERSECOUPON2 (const VECTOR<double>& date,
									  const VECTOR<double>& strike,
									  const VECTOR<double>& power,
									  const VECTOR<double>& callput,
									  double xo,
									  ARM_result& result,
									  long objId = -1);



#endif	// ARM_REV_H

// EOF %M%