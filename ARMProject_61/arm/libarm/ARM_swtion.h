#ifndef ARM_SWTION_H
#define ARM_SWTION_H



#include "ARM_result.h"



extern long ARM_SWAPTION (long swapId,
					      long isRecOrPay,
					      double strike,
					      double maturity,
					      long exerciseType,
						  ARM_result& result,
					      long objId = -1);

extern long ARM_GEN_SWAPTION (double swapTerm,
							  double optionExpiry,
							  long liborType,
							  long isRecOrPay,
							  const CCString& ccy,
							  ARM_result& result,
							  long objId = -1);

extern long ARM_LIBORSWAPTION (double startDate,
							   double endDate,
				               long receiveOrPay,
							   double strike,
						       double maturity,
							   long liborType,
						       double spread,
							   long exerciseType,
						       long resetFreq,
							   long payFreq,
					           long CcyId,
							   ARM_result& result,
							   long objId = -1);

extern long ARM_EXOSWAPTION (long swapId,
			                 long isRecOrPay,
				             long xStyleId,
				             long kRefValId,
				             double swapYearTerm,
				             ARM_result& result,
					         long objId = -1);

extern long ARM_EXOCFSWAPTION (long swapId,
							   long isRecOrPay,
							   long isCapOrFloor,
							   long xStyleId,
							   long kSptionRefValId,
							   long kCFloorRefValId,
							   double cFloorPosition,
                               long IsBarrierCF,
							   ARM_result& result,
							   long objId = -1);

extern long ARM_VARFIXSWAPTION (double startDate,
								double endDate,
								long spreadsId,
								long exStyleId,
								long receiveOrPay,
								double strike,
								double maturity,
								long liborType,
								double spread,
								long resetFreq,
								long payFreq,
								long ccyId,
								ARM_result& result,
								long objId = -1);
 
extern long ARM_ARM_OPTIONALACCRUALZCBOND (double startDate,
										   double endDate,
										   double strike,
										   long nbCurPerforAcc,
										   long payFreqId,
										   long ccyId,
										   ARM_result& result,
										   long objId = -1);

extern long ARM_ARM_FlexAccretSwaption (double startDate,
										double endDate,
										double fixedRate,
										long nbCurPerforAcc,
										long receiveOrPay,
										long freqId,
										long liborTypeId,
										double spread,
										long exerciseTypeId,
										long ccyId,
										ARM_result& result,
										long objId = -1);

					    
#endif	// ARM_SWTION_H



// EOF %M%