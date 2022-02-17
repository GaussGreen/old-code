#ifndef ARM_LOCAL_SWTION_H
#define ARM_LOCAL_SWTION_H


class ARM_result;


extern long ARMLOCAL_SWAPTION (long swapId,
							   long isRecOrPay,
							   long strikeType,
							   double strike,
							   double maturity,
							   long exerciseType,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_SwaptionFromExpiry(CCString& optionExpiry,
										 CCString& swapTerm,
										 long liborType,
										 long receiveOrPay,
										 double strike,
										 long spreadType,
										 double spread,
										 bool ccyIsObject,
										 const CCString& CcyName,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_LIBORSWAPTION (double startDate,
									double endDate,
									long receiveOrPay,
									double strike,
									double maturity,
									long liborType,
									long spreadType,
									double spread,
									long exerciseType,
									long resetFreq,
									long payFreq,
									bool ccyIsObject,
									const CCString& CcyName,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_EXOSWAPTION (long swapId,
								  long isRecOrPay,
								  long xStyleId,
								  long kRefValId,
								  double swapYearTerm,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_VARFIXSWAPTION (double startDate,
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

extern long ARMLOCAL_ARM_OPTIONALACCRUALZCBOND (double startDate,
												double endDate,
												double strike,
												long nbCurPerforAcc,
												long payFreqId,
												long ccyId,
												ARM_result& result,
												long objId = -1);

extern long ARMLOCAL_EXOCFSWAPTION (long swapId,
									long isRecOrPay,
									long isCapOrFloor,
									long xStyleId,
									long kSptionRefValId,
									long kCFloorRefValId,
									double cFloorPosition,
									long IsBarrierCF,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_FlexAccretSwaption (double startDate,
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
//-----------------------------------
extern double ARMLOCAL_SwaptionStickyDelta(long swaptionId,
										   long modelId,
										   long perturbeDiscountCurvId,
								           ARM_result& result);

#endif /* ARM_LOCAL_SWTION_H */
