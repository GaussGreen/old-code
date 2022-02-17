#ifndef ARM_CAPFL_H
#define ARM_CAPFL_H



#include "ARM_result.h"



extern long ARM_CAPFLOOR (long swapLegId,
						  long capOrFloor,
						  long strikeType,
						  double strike,
						  ARM_result& result,
						  long objId = -1);

extern long ARM_MATCAPFLOOR (long swapLegId,
							 double annuity,
							 double initNominal,
							 long isTRI,
							 long capOrFloor,
							 double coeff,
							 double firstTRIstrike,
							 long minStrikesId,
							 long isDigitalPayoff,
							 double increasingCoef,
							 double maxMatDate,
                             ARM_result& result,
							 long objId = -1);

extern long ARM_GEN_CAPFLOOR (long capOrFloor,
							  const CCString& delivery,
							  long liborType,
							  const CCString& currency,
							  ARM_result& result,
							  long objId = -1);

extern long ARM_GEN_CAP (const CCString& delivery,
						 long liborType,
						 const CCString& currency,
						 ARM_result& result,
						 long objId = -1);

extern long ARM_GEN_FLOOR (const CCString& delivery,
						   long liborType,
						   const CCString& currency,
						   ARM_result& result,
						   long objId = -1);

extern long ARM_LIBORCF (double startDate,
						 double endDate,
						 long isItCapOrFloor,
						 long strikeType,
						 double strike,
						 long liborType,
						 double spread,
						 long resetFreq,
						 long payFreq,
						 long currencyId,
						 ARM_result& result,
						 long objId = -1);

extern long ARM_FLEXCF (long swapLegId,
						long isItCapOrFloor,
						double strike,
						long nbEx,
						long exerciseType,
						ARM_result& result,
						long objId = -1);

extern long ARM_LIBORFLEXCF (double startDate,
							 double endDate,
							 long isItCapOrFloor,
							 double strike,
							 long nbEx,
							 long exerciseType,
							 long liborType,
							 double spread,
							 long resetFreq,
							 long payFreq,
							 long currencyId,
							 ARM_result& result,
							 long objId = -1);

extern long ARM_EXOFLEXCF (long swapLegId,
						   long isItCapOrFloor,
						   long kRefValId,
						   long nbEx,
						   long exerciseType,
						   ARM_result& result,
						   long objId = -1);

extern long ARM_CapLetPrice (long secId,
							 long modId,
							 long numEx,
						     ARM_result& result);

extern long ARM_STICKY (long swapLegId,
						long capOrFloor,
						double strike,
						const VECTOR<double>& spreadDates,
						const VECTOR<double>& spreadValues,
						long kRefValId,
						ARM_result& result,
						long objId = -1);

extern long ARM_RATCHET (long swapLegId,
						 long capOrFloor,
						 double strike,
						 const VECTOR<double>& spreadDates,
						 const VECTOR<double>& spreadValues,
						 const VECTOR<double>& correlDates,
						 const VECTOR<double>& correlValues,
						 const VECTOR<double>& fwdVolsDates,
						 const VECTOR<double>& fwdVolsValues,
						 ARM_result& result,
						 long objId = -1);

extern long ARM_ARM_SPREADOPTION (double startDate,
								  double endDate,
								  long capOrFloorId,
								  long strike_type,
								  double strike,
								  long liborType1Id,
								  long liborType2Id,
								  double weight1,
								  double weight2,
								  long dayCountId,
								  long resetFreqId,
								  long payFreqId,
								  long resetTimingId,
								  long payTimingId,
								  long ccyId,
								  ARM_result& result,
								  long objId = -1);

extern long ARM_ARM_DIGITAL (long swapLegId,
							 long isItCapOrFloorId,
							 long strikeType,
							 double strike,
							 double spread1,
							 double spread2,
							 long payoffType,
							 double payoff,
							 ARM_result& result,
							 long objId = -1);


#endif	// ARM_SWAP_H

// EOF %M%