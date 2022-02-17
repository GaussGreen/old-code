#ifndef ARM_RGNOTE_H
#define ARM_RGNOTE_H



#include "ARM_result.h"



long ARM_RNGNOTE (long swapLegId,
				  double lowerBound,
				  double upperBound,
				  long recordFreq,
				  double accruedRate,
				  ARM_result& result,
				  long objId = -1);

long ARM_LIBORRNGNOTE (double startDate,
					   double endDate,
					   double lowerBound,
					   double upperBound,
					   long recordFreq,
					   long liborType,
					   double accruedRate,
					   double spread,
					   long resetFreq,
					   long payFreq,
					   long ccyId,
					   ARM_result& result,
					   long objId = -1);



#endif	// ARM_RGNOTE_H

// EOF %M%