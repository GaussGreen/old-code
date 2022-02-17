#ifndef ARM_BARRIER_H
#define ARM_BARRIER_H



#include "ARM_result.h"



extern long ARM_CONSTBARRIER (long underlyingId,
							  long tAssetId,
							  double maturity,
							  double barrier,
							  long upDown,
							  long inOut,
							  long triggerVar,
							  double rebate,
							  double firstX,
							  ARM_result& result,
							  long objId = -1);

extern long ARM_BARRIER (long underlyingId,
					     long tAssetId,
						 long xStyleId,
						 long refValId,
						 long upDown,
						 long inOut,
						 long triggerVar,
						 double rebate,
						 ARM_result& result,
						 long objId = -1);



#endif	// ARM_BARRIER_H

// EOF %M%