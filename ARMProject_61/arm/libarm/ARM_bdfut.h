#ifndef ARM_BDFUT_H
#define ARM_BDFUT_H



#include "ARM_result.h"



extern long ARM_BDFUT (double delivery,
					   long underIsBd,
					   long underId,
					   double coupon,
					   double convFactor,
					   ARM_result& result,
					   long objId = -1);

extern long ARM_GetConversionFactor (long bdFutId,
									 long factId,
									 ARM_result& result);

extern long ARM_GetCheapest (long bdFutId, ARM_result& result);

extern long ARM_GILT_NOTIONNAL_BUND (const CCString &delivery,
									 long underId,
									 long notioOrGilt,
									 long market,
									 ARM_result& result,
									 long objId = -1);

extern long ARM_GILT (const CCString &delivery,
					  long underId,
	  				  ARM_result& result,
					  long objId = -1);

extern long ARM_NOTIONNAL (const CCString &delivery,
					       long underId,
	  				       ARM_result& result,
					       long objId = -1);

extern long ARM_BUND_LIFFE (const CCString &delivery,
					        long underId,
	  				        ARM_result& result,
					        long objId = -1);

extern long ARM_BUND_DTB (const CCString &delivery,
					      long underId,
	  				      ARM_result& result,
					      long objId = -1);

extern long ARM_BDFUTBASKET (double delivery,
						 long underIsBd,
						 long underId,
						 double coupon,
						 double convFactor,
						 ARM_result& result,
					     long objId = -1);


					    
#endif	// ARM_BDFUT_H



// EOF %M%