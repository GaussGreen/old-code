#ifndef ARMLOCAL_BDFUT_H
#define ARMLOCAL_BDFUT_H



#include "ARM_result.h"



extern long ARMLOCAL_BDFUT (double delivery,
					   long underIsBd,
					   long underId,
					   double coupon,
					   double convFactor,
					   ARM_result& result,
					   long objId = -1);

extern long ARMLOCAL_GetConversionFactor (long bdFutId,
									 long factId,
									 ARM_result& result);

extern long ARMLOCAL_GetCheapest (long bdFutId, ARM_result& result);

extern long ARMLOCAL_GILT_NOTIONNAL_BUND (const CCString &delivery,
									 long underId,
									 long notioOrGilt,
									 long market,
									 ARM_result& result,
									 long objId = -1);

extern long ARMLOCAL_GILT (const CCString &delivery,
					  long underId,
	  				  ARM_result& result,
					  long objId = -1);

extern long ARMLOCAL_NOTIONNAL (const CCString &delivery,
					       long underId,
	  				       ARM_result& result,
					       long objId = -1);

extern long ARMLOCAL_BUND_LIFFE (const CCString &delivery,
					        long underId,
	  				        ARM_result& result,
					        long objId = -1);

extern long ARMLOCAL_BUND_DTB (const CCString &delivery,
					      long underId,
	  				      ARM_result& result,
					      long objId = -1);

extern long ARMLOCAL_BDFUTBASKET (double delivery,
						 long underIsBd,
						 long underId,
						 double coupon,
						 double convFactor,
						 ARM_result& result,
					     long objId = -1);


					    
#endif	// ARMLOCAL_BDFUT_H



// EOF %M%