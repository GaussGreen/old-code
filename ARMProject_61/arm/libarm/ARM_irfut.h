#ifndef ARM_IRFUT_H
#define ARM_IRFUT_H



#include "ARM_result.h"



extern long ARM_IRFUT (double delivery, long idUnderlying,
					   ARM_result& result, long objId = -1);

extern long ARM_THREE_MONTH_FUT (const CCString& delivery, long market, const CCString& ccy,
							     ARM_result& result, long objId = -1);

extern long ARM_FUT_PIBOR (const CCString& delivery, ARM_result& result, long objId = -1);

extern long ARM_FUT_SHORTSTERLING (const CCString& delivery, ARM_result& result, long objId = -1);



#endif	// ARM_IRFUT_H

// EOF %M%