#ifndef ARM_CCY_H
#define ARM_CCY_H


#include "ARM_result.h"


extern long ARM_CCY (const CCString& name,
					 long idCurve,
					 double crossValue,
					 long daycount,
					 ARM_result& result,
					 long objId = -1);

extern long ARM_ISOCCY (const CCString& name,
						ARM_result& result,
						long objId = -1);

extern long ARM_GetDefaultIndexFromCurrency(const CCString& name,
											ARM_result& result);

extern long ARM_GetPayCalName (const CCString& name,
							   long idxtype,
							   ARM_result& result);

extern long ARM_GetCcyName (const long ccyid,
							ARM_result& result);

extern long ARM_GetSpotDays (const CCString&  ccy,
							 ARM_result& result);

#endif	// ARM_CCY_H

// EOF %M%