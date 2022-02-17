#ifndef ARM_IASEC_H
#define ARM_IASEC_H



#include "ARM_result.h"



extern long ARM_IASEC (long underlyingId,
					   long iaCtrlId,
					   long refValId,
					   long iaCtrlType,
					   ARM_result& result,
					   long objId = -1);



#endif	// ARM_IASEC_H

// EOF %M%