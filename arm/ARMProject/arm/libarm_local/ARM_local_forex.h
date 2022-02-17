#ifndef ARM_LOCAL_FOREX_H
#define ARM_LOCAL_FOREX_H

class ARM_result;

extern long ARMLOCAL_FOREX (long LeftCcyId,
							long RightCcyId,
                            double SpotValue,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_FOREX ( const CCString& LeftCcyName,
							 const CCString& RightCcyName,
							 double SpotValue,
							 ARM_result& result,
							 long objId = -1);

#endif	// ARM_LOCAL_FOREX_H

// EOF %M%