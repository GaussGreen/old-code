#ifndef ARM_LOCAL_UTILITIES_H
#define ARM_LOCAL_UTILITIES_H

#include <vector>
#include <ARM/libarm/ARM_result.h>

#include <GP_Base/gpbase/curve.h>
#include <GP_Base/gpbase/curvetypedef.h>

extern long ARMLOCAL_CQSOCREATE(double dStartDate, 
						 double dEndDate,
						 int iFrequency,
						 int iDayCounter,
						 int iIsAdjusted,
						 int iResetType,
						 int iResetGap,
						 CCString sResetCalendar,
						 CCString sPaymtCalendar,						
						 CCString domCurrency,
						 CCString forCurrency,
						 CCString index1,
						 CCString index2,						 
						 const ARM::ARM_Curve& nominal, 
						 const ARM::ARM_Curve& strikes, 
						 const ARM::ARM_Curve& leveragesShort, 
						 const ARM::ARM_Curve& leveragesLong, 
						 const ARM::ARM_Curve& cpnMin, 
						 const ARM::ARM_Curve& cpnMax, 
						 const ARM::ARM_Curve& margins,
						 const ARM::ARM_Curve& fees,
						 ARM_result& result,
						 long& objId
								);


#endif /* ARM_LOCAL_UTILITIES_H */