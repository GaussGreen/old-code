#ifndef ARM_BOND_H
#define ARM_BOND_H



#include "ARM_result.h"



extern long ARM_GetBondFromSummit (const CCString& rga, ARM_result& result, 
                                   long objId = -1);

extern long ARM_bond (double issueDate, double maturityDate,
					  double firstCouponDate, double couponRate,
					  double redemptionPrice, long periodicity,
					  long dayCount, long settleGap, long couponDateFlag,
					  ARM_result& result, long objId = -1);

extern long ARM_YTOPRICE (long bondId, double settlement, double yield, ARM_result& result);

extern long ARM_PTOYIELD (long bondId, double settlement, double price, ARM_result& result);

extern long ARM_SetYield (long bondId, double yield, ARM_result& result);

extern long ARM_CalcDuration (long bondId, double settlement, double yield, long flagCpn, ARM_result& result);

extern long ARM_YTOCONVEXITY (long bondId, double settlement, double actuRate, ARM_result& result);

extern long ARM_YTODURATION (long bondId, double settlement, double actuRate, long flagCpn, ARM_result& result);

extern long ARM_BDFAPRICE (long bondId, double settlement, double actuPrice, double forwardDate, 
					       double repoRate, ARM_result& result);

extern long ARM_BDREPORATE (long bondId, double settlement, 
					        double actuPrice, double forwardDate,
					        double forwardPrice,
					        ARM_result& result);


					    
#endif	// ARM_BOND_H

// EOF %M%