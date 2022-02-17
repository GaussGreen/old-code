#ifndef ARM_XSTYLE_H
#define ARM_XSTYLE_H



#include "ARM_result.h"



extern long ARM_EUROPEANXSTYLE (double xdate,
						        ARM_result& result,
						        long objId = -1);

extern long ARM_AMERICANXSTYLE (double xStartDate,
						        double xEndDate,
						        ARM_result& result,
						        long objId = -1);

extern long ARM_BERMUDANXSTYLE (VECTOR<double>& xDates,
						        ARM_result& result,
						        long objId = -1);

extern long ARM_CUSTOMXSTYLE (VECTOR<double>& xStartDates,
							  VECTOR<double>& xEndDates,
						      ARM_result& result,
						      long objId = -1);


					    
#endif	// ARM_XSTYLE_H



// EOF %M%