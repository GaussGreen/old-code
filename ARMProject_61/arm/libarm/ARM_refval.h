#ifndef ARM_REFVAL_H
#define ARM_REFVAL_H


#include "ARM_result.h"



extern long ARM_CONSTREFVALUE (double value,
							   ARM_result& result,
							   long objId = -1);
                                   
extern long ARM_REFVALUE (VECTOR<double>& dates,
						  VECTOR<double>& values,
						  VECTOR<double>& values2,
						  long valueType,
						  long conversion,
						  long calcMethod,
						  ARM_result& result,
						  long objId = -1);

extern long ARM_IATHREELEVREFVAL (double value,
								  double level0,
								  double amort0,
								  double level1,
								  double amort1,
								  double level2,
								  double amort2,
								  ARM_result& result, long objId = -1);


					    
#endif	// ARM_REFVAL_H

// EOF %M%