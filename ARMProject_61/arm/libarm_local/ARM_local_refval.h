#ifndef ARM_LOCAL_REFVAL_H
#define ARM_LOCAL_REFVAL_H


extern long ARMLOCAL_REFVALUE (VECTOR<double>& dates,
							   VECTOR<double>& values,
							   VECTOR<double>& values2,
							   long valueType,
							   long conversion,
							   long calcMethod,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_CONSTREFVALUE (double value,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_IATHREELEVREFVAL (double value,
									   double level0,
									   double amort0,
									   double level1,
									   double amort1,
									   double level2,
									   double amort2,
									   ARM_result& result,
									   long objId = -1);

extern long ARMLOCAL_CptRefValue (long refvalId,
								  double date,
								  ARM_result& result);

extern long ARMLOCAL_SumRefValue (long refval1Id,
								  long refval2Id,
								  double coef,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_DisplayRefValue(long refvalId,
									 bool isDate,
									 ARM_result& result);

#endif /* ARM_LOCAL_REFVAL_H */

