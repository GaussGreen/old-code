#ifndef ARM_LOCAL_IRFUT_H
#define ARM_LOCAL_IRFUT_H


extern long ARMLOCAL_IRFUT (double delivery, long idUnderlying,
					   ARM_result& result, long objId = -1);

extern long ARMLOCAL_THREE_MONTH_FUT (const CCString& delivery, long market, const CCString& ccy,
							     ARM_result& result, long objId = -1);


#endif /* ARM_LOCAL_IRFUT_H */