#ifndef ARM_LOCAL_CCY_H
#define ARM_LOCAL_CCY_H

class ARM_result;

extern long ARMLOCAL_ISOCCY (const CCString& cname,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_CCY (const CCString& name,
						  long idCurve,
						  double crossValue,
						  long daycount,
						  ARM_result& result,
						  long objId = -1);

extern long ARMLOCAL_GetSpotDays (const CCString&  ccy,
								  ARM_result& result);

extern long ARMLOCAL_GetInfoFromCcy(long ccyId,
									const CCString& type,
									ARM_result& result);

#endif	// ARM_LOCAL_CCY_H

// EOF %M%