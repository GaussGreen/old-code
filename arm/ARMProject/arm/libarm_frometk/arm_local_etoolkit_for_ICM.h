#ifndef ARM_LOCAL_ETOOLKIT_FOR_ICM_H
#define ARM_LOCAL_ETOOLKIT_FOR_ICM_H

#include <libCCTools++/CCString.h>
// #include "ARM_local_eToolkit.h"


// A cause de windows.h
#ifdef GetObject
#undef GetObject
#endif


#include <libCCTools++/CCString.h>

class ARM_Date;
class ARM_result;


extern long etoolkit_GetDefProbCurve (const CCString& issuer,
							   ARM_Date asof,
							   const CCString& CurveId,
							   CCString& xmlResponse,
							   CCString& messageList);



extern long etoolkit_GetVolFromSummit (ARM_Date asof,
							   const CCString& IssuerId,
							   const CCString& curveName,
							   const CCString& currency,
							   CCString& xmlResponse,
							   CCString& messageList);


extern long etoolkit_GetCorrFromSummit (ARM_Date asof,
							   const CCString& issuerA,
							   const CCString& issuerB,
							   const CCString& curveId,
							   CCString& xmlResponse,
							   CCString& messageList);

extern long etoolkit_GetRecoveryCurve (const CCString& SecId,
							   ARM_Date asof,
							   const CCString& CurveId,
							   CCString& xmlResponse,
							   CCString& messageList);

extern long etoolkit_getXML_ZC_DP_FromSummit(const CCString& Issuer,
											const CCString& currency,
											const CCString& cvName,
											ARM_Date aSdate,
											CCString& xmlResponse,
											CCString& messageList);

extern long  ARMLOCAL_ParseDefProbCurveMktData (const char* chaineXML,
											 ARM_Date CurveDate, 
											 ARM_Currency* ccy,
										     CCString& SecId);

						   
#endif
