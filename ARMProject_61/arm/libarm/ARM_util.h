#ifndef ARM_UTIL_H
#define ARM_UTIL_H



#include "ARM_result.h"




extern long ARM_GetPrtyByName (long objId, const CCString& prtByName, 
                               long nbArg, ARM_result& result);


extern long ARM_GetVolOrRatesRange(double date1,
       					           double date2,
					               const CCString& ccy,
					               const CCString& index,
					               const CCString& cvName,
					               const CCString& expiry,
					               const CCString& matu,
					               long yieldOrValId,
					               long calcModId,
					               const CCString& volType,
					               const CCString& outFile,
					               ARM_result& result);

extern long ARM_GetHistoricalVol(double date1,
					             double date2,
					             const CCString& ccy,
					             const CCString& index,
					             const CCString& cvName,
					             const CCString& expiry,
					             const CCString& matu,
					             long yieldOrValId,
					             long calcModId,
					             const CCString& volType,
					             ARM_result& result);

extern long ARM_GetAsOfVolOrRate(double date1,
  					             const CCString& ccy,
					             const CCString& index,
					             const CCString& cvName,
					             const CCString& expiry,
					             const CCString& matu,
					             long yieldOrValId,
					             long calcModId,
					             const CCString& volType,
					             ARM_result& result);

extern long ARM_ComputeHistoCorrel(double date1,
 					               double date2,
					               const CCString& ccy1,
					               const CCString& index1,
					               const CCString& cvName1,
					               const CCString& expiry1,
					               const CCString& matu1,
					               long yieldOrValId1,
					               long calcModId1,
					               const CCString& volType1,
                                   const CCString& ccy2,
					               const CCString& index2,
					               const CCString& cvName2,
					               const CCString& expiry2,
					               const CCString& matu2,
					               long yieldOrValId2,
					               long calcModId2,
					               const CCString& volType2,
					               ARM_result& result);

extern long ARM_GetHistoFwdRates (double date1,
		        			      double date2,
					              const CCString& expiry,
					              const CCString& matu,
                                  const CCString& ccy,
					              const CCString& outFile,
					              ARM_result& result);

extern long ARM_GetFwdRatesMatrix (double date,
		        			      const CCString& ccy,
					              ARM_result& result);


#endif	// ARM_UTIL_H

// EOF %M%