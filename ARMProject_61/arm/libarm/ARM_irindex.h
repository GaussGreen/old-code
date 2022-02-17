#ifndef ARM_IRINDEX_H
#define ARM_IRINDEX_H



#include "ARM_result.h"



extern long ARM_IRINDEX (long dayCountId,
						 long frequencyId,
						 double maturity,
						 long compMethId,
						 long fwdRuleId,
						 long resetTimingId,
						 long resetGap,
						 long payTimingId,
						 long payGap,
						 long ccyId,
						 long indexType,
						 long decompFreq,
						 ARM_result& result,
						 long objId = -1);

extern long ARM_IRINDEX_MONEY_MARKET (const CCString& mmTerm,
									  const CCString& ccy,
									  ARM_result& result,
									  long objId = -1);

extern long ARM_LIBOR (long liborTypeId,
					   long ccyId,
					   long resetFreqId,
					   long payFreqId,
					   ARM_result& result,
					   long objId = -1);

extern long ARM_CMS(long CMSType,
					long liborTypeId,
					long ccyId,
	        		ARM_result& result,
					long objId = -1);

extern long ARM_FixedIndex(long dayCountId,
						   const CCString& ccy,
						   ARM_result& result,
						   long objId = -1);



#endif	// ARM_IRINDEX_H

// EOF %M%