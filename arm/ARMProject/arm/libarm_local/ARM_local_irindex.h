#ifndef ARM_LOCAL_IRINDEX_H
#define ARM_LOCAL_IRINDEX_H


extern long ARMLOCAL_LIBOR (long liborTypeId,
							bool ccyIsObject,
							const CCString& ccyName,
							long resetFreqId,
							long payFreqId,
                            long daycount,
							long intRuleId,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_IRINDEX (long dayCountId,
							  long payFreqId,
							  double maturity,
							  long compMethId,
							  long fwdRuleId,
							  long resetTimingId,
							  long resetGap,
							  long payTimingId,
							  long payGap,
							  bool ccyIsObject,
							  const CCString& ccyName,
							  long indexType,
							  long decompFreq,
							  long intRuleId,
							  long resetFreqId,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_MultiIrindex(VECTOR<CCString>& irindexVect,
								  VECTOR<double>& weightVect,
								  ARM_result& result, 
								  long objId = -1);

extern long ARMLOCAL_FixedIndex(long dayCountId,
								const CCString& ccy,
								ARM_result& result,
								long objId = -1);

extern long ARMLOCAL_IRINDEX_MONEY_MARKET (const CCString& mmTerm,
										   const CCString& ccy,
										   ARM_result& result,
										   long objId = -1);

extern long ARMLOCAL_CMS(long CMSType,
						 long liborTypeId,
						 bool ccyIsObject,
						 const CCString& ccyName,
						 ARM_result& result,
						 long objId = -1);

#endif /* ARM_LOCAL_IRINDEX_H */