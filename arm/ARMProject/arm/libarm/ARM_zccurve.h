#ifndef ARM_ZCCURVE_H
#define ARM_ZCCURVE_H

#include "ARM_result.h"



extern long ARM_zclint (const VECTOR<double>& matu,
						const VECTOR<double>& rate,
						long meth,
					    double aDate,
						const CCString& ccy,
						long interpMethId,
						ARM_result& result,
						long objId = -1);

extern long ARM_GetZCFromSummit (const CCString& index,
								 const CCString& currency,
								 const CCString& cvName,
								 double aSdate,
								 ARM_result& result,
								 long objId = -1);

extern long ARM_DiscountPrice (long idCurve,
							   double matu,
							   ARM_result& result);

extern long ARM_DiscountYield (long idCurve,
							   double matu,
							   long meth,
							   ARM_result& result);

extern long ARM_CreateZCSwapFutInt (double Date,
									VECTOR<CCString>& matu,
									VECTOR<double>& rate,
									long MMVsFut,
									long SwapVsFut,
									long Raw,
									long interp,
									const CCString& Ccy,
									ARM_result& result,
									long objId = -1);

extern long ARM_CreateZCSwapInt (double Date,
								 VECTOR<CCString>& matu,
								 VECTOR<double>& rate,
								 long MMVsFut,
								 long SwapVsFut,
								 long Raw,
								 long interp,
								 const CCString& Ccy,
								 ARM_result& result,
								 long objId = -1);

extern long ARM_C_CreateZCSwapInt (double Date, 
								   CCString* C_matu,
								   long matu_size,
								   double* C_rate,
								   long rate_size,
								   long MMVsFut,
								   long SwapVsFut,
								   long Raw, 
								   long interp,
								   const CCString& Ccy,
								   ARM_result& result,
								   long objId);

extern long ARM_ForwardYield (long idCurve,
							  double matu1,
							  double matu2,
							  long meth,
							  long adjId,
							  ARM_result& result);

extern long ARM_ForwardPrice (long idCurve,
							  double matu1,
							  double matu2,
							  ARM_result& result);

extern long ARM_zcvsk (const VECTOR<double>& param,
					   double date,
                       ARM_result& result,
					   long objId = -1);

extern long ARM_zcflat (const double zeroFlat,
						double date, 
                        ARM_result& result,
						long objId = -1);

extern long ARM_zcspl (const VECTOR<double>& param,
					   double date, 
                       ARM_result& result,
					   long objId = -1);

extern long ARM_CreateZCTAMInt (double date,
								VECTOR<CCString>& matu, 
								VECTOR<double>& rate,
								double mean_rates, 
								long raw,
								long interp, 
								long lastBucketInt,
								const CCString& Ccy,
								ARM_result& result,
								long objId = -1);

extern long ARM_CreateZCCashInt (double date,
								 VECTOR<CCString>& matu,
								 VECTOR<double>& rate,
								 VECTOR<long>& bondsId,
								 VECTOR<double>& yields,
								 long MMVsFut,
								 const CCString& Ccy,
								 ARM_result& result,
								 long objId = -1);

extern long ARM_zcsplicub (VECTOR<double>& matu,
						   VECTOR<double>& rate,
					       long meth,
						   double date,
						   long lastBucket,
					       ARM_result& result,
						   long objId = -1);

extern long ARM_zccubdiff (VECTOR<double>& matu,
						   VECTOR<double>& rate,
						   long meth,
						   double date,
						   ARM_result& result,
						   long objId = -1);

extern long ARM_zcswapcubdiff (double date,
							   VECTOR<CCString>& matu,
							   VECTOR<double>& rate,
							   long mmVsFut,
							   long swapVsFut,
							   long raw,
							   long interp,
							   const CCString& ccy,
							   ARM_result& result,
							   long objId = -1);

extern long ARM_zcspreaded (long zcSprId,
							long zcInitId,
							double date,
							long MMFreq,
							long SwapFreq,
							long ccyId,
							ARM_result& result,
							long objId = -1);

extern long ARM_GetInitialCurveFromSummit (const CCString& index,
										   const CCString& currency, 
										   const CCString& cvName,
										   double aSdate,
										   long AdjOrNotId,
										   const CCString& outFile,
										   ARM_result& result);

extern long ARM_CreateZCSwapIntSmooth (double Date,
									   VECTOR<CCString>& matu, 
									   VECTOR<double>& rate,
									   long MMVsFut,
									   long SwapVsFut,
									   long Raw, 
									   long interp,
									   const CCString& Ccy,
									   double lambda,
									   long prec,
									   ARM_result& result,
									   long objId = -1);

extern long ARM_CreateZCSwapFutIntSmooth (double Date,
										  VECTOR<CCString>& matu, 
										  VECTOR<double>& rate,
										  long MMVsFut,
										  long SwapVsFut,
										  long Raw, 
										  long interp,
										  const CCString& Ccy,
										  double lambda,
										  long prec,
										  ARM_result& result,
										  long objId = -1);

extern long ARM_ZCINTSMOOTH (const VECTOR<double>& matu,
							 const VECTOR<double>& rate,
							 double aDate,
							 long meth,
							 double lambda,
							 long prec,
							 ARM_result& result,
							 long objId = -1);

extern long ARM_ZCINTSMOOTH (long inCvId,
							 double lambda,
							 long prec,
							 ARM_result& result,
							 long objId = -1);

extern long ARM_CreateTOYNYZCSwapInt (double Date,
									  VECTOR<CCString>& matu,
									  VECTOR<double>& rate,
									  long MMVsFut,
									  long SwapVsFut,
									  long Raw,
									  long interp,
									  const CCString& Ccy,
									  long frq,
									  ARM_result& result,
									  long objId = -1);

extern long ARM_ARM_BumpCurve (long cvId,
							   VECTOR<CCString>& matu,
							   VECTOR<double>& epsilon,
							   ARM_result& result,
							   long objId = -1);

#endif	// ARM_ZCCURVE_H

// EOF %M%
