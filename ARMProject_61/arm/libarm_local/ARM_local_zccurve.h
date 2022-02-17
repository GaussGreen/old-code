#ifndef ARM_LOCAL_ZCCURVE_H
#define ARM_LOCAL_ZCCURVE_H

#include "armglob.h"
#include <gpbase\datestrip.h>

class ARM_result;
class ARM_ZeroLInterpol;
class ARM_Date;

using namespace ARM;

extern long ARMLOCAL_CreateZCSwapInt (double Date,
									  VECTOR<CCString>& matu,
									  VECTOR<double>& rate,
									  long MMVsFut,
									  long SwapVsFut,
									  long Raw,
									  long interp,
									  const CCString& Ccy,
									  long swapFrqId,
									  long fixDayCount,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_zclint(const VECTOR<double>& matu,
							const VECTOR<double>& rate,
							long meth,
							double aDate,
							const CCString& sCcy,
							long interpMeth,
                            int    matuAreDoubles, // 1: YES, O: NO
                            ARM_CRV_TERMS& sTerms,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_shiftzclint(double value,
                                 long nbplot,
                                 const VECTOR<double>& matu,
					             const VECTOR<double>& rate,
					             long meth,
					             double aDate,
					             const CCString& sCcy,
					             long interpMeth,
					             ARM_result& result,
					             long objId = -1);

extern long ARMLOCAL_zcflat (const double zeroFlat,
							 double date,
							 const CCString& ccy,
                             ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_DiscountPrice (long idCurve,
									double matu,
									ARM_result& result);
extern long ARMLOCAL_DiscountYield (long idCurve,
									double matu,
									long meth,
									ARM_result& result);

extern long ARMLOCAL_ForwardPrice (long idCurve,
								   double matu1,
								   double matu2,
								   ARM_result& result);
extern long ARMLOCAL_ForwardYield (long idCurve,
								   double matu1,
								   double matu2,
								   long meth,
								   long adjDaycountId,
								   long decompFreqId,
								   long daycountId,
								   ARM_result& result);

extern long ARMLOCAL_zcspreaded (long zcSprId,
								 long zcInitId,
								 double date,
								 long MMFreq,
								 long SwapFreq,
								 bool ccyIsObject,
								 const CCString& ccyName,
								 ARM_result& result,
								 long objId = -1);

extern long ARMLOCAL_CreateZCSwapIntSmooth (double Date,
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

extern long ARMLOCAL_CreateZCSwapFutInt (double Date,
										 VECTOR<CCString>& matu,
										 VECTOR<double>& rate,
										 long MMVsFut,
										 long SwapVsFut,
										 long Raw,
										 long interp,
										 const CCString& Ccy,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_CreateZCSwapFutIntSmooth (double Date,
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

extern long ARMLOCAL_ZCINTSMOOTH (long inCvId,
								  double lambda,
								  long prec,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_ZCINTSMOOTH (const VECTOR<double>& matu,
								  const VECTOR<double>& rate,
								  double aDate,
								  long meth,
								  double lambda,
								  long prec,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_zcswapcubdiff (double date,
									VECTOR<CCString>& matu,
									VECTOR<double>& rate,
									long mmVsFut,
									long swapVsFut,
									long raw,
									long interp,
									const CCString& ccy,
									ARM_result& result,
									long objId = -1);

extern ARM_ZeroLInterpol* GetZCFromSummitNoETK(const CCString& index,
							                   const CCString& currency,
							                   const CCString& cvName,
							                   ARM_Date asof,
											   long interpId,
                                               ARM_result& result);

long ARMLOCAL_GetZCFromSummit (			const CCString& index,
										const CCString& currency,
										const CCString& cvName,
										double aSdate,
										long interpId,
										ARM_result& result,
										long objId = -1);

long ARMLOCAL_GetZCFromCalypso (const ARM_Date& AsOf,
							const std::string& index,
							const std::string& ccy,
							const std::string& term,
							const std::string& PricingEnv,
							const std::string& ForceCurveName,
							long interpId,
							const std::string& xmlFileName,
							long objId =-1); 

extern long ARMLOCAL_GetInitialCurveFromSummit (const CCString& index,
												const CCString& currency,
												const CCString& cvName,
												double aSdate,
												long adjOrNotId,
												VECTOR<CCString>* matu,
												VECTOR<double>* yield,
												ARM_result& result);

void  ARMLOCAL_GetInitialCurveFromCalypso (const ARM_Date& AsOf,
										 const std::string & pricingEnv,
										 const std::string & index,
										 const std::string & currency,
										 const std::string & term,
										 const std::string & forceCurveName,
										 const std::string & xmlFileName,
										 bool doAdj,
											std::vector<std::string >& matu,
											std::vector<double>& yield);

extern long ARMLOCAL_CreateZCFromSummit (const CCString& index,
										 const CCString& currency,
										 const CCString& cvName,
										 double aSdate,
										 long adjOrNotId,
										 const CCString& raw,
										 long swapFrqId,
										 ARM_result& result,
										 long objId = -1);

long ARMLOCAL_CreateZCSpreadedFromSummit(const CCString& index,
										 const CCString& currency,
										 const CCString& cvName,
										 double aSdate,
										 long adjOrNotId,
										 const CCString& raw,
										 long swapFrqId,
										 long mmFrqId,
										 long interpId,
										 long zcInitId,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_zcvsk (const VECTOR<double>& param,
							double date,
							ARM_result& result,
							long objId=-1);

extern long ARMLOCAL_zcsplicub (VECTOR<double>& matu,
								VECTOR<double>& rate,
								long meth,
								double date,
								long lastBucket,
								ARM_result& result,
								long objId = -1);

extern long ARMLOCAL_zcspl (const VECTOR<double>& param,
							double date,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_CreateTOYNYZCSwapInt (double Date,
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

extern long ARMLOCAL_BumpCurve (long cvId,
								VECTOR<CCString>& matu,
								VECTOR<double>& epsilon,
								long meth,
								ARM_result& result,
								long objId = -1);

extern long ARMLOCAL_BumpSpreadedCurve (long cvId,
										VECTOR<CCString>& matu,
										VECTOR<double>& epsilon,
										CCString& curveToBump,
										long meth,
										ARM_result& result,
										long objId = -1);

extern long ARMLOCAL_CreateZCCashInt (double date,
									  VECTOR<CCString>& matu,
									  VECTOR<double>& rate,
									  VECTOR<long>& bondsId,
									  VECTOR<double>& yields,
									  long MMVsFut,
									  const CCString& Ccy,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_GenerateBasisAdjCurve(long Crv1Id,
										   long BSCrv1Id,
										   long Crv2Id,
										   long BSCrv2Id,
										   long flagInputAsSprds,
										   long flagRetSprds,
										   VECTOR<CCString>& matu,
										   ARM_result& result,
										   long objId = -1);

ARM_ZeroLInterpol* obj_getZcFromSummit (const CCString& index,
										const CCString& currency,
										const CCString& cvName,
										double aSdate,
										long interpId,
										ARM_result& result);

ARM_ZeroLInterpol* obj_getZcFromCalypso(const ARM_Date& AsOf,
										const std::string & index,
										const std::string & ccy,
										const std::string & term,
										const std::string & PricingEnv,
										const std::string & ForceCurveName,
										long interpId,
										const std::string & xmlFileName) ;

extern long ARMLOCAL_GenForwardYield (long idCurve,
									  double matu1,
									  double matu2,
									  long isSwapRate,
									  long decompFreqId,
									  long daycountId,
									  ARM_result& result);

extern long ARMLOCAL_GetMaturitiesFromZC(ARM_result& result,
		    							 long objId);


extern long ARMLOCAL_OldZcCurve(const long& zcId, const double& asOf, ARM_result& result, long objId);


extern long ARMLOCAL_zcswapsplsum (double date,
							 VECTOR<CCString>& matu,
							 VECTOR<double>& rate,
							 long mmVsFut,
							 long swapVsFut,
							 long raw,
							 long interp,
							 const CCString& ccy,
							 ARM_result& result,
							 long objId = -1);


extern long ARMLOCAL_FixingSched(const double&		asOf,
								 VECTOR<string>		LiborKeys,
								 VECTOR<string>		FXKeys,
								 VECTOR<long>		LiborCurveId,
								 VECTOR<long>		FXCurveId,
								 ARM_result&		result,
								 long				objId=-1);


extern long ARMLOCAL_GetFixingSchedFromSummit(VECTOR<string>&		ListOfKeys,
											  const double&			AsOf,
											  const CCString&		Source,
											  long&					DateStrip,
											  ARM_result&			result,
											  long					objId=-1);



#endif	// ARM_LOCAL_ZCCURVE_H

// EOF %M%