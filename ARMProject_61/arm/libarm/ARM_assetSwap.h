#ifndef ARM_ASSETSWAP_H
#define ARM_ASSETSWAP_H



#include "ARM_result.h"



#define BasisSwap_Method_Num		0
#define BasisSwap_Method_Alg		1
#define BasisSwap_Method_Num_Spr	2
#define BasisSwap_Method_Alg_Spr	3



extern long ARM_BasisSwap (double asOfDate,
						   double delivery,
						   double maturity,
						   double spread1,
						   const CCString& ccy1,
						   long liborType1Id,
						   long forwardCurve1Id,
						   long discountCurve1Id,
						   const CCString& ccy2,
						   long liborType2Id,
						   long forwardCurve2Id,
						   long discountCurve2Id,
						   long amortizationId,
						   long solve,
						   ARM_result& result);

extern double ASP_getStartXNL1 ();
extern double ASP_getEndXNL1 ();
extern double ASP_getStartXNL2 ();
extern double ASP_getEndXNL2 ();
extern double ASP_getPrice1 ();
extern double ASP_getPrice2 ();
extern double ASP_getMinimum ();
extern double ASP_getNbIter ();

extern long ARM_FRNPrice (double asOfDate,
					      double delivery,
					      double maturity,
						  const CCString& ccy1,
						  long liborType1Id,
					      long forwardCurve1Id,
						  long discountCurve1Id,
						  double facialMargin,
						  double valoMargin,
						  long frequencyId,
					      const CCString& ccy2,
					      long liborType2Id,
						  long forwardCurve2Id,	
						  long discountCurve2Id,
						  long amortizationId,
						  long frequencyId2,
					      double fixing,
						  double spread,
						  long solve,
						  ARM_result& result);


extern long ARM_ASWMargin (double bondMaturity,
						   double bondCoupon,
						   long bondFrequency,
						   long bondBase,
						   double bondPrice,
						   double bondRedemptionPrice,
						   double asOfDate,
						   double delivery,
						   long fixDecompFrequency,
						   long floatResetFreq,
						   long floatPayFreq,
						   const CCString& ccy1,
						   long liborType1Id,
						   long forwardCurve1Id,
						   long discountCurve1Id,
						   const CCString& ccy2,
						   long liborType2Id,
						   long forwardCurve2Id,
						   long discountCurve2Id,
						   long amortizationId,
						   long solve,
						   long viewFlag,
						   const CCString& id,
						   ARM_result& result);

extern long ARM_ASWPrice (double bondMaturity,
						  double bondCoupon,
						  long bondFrequency,
						  long bondBase,
						  double bondMargin,
						  double bondRedemptionPrice,
						  double asOfDate,
						  double delivery,
						  long fixDecompFrequency,
						  long floatResetFreq,
						  long floatPayFreq,
						  const CCString& ccy1,
						  long liborType1Id,
						  long forwardCurve1Id,
						  long discountCurve1Id,
						  const CCString& ccy2,
						  long liborType2Id,
						  long forwardCurve2Id,
						  long discountCurve2Id,
						  long amortizationId,
						  long solve,
						  double minValue,
						  double maxValue,
						  ARM_result& result);

extern long ARM_ARM_NextCpnDate (double asOfDate,
								 double maturity,
								 long frequencyId,
								 long ruleId,
								 const CCString& ccy,
								 ARM_result& result,
								 long intruleId = K_UNADJUSTED);

extern long ARM_ARM_PrevCpnDate (double asOfDate,
								 double maturity,
								 long frequencyId,
								 long ruleId,
								 const CCString& ccy,
								 ARM_result& result,
								 long intruleId = K_UNADJUSTED);

extern long ARM_CptBPV (double asOfDate,
						double delivery,
						double maturity,
						long zcId,
						long frequencyId,
						long dayCountId,
						const CCString& ccy,
						long amortizationId,
						ARM_result& result);

extern long ARM_LiborAssetSwapPriceAlg (long zcId,
										double startDate,
										double delivery,
										double endDate,
										double fixedRate,
										long fixDayCountId,
										long fixReceiveOrPayId, 
										long fixFrequencyId,
										long indexFrequencyId,
										long fixIntRuleId,
										double margin,
										const CCString& discountCcy,
										long amortizationId,
										double redemptionPrice,
										ARM_result& result);

extern long ARM_BasisSwapNum (double asOfDate,
						      double delivery,
						      double maturity,
						      double spread1,
						      const CCString& ccy1,
						      long liborType1Id,
						      long forwardCurve1Id,
						      long discountCurve1Id,
						      const CCString& ccy2,
						      long liborType2Id,
						      long forwardCurve2Id,
						      long discountCurve2Id,
							  long amortizationId,
						      ARM_result& result);

/*
extern long ARM_LiborAssetSwapPrice (long modelId, double startDate, double endDate,
							  double fixedRate, long fixReceiveOrPayId, 
							  long fixDayCountId, long fixFrequencyId,
							  long fixDecompFrequencyId, long fixPayTimingId, long fixIntRuleId,
							  long liborTypeId, double spread, long floatResetFreqId,
							  long floatPayFreqId, long assetGap, long vFlag,
							  double margin, const CCString& discountCcy,
							  double redemptionPrice, double supplFee,
							  long solve, const CCString& id,
							  double minValue, double maxValue,
							  ARM_result& result);
*/

extern long ARM_LiborAssetSwapPriceNew (long ycModId,
							  long forwardCurve1Id,
							  long fixedLeg1Id,
							  long fixedLeg2Id,
							  long fixedLeg3Id,
							  double asOfDate,
							  double delivery,
							  double bondMaturity,
							  double bondCoupon,
							  long fixReceiveOrPayId,
							  long floatPayFreq, long fixDayCountId, long fixFrequencyId,
							  long fixIntRuleId, double bondMargin, long ccy1Id,
							  long amortizationId, double bondRedemptionPrice,
							  double minValue, double maxValue, ARM_result& result);


extern long ARM_FRNPriceNew (double asOfDate,
					      double delivery,
					      double maturity,
						  const CCString& ccy1,
						  long liborType1Id,
					      long forwardCurve1Id,
						  long discountCurve1Id,
						  double facialMargin,
						  double valoMargin,
						  long frequencyId,
					      const CCString& ccy2,
					      long liborType2Id,
						  long forwardCurve2Id,	
						  long discountCurve2Id,
						  long amortizationId,
						  long frequencyId2,
					      double fixing,
						  double spread,
						  long solve,
						  ARM_result& result);

extern long ARM_FRNMarginNew (double asOfDate,
						   double delivery,
						   double maturity,
						   const CCString& ccy1,
						   long liborType1Id,
						   long forwardCurve1Id,
						   long discountCurve1Id,
						   double facialMargin,
						   double price,
						   long frequencyId,
						   const CCString& ccy2,
						   long liborType2Id,
						   long forwardCurve2Id,
						   long discountCurve2Id,
						   long amortizationId,
						   long frequencyId2,
						   double fixing,
						   double spread,
						   long solve,
						   ARM_result& result);

extern long ARM_ASWPriceNew (double bondMaturity,
						  double bondCoupon,
						  long bondFrequency,
						  long bondBase,
						  double bondMargin,
						  double bondRedemptionPrice,
						  double asOfDate,
						  double delivery,
						  long fixDecompFrequency,
						  long floatResetFreq,
						  long floatPayFreq,
						  const CCString& ccy1,
						  long liborType1Id,
						  long forwardCurve1Id,
						  long discountCurve1Id,
						  const CCString& ccy2,
						  long liborType2Id,
						  long forwardCurve2Id,
						  long discountCurve2Id,
						  long amortizationId,
						  long solve,
						  double minValue,
						  double maxValue,
						  ARM_result& result);

extern long ARM_ASWMarginNew (double bondMaturity,
					double bondCoupon,
					long bondFrequency,
					long bondBase,
					double bondPrice,
					double bondRedemptionPrice,
					double asOfDate,
					double delivery,
					long fixDecompFrequency,
					long floatResetFreq,
					long floatPayFreq,
					const CCString& ccy1,
					long liborType1Id,
					long forwardCurve1Id,
					long discountCurve1Id,
					const CCString& ccy2,
					long liborType2Id,
					long forwardCurve2Id,
					long discountCurve2Id,
					long amortizationId,
					long solve,
					double minValue,
					double maxValue,
					ARM_result& result);

extern long ARM_BasisSwapNew (double asOfDate,
					double delivery,
					double maturity,
					double margin1,
					long ccy1Id,
					long liborType1Id,
					long forwardCurve1Id,
					long discountCurve1Id,
					long ccy2Id,
					long liborType2Id,
					long forwardCurve2Id,
					long discountCurve2Id,
					long amortizationId,
					long solve,
					const CCString & ccy1,
					const CCString & ccy2,
					ARM_result& result);



long ARM_LiborAssetSwapPriceAlgNew (long modelId,
										 long zcId,
										 long fixedLeg1Id,
										 long fixedLeg2Id,
										 long fixedLeg3Id,
										 double startDate,
										 double delivery,
										 double endDate,
										 long indexFrequencyId,
										 double margin,
										 long discountCcyId,
										 long amortizationId,
										 double redemptionPrice,
										 ARM_result& result);


long ARM_LiborAssetSwapMarginNumNew (long ycModId,
										  long forwardCurve1Id,
										  long fixedLeg1Id,
										  long fixedLeg2Id,
										  long fixedLeg3Id,
										  double asOfDate,
										  double delivery,
										  long bondMaturity,
										  long floatPayFreq,
										  double bondPrice,
										  const CCString & ccy1,
										  long amortizationId,
										  double bondRedemptionPrice,
										  double minValue,
										  double maxValue,
										  ARM_result& result);

long ARM_LiborAssetSwapMarginAlgNew (long ycmodId,
										  long zcId,
										  long fixedLeg1Id,
										  long fixedLeg2Id,
										  long fixedLeg3Id,
										  double startDate,
										  double delivery,
										  double endDate,
										  double fixedRate,
										  long fixDayCountId,
										  long fixReceiveOrPayId,
										  long fixFrequencyId,
										  long indexFrequencyId,
										  long fixIntRuleId,
										  double price,
										  long discountCcyId,
										  long amortizationId,
										  double redemptionPrice,
										  ARM_result& result);

extern long ARM_CptBPVNew (double asOfDate,
						double delivery,
						double maturity,
						long zcId,
						long frequencyId,
						long dayCountId,
						long ccyId,
						long amortizationId,
						ARM_result& result);

extern long ARM_BasisSwapNumNew (double asOfDate,
						      double delivery,
						      double maturity,
						      double spread1,
						      const CCString& ccy1,
						      long liborType1Id,
						      long forwardCurve1Id,
						      long discountCurve1Id,
						      const CCString& ccy2,
						      long liborType2Id,
						      long forwardCurve2Id,
						      long discountCurve2Id,
							  long amortizationId,
						      ARM_result& result);

extern long ARM_BasisSwapAlgNew (double asOfDate,
					   double delivery,
					   double maturity,
					   double margin1,
					   const CCString& ccy1,
					   long liborType1Id,
					   long forwardCurve1Id,
					   long discountCurve1Id,
					   const CCString& ccy2,
					   long liborType2Id,
				 	   long forwardCurve2Id,
					   long discountCurve2Id,
					   long amortizationId,
					   long mode,
					   ARM_result& result);


extern long ARM_DiscountPriceRefvalue(long zcId, long refvalId, long discountCcyId,
									double startDate, double endDate, ARM_result& result);;
						

long Create1IdCCY(const CCString & ccy1,long & ccy1Id,long CurveId);
long Create2IdCCY(const CCString & ccy1,const CCString & ccy2,long & ccy1Id, long & ccy2Id ,long CurveId ,bool & oneCurrency );

#endif	// ARM_ASSETSWAP_H

// EOF %M%