#ifndef ARM_LOCAL_BOND_H
#define ARM_LOCAL_BOND_H



extern long ARMLOCAL_bond (double issueDate,
						   double maturityDate,
						   double firstCouponDate,
						   double couponRate,
						   double redemptionPrice,
						   long periodicity,
						   long dayCount,
						   long settleGap,
						   long couponDateFlag,
						   long ccy,
						   ARM_result& result,
						   long objId = -1);


extern long ARMLOCAL_RiskyBond (double issueDate,
						   double maturityDate,
						   double firstCouponDate,
						   double couponRate,
						   double redemptionPrice,
						   long periodicity,
						   long dayCount,
						   long settleGap,
						   long couponDateFlag,
						   long ccy,
						   double srepo,
						   double ssl,
						   double recoveryRate,
						   ARM_result& result,
						   long objId = -1);

extern long ARMLOCAL_RiskyBondWithCF (double asOfDate,
									  double redemptionPrice,
									  long periodicity,
									  long dayCount,
									  long settleGap,
									  long couponDateFlag,
									  long ccy,
									  VECTOR<double>& yearTerms, 
									  VECTOR<double>& cashFlows,
									  double srepo,
									  double ssl,
									  double recoveryRate,
									  ARM_result& result,
									  long objId = -1);
 

extern long ARMLOCAL_YTOPRICE (long bondId,
							   double settlement,
							   double yield,
							   ARM_result& result);

extern long ARMLOCAL_PTOYIELD (long bondId,
							   double settlement,
							   double price,
							   ARM_result& result);

extern long ARMLOCAL_BDFAPRICE (long bondId,
								double settlement,
								double actuPrice,
								double forwardDate,
								double repoRate,
								ARM_result& result);

extern long ARMLOCAL_BDREPORATE (long bondId,
								 double settlement,
								 double actuPrice,
								 double forwardDate,
								 double forwardPrice,
								 ARM_result& result);

extern long ARMLOCAL_bondTEC (double issueDate,
							  double maturityDate,
							  double firstCouponDate,
							  double couponRate,
							  double redemptionPrice,
							  long periodicity,
							  long dayCount,
							  long settleGap,
							  long couponDateFlag,
							  long ccyId,
							  double tec,
							  long pfTECId,
							  long ModtecId,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_YTODURATION (long bondId,
								  double settlement,
								  double actuRate,
								  long flagCpn,
								  ARM_result& result);

extern long ARMLOCAL_YTOCONVEXITY (long bondId,
								   double settlement,
								   double actuRate,
								   ARM_result& result);

#endif	// ARM_LOCAL_BOND_H
