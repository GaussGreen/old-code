#ifndef ARM_OPTION_H
#define ARM_OPTION_H



#include "ARM_result.h"



extern long ARM_OPTION (long underId,
					    double maturityDate,
					    double strike,
					    long optionType,
					    long exerciseType,
						long strikeType,
						double FstXDate,
					    ARM_result& result,
					    long objId = -1);

extern long ARM_EXOPTION (long underId,
						  long optionType,
						  long styleId,
						  long kRefValId,
						  ARM_result& result,
						  long objId = -1);

extern long ARM_VolImp (long secId,
						long modId,
						double price,
						ARM_result& result);

extern long ARM_GetUnderPrice (long secId,
							   ARM_result& result);


extern long ARM_bsOption (double spot,
						  double strike,
						  double volatility,
						  double dividend,
						  double discountRate,
						  double maturity,
						  long CallPut,
						  ARM_result& result);

extern long ARM_bsDelta (double spot,
						 double strike,
						 double volatility,
						 double dividend,
						 double discountRate,
						 double maturity,
						 long CallPut,
						 ARM_result& result);

extern long ARM_bsVega (double spot,
						double strike,
						double volatility,
						double dividend,
						double discountRate,
						double maturity,
						long CallPut,
						ARM_result& result);

extern long ARM_bsTheta (double spot,
						 double strike,
						 double volatility,
						 double dividend,
						 double discountRate,
						 double maturity,
						 long CallPut,
						 ARM_result& result);

extern long ARM_bsGamma (double spot,
						 double strike,
						 double volatility,
						 double dividend,
						 double discountRate,
						 double maturity,
						 long CallPut,
						 ARM_result& result);

extern long ARM_ARM_OPTIONPORTFOLIO (long portfolioId,
									 long styleId,
									 long kRefValId,
									 long optionType,
									 ARM_result& result,
									 long objId = -1);
					    
#endif	// ARM_OPTION_H



// EOF %M%