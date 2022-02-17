#ifndef ICM_LOCAL_SUMMIT_H
#define ICM_LOCAL_SUMMIT_H

#include <ARM\libarm\ARM_result.h>
#include <icmkernel\util\icm_qmatrix.h>

/** 
long ICMLOCAL_GetModelFromSummit (const long DiscountCurveId,
									  const CCString& idSummit,
									  const CCString& type,
									  const CCString& CurveId,
									  const CCString& CorrCurveId,
									  ARM_result& result,
									  long objId = -1);
									  **/ 

void 
ICMLOCAL_GetBasketCorrelMkDataFromCalypso(const std::string& pricingEnv,const ARM_Date& date,
										  const std::string& forceCurveName,const std::string& xmlFilename,
										  std::vector<std::string>& matus,
										  std::vector<double>&tenors,
										  ICM_QMatrix<double>& correls) ;

long ICMLOCAL_CreateBasketCorrelMkDataFromCalypso(const std::string& pricingEnv,const ARM_Date& date,
												  const std::string& forceCurveName,const std::string Ccy ,const std::string& xmlFilename,
										          long IndexId ,long objId=-1);
#endif 
