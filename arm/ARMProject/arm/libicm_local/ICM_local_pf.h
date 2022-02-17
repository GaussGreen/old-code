#ifndef ICM_LOCAL_PORTFOLIO_H
#define ICM_LOCAL_PORTFOLIO_H

#include "ICMKernel/util/icm_vector.h"

class ARM_result;


long ICMLOCAL_Portfolio (const VECTOR<long>& SecuritiesID,
								const long& CashFlowId,
								ARM_result& result,
								long objId = -1);


long ICMLOCAL_Leg_Basket (const VECTOR<long>& SecuritiesID,
								const long& CashFlowId,
								ARM_result& result,
								long objId = -1); 


long ICMLOCAL_Collateral (const vector<std::string>& Labels,
						  const ARM_Vector & Notionals,
						  // const vector<bool>& IsInDefault,
						  long objId = -1); 

long ICMLOCAL_VariableCollateral (const vector<std::string>& Labels,
						  const ICM_Vector<long>& NotionalsIds,
						  long objId = -1); 


long ICMLOCAL_CDO_SQUARE_GEN(const int& CdsId,
								const double& SubAmount,
								const int& PortfolioId,
								const double& Binary,
								const double& RcvFee,
								ARM_result& result,
								long	objId = -1);


#endif	

/*----End Of File ----*/
// EOF %M%