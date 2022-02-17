#ifndef ICM_LOCAL_CASH_FLOW_H
#define ICM_LOCAL_CASH_FLOW_H

class ARM_result;


extern long ICMLOCAL_CashFlows (const VECTOR<CCString>& matrice,
								int nbrows,
								int nbcolumns,
								ARM_result& result,
								long objId = -1);

extern long ICMLOCAL_Parameters (const VECTOR<CCString>& matrice,
								int nbrows,
								int nbcolumns,
								ARM_result& result,
								long objId = -1);

#endif	

/*----End Of File ----*/
// EOF %M%