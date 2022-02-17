#ifndef ICM_LOCAL_PRICER_H
#define ICM_LOCAL_PRICER_H

#include <ARM\libarm\ARM_result.h>


long ICMLOCAL_Pricer (const ARM_Date*asof,	/** optional **/ 
					  long idSecurity, 
					  long idModel,
					  int PricerType,
					  int nbpaths,
					  long idParameters,
					  // double AsOfdate,
					  // long idMktDataManager,
					  // ARM_result& result, 
					  long objId = -1);


extern long ICMLOCAL_CacheOption (long pricerId, 
								  int Option, 
								  ARM_result& result);

/** 
extern double ICMLOCAL_Debug_Function (long pricerId,
									   double Double,
									   VECTOR<double>& Data,
									   ARM_result& result);
									   **/ 

long ICMLOCAL_PricerDefaultCdsNew (const ARM_Date*  asof /** optionnal */ ,
								   long idSecurity, 
								long idModel,
								// ARM_result& result, 
								long objId = -1);

 

extern long ICMLOCAL_SetVolatility(long pricerId,
								   long idvolcurve,
								   ARM_result& result);

extern long ICMLOCAL_GetPricer_DataFromLabel(
							long				PricerId,
							CCString			DataLabel,
							ARM_result&			result);



extern long ICMLOCAL_GenerateImpliedCurve (long pricerId, ARM_result& result,long objId=-1);

#endif	// ARM_LOCAL_MOD_H

// EOF %M%
