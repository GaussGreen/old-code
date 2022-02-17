#ifndef ICM_LOCAL_MOD_H
#define ICM_LOCAL_MOD_H

class ARM_Date;

#include <ARM\libarm\ARM_result.h>



extern long ICMLOCAL_ModelMultiCurves ( int NbDefCurves, 
									VECTOR<long> DefCurvesID,
									long DiscCurveId,
									VECTOR<double> LossRates ,
									long CorrelationId,	
									long idvolcurve,
									bool  CloneOrNot,
									long idCpnInfCurve,
									long idCpnIRCurve,
									ARM_result& result, 
									long objId = -1);


extern long ICMLOCAL_ModelMultiCurves ( int NbDefCurves, 
									VECTOR<long> DefCurvesID,
									long DiscCurveId,
									VECTOR<double> LossRates ,
									long CorrelationId,
									long marketDataMngId,
									long idvolcurve,
									bool  CloneOrNot,
									long objId = -1);


/*extern long ICMLOCAL_SetPropMixCopule (long ModelId, 
									   double PropIndep,
									   double PropFullCorrel, 
									   ARM_result& result);
*/
extern long ICMLOCAL_DefProbModel ( long idDefProb, 
										long idIRcurve,
										long idvolcurve,
										ARM_result& result, 
										long objId = -1);

extern long ICMLOCAL_MetaModel (int Nbmodel, 
								 const VECTOR<long>& ModelID,
								 VECTOR<int>& PricerType,
								 ARM_result& result, 
								 long objId = -1);


extern long ICMLOCAL_MarketDataMng (const	VECTOR<long>& Data,
									long objId = -1);

extern long ICMLOCAL_MarkerDataMng (const	VECTOR<long>& Data,
									long objId = -1);

extern long ICMLOCAL_Customized_Credit_MultiCurves(
									long				DiscountCurveId,
									VECTOR<CCString>	labels,
									vector<double>		spreads,
									VECTOR<CCString>	maturities,								
									long				Data_DescriptionId,
									long				Market_ParametersId,
									long				CDO_Square_ParametersId,
									long				CDO_Square_DataId,
									long				CorrelationId,
									long				ModelMultiCurvesId,
									ARM_result& result, 
									long objId = -1);

extern long ICMLOCAL_SetVolCurve(long ModelId,
								   long idvolcurve,
								   ARM_result& result);
#endif	// ARM_LOCAL_MOD_H

// EOF %M%

