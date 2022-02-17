#ifndef ICM_LOCAL_CAIR_H
#define ICM_LOCAL_CAIR_H

#include <vector>
#include <ARM\libarm\ARM_result.h>
//#include "CCstring.h"
using namespace std;


extern	long ICMLOCAL_SetCreditManager_MarketData(
							long			CreditManagerId,
							const double&	AsOfDate,
							CCString		CurrencyLabel,
							int				IR_SummitCurveId,
							long			MarketDate_ParametersId,
							long			Imposed_IR_ValuesId,
							ARM_result&		result);


extern	long ICMLOCAL_SetCreditManager_CreditData(
							long			CreditManagerId,
							VECTOR<CCString>	labels,
							long				CD_DescriptionId,
							vector<double>		MatIssuersSpread,
							VECTOR<CCString>	maturities,
							long				CD_ParametersId,
//							const VECTOR<long>& DefCurvesID,
							CCString			HedgesCDSMaturity,
							long				CD_CDOSquare_ParametersId,
							long				CD_CDOSquare_DataId,
							ARM_result&			result);


extern long ICMLOCAL_SetCreditManager_CreditModel(
							long			CreditManagerId,
							long				CM_ParameterId,
							double				correlation_value,
							VECTOR<double>		beta_vector,
							VECTOR<double>		base_correlation_strikes,
							VECTOR<double>		base_correlation_values,
							long				CorrelationMatrixId,
							long				CorrelationId,
							long				FactorLoading_ParameterId,
							ARM_result&			result);


extern long ICMLOCAL_SetCreditManager_CreditProduct(
							long			CreditManagerId,
							long				CP_DefaultLegId,
							long				CP_PremiumLegId,
							long				CP_PricingParametersId,
							ARM_result&			result);


extern long ICMLOCAL_GetCreditManager_DataFromLabel(
							long			CreditManagerId,
							CCString			DataLabel,
							ARM_result&			result);

extern long ICMLOCAL_GetCreditManager_DataMatrixFromLabel(
							long				CreditManagerId,
							CCString			DataLabel,
							VECTOR<double*>&	OutputMatrix,
							VECTOR<CCString>&	OutputLabels,
							int&				OutputNbRows,
							int&				OutputNbCols,
							ARM_result&			result);

extern long ICMLOCAL_SetCreditManager_CreditCalibrator(
							long			CreditManagerId,
							long				CC_DescriptionId,
							VECTOR<CCString>	maturities,
							ARM_result&			result);


extern long ICMLOCAL_SetCreditManager_CorrelationCalibrator(
							long				CreditManagerId,
							long				CC_DescriptionId,
							long				CC_ParametersId,
							VECTOR<long>		TranchesId,
							ARM_result&			result);


extern long ICMLOCAL_CreditManager(ARM_result&	result,
							long objId	=	-1);

#endif
