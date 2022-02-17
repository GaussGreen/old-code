#ifndef _CRIF_H_
#define _CRIF_H_
// ------------------------------------------------------------------------------------------------------------------ //
//
// CallableRatchetInverseFloater.h
//
//

#include "CRIFProdStruct.h"
#include "CRIFUtil.h"

/* Function to price the underlying */
Err	crif_calc_mkt_value(	CRIF_DEAL			sDeal,
							CMSCTS_MARKET			sMarket,
							CRIF_PRICING_PARAMS	sPricingParams,
							CMSCTS_CALIB			sCalibration);

Err	crif_calc_cv_value(		CMSCTS_MARKET				sMarket,
							CRIF_DEAL				sDeal,
							CRIF_SIMULARG			sSimulArg,
							CRIF_PRICING_PARAMS		sPricingParams,
							CRIF_OUTPUTS			sOutputs);


#endif 