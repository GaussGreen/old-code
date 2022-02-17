#ifndef _CMSCTS_CALLER_H_
#define _CMSCTS_CALLER_H_

#include "CMSCTSProdStruct.h"

Err	CMSCTSAutocal(CMSCTS_MARKET			sMarket,
				  CMSCTS_DEAL			sDeal,
				  CMSCTS_CALIB			sCalibration,
				  CMSCTS_PRICING_PARAMS	sPricingParams,
				  CMSCTS_OUTPUTS		sOutputs);

void cmscts_market_set_default(CMSCTS_MARKET sMarket);
void cmscts_deal_set_default(CMSCTS_DEAL sDeal);
void cmscts_calibration_set_default(CMSCTS_CALIB sCalibration);
void cmscts_pricingparams_set_default(CMSCTS_PRICING_PARAMS sPricingParams);
void cmscts_outputs_set_default(CMSCTS_OUTPUTS sOutputs);

#endif