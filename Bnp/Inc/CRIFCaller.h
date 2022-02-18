#ifndef _CRIF_CALLER_H_
#define _CRIF_CALLER_H_

#include "CRIFProdStruct.h"

Err CRIFAutocal(
    CMSCTS_MARKET       sMarket,
    CRIF_DEAL           sDeal,
    CMSCTS_CALIB        sCalibration,
    CRIF_PRICING_PARAMS sPricingParams,
    CRIF_OUTPUTS        sOutputs);

void crif_deal_set_default(CRIF_DEAL sDeal);
void crif_calibration_set_default(CMSCTS_CALIB sCalibration);
void crif_pricingparams_set_default(CRIF_PRICING_PARAMS sPricingParams);
void crif_outputs_set_default(CRIF_OUTPUTS sOutputs);

#endif