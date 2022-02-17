#ifndef _CMSCTS_PRICING_H_
#define _CMSCTS_PRICING_H_

#include "CMSCTSProdStruct.h"

/* Function to price the underlying */
Err cmscts_calc_mkt_value(CMSCTS_DEAL sDeal, CMSCTS_MARKET sMarket,
                          CMSCTS_PRICING_PARAMS sPricingParams,
                          CMSCTS_CALIB sCalibration);

#endif