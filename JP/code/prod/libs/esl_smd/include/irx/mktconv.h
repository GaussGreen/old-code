#ifndef IRX_MKTCONV_H
#define IRX_MKTCONV_H

#include "irxflow.h"

#ifdef __cplusplus
extern "C" 
{
#endif

/**
 * Validates that a market conv object makes sense.
 */
int irxMarketConvValidate (const IrxTMarketConv *marketConv);

/**
 * Defines the market conventions for a particular currency code.
 *
 * These are hard-coded within this function and is provided as a
 * convenience function.
 */
IrxTMarketConv* irxCcyMarketConv (const char* ccy);

#ifdef __cplusplus
}
#endif

#endif
