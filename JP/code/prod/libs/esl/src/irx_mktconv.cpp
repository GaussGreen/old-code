#include "irx/mktconv.h"
#include "irx/swap.h"
#include <irx/macros.h>

/*
 * Validates that a market conv object makes sense.
 */
int irxMarketConvValidate (const IrxTMarketConv *marketConv)
{
/* potentially down the line we could have market conventions which do
   not refer to interest rates */
    return irxSwapMarketConvValidate (marketConv);
}

/*
 * Defines the market conventions for a particular currency code.
 *
 * These are hard-coded within this function and is provided as a
 * convenience function.
 */
IrxTMarketConv* irxCcyMarketConv (const char* ccy)
{
    static char routine[] = "irxCcyMarketConv";

    static IrxTMarketConv USD_MKT_CONV =
    {
        {6, 'M', FALSE},      /* fixedIvl */
        {3, 'M', FALSE},      /* floatIvl */
        {6, 'M', FALSE},      /* cbsIvl */
        IRX_B30_360,          /* fixedDcc */
        IRX_ACT_360,          /* floatDcc */
        IRX_B30_360,          /* cbsDcc */
        IRX_ACT_360,          /* mmDcc */
        IRX_BAD_DAY_MODIFIED, /* paymentBdc */
        IRX_BAD_DAY_MODIFIED, /* accrualBdc */
        IRX_BAD_DAY_MODIFIED, /* resetBdc */
        IRX_BAD_DAY_MODIFIED, /* mmBdc */
        2                     /* daysToSpot */
    };

    static IrxTMarketConv EUR_MKT_CONV =
    {
        {1, 'Y', FALSE},      /* fixedIvl */
        {6, 'M', FALSE},      /* floatIvl */
        {1, 'Y', FALSE},      /* cbsIvl */
        IRX_B30_360,          /* fixedDcc */
        IRX_ACT_360,          /* floatDcc */
        IRX_B30_360,          /* cbsDcc */
        IRX_ACT_360,          /* mmDcc */
        IRX_BAD_DAY_MODIFIED, /* paymentBdc */
        IRX_BAD_DAY_MODIFIED, /* accrualBdc */
        IRX_BAD_DAY_MODIFIED, /* resetBdc */
        IRX_BAD_DAY_MODIFIED, /* mmBdc */
        2                     /* daysToSpot */
    };

    static IrxTMarketConv JPY_MKT_CONV =
    {
        {6, 'M', FALSE},      /* fixedIvl */
        {6, 'M', FALSE},      /* floatIvl */
        {6, 'M', FALSE},      /* cbsIvl */
        IRX_ACT_365F,         /* fixedDcc */
        IRX_ACT_360,          /* floatDcc */
        IRX_ACT_365F,         /* cbsDcc */
        IRX_ACT_360,          /* mmDcc */
        IRX_BAD_DAY_MODIFIED, /* paymentBdc */
        IRX_BAD_DAY_MODIFIED, /* accrualBdc */
        IRX_BAD_DAY_MODIFIED, /* resetBdc */
        IRX_BAD_DAY_MODIFIED, /* mmBdc */
        2                     /* daysToSpot */
    };
    
    REQUIRE (ccy != NULL);

    switch (*ccy)
    {
    case 'U':
        if (strcmp (ccy, "USD") == 0)
            return irxMarketConvCopy (&USD_MKT_CONV);
    case 'E':
        if (strcmp (ccy, "EUR") == 0)
            return irxMarketConvCopy (&EUR_MKT_CONV);
    case 'J':
        if (strcmp (ccy, "JPY") == 0)
            return irxMarketConvCopy (&JPY_MKT_CONV);
    default:
        break;
    }

    irxError ("%s: Unknown currency code %s\n", routine, ccy);

 RETURN: /* failure only */

    return NULL;
}
        
    

