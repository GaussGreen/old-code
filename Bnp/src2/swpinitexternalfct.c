/* ===================================================================================

   FILENAME:   SwpInitExternalFct.c

   PURPOSE:    When using the library in an external environment, this function
                           sets up the official functions required from the outside by the library,
                           such as;
                                 - discount factor
                                         - volatility interpolation
                                         - spread calculation
                                         - reference rate details
                                         - historical fixings
                                         - ...
                           to provide OFFICIAL values for the pricing of vanilla (swaps and
                           simple options)

  =================================================================================== */

#include "SwpAccess.h"

Err SwpInitExternalFct(
    DiscFuncType     disc_func,
    SpreadFuncType   spread_func,
    RateInfoFuncType rate_info_func,
    VolFuncType      vol_func,
    SABRVolFuncType  sabrvol_func,
    FixingFuncType   fixing_func)
{
    Err err = NULL;

    err = swp_f_SetDiscFunc(disc_func);

    err = swp_f_SetSpreadFunc(spread_func);

    err = swp_f_SetRateInfoFunc(rate_info_func);

    err = swp_f_SetVolFunc(vol_func);

    err = swp_f_SetSABRVolFunc(sabrvol_func);

    err = swp_f_SetFixingFunc(fixing_func);

    return NULL;

} /* END Err SwpInitExternalFct(...) */

/* ---------------------------------------------------------------------------------- */
