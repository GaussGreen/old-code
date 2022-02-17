/* -------------------------------------------------------------------------

   FILE NAME	: swp_f_external_fct.c
	
   PURPOSE		: Sets a global variable that stores all the external
                  functions that have to be set before calling any SRT
				  financial function
				  This include:
					- DiscFunc (df_curve_id/ccy, start, end, *df)
					- SpreadFunc (ref_rate_name, start, end, *spread)
					- GetVolFunc ( vol_curve_id, start, end, strike, *vol )
					- FixingFunc ()
					- ...
				  Gives the functions to communicate with this structure
   ------------------------------------------------------------------------- */

#include "swp_h_all.h"
#include "swp_h_external_fct.h"

/* -----------------------------------------------------------------------------
   This is the very important static declaration: this is the global structure
   to access these external functions
   ----------------------------------------------------------------------------- */

static SrtExternalFunctions  _external_functions ;

/* static void *_external_market; */


/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the External Functions to NULL : this has to be called when openeing SRT 
   ----------------------------------------------------------------------------- */

Err swp_f_ResetExternalFunctions ( )
{
	_external_functions.disc_func = NULL;
	
	_external_functions.spread_func = NULL;

	_external_functions.rate_info_func = NULL;
	
	_external_functions.vol_func = NULL;

	_external_functions.sabrvol_func = NULL;

	_external_functions.fixing_func = NULL;

	return NULL;
}

/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the Official DiscountFactor Function 
   ----------------------------------------------------------------------------- */

Err swp_f_SetDiscFunc ( DiscFuncType disc_func)
{
	_external_functions.disc_func = disc_func ;
	
	return NULL;
}

/* -----------------------------------------------------------------------------
   Gets the Official DiscountFactor Function 
   ----------------------------------------------------------------------------- */

Err swp_f_GetDiscFunc ( DiscFuncType *disc_func)
{
	*disc_func = _external_functions.disc_func;
	
	if ( !(*disc_func) )
		return serror ("No Discount Function available");

	return NULL;
}

/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the Official Spread Function 
   ----------------------------------------------------------------------------- */

Err swp_f_SetSpreadFunc ( SpreadFuncType spread_func)
{
	_external_functions.spread_func = spread_func ;
	
	return NULL;
}

/* -----------------------------------------------------------------------------
   Gets the Official Spread Function 
   ----------------------------------------------------------------------------- */

Err swp_f_GetSpreadFunc ( SpreadFuncType *spread_func)
{
	*spread_func = _external_functions.spread_func;

	if (!(*spread_func))
		return serror ("No Spread Function available");

	return NULL;
}


/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the Official Rate Information Function used in the pricing 
   ----------------------------------------------------------------------------- */
Err swp_f_SetRateInfoFunc ( RateInfoFuncType rate_info_func)
{
	_external_functions.rate_info_func = rate_info_func ;
	
	return NULL;
}

Err swp_f_GetRateInfoFunc ( RateInfoFuncType *rate_info_func)
{
	*rate_info_func = _external_functions.rate_info_func;

	if (!(*rate_info_func))
		return serror ("No Rate Information Function available");

	return NULL;
}

/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the Official Volatility Interpolation Function 
   ----------------------------------------------------------------------------- */

Err swp_f_SetVolFunc ( VolFuncType vol_func)
{
	_external_functions.vol_func = vol_func;

	return NULL;
}

/* -----------------------------------------------------------------------------
   Gets the Official Volatility Interpolation Function 
   ----------------------------------------------------------------------------- */

Err swp_f_GetVolFunc ( VolFuncType *vol_func)
{
	*vol_func = _external_functions.vol_func;

	if (!(*vol_func))
		return serror ("No Volatility Interpolation Function available");

	return NULL;
}


/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the Official SABRVolatility Function 
   ----------------------------------------------------------------------------- */

Err swp_f_SetSABRVolFunc ( SABRVolFuncType sabrvol_func)
{
	_external_functions.sabrvol_func = sabrvol_func;

	return NULL;
}

/* -----------------------------------------------------------------------------
   Gets the Official Volatility Interpolation Function 
   ----------------------------------------------------------------------------- */

Err swp_f_GetSABRVolFunc ( SABRVolFuncType *sabrvol_func)
{
	*sabrvol_func = _external_functions.sabrvol_func;

	if (!(*sabrvol_func))
		return serror ("No SABRVolatility Function available");

	return NULL;
}

/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   Sets the Official Historical Fixing  Function used in the pricing 
   ----------------------------------------------------------------------------- */

Err swp_f_SetFixingFunc ( FixingFuncType fixing_func)
{
	_external_functions.fixing_func = fixing_func;

	return NULL;
}

/* -----------------------------------------------------------------------------
   Gets the Official Historical Fixing Function used in the pricing 
   ----------------------------------------------------------------------------- */

Err swp_f_GetFixingFunc ( FixingFuncType *fixing_func)
{
	*fixing_func = _external_functions.fixing_func;

	if (!(*fixing_func))
		return serror ("No Historical Fixing Function available");

	return NULL;
}

/* ----------------------------------------------------------------------------- */
