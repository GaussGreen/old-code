/* -------------------------------------------------------------------------

   FILE NAME	: swp_f_spread.c
	
   PURPOSE		: compute spreads via a pointer to a function SpreadFunc
                  passed to the SrtExternalFunctions static 
				  ( see srt_f_external_fct.c )
   ------------------------------------------------------------------------- */

#include "swp_h_all.h"
#include "swp_h_external_fct.h"
#include "swp_h_df.h"
#include "swp_h_spread.h"


/* --------------------------------------------------------------------------- 
	Calculate spreads , using the SpreadFunc that has been attached to the 
	library (see srt_f_external_fct.c)
  --------------------------------------------------------------------------- */

double 	swp_f_spread( Ddate start, Ddate end, String ref_rate )
{
double     spread;
Err        err;
Err        (*spread_func)(char *name, double start, double end, double *spread);

/* If the reference rate is cash: no need to do anything : spread is 0.0 */
	if (!strcmp(ref_rate,"CASH"))
		return 0.0;
		
/* Gets spread function attached to the library */
	err = swp_f_GetSpreadFunc(&spread_func);
	if (err)
		return SRT_SPREAD_ERROR;

/* Computes the spread for the */
	err = spread_func(ref_rate, start, end, &spread);
	if (err)
		return SRT_SPREAD_ERROR;
	
	return spread;

} /* END srt_f_spread(...) */

/* ----------------------------------------------------------------------------- */

/* When dealing with spreads in swaps or fra's, it is useful to get floating 
   compd and basis for the reference rate */ 
        
Err swp_f_get_ref_rate_details(
			String          ref_rate_name,
			SrtBasisCode    *float_basis,
			SrtCompounding  *float_compounding)
{
Err                 err;
RateInfoFuncType    rate_info_func;
char                family_code[32];
char                ccy[32];
char                tenor[32];
char                basis[32];
int                 frequency;
char                calc_method[32];
int                 number_of_period;

/* If there is no reference rate mentionned, do not do anything */
	if (!ref_rate_name)
		return "No reference rate mentionned in srt_f_get_ref_rate_details";

/* If the ref rate is CASH, make it quick */
	if (!strcmp (ref_rate_name,"CASH") ) 
	{
		*float_basis  = BASIS_ACT_360;
		*float_compounding = SRT_QUARTERLY;
		return NULL;
	}
	
/* Get the official reference rate properties */
	err = swp_f_GetRateInfoFunc ( &rate_info_func);
	if (err) 
		return err;
	err =   rate_info_func ( ref_rate_name, family_code, ccy, tenor, basis,  
				&frequency, calc_method, &number_of_period);
	if (err) 
		return err;

/* Turns the compounding and basis strings into types */
	err = interp_basis(basis, float_basis);
	if (err)
		return err;
	*float_compounding = (SrtCompounding)frequency;

	return NULL;

}/* END srt_f_get_ref_rate_details

/* ----------------------------------------------------------------------------- */

Err srt_f_get_spot_lag_from_refrate(
			String          ref_rate_name,
			int             *spot_lag)
{
Err                 err;
RateInfoFuncType    rate_info_func;
char                family_code[32];
char                ccy[32];
char                tenor[32];
char                basis[32];
int                 frequency;
char                calc_method[32];
int                 number_of_period;
SrtCcyParam			*ccy_param;

/* If there is no reference rate mentionned, do not do anything */
	if (!ref_rate_name)
		return "No reference rate mentionned in srt_f_get_spot_lag_from_refrate";

/* Get the official reference rate properties */
	err = swp_f_GetRateInfoFunc ( &rate_info_func);
	if (err) 
		return err;
	err =   rate_info_func ( ref_rate_name, family_code, ccy, tenor, basis,  
				&frequency, calc_method, &number_of_period);
	if (err) 
		return err;

   err = swp_f_get_CcyParam_from_CcyStr(ccy, &ccy_param);
   if (err) return err;

   *spot_lag = ccy_param->spot_lag;

	return NULL;

}/* END srt_f_get_spot_lag_from_refrate

/* ----------------------------------------------------------------------------- */

/* This function returns the short forward rate, adjusted by the cash spread  */ 

double 	swp_f_fra( Ddate start, Ddate end, BasisCode b, char *crvname, String ref_rate )
{
double  fwd_cash;
double  spread;
double  fra;
	
/* Computes the Cash forward rate (Fra) */
	fwd_cash = swp_f_fwdcash( start, end, b, crvname);

/* Computes the spread for the relevant period */
	spread = swp_f_spread(start, end, ref_rate );
	
/* The Fra is the sum of the cash forwar and the spread */
	fra = fwd_cash + spread;

	return fra;

} /* END double swp_f_fra_from_crv(...) */

/* ----------------------------------------------------------------------- */


/* ======================================================================= */


