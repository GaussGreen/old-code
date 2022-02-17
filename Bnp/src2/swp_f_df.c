/* --------------------------------------------------------------------------- 
   FILENAME:     swp_f_df.c
   
   PURPOSE		: compute a discount factor via a pointer to a function DiscFunc
                  passed to the SrtExternalFunctions static 
				  ( see srt_f_external_fct.c )
				  Also available, a function to compute a CashFwd 
   --------------------------------------------------------------------------- */
#include "swp_h_all.h"
#include "swp_h_external_fct.h"
#include "swp_h_df.h"
#include "math.h"

/* --------------------------------------------------------------------------- */

/* Note that the crv_name could also be the Currency (3 characters long) */

double swp_f_disc_or_zrate(
			Ddate 		  dstart ,
			Ddate 		  dend,
			char         *crv_name, 
			int 		  rate_or_disc)   
{
Err                  err = NULL;
double               zr;
double               df;
Err                  (*df_func)(char *crvname, double start, double end, double *df);
	
/* Gets the discount factor function attached to the library */
	err = swp_f_GetDiscFunc (&df_func);
	if (err) return SRT_SPREAD_ERROR;

/* Computes the discount factor */
	if (dend == dstart)
	{
		df = 1.0;
	}
	else
	if (dend > dstart)
	{
		err = df_func( crv_name, dstart, dend, &df);
	}
	else
	{
		err = df_func( crv_name, dend, dstart, &df);
		df = 1.0 / df;
	}
	
	if (err) return SRT_DF_ERROR;
	
/* Return the rate or the discount factor */
	if (rate_or_disc == 1)
	{
		return df;
	}
	else
	if (dstart == dend)
	{
		zr = 1.0;
	}
	else
	{
		zr = log (df) / (dstart-dend) / YEARS_IN_DAY;
	}
	return zr;

}/* END disc_crv(...) */ 

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
   This is a CASH fwd rate ( no reference rate)
   ----------------------------------------------------------------------- */
double 	swp_f_fwdcash( Ddate start, Ddate end, BasisCode b, char *crvname)
{
double               df;
double               fra;
double               cvg;

 /* If end = start, increment end by one day to calc IFR */
	if ( end == start) 
		end += 1.0;
	
/* Computes the fra (cash) through discount factors */
	df = swp_f_df( start, end, crvname);
	if (df == SRT_DF_ERROR)
		return SRT_DF_ERROR;

/* Computes the coverage */
	cvg = coverage( DTOL( start), DTOL( end), b);

/* Computes the Fra */
	fra = ( 1.0 / df ) - 1.0;
	fra = fra / cvg;

	return fra;

} /* END double 	swp_f_fwdcash(...) */

/*-----------------------------------------------------------------------------------
This is a proper cash rate taking into account the fixing date, the start and end dates
 ----------------------------------------------------------------------- */
double 	swp_f_fwdcashrate( Ddate fixing, Ddate start, Ddate end, BasisCode b, char *crvname)
{
double               df1, df2;
double               fra;
double               cvg;

 /* If end = start, increment end by one day to calc IFR */
	if ( end == start) 
		end += 1.0;
	
/* Computes the fra (cash) through discount factors */
	df1 = swp_f_df(fixing, start, crvname);
	df2 = swp_f_df(fixing, end, crvname);
	if ((df1 == SRT_DF_ERROR)||(df2 == SRT_DF_ERROR))
		return SRT_DF_ERROR;

/* Computes the coverage */
	cvg = coverage( DTOL( start), DTOL( end), b);

/* Computes the Fra */
	fra = ( df1 / df2 ) - 1.0;
	fra = fra / cvg;

	return fra;

} /* END double 	swp_f_fwdcashrate(...) */

/*-----------------------------------------------------------------------------------
This is a swap cash rate: It uses a Theo End date.
 ----------------------------------------------------------------------- */
double 	swp_f_swapcashrate( long start, long theoend, SrtBasisCode Basis, SrtCompounding Comp, char  *ycName, char  *refRateCode)
{
double               df1, df2;
long today;
String basisStr;
String compStr;
SrtCrvPtr yldcrv;
double level;
double dswap;
long end;


yldcrv = lookup_curve(ycName);
today = get_clcndate_from_yldcrv(yldcrv);
translate_basis(&basisStr, Basis);
translate_compounding(&compStr, Comp);

end = add_unit(theoend,0,SRT_BDAY, MODIFIED_SUCCEEDING);

/*Computes the discount factors for the numerator*/
	df1 = swp_f_df(today,start, ycName);
	df2 = swp_f_df(today, end, ycName);

	if ((df1 == SRT_DF_ERROR)||(df2 == SRT_DF_ERROR))
		return SRT_DF_ERROR;

/* Computes the Level */
swp_f_LevelPayment(start,theoend,compStr,basisStr,ycName,refRateCode,&level);
	
/* Computes the Swaprate */
dswap = (df1-df2)/level;
	
	return dswap;

} /* END double 	swp_f_swapcashrate(...) */

/* ----------------------------------------------------------------------- */

