/* ===============================================================================

   FILENAME	:	swp_f_curve_fct.c
	
   PURPOSE  :   macros and functions to access curves information from the 
                curve object stored in the double linked list :
                SrtCurve = SrtCurveDesc

   =============================================================================== */

#include "swp_h_all.h"
#include "swp_h_curve_fct.h"
#include "srt_h_fwdobj.h"

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
                      Functions to get today == clcn_date 
   ------------------------------------------------------------------------ */

Err get_curve_clcndate(SrtCurvePtr crv, long *clcn_date)
{
Err     err  = NULL;

	if (!crv)
		return serror("Empty Curve in get_curve_clcndate");

	switch (crv->curve_type)
	{
		case YIELD_CURVE:
			*clcn_date = get_clcndate_from_yldcrv(crv);	
		break;
		case DVD_CURVE:
			*clcn_date = get_clcndate_from_dvdcrv(crv);	
		break;
		case REPO_CURVE:
			*clcn_date = get_clcndate_from_repocrv(crv);	
		break;
		case CMT_CURVE:
			*clcn_date = get_clcndate_from_cmtcrv(crv);	
		break;
		default:
			*clcn_date = 0;
		break;
	}

	return err;

} /* Err get_curve_clcndate(...) */

/* ------------------------------------------------------------------------ */

long get_clcndate_from_curve(SrtCurvePtr crv)
{
Err    err = NULL;
long   today;

	err = get_curve_clcndate(crv,&today);
	if (err)
		return 0;
	else
		return today;

} /* long get_clcndate_from_curve(...) */

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
     Functions to get the name of the yield curve asociated to a curve  
   ------------------------------------------------------------------------ */
Err get_curve_ycname(SrtCurvePtr crv, String *yc_name)
{
Err err  = NULL;  

	if (!crv)
		return serror("Empty Curve in get_curve_ycname");

	switch (crv->curve_type)
	{
		case DVD_CURVE:
			(*yc_name) = get_curve_name(crv);
		break;
		case REPO_CURVE:
			(*yc_name) = get_curve_name(crv);
		break;
		case CMT_CURVE:
			(*yc_name) = get_ycname_from_cmtcrv(crv);
		case YIELD_CURVE:
			(*yc_name) = get_curve_name(crv);
		break;
	}
	
	return NULL;

} /* END Err get_curve_ycname(...) */

/* ------------------------------------------------------------------------ */

char *get_ycname_from_curve(SrtCurvePtr crv)
{
String   yc_name = NULL;
Err      err;	

	err = get_curve_ycname(crv, &yc_name);
	if (err)
		return NULL;
	else
		return yc_name;

} /* END char *get_ycname_from_curve(...) */


/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
                      Functions to get the spot lag 
   ------------------------------------------------------------------------ */

Err get_curve_spotlag(SrtCurvePtr crv, int *spotlag)
{
Err          err  = NULL;
SrtCcyParam	*ccy_param;
String		 ccy_str;

	if (!crv)
		return serror("Empty Curve in get_curve_clcndate");


	switch (crv->curve_type)
	{
		case DVD_CURVE:
			*spotlag = 0;
		break;
		case REPO_CURVE:
			*spotlag = 0;
		break;
		case CMT_CURVE:
		case YIELD_CURVE:
			ccy_param = get_ccyparam_from_yldcrv(crv);
			if (!ccy_param)
			{
				ccy_str = get_curve_ccy(crv);
				err = swp_f_get_CcyParam_from_CcyStr(ccy_str, &ccy_param );
				if (err)
					return err;
			}
			*spotlag = ccy_param->spot_lag;
		break;
	}
	
	return NULL;

} /* END Err get_curve_spotlag(...) */

/* ------------------------------------------------------------------------------- */

int get_spotlag_from_curve(SrtCurvePtr crv)
{
int   spotlag;
Err    err      = NULL;

	err = get_curve_spotlag(crv, &spotlag);
	if (err)
		return 0;
	else
		return spotlag;

} /* long get_spotlag_from_curve(...) */

/* ------------------------------------------------------------------------------- */

