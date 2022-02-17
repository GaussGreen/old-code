/*******************************************************************************
*                                                                           
* PURPOSE      	: iniatializes the shift parameters for greek computation
*                                                                                                          
*******************************************************************************/

#include        <OPFNCTNS.H>
#include        <UTCONST.H>

/* global variable for options tools */

GlobVarOptionTool GVOPT = {VOL_ADD,VOL_SHIFT,SABRVOL_SHIFT,STRIKE_SHIFT,DELTA_SHIFT,VANNA_SHIFT, YEARS_IN_DAY,ALPHA_SHIFT,BETA_SHIFT,RHO_SHIFT,ZETA_SHIFT};

  
Err srt_f_defaultgreekshift(void)
{
	GVOPT.vol_add =			VOL_ADD;
	GVOPT.vol_shift =		VOL_SHIFT;
	GVOPT.sabrvol_shift =	SABRVOL_SHIFT;
	GVOPT.strike_shift =	STRIKE_SHIFT;
	GVOPT.delta_shift =		DELTA_SHIFT;
	GVOPT.theta_shift =		YEARS_IN_DAY;
	GVOPT.alpha_shift =		ALPHA_SHIFT;
	GVOPT.beta_shift =		BETA_SHIFT;
	GVOPT.rho_shift =		RHO_SHIFT;
	GVOPT.zeta_shift =		ZETA_SHIFT;
	GVOPT.vanna_shift =		VANNA_SHIFT;

	return NULL;
}

Err srt_f_initgreekshift(double vol_add, double vol_shift,double sabrvol_shift, double strike_shift,double delta_shift,
						 double theta_shift,double alpha_shift,double beta_shift,double rho_shift,double zeta_shift)
{
	GVOPT.vol_add =			vol_add;
	GVOPT.vol_shift =		vol_shift;
	GVOPT.sabrvol_shift =	sabrvol_shift;
	GVOPT.strike_shift =	strike_shift;
	GVOPT.delta_shift =		delta_shift;
	GVOPT.theta_shift =		theta_shift;
	GVOPT.alpha_shift =		alpha_shift;
	GVOPT.beta_shift =		beta_shift;
	GVOPT.rho_shift =		rho_shift;
	GVOPT.zeta_shift =		zeta_shift;
	GVOPT.vanna_shift =		delta_shift;
	return NULL;
}
Err srt_f_setsabrvolshift(double sabrvol_shift)
{
	GVOPT.sabrvol_shift = sabrvol_shift;

	return NULL;
}

double srt_f_getsabrvolshift(void)
{
	return(GVOPT.sabrvol_shift);
}

Err srt_f_setstrikeshift(double strike_shift)
{
	GVOPT.strike_shift = strike_shift;

	return NULL;
}

double srt_f_getstrikeshift(void)
{
	return(GVOPT.strike_shift);
}


Err srt_f_setdeltashift(double delta_shift)
{
	GVOPT.delta_shift = delta_shift;

	return NULL;
}

double srt_f_getdeltashift(void)
{
	return(GVOPT.delta_shift);
}


Err srt_f_setthetashift(double theta_shift)
{
	GVOPT.theta_shift = theta_shift;

	return NULL;
}

double srt_f_getthetashift(void)
{
	return(GVOPT.theta_shift);
}


Err srt_f_setvegashift(double vol_add)
{
	GVOPT.vol_add = vol_add;

	return NULL;
}

double srt_f_getvegashift(void)
{
	return(GVOPT.vol_add);
}

Err srt_f_setvannashift(double rates_add)
{
	GVOPT.vanna_shift = rates_add;

	return NULL;
}

double srt_f_getvannashift(void)
{
	return(GVOPT.vanna_shift);
}

Err srt_f_setalphashift(double alpha_shift)
{
	GVOPT.alpha_shift = alpha_shift;

	return NULL;
}

double srt_f_getalphashift(void)
{
	return(GVOPT.alpha_shift);
}

Err srt_f_setbetashift(double beta_shift)
{
	GVOPT.beta_shift = beta_shift;

	return NULL;
}

double srt_f_getbetashift(void)
{
	return(GVOPT.beta_shift);
}

Err srt_f_setrhoshift(double rho_shift)
{
	GVOPT.rho_shift = rho_shift;

	return NULL;
}

double srt_f_getrhoshift(void)
{
	return(GVOPT.rho_shift);
}

Err srt_f_setzetashift(double zeta_shift)
{
	GVOPT.zeta_shift = zeta_shift;

	return NULL;
}

double srt_f_getzetashift(void)
{
	return(GVOPT.zeta_shift);
}
