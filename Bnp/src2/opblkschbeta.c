

/******************************************************************************/
/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optblkschbeta(...)                                   
*                                                                       
* PURPOSE      	: Premium of a European call option for a single risk factor 
*				  using dF=a*F^{beta}dW    
*                                                                      
* DESCRIPTION  	: Black-Scholes: calculates the option premium in which the
*		          underlying can vary between normal to lognormal by varying 0<=beta<1                       
*                                                                   
* CALLS		: I_nu(...) in opblkschbeta.c
*			  integrand(...) in opblkschbeta.c
*			  safe_exp(...) in opblkschbeta.c		
*		                                    
* PARAMETERS                                                           
*	INPUT	    : fwd_price	    - forward underlying price    	 
*              	: strike      	- strike price		      	  
*              	: vol         	- annual volatility   	      	  
*              	: mat         	- initial time, in years  	      	   	      
*				: beta          - dF = a*F^{beta}dW
*              	: call_put      - type of option: 0 call, 1 put 
*				: greek			- info wanted (premium, greeks...)  
*
*                                                      
* RETURNS      	: premium       - option premium		      	 
*                                                                       
*******************************************************************************/

/*******************************************************************************
**                      Include Files
*******************************************************************************/


#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"


/* ========================================================================
   FUNC: integrand  (static)
   DESC: Function under numerical integration for the B-S beta.
   MODIFIES:
   DECLARATION:
   ======================================================================== */

static double safe_exp(double x)
{
	if (x < -27)
	  return 0.;
	else
	  return exp(x);
}

static double integranda(double y,va_list argptr)
{/* exponent < .5 */
	
	double nu;
	double x_0;
	double k_bar;
	double ans;
	double I_minusnu,I_posnu;
	Err	   err = NULL;

	nu = va_arg(argptr, double);
	x_0 = va_arg(argptr, double);
	k_bar = va_arg(argptr, double);

	err = I_nu(-nu,(2*x_0*x_0*y),&I_minusnu);
	err = I_nu( nu,(2*x_0*x_0*y),&I_posnu);
	
	ans = (pow(y,2*nu) - pow(k_bar,2*nu))*pow(y,.5-nu)*
		safe_exp(-x_0*x_0*((y-1)*(y-1)))*(2*x_0*x_0) * sqrt(y)*
		(.5*I_posnu+.5*I_minusnu);


	return ans;
}

static double integrandb(double y,va_list argptr)
{/* exponent >= .5 */
	
	double safe_exp(double x);

	double nu;
	double x_0;
	double k_bar;
	double ans;
	double I_posnu;
	Err	   err = NULL;
	
	nu = va_arg(argptr, double);
	x_0 = va_arg(argptr, double);
	k_bar = va_arg(argptr, double);

	err = I_nu( nu,(2*x_0*x_0*y),&I_posnu);

	ans = (pow(y,2*nu) - pow(k_bar,2*nu))*pow(y,.5-nu)*
		safe_exp(-x_0*x_0*((y-1)*(y-1)))*(2*x_0*x_0) * 
	    sqrt(y)*I_posnu;

	return ans;
}


/* ----------------------------------------------------------------------------------- */
/*   NEW CODE WITH EXACT FORMULA                                                       */
/*   NOW ABLE THE BETA TO BE GREATER THAN ZERO                                         */
/* ----------------------------------------------------------------------------------- */

double srt_f_optblkschbeta(
		double            fwd_price,							
		double            strike, 
		double            vol, 
		double            mat,  
		double            beta,
		double            disc,
		SrtCallPutType    call_put, 
		SrtGreekType      greek)
{
/* parameters */
double b;	/* modified bessel I_{nu} */
double a;	/* change of var:   
			 {K^2{1-beta}}{sigma^2*(T-t)*(1-beta)^2}
			 */
double c;	/* change of var:
			{F^2{1-beta}}{sigma^2*(T-t)*(1-beta)^2}
			 */ 
double varc;   
double ans = 0.; 
double shift;
double premium;
double term1;
double term2;
double  precision1;
double  precision2;
Err err=NULL;
					/*** VARIABLES CHECK ***/

	if ((mat <=0) || (vol <=0) || (pow(fwd_price, beta-1) * vol * sqrt(mat) <= 1e-4))
	{
		if (call_put == SRT_PUT)
				premium = FMAX(strike - fwd_price,0);
		else
				premium = FMAX(fwd_price - strike,0);
	}
	else
	{
		if (fabs(beta - 1.0) <= 1E-2)
		{
			premium = 	srt_f_optblksch(fwd_price, strike, vol, mat, disc, call_put, greek);
		}
		else
		if (beta < 0)
		{
			premium = 	srt_f_optblknrm(fwd_price, strike, vol, mat, disc, call_put, greek);
		}
		else
		{
							/*** Change of variables ***/
			/* level 1 */
			b  = 1./(1.-beta);
			
			/* cumulative variance */
			varc = vol*vol*mat;

			/* level 2 */
			c = pow(fwd_price,2. - 2.*beta)/(varc*(1. - beta)*(1.-beta));
			a = pow(strike,2. - 2.*beta)/(varc*(1. - beta)*(1.-beta));

			if (beta <  1)
			{
			 err = cumchn(a,2.+b,c,&term1,&precision1);
			 if (err)
				 return -1;

			 err = cumchn(c,b,a,&term2,&precision2);
			 if (err)
				 return -1;

				premium = disc*(fwd_price*(1.-term1)-strike*term2);
			}
			else
			{
			 err = cumchn(c,-b,a,&term1,&precision1);
			 if (err)
				 return -1;
			 err = cumchn(a,2.-b,c,&term2,&precision2);
			 if (err)
				 return -1;

				premium = disc*(fwd_price*(1.-term1)-strike*term2);
			}
			/* Use put-call parity for computing put */
			if ( call_put == SRT_PUT )
					premium = premium + disc * (strike - fwd_price);
		}

	}


switch (greek)
{
	 case PREMIUM:    
		ans = premium;
		break;     

	case DELTA :
		shift = fwd_price / 10000;
		ans = ( srt_f_optblkschbeta(	
					fwd_price  + shift, 
		 			strike, 	
		 			vol,
		 			mat, 
					beta,
					disc,
					call_put,
					PREMIUM) 
			- premium ) / shift;  
		break;         

	case GAMMA :
		shift = fwd_price/ 10000;
		ans =  srt_f_optblkschbeta(	
					fwd_price +  shift, 
		 			strike, 	
		 			vol,
		 			mat,
					beta,
					disc,
					call_put,
					PREMIUM);  
		ans +=  srt_f_optblkschbeta(	
					fwd_price - shift, 
		 			strike, 	
		 			vol,
		 			mat,
					beta,
					disc,
					call_put,
					PREMIUM);  
		ans -= 2*premium;
		ans /= shift*shift; 
                break;

	case VEGA:
		shift = GVOPT.vol_add;
		ans = ( srt_f_optblkschbeta(	
					fwd_price , 
		 			strike, 	
		 			vol+ shift,
		 			mat, 
					beta,
					disc,
					call_put,
					PREMIUM) 
			- premium ) / shift;  
		break;


	case THETA:
		shift = YEARS_IN_DAY;
		ans =  srt_f_optblkschbeta(	
					fwd_price, 
		 			strike, 	
		 			vol,
		 			mat - shift, 
					beta,
					disc*exp(-shift*log(disc)/mat),
					call_put,
					PREMIUM) 
			- premium;  
		break;

	default:
		ans = UNKNOWN_GREEK;
}

	return ans;
}

