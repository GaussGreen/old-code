/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                           
* FUNCTION     	: srt_f_optimpvolbeta()                                             
*                                                                           
* PURPOSE      	: Calculates the volatility implied from a premium           
*		  	computed using the Black-Scholes Beta formula     
*                                                                           
* DESCRIPTION  	: XX                                                         
*		  							    
* CALLS		: srt_f_optblkschbeta()                                      
*                                                                           
* PARAMETERS                                                                
*	INPUT	: premium 	- premium on European option                
*              	: fwd_price   	- forward price of underlying   	    
*              	: strike        - strike price               	      	    
*              	: vol_init      - initial volatility of underlying 	    
*              	: mat           - maturity, in years 
*				: beta			- dF = a*F^{beta}dW	      		    
*              	: disc        	- discount factor to expiry 	      	    
*              	: call_put      - type of option: 0 call, 1 put	      	    
*                                                                           
* RETURNS      	: vol_guess     - guess of implied volatility		    
*                                                                           
*******************************************************************************/

static double optvalbeta(double fwd,double strike,double vol,double mat, double beta, double disc_fact,
					SrtCallPutType call_put, SrtGreekType greek)
{
	
		return srt_f_optblkschbeta(fwd,strike,vol,mat,beta,disc_fact,call_put,greek);  

}

static double optval(double fwd,double strike,double vol,double mat, double disc_fact,
					SrtCallPutType call_put, SrtGreekType greek)
{
	
		return srt_f_optblksch(fwd,strike,vol,mat,disc_fact,call_put,greek);  

}

double 	srt_f_optimpvolbeta	(
			double 	premium, 
			double 	fwd_price, 
			double 	strike,
			double 	mat, 
			double  beta,
			double 	disc, 
			SrtCallPutType 		call_put
			)
{ 	
	int 	i;
	double  vol_init = VOL_INIT;
	double  vol_shift = VOL_SHIFT;
	double 	vol_guess;
	double	tmp_vol;
	double 	deriv;
	double 	prem_shift;
	double	prem_new;


/*  nab  */
	vol_init = beta + .01;
	vol_shift = .001*(beta + .1);

/*
if normal vols then1
	vol_init *= .1;
	vol_shift *= .1;
*/

vol_guess = vol_init;

#if 1

/* ========================================================================== */
/* XXX       	      	      		      */
/* ========================================================================== */
   
	

prem_new  = optval(	fwd_price,
				strike,
				NULL_VOL,
				mat, 
				disc, 
				call_put, 
				PREMIUM
				);  

/* ========================================================================== */
/* XXX       	      	      		      */
/* ========================================================================== */
                          
if (  ( fabs( prem_new - premium ) < premium * PREM_TOL_PCT )
   || ( premium < PREM_TOL_PCT/10.0) ) 
{
	vol_guess = NULL_VOL;
}
else

#endif

{
	/* ------------------------------------------------------------------ */
	/* XXX  	      */
	/* ------------------------------------------------------------------ */

	for (i = 0; i <= MAX_ITER; i++)
    {

		/* ---------------------------------------------------------- */
		/* XXX  	      */
		/* ---------------------------------------------------------- */

		prem_new = optvalbeta(fwd_price,
						strike,
						vol_guess,
						mat, 
						beta,
						disc, 
						call_put, 
						PREMIUM
						);  

		/* ---------------------------------------------------------- */
		/* XXX  	      */
		/* ---------------------------------------------------------- */

	    if ( fabs( prem_new - premium ) < PREM_TOL ) 
		{
			break;
        }

		/* ---------------------------------------------------------- */
		/* XXX  	      */
		/* ---------------------------------------------------------- */

    	prem_shift = optvalbeta(fwd_price,
						strike,
						vol_guess + vol_shift,
						mat, 
						beta,
						disc, 
						call_put, 
						PREMIUM
						);

		deriv = (prem_shift - prem_new) / vol_shift;
	     	
		/* ---------------------------------------------------------- */
		/* XXX  	      */
		/* ---------------------------------------------------------- */

		if ( deriv != 0 ) 
		{	
			/* -------------------------------------------------- */
			/* XXX  	      */
			/* -------------------------------------------------- */

	        tmp_vol = vol_guess - ( (prem_new-premium) / deriv );

		  	if ( tmp_vol > (vol_guess * 10) )
				vol_guess = vol_guess * 2.0;
		  	else if ( tmp_vol < (vol_guess / 10) ) 
				vol_guess = vol_guess / 2.0;	
		  	else
				vol_guess = tmp_vol;
		}	
	    else 
		{
			/* -------------------------------------------------- */
			/* XXX  	      */
			/* -------------------------------------------------- */

		  	if ( prem_new < premium )
				vol_guess = vol_guess * 2.0;
		  	else
				vol_guess = vol_guess / 1.5;
		}		
	}
}

/* ========================================================================== */
/* Return implied volatility guess       	      	      		      */
/* ========================================================================== */

return( vol_guess );

} /* srt_f_optimpvolbeta() */

/******************************************************************************/
