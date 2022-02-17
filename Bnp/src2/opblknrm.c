/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optblknrm(...)                                   
*                                                                       
* PURPOSE      	: European option pricing for a single risk factor     
*                                                                      
* DESCRIPTION  	: Black-Scholes: calculates the option premium and     
*		  theoretical delta, gamma, vega and theta for an 	
*		  underlying which is lognormal                        
*                                                                   
* CALLS		: gauss(...) in gen_math.c
*		  norm(...)  in gen_math.c	      
*		                                    
* PARAMETERS                                                           
*	INPUT	: fwd_price	- forward underlying price    	 
*              	: strike      	- strike price		      	  
*              	: vol         	- annual volatility   	      	  
*              	: mat         	- maturity, in years  	      	       
*              	: disc        	- discount factor to expiry   		 
*              	: call_put      - type of option: 0 call, 1 put 
*		: greek		- info wanted (premium, greeks...)
*                                                                        
* RETURNS      	: premium       - option premium		      	 
*                                                                       
*******************************************************************************/

double 	srt_f_optblknrm(	
		double 				fwd_price,
		double 				strike,
		double 				vol,
		double 				mat, 
		double 	  			disc, 
		SrtCallPutType		call_put, 
		SrtGreekType		greek
		)
{			
	double 	d1;
	double  answer;
	double 	nd1;
	double 	gd1;
	double 	cp;
	double  vsqrtt;
 
	if (mat < 0) return 0;

	if (call_put == SRT_STRADDLE)
	{
		answer = srt_f_optblknrm(fwd_price, strike,vol,mat,disc,SRT_CALL,greek)+srt_f_optblknrm(fwd_price, strike,vol,mat,disc,SRT_PUT,greek);
		return answer;
	}
	
	cp = (call_put == SRT_CALL) ? 1 : -1 ;
        
	vsqrtt = vol * sqrt(mat);

/* Set-up d1 and d1 (see Hull, p. 224) to avoid divide-by-zeros 	      */

if ((vol != 0)  && (mat != 0)) 
{
	/* Expected  case	      				      	      */

	d1  = (fwd_price-strike)/vsqrtt; 
	nd1 = norm(cp*d1);
    gd1 = gauss(d1);
}
else 
{
	/* Allow intrinsic value of underlying to be returned if required     */

	if (cp * (fwd_price - strike) > 0.0)
	{                     
		d1  = INFINITY; 			/* Equivalent of +oo */
 	    nd1 = 1;
	}
	else 
	{
		d1  = -INFINITY;			/* Equivalent of -oo */
 	    nd1 = 0;
	}
    gd1 = 0.0;
}		 			
		
/* Evaluate option-type dependent greeks (spot-deltas, theta) and premium     */


if ((mat == 0) && (greek!=PREMIUM))
	answer = 0.0;
else
{
	switch (greek)
	{
	case DELTA: 		/*** DELTA ***/
	case DELTA_FWD: 	/*** FWD DELTA ***/
		answer = cp * nd1 * disc;
		break;

	case GAMMA : 		/*** GAMMA ***/
	case GAMMA_FWD : 	/*** FWD GAMMA ***/
		answer = gd1/vsqrtt;
		break ;

	case VEGA : 		/*** VEGA ***/
		answer =  0.01 * disc * sqrt(mat) * gd1;
		break ;

	case THETA : 		/*** THETA ***/
		answer = - (1 / 365.00) *
			(0.5 * disc * vol/sqrt(mat) * gd1 + 
			  disc * log(disc)/mat);
      	break ;

	case PREMIUM : 		/*** PREMIUM ***/   
		answer = disc * (vsqrtt*gd1 + cp * (fwd_price-strike) * nd1);
		break;

	case RHO :			/*** RHO ***/
		answer =  disc * cp * strike * mat * nd1;   
		break;

	case VANNA :
		answer = - disc * d1 * gauss(d1) / vol;
		break;

	case VOLGA :
		answer = disc * sqrt(mat) * d1 * d1 * gauss(d1) / vol;		
		break;

	default:
		answer = UNKNOWN_GREEK;
		break;
	}
}
return (answer);

} /* END srt_f_optblknrm() */

/* ---------------------------------------------------------------------------------------- */

/* A function that checks the Normal Black-Scholes inputs */

Err srt_f_optchkblknrm(
		double 			vol,
		double 			mat)
{
	if (mat <0)
		return serror("Negative Maturity in Normal BlackScholes inputs");

	if (vol <0)
		return serror("Negative Volatility in Normal BlackScholes inputs");

	return NULL;
}


/* ---------------------------------------------------------------------------------------- */

/* A Higher level interface function to call Normal Black Scholes */
Err OptBlkNrm(	
		double 				fwd_price,
		double 				strike,
		double 				vol,
		double 				mat, 
		double 	  			disc, 
		char               *call_put_str, 
		char               *greek_str,
		double             *result
		)
{
Err                 err             = NULL;
SrtCallPutType      call_put_type;
SrtGreekType        greek_type      = PREMIUM;


/* Transforms the call_put string into a type */
	err = interp_call_put(call_put_str, &call_put_type);
	if (err)
		return err;

/* Transforms the greek string into a type */
	if (strcmp(greek_str, "") != 0)
	{
		err = interp_greeks(greek_str, &greek_type);
		if (err)
			return err;
	}

/* Check the inputs */
	err = srt_f_optchkblknrm(vol, mat);
	if (err)
		return err;

/* Compute the requested greek */
	*result = srt_f_optblknrm(fwd_price, strike, vol, mat,
		disc, call_put_type, greek_type);

/* Return a success message */
	return NULL;

}

double 	srt_f_optblknrm_accurate(	
		double 				fwd_price,
		double 				strike,
		double 				vol,
		double 				mat, 
		double 	  			disc, 
		SrtCallPutType		call_put, 
		SrtGreekType		greek
		)
{			
	double 	d1;
	double  answer;
	double 	nd1;
	double 	gd1;
	double 	cp;
	double  vsqrtt;
 
	if (mat < 0) return 0;

	if (call_put == SRT_STRADDLE)
	{
		answer = srt_f_optblknrm_accurate(fwd_price, strike,vol,mat,disc,SRT_CALL,greek)+srt_f_optblknrm(fwd_price, strike,vol,mat,disc,SRT_PUT,greek);
		return answer;
	}
	
	cp = (call_put == SRT_CALL) ? 1 : -1 ;
        
	vsqrtt = vol * sqrt(mat);

/* Set-up d1 and d1 (see Hull, p. 224) to avoid divide-by-zeros 	      */

if ((vol != 0)  && (mat != 0)) 
{
	/* Expected  case	      				      	      */

	d1  = (fwd_price-strike)/vsqrtt; 
	nd1 = norm_accurate(cp*d1);
    gd1 = gauss(d1);
}
else 
{
	/* Allow intrinsic value of underlying to be returned if required     */

	if (cp * (fwd_price - strike) > 0.0)
	{                     
		d1  = INFINITY; 			/* Equivalent of +oo */
 	    nd1 = 1;
	}
	else 
	{
		d1  = -INFINITY;			/* Equivalent of -oo */
 	    nd1 = 0;
	}
    gd1 = 0.0;
}		 			
		
/* Evaluate option-type dependent greeks (spot-deltas, theta) and premium     */


if ((mat == 0) && (greek!=PREMIUM))
	answer = 0.0;
else
{
	switch (greek)
	{
	case DELTA: 		/*** DELTA ***/
	case DELTA_FWD: 	/*** FWD DELTA ***/
		answer = cp * nd1 * disc;
		break;

	case GAMMA : 		/*** GAMMA ***/
	case GAMMA_FWD : 	/*** FWD GAMMA ***/
		answer = gd1/vsqrtt;
		break ;

	case VEGA : 		/*** VEGA ***/
		answer =  0.01 * disc * sqrt(mat) * gd1;
		break ;

	case THETA : 		/*** THETA ***/
		answer = - (1 / 365.00) *
			(0.5 * disc * vol/sqrt(mat) * gd1 + 
			  disc * log(disc)/mat);
      	break ;

	case PREMIUM : 		/*** PREMIUM ***/   
		answer = disc * (vsqrtt*gd1 + cp * (fwd_price-strike) * nd1);
		break;

	case RHO :			/*** RHO ***/
		answer =  disc * cp * strike * mat * nd1;   
		break;

	case VANNA :
		answer = - disc * d1 * gauss(d1) / vol;
		break;

	case VOLGA :
		answer = disc * sqrt(mat) * d1 * d1 * gauss(d1) / vol;		
		break;

	default:
		answer = UNKNOWN_GREEK;
		break;
	}
}
return (answer);

} /* END srt_f_optblknrm_accurate() */
