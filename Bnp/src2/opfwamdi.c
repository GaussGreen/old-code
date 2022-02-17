/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optfwdamdig( 7)                                   
*                                                                       
* PURPOSE      	: Forward American Digital, with payment at maturity     
*                                                                      
* DESCRIPTION  	: ??   
*		                                    
* PARAMETERS      ???     
*	INPUT	: fwd_price1	- forward underlying price    	 
		: fwd_price2
*              	: strike      	- strike price		      	  
*              	: vol1         	- annual volatility   	      	  
*              	: vol2         	- annual volatility   	      	  
*              	: mat1         	- maturity, in years  	      	       
*		: mat2
*              	: disc        	- discount factor to expiry   		 
*              	: below_above   - position of barrier       	 
*              	: min_max  	-       	 
*	OUTPUT	: greeks        - greeks for given option pricing     	  
*                                                                        
* RETURNS      	: premium       - option premium		      	 
*                                                                       
*******************************************************************************/

double 	srt_f_optfwdamdig(	
		double 			fwd1,
		double 			fwd2,
		double			barrier,
		double 			vol1,
		double			vol2,
		double 			mat1,
		double 			mat2, 
		double 			disc, 
		SrtBarrierType		below_above,
		SrtMinmaxType 		min_max,
		SrtGreekType		greek)
{			
	double	premium          = 0.0;
	double	answer          = 0.0;
	double	shift          = 0.0;


	premium = srt_f_optamesla (	fwd1,
								fwd2,
								barrier,
								vol1,
								vol2,
								mat1,
								mat2,
								disc,
								below_above,
								PREMIUM);

	premium +=  srt_f_opteurdig(fwd1,
								barrier,
								vol1,
								mat1,
								disc,
								below_above,
								PREMIUM);


switch (greek)
{
	case PREMIUM : 	/*** PREMIUM ***/  
		answer = premium;
		break;

	case DELTA_FWD1 : 	/*** DELTA FWD ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_optfwdamdig(	fwd1 + shift,
										fwd2,
				        				barrier,
										vol1,
										vol2,
										mat1,
										mat2,
										disc,
										below_above,
										min_max,
										PREMIUM)
				- premium ) /shift;      		 
		break;

	case DELTA_FWD2 : 	/*** DELTA FWD 2 ***/  
		shift = fwd2 / 10000;
		answer = ( srt_f_optfwdamdig(	fwd1,
										fwd2 + shift,
				        				barrier,
										vol1,
										vol2,
										mat1,
										mat2,
										disc,
										below_above,
										min_max,
										PREMIUM)
					- premium ) /shift;      		 
		break;

	case DELTA :	/*** DELTA SPOT + FWD ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_optfwdamdig(	fwd1 * (1 + shift/fwd1),
										fwd2 * (1 + shift/fwd1),
				        				barrier,
										vol1,
										vol2,
										mat1,
										mat2,
										disc,
										below_above,
										min_max,
										PREMIUM)
				- premium ) /shift;      		 
		break;

	case GAMMA :	/*** GAMMA ***/
		shift = fwd1 / 10000;
		answer =  srt_f_optfwdamdig(fwd1 * (1 + shift/fwd1), 
									fwd2 * (1 + shift/fwd1), 
				        			barrier,
									vol1,
									vol2,
									mat1,
									mat2,
									disc,
									below_above,
									min_max,
									PREMIUM);
		answer += srt_f_optfwdamdig(fwd1 * (1 - shift/fwd1),
									fwd2 * (1 - shift/fwd1), 
				        			barrier,
									vol1,
									vol2,
									mat1,
									mat2,
									disc,
									below_above,
									min_max,
									PREMIUM);
		answer -= 2 * premium;
		answer /= shift * shift;
		break ;

	case VEGA1 : 	/*** VEGA ***/
		shift = GVOPT.vol_add;
		answer =  ( srt_f_optfwdamdig(	fwd1 ,
										fwd2,
				        				barrier,
										vol1 + shift,
										vol2,
										mat1,
										mat2,
										disc,
										below_above,
										min_max,
										PREMIUM)
				- premium )  / shift;
		break;

	case VEGA2 : 	/*** VEGA ***/
		shift = GVOPT.vol_add;
		answer =  ( srt_f_optfwdamdig(	fwd1 ,
										fwd2,
				        				barrier,
										vol1 ,
										vol2 + shift,
										mat1,
										mat2,
										disc,
										below_above,
										min_max,
										PREMIUM)
				- premium )  / shift;
		break;

	case THETA:		/*** THETA  ***/
		shift = YEARS_IN_DAY;
		answer =  srt_f_optfwdamdig(fwd1,
									fwd2,
				        			barrier,
									vol1 ,
									vol2, 
									mat1 - shift,
									mat2 - shift,
									disc
									 * exp(-shift*log(disc)/mat2),
									below_above,
									min_max,
									PREMIUM)
				- premium ;
		break;

	default:
		answer = UNKNOWN_GREEK;
		break;
}

return (answer);            

}  /*  END srt_f_opfwdamdig() */

