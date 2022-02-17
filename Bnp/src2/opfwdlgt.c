/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                        
* FUNCTION     	: srt_f_optfwdlgt( 11 )                                
*                                                                        
* PURPOSE      	: Forward Lightable Option                                    
*                                                                         
* DESCRIPTION  	: lightable option where the lightable period (mat1-mat2) is    
*		  starts after the option start              
*                                                                         
* CALLS		: srt_f_optblksch()                                          
*		: srt_f_optfwdext()                               	  
*                                                                         
* PARAMETERS   	: fwd1     	- forward price at T1	  
*              	: fwd2         	- forward price at T2	  
*              	: strike       	- strike price of ??           	   	  
*              	: barrier      	- barrier level           	   	  
*              	: vol1          - volatility of ??           	   	  
*              	: vol2         	- volatility of ??           	   	  
*              	: mat1        	- maturity of 1st extinguishable period   
*              	: mat2        	- maturity of 2nd extinguishable period   
*              	: disc        	- discount factor            		  
*              	: call_put      - ??            	 		  
*              	: down_up      	- ??                		      	  
*                                                                         
* RETURNS      	: ??          	- ??                                      
*                                                                         
*******************************************************************************/
	   	   
double 	srt_f_optfwdlgt	(	
				double 	fwd1, 
				double 	fwd2, 
				double 	strike,
				double 	barrier, 
				double 	vol1, 
				double 	vol2, 
				double 	mat1, 
				double 	mat2, 
				double 	disc, 
				SrtCallPutType 	call_put, 
				SrtBarrierType 	down_up,
				SrtGreekType 	greek
				)
{ 
double shift;
double premium, answer;

/* A forward extinguishable and a forward lightable is equal to an european */
premium = srt_f_optblksch(fwd2,
				strike,
				vol2,
				mat2,
				disc,
				call_put,
				PREMIUM)  
          - srt_f_optfwdext(fwd1,
				fwd2,
				strike,
				barrier,
				vol1,
				vol2,
				mat1,
				mat2,
		  		disc,
				call_put,
				down_up,
				PREMIUM);

switch (greek)
{
	case PREMIUM : 		/*** PREMIUM ***/  
		answer = premium;
		break;

	case DELTA_FWD1 : 	/*** DELTA FWD1 ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_optfwdlgt(	fwd1 + shift,
						fwd2,
						strike,
				        barrier,
						vol1,
						vol2,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM)
			- premium ) /shift;      		 
		break;

	case DELTA_FWD2 : 	/*** DELTA FWD2 ***/  
		shift = fwd2 / 10000;
		answer = ( srt_f_optfwdlgt(	fwd1 ,
						fwd2 + shift,
						strike,
				        barrier,
						vol1,
						vol2,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM)
			- premium ) /shift;      		 
		break;

	case DELTA :	/*** DELTA SPOT + FWDS ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_optfwdlgt(	fwd1 + shift,
						fwd2 * (1 + shift/fwd1),
						strike,
				        barrier,
						vol1,
						vol2,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM)
			- premium ) /shift;      		 
		break;

	case GAMMA :	/*** GAMMA ***/
		shift = fwd1 / 10000;
		answer =  srt_f_optfwdlgt(	fwd1 + shift,
						fwd2 * (1 + shift/fwd1),  
						strike,
				        barrier,
						vol1,
						vol2,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM);
		answer += srt_f_optfwdlgt(	fwd1 - shift,
						fwd2 * (1 - shift/fwd1), 
						strike,
				        barrier,
						vol1,
						vol2,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM);
		answer -= 2 * premium;
		answer /= shift * shift;
		break ;

	case VEGA1 : 	/*** VEGA1 ***/
		shift = GVOPT.vol_add;
		answer =  ( srt_f_optfwdlgt(	fwd1 ,
						fwd2,
						strike,
				        barrier,
						vol1 + shift,
						vol2,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM)
			- premium )  / shift;
		break;

	case VEGA2 : 	/*** VEGA2 ***/
		shift = GVOPT.vol_add;
		answer =  ( srt_f_optfwdlgt(	fwd1 ,
						fwd2,
						strike,
				        barrier,
						vol1 ,
						vol2 + shift,
						mat1,
						mat2,
						disc,
						call_put,
						down_up,
						PREMIUM)
			- premium )  / shift;
		break;

	case THETA:		/*** THETA  ***/
		shift = YEARS_IN_DAY;
		answer =  srt_f_optfwdlgt(	fwd1 ,
						fwd2 ,
						strike,
				        barrier,
						vol1,
						vol2,
						mat1 - shift,
						mat2 - shift,
						disc
							* exp(-shift*log(disc)/mat2),
						call_put,
						down_up,
						PREMIUM)
			- premium ;
		break;

	default:
		answer = UNKNOWN_GREEK;
		break;
}
return( answer );

} /* END srt_f_optfwdlgt() */

/******************************************************************************/
