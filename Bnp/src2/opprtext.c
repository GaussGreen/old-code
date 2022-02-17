/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/

/******************************************************************************
 	computes the following probabilities:
		Prob ( x(T2) > k ; m(T1) > b ): 	cp = 1; du = 1 
		Prob ( x(T2) < k ; m(T1) > b ): 	cp = -1; du = 1 
		Prob ( x(T2) > k ; M(T1) < b ): 	cp = 1; du = -1 
		Prob ( x(T2) < k ; M(T1) < b ): 	cp = -1; du = -1 
******************************************************************************/
static double proba_joint_spot_part_min_max(
				double mu1,
				double mu2,
				double sig1,
				double sig2,
				double mat1,
				double mat2,
				double k,
				double b,
				double cp,   /*1 if call:S>K  -1 if put:S<K */
				double du    /*1 if down:min>B  -1 if up:Max<B*/
                                )
{                                      
double part1,part2,scalar,prob;
double sig_sqrt1 = sig1 * sqrt(mat1);
double sig_sqrt2 = sig2 * sqrt(mat2);
double corr		 = sig_sqrt1/sig_sqrt2;	

	part1 	= bivar( du * ( -b + mu1*mat1 ) / sig_sqrt1,     
			 cp * ( -k + mu2*mat2 ) / sig_sqrt2,
			 du * cp * corr );       

	part2 	= bivar( du * (b + mu1*mat1) / sig_sqrt1,
			 cp * (-k + 2*b + mu2*mat2) / sig_sqrt2,  
			 du * cp * corr);

	scalar 	= exp(2*mu1*b/(sig1*sig1));       

	prob 	= part1 - scalar * part2;

	return(prob);
}


/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optprtext( 11 )                                
*                                                                         
* PURPOSE      	: Extinguishable option                                    
*                                                                         
* DESCRIPTION  	: Extinguishable option where the extinguishable period	  
*		  	(mat1) is shorter than the option maturity (mat2) 
*
* CALLS		: bivar() in gen_math.c
*		: srt_f_optblksch()
*                                                                         
* PARAMETERS   	: fwd1     	- forward price of 1st underlying    	  
*              	: fwd2       	- forward price for 2nd underlying     	  
*              	: spot       	- spot price of underlying            	  
*              	: strike     	- strike price of option          	  
*              	: barrier     	- barrier level        	      		  
*              	: vol         	- volatility of underlying           	  
*              	: mat1        	- maturity of 1st extinguishable period   
*              	: mat2        	- maturity of 2nd extinguishable period   
*              	: disc         	- discount factor            	      	  
*              	: call_put      - ??            	      		  
*              	: down_up    	- ??                    		  
* 		: greek	- ??
*                                                                         
* RETURNS      	: ??          	- ??                                      
*                                                                         
*******************************************************************************/

double 	srt_f_optprtext  	(
				double 	fwd1, 
				double 	fwd2, 
				double 	spot, 
				double 	strike,
				double 	barrier, 
				double 	vol1, 
				double 	vol2, 
				double 	mat1, 
				double 	mat2, 
				double 	disc, 
				SrtCallPutType 		call_put, 
				SrtBarrierType 		down_up, 
				SrtGreekType 		greek)

{ /* BEGIN srt_f_optprtext(11) */

double premium;
double answer;
double shift;
double mu_1, mu_1p, mu_1m, mu_2, mu_2p, mu_2m;
double k,b;
double du, cp;


	cp = ( call_put == SRT_CALL ) ? 1 : -1 ;
	du = ( down_up == SRT_DOWN ) ? 1 : -1 ;

	/* Spot should better be above the barrier for a Down&Out*/
	/* and better be below the barrier for an Up&Out */
	if ( (spot - barrier) * du < 0 )
		premium = 0.0;
	else
	/*** PARTIAL EXTINGUISHABLE => EXTINGUISHABLE ***/
	if ( mat1 >= mat2 )   
	{
		premium = srt_f_optexting(	fwd2,
									spot,
									strike,
									barrier,
									vol2,
									mat2,
									disc,
									call_put,
									down_up,
									PREMIUM);
	}                 
	else	
	if (mat2 == 0) /* so we don't care about mat1 */
	{
		if ( (spot - strike) * cp > 0 ) 
			premium = cp * (spot - strike) ;
		else
			premium =  0.0;
	}
	else 
	if (mat1 <= 0)   /* Standard European  */
	{
		premium = srt_f_optblksch(	fwd2,
									strike,
									vol2,
									mat2,
									disc,
									call_put,
									PREMIUM);
	}
	else
	if ( vol1 == 0 )
	{
		if ( (fwd1 - barrier) * du > 0 )   
			premium = srt_f_optblksch(	fwd2,
										strike,
										vol2,
										mat2,
										disc,
										call_put,
										PREMIUM);
		else 
			premium = 0.0;
	}
	else
	{
		k = log(strike/spot);
		b = log(barrier/spot);

		mu_1 = log(fwd1/spot) / mat1;                                
		mu_2 = log(fwd2/spot) / mat2;                                

		mu_1p = mu_1 + 0.5 * vol1 * vol1;
		mu_1m = mu_1 - 0.5 * vol1 * vol1; 
		mu_2p = mu_2 + 0.5 * vol2 * vol2;
		mu_2m = mu_2 - 0.5 * vol2 * vol2; 

		premium = 	fwd2 * 
					proba_joint_spot_part_min_max( 	mu_1p,
													mu_2p,
													vol1,
													vol2,
													mat1,
													mat2,
													k,
													b,
													cp,   
													du);

		premium -=  strike * 
					proba_joint_spot_part_min_max(	mu_1m,
							        				mu_2m,
													vol1,
													vol2,
													mat1,
													mat2,
													k,
													b,
													cp,   
													du);
		premium *= cp * disc;
	}

	switch (greek)
	{
		case PREMIUM : 		/*** PREMIUM ***/  
			answer = premium;
			break;

		case DELTA_FWD1 : 	/*** DELTA FWD1 ***/  
			shift = fwd1 / 10000;
			answer = (srt_f_optprtext(	fwd1+ shift,
										fwd2,
			   							spot,
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
						- premium ) / shift;
			break;

		case DELTA_FWD2 : 	/*** DELTA FWD2 ***/  
			shift = fwd2 / 10000;
			answer = (srt_f_optprtext(	fwd1,
										fwd2+ shift,
			   							spot,
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
					- premium ) / shift;
			break;

		case DELTA : 		/*** DELTA FWD1 + FWD2 + SPOT ***/  
			shift = spot / 10000;
            answer = (srt_f_optprtext(	fwd1 *(1 + shift/spot),
										fwd2 *(1 + shift/spot),
			   							spot + shift,
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
					- premium ) / shift;
			break;

		case GAMMA : 		/*** GAMMA ***/  
			shift = spot / 10000;
            answer = srt_f_optprtext(	fwd1 * (1 + shift/spot),
										fwd2 * (1 + shift/spot),
			   							spot + shift,
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
			answer += srt_f_optprtext(	fwd1 * (1 - shift/spot),
										fwd2 * (1 - shift/spot),
			   							spot - shift,
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
			break;

		case VEGA1: 		/*** VEGA ***/
			shift = GVOPT.vol_add;
            answer = (srt_f_optprtext(	fwd1 ,
										fwd2 ,
			   							spot,
										strike,
										barrier,
										vol1 + shift,
										vol2 ,
										mat1,
										mat2,
										disc,
										call_put,
										down_up,
										PREMIUM)
					- premium ) / shift;
			break;

		case VEGA2: 		/*** VEGA ***/
			shift = GVOPT.vol_add;
			answer = (srt_f_optprtext(	fwd1 ,
										fwd2 ,
			   							spot ,
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
					- premium ) / shift;
			break;

		case THETA: 		/*** THETA ***/
			shift = YEARS_IN_DAY;
			answer = (srt_f_optprtext(	fwd1 ,
										fwd2 ,
			   							spot ,
										strike,
										barrier,
										vol1 ,
										vol2 ,
										mat1 -  shift,
										mat2 - shift,
										disc
											*exp(-shift*log(disc)/mat2),
										call_put,
										down_up,
										PREMIUM)
						- premium ) / shift;
			break;

		default:
			answer = UNKNOWN_GREEK;
	}                 

	return(answer);

} /* END   srt_f_optprtext() */ 
