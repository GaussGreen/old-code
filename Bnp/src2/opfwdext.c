/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/

static double 	optfwdext_fn(
			double mu1,
			double mu2,
			double sig1,
			double sig2,
			double mat1,
			double mat2,
			double k,
			double l,
			int    du );      

static double 	srt_f_optfwdext_price  	(
				double 	fwd1, 
				double 	fwd2, 
				double 	strike,
				double 	barrier, 
				double 	vol1, 
				double 	vol2, 
				double 	mat1, 
				double 	mat2, 
				double 	disc, 
				SrtCallPutType	call_put, 
				SrtBarrierType	down_up);

/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optfwdext()                               
*                                                                         
* PURPOSE      	: Forward extinguishable option  
*                                                                         
* DESCRIPTION  	: Forward extinguishable option where the extinguishable  
*		  period is from mat1 to mat2. 
*
* CALLS		: bivar() in gen_math.c
*                                                                         
* PARAMETERS   	: fwd1     	- forward price of underlying at T1   	  
*              	: fwd2       	- forward price of underlying at T2    	  
*              	: strike     	- strike price of option          	  
*              	: barrier     	- barrier level        	      		  
*              	: vol         	- volatility of underlying           	  
*              	: mat1        	- maturity of 1st extinguishable period   
*              	: mat2        	- maturity of 2nd extinguishable period   
*              	: disc         	- discount factor            	      	  
*              	: call_put        	- 0/1 - call/put          	  
*              	: down_up    	- 0/1 - barrier down or up	  
*                                                                         
* RETURNS      	: premium       - of option        
*               : deltas,
*		: gamma,
*		: vegas,
*		: theta                                                           
*******************************************************************************/

double 	srt_f_optfwdext  	(
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
				SrtGreekType 	greek) 

{ /* BEGIN srt_f_optfwdext() */
double	answer;
double	premium;
double	shift;

premium = srt_f_optfwdext_price (
				 fwd1, 
				 fwd2, 
				 strike,
				 barrier, 
				 vol1, 
				 vol2, 
				 mat1, 
				 mat2, 
				 disc, 
				 call_put, 
				 down_up ) ;
        
switch (greek)
{
	case PREMIUM : 	/*** PREMIUM ***/  
		answer = premium;
		break;

	case DELTA_FWD1 : 	/*** DELTA FWD ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_optfwdext(	fwd1+ shift, 
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
	case DELTA_FWD2 : 	/*** DELTA FWD ***/  
		shift = fwd2 / 10000;
		answer = ( srt_f_optfwdext(	fwd1, 
						fwd2+ shift,
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

	case DELTA :	/*** DELTA FWD1 + FWD2 ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_optfwdext(	fwd1 * (1 + shift/fwd1),
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
		answer =  srt_f_optfwdext(	fwd1 * (1 + shift/fwd1),
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
		answer += srt_f_optfwdext(	fwd1 * (1 - shift/fwd1),
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
	case VEGA1 : 	/*** VEGA AT T1***/
		shift = GVOPT.vol_add;
		answer =  ( srt_f_optfwdext(	fwd1,
						fwd2 ,
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
			- premium )  / shift;
		break;
	case VEGA2 : 	/*** VEGA AT T2***/
		shift = GVOPT.vol_add;
		answer =  ( srt_f_optfwdext(	fwd1,
						fwd2 ,
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
		answer =  srt_f_optfwdext(	fwd1,
						fwd2 ,
						strike,
				        barrier,
						vol1,
						vol2 ,
						mat1 - shift,
						mat2 - shift,
						disc
							*exp(shift*log(disc)/mat2),
						call_put,
						down_up,
						PREMIUM)
			- premium ;
		break;

	default:
		answer = UNKNOWN_GREEK;
		break;
}
return(answer) ;
} /* END   srt_f_optfwdext() */ 
/************************************************************************/


static double 	srt_f_optfwdext_price  	(
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
				SrtBarrierType 	down_up)

{ /* BEGIN srt_f_optfwdext_price */
double	premium = 0.0;
double 	measure = 0.0;
double 	mu, mu1, mu2,k, l, sig1, sig2;

int	du, cp, flag;


	cp 	 = ( call_put == SRT_CALL ) ? 1 : -1 ;   
	du 	 = ( down_up == SRT_DOWN ) ? 1 : -1 ;
	flag = du*cp;

	/* the spot will cross the barrier if we want a payoff different from 0 */
	/* so premium = 0 */
    if ( flag == -1)
	{ 
		if ( du*(strike-barrier ) < 0.0  ) 
			return(0.0);
    }

	/* Maturity 2 should be greater than maturity 1 */
	/* otherwise we have an european */ 
	if ( mat1 >= mat2 )
	{
		if (flag == -1) /* Call Up or Put Down */
		{
			/* It is an european with strike K minus an european with strike L */
			/* minus a digital multiply by du*(strike-barrier) */
			premium = srt_f_optblksch(	fwd1,
										strike,
										vol1,
										mat1,
										disc,
										call_put,
										PREMIUM);

			premium -= srt_f_optblksch(	fwd1,
										barrier,
										vol1,
										mat1,
										disc,
										call_put,
										PREMIUM);

			premium -= du * (strike-barrier) *
					   srt_f_opteurdig(	fwd1,
										barrier,
										vol1,
										mat1,
										disc,
										down_up,
										PREMIUM);


		}
		else /* Call Down or Put Up */
		{
			/* It is an european */
			premium =  srt_f_optblksch(	fwd1, 
										strike,
										vol1,
										mat1,
										disc,
										call_put,
										PREMIUM);
		}

		return premium;
	}


	k 		= log( strike / fwd1 );
	l		= log( barrier / fwd1 );
    sig1 	= vol1;
	sig2 	= sqrt( (vol2*vol2*mat2 - vol1*vol1*mat1)/(mat2-mat1) ) ;
                   
	mu		= log(fwd2/fwd1)/(mat2-mat1);
	

	/* First part of the formula */
	mu1 	= sig1 * sig1 / 2;     
	mu2 	= mu + sig2 * sig2 / 2;

	measure = 0.0;

	if ( flag == -1)
	{         	
		measure = optfwdext_fn(	mu1 , mu2,
								sig1, sig2,
								mat1,mat2,
								l,l, du ); 
	}                

	measure += flag * optfwdext_fn(	mu1, mu2,
									sig1, sig2,
									mat1, mat2,
									k,l, du );      
                  
	premium = fwd2 * measure;



	/* Second part of the formula */
	measure  = 0.0 ;

	mu1 	= ( - sig1*sig1 / 2 ) ;     
	mu2 	= mu - sig2*sig2 / 2 ;
    
	if ( flag == -1)
	{ 	
		measure = optfwdext_fn(	mu1 , mu2,
								sig1, sig2,
								mat1,mat2,
								l,l, du );        
	}

	measure += flag * optfwdext_fn(	mu1 , mu2,
									sig1, sig2,
									mat1, mat2,
									k,l, du );

  	premium -= strike * measure;

	/* Price */
	premium *= disc * cp;

	return (premium);

} /* END   srt_f_fwdext_price */ 


/******************************************************************************/
static double 	optfwdext_fn (
			double	mu1,
			double	mu2,
			double	sig1,
			double	sig2,
			double	mat1,
			double	mat2,
			double	k,
			double	l,
			int		du )      

/* calculates the prob pr( X(t2)<K, min(t1,t2) X > l )	*/

{	
double	part1, part2,result;
double	du1, a, b, a1, b1, alpha, mat12;

	du1 	= (-1.0 * du);		
	mat12	= mat2 - mat1;

	a	= du1 * (k - l - mu2 * mat12)  ;
    b	= sig2 * sqrt( mat12 ) ;
    a1	= du1 * (l - mu1 * mat1) ;          
    b1	= sig1 * sqrt( mat1 ) ;

	part1	= srt_f_intnrm_i1( a, b, a1, b1 );

        
   	alpha 	= du1 * 2.0 * mu2 / (sig2*sig2) ;

	part2	= srt_f_intnrm_j2(alpha, a, b, a1, b1);

	result 	= part1 -  part2; 

	return result;

} /* END optfwdext_fn  */
