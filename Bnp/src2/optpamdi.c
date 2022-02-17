/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


static double 	srt_f_opttpamdi_price  	(
				double 			fwd1, 
				double 			fwd2, 
				double 			barrier, 
				double 			vol1, 
				double 			vol2, 
				double 			mat1, 
				double 			mat2, 
				double 			disc1, 
				double 			disc2, 
				SrtBarrierType 	below_above);


/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_opttpamdi()                               
*                                                                         
* PURPOSE      	: "touch and pay" american digital option  
*                                                                         
* DESCRIPTION  	: American digital option where the payoff ($1) is paid
*		  when the barrier is touched  
*
* CALLS		: bivar() in gen_math.c
*                                                                         
* PARAMETERS   	: fwd1     	- forward price of underlying at T1    	  
*              	: fwd2       	- forward price of underlying at T2    	  
*              	: barrier     	- barrier level        	      		  
*              	: vol1         	- volatility of underlying(BS) from 0 to T1
*              	: vol2         	- volatility of underlying(BS) from 0 to T2
*              	: mat1        	- first barrier date   
*              	: mat2        	- maturity of option (last barrier date)   
*              	: disc1         - discount factor at T1	  
*              	: disc2         - discount factor at T2	  
*              	: below_above    	- 0/1 - barrier down or up	  
*                                                                         
* RETURNS      	: premium       - of option        
*               : deltas,
*		: gamma,
*		: vegas,
*		: theta                                                           
*******************************************************************************/

double 	srt_f_opttpamdi  	(
				double 	fwd1, 
				double 	fwd2, 
				double 	barrier, 
				double 	vol1, 
				double 	vol2, 
				double 	mat1, 
				double 	mat2, 
				double 	disc1, 
				double 	disc2, 
				SrtBarrierType	below_above,
				SrtGreekType 	greek) 

{ /* BEGIN srt_f_opttpamdi() */
	double	answer ;
	double	premium;
	double	shift;

premium = srt_f_opttpamdi_price (
				 fwd1, 
				 fwd2, 
				 barrier, 
				 vol1, 
				 vol2, 
				 mat1, 
				 mat2, 
				 disc1, 
				 disc2, 
				 below_above ) ;
        
switch (greek)
{
	case PREMIUM : 	/*** PREMIUM ***/  
		answer = premium;
		break;

	case DELTA_FWD : 	/*** DELTA FWD ***/  
		shift = fwd2 / 10000;
		answer = ( srt_f_opttpamdi(	fwd1 ,
									fwd2 + shift,
									barrier,
									vol1,
									vol2,
									mat1,
									mat2,
									disc1,
									disc2,
									below_above,
									PREMIUM)
				- premium ) /shift;      		 
		break;

	case DELTA :	/*** DELTA FWD1 + FWD2 ***/  
		shift = fwd1 / 10000;
		answer = ( srt_f_opttpamdi(	fwd1 * (1 + shift/fwd1),
									fwd2 * (1 + shift/fwd2),
									barrier,
									vol1,
									vol2,
									mat1,
									mat2,
									disc1,
									disc2,
									below_above,
									PREMIUM)
				- premium ) /shift;      		 
		break;

	case GAMMA :	/*** GAMMA ***/
		shift = fwd1 / 10000;
		answer = srt_f_opttpamdi(	fwd1 * (1 + shift/fwd1),
									fwd2 * (1 + shift/fwd2),
									barrier,
									vol1,
									vol2,
									mat1,
									mat2,
									disc1,
									disc2,
									below_above,
									PREMIUM );
		answer += srt_f_opttpamdi(	fwd1 * (1 + shift/fwd1),
									fwd2 * (1 + shift/fwd2),
									barrier,
									vol1,
									vol2,
									mat1,
									mat2,
									disc1,
									disc2,
									below_above,
									PREMIUM);
		answer -= 2 * premium;
		answer /= shift * shift;
		break ;

	case VEGA: 	/*** VEGA ***/
		shift = GVOPT.vol_add;
		answer = ( srt_f_opttpamdi(	fwd1 ,
									fwd2 ,
									barrier,
									vol1 + shift,
									vol2 + shift,
									mat1,
									mat2,
									disc1,
									disc2,
									below_above,
									PREMIUM)
				- premium )  / shift;
		break;

	case THETA:		/*** THETA  ***/
		shift = YEARS_IN_DAY;
		answer = srt_f_opttpamdi(	fwd1 ,
									fwd2 ,
									barrier,
									vol1 ,
									vol2 ,
									mat1 - shift,
									mat2 - shift,
									disc1
										* exp(-shift*log(disc1)/mat1),
									disc2
										* exp(-shift*log(disc2)/mat2),
									below_above,
									PREMIUM)
				- premium ;
		break;

	default:
		answer = UNKNOWN_GREEK;
		break;
}

return(answer);

} /* END   srt_f_opttpamdi() */ 
/************************************************************************/


static double 	srt_f_opttpamdi_price  	(
				double 			fwd1, 
				double 			fwd2, 
				double 			barrier, 
				double 			vol1, 
				double 			vol2, 
				double 			mat1, 
				double 			mat2, 
				double 			disc1, 
				double 			disc2, 
				SrtBarrierType 	below_above)
{ /* BEGIN srt_f_opttpamdi_price */
double	premium=0;
double 	measure=0;
double 	mu12, mu12_p, r, 
		b , sig1 , sig1_sqrt, sig2, sig2_sqrt ,
		alpha, a0, b0, a1, b1 ;
int	du;

	/* du = 1 if DOWN  ;   du = -1 if UP */
    du 	= ((below_above == SRT_DOWN) ? 1 : -1);

 	if ( mat1 > mat2 )  
		return(0);
	else
	if (mat1 == mat2)
	{
		if (du * (fwd1 - barrier) < 0)
			return (0);
		else
			return (1);
	}
	
	sig1 = vol1;
	sig2_sqrt = vol2*vol2*mat2 - vol1*vol1*mat1 ;
	if (sig2_sqrt <0 ) 
		return 0;
	sig2_sqrt = sqrt(sig2_sqrt);
    sig2	  = sig2_sqrt/sqrt(mat2-mat1);
                   
	b	= log( barrier / fwd1 ) + sig1 * sig1 * mat1/ 2;

	if (disc1 < disc2)
		return 0;
	r 	= -log(disc2/disc1)/(mat2-mat1);

	mu12	= log(fwd2/fwd1)/(mat2-mat1) - sig2*sig2 / 2 ;
	mu12_p = sqrt( 2*r*sig2*sig2 + mu12*mu12 );     

	sig1	  = vol1 ;
	sig1_sqrt = vol1 * sqrt( mat1 );

	if (mat1 > 0 )
		premium = norm( du * b / sig1_sqrt );
	else
	{
        if ( du * b > 0 )
			premium = 1;
		else
			premium = 0.0;
	}

	alpha 	= du * (mu12_p - mu12) / (sig2*sig2);
	a0 	= -du * mu12_p * (mat2-mat1);
	b0 	= sig2_sqrt;
	a1 	= - du * b;
	b1 	= sig1_sqrt;

	premium += srt_f_intnrm_j2(alpha, a0, b0, a1, b1);

	alpha 	= -du * (mu12_p + mu12) / (sig2*sig2);
	a0 	= du * mu12_p * (mat2-mat1);
	b0 	= sig2_sqrt;
	a1 	= - du * b;
	b1 	= sig1_sqrt;

	premium += srt_f_intnrm_j2(alpha, a0, b0, a1, b1);

	premium *= disc1;

	return (premium);

} /* END   srt_f_optamdi_price */ 


