/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_opteurscdout                                         
*                                                                         
* PURPOSE      	: European SCUD Option                                              
*                                                                         
* DESCRIPTION  	: European Barrier on a secondary underlying                        
*                                                                         
* CALLS		: bivar()                                                      
*                                                                         
* PARAMETERS   	: fwdx     	    - forward price of 1st underlying x       
*              	: fwdy          - forward price of 2nd underlying y        
*              	: strike        - strike price			          
*              	: barrier       - barrier level 			       
*              	: sigx          - vol of 1st underlying x               	  
*              	: sigy          - vol of 2nd underlying y               	  
*              	: rho           - correlation                    		  
*              	: matx           - maturity of option, in years             
*              	: maty           - time until which Y can extinguish
				 disc          - discount factor                  	  
*              	: call_put      - type of option: 0 call, 1 put	     	  
*              	: down_up       - type of scud: ??                    	  
*              	: greek        	- ??                      		  
*                                                                         
* RETURNS      	: ??          	- ??                                      
*                                                                         
*******************************************************************************/

double 	srt_f_opteurscdout	(
			double 	         fwdx, 
			double 	         fwdy, 
			double 	         strike,
			double 	         barrier, 
			double 	         sigx, 
			double 	         sigy, 
			double 	         rho, 
			double 	         matx,
			double 	         maty,
   			double 	         disc, 
			SrtCallPutType 	 call_put,
			SrtBarrierType 	 down_up,
			SrtGreekType 	 greek
			)
{ 
	double 	mux;
	double	muy;

	double  varx;
	double  vary;
	double  stdevx;
	double  stdevy;

	double 	b;
	double	k;
	double  corr;

	double  dx;
	double  dy;
	
	double  prob_QT;
	double  prob_QX;
	
	double	premium;
	double	result;

	double 	shift;
	double	shiftx;
	double	shifty;

	int 	cp;
	int		du;

if (matx <= 0)
{
	premium = srt_f_optblksch(	fwdx,
								strike,
								sigx,
								matx,
								disc,
								call_put,
								PREMIUM);
}
else
{
	/* Relative positions of barriers and strikes (trigger levels) */
	k     = log( fwdx / strike );
	b     = log( fwdy / barrier );

	/*Variances and standard deviations */
	varx = sigx * sigx * matx;
	vary = sigy * sigy * maty;
	stdevx = sqrt(varx);
	stdevy = sqrt(vary);
	
	/* Time corrected correlation */
	corr = rho * sqrt(maty / matx);

	/* Changes in sign for call/put, up&out/down&out */
	cp    = ( call_put == SRT_CALL ) ? 1 : -1 ;
	du = ( down_up == SRT_DOWN ) ? 1 : -1;


/* Probability under QT */
	/* Drifts corrections*/
	mux   = - 0.5 * stdevx;
	muy   = - 0.5 * stdevy;

	/* Trigger levels for the normalised Gaussians */
	dx = k / stdevx + mux;
	dy = b / stdevy + muy ;

	/* Probability [ ( X <> K ) && ( Y <>B ) ] */
	prob_QT = bivar( 	cp * dx					 , 
						du * dy					 , 
						cp * du * corr 					 );

/* Probability under QX */	
	/* Drifts corrections*/
	mux   = + 0.5 * stdevx;
	muy   = - 0.5 * stdevy + rho * sigx * sqrt(maty) ;

	/* Trigger levels for the normalised Gaussians */
	dx = k / stdevx + mux;
	dy = b / stdevy + muy ;


	/* Probability [ ( X <> K ) && ( Y <>B ) ] */
	prob_QX = bivar( 	cp * dx					 , 
						du * dy					 , 
						cp * du * corr 					 );
	
	
/* The Premium */
	premium  = fwdx * prob_QX - strike * prob_QT;
	premium *= disc * cp;
}


switch ( greek ) 
{
	case PREMIUM :		
		return ( premium );        
		break;
	
	case DELTAX :
		shift  = fwdx / 10000;
		result = ( srt_f_opteurscdout(	fwdx + shift,
									fwdy,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM)
			   - premium ) / shift;
		return ( result );
		break;

	case DELTAY :
		shift  = fwdy / 10000;
		result = ( srt_f_opteurscdout(	fwdx,
									fwdy + shift,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM)
			   - premium ) / shift;
		return ( result );
		break;

	case GAMMAX :
		shift   = fwdx / 1000;
		result  = srt_f_opteurscdout(	fwdx + shift,
									fwdy,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
		result += srt_f_opteurscdout(	fwdx - shift,
									fwdy,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
        result -= 2 * premium;
		result /= shift * shift;
		return( result );
		break;

	case GAMMAY :
		shift   = fwdy / 1000; 
		result  = srt_f_opteurscdout(	fwdx,
									fwdy + shift,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
		result += srt_f_opteurscdout(	fwdx,
									fwdy - shift,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
        result -= 2 * premium;
		result /= shift * shift;
		return( result );
		break;

	case GAMMAXY :
		shifty   = fwdy / 1000; 
		shiftx   = fwdx / 1000;
		result  = srt_f_opteurscdout(	fwdx + shiftx,
									fwdy + shifty,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
		result  += srt_f_opteurscdout(	fwdx - shiftx,
									fwdy - shifty,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
		result  -= srt_f_opteurscdout(	fwdx - shiftx,
									fwdy + shifty,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
		result  -= srt_f_opteurscdout(	fwdx + shiftx,
									fwdy - shifty,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM);
		result /= 4 * shiftx * shifty;
		return( result );
		break;

	case VEGAX :
		shift = GVOPT.vol_add;
		result = ( srt_f_opteurscdout(	fwdx,
									fwdy,
									strike,
									barrier,
									sigx + shift,
									sigy,
									rho,
									maty,
									maty,
									disc,
									call_put,		
									down_up,
									PREMIUM )
				- premium ) / shift;
		return( result );
		break;

	case VEGAY :
		shift = GVOPT.vol_add;
		result = ( srt_f_opteurscdout(	fwdx,
									fwdy,
									strike,
									barrier,
									sigx,
									sigy + shift,
									rho,
									maty,
									maty,
									disc,
									call_put,
									down_up,
									PREMIUM )
				- premium ) / shift;
		return( result );
		break;

	case THETA :
		shift = YEARS_IN_DAY;
		result = srt_f_opteurscdout(	fwdx,
									fwdy,
									strike,
									barrier,
									sigx,
									sigy,
									rho,
									maty - shift,
									maty - shift,
									disc
										* exp(-shift*log(disc)/maty),
									call_put,
									down_up,
									PREMIUM )
				- premium ;
		return( result );
		break;

	default:
		return(UNKNOWN_GREEK);
		break;
	}

} /* END srt_f_opteurscdout() */

/******************************************************************************/
