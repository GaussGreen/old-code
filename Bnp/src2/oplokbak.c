/*******************************************************************************
**                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optlokbak()                               	  
*                                                                         
* PURPOSE     	: lookback options for stocks, bonds & swaps	  
*                                                                         
* DESCRIPTION  	: XX 
*
* CALLS		: norm() in gen_math.c                            
*                                                                         
* PARAMETERS  	: fwd_price 	- forward price of underlying        	  
*             	: spot     	- spot price of underlying                
*             	: hist_extr   	- historical extrema (for hedging purposes)   
*             	: vol         	- volatility of underlying                
*             	: mat        	- maturity of option                  
*             	: disc       	- discount factor	                  
*             	: call_put      - type of option: 0 call, 1 put      	  
*             	: greek   	- premium or greek wanted       	  
*                                                                         
* RETURNS      	: premium    	- price of the option   
*                                                                         
*******************************************************************************/
     
double 	srt_f_optlokbak	(	
			double 	fwd_price, 
			double 	spot, 
			double 	hist_extr, 
			double 	vol, 
			double 	mat, 
			double 	disc, 
			SrtCallPutType 	call_put, 
			SrtGreekType	greek
			)
{  /* BEGIN srt_f_optlokbak */
                                

/****   LOOK_BACK payoff :  				      	 	     					****/
/****		CALLS (call_put=0) 					    		   	   			****/
/**** 			S(T) - min S(t) 	       	       	   			****/
/**** 		PUTS (call_put=1) 				       	       	  					****/
/**** 			Max S(t) - S(T)  (t<T)	       	       				****/
/****  this function requires the historical extremum between start and t 	****/
	


	
	double answer, shift;
	double l_ext; 		/** = log(MIN..MAX /SPOT) **/
	double mu;	       	/** dS/S = mu dt + vol dW         **/
	double d_ext; 		/** used in BS call_put formulae **/
	double var;			/** = vol*vol*mat		  **/	
	double nu;			/** = 2*mu/(vol*vol)		  **/

	double premium=0.0;

/** 		call_put   	0=CALL	 	1=PUT
	----------------------------------------------
		0=STD		S-min		Max-S
	----------------------------------------------
      
**/                                                             

/**IN ORER TO PREVENT ERRORS IN EXTREMA **/
	
if ( ( call_put == SRT_CALL ) && ( hist_extr > spot) )
	hist_extr = spot;
if ( ( call_put == SRT_PUT  ) && ( hist_extr < spot) )
	hist_extr = spot;
	
/** DEFINE VARIABLES USED IN PRICES **/
if ( (mat != 0) && (vol != 0) )
{ 
	l_ext = log( hist_extr/spot );

	mu = log( fwd_price/spot ) / mat ;
	if (fabs(mu)<0.000001)
		mu = 0;
	
	var = vol*vol*mat;
	d_ext = ( -l_ext + mu*mat + 0.5*var ) / sqrt(var);
	nu = 2*mu / (vol*vol);
}

if ( call_put == SRT_CALL ) 
{ 
	/** CALL S(T) - Min S(s)  **/
	if ( mat <= 0 )
	{
		premium = spot-hist_extr;
	}
	else
	if ( vol == 0 )
	{
		if (fwd_price > hist_extr)
		{
			premium = fwd_price-hist_extr;
			premium *= disc;
		}
		else 
			premium = 0.0;
	}
	else
    if ( mu != 0 ) 
	{
		premium =  exp(nu*l_ext) * norm(-d_ext+ 2*mu*sqrt(mat)/vol);
		premium -=  exp( mu*mat ) * norm( -d_ext );
		premium *= spot / nu;
	}
	else
	{
		premium = spot*sqrt(var)*( -d_ext*norm(-d_ext)+gauss(-d_ext) );
	
		premium += fwd_price * norm( d_ext);
		premium -= hist_extr * norm( d_ext - sqrt(var) );
                   
  		premium *= disc;
	}	
}
else       
if ( call_put == SRT_PUT ) 
{
	/** PUT Max S(s) - S(T)     **/
    if ( mat <= 0 )
	{	
		premium = hist_extr-spot;
	}
	else
	if ( vol == 0 )
	{
		if (fwd_price< hist_extr)
		{
			premium = hist_extr-fwd_price;
			premium *= disc;
		}
		else 
			premium = 0.0;
	}
	else
	if ( mu != 0 ) 
	{
		premium = -exp(nu*l_ext)*norm(d_ext-2*mu*sqrt(mat)/vol);
		premium +=  exp( mu*mat ) * norm( d_ext );
		premium *= spot / nu;
	}
	else
	{
		premium = spot*sqrt(var)*( d_ext*norm(d_ext)+gauss(d_ext) ); 
   		premium -= fwd_price * norm( -d_ext);
		premium += hist_extr * norm( -d_ext + sqrt(var) );   
		premium *= disc;	
	}
}


if ( greek == PREMIUM)      
{
	answer = premium; 
}
else if ( greek == DELTA_FWD )
{		                        /* DELTA FWD UP */
	shift = fwd_price/ 10000;
	answer = ( srt_f_optlokbak(	fwd_price+shift,
								spot,
								hist_extr,
								vol,
								mat,
								disc,
								call_put,
								PREMIUM) 
			- premium ) / shift;  
}         
else if ( greek == DELTA )
{                                  /* DELTA */
	shift = spot/ 10000;
	answer = ( srt_f_optlokbak(	fwd_price * ( 1 + shift/spot),
								spot+shift,
								hist_extr,
								vol,
								mat,
								disc,
								call_put,
								PREMIUM) 
			- premium ) / shift;  
}
/*
else if ( greek == GAMMA_FWD )
{                                  
	shift = fwd_price/ 10000;
	answer = ( srt_f_optlokbak(	fwd_price+shift,
					spot,
					hist_extr,
					vol,
					mat,
					disc,
					call_put,
					PREMIUM)
		+ srt_f_optlokbak(	fwd_price-shift,
					spot,
					hist_extr,
					vol,
					mat,
					disc,
					call_put,
					PREMIUM) 
		- 2*premium) / (shift*shift); 
}
*/
else if ( greek == GAMMA )
{                                  /*   GAMMA SPOT */
	shift = spot/ 10000;
	answer = ( srt_f_optlokbak(	fwd_price * ( 1 + shift/spot),
								spot + shift,
								hist_extr,
								vol,
								mat,
								disc,
								call_put,
								PREMIUM)
			+ srt_f_optlokbak(	fwd_price * ( 1 - shift/spot),
								spot-shift,
								hist_extr,
								vol,
								mat,
								disc,
				       	 		call_put,
								PREMIUM) 
			- 2*premium) / (shift*shift); 
}
else if ( greek == VEGA )
{                                  /*   VEGA    */
	shift = GVOPT.vol_add;
	answer = ( srt_f_optlokbak(	fwd_price,
								spot,
								hist_extr,
								vol + shift,
								mat,
								disc,
								call_put,
								PREMIUM) 
			- premium ) / shift;  
}          
else if ( greek == THETA )
{                                  /*   THETA   */
	shift = YEARS_IN_DAY;
	answer =  ( srt_f_optlokbak(fwd_price,
								spot,
								hist_extr,
								vol,
		     					mat - shift,
								disc
									* exp(-shift*log(disc)/mat),
								call_put,
								PREMIUM)  
		     - premium);  	
}              
else
	answer = UNKNOWN_GREEK;

return(answer);   

} /* END srt_f_optlokbak() */         

/******************************************************************************/
