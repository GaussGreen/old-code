
/*****************************************************************************
*                      Include Files                               
*******************************************************************************/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/******************************************************************************/

 
/******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_opttriscud                                        
*                                                                         
* PURPOSE      	: DOUBLE BARRIER SCUD Option -BARRIER ON SECONDARY UND      
*
* DESCRIPTION  	: Barrier on a secondary underlying                        
*                                                                         
* AUTHORS	: O Van Eyseren & M Baker                        
*******************************************************************************/

static double  prob_tri_scud(      
		double 		y, 
		double 		L, 
		double 		U, 
		double 		A, 
		double 		mux,
		double		muy,
		double 		sigx,
		double 		sigy,
		double 		rho,
		double 		mat,   
		int		du,
		int		cp,
		int		nterms ) ;


static double 	triscud_price	(
			double 	fwdy, 
			double 	fwdx, 
			double 	spoty, 
			double 	spotx, 
			double 	strikey,
			double 	barriery_up, 
			double 	barriery_dwn, 
			double 	barrierx, 
			double 	sigy, 
			double 	sigx, 
			double 	rho, 
			double 	mat,
   			double 	disc, 
			SrtCallPutType 	call_put,     
			int 	scud_type,   
			int	nterms ) ;

             
double srt_f_opttriscud( double fwdy,    
			 double fwdx,    
			 double spoty,    
			 double spotx,    
			 double strikey,     
			 double barriery_up,  
			 double barriery_dwn, 
			 double barrierx,   
			 double sigy,       
			 double sigx,       
			 double rho,        
			 double mat,        
			 double disc,    
			 SrtCallPutType	call_put,    	    /*ca_pu in y */ 
			 int	scud_type,	    /*up_do in x */    
			 int	nterms,   
			 SrtGreekType	greek)   
{ 
        double  f0, f1, f2, f3;
	double	temp, temp1, dx, dx1;

   	double  answer;

{
	f0  =  triscud_price(fwdy,
				fwdx,
				spoty,
				spotx,
				strikey,
				barriery_up,
				barriery_dwn,
				barrierx,
				sigy,
				sigx,
				rho,
				mat,
				disc,
				call_put,
				scud_type,
				nterms );  }

/* greek 	= 0 	premium
		= 1     fwd delta y
		= 2     fwd delta x
		= 3 	spot delta y
		= 4     spot delta x
		= 5     fwd gamma y
		= 6     fwd gamma x    
		= 7     cross fwd gamma
		= 8     vega y 
		= 9     vega x 
		= 10    cross vega - change in price wrt correlation 
		= 11    theta 
		= 21    total delta y 
		= 22    total delta x  */

if (greek != PREMIUM)   
{
		if ((greek==DELTAY) || (greek==GAMMAY))
		{       
			temp 	= fwdy;			/* fwd delta y */ 
			dx   	= .01 *temp;           	/* fwd gamma y */
			fwdy 	= temp + dx;            /* cross gamma */
			spoty 	= temp+ dx; 
		}
		if ((greek==DELTAX) || (greek==GAMMAX))
		{       
			temp 	= fwdx;                 /* fwd delta x */ 
			dx   	= .01 *temp;  
			fwdx 	= temp+ dx; 
			spotx	= temp+ dx; 
		}
		if (greek==GAMMAXY)
		{       
			temp 	= fwdy ;                 
			dx   	= .001 *temp;  
			fwdy 	= temp + dx;  
			spoty 	= temp+ dx; 
	      		temp1 	= fwdx;                 
			dx1   	= .001 *temp1;           
    			fwdx 	= temp1 + dx1;           
			spotx	= temp+ dx; 
		}
		if (greek==VEGAY)
		{       
			temp 	= sigy ;                /* vega y */   
			dx   	= .01   ;  
			sigy 	= temp + dx ; 
		}
		if (greek==VEGAX)
		{       
			temp 	= sigx ;                /* vega x */   
			dx   	= .01  ;  
			sigx 	= temp + dx ; 
		}
		if (greek==VEGAXY)
		{       
			temp 	= rho ;                 /* cross vega */   
			dx   	= .01 *temp ;  
			rho 	= temp+ dx ; 
		}
		if (greek==THETA)
		{       
			temp 	= mat ;                 /* theta */   
			dx   	= - YEARS_IN_DAY;  
			mat 	= temp + dx ; 
		}
              
              
		f1 =   triscud_price(fwdy,
				fwdx,
				spoty,
				spotx,
				strikey,
				barriery_up,
				barriery_dwn,
				barrierx,
				sigy,
				sigx,
				rho,
				mat,
				disc,
				call_put,
				scud_type,
				nterms );  
               
		if ((greek==DELTAY) || (greek==GAMMAY))
		{      	
			fwdy 	= temp - dx ; 
			spoty 	= temp - dx ; 
		}

		if ((greek==DELTAX) || (greek==GAMMAX))
		{                                    
			fwdx 	= temp - dx ; 
		 	spotx 	= temp - dx ; 
		}
                                         

		if (greek==GAMMAXY)
		{	
			fwdy  	=  temp + dx ;   
			fwdx  	=  temp1 ;   
			spoty 	= temp - dx ; 
		 	spotx 	= temp - dx ; 
		}
                                                                
		if (greek==VEGAY)
		{	
			sigy 	= temp - dx ;
			dx   	= 1.0 ; 
		}

		if (greek==VEGAX)
		{	
			sigx 	= temp - dx ;
			dx   	= 1.0 ; 
		}

		if (greek==VEGAXY)
		{	
			rho 	= temp - dx ;
			dx  	= 1.0 ; 
		}

		if (greek==THETA)
		{ 	
			mat 	= temp ;
			dx  	= 2.0 ; 
		}


		f2 =    triscud_price(fwdy,
				fwdx,
				spoty,
				spotx,
				strikey,
				barriery_up,
				barriery_dwn,
				barrierx,
				sigy,
				sigx,
				rho,
				mat,
				disc,
				call_put,
				scud_type,
				nterms );      
	}

	if (greek==PREMIUM) 
	{
		answer = f0; 
	}
	else if ((greek==DELTAX) || (greek==DELTAY))      
	{
		answer = 0.5 * ( f1 - f2 ) / dx; 
	}
	else if ((greek==GAMMAY) || (greek==GAMMAX)) 
	{
		answer =  ( f1 - 2.0 * f0 + f2 ) / (dx*dx); 
	}
	else if (greek==THETA) 
	{
		answer =  f1 - f2; 
	}
	else if (greek==GAMMAXY) 
	{                                   
		fwdy = temp ;
		fwdx = temp1 + dx1 ;

		f3  =  triscud_price(fwdy,
				fwdx,
				spoty,
		 		spotx,
				strikey,
				barriery_up,
			 	barriery_dwn,
				barrierx,
				sigy,
				sigx,
				rho,
				mat,
				disc,
				call_put,
				scud_type,
				nterms );

		answer = ( ( f1 - f2 ) - ( f3 - f0 ) )/ (dx* dx1) ;
	}

	return( answer );

} /* END srt_f_opttriscud() */     


/****************************************************************************/
static double triscud_price (
			double 	fwdy, 
			double 	fwdx, 
			double 	spoty, 
			double 	spotx, 
			double 	strikey,
			double 	barriery_up, 
			double 	barriery_dwn, 
			double 	barrierx, 
			double 	sigy, 
			double 	sigx, 
			double 	rho, 
			double 	mat,
   			double 	disc, 
			SrtCallPutType 	call_put,     /*ca_pu in y */
			int 	scud_type,    /*up_do in x */
			int	nterms )
{ 
	double 	mux;
	double	muy;

	double 	l_yu;
	double 	l_yd;
	double 	l_x;
	double 	l_k;

	double	premium;

	int 	cp;
	int	du;

/*	note double barrier & payoff on  y - extinguish on x 		*/
/*		notation has reversed from the call routine  		*/
	 
mux   	= log( fwdx / spotx ) / mat - ( sigx * sigx / 2);
muy   	= log( fwdy / spoty ) / mat - ( sigy * sigy / 2);

l_yu    = log( barriery_up  / spoty );
l_yd    = log( barriery_dwn / spoty );
l_x     = log( barrierx / spotx );

l_k     = log( strikey / spoty );

cp    	= ( call_put == SRT_CALL ) ? 1 : -1 ;     /* 1 for a call */
du    	= 1 - ( 2 * scud_type );       /* 1 for down   */  

if (l_x * du > 0) return 0;
if (l_yu < 0 ) return 0;
if (l_yd > 0 ) return 0;

if ( cp==1 ) 
{ 	if ( l_k > l_yu ) {	 return 0 ; } 
}
else if ( cp == -1 ) 
{ 	if ( l_k < l_yd ) {	 return 0 ; } 
}
   
premium =  fwdy 
	*  prob_tri_scud ( 
			  l_k, 
			  l_yd, l_yu,
			  l_x, 
			  mux+rho*sigx*sigy, muy+sigy*sigy, 
			  sigx, sigy, 
			  rho, 
			  mat,
			  du,
			  cp,
			  nterms ); 
                       
premium -=  strikey
	 * prob_tri_scud ( 
			  l_k, 
			  l_yd, l_yu, 
			  l_x, 
			  mux, muy, 
			  sigx, sigy, 
			  rho, 
			  mat,
			  du,
			  cp,
			  nterms); 

        premium *= disc;
        premium *= cp;

	return premium;
} /* END triscud_price() */     


/****************************************************************************/
static double  prob_tri_scud(      
		double 		y, 
		double 		L, 
		double 		U, 
		double 		A, 
		double 		mux,
		double		muy,
		double 		sigx,
		double 		sigy,
		double 		rho,
		double 		mat,   
		int		du,
		int		cp,
		int		nterms )
              
{

	double 		result = 0.0;
	double 		B;
	double 		bn; 
	double		xtest, term1, term2 ;
	double 		sigxrt,
			sigyrt,
			sigx_2,
			sigy_2 ;
        int n;

/*	payoff in y dble barrier in y;  scud one barrier in x	*/
	sigx_2 = sigx * sigx ;
	sigy_2 = sigy * sigy ;

	sigxrt = sigx * sqrt(mat);
	sigyrt = sigy * sqrt(mat);

	B = (cp == 1) ?  U : L;	

	for( n= - nterms; n<=nterms; n++)
	{ 

        bn	= 2.0 * n * (U-L) ;

	term1 = ( bivar( du*(-A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			 cp* (-y + bn + muy*mat) /sigyrt,
			 cp * du *rho )
		- bivar( du*(-A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			 cp *(-B + bn + muy*mat) /sigyrt,
			 cp * du * rho ) ) ;



	xtest =  2*mux*A/sigx_2 + 2*A*bn*rho/(sigx*sigy*mat) ;

	term1 -= exp( xtest ) *( 
		 bivar( du * ( A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
		       cp * (-y + bn + 2*rho*sigy*A/sigx + muy*mat) /sigyrt, 
		       du * cp * rho )
		-bivar( du * ( A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			cp * (-B + bn +2*rho*sigy*A/sigx +muy*mat) /sigyrt, 
			cp * du * rho )	 );

	result +=  exp(muy*bn /sigy_2) * term1 ;

        bn	= 2.0 * U - bn ;

	term2 = ( bivar( du*(-A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			 cp* (-y + bn + muy*mat) /sigyrt,
			 cp * du *rho )
		- bivar( du*(-A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			 cp *(-B + bn + muy*mat) /sigyrt,
			 cp * du * rho )	 ) ; 


	xtest = 2*mux*A/sigx_2 + 2*A*bn*rho/(sigx*sigy*mat) ;


	term2 -= exp( xtest ) *(
		 bivar( du * ( A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			cp * (-y + bn + 2*rho*sigy*A/sigx + muy*mat) /sigyrt, 
			du * cp * rho )
		 -bivar( du * ( A + rho*sigx*bn/sigy + mux*mat) /sigxrt,
			 cp * (-B + bn + 2*rho*sigy*A/sigx + muy*mat) /sigyrt, 
			 cp * du * rho )  	);

	result -= exp(muy*bn /sigy_2) * term2 ;   

        }
	return result;
}
