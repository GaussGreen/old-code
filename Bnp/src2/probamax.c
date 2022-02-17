/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include        "utallhdr.h"
#include        <num_h_proba.h"
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                         
* FUNCTION     	: proba_max_fct(...)                                   
*                                                                                                                                            
* DESCRIPTION  	: Returns the probability to get 1 unit of currency at maturity
*                 if the spot stays under the barrier during mat.                    
*                                                                   
* CALLS			: gauss(...) in gen_math.c
*				  norm(...)  in gen_math.c	      
*		                                    
* PARAMETERS                                                           
*	INPUT		: spot			- spot underlying price    	 
*				: barrier      	- barrier level to stay under		      	  
*				: vol         	- annual volatility   	      	  
*				: mat         	- maturity, in years  	      	       
*				: disc        	- discount factor for the domestic currency   		 
*                                                                        
* RETURNS		: proba			- probability		      	 
*                                                                       
*******************************************************************************/

double proba_max_fct(double spot,
					 double barrier,
					 double vol, 
					 double mat,
					 double disc)
{
double mu,u,d1,d2,price = 0.0,coeff;

	if (spot>=barrier) return(0);
		
	mu = -log(disc) / mat - vol*vol/2;
	
	u=log(barrier/spot);
	
	coeff=exp(2*mu*u/vol/vol);

	d1=(u-mu*mat)/vol/sqrt(mat);
	d2=(-u-mu*mat)/vol/sqrt(mat);

	price = (norm(d1)-coeff*norm(d2));
	price *= disc;

	return(price);
}
