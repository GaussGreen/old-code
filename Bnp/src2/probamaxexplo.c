/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include        <num_h_allhdr.h"
#include		<math.h"

/******************************************************************************/


/*******************************************************************************
*                                                                         
* FUNCTION     	: proba_max_explo_fct(...)                                   
*                                                                                                                                            
* DESCRIPTION  	: Returns the probability to get 1 unit of currency as soon as
*                 the spot becomes greater than the barrier during mat.                    
*                                                                   
* CALLS			: gauss(...) in gen_math.c
*				  norm(...)  in gen_math.c	      
*		                                    
* PARAMETERS                                                           
*	INPUT		: spot			- spot underlying price    	 
*				: barrier      	- barrier level to stay under		      	  
*				: vol         	- annual volatility   	      	  
*				: mat         	- maturity, in years  	      	       
*				: disc       	- discount factor for the domestic currency   		 
*                                                                        
* RETURNS		: proba			- probability		      	 
*                                                                       
*******************************************************************************/

double proba_max_explo_fct(double spot,							
						   double barrier,					 
						   double vol, 					 
						   double mat,					 
						   double disc)
{
double mu,mu1,mu2,u,d1,d2,price=0.0,coeff;

	if (spot>=barrier) return(1);

	mu1 = - log(disc) / mat - vol*vol/2;
	mu2 = sqrt(mu1*mu1-2*vol*vol*log(disc)/mat);
	mu=mu1-mu2;
	
	u=log(barrier/spot);
	
	coeff=exp(2*mu2*u/vol/vol);
	
	d1=(u-mu2*mat)/vol/sqrt(mat);
	d2=(-u-mu2*mat)/vol/sqrt(mat);

	price = 1-(norm(d1)-coeff*norm(d2));
	price *= exp(u*mu/vol/vol);

	return(price);
}
