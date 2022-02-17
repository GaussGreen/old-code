/*------------------------------------------------------------------
PURPOSE:       Computation of option premium 
		where option has payoff on a variable 
		and is subject
	       to a corridor condition on another.
AUTHOR:        Bing-Le Wu
-------------------------------------------------------------------------*/
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <math.h>

#include "gtobf.h"
#include "binorm.h"


#define ERROR	(0.0001)


/*f--------------------------------------------------------------
 * Options formulae : double KO bivariate.
 *
 * <br><br>
 * Computes the premium of bivariate double KO option.
 * Written by Bing Le Wu (aka "Da Bing Machine").
 */

DLL_EXPORT(int)
DrlOptionDoubleKO(
	double expiration,	/* (I) time to expiration */
	double forward2,	/* (I) payoff parameters*/
	double strike,
	double vol2,
	double spot1,		/* (I) barrier variable parameters*/
	double forward1,
	double vol1,
	double upperBarrier,
	double lowerBarrier,
	double correlation,	/* (I) correlation between var1 and var2*/
	char   optionType,	/* (I) C for call and P for put*/
	double *premium)	/* (O) output*/

{  
	double A=0;
	double B=0;
	double d[10];
	double n[8];
	double temp;
	int    type;
	int    k;
	int    N;
	
	
		
	if(optionType=='C'|| optionType=='c')
		type=1;
	if(optionType=='P'|| optionType=='p')
		type=-1;

	temp=2*log(upperBarrier/lowerBarrier)*log(upperBarrier/lowerBarrier)
				/(vol1*vol1*expiration);
	N=0;
	while(exp(-N*N*temp)*8*forward2>ERROR)
		N++;
	N=N+10;       
	d[0]=(log(lowerBarrier/forward1)+0.5*vol1*vol1*expiration
			-correlation*vol1*vol2*expiration)*
					1/(vol1*sqrt(expiration));
	d[1]=(log(upperBarrier/forward1)+0.5*vol1*vol1*expiration
			-correlation*vol1*vol2*expiration)*  
					1/(vol1*sqrt(expiration));    
	d[2]=log(upperBarrier/lowerBarrier)/(vol1*sqrt(expiration));
	d[3]=(log(forward2/strike)+0.5*vol2*vol2*expiration)
			/(vol2*sqrt(expiration));
	d[4]=(log(lowerBarrier/forward1)+0.5*vol1*vol1*expiration)
			/(vol1*sqrt(expiration));
	d[5]=(log(upperBarrier/forward1)+0.5*vol1*vol1*expiration)
			/(vol1*sqrt(expiration));   
	d[6]=(log(forward2/strike)-0.5*vol2*vol2*expiration)
			/(vol2*sqrt(expiration));
	d[7]=log(lowerBarrier/spot1)/(vol1*sqrt(expiration));
	d[8]=2/(vol1*vol1)*(log(forward1/spot1)/expiration- 0.5*vol1*vol1);
	d[9]=2*correlation*vol2/vol1;
	
	for(k=0; k<N+1;  k++) 
	{
	   GtoBiNormalCum(d[0]+2*k*d[2], type*(d[3]-2*k*correlation*d[2]), 
			-type*correlation, &n[0]);
	   GtoBiNormalCum(d[1]+2*k*d[2], type*(d[3]-2*k*correlation*d[2]), 
			-type*correlation, &n[1]);
	   GtoBiNormalCum(d[4]+2*k*d[2], type*(d[6]-2*k*correlation*d[2]), 
			-type*correlation, &n[2]);
	   GtoBiNormalCum(d[5]+2*k*d[2], type*(d[6]-2*k*correlation*d[2]), 
			-type*correlation, &n[3]);

	   GtoBiNormalCum(d[0]+2*k*d[2]-2*d[7], 
			type*(d[3]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[4]);
	   GtoBiNormalCum(d[1]+2*k*d[2]-2*d[7], 
			type*(d[3]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[5]);
	   GtoBiNormalCum(d[4]+2*k*d[2]-2*d[7], 
			type*(d[6]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[6]);
	   GtoBiNormalCum(d[5]+2*k*d[2]-2*d[7], 
			type*(d[6]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[7]);


	   A+=type*pow(lowerBarrier/upperBarrier, k*d[8])*
			(forward2*pow(lowerBarrier/upperBarrier,k*d[9])*(n[1]-n[0])
			-strike*(n[3]-n[2]));
	   B+=type*pow(pow(lowerBarrier, k+1)/(pow(upperBarrier, k)*spot1), d[8])*
		(forward2*pow(pow(lowerBarrier, k+1)/(pow(upperBarrier, k)*spot1), d[9])*(n[5]-n[4])-strike*(n[7]-n[6])); 
					
	}

	for(k=-1; k>-N-1;  k--) 
	{
	   if((d[0]+2*k*d[2]<-7.0) &&
		(d[1]+2*k*d[2]<-7.0) &&
		(d[4]+2*k*d[2]<-7.0) &&
		(d[5]+2*k*d[2]<-7.0) )
		A+=0.0;
	   else{
		GtoBiNormalCum(d[0]+2*k*d[2], 
			type*(d[3]-2*k*correlation*d[2]), 
				-type*correlation, &n[0]);
		GtoBiNormalCum(d[1]+2*k*d[2], 
			type*(d[3]-2*k*correlation*d[2]), 
				-type*correlation, &n[1]);
		GtoBiNormalCum(d[4]+2*k*d[2], 
			type*(d[6]-2*k*correlation*d[2]), 
				-type*correlation, &n[2]);
		GtoBiNormalCum(d[5]+2*k*d[2], 
			type*(d[6]-2*k*correlation*d[2]), 
				-type*correlation, &n[3]);
		A+=type*pow(lowerBarrier/upperBarrier, k*d[8])*
			(forward2*pow(lowerBarrier/upperBarrier,k*d[9])*(n[1]-n[0])
				-strike*(n[3]-n[2]));
	   }
	   if((d[0]+2*k*d[2]-2*d[7]<-7.0) &&
		(d[1]+2*k*d[2]-2*d[7]<-7.0) &&
		(d[4]+2*k*d[2]-2*d[7]<-7.0) &&
		(d[5]+2*k*d[2]-2*d[7]<-7.0) )
		B+=0.0;
	   else{
		GtoBiNormalCum(d[0]+2*k*d[2]-2*d[7], 
			type*(d[3]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[4]);
		GtoBiNormalCum(d[1]+2*k*d[2]-2*d[7], 
			type*(d[3]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[5]);
		GtoBiNormalCum(d[4]+2*k*d[2]-2*d[7], 
			type*(d[6]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[6]);
		GtoBiNormalCum(d[5]+2*k*d[2]-2*d[7], 
			type*(d[6]-2*correlation*(k*d[2]-d[7])), 
				-type*correlation, &n[7]);
		B+=type*pow(pow(lowerBarrier, k+1)/(pow(upperBarrier, k)*spot1), d[8])*
		 (forward2*pow(pow(lowerBarrier, k+1)/(pow(upperBarrier, k)*spot1), d[9])*(n[5]-n[4])-strike*(n[7]-n[6]));
	  }
	}
	if(A-B<ERROR/10.0)
		*premium=0;     
	else
		*premium=A-B;

	return(SUCCESS);
}
