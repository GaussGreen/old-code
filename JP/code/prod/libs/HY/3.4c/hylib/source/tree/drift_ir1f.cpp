/****************************************************************************/
/*	
/*	drift1.cpp
/*
/****************************************************************************/
#include "tree_ir1f.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>



/*****  FindDrift   ********************************************************/
/*
*       Calculate the interest rate drift. 
*/
void	IRTree1F::drift(void)
{
	double	jump,expJump;     
	double	d, du, grid;       
	double	qu, qd, q0, Qi;    
	double	p, temp;
	double	zeroPrice;
	double	timeTillTp;
	double	timeTillNextTp;
	double	P[3];
	int		dim, nextDim; 
	int		i, m;
	int		iu, i0, id;
	double	beta = _process.get_mean_reversion()[0];

	/* the first drift*/
	KTimeLine::iterator	iter, pIter, nIter;

	/*find the last node, and get the max memory needed */
	iter = tl_end();
	iter--;
	TimeSlice1D<double>	statePr(iter->second[__irTreeMaxLimit1F]);
	TimeSlice1D<double>	newStatePr(iter->second[__irTreeMaxLimit1F]);

	
	iter = tl_begin();
	iter->second[__irTreeDrift] = iter->second[__fwdRateTillNextTp];
	/* update state price for the second node*/
	iter++;
	zeroPrice = iter->second[__zeroPrice];
	statePr[1] = 1./(2*K_JUMPCOEFF) * zeroPrice;
	statePr[0] = (1 - 1./K_JUMPCOEFF)*zeroPrice;
	statePr[-1] = 1./(2*K_JUMPCOEFF) * zeroPrice;
	
	while(iter != tl_end())
	{									
		nIter = pIter = iter;
		pIter--;
		nIter++;
		if(nIter == tl_end())
			break;
		/*Calibrate drift*/
		timeTillTp = pIter->second[__timeTillNextTp];
		timeTillNextTp = iter->second[__timeTillNextTp];

		du = sqrt (K_JUMPCOEFF * timeTillTp);      
		jump = pIter->second[__irTreeSpotVol1F] * du;   
		expJump = exp (jump);

      	P[0] = P[1] = P[2] = 0.;
		dim = iter->second[__irTreeLimit1F];
     	grid = pIter->second[__irTreeDrift] * iter->second[__fwdRateTillNextTp]
				/pIter->second[__fwdRateTillNextTp] * exp (-dim*jump);
		for (i = -dim; i <= dim; i++)	
		{
			temp =  statePr[i]/(1. + grid);
			P[0] += temp;
			temp *= grid / (1. + grid);
			P[1] += temp;
			temp *= grid / (1. + grid);
			P[2] += temp;
			grid *= expJump;
		}   
		zeroPrice = nIter->second[__zeroPrice];
		P[0] -= zeroPrice;  
		iter->second[__irTreeDrift] =  pIter->second[__irTreeDrift] * iter->second[__fwdRateTillNextTp]
										/pIter->second[__fwdRateTillNextTp]
										*(1. - (-P[1] + sqrt(P[1]*P[1] - 4* P[0]*P[2]))/(2.*P[2]));								


		/*Update State Prices*/
		nextDim = nIter->second[__irTreeLimit1F];
		for (i = -nextDim; i <= nextDim; i++)	
			newStatePr[i] = 0.;

		d  = (1. - beta * timeTillNextTp)  
				 *pIter->second[__irTreeSpotVol1F] * sqrt(timeTillTp)
				 /(iter->second[__irTreeSpotVol1F] * sqrt(timeTillNextTp)) - 1.;
	
		grid = iter->second[__irTreeDrift] * exp (-dim*jump);
		for (i = -dim; i <= dim; i++)		
		{
			Qi = d * i;
			if(Qi > 0)
				m = (int)(Qi + 0.5);
			else		
				m = (int)(Qi - 0.5);
			Qi -= m;
			qu = .5 * (1./K_JUMPCOEFF + Qi + Qi * Qi);
			qd = qu - Qi;
			q0 = 1. - qu - qd;

            temp = statePr[i] /(1. + grid);
			iu = Min(nextDim, Max(-nextDim, i+1+m));
			i0 = Min(nextDim, Max(-nextDim, i+m  ));
			id = Min(nextDim, Max(-nextDim, i-1+m));
			newStatePr[iu] += qu * temp;  
			newStatePr[i0] += q0 * temp; 	
			newStatePr[id] += qd * temp;

    		grid *= expJump;
		}  

		p = 0.;	
		for (i =-nextDim; i <= nextDim; i++)
		{                                                               
			statePr[i] = newStatePr[i];
		  	p += statePr[i];
	    }
		p /= zeroPrice;			
		if(p < 1. - 0.00001)
			throw KException("IRTree1F::drift Failed!");
		//next iter
		iter++;
	}
}

