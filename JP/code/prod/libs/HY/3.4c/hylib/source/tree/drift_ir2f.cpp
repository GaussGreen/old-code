/****************************************************************************/
/*      Building of the interest rate tree: calculation of forward rates,   */
/*      interest rate irDrift at each time step.          		    */
/****************************************************************************/
/*      DRIFT.c                                                             */
/****************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "tree_ir2f.h"

/*****  FindDrift   ********************************************************/
/*
*       Calculate the interest rate irDrift. 
*/
void	IRTree2F::drift()
{
	double		d[2], P[3];
	double		rho;
	double		jump[2],expJump[2];     
	double		du, grid[2];       
	double		Pi, Qi, Qij;    
	double		p, temp;
	double		pu, pd, p0;               
	double		qu, qd, q0;
	double		zeroPrice;
	double		timeTillTp;
	double		timeTillNextTp;
	int			dim1, dim2, nextDim1, nextDim2;
	int			i, j, l, m;
	int			iu, i0, id, ju, j0, jd;
	KValarray<double>		&beta =  _process.get_mean_reversion();
	rho = _process.get_correlation()[0][1];

	KTimeLine::iterator	iter, pIter, nIter;
	/*find the last node, and get the max memory needed */
	iter = tl_end();
	iter--;
	TimeSlice2D<double>	statePr(iter->second[__irTreeMaxLimit1F], iter->second[__irTreeMaxLimit2F]);
	TimeSlice2D<double>	newStatePr(iter->second[__irTreeMaxLimit1F], iter->second[__irTreeMaxLimit2F]);
	/* the first irDrift*/
	iter = tl_begin();
	iter->second[__irTreeDrift] = iter->second[__fwdRateTillNextTp];
	iter++;
	/* update state price for the second node*/
	zeroPrice = iter->second[__zeroPrice];
	statePr( 1,  1) = 1./(2*K_JUMPCOEFF) * 1./(2*K_JUMPCOEFF) * zeroPrice;
	statePr( 1, -1) = statePr(1,  1);
	statePr(-1, -1) = statePr(1,  1);
	statePr(-1,  1) = statePr(1,  1);
	statePr( 1,  0) = 1./(2*K_JUMPCOEFF) * (1. - 1./K_JUMPCOEFF) * zeroPrice;
	statePr( 0, -1) = statePr(1, 0);
	statePr( 0,  1) = statePr(1, 0);
	statePr(-1,  0) = statePr(1, 0);
	statePr( 0,  0) = (1. - 1./K_JUMPCOEFF) * (1. - 1./K_JUMPCOEFF) * zeroPrice;

	while(iter != tl_end())
	{									
		nIter = pIter = iter;
		pIter--;
		nIter++;
		if(nIter == tl_end())
			break;

		timeTillTp = pIter->second[__timeTillNextTp];
		timeTillNextTp = iter->second[__timeTillNextTp];
		du = sqrt (K_JUMPCOEFF * timeTillTp);      
		jump[0] = (pIter->second[__irTreeSpotVol1F] + rho * pIter->second[__irTreeSpotVol2F]) * du;   
    	jump[1] = sqrt(1. - rho * rho) * pIter->second[__irTreeSpotVol2F] * du;   
		expJump[0] = exp (jump[0]);
    	expJump[1] = exp (jump[1]);

		dim1 = iter->second[__irTreeLimit1F];
		dim2 = iter->second[__irTreeLimit2F];
		nextDim1 = nIter->second[__irTreeLimit1F];
		nextDim2 = nIter->second[__irTreeLimit2F];
	
		grid[0]  =  pIter->second[__irTreeDrift] * iter->second[__fwdRateTillNextTp]
					/pIter->second[__fwdRateTillNextTp] * exp (-dim1 * jump[0]);	
 		for(i=0; i<3; i++)
			P[i] = 0.; 
		for (i = -dim1; i <= dim1; i++)       	
		{                          
			grid[1] = grid[0] * exp (-dim2*jump[1]);
			for (j = -dim2; j <= dim2; j++)	
			{
				temp =  statePr(i, j)/(1. + grid[1]);
				P[0] += temp;
				temp *= grid[1] / (1. + grid[1]);
				P[1] += temp;
				temp *= grid[1] / (1. + grid[1]);
				P[2] += temp;
				grid[1] *= expJump[1];
     		}   
			grid[0] *= expJump[0];
		}
		P[0] -= nIter->second[__zeroPrice];  
		iter->second[__irTreeDrift] =  pIter->second[__irTreeDrift] *iter->second[__fwdRateTillNextTp]
										/pIter->second[__fwdRateTillNextTp] 
										*(1. - (-P[1] + sqrt(P[1]*P[1] - 4* P[0]*P[2]))/(2.*P[2]));								

		/*Update State Prices*/
		for (i = -nextDim1; i <= nextDim1; i++)     
			for (j = -nextDim2; j <= nextDim2; j++)
				newStatePr(i, j) = 0.;

		d[0]  = (1. - beta[0] * timeTillNextTp)*pIter->second[__irTreeSpotVol1F] * sqrt(timeTillTp)
				 /(iter->second[__irTreeSpotVol1F]  * sqrt(timeTillNextTp)) - 1.;
		d[1]  = (1. - beta[1] * timeTillNextTp)*pIter->second[__irTreeSpotVol1F] * sqrt(timeTillTp)
				 /(iter->second[__irTreeSpotVol1F]  * sqrt(timeTillNextTp)) - 1.;
		
		grid[0] = iter->second[__irTreeDrift] * exp(-dim1*jump[0]);
		for (i = -dim1; i <= dim1; i++)      
		{   
			Pi = d[0] * i;	
			Qi = rho / sqrt(1. - rho * rho) 
				* (d[1] - d[0]) * i;
			if(Pi > 0)
				l = (int)(Pi + 0.5);
			else		
				l = (int)(Pi - 0.5);
			Pi -= l;
			pu = .5 * (1./K_JUMPCOEFF + Pi + Pi * Pi);
			pd = pu - Pi;
			p0 = 1. - pu - pd;

			iu = Min(nextDim1, Max(-nextDim1, i+1+l));
			i0 = Min(nextDim1, Max(-nextDim1, i+l  ));
			id = Min(nextDim1, Max(-nextDim1, i-1+l));

			grid[1] = grid[0] * exp (-dim2*jump[1]);
			for (j = -dim2; j <= dim2; j++)		
			{
				Qij = Qi + d[1] * j;
				if(Qij > 0)
					m = (int)(Qij + 0.5);
				else		
					m = (int)(Qij - 0.5);
				Qij -= m;
				qu = .5 * (1./K_JUMPCOEFF + Qij + Qij * Qij);
				qd = qu - Qij;
				q0 = 1. - qu - qd;
				ju = Min(nextDim2, Max(-nextDim2, j+1+m));
				j0 = Min(nextDim2, Max(-nextDim2, j+m  ));
				jd = Min(nextDim2, Max(-nextDim2, j-1+m));
            	temp = statePr(i, j)/(1. + grid[1]);
				newStatePr(iu, ju) += pu *qu * temp;  
				newStatePr(iu, j0) += pu *q0 * temp; 	
				newStatePr(iu, jd) += pu *qd * temp;
				newStatePr(i0, ju) += p0 *qu * temp;
				newStatePr(i0, j0) += p0 *q0 * temp;
				newStatePr(i0, jd) += p0 *qd * temp;
				newStatePr(id, ju) += pd *qu * temp;
				newStatePr(id, j0) += pd *q0 * temp;
				newStatePr(id, jd) += pd *qd * temp;

				grid[1] *= expJump[1];
			}  
			grid[0] *= expJump[0];
		}
		p = 0.;	
		for (i = -nextDim1; i <= nextDim1; i++)     
			for (j = -nextDim2; j <= nextDim2; j++)
			{                                                               
				statePr(i, j) = newStatePr(i,j);
		  		p += statePr(i,j);
	    	}
		p /= nIter->second[__zeroPrice];
		if(p < 1. - 0.00001)
			throw KException("IRTree2F::irDrift Failed!");
		//next iter
		iter++;
	}

}

