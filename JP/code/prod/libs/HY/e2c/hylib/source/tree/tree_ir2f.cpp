/**********************************************************
*
*	tree.cpp
*
***********************************************************/
#include "kplatform.h"  
#include "tree_ir2f.h"
#include <iomanip>
#include <math.h>

#include "kcurve.h"
#include "kplatdep.h"


KValarray<int>	IRTree2F::get_curr_limit(void)const
{
	KValarray<int>	lim(2);
	KTimePoint &tp = get_curr_tp();
	lim[0] = tp[__irTreeLimit1F]; 
	lim[1] = tp[__irTreeLimit2F]; 
	return  lim;
}
KValarray<int>	IRTree2F::get_next_limit(void)const
{
	KValarray<int>	lim(2);
	KTimePoint &tp = get_next_tp();
	lim[0] = tp[__irTreeLimit1F]; 
	lim[1] = tp[__irTreeLimit2F]; 
	return  lim;
}

KValarray<int>	IRTree2F::get_max_limit(void)const
{
	KValarray<int>	lim(2);
	KTimePoint &tp = get_curr_tp();
	lim[0] = tp[__irTreeMaxLimit1F];
	lim[1] = tp[__irTreeMaxLimit2F];
	return  lim;
}


void	IRTree2F::treeSpotVol() 
{
	KTimeLine::iterator iter, nIter;

	for(iter = tl_begin();iter != tl_end(); iter++)
	{
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;

		KValarray<double> &spotVol = _process.get_avg_vol(iter->first, nIter->first);
		iter->second[__irTreeSpotVol1F] = spotVol[0];
		iter->second[__irTreeSpotVol2F] = spotVol[1];
	} 

	/* Adjustment to the spot volatility as FwdRate is simple rate */
	double  fwd;
	KValarray<double>	&beta = _process.get_mean_reversion();
	for(iter = tl_begin();iter != tl_end(); iter++)
	{	
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;
		fwd = nIter->second[__fwdRateTillNextTp];
		iter->second[__irTreeSpotVol1F]  *= (1+fwd)*log(1+fwd)/fwd * A (beta[0], iter->second[__timeTillNextTp]);
		iter->second[__irTreeSpotVol2F]  *= (1+fwd)*log(1+fwd)/fwd * A (beta[1], iter->second[__timeTillNextTp]);

	}
	//the last
	nIter = tl_end();
	nIter--;
	iter= nIter;
	iter--;
	nIter->second[__irTreeSpotVol1F] = iter->second[__irTreeSpotVol1F];
	nIter->second[__irTreeSpotVol2F] = iter->second[__irTreeSpotVol2F];
} 


void	IRTree2F::limits (void)
{
	double	currentMax;						
	double	vol[2], u, d[2], Q;
	KValarray<double>	&beta = _process.get_mean_reversion();
	double	corr = _process.get_correlation()[0][1];
	
	KTimeLine::iterator	iter, nIter;
	iter = tl_begin();
	/* limits for the first node*/
	iter->second[__irTreeLimit1F] = 0;
	iter->second[__irTreeLimit2F] = 0;
	iter->second[__irTreeMaxLimit1F] = 0;
	iter->second[__irTreeMaxLimit2F] = 0;

	vol[0] = vol[1] = 0.;
	while(iter != tl_end())
	{   
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;

		u = iter->second[__irTreeSpotVol1F] * iter->second[__irTreeSpotVol1F] * iter->second[__timeTillNextTp];
		vol[0] *= exp (- 2. * beta[0] * iter->second[__timeTillNextTp]);  
		vol[0] += u	* A (2. * beta[0], iter->second[__timeTillNextTp]);
		u *= K_JUMPCOEFF;
		currentMax  = floor(1 + K_NUM_SIGMA * sqrt(vol[0] / u));
		d[0] = (1. - beta[0]  * nIter->second[__timeTillNextTp]);  
		d[0] *= iter->second[__irTreeSpotVol1F] * sqrt(iter->second[__timeTillNextTp]);
		d[0] /= nIter->second[__irTreeSpotVol1F] * sqrt(nIter->second[__timeTillNextTp]);
		d[0] -= 1.;
		nIter->second[__irTreeLimit1F] = Min(currentMax, iter->second[__irTreeLimit1F] + 1 
										+ int(fabs(d[0]) * iter->second[__irTreeLimit1F] + .5));

		u = iter->second[__irTreeSpotVol2F] * iter->second[__irTreeSpotVol2F] * iter->second[__timeTillNextTp];
		vol[1] *= exp (- 2. * beta[1] * iter->second[__timeTillNextTp]);  
		vol[1] += u	* A (2. * beta[1], iter->second[__timeTillNextTp]);
		u *= K_JUMPCOEFF;
		currentMax  = floor(1 + K_NUM_SIGMA * sqrt(vol[1] / u));
		d[1] = (1. - beta[1]  * nIter->second[__timeTillNextTp]);  
		d[1] *= iter->second[__irTreeSpotVol2F] * sqrt(iter->second[__timeTillNextTp]);
		d[1] /= nIter->second[__irTreeSpotVol2F] * sqrt(nIter->second[__timeTillNextTp]);
		d[1] -= 1.;
		Q = fabs(corr/sqrt(1. - corr * corr)*(d[1] - d[0]))*iter->second[__irTreeLimit1F];
		nIter->second[__irTreeLimit2F] = Min(currentMax, iter->second[__irTreeLimit2F] + 1 
										+ int(Q + fabs(d[1]) * iter->second[__irTreeLimit2F] + .5));
		nIter->second[__irTreeMaxLimit1F] = Max(iter->second[__irTreeMaxLimit1F] , nIter->second[__irTreeLimit1F]);
		nIter->second[__irTreeMaxLimit2F] = Max(iter->second[__irTreeMaxLimit2F] , nIter->second[__irTreeLimit2F]);

		iter++;
	}
}


/*****  Lattice   ***********************************************************/
/*
*       Position the nodes in the lattice: calculate the one period discount
*       factor and the probabilities at each node in the tree.
*/
void	IRTree2F::buildLattice (void)
{
	double	timeTillTp,
			timeTillNextTp,
			Pi, Qi, Qij,
			pu, pd, p0, qu, qd, q0,
			jump[2], expJump[2],   
			du, grid[2];
	int		i, j,l, m, dim1, dim2, 
			nextDim1, nextDim2;		
	double	d[2];
	KValarray<double>	&beta = _process.get_mean_reversion();
	double	rho = _process.get_correlation()[0][1];

	KVector(double)	prob(9);
	KValarray<int>	shift(2);
	
	KTimePoint	&tp = get_curr_tp();
	/* if it is the first time point*/
	if (tp[__timeFromValueDate] == 0.)
	{
		IR2FTreeParam &param = _irParam(0, 0);
		param._idx1[0] = -1;
		param._idx1[1] =  0;
		param._idx1[2] =  1;
		param._idx2[0] = -1;
		param._idx2[1] =  0;
		param._idx2[2] =  1;
		param._prob[2][2] = 1/(4*K_JUMPCOEFF*K_JUMPCOEFF);
		param._prob[2][1] = (1. - 1/K_JUMPCOEFF)/(2*K_JUMPCOEFF);
		param._prob[2][0] = 1/(4*K_JUMPCOEFF*K_JUMPCOEFF);
		param._prob[1][2] = (1. - 1/K_JUMPCOEFF)/(2*K_JUMPCOEFF);
		param._prob[1][1] = (1. - 1/K_JUMPCOEFF) * (1. - 1/K_JUMPCOEFF);
		param._prob[1][0] = (1. - 1/K_JUMPCOEFF)/(2*K_JUMPCOEFF);
		param._prob[0][2] = 1/(4*K_JUMPCOEFF*K_JUMPCOEFF);
		param._prob[0][1] = (1. - 1/K_JUMPCOEFF)/(2*K_JUMPCOEFF);
		param._prob[0][0] = 1/(4*K_JUMPCOEFF*K_JUMPCOEFF);
		param._idxDisc = 1. / (1. + tp[__fwdRateTillNextTp]);
	}
	else
	{
		KTimePoint	&pTp = get_prev_tp();
		KTimePoint	&nTp = get_next_tp();

		timeTillTp = pTp[__timeTillNextTp];
		timeTillNextTp =  tp[__timeTillNextTp];
		du = sqrt (K_JUMPCOEFF * timeTillTp);      
		
		jump[0] = (pTp[__irTreeSpotVol1F] + rho * pTp[__irTreeSpotVol2F]) * du;   
		jump[1] = sqrt(1. - rho*rho) * pTp[__irTreeSpotVol2F] * du;   
		expJump[0] = exp (jump[0]);
		expJump[1] = exp (jump[1]);

		dim1 = tp[__irTreeLimit1F];
		dim2 = tp[__irTreeLimit2F];
		nextDim1 = nTp[__irTreeLimit1F];
		nextDim2 = nTp[__irTreeLimit2F];

		d[0] = (1. - beta[0] * timeTillNextTp)  
				*pTp[__irTreeSpotVol1F] * sqrt(timeTillTp)
				 /(tp[__irTreeSpotVol1F] * sqrt(timeTillNextTp)) - 1.;
		d[1] = (1. - beta[1] * timeTillNextTp)  
				*pTp[__irTreeSpotVol2F] * sqrt(timeTillTp)
				 /(tp[__irTreeSpotVol2F] * sqrt(timeTillNextTp)) - 1.;

		grid[0]  =	tp[__irTreeDrift] * exp (-dim1*jump[0]);
		for (i = -dim1; i <= dim1; i++)                   	
		{   
			Pi = d[0] * i;	
			Qi = rho/sqrt(1. - rho * rho) * (d[1] - d[0]) * i;
			if(Pi > 0)
				l = (int)(Pi + 0.5);
			else		
				l = (int)(Pi - 0.5);

			Pi -= l;
			pu = .5 * (1./K_JUMPCOEFF + Pi + Pi * Pi);
			pd = pu - Pi;
			p0 = 1. - pu - pd;

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

				IR2FTreeParam &param = _irParam(i, j);
				param._idx1[0] = Min(nextDim1, Max(-nextDim1, i+l-1));
				param._idx1[1] = Min(nextDim1, Max(-nextDim1, i+l  ));
				param._idx1[2] = Min(nextDim1, Max(-nextDim1, i+l+1));
				param._idx2[0] = Min(nextDim2, Max(-nextDim2, j+m-1));
				param._idx2[1] = Min(nextDim2, Max(-nextDim2, j+m  ));
				param._idx2[2] = Min(nextDim2, Max(-nextDim2, j+m+1));
				param._prob[2][2] = pu * qu;
				param._prob[2][1] = pu * q0;
				param._prob[2][0] = pu * qd;
				param._prob[1][2] = p0 * qu;
				param._prob[1][1] = p0 * q0;
				param._prob[1][0] = p0 * qd;
				param._prob[0][2] = pd * qu;
				param._prob[0][1] = pd * q0;
				param._prob[0][0] = pd * qd;
				param._idxDisc = 1. / (1. + grid[1]);
				grid[1] *= expJump[1];
			}
			grid[0] *= expJump[0];
		}  
	}
}




