/**********************************************************
*
*	treeIR1F.cpp
*
***********************************************************/
#include <iomanip>
#include <math.h>
#include "tree_ir1f.h"
#include "kcurve.h"

KValarray<int>	IRTree1F::get_curr_limit(void)const
{
	KValarray<int>	lim(1);
	KTimePoint &tp = get_curr_tp();
	lim[0] = tp[__irTreeLimit1F];
	return  lim;
}
KValarray<int>	IRTree1F::get_next_limit(void)const
{
	KValarray<int>	lim(1);
	KTimePoint &tp = get_next_tp();
	lim[0] = tp[__irTreeLimit1F];
	return  lim;
}
KValarray<int>	IRTree1F::get_max_limit(void)const
{
	KValarray<int>	lim(1);
	KTimePoint &tp = get_curr_tp();
	lim[0] = tp[__irTreeMaxLimit1F];
	return  lim;
}


void	IRTree1F::treeSpotVol(void) 
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
	} 

	double	fwd;
	double	beta = _process.get_mean_reversion()[0];
	/* Adjustment to the spot volatility as FwdRate is simple rate */
	for(iter = tl_begin();iter != tl_end(); iter++)
	{	
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;
		fwd = nIter->second[__fwdRateTillNextTp];
		iter->second[__irTreeSpotVol1F]  *= (1.+fwd) * log(1+fwd)/fwd
											* A (beta, iter->second[__timeTillNextTp]);
	}  
	//the last
	nIter = tl_end();
	nIter--;
	iter= nIter;
	iter--;
	nIter->second[__irTreeSpotVol1F] = iter->second[__irTreeSpotVol1F];

} 



void	IRTree1F::limits (void)
{
	double	currentMax;
	double	vol, u, d;
	double	beta = _process.get_mean_reversion()[0];

	KTimeLine::iterator	iter, nIter;
	iter = tl_begin();
	/* limits for the first node*/
	iter->second[__irTreeLimit1F] = 0;
	iter->second[__irTreeMaxLimit1F] = 0;
	vol = 0.;
	while(iter != tl_end())
	{   
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;

		u = iter->second[__irTreeSpotVol1F] * iter->second[__irTreeSpotVol1F] * iter->second[__timeTillNextTp];                  
		vol *= exp (- 2. * beta * iter->second[__timeTillNextTp]);  
		vol += u * A (2. * beta, iter->second[__timeTillNextTp]);
		u *= K_JUMPCOEFF; 
		currentMax = floor(1 + K_NUM_SIGMA * sqrt(vol / u));
		
		d = (1. - beta  * nIter->second[__timeTillNextTp]);  
		d *= iter->second[__irTreeSpotVol1F] * sqrt(iter->second[__timeTillNextTp]);
		d /= nIter->second[__irTreeSpotVol1F] * sqrt(nIter->second[__timeTillNextTp]);
		d -= 1.;
		nIter->second[__irTreeLimit1F] = Min(currentMax, iter->second[__irTreeLimit1F] + 1 
											+ int(fabs(d) * iter->second[__irTreeLimit1F] + .5));
		nIter->second[__irTreeMaxLimit1F] = Max(iter->second[__irTreeMaxLimit1F], nIter->second[__irTreeLimit1F]);
		iter++;
	}

}



/*****  Lattice   ***********************************************************/
/*
*       Position the nodes in the lattice: calculate the one period discount
*       factor and the probabilities at each node in the tree.
*/
void	IRTree1F::buildLattice (void)
{
	double	timeTillTp,
			timeTillNextTp,
			Qi,jump, expJump,   
			du, d, grid;
	int		dim,nextDim,
			i, m;
	double	beta = _process.get_mean_reversion()[0];

	KTimePoint	&tp = get_curr_tp();
	/* if it is the first time point*/
	if (tp[__timeFromValueDate] == 0.)
	{
		IR1FTreeParam	&param = _irParam[0];
		param._idx[0] = -1;
		param._idx[1] =  0;
		param._idx[2] =  1;
		param._prob[0] = 1./(2*K_JUMPCOEFF);
		param._prob[1] = 1. - 1./K_JUMPCOEFF;
		param._prob[2] = 1./(2*K_JUMPCOEFF);
		param._idxDisc = 1. / (1. + tp[__fwdRateTillNextTp]);
	}
	else
	{
		KTimePoint	&pTp = get_prev_tp();
		KTimePoint	&nTp = get_next_tp();

		timeTillTp = pTp[__timeTillNextTp];
		timeTillNextTp =  tp[__timeTillNextTp];
		du = sqrt (K_JUMPCOEFF * timeTillTp);      
		jump  = pTp[__irTreeSpotVol1F] * du;   
		expJump  = exp (jump );

		d = (1. - beta  * timeTillNextTp)  
			*pTp[__irTreeSpotVol1F] * sqrt(timeTillTp)
				 /(tp[__irTreeSpotVol1F] * sqrt(timeTillNextTp)) - 1.;
		
		dim = tp[__irTreeLimit1F];
		nextDim = nTp[__irTreeLimit1F];
	
		grid = tp[__irTreeDrift] * exp (-dim*jump);
		for (i = -dim; i <= dim; i++)			
		{		
			Qi = d * i;
			if(Qi > 0)
				m = (int)(Qi + 0.5);
			else		
				m = (int)(Qi - 0.5);
		
			IR1FTreeParam	&param = _irParam[i];	

			param._idx[0]  = Min(nextDim, Max(-nextDim, i+m-1));;                                                  
			param._idx[1]  = Min(nextDim, Max(-nextDim, i+m  ));;                                                  
			param._idx[2]  = Min(nextDim, Max(-nextDim, i+m+1));;                                                  
			
			Qi -= m;
			param._prob[2] = .5 * (1./K_JUMPCOEFF + Qi + Qi * Qi);
			param._prob[0] = param._prob[2] - Qi;
			param._prob[1] = 1. - param._prob[0] - param._prob[2];
			param._idxDisc = 1. / (1. + grid);
			grid *= expJump ;
		}
	}
}




