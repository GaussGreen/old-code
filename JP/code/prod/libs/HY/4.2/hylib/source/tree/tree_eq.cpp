//TreeEq.cpp
#include "tree_eq.h"
#include <math.h>

KValarray<int>	EqTree::get_curr_limit(void)const
{
	KValarray<int> lim(1);
	KTimePoint	&tp = get_curr_tp();
	lim[0] = tp[__eqTreeLimit];
	return lim;
}
KValarray<int>	EqTree::get_next_limit(void)const
{
	KValarray<int> lim(1);
	KTimePoint	&tp = get_next_tp();
	lim[0] = tp[__eqTreeLimit];
	return lim;
}

KValarray<int>	EqTree::get_max_limit(void)const
{
	KValarray<int> lim(1);
	KTimePoint	&tp = get_curr_tp();
	lim[0] = tp[__eqTreeMaxLimit];
	return lim;
}



void	EqTree::set_process(const KEqProcess &p) 
{
	_process = p;
	KTimeLine	&tl = p.get_timeline();
	//insert time line
	for(KTimeLine::iterator iter = tl.begin(); iter!= tl.end(); iter++)
		insert_tp(iter->first, iter->second);
}


void	EqTree::drift(void)
{
	double	variance = 0;
	KRateCurve	&zr = get_index_zero();

	for(KTimeLine::iterator iter = tl_begin(); iter!=  tl_end(); iter++)
	{
		iter->second[__eqFwdPrice] = _process.get_fwd_price(iter->first, zr);
		iter->second[__eqTreeMidNode] = iter->second[__eqFwdPrice]*exp(-.5 *variance);
		variance += iter->second[__eqTreeSpotVol] * iter->second[__eqTreeSpotVol] 
				*iter->second[__timeTillNextTp];
	}  
}

// Spot Vol interpolation
void	EqTree::treeSpotVol(void) 
{
	KTimeLine::iterator iter, nIter;
	for(iter = tl_begin(); iter != tl_end(); iter++)					
	{
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;
		KValarray<double>	&spotVol = _process.get_avg_vol(iter->first, nIter->first);
		iter->second[__eqTreeSpotVol] = spotVol[0];
	}
	/*the last spot vol is irrelevant*/
	nIter = tl_end();
	nIter--;
	iter = nIter;
	iter--;
	nIter->second[__eqTreeSpotVol] = iter->second[__eqTreeSpotVol];
}

void	EqTree::limits (void)
{
	double	currentMax;
	double	vol,u, d;

	KTimeLine::iterator	iter, nIter;

	iter = tl_begin();
	/* limits for the first node*/
	iter->second[__eqTreeLimit] = 0;
	iter->second[__eqTreeMaxLimit] = 0;
	vol = 0;
	while(iter != tl_end())
	{   
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;
		u = iter->second[__eqTreeSpotVol] * iter->second[__eqTreeSpotVol]  * iter->second[__timeTillNextTp];
		vol += u;
		u *= K_JUMPCOEFF;
		currentMax = floor(1 + K_NUM_SIGMA * sqrt(vol/ u));
		d = iter->second[__eqTreeSpotVol] * sqrt(iter->second[__timeTillNextTp]);
		d /= nIter->second[__eqTreeSpotVol] * sqrt(nIter->second[__timeTillNextTp]);
		d -= 1.;
		nIter->second[__eqTreeLimit] = Min(currentMax, iter->second[__eqTreeLimit] + 1 
											+ int(fabs(d) * iter->second[__eqTreeLimit] + .5));
		nIter->second[__eqTreeMaxLimit] = Max(iter->second[__eqTreeMaxLimit], nIter->second[__eqTreeLimit]);
		iter++;
	}
}



DTimeSlice	EqTree::get_stock_price(void)const 
{
	DTimeSlice	price(this);
	int	dim = get_curr_limit()[0];
	int	i;

	for(i=-dim; i<= dim; i++)
		price[i] = _eqParam[i]._under;

	return price;

}



void	EqTree::buildLattice (void)
{
	double	timeTillTp,timeTillNextTp,
			d,jump, expJump, grid;
	double	Pi;
	int		dim, nextDim,i, l;
	KTimePoint	&tp = get_curr_tp();
	/* if it is the first time point*/
	if (tp[__timeFromValueDate] == 0.)
	{
		EqTreeParam	&param = _eqParam[0];
		param._idx[0] = -1;
		param._idx[1] =  0;
		param._idx[2] =  1;
		param._prob[0] = 1/(2*K_JUMPCOEFF);
		param._prob[1] = 1. - 1/K_JUMPCOEFF;
		param._prob[2] = 1/(2*K_JUMPCOEFF);
		param._under = tp[__eqTreeMidNode];
	}
	else
	{
		KTimePoint	&pTp = get_prev_tp();
		KTimePoint	&nTp = get_next_tp();
		timeTillTp = pTp[__timeTillNextTp];
		timeTillNextTp =  tp[__timeTillNextTp];

		jump     = pTp[__eqTreeSpotVol] * sqrt (K_JUMPCOEFF * timeTillTp);   
		expJump = exp (jump);

		dim = tp[__eqTreeLimit];
		nextDim = nTp[__eqTreeLimit];
		
		d  = pTp[__eqTreeSpotVol] * sqrt(timeTillTp)
				/(tp[__eqTreeSpotVol] * sqrt(timeTillNextTp)) - 1.;
		grid = tp[__eqTreeMidNode] * exp (-dim*jump);
		for (i = -dim; i <= dim; i++)                   	
		{   
			Pi =  d * i;
			if(Pi > 0)
				l = (int)(Pi + 0.5);
			else		
				l = (int)(Pi - 0.5);
			
			Pi -= l;
			
			EqTreeParam	&param = _eqParam[i];
			param._idx[0] = Min(nextDim, Max(-nextDim, i+l-1));
			param._idx[1] = Min(nextDim, Max(-nextDim, i+l  ));
			param._idx[2] = Min(nextDim, Max(-nextDim, i+l+1));
	
			param._prob[2] = .5 * (1./K_JUMPCOEFF + Pi + Pi * Pi);
			param._prob[0] = param._prob[2] - Pi;
			param._prob[1] = 1 - param._prob[0]-param._prob[2];
			param._under = grid;
			
			grid *= expJump;
		}  
	}

}