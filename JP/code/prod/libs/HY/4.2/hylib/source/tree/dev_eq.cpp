#include "tree_eq.h"



/*****  Dev  ***************************************************************/
/*
*       Discounted expected value. 
*/
void	EqTree::dev(DTimeSlice &ts)	
{
	int		i, iu, i0, id;
	double	disc;
	KTimePoint	&tp = get_curr_tp();

	int	dim  = get_curr_limit()[0];
	disc = 1./(1 + tp[__fwdRateTillNextTp]);

	for (i = -dim; i <= dim; i++)
	{
		EqTreeParam	&param = _eqParam[i];
		iu = param._idx[2];
		i0 = param._idx[1];
		id = param._idx[0];
		_temp[i]  =  ( param._prob[2] * ts[iu]+ 
						param._prob[1] * ts[i0]+ 
						param._prob[0] * ts[id] ) * disc;
	}
	for(i=-dim; i<= dim; i++)
		ts[i] = _temp[i];
	
}

void	EqTree::dev(FTimeSlice &ts)	
{
	int		i, iu, i0, id;
	double	disc;
	KTimePoint	&tp = get_curr_tp();

	int	dim  = get_curr_limit()[0];
	disc = 1./(1 + tp[__fwdRateTillNextTp]);

	for (i = -dim; i <= dim; i++)
	{
		EqTreeParam	&param = _eqParam[i];
		iu = param._idx[2];
		i0 = param._idx[1];
		id = param._idx[0];
		_tempFunc[i]  =  ( param._prob[2] * ts[iu]+ 
						param._prob[1] * ts[i0]+ 
						param._prob[0] * ts[id] ) * disc;
	}
	for(i=-dim; i<= dim; i++)
		ts[i] = _tempFunc[i];
	
}
