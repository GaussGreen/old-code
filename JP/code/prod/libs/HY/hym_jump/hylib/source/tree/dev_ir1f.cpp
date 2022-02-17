
#include "tree_ir1f.h"


void	IRTree1F::dev(DTimeSlice &ts)
{
	int		i, iu, i0, id;
	int	dim  = get_curr_limit()[0];
	for (i = -dim; i <= dim; i++)			
	{
		IR1FTreeParam	&param = _irParam[i];
		iu = param._idx[2];
		i0 = param._idx[1];
		id = param._idx[0];
		_temp[i] = param._idxDisc * 
					(	param._prob[2] * ts[iu] + 
						param._prob[1] * ts[i0] + 
						param._prob[0] * ts[id]  );

	} 
	for (i = -dim; i <= dim; i++)			
		ts[i] = _temp[i];
}



void	IRTree1F::dev(FTimeSlice &ts)
{
	int	i, iu, i0, id;
	int	dim  = get_curr_limit()[0];
	for (i = -dim; i <= dim; i++)			
	{
		IR1FTreeParam	&param = _irParam[i];
		iu = param._idx[2];
		i0 = param._idx[1];
		id = param._idx[0];
		_tempFunc[i] =	param._idxDisc *
						(param._prob[2] * ts[iu] + 
						param._prob[1] * ts[i0] + 
						param._prob[0] * ts[id]  );

	} 
	for (i = -dim; i <= dim; i++)			
		ts[i] = _tempFunc[i];
}



