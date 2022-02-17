#include "tree_ir2f.h"

void	IRTree2F::dev(DTimeSlice &ts)
{
	int		i, j;
	int		ku, k0, kd;
	int		lu, l0, ld;

	int	dim1  = get_curr_limit()[0];
	int	dim2  = get_curr_limit()[1];
	for (i = -dim1; i <= dim1; i++)			
    	for (j = -dim2; j <= dim2; j++)		
		{  
			IR2FTreeParam	&param = _irParam(i, j);
		
			ku = param._idx1[2];
			k0 = param._idx1[1];
			kd = param._idx1[0];
			lu = param._idx2[2];
			l0 = param._idx2[1];
			ld = param._idx2[0];
			_temp(i, j) = param._idxDisc *
							(param._prob[2][2] * ts(ku, lu) +
							 param._prob[2][1] * ts(ku, l0) +
							 param._prob[2][0] * ts(ku, ld) +
							 param._prob[1][2] * ts(k0, lu) +
							 param._prob[1][1] * ts(k0, l0) +
							 param._prob[1][0] * ts(k0, ld) +
							 param._prob[0][2] * ts(kd, lu) +
							 param._prob[0][1] * ts(kd, l0) +
							 param._prob[0][0] * ts(kd, ld));
		}
			
	for (i = -dim1; i <= dim1; i++)
		for (j = -dim2; j <= dim2; j++)
			ts(i, j) = _temp(i, j);
}

void	IRTree2F::dev(FTimeSlice &ts)
{
	int		i, j;
	int		ku, k0, kd;
	int		lu, l0, ld;

	int	dim1  = get_curr_limit()[0];
	int	dim2  = get_curr_limit()[1];
	for (i = -dim1; i <= dim1; i++)			
    	for (j = -dim2; j <= dim2; j++)		
		{  
			IR2FTreeParam	&param = _irParam(i, j);
		
			ku = param._idx1[2];
			k0 = param._idx1[1];
			kd = param._idx1[0];
			lu = param._idx2[2];
			l0 = param._idx2[1];
			ld = param._idx2[0];

			_tempFunc(i, j) = param._idxDisc*
							(param._prob[2][2] * ts(ku, lu) +
							 param._prob[2][1] * ts(ku, l0) +
							 param._prob[2][0] * ts(ku, ld) +
							 param._prob[1][2] * ts(k0, lu) +
							 param._prob[1][1] * ts(k0, l0) +
							 param._prob[1][0] * ts(k0, ld) +
							 param._prob[0][2] * ts(kd, lu) +
							 param._prob[0][1] * ts(kd, l0) +
							 param._prob[0][0] * ts(kd, ld));
		}
	for (i = -dim1; i <= dim1; i++)
		for (j = -dim2; j <= dim2; j++)
			ts(i, j) = _tempFunc(i, j);
}


