#include "gpbase\gptrigomatrix.h"
#include "gpbase\numericconstant.h"
#include <math.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////
///	Utility function
///	Routine: TrigoMatrix
///	Returns: ARM_GP_Matrix
///	Action : Compute matrix with trigonometric 
/// coefficient cij = cos(alpha*PI*(i-j)/N)
////////////////////////////////////////////////////
ARM_GP_Matrix* TrigoMatrix(size_t n, double alpha)
{
	ARM_GP_Matrix* matrix = new ARM_GP_Matrix(n,n);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (n > 1)
			{
				(*matrix)(i,j) = cos(ARM_NumericConstants::ARM_PI*alpha*(i-j)/(n-1));
			}
			else
			{
				(*matrix)(i,j) = 0.0;
			}
		}

	return matrix;
}

CC_END_NAMESPACE()
