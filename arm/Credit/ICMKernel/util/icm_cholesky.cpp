#include "math.h"
#include <vector>
#include "ICMKernel\util\icm_cholesky.h"
#include "ICMKernel\util\icm_qmatrix.h"

void cholesky::init(ICM_QMatrix<double>& matrice,int size)
{
	int i=0,j=0;
	m_matrix= ICM_QMatrix<double>(size,size);

	for(i=0;i<size;i++)
		{
			for(j=0;j<size;j++)
				m_matrix.SetValue(i,j, matrice.Getvalue(i,j));
		}
}



int
cholesky::choldc(int n)
{
	int     i = 0, j = 0, k = 0;
	double  sum = 0.0,tmp=1.0;
	m_eigenvalues.resize(n);

     
    {	

	for (i = 0; i < n; i++)
	{
		for (j = i; j < n; j++)
		{
			for (sum = m_matrix.Getvalue(i,j), k = i - 1; k >= 0; k--)
				sum -= m_matrix.Getvalue(i,k) * m_matrix.Getvalue(j,k);
			if (i == j)
			{
				if (sum <= 0.0)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR : The Correlation Matrix is not define >0 !");

				m_eigenvalues[i] = sqrt(sum);
				tmp*=m_eigenvalues[i];
			}
			else
            {
				m_matrix.SetValue(j,i, sum / m_eigenvalues[i]);
            }
		}
	}
    }	
    

	return 1;
}

void
cholesky::getindex()

{
	int n=m_eigenvalues.size();
	m_index.resize(n);

	//sort_descending(m_eigenvalues, m_index);
}