#ifndef _CHOLESKY_TOOLS
#define _CHOLESKY_TOOLS

#include <vector>
//#include "basicTemplateFunctions.h"
#include "icm_qmatrix.h"

class cholesky
{
public:
	cholesky() {}

	virtual ~cholesky() {};

	void init (ICM_QMatrix<double>& matrice,int size);
	int	choldc(int n);
	void getindex();

	ICM_QMatrix<double> m_matrix;
	std::vector<double> m_eigenvalues;
	std::vector<int> m_index;  // stocke les indices des valeurs propres rangées dans l'ordre décroissant

};

#endif