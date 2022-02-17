/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Integrator_Dim_2_DIM_2.CPP
	PROJECT:	UTIL
	
	DESCRIPTION:	Numerical Integrator

  -----------------------------------------------------------------

 	CREATION:	February 16, 2006

	LAST MODIF:	February 16, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */


# include "ICM_Integrator_Dim_2.h"

// # include <stdlib.h>
# include <memory.h>
# include <vector>
# include <map>

# include "ICMKernel\glob\icm_constants.h"
# include "ARMKernel\glob\expt.h"
# include "ICMKernel\util\icm_qmatrix.h"

# include "ICMKernel\glob\icm_maths.h"


using namespace std;


void ICM_Integrator_Dim_2::Init()
{
	itsIntegrationStep_2	= 20;
	itsLowBound_2		=	0.0;
	itsUpBound_2		=	0.0;

	itsIntegratorType_2 = qGAUSS_HERMITE;

	itsX_2	=	NULL;
	itsW_2	=	NULL;
}


void ICM_Integrator_Dim_2::BitwiseCopy(const ICM_Integrator_Dim_2* src)
{
	itsIntegrationStep_2	=	src->itsIntegrationStep_2;
	itsLowBound_2			=	src->itsLowBound_2;
	itsUpBound_2			=	src->itsUpBound_2;

	itsIntegratorType_2		=	src->itsIntegratorType_2;

	if (itsX_2)
		delete[] itsX_2;
	if (src->itsX_2)
	{
		itsX_2 = (double *)malloc(itsIntegrationStep_2 * sizeof(double));
		if (itsX_2 != NULL)
			memcpy(itsX_2, src->itsX_2, itsIntegrationStep_2 * sizeof(double));
	}

	if (itsW_2)
		delete[] itsW_2;
	if (src->itsW_2)
	{
		itsW_2 = (double *)malloc(itsIntegrationStep_2 * sizeof(double));
		if (itsW_2 != NULL)
			memcpy(itsW_2, src->itsW_2, itsIntegrationStep_2 * sizeof(double));
	}
}


void ICM_Integrator_Dim_2::Copy(const ICM_Integrator_Dim_2* src)
{
    ICM_Integrator::Copy(src);
	
	if (src != this)
		BitwiseCopy(src);
}


ICM_Integrator_Dim_2* ICM_Integrator_Dim_2 :: Clone() const
{
    ICM_Integrator_Dim_2* theClone = new ICM_Integrator_Dim_2();

    theClone->Copy(this);

    return (theClone);
}

void ICM_Integrator_Dim_2::Reset()
{
	if (itsX_2)
		delete[] itsX_2;
	if (itsW_2)
		delete[] itsW_2;
}


void	ICM_Integrator_Dim_2 :: GetAbscissaVector_2(DoubleVector& data)
{
	data.clear();
	data.resize(itsIntegrationStep_2);
	for (int i=0; i<itsIntegrationStep_2;i++)
		data[i]	=	itsX_2[i];
}

void	ICM_Integrator_Dim_2 :: GetWeightVector_2(DoubleVector& data)
{
	data.clear();
	data.resize(itsIntegrationStep_2);
	for (int i=0; i<itsIntegrationStep_2;i++)
		data[i]	=	itsW_2[i];
}



//----------------------------------------------------------------------------


void  ICM_Integrator_Dim_2::Integrate(	const double& x1,
									const double& x2,
									const double& y1,
									const double& y2,
									void (*fct)(void*, double, double,  double&),
									AddressVector* param,
									double& res)
{
	int	i, j;
	double fxy;
	std::map<int,ICM_QMatrix<double> >::const_iterator it;
	ICM_QMatrix<double>* matrix=NULL;
	std::map<int,ICM_QMatrix<double> >::const_iterator it_2;
	ICM_QMatrix<double>* matrix_2=NULL;
	
	// qTrapeze Data
	int NbStepTrapeze	= itsIntegrationStep - 1;
	int NbStepTrapeze_2 = itsIntegrationStep_2 - 1;

	if (fct == NULL)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
			"ERROR: Null function for Integrator!");
	}

	res	=	0.0;
	double	xi;
	double	size	=	x2 - x1;
	double	yj;
	double	size_2	=	y2 - y1;
	
	double	sum_2;

	int	n	=	itsIntegrationStep;

	//On load les valeurs dans la map en fonction du step d'intégration
	LoadValues(itsIntegrationStep);
	if (itsIntegrationStep_2 != itsIntegrationStep)
		LoadValues(itsIntegrationStep_2);

	//On verifie que la méthode numérique est compatible avec le nombre de pas employé
	CheckNumericalMethod(itsIntegrationStep, itsIntegratorType);
	CheckNumericalMethod(itsIntegrationStep_2, itsIntegratorType_2);

	if (itsIntegratorType != itsIntegratorType_2)
		ICMTHROW(ERR_INVALID_DATA,"Integrator Dim 2 deals only with ONE Numerical Integrator Type!");

	int* nointegrationcoef_x = (int*)param->GetLastMinusOne();
	int* nointegrationcoef_y = (int*)param->GetLast();

	switch (itsIntegratorType)
	{
	case qGAUSS_LEGENDRE:
		{
		it = itsMap.find(itsIntegrationStep); 
		if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep];

		if (itsIntegrationStep == itsIntegrationStep_2)
			matrix_2	=	matrix;
		else
		{
			it_2 = itsMap.find(itsIntegrationStep_2); 
			if (it_2!=itsMap.end()) matrix = &itsMap[itsIntegrationStep_2];
		}

		if (matrix == NULL)
		{
			ComputeAbscissasAndWeightsGaussLegendre(0.0, 1.0, itsIntegrationStep, itsX, itsW);

			ComputeAbscissasAndWeightsGaussLegendre(0.0, 1.0, itsIntegrationStep_2, itsX_2, itsW_2);
		}
		
		for (i=0;i<itsIntegrationStep;i++)
		{
			*nointegrationcoef_x = i;
			xi	=	x1 + size * (*matrix)(i,0);
			
			sum_2	=	0.0;

			for (j=0; j <itsIntegrationStep_2; j++)
			{
				*nointegrationcoef_y	=	j;
				yj	=	y1 + size_2 * (*matrix_2)(j,0);

				(*fct)(param, xi, yj, fxy);

				sum_2 += fxy * (*matrix_2)(j,1) * size_2 * exp(-0.5 * yj * yj);
			}

			res += sum_2 * (*matrix)(i,1) * size * exp(-0.5 * xi * xi);
		}

		res	/=	TWOPI;

		break;
		
		}
	case qGAUSS_HERMITE:
		{
		it = itsMap.find(itsIntegrationStep); 
		if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep];

		if (itsIntegrationStep == itsIntegrationStep_2)
			matrix_2	=	matrix;
		else
		{
			it_2 = itsMap.find(itsIntegrationStep_2); 
			if (it_2!=itsMap.end()) matrix = &itsMap[itsIntegrationStep_2];
		}

		if (matrix == NULL)
		{
			ComputeAbscissasAndWeightsGaussHermite(itsIntegrationStep, itsX, itsW);

			ComputeAbscissasAndWeightsGaussHermite(itsIntegrationStep_2, itsX_2, itsW_2);
		}

		for (i=0;i<itsIntegrationStep;i++)
		{
			*nointegrationcoef_x = i;
			xi	=	(*matrix)(i,0);
			
			sum_2	=	0.0;

			for (j=0; j <itsIntegrationStep_2; j++)
			{
				*nointegrationcoef_y	=	j;
				yj	=	(*matrix_2)(j,0);

				(*fct)(param, xi, yj, fxy);

				sum_2 += fxy * (*matrix_2)(j,1);
			}

			res += sum_2 * (*matrix)(i,1);
		}

		res	*=	ONEOVERPI;

		break;
		}
	default:
/*		
		// TRAPEZIUM Case
		double AbsTemp		=	x2;
		double AbsTemp_2	=	y2;

		double Pas		=	size / NbStepTrapeze;
		double Pas_2	=	size_2 / NbStepTrapeze_2;
		double coef = 0.;

			// On calcule les valeurs aux bornes de l'intervalle d'intégration
			// Borne Sup
			(*fct)(param, AbsTemp, AbsTemp_2, fxy);
				coef = exp(-0.5 * AbsTemp * AbsTemp)/2.;
				res1 += fx1*coef;
			
			// Borne Inf
				AbsTemp = x1;
				coef = exp(-0.5 * AbsTemp * AbsTemp)/2.;
				(*fct)(param, AbsTemp, AbsTemp_2, fx1);
				res1 += fx1*coef;
			
			// On calcule toutes les valeurs intermediaires
			for (i=1;i<NbStepTrapeze;i++)
			{
				for (j=1;j<NbStepTrapeze_2;j++)
				{
					AbsTemp += Pas;
					coef = exp(- 0.5 * AbsTemp * AbsTemp);
					(*fct)(param, AbsTemp, AbsTemp_2, fx1, fx2);
					res1 += fx1*coef;
					res2 += fx2*coef;
				}
			}


			res1	*=	Pas / TWO_PI;
*/
			break;
	}

}

void  ICM_Integrator_Dim_2::Integrate(
							const double& x1,
							const double& x2,
								const double& y1,
								const double& y2,
							void (*fct)(void*, double, double, double&, double&),
							AddressVector* param,
							double& res1,
							double& res2)
{
/*	int	i;
	double fx1, fx2;
	
	// Data Trapeze
	int NbStepTrapeze = itsIntegrationStep - 1;
	std::map<int,ICM_QMatrix<double> >::const_iterator it;
	ICM_QMatrix<double>* matrix=NULL;

	if (fct == NULL)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
			"ERROR: Null function for Integrator!");
	}

	int	n	=	itsIntegrationStep;

	//On load les valeurs dans la map en fonction du step d'intégration
	LoadValues(itsIntegrationStep);

	//On verifie que la méthode numérique est compatible avec le nombre de pas employé
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegratorType == qGAUSS_LEGENDRE))
		itsIntegratorType = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegratorType == qGAUSS_HERMITE))
		itsIntegratorType = qGAUSS_LEGENDRE;

	double	size	=	x2 - x1;
	double	xi =0., coef =0.;
	res1	=	0.0;
	res2	=	0.0;

	int* nointegrationcoef = (int*)param->GetLast();
	
	switch (itsIntegratorType)
	{
	case qGAUSS_LEGENDRE:
		{
		it = itsMap.find(itsIntegrationStep); 
		if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep];
		
		if (matrix)
		{
			for (i=0;i<itsIntegrationStep;i++)
			{
				*nointegrationcoef = i;
				xi	=	x1 + size * (*matrix)(i,0);
				(*fct)(param, xi, fx1, fx2);
			
				coef	=	size * exp(- 0.5 * xi * xi);
				res1 += fx1 * (*matrix)(i,1) * coef;
				res2 += fx2 * (*matrix)(i,1) * coef;
			}
	
			res1	/=	SQRT2PI;
			res2	/=	SQRT2PI;
		}
		else
		{		
		
			ComputeAbscissasAndWeightsGaussLegendre(0.0, 1.0, itsIntegrationStep);
			for (i=0;i<n;i++)
			{
				*nointegrationcoef = i;
				xi	=	x1 + size * itsX[i];
				(*fct)(param, xi, fx1, fx2);
			
				coef	=	size * exp(- 0.5 * xi * xi);
				res1 += fx1 * itsW[i] * coef; 
				res2 += fx2 * itsW[i] * coef;
			}

			res1	/=	SQRT2PI;
			res2	/=	SQRT2PI;

		}
		break;
		}

	case qGAUSS_HERMITE:
		{
		it = itsMap.find(itsIntegrationStep); 
		if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep];
		
		if (matrix)
		{
			for (i=0;i<itsIntegrationStep;i++)
			{
				*nointegrationcoef = i;
				(*fct)(param, (*matrix)(i,0), fx1, fx2);
		
				res1 += fx1 * (*matrix)(i,1);
				res2 += fx2 * (*matrix)(i,1);
			}

			res1	*=	ONEOVERSQRTPI;
			res2	*=	ONEOVERSQRTPI;
		}
		else
		{
	
			ComputeAbscissasAndWeightsGaussHermite(itsIntegrationStep);

			for (i=0;i<n;i++)
			{
				*nointegrationcoef = i;
				xi	=	itsX[i];
				(*fct)(param, xi, fx1, fx2);
			
				res1 += fx1 * itsW[i]; 
				res2 += fx2 * itsW[i];
			}

			res1	*=	ONEOVERSQRTPI;
			res2	*=	ONEOVERSQRTPI;
			break;
		}
		break;
		}
	
	default:
			double AbsTemp = x2;
			double Pas = size/NbStepTrapeze;
			double coef = 0.;

			// On calcule les valeurs aux bornes de l'intervalle d'intégration
			// Borne Sup
				(*fct)(param, AbsTemp, fx1, fx2);
				coef = exp(-0.5 * AbsTemp * AbsTemp)/2.;
				res1 += fx1*coef;
				res2 += fx2*coef;
			
			// Borne Inf
				AbsTemp = x1;
				coef = exp(-0.5 * AbsTemp * AbsTemp)/2.;
				(*fct)(param, AbsTemp, fx1, fx2);
				res1 += fx1*coef;
				res2 += fx2*coef;
			
			// On calcule toutes les valeurs intermediaires
			for (i=1;i<NbStepTrapeze;i++)
			{
				AbsTemp += Pas;
				coef = exp(- 0.5 * AbsTemp * AbsTemp);
				(*fct)(param, AbsTemp, fx1, fx2);
				res1 += fx1*coef;
				res2 += fx2*coef;
			}


			res1	*=	Pas/SQRT2PI;
			res2	*=	Pas/SQRT2PI;

			break;
		}
		*/
}


//Integration à partir d'une fonction donnée par un vecteur de points
//La taille du vecteur doit correspondre à nb de points de la méthode 
//Intégration effectuée des membre LBound a Ubound si Hors Gauss Hermite
/*
double  ICM_Integrator_Dim_2::IntegrateVector(const std::vector<double> fx)
{
	//Vector Size Test
	if (fx.size() != itsIntegrationStep)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
			"ERROR: Vector Size different from Integrator Size");
	}

	std::map<int,ICM_QMatrix<double> >::const_iterator it;
	ICM_QMatrix<double>* matrix=NULL;

	int	i;
	double res=0.0;
	double xi;
	double size	= itsUpBound - itsLowBound;
	int	n =	itsIntegrationStep;

	//On load les valeurs dans la map en fonction du step d'intégration
	LoadValues(itsIntegrationStep);

	//On verifie que la méthode numérique est compatible avec le nombre de pas employé
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegratorType == qGAUSS_LEGENDRE))
		itsIntegratorType = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegratorType == qGAUSS_HERMITE))
		itsIntegratorType = qGAUSS_LEGENDRE;

	switch (itsIntegratorType)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	it = itsMap.find(itsIntegrationStep); 
	if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep];
	break;
	default:;
	}

	switch (itsIntegratorType)
	{
	case qGAUSS_LEGENDRE:
		{
		if (matrix)
		{
			for (i=0;i<itsIntegrationStep;i++)
				{
					xi	=	itsLowBound + size * (*matrix)(i,0);
					res += fx[i] * (*matrix)(i,1) * size * exp(-0.5 * xi * xi);
				}

				res	/=	SQRT2PI;
		}
		else
		{
				ComputeAbscissasAndWeightsGaussLegendre(0.0, 1.0, itsIntegrationStep);

				for (i=0;i<n;i++)
				{
					xi	=	itsLowBound + size * itsX[i];
					res += fx[i] * itsW[i] * size * exp(-0.5 * xi * xi);
				}

				res	/=	SQRT2PI;
		}
		break;
		}
	case qGAUSS_HERMITE:
		{
		if (matrix)
		{
				for (i=0;i<itsIntegrationStep;i++)
				{
					res += fx[i] * (*matrix)(i,1);
				}
				res	*=	ONEOVERSQRTPI;
		}
		else	
		{
				ComputeAbscissasAndWeightsGaussHermite(itsIntegrationStep);

				for (i=0;i<n;i++)
				{
					xi	=	itsX[i];
					res += fx[i] * itsW[i] ;
				}

				res	*=	ONEOVERSQRTPI;
		}
		break;
		}
	default: // Cas TRAPEZE
		// -------- Rajout méthode TRAPEZE
		// onin le 14/03/05
		int NbStepTrapeze = itsIntegrationStep - 1;
		double AbsTemp = 0.;
		double coef = size/NbStepTrapeze;

		// On calcule les valeurs aux bornes de l'intervalle d'intégration
		// Borne sup
		AbsTemp = itsUpBound;
		res += fx[NbStepTrapeze-1]*exp(-0.5 * AbsTemp * AbsTemp)/2.;

		// Borne Inf
		AbsTemp = itsLowBound;
		res += fx[0]*exp(-0.5 * AbsTemp * AbsTemp)/2.;	
		
		// On calcule toutes les valeurs intermediaires
		for (i=1;i<NbStepTrapeze;i++)
		{
			AbsTemp += coef;
			res += fx[i]*exp(-0.5 * AbsTemp * AbsTemp);
		}

		res	*= coef/SQRT2PI;
		break;
	}
	return res;
}
*/

//Retourne la vraie valeur de l'abcisse 
//En fonction des bornes Lbound et Ubound si GL ou TR
double ICM_Integrator_Dim_2::GetRealAbscissa_2(const int& index)
{
	double res = 0., size =0., coef =0.;
	switch (itsIntegratorType_2)
	{
	//Gauss Legendre
	case qGAUSS_LEGENDRE:
		size	= itsUpBound_2 - itsLowBound_2;
		res	=	itsLowBound_2 + size * GetAbscissa_2(index);
		break;
	//Gauss Hermite
	case qGAUSS_HERMITE:
		res = GetAbscissa_2(index) * SQRT2;
		break;
	//Trapeze
	default:
		int NbStepTrapeze = itsIntegrationStep_2 - 1;
		size	= itsUpBound_2 - itsLowBound_2;
		double coef = size/NbStepTrapeze;
		if (index <0 || index>=itsIntegrationStep_2)
			res = 0.;
		else
			res = itsLowBound_2 + ((double) index * coef);
		break;
	}
	return res;
}

double ICM_Integrator_Dim_2::GetAbscissa_2(const int& index)
{
	std::map<int,ICM_QMatrix<double> >::const_iterator it;
	ICM_QMatrix<double>* matrix=NULL;

	//On load les valeurs dans la map en fonction du step d'intégration
	LoadValues(itsIntegrationStep_2);

	switch (itsIntegratorType_2)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	it = itsMap.find(itsIntegrationStep_2); 
	if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep_2];
	break;
	default:;
	}

	switch (itsIntegratorType_2)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	
		if (matrix)
		{
			if (index <0 || index>=itsIntegrationStep_2)
				return 0.;
			else
				return (*matrix)(index, 0);
		}
		break;
	default:
			
			if (index <0 || index>=itsIntegrationStep_2)
				return 0.;
			else
				return itsX_2[index];
	}

	return 0.;
}
	
double ICM_Integrator_Dim_2::GetWeight_2(const int& index)
{
	std::map<int,ICM_QMatrix<double> >::const_iterator it;
	ICM_QMatrix<double>* matrix=NULL;

	//On load les valeurs dans la map en fonction du step d'intégration
	LoadValues(itsIntegrationStep_2);

	switch (itsIntegratorType_2)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	it = itsMap.find(itsIntegrationStep_2); 
	if (it!=itsMap.end()) matrix = &itsMap[itsIntegrationStep_2];
	break;
	default:;
	}

	switch (itsIntegratorType_2)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	
		if (matrix)
		{
			if (index <0 || index>=itsIntegrationStep_2)
				return 0.;
			else
				return (*matrix)(index, 1);
		}
		break;
	default:
			
			if (index <0 || index>=itsIntegrationStep_2)
				return 0.;
			else
				return 0.;
	}

	return 0.;

}


ReturnCode  ICM_Integrator_Dim_2 :: ComputeAbscissasAndWeightsGaussLegendre(const double& x1, 
																	  const double& x2, 
																	  const int& N,
																	  double*& TheX,
																	  double*& TheW
																	  )
/* Given the lower and upper limits of integration x1 and x2, and given N, this routine returns
arrays x[1..N] and w[1..N] of length N, containing the abscissas and weights of the Gauss-
Legendre N-point quadrature formula.*/

// int gauleg(double x1, double x2, double *x, double *w, int N)
{
	int m, j, i;
	double z1, z, xm, xl, pp, p3, p2, p1; //High precision is a good idea for this routine.

	if (N < 1) return RetNotOk;

	ICM_QMatrix<double> Matrix(N,2);

	if (TheX)
		delete[] TheX;
	TheX = (double *)malloc(N*sizeof(double));
	if (TheW)
		delete[] TheW;
	TheW = (double *)malloc(N*sizeof(double));

	m = (N+1)/2; //The roots are symmetric in the interval, so we only have to nd half of them. 
	xm = 0.5 * (x2 + x1);
	xl = 0.5 * (x2 - x1);

	//Loop over the desired roots.
	for (i=1;i<=m;i++) 
	{ 
		z=cos(3.141592654*(i-0.25)/(N+0.5));
		//Starting with the above approximation to the ith root, we enter the main loop of re nement by Newton's method.
		do 
		{
			p1=1.0;
			p2=0.0;

			//Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
			for (j=1;j<=N;j++) 
			{ 
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j-1.0) * z * p2 - (j - 1.0) * p3) / j;
			}

			//p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp = N * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp; //Newton's method.
		} 
		while (fabs(z-z1) > ZEPS);

		Matrix(i-1,0) = TheX[i-1] = xm - xl * z;								//Scale the root to the desired interval,
		Matrix(N-i,0) = TheX[N-i] = xm + xl * z;							//and put in its symmetric counterpart.
		Matrix(i-1,1) = TheW[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);		//Compute the weight
		Matrix(N-i,1) = TheW[N-i] = TheW[i-1];									//and its symmetric counterpart.
	}

	itsMap[N]=Matrix;	

	return RetOk;
}


ReturnCode  ICM_Integrator_Dim_2 :: ComputeAbscissasAndWeightsGaussHermite(const int& N,
																	  double*& TheX,
																	  double*& TheW)
{
	int i,its,j,m;
	double p1,p2,p3,pp,z,z1;	// High precision is a good idea for this routine.

	if (N < 1) return RetNotOk;

	ICM_QMatrix<double> Matrix(N,2);

	if (TheX)
		delete[] TheX;
	TheX = (double *)malloc(N*sizeof(double));
	if (TheW)
		delete[] TheW;
	TheW = (double *)malloc(N*sizeof(double));

	m=(N+1)/2;
	// The roots are symmetric about the origin, so we have to .nd only half of them.
	for (i=1;i<=m;i++)
	{			
		//	Loop over the desired roots.
		if (i == 1)
		{
			// Initial guess for the largest root.
			z=sqrt((double)(2*N+1))-1.85575*pow((double)(2*N+1),-0.16667);
		}
		else if (i == 2)
		{
			// Initial guess for the second largest root.
			z -= 1.14*pow((double)N,0.426)/z;
		}
		else if (i == 3)
		{
			// Initial guess for the third largest root.
			z=1.86*z-0.86*itsX[0];
		}
		else if (i == 4)
		{
			// Initial guess for the fourth largest root.
			z=1.91*z-0.91*itsX[1];
		}
		else
		{
			// Initial guess for the other roots.
			z=2.0*z-itsX[i-3];
		}

		for (its=1;its<=MAXIT_GH;its++)
		{
			// Refinement by Newton’s method.
			p1=PIM4_GH;
			p2=0.0;
			for (j=1;j<=N;j++)
			{
				// Loop up the recurrence relation to get
				// the Hermite polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
			}
			// p1 is now the desired Hermite polynomial. We next compute pp, its derivative, by
			// the relation (4.5.21) using p2, the polynomial of one lower order.
			pp=sqrt((double)2*N)*p2;
			z1=z;
			z=z1-p1/pp;		// Newton’s formula.
			if (fabs(z-z1) <= EPS_GH) break;
		}

		if (its > MAXIT_GH) return NULL;	// too many iterations

		Matrix(i-1,0) = TheX[i-1]=z;					// Store the root
		Matrix(N-i,0) = TheX[N-i] = -z;			// and its symmetric counterpart.
		Matrix(i-1,1) = TheW[i-1]=2.0/(pp*pp);		// Compute the weight
		Matrix(N-i,1) = TheW[N-i]=TheW[i-1];			// and its symmetric counterpart.
	}
	
	itsMap[N]=Matrix;	

	return RetOk;
}

