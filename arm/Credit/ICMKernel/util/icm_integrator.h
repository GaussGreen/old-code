/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_INTEGRATOR.H
	PROJECT:	UTIL
	
	DESCRIPTION:	Numerical Integrator

  -----------------------------------------------------------------

 	CREATION:	October 8, 2004

	LAST MODIF:	October 8, 2004
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

/*
	Modif CN : LongShort CDO
	double GetAbscissa(const int& index) : point
	double GetWeight(const int& index) : weight
*/

/*
	Modif CN : Integration à partir d'une fonction donnée par une vecteur
	double IntegrateVector(const vector<double> fx) : point
	
*/

#ifndef _ICM_INTEGRATOR_H
#define _ICM_INTEGRATOR_H

# include <math.h> 

# include "ICMKernel\cair\types.h"
# include "ICMKernel\glob\icm_enums.h"
# include "ICMKernel\util\icm_qmatrix.h"
#include "ICMKernel\glob\icm_addressvector.h"
#include <map>

class ICM_Integrator
{
	// ------------------------------------------------------------------------------
	// DATA
	// ------------------------------------------------------------------------------

	// Recall Gauss-Hermite Formula is given with the integration formula
	// I(-inf ; +inf) exp(-x^2) * f(x) dx = Sum(i=1..N) w(i) * f(x(i))

protected:

	typedef std::map<int,ICM_QMatrix<double> > map_integrator; 
	static	map_integrator itsMap;

	int		itsIntegrationStep;
	double	itsLowBound;
	double	itsUpBound;

	qIntegratorChoice	itsIntegratorType;

	double	*itsX;		// abscissas
	double	*itsW;		// weigths

public:

	void LoadValues(const int& nb);
	void FullLoadValues();

	// Methods to handle with the Data
	inline void	SetIntegrationStep(const int& data) {itsIntegrationStep = data; }
	inline int	GetIntegrationStep(void) { return itsIntegrationStep;}

	inline void	SetLowBound(const double& data) {itsLowBound = data; }
	inline double GetLowBound(void) { return itsLowBound;}

	inline void	SetUpBound(const double& data) {itsUpBound = data; }
	inline double GetUpBound(void) { return itsUpBound;}

	inline void	SetIntegrationType(const qIntegratorChoice& data) {itsIntegratorType = data; }
	inline qIntegratorChoice GetIntegrationType(void) { return itsIntegratorType;}

	void	GetAbscissaVector(DoubleVector& data);
	void	GetWeightVector(DoubleVector& data);

	// ------------------------------------------------------------------------------
	// Constructors & so on
	// ------------------------------------------------------------------------------
	public :

		ICM_Integrator() 
		{Init();}	

		void BitwiseCopy(const ICM_Integrator* src);
		void Copy(const ICM_Integrator* srcZc);
        ICM_Integrator* Clone() const;

		~ICM_Integrator()	{Reset();}

		void	Init();
		void	Reset();

public :
	
	ReturnCode	ComputeAbscissasAndWeightsGaussLegendre(const double& x1, const double& x2, const int& n);

	ReturnCode	ComputeAbscissasAndWeightsGaussHermite(const int& n);

	double	Integrate();

	double GetAbscissa(const int& index);
	double GetWeight(const int& index);

public:
	
	void  Integrate(const double& x1,
					const double& x2,
					void (*fct)(void*, double, double&),
					AddressVector* param,
					double& res);

	void	Integrate(const double&	x1,
					  const double&	x2,
					  void (*fct)(void*, double, double&, double&),
					  AddressVector* param,
					  double& res1,
					  double& res2);
	//Integration
	double IntegrateVector(	const std::vector<double> fx);	
	double GetRealAbscissa(const int& index);

	ICM_QMatrix<double>* GetCoefMatrix(const int& nbstep)
	{	LoadValues(nbstep);
		std::map<int,ICM_QMatrix<double> >::const_iterator it;
		it = itsMap.find(itsIntegrationStep);
		if (it!=itsMap.end()) return &itsMap[itsIntegrationStep];

		return NULL; }

	// UTILS
	void	CheckNumericalMethod(int TheIntegrationStep, qIntegratorChoice& TheIntegratorType);


protected:

	// GAUSS LEGENDRE
	static double GaussLegendreCoeff_11[11][2];
	static double GaussLegendreCoeff_21[21][2];
	static double GaussLegendreCoeff_51[51][2];
	static double GaussLegendreCoeff_101[101][2];
	static double GaussLegendreCoeff_201[201][2];
	static double GaussLegendreCoeff_501[501][2];
	static double GaussLegendreCoeff_901[901][2];

	// GAUSS HERMITE
	static double HermiteCoeff_20[20][2];
	static double HermiteCoeff_40[40][2];
	static double HermiteCoeff_60[60][2];
	static double HermiteCoeff_100[100][2];

};

# endif

