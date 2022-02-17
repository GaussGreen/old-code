/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Integrator_Dim_2_DIM_2.H
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

/*
	Modif CN : LongShort CDO
	double GetAbscissa(const int& index) : point
	double GetWeight(const int& index) : weight
*/

/*
	Modif CN : Integration à partir d'une fonction donnée par une vecteur
	double IntegrateVector(const vector<double> fx) : point
	
*/

#ifndef _ICM_Integrator_Dim_2_H
#define _ICM_Integrator_Dim_2_H

# include <math.h> 

# include "ICMKernel\cair\types.h"
# include "ICMKernel\glob\icm_enums.h"
# include "ICMKernel\util\icm_qmatrix.h"
# include "ICMKernel\glob\icm_addressvector.h"

# include "ICMKernel\util\icm_integrator.h"

class ICM_Integrator_Dim_2 : public ICM_Integrator
{
	// ------------------------------------------------------------------------------
	// DATA
	// ------------------------------------------------------------------------------

	// Recall Gauss-Hermite Formula is given with the integration formula
	// I(-inf ; +inf) exp(-x^2) * f(x) dx = Sum(i=1..N) w(i) * f(x(i))

	// REMARK, for an Integration in Dimension n, play with vectors<...>
protected:

	int		itsIntegrationStep_2;
	double	itsLowBound_2;
	double	itsUpBound_2;

	qIntegratorChoice	itsIntegratorType_2;

	// Second Integrand
	double	*itsX_2;		// abscissas
	double	*itsW_2;		// weigths

public:

	// Methods to handle with the Data
	inline void	SetIntegrationStep_2(const int& data) {itsIntegrationStep_2 = data;}
	inline int	GetIntegrationStep_2(void) { return itsIntegrationStep_2;}

	inline void	SetLowBound_2(const double& data) {itsLowBound_2 = data;}
	inline double GetLowBound_2(void) { return itsLowBound_2;}

	inline void	SetUpBound_2(const double& data) {itsUpBound_2 = data;}
	inline double GetUpBound_2(void) { return itsUpBound_2;}

	inline void	SetIntegrationType_2(const qIntegratorChoice& data) {itsIntegratorType_2 = data;}
	inline qIntegratorChoice GetIntegrationType_2(void) { return itsIntegratorType_2;}

	void	GetAbscissaVector_2(DoubleVector& data);
	void	GetWeightVector_2(DoubleVector& data);

	// ------------------------------------------------------------------------------
	// Constructors & so on
	// ------------------------------------------------------------------------------
	public :

		ICM_Integrator_Dim_2() 
		{Init();}	

		void BitwiseCopy(const ICM_Integrator_Dim_2* src);
		void Copy(const ICM_Integrator_Dim_2* srcZc);
        ICM_Integrator_Dim_2* Clone() const;

		~ICM_Integrator_Dim_2()	{Reset();}

		void	Init();
		void	Reset();

public :
	
	double	Integrate() {return -1.0;}

	ReturnCode  ComputeAbscissasAndWeightsGaussLegendre(const double& x1, 
													  const double& x2, 
													  const int& N,
													  double*& TheX,
													  double*& TheW
													  );
	
	ReturnCode  ComputeAbscissasAndWeightsGaussHermite(const int& N,
													  double*& TheX,
													  double*& TheW);

	double GetAbscissa_2(const int& index);
	double GetWeight_2(const int& index);

public:
	
	void	Integrate(const double& x1,
					const double& x2,
					const double& y1,
					const double& y2,
					void (*fct)(void*, double, double, double&),
					AddressVector* param,
					double& res);

	void	Integrate(const double&	x1,
					  const double&	x2,
						const double& y1,
						const double& y2,
					  void (*fct)(void*, double, double, double&, double&),
					  AddressVector* param,
					  double& res1,
					  double& res2);
	//Integration
//	double IntegrateVector(const std::vector<double> fx);	
	double GetRealAbscissa_2(const int& index);

	ICM_QMatrix<double>* GetCoefMatrix(const int& nbstep)
	{	LoadValues(nbstep);
		std::map<int,ICM_QMatrix<double> >::const_iterator it;
		it = itsMap.find(itsIntegrationStep_2);
		if (it!=itsMap.end()) return &itsMap[itsIntegrationStep_2];

		return NULL; }


};

# endif

