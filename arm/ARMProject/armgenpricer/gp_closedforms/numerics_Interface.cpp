/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file numerics_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/bivariate_normal.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/incompletebeta.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/lambert_function.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/heston_calibration.h"
#include "gpclosedforms/smile_shiftedlognormal.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/smile_bisabr.h"
#include "gpclosedforms/student_copula.h"
#include "gpclosedforms/whittaker.h"
#include "gpclosedforms/bessel.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/tridiagonalsolve.h"
#include "gpclosedforms/eigenvalues.h"
#include "gpclosedforms/smile_2logsmiled.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/typedef.h"
#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/heston_calibration.h"


using namespace std;

CC_BEGIN_NAMESPACE(ARM)

double Export_Bivariate(double x, double y, double rho, int p, int q)
{
	if (p == 0 && q == 0)
		return NormalCDF( x, y,  rho);
	else if (p==1 && q==0)
		return Normal_X_Expectation(x, y, rho);
	else if (p==0 && q==1)
		return Normal_Y_Expectation(x, y, rho);
	else if (p==1 && q==1)
		return Normal_XY_Expectation(x, y, rho);
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Bivariate:  ");

}

double Export_Gamma(double a)
{
	return gamma(a);
}

double Export_Lambert_Function(double z)
{
		if (z <= -ARM_NumericConstants::ARM_INVERSE_E) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"z should be > -1/e in lambertfunc" );

	return lambertfunc(z);
}


double Export_IncompleteBeta(const double& a, const double& b, const double& x)
{
	return IncompleteBeta(a,  b,  x);
}

double Export_InverseIncompleteBeta(const double& a, const double& b, const double& x0,const double& x)
{
	return IncompleteBetaInverse(a,  b,  x0,x);
}


double  Export_Hypergeometric2F1( double  a, double  b, double  c, double  z)
{
	if (z>=1.0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Hypergeometric2F1  : abs(z)>1 !  ");
	}
		return Hypergeometric2F1(a,b,c,z);
}


GaussLegendre_Coefficients* Export_GaussLegendre_Coefficients(double a,double b,int n)
{
	GaussLegendre_Coefficients* ptr=new GaussLegendre_Coefficients(n,a,b);
	return ptr;
}

GaussHermite_Coefficients* Export_GaussHermite_Coefficients(int n)
{
	GaussHermite_Coefficients* ptr=new GaussHermite_Coefficients(n);
	return ptr;
}

GaussLaguerre_Coefficients* Export_GaussLaguerre_Coefficients(double alpha ,int n)
{
	GaussLaguerre_Coefficients* ptr=new GaussLaguerre_Coefficients(alpha,n);
	return ptr;
}



GeneralizedHeston_ParameterSet*   Export_GeneralizedHeston_Model_Calibrate(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,double f,double t,int nbsteps,int algorithm,
												double V0,double omega,double theta,double ksi,double rho,double muJ,double sigmaJ,double lambda)
{
	int arg1=K_Vec->size();
	int arg2=ImpVol_Vec->size();
	if (arg1!=arg2)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_GeneralizedHeston_Model_Calibrate::value  : nb element of K_Vec and ImpVol_Vec should be the same !");
	}
	GeneralizedHeston_ParameterSet* setptr= GeneralizedHeston_CalibrateToSmile(K_Vec,ImpVol_Vec,f,t,nbsteps,algorithm,
												  V0,omega, theta, ksi, rho, muJ, sigmaJ, lambda);
	return setptr;

}

GeneralizedHeston_ParameterSet*   Export_GeneralizedHeston_Model_Calibrate_NoJump(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,double f,double t,double muJ0,double sigmaJ0,double lambda0,int nbsteps,int algorithm,
												double V0,double omega,double theta,double ksi,double rho)
{
	int arg1=K_Vec->size();
	int arg2=ImpVol_Vec->size();
	if (arg1!=arg2)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_GeneralizedHeston_Model_Calibrate::value  : nb element of K_Vec and ImpVol_Vec should be the same !");
	}
	GeneralizedHeston_ParameterSet* setptr= GeneralizedHeston_CalibrateToSmile(K_Vec,ImpVol_Vec,f,t,muJ0,sigmaJ0,lambda0,nbsteps,algorithm,
												  V0,omega, theta, ksi, rho);
	return setptr;

}



double Export_ImcompleteBeta_Inverse(double a, double b, double x)
{
	return IncompleteBetaInverse(a,b,x);
}

double Export_Student_QIntegral(double rank, double b, double x)
{
	ArgumentList a(0.5,rank);
	return StudentCopula::Q_Integral(a,b,x);
}



double  Export_Hypergeometric_Whittaker_W (double  a,double  b,double  z)
{
	return Hypergeometric_Whittaker_W (  a,  b,  z);
}

double  Export_Hypergeometric_Whittaker_M (double  a,double  b,double  z)
{
	return Hypergeometric_Whittaker_M (  a,  b,  z);
}


double Export_Bessel_I(double nu, double x)
{
	return bessel_fractional_I( nu,  x);
}


double Export_Bessel_J(double nu, double x)
{
	return bessel_fractional_J( nu,  x);
}

double Export_Bessel_K(double nu, double x)
{
	return bessel_fractional_K( nu,  x);
}

double Export_Bessel_Y(double nu, double x)
{
	return  bessel_fractional_Y( nu,  x);
}


double Export_Hypergeometric_Appell(double a,double b1,double b2,double c ,double x,double xim, double y,double yim,int nb)
{
	complex<double> ca(a,0);
	complex<double> cb1(b1,0);
	complex<double> cb2(b2,0);
	complex<double> cc(c,0);
	complex<double> cx(x,xim);
	complex<double> cy(y,yim);
	
	complex<double> res=HypergeometricAppellF1(ca,cb1,cb2,cc,cx,cy,nb);
	return real(res);;
}


/*
a :diagonale inferieure
b: diagonale
c: diagonale superieure
r: vecteur constant
u: resultat
gam: vecteur tampon
  */
vector<double>* Export_Util_Trigonal_Solve(ARM_GP_T_Vector<double>* A_Vec,ARM_GP_T_Vector<double>* B_Vec,ARM_GP_T_Vector<double>* C_Vec, ARM_GP_T_Vector<double>* R_Vec)
{
	int n=B_Vec->size();
	vector<double>* result=new vector<double>(n);
	vector<double> tampon(n);
	const vector<double> a_Vec=A_Vec->GetValues();
	const vector<double> b_Vec=B_Vec->GetValues();
	const vector<double> c_Vec=C_Vec->GetValues();
	const vector<double> r_Vec=R_Vec->GetValues();
	tridiagonalsolve(a_Vec,b_Vec,c_Vec,r_Vec,*result,tampon);
	return result;
}


void Export_Eigenvalues4(double rho12,double rho13, double rho14,double rho23,double rho24,double rho34,double* e1,double* e2,double* e3,double* e4)
{
	eigenvalues4( rho12, rho13,  rho14, rho23, rho24, rho34, e1, e2, e3, e4);
	return;
}

void Export_Eigenvalues3(double rho12,double rho13,double rho23,double* e1,double* e2,double* e3)
{
	eigenvalues3( rho12, rho13, rho23, e1, e2, e3);
	return;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
