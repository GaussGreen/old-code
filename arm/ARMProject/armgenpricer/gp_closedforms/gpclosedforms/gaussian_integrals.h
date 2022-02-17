/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_integrals.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_GAUSSIAN_INTEGRALS_H
#define _GP_CF_GAUSSIAN_INTEGRALS_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/long_double.h"
#include <complex>
#include <vector>
#include <map>

#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"
#include "gpclosedforms/inverse.h"

using std::complex;

using std::vector;

CC_BEGIN_NAMESPACE(ARM)

long double atanh(long double x);

complex<long double> atanh(complex<long double> y);

complex<long double> arctan(complex<long double> y);


long double acoth(long double x);

class Orthogonal_Polynomial_Coefficients
{
public:
	Orthogonal_Polynomial_Coefficients(int n);
	virtual ~Orthogonal_Polynomial_Coefficients() {};

	inline std::vector<double> get_weights() const { return weigths;}
	inline std::vector<double> get_points() const{return points;}
	inline long double get_weight(int n) const {return weigths[n];}
	inline long double get_point(int n) const {return points[n];}
	inline int get_order() const { return points.size(); }

	inline void set_points(int i,double x){points[i]=x;}
	inline void set_weigths(int i,double x){weigths[i]=x;}
	inline void set_points(std::vector<double> x){weigths=x;}
	inline void set_weigths(std::vector<double> x){points=x;}

	void fill( Orthogonal_Polynomial_Coefficients coeff){ weigths= coeff.get_weights(); points= coeff.get_points();}

	protected:
	std::vector<double> weigths;
	std::vector<double> points;
	Orthogonal_Polynomial_Coefficients* clone();
};
/////////////////////////////////////////////////////////////////////////
///
/// Gauss hermite coefficients: can be used to integrales like
/// Integral{-infinity, +infinity} of exp(-x^2)*f(x)
///
/// To compute integrals like exp(-x^2/2)/Sqrt[2Pi] * f(x) , you must multiply
///  evry point by sqrt(2)  and the final summation by 1/sqrt(Pi)
///
/////////////////////////////////////////////////////////////////////////

class GaussHermite_Coefficients_repository
{
public:
	typedef std::map<int,Orthogonal_Polynomial_Coefficients*> GaussMap;
	static	GaussHermite_Coefficients_repository*	get_instance()
	{
		if(!instance) instance=new GaussHermite_Coefficients_repository;
		return instance;
	}
	Orthogonal_Polynomial_Coefficients* get_coefficients(int n)
	{
		GaussMap::iterator iter=repository.find(n);
		if(iter != repository.end())
		{
			return iter->second;
		}
		else
		{
			return NULL;
		}
	}
	GaussMap* get_repository()
	{ return &repository;}
private:
	static	GaussHermite_Coefficients_repository* instance;
	GaussHermite_Coefficients_repository()
	{}

	~GaussHermite_Coefficients_repository()
	{}

	GaussMap repository;
	
	
};


class GaussHermite_Coefficients : public Orthogonal_Polynomial_Coefficients
{
	public:
	GaussHermite_Coefficients(int n) ;
};

class GaussLegendre_Coefficients;
/////////////////////////////////////////////////////////////////////////
///
/// ReducedGauss hermite coefficients: can be used to integrales like
/// Integral{-infinity, +infinity} of exp(-x^2/2) f(x)
///
/// you must multiply
///  evry point by sqrt(2)  and the final summation by 1/sqrt(2 Pi)
///
///  with the constructor (GaussLegendre_Coefficients* glcoeffs_ptr, double x1,double x2)
///  we compute integrals like Integral{a, b} of exp(-x^2/2) f(x)
/////////////////////////////////////////////////////////////////////////

class ReducedGaussHermite_Coefficients : public Orthogonal_Polynomial_Coefficients
{
	public:
	ReducedGaussHermite_Coefficients(int n) ;
	ReducedGaussHermite_Coefficients(GaussLegendre_Coefficients* glcoeffs_ptr, double x1,double x2);
};

/////////////////////////////////////////////////////////////////////////
///
/// Gauss laguerre coefficients: can be used to integrales like
/// Integral{0, +infinity} of x^alpha*exp(-x)*f(x)
///
/////////////////////////////////////////////////////////////////////////

class GaussLaguerre_Coefficients : public Orthogonal_Polynomial_Coefficients
{
	public:
	GaussLaguerre_Coefficients(long double alpha, int n) ;
};


/////////////////////////////////////////////////////////////////////////
///
/// Stratified (Hermite-Legendre) coefficients: can be used to integrales like
/// Integral{a, b} of exp(-x^2/2)/sqrt(2pi) *f(x)
///
/////////////////////////////////////////////////////////////////////////


class GaussStratifiedHermiteLegendre_Coefficients: public Orthogonal_Polynomial_Coefficients
{
	public:
		GaussStratifiedHermiteLegendre_Coefficients(GaussLegendre_Coefficients* coefs,int p,double a,double b);
};


/////////////////////////////////////////////////////////////////////////
///
/// Gauss legendre coefficients: can be used to integrales like
/// Integral{-1, +1} of f(x)
///
/////////////////////////////////////////////////////////////////////////

class GaussLegendre_Coefficients_repository
{
public:
	typedef std::map<int,Orthogonal_Polynomial_Coefficients*> GaussMap;
	static	GaussLegendre_Coefficients_repository*	get_instance()
	{
		if(!instance) instance=new GaussLegendre_Coefficients_repository;
		return instance;
	}
	Orthogonal_Polynomial_Coefficients* get_coefficients(int n)
	{
		GaussMap::iterator iter=repository.find(n);
		if(iter != repository.end())
		{
			return iter->second;
		}
		else
		{
			return NULL;
		}
	}
	GaussMap* get_repository()
	{ return &repository;}
private:
	static	GaussLegendre_Coefficients_repository* instance;
	GaussLegendre_Coefficients_repository()
	{}

	~GaussLegendre_Coefficients_repository()
	{
		GaussMap::iterator it;

		for (it =  repository.begin(); it !=  repository.end(); it++)
			if( it->second ) { delete ( it->second); it->second= NULL;
		}
	}

	GaussMap repository;
	
	
};

class GaussLegendre_Coefficients : public Orthogonal_Polynomial_Coefficients
{
	public:
	GaussLegendre_Coefficients(int n) ;
	virtual ~GaussLegendre_Coefficients() {};

	GaussLegendre_Coefficients(int n,double x1,double x2) ;
	GaussLegendre_Coefficients(GaussLegendre_Coefficients* coefs,double x1,double x2);
	void AdaptCoefficients(GaussLegendre_Coefficients* coefs2,double x1,double x2) ;

	
	
};

/////////////////////////////////////////////////////////////////////////
///
/// ConstantSteps: can be used to integrales like
/// Integral{0, Xmax} of f(x)
///
/////////////////////////////////////////////////////////////////////////

class ConstantSteps_Coefficients : public Orthogonal_Polynomial_Coefficients
{
	public:
	ConstantSteps_Coefficients(double xmax,int n) ;
};


double LegendreInt(const DoubleToDoubleFunc& func, double a, double b, int nbpts);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

