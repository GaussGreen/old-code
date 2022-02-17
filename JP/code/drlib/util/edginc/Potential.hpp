//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : Potential.hpp
//
//   Description : Class that can be used for minimisation of entropy of density associated to potential
//
//   Author      : Matthias Arnsdorf
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/AtomicArray.hpp"
#include "edginc/Function.hpp"
#include "edginc/Optimizer.hpp"


#ifndef QLIB_POTENTIAL
#define QLIB_POTENTIAL


DRLIB_BEGIN_NAMESPACE



class UTIL_DLL Potential : public MFunctionND
{
public:
	Potential(
			const DoubleArray & c,				// constraints
			const DoubleArrayArray & h,			// constraint functional h[numConstraints][gridPoints]. Assumed PIECEWISE LINEAR between grid points
			const DoubleArray& integrationGrid,	// grid points for integration. Both h and prior need to be defined at the grid points
			const DoubleArray& prior,			// Prior function. Assumed PIECEWISE LINEAR between grid points
			int numConstraints					// number of constraints
			); 

	virtual ~Potential(){};

	/** Potential function */
	virtual void operator()(const DoubleArray&  x,	// (I) input array (n-D)
                            DoubleArray&        f	// (O) function values (1-D)	
							) const;


	/** integrates functions given by array of function values y defined at grid points x
	Assumes y is piecewise linear between grid points (i.e. use Trapeze rule)
	
	inline for performance 
	no validation
	*/
	static double integrate(const DoubleArray & y, 
							const DoubleArray & x, 
							int n // number of points in array to avoid size() call
 							)
							{
		double result = 0;
		for(int i = 1;i<n;i++)	result += (y[i]+y[i-1])/2*(x[i]-x[i-1]);
		return result;
	}



	/** Integrates function f(x) against potential
		i.e. return outIntegral(x) = \int_0^x Phi(x)f(x) dx
		aussume that f(x) is defined at integrationGrid points and that f(x) is piecewise linear between points

		Phi defined by lagrange multipliers

		output Integral is also defined at integrationGrid points */
	void integral_PhiXfX(
		const DoubleArray & lagrangeMults,		// lagrange multipliers
		const DoubleArray & f,					// function to integrat. Piecewise linear between integrationGrid points
		double norm,							// norm of potential
		DoubleArray & outIntegral				// (O) output integral defined at integrationGrid points
		) const;


	

	/** integrates functions of form a(t)Exp(b(t)] where a(t) and b(t) are linear
		Inputs are time interval and values of a and b at end points */
	static double integrate_AtexpBt(double t1, double t2, double a1, double a2, double b1, double b2);
	
	/** integrates functions of form a(t)b(t)Exp(c(t)] where a(t), b(t) and c(t) are linear
		Inputs are time interval and values of a, b and c at end points */
	static double integrate_AtBtexpCt(
		double t1, double t2, 
		double a1, double a2, 
		double b1, double b2,
		double c1, double c2);
	

	/** density corresponding to potential defined at integrationGrid points*/
	DoubleArraySP PDF(const DoubleArray & lag) const;

	/** intergral of density corresponding to potential defined at integrationGrid points*/
	DoubleArraySP CDF(const DoubleArray & lag) const;

	/** returns diff c_n - \int PDF(x)h_n(x) dx */
	DoubleArraySP check(const DoubleArray & lagrangeMults, double norm) const; 

	/** returns norm of potential needed to calculate integral */
	double norm(const DoubleArray & lag) const;

private:

	static const double MyTINY; // A small number.


	// FIELDS 
	
	// constraints are defined by:
	// int \sum_i phi(x)h[i][x] dx = c[i]
	// i labels the constraints and x an element of the integrationGrid
	
	/** constraints */
	DoubleArray c;
	/** constraint functionals */
	DoubleArrayArray  h; //h[constraint][gridPoint]

	/** prior density */
	DoubleArray prior;
	
	/** number of constraints = number of lagrange multipliers */
	int numConstraints;
	/** size of integration grid */
	int gridSize;
	/** integrationGrid */
	DoubleArray integrationGrid;

};
 

DRLIB_END_NAMESPACE

#endif
