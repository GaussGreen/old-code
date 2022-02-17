//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : Potential.cpp
//
//   Description : Potential function that can be minimised to solve constrained
//	               minimisation problem with lagrange multipliers
//
//   Author      : Matthias Arnsdorf
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/Potential.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
const double Potential::MyTINY = 1.0e-10; // A small number.

/** Constructor */
Potential::Potential(
			const DoubleArray & c,				// constraints
			const DoubleArrayArray & h,			// constraint functional h[numConstraints][integrationGridPoints]. Assumed PIECEWISE LINEAR between integrationGrid points
			const DoubleArray& integrationGrid,	// integrationGrid points for integration. Both h and prior need to be defined at the integrationGrid points
			const DoubleArray& _prior,			// Prior function. Assumed PIECEWISE LINEAR between integrationGrid points
			int numConstraints
			) :
MFunctionND(numConstraints,1,*InfiniteRange::createInfiniteRangeArray(numConstraints)), 
			c(c), 
			h(h), 
			prior(_prior), 
			numConstraints(numConstraints),
			integrationGrid(integrationGrid)
{
	static const string method = "Potential::Potential";
	try
	{
		int i; //counter
		gridSize = integrationGrid.size();

		if(gridSize < 1) throw ModelException("gridSize <1");

		/** validation */
		if(h.size() != numConstraints) throw ModelException("Sizes of h != numConstraints");
		if(c.size() != numConstraints) throw ModelException("Size of c != numConstraints");
		if(prior.size() != gridSize) throw ModelException("Size of Prior != gridSize");
		for(i = 0;i<numConstraints;i++) 
		{
			if(h[i].size() != gridSize) 
			{
				throw ModelException("h[" + Format::toString(i) +"] != gridSize");
			}
		}
		
		// calculate norm of prior
		double priorNorm = integrate(prior,integrationGrid,gridSize);
		if(priorNorm <= 0.0) throw ModelException("priorNorm <= 0");

		// check if norm is one otherwise normalise
		if(!Maths::areEqualWithinTol(priorNorm,1.0,MyTINY))
		{
			for(i = 0;i<gridSize;i++) prior[i] /= priorNorm;
		}
		


	} catch(exception& e){
		throw ModelException(e,method);
	}
	

}


//////////////////////////////////////////////////////////////////////////////////////////////////
/** overload () 
	This defines potential as function of the lagrange multipliers x[j]
	The potential is given by:

			\int_x prior[x_i]*[exp(\sum_j x[j](c[j] - h[j][x_i]))  dx

	We are assuming that both h and prior are piecewise linear between integrationGrid points.
	Hence the above is a sum of integrals of a linear function times an exponential of a linear function
	This can be done analytically.
*/
//////////////////////////////////////////////////////////////////////////////////////////////
void Potential::operator()(const DoubleArray&  x,
                    DoubleArray&        f) const
{
	int i,j; //counter
	
	// check size of x
	if(x.size() != numConstraints) throw ModelException("x size != numConstraints");
	// check size of f
	if(f.size() < 1) throw ModelException("f size < 1");
	
	

	// initialise integral
	double integral = 0.0;
	
	/** exponent at start of integrationGrid interval */
	double exponent1 = 0.0;
	for(j = 0;j<numConstraints;j++)
	{
		exponent1 += x[j]*(c[j]-h[j][0]);
	}
	
	/** exponent at end of time grid interval */
	double exponent2 = 0.0;

	// loop through integrationGrid
	for(i = 0; i < gridSize - 1;i++)
	{
		exponent2 = 0.0;
		for(j = 0;j<numConstraints;j++)
		{
			exponent2 += x[j]*(c[j]-h[j][i+1]);
		}

		// update integral
		integral += integrate_AtexpBt(
			integrationGrid[i],		integrationGrid[i+1], 
			prior[i],				prior[i+1], 
			exponent1,				exponent2);

		exponent1 = exponent2; // assign old value of exponent2
		
	}
	
	// return output	
	f[0] = integral;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
/** 
	Integrates function f(x) against density Phi(x)
	i.e. return outIntegral(x) = \int_0^x Phi(x)f(x) dx
	aussume that f(x) is defined at integration grid points and that f(x) is piecewise linear between points

	density Phi defined by lagrange multipliers
 */
//////////////////////////////////////////////////////////////////////////////////////////////////
void Potential::integral_PhiXfX(const DoubleArray & lag,
								const DoubleArray & f, 
								double norm, 
								DoubleArray & outIntegral) const
{
	static const string method = "Potential::integral_PhiXfX";
	try
	{
		int i,j;

		if(lag.size() != numConstraints) throw ("Number of lagrange multipliers not equal to number of constraints");
		if(f.size() != gridSize) throw ("size of f not same as grid size");

		if(outIntegral.size() != gridSize) throw ModelException("outIntegral has incorrect size");
				
		
		
		// calculate outIntegral[x] = 1/norm * \int f(x)prior(x)exp(-\sum_j lag[j]*(h[j][x]-c[j]) dx
		outIntegral[0] = 0.0;

		/** exponent at start of grid interval */
		double exponent1 = 0.0;
		for(j = 0;j<numConstraints;j++)
		{
			exponent1 -= lag[j]*(h[j][0]-c[j]);
		}
		/** exponent at end of time grid interval */
		double exponent2 = 0.0;

		for(i = 0; i < gridSize - 1;i++)
		{
			exponent2 = 0.0;
			for(j = 0;j<numConstraints;j++)
			{
				exponent2 -= lag[j]*(h[j][i+1]-c[j]);
			}


			outIntegral[i+1] = outIntegral[i];

			outIntegral[i+1] += integrate_AtBtexpCt(
				integrationGrid[i], integrationGrid[i+1], 
				prior[i], prior[i+1],
				f[i], f[i+1],
				exponent1, exponent2)/norm;

			exponent1 = exponent2; // assign old value of exponent2
			
		}
		

	} catch(exception& e){
		throw ModelException(e,method);
	}
}


/////////////////////////////////////////////////////
/** return norm */
////////////////////////////////////////////////////
double Potential::norm(const DoubleArray & lag) const
{
	static const string method = "Potential::norm";
	try
	{
		if(lag.size() != numConstraints) throw ("Number of lagrange multipliers not equal to number of constraints");
	
		// calculate norm
		DoubleArray norm(1);
		operator()(lag,norm);

		if(norm[0] < MyTINY) throw ModelException("Normalisation factor = 0");

		return norm[0];	

	} catch(exception& e){
		throw ModelException(e,method);
	}

}

/////////////////////////////////////////////////////////////////////////////////////
/** integrate linear exponentials 
	i.e. functions of form a(t)exp[b(t)]
	where a(t) and b(t) are linear	*/
/////////////////////////////////////////////////////////////////////////////////////
double Potential::integrate_AtexpBt(double t1, double t2, double a1, double a2, double b1, double b2)
{
	static const string method = "Potential::integrate_AtexpBt";
	try{
		if(t2 == t1) return 0;

		

		// First calculate linear coefficients for a(t) and b(t) from given endpoints
		// a(t) = constA + linA*t
		// b(t)= constB + linB*t
		double linA = (a2 - a1)/(t2-t1);
		double linB = (b2 - b1)/(t2-t1);

		double constA = a1 - linA*t1;
		double constB = b1 - linB*t1;

		
		// Want integral of (constA + linA*t)*exp(constB + linB*t) --------------------
		double integral = 0;
		
		// Have special case when linB  is small 
		if(abs(linB) < sqrt(MyTINY))
		{
			double dt1 = t2-t1;
			double dt2 = (t2*t2-t1*t1)/2;
			double A = exp(constB)*constA;
			double B = exp(constB)*linA;
			// first order approximation of integral
			// integral = exp(constB)*constA(t2-t1) + exp(constB)*linA/2*(t2^2-t1^2)
			integral = A*dt1 + B*dt2;
			
			// higher orders
			if(linB != 0) 
			{
				// add second order approximation if linB != 0
				// need to add integral of (constA + linA*t)*exp(constB)*linB*t
				// i.e exp(constB)*constA*linB(t2^2-t1^2)/2 + exp(constB)*linA*linB*(t2^3-t1^3)/3
				double dt3 = (t2*t2*t2-t1*t1*t1)/3;
				integral += linB*(A*dt2 + B*dt3);

				// third order
				double dt4 = (t2*t2*t2*t2-t1*t1*t1*t1)/4;
				integral += linB*linB/2*(A*dt3 + B*dt4);

				// fourth order
				double dt5 = (t2*t2*t2*t2*t2-t1*t1*t1*t1*t1)/5;
				integral += linB*linB*linB/6*(A*dt4+B*dt5);
			}
		}
		// now consider full calculation
		else
		{	
			// write as constA\int exp(b(t2)) dt + linA\int t*exp(b(t1)) dt
			//
			// First term is:
			// constA/linB*[exp(b(t2))-exp(b(t1))]
			//
			// second term is:
			// linA/linB*[(t2*exp(b(t2))-t1*exp(b(t1)) - (exp(b(t2))-exp(b(t1)))/linB]
			
			// check that b1 and b2 don't exceed bounds
			const double maxExponent = 700.0;
			if(b1 > maxExponent) throw ModelException("Internal error: b1 exceeds max bound");
			if(b2 > maxExponent) throw ModelException("Internal error: b2 exceeds max bound");


			double expB1 = exp(b1);
			double expB2 = exp(b2);
			
			double term1 = constA/linB*(expB2 - expB1);

			double lalb = (a2-a1)/(b2-b1); //= linA/linB
			double lalblb = lalb/linB;

			//double term2 = linA/linB*((t2-1/linB)*expB2 - (t1-1/linB)*expB1);
			double term2 = (t2*expB2 - t1*expB1)*lalb + (expB1 - expB2)*lalblb;

			integral = term1 + term2;
		}

		return integral;

	} catch(exception& e){
		throw ModelException(e,method);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
/** integrate linear exponentials 
	i.e. functions of form a(t)b(t)exp[c(t)]
	where a(t), b(t) and c(t) are linear	*/
/////////////////////////////////////////////////////////////////////////////////////
double Potential::integrate_AtBtexpCt(double t1, double t2, 
									  double a1, double a2, 
									  double b1, double b2,
									  double c1, double c2)
{
	static const string method = "Potential::integrate_AtBtexpCt";
	try{
		if(t2 == t1) return 0.0;

		// First caluclate linear coefficients for a(t), b(t), c(t) from given endpoints
		// c(t) = constC + linC*t
		double linA = (a2 - a1)/(t2-t1);
		double linB = (b2 - b1)/(t2-t1);
		double linC = (c2 - c1)/(t2-t1);
		
		double constA = a1 - linA*t1;
		double constB = b1 - linB*t1;
		double constC = c1 - linC*t1;
	
		// Now write a(t)b(t) as quadratic p2*t^2 + p1*t+ p0
		double p2 = linA*linB;
		double p1 = linA*constB + linB*constA;
		double p0 = constA*constB;



		// Want integral of (constA + linA*t)*exp(constB + linB*t) --------------------
		double integral = 0.0;
		
		// when linC = 0 then are just integrating quadratic
		if(linC == 0.0) 
		{
			integral = p2/3*(t2*t2*t2-t1*t1*t1) + p1/2*(t2*t2-t1*t1) + p0*(t2-t1);
			integral *= exp(constC);
		}
		else
		{
			// split integral into linear and quadratic part:
			// integral = \int p2t^2*exp(linC*t + constC) + \int (p1*t+p0)*exp(c(t))

			// first do linear part
			integral = integrate_AtexpBt(t1, t2, (p1*t1+p0), (p1*t2 + p0), c1, c2);

			// quadratic part given by sum of three terms
			
			// term1 = p2/linC(t2^2exp(c(t2)) - t1^2exp(c(t1)))
			// term2 = -2p2/linC^2(t2*exp(c(t2)) - t1*exp(c(t1)))
			// term2 = 2p2/linC^3(exp(c(t2)) - exp(c(t1)))

			double expC2 = exp(c2);
			double expC1 = exp(c1);

			integral += p2/linC*(t2*t2*expC2 - t1*t1*expC1);
			integral -= 2*p2/(linC*linC)*(t2*expC2 - t1*expC1);
			integral += 2*p2/(linC*linC*linC)*(expC2 - expC1);

		}
		

		return integral;

	} catch(exception& e){
		throw ModelException(e,method);
	}
}
//////////////////////////////////////////////////////////////////////////////
/**  density 

  l are the lagrange multipliers
  the density phi(x) is defined at grid points and given by:
  phi(x) = prior[x]*exp(sum_i l[i]*c[i]-h[i](t))/norm
  norm = Potential(l)

 */
///////////////////////////////////////////////////////////////////////////
DoubleArraySP Potential::PDF(const DoubleArray & lag) const
{
	static const string method = "Potential::PDF";
	try{
		int i,j; //counters
		if(lag.size() != numConstraints) throw ("Number of lagrange multipliers not equal to number of constraints");
		
		DoubleArraySP phi(new DoubleArray(gridSize));
	
		/** normalisation factor */
		double normFac = norm(lag);

		double exponent;

		for(i = 0; i<gridSize;i++)
		{
			exponent = 0.0;
			for(j = 0;j<numConstraints;j++)	exponent += lag[j]*(c[j]-h[j][i]);
			(*phi)[i] = prior[i]*exp(exponent)/normFac;
		}
		
		
		return phi;


	} catch(exception& e){
		throw ModelException(e,method);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
/** cummulative density function */
////////////////////////////////////////////////////////////////////////////////////////////
DoubleArraySP Potential::CDF(const DoubleArray & lag) const
{
	static const string method = "Potential::CDF";
	try{
		
		if(lag.size() != numConstraints) throw ("Number of lagrange multipliers not equal to number of constraints");
		
		DoubleArraySP cummDist(new DoubleArray(gridSize));
	
		// integrate PDF with unit test function
		DoubleArray testFunction(gridSize, 1.0);

		integral_PhiXfX(lag, testFunction, norm(lag), *cummDist);
		
		
		return cummDist;
		
	
	} catch(exception& e){
		throw ModelException(e,method);
	}

}

////////////////////////////////////////////////////////////////////////
/** checks constraint fit 
	last check element is normalisation constraint 
	can be used for testing purposes */
////////////////////////////////////////////////////////////////////////
DoubleArraySP Potential::check(const DoubleArray & lag, double norm) const
{
	static const string method = "Potential::check";
	try{
		int j; //counters
		
		DoubleArraySP diff(new DoubleArray(c.size()+1));
		
		DoubleArray integral(gridSize);

		for(j = 0;j<c.size();j++)
		{
			integral_PhiXfX(lag, h[j], norm, integral); 
			
			//(*diff)[j] = (integral.back() - c[j])/c[j];	// relative error
			(*diff)[j] = (integral.back() - c[j]);			// absolute error	
		}


		// check normalisation
		(*diff)[c.size()] = CDF(lag)->back()-1;
		
		return diff;
		

	} catch(exception& e){
		throw ModelException(e,method);
	}

}
DRLIB_END_NAMESPACE
