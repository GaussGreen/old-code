

#include	"stdafx.h"

#include "MCInt\INTMiscIntegrands.h"



//----------------------------------------------------------------------------------------------

double CZaremba1integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CZaremba2integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CZaremba3integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

CHaselgroveWorstCaseintegrand::CHaselgroveWorstCaseintegrand(int s) : _s(s) 
{ 

  if ((s != 2) && (s != 4)) 
    {
      ERROR("In CHaselgroveWorstCaseintegrand::CHaselgroveWorstCaseintegrand: s must be in {2, 4}. ") 
    }

  _S = string("f-hwc(") + _s + ")"; 
}

double CHaselgroveWorstCaseintegrand::evaluate(const STLDoubleVector& x) const
{
  double result = 1.0;
  int k = x.size();
  int i;
  double factor;

  switch(_s)
    {
    case 2:
      for (i = 0; i < k; i++)
	        {
	          factor = 1.0 - fabs(x[i]);
	          factor *= factor;
	          result *= factor;
	        }
      break;
    case 4:
      for (i = 0; i < k; i++)
	{
	  factor = 1.0 - fabs(x[i]);
	  factor *= factor;
	  result *= (2.0 - factor)*factor;
	}
      break;
    }

  return result;
}

double CHaselgroveintegrand::evaluate(const STLDoubleVector& x) const
{
	int k = x.size();
	double xprod = 1.0;

	for (int i = 0; i < k; i++)
	{
		xprod *= x[i];
	}

	return exp(-xprod);
}

//----------------------------------------------------------------------------------------------

double CTsudaintegrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

double CRoosArnold1integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CRoosArnold2integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CRoosArnold3integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

double CGenz1integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CGenz2integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CGenz3integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

double CGenz4integrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

double CSchmitzbergerintegrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

double CCapstickKeisterintegrand::evaluate(const STLDoubleVector& x) const
{
	int k = x.size();
	double sum = 0.0;
	double y;

	for (int i=0; i<k; i++)
	{
		y = Nicdf((x[i] == 0.0) ? 0.00001 : x[i]);
		sum += y*y;
	}
	sum /= 2.0;

	return sqrt(pow(PI, k))*cos(sqrt(sum));
//	return sqrt(pow<double>(PI, k))*cos(sqrt(sum));
}

// Note: analytic integral over [0, 1)^k only available for k = 1, 25.
double CCapstickKeisterintegrand::I(int k) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

double CHellekalekintegrand::evaluate(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

