#include	"stdafx.h"



#include "INTIntegrator.h"



//----------------------------------------------------------------------------------------------

// Compute integral of f over domain D using N function samples
double CMCintegrator::integrate(const CIntegrand& f, const CDomain& D, long N)
{
  double integral = 0.0;
  STLDoubleVector x(D._k);

  for (int m = 1; m <= N; m++)
    {
      randomU01k(x);
      integral += f.evaluate(x, D);
    }
 
   return integral/(double)N;
}

//----------------------------------------------------------------------------------------------

// Compute integral of f over domain D using N function samples and specified set of parameters
double CNTintegrator::integrate(const CIntegrand& f, const CDomain& D, long N, const CNTparameters& par)
{
  CAlpha* a = ((par._a != NULL) ? par._a : _defaultpar._a);
  CWeights* w = ((par._w != NULL) ? par._w : _defaultpar._w);
  CPeriodization* p = ((par._p != NULL) ? par._p : _defaultpar._p);

  if (a == NULL) ERROR("In CNTintegrator::integrate: No alpha available. ")
  if (w == NULL) ERROR("In CNTintegrator::integrate: No weights available. ")
  if (p == NULL) ERROR("In CNTintegrator::integrate: No periodization available. ")

  STLDoubleVector alpha;
  a->get(alpha, D._k);

  long Neff = w->Neff(N);
  EQSP_CC_STL_PAIR(long, long) Nrange;
  w->Nrange(Nrange, Neff);
  double C = w->C(Neff);
  double xshift = w->xshift();

  double integral = 0.0;
  double Am;
  STLDoubleVector x;

  for (int m = Nrange.first; m <= Nrange.second; m++)
    {
      Am = w->A(Neff, m);

      x = alpha;
      x *= (double)m;
      x += xshift;

      integral += Am * f.evaluate(x, D, p);
    }

  return C*integral;
}




void CNTintegrator::integrate(vector<double>& integral,const vector<CIntegrand>& f, const CDomain& D, long N, const CNTparameters& par)
{
  CAlpha* a = ((par._a != NULL) ? par._a : _defaultpar._a);
  CWeights* w = ((par._w != NULL) ? par._w : _defaultpar._w);
  CPeriodization* p = ((par._p != NULL) ? par._p : _defaultpar._p);

  if (a == NULL) ERROR("In CNTintegrator::integrate: No alpha available. ")
  if (w == NULL) ERROR("In CNTintegrator::integrate: No weights available. ")
  if (p == NULL) ERROR("In CNTintegrator::integrate: No periodization available. ")

  vector<double> alpha;
  a->get(alpha, D._k);

  long Neff = w->Neff(N);
  pair<long, long> Nrange;
  w->Nrange(Nrange, Neff);
  double C = w->C(Neff);
  double xshift = w->xshift();


  double Am;
  vector<double> x;

  int nfunctions = f.size();

  for (int m = Nrange.first; m <= Nrange.second; m++)
    {
      Am = w->A(Neff, m);

      x = alpha;
      x *= (double)m;
      x += xshift;

	  for (int i = 0 ; i < nfunctions; i++ )
	  {
	      integral[i] += Am * C*f[i].evaluate(x, D, p);
	  }
    }
}



void CNTintegrator::getIntegrationPoint( int m, vector<double>& x, const CDomain& D,  const CNTparameters& par)
{
  CAlpha* a = ((par._a != NULL) ? par._a : _defaultpar._a);
  CWeights* w = ((par._w != NULL) ? par._w : _defaultpar._w);

  if (a == NULL) ERROR("In CNTintegrator::integrate: No alpha available. ")
  if (w == NULL) ERROR("In CNTintegrator::integrate: No weights available. ")

  vector<double> alpha;
  a->get(alpha, D._k);

  double xshift = w->xshift();

   x = alpha;
   x *= (double)m;
   x += xshift;

   modto01k(x);	 
}



void CNTintegrator::getSlicedIntegrationPoints(vector<double>& x,int islice, const CDomain& D, long N, const CNTparameters& par)
{

// size of vector is npaths

  CAlpha* a = ((par._a != NULL) ? par._a : _defaultpar._a);
  CWeights* w = ((par._w != NULL) ? par._w : _defaultpar._w);
  CPeriodization* p = ((par._p != NULL) ? par._p : _defaultpar._p);

  if (a == NULL) ERROR("In CNTintegrator::integrate: No alpha available. ")
  if (w == NULL) ERROR("In CNTintegrator::integrate: No weights available. ")
  if (p == NULL) ERROR("In CNTintegrator::integrate: No periodization available. ")

  vector<double> alpha;
  a->get(alpha, D._k);

  long Neff = w->Neff(N);
  pair<long, long> Nrange;
  w->Nrange(Nrange, Neff);
  double C = w->C(Neff);
  double xshift = w->xshift();

  double integral = 0.0;
  double Am;

  int  ndim = Nrange.second-Nrange.first+1;

  x.resize(ndim); 

  int k;

  for (int m = Nrange.first; m <= Nrange.second; m++)
    {
      Am = w->A(Neff, m);

	  k	 = m-Nrange.first;

      x[k] = alpha[islice];

      x[k] *= (double)m;
      x[k] += xshift;
	  modto01(x[k]);	
    }

 
}

void CNTintegrator::getSlicedIntegrationPoints(CVector& x,int islice, const CDomain& D, long N, const CNTparameters& par)
{

// size of vector is npaths

  CAlpha* a = ((par._a != NULL) ? par._a : _defaultpar._a);
  CWeights* w = ((par._w != NULL) ? par._w : _defaultpar._w);
  CPeriodization* p = ((par._p != NULL) ? par._p : _defaultpar._p);

  if (a == NULL) ERROR("In CNTintegrator::integrate: No alpha available. ")
  if (w == NULL) ERROR("In CNTintegrator::integrate: No weights available. ")
  if (p == NULL) ERROR("In CNTintegrator::integrate: No periodization available. ")

  vector<double> alpha;
  a->get(alpha, D._k);

  long Neff = w->Neff(N);
  pair<long, long> Nrange;
  w->Nrange(Nrange, Neff);
  double C = w->C(Neff);
  double xshift = w->xshift();

  double integral = 0.0;
  double Am;

  int  ndim = Nrange.second-Nrange.first+1;

//  x.resize(ndim); 

  int k;

  for (int m = Nrange.first; m <= Nrange.second; m++)
    {
      Am = w->A(Neff, m);

	  k	 = m-Nrange.first;

      x[k] = alpha[islice];

      x[k] *= (double)m;
      x[k] += xshift;
	  modto01(x[k]);	
    }

 
}















// ALEXL
class CMyintegrand: public CIntegrand
{

	public:

	inline double evaluate(const vector<double>& x)const
	{

		double value;
		value = x[0]+x[1]+x[2]+x[3];
		return value;
	}

};

double stupidAlextest()
{

	long idum = -1234;

	ran2init(idum);

	CMyintegrand f;

	int ndim = 4;
	CUnitcube D(4);// CGeneralCube,Cube see testintegrands

	CNTintegrator NTintegrator;

//	select alpha

	CNiederreiteralpha alpha;

	CArithmeticmeanweights weights;// weight

	int d = 3;// d >= r 
	CPolynomialperiodization periodization(d);

	CNTparameters par(&alpha,&weights,&periodization);

	int npoints = 1000;

	double res;

	res = NTintegrator.integrate(f,D,npoints,par);

	return res;
};












