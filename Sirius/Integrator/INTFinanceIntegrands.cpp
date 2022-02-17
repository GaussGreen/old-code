#include "stdafx.h"


#include ".\mcint\INTFinanceIntegrands.h"

static double xMax(double x,double y)
{
	if ( x > y )
		return x;
	else
		return y;
}
	//----------------------------------------------------------------------------------------------

/*CEuropeanintegrand::CEuropeanintegrand(const STLDoubleVector& t, int type, double S0, 
                                       double K, double rho, double T, double vol) :
  _t(t), _type(type), _S0(S0), _K(K), _rho(rho), _T(T), _vol(vol)
{
  // Statistics from sequence of intermediate times (t_i), required for evaluation of payoff
  int k = t.size() - 1;
  _c.resize(k);

  for (int i=0; i<k; i++)
    {
      _c[i] = sqrt(t[i+1] - t[i]);
    }

  _S = "f-eo";
}
*/
CEuropeanintegrand::CEuropeanintegrand(const STLDoubleVector& t, int type, double S0, 
                                       double K, double rho, double T, double vol) 
{}
//----------------------------------------------------------------------------------------------

double CEuropeanintegrand::evaluate(const STLDoubleVector& x) const
{
	int k = x.size();
	
	double sum = 0.0;
	for (int i = 0; i < k; i++)
    {
		sum += _c[i]*Nicdf(x[i]);
    }
	
	double ST = _S0*exp((_rho - 0.5*_vol*_vol)*_T + _vol*sum);
	
	switch(_type)
    {
	case 0:
	default:
		return xMax(ST - _K, 0.0);
		
    }
	
	return 0.0;
}

//----------------------------------------------------------------------------------------------

double CEuropeanintegrand::I(int k) const 
{ 
  return BS(_type, _S0, _K, _rho, _T, _vol);
}

//----------------------------------------------------------------------------------------------

void CEuropeanintegrand::report(void) const
{
	std::cout << "European option: " << std::endl;
	std::cout << "type = ";
	switch(_type)
    {
	case 0:
	default:
		std::cout << "call, ";
		break;
		
    }
	std::cout << "S0 = " << _S0 << ", K = " << _K << ", rho = " << _rho << ", ";
	std::cout << "T = " << _T << ", vol = " << _vol << ", " << std::endl;
	std::cout << "with intermediate times t = " << _t;
}

//----------------------------------------------------------------------------------------------

// Asian option using geometric average of the spot price at control times 0 = t_0 <= t_1 <= ... <= t_k = T
/*CAsianintegrand::CAsianintegrand(const STLDoubleVector& t, int type, double S0, 
                                 double K, double rho, double T, double vol) :
  _t(t), _type(type), _S0(S0), _K(K), _rho(rho), _T(T), _vol(vol)
{
  int k = t.size() - 1;

  _c1 = 0.0;
  _c2.resize(k);
  _c3 = 0.0;

  for (int i = 0; i < k; i++)
    {
      _c1 += t[i+1];
      _c2[i] = (k - i)*sqrt(t[i+1] - t[i]);
      _c3 += _c2[i]*_c2[i];
    }

  _c1 /= (double)k;
  _c3 = sqrt(_c3)/(double)k;

  _S = "f-ao";
}
*/


CAsianintegrand::CAsianintegrand(const STLDoubleVector& t, int type, double S0, 
                                 double K, double rho, double T, double vol) 
 {
  int k = t.size() - 1;

  _c1 = 0.0;
  _c2.resize(k);
  _c3 = 0.0;

  for (int i = 0; i < k; i++)
    {
      _c1 += t[i+1];
      _c2[i] = (k - i)*sqrt(t[i+1] - t[i]);
      _c3 += _c2[i]*_c2[i];
    }

  _c1 /= (double)k;
  _c3 = sqrt(_c3)/(double)k;

  _S = "f-ao";
}

//----------------------------------------------------------------------------------------------

double CAsianintegrand::evaluate(const STLDoubleVector& x) const
{
	int k = x.size();
	
	double sum = 0.0;
	for (int i = 0; i < k; i++)
    {
		sum += _c2[i]*Nicdf(x[i]);
    }
	
	double ST = _S0 * exp((_rho - 0.5*_vol*_vol)*_c1 + _vol/(double)k*sum);
	
	switch(_type)
    {
	case 0:
	default:
		return xMax(ST - _K, 0.0);
		
    }
	
	return 0.0;  
}

//----------------------------------------------------------------------------------------------

double CAsianintegrand::I(int k) const
{
  double volhat = _vol*_c3/sqrt(_T);
  double S0hat = _S0*exp((_rho - 0.5*_vol*_vol)*_c1 + 0.5*volhat*volhat*_T);

  return BS(_type, S0hat, _K, 0.0, _T, volhat);
}

//----------------------------------------------------------------------------------------------

void CAsianintegrand::report(void) const
{
	std::cout << "Asian option: " << std::endl;
	std::cout << "type = ";
	switch(_type)
    {
	case 0:
	default:
		std::cout << "call, ";
		
    }
	std::cout << "S0 = " << _S0 << ", K = " << _K << ", rho = " << _rho << ", ";
	std::cout << "T = " << _T << ", vol = " << _vol << ", " << std::endl;
	std::cout << "with control times t = " << _t;
}

//----------------------------------------------------------------------------------------------

CBasketintegrand::CBasketintegrand(const STLDoubleVector& S0, const STLDoubleVector& vol, const STLDoubleVectorVector& corr, 
                                   int type, double K, double rho, double T) 
//  _S0(S0), _vol(vol), _corr(corr), _type(type), _K(K), _rho(rho), _T(T) 
{
  int k = S0.size();
  if (vol.size() != k) ERROR("In CBasketintegrand::CBasketintegrand: Dimensions of vol and S0 incompatible. ")
  if (corr.size() != k) ERROR("In CBasketintegrand::CBasketintegrand: Dimensions of corr and S0 incompatible. ")

  STLDoubleVectorVector L;
  cholesky(L, corr);

  _c1 = 1.0;
  _c2 = 0.0;
  _c3.resize(k);

  for (int i = 0; i < k; i++)
    {
      _c1 *= _S0[i];
      _c2 += _vol[i]*_vol[i];

      _c3[i] = 0.0;
      for (int j = i; j < k; j++)
        {
          _c3[i] += _vol[j]*L[j][i];
        }
    }

  _c1 = exp(log(_c1)/(double)k);
  _c2 /= (double)k;

  _S = "f-bo";
}

//----------------------------------------------------------------------------------------------

double CBasketintegrand::evaluate(const STLDoubleVector& x) const        
{
	int k = x.size();
	
	double sum = 0.0;
	for (int i = 0; i < k; i++)
    {
		sum += _c3[i]*Nicdf(x[i]);
    }
	
	double ST = _c1*exp((_rho - 0.5*_c2)*_T + (sqrt(_T)/(double)k)*sum);
	
	switch(_type)
    {
	case 0:
	default:
		return xMax(ST - _K, 0.0);
    }
	
	return 0.0;  
}

//----------------------------------------------------------------------------------------------

double CBasketintegrand::I(int k) const
{
  double sum = 0.0;
  for (int i = 0; i < k; i++)
    {
      sum += _c3[i]*_c3[i];
    }

  double volhat = (1.0/(double)k)*sqrt(sum);
  double S0hat = _c1*exp((_rho - 0.5*(_c2 - volhat*volhat))*_T);

  return BS(_type, S0hat, _K, 0.0, _T, volhat);
}

//----------------------------------------------------------------------------------------------

void CBasketintegrand::report(void) const
{
	std::cout << "Basket option: " << std::endl;
	std::cout << "type = ";
	switch(_type)
    {
	case 0:
	default:
		std::cout << "call, ";
		
    }
	std::cout << "K = " << _K << ", rho = " << _rho << ", T = " << _T << ", with " << std::endl;
	std::cout << "SO = " << _S0 << "," << std::endl;
	std::cout << "vol = " << _vol << " and correlation matrix" << std::endl;
	std::cout << "C = " << _corr;
}

//----------------------------------------------------------------------------------------------

