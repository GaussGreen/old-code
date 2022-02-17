#include	"stdafx.h"



#include "MCInt\INTSequence.h"


//----------------------------------------------------------------------------------------------

// Retrieve i-th sequence element x_i, i >= 1
double CSequence::get(int i)
{
  if (i < 1) ERROR("In CSequence::get: i must be positive. ")

  for (int j = _n; j < i; j++)
    {
      // Compute missing x_i on demand
      computenext();
    }

  return _x[i-1];
}

//----------------------------------------------------------------------------------------------

// Retrieve vector (x_1, ..., x_k) of k leading sequence elements
STLDoubleVector& CSequence::get(STLDoubleVector& x)
{
  int k = x.size();
  int i;

  for (i = _n; i < k; i++)
    {
      // Compute missing x_i on demand
      computenext();
    }

  for (i = 0; i < k; i++)
    {
      x[i] = _x[i];
    }

  return x;
}

//----------------------------------------------------------------------------------------------

// Compute and buffer next rational r_{n+1}
void CRationalsequence::computenext(bool verbose)
{
  do
    {
      if (--_num == 0)
	      {
	        _d++;
	        _num = (_d%2 == 0) ? _d/2 - 1 : _d/2;
	        _den = _d - _num;
	      }
      else
	      {
	        _den++;
	      }
    } 
  while (hcf(_num, _den) > 1);

  _n++;
  _x.push_back((double)_num/(double)_den);

  if (verbose)
    {
      std::cout << _n << "th rational = " << _num << "/" << _den << std::endl;
    }
}

//----------------------------------------------------------------------------------------------

// Compute and buffer next prime p_{n+1}
void CPrimesequence::computenext(bool verbose)
{
  int p = (_n == 0 ? 1 : _x[_n-1]);
  while(!isprime(++p)) {};

  _n++;
  _x.push_back(p);

  if (verbose)
    {
      std::cout << _n << "th prime = " << p << std::endl;
    }
}

//----------------------------------------------------------------------------------------------

