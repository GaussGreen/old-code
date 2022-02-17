#include	"stdafx.h"


#include "MCInt\INTEnumerator.h"


//----------------------------------------------------------------------------------------------

void CEnumerator::init(int k, const STLIntegerVector& n, const STLIntegerVector& margin)
{
  _k = k;

  if (n.size() != k) ERROR("In CEnumerator::CEnumerator: n has wrong size. ")
  _n = n;

  if (margin.size() == 0) _margin.resize(k, 0);
  else if (margin.size() != k) ERROR("In CEnumerator::CEnumerator: margin has wrong size. ")
  else _margin = margin;

  _x.resize(k);

  reset();
}

//----------------------------------------------------------------------------------------------

bool CEnumerator::operator++(void)
{
  do
    {
      ++_x[_i];

      if (_x[_i] == _n[_i] - _margin[_i])
        {
          _x[_i] = _margin[_i] - 1;
          _i--;
        }
      else
        {
          if (_i == _k - 1) return true;
          else _i++;
        }
    }
  while (_i >= 0);

  return false;
}

//----------------------------------------------------------------------------------------------

void CEnumerator::reset(void)
{
  for (int i = 0; i < _k; i++) _x[i] = _margin[i] - 1;
  _i = 0;

  if (!++(*this)) ERROR("In CEnumerator::reset: No configuration satisfies requirements. ");
}

//----------------------------------------------------------------------------------------------

STLDoubleVector& CEnumerator::xprop(STLDoubleVector& x)
{
  x.resize(_k);
  
  for (int i = 0; i < _k; i++)
    {
      x[i] = (double)_x[i]/(double)(_n[_i] - 1);
    }

  return x;
}

//----------------------------------------------------------------------------------------------

void CEnumerator::report(void) const
{
  std::cout << "Enumerator for cubic lattice in " << _k << " dimensions, " << std::endl;
  std::cout << "of size = " << _n << " with margins " << _margin << ";" << std::endl;
  std::cout << "current position x = " << _x;
}

//----------------------------------------------------------------------------------------------

bool CPermutator::ispermutation(int n, const STLIntegerVector& pi)
{
  if (pi.size() != n) ERROR("In CPermutator::ispermutation: vector pi has wrong size. ")
  
  STLBoolVector visited(n, false);
  for (int x = 0; x < n; x++)
    {
      // pi ist not injective
      if (visited[pi[x]]) return false;
        
      visited[pi[x]] = true;
    }
  
  return true;
}

//----------------------------------------------------------------------------------------------

void CPermutator::init(int k, const STLIntegerVector& n, const STLIntegerVector& margin)
{
//  _n_max = max(n); // max of a vector
	_n_max = * (max_element(n.begin(), n.end() ) ); // using proper stl
  _used.resize(_n_max, false);

  CEnumerator::init(k, n, margin);
}

//----------------------------------------------------------------------------------------------

bool CPermutator::operator++(void)
{
  do
    {
      if(_x[_i] != _margin[_i] - 1) _used[_x[_i]] = false;
      do
        {
          ++_x[_i];
          if (_x[_i] == _n[_i] - _margin[_i]) break;
        } 
      while (_used[_x[_i]]);

      if (_x[_i] == _n[_i] - _margin[_i])
        {
          _x[_i]= _margin[_i] - 1;
          _i--;
        }
      else
        {
          _used[_x[_i]] = true;
          if (_i == _k - 1) return true;
          else _i++;
        }
    }
  while (_i >= 0);

  return false;
}

//----------------------------------------------------------------------------------------------

void CPermutator::reset(void)
{
  for (int i = 0; i < _k; i++) _x[i] = _margin[i] - 1;
  for (int x = 0; x < _n_max; x++) _used[x] = false;
  _i = 0;

  if (!++(*this)) ERROR("In CPermutator::reset: No configuration satisfies requirements. ");
}

//----------------------------------------------------------------------------------------------

void CPermutator::report(void) const
{
  std::cout << "Permutator for cubic lattice in " << _k << " dimensions, " << std::endl;
  std::cout << "of size = " << _n << " with margins " << _margin << ";" << std::endl;
  std::cout << "current position x = " << _x;
}

//----------------------------------------------------------------------------------------------


