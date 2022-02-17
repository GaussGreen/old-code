#include	"stdafx.h"

#include "MCInt\INTDomain.h"


//----------------------------------------------------------------------------------------------

CParallelepiped::CParallelepiped(int k, const STLDoubleVectorVector& M) : CLineardomain(k), _M(M) 
{ 
  _J = fabs(det(M));
  _S = string("D-parepip(") + _k + ")";
}

//----------------------------------------------------------------------------------------------

// Domain transformation for parallelepiped
STLDoubleVector& CParallelepiped::psi(STLDoubleVector& x) const
{
  if (x.size() != _k) ERROR("In CParallelepiped::psi: vector x has wrong size. ")

  STLDoubleVector y(x);
  for (int i = 0; i < _k; i++)
    {
      x[i] = 0.0;
      for (int j = 0; j < _k; j++)
        {
          x[i] += _M[i][j]*y[j];
        }
    }

  return x;
}

//----------------------------------------------------------------------------------------------

// Report parallelepiped
void CParallelepiped::report(void) const
{ 
  std::cout << "Parallelepiped in " << _k << " dimensions with transformation matrix M, " << std::endl;
  std::cout << "M = " << _M << std::endl;
};

//----------------------------------------------------------------------------------------------

CGeneralcube::CGeneralcube(int k, ...) : CRectilineardomain(k)
{
  _L.resize(k);
  int i, k2 = k*k;

  // Initialise limits L
  va_list ap;
  va_start(ap, k);
  for (i = 0; i < k; i++)
    {
      _L[i] = make_pair(va_arg(ap, double), va_arg(ap, double));
    }
  va_end(ap);
    
  // Compute constant Jacobian J
  _J = 1.0;
  for (i = 0; i < _k; i ++)
  {
	  _J *= (_L[i].second - _L[i].first);
  }

  _S = string("D-gencube(") + _k + ")"; 
};

//----------------------------------------------------------------------------------------------

// Domain transformation for general cube
STLDoubleVector& CGeneralcube::psi(STLDoubleVector& x) const
{
  if (x.size() != _k) ERROR("In CGeneralcube::psi: vector x has wrong size. ")

  for (int i = 0; i < _k; i++)
    {
      x[i] = (_L[i].second - _L[i].first)*x[i] + _L[i].first;
    }

  return x;
}

//----------------------------------------------------------------------------------------------

// Report general cube
void CGeneralcube::report(void) const
{ 
  std::cout << "General cube ";
  for (int i = 0; i < _k; i++)
    {
      if (i > 0) std::cout << " x ";
      std::cout << "[" << _L[i].first << ", " << _L[i].second << ")";
    }
  std::cout << " in " << _k << " dimensions";
};

//----------------------------------------------------------------------------------------------

CCube::CCube(int k, double a, double b) : CRectilineardomain(k)
{
  _L.first = a;
  _L.second = b;

  _J = pow(b - a, k);

  _S = string("D-[") + _L.first + ", " + _L.second + ")^" + _k; 
}

//----------------------------------------------------------------------------------------------

// Domain transformation for cube
STLDoubleVector& CCube::psi(STLDoubleVector& x) const
{
  if (x.size() != _k) ERROR("In CCube::psi: vector x has wrong size. ")

  for (int i = 0; i < _k; i++)
    {
      x[i] = (_L.second - _L.first)*x[i] + _L.first;
    }

  return x;
}

//----------------------------------------------------------------------------------------------

// Report cube
void CCube::report(void) const
{ 
  std::cout << "Cube ";
  std::cout << "[" << _L.first << ", " << _L.second << ")^" << _k; 
  std::cout << " in " << _k << " dimensions";
};

//----------------------------------------------------------------------------------------------

// Report unit cube
void CUnitcube::report(void) const
{ 
  std::cout << "Unit cube [0, 1)^" << _k;
  std::cout << " in " << _k << " dimensions";
};

//----------------------------------------------------------------------------------------------

// Domain transformation for ellipsoidal domain
STLDoubleVector& CEllipsoidaldomain::psi(STLDoubleVector& x) const
{
  // TODO
  return x;  
}

//----------------------------------------------------------------------------------------------

// Jacobian for ellipsoidal domain
double CEllipsoidaldomain::J(const STLDoubleVector& x) const
{
  // TODO
  return 0.0;
}

//----------------------------------------------------------------------------------------------

CEllipsoid::CEllipsoid(int k, STLDoubleVector& r, STLDoubleVector& C) : CEllipsoidaldomain(k, C), _r(r)
{
  _S = string("D-ellip(") + _k + ")"; 
}

//----------------------------------------------------------------------------------------------

// Report ellipsoid
void CEllipsoid::report(void) const
{ 
  // TODO
};

//----------------------------------------------------------------------------------------------

CSphere::CSphere(int k, double r, STLDoubleVector& C) : CEllipsoidaldomain(k, C), _r(r)
{
  _S = string("D-sphere(") + _r + ", " + _k + ")";   
}

//----------------------------------------------------------------------------------------------

// Report sphere
void CSphere::report(void) const
{ 
  // TODO
};

//----------------------------------------------------------------------------------------------

CUnitsphere::CUnitsphere(int k, STLDoubleVector& C) : CEllipsoidaldomain(k, C) 
{  
  _S = string("D-sphere(1, ") + _k + ")"; 
}

//----------------------------------------------------------------------------------------------

// Report unit sphere
void CUnitsphere::report(void) const
{ 
  // TODO
};

//----------------------------------------------------------------------------------------------

