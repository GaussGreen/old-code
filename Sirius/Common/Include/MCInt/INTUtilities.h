
#ifndef __INTUTILITIES_H__
#define __INTUTILITIES_H__

#include "INTAlgebra.h"
#include "INTCombinatorics.h"
#include "INTRng.h"
#include "INTNormal.h"

// 1/x with overflow protection
inline  double oneby(double x) { return (fabs(x) < 0.000001) ? ((x < 0) ? -1.0 : 1.0)*1000000.0 : 1/x; }

// Fractional part {x} = x - floor(x) = x mod 1 of x in |R
inline  double& modto01(double& x) { x -= floor(x); return x; }


//----------------------------------------------------------------------------------------------

// Vector and matrix manipulation


template<class T> class  CreateVectorClass
{

 public:
  // Create a vector x = (x_1, ..., x_k)
  static vector < T >& createvector(vector < T >& x, int k, ...)
    {
      x.resize(k);
      
      va_list ap;
      va_start(ap, k);
      for (int i = 0; i < k; i++) 
        x[i] = va_arg(ap, T);
      
      va_end(ap);
      
      return x;
    }
  
  // Create a pseumatrix M = ((m_11, ..., m_1k), ..., (m_k1, ..., m_kk))
  static vector < vector <T> >& creatematrix(vector <vector <T> >& M, int k, ...)
    {
      M.clear();
      M.resize(k, vector <T>(k));
      int k2 = k*k;
      
      va_list ap;
      va_start(ap, k);
      for (int i = 0; i < k2; i++)
        M[i/k][i%k] = va_arg(ap, T);
      
      va_end(ap);
      
      return M;
    }
};



// Print a vector x = (x_1, ..., x_k)
template<class T>  std::ostream& operator<<(std::ostream& os, const vector <T>& x)
{
  int k = x.size();
  os << "(";
  for (int i = 0; i < k; i++)
    {
      if (i > 0) os << ", ";
      os << x[i];
    }
  os << ")";

  return os;
}


//----------------------------------------------------------------------------------------------

// Vector operations

// Vector addition (x + y)_i := x_i + y_i
template<class T>  vector <T>& operator+=(vector <T>& x, const vector <T>& y)
{
  int k = x.size();
  if (y.size() != k) ERROR("In operator+=: vector y has wrong size. ")

  for (int i = 0; i < k; i++)
      x[i] += y[i];

  return x; 
}

// Scalar addition (x + lambda)_i := x_i + lambda
template<class T>  vector <T>& operator+=(vector <T>& x, T lambda)
{
  int k = x.size();
  for (int i = 0; i < k; i++)
      x[i] += lambda;

  return x;
}

// Scalar multiplication (x * lambda)_i := x_i * lambda
template<class T>  vector <T>& operator*=(vector <T>& x, T lambda)
{
  int k = x.size();
  for (int i = 0; i < k; i++)
      x[i] *= lambda;

  return x;
}

// Fractional parts {x}_i := x_i - floor(x_i) = x_i mod 1 of x in |R^k
 STLDoubleVector& modto01k(STLDoubleVector& x);

//----------------------------------------------------------------------------------------------

// String operations

// Int to string conversion
 string i2s(int i);

// Long to string conversion
 string l2s(long l);

// Double to string conversion
 string d2s(double d );


// String concatenation
 string operator+(const string& s1, const string& s2);

// Const char * to string concatenation
 string operator+(const string& s1, const char* s2);

// String to const char * concatenation
 string operator+(const char* s1, const string& s2);

// Int to string concatenation
 string operator+(const string& s, int i);

// Long to string concatenation
 string operator+(const string& s, long l);

// Double to string concatenation
 string operator+(const string& s, double d);

//----------------------------------------------------------------------------------------------

// Algebra

// Highest common factor
 int hcf(int a, int b, int q = 1, int r = 0, int s = 0, int t = 1, int u = 0);

// Primality test
 bool isprime(int p);

// Cholesky decomposition of a symmetric positive definite square matrix M = L*L^t
#define cholesky choldc

// Determinant of a non-singular square matrix M
 double det(const STLDoubleVectorVector& M);

//----------------------------------------------------------------------------------------------

// Combinatorics

// Kronecker delta function
inline  int delta(int i, int j) { return (i == j) ? 1 : 0; }

// Factorial n! = n*(n-1)*...*1, n in |N
inline  double factorial(int n) { return factr(n); }

// Binomial coefficient (n choose k) = n! / (k! * (n - k)!), n >= k in |N
template<class T> inline  double binomial(T n, T k) { return bico(n, k); };

// (n(n-2)...)/((n-1)(n-3)...)
template<class T>  double specialbinomial(T n)
{
  T i = 0;
  double result = 1.0;

  while (i < n)
    {
      if (i%2 == 0) result *= (double)(n-i);
      else result /= (double)(n-i);
      i++;
    }
        
  return result;
}

//----------------------------------------------------------------------------------------------

// Statistics

// Uniform deviate U ~ U[0, 1)
#define randomU01 RAN2

// Uniform deviate U ~ U[a, b), a < b
inline  double randomUab(double a, double b) { return (b - a)*randomU01() + a; }

// Vector of IID uniform deviates U_i ~ U[0, 1), 1 <= i <= k
 STLDoubleVector& randomU01k(STLDoubleVector& U);

// Vector of IID uniform deviates U_i ~ U[a, b), 1 <= i <= k, a < b
 STLDoubleVector& randomUabk(STLDoubleVector& U, double a, double b);

// Normal cumulative density function
#define Ncdf Ncdf2

// Normal inverse cumulative density function
#define Nicdf Nicdf2

//----------------------------------------------------------------------------------------------

// Plotting

// Plot table Z with row labels X, column labels Y and title t to the specified file
template<class S, class T, class U, class V> 
 void plotXYZ(const vector <S>& X, const vector <T>& Y, const vector <U>& Z, 
								  const V& t, const string& filename)
{
  int R = X.size();
  int r;

  int C = Y.size();
  int c;

  // Validate input
  if ((R <= 0) || (C <= 0))						ERROR("In plotXYZ: empty plot. ")
  if (Z.size() != R)							ERROR("In plotXYZ: table Z has wrong number of rows. ")
  if (Z[0].size() != C)							ERROR("In plotXYZ: table Z has wrong number of columns. ")
  for (r = 1; r < R; r++) if (Z[r].size() != C) ERROR("In plotXYZ: table Z not an NxM square matrix. ")

  // Open output file
  std::ofstream os(filename.c_str());

  // Plot column labels
  os << t << "\t";
  for (c = 0; c < C; c++) os << Y[c] << "\t";
  os << std::endl;

  // Plot rows, with leading row label
  for (r = 0; r < R; r++)
    {
      os << X[r] << "\t";
      for (c = 0; c < C; c++) 
          os << Z[r][c] << "\t";
      os << endl;
    }

  os.close();
}

//----------------------------------------------------------------------------------------------


#endif
