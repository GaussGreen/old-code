#include	"stdafx.h"


#include "MCInt\INTUtilities.h"



// Fractional parts {x}_i := x_i - floor(x_i) = x_i mod 1
STLDoubleVector& modto01k(STLDoubleVector& x)
{
  int k = x.size();
  for (int i = 0; i < k; i++)
      modto01(x[i]);
  return x;
}

//----------------------------------------------------------------------------------------------

// String operations

// Int to string conversion
string i2s(int i)
{
  char buffer[512];
  //  _itoa(i, buffer, 10); // does not exists in Solaris
  sprintf( buffer, "%d",i);
  return string(buffer);
}

// Double to string conversion
string d2s(double d)
{
  char buffer[64];
  //  _gcvt(d, digits, buffer);
  sprintf( buffer, "%.7f",d);
  return string(buffer);
}

// Long to string conversion
string l2s(long l)
{
  char buffer[64];
  //  _ltoa(l, buffer, 10);
  sprintf( buffer, "%.7f",l);
  return string(buffer);
}



// String concatenation
string operator+(const string& s1, const string& s2)
{
  string s1s2 = s1;
  s1s2 += s2;
  return s1s2;
}

// Const char * to string concatenation
string operator+(const string& s1, const char* s2)
{
  return s1 + string(s2);  
}

// String to const char * concatenation
string operator+(const char* s1, const string& s2)
{
  return string(s1) + s2;  
}

// Int to string concatenation
string operator+(const string& s, int i)
{
  return s + i2s(i);
}

// Long to string concatenation
string operator+(const string& s, long l)
{
  return s + l2s(l);
}

// Double to string concatenation
string operator+(const string& s, double d)
{
  return s + d2s(d);
}

//----------------------------------------------------------------------------------------------

// Algebra

// Highest common factor (Euclid's algorithm)
int hcf(int a, int b, int q, int r, int s, int t, int u)
{
  if (b == 0)
      return a;
  else
      return hcf(b, a%b, s, t, q - u*s, r - u*t, a/b);
}

// Primality test
bool isprime(int p)
{
  int div = 1;
  int div_max = sqrt(p);
        
  while(++div <= div_max)
      if (hcf(div, p) > 1) return false;

  return true;
}

// Determinant of a non-singular square matrix M
double det(const STLDoubleVectorVector& M)
{
  STLDoubleVectorVector A(M);
  int k = M.size();
  STLIntegerVector indx(k);
  double det;

  ludcmp(A, k, indx, &det);

  for (int i = 0; i < k; i++) det *= A[i][i];

  return det;
}

//----------------------------------------------------------------------------------------------

// Statistics

// Vector of IID uniform deviates U_i ~ U[0,1), 1 <= i <= k
STLDoubleVector& randomU01k(STLDoubleVector& U)
{
  int k = U.size();
  for (int i = 0; i < k; i++) U[i] = randomU01();

  return U;
}

// Vector of IID uniform deviates U_i ~ U[a, b), 1 <= i <= k, a < b
STLDoubleVector& randomUabk(STLDoubleVector& U, double a, double b)
{
  randomU01k(U);
  U *= b - a;
  U += a;

  return U;
}

//----------------------------------------------------------------------------------------------

