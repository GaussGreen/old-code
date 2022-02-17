#include	"stdafx.h"



#include "MCInt\INTEnumerator.h"
#include "MCInt\INTSequence.h"

#include "MCInt\INTTestUtilities.h"



//----------------------------------------------------------------------------------------------

void testutilities(void)
{
  //int i;

/*  cout << "Arithmetics: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  cout << "Comparison: " << std::endl;
  cout << "min(3, 5) = min(5, 3) = " << CC_MIN(3, 5) << " = " << CC_MIN(5, 3) << std::endl;
  cout << "max(3, 5) = max(5, 3) = " << CC_MAX(3, 5) << " = " << CC_MAX(5, 3) << std::endl;
  cout << std::endl;

  // ----------

  cout << "Exponentiation: " << std::endl;
  cout << "pow(2, 16) = " << pow(2, 16) << std::endl;
  cout << "pow(pi, 2) = " << pow(PI, 2) << std::endl;
  cout << std::endl;

  // ----------

  cout << "Absolute value: " << std::endl;
  cout << "abs(-32767) = abs(32767) = " << abs(-32767) << " = " << abs(32767) << std::endl;
  cout << "abs(-2147483647) = abs(2147483647) = " << abs(-2147483647L) << " = " << abs(2147483647L) << std::endl;
  cout << "fabs(-pi) = fabs(pi) = " << fabs(-PI) << " = " << fabs(PI) << std::endl;
  cout << std::endl;

  // ----------

  cout << "1/x with overflow protection: " << std::endl;
  cout << "1/2 = " << oneby(2) << ", 1/(-2) = " << oneby(-2) << std::endl;
  cout << "1/0 = " << oneby(0) << ", 1/(-0.0000005) = " << oneby(-0.0000005) << std::endl;
  cout << std::endl;

  // ----------

  cout << "Fractional part {x} = x - floor(x) = x mod 1: " << std::endl;
  double x;
  x = 1.0/2.0; cout << "{1/2} = " << modto01(x);
  x = 13.0/5.0; cout << ", {13/5} = " << modto01(x) << std::endl;
  x = -1.0/2.0; cout << "{-1/2} = " << modto01(x);
  x = -13.0/5.0; cout << ", {-13/5} = " << modto01(x) << std::endl;
  cout << std::endl;
  
  // ----------

  cout << "Vector and matrix manipulation: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  STLIntegerVector X;
  CreateVectorClass<int>::createvector(X, 5, 1,2,3,4,5);
  cout << "Create and print an int vector X = (x_1, ..., x_k) " << std::endl;
  cout << "X = " << X << std::endl;
  cout << std::endl;

  // ----------

  STLDoubleVector Y;
  CreateVectorClass<double>::createvector(Y, 3, PI, 2.0*PI, 3.0*PI);
  cout << "Create and print a double vector Y = (y_1, ..., y_k) " << std::endl;
  cout << "Y = " << Y << std::endl;
  cout << std::endl;

  // ----------

  STLIntegerVectorVector M;
  CreateVectorClass<int>::creatematrix(M, 4, 1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16);
  cout << "Create and print an int matrix M = ((m_11, ..., m_1k), ..., (m_k1, ..., m_kk)) " << std::endl;
  cout << "M = " << M << std::endl;
  cout << std::endl;

  // ----------

  STLDoubleVectorVector N;
  CreateVectorClass<double>::creatematrix(N, 3, 1.0,0.1,0.2, 0.1,1.0,0.1, 0.2,0.1,1.0);
  cout << "Create and print an double matrix N = ((n_11, ..., n_1k), ..., (n_k1, ..., n_kk)) " << std::endl;
  cout << "N = " << N << std::endl;
  cout << std::endl;

  cout << "Vector operations: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  CreateVectorClass<double>::createvector(Y, 4, exp(1.0),2.0*exp(1.0),3.0*exp(1.0),4.0*exp(1.0));
  cout << "Y + Y = " << Y << " + " << Y << std::endl;
  cout << "      = " << (Y += Y) << std::endl;
  cout << std::endl;

  // ----------

  int lambda = 2;
  CreateVectorClass<int>::createvector(X, 5, 1,2,3,4,5);
  cout << "X * lambda = " << X << " * " << lambda << std::endl;
  cout << "           = " << (X *= lambda) << std::endl;
  cout << std::endl;

  // ----------

  double kappa = 0.5;
  CreateVectorClass<double>::createvector(Y, 3, PI, 2.0*PI, 3.0*PI);
  cout << "Y + kappa = " << Y << " + " << kappa << std::endl;
  cout << "          = " << (Y += kappa) << std::endl;
  cout << std::endl;
  
  // ----------

  cout << "{Y + kappa} = " << modto01k(Y) << std::endl;
  cout << std::endl;

  cout << "String operations: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  string astr = "a", bstr = "b";
  cout << "String '" << astr.c_str() << "' + string '" << bstr.c_str() << "' = '" << (astr + bstr).c_str() << "'" << std::endl;
  const char *b = "b";
  cout << "String '" << astr.c_str() << "' + const char* '" << b << "' = '" << (astr + b).c_str() << "'" << std::endl;
  const char *a = "a";
  cout << "Const char* '" << a << "' + string '" << bstr.c_str() << "' = '" << (a + bstr).c_str() << "'" << std::endl;
  i = 32767;
  cout << "String '" << astr.c_str() << "' + int '" << i << "' = '" << (astr + i).c_str() << "'" << std::endl;
  long l = 2147483647L;
  cout << "String '" << astr.c_str() << "' + long '" << l << "' = '" << (astr + l).c_str() << "'" << std::endl;
  double d = PI;
  cout << "String '" << astr.c_str() << "' + double '" << d << "' = '" << (astr + d).c_str() << "'" << std::endl;
  cout << std::endl;

  cout << "Algebra: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  cout << "The highest common factor of 54 and 27 is " << hcf(54, 27) << "." << std::endl;
  cout << "The highest common factor of 48 and 64 is " << hcf(48, 64) << "." << std::endl;
  cout << std::endl;

  // ----------

  cout << "137 is " << (isprime(137) ? "" : "not ") << "prime. " << std::endl;
  cout << "138 is " << (isprime(138) ? "" : "not ") << "prime. " << std::endl;
  cout << std::endl;

  // ----------

  cout << "Cholesky decomposition of a symmetric positive definite matrix N = L*L^t" << std::endl;
  STLDoubleVectorVector L;
  cholesky(L, N);
  cout << "N = " << N << std::endl; 
  cout << "L = " << L << std::endl;
  cout << std::endl;

  N.clear();
  N.resize(10, STLDoubleVector(10, 0.0));
  for (int r = 0; r < 10; r++) N[r][r] = 1.0;
  cout << N << std::endl;
  cholesky(L, N);
  cout << "N = " << N << std::endl; 
  cout << "L = " << L << std::endl;
  cout << std::endl;

  // ----------

  STLDoubleVectorVector A;
  CreateVectorClass<double>::creatematrix(A, 4, 7.0,2.0,1.0,5.0, 3.0,1.0,2.0,1.0, 1.0,6.0,1.0,1.0, 2.0,1.0,2.0,6.0);
  cout << "Determinant of a non-singular square matrix A, " << std::endl;
  cout << "A = " << A << std::endl;
  cout << "det(A) = " << det(A) << std::endl;
  cout << std::endl;

  cout << "Combinatorics: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  cout << "Kronecker delta function: " << std::endl;
  cout << "delta(1, -1) = " << delta(1, -1) << ", delta(1, 1) = " << delta(1, 1) << std::endl;
  cout << std::endl;

  // ----------

  cout << "Factorial n! = n*(n-1)*...*1, n in |N: " << std::endl;
  cout << "5! = " << factorial(5) << ", 33! = " << factorial(33) << std::endl;
  cout << std::endl;

  // ----------

  cout << "Binomial coefficient (n choose k) = n! / (k! * (n - k)!), n >= k in |N: " << std::endl;
  cout << "(6 choose 4) = " << binomial(6, 4) << std::endl;
  cout << std::endl;

  // ----------

  cout << "(n(n-2)...)/((n-1)(n-3)...): " << std::endl;
  cout << "n = 5 odd: " << specialbinomial(5) << std::endl;
  cout << "n = 6 even: " << specialbinomial(6) << std::endl;
  cout << std::endl;

  cout << "Statistics: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  cout << "Uniform deviates U ~ U[0,1): " << std::endl;
  for (i = 0; i < 10; i++) 
    {
      if (i > 0) cout << ", ";
      cout << randomU01();
    }
  cout << std::endl << std::endl;

  // ----------

  cout << "Vectors of IID uniform deviates U_i ~ U[0,1), 1 <= i <= k: " << std::endl;
  STLDoubleVector U(5);
  for (i = 0; i < 10; i++) 
    {
      if (i > 0) cout << ", " << std::endl;
      cout << randomU01k(U);
    }
  cout << std::endl << std::endl;

  // ----------

  cout << "Normal cumulative density function: " << std::endl;
  cout << "Ncdf(-10.0) = " << Ncdf(-10.0) << ", Ncdf(10.0) = " << Ncdf(10.0) << ", " << std::endl;
  cout << "Ncdf(0) = " << Ncdf(0.0) << ", ";
  cout << "Ncdf(-0.5) = " << Ncdf(-0.5) << ", Ncdf(0.5) = " << Ncdf(0.5) << std::endl;
  cout << std::endl;

  // ----------

  cout << "Normal inverse cumulative density function: " << std::endl;
  cout << "Nicdf(0.01) = " << Nicdf(0.01) << ", Nicdf(0.99) = " << Nicdf(0.99) << ", " << std::endl;
  cout << "Nicdf(0.5) = " << Nicdf(0.5) << ", ";
  cout << "Nicdf(0.4) = " << Nicdf(0.4) << ", Nicdf(0.6) = " << Nicdf(0.6) << std::endl;
  cout << "Nicdf(Ncdf(0.3)) = " << Nicdf(Ncdf(0.3)) << ", Ncdf(Nicdf(0.7)) = " << Ncdf(Nicdf(0.7)) << std::endl;
  cout << std::endl;

  cout << "Sequences: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  bool verbose = true;

  cout << "Sequence of rationals (r_i) in |Q intersect [0,1): " << std::endl;
  CRationalsequence R;
  
  for (i = 1; i < 4; i++) R.computenext(verbose);
  cout << std::endl; 
  
  Y.resize(10);
  cout << R.get(Y) << std::endl;
  cout << std::endl; 

  // ----------

  cout << "Sequence of primes (p_i): " << std::endl;
  CPrimesequence P;

  for (i = 1; i < 6; i++) P.computenext(verbose);
  cout << std::endl;
  
  Y.resize(10);
  cout << P.get(Y) << std::endl;
  cout << std::endl; 

  cout << "Lattice: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  int k = 3;
  STLIntegerVector res, margin;
  CreateVectorClass<int>::createvector(res, k, 5, 5, 6);
  CreateVectorClass<int>::createvector(margin, k, 2, 1, 1);
  CEnumerator enumerator(k, res, margin);

  enumerator.report(); cout << std::endl;
  cout << std::endl;

  i = 0;
  do { cout << i << ". " << enumerator._x << std::endl; i++; } while(++enumerator);
  cout << std::endl;

  cout << "Permutator: " << std::endl;
  cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  CPermutator permutator(k, res, margin);

  permutator.report(); cout << std::endl;
  cout << std::endl;

  i = 0;
  do { cout << i << ". " << permutator._x << std::endl; i++; } while(++permutator);
  cout << std::endl;
*/
}

//----------------------------------------------------------------------------------------------

