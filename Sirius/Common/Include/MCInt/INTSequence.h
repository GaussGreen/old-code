
#ifndef __INTSEQUENCE_H__
#define __INTSEQUENCE_H__

#include "INTSTL.h"
#include "INTUtilities.h"

//----------------------------------------------------------------------------------------------

// Base class for sequences
class  CSequence
{
 protected:
  // Buffered sequence elements x
  STLDoubleVector _x;

  // Current size n of buffer
  int _n;

  CSequence() : _n(0) {}

 public:
  // Retrieve i-th sequence element x_i, i >= 1
  double get(int i);

  // Retrieve vector (x_1, ..., x_k) of k leading sequence elements
  STLDoubleVector& get(STLDoubleVector& x);

  // Compute and buffer next sequence element x_{n+1}
  virtual void computenext(bool verbose = false) = 0;
};

//----------------------------------------------------------------------------------------------

// Sequence of rationals (r_i) in |Q intersect [0,1) (for Baker alpha: alpha_i = exp(m*r_i))
class  CRationalsequence : public CSequence
{
 private:
  // Current distance from origin
  int _d;

  // Current numerator and denominator
  int _num;
  int _den;

 public:
   CRationalsequence(void) : _d(2), _num(1), _den(1) {}

  // Compute and buffer next rational r_{n+1}
  void computenext(bool verbose = false);
};

//----------------------------------------------------------------------------------------------

// Sequence of primes (p_i) (for prime root alpha: alpha_i = sqrt(p_i))
class  CPrimesequence : public CSequence
{
 public:
  // Compute and buffer next prime p_{n+1}
  void computenext(bool verbose = false);
};

//----------------------------------------------------------------------------------------------


#endif
