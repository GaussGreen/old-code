
#ifndef __INTNUMERATOR_H__
#define __INTNUMERATOR_H__


#include "INTSTL.h"
#include "INTUtilities.h"

//----------------------------------------------------------------------------------------------

class  CEnumerator
{
 protected:
  int _k;
  int _i;

  STLIntegerVector _n;
  STLIntegerVector _margin;

  CEnumerator() {}

  virtual void init(int k, const STLIntegerVector& n, const STLIntegerVector& margin);

 public:
  STLIntegerVector _x;

  CEnumerator(int k, const STLIntegerVector& n, const STLIntegerVector& margin) 
    { 
      init(k, n, margin); 
    }
  CEnumerator(int k, const STLIntegerVector& n, int margin = 0) 
    { 
      STLIntegerVector margink(k, margin);
      init(k, n, margink); 
    }
  CEnumerator(int k, int n, int margin = 0)
    { 
      STLIntegerVector nk(k, n);
      STLIntegerVector margink(k, margin);
      init(k, nk, margink); 
    }

  STLDoubleVector& xprop(STLDoubleVector& x);

  virtual bool operator++(void);

  virtual void reset(void);

  virtual void report(void) const;
};

//----------------------------------------------------------------------------------------------

class  CPermutator : public CEnumerator
{
 private:  
  int _n_max;
  vector<bool> _used;

  void init(int k, const STLIntegerVector& n, const STLIntegerVector& margin);

 public:
  CPermutator(int k, const STLIntegerVector& n, const STLIntegerVector& margin) 
    { 
      init(k, n, margin); 
    }
  CPermutator(int k, const STLIntegerVector& n, int margin = 0) 
    { 
      STLIntegerVector margink(k, margin);
      init(k, n, margink); 
    }
  CPermutator(int k, int n, int margin = 0)
    { 
      STLIntegerVector nk(k, n);
      STLIntegerVector margink(k, margin);
      init(k, nk, margink); 
    }
  CPermutator(int n)
    { 
      STLIntegerVector nn(n, n);
      STLIntegerVector marginn(n, 0);
      init(n, nn, marginn); 
    }

  static bool ispermutation(int n, const STLIntegerVector& x);

  bool operator++(void);

  void reset(void);

  void report(void) const;
};

//----------------------------------------------------------------------------------------------


#endif
