#ifndef __konedimsolver_h
#define __konedimsolver_h

#include "kexception.h"

/** Small double constant close to machine precision */
const double SOLVER_DOUBLE_EPS = 3e-15;

/** One dimensional solver interface. 
@see KSecantRoot */
template <class F>
class KOneDimSolver {
public:
/** Perform solving algorithm

  Attempts to find solution f(x) = target, starting with x = guess.
  @param guess Initial guess
  @param target Target value
	*/
	virtual double solve (F& f, double guess, double target) = 0;
	
	/**	Number of calls to f in last solver run. */
	int numFCalls() {return _numFCalls;}
	/**	f(x) where x (the solution) and f are from the last solver run. */
	double fValue() {return _fValue;}
protected:
	KOneDimSolver () : _numFCalls(-999999), _fValue (-999999){}
	int _numFCalls;
	double _fValue;
};

/** Secant root solver.

  Finds root using secant search (see NRC2)
  Secant search works well for near-linear functions
  For HIGHLY non-linear functions, you probably want to use Brent.
  
	F can either be a functional or a function reference.
	They both work, since they have the same syntax with respect to operator ()
	
	  Functionals have the method:
	  
		<pre><dir>	
		@@double operator() (double);
		</dir></pre>  
		
		  Our goal is to find x such that f(x) = target.
		  
			Examples:
			<pre><dir>
			@@	typedef double func (double);
			@@	double foo (double x) {return x+4;}	
			@@
			@@	KSecantRoot<func> solver (1); // Initial start points will be guess and guess + 1
			@@	double ans = solver.solve(foo, 0, 2); // find x such that foo(x) = 2 starting with x = 0
			@@	int numFooCalls = solver.numFCalls;
			@@
			@@	class Bar {
			@@	public: 
			@@	  bar(double x) : _x(x) {}
			@@	  double _x;
			@@	  double operator () (double y) {return y + _x;}
			@@	};
			@@
			@@	Bar bar (4);
			@@	KSecantRoot<Bar> solver2 (1);
			@@	double ans = solver2.solve(bar, 0, 2);
			</dir></pre>
			
			  The previous examples are equivalent, but illustrate the two ways the solver can be called.
*/


template <class F>
class KSecantRoot : public KOneDimSolver<F> {
public:
/** Constructor 
@param initDx  First two values passed to secant solver are guess and guess + initDx.
@param xacc  Solver runs until the two x values are within xacc of each other (except for initial two points).
@param facc  Solver runs until fabs(f(x)-target) &lt facc.
@param maxIters  Maximum number of iterations. Returns exception if exceeded.
	*/
	KSecantRoot (double initDx, double xacc = SOLVER_DOUBLE_EPS, double facc = SOLVER_DOUBLE_EPS, int maxIters = 30) 
		: _initDx(initDx),_xacc(xacc), _facc(facc), _maxIters(maxIters) {}
	
	/** Performs solver to find x such that f(x) = target starting with x = guess. */
	double solve (F& f, double guess, double target)
	{
		_numFCalls = 0;
		double x1 = guess;
		double x2 = guess + _initDx;
		double f1, f2, dx2;
		
		f1 = f(x1) - target; _numFCalls++;
		f2 = f(x2) - target; _numFCalls++;
		
		if (fabs(f1) < fabs(f2)) {	// swaps x1, x2 and f1, f2 if f1 is closer to target
			double swap;
			swap = x1;
			x1 = x2;
			x2 = swap;
			swap = f1;
			f1 = f2;
			f2 = swap;
		}
		
		// x1 is the first value, x2 is the second value
		for (int i = 0; i < _maxIters; i++)
		{
			if (f1 == f2) 
				throw KException ("KSecantRoot has x1, x2 such that f(x1)=f(x2)");
			
			dx2 = (x2-x1) / (f2-f1) * (-f2);
			x1 = x2; 
			f1 = f2;
			x2 += dx2;
			f2 = f(x2) - target; _numFCalls++;
			
			if (fabs (dx2) < _xacc || fabs(f2) < _facc) {
				_fValue = f2;
				return x2;
			}
		}
		
		throw KException ("Exceeded maximum number of iterations");
		return 0;
	}

private:
	double _initDx, _xacc, _facc;
	int _maxIters;
};

#endif
