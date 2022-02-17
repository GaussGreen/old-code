#include "drplatdep.h"
#include "drexception.h"

//f Finds root using secant search (see NRC2)
//f Secant search works well for near-linear functions
//f For highly non-linear functions, you probably want to use Brent.
//f
//f F can either be a functional or a function reference.
//f They both work, since they have the same syntax with respect to operator ()
//f
//f Functionals have the method:
//f		double operator() (double)
//f 
//f x1, x2 are the two starting x points.
//f Our goal is to find x such that f(x) = target and we will run till
//f the iterative change in x is less than xacc.
//f
//f If we exceed MAX_ITERS, we will throw an exception.

template <class F>
double DRSecantRoot (F& f, double x1, double x2, double xacc, 
					 double target = 0, int MAX_ITERS = 30)
{
	double f1, f2, dx2;

	f1 = f(x1) - target;
	f2 = f(x2) - target;

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
	for (int i = 0; i < MAX_ITERS; i++)
	{
		if (f1 == f2) 
			throw DRException ("DRSecantRoot has x1, x2 such that f(x1)=f(x2)");

		dx2 = (x2-x1) / (f2-f1) * (-f2);
		x1 = x2; 
		f1 = f2;
		x2 += dx2;
		f2 = f(x2) - target;

		if (fabs (dx2) < xacc || f2 == 0) return x2;
	}

	throw DRException ("Exceeded maximum number of iterations");
	return 0;
}





