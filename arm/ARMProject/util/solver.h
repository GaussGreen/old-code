/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: solver.h,v $
 * Revision 1.11  2004/03/15 11:21:36  jmprie
 * ajout SetGuessBounds()
 *
 * Revision 1.10  2004/01/21 12:05:23  emezzine
 * added boundaries into solver.
 *
 * Revision 1.8  2004/01/07 07:41:54  jmprie
 * ajout du zero Brent solver
 *
 * Revision 1.7  2003/12/09 17:27:12  emezzine
 * improvment
 *
 * Revision 1.6  2003/11/17 13:48:27  emezzine
 * New class ModifiedNRSolver.
 *
 * Revision 1.5  2003/10/17 14:03:27  jmprie
 *  chgt de fonctionnement du NewtonRaphsonSolver pour pouvoir modifier
 * le contexte de la f(x) sans toucher au solveur
 *
 * Revision 1.4  2003/10/13 07:56:39  jmprie
 *  ajout pour Newton, du SetInitialGuess()
 *
 * Revision 1.3  2003/08/18 11:16:32  ebenhamou
 * exception is sent in the proper format
 *
 */

           
#ifndef SOLVER_H
#define SOLVER_H

#include "function.h"


#define SOLVE_PRECISION 0.0001
#define FIN_MAX_ITER_SOLVE 2000


class Solver{
public:
  virtual ~Solver(){};
  virtual double Solve() = 0;
};

class DichotomySolver: public Solver
{
public:
	DichotomySolver( const RealValuedFunction& f_, 
		             double minorant_,
		             double majorant_,
		             double right_member=0.0,
		             double precision_=SOLVE_PRECISION, 
		             long max_iters_=FIN_MAX_ITER_SOLVE );
	double Solve();		// virtual function
	
private:
	const RealValuedFunction& f;
	double minorant, majorant;
	double b;//on cherche les solutions de f(x)=b
	double precision;
	long max_iters;
};



class NewtonRaphsonSolver: public Solver
{
public:
  NewtonRaphsonSolver( const FunctionWithKnownDerivative* f_, 
		               double initial_guess_,
		               double right_member=0.0,
		               double precision_=SOLVE_PRECISION, 
		               long max_iters_=FIN_MAX_ITER_SOLVE );
  ~NewtonRaphsonSolver();
  virtual double Solve();
  void setInitialGuess(double guess);

protected:
    // Pointers to manage outside the solver the data
    // context of f and df
  const RealValuedFunction* f;
  const RealValuedFunction* fprime;
  double initial_guess;
  double b;
  double precision;
  long max_iters;
};


class ModifiedNRSolver: public Solver
{
public:
  ModifiedNRSolver( const FunctionWithKnownDerivative* f_, 
		               double initial_guess_,
		               double right_member=0.0,
		               double precision_=SOLVE_PRECISION, 
		               long max_iters_=FIN_MAX_ITER_SOLVE );
  ~ModifiedNRSolver();
  virtual double Solve();
  void SetInitialParams(double guess,
                        double lowerbound = LOW_INFINITE_BOUND,
                        double upperbound =UPPER_INFINITE_BOUND);

protected:
    // Pointers to manage outside the solver the data
    // context of f and df
  const RealValuedFunction* itsFunction;
  const RealValuedFunction* itsFunctionDerivative;
  double initial_guess;
  double lower_bound;
  double upper_bound;
  double b;
  double precision;
  long max_iters;
};


class BrentSolver: public Solver
{
public:
  BrentSolver( const RealValuedFunction& f_, 
		       double minGuess_,
		       double maxGuess_,
		       double target_=0.0,
		       double precision_=SOLVE_PRECISION, 
		       long max_iters_=FIN_MAX_ITER_SOLVE );
  ~BrentSolver();

  void SetGuessBounds(const double minVal,const double maxVal) {minGuess=minVal; maxGuess=maxVal;}

  virtual double Solve();

protected:
private:
	const RealValuedFunction& f;
    double minGuess,maxGuess;
	double target;// f(x)=target
	double precision;
	long max_iters;
};



#undef SOLVE_PRECISION 
#undef FIN_MAX_ITER_SOLVE

#endif

