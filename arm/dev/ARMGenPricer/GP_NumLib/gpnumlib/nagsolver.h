/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *	\file solver.h
 *
 *  \brief 
 *	\author  A. E.EZZINE
 *	\version 1.0
 *	\date June 2005
 */

#ifndef _INGPNUMLIB_NAGSOLVER_H
#define _INGPNUMLIB_NAGSOLVER_H


/// nag headers
#include "gpbase/removenagwarning.h"
#include "gpbase/functor.h"
#include "nag.h"
#include "nage04.h"
#include "gpnumlib/numfunction.h"

#include "solver.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/port.h"
#include "gpbase/warningkeeper.h"
#include "mrgk5.h"
#include <glob/expt.h>
#include <cmath>
#include <iomanip> /// for setprecision()
#include <sstream>

using namespace std;

CC_BEGIN_NAMESPACE( ARM )

inline void ProcessNagWarningOnObj( NagError fail, ARM_WarningKeeper& wk )
{
    switch( fail.code )
	{
	case NE_USER_STOP:
		wk.AddWarningMessage( "User requested termination\n" );
		break;
	case NW_COND_MIN:
		wk.AddWarningMessage( "Conditional minumum: A Solution has been found, but may be it's not the best one!\n" );
		break;
	case NE_INT_ARG_LT:
		wk.AddWarningMessage( "On entry, n must not be less than 1\n" );
		break;
	case NE_BOUND:
		wk.AddWarningMessage( "The lower bound is greater than the upper bound\n" );
		break;
	case NE_DERIV_ERRORS:
		wk.AddWarningMessage( "Large errors were found in the derivatives of the objective function\n" );
		break;
     
	case NE_OPT_NOT_INIT:
		wk.AddWarningMessage( "Options structure not initialized\n" );
		break;
	case NE_BAD_PARAM:
		wk.AddWarningMessage( "On entry parameter bound had an illegal value\n" );
		break;
	case NE_2_REAL_ARG_LT:
		wk.AddWarningMessage( "On entry, options.step_max mispecified with options.optim_tol\n" );
		break;
    case NE_INVALID_INT_RANGE_1:
		wk.AddWarningMessage( "Value given to options.max_iter not valid\n" );
		break;
	case NE_INVALID_REAL_RANGE_EF:
		wk.AddWarningMessage( "Value given to options.optim_tol not valid\n" );
		break;
	case NE_INVALID_REAL_RANGE_FF:
		wk.AddWarningMessage( "Value given to options.linesearch_tol not valid\n" );
		break;
	case NE_INIT_MEM:
		wk.AddWarningMessage( "Option init_state incorrect\n" );
		break;
	case NE_NO_MEM:
		wk.AddWarningMessage( "Option init_state incorrect\n" );
		break;
	case NE_HESD:
		wk.AddWarningMessage( "initial values of the supplied options.hesd incorrect\n" );
		break;
	case NE_ALLOC_FAIL:
		wk.AddWarningMessage( "Memory allocation failed\n" );
		break;
	case NW_TOO_MANY_ITER:
		wk.AddWarningMessage( "Maximum number of iterations reached\n" );
		break;
	case NE_CHOLESKY_OVERFLOW:
		wk.AddWarningMessage( "Cholesky overflow\n" );
		break;
	case NW_LOCAL_SEARCH:
		wk.AddWarningMessage( "local search has failed to find a feasible point\n" );
		break;
	case NE_NOT_APPEND_FILE:
		wk.AddWarningMessage( "Cannot open file \n" );
		break;
	case NE_WRITE_ERROR:
		wk.AddWarningMessage( "Cannot write in file \n" );
		break;
	case NE_NOT_CLOSE_FILE:
		wk.AddWarningMessage( "Cannot close file \n" );
		break;
	case NE_NOERROR : case NW_KT_CONDITIONS:
		break;
	default:
		std::ostringstream os;
		os << "Exited with an unknown flag error " << fail.code << "\n";
		wk.AddWarningMessage( os.str() );
		break;
	}
}



//////////////////////////
/// Nag solver
/// a function should support T operator(T);
//////////////////////////
template <typename Func, typename T=double> class T_NagSolver : public T_SolverWithInitialGuess<Func,T>
{
public:
	T_NagSolver(const Func& f,
		const T target	= 0,
		const T ftolerance	= SolverConstant::DefaultFxTolerance,
		const T xtolerance	= SolverConstant::DefaultFxTolerance,
		const size_t max_iter	= SolverConstant::DefaultMax_Iter,
		const bool printLevel	= false,
		ARM_WarningKeeper* wKeeper = 0 )
	:
		T_SolverWithInitialGuess<Func,T>(f,
			target,
			max_iter,
			ftolerance,
			xtolerance,
			printLevel,
			wKeeper)
	{}
	~T_NagSolver() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const { return new T_NagSolver<Func,T>(*this); } 
};


///////////////////////////////////
//// implementation
///////////////////////////////////

////////////////////////////////////////////////////
///	Routine: WeightedSquared
///	Returns: void
///	Action : global function used to call Nag optimizer 
////////////////////////////////////////////////////
void NAG_CALL WeightedSquared(double x, 
                              double* fx,
                              Nag_Comm *comm);

////////////////////////////////////////////////////
///	Class  : T_NagSolver
///	Routine: Solve()
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
template <typename Func, typename T> T T_NagSolver<Func,T>::Solve()
{
    T bl = GetLowerBound();
	T bu = GetUpperBound();
	T x = GetInitialGuess();

    Nag_Comm comm;
    comm.p = (void*) &GetF();
    static NagError fail;
    fail.print = FALSE;
	
	T e1 =GetFxTolerance();
	T e2 = GetXTolerance();
	int max_func = 3*GetMaxIters();
	/// (e04abc) searches for a minimum, in a given finite interval
	/// of a continous function of a single variable,using function values only
	/// The method (based on quadratic interpolation) is intended for functions
	/// which have a continous first derivative.
	T err;
    nag_opt_one_var_no_deriv(&WeightedSquared,e1,e2,&bl,&bu,max_func,&x, &err,&comm, &fail);

	if( GetWarningKeeper() )
		ProcessNagWarningOnObj( fail, *GetWarningKeeper() );
	T result = x;
	return result;
}

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
