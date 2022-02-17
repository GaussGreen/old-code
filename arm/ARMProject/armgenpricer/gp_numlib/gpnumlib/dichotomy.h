/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *	\file dichotomy.h
 *
 *  \brief 
 *	\author  A. TRIKI
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPNUMLIB_DICHOTOMY_H
#define _INGPNUMLIB_DICHOTOMY_H

#include "solver.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/port.h"
#include "mrgk5.h"
#include "expt.h"
#include <cmath>

using namespace std;

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////
/// dichotomy solver
/// a function should support T operator(T);
//////////////////////////

template <typename Func, typename T=double> class T_DichotomySolver : public T_SolverWithInitialGuess<Func,T>
{
private:
	///Test if root exist 
	void IsRootExisting() const;
public:
	T_DichotomySolver( const Func& f,
		const T target	= 0,
		const T ftolerance	= SolverConstant::DefaultFxTolerance,
		const T xtolerance	= SolverConstant::DefaultFxTolerance,
		const size_t max_iter	= SolverConstant::DefaultMax_Iter,
		const bool printLevel	= false,
		ARM_WarningKeeper* wKeeper = NULL ,
        ModifiedNRSolverPtr& solver = ModifiedNRSolverPtr(NULL))
	:
		T_SolverWithInitialGuess<Func,T>(f,
			target,
			max_iter,
			ftolerance,
			xtolerance,
			printLevel,
			wKeeper,
            solver)
	{}
	/// virtual functions
	~T_DichotomySolver() {}
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const	{ return new T_DichotomySolver<Func,T>(*this);	}
};


///////////////////////////////////////////////////////////////////////////
///	Class  : template T_DichotomySolver
///	Routine: IsRootExisting
///	Returns: nothing
///	Action : We use bounds to test if root exist
///////////////////////////////////////////////////////////////////////////
template <typename Func, typename T> void T_DichotomySolver<Func,T>::IsRootExisting() const
{
	T fmin=GetF()(GetLowerBound())- GetTarget();
	T fmax=GetF()(GetUpperBound())- GetTarget();
	bool RootState = (fmin*fmax<=0.0);
	SetRootState(RootState);
}

///////////////////////////////////////////////////////////////////////////
///	Class  : template T_DichotomySolver
///	Routine: Solve
///	Returns: T
///	Action : Dichotomy routine
///////////////////////////////////////////////////////////////////////////
template <typename Func, typename T> T T_DichotomySolver<Func,T>::Solve()
{
    T result=0.0;;
	IsRootExisting();
	if (!GetRootState())
	{
        CC_Ostringstream Os;
        Os << "Could not bracket any solution between"<<GetLowerBound()<<"'"<<GetUpperBound();
        Os <<"in DichotomySolver::Solve"<< "\n";
        if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( Os.str());

		return result;
	}

	T fmin=GetF()(GetLowerBound())- GetTarget();
	
	T dx,xmid;
	result = fmin < 0.0 ? (dx=GetUpperBound()-GetLowerBound(),GetLowerBound()) : (dx=GetLowerBound()-GetUpperBound(),GetUpperBound());
	for(size_t i=1; i<=GetMaxIters(); i++) 
	{
		dx*=0.5;
		xmid=result+dx;
		T fmid=GetF()(xmid)-GetTarget();
		if (fmid <= 0.0) 
			result=xmid;
		if (fabs(dx) < GetXTolerance() || fabs(fmid) < GetFxTolerance() )
            break;
	}
	
	if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( "problem in Dichotomy: reach max iteration\n" );

    /// call solver to find root with high precision
    if(GetSolver()!= ModifiedNRSolverPtr(NULL))
    {
        GetSolver()->setInitialGuess(result);
	    /// Then we call the final solver
	    return GetSolver()->Solve();
    }
    
    return result;
}

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
