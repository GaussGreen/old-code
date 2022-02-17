/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *	\file newtonraphson.h
 *
 *  \brief 
 *	\author  A. TRIKI
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPNUMLIB_NEWTONRAPHSON_H
#define _INGPNUMLIB_NEWTONRAPHSON_H

#include "solver.h"

///gpbase
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/port.h"
#include "gpbase/numericconstant.h"

#include "mrgk5.h"
#include "expt.h"
#include <cmath>
using namespace std;

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////
/// standard newton raphson solver
///		a functon with known derivative
///		should support T operator(T)
///		and a function T Derivative() that return a pointor to a function!
//////////////////////////////////////////

template <typename Func, typename T=double> class T_NewtonRaphsonSolver : public T_SolverWithInitialGuess<Func,T>
{
public:
	T_NewtonRaphsonSolver(const Func& f,
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
	
	~T_NewtonRaphsonSolver() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const	{ return new T_NewtonRaphsonSolver<Func,T>(*this); }
};


//////////////////////////////////////////
/// standard newton raphson solver (with no throw)
///		a functon with known derivative
///		should support T operator(T)
///		and a function T Derivative() that return a pointor to a function!
//////////////////////////////////////////

template <typename Func, typename T=double> class T_NewtonRaphsonSolverNoThrow : public T_SolverWithInitialGuess<Func,T>
{
public:
	T_NewtonRaphsonSolverNoThrow(const Func& f,
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
	
	~T_NewtonRaphsonSolverNoThrow() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const	{ return new T_NewtonRaphsonSolverNoThrow<Func,T>(*this); }
};

//////////////////////////////////////////
/// standard newton raphson solver (with no throw and bounded lenght of dx)
///		a functon with known derivative
///		should support T operator(T)
///		and a function T Derivative() that return a pointor to a function!
//////////////////////////////////////////

template <typename Func, typename T=double> 
class T_NewtonRaphsonSolverBoundedLenght : public T_SolverWithInitialGuess<Func,T>{
public:
	T_NewtonRaphsonSolverBoundedLenght(	const Func&				f,
										const T target			= 0,
										const T length			= 1.0,
										const T ftolerance		= SolverConstant::DefaultFxTolerance,
										const T xtolerance		= SolverConstant::DefaultFxTolerance,
										const size_t max_iter	= SolverConstant::DefaultMax_Iter,
										const bool printLevel	= false,
										ARM_WarningKeeper* wKeeper = 0 ):
											T_SolverWithInitialGuess<Func,T>(	f,
																				target,
																				max_iter,
																				ftolerance,
																				xtolerance,
																				printLevel,
																				wKeeper)	{}
	
	~T_NewtonRaphsonSolverBoundedLenght() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const	{ return new T_NewtonRaphsonSolverBoundedLenght<Func,T>(*this); }
	void SetLength(  const T & length) { itsLength= length; }
private:
	T itsLength;

};

//////////////////////////////////////////
/// standard newton raphson solver with retrial
///		a functon with known derivative
///		should support T operator(T)
///		and a function T Derivative() that return a pointor to a function!
//////////////////////////////////////////

template <typename Func,typename T=double,typename RandomGen=ARM_RandUniform_MRGK5> class T_NewtonRaphsonSolverRetrial : public T_SolverWithInitialGuess<Func,T>
{
public:
	T_NewtonRaphsonSolverRetrial(const Func& f,
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

	~T_NewtonRaphsonSolverRetrial() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const	{ return new T_NewtonRaphsonSolverRetrial<Func,T,RandomGen>(*this); }
};

//////////////////////////////////////////
/// smooth newton raphson solver
///		a functon with known derivative
///		should support T operator(T)
///		and a function T Derivative() that return a pointor to a function!
/// compared to the standard newton raphson solver
/// it can cope with wrong initial guess by taking a random guess
///	to find another initial point!
//////////////////////////////////////////

template <typename Func, typename T=double,typename RandomGen=ARM_RandUniform_MRGK5> class T_SmoothNewtonRaphsonSolver : public T_SolverWithInitialGuess<Func,T>
{
public:
	T_SmoothNewtonRaphsonSolver(const Func& f,
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
	~T_SmoothNewtonRaphsonSolver() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const	{ return new T_SmoothNewtonRaphsonSolver<Func,T,RandomGen>(*this); }
};

///////////////////////////////////
//// implementation
///////////////////////////////////
//////////////////////////////////////
/// standard Newton Rhapson
//////////////////////////////////////

template <typename Func, typename T> T T_NewtonRaphsonSolver<Func,T>::Solve()
{
	T df_x,dx,f_x,result;
	result=GetInitialGuess();
	for (size_t i=0;i<GetMaxIters();++i)
	{
		f_x=GetF()(result)-GetTarget();
		df_x=(*GetF().Derivative())(result);
		
		if(fabs(df_x)<SolverConstant::DefaultGradTolerance)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Slope to low to find root in <NewtonRaphsonSolver>");
		
		GetSolverDetails()->StoreSolverDetails(i, result ,f_x, df_x);
		dx=f_x/df_x;
		result -= dx;
		if (fabs(dx) < GetFxTolerance()) 
		{
			/// to keep track of things, we evaluate with the resulting value
			GetF()(result);
			GetSolverDetails()->StoreFinalSolution(result);
			return result;
		}
	}
	
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"reached max iteration in <NewtonRaphsonSolver>");
}



///////////////////////////////////
//// implementation
///////////////////////////////////
//////////////////////////////////////
/// standard Newton Rhapson
//////////////////////////////////////

template <typename Func, typename T> T T_NewtonRaphsonSolverNoThrow<Func,T>::Solve()
{


	T df_x,dx,f_x,result;
	result=GetInitialGuess();

	for (size_t i=0;i<GetMaxIters();++i)
	{
		f_x=GetF()(result)-GetTarget();
		df_x=(*GetF().Derivative())(result);	
		GetSolverDetails()->StoreSolverDetails(i, result ,f_x, df_x);
		if ( fabs(df_x) <SolverConstant::DefaultGradTolerance )
			return result;
		dx=f_x/df_x;
		result -= dx;

		if ( fabs(dx) < GetFxTolerance()) 
		{
			/// to keep track of things, we evaluate with the resulting value
			GetF()(result);
			GetSolverDetails()->StoreFinalSolution(result);
			return result;
		}
	}
	if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( "problem in NewtonRaphsonSolverNoThrow: reach max iteration\n" );
	return result;
}	


///////////////////////////////////
//// implementation
///////////////////////////////////
//////////////////////////////////////
/// standard Newton Rhapson
//////////////////////////////////////

template <typename Func, typename T> T T_NewtonRaphsonSolverBoundedLenght<Func,T>::Solve()
{


	T df_x,dx,f_x,result;
	result=GetInitialGuess();

	for (size_t i=0;i<GetMaxIters();++i)
	{
		f_x=GetF()(result)-GetTarget();
		df_x=(*GetF().Derivative())(result);	
		GetSolverDetails()->StoreSolverDetails(i, result ,f_x, df_x);
		if ( fabs(df_x) <SolverConstant::DefaultGradTolerance )
			return result;
		dx=f_x/df_x;
		dx = abs(dx)>itsLength? (dx>0?1.0:-1.0)*itsLength: dx;

		result -= dx;

		if ( fabs(dx) < GetFxTolerance() ||  fabs(f_x) < GetFxTolerance() ) 
		{
			/// to keep track of things, we evaluate with the resulting value
			GetF()(result);
			GetSolverDetails()->StoreFinalSolution(result);
			return result;
		}
	}
	if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( "problem in NewtonRaphsonSolverNoThrow: reach max iteration\n" );
	return result;
}


//////////////////////////////////////
///Newton Rhapson Solver Retrial
//////////////////////////////////////

template <typename Func, typename T, typename RandomGen>  T T_NewtonRaphsonSolverRetrial<Func,T,RandomGen>::Solve()
{
	T df_x,dx,f_x,result;
	result=GetInitialGuess();
	size_t iter = 0.0;
	RandomGen unifGen;
	int rtnSup=0;

	for (size_t i=0;i<GetMaxIters();++i)
	{
		f_x=GetF()(result)-GetTarget();
		df_x=(*GetF().Derivative())(result);
		if(fabs(df_x) < SolverConstant::DefaultGradTolerance)
		{
			++iter;
			double unif = unifGen();
			result=(1.0-unif)*GetLowerBound()+unif*GetUpperBound();
			if (iter == 5)
			{
				if( GetWarningKeeper() )
					GetWarningKeeper()->AddWarningMessage( "problem in NewtonRaphsonSolverRetrial: use initial solution\n" );
				return GetInitialGuess();
			}
		}
		else
		{
			GetSolverDetails()->StoreSolverDetails(i, result ,f_x, df_x);
			dx=f_x/df_x;
			result -= dx;
			if(result>GetUpperBound())
			{
				if( GetWarningKeeper() )
					GetWarningKeeper()->AddWarningMessage( "problem in NewtonRaphsonSolverRetrial: found result above upper boundary\n" );
				result=GetUpperBound();
			}

			if(result<GetLowerBound())
			{	
				if( GetWarningKeeper() )
					GetWarningKeeper()->AddWarningMessage( "problem in NewtonRaphsonSolverRetrial: found result above lower boundary\n" );
				result=GetLowerBound();
			}
			if (fabs(dx) < GetFxTolerance()) 
			{
				/// to keep track of things, we evaluate with the resulting value
				GetF()(result);
				GetSolverDetails()->StoreFinalSolution(result);
				return result;
			}
		}
	}
	if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( "problem in NewtonRaphsonSolverRetrial: reach max iteration\n" );
	//Updates Values in Curve Model Param
	GetF()(result);
	GetSolverDetails()->StoreFinalSolution(result);
	return result;
}



//////////////////////////////////////
/// smooth Newton Rhapson
/// with 
//////////////////////////////////////
template <typename Func, typename T, typename RandomGen> T T_SmoothNewtonRaphsonSolver<Func,T,RandomGen>::Solve()
{
	T df_x,dx,f_x;  
	T result=GetInitialGuess();
	T x0=result;
	RandomGen unifGen;
	// to store  details solving
	GetSolverDetails()->StoreTitleDetails();
	
	int iter = 0.0;
	for (size_t i=1;i<=GetMaxIters();i++)
	{
		f_x=GetF()(result)-GetTarget();
		df_x=(*GetF().Derivative())(result);
		T rtn0=result;
        GetSolverDetails()->StoreSolverDetails(i, result ,f_x, df_x);
		if(fabs(df_x) < SolverConstant::DefaultGradTolerance)
		{
			++iter;
			double unif = unifGen();
			result=(1.0-unif)*GetLowerBound()+unif*GetUpperBound();
			if (iter == 5)
			{
				if( GetWarningKeeper() )
					GetWarningKeeper()->AddWarningMessage( "problem in SmoothNewtonRaphsonSolver: use initial value\n" );
				return x0;
			}
		}
		else
		{
			dx=f_x/df_x;
			result-= dx;
			iter=0;
			x0=result;
			
			if (fabs(dx)<GetFxTolerance()) 
			{
				/// to keep track of things, we evaluate with the resulting value
				GetF()(result);
				GetSolverDetails()->StoreFinalSolution(result);
				return result;
			}
		} 
		
	    if (result > 1.25*rtn0)
			result = 1.25*rtn0;     
		else if (result < 0.75*rtn0)
			result = 0.75*rtn0;
	}

	if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( "problem in SmoothNewtonRaphsonSolver: reach max iteration\n" );
	//Updates Values in Curve Model Param
	GetF()(result);
	GetSolverDetails()->StoreFinalSolution(result);
	return result;
}

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
