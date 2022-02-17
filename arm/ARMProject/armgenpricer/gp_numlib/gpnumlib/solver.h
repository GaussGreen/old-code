/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file solver.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPNUMLIB_SOLVER_H
#define _INGPNUMLIB_SOLVER_H

///GP Base
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/ostringstream.h"
#include "gpbase/warningkeeper.h"

#include "mrgk5.h"
#include "expt.h"
#include <cmath>
#include <iomanip> /// for setprecision()
using namespace std;

CC_BEGIN_NAMESPACE( ARM )

/// general constant
struct SolverConstant
{
	static const double DefaultXTolerance;
	static const double DefaultFxTolerance; 
	static const double DefaultGradTolerance;
	static const double DoubleTolerance;
	static const size_t DefaultMax_Iter;
	static const double DefaultDichoXTolerance;
	static const size_t DefaultDichoMax_Iter;
	static const double UpperInfiniteBound;
	static const double LowerInfiniteBound;
	static const double DefaultStepMax;
	static const bool   DefaultPrint;

};

///////////////////////////////////////////////////
///	Struct : SolverDetails 
///	Action : Functor whose responsability is to 
///			 store details in the file Method Description
////////////////////////////////////////////////////

class SolverDetails
{
	protected:
	CC_Ostringstream os;
	public:
// FIXMEFRED: mig.vc8 (29/05/2007 09:24:49): missing default ctor leads to intern C++ ctor (basic_ios error)
	SolverDetails()
	{}

// FIXMEFRED: mig.vc8 (29/05/2007 09:24:49): missing default ctor leads to intern C++ ctor (basic_ios error)
	SolverDetails(const SolverDetails& cur)
	{}

	virtual inline const CC_Ostringstream& GetOs() const	{	return os;	}
	virtual void StoreTitleDetails() =0;
	virtual void StoreSolverDetails(int iter,double x,double f, double df)=0;
	virtual void StoreFinalSolution(double result)=0;
	virtual SolverDetails* Clone() const= 0;
};

///////////////////////////////////////////////////
///	Struct : SolverDetailsOn 
///	Action : 
////////////////////////////////////////////////////

class SolverDetailsOn:public SolverDetails
{
	virtual SolverDetails* Clone() const	{ return new SolverDetailsOn(*this); }
	virtual void StoreTitleDetails()
	{
		const string indent="\t";
		const string dbIndent=indent+indent;
		os <<"Iter"<<indent<<"X"<<dbIndent<<"F"<<dbIndent<<"dF\n";
	}
	virtual void StoreSolverDetails(int iter, double x,double f, double df)
	{
		const string indent="\t";
		os <<fixed<<setprecision(6)<<iter<<indent<<fixed<<setprecision(12)<<x<<indent<<fixed<<setprecision(12)<<f<<indent<<fixed<<setprecision(12)<<df<<"\n";
	}
	virtual void StoreFinalSolution(double result)
	{
		os<<"\n"<<"Final Solution : \t"<<fixed<<setprecision(12)<<result<<"\n";
	}
};

///////////////////////////////////////////////////
///	Struct : SolverDetailsOff 
///	Action : 
////////////////////////////////////////////////////

class SolverDetailsOff:public SolverDetails
{
	virtual SolverDetails* Clone() const	{ return new SolverDetailsOff(*this); }
	virtual void StoreTitleDetails(){}
	virtual void StoreSolverDetails(int iter, double x,double f, double df)	{}
	virtual void StoreFinalSolution(double result){}
};

///////////////////////////////////////////
/// \brief
///	in this file, you can find the following solvers:
///		-dichotomy solver
///		-brent solver
///		-standard newton raphson
///		-smooth newton raphson
///////////////////////////////////////////


template <typename Func, typename T=double> class T_SolverWithInitialGuess
{
private:
	/// function supplied by the user, 
	/// must calculate the value F(x) 
	/// at any point x  in [a,b]
	const Func& itsF;

	/// target to reach
	T itsTarget;

	/// the point at with the value of F is required
	T itsInitialGuess;

	/// the lower bound of the interval, a
	T itsLowerBound;

	/// the upper bound of the interval, b
	T itsUpperBound;

	/// max of iterations
	size_t itsMaxIters;

	/// the value such that |F(x)|< ftol,
	/// x is accepted as the zero
	T itsFxTolerance;

	/// the absolute tolerance to which the zerao is required
	T itsXTolerance;

	/// to return if or not zero were required
	CC_IS_MUTABLE bool itsRootState;

	/// to controle warning
	ARM_WarningKeeper* itsWarningKeeper;

	/// to store all details
	SolverDetails* itsFtor;

    /// to pre-optimize with another solver
     ModifiedNRSolverPtr itsSolver;

public:
	T_SolverWithInitialGuess(const Func& F,
		const T target = 0, 
		const size_t max_iter = SolverConstant::DefaultMax_Iter,
		const T fxTolerance	= SolverConstant::DefaultFxTolerance,
		const T xTolerance    = SolverConstant::DefaultXTolerance,
		bool printLevel = false,
		ARM_WarningKeeper* wKeeper = NULL,
        ModifiedNRSolverPtr& solver = ModifiedNRSolverPtr(NULL))
	:
	itsF(F),
	itsTarget(target),
	itsInitialGuess(0.0),
	itsLowerBound(SolverConstant::LowerInfiniteBound),
	itsUpperBound(SolverConstant::UpperInfiniteBound),
	itsMaxIters(max_iter),
	itsFxTolerance(fxTolerance),
	itsXTolerance(xTolerance),
	itsWarningKeeper(wKeeper),
    itsSolver(solver),
	itsRootState(false)
	{
		if(printLevel)
			itsFtor =new SolverDetailsOn;
		else
			itsFtor = new SolverDetailsOff;
	}
	T_SolverWithInitialGuess(const T_SolverWithInitialGuess<Func,T>& rhs)
	:
		itsTarget(rhs.itsTarget),
		itsInitialGuess(rhs.itsInitialGuess),
		itsLowerBound(rhs.itsLowerBound),
		itsUpperBound(rhs.itsUpperBound),
		itsMaxIters(rhs.itsMaxIters),
		itsFxTolerance(rhs.itsFxTolerance),
		itsXTolerance(rhs.itsXTolerance),
		itsRootState(rhs.itsRootState),
		itsF(rhs.itsF)
	{
		itsWarningKeeper = rhs.itsWarningKeeper ? new ARM_WarningKeeper( *rhs.itsWarningKeeper): NULL;
		itsFtor= rhs.itsFtor ? (SolverDetails*)rhs.itsFtor->Clone(): NULL;
        itsSolver= ModifiedNRSolverPtr(rhs.itsSolver != ModifiedNRSolverPtr(NULL) ? (ModifiedNRSolver*)rhs.itsFtor->Clone(): NULL);
	}

	T_SolverWithInitialGuess<Func,T>& operator = (const T_SolverWithInitialGuess<Func,T>& rhs)
	{
		if(this != rhs)
		{
			itsTarget		= rhs.itsTarget;
			itsInitialGuess = rhs.itsInitialGuess;
			itsLowerBound	= rhs.itsLowerBound;
			itsUpperBound	= rhs.itsUpperBound;
			itsMaxIters		= rhs.itsMaxIters;
			itsFxTolerance	= rhs.itsFxTolerance;
			itsXTolerance	= rhs.itsXTolerance;
			itsRootState	= rhs.itsRootState;
			itsF			= rhs.itsF;
			delete itsWarningKeeper;
			itsWarningKeeper = rhs.itsWarningKeeper ? new ARM_WarningKeeper( *rhs.itsWarningKeeper): NULL;
			delete itsFtor;
			itsFtor= rhs.itsFtor ? (SolverDetails*)rhs.itsFtor->Clone(): NULL;
            itsSolver= ModifiedNRSolverPtr(rhs.itsSolver != ModifiedNRSolverPtr(NULL) ? (ModifiedNRSolver*)rhs.itsFtor->Clone(): NULL);
		}
	}

	virtual ~T_SolverWithInitialGuess() { delete itsFtor;}
	virtual T Solve()= 0;
	void setInitialGuess(T initialGuess, 
		T lowerBound = SolverConstant::LowerInfiniteBound,
		T upperBound = SolverConstant::UpperInfiniteBound )
	{
		itsInitialGuess	= initialGuess;
		itsLowerBound=lowerBound;
		itsUpperBound=upperBound;
	}

	virtual T_SolverWithInitialGuess<Func,T>* Clone() const= 0;

	/// Accessors
	inline SolverDetails* GetSolverDetails() const		{ return itsFtor;			}
	inline bool GetRootState() const					{ return itsRootState;		}
	inline void SetRootState(bool rootState) const		{ itsRootState = rootState;	}
	inline const Func&  GetF() const 					{ return itsF;				} 
	inline T GetTarget() const							{ return itsTarget;			}
	inline T GetInitialGuess()  const					{ return itsInitialGuess;   }
	inline T GetLowerBound() const						{ return itsLowerBound;     }
	inline T GetUpperBound() const						{ return itsUpperBound;     }
	inline size_t GetMaxIters() const					{ return itsMaxIters;		}
	inline T GetFxTolerance() const						{ return itsFxTolerance;	}
	inline T GetXTolerance() const						{ return itsXTolerance;		}
	inline ARM_WarningKeeper* GetWarningKeeper() const  { return itsWarningKeeper;	}
    inline ModifiedNRSolverPtr GetSolver() const        { return itsSolver;         }

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
