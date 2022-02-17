/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 * $Log: solver.cpp,v $
 * Revision 1.1  2004/09/22 10:15:09  aschauly
 * Initial revision
 *
 *
 *
 */


/*! \file optimizernd.cpp
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date December 2005
 */

#include "gpbase/removenagwarning.h"
#include "gpnumlib/optimizernd.h"
#include "nag.h"
#include "nage04.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
/// Class: NagNDOptimizer
///	Routine: WeightedSquaredND
///	Returns: void
///	Action : global function used to call Nag optimizer 
////////////////////////////////////////////////////
NagNDOptimizer::NagNDOptimizer(
NDFunc* func,
const size_t max_iter,
double tolerance,
double stepMax,
bool localSearch)
:
itsFunc(func),
itsMaxIter(max_iter),
itsTolerance(tolerance),
itsStepMax(stepMax),
itsLocalSearch(localSearch)
{
}

////////////////////////////////////////////////////
/// Class: NagNDOptimizer
///	Routine: operator=
///	Returns: void
///	Action : global function used to call Nag optimizer 
////////////////////////////////////////////////////
NagNDOptimizer::NagNDOptimizer(const NagNDOptimizer& rhs)
:
itsMaxIter(rhs.itsMaxIter),
itsTolerance(rhs.itsTolerance),
itsStepMax(rhs.itsStepMax),
itsLocalSearch(rhs.itsLocalSearch)
{
}


////////////////////////////////////////////////////
///	Routine: WeightedSquaredND
///	Returns: void
///	Action : global function used to call Nag optimizer 
////////////////////////////////////////////////////
void NAG_CALL NAGWeightedSquaredND(
	Integer n,
	double* x, 
    double* fx,
	double* g,
    Nag_Comm *comm )
{
    NagNDOptimizer* Optimise = (NagNDOptimizer*)(comm->p);
	ARM_GP_Vector var(n,x);
	*fx = (*Optimise->GetFunction())(var);
}

////////////////////////////////////////////////////
/// Class: NagNDOptimizer
///	Routine: ARM_GP_Vector
///	Returns: void
///	Action : Compute the optimal point
////////////////////////////////////////////////////
ARM_GP_Vector NagNDOptimizer::Optimize(
const ARM_GP_Vector& initialGuess,
const ARM_GP_Vector& lowerBound,
const ARM_GP_Vector& upperBound)
{
	double err;
	Nag_BoundType bound;
	bound = Nag_Bounds;
	int NbVar   = initialGuess.size();  
	double* x	= new double [NbVar];
	double* g	= new double [NbVar];
	double* bl	= new double [NbVar];
	double* bu	= new double [NbVar];

	size_t i;
	for (i = 0; i<NbVar; ++i)
	{
		bl[i] = lowerBound[i];
		bu[i] = upperBound[i];
		x[i] = initialGuess[i];
	}

	/// initialise NAG (e04xxc)
	Nag_E04_Opt options;	
	nag_opt_init(&options);

	Nag_Comm comm;
	comm.p = this;
	static NagError fail;

	fail.print					= FALSE;
    options.list				= FALSE;
    options.print_level			= Nag_NoPrint;
    options.output_level		= Nag_NoOutput;
	options.minor_print_level	= Nag_NoPrint;
	options.print_deriv			= Nag_D_NoPrint;


	options.max_iter	= itsMaxIter;
	options.optim_tol	= itsTolerance;
	options.step_max	= itsStepMax;
	options.local_search= itsLocalSearch;

	/// optimisation with boundary using no derivatives (e04jbc)
	nag_opt_bounds_no_deriv(NbVar,&NAGWeightedSquaredND,bound,bl,bu,x, &err, g,&options, &comm, &fail);

	/// free the option (e04xzc)
	nag_opt_free(&options,"all",&fail);

	ARM_GP_Vector ret(NbVar);

	for (i = 0; i<NbVar; ++i)
		ret[i] = x[i];

	/// free memory
	delete g;
	delete bl;
	delete bu;
	delete x;

	return ret;
}

CC_END_NAMESPACE()