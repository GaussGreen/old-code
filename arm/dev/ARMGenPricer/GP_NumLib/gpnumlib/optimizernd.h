/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *	\file optimizernd.h
 *
 *  \brief 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date December 2005
 */

#ifndef _INGPNUMLIB_OPTIMIZERND_H
#define _INGPNUMLIB_OPTIMIZERND_H

#include "gpbase\assignop.h"
#include "gpbase\numericconstant.h"
#include "gpbase\functor.h"
#include "gpbase\gpvector.h"
#include "gpbase\typedef.h"

CC_BEGIN_NAMESPACE( ARM )

// This class is a wrapper of nag optimizer for multidimensional functions
class NagNDOptimizer
{
public:

	// Multi dimension function type
	typedef ARM_GP::UnaryFunc<ARM_GP_Vector, double> NDFunc;

	// Main constructor
	NagNDOptimizer(
		NDFunc* func,
        const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER,
		double tolerance						= 1.0e-3,
		double stepMax							= 1.0e+002,
		bool localSearch						= false);

	// Copy constructor, assignment operator and Clone Function
	NagNDOptimizer(const NagNDOptimizer& rhs);

	ASSIGN_OPERATOR(NagNDOptimizer)

	/*virtual NagNDOptimizer* Clone() const
	{
		return new(*this);
	}*/

	// Compute the optimal point
	virtual ARM_GP_Vector Optimize(
		const ARM_GP_Vector& initialGuess,
		const ARM_GP_Vector& lowerBound,
		const ARM_GP_Vector& upperBound);

	// Accessor
	NDFunc* GetFunction() const {return itsFunc; }

private:
	size_t itsMaxIter;
	bool itsGetDetails;
	double itsTolerance;
	double itsStepMax;
	bool itsLocalSearch;

	ARM_GP_Vector itsInitialGuess;
	NDFunc* itsFunc;
};


CC_END_NAMESPACE()

#endif _INGPNUMLIB_OPTIMIZERND_H