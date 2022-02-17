/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amc_exercfrcalc.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPINFRA_EXERCFRCALC_H
#define _INGPINFRA_EXERCFRCALC_H

#include "gpbase/port.h"
#include "gpinfra/exerciseboundary.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_ExerciseBoundaryCalc : public ARM_RootObject
{
protected: 
	size_t itsItersNb;

public: 
	/// Constructor - Destructor
	ARM_ExerciseBoundaryCalc() {};
	~ARM_ExerciseBoundaryCalc() {};

	/// part to compute an exercise boundary
	virtual ARM_ExerciseBoundary * ComputeExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, 
		const ARM_GP_MatrixPtr& StatesVector  ) = 0;
	/// part to tell if needs default argument for the exercise function
	virtual bool NeedToCreateDefaultArgument() const = 0;
	// Do we use the model states for the exercise boundary calculation
	virtual bool IsAutomatic() const { return false; };

	/// accessors
	inline size_t getItsItersNb() { return itsItersNb; };

	// the root name
	virtual ARM_CLASS_NAME GetRootName() = NULL;

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
