/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amc_andersen.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPINFRA_AMCANDERSEN_H
#define _INGPINFRA_AMCANDERSEN_H

#include "amc_exercboundcalc.h"
#include "gpbase/port.h"
#include "gpinfra/exerciseboundary.h"

#include <algorithm>
#include <functional>

using std::pair;

CC_BEGIN_NAMESPACE( ARM )

class ARM_AMCAndersen : public ARM_ExerciseBoundaryCalc
{
private: 
	bool itsUseSortedMaximisation;
	pair<double,bool> ComputeStdExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_VectorPtr& trigValue);
	pair<double,bool> ComputeSortedExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_VectorPtr& trigValue);

public: 
	/// the central function
	virtual ARM_ExerciseBoundary * ComputeExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt,
		const ARM_GP_MatrixPtr& StatesVector );

	/// Constructors/Destructor
	ARM_AMCAndersen( size_t ItersNb, bool sortedMaximization =false );
	ARM_AMCAndersen(const ARM_AMCAndersen& rhs)
	: ARM_ExerciseBoundaryCalc(rhs), itsUseSortedMaximisation(rhs.itsUseSortedMaximisation) {};

	ARM_AMCAndersen& operator=(const ARM_AMCAndersen& rhs )
	{
		if( this != &rhs )
		{
			ARM_ExerciseBoundaryCalc::operator =( rhs );
			itsUseSortedMaximisation = rhs.itsUseSortedMaximisation;
		}
		return *this;
	}
	~ARM_AMCAndersen() {};

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// part to tell if needs default argument for the exercise function
	virtual bool NeedToCreateDefaultArgument() const { return false; }

	/// the root name
	virtual ARM_CLASS_NAME GetRootName() { return ARM_AMCANDERSEN; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
