/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amc_LongstaffSchwartz.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPINFRA_AMCLongstaffSchwartz_H
#define _INGPINFRA_AMCLongstaffSchwartz_H

#include "gpbase/port.h"
#include "amc_exercboundcalc.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_AMCLongstaffSchwartz : public ARM_ExerciseBoundaryCalc
{
private: 
	ARM_VectorPtr& ComputeRegression( const ARM_VectorPtrDbleFuncPtrVector& Fcts, 
		const ARM_VectorPtrVector& X, const ARM_VectorPtr& Y );

	void CopyNoCleanUp(const ARM_AMCLongstaffSchwartz& rhs);
    void CleanUp();
public:
	ARM_ExerciseBoundary * ComputeExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt,
		const ARM_GP_MatrixPtr& StatesVector);
	 
	ARM_AMCLongstaffSchwartz( size_t ItersNb, ARM_Regression::RegressionMode regressionMode = ARM_Regression::LS, double span = 0.0, int useModelStates = false, int degree = 3);
	~ARM_AMCLongstaffSchwartz() {};

	/// part to tell if needs default argument for the exercise function
	virtual bool NeedToCreateDefaultArgument() const { return true; }
	virtual bool IsAutomatic() const { return itsIsAutomatic==1; }
	virtual int Degree() const { return itsDegree; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// the root name
	virtual ARM_CLASS_NAME GetRootName() { return ARM_AMCLONGSTAFFSCHWARTZ; }

private:
	ARM_Regression::RegressionMode itsRegressionMode;
	double itsSpan;
	int itsIsAutomatic;
	int itsDegree;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
