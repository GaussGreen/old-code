/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amc_ls.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/amc_ls.h"

/// gpbase
#include "gpbase/typedef.h"
#include "gpbase/gpmatrixlinalg.h"

/// gpnummethods
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/amcmethod.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"
#include "gpnumlib/regression.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_AMCLongstaffSchwartz
///	Routine: ARM_AMCLongstaffSchwartz
///	Returns:
///	Action : Default constructor
////////////////////////////////////////////////////
ARM_AMCLongstaffSchwartz::ARM_AMCLongstaffSchwartz( size_t ItersNb, ARM_Regression::RegressionMode regressionMode, double span, int isAutomatic, int degree )
{
	itsItersNb = ItersNb;
	itsRegressionMode = regressionMode;
	itsSpan = span;
	itsIsAutomatic = isAutomatic;
	itsDegree = degree;
	SetName(ARM_AMCLONGSTAFFSCHWARTZ);
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCLongstaffSchwartz
///	Routine: Clone,View, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_AMCLongstaffSchwartz::Clone() const
{
	return new ARM_AMCLongstaffSchwartz(*this);
}


string ARM_AMCLongstaffSchwartz::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> Longstaff Schwartz METHOD <====== " << CC_NS(std,endl);
	os << indent << itsItersNb << "Iterations for exercise boundary calculation" << CC_NS(std,endl);
	os << indent << " No Comments yet" << CC_NS(std,endl);

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCLongstaffSchwartz
///	Routine: ComputeExerciseBoundary
///	Returns: ARM_ExerciseBoundaryPtr
///	Action : Does Nothing yet
////////////////////////////////////////////////////

ARM_ExerciseBoundary * ARM_AMCLongstaffSchwartz::ComputeExerciseBoundary(	const ARM_VectorPtr& payoff, 
																			const ARM_VectorPtr& contOpt,
																			const ARM_GP_MatrixPtr& StatesVector )
{
	ARM_LSExerciseBoundary * lseb;
	if (!payoff.IsNull()&&!contOpt.IsNull()&&!StatesVector.IsNull())
	{
		ARM_RegressionPtr regression;
		ARM_GP_VectorPtr tmp = ARM_GP_VectorPtr(new ARM_GP_Vector(*contOpt));
		if (itsRegressionMode == ARM_Regression::LS)
			regression = ARM_RegressionPtr(new ARM_LSRegression(tmp, StatesVector ));
		else if (itsRegressionMode == ARM_Regression::LOESS)
			regression = ARM_RegressionPtr(new ARM_LOESSRegression(tmp, StatesVector, itsSpan ));
		lseb = new ARM_LSExerciseBoundary( regression, IsAutomatic(), Degree() );
	}
	else
	{
		lseb = new ARM_LSExerciseBoundary( ARM_RegressionPtr(NULL), IsAutomatic(), Degree() );
	}

	return lseb;
}
////////////////////////////////////////////////////
///	Class  : ARM_AMCLongstaffSchwartz
///	Routine: ComputeRegression
///	Returns: ARM_GP_VectorPtr B = (X'X)^(-1)X'Y
///	Action : Does Mean Square Regression on a functions space
////////////////////////////////////////////////////

ARM_VectorPtr& ARM_AMCLongstaffSchwartz::ComputeRegression( const ARM_VectorPtrDbleFuncPtrVector& Fcts, 
															const ARM_VectorPtrVector& X, 
															const ARM_VectorPtr& Y )
{
	size_t Rows = Fcts.size();
	size_t Cols = X.size();
	size_t i,j;
	ARM_GP_Matrix Mat( Rows, Cols );

	for ( i = 0 ; i <Rows ; i++ )
		for ( j = 0 ; j< Cols ; j++ )
			Mat(i,j) = (*(Fcts[i]))(X[j]);

	std::vector<double> * A = LeastSquareRegression( Mat, *Y );
	ARM_VectorPtr B = ARM_VectorPtr( new std::vector<double>( *A ) );

	return B;
}


CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
