/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pde1FExplicit.cpp
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 *
 *****  Remarks on the bottom of the file
 */

#include "gpnummethods/pde1Fnumericalschemes.h"
#include "gpinfra/pricingmodel.h"
#include "gpbase/vectormanip.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpinfra/pricingstates.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpnummethods/pdetruncators.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FExplicitNumericalScheme
///	Routine: BuildTransitionMatrixes
///	Returns: void
///	Action : Build Transisition Matrixes for explicit scheme
////////////////////////////////////////////////////
void ARM_PDE1FExplicitNumericalScheme::BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart, double* timeStepEnd )
{
	if (timeStepStart || timeStepEnd)
		ARM_THROW( ERR_INVALID_ARGUMENT, "BuildTransitionMatrixes : timeStepStart/timeStepEnd not implemented " );


	/// ***** Builds matrixes for numerical integration ****
	/// Extract the drift part of the Model
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
    model.EulerLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);

	/// Extract the volatility and correlation part of the Model
    ARM_GP_MatrixPtr vols,d1Vols,correls;
    model.VolatilitiesAndCorrelations(*timeSteps,vols,d1Vols,correls,false/*flatvol*/);

	/// Builds the transition matrixes with the vols and Drifts
	setTimeIndexesToMatrixIndexes(  ARM_GP_VectorPtr( new ARM_GP_Vector( timeSteps->size() ) ) );


	size_t TimeStepSize, TimeIndex;
	size_t StatesSize, diagIndex;
	double DiscretizationStep = getSpaceDiscretizationStep();
	double DiscretizationStepSquare = DiscretizationStep*DiscretizationStep;
	double volSquare = 0;
	double absoluteDrift = 0;
	double relativeDrift = 0;
	double dTdX2 = 0; // dT/dX2
	double dTdX = 0;  // dT/dX
	double dT = 0;
	double state;
	double LastDT;
	double LastVolSquare;

	TimeStepSize = timeSteps->size();
	StatesSize = states->size();

	/// Reinitialising matrixVector
	itsDiagonalElts = ARM_VectorPtrVector(0);
	itsUpperDiagonalElts = ARM_VectorPtrVector(0);
	itsLowerDiagonalElts = ARM_VectorPtrVector(0);
	
	/* For truncation
	itsULLimitConditions = ARM_VectorPtrVector(0);
	itsDRLimitConditions = ARM_VectorPtrVector(0);
	ARM_GP_VectorPtr ULLimitConditions( new ARM_GP_Vector(3,0.) );
	ARM_GP_VectorPtr DRLimitConditions( new ARM_GP_Vector(3,0.) ); */

	LastDT = 0.0;
	LastVolSquare = 0.0;
	size_t MatrixIndex = 0;

	for( TimeIndex = 0 ; TimeIndex < TimeStepSize-1 ; TimeIndex++ )
	{
		/// Compute Vol Square and deltaT
		volSquare = (*vols)(0,TimeIndex);
		volSquare *= volSquare;
		dT = ((*timeSteps)[TimeIndex+1]-(*timeSteps)[TimeIndex])/K_YEAR_LEN;

		/// We build a new matrix only if volsquare or deltaT have changed
		if( (LastVolSquare-volSquare)>K_DOUBLE_TOL || (LastDT-dT)>K_DOUBLE_TOL || (LastVolSquare-volSquare)<-K_DOUBLE_TOL || (LastDT-dT)<-K_DOUBLE_TOL)
		{
			/// Those Vectors will be enqueued
			ARM_GP_VectorPtr DiagElt( new ARM_GP_Vector(StatesSize) );
			ARM_GP_VectorPtr UpperDiagElt( new ARM_GP_Vector(StatesSize-1) );
			ARM_GP_VectorPtr LowerDiagElt( new ARM_GP_Vector(StatesSize-1) );

			absoluteDrift = dT*(*absoluteDrifts)(TimeIndex,0);
			relativeDrift = (*relativeDrifts)(TimeIndex,0);
			dTdX2 = dT/DiscretizationStepSquare;
			dTdX = dT/DiscretizationStep;

			/// Central Part of the matrix
			for( diagIndex = 1 ; diagIndex < StatesSize-1 ; diagIndex ++ )
			{
				state = (*(states->GetNumMethodStates()))(0,diagIndex);
				(*DiagElt)[diagIndex] = 1-volSquare*dTdX2;
				(*LowerDiagElt)[diagIndex-1] = (-0.5*relativeDrift*state+ absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2;
				(*UpperDiagElt)[diagIndex] = (0.5*relativeDrift*state + absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2;
			}
	

			/// Limit conditions
			/// Von Neumann limit conditions
			state = (*(states->GetNumMethodStates()))(0,0);
			(*DiagElt)[0] = 1 + (absoluteDrift- relativeDrift*state)/DiscretizationStep;
			(*UpperDiagElt)[0] = relativeDrift*state/DiscretizationStep;
			state = (*(states->GetNumMethodStates()))(0,StatesSize-1);
			(*DiagElt)[StatesSize-1] = 1 + (absoluteDrift+ relativeDrift*state)/DiscretizationStep ;
			(*LowerDiagElt)[StatesSize-2] = - relativeDrift*state/DiscretizationStep;

			/// Enqueue it all
			itsDiagonalElts.push_back( DiagElt );
			itsUpperDiagonalElts.push_back( UpperDiagElt );
			itsLowerDiagonalElts.push_back( LowerDiagElt );

			if( TimeIndex != 0 )
				MatrixIndex++;

			LastVolSquare = volSquare;
			LastDT = dT;
		}
		setTimeToMatrixIndex(TimeIndex, MatrixIndex);

		LastVolSquare = volSquare;
		LastDT = dT;

	}

/*		For truncation
	if( TruncationSizeVector.size() != 0 )
	{
		for( TimeIndex = 0 ; TimeIndex<TimeStepSize-1 ; ++TimeIndex )
		{
			if( TruncationSizeVector[TimeIndex] != TruncationSizeVector[TimeIndex+1] )
			{
				ULLimitConditions = ARM_GP_VectorPtr( new ARM_GP_Vector(3) );
				DRLimitConditions = ARM_GP_VectorPtr( new ARM_GP_Vector(3) );

				dT = ((*timeSteps)[TimeIndex+1]-(*timeSteps)[TimeIndex])/K_YEAR_LEN;
				dTdX2 = dT/DiscretizationStepSquare;

				absoluteDrift = dT*(*absoluteDrifts)(TimeIndex,0);
				relativeDrift = (*relativeDrifts)(TimeIndex,0);
				state = (*(states->GetNumMethodStates()))(0, (StatesSize-TruncationSizeVector[TimeIndex])/2 );
				/// UpperDiag
				ULLimitConditions->Elt(0) = relativeDrift*state/DiscretizationStep;
				/// Diag
				ULLimitConditions->Elt(1) = 1 + (absoluteDrift- relativeDrift*state)/DiscretizationStep;
				state = (*(states->GetNumMethodStates()))(0, (StatesSize-TruncationSizeVector[TimeIndex])/2+1 );
				/// LowerDiag
				ULLimitConditions->Elt(2) = (-0.5*relativeDrift*state+ absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2;

				state = (*(states->GetNumMethodStates()))(0, (StatesSize+TruncationSizeVector[TimeIndex])/2-2 );
				/// UpperDiag
				DRLimitConditions->Elt(0) = (0.5*relativeDrift*state + absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2;;
				state = (*(states->GetNumMethodStates()))(0, (StatesSize+TruncationSizeVector[TimeIndex])/2-1 );
				/// Diag
				DRLimitConditions->Elt(1) = 1 + (absoluteDrift+ relativeDrift*state)/DiscretizationStep ;
				/// LowerDiag
				DRLimitConditions->Elt(2) = - relativeDrift*state/DiscretizationStep;
			}
			itsULLimitConditions.push_back(ULLimitConditions);
			itsDRLimitConditions.push_back(DRLimitConditions);
		}
	}
	*/

}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FExplicitNumericalScheme
///	Routine: Induct
///	Returns: void
///	Action : Inducts payoffs from one TimeIDx to the next(previous one)
////////////////////////////////////////////////////
void ARM_PDE1FExplicitNumericalScheme::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	/// finds right matrixes
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
	ARM_GP_VectorPtr DiagElt = itsDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr UpperDiagElt = itsUpperDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr LowerDiagElt = itsLowerDiagonalElts[MatrixIndex];
	size_t i;

	/// Truncate Matrixes and states
	// getTruncator()->TruncateTranstionMatrixes( toTimeIdx, UpperDiagElt, LowerDiagElt, DiagElt, itsULLimitConditions[toTimeIdx], itsDRLimitConditions[toTimeIdx]);
	// getTruncator()->TruncatePricingStates( toTimeIdx, states );
	
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();

	/// Update PayoffStates
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);

		ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
		UpdateVectorWithMatrixByProduct( DiagElt, UpperDiagElt, LowerDiagElt, intermediatePayoffs );

		ARM_GP_Vector vec = payoffStates.GetPayoffs();
		UpdateVectorWithMatrixByProduct( DiagElt, UpperDiagElt, LowerDiagElt, vec );
		payoffStates.SetPayoffs( vec );
	}

	/// Update Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_VectorPtr vec = states->GetPayoffVec(i);
		UpdateVectorWithMatrixByProduct( DiagElt, UpperDiagElt, LowerDiagElt, *vec );

		for( int j=0; j<vec->size() ; ++j)
			states->SetPayoff( j, i, vec->Elt(j) );

	}

}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FExplicitNumericalScheme
///	Routine: toString
///	Returns: string
///	Action : As usual for a toString method
////////////////////////////////////////////////////
string ARM_PDE1FExplicitNumericalScheme::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " 1F Explicit Numerical Scheme " << CC_NS(std,endl);
	os << indent << "     Space Discretization : " << getSpaceDiscretizationPointsNb() << " Points" << CC_NS(std,endl);
	os << indent << "     Operture         : " << getStdDevNb() << " Standard Deviations" << CC_NS(std,endl);

#if defined(__GP_STRICT_VALIDATION)
	os << indent << CC_NS(std,endl) << indent << "Used Instant Vols:" << CC_NS(std,endl);
	os << StoredVols->toString();
	os << CC_NS(std,endl);
#endif

	return os.str();
}

CC_END_NAMESPACE()

/*
* The implementation has been made as follows
*  moving from one state to another relies on multiplication by a tridiagonal matrix
*  tridiagonal matrixes are stored in three vectors, of size nbstates for itsDiagonal, and 
*  nbstates-1 for the upper and lower lines. 
* 
*  Basically, 
*						DiagonalElt(0)			  itsUpperDiagonalElt(0)   0
*  TransitionMatrix =   itsLowerDiagonalElt(0)    DiagonalElt(1)			itsUpperDiagonalElt(1)	0
*								0					itsLowerDiagonalElt(1)	DiagonalElt(2)			itsUpperDiagonalElt(2)	0
*								...							...					...							...					...
*								...							...					...							...					...			itsUpperDiagonalElt(matrixSize-2)
*																									itsLowerDiagonalElt(matrixSize-2)		DiagonalElt(matrixSize-1)
*  This explains how multiplications are performed (in UpdateVectorWithMatrixByProduct).
*
*	To get those elements from date 
*   size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
*	itsDiagonalElts[MatrixIndex]
*	itsUpperDiagonalElts[MatrixIndex]
*	itsLowerDiagonalElts[MatrixIndex]
*/