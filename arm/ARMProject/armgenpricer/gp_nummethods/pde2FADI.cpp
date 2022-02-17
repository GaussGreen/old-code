/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pde2FExplicit.cpp
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 *
 *****  Remarks on the bottom of the file
 */

#include "gpnummethods/pde2Fnumericalschemes.h"
#include "gpinfra/pricingmodel.h"
#include "gpbase/vectormanip.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpinfra/pricingstates.h"
#include "gpbase/gpmatrixlinalg.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: toString
///	Returns: string
///	Action : toString
////////////////////////////////////////////////////
string ARM_PDE2FADINumericalScheme::toString(const string& indent, const string&nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " 2F Crank Nicholson Numerical Scheme " << CC_NS(std,endl);
	os << indent << "     Space Discretization : " << getSpaceDiscretizationPointsNb() << " Points" << CC_NS(std,endl);
	os << indent << "     Operture         : " << getStdDevNb() << " Standard Deviations" << CC_NS(std,endl);

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: ~ARM_PDE2FADINumericalScheme
///	Returns: nothing
///	Action : destructor
////////////////////////////////////////////////////
ARM_PDE2FADINumericalScheme::~ARM_PDE2FADINumericalScheme()
{
	size_t diagIndex;
	size_t matrixSize = getSpaceDiscretizationPointsNb();

	if( vector1 && vector2 )
	{
		for( diagIndex = 0 ; diagIndex < matrixSize ; diagIndex++)
		{
			free( vector1[diagIndex] );
			free( vector2[diagIndex] );
		}
		free( vector1 );
		free( vector2 );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: BuildTransitionMatrixes
///	Returns: viod
///	Action : Builds Transition Matrixes used in Induct
////////////////////////////////////////////////////
void ARM_PDE2FADINumericalScheme::BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart, double* timeStepEnd)
{
	if (timeStepStart || timeStepEnd)
		ARM_THROW( ERR_INVALID_ARGUMENT, "BuildTransitionMatrixes : timeStepStart/timeStepEnd not implemented " );

	/// ***** Builds matrixes for numerical integration ****
	/// Extract the drift part of the Model
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
    model.EulerLocalDrifts(timeSteps->GetValues(),relativeDrifts,absoluteDrifts);

	/// Extract the volatility and correlation part of the Model
    ARM_GP_MatrixPtr vols,d1Vols,correls;
    model.VolatilitiesAndCorrelations(timeSteps->GetValues(),vols,d1Vols,correls,false/*flatvol*/);

	/// Builds the transition matrixes with the vols and Drifts
	setTimeIndexesToMatrixIndexes(  ARM_GP_VectorPtr( new ARM_GP_Vector( timeSteps->size() ) ) );


	size_t TimeStepSize, TimeIndex;
	size_t StatesSize, diagIndex;
	double DiscretizationStepX = getSpaceDiscretizationStepX();
	double DiscretizationStepXSquare = DiscretizationStepX*DiscretizationStepX;
	double DiscretizationStepY = getSpaceDiscretizationStepY();
	double DiscretizationStepYSquare = DiscretizationStepY*DiscretizationStepY;
	double volSquareX = 0;
	double volSquareY = 0;
	double absoluteDriftX = 0;
	double relativeDriftX = 0;
	double absoluteDriftY = 0;
	double relativeDriftY = 0;
	double dTdX2 = 0;
	double dTdX = 0;
	double dTdY2 = 0;
	double dTdY = 0;
	double dT = 0;
	double OneOndX2 = 0;
	double OneOndY2 = 0;
	double OneOndT = 0;
	double stateX,stateY;
	double LastDT;
	double LastVolSquareX = 0;
	double LastVolSquareY = 0;
	size_t matrixSize;
	size_t MatrixIndex = 0;

	TimeStepSize = timeSteps->size();
	StatesSize = states->size();
	matrixSize = getSpaceDiscretizationPointsNb();

	itsTimeSteps = timeSteps;

	/// Initialises vectors that contains the tridiag matrixes
	itsRBYUOperators =ARM_VectorPtrVector(0);
	itsLAXUOperators = ARM_VectorPtrVector(0);
	itsLBYUOperators = ARM_VectorPtrVector(0);

	itsRBYDOperators =ARM_VectorPtrVector(0);
	itsLAXDOperators = ARM_VectorPtrVector(0);
	itsLBYDOperators = ARM_VectorPtrVector(0);

	itsRBYDiagOperators =ARM_VectorPtrVector(0);
	itsLAXDiagOperators = ARM_VectorPtrVector(0);
	itsLBYDiagOperators = ARM_VectorPtrVector(0);

	itsCorrelFactors = ARM_GP_VectorPtr( new ARM_GP_Vector(0) );

	/// building intermediate vectors in the form of matrixes
	vector1 = (double**) malloc( matrixSize* sizeof(double* ) );
	vector2 = (double**) malloc( matrixSize* sizeof(double* ) );

	for( diagIndex = 0 ; diagIndex < matrixSize ; diagIndex++)
	{
		vector1[diagIndex] = (double*) malloc( matrixSize*sizeof( double ) );
		vector2[diagIndex] = (double*) malloc( matrixSize*sizeof( double ) );
	}

	for( TimeIndex = 0 ; TimeIndex < TimeStepSize-1 ; ++TimeIndex )
	{
		volSquareX = (*vols)(0,TimeIndex);
		volSquareX *= volSquareX;
		volSquareY = (*vols)(1,TimeIndex);
		volSquareY *= volSquareY;
		dT = ((*timeSteps)[TimeIndex+1]-(*timeSteps)[TimeIndex])/K_YEAR_LEN;
		OneOndT = 1/dT;
		OneOndX2 = 1/DiscretizationStepXSquare;
		OneOndY2 = 1/DiscretizationStepYSquare;

		if( (LastVolSquareX-volSquareX)>K_DOUBLE_TOL || (LastDT-dT)>K_DOUBLE_TOL || (LastVolSquareX-volSquareX)<-K_DOUBLE_TOL || (LastDT-dT)<-K_DOUBLE_TOL)
		{
			/// Tridiagonal matrixes are enqueued in the form of three vectors. 
			ARM_GP_VectorPtr RBYUOperator( new ARM_GP_Vector( matrixSize-1 ) );
			ARM_GP_VectorPtr LAXUOperator( new ARM_GP_Vector( matrixSize-1 ) );
			ARM_GP_VectorPtr LBYUOperator( new ARM_GP_Vector( matrixSize-1 ) );

			ARM_GP_VectorPtr RBYDOperator( new ARM_GP_Vector( matrixSize-1 ) );
			ARM_GP_VectorPtr LAXDOperator( new ARM_GP_Vector( matrixSize-1 ) );
			ARM_GP_VectorPtr LBYDOperator( new ARM_GP_Vector( matrixSize-1 ) );

			ARM_GP_VectorPtr RBYDiagOperator( new ARM_GP_Vector( matrixSize ) );
			ARM_GP_VectorPtr LAXDiagOperator( new ARM_GP_Vector( matrixSize ) );
			ARM_GP_VectorPtr LBYDiagOperator( new ARM_GP_Vector( matrixSize ) );

			double crossElt=0;

			/// Precomputing stuff, for faster matrix building
			/// because Euler Local drifts are integrated in time, we divide by dT
			absoluteDriftX = 0.25*(*absoluteDrifts)(TimeIndex,0)*OneOndT/DiscretizationStepX;
			relativeDriftX = 0.25*(*relativeDrifts)(TimeIndex,0)*OneOndT/DiscretizationStepX;
			absoluteDriftY = 0.25*(*absoluteDrifts)(TimeIndex,1)*OneOndT/DiscretizationStepY;
			relativeDriftY = 0.25*(*relativeDrifts)(TimeIndex,1)*OneOndT/DiscretizationStepY;

			/// Diagonal terms of the matrix
			double RBYDiagTerm = OneOndT-0.5*volSquareY*OneOndY2;
			double LBYDiagTerm = OneOndT+0.5*volSquareY*OneOndY2;
			double LAXDiagTerm = OneOndT+0.5*volSquareX*OneOndX2;

			/// Constant term of their upper diagonals
			double constantTermLX = -0.25*volSquareX*OneOndX2;
			double constantTermLY = -0.25*volSquareY*OneOndY2;
			double constantTermRY = -constantTermLY;

			double relDriftX = 0;
			double relDriftY = 0;

			/// Central part of the matrixes
			for( diagIndex = 1 ; diagIndex < matrixSize-1 ; ++diagIndex )
			{
				stateX = (*(states->GetNumMethodStates()))(0,diagIndex*matrixSize);
				stateY = (*(states->GetNumMethodStates()))(1,diagIndex);

				relDriftX = stateX*relativeDriftX;
				relDriftY = stateY*relativeDriftY;

				RBYUOperator->Elt(diagIndex) = constantTermRY+relDriftY+absoluteDriftY;
				RBYDOperator->Elt(diagIndex-1) = constantTermRY-relDriftY-absoluteDriftY;
				RBYDiagOperator->Elt(diagIndex) = RBYDiagTerm;

				LBYUOperator->Elt(diagIndex) = constantTermLY-relDriftY-absoluteDriftY;
				LBYDOperator->Elt(diagIndex-1) = constantTermLY+relDriftY+absoluteDriftY;
				LBYDiagOperator->Elt(diagIndex) = LBYDiagTerm;

				LAXUOperator->Elt(diagIndex) = constantTermLX-relDriftX-absoluteDriftX;
				LAXDOperator->Elt(diagIndex-1) = constantTermLX+relDriftX+absoluteDriftX;
				LAXDiagOperator->Elt(diagIndex) = LAXDiagTerm;
			}

			/// Limit conditions, Xmax, Ymax
			stateX = (*(states->GetNumMethodStates()))(0,diagIndex*matrixSize);
			stateY = (*(states->GetNumMethodStates()))(1,diagIndex);

			absoluteDriftX*=2;
			relativeDriftX*=2;
			absoluteDriftY*=2;
			relativeDriftY*=2;

			relDriftX = stateX*relativeDriftX;
			relDriftY = stateY*relativeDriftY;

			RBYDOperator->Elt(diagIndex-1) = -relDriftY-absoluteDriftY;
			RBYDiagOperator->Elt(diagIndex) = OneOndT+relDriftY+absoluteDriftY;

			LBYDOperator->Elt(diagIndex-1) = relDriftY+absoluteDriftY;
			LBYDiagOperator->Elt(diagIndex) = OneOndT-relDriftY-absoluteDriftY;

			LAXDOperator->Elt(diagIndex-1) = relDriftX+absoluteDriftX;
			LAXDiagOperator->Elt(diagIndex) = OneOndT-relDriftX-absoluteDriftX;

			/// Limit conditions Xmin, Ymin
			stateX = (*(states->GetNumMethodStates()))(0,0);
			stateY = (*(states->GetNumMethodStates()))(1,0);

			relDriftX = stateX*relativeDriftX;
			relDriftY = stateY*relativeDriftY;

			RBYUOperator->Elt(0) = relDriftY+absoluteDriftY;
			RBYDiagOperator->Elt(0) = OneOndT-relDriftY-absoluteDriftY;

			LBYUOperator->Elt(0) = -relDriftY-absoluteDriftY;
			LBYDiagOperator->Elt(0) = OneOndT+relDriftY+absoluteDriftY;

			LAXUOperator->Elt(0) = -relDriftX-absoluteDriftX;
			LAXDiagOperator->Elt(0) = OneOndT+relDriftX+absoluteDriftX;

			
			if( TimeIndex != 0 )
				MatrixIndex++;

			/// Enqueueing it all
			itsRBYUOperators.push_back(	RBYUOperator );
			itsLAXUOperators.push_back(	LAXUOperator );
			itsLBYUOperators.push_back( LBYUOperator );

			itsRBYDOperators.push_back(	RBYDOperator );
			itsLAXDOperators.push_back( LAXDOperator );
			itsLBYDOperators.push_back( LBYDOperator );

			itsRBYDiagOperators.push_back( RBYDiagOperator );
			itsLAXDiagOperators.push_back( LAXDiagOperator );
			itsLBYDiagOperators.push_back( LBYDiagOperator );

			itsCorrelFactors->push_back( correls->Elt(0,TimeIndex)*sqrt( volSquareX*volSquareY )/(16.0*DiscretizationStepX*DiscretizationStepY ) );

			LastVolSquareX = volSquareX;
			LastDT = dT;
		}

		setTimeToMatrixIndex(TimeIndex, MatrixIndex);

		LastVolSquareX = volSquareX;
		LastDT = dT;
	}

	nextPayoffStates = ARM_MatrixPtrVector(0);
	nextPayoffs = ARM_VectorPtrVector(0);

}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: Induct
///	Returns: viod
///	Action : Updates pricing states from one date to another
////////////////////////////////////////////////////
void ARM_PDE2FADINumericalScheme::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	/// gets matrixes needed for Induction
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
	ARM_GP_VectorPtr LAXU = itsLAXUOperators[MatrixIndex];
	ARM_GP_VectorPtr LAXD = itsLAXDOperators[MatrixIndex];
	ARM_GP_VectorPtr LAXDiag = itsLAXDiagOperators[MatrixIndex];
	ARM_GP_VectorPtr LBYU = itsLBYUOperators[MatrixIndex];
	ARM_GP_VectorPtr LBYD = itsLBYDOperators[MatrixIndex];
	ARM_GP_VectorPtr LBYDiag = itsLBYDiagOperators[MatrixIndex];

	size_t i;
	
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();

	/// Precomputation in the Thomas algorithm
	PrecomputeMatrixesForInversion( LAXD, LAXDiag, LAXU, itsALNormalizationTerms, itsALOtherCoeffs );
	PrecomputeMatrixesForInversion( LBYD, LBYDiag, LBYU, itsBLNormalizationTerms, itsBLOtherCoeffs );

	/// Precomputation of the crossed part in the predictor/corrector stuff
	PrecomputeCrossedTerms( states );

	/// PayoffStatesSize
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);

		ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
		UpdatePayoffs( intermediatePayoffs, i, toTimeIdx );

		std::vector<double> vec = payoffStates.GetPayoffs();
		UpdatePayoffs( vec, diffPayoffStatesPayoffs[i]->GetValues(), toTimeIdx  );
		payoffStates.SetPayoffs( vec );
	}

	/// Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_VectorPtr vec = states->GetPayoffVec(i);
		UpdatePayoffs(*vec, diffPayoffs[i]->GetValues(), toTimeIdx  );

		for( int j=0; j<vec->size() ; ++j)
			states->SetPayoff( j, i,(*vec)[j] );

	}

	if( toTimeIdx == 0)
	{
		if( vector1 && vector2 )
		{
			int matrixSize = LAXDiag->size();
			for( int diagIndex = 0 ; diagIndex < matrixSize ; diagIndex++)
			{
				free( vector1[diagIndex] );
				free( vector2[diagIndex] );
			}
			free( vector1 );
			free( vector2 );

			vector1 = NULL;
			vector2 = NULL;
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: UpdatePayoffs
///	Returns: viod
///	Action : Updates payoffs from one date to another
////////////////////////////////////////////////////
void ARM_PDE2FADINumericalScheme::UpdatePayoffs(std::vector<double>& vec,const std::vector<double>& PredictorCorrector, int toTimeIdx)
{
	size_t StatesSize = vec.size();
	int i,j;
	size_t matrixSize = getSpaceDiscretizationPointsNb();
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
	ARM_GP_VectorPtr CorrelTerm( new ARM_GP_Vector(StatesSize) );
	std::vector<double>::iterator iter, iterEnd;
	std::vector<double>::iterator DiagIter, UIter, LIter;
	double *matrixPointer1, *matrixPointer2;

	iterEnd = vec.end();
	iter = vec.begin();

	/// get matrixes to update payoffs
	ARM_GP_VectorPtr RBYDiagOperator = itsRBYDiagOperators[MatrixIndex];
	ARM_GP_VectorPtr RBYUOperator = itsRBYUOperators[MatrixIndex];
	ARM_GP_VectorPtr RBYLOperator = itsRBYDOperators[MatrixIndex];

	ARM_GP_VectorPtr LAXDOperator = itsLAXDOperators[MatrixIndex];
	ARM_GP_VectorPtr LBYDOperator = itsLBYDOperators[MatrixIndex];
	double correlTerm = itsCorrelFactors->Elt(MatrixIndex);

	for( i = 0 ; i< matrixSize ; ++i )
	{
		matrixPointer1 = vector1[i];
		DiagIter = RBYDiagOperator->begin();
		UIter = RBYUOperator->begin();
		LIter = RBYLOperator->begin();

		(*matrixPointer1) = (*iter)*(*DiagIter)
			+(*(iter+1))*(*UIter);

		iter++;
		DiagIter++;
		UIter++;
		matrixPointer1++;

		for( j=1 ; j < matrixSize-1 ; ++j )
		{
			(*matrixPointer1) = (*iter)*(*DiagIter)
				+(*(iter-1))*(*LIter)
				+(*(iter+1))*(*UIter);

			matrixPointer1++;
			iter++;
			UIter++;
			LIter++;
			DiagIter++;
		}

		(*matrixPointer1) = (*iter)*(*DiagIter)
			+(*(iter-1))*(*LIter);

		DiagIter++;
		LIter++;
		iter++;

#if defined(__GP_STRICT_VALIDATION)
		if( DiagIter != RBYDiagOperator->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "BYDiagIter not at then end!" );

		if( LIter!= RBYLOperator->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "LIter not at then end!" );

		if( UIter!= RBYUOperator->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "UIter not at then end!" );

		if( (*matrixPointer1) != vector1[i][matrixSize-1] )
			ARM_THROW( ERR_INVALID_ARGUMENT, "wrong value vector[i][end]" );

		if( iter != (vec.begin()+(i+1)*matrixSize) )
			ARM_THROW( ERR_INVALID_ARGUMENT, "iter ill placed!" );
#endif

	}

	/// Now vector1 contains BR*vec

	i=0;
	j=0;

	matrixPointer2 = vector2[i];

	(*matrixPointer2) = vector1[j][i]+
				4*correlTerm*( PredictorCorrector[(i+1)*matrixSize+j+1]+PredictorCorrector[(i)*matrixSize+j]
				- PredictorCorrector[(i)*matrixSize+j+1] - PredictorCorrector[(i+1)*matrixSize+j] );
	matrixPointer2++;

	for( j=1 ; j<matrixSize-1 ; ++j )
	{
		(*matrixPointer2) = vector1[j][i]+
				2*correlTerm*( PredictorCorrector[(i+1)*matrixSize+j+1]+PredictorCorrector[(i)*matrixSize+j-1]
				- PredictorCorrector[(i)*matrixSize+j+1] - PredictorCorrector[(i+1)*matrixSize+j-1] );
		matrixPointer2++;
	}

	(*matrixPointer2) = vector1[j][i]+
			4*correlTerm*( PredictorCorrector[(i+1)*matrixSize+j]+PredictorCorrector[(i)*matrixSize+j-1]
			- PredictorCorrector[(i)*matrixSize+j] - PredictorCorrector[(i+1)*matrixSize+j-1] );
	matrixPointer2++;

	for( i=1 ; i< matrixSize-1 ; ++i )
	{
		j=0;
		matrixPointer2 = vector2[i];

		(*matrixPointer2) = vector1[j][i]+
				2*correlTerm*( PredictorCorrector[(i+1)*matrixSize+j+1]+PredictorCorrector[(i-1)*matrixSize+j]
				- PredictorCorrector[(i-1)*matrixSize+j+1] - PredictorCorrector[(i+1)*matrixSize+j] );
		matrixPointer2++;

 		for( j=1 ; j<matrixSize-1 ; ++j )
		{
			(*matrixPointer2) = vector1[j][i]+
				correlTerm*( PredictorCorrector[(i+1)*matrixSize+j+1]+PredictorCorrector[(i-1)*matrixSize+j-1]
				- PredictorCorrector[(i-1)*matrixSize+j+1] - PredictorCorrector[(i+1)*matrixSize+j-1] );
			matrixPointer2++;
		}

		(*matrixPointer2) = vector1[j][i]+
				2*correlTerm*( PredictorCorrector[(i+1)*matrixSize+j]+PredictorCorrector[(i-1)*matrixSize+j-1]
				- PredictorCorrector[(i-1)*matrixSize+j] - PredictorCorrector[(i+1)*matrixSize+j-1] );
		matrixPointer2++;

	}

	j=0;
	i=matrixSize-1;

	matrixPointer2 = vector2[i];

	(*matrixPointer2) = vector1[j][i]+
		4*correlTerm*( PredictorCorrector[(i)*matrixSize+j+1]+PredictorCorrector[(i-1)*matrixSize+j]
		- PredictorCorrector[(i-1)*matrixSize+j+1] - PredictorCorrector[(i)*matrixSize+j] );
	matrixPointer2++;

	for( j=1 ; j<matrixSize-1 ; ++j )
	{
		(*matrixPointer2) = vector1[j][i]+
			2*correlTerm*( PredictorCorrector[(i)*matrixSize+j+1]+PredictorCorrector[(i-1)*matrixSize+j-1]
			- PredictorCorrector[(i-1)*matrixSize+j+1] - PredictorCorrector[(i)*matrixSize+j-1] );
		matrixPointer2++;
	}

	(*matrixPointer2) = vector1[j][i]+
			4*correlTerm*( PredictorCorrector[(i)*matrixSize+j]+PredictorCorrector[(i-1)*matrixSize+j-1]
			- PredictorCorrector[(i-1)*matrixSize+j] - PredictorCorrector[(i)*matrixSize+j-1] );
	matrixPointer2++;


	/// Now vector2 Contains CrossedTerm + Transposed Vector1. Vector2 is organised X;
	/// Computing A^-1*vector1

	for( i=0 ; i< matrixSize ; ++i )
	{
		matrixPointer2 = vector2[i];	

		DiagIter = itsALNormalizationTerms->begin();
		(*matrixPointer2) /= (*DiagIter);
		LIter = LAXDOperator->begin();
		DiagIter++;
		matrixPointer2++;
		
		for( j=1 ; j<matrixSize; ++j )
		{
			(*matrixPointer2) = ((*matrixPointer2)-(*LIter)*(*(matrixPointer2-1)))/(*DiagIter);
			matrixPointer2++;
			DiagIter++;
			LIter++;
		}

#if defined(__GP_STRICT_VALIDATION)
		if( DiagIter != itsALNormalizationTerms->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "DiagIter not at the end!" );

		if( LIter != LAXDOperator->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "LIter not at the end!" );
#endif

		matrixPointer2--;
		matrixPointer2--;

#if defined(__GP_STRICT_VALIDATION)
		if( (*matrixPointer2) != vector2[i][matrixSize-2] )
			ARM_THROW( ERR_INVALID_ARGUMENT, "matrixPointer2 ill positioned!" );
#endif

		DiagIter = itsALOtherCoeffs->end()-1;

		for( j=matrixSize-2 ; j>0; --j )
		{
			(*matrixPointer2) = (*matrixPointer2)-(*DiagIter)*(*(matrixPointer2+1));
			matrixPointer2--;
			DiagIter--;
		}

		(*matrixPointer2) = (*matrixPointer2)-(*DiagIter)*(*(matrixPointer2+1));

#if defined(__GP_STRICT_VALIDATION)
		if( DiagIter != itsALOtherCoeffs->begin() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "DiagIter ill positioned!" );

		if( (*matrixPointer2) != vector2[i][0] )
			ARM_THROW( ERR_INVALID_ARGUMENT, "matrixPointer2 ill positioned!" );
#endif

	}


	double TwoOnDT = 2.0/(((*itsTimeSteps)[toTimeIdx+1]-(*itsTimeSteps)[toTimeIdx])/K_YEAR_LEN );

	for( i=0 ; i< matrixSize ; ++i )
	{
		matrixPointer1 = vector1[i];

		for( j=0 ; j<matrixSize ; ++j )
		{
			(*matrixPointer1) = TwoOnDT*(vector2[j][i])-(*matrixPointer1);
			matrixPointer1++;
		}
	}

	/// Now vector1 contains BL*Solution. It is organised Y;

	iter = iterEnd-1;

	for( i=matrixSize-1 ; i>=0 ; --i )
	{
		matrixPointer1 = vector1[i];

		DiagIter = itsBLNormalizationTerms->begin();
		(*matrixPointer1) /= (*DiagIter);
		LIter = LBYDOperator->begin();
		DiagIter++;
		matrixPointer1++;
		
		for( j=1 ; j<matrixSize; ++j )
		{
			(*matrixPointer1) = ((*matrixPointer1)-(*LIter)*(*(matrixPointer1-1)))/(*DiagIter);
			matrixPointer1++;
			DiagIter++;
			LIter++;
		}

#if defined(__GP_STRICT_VALIDATION)
		if( DiagIter != itsBLNormalizationTerms->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "DiagIter not at the end!" );

		if( LIter != LBYDOperator->end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "LIter not at the end!" );
#endif

		matrixPointer1--;

		(*iter) = (*matrixPointer1);
		iter--;

		matrixPointer1--;

#if defined(__GP_STRICT_VALIDATION)
		if( (*matrixPointer1) != vector1[i][matrixSize-2] )
			ARM_THROW( ERR_INVALID_ARGUMENT, "matrixPointer1 ill positioned!" );
#endif

		DiagIter = itsBLOtherCoeffs->end()-1;

		for( j=matrixSize-2 ; j>0; --j )
		{
			(*iter) = (*matrixPointer1)-(*DiagIter)*(*(iter+1));
			matrixPointer1--;
			DiagIter--;
			iter--;
		}

		(*iter) = (*matrixPointer1)-(*DiagIter)*(*(iter+1));
		iter--;

#if defined(__GP_STRICT_VALIDATION)
		if( DiagIter != itsBLOtherCoeffs->begin() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "DiagIter ill positioned!" );

		if( (*matrixPointer1) != vector1[i][0] )
			ARM_THROW( ERR_INVALID_ARGUMENT, "matrixPointer1 ill positioned!" );
#endif

	}

#if defined(__GP_STRICT_VALIDATION)
	if( iter+1 != vec.begin() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Iter ill positioned!" );
#endif

}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: UpdatePayoffs
///	Returns: viod
///	Action : Same thing as above. Uses UpdatePayoffs from above for 
///  each line of the matrix
////////////////////////////////////////////////////
void ARM_PDE2FADINumericalScheme::UpdatePayoffs(ARM_GP_Matrix& vec, int payoffStateIdx, int toTimeIdx)
{
	std::vector<double> vec2(vec.cols());

	for( size_t j=0 ; j<vec.rows() ; ++j )
	{
		for( size_t i=0 ; i<vec.cols() ; ++i)
			vec2[i] = vec(j,i);

		UpdatePayoffs( vec2, (diffPayoffStates[payoffStateIdx])->GetRow(j)->GetValues(), toTimeIdx );
	
		for(int i=0 ; i<vec.cols() ; ++i)
			vec(j,i) = vec2[i];
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: PrecomputeMatrixesForInversion
///	Returns: viod
///	Action : Precomputes inverses according to the Thomas Alogorithm and stores them
////////////////////////////////////////////////////
void ARM_PDE2FADINumericalScheme::PrecomputeMatrixesForInversion( const ARM_GP_VectorPtr& LowerTerm, const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, 
		ARM_GP_VectorPtr& NormalizationTerms, ARM_GP_VectorPtr& OtherCoeffs)
{
	size_t MatrixSize = DiagTerm->size();

	if( NormalizationTerms == ARM_GP_VectorPtr(NULL) )
		NormalizationTerms = ARM_GP_VectorPtr( new ARM_GP_Vector( MatrixSize ) );

	if( OtherCoeffs == ARM_GP_VectorPtr(NULL) )
		OtherCoeffs = ARM_GP_VectorPtr( new ARM_GP_Vector( MatrixSize-1 ) );

	size_t i;
	

	double Normalisation = DiagTerm->Elt(0);
	double OtherCoeff = UpperTerm->Elt(0) / Normalisation;

	NormalizationTerms->Elt(0) = Normalisation;
	OtherCoeffs->Elt(0) = OtherCoeff;

	for( i=1 ; i<MatrixSize-1 ; ++i )
	{
		Normalisation = DiagTerm->Elt(i)-LowerTerm->Elt(i-1)*OtherCoeff;
		OtherCoeff = UpperTerm->Elt(i)/Normalisation;
		NormalizationTerms->Elt(i) = Normalisation;
		OtherCoeffs->Elt(i) = OtherCoeff;
	}

	Normalisation = DiagTerm->Elt(MatrixSize-1)-LowerTerm->Elt(MatrixSize-2)*OtherCoeff;
	NormalizationTerms->Elt(MatrixSize-1) = Normalisation;
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FADINumericalScheme
///	Routine: PrecomputeCrossedTerms
///	Returns: viod
///	Action : Precomputes the cross matrixes for payoffs, intermediatepayoffs
///   This is the crossterm in the predictor/corrector stuff
////////////////////////////////////////////////////
void ARM_PDE2FADINumericalScheme::PrecomputeCrossedTerms( const ARM_PricingStatesPtr& states )
{
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();
	size_t i;
	std::vector<double>::iterator storedIter, otherIter, endIter, diffIter;
	ARM_GP_Matrix::iterator mStoredIter, mEndIter, mDiffIter;
	ARM_GP_Matrix::const_iterator mOtherIter;

	//// If storage matrixes do not exist ...
	if( nextPayoffs.size() == 0 )
	{
		for( i=0 ; i<PayoffsSize ; ++i)
		{
			ARM_GP_VectorPtr vec = ARM_GP_VectorPtr(new ARM_GP_Vector(*states->GetPayoffVec(i)));
			nextPayoffs.push_back(vec);
			diffPayoffs.push_back(vec);
		}
	}

	if( nextPayoffStates.size() == 0 )
	{
		for( i=0 ; i<PayoffStatesSize ; ++i)
		{
			ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);
			const ARM_GP_Matrix& intermediatePayoffs = payoffStates.GetIntermediatePayoffs();

			nextPayoffStates.push_back( ARM_GP_MatrixPtr( static_cast<ARM_GP_Matrix*>( intermediatePayoffs.Clone() ) ) );
			diffPayoffStates.push_back( ARM_GP_MatrixPtr( static_cast<ARM_GP_Matrix*>( intermediatePayoffs.Clone() ) ) );

			std::vector<double>& vec = const_cast<std::vector<double>&> (payoffStates.GetPayoffs());
			nextPayoffStatesPayoffs.push_back(ARM_GP_VectorPtr(new ARM_GP_Vector(vec)));
			diffPayoffStatesPayoffs.push_back(ARM_GP_VectorPtr(new ARM_GP_Vector(vec)));
		}
	}

	size_t PayoffStatesRows = nextPayoffStates[0]->rows(), PayoffStatesCols = nextPayoffStates[0]->cols();

	/// If storage matrixes do not have the right size.
	if( states->GetPayoffStates(0).GetIntermediatePayoffs().rows() != PayoffStatesRows )
	{
		for( i=0 ; i<PayoffStatesSize ; ++i)
		{
			nextPayoffStates[i]->resize( PayoffStatesRows+1, PayoffStatesCols);
			diffPayoffStates[i]->resize( PayoffStatesRows+1, PayoffStatesCols);
			ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);
			const ARM_GP_Matrix& intermediatePayoffs = payoffStates.GetIntermediatePayoffs();
			ARM_GP_Matrix::iterator PayoffStatesPayoffsIterator = nextPayoffStatesPayoffs[i]->begin();

			ARM_GP_Matrix::const_iterator PayoffStatesIterator = intermediatePayoffs.begin();
			PayoffStatesIterator += PayoffStatesRows*PayoffStatesCols;
			ARM_GP_Matrix::const_iterator PayoffStatesIteratorEnd = intermediatePayoffs.end();
			ARM_GP_Matrix::iterator StoredPayoffsStatesIterator = nextPayoffStates[i]->begin() + PayoffStatesRows*PayoffStatesCols;

			for( ; PayoffStatesIterator != PayoffStatesIteratorEnd ; ++PayoffStatesIterator, ++StoredPayoffsStatesIterator, ++PayoffStatesPayoffsIterator)
			{
				(*StoredPayoffsStatesIterator) = (*PayoffStatesIterator);
				(*PayoffStatesPayoffsIterator) += (*PayoffStatesIterator);
			}

		}
	}

	/// Filling payoffstates/payofsf
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_VectorPtr vec = states->GetPayoffVec(i);
		ARM_VectorPtr nextvec = nextPayoffs[i];
		otherIter = vec->begin();
		storedIter = nextvec->begin();
		endIter = nextvec->end();
		diffIter = diffPayoffs[i]->begin();

		for( ; storedIter != endIter ; storedIter++, otherIter++, diffIter++ )
		{
			(*diffIter) = 3*(*otherIter) - (*storedIter);
			(*storedIter) = (*otherIter);
		}
	}

	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);
		const ARM_GP_Matrix& intermediatePayoffs = payoffStates.GetIntermediatePayoffs();
		mStoredIter = nextPayoffStates[i]->begin();
		mOtherIter = intermediatePayoffs.begin();
		mEndIter = nextPayoffStates[i]->end();
		mDiffIter = diffPayoffStates[i]->begin();

		for( ; mStoredIter != mEndIter ; mStoredIter++, mOtherIter++, mDiffIter++ )
		{
			(*mDiffIter) = 3*(*mOtherIter) - (*mStoredIter);
			(*mStoredIter) = (*mOtherIter);
		}


		std::vector<double> vec = payoffStates.GetPayoffs();
		ARM_VectorPtr nextvec = nextPayoffStatesPayoffs[i];
		otherIter = vec.begin();
		storedIter = nextvec->begin();
		endIter = nextvec->end();
		diffIter = diffPayoffStatesPayoffs[i]->begin();

		for( ; storedIter != endIter ; storedIter++, otherIter++, diffIter++ )
		{
			(*diffIter) = 3*(*otherIter) - (*storedIter);
			(*storedIter) = (*otherIter);
		}
	}
}


CC_END_NAMESPACE()

/*
* The implementation has been made as follows
*
* The numerical scheme works this way: 
*  it is an ADI with crossed terms using predictor/corrector techniques. 
*
*  as usual, we want to move from C^n+1 to C^n
*   U = V^n + 1/2*(V^n - V^n+1)
*   (1/dT I - 1/2 A)V = (1/dT I + 1/2 B)V^n + 1/2 D U
*	(1/dT I - 1/2 B) C^n-1 = (1/dT I + 1/2 A)V + 1/2 D U
*
*  We write: 
*  1/dT I - 1/2 A = AL
*  1/dT I + 1/2 B = BR
*  1/dT I - 1/2 B = BL
*  Note that (1/dT I + 1/2 A)V = 2/dT V - (1/dT I - 1/2 A)V
*
*   U = V^n + 1/2*(V^n - V^n+1) is computed and V^n+1 is stored for every payoffstates, etc
*  in PrecomputeCrossedTerms. 
*  The difference are stored ni diffPayoffStates, diffPayoffStatesPayoffs, diffPayoffs. 
*  The V^n+1 (previous values) are stored in nextPayoffStates, nextPayoffStatesPayoffs, nextPayoffs. 
*   in fact, 3*V^n-V^n+1 is stored. 
* 
*  Then we compute BR*V^n+1 ; it is stored in vector1[i][j]. where i corresponds to an X index, and j to Y. 
*   Then we compute vector2[i][j] = vector1[j][i] + 1/2 D U
*    vector2[i][j] : i -> Y, j-> X. 
*   then vector1[i][j] = AL^-1 * vector2[i][j]
*     vector1[i][j] = 2/dT vector1[i][j] - vector2[j][i]  (according to the previous remark). 
*    V^n-1 = BL^-1 vector1[i][j]. 
*
*     BR
*	itsRBYUOperators;
*	itsRBYDOperators;
*	itsRBYDiagOperators;
*
*     AL
*	itsLAXUOperators;
*	itsLAXDOperators;
*	itsLAXDiagOperators;

      BL
*	itsLBYUOperators;
*	itsLBYDOperators;
*	itsLBYDiagOperators;
*
*
*
     Coeffs for inversion are stored in: 
*	itsALNormalizationTerms
*	itsALOtherCoeffs;
*	itsBLNormalizationTerms
*	itsBLOtherCoeffs;
 

*/
