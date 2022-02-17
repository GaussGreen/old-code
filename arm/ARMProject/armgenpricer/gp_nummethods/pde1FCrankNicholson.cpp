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
///	Class  : ARM_PDE1FCrankNicholsonNumericalScheme
///	Routine: BuildTransitionMatrixes
///	Returns: void
///	Action : Build Transisition Matrixes for explicit scheme
////////////////////////////////////////////////////
void ARM_PDE1FCrankNicholsonNumericalScheme::BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart, double* timeStepEnd)
{
	
	/// ***** Builds matrixes for numerical integration ****
	/// Extract the drift part of the Model
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
	model.EulerLocalDrifts(timeSteps->GetValues(),relativeDrifts,absoluteDrifts);

	/// Extract the volatility and correlation part of the Model
    ARM_GP_MatrixPtr vols,d1VolsNotUsed,correlsNotUsed;
    model.VolatilitiesAndCorrelations(timeSteps->GetValues(),vols,d1VolsNotUsed,correlsNotUsed, false/*flatvol*/);

	/// Get model time steps to know when the vol changes
	ARM_GP_VectorPtr modelTimeSteps = model.VolatilitiesAndCorrelationTimesSteps();

	
	size_t TimeStepSize, TimeIndex;
	size_t StatesSize;
	double DiscretizationStep = getSpaceDiscretizationStep();
	double DiscretizationStepSquare = DiscretizationStep*DiscretizationStep;
	double volSquare = 0;
	double absoluteDrift = 0;
	double relativeDrift = 0;
	double dTdX2 = 0;
	double dT = 0;
	double state;
	double LastDT;
	double halfVolSquareDTdX2;

	TimeStepSize = timeSteps->size();
	StatesSize   = states->size();

	/// check if transition matrixes have already been computed
	/// in HK case, we may call this method many times for bootstrap calibration (as the model vol is changing)
	bool Updating = (timeStepStart && timeStepEnd);
		
	if ( !Updating )
	{
		itsNumDiagonalElts		= ARM_VectorPtrVector(0);
		itsNumUpperDiagonalElts = ARM_VectorPtrVector(0);
		itsNumLowerDiagonalElts = ARM_VectorPtrVector(0);
		itsDenomDiagonalElts	= ARM_VectorPtrVector(0);
		itsDenomUpperDiagonalElts = ARM_VectorPtrVector(0);
		itsDenomLowerDiagonalElts = ARM_VectorPtrVector(0);
	}

	/// a little sanity test
	if ( Updating )
	{
		if ( getTimeIndexesToMatrixIndexes()->size() != TimeStepSize )
			ARM_THROW( ERR_INVALID_ARGUMENT, "BuildTransitionMatrixes : unresolved problem with time discretization" );

		// if ( itsLastMatrixIndex != itsNumDiagonalElts.size() )
		//	ARM_THROW( ERR_INVALID_ARGUMENT, "BuildTransitionMatrixes : unresolved problem with time discretization" );
	}
	else
		setTimeIndexesToMatrixIndexes(  ARM_GP_VectorPtr( new ARM_GP_Vector( timeSteps->size() ) ) );


	/// Define start and end of build process
	int startIndex = 0;
	int endIndex   = TimeStepSize - 1;

	/// Manage case where timeStepStart or timeStepEnd is provided
	if (timeStepStart)
	{
		startIndex = -1;
		/// find corresponding index
		for (size_t i(0); i<TimeStepSize; i++)
		{	if ( fabs( (*timeStepStart) -  timeSteps->Elt(i) ) < 1.0e-10 ) 
			{
				startIndex = i;
				break;
			}
		}

		if (startIndex == -1)
			ARM_THROW( ERR_INVALID_ARGUMENT, "BuildTransitionMatrixes : timeStepStart does not belong to the PDE time discretization" );


	}

	if (timeStepEnd)
	{
		endIndex = -1;
		/// find corresponding index
		for (size_t i = startIndex; i<TimeStepSize; i++)
		{	if ( fabs( (*timeStepEnd) -  timeSteps->Elt(i) ) < 1.0e-10 ) 
			{
				endIndex = i;
				break;
			}
		}

		if (endIndex == -1)
			ARM_THROW( ERR_INVALID_ARGUMENT, "BuildTransitionMatrixes : timeStepEnd does not belong to the PDE time discretization" );
		
	}


	///	for truncation
	/*-----------------------------------------------
	itsNumULLimitConditions = ARM_VectorPtrVector(0);
	itsNumDRLimitConditions = ARM_VectorPtrVector(0);
	itsDenomULLimitConditions = ARM_VectorPtrVector(0);
	itsDenomDRLimitConditions = ARM_VectorPtrVector(0);

	ARM_GP_VectorPtr NumULLimitConditions( new ARM_GP_Vector(3,0.0) );
	ARM_GP_VectorPtr NumDRLimitConditions( new ARM_GP_Vector(3,0.0) );
	ARM_GP_VectorPtr DenomULLimitConditions( new ARM_GP_Vector(3,0.0) );
	ARM_GP_VectorPtr DenomDRLimitConditions( new ARM_GP_Vector(3,0.0) );
	-----------------------------------------------*/

	LastDT = 0.0;
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(startIndex);
	
	/// find ModelIndex
	size_t ModelIndex  = 0 ;
	for (size_t i(0); i<modelTimeSteps->size(); i++)
	{
		if ( modelTimeSteps->Elt(i) == timeSteps->Elt(startIndex)  )
		{
			ModelIndex = i ;
			break;
		}
	}
	
	for( TimeIndex = startIndex ; TimeIndex < endIndex ; TimeIndex++ )
	{
		/// Compute Vol Square and deltaT
		volSquare = (*vols)(0,TimeIndex);
		volSquare *= volSquare;
		dT = (timeSteps->Elt(TimeIndex+1)-timeSteps->Elt(TimeIndex))/K_YEAR_LEN;

		/// We build a new matrix only if model vol or deltaT have changed
		if( fabs(modelTimeSteps->Elt(ModelIndex) - timeSteps->Elt(TimeIndex))<K_DOUBLE_TOL || fabs(LastDT-dT)>K_DOUBLE_TOL )
		{
			ARM_GP_VectorPtr NumDiagElt( new ARM_GP_Vector(StatesSize) );
			ARM_GP_VectorPtr NumUpperDiagElt( new ARM_GP_Vector(StatesSize-1) );
			ARM_GP_VectorPtr NumLowerDiagElt( new ARM_GP_Vector(StatesSize-1) );
			ARM_GP_VectorPtr DenomDiagElt( new ARM_GP_Vector(StatesSize) );
			ARM_GP_VectorPtr DenomUpperDiagElt( new ARM_GP_Vector(StatesSize-1) );
			ARM_GP_VectorPtr DenomLowerDiagElt( new ARM_GP_Vector(StatesSize-1) );

			absoluteDrift = dT*(*absoluteDrifts)(TimeIndex,0);
			relativeDrift = (*relativeDrifts)(TimeIndex,0);
			dTdX2 = dT / DiscretizationStepSquare;
			halfVolSquareDTdX2 = 0.5 * volSquare * dTdX2;

			absoluteDrift /= DiscretizationStep;
			relativeDrift /= DiscretizationStep;

			std::vector<double>::iterator 
				iterNumDiag = NumDiagElt->begin()+1,
				iterNumLow  = NumLowerDiagElt->begin(),
				iterNumUp   = NumUpperDiagElt->begin() + 1,
				iterDenDiag = DenomDiagElt->begin() + 1,
				iterDenLow  = DenomLowerDiagElt->begin(),
				iterDenUp   = DenomUpperDiagElt->begin() + 1,

				end		   = NumDiagElt->end()-1;

			ARM_GP_Matrix::iterator iterState = states->GetNumMethodStates()->begin() + 1;

			/// Building the central part of the matrix
			for (; iterNumDiag != end; ++iterNumDiag, ++iterNumLow, ++iterNumUp, ++iterDenDiag, ++iterDenLow, ++iterDenUp, ++iterState)
			{				
				(*iterNumDiag) = 1 - halfVolSquareDTdX2;
				(*iterNumLow)  = 0.5 * (-0.5 * relativeDrift * (*iterState) + absoluteDrift + halfVolSquareDTdX2);
				(*iterNumUp)   = 0.5 * ( 0.5 * relativeDrift * (*iterState) + absoluteDrift + halfVolSquareDTdX2);

				(*iterDenDiag) = 1.0 + halfVolSquareDTdX2;
				(*iterDenLow)  = -0.5 * (-0.5 * relativeDrift * (*iterState)  + absoluteDrift + halfVolSquareDTdX2);
				(*iterDenUp)   = -0.5 * ( 0.5 * relativeDrift * (*iterState)  + absoluteDrift + halfVolSquareDTdX2);
			}
	

			/// Building matrix limit conditions
			/// Von Neumann limit conditions
			state					= (*(states->GetNumMethodStates()))(0,0);
			(*NumDiagElt)[0]		= 1 + 0.5 * (absoluteDrift - relativeDrift * state);
			(*NumUpperDiagElt)[0]	= 0.5 * relativeDrift * state;
			(*DenomDiagElt)[0]		= 1 - 0.5 * (absoluteDrift- relativeDrift * state);
			(*DenomUpperDiagElt)[0] = -0.5 * relativeDrift * state;
			
			state								= (*(states->GetNumMethodStates()))(0,StatesSize-1);
			(*NumDiagElt)[StatesSize-1]			= 1 + 0.5 * (absoluteDrift + relativeDrift*state);
			(*NumLowerDiagElt)[StatesSize-2]	= - 0.5 * relativeDrift * state;
			(*DenomDiagElt)[StatesSize-1]		= 1 - 0.5 * (absoluteDrift + relativeDrift*state);
			(*DenomLowerDiagElt)[StatesSize-2]	= 0.5 * relativeDrift * state;


			if( TimeIndex != startIndex )
				MatrixIndex++;

			if ( !Updating )
			{
				/// Enqueueing numerator
				itsNumDiagonalElts.push_back( NumDiagElt );
				itsNumUpperDiagonalElts.push_back( NumUpperDiagElt );
				itsNumLowerDiagonalElts.push_back( NumLowerDiagElt );

				/// Enqueueing denominator
				itsDenomDiagonalElts.push_back( DenomDiagElt );
				itsDenomUpperDiagonalElts.push_back( DenomUpperDiagElt );
				itsDenomLowerDiagonalElts.push_back( DenomLowerDiagElt );
			}
			else
			{
				/// Enqueueing numerator
				itsNumDiagonalElts[MatrixIndex]		 = NumDiagElt ;
				itsNumUpperDiagonalElts[MatrixIndex] = NumUpperDiagElt;
				itsNumLowerDiagonalElts[MatrixIndex] = NumLowerDiagElt;

				/// Enqueueing denominator
				itsDenomDiagonalElts[MatrixIndex]		= DenomDiagElt;
				itsDenomUpperDiagonalElts[MatrixIndex]	= DenomUpperDiagElt ;
				itsDenomLowerDiagonalElts[MatrixIndex]	= DenomLowerDiagElt ;

			}
			
			LastDT = dT;

			if( fabs(modelTimeSteps->Elt(ModelIndex) -timeSteps->Elt(TimeIndex))<K_DOUBLE_TOL && (ModelIndex<modelTimeSteps->size()-1) )
				ModelIndex ++;
		}
		setTimeToMatrixIndex(TimeIndex, MatrixIndex);
	
		LastDT = dT;
	}

	if ( !Updating )
		itsLastMatrixIndex = MatrixIndex+1;

	/// Allocates matrixes for inversion
	itsOtherCoeff = ARM_GP_VectorPtr( new ARM_GP_Vector( StatesSize-1 ) );
	itsNormalizationTerms = ARM_GP_VectorPtr( new ARM_GP_Vector( StatesSize ) );

	/// Uncomment this if you are brave (Truncation)
	/*-----------------------------------------------------------------------
	if( TruncationSizeVector.size() != 0 )
	{
		for( TimeIndex = 0 ; TimeIndex<TimeStepSize-1 ; ++TimeIndex )
		{
			if( TruncationSizeVector[TimeIndex] != TruncationSizeVector[TimeIndex+1] )
			{
				NumULLimitConditions = ARM_GP_VectorPtr( new ARM_GP_Vector(3) );
				NumDRLimitConditions = ARM_GP_VectorPtr( new ARM_GP_Vector(3) );
				DenomULLimitConditions = ARM_GP_VectorPtr( new ARM_GP_Vector(3) );
				DenomDRLimitConditions = ARM_GP_VectorPtr( new ARM_GP_Vector(3) );

				dT = (timeSteps->Elt(TimeIndex+1)-timeSteps->Elt(TimeIndex))/K_YEAR_LEN;
				dTdX2 = dT/DiscretizationStepSquare;
				absoluteDrift = dT*(*absoluteDrifts)(TimeIndex,0);
				relativeDrift = (*relativeDrifts)(TimeIndex,0);


				state = (*(states->GetNumMethodStates()))(0, (StatesSize-TruncationSizeVector[TimeIndex])/2 );
				/// UpperDiag
				NumULLimitConditions->Elt(0) = 0.5*relativeDrift*state/DiscretizationStep;
				DenomULLimitConditions->Elt(0) = -0.5*relativeDrift*state/DiscretizationStep;
				/// Diag
				NumULLimitConditions->Elt(1) = 1 + 0.5*(absoluteDrift- relativeDrift*state)/DiscretizationStep;
				DenomULLimitConditions->Elt(1) = 1 - 0.5*(absoluteDrift- relativeDrift*state)/DiscretizationStep;
				state = (*(states->GetNumMethodStates()))(0, (StatesSize-TruncationSizeVector[TimeIndex])/2+1 );
				/// LowerDiag
				NumULLimitConditions->Elt(2) = 0.5*( (-0.5*relativeDrift*state+ absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2 );
				DenomULLimitConditions->Elt(2) = -0.5*( (-0.5*relativeDrift*state+ absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2 );



				state = (*(states->GetNumMethodStates()))(0, (StatesSize+TruncationSizeVector[TimeIndex])/2-2 );
				/// UpperDiag
				NumDRLimitConditions->Elt(0) = 0.5* ( (0.5*relativeDrift*state + absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2 );
				DenomDRLimitConditions->Elt(0) = -0.5* ( (0.5*relativeDrift*state + absoluteDrift )/DiscretizationStep + 0.5*volSquare*dTdX2 );
				state = (*(states->GetNumMethodStates()))(0, (StatesSize+TruncationSizeVector[TimeIndex])/2-1 );
				/// Diag
				NumDRLimitConditions->Elt(1) = 1+0.5*(absoluteDrift+ relativeDrift*state)/DiscretizationStep ;
				DenomDRLimitConditions->Elt(1) = 1-0.5*(absoluteDrift+ relativeDrift*state)/DiscretizationStep ;
				/// LowerDiag
				NumDRLimitConditions->Elt(2) = -0.5*relativeDrift*state/DiscretizationStep;
				DenomDRLimitConditions->Elt(2) = 0.5*relativeDrift*state/DiscretizationStep;
			}
			itsNumULLimitConditions.push_back(NumULLimitConditions);
			itsNumDRLimitConditions.push_back(NumDRLimitConditions);
			itsDenomULLimitConditions.push_back(DenomULLimitConditions);
			itsDenomDRLimitConditions.push_back(DenomDRLimitConditions);
		}
	}
	-----------------------------------------------------------------------*/
	/// End Truncation

}


////////////////////////////////////////////////////
///	Class  : ARM_PDE1FCrankNicholsonNumericalScheme
///	Routine: Induct
///	Returns: void
///	Action : Inducts payoffs from one TimeIDx to the next(previous one)
////////////////////////////////////////////////////
void ARM_PDE1FCrankNicholsonNumericalScheme::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	/// Gets the matrixes needed for induct
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
	/// Numerator Matrixes
	ARM_GP_VectorPtr NumDiagElt = itsNumDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr NumUpperDiagElt = itsNumUpperDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr NumLowerDiagElt = itsNumLowerDiagonalElts[MatrixIndex];
	/// Denominator Matrixes
	ARM_GP_VectorPtr DenomDiagElt = itsDenomDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr DenomUpperDiagElt = itsDenomUpperDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr DenomLowerDiagElt = itsDenomLowerDiagonalElts[MatrixIndex];
	
	size_t i;

	// Truncation part
	/*getTruncator()->TruncateTranstionMatrixes( toTimeIdx, NumUpperDiagElt, NumLowerDiagElt, NumDiagElt, itsNumULLimitConditions[toTimeIdx], itsNumDRLimitConditions[toTimeIdx]);
	getTruncator()->TruncateTranstionMatrixes( toTimeIdx, DenomUpperDiagElt, DenomLowerDiagElt, DenomDiagElt, itsDenomULLimitConditions[toTimeIdx], itsDenomDRLimitConditions[toTimeIdx]);
	getTruncator()->TruncatePricingStates( toTimeIdx, states );*/

  	if( DenomDiagElt->size() != itsNormalizationTerms->size() )
	{
		itsNormalizationTerms->resize( DenomDiagElt->size() );
		itsOtherCoeff->resize( DenomUpperDiagElt->size() );
		UpdateDiagsForInversion( MatrixIndex );
	}
	
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();
	size_t ProbaChangeSize = states->GetProbaChangesSize();

	/// Precomputes inverse matrix only if retropropagation matrix has changed
	if( true )  /// if( MatrixIndex != itsLastMatrixIndex ) /// commented to be sure there is no problem when induncting several times on same time interval
		UpdateDiagsForInversion( MatrixIndex );

	itsLastMatrixIndex = MatrixIndex;
 
	/// PayoffStatesSize
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);

		ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix& >(payoffStates.GetIntermediatePayoffs() );
		UpdateVectorWithMatrixByProduct( NumDiagElt, NumUpperDiagElt, NumLowerDiagElt, intermediatePayoffs );
		UpdateVectorWithMatrixByInverse( DenomDiagElt, DenomUpperDiagElt, DenomLowerDiagElt, intermediatePayoffs );

		std::vector<double>& vec = const_cast<std::vector<double>&> (payoffStates.GetPayoffs());
		UpdateVectorWithMatrixByProduct( NumDiagElt, NumUpperDiagElt, NumLowerDiagElt, vec );
		UpdateVectorWithMatrixByInverse( DenomDiagElt, DenomUpperDiagElt, DenomLowerDiagElt, vec );
//		payoffStates.SetPayoffs( vec );
	}

	/// Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_MemPool_Matrix::iterator vecIter = states->payoffsBeginIterator(i);
		UpdateVectorWithMatrixByProduct( NumDiagElt, NumUpperDiagElt, NumLowerDiagElt, vecIter );
		UpdateVectorWithMatrixByInverse( DenomDiagElt, DenomUpperDiagElt, DenomLowerDiagElt, vecIter );
	}

	for (i=0;i<ProbaChangeSize;++i)
	{
		ARM_VectorPtr vec = states->GetProbaChangeVec(i);
		UpdateVectorWithMatrixByProduct( NumDiagElt, NumUpperDiagElt, NumLowerDiagElt, *vec );
		UpdateVectorWithMatrixByInverse( DenomDiagElt, DenomUpperDiagElt, DenomLowerDiagElt, *vec );

		for( int j=0; j<vec->size() ; ++j)
			states->SetProbaChange( j, i, (*vec)[j] );
	}

}


////////////////////////////////////////////////////
///	Class  : ARM_PDE1FCrankNicholsonNumericalScheme
///	Routine: UpdateVectorWithMatrixByInverse
///	Returns: void
///	Action : updates payoffs in TransposedVector 
///  using tridiag matrix given as parameter. Inverse
///  has been precomputed
////////////////////////////////////////////////////
void ARM_PDE1FCrankNicholsonNumericalScheme::UpdateVectorWithMatrixByInverse( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_GP_Matrix& TransposedVector )
{
	size_t j, linesNb = TransposedVector.rows(); 
	size_t i;
	size_t StateSize = TransposedVector.cols();
	std::vector<double>::iterator normIter, lowerIter, vecIter;

	vecIter = TransposedVector.begin();

	for( j=0 ; j<linesNb ; ++j, vecIter+=StateSize )
	{
		normIter = itsNormalizationTerms->begin();
		double CurrentState = (*vecIter)/(*normIter);

		(*vecIter) = CurrentState;
		++vecIter;
		++normIter;
		lowerIter = LowerTerm->begin();

		for( i=1 ; i< StateSize ; ++i, ++vecIter, ++lowerIter, ++normIter)
		{
			CurrentState = ((*vecIter)-(*lowerIter)*CurrentState)/(*normIter);
			(*vecIter) = CurrentState;
		}

#if defined(__GP_STRICT_VALIDATION)
		if( normIter != itsNormalizationTerms->end()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "normalization not at the end" );
		if( lowerIter!= LowerTerm->end()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "lowerTerm not at the end" );
		if( vecIter != TransposedVector.begin()+(1+j)*StateSize )
			ARM_THROW( ERR_INVALID_ARGUMENT, "vecIter ill placed" );
#endif

		--vecIter;
		--vecIter;
		normIter = itsOtherCoeff->end()-1;

		for( i = StateSize-2 ; i>0 ; --i, --vecIter, --normIter )
		{
			CurrentState = (*vecIter)-(*normIter)*CurrentState;
			(*vecIter) = CurrentState;
		}

#if defined(__GP_STRICT_VALIDATION)
		if( normIter != itsOtherCoeff->begin()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "itsOtherCoeff not at the beginning" );
		if( vecIter != TransposedVector.begin()+j*StateSize )
			ARM_THROW( ERR_INVALID_ARGUMENT, "vecIter ill placed" );
#endif

		/// Smoothing
		TransposedVector(j,0) = 2*TransposedVector(j,1) - TransposedVector(j,2);
		TransposedVector(j,StateSize-1) = 2*TransposedVector(j,StateSize-2)-TransposedVector(j,StateSize-3);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FCrankNicholsonNumericalScheme
///	Routine: UpdateVectorWithMatrixByInverse
///	Returns: void
///	Action : updates payoffs in vec 
///  using tridiag matrix given as parameter
///    Same thing as above but for solely vector
////////////////////////////////////////////////////
void ARM_PDE1FCrankNicholsonNumericalScheme::UpdateVectorWithMatrixByInverse( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, std::vector<double>& vec )
{
	size_t i;
	size_t StateSize = vec.size();
	std::vector<double>::iterator normIter, lowerIter, vecIter;

	vecIter = vec.begin();
	normIter = itsNormalizationTerms->begin();
	double CurrentState = (*vecIter)/(*normIter);

	(*vecIter) = CurrentState;
	++vecIter;++normIter;
	lowerIter = LowerTerm->begin();

	for( i=1 ; i< StateSize ; ++i, ++vecIter, ++lowerIter, ++normIter)
	{
		CurrentState = ((*vecIter)-(*lowerIter)*CurrentState)/(*normIter);
		(*vecIter) = CurrentState;
	}

#if defined(__GP_STRICT_VALIDATION)
		if( normIter != itsNormalizationTerms->end()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "normalization not at the end" );
		if( lowerIter!= LowerTerm->end()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "lowerTerm not at the end" );
		if( vecIter != vec.end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "vecIter ill placed" );
#endif

	--vecIter;
	--vecIter;

	normIter = itsOtherCoeff->end()-1;

	for( i = StateSize-2 ; i>0 ; --i, --vecIter, --normIter )
	{
		CurrentState = (*vecIter)-(*normIter)*CurrentState;
		(*vecIter) = CurrentState;
	}

#if defined(__GP_STRICT_VALIDATION)
	if( normIter != itsOtherCoeff->begin() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "itsOtherCoeff iter != begin" );
	if( vecIter != vec.begin() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "vecIterator != begin" );
#endif

	/// Smoothing
	vec[0] = 2*vec[1]-vec[2];
	vec[StateSize-1] = 2*vec[StateSize-2]-vec[StateSize-3];
}

void ARM_PDE1FCrankNicholsonNumericalScheme::UpdateVectorWithMatrixByInverse( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_MemPool_Matrix::iterator vecBegin )
{
	size_t i;
	size_t StateSize = DiagTerm->size();
	std::vector<double>::iterator normIter, lowerIter;

	ARM_MemPool_Matrix::iterator vecIter = vecBegin;
	normIter = itsNormalizationTerms->begin();
	double CurrentState = (*vecIter)/(*normIter);

	(*vecIter) = CurrentState;
	++vecIter;++normIter;
	lowerIter = LowerTerm->begin();

	for( i=1 ; i< StateSize ; ++i, ++vecIter, ++lowerIter, ++normIter)
	{
		CurrentState = ((*vecIter)-(*lowerIter)*CurrentState)/(*normIter);
		(*vecIter) = CurrentState;
	}

#if defined(__GP_STRICT_VALIDATION)
		if( normIter != itsNormalizationTerms->end()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "normalization not at the end" );
		if( lowerIter!= LowerTerm->end()  )
			ARM_THROW( ERR_INVALID_ARGUMENT, "lowerTerm not at the end" );
#endif

	--vecIter;
	--vecIter;

	normIter = itsOtherCoeff->end()-1;

	for( i = StateSize-2 ; i>0 ; --i, --vecIter, --normIter )
	{
		CurrentState = (*vecIter)-(*normIter)*CurrentState;
		(*vecIter) = CurrentState;
	}

#if defined(__GP_STRICT_VALIDATION)
	if( normIter != itsOtherCoeff->begin() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "itsOtherCoeff iter != begin" );
#endif

	/// Smoothing
	(*vecIter) = 2* (*(vecIter+1))-(*(vecIter+2));
	*(vecIter+StateSize-1) = 2* (*(vecIter+StateSize-2)) - *(vecIter+StateSize-3);
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FCrankNicholsonNumericalScheme
///	Routine: UpdateDiagsForInversion
///	Returns: void
///	Action : precomputes and stores stuff for 
///   enhanced matrix inversion in UpdateVectorWithMatrixByInverse
///    See thomas Algortihm
////////////////////////////////////////////////////
void ARM_PDE1FCrankNicholsonNumericalScheme::UpdateDiagsForInversion( size_t MatrixIndex )
{
	size_t i;
	ARM_GP_VectorPtr LowerTerm = itsDenomLowerDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr DiagTerm = itsDenomDiagonalElts[MatrixIndex];
	ARM_GP_VectorPtr UpperTerm = itsDenomUpperDiagonalElts[MatrixIndex];

	size_t MatrixSize = DiagTerm->size();

	double Normalisation = DiagTerm->Elt(0);
	double OtherCoeff = UpperTerm->Elt(0) / Normalisation;

	(*itsNormalizationTerms)[0] = Normalisation;
	(*itsOtherCoeff)[0] = OtherCoeff;

	for( i=1 ; i<MatrixSize-1 ; ++i )
	{
		Normalisation = DiagTerm->Elt(i)-LowerTerm->Elt(i-1)*OtherCoeff;
		OtherCoeff = UpperTerm->Elt(i)/Normalisation;
		(*itsNormalizationTerms)[i] = Normalisation;
		(*itsOtherCoeff)[i] = OtherCoeff;
	}

	Normalisation = DiagTerm->Elt(MatrixSize-1)-LowerTerm->Elt(MatrixSize-2)*OtherCoeff;
	(*itsNormalizationTerms)[MatrixSize-1] = Normalisation;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FCrankNicholsonNumericalScheme
///	Routine: toString
///	Returns: string
///	Action : As usual for a toString method
////////////////////////////////////////////////////
string ARM_PDE1FCrankNicholsonNumericalScheme::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " 1F Crank Nicholson Numerical Scheme " << CC_NS(std,endl);
	os << indent << "     Space Discretization : " << getSpaceDiscretizationPointsNb() << " Points" << CC_NS(std,endl);
	os << indent << "     Operture         : " << getStdDevNb() << " Standard Deviations" << CC_NS(std,endl);

#if defined(__GP_STRICT_VALIDATION)
	if (!StoredVols.IsNull())
	{
		os << indent << CC_NS(std,endl) << indent << "Used Instant Vols:" << CC_NS(std,endl);
		os << StoredVols->toString();
		os << CC_NS(std,endl);
	}
#endif

	return os.str();
}
CC_END_NAMESPACE()

/*
* Once again, all matrixes are tridiagonal. For their implentation, see explicit scheme. 
*  (I+M) C^n+1 = (I-M) C^n
* Which we will also write
*  D C^n = N C^n+1
*
* where D and N are stored (respectively) in
*  
*	itsNumDiagonalElts, itsNumUpperDiagonalElts, itsNumLowerDiagonalElts;
* and 
*   itsDenomDiagonalElts, itsDenomUpperDiagonalElts, itsDenomLowerDiagonalElts;
*
*
*  As far as matrix inversion is concerned. The last used matrix index is stored in itsLastMatrixIndex. 
*  Coefficients of the thomas algorithms are stored in itsNormalizationTerms and itsOtherCoeff. 
*  if itsLastMatrixIndex != getMatrixIdxFromTimeIdx(toTimeIdx), then itsNormalizationTerms and itsOtherCoeff
*  are updated (using UpdateDiagsForInversion). 
*
*/
