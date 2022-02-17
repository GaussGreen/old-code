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
#include "gpnummethods/pdetruncators.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: BuildTransitionMatrixes
///	Returns: void
///	Action : Builds Transition Matrixes used in Induct
////////////////////////////////////////////////////
void ARM_PDE2FExplicitNumericalScheme::BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart, double* timeStepEnd)
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
	double stateX,stateY;
	double LastDT;
	double LastVolSquareX;
	double LastVolSquareY;
	size_t matrixSize;

	TimeStepSize = timeSteps->size();
	StatesSize = states->size();
	matrixSize = getSpaceDiscretizationPointsNb();


	/// Predictor and Corrector matrixes
	itsXPDiagElts= ARM_VectorPtrVector(0);
	itsYPDiagElts= ARM_VectorPtrVector(0);

	itsXCDiagElts= ARM_VectorPtrVector(0);
	itsYCDiagElts= ARM_VectorPtrVector(0);

	itsXPOperatorUElts= ARM_VectorPtrVector(0);
	itsXPOperatorLElts= ARM_VectorPtrVector(0);
	itsYPOperatorUElts= ARM_VectorPtrVector(0);
	itsYPOperatorDElts= ARM_VectorPtrVector(0);

	itsXCOperatorUElts= ARM_VectorPtrVector(0);
	itsXCOperatorLElts= ARM_VectorPtrVector(0);
	itsYCOperatorUElts= ARM_VectorPtrVector(0);
	itsYCOperatorDElts= ARM_VectorPtrVector(0);

	itsCrossCoeffs = ARM_GP_VectorPtr( new ARM_GP_Vector(0) );

	/// Boundary limits in case of necessity (truncation)
	/// NullVector = ARM_GP_VectorPtr ( new ARM_GP_Vector(3,0.0) );

	LastDT = 0.0;
	LastVolSquareX = 0.0;
	LastVolSquareY = 0.0;
	size_t MatrixIndex = 0;

	for( TimeIndex = 0 ; TimeIndex < TimeStepSize-1 ; TimeIndex++ )
	{
		/// Precomputes val square for second derivative terms
		volSquareX = (*vols)(0,TimeIndex);
		volSquareX *= volSquareX;
		volSquareY = (*vols)(1,TimeIndex);
		volSquareY *= volSquareY;
		dT = ((*timeSteps)[TimeIndex+1]-(*timeSteps)[TimeIndex])/K_YEAR_LEN;

		if( (LastVolSquareX-volSquareX)>K_DOUBLE_TOL || (LastDT-dT)>K_DOUBLE_TOL || (LastVolSquareX-volSquareX)<-K_DOUBLE_TOL || (LastDT-dT)<-K_DOUBLE_TOL)
		{
			ARM_GP_VectorPtr XPDiagElt( new ARM_GP_Vector(matrixSize,0.0) );
			ARM_GP_VectorPtr YPDiagElt( new ARM_GP_Vector(matrixSize,0.0) );

			ARM_GP_VectorPtr XCDiagElt( new ARM_GP_Vector(matrixSize,0.0) );
			ARM_GP_VectorPtr YCDiagElt( new ARM_GP_Vector(matrixSize,0.0) );

			ARM_GP_VectorPtr XUPElt( new ARM_GP_Vector(matrixSize-1,0.0) );
			ARM_GP_VectorPtr XLPElt( new ARM_GP_Vector(matrixSize-1,0.0) );
			ARM_GP_VectorPtr YUPElt( new ARM_GP_Vector(matrixSize-1,0.0) );
			ARM_GP_VectorPtr YLPElt( new ARM_GP_Vector(matrixSize-1,0.0) );

			ARM_GP_VectorPtr XUCElt( new ARM_GP_Vector(matrixSize-1,0.0) );
			ARM_GP_VectorPtr XLCElt( new ARM_GP_Vector(matrixSize-1,0.0) );
			ARM_GP_VectorPtr YUCElt( new ARM_GP_Vector(matrixSize-1,0.0) );
			ARM_GP_VectorPtr YLCElt( new ARM_GP_Vector(matrixSize-1,0.0) );

			double crossElt=0;

			absoluteDriftX = dT*(*absoluteDrifts)(TimeIndex,0);
			relativeDriftX = (*relativeDrifts)(TimeIndex,0);
			absoluteDriftY = dT*(*absoluteDrifts)(TimeIndex,1);
			relativeDriftY = (*relativeDrifts)(TimeIndex,1);
			dTdX2 = dT/DiscretizationStepXSquare;
			dTdX = dT/DiscretizationStepX;
			dTdY2 = dT/DiscretizationStepYSquare;
			dTdY = dT/DiscretizationStepY;

			/// Precomputes diagonal terms
			double DiagTermPX = 1.0-dTdX2*volSquareX-absoluteDriftX/DiscretizationStepX;
			double DiagTermPY = -dTdY2*volSquareY-absoluteDriftY/DiscretizationStepY;
			double DiagTermCX = 1.0-dTdX2*volSquareX+absoluteDriftX/DiscretizationStepX;
			double DiagTermCY = -dTdY2*volSquareY+absoluteDriftY/DiscretizationStepY;

			/// Precomputes constant terms on the upper and upper diagonals
			double constantTermPXU = 0.5*volSquareX*dTdX2+absoluteDriftX/DiscretizationStepX;
			double constantTermPXD = 0.5*volSquareX*dTdX2;
			double constantTermPYU = 0.5*volSquareY*dTdY2+absoluteDriftY/DiscretizationStepY;
			double constantTermPYD = 0.5*volSquareY*dTdY2;

			double constantTermCXU = 0.5*volSquareX*dTdX2;
			double constantTermCXD = 0.5*volSquareX*dTdX2-absoluteDriftX/DiscretizationStepX;
			double constantTermCYU = 0.5*volSquareY*dTdY2;
			double constantTermCYD = 0.5*volSquareY*dTdY2-absoluteDriftY/DiscretizationStepY;

			double relDriftX = 0;
			double relDriftY = 0;


			/// Loop to fill matrixes (hum, vectors...) ; center part of the matrix
			for( diagIndex=1 ; diagIndex < matrixSize-1 ; diagIndex++ )
			{
				stateX = (*(states->GetNumMethodStates()))(0,diagIndex*matrixSize);
				stateY = (*(states->GetNumMethodStates()))(1,diagIndex);

				relDriftX = stateX*relativeDriftX/DiscretizationStepX;
				relDriftY = stateY*relativeDriftY/DiscretizationStepY;

				XPDiagElt->Elt(diagIndex) = DiagTermPX-relDriftX;
				YPDiagElt->Elt(diagIndex) = DiagTermPY-relDriftY;
				XCDiagElt->Elt(diagIndex) = DiagTermCX+relDriftX;
				YCDiagElt->Elt(diagIndex) = DiagTermCY+relDriftY;

				XUPElt->Elt(diagIndex)=constantTermPXU+relDriftX;
				XLPElt->Elt(diagIndex-1)=constantTermPXD;
				YUPElt->Elt(diagIndex)=constantTermPYU+relDriftY;
				YLPElt->Elt(diagIndex-1)=constantTermPYD;

				XUCElt->Elt(diagIndex)=constantTermCXU;
				XLCElt->Elt(diagIndex-1)=constantTermCXD-relDriftX;
				YUCElt->Elt(diagIndex)=constantTermCYU;
				YLCElt->Elt(diagIndex-1)=constantTermCYD-relDriftY;
			}

			double crossTerm = 0.25*correls->Elt(0,TimeIndex)*sqrt(volSquareX*volSquareY)*dTdX/DiscretizationStepY;

			/// Enqueueing Matrixes
			itsXPDiagElts.push_back(XPDiagElt);
			itsYPDiagElts.push_back(YPDiagElt);
			itsXPOperatorUElts.push_back(XUPElt);
			itsXPOperatorLElts.push_back(XLPElt);
			itsYPOperatorUElts.push_back(YUPElt);
			itsYPOperatorDElts.push_back(YLPElt);

			itsXCDiagElts.push_back(XCDiagElt);
			itsYCDiagElts.push_back(YCDiagElt);
			itsXCOperatorUElts.push_back(XUCElt);
			itsXCOperatorLElts.push_back(XLCElt);
			itsYCOperatorUElts.push_back(YUCElt);
			itsYCOperatorDElts.push_back(YLCElt);

			itsCrossCoeffs->push_back(crossTerm);

			if( TimeIndex != 0 )
				MatrixIndex++;

			LastVolSquareX = volSquareX;
			LastDT = dT;
		}
		setTimeToMatrixIndex(TimeIndex, MatrixIndex);

		LastVolSquareX = volSquareX;
		LastDT = dT;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: toString
///	Returns: string
///	Action : toString
////////////////////////////////////////////////////
string ARM_PDE2FExplicitNumericalScheme ::toString(const string& indent, const string&nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " 2F Explicit Numerical Scheme " << CC_NS(std,endl);
	os << indent << "     Space Discretization : " << getSpaceDiscretizationPointsNb() << " Points" << CC_NS(std,endl);
	os << indent << "     Operture         : " << getStdDevNb() << " Standard Deviations" << CC_NS(std,endl);

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: TruncateMatrixesAndStates
///	Returns: void
///	Action : truncates nummethodstates 
////////////////////////////////////////////////////
void ARM_PDE2FExplicitNumericalScheme ::TruncateMatrixesAndStates( ARM_PricingStatesPtr& states, int toTimeIdx )
{
	ARM_GP_VectorPtr NullVector; /// To be removed if you want to use Truncation
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
	getTruncator()->TruncateTranstionMatrixes( toTimeIdx, itsXPOperatorUElts[MatrixIndex], itsXPOperatorLElts[MatrixIndex], itsXPDiagElts[MatrixIndex], NullVector, NullVector);
	getTruncator()->TruncateTranstionMatrixes( toTimeIdx, itsXCOperatorUElts[MatrixIndex], itsXCOperatorLElts[MatrixIndex], itsXCDiagElts[MatrixIndex], NullVector, NullVector);
	getTruncator()->TruncateTranstionMatrixes( toTimeIdx, itsYPOperatorUElts[MatrixIndex], itsYPOperatorDElts[MatrixIndex], itsYPDiagElts[MatrixIndex], NullVector, NullVector);
	getTruncator()->TruncateTranstionMatrixes( toTimeIdx, itsYCOperatorUElts[MatrixIndex], itsYCOperatorDElts[MatrixIndex], itsYCDiagElts[MatrixIndex], NullVector, NullVector);
	getTruncator()->TruncatePricingStates( toTimeIdx, states );
	size_t newMatrixSize = itsXPDiagElts[MatrixIndex]->size();
	setSpaceDiscretizationPointsNb( newMatrixSize );
	newMatrixSize *= newMatrixSize;
	if( itsPayoffs != ARM_GP_VectorPtr(NULL) )
	{
		itsPayoffs->resize( newMatrixSize );
		itsPayoffs2->resize( newMatrixSize );
	}

}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: Induct
///	Returns: void
///	Action : induct nummethodstates/prices from one date to
///  another
////////////////////////////////////////////////////
void ARM_PDE2FExplicitNumericalScheme ::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	size_t MatrixIndex = getMatrixIdxFromTimeIdx(toTimeIdx);
	
	size_t i;

	//// For truncation (removed yet)
	//// TruncateMatrixesAndStates( states, toTimeIdx );
	
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();

	/// PayoffStatesSize
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);

		ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
		UpdatePayoffs( intermediatePayoffs, toTimeIdx );

		ARM_GP_Vector vec = payoffStates.GetPayoffs(0);
		UpdatePayoffs( vec, toTimeIdx  );
		payoffStates.SetPayoffs( vec );
	}

	/// Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_GP_VectorPtr vec = states->GetPayoffVec(i);
		UpdatePayoffs( *vec, toTimeIdx  );

		for( int j=0; j<vec->size() ; ++j)
			states->SetPayoff( j, i,(*vec)[j] );

	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: UpdatePayoffs
///	Returns: void
///	Action : updates payoffs from one date to another
////////////////////////////////////////////////////
void ARM_PDE2FExplicitNumericalScheme::UpdatePayoffs(ARM_GP_Vector& vec, int toTimeIdx)
{
	/// Gets matrixes needed for update
	size_t matIdx = getMatrixIdxFromTimeIdx(toTimeIdx);
	size_t matrixSize = getSpaceDiscretizationPointsNb()*getSpaceDiscretizationPointsNb();

	ARM_GP_VectorPtr XPDiagElt = itsXPDiagElts[matIdx];
	ARM_GP_VectorPtr YPDiagElt = itsYPDiagElts[matIdx];
	ARM_GP_VectorPtr XPUElt = itsXPOperatorUElts[matIdx];
	ARM_GP_VectorPtr XPLElt = itsXPOperatorLElts[matIdx];
	ARM_GP_VectorPtr YPUElt = itsYPOperatorUElts[matIdx];
	ARM_GP_VectorPtr YPLElt = itsYPOperatorDElts[matIdx];
	ARM_GP_VectorPtr XCDiagElt = itsXCDiagElts[matIdx];
	ARM_GP_VectorPtr YCDiagElt = itsYCDiagElts[matIdx];
	ARM_GP_VectorPtr XCUElt = itsXCOperatorUElts[matIdx];
	ARM_GP_VectorPtr XCLElt = itsXCOperatorLElts[matIdx];
	ARM_GP_VectorPtr YCUElt = itsYCOperatorUElts[matIdx];
	ARM_GP_VectorPtr YCLElt = itsYCOperatorDElts[matIdx];
	double crossTerm = (*itsCrossCoeffs)[matIdx];

	/// Building intermediate vectors 
	if( itsPayoffs == ARM_GP_VectorPtr(NULL) )
		itsPayoffs = ARM_GP_VectorPtr( new ARM_GP_Vector( matrixSize ) );

	if( itsPayoffs2 == ARM_GP_VectorPtr(NULL) )
		itsPayoffs2 = ARM_GP_VectorPtr( new ARM_GP_Vector( matrixSize ) );

	/// Predictor Part
	UpdateVectorWithVectors( vec, *itsPayoffs, XPDiagElt, YPDiagElt, XPUElt, XPLElt, YPUElt, YPLElt, crossTerm ); 
	/// Corrector part
	UpdateVectorWithVectors( *itsPayoffs, *itsPayoffs2, XCDiagElt, YCDiagElt, XCUElt, XCLElt, YCUElt, YCLElt, crossTerm ); 

	/// Mean of the Corrector and Predictor Part
	for( size_t i=0 ; i < matrixSize ; ++i)
		vec[i] = 0.5*(itsPayoffs2->Elt(i) + vec[i]);
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: UpdatePayoffs
///	Returns: void
///	Action : updates payoff from one date to another matrix part. 
///  uses the other one above
////////////////////////////////////////////////////
void ARM_PDE2FExplicitNumericalScheme::UpdatePayoffs(ARM_GP_Matrix& vec, int toTimeIdx)
{
	ARM_GP_Vector vec2(vec.cols());

	for( size_t j=0 ; j<vec.rows() ; ++j )
	{
		for( size_t i=0 ; i<vec.cols() ; ++i)
			vec2[i] = vec(j,i);

		UpdatePayoffs( vec2, toTimeIdx );
	
		for(size_t i=0 ; i<vec.cols() ; ++i)
			vec(j,i) = vec2[i];
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE2FExplicitNumericalScheme 
///	Routine: UpdateVectorWithVectors
///	Returns: void
///	Action : basic matrix multiplication, using vectors
////////////////////////////////////////////////////
void ARM_PDE2FExplicitNumericalScheme::UpdateVectorWithVectors( ARM_GP_Vector& prevVec, ARM_GP_Vector& nextVec, ARM_GP_VectorPtr& XDiag, ARM_GP_VectorPtr& YDiag, ARM_GP_VectorPtr& UU, ARM_GP_VectorPtr& LL, ARM_GP_VectorPtr& U, ARM_GP_VectorPtr& L, double CorrelTerm )
{
	std::vector<double>::iterator prevXY, prevXprevY, prevXnextY, nextXprevY, nextXY, nextXnextY, XnextY, XprevY, curPos, XY;
	size_t matrixSize = getSpaceDiscretizationPointsNb();
	size_t i,j;

	curPos = nextVec.begin();
	XY = prevVec.begin();
	nextXY = XY+matrixSize;
	nextXnextY = nextXY + 1;
	XnextY = XY+1;
	nextXprevY = nextXY-1;

	i=0;j=0;
	XprevY = curPos;

	/// Limit conditions = constant

	/// Xmin
	/// Ymin

	(*curPos) = (*XY);

	XprevY = XY;
	curPos ++; XY++;
	nextXY++;XnextY++;
	nextXnextY++;nextXprevY++;

	/// Middle of the matrix

	for( j=1 ; j<matrixSize ; ++j)
	{
		(*curPos) = (*XY);

		curPos ++; XY++;
		nextXY++;XprevY++;XnextY++;
		nextXnextY++;nextXprevY++;
	}

	/// Xmin YMax
	prevXY = prevVec.begin();
	prevXprevY = prevXY-1;
	prevXnextY = prevXY+1;

	//// Core stuff

	for( i=1; i<matrixSize-1 ; ++i )
	{
		/// Limit condition X=Xmin
		(*curPos) = (*XY);

		curPos ++; XY++;
		prevXY++;nextXY++;XprevY++;XnextY++;
		nextXnextY++;nextXprevY++;prevXnextY++;
		prevXprevY++;

		/// Middle of the matrix
		for( j=1 ; j<matrixSize-1 ; ++j)
		{
			(*curPos) = ( XDiag->Elt(i) + YDiag->Elt(j) ) * (*XY)
				+ (*prevXY) * LL->Elt(i-1)
				+ (*nextXY) * UU->Elt(i)
				+ (*XprevY) * L->Elt(j-1)
				+ (*XnextY) * U->Elt(j)
				+ CorrelTerm*( (*nextXnextY) + (*prevXprevY) - (*nextXprevY) - (*prevXnextY) );

			curPos ++; XY++;
			prevXY++;nextXY++;XprevY++;XnextY++;
			nextXnextY++;prevXprevY++;nextXprevY++;prevXnextY++;
		}

		/// Limit condition X=Xmax
		(*curPos) = (*XY);
		curPos ++; XY++;
		prevXY++;nextXY++;XprevY++;XnextY++;
		nextXnextY++;prevXprevY++;nextXprevY++;prevXnextY++;
	}

	/// Limit conditions = constant
	/// Xmax terms of the matrix
	for( j=0 ; j<matrixSize-1 ; ++j)
	{
		(*curPos) = (*XY);

		curPos ++; XY++;
		prevXY++;XprevY++;XnextY++;
		prevXprevY++;prevXnextY++;
	}

	/// Xmax
	/// Ymax

	(*curPos) = (*XY);

#if defined(__GP_STRICT_VALIDATION)
	curPos++;
	XY++;

	prevXY++;XprevY++;
	prevXprevY++;prevXnextY++;
	prevXnextY;nextXprevY++;

	if( curPos != nextVec.end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "curPos not at the end" );

	if( XY != prevVec.end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "XY not at the end" );

	if( prevXY != XY-matrixSize )
		ARM_THROW( ERR_INVALID_ARGUMENT, "prevXY misplaced" );

	if( prevXprevY != prevXY-1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "prevXprevY misplaced" );

	if( prevXnextY != prevXY+1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "prevXnextY misplaced" );

	if( XprevY != XY-1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "XprevY misplaced" );

	if( XnextY != prevVec.end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "XnextY misplaced" );

	if( nextXnextY != XY+1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "nextXnextY misplaced" );

	if( nextXY != XY )
		ARM_THROW( ERR_INVALID_ARGUMENT, "nextXY misplaced" );

	if( nextXprevY != XY )
		ARM_THROW( ERR_INVALID_ARGUMENT, "nextXnextY misplaced" );

#endif

}

CC_END_NAMESPACE()

/*
*
*	The same way as before, everything is based on tridiagonal matrixes. 
*    Computation is made:
*      P = (AP+BP+CP) C^n+1
*      C = (AC+BC+CC) P
*      C^n = 1/2*(P+C)
*
*	C^n is stored as a statesSize*statesSize vector. 
*    C(Xmin,Ymin), C(Xmin,Ymin+1),  C(Xmin,Ymin+2), ..., C(Xmin+1,Ymin), C(Xmin+1,Ymin+1), ..., C(Xmax,Ymax)
*
*	AP is a tridiagonal matrix that does the X derivatives part in X. It is applied statesSize times on statesSize sized vectors. *
*   Same thing for BP, but index manipulation is somewhat more complicated.
*   CP is stored as a double, (because its value does not depend on X and Y). 
*
*	Same thing for AC,BC and CC. 
*
*
*    AP: 
* itsXPDiagElts;
* itsXPOperatorUElts;
* itsXPOperatorLElts;
* 
*	BP: 
* itsYPDiagElts;
* itsYPOperatorUElts;
* itsYPOperatorDElts;
*
*	AC:
* itsXCDiagElts;
* itsXCOperatorUElts
* itsXCOperatorLElts
*
*   BC: 
* itsYCDiagElts;
* itsYCOperatorUElts
* itsYCOperatorDElts
*
* CC and CP: 
* itsCrossCoeffs
*
*
* P and C are respectively stored in itsPayoffs and itsPayoffs2. 
*
*/ 
