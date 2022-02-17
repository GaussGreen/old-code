/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pdenumericalschemes.cpp
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */

#include "gpnummethods/pdenumericalschemes.h"
#include "gpnummethods/pde1Fnumericalschemes.h"
#include "gpnummethods/pde2Fnumericalschemes.h"
#include "gpnummethods/pde3Fnumericalschemes.h"
#include "gpnummethods/pdetruncators.h"
#include "gpinfra/pricingmodel.h"
#include "gpbase/vectormanip.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpinfra/pricingstates.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/gplinalgconvert.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Checks if the model is compatible with the 
///   current numerical scheme. Basically checks that model is 1F
////////////////////////////////////////////////////
ARM_PDENumericalScheme* ARM_PDENumericalScheme::getNumericalSchemeInstanceById( int NumSchemeId, size_t NX, size_t NY, size_t NZ, double Theta1, double Theta2, double Theta3, int BoundConditionId, double lambda,  int gridType, ARM_GP_Matrix gridData, ARM_GP_Matrix schedulerData)
{
	switch( NumSchemeId )
	{
	case ARM_PDENumericalScheme::None :
		return NULL;
		break;
	case ARM_PDENumericalScheme::Explicit1F :
		return new ARM_PDE1FExplicitNumericalScheme();
		break;
	case ARM_PDENumericalScheme::CN1F :
		return new ARM_PDE1FCrankNicholsonNumericalScheme();
		break;
	case ARM_PDENumericalScheme::Explicit2F :
		return new ARM_PDE2FExplicitNumericalScheme();
		break;
	case ARM_PDENumericalScheme::ADI2F :
		return new ARM_PDE2FADINumericalScheme();
		break;
	case ARM_PDENumericalScheme::CS3F :
		return new ARM_PDE3FCraigSneydNumericalScheme(NX,NY,NZ,Theta1,Theta2,Theta3,(ARM_PDE3FCraigSneydNumericalScheme::BoundConditionType)BoundConditionId,lambda,(ARM_PDE3FNumericalScheme::GridType)gridType,gridData, schedulerData);
		break;

	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, "Unknown Numerical Scheme" );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDENumericalScheme
///	Routine: ~ARM_PDENumericalScheme
///	Returns: ...
///	Action : destructor
////////////////////////////////////////////////////
ARM_PDENumericalScheme::~ARM_PDENumericalScheme()
{ 
	delete itsTruncator;
	DeletePointorVector<ARM_GP_Matrix>(itsGlobalVariances);
}

////////////////////////////////////////////////////
///	Class  : ARM_PDENumericalScheme
///	Routine: ARM_PDENumericalScheme
///	Returns: ...///	Returns: ...
///	Action : constructor
////////////////////////////////////////////////////
ARM_PDENumericalScheme::ARM_PDENumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator) : itsSpaceDiscretizationPointsNb(DiscretizationPointsNb), itsStdDevNb( StdDevNb ),itsSpaceDiscretizationStep(0), itsTruncator(Truncator) 
{
	if (!itsTruncator)
		itsTruncator = new PDE_DummyTruncator();
}

////////////////////////////////////////////////////
///	Class  : ARM_PDENumericalScheme
///	Routine: ARM_PDENumericalScheme
///	Returns: ...
///	Action : cpoy constructor
////////////////////////////////////////////////////
ARM_PDENumericalScheme::ARM_PDENumericalScheme( const ARM_PDENumericalScheme& rhs ) : 
itsSpaceDiscretizationPointsNb(rhs.itsSpaceDiscretizationPointsNb), 
itsStdDevNb(rhs.itsStdDevNb),
itsSpaceDiscretizationStep(rhs.itsSpaceDiscretizationStep), 
itsTruncator(NULL)
{
	itsTruncator = static_cast<PDE_Truncator*> (rhs.itsTruncator->Clone() );
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Checks if the model is compatible with the 
///   current numerical scheme. Basically checks that model is 1F
////////////////////////////////////////////////////
bool ARM_PDE1FNumericalScheme::CheckCompatibilityWithModel( const ARM_PricingModel& model )
{
	return ( (model.FactorCount() == 1) ? true : 0);
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Initialiation of the NumericalScheme
///  builds nummethodstates and discretization schemes. 
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PDE1FNumericalScheme::Init( const ARM_PricingModel& model )
{
	ARM_GP_VectorPtr timeSteps( static_cast<ARM_GP_Vector*> (model.GetNumMethod()->GetTimeSteps()->Clone() ));
	itsGlobalVariances = ARM_MatrixVector(0);
	model.NumMethodStateGlobalVariances( *timeSteps, itsGlobalVariances );

	int SpaceDiscretizationPointsNb = getSpaceDiscretizationPointsNb();

	//////////////////////////////////////////////////////////////////////////////////
	///// For testing purpose ONLY (truncation)
	/// ARM_GP_T_Vector<size_t> TruncationVector(timeSteps->size());
	/// for (int k=0 ; k < timeSteps->size() ; ++k)
	/// TruncationVector[k] = ((SpaceDiscretizationPointsNb+(SpaceDiscretizationPointsNb*k)/timeSteps->size()) % 2) ? (SpaceDiscretizationPointsNb+(SpaceDiscretizationPointsNb*k)/timeSteps->size()) : (SpaceDiscretizationPointsNb+(SpaceDiscretizationPointsNb*k)/timeSteps->size()-1);
	/// setTruncator( new PDE1F_Truncator() );
	/// getTruncator()->setToTimeSizes( TruncationVector );
	//////////////////////////////////////////////////////////////////////////////////

	/// Price Will be in the middle of the Grid
	int PriceIndex = (SpaceDiscretizationPointsNb-1)/2;
	setPriceIndex( PriceIndex );

	//////////////////////////////////////////////////////////////////////////////////
	///////////TEST
	/// SpaceDiscretizationPointsNb*=2;
	/// SpaceDiscretizationPointsNb-=1;
	/// setSpaceDiscretizationPointsNb(SpaceDiscretizationPointsNb);
	//////////////////////////////////////////////////////////////////////////////////

	/// Computation of DeltaX
	double DiscretizationStep = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/SpaceDiscretizationPointsNb;
	setSpaceDiscretizationStep( DiscretizationStep );

	/// *****  building NumMethodStates *******

	/// Building Pricing States ; its size is equal to the number of DiscreizationPoints
	ARM_PricingStatesPtr states( model.FirstPricingStates( SpaceDiscretizationPointsNb ) );
	
	/// Building NumMethodSates
	ARM_GP_Matrix * matrix = new ARM_GP_Matrix( 1, SpaceDiscretizationPointsNb );

	for( int i=0 ; i  <SpaceDiscretizationPointsNb ; ++i )
		(*matrix)(0,i) = DiscretizationStep*(i-(SpaceDiscretizationPointsNb-1)/2);

	states->SetNumMethodStates(	ARM_GP_MatrixPtr(matrix) );
	
	/// build transition matrixes
	BuildTransitionMatrixes( model, timeSteps, states );

	/// Returns the first pricing states
	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FNumericalScheme
///	Routine: UpdateVectorWithMatrixByProduct
///	Returns: void
///	Action : updates payoffs in TransposedVector 
///  using tridiag matrix given as parameter
////////////////////////////////////////////////////
void ARM_PDE1FNumericalScheme::UpdateVectorWithMatrixByProduct( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_GP_Matrix& TransposedVector )
{
	size_t linesNb = TransposedVector.rows();
	size_t size = TransposedVector.cols();
	size_t i = 0;
	ARM_GP_Vector::iterator vecIter, diagIter, uIter, lIter;
	vecIter = TransposedVector.begin();

	for( size_t j=0 ; j< linesNb ; ++j)
	{
		double N,NM1,NP1;
		N=(*vecIter);
		NP1=*(vecIter+1);
		NM1=0;
		diagIter = DiagTerm->begin();
		uIter = UpperTerm->begin();

		(*vecIter) = (*diagIter)*N + (*uIter) * NP1;
		NM1=N;
		N=NP1;

		++vecIter;
		++diagIter;
		++uIter;
		lIter = LowerTerm->begin();

		for( i=1; i<size-1; ++i, ++vecIter, ++lIter, ++uIter, ++diagIter)
		{
			NP1 = *(vecIter+1);
			(*vecIter) = (*lIter)*NM1 + (*diagIter)*N + (*uIter) * NP1;
			NM1 = N;
			N=NP1;
		}
		(*vecIter) = (*lIter)*NM1 + (*diagIter)*N;
		++vecIter;

#if defined(__GP_STRICT_VALIDATION)
	++diagIter;
	++lIter;
	if( vecIter != TransposedVector.begin()+(j+1)*size )
		ARM_THROW( ERR_INVALID_ARGUMENT, "veciter ill palced" );
	if( diagIter != DiagTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "diagIter =! end" );
	if( uIter != UpperTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "upperIter =! end" );
	if( lIter != LowerTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "lowerIter =! end" );
#endif
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE1FNumericalScheme
///	Routine: UpdateVectorWithMatrixByProduct
///	Returns: void
///	Action : updates payoffs in vec from one date 
///   using tridiag matrix given as parameter
////////////////////////////////////////////////////
void ARM_PDE1FNumericalScheme::UpdateVectorWithMatrixByProduct( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_GP_Vector& vec )
{
	size_t size = vec.size();
	size_t i = 0;
	double N,NM1,NP1;
	ARM_GP_Vector::iterator vecIter, lIter, uIter, diagIter;
	vecIter = vec.begin();
	diagIter = DiagTerm->begin();
	uIter = UpperTerm->begin();
	lIter = LowerTerm->begin();
	N=(*vecIter);
	NP1=*(vecIter+1);
	NM1=0;

	(*vecIter) = (*diagIter)*N + (*uIter) * NP1;
	++uIter;++diagIter;	++vecIter;
	NM1=N;
	N=NP1;

	for( i=1; i<size-1; ++i, ++uIter, ++diagIter, ++lIter, ++vecIter)
	{
		NP1 = *(vecIter+1);
		(*vecIter) = (*lIter)*NM1 + (*diagIter)*N + (*uIter) * NP1;
		NM1 = N;
		N=NP1;
	}
	(*vecIter) = (*lIter)*NM1 + (*diagIter)*N;

#if defined(__GP_STRICT_VALIDATION)
	++vecIter;
	++diagIter;
	++lIter;
	if( vecIter != vec.end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "veciter ill palced" );
	if( diagIter != DiagTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "diagIter =! end" );
	if( uIter != UpperTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "upperIter =! end" );
	if( lIter != LowerTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "lowerIter =! end" );
#endif
}

void ARM_PDE1FNumericalScheme::UpdateVectorWithMatrixByProduct( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_MemPool_Matrix::iterator vecBegin )
{
	size_t size = DiagTerm->size();
	size_t i = 0;
	double N,NM1,NP1;
	ARM_GP_Vector::iterator lIter, uIter, diagIter;
	ARM_MemPool_Matrix::iterator vecIter = vecBegin;
	diagIter = DiagTerm->begin();
	uIter = UpperTerm->begin();
	lIter = LowerTerm->begin();
	N=(*vecIter);
	NP1=*(vecIter+1);
	NM1=0;

	(*vecIter) = (*diagIter)*N + (*uIter) * NP1;
	++uIter;++diagIter;	++vecIter;
	NM1=N;
	N=NP1;

	for( i=1; i<size-1; ++i, ++uIter, ++diagIter, ++lIter, ++vecIter)
	{
		NP1 = *(vecIter+1);
		(*vecIter) = (*lIter)*NM1 + (*diagIter)*N + (*uIter) * NP1;
		NM1 = N;
		N=NP1;
	}
	(*vecIter) = (*lIter)*NM1 + (*diagIter)*N;

#if defined(__GP_STRICT_VALIDATION)
	++vecIter;
	++diagIter;
	++lIter;
	if( diagIter != DiagTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "diagIter =! end" );
	if( uIter != UpperTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "upperIter =! end" );
	if( lIter != LowerTerm->end() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "lowerIter =! end" );
#endif
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FNumericalScheme
///	Routine: CheckCompatibilityWithModel
///	Returns: ARM_PricingStatesPtr
///	Action : Checks if the model is compatible with the 
///   current numerical scheme. Basically checks that model is 2F
////////////////////////////////////////////////////
bool ARM_PDE2FNumericalScheme::CheckCompatibilityWithModel( const ARM_PricingModel& model )
{
	return ( (model.FactorCount() == 2) ? true : 0);
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE2FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Initialiation of the NumericalScheme
///  builds nummethodstates and discretization schemes. 
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PDE2FNumericalScheme::Init( const ARM_PricingModel& model )
{
	ARM_GP_VectorPtr timeSteps( static_cast<ARM_GP_Vector*> (model.GetNumMethod()->GetTimeSteps()->Clone() ));
	model.NumMethodStateGlobalVariances( *timeSteps, itsGlobalVariances );

	size_t SpaceDiscretizationPointsNb = getSpaceDiscretizationPointsNb();
	SpaceDiscretizationPointsNb = (SpaceDiscretizationPointsNb%2 ) ? SpaceDiscretizationPointsNb : SpaceDiscretizationPointsNb+1;
	setSpaceDiscretizationPointsNb( SpaceDiscretizationPointsNb );
	size_t SpaceDiscretizationPointsNbSquare = SpaceDiscretizationPointsNb*SpaceDiscretizationPointsNb;


	//////////////////////////////////////////////////////////////////////////////////
	///// For testing purpose ONLY (Truncation)
	/// ARM_GP_T_Vector<size_t> TruncationVector(timeSteps->size());
	/// for (int k=0 ; k < timeSteps->size() ; ++k)
    /// TruncationVector[k] = ((SpaceDiscretizationPointsNb+(SpaceDiscretizationPointsNb*k)/timeSteps->size()) % 2) ? (SpaceDiscretizationPointsNb+(SpaceDiscretizationPointsNb*k)/timeSteps->size()) : (SpaceDiscretizationPointsNb+(SpaceDiscretizationPointsNb*k)/timeSteps->size()-1);
	/// setTruncator( new PDE2F_Truncator() );
	/// getTruncator()->setToTimeSizes( TruncationVector );
	//////////////////////////////////////////////////////////////////////////////////

	/// Price Will be in the middle of the Grid
	int PriceIndex = (SpaceDiscretizationPointsNbSquare-1)/2;
	setPriceIndex( PriceIndex );
	int midIndex = (SpaceDiscretizationPointsNb-1)/2;

	//////////////////////////////////////////////////////////////////////////////////
	///////////TEST (Truncation)
	// SpaceDiscretizationPointsNb*=2;
	// SpaceDiscretizationPointsNb-=1;
	// setSpaceDiscretizationPointsNb(SpaceDiscretizationPointsNb);
	// SpaceDiscretizationPointsNbSquare = SpaceDiscretizationPointsNb*SpaceDiscretizationPointsNb;
	// midIndex = (SpaceDiscretizationPointsNb-1)/2;
	//////////////////////////////////////////////////////////////////////////////////

	/// Computation of DeltaX
	double DiscretizationStepX = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/SpaceDiscretizationPointsNb;
	setSpaceDiscretizationStepX( DiscretizationStepX );
	double DiscretizationStepY = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/SpaceDiscretizationPointsNb;
	setSpaceDiscretizationStepY( DiscretizationStepY );

	/// ******** building NumMethodStates ********
	/// Building Pricing States ; its size is equal to the number of DiscreizationPoints
	ARM_PricingStatesPtr states( model.FirstPricingStates( SpaceDiscretizationPointsNbSquare  ) );
	
	/// Building NumMethodSates
	ARM_GP_Matrix * matrix = new ARM_GP_Matrix( 2, SpaceDiscretizationPointsNbSquare  );

	int i,j;
	int rowIndex = 0;
	double XValue = 0;
	for( i=0 ; i  <SpaceDiscretizationPointsNb ; ++i )
	{
		XValue = DiscretizationStepX*(i-midIndex);
		for( j=0 ; j<SpaceDiscretizationPointsNb ; ++j )
		{
			(*matrix)(0,rowIndex+j) = XValue;
			(*matrix)(1,rowIndex+j) = DiscretizationStepY*(j-midIndex);
		}
		rowIndex += SpaceDiscretizationPointsNb;
	}

	states->SetNumMethodStates(	ARM_GP_MatrixPtr(matrix) );

	/// build transition matrixes
	BuildTransitionMatrixes( model, timeSteps , states );

	/// Returns the first pricing states
	return states;
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: constructor
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

ARM_PDE3FNumericalScheme::ARM_PDE3FNumericalScheme(
size_t NX,
size_t NY,
size_t NZ,
GridType gridType,
const ARM_GP_Matrix& gridData,
const ARM_GP_Matrix& schedulerData) 
: ARM_PDENumericalScheme( 0, 0, NULL  ),
itsNX(NX),
itsNY(NY),
itsNZ(NZ),
itsGridType(gridType),
itsGridData(gridData),
itsSchedulerData(schedulerData),
itsDeltaX(NULL),
itsDeltaY(NULL),
itsDeltaZ(NULL),
itsXGrid(NULL),
itsYGrid(NULL),
itsZGrid(NULL),
itsTimeSteps(NULL),
itsModel(NULL)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: copy constructor
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

ARM_PDE3FNumericalScheme::ARM_PDE3FNumericalScheme(
const ARM_PDE3FNumericalScheme& rhs) 
: ARM_PDENumericalScheme(rhs),
itsNX(rhs.itsNX),
itsNY(rhs.itsNY),
itsNZ(rhs.itsNZ),
itsDeltaX(rhs.itsDeltaX),
itsDeltaY(rhs.itsDeltaY),
itsDeltaZ(rhs.itsDeltaZ),
itsTimeSteps(rhs.itsTimeSteps),
itsGridType(rhs.itsGridType),
itsGridData(rhs.itsGridData),
itsSchedulerData(rhs.itsSchedulerData),
itsModel(rhs.itsModel)
{
}


////////////////////////////////////////////////////
///	Class  : ARM _PDE3FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Initialiation of the NumericalScheme
///  builds nummethodstates and discretization schemes. 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_PDE3FNumericalScheme::Init( const ARM_PricingModel& model, double lambda )
{
	//TimeStep
	itsModel = &model;
	ARM_GP_VectorPtr timeSteps( static_cast<ARM_GP_Vector*> (model.GetNumMethod()->GetTimeSteps()->Clone() ));
	
	if(itsSchedulerData.rows()>1)
	{	
		if(itsSchedulerData.Elt(1,0)==1)//TimeStepPerYearScheduler
		{	
			int timeStepPerYear = itsSchedulerData.Elt(0,1);
			ARM_TimeStepPerYearScheduler timeStepPerYearSchedul(timeStepPerYear);
			timeSteps=timeStepPerYearSchedul.ComputeTimeSteps(model);	
		}
		else if(itsSchedulerData.Elt(1,0)==2)//ConstantVarianceScheduler
		{		
			size_t minStepBefore1stEventDate=itsSchedulerData.Elt(0,2);
			double minStdDev =	itsSchedulerData.Elt(1,2);

			ARM_ConstantVarianceMeanRevertingScheduler cstVarMeanRevertingSchedul(minStepBefore1stEventDate,minStdDev);
			timeSteps=cstVarMeanRevertingSchedul.ComputeTimeSteps(model);	
		}
		else if(itsSchedulerData.Elt(1,0)==3)//ConstantVarianceMeanRevertingScheduler
		{	
			size_t minStepBefore1stEventDate        = itsSchedulerData.Elt(0,3);
			size_t maxStepBefore1stEventDate        = itsSchedulerData.Elt(1,3);
			size_t stepNbPerYearBefore1stEventDate  = itsSchedulerData.Elt(2,3);
			double optimalYearFracTime              = itsSchedulerData.Elt(3,3);
			size_t stepNbPerYearBeforeOptimalTime   = itsSchedulerData.Elt(4,3);
			size_t stepNbPerYearAfterOptimalTime    = itsSchedulerData.Elt(5,3);	

			ARM_MultiRegimeScheduler multiRegimeSchedul( 
				minStepBefore1stEventDate,
				maxStepBefore1stEventDate,
				stepNbPerYearBefore1stEventDate,
				optimalYearFracTime,
				stepNbPerYearBeforeOptimalTime,
				stepNbPerYearAfterOptimalTime);

			timeSteps=multiRegimeSchedul.ComputeTimeSteps(model);	
		}
		else
		{
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[1,0] of Scheduler Data must be equals to 1 or 2 or 3");
		}
	}
	itsTimeSteps = ARM_GP_VectorPtr(timeSteps);

	//update of the timestep of the nummethod
	model.GetNumMethod()->SetTimeSteps(*timeSteps);
	model.GetNumMethod()->SetLastTimeIdx(timeSteps->size()-1);

	//XGrid
	size_t XGridDiscretizationPointsNb = GetNX();
	XGridDiscretizationPointsNb = (XGridDiscretizationPointsNb%2 ) ? XGridDiscretizationPointsNb : XGridDiscretizationPointsNb+1;
	//particular case for degnenerate: to be sure it works,the minimum of point in the x direction because of the discounting which is a function of the first factor
	if (XGridDiscretizationPointsNb <5)
	{
		XGridDiscretizationPointsNb =5;
	}
	SetNX( XGridDiscretizationPointsNb  );
	//YGrid
	size_t YGridDiscretizationPointsNb = GetNY();
	YGridDiscretizationPointsNb = (YGridDiscretizationPointsNb%2 ) ? YGridDiscretizationPointsNb : YGridDiscretizationPointsNb+1;
	SetNY( YGridDiscretizationPointsNb );
	//ZGrid
	size_t ZGridDiscretizationPointsNb = GetNZ();
	ZGridDiscretizationPointsNb = (ZGridDiscretizationPointsNb%2 ) ? ZGridDiscretizationPointsNb : ZGridDiscretizationPointsNb+1;
	SetNZ( ZGridDiscretizationPointsNb );

	//TotalGrid
	size_t SpaceDiscretizationPointsNb = XGridDiscretizationPointsNb*YGridDiscretizationPointsNb*ZGridDiscretizationPointsNb;

	/// Price Will be in the middle of the Grid
	int PriceIndex = (SpaceDiscretizationPointsNb-1)/2;
	setPriceIndex( PriceIndex );
	/// MedianIndex for each direction
	int XmidIndex = (XGridDiscretizationPointsNb-1)/2;
	int YmidIndex = (YGridDiscretizationPointsNb-1)/2;
	int ZmidIndex = (ZGridDiscretizationPointsNb-1)/2;

	double DiscretizationStepX; 
	double DiscretizationStepY;
	double DiscretizationStepZ;

		/// ******** building NumMethodStates ********
	/// Building Pricing States ; its size is equal to the number of DiscretizationPoints
	ARM_PricingStatesPtr states( model.FirstPricingStates( SpaceDiscretizationPointsNb  ) );

	/// Building NumMethodSates
	ARM_GP_Matrix * matrix = new ARM_GP_Matrix( 3, SpaceDiscretizationPointsNb  );

	int i,j,k,I,K;
	//size_t YZGridDiscretizationPointsNb=YGridDiscretizationPointsNb*ZGridDiscretizationPointsNb;
	size_t YXGridDiscretizationPointsNb=YGridDiscretizationPointsNb*XGridDiscretizationPointsNb;
	//vectors for the setting of the Grids (using the hypothesis of constant step)
	ARM_GP_VectorPtr XGrid( new ARM_GP_Vector(SpaceDiscretizationPointsNb) );
	ARM_GP_VectorPtr YGrid( new ARM_GP_Vector(SpaceDiscretizationPointsNb) );
	ARM_GP_VectorPtr ZGrid( new ARM_GP_Vector(SpaceDiscretizationPointsNb) );
	
	switch(itsGridType)
	{
		case StdDev:
		{
			//ARM_PricingStatesPtr states( NULL );
			itsGlobalVariances = ARM_MatrixVector(0);
			ARM_MatrixVector& localVariances  = ARM_MatrixVector(0);
			model.NumMethodStateLocalGlobalVariances( *timeSteps, localVariances, itsGlobalVariances );
			//Degenrescences
			//Denegerescence 1F
			if(itsGlobalVariances[timeSteps->size()-1]->rows()==1)
			{
				DiscretizationStepX = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/(XGridDiscretizationPointsNb-1);
				DiscretizationStepY= 0.001;
				DiscretizationStepZ= 0.001;
			}
			//Denegerescence 2F
			else if(itsGlobalVariances[timeSteps->size()-1]->rows()==2)
			{
				DiscretizationStepX = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/(XGridDiscretizationPointsNb-1);
				DiscretizationStepY = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/(YGridDiscretizationPointsNb-1);
				DiscretizationStepZ= 0.001;
			}
			else
			{
				DiscretizationStepX = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/(XGridDiscretizationPointsNb-1);
				DiscretizationStepY = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/(YGridDiscretizationPointsNb-1);
				DiscretizationStepZ = getStdDevNb()*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(2,2))/(ZGridDiscretizationPointsNb-1);
			}
			//remplissage des vecteurs 
			/// vectors delta
			ARM_GP_VectorPtr deltaX(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepX));
			ARM_GP_VectorPtr deltaY(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepY));
			ARM_GP_VectorPtr deltaZ(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepZ));
			///Setting the delta
			SetDeltaX( deltaX );
			SetDeltaY( deltaY );
			SetDeltaZ( deltaZ );
			if(abs(lambda)<3)//Qmapping, LogFxSpot
			{
				for( K=0 ; K<SpaceDiscretizationPointsNb ; ++K )
				{
					k=K/YXGridDiscretizationPointsNb;
					I=K%YXGridDiscretizationPointsNb;
					j=I/XGridDiscretizationPointsNb;
					i=I%XGridDiscretizationPointsNb;
					(*matrix)(0,K) = DiscretizationStepX*(i-XmidIndex);
					(*matrix)(1,K) = DiscretizationStepY*(j-YmidIndex);
					(*matrix)(2,K) = DiscretizationStepZ*(k-ZmidIndex);
					//a verifier
					(*XGrid)[K]=(*matrix)(0,K);
					(*YGrid)[K]=(*matrix)(1,K);
					(*ZGrid)[K]=(*matrix)(2,K);
				}
			}
			if(abs(lambda)==3)//SpotFX
			{
				for( K=0 ; K<SpaceDiscretizationPointsNb ; ++K )
				{
					k=K/YXGridDiscretizationPointsNb;
					I=K%YXGridDiscretizationPointsNb;
					j=I/XGridDiscretizationPointsNb;
					i=I%XGridDiscretizationPointsNb;
					(*matrix)(0,K) = DiscretizationStepX*(i-XmidIndex);
					(*matrix)(1,K) = DiscretizationStepY*(j-YmidIndex);
					(*matrix)(2,K) = (1+DiscretizationStepZ*(k-ZmidIndex));
					//a verifier
					(*XGrid)[K]=(*matrix)(0,K);
					(*YGrid)[K]=(*matrix)(1,K);
					(*ZGrid)[K]=(*matrix)(2,K);
				}
			}
		}
		break;
		case Fixed:
		{
			/// vectors delta
			ARM_GP_VectorPtr deltaX(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepX));
			ARM_GP_VectorPtr deltaY(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepY));
			ARM_GP_VectorPtr deltaZ(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepZ));
			///Setting the delta
			SetDeltaX( deltaX );
			SetDeltaY( deltaY );
			SetDeltaZ( deltaZ );
			//Degenrescences
			//Denegerescence 1F
			if(itsGridData.cols()==1)
			{
				DiscretizationStepX = 2*itsGridData.Elt(0,0)/(XGridDiscretizationPointsNb-1);
				DiscretizationStepY= 0.001;
				DiscretizationStepZ= 0.001;
			}
			//Denegerescence 2F
			else if(itsGridData.cols()==2)
			{
				DiscretizationStepX = 2*itsGridData.Elt(0,0)/(XGridDiscretizationPointsNb-1);
				DiscretizationStepY = 2*itsGridData.Elt(0,1)/(YGridDiscretizationPointsNb-1);
				DiscretizationStepZ= 0.001;
			}
			else
			{
				DiscretizationStepX = 2*itsGridData.Elt(0,0)/(XGridDiscretizationPointsNb-1);
				DiscretizationStepY = 2*itsGridData.Elt(0,1)/(YGridDiscretizationPointsNb-1);
				DiscretizationStepZ = 2*itsGridData.Elt(0,2)/(ZGridDiscretizationPointsNb-1);
			}
			//remplissage des vecteurs 
			if(abs(lambda)<3)//Qmapping, LogFxSpot
			{
				for( K=0 ; K<SpaceDiscretizationPointsNb ; ++K )
				{
					k=K/YXGridDiscretizationPointsNb;
					I=K%YXGridDiscretizationPointsNb;
					j=I/XGridDiscretizationPointsNb;
					i=I%XGridDiscretizationPointsNb;
					(*matrix)(0,K) = DiscretizationStepX*(i-XmidIndex);
					(*matrix)(1,K) = DiscretizationStepY*(j-YmidIndex);
					(*matrix)(2,K) = DiscretizationStepZ*(k-ZmidIndex);
					//a verifier
					(*XGrid)[K]=(*matrix)(0,K);
					(*YGrid)[K]=(*matrix)(1,K);
					(*ZGrid)[K]=(*matrix)(2,K);
				}
			}
			if(abs(lambda)==3)//SpotFX
			{
				for( K=0 ; K<SpaceDiscretizationPointsNb ; ++K )
				{
					k=K/YXGridDiscretizationPointsNb;
					I=K%YXGridDiscretizationPointsNb;
					j=I/XGridDiscretizationPointsNb;
					i=I%XGridDiscretizationPointsNb;
					(*matrix)(0,K) = DiscretizationStepX*(i-XmidIndex);
					(*matrix)(1,K) = DiscretizationStepY*(j-YmidIndex);
					(*matrix)(2,K) = (1+DiscretizationStepZ*(k-ZmidIndex));
					//a verifier
					(*XGrid)[K]=(*matrix)(0,K);
					(*YGrid)[K]=(*matrix)(1,K);
					(*ZGrid)[K]=(*matrix)(2,K);
				}
			}
		}
		break;
		case BiReg:
		{
			//ARM_PricingStatesPtr states( NULL );
			itsGlobalVariances = ARM_MatrixVector(0);
			ARM_MatrixVector& localVariances  = ARM_MatrixVector(0);
			model.NumMethodStateLocalGlobalVariances( *timeSteps, localVariances, itsGlobalVariances );
			//alpha, alphaX is the extreme standard dev in the X direction
			double alphaX=itsGridData.Elt(0,0);
			double alphaY=itsGridData.Elt(0,1);
			double alphaZ=itsGridData.Elt(0,2);
			//beta, betaX is the intern standard dev in the X direction, must be < alphaX
			double betaX=itsGridData.Elt(1,0);
			double betaY=itsGridData.Elt(1,1);
			double betaZ=itsGridData.Elt(1,2);
			if( betaX>alphaX)
			{
			    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[1,0] of the Grid Data must be inferior to Elt[0,0] ");
			}
			if( betaY>alphaY)
			{
			    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[1,1] of the Grid Data must be inferior to Elt[0,1] ");
			}
			if( betaZ>alphaZ)
			{
			    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[1,2] of the Grid Data must be inferior to Elt[0,2] ");
			}
			//Nint, NintX is the number of point inn the interior grid, 0<NintX< NtotX
			int NintX=itsGridData.Elt(2,0);
			int NintY=itsGridData.Elt(2,1);
			int NintZ=itsGridData.Elt(2,2);

			NintX = (NintX%2 ) ? NintX : NintX+1;
			NintY = (NintY%2 ) ? NintY : NintY+1;
			NintZ = (NintZ%2 ) ? NintZ : NintZ+1;

			/// intermediaire Index for each direction
			int XintIndex=-1;
			int YintIndex=-1;
			int ZintIndex=-1;
			if ( (NintX>=0)&&(NintX<=XGridDiscretizationPointsNb) )
			{
				XintIndex = XmidIndex-(NintX-1)/2;
			}
			else
			{
			    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[2,0] of the Grid Data must be an integer between 0 and NX ");
			}
			if ( (NintY>=0)&&(NintY<=YGridDiscretizationPointsNb) )
			{
				YintIndex = YmidIndex-(NintY-1)/2;
			}
			else
			{
			    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[2,1] of the Grid Data must be an integer between 0 and NXY ");
			}
			if ( (NintZ>=0)&&(NintZ<=ZGridDiscretizationPointsNb) )
			{
				ZintIndex = ZmidIndex-(NintZ-1)/2;
			}
			else
			{
			    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Elt[2,2] of the Grid Data must be an integer between 0 and NZ ");
			}
			//calculation of the delta, delta1 is the delta of the intern grid, delta2 is the one of the extern grid
			double DeltaX1=0.0;
			double DeltaX2=0.0;
			double DeltaY1=0.0;
			double DeltaY2=0.0;
			double DeltaZ1=0.0;
			double DeltaZ2=0.0;

			double denomDeltaX2=1.0;
			double denomDeltaY2=1.0;
			double denomDeltaZ2=1.0;

			if(XGridDiscretizationPointsNb-NintX>0)
			{
				denomDeltaX2=XGridDiscretizationPointsNb-NintX;
			}
			if(YGridDiscretizationPointsNb-NintY>0)
			{
				denomDeltaY2=YGridDiscretizationPointsNb-NintY;
			}
			if(ZGridDiscretizationPointsNb-NintZ>0)
			{
				denomDeltaZ2=ZGridDiscretizationPointsNb-NintZ;
			}


			//Degenrescences
			//Denegerescence 1F
			if(itsGlobalVariances[timeSteps->size()-1]->rows()==1)
			{
				DeltaX1 = betaX*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/(NintX-1);
				DeltaX2 = (alphaX-betaX)*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/denomDeltaX2;
				DeltaY1= 0.001;
				DeltaY2= 0.001;
				DeltaZ1= 0.001;
				DeltaZ2= 0.001;
			}
			//Degenerescence 2F
			else if(itsGlobalVariances[timeSteps->size()-1]->rows()==2)
			{
				DeltaX1 = betaX*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/(NintX-1);
				DeltaX2 = (alphaX-betaX)*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/denomDeltaX2;
				DeltaY1 = betaY*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/(NintY-1);
				DeltaY2 = (alphaY-betaY)*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/denomDeltaY2;
				DeltaZ1= 0.001;
				DeltaZ2= 0.001;
			}
			else
			{
				DeltaX1 = betaX*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/(NintX-1);
				DeltaX2 = (alphaX-betaX)*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(0,0))/denomDeltaX2;
				DeltaY1 = betaY*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/(NintY-1);
				DeltaY2 = (alphaY-betaY)*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(1,1))/denomDeltaY2;
				DeltaZ1 = betaZ*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(2,2))/(NintZ-1);
				DeltaZ2 = (alphaZ-betaZ)*sqrt((*(itsGlobalVariances[timeSteps->size()-1]))(2,2))/denomDeltaZ2;
			}
			//remplissage des vecteurs 
			/// vectors delta
			ARM_GP_VectorPtr deltaX(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepX));
			ARM_GP_VectorPtr deltaY(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepY));
			ARM_GP_VectorPtr deltaZ(new ARM_GP_Vector(SpaceDiscretizationPointsNb,DiscretizationStepZ));
			if(abs(lambda)<3)//Qmapping, LogFxSpot
			{
				for( K=0 ; K<SpaceDiscretizationPointsNb ; ++K )
				{
					k=K/YXGridDiscretizationPointsNb;
					I=K%YXGridDiscretizationPointsNb;
					j=I/XGridDiscretizationPointsNb;
					i=I%XGridDiscretizationPointsNb;
					//construction de deltaX et Xgrid 
					if( i<XintIndex )
					{
						(*deltaX)[K]=DeltaX2;
						(*matrix)(0,K) = DeltaX2*(i-XintIndex)+DeltaX1*(XintIndex-XmidIndex);
						(*XGrid)[K]=(*matrix)(0,K);
					}
					if( (i>=XintIndex)&&(i<=XGridDiscretizationPointsNb-1-XintIndex) )
					{
						(*deltaX)[K]=DeltaX1;
						(*matrix)(0,K) = DeltaX1*(i-XmidIndex);
						(*XGrid)[K]=(*matrix)(0,K);
					}
					if( (i>XGridDiscretizationPointsNb-1-XintIndex) )
					{
						(*deltaX)[K]=DeltaX2;
						(*matrix)(0,K) = DeltaX2*(i-(XGridDiscretizationPointsNb-1-XintIndex))-DeltaX1*(XintIndex-XmidIndex);
						(*XGrid)[K]=(*matrix)(0,K);
					}
					//construction de deltaY et Ygrid
					if( j<YintIndex )
					{
						(*deltaY)[K]=DeltaY2;
						(*matrix)(1,K) = DeltaY2*(j-YintIndex)+DeltaY1*(YintIndex-YmidIndex);
						(*YGrid)[K]=(*matrix)(1,K);
					}
					if( (j>=YintIndex)&&(j<=YGridDiscretizationPointsNb-1-YintIndex) )
					{
						(*deltaY)[K]=DeltaY1;
						(*matrix)(1,K) = DeltaY1*(j-YmidIndex);
						(*YGrid)[K]=(*matrix)(1,K);
					}
					if( j>YGridDiscretizationPointsNb-1-YintIndex )
					{
						(*deltaY)[K]=DeltaY2;
						(*matrix)(1,K) = DeltaY2*(j-(YGridDiscretizationPointsNb-1-YintIndex))-DeltaY1*(YintIndex-YmidIndex);
						(*YGrid)[K]=(*matrix)(1,K);
					}
					//construction de deltaZ et Zgrid
					if( k<ZintIndex ) 
					{
						(*deltaZ)[K]=DeltaZ2;
						(*matrix)(2,K) = DeltaZ2*(k-ZintIndex)+DeltaZ1*(ZintIndex-ZmidIndex);
						(*ZGrid)[K]=(*matrix)(2,K);
					}
					if( (k>=ZintIndex)&&(k<=ZGridDiscretizationPointsNb-1-ZintIndex) )
					{
						(*deltaZ)[K]=DeltaZ1;
						(*matrix)(2,K) = DeltaZ1*(k-ZmidIndex);
						(*ZGrid)[K]=(*matrix)(2,K);
					}
					if( (k>ZGridDiscretizationPointsNb-1-ZintIndex) )
					{
						(*deltaZ)[K]=DeltaZ2;
						(*matrix)(2,K) = DeltaZ2*(k-(ZGridDiscretizationPointsNb-1-ZintIndex))-DeltaZ1*(ZintIndex-ZmidIndex);
						(*ZGrid)[K]=(*matrix)(2,K);
					}
				}
			}
			///Setting the delta
			SetDeltaX( deltaX );
			SetDeltaY( deltaY );
			SetDeltaZ( deltaZ );
		}
		break;
	}
		///Setting of the Grids
	SetXGrid(XGrid);
	SetYGrid(YGrid);
	SetZGrid(ZGrid);

	//return of the states
	states->SetNumMethodStates(	ARM_GP_MatrixPtr(matrix) );

	/// Returns the first pricing states
	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Checks if the model is compatible with the 
/// current numerical scheme. Basically checks that model 
/// less or equal than 3F
////////////////////////////////////////////////////

bool ARM_PDE3FNumericalScheme::CheckCompatibilityWithModel( const ARM_PricingModel& model )
{
	return ( (model.FactorCount() <= 3) ? true : 0);

}

CC_END_NAMESPACE()
