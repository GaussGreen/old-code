/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarkovFunctional.cpp
 *
 *  \brief Markov Functional Model 1 Factor
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date August 2005
 */


#include "gpmodels/MFDFMap.h"
#include "gpinfra/pricingstates.h"
#include "gpbase/mempoolmatrix.h"
#include "gpbase/cloneutilityfunc.h"
#include <math.h>

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: constuctor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_DiscountFactorMap::ARM_DiscountFactorMap()
:	itsNumMethodStates(0),
	itsDiscountFactors(0),
	itsResetTimes(NULL),
	itsStartTimes(NULL),
	itsEndTimes(NULL) 
{
}

////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_DiscountFactorMap::ARM_DiscountFactorMap( const ARM_DiscountFactorMap& rhs) 
:	itsResetTimes	(CreateClonedPtr(&*rhs.itsResetTimes)),
	itsStartTimes	(CreateClonedPtr(&*rhs.itsStartTimes)),
	itsEndTimes		(CreateClonedPtr(&*rhs.itsEndTimes))
{
	DuplicateCloneablePtrVectorInPlace<ARM_GP_Matrix> (rhs.itsNumMethodStates, itsNumMethodStates);
	DuplicateCloneablePtrVectorInPlace<ARM_GP_Matrix> (rhs.itsDiscountFactors, itsDiscountFactors);
}



////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: IdxFromTime
///	Returns: size_t
///	Action : returns the position of date into timeVector
///  throw exception if not found
////////////////////////////////////////////////////
size_t ARM_DiscountFactorMap::IdxFromTime( const ARM_GP_VectorPtr& timeVector, double time ) const
{
	size_t k=0;
	while( k<timeVector->size() && timeVector->Elt(k) < time )
		++k;

	if( k >= timeVector->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": date not found in the map!");

	return k;
}

////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: NextIdxFromTime
///	Returns: size_t
///	Action : returns the position of date into timeVector
///  throws no exception if not found
////////////////////////////////////////////////////
size_t ARM_DiscountFactorMap::NextIdxFromTime( const ARM_GP_VectorPtr& timeVector, double time ) const
{
	size_t k=0;
	while( k<timeVector->size() && timeVector->Elt(k) < time )
		++k;

	return k;
}


////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: existsInVector
///	Returns: bool
///	Action : tells if time is in timeVector
////////////////////////////////////////////////////
bool ARM_DiscountFactorMap::existsInVector( const ARM_GP_VectorPtr& timeVector, double time ) const
{
	size_t k=0;
	while( k<timeVector->size() && timeVector->Elt(k) < time )
		++k;

	return (k < timeVector->size() && abs( time - timeVector->Elt(k) ) < K_DOUBLE_TOL ) ? true : false;
}



////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: StoreDiscountFactors
///	Returns: void
///	Action : stores discount factors and nummethodstates
////////////////////////////////////////////////////
void ARM_DiscountFactorMap::StoreDiscountFactors( double StorageTime, const ARM_PricingStatesPtr& states  )
{
	/// Get position in matrixes for storage
	size_t storageIndex = IdxFromTime( itsResetTimes, StorageTime );
	/// Stores Nummethod States
	itsNumMethodStates[storageIndex] = states->GetNumMethodStates(); //ARM_GP_MatrixPtr( static_cast<ARM_GP_Matrix*> (states->GetNumMethodStates()->Clone()) );

	/// Payoffs are copied and stored into a gp_matrixptr
	const ARM_MemPool_Matrix pf = states->GetPayoffs();
	size_t k,l, rowsNb, colsNb;	
	rowsNb = pf.GetRowsNb();
	colsNb = pf.GetColsNb();
	ARM_GP_MatrixPtr discountFactors(new ARM_GP_Matrix( rowsNb, colsNb ));
	ARM_GP_Matrix::iterator storageIter = discountFactors->begin();
	ARM_MemPool_Matrix::iterator statesIterator;

	for(k=0 ; k<rowsNb; ++k)
	{
		statesIterator = states->payoffsBeginIterator(k);
		for( l=0 ; l<colsNb ; ++l, ++statesIterator, ++storageIter)
			(*storageIter)=(*statesIterator);
	}

	/// Stores discount factors
	itsDiscountFactors[storageIndex] = discountFactors;
}


////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: DiscountFactorSorted
///	Returns: ARM_GP_VectorPtr
///	Action : returns DFs if nummethod states are sorted
///  useful for instance if pricing is made with another discretization
///  than calibration; interpolation made on DFs
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DiscountFactorMap::DiscountFactorSorted( double ResetTime, double EndTime, const ARM_PricingStatesPtr& states ) const
{
	ARM_GP_VectorPtr dfs = ReturnCalibrationDiscountFactors(ResetTime, EndTime);
	size_t storageIndex  = NDaysIdxFromTime( itsResetTimes, ResetTime, 1. );
	ARM_GP_MatrixPtr storageNumMethodStates = itsNumMethodStates[storageIndex];

	/// we need to ask for model states instead of nummethod states
	/// this is because of AMC (2ème passe met les num method states à 0)
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();
	size_t i,j;
	size_t storageSize  = storageNumMethodStates->cols();
	size_t newSize		= newNumMethodStates->cols();

	double curNumMethState, prevNumMethState, nextNumMethState;

	std::vector<double>& result = new std::vector<double>( newSize );
	
	///
	/// FLAT extrapol if new states are wider than storage states
	/// should send a warning ....
	i = 0;
	while (newNumMethodStates->Elt(0,i) <= storageNumMethodStates->Elt(0,0))
	{	
		result->Elt(i) = dfs->Elt(0); 
		i++;
	}
	size_t firstIdx = i;

	i = newSize - 1;
	while (newNumMethodStates->Elt(0,i) >= storageNumMethodStates->Elt(0,storageSize-1))
	{	
		result->Elt(i) = dfs->Elt(storageSize-1); 
		i--;
	}
	size_t lastIdx = i;
	
	////
	/// std linear interp
	j = 0;
	for( i=firstIdx; i<=lastIdx; ++i)
	{
		curNumMethState = newNumMethodStates->Elt(0,i);

		while (storageNumMethodStates->Elt(0,j) < curNumMethState)
			++j;
		
		nextNumMethState = storageNumMethodStates->Elt(0,j);
		prevNumMethState = storageNumMethodStates->Elt(0,j-1);

		result->Elt(i) = (  (curNumMethState - prevNumMethState) * dfs->Elt(j)
						  + (nextNumMethState - curNumMethState) * dfs->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
	}
	
	return ARM_GP_VectorPtr(result);

}

////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: DiscountFactorUnSorted
///	Returns: ARM_GP_VectorPtr
///	Action : returns DFs if nummethod states are not sorted
///  useful for pricing with Monte Carlo
///  --------    Could be enhanced using iterators ?
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DiscountFactorMap::DiscountFactorUnSorted( double ResetTime, double EndTime, const ARM_PricingStatesPtr& states ) const
{	
	ARM_GP_VectorPtr dfs = ReturnCalibrationDiscountFactors(ResetTime, EndTime);
	size_t storageIndex  = NDaysIdxFromTime( itsResetTimes, ResetTime, 1. );
	ARM_GP_MatrixPtr storageNumMethodStates = itsNumMethodStates[storageIndex];
	
	/// we need to ask for model states instead of nummethod states
	/// this is because of AMC (2ème passe met les num method states à 0)
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();
	size_t i, j;
	size_t storageSize  = storageNumMethodStates->cols();
	size_t newSize		= newNumMethodStates->cols();

	double curNumMethState, prevNumMethState, nextNumMethState;
	
	double storageXmin = storageNumMethodStates->Elt(0,0) ;
	double storageXmax = storageNumMethodStates->Elt(0,storageSize-1) ;
		
	double storageDx = (storageXmax - storageXmin) / (storageSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	const double tolerance = 1.e-8;

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(0,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  dfs->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  dfs->Elt(storageSize-1);
		/// general case
		else
		{
			/// CAUTION : we assume that x(i)- x(i-1) is constant
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			if (j>storageSize-1) j = storageSize-1;
			nextNumMethState = storageNumMethodStates->Elt(0,j);
			prevNumMethState = storageNumMethodStates->Elt(0,j-1);
			
			//// consistency test
			/// pas top pour les perfs
			if (curNumMethState>nextNumMethState + tolerance || curNumMethState<prevNumMethState - tolerance)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "DiscountFactorUnSorted: storage states don't seem to be regularly spaced");
						

			result->Elt(i) = (  (curNumMethState - prevNumMethState) * dfs->Elt(j)
							  + (nextNumMethState - curNumMethState) * dfs->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
		}
	}

	return ARM_GP_VectorPtr(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: CalibrationDiscountFactorRatiosInterpolate
///	Returns: ARM_GP_VectorPtr
///	Action : returns B(Treset, Tend)/B(Treset, Tn+1) interpolated on states
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DiscountFactorMap::DiscountFactorInterpolate( double ResetTime, double EndTime, const std::vector<double>& states, bool dividedByTerminalDf) const
{
	ARM_GP_VectorPtr dfs;
	
	if (dividedByTerminalDf)
		dfs = ReturnCalibrationDiscountFactorRatios(ResetTime, EndTime);
	else
		dfs = ReturnCalibrationDiscountFactors(ResetTime, EndTime);

	size_t storageIndex  = NDaysIdxFromTime( itsResetTimes, ResetTime, 1. );
	ARM_GP_MatrixPtr storageNumMethodStates = itsNumMethodStates[storageIndex];
	size_t i, j;
	size_t storageSize  = storageNumMethodStates->cols();
	size_t newSize		= states.size();

	double curNumMethState, prevNumMethState, nextNumMethState;
	
	double storageXmin = storageNumMethodStates->Elt(0,0) ;
	double storageXmax = storageNumMethodStates->Elt(0,storageSize-1) ;
	
	double storageDx = (storageXmax - storageXmin) / (storageSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = states[i];
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  dfs->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  dfs->Elt(storageSize-1);
		/// general case
		else
		{
			/// CAUTION : we assume that x(i)- x(i-1) is constant
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = storageNumMethodStates->Elt(0,j);
			prevNumMethState = storageNumMethodStates->Elt(0,j-1);
			
			//// consistency test
			/// pas top pour les perfs
			if (curNumMethState>nextNumMethState || curNumMethState<prevNumMethState)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "DiscountFactorUnSorted: storage states don't seem to be regularly spaced");
						

			result->Elt(i) = (  (curNumMethState - prevNumMethState) * dfs->Elt(j)
							  + (nextNumMethState - curNumMethState) * dfs->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
		}
	}

	return ARM_GP_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: DiscountFactor
///	Returns: ARM_GP_VectorPtr
///	Action : returns DFs; uses DFSorted or UnSorted depending on case
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DiscountFactorMap::DiscountFactor( double ResetTime, double EndTime, const ARM_PricingStatesPtr& states, bool statesAreSorted ) const
{
	/// ResetTime == EndTime
	if( abs ( ResetTime - EndTime ) < 1. )
		return ARM_GP_VectorPtr( new std::vector<double>( states->size(), 1. ) );

	if( !DoesResetTimeExist( ResetTime ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ResetTime not in schedule ; maybe this case will work later");
	
	if( statesAreSorted ) 
		return DiscountFactorSorted( ResetTime, EndTime, states );

	return DiscountFactorUnSorted( ResetTime, EndTime, states );
}


////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: ReturnCalibrationDiscountFactors
///	Returns: ARM_GP_VectorPtr
///	Action : returns DFs computed during calibration
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DiscountFactorMap::ReturnCalibrationDiscountFactors( double ResetTime, double EndTime ) const
{
	/// Sanity Checks
	if( !DoesResetTimeExist( ResetTime ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ResetTime not in schedule ; maybe this case will work later");
	if( EndTime > itsEndTimes->Elt(itsEndTimes->size() -1 ) + 30)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": EndTime too far beyond numeraire date");
	
	size_t storageIndex = IdxFromTime( itsResetTimes, ResetTime );
	ARM_GP_MatrixPtr discountFactors = itsDiscountFactors[storageIndex];

	/// ResetTime == EndTime
	if( abs ( ResetTime - EndTime ) < 1. )
		return ARM_GP_VectorPtr( new std::vector<double>( discountFactors->cols(), 1. ) );

	size_t ResetIdx = NDaysIdxFromTime( itsResetTimes, ResetTime, 1. );
	size_t EndIdx = 0;
	size_t Lastline = discountFactors->rows()-1;

	std::vector<double> * result;
	
	/// We want terminal zero coupon
	if( abs (EndTime - itsEndTimes->Elt(itsEndTimes->size() -1 ) ) < K_DOUBLE_TOL )
	{
		result = discountFactors->GetRow(Lastline);
		for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter )
			(*iter) = 1.0/(*iter);

		return ARM_GP_VectorPtr(result);
	}

	/// EndTime is in scheduler
	if( DoesStartTimeExist( EndTime ) )
	{
		EndIdx = NDaysIdxFromTime( itsStartTimes, EndTime, 2. );
		result = discountFactors->GetRow(Lastline-EndIdx+ResetIdx);
		ARM_GP_Matrix::iterator miter = discountFactors->begin() + Lastline * discountFactors->cols();
		for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter,++miter)
			(*iter) = (*iter)/(*miter);

		return ARM_GP_VectorPtr(result);
	}

	/// EndTime if not in scheduler
	/// -> interpolation (linear on zc rates)
	else
	{	
		if (EndTime > itsEndTimes->Elt(itsEndTimes->size()-1) ) 
			EndIdx = itsEndTimes->size()-1;
		else
			while( itsEndTimes->Elt(EndIdx) < EndTime )
				EndIdx ++;
			
		size_t ReverseEndIdx = ResetIdx + Lastline - EndIdx;

		/// case where EndTime is between reset and first start date
		if (itsResetTimes->Elt(ResetIdx) > itsStartTimes->Elt(EndIdx))
		{
			EndIdx ++;
			ReverseEndIdx --;
		}

		ARM_GP_Matrix::iterator termIter = discountFactors->begin() + Lastline*discountFactors->cols();
		
		result = discountFactors->GetRow(ReverseEndIdx);
		
		double nextWeight = (EndTime-itsStartTimes->Elt(EndIdx))/(itsEndTimes->Elt(EndIdx)-itsStartTimes->Elt(EndIdx)); 
		double prevWeight = 1.0 - nextWeight; 
		double prevZcRate, nextZcRate, currentZcRate;
				
		if (ReverseEndIdx != 0)
		{
			ARM_GP_Matrix::iterator nextIter = discountFactors->begin() + (ReverseEndIdx-1)*(discountFactors->cols());

			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter,++nextIter,++termIter)
			{			
				prevZcRate = - log( (*iter)/(*termIter)  );
				nextZcRate = - log( (*nextIter)/(*termIter) );
				currentZcRate = prevWeight * prevZcRate + nextWeight * nextZcRate;
				(*iter) = exp(-currentZcRate);
			}
		}
		else
		{			
			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter, ++termIter)
			{			
				prevZcRate = - log( (*iter)/(*termIter)  );
				nextZcRate = - log( 1./(*termIter) );
				currentZcRate = prevWeight * prevZcRate + nextWeight * nextZcRate;
				(*iter) = exp(-currentZcRate);
			}
		}


		return ARM_GP_VectorPtr(result);

		/***** old interpolation (linear rate)
		double delta = (EndTime-itsStartTimes->Elt(EndIdx))/(itsEndTimes->Elt(EndIdx)-itsStartTimes->Elt(EndIdx));

		if (ReverseEndIdx != 0)
		{
			ARM_GP_Matrix::iterator nextIter = discountFactors->begin() + (ReverseEndIdx-1)*(discountFactors->cols());

			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter,++nextIter,++termIter)
				(*iter) = (*iter)/((*termIter)*(1.+delta*((*iter)/(*nextIter)-1.)) );
		}
		else
		{
			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter, ++termIter)
				(*iter) = (*iter)/((*termIter)*(1.+delta*((*iter)-1.)) );
		}


		return ARM_GP_VectorPtr(result);
		******/
	}

}

////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: ReturnCalibrationDiscountFactorRatios
///	Returns: ARM_GP_VectorPtr
///	Action : returns DF Ratios computed during calibration
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DiscountFactorMap::ReturnCalibrationDiscountFactorRatios( double ResetTime, double EndTime ) const
{
	/// Sanity Checks
	if( !DoesResetTimeExist( ResetTime ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ResetTime not in schedule ; maybe this case will work later");
	if( EndTime > itsEndTimes->Elt(itsEndTimes->size() -1 ) + 30)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": EndTime too far beyond numeraire date");

	size_t storageIndex = IdxFromTime( itsResetTimes, ResetTime );
	ARM_GP_MatrixPtr discountFactors = itsDiscountFactors[storageIndex];

	/// ResetTime == EndTime
	if( abs ( ResetTime - EndTime ) < 1. )
		return ARM_GP_VectorPtr( new std::vector<double>( discountFactors->cols(), 1. ) );

	size_t ResetIdx = NDaysIdxFromTime( itsResetTimes, ResetTime, 1. );
	size_t EndIdx = 0;
	size_t Lastline = discountFactors->rows()-1;

	std::vector<double> * result;
	
	/// We want terminal zero coupon
	if( abs (EndTime - itsEndTimes->Elt(itsEndTimes->size() -1 ) ) < K_DOUBLE_TOL )
	{
		result = discountFactors->GetRow(Lastline);
		return ARM_GP_VectorPtr(new ARM_GP_Vector(discountFactors->cols(), 1.0));
	}

	/// EndTime is in scheduler
	if( DoesStartTimeExist( EndTime ) )
	{
		EndIdx = NDaysIdxFromTime( itsStartTimes, EndTime, 2. );
		result = discountFactors->GetRow(Lastline-EndIdx+ResetIdx);
		return ARM_GP_VectorPtr(result);
	}

	/// EndTime if not in scheduler
	/// interpolation (linear on zc rates)
	else
	{
		if (EndTime > itsEndTimes->Elt(itsEndTimes->size()-1) ) 
			EndIdx = itsEndTimes->size()-1;
		else
			while( itsEndTimes->Elt(EndIdx) < EndTime )
				EndIdx ++;
				
		size_t ReverseEndIdx = ResetIdx + Lastline - EndIdx;

		/// case where EndTime is between reset and first start date
		if (itsResetTimes->Elt(ResetIdx) > itsStartTimes->Elt(EndIdx))
		{
			EndIdx ++;
			ReverseEndIdx --;
		}

		result = discountFactors->GetRow(ReverseEndIdx);

		double nextWeight = (EndTime-itsStartTimes->Elt(EndIdx))/(itsEndTimes->Elt(EndIdx)-itsStartTimes->Elt(EndIdx)); 
		double prevWeight = 1.0 - nextWeight; 
		double prevZcRate, nextZcRate, currentZcRate;
				
		if (ReverseEndIdx != 0)
		{
			ARM_GP_Matrix::iterator nextIter = discountFactors->begin() + (ReverseEndIdx-1)*(discountFactors->cols());

			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter,++nextIter)
			{			
				prevZcRate = - log( (*iter) );
				nextZcRate = - log( (*nextIter) );
				currentZcRate = prevWeight * prevZcRate + nextWeight * nextZcRate;
				(*iter) = exp(-currentZcRate);
			}
		}
		else
		{			
			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter)
			{			
				prevZcRate = - log( (*iter) );
				nextZcRate = 0.0;
				currentZcRate = prevWeight * prevZcRate + nextWeight * nextZcRate;
				(*iter) = exp(-currentZcRate);
			}
		}

		return ARM_GP_VectorPtr(result);

		/***** old interpolation (linear rate)
		double delta = (EndTime-itsStartTimes->Elt(EndIdx))/(itsEndTimes->Elt(EndIdx)-itsStartTimes->Elt(EndIdx));

		if (ReverseEndIdx != 0)
		{
			ARM_GP_Matrix::iterator nextIter = discountFactors->begin() + (ReverseEndIdx-1)*(discountFactors->cols());

			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter,++nextIter)
				(*iter) = (*iter)/(1.+delta*((*iter)/(*nextIter)-1.)) ;
		}
		else
		{
			for( std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter)
				(*iter) = (*iter)/(1.+delta*((*iter)-1.)) ;
		}


		return ARM_GP_VectorPtr(result);
		*****/

	}

}


///////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: NDaysIdxFromTime
///	Returns: size_t
///	Action : 
////////////////////////////////////////////////////
size_t ARM_DiscountFactorMap::NDaysIdxFromTime( const ARM_GP_VectorPtr& timeVector, double time, double daysNb ) const
{
	int i=-1, N = timeVector->size();

	while( (abs( time - timeVector->Elt(++i) ) > daysNb ) && i<N ) {}

	if( i == N )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": date not found in the map!");

	return i;
}

///////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: isNDaysIdxFromTime
///	Returns: bool
///	Action : 
////////////////////////////////////////////////////
bool ARM_DiscountFactorMap::isNDaysIdxFromTime( const ARM_GP_VectorPtr& timeVector, double time, double daysNb ) const
{
	int i=-1, N = timeVector->size();

	while( (abs( time - timeVector->Elt(++i) ) > daysNb ) && i<N ) {}

	if( i == N )
		return false;

	return true;

}

///////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: toString
///	Returns: string
///	Action : 
////////////////////////////////////////////////////

string ARM_DiscountFactorMap::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_DiscountFactorMap\n";
    os << indent << "-----------------------\n";

	os << CC_NS(std,endl);
	
	size_t resetIndex = 0;
	double dfMax (0.0);
	double resetIndexMax;
	double maturityIndexMax;
	ARM_MatrixPtrVector::const_iterator iterMax;
	size_t stateIndexMax;
	
	for( ARM_MatrixPtrVector::const_iterator iter = itsDiscountFactors.begin() ; iter != itsDiscountFactors.end() ; ++iter, ++resetIndex )
	{				
		os << indent << CC_NS(std,fixed)    << CC_NS(std,setw)(0) << CC_NS(std,setprecision)(0) <<  "DFs stored for reset date : " << itsResetTimes->Elt(resetIndex) << CC_NS(std,endl);
		os << "-------------------------------------"<< CC_NS(std,endl);

		for(int j=0; j<(*iter)->GetColsNb(); ++j)
		{
			if (j%20==0)
			{
				os	<< "state #" << j << "\t\t";
			}
		}
		os << CC_NS(std,endl);

		for(size_t i=0; i<(*iter)->GetRowsNb(); ++i)
		{
			for(int j=0; j<(*iter)->GetColsNb(); ++j)
			{
				if ((*iter)->Elt(i, j)>dfMax)
				{
					dfMax			 = (*iter)->Elt(i, j);
					resetIndexMax	 = resetIndex;
					maturityIndexMax = i;
					iterMax			 = iter;
					stateIndexMax	 = j;
				}

				if (j%20==0)
				{
					if ((*iter)->Elt(i, j) < 1.e5)
						os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
							<< CC_NS(std,setw)(8) << (*iter)->Elt(i, j) << "\t\t";
					else
						os	<< CC_NS(std,scientific)<< CC_NS(std,setw)(4) 
							<< CC_NS(std,setprecision)(4) << (*iter)->Elt(i, j) << "\t\t";
				}
			}
			os << CC_NS(std,endl);
		}

		os << CC_NS(std,endl);
	}

	if (!itsResetTimes.IsNull())
	{
		os << "Max Discount Factor :"  << CC_NS(std,endl);
		os << "----------------------" << CC_NS(std,endl);

		os << CC_NS(std,fixed)     << CC_NS(std,setw)(0) << CC_NS(std,setprecision)(0) << "reset date \t#" << resetIndexMax	 << "\t("<< itsResetTimes->Elt(resetIndexMax) << ")"<< CC_NS(std,endl);
		os << CC_NS(std,fixed)     << CC_NS(std,setw)(0) << CC_NS(std,setprecision)(0) << "maturity \t#"   << maturityIndexMax << "\t("<< itsStartTimes->Elt(resetIndexMax+maturityIndexMax) << ")"<< CC_NS(std,endl);
		os << CC_NS(std,fixed)     << CC_NS(std,setw)(0) << CC_NS(std,setprecision)(0) << "state \t#"	   << stateIndexMax << CC_NS(std,endl);
		os << CC_NS(std,scientific)<< CC_NS(std,setw)(4) << CC_NS(std,setprecision)(4) << "DF = \t\t"	   << (*iterMax)->Elt(maturityIndexMax, stateIndexMax) <<  CC_NS(std,endl);
	}

    return os.str();	
}


CC_END_NAMESPACE()