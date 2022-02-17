/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MFDFMap.h
 *
 *  \brief 
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\time September 2005
 */

#ifndef _INGPMODELS_MFDF_H
#define _INGPMODELS_MFDF_H

#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/gpmatrix.h"
#include "gpinfra/pricingstates.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_DiscountFactorMap : public ARM_RootObject
{
private:
	ARM_GP_VectorPtr itsResetTimes;
	ARM_GP_VectorPtr itsStartTimes; // to be completed (not filled yet)
	ARM_GP_VectorPtr itsEndTimes;
	
	ARM_MatrixPtrVector itsNumMethodStates;
	ARM_MatrixPtrVector itsDiscountFactors;

	size_t IdxFromTime( const ARM_GP_VectorPtr& timeVector, double time ) const;
	size_t NDaysIdxFromTime( const ARM_GP_VectorPtr& timeVector, double time, double daysNb ) const;
	size_t NextIdxFromTime( const ARM_GP_VectorPtr& timeVector, double time ) const;
	bool existsInVector( const ARM_GP_VectorPtr& timeVector, double time ) const;
	bool isNDaysIdxFromTime( const ARM_GP_VectorPtr& timeVector, double time, double daysNb ) const;


	/// Two private DiscountFactor methods for both sorted and unsorted
	ARM_GP_VectorPtr DiscountFactorSorted( double ResetTime, double EndTime, const ARM_PricingStatesPtr& states ) const;
	ARM_GP_VectorPtr DiscountFactorUnSorted( double ResetTime, double EndTime, const ARM_PricingStatesPtr& states ) const;

public:
	/// Constructors/destructor
	ARM_DiscountFactorMap();
	ARM_DiscountFactorMap( const ARM_DiscountFactorMap& rhs) ;
	virtual ~ARM_DiscountFactorMap() {}
	
	ASSIGN_OPERATOR(ARM_DiscountFactorMap)

	/// Storage
	void StoreDiscountFactors( double StorageTime, const ARM_PricingStatesPtr& states ); 
	
	/// Calibration discount factors
	ARM_GP_VectorPtr ReturnCalibrationDiscountFactors( double EvalTime, double EndTime ) const;
	ARM_GP_VectorPtr ReturnCalibrationDiscountFactorRatios( double EvalTime, double EndTime ) const;
	ARM_GP_VectorPtr DiscountFactorInterpolate( double ResetTime, double EndTime, const ARM_GP_Vector& states, bool dividedByTerminalDf = true ) const;
	
	/// Forward discount factor and discountfactor
	ARM_GP_VectorPtr DiscountFactor( double EvalTime, double EndTime, const ARM_PricingStatesPtr& states, bool statesAreSorted = false ) const; 

	/// Bullshit
	inline double	DFRatioInterpolateFromStates( const ARM_PricingStatesPtr& states, double state, size_t startSearchIdx, size_t EndIdx ) const;
	//inline double LevelOnLastDfInterpolateFromStates( const ARM_PricingStatesPtr& states, double state, size_t startSearchIdx, size_t firstEndIdx, size_t nbPeriods, const ARM_GP_Vector& interestTerms) const;
	inline double LevelOnLastDfInterpolateFromStates( const ARM_PricingStatesPtr& states, double state, size_t startSearchIdx, size_t firstEndIdx, const ARM_IntVector& relIdx, const ARM_GP_Vector& interestTerms) const;
	/// Accessors
	inline void setResetTimes( const ARM_GP_VectorPtr& resetTimes ) { itsResetTimes = resetTimes; itsNumMethodStates.resize(resetTimes->size()); itsDiscountFactors.resize(resetTimes->size());}
	inline ARM_GP_VectorPtr getResetTimes() const { return itsResetTimes; }
	inline void setEndTimes( const ARM_GP_VectorPtr& endTimes ) { itsEndTimes = endTimes; }
	inline ARM_GP_VectorPtr getEndTimes() const { return itsEndTimes; }
	inline void setStartTimes( const ARM_GP_VectorPtr& startTimes ) { itsStartTimes = startTimes; }
	inline ARM_GP_VectorPtr getStartTimes() const { return itsStartTimes; }
	
	
	inline size_t getResetTimeIdx( double time ) const { return IdxFromTime( itsResetTimes, time ); }
	inline size_t getEndTimeIdx( double time ) const { return IdxFromTime( itsEndTimes, time ); }
	
	inline size_t dfsNbAtTime( double time ) const { return itsResetTimes->size() - getResetTimeIdx( time ); }

	inline size_t NDaysFromStartTimeIdx( double time, double daysNb ) const { return NDaysIdxFromTime( itsStartTimes, time, daysNb ); }
	inline size_t NDaysFromEndTimeIdx( double time, double daysNb ) const { return NDaysIdxFromTime( itsEndTimes, time, daysNb ); }
	inline size_t NDaysFromResetTimeIdx( double time, double daysNb ) const { return NDaysIdxFromTime( itsResetTimes, time, daysNb ); }

	inline double getTerminalTime() const { return itsEndTimes->Elt( itsEndTimes->size() -1 ); }
	inline double getLastResetTime() const { return itsResetTimes->Elt( itsResetTimes->size() -1 ); }
	inline double getResetTimesSize() const { return itsResetTimes->size(); }

	inline double getResetTime( size_t Idx ) const { return itsResetTimes->Elt(Idx); }
	inline double getStartTime( size_t Idx ) const { return itsStartTimes->Elt(Idx); }
	inline double getEndTime( size_t Idx ) const { return itsEndTimes->Elt(Idx); }
		
	inline double getNumMethStateMax( size_t StorageIdx ) const { return itsNumMethodStates[StorageIdx]->Elt(0,itsNumMethodStates[StorageIdx]->size()-1); }
	inline double getNumMethStateMin( size_t StorageIdx ) const { return itsNumMethodStates[StorageIdx]->Elt(0,0); }
	inline double getNumMethStateMax( double evalTime ) const 
	{
		size_t resetTimeIdx = getResetTimeIdx(evalTime);
		return getNumMethStateMax(resetTimeIdx);
	}
	inline double getNumMethStateMin( double evalTime ) const 
	{
		size_t resetTimeIdx = getResetTimeIdx(evalTime);
		return getNumMethStateMin(resetTimeIdx);
	}


	inline bool DoesResetTimeExist( double time ) const { return existsInVector( itsResetTimes, time );}
	inline bool DoesStartTimeExist( double time ) const { return existsInVector( itsStartTimes, time ); }
	inline bool DoesEndTimeExist( double time ) const { return existsInVector( itsEndTimes, time ); }

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_DiscountFactorMap(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};



////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: DFRatioInterpolateFromStates
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
inline double ARM_DiscountFactorMap::DFRatioInterpolateFromStates( const ARM_PricingStatesPtr& states, double state, size_t startSearchIdx, size_t EndIdx ) const
{
	ARM_GP_MatrixPtr NumMethodStates = states->GetNumMethodStates();

	size_t i=startSearchIdx, N = NumMethodStates->cols();
	while( NumMethodStates->Elt(0,i) < state && i < N-1 )
		++i;

	double MinState = state - NumMethodStates->Elt(0,i-1);
	double MaxState = NumMethodStates->Elt(0,i) - state;

	return (states->GetPayoff(i-1,EndIdx) * MaxState + states->GetPayoff(i,EndIdx) * MinState)/(MaxState+MinState);
}

////////////////////////////////////////////////////
///	Class  : ARM_DiscountFactorMap
///	Routine: LevelOnLastDfInterpolateFromStates
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_DiscountFactorMap::LevelOnLastDfInterpolateFromStates( const ARM_PricingStatesPtr& states, double state, size_t startSearchIdx, size_t firstEndIdx, const ARM_IntVector& relIdx, const ARM_GP_Vector& interestTerms) const
{
	ARM_GP_MatrixPtr NumMethodStates = states->GetNumMethodStates();

	size_t	i = startSearchIdx, 
			N = NumMethodStates->cols();
	
	while( NumMethodStates->Elt(0,i) < state && i < N-1 )
		++i;

	if (i==0) i++;

	double MinState = state - NumMethodStates->Elt(0,i-1);
	double MaxState = NumMethodStates->Elt(0,i) - state;
	
	/// result assumed to be of size nbPeriods
	double prevPayoff(0);
	double curPayoff(0);
	double delta;
			
	int	idx, 
		size = interestTerms.size();

	for (int k(0); k<size; k++)
	{
		delta = interestTerms[k];
		idx   = firstEndIdx - relIdx[k] ;
			
		if (idx >= 0)
		{
			prevPayoff += delta * states->GetPayoff(i-1,idx);
			curPayoff  += delta * states->GetPayoff(i,  idx);
		}
		/// if last rate date is terminal date
		/// df ratio = 1 (not stored in payoffs)
		else
		{
			prevPayoff += delta ;
			curPayoff  += delta ;
		}

	}

	return (prevPayoff * MaxState + curPayoff * MinState)/(MaxState+MinState);
}





CC_END_NAMESPACE()

#endif