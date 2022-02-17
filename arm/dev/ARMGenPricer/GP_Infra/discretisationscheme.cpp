/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file discretisationscheme.cpp
 *  \brief implements various discretisation scheme
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#include "gpinfra/discretisationscheme.h"

#include "gpinfra/pricingmodel.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/timeinfo.h"

#include "gpbase/ostringstream.h"
#include "gpbase/env.h"
#include "gpbase/vectormanip.h"

#include <memory>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_EventTime
///	Routine: ModelTimesFromTimesInfo
///	Returns: 
///	Action : returns zero plus the event time from time Info
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_EventTime::ModelTimesFromTimesInfo(
	const ARM_TimeInfoPtrVector& timeInfos,
	ARM_PricingModel&	model) const 
{
	size_t nbEvents=timeInfos.size();
	
#if defined(__GP_STRICT_VALIDATION)
	if( model.GetNumMethod() == ARM_NumMethodPtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		ARM_USERNAME + ": no numerical method set for the model!" );
#endif
	
	/// an ugly switch... this may be revisited to split the 
	/// numMethod in the type BackwardLooking, ForwardLooking and so on..
	/// and have derived class manage everything!
	ARM_NumMethod::GP_PricingDirection pricingDirection = model.GetNumMethod()->GetPricingDirection();
	switch(pricingDirection)
	{
	/// backward and forward looking add zero to the model times!
	case ARM_NumMethod::GP_BCKWDLOOKING:
	case ARM_NumMethod::GP_FWDLOOKING:
	case ARM_NumMethod::GP_FWDBCKWDLOOKING:
		{
			/// Initialise a basic schedule event dates
			/// Only eventTime >= 0 were passed in increasing order
			/// Don't insert 0 if already an event time
			ARM_GP_Vector* timeSteps;
			size_t i,j;

			if(timeInfos[0]->GetEventTime()< K_NEW_DOUBLE_TOL  && pricingDirection == ARM_NumMethod::GP_FWDLOOKING)  
				timeInfos[0]->SetEventTime(K_NEW_DOUBLE_TOL);
			
			if(timeInfos[0]->GetEventTime()>= K_NEW_DOUBLE_TOL )
			{
				timeSteps=new ARM_GP_Vector(nbEvents+1);
				(*timeSteps)[0]=0.0;
				j=1;
			}
			else
			{
				timeSteps=new ARM_GP_Vector(nbEvents);
				j=0;
			}

			for(i=0;i<nbEvents;++i,++j)
				(*timeSteps)[j]= timeInfos[i]->GetEventTime();

			return timeSteps;
		}
		
	default:
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": non supported numethod pricing direction " << ARM_NumMethod::GP_PricingDirectionTxt[pricingDirection];
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
}





////////////////////////////////////////////////////
///	Struct : ARM_EventAndModelTime
///	Routine: ModelTimesFromTimesInfo
///	Returns: 
///	Action : returns eventTime and model fixing
///			 if the model support analytical marginal
///		     there is no reason to insert model times
///			 and we only return the event times
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_EventAndModelTime::ModelTimesFromTimesInfo(
	const ARM_TimeInfoPtrVector& timeInfos,
	ARM_PricingModel& model ) const 
{
	/// first computes eventTimes
	CC_NS(std,auto_ptr)<ARM_GP_Vector> pEventTimes( ARM_EventTime().ModelTimesFromTimesInfo( timeInfos, model ) );
	
	/// get model times only if appropriate
	if( !model.SupportAnalyticMarginal() )
	{
		CC_NS(std,auto_ptr)<ARM_GP_Vector> pModelTimes( model.ComputeModelTimes( timeInfos ) );

		return MergeSortedVectorNoDuplicates( *pEventTimes, *pModelTimes );
	}
	else
		return VectorUnique( *pEventTimes );
}


////////////////////////////////////////////////////
///	Struct : ARM_FixStepTimeInserter
///	Routine: AddTimes
///	Returns: 
///	Action : returns times after inserting fix step times between two times
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_FixStepTimeInserter::AddTimes(const ARM_GP_Vector& initialTimes) const
{
	/// use vector<double> because of its dynamic growth property!
	vector<double> tmpResult;
	int initialTimesSize=initialTimes.size();	/// int because can be equal to zero potentially!
	tmpResult.reserve(initialTimesSize);
    
	double currentTime,nextTime; /// for better readibility
	int i,j;

	for(i=0;i<initialTimesSize-1;++i)
	{
		currentTime = initialTimes[i];
		nextTime	= initialTimes[i+1];
		
		tmpResult.push_back(currentTime);

		for(j=1;j<(nextTime-currentTime)/itsFixStep; ++j)
		{
			/// should we include the last model time
			if(currentTime+j*itsFixStep<nextTime - itsSensitivity )
				tmpResult.push_back(currentTime+j*itsFixStep);
		}
	}

	/// include the last time!
	tmpResult.push_back(nextTime);

	return new ARM_GP_Vector(tmpResult.size(),&tmpResult[0]);
}



////////////////////////////////////////////////////
///	Struct : ARM_FixPointNbTimeInserter
///	Routine: AddTimes
///	Returns: 
///	Action : returns times after inserting fix nb of points between two times
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_FixPointNbTimeInserter::AddTimes(const ARM_GP_Vector& initialTimes) const
{
	/// use vector<double> because of its dynamic growth property!
	vector<double> tmpResult;
	size_t initialTimesSize=initialTimes.size();
	tmpResult.reserve(initialTimesSize);
    
	double currentTime,nextTime; /// for better readibility
	size_t i,j;

	for(i=0;i<initialTimesSize-1;++i)
	{
		currentTime = initialTimes[i];
		nextTime	= initialTimes[i+1];
		
		int increment=(nextTime-currentTime)/itsFixPointNb;
		for(j=0;j<itsFixPointNb; ++j)
		/// should we include the last model time
		if(currentTime+j*increment<nextTime)
			tmpResult.push_back(currentTime+j*increment);
	}
	// include the last one!
	tmpResult.push_back(nextTime);

	return new ARM_GP_Vector(tmpResult.size(),&tmpResult[0]);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

