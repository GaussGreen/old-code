/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file schedulerfactory.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#include "gpnummethods/schedulerfactory.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"


/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/scheduler.h"
#include "gpnummethods/argconvdefault.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SchedulerFactoryData
///	Routine: toString
///	Returns: 
///	Action : Object dump
////////////////////////////////////////////////////
string ARM_SchedulerFactoryData::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

	os << indent << "Scheduler Factory Data\n";
	os << indent << "--------------------\n";
	os << indent << "\nScheduler=" << ARM_ArgConvReverse_SchedulerType.GetString(itsSchedulerType) << "\n";
    if(itsSchedulerDatas.size()>1)
        os << indent << "   NbSteps="<<itsSchedulerDatas[0] << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_SchedulerFactoryImp
///	Routine: CreateScheduler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SchedulerBase* ARM_SchedulerFactoryImp::CreateScheduler(const ARM_SchedulerFactoryData& schedulerFactoryData, int& nbSteps)
{
    return CreateScheduler(
		schedulerFactoryData.GetSchedulerType(), 
		schedulerFactoryData.GetSchedulerDatas(),
		nbSteps);
}

////////////////////////////////////////////////////
///	Class  : ARM_SchedulerFactoryImp
///	Routine: CreateScheduler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SchedulerBase* ARM_SchedulerFactoryImp::CreateScheduler(
        int schedulerType, 
		const ARM_GP_Vector& schedulerDatas,
		int& nbSteps)
{
    ARM_SchedulerBase* scheduler        = NULL;

    /// Default value
    nbSteps=25;
    size_t nbDatas = schedulerDatas.size();

    /// Scheduler creation
    switch(schedulerType)
    {
    case ARM_SchedulerBase::ConstantVariance:
        if(nbDatas>1)
        {
            nbSteps     = schedulerDatas[0];
            scheduler   = new ARM_ConstantVarianceScheduler(schedulerDatas[1]);
        }
        else
        {
            scheduler   = new ARM_ConstantVarianceScheduler;
            if(nbDatas>0)
                nbSteps = schedulerDatas[0];
        }
        break;

    case ARM_SchedulerBase::ConstantVarianceMeanReverting:
        if(nbDatas>2)
        {
            nbSteps     = schedulerDatas[0];
            scheduler   = new ARM_ConstantVarianceMeanRevertingScheduler(schedulerDatas[1],schedulerDatas[2]);
        }
        else
        {
            scheduler   = new ARM_ConstantVarianceMeanRevertingScheduler;
            if(nbDatas>0)
                nbSteps = schedulerDatas[0];
        }
        break;

    case ARM_SchedulerBase::MultiRegime:
        nbSteps = 0; // not set because computed by the multi regime scheduler itself
        if(nbDatas == 6)
            scheduler   = new ARM_MultiRegimeScheduler(schedulerDatas[0],
			schedulerDatas[1],
			schedulerDatas[2],
			schedulerDatas[3],
			schedulerDatas[4],
			schedulerDatas[5]);
        else if(nbDatas == 8)
            scheduler   = new ARM_MultiRegimeScheduler(schedulerDatas[0],
			schedulerDatas[1],
			schedulerDatas[2],
			schedulerDatas[3],
			schedulerDatas[4],
			schedulerDatas[5],
			schedulerDatas[6],	// second optimal year
			schedulerDatas[7]); // nb of steps before second optimal year*/
        else
            scheduler   = new ARM_MultiRegimeScheduler;
        break;

	case ARM_SchedulerBase::TimeStepPerYear:
		if(nbDatas>0)
			scheduler = new ARM_TimeStepPerYearScheduler(schedulerDatas[0]);
		else
			scheduler = new ARM_TimeStepPerYearScheduler;
		break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown scheduler type");
    }

	/// return the result
    return scheduler;
}


ARM_SingletonHolder<ARM_SchedulerFactoryImp> ARM_SchedulerFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/