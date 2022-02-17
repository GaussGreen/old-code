/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file schedulerbase.cpp
 *
 *  \brief scheduler object are used to discretise in time
 *
 *	\author  JM Prie, E Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */

#include "gpinfra/schedulerbase.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SchedulerBase
///	Routine: ComputeTimeSteps
///	Returns: ARM_GP_VectorPtr
///	Action : Default method for the computation of the time steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_SchedulerBase::ComputeTimeSteps(const ARM_PricingModel& model) const
{
    /// Build a default schedule with the required number of steps
    const ARM_GP_Vector& modelSchedule = *(model.GetNumMethod()->GetTimeSteps());
    int nbModelSteps=modelSchedule.size()-1; // [0]=0

    int i,nbSteps = model.GetNumMethod()->GetNbSteps();

    int nbToAdd=nbSteps-nbModelSteps;
    if(nbToAdd <= 0)
        return ARM_GP_VectorPtr(new ARM_GP_Vector(modelSchedule)); /// keep the model schedule

    ARM_GP_Vector* timeSteps = new ARM_GP_Vector(nbSteps+1);
    timeSteps[0]=0.0;
    int j,idx=1;

    int nbInSteps = nbToAdd/nbModelSteps;
    int nbExtraInSteps = nbToAdd - nbInSteps*nbModelSteps;
    int jmax;
    double lastModelTime=0.0,modelTime;
    double step,time;
    for(i=1;i<=nbModelSteps;++i)
    {
        modelTime = modelSchedule[i];
        jmax = (i-1<nbExtraInSteps ? nbInSteps+1 : nbInSteps);
        step = (modelTime-lastModelTime)/(jmax+1);
        for(j=0,time=lastModelTime+step;j<jmax;++j,time+=step,++idx)
            (*timeSteps)[idx]=time;
        (*timeSteps)[idx]=modelTime; // to avoid rounding pbs !
        ++idx;
        lastModelTime = modelTime;
    }
    return ARM_GP_VectorPtr(timeSteps);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/