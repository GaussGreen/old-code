/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file scheduler.cpp
 *
 *  \brief scheduler object are used to discretise in time
 *
 *	\author  JM Prie, E Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#include "gpnummethods/scheduler.h"

/// gpbase
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/discretisationscheme.h"

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
///	Class  : ARM_MultiRegimeScheduler
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_MultiRegimeScheduler
///	Routine: InsertCstSteps
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_MultiRegimeScheduler::InsertCstSteps(ARM_GP_Vector* timeSteps,
	int startIdx,
	int endIdx,
	const ARM_GP_Vector& modelSchedule,
	size_t nbStepPerYear) const
{
    double time,eventTime = modelSchedule[startIdx];
    double step,nextEventTime;
    int i,idx,nbSteps;
    for(idx=startIdx;idx<endIdx;++idx)
    {
        nextEventTime = modelSchedule[idx+1];
        nbSteps=static_cast<int>( floor((nextEventTime-eventTime)/K_YEAR_LEN*nbStepPerYear+0.5) );
        if(nbSteps>0)
        {
            step = (nextEventTime-eventTime)/nbSteps;
            for(i=0,time=eventTime+step;i<nbSteps-1;++i,time += step)
                timeSteps->push_back(time);
        }
        timeSteps->push_back(nextEventTime); // to avoid rounding pbs !
        eventTime = nextEventTime;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiRegimeScheduler 
///	Routine: ComputeTimeSteps
///	Returns: 
///	Action : Build a schedule based on different number
///          of steps per year for several periods :
///             - before 1st event date
///             - between the 1st event date and an "optimal" date
///             - after the "optimal" date
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MultiRegimeScheduler::ComputeTimeSteps(const ARM_PricingModel& model) const
{
    const ARM_GP_Vector& modelSchedule = *(model.GetNumMethod()->GetTimeSteps());
    size_t nbEvents = modelSchedule.size();
	
#ifdef __GP_STRICT_VALIDATION
    if(nbEvents < 2 || modelSchedule[0] != 0.0)
        ARM_THROW( ERR_INVALID_ARGUMENT,ARM_USERNAME+ "Model schedule must begin at spot time (t=0)" );
#endif
	
	
    ARM_GP_Vector* timeSteps = new ARM_GP_Vector;
	
    /// Compute steps before 1st event time
    double firstEventTime=modelSchedule[1];
//    int nbSteps=static_cast<int>( floor(firstEventTime/K_YEAR_LEN*itsStepNbPerYearBefore1stEventDate+0.5) );
    int nbSteps=static_cast<int>( floor(firstEventTime/K_YEAR_LEN*itsStepNbPerYearBefore1stEventDate) ); // integer part to compute as Tree3F does

    if(nbSteps < itsMinStepBefore1stEventDate )
        nbSteps = itsMinStepBefore1stEventDate;
    else if(nbSteps > itsMaxStepBefore1stEventDate )
        nbSteps = itsMaxStepBefore1stEventDate;
	
    double time,step = firstEventTime/nbSteps;
    int i;

	/// Force to insert 0 in any cases!
	timeSteps->push_back(0.0);
	for(i=0,time=step;i<nbSteps-1;++i,time += step)
        timeSteps->push_back(time);
    timeSteps->push_back(firstEventTime);

    /// If only one call the schedule is done
    if(nbEvents < 3)
        return ARM_GP_VectorPtr(timeSteps);

    bool isGPStyle=false;
    if(isGPStyle)
    {
        /// GP style for lattice schedule building

        /// Locate event time just before optimal time
        int optIdx = CC_NS(std,lower_bound)(modelSchedule.begin(), modelSchedule.end(), itsOptimalYearFracTime * K_YEAR_LEN ) - modelSchedule.begin()-1;
	    if( optIdx < 1 ) optIdx = 1;
        double optimalEventTime = modelSchedule[optIdx];

        /// Compute steps between 1st event time & optimal time
        if(optimalEventTime > firstEventTime)
            InsertCstSteps(timeSteps,1,optIdx,modelSchedule,itsStepNbPerYearBeforeOptimalTime);

        /// Compute steps after optimal time
        InsertCstSteps(timeSteps,optIdx,nbEvents-1,modelSchedule,itsStepNbPerYearAfterOptimalTime);
    }
    else
    {
        /// Summit production style for lattice schedule building
        double lastEventTime=modelSchedule[nbEvents-1];

        double optimalTime = itsOptimalYearFracTime*K_YEAR_LEN;
		optimalTime = lastEventTime < optimalTime ? lastEventTime : optimalTime;
        size_t nbStepBeforeOptimal = static_cast<size_t>( floor((optimalTime-firstEventTime)/K_YEAR_LEN*itsStepNbPerYearBeforeOptimalTime + 0.5) );
        double stepBeforeOptimal = (optimalTime - firstEventTime)/nbStepBeforeOptimal;
		
		double optimal2Time = its2OptimalYearFracTime*K_YEAR_LEN;
		optimal2Time = lastEventTime < optimal2Time ? lastEventTime : optimal2Time;
        size_t nbStepBetweenOptimal = static_cast<size_t>( floor((optimal2Time-optimalTime)/K_YEAR_LEN*itsStepNbPerYearAfterOptimalTime + 0.5) );
        double stepBetweenOptimal = (optimal2Time - optimalTime)/nbStepBetweenOptimal;

        size_t nbStepAfterOptimal2 = static_cast<size_t>( floor((lastEventTime-optimal2Time)/K_YEAR_LEN*itsStepNbPerYearAfter2OptimalTime + 0.5) );
        double stepAfterOptimal2 = (lastEventTime-optimal2Time)/nbStepAfterOptimal2;

        double time,nextTime;
        for(size_t eventIdx=1;eventIdx+1 < nbEvents;++eventIdx)
        {
            time = modelSchedule[eventIdx];
            /// Insert steps between each event using summit rule
            if(time < optimalTime)
            {
                /// Use the "after 1st notice & before optimal date gap"
                while((time = *(timeSteps->end()-1)) < modelSchedule[eventIdx+1] - K_NEW_DOUBLE_TOL)
                {
                    nextTime = time + stepBeforeOptimal;
                    if(nextTime < modelSchedule[eventIdx+1] - 0.5*stepBeforeOptimal)
                        timeSteps->push_back(nextTime);
                    else
                        timeSteps->push_back(modelSchedule[eventIdx+1]);
                }
            }
			else if(time < optimal2Time)
            {
                /// Use the "after 1st notice & before optimal date gap"
                while((time = *(timeSteps->end()-1)) < modelSchedule[eventIdx+1] - K_NEW_DOUBLE_TOL)
                {
                    nextTime = time + stepBetweenOptimal;
                    if(nextTime < modelSchedule[eventIdx+1] - 0.5*stepBetweenOptimal)
                        timeSteps->push_back(nextTime);
                    else
                        timeSteps->push_back(modelSchedule[eventIdx+1]);
                }
            }
            else
            {
                /// Use the "after optimal date gap"
                while((time = *(timeSteps->end()-1)) < modelSchedule[eventIdx+1] - K_NEW_DOUBLE_TOL)
                {
                    nextTime = time + stepAfterOptimal2;
                    if(nextTime < modelSchedule[eventIdx+1] - 0.5*stepAfterOptimal2)
                        timeSteps->push_back(nextTime);
                    else
                        timeSteps->push_back(modelSchedule[eventIdx+1]);
                }
            }
        } // for eventTimes
    }
    return ARM_GP_VectorPtr(timeSteps);
}


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
///	Class  : ARM_ConstantVarianceMeanRevertingScheduler
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_ConstantVarianceMeanRevertingScheduler
///	Routine: YtoX
///	Returns: double
///	Action : Function to invert from Y to X
////////////////////////////////////////////////////

double ARM_ConstantVarianceMeanRevertingScheduler::YtoX(ARM_GP_Vector& X,ARM_GP_Vector& Y,double y,int& lastIdx,int nextIdx) const
{
    /// Search x such that Y(x)= y, Y() is linearly interpolated
    /// Y is strictly increasing and we have Y[lastIdx] < y <= Y[nextIdx]
    for(int idx=lastIdx+1;idx<=nextIdx;++idx)
    {
        if(y<=Y[idx])
            break;
    }

    lastIdx=idx-1; /// for next call by increasing y
    double xinf=X[lastIdx];
    double yinf=Y[lastIdx];
    double xsup=X[idx];
    double ysup=Y[idx];
    return xinf + (y-yinf)*(xsup-xinf)/(ysup-yinf);
}



////////////////////////////////////////////////////
///	Class  : ARM_ConstantVarianceMeanRevertingScheduler
///	Routine: ComputeTimeSteps
///	Returns: ARM_GP_VectorPtr
///	Action : Computes the time steps. The logic is to have
///             -a minimum of steps before the first notice date
///             -local variance almost constant with the constraint to have the event date included
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_ConstantVarianceMeanRevertingScheduler::ComputeTimeSteps(const ARM_PricingModel& model) const
{
    const ARM_GP_Vector* const modelSchedule = model.GetNumMethod()->GetTimeSteps();
    int nbModelSteps= modelSchedule->size()-1; /// Actual model steps
    int nbSteps = model.GetNumMethod()->GetNbSteps();

    if(nbSteps-nbModelSteps <= 0)
        return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(modelSchedule->Clone()) );

    ARM_MatrixVector localVCV;
    model.NumMethodStateLocalVariances(*modelSchedule,localVCV);

    ARM_GP_Vector varIncr(localVCV.size());
    double varIncrSum=0.0;
    int i,j;
    for(i=0;i<localVCV.size();++i)
    {
        varIncr[i] = localVCV[i]->trace();
        varIncrSum += varIncr[i];
    }

    if(varIncrSum < itsMinStdDev*itsMinStdDev*(*modelSchedule)[nbModelSteps]/K_YEAR_LEN)
    {
        /// Variance too low, switch to a time
        /// uniform spacing
        return ARM_SchedulerBase::ComputeTimeSteps( model);
    }

    
    /// 1st step : set steps before the 1st model time and
    /// fulfill the constraint on the minimum number
    int nbInSteps = (int)(floor(varIncr[0]/varIncrSum*nbSteps));
    if(nbInSteps<1)
        nbInSteps=1; // current event date
    double time=0.0;
    double lastModelTime=0.0,modelTime;
    double var = 0.0;
    int ifirst=1,nbAdded=0;
    double step;

    ARM_GP_Vector timeSteps(nbSteps+nbModelSteps);
    timeSteps[0]=0.0;
    int idx=1;

    ARM_IntVector modelTimeIdx(nbModelSteps+1);
    modelTimeIdx[0]=0;

    if(nbInSteps <= itsMinStepBefore1stEventDate)
    {
        /// Force itsMinStepBefore1stEventDate steps between 0 and
        /// the first model time

        var += varIncr[ifirst-1];
        modelTime = (*modelSchedule)[ifirst];
        step = modelTime/itsMinStepBefore1stEventDate;
        for(;idx<itsMinStepBefore1stEventDate;++idx)
        {
            time += step;
            timeSteps[idx]=time;
        }
        time = modelTime;
        timeSteps[idx]=time;  // to avoid rounding pbs !
        modelTimeIdx[ifirst]=idx;
        ++idx;
        nbAdded = itsMinStepBefore1stEventDate;
        lastModelTime = modelTime;
        ifirst = 2;
    }

    /// 2nd step : insert steps between all model times
    for(i=ifirst;i<=nbModelSteps;++i)
    {
        var += varIncr[i-1] ;
        nbInSteps = (int)(floor(var/varIncrSum*nbSteps)) - nbAdded;
        if(nbInSteps<1)
            nbInSteps=1; // current event date
        if(nbAdded+nbInSteps > timeSteps.size())
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                " : Inconsistency in the mix time/var tree schedule initialisation" );
        }
        modelTime = (*modelSchedule)[i];
        step = (modelTime-lastModelTime)/nbInSteps;
        for(j=0;j<nbInSteps-1;++j,++idx)
        {
            time += step;
            timeSteps[idx]=time;
        }
        time = modelTime;
        timeSteps[idx]=time;  // to avoid rounding pbs !
        modelTimeIdx[i]=idx;
        ++idx;
        nbAdded += nbInSteps;
        lastModelTime = modelTime;
    }

    /// 3rd step : review this time schedule using the
    /// local variance sum profile which is now
    /// strictly increasing (but not the variance profile)
    nbSteps=idx;
    ARM_GP_Vector X(nbSteps),Y(nbSteps);
    Y[0]=0.0;
    for(j=0;j<nbSteps;++j)
        X[j]=timeSteps[j];

	/// management of simplematrix to triangular matrix
    ARM_MatrixVector lastLocalVCV;
	model.NumMethodStateLocalVariances(X,lastLocalVCV);
    
    for(j=1;j<nbSteps;++j)
        Y[j] = Y[j-1] + lastLocalVCV[j-1]->trace();

    ARM_GP_Vector* finalTimeSteps = new ARM_GP_Vector(nbSteps);
    (*finalTimeSteps)[0]=0.0;
    idx=1;
    var=0.0;
    int nextIdx,lastIdx=0;
    double nextVar,varStep;
    for(int nextModelIdx=1;nextModelIdx<=nbModelSteps;++nextModelIdx)
    {
        nextIdx=modelTimeIdx[nextModelIdx];
        nextVar=Y[nextIdx];
        nbSteps=nextIdx-lastIdx;
        varStep=(nextVar-var)/nbSteps;
        var += varStep;
        for(int i=0;i<nbSteps-1;++i)
        {
            time=YtoX(X,Y,var,lastIdx,nextIdx);
            (*finalTimeSteps)[idx]=time;
            var += varStep;
            ++idx;
        }
        if(idx != nextIdx)
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                " : Inconsistency in the final rescheduling in constant variance" );
        }
        time=(*modelSchedule)[nextModelIdx];
        (*finalTimeSteps)[idx]=time;
        lastIdx=nextIdx;
        var=nextVar;
        ++idx;
    }

    /// Free memory
	DeletePointorVector<ARM_GP_Matrix>( localVCV );
	DeletePointorVector<ARM_GP_Matrix>( lastLocalVCV );

    return static_cast<ARM_GP_VectorPtr>(finalTimeSteps);
}





////////////////////////////////////////////////////
///	Class  : ARM_ConstantVarianceScheduler
///	Routine: ComputeTimeSteps
///	Returns: ARM_GP_VectorPtr
///	Action : Computes the time steps. The logic is to have
///             -a minimum of steps before the first notice date
///             -constant variance with the constraint to have the event date included
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_ConstantVarianceScheduler::ComputeTimeSteps(const ARM_PricingModel& model) const
{
    const ARM_GP_Vector* const modelSchedule = model.GetNumMethod()->GetTimeSteps();
    int nbModelSteps= modelSchedule->size()-1; /// Actual model steps
    int nbSteps     = model.GetNumMethod()->GetNbSteps();
    int nbToAdd=nbSteps-(nbModelSteps); 
    if(nbToAdd <= 0)
		return ARM_GP_VectorPtr(new ARM_GP_Vector(*modelSchedule));
       
	/// management of simplematrix to triangular matrix
    ARM_MatrixVector localVCV;
	ARM_MatrixVector globalVCV;
   
    model.NumMethodStateLocalGlobalVariances(*modelSchedule,localVCV,globalVCV);
	double firstModelTime=(*modelSchedule)[1];
    double firstVar = globalVCV[1]->trace();
    double lastModelTime=(*modelSchedule)[nbModelSteps];
    double lastVar = globalVCV[nbModelSteps]->trace();

    ARM_GP_Vector timeSteps(nbSteps+1); // +1 for [0]=0
    timeSteps[0]=0.0;
    int idx=1;

    double time,timeErr;
    double varErrRel=0.0001;

    double var;

    double varStep=lastVar/nbToAdd;
    int nbStepBefore=(int)(floor(firstVar/varStep));
    time=model.VarianceToTime(nbStepBefore*varStep,0.0,firstModelTime);

    if(time < firstModelTime - varStep/10.0)
    {
        /// Don't merge additional step
        nbStepBefore++;
    }

    ARM_IntVector modelTimeIdx(nbModelSteps+1);
    modelTimeIdx[0]=0;

    /// 1st stage : set steps before the 1st model time with
    /// fulfilling the constraint on the minimum number
    if(nbStepBefore <= itsMinStepBefore1stEventDate)
    {
        /// Force itsMinStepBefore1stEventDate steps between 0 and
        /// the first model time

        double varStepBefore=firstVar/itsMinStepBefore1stEventDate;
        var=varStepBefore;
        time=model.VarianceToTime(var,0.0,firstModelTime);
        for(idx=1;idx<itsMinStepBefore1stEventDate;++idx)
        {
            timeSteps[idx]=time;
            var+=varStepBefore;
            time=model.VarianceToTime(var,time,firstModelTime);
        }
        timeErr=model.VarianceToTime(var*(1+varErrRel),time,time+1.0)-time;
        if(time > firstModelTime + timeErr)
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                " : Inconsistency before first model time in the tree schedule initialisation" );
        }

        varStep=(lastVar-firstVar)/(nbToAdd-itsMinStepBefore1stEventDate+1);
    }
    else if(nbStepBefore > 1) 
    {
        /// Keep the same spacing in variance before and
        /// after the first model time
        var=varStep;
        time=model.VarianceToTime(var,0.0,firstModelTime);
        for(idx=1;idx<nbStepBefore;++idx)
        {
            timeSteps[idx]=time;
            var+=varStep;
            time=model.VarianceToTime(var,time,firstModelTime);
            timeErr=model.VarianceToTime(var*(1+varErrRel),time,firstModelTime)-time;
        }
        if(time <= firstModelTime + timeErr)
        {
            /// Next time step to add is the first model time step
            var=firstVar;
            time=firstModelTime;
        }
    }


    /// 2nd stage : insert steps between all model times
    /// using the global constant variance step
    double nextVar=firstVar;
    double nextModelTime=firstModelTime;
    int nextModelIdx=1;
    timeErr=model.VarianceToTime(firstVar*(1+varErrRel),firstModelTime,firstModelTime+1.0)-
        firstModelTime;
    for(;idx<=nbSteps;++idx)
    {
        if(nextModelTime - timeErr <= time)
        {
            // After next model time then stop just on it
            timeSteps[idx]=nextModelTime;
            modelTimeIdx[nextModelIdx]=idx;

            if(nextModelIdx < nbModelSteps)
            {
                if(time <= nextModelTime + timeErr)
                {
                    // Merge with the current model time in fact
                    var = nextVar + varStep;
                    time=model.VarianceToTime(var,time,(*modelSchedule)[nextModelIdx+1]);
                    timeErr=model.VarianceToTime(var*(1+varErrRel),time,time+1.0)-
                        time;
                }

                // Next model time is the new target to reach
                nextModelIdx++;
                nextModelTime=(*modelSchedule)[nextModelIdx];
                nextVar=globalVCV[nextModelIdx]->trace();
            }
            else
            {
                // Last model time is reached
				nextModelIdx++;
				idx++;
                break;
            }
        }
        else
        {
            // Strictly before next model time => a step is added
            timeSteps[idx]=time;
            var+=varStep;
			if (var >= lastVar-K_NEW_DOUBLE_TOL)
			{
				idx++;
				break;
			}
            time=model.VarianceToTime(var,time,nextModelTime);
            timeErr=model.VarianceToTime(var*(1+varErrRel),time,time+1.0)-time;
        }
    } // idx <= nbSteps

	// If we go out of the loop we need at least to had the model timesteps.
	for (;nextModelIdx < nbModelSteps+1;++nextModelIdx,++idx)
		timeSteps[idx]=(*modelSchedule)[nextModelIdx];


    /// Copy computed schedule to the final one (skipping the 3rd uniformisation stage)
    ARM_GP_Vector* finalTimeSteps = new ARM_GP_Vector(idx<nbSteps+1?idx:nbSteps+1);
    for(idx=0;idx<finalTimeSteps->size();++idx)
        (*finalTimeSteps)[idx]=timeSteps[idx];


/********** commented to be consistent with SFRM tree scheduler **********
    /// 3rd stage : reschedule evenly in variance between each
    /// model time with the number of time steps computed in 2nd stage
    int lastIdx=0;
    var=0.0;
    time=0.0;
    (*finalTimeSteps)[lastIdx]=0.0;
    idx=1;
    for(nextModelIdx=1;nextModelIdx<=nbModelSteps;++nextModelIdx)
    {
        nextVar=globalVCV[nextModelIdx]->trace();
        nbSteps=modelTimeIdx[nextModelIdx];
        nbSteps -= lastIdx;
        varStep=(nextVar-var)/nbSteps;
        nextModelTime=(*modelSchedule)[nextModelIdx];
        var += varStep;
        for(int i=0;i<nbSteps-1;++i)
        {
            time=model.VarianceToTime(var,time,nextModelTime);
            (*finalTimeSteps)[idx]=time;
            var += varStep;
            ++idx;
        }
        if(idx != modelTimeIdx[nextModelIdx])
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                " : Inconsistency in the final rescheduling in constant variance" );
        }
        time=nextModelTime;
        (*finalTimeSteps)[idx]=time;
        lastIdx=idx;
        var=nextVar;
        ++idx;
    }

********** commented to be consistent with SFRM tree scheduler **********/

    /// Free memory
	DeletePointorVector<ARM_GP_Matrix>( localVCV );
	DeletePointorVector<ARM_GP_Matrix>( globalVCV );

    return static_cast<ARM_GP_VectorPtr>(finalTimeSteps);
}

////////////////////////////////////////////////////
///	Class  : ARM_TimeStepPerYearScheduler
///	Routine: ComputeTimeSteps
///	Returns: ARM_GP_VectorPtr
///	Action : Computes the time steps. The logic is to have
///             n time steps per year
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_TimeStepPerYearScheduler::ComputeTimeSteps(const ARM_PricingModel& model) const
{
	if( itsTimeStepPerYear >= 0 )
	{
		int fixStep = model.ModelFixTimeStep( itsTimeStepPerYear );
		ARM_GP_Vector& modelSchedule = *(model.GetNumMethod()->GetTimeSteps());
		ARM_FixStepTimeInserter timeInserter(fixStep);
		ARM_GP_VectorPtr pNewTimeSteps( timeInserter.AddTimes(modelSchedule) );
		return pNewTimeSteps;
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                " : TimeStepPerYear should be positive." );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_TimeStepPerYearScheduler
///	Routine: toString
///	Returns: string
///	Action : Display the contents
////////////////////////////////////////////////////

string ARM_TimeStepPerYearScheduler::toString(const string& indent,const string& nextIndent) const
{ 
	CC_Ostringstream os;

	os << "Time Step Per Year Scheduler: " << itsTimeStepPerYear;

	return os.str();
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/