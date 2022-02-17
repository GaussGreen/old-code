/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file scheduler.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPNUMMETHODS_SCHEDULER_H
#define _INGPNUMMETHODS_SCHEDULER_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpinfra/typedef.h"
#include "gpinfra/schedulerbase.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;

///////////////////////////////////////////////////////////////
/// scheduler used originally in the markovian drift sampler
///////////////////////////////////////////////////////////////
struct ARM_MultiRegimeScheduler : public ARM_SchedulerBase
{
private:
    size_t itsMinStepBefore1stEventDate;
    size_t itsMaxStepBefore1stEventDate;
	size_t itsStepNbPerYearBefore1stEventDate;
    double itsOptimalYearFracTime;
    size_t itsStepNbPerYearBeforeOptimalTime;
	double its2OptimalYearFracTime;
	size_t itsStepNbPerYearAfterOptimalTime;
    size_t itsStepNbPerYearAfter2OptimalTime;
    void InsertCstSteps(ARM_GP_Vector* timeSteps,int startIdx,int endIdx,const ARM_GP_Vector& modelSchedule,size_t nbStepPerYear) const;

public:
	ARM_MultiRegimeScheduler( 
		size_t minStepBefore1stEventDate        = 1,
		size_t maxStepBefore1stEventDate        = 10,
		size_t stepNbPerYearBefore1stEventDate  = 4,
		double optimalYearFracTime              = 7.0,
		size_t stepNbPerYearBeforeOptimalTime   = 2,
		size_t stepNbPerYearAfterOptimalTime    = 1,
		double optimal2YearFracTime             = 1000.0,
		size_t stepNbPerYearAfter2OptimalTime	= 1)
	:
		ARM_SchedulerBase(),
		itsMinStepBefore1stEventDate( minStepBefore1stEventDate ),
		itsMaxStepBefore1stEventDate( maxStepBefore1stEventDate ),
		itsStepNbPerYearBefore1stEventDate( stepNbPerYearBefore1stEventDate ),
		itsOptimalYearFracTime( optimalYearFracTime ),
		itsStepNbPerYearBeforeOptimalTime( stepNbPerYearBeforeOptimalTime ),
		itsStepNbPerYearAfterOptimalTime( stepNbPerYearAfterOptimalTime ),
		its2OptimalYearFracTime( optimal2YearFracTime ),
		itsStepNbPerYearAfter2OptimalTime( stepNbPerYearAfter2OptimalTime )
	{
	}

	/// function to compute time steps
	virtual ARM_GP_VectorPtr ComputeTimeSteps(const ARM_PricingModel& model ) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_MultiRegimeScheduler(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MultiRegimeScheduler";}
};

///////////////////////////////////////////////////////////////
/// scheduler used originally in the mean reverting scheduler
///////////////////////////////////////////////////////////////
struct ARM_ConstantVarianceMeanRevertingScheduler : public ARM_SchedulerBase
{
private:
    size_t itsMinStepBefore1stEventDate;
    double itsMinStdDev;
	double YtoX(std::vector<double>& X,std::vector<double>& Y,double y,int& lastIdx,int nextIdx) const;

public:
	ARM_ConstantVarianceMeanRevertingScheduler(
        size_t minStepBefore1stEventDate    = 1,
        double minStdDev                    = 1.0e-3 ) 
	: itsMinStepBefore1stEventDate( minStepBefore1stEventDate ), itsMinStdDev( minStdDev ) 
	{}

	/// function to compute time steps
	virtual ARM_GP_VectorPtr ComputeTimeSteps(const ARM_PricingModel& model ) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_ConstantVarianceMeanRevertingScheduler(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_ConstantVarianceMeanRevertingScheduler";}
};


///////////////////////////////////////////////////////////////
/// scheduler used originally in the normal centered sampler
///////////////////////////////////////////////////////////////
struct ARM_ConstantVarianceScheduler : public ARM_SchedulerBase
{
private:
    size_t itsMinStepBefore1stEventDate;
public:
	ARM_ConstantVarianceScheduler( size_t minStepBefore1stEventDate = 1 ) : itsMinStepBefore1stEventDate(minStepBefore1stEventDate) {}

	/// function to compute time steps
	virtual ARM_GP_VectorPtr ComputeTimeSteps(const ARM_PricingModel& model ) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_ConstantVarianceScheduler(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_ConstantVarianceScheduler";}
};


///////////////////////////////////////////////////////////////
/// scheduler used for MC or PDE
/// It creates some time steps per year
///////////////////////////////////////////////////////////////
struct ARM_TimeStepPerYearScheduler : public ARM_SchedulerBase
{
private:
	size_t itsTimeStepPerYear;
public:
	ARM_TimeStepPerYearScheduler( size_t timeStepPerYear = 1 ) : itsTimeStepPerYear(timeStepPerYear) {}

	/// function to compute time steps
	virtual ARM_GP_VectorPtr ComputeTimeSteps(const ARM_PricingModel& model ) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_TimeStepPerYearScheduler(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

