/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file schedulerfactory.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPNUMMETHODS_SCHEDULERFACTORY_H
#define _INGPNUMMETHODS_SCHEDULERFACTORY_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
struct ARM_SchedulerBase;

class ARM_SchedulerFactoryData : public ARM_RootObject
{
private:
    int itsSchedulerType;
    std::vector<double> itsSchedulerDatas;

public:
    ARM_SchedulerFactoryData(int schedulerType, const std::vector<double>& schedulerDatas)
        : itsSchedulerType(schedulerType), itsSchedulerDatas(schedulerDatas)
    {}

    virtual ~ARM_SchedulerFactoryData() {}
	virtual ARM_Object* Clone() const { return new ARM_SchedulerFactoryData(*this); }

    virtual string toString(const string& indent="",const string& nextIndent="") const;

    int GetSchedulerType() const { return itsSchedulerType; }
    const std::vector<double>& GetSchedulerDatas() const { return itsSchedulerDatas; }
};


struct ARM_SchedulerFactoryImp
{
    ARM_SchedulerBase* CreateScheduler(
		int schedulerType, 
		const std::vector<double>& schedulerDatas,
		// It need to be return to be set in the numerical method
		int& nbSteps);

    ARM_SchedulerBase* CreateScheduler(const ARM_SchedulerFactoryData& schedulerData, int& nbSteps);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_SchedulerFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_SchedulerFactoryImp>;
};

extern ARM_SingletonHolder<ARM_SchedulerFactoryImp> ARM_SchedulerFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

