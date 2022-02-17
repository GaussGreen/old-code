/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file scheduler.h
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */

#ifndef _INGPNUMMETHODS_SCHEDULERBASE_H
#define _INGPNUMMETHODS_SCHEDULERBASE_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_SchedulerBase : public ARM_RootObject
{
    enum SchedulerType
	{
		ConstantVariance = 0,
		ConstantVarianceMeanReverting,
        MultiRegime,
		TimeStepPerYear
    };

    virtual ARM_GP_VectorPtr ComputeTimeSteps(const ARM_PricingModel& model ) const;
};


CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/